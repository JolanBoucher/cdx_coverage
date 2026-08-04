"""Standalone circular coverage plot renderer, invoked as a subprocess by cdx_coverage.

All numeric preparation (downsampling, smoothing, log transform, tick math,
grid layout) is done in C++ (src/python_circular_plot.cpp), reusing the same
numeric core as the linear (Cairo) backend - that work is CPU-bound and
scales with full-resolution coverage length, so it doesn't belong in a
subprocess. This script receives already-processed points and ticks and is
responsible only for the pycirclize/matplotlib rendering itself: coordinate
re-projection onto the polar track, segment/gap drawing, axis labels, title,
and PNG export - ported from the Python prototype's
output_tools.py:_plotCircularGraph.

Only circular graphs go through Python. Linear graphs remain pure C++/Cairo
(no subprocess, no Python dependency) in src/cairo_plot.cpp.

Global mode (many components) renders each panel in its own process
(multiprocessing.Pool, one panel ~1s of matplotlib/pycirclize work that a
single process can't parallelize across cores due to the GIL) into a
fixed-size PNG, then composites them into the final grid via matplotlib's
own PNG reader/writer (matplotlib.image) - no Pillow dependency needed for
that. This cut a 17-component render from ~19s to a few seconds on a
multi-core machine. Query mode (a single panel) renders directly, no pool.

--------------------------------------------------------------------------
Binary input format written by src/python_circular_plot.cpp (little-endian,
format version 2):

Header (common to both modes):
    magic          : 4 bytes, ASCII "CXCP"
    version        : uint8   (=2)
    mode           : uint8   (0 = query, 1 = global)
    dpi            : uint32
    fig_width      : float64
    fig_height     : float64
    line_color     : uint16 length-prefix + UTF-8 bytes
    fill_color     : uint16 length-prefix + UTF-8 bytes

Query payload (mode = 0): one Package (see below).

Global payload (mode = 1):
    subplot_width   : float64
    subplot_height  : float64
    rows            : uint32
    columns         : uint32
    component_count : uint32
    packages        : component_count * Package

Package (one component or one query, already fully processed):
    component_name    : uint16 length-prefix + UTF-8 bytes
    compo_length       : uint64
    query_start         : uint64
    query_end           : uint64
    full_component      : uint8
    crosses_origin       : uint8
    visible              : uint8   (0 => nothing to draw, hide the axis)
    logarithmic          : uint8
    log_base              : int32
    y_upper_limit         : float64
    point_count            : uint64
    local_x                : point_count * uint64   (indices into the traversal)
    plot_y                 : point_count * float64  (may contain NaN = gap)
    tick_count              : uint32
    tick_raw_values         : tick_count * float64
    tick_display_values     : tick_count * float64
--------------------------------------------------------------------------
"""

import os
import struct
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from pycirclize import Circos


# ---------------------------------------------------------------------------
# Binary reader
# ---------------------------------------------------------------------------

class _Reader:
    """Small sequential binary reader matching the format documented above."""

    def __init__(self, f):
        self.f = f

    def u8(self) -> int:
        return struct.unpack("<B", self.f.read(1))[0]

    def i32(self) -> int:
        return struct.unpack("<i", self.f.read(4))[0]

    def u32(self) -> int:
        return struct.unpack("<I", self.f.read(4))[0]

    def u64(self) -> int:
        return struct.unpack("<Q", self.f.read(8))[0]

    def f64(self) -> float:
        return struct.unpack("<d", self.f.read(8))[0]

    def string(self) -> str:
        length = struct.unpack("<H", self.f.read(2))[0]
        return self.f.read(length).decode("utf-8") if length else ""

    def u64_array(self, count: int) -> np.ndarray:
        return np.fromfile(self.f, dtype="<u8", count=count)

    def f64_array(self, count: int) -> np.ndarray:
        return np.fromfile(self.f, dtype="<f8", count=count)

    def package(self) -> dict:
        component_name = self.string()
        compo_length = self.u64()
        query_start = self.u64()
        query_end = self.u64()
        full_component = bool(self.u8())
        crosses_origin = bool(self.u8())
        visible = bool(self.u8())
        logarithmic = bool(self.u8())
        log_base = self.i32()
        y_upper_limit = self.f64()

        point_count = self.u64()
        local_x = self.u64_array(point_count)
        plot_y = self.f64_array(point_count)

        tick_count = self.u32()
        tick_raw_values = self.f64_array(tick_count)
        tick_display_values = self.f64_array(tick_count)

        return dict(
            component_name=component_name,
            compo_length=compo_length,
            query_start=query_start,
            query_end=query_end,
            full_component=full_component,
            crosses_origin=crosses_origin,
            visible=visible,
            logarithmic=logarithmic,
            log_base=log_base,
            y_upper_limit=y_upper_limit,
            local_x=local_x,
            plot_y=plot_y,
            tick_raw_values=tick_raw_values,
            tick_display_values=tick_display_values,
        )


def _load_request(path) -> dict:
    with open(path, "rb") as fh:
        r = _Reader(fh)

        magic = fh.read(4)
        if magic != b"CXCP":
            raise ValueError(f"Invalid magic number: {magic!r}, expected b'CXCP'")

        version = r.u8()
        if version != 2:
            raise ValueError(f"Unsupported binary format version: {version} (expected 2)")

        mode = r.u8()
        dpi = r.u32()
        fig_width = r.f64()
        fig_height = r.f64()
        line_color = r.string()
        fill_color = r.string()

        common = dict(
            dpi=dpi,
            fig_width=fig_width,
            fig_height=fig_height,
            line_color=line_color,
            fill_color=fill_color,
        )

        if mode == 0:
            return dict(mode="query", package=r.package(), **common)

        if mode == 1:
            subplot_width = r.f64()
            subplot_height = r.f64()
            rows = r.u32()
            columns = r.u32()
            component_count = r.u32()
            packages = [r.package() for _ in range(component_count)]

            return dict(
                mode="global",
                subplot_width=subplot_width,
                subplot_height=subplot_height,
                rows=rows,
                columns=columns,
                packages=packages,
                **common,
            )

        raise ValueError(f"Unknown mode byte: {mode}")


# ---------------------------------------------------------------------------
# Core renderer - pycirclize/matplotlib only, all numeric prep already done.
# Ported from output_tools.py:_plotCircularGraph, steps 6-8 (canvas/track
# setup, coordinate re-projection & drawing, annotations & ticks); steps 1-5
# (traversal extraction, masking, smoothing, scale transform, tick math) now
# happen in src/python_circular_plot.cpp:prepareCircularPlotPackage.
# ---------------------------------------------------------------------------

def _plot_circular_graph(axis, pkg: dict, query_mode: bool, line_color: str, fill_color: str) -> None:
    if not pkg["visible"]:
        axis.set_visible(False)
        return

    component_id = pkg["component_name"]
    compo_start, compo_end = 0, int(pkg["compo_length"])
    compo_length = compo_end - compo_start

    query_start = int(pkg["query_start"])
    query_end = int(pkg["query_end"])
    full_component = pkg["full_component"]
    crosses_origin = pkg["crosses_origin"]

    local_x = pkg["local_x"].astype(np.float64, copy=False)
    plot_y = pkg["plot_y"]
    coverage_upper_limit = pkg["y_upper_limit"]
    raw_ticks = pkg["tick_raw_values"]
    display_ticks = pkg["tick_display_values"]
    log_base = pkg["log_base"] if pkg["logarithmic"] else None

    if local_x.size == 0 or np.all(np.isnan(plot_y)):
        axis.set_visible(False)
        return

    circos = Circos(sectors={str(component_id): (compo_start, compo_end)})
    sector = circos.sectors[0]

    inner_radius, outer_radius = 40.0, 95.0
    radial_size = outer_radius - inner_radius

    track = sector.add_track((inner_radius, outer_radius), r_pad_ratio=0.0)
    track.axis(fc="#F8FAFC", ec="#CBD5E1", lw=0.5)

    def draw_segment(x: np.ndarray, y: np.ndarray) -> None:
        """Draw contiguous finite segments to prevent rendering gaps across NaNs."""
        if x.size == 0 or y.size == 0:
            return
        finite_mask = np.isfinite(y)
        if not np.any(finite_mask):
            return

        finite_indices = np.flatnonzero(finite_mask)
        split_positions = np.where(np.diff(finite_indices) > 1)[0] + 1

        for block in np.split(finite_indices, split_positions):
            if block.size == 0:
                continue
            segment_x = x[block]
            segment_y = y[block]

            track.fill_between(
                segment_x, segment_y, 0,
                vmin=0.0, vmax=coverage_upper_limit,
                color=fill_color, alpha=0.35,
            )
            track.line(
                segment_x, segment_y,
                vmin=0.0, vmax=coverage_upper_limit,
                color=line_color, lw=1.0,
            )

    if full_component:
        display_x = np.append(local_x + compo_start, compo_end)
        display_y = np.append(plot_y, plot_y[0])
        draw_segment(display_x, display_y)
    elif not crosses_origin:
        display_x = np.append(local_x + query_start, query_end + 1)
        display_y = np.append(plot_y, plot_y[-1])
        draw_segment(display_x, display_y)
    else:
        first_segment_length = compo_end - query_start
        first_mask = local_x < first_segment_length
        second_mask = ~first_mask

        first_x = local_x[first_mask] + query_start
        first_y = plot_y[first_mask]
        second_x = local_x[second_mask] - first_segment_length + compo_start
        second_y = plot_y[second_mask]

        finite_mask = np.isfinite(plot_y)
        if np.count_nonzero(finite_mask) >= 2:
            origin_value = float(np.interp(first_segment_length, local_x[finite_mask], plot_y[finite_mask]))
        else:
            origin_value = np.nan

        if first_x.size:
            first_x = np.append(first_x, compo_end)
            first_y = np.append(first_y, origin_value)
            draw_segment(first_x, first_y)

        if second_x.size:
            second_x = np.insert(second_x, 0, compo_start)
            second_y = np.insert(second_y, 0, origin_value)
            second_x = np.append(second_x, query_end + 1)
            second_y = np.append(second_y, second_y[-1])
            draw_segment(second_x, second_y)

    tick_interval = max(1, compo_length // 6)
    tick_positions = np.arange(compo_start, compo_end, tick_interval, dtype=np.int64)

    track.xticks(
        tick_positions,
        [f"{position:,}" for position in tick_positions],
        label_size=7,
        label_orientation="horizontal",
    )

    if full_component:
        query_label = "Complete component"
    elif crosses_origin:
        query_label = f"Query: {query_start:,}-{compo_end - 1:,}; {compo_start:,}-{query_end:,} bp"
    else:
        query_label = f"Query: {query_start:,}-{query_end:,} bp"

    scale_label = f"Log scale: base {log_base}" if log_base is not None else "Linear scale"

    if query_mode:
        center_text = f"Component {component_id}\nLength: {compo_length:,} bp\n{query_label}\n"
    else:
        center_text = f"Component {component_id}\n{compo_length:,} bp\n{scale_label}"

    circos.text(center_text, r=0, size=8, weight="bold")
    circos.plotfig(ax=axis)

    angles = np.linspace(0.0, 2.0 * np.pi, 720)
    label_angle = np.deg2rad(337.5)

    for raw_tick, display_tick in zip(raw_ticks, display_ticks):
        if log_base is not None and raw_tick == 0:
            continue

        radius = inner_radius + (display_tick / coverage_upper_limit) * radial_size

        axis.plot(angles, np.full(angles.shape, radius), color="#94A3B8", lw=0.45, ls="--", alpha=0.60, zorder=1)
        axis.text(
            label_angle, radius, f"{raw_tick:,.0f}x",
            fontsize=6, color="#334155", ha="left", va="center", clip_on=False, zorder=10,
        )


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def _render_query(request: dict, output_png: Path) -> None:
    figure = plt.figure(figsize=(request["fig_width"], request["fig_height"]), dpi=request["dpi"])
    axis = figure.add_subplot(1, 1, 1, projection="polar")

    _plot_circular_graph(
        axis=axis,
        pkg=request["package"],
        query_mode=True,
        line_color=request["line_color"],
        fill_color=request["fill_color"],
    )

    figure.tight_layout()
    try:
        figure.savefig(output_png, dpi=request["dpi"], bbox_inches="tight", facecolor="white")
    finally:
        plt.close(figure)


def _render_panel_worker(args: tuple) -> str | None:
    """Renders ONE component's polar panel to its own fixed-size PNG.

    Runs in a worker process (multiprocessing.Pool): each pycirclize/
    matplotlib panel is independent CPU-bound work (~1s/panel measured),
    and Python's GIL means a single process can't parallelize that across
    cores even with threads - hence a process pool. Must be a module-level
    function (not a closure) to stay picklable under the "spawn" start
    method (the default on macOS).

    No `bbox_inches="tight"`: every panel must come out at EXACTLY
    `panel_px x panel_px` pixels for the grid composite in _render_global to
    line cells up correctly. This also happens to be faster (tight bbox
    computation requires an extra full render pass).
    """
    pkg, line_color, fill_color, dpi, panel_size_in, panel_path = args

    if not pkg["visible"]:
        return None

    figure = plt.figure(figsize=(panel_size_in, panel_size_in), dpi=dpi)
    axis = figure.add_subplot(1, 1, 1, projection="polar")

    _plot_circular_graph(
        axis=axis,
        pkg=pkg,
        query_mode=False,
        line_color=line_color,
        fill_color=fill_color,
    )

    try:
        figure.savefig(panel_path, dpi=dpi, facecolor="white")
    finally:
        plt.close(figure)

    return panel_path


def _render_global(request: dict, output_png: Path, work_dir: Path) -> None:
    packages = request["packages"]
    rows, columns = request["rows"], request["columns"]
    component_count = len(packages)
    dpi = request["dpi"]

    subplot_size = max(request["subplot_width"], request["subplot_height"], 5.5)
    panel_px = max(1, round(subplot_size * dpi))

    tasks = [
        (pkg, request["line_color"], request["fill_color"], dpi, subplot_size, str(work_dir / f"panel_{i}.png"))
        for i, pkg in enumerate(packages)
    ]

    # cpu_count() can return None on some platforms; component_count caps
    # worker count so we never spin up more processes than there is work.
    worker_count = max(1, min(component_count, os.cpu_count() or 1))

    if worker_count == 1:
        # Skip Pool entirely for a single panel/core: avoids paying a second
        # process's interpreter+import startup cost for zero parallel gain.
        panel_paths = [_render_panel_worker(task) for task in tasks]
    else:
        import multiprocessing
        with multiprocessing.Pool(processes=worker_count) as pool:
            panel_paths = pool.map(_render_panel_worker, tasks)

    # Composite panels into the final grid via matplotlib's own PNG
    # reader/writer (no Pillow dependency needed - PNG doesn't require it).
    # White background so empty trailing grid cells and invisible/masked
    # components render as blank cells, exactly like set_visible(False) did
    # in the previous single-figure approach.
    canvas = np.ones((panel_px * rows, panel_px * columns, 3), dtype=np.float32)

    for i, panel_path in enumerate(panel_paths):
        if panel_path is None:
            continue
        panel = mpimg.imread(panel_path)
        if panel.shape[2] == 4:
            panel = panel[:, :, :3]
        height, width = panel.shape[:2]
        row, col = divmod(i, columns)
        top, left = row * panel_px, col * panel_px
        canvas[top:top + height, left:left + width, :] = panel

    mpimg.imsave(output_png, canvas)


def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.bin> <output.png>", file=sys.stderr)
        return 2

    input_path = Path(sys.argv[1])
    output_png = Path(sys.argv[2])

    request = _load_request(input_path)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    if request["mode"] == "query":
        _render_query(request, output_png)
    else:
        _render_global(request, output_png, input_path.parent)

    return 0


if __name__ == "__main__":
    sys.exit(main())
