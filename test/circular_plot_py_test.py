#!/usr/bin/env python3
"""Standalone unit tests for python_script/circular_plot.py, independent of the C++ suite.

circular_plot.py does `from pycirclize import Circos` at import time, so it cannot even be
imported without the real pycirclize package installed. To test it in isolation - without
requiring pycirclize (which may not be installed, e.g. no bundled venv has been provisioned
yet) - a minimal fake "pycirclize" module is injected into sys.modules *before* importing
circular_plot. matplotlib and numpy, by contrast, are real imports: circular_plot.py actually
uses their APIs (np.isfinite, np.interp, np.split, ...) in ways worth exercising for real, and
they are lightweight/commonly available. If either is missing, every test class here is
skipped (rather than erroring) via _DEPENDENCIES_AVAILABLE.

Coverage, by test class:
  - ReaderTest: byte-level round trip of every _Reader primitive (u8/i32/u32/u64/f64/string/
    u64_array/f64_array) against hand-built buffers.
  - LoadRequestTest: _load_request() for query mode and global mode against hand-built request
    buffers matching the documented binary format, plus the three validation error paths (bad
    magic, unsupported version, unknown mode byte).
  - PlotCircularGraphEarlyExitTest: the two early-return branches of _plot_circular_graph()
    (pkg["visible"] is False; the downsampled series is empty or entirely NaN) - verifies
    axis.set_visible(False) is called and that Circos is never constructed (no rendering work
    is attempted).
  - DrawSegmentSplittingTest: the NaN-gap segment-splitting logic (nested draw_segment() inside
    _plot_circular_graph), exercised for all three traversal topologies - full component,
    contiguous sub-range, and origin-crossing - via a mocked Circos/sector/track chain so no
    real polar rendering happens. All expected coordinates below were hand-derived by tracing
    the exact algorithm in circular_plot.py (see comments on each test).
  - RadialTicksTest: the radial coverage-tick loop - specifically that the raw_tick == 0 entry
    is skipped on a logarithmic scale (it isn't representable there) but kept on a linear scale.

Run directly with:
    python3 -m unittest test/circular_plot_py_test.py -v
or let CMake/ctest run it (see the optional add_test() in test/CMakeLists.txt, gated on
Python3 being found at configure time).
"""

from __future__ import annotations

import io
import struct
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

# ---------------------------------------------------------------------------
# Make circular_plot.py importable without a real pycirclize installation.
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent.parent / "python_script"
sys.path.insert(0, str(_SCRIPT_DIR))

_fake_pycirclize = types.ModuleType("pycirclize")
_fake_pycirclize.Circos = mock.MagicMock(name="Circos_module_default")
# Overwrite (not setdefault): force our fake to be used even if a real
# pycirclize happens to be installed on this machine, so these tests behave
# identically and deterministically everywhere.
sys.modules["pycirclize"] = _fake_pycirclize

_DEPENDENCIES_AVAILABLE = True
_SKIP_REASON = ""
try:
    import circular_plot as cp  # noqa: E402  (import after sys.path/sys.modules setup above)
except ModuleNotFoundError as exc:  # pragma: no cover - only hit on machines missing numpy/matplotlib
    _DEPENDENCIES_AVAILABLE = False
    _SKIP_REASON = f"circular_plot.py could not be imported (missing dependency: {exc})"
    cp = None  # type: ignore[assignment]


def _skip_unless_deps():
    return unittest.skipUnless(_DEPENDENCIES_AVAILABLE, _SKIP_REASON)


# ---------------------------------------------------------------------------
# _Reader primitives
# ---------------------------------------------------------------------------

@_skip_unless_deps()
class ReaderTest(unittest.TestCase):
    """Byte-level round trip of every _Reader primitive against hand-built buffers."""

    def test_u8_reads_single_unsigned_byte(self):
        reader = cp._Reader(io.BytesIO(bytes([5])))
        self.assertEqual(reader.u8(), 5)

    def test_u8_reads_max_value_255(self):
        reader = cp._Reader(io.BytesIO(bytes([255])))
        self.assertEqual(reader.u8(), 255)

    def test_i32_reads_negative_value(self):
        reader = cp._Reader(io.BytesIO(struct.pack("<i", -123)))
        self.assertEqual(reader.i32(), -123)

    def test_i32_reads_positive_value(self):
        reader = cp._Reader(io.BytesIO(struct.pack("<i", 2_000_000_000)))
        self.assertEqual(reader.i32(), 2_000_000_000)

    def test_u32_reads_value_beyond_int32_range(self):
        reader = cp._Reader(io.BytesIO(struct.pack("<I", 4_000_000_000)))
        self.assertEqual(reader.u32(), 4_000_000_000)

    def test_u64_reads_large_value(self):
        reader = cp._Reader(io.BytesIO(struct.pack("<Q", 12_345_678_901_234)))
        self.assertEqual(reader.u64(), 12_345_678_901_234)

    def test_f64_reads_finite_value(self):
        reader = cp._Reader(io.BytesIO(struct.pack("<d", 3.14159)))
        self.assertAlmostEqual(reader.f64(), 3.14159)

    def test_f64_round_trips_nan(self):
        # plot_y arrays legitimately contain NaN (masked/out-of-query gaps),
        # so the reader must round-trip it rather than choke on it.
        reader = cp._Reader(io.BytesIO(struct.pack("<d", float("nan"))))
        self.assertTrue(np.isnan(reader.f64()))

    def test_string_reads_length_prefixed_utf8(self):
        payload = "chr1".encode("utf-8")
        buffer = struct.pack("<H", len(payload)) + payload
        reader = cp._Reader(io.BytesIO(buffer))
        self.assertEqual(reader.string(), "chr1")

    def test_string_reads_empty_string(self):
        reader = cp._Reader(io.BytesIO(struct.pack("<H", 0)))
        self.assertEqual(reader.string(), "")

    def test_u64_array_reads_zero_count_as_empty(self):
        reader = cp._Reader(io.BytesIO(b""))
        result = reader.u64_array(0)
        self.assertEqual(result.size, 0)
        self.assertEqual(result.dtype, np.uint64)

    def test_u64_array_reads_multiple_values_in_order(self):
        values = [0, 1, 2 ** 63, 2 ** 64 - 1]
        buffer = struct.pack("<4Q", *values)
        reader = cp._Reader(io.BytesIO(buffer))
        np.testing.assert_array_equal(reader.u64_array(4), np.array(values, dtype=np.uint64))

    def test_f64_array_reads_zero_count_as_empty(self):
        reader = cp._Reader(io.BytesIO(b""))
        result = reader.f64_array(0)
        self.assertEqual(result.size, 0)
        self.assertEqual(result.dtype, np.float64)

    def test_f64_array_reads_multiple_values_including_nan(self):
        values = [1.5, float("nan"), -2.25, 0.0]
        buffer = struct.pack("<4d", *values)
        reader = cp._Reader(io.BytesIO(buffer))
        result = reader.f64_array(4)
        np.testing.assert_array_equal(result[[0, 2, 3]], np.array([1.5, -2.25, 0.0]))
        self.assertTrue(np.isnan(result[1]))


# ---------------------------------------------------------------------------
# _load_request
# ---------------------------------------------------------------------------

def _pack_string(value: str) -> bytes:
    payload = value.encode("utf-8")
    return struct.pack("<H", len(payload)) + payload


def _pack_package(
        *,
        component_name="chr1",
        compo_length=100,
        query_start=0,
        query_end=99,
        full_component=True,
        crosses_origin=False,
        visible=True,
        logarithmic=False,
        log_base=0,
        y_upper_limit=1.0,
        local_x=(0, 1, 2),
        plot_y=(1.0, 2.0, 3.0),
        tick_raw_values=(0.0, 1.0),
        tick_display_values=(0.0, 1.0),
) -> bytes:
    """Builds the exact byte layout BinaryRequestWriter::package() produces (see the format
    docstring at the top of circular_plot.py), so _load_request can be exercised without any
    C++ involvement."""
    out = bytearray()
    out += _pack_string(component_name)
    out += struct.pack("<Q", compo_length)
    out += struct.pack("<Q", query_start)
    out += struct.pack("<Q", query_end)
    out += struct.pack("<B", 1 if full_component else 0)
    out += struct.pack("<B", 1 if crosses_origin else 0)
    out += struct.pack("<B", 1 if visible else 0)
    out += struct.pack("<B", 1 if logarithmic else 0)
    out += struct.pack("<i", log_base)
    out += struct.pack("<d", y_upper_limit)
    out += struct.pack("<Q", len(local_x))
    out += struct.pack(f"<{len(local_x)}Q", *local_x)
    out += struct.pack(f"<{len(plot_y)}d", *plot_y)
    out += struct.pack("<I", len(tick_raw_values))
    out += struct.pack(f"<{len(tick_raw_values)}d", *tick_raw_values)
    out += struct.pack(f"<{len(tick_display_values)}d", *tick_display_values)
    return bytes(out)


def _pack_header(*, mode, dpi=300, fig_width=7.0, fig_height=4.5,
                  line_color="#1E3A8A", fill_color="#93C5FD", version=2) -> bytes:
    out = bytearray()
    out += b"CXCP"
    out += struct.pack("<B", version)
    out += struct.pack("<B", mode)
    out += struct.pack("<I", dpi)
    out += struct.pack("<d", fig_width)
    out += struct.pack("<d", fig_height)
    out += _pack_string(line_color)
    out += _pack_string(fill_color)
    return bytes(out)


@_skip_unless_deps()
class LoadRequestTest(unittest.TestCase):
    """_load_request() against hand-built request buffers, written to real temp files (the
    function takes a path and opens the file itself, mirroring how src/python_circular_plot.cpp
    invokes it)."""

    def _write_and_load(self, tmp_path, payload: bytes) -> dict:
        tmp_path.write_bytes(payload)
        return cp._load_request(tmp_path)

    def test_query_mode_round_trips_all_fields(self, ):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "request.bin"
            payload = _pack_header(mode=0, dpi=150, fig_width=5.0, fig_height=3.0,
                                    line_color="#010203", fill_color="#040506")
            payload += _pack_package(
                component_name="chrX", compo_length=1000, query_start=10, query_end=20,
                full_component=False, crosses_origin=False, visible=True,
                logarithmic=True, log_base=2, y_upper_limit=42.5,
                local_x=(0, 5, 10), plot_y=(1.0, float("nan"), 3.0),
                tick_raw_values=(0.0, 10.0, 20.0), tick_display_values=(0.0, 1.0, 2.0),
            )
            request = self._write_and_load(path, payload)

        self.assertEqual(request["mode"], "query")
        self.assertEqual(request["dpi"], 150)
        self.assertEqual(request["fig_width"], 5.0)
        self.assertEqual(request["fig_height"], 3.0)
        self.assertEqual(request["line_color"], "#010203")
        self.assertEqual(request["fill_color"], "#040506")

        pkg = request["package"]
        self.assertEqual(pkg["component_name"], "chrX")
        self.assertEqual(pkg["compo_length"], 1000)
        self.assertEqual(pkg["query_start"], 10)
        self.assertEqual(pkg["query_end"], 20)
        self.assertFalse(pkg["full_component"])
        self.assertFalse(pkg["crosses_origin"])
        self.assertTrue(pkg["visible"])
        self.assertTrue(pkg["logarithmic"])
        self.assertEqual(pkg["log_base"], 2)
        self.assertEqual(pkg["y_upper_limit"], 42.5)
        np.testing.assert_array_equal(pkg["local_x"], np.array([0, 5, 10], dtype=np.uint64))
        np.testing.assert_array_equal(pkg["plot_y"][[0, 2]], np.array([1.0, 3.0]))
        self.assertTrue(np.isnan(pkg["plot_y"][1]))
        np.testing.assert_array_equal(pkg["tick_raw_values"], np.array([0.0, 10.0, 20.0]))
        np.testing.assert_array_equal(pkg["tick_display_values"], np.array([0.0, 1.0, 2.0]))

    def test_global_mode_round_trips_layout_and_multiple_packages(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "request.bin"
            payload = _pack_header(mode=1)
            payload += struct.pack("<d", 6.0)  # subplot_width
            payload += struct.pack("<d", 6.0)  # subplot_height
            payload += struct.pack("<I", 2)  # rows
            payload += struct.pack("<I", 3)  # columns
            payload += struct.pack("<I", 2)  # component_count
            payload += _pack_package(component_name="chr1", visible=True)
            payload += _pack_package(component_name="chr2", visible=False, local_x=(), plot_y=(),
                                      tick_raw_values=(), tick_display_values=())
            request = self._write_and_load(path, payload)

        self.assertEqual(request["mode"], "global")
        self.assertEqual(request["subplot_width"], 6.0)
        self.assertEqual(request["subplot_height"], 6.0)
        self.assertEqual(request["rows"], 2)
        self.assertEqual(request["columns"], 3)
        self.assertEqual(len(request["packages"]), 2)
        self.assertEqual(request["packages"][0]["component_name"], "chr1")
        self.assertTrue(request["packages"][0]["visible"])
        self.assertEqual(request["packages"][1]["component_name"], "chr2")
        self.assertFalse(request["packages"][1]["visible"])
        self.assertEqual(request["packages"][1]["local_x"].size, 0)

    def test_invalid_magic_raises_value_error(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "request.bin"
            payload = b"XXXX" + _pack_header(mode=0)[4:] + _pack_package()
            path.write_bytes(payload)
            with self.assertRaisesRegex(ValueError, "Invalid magic number"):
                cp._load_request(path)

    def test_unsupported_version_raises_value_error(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "request.bin"
            payload = _pack_header(mode=0, version=1) + _pack_package()
            path.write_bytes(payload)
            with self.assertRaisesRegex(ValueError, "Unsupported binary format version"):
                cp._load_request(path)

    def test_unknown_mode_byte_raises_value_error(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "request.bin"
            payload = _pack_header(mode=7) + _pack_package()
            path.write_bytes(payload)
            with self.assertRaisesRegex(ValueError, "Unknown mode byte"):
                cp._load_request(path)


# ---------------------------------------------------------------------------
# _plot_circular_graph: shared helpers
# ---------------------------------------------------------------------------

def _make_pkg(**overrides) -> dict:
    """A minimal-but-complete package dict, with every field _plot_circular_graph reads."""
    pkg = dict(
        component_name="chr1",
        compo_length=10,
        query_start=0,
        query_end=9,
        full_component=True,
        crosses_origin=False,
        visible=True,
        logarithmic=False,
        log_base=0,
        y_upper_limit=10.0,
        local_x=np.array([0, 1, 2, 3, 4], dtype=np.uint64),
        plot_y=np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
        tick_raw_values=np.array([0.0, 5.0, 10.0]),
        tick_display_values=np.array([0.0, 5.0, 10.0]),
    )
    pkg.update(overrides)
    return pkg


def _make_circos_mock():
    """Builds a Circos(...)/sector/track mock chain deep enough to satisfy every call
    _plot_circular_graph makes when it does attempt to render, without doing any real
    plotting work."""
    mock_track = mock.MagicMock(name="track")
    mock_sector = mock.MagicMock(name="sector")
    mock_sector.add_track.return_value = mock_track
    mock_circos_instance = mock.MagicMock(name="circos_instance")
    mock_circos_instance.sectors = [mock_sector]
    mock_circos_class = mock.MagicMock(name="Circos", return_value=mock_circos_instance)
    return mock_circos_class, mock_circos_instance, mock_sector, mock_track


@_skip_unless_deps()
class PlotCircularGraphEarlyExitTest(unittest.TestCase):
    """The two early-return branches of _plot_circular_graph: neither should attempt to
    construct a Circos instance at all, since there is nothing drawable."""

    def test_invisible_package_hides_axis_without_constructing_circos(self):
        pkg = _make_pkg(visible=False)
        mock_circos_class, *_ = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=True,
                                     line_color="#000000", fill_color="#ffffff")

        axis.set_visible.assert_called_once_with(False)
        mock_circos_class.assert_not_called()

    def test_empty_local_x_hides_axis_without_constructing_circos(self):
        pkg = _make_pkg(local_x=np.array([], dtype=np.uint64), plot_y=np.array([]))
        mock_circos_class, *_ = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=True,
                                     line_color="#000000", fill_color="#ffffff")

        axis.set_visible.assert_called_once_with(False)
        mock_circos_class.assert_not_called()

    def test_all_nan_plot_y_hides_axis_without_constructing_circos(self):
        pkg = _make_pkg(
            local_x=np.array([0, 1, 2], dtype=np.uint64),
            plot_y=np.array([float("nan"), float("nan"), float("nan")]),
        )
        mock_circos_class, *_ = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=True,
                                     line_color="#000000", fill_color="#ffffff")

        axis.set_visible.assert_called_once_with(False)
        mock_circos_class.assert_not_called()


@_skip_unless_deps()
class DrawSegmentSplittingTest(unittest.TestCase):
    """The NaN-gap segment-splitting logic (draw_segment, nested inside
    _plot_circular_graph), for all three traversal topologies. Expected coordinates were
    hand-derived by tracing circular_plot.py's actual algorithm; see the comment above each
    test for the derivation.
    """

    @staticmethod
    def _segments_drawn(track_mock) -> list:
        """Extracts the (x, y) arrays passed to every track.line(...) call, in call order.
        (track.fill_between is called with the same segment_x/segment_y once per segment
        too; asserting on .line is enough to pin down the splitting behavior, and a
        dedicated test below checks fill_between/line stay in lockstep.)"""
        return [(call.args[0], call.args[1]) for call in track_mock.line.call_args_list]

    def test_full_component_splits_on_two_internal_gaps(self):
        # compo_length=5, local_x=[0,1,2,3,4] (no downsampling), plot_y has NaNs at
        # traversal indices 1 and 4.
        #   display_x = local_x + 0, then append compo_end(5)        -> [0,1,2,3,4,5]
        #   display_y = plot_y,      then append plot_y[0]  (=1.0)   -> [1,nan,3,4,nan,1]
        #   finite_indices = [0,2,3,5]; diff = [2,1,2] -> splits at [1,3]
        #   blocks -> [0], [2,3], [5]
        pkg = _make_pkg(
            compo_length=5, query_start=0, query_end=4, full_component=True, crosses_origin=False,
            local_x=np.array([0, 1, 2, 3, 4], dtype=np.uint64),
            plot_y=np.array([1.0, float("nan"), 3.0, 4.0, float("nan")]),
        )
        mock_circos_class, _instance, _sector, mock_track = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=False,
                                     line_color="#000000", fill_color="#ffffff")

        segments = self._segments_drawn(mock_track)
        self.assertEqual(len(segments), 3)
        np.testing.assert_array_equal(segments[0][0], [0]);        np.testing.assert_array_equal(segments[0][1], [1.0])
        np.testing.assert_array_equal(segments[1][0], [2, 3]);     np.testing.assert_array_equal(segments[1][1], [3.0, 4.0])
        np.testing.assert_array_equal(segments[2][0], [5]);        np.testing.assert_array_equal(segments[2][1], [1.0])

    def test_subrange_query_splits_on_single_internal_gap(self):
        # query_start=50, query_end=55 (contiguous, not full, not crossing origin), a
        # downsampled traversal local_x=[0,2,4] with a NaN at the middle point.
        #   display_x = local_x + query_start, then append query_end+1(56) -> [50,52,54,56]
        #   display_y = plot_y,     then append plot_y[-1] (=30.0)         -> [10,nan,30,30]
        #   finite_indices = [0,2,3]; diff = [2,1] -> split at [1]
        #   blocks -> [0], [2,3]
        pkg = _make_pkg(
            compo_length=1000, query_start=50, query_end=55, full_component=False, crosses_origin=False,
            local_x=np.array([0, 2, 4], dtype=np.uint64),
            plot_y=np.array([10.0, float("nan"), 30.0]),
        )
        mock_circos_class, _instance, _sector, mock_track = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=True,
                                     line_color="#000000", fill_color="#ffffff")

        segments = self._segments_drawn(mock_track)
        self.assertEqual(len(segments), 2)
        np.testing.assert_array_equal(segments[0][0], [50]);       np.testing.assert_array_equal(segments[0][1], [10.0])
        np.testing.assert_array_equal(segments[1][0], [54, 56]);   np.testing.assert_array_equal(segments[1][1], [30.0, 30.0])

    def test_origin_crossing_query_splits_both_sides_and_interpolates_the_join(self):
        # compo_length=20, query_start=18, query_end=2 -> first_segment_length = 20-18 = 2.
        # local_x=[0,1,3,4]: traversal indices 0,1 (< 2) belong to the "first" side
        # (positions 18,19), indices 3,4 (>= 2) belong to the "second" side (positions 1,2).
        # plot_y=[10,20,nan,50] (a gap on the second side only).
        #
        # origin_value = np.interp(first_segment_length=2, local_x[finite]=[0,1,4],
        #                           plot_y[finite]=[10,20,50]) = 20 + (2-1)/(4-1)*(50-20) = 30.0
        # (chosen so this lands on an exact value, no floating-point fuzz needed).
        #
        # first_x  = [0,1]+18 = [18,19], then append compo_end(20)         -> [18,19,20]
        # first_y  = [10,20],            then append origin_value(30)     -> [10,20,30]
        #   -> no NaNs -> one segment: x=[18,19,20], y=[10,20,30]
        #
        # second_x = [3,4]-2+0 = [1,2], insert compo_start(0) at front     -> [0,1,2]
        # second_y = [nan,50],          insert origin_value(30) at front  -> [30,nan,50]
        #   then append query_end+1(3) to x, and current-last(50) to y:
        #            second_x -> [0,1,2,3]; second_y -> [30,nan,50,50]
        #   finite_indices = [0,2,3]; diff = [2,1] -> split at [1] -> blocks [0], [2,3]
        #   -> two segments: x=[0],y=[30]  and  x=[2,3],y=[50,50]
        pkg = _make_pkg(
            compo_length=20, query_start=18, query_end=2, full_component=False, crosses_origin=True,
            local_x=np.array([0, 1, 3, 4], dtype=np.uint64),
            plot_y=np.array([10.0, 20.0, float("nan"), 50.0]),
        )
        mock_circos_class, _instance, _sector, mock_track = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=True,
                                     line_color="#000000", fill_color="#ffffff")

        segments = self._segments_drawn(mock_track)
        self.assertEqual(len(segments), 3)
        np.testing.assert_array_equal(segments[0][0], [18, 19, 20]); np.testing.assert_array_equal(segments[0][1], [10.0, 20.0, 30.0])
        np.testing.assert_array_equal(segments[1][0], [0]);          np.testing.assert_array_equal(segments[1][1], [30.0])
        np.testing.assert_array_equal(segments[2][0], [2, 3]);       np.testing.assert_array_equal(segments[2][1], [50.0, 50.0])

    def test_fill_between_and_line_are_called_once_per_segment_in_lockstep(self):
        pkg = _make_pkg(
            compo_length=5, query_start=0, query_end=4, full_component=True, crosses_origin=False,
            local_x=np.array([0, 1, 2, 3, 4], dtype=np.uint64),
            plot_y=np.array([1.0, float("nan"), 3.0, 4.0, float("nan")]),
        )
        mock_circos_class, _instance, _sector, mock_track = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=False,
                                     line_color="#123456", fill_color="#abcdef")

        self.assertEqual(mock_track.fill_between.call_count, mock_track.line.call_count)
        self.assertEqual(mock_track.line.call_count, 3)
        # Spot-check that the configured colors and y-axis bounds reach the draw calls.
        _, kwargs = mock_track.line.call_args_list[0]
        self.assertEqual(kwargs["color"], "#123456")
        self.assertEqual(kwargs["vmax"], pkg["y_upper_limit"])
        _, kwargs = mock_track.fill_between.call_args_list[0]
        self.assertEqual(kwargs["color"], "#abcdef")

    def test_no_drawable_segments_when_query_mode_sub_range_is_all_nan_but_not_flagged_invisible(self):
        # Defensive test: even if upstream preprocessing somehow left visible=True with an
        # entirely-NaN plot_y that still has local_x entries, _plot_circular_graph's own
        # `np.all(np.isnan(plot_y))` guard (checked before touching Circos at all) must catch
        # it - this duplicates one PlotCircularGraphEarlyExitTest case but from the
        # draw-segment test's angle (asserting zero draw calls, not just axis.set_visible).
        pkg = _make_pkg(
            local_x=np.array([0, 1, 2], dtype=np.uint64),
            plot_y=np.array([float("nan")] * 3),
        )
        mock_circos_class, _instance, _sector, mock_track = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=True,
                                     line_color="#000000", fill_color="#ffffff")

        mock_track.line.assert_not_called()
        mock_track.fill_between.assert_not_called()


@_skip_unless_deps()
class RadialTicksTest(unittest.TestCase):
    """The radial coverage-tick loop: raw_tick == 0 is skipped on a log scale (not
    representable there) but kept on a linear scale. Uses axis.plot/axis.text call counts,
    since those two calls are made ONLY by this tick loop (track.xticks is a separate,
    unrelated call for the genomic-position ticks)."""

    def test_zero_tick_is_skipped_on_logarithmic_scale(self):
        pkg = _make_pkg(
            logarithmic=True, log_base=2,
            tick_raw_values=np.array([0.0, 10.0, 20.0]),
            tick_display_values=np.array([0.0, 1.0, 2.0]),
        )
        mock_circos_class, *_ = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=False,
                                     line_color="#000000", fill_color="#ffffff")

        # 3 ticks provided, 1 skipped (raw_tick == 0 under log scale) -> 2 draw calls each.
        self.assertEqual(axis.plot.call_count, 2)
        self.assertEqual(axis.text.call_count, 2)

    def test_zero_tick_is_kept_on_linear_scale(self):
        pkg = _make_pkg(
            logarithmic=False, log_base=0,
            tick_raw_values=np.array([0.0, 10.0, 20.0]),
            tick_display_values=np.array([0.0, 10.0, 20.0]),
        )
        mock_circos_class, *_ = _make_circos_mock()
        axis = mock.MagicMock(name="axis")

        with mock.patch.object(cp, "Circos", mock_circos_class):
            cp._plot_circular_graph(axis=axis, pkg=pkg, query_mode=False,
                                     line_color="#000000", fill_color="#ffffff")

        # All 3 ticks kept on a linear scale, including the zero one.
        self.assertEqual(axis.plot.call_count, 3)
        self.assertEqual(axis.text.call_count, 3)


if __name__ == "__main__":
    unittest.main()
