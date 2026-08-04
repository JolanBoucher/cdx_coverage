# cdx_coverage

A high-performance C++ tool for computing per-base sequencing coverage over
a linearized pangenome graph, from a custom binary graph index (**CDX**)
and read alignments in VG's **GAM** format.

`cdx_coverage` is designed for large pangenome graphs: coverage extraction,
projection to genomic coordinates, and coverage-graph rendering all run
in-process and are multi-threaded, with no dependency on the graph
toolkit's own `vg` binary at runtime.

## Overview

Given a CDX index and a GAM alignment file built from the *same* pangenome
graph, `cdx_coverage` computes, for every graph component (chromosome/
contig) or a user-selected sub-region:

1. per-base coverage depth, projected from node-level alignment counts to
   linear genomic coordinates;
2. summary statistics (mapping rate, breadth/depth of coverage, coverage
   distribution quantiles); and
3. a coverage plot, rendered either as a conventional linear track or as a
   circular (Circos-style) genome plot.

> **Important:** the GAM file must have been produced from alignments to
> the exact same pangenome graph used to build the CDX index — the node
> identifiers referenced in the GAM must match those stored in the index.
> `cdx_coverage` does not verify graph identity beyond structural checks.

## Features

- **Two query modes**: whole-graph coverage (all components at once) or a
  single component/sub-region (`-q/--query`), by name or numeric ID.
- **Two plot styles**: linear tracks (`--component-type linear`, default)
  or circular genome plots (`--component-type circular`), each with a
  single-panel view for one component/region and a multi-panel grid view
  for the whole graph.
- **Linear, in-process rendering.** Linear plots are rasterized directly
  via Cairo in the same process — no subprocess, no external interpreter.
- **Multi-core throughout**: GAM decompression, coordinate projection, and
  multi-panel plot rendering (both linear and circular) are parallelized.
- **Linear or logarithmic coverage scale** (`--log`), configurable
  smoothing, downsampling, resolution, figure size, and colors.
- **CDX inspection mode** (`-i/--inspect`) to list or describe graph
  components without processing any alignments.

## Requirements

### Build-time dependencies

- A C++17 compiler (GCC or Clang)
- [CMake](https://cmake.org/) ≥ 3.16
- [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
- [Protobuf](https://github.com/protocolbuffers/protobuf)
- [HTSlib](https://github.com/samtools/htslib) ≥ 1.10
- [Abseil](https://abseil.io/) (required by modern Protobuf)
- OpenMP
- POSIX threads
- [Cairo](https://www.cairographics.org/) — in-process linear plot rendering
- [CLI11](https://github.com/CLIUtils/CLI11) — command-line parsing
  (single-header library; see [Installing dependencies](#installing-dependencies))
- [libvgio](https://github.com/vgteam/libvgio) — built from the
  `deps/libvgio` git submodule; deserializes GAM files and exposes the VG
  Protobuf alignment types (`vg::Alignment`, `vg::Mapping`, `vg::Position`,
  `vg::Edit`) that `cdx_coverage` reads
- **`cdx_lib`** — a sibling library providing the CDX format types and I/O.
  `CMakeLists.txt` expects it at `../cdx_lib` relative to this repository
  (i.e. checked out next to `cdx_coverage`, not inside it)

### Runtime-only dependency (circular graphs)

Circular plots are rendered by a small Python script invoked as a
subprocess (all coverage/smoothing/statistics computation still happens in
C++; Python only handles the polar plot drawing via
[pycirclize](https://github.com/moshi4/pyCirclize)). This is **not**
required to build or run `cdx_coverage` with linear plots.

By default, the build automatically provisions a private Python
environment for this — see **[CIRCULAR_PLOT_SETUP.md](CIRCULAR_PLOT_SETUP.md)**
for details, manual setup, and troubleshooting.

## Installing dependencies

### Clone with submodules

```bash
git clone --recurse-submodules <repository-url> cdx_coverage
cd cdx_coverage
```

If already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

`cdx_lib` is **not** a submodule of this repository: clone it separately as
a sibling directory (`../cdx_lib` relative to `cdx_coverage`).

### macOS (Apple Silicon or Intel, via Homebrew)

```bash
brew install cmake pkg-config protobuf abseil htslib libomp cairo cli11
```

### Ubuntu 22.04 or newer

```bash
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    git \
    pkg-config \
    protobuf-compiler \
    libprotobuf-dev \
    libabsl-dev \
    libhts-dev \
    libomp-dev \
    libcairo2-dev \
    zlib1g-dev
```

Package availability for CLI11 and Abseil varies by Ubuntu release. If
`libcli11-dev` isn't available in your repositories, install CLI11 as a
single header instead:

```bash
mkdir -p include/CLI
curl -L -o include/CLI/CLI.hpp \
    https://github.com/CLIUtils/CLI11/releases/latest/download/CLI11.hpp
```

## Building

```bash
cmake -B build -S .
cmake --build build -j
```

The compiled binary is `build/cdx_coverage`. If a Python 3 interpreter was
found during configuration, a private virtual environment for circular
graph rendering is also provisioned automatically next to it (best-effort,
never fails the build — see
[CIRCULAR_PLOT_SETUP.md](CIRCULAR_PLOT_SETUP.md)).

To install system-wide:

```bash
cmake --install build
```

## Usage

```
cdx_coverage <CDX> [GAM] [OPTIONS]
```

### Inspect a CDX index

List all components in the index, without processing any alignments:

```bash
cdx_coverage graph.cdx --inspect
```

Describe a single component:

```bash
cdx_coverage graph.cdx --inspect chr1
```

### Whole-graph coverage

Compute coverage for every component and render a multi-panel grid plot:

```bash
cdx_coverage graph.cdx reads.gam -o results/
```

### Single component or sub-region

```bash
# Entire component, by name or numeric ID
cdx_coverage graph.cdx reads.gam -q chr1 -o results/
cdx_coverage graph.cdx reads.gam -q 0    -o results/

# 0-based sub-region
cdx_coverage graph.cdx reads.gam -q "chr1 1000-5000" -o results/
```

### Circular plots

```bash
cdx_coverage graph.cdx reads.gam -q chr1 -c circular -o results/
```

### Logarithmic scale, custom styling

```bash
cdx_coverage graph.cdx reads.gam \
    -q chr1 --log 10 \
    --color-line "#1E3A8A" --color-filling "#93C5FD" \
    --fig-size 7x4.5 --dpi 300 \
    -o results/
```

## Command-line reference

| Option | Description |
|---|---|
| `<CDX>` | Path to the binary CDX graph index (required). |
| `[GAM]` | Path to the GAM alignment file (required unless `-i/--inspect`). |
| `-q, --query TEXT` | Scope the computation to one component: `COMPONENT` or `"COMPONENT START-END"` (0-based), by name or numeric ID. Omit to process the whole graph. |
| `-c, --component-type` | `linear`/`l` (default) or `circular`/`c`. |
| `-i, --inspect [COMPONENT]` | Print CDX index contents and exit; no value lists all components. |
| `-o, --output PATH` | Output directory. Default: `.` |
| `--no-graph` | Skip coverage graph generation. |
| `--no-stats` | Skip the statistics report. |
| `--no-table` | Skip the per-base TSV table. |
| `--log [BASE]` | Logarithmic coverage scale; base defaults to 10. |
| `--smoothing FLOAT` | Moving-average window, as a fraction of length, in `[0.0, 1.0]`. Default: `0.01`. |
| `--max-point, --max-points N` | Maximum points passed to the plotting backend; `0` disables downsampling (full resolution). Default: `10000`. |
| `--dpi N` | Output graph resolution. Default: `300`. |
| `--fig-size WIDTHxHEIGHT` | Figure size in inches, e.g. `7x4.5`. |
| `--color-line HEX` | Coverage line color. Default: `#1E3A8A`. |
| `--color-filling HEX` | Coverage area fill color. Default: `#93C5FD`. |
| `-t, --worker-threads N` | Threads used for computation. Default: `auto` (all cores). |
| `-T, --decompression-threads N` | Threads used for GAM decompression. Default: `auto` (half the cores). |

Run `cdx_coverage --help` for the authoritative, up-to-date list.

## Output files

Written to the directory given by `-o/--output` (unless disabled via
`--no-table`/`--no-stats`/`--no-graph`):

| File | Description |
|---|---|
| `coverage_profile.tsv` | Per-base coverage table: `component_name`, `position`, `coverage`. |
| `coverage_stats.txt` | Mapping statistics (total/mapped/unmapped reads) and coverage statistics (breadth, mean, median, standard deviation, quartiles, min/max) per component or for the queried region. |
| `coverage_graph.png` | The rendered coverage plot (linear or circular, single panel or multi-panel grid, depending on the query and `-c/--component-type`). |

## Project layout

```
src/                 C++ source (CLI, CDX/GAM I/O, coverage projection,
                      statistics, linear plotting, circular plotting driver)
python_script/       circular_plot.py — pycirclize rendering backend for
                      circular graphs, invoked as a subprocess
cmake/                build-time helper scripts (Python env provisioning)
deps/libvgio/         libvgio git submodule
```

## License

No license has been specified yet for this project.
