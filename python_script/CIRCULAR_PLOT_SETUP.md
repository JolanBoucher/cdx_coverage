# Circular graph runtime setup

`cdx_coverage` renders circular coverage graphs (`--component-type circular`
/ `-c circular`) by invoking `python3` as a subprocess on a small,
standalone rendering script, `python_script/circular_plot.py`, which draws
Circos-style polar plots via [pycirclize](https://github.com/moshi4/pyCirclize)
on top of matplotlib.

This is the **only** part of `cdx_coverage` that touches Python, and only
at runtime:

- All coverage extraction, downsampling, smoothing, log-scale transforms,
  and tick math run in C++ (`src/python_circular_plot.cpp`), reusing the
  same numeric core as the linear (Cairo) backend. `circular_plot.py`
  receives already-computed points and is responsible only for drawing them.
- Linear graphs (`--component-type linear`, the default) never invoke
  Python at all — they are rendered in-process via Cairo
  (`src/cairo_plot.cpp`).
- No Python is required to *compile* `cdx_coverage`.

## Automatic setup (default)

Building `cdx_coverage` normally is enough — there is no manual venv to
create or activate. During the build:

1. `circular_plot.py` is copied next to the built executable.
2. If a Python 3 interpreter was found when CMake was configured, a
   `cmake/setup_circular_env.sh` build step provisions a private virtual
   environment (`pyenv/`, next to the executable) and installs
   `numpy`, `matplotlib`, and `pycirclize` into it
   (`python_script/requirements.txt`).

At runtime, `cdx_coverage` prefers this bundled environment automatically
(see [Interpreter resolution](#interpreter-resolution) below) — nothing
further to configure.

This step requires network access at build time (to `pip install`) and is
**best-effort**: it never fails the build. On an offline machine, or if the
`venv` module isn't available, you'll see a warning in the build log and
circular graphs won't work until dependencies are made available some
other way (see [Manual setup](#manual-setup)); everything else (compiling,
linear graphs) is unaffected.

To force re-provisioning (e.g. after editing `python_script/requirements.txt`),
delete the `pyenv/` directory next to the executable and rebuild.

## Manual setup

Useful for offline builds, custom environments, or CI:

```bash
python3 -m pip install numpy matplotlib pycirclize
```

Python 3.9+ is required; no other version is pinned. Verify with:

```bash
python3 -c "import numpy, matplotlib, pycirclize; print('ok')"
```

If this isn't the interpreter `cdx_coverage` would otherwise find, point it
there explicitly with `CDX_PYTHON_EXECUTABLE` (below).

## Interpreter resolution

`cdx_coverage` picks which `python3` to invoke in this order:

1. **`CDX_PYTHON_EXECUTABLE`** environment variable, if set — an explicit
   override pointing at any interpreter binary.
2. The bundled venv provisioned at build time (`pyenv/bin/python3`, next to
   the executable).
3. Plain `python3` resolved from the current process's `PATH`, as a last
   resort.

### Running from an IDE (CLion, etc.)

If step 2 didn't run (e.g. an older build directory, or provisioning
failed — check the build log) `cdx_coverage` falls back to whatever
`python3` is in `PATH`. That's the **process's own** `PATH`, not a login
shell's: a virtualenv activated only via a shell profile
(`~/.zshrc`, `~/.bashrc`) or `source venv/bin/activate` is invisible to a
binary launched directly by an IDE — it runs fine from a terminal but fails
with `ModuleNotFoundError: No module named 'pycirclize'` when run from
CLion (or similar).

Fix: rebuild so the bundled venv gets provisioned (simplest), or set
`CDX_PYTHON_EXECUTABLE` explicitly in the IDE's run configuration:

```
CDX_PYTHON_EXECUTABLE=/path/to/venv/bin/python3
```

In CLion: *Run > Edit Configurations… > Environment variables* for the
`cdx_coverage` target.

## Script resolution

`cdx_coverage` locates `circular_plot.py` in this order:

1. **`CDX_CIRCULAR_PLOT_SCRIPT`** environment variable, if set — an
   explicit override pointing at any path.
2. `circular_plot.py` next to the running executable (where CMake places
   it, both in the build directory and after `cmake --install`).
3. `python_script/circular_plot.py` in the source tree — a development
   convenience for running the binary directly from a build directory
   before any install step.

If none of these resolve, `cdx_coverage` raises an error naming
`CDX_CIRCULAR_PLOT_SCRIPT` as the fastest fix.

## Performance notes

For a whole-graph view (`--component-type circular` without `-q/--query`),
each component's panel is rendered independently and — since a single
Python process can't parallelize CPU-bound matplotlib work across cores
due to the GIL — in its own worker process (`multiprocessing.Pool`, sized
to the number of available cores), then composited into the final grid.
For a graph with many components, this is substantially faster than
rendering panels one at a time.

To see where time is spent on a specific run, set
`CDX_CIRCULAR_PLOT_TIMING=1`. `circular_plot.py` will print a per-stage
timing breakdown to stderr (interpreter startup/imports, request loading,
per-panel drawing, grid compositing, PNG export), forwarded by
`cdx_coverage` even on success (normally only shown on failure):

```bash
CDX_CIRCULAR_PLOT_TIMING=1 cdx_coverage graph.cdx reads.gam -c circular -o results/
```

Temporary working files (the binary request passed to Python, and — in
whole-graph mode — the individual panel PNGs before compositing) are
written next to the output file, under `<name>_circular_tmp/`, and removed
automatically once the final PNG has been written successfully.

## Troubleshooting

- **`Failed to launch 'python3' ...`**: no working `python3` was found by
  any of the steps in [Interpreter resolution](#interpreter-resolution).
  Install Python 3, or set `CDX_PYTHON_EXECUTABLE`.
- **`Circular plot rendering failed (exit status ...)`**: the captured
  Python traceback is included directly in the error message. The usual
  cause is a missing package (`numpy`, `matplotlib`, `pycirclize`) in
  whichever interpreter was resolved, or a coverage array that violates an
  internal invariant (e.g. mismatched lengths) — please report the latter
  as a bug.
- **`Cannot locate circular_plot.py`**: the script wasn't found next to the
  executable and no source-tree fallback applies (e.g. the binary was
  moved without copying the script alongside it). Set
  `CDX_CIRCULAR_PLOT_SCRIPT` to its full path.
