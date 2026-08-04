# Circular graph runtime dependencies

`cdx_coverage` renders circular coverage graphs (`--component-type circular`)
by invoking `python3` as a subprocess (via `popen`, in
`src/python_circular_plot.cpp`) on a small standalone script,
`python_script/circular_plot.py`, which uses
[pycirclize](https://github.com/moshi4/pyCirclize) on top of matplotlib.
This is a **runtime-only** dependency: no Python is required to compile
`cdx_coverage`, and linear graphs (`--component-type linear`) never touch
Python at all (pure C++/Cairo, in-process).

`circular_plot.py` only handles the pycirclize/matplotlib rendering itself:
all numeric preparation (downsampling, smoothing, log transform, tick math)
runs in C++ (`src/python_circular_plot.cpp`), so the Python side stays thin
and the CPU-bound work happens at C++ speed.

## Install: automatic (default)

**Building `cdx_coverage` normally is enough** - you don't need to create or
activate a venv yourself. If a Python 3 interpreter is found at CMake
configure time, a `POST_BUILD` step (`cmake/setup_circular_env.sh`)
provisions a private venv at `pyenv/` next to the built executable and
installs `numpy`/`matplotlib`/`pycirclize` (see `python_script/requirements.txt`)
into it automatically. `circular_plot.py` itself is copied next to the
executable the same way. `cdx_coverage` prefers this bundled venv at
runtime automatically (see "How the interpreter is resolved" below) - no
environment variable, no manual step.

This step requires network access at build time (to `pip install`). It is
**best-effort and never fails the build**: on an offline machine, or if
`python3 -m venv` isn't available, you'll see a warning during the build and
circular graphs simply won't work until dependencies are installed some
other way (see "Install: manual" below). Everything else (compiling, linear
graphs) is unaffected.

To force re-provisioning (e.g. after editing `requirements.txt`), delete the
`pyenv/` directory next to the executable and rebuild.

## Install: manual (offline machines, custom setups, CI)

```bash
python3 -m pip install numpy matplotlib pycirclize
```

Any reasonably recent Python 3 (3.9+) works. There is no pinned version
requirement. Verify with:

```bash
python3 -c "import numpy, matplotlib, pycirclize; print('ok')"
```

If this `python3` isn't the one `cdx_coverage` would otherwise pick, point
it there explicitly with `CDX_PYTHON_EXECUTABLE` (see below).

## How the interpreter is resolved

`resolvePythonExecutable()` (in `src/python_circular_plot.cpp`) checks, in
order:

1. The `CDX_PYTHON_EXECUTABLE` environment variable, if set (explicit
   override - point it at any interpreter binary).
2. The bundled venv CMake provisions automatically (`pyenv/bin/python3`
   next to the executable, see above) - this is what makes circular graphs
   work out of the box after a normal build, no manual step.
3. Plain `python3` resolved from the current process's `PATH`, as a last
   resort (e.g. if provisioning was skipped or failed).

## Running from an IDE (CLion, etc.)

If step 2 above didn't run (e.g. you're using an older build directory from
before this feature existed, or provisioning failed - see the build log),
`cdx_coverage` falls back to whatever `python3` is in `PATH`. That's the
**process's own** `PATH`, not a login shell's: if you activate a virtualenv
only through your shell profile (`~/.zshrc`, `~/.bashrc`) or
`source venv/bin/activate` manually, that's invisible to a binary launched
directly by an IDE - it runs fine from a terminal but fails with
`ModuleNotFoundError: No module named 'pycirclize'` when run from CLion (or
similar).

Fix: either rebuild so the bundled venv gets provisioned (simplest), or set
`CDX_PYTHON_EXECUTABLE` explicitly in the IDE's run configuration:

```
CDX_PYTHON_EXECUTABLE=/path/to/venv/bin/python3
```

In CLion: *Run > Edit Configurations... > Environment variables* for the
`cdx_coverage` target. This bypasses both the bundled venv and `PATH`
resolution, and always uses that exact interpreter.

## How the script is located at runtime

`resolveCircularPlotScript()` (in `src/python_circular_plot.cpp`) checks, in
order:

1. The `CDX_CIRCULAR_PLOT_SCRIPT` environment variable, if set (explicit
   override — point it at any path).
2. `circular_plot.py` next to the running executable (where CMake places
   it, both in the build directory and after `install`).
3. `python_script/circular_plot.py` in the source tree (dev convenience,
   compiled in via `CDX_CIRCULAR_PLOT_SCRIPT_SOURCE_DIR`, for running the
   binary directly out of a build directory before any install step).

If none of these resolve, `cdx_coverage` raises a clear error naming the
environment variable as the fastest fix.

## Troubleshooting

- **`Failed to launch python3 ...`**: `python3` is not in the `PATH` of the
  shell/process running `cdx_coverage`. Install Python 3 or adjust `PATH`.
- **`Circular plot rendering failed (exit status ...)`**: the captured
  Python traceback is included directly in the exception message — the
  usual cause is a missing package (`numpy`, `matplotlib`, or `pycirclize`)
  or a coverage array that violates an invariant (e.g. mismatched lengths).
- **`Cannot locate circular_plot.py`**: the script wasn't found next to the
  executable and no source-tree fallback applies (e.g. you moved the binary
  without also copying the script). Set `CDX_CIRCULAR_PLOT_SCRIPT` to its
  full path.

## Why Python only for circular graphs, not linear?

Circos (the previous external tool for circular rendering) required a
notoriously fragile install (Perl GD/libgd with PNG support compiled in,
plus a separate `librsvg` conversion step — see the retired
`CIRCOS_SETUP.md`) and, after many rounds of configuration fixes, still had
persistent radial axis/gridline misalignment. pycirclize is a small,
well-behaved Python library that already renders exactly the intended look
(the original Python prototype used it), so a lightweight subprocess call
is both simpler and more reliable than continuing to hand-roll Circos-style
circular rendering in raw Cairo. Linear graphs don't have this problem —
Cairo renders them directly, in-process, with no subprocess at all.
