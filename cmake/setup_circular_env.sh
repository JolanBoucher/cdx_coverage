#!/usr/bin/env bash
# Best-effort provisioning of a private Python venv for circular graph
# rendering (numpy/matplotlib/pycirclize), run as a build step so end users
# never have to create/activate a venv by hand - see CIRCULAR_PLOT_SETUP.md.
#
# Invoked by CMake as a POST_BUILD step on the cdx_coverage target. Must
# NEVER fail the build: any problem here (offline machine, no `venv` module,
# etc.) is reported as a warning and the script still exits 0. Circular
# graphs simply won't work until dependencies are available (same as if this
# script didn't exist at all) - see resolvePythonExecutable() in
# src/python_circular_plot.cpp, which falls back to the system "python3" (or
# CDX_PYTHON_EXECUTABLE) when no bundled venv is found.
#
# Usage: setup_circular_env.sh <pyenv_dir> <requirements_file> <base_python3>
set -u

PYENV_DIR="$1"
REQUIREMENTS_FILE="$2"
BASE_PYTHON="${3:-python3}"

PYTHON_EXE="$PYENV_DIR/bin/python3"
LOG_PREFIX="[cdx_coverage]"

# Fast path: already provisioned and satisfies requirements.
if [ -x "$PYTHON_EXE" ] && "$PYTHON_EXE" -c "import numpy, matplotlib, pycirclize" >/dev/null 2>&1; then
    exit 0
fi

echo "$LOG_PREFIX Provisioning Python environment for circular graphs at: $PYENV_DIR"

if ! "$BASE_PYTHON" -m venv "$PYENV_DIR" >/tmp/cdx_pyenv_setup.log 2>&1; then
    echo "$LOG_PREFIX WARNING: could not create a Python venv at $PYENV_DIR." >&2
    cat /tmp/cdx_pyenv_setup.log >&2
    echo "$LOG_PREFIX Circular graphs (--component-type circular) will not work until you" >&2
    echo "$LOG_PREFIX install numpy/matplotlib/pycirclize yourself - see CIRCULAR_PLOT_SETUP.md." >&2
    exit 0
fi

"$PYTHON_EXE" -m pip install --quiet --upgrade pip >/tmp/cdx_pyenv_setup.log 2>&1 || true

if ! "$PYTHON_EXE" -m pip install --quiet -r "$REQUIREMENTS_FILE" >/tmp/cdx_pyenv_setup.log 2>&1; then
    echo "$LOG_PREFIX WARNING: could not install numpy/matplotlib/pycirclize into $PYENV_DIR." >&2
    echo "$LOG_PREFIX This is usually a network issue during the build (no internet access)." >&2
    cat /tmp/cdx_pyenv_setup.log >&2
    echo "$LOG_PREFIX Circular graphs will not work until dependencies are installed - see" >&2
    echo "$LOG_PREFIX CIRCULAR_PLOT_SETUP.md (or re-run this build with network access)." >&2
    exit 0
fi

echo "$LOG_PREFIX Python environment ready: $PYENV_DIR"
exit 0
