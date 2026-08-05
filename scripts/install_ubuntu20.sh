#!/usr/bin/env bash
#
# Bootstrap a fresh Ubuntu machine to build cdx_coverage, end to end.
# Ubuntu 20.04 is the primary target (hence the filename), but this script
# also runs cleanly on newer releases (e.g. 24.04): it detects the Ubuntu
# version and skips the GCC-11 PPA step where the default compiler is
# already new enough.
#
# .github/workflows/ubuntu.yml calls this script directly (for both its
# ubuntu-20.04 and ubuntu-24.04 matrix legs) instead of duplicating these
# steps, so this script and CI can never silently drift apart. If you need
# to change the dependency list or Abseil version, change it here only.
#
# Assumed layout when this script is run:
#   <parent>/cdx_coverage   <- this repo, submodules NOT required to be
#                              initialized yet (the script does it)
#
# cdx_lib is neither a submodule nor a sibling checkout: CMake's
# FetchContent fetches and builds it automatically the first time you run
# `cmake -S . -B build` (see CMakeLists.txt), so this script doesn't need to
# do anything for it - a working internet connection at configure time is
# all that's required.
#
# Usage:
#   cd cdx_coverage
#   ./scripts/install_ubuntu20.sh              # installs deps, configures,
#                                               # builds, and runs the test
#                                               # suite (ctest)
#   ./scripts/install_ubuntu20.sh --no-tests    # same, but skips building
#                                               # and running the tests
#
# Safe to re-run: every step is idempotent (apt install of already-installed
# packages is a no-op, Abseil/CLI11 steps are skipped if already present,
# cmake --build only rebuilds what changed).

set -euo pipefail

RUN_TESTS=1
for arg in "$@"; do
    case "$arg" in
        --no-tests | --no-test | --skip-tests)
            RUN_TESTS=0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Usage: $0 [--no-tests]" >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CDX_COVERAGE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${CDX_COVERAGE_DIR}/build"

ABSEIL_TAG="20240722.0"   # Kept in sync with .github/workflows/ubuntu.yml -
                          # this is the version CI actually tests, not just
                          # the ">=20230802" floor CMakeLists.txt requires.

# Root containers (common for Docker-based Ubuntu images/CI) typically don't
# have `sudo` installed at all - it would be a no-op anyway since root
# already has every permission `sudo` would grant. Detect that up front
# instead of hardcoding `sudo` in every privileged command below.
if [ "$(id -u)" -eq 0 ]; then
    SUDO=""
else
    SUDO="sudo"
fi

echo "==> [1/6] Installing system packages (apt)"
$SUDO apt update
$SUDO apt install -y \
    build-essential \
    cmake \
    git \
    pkg-config \
    protobuf-compiler \
    libprotobuf-dev \
    libhts-dev \
    libomp-dev \
    libcairo2-dev \
    libjansson-dev \
    zlib1g-dev \
    curl \
    software-properties-common \
    python3 \
    python3-venv \
    python3-pip
# NOTE: no python3-dev / python3.8-dev here on purpose. CMakeLists.txt only
# requests find_package(Python3 COMPONENTS Interpreter) - no Python.h is
# ever included, the circular-plot integration is subprocess-based, not
# embedded. python3-venv IS required though: cmake/setup_circular_env.sh
# runs `python3 -m venv` as a post-build step to provision numpy/matplotlib/
# pycirclize; without python3-venv that step fails silently and circular
# graphs won't render.

UBUNTU_VERSION_ID="$(. /etc/os-release && echo "${VERSION_ID:-unknown}")"

if [ "${UBUNTU_VERSION_ID}" = "20.04" ]; then
    echo "==> [2/6] Installing GCC 11 (Ubuntu 20.04's default GCC 9.4 is too"
    echo "    old for a robust C++17 baseline)"
    $SUDO add-apt-repository -y ppa:ubuntu-toolchain-r/test
    $SUDO apt update
    $SUDO apt install -y gcc-11 g++-11
    export CC=gcc-11
    export CXX=g++-11
else
    echo "==> [2/6] Ubuntu ${UBUNTU_VERSION_ID}: default compiler is new"
    echo "    enough, skipping the GCC-11 PPA step"
    export CC="${CC:-gcc}"
    export CXX="${CXX:-g++}"
fi
echo "    Using CC=$CC CXX=$CXX for the rest of this script."
echo "    Export these yourself before running cmake by hand later:"
echo "      export CC=$CC CXX=$CXX"

echo "==> [3/6] Installing CLI11 (header-only)"
cd "${CDX_COVERAGE_DIR}"
if [ -f include/CLI/CLI.hpp ]; then
    echo "    include/CLI/CLI.hpp already present, skipping"
elif $SUDO apt install -y libcli11-dev 2>/dev/null; then
    echo "    Installed libcli11-dev from apt"
else
    echo "    libcli11-dev not available on this Ubuntu release, fetching"
    echo "    the single-header release instead"
    mkdir -p include/CLI
    curl -L -o include/CLI/CLI.hpp \
        https://github.com/CLIUtils/CLI11/releases/latest/download/CLI11.hpp
fi

echo "==> [4/6] Building and installing Abseil ${ABSEIL_TAG} from source"
echo "    (Ubuntu's packaged libabsl-dev predates the Logging library"
echo "    symbols - absl::log_internal_message, absl::log_internal_check_op"
echo "    - that CMakeLists.txt links against)"
if ldconfig -p | grep -q libabsl_log_internal_message; then
    echo "    A matching Abseil already appears to be installed, skipping build"
else
    ABSEIL_SRC_DIR="$(mktemp -d)/abseil-cpp"
    git clone --depth 1 --branch "${ABSEIL_TAG}" \
        https://github.com/abseil/abseil-cpp.git "${ABSEIL_SRC_DIR}"
    cmake -S "${ABSEIL_SRC_DIR}" -B "${ABSEIL_SRC_DIR}/build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DCMAKE_CXX_STANDARD=17 \
        -DABSL_PROPAGATE_CXX_STD=ON \
        -DABSL_BUILD_TESTING=OFF
    cmake --build "${ABSEIL_SRC_DIR}/build" -j"$(nproc)"
    $SUDO cmake --install "${ABSEIL_SRC_DIR}/build"
    $SUDO ldconfig
fi

echo "==> [5/6] Initializing cdx_coverage Git submodules (deps/libvgio +"
echo "    nested deps/libhandlegraph)"
cd "${CDX_COVERAGE_DIR}"
git submodule update --init --recursive

echo "==> [6/6] Configuring, building"
if [ "${RUN_TESTS}" -eq 1 ]; then
    echo "    and testing cdx_coverage"
else
    echo "    cdx_coverage (--no-tests given: skipping the test suite)"
fi

CDX_CMAKE_CONFIGURE_ARGS=(-DCMAKE_BUILD_TYPE=Release)
if [ "${RUN_TESTS}" -eq 0 ]; then
    # Also skip fetching/building GoogleTest entirely (CDX_BUILD_TESTS, see
    # top-level CMakeLists.txt) - not just skip running ctest below - so this
    # path works offline and doesn't waste time building test binaries no
    # one is about to run.
    CDX_CMAKE_CONFIGURE_ARGS+=(-DCDX_BUILD_TESTS=OFF)
fi

cmake -S "${CDX_COVERAGE_DIR}" -B "${BUILD_DIR}" "${CDX_CMAKE_CONFIGURE_ARGS[@]}"
cmake --build "${BUILD_DIR}" -j"$(nproc)"

if [ "${RUN_TESTS}" -eq 1 ]; then
    echo "==> Running test suite (ctest)"
    # `ctest --test-dir` needs CMake 3.20+; Ubuntu 20.04's apt cmake is
    # 3.16.3, where that flag is silently ignored (ctest then looks for
    # CTestTestfile.cmake in the current directory instead of in
    # ${BUILD_DIR}, finds none, and reports "No tests were found!!!" without
    # error). `cd` into the build dir instead - works on every CMake version.
    (cd "${BUILD_DIR}" && ctest --output-on-failure)
fi

cat <<EOF

==> Done.

    "${BUILD_DIR}/cdx_coverage" --help

EOF
