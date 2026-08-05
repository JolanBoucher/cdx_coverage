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
#   ./scripts/install_ubuntu20.sh
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j"$(nproc)"
#   ./build/cdx_coverage --help
#
# Safe to re-run: every step is idempotent (apt install of already-installed
# packages is a no-op, Abseil/CLI11 steps are skipped if already present).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CDX_COVERAGE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

ABSEIL_TAG="20240722.0"   # Kept in sync with .github/workflows/ubuntu.yml -
                          # this is the version CI actually tests, not just
                          # the ">=20230802" floor CMakeLists.txt requires.

echo "==> [1/5] Installing system packages (apt)"
sudo apt update
sudo apt install -y \
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
    echo "==> [2/5] Installing GCC 11 (Ubuntu 20.04's default GCC 9.4 is too"
    echo "    old for a robust C++17 baseline)"
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt update
    sudo apt install -y gcc-11 g++-11
    export CC=gcc-11
    export CXX=g++-11
else
    echo "==> [2/5] Ubuntu ${UBUNTU_VERSION_ID}: default compiler is new"
    echo "    enough, skipping the GCC-11 PPA step"
    export CC="${CC:-gcc}"
    export CXX="${CXX:-g++}"
fi
echo "    Using CC=$CC CXX=$CXX for the rest of this script."
echo "    Export these yourself before running cmake by hand later:"
echo "      export CC=$CC CXX=$CXX"

echo "==> [3/5] Installing CLI11 (header-only)"
cd "${CDX_COVERAGE_DIR}"
if [ -f include/CLI/CLI.hpp ]; then
    echo "    include/CLI/CLI.hpp already present, skipping"
elif sudo apt install -y libcli11-dev 2>/dev/null; then
    echo "    Installed libcli11-dev from apt"
else
    echo "    libcli11-dev not available on this Ubuntu release, fetching"
    echo "    the single-header release instead"
    mkdir -p include/CLI
    curl -L -o include/CLI/CLI.hpp \
        https://github.com/CLIUtils/CLI11/releases/latest/download/CLI11.hpp
fi

echo "==> [4/5] Building and installing Abseil ${ABSEIL_TAG} from source"
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
    sudo cmake --install "${ABSEIL_SRC_DIR}/build"
    sudo ldconfig
fi

echo "==> [5/5] Initializing cdx_coverage Git submodules (deps/libvgio +"
echo "    nested deps/libhandlegraph)"
cd "${CDX_COVERAGE_DIR}"
git submodule update --init --recursive

cat <<EOF

==> Done. Next steps:

    export CC=$CC CXX=$CXX
    cd "${CDX_COVERAGE_DIR}"
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j"\$(nproc)"
    ./build/cdx_coverage --help

EOF
