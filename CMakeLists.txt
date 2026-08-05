cmake_minimum_required(VERSION 3.16)

# Fix policy warnings BEFORE project()
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

project(
        cdx_coverage
        VERSION 0.1.0
        DESCRIPTION "Compute coverage from a pangenome using a CDX index and GAM alignments"
        LANGUAGES CXX
)

# ---------------------------------------------------------------------------
# Project configuration
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---------------------------------------------------------------------------
# Verify libvgio system dependencies
# ---------------------------------------------------------------------------
find_package(PkgConfig REQUIRED)
find_package(Protobuf REQUIRED)
find_package(Threads REQUIRED)

# Find Abseil (required by modern Protobuf on Homebrew / macOS)
find_package(absl REQUIRED)

# Homebrew installs libomp as a keg-only package on Apple Silicon / Intel macOS
if(APPLE)
    if(EXISTS "/opt/homebrew/opt/libomp")
        set(LIBOMP_PREFIX "/opt/homebrew/opt/libomp")
    elseif(EXISTS "/usr/local/opt/libomp")
        set(LIBOMP_PREFIX "/usr/local/opt/libomp")
    endif()

    if(LIBOMP_PREFIX)
        set(OpenMP_ROOT "${LIBOMP_PREFIX}")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
    endif()
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)

# libvgio requires HTSlib 1.10 or newer.
pkg_check_modules(
        HTSLIB
        REQUIRED
        IMPORTED_TARGET
        htslib>=1.10
)

# Cairo is used for in-process linear coverage plot rendering (no subprocess).
pkg_check_modules(
        CAIRO
        REQUIRED
        IMPORTED_TARGET
        cairo
)

# Circular coverage plots are rendered by a Python (pycirclize) subprocess
# (src/python_circular_plot.cpp + python_script/circular_plot.py) - see
# CIRCULAR_PLOT_SETUP.md. Python is only needed at RUNTIME for that one
# feature (--component-type circular); this project must still compile with
# no Python present. This lookup is purely informational (a configure-time
# warning), never REQUIRED. If found, a POST_BUILD step further below
# auto-provisions a private venv with numpy/matplotlib/pycirclize next to
# the built executable, so end users don't have to set one up by hand.
find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_FOUND)
    message(STATUS
            "Python 3 interpreter not found at configure time. This is fine: "
            "cdx_coverage compiles without it. It is only needed at RUNTIME to "
            "render circular graphs (--component-type circular). A Python 3 "
            "install would also let this build auto-provision numpy/matplotlib/"
            "pycirclize for you - see CIRCULAR_PLOT_SETUP.md."
    )
endif()

# ---------------------------------------------------------------------------
# Build libvgio as part of this project
# ---------------------------------------------------------------------------
set(LIBVGIO_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/libvgio")

if(NOT EXISTS "${LIBVGIO_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
            "libvgio was not found at:\n"
            "  ${LIBVGIO_SOURCE_DIR}\n\n"
            "Initialize the Git submodules with:\n"
            "  git submodule update --init --recursive"
    )
endif()

# Silence policy warning for CMP0074 in subdirectories
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

# Prevent libvgio from invoking find_package(Protobuf) a second time
set(CMAKE_DISABLE_FIND_PACKAGE_Protobuf TRUE)

add_subdirectory(
        "${LIBVGIO_SOURCE_DIR}"
        "${CMAKE_CURRENT_BINARY_DIR}/deps/libvgio"
)

# Re-enable find_package(Protobuf) for the rest of your project
set(CMAKE_DISABLE_FIND_PACKAGE_Protobuf FALSE)

# ---------------------------------------------------------------------------
# Normalize the libvgio target name & link Abseil
# ---------------------------------------------------------------------------
if(TARGET VGio::VGio)
    set(CDX_VGIO_TARGET VGio::VGio)
elseif(TARGET vgio)
    set(CDX_VGIO_TARGET vgio)
elseif(TARGET libvgio)
    set(CDX_VGIO_TARGET libvgio)
else()
    message(FATAL_ERROR "libvgio was configured, but no recognized CMake target was found.")
endif()

message(STATUS "Using libvgio target: ${CDX_VGIO_TARGET}")

# ---------------------------------------------------------------------------
# Portability fix: make sure "vg/vg.pb.h" resolves no matter which Protobuf
# CMake integration the host machine's libvgio build picked (see
# cmake/fix_libvgio_proto_symlink.cmake for the full explanation). This is
# self-healing on every machine/OS - no manual submodule patching needed.
# ---------------------------------------------------------------------------
set(LIBVGIO_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/deps/libvgio")

add_custom_command(
        OUTPUT "${LIBVGIO_BINARY_DIR}/.vg_proto_fixup_stamp"
        COMMAND "${CMAKE_COMMAND}"
                -DLIBVGIO_BINARY_DIR=${LIBVGIO_BINARY_DIR}
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/fix_libvgio_proto_symlink.cmake"
        COMMAND "${CMAKE_COMMAND}" -E touch "${LIBVGIO_BINARY_DIR}/.vg_proto_fixup_stamp"
        DEPENDS run_protoc
        COMMENT "Ensuring vg/vg.pb.h is reachable regardless of the host's Protobuf CMake integration"
        VERBATIM
)
add_custom_target(cdx_fix_libvgio_proto_include ALL
        DEPENDS "${LIBVGIO_BINARY_DIR}/.vg_proto_fixup_stamp")
# Must run after libvgio's own "vg -> ." symlink target so ours wins in all cases.
add_dependencies(cdx_fix_libvgio_proto_include run_protoc link_target)
add_dependencies(${CDX_VGIO_TARGET} cdx_fix_libvgio_proto_include)

# Link Abseil components required by Protobuf into vgio
target_link_libraries(${CDX_VGIO_TARGET} PRIVATE
        absl::log_internal_message
        absl::log_internal_check_op
        absl::hash
        absl::status
        absl::raw_logging_internal
)

# ---------------------------------------------------------------------------
# cdx_coverage executable
# ---------------------------------------------------------------------------
add_executable(cdx_coverage
        src/main.cpp
        src/gam_io.cpp
        src/gam_io.h
        src/config.h
        src/cdx_loader.cpp
        src/cdx_loader.h
        src/cov_projection.cpp
        src/cov_projection.h
        src/output_stats.cpp
        src/output_stats.h
        src/output_coverage.cpp
        src/output_coverage.h
        src/cli.cpp
        src/cli.hpp
        src/query_resolver.cpp
        src/query_resolver.h
        src/linear_plot.cpp
        src/output_plot.h
        src/python_circular_plot.cpp
        src/cairo_plot.cpp
)

target_compile_features(cdx_coverage PRIVATE cxx_std_17)

target_include_directories(cdx_coverage PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Dev-time fallback so python_circular_plot.cpp can find circular_plot.py
# straight from the source tree when running the binary directly out of a
# build directory (before any install step ran). At runtime, the script is
# looked up next to the executable first (see the copy/install rules below);
# this compile definition is only consulted if that lookup fails.
target_compile_definitions(cdx_coverage PRIVATE
        CDX_CIRCULAR_PLOT_SCRIPT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/python_script"
)

add_subdirectory(
        ../cdx_lib
        cdx_lib_build
)

target_link_libraries(cdx_coverage PRIVATE
        cdx_lib
        ${CDX_VGIO_TARGET}
        protobuf::libprotobuf
        absl::log_internal_message
        absl::log_internal_check_op
        absl::hash
        absl::status
        absl::raw_logging_internal
        PkgConfig::CAIRO
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(cdx_coverage PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()

# Copy circular_plot.py next to the built executable so it is found there at
# runtime (see resolveCircularPlotScript in src/python_circular_plot.cpp)
# without requiring an install step during development.
#
# NOT a `TARGET cdx_coverage POST_BUILD` command: those only fire when the
# cdx_coverage TARGET itself gets rebuilt/relinked. Editing only
# circular_plot.py (no C++ change) would then leave a stale copy next to the
# executable - `cmake --build` sees cdx_coverage is already up to date and
# skips it entirely, silently running old Python logic. Tying this to its
# own OUTPUT/DEPENDS instead makes `cmake --build` re-copy it on its own
# whenever the source .py file's content changes, independent of the C++
# target's build state.
set(CDX_CIRCULAR_PLOT_SCRIPT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/python_script/circular_plot.py")
set(CDX_CIRCULAR_PLOT_SCRIPT_DEST "${CMAKE_CURRENT_BINARY_DIR}/circular_plot.py")

add_custom_command(
        OUTPUT "${CDX_CIRCULAR_PLOT_SCRIPT_DEST}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CDX_CIRCULAR_PLOT_SCRIPT_SRC}"
        "${CDX_CIRCULAR_PLOT_SCRIPT_DEST}"
        DEPENDS "${CDX_CIRCULAR_PLOT_SCRIPT_SRC}"
        COMMENT "Copying circular_plot.py next to cdx_coverage"
)
add_custom_target(cdx_coverage_circular_plot_script ALL
        DEPENDS "${CDX_CIRCULAR_PLOT_SCRIPT_DEST}"
)
add_dependencies(cdx_coverage cdx_coverage_circular_plot_script)

# Auto-provision a private venv (numpy/matplotlib/pycirclize) next to the
# built executable, so building cdx_coverage on a fresh machine is enough to
# get working circular graphs too - no one has to create/activate a venv by
# hand. Only attempted if a Python 3 interpreter was found at configure time.
# NEVER fails the build (see cmake/setup_circular_env.sh: any problem, e.g.
# no network access, is a warning, not a build error). At runtime,
# resolvePythonExecutable() (src/python_circular_plot.cpp) prefers this
# bundled venv over the system "python3" when present, and
# CDX_PYTHON_EXECUTABLE always overrides both.
if(Python3_FOUND)
    add_custom_command(TARGET cdx_coverage POST_BUILD
            COMMAND bash
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake/setup_circular_env.sh"
            "$<TARGET_FILE_DIR:cdx_coverage>/pyenv"
            "${CMAKE_CURRENT_SOURCE_DIR}/python_script/requirements.txt"
            "${Python3_EXECUTABLE}"
            COMMENT "Provisioning bundled Python environment for circular graphs (best-effort)"
    )
endif()

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------
install(TARGETS cdx_coverage DESTINATION bin)
install(FILES python_script/circular_plot.py DESTINATION bin)

# Bundle the provisioned venv with `make install` too, if it was built
# successfully (OPTIONAL: silently skipped if provisioning failed/was
# skipped, e.g. no network access during the build).
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pyenv" DESTINATION bin OPTIONAL)

# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------
# CDX_BUILD_TESTS lets the tests (and their GoogleTest download/build) be
# skipped entirely, e.g. for a minimal/offline build of cdx_coverage itself.
option(CDX_BUILD_TESTS "Build unit tests (requires network access to fetch GoogleTest on first configure)" ON)

if(CDX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()