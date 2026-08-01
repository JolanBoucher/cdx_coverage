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
)

target_compile_features(cdx_coverage PRIVATE cxx_std_17)

target_include_directories(cdx_coverage PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

add_subdirectory(
        ../cdx_lib
        cdx_lib_build
)

target_link_libraries(
        cdx_coverage
        PRIVATE
        cdx_lib
)

target_link_libraries(cdx_coverage PRIVATE
        ${CDX_VGIO_TARGET}
        protobuf::libprotobuf
        absl::log_internal_message
        absl::log_internal_check_op
        absl::hash
        absl::status
        absl::raw_logging_internal
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(cdx_coverage PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()