# Portability fix for libvgio's protoc-generated headers.
#
# libvgio/CMakeLists.txt tries two ways to invoke Protobuf's CMake
# integration:
#   1) find_package(Protobuf CONFIG) + protobuf_generate() -> writes
#      vg.pb.h FLAT at <libvgio_build_dir>/vg.pb.h, and libvgio's own
#      "vg -> ." symlink makes "vg/vg.pb.h" resolve correctly.
#   2) fallback find_package(Protobuf) (legacy FindProtobuf.cmake module)
#      + protobuf_generate_cpp() -> PRESERVES the "deps/" subpath, writing
#      vg.pb.h at <libvgio_build_dir>/deps/vg.pb.h instead.
#
# Which path is taken depends on whether the system's installed Protobuf
# ships a protobuf-config.cmake (varies by distro/version - e.g. present on
# Ubuntu 24.04's libprotobuf-dev, absent on Ubuntu 20.04's). In case (2),
# libvgio's flat "vg -> ." symlink does NOT cover "vg/vg.pb.h", and any
# source doing #include "vg/vg.pb.h" fails to compile.
#
# This script runs after libvgio's protoc step, locates wherever vg.pb.h
# actually landed, and (re)builds "<libvgio_build_dir>/vg/" as a small
# directory of symlinks pointing at the real generated files - so the
# include always resolves, regardless of which Protobuf CMake integration
# the host machine happens to have.

if(NOT DEFINED LIBVGIO_BINARY_DIR)
    message(FATAL_ERROR "LIBVGIO_BINARY_DIR must be defined")
endif()

file(GLOB_RECURSE _candidate_headers "${LIBVGIO_BINARY_DIR}/*vg.pb.h")
# Exclude anything already sitting under the "vg/" fixup dir from a previous run.
list(FILTER _candidate_headers EXCLUDE REGEX "/vg/vg\\.pb\\.h$")

if(NOT _candidate_headers)
    message(FATAL_ERROR
        "Could not locate a protoc-generated vg.pb.h anywhere under "
        "${LIBVGIO_BINARY_DIR}. Did the run_protoc step succeed?")
endif()

list(GET _candidate_headers 0 _hdr_path)
string(REGEX REPLACE "\\.h$" ".cc" _src_path "${_hdr_path}")

# libvgio's own link_target step may have already created "vg" as a symlink
# pointing to "." (i.e. to LIBVGIO_BINARY_DIR itself, self-referentially).
# file(REMOVE_RECURSE) follows symlinks-to-directories, so calling it on that
# self-referential symlink recurses into itself forever and crashes CMake.
# Only ever remove it as a plain (non-recursing) unlink when it's a symlink;
# only recurse when it's a genuine directory.
if(IS_SYMLINK "${LIBVGIO_BINARY_DIR}/vg")
    file(REMOVE "${LIBVGIO_BINARY_DIR}/vg")
elseif(IS_DIRECTORY "${LIBVGIO_BINARY_DIR}/vg")
    file(REMOVE_RECURSE "${LIBVGIO_BINARY_DIR}/vg")
endif()
file(MAKE_DIRECTORY "${LIBVGIO_BINARY_DIR}/vg")
file(CREATE_LINK "${_hdr_path}" "${LIBVGIO_BINARY_DIR}/vg/vg.pb.h" SYMBOLIC)

if(EXISTS "${_src_path}")
    file(CREATE_LINK "${_src_path}" "${LIBVGIO_BINARY_DIR}/vg/vg.pb.cc" SYMBOLIC)
endif()

message(STATUS "libvgio proto fixup: vg/vg.pb.h -> ${_hdr_path}")
