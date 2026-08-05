# Portability fix for libvgio's protoc-generated headers.
#
# libvgio/CMakeLists.txt tries two ways to invoke Protobuf's CMake
# integration:
#   1) find_package(Protobuf CONFIG) + protobuf_generate() -> writes
#      vg.pb.h FLAT at <libvgio_build_dir>/vg.pb.h, and libvgio's own
#      "vg -> ." symlink (its own "link_target" custom target) makes
#      "vg/vg.pb.h" resolve correctly.
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
# IMPORTANT: this script deliberately does NOT touch libvgio's own "vg"
# symlink/entry. That target ("link_target") is an `add_custom_target(...
# ALL COMMAND ... create_symlink . vg)` with no declared OUTPUT, so Ninja/
# Make re-runs it on *every* build, unconditionally recreating "vg" as a
# plain symlink pointing at ".". An earlier version of this script replaced
# "vg" with a real directory the first time it ran, which then made every
# subsequent build's "create_symlink . vg" step fail outright (fatal on
# macOS: "Operation not permitted" trying to unlink a non-empty directory;
# merely fragile elsewhere). Fighting over the same path across builds is
# not fixable while link_target keeps unconditionally recreating it - so
# instead this script builds its OWN parallel "vg/vg.pb.h" under a
# different directory (cdx_vg_include/vg/...), which the parent
# CMakeLists.txt adds to the relevant targets' include paths. libvgio's
# native "vg -> ." symlink is left completely alone, whether it's correct
# (case 1, flat) or wrong (case 2, deps/) - it simply doesn't matter
# anymore, because our own directory covers "vg/vg.pb.h" independently of
# it.

if(NOT DEFINED LIBVGIO_BINARY_DIR)
    message(FATAL_ERROR "LIBVGIO_BINARY_DIR must be defined")
endif()

set(_fixup_dir "${LIBVGIO_BINARY_DIR}/cdx_vg_include")

file(GLOB_RECURSE _candidate_headers "${LIBVGIO_BINARY_DIR}/*vg.pb.h")
# Exclude anything already sitting under our own fixup dir from a previous run.
list(FILTER _candidate_headers EXCLUDE REGEX "/cdx_vg_include/vg/vg\\.pb\\.h$")

if(NOT _candidate_headers)
    message(FATAL_ERROR
        "Could not locate a protoc-generated vg.pb.h anywhere under "
        "${LIBVGIO_BINARY_DIR}. Did the run_protoc step succeed?")
endif()

list(GET _candidate_headers 0 _hdr_path)
string(REGEX REPLACE "\\.h$" ".cc" _src_path "${_hdr_path}")

# Our own directory - never touched by libvgio's own build rules, so it's
# safe to remove/recreate on every run without racing anything.
file(REMOVE_RECURSE "${_fixup_dir}")
file(MAKE_DIRECTORY "${_fixup_dir}/vg")
file(CREATE_LINK "${_hdr_path}" "${_fixup_dir}/vg/vg.pb.h" SYMBOLIC)

if(EXISTS "${_src_path}")
    file(CREATE_LINK "${_src_path}" "${_fixup_dir}/vg/vg.pb.cc" SYMBOLIC)
endif()

message(STATUS "libvgio proto fixup: vg/vg.pb.h -> ${_hdr_path} (via ${_fixup_dir})")
