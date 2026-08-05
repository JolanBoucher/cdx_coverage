# ---------------------------------------------------------------------------
# Unit tests (GoogleTest, fetched from source)
# ---------------------------------------------------------------------------
# Each source module under src/ that has no compile-time dependency on
# another module gets its own standalone test executable, compiling that
# module's single .cpp file directly alongside its test file (see
# CMakeLists.txt discussion: no shared internal library, since the tested
# modules don't depend on one another).

include(FetchContent)

# Pin to a specific release so builds are reproducible; bump deliberately.
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
)

# Keep GoogleTest's own install/build options out of this project's install
# rules and avoid rebuilding it against a different runtime than the rest
# of cdx_coverage.
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# ---------------------------------------------------------------------------
# cov_projection
# ---------------------------------------------------------------------------
add_executable(cov_projection_test
        cov_projection_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/cov_projection.cpp"
)

target_compile_features(cov_projection_test PRIVATE cxx_std_17)

target_link_libraries(cov_projection_test PRIVATE
        cdx_lib
        GTest::gtest_main
)

gtest_discover_tests(cov_projection_test)

# ---------------------------------------------------------------------------
# query_resolver
# ---------------------------------------------------------------------------
add_executable(query_resolver_test
        query_resolver_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/query_resolver.cpp"
)

target_compile_features(query_resolver_test PRIVATE cxx_std_17)

# No cdx_lib dependency here: query_resolver.h/.cpp are self-contained
# (only <string>/<unordered_map>/... from the standard library).
target_link_libraries(query_resolver_test PRIVATE
        GTest::gtest_main
)

gtest_discover_tests(query_resolver_test)

# ---------------------------------------------------------------------------
# output_stats
# ---------------------------------------------------------------------------
add_executable(output_stats_test
        output_stats_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/output_stats.cpp"
)

target_compile_features(output_stats_test PRIVATE cxx_std_17)

# Needed for cdx_types.h (used by output_stats.h/gam_io.h), pulled in via
# cdx_lib's propagated include directories - no compiled cdx_lib code is
# actually linked into or called by this test.
target_link_libraries(output_stats_test PRIVATE
        cdx_lib
        GTest::gtest_main
)

gtest_discover_tests(output_stats_test)

# ---------------------------------------------------------------------------
# output_coverage
# ---------------------------------------------------------------------------
add_executable(output_coverage_test
        output_coverage_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/output_coverage.cpp"
)

target_compile_features(output_coverage_test PRIVATE cxx_std_17)

# output_coverage.cpp uses std::thread internally (DoubleBufferedTsvWriter).
find_package(Threads REQUIRED)

target_link_libraries(output_coverage_test PRIVATE
        cdx_lib
        Threads::Threads
        GTest::gtest_main
)

gtest_discover_tests(output_coverage_test)

# ---------------------------------------------------------------------------
# cdx_loader
# ---------------------------------------------------------------------------
add_executable(cdx_loader_test
        cdx_loader_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/cdx_loader.cpp"
)

target_compile_features(cdx_loader_test PRIVATE cxx_std_17)

# Needs cdx_lib for cdx_types.h/cdx_format.h/cdx_IO.h (real header/component
# I/O), and links the compiled cdx_lib code: test fixtures serialize valid
# binary .cdx files on disk using cdx::CdxFormat::pack_* directly.
target_link_libraries(cdx_loader_test PRIVATE
        cdx_lib
        GTest::gtest_main
)

gtest_discover_tests(cdx_loader_test)

# ---------------------------------------------------------------------------
# gam_io
# ---------------------------------------------------------------------------
add_executable(gam_io_test
        gam_io_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/gam_io.cpp"
)

target_compile_features(gam_io_test PRIVATE cxx_std_17)

# Fixtures are built with the real vg::io serialization stack (ProtobufEmitter<Alignment>),
# so this target needs the same dependencies as gam_io.cpp itself / the main cdx_coverage
# executable: libvgio, Protobuf, the Abseil components Protobuf's generated code pulls in,
# and OpenMP (process_gam is parallelized with #pragma omp).
target_link_libraries(gam_io_test PRIVATE
        cdx_lib
        ${CDX_VGIO_TARGET}
        protobuf::libprotobuf
        absl::log_internal_message
        absl::log_internal_check_op
        absl::hash
        absl::status
        absl::raw_logging_internal
        OpenMP::OpenMP_CXX
        GTest::gtest_main
)

gtest_discover_tests(gam_io_test)

# ---------------------------------------------------------------------------
# cli
# ---------------------------------------------------------------------------
add_executable(cli_test
        cli_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/cli.cpp"
)

target_compile_features(cli_test PRIVATE cxx_std_17)

# cli.cpp has no cdx_lib dependency, but it does need CLI11 - resolved
# portably via the CLI11::CLI11 imported target set up in the root
# CMakeLists.txt (same one the main cdx_coverage executable links against),
# so this works on any machine/CI without a hardcoded path.
target_link_libraries(cli_test PRIVATE
        CLI11::CLI11
        GTest::gtest_main
)

gtest_discover_tests(cli_test)

# ---------------------------------------------------------------------------
# linear_plot
# ---------------------------------------------------------------------------
# Pure numerical preparation (downsampling, smoothing, tick math, log
# transform) - no Cairo/PNG/subprocess dependency, so this gets full unit
# test coverage rather than the lighter integration treatment used for the
# actual rendering backends (cairo_plot.cpp, python_circular_plot.cpp).
add_executable(linear_plot_test
        linear_plot_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/linear_plot.cpp"
)

target_compile_features(linear_plot_test PRIVATE cxx_std_17)

target_link_libraries(linear_plot_test PRIVATE
        cdx_lib
        GTest::gtest_main
)

gtest_discover_tests(linear_plot_test)

# ---------------------------------------------------------------------------
# cairo_plot
# ---------------------------------------------------------------------------
# Light integration tests: cairo_plot.cpp's internals are all in an anonymous
# namespace, so only the two public entry points (writeLinearPlotQuery,
# writeLinearPlotGlobal) are exercised, by checking the PNG files they
# produce (signature + IHDR-declared dimensions) rather than pixel content.
# Both writeLinearPlotQuery/Global call prepareLinearPlotPackage() and
# chooseGlobalGraphGrid(), defined in linear_plot.cpp, so that source file
# is compiled into this target too - it is not linked from a shared library.
add_executable(cairo_plot_test
        cairo_plot_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/cairo_plot.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/linear_plot.cpp"
)

target_compile_features(cairo_plot_test PRIVATE cxx_std_17)

target_link_libraries(cairo_plot_test PRIVATE
        cdx_lib
        PkgConfig::CAIRO
        GTest::gtest_main
)

gtest_discover_tests(cairo_plot_test)

# ---------------------------------------------------------------------------
# python_circular_plot (C++ side: light integration tests)
# ---------------------------------------------------------------------------
# Same anonymous-namespace situation as cairo_plot.cpp: only writeCircularPlotQuery/
# writeCircularPlotGlobal are reachable. Rather than depending on the real Python
# pycirclize/matplotlib rendering stack (heavyweight, optional at runtime, not always
# installed), most tests here point CDX_CIRCULAR_PLOT_SCRIPT at tiny stub .py scripts
# (always-succeeds / always-fails) written by the test itself, so the C++ orchestration
# (temp dirs, binary request serialization, subprocess invocation, error propagation,
# cleanup) is exercised end-to-end with a real `python3` subprocess, without ever
# requiring pycirclize. See circular_plot_test.cpp's file docstring for the full
# breakdown. Needs linear_plot.cpp too, for chooseGlobalGraphGrid()/prepareCoverageForPlot().
add_executable(circular_plot_test
        circular_plot_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/python_circular_plot.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/linear_plot.cpp"
)

target_compile_features(circular_plot_test PRIVATE cxx_std_17)

target_link_libraries(circular_plot_test PRIVATE
        cdx_lib
        GTest::gtest_main
)

gtest_discover_tests(circular_plot_test)

# ---------------------------------------------------------------------------
# query_plot_slice
# ---------------------------------------------------------------------------
# sliceLinearQueryCoverage() was extracted out of main.cpp specifically to be unit-testable
# (see src/query_plot_slice.h/.cpp) - pure computation, only needs cdx_types.h/config.h.
add_executable(query_plot_slice_test
        query_plot_slice_test.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../src/query_plot_slice.cpp"
)

target_compile_features(query_plot_slice_test PRIVATE cxx_std_17)

target_link_libraries(query_plot_slice_test PRIVATE
        cdx_lib
        GTest::gtest_main
)

gtest_discover_tests(query_plot_slice_test)

# ---------------------------------------------------------------------------
# main.cpp end-to-end smoke tests (run the real, fully-built cdx_coverage binary)
# ---------------------------------------------------------------------------
# main.cpp itself is pure orchestration of already-unit-tested modules and has no reachable
# public API, so it's verified by actually running the compiled cdx_coverage executable
# end-to-end against small real CDX+GAM fixtures - see main_e2e_test.cpp's file docstring.
# This needs the same fixture-building dependencies as cdx_loader_test/gam_io_test
# (cdx_lib + libvgio/Protobuf/Abseil), plus a hard build-order dependency on cdx_coverage
# itself and its executable path passed in as a compile definition.
add_executable(main_e2e_test
        main_test.cpp
)

target_compile_features(main_e2e_test PRIVATE cxx_std_17)

target_link_libraries(main_e2e_test PRIVATE
        cdx_lib
        ${CDX_VGIO_TARGET}
        protobuf::libprotobuf
        absl::log_internal_message
        absl::log_internal_check_op
        absl::hash
        absl::status
        absl::raw_logging_internal
        GTest::gtest_main
)

add_dependencies(main_e2e_test cdx_coverage)
target_compile_definitions(main_e2e_test PRIVATE
        CDX_COVERAGE_EXECUTABLE_PATH="$<TARGET_FILE:cdx_coverage>"
)

gtest_discover_tests(main_e2e_test)

# ---------------------------------------------------------------------------
# circular_plot.py (pure-Python unit tests, independent of the C++ suite)
# ---------------------------------------------------------------------------
# Standalone unittest suite for python_script/circular_plot.py itself (binary
# request parsing, NaN-gap segment splitting, radial tick logic). pycirclize
# is mocked inside the test file, so it only needs a Python3 interpreter with
# numpy/matplotlib - both real imports in circular_plot.py - available; if
# either is missing, every test in the file self-skips rather than erroring
# (see _DEPENDENCIES_AVAILABLE in circular_plot_py_test.py).
if(Python3_FOUND)
    add_test(
            NAME circular_plot_py_test
            COMMAND ${Python3_EXECUTABLE} -m unittest circular_plot_py_test -v
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
endif()
