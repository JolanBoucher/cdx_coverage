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
