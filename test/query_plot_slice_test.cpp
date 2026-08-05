/**
 * @file query_plot_slice_test.cpp
 * @brief Unit tests for query_plot_slice.cpp (sliceLinearQueryCoverage).
 *
 * sliceLinearQueryCoverage() was extracted out of main.cpp's runQueryPipeline() (the linear
 * graph branch) specifically so it could get real unit tests: unlike the rest of main.cpp, this
 * one piece of logic (origin-crossing rejection, out-of-range validation, sub-vector slicing to
 * align the plotted X-axis with the requested query range) has no I/O and no dependency on the
 * rest of the pipeline, so it's exercised directly here rather than only implicitly through a
 * full end-to-end run of the real binary (see test/main_e2e_test.cpp for that).
 *
 * sliceLinearQueryCoverage() is a free function (not in an anonymous namespace), so it's called
 * directly rather than through some other public entry point.
 */

#include <gtest/gtest.h>

#include "../src/query_plot_slice.h"
#include "../src/config.h"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace {
    std::vector<cdx::Coverage> makeCoverage(const std::size_t n) {
        std::vector<cdx::Coverage> coverage;
        coverage.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            coverage.push_back(static_cast<cdx::Coverage>(i * 10));
        }
        return coverage;
    }
} // namespace

TEST(SliceLinearQueryCoverageTest, MiddleRangeReturnsExactInclusiveSlice) {
    const auto coverage = makeCoverage(10); // {0,10,20,...,90}
    const auto result = output::sliceLinearQueryCoverage(coverage, {3, 6});

    EXPECT_EQ(result, (std::vector<cdx::Coverage>{30, 40, 50, 60}));
}

TEST(SliceLinearQueryCoverageTest, FullComponentRangeReturnsWholeVectorUnchanged) {
    const auto coverage = makeCoverage(5);
    const auto result = output::sliceLinearQueryCoverage(coverage, {0, 4});

    EXPECT_EQ(result, coverage);
}

TEST(SliceLinearQueryCoverageTest, SingleBasePairRangeReturnsOneElement) {
    const auto coverage = makeCoverage(10);
    const auto result = output::sliceLinearQueryCoverage(coverage, {5, 5});

    EXPECT_EQ(result, (std::vector<cdx::Coverage>{50}));
}

TEST(SliceLinearQueryCoverageTest, FirstElementOnlyRangeWorks) {
    const auto coverage = makeCoverage(10);
    const auto result = output::sliceLinearQueryCoverage(coverage, {0, 0});

    EXPECT_EQ(result, (std::vector<cdx::Coverage>{0}));
}

TEST(SliceLinearQueryCoverageTest, LastValidEndAtSizeMinusOneWorks) {
    const auto coverage = makeCoverage(10);
    const auto result = output::sliceLinearQueryCoverage(coverage, {9, 9});

    EXPECT_EQ(result, (std::vector<cdx::Coverage>{90}));
}

TEST(SliceLinearQueryCoverageTest, StartGreaterThanEndThrowsInvalidArgument) {
    const auto coverage = makeCoverage(10);
    EXPECT_THROW(output::sliceLinearQueryCoverage(coverage, {6, 3}), std::invalid_argument);
}

TEST(SliceLinearQueryCoverageTest, StartOneBeyondEndThrowsInvalidArgument) {
    // The smallest possible "crossing" case: start is exactly one past end.
    const auto coverage = makeCoverage(10);
    EXPECT_THROW(output::sliceLinearQueryCoverage(coverage, {4, 3}), std::invalid_argument);
}

TEST(SliceLinearQueryCoverageTest, EndEqualToSizeThrowsOutOfRange) {
    const auto coverage = makeCoverage(10); // valid indices [0,9]
    EXPECT_THROW(output::sliceLinearQueryCoverage(coverage, {0, 10}), std::out_of_range);
}

TEST(SliceLinearQueryCoverageTest, EndFarBeyondSizeThrowsOutOfRange) {
    const auto coverage = makeCoverage(10);
    EXPECT_THROW(output::sliceLinearQueryCoverage(coverage, {0, 1000}), std::out_of_range);
}

TEST(SliceLinearQueryCoverageTest, EmptyCoverageWithZeroZeroThrowsOutOfRange) {
    const std::vector<cdx::Coverage> coverage; // size 0
    // {0,0} is not a crossing range (start<=end), but end(0) >= size(0).
    EXPECT_THROW(output::sliceLinearQueryCoverage(coverage, {0, 0}), std::out_of_range);
}

TEST(SliceLinearQueryCoverageTest, SameStartAndEndBeyondSizeThrowsOutOfRangeNotInvalidArgument) {
    // Regression guard: start==end (not a crossing case) must be reported via out_of_range
    // when out of bounds, never misclassified as the origin-crossing invalid_argument case.
    const auto coverage = makeCoverage(5);
    EXPECT_THROW(output::sliceLinearQueryCoverage(coverage, {100, 100}), std::out_of_range);
}

TEST(SliceLinearQueryCoverageTest, PreservesValuesAndOrderExactly) {
    const std::vector<cdx::Coverage> coverage{7, 3, 9, 1, 4, 8};
    const auto result = output::sliceLinearQueryCoverage(coverage, {1, 4});

    EXPECT_EQ(result, (std::vector<cdx::Coverage>{3, 9, 1, 4}));
}

TEST(SliceLinearQueryCoverageTest, SentinelValuesArePassedThroughUnfiltered) {
    // sliceLinearQueryCoverage only slices; sentinel filtering (NOT_IN_QUERY / NOT_IN_COMPO)
    // is prepareLinearPlotPackage's responsibility downstream, not this function's.
    const std::vector<cdx::Coverage> coverage{1, cfg::NOT_IN_QUERY, 3, cfg::NOT_IN_COMPO, 5};
    const auto result = output::sliceLinearQueryCoverage(coverage, {0, 4});

    EXPECT_EQ(result, coverage);
}

TEST(SliceLinearQueryCoverageTest, LargeVectorMiddleSliceIsExact) {
    const auto coverage = makeCoverage(10000);
    const auto result = output::sliceLinearQueryCoverage(coverage, {4999, 5001});

    EXPECT_EQ(result, (std::vector<cdx::Coverage>{49990, 50000, 50010}));
}
