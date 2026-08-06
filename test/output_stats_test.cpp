/**
 * @file output_stats_test.cpp
 * @brief Unit tests for src/output_stats.cpp.
 *
 * Three suites, one per public responsibility of the module:
 *   - ComputeCoverageStatsTest exercises output::computeCoverageStats(), the
 *     pure statistics-computation function, with exact reference values
 *     computed independently (see the Python simulation used to derive
 *     them - the dense-histogram and nth_element algorithms are each
 *     reimplemented faithfully in that script, not just "expected output
 *     copied from a first run of the C++ code").
 *   - WriteStatsReportQueryTest / WriteStatsReportGlobalTest exercise the
 *     validation logic and file-writing side effects of the two report
 *     generators. Report *content* is checked at the substring level
 *     (section headers, component names, absence of nan/inf) rather than
 *     matched character-for-character, since exact formatting is a
 *     presentation detail already indirectly covered by the numeric tests
 *     in ComputeCoverageStatsTest; a full-text match would be fragile and
 *     wouldn't add real correctness coverage.
 *
 * Note: computeCoverageStatsNthElement() is a file-local (anonymous
 * namespace) implementation detail, not part of the public API, so it
 * cannot be invoked directly nor directly cross-checked against the
 * dense-histogram path on the same input from outside this translation
 * unit. It is instead exercised on its own through computeCoverageStats()
 * by supplying data whose maximum forces the fallback dispatch.
 */

#include "../src/output_stats.h"
#include "../src/config.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {
    // Absolute tolerance for floating-point statistic comparisons.
    constexpr double kEps = 1e-9;

    // Unique temp file path per test, under the system temp directory;
    // callers are responsible for removing it (see TearDownFile below).
    std::filesystem::path tempReportPath(const std::string& test_name) {
        return std::filesystem::temp_directory_path() / ("cdx_output_stats_test_" + test_name + ".txt");
    }

    void removeIfExists(const std::filesystem::path& path) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    std::string readFile(const std::filesystem::path& path) {
        std::ifstream in(path);
        std::ostringstream buffer;
        buffer << in.rdbuf();
        return buffer.str();
    }
}

// =============================================================================
// computeCoverageStats
//
// Computes summary statistics (breadth, mean, median, stddev, cv,
// quartiles, extrema) from a coverage vector, automatically dispatching
// between a dense-histogram implementation and an nth_element-based
// fallback for coverage ranges exceeding MAX_DENSE_HISTOGRAM (1,000,000).
// Coverage values >= cfg::NOT_IN_QUERY are sentinels and excluded from all
// calculations.
// =============================================================================

// An empty coverage vector produces an all-zero result, no NaN/division
// by zero.
TEST(ComputeCoverageStatsTest, EmptyCoverageReturnsAllZero) {
    const std::vector<cdx::Coverage> coverage;

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 0u);
    EXPECT_EQ(stats.covered_positions, 0u);
    EXPECT_EQ(stats.min, 0u);
    EXPECT_EQ(stats.max, 0u);
    EXPECT_NEAR(stats.mean, 0.0, kEps);
    EXPECT_NEAR(stats.breadth, 0.0, kEps);
}

// A coverage vector containing only sentinel values behaves identically to
// an empty one: nothing is a valid position.
TEST(ComputeCoverageStatsTest, AllSentinelCoverageIsEquivalentToEmpty) {
    const std::vector<cdx::Coverage> coverage = {cfg::NOT_IN_QUERY, cfg::NOT_IN_COMPO};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 0u);
    EXPECT_EQ(stats.min, 0u);
    EXPECT_EQ(stats.max, 0u);
}

// All-zero (but valid, non-sentinel) coverage: real positions with no
// reads, distinct from the "no valid positions at all" case above.
TEST(ComputeCoverageStatsTest, AllZeroCoverage) {
    const std::vector<cdx::Coverage> coverage = {0, 0, 0};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 3u);
    EXPECT_EQ(stats.covered_positions, 0u);
    EXPECT_NEAR(stats.breadth, 0.0, kEps);
    EXPECT_NEAR(stats.mean, 0.0, kEps);
    EXPECT_NEAR(stats.stddev, 0.0, kEps);
    EXPECT_NEAR(stats.cv, 0.0, kEps); // mean == 0, so cv stays 0 (no division).
}

// Uniform non-zero coverage: zero spread, cv well-defined (mean > 0).
TEST(ComputeCoverageStatsTest, UniformNonzeroCoverageHasZeroSpread) {
    const std::vector<cdx::Coverage> coverage = {7, 7, 7, 7};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 4u);
    EXPECT_EQ(stats.covered_positions, 4u);
    EXPECT_NEAR(stats.breadth, 100.0, kEps);
    EXPECT_NEAR(stats.mean, 7.0, kEps);
    EXPECT_NEAR(stats.stddev, 0.0, kEps);
    EXPECT_NEAR(stats.cv, 0.0, kEps);
    EXPECT_NEAR(stats.q1, 7.0, kEps);
    EXPECT_NEAR(stats.median, 7.0, kEps);
    EXPECT_NEAR(stats.q3, 7.0, kEps);
    EXPECT_EQ(stats.min, 7u);
    EXPECT_EQ(stats.max, 7u);
}

// N = 5: quartile ranks land exactly on integers (no interpolation
// fraction), a useful "obviously correct" reference case.
TEST(ComputeCoverageStatsTest, ExactIntegerQuantileRanks) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4, 5};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 5u);
    EXPECT_EQ(stats.covered_positions, 5u);
    EXPECT_NEAR(stats.breadth, 100.0, kEps);
    EXPECT_NEAR(stats.mean, 3.0, kEps);
    EXPECT_NEAR(stats.stddev, std::sqrt(2.0), kEps);
    EXPECT_NEAR(stats.cv, std::sqrt(2.0) / 3.0, kEps);
    EXPECT_NEAR(stats.q1, 2.0, kEps);
    EXPECT_NEAR(stats.median, 3.0, kEps);
    EXPECT_NEAR(stats.q3, 4.0, kEps);
    EXPECT_EQ(stats.min, 1u);
    EXPECT_EQ(stats.max, 5u);
}

// N = 4: quartile ranks are fractional (0.75 / 1.5 / 2.25), genuinely
// exercising the linear-interpolation branch between two order statistics.
TEST(ComputeCoverageStatsTest, FractionalQuantileRanksInterpolate) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 4u);
    EXPECT_NEAR(stats.mean, 2.5, kEps);
    EXPECT_NEAR(stats.stddev, 1.118033988749895, 1e-12);
    EXPECT_NEAR(stats.cv, 0.447213595499958, 1e-12);
    EXPECT_NEAR(stats.q1, 1.75, kEps);
    EXPECT_NEAR(stats.median, 2.5, kEps);
    EXPECT_NEAR(stats.q3, 3.25, kEps);
}

// Sentinel values interleaved among valid ones are excluded from every
// statistic; only the valid subset {10, 20, 30} should be reflected.
TEST(ComputeCoverageStatsTest, MixedSentinelAndValidValuesExcludesSentinels) {
    const std::vector<cdx::Coverage> coverage = {10, cfg::NOT_IN_QUERY, 20, cfg::NOT_IN_COMPO, 30};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 3u);
    EXPECT_EQ(stats.covered_positions, 3u);
    EXPECT_NEAR(stats.mean, 20.0, kEps);
    EXPECT_NEAR(stats.stddev, 8.16496580927726, 1e-11);
    EXPECT_NEAR(stats.q1, 15.0, kEps);
    EXPECT_NEAR(stats.median, 20.0, kEps);
    EXPECT_NEAR(stats.q3, 25.0, kEps);
    EXPECT_EQ(stats.min, 10u);
    EXPECT_EQ(stats.max, 30u);
}

// A single value exactly at MAX_DENSE_HISTOGRAM (1,000,000) must remain on
// the dense-histogram path (the check is a strict '>').
TEST(ComputeCoverageStatsTest, DenseHistogramThresholdBoundaryStaysDense) {
    const std::vector<cdx::Coverage> coverage = {1'000'000};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 1u);
    EXPECT_NEAR(stats.mean, 1'000'000.0, kEps);
    EXPECT_NEAR(stats.stddev, 0.0, kEps);
    EXPECT_EQ(stats.min, 1'000'000u);
    EXPECT_EQ(stats.max, 1'000'000u);
}

// One unit above the threshold switches to the nth_element fallback; the
// result must still be correct, not just "doesn't crash".
TEST(ComputeCoverageStatsTest, FallbackThresholdBoundaryTriggersFallback) {
    const std::vector<cdx::Coverage> coverage = {1'000'001};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 1u);
    EXPECT_NEAR(stats.mean, 1'000'001.0, kEps);
    EXPECT_NEAR(stats.stddev, 0.0, kEps);
    EXPECT_EQ(stats.min, 1'000'001u);
    EXPECT_EQ(stats.max, 1'000'001u);
}

// A small (N=3) dataset with one very large value: exercises the
// nth_element fallback path end-to-end with genuine interpolation.
TEST(ComputeCoverageStatsTest, NthElementFallbackProducesCorrectStats) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 2'000'000};

    const auto stats = output::computeCoverageStats(coverage);

    EXPECT_EQ(stats.region_length, 3u);
    EXPECT_EQ(stats.covered_positions, 3u);
    EXPECT_NEAR(stats.mean, 666667.6666666666, 1e-6);
    EXPECT_NEAR(stats.stddev, 942808.3344753706, 1e-6);
    EXPECT_NEAR(stats.cv, 1.4142103803974855, 1e-9);
    EXPECT_NEAR(stats.q1, 1.5, kEps);
    EXPECT_NEAR(stats.median, 2.0, kEps);
    EXPECT_NEAR(stats.q3, 1000001.0, 1e-6);
    EXPECT_EQ(stats.min, 1u);
    EXPECT_EQ(stats.max, 2'000'000u);
}

// =============================================================================
// writeStatsReportQuery
//
// Validates internal consistency of GamMappingStats, computes coverage
// statistics for a single component, and writes a formatted text report.
// =============================================================================

// mapped > total is an inconsistent mapping-statistics state.
TEST(WriteStatsReportQueryTest, MappedExceedsTotalThrows) {
    const GamMappingStats stats{/*total*/ 10, /*mapped*/ 20, /*mapped_to_query*/ 0, /*unmapped*/ 0};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempReportPath("mapped_exceeds_total");

    EXPECT_THROW(
        output::writeStatsReportQuery(path, stats, coverage, "chr1"),
        std::invalid_argument);
    removeIfExists(path);
}

// unmapped > total is likewise inconsistent.
TEST(WriteStatsReportQueryTest, UnmappedExceedsTotalThrows) {
    const GamMappingStats stats{/*total*/ 10, /*mapped*/ 5, /*mapped_to_query*/ 0, /*unmapped*/ 20};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempReportPath("unmapped_exceeds_total");

    EXPECT_THROW(
        output::writeStatsReportQuery(path, stats, coverage, "chr1"),
        std::invalid_argument);
    removeIfExists(path);
}

// mapped_to_query cannot exceed mapped (a subset relationship).
TEST(WriteStatsReportQueryTest, MappedToQueryExceedsMappedThrows) {
    const GamMappingStats stats{/*total*/ 10, /*mapped*/ 5, /*mapped_to_query*/ 8, /*unmapped*/ 5};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempReportPath("mapped_to_query_exceeds_mapped");

    EXPECT_THROW(
        output::writeStatsReportQuery(path, stats, coverage, "chr1"),
        std::invalid_argument);
    removeIfExists(path);
}

// mapped + unmapped must sum to exactly total.
TEST(WriteStatsReportQueryTest, MappedPlusUnmappedMustEqualTotalThrows) {
    const GamMappingStats stats{/*total*/ 10, /*mapped*/ 5, /*mapped_to_query*/ 0, /*unmapped*/ 4};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempReportPath("mapped_plus_unmapped_mismatch");

    EXPECT_THROW(
        output::writeStatsReportQuery(path, stats, coverage, "chr1"),
        std::invalid_argument);
    removeIfExists(path);
}

// total == 0 must not produce a division-by-zero (NaN/inf) anywhere in the
// percentage calculations or the written report.
TEST(WriteStatsReportQueryTest, ZeroTotalReadsProducesNoDivisionByZero) {
    const GamMappingStats stats{}; // all zero, still internally consistent (0+0==0)
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempReportPath("zero_total_reads");

    output::writeStatsReportQuery(path, stats, coverage, "chr1");
    const std::string content = readFile(path);

    EXPECT_EQ(content.find("nan"), std::string::npos);
    EXPECT_EQ(content.find("-nan"), std::string::npos);
    EXPECT_EQ(content.find("inf"), std::string::npos);
    removeIfExists(path);
}

// An unopenable output path (nonexistent parent directory) must raise
// std::runtime_error rather than fail silently.
TEST(WriteStatsReportQueryTest, UnopenableOutputPathThrowsRuntimeError) {
    const GamMappingStats stats{10, 8, 5, 2};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const std::filesystem::path bad_path = "/cdx_this_directory_should_not_exist_xyz/report.txt";

    EXPECT_THROW(
        output::writeStatsReportQuery(bad_path, stats, coverage, "chr1"),
        std::runtime_error);
}

// The written report contains the component name and both expected
// section headers.
TEST(WriteStatsReportQueryTest, ReportContainsComponentNameAndSections) {
    const GamMappingStats stats{10, 8, 5, 2};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4, 5};
    const auto path = tempReportPath("report_sections");

    output::writeStatsReportQuery(path, stats, coverage, "chr1");
    const std::string content = readFile(path);

    EXPECT_NE(content.find("chr1"), std::string::npos);
    EXPECT_NE(content.find("Mapping Statistics"), std::string::npos);
    EXPECT_NE(content.find("Coverage Statistics"), std::string::npos);
    removeIfExists(path);
}

// The default precision parameter (CoveragePrecision::Node, kept for
// backward compatibility with pre-existing callers/tests) is recorded as
// "node" in the report header.
TEST(WriteStatsReportQueryTest, DefaultPrecisionLabelledNode) {
    const GamMappingStats stats{10, 8, 5, 2};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempReportPath("default_precision_label");

    output::writeStatsReportQuery(path, stats, coverage, "chr1");
    const std::string content = readFile(path);

    EXPECT_NE(content.find("Coverage precision: node"), std::string::npos);
    removeIfExists(path);
}

// An explicit CoveragePrecision::Base is recorded as "base" instead.
TEST(WriteStatsReportQueryTest, ExplicitBasePrecisionLabelledBase) {
    const GamMappingStats stats{10, 8, 5, 2};
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempReportPath("explicit_base_precision_label");

    output::writeStatsReportQuery(path, stats, coverage, "chr1", CoveragePrecision::Base);
    const std::string content = readFile(path);

    EXPECT_NE(content.find("Coverage precision: base"), std::string::npos);
    removeIfExists(path);
}

// =============================================================================
// writeStatsReportGlobal
//
// Validates that component_offsets/component_names form a complete,
// consistent partition of the flattened coverage table, computes global
// and per-component statistics (in parallel when OpenMP is available), and
// writes one report section per component plus a global summary.
// =============================================================================

// component_offsets must contain at least two boundaries.
TEST(WriteStatsReportGlobalTest, ComponentOffsetsMustContainTwoBoundariesThrows) {
    const std::vector<cdx::Coverage> coverage;
    const std::vector<cdx::PosBp> offsets = {0};
    const std::vector<std::string> names;
    const GamMappingStats stats{};
    const auto path = tempReportPath("offsets_too_short");

    EXPECT_THROW(
        output::writeStatsReportGlobal(path, coverage, offsets, names, stats, 1),
        std::invalid_argument);
    removeIfExists(path);
}

// component_names.size() must equal component_offsets.size() - 1.
TEST(WriteStatsReportGlobalTest, ComponentNameCountMismatchThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 2, 4};
    const std::vector<std::string> names = {"chr1"}; // should have 2 entries
    const GamMappingStats stats{};
    const auto path = tempReportPath("name_count_mismatch");

    EXPECT_THROW(
        output::writeStatsReportGlobal(path, coverage, offsets, names, stats, 1),
        std::invalid_argument);
    removeIfExists(path);
}

// The first offset must be zero.
TEST(WriteStatsReportGlobalTest, FirstOffsetMustBeZeroThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {1, 4};
    const std::vector<std::string> names = {"chr1"};
    const GamMappingStats stats{};
    const auto path = tempReportPath("first_offset_nonzero");

    EXPECT_THROW(
        output::writeStatsReportGlobal(path, coverage, offsets, names, stats, 1),
        std::invalid_argument);
    removeIfExists(path);
}

// The final offset must match the flattened coverage table size.
TEST(WriteStatsReportGlobalTest, FinalOffsetMustMatchCoverageSizeThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 99};
    const std::vector<std::string> names = {"chr1"};
    const GamMappingStats stats{};
    const auto path = tempReportPath("final_offset_mismatch");

    EXPECT_THROW(
        output::writeStatsReportGlobal(path, coverage, offsets, names, stats, 1),
        std::invalid_argument);
    removeIfExists(path);
}

// Offsets must be non-decreasing (a well-formed partition).
TEST(WriteStatsReportGlobalTest, OffsetsMustBeNonDecreasingThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 3, 2, 4};
    const std::vector<std::string> names = {"chr1", "chr2"};
    const GamMappingStats stats{};
    const auto path = tempReportPath("offsets_not_nondecreasing");

    EXPECT_THROW(
        output::writeStatsReportGlobal(path, coverage, offsets, names, stats, 1),
        std::invalid_argument);
    removeIfExists(path);
}

// An unopenable output path must raise std::runtime_error.
TEST(WriteStatsReportGlobalTest, UnopenableOutputPathThrowsRuntimeError) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 4};
    const std::vector<std::string> names = {"chr1"};
    const GamMappingStats stats{};
    const std::filesystem::path bad_path = "/cdx_this_directory_should_not_exist_xyz/report.txt";

    EXPECT_THROW(
        output::writeStatsReportGlobal(bad_path, coverage, offsets, names, stats, 1),
        std::runtime_error);
}

// A single-component graph is the degenerate case of the multi-component
// report; the file must still contain both a global and a per-component
// section.
TEST(WriteStatsReportGlobalTest, SingleComponentReportContainsSections) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4, 5};
    const std::vector<cdx::PosBp> offsets = {0, 5};
    const std::vector<std::string> names = {"chr1"};
    const GamMappingStats stats{10, 8, 5, 2};
    const auto path = tempReportPath("single_component");

    output::writeStatsReportGlobal(path, coverage, offsets, names, stats, 1);
    const std::string content = readFile(path);

    EXPECT_NE(content.find("Global Coverage Statistics"), std::string::npos);
    EXPECT_NE(content.find("chr1 Coverage Statistics"), std::string::npos);
    removeIfExists(path);
}

// One component small enough for the dense-histogram path, one component
// whose max forces the nth_element fallback: global statistics must be
// recomputed from the full table (global_fallback == true) rather than
// merging incompatible accumulators, and the report must still contain a
// section for every component without crashing or emitting nan/inf.
TEST(WriteStatsReportGlobalTest, MixedDenseAndFallbackComponentsProducesSaneReport) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 2'000'000, 1, 1};
    const std::vector<cdx::PosBp> offsets = {0, 3, 6};
    const std::vector<std::string> names = {"chr1", "chr2"};
    const GamMappingStats stats{10, 8, 5, 2};
    const auto path = tempReportPath("mixed_dense_fallback");

    output::writeStatsReportGlobal(path, coverage, offsets, names, stats, 1);
    const std::string content = readFile(path);

    EXPECT_NE(content.find("Global Coverage Statistics"), std::string::npos);
    EXPECT_NE(content.find("chr1 Coverage Statistics"), std::string::npos);
    EXPECT_NE(content.find("chr2 Coverage Statistics"), std::string::npos);
    EXPECT_EQ(content.find("nan"), std::string::npos);
    EXPECT_EQ(content.find("inf"), std::string::npos);
    removeIfExists(path);
}

// Same precision-label check as WriteStatsReportQueryTest above, exercising
// the trailing defaulted parameter and an explicit override.
TEST(WriteStatsReportGlobalTest, PrecisionLabelDefaultsToNodeAndHonoursOverride) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4, 5};
    const std::vector<cdx::PosBp> offsets = {0, 5};
    const std::vector<std::string> names = {"chr1"};
    const GamMappingStats stats{10, 8, 5, 2};

    const auto default_path = tempReportPath("global_default_precision_label");
    output::writeStatsReportGlobal(default_path, coverage, offsets, names, stats, 1);
    EXPECT_NE(readFile(default_path).find("Coverage precision: node"), std::string::npos);
    removeIfExists(default_path);

    const auto base_path = tempReportPath("global_explicit_base_precision_label");
    output::writeStatsReportGlobal(base_path, coverage, offsets, names, stats, 1, CoveragePrecision::Base);
    EXPECT_NE(readFile(base_path).find("Coverage precision: base"), std::string::npos);
    removeIfExists(base_path);
}
