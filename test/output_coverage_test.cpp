/**
 * @file output_coverage_test.cpp
 * @brief Unit tests for src/output_coverage.cpp.
 *
 * Two suites, one per public TSV-writing function. Both write a fixed,
 * deterministic tab-separated format, so tests compare exact file content
 * rather than substrings.
 *
 * The double-buffered background-thread writer (DoubleBufferedTsvWriter) is
 * an implementation detail declared only in output_coverage.cpp (external
 * linkage, but no header declares it), so it cannot be instantiated or
 * driven directly from this file - see the discussion in this session for
 * why. It is instead exercised indirectly, end-to-end, through the public
 * functions: LargeVolumeCrossesBufferBoundary below writes enough rows to
 * force multiple buffer submissions/flips through the background writer
 * thread and verifies the output is complete and uncorrupted, not just
 * "didn't crash".
 */

#include "../src/output_coverage.h"
#include "../src/config.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {
    std::filesystem::path tempTsvPath(const std::string& test_name) {
        return std::filesystem::temp_directory_path() / ("cdx_output_coverage_test_" + test_name + ".tsv");
    }

    void removeIfExists(const std::filesystem::path& path) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    std::string readFile(const std::filesystem::path& path) {
        std::ifstream in(path, std::ios::binary);
        std::ostringstream buffer;
        buffer << in.rdbuf();
        return buffer.str();
    }

    std::size_t countLines(const std::string& content) {
        return static_cast<std::size_t>(std::count(content.begin(), content.end(), '\n'));
    }
}

// =============================================================================
// writeCoverageTsvQuery
//
// Writes one TSV row per base-pair position for a single component:
// "component_name\tposition\tcoverage\n", skipping positions whose
// coverage is >= cfg::NOT_IN_QUERY. Position numbers are the raw index
// into bp_cov_table (skipped positions leave gaps, not renumbering).
// =============================================================================

// Header plus one exact row per position, in order.
TEST(WriteCoverageTsvQueryTest, HeaderAndRowsMatchInputExactly) {
    const std::vector<cdx::Coverage> coverage = {10, 20, 30};
    const auto path = tempTsvPath("query_exact_rows");

    output::writeCoverageTsvQuery(path, coverage, "chr1");
    const std::string content = readFile(path);

    const std::string expected =
        "component_name\tposition\tcoverage\n"
        "chr1\t0\t10\n"
        "chr1\t1\t20\n"
        "chr1\t2\t30\n";
    EXPECT_EQ(content, expected);
    removeIfExists(path);
}

// Sentinel positions are omitted, but surviving rows keep their original
// (non-renumbered) index.
TEST(WriteCoverageTsvQueryTest, SentinelPositionsAreSkippedWithoutRenumbering) {
    const std::vector<cdx::Coverage> coverage = {10, cfg::NOT_IN_QUERY, 20, cfg::NOT_IN_COMPO, 30};
    const auto path = tempTsvPath("query_sentinels_skipped");

    output::writeCoverageTsvQuery(path, coverage, "chr1");
    const std::string content = readFile(path);

    const std::string expected =
        "component_name\tposition\tcoverage\n"
        "chr1\t0\t10\n"
        "chr1\t2\t20\n"
        "chr1\t4\t30\n";
    EXPECT_EQ(content, expected);
    removeIfExists(path);
}

// The largest valid value, one below the sentinel range, is still written.
TEST(WriteCoverageTsvQueryTest, MaxValidValueBelowSentinelIncluded) {
    const std::vector<cdx::Coverage> coverage = {cfg::NOT_IN_QUERY - 1};
    const auto path = tempTsvPath("query_max_valid_value");

    output::writeCoverageTsvQuery(path, coverage, "chr1");
    const std::string content = readFile(path);

    const std::string expected =
        "component_name\tposition\tcoverage\n"
        "chr1\t0\t" + std::to_string(cfg::NOT_IN_QUERY - 1) + "\n";
    EXPECT_EQ(content, expected);
    removeIfExists(path);
}

// An empty component name is rejected before any file is touched.
TEST(WriteCoverageTsvQueryTest, EmptyComponentNameThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const auto path = tempTsvPath("query_empty_name");

    EXPECT_THROW(output::writeCoverageTsvQuery(path, coverage, ""), std::invalid_argument);
    removeIfExists(path);
}

// An empty coverage table still produces a valid file: header only.
TEST(WriteCoverageTsvQueryTest, EmptyCoverageTableWritesHeaderOnly) {
    const std::vector<cdx::Coverage> coverage;
    const auto path = tempTsvPath("query_empty_coverage");

    output::writeCoverageTsvQuery(path, coverage, "chr1");
    const std::string content = readFile(path);

    EXPECT_EQ(content, "component_name\tposition\tcoverage\n");
    removeIfExists(path);
}

// A coverage table made entirely of sentinel values behaves like an empty
// one from the writer's perspective: header only, no rows, no error -
// distinct from an actually-empty vector (there are elements, just none
// of them valid).
TEST(WriteCoverageTsvQueryTest, AllSentinelCoverageWritesHeaderOnly) {
    const std::vector<cdx::Coverage> coverage = {cfg::NOT_IN_QUERY, cfg::NOT_IN_COMPO, cfg::NOT_IN_QUERY};
    const auto path = tempTsvPath("query_all_sentinel");

    output::writeCoverageTsvQuery(path, coverage, "chr1");
    const std::string content = readFile(path);

    EXPECT_EQ(content, "component_name\tposition\tcoverage\n");
    removeIfExists(path);
}

// A single row whose component-name prefix alone exceeds the writer's
// internal buffer capacity cannot be formatted and must raise
// std::length_error, rather than silently truncating or corrupting output.
TEST(WriteCoverageTsvQueryTest, ComponentNameTooLongThrowsLengthError) {
    const std::string huge_name(5 * 1024 * 1024, 'x'); // > default 4MB buffer
    const std::vector<cdx::Coverage> coverage = {1};
    const auto path = tempTsvPath("query_name_too_long");

    EXPECT_THROW(output::writeCoverageTsvQuery(path, coverage, huge_name), std::length_error);
    removeIfExists(path);
}

// An unopenable output path raises std::runtime_error.
TEST(WriteCoverageTsvQueryTest, UnopenableOutputPathThrowsRuntimeError) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3};
    const std::filesystem::path bad_path = "/cdx_this_directory_should_not_exist_xyz/report.tsv";

    EXPECT_THROW(output::writeCoverageTsvQuery(bad_path, coverage, "chr1"), std::runtime_error);
}

// A large enough row count to force multiple internal buffer
// submissions/flips through the background writer thread: verifies
// completeness (no lost/duplicated rows) and correctness at both ends of
// the file, not just "the process didn't crash".
TEST(WriteCoverageTsvQueryTest, LargeVolumeCrossesBufferBoundary) {
    constexpr std::size_t kRowCount = 700'000; // several multiples of the 4MB default buffer
    std::vector<cdx::Coverage> coverage(kRowCount);
    for (std::size_t i = 0; i < kRowCount; ++i) {
        coverage[i] = static_cast<cdx::Coverage>(i % 1000);
    }
    const auto path = tempTsvPath("query_large_volume");

    output::writeCoverageTsvQuery(path, coverage, "big");
    const std::string content = readFile(path);

    EXPECT_EQ(countLines(content), kRowCount + 1); // + header
    EXPECT_EQ(content.substr(0, 33), "component_name\tposition\tcoverage\n");
    EXPECT_NE(content.find("big\t0\t0\n"), std::string::npos);
    EXPECT_NE(content.find("big\t699999\t999\n"), std::string::npos);
    // The file must end with exactly this last row (no trailing garbage,
    // no partially-written final buffer).
    const std::string last_row = "big\t699999\t999\n";
    EXPECT_EQ(content.compare(content.size() - last_row.size(), last_row.size(), last_row), 0);
    removeIfExists(path);
}

// =============================================================================
// writeCoverageTsvGlobal
//
// Writes rows for every component in a flattened multi-component coverage
// table, converting flattened positions to component-relative coordinates
// (position resets to 0 at the start of each component). Validates that
// bp_component_offsets/component_names form a complete, well-formed
// partition of the flattened table before writing anything.
// =============================================================================

// Position resets to 0 at each component boundary (component-relative,
// not the flattened/global coordinate).
TEST(WriteCoverageTsvGlobalTest, PositionResetsPerComponent) {
    const std::vector<cdx::Coverage> coverage = {10, 20, 30, 40, 50};
    const std::vector<cdx::PosBp> offsets = {0, 2, 5};
    const std::vector<std::string> names = {"chr1", "chr2"};
    const auto path = tempTsvPath("global_position_reset");

    output::writeCoverageTsvGlobal(path, coverage, offsets, names);
    const std::string content = readFile(path);

    const std::string expected =
        "component_name\tposition\tcoverage\n"
        "chr1\t0\t10\n"
        "chr1\t1\t20\n"
        "chr2\t0\t30\n"
        "chr2\t1\t40\n"
        "chr2\t2\t50\n";
    EXPECT_EQ(content, expected);
    removeIfExists(path);
}

// A zero-length component (two equal consecutive offsets) contributes no
// rows at all, but does not break the components after it.
TEST(WriteCoverageTsvGlobalTest, ZeroLengthComponentProducesNoRows) {
    const std::vector<cdx::Coverage> coverage = {10, 20, 30};
    const std::vector<cdx::PosBp> offsets = {0, 0, 3};
    const std::vector<std::string> names = {"empty", "chr1"};
    const auto path = tempTsvPath("global_zero_length_component");

    output::writeCoverageTsvGlobal(path, coverage, offsets, names);
    const std::string content = readFile(path);

    const std::string expected =
        "component_name\tposition\tcoverage\n"
        "chr1\t0\t10\n"
        "chr1\t1\t20\n"
        "chr1\t2\t30\n";
    EXPECT_EQ(content, expected);
    removeIfExists(path);
}

// Sentinel values are skipped independently within each component, again
// without renumbering the surviving positions.
TEST(WriteCoverageTsvGlobalTest, SentinelValuesSkippedAcrossComponents) {
    const std::vector<cdx::Coverage> coverage = {10, cfg::NOT_IN_QUERY, 20, cfg::NOT_IN_COMPO, 30};
    const std::vector<cdx::PosBp> offsets = {0, 2, 5};
    const std::vector<std::string> names = {"chr1", "chr2"};
    const auto path = tempTsvPath("global_sentinels_skipped");

    output::writeCoverageTsvGlobal(path, coverage, offsets, names);
    const std::string content = readFile(path);

    const std::string expected =
        "component_name\tposition\tcoverage\n"
        "chr1\t0\t10\n"
        "chr2\t0\t20\n"
        "chr2\t2\t30\n";
    EXPECT_EQ(content, expected);
    removeIfExists(path);
}

// A single-component graph is the degenerate case of the multi-component
// writer; rows and positions must still be correct.
TEST(WriteCoverageTsvGlobalTest, SingleComponentWritesRelativePositions) {
    const std::vector<cdx::Coverage> coverage = {10, 20, 30};
    const std::vector<cdx::PosBp> offsets = {0, 3};
    const std::vector<std::string> names = {"chr1"};
    const auto path = tempTsvPath("global_single_component");

    output::writeCoverageTsvGlobal(path, coverage, offsets, names);
    const std::string content = readFile(path);

    const std::string expected =
        "component_name\tposition\tcoverage\n"
        "chr1\t0\t10\n"
        "chr1\t1\t20\n"
        "chr1\t2\t30\n";
    EXPECT_EQ(content, expected);
    removeIfExists(path);
}

// All-sentinel coverage across every component: header only, no rows.
TEST(WriteCoverageTsvGlobalTest, AllSentinelCoverageWritesHeaderOnly) {
    const std::vector<cdx::Coverage> coverage = {cfg::NOT_IN_QUERY, cfg::NOT_IN_COMPO, cfg::NOT_IN_QUERY};
    const std::vector<cdx::PosBp> offsets = {0, 1, 3};
    const std::vector<std::string> names = {"chr1", "chr2"};
    const auto path = tempTsvPath("global_all_sentinel");

    output::writeCoverageTsvGlobal(path, coverage, offsets, names);
    const std::string content = readFile(path);

    EXPECT_EQ(content, "component_name\tposition\tcoverage\n");
    removeIfExists(path);
}

// bp_component_offsets must contain at least two boundaries.
TEST(WriteCoverageTsvGlobalTest, ComponentOffsetsMustContainTwoBoundariesThrows) {
    const std::vector<cdx::Coverage> coverage;
    const std::vector<cdx::PosBp> offsets = {0};
    const std::vector<std::string> names;
    const auto path = tempTsvPath("global_offsets_too_short");

    EXPECT_THROW(output::writeCoverageTsvGlobal(path, coverage, offsets, names), std::invalid_argument);
    removeIfExists(path);
}

// component_names.size() must equal bp_component_offsets.size() - 1.
TEST(WriteCoverageTsvGlobalTest, ComponentNameCountMismatchThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 2, 4};
    const std::vector<std::string> names = {"chr1"}; // should have 2 entries
    const auto path = tempTsvPath("global_name_count_mismatch");

    EXPECT_THROW(output::writeCoverageTsvGlobal(path, coverage, offsets, names), std::invalid_argument);
    removeIfExists(path);
}

// The first offset must be zero.
TEST(WriteCoverageTsvGlobalTest, FirstOffsetMustBeZeroThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {1, 4};
    const std::vector<std::string> names = {"chr1"};
    const auto path = tempTsvPath("global_first_offset_nonzero");

    EXPECT_THROW(output::writeCoverageTsvGlobal(path, coverage, offsets, names), std::invalid_argument);
    removeIfExists(path);
}

// The final offset must match the flattened coverage-table size.
TEST(WriteCoverageTsvGlobalTest, FinalOffsetMustMatchTableSizeThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 99};
    const std::vector<std::string> names = {"chr1"};
    const auto path = tempTsvPath("global_final_offset_mismatch");

    EXPECT_THROW(output::writeCoverageTsvGlobal(path, coverage, offsets, names), std::invalid_argument);
    removeIfExists(path);
}

// Offsets must be non-decreasing.
TEST(WriteCoverageTsvGlobalTest, OffsetsMustBeNonDecreasingThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 3, 2, 4};
    const std::vector<std::string> names = {"chr1", "chr2"};
    const auto path = tempTsvPath("global_offsets_not_nondecreasing");

    EXPECT_THROW(output::writeCoverageTsvGlobal(path, coverage, offsets, names), std::invalid_argument);
    removeIfExists(path);
}

// Every component name must be non-empty, since it is written verbatim
// into the TSV.
TEST(WriteCoverageTsvGlobalTest, EmptyComponentNameThrows) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 2, 4};
    const std::vector<std::string> names = {"chr1", ""};
    const auto path = tempTsvPath("global_empty_component_name");

    EXPECT_THROW(output::writeCoverageTsvGlobal(path, coverage, offsets, names), std::invalid_argument);
    removeIfExists(path);
}

// An unopenable output path raises std::runtime_error.
TEST(WriteCoverageTsvGlobalTest, UnopenableOutputPathThrowsRuntimeError) {
    const std::vector<cdx::Coverage> coverage = {1, 2, 3, 4};
    const std::vector<cdx::PosBp> offsets = {0, 4};
    const std::vector<std::string> names = {"chr1"};
    const std::filesystem::path bad_path = "/cdx_this_directory_should_not_exist_xyz/report.tsv";

    EXPECT_THROW(output::writeCoverageTsvGlobal(bad_path, coverage, offsets, names), std::runtime_error);
}
