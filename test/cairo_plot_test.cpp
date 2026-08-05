/**
 * @file cairo_plot_test.cpp
 * @brief Light integration tests for cairo_plot.cpp (in-process libcairo PNG rendering).
 *
 * cairo_plot.cpp's internals (parseHexColor, niceAxisTicks, drawPanel, buildPanelSpec,
 * writePngOrThrow) all live in an anonymous namespace, so only the two public entry points
 * declared in output_plot.h - writeLinearPlotQuery() and writeLinearPlotGlobal() - are
 * reachable from a test. These are therefore "light" integration tests: rather than
 * verifying pixel content (which would require decoding the PNG image data), each test
 * checks that:
 *   - a syntactically valid PNG file is produced (correct 8-byte signature + a well-formed
 *     IHDR chunk right after it),
 *   - the PNG's declared width/height (read directly out of the IHDR chunk, big-endian,
 *     bytes [16,20) and [20,24)) match what the figure_width/figure_height/dpi/grid
 *     configuration should produce,
 *   - error paths (bad component offsets, undependable output directories) throw the
 *     expected exception type instead of crashing or silently producing a broken file.
 *
 * NOTE: this file could not be compiled or executed in the sandbox this was written in
 * (no libcairo available there), unlike linear_plot_test.cpp. It has been written and
 * reviewed carefully against the real cairo_plot.cpp/output_plot.h sources, but has not
 * been locally verified - please build it via `ctest` (or run the cairo_plot_test binary
 * directly) on your machine and report back any failures.
 *
 * Because cairo_plot.cpp calls prepareLinearPlotPackage()/chooseGlobalGraphGrid() (defined
 * in linear_plot.cpp, declared in output_plot.h), this test target compiles both .cpp files
 * together - see the CMakeLists.txt target below.
 */

#include <gtest/gtest.h>

#include "../src/output_plot.h"
#include "../src/config.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {
    namespace fs = std::filesystem;

    // ---------------------------------------------------------------------
    // Test fixtures / helpers
    // ---------------------------------------------------------------------

    /** @brief RAII temporary directory, recursively removed on destruction. */
    class TempDir {
    public:
        TempDir() {
            path_ = fs::temp_directory_path() /
                    ("cairo_plot_test_" + std::to_string(reinterpret_cast<std::uintptr_t>(this)));
            fs::create_directories(path_);
        }

        ~TempDir() {
            std::error_code ec;
            fs::remove_all(path_, ec);
        }

        TempDir(const TempDir &) = delete;

        TempDir &operator=(const TempDir &) = delete;

        [[nodiscard]] const fs::path &path() const { return path_; }

    private:
        fs::path path_;
    };

    /** @brief Width/height declared in a PNG file's IHDR chunk. */
    struct PngDims {
        std::uint32_t width = 0;
        std::uint32_t height = 0;
    };

    /**
     * @brief Reads and validates the PNG signature + IHDR chunk header of a file on disk,
     *        returning the width/height it declares.
     *
     * Fails the current test (via GTEST_FAIL / ADD_FAILURE, through the caller's EXPECT/ASSERT
     * usage) rather than throwing, so callers should check the returned optional-like validity
     * via the companion isValidPng() before trusting the dimensions, or just call this after
     * already asserting isValidPng().
     */
    PngDims readPngDims(const fs::path &png_path) {
        std::ifstream file(png_path, std::ios::binary);
        EXPECT_TRUE(file.is_open()) << "Could not open PNG file: " << png_path;

        unsigned char header[24] = {};
        file.read(reinterpret_cast<char *>(header), sizeof(header));
        EXPECT_EQ(file.gcount(), static_cast<std::streamsize>(sizeof(header)))
            << "PNG file shorter than the signature + IHDR header: " << png_path;

        static constexpr unsigned char kPngSignature[8] = {
            0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A
        };
        for (int i = 0; i < 8; ++i) {
            EXPECT_EQ(header[i], kPngSignature[i]) << "Invalid PNG signature byte at index " << i;
        }

        // Bytes [8,12) = chunk length (should be 13 for IHDR), [12,16) = "IHDR" tag.
        EXPECT_EQ(header[12], 'I');
        EXPECT_EQ(header[13], 'H');
        EXPECT_EQ(header[14], 'D');
        EXPECT_EQ(header[15], 'R');

        PngDims dims;
        dims.width = (static_cast<std::uint32_t>(header[16]) << 24) |
                     (static_cast<std::uint32_t>(header[17]) << 16) |
                     (static_cast<std::uint32_t>(header[18]) << 8) |
                     static_cast<std::uint32_t>(header[19]);
        dims.height = (static_cast<std::uint32_t>(header[20]) << 24) |
                      (static_cast<std::uint32_t>(header[21]) << 16) |
                      (static_cast<std::uint32_t>(header[22]) << 8) |
                      static_cast<std::uint32_t>(header[23]);
        return dims;
    }

    /** @brief Builds a plausible-looking coverage vector of length @p n (simple ramp/wave). */
    std::vector<cdx::Coverage> makeCoverage(const std::size_t n) {
        std::vector<cdx::Coverage> coverage;
        coverage.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            coverage.push_back(static_cast<cdx::Coverage>((i * 7) % 23));
        }
        return coverage;
    }

    /** @brief A minimal, otherwise-default PlotConfig with small dimensions/dpi for fast tests. */
    output::PlotConfig smallConfig() {
        output::PlotConfig config;
        config.dpi = 72;
        config.figure_width = 3.0;
        config.figure_height = 2.0;
        return config;
    }
} // namespace

// ===========================================================================
// writeLinearPlotQuery
// ===========================================================================

TEST(WriteLinearPlotQueryTest, ValidCoverageProducesPngWithExpectedDimensions) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    output::writeLinearPlotQuery(out, makeCoverage(200), "chr1", 0, config);

    ASSERT_TRUE(fs::exists(out));
    EXPECT_GT(fs::file_size(out), 0u);

    const PngDims dims = readPngDims(out);
    const auto expected_w = static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi));
    const auto expected_h = static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi));
    EXPECT_EQ(dims.width, expected_w);
    EXPECT_EQ(dims.height, expected_h);
}

TEST(WriteLinearPlotQueryTest, DifferentDpiAndFigureSizeChangeDimensions) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";

    output::PlotConfig config;
    config.dpi = 150;
    config.figure_width = 5.0;
    config.figure_height = 2.5;

    output::writeLinearPlotQuery(out, makeCoverage(50), "chr2", 1000, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);
    EXPECT_EQ(dims.width, static_cast<std::uint32_t>(std::llround(5.0 * 150)));
    EXPECT_EQ(dims.height, static_cast<std::uint32_t>(std::llround(2.5 * 150)));
}

TEST(WriteLinearPlotQueryTest, EmptyCoverageStillProducesValidPng) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    output::writeLinearPlotQuery(out, {}, "empty_region", 0, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);
    EXPECT_EQ(dims.width, static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi)));
    EXPECT_EQ(dims.height, static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi)));
}

TEST(WriteLinearPlotQueryTest, SentinelValuesAreFilteredWithoutCrashing) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    std::vector<cdx::Coverage> coverage{5, 10, cfg::NOT_IN_QUERY, 15, cfg::NOT_IN_COMPO, 20, 0, 3};

    output::writeLinearPlotQuery(out, coverage, "sentinel_region", 0, config);

    ASSERT_TRUE(fs::exists(out));
    EXPECT_GT(fs::file_size(out), 0u);
    (void) readPngDims(out); // just validate it's well-formed
}

TEST(WriteLinearPlotQueryTest, LogarithmicScaleProducesValidPng) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    auto config = smallConfig();
    config.log_base = 2;

    output::writeLinearPlotQuery(out, makeCoverage(300), "log_region", 0, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);
    EXPECT_EQ(dims.width, static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi)));
    EXPECT_EQ(dims.height, static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi)));
}

TEST(WriteLinearPlotQueryTest, NestedOutputDirectoryIsCreatedAutomatically) {
    TempDir dir;
    const auto out = dir.path() / "a" / "b" / "c" / "plot.png";
    ASSERT_FALSE(fs::exists(out.parent_path()));

    output::writeLinearPlotQuery(out, makeCoverage(20), "nested", 0, smallConfig());

    EXPECT_TRUE(fs::exists(out));
}

TEST(WriteLinearPlotQueryTest, UncreatableParentDirectoryThrows) {
    TempDir dir;
    const auto blocker = dir.path() / "blocker";
    {
        std::ofstream blocker_file(blocker);
        blocker_file << "not a directory";
    }
    // "blocker" exists as a regular file, so create_directories() on
    // "blocker/plot.png"'s parent ("blocker") must fail.
    const auto out = blocker / "plot.png";

    EXPECT_THROW(
        output::writeLinearPlotQuery(out, makeCoverage(10), "x", 0, smallConfig()),
        fs::filesystem_error
    );
}

// ===========================================================================
// writeLinearPlotGlobal
// ===========================================================================

TEST(WriteLinearPlotGlobalTest, FewerThanTwoOffsetsThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const std::vector<cdx::Coverage> coverage(10, 5);
    const std::vector<cdx::PosBp> offsets{0}; // only one boundary => zero components
    const std::vector<std::string> names{"chr1"};

    EXPECT_THROW(
        output::writeLinearPlotGlobal(out, coverage, offsets, names, smallConfig()),
        std::invalid_argument
    );
}

TEST(WriteLinearPlotGlobalTest, EmptyOffsetsThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const std::vector<cdx::Coverage> coverage(10, 5);

    EXPECT_THROW(
        output::writeLinearPlotGlobal(out, coverage, {}, {}, smallConfig()),
        std::invalid_argument
    );
}

TEST(WriteLinearPlotGlobalTest, EndBeforeStartOffsetThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const std::vector<cdx::Coverage> coverage(10, 5);
    // A single component whose end (3) is before its start (5).
    const std::vector<cdx::PosBp> offsets{5, 3};
    const std::vector<std::string> names{"chr1"};

    EXPECT_THROW(
        output::writeLinearPlotGlobal(out, coverage, offsets, names, smallConfig()),
        std::invalid_argument
    );
}

TEST(WriteLinearPlotGlobalTest, EndBeyondCoverageSizeThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const std::vector<cdx::Coverage> coverage(10, 5);
    // A single component whose end (100) is beyond flat_bp_coverage.size() (10).
    const std::vector<cdx::PosBp> offsets{0, 100};
    const std::vector<std::string> names{"chr1"};

    EXPECT_THROW(
        output::writeLinearPlotGlobal(out, coverage, offsets, names, smallConfig()),
        std::invalid_argument
    );
}

TEST(WriteLinearPlotGlobalTest, SingleComponentProducesOnePanelGrid) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    const auto coverage = makeCoverage(120);
    const std::vector<cdx::PosBp> offsets{0, 120};
    const std::vector<std::string> names{"chr1"};

    output::writeLinearPlotGlobal(out, coverage, offsets, names, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);

    const auto grid = output::chooseGlobalGraphGrid(1);
    const auto panel_w = static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi));
    const auto panel_h = static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi));

    EXPECT_EQ(dims.width, panel_w * static_cast<std::uint32_t>(grid.columns));
    EXPECT_EQ(dims.height, panel_h * static_cast<std::uint32_t>(grid.rows));
}

TEST(WriteLinearPlotGlobalTest, MultipleComponentsMatchChooseGlobalGraphGridLayout) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    // 4 components of 50 bp each, contiguous in the flattened coverage vector.
    constexpr std::size_t per_component = 50;
    constexpr std::size_t component_count = 4;
    const auto coverage = makeCoverage(per_component * component_count);

    std::vector<cdx::PosBp> offsets;
    for (std::size_t i = 0; i <= component_count; ++i) {
        offsets.push_back(static_cast<cdx::PosBp>(i * per_component));
    }
    const std::vector<std::string> names{"chr1", "chr2", "chr3", "chr4"};

    output::writeLinearPlotGlobal(out, coverage, offsets, names, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);

    const auto grid = output::chooseGlobalGraphGrid(component_count);
    const auto panel_w = static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi));
    const auto panel_h = static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi));

    EXPECT_EQ(dims.width, panel_w * static_cast<std::uint32_t>(grid.columns));
    EXPECT_EQ(dims.height, panel_h * static_cast<std::uint32_t>(grid.rows));
}

TEST(WriteLinearPlotGlobalTest, TwelveComponentsMatchChooseGlobalGraphGridLayout) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    constexpr std::size_t per_component = 30;
    constexpr std::size_t component_count = 12;
    const auto coverage = makeCoverage(per_component * component_count);

    std::vector<cdx::PosBp> offsets;
    for (std::size_t i = 0; i <= component_count; ++i) {
        offsets.push_back(static_cast<cdx::PosBp>(i * per_component));
    }
    std::vector<std::string> names;
    for (std::size_t i = 0; i < component_count; ++i) {
        names.push_back("chr" + std::to_string(i + 1));
    }

    output::writeLinearPlotGlobal(out, coverage, offsets, names, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);

    const auto grid = output::chooseGlobalGraphGrid(component_count);
    const auto panel_w = static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi));
    const auto panel_h = static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi));

    EXPECT_EQ(dims.width, panel_w * static_cast<std::uint32_t>(grid.columns));
    EXPECT_EQ(dims.height, panel_h * static_cast<std::uint32_t>(grid.rows));
}

TEST(WriteLinearPlotGlobalTest, ComponentNamesShorterThanCountFallsBackToIndex) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    constexpr std::size_t per_component = 40;
    constexpr std::size_t component_count = 3;
    const auto coverage = makeCoverage(per_component * component_count);

    std::vector<cdx::PosBp> offsets;
    for (std::size_t i = 0; i <= component_count; ++i) {
        offsets.push_back(static_cast<cdx::PosBp>(i * per_component));
    }
    // Only one name for three components: components 1 and 2 must fall back
    // to std::to_string(component_id) internally instead of crashing.
    const std::vector<std::string> names{"only_one_name"};

    output::writeLinearPlotGlobal(out, coverage, offsets, names, config);

    ASSERT_TRUE(fs::exists(out));
    EXPECT_GT(fs::file_size(out), 0u);
}

TEST(WriteLinearPlotGlobalTest, EmptyComponentNamesFallsBackToIndexForAll) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    constexpr std::size_t per_component = 25;
    constexpr std::size_t component_count = 2;
    const auto coverage = makeCoverage(per_component * component_count);

    std::vector<cdx::PosBp> offsets{0, per_component, per_component * 2};

    output::writeLinearPlotGlobal(out, coverage, offsets, {}, config);

    ASSERT_TRUE(fs::exists(out));
    EXPECT_GT(fs::file_size(out), 0u);
}

TEST(WriteLinearPlotGlobalTest, LogarithmicScaleInGlobalModeProducesValidPng) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    auto config = smallConfig();
    config.log_base = 10;

    constexpr std::size_t per_component = 60;
    constexpr std::size_t component_count = 3;
    const auto coverage = makeCoverage(per_component * component_count);

    std::vector<cdx::PosBp> offsets;
    for (std::size_t i = 0; i <= component_count; ++i) {
        offsets.push_back(static_cast<cdx::PosBp>(i * per_component));
    }
    const std::vector<std::string> names{"a", "b", "c"};

    output::writeLinearPlotGlobal(out, coverage, offsets, names, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);

    const auto grid = output::chooseGlobalGraphGrid(component_count);
    const auto panel_w = static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi));
    const auto panel_h = static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi));

    EXPECT_EQ(dims.width, panel_w * static_cast<std::uint32_t>(grid.columns));
    EXPECT_EQ(dims.height, panel_h * static_cast<std::uint32_t>(grid.rows));
}

TEST(WriteLinearPlotGlobalTest, NestedOutputDirectoryIsCreatedAutomatically) {
    TempDir dir;
    const auto out = dir.path() / "nested" / "dir" / "plot.png";
    ASSERT_FALSE(fs::exists(out.parent_path()));

    const auto coverage = makeCoverage(60);
    const std::vector<cdx::PosBp> offsets{0, 60};
    const std::vector<std::string> names{"chr1"};

    output::writeLinearPlotGlobal(out, coverage, offsets, names, smallConfig());

    EXPECT_TRUE(fs::exists(out));
}

TEST(WriteLinearPlotGlobalTest, ZeroLengthComponentIsHandledWithoutCrashing) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto config = smallConfig();

    // Component 0: 40 bp. Component 1: zero-length (start == end).
    const auto coverage = makeCoverage(40);
    const std::vector<cdx::PosBp> offsets{0, 40, 40};
    const std::vector<std::string> names{"chr1", "chr2_empty"};

    output::writeLinearPlotGlobal(out, coverage, offsets, names, config);

    ASSERT_TRUE(fs::exists(out));
    const PngDims dims = readPngDims(out);

    const auto grid = output::chooseGlobalGraphGrid(2);
    const auto panel_w = static_cast<std::uint32_t>(std::llround(config.figure_width * config.dpi));
    const auto panel_h = static_cast<std::uint32_t>(std::llround(config.figure_height * config.dpi));

    EXPECT_EQ(dims.width, panel_w * static_cast<std::uint32_t>(grid.columns));
    EXPECT_EQ(dims.height, panel_h * static_cast<std::uint32_t>(grid.rows));
}
