/**
 * @file cli_test.cpp
 * @brief Unit tests for cli.cpp (parse_args and its internal parsing helpers).
 *
 * parse_args() is exercised entirely through its public, real CLI11-based entry point, using
 * synthetic argv arrays and small placeholder files on disk (CLI::ExistingFile only checks that
 * the path exists, so empty files are sufficient fixtures for the mandatory `cdx`/`gam`
 * positionals). This indirectly covers every internal anonymous-namespace helper
 * (parse_query, parse_component_type, parse_fig_size, parse_hex_color, trim).
 *
 * As agreed, every CLI option is tested with both valid values and inputs designed to force a
 * parsing/validation failure (wrong type, out of range, malformed strings, etc.) - this includes
 * options bound directly to a typed CLI11 value (--dpi is `int`, --smoothing is `double`,
 * --max-points is `size_t`, --log is `int`) and options that are parsed manually as strings
 * after CLI11 hands them over as-is (-t/-T, --query, --component-type, --fig-size, the colors).
 *
 * parse_args() was refactored (see cli.cpp) so that every parsing/validation failure - both
 * CLI11-level (CLI::ParseError, e.g. a non-numeric value for a typed int/double option) and the
 * hand-written post-parse checks - now throws std::runtime_error instead of calling std::exit(1)
 * directly, specifically so this module could be unit tested in-process. The sole remaining
 * std::exit() call is for --help (exit code 0, not an error), which two dedicated death tests at
 * the bottom of this file cover with EXPECT_EXIT.
 *
 * NOTE ON VERIFICATION: like gam_io_test.cpp, this file could not be locally compiled in the
 * sandbox used to author it - CLI11 is not vendored in this repository (`include/` is empty here)
 * and is not fetched by CMake, so there was nothing to build against. It was written from a
 * careful line-by-line reading of cli.cpp/cli.hpp and CLI11's documented option semantics
 * (required(), ->check(CLI::ExistingFile), ->type_size(0, 1), ->default_str(), ->check(CLI::Range),
 * ->check(CLI::PositiveNumber)) rather than the usual compile+mutate verification loop. Please
 * build and run it and report back anything that doesn't match.
 */

#include "../src/cli.hpp"

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
    // =====================================================================
    // Test infrastructure.
    // =====================================================================

    /** @brief Creates an empty placeholder file (only its existence matters to CLI::ExistingFile). */
    std::filesystem::path makeTempFile(const std::string &suffix) {
        static int counter = 0;
        const std::filesystem::path path =
                std::filesystem::temp_directory_path() / ("cli_test_" + std::to_string(++counter) + suffix);
        std::ofstream(path).close();
        return path;
    }

    /** @brief Shared placeholder cdx/gam files, created once and removed at process exit. */
    struct FixturePaths {
        std::filesystem::path cdx = makeTempFile(".cdx");
        std::filesystem::path gam = makeTempFile(".gam");

        ~FixturePaths() {
            std::error_code ec;
            std::filesystem::remove(cdx, ec);
            std::filesystem::remove(gam, ec);
        }
    };

    const FixturePaths &fixtures() {
        static const FixturePaths paths;
        return paths;
    }

    /**
     * @brief Builds a stable argv array (argv[0] plus every supplied token) and calls parse_args.
     *        Owns the underlying string storage for the lifetime of the call.
     */
    CliArgs parseWith(const std::vector<std::string> &tokens) {
        std::vector<std::string> storage;
        storage.emplace_back("cdx_coverage"); // argv[0]
        storage.insert(storage.end(), tokens.begin(), tokens.end());

        std::vector<char *> argv;
        argv.reserve(storage.size());
        for (auto &token: storage) {
            argv.push_back(token.data());
        }

        return parse_args(static_cast<int>(argv.size()), argv.data());
    }

    /** @brief Convenience: cdx + gam + extra tokens (the common non-inspect case). */
    CliArgs parseWithCdxGam(std::initializer_list<std::string> extra = {}) {
        std::vector<std::string> tokens{fixtures().cdx.string(), fixtures().gam.string()};
        tokens.insert(tokens.end(), extra.begin(), extra.end());
        return parseWith(tokens);
    }

    /** @brief Convenience: cdx only + extra tokens (for inspect-mode tests, no gam required). */
    CliArgs parseWithCdxOnly(std::initializer_list<std::string> extra = {}) {
        std::vector<std::string> tokens{fixtures().cdx.string()};
        tokens.insert(tokens.end(), extra.begin(), extra.end());
        return parseWith(tokens);
    }
} // anonymous namespace

// =============================================================================
// Basic positional arguments and required-file handling.
// =============================================================================
namespace {
    TEST(CliBasicArgsTest, MinimalValidArgsPopulateDefaults) {
        const CliArgs args = parseWithCdxGam();

        EXPECT_EQ(args.cdx_file, fixtures().cdx.string());
        EXPECT_EQ(args.gam_file, fixtures().gam.string());
        EXPECT_FALSE(args.query.has_value());
        EXPECT_EQ(args.component_type, ComponentType::Linear);
        EXPECT_FALSE(args.inspectMode());
        EXPECT_EQ(args.output_directory, ".");
        EXPECT_TRUE(args.generateGraph());
        EXPECT_TRUE(args.generateStats());
        EXPECT_TRUE(args.generateTable());
        EXPECT_GE(args.worker_threads, 1);
        EXPECT_GE(args.decompression_threads, 1);
        EXPECT_DOUBLE_EQ(args.smoothing, 0.01);
        EXPECT_EQ(args.max_plot_points, 10000u);
        EXPECT_EQ(args.dpi, 300);
        EXPECT_FALSE(args.log_base.has_value());
        EXPECT_EQ(args.line_color, "#1E3A8A");
        EXPECT_EQ(args.fill_color, "#93C5FD");
        // No query, linear -> standard overview canvas
        EXPECT_DOUBLE_EQ(args.figure_width, 5.5);
        EXPECT_DOUBLE_EQ(args.figure_height, 3.5);
    }

    TEST(CliBasicArgsTest, MissingRequiredCdxPositionalThrows) {
        EXPECT_THROW(parseWith({}), std::runtime_error);
    }

    TEST(CliBasicArgsTest, NonExistentCdxFileThrows) {
        const std::filesystem::path missing = std::filesystem::temp_directory_path() / "does_not_exist.cdx";
        EXPECT_THROW(parseWith({missing.string(), fixtures().gam.string()}), std::runtime_error);
    }

    TEST(CliBasicArgsTest, NonExistentGamFileThrows) {
        const std::filesystem::path missing = std::filesystem::temp_directory_path() / "does_not_exist.gam";
        EXPECT_THROW(parseWith({fixtures().cdx.string(), missing.string()}), std::runtime_error);
    }

    TEST(CliBasicArgsTest, MissingGamFileWithoutInspectThrows) {
        EXPECT_THROW(parseWithCdxOnly(), std::runtime_error);
    }

    TEST(CliBasicArgsTest, InspectModeDoesNotRequireGamFile) {
        const CliArgs args = parseWithCdxOnly({"-i"});
        EXPECT_TRUE(args.inspectMode());
        EXPECT_TRUE(args.gam_file.empty());
    }

    TEST(CliBasicArgsTest, UnknownOptionThrows) {
        EXPECT_THROW(parseWithCdxGam({"--this-flag-does-not-exist"}), std::runtime_error);
    }
} // anonymous namespace

// =============================================================================
// --query / -q
// =============================================================================
namespace {
    TEST(CliQueryOptionTest, ComponentNameOnly) {
        const CliArgs args = parseWithCdxGam({"-q", "chr1"});
        ASSERT_TRUE(args.query.has_value());
        EXPECT_EQ(args.query->component, "chr1");
        EXPECT_FALSE(args.query->range.has_value());
    }

    TEST(CliQueryOptionTest, ComponentIdOnly) {
        const CliArgs args = parseWithCdxGam({"-q", "0"});
        ASSERT_TRUE(args.query.has_value());
        EXPECT_EQ(args.query->component, "0");
        EXPECT_FALSE(args.query->range.has_value());
    }

    TEST(CliQueryOptionTest, ComponentWithPositiveRange) {
        const CliArgs args = parseWithCdxGam({"-q", "chr1 1000:5000"});
        ASSERT_TRUE(args.query.has_value());
        EXPECT_EQ(args.query->component, "chr1");
        ASSERT_TRUE(args.query->range.has_value());
        EXPECT_EQ(args.query->range->start, 1000);
        EXPECT_EQ(args.query->range->end, 5000);
    }

    TEST(CliQueryOptionTest, ComponentWithNegativeRange) {
        const CliArgs args = parseWithCdxGam({"-q", "chrX -10:-100"});
        ASSERT_TRUE(args.query.has_value());
        ASSERT_TRUE(args.query->range.has_value());
        EXPECT_EQ(args.query->range->start, -10);
        EXPECT_EQ(args.query->range->end, -100);
    }

    // The help text used to say "chr1 1000-5000" (dash) and "COMPONENT:START-END"
    // (colon before the range) - neither matches the actual regex, which requires
    // a space before the range and a colon between start/end. Both wrong forms
    // must be rejected.
    TEST(CliQueryOptionTest, DashSeparatedRangeThrows) {
        EXPECT_THROW(parseWithCdxGam({"-q", "chr1 1000-5000"}), std::runtime_error);
    }

    // Not an error: the regex's component group ([^\s]+) matches any run of
    // non-whitespace characters, colons included, and the range group only
    // ever activates after a *space*. With no space anywhere in the string,
    // the whole thing - colons and digits together - is swallowed as one
    // literal component name, with no range attached.
    TEST(CliQueryOptionTest, ColonBetweenComponentAndRangeIsTakenAsLiteralComponentName) {
        const CliArgs args = parseWithCdxGam({"-q", "chr1:1000:5000"});
        ASSERT_TRUE(args.query.has_value());
        EXPECT_EQ(args.query->component, "chr1:1000:5000");
        EXPECT_FALSE(args.query->range.has_value());
    }

    TEST(CliQueryOptionTest, EmptyStringIsTreatedAsNoQuery) {
        const CliArgs args = parseWithCdxGam({"-q", ""});
        EXPECT_FALSE(args.query.has_value());
    }

    TEST(CliQueryOptionTest, WhitespaceOnlyThrows) {
        EXPECT_THROW(parseWithCdxGam({"-q", "   "}), std::runtime_error);
    }
} // anonymous namespace

// =============================================================================
// --component-type / -c
// =============================================================================
namespace {
    TEST(CliComponentTypeOptionTest, DefaultIsLinear) {
        EXPECT_EQ(parseWithCdxGam().component_type, ComponentType::Linear);
    }

    TEST(CliComponentTypeOptionTest, LinearAliasesAccepted) {
        for (const std::string &value: {"l", "linear", "L", "LINEAR", "Linear"}) {
            EXPECT_EQ(parseWithCdxGam({"-c", value}).component_type, ComponentType::Linear) << "value=" << value;
        }
    }

    TEST(CliComponentTypeOptionTest, CircularAliasesAccepted) {
        for (const std::string &value: {"c", "circular", "C", "CIRCULAR", "Circular"}) {
            EXPECT_EQ(parseWithCdxGam({"-c", value}).component_type, ComponentType::Circular) << "value=" << value;
        }
    }

    TEST(CliComponentTypeOptionTest, InvalidValueThrows) {
        EXPECT_THROW(parseWithCdxGam({"-c", "triangle"}), std::runtime_error);
    }
} // anonymous namespace

// =============================================================================
// --coverage-precision / -p
//
// Mirrors CliComponentTypeOptionTest above (parse_coverage_precision() is a
// direct copy of parse_component_type()'s structure), with one addition:
// the default is Base (not the "cheap/original" option), unlike
// --component-type defaulting to Linear - see coverage_precision.h for why.
// =============================================================================
namespace {
    TEST(CliCoveragePrecisionOptionTest, DefaultIsBase) {
        EXPECT_EQ(parseWithCdxGam().coverage_precision, CoveragePrecision::Base);
    }

    TEST(CliCoveragePrecisionOptionTest, BaseAliasesAccepted) {
        for (const std::string &value: {"b", "base", "B", "BASE", "Base"}) {
            EXPECT_EQ(parseWithCdxGam({"-p", value}).coverage_precision, CoveragePrecision::Base) << "value=" << value;
        }
    }

    TEST(CliCoveragePrecisionOptionTest, NodeAliasesAccepted) {
        for (const std::string &value: {"n", "node", "N", "NODE", "Node"}) {
            EXPECT_EQ(parseWithCdxGam({"-p", value}).coverage_precision, CoveragePrecision::Node) << "value=" << value;
        }
    }

    TEST(CliCoveragePrecisionOptionTest, InvalidValueThrows) {
        EXPECT_THROW(parseWithCdxGam({"-p", "ultra"}), std::runtime_error);
    }

    TEST(CliCoveragePrecisionOptionTest, LongFormOptionNameAccepted) {
        EXPECT_EQ(parseWithCdxGam({"--coverage-precision", "node"}).coverage_precision, CoveragePrecision::Node);
    }
} // anonymous namespace

// =============================================================================
// -t/--worker-threads and -T/--decompression-threads
// =============================================================================
namespace {
    TEST(CliThreadOptionTest, AutoResolvesToAtLeastOneThread) {
        const CliArgs args = parseWithCdxGam();
        EXPECT_GE(args.worker_threads, 1);
        EXPECT_GE(args.decompression_threads, 1);
    }

    // Any sufficiently large explicit request must be capped to the same
    // machine-dependent maximum, regardless of exactly how large it is.
    TEST(CliThreadOptionTest, ExplicitWorkerThreadsCapToHardwareMaximum) {
        const CliArgs a = parseWithCdxGam({"-t", "100000"});
        const CliArgs b = parseWithCdxGam({"-t", "200000"});
        EXPECT_EQ(a.worker_threads, b.worker_threads);
        EXPECT_LT(a.worker_threads, 100000);
    }

    TEST(CliThreadOptionTest, WorkerThreadsZeroThrows) {
        EXPECT_THROW(parseWithCdxGam({"-t", "0"}), std::runtime_error);
    }

    TEST(CliThreadOptionTest, WorkerThreadsNegativeThrows) {
        EXPECT_THROW(parseWithCdxGam({"-t", "-5"}), std::runtime_error);
    }

    // -t is parsed as a raw string by CLI11 (to allow "auto"), so a non-numeric
    // value is not caught by CLI11's own type system - it fails later in the
    // manual std::stoi conversion instead. Either way it must still throw.
    TEST(CliThreadOptionTest, WorkerThreadsNonNumericStringThrows) {
        EXPECT_THROW(parseWithCdxGam({"-t", "abc"}), std::runtime_error);
    }

    TEST(CliThreadOptionTest, DecompressionThreadsZeroThrows) {
        EXPECT_THROW(parseWithCdxGam({"-T", "0"}), std::runtime_error);
    }

    TEST(CliThreadOptionTest, DecompressionThreadsNegativeThrows) {
        EXPECT_THROW(parseWithCdxGam({"-T", "-1"}), std::runtime_error);
    }

    TEST(CliThreadOptionTest, DecompressionThreadsNonNumericStringThrows) {
        EXPECT_THROW(parseWithCdxGam({"-T", "xyz"}), std::runtime_error);
    }
} // anonymous namespace

// =============================================================================
// -i/--inspect
// =============================================================================
namespace {
    TEST(CliInspectOptionTest, FlagAloneEnablesWithNoComponentSelected) {
        const CliArgs args = parseWithCdxOnly({"-i"});
        EXPECT_TRUE(args.inspectMode());
        EXPECT_FALSE(args.inspect.component.has_value());
    }

    TEST(CliInspectOptionTest, FlagWithComponentValue) {
        const CliArgs args = parseWithCdxOnly({"-i", "chr1"});
        EXPECT_TRUE(args.inspectMode());
        ASSERT_TRUE(args.inspect.component.has_value());
        EXPECT_EQ(*args.inspect.component, "chr1");
    }

    TEST(CliInspectOptionTest, AbsentByDefault) {
        EXPECT_FALSE(parseWithCdxGam().inspectMode());
    }
} // anonymous namespace

// =============================================================================
// -o/--output and --no-graph/--no-stats/--no-table
// =============================================================================
namespace {
    TEST(CliOutputOptionsTest, DefaultOutputDirectoryIsCurrentDir) {
        EXPECT_EQ(parseWithCdxGam().output_directory, ".");
    }

    TEST(CliOutputOptionsTest, CustomOutputDirectory) {
        EXPECT_EQ(parseWithCdxGam({"-o", "/tmp/some/dir"}).output_directory, "/tmp/some/dir");
    }

    TEST(CliOutputOptionsTest, NoGraphFlag) {
        const CliArgs args = parseWithCdxGam({"--no-graph"});
        EXPECT_TRUE(args.no_graph);
        EXPECT_FALSE(args.generateGraph());
    }

    TEST(CliOutputOptionsTest, NoStatsFlag) {
        const CliArgs args = parseWithCdxGam({"--no-stats"});
        EXPECT_TRUE(args.no_stats);
        EXPECT_FALSE(args.generateStats());
    }

    TEST(CliOutputOptionsTest, NoTableFlag) {
        const CliArgs args = parseWithCdxGam({"--no-table"});
        EXPECT_TRUE(args.no_table);
        EXPECT_FALSE(args.generateTable());
    }

    TEST(CliOutputOptionsTest, TwoDisabledOutputsStillValid) {
        const CliArgs args = parseWithCdxGam({"--no-graph", "--no-stats"});
        EXPECT_TRUE(args.generateTable());
    }

    TEST(CliOutputOptionsTest, AllThreeOutputsDisabledThrows) {
        EXPECT_THROW(parseWithCdxGam({"--no-graph", "--no-stats", "--no-table"}), std::runtime_error);
    }
} // anonymous namespace

// =============================================================================
// Graph rendering options: --log, --smoothing, --max-points, --dpi, --fig-size,
// --color-line, --color-filling.
// =============================================================================
namespace {
    TEST(CliGraphOptionsTest, LogFlagAbsentByDefault) {
        EXPECT_FALSE(parseWithCdxGam().log_base.has_value());
    }

    TEST(CliGraphOptionsTest, LogFlagAloneDefaultsToBase10) {
        const CliArgs args = parseWithCdxGam({"--log"});
        ASSERT_TRUE(args.log_base.has_value());
        EXPECT_EQ(*args.log_base, 10);
    }

    TEST(CliGraphOptionsTest, LogFlagWithExplicitBase) {
        const CliArgs args = parseWithCdxGam({"--log", "2"});
        ASSERT_TRUE(args.log_base.has_value());
        EXPECT_EQ(*args.log_base, 2);
    }

    TEST(CliGraphOptionsTest, LogFlagBelowRangeThrows) {
        EXPECT_THROW(parseWithCdxGam({"--log", "1"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, LogFlagAboveRangeThrows) {
        EXPECT_THROW(parseWithCdxGam({"--log", "100000"}), std::runtime_error);
    }

    // --log is bound directly to a CLI11-typed `int` - this is a genuine
    // CLI11-level type-conversion failure, not a manual std::sto* one.
    TEST(CliGraphOptionsTest, LogFlagNonIntegerThrows) {
        EXPECT_THROW(parseWithCdxGam({"--log", "abc"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, SmoothingDefault) {
        EXPECT_DOUBLE_EQ(parseWithCdxGam().smoothing, 0.01);
    }

    TEST(CliGraphOptionsTest, SmoothingCustomValue) {
        EXPECT_DOUBLE_EQ(parseWithCdxGam({"--smoothing", "0.5"}).smoothing, 0.5);
    }

    TEST(CliGraphOptionsTest, SmoothingBelowRangeThrows) {
        EXPECT_THROW(parseWithCdxGam({"--smoothing", "-0.1"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, SmoothingAboveRangeThrows) {
        EXPECT_THROW(parseWithCdxGam({"--smoothing", "1.5"}), std::runtime_error);
    }

    // --smoothing is bound directly to a CLI11-typed `double`.
    TEST(CliGraphOptionsTest, SmoothingNonNumericThrows) {
        EXPECT_THROW(parseWithCdxGam({"--smoothing", "abc"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, MaxPointsDefault) {
        EXPECT_EQ(parseWithCdxGam().max_plot_points, 10000u);
    }

    TEST(CliGraphOptionsTest, MaxPointsCustomValue) {
        EXPECT_EQ(parseWithCdxGam({"--max-points", "500"}).max_plot_points, 500u);
    }

    // Documented as the way to disable downsampling - must not throw.
    TEST(CliGraphOptionsTest, MaxPointsZeroIsAllowed) {
        EXPECT_EQ(parseWithCdxGam({"--max-points", "0"}).max_plot_points, 0u);
    }

    // --max-points is bound directly to a CLI11-typed `size_t`.
    TEST(CliGraphOptionsTest, MaxPointsNonIntegerThrows) {
        EXPECT_THROW(parseWithCdxGam({"--max-points", "abc"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, MaxPointsNegativeThrows) {
        EXPECT_THROW(parseWithCdxGam({"--max-points", "-1"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, DpiDefault) {
        EXPECT_EQ(parseWithCdxGam().dpi, 300);
    }

    TEST(CliGraphOptionsTest, DpiCustomValue) {
        EXPECT_EQ(parseWithCdxGam({"--dpi", "150"}).dpi, 150);
    }

    TEST(CliGraphOptionsTest, DpiZeroThrows) {
        EXPECT_THROW(parseWithCdxGam({"--dpi", "0"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, DpiNegativeThrows) {
        EXPECT_THROW(parseWithCdxGam({"--dpi", "-5"}), std::runtime_error);
    }

    // --dpi is bound directly to a CLI11-typed `int`.
    TEST(CliGraphOptionsTest, DpiNonIntegerThrows) {
        EXPECT_THROW(parseWithCdxGam({"--dpi", "abc"}), std::runtime_error);
    }

    // --- --fig-size ---

    TEST(CliGraphOptionsTest, FigSizeDefault_NoQuery_Linear) {
        const CliArgs args = parseWithCdxGam();
        EXPECT_DOUBLE_EQ(args.figure_width, 5.5);
        EXPECT_DOUBLE_EQ(args.figure_height, 3.5);
        EXPECT_FALSE(args.custom_figure_size.has_value());
    }

    TEST(CliGraphOptionsTest, FigSizeDefault_NoQuery_Circular) {
        const CliArgs args = parseWithCdxGam({"-c", "circular"});
        EXPECT_DOUBLE_EQ(args.figure_width, 4.0);
        EXPECT_DOUBLE_EQ(args.figure_height, 4.0);
    }

    TEST(CliGraphOptionsTest, FigSizeDefault_WithQuery_Linear) {
        const CliArgs args = parseWithCdxGam({"-q", "chr1"});
        EXPECT_DOUBLE_EQ(args.figure_width, 7.0);
        EXPECT_DOUBLE_EQ(args.figure_height, 4.5);
    }

    TEST(CliGraphOptionsTest, FigSizeDefault_WithQuery_Circular) {
        const CliArgs args = parseWithCdxGam({"-q", "chr1", "-c", "circular"});
        EXPECT_DOUBLE_EQ(args.figure_width, 7.0);
        EXPECT_DOUBLE_EQ(args.figure_height, 7.0);
    }

    TEST(CliGraphOptionsTest, FigSizeCustomOverridesDefault) {
        const CliArgs args = parseWithCdxGam({"--fig-size", "7x4.5"});
        EXPECT_DOUBLE_EQ(args.figure_width, 7.0);
        EXPECT_DOUBLE_EQ(args.figure_height, 4.5);
        ASSERT_TRUE(args.custom_figure_size.has_value());
        EXPECT_DOUBLE_EQ(args.custom_figure_size->first, 7.0);
        EXPECT_DOUBLE_EQ(args.custom_figure_size->second, 4.5);
    }

    TEST(CliGraphOptionsTest, FigSizeAcceptsUppercaseXSeparator) {
        const CliArgs args = parseWithCdxGam({"--fig-size", "10X8"});
        EXPECT_DOUBLE_EQ(args.figure_width, 10.0);
        EXPECT_DOUBLE_EQ(args.figure_height, 8.0);
    }

    TEST(CliGraphOptionsTest, FigSizeInvalidFormatThrows) {
        EXPECT_THROW(parseWithCdxGam({"--fig-size", "not-a-size"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, FigSizeZeroWidthThrows) {
        EXPECT_THROW(parseWithCdxGam({"--fig-size", "0x5"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, FigSizeNegativeHeightThrows) {
        // The regex itself only accepts unsigned digits, so a leading '-' simply
        // fails to match the whole pattern - still must throw.
        EXPECT_THROW(parseWithCdxGam({"--fig-size", "5x-3"}), std::runtime_error);
    }

    // --- colors ---

    TEST(CliGraphOptionsTest, ColorLineDefault) {
        EXPECT_EQ(parseWithCdxGam().line_color, "#1E3A8A");
    }

    TEST(CliGraphOptionsTest, ColorLineLowercaseIsNormalizedToUppercase) {
        EXPECT_EQ(parseWithCdxGam({"--color-line", "#1e3a8a"}).line_color, "#1E3A8A");
    }

    TEST(CliGraphOptionsTest, ColorLineShorthandFormAccepted) {
        EXPECT_EQ(parseWithCdxGam({"--color-line", "#fff"}).line_color, "#FFF");
    }

    TEST(CliGraphOptionsTest, ColorLineWithAlphaAccepted) {
        EXPECT_EQ(parseWithCdxGam({"--color-line", "#1e3a8aff"}).line_color, "#1E3A8AFF");
    }

    TEST(CliGraphOptionsTest, ColorLineMissingHashThrows) {
        EXPECT_THROW(parseWithCdxGam({"--color-line", "1E3A8A"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, ColorLineNonHexCharactersThrow) {
        EXPECT_THROW(parseWithCdxGam({"--color-line", "#GGGGGG"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, ColorLineWrongLengthThrows) {
        EXPECT_THROW(parseWithCdxGam({"--color-line", "#12345"}), std::runtime_error);
    }

    TEST(CliGraphOptionsTest, ColorFillDefault) {
        EXPECT_EQ(parseWithCdxGam().fill_color, "#93C5FD");
    }

    TEST(CliGraphOptionsTest, ColorFillCustomValue) {
        EXPECT_EQ(parseWithCdxGam({"--color-filling", "#00FF00"}).fill_color, "#00FF00");
    }

    TEST(CliGraphOptionsTest, ColorFillInvalidThrows) {
        EXPECT_THROW(parseWithCdxGam({"--color-filling", "not-a-color"}), std::runtime_error);
    }
} // anonymous namespace

// =============================================================================
// --help: the one remaining std::exit() path (exit code 0, not an error).
// Uses GoogleTest death tests since the process legitimately terminates.
// =============================================================================
namespace {
    void invokeHelp() {
        const CliArgs args = parseWithCdxGam({"--help"});
        (void) args; // unreachable: --help exits the process before returning
    }

    TEST(CliHelpTest, HelpFlagExitsWithCodeZero) {
        EXPECT_EXIT(invokeHelp(), ::testing::ExitedWithCode(0), "");
    }

    void invokeShortHelp() {
        const CliArgs args = parseWithCdxGam({"-h"});
        (void) args;
    }

    TEST(CliHelpTest, ShortHelpFlagExitsWithCodeZero) {
        EXPECT_EXIT(invokeShortHelp(), ::testing::ExitedWithCode(0), "");
    }
} // anonymous namespace
