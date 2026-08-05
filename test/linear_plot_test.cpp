/**
 * @file linear_plot_test.cpp
 * @brief Unit tests for linear_plot.cpp.
 *
 * Unlike cairo_plot.cpp (actual Cairo/PNG rendering) and python_circular_plot.cpp (subprocess
 * invocation), everything in linear_plot.cpp is pure numerical computation with no I/O, so it is
 * fully unit tested here rather than treated as a lighter integration module. Only
 * `prepareCoverageForPlotImpl` lives in an anonymous namespace; it is exercised through its two
 * public overloads, `prepareCoverageForPlot(vector<cdx::Coverage>, ...)` and
 * `prepareCoverageForPlot(vector<double>, ...)`, which both delegate to it directly.
 */

#include "../src/output_plot.h"
#include "../src/config.h"
#include "cdx_types.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

// =============================================================================
// chooseGlobalGraphGrid
//
// Covers the 9 hand-tuned layouts, the sqrt-based formula for 10-30, and the
// [1, 30] validity range.
// =============================================================================
namespace {
    TEST(ChooseGlobalGraphGridTest, ZeroThrows) {
        EXPECT_THROW(output::chooseGlobalGraphGrid(0), std::invalid_argument);
    }

    TEST(ChooseGlobalGraphGridTest, AboveThirtyThrows) {
        EXPECT_THROW(output::chooseGlobalGraphGrid(31), std::invalid_argument);
    }

    TEST(ChooseGlobalGraphGridTest, HandTunedLayoutsOneThroughNine) {
        const std::vector<output::GridLayout> expected = {
            {1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {2, 3}, {3, 3}, {2, 4}, {3, 3}
        };
        for (std::size_t n = 1; n <= 9; ++n) {
            const output::GridLayout got = output::chooseGlobalGraphGrid(n);
            EXPECT_EQ(got.rows, expected[n - 1].rows) << "n=" << n;
            EXPECT_EQ(got.columns, expected[n - 1].columns) << "n=" << n;
        }
    }

    TEST(ChooseGlobalGraphGridTest, AutoFormulaTen) {
        const output::GridLayout grid = output::chooseGlobalGraphGrid(10);
        EXPECT_EQ(grid.columns, 4u); // ceil(sqrt(10)) = 4
        EXPECT_EQ(grid.rows, 3u); // ceil(10/4) = 3
    }

    TEST(ChooseGlobalGraphGridTest, AutoFormulaPerfectSquareSixteen) {
        const output::GridLayout grid = output::chooseGlobalGraphGrid(16);
        EXPECT_EQ(grid.columns, 4u);
        EXPECT_EQ(grid.rows, 4u);
    }

    TEST(ChooseGlobalGraphGridTest, AutoFormulaSeventeen) {
        const output::GridLayout grid = output::chooseGlobalGraphGrid(17);
        EXPECT_EQ(grid.columns, 5u); // ceil(sqrt(17)) = 5
        EXPECT_EQ(grid.rows, 4u); // ceil(17/5) = 4
    }

    TEST(ChooseGlobalGraphGridTest, AutoFormulaPerfectSquareTwentyFive) {
        const output::GridLayout grid = output::chooseGlobalGraphGrid(25);
        EXPECT_EQ(grid.columns, 5u);
        EXPECT_EQ(grid.rows, 5u);
    }

    // Upper boundary of the supported range (recently raised from 25 to 30).
    TEST(ChooseGlobalGraphGridTest, AutoFormulaUpperBoundaryThirty) {
        const output::GridLayout grid = output::chooseGlobalGraphGrid(30);
        EXPECT_EQ(grid.columns, 6u); // ceil(sqrt(30)) = 6
        EXPECT_EQ(grid.rows, 5u); // ceil(30/6) = 5
    }
} // anonymous namespace

// =============================================================================
// prepareCoverageForPlot (downsampling + moving-average smoothing)
// =============================================================================
namespace {
    TEST(PrepareCoverageForPlotTest, EmptyInputReturnsDefaultPlotData) {
        const output::PlotData data = output::prepareCoverageForPlot(
            std::vector<double>{}, 0.1, 100, output::Topology::Linear
        );
        EXPECT_TRUE(data.x.empty());
        EXPECT_TRUE(data.y.empty());
        EXPECT_EQ(data.window_size, 1u); // default member value, untouched
    }

    TEST(PrepareCoverageForPlotTest, NonFiniteSmoothingThrows) {
        const std::vector<double> coverage{1, 2, 3};
        EXPECT_THROW(
            output::prepareCoverageForPlot(coverage, std::numeric_limits<double>::quiet_NaN(), 0,
                output::Topology::Linear),
            std::invalid_argument
        );
    }

    TEST(PrepareCoverageForPlotTest, SmoothingAboveOneThrows) {
        const std::vector<double> coverage{1, 2, 3};
        EXPECT_THROW(
            output::prepareCoverageForPlot(coverage, 1.5, 0, output::Topology::Linear),
            std::invalid_argument
        );
    }

    // Negative smoothing is not explicitly rejected by the validation checks
    // (only non-finite or >1.0 is) - it simply behaves like smoothing=0 (raw
    // pass-through), since the only other check is `smoothing <= 0.0`.
    TEST(PrepareCoverageForPlotTest, NegativeSmoothingBehavesLikeDisabled) {
        const std::vector<double> coverage{1, 2, 3, 4, 5};
        const output::PlotData data = output::prepareCoverageForPlot(coverage, -0.5, 0, output::Topology::Linear);
        EXPECT_EQ(data.window_size, 1u);
        EXPECT_EQ(data.y, coverage);
    }

    // max_plot_points == 0 means full resolution: every point is kept.
    TEST(PrepareCoverageForPlotTest, MaxPlotPointsZeroMeansNoDownsampling) {
        const std::vector<double> coverage{10, 20, 30, 40, 50};
        const output::PlotData data = output::prepareCoverageForPlot(coverage, 0.0, 0, output::Topology::Linear);
        EXPECT_EQ(data.x, (std::vector<std::size_t>{0, 1, 2, 3, 4}));
        EXPECT_EQ(data.y, coverage);
        EXPECT_EQ(data.window_size, 1u);
    }

    // Downsampling forces the final index to always be included, whether it
    // needs to be appended (list not yet at capacity) or must overwrite the
    // last sampled index (list already at capacity).
    TEST(PrepareCoverageForPlotTest, DownsamplingAlwaysIncludesLastIndex) {
        // n=10, max_plot_points=3 -> step=ceil(10/3)=4 -> raw samples at 0,4,8;
        // list already has 3 elements (== effective_max_points), so index 8 is
        // overwritten by the forced last index 9, not appended.
        std::vector<double> coverage(10);
        for (std::size_t i = 0; i < 10; ++i) coverage[i] = static_cast<double>(i);

        const output::PlotData data = output::prepareCoverageForPlot(coverage, 0.0, 3, output::Topology::Linear);
        EXPECT_EQ(data.x, (std::vector<std::size_t>{0, 4, 9}));
    }

    TEST(PrepareCoverageForPlotTest, DownsamplingAppendsLastIndexWhenBelowCapacity) {
        // n=9, max_plot_points=5 -> step=ceil(9/5)=2 -> raw samples at 0,2,4,6,8;
        // last sample (8) already equals n-1, so no adjustment needed at all
        // here - use a case where it doesn't land exactly instead.
        // n=7, max_plot_points=5 -> step=ceil(7/5)=2 -> samples 0,2,4,6; size=4 < 5
        // (capacity), so the missing last index (6) ... already included. Use
        // n=8, max_plot_points=5 -> step=ceil(8/5)=2 -> samples 0,2,4,6; size=4<5,
        // last index is 7, not yet present -> appended.
        std::vector<double> coverage(8);
        for (std::size_t i = 0; i < 8; ++i) coverage[i] = static_cast<double>(i);

        const output::PlotData data = output::prepareCoverageForPlot(coverage, 0.0, 5, output::Topology::Linear);
        EXPECT_EQ(data.x, (std::vector<std::size_t>{0, 2, 4, 6, 7}));
    }

    // window_size rounds to 0 (tiny smoothing fraction on a large dataset) -> raw pass-through.
    TEST(PrepareCoverageForPlotTest, TinySmoothingCollapsesToRawValues) {
        std::vector<double> coverage(100);
        for (std::size_t i = 0; i < 100; ++i) coverage[i] = static_cast<double>(i);

        const output::PlotData data = output::prepareCoverageForPlot(coverage, 0.001, 0, output::Topology::Linear);
        EXPECT_EQ(data.window_size, 1u);
        EXPECT_EQ(data.y, coverage);
    }

    // Linear topology: window clipped at the boundaries (edge points average
    // over fewer than window_size elements, normalized by actual length).
    TEST(PrepareCoverageForPlotTest, LinearTopologySmoothingWithClippedEdges) {
        const std::vector<double> coverage{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        // smoothing=0.3 -> window_size = round(10*0.3) = 3 (already odd).
        const output::PlotData data = output::prepareCoverageForPlot(coverage, 0.3, 0, output::Topology::Linear);

        ASSERT_EQ(data.window_size, 3u);
        const std::vector<double> expected{1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.5};
        ASSERT_EQ(data.y.size(), expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) {
            EXPECT_DOUBLE_EQ(data.y[i], expected[i]) << "i=" << i;
        }
    }

    // Circular topology: the same window wraps around the array boundary
    // instead of being clipped, producing different results at the edges.
    TEST(PrepareCoverageForPlotTest, CircularTopologySmoothingWrapsAround) {
        const std::vector<double> coverage{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        const output::PlotData data = output::prepareCoverageForPlot(coverage, 0.3, 0, output::Topology::Circular);

        ASSERT_EQ(data.window_size, 3u);
        const std::vector<double> expected{
            13.0 / 3.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 20.0 / 3.0
        };
        ASSERT_EQ(data.y.size(), expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) {
            EXPECT_NEAR(data.y[i], expected[i], 1e-9) << "i=" << i;
        }
    }

    // An even computed window size below n_points is bumped up by one to stay odd.
    TEST(PrepareCoverageForPlotTest, EvenWindowBelowNPointsIsIncrementedToOdd) {
        std::vector<double> coverage(8, 0.0);
        // smoothing=0.5 -> round(8*0.5) = 4 (even, < 8) -> becomes 5.
        const output::PlotData data = output::prepareCoverageForPlot(coverage, 0.5, 0, output::Topology::Linear);
        EXPECT_EQ(data.window_size, 5u);
    }

    // An even computed window size equal to n_points is decremented instead
    // (incrementing would exceed the dataset size).
    TEST(PrepareCoverageForPlotTest, EvenWindowEqualToNPointsIsDecrementedToOdd) {
        std::vector<double> coverage(4, 0.0);
        // smoothing=1.0 -> round(4*1.0) = 4 (even, == 4) -> becomes 3.
        const output::PlotData data = output::prepareCoverageForPlot(coverage, 1.0, 0, output::Topology::Linear);
        EXPECT_EQ(data.window_size, 3u);
    }

    // The cdx::Coverage (uint32) overload must simply widen to double and
    // delegate, producing identical results to the double overload.
    TEST(PrepareCoverageForPlotTest, IntegerOverloadMatchesDoubleOverload) {
        const std::vector<cdx::Coverage> int_coverage{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        const std::vector<double> double_coverage{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        const output::PlotData from_int = output::prepareCoverageForPlot(int_coverage, 0.3, 0, output::Topology::Linear);
        const output::PlotData from_double = output::prepareCoverageForPlot(
            double_coverage, 0.3, 0, output::Topology::Linear
        );

        EXPECT_EQ(from_int.x, from_double.x);
        EXPECT_EQ(from_int.window_size, from_double.window_size);
        ASSERT_EQ(from_int.y.size(), from_double.y.size());
        for (std::size_t i = 0; i < from_int.y.size(); ++i) {
            EXPECT_DOUBLE_EQ(from_int.y[i], from_double.y[i]);
        }
    }
} // anonymous namespace

// =============================================================================
// calculateCoverageTicks (linear axis "nice" tick positions)
// =============================================================================
namespace {
    TEST(CalculateCoverageTicksTest, NonFiniteMaximumThrows) {
        EXPECT_THROW(
            output::calculateCoverageTicks(std::numeric_limits<double>::infinity()),
            std::invalid_argument
        );
    }

    TEST(CalculateCoverageTicksTest, TargetTickCountBelowTwoThrows) {
        EXPECT_THROW(output::calculateCoverageTicks(100.0, 1), std::invalid_argument);
    }

    TEST(CalculateCoverageTicksTest, ZeroOrNegativeMaximumReturnsTrivialBaseline) {
        const output::CoverageTicks ticks = output::calculateCoverageTicks(0.0);
        EXPECT_EQ(ticks.values, (std::vector<double>{1.0}));
        EXPECT_DOUBLE_EQ(ticks.upper_limit, 1.0);

        const output::CoverageTicks negative_ticks = output::calculateCoverageTicks(-5.0);
        EXPECT_EQ(negative_ticks.values, (std::vector<double>{1.0}));
    }

    TEST(CalculateCoverageTicksTest, RoundsUpToNiceStepOfTen) {
        // raw_step = 37/4 = 9.25 -> magnitude 1 -> normalized 9.25 -> nearest
        // multiplier >= 9.25 is 10 -> tick_step = 10, upper_limit = ceil(3.7)*10 = 40.
        const output::CoverageTicks ticks = output::calculateCoverageTicks(37.0, 4);
        EXPECT_EQ(ticks.values, (std::vector<double>{10.0, 20.0, 30.0, 40.0}));
        EXPECT_DOUBLE_EQ(ticks.upper_limit, 40.0);
    }

    TEST(CalculateCoverageTicksTest, RoundsUpToNiceStepOfTwoPointFive) {
        // raw_step = 1234/5 = 246.8 -> magnitude 100 -> normalized 2.468 ->
        // nearest multiplier >= 2.468 is 2.5 -> tick_step = 250,
        // upper_limit = ceil(4.936)*250 = 1250.
        const output::CoverageTicks ticks = output::calculateCoverageTicks(1234.0, 5);
        EXPECT_EQ(ticks.values, (std::vector<double>{250.0, 500.0, 750.0, 1000.0, 1250.0}));
        EXPECT_DOUBLE_EQ(ticks.upper_limit, 1250.0);
    }
} // anonymous namespace

// =============================================================================
// logTransformCoverage
// =============================================================================
namespace {
    TEST(LogTransformCoverageTest, LogBaseOfOneOrLessThrows) {
        EXPECT_THROW(output::logTransformCoverage({1.0}, 1), std::invalid_argument);
        EXPECT_THROW(output::logTransformCoverage({1.0}, 0), std::invalid_argument);
        EXPECT_THROW(output::logTransformCoverage({1.0}, -3), std::invalid_argument);
    }

    TEST(LogTransformCoverageTest, NegativeValueThrows) {
        EXPECT_THROW(output::logTransformCoverage({-1.0}, 10), std::invalid_argument);
    }

    TEST(LogTransformCoverageTest, ValuesAtOrBelowOneStayLinear) {
        const std::vector<double> result = output::logTransformCoverage({0.0, 0.5, 1.0}, 10);
        ASSERT_EQ(result.size(), 3u);
        EXPECT_DOUBLE_EQ(result[0], 0.0);
        EXPECT_DOUBLE_EQ(result[1], 0.5);
        EXPECT_DOUBLE_EQ(result[2], 1.0);
    }

    TEST(LogTransformCoverageTest, ValuesAboveOneAreLogTransformed) {
        // 1 + log_10(100) = 1 + 2 = 3 ; 1 + log_2(8) = 1 + 3 = 4
        const std::vector<double> base10 = output::logTransformCoverage({100.0}, 10);
        EXPECT_NEAR(base10[0], 3.0, 1e-9);

        const std::vector<double> base2 = output::logTransformCoverage({8.0}, 2);
        EXPECT_NEAR(base2[0], 4.0, 1e-9);
    }

    TEST(LogTransformCoverageTest, MixedVector) {
        const std::vector<double> result = output::logTransformCoverage({0.0, 0.5, 1.0, 2.0, 100.0}, 10);
        ASSERT_EQ(result.size(), 5u);
        EXPECT_DOUBLE_EQ(result[0], 0.0);
        EXPECT_DOUBLE_EQ(result[1], 0.5);
        EXPECT_DOUBLE_EQ(result[2], 1.0);
        EXPECT_NEAR(result[3], 1.0 + std::log10(2.0), 1e-9);
        EXPECT_NEAR(result[4], 3.0, 1e-9);
    }
} // anonymous namespace

// =============================================================================
// logCoverageTicks
// =============================================================================
namespace {
    TEST(LogCoverageTicksTest, LogBaseOfOneOrLessThrows) {
        EXPECT_THROW(output::logCoverageTicks(100.0, 1), std::invalid_argument);
    }

    TEST(LogCoverageTicksTest, ZeroOrNegativeMaximumReturnsTrivialBaseline) {
        const output::LogCoverageTicks ticks = output::logCoverageTicks(0.0, 10);
        EXPECT_EQ(ticks.raw_values, (std::vector<double>{0.0}));
        EXPECT_EQ(ticks.display_values, (std::vector<double>{0.0}));
        EXPECT_DOUBLE_EQ(ticks.display_upper_limit, 1.0);
    }

    // maximum_coverage == 1 exactly: log(1) == 0, so max_exponent collapses to
    // 0 and the power-of-base loop never runs.
    TEST(LogCoverageTicksTest, MaximumOfOneProducesNoExponentTicks) {
        const output::LogCoverageTicks ticks = output::logCoverageTicks(1.0, 10);
        EXPECT_EQ(ticks.raw_values, (std::vector<double>{0.0, 1.0}));
        EXPECT_EQ(ticks.display_values, (std::vector<double>{0.0, 1.0}));
        EXPECT_DOUBLE_EQ(ticks.display_upper_limit, 1.0);
    }

    TEST(LogCoverageTicksTest, Base10UpToFifty) {
        // ceil(log_10(50)) = ceil(1.69897) = 2 -> raw ticks 0, 1, 10, 100.
        const output::LogCoverageTicks ticks = output::logCoverageTicks(50.0, 10);
        EXPECT_EQ(ticks.raw_values, (std::vector<double>{0.0, 1.0, 10.0, 100.0}));

        ASSERT_EQ(ticks.display_values.size(), 4u);
        EXPECT_DOUBLE_EQ(ticks.display_values[0], 0.0);
        EXPECT_DOUBLE_EQ(ticks.display_values[1], 1.0);
        EXPECT_NEAR(ticks.display_values[2], 2.0, 1e-9); // 1 + log10(10)
        EXPECT_NEAR(ticks.display_values[3], 3.0, 1e-9); // 1 + log10(100)
        EXPECT_NEAR(ticks.display_upper_limit, 3.0, 1e-9);
    }

    TEST(LogCoverageTicksTest, Base2UpToTen) {
        // ceil(log_2(10)) = ceil(3.3219) = 4 -> raw ticks 0, 1, 2, 4, 8, 16.
        const output::LogCoverageTicks ticks = output::logCoverageTicks(10.0, 2);
        EXPECT_EQ(ticks.raw_values, (std::vector<double>{0.0, 1.0, 2.0, 4.0, 8.0, 16.0}));
        EXPECT_NEAR(ticks.display_upper_limit, 5.0, 1e-9); // 1 + log2(16)
    }
} // anonymous namespace

// =============================================================================
// prepareLinearPlotPackage
// =============================================================================
namespace {
    // Sentinel filtering (NOT_IN_QUERY -> 0) and standard linear-scale headroom.
    TEST(PrepareLinearPlotPackageTest, SentinelFilteringAndLinearScaleHeadroom) {
        const std::vector<cdx::Coverage> coverage{0, 5, 10, cfg::NOT_IN_QUERY, 3};
        output::PlotConfig config; // defaults: smoothing=0.01, max_plot_points=10000, no log_base

        const cdx::LinearPlotPackageBin pkg = output::prepareLinearPlotPackage(coverage, "chrT", 100, config);

        EXPECT_EQ(pkg.component_name, "chrT");
        EXPECT_EQ(pkg.query_start, 100u);
        EXPECT_EQ(pkg.query_end, 104u); // 100 + 5 - 1
        EXPECT_DOUBLE_EQ(pkg.x_start, 100.0);
        EXPECT_DOUBLE_EQ(pkg.x_step, 1.0);
        EXPECT_FALSE(pkg.logarithmic);
        EXPECT_EQ(pkg.log_base, 0);

        // window_size collapses to 1 here (round(5*0.01)=0 -> clamped to 1 ->
        // "<=1" -> raw pass-through), so sentinel filtering is directly visible.
        EXPECT_EQ(pkg.y, (std::vector<double>{0.0, 5.0, 10.0, 0.0, 3.0}));

        // 10% headroom over the peak raw value (10).
        EXPECT_DOUBLE_EQ(pkg.y_upper_limit, 11.0);
    }

    TEST(PrepareLinearPlotPackageTest, EmptyCoverageFallsBackToUnitUpperLimit) {
        output::PlotConfig config;
        const cdx::LinearPlotPackageBin pkg = output::prepareLinearPlotPackage({}, "chrEmpty", 50, config);

        EXPECT_EQ(pkg.query_start, 50u);
        EXPECT_EQ(pkg.query_end, 50u); // empty coverage: query_end falls back to offset
        EXPECT_TRUE(pkg.y.empty());
        EXPECT_FALSE(pkg.logarithmic);
        EXPECT_DOUBLE_EQ(pkg.y_upper_limit, 1.0);
    }

    // Logarithmic scale: y is log-transformed, ticks/labels come from logCoverageTicks.
    TEST(PrepareLinearPlotPackageTest, LogarithmicScaleTransformsYAndBuildsTicks) {
        const std::vector<cdx::Coverage> coverage{0, 5, 10, cfg::NOT_IN_QUERY, 3};
        output::PlotConfig config;
        config.log_base = 2;

        const cdx::LinearPlotPackageBin pkg = output::prepareLinearPlotPackage(coverage, "chrT", 0, config);

        EXPECT_TRUE(pkg.logarithmic);
        EXPECT_EQ(pkg.log_base, 2);

        // Raw filtered values [0, 5, 10, 0, 3] log2-transformed:
        // 0 -> 0 ; 5 -> 1+log2(5) ; 10 -> 1+log2(10) ; 0 -> 0 ; 3 -> 1+log2(3)
        ASSERT_EQ(pkg.y.size(), 5u);
        EXPECT_NEAR(pkg.y[0], 0.0, 1e-9);
        EXPECT_NEAR(pkg.y[1], 1.0 + std::log2(5.0), 1e-9);
        EXPECT_NEAR(pkg.y[2], 1.0 + std::log2(10.0), 1e-9);
        EXPECT_NEAR(pkg.y[3], 0.0, 1e-9);
        EXPECT_NEAR(pkg.y[4], 1.0 + std::log2(3.0), 1e-9);

        // Ticks are derived from the peak *raw* value (10), not the
        // already-log-transformed y: raw ticks 0,1,2,4,8,16 (see
        // LogCoverageTicksTest.Base2UpToTen), formatted as "Nx".
        ASSERT_EQ(pkg.tick_labels.size(), 6u);
        EXPECT_EQ(pkg.tick_labels, (std::vector<std::string>{"0x", "1x", "2x", "4x", "8x", "16x"}));
        EXPECT_NEAR(pkg.y_upper_limit, 5.0, 1e-9); // 1 + log2(16)
    }

    // All-zero coverage in log mode: maximum_raw stays 0, hitting the
    // fallback {0x, 1x} ticks instead of logCoverageTicks.
    TEST(PrepareLinearPlotPackageTest, LogarithmicScaleWithAllZeroCoverageFallsBackToUnitTicks) {
        const std::vector<cdx::Coverage> coverage{0, 0, 0};
        output::PlotConfig config;
        config.log_base = 5;

        const cdx::LinearPlotPackageBin pkg = output::prepareLinearPlotPackage(coverage, "chrZero", 0, config);

        EXPECT_TRUE(pkg.logarithmic);
        EXPECT_EQ(pkg.tick_positions, (std::vector<double>{0.0, 1.0}));
        EXPECT_EQ(pkg.tick_labels, (std::vector<std::string>{"0x", "1x"}));
        EXPECT_DOUBLE_EQ(pkg.y_upper_limit, 1.0);
    }
} // anonymous namespace
