/**
 * @file linear_plot.cpp
 * @brief Numeric preparation shared by the linear coverage plot renderers
 *        (downsampling, smoothing, log transform, grid layout, tick math).
 *
 * @note Le rendu proprement dit (Cairo, PNG) vit dans cairo_plot.cpp. Ce
 *       fichier ne produit que des structures de données prêtes à dessiner.
 */

#include "output_plot.h"
#include "config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

namespace output {
    // Computes the grid layout for a global graph components.
    GridLayout chooseGlobalGraphGrid(const std::size_t component_count) {
        // Validate that component count falls within the supported range [1, 30]
        if (component_count < 1 || component_count > 30) {
            throw std::invalid_argument("The global graph supports between 1 and 30 components.");
        }

        // Optimal grid configurations for small component counts (1 to 9)
        switch (component_count) {
            case 1: return {1, 1};
            case 2: return {1, 2};
            case 3: return {1, 3};
            case 4: return {2, 2};
            case 5: return {2, 3};
            case 6: return {2, 3};
            case 7: return {3, 3};
            case 8: return {2, 4};
            case 9: return {3, 3};
            default:
                break;
        }

        // Automatically calculate balanced row/column dimensions for counts from 10 to 25
        const auto columns = static_cast<std::size_t>(
            std::ceil(std::sqrt(static_cast<double>(component_count)))
        );

        const auto rows = static_cast<std::size_t>(
            std::ceil(static_cast<double>(component_count) / static_cast<double>(columns))
        );

        return {rows, columns};
    }

    namespace {
        /**
         * @brief Prepares and optimizes coverage data for plotting by performing downsampling and moving-average smoothing.
         *
         * This function serves as the shared backend for formatting coverage datasets. It downsamples
         * high-resolution datasets to fit within specified point budgets (`max_plot_points`) and computes
         * sliding-window moving averages using prefix sum arrays (cumulative sums) for $O(1)$ window queries.
         * It properly handles boundary conditions for both **Linear** and **Circular** topologies (such as wrapping around
         * circular chromosomes).
         *
         * @param coverage Vector of double-precision coverage values (supporting NaN/sentinel representations).
         * @param smoothing Moving-average window size expressed as a fraction of total points [0.0, 1.0].
         * @param max_plot_points Maximum number of points allowed in the downsampled output (0 = full resolution).
         * @param topology Graph topology mode (Topology::Linear or Topology::Circular).
         * @return PlotData Processed structure containing downsampled X coordinates, smoothed Y values, and the applied window size.
         *
         * @throws std::invalid_argument If @p smoothing is non-finite or outside the range [0.0, 1.0].
         */
        PlotData prepareCoverageForPlotImpl(
            const std::vector<double> &coverage,
            const double smoothing,
            const std::size_t max_plot_points,
            const Topology topology
        ) {
            // Validate smoothing parameters
            if (!std::isfinite(smoothing)) {
                throw std::invalid_argument("smoothing must be finite.");
            }
            if (smoothing > 1.0) {
                throw std::invalid_argument("smoothing must be between 0 and 1.");
            }

            const std::size_t n_points = coverage.size();

            // Return empty dataset immediately if input coverage is empty
            if (n_points == 0) {
                return {};
            }

            // max_plot_points == 0 means "no downsampling": preserve full-resolution
            // dataset plotting (effectively mapping max_plot_points = n_points).
            const std::size_t effective_max_points = (max_plot_points == 0) ? n_points : max_plot_points;

            // Calculate step stride to sample points uniformly across the dataset
            const std::size_t step = std::max<std::size_t>(
                1, (n_points + effective_max_points - 1) / effective_max_points
            );

            PlotData result;
            result.x.reserve(std::min(effective_max_points, n_points));

            // Populate downsampled X-axis coordinate indices
            for (std::size_t i = 0; i < n_points; i += step) {
                result.x.push_back(i);
            }

            // Ensure the final data point (last index) is always included for accurate trailing bounds
            if (result.x.back() != n_points - 1) {
                if (result.x.size() < effective_max_points) {
                    result.x.push_back(n_points - 1);
                } else {
                    result.x.back() = n_points - 1;
                }
            }

            // Helper lambda to bypass smoothing and populate raw data points directly
            auto fill_raw_values = [&]() -> PlotData {
                result.y.reserve(result.x.size());
                for (const std::size_t idx: result.x) {
                    result.y.push_back(coverage[idx]);
                }
                result.window_size = 1;
                return result;
            };

            // If smoothing is disabled (<= 0.0), return raw values straightaway
            if (smoothing <= 0.0) return fill_raw_values();

            // Compute smoothing window size based on the dataset fraction
            std::size_t window_size =
                    std::max<std::size_t>(1, std::llround(static_cast<double>(n_points) * smoothing));
            window_size = std::min(window_size, n_points);

            // Force window size to be odd for symmetric centering, unless it matches total points
            if (window_size > 1 && window_size % 2 == 0) {
                if (window_size < n_points) ++window_size;
                else --window_size;
            }
            if (window_size <= 1) return fill_raw_values();

            // Build a prefix sum (cumulative sum) array for O(1) range query performance
            std::vector<double> cumulative_sum(n_points + 1, 0.0);
            for (std::size_t i = 0; i < n_points; ++i) {
                cumulative_sum[i + 1] = cumulative_sum[i] + coverage[i];
            }

            const std::size_t half_window = window_size / 2;
            result.y.reserve(result.x.size());

            // Apply moving average calculation depending on structural topology
            if (topology == Topology::Circular) {
                // Circular topology: handles wrapping past the boundaries of the genome/component
                for (const std::size_t center: result.x) {
                    const std::size_t start = (center + n_points - half_window) % n_points;
                    const std::size_t end = start + window_size;

                    double sum = 0.0;
                    if (end <= n_points) {
                        sum = cumulative_sum[end] - cumulative_sum[start];
                    } else {
                        // Handle wrap-around past the end back to the beginning of the sequence
                        const std::size_t wrapped_end = end - n_points;
                        sum = cumulative_sum[n_points] - cumulative_sum[start] + cumulative_sum[wrapped_end];
                    }

                    result.y.push_back(sum / static_cast<double>(window_size));
                }
            } else {
                // Linear topology: clips windows to absolute structural boundaries
                for (const std::size_t center: result.x) {
                    const std::size_t start = center > half_window ? center - half_window : 0;
                    const std::size_t end = std::min(center + half_window + 1, n_points);

                    const double sum = cumulative_sum[end] - cumulative_sum[start];
                    const std::size_t length = end - start;

                    // Normalize by actual window length (important near edges where window is truncated)
                    result.y.push_back(sum / static_cast<double>(length));
                }
            }

            result.window_size = window_size;
            return result;
        }
    } // namespace

    // Adapter overload of prepareCoverageForPlotImpl converts raw integer coverage data to double
    // used for plotting linear graph
    PlotData prepareCoverageForPlot(
        const std::vector<cdx::Coverage> &coverage,
        const double smoothing,
        const std::size_t max_plot_points,
        const Topology topology
    ) {
        std::vector<double> widened(coverage.size());
        for (std::size_t i = 0; i < coverage.size(); ++i) {
            widened[i] = static_cast<double>(coverage[i]);
        }
        return prepareCoverageForPlotImpl(widened, smoothing, max_plot_points, topology);
    }


    // Adapter overload of prepareCoverageForPlotImpl without doing anything
    // used for plotting circular graph
    PlotData prepareCoverageForPlot(
        const std::vector<double> &coverage,
        const double smoothing,
        const std::size_t max_plot_points,
        const Topology topology
    ) {
        return prepareCoverageForPlotImpl(coverage, smoothing, max_plot_points, topology);
    }

    // Compute Y-axis label for logarithmic log coverage
    LogCoverageTicks logCoverageTicks(
        const double maximum_coverage,
        const int log_base
    ) {
        if (log_base <= 1) {
            throw std::invalid_argument("log_base must be greater than 1.");
        }
        // Return baseline trivial ticks if maximum coverage is zero or negative
        if (maximum_coverage <= 0.0) {
            return {{0.0}, {0.0}, 1.0};
        }

        // Determine the maximum exponent needed to span up to the maximum coverage value
        const auto max_exponent = std::max(
            0,
            static_cast<int>(std::ceil(std::log(maximum_coverage) / std::log(static_cast<double>(log_base))))
        );

        // Build raw tick values starting from 0.0, 1.0, and scaling up by powers of log_base
        std::vector<double> raw_ticks;
        raw_ticks.reserve(static_cast<std::size_t>(max_exponent) + 2);
        raw_ticks.push_back(0.0);
        raw_ticks.push_back(1.0);
        double value = log_base;

        for (int exponent = 1; exponent <= max_exponent; ++exponent) {
            raw_ticks.push_back(value);
            value *= static_cast<double>(log_base);
        }

        // Transform raw ticks into their corresponding log-scale plotting positions
        std::vector<double> display_ticks = logTransformCoverage(raw_ticks, log_base);

        // Important: Capture the upper limit value BEFORE moving display_ticks into
        // the return list. Moving display_ticks invalidates its contents; calling
        // .back() afterwards would read from a moved-from vector (undefined behavior -> SIGSEGV).
        const double display_upper_limit = display_ticks.back();

        return {
            std::move(raw_ticks),
            std::move(display_ticks),
            display_upper_limit
        };
    }

    // chose nice label for the Y-axis
    CoverageTicks calculateCoverageTicks(
        const double maximum_coverage,
        const std::size_t target_tick_count
    ) {
        // Validate inputs for finiteness and boundary constraints
        if (!std::isfinite(maximum_coverage)) {
            throw std::invalid_argument("maximum_coverage must be finite.");
        }
        if (target_tick_count < 2) {
            throw std::invalid_argument("target_tick_count must be at least 2.");
        }

        // Return baseline trivial ticks if maximum coverage is zero or negative
        if (maximum_coverage <= 0.0) {
            return {{1.0}, 1.0};
        }

        // Determine the raw step size and extract its order of magnitude (base-10 exponent)
        const double raw_step = maximum_coverage / static_cast<double>(target_tick_count);
        const double magnitude = std::pow(10.0, std::floor(std::log10(raw_step)));
        const double normalized_step = raw_step / magnitude;

        // Standard multipliers used for selecting "nice" human-readable axis steps
        constexpr std::array<double, 7> multipliers{
            1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0
        };

        double multiplier = 10.0;

        // Find the smallest nice multiplier greater than or equal to the normalized step
        for (const double candidate: multipliers) {
            if (candidate >= normalized_step) {
                multiplier = candidate;
                break;
            }
        }

        // Compute the final clean tick step and rounded upper limit boundary
        const double tick_step = multiplier * magnitude;
        const double upper_limit = std::ceil(maximum_coverage / tick_step) * tick_step;
        std::vector<double> tick_values;
        const auto tick_count = static_cast<std::size_t>(std::llround(upper_limit / tick_step));

        // Generate individual tick positions, including a tolerance margin to avoid floating-point drift
        tick_values.reserve(tick_count);
        for (double value2 = tick_step;
             value2 <= upper_limit + tick_step * 0.5;
             value2 += tick_step) {
            tick_values.push_back(value2);
        }

        return {std::move(tick_values), upper_limit};
    }


    // Transforms a vector of coverage values onto a custom logarithmic scale
    std::vector<double> logTransformCoverage(
        const std::vector<double> &values,
        const int log_base
    ) {
        // Validate that the log base is greater than 1
        if (log_base <= 1) {
            throw std::invalid_argument("log_base must be greater than 1.");
        }

        // Precompute the natural logarithm of the base for change-of-base calculations
        const double log_base_factor = std::log(static_cast<double>(log_base));

        std::vector<double> transformed_val;
        transformed_val.reserve(values.size());

        for (const double value: values) {
            // Ensure coverage values are non-negative
            if (value < 0.0) {
                throw std::invalid_argument("Coverage values cannot be negative.");
            }

            // Keep values <= 1.0 linear to ensure smooth handling near zero/one bounds
            if (value <= 1.0) {
                transformed_val.push_back(value);
            } else {
                // Apply custom log base transformation formula: 1.0 + log_base(value)
                transformed_val.push_back(1.0 + std::log(value) / log_base_factor);
            }
        }

        return transformed_val;
    }


    // Prepares and formats a linear coverage plot package for rendering and visualization
    cdx::LinearPlotPackageBin prepareLinearPlotPackage(
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name,
        std::size_t offset,
        const PlotConfig &config
    ) {
        cdx::LinearPlotPackageBin pkg;

        // Set basic structural metadata and query boundaries
        pkg.component_name = component_name;
        pkg.query_start = offset;
        pkg.query_end = coverage.empty() ? offset : (offset + coverage.size() - 1);

        // =====================================================================
        // 1. Downsampling & Smoothing Preparation
        // =====================================================================
        // Filter out sentinel values representing positions outside the query scope, replacing them with 0
        std::vector<cdx::Coverage> filtered_coverage;
        filtered_coverage.reserve(coverage.size());

        for (const auto value: coverage) {
            if (value >= cfg::NOT_IN_QUERY) {
                filtered_coverage.push_back(0);
            } else {
                filtered_coverage.push_back(value);
            }
        }

        // Execute core downsampling and smoothing routines for linear topology
        const PlotData plot_data = prepareCoverageForPlot(
            filtered_coverage,
            config.smoothing,
            config.max_plot_points,
            Topology::Linear
        );

        // Configure linear X-axis coordinate mapping steps
        // (max_plot_points == 0 => no downsampling, render at full resolution)
        const std::size_t n_points = coverage.size();
        const std::size_t effective_max_points = (config.max_plot_points == 0) ? n_points : config.max_plot_points;
        const std::size_t step = (n_points == 0 || effective_max_points == 0)
                                     ? 1
                                     : std::max<std::size_t>(
                                         1, (n_points + effective_max_points - 1) / effective_max_points);

        pkg.x_start = static_cast<double>(offset);
        pkg.x_step = static_cast<double>(step);

        pkg.y.reserve(plot_data.y.size());
        double maximum_raw = 0.0;

        // =====================================================================
        // 2. Map Processed Y Values and Track Peak Raw Coverage
        // =====================================================================
        for (std::size_t i = 0; i < plot_data.x.size(); ++i) {
            const double value = plot_data.y[i];
            pkg.y.push_back(value);
            maximum_raw = std::max(maximum_raw, value);
        }

        // =====================================================================
        // 3. Vertical Scale Configuration (Logarithmic vs. Linear)
        // =====================================================================
        if (config.log_base.has_value()) {
            const int log_base = *config.log_base;

            pkg.logarithmic = true;
            pkg.log_base = log_base;

            // Apply logarithmic transformation to the Y values
            pkg.y = logTransformCoverage(pkg.y, log_base);

            if (maximum_raw > 0.0) {
                // Compute tailored tick positions and upper limits for log scale
                const LogCoverageTicks ticks = logCoverageTicks(maximum_raw, log_base);
                pkg.tick_positions = ticks.display_values;
                pkg.y_upper_limit = std::max(ticks.display_upper_limit, 1.0);

                // Generate formatted string labels for each log tick (e.g., "10x", "100x")
                pkg.tick_labels.reserve(ticks.raw_values.size());
                for (const double raw: ticks.raw_values) {
                    std::ostringstream s;
                    s << std::fixed << std::setprecision(0) << raw << 'x';
                    pkg.tick_labels.push_back(s.str());
                }
            } else {
                // Fallback default ticks if coverage is entirely zero
                pkg.tick_positions = {0.0, 1.0};
                pkg.tick_labels = {"0x", "1x"};
                pkg.y_upper_limit = 1.0;
            }
        } else {
            // Standard linear scale configuration with a 10% headroom buffer
            pkg.logarithmic = false;
            pkg.log_base = 0;
            pkg.y_upper_limit = (maximum_raw > 0.0) ? (maximum_raw * 1.10) : 1.0;
        }

        return pkg;
    }
} // namespace output
