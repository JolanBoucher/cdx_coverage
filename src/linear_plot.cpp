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
    GridLayout chooseGlobalGraphGrid(const std::size_t component_count) {
        if (component_count < 1 || component_count > 25) {
            throw std::invalid_argument("The global graph supports between 1 and 25 components.");
        }

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
         * @brief Shared implementation for both prepareCoverageForPlot
         *        overloads. `coverage` is already a double vector at this
         *        point (either widened from cdx::Coverage, or the caller's
         *        own NaN-capable array, e.g. the circular backend's masked
         *        sentinels - cf. output_plot.h).
         */
        PlotData prepareCoverageForPlotImpl(
            const std::vector<double> &coverage,
            const double smoothing,
            const std::size_t max_plot_points,
            const Topology topology
        ) {
            if (!std::isfinite(smoothing)) {
                throw std::invalid_argument("smoothing must be finite.");
            }
            if (smoothing > 1.0) {
                throw std::invalid_argument("smoothing must be between 0 and 1.");
            }

            const std::size_t n_points = coverage.size();

            if (n_points == 0) {
                return {};
            }

            // max_plot_points == 0 signifie "pas de réduction de points" : on trace
            // la résolution complète (équivalent à max_plot_points = n_points).
            const std::size_t effective_max_points = (max_plot_points == 0) ? n_points : max_plot_points;

            const std::size_t step = std::max<std::size_t>(
                1, (n_points + effective_max_points - 1) / effective_max_points
            );

            PlotData result;
            result.x.reserve(std::min(effective_max_points, n_points));

            for (std::size_t i = 0; i < n_points; i += step) {
                result.x.push_back(i);
            }

            if (result.x.back() != n_points - 1) {
                if (result.x.size() < effective_max_points) {
                    result.x.push_back(n_points - 1);
                } else {
                    result.x.back() = n_points - 1;
                }
            }

            auto fill_raw_values = [&]() -> PlotData {
                result.y.reserve(result.x.size());
                for (const std::size_t idx: result.x) {
                    result.y.push_back(coverage[idx]);
                }
                result.window_size = 1;
                return result;
            };

            if (smoothing <= 0.0) return fill_raw_values();

            std::size_t window_size =
                    std::max<std::size_t>(1, std::llround(static_cast<double>(n_points) * smoothing));
            window_size = std::min(window_size, n_points);

            if (window_size > 1 && window_size % 2 == 0) {
                if (window_size < n_points) ++window_size;
                else --window_size;
            }
            if (window_size <= 1) return fill_raw_values();

            std::vector<double> cumulative_sum(n_points + 1, 0.0);
            for (std::size_t i = 0; i < n_points; ++i) {
                cumulative_sum[i + 1] = cumulative_sum[i] + coverage[i];
            }

            const std::size_t half_window = window_size / 2;
            result.y.reserve(result.x.size());

            if (topology == Topology::Circular) {
                for (const std::size_t center: result.x) {
                    const std::size_t start = (center + n_points - half_window) % n_points;
                    const std::size_t end = start + window_size;

                    double sum = 0.0;
                    if (end <= n_points) {
                        sum = cumulative_sum[end] - cumulative_sum[start];
                    } else {
                        const std::size_t wrapped_end = end - n_points;
                        sum = cumulative_sum[n_points] - cumulative_sum[start] + cumulative_sum[wrapped_end];
                    }

                    result.y.push_back(sum / static_cast<double>(window_size));
                }
            } else {
                for (const std::size_t center: result.x) {
                    const std::size_t start = center > half_window ? center - half_window : 0;
                    const std::size_t end = std::min(center + half_window + 1, n_points);

                    const double sum = cumulative_sum[end] - cumulative_sum[start];
                    const std::size_t length = end - start;

                    result.y.push_back(sum / static_cast<double>(length));
                }
            }

            result.window_size = window_size;
            return result;
        }
    } // namespace

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

    PlotData prepareCoverageForPlot(
        const std::vector<double> &coverage,
        const double smoothing,
        const std::size_t max_plot_points,
        const Topology topology
    ) {
        return prepareCoverageForPlotImpl(coverage, smoothing, max_plot_points, topology);
    }

    LogCoverageTicks logCoverageTicks(
        const double maximum_coverage,
        const int log_base
    ) {
        if (log_base <= 1) {
            throw std::invalid_argument("log_base must be greater than 1.");
        }

        if (maximum_coverage <= 0.0) {
            return {{0.0}, {0.0}, 1.0};
        }

        const auto max_exponent = std::max(
            0,
            static_cast<int>(std::ceil(std::log(maximum_coverage) / std::log(static_cast<double>(log_base))))
        );

        std::vector<double> raw_ticks;
        raw_ticks.reserve(static_cast<std::size_t>(max_exponent) + 2);
        raw_ticks.push_back(0.0);
        raw_ticks.push_back(1.0);
        double value = log_base;

        for (int exponent = 1; exponent <= max_exponent; ++exponent) {
            raw_ticks.push_back(value);
            value *= static_cast<double>(log_base);
        }
        std::vector<double> display_ticks = logTransformCoverage(raw_ticks, log_base);

        // Important : capturer la valeur AVANT de déplacer display_ticks dans
        // la liste d'initialisation. std::move(display_ticks) ci-dessous vide
        // le vecteur ; appeler .back() après coup lirait un vecteur déplacé
        // (comportement indéfini -> SIGSEGV observé avec --log).
        const double display_upper_limit = display_ticks.back();

        return {
            std::move(raw_ticks),
            std::move(display_ticks),
            display_upper_limit
        };
    }

    CoverageTicks calculateCoverageTicks(
        const double maximum_coverage,
        const std::size_t target_tick_count
    ) {
        if (!std::isfinite(maximum_coverage)) {
            throw std::invalid_argument("maximum_coverage must be finite.");
        }

        if (target_tick_count < 2) {
            throw std::invalid_argument("target_tick_count must be at least 2.");
        }

        if (maximum_coverage <= 0.0) {
            return {{1.0}, 1.0};
        }

        const double raw_step = maximum_coverage / static_cast<double>(target_tick_count);
        const double magnitude = std::pow(10.0, std::floor(std::log10(raw_step)));
        const double normalized_step = raw_step / magnitude;
        constexpr std::array<double, 7> multipliers{
            1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0
        };

        double multiplier = 10.0;

        for (const double candidate: multipliers) {
            if (candidate >= normalized_step) {
                multiplier = candidate;
                break;
            }
        }

        const double tick_step = multiplier * magnitude;
        const double upper_limit = std::ceil(maximum_coverage / tick_step) * tick_step;
        std::vector<double> tick_values;
        const auto tick_count = static_cast<std::size_t>(std::llround(upper_limit / tick_step));

        tick_values.reserve(tick_count);
        for (double value2 = tick_step;
             value2 <= upper_limit + tick_step * 0.5;
             value2 += tick_step) {
            tick_values.push_back(value2);
        }

        return {std::move(tick_values), upper_limit};
    }

    std::vector<double> logTransformCoverage(
        const std::vector<double> &values,
        const int log_base
    ) {
        if (log_base <= 1) {
            throw std::invalid_argument("log_base must be greater than 1.");
        }

        const double log_base_factor = std::log(static_cast<double>(log_base));

        std::vector<double> transformed_val;
        transformed_val.reserve(values.size());

        for (const double value: values) {
            if (value < 0.0) {
                throw std::invalid_argument("Coverage values cannot be negative.");
            }
            if (value <= 1.0) {
                transformed_val.push_back(value);
            } else {
                transformed_val.push_back(1.0 + std::log(value) / log_base_factor);
            }
        }

        return transformed_val;
    }

    cdx::LinearPlotPackageBin prepareLinearPlotPackage(
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name,
        std::size_t offset,
        const PlotConfig &config
    ) {
        cdx::LinearPlotPackageBin pkg;

        pkg.component_name = component_name;
        pkg.query_start = offset;
        pkg.query_end = coverage.empty() ? offset : (offset + coverage.size() - 1);

        // ---------------------------------------------------------------------
        // 1. Downsampling & Lissage
        // ---------------------------------------------------------------------
        std::vector<cdx::Coverage> filtered_coverage;
        filtered_coverage.reserve(coverage.size());

        for (const auto value: coverage) {
            if (value >= cfg::NOT_IN_QUERY) {
                filtered_coverage.push_back(0);
            } else {
                filtered_coverage.push_back(value);
            }
        }

        const PlotData plot_data = prepareCoverageForPlot(
            filtered_coverage,
            config.smoothing,
            config.max_plot_points,
            Topology::Linear
        );

        // Configuration du repère linéaire X
        // (max_plot_points == 0 => pas de réduction, résolution complète)
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

        // ---------------------------------------------------------------------
        // 2. Remplacement des sentinelles Y par 0.0
        // ---------------------------------------------------------------------
        for (std::size_t i = 0; i < plot_data.x.size(); ++i) {
            const double value = plot_data.y[i];
            pkg.y.push_back(value);
            maximum_raw = std::max(maximum_raw, value);
        }

        // ---------------------------------------------------------------------
        // 3. Échelle verticale (Log vs Linéaire)
        // ---------------------------------------------------------------------
        if (config.log_base.has_value()) {
            const int log_base = *config.log_base;

            pkg.logarithmic = true;
            pkg.log_base = log_base;

            pkg.y = logTransformCoverage(pkg.y, log_base);

            if (maximum_raw > 0.0) {
                const LogCoverageTicks ticks = logCoverageTicks(maximum_raw, log_base);
                pkg.tick_positions = ticks.display_values;
                pkg.y_upper_limit = std::max(ticks.display_upper_limit, 1.0);

                pkg.tick_labels.reserve(ticks.raw_values.size());
                for (const double raw: ticks.raw_values) {
                    std::ostringstream s;
                    s << std::fixed << std::setprecision(0) << raw << 'x';
                    pkg.tick_labels.push_back(s.str());
                }
            } else {
                pkg.tick_positions = {0.0, 1.0};
                pkg.tick_labels = {"0x", "1x"};
                pkg.y_upper_limit = 1.0;
            }
        } else {
            pkg.logarithmic = false;
            pkg.log_base = 0;
            pkg.y_upper_limit = (maximum_raw > 0.0) ? (maximum_raw * 1.10) : 1.0;
        }

        return pkg;
    }
} // namespace output
