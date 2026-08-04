/**
 * @file output_stats.cpp
 * @brief Coverage statistics calculation and text report generation.
 */

#include "output_stats.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
    constexpr cdx::Coverage MAX_DENSE_HISTOGRAM = 1'000'000;

    /* @brief Intermediate accumulator aggregating frequencies and statistical moments. */
    struct CoverageAccumulator {
        std::size_t region_length = 0;              // Total number of positions contributing to this accumulator.
        std::size_t covered_positions = 0;          // Number of positions with non-zero coverage.
        std::uint64_t sum = 0;                      // Sum of all coverage values across the region.
        double sum_squares = 0.0;                   // Sum of squared coverage values
        cdx::Coverage min =                         // Lowest observed coverage value.
            std::numeric_limits<cdx::Coverage>::max();
        cdx::Coverage max = 0;                      // Highest observed coverage value.
        bool requires_fallback = false;             // Indicates that exact statistics could not be computed directly
        std::vector<std::size_t> histogram;         // Coverage frequency table

        /**
         * @brief Merge coverage statistics from another accumulator.
         *
         * Adds the contents of @p other to this accumulator, combining aggregate
         * coverage metrics, histogram counts, and minimum/maximum coverage values.
         * If @p other represents an empty region (region_length == 0), no changes
         * are made.
         *
         * The histogram is automatically expanded as needed to accommodate the
         * coverage range present in @p other.
         *
         * @param other Accumulator whose statistics will be incorporated into this
         *        accumulator.
         *
         * @throws std::overflow_error In debug builds if the merged coverage sum
         *         would exceed the capacity of std::uint64_t.
         */
        void merge(const CoverageAccumulator &other) {
            if (other.region_length == 0) return;

#ifndef NDEBUG
            if (other.sum > std::numeric_limits<std::uint64_t>::max() - sum) {
                throw std::overflow_error("Coverage sum exceeds uint64_t capacity.");
            }
#endif

            // filling the member with the struct
            region_length += other.region_length;
            covered_positions += other.covered_positions;
            sum += other.sum;
            sum_squares += other.sum_squares;

            min = std::min(min, other.min);
            max = std::max(max, other.max);

            if (histogram.size() < other.histogram.size()) {
                histogram.resize(other.histogram.size(), 0);
            }

            for (std::size_t i = 0; i < other.histogram.size(); ++i) {
                histogram[i] += other.histogram[i];
            }
        }
    };

    /**
     * @brief Format an unsigned integer with comma thousands separators.
     *
     * Converts the supplied integer to its decimal string representation and
     * inserts commas every three digits, counting from the rightmost digit.
     *
     * Examples:
     *   42       -> "42"
     *   1234     -> "1,234"
     *   1234567  -> "1,234,567"
     *
     * @param value Integer value to format.
     * @return The formatted string representation of @p value.
     */
    [[nodiscard]]
    std::string formatInteger(const std::size_t value) {
        std::string result = std::to_string(value);
        for (std::ptrdiff_t position = static_cast<std::ptrdiff_t>(result.size()) - 3; position > 0; position -= 3) {
            result.insert(static_cast<std::size_t>(position), 1, ',');
        }
        return result;
    }

    /**
     * @brief Write a formatted coverage-statistics report section.
     *
     * Emits a human-readable summary of coverage metrics for a genomic region,
     * component, or dataset. The section includes region size, coverage breadth,
     * central-tendency statistics, dispersion metrics, quartiles, and coverage
     * extrema, formatted as aligned key/value pairs.
     *
     * @param stream Output stream receiving the formatted report section.
     * @param title Section title displayed in the report header.
     * @param stats Coverage statistics to display.
     */
    void writeCoverageStatsSection(
        std::ostream &stream,
        const std::string &title,
        const output::CoverageStats &stats
    ) {
        // Write section header and separator line.
        stream << "[ " << title << " ]\n" << std::string(70, '-') << '\n';
        stream << std::left << std::setw(30) << "Indexed positions" << " : " << formatInteger(stats.region_length) << '\n';
        stream << std::left << std::setw(30) << "Covered positions" << " : " << formatInteger(stats.covered_positions) << '\n';

        // Use fixed-point formatting for percentage and coverage statistics.
        stream << std::fixed << std::setprecision(2);
        stream << std::left << std::setw(30) << "Breadth of coverage" << " : " << stats.breadth << "%\n";
        stream << std::left << std::setw(30) << "Mean coverage" << " : " << stats.mean << "x\n";
        stream << std::left << std::setw(30) << "Median coverage" << " : " << stats.median << "x\n";
        stream << std::left << std::setw(30) << "Standard deviation" << " : " << stats.stddev << '\n';

        // Display coefficient of variation with additional precision since
        // values are often small and sensitive to rounding.
        stream << std::setprecision(4) << std::left << std::setw(30) << "Coefficient of variation" << " : " << stats.cv << '\n';
        stream << std::setprecision(2) << std::left << std::setw(30) << "Q1" << " : " << stats.q1 << "x\n";
        stream << std::left << std::setw(30) << "Q3" << " : " << stats.q3 << "x\n";

        // Switch back to default formatting for integer-like coverage extrema.
        stream << std::defaultfloat << std::left << std::setw(30) << "Minimum coverage" << " : " << stats.min << "x\n";
        stream << std::left << std::setw(30) << "Maximum coverage" << " : " << stats.max << "x\n";
    }

    /**
     * @brief Compute coverage statistics using selection-based quantile estimation.
     *
     * This fallback implementation is used when the coverage range is too large
     * for efficient dense-histogram processing (i.e. maximum coverage exceeds
     * MAX_DENSE_HISTOGRAM). Valid coverage values are copied into a temporary
     * container, after which summary statistics and quantiles are computed
     * directly from the observed values.
     *
     * Coverage entries greater than or equal to cfg::NOT_IN_QUERY are ignored.
     *
     * Quartiles and the median are computed using repeated calls to
     * std::nth_element, avoiding the memory cost of a dense coverage histogram
     * while maintaining linear-average-time selection.
     *
     * @tparam Iterator Iterator type over a sequence of coverage values.
     * @param begin Iterator to the first coverage value.
     * @param end Iterator one past the last coverage value.
     * @return Computed coverage statistics for all valid positions in the range.
     */
    template<typename Iterator>
    output::CoverageStats computeCoverageStatsNthElement(Iterator begin, Iterator end) {
        output::CoverageStats stats;
        std::vector<cdx::Coverage> valid_cov;

        // Collect coverage values corresponding to indexed positions.
        for (auto iterator = begin; iterator != end; ++iterator) {
            const cdx::Coverage value = *iterator;
            if (value < cfg::NOT_IN_QUERY) {
                valid_cov.push_back(value);
            }
        }

        stats.region_length = valid_cov.size();
        if (valid_cov.empty()) return stats;

        stats.min = valid_cov.front();
        stats.max = valid_cov.front();

        std::uint64_t sum = 0;
        double sum_squares = 0.0;

        // Accumulate coverage moments and extrema used by summary statistics.
        for (const auto value: valid_cov) {
            if (value != 0) ++stats.covered_positions;
            if (value < stats.min) stats.min = value;
            if (value > stats.max) stats.max = value;

            sum += value;
            const auto d_val = static_cast<double>(value);
            sum_squares += d_val * d_val;
        }

        const auto N = static_cast<double>(stats.region_length);
        const double mean = static_cast<double>(sum) / N;
        const double variance = sum_squares / N - mean * mean;

        stats.mean = mean;
        stats.stddev = std::sqrt(std::max(0.0, variance));
        stats.cv = stats.mean > 0.0 ? stats.stddev / stats.mean : 0.0;
        stats.breadth = 100.0 * static_cast<double>(stats.covered_positions) / N;

        // Compute an interpolated quantile using nth_element instead of sorting
        // the entire coverage vector.
        auto get_quantile = [&](const double percentile) -> double {
            const double pos = percentile * (N - 1.0);
            const auto lower_idx = static_cast<std::size_t>(std::floor(pos));
            const auto upper_idx = static_cast<std::size_t>(std::ceil(pos));

            std::nth_element(valid_cov.begin(), valid_cov.begin() + static_cast<std::ptrdiff_t>(lower_idx),
                             valid_cov.end());
            const double lower_val = valid_cov[lower_idx];

            if (lower_idx == upper_idx) return lower_val;

            std::nth_element(valid_cov.begin() + static_cast<std::ptrdiff_t>(lower_idx) + 1,
                             valid_cov.begin() + static_cast<std::ptrdiff_t>(upper_idx),
                             valid_cov.end());
            const double upper_val = valid_cov[upper_idx];

            return lower_val + (pos - static_cast<double>(lower_idx)) * (upper_val - lower_val);
        };

        // get put the quantile in the structure
        stats.q1 = get_quantile(0.25);
        stats.median = get_quantile(0.50);
        stats.q3 = get_quantile(0.75);

        return stats;
    }

    /**
     * @brief Convert an accumulated coverage histogram into final coverage statistics.
     *
     * Computes summary statistics, including breadth of coverage, mean coverage,
     * standard deviation, coefficient of variation, median, and quartiles from
     * a CoverageAccumulator. Quantiles are derived directly from the coverage
     * histogram without materializing or sorting per-position coverage values.
     *
     * This implementation is intended for accumulators whose coverage range is
     * small enough to be represented efficiently by a dense histogram.
     *
     * @param accumulator Coverage accumulator containing aggregated moments,
     *        extrema, and coverage frequencies.
     * @return Fully populated CoverageStats structure. If the accumulator
     *         contains no positions, all statistics remain at their default
     *         values except min and max, which are set to zero.
     */
    output::CoverageStats finalizeStatsFromAccumulator(const CoverageAccumulator &accumulator) {
        output::CoverageStats stats;
        stats.region_length = accumulator.region_length;

        if (accumulator.region_length == 0) {
            stats.min = 0;
            stats.max = 0;
            return stats;
        }

        stats.min = accumulator.min;
        stats.max = accumulator.max;
        stats.covered_positions = accumulator.covered_positions;

        // Compute population mean and variance from accumulated moments.
        const auto N = static_cast<double>(accumulator.region_length);
        const double mean = static_cast<double>(accumulator.sum) / N;
        const double variance = accumulator.sum_squares / N - mean * mean;

        stats.mean = mean;
        stats.stddev = std::sqrt(std::max(0.0, variance));
        stats.cv = stats.mean > 0.0 ? stats.stddev / stats.mean : 0.0;
        stats.breadth = 100.0 * static_cast<double>(accumulator.covered_positions) / N;

        const double rank_q1 = 0.25 * (N - 1.0);
        const double rank_med = 0.50 * (N - 1.0);
        const double rank_q3 = 0.75 * (N - 1.0);

        // Determine the lower and upper order-statistic indices surrounding each
        // quantile rank. When a quantile falls between two positions, these bounds
        // are later used for linear interpolation.
        const auto r_l_q1  = static_cast<std::size_t>(std::floor(rank_q1));
        const auto r_u_q1  = static_cast<std::size_t>(std::ceil(rank_q1));
        const auto r_l_med = static_cast<std::size_t>(std::floor(rank_med));
        const auto r_u_med = static_cast<std::size_t>(std::ceil(rank_med));
        const auto r_l_q3  = static_cast<std::size_t>(std::floor(rank_q3));
        const auto r_u_q3  = static_cast<std::size_t>(std::ceil(rank_q3));

        std::size_t cumulative = 0;
        cdx::Coverage val_l_q1 = 0, val_u_q1 = 0;
        cdx::Coverage val_l_med = 0, val_u_med = 0;
        cdx::Coverage val_l_q3 = 0, val_u_q3 = 0;

        bool f_l_q1 = false, f_u_q1 = false;
        bool f_l_med = false, f_u_med = false;
        bool f_l_q3 = false, f_u_q3 = false;

        for (std::size_t value = accumulator.min; value <= accumulator.max; ++value) {
            if (value >= accumulator.histogram.size()) break;

            cumulative += accumulator.histogram[value];
            if (!f_l_q1 && cumulative > r_l_q1) {
                val_l_q1 = static_cast<cdx::Coverage>(value);
                f_l_q1 = true;
            }
            if (!f_u_q1 && cumulative > r_u_q1) {
                val_u_q1 = static_cast<cdx::Coverage>(value);
                f_u_q1 = true;
            }
            if (!f_l_med && cumulative > r_l_med) {
                val_l_med = static_cast<cdx::Coverage>(value);
                f_l_med = true;
            }
            if (!f_u_med && cumulative > r_u_med) {
                val_u_med = static_cast<cdx::Coverage>(value);
                f_u_med = true;
            }
            if (!f_l_q3 && cumulative > r_l_q3) {
                val_l_q3 = static_cast<cdx::Coverage>(value);
                f_l_q3 = true;
            }
            if (cumulative > r_u_q3) {
                val_u_q3 = static_cast<cdx::Coverage>(value);
                break;
            }
        }

        // Linearly interpolate a quantile value between the lower and upper order statistics
        // that bound the requested rank.
        auto interp = [](
            const double rank,
            const std::size_t r_l,
            const cdx::Coverage v_l,
            const cdx::Coverage v_u
        ) {
            const double fraction = rank - static_cast<double>(r_l);
            return static_cast<double>(v_l) + fraction * static_cast<double>(v_u - v_l);
        };

        stats.q1 = interp(rank_q1, r_l_q1, val_l_q1, val_u_q1);
        stats.median = interp(rank_med, r_l_med, val_l_med, val_u_med);
        stats.q3 = interp(rank_q3, r_l_q3, val_l_q3, val_u_q3);

        return stats;
    }

    /**
     * @brief Accumulate coverage statistics over a range of coverage values.
     *
     * Performs two passes over the input range. The first pass determines the
     * minimum and maximum valid coverage values and checks whether a dense
     * histogram can be used efficiently. The second pass populates the histogram
     * and accumulates the statistical moments required for downstream summary
     * statistics.
     *
     * Coverage values greater than or equal to cfg::NOT_IN_QUERY are treated as
     * sentinel values and are excluded from all calculations.
     *
     * If the observed maximum coverage exceeds MAX_DENSE_HISTOGRAM, histogram
     * construction is skipped and the returned accumulator is marked for fallback
     * processing via CoverageAccumulator::requires_fallback.
     *
     * @tparam Iterator Iterator type over a sequence of coverage values.
     * @param begin Iterator to the first coverage value.
     * @param end Iterator one past the last coverage value.
     * @return Accumulated coverage statistics and histogram information.
     *
     * @throws std::overflow_error In debug builds if the accumulated coverage
     *         sum exceeds the capacity of std::uint64_t.
     */
    template<typename Iterator>
    [[nodiscard]]
    CoverageAccumulator accumulateRange(
        Iterator begin,
        Iterator end
    ) {
        CoverageAccumulator acc;

        // Pass 1: identify the valid coverage range. This determines both the
        // histogram size and whether the dense-histogram approach is feasible.
        for (auto iterator = begin; iterator != end; ++iterator) {
            const cdx::Coverage value = *iterator;
            if (value >= cfg::NOT_IN_QUERY) {
                continue;
            }

            acc.min = std::min(acc.min, value);
            acc.max = std::max(acc.max, value);
        }

        // Empty range or range containing only sentinels.
        if (acc.min > acc.max) {
            acc.min = 0;
            acc.max = 0;
            return acc;
        }

        // Avoid allocating an excessively large dense histogram.
        if (acc.max > MAX_DENSE_HISTOGRAM) {
            acc.requires_fallback = true;
            return acc;
        }

        acc.histogram.assign(static_cast<std::size_t>(acc.max) + 1, 0);

        // Pass 2: populate the histogram and accumulate the moments needed to compute mean and variance.
        for (auto iterator = begin; iterator != end; ++iterator) {
            const cdx::Coverage value = *iterator;
            if (value >= cfg::NOT_IN_QUERY) {
                continue;
            }

            ++acc.region_length;
            if (value != 0) {
                ++acc.covered_positions;
            }

            // histogram[i] stores the number of positions observed with coverage depth i.
            ++acc.histogram[static_cast<std::size_t>(value)];

            const auto unsigned_value = static_cast<std::uint64_t>(value);
#ifndef NDEBUG
            if (unsigned_value > std::numeric_limits<std::uint64_t>::max() - acc.sum) {
                throw std::overflow_error("Coverage sum exceeds uint64_t capacity.");
            }
#endif
            acc.sum += unsigned_value;

            // Accumulate the second moment for later variance calculation.
            const double numeric_value = value;
            acc.sum_squares += numeric_value * numeric_value;
        }

        return acc;
    }
} // anonymous namespace

namespace output {

    // Compute coverage statistics for a coverage vector
    CoverageStats computeCoverageStats(
        const std::vector<cdx::Coverage> &coverage
    ) {
        const CoverageAccumulator accumulator = accumulateRange(
            coverage.begin(),
            coverage.end()
        );

        // Large coverage ranges are handled by a selection-based implementation
        // to avoid constructing an excessively large dense histogram.
        if (accumulator.requires_fallback) {
            return computeCoverageStatsNthElement(
                coverage.begin(),
                coverage.end()
            );
        }

        return finalizeStatsFromAccumulator(accumulator);
    }

    // Generate a coverage report for a single query component.
    void writeStatsReportQuery(
        const std::filesystem::path &output_path,
        const GamMappingStats &mapping_stats,
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name
    ) {

        // Validate basic accounting relationships between mapping categories.
        if (mapping_stats.mapped > mapping_stats.total) {
            throw std::invalid_argument("Mapped reads cannot exceed total reads.");
        }
        if (mapping_stats.unmapped > mapping_stats.total) {
            throw std::invalid_argument("Unmapped reads cannot exceed total reads.");
        }
        if (mapping_stats.mapped_to_query > mapping_stats.mapped) {
            throw std::invalid_argument("Reads mapped to query cannot exceed mapped reads.");
        }
        if (mapping_stats.mapped + mapping_stats.unmapped != mapping_stats.total) {
            throw std::invalid_argument("Mapped and unmapped reads must sum to total reads.");
        }

        const CoverageStats cov_stats = computeCoverageStats(coverage);

        const double pct_mapped = mapping_stats.total
                ? 100.0 * static_cast<double>(mapping_stats.mapped) / static_cast<double>(mapping_stats.total)
                : 0.0;

        const double pct_query_total = mapping_stats.total
                ? 100.0 * static_cast<double>(mapping_stats.mapped_to_query) / static_cast<double>(mapping_stats.total)
                : 0.0;

        const double pct_query_mapped = mapping_stats.mapped
                ? 100.0 * static_cast<double>(mapping_stats.mapped_to_query) / static_cast<double>(mapping_stats.mapped)
                : 0.0;

        const double pct_unmapped = mapping_stats.total
                ? 100.0 * static_cast<double>(mapping_stats.unmapped) / static_cast<double>(
                mapping_stats.total)
                : 0.0;

        std::ofstream report(output_path);
        if (!report) {
            throw std::runtime_error("Unable to open report file: " + output_path.string());
        }

        // Write report header and component metadata.
        report << std::string(70, '=') << '\n';
        report << "COVERAGE REPORT\n";
        report << "Component: " << component_name << '\n';
        report << std::string(70, '=') << "\n\n";

        report << "[ Mapping Statistics ]\n";
        report << std::string(70, '-') << '\n';
        report << std::left << std::setw(30) << "Total reads" << " : " << formatInteger(mapping_stats.total) << '\n';
        report << std::left << std::setw(30) << "Mapped reads" << " : " << formatInteger(mapping_stats.mapped)
                << " (" << std::fixed << std::setprecision(2) << pct_mapped << "%)\n";
        report << std::left << std::setw(30) << "Reads mapped to query" << " : " << formatInteger(
            mapping_stats.mapped_to_query) << '\n';
        report << std::left << std::setw(30) << "  Percentage of total" << " : " << pct_query_total << "%\n";
        report << std::left << std::setw(30) << "  Percentage of mapped" << " : " << pct_query_mapped << "%\n";
        report << std::left << std::setw(30) << "Unmapped reads" << " : " << formatInteger(mapping_stats.unmapped)
                << " (" << pct_unmapped << "%)\n\n";

        writeCoverageStatsSection(report, "Coverage Statistics", cov_stats);
    }

    // Generate a coverage report for a multi-component coverage table
    void writeStatsReportGlobal(
        const std::filesystem::path &output_path,
        const std::vector<cdx::Coverage> &coverage,
        const std::vector<cdx::PosBp> &component_offsets,
        const std::vector<std::string> &component_names,
        const GamMappingStats &mapping_stats,
        int threads
    ) {
        // Validate that component_offsets describes a complete partition of the flattened coverage table.
        if (component_offsets.size() < 2) {
            throw std::invalid_argument("bp_component_offsets must contain at least two boundaries.");
        }
        if (component_names.size() + 1 != component_offsets.size()) {
            throw std::invalid_argument("Component name count does not match component count.");
        }

        if (component_offsets.front() != 0) {
            throw std::invalid_argument("The first bp component offset must be zero.");
        }
        if (component_offsets.back() != static_cast<cdx::PosBp>(coverage.size())) {
            throw std::invalid_argument("The final bp component offset must match the flattened coverage table size.");
        }
        for (std::size_t i = 1; i < component_offsets.size(); ++i) {
            if (component_offsets[i] < component_offsets[i - 1]) {
                throw std::invalid_argument("bp_component_offsets must be non-decreasing.");
            }
        }

        std::ofstream report(output_path);
        if (!report) {
            throw std::runtime_error("Unable to open report file: " + output_path.string());
        }

        // Write report header and mapping statistics shared by all components.
        report << std::string(70, '=') << '\n';
        report << "COVERAGE REPORT\n";
        report << std::string(70, '=') << "\n\n";

        report << "[ Mapping Statistics ]\n";
        report << std::string(70, '-') << '\n';
        report << std::left << std::setw(30) << "Total reads" << " : " << formatInteger(mapping_stats.total) << '\n';
        report << std::left << std::setw(30) << "Mapped reads" << " : " << formatInteger(mapping_stats.mapped) << '\n';
        report << std::left << std::setw(30) << "Reads mapped to query" << " : " << formatInteger(
            mapping_stats.mapped_to_query) << "\n\n";

        const std::size_t component_count = component_offsets.size() - 1;

        // Store intermediate accumulators and final statistics for each component.
        std::vector<CoverageAccumulator> component_accumulators(component_count);
        std::vector<CoverageStats> component_stats(component_count);
        std::vector<std::uint8_t> component_uses_fallback(component_count, 0);

        {
// Compute component statistics independently; components can be processed
// in parallel because each iteration writes only to its own index.
#if defined(_OPENMP)

#pragma omp parallel for schedule(dynamic) num_threads(threads)

#endif
            for (std::ptrdiff_t component_index = 0;
                 component_index < static_cast<std::ptrdiff_t>(component_count);
                 ++component_index) {
                const auto component_id = static_cast<std::size_t>(component_index);
                const auto start_index = static_cast<std::size_t>(component_offsets[component_id]);
                const auto end_index = static_cast<std::size_t>(component_offsets[component_id + 1]);

                // Convert component boundary offsets into iterators spanning the
                // corresponding subrange of the flattened coverage table.
                const auto begin_iterator
                        = coverage.begin() + static_cast<std::ptrdiff_t>(start_index);
                const auto end_iterator
                        = coverage.begin() + static_cast<std::ptrdiff_t>(end_index);

                component_accumulators[component_id] = accumulateRange(
                    begin_iterator,
                    end_iterator
                );

                // Large coverage ranges are processed using the selection-based fallback
                // implementation instead of a dense histogram.
                if (component_accumulators[component_id].requires_fallback) {
                    component_uses_fallback[component_id] = 1;
                    component_stats[component_id] = computeCoverageStatsNthElement(
                        begin_iterator,
                        end_iterator
                    );
                } else {
                    component_stats[component_id] = finalizeStatsFromAccumulator(
                        component_accumulators[component_id]
                    );
                }
            }
        }

        CoverageAccumulator global_accumulator;
        bool global_fallback = false;

        // If all components used histogram-based accumulation, merge their accumulators
        // to construct global statistics without rescanning the coverage table.
        for (std::size_t component_id = 0; component_id < component_count; ++component_id) {
            if (component_uses_fallback[component_id] != 0) {
                global_fallback = true;
                break;
            }
            global_accumulator.merge(component_accumulators[component_id]);
        }

        CoverageStats global_stats;
        // If any component required fallback processing,
        // recompute global statistics directly from the full coverage table.
        if (global_fallback) {
            global_stats = computeCoverageStatsNthElement(
                coverage.begin(),
                coverage.end()
            );
        } else {
            global_stats = finalizeStatsFromAccumulator(global_accumulator);
        }

        // Write overall statistics followed by one section per component.
        writeCoverageStatsSection(report, "Global Coverage Statistics", global_stats);

        for (std::size_t component_id = 0; component_id < component_count; ++component_id) {
            writeCoverageStatsSection(
                report,
                component_names[component_id] + " Coverage Statistics",
                component_stats[component_id]
            );
        }
    }
} // namespace output