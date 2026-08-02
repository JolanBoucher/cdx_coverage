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

    /**
     * @brief Intermediate accumulator aggregating frequencies and statistical moments.
     */
    struct CoverageAccumulator {
        std::size_t region_length = 0;
        std::size_t covered_positions = 0;

        std::uint64_t sum = 0;
        double sum_squares = 0.0;

        cdx::Coverage min = std::numeric_limits<cdx::Coverage>::max();
        cdx::Coverage max = 0;

        bool requires_fallback = false;

        std::vector<std::size_t> histogram;

        /**
         * @brief Merges a component accumulator into the global accumulator.
         */
        void merge(const CoverageAccumulator &other) {
            if (other.region_length == 0) return;

#ifndef NDEBUG
            if (other.sum > std::numeric_limits<std::uint64_t>::max() - sum) {
                throw std::overflow_error("Coverage sum exceeds uint64_t capacity.");
            }
#endif

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
     * @brief Formats an unsigned integer with thousands separators (commas).
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
     * @brief Writes one formatted coverage-statistics section to an output stream.
     */
    void writeCoverageStatsSection(
        std::ostream &stream,
        const std::string &title,
        const output::CoverageStats &stats
    ) {
        stream << "[ " << title << " ]\n"
                << std::string(70, '-') << '\n';

        stream << std::left << std::setw(30) << "Indexed positions" << " : "
                << formatInteger(stats.region_length) << '\n';

        stream << std::left << std::setw(30) << "Covered positions" << " : "
                << formatInteger(stats.covered_positions) << '\n';

        stream << std::fixed << std::setprecision(2);
        stream << std::left << std::setw(30) << "Breadth of coverage" << " : " << stats.breadth << "%\n";
        stream << std::left << std::setw(30) << "Mean coverage" << " : " << stats.mean << "x\n";
        stream << std::left << std::setw(30) << "Median coverage" << " : " << stats.median << "x\n";
        stream << std::left << std::setw(30) << "Standard deviation" << " : " << stats.stddev << '\n';
        stream << std::setprecision(4) << std::left << std::setw(30)
                << "Coefficient of variation" << " : " << stats.cv << '\n';
        stream << std::setprecision(2) << std::left << std::setw(30) << "Q1" << " : " << stats.q1 << "x\n";
        stream << std::left << std::setw(30) << "Q3" << " : " << stats.q3 << "x\n";

        stream << std::defaultfloat << std::left << std::setw(30) << "Minimum coverage" << " : " << stats.min << "x\n";
        stream << std::left << std::setw(30) << "Maximum coverage" << " : " << stats.max << "x\n";
    }

    /**
     * @brief Fallback calculation using std::nth_element if max > MAX_DENSE_HISTOGRAM.
     */
    template<typename Iterator>
    output::CoverageStats computeCoverageStatsNthElement(Iterator begin, Iterator end) {
        output::CoverageStats stats;
        std::vector<cdx::Coverage> valid_cov;

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

        stats.q1 = get_quantile(0.25);
        stats.median = get_quantile(0.50);
        stats.q3 = get_quantile(0.75);

        return stats;
    }

    /**
     * @brief Computes final stats (Mean, StdDev, Quantiles) directly from an accumulator.
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

        const auto r_l_q1 = static_cast<std::size_t>(std::floor(rank_q1));
        const auto r_u_q1 = static_cast<std::size_t>(std::ceil(rank_q1));
        const auto r_l_med = static_cast<std::size_t>(std::floor(rank_med));
        const auto r_u_med = static_cast<std::size_t>(std::ceil(rank_med));
        const auto r_l_q3 = static_cast<std::size_t>(std::floor(rank_q3));
        const auto r_u_q3 = static_cast<std::size_t>(std::ceil(rank_q3));

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

        auto interp = [](const double rank, const std::size_t r_l, const cdx::Coverage v_l, const cdx::Coverage v_u) {
            const double fraction = rank - static_cast<double>(r_l);
            return static_cast<double>(v_l) + fraction * static_cast<double>(v_u - v_l);
        };

        stats.q1 = interp(rank_q1, r_l_q1, val_l_q1, val_u_q1);
        stats.median = interp(rank_med, r_l_med, val_l_med, val_u_med);
        stats.q3 = interp(rank_q3, r_l_q3, val_l_q3, val_u_q3);

        return stats;
    }

    /**
     * @brief Accumulates a coverage range using two sequential passes.
     */
    template<typename Iterator>
    [[nodiscard]]
    CoverageAccumulator accumulateRange(
        Iterator begin,
        Iterator end
    ) {
        CoverageAccumulator acc;

        // Pass 1: determine the valid coverage range.
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

        // Pass 2: histogram and numerical moments.
        for (auto iterator = begin; iterator != end; ++iterator) {
            const cdx::Coverage value = *iterator;
            if (value >= cfg::NOT_IN_QUERY) {
                continue;
            }

            ++acc.region_length;
            if (value != 0) {
                ++acc.covered_positions;
            }

            ++acc.histogram[static_cast<std::size_t>(value)];

            const auto unsigned_value = static_cast<std::uint64_t>(value);
#ifndef NDEBUG
            if (unsigned_value > std::numeric_limits<std::uint64_t>::max() - acc.sum) {
                throw std::overflow_error("Coverage sum exceeds uint64_t capacity.");
            }
#endif
            acc.sum += unsigned_value;

            const double numeric_value = value;
            acc.sum_squares += numeric_value * numeric_value;
        }

        return acc;
    }
} // anonymous namespace

namespace output {
    CoverageStats computeCoverageStats(
        const std::vector<cdx::Coverage> &coverage
    ) {
        const CoverageAccumulator accumulator = accumulateRange(
            coverage.begin(),
            coverage.end()
        );

        if (accumulator.requires_fallback) {
            return computeCoverageStatsNthElement(
                coverage.begin(),
                coverage.end()
            );
        }

        return finalizeStatsFromAccumulator(accumulator);
    }

    void writeStatsReportQuery(
        const std::filesystem::path &output_txt,
        const std::map<std::string, std::uint64_t> &mapping,
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name
    ) {
        constexpr std::array<const char *, 4> required_keys{"total", "mapped", "mapped_to_query", "unmapped"};
        for (const auto *key: required_keys) {
            if (mapping.find(key) == mapping.end()) {
                throw std::invalid_argument("Missing mapping statistic: " + std::string(key));
            }
        }

        const std::uint64_t total = mapping.at("total");
        const std::uint64_t mapped = mapping.at("mapped");
        const std::uint64_t mapped_to_query = mapping.at("mapped_to_query");
        const std::uint64_t unmapped = mapping.at("unmapped");

        if (mapped > total) throw std::invalid_argument("Mapped reads cannot exceed total reads.");
        if (unmapped > total) throw std::invalid_argument("Unmapped reads cannot exceed total reads.");
        if (mapped_to_query > mapped) throw std::invalid_argument("Reads mapped to query cannot exceed mapped reads.");
        if (mapped + unmapped != total) {
            throw std::invalid_argument("Mapped and unmapped reads must sum to total reads.");
        }

        const CoverageStats coverage_stats = computeCoverageStats(coverage);

        const double pct_mapped = total
                                      ? 100.0 * static_cast<double>(mapped) / static_cast<double>(total)
                                      : 0.0;
        const double pct_query_total = total
                                           ? 100.0 * static_cast<double>(mapped_to_query) / static_cast<double>(total)
                                           : 0.0;
        const double pct_query_mapped = mapped
                                            ? 100.0 * static_cast<double>(mapped_to_query) / static_cast<double>(mapped)
                                            : 0.0;
        const double pct_unmapped = total
                                        ? 100.0 * static_cast<double>(unmapped) / static_cast<double>(total)
                                        : 0.0;

        std::ofstream report(output_txt);
        if (!report) {
            throw std::runtime_error("Unable to open report file: " + output_txt.string());
        }

        report << std::string(70, '=') << '\n';
        report << "COVERAGE REPORT\n";
        report << "Component: " << component_name << '\n';
        report << std::string(70, '=') << "\n\n";

        report << "[ Mapping Statistics ]\n";
        report << std::string(70, '-') << '\n';
        report << std::left << std::setw(30) << "Total reads" << " : " << formatInteger(total) << '\n';
        report << std::left << std::setw(30) << "Mapped reads" << " : " << formatInteger(mapped)
                << " (" << std::fixed << std::setprecision(2) << pct_mapped << "%)\n";
        report << std::left << std::setw(30) << "Reads mapped to query" << " : " << formatInteger(mapped_to_query) <<
                '\n';
        report << std::left << std::setw(30) << "  Percentage of total" << " : " << pct_query_total << "%\n";
        report << std::left << std::setw(30) << "  Percentage of mapped" << " : " << pct_query_mapped << "%\n";
        report << std::left << std::setw(30) << "Unmapped reads" << " : " << formatInteger(unmapped)
                << " (" << pct_unmapped << "%)\n\n";

        writeCoverageStatsSection(report, "Coverage Statistics", coverage_stats);
    }

    void writeStatsReportGlobal(
        const std::filesystem::path &output_txt,
        const std::vector<cdx::Coverage> &flat_bp_cov_table,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names,
        const std::map<std::string, std::uint64_t> *mapping_stats
    ) {
        if (bp_component_offsets.size() < 2) {
            throw std::invalid_argument("bp_component_offsets must contain at least two boundaries.");
        }
        if (component_names.size() + 1 != bp_component_offsets.size()) {
            throw std::invalid_argument("Component name count does not match component count.");
        }

        if (bp_component_offsets.front() != 0) {
            throw std::invalid_argument("The first bp component offset must be zero.");
        }
        if (bp_component_offsets.back() != static_cast<cdx::PosBp>(flat_bp_cov_table.size())) {
            throw std::invalid_argument("The final bp component offset must match the flattened coverage table size.");
        }
        for (std::size_t i = 1; i < bp_component_offsets.size(); ++i) {
            if (bp_component_offsets[i] < bp_component_offsets[i - 1]) {
                throw std::invalid_argument("bp_component_offsets must be non-decreasing.");
            }
        }

        std::ofstream report(output_txt);
        if (!report) {
            throw std::runtime_error("Unable to open report file: " + output_txt.string());
        }

        report << std::string(70, '=') << '\n';
        report << "COVERAGE REPORT\n";
        report << std::string(70, '=') << '\n';

        if (mapping_stats != nullptr) {
            report << "[ Mapping Statistics ]\n";
            report << std::string(70, '-') << '\n';
            for (const auto &[key, value]: *mapping_stats) {
                report << std::left << std::setw(30) << key << " : " << value << '\n';
            }
            report << '\n';
        }

        const std::size_t component_count = bp_component_offsets.size() - 1;
        std::vector<CoverageAccumulator> component_accumulators(component_count);
        std::vector<CoverageStats> component_stats(component_count);
        std::vector<std::uint8_t> component_uses_fallback(component_count, 0);

        {
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
            for (std::ptrdiff_t component_index = 0;
                 component_index < static_cast<std::ptrdiff_t>(component_count);
                 ++component_index) {
                const auto component_id = static_cast<std::size_t>(component_index);
                const auto start_index = static_cast<std::size_t>(bp_component_offsets[component_id]);
                const auto end_index = static_cast<std::size_t>(bp_component_offsets[component_id + 1]);

                const auto begin_iterator
                        = flat_bp_cov_table.begin() + static_cast<std::ptrdiff_t>(start_index);
                const auto end_iterator
                        = flat_bp_cov_table.begin() + static_cast<std::ptrdiff_t>(end_index);

                component_accumulators[component_id] = accumulateRange(
                    begin_iterator,
                    end_iterator
                );

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

        for (std::size_t component_id = 0; component_id < component_count; ++component_id) {
            if (component_uses_fallback[component_id] != 0) {
                global_fallback = true;
                break;
            }
            global_accumulator.merge(component_accumulators[component_id]);
        }

        CoverageStats global_stats;
        if (global_fallback) {
            global_stats = computeCoverageStatsNthElement(
                flat_bp_cov_table.begin(),
                flat_bp_cov_table.end()
            );
        } else {
            global_stats = finalizeStatsFromAccumulator(global_accumulator);
        }

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
