/**
 * @file output_stats.cpp
 * @brief Coverage statistics calculation and text report generation.
 */

#include "output_stats.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.h"

namespace {
    /**
     * @brief Computes a linearly interpolated percentile on sorted values.
     */
    [[nodiscard]]
    double linearPercentile(
        const std::vector<cdx::Coverage> &sorted_coverage,
        const double percentile
    ) {
        if (sorted_coverage.empty()) {
            return 0.0;
        }
        if (sorted_coverage.size() == 1) {
            return sorted_coverage.front();
        }

        const double position = percentile * static_cast<double>(sorted_coverage.size() - 1);
        const auto lower_index = static_cast<std::size_t>(std::floor(position));
        const auto upper_index = static_cast<std::size_t>(std::ceil(position));
        const double lower_value = sorted_coverage[lower_index];

        if (lower_index == upper_index) {
            return lower_value;
        }

        const double upper_value = sorted_coverage[upper_index];
        const double interpolation_weight = position - static_cast<double>(lower_index);

        return lower_value + interpolation_weight * (upper_value - lower_value);
    }

    /**
     * @brief Add comma thousands separators to an unsigned integer.
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
     * @brief Writes one formatted coverage-statistics section.
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
} // anonymous namespace

namespace output {
    CoverageStats computeCoverageStats(
        const std::vector<cdx::Coverage> &coverage
    ) {
        CoverageStats stats;
        stats.region_length = coverage.size();

        if (coverage.empty()) {
            return stats;
        }

        stats.min = coverage.front();
        stats.max = coverage.front();

        double running_mean = 0.0;
        double m2 = 0.0;
        std::size_t processed_count = 0;

        for (const cdx::Coverage value: coverage) {
            ++processed_count;
            if (value != 0) {
                ++stats.covered_positions;
            }

            if (value < stats.min) stats.min = value;
            if (value > stats.max) stats.max = value;

            const double numeric_value = value;
            const double delta = numeric_value - running_mean;
            running_mean += delta / static_cast<double>(processed_count);
            const double updated_delta = numeric_value - running_mean;

            m2 += delta * updated_delta;
        }

        stats.mean = running_mean;

        stats.breadth = 100.0 * static_cast<double>(stats.covered_positions) /
                        static_cast<double>(stats.region_length);

        if (stats.region_length > 0) {
            const double population_variance = m2 / static_cast<double>(stats.region_length);
            stats.stddev = std::sqrt(std::max(population_variance, 0.0));
        }

        stats.cv = stats.mean > 0.0 ? stats.stddev / stats.mean : 0.0;

        // Tri simple et lisible
        std::vector<cdx::Coverage> sorted_coverage = coverage;
        std::sort(sorted_coverage.begin(), sorted_coverage.end());

        stats.q1 = linearPercentile(sorted_coverage, 0.25);
        stats.median = linearPercentile(sorted_coverage, 0.50);
        stats.q3 = linearPercentile(sorted_coverage, 0.75);

        return stats;
    }

    void writeStatsReportQuery(
        const std::filesystem::path &output_txt,
        const std::map<std::string, std::uint64_t> &mapping,
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name
    ) {
        // 1. Validation unique des clés requises
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

        // 2. Vérifications de cohérence
        if (mapped > total) throw std::invalid_argument("Mapped reads cannot exceed total reads.");
        if (unmapped > total) throw std::invalid_argument("Unmapped reads cannot exceed total reads.");
        if (mapped_to_query > mapped) throw std::invalid_argument("Reads mapped to query cannot exceed mapped reads.");
        if (mapped + unmapped != total) throw std::invalid_argument(
            "Mapped and unmapped reads must sum to total reads.");

        // 3. Filtrage du vecteur de couverture
        std::vector<cdx::Coverage> valid_coverage;
        valid_coverage.reserve(coverage.size());
        for (const auto value: coverage) {
            if (value < cfg::NOT_IN_QUERY) {
                valid_coverage.push_back(value);
            }
        }

        const CoverageStats coverage_stats = computeCoverageStats(valid_coverage);

        // 4. Calcul des pourcentages
        const double pct_mapped = total ? 100.0 * static_cast<double>(mapped) / static_cast<double>(total) : 0.0;
        const double pct_query_total = total
                                           ? 100.0 * static_cast<double>(mapped_to_query) / static_cast<double>(total)
                                           : 0.0;
        const double pct_query_mapped = mapped
                                            ? 100.0 * static_cast<double>(mapped_to_query) / static_cast<double>(mapped)
                                            : 0.0;
        const double pct_unmapped = total ? 100.0 * static_cast<double>(unmapped) / static_cast<double>(total) : 0.0;

        // 5. Écriture unique du fichier
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
        // Validation stricte de la cohérence des noms et des offsets
        if (bp_component_offsets.size() < 2) {
            throw std::invalid_argument("bp_component_offsets must contain at least two boundaries.");
        }
        if (component_names.size() + 1 != bp_component_offsets.size()) {
            throw std::invalid_argument("Component name count does not match component count.");
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

        // Statistiques globales
        std::vector<cdx::Coverage> global_coverage;
        global_coverage.reserve(flat_bp_cov_table.size());
        for (const cdx::Coverage value: flat_bp_cov_table) {
            if (value < cfg::NOT_IN_QUERY) {
                global_coverage.push_back(value);
            }
        }

        const CoverageStats global_stats = computeCoverageStats(global_coverage);
        writeCoverageStatsSection(report, "Global Coverage Statistics", global_stats);

        // Statistiques par composante avec réutilisation du buffer mémoire
        std::vector<cdx::Coverage> component_cov;

        for (std::size_t component_id = 0; component_id + 1 < bp_component_offsets.size(); ++component_id) {
            const std::string &component_name = component_names[component_id];
            const auto component_start = static_cast<std::size_t>(bp_component_offsets[component_id]);
            const auto component_end = static_cast<std::size_t>(bp_component_offsets[component_id + 1]);

#ifndef NDEBUG
            if (component_start > component_end || component_end > flat_bp_cov_table.size()) {
                throw std::out_of_range("Invalid component offset boundaries.");
            }
#endif

            component_cov.clear();
            component_cov.reserve(component_end - component_start);

            for (std::size_t pos = component_start; pos < component_end; ++pos) {
                const cdx::Coverage value = flat_bp_cov_table[pos];
                if (value < cfg::NOT_IN_QUERY) {
                    component_cov.push_back(value);
                }
            }

            const CoverageStats component_stats = computeCoverageStats(component_cov);

            writeCoverageStatsSection(
                report,
                component_name + " Coverage Statistics ",
                component_stats
            );
        }
    }
} // namespace output
