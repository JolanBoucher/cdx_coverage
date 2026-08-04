/**
* @file output_stats.h
 * @brief Coverage-statistics computation and report-generation utilities.
 *
 * This module provides data structures and functions for computing summary
 * statistics from per-base coverage data and generating human-readable
 * coverage reports. It supports both single-component and multi-component
 * coverage tables, including calculation of coverage breadth, mean,
 * median, standard deviation, coefficient of variation, quartiles, and
 * coverage extrema.
 *
 * Statistics are computed using either a dense-histogram approach or an
 * nth_element-based fallback for large coverage ranges. Report-generation
 * utilities produce formatted summaries of mapping and coverage metrics for
 * downstream analysis and quality control.
 */

#pragma once

#include "cdx_types.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "gam_io.h"

namespace output {

    /**
    * @brief Summary statistics describing coverage over a genomic region.
    *
    * Stores coverage breadth, central-tendency metrics, dispersion measures,
    * quartiles, extrema, and position counts computed from a coverage dataset.
    */
    struct CoverageStats {
        std::size_t region_length = 0;
        std::size_t covered_positions = 0;
        double breadth = 0.0;
        double mean = 0.0;
        double median = 0.0;
        double stddev = 0.0;
        double cv = 0.0;
        double q1 = 0.0;
        double q3 = 0.0;
        cdx::Coverage min = 0;
        cdx::Coverage max = 0;
    };

    /**
    * @brief Compute coverage statistics for a coverage vector.
    *
    * Accumulates coverage information and selects the most appropriate
    * statistics-computation strategy based on the observed coverage range.
    *
    * For typical coverage ranges, statistics are computed from a dense
    * histogram for efficient quantile calculation. If the coverage range is
    * too large to represent efficiently as a dense histogram, the function
    * falls back to a std::nth_element-based implementation.
    *
    * @param coverage Per-position coverage values.
    * @return Computed coverage statistics.
    */
    [[nodiscard]]
    CoverageStats computeCoverageStats(
        const std::vector<cdx::Coverage> &coverage
    );


    /**
    * @brief Generate a coverage report for a single query component.
    *
    * Validates mapping statistics, computes coverage statistics from the
    * supplied coverage vector, and writes a human-readable report containing
    * mapping metrics and coverage summary statistics.
    *
    * The report includes read counts, mapping percentages, coverage breadth,
    * central-tendency statistics, dispersion metrics, quartiles, and coverage
    * extrema for the specified component.
    *
    * @param output_path Path of the report file to create.
    * @param mapping_stats Read-mapping statistics associated with the query.
    * @param coverage Per-position coverage values for the query component.
    * @param component_name Name of the query component being reported.
    *
    * @throws std::invalid_argument If the mapping statistics are internally inconsistent.
    * @throws std::runtime_error If the report file cannot be opened.
    */
    void writeStatsReportQuery(
        const std::filesystem::path &output_path,
        const GamMappingStats &mapping_stats,
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name
    );

    /**
     * @brief Generate a coverage report for a multi-component coverage table.
     *
     * Computes and writes both global and per-component coverage statistics for
     * a flattened coverage table. Component boundaries are defined by
     * component_offsets, which partition the coverage vector into contiguous
     * component regions described by component_names.
     *
     * Coverage statistics are computed independently for each component and,
     * when possible, merged into a global accumulator to efficiently derive
     * overall statistics. Components whose coverage range exceeds
     * MAX_DENSE_HISTOGRAM are processed using the nth_element-based fallback
     * implementation.
     *
     * Component statistics may be computed in parallel when OpenMP support is
     * available.
     *
     * @param output_path Path of the report file to create.
     * @param coverage Flattened per-position coverage table.
     * @param component_offsets Component boundaries within the flattened
     *        coverage table. Must begin at 0 and terminate at coverage.size().
     * @param component_names Names associated with each component interval.
     * @param mapping_stats Read-mapping statistics to include in the report.
     * @param threads Number of worker threads to use when OpenMP is enabled.
     *
     * @throws std::invalid_argument If component boundaries or component names
     *         are inconsistent with the coverage table.
     * @throws std::runtime_error If the report file cannot be opened.
     */
    void writeStatsReportGlobal(
        const std::filesystem::path &output_path,
        const std::vector<cdx::Coverage> &coverage,
        const std::vector<cdx::PosBp> &component_offsets,
        const std::vector<std::string> &component_names,
        const GamMappingStats &mapping_stats,
        int threads
    );
} // namespace output
