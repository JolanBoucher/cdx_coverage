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

    [[nodiscard]]
    CoverageStats computeCoverageStats(
        const std::vector<cdx::Coverage> &coverage
    );


    void writeStatsReportQuery(
        const std::filesystem::path &output_path,
        const GamMappingStats &mapping_stats,
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name
    );

    void writeStatsReportGlobal(
        const std::filesystem::path &output_path,
        const std::vector<cdx::Coverage> &coverage,
        const std::vector<cdx::PosBp> &component_offsets,
        const std::vector<std::string> &component_names,
        const GamMappingStats &mapping_stats,
        int threads
    );
} // namespace output
