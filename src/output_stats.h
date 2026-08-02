#pragma once

#include "cdx_types.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace output {

    struct CoverageStats
    {
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
        const std::vector<cdx::Coverage>& coverage
    );


    void writeStatsReportQuery(
        const std::filesystem::path& output_txt,
        const std::map<std::string, std::uint64_t>& mapping,
        const std::vector<cdx::Coverage>& coverage,
        const std::string& component_name
    );

    void writeStatsReportGlobal(
        const std::filesystem::path& output_txt,
        const std::vector<cdx::Coverage>& flat_bp_cov_table,
        const std::vector<cdx::PosBp>& bp_component_offsets,
        const std::vector<std::string>& component_names,
        const std::map<std::string, std::uint64_t>* mapping_stats = nullptr
    );

} // namespace output