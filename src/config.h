#pragma once
#ifndef CDX_COVERAGE_CONFIG_H
#define CDX_COVERAGE_CONFIG_H

#include <cstdint>
#include <string_view>

namespace cfg {
    // sentinels
    inline constexpr std::uint32_t NOT_IN_COMPO = 0xFFFFFFFF;
    inline constexpr std::uint32_t NOT_IN_QUERY = 0xFFFFFFFF - 1;
    inline constexpr std::uint64_t INVALID_FLAT_IDX = 0xFFFFFFFFFFFFFFFF;
    inline constexpr std::uint32_t INVALID_NODE = 0xFFFFFFFF;

    // default file name
    inline constexpr std::string_view NAME_TSV_FILE = "coverage_profile";
    inline constexpr std:: string_view NAME_STATS_FILE = "coverage_stats";
    inline constexpr std::string_view NAME_GRAPH_FILE = "coverage_graph";

}
#endif //CDX_COVERAGE_CONFIG_H