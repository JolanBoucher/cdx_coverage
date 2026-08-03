#pragma once
#ifndef CDX_COVERAGE_CONFIG_H
#define CDX_COVERAGE_CONFIG_H

#include <cstdint>
#include <iostream>
#include <string_view>

namespace cfg {
    // sentinels
    inline constexpr std::uint32_t NOT_IN_COMPO = 0xFFFFFFFF;
    inline constexpr std::uint32_t NOT_IN_QUERY = 0xFFFFFFFF - 1;
    inline constexpr std::uint64_t INVALID_FLAT_IDX = 0xFFFFFFFFFFFFFFFF;
    inline constexpr std::uint32_t INVALID_NODE = 0xFFFFFFFF;

    // default file name
    inline constexpr std::string_view NAME_TSV_FILE = "coverage_profile.txt";
    inline constexpr std:: string_view NAME_STATS_FILE = "coverage_stats.tsv";
    inline constexpr std::string_view NAME_GRAPH_FILE = "coverage_graph.png";

    class ScopedTimer {
    public:
        explicit ScopedTimer(std::string name)
            : name_(std::move(name)),
              start_(std::chrono::steady_clock::now()) {}

        ~ScopedTimer() {
            const auto end = std::chrono::steady_clock::now();
            const double seconds = std::chrono::duration<double>(end - start_).count();
            std::cout << "[TIME] " << name_ << ": " << seconds << " s\n";
        }

    private:
        std::string name_;
        std::chrono::steady_clock::time_point start_;
    };
}
#endif //CDX_COVERAGE_CONFIG_H