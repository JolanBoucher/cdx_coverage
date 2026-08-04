#pragma once
#ifndef CDX_COVERAGE_CONFIG_H
#define CDX_COVERAGE_CONFIG_H

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <iomanip>

namespace cfg {
    // sentinels
    inline constexpr std::uint32_t NOT_IN_COMPO = 0xFFFFFFFF;
    inline constexpr std::uint32_t NOT_IN_QUERY = 0xFFFFFFFF - 1;
    inline constexpr std::uint64_t INVALID_FLAT_IDX = 0xFFFFFFFFFFFFFFFF;
    inline constexpr std::uint32_t INVALID_NODE = 0xFFFFFFFF;

    // default file name
    inline constexpr std::string_view NAME_TSV_FILE = "coverage_profile.tsv";
    inline constexpr std::string_view NAME_STATS_FILE = "coverage_stats.txt";
    inline constexpr std::string_view NAME_GRAPH_FILE = "coverage_graph.png";

    /** @brief RAII timer measuring execution duration with steady_clock and exception tracking. */
    class ScopedTimer {
public:
    using Clock = std::chrono::steady_clock;

    /**
     * @brief Construct timer and capture current clock time and active exception count.
     * @param step_name Label displayed upon timer completion or failure.
     */
    explicit ScopedTimer(std::string step_name)
        : name_(std::move(step_name)),
          start_(Clock::now()),
          uncaught_on_entry_(std::uncaught_exceptions()) {}
    /**
     * @brief Dynamically updates the step label displayed upon timer destruction.
     * @param new_name The new string description for the timed step.
     */
    void update_name(std::string new_name) {
        name_ = std::move(new_name);
    }

    /**
     * @brief Destructor that computes elapsed time and logs completion or failure status.
     * @note Marked noexcept to prevent throwing during unwinding.
     */
    ~ScopedTimer() noexcept {
        try {
            const auto end = Clock::now();
            const double seconds = std::chrono::duration<double>(end - start_).count();
            const bool failed = std::uncaught_exceptions() > uncaught_on_entry_;
            print_step_time(name_, seconds, failed);
        } catch (...) {
            // Guarantee noexcept behavior in destructor
        }
    }

    // Prevent copy/move to avoid duplicate measurement prints
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
    ScopedTimer(ScopedTimer&&) = delete;
    ScopedTimer& operator=(ScopedTimer&&) = delete;

private:
    std::string name_;        ///< Name of the pipeline step being measured.
    Clock::time_point start_; ///< Start timestamp recorded on construction.
    int uncaught_on_entry_;   ///< Number of uncaught exceptions present at construction.

    /**
     * @brief Formats and prints the timing output to std::cerr.
     * @param step_name Label of the executed step.
     * @param seconds Elapsed duration in seconds.
     * @param failed True if step failed due to an uncaught exception.
     */
    static void print_step_time(const std::string& step_name, const double seconds, const bool failed) {
        std::cerr << "  - " << std::left << std::setw(50) << step_name
                  << (failed ? " Failed after   " : " Completed in ")
                  << std::right << std::setw(7) << std::fixed
                  << std::setprecision(4) << seconds << " s\n";
    }
};

   // tiny formating function that transform interger like this 1534430 to this 1,534,430
    template<typename Integer> [[nodiscard]]
    std::string formatInteger(const Integer value) {

        std::string text = std::to_string(value);
        const std::size_t sign_offset = !text.empty() && text.front() == '-' ? 1 : 0;

        for (std::size_t position = text.size(); position > sign_offset + 3; position -= 3)
            text.insert(position - 3, 1, ',');

        return text;
    }
}
#endif //CDX_COVERAGE_CONFIG_H
