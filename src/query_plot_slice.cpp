/**
 * @file query_plot_slice.cpp
 * @brief Implementation of sliceLinearQueryCoverage() - see query_plot_slice.h.
 */

#include "query_plot_slice.h"

#include <stdexcept>

namespace output {
    std::vector<cdx::Coverage> sliceLinearQueryCoverage(
        const std::vector<cdx::Coverage> &bp_coverage,
        const std::pair<cdx::PosBp, cdx::PosBp> query_range_bp
    ) {
        const auto [query_start_bp, query_end_bp] = query_range_bp;

        if (query_start_bp > query_end_bp) {
            throw std::invalid_argument(
                "A linear graph cannot represent an origin-crossing query."
            );
        }
        if (query_end_bp >= bp_coverage.size()) {
            throw std::out_of_range("Query range exceeds component coverage bounds.");
        }

        return std::vector<cdx::Coverage>(
            bp_coverage.begin() + static_cast<std::ptrdiff_t>(query_start_bp),
            bp_coverage.begin() + static_cast<std::ptrdiff_t>(query_end_bp) + 1
        );
    }
} // namespace output
