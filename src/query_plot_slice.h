/**
 * @file query_plot_slice.h
 * @brief Small helper extracted from main.cpp's runQueryPipeline(): trims a full-component
 *        base-pair coverage vector down to the exact requested query interval, for the linear
 *        (Cairo) graph backend only.
 *
 * writeLinearPlotQuery() plots whatever vector it is handed as spanning [offset, offset +
 * vector.size()); to make the graph's X-axis match the user-requested query range exactly
 * (e.g. `-q "chr1 0:100"` should produce an axis from 0 to 100, not from the component's
 * start to its end), the component-wide bp_coverage vector produced by
 * expandPosCovQuery()/trimCoverageToQuery() must first be sliced down to [query_start,
 * query_end]. This is pure computation with no I/O, so unlike the rest of main.cpp it can be
 * pulled out into its own testable module (see query_plot_slice_test.cpp).
 *
 * The circular backend does NOT use this: writeCircularPlotQuery() needs the full,
 * un-sliced component coverage (circular smoothing/origin-crossing queries need context on
 * both sides of the query boundary) and extracts its own traversal internally.
 */

#ifndef CDX_COVERAGE_QUERY_PLOT_SLICE_H
#define CDX_COVERAGE_QUERY_PLOT_SLICE_H

#include "cdx_types.h"

#include <utility>
#include <vector>

namespace output {
    /**
     * @brief Extracts the coverage sub-vector to plot for a linear query graph, so the
     *        rendered X-axis spans exactly [query_start, query_end] rather than the whole
     *        component.
     *
     * @param bp_coverage Full-component coverage vector, indexed by base-pair position
     *        (0-based, one entry per position in the component).
     * @param query_range_bp Inclusive [query_start, query_end] base-pair interval requested by
     *        the user.
     * @return std::vector<cdx::Coverage> The sub-vector bp_coverage[query_start..query_end]
     *         (inclusive), preserved in original order/values (no filtering or transformation).
     *
     * @throws std::invalid_argument If query_start > query_end (an origin-crossing range) -
     *         a linear graph cannot represent that; only the circular backend can.
     * @throws std::out_of_range If query_end is beyond bp_coverage's bounds.
     */
    [[nodiscard]]
    std::vector<cdx::Coverage> sliceLinearQueryCoverage(
        const std::vector<cdx::Coverage> &bp_coverage,
        std::pair<cdx::PosBp, cdx::PosBp> query_range_bp
    );
} // namespace output

#endif // CDX_COVERAGE_QUERY_PLOT_SLICE_H
