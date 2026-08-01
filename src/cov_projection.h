#pragma once

/**
 * @file cov_projection.h
 * @brief Coverage projection and base-pair expansion routines.
 */

#include "cdx_types.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>


/** @brief Projects query coverage from nid-space into local idx-space. */
[[nodiscard]] std::vector<cdx::Coverage> projectCov2IdxQuery(
    const std::vector<cdx::Coverage> &cov_table,
    const std::vector<cdx::Idx> &nid2idx,
    std::size_t component_size
);

/** @brief Expands node-level query coverage into bp-space. */
[[nodiscard]] std::vector<cdx::Coverage> expandPosCovQuery(
    const std::vector<cdx::Coverage> &idx_cov_table,
    const std::vector<cdx::PosBp> &idx2bp
);

/** @brief Projects global coverage from nid-space into flat idx-space. */
[[nodiscard]] std::vector<cdx::Coverage> projectCov2IdxGlobal(
    const std::vector<cdx::Coverage> &cov_table,
    const std::vector<cdx::FlatIdx> &nid2flat_idx,
    const std::vector<cdx::RecordCount> &component_offsets
);

/**
 * @brief Expands flattened node coverage into flattened bp-space.
 *
 * Returns:
 *  - bp coverage table
 *  - component bp offsets
 */
[[nodiscard]] std::pair<std::vector<cdx::Coverage>, std::vector<cdx::PosBp> >
expandPosCovGlobal(
    const std::vector<cdx::Coverage> &flat_idx_cov_table,
    const std::vector<cdx::RecordCount> &component_offsets,
    const std::vector<cdx::PosBp> &idx2bp,
    const std::vector<cdx::RecordCount> &idx2bp_offsets
);

/** @brief Masks positions outside a query interval. */
[[nodiscard]] std::vector<cdx::Coverage> trimCoverageToQuery(
    const std::vector<cdx::Coverage> &bp_cov_table,
    const std::pair<cdx::PosBp, cdx::PosBp> &query_bound
);
