#pragma once

/**
 * @file cov_projection.h
 * @brief Coverage projection and base-pair expansion algorithms.
 *
 * Provides utilities for transforming coverage information between the
 * different coordinate systems used throughout the project:
 *   - relative node coverage
 *   - component-local index space
 *   - graph-wide flattened index space
 *   - base-pair level coverage
 *
 * The routines in this file are responsible for projecting node-associated
 * coverage values onto component or global graph layouts, expanding node
 * coverages into per-base-pair depth arrays, concatenating component-level
 * coverage tracks, and trimming results to query-specific genomic intervals.
 */

#include "cdx_types.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>


/**
 * @brief Maps relative node coverage values into a component's local topological index space.
 *
 * Takes node-level coverage indexed by relative offset and redistributes values into a local
 * array ordered by component rank (`local_idx`). Unmapped nodes or nodes matching sentinels are skipped.
 */
[[nodiscard]] std::vector<cdx::Coverage> projectCov2IdxQuery(
    const std::vector<cdx::Coverage> &cov_table,
    const std::vector<cdx::Idx> &nid2idx,
    std::size_t component_size
);


/**
 * @brief Expands node-level coverage into per-base-pair depth across a single component.
 *
 * Uses cumulative prefix-sum array `idx2bp` to calculate each node's [start, end) base-pair interval,
 * filling those intervals with the node's assigned depth.
 */
[[nodiscard]] std::vector<cdx::Coverage> expandPosCovQuery(
    const std::vector<cdx::Coverage> &idx_cov_table,
    const std::vector<cdx::PosBp> &idx2bp
);

/**
 * @brief Projects node coverage across the whole graph into a single flattened index array.
 *
 * Maps relative node offsets to global flat indices (`nid2flat_idx`). Nodes unmapped in the graph
 * or assigned `cfg::INVALID_FLAT_IDX` are skipped.
 */
[[nodiscard]] std::vector<cdx::Coverage> projectCov2IdxGlobal(
    const std::vector<cdx::Coverage> &cov_table,
    const std::vector<cdx::FlatIdx> &nid2flat_idx,
    const std::vector<cdx::RecordCount> &component_offsets
);


/**
 * @brief Concatenates per-component base-pair coverages into a single graph-wide buffer.
 *
 * Sequentially expands node coverage entries across all components, keeping track of global base-pair
 * boundary offsets.
 *
 * @return std::pair containing:
 *         - first:  Flattened base-pair coverage array for the whole graph.
 *         - second: Cumulative base-pair offset boundaries per component (size = num_components + 1).
 */
[[nodiscard]] std::pair<std::vector<cdx::Coverage>, std::vector<cdx::PosBp> >
expandPosCovGlobal(
    const std::vector<cdx::Coverage> &flat_idx_cov_table,
    const std::vector<cdx::RecordCount> &component_offsets,
    const std::vector<cdx::PosBp> &idx2bp,
    const std::vector<cdx::RecordCount> &idx2bp_offsets
);

/**
 * @brief Masks out-of-query regions of a base-pair coverage array using sentinel values.
 *
 * Preserves total vector dimensions while setting positions outside the requested query
 * interval to `cfg::NOT_IN_QUERY` (if currently less than the sentinel value). Supports both standard
 * [start, end] and inverted range conditions.
 */
[[nodiscard]] std::vector<cdx::Coverage> trimCoverageToQuery(
    const std::vector<cdx::Coverage> &bp_cov_table,
    const std::pair<cdx::PosBp, cdx::PosBp> &query_bound
);
