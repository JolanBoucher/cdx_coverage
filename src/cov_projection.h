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
#include "coverage_gaps.h"

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


// ---------------------------------------------------------------------------
// Base-pair-precision coverage gaps.
//
// The functions below are the base-pair-precision counterpart to the
// node-level pipeline above. They operate strictly *after* the expand step
// (expandPosCovQuery/expandPosCovGlobal) and *before* any query-specific
// trimming (trimCoverageToQuery, the circular backend's own extraction,
// etc.), so a single application covers the linear, circular, range-
// restricted and whole-graph query modes uniformly - none of them need to
// know that base-pair correction exists.
// ---------------------------------------------------------------------------

/**
 * @brief Builds a per-node base-pair length lookup table for a single-
 *        component query, indexed identically to `QueryData::node_coverage`
 *        (i.e. by raw node-ID offset from `nid_min`, the same space
 *        `process_gam`'s `target` vector uses) rather than by topological
 *        index.
 *
 * This is the table `process_gam` needs (as its `node_lengths` parameter)
 * to compute base-pair-precision coverage gaps: node length isn't otherwise
 * available in that raw nid-offset coordinate space, only via `idx2bp`
 * (indexed by topological index) combined with `nid2idx` (which translates
 * between the two).
 *
 * @param nid2idx Maps a node identifier offset to its local topological index.
 * @param idx2bp Cumulative base-pair prefix-sum array (size = component
 *        node count + 1); consecutive entries give each node's length.
 * @return Per-node lengths, same size and indexing as @p nid2idx. Nodes
 *         mapped to a sentinel index (`cfg::NOT_IN_QUERY`/`cfg::NOT_IN_COMPO`)
 *         get a length of 0, since `process_gam` never looks at their length
 *         anyway (they are filtered out by `valid_nodes` before gap detection
 *         would run).
 *
 * @throws std::invalid_argument If idx2bp has fewer than 2 entries.
 * @throws std::out_of_range If a non-sentinel local index exceeds the
 *         component's node count implied by idx2bp.
 * @throws std::runtime_error If idx2bp is non-monotonic for some node.
 */
[[nodiscard]] std::vector<cdx::SeqLen> buildNodeLengthsQuery(
    const std::vector<cdx::Idx> &nid2idx,
    const std::vector<cdx::PosBp> &idx2bp
);

/**
 * @brief Global-graph counterpart to buildNodeLengthsQuery(), indexed by raw
 *        node-ID offset from the graph's minimum node ID (matching
 *        `GlobalData::node_coverage`'s coordinate space).
 *
 * @param nid2flat_idx Maps a node identifier offset to its global flattened index.
 * @param component_offsets Cumulative node-count offset of each component
 *        (size = component count + 1), partitioning the flat index space.
 * @param idx2bp Concatenated per-component base-pair prefix-sum array.
 * @param idx2bp_offsets Per-component starting offset into idx2bp (size =
 *        component count + 1).
 * @return Per-node lengths, same size and indexing as @p nid2flat_idx. Nodes
 *         mapped to `cfg::INVALID_FLAT_IDX` get a length of 0.
 *
 * @throws std::invalid_argument If the component/idx2bp boundary tables are
 *         inconsistently sized.
 * @throws std::out_of_range If a flat index or its derived idx2bp position
 *         falls outside the tables' bounds.
 */
[[nodiscard]] std::vector<cdx::SeqLen> buildNodeLengthsGlobal(
    const std::vector<cdx::FlatIdx> &nid2flat_idx,
    const std::vector<cdx::RecordCount> &component_offsets,
    const std::vector<cdx::PosBp> &idx2bp,
    const std::vector<cdx::RecordCount> &idx2bp_offsets
);

/**
 * @brief Subtracts a set of base-pair-precision coverage gaps from an
 *        already-expanded, single-component base-pair coverage array.
 *
 * For each gap, the node it refers to (identified by raw nid-offset) is
 * translated into a topological index via @p nid2idx, then into its
 * component-local base-pair start via @p idx2bp; the gap's node-local
 * [start, end) range is offset by that and decremented by 1 in @p
 * bp_coverage. Gaps referring to a node outside the active
 * query/component (a sentinel `nid2idx` entry) are silently skipped - they
 * describe a read that touched a node process_gam already excluded from
 * `target`, so there is nothing in @p bp_coverage to correct.
 *
 * @param bp_coverage Component-wide base-pair coverage array, modified in place.
 * @param gaps Coverage gaps collected by process_gam() in CoveragePrecision::Base mode.
 * @param nid2idx Maps a node identifier offset to its local topological index.
 * @param idx2bp Cumulative base-pair prefix-sum array for the component.
 *
 * @throws std::out_of_range If a gap's translated position falls outside
 *         @p bp_coverage's bounds, or its nid_offset exceeds @p nid2idx.
 * @throws std::logic_error If applying a gap would decrement a position
 *         already at 0 - this can only happen if the gap geometry produced
 *         by process_gam is inconsistent with the node-level coverage it is
 *         meant to refine (a bug), not from normal/malformed GAM content.
 */
void applyBpGapsQuery(
    std::vector<cdx::Coverage> &bp_coverage,
    const std::vector<BpGap> &gaps,
    const std::vector<cdx::Idx> &nid2idx,
    const std::vector<cdx::PosBp> &idx2bp
);

/**
 * @brief Global-graph counterpart to applyBpGapsQuery(), operating on the
 *        flattened multi-component base-pair coverage array produced by
 *        expandPosCovGlobal().
 *
 * Each gap's node is translated: raw nid-offset -> global flat index (via
 * @p nid2flat_idx) -> owning component (located within @p component_offsets
 * by binary search) -> component-local base-pair start (via @p idx2bp /
 * @p idx2bp_offsets) -> final flattened position (shifted by @p
 * component_bp_offsets, the per-component base-pair boundaries returned
 * alongside expandPosCovGlobal()'s coverage array). Gaps referring to a node
 * with no flat index (`cfg::INVALID_FLAT_IDX`) are silently skipped, for the
 * same reason as in applyBpGapsQuery().
 *
 * @param flat_bp_coverage Flattened whole-graph base-pair coverage array, modified in place.
 * @param gaps Coverage gaps collected by process_gam() in CoveragePrecision::Base mode.
 * @param nid2flat_idx Maps a node identifier offset to its global flattened index.
 * @param component_offsets Cumulative node-count offset of each component.
 * @param idx2bp Concatenated per-component base-pair prefix-sum array.
 * @param idx2bp_offsets Per-component starting offset into idx2bp.
 * @param component_bp_offsets Per-component starting offset into the final
 *        flattened base-pair array (the second element of
 *        expandPosCovGlobal()'s return value).
 *
 * @throws std::invalid_argument If the component/idx2bp/component_bp_offsets
 *         boundary tables are inconsistently sized.
 * @throws std::out_of_range If a gap's translated position falls outside
 *         @p flat_bp_coverage's bounds, its nid_offset exceeds
 *         @p nid2flat_idx, or its flat index cannot be located within
 *         @p component_offsets.
 * @throws std::logic_error Same invariant as applyBpGapsQuery().
 */
void applyBpGapsGlobal(
    std::vector<cdx::Coverage> &flat_bp_coverage,
    const std::vector<BpGap> &gaps,
    const std::vector<cdx::FlatIdx> &nid2flat_idx,
    const std::vector<cdx::RecordCount> &component_offsets,
    const std::vector<cdx::PosBp> &idx2bp,
    const std::vector<cdx::RecordCount> &idx2bp_offsets,
    const std::vector<cdx::PosBp> &component_bp_offsets
);
