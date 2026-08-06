/**
 * @file gam_io.h
 * @brief Header file defining structures and functions for processing GAM alignment files[cite: 1].
 *
 * This module provides the core parallel engine (`process_gam`) for streaming and analyzing
 * GAM alignments against a CDX graph index, alongside tracking structures (`GamMappingStats`)[cite: 1]
 * for performance metrics and mapping statistics.
 *
 * Two coverage-computation granularities are supported, selected via `CoveragePrecision`
 * (see coverage_precision.h):
 *   - `Node`: the original behaviour - every read touching a node increments that node's
 *     coverage by exactly one. `Edit` submessages are never inspected.
 *   - `Base`: additionally collects `BpGap` entries (see coverage_gaps.h) describing, per
 *     read, which base-pair ranges of a touched node were *not* actually traversed at
 *     base-pair resolution (deletions, and read-boundary under/over-reach). These gaps are
 *     only geometry expressed in raw node-ID-offset space; turning them into a corrected
 *     base-pair coverage array is the job of `applyBpGapsQuery`/`applyBpGapsGlobal` in
 *     cov_projection.h, applied downstream after node coverage has been expanded into
 *     genomic coordinates.
 */

#ifndef GAM_IO_H
#define GAM_IO_H

#include "cdx_types.h"
#include "coverage_gaps.h"
#include "coverage_precision.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

/**
 * @brief Statistics collected during GAM file parsing and alignment processing.
 *
 * The primary counters track alignment-level mapping metrics, while the compound
 * assignment operator allows safe accumulation of metrics across multithreaded workers.
 */
struct GamMappingStats {
    std::uint64_t total = 0;            // Total number of successfully deserialized alignments.
    std::uint64_t mapped = 0;           // Number of alignments containing at least one valid node ID.
    std::uint64_t mapped_to_query = 0;  // Number of alignments overlapping at least one node in the active query.
    std::uint64_t unmapped = 0;         // Number of alignments containing no valid node IDs.

    /**
     * @brief Accumulates statistics from another instance in-place.
     *
     * @param other The GamMappingStats instance to add.
     * @return Reference to this updated GamMappingStats object.
     */
    GamMappingStats& operator+=(const GamMappingStats& other) noexcept {
        total += other.total;
        mapped += other.mapped;
        mapped_to_query += other.mapped_to_query;
        unmapped += other.unmapped;
        return *this;
    }
};

/**
 * @brief Main parallel coverage engine for processing and analyzing GAM alignments.
 *
 * This function reads and parses a GAM alignment file in parallel using OpenMP tasks and
 * `vg::io::MessageIterator`. It streams alignments in batches, using Protocol Buffers
 * memory arenas for efficient bump allocation. It tracks coverage per node relative to
 * a CDX graph index, restricts calculations to valid query bounds, handles thread-local
 * aggregations, and safely reduces final node coverages while preventing integer overflows.
 *
 * When `precision == CoveragePrecision::Base`, it additionally scans each mapping's `Edit`
 * list (and its `Position::offset`/`is_reverse` fields) to record `BpGap` entries - see
 * coverage_gaps.h - into `out_gaps`. This is strictly additive: the node-level `target`
 * vector is filled identically in both modes. In `CoveragePrecision::Node` mode, `Edit`
 * submessages are never even read, so processing cost and memory stay exactly what they
 * were before base-pair precision existed.
 *
 * @param gam_file Path to the input GAM alignment file.
 * @param target Reference to the output coverage vector where calculated per-node values are stored.
 * @param nid_min Minimum node ID offset defining the start of the local coordinate range.
 * @param batch_size Number of alignments grouped into a single processing task batch.
 * @param decompression_threads Number of dedicated threads used for parallel GAM stream decompression.
 * @param worker_threads Number of active worker threads allocated for parallel alignment parsing and counting.
 * @param precision Coverage computation granularity. Defaults to `CoveragePrecision::Node` so
 *        existing callers that omit this parameter keep their exact original behaviour.
 * @param node_lengths Per-node base-pair length, indexed identically to `target` (i.e.
 *        `node_lengths->at(node_id - nid_min)`). Required (non-null, same size as `target`)
 *        when `precision == CoveragePrecision::Base`; ignored otherwise, so callers running in
 *        `Node` mode don't need to build this lookup table at all.
 * @param out_gaps Output vector of coverage gaps. Cleared and populated when
 *        `precision == CoveragePrecision::Base` (required, must not be null in that case);
 *        left completely untouched in `Node` mode.
 * @return GamMappingStats Aggregated mapping statistics (total, mapped, unmapped, and query-mapped counts).
 *
 * @throws std::invalid_argument If @p batch_size is zero, @p decompression_threads is non-positive,
 *         @p worker_threads is non-positive, @p target is empty, or (in `Base` mode) @p node_lengths /
 *         @p out_gaps are null or @p node_lengths does not match @p target in size.
 * @throws std::runtime_error If the GAM file fails to open or encounters stream reading/parsing errors.
 * @throws std::logic_error If internal statistical invariants are violated.
 */
[[nodiscard]]
GamMappingStats process_gam(
    const std::string& gam_file,
    std::vector<cdx::Coverage>& target,
    cdx::Nid nid_min,
    std::size_t batch_size,
    int decompression_threads,
    int worker_threads,
    CoveragePrecision precision = CoveragePrecision::Node,
    const std::vector<cdx::SeqLen>* node_lengths = nullptr,
    std::vector<BpGap>* out_gaps = nullptr
);

#endif // GAM_IO_H
