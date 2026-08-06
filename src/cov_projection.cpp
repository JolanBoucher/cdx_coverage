/**
 * @file cov_projection.cpp
 * @brief Coverage projection and base-pair expansion routines.
 */

#include "cov_projection.h"
#include "config.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
    /**
     * @brief Decrements every position in [start, end) of @p bp_coverage by 1,
     *        throwing rather than silently underflowing if any position is
     *        already at 0.
     *
     * `cdx::Coverage` is unsigned, so an unchecked `--bp_coverage[pos]` on a
     * value already at 0 would wrap around to the type's maximum value
     * instead of going negative - silently turning a correction into a
     * spurious coverage spike. That should never legitimately happen: every
     * BpGap is produced by exactly one read's mapping, and that same read
     * also contributed the +1 node-level credit being refined (see
     * gam_io.cpp's per-mapping gap-detection block), so the total number of
     * gaps ever touching a given position can never exceed that position's
     * coverage value. Hitting 0 here therefore indicates a real bug in the
     * gap geometry (coverage_gaps.cpp) or in how it was wired up - not
     * malformed GAM content - hence throwing instead of clamping.
     *
     * @throws std::out_of_range If @p end exceeds bp_coverage's size.
     * @throws std::logic_error If a position within the range is already 0.
     */
    void decrementRangeOrThrow(
        std::vector<cdx::Coverage> &bp_coverage,
        const cdx::PosBp start,
        const cdx::PosBp end
    ) {
        if (end > bp_coverage.size()) {
            throw std::out_of_range(
                "Coverage gap range [" + std::to_string(start) + ", " + std::to_string(end) +
                ") exceeds bp_coverage bounds (" + std::to_string(bp_coverage.size()) + ")."
            );
        }

        for (cdx::PosBp pos = start; pos < end; ++pos) {
            if (bp_coverage[pos] == 0) {
                throw std::logic_error(
                    "Coverage gap would decrement bp_coverage[" + std::to_string(pos) +
                    "] below zero - this indicates a bug in coverage-gap geometry rather than "
                    "malformed input, since every gap is bounded by the read-derived node-level "
                    "credit it refines."
                );
            }
            --bp_coverage[pos];
        }
    }

    /**
     * @brief Locates the component owning a global flattened node index.
     *
     * @p component_offsets partitions the flat index space [0, total_nodes)
     * into contiguous per-component ranges; the component owning @p flat_idx
     * is the last boundary that is <= flat_idx.
     *
     * @throws std::out_of_range If @p flat_idx does not fall within any
     *         component's range described by @p component_offsets.
     */
    [[nodiscard]] std::size_t findComponentForFlatIdx(
        const std::vector<cdx::RecordCount> &component_offsets,
        const cdx::FlatIdx flat_idx
    ) {
        // upper_bound finds the first boundary strictly greater than
        // flat_idx; the owning component is the one just before it.
        const auto it = std::upper_bound(
            component_offsets.begin(), component_offsets.end(),
            static_cast<cdx::RecordCount>(flat_idx)
        );

        if (it == component_offsets.begin() || it == component_offsets.end()) {
            throw std::out_of_range("Flat index does not fall within any component's node range.");
        }

        return static_cast<std::size_t>(std::distance(component_offsets.begin(), it)) - 1;
    }
} // anonymous namespace

// Projects query coverage from nid-space into local idx-space.
[[nodiscard]]
std::vector<cdx::Coverage> projectCov2IdxQuery(
    const std::vector<cdx::Coverage>& cov_table,
    const std::vector<cdx::Idx>& nid2idx,
    const std::size_t component_size
) {
    // 1. Structural check: verify 1:1 mapping length between coverage table and index map
    if (cov_table.size() != nid2idx.size()) {
        throw std::invalid_argument(
            "Dimension mismatch: cov_table size (" + std::to_string(cov_table.size()) +
            ") must match nid2idx size (" + std::to_string(nid2idx.size()) + ")."
        );
    }

    // 2. Allocate zero-initialized destination array sized to local component node count
    std::vector<cdx::Coverage> idx_cov_table(component_size, 0);

    // 3. Iterate over input offsets and remap coverage entries
    for (std::size_t nid_offset = 0; nid_offset < cov_table.size(); ++nid_offset) {
        const cdx::Idx local_idx = nid2idx[nid_offset];

        // Filter out nodes not included in query/component space. Both
        // sentinels are checked explicitly (rather than relying on
        // NOT_IN_COMPO incidentally overflowing the bounds check below)
        // so the skip condition documents its own intent.
        if (local_idx == cfg::NOT_IN_QUERY || local_idx == cfg::NOT_IN_COMPO) {
            continue;
        }

        // Always validated (not just in debug builds): an out-of-range
        // local_idx here indicates malformed input data, not an internal
        // invariant violation, and callers rely on this being caught.
        if (static_cast<std::size_t>(local_idx) >= component_size) {
            throw std::out_of_range("Local index exceeds component size bounds.");
        }

        // Copy coverage value into the local topological slot
        const cdx::Coverage coverage = cov_table[nid_offset];
        if (coverage != 0) {
            idx_cov_table[static_cast<std::size_t>(local_idx)] = coverage;
        }
    }

    return idx_cov_table;
}

// Expands node-level query coverage into bp-space.
[[nodiscard]]
std::vector<cdx::Coverage> expandPosCovQuery(
    const std::vector<cdx::Coverage>& idx_cov_table,
    const std::vector<cdx::PosBp>& idx2bp
) {
    // 1. Ensure prefix-sum array contains at least origin (0) and total length
    if (idx2bp.size() < 2) {
        throw std::invalid_argument("idx2bp array must contain at least origin and final position.");
    }

    // 2. Validate prefix-sum sizing invariant: len(idx2bp) == num_nodes + 1
    const std::size_t node_count = idx_cov_table.size();
    if (idx2bp.size() != node_count + 1) {
        throw std::invalid_argument(
            "Dimension mismatch: idx2bp size (" + std::to_string(idx2bp.size()) +
            ") must equal idx_cov_table size + 1 (" + std::to_string(node_count + 1) + ")."
        );
    }

    // Total base pairs in component is recorded at the last slot of the prefix-sum vector
    const cdx::PosBp component_length_bp = idx2bp.back();

    // 3. Allocate per-base-pair array initialized with sentinel default values
    // 4. Fix: Explicit vector type instantiation
    std::vector bp_cov_table(component_length_bp, cfg::INVALID_NODE);

    // 4. Expand coverage depth across each node's base-pair range
    for (std::size_t local_idx = 0; local_idx < node_count; ++local_idx) {
        const cdx::PosBp start_bp = idx2bp[local_idx];
        const cdx::PosBp end_bp = idx2bp[local_idx + 1];

        // Always validated: a non-monotonic idx2bp reflects malformed
        // input data, not just an internal debug invariant.
        if (start_bp > end_bp) {
            throw std::runtime_error("Non-monotonic base-pair coordinates detected in idx2bp.");
        }

        const cdx::Coverage coverage = idx_cov_table[local_idx];

        // Fill memory range [start_bp, end_bp) with node depth
        std::fill(
            bp_cov_table.begin() + static_cast<std::ptrdiff_t>(start_bp),
            bp_cov_table.begin() + static_cast<std::ptrdiff_t>(end_bp),
            coverage
        );
    }

    return bp_cov_table;
}

// Projects global coverage from nid-space into flat idx-space.
[[nodiscard]]
std::vector<cdx::Coverage> projectCov2IdxGlobal(
    const std::vector<cdx::Coverage>& cov_table,
    const std::vector<cdx::FlatIdx>& nid2flat_idx,
    const std::vector<cdx::RecordCount>& component_offsets
) {
    // 1. Verify input vector dimensions match
    if (cov_table.size() != nid2flat_idx.size()) {
        throw std::invalid_argument(
            "Dimension mismatch: cov_table size (" + std::to_string(cov_table.size()) +
            ") must match nid2flat_idx size (" + std::to_string(nid2flat_idx.size()) + ")."
        );
    }

    if (component_offsets.size() < 2) {
        throw std::invalid_argument("component_offsets must contain at least two boundaries.");
    }

    // 2. Allocate flat global table sized to total node count (component_offsets.back())
    std::vector<cdx::Coverage> flat_idx_cov_table( (component_offsets.back()),0);

    // 3. Remap non-zero coverages to flat indices
    for (std::size_t node_offset = 0; node_offset < cov_table.size(); ++node_offset) {
        const cdx::FlatIdx flat_idx = nid2flat_idx[node_offset];

        if (flat_idx == cfg::INVALID_FLAT_IDX) {
            continue;
        }

        // Always validated: an out-of-range flat_idx indicates malformed
        // input data, not just an internal debug invariant.
        if (static_cast<std::size_t>(flat_idx) >= flat_idx_cov_table.size()) {
            throw std::out_of_range("Flat index exceeds global coverage table bounds.");
        }

        const cdx::Coverage coverage = cov_table[node_offset];
        if (coverage != 0) {
            flat_idx_cov_table[static_cast<std::size_t>(flat_idx)] = coverage;
        }
    }

    return flat_idx_cov_table;
}

// Expands flattened global node coverage into flattened bp-space.
[[nodiscard]]
std::pair<std::vector<cdx::Coverage>, std::vector<cdx::PosBp>> expandPosCovGlobal(
    const std::vector<cdx::Coverage>& flat_idx_cov_table,
    const std::vector<cdx::RecordCount>& component_offsets,
    const std::vector<cdx::PosBp>& idx2bp,
    const std::vector<cdx::RecordCount>& idx2bp_offsets
) {
    // 1. Boundary array validation
    if (component_offsets.size() < 2 || idx2bp_offsets.size() != component_offsets.size()) {
        throw std::invalid_argument(
            "Component boundary offset tables must contain identical sizes of at least 2 entries.");
    }
    if (component_offsets.front() != 0 || idx2bp_offsets.front() != 0) {
        throw std::invalid_argument("The first component and idx2bp offsets must be zero.");
    }
    if (component_offsets.back() != flat_idx_cov_table.size()) {
        throw std::invalid_argument(
            "Final component offset (" + std::to_string(component_offsets.back()) +") "
            "must match flat_idx_cov_table length (" + std::to_string(flat_idx_cov_table.size()) + ").");
    }
    if (idx2bp_offsets.back() != idx2bp.size()) {
        throw std::invalid_argument(
            "Final idx2bp offset (" + std::to_string(idx2bp_offsets.back()) +
            ") must match total idx2bp array size (" + std::to_string(idx2bp.size()) + ").");
    }

    const std::size_t component_count = component_offsets.size() - 1;

    // 2. Pre-pass: Compute total graph base-pair length across all components
    cdx::PosBp running_bp = 0;
    for (std::size_t comp_id = 0; comp_id < component_count; ++comp_id) {

        const std::size_t pos_end = idx2bp_offsets[comp_id + 1];

        if (pos_end > 0 && pos_end <= idx2bp.size()) {
            running_bp += idx2bp[pos_end - 1];
        }
    }

    // Guard against allocations exceeding machine size_t limits
    if (running_bp > static_cast<cdx::PosBp>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("BP coverage table size exceeds system size_t memory capacity.");
    }

    // 3. Allocate combined global base-pair coverage table
    std::vector flat_bp_cov_table(running_bp, cfg::INVALID_NODE);
    std::vector<cdx::PosBp> component_bp_offsets;

    component_bp_offsets.reserve(component_count + 1);
    component_bp_offsets.push_back(0);

    cdx::PosBp global_bp_offset = 0;

    // 4. Expand nodes component by component into contiguous base-pair output
    for (std::size_t comp_id = 0; comp_id < component_count; ++comp_id) {

        const cdx::RecordCount node_start = component_offsets[comp_id];
        const cdx::RecordCount node_end = component_offsets[comp_id + 1];
        const std::size_t node_count = node_end - node_start;
        const std::size_t pos_offset = idx2bp_offsets[comp_id];

        for (std::size_t local_idx = 0; local_idx < node_count; ++local_idx) {
            const std::size_t local_pos_idx = pos_offset + local_idx;

            // Always validated: inconsistent component/idx2bp offset
            // tables indicate malformed input data, not just an internal
            // debug invariant.
            if (local_pos_idx + 1 >= idx2bp.size()) {
                throw std::out_of_range("idx2bp access exceeds bounds.");
            }

            // Local node base-pair coordinates within component
            const cdx::PosBp local_start_bp = idx2bp[local_pos_idx];
            const cdx::PosBp local_end_bp = idx2bp[local_pos_idx + 1];

            // Map local offsets to global concatenated base-pair coordinates
            const cdx::PosBp flat_start_bp = global_bp_offset + local_start_bp;
            const cdx::PosBp flat_end_bp = global_bp_offset + local_end_bp;

            const cdx::Coverage coverage = flat_idx_cov_table[static_cast<std::size_t>(node_start) + local_idx];

            std::fill(
                flat_bp_cov_table.begin() + static_cast<std::ptrdiff_t>(flat_start_bp),
                flat_bp_cov_table.begin() + static_cast<std::ptrdiff_t>(flat_end_bp),
                coverage
            );
        }

        // Store component global base-pair boundary offset
        const std::size_t comp_last_pos_idx = idx2bp_offsets[comp_id + 1] - 1;
        const cdx::PosBp comp_length_bp = idx2bp[comp_last_pos_idx];
        global_bp_offset += comp_length_bp;

        component_bp_offsets.push_back(global_bp_offset);
    }

    return {std::move(flat_bp_cov_table), std::move(component_bp_offsets)};
}

// Masks positions outside a query interval.
[[nodiscard]]
std::vector<cdx::Coverage> trimCoverageToQuery(
    const std::vector<cdx::Coverage>& bp_cov_table,
    const std::pair<cdx::PosBp, cdx::PosBp>& query_bound
){
    // Early exit on empty input array
    if (bp_cov_table.empty()) return bp_cov_table;

    std::vector<cdx::Coverage> trimmed = bp_cov_table;
    const auto component_length = trimmed.size();
    const auto [query_start, query_end] = query_bound;

    // Validate query boundaries stay within array dimensions
    if (query_start >= component_length) {
        throw std::out_of_range("query_start outside component bounds.");
    }
    if (query_end >= component_length) {
        throw std::out_of_range("query_end outside component bounds.");
    }

    if (query_start <= query_end){
        // Standard range: Mask prefix [0, query_start)
        for (std::size_t pos = 0; pos < query_start; ++pos){
            if (trimmed[pos] < cfg::NOT_IN_QUERY) {
                trimmed[pos] = cfg::NOT_IN_QUERY;
            }
        }

        // Standard range: Mask suffix (query_end, component_length)
        for (std::size_t pos = query_end + 1; pos < component_length; ++pos){
            if (trimmed[pos] < cfg::NOT_IN_QUERY) {
                trimmed[pos] = cfg::NOT_IN_QUERY;
            }
        }
    }
    else {
        // Inverted range: Mask region between query_end and query_start
        for (std::size_t pos = query_end + 1; pos < query_start; ++pos){
            if (trimmed[pos] < cfg::NOT_IN_QUERY) {
                trimmed[pos] = cfg::NOT_IN_QUERY;
            }
        }
    }
    return trimmed;
}

// Builds a raw-nid-offset-indexed node length lookup for a single component.
[[nodiscard]]
std::vector<cdx::SeqLen> buildNodeLengthsQuery(
    const std::vector<cdx::Idx> &nid2idx,
    const std::vector<cdx::PosBp> &idx2bp
) {
    if (idx2bp.size() < 2) {
        throw std::invalid_argument("idx2bp array must contain at least origin and final position.");
    }

    const std::size_t component_size = idx2bp.size() - 1;
    std::vector<cdx::SeqLen> node_lengths(nid2idx.size(), 0);

    for (std::size_t nid_offset = 0; nid_offset < nid2idx.size(); ++nid_offset) {
        const cdx::Idx local_idx = nid2idx[nid_offset];

        // Sentinel-mapped nodes are outside the active query/component -
        // process_gam filters them out via valid_nodes before gap detection
        // would ever consult their length, so a length of 0 is never read.
        if (local_idx == cfg::NOT_IN_QUERY || local_idx == cfg::NOT_IN_COMPO) {
            continue;
        }

        if (static_cast<std::size_t>(local_idx) >= component_size) {
            throw std::out_of_range("Local index exceeds component size bounds.");
        }

        const cdx::PosBp start_bp = idx2bp[local_idx];
        const cdx::PosBp end_bp = idx2bp[static_cast<std::size_t>(local_idx) + 1];

        if (start_bp > end_bp) {
            throw std::runtime_error("Non-monotonic base-pair coordinates detected in idx2bp.");
        }

        const cdx::PosBp length = end_bp - start_bp;

        // Defensive clamp against SeqLen (32-bit) overflow: no real node
        // approaches 4 billion base pairs, but this keeps the cast total.
        node_lengths[nid_offset] = length > std::numeric_limits<cdx::SeqLen>::max()
                ? std::numeric_limits<cdx::SeqLen>::max()
                : static_cast<cdx::SeqLen>(length);
    }

    return node_lengths;
}

// Builds a raw-nid-offset-indexed node length lookup for the whole graph.
[[nodiscard]]
std::vector<cdx::SeqLen> buildNodeLengthsGlobal(
    const std::vector<cdx::FlatIdx> &nid2flat_idx,
    const std::vector<cdx::RecordCount> &component_offsets,
    const std::vector<cdx::PosBp> &idx2bp,
    const std::vector<cdx::RecordCount> &idx2bp_offsets
) {
    if (component_offsets.size() < 2 || idx2bp_offsets.size() != component_offsets.size()) {
        throw std::invalid_argument(
            "Component boundary offset tables must contain identical sizes of at least 2 entries.");
    }

    std::vector<cdx::SeqLen> node_lengths(nid2flat_idx.size(), 0);

    for (std::size_t nid_offset = 0; nid_offset < nid2flat_idx.size(); ++nid_offset) {
        const cdx::FlatIdx flat_idx = nid2flat_idx[nid_offset];

        if (flat_idx == cfg::INVALID_FLAT_IDX) {
            continue;
        }

        const std::size_t comp_id = findComponentForFlatIdx(component_offsets, flat_idx);
        const cdx::FlatIdx local_idx_in_component = flat_idx - static_cast<cdx::FlatIdx>(component_offsets[comp_id]);
        const std::size_t pos_idx = static_cast<std::size_t>(idx2bp_offsets[comp_id]) +
                                     static_cast<std::size_t>(local_idx_in_component);

        if (pos_idx + 1 >= idx2bp.size()) {
            throw std::out_of_range("idx2bp access exceeds bounds while building global node lengths.");
        }

        const cdx::PosBp start_bp = idx2bp[pos_idx];
        const cdx::PosBp end_bp = idx2bp[pos_idx + 1];

        if (start_bp > end_bp) {
            throw std::runtime_error("Non-monotonic base-pair coordinates detected in idx2bp.");
        }

        const cdx::PosBp length = end_bp - start_bp;
        node_lengths[nid_offset] = length > std::numeric_limits<cdx::SeqLen>::max()
                ? std::numeric_limits<cdx::SeqLen>::max()
                : static_cast<cdx::SeqLen>(length);
    }

    return node_lengths;
}

// Applies base-pair-precision coverage gaps to a single-component array.
void applyBpGapsQuery(
    std::vector<cdx::Coverage> &bp_coverage,
    const std::vector<BpGap> &gaps,
    const std::vector<cdx::Idx> &nid2idx,
    const std::vector<cdx::PosBp> &idx2bp
) {
    for (const BpGap &gap: gaps) {
        if (gap.nid_offset >= nid2idx.size()) {
            throw std::out_of_range("BpGap nid_offset exceeds nid2idx bounds.");
        }

        const cdx::Idx local_idx = nid2idx[gap.nid_offset];

        // The gap's node was excluded from the active query/component by
        // process_gam itself (same sentinel check as projectCov2IdxQuery) -
        // there is nothing in bp_coverage to correct for it.
        if (local_idx == cfg::NOT_IN_QUERY || local_idx == cfg::NOT_IN_COMPO) {
            continue;
        }

        if (static_cast<std::size_t>(local_idx) + 1 >= idx2bp.size()) {
            throw std::out_of_range("Local index exceeds idx2bp bounds while applying coverage gaps.");
        }

        const cdx::PosBp node_start_bp = idx2bp[local_idx];
        const cdx::PosBp global_start = node_start_bp + gap.range.start;
        const cdx::PosBp global_end = node_start_bp + gap.range.end;

        decrementRangeOrThrow(bp_coverage, global_start, global_end);
    }
}

// Applies base-pair-precision coverage gaps to the flattened global array.
void applyBpGapsGlobal(
    std::vector<cdx::Coverage> &flat_bp_coverage,
    const std::vector<BpGap> &gaps,
    const std::vector<cdx::FlatIdx> &nid2flat_idx,
    const std::vector<cdx::RecordCount> &component_offsets,
    const std::vector<cdx::PosBp> &idx2bp,
    const std::vector<cdx::RecordCount> &idx2bp_offsets,
    const std::vector<cdx::PosBp> &component_bp_offsets
) {
    if (component_offsets.size() < 2
        || idx2bp_offsets.size() != component_offsets.size()
        || component_bp_offsets.size() != component_offsets.size()) {
        throw std::invalid_argument(
            "Component boundary offset tables must contain identical sizes of at least 2 entries.");
    }

    for (const BpGap &gap: gaps) {
        if (gap.nid_offset >= nid2flat_idx.size()) {
            throw std::out_of_range("BpGap nid_offset exceeds nid2flat_idx bounds.");
        }

        const cdx::FlatIdx flat_idx = nid2flat_idx[gap.nid_offset];

        // Same reasoning as applyBpGapsQuery(): a node with no flat index
        // was excluded from GlobalData::node_coverage entirely.
        if (flat_idx == cfg::INVALID_FLAT_IDX) {
            continue;
        }

        const std::size_t comp_id = findComponentForFlatIdx(component_offsets, flat_idx);
        const cdx::FlatIdx local_idx_in_component = flat_idx - static_cast<cdx::FlatIdx>(component_offsets[comp_id]);
        const std::size_t pos_idx = static_cast<std::size_t>(idx2bp_offsets[comp_id]) +
                                     static_cast<std::size_t>(local_idx_in_component);

        if (pos_idx + 1 >= idx2bp.size()) {
            throw std::out_of_range("idx2bp access exceeds bounds while applying coverage gaps.");
        }

        const cdx::PosBp node_start_bp = idx2bp[pos_idx];
        const cdx::PosBp component_offset_bp = component_bp_offsets[comp_id];

        const cdx::PosBp global_start = component_offset_bp + node_start_bp + gap.range.start;
        const cdx::PosBp global_end = component_offset_bp + node_start_bp + gap.range.end;

        decrementRangeOrThrow(flat_bp_coverage, global_start, global_end);
    }
}