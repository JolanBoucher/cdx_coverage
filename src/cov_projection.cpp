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

        // Filter out nodes not included in query/component space
        // 3. Fix: Use correct index sentinel (cfg::INVALID_IDX) rather than coverage sentinel
        if (local_idx == cfg::NOT_IN_QUERY) {
            continue;
        }

#ifndef NDEBUG
        // Sanity check: Ensure local index does not write outside component target vector
        if (static_cast<std::size_t>(local_idx) >= component_size) {
            throw std::out_of_range("Local index exceeds component size bounds.");
        }
#endif

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

#ifndef NDEBUG
        // Assert prefix sum monotonicity
        if (start_bp > end_bp) {
            throw std::runtime_error("Non-monotonic base-pair coordinates detected in idx2bp.");
        }
#endif

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

#ifndef NDEBUG
        if (static_cast<std::size_t>(flat_idx) >= flat_idx_cov_table.size()) {
            throw std::out_of_range("Flat index exceeds global coverage table bounds.");
        }
#endif

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

#ifndef NDEBUG
            if (local_pos_idx + 1 >= idx2bp.size()) {
                throw std::out_of_range("idx2bp access exceeds bounds.");
            }
#endif

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