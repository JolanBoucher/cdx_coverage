#include "gam_io.h"
#include "cdx_types.h"
#include "config.h"
#include "coverage_gaps.h"
#include "coverage_precision.h"

#include <vg/vg.pb.h>
#include <vg/io/edit.hpp>
#include <vg/io/message_iterator.hpp>
#include <vg/io/registry.hpp>
#include <google/protobuf/arena.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

// Processing and analyzing GAM alignments in parallel
GamMappingStats process_gam(
    const std::string& gam_file,
    std::vector<cdx::Coverage>& target,
    cdx::Nid nid_min,
    std::size_t batch_size,
    int decompression_threads,
    int worker_threads,
    CoveragePrecision precision,
    const std::vector<cdx::SeqLen>* node_lengths,
    std::vector<BpGap>* out_gaps
) {
    // 1. Validate mandatory function arguments
    if (batch_size == 0) {
        throw std::invalid_argument("batch_size must be strictly greater than zero.");
    }
    if (decompression_threads <= 0) {
        throw std::invalid_argument("decompression_threads must be strictly greater than zero.");
    }
    if (worker_threads <= 0) {
        throw std::invalid_argument("worker_threads must be strictly greater than zero.");
    }
    if (target.empty()) {
        throw std::invalid_argument("Target coverage vector must not be empty.");
    }

    // Base-precision mode needs a per-node length lookup (to know where a
    // node ends, for deletion/boundary-gap geometry) and somewhere to write
    // the gaps it finds. Both are optional (nullable) parameters precisely
    // so that Node-mode callers don't have to build/own them at all.
    const bool base_precision = (precision == CoveragePrecision::Base);
    if (base_precision) {
        if (node_lengths == nullptr || node_lengths->size() != target.size()) {
            throw std::invalid_argument(
                "node_lengths must be provided and match target's size when "
                "precision == CoveragePrecision::Base."
            );
        }
        if (out_gaps == nullptr) {
            throw std::invalid_argument("out_gaps must be provided when precision == CoveragePrecision::Base.");
        }
        out_gaps->clear();
    }

    std::ifstream gam_stream(gam_file, std::ios::in | std::ios::binary);

    if (!gam_stream) {
        throw std::runtime_error("Failed to open GAM file: " + gam_file);
    }

    const std::size_t coverage_size = target.size();

    // Precompute a fast lookup mask to preserve local nodes belonging strictly to the query region
    std::vector<std::uint8_t> valid_nodes(coverage_size, 0);

    for (std::size_t node_offset = 0; node_offset < coverage_size; ++node_offset) {
        if (target[node_offset] < cfg::NOT_IN_QUERY) {
            valid_nodes[node_offset] = 1;
        }
    }

    const int active_threads = worker_threads;

    // Allocate thread-local structures to prevent false sharing and race conditions during computation
    std::vector local_coverages(static_cast<std::size_t>(active_threads), std::vector<std::uint32_t>(coverage_size, 0));
    std::vector<GamMappingStats> local_stats(static_cast<std::size_t>(active_threads));

    // One gap list per thread, mirroring local_coverages/local_stats above -
    // each thread only ever appends to its own slot, so no synchronization
    // is needed. Left as empty, unused vectors (negligible cost: no heap
    // allocation happens until the first push_back) when base_precision is
    // false, so Node-mode runs pay nothing extra here.
    std::vector<std::vector<BpGap> > local_gaps(static_cast<std::size_t>(active_threads));

    std::exception_ptr reader_exception;

    try {
        // Initialize message iterator with built-in multithreaded GAM decompression
        vg::io::MessageIterator message_it(gam_stream, false, decompression_threads);

#pragma omp parallel num_threads(active_threads)
        {
#pragma omp single
            {
                try {
                    // Main single-threaded loop feeding batches to the OpenMP task pool
                    while (message_it.has_current()) {
                        auto batch = std::make_unique<std::vector<std::string>>();
                        batch->reserve(batch_size);

                        // Group serialized alignments into batches to optimize task granularity
                        while (message_it.has_current() && batch->size() < batch_size) {
                            auto tag_and_data = std::move(message_it.take());

                            if (!vg::io::Registry::check_protobuf_tag<vg::Alignment>(tag_and_data.first)) {
                                continue;
                            }
                            if (tag_and_data.second) {
                                batch->push_back(std::move(*tag_and_data.second));
                            }
                        }

                        if (batch->empty()) continue;
                        auto *raw_batch_ptr = batch.release();

                        // Dispatch each alignment batch as an asynchronous parallel task
#pragma omp task firstprivate(raw_batch_ptr)
                        {
                            std::unique_ptr<std::vector<std::string>> owned_batch(raw_batch_ptr);
                            const auto tid = static_cast<std::size_t>(omp_get_thread_num());

                            std::uint32_t *const coverage = local_coverages[tid].data();

                            GamMappingStats task_stats;

                            // One arena per task/batch: every vg::Alignment
                            // parsed below is bump-allocated from it (its
                            // deeply nested Path -> repeated Mapping ->
                            // Position/Edit submessages included) instead of
                            // going through malloc/free individually, and
                            // everything is released in one bulk sweep when
                            // the arena is destroyed at the end of the task -
                            // no per-message Clear() traversal, no per-field
                            // free(). batch_size bounds how much a single arena
                            // grows to before being freed.
                            google::protobuf::Arena arena;

                            for (const std::string &serialized: *owned_batch) {
                                if (serialized.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                                    continue;
                                }

                                auto *alignment = google::protobuf::Arena::Create<vg::Alignment>(&arena);

                                if (!alignment->ParseFromArray(serialized.data(), static_cast<int>(serialized.size()))) {
                                    continue;
                                }

                                ++task_stats.total;

                                bool read_is_mapped = false;
                                bool read_maps_to_query = false;

                                const vg::Path &path = alignment->path();

                                for (int i = 0; i < path.mapping_size(); ++i) {

                                    const std::int64_t signed_node_id = path.mapping(i).position().node_id();

                                    if (signed_node_id <= 0) continue;

                                    // Mark read as mapped upon encountering at least one positive node ID
                                    read_is_mapped = true;

                                    const auto node_id = static_cast<cdx::Nid>(signed_node_id);

                                    // Skip nodes located upstream of the local range
                                    if (node_id < nid_min) {
                                        continue;
                                    }

                                    const cdx::Nid raw_offset = node_id - nid_min;

                                    // Skip nodes located downstream of the local range
                                    if (raw_offset >= static_cast<cdx::Nid>(coverage_size)) {
                                        continue;
                                    }

                                    const std::size_t node_offset = static_cast<std::size_t>(raw_offset);

                                    // Skip nodes within the component but excluded from the specific query scope
                                    if (!valid_nodes[node_offset]) {
                                        continue;
                                    }

                                    ++coverage[node_offset];
                                    read_maps_to_query = true;

                                    // ---- Base-precision gap detection ----
                                    // The increment above already gave this
                                    // node a full +1 node-level credit for
                                    // this read. Everything below is a pure
                                    // refinement of that credit: it locates
                                    // which base-pair ranges of the node
                                    // this specific mapping did NOT actually
                                    // traverse (deletions, and - only for
                                    // the first/last mapping of the path -
                                    // read-boundary under-reach), and
                                    // records them as BpGap entries for
                                    // later subtraction once coverage has
                                    // been expanded into base-pair space
                                    // (see applyBpGapsQuery/Global in
                                    // cov_projection.h). It never modifies
                                    // `coverage` itself.
                                    if (base_precision) {
                                        const cdx::SeqLen node_length = (*node_lengths)[node_offset];

                                        // A zero-length node lookup can only happen on malformed/
                                        // inconsistent input (every real node has >= 1 base pair) -
                                        // skip rather than feed a degenerate length into the geometry
                                        // helpers below.
                                        if (node_length > 0) {
                                            const vg::Mapping &mapping = path.mapping(i);
                                            const vg::Position &mapping_position = mapping.position();

                                            // Clamp defensively: a corrupt/adversarial GAM could in
                                            // principle claim an offset past the node's own length.
                                            // Named distinctly from the outer `raw_offset` (the node's
                                            // nid-offset, used below when emitting BpGap entries) to
                                            // avoid shadowing it - this one is the mapping's intra-node
                                            // Position::offset() instead.
                                            const std::int64_t raw_mapping_offset = mapping_position.offset();
                                            const cdx::SeqLen offset = raw_mapping_offset > 0
                                                ? static_cast<cdx::SeqLen>(std::min<std::int64_t>(raw_mapping_offset, node_length))
                                                : cdx::SeqLen{0};
                                            const bool is_reverse = mapping_position.is_reverse();

                                            const bool is_first_mapping = (i == 0);
                                            const bool is_last_mapping = (i == path.mapping_size() - 1);

                                            // Bases consumed so far, expressed in this mapping's own
                                            // walk order (0 = the first base the mapping consumes,
                                            // regardless of strand) - fed to walkSpanToForwardRange /
                                            // trailingUncoveredRange, which handle the forward/reverse
                                            // orientation conversion.
                                            cdx::SeqLen walk_cursor = 0;

                                            for (int edit_idx = 0; edit_idx < mapping.edit_size(); ++edit_idx) {
                                                const vg::Edit &edit = mapping.edit(edit_idx);
                                                const std::int32_t raw_from_length = edit.from_length();
                                                const cdx::SeqLen edit_from_length = raw_from_length > 0
                                                    ? static_cast<cdx::SeqLen>(raw_from_length)
                                                    : cdx::SeqLen{0};

                                                // Only deletions leave a base-pair gap: they consume
                                                // "from" (reference/node) length without the read ever
                                                // being sequenced against it. Matches AND mismatches/
                                                // substitutions both consume from_length while genuinely
                                                // being covered by a read base, so they are deliberately
                                                // left uncorrected - this tool reports raw read depth,
                                                // not concordance with the reference. Insertions consume
                                                // no "from" length at all, so they can never produce a
                                                // gap either (there is no node position to subtract from).
                                                if (vg::io::edit_is_deletion(edit)) {
                                                    const ForwardRange gap_range = walkSpanToForwardRange(
                                                        node_length, offset, is_reverse,
                                                        walk_cursor, walk_cursor + edit_from_length
                                                    );
                                                    if (!gap_range.empty()) {
                                                        local_gaps[tid].push_back(BpGap{raw_offset, gap_range});
                                                    }
                                                }

                                                walk_cursor += edit_from_length;
                                            }

                                            // Boundary under-coverage: only the very first and last
                                            // mapping of a read's path can start/end partway through
                                            // their node - internal mappings always walk their node
                                            // edge to edge, so these are no-ops for them by construction.
                                            if (is_first_mapping) {
                                                if (const auto gap = leadingUncoveredRange(node_length, offset, is_reverse)) {
                                                    local_gaps[tid].push_back(BpGap{raw_offset, *gap});
                                                }
                                            }
                                            if (is_last_mapping) {
                                                if (const auto gap = trailingUncoveredRange(node_length, offset, is_reverse, walk_cursor)) {
                                                    local_gaps[tid].push_back(BpGap{raw_offset, *gap});
                                                }
                                            }
                                        }
                                    }
                                }

                                if (read_is_mapped) {
                                    ++task_stats.mapped;
                                } else {
                                    ++task_stats.unmapped;
                                }

                                if (read_maps_to_query) {
                                    ++task_stats.mapped_to_query;
                                }
                            }
                            local_stats[tid] += task_stats;
                        }
                    }
                } catch (...) {
                    reader_exception = std::current_exception();
                }

#pragma omp taskwait
            }
        }
    } catch (...) {
        reader_exception = std::current_exception();
    }

    // Propagate any exceptions caught during asynchronous stream processing
    if (reader_exception) {
        try {
            std::rethrow_exception(reader_exception);
        } catch (const std::exception &error) {
            throw std::runtime_error("Failed to parse GAM file '" + gam_file + "': " + error.what());
        }
    }

    // Aggregate thread-local statistics into global metrics
    GamMappingStats stats;
    for (const GamMappingStats &thread_stats: local_stats) {
        stats += thread_stats;
    }

    // Validate fundamental statistical invariants
    if (stats.total != stats.mapped + stats.unmapped) {
        throw std::logic_error("Inconsistent GAM mapping statistics: total != mapped + unmapped.");
    }

    if (stats.mapped_to_query > stats.mapped) {
        throw std::logic_error("Inconsistent GAM mapping statistics: mapped_to_query exceeds mapped.");
    }

#ifndef NDEBUG
    std::cerr
            << "[DEBUG] GAM alignments total: "
            << stats.total
            << '\n'
            << "[DEBUG] GAM alignments mapped: "
            << stats.mapped
            << '\n'
            << "[DEBUG] GAM alignments mapped to query: "
            << stats.mapped_to_query
            << '\n'
            << "[DEBUG] GAM alignments unmapped: "
            << stats.unmapped
            << '\n';
#endif

    // Parallel reduction of thread-local coverage vectors into the final target vector
#pragma omp parallel num_threads(active_threads)
    {
        const std::size_t worker = static_cast<std::size_t>(omp_get_thread_num());
        const std::size_t workers = static_cast<std::size_t>(omp_get_num_threads());
        const std::size_t block = (coverage_size + workers - 1) / workers;
        const std::size_t begin = worker * block;
        const std::size_t end = std::min(begin + block, coverage_size);

        for (std::size_t node_offset = begin;
             node_offset < end;
             ++node_offset) {
            if (!valid_nodes[node_offset]) continue;

            std::uint64_t total_coverage = 0;
            for (const auto &source: local_coverages) {
                total_coverage += source[node_offset];
            }

            // Clamp coverage values to prevent integer overflow past cdx::Coverage capacity
            if (total_coverage > static_cast<std::uint64_t>(std::numeric_limits<cdx::Coverage>::max())) {
                target[node_offset] = std::numeric_limits<cdx::Coverage>::max();
            } else {
                target[node_offset] = static_cast<cdx::Coverage>(total_coverage);
            }
        }
    }

    // Concatenate every thread's gap list into the caller-supplied output
    // vector. This is a single-threaded linear pass, but it operates on
    // however many gaps were actually found (proportional to read/edit
    // count), not on graph or genome size - the same cost class as
    // `expandPosCovQuery`'s single-threaded fill, not the parallel scan
    // above, so it stays cheap even for very large runs.
    if (base_precision) {
        std::size_t total_gap_count = 0;
        for (const auto &thread_gap_list: local_gaps) {
            total_gap_count += thread_gap_list.size();
        }

        out_gaps->reserve(out_gaps->size() + total_gap_count);
        for (auto &thread_gap_list: local_gaps) {
            out_gaps->insert(
                out_gaps->end(),
                std::make_move_iterator(thread_gap_list.begin()),
                std::make_move_iterator(thread_gap_list.end())
            );
        }
    }

    return stats;
}