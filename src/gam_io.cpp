#include "gam_io.h"
#include "cdx_types.h"
#include "config.h"

#include <vg/vg.pb.h>
#include <vg/io/message_iterator.hpp>
#include <vg/io/registry.hpp>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

GamMappingStats process_gam(
    const std::string& gam_file,
    std::vector<cdx::Coverage>& target,
    cdx::Nid nid_min,
    std::size_t batch_size,
    int decompression_threads,
    int worker_threads
) {
    // 1. Validation des arguments
    if (batch_size == 0) {
        throw std::invalid_argument(
            "batch_size must be strictly greater than zero."
        );
    }

    if (decompression_threads <= 0) {
        throw std::invalid_argument(
            "decompression_threads must be strictly greater than zero."
        );
    }

    if (target.empty()) {
        throw std::invalid_argument(
            "Target coverage vector must not be empty."
        );
    }

    std::ifstream gam_stream(gam_file, std::ios::in | std::ios::binary);

    if (!gam_stream) {
        throw std::runtime_error("Failed to open GAM file: " + gam_file);
    }

    const std::size_t coverage_size = target.size();

    /*
     * Préserve les nœuds locaux appartenant réellement à la query.
     */
    std::vector<std::uint8_t> valid_nodes(coverage_size, 0);

    for (std::size_t node_offset = 0; node_offset < coverage_size; ++node_offset) {
        if (target[node_offset] < cfg::NOT_IN_QUERY) {
            valid_nodes[node_offset] = 1;
        }
    }

    // ajoute les worker threads passer en args
    const int active_threads = worker_threads;

    std::vector local_coverages(static_cast<std::size_t>(active_threads), std::vector<std::uint32_t>(coverage_size, 0));
    std::vector<GamMappingStats> local_stats(static_cast<std::size_t>(active_threads));
    std::exception_ptr reader_exception;

    try {
        vg::io::MessageIterator message_it(gam_stream, false, decompression_threads);
#pragma omp parallel num_threads(active_threads)
        {
#pragma omp single
            {
                try {
                    while (message_it.has_current()) {
                        auto batch = std::make_unique<std::vector<std::string> >();
                        batch->reserve(batch_size);

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

#pragma omp task firstprivate(raw_batch_ptr)
                        {
                            std::unique_ptr<std::vector<std::string> > owned_batch(raw_batch_ptr);
                            const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());

                            std::uint32_t *const coverage = local_coverages[tid].data();

                            GamMappingStats task_stats;
                            vg::Alignment alignment;

                            for (const std::string &serialized: *owned_batch) {
                                if (serialized.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                                    continue;
                                }
                                alignment.Clear();

                                if (!alignment.ParseFromArray(serialized.data(), static_cast<int>(serialized.size()))) {
                                    continue;
                                }

                                ++task_stats.total;

                                bool read_is_mapped = false;
                                bool read_maps_to_query = false;

                                const vg::Path &path = alignment.path();

                                for (int i = 0; i < path.mapping_size(); ++i) {

                                    const std::int64_t signed_node_id = path.mapping(i).position().node_id();

                                    if (signed_node_id <= 0) continue;

                                    // Le read est considéré comme mapped dès qu'il contient au moins un node ID positif
                                    read_is_mapped = true;

                                    const auto node_id = static_cast<cdx::Nid>(signed_node_id);

                                    // Nœud situé avant la plage locale.
                                    if (node_id < nid_min) {
                                        continue;
                                    }

                                    const cdx::Nid raw_offset = node_id - nid_min;

                                    // Nœud situé après la plage locale.
                                    if (raw_offset >= static_cast<cdx::Nid>(coverage_size)) {
                                        continue;
                                    }

                                    const std::size_t node_offset = static_cast<std::size_t>(raw_offset);

                                    // Nœud dans la composante, mais exclu de la query.
                                    if (!valid_nodes[node_offset]) {
                                        continue;
                                    }

                                    ++coverage[node_offset];
                                    read_maps_to_query = true;
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

    if (reader_exception) {
        try {
            std::rethrow_exception(reader_exception);
        } catch (const std::exception &error) {
            throw std::runtime_error("Failed to parse GAM file '" + gam_file + "': " + error.what());
        }
    }

    // Agrégation des statistiques thread-local
    GamMappingStats stats;
    for (const GamMappingStats &thread_stats: local_stats) {
        stats += thread_stats;
    }

    // Vérification des invariants de statistiques
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

    /*
     * Réduction du vecteur de couverture local vers le vecteur target final.
     */
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

            if (total_coverage > static_cast<std::uint64_t>(std::numeric_limits<cdx::Coverage>::max())) {
                target[node_offset] = std::numeric_limits<cdx::Coverage>::max();
            } else {
                target[node_offset] = static_cast<cdx::Coverage>(total_coverage);
            }
        }
    }

    return stats;
}