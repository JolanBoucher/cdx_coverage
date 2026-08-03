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
    const std::string &gam_file,
    std::vector<cdx::Coverage> &target,
    const cdx::Nid nid_min,
    const std::size_t batch_size,
    const int decompression_threads
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

    // Limite les workers à un maximum de 8 threads.
    const int available_threads = omp_get_max_threads();
    const int active_threads = std::max(1, std::min(available_threads, 8));

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

/**
 * Lit un fichier GAM et affiche les informations du premier alignment trouvé.
 */
bool inspect_first_alignment(const std::string &gam_file) {
    std::ifstream gam_stream(gam_file, std::ios::in | std::ios::binary);
    if (!gam_stream.is_open()) {
        std::cerr << "Cannot open GAM file: " << gam_file << '\n';
        return false;
    }

    try {
        vg::io::MessageIterator it(gam_stream);

        while (it.has_current()) {
            auto tag_and_data = it.take();

            if (vg::io::Registry::check_protobuf_tag<vg::Alignment>(tag_and_data.first)) {
                if (!tag_and_data.second || tag_and_data.second->empty()) {
                    continue;
                }

                vg::Alignment alignment;
                if (!alignment.ParseFromString(*tag_and_data.second)) {
                    std::cerr << "Warning: Could not parse alignment body.\n";
                    continue;
                }

                std::cout << "Read name: " << alignment.name() << '\n';
                std::cout << "Read length: " << alignment.sequence().size() << '\n';
                std::cout << "Mapping quality: " << alignment.mapping_quality() << '\n';
                std::cout << "Alignment score: " << alignment.score() << '\n';
                std::cout << "Secondary: " << (alignment.is_secondary() ? "yes" : "no") << '\n';

                const vg::Path &path = alignment.path();
                std::cout << "Number of mappings: " << path.mapping_size() << "\n\n";

                for (int i = 0; i < path.mapping_size(); ++i) {
                    const vg::Mapping &mapping = path.mapping(i);
                    const vg::Position &position = mapping.position();

                    std::cout << "Mapping " << i << '\n';
                    std::cout << "  node_id: " << position.node_id() << '\n';
                    std::cout << "  orientation: " << (position.is_reverse() ? "reverse" : "forward") << '\n';
                    std::cout << "  node_offset: " << position.offset() << '\n';
                    std::cout << "  rank: " << mapping.rank() << '\n';
                    std::cout << "  edits: " << mapping.edit_size() << '\n';

                    for (int j = 0; j < mapping.edit_size(); ++j) {
                        const vg::Edit &edit = mapping.edit(j);
                        std::cout << "    Edit " << j << '\n';
                        std::cout << "      from_length: " << edit.from_length() << '\n';
                        std::cout << "      to_length: " << edit.to_length() << '\n';
                        if (!edit.sequence().empty()) {
                            std::cout << "      sequence: " << edit.sequence() << '\n';
                        }
                    }
                    std::cout << '\n';
                }

                std::cout << "Node path: ";
                for (int i = 0; i < path.mapping_size(); ++i) {
                    const vg::Position &position = path.mapping(i).position();
                    std::cout << (position.is_reverse() ? '<' : '>') << position.node_id();
                }
                std::cout << '\n';

                return true;
            }
        }
    } catch (const std::exception &error) {
        std::cerr << "Error while reading GAM file: " << error.what() << '\n';
        return false;
    }

    std::cerr << "The GAM file contains no alignments.\n";
    return false;
}
