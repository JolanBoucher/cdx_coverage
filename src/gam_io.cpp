#include "gam_io.h"

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

#include <vg/vg.pb.h>
#include <vg/io/message_iterator.hpp>
#include <vg/io/registry.hpp>

void process_gam_fast(
    const std::string& gam_file,
    std::vector<uint32_t>& global_coverage,
    uint64_t& read_count,
    std::size_t max_node_id,
    std::size_t batch_size,
    int decompression_threads
) {
    // 1. Argument validation
    if (batch_size == 0) {
        throw std::invalid_argument("batch_size must be strictly greater than zero.");
    }

    if (max_node_id == std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("max_node_id is too large and would cause integer overflow.");
    }

    read_count = 0;

    std::ifstream gam_stream(gam_file, std::ios::in | std::ios::binary);
    if (!gam_stream) {
        throw std::runtime_error("Failed to open GAM file: " + gam_file);
    }

    const std::size_t coverage_size = max_node_id + 1;
    global_coverage.assign(coverage_size, 0);

    // Caps workers to a maximum of 8 threads to prevent thread contention on small tasks/Apple Silicon
    const int available_threads = omp_get_max_threads();
    const int active_threads = std::min(available_threads, 8);
    omp_set_num_threads(active_threads);

    // Thread-local accumulation storage to prevent atomic locks/data races during processing
    std::vector<std::vector<uint32_t>> local_coverages(
        static_cast<std::size_t>(active_threads),
        std::vector<uint32_t>(coverage_size, 0)
    );

    std::vector<uint64_t> local_read_counts(static_cast<std::size_t>(active_threads), 0);
    std::vector<uint64_t> local_out_of_range_counts(static_cast<std::size_t>(active_threads), 0);

    std::exception_ptr reader_exception;

    try {
        // Asynchronous BGZF block decompression stream iterator
        vg::io::MessageIterator message_it(gam_stream, false, decompression_threads);

        #pragma omp parallel
        {
            #pragma omp single
            {
                try {
                    // Producer loop: Main thread extracts Protobuf payloads and delegates work via tasks
                    while (message_it.has_current()) {
                        auto batch = std::make_unique<std::vector<std::string>>();
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

                        if (batch->empty()) {
                            continue;
                        }

                        auto raw_batch_ptr = batch.release();

                        // Consumer task: Workers deserialize Protobuf objects and calculate coverage
                        #pragma omp task firstprivate(raw_batch_ptr)
                        {
                            std::unique_ptr<std::vector<std::string>> owned_batch(raw_batch_ptr);
                            const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());

                            uint32_t* const coverage = local_coverages[tid].data();
                            uint64_t parsed_reads = 0;
                            uint64_t out_of_range = 0;

                            vg::Alignment alignment;

                            for (const std::string& serialized : *owned_batch) {
                                if (serialized.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                                    continue;
                                }

                                if (!alignment.ParseFromArray(serialized.data(), static_cast<int>(serialized.size()))) {
                                    continue;
                                }

                                ++parsed_reads;
                                const vg::Path& path = alignment.path();

                                for (int i = 0; i < path.mapping_size(); ++i) {
                                    const int64_t signed_node_id = path.mapping(i).position().node_id();
                                    if (signed_node_id <= 0) {
                                        continue;
                                    }

                                    const auto node_id = static_cast<std::size_t>(signed_node_id);
                                    if (node_id < coverage_size) {
                                        ++coverage[node_id];
                                    } else {
                                        ++out_of_range;
                                    }
                                }
                            }

                            local_read_counts[tid] += parsed_reads;
                            local_out_of_range_counts[tid] += out_of_range;
                        }
                    }
                } catch (...) {
                    reader_exception = std::current_exception();
                }

                // Wait for all worker tasks to complete before exiting the single region
                #pragma omp taskwait
            }
        }
    } catch (...) {
        reader_exception = std::current_exception();
    }

    if (reader_exception) {
        try {
            std::rethrow_exception(reader_exception);
        } catch (const std::exception& error) {
            throw std::runtime_error("Failed to parse GAM file '" + gam_file + "': " + error.what());
        }
    }

    // Aggregate thread-local read counters and check for index boundary errors
    uint64_t out_of_range_count = 0;
    for (std::size_t tid = 0; tid < local_read_counts.size(); ++tid) {
        read_count += local_read_counts[tid];
        out_of_range_count += local_out_of_range_counts[tid];
    }

    if (out_of_range_count > 0) {
        throw std::runtime_error(
            std::to_string(out_of_range_count) +
            " occurrence(s) of node_id exceeded max_node_id (" + std::to_string(max_node_id) + ")"
        );
    }

    // 2. Parallel SIMD Reduction: Combine thread-local coverage vectors into global_coverage
    #pragma omp parallel
    {
        const std::size_t worker = static_cast<std::size_t>(omp_get_thread_num());
        const std::size_t workers = static_cast<std::size_t>(omp_get_num_threads());

        const std::size_t block = (coverage_size + workers - 1) / workers;
        const std::size_t begin = worker * block;
        const std::size_t end = std::min(begin + block, coverage_size);

        if (begin < end) {
            uint32_t* const destination = global_coverage.data();

            for (const auto& source : local_coverages) {
                const uint32_t* const input = source.data();

                #pragma omp simd
                for (std::size_t i = begin; i < end; ++i) {
                    destination[i] += input[i];
                }
            }
        }
    }
}

/**
 * Lit un fichier GAM et affiche les informations du premier alignment trouvé.
 *
 * @param gam_file Le chemin vers le fichier GAM.
 * @return True si au moins un alignment a été lu et affiché avec succès, false sinon.
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
                    continue; // Ignore les messages d'en-tête/métadonnées vides
                }

                vg::Alignment alignment;
                if (!alignment.ParseFromString(*tag_and_data.second)) {
                    std::cerr << "Warning: Could not parse alignment body.\n";
                    continue;
                }

                // Premier Alignment valide trouvé : affichage des détails
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

                return true; // Succès : traitement terminé dès le premier alignment
            }
        }
    } catch (const std::exception &error) {
        std::cerr << "Error while reading GAM file: " << error.what() << '\n';
        return false;
    }

    std::cerr << "The GAM file contains no alignments.\n";
    return false;
}
