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

void process_gam(
    const std::string& gam_file,
    std::vector<cdx::Coverage>& target,
    const cdx::Nid nid_min,
    std::uint64_t& read_count,
    const std::size_t batch_size,
    const int decompression_threads
) {
    // 1. Argument validation
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

    read_count = 0;

    std::ifstream gam_stream(
        gam_file,
        std::ios::in | std::ios::binary
    );

    if (!gam_stream) {
        throw std::runtime_error(
            "Failed to open GAM file: " + gam_file
        );
    }

    const std::size_t coverage_size = target.size();

    if (coverage_size == 0) {
        return;
    }

    /*
     * Preserve which local nodes belong to the query.
     *
     * We cannot use target directly from the worker tasks because target
     * will receive the final coverage during the reduction.
     */
    std::vector<std::uint8_t> valid_nodes(coverage_size, 0);

    for (std::size_t node_offset = 0;
         node_offset < coverage_size;
         ++node_offset) {

        if (target[node_offset] < cfg::NOT_IN_QUERY) {
            valid_nodes[node_offset] = 1;
        }
    }

    // Cap workers to a maximum of 8 threads.
    const int available_threads = omp_get_max_threads();
    const int active_threads = std::max(
        1,
        std::min(available_threads, 8)
    );

    /*
     * Avoid changing the global OpenMP thread setting with
     * omp_set_num_threads(). The num_threads clause only applies to
     * this region.
     */

    std::vector<std::vector<std::uint32_t>> local_coverages(
        static_cast<std::size_t>(active_threads),
        std::vector<std::uint32_t>(coverage_size, 0)
    );

    std::vector<std::uint64_t> local_read_counts(
        static_cast<std::size_t>(active_threads),
        0
    );

    /*
     * These counters are informational. Nodes outside the local component
     * or outside the query are expected and are not errors.
     */
    std::vector<std::uint64_t> local_outside_range_counts(
        static_cast<std::size_t>(active_threads),
        0
    );

    std::vector<std::uint64_t> local_outside_query_counts(
        static_cast<std::size_t>(active_threads),
        0
    );

    std::exception_ptr reader_exception;

    try {
        vg::io::MessageIterator message_it(
            gam_stream,
            false,
            decompression_threads
        );

        #pragma omp parallel num_threads(active_threads)
        {
            #pragma omp single
            {
                try {
                    while (message_it.has_current()) {
                        auto batch =
                            std::make_unique<std::vector<std::string>>();

                        batch->reserve(batch_size);

                        while (
                            message_it.has_current() &&
                            batch->size() < batch_size
                        ) {
                            auto tag_and_data =
                                std::move(message_it.take());

                            if (!vg::io::Registry::check_protobuf_tag<
                                    vg::Alignment
                                >(tag_and_data.first)) {
                                continue;
                            }

                            if (tag_and_data.second) {
                                batch->push_back(
                                    std::move(*tag_and_data.second)
                                );
                            }
                        }

                        if (batch->empty()) {
                            continue;
                        }

                        auto* raw_batch_ptr = batch.release();

                        #pragma omp task firstprivate(raw_batch_ptr)
                        {
                            std::unique_ptr<std::vector<std::string>>
                                owned_batch(raw_batch_ptr);

                            const std::size_t tid =
                                static_cast<std::size_t>(
                                    omp_get_thread_num()
                                );

                            std::uint32_t* const coverage =
                                local_coverages[tid].data();

                            std::uint64_t parsed_reads = 0;
                            std::uint64_t outside_range = 0;
                            std::uint64_t outside_query = 0;

                            vg::Alignment alignment;

                            for (const std::string& serialized :
                                 *owned_batch) {

                                if (
                                    serialized.size() >
                                    static_cast<std::size_t>(
                                        std::numeric_limits<int>::max()
                                    )
                                ) {
                                    continue;
                                }

                                alignment.Clear();

                                if (!alignment.ParseFromArray(
                                        serialized.data(),
                                        static_cast<int>(
                                            serialized.size()
                                        )
                                    )) {
                                    continue;
                                }

                                ++parsed_reads;

                                const vg::Path& path = alignment.path();

                                for (int i = 0;
                                     i < path.mapping_size();
                                     ++i) {

                                    const std::int64_t signed_node_id =
                                        path.mapping(i)
                                            .position()
                                            .node_id();

                                    if (signed_node_id <= 0) {
                                        continue;
                                    }

                                    const auto node_id =
                                        static_cast<cdx::Nid>(
                                            signed_node_id
                                        );

                                    /*
                                     * The GAM uses global node IDs.
                                     * target uses local offsets beginning
                                     * at nid_min.
                                     */
                                    if (node_id < nid_min) {
                                        ++outside_range;
                                        continue;
                                    }

                                    const cdx::Nid raw_offset =
                                        node_id - nid_min;

                                    if (
                                        raw_offset >=
                                        static_cast<cdx::Nid>(
                                            coverage_size
                                        )
                                    ) {
                                        ++outside_range;
                                        continue;
                                    }

                                    const std::size_t node_offset =
                                        static_cast<std::size_t>(
                                            raw_offset
                                        );

                                    /*
                                     * The node belongs to the component,
                                     * but not necessarily to the requested
                                     * genomic interval.
                                     */
                                    if (!valid_nodes[node_offset]) {
                                        ++outside_query;
                                        continue;
                                    }

                                    ++coverage[node_offset];
                                }
                            }

                            local_read_counts[tid] += parsed_reads;
                            local_outside_range_counts[tid] +=
                                outside_range;
                            local_outside_query_counts[tid] +=
                                outside_query;
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
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Failed to parse GAM file '" +
                gam_file +
                "': " +
                error.what()
            );
        }
    }

    std::uint64_t outside_range_count = 0;
    std::uint64_t outside_query_count = 0;

    for (std::size_t tid = 0;
         tid < local_read_counts.size();
         ++tid) {

        read_count += local_read_counts[tid];

        outside_range_count +=
            local_outside_range_counts[tid];

        outside_query_count +=
            local_outside_query_counts[tid];
    }

    /*
     * Do not throw for nodes outside the local component.
     *
     * A GAM may legitimately contain alignments from every component,
     * while target represents only one component or one query.
     */
#ifdef CDX_GAM_DEBUG
    std::cerr
        << "[DEBUG] GAM mappings outside local component: "
        << outside_range_count
        << '\n'
        << "[DEBUG] GAM mappings inside component but outside query: "
        << outside_query_count
        << '\n';
#endif

    /*
     * Combine thread-local vectors into target.
     *
     * Only valid query nodes are overwritten. NOT_IN_QUERY values remain
     * unchanged.
     */
    #pragma omp parallel num_threads(active_threads)
    {
        const std::size_t worker =
            static_cast<std::size_t>(omp_get_thread_num());

        const std::size_t workers =
            static_cast<std::size_t>(omp_get_num_threads());

        const std::size_t block =
            (coverage_size + workers - 1) / workers;

        const std::size_t begin = worker * block;
        const std::size_t end =
            std::min(begin + block, coverage_size);

        for (std::size_t node_offset = begin;
             node_offset < end;
             ++node_offset) {

            if (!valid_nodes[node_offset]) {
                continue;
            }

            std::uint64_t total_coverage = 0;

            for (const auto& source : local_coverages) {
                total_coverage += source[node_offset];
            }

            /*
             * Preserve the behavior of the existing uint32_t storage,
             * but detect an actual coverage overflow.
             */
            if (
                total_coverage >
                static_cast<std::uint64_t>(
                    std::numeric_limits<cdx::Coverage>::max()
                )
            ) {
                /*
                 * An exception cannot safely escape an OpenMP region.
                 * Saturating here is one possible policy.
                 *
                 * If saturation is undesirable, collect an overflow flag
                 * and throw after the parallel region.
                 */
                target[node_offset] =
                    std::numeric_limits<cdx::Coverage>::max();
            } else {
                target[node_offset] =
                    static_cast<cdx::Coverage>(total_coverage);
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
