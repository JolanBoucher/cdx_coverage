//
// Created by Jolan on 2026-07-31.
//

#ifndef GAM_IO_H
#define GAM_IO_H

#include "cdx_types.h"
#include <string>
#include <vector>
#include <cstdint>


/**
 * @brief Rapidly computes node coverage from a GAM (Genomic Alignment Map) file.
 *
 * Uses a Producer-Consumer pattern with OpenMP tasks:
 * - Thread 0 (Producer) streams and deserializes Protobuf batches from libvgio.
 * - Worker Threads (Consumers) parse alignment paths into thread-local coverage vectors.
 * - Thread-local vectors are reduced in parallel at the end to prevent mutex contention.
 *
 * @param gam_file Path to the input .gam file.
 * @param target
 * @param nid_min
 * @param read_count Output reference updated with total reads processed.
 * @param batch_size Number of GAM alignments packaged per OpenMP task (default: 512).
 * @param decompression_threads BGZF decompression threads for libvgio (default: 4).
 */
void process_gam(
    const std::string& gam_file,
    std::vector<cdx::Coverage>& target,
    cdx::Nid nid_min,
    std::uint64_t& read_count,
    std::size_t batch_size,
    int decompression_threads
);

/**
 * Inspecte et affiche dans la console les détails du premier alignment d'un fichier GAM.
 *
 * @param gam_file Chemin vers le fichier GAM à inspecter.
 * @return true si un alignment a été lu et affiché avec succès, false sinon.
 */
bool inspect_first_alignment(const std::string& gam_file);

#endif // GAM_IO_H