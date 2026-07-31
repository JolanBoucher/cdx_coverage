#include "gam_io.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
#include <omp.h>

int main(int argc, char** argv) {
    const std::string gam_file = (argc > 1) ? argv[1] : "./bio-file/messy.gam";

    // Bornes connues pour le graphe mitochondrial (ex: 35 000 max)
    // À ajuster selon le pangenome (ex: 1 000 000 pour des graphes plus grands)
    const std::size_t max_node_id = 550000;

    // taille de batch le plus rapide
    //      512 sur petit gam de <500k read
    //      2048 sur gam à grand 500k-25M
    const std::size_t batch_size = 2048;
    const int nb_runs = 5;
    int decompression_threads = 4;

    std::cout << "=== Benchmarking process_gam_fast (Optimisé) ===\n";
    std::cout << "Fichier          : " << gam_file << "\n";
    std::cout << "Threads OpenMP   : " << omp_get_max_threads() << " / " << omp_get_num_procs() << " cœurs\n";
    std::cout << "Max Node ID      : " << max_node_id << "\n";
    std::cout << "Batch Size       : " << batch_size << "\n";
    std::cout << "Répétitions      : " << nb_runs << " exécution(s)\n";

    std::vector<uint32_t> global_coverage;
    uint64_t read_count = 0;

    // Run de chauffe (Warmup)
    std::cout << "Exécution de chauffe (warmup)... " << std::flush;
    try {
        process_gam_fast(gam_file, global_coverage, read_count, max_node_id, batch_size);
        std::cout << "Terminé.\n";
    } catch (const std::exception& e) {
        std::cerr << "\n[Erreur Warmup] " << e.what() << "\n";
        return 1;
    }

    // Benchmark
    std::cout << "Lancement du benchmark sur " << nb_runs << " itération(s)...\n";
    std::vector<double> timings_ms;
    timings_ms.reserve(nb_runs);

    for (int i = 0; i < nb_runs; ++i) {
        global_coverage.clear();
        read_count = 0;

        auto start = std::chrono::high_resolution_clock::now();

        process_gam_fast(
            gam_file,
            global_coverage,
            read_count,
            max_node_id,
            batch_size,
            decompression_threads // Ou laisser par défaut si non spécifié
        );

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;
        timings_ms.push_back(duration.count());
    }

    // Analyse des métriques
    double total_time_ms = std::accumulate(timings_ms.begin(), timings_ms.end(), 0.0);
    double avg_time_ms = total_time_ms / nb_runs;

    std::size_t covered_nodes = 0;
    uint64_t total_hits = 0;
    uint32_t max_cov = 0;
    std::size_t max_cov_node_id = 0;

    for (std::size_t node_id = 0; node_id < global_coverage.size(); ++node_id) {
        uint32_t cov = global_coverage[node_id];
        if (cov > 0) {
            covered_nodes++;
            total_hits += cov;
            if (cov > max_cov) {
                max_cov = cov;
                max_cov_node_id = node_id;
            }
        }
    }

    std::cout << "\n--- Résultats ---\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Temps total (" << nb_runs << " runs) : " << total_time_ms << " ms (" << total_time_ms / 1000.0 << " s)\n";
    std::cout << "Temps moyen / run  : " << avg_time_ms << " ms (" << avg_time_ms / 1000.0 << " s)\n";
    std::cout << "Reads parsés / run : " << read_count << "\n";
    std::cout << "Nœuds couverts     : " << covered_nodes << "\n";
    std::cout << "Total des hits     : " << total_hits << "\n";
    std::cout << "Couverture max     : " << max_cov << " (Nœud ID: " << max_cov_node_id << ")\n";
    std::cout << "=====================================\n";

    return 0;
}