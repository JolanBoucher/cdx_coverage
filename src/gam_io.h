#ifndef GAM_IO_H
#define GAM_IO_H

#include "cdx_types.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

/**
 * @brief Statistiques produites pendant la lecture d'un fichier GAM.
 *
 * Les quatre premiers compteurs sont calculés par alignment/read.
 * Les deux derniers sont des compteurs diagnostiques calculés par mapping.
 */
struct GamMappingStats {
    std::uint64_t total = 0;                       // Nombre total d'alignements correctement désérialisés.
    std::uint64_t mapped = 0;                      // Nombre d'alignements contenant au moins un node ID valide.
    std::uint64_t mapped_to_query = 0;             // Nombre d'alignements touchant au moins un nœud de la query active.
    std::uint64_t unmapped = 0;                    // Nombre d'alignements ne contenant aucun node ID valide.

    GamMappingStats& operator+=(const GamMappingStats& other) noexcept {
        total += other.total;
        mapped += other.mapped;
        mapped_to_query += other.mapped_to_query;
        unmapped += other.unmapped;
        return *this;
    }
};

/**
 * @brief Calcule la couverture locale et les statistiques d'un fichier GAM.
 *
 * @param gam_file Chemin vers le fichier GAM.
 * @param target Vecteur de couverture dans l'espace local CDX.
 * @param nid_min Premier node ID global représenté par target.
 * @param batch_size Nombre maximal d'alignements par tâche OpenMP.
 * @param decompression_threads Nombre de threads de décompression BGZF.
 *
 * @return Statistiques calculées pendant le traitement du GAM.
 */
[[nodiscard]]
GamMappingStats process_gam(
    const std::string& gam_file,
    std::vector<cdx::Coverage>& target,
    cdx::Nid nid_min,
    std::size_t batch_size,
    int decompression_threads
);

/**
 * @brief Inspecte le premier alignement d'un fichier GAM.
 */
[[nodiscard]] bool inspect_first_alignment(const std::string& gam_file);

#endif // GAM_IO_H