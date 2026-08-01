// contain only the data structure associated with the cdx domain

#ifndef CDX_COVERAGE_CDX_TYPES_H
#define CDX_COVERAGE_CDX_TYPES_H

#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace cdx {

    // defining some types for this namespace
    using Nid = std::uint32_t;
    using Idx = std::uint32_t;
    using FlatIdx = std::uint64_t;
    using PosBp = std::uint64_t;
    using Coverage = std::uint32_t;

    //
    // Bitvector + index de rang
    //
    struct PositionIndex {
        std::vector<std::uint64_t> bitvector;
        std::vector<std::uint32_t> rank;
    };


    //
    // Informations globales calculées avant la construction
    //
    struct GraphLayout {
        Nid graph_nid_min = 0;                                  // min nid du graphe
        Nid graph_nid_max = 0;                                  // max nid du graph
        Nid total_nodes = 0;                                    // nb_node du graphe
        std::size_t component_count = 0;                      // nb_compo du graphe
        std::vector<PosBp> component_offsets;           // offset de chaque compo
    };

    struct ComponentInfo {
        std::size_t compo_id; // à voir si c'est la bonne taille
        std:: string compo_name;
        Nid nid_min;
        Nid nid_max;
        Nid nb_nodes;
        PosBp component_length;
        std::uint64_t payload_offset;
        std::uint64_t payload_size;
    };


    //
    // Résultat construit par loadQuery()
    //
    struct QueryData {
        std::vector<Idx> nid2idx;                                   // nid -> idx
        std::vector<PosBp> idx2pb;                                  // idx -> position cumulative
        std::vector<Coverage> node_coverage;                        // table de couverture locale
        PositionIndex position_index;                               // projection bp -> idx
        std::pair<PosBp, PosBp>  query_range_bp {0, 0};             // bornes bp
        std::pair<Idx, Idx> query_range_idx {0, 0};                 // bornes idx
        Nid nid_min = 0;                                            // plus petit nid du composant
        Nid node_count = 0;                                         // nb_node du composant
        PosBp component_length = 0;                                 // longueur totale du composant
    };


    //
    // Résultat construit par loadGlobal()
    //
    struct GlobalData {
        std::vector<FlatIdx> nid2flat_idx;              // nid -> flat_idx
        std::vector<Coverage> node_coverage;            // couverture globale
        std::vector<PosBp> idx2pb;                      // flat_idx -> position cumulée
        std::vector<PosBp> component_offsets;           // offsets des compo
        std::vector<PosBp> component_lengths;           // longueur cumulative des composants
        Nid graph_nid_min = 0;                          // min nid du graphe
        Nid graph_nid_max = 0;                          // max nid du graphe
        Nid total_nodes = 0;                            // nb_node du graphe
    };
} // namespace cdx

#endif //CDX_COVERAGE_CDX_TYPES_H