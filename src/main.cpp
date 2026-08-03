#include "cdx_loader.h"
#include "cli.hpp"
#include "config.h"
#include "cov_projection.h"
#include "gam_io.h"
#include "output_coverage.h"
#include "output_stats.h"
#include "query_resolver.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using MappingStats = std::map<std::string, std::uint64_t>;

constexpr std::size_t GAM_BATCH_SIZE = 2048;
constexpr int GAM_DECOMPRESSION_THREADS = 4;

/**
 * @brief Calcule directement la couverture GAM dans l'espace local CDX.
 *
 * @param gam_path Chemin vers le fichier GAM.
 * @param nid_min Premier node ID global représenté par target.
 * @param target Vecteur local de couverture CDX. Les valeurs NOT_IN_QUERY
 *               sont conservées et ignorées pendant le calcul.
 * @param mapping_stats Statistiques de mapping.
 */
void processGam(
    const std::string& gam_path,
    const cdx::Nid nid_min,
    std::vector<cdx::Coverage>& target,
    MappingStats& mapping_stats
) {
    std::uint64_t read_count = 0;

    process_gam(
        gam_path,
        target,
        nid_min,
        read_count,
        GAM_BATCH_SIZE,
        GAM_DECOMPRESSION_THREADS
    );

    mapping_stats["total"] = read_count;
    mapping_stats["mapped"] = read_count;
    mapping_stats["mapped_to_query"] = 0;
    mapping_stats["unmapped"] = 0;
}

void runGlobalPipeline(const CliArgs& args, MappingStats& mapping_stats) {
    cdx::GlobalData data;
    {
        cfg::ScopedTimer timer("loadGlobal");
        data = cdx::loadGlobal(args.cdx_file);
    }

    {
        cfg::ScopedTimer timer("processGam");
        processGam(
            args.gam_file,
            data.layout.graph_nid_min,
            data.node_coverage,
            mapping_stats
        );
    }

    std::vector<cdx::Coverage> flat_idx_coverage;
    {
        cfg::ScopedTimer timer("projectCov2IdxGlobal");
        flat_idx_coverage = projectCov2IdxGlobal(
            data.node_coverage,
            data.nid2flat_idx,
            data.layout.component_offsets
        );
    }

    std::vector<cdx::Coverage> flat_bp_coverage;
    std::vector<cdx::PosBp> bp_component_offsets;
    {
        cfg::ScopedTimer timer("expandPosCovGlobal");
        std::tie(flat_bp_coverage, bp_component_offsets) = expandPosCovGlobal(
            flat_idx_coverage,
            data.layout.component_offsets,
            data.idx2bp,
            data.idx2bp_offsets
        );
    }

    const std::filesystem::path out_dir(args.output_directory);

    if (args.generateTable()) {
        cfg::ScopedTimer timer("writeCoverageTsvGlobal");
        output::writeCoverageTsvGlobal(
            out_dir / cfg::NAME_TSV_FILE,
            flat_bp_coverage,
            bp_component_offsets,
            data.layout.component_names
        );
    }

    if (args.generateStats()) {
        cfg::ScopedTimer timer("writeStatsReportGlobal");
        output::writeStatsReportGlobal(
            out_dir / cfg::NAME_STATS_FILE,
            flat_bp_coverage,
            bp_component_offsets,
            data.layout.component_names,
            &mapping_stats
        );
    }
}

void runQueryPipeline(
    const CliArgs& args,
    std::size_t target_cid,
    MappingStats& mapping_stats
) {
    // Conversion de la plage optionnelle issue des arguments CLI
    std::optional<std::pair<std::int64_t, std::int64_t>> query_range = std::nullopt;
    if (args.query && args.query->range) {
        query_range = std::make_pair(
            args.query->range->start,
            args.query->range->end
        );
    }

    if (query_range) {
        std::cerr << "[DEBUG] - Query range: "
                  << query_range->first << ':' << query_range->second << '\n';
    } else {
        std::cerr << "[DEBUG] - Query range: entire component\n";
    }

    const bool is_circular = (args.component_type == ComponentType::Circular);

    cdx::QueryData data;
    {
        cfg::ScopedTimer timer("loadQuery");
        data = cdx::loadQuery(
            args.cdx_file,
            static_cast<cdx::Cid>(target_cid),
            query_range,
            is_circular
        );
    }

    {
        cfg::ScopedTimer timer("processGam");
        processGam(
            args.gam_file,
            data.component.nid_min,
            data.node_coverage,
            mapping_stats
        );
    }

    std::vector<cdx::Coverage> idx_coverage;
    {
        cfg::ScopedTimer timer("projectCov2IdxQuery");
        idx_coverage = projectCov2IdxQuery(
            data.node_coverage,
            data.nid2idx,
            data.idx2bp.size() - 1
        );
    }

    std::vector<cdx::Coverage> bp_coverage;
    {
        cfg::ScopedTimer timer("expandPosCovQuery");
        bp_coverage = expandPosCovQuery(idx_coverage, data.idx2bp);
    }

    {
        cfg::ScopedTimer timer("trimCoverageToQuery");
        bp_coverage = trimCoverageToQuery(bp_coverage, data.query_range_bp);
    }

    const std::filesystem::path out_dir(args.output_directory);

    if (args.generateTable()) {
        cfg::ScopedTimer timer("writeCoverageTsvQuery");
        output::writeCoverageTsvQuery(
            out_dir / cfg::NAME_TSV_FILE,
            bp_coverage,
            data.component.compo_name
        );
    }

    if (args.generateStats()) {
        cfg::ScopedTimer timer("writeStatsReportQuery");
        output::writeStatsReportQuery(
            out_dir / cfg::NAME_STATS_FILE,
            mapping_stats,
            bp_coverage,
            data.component.compo_name
        );
    }
}

} // anonymous namespace

int main(int argc, char** argv) {
    try {
        cfg::ScopedTimer total_timer("TOTAL PROGRAM");

        // 1. Parsing CLI
        CliArgs args = parse_args(argc, argv);

        std::filesystem::create_directories(args.output_directory);

        // 2. Chargement léger de la métadonnée CDX pour alimenter le resolver
        ComponentResolver resolver;
        {
            cdx::GlobalData meta_data = cdx::loadGlobal(args.cdx_file);
            for (std::size_t cid = 0; cid < meta_data.layout.component_names.size(); ++cid) {
                resolver.register_component(cid, meta_data.layout.component_names[cid]);
            }
        }

        // 3. Traitement du MODE INSPECT
        if (args.inspectMode()) {
            if (args.inspect.component.has_value()) {
                ResolvedComponent resolved = resolver.resolve(*args.inspect.component);
                std::cout << "Component Inspection:\n"
                          << "  CID  : " << resolved.cid << "\n"
                          << "  Name : " << resolved.name << "\n";
            } else {
                std::cout << "CDX Index contains " << resolver.size() << " components.\n";
                for (std::size_t cid = 0; cid < resolver.size(); ++cid) {
                    std::cout << "  - CID " << cid << " : " << resolver.get_name(cid) << "\n";
                }
            }
            return 0;
        }

        // 4. Traitement du PIPELINE DE COUVERTURE
        MappingStats mapping_stats{
            {"total", 0},
            {"mapped", 0},
            {"mapped_to_query", 0},
            {"unmapped", 0}
        };

        if (args.query.has_value()) {
            ResolvedComponent resolved = resolver.resolve(args.query->component);

            std::clog << "[INFO] Query target resolved: '" << args.query->component
                      << "' -> CID " << resolved.cid << " (" << resolved.name << ")\n";

            runQueryPipeline(args, resolved.cid, mapping_stats);
        } else {
            runGlobalPipeline(args, mapping_stats);
        }

        std::cout << "Processing completed successfully.\n";
        return 0;

    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}