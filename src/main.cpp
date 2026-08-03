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
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
    constexpr std::size_t GAM_BATCH_SIZE = 2048;
    constexpr int GAM_DECOMPRESSION_THREADS = 4;

    /**
     * @brief Calcule la couverture GAM et retourne les statistiques de mapping.
     */
    [[nodiscard]]
    GamMappingStats processGam(
        const std::string &gam_path,
        const cdx::Nid nid_min,
        std::vector<cdx::Coverage> &target
    ) {
        return process_gam(
            gam_path,
            target,
            nid_min,
            GAM_BATCH_SIZE,
            GAM_DECOMPRESSION_THREADS
        );
    }

    [[nodiscard]]
    GamMappingStats runGlobalPipeline(const CliArgs &args) {
        cdx::GlobalData data;
        {
            cfg::ScopedTimer timer("loadGlobal");
            data = cdx::loadGlobal(args.cdx_file);
        }

        GamMappingStats mapping_stats;
        {
            cfg::ScopedTimer timer("processGam");
            mapping_stats = processGam(
                args.gam_file,
                data.layout.graph_nid_min,
                data.node_coverage
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
                mapping_stats
            );
        }

        return mapping_stats;
    }

    [[nodiscard]]
    GamMappingStats runQueryPipeline(
        const CliArgs &args,
        const std::size_t target_cid
    ) {
        // Conversion de la plage optionnelle issue des arguments CLI
        std::optional<std::pair<std::int64_t, std::int64_t> > query_range = std::nullopt;
        if (args.query && args.query->range) {
            query_range = std::make_pair(
                args.query->range->start,
                args.query->range->end
            );
        }

        const bool is_circular = args.component_type == ComponentType::Circular;

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

        GamMappingStats mapping_stats;
        {
            cfg::ScopedTimer timer("processGam");
            mapping_stats = processGam(
                args.gam_file,
                data.component.nid_min,
                data.node_coverage
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
        return mapping_stats;
    }
} // anonymous namespace

int main(int argc, char **argv) {
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
            std::optional<cdx::Cid> component_id = std::nullopt;

            if (args.inspect.component.has_value()) {
                const ResolvedComponent resolved =
                        resolver.resolve(*args.inspect.component);

                if (resolved.cid > static_cast<std::size_t>(std::numeric_limits<cdx::Cid>::max())) {
                    throw std::overflow_error(
                        "Resolved component ID exceeds cdx::Cid capacity."
                    );
                }

                component_id = static_cast<cdx::Cid>(resolved.cid);
            }

            cdx::inspectComponent(
                std::filesystem::path(args.cdx_file),
                component_id
            );
            return 0;
        }

        // 4. Traitement du PIPELINE DE COUVERTURE
        GamMappingStats mapping_stats;

        if (args.query.has_value()) {
            const ResolvedComponent resolved = resolver.resolve(args.query->component);
            mapping_stats = runQueryPipeline(args, resolved.cid);
        } else {
            mapping_stats = runGlobalPipeline(args);
        }

        std::cout << "Processing completed successfully.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
