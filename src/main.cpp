/**
 * @file main.cpp
 * @brief Main entry point for the CDX Coverage application.
 *
 * CDX Coverage computes per-position coverage statistics from GAM
 * alignments projected onto a graph indexed by a CDX file. Depending on
 * the command-line options, the program can:
 *
 *   - Analyze the entire pangenome (global mode).
 *   - Analyze a single component or genomic interval (query mode).
 *   - Inspect indexed components and metadata (inspection mode).
 *   - Generate coverage tables, statistical reports, and graphical
 *     visualizations.
 *
 * This file contains the top-level application logic responsible for
 * argument parsing, component resolution, pipeline dispatch, execution
 * timing, status reporting, and fatal-error handling.
 */

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
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "output_plot.h"
#include "query_plot_slice.h"

namespace {
    constexpr std::size_t GAM_BATCH_SIZE = 2048;
    constexpr std::size_t TERMINAL_WIDTH = 77;

   /**
    * @brief Execute the complete global coverage-analysis pipeline.
    *
    * Runs the full pangenome workflow from input loading through coverage
    * generation and output production. The pipeline:
    *
    * 1. Loads the global CDX index and component metadata.
    * 2. Processes GAM alignments and accumulates node-level coverage.
    * 3. Projects coverage into index space and expands it into genomic coordinate space.
    * 4. Generates the user-requested outputs:
    *    - Coverage table (TSV)
    *    - Coverage statistics report
    *    - Coverage graph (linear or circular)
    *
    * Progress information and execution timings are reported for each stage.
    * Output generation steps are executed only when enabled through the
    * corresponding command-line arguments.
    *
    * @param args Command-line configuration controlling input files,
    *   processing parameters, output selection, and rendering options.
    *
    * @throws std::exception Propagates any exception raised during data
    *   loading, coverage computation, output generation, or rendering.
    */
    void runGlobalPipeline(const CliArgs &args) {
        const int output_steps =
                static_cast<int>(args.generateTable()) +
                static_cast<int>(args.generateStats()) +
                static_cast<int>(args.generateGraph());

        const int total_steps = 3 + output_steps;
        int current_step = 1;

        // ============================================================
        // STEP 1: CDX Loading
        // ============================================================
        std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] CDX Loading\n";
        cdx::GlobalData data;

        {
            cfg::ScopedTimer timer("Global CDX data loaded");
            data = cdx::loadGlobal(args.cdx_file);
            timer.update_name("Whole pangenome loaded");
        }

        // Base-precision mode needs a per-node length lookup before
        // process_gam can be called (see coverage_gaps.h/gam_io.h) - built
        // only when actually needed, so Node-mode runs never pay for it.
        const bool base_precision = (args.coverage_precision == CoveragePrecision::Base);
        std::vector<cdx::SeqLen> node_lengths;
        std::vector<BpGap> coverage_gaps;

        if (base_precision) {
            cfg::ScopedTimer timer("Node length lookup built");
            node_lengths = buildNodeLengthsGlobal(
                data.nid2flat_idx,
                data.layout.component_offsets,
                data.idx2bp,
                data.idx2bp_offsets
            );
        }

        // ============================================================
        // STEP 2: GAM Processing
        // ============================================================
        std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] GAM Processing\n";
        GamMappingStats mapping_stats;

        {
            cfg::ScopedTimer timer("GAM alignments processed");

            mapping_stats = process_gam(
                args.gam_file,
                data.node_coverage,
                data.layout.graph_nid_min,
                GAM_BATCH_SIZE,
                args.decompression_threads,
                args.worker_threads,
                args.coverage_precision,
                base_precision ? &node_lengths : nullptr,
                base_precision ? &coverage_gaps : nullptr
            );

            timer.update_name("Processed " + cfg::formatInteger(mapping_stats.total) + " reads");
        }

        // ============================================================
        // STEP 3: Coverage Projection
        // ============================================================
        std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Genomic Coordinate Projection\n";
        std::vector<cdx::Coverage> flat_idx_coverage;

        {
            cfg::ScopedTimer timer("Node coverage projected to component index space");

            flat_idx_coverage = projectCov2IdxGlobal(
                data.node_coverage,
                data.nid2flat_idx,
                data.layout.component_offsets
            );
        }

        std::vector<cdx::Coverage> flat_bp_coverage;
        std::vector<cdx::PosBp> bp_component_offsets;

        {
            cfg::ScopedTimer timer("Index coverage expanded to genomic coordinates");

            std::tie(flat_bp_coverage, bp_component_offsets) =
                    expandPosCovGlobal(
                        flat_idx_coverage,
                        data.layout.component_offsets,
                        data.idx2bp,
                        data.idx2bp_offsets
                    );
        }

        // Refine the uniform per-node fill above down to base-pair
        // precision by subtracting every collected coverage gap. Applied
        // here - after expansion, before any output is written - so every
        // output (TSV, stats, graph) automatically reflects base-pair
        // precision without needing to know it exists.
        if (base_precision) {
            cfg::ScopedTimer timer("Coverage gaps applied at base-pair precision");
            applyBpGapsGlobal(
                flat_bp_coverage,
                coverage_gaps,
                data.nid2flat_idx,
                data.layout.component_offsets,
                data.idx2bp,
                data.idx2bp_offsets,
                bp_component_offsets
            );
        }

        const std::filesystem::path output_directory(args.output_directory);

        // ============================================================
        // OPTIONAL OUTPUT: TSV
        // ============================================================
        if (args.generateTable()) {
            std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Coverage Table Export\n";

            {
                const std::filesystem::path output_path = output_directory / cfg::NAME_TSV_FILE;
                cfg::ScopedTimer timer("Coverage TSV table written");

                output::writeCoverageTsvGlobal(
                    output_path,
                    flat_bp_coverage,
                    bp_component_offsets,
                    data.layout.component_names
                );

                timer.update_name("Coverage table saved (" + output_path.filename().string() + ')');
            }
        }

        // ============================================================
        // OPTIONAL OUTPUT: Statistics
        // ============================================================
        if (args.generateStats()) {
            std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Statistics Report\n";

            {
                const std::filesystem::path output_path = output_directory / cfg::NAME_STATS_FILE;
                cfg::ScopedTimer timer("Coverage statistics report written");

                output::writeStatsReportGlobal(
                    output_path,
                    flat_bp_coverage,
                    bp_component_offsets,
                    data.layout.component_names,
                    mapping_stats,
                    args.worker_threads,
                    args.coverage_precision
                );

                timer.update_name("Statistics report saved (" + output_path.filename().string() + ')');
            }
        }

        // ============================================================
        // OPTIONAL OUTPUT: Graph
        // ============================================================
        if (args.generateGraph()) {
            std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Coverage Graph Generation\n";

            cfg::ScopedTimer timer("Coverage graph");

            output::PlotConfig plot_config;

            plot_config.smoothing = static_cast<double>(args.smoothing);
            plot_config.max_plot_points = static_cast<std::size_t>(args.max_plot_points);
            plot_config.dpi = args.dpi;
            plot_config.figure_width = args.figure_width;
            plot_config.figure_height = args.figure_height;
            plot_config.line_color = args.line_color;
            plot_config.fill_color = args.fill_color;
            if (args.log_base.has_value()) {
                plot_config.log_base = args.log_base.value();
            }

            if (args.component_type == ComponentType::Circular) {
                output::writeCircularPlotGlobal(
                    output_directory / cfg::NAME_GRAPH_FILE,
                    flat_bp_coverage,
                    bp_component_offsets,
                    data.layout.component_names,
                    plot_config
                );
            } else {
                output::writeLinearPlotGlobal(
                    output_directory / cfg::NAME_GRAPH_FILE,
                    flat_bp_coverage,
                    bp_component_offsets,
                    data.layout.component_names,
                    plot_config
                );
            }

            timer.update_name("Coverage graph saved (" + std::string(cfg::NAME_GRAPH_FILE) + ')');
        }
    }

    /**
     * @brief Execute the query-level coverage-analysis pipeline.
     *
     * Runs the complete workflow for a single component or component subregion.
     * The pipeline:
     *
     *   1. Loads the requested component (or query interval) from the CDX index.
     *   2. Processes GAM alignments and accumulates node-level coverage.
     *   3. Projects coverage into index space, expands it into genomic
     *      coordinates, and restricts the result to the requested query region.
     *   4. Generates the user-requested outputs:
     *        - Coverage table (TSV)
     *        - Coverage statistics report
     *        - Coverage graph (linear or circular)
     *
     * Query ranges may represent either a complete component or a specific
     * genomic interval. Circular queries may additionally span the component
     * origin, in which case the circular plotting backend preserves the full
     * component coverage context required for wrap-around smoothing and
     * visualization.
     *
     * Progress information and execution timings are reported for each stage.
     * Output-generation steps are executed only when enabled through the
     * corresponding command-line arguments.
     *
     * @param args Command-line configuration controlling input files,
     *        processing parameters, output selection, and rendering options.
     * @param target_cid Component identifier selected for analysis.
     * @param target_name Human-readable name of the selected component.
     *
     * @throws std::exception Propagates any exception raised during data
     *         loading, coverage computation, output generation, or rendering.
     */
    void runQueryPipeline(
        const CliArgs &args,
        const std::size_t target_cid,
        const std::string &target_name
    ) {
        const int output_steps =
                static_cast<int>(args.generateTable()) +
                static_cast<int>(args.generateStats()) +
                static_cast<int>(args.generateGraph());

        const int total_steps = 3 + output_steps;
        int current_step = 1;

        std::optional<std::pair<std::int64_t, std::int64_t> > query_range = std::nullopt;

        if (args.query && args.query->range) {
            query_range = std::make_pair(args.query->range->start, args.query->range->end);
        }

        const bool is_circular = args.component_type == ComponentType::Circular;

        // ============================================================
        // STEP 1: Query CDX Loading
        // ============================================================
        std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] CDX Query Loading\n";
        cdx::QueryData data;

        {
            cfg::ScopedTimer timer("Query component loaded");

            data = cdx::loadQuery(
                args.cdx_file,
                static_cast<cdx::Cid>(target_cid),
                query_range,
                is_circular
            );

            timer.update_name("Pangenome component " + target_name + " loaded");
        }

        // Base-precision mode needs a per-node length lookup before
        // process_gam can be called (see coverage_gaps.h/gam_io.h) - built
        // only when actually needed, so Node-mode runs never pay for it.
        const bool base_precision = (args.coverage_precision == CoveragePrecision::Base);
        std::vector<cdx::SeqLen> node_lengths;
        std::vector<BpGap> coverage_gaps;

        if (base_precision) {
            cfg::ScopedTimer timer("Node length lookup built");
            node_lengths = buildNodeLengthsQuery(data.nid2idx, data.idx2bp);
        }

        // ============================================================
        // STEP 2: GAM Processing
        // ============================================================
        std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] GAM Processing\n";
        GamMappingStats mapping_stats;

        {
            cfg::ScopedTimer timer("GAM alignments processed");

            mapping_stats = process_gam(
                args.gam_file,
                data.node_coverage,
                data.component.nid_min,
                GAM_BATCH_SIZE,
                args.decompression_threads,
                args.worker_threads,
                args.coverage_precision,
                base_precision ? &node_lengths : nullptr,
                base_precision ? &coverage_gaps : nullptr
            );

            timer.update_name("Processed " + std::to_string(mapping_stats.total) + " GAM reads");
        }

        // ============================================================
        // STEP 3: Coverage Projection
        // ============================================================
        std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Genomic Coordinate Projection\n";
        std::vector<cdx::Coverage> idx_coverage;

        {
            cfg::ScopedTimer timer("Node coverage projected to component index space");

            idx_coverage = projectCov2IdxQuery(
                data.node_coverage,
                data.nid2idx,
                data.idx2bp.size() - 1
            );
        }

        std::vector<cdx::Coverage> bp_coverage;

        {
            cfg::ScopedTimer timer("Index coverage expanded to genomic coordinates");
            bp_coverage = expandPosCovQuery(idx_coverage, data.idx2bp);
        }

        // Refine the uniform per-node fill above down to base-pair
        // precision, before the query-interval trim below - trimming only
        // masks positions outside the requested range with a sentinel, it
        // never changes indices, so applying gaps first or after would give
        // the same result inside the query interval; doing it first keeps
        // this step symmetric with the global pipeline (gaps are always
        // applied right after expansion, before any query/output-specific
        // step).
        if (base_precision) {
            cfg::ScopedTimer timer("Coverage gaps applied at base-pair precision");
            applyBpGapsQuery(bp_coverage, coverage_gaps, data.nid2idx, data.idx2bp);
        }

        {
            cfg::ScopedTimer timer("Coverage trimmed to requested query interval");
            bp_coverage = trimCoverageToQuery(bp_coverage, data.query_range_bp);
        }

        const std::filesystem::path output_directory(args.output_directory);

        // ============================================================
        // OPTIONAL OUTPUT: TSV
        // ============================================================
        if (args.generateTable()) {
            std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Coverage Table Export\n";

            {
                const std::filesystem::path output_path = output_directory / cfg::NAME_TSV_FILE;
                cfg::ScopedTimer timer("Coverage TSV table written");

                output::writeCoverageTsvQuery(
                    output_path,
                    bp_coverage,
                    data.component.compo_name
                );

                timer.update_name("Coverage table saved (" + output_path.filename().string() + ')');
            }
        }

        // ============================================================
        // OPTIONAL OUTPUT: Statistics
        // ============================================================
        if (args.generateStats()) {
            std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Statistics Report\n";

            {
                const std::filesystem::path output_path = output_directory / cfg::NAME_STATS_FILE;
                cfg::ScopedTimer timer("Coverage statistics report written");

                output::writeStatsReportQuery(
                    output_path,
                    mapping_stats,
                    bp_coverage,
                    data.component.compo_name,
                    args.coverage_precision
                );

                timer.update_name("Statistics report saved (" + output_path.filename().string() + ')');
            }
        }

        // ============================================================
        // OPTIONAL OUTPUT: Graph
        // ============================================================
        if (args.generateGraph()) {
            std::cerr << "\n[STEP " << current_step++ << '/' << total_steps << "] Coverage Graph Generation\n";

            cfg::ScopedTimer timer("Coverage graph");

            output::PlotConfig plot_config;

            plot_config.smoothing = static_cast<double>(args.smoothing);
            plot_config.max_plot_points = static_cast<std::size_t>(args.max_plot_points);
            plot_config.dpi = args.dpi;
            plot_config.figure_width = args.figure_width;
            plot_config.figure_height = args.figure_height;
            plot_config.line_color = args.line_color;
            plot_config.fill_color = args.fill_color;
            if (args.log_base.has_value()) {
                plot_config.log_base = args.log_base.value();
            }

            if (is_circular) {
                /*
                 * Contrairement au backend linéaire, on NE découpe PAS bp_coverage :
                 * le lissage circulaire (frontière position 0 / position compo_end
                 * voisines sur la molécule) et les requêtes traversant l'origine ont
                 * besoin du tableau complet de la composante. writeCircularPlotQuery
                 * gère lui-même l'extraction du parcours à partir de query_range_bp.
                 */
                output::writeCircularPlotQuery(
                    output_directory / cfg::NAME_GRAPH_FILE,
                    bp_coverage,
                    data.component.compo_name,
                    data.getComponentLength(),
                    data.query_range_bp,
                    plot_config
                );
            } else {
                /*
                 * bp_coverage couvre toute la composante (trimCoverageToQuery ne
                 * fait que masquer les positions hors requête, il ne redimensionne
                 * pas le vecteur). Pour que l'axe X du graphe corresponde
                 * exactement à la plage demandée (ex. -q "chr1 0:100" -> axe 0-100
                 * et non compo_start-compo_end), on découpe ici un sous-vecteur
                 * borné à [query_start, query_end] (voir query_plot_slice.h/.cpp).
                 */
                const std::vector<cdx::Coverage> plot_coverage =
                        output::sliceLinearQueryCoverage(bp_coverage, data.query_range_bp);

                output::writeLinearPlotQuery(
                    output_directory /cfg::NAME_GRAPH_FILE,
                    plot_coverage,
                    data.component.compo_name,
                    data.query_range_bp.first,
                    plot_config
                );
            }
        }
    }
} // namespace

/**
* @brief Entry point of the CDX coverage analysis application.
*
* Parses command-line arguments, initializes component-resolution metadata,
* and dispatches execution to the appropriate operating mode:
*
*   - Inspection mode: display component information from the CDX index.
*   - Query mode: analyze a single component or query interval.
*   - Global mode: analyze the entire pangenome.
*
* The function also handles output-directory creation, runtime reporting,
* component-name resolution, and top-level exception handling.
*
* @param argc Number of command-line arguments.
* @param argv Command-line argument array.
* @return EXIT_SUCCESS on successful completion, EXIT_FAILURE if an error occurs.
*/
int main(int argc, char **argv) {
    using Clock = std::chrono::steady_clock;
    const auto total_start = Clock::now();

    try {
        const CliArgs args = parse_args(argc, argv);

        // Inspection mode does not generate output files, so there is
        // no need to create the output directory.
        if (!args.inspectMode())
            std::filesystem::create_directories(args.output_directory);

        ComponentResolver resolver;

        // This metadata load is used only to build the mapping between component names and component IDs.
        {
            const cdx::GlobalData metadata = cdx::loadGlobal(args.cdx_file);

            for (std::size_t cid = 0; cid < metadata.layout.component_names.size(); ++cid) {
                resolver.register_component(cid, metadata.layout.component_names[cid]);
            }
        }

        // ========================================================
        // INSPECTION MODE
        // ========================================================
        if (args.inspectMode()) {
            std::optional<cdx::Cid> component_id = std::nullopt;

            if (args.inspect.component) {
                const ResolvedComponent resolved = resolver.resolve(*args.inspect.component);

                if (resolved.cid > static_cast<std::size_t>(std::numeric_limits<cdx::Cid>::max())) {
                    throw std::overflow_error("Resolved component ID exceeds cdx::Cid capacity.");
                }

                component_id = static_cast<cdx::Cid>(resolved.cid);
            }

            cdx::inspectComponent(
                std::filesystem::path(args.cdx_file),
                component_id
            );

            return EXIT_SUCCESS;
        }

        // Display program banner.
        std::cerr
                << std::string(TERMINAL_WIDTH, '=') << '\n'
                << std::setw(46) << "CDX COVERAGE" << '\n'
                << std::string(TERMINAL_WIDTH, '=') << '\n';


        // Dispatch execution according to the requested analysis mode.
        if (args.query) {

            // Resolve the user-specified component name or CID before running the query-level pipeline.
            const ResolvedComponent resolved = resolver.resolve(args.query->component);


            runQueryPipeline(args, resolved.cid, std::string(resolved.name));
        } else {
            runGlobalPipeline(args);
        }

        const double total_seconds = std::chrono::duration<double>(Clock::now() - total_start).count();
        std::cerr
                << '\n' << std::string(TERMINAL_WIDTH, '=') << '\n'
                << " [SUCCESS] Execution finished in " << std::fixed << std::setprecision(3) << total_seconds << " s\n"
                << std::string(TERMINAL_WIDTH, '=') << '\n';

        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        const double total_seconds = std::chrono::duration<double>(Clock::now() - total_start).count();
        std::cerr
                << '\n' << std::string(TERMINAL_WIDTH, '=') << '\n'
                << " [FATAL ERROR] " << error.what() << '\n'
                << " Execution failed after " << std::fixed << std::setprecision(3) << total_seconds << " s\n"
                << std::string(TERMINAL_WIDTH, '=') << '\n';

        return EXIT_FAILURE;
    } catch (...) {
        const double total_seconds = std::chrono::duration<double>(Clock::now() - total_start).count();
        std::cerr
                << '\n' << std::string(TERMINAL_WIDTH, '=') << '\n'
                << " [FATAL ERROR] Unknown non-standard error occurred.\n"
                << " Execution failed after " << std::fixed << std::setprecision(3) << total_seconds << " s\n"
                << std::string(TERMINAL_WIDTH, '=') << '\n';
        return EXIT_FAILURE;
    }
}
