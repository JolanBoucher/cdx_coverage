#include "cdx_loader.h"
#include "config.h"
#include "cov_projection.h"
#include "gam_io.h"
#include "output_coverage.h"
#include "output_stats.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

    class ScopedTimer {
    public:
        explicit ScopedTimer(std::string name)
            : name_(std::move(name)),
              start_(std::chrono::steady_clock::now()) {}

        ~ScopedTimer() {
            const auto end = std::chrono::steady_clock::now();

            const double seconds =
                std::chrono::duration<double>(
                    end - start_
                ).count();

            std::cout
                << "[TIME] "
                << name_
                << ": "
                << seconds
                << " s\n";
        }

    private:
        std::string name_;
        std::chrono::steady_clock::time_point start_;
    };

// ================================================================
// MANUAL CONFIGURATION
// ================================================================

const std::filesystem::path CDX_PATH = "./bio-file/messy.cdx";
const std::filesystem::path GAM_PATH = "./bio-file/messy.gam";
const std::filesystem::path OUTPUT_DIRECTORY = "./output";

// Set to std::nullopt for whole-graph mode.
// Set to a component ID, for example cdx::Cid{0}, for query mode.
constexpr std::optional<cdx::Cid> COMPONENT_ID = std::nullopt;

// Inclusive query coordinates used only in query mode.
// std::nullopt means the complete component.
constexpr std::optional<std::pair<std::int64_t, std::int64_t>> QUERY_RANGE = std::nullopt;

// Allows start > end for a circular origin-crossing query.
constexpr bool CIRCULAR = false;

constexpr std::size_t GAM_BATCH_SIZE = 2048;
constexpr int GAM_DECOMPRESSION_THREADS = 4;

const std::filesystem::path TSV_OUTPUT = OUTPUT_DIRECTORY / "tsv_coverage.tsv";
const std::filesystem::path STATS_OUTPUT = OUTPUT_DIRECTORY / "stats_coverage.txt";

using MappingStats = std::map<std::string, std::uint64_t>;

/**
 * @brief Copies absolute-node GAM coverage into a relative nid-space CDX table.
 *
 * Sentinel entries in the target table are preserved. Only active entries whose
 * value is below cfg::NOT_IN_QUERY receive coverage values.
 */
void transferGamCoverage(
    const std::vector<std::uint32_t>& absolute_coverage,
    const cdx::Nid nid_min,
    std::vector<cdx::Coverage>& target
) {
    for (std::size_t node_offset = 0;
         node_offset < target.size();
         ++node_offset) {
        if (target[node_offset] >= cfg::NOT_IN_QUERY) {
            continue;
        }

        const cdx::Nid node_id =
            nid_min + static_cast<cdx::Nid>(node_offset);

        if (node_id >=
            static_cast<cdx::Nid>(absolute_coverage.size())) {
            continue;
        }

        target[node_offset] =
            absolute_coverage[static_cast<std::size_t>(node_id)];
    }
}

/**
 * @brief Processes the GAM file and returns coverage indexed by absolute node ID.
 */
[[nodiscard]]
std::vector<std::uint32_t> processGam(
    const cdx::Nid maximum_node_id,
    MappingStats& mapping_stats
) {
    if (maximum_node_id >
        static_cast<cdx::Nid>(
            std::numeric_limits<std::size_t>::max()
        )) {
        throw std::overflow_error(
            "Maximum node ID exceeds size_t capacity."
        );
    }

    std::vector<std::uint32_t> gam_coverage;
    std::uint64_t read_count = 0;

    process_gam_fast(
        GAM_PATH.string(),
        gam_coverage,
        read_count,
        maximum_node_id,
        GAM_BATCH_SIZE,
        GAM_DECOMPRESSION_THREADS
    );

    // Temporary interpretation of process_gam_fast's read_count.
    // Adapt these fields if gam_io later exposes mapped/unmapped counters.
    mapping_stats["total"] = read_count;
    mapping_stats["mapped"] = read_count;
    mapping_stats["mapped_to_query"] = 0;
    mapping_stats["unmapped"] = 0;

    return gam_coverage;
}

void runGlobalPipeline(MappingStats& mapping_stats) {
    cdx::GlobalData data;
    {
        ScopedTimer timer("loadGlobal");
        data = cdx::loadGlobal(CDX_PATH);
    }

    std::vector<std::uint32_t> gam_coverage;
    {
        ScopedTimer timer("processGam");
        gam_coverage = processGam(data.layout.graph_nid_max, mapping_stats);
    }

    {
        ScopedTimer timer("transferGamCoverage");
        transferGamCoverage(
            gam_coverage,
            data.layout.graph_nid_min,
            data.node_coverage
        );
    }

    std::vector<cdx::Coverage> flat_idx_coverage;
    {
        ScopedTimer timer("projectCov2IdxGlobal");
        flat_idx_coverage = projectCov2IdxGlobal(
            data.node_coverage,
            data.nid2flat_idx,
            data.layout.component_offsets
        );
    }

    std::vector<cdx::Coverage> flat_bp_coverage;
    std::vector<cdx::PosBp> bp_component_offsets;
    {
        ScopedTimer timer("expandPosCovGlobal");
        std::tie(flat_bp_coverage, bp_component_offsets) = expandPosCovGlobal(
            flat_idx_coverage,
            data.layout.component_offsets,
            data.idx2bp,
            data.idx2bp_offsets
        );
    }

    {
        ScopedTimer timer("writeCoverageTsvGlobal");
        output::writeCoverageTsvGlobal(
            TSV_OUTPUT,
            flat_bp_coverage,
            bp_component_offsets,
            data.layout.component_names
        );
    }

    {
        ScopedTimer timer("writeStatsReportGlobal");
        output::writeStatsReportGlobal(
            STATS_OUTPUT,
            flat_bp_coverage,
            bp_component_offsets,
            data.layout.component_names,
            &mapping_stats
        );
    }
}

void runQueryPipeline(MappingStats& mapping_stats) {
    if (!COMPONENT_ID) {
        throw std::logic_error(
            "Query mode requires COMPONENT_ID."
        );
    }

    cdx::QueryData data;
    {
        ScopedTimer timer("loadQuery");
        data = cdx::loadQuery(
            CDX_PATH,
            *COMPONENT_ID,
            QUERY_RANGE,
            CIRCULAR
        );
    }

    std::vector<std::uint32_t> gam_coverage;
    {
        ScopedTimer timer("processGam");
        gam_coverage = processGam(data.component.nid_max, mapping_stats);
    }

    {
        ScopedTimer timer("transferGamCoverage");
        transferGamCoverage(
            gam_coverage,
            data.component.nid_min,
            data.node_coverage
        );
    }

    std::vector<cdx::Coverage> idx_coverage;
    {
        ScopedTimer timer("projectCov2IdxQuery");
        idx_coverage = projectCov2IdxQuery(
            data.node_coverage,
            data.nid2idx,
            data.idx2bp.size() - 1
        );
    }

    std::vector<cdx::Coverage> bp_coverage;
    {
        ScopedTimer timer("expandPosCovQuery");
        bp_coverage = expandPosCovQuery(
            idx_coverage,
            data.idx2bp
        );
    }

    {
        ScopedTimer timer("trimCoverageToQuery");
        bp_coverage = trimCoverageToQuery(
            bp_coverage,
            data.query_range_bp
        );
    }

    {
        ScopedTimer timer("writeCoverageTsvQuery");
        output::writeCoverageTsvQuery(
            TSV_OUTPUT,
            bp_coverage,
            data.component.compo_name
        );
    }

    {
        ScopedTimer timer("writeStatsReportQuery");
        output::writeStatsReportQuery(
            STATS_OUTPUT,
            mapping_stats,
            bp_coverage,
            data.component.compo_name
        );
    }
}

} // anonymous namespace

int main() {
    try {
        ScopedTimer total_timer("TOTAL PROGRAM");

        std::filesystem::create_directories(
            OUTPUT_DIRECTORY
        );

        if (!std::filesystem::is_directory(
                OUTPUT_DIRECTORY)) {
            throw std::runtime_error(
                "Output path is not a directory: " +
                OUTPUT_DIRECTORY.string()
            );
        }

        MappingStats mapping_stats{
            {"total", 0},
            {"mapped", 0},
            {"mapped_to_query", 0},
            {"unmapped", 0}
        };

        if (COMPONENT_ID) {
            runQueryPipeline(mapping_stats);
        }
        else {
            runGlobalPipeline(mapping_stats);
        }

        std::cout
            << "Coverage TSV written to: "
            << TSV_OUTPUT
            << '\n'
            << "Coverage statistics written to: "
            << STATS_OUTPUT
            << '\n';

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "Error: "
            << error.what()
            << '\n';

        return 1;
    }
}