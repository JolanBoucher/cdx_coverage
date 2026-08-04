/**
 * @file cdx_loader.h
 * @brief High-level CDX loading API.
 *
 * Provides:
 *  - component loading
 *  - query initialization
 *  - graph-wide loading
 *  - component inspection
 */

#ifndef CDX_COVERAGE_CDX_LOADER_H
#define CDX_COVERAGE_CDX_LOADER_H

#include "cdx_types.h"
#include <filesystem>
#include <optional>

namespace cdx {
    /**
     * @brief Load, build, and initialize all data structures required for a query into a consolidated structure.
     *
     * Opens the binary CDX file, retrieves component boundaries, parses records, resolves query ranges,
     * builds topological position-to-index translation maps, and populates the dense node coverage table.
     *
     * @param cdx_path Filesystem path to the binary CDX index file.
     * @param component_id Numeric identifier of the component to load.
     * @param query_range Optional inclusive (start_bp, end_bp) coordinate pair. std::nullopt defaults to full component.
     * @param circular Boolean flag indicating whether the component supports circular coordinate wrapping.
     * @return cdx::QueryData Consolidated structure containing initialized query maps and metadata.
     */
    [[nodiscard]] QueryData loadQuery(
        const std::filesystem::path &cdx_path,
        std::size_t component_id,
        std::optional<std::pair<int64_t, int64_t> > query_range = std::nullopt,
        bool circular = false
    );


    /**
     * @brief Loads and initializes global pangenome graph metadata, coordinate lookup tables,
     *        and dense coverage structures from a binary CDX archive.
     *
     * Performs stream validation, parses total layout bounds, builds global node-ID to flat index mappings,
     * constructs position prefix-sum arrays, and derives per-component base-pair lengths.
     *
     * @param cdx_path Filesystem path to the target uncompressed binary CDX file.
     * @return GlobalData Fully populated global lookup tables and graph layout metadata.
     *
     * @throws std::runtime_error If file cannot be opened or underlying reader passes fail.
     */
    [[nodiscard]] GlobalData loadGlobal(
        const std::filesystem::path &cdx_path
    );


    /**
     * @brief Inspects CDX file structure and prints a formatted summary table to stdout.
     *
     * Depending on whether `component_id` is specified, prints metadata (Component ID, name, sequence length,
     * node count, node ID range) for either a single component or all components in the archive.
     *
     * @param cdx_path Filesystem path to the binary CDX index file.
     * @param component_id Optional specific component ID to inspect. std::nullopt inspects all components.
     *
     * @throws std::runtime_error If the file cannot be opened or contains no components.
     * @throws std::out_of_range If the requested component ID is invalid.
     */
    void inspectComponent(
        const std::filesystem::path &cdx_path,
        std::optional<Cid> component_id = std::nullopt
    );
}

#endif //CDX_COVERAGE_CDX_LOADER_H
