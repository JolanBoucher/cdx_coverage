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
    [[nodiscard]] QueryData loadQuery(
        const std::filesystem::path &cdx_path,
        std::size_t component_id,
        std::optional<std::pair<int64_t, int64_t> > query_range = std::nullopt,
        bool circular = false
    );

    [[nodiscard]] GlobalData loadGlobal(
        const std::filesystem::path &cdx_path
    );

    void inspectComponent(
        const std::filesystem::path &cdx_path,
        std::optional<Cid> component_id = std::nullopt
    );
}

#endif //CDX_COVERAGE_CDX_LOADER_H
