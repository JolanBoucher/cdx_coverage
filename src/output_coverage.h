#ifndef CDX_COVERAGE_OUTPUT_COVERAGE_H
#define CDX_COVERAGE_OUTPUT_COVERAGE_H

#include "cdx_types.h"

#include <filesystem>
#include <string>
#include <vector>

namespace output {

    void writeCoverageTsvQuery(
        const std::filesystem::path& output_tsv,
        const std::vector<cdx::Coverage>& bp_cov_table,
        const std::string& component_name
    );

    void writeCoverageTsvGlobal(
        const std::filesystem::path& output_tsv,
        const std::vector<cdx::Coverage>& flat_bp_cov_table,
        const std::vector<cdx::PosBp>& bp_component_offsets,
        const std::vector<std::string>& component_names
    );

} // namespace output

#endif // CDX_COVERAGE_OUTPUT_COVERAGE_H