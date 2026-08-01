//
// Created by Jolan on 2026-07-31.
//

#include "cdx_loader.h"

#include "cdx_types.h"
#include "cdx_format.h"
#include "config.h"

#include <fstream>
#include <stdexcept>

namespace {
    // this namespace for the private function used by loadQuery

    std::vector<cdx::PosBp> buildIdx2Pos(
        const std::vector<cdx::NodeRecord> &records,
        const std::size_t node_count
    ) {
        // initialize the output vector
        std::vector<cdx::PosBp> idx2bp(node_count + 1, 0);

        // Placing the length in topological order
        for (const auto &record: records) {
#ifndef NDEBUG
            if (record.idx >= node_count) {
                throw std::runtime_error("Node record for " + std::to_string(record.node_id) + " out of bounds");
            }
#endif
            idx2bp[record.idx + 1] = record.seq_len;
        }

        // transform the topological order in cumulative position table
        for (std::size_t idx = 1; idx < idx2bp.size(); ++idx){
            idx2bp[idx] += idx2bp[idx - 1];
        }
    };

    /**
     * @brief Builds a rank/select-like position index to map genomic base pairs (bp) to node indices in O(1) time.
     *
     * Takes a sorted vector of base pair coordinates representing node boundaries and constructs
     * a bitvector marking node starts, along with a cumulative rank lookup table.
     *
     * @param idx2bp Sorted vector containing the starting base pair position for each node.
     * @return cdx::PositionIndex Structure containing the compressed bitvector and precomputed rank indices.
     */
    [[nodiscard]]
    cdx::PositionIndex buildPosRankIndex(
        const std::vector<cdx::PosBp> &idx2bp
    ) {
        constexpr std::uint64_t WORD_SHIFT = 6;
        constexpr std::uint64_t WORD_MASK = 63;

        cdx::PositionIndex index;
        if (idx2bp.empty()) {
            return index;
        }

        // --- Debug safeguard: Ensure input base pair coordinates are monotonically sorted ---
#ifndef NDEBUG
        for (std::size_t idx = 1; idx < idx2bp.size(); ++idx) {
            if (idx2bp[idx] < idx2bp[idx - 1]) {
                throw std::runtime_error("idx2bp must be sorted");
            }
        }
#endif

        const cdx::PosBp component_length_bp = idx2bp.back();
        const std::size_t node_count = idx2bp.size() - 1;

        // Calculate total 64-bit words required to cover the component length
        const std::size_t word_count = (component_length_bp + WORD_MASK) >> WORD_SHIFT;
        index.bitvector.assign(word_count, 0);

        // --- Mark every node start bit inside the bitvector ---
        for (std::size_t idx = 0; idx < node_count; ++idx) {
            const cdx::PosBp start_bp = idx2bp[idx];
            const std::size_t word_idx = start_bp >> WORD_SHIFT;
            const std::size_t bit_idx = start_bp & WORD_MASK;
            index.bitvector[word_idx] |= 1ULL << bit_idx;
        }

        // --- Build the cumulative rank index table ---
        index.rank_index.resize(word_count + 1, 0);

        std::uint32_t cumulative_rank = 0;

        for (std::size_t word_idx = 0; word_idx < word_count; ++word_idx) {
            index.rank_index[word_idx] = cumulative_rank;

            // count the number of 1 bit in a word
            cumulative_rank += __builtin_popcountll(index.bitvector[word_idx]);
        }
        index.rank_index[word_count] = cumulative_rank;
        return index;
    }

    /**
     * @brief Builds a user-friendly origin-crossing query error message.
     *
     * @param component_id Identifier of the graph component being queried.
     * @param original_start Original signed start coordinate supplied by the user.
     * @param original_end Original signed end coordinate supplied by the user.
     * @param start Resolved absolute start base pair position.
     * @param end Resolved absolute end base pair position.
     * @return std::string Formatted error message string.
     */
    std::string formatOriginCrossingError(
        const std::size_t component_id,
        const cdx::PosBp original_start,
        const cdx::PosBp original_end,
        const cdx::PosBp start,
        const cdx::PosBp end
    ) {
        const cdx::PosBp linear_start = std::min(start, end);
        const cdx::PosBp linear_end = std::max(start, end);

        return (
            "Query " + std::to_string(component_id) + ":" + std::to_string(original_start) + "-" +
            std::to_string(original_end) + " resolves to " + std::to_string(start) + "-" + std::to_string(end)
            + " and crosses the component origin.\n\n" +

            "If this component is circular, rerun with:\n  --component-type circular\n\n" +

            "If you intended the linear interval between both positions, use:\n" +
            "  " + std::to_string(component_id) + ":" + std::to_string(linear_start) + "-" + std::to_string(linear_end)
        );
    }


    /**
     * @brief Validates, normalizes, and resolves user query coordinate ranges.
     *
     * Handles negative indexing by wrapping coordinates relative to the component length,
     * validates absolute bounds, and verifies origin-crossing rules against the component's
     * circularity status.
     *
     * @param component_id Identifier of the graph component being queried.
     * @param component_length Total length of the component in base pairs.
     * @param query_range Optional pair of signed user coordinates (start, end).
     * @param circular Boolean flag indicating whether the component supports circular coordinate wrapping.
     * @return std::pair<cdx::PosBp, cdx::PosBp> The validated and normalized absolute base pair boundaries.
     */
    [[nodiscard]] std::pair<cdx::PosBp, cdx::PosBp> resolveQueryRange(
        const std::size_t component_id,
        const cdx::PosBp component_length,
        const std::optional<std::pair<std::int64_t, std::int64_t> > &query_range,
        const bool circular
    ) {
        // Guard against division-by-zero or malformed components
        if (component_length == 0) {
            throw std::invalid_argument("Component " + std::to_string(component_id) + " has invalid length 0.");
        }

        // Default to the full span of the component if no query range is provided
        if (!query_range) {
            return {0, component_length - 1};
        }

        const auto [original_start, original_end] = *query_range;

        // Normalize negative coordinates by wrapping them around the component length (Option B behavior)
        const std::int64_t resolved_start = (original_start < 0)
                                                ? original_start + static_cast<std::int64_t>(component_length)
                                                : original_start;

        const std::int64_t resolved_end = (original_end < 0)
                                              ? original_end + static_cast<std::int64_t>(component_length)
                                              : original_end;

        // Validate that resolved boundaries fit strictly within component coordinate space
        if (resolved_start < 0 || resolved_start >= static_cast<std::int64_t>(component_length)) {
            throw std::out_of_range("Query start coordinate out of bounds.");
        }

        if (resolved_end < 0 || resolved_end >= static_cast<std::int64_t>(component_length)) {
            throw std::out_of_range("Query end coordinate out of bounds.");
        }

        const auto query_start = static_cast<cdx::PosBp>(resolved_start);
        const auto query_end = static_cast<cdx::PosBp>(resolved_end);

        // Check if the coordinate range wraps around the graph origin (start > end)
        const bool crosses_origin = query_start > query_end;

        if (crosses_origin && !circular) {
            throw std::invalid_argument(formatOriginCrossingError(
                component_id,
                original_start,
                original_end,
                query_start,
                query_end));
        }
        return {
            query_start, query_end
        };
    }

    /**
     * @brief Projects a genomic coordinate to the node containing that position.
     *
     * Uses precomputed rank tables and a bitmask combined with compiler intrinsics (__builtin_popcountll)
     * to count active node-start bits up to the target coordinate position in O(1) time.
     *
     * @param bp Absolute base pair position to project.
     * @param index Reference to the PositionIndex containing the bitvector and rank table.
     * @return cdx::Idx The mapped node index containing the target base pair.
     */
    [[nodiscard]] cdx::Idx projectBpToIdx(
        const cdx::PosBp bp,
        const cdx::PositionIndex &index
    ) {
        constexpr std::uint64_t WORD_SHIFT = 6;
        constexpr std::uint64_t WORD_MASK = 63;

        // Locate the 64-bit word index and specific bit offset for the target base pair
        const std::size_t word_idx = bp >> WORD_SHIFT;
        const std::size_t bit_idx = bp & WORD_MASK;

        // Guard against undefined behavior (1ULL << 64) when bit_idx reaches 63
        const std::uint64_t mask = (bit_idx == 63) ? UINT64_MAX : ((1ULL << (bit_idx + 1)) - 1ULL);

        // Rank equals the number of node starts up to and including bp
        const std::uint32_t rank = index.rank_index[word_idx] +
                                   __builtin_popcountll(index.bitvector[word_idx] & mask);

        // Convert rank to 0-based node index
        return rank - 1;
    }

    /**
     * @brief Projects a genomic query range (start and end base pairs) into corresponding node indices.
     *
     * Eliminates code duplication by delegating single coordinate mapping to projectBpToIdx,
     * cleanly transforming bp intervals into structural index boundaries.
     *
     * @param query_range_bp Pair of absolute base pair coordinates (start, end).
     * @param index Reference to the PositionIndex containing the bitvector and rank lookup structures.
     * @return std::pair<cdx::Idx, cdx::Idx> The mapped start and end node indices.
     */
    [[nodiscard]] std::pair<cdx::Idx, cdx::Idx> projectPos2IdxQuery(
        const std::pair<cdx::PosBp, cdx::PosBp> &query_range_bp,
        const cdx::PositionIndex &index
    ) {
        return {
            projectBpToIdx(query_range_bp.first, index),
            projectBpToIdx(query_range_bp.second, index)
        };
    }

    /**
     * @brief Builds a dense node-id to topological-index lookup table.
     *
     * The returned vector is indexed by  node_id - nid_min  allowing constant-time translation:
     * local_idx = `nid2idx[node_id - nid_min]`;
     *
     * Entries corresponding to node identifiers absent from the component
     * are initialized with cfg::NOT_IN_COMPO.
     *
     * @param records Component node records.
     * @param nid_min Minimum node identifier in the component.
     * @param nid_max Maximum node identifier in the component.
     *
     * @return A dense nid lookup table.
     */
    [[nodiscard]] std::vector<cdx::Idx> buildNid2IdxQuery(
        const std::vector<cdx::NodeRecord> &records,
        const cdx::Nid nid_min,
        const cdx::Nid nid_max
    ) {
        // Compute total dense table size covering the node ID span
        const std::size_t component_size = nid_max - nid_min + 1;

        // Initialize the lookup table with the default 'not in component' sentinel value
        std::vector<cdx::Idx> nid2idx(component_size, cfg::NOT_IN_COMPO);

        // Populate the dense array using direct offset indexing for O(1) lookups
        for (const auto &record: records) {
            // --- Debug safeguard: Ensure record node IDs fall strictly within bounds ---
#ifndef NDEBUG
            if (record.node_id < nid_min || record.node_id > nid_max) {
                throw std::runtime_error("NodeRecord node_id outside component range");
            }
#endif

            const std::size_t node_offset = record.node_id - nid_min;
            nid2idx[node_offset] = record.idx;
        }
        return nid2idx;
    }

    /**
      * @brief Initialize and populate the dense node coverage table for a target query.
      *
      * Pre-allocates a vector indexed by relative node offset (`node_id - nid_min`)
      * initialized to `cfg::NOT_IN_COMPO`. As records are processed, each node's topological
      * index is evaluated against the query boundaries:
      *   - Active query nodes are assigned `0` (ready for coverage accumulation).
      *   - In-component but out-of-query nodes are assigned `cfg::NOT_IN_QUERY`.
      *
      * @param records Component node records.
      * @param nid_min Minimum node ID present in the component.
      * @param nid_max Maximum node ID present in the component.
      * @param query_range_idx Inclusive `(query_start, query_end)` topological index bounds.
      * @return std::vector<cdx::Coverage> Dense coverage table vector.
      */
    [[nodiscard]] std::vector<cdx::Coverage> buildCovTableQuery(
        const std::vector<cdx::NodeRecord> &records,
        const cdx::Nid nid_min,
        const cdx::Nid nid_max,
        const std::pair<cdx::Idx, cdx::Idx> &query_range_idx
    ) {
        // Compute total dense table size covering the node ID span
        const std::size_t component_size = nid_max - nid_min + 1;

        // Initialize the lookup table with the default 'not in component' sentinel value
        std::vector<cdx::Coverage> cov_table(component_size, cfg::NOT_IN_COMPO);

        const auto [query_start, query_end] = query_range_idx;

        // Detect whether the query wraps around a circular component boundary
        const bool crosses_origin = query_start > query_end;

        // Process all node records from the component chunk
        for (const auto &record: records) {
            // --- Debug safeguard: Ensure record node IDs fall strictly within bounds ---
#ifndef NDEBUG
            if (record.node_id < nid_min || record.node_id > nid_max) {
                throw std::runtime_error("NodeRecord node_id outside component range during coverage table build");
            }
#endif

            const std::size_t node_offset = record.node_id - nid_min;

            // Evaluate topological index membership based on query topology (handling circular wrap-around if needed)
            const bool in_query = crosses_origin
                                      ? (record.idx >= query_start || record.idx <= query_end)
                                      : (record.idx >= query_start && record.idx <= query_end);

            // Mark active query nodes with 0, and out-of-query component nodes with sentinels
            cov_table[node_offset] = in_query ? 0 : cfg::NOT_IN_QUERY;
        }
        return cov_table;
    }
} // anonymous namespace

namespace cdx {
    /**
 * @brief Load, build, and initialize all data structures required for a query into a consolidated structure.
 *
 * Opens the binary CDX file, retrieves component boundaries, parses records, resolves query ranges,
 * builds topological position-to-index translation maps, and populates the dense node coverage table.
 *
 * @param cdx_path Filesystem path to the binary CDX index file.
 * @param component_id Numeric identifier of the component to load.
 * @param query_range Optional inclusive (start_bp, end_bp) coordinate pair. std::nullopt defaults to the full component.
 * @param circular Boolean flag indicating whether the component supports circular coordinate wrapping.
 * @return cdx::QueryData Consolidated structure containing all initialized query maps and metadata.
 */
    [[nodiscard]] QueryData loadQuery(
        const std::filesystem::path &cdx_path,
        std::size_t component_id,
        std::optional<std::pair<std::int64_t, std::int64_t> > query_range,
        bool circular
    ) {
        QueryData query_data;

        // 1. Open binary CDX index stream
        std::ifstream cdx_stream(cdx_path, std::ios::binary);

        if (!cdx_stream) {
            throw std::runtime_error("Unable to open CDX file: " + cdx_path.string());
        }

        // 2. Locate component boundaries and metadata
        query_data.component = cdx::seekComponent(cdx_stream, component_id);

        // 3. Read node records payload from the binary stream
        std::vector<cdx::NodeRecord> records = cdx::readComponentPayload(
            cdx_stream,
            query_data.component.nb_nodes
        );

        // 4. Build idx-space to genomic position-space prefix sum array
        query_data.idx2bp = buildIdx2Pos(records, query_data.component.nb_nodes);

        if (query_data.idx2bp.empty()) {
            throw std::runtime_error("Component contains no nodes.");
        }

        const cdx::PosBp component_length = query_data.getComponentLength();

        // 5. Resolve and project user query boundaries or default to full component span
        if (query_range) {
            query_data.query_range_bp = resolveQueryRange(
                component_id,
                component_length,
                query_range,
                circular
            );

            query_data.position_index = buildPosRankIndex(query_data.idx2bp);

            query_data.query_range_idx = projectPos2IdxQuery(
                query_data.query_range_bp,
                query_data.position_index
            );
        } else {
            query_data.query_range_bp = {0, component_length - 1};
            query_data.query_range_idx = {0, static_cast<Idx>(query_data.component.nb_nodes - 1)};
        }

        // 6. Build dense node-id to topological-index translation table
        query_data.nid2idx = buildNid2IdxQuery(
            records,
            query_data.component.nid_min,
            query_data.component.nid_max
        );

        // 7. Initialize and activate the dense node coverage table
        query_data.node_coverage = buildCovTableQuery(
            records,
            query_data.component.nid_min,
            query_data.component.nid_max,
            query_data.query_range_idx
        );

        return query_data;
    }
}

//TODO ajouté et complété cdx_io à cdx_lib (anciennement shared_cdx)