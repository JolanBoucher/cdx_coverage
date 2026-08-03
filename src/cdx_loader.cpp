//
// Created by Jolan on 2026-07-31.
//

#include "cdx_loader.h"

#include "cdx_format.h"
#include "cdx_IO.h"
#include "config.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace {
    // this namespace for the private function used by loadQuery

    /**
     * @brief Constructs the cumulative base-pair position prefix sum array from node records.
     *
     * @param records Collection of node records containing local topological indices and sequence lengths.
     * @param node_count Total number of nodes in the component.
     * @return std::vector<cdx::PosBp> Array of size (node_count + 1) mapping local indices to cumulative base-pair positions.
     */
    [[nodiscard]] std::vector<cdx::PosBp> buildIdx2Pos(
        const std::vector<cdx::NodeRecord> &records,
        const std::size_t node_count
    ) {
        // Allocate node_count + 1 entries initialized to zero (idx2bp[0] serves as origin 0)
        std::vector<cdx::PosBp> idx2bp(node_count + 1, 0);

        // Populate element sequence lengths offset by 1 to reserve index 0
        for (const cdx::NodeRecord &record: records) {
#ifndef NDEBUG
            if (record.idx >= node_count) {
                throw std::runtime_error(
                    "Node record for node " + std::to_string(record.node_id) +
                    " has an out-of-bounds local index " + std::to_string(record.idx) + "."
                );
            }
#endif

            idx2bp[static_cast<std::size_t>(record.idx) + 1] = record.seq_len;
        }
        // Convert individual sequence lengths into prefix sums (e.g., [10, 25, 15] -> [0, 10, 35, 50])
        for (std::size_t idx = 1; idx < idx2bp.size(); ++idx) {
            idx2bp[idx] += idx2bp[idx - 1];
        }
        return idx2bp;
    }

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

        // --- Mark every node start a bit inside the bitvector ---
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
     * @param resolved_start Resolved absolute start base pair position.
     * @param resolved_end Resolved absolute end base pair position.
     * @return std::string Formatted error message string.
     */
    std::string formatOriginCrossingError(
        const std::size_t component_id,
        const std::int64_t original_start,
        const std::int64_t original_end,
        const cdx::PosBp resolved_start,
        const cdx::PosBp resolved_end
    ) {
        const cdx::PosBp linear_start = std::min(resolved_start, resolved_end);
        const cdx::PosBp linear_end = std::max(resolved_start, resolved_end);

        return "Query " + std::to_string(component_id) + ":" + std::to_string(original_start) + "-" +
               std::to_string(original_end) + " resolves to " + std::to_string(resolved_start) + "-" + std::to_string(
                   resolved_end)
               + " and crosses the component origin.\n\n" +

               "If this component is circular, rerun with:\n  --component-type circular\n\n" +

               "If you intended the linear interval between both positions, use:\n" +
               "  " + std::to_string(component_id) + ":" + std::to_string(linear_start) + "-" + std::to_string(
                   linear_end);
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
        const std::int64_t resolved_start = original_start < 0
                                                ? original_start + static_cast<std::int64_t>(component_length)
                                                : original_start;

        const std::int64_t resolved_end = original_end < 0
                                              ? original_end + static_cast<std::int64_t>(component_length)
                                              : original_end;

        // Validate that both resolved coordinates lie inside the component.
        if (resolved_start < 0 || resolved_start >= static_cast<std::int64_t>(component_length)) {
            throw std::out_of_range("Query start coordinate out of bounds.");
        }

        if (resolved_end < 0 || resolved_end >= static_cast<std::int64_t>(component_length)) {
            throw std::out_of_range("Query end coordinate out of bounds.");
        }

        // A start greater than the end represents an interval wrapping
        // around the origin of a circular component.
        const bool crosses_origin = resolved_start > resolved_end;

        if (crosses_origin && !circular) {
            throw std::invalid_argument(
                formatOriginCrossingError(
                    component_id,
                    original_start,
                    original_end,
                    static_cast<cdx::PosBp>(resolved_start),
                    static_cast<cdx::PosBp>(resolved_end)
                )
            );
        }

        // The validation above guarantees that these signed values are
        // non-negative and safely representable as PosBp.
        return {
            static_cast<cdx::PosBp>(resolved_start),
            static_cast<cdx::PosBp>(resolved_end)
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
        const std::uint64_t mask = bit_idx == 63 ? UINT64_MAX : (1ULL << (bit_idx + 1)) - 1ULL;

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
        std::vector nid2idx(component_size, cfg::NOT_IN_COMPO);

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
        std::vector cov_table(component_size, cfg::NOT_IN_COMPO);

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
                                      ? record.idx >= query_start || record.idx <= query_end
                                      : record.idx >= query_start && record.idx <= query_end;

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
        query_data.component = seekComponent(cdx_stream, component_id);

        // 3. Read node records payload from the binary stream
        std::vector<NodeRecord> records = readComponentPayload(
            cdx_stream,
            query_data.component.nb_nodes
        );

        // 4. Build idx-space to genomic position-space prefix sum array
        query_data.idx2bp = buildIdx2Pos(records, query_data.component.nb_nodes);

        if (query_data.idx2bp.empty()) {
            throw std::runtime_error("Component contains no nodes.");
        }

        const PosBp component_length = query_data.getComponentLength();

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

namespace {
    // this namespace used for private functions of loadGlobal

    [[nodiscard]]
    cdx::GraphLayout computeGraphLayout(
        std::istream &cdx
    ) {
        cdx.clear();
        cdx.seekg(0, std::ios::beg);

        if (!cdx) {
            throw std::runtime_error("Unable to seek to beginning of CDX stream.");
        }

        const cdx::FileHeader file_header = cdx::readGlobalHeader(cdx);
        cdx::GraphLayout layout;

        layout.component_count = static_cast<std::size_t>(file_header.n_components);
        layout.component_offsets.reserve(layout.component_count + 1);
        layout.component_names.reserve(layout.component_count);
        layout.component_offsets.push_back(0);
        bool first_component = true;

        for (cdx::Cid cid = 0; cid < file_header.n_components; ++cid) {
            const cdx::ComponentInfo component = cdx::readComponentHeader(cdx, cid);
            if (component.nb_nodes == 0) {
                throw std::runtime_error("Component " + std::to_string(cid) + " contains zero records.");
            }

            // Preserve the component name returned by readComponentHeader().
            layout.component_names.push_back(component.compo_name);

            if (first_component) {
                layout.graph_nid_min = component.nid_min;
                layout.graph_nid_max = component.nid_max;
                first_component = false;
            } else {
                layout.graph_nid_min = std::min(layout.graph_nid_min, component.nid_min);
                layout.graph_nid_max = std::max(layout.graph_nid_max, component.nid_max);
            }

            layout.total_nodes += component.nb_nodes;
            layout.component_offsets.push_back(layout.total_nodes);

            /*
             * readComponentHeader() has consumed the fixed header and
             * the variable-length component name. The stream is now at
             * the beginning of the NodeRecord payload.
             */
            cdx.seekg(component.payload_size, std::ios::cur);
            if (!cdx) {
                throw std::runtime_error("Unable to skip payload of component " + std::to_string(cid));
            }
        }

        return layout;
    }

    [[nodiscard]] std::vector<cdx::FlatIdx> buildNid2FlatIdxGlobal(
        std::istream &cdx,
        const cdx::Nid graph_nid_min,
        const cdx::Nid graph_nid_max,
        const std::vector<cdx::RecordCount> &component_offsets
    ) {
        const std::size_t component_count = component_offsets.size() - 1;
        const std::size_t graph_nid_range = graph_nid_max - graph_nid_min + 1;

        // Initialize a dense nid-space lookup
        std::vector nid2flat_idx(graph_nid_range, cfg::INVALID_FLAT_IDX);

        // Reset stream
        cdx.clear();
        cdx.seekg(0, std::ios::beg);

        if (!cdx) {
            throw std::runtime_error("Unable to seek to beginning of CDX stream.");
        }

        cdx::readGlobalHeader(cdx); // we just use to place the cursor after the global header

        // Process all components sequentially
        for (cdx::Cid component_id = 0; component_id < component_count; ++component_id) {
            const cdx::ComponentInfo component = cdx::readComponentHeader(cdx, component_id);
            const auto records = cdx::readComponentPayload(cdx, component.nb_nodes);
            const auto component_offset = static_cast<cdx::FlatIdx>(component_offsets[component_id]);

            for (const auto &record: records) {
                const std::size_t node_offset = record.node_id - graph_nid_min;
#ifndef NDEBUG
                if (record.node_id < graph_nid_min || record.node_id > graph_nid_max) {
                    throw std::runtime_error("Node ID outside graph bounds.");
                }
#endif
                nid2flat_idx[node_offset] = component_offset + static_cast<cdx::FlatIdx>(record.idx);
            }
        }

        return nid2flat_idx;
    }

    [[nodiscard]] std::vector<cdx::Coverage> buildCovTableGlobal(
        const std::vector<cdx::FlatIdx> &nid2flat_idx
    ) {
        // initialize global cov_table
        std::vector cov_table(nid2flat_idx.size(), cfg::NOT_IN_COMPO);
        for (std::size_t node_offset = 0; node_offset < nid2flat_idx.size(); ++node_offset) {
            // initialize the vector valid slots
            if (nid2flat_idx[node_offset] != cfg::INVALID_FLAT_IDX) {
                cov_table[node_offset] = 0;
            }
        }
        return cov_table;
    }


    /**
     * @brief Builds a global concatenated idx2pos map across all components in a CDX stream.
     *
     * Resets the input stream, verifies the global header, and constructs contiguous local position
     * prefix-sum tables while recording component start/end boundary indices inside the array.
     *
     * @param input Open seekable binary CDX stream.
     * @param component_count Total expected number of components to process.
     *
     * @return std::pair<std::vector<PosBp>, std::vector<Idx>>
     *         - first:  Concatenated prefix-sum arrays across all processed components.
     *         - second: Boundary element indices in idx2pos for each component (size component_count + 1).
     *
     * @throws std::runtime_error If stream seeking or payload reading fails.
     * @throws std::out_of_range If requested component_count exceeds archive contents.
     */
    [[nodiscard]] std::pair<std::vector<cdx::PosBp>, std::vector<cdx::RecordCount> > buildIdx2PosGlobal(
        std::istream &input,
        cdx::Cid component_count
    ) {
        // 1. Reset stream to origin and validate global header
        input.clear();
        input.seekg(0, std::ios::beg);

        if (!input) {
            throw std::runtime_error("Unable to seek to beginning of CDX stream.");
        }

        const cdx::FileHeader global_header = cdx::readGlobalHeader(input);
        if (component_count > global_header.n_components) {
            throw std::out_of_range("Requested component count (" + std::to_string(component_count) +
                                    ") exceeds archive capacity (" + std::to_string(global_header.n_components) + ").");
        }

        std::vector<cdx::PosBp> idx2pos;
        std::vector<cdx::RecordCount> idx2pos_offsets;
        idx2pos_offsets.reserve(static_cast<std::size_t>(component_count) + 1);
        idx2pos_offsets.push_back(0);

        // 2. Pass 1: compute exact element count to pre-allocate memory
        const std::streampos payload_start_pos = input.tellg();
        std::size_t total_elements = 0;

        for (cdx::Cid id = 0; id < component_count; ++id) {
            const cdx::ComponentInfo comp_info = cdx::readComponentHeader(input, id);
            total_elements += static_cast<std::size_t>(comp_info.nb_nodes) + 1;
            input.seekg(comp_info.payload_size, std::ios::cur);
        }

        idx2pos.reserve(total_elements);

        // 3. Pass 2: rewind to headers and populate payload data
        input.clear();
        input.seekg(payload_start_pos, std::ios::beg);

        for (cdx::Cid id = 0; id < component_count; ++id) {
            const cdx::ComponentInfo component = cdx::readComponentHeader(input, id);
            std::vector<cdx::NodeRecord> records;
            records = cdx::readComponentPayload(input, component.nb_nodes);

            const std::vector<cdx::PosBp> comp_idx2pos = buildIdx2Pos(records, component.nb_nodes);

            // Append component-local prefix sums into the contiguous vector
            idx2pos.insert(idx2pos.end(), comp_idx2pos.begin(), comp_idx2pos.end());

            // Store exclusive end boundary in index space
            idx2pos_offsets.push_back(static_cast<cdx::Idx>(idx2pos.size()));
        }

        return {std::move(idx2pos), std::move(idx2pos_offsets)};
    }
} // end private function supporting loadGlobal()

namespace cdx {
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
     * @throws std::runtime_error If the file cannot be opened or if underlying reader passes fail.
     */
    [[nodiscard]] GlobalData loadGlobal(const std::filesystem::path &cdx_path) {
        GlobalData global_data;

        // 1. Open the binary CDX stream.
        std::ifstream cdx_stream(cdx_path, std::ios::binary);
        if (!cdx_stream) {
            throw std::runtime_error("Unable to open CDX file: " + cdx_path.string());
        }

        // 2. Extract graph-wide metadata, component boundaries, and component names.
        GraphLayout layout = computeGraphLayout(cdx_stream);

        // 3. Build the dense node-ID to flat-index translation table.
        global_data.nid2flat_idx = buildNid2FlatIdxGlobal(
            cdx_stream,
            layout.graph_nid_min,
            layout.graph_nid_max,
            layout.component_offsets
        );

        // 4. Initialize the global coverage table in relative nid-space.
        global_data.node_coverage = buildCovTableGlobal(global_data.nid2flat_idx);

        // 5. Build concatenated component-local idx-to-bp tables and store them directly.
        auto [idx2bp, idx2bp_offsets] = buildIdx2PosGlobal(
            cdx_stream,
            static_cast<Cid>(layout.component_count)
        );

        global_data.idx2bp = std::move(idx2bp);
        global_data.idx2bp_offsets = std::move(idx2bp_offsets);

        // Validation post-move
        if (global_data.idx2bp_offsets.size() != layout.component_count + 1) {
            throw std::runtime_error("idx2bp offset count does not match the number of CDX components.");
        }

        // 6. Extract each component's total bp length using global_data fields.
        global_data.component_lengths.reserve(layout.component_count);

        for (std::size_t cid = 0; cid < layout.component_count; ++cid) {
            const RecordCount end_offset = global_data.idx2bp_offsets[cid + 1];

            if (end_offset == 0 || end_offset > static_cast<RecordCount>(global_data.idx2bp.size())) {
                throw std::runtime_error(
                    "Invalid idx2bp boundary for component " + std::to_string(cid) + "."
                );
            }

            const auto end_index = static_cast<std::size_t>(end_offset);

            // The last element of each local prefix-sum table stores the total component length.
            global_data.component_lengths.push_back(global_data.idx2bp[end_index - 1]);
        }

        // 7. Validate layout consistency.
        if (layout.component_names.size() != layout.component_count) {
            throw std::runtime_error("Component name count does not match the number of CDX components.");
        }
        if (global_data.component_lengths.size() != layout.component_count) {
            throw std::logic_error("Component length count does not match the number of CDX components.");
        }

        // 8. Move layout to data structure and return.
        global_data.layout = std::move(layout);
        return global_data;
    }

    void inspectComponent(
        const std::filesystem::path &cdx_path,
        const std::optional<Cid> component_id
    ) {
        std::ifstream cdx_stream(cdx_path, std::ios::binary);

        if (!cdx_stream) {
            throw std::runtime_error("Unable to open CDX file: " + cdx_path.string());
        }

        const FileHeader header = readGlobalHeader(cdx_stream);
        const auto component_count = header.n_components;

        // Validate requested component.
        if (component_id && *component_id >= component_count) {
            throw std::out_of_range("Component " + std::to_string(*component_id) + " does not exist. "
                                    "Valid range: [0-" + std::to_string(component_count - 1) + "].");
        }

        std::cout
                << "=========================================================================\n"
                << "CDX COMPONENT SUMMARY\n"
                << "File: " + cdx_path.filename().string() + "\n"
                << "=========================================================================\n\n";

        std::cout << std::left << std::setw(12) << "Component" << std::right << std::setw(18) << "Length(bp)"
                << std::setw(15) << "Nodes" << std::setw(25) << "NodeID Range" << '\n';

        std::cout << std::string(73, '-') << '\n';

        // Single-component inspection.
        if (component_id) {
            const auto component = seekComponent(cdx_stream, *component_id);
            const auto records = readComponentPayload(cdx_stream, component.nb_nodes);
            const auto idx2bp = buildIdx2Pos(records, component.nb_nodes);
            const PosBp component_length = idx2bp.back();

            std::ostringstream nid_range;

            nid_range << component.nid_min << "-" << component.nid_max;
            std::cout << std::left << std::setw(12) << *component_id
                    << std::right << std::setw(18) << component_length
                    << std::setw(15) << component.nb_nodes
                    << std::setw(25) << nid_range.str() << '\n';
            return;
        }

        // Full inspection.
        RecordCount total_nodes = 0;
        PosBp total_length_bp = 0;

        for (Cid compo_id = 0; compo_id < component_count; ++compo_id) {
            const auto component = readComponentHeader(cdx_stream, compo_id);
            const auto records = readComponentPayload(cdx_stream, component.nb_nodes);
            const auto idx2bp = buildIdx2Pos(records, component.nb_nodes);
            const PosBp component_length = idx2bp.back();

            total_nodes += component.nb_nodes;
            total_length_bp += component_length;

            std::ostringstream nid_range;

            nid_range << component.nid_min << "-" << component.nid_max;
            std::cout << std::left << std::setw(12) << compo_id
                    << std::right << std::setw(18) << component_length
                    << std::setw(15) << component.nb_nodes
                    << std::setw(25) << nid_range.str() << '\n';
        }

        std::cout << std::string(73, '-') << '\n';
        std::cout << std::left << std::setw(12) << "Total"
                << std::right << std::setw(18) << total_length_bp
                << std::setw(15) << total_nodes
                << std::setw(25) << "-" << '\n';
    }
}
