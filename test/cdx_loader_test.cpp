/**
 * @file cdx_loader_test.cpp
 * @brief Unit tests for cdx_loader.cpp (loadQuery, loadGlobal, inspectComponent).
 *
 * All logic under test lives in an anonymous namespace inside cdx_loader.cpp, so it is only
 * reachable through the three public entry points. Each test therefore builds a real, valid
 * binary .cdx fixture on disk (via the CdxFileBuilder helper below, which uses the same
 * cdx::CdxFormat::pack_* routines the production reader/writer rely on) and exercises the
 * loader through its public API.
 *
 * Test fixtures largely mirror the pre-existing Python reference test suite, translated to
 * the public C++ surface (the internal helpers such as _resolveQueryRange, _buildIdx2Pos,
 * _buildCovTableQuery, _computeGraphLayout, etc. are all indirectly exercised through
 * loadQuery/loadGlobal/inspectComponent). A dedicated section is also added for "funky" query
 * coordinate combinations (negative/huge/crossing endpoints) per explicit request, since these
 * are historically the most bug-prone inputs to loadQuery.
 */

#include "../src/cdx_loader.h"
#include "cdx_format.h"
#include "cdx_types.h"
#include "../src/config.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {
    // =====================================================================
    // Test infrastructure: builds real binary .cdx fixture files on disk.
    // =====================================================================

    /**
     * @brief Describes one component to serialize into a test .cdx fixture.
     *
     * nid_min/nid_max are supplied explicitly (rather than inferred from
     * `records`) so that tests can deliberately construct sparse node-ID
     * ranges (gaps) or edge-case component headers.
     */
    struct ComponentSpec {
        std::string name;
        cdx::Nid nid_min;
        cdx::Nid nid_max;
        // (node_id, local_idx, seq_len) triples, may be supplied out of order.
        std::vector<std::tuple<cdx::Nid, cdx::Idx, cdx::SeqLen> > records;
    };

    /**
     * @brief Writes a valid binary CDX file built from a list of ComponentSpec,
     *        using the production cdx::CdxFormat pack_* routines, and removes
     *        it automatically when the object goes out of scope.
     */
    class CdxFileFixture {
    public:
        explicit CdxFileFixture(const std::vector<ComponentSpec> &components) {
            static int counter = 0;
            path_ = std::filesystem::temp_directory_path() /
                    ("cdx_loader_test_" + std::to_string(++counter) + ".cdx");

            std::ofstream out(path_, std::ios::binary);

            std::array<char, cdx::CdxFormat::FILE_HEADER_SIZE> file_header_buf{};
            cdx::CdxFormat::pack_file_header(file_header_buf.data(), static_cast<uint32_t>(components.size()));
            out.write(file_header_buf.data(), static_cast<std::streamsize>(file_header_buf.size()));

            for (const auto &component: components) {
                std::array<char, cdx::CdxFormat::COMPONENT_HEADER_SIZE> comp_header_buf{};
                cdx::CdxFormat::pack_component_header(
                    comp_header_buf.data(),
                    component.records.size(),
                    component.nid_min,
                    component.nid_max,
                    static_cast<uint32_t>(component.name.size())
                );
                out.write(comp_header_buf.data(), static_cast<std::streamsize>(comp_header_buf.size()));
                out.write(component.name.data(), static_cast<std::streamsize>(component.name.size()));

                for (const auto &[node_id, local_idx, seq_len]: component.records) {
                    std::array<char, cdx::CdxFormat::RECORD_SIZE> record_buf{};
                    cdx::CdxFormat::pack_node_record(record_buf.data(), node_id, local_idx, seq_len);
                    out.write(record_buf.data(), static_cast<std::streamsize>(record_buf.size()));
                }
            }
        }

        ~CdxFileFixture() {
            std::error_code ec;
            std::filesystem::remove(path_, ec);
        }

        CdxFileFixture(const CdxFileFixture &) = delete;

        CdxFileFixture &operator=(const CdxFileFixture &) = delete;

        [[nodiscard]] const std::filesystem::path &path() const { return path_; }

    private:
        std::filesystem::path path_;
    };

    /**
     * @brief Convenience helper: builds a single-component fixture with node IDs
     *        inferred from `records` (min/max taken from the record list itself).
     */
    std::vector<ComponentSpec> singleComponent(
        const std::string &name,
        const std::vector<std::tuple<cdx::Nid, cdx::Idx, cdx::SeqLen> > &records
    ) {
        cdx::Nid nid_min = std::get<0>(records.front());
        cdx::Nid nid_max = nid_min;
        for (const auto &[nid, idx, len]: records) {
            nid_min = std::min(nid_min, nid);
            nid_max = std::max(nid_max, nid);
        }
        return {ComponentSpec{name, nid_min, nid_max, records}};
    }

    /**
     * @brief RAII helper redirecting std::cout into a string buffer for the
     *        lifetime of the object, restoring the original buffer on destruction.
     */
    class CoutCapture {
    public:
        CoutCapture() : old_buf_(std::cout.rdbuf(buffer_.rdbuf())) {
        }

        ~CoutCapture() { std::cout.rdbuf(old_buf_); }

        [[nodiscard]] std::string str() const { return buffer_.str(); }

    private:
        std::ostringstream buffer_;
        std::streambuf *old_buf_;
    };

    // Fixture component reused across the "funky query" section: 3 nodes of
    // 100bp each -> idx2bp = [0, 100, 200, 300], component_length = 300.
    std::vector<ComponentSpec> threeNodeFixture() {
        return singleComponent("chrX", {
                                    {10, 0, 100},
                                    {11, 1, 100},
                                    {12, 2, 100}
                                });
    }
} // anonymous namespace

// =============================================================================
// loadQuery
//
// Covers: default full-component span, partial ranges, negative-coordinate
// wrapping, node-boundary touching, circular vs. linear origin crossing,
// sparse node IDs, unordered on-disk records, and error propagation
// (missing file, out-of-range component, corrupted local index).
// =============================================================================
namespace {
    // Default query_range (nullopt) must resolve to the full component span.
    TEST(LoadQueryTest, DefaultRangeReturnsFullComponentSpan) {
        CdxFileFixture fixture(singleComponent("c0", {{10, 0, 100}, {11, 1, 200}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::nullopt, false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 299}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 1}));
        EXPECT_EQ(data.component.nid_min, 10u);
        EXPECT_EQ(data.idx2bp, (std::vector<cdx::PosBp>{0, 100, 300}));
        EXPECT_EQ(data.nid2idx, (std::vector<cdx::Idx>{0, 1}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, 0}));
    }

    // A sub-range is resolved, projected to idx space, and activates only the
    // covered node(s) in node_coverage.
    TEST(LoadQueryTest, PartialRangeActivatesOnlyCoveredNode) {
        CdxFileFixture fixture(singleComponent("c0", {{100, 0, 10}, {101, 1, 20}, {102, 2, 30}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(15, 25), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{15, 25}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{1, 1}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY, 0, cfg::NOT_IN_QUERY}));
    }

    // Circular components may specify a query range that crosses the origin.
    TEST(LoadQueryTest, CircularCrossOriginQueryIsAccepted) {
        CdxFileFixture fixture(singleComponent("c0", {{5, 0, 50}, {6, 1, 50}, {7, 2, 50}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(140, 10), true);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{140, 10}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{2, 0}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, cfg::NOT_IN_QUERY, 0}));
    }

    // The same origin-crossing query must be rejected with a helpful message
    // when the component is not marked circular.
    TEST(LoadQueryTest, LinearCrossOriginQueryThrowsWithHelpfulMessage) {
        CdxFileFixture fixture(singleComponent("c0", {{5, 0, 50}, {6, 1, 50}, {7, 2, 50}}));

        try {
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(140, 10), false);
            FAIL() << "Expected std::invalid_argument";
        } catch (const std::invalid_argument &e) {
            const std::string message = e.what();
            EXPECT_NE(message.find("crosses the component origin"), std::string::npos);
            EXPECT_NE(message.find("--component-type circular"), std::string::npos);
        }
    }

    // A single-node component must be handled consistently across every returned structure.
    TEST(LoadQueryTest, SingleNodeComponent) {
        CdxFileFixture fixture(singleComponent("c0", {{42, 0, 500}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::nullopt, false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 499}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 0}));
        EXPECT_EQ(data.idx2bp, (std::vector<cdx::PosBp>{0, 500}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0}));
    }

    // Negative coordinates must resolve relative to the component's end.
    TEST(LoadQueryTest, NegativeCoordinatesResolveRelativeToEnd) {
        CdxFileFixture fixture(singleComponent("c0", {{10, 0, 100}, {11, 1, 100}, {12, 2, 100}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-50, -1), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{250, 299}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{2, 2}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY, 0}));
    }

    // A query touching an exact node boundary must project to that node only.
    TEST(LoadQueryTest, QueryTouchingNodeBoundary) {
        CdxFileFixture fixture(singleComponent("c0", {{10, 0, 100}, {11, 1, 100}, {12, 2, 100}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(100, 199), false);

        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{1, 1}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY, 0, cfg::NOT_IN_QUERY}));
    }

    // Any node partially touched by the query must be fully included.
    TEST(LoadQueryTest, QueryOverlappingTwoNodesIncludesBoth) {
        CdxFileFixture fixture(singleComponent("c0", {{10, 0, 100}, {11, 1, 100}, {12, 2, 100}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(50, 150), false);

        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 1}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, 0, cfg::NOT_IN_QUERY}));
    }

    // The CDX format guarantees records are written in node_id-ascending order
    // (enforced by cdx_builder, and checked by cdx_lib in debug builds), but it
    // does NOT guarantee that the topological `idx` field matches that node_id
    // order (idx reflects graph layout order, node_id is just an identifier).
    // buildIdx2Pos/buildNid2IdxQuery must place every record by its `idx` field
    // regardless of the on-disk (node_id) order.
    TEST(LoadQueryTest, TopologicalIdxOrderMayDifferFromNodeIdOrder) {
        CdxFileFixture fixture(singleComponent("c0", {{10, 2, 30}, {11, 0, 10}, {12, 1, 20}}));

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::nullopt, false);

        EXPECT_EQ(data.idx2bp, (std::vector<cdx::PosBp>{0, 10, 30, 60}));
        EXPECT_EQ(data.nid2idx, (std::vector<cdx::Idx>{2, 0, 1})); // node 10->idx2, node 11->idx0, node 12->idx1
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, 0, 0}));
    }

    // Node IDs missing from a sparse component must read back as NOT_IN_COMPO
    // in both nid2idx and node_coverage.
    TEST(LoadQueryTest, SparseNodeIdsProduceNotInCompoSentinel) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 10, 20, {{10, 0, 10}, {15, 1, 10}, {20, 2, 10}}}
        };
        CdxFileFixture fixture(components);

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::nullopt, false);

        EXPECT_EQ(data.nid2idx[0], 0u); // node 10 (offset 0), present
        EXPECT_EQ(data.nid2idx[5], 1u); // node 15 (offset 5), present with local idx 1
        EXPECT_EQ(data.nid2idx[10], 2u); // node 20 (offset 10), present with local idx 2
        EXPECT_EQ(data.nid2idx[1], cfg::NOT_IN_COMPO); // node 11 (offset 1), absent
        EXPECT_EQ(data.node_coverage[1], cfg::NOT_IN_COMPO); // node 11 (absent)
    }

    // component_id beyond the archive's component count must propagate std::out_of_range.
    TEST(LoadQueryTest, OutOfRangeComponentIdThrowsOutOfRange) {
        CdxFileFixture fixture(singleComponent("c0", {{10, 0, 100}}));

        EXPECT_THROW(cdx::loadQuery(fixture.path(), 5, std::nullopt, false), std::out_of_range);
    }

    // A missing file must raise a clear runtime_error rather than crash.
    TEST(LoadQueryTest, MissingFileThrowsRuntimeError) {
        const std::filesystem::path missing = std::filesystem::temp_directory_path() / "does_not_exist_cdx_loader.cdx";

        EXPECT_THROW(cdx::loadQuery(missing, 0, std::nullopt, false), std::runtime_error);
    }

    // A corrupted record whose `idx` is out of bounds must now unconditionally
    // throw (previously guarded by #ifndef NDEBUG in buildIdx2Pos).
    TEST(LoadQueryTest, CorruptedRecordWithOutOfBoundsIdxThrows) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 100, 200, {{100, 0, 400}, {200, 5, 300}}} // idx 5 with only 2 records
        };
        CdxFileFixture fixture(components);

        EXPECT_THROW(cdx::loadQuery(fixture.path(), 0, std::nullopt, false), std::runtime_error);
    }

    // =========================================================================
    // "Funky" query coordinate combinations.
    //
    // Fixed fixture: 3 nodes of 100bp each on component "chrX"
    // -> idx2bp = [0, 100, 200, 300], component_length = 300.
    // =========================================================================

    // chrX 0:-100 -> resolved (0, 200): negative end wraps but the range does
    // not cross the origin, so it must succeed under both circular=true/false.
    TEST(LoadQueryTest, FunkyQuery_ZeroToNegative100_ResolvesNonCrossing) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(0, -100), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 200}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 2}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, 0, 0}));
    }

    // chrX -10:-100 -> resolved (290, 200): crosses the origin at the bp level.
    // Rejected in linear mode...
    TEST(LoadQueryTest, FunkyQuery_Neg10ToNeg100_LinearRejected) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-10, -100), false),
            std::invalid_argument
        );
    }

    // ...but accepted in circular mode. Interestingly, both resolved bp
    // endpoints (290 and 200) fall inside the *same* last node, so the
    // idx-level projection does NOT cross the origin even though the bp-level
    // range does: query_range_idx collapses to (2, 2), activating only that node.
    TEST(LoadQueryTest, FunkyQuery_Neg10ToNeg100_CircularCollapsesToSingleNode) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-10, -100), true);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{290, 200}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{2, 2}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY, 0}));
    }

    // chrX 10:-100 -> resolved (10, 200): positive start unchanged, negative
    // end wraps; does not cross the origin, must succeed regardless of circular.
    TEST(LoadQueryTest, FunkyQuery_10ToNegative100_ResolvesNonCrossing) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(10, -100), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{10, 200}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 2}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, 0, 0}));
    }

    // chrX 20:10 -> both positive, start > end: crosses the origin at the bp
    // level. Rejected in linear mode...
    TEST(LoadQueryTest, FunkyQuery_20To10_LinearRejected) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(20, 10), false),
            std::invalid_argument
        );
    }

    // ...but accepted in circular mode. Both 20bp and 10bp land inside node 0
    // ([0, 100)), so again the idx-level projection does NOT cross the origin:
    // query_range_idx collapses to (0, 0).
    TEST(LoadQueryTest, FunkyQuery_20To10_CircularCollapsesToSingleNode) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(20, 10), true);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{20, 10}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 0}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY}));
    }

    // chrX 0:0 -> trivial single base-pair query.
    TEST(LoadQueryTest, FunkyQuery_ZeroToZero_SingleBasePair) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(0, 0), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 0}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 0}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY}));
    }

    // chrX 999999999999:9999999999 -> start is checked (and rejected) first,
    // regardless of the (also out-of-range) end value.
    TEST(LoadQueryTest, FunkyQuery_HugeStartOutOfRangeThrowsWithStartMessage) {
        CdxFileFixture fixture(threeNodeFixture());

        try {
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(999999999999, 9999999999), false);
            FAIL() << "Expected std::out_of_range";
        } catch (const std::out_of_range &e) {
            EXPECT_NE(std::string(e.what()).find("start"), std::string::npos);
        }
    }

    // A valid start with a huge, out-of-range end must trip the end-specific message.
    TEST(LoadQueryTest, FunkyQuery_HugeEndOutOfRangeThrowsWithEndMessage) {
        CdxFileFixture fixture(threeNodeFixture());

        try {
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(10, 999999999999), false);
            FAIL() << "Expected std::out_of_range";
        } catch (const std::out_of_range &e) {
            EXPECT_NE(std::string(e.what()).find("end"), std::string::npos);
        }
    }

    // Extreme int64 endpoints must be rejected cleanly (no crash / UB), even
    // though the negative-wrap arithmetic (original + component_length) is
    // evaluated on the minimum representable int64 value.
    TEST(LoadQueryTest, FunkyQuery_Int64MinStartThrowsCleanly) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(
                fixture.path(), 0,
                std::make_pair<int64_t, int64_t>(std::numeric_limits<int64_t>::min(), 10),
                false
            ),
            std::out_of_range
        );
    }

    TEST(LoadQueryTest, FunkyQuery_Int64MaxEndThrowsCleanly) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(
                fixture.path(), 0,
                std::make_pair<int64_t, int64_t>(10, std::numeric_limits<int64_t>::max()),
                false
            ),
            std::out_of_range
        );
    }

    TEST(LoadQueryTest, FunkyQuery_Int64MinAndMaxThrowsCleanly) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(
                fixture.path(), 0,
                std::make_pair<int64_t, int64_t>(
                    std::numeric_limits<int64_t>::min(),
                    std::numeric_limits<int64_t>::max()
                ),
                false
            ),
            std::out_of_range
        );
    }

    // chrX 0:-0 -> literal negative zero. int64_t has no signed-zero distinction
    // so -0 == 0, meaning this must resolve identically to (0, 0). Included
    // explicitly since it looks alarming on a command line even though it is
    // numerically inert.
    TEST(LoadQueryTest, FunkyQuery_ZeroToNegativeZero_IsIdenticalToZeroZero) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(0, -0), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 0}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 0}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY}));
    }

    // chrX 299:299 -> the very last valid bp queried as a single position.
    TEST(LoadQueryTest, FunkyQuery_LastValidBpAsSinglePosition) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(299, 299), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{299, 299}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{2, 2}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY, 0}));
    }

    // chrX -300:299 -> negative start wraps to exactly 0 (-component_length),
    // equivalent to the default full-component span.
    TEST(LoadQueryTest, FunkyQuery_NegativeStartWrapsToExactZero) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-300, 299), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 299}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 2}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, 0, 0}));
    }

    // chrX -301:0 -> one past the negative wrap boundary (-component_length - 1)
    // resolves to -1, which must be rejected.
    TEST(LoadQueryTest, FunkyQuery_NegativeStartOneBeyondWrapBoundaryThrows) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-301, 0), false),
            std::out_of_range
        );
    }

    // chrX 300:300 -> start exactly equal to component_length (one past the
    // last valid index) must be rejected, not silently clamped.
    TEST(LoadQueryTest, FunkyQuery_StartEqualToComponentLengthThrows) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(300, 300), false),
            std::out_of_range
        );
    }

    // chrX 299:0, circular -> adjacent-boundary origin crossing covering just
    // the last and first node (2-node span, not a "collapse" like 20:10 above:
    // here the crossing genuinely survives at the idx level).
    TEST(LoadQueryTest, FunkyQuery_AdjacentBoundaryCrossOrigin_CircularSpansEnds) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(299, 0), true);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{299, 0}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{2, 0}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, cfg::NOT_IN_QUERY, 0}));
    }

    TEST(LoadQueryTest, FunkyQuery_AdjacentBoundaryCrossOrigin_LinearRejected) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(299, 0), false),
            std::invalid_argument
        );
    }

    // chrX -50:-50 -> both endpoints negative and identical; must resolve to a
    // single valid position deep inside the component (250), not be mistaken
    // for an empty/crossing range.
    TEST(LoadQueryTest, FunkyQuery_BothEndpointsSameNegativeValue) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-50, -50), false);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{250, 250}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{2, 2}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY, 0}));
    }

    // chrX -10:5 -> negative start resolves to 290, crosses the origin against
    // a small positive end. Rejected in linear mode...
    TEST(LoadQueryTest, FunkyQuery_Neg10To5_LinearRejected) {
        CdxFileFixture fixture(threeNodeFixture());

        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-10, 5), false),
            std::invalid_argument
        );
    }

    // ...and accepted (surviving as a genuine idx-level crossing, unlike the
    // -10:-100 case above) in circular mode.
    TEST(LoadQueryTest, FunkyQuery_Neg10To5_CircularAccepted) {
        CdxFileFixture fixture(threeNodeFixture());

        const cdx::QueryData data = cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-10, 5), true);

        EXPECT_EQ(data.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{290, 5}));
        EXPECT_EQ(data.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{2, 0}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, cfg::NOT_IN_QUERY, 0}));
    }

    // A single-bp component (length 1) exercises the tightest possible bounds:
    // both None and (0, 0) must resolve identically, and anything else is out of range.
    TEST(LoadQueryTest, SingleBpComponent_DefaultAndExplicitZeroZeroAgree) {
        CdxFileFixture fixture(singleComponent("c0", {{5, 0, 1}}));

        const cdx::QueryData default_range = cdx::loadQuery(fixture.path(), 0, std::nullopt, false);
        EXPECT_EQ(default_range.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 0}));

        const cdx::QueryData explicit_range = cdx::loadQuery(
            fixture.path(), 0, std::make_pair<int64_t, int64_t>(0, 0), false
        );
        EXPECT_EQ(explicit_range.query_range_bp, (std::pair<cdx::PosBp, cdx::PosBp>{0, 0}));
        EXPECT_EQ(explicit_range.query_range_idx, (std::pair<cdx::Idx, cdx::Idx>{0, 0}));
    }

    TEST(LoadQueryTest, SingleBpComponent_AnyOtherRangeThrows) {
        CdxFileFixture fixture(singleComponent("c0", {{5, 0, 1}}));

        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(1, 1), false),
            std::out_of_range
        );
        EXPECT_THROW(
            cdx::loadQuery(fixture.path(), 0, std::make_pair<int64_t, int64_t>(-2, 0), false),
            std::out_of_range
        );
    }
} // anonymous namespace

// =============================================================================
// loadGlobal
//
// Covers: single/multiple-component aggregation, sparse global node-ID space,
// component-local coordinate independence, cross-structure consistency, and
// zero-record component rejection.
// =============================================================================
namespace {
    TEST(LoadGlobalTest, SingleComponentGraph) {
        CdxFileFixture fixture(singleComponent("c0", {{50, 0, 40}, {51, 1, 60}}));

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        EXPECT_EQ(data.layout.graph_nid_min, 50u);
        EXPECT_EQ(data.layout.component_offsets, (std::vector<cdx::RecordCount>{0, 2}));
        EXPECT_EQ(data.idx2bp_offsets, (std::vector<cdx::RecordCount>{0, 3}));
        EXPECT_EQ(data.idx2bp, (std::vector<cdx::PosBp>{0, 40, 100}));
        EXPECT_EQ(data.nid2flat_idx, (std::vector<cdx::FlatIdx>{0, 1}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0, 0}));
        EXPECT_EQ(data.component_lengths, (std::vector<cdx::PosBp>{100}));
    }

    TEST(LoadGlobalTest, MultipleComponentsAndSparseNids) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 10, 11, {{10, 0, 100}, {11, 1, 100}}},
            ComponentSpec{"c1", 20, 22, {{20, 0, 50}, {22, 1, 50}}}
        };
        CdxFileFixture fixture(components);

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        EXPECT_EQ(data.layout.graph_nid_min, 10u);
        EXPECT_EQ(data.layout.component_offsets, (std::vector<cdx::RecordCount>{0, 2, 4}));
        EXPECT_EQ(data.idx2bp_offsets, (std::vector<cdx::RecordCount>{0, 3, 6}));
        EXPECT_EQ(data.idx2bp, (std::vector<cdx::PosBp>{0, 100, 200, 0, 50, 100}));

        EXPECT_EQ(data.nid2flat_idx[10 - 10], 0u);
        EXPECT_EQ(data.nid2flat_idx[11 - 10], 1u);
        EXPECT_EQ(data.nid2flat_idx[15 - 10], cfg::INVALID_FLAT_IDX);
        EXPECT_EQ(data.nid2flat_idx[20 - 10], 2u);
        EXPECT_EQ(data.nid2flat_idx[21 - 10], cfg::INVALID_FLAT_IDX);
        EXPECT_EQ(data.nid2flat_idx[22 - 10], 3u);
    }

    // Cross-validates that consecutive component boundaries stay consistent
    // between component_offsets and idx2bp_offsets (idx2bp always has one more
    // element per component than there are node records).
    TEST(LoadGlobalTest, InternalConsistencyBetweenOffsetArrays) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 100, 101, {{100, 0, 30}, {101, 1, 40}}},
            ComponentSpec{"c1", 200, 202, {{200, 0, 10}, {202, 1, 20}}}
        };
        CdxFileFixture fixture(components);

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        ASSERT_EQ(data.layout.component_offsets.size(), data.idx2bp_offsets.size());

        for (std::size_t i = 0; i + 1 < data.layout.component_offsets.size(); ++i) {
            const auto n_records_compo = data.layout.component_offsets[i + 1] - data.layout.component_offsets[i];
            const auto idx2pos_len = data.idx2bp_offsets[i + 1] - data.idx2bp_offsets[i];
            EXPECT_EQ(idx2pos_len, n_records_compo + 1);
        }
    }

    TEST(LoadGlobalTest, SingleNodeGraph) {
        CdxFileFixture fixture(singleComponent("c0", {{42, 0, 500}}));

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        EXPECT_EQ(data.layout.graph_nid_min, 42u);
        EXPECT_EQ(data.nid2flat_idx, (std::vector<cdx::FlatIdx>{0}));
        EXPECT_EQ(data.node_coverage, (std::vector<cdx::Coverage>{0}));
    }

    // See LoadQueryTest.TopologicalIdxOrderMayDifferFromNodeIdOrder: node_id
    // order is guaranteed on disk, but the topological `idx` order is not.
    TEST(LoadGlobalTest, TopologicalIdxOrderMayDifferFromNodeIdOrder) {
        CdxFileFixture fixture(singleComponent("c0", {{10, 2, 30}, {11, 0, 10}, {12, 1, 20}}));

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        EXPECT_EQ(data.idx2bp, (std::vector<cdx::PosBp>{0, 10, 30, 60}));
        EXPECT_EQ(data.nid2flat_idx, (std::vector<cdx::FlatIdx>{2, 0, 1}));
    }

    TEST(LoadGlobalTest, LargeNidGapBetweenComponentsPreservesSentinels) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 1, 2, {{1, 0, 10}, {2, 1, 10}}},
            ComponentSpec{"c1", 1000, 1001, {{1000, 0, 10}, {1001, 1, 10}}}
        };
        CdxFileFixture fixture(components);

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        EXPECT_EQ(data.nid2flat_idx[0], 0u);
        EXPECT_EQ(data.nid2flat_idx[1], 1u);
        EXPECT_EQ(data.nid2flat_idx[500], cfg::INVALID_FLAT_IDX);
        EXPECT_EQ(data.nid2flat_idx[999], 2u);
        EXPECT_EQ(data.nid2flat_idx[1000], 3u);
    }

    // Each component keeps an independent, zero-based local coordinate system
    // inside the concatenated idx2bp array.
    TEST(LoadGlobalTest, ComponentLocalCoordinatesRestartAtZero) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 10, 11, {{10, 0, 100}, {11, 1, 200}}},
            ComponentSpec{"c1", 20, 21, {{20, 0, 10}, {21, 1, 20}}}
        };
        CdxFileFixture fixture(components);

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        const std::vector compo0(
            data.idx2bp.begin() + static_cast<long>(data.idx2bp_offsets[0]),
            data.idx2bp.begin() + static_cast<long>(data.idx2bp_offsets[1])
        );
        const std::vector compo1(
            data.idx2bp.begin() + static_cast<long>(data.idx2bp_offsets[1]),
            data.idx2bp.begin() + static_cast<long>(data.idx2bp_offsets[2])
        );

        EXPECT_EQ(compo0, (std::vector<cdx::PosBp>{0, 100, 300}));
        EXPECT_EQ(compo1, (std::vector<cdx::PosBp>{0, 10, 30}));
    }

    // Every valid flat_idx must correspond to an activated (0) coverage entry,
    // and every sentinel flat_idx to a NOT_IN_COMPO coverage entry.
    TEST(LoadGlobalTest, CovTableMatchesFlatIdxPresence) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 10, 11, {{10, 0, 10}, {11, 1, 10}}},
            ComponentSpec{"c1", 20, 20, {{20, 0, 10}}}
        };
        CdxFileFixture fixture(components);

        const cdx::GlobalData data = cdx::loadGlobal(fixture.path());

        for (std::size_t offset = 0; offset < data.nid2flat_idx.size(); ++offset) {
            if (data.nid2flat_idx[offset] == cfg::INVALID_FLAT_IDX) {
                EXPECT_EQ(data.node_coverage[offset], cfg::NOT_IN_COMPO);
            } else {
                EXPECT_EQ(data.node_coverage[offset], 0u);
            }
        }
    }

    // A component with zero node records must be rejected (guaranteed elsewhere
    // that real archives never contain one, but the reader must fail cleanly
    // rather than divide-by-zero / produce a bogus empty component length).
    TEST(LoadGlobalTest, ZeroRecordComponentThrows) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 10, 10, {{10, 0, 50}}},
            ComponentSpec{"c1", 20, 20, {}}
        };
        CdxFileFixture fixture(components);

        EXPECT_THROW(cdx::loadGlobal(fixture.path()), std::runtime_error);
    }

    TEST(LoadGlobalTest, MissingFileThrowsRuntimeError) {
        const std::filesystem::path missing = std::filesystem::temp_directory_path() / "does_not_exist_global.cdx";

        EXPECT_THROW(cdx::loadGlobal(missing), std::runtime_error);
    }
} // anonymous namespace

// =============================================================================
// inspectComponent
//
// Prints a formatted summary to stdout; assertions check substrings/columns
// rather than pinning the exact layout, to stay robust to minor formatting
// tweaks.
// =============================================================================
namespace {
    std::vector<ComponentSpec> twoComponentFixture() {
        return {
            ComponentSpec{"c0", 10, 11, {{10, 0, 100}, {11, 1, 200}}},
            ComponentSpec{"c1", 20, 21, {{20, 0, 10}, {21, 1, 20}}}
        };
    }

    // Inspecting a single component prints only that component's row, no total.
    TEST(InspectComponentTest, SingleComponentPrintsOwnMetricsOnly) {
        CdxFileFixture fixture(twoComponentFixture());
        CoutCapture capture;

        cdx::inspectComponent(fixture.path(), 0);
        const std::string output = capture.str();

        EXPECT_NE(output.find("CDX COMPONENT SUMMARY"), std::string::npos);
        EXPECT_NE(output.find("300"), std::string::npos); // component length (100 + 200)
        EXPECT_NE(output.find("10-11"), std::string::npos);
        EXPECT_EQ(output.find("20-21"), std::string::npos);
        EXPECT_EQ(output.find("Total"), std::string::npos);
    }

    // Inspecting all components prints every row plus a trailing "Total" row
    // whose bp/node counts are the exact sum across components.
    TEST(InspectComponentTest, AllComponentsPrintEveryRowAndExactTotal) {
        CdxFileFixture fixture(twoComponentFixture());
        CoutCapture capture;

        cdx::inspectComponent(fixture.path(), std::nullopt);
        const std::string output = capture.str();

        EXPECT_NE(output.find("10-11"), std::string::npos);
        EXPECT_NE(output.find("20-21"), std::string::npos);
        EXPECT_NE(output.find("Total"), std::string::npos);
        EXPECT_NE(output.find("330"), std::string::npos); // 300 + 30 total bp
    }

    // An out-of-range component ID must propagate std::out_of_range with the
    // valid range spelled out in the message.
    TEST(InspectComponentTest, InvalidComponentIdThrowsOutOfRangeWithValidRange) {
        CdxFileFixture fixture(twoComponentFixture());

        try {
            cdx::inspectComponent(fixture.path(), 99);
            FAIL() << "Expected std::out_of_range";
        } catch (const std::out_of_range &e) {
            const std::string message = e.what();
            EXPECT_NE(message.find("99"), std::string::npos);
            EXPECT_NE(message.find("[0-1]"), std::string::npos);
        }
    }

    // Component names longer than the display column must be truncated with a
    // trailing '~' marker rather than corrupting the table layout.
    TEST(InspectComponentTest, OverlongComponentNameIsTruncatedWithTilde) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"this_name_is_way_too_long_for_the_column", 10, 10, {{10, 0, 50}}}
        };
        CdxFileFixture fixture(components);
        CoutCapture capture;

        cdx::inspectComponent(fixture.path(), 0);
        const std::string output = capture.str();

        EXPECT_NE(output.find('~'), std::string::npos);
    }

    // Large lengths/node IDs must use thousands separators (cfg::formatInteger).
    TEST(InspectComponentTest, LargeValuesUseThousandsSeparators) {
        const std::vector<ComponentSpec> components = {
            ComponentSpec{"c0", 10000, 10001, {{10000, 0, 1000}, {10001, 1, 1000}}}
        };
        CdxFileFixture fixture(components);
        CoutCapture capture;

        cdx::inspectComponent(fixture.path(), 0);
        const std::string output = capture.str();

        EXPECT_NE(output.find("2,000"), std::string::npos);
        EXPECT_NE(output.find("10,000-10,001"), std::string::npos);
    }

    // An empty CDX archive (0 components) cannot be constructed by
    // CdxFileFixture through the normal builder path used above, but the
    // reader (cdx_lib's readGlobalHeader) already rejects it upstream; verify
    // that inspectComponent surfaces that failure as a runtime_error.
    TEST(InspectComponentTest, EmptyArchiveThrowsRuntimeError) {
        CdxFileFixture fixture(std::vector<ComponentSpec>{});

        EXPECT_THROW(cdx::inspectComponent(fixture.path(), std::nullopt), std::runtime_error);
    }
} // anonymous namespace
