/**
 * @file cov_projection_test.cpp
 * @brief Unit tests for src/cov_projection.cpp.
 *
 * One GoogleTest suite per function declared in cov_projection.h. Each
 * suite is preceded by a short block comment describing what that function
 * is responsible for and the angle its test cases cover; each individual
 * TEST() is preceded by a one-line comment describing that specific case.
 *
 * Note: cfg::NOT_IN_COMPO and cfg::INVALID_NODE are numerically identical
 * (both 0xFFFFFFFF, see config.h) even though they document distinct
 * concepts (an input coverage sentinel vs. an output/default sentinel).
 * Several cases below rely on that identity being preserved verbatim by
 * straight pass-through, not on any special-cased conversion.
 */

#include "../src/cov_projection.h"
#include "../src/config.h"
#include "cdx_types.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <utility>
#include <vector>

// =============================================================================
// projectCov2IdxQuery
//
// Remaps node-level coverage (indexed by relative nid offset) into a
// component-local index space via nid2idx. Two independent sentinel
// concerns are exercised here:
//   - sentinel *index* values in nid2idx (cfg::NOT_IN_QUERY / NOT_IN_COMPO)
//     mark a node as excluded from the local index space and must be
//     skipped without writing anything;
//   - sentinel *coverage* values in cov_table are ordinary payload as far
//     as this function is concerned and must simply pass through unchanged
//     at whatever local_idx they get mapped to.
// Bounds violations on local_idx are validated unconditionally (not just
// in debug builds), since they reflect malformed input, not an internal
// invariant.
// =============================================================================

// Identity mapping: nid offset order already matches local_idx order.
TEST(ProjectCov2IdxQueryTest, LinearIdentityMapping) {
    const std::vector<cdx::Coverage> cov = {10, 20, 30};
    const std::vector<cdx::Idx> nid2idx = {0, 1, 2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{10, 20, 30}));
}

// Coverage must follow nid2idx, not physical array order.
TEST(ProjectCov2IdxQueryTest, UnorderedIdxMapping) {
    const std::vector<cdx::Coverage> cov = {100, 200, 300};
    const std::vector<cdx::Idx> nid2idx = {2, 0, 1};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{200, 300, 100}));
}

// cov_table and nid2idx must share identical sizes.
TEST(ProjectCov2IdxQueryTest, DimensionMismatchThrows) {
    const std::vector<cdx::Coverage> cov = {1, 2};
    const std::vector<cdx::Idx> nid2idx = {0, 1, 2};

    EXPECT_THROW(projectCov2IdxQuery(cov, nid2idx, 3), std::invalid_argument);
}

// A NOT_IN_COMPO coverage *value* is opaque payload here: it must survive
// unchanged, since only nid2idx entries are inspected for sentinels.
TEST(ProjectCov2IdxQueryTest, SentinelCoverageValueNotInCompoIsPreserved) {
    const std::vector<cdx::Coverage> cov = {5, cfg::NOT_IN_COMPO, 7};
    const std::vector<cdx::Idx> nid2idx = {0, 1, 2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{5, cfg::NOT_IN_COMPO, 7}));
}

// Same as above for NOT_IN_QUERY appearing as a coverage value.
TEST(ProjectCov2IdxQueryTest, SentinelCoverageValueNotInQueryIsPreserved) {
    const std::vector<cdx::Coverage> cov = {5, cfg::NOT_IN_QUERY, 7};
    const std::vector<cdx::Idx> nid2idx = {0, 1, 2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{5, cfg::NOT_IN_QUERY, 7}));
}

// Both coverage-value sentinels together, still simple pass-through.
TEST(ProjectCov2IdxQueryTest, BothSentinelCoverageValuesArePreserved) {
    const std::vector<cdx::Coverage> cov = {cfg::NOT_IN_COMPO, cfg::NOT_IN_QUERY, 12};
    const std::vector<cdx::Idx> nid2idx = {0, 1, 2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{cfg::NOT_IN_COMPO, cfg::NOT_IN_QUERY, 12}));
}

// A NOT_IN_QUERY *index* in nid2idx must be skipped: nothing is written
// for that node, leaving the destination slot at its zero default.
TEST(ProjectCov2IdxQueryTest, NidIdxSentinelNotInQueryIsSkipped) {
    const std::vector<cdx::Coverage> cov = {10, 20, 30};
    const std::vector<cdx::Idx> nid2idx = {0, cfg::NOT_IN_QUERY, 2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{10, 0, 30}));
}

// A NOT_IN_COMPO *index* in nid2idx must likewise be skipped.
TEST(ProjectCov2IdxQueryTest, NidIdxSentinelNotInCompoIsSkipped) {
    const std::vector<cdx::Coverage> cov = {10, 20, 30};
    const std::vector<cdx::Idx> nid2idx = {0, cfg::NOT_IN_COMPO, 2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{10, 0, 30}));
}

// Large nid-space gaps must not affect the projection.
TEST(ProjectCov2IdxQueryTest, LargeSparseIdxSpace) {
    const std::vector<cdx::Coverage> cov = {50, 60, 70};
    const std::vector<cdx::Idx> nid2idx = {99, 0, 50};

    const auto result = projectCov2IdxQuery(cov, nid2idx, 100);

    EXPECT_EQ(result[0], 60);
    EXPECT_EQ(result[50], 70);
    EXPECT_EQ(result[99], 50);
}

// 0 is a legitimate coverage value, not a sentinel; must round-trip.
TEST(ProjectCov2IdxQueryTest, CoverageZeroIsValid) {
    const std::vector<cdx::Coverage> cov = {0, 0, 5};
    const std::vector<cdx::Idx> nid2idx = {0, 1, 2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{0, 0, 5}));
}

// The largest coverage value just below the sentinel range must survive.
TEST(ProjectCov2IdxQueryTest, MaxValidCoverageBelowSentinel) {
    const std::vector<cdx::Coverage> cov = {cfg::NOT_IN_QUERY - 1};
    const std::vector<cdx::Idx> nid2idx = {0};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 1),
              (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY - 1}));
}

// Degenerate one-node component.
TEST(ProjectCov2IdxQueryTest, SingleNodeComponent) {
    const std::vector<cdx::Coverage> cov = {42};
    const std::vector<cdx::Idx> nid2idx = {0};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 1), (std::vector<cdx::Coverage>{42}));
}

// Output size must always equal the requested component_size, independent
// of how many entries were actually written.
TEST(ProjectCov2IdxQueryTest, OutputLengthMatchesComponentSize) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3};
    const std::vector<cdx::Idx> nid2idx = {1, 2, 0};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3).size(), 3u);
}

// Repeated calls with identical inputs must produce identical outputs.
TEST(ProjectCov2IdxQueryTest, RepeatedExecutionIsDeterministic) {
    const std::vector<cdx::Coverage> cov = {7, 8, 9};
    const std::vector<cdx::Idx> nid2idx = {2, 0, 1};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3), projectCov2IdxQuery(cov, nid2idx, 3));
}

// local_idx exceeding component_size must raise, unconditionally (not
// just in debug builds) since it reflects malformed input data.
TEST(ProjectCov2IdxQueryTest, LocalIdxOutOfBoundsThrows) {
    const std::vector<cdx::Coverage> cov = {10};
    const std::vector<cdx::Idx> nid2idx = {5};

    EXPECT_THROW(projectCov2IdxQuery(cov, nid2idx, 3), std::out_of_range);
}

// Boundary case: local_idx == component_size - 1 is the largest valid
// index and must succeed (distinguishes strict bounds handling from an
// accidental off-by-one).
TEST(ProjectCov2IdxQueryTest, ExactUpperBoundaryIndexSucceeds) {
    const std::vector<cdx::Coverage> cov = {42};
    const std::vector<cdx::Idx> nid2idx = {2};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{0, 0, 42}));
}

// Two nodes mapping to the same local_idx: the later one wins.
TEST(ProjectCov2IdxQueryTest, DuplicateLocalIdxOverwrites) {
    const std::vector<cdx::Coverage> cov = {10, 20};
    const std::vector<cdx::Idx> nid2idx = {1, 1};

    EXPECT_EQ(projectCov2IdxQuery(cov, nid2idx, 3),
              (std::vector<cdx::Coverage>{0, 20, 0}));
}

// Empty inputs with a zero-sized component return an empty vector.
TEST(ProjectCov2IdxQueryTest, EmptyInputTablesReturnEmpty) {
    const std::vector<cdx::Coverage> cov;
    const std::vector<cdx::Idx> nid2idx;

    EXPECT_TRUE(projectCov2IdxQuery(cov, nid2idx, 0).empty());
}

// =============================================================================
// expandPosCovQuery
//
// Expands node-level coverage into per-base-pair depth for a single
// component, using idx2bp as a cumulative base-pair prefix-sum table
// (size == node_count + 1). A component is guaranteed by construction of
// the CDX index to contain at least one node, so a degenerate idx2bp of
// size < 2 is treated as malformed input, not a valid empty result.
// =============================================================================

// A single node expands to its full bp length.
TEST(ExpandPosCovQueryTest, SingleNode) {
    const std::vector<cdx::Coverage> cov = {20};
    const std::vector<cdx::PosBp> idx2bp = {0, 5};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp),
              (std::vector<cdx::Coverage>{20, 20, 20, 20, 20}));
}

// A one-bp node expands to exactly one bp of coverage.
TEST(ExpandPosCovQueryTest, SingleBpNode) {
    const std::vector<cdx::Coverage> cov = {7};
    const std::vector<cdx::PosBp> idx2bp = {0, 1};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp), (std::vector<cdx::Coverage>{7}));
}

// Nodes with heterogeneous lengths expand correctly and contiguously.
TEST(ExpandPosCovQueryTest, MultipleNodesVariableLengths) {
    const std::vector<cdx::Coverage> cov = {20, 15, 33};
    const std::vector<cdx::PosBp> idx2bp = {0, 5, 8, 13};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp),
              (std::vector<cdx::Coverage>{
                  20, 20, 20, 20, 20,
                  15, 15, 15,
                  33, 33, 33, 33, 33,
              }));
}

// All-zero node coverage remains zero after expansion.
TEST(ExpandPosCovQueryTest, AllZeroCoverage) {
    const std::vector<cdx::Coverage> cov = {0, 0, 0};
    const std::vector<cdx::PosBp> idx2bp = {0, 5, 10, 15};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp), (std::vector<cdx::Coverage>(15, 0)));
}

// Zero-coverage nodes interleaved with covered nodes expand independently.
TEST(ExpandPosCovQueryTest, MixedZeroAndNonzeroCoverage) {
    const std::vector<cdx::Coverage> cov = {10, 0, 30};
    const std::vector<cdx::PosBp> idx2bp = {0, 3, 6, 9};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp),
              (std::vector<cdx::Coverage>{10, 10, 10, 0, 0, 0, 30, 30, 30}));
}

// A large uncovered region in the middle stays fully zero-filled.
TEST(ExpandPosCovQueryTest, LargeHoleOfZeroCoverage) {
    const std::vector<cdx::Coverage> cov = {50, 0, 0, 0, 75};
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 4, 6, 8, 10};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp),
              (std::vector<cdx::Coverage>{50, 50, 0, 0, 0, 0, 0, 0, 75, 75}));
}

// Dense component with identical single-bp coverage everywhere.
TEST(ExpandPosCovQueryTest, DenseUniformCoverage) {
    const std::vector<cdx::Coverage> cov(100, 42);
    std::vector<cdx::PosBp> idx2bp(101);
    for (std::size_t i = 0; i < idx2bp.size(); ++i) {
        idx2bp[i] = static_cast<cdx::PosBp>(i);
    }

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp), (std::vector<cdx::Coverage>(100, 42)));
}

// Only one node in the middle contributes non-zero coverage.
TEST(ExpandPosCovQueryTest, SparseSingleActiveNode) {
    const std::vector<cdx::Coverage> cov = {0, 0, 99, 0, 0};
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 4, 6, 8, 10};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp),
              (std::vector<cdx::Coverage>{0, 0, 0, 0, 99, 99, 0, 0, 0, 0}));
}

// Coverage must change exactly at idx2bp boundaries, no off-by-one drift.
TEST(ExpandPosCovQueryTest, NodeBoundariesAreExact) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 5, 9};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp),
              (std::vector<cdx::Coverage>{1, 1, 2, 2, 2, 3, 3, 3, 3}));
}

// Largest valid coverage value, just below the sentinel range, survives.
TEST(ExpandPosCovQueryTest, MaxValidCoverageValue) {
    const std::vector<cdx::Coverage> cov = {cfg::NOT_IN_QUERY - 1};
    const std::vector<cdx::PosBp> idx2bp = {0, 4};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp),
              (std::vector<cdx::Coverage>(4, cfg::NOT_IN_QUERY - 1)));
}

// idx2bp must contain exactly one more entry than idx_cov_table.
TEST(ExpandPosCovQueryTest, DimensionMismatchThrows) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 5, 10};

    EXPECT_THROW(expandPosCovQuery(cov, idx2bp), std::invalid_argument);
}

// A degenerate idx2bp of size < 2 (0-node component) is rejected: the CDX
// format guarantees every component has at least one node, so this is
// treated as malformed input rather than a valid "empty" result.
TEST(ExpandPosCovQueryTest, DegenerateEmptyComponentThrows) {
    const std::vector<cdx::Coverage> cov;
    const std::vector<cdx::PosBp> idx2bp = {0};

    EXPECT_THROW(expandPosCovQuery(cov, idx2bp), std::invalid_argument);
}

// Output length must equal the component's total bp length.
TEST(ExpandPosCovQueryTest, OutputLengthEqualsComponentLength) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 3, 7, 15};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp).size(), 15u);
}

// Position space need not start at zero.
TEST(ExpandPosCovQueryTest, IdxOffsetStartNotZero) {
    const std::vector<cdx::Coverage> cov = {10, 20};
    const std::vector<cdx::PosBp> idx2bp = {10, 15, 20};

    const auto result = expandPosCovQuery(cov, idx2bp);

    ASSERT_EQ(result.size(), 20u);
    EXPECT_EQ(std::vector<cdx::Coverage>(result.begin() + 10, result.begin() + 15),
              (std::vector<cdx::Coverage>(5, 10)));
    EXPECT_EQ(std::vector<cdx::Coverage>(result.begin() + 15, result.begin() + 20),
              (std::vector<cdx::Coverage>(5, 20)));
}

// Each singleton (1-bp) node occupies exactly one output position.
TEST(ExpandPosCovQueryTest, ManySingleBpNodes) {
    std::vector<cdx::Coverage> cov(64);
    std::vector<cdx::PosBp> idx2bp(65);
    for (std::size_t i = 0; i < cov.size(); ++i) {
        cov[i] = static_cast<cdx::Coverage>(i + 1);
    }
    for (std::size_t i = 0; i < idx2bp.size(); ++i) {
        idx2bp[i] = static_cast<cdx::PosBp>(i);
    }

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp), cov);
}

// Repeated calls with identical inputs must produce identical outputs.
TEST(ExpandPosCovQueryTest, RepeatedExecutionIsDeterministic) {
    const std::vector<cdx::Coverage> cov = {5, 10, 15};
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 5, 8};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp), expandPosCovQuery(cov, idx2bp));
}

// Total expanded bp coverage must equal sum(node_coverage * node_length).
TEST(ExpandPosCovQueryTest, WeightedCoverageSumInvariant) {
    const std::vector<cdx::Coverage> cov = {10, 0, 25, 40};
    const std::vector<cdx::PosBp> idx2bp = {0, 100, 250, 300, 1000};

    const auto result = expandPosCovQuery(cov, idx2bp);

    std::uint64_t expected = 0;
    for (std::size_t i = 0; i < cov.size(); ++i) {
        expected += static_cast<std::uint64_t>(cov[i]) * (idx2bp[i + 1] - idx2bp[i]);
    }
    std::uint64_t actual = 0;
    for (const auto v : result) actual += v;

    EXPECT_EQ(actual, expected);
}

// The last base pair of the component must be written, not left off by
// one past the final fill.
TEST(ExpandPosCovQueryTest, LastBpIsWritten) {
    const std::vector<cdx::Coverage> cov = {5};
    const std::vector<cdx::PosBp> idx2bp = {0, 4};

    EXPECT_EQ(expandPosCovQuery(cov, idx2bp), (std::vector<cdx::Coverage>{5, 5, 5, 5}));
}

// A non-monotonic idx2bp (start > end for some node) reflects malformed
// input and must raise unconditionally.
TEST(ExpandPosCovQueryTest, NonMonotonicPrefixSumThrows) {
    const std::vector<cdx::Coverage> cov = {1, 2};
    const std::vector<cdx::PosBp> idx2bp = {0, 5, 3};

    EXPECT_THROW(expandPosCovQuery(cov, idx2bp), std::runtime_error);
}

// =============================================================================
// projectCov2IdxGlobal
//
// Same remapping idea as projectCov2IdxQuery but graph-wide: nid2flat_idx
// maps node offsets into a single flattened index space spanning every
// component, bounded by component_offsets. INVALID_FLAT_IDX marks a node
// as unmapped and must be skipped; other sentinel values only ever appear
// as coverage payload here and pass through unchanged.
// =============================================================================

// Identity flat mapping within a single component.
TEST(ProjectCov2IdxGlobalTest, IdentityFlatMapping) {
    const std::vector<cdx::Coverage> cov = {10, 20, 30};
    const std::vector<cdx::FlatIdx> nid2flat = {0, 1, 2};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{10, 20, 30}));
}

// Coverage must follow nid2flat_idx order, not physical array order.
TEST(ProjectCov2IdxGlobalTest, UnorderedFlatMapping) {
    const std::vector<cdx::Coverage> cov = {100, 200, 300};
    const std::vector<cdx::FlatIdx> nid2flat = {2, 0, 1};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{200, 300, 100}));
}

// Two components decode into disjoint, correctly bounded output ranges.
TEST(ProjectCov2IdxGlobalTest, MultipleComponentsDecodingBoundaries) {
    const std::vector<cdx::Coverage> cov = {10, 20, 30, 40};
    const std::vector<cdx::FlatIdx> nid2flat = {0, 1, 2, 3};
    const std::vector<cdx::RecordCount> comp = {0, 2, 4};

    const auto result = projectCov2IdxGlobal(cov, nid2flat, comp);

    EXPECT_EQ(std::vector<cdx::Coverage>(result.begin(), result.begin() + 2),
              (std::vector<cdx::Coverage>{10, 20}));
    EXPECT_EQ(std::vector<cdx::Coverage>(result.begin() + 2, result.begin() + 4),
              (std::vector<cdx::Coverage>{30, 40}));
}

// Writing into one component's flat range must not disturb another's.
TEST(ProjectCov2IdxGlobalTest, ComponentIsolationIntegrity) {
    const std::vector<cdx::Coverage> cov = {50, 99};
    const std::vector<cdx::FlatIdx> nid2flat = {2, 3};
    const std::vector<cdx::RecordCount> comp = {0, 2, 4};

    const auto result = projectCov2IdxGlobal(cov, nid2flat, comp);

    EXPECT_EQ(std::vector<cdx::Coverage>(result.begin(), result.begin() + 2),
              (std::vector<cdx::Coverage>{0, 0}));
    EXPECT_EQ(std::vector<cdx::Coverage>(result.begin() + 2, result.begin() + 4),
              (std::vector<cdx::Coverage>{50, 99}));
}

// A large, sparsely populated flat index space still projects correctly.
TEST(ProjectCov2IdxGlobalTest, LargeSparseFlatSpace) {
    const std::vector<cdx::Coverage> cov = {50, 60, 70};
    const std::vector<cdx::FlatIdx> nid2flat = {99, 0, 50};
    const std::vector<cdx::RecordCount> comp = {0, 100};

    const auto result = projectCov2IdxGlobal(cov, nid2flat, comp);

    EXPECT_EQ(result[0], 60);
    EXPECT_EQ(result[50], 70);
    EXPECT_EQ(result[99], 50);
}

// An empty graph (no nodes, no components) returns an empty vector.
TEST(ProjectCov2IdxGlobalTest, EmptyGraph) {
    const std::vector<cdx::Coverage> cov;
    const std::vector<cdx::FlatIdx> nid2flat;
    const std::vector<cdx::RecordCount> comp = {0, 0};

    EXPECT_TRUE(projectCov2IdxGlobal(cov, nid2flat, comp).empty());
}

// A NOT_IN_COMPO coverage *value* at a validly-mapped flat index is inert
// payload and passes through unchanged (it is numerically INVALID_NODE,
// see file header comment).
TEST(ProjectCov2IdxGlobalTest, SentinelCoverageValuePassesThroughAtValidFlatIdx) {
    const std::vector<cdx::Coverage> cov = {11, cfg::NOT_IN_COMPO, 33};
    const std::vector<cdx::FlatIdx> nid2flat = {0, 1, 2};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{11, cfg::INVALID_NODE, 33}));
}

// A node whose nid2flat_idx is INVALID_FLAT_IDX is unmapped and must be
// skipped entirely (destination slot stays at its zero default).
TEST(ProjectCov2IdxGlobalTest, InvalidFlatIdxIsSkipped) {
    const std::vector<cdx::Coverage> cov = {10, 99, 20};
    const std::vector<cdx::FlatIdx> nid2flat = {0, cfg::INVALID_FLAT_IDX, 2};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{10, 0, 20}));
}

// Zero coverage is a no-op write onto an already zero-initialized slot.
TEST(ProjectCov2IdxGlobalTest, ZeroCoverageIsSkipped) {
    const std::vector<cdx::Coverage> cov = {10, 0, 30};
    const std::vector<cdx::FlatIdx> nid2flat = {0, 1, 2};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{10, 0, 30}));
}

// All-zero coverage input.
TEST(ProjectCov2IdxGlobalTest, AllZeroCoverage) {
    const std::vector<cdx::Coverage> cov = {0, 0, 0};
    const std::vector<cdx::FlatIdx> nid2flat = {0, 1, 2};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{0, 0, 0}));
}

// Two nodes mapping to the same flat index: the later one wins.
TEST(ProjectCov2IdxGlobalTest, DuplicateFlatIdxOverwrites) {
    const std::vector<cdx::Coverage> cov = {10, 20};
    const std::vector<cdx::FlatIdx> nid2flat = {1, 1};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{0, 20, 0}));
}

// Total coverage sum is preserved across a pure permutation (no overwrite,
// no unmapped nodes).
TEST(ProjectCov2IdxGlobalTest, TotalCoverageSumPreserved) {
    const std::vector<cdx::Coverage> cov = {100, 200, 300};
    const std::vector<cdx::FlatIdx> nid2flat = {2, 0, 1};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    const auto result = projectCov2IdxGlobal(cov, nid2flat, comp);

    std::uint64_t expected = 0, actual = 0;
    for (const auto v : cov) expected += v;
    for (const auto v : result) actual += v;
    EXPECT_EQ(actual, expected);
}

// Largest valid coverage value, just below the sentinel range, survives.
TEST(ProjectCov2IdxGlobalTest, MaxValidCoverageValue) {
    const std::vector<cdx::Coverage> cov = {cfg::NOT_IN_QUERY - 1};
    const std::vector<cdx::FlatIdx> nid2flat = {0};
    const std::vector<cdx::RecordCount> comp = {0, 1};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              (std::vector<cdx::Coverage>{cfg::NOT_IN_QUERY - 1}));
}

// Output length must equal the final component offset (total node count).
TEST(ProjectCov2IdxGlobalTest, OutputLengthMatchesLastComponentOffset) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3};
    const std::vector<cdx::FlatIdx> nid2flat = {1, 2, 0};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp).size(), comp.back());
}

// Repeated calls with identical inputs must produce identical outputs.
TEST(ProjectCov2IdxGlobalTest, RepeatedExecutionIsDeterministic) {
    const std::vector<cdx::Coverage> cov = {7, 8, 9};
    const std::vector<cdx::FlatIdx> nid2flat = {2, 0, 1};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_EQ(projectCov2IdxGlobal(cov, nid2flat, comp),
              projectCov2IdxGlobal(cov, nid2flat, comp));
}

// cov_table and nid2flat_idx must share identical sizes.
TEST(ProjectCov2IdxGlobalTest, DimensionMismatchThrows) {
    const std::vector<cdx::Coverage> cov = {1, 2};
    const std::vector<cdx::FlatIdx> nid2flat = {0, 1, 2};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_THROW(projectCov2IdxGlobal(cov, nid2flat, comp), std::invalid_argument);
}

// component_offsets must contain at least two boundaries.
TEST(ProjectCov2IdxGlobalTest, ComponentOffsetsMustContainTwoBoundariesThrows) {
    const std::vector<cdx::Coverage> cov;
    const std::vector<cdx::FlatIdx> nid2flat;
    const std::vector<cdx::RecordCount> comp = {0};

    EXPECT_THROW(projectCov2IdxGlobal(cov, nid2flat, comp), std::invalid_argument);
}

// A flat_idx exceeding the global table size must raise, unconditionally.
TEST(ProjectCov2IdxGlobalTest, FlatIdxOutOfBoundsThrows) {
    const std::vector<cdx::Coverage> cov = {10};
    const std::vector<cdx::FlatIdx> nid2flat = {5};
    const std::vector<cdx::RecordCount> comp = {0, 3};

    EXPECT_THROW(projectCov2IdxGlobal(cov, nid2flat, comp), std::out_of_range);
}

// =============================================================================
// expandPosCovGlobal
//
// Concatenates per-component base-pair expansions into one graph-wide
// buffer, tracking cumulative bp offsets per component. Exercises the
// structural validation of the four offset/length invariants between
// component_offsets, idx2bp and idx2bp_offsets, plus the arithmetic
// correctness of the expansion and concatenation itself.
// =============================================================================

// Two components with the same local bp positions must still land at
// different global bp offsets.
TEST(ExpandPosCovGlobalTest, TwoComponentsSameLocalPositions) {
    const std::vector<cdx::Coverage> cov = {10, 20, 30, 40};
    const std::vector<cdx::RecordCount> comp = {0, 2, 4};
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 5, 0, 3, 4};
    const std::vector<cdx::RecordCount> posoff = {0, 3, 6};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);

    EXPECT_EQ(bp_off, (std::vector<cdx::PosBp>{0, 5, 9}));
    EXPECT_EQ(bp_cov, (std::vector<cdx::Coverage>{10, 10, 20, 20, 20, 30, 30, 30, 40}));
}

// The base pair immediately before/after a component boundary must belong
// to the correct, distinct component.
TEST(ExpandPosCovGlobalTest, ComponentBoundaryIsExact) {
    const std::vector<cdx::Coverage> cov = {10, 20};
    const std::vector<cdx::RecordCount> comp = {0, 1, 2};
    const std::vector<cdx::PosBp> idx2bp = {0, 5, 0, 3};
    const std::vector<cdx::RecordCount> posoff = {0, 2, 4};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);

    const auto split = bp_off[1];
    EXPECT_EQ(bp_cov[split - 1], 10);
    EXPECT_EQ(bp_cov[split], 20);
}

// Each component's bp offset delta must equal its own total bp length.
TEST(ExpandPosCovGlobalTest, BpComponentOffsetsConsistency) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3};
    const std::vector<cdx::RecordCount> comp = {0, 1, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 10, 0, 4, 7};
    const std::vector<cdx::RecordCount> posoff = {0, 2, 5};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);
    (void)bp_cov;

    EXPECT_EQ(bp_off, (std::vector<cdx::PosBp>{0, 10, 17}));
}

// Cross-check bp_off deltas against idx2bp's own per-component totals.
TEST(ExpandPosCovGlobalTest, ComponentBpLengthsMatchOffsets) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3, 4};
    const std::vector<cdx::RecordCount> comp = {0, 2, 4};
    const std::vector<cdx::PosBp> idx2bp = {0, 10, 25, 0, 8, 12};
    const std::vector<cdx::RecordCount> posoff = {0, 3, 6};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);
    (void)bp_cov;

    const std::size_t component_count = comp.size() - 1;
    for (std::size_t c = 0; c < component_count; ++c) {
        const auto actual_bp_len = bp_off[c + 1] - bp_off[c];
        const auto expected_bp_len = idx2bp[posoff[c + 1] - 1];
        EXPECT_EQ(actual_bp_len, expected_bp_len);
    }
}

// The concatenated output length must equal the final global bp offset.
TEST(ExpandPosCovGlobalTest, ConcatenationInvariantMultiComponents) {
    const std::vector<cdx::Coverage> cov = {5, 5, 5};
    const std::vector<cdx::RecordCount> comp = {0, 1, 2, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 10, 0, 20, 0, 30};
    const std::vector<cdx::RecordCount> posoff = {0, 2, 4, 6};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);

    EXPECT_EQ(bp_cov.size(), bp_off.back());
    EXPECT_EQ(bp_cov.size(), 60u);
}

// Total expanded bp coverage must equal sum(node_coverage * node_length).
TEST(ExpandPosCovGlobalTest, WeightedCoverageSumInvariant) {
    const std::vector<cdx::Coverage> cov = {10, 0, 25, 40};
    const std::vector<cdx::RecordCount> comp = {0, 4};
    const std::vector<cdx::PosBp> idx2bp = {0, 100, 250, 300, 1000};
    const std::vector<cdx::RecordCount> posoff = {0, 5};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);
    (void)bp_off;

    std::uint64_t expected = 0;
    for (std::size_t i = 0; i < cov.size(); ++i) {
        expected += static_cast<std::uint64_t>(cov[i]) * (idx2bp[i + 1] - idx2bp[i]);
    }
    std::uint64_t actual = 0;
    for (const auto v : bp_cov) actual += v;

    EXPECT_EQ(actual, expected);
}

// A zero-length node (idx2bp start == end) contributes no output bp.
TEST(ExpandPosCovGlobalTest, ZeroLengthNode) {
    const std::vector<cdx::Coverage> cov = {10, 99, 20};
    const std::vector<cdx::RecordCount> comp = {0, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 3, 3, 6};
    const std::vector<cdx::RecordCount> posoff = {0, 4};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);
    (void)bp_off;

    EXPECT_EQ(bp_cov, (std::vector<cdx::Coverage>{10, 10, 10, 20, 20, 20}));
}

// All-zero node coverage across the whole graph.
TEST(ExpandPosCovGlobalTest, AllZeroCoverage) {
    const std::vector<cdx::Coverage> cov = {0, 0, 0};
    const std::vector<cdx::RecordCount> comp = {0, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 5, 9};
    const std::vector<cdx::RecordCount> posoff = {0, 4};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);

    EXPECT_EQ(bp_cov, (std::vector<cdx::Coverage>(9, 0)));
    EXPECT_EQ(bp_off, (std::vector<cdx::PosBp>{0, 9}));
}

// A graph with a single component is the degenerate/base case of the
// multi-component concatenation logic.
TEST(ExpandPosCovGlobalTest, SingleComponentGraph) {
    const std::vector<cdx::Coverage> cov = {5, 10, 15};
    const std::vector<cdx::RecordCount> comp = {0, 3};
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 5, 7};
    const std::vector<cdx::RecordCount> posoff = {0, 4};

    const auto [bp_cov, bp_off] = expandPosCovGlobal(cov, comp, idx2bp, posoff);

    EXPECT_EQ(bp_off, (std::vector<cdx::PosBp>{0, 7}));
    EXPECT_EQ(bp_cov, (std::vector<cdx::Coverage>{5, 5, 10, 10, 10, 15, 15}));
}

// component_offsets must contain at least two boundaries, and match
// idx2bp_offsets in size.
TEST(ExpandPosCovGlobalTest, LessThanTwoComponentBoundariesThrows) {
    const std::vector<cdx::Coverage> cov;
    const std::vector<cdx::RecordCount> comp = {0};
    const std::vector<cdx::PosBp> idx2bp;
    const std::vector<cdx::RecordCount> posoff = {0};

    EXPECT_THROW(expandPosCovGlobal(cov, comp, idx2bp, posoff), std::invalid_argument);
}

// component_offsets and idx2bp_offsets must have identical sizes.
TEST(ExpandPosCovGlobalTest, ComponentOffsetLengthMismatchThrows) {
    const std::vector<cdx::Coverage> cov = {1, 2};
    const std::vector<cdx::RecordCount> comp = {0, 2, 4};
    const std::vector<cdx::PosBp> idx2bp = {0, 5, 10};
    const std::vector<cdx::RecordCount> posoff = {0, 3};

    EXPECT_THROW(expandPosCovGlobal(cov, comp, idx2bp, posoff), std::invalid_argument);
}

// Both offset tables must start at zero.
TEST(ExpandPosCovGlobalTest, FirstOffsetsMustBeZeroThrows) {
    const std::vector<cdx::Coverage> cov = {1};
    const std::vector<cdx::RecordCount> comp = {1, 2};
    const std::vector<cdx::PosBp> idx2bp = {0, 5};
    const std::vector<cdx::RecordCount> posoff = {0, 2};

    EXPECT_THROW(expandPosCovGlobal(cov, comp, idx2bp, posoff), std::invalid_argument);
}

// component_offsets.back() must equal flat_idx_cov_table.size().
TEST(ExpandPosCovGlobalTest, FinalComponentOffsetMismatchThrows) {
    const std::vector<cdx::Coverage> cov = {1, 2};
    const std::vector<cdx::RecordCount> comp = {0, 99};
    const std::vector<cdx::PosBp> idx2bp = {0, 5, 10};
    const std::vector<cdx::RecordCount> posoff = {0, 3};

    EXPECT_THROW(expandPosCovGlobal(cov, comp, idx2bp, posoff), std::invalid_argument);
}

// idx2bp_offsets.back() must equal idx2bp.size().
TEST(ExpandPosCovGlobalTest, FinalIdx2bpOffsetMismatchThrows) {
    const std::vector<cdx::Coverage> cov = {1};
    const std::vector<cdx::RecordCount> comp = {0, 1};
    const std::vector<cdx::PosBp> idx2bp = {0, 5};
    const std::vector<cdx::RecordCount> posoff = {0, 999};

    EXPECT_THROW(expandPosCovGlobal(cov, comp, idx2bp, posoff), std::invalid_argument);
}

// Aggregate checks (sizes, first/last totals) can all pass while the
// per-component node count implied by component_offsets still disagrees
// with what idx2bp_offsets/idx2bp actually provide for that component;
// the inner per-node bounds check must catch that and throw.
TEST(ExpandPosCovGlobalTest, InconsistentIdx2bpOffsetsThrows) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3, 4, 5};
    const std::vector<cdx::RecordCount> comp = {0, 5};   // claims 5 nodes
    const std::vector<cdx::PosBp> idx2bp = {0, 2, 4};    // only room for 2 nodes
    const std::vector<cdx::RecordCount> posoff = {0, 3};

    EXPECT_THROW(expandPosCovGlobal(cov, comp, idx2bp, posoff), std::out_of_range);
}

// =============================================================================
// trimCoverageToQuery
//
// Masks bp positions outside a query interval with cfg::NOT_IN_QUERY,
// supporting both standard [start, end] and circular/inverted ranges. The
// C++ implementation always returns a new vector and never mutates its
// input (a stronger, deliberately different contract than the Python
// reference implementation, which trims in place and returns the same
// object).
// =============================================================================

// Empty input short-circuits to an empty result without touching bounds.
TEST(TrimCoverageToQueryTest, EmptyArrayReturnsEmpty) {
    EXPECT_TRUE(trimCoverageToQuery({}, {0, 0}).empty());
}

// A query spanning the whole component leaves every value untouched.
TEST(TrimCoverageToQueryTest, FullComponentQueryLeavesArrayUnchanged) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3, 4, 5};

    EXPECT_EQ(trimCoverageToQuery(cov, {0, 4}), cov);
}

// Positions after the query end are masked.
TEST(TrimCoverageToQueryTest, TrimRightSide) {
    const std::vector<cdx::Coverage> cov = {3, 5, 6, 6, 7, 5, 3};
    const std::vector<cdx::Coverage> expected = {3, 5, 6, 6, 7, cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY};

    EXPECT_EQ(trimCoverageToQuery(cov, {0, 4}), expected);
}

// Positions before the query start are masked.
TEST(TrimCoverageToQueryTest, TrimLeftSide) {
    const std::vector<cdx::Coverage> cov = {3, 5, 6, 6, 7, 5, 3};
    const std::vector<cdx::Coverage> expected = {cfg::NOT_IN_QUERY, cfg::NOT_IN_QUERY, 6, 6, 7, 5, 3};

    EXPECT_EQ(trimCoverageToQuery(cov, {2, 6}), expected);
}

// Both prefix and suffix are masked, leaving only the query interval.
TEST(TrimCoverageToQueryTest, TrimBothSides) {
    const std::vector<cdx::Coverage> cov = {3, 5, 6, 6, 7, 5, 3};
    const cdx::Coverage n = cfg::NOT_IN_QUERY;
    const std::vector<cdx::Coverage> expected = {n, n, 6, 6, 7, n, n};

    EXPECT_EQ(trimCoverageToQuery(cov, {2, 4}), expected);
}

// Existing NOT_IN_QUERY sentinels inside the query range must be left as-is.
TEST(TrimCoverageToQueryTest, PreserveExistingNotInQuery) {
    const cdx::Coverage n = cfg::NOT_IN_QUERY;
    const std::vector<cdx::Coverage> cov = {n, n, 6, 7, n};

    EXPECT_EQ(trimCoverageToQuery(cov, {2, 3}), (std::vector<cdx::Coverage>{n, n, 6, 7, n}));
}

// Existing NOT_IN_COMPO sentinels must never be downgraded/overwritten.
TEST(TrimCoverageToQueryTest, PreserveExistingNotInCompo) {
    const cdx::Coverage c = cfg::NOT_IN_COMPO;
    const std::vector<cdx::Coverage> cov = {c, 4, 5, c};

    EXPECT_EQ(trimCoverageToQuery(cov, {1, 2}), (std::vector<cdx::Coverage>{c, 4, 5, c}));
}

// NOT_IN_COMPO sentinels *inside* the query range must survive untouched
// (only masking with NOT_IN_QUERY happens, and only outside the range).
TEST(TrimCoverageToQueryTest, PreserveNotInCompoInsideQueryRange) {
    const cdx::Coverage c = cfg::NOT_IN_COMPO;
    const cdx::Coverage n = cfg::NOT_IN_QUERY;
    const std::vector<cdx::Coverage> cov = {10, c, 20, c, 30};

    EXPECT_EQ(trimCoverageToQuery(cov, {1, 3}), (std::vector<cdx::Coverage>{n, c, 20, c, n}));
}

// A circular query (start > end) masks the middle region instead of the
// prefix/suffix.
TEST(TrimCoverageToQueryTest, CircularQuery) {
    const cdx::Coverage n = cfg::NOT_IN_QUERY;
    const std::vector<cdx::Coverage> cov = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    const std::vector<cdx::Coverage> expected = {10, 11, 12, n, n, n, n, n, 18, 19};

    EXPECT_EQ(trimCoverageToQuery(cov, {8, 2}), expected);
}

// A circular query that wraps the entire component masks nothing.
TEST(TrimCoverageToQueryTest, CircularQueryWrapsWholeComponent) {
    const std::vector<cdx::Coverage> cov = {10, 20, 30, 40, 50};

    EXPECT_EQ(trimCoverageToQuery(cov, {3, 2}), cov);
}

// A single-base query masks everything except that one position.
TEST(TrimCoverageToQueryTest, SingleBaseQuery) {
    const cdx::Coverage n = cfg::NOT_IN_QUERY;
    const std::vector<cdx::Coverage> cov = {1, 2, 3, 4, 5};

    EXPECT_EQ(trimCoverageToQuery(cov, {2, 2}), (std::vector<cdx::Coverage>{n, n, 3, n, n}));
}

// query_start/query_end at or beyond the component length must raise.
TEST(TrimCoverageToQueryTest, InvalidBoundsThrow) {
    const std::vector<cdx::Coverage> cov = {1, 2, 3};

    EXPECT_THROW(trimCoverageToQuery(cov, {0, 3}), std::out_of_range);
    EXPECT_THROW(trimCoverageToQuery(cov, {3, 1}), std::out_of_range);
}

// The C++ contract never mutates the caller's input vector, unlike the
// Python reference implementation (see suite docstring above).
TEST(TrimCoverageToQueryTest, DoesNotMutateInputArgument) {
    const std::vector<cdx::Coverage> original = {1, 2, 3};
    std::vector<cdx::Coverage> cov = original;

    const auto result = trimCoverageToQuery(cov, {1, 2});

    EXPECT_EQ(cov, original);
    EXPECT_NE(&result, &cov);
}
