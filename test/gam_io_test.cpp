/**
 * @file gam_io_test.cpp
 * @brief Unit tests for gam_io.cpp (process_gam).
 *
 * process_gam streams real vg::Alignment Protobuf messages framed as a GAM file (via
 * vg::io::MessageIterator) and accumulates per-node coverage in parallel with OpenMP. Fixtures
 * here are built with the production serialization stack itself (vg::io::ProtobufEmitter<Alignment>,
 * which uses the same Registry-assigned "GAM" tag process_gam checks against) rather than hand-rolled
 * bytes, so they are indistinguishable from what a real aligner would produce.
 *
 * Only `path.mapping[].position.node_id` is read by process_gam; every other Alignment field is left
 * at its default value in these fixtures.
 *
 * NOTE ON VERIFICATION: unlike the other modules in this test suite, this file could not be locally
 * compiled/run in the sandbox used to author it (no cmake, no protobuf, no htslib available there). It
 * was written from a careful line-by-line reading of gam_io.cpp, protobuf_emitter.hpp, and
 * message_iterator.hpp rather than the usual compile+mutate verification loop - please build and run
 * it (ctest / CLion) and report back anything that doesn't match.
 */

#include "../src/gam_io.h"
#include "../src/config.h"
#include "cdx_types.h"

#include <vg/vg.pb.h>
#include <vg/io/protobuf_emitter.hpp>

#include <gtest/gtest.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace {
    // =====================================================================
    // Test infrastructure: builds real, valid .gam fixture files on disk
    // using vg::io's own serialization stack.
    // =====================================================================

    /**
     * @brief Describes one alignment purely by the node IDs its path visits,
     *        in order. An empty list represents an unmapped read (no Path
     *        mappings at all) - the most common real-world "unmapped" shape.
     */
    using AlignmentSpec = std::vector<std::int64_t>;

    class GamFileFixture {
    public:
        explicit GamFileFixture(const std::vector<AlignmentSpec> &alignments) {
            static int counter = 0;
            path_ = std::filesystem::temp_directory_path() /
                    ("gam_io_test_" + std::to_string(++counter) + ".gam");

            std::ofstream out(path_, std::ios::binary);
            {
                // Uncompressed (compress=false): simpler and faster for small
                // fixtures; process_gam's MessageIterator handles both.
                vg::io::ProtobufEmitter<vg::Alignment> emitter(out, false);

                for (const auto &node_ids: alignments) {
                    vg::Alignment alignment;
                    for (const std::int64_t node_id: node_ids) {
                        alignment.mutable_path()->add_mapping()->mutable_position()->set_node_id(node_id);
                    }
                    emitter.write(std::move(alignment));
                }
                // emitter destructor flushes the final group here, before `out` closes below.
            }
        }

        ~GamFileFixture() {
            std::error_code ec;
            std::filesystem::remove(path_, ec);
        }

        GamFileFixture(const GamFileFixture &) = delete;

        GamFileFixture &operator=(const GamFileFixture &) = delete;

        [[nodiscard]] const std::filesystem::path &path() const { return path_; }

    private:
        std::filesystem::path path_;
    };
} // anonymous namespace

// =============================================================================
// Argument validation.
//
// All four checks happen before the GAM file is even opened, so a
// placeholder (non-existent) path is used throughout this section.
// =============================================================================
namespace {
    const std::filesystem::path kUnusedPath = "/nonexistent/path/should/not/be/opened.gam";

    TEST(GamIoArgumentValidationTest, BatchSizeZeroThrows) {
        std::vector<cdx::Coverage> target(1, 0);
        EXPECT_THROW(process_gam(kUnusedPath.string(), target, 0, 0, 1, 1), std::invalid_argument);
    }

    TEST(GamIoArgumentValidationTest, DecompressionThreadsZeroThrows) {
        std::vector<cdx::Coverage> target(1, 0);
        EXPECT_THROW(process_gam(kUnusedPath.string(), target, 0, 100, 0, 1), std::invalid_argument);
    }

    TEST(GamIoArgumentValidationTest, DecompressionThreadsNegativeThrows) {
        std::vector<cdx::Coverage> target(1, 0);
        EXPECT_THROW(process_gam(kUnusedPath.string(), target, 0, 100, -1, 1), std::invalid_argument);
    }

    // worker_threads<=0 previously went unvalidated and fed straight into
    // std::vector sizing / #pragma omp parallel num_threads(...); now rejected
    // up front just like decompression_threads.
    TEST(GamIoArgumentValidationTest, WorkerThreadsZeroThrows) {
        std::vector<cdx::Coverage> target(1, 0);
        EXPECT_THROW(process_gam(kUnusedPath.string(), target, 0, 100, 1, 0), std::invalid_argument);
    }

    TEST(GamIoArgumentValidationTest, WorkerThreadsNegativeThrows) {
        std::vector<cdx::Coverage> target(1, 0);
        EXPECT_THROW(process_gam(kUnusedPath.string(), target, 0, 100, 1, -4), std::invalid_argument);
    }

    TEST(GamIoArgumentValidationTest, EmptyTargetVectorThrows) {
        std::vector<cdx::Coverage> target; // empty
        EXPECT_THROW(process_gam(kUnusedPath.string(), target, 0, 100, 1, 1), std::invalid_argument);
    }
} // anonymous namespace

// =============================================================================
// File handling.
// =============================================================================
namespace {
    TEST(GamIoFileHandlingTest, MissingFileThrowsRuntimeError) {
        const std::filesystem::path missing = std::filesystem::temp_directory_path() / "does_not_exist_gam_io.gam";
        std::vector<cdx::Coverage> target(4, 0);

        EXPECT_THROW(process_gam(missing.string(), target, 100, 100, 1, 1), std::runtime_error);
    }

    // A syntactically valid GAM file with zero alignments (ProtobufEmitter
    // always writes at least the type tag) must produce all-zero stats.
    //
    // IMPORTANT: process_gam does NOT accumulate onto whatever was already in
    // `target` - the final parallel-reduction step unconditionally overwrites
    // every valid (non-sentinel) slot with the total computed for *this*
    // call (`target[node_offset] = total_coverage;`, not `+=`). With zero
    // alignments that total is 0, so a valid slot that started at a nonzero
    // value (5 here) still gets reset to 0. Only the two sentinel slots
    // (NOT_IN_QUERY / NOT_IN_COMPO) are left completely untouched, since
    // valid_nodes[] is false for them and the write is skipped entirely.
    TEST(GamIoFileHandlingTest, EmptyGamFileYieldsZeroStatsAndResetsValidSlotsOnly) {
        const GamFileFixture fixture({});
        std::vector<cdx::Coverage> target{0, cfg::NOT_IN_QUERY, 5, cfg::NOT_IN_COMPO};

        const GamMappingStats stats = process_gam(fixture.path().string(), target, 100, 100, 1, 1);

        EXPECT_EQ(stats.total, 0u);
        EXPECT_EQ(stats.mapped, 0u);
        EXPECT_EQ(stats.unmapped, 0u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, (std::vector<cdx::Coverage>{0, cfg::NOT_IN_QUERY, 0, cfg::NOT_IN_COMPO}));
    }
} // anonymous namespace

// =============================================================================
// Per-mapping filtering logic.
//
// Fixed setup used throughout: nid_min = 100, a 5-node local coverage range
// (nodes 100..104, offsets 0..4), with node offset 2 (node 102) marked
// out-of-query via the NOT_IN_QUERY sentinel and every other offset active
// (starting coverage 0).
// =============================================================================
namespace {
    std::vector<cdx::Coverage> baseTarget() {
        return {0, 0, cfg::NOT_IN_QUERY, 0, 0};
    }

    constexpr cdx::Nid kNidMin = 100;

    // A read whose only mapping has a non-positive node_id is never marked
    // mapped at all (the `continue` happens before read_is_mapped is set).
    TEST(GamIoMappingFilterTest, NonPositiveNodeIdCountsAsUnmapped) {
        const GamFileFixture fixture({{0}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.total, 1u);
        EXPECT_EQ(stats.mapped, 0u);
        EXPECT_EQ(stats.unmapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, baseTarget());
    }

    TEST(GamIoMappingFilterTest, NegativeNodeIdCountsAsUnmapped) {
        const GamFileFixture fixture({{-5}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.mapped, 0u);
        EXPECT_EQ(stats.unmapped, 1u);
    }

    // An alignment with an entirely empty Path (no mappings at all) - the
    // most common real-world representation of an unmapped read.
    TEST(GamIoMappingFilterTest, EmptyPathCountsAsUnmapped) {
        const GamFileFixture fixture({{}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.total, 1u);
        EXPECT_EQ(stats.mapped, 0u);
        EXPECT_EQ(stats.unmapped, 1u);
        EXPECT_EQ(target, baseTarget());
    }

    // A positive node_id below nid_min still marks the read "mapped" (the
    // flag is set before the range check), but must not touch coverage or
    // count toward mapped_to_query.
    TEST(GamIoMappingFilterTest, NodeIdBelowNidMinCountsMappedButNotQuery) {
        const GamFileFixture fixture({{kNidMin - 1}}); // node 99, upstream of the local range
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.mapped, 1u);
        EXPECT_EQ(stats.unmapped, 0u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, baseTarget());
    }

    // A positive node_id at/beyond nid_min + coverage_size is downstream of
    // the local range: same "mapped but not query" outcome.
    TEST(GamIoMappingFilterTest, NodeIdBeyondCoverageSizeCountsMappedButNotQuery) {
        const GamFileFixture fixture({{kNidMin + 5}}); // offset 5, target has only offsets 0..4
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.mapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, baseTarget());
    }

    // A node inside the local range but excluded from the active query
    // (NOT_IN_QUERY sentinel, offset 2 / node 102): mapped, but not to query.
    TEST(GamIoMappingFilterTest, NodeOutsideActiveQueryCountsMappedButNotQuery) {
        const GamFileFixture fixture({{kNidMin + 2}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.mapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, baseTarget()); // untouched: still NOT_IN_QUERY, not incremented
    }

    // A node inside the active query is the only case that increments coverage.
    TEST(GamIoMappingFilterTest, NodeInsideActiveQueryIncrementsCoverage) {
        const GamFileFixture fixture({{kNidMin + 1}}); // offset 1, active
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.total, 1u);
        EXPECT_EQ(stats.mapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 1u);
        EXPECT_EQ(target[1], 1u);
        EXPECT_EQ(target[0], 0u);
        EXPECT_EQ(target[3], 0u);
        EXPECT_EQ(target[4], 0u);
    }

    // Every mapping in a single alignment's path contributes its own
    // coverage increment, but the read is only counted once in the stats.
    TEST(GamIoMappingFilterTest, MultipleMappingsEachIncrementCoverageOnceReadCountedOnce) {
        const GamFileFixture fixture({{kNidMin + 0, kNidMin + 1, kNidMin + 3, kNidMin + 4}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.total, 1u);
        EXPECT_EQ(stats.mapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 1u);
        EXPECT_EQ(target, (std::vector<cdx::Coverage>{1, 1, cfg::NOT_IN_QUERY, 1, 1}));
    }

    // A path that revisits the same node twice (e.g. a loop) accumulates
    // coverage twice from a single read.
    TEST(GamIoMappingFilterTest, RevisitingSameNodeAccumulatesCoverageTwice) {
        const GamFileFixture fixture({{kNidMin + 0, kNidMin + 1, kNidMin + 0}});
        std::vector<cdx::Coverage> target = baseTarget();

        process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(target[0], 2u);
        EXPECT_EQ(target[1], 1u);
    }
} // anonymous namespace

// =============================================================================
// Aggregate stats bookkeeping across a mixed batch of reads.
// =============================================================================
namespace {
    TEST(GamIoStatsBookkeepingTest, MixedBatchProducesExactCounts) {
        const GamFileFixture fixture({
            {}, // unmapped (empty path)
            {0}, // unmapped (non-positive node_id)
            {kNidMin - 1}, // mapped, not query (upstream)
            {kNidMin + 5}, // mapped, not query (downstream)
            {kNidMin + 2}, // mapped, not query (excluded by sentinel)
            {kNidMin + 1}, // mapped, to query
            {kNidMin + 3, kNidMin + 4}, // mapped, to query (2 mappings, 1 read)
        });
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.total, 7u);
        EXPECT_EQ(stats.unmapped, 2u);
        EXPECT_EQ(stats.mapped, 5u);
        EXPECT_EQ(stats.mapped_to_query, 2u);
        EXPECT_EQ(target, (std::vector<cdx::Coverage>{0, 1, cfg::NOT_IN_QUERY, 1, 1}));
    }
} // anonymous namespace

// =============================================================================
// Parallel reduction correctness: the final coverage/stats must not depend on
// how many worker threads (nor how the stream is chopped into batches) were
// used, since accumulation is a commutative sum across threads/batches.
// =============================================================================
namespace {
    std::vector<AlignmentSpec> manyAlignmentsOverTenNodes() {
        // 200 reads, each hitting 1-3 nodes within offsets [0, 9] (nid_min..nid_min+9).
        std::vector<AlignmentSpec> alignments;
        alignments.reserve(200);
        for (int i = 0; i < 200; ++i) {
            AlignmentSpec spec;
            const int base = i % 8; // keeps every read fully inside [0, 9]
            spec.push_back(kNidMin + base);
            spec.push_back(kNidMin + base + 1);
            if (i % 3 == 0) {
                spec.push_back(kNidMin + base + 2);
            }
            alignments.push_back(spec);
        }
        return alignments;
    }

    TEST(GamIoThreadingTest, MultipleWorkerThreadsMatchSingleThreadedResult) {
        const GamFileFixture fixture(manyAlignmentsOverTenNodes());

        std::vector<cdx::Coverage> target_single(10, 0);
        const GamMappingStats stats_single = process_gam(fixture.path().string(), target_single, kNidMin, 7, 1, 1);

        std::vector<cdx::Coverage> target_multi(10, 0);
        const GamMappingStats stats_multi = process_gam(fixture.path().string(), target_multi, kNidMin, 7, 1, 4);

        EXPECT_EQ(target_single, target_multi);
        EXPECT_EQ(stats_single.total, stats_multi.total);
        EXPECT_EQ(stats_single.mapped, stats_multi.mapped);
        EXPECT_EQ(stats_single.unmapped, stats_multi.unmapped);
        EXPECT_EQ(stats_single.mapped_to_query, stats_multi.mapped_to_query);
        EXPECT_EQ(stats_single.total, 200u);
    }
} // anonymous namespace
