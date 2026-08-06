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
#include "../src/coverage_gaps.h"
#include "../src/coverage_precision.h"
#include "cdx_types.h"

#include <vg/vg.pb.h>
#include <vg/io/protobuf_emitter.hpp>

#include <gtest/gtest.h>

#include <algorithm>
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

    // -------------------------------------------------------------------
    // Base-precision fixtures: unlike GamFileFixture above (which only sets
    // node_id, everything else defaulted), the CoveragePrecision::Base tests
    // below need full control over Position::offset/is_reverse and each
    // mapping's Edit list, so alignments are built directly with the real
    // Protobuf setters and handed in ready-made.
    // -------------------------------------------------------------------

    /** @brief Same on-disk fixture mechanics as GamFileFixture, but takes fully-built Alignments. */
    class DetailedGamFileFixture {
    public:
        explicit DetailedGamFileFixture(const std::vector<vg::Alignment> &alignments) {
            static int counter = 0;
            path_ = std::filesystem::temp_directory_path() /
                    ("gam_io_test_detailed_" + std::to_string(++counter) + ".gam");

            std::ofstream out(path_, std::ios::binary);
            {
                vg::io::ProtobufEmitter<vg::Alignment> emitter(out, false);
                for (const vg::Alignment &alignment: alignments) {
                    vg::Alignment copy = alignment;
                    emitter.write(std::move(copy));
                }
            }
        }

        ~DetailedGamFileFixture() {
            std::error_code ec;
            std::filesystem::remove(path_, ec);
        }

        DetailedGamFileFixture(const DetailedGamFileFixture &) = delete;

        DetailedGamFileFixture &operator=(const DetailedGamFileFixture &) = delete;

        [[nodiscard]] const std::filesystem::path &path() const { return path_; }

    private:
        std::filesystem::path path_;
    };

    /** @brief Appends a mapping (node_id/offset/is_reverse) to a Path and returns it for edit-building. */
    vg::Mapping *addMapping(
        vg::Path &path,
        const std::int64_t node_id,
        const std::int64_t offset = 0,
        const bool is_reverse = false
    ) {
        vg::Mapping *mapping = path.add_mapping();
        mapping->mutable_position()->set_node_id(node_id);
        mapping->mutable_position()->set_offset(offset);
        mapping->mutable_position()->set_is_reverse(is_reverse);
        return mapping;
    }

    /** @brief Appends a match edit (from_length == to_length, no sequence) - edit_is_match(). */
    void addMatchEdit(vg::Mapping &mapping, const std::int32_t length) {
        vg::Edit *edit = mapping.add_edit();
        edit->set_from_length(length);
        edit->set_to_length(length);
    }

    /** @brief Appends a deletion edit (from_length > 0, to_length == 0) - edit_is_deletion(). */
    void addDeletionEdit(vg::Mapping &mapping, const std::int32_t length) {
        vg::Edit *edit = mapping.add_edit();
        edit->set_from_length(length);
        edit->set_to_length(0);
    }

    /** @brief Appends an insertion edit (from_length == 0, to_length > 0) - edit_is_insertion(). */
    void addInsertionEdit(vg::Mapping &mapping, const std::int32_t length) {
        vg::Edit *edit = mapping.add_edit();
        edit->set_from_length(0);
        edit->set_to_length(length);
        edit->set_sequence(std::string(static_cast<std::size_t>(length), 'A'));
    }

    /** @brief Appends a substitution/mismatch edit (from_length == to_length, non-empty sequence) - edit_is_sub(). */
    void addSubstitutionEdit(vg::Mapping &mapping, const std::int32_t length) {
        vg::Edit *edit = mapping.add_edit();
        edit->set_from_length(length);
        edit->set_to_length(length);
        edit->set_sequence(std::string(static_cast<std::size_t>(length), 'A'));
    }

    /** @brief Finds the (single) gap for a given nid_offset, or nullptr if none was recorded. */
    const BpGap *findGap(const std::vector<BpGap> &gaps, const cdx::Nid nid_offset) {
        const auto it = std::find_if(gaps.begin(), gaps.end(), [nid_offset](const BpGap &g) {
            return g.nid_offset == nid_offset;
        });
        return it == gaps.end() ? nullptr : &*it;
    }
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
        const GamFileFixture fixture(std::vector<AlignmentSpec>{});
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
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{0}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.total, 1u);
        EXPECT_EQ(stats.mapped, 0u);
        EXPECT_EQ(stats.unmapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, baseTarget());
    }

    TEST(GamIoMappingFilterTest, NegativeNodeIdCountsAsUnmapped) {
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{-5}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.mapped, 0u);
        EXPECT_EQ(stats.unmapped, 1u);
    }

    // An alignment with an entirely empty Path (no mappings at all) - the
    // most common real-world representation of an unmapped read.
    TEST(GamIoMappingFilterTest, EmptyPathCountsAsUnmapped) {
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{}});
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
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{kNidMin - 1}}); // node 99, upstream of the local range
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
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{kNidMin + 5}}); // offset 5, target has only offsets 0..4
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.mapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, baseTarget());
    }

    // A node inside the local range but excluded from the active query
    // (NOT_IN_QUERY sentinel, offset 2 / node 102): mapped, but not to query.
    TEST(GamIoMappingFilterTest, NodeOutsideActiveQueryCountsMappedButNotQuery) {
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{kNidMin + 2}});
        std::vector<cdx::Coverage> target = baseTarget();

        const GamMappingStats stats = process_gam(fixture.path().string(), target, kNidMin, 100, 1, 1);

        EXPECT_EQ(stats.mapped, 1u);
        EXPECT_EQ(stats.mapped_to_query, 0u);
        EXPECT_EQ(target, baseTarget()); // untouched: still NOT_IN_QUERY, not incremented
    }

    // A node inside the active query is the only case that increments coverage.
    TEST(GamIoMappingFilterTest, NodeInsideActiveQueryIncrementsCoverage) {
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{kNidMin + 1}}); // offset 1, active
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
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{kNidMin + 0, kNidMin + 1, kNidMin + 3, kNidMin + 4}});
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
        const GamFileFixture fixture(std::vector<AlignmentSpec>{{kNidMin + 0, kNidMin + 1, kNidMin + 0}});
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
        const GamFileFixture fixture(std::vector<AlignmentSpec>{
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

// =============================================================================
// CoveragePrecision::Base - argument validation.
//
// node_lengths/out_gaps are only meaningful (and required) in Base mode; in
// Node mode (the default, exercised by every test above) they are ignored
// entirely, which is what lets Node-mode callers skip building them.
// =============================================================================
namespace {
    TEST(GamIoBasePrecisionValidationTest, MissingNodeLengthsThrows) {
        std::vector<cdx::Coverage> target(1, 0);
        std::vector<BpGap> gaps;

        EXPECT_THROW(
            process_gam(kUnusedPath.string(), target, 0, 100, 1, 1, CoveragePrecision::Base, nullptr, &gaps),
            std::invalid_argument
        );
    }

    TEST(GamIoBasePrecisionValidationTest, MissingOutGapsThrows) {
        std::vector<cdx::Coverage> target(1, 0);
        std::vector<cdx::SeqLen> node_lengths(1, 10);

        EXPECT_THROW(
            process_gam(kUnusedPath.string(), target, 0, 100, 1, 1, CoveragePrecision::Base, &node_lengths, nullptr),
            std::invalid_argument
        );
    }

    TEST(GamIoBasePrecisionValidationTest, NodeLengthsSizeMismatchThrows) {
        std::vector<cdx::Coverage> target(3, 0);
        std::vector<cdx::SeqLen> node_lengths(2, 10); // wrong size vs. target
        std::vector<BpGap> gaps;

        EXPECT_THROW(
            process_gam(kUnusedPath.string(), target, 0, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps),
            std::invalid_argument
        );
    }
} // anonymous namespace

// =============================================================================
// CoveragePrecision::Base - gap detection.
//
// Fixed setup used throughout: nid_min = 100 (so node 100 is nid_offset 0),
// a single node of length 10bp unless noted otherwise.
// =============================================================================
namespace {
    constexpr cdx::Nid kBaseNidMin = 100;

    // A deletion in the middle of an otherwise full-length mapping produces
    // exactly one gap, spanning just the deleted bases.
    TEST(GamIoBasePrecisionGapTest, MiddleDeletionProducesExactGap) {
        vg::Alignment alignment;
        vg::Mapping *mapping = addMapping(*alignment.mutable_path(), kBaseNidMin, 0, false);
        addMatchEdit(*mapping, 3);
        addDeletionEdit(*mapping, 2);
        addMatchEdit(*mapping, 5); // 3+2+5 = 10 = full node length -> no boundary gap

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target(1, 0);
        std::vector<cdx::SeqLen> node_lengths(1, 10);
        std::vector<BpGap> gaps;

        process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps);

        EXPECT_EQ(target[0], 1u); // node-level credit unaffected
        ASSERT_EQ(gaps.size(), 1u);
        EXPECT_EQ(gaps[0].nid_offset, 0u);
        EXPECT_EQ(gaps[0].range.start, 3u);
        EXPECT_EQ(gaps[0].range.end, 5u);
    }

    // Matches and mismatches/substitutions never produce gaps: this tool
    // reports raw read depth, not concordance with the reference (see the
    // discussion this feature was designed from).
    TEST(GamIoBasePrecisionGapTest, SubstitutionsAndMatchesProduceNoGaps) {
        vg::Alignment alignment;
        vg::Mapping *mapping = addMapping(*alignment.mutable_path(), kBaseNidMin, 0, false);
        addMatchEdit(*mapping, 4);
        addSubstitutionEdit(*mapping, 2);
        addMatchEdit(*mapping, 4); // 4+2+4 = 10 = full node length

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target(1, 0);
        std::vector<cdx::SeqLen> node_lengths(1, 10);
        std::vector<BpGap> gaps;

        process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps);

        EXPECT_EQ(target[0], 1u);
        EXPECT_TRUE(gaps.empty());
    }

    // An insertion consumes no "from" (node) length at all, so it can never
    // produce a gap either - there is no node position to subtract from.
    TEST(GamIoBasePrecisionGapTest, InsertionProducesNoGap) {
        vg::Alignment alignment;
        vg::Mapping *mapping = addMapping(*alignment.mutable_path(), kBaseNidMin, 0, false);
        addMatchEdit(*mapping, 5);
        addInsertionEdit(*mapping, 3); // consumes 0 "from" length
        addMatchEdit(*mapping, 5); // 5+0+5 = 10 = full node length

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target(1, 0);
        std::vector<cdx::SeqLen> node_lengths(1, 10);
        std::vector<BpGap> gaps;

        process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps);

        EXPECT_EQ(target[0], 1u);
        EXPECT_TRUE(gaps.empty());
    }

    // A single-mapping read that starts partway through its node (offset >
    // 0) and whose edits fall short of the node's far end produces two
    // boundary gaps: one before the walk starts, one after it ends.
    TEST(GamIoBasePrecisionGapTest, PartialSingleMappingProducesLeadingAndTrailingGaps) {
        vg::Alignment alignment;
        vg::Mapping *mapping = addMapping(*alignment.mutable_path(), kBaseNidMin, 3, false);
        addMatchEdit(*mapping, 5); // walk covers forward [3, 8) out of a 10bp node

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target(1, 0);
        std::vector<cdx::SeqLen> node_lengths(1, 10);
        std::vector<BpGap> gaps;

        process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps);

        EXPECT_EQ(target[0], 1u);
        ASSERT_EQ(gaps.size(), 2u);

        // Order isn't guaranteed (thread-local lists are concatenated), so
        // find each gap by its range instead of assuming positions [0]/[1].
        const bool has_leading = std::any_of(gaps.begin(), gaps.end(), [](const BpGap &g) {
            return g.range.start == 0 && g.range.end == 3;
        });
        const bool has_trailing = std::any_of(gaps.begin(), gaps.end(), [](const BpGap &g) {
            return g.range.start == 8 && g.range.end == 10;
        });
        EXPECT_TRUE(has_leading);
        EXPECT_TRUE(has_trailing);
    }

    // Reverse-strand mappings must mirror gap positions relative to the
    // node's forward orientation, not reuse the forward-strand walk order
    // directly - this is the main correctness risk in the whole feature.
    TEST(GamIoBasePrecisionGapTest, ReverseStrandDeletionIsMirrored) {
        vg::Alignment alignment;
        // is_reverse, offset=0: walk starts at forward position 9, moving down.
        vg::Mapping *mapping = addMapping(*alignment.mutable_path(), kBaseNidMin, 0, true);
        addMatchEdit(*mapping, 4); // walk positions [0,4) -> forward [6,10)
        addDeletionEdit(*mapping, 1); // walk positions [4,5) -> forward [5,6)
        addMatchEdit(*mapping, 5); // walk positions [5,10) -> forward [0,5)

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target(1, 0);
        std::vector<cdx::SeqLen> node_lengths(1, 10);
        std::vector<BpGap> gaps;

        process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps);

        ASSERT_EQ(gaps.size(), 1u);
        EXPECT_EQ(gaps[0].range.start, 5u);
        EXPECT_EQ(gaps[0].range.end, 6u);
    }

    // In a multi-node path, only the first and last mapping can produce
    // boundary gaps - an internal mapping always walks its node edge to
    // edge, by construction of how aligners emit paths.
    TEST(GamIoBasePrecisionGapTest, OnlyFirstAndLastMappingCanHaveBoundaryGaps) {
        vg::Alignment alignment;
        vg::Path &path = *alignment.mutable_path();

        // First node: starts at offset 2 (leading gap [0,2)).
        vg::Mapping *first = addMapping(path, kBaseNidMin, 2, false);
        addMatchEdit(*first, 8); // covers [2,10) of a 10bp node

        // Middle node: always offset 0, full length, no gap possible.
        vg::Mapping *middle = addMapping(path, kBaseNidMin + 1, 0, false);
        addMatchEdit(*middle, 10);

        // Last node: full offset 0 but only consumes 6 of 10bp (trailing gap [6,10)).
        vg::Mapping *last = addMapping(path, kBaseNidMin + 2, 0, false);
        addMatchEdit(*last, 6);

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target(3, 0);
        std::vector<cdx::SeqLen> node_lengths = {10, 10, 10};
        std::vector<BpGap> gaps;

        process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps);

        EXPECT_EQ(target, (std::vector<cdx::Coverage>{1, 1, 1}));

        EXPECT_EQ(findGap(gaps, 1), nullptr) << "middle node must never have a gap";

        const BpGap *first_gap = findGap(gaps, 0);
        ASSERT_NE(first_gap, nullptr);
        EXPECT_EQ(first_gap->range.start, 0u);
        EXPECT_EQ(first_gap->range.end, 2u);

        const BpGap *last_gap = findGap(gaps, 2);
        ASSERT_NE(last_gap, nullptr);
        EXPECT_EQ(last_gap->range.start, 6u);
        EXPECT_EQ(last_gap->range.end, 10u);
    }

    // A node excluded from the active query (NOT_IN_QUERY sentinel) never
    // reaches gap detection at all, even if its mapping carries a deletion -
    // the same gating process_gam already applies to the coverage increment.
    TEST(GamIoBasePrecisionGapTest, NodeOutsideActiveQueryProducesNoGap) {
        vg::Alignment alignment;
        vg::Mapping *mapping = addMapping(*alignment.mutable_path(), kBaseNidMin, 0, false);
        addDeletionEdit(*mapping, 5);
        addMatchEdit(*mapping, 5);

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target = {cfg::NOT_IN_QUERY};
        std::vector<cdx::SeqLen> node_lengths(1, 10);
        std::vector<BpGap> gaps;

        process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Base, &node_lengths, &gaps);

        EXPECT_TRUE(gaps.empty());
    }

    // Node mode (the default used throughout this file) must ignore Edit
    // data entirely: a deletion in the fixture has zero effect on the
    // node-level result, matching plain node-presence counting.
    TEST(GamIoBasePrecisionGapTest, NodeModeIgnoresEditsEvenWhenPresent) {
        vg::Alignment alignment;
        vg::Mapping *mapping = addMapping(*alignment.mutable_path(), kBaseNidMin, 3, true);
        addDeletionEdit(*mapping, 4);

        const DetailedGamFileFixture fixture({alignment});
        std::vector<cdx::Coverage> target(1, 0);

        // Explicit CoveragePrecision::Node with null node_lengths/out_gaps -
        // exactly what a Node-mode caller does, must not throw or crash
        // despite the fixture carrying real edit/offset/strand data.
        const GamMappingStats stats =
                process_gam(fixture.path().string(), target, kBaseNidMin, 100, 1, 1, CoveragePrecision::Node);

        EXPECT_EQ(stats.mapped_to_query, 1u);
        EXPECT_EQ(target[0], 1u);
    }
} // anonymous namespace
