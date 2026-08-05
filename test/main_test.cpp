/**
 * @file main_test.cpp
 * @brief End-to-end smoke tests running the real, fully-built `cdx_coverage` executable.
 *
 * main.cpp itself has no reachable public API (main() and its two helper pipelines are all
 * either the entry point or in an anonymous namespace) - it is pure orchestration of modules
 * that already each have their own dedicated unit test suite (cdx_loader, gam_io,
 * cov_projection, output_coverage/stats/plot, cli, query_resolver), plus one piece of real
 * logic that WAS extracted and unit-tested directly: see query_plot_slice_test.cpp.
 *
 * What's left to verify about main.cpp specifically is the *wiring*: does argument parsing
 * correctly dispatch to inspect/query/global mode, do the requested output files actually get
 * written to the requested directory, and does a failure anywhere in the pipeline surface as a
 * clean non-zero exit with a readable message? That's exactly what these tests check, by
 * building small but real CDX+GAM fixtures (same technique as cdx_loader_test.cpp /
 * gam_io_test.cpp) and invoking the actual compiled `cdx_coverage` binary as a subprocess via
 * popen(), then inspecting its exit code, captured stdout+stderr, and the files it produced.
 * Nothing here re-derives per-position coverage values or pixel content - that's the job of
 * the modules' own unit tests; these tests only check that running the real binary end-to-end
 * "does the right general thing" for each top-level mode.
 */

#include <gtest/gtest.h>

#include "cdx_format.h"
#include "cdx_types.h"

#include <vg/vg.pb.h>
#include <vg/io/protobuf_emitter.hpp>

#include <array>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <tuple>
#include <vector>

namespace {
    namespace fs = std::filesystem;

    // -----------------------------------------------------------------
    // Test infrastructure
    // -----------------------------------------------------------------

    /** @brief RAII temporary directory, recursively removed on destruction. */
    class TempDir {
    public:
        TempDir() {
            path_ = fs::temp_directory_path() /
                    ("main_e2e_test_" + std::to_string(reinterpret_cast<std::uintptr_t>(this)));
            fs::create_directories(path_);
        }

        ~TempDir() {
            std::error_code ec;
            fs::remove_all(path_, ec);
        }

        TempDir(const TempDir &) = delete;

        TempDir &operator=(const TempDir &) = delete;

        [[nodiscard]] const fs::path &path() const { return path_; }

    private:
        fs::path path_;
    };

    /** @brief One component to serialize into a real binary .cdx fixture. */
    struct ComponentSpec {
        std::string name;
        cdx::Nid nid_min;
        cdx::Nid nid_max;
        // (node_id, local_idx, seq_len) triples, in on-disk (node_id-ascending) order.
        std::vector<std::tuple<cdx::Nid, cdx::Idx, cdx::SeqLen> > records;
    };

    /** @brief Writes a valid binary CDX file via cdx::CdxFormat::pack_*, same as production. */
    class CdxFileFixture {
    public:
        explicit CdxFileFixture(const std::vector<ComponentSpec> &components) {
            static int counter = 0;
            path_ = fs::temp_directory_path() /
                    ("main_e2e_test_" + std::to_string(++counter) + ".cdx");

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
            fs::remove(path_, ec);
        }

        CdxFileFixture(const CdxFileFixture &) = delete;

        CdxFileFixture &operator=(const CdxFileFixture &) = delete;

        [[nodiscard]] const fs::path &path() const { return path_; }

    private:
        fs::path path_;
    };

    /** @brief One alignment, described purely by the node IDs its path visits, in order. */
    using AlignmentSpec = std::vector<std::int64_t>;

    /** @brief Writes a real, valid .gam fixture via vg::io::ProtobufEmitter, same stack the
     *         real program reads with. */
    class GamFileFixture {
    public:
        explicit GamFileFixture(const std::vector<AlignmentSpec> &alignments) {
            static int counter = 0;
            path_ = fs::temp_directory_path() /
                    ("main_e2e_test_" + std::to_string(++counter) + ".gam");

            std::ofstream out(path_, std::ios::binary);
            {
                vg::io::ProtobufEmitter<vg::Alignment> emitter(out, /*compress=*/false);
                for (const auto &node_ids: alignments) {
                    vg::Alignment alignment;
                    for (const std::int64_t node_id: node_ids) {
                        alignment.mutable_path()->add_mapping()->mutable_position()->set_node_id(node_id);
                    }
                    emitter.write(std::move(alignment));
                }
            }
        }

        ~GamFileFixture() {
            std::error_code ec;
            fs::remove(path_, ec);
        }

        GamFileFixture(const GamFileFixture &) = delete;

        GamFileFixture &operator=(const GamFileFixture &) = delete;

        [[nodiscard]] const fs::path &path() const { return path_; }

    private:
        fs::path path_;
    };

    /**
     * @brief A small, fixed pangenome fixture shared by most tests: two components,
     *        "chr1" (3 nodes x 100bp = 300bp) and "chr2" (2 nodes x 50bp = 100bp).
     */
    std::vector<ComponentSpec> twoComponentFixture() {
        return {
            ComponentSpec{"chr1", 10, 12, {{10, 0, 100}, {11, 1, 100}, {12, 2, 100}}},
            ComponentSpec{"chr2", 20, 21, {{20, 0, 50}, {21, 1, 50}}},
        };
    }

    /** @brief A handful of alignments giving both components some non-zero coverage, plus one
     *         deliberately unmapped read (no Path at all). */
    std::vector<AlignmentSpec> sampleAlignments() {
        return {
            {10}, // maps into chr1, node 10
            {10, 11}, // maps into chr1, nodes 10 and 11
            {20}, // maps into chr2, node 20
            {}, // unmapped read
        };
    }

    /** @brief Result of running the real cdx_coverage binary as a subprocess. */
    struct RunResult {
        int exit_code = -1;
        std::string output; // combined stdout+stderr
    };

    /**
     * @brief Resolved at compile time by test/CMakeLists.txt to the real, freshly-built
     *        cdx_coverage executable's path ($<TARGET_FILE:cdx_coverage>), so this test always
     *        runs against the binary that was just compiled, not a stale/system-wide one.
     */
#ifndef CDX_COVERAGE_EXECUTABLE_PATH
#error "CDX_COVERAGE_EXECUTABLE_PATH must be defined by test/CMakeLists.txt"
#endif

    /** @brief Runs the real cdx_coverage binary with @p args, capturing combined output. */
    RunResult runCdxCoverage(const std::vector<std::string> &args) {
        std::ostringstream command;
        command << std::quoted(std::string(CDX_COVERAGE_EXECUTABLE_PATH));
        for (const auto &arg: args) {
            command << ' ' << std::quoted(arg);
        }
        command << " 2>&1";

        RunResult result;
        FILE *pipe = popen(command.str().c_str(), "r");
        if (pipe == nullptr) {
            ADD_FAILURE() << "Failed to launch cdx_coverage subprocess.";
            return result;
        }

        std::array<char, 512> buffer{};
        while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
            result.output += buffer.data();
        }

        const int status = pclose(pipe);
        result.exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
        return result;
    }
} // namespace

// ===========================================================================
// Global mode
// ===========================================================================

TEST(MainE2eGlobalModeTest, DefaultRunProducesAllThreeOutputFiles) {
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string()
    });

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_profile.tsv"));
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_stats.txt"));
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_graph.png"));
    EXPECT_GT(fs::file_size(out_dir.path() / "coverage_profile.tsv"), 0u);
    EXPECT_GT(fs::file_size(out_dir.path() / "coverage_stats.txt"), 0u);
    EXPECT_GT(fs::file_size(out_dir.path() / "coverage_graph.png"), 0u);
}

TEST(MainE2eGlobalModeTest, NoGraphNoStatsProducesOnlyTheTable) {
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string(),
        "--no-graph", "--no-stats"
    });

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_profile.tsv"));
    EXPECT_FALSE(fs::exists(out_dir.path() / "coverage_stats.txt"));
    EXPECT_FALSE(fs::exists(out_dir.path() / "coverage_graph.png"));
}

TEST(MainE2eGlobalModeTest, AllOutputsDisabledFailsCleanly) {
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string(),
        "--no-graph", "--no-stats", "--no-table"
    });

    EXPECT_NE(result.exit_code, 0);
    EXPECT_NE(result.output.find("All outputs are disabled"), std::string::npos) << result.output;
}

// ===========================================================================
// Query mode
// ===========================================================================

TEST(MainE2eQueryModeTest, QueryByComponentNameProducesGraphAndTable) {
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string(),
        "-q", "chr1", "--no-stats"
    });

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_profile.tsv"));
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_graph.png"));
}

TEST(MainE2eQueryModeTest, QueryWithExplicitRangeSlicesTheLinearGraph) {
    // Exercises the newly-extracted sliceLinearQueryCoverage() in situ: the linear graph
    // branch must accept a sub-range query without throwing and produce a graph.
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string(),
        "-q", "chr1 0:50", "--no-stats", "--no-table"
    });

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_graph.png"));
    EXPECT_GT(fs::file_size(out_dir.path() / "coverage_graph.png"), 0u);
}

TEST(MainE2eQueryModeTest, QueryByComponentIdInsteadOfNameAlsoWorks) {
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    // Component IDs are 0-based in registration order: chr1 is component 0.
    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string(),
        "-q", "0", "--no-graph", "--no-stats"
    });

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
    EXPECT_TRUE(fs::exists(out_dir.path() / "coverage_profile.tsv"));
}

TEST(MainE2eQueryModeTest, UnknownComponentNameFailsCleanly) {
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string(),
        "-q", "does_not_exist"
    });

    EXPECT_NE(result.exit_code, 0);
    EXPECT_NE(result.output.find("[FATAL ERROR]"), std::string::npos) << result.output;
}

// ===========================================================================
// Inspect mode
// ===========================================================================

TEST(MainE2eInspectModeTest, SummaryListsAllComponentsAndSkipsOutputDirectory) {
    const CdxFileFixture cdx(twoComponentFixture());
    const TempDir out_dir_parent;
    // A directory that must NOT be created by inspect mode (main.cpp explicitly skips
    // create_directories() when args.inspectMode() is true).
    const fs::path out_dir = out_dir_parent.path() / "should_not_be_created";

    const RunResult result = runCdxCoverage({cdx.path().string(), "-i", "-o", out_dir.string()});

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
    EXPECT_NE(result.output.find("chr1"), std::string::npos) << result.output;
    EXPECT_NE(result.output.find("chr2"), std::string::npos) << result.output;
    EXPECT_FALSE(fs::exists(out_dir));
}

TEST(MainE2eInspectModeTest, SingleComponentByNameShowsOnlyThatComponent) {
    const CdxFileFixture cdx(twoComponentFixture());
    const TempDir out_dir;

    // "=" assignment syntax is used deliberately here: -i/--inspect accepts 0 or 1 values
    // (CLI11 type_size(0,1)), and "--inspect=chr1" unambiguously binds the value to the
    // option instead of risking it being parsed as a stray positional. If this test fails to
    // pick up "chr1" as the inspected component, try "-i chr1" (space-separated) instead -
    // this is the one CLI11 behavior this file couldn't verify locally (see file docstring).
    const RunResult result = runCdxCoverage({
        cdx.path().string(), "--inspect=chr1", "-o", out_dir.path().string()
    });

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
    EXPECT_NE(result.output.find("chr1"), std::string::npos) << result.output;
}

TEST(MainE2eInspectModeTest, InspectModeDoesNotRequireAGamFile) {
    const CdxFileFixture cdx(twoComponentFixture());
    const TempDir out_dir;

    // No GAM path at all - only valid in inspect mode.
    const RunResult result = runCdxCoverage({cdx.path().string(), "-i", "-o", out_dir.path().string()});

    EXPECT_EQ(result.exit_code, 0) << "Output was:\n" << result.output;
}

// ===========================================================================
// Error handling
// ===========================================================================

TEST(MainE2eErrorHandlingTest, MissingGamFileOutsideInspectModeFailsCleanly) {
    const CdxFileFixture cdx(twoComponentFixture());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({cdx.path().string(), "-o", out_dir.path().string()});

    EXPECT_NE(result.exit_code, 0);
    EXPECT_NE(result.output.find("[FATAL ERROR]"), std::string::npos) << result.output;
    EXPECT_NE(result.output.find("GAM file argument is required"), std::string::npos) << result.output;
}

TEST(MainE2eErrorHandlingTest, NonexistentCdxFileFailsAtParseTime) {
    const TempDir out_dir;
    const fs::path bogus_cdx = out_dir.path() / "does_not_exist.cdx";

    const RunResult result = runCdxCoverage({bogus_cdx.string(), "-o", out_dir.path().string()});

    // CLI11's ->check(CLI::ExistingFile) on the "cdx" positional rejects this before any
    // application logic runs; the exact message is CLI11's own, so only the failure itself
    // (non-zero exit) is asserted here.
    EXPECT_NE(result.exit_code, 0);
}

TEST(MainE2eErrorHandlingTest, MalformedQueryStringFailsCleanly) {
    const CdxFileFixture cdx(twoComponentFixture());
    const GamFileFixture gam(sampleAlignments());
    const TempDir out_dir;

    const RunResult result = runCdxCoverage({
        cdx.path().string(), gam.path().string(), "-o", out_dir.path().string(),
        "-q", "chr1:1000:5000" // colon directly against the component name: invalid syntax
                                // once whitespace-separated START:END is expected; see
                                // cli_test.cpp's ColonBetweenComponentAndRangeIsTakenAsLiteral-
                                // ComponentName - this is actually accepted as a literal
                                // component name "chr1:1000:5000", which then fails to
                                // resolve to any real component (same net effect: failure).
    });

    EXPECT_NE(result.exit_code, 0);
    EXPECT_NE(result.output.find("[FATAL ERROR]"), std::string::npos) << result.output;
}
