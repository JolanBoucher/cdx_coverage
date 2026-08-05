/**
 * @file circular_plot_test.cpp
 * @brief Light integration tests for python_circular_plot.cpp (C++ side of circular plots).
 *
 * python_circular_plot.cpp's internals (prepareCircularPlotPackage, BinaryRequestWriter,
 * resolveCircularPlotScript, resolvePythonExecutable, runCircularPlotScriptOrThrow, ...) are
 * all in an anonymous namespace, so only the two public entry points declared in
 * output_plot.h - writeCircularPlotQuery() and writeCircularPlotGlobal() - are reachable.
 *
 * The real rendering backend (python_script/circular_plot.py) depends on pycirclize, an
 * optional runtime-only dependency (see CIRCULAR_PLOT_SETUP.md) that may not be installed on
 * every machine that builds/tests this project. Rather than skip most of this module's C++
 * orchestration logic for that reason, most tests below point the CDX_CIRCULAR_PLOT_SCRIPT
 * environment variable at tiny *stub* Python scripts (written by this test file itself) that
 * either always succeed (writing a trivially valid PNG, ignoring the actual request content)
 * or always fail (non-zero exit + a distinctive stderr marker). Since the C++ side only cares
 * about (a) the subprocess's exit status and (b) whether the expected PNG file exists
 * afterwards, this lets the *entire* C++ pipeline - temp directory creation, traversal
 * extraction (including the tricky origin-crossing case), sentinel masking, binary request
 * serialization, real subprocess invocation via popen/pclose, stdout/stderr capture, error
 * propagation, and work-dir cleanup - be verified end-to-end with a real `python3` process,
 * without ever requiring pycirclize/matplotlib/numpy to be installed for the *stub*. (The
 * *real* circular_plot.py script has its own independent test suite:
 * test/circular_plot_py_test.py.)
 *
 * Test categories:
 *   A - Pure C++ validation, no Python involved at all (errors thrown before any subprocess
 *       is ever considered).
 *   C - "Stub always succeeds" scripts: full happy-path orchestration, for all three query
 *       traversal topologies (full component / contiguous sub-range / origin-crossing), plus
 *       logarithmic scale, global mode, and directory auto-creation.
 *   D - "Stub always fails" script: verifies runtime_error propagation, including the
 *       captured stderr diagnostic text, and that the work directory survives a failed run
 *       (left behind for debugging - cleanupWorkDir() is only reached on success).
 *   E - Script-resolution failure: no CDX_CIRCULAR_PLOT_SCRIPT override, and no
 *       circular_plot.py sitting next to this test binary (this test target is deliberately
 *       NOT compiled with CDX_CIRCULAR_PLOT_SCRIPT_SOURCE_DIR), so resolution must fail with a
 *       clear message. Also verifies a *nonexistent* override path is correctly ignored
 *       (falls through) rather than being blindly trusted.
 *   F - CDX_PYTHON_EXECUTABLE override plumbing: pointing it at a bogus interpreter path must
 *       make the run fail, proving the override is genuinely used (not silently ignored).
 *
 * Categories C/D/F need a real `python3` on PATH; each self-skips via GTEST_SKIP() if none is
 * found, rather than failing.
 */

#include <gtest/gtest.h>

#include "../src/output_plot.h"
#include "../src/config.h"

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {
    namespace fs = std::filesystem;

    // -----------------------------------------------------------------
    // Shared test infrastructure
    // -----------------------------------------------------------------

    /** @brief RAII temporary directory, recursively removed on destruction. */
    class TempDir {
    public:
        TempDir() {
            path_ = fs::temp_directory_path() /
                    ("circular_plot_test_" + std::to_string(reinterpret_cast<std::uintptr_t>(this)));
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

    /**
     * @brief RAII environment variable override, restoring whatever value (or absence) the
     *        variable had beforehand. Needed because gtest binaries invoked without
     *        --gtest_filter run every test in one process, so env var state must not leak
     *        across tests even though ctest normally isolates each case into its own process.
     */
    class EnvVarGuard {
    public:
        /** @brief Sets @p name to @p value. */
        EnvVarGuard(std::string name, const std::string &value) : name_(std::move(name)) {
            capturePrevious();
            ::setenv(name_.c_str(), value.c_str(), /*overwrite=*/1);
        }

        /** @brief Forcibly unsets @p name for the guard's lifetime. */
        explicit EnvVarGuard(std::string name) : name_(std::move(name)) {
            capturePrevious();
            ::unsetenv(name_.c_str());
        }

        ~EnvVarGuard() {
            if (had_previous_) {
                ::setenv(name_.c_str(), previous_value_.c_str(), 1);
            } else {
                ::unsetenv(name_.c_str());
            }
        }

        EnvVarGuard(const EnvVarGuard &) = delete;

        EnvVarGuard &operator=(const EnvVarGuard &) = delete;

    private:
        void capturePrevious() {
            if (const char *old = std::getenv(name_.c_str())) {
                had_previous_ = true;
                previous_value_ = old;
            }
        }

        std::string name_;
        bool had_previous_ = false;
        std::string previous_value_;
    };

    /** @brief Writes @p content to @p path as a plain text file (used for the stub scripts). */
    void writeTextFile(const fs::path &path, const std::string &content) {
        std::ofstream out(path, std::ios::binary);
        out << content;
    }

    /** @brief True if @p path starts with a well-formed 8-byte PNG signature. */
    bool hasValidPngSignature(const fs::path &path) {
        std::ifstream file(path, std::ios::binary);
        if (!file.is_open()) return false;
        unsigned char sig[8] = {};
        file.read(reinterpret_cast<char *>(sig), sizeof(sig));
        if (file.gcount() != static_cast<std::streamsize>(sizeof(sig))) return false;
        static constexpr unsigned char kPngSignature[8] = {
            0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A
        };
        for (int i = 0; i < 8; ++i) {
            if (sig[i] != kPngSignature[i]) return false;
        }
        return true;
    }

    /** @brief True if a `python3` interpreter can actually be invoked on this machine. */
    bool pythonAvailable() {
        const int rc = std::system("python3 -c \"import sys\" >/dev/null 2>&1");
        return rc == 0;
    }

    /**
     * @brief A stub renderer that ignores its input entirely and always writes a trivially
     *        valid (1x1) PNG to argv[2], exiting 0. Used to exercise the full C++ success
     *        path without depending on pycirclize/matplotlib.
     */
    constexpr const char *kStubSuccessScript = R"PY(
import sys
import base64
from pathlib import Path

# A minimal valid 1x1 PNG, hardcoded, entirely independent of the input request's content.
_PNG_B64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk"
    "+A8AAQUBAScY42YAAAAASUVORK5CYII="
)

def main() -> int:
    if len(sys.argv) != 3:
        return 2
    out = Path(sys.argv[2])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(base64.b64decode(_PNG_B64))
    return 0

if __name__ == "__main__":
    sys.exit(main())
)PY";

    /**
     * @brief A stub renderer that always fails: prints a distinctive marker to stderr and
     *        exits with a non-zero status. Used to exercise error propagation.
     */
    constexpr const char *kStubFailureScript = R"PY(
import sys
print("STUB_FAILURE_MARKER_9f3c21", file=sys.stderr)
sys.exit(3)
)PY";

    /** @brief Builds a plausible-looking coverage vector of length @p n (simple ramp/wave). */
    std::vector<cdx::Coverage> makeCoverage(const std::size_t n) {
        std::vector<cdx::Coverage> coverage;
        coverage.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            coverage.push_back(static_cast<cdx::Coverage>((i * 5) % 17));
        }
        return coverage;
    }

    /** @brief A coverage vector of length @p n where every value is the NOT_IN_QUERY sentinel. */
    std::vector<cdx::Coverage> makeAllSentinelCoverage(const std::size_t n) {
        return std::vector<cdx::Coverage>(n, cfg::NOT_IN_QUERY);
    }

    output::PlotConfig defaultConfig() {
        output::PlotConfig config;
        config.dpi = 72;
        config.figure_width = 3.0;
        config.figure_height = 3.0;
        return config;
    }
} // namespace

// ===========================================================================
// Category A - pure C++ validation, no Python involved
// ===========================================================================

TEST(CircularPlotQueryValidationTest, ZeroComponentLengthThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";

    EXPECT_THROW(
        output::writeCircularPlotQuery(out, {}, "chr1", /*compo_length=*/0, {0, 0}, defaultConfig()),
        std::invalid_argument
    );
}

TEST(CircularPlotQueryValidationTest, CoverageSizeMismatchThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(5); // compo_length says 10, vector only has 5

    EXPECT_THROW(
        output::writeCircularPlotQuery(out, coverage, "chr1", 10, {0, 4}, defaultConfig()),
        std::invalid_argument
    );
}

TEST(CircularPlotQueryValidationTest, QueryStartBeyondComponentThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(10);

    EXPECT_THROW(
        output::writeCircularPlotQuery(out, coverage, "chr1", 10, {10, 5}, defaultConfig()),
        std::out_of_range
    );
}

TEST(CircularPlotQueryValidationTest, QueryEndBeyondComponentThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(10);

    EXPECT_THROW(
        output::writeCircularPlotQuery(out, coverage, "chr1", 10, {5, 10}, defaultConfig()),
        std::out_of_range
    );
}

TEST(CircularPlotGlobalValidationTest, FewerThanTwoOffsetsThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(10);

    EXPECT_THROW(
        output::writeCircularPlotGlobal(out, coverage, {0}, {"chr1"}, defaultConfig()),
        std::invalid_argument
    );
}

TEST(CircularPlotGlobalValidationTest, FinalOffsetBeyondCoverageSizeThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(10);
    const std::vector<cdx::PosBp> offsets{0, 20};

    EXPECT_THROW(
        output::writeCircularPlotGlobal(out, coverage, offsets, {"chr1"}, defaultConfig()),
        std::invalid_argument
    );
}

TEST(CircularPlotGlobalValidationTest, PerComponentEndBeforeStartThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(10);
    // Component 0: (0,5) valid. Component 1: (5,3), end before start.
    const std::vector<cdx::PosBp> offsets{0, 5, 3};

    EXPECT_THROW(
        output::writeCircularPlotGlobal(out, coverage, offsets, {"chr1", "chr2"}, defaultConfig()),
        std::invalid_argument
    );
}

TEST(CircularPlotGlobalValidationTest, PerComponentEndBeyondCoverageSizeThrows) {
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(60);
    // Component 0: (0,100), end far beyond coverage.size()=60, even though the *final*
    // offset (50) is within bounds - this specifically targets the per-worker check rather
    // than the upfront bp_component_offsets.back() check (FinalOffsetBeyondCoverageSizeThrows).
    const std::vector<cdx::PosBp> offsets{0, 100, 50};

    EXPECT_THROW(
        output::writeCircularPlotGlobal(out, coverage, offsets, {"chr1", "chr2"}, defaultConfig()),
        std::invalid_argument
    );
}

// ===========================================================================
// Category E - script resolution failure (no override, no script next to binary)
// ===========================================================================

TEST(CircularPlotScriptResolutionTest, QueryModeThrowsWhenScriptCannotBeLocated) {
    const EnvVarGuard clear_script("CDX_CIRCULAR_PLOT_SCRIPT");
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(20);

    EXPECT_THROW(
        {
            try {
                output::writeCircularPlotQuery(out, coverage, "chr1", 20, {0, 19}, defaultConfig());
            } catch (const std::runtime_error &e) {
                EXPECT_NE(std::string(e.what()).find("Cannot locate circular_plot.py"), std::string::npos);
                throw;
            }
        },
        std::runtime_error
    );
}

TEST(CircularPlotScriptResolutionTest, GlobalModeThrowsWhenScriptCannotBeLocated) {
    const EnvVarGuard clear_script("CDX_CIRCULAR_PLOT_SCRIPT");
    TempDir dir;
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(30);
    const std::vector<cdx::PosBp> offsets{0, 15, 30};

    EXPECT_THROW(
        output::writeCircularPlotGlobal(out, coverage, offsets, {"chr1", "chr2"}, defaultConfig()),
        std::runtime_error
    );
}

TEST(CircularPlotScriptResolutionTest, NonexistentOverridePathFallsThroughInsteadOfBeingTrusted) {
    TempDir dir;
    const EnvVarGuard override_script(
        "CDX_CIRCULAR_PLOT_SCRIPT", (dir.path() / "does_not_exist.py").string());
    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(20);

    // The override path doesn't exist on disk, so resolution must fall through to the other
    // lookup locations - none of which exist for this test binary either - and still fail
    // with the same "cannot locate" error, rather than blindly trusting the override.
    EXPECT_THROW(
        {
            try {
                output::writeCircularPlotQuery(out, coverage, "chr1", 20, {0, 19}, defaultConfig());
            } catch (const std::runtime_error &e) {
                EXPECT_NE(std::string(e.what()).find("Cannot locate circular_plot.py"), std::string::npos);
                throw;
            }
        },
        std::runtime_error
    );
}

// ===========================================================================
// Category C - stub script that always succeeds: full C++ orchestration path
// ===========================================================================

class CircularPlotStubSuccessTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!pythonAvailable()) {
            GTEST_SKIP() << "No python3 interpreter available on PATH; skipping.";
        }
        script_path_ = dir_.path() / "stub_success.py";
        writeTextFile(script_path_, kStubSuccessScript);
    }

    TempDir dir_;
    fs::path script_path_;
};

TEST_F(CircularPlotStubSuccessTest, FullComponentQueryProducesPngAndCleansUpWorkDir) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeCoverage(50);

    output::writeCircularPlotQuery(out, coverage, "chr1", 50, {0, 49}, defaultConfig());

    EXPECT_TRUE(fs::exists(out));
    EXPECT_TRUE(hasValidPngSignature(out));
    EXPECT_FALSE(fs::exists(dir_.path() / "plot_circular_tmp"));
}

TEST_F(CircularPlotStubSuccessTest, ContiguousSubRangeQueryProducesPng) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeCoverage(200);

    // Sub-range: not full component, not crossing the origin.
    output::writeCircularPlotQuery(out, coverage, "chr2", 200, {50, 120}, defaultConfig());

    EXPECT_TRUE(fs::exists(out));
    EXPECT_TRUE(hasValidPngSignature(out));
}

TEST_F(CircularPlotStubSuccessTest, OriginCrossingQueryProducesPng) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeCoverage(100);

    // query_start > query_end => wraps around the circular origin; exercises the trickiest
    // traversal-building branch of prepareCircularPlotPackage.
    output::writeCircularPlotQuery(out, coverage, "chr3", 100, {90, 10}, defaultConfig());

    EXPECT_TRUE(fs::exists(out));
    EXPECT_TRUE(hasValidPngSignature(out));
}

TEST_F(CircularPlotStubSuccessTest, LogarithmicScaleQueryProducesPng) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    auto config = defaultConfig();
    config.log_base = 2;
    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeCoverage(80);

    output::writeCircularPlotQuery(out, coverage, "chr4", 80, {0, 79}, defaultConfig());
    output::writeCircularPlotQuery(out, coverage, "chr4", 80, {0, 79}, config);

    EXPECT_TRUE(fs::exists(out));
    EXPECT_TRUE(hasValidPngSignature(out));
}

TEST_F(CircularPlotStubSuccessTest, AllSentinelQueryIsInvisibleButStillProducesPng) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeAllSentinelCoverage(40);

    // Every value masked to NaN -> prepareCircularPlotPackage marks the package invisible,
    // but the C++ pipeline still serializes the request and invokes the renderer (visibility
    // is a rendering-time decision made in circular_plot.py, not a reason to skip the
    // subprocess call on the C++ side).
    output::writeCircularPlotQuery(out, coverage, "chr5", 40, {0, 39}, defaultConfig());

    EXPECT_TRUE(fs::exists(out));
    EXPECT_TRUE(hasValidPngSignature(out));
}

TEST_F(CircularPlotStubSuccessTest, GlobalModeWithNameFallbackProducesPng) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeCoverage(90);
    const std::vector<cdx::PosBp> offsets{0, 30, 60, 90};
    // Only one name provided for three components: components 1 and 2 must fall back to
    // std::to_string(component_id) internally instead of crashing.
    const std::vector<std::string> names{"only_one_name"};

    output::writeCircularPlotGlobal(out, coverage, offsets, names, defaultConfig());

    EXPECT_TRUE(fs::exists(out));
    EXPECT_TRUE(hasValidPngSignature(out));
}

TEST_F(CircularPlotStubSuccessTest, GlobalModeWithMixedVisibilityProducesPng) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto normal = makeCoverage(40);
    const auto sentinel = makeAllSentinelCoverage(40);

    std::vector<cdx::Coverage> coverage;
    coverage.insert(coverage.end(), normal.begin(), normal.end());
    coverage.insert(coverage.end(), sentinel.begin(), sentinel.end());

    const std::vector<cdx::PosBp> offsets{0, 40, 80};
    const std::vector<std::string> names{"visible_component", "invisible_component"};

    output::writeCircularPlotGlobal(out, coverage, offsets, names, defaultConfig());

    EXPECT_TRUE(fs::exists(out));
    EXPECT_TRUE(hasValidPngSignature(out));
}

TEST_F(CircularPlotStubSuccessTest, NestedOutputDirectoryIsCreatedAutomatically) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "nested" / "sub" / "plot.png";
    ASSERT_FALSE(fs::exists(out.parent_path()));

    const auto coverage = makeCoverage(30);
    output::writeCircularPlotQuery(out, coverage, "chr1", 30, {0, 29}, defaultConfig());

    EXPECT_TRUE(fs::exists(out));
}

// ===========================================================================
// Category D - stub script that always fails: error propagation
// ===========================================================================

class CircularPlotStubFailureTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!pythonAvailable()) {
            GTEST_SKIP() << "No python3 interpreter available on PATH; skipping.";
        }
        script_path_ = dir_.path() / "stub_failure.py";
        writeTextFile(script_path_, kStubFailureScript);
    }

    TempDir dir_;
    fs::path script_path_;
};

TEST_F(CircularPlotStubFailureTest, NonZeroExitBecomesRuntimeErrorWithCapturedDiagnostics) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeCoverage(20);

    bool threw = false;
    try {
        output::writeCircularPlotQuery(out, coverage, "chr1", 20, {0, 19}, defaultConfig());
    } catch (const std::runtime_error &e) {
        threw = true;
        const std::string message = e.what();
        EXPECT_NE(message.find("Circular plot rendering failed"), std::string::npos);
        EXPECT_NE(message.find("exit status"), std::string::npos);
        // The subprocess's stderr output is captured and folded into the exception message.
        EXPECT_NE(message.find("STUB_FAILURE_MARKER_9f3c21"), std::string::npos);
    }
    EXPECT_TRUE(threw);

    EXPECT_FALSE(fs::exists(out));
}

TEST_F(CircularPlotStubFailureTest, WorkDirSurvivesAFailedRunForDebugging) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto work_dir = dir_.path() / "plot_circular_tmp";
    const auto coverage = makeCoverage(20);

    EXPECT_THROW(
        output::writeCircularPlotQuery(out, coverage, "chr1", 20, {0, 19}, defaultConfig()),
        std::runtime_error
    );

    // cleanupWorkDir() is only reached after a successful render; on failure the request.bin
    // (and the temp work dir containing it) should still be on disk for debugging.
    EXPECT_TRUE(fs::exists(work_dir));
    EXPECT_TRUE(fs::exists(work_dir / "request.bin"));
}

TEST_F(CircularPlotStubFailureTest, GlobalModeAlsoPropagatesRuntimeError) {
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "python3");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path_.string());

    const auto out = dir_.path() / "plot.png";
    const auto coverage = makeCoverage(40);
    const std::vector<cdx::PosBp> offsets{0, 20, 40};

    EXPECT_THROW(
        output::writeCircularPlotGlobal(out, coverage, offsets, {"chr1", "chr2"}, defaultConfig()),
        std::runtime_error
    );
}

// ===========================================================================
// Category F - CDX_PYTHON_EXECUTABLE override plumbing
// ===========================================================================

TEST(CircularPlotPythonExecutableOverrideTest, BogusInterpreterPathCausesFailure) {
    if (!pythonAvailable()) {
        GTEST_SKIP() << "No python3 interpreter available on PATH; skipping.";
    }

    TempDir dir;
    const auto script_path = dir.path() / "stub_success.py";
    writeTextFile(script_path, kStubSuccessScript);

    // A script that would succeed under a real interpreter, but CDX_PYTHON_EXECUTABLE is
    // deliberately pointed at a nonexistent binary - proving the override is genuinely
    // plumbed into the invoked command rather than silently ignored in favor of a real
    // python3 found elsewhere.
    const EnvVarGuard python_override("CDX_PYTHON_EXECUTABLE", "/nonexistent/interpreter_xyz_9f3c");
    const EnvVarGuard script_override("CDX_CIRCULAR_PLOT_SCRIPT", script_path.string());

    const auto out = dir.path() / "plot.png";
    const auto coverage = makeCoverage(20);

    EXPECT_THROW(
        output::writeCircularPlotQuery(out, coverage, "chr1", 20, {0, 19}, defaultConfig()),
        std::runtime_error
    );
    EXPECT_FALSE(fs::exists(out));
}
