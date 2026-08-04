/**
 * @file python_circular_plot.cpp
 * @brief Circular coverage plot rendering: all numeric preparation
 *        (downsampling, smoothing, log transform, tick math, grid layout)
 *        happens here in C++; a Python (pycirclize) subprocess only turns
 *        the already-computed points into a picture.
 *
 * @note Circular rendering used to go through the external Circos tool
 *       (Perl) plus rsvg-convert for SVG->PNG conversion. That pipeline
 *       required a fragile, hard-to-install dependency chain and produced
 *       persistent visual alignment issues (axis/gridline mismatch) that
 *       many rounds of Circos config patches never fully fixed. It was
 *       replaced by a Python (pycirclize) subprocess, since pycirclize
 *       already renders exactly the look the Python prototype
 *       (output_tools.py: `_plotCircularGraph`) established.
 *
 *       Initially, that subprocess also duplicated the prototype's numeric
 *       preparation (smoothing, downsampling, log transform, tick math) in
 *       Python/NumPy. That work is CPU-bound and scales with full-resolution
 *       coverage length (potentially millions of points before
 *       downsampling), so it has been moved here, reusing the exact same
 *       numeric core the linear (Cairo) backend already relies on
 *       (chooseGlobalGraphGrid, prepareCoverageForPlot, logTransformCoverage,
 *       logCoverageTicks, calculateCoverageTicks - see linear_plot.cpp).
 *       The Python side (python_script/circular_plot.py) now only receives
 *       already-downsampled/smoothed/log-transformed points plus tick
 *       values, and is responsible solely for the pycirclize/matplotlib
 *       rendering itself (coordinate re-projection onto the polar track,
 *       segment/gap drawing, axis labels, title, PNG export).
 *
 *       Python is a RUNTIME-only dependency of this one feature (circular
 *       graphs). The project still compiles with no Python present; only
 *       running with `--component-type circular` requires python3 plus
 *       numpy/matplotlib/pycirclize to be installed (see
 *       CIRCULAR_PLOT_SETUP.md). Linear graphs remain pure C++/Cairo,
 *       in-process, no subprocess (src/cairo_plot.cpp).
 */

#include "output_plot.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

namespace output {
    namespace {
        // -------------------------------------------------------------
        // Numeric preparation (mirrors what the Python prototype's
        // _plotCircularGraph used to compute in NumPy, steps 1-5: traversal
        // extraction, sentinel masking, invariant-window smoothing +
        // downsampling, log transform, tick math).
        // -------------------------------------------------------------

        /**
         * @brief Everything a circular plot needs to draw ONE component (or
         *        one sub-query of it), already downsampled/smoothed/
         *        log-transformed: pycirclize/matplotlib rendering is all
         *        that remains to be done, in Python.
         *
         * `local_x` mirrors the Python prototype's `plot_x`: indices into
         * the *traversal* (the extracted, possibly origin-crossing-
         * concatenated sub-range), NOT absolute genomic coordinates. The
         * Python side still does the small (post-downsampling, so cheap)
         * arithmetic that turns these into display coordinates on the
         * polar track, exactly as _plotCircularGraph did - only the
         * expensive full-resolution numeric work moved here.
         */
        struct CircularPlotPackage {
            std::string component_name;
            cdx::PosBp compo_length = 0;
            cdx::PosBp query_start = 0;
            cdx::PosBp query_end = 0;
            bool full_component = true;
            bool crosses_origin = false;
            bool visible = true; // false => nothing to draw (empty/fully-masked range)
            bool logarithmic = false;
            int log_base = 0;
            double y_upper_limit = 1.0;

            std::vector<std::size_t> local_x;
            std::vector<double> plot_y; // may contain NaN (masked/out-of-query gaps)

            std::vector<double> tick_raw_values;
            std::vector<double> tick_display_values;
        };

        /**
         * @brief Prepares a circular plot package: extraction of the
         *        traversal (handling origin-crossing queries), sentinel
         *        masking to NaN (preserved as visual gaps, unlike the
         *        linear backend which zero-fills - cf. prepareLinearPlotPackage),
         *        invariant-window smoothing/downsampling, and radial scale
         *        (log or linear) tick computation.
         *
         * `component_coverage` must span the WHOLE component (never
         * pre-sliced): circular wrap-around smoothing and origin-crossing
         * queries both need the full context, exactly as documented on
         * writeCircularPlotQuery/Global.
         */
        CircularPlotPackage prepareCircularPlotPackage(
            const std::vector<cdx::Coverage> &component_coverage,
            const std::string &component_name,
            const cdx::PosBp compo_length,
            const std::pair<cdx::PosBp, cdx::PosBp> &query_bound,
            const PlotConfig &config
        ) {
            if (compo_length == 0) {
                throw std::invalid_argument("compo_length must be greater than zero.");
            }
            if (component_coverage.size() != static_cast<std::size_t>(compo_length)) {
                throw std::invalid_argument(
                    "component_coverage size must equal compo_length "
                    "(the array must NOT be pre-sliced for circular plots)."
                );
            }

            CircularPlotPackage pkg;
            pkg.component_name = component_name;
            pkg.compo_length = compo_length;
            pkg.query_start = query_bound.first;
            pkg.query_end = query_bound.second;

            if (pkg.query_start >= compo_length || pkg.query_end >= compo_length) {
                throw std::out_of_range("query_bound is outside component bounds.");
            }

            pkg.full_component = (pkg.query_start == 0 && pkg.query_end == compo_length - 1);
            pkg.crosses_origin = !pkg.full_component && (pkg.query_start > pkg.query_end);

            // --- 1. Extraction du parcours (traversal) -----------------------
            std::vector<double> traversal;

            if (pkg.full_component) {
                traversal.resize(component_coverage.size());
                for (std::size_t i = 0; i < component_coverage.size(); ++i) {
                    traversal[i] = static_cast<double>(component_coverage[i]);
                }
            } else if (pkg.crosses_origin) {
                const auto first_len = static_cast<std::size_t>(compo_length - pkg.query_start);
                const auto second_len = static_cast<std::size_t>(pkg.query_end + 1);
                traversal.reserve(first_len + second_len);

                for (cdx::PosBp bp = pkg.query_start; bp < compo_length; ++bp) {
                    traversal.push_back(static_cast<double>(component_coverage[static_cast<std::size_t>(bp)]));
                }
                for (cdx::PosBp bp = 0; bp <= pkg.query_end; ++bp) {
                    traversal.push_back(static_cast<double>(component_coverage[static_cast<std::size_t>(bp)]));
                }
            } else {
                const auto start = static_cast<std::size_t>(pkg.query_start);
                const auto end = static_cast<std::size_t>(pkg.query_end);

                traversal.resize(end - start + 1);
                for (std::size_t i = 0; i < traversal.size(); ++i) {
                    traversal[i] = static_cast<double>(component_coverage[start + i]);
                }
            }

            // --- 2. Sentinelles -> NaN (trous visuels, fidèle au prototype) --
            constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
            bool any_finite = false;
            for (double &value: traversal) {
                if (value >= static_cast<double>(cfg::NOT_IN_QUERY)) {
                    value = kNaN;
                } else {
                    any_finite = true;
                }
            }

            if (!any_finite) {
                pkg.visible = false;
                return pkg;
            }

            // --- 3. Lissage à fenêtre invariante en bp, topologie adaptée ----
            const double adjusted_smoothing =
                (config.smoothing <= 0.0)
                    ? config.smoothing
                    : config.smoothing * static_cast<double>(compo_length) / static_cast<double>(traversal.size());

            const PlotData plot_data = prepareCoverageForPlot(
                traversal,
                adjusted_smoothing,
                config.max_plot_points,
                pkg.full_component ? Topology::Circular : Topology::Linear
            );

            if (plot_data.y.empty()) {
                pkg.visible = false;
                return pkg;
            }

            double maximum_raw = 0.0;
            bool any_finite_plot_y = false;
            for (const double value: plot_data.y) {
                if (std::isfinite(value)) {
                    maximum_raw = std::max(maximum_raw, value);
                    any_finite_plot_y = true;
                }
            }

            if (!any_finite_plot_y) {
                pkg.visible = false;
                return pkg;
            }

            pkg.local_x = plot_data.x;

            // --- 4. Échelle radiale (log vs linéaire) ------------------------
            if (config.log_base.has_value()) {
                const int log_base = *config.log_base;

                pkg.logarithmic = true;
                pkg.log_base = log_base;
                pkg.plot_y = logTransformCoverage(plot_data.y, log_base);

                const LogCoverageTicks ticks = logCoverageTicks(maximum_raw, log_base);
                pkg.tick_raw_values = ticks.raw_values;
                pkg.tick_display_values = ticks.display_values;
                pkg.y_upper_limit = std::max(ticks.display_upper_limit, 1.0);
            } else {
                pkg.logarithmic = false;
                pkg.log_base = 0;
                pkg.plot_y = plot_data.y;

                const CoverageTicks ticks = calculateCoverageTicks(maximum_raw);
                pkg.tick_raw_values = ticks.values;
                pkg.tick_display_values = ticks.values;
                pkg.y_upper_limit = std::max(ticks.upper_limit, 1.0);
            }

            return pkg;
        }

        // -------------------------------------------------------------
        // Binary request writer (little-endian), matching the format
        // documented at the top of python_script/circular_plot.py.
        // -------------------------------------------------------------
        class BinaryRequestWriter {
        public:
            explicit BinaryRequestWriter(const std::filesystem::path &path)
                : out_(path, std::ios::binary) {
                if (!out_) {
                    throw std::runtime_error("Cannot open circular plot request file: " + path.string());
                }
            }

            void bytes(const void *data, const std::size_t size) {
                out_.write(static_cast<const char *>(data), static_cast<std::streamsize>(size));
            }

            void u8(const std::uint8_t value) { bytes(&value, sizeof(value)); }
            void i32(const std::int32_t value) { bytes(&value, sizeof(value)); }
            void u32(const std::uint32_t value) { bytes(&value, sizeof(value)); }
            void u64(const std::uint64_t value) { bytes(&value, sizeof(value)); }
            void f64(const double value) { bytes(&value, sizeof(value)); }

            void string(const std::string &value) {
                const auto length = static_cast<std::uint16_t>(std::min<std::size_t>(value.size(), 0xFFFF));
                u16(length);
                bytes(value.data(), length);
            }

            void u64_array(const std::vector<std::size_t> &values) {
                if (values.empty()) return;

                // Bulk-write when std::size_t already IS a 64-bit value (true
                // on every platform this project targets: macOS/Linux, 64-bit
                // only, per the libvgio/Homebrew/apt setup docs) - avoids one
                // iostream call per point, which mattered once local_x arrays
                // routinely held tens of thousands of points across many
                // global-mode components. Falls back to a per-element loop
                // on any platform where that assumption doesn't hold.
                if constexpr (sizeof(std::size_t) == sizeof(std::uint64_t)) {
                    bytes(values.data(), values.size() * sizeof(std::uint64_t));
                } else {
                    for (const std::size_t value: values) {
                        u64(static_cast<std::uint64_t>(value));
                    }
                }
            }

            void f64_array(const std::vector<double> &values) {
                if (!values.empty()) {
                    bytes(values.data(), values.size() * sizeof(double));
                }
            }

            void package(const CircularPlotPackage &pkg) {
                string(pkg.component_name);
                u64(static_cast<std::uint64_t>(pkg.compo_length));
                u64(static_cast<std::uint64_t>(pkg.query_start));
                u64(static_cast<std::uint64_t>(pkg.query_end));
                u8(pkg.full_component ? 1 : 0);
                u8(pkg.crosses_origin ? 1 : 0);
                u8(pkg.visible ? 1 : 0);
                u8(pkg.logarithmic ? 1 : 0);
                i32(pkg.log_base);
                f64(pkg.y_upper_limit);

                u64(static_cast<std::uint64_t>(pkg.local_x.size()));
                u64_array(pkg.local_x);
                f64_array(pkg.plot_y);

                u32(static_cast<std::uint32_t>(pkg.tick_raw_values.size()));
                f64_array(pkg.tick_raw_values);
                f64_array(pkg.tick_display_values);
            }

            void finish() {
                out_.flush();
                if (!out_.good()) {
                    throw std::runtime_error("Error writing circular plot request file.");
                }
            }

        private:
            void u16(const std::uint16_t value) { bytes(&value, sizeof(value)); }

            std::ofstream out_;
        };

        void writeCommonHeader(
            BinaryRequestWriter &writer,
            const std::uint8_t mode,
            const PlotConfig &config
        ) {
            writer.bytes("CXCP", 4);
            writer.u8(2); // version: 2 = pre-processed packages (numeric prep done in C++)
            writer.u8(mode);
            writer.u32(static_cast<std::uint32_t>(config.dpi));
            writer.f64(config.figure_width);
            writer.f64(config.figure_height);
            writer.string(config.line_color);
            writer.string(config.fill_color);
        }

        /**
         * @brief Locates the directory containing the running executable, so
         *        circular_plot.py (copied/installed next to it by CMake) can
         *        be found without relying on the current working directory.
         */
        std::filesystem::path executableDirectory() {
#if defined(__APPLE__)
            std::array<char, 4096> buffer{};
            std::uint32_t size = static_cast<std::uint32_t>(buffer.size());
            if (_NSGetExecutablePath(buffer.data(), &size) == 0) {
                std::error_code error_code;
                auto resolved = std::filesystem::canonical(std::filesystem::path(buffer.data()), error_code);
                if (!error_code) {
                    return resolved.parent_path();
                }
            }
#elif defined(__linux__)
            std::array<char, 4096> buffer{};
            const ssize_t length = readlink("/proc/self/exe", buffer.data(), buffer.size() - 1);
            if (length != -1) {
                buffer[static_cast<std::size_t>(length)] = '\0';
                return std::filesystem::path(buffer.data()).parent_path();
            }
#endif
            return std::filesystem::current_path();
        }

        /**
         * @brief Resolves the path to circular_plot.py, trying (in order): an
         *        explicit environment override, next to the running
         *        executable (build output dir / installed bin dir, where
         *        CMake copies/installs it), then the source tree location
         *        (dev convenience when running straight from a build dir).
         */
        std::filesystem::path resolveCircularPlotScript() {
            if (const char *override_path = std::getenv("CDX_CIRCULAR_PLOT_SCRIPT")) {
                if (std::filesystem::exists(override_path)) {
                    return override_path;
                }
            }

            const std::filesystem::path next_to_binary = executableDirectory() / "circular_plot.py";
            if (std::filesystem::exists(next_to_binary)) {
                return next_to_binary;
            }

#ifdef CDX_CIRCULAR_PLOT_SCRIPT_SOURCE_DIR
            const std::filesystem::path source_tree_copy =
                    std::filesystem::path(CDX_CIRCULAR_PLOT_SCRIPT_SOURCE_DIR) / "circular_plot.py";
            if (std::filesystem::exists(source_tree_copy)) {
                return source_tree_copy;
            }
#endif

            throw std::runtime_error(
                "Cannot locate circular_plot.py. Set the CDX_CIRCULAR_PLOT_SCRIPT "
                "environment variable to its full path, or ensure it was installed "
                "next to the cdx_coverage executable (see CIRCULAR_PLOT_SETUP.md)."
            );
        }

        /**
         * @brief Resolves which python3 binary to invoke, in order:
         *
         *        1. CDX_PYTHON_EXECUTABLE, if set - explicit override,
         *           useful e.g. when an IDE run configuration doesn't
         *           inherit the shell PATH that would normally expose a
         *           manually-activated venv (see CIRCULAR_PLOT_SETUP.md).
         *        2. The bundled venv CMake provisions automatically next to
         *           the executable at build time (cmake/setup_circular_env.sh,
         *           "pyenv/bin/python3") - so a normal build is enough to get
         *           working circular graphs, with no manual venv setup by
         *           the end user. Skipped if that provisioning step failed
         *           or was never run (e.g. no Python3 found at configure
         *           time, or no network access during the build).
         *        3. Plain "python3" resolved from the current process's
         *           PATH, as a last resort.
         */
        std::string resolvePythonExecutable() {
            if (const char *override_path = std::getenv("CDX_PYTHON_EXECUTABLE")) {
                if (*override_path != '\0') {
                    return override_path;
                }
            }

            const std::filesystem::path bundled_venv_python =
                    executableDirectory() / "pyenv" / "bin" / "python3";
            if (std::filesystem::exists(bundled_venv_python)) {
                return bundled_venv_python.string();
            }

            return "python3";
        }

        /**
         * @brief Invokes `python3 circular_plot.py <request.bin> <output.png>`
         *        and captures stdout/stderr for diagnostics on failure.
         */
        void runCircularPlotScriptOrThrow(
            const std::filesystem::path &script_path,
            const std::filesystem::path &request_path,
            const std::filesystem::path &output_png
        ) {
            const std::string python_executable = resolvePythonExecutable();

            std::ostringstream command;
            command << std::quoted(python_executable) << ' '
                    << std::quoted(script_path.string())
                    << ' ' << std::quoted(request_path.string())
                    << ' ' << std::quoted(output_png.string())
                    << " 2>&1";

            FILE *pipe = popen(command.str().c_str(), "r");
            if (pipe == nullptr) {
                throw std::runtime_error(
                    "Failed to launch '" + python_executable + "' to render the circular graph. "
                    "Is Python 3 installed and available in PATH? If running from an IDE that "
                    "doesn't inherit your shell's venv, set CDX_PYTHON_EXECUTABLE to the exact "
                    "interpreter path (see CIRCULAR_PLOT_SETUP.md)."
                );
            }

            std::string captured_output;
            std::array<char, 512> buffer{};

            while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
                captured_output += buffer.data();
            }

            const int status = pclose(pipe);

            if (status != 0) {
                throw std::runtime_error(
                    "Circular plot rendering failed (exit status " + std::to_string(status) + ").\n"
                    "Script: " + script_path.string() + "\n"
                    "Request file: " + request_path.string() + "\n"
                    "--- python output ---\n" + captured_output
                );
            }

            if (!std::filesystem::exists(output_png)) {
                throw std::runtime_error(
                    "circular_plot.py completed successfully but the expected output "
                    "file was not found: " + output_png.string()
                );
            }
        }

        /**
         * @brief Best-effort removal of the per-graph work directory
         *        (request.bin) once the PNG has been produced successfully.
         *        Never throws: a leftover temp directory is a cosmetic issue,
         *        not worth failing an otherwise-successful render over.
         */
        void cleanupWorkDir(const std::filesystem::path &work_dir) {
            std::error_code error_code;
            std::filesystem::remove_all(work_dir, error_code);
        }
    } // namespace

    void writeCircularPlotQuery(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &component_coverage,
        const std::string &component_name,
        const cdx::PosBp compo_length,
        const std::pair<cdx::PosBp, cdx::PosBp> query_bound,
        const PlotConfig &config
    ) {
        const std::filesystem::path image_dir =
                output_png.parent_path().empty() ? std::filesystem::path(".") : output_png.parent_path();
        std::filesystem::create_directories(image_dir);

        const std::filesystem::path work_dir = image_dir / (output_png.stem().string() + "_circular_tmp");
        std::filesystem::create_directories(work_dir);
        const std::filesystem::path request_path = work_dir / "request.bin";

        const CircularPlotPackage pkg = prepareCircularPlotPackage(
            component_coverage, component_name, compo_length, query_bound, config
        );

        {
            BinaryRequestWriter writer(request_path);
            writeCommonHeader(writer, /*mode=*/0, config);
            writer.package(pkg);
            writer.finish();
        }

        const std::filesystem::path script_path = resolveCircularPlotScript();
        runCircularPlotScriptOrThrow(script_path, request_path, output_png);
        cleanupWorkDir(work_dir);
    }

    void writeCircularPlotGlobal(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &flat_bp_coverage,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names,
        const PlotConfig &config
    ) {
        if (bp_component_offsets.size() < 2) {
            throw std::invalid_argument("bp_component_offsets must contain at least two offsets.");
        }
        if (bp_component_offsets.back() > static_cast<cdx::PosBp>(flat_bp_coverage.size())) {
            throw std::invalid_argument("The final component offset exceeds the coverage size.");
        }

        const std::filesystem::path image_dir =
                output_png.parent_path().empty() ? std::filesystem::path(".") : output_png.parent_path();
        std::filesystem::create_directories(image_dir);

        const std::filesystem::path work_dir = image_dir / (output_png.stem().string() + "_circular_tmp");
        std::filesystem::create_directories(work_dir);
        const std::filesystem::path request_path = work_dir / "request.bin";

        const std::size_t component_count = bp_component_offsets.size() - 1;
        const GridLayout grid = chooseGlobalGraphGrid(component_count);

        // Numeric prep per component is independent CPU-bound work (same
        // pattern as the old per-component Circos rendering, and the linear
        // backend's multi-panel prep): parallelize it, then serialize the
        // results sequentially (std::ofstream is not thread-safe) and hand
        // the whole grid to ONE Python invocation - not one per component -
        // since process-launch overhead was the actual bottleneck reported.
        std::vector<CircularPlotPackage> packages(component_count);
        std::exception_ptr worker_exception;
        std::mutex exception_mutex;

        auto prepare_component = [&](const std::size_t component_id) {
            try {
                const auto start = static_cast<std::size_t>(bp_component_offsets[component_id]);
                const auto end = static_cast<std::size_t>(bp_component_offsets[component_id + 1]);

                if (end < start || end > flat_bp_coverage.size()) {
                    throw std::invalid_argument("bp_component_offsets is inconsistent with flat_bp_coverage size.");
                }

                const std::vector<cdx::Coverage> component_coverage(
                    flat_bp_coverage.begin() + static_cast<std::ptrdiff_t>(start),
                    flat_bp_coverage.begin() + static_cast<std::ptrdiff_t>(end)
                );

                const auto compo_length = static_cast<cdx::PosBp>(component_coverage.size());
                const std::string name = (component_id < component_names.size())
                                              ? component_names[component_id]
                                              : std::to_string(component_id);

                packages[component_id] = prepareCircularPlotPackage(
                    component_coverage, name, compo_length, {0, compo_length - 1}, config
                );
            } catch (...) {
                const std::lock_guard<std::mutex> lock(exception_mutex);
                if (!worker_exception) {
                    worker_exception = std::current_exception();
                }
            }
        };

        {
            const unsigned int hardware_threads = std::thread::hardware_concurrency();
            const std::size_t thread_count = std::min(
                component_count,
                static_cast<std::size_t>(std::max(1u, hardware_threads))
            );

            std::atomic<std::size_t> next_index{0};
            std::vector<std::thread> workers;
            workers.reserve(thread_count);

            for (std::size_t t = 0; t < thread_count; ++t) {
                workers.emplace_back([&]() {
                    for (std::size_t idx; (idx = next_index.fetch_add(1)) < component_count;) {
                        prepare_component(idx);
                    }
                });
            }

            for (auto &worker: workers) {
                worker.join();
            }
        }

        if (worker_exception) {
            std::rethrow_exception(worker_exception);
        }

        {
            BinaryRequestWriter writer(request_path);
            writeCommonHeader(writer, /*mode=*/1, config);

            writer.f64(config.figure_width);
            writer.f64(config.figure_height);
            writer.u32(static_cast<std::uint32_t>(grid.rows));
            writer.u32(static_cast<std::uint32_t>(grid.columns));

            writer.u32(static_cast<std::uint32_t>(component_count));
            for (const CircularPlotPackage &pkg: packages) {
                writer.package(pkg);
            }

            writer.finish();
        }

        const std::filesystem::path script_path = resolveCircularPlotScript();
        runCircularPlotScriptOrThrow(script_path, request_path, output_png);
        cleanupWorkDir(work_dir);
    }
} // namespace output
