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
        * @brief Preprocessed data required to render a single circular coverage track.
        *
        * Contains all information needed by the Python rendering layer to draw one
        * component or query sub-range in a circular plot. Coverage values have
        * already been extracted, filtered, optionally transformed, and reduced to
        * plotting resolution; the Python side is only responsible for converting
        * the provided coordinates into display-space coordinates and rendering the
        * figure.
        *
        * Coordinates in local_x are traversal-relative positions within the plotted
        * region rather than absolute component coordinates. This allows the package
        * to represent both contiguous ranges and origin-crossing ranges using a
        * single linear traversal.
        */
        struct CircularPlotPackage {
            // Component and query-region metadata.
            std::string component_name;
            cdx::PosBp compo_length = 0;
            cdx::PosBp query_start = 0;
            cdx::PosBp query_end = 0;

            // Traversal characteristics.
            bool full_component = true; // Entire component is plotted.
            bool crosses_origin = false; // Query interval wraps around the origin.
            bool visible = true; // False == no drawable data remains.

            // Y-axis transformation and scaling parameters.
            bool logarithmic = false;
            int log_base = 0;
            double y_upper_limit = 1.0;

            // Downsampled plotting coordinates. local_x contains positions in the
            // extracted traversal, while plot_y contains the corresponding values.
            std::vector<std::size_t> local_x;
            std::vector<double> plot_y; // may contain NaN (masked/out-of-query gaps)

            // Tick labels for the radial axis. Raw values are used for placement,
            // while display values contain the transformed values shown to users.
            std::vector<double> tick_raw_values;
            std::vector<double> tick_display_values;
        };

        /**
         * @brief Build a render-ready circular coverage plot package.
         *
         * Extracts the requested query interval from a component coverage vector,
         * handling full-component, linear sub-range, and origin-crossing circular
         * traversals. Coverage values are then converted into plotting data through
         * sentinel masking, smoothing, downsampling, and optional logarithmic
         * transformation.
         *
         * Unlike the linear plotting pipeline, positions outside the query are
         * preserved as NaN values so that the renderer displays visual gaps rather
         * than zero-filled regions.
         *
         * The returned package contains all data required by the Python rendering
         * layer, including plotting coordinates, axis scaling information, and
         * radial tick metadata.
         *
         * @param component_coverage Coverage values for the entire component.
         * @param component_name Name of the component being plotted.
         * @param compo_length Total component length in base pairs.
         * @param query_bound Inclusive query interval within the component.
         * @param config Plot configuration controlling smoothing, downsampling, and axis scaling.
         *
         * @return A fully prepared CircularPlotPackage.
         *
         * @throws std::invalid_argument If component dimensions are inconsistent.
         * @throws std::out_of_range If query_bound falls outside component bounds.
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

            // Determine whether the requested interval spans
            // the whole component, a contiguous sub-range, or wraps around the circular origin.
            pkg.full_component = (pkg.query_start == 0 && pkg.query_end == compo_length - 1);
            pkg.crosses_origin = !pkg.full_component && (pkg.query_start > pkg.query_end);


            // Step 1: Construct a linear traversal corresponding to the query region.
            // Origin-crossing intervals are represented by concatenating the tail and
            // head of the circular component into a single traversal.
            std::vector<double> traversal;

            // Full-component view: preserve the entire coverage vector.
            if (pkg.full_component) {
                traversal.resize(component_coverage.size());
                for (std::size_t i = 0; i < component_coverage.size(); ++i) {
                    traversal[i] = static_cast<double>(component_coverage[i]);
                }

                // Origin-crossing interval: traverse [query_start, end] followed by
                // [0, query_end] so that the wrapped query becomes a single continuous
                // sequence for smoothing, downsampling, and plotting.
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

                // Standard contiguous interval that does not cross the component origin.
            } else {
                const auto start = static_cast<std::size_t>(pkg.query_start);
                const auto end = static_cast<std::size_t>(pkg.query_end);

                traversal.resize(end - start + 1);
                for (std::size_t i = 0; i < traversal.size(); ++i) {
                    traversal[i] = static_cast<double>(component_coverage[start + i]);
                }
            }

            // Step 2: Convert out-of-query sentinel values to NaN so the renderer
            // produces visible gaps instead of artificial zero-coverage segments.
            constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
            bool any_finite = false;
            for (double &value: traversal) {
                if (value >= static_cast<double>(cfg::NOT_IN_QUERY)) {
                    value = kNaN;
                } else {
                    any_finite = true;
                }
            }

            // No finite coverage values remain after masking.
            if (!any_finite) {
                pkg.visible = false;
                return pkg;
            }

            // Step 3 : Preserve an approximately constant smoothing window in base-pair units.
            // When plotting only a subrange, the smoothing parameter is rescaled from component
            // coordinates to traversal coordinates.
            const double adjusted_smoothing =
                    (config.smoothing <= 0.0)
                        ? config.smoothing
                        : config.smoothing * static_cast<double>(compo_length) / static_cast<double>(traversal.size());

            // Smooth and downsample the traversal. Full-component plots use circular
            // topology so values near the origin interact correctly across the wrap;
            // subrange plots use linear topology.
            const PlotData plot_data = prepareCoverageForPlot(
                traversal,
                adjusted_smoothing,
                config.max_plot_points,
                pkg.full_component ? Topology::Circular : Topology::Linear
            );

            // Downsampling and masking can eliminate all drawable points.
            if (plot_data.y.empty()) {
                pkg.visible = false;
                return pkg;
            }

            // Find the largest finite plotted value for axis scaling and tick
            // generation. NaN gap values are ignored.
            double maximum_raw = 0.0;
            bool any_finite_plot_y = false;
            for (const double value: plot_data.y) {
                if (std::isfinite(value)) {
                    maximum_raw = std::max(maximum_raw, value);
                    any_finite_plot_y = true;
                }
            }

            // The plotted series contains only NaN values.
            if (!any_finite_plot_y) {
                pkg.visible = false;
                return pkg;
            }

            pkg.local_x = plot_data.x;

            // Step 4: Configure radial-axis scaling and prepare the values that
            // will be displayed on the circular coverage track.

            // Transform coverage values into log space and generate both the raw
            // coverage tick positions and their displayed logarithmic labels.
            if (config.log_base.has_value()) {
                const int log_base = *config.log_base;

                pkg.logarithmic = true;
                pkg.log_base = log_base;
                pkg.plot_y = logTransformCoverage(plot_data.y, log_base);

                const LogCoverageTicks ticks = logCoverageTicks(maximum_raw, log_base);
                pkg.tick_raw_values = ticks.raw_values;
                pkg.tick_display_values = ticks.display_values;
                pkg.y_upper_limit = std::max(ticks.display_upper_limit, 1.0);

                // Keep coverage values unchanged and generate standard linear-axis ticks.
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

        /**
        * @brief Writes preprocessed circular-plot requests in the binary format
        * consumed by circular_plot.py.
        *
        * This class serializes CircularPlotPackage objects and plot configuration
        * metadata into a compact little-endian binary representation. The generated
        * file contains all numerical preprocessing results, allowing the Python
        * rendering layer to focus solely on visualization.
        */
        class BinaryRequestWriter {
        public:
            // Open the binary request file that will be consumed by the Python
            // rendering backend.
            explicit BinaryRequestWriter(const std::filesystem::path &path)
                : out_(path, std::ios::binary) {
                if (!out_) {
                    throw std::runtime_error("Cannot open circular plot request file: " + path.string());
                }
            }

            // Low-level primitive used by all typed writers.
            void bytes(const void *data, const std::size_t size) {
                out_.write(static_cast<const char *>(data), static_cast<std::streamsize>(size));
            }

            void u8(const std::uint8_t value) { bytes(&value, sizeof(value)); }
            void i32(const std::int32_t value) { bytes(&value, sizeof(value)); }
            void u32(const std::uint32_t value) { bytes(&value, sizeof(value)); }
            void u64(const std::uint64_t value) { bytes(&value, sizeof(value)); }
            void f64(const double value) { bytes(&value, sizeof(value)); }

            // Strings are stored as: [u16 length][raw bytes].
            // Length is capped at UINT16_MAX to match the on-disk format specification.
            void string(const std::string &value) {
                const auto length = static_cast<std::uint16_t>(std::min<std::size_t>(value.size(), 0xFFFF));
                u16(length);
                bytes(value.data(), length);
            }

            // Fast path: bulk-write the entire buffer when std::size_t is already
            // 64-bit and layout-compatible with the serialized format. Otherwise,
            // fall back to element-by-element conversion.
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

            // Arrays of doubles are already stored in the required binary layout.
            void f64_array(const std::vector<double> &values) {
                if (!values.empty()) {
                    bytes(values.data(), values.size() * sizeof(double));
                }
            }

            /**
            * @brief Serialize a single CircularPlotPackage.
            *
            * Fields are written in the exact order expected by the Python reader.
            * Any modification to this layout must be mirrored in circular_plot.py.
            */
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

                // Tick positions and display labels are written separately because they
                // differ when logarithmic scaling is enabled.
                u32(static_cast<std::uint32_t>(pkg.tick_raw_values.size()));
                f64_array(pkg.tick_raw_values);
                f64_array(pkg.tick_display_values);
            }

            // Ensure all buffered data reaches disk and detect write failures.
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


        /**
        * @brief Write the binary file header shared by all circular-plot requests.
        *
        * The header identifies the file format, version, rendering mode, and
        * global plotting configuration. It must be read before any serialized
        * CircularPlotPackage objects.
        *
        * @param writer Binary request writer receiving the encoded header.
        * @param mode Plot-generation mode understood by the Python backend.
        * @param config Rendering configuration shared by all plot packages.
        */
        void writeCommonHeader(
            BinaryRequestWriter &writer,
            const std::uint8_t mode,
            const PlotConfig &config
        ) {
            // File magic used to validate that the payload is a circular-plot request.
            writer.bytes("CXCP", 4);

            // version: 2 = pre-processed packages (numeric prep done in C++)
            writer.u8(2);
            writer.u8(mode);

            // Global rendering parameters shared by every package in the request.
            writer.u32(static_cast<std::uint32_t>(config.dpi));
            writer.f64(config.figure_width);
            writer.f64(config.figure_height);
            writer.string(config.line_color);
            writer.string(config.fill_color);
        }

        /**
         * @brief Determine the directory containing the currently running executable.
         *
         * Resolves the executable location using platform-specific mechanisms and
         * returns its parent directory. This allows runtime resources installed
         * alongside the executable, such as circular_plot.py, to be located
         * independently of the current working directory.
         *
         * On failure, falls back to std::filesystem::current_path().
         *
         * @return Filesystem path to the executable directory, or the current
         *         working directory if the executable location cannot be resolved.
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
         * @brief Locate the Python circular-plot rendering script.
         *
         * Searches for circular_plot.py using several locations in order of
         * preference:
         *   1. The path specified by the CDX_CIRCULAR_PLOT_SCRIPT environment variable.
         *   2. The directory containing the running executable.
         *   3. The source-tree copy configured at build time (development fallback).
         *
         * This allows both installed deployments and local development builds to
         * locate the rendering script without relying on the current working
         * directory.
         *
         * @return Absolute path to circular_plot.py.
         *
         * @throws std::runtime_error If the script cannot be found in any of the
         *         supported locations.
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
                "next to the cdx_coverage executable (see python_script/CIRCULAR_PLOT_SETUP.md)."
            );
        }

        /**
         * @brief Determine which Python interpreter should be used for circular-plot rendering.
         *
         * The interpreter is resolved in the following order:
         *   1. The path specified by the CDX_PYTHON_EXECUTABLE environment variable.
         *   2. A bundled virtual environment installed alongside the executable.
         *   3. A python3 executable available on the system PATH.
         *
         * This resolution strategy allows explicit user overrides, supports
         * self-contained deployments, and falls back to standard system Python
         * installations when no dedicated environment is available.
         *
         * @return Path or command name of the Python interpreter to invoke.
         */
        std::string resolvePythonExecutable() {
            // try to find the executable in the environment variable
            if (const char *override_path = std::getenv("CDX_PYTHON_EXECUTABLE")) {
                if (*override_path != '\0') {
                    return override_path;
                }
            }

            // use the bundled virtual environment installed alongside the executable
            // (usually done automatically at compilation)
            const std::filesystem::path bundled_venv_python =
                    executableDirectory() / "pyenv" / "bin" / "python3";
            if (std::filesystem::exists(bundled_venv_python)) {
                return bundled_venv_python.string();
            }

            return "python3";
        }

        /**
         * @brief Execute the Python circular-plot renderer and verify output creation.
         *
         * Launches circular_plot.py with the specified request file and output path,
         * captures the renderer's stdout/stderr streams, and reports any execution
         * failures with diagnostic output.
         *
         * The function additionally verifies that the expected PNG file was created
         * before returning successfully.
         *
         * @param script_path Path to circular_plot.py.
         * @param request_path Path to the serialized plot request.
         * @param output_png Expected output image path.
         *
         * @throws std::runtime_error If the Python process cannot be started,
         * terminates with a non-zero exit status, or fails to produce
         * the expected output image.
         */
        void runCircularPlotScriptOrThrow(
            const std::filesystem::path &script_path,
            const std::filesystem::path &request_path,
            const std::filesystem::path &output_png
        ) {
            const std::string python_executable = resolvePythonExecutable();


            // Build a shell command with proper path quoting and redirect stderr
            // into stdout so all diagnostic output can be captured through one pipe.
            std::ostringstream command;
            command << std::quoted(python_executable) << ' '
                    << std::quoted(script_path.string())
                    << ' ' << std::quoted(request_path.string())
                    << ' ' << std::quoted(output_png.string())
                    << " 2>&1";

            // Launch the renderer and open a pipe for consuming its output.
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

            // Collect all renderer output so it can be included in exception messages if the script fails.
            while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
                captured_output += buffer.data();
            }

            const int status = pclose(pipe);

            // Return Python-side errors together with any emitted diagnostics.
            if (status != 0) {
                throw std::runtime_error(
                    "Circular plot rendering failed (exit status " + std::to_string(status) + ").\n"
                    "Script: " + script_path.string() + "\n"
                    "Request file: " + request_path.string() + "\n"
                    "--- python output ---\n" + captured_output
                );
            }

            // Successful process completion should always produce the expected image.
            if (!std::filesystem::exists(output_png)) {
                throw std::runtime_error(
                    "circular_plot.py completed successfully but the expected output "
                    "file was not found: " + output_png.string()
                );
            }
        }

        /**
         * @brief Remove a temporary circular-plot working directory.
         *
         * Attempts to recursively delete the specified working directory and its
         * contents. Any filesystem errors are intentionally ignored so that cleanup
         * failures do not affect an otherwise successful plot-generation operation.
         *
         * @param work_dir Temporary working directory to remove.
         */
        void cleanupWorkDir(const std::filesystem::path &work_dir) {
            std::error_code error_code;
            std::filesystem::remove_all(work_dir, error_code);
        }
    } // namespace

    // Generate a circular coverage plot for a single component or query.
    void writeCircularPlotQuery(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &component_coverage,
        const std::string &component_name,
        const cdx::PosBp compo_length,
        const std::pair<cdx::PosBp, cdx::PosBp> query_bound,
        const PlotConfig &config
    ) {
        // Ensure the output directory exists.
        const std::filesystem::path image_dir =
                output_png.parent_path().empty()
                    ? std::filesystem::path(".")
                    : output_png.parent_path();
        std::filesystem::create_directories(image_dir);

        // Create a temporary working directory used to exchange data with the
        // Python rendering backend.
        const std::filesystem::path work_dir =
                image_dir / (output_png.stem().string() + "_circular_tmp");
        std::filesystem::create_directories(work_dir);

        const std::filesystem::path request_path = work_dir / "request.bin";

        // Precompute and package all plot data required for rendering.
        const CircularPlotPackage pkg = prepareCircularPlotPackage(
            component_coverage,
            component_name,
            compo_length,
            query_bound,
            config
        );

        // Serialize the plotting request in the binary format understood by
        // circular_plot.py.
        {
            BinaryRequestWriter writer(request_path);
            writeCommonHeader(writer, /*mode=*/0, config);
            writer.package(pkg);
            writer.finish();
        }

        // Execute the Python renderer and generate the final PNG.
        const std::filesystem::path script_path = resolveCircularPlotScript();
        runCircularPlotScriptOrThrow(script_path, request_path, output_png);

        // Remove intermediate files once rendering succeeds.
        cleanupWorkDir(work_dir);
    }

    // Generate a multi-panel circular coverage plot for all components.
    void writeCircularPlotGlobal(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &flat_bp_coverage,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names,
        const PlotConfig &config
    ) {
        // Validate that the component offsets describe a valid partition of the flattened coverage table.
        if (bp_component_offsets.size() < 2) {
            throw std::invalid_argument("bp_component_offsets must contain at least two offsets.");
        }
        if (bp_component_offsets.back() > static_cast<cdx::PosBp>(flat_bp_coverage.size())) {
            throw std::invalid_argument("The final component offset exceeds the coverage size.");
        }

        // Create the output and temporary working directories
        // used for exchanging data with the Python rendering backend.
        const std::filesystem::path image_dir =
                output_png.parent_path().empty() ? std::filesystem::path(".") : output_png.parent_path();
        std::filesystem::create_directories(image_dir);

        const std::filesystem::path work_dir = image_dir / (output_png.stem().string() + "_circular_tmp");
        std::filesystem::create_directories(work_dir);
        const std::filesystem::path request_path = work_dir / "request.bin";

        // Determine a suitable subplot layout for the number of components.
        const std::size_t component_count = bp_component_offsets.size() - 1;
        const GridLayout grid = chooseGlobalGraphGrid(component_count);

        // Component preprocessing is independent CPU-bound work and can therefore
        // be parallelized. Serialization remains single-threaded because the binary
        // output stream is shared.
        std::vector<CircularPlotPackage> packages(component_count);
        std::exception_ptr worker_exception;
        std::mutex exception_mutex;

        // Prepare one component package from its coverage subrange. Any exception
        // is captured and forwarded to the main thread after all workers complete.
        auto prepare_component = [&](const std::size_t component_id) {
            try {
                const auto start = static_cast<std::size_t>(bp_component_offsets[component_id]);
                const auto end = static_cast<std::size_t>(bp_component_offsets[component_id + 1]);

                if (end < start || end > flat_bp_coverage.size()) {
                    throw std::invalid_argument("bp_component_offsets is inconsistent with flat_bp_coverage size.");
                }

                // Extract the component-specific coverage range from the flattened table.
                const std::vector<cdx::Coverage> component_coverage(
                    flat_bp_coverage.begin() + static_cast<std::ptrdiff_t>(start),
                    flat_bp_coverage.begin() + static_cast<std::ptrdiff_t>(end)
                );

                const auto compo_length = static_cast<cdx::PosBp>(component_coverage.size());

                // Fall back to the component index if no name was provided.
                const std::string name = (component_id < component_names.size())
                                             ? component_names[component_id]
                                             : std::to_string(component_id);

                // Generate a full-component circular plot package.
                packages[component_id] = prepareCircularPlotPackage(
                    component_coverage, name, compo_length, {0, compo_length - 1}, config
                );

            // Preserve the first worker exception and rethrow it once all threads have completed.
            } catch (...) {
                const std::lock_guard<std::mutex> lock(exception_mutex);
                if (!worker_exception) {
                    worker_exception = std::current_exception();
                }
            }
        };

        {
            // Use at most one worker per component
            // and never exceed the number of hardware threads reported by the system.
            const unsigned int hardware_threads = std::thread::hardware_concurrency();
            const std::size_t thread_count = std::min(
                component_count,
                static_cast<std::size_t>(std::max(1u, hardware_threads))
            );

            // Components are assigned dynamically to workers through an atomic index to keep workload balanced.
            std::atomic<std::size_t> next_index{0};
            std::vector<std::thread> workers;
            workers.reserve(thread_count);

            // Claim the next available component and process it.
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

        // Throw any failure that occurred during parallel preprocessing.
        if (worker_exception) {
            std::rethrow_exception(worker_exception);
        }

        {
            // Serialize the complete multi-component plot request.
            BinaryRequestWriter writer(request_path);
            writeCommonHeader(writer, /*mode=*/1, config);

            // Figure size metadata
            writer.f64(config.figure_width);
            writer.f64(config.figure_height);

            // Global-layout metadata consumed by the Python renderer.
            writer.u32(static_cast<std::uint32_t>(grid.rows));
            writer.u32(static_cast<std::uint32_t>(grid.columns));

            // Serialize all component plot packages.
            writer.u32(static_cast<std::uint32_t>(component_count));
            for (const CircularPlotPackage &pkg: packages) {
                writer.package(pkg);
            }

            writer.finish();
        }

        // Render the complete plot grid using a single Python process.
        const std::filesystem::path script_path = resolveCircularPlotScript();
        runCircularPlotScriptOrThrow(script_path, request_path, output_png);

        // Remove temporary request files once rendering succeeds.
        cleanupWorkDir(work_dir);
    }
} // namespace output
