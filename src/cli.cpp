/**
 * @file cli.cpp
 * @brief Command-line interface parser and validation implementation for cdx_coverage.
 *
 * Provides argument parsing via CLI11 for input files (CDX index, GAM alignments),
 * execution modes (inspect vs. coverage calculation), thread dispatching configuration,
 * target region/component query parsing, and custom plot styling options.
 */

#include "cli.hpp"
#include <CLI/CLI.hpp>

#if defined(_OPENMP)
#include <omp.h>
#else
#include <thread>
#endif

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
    [[nodiscard]] std::string trim(const std::string &str) {
        static constexpr auto whitespace = " \t\n\r\v\f";

        const auto first = str.find_first_not_of(whitespace);
        if (first == std::string::npos) {
            return {};
        }

        const auto last = str.find_last_not_of(whitespace);
        return str.substr(first, last - first + 1);
    }

    [[nodiscard]] int get_hardware_threads() noexcept {
#if defined(_OPENMP)
        return omp_get_max_threads();
#else
        const unsigned int threads = std::thread::hardware_concurrency();
        return threads > 0 ? static_cast<int>(threads) : 1;
#endif
    }

    /**
     * @brief Custom parser for the `--query` command-line argument.
     *
     * Extracts a component identifier (component name or numerical component ID)
     * and an optional 0-based index range [start, end] from a string.
     *
     * @param text The raw input string supplied by the user (e.g., "chr1", "0", "chrM 1000:5000").
     * @return QuerySelection A populated struct containing the component target and optional range.
     * @throws std::invalid_argument If the input text does not match any valid query format.
     */
    QuerySelection parse_query(const std::string &text) {
        // Regex matching either:
        // Group 1: Component identifier (name or ID)
        // Group 2 & 3 (Optional): Range start and end coordinates separated by ':'
        static const std::regex query_regex(
            R"(^\s*([^\s]+)(?:\s+(-?\d+):(-?\d+))?\s*$)",
            std::regex::ECMAScript
        );

        std::smatch match;

        // Validate the input string against expected formats
        if (!std::regex_match(text, match, query_regex)) {
            throw std::invalid_argument(
                "Invalid query format: '" + text + "'.\n"
                "Expected formats:\n"
                "  chr1           (Entire component by name)\n"
                "  0              (Entire component by compoID)\n"
                "  chrM 1000:5000 (Range selection by name)\n"
                "  5 100:1000     (Range selection by compoID)"
            );
        }

        QuerySelection sel;
        // Capture mandatory component name or ID (Group 1)
        sel.component = match[1].str();

        // Check if an optional range was captured (Groups 2 and 3)
        if (match[2].matched && match[3].matched) {
            QueryRange r;
            // Parse string coordinates to 64-bit signed integers
            r.start = std::stoll(match[2].str());
            r.end = std::stoll(match[3].str());
            sel.range = r;
        }

        return sel;
    }


    /**
 * @brief Parses and normalizes the graph component type from a string.
 *
 * Trims input whitespace, converts characters to lower case, and maps
 * recognized flags to the corresponding `ComponentType` enum value.
 *
 * @param value Raw string input representing component type (e.g., "linear", "l", "CIRCULAR", "c").
 * @return ComponentType Enum value indicating `ComponentType::Linear` or `ComponentType::Circular`.
 * @throws std::invalid_argument If the normalized string is not a valid component type alias.
 */
    ComponentType parse_component_type(std::string value) {
        value = trim(value); // Strip leading/trailing whitespace

        // Convert string to lower case for case-insensitive matching
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });

        // Map recognized aliases to enum values
        if (value == "l" || value == "linear") return ComponentType::Linear;
        if (value == "c" || value == "circular") return ComponentType::Circular;

        throw std::invalid_argument("Component type must be 'linear'/'l' or 'circular'/'c'.");
    }


    /**
     * @brief Parses and normalizes the coverage-precision mode from a string.
     *
     * Trims input whitespace, converts characters to lower case, and maps
     * recognized flags to the corresponding `CoveragePrecision` enum value.
     * Mirrors `parse_component_type()` above.
     *
     * @param value Raw string input (e.g., "base", "b", "NODE", "n").
     * @return CoveragePrecision Enum value indicating `CoveragePrecision::Node`
     *         or `CoveragePrecision::Base`.
     * @throws std::invalid_argument If the normalized string is not a valid alias.
     */
    CoveragePrecision parse_coverage_precision(std::string value) {
        value = trim(value); // Strip leading/trailing whitespace

        // Convert string to lower case for case-insensitive matching
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });

        // Map recognized aliases to enum values
        if (value == "n" || value == "node") return CoveragePrecision::Node;
        if (value == "b" || value == "base") return CoveragePrecision::Base;

        throw std::invalid_argument("Coverage precision must be 'node'/'n' or 'base'/'b'.");
    }


    /**
     * @brief Parses plot figure dimensions from a string formatted as WIDTHxHEIGHT.
     *
     * Supports standard 'x' or 'X' separators as well as the UTF-8 multiplication symbol ('×').
     *
     * @param value String containing figure dimensions (e.g., "7x4.5", "10X8").
     * @return std::pair<double, double> A pair containing {width, height} in inches.
     * @throws std::invalid_argument If the string format is invalid or dimensions are <= 0.
     */
    std::pair<double, double> parse_fig_size(const std::string &value) {
        // Regex matching:
        // Group 1: Width (positive integer or floating-point)
        // Separator: 'x', 'X', or UTF-8 '×' (\xC3\x97)
        // Group 2: Height (positive integer or floating-point)
        static const std::regex size_regex(
            R"(\s*(\d+(?:\.\d+)?)\s*[xX\xC3\x97]\s*(\d+(?:\.\d+)?)\s*)",
            std::regex::ECMAScript
        );

        std::smatch match;
        // Validate string layout against the regex pattern
        if (!std::regex_match(value, match, size_regex)) {
            throw std::invalid_argument("Figure size must use WIDTHxHEIGHT format, e.g., '7x4.5'.");
        }

        // Convert matched numerical string components to double values
        const double w = std::stod(match[1].str());
        const double h = std::stod(match[2].str());

        // Ensure figure dimensions are strictly positive
        if (w <= 0.0 || h <= 0.0) {
            throw std::invalid_argument("Figure width and height must be strictly positive.");
        }

        return {w, h};
    }

    /**
 * @brief Parses and validates a hexadecimal color code string.
 *
 * Trims leading/trailing whitespace, checks against standard hex color formats
 * (including short-hand, RGB, RGBA, and ARGB variants), and normalizes the
 * output string to uppercase.
 *
 * @param value Raw color string provided by the user (e.g., "#1e3a8a", "#fff", "#1E3A8AFF").
 * @return std::string Standardized, uppercase hexadecimal color string (e.g., "#1E3A8A").
 * @throws std::invalid_argument If the string does not match valid hexadecimal color patterns.
 */
    std::string parse_hex_color(std::string value) {
        value = trim(value); // Strip leading and trailing whitespace


        // Regex matching valid hex color formats starting with '#':
        // - 3 or 4 digits: #RGB or #RGBA (shorthand)
        // - 6 digits:      #RRGGBB
        // - 8 digits:      #RRGGBBAA
        static const std::regex color_regex(
            R"(#(?:[0-9A-Fa-f]{3,4}|[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8}))",
            std::regex::ECMAScript
        );

        if (!std::regex_match(value, color_regex)) {
            throw std::invalid_argument("Color must be hexadecimal, for example #1E3A8A.");
        }

        // Convert hex letters (a-f) to uppercase for uniform color code representation
        std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char c) {
                           return static_cast<char>(std::toupper(c));
                       }
        );

        return value;
    }
} // namespace anonyme

// @brief Parses command-line options and validates inputs using CLI11.
CliArgs parse_args(const int argc, char **argv) {
    CliArgs args;

    // Initialize CLI11 application with description and custom usage format
    CLI::App app{"Calculate per-node coverage from a GAM file using a CDX index."};
    app.usage("cdx_coverage <CDX> [GAM] [OPTIONS]");

    // ============================================================
    // 1. POSITIONAL INPUT FILES
    // ============================================================
    // CDX file is mandatory and must exist on disk
    app.add_option("cdx", args.cdx_file, "Path to the binary CDX graph index.")
            ->required()
            ->check(CLI::ExistingFile);

    // GAM file is optional at initial parse (validated post-parse based on mode)
    app.add_option("gam", args.gam_file, "Path to the GAM alignment file.")
            ->check(CLI::ExistingFile);

    // ============================================================
    // 2. QUERY SELECTION GROUP
    // ============================================================
    auto *group_query = app.add_option_group("QUERY");

    std::string query_raw;
    group_query->add_option(
        "-q,--query",
        query_raw,
        "Scope of the coverage calculation (0-based coordinates).\n"
        "Formats: COMPONENT or 'COMPONENT START:END' (accepts name or ComponentID)\n"
        "Examples: -q chr1, -q 0, -q \"chr1 1000:5000\""
    );

    std::string comp_type_raw;
    group_query->add_option(
        "-c,--component-type",
        comp_type_raw,
        "Graph coordinate mapping structure:\n"
        "'linear'/'l' or 'circular'/'c'.\n"
        "Default: linear."
    );

    std::string coverage_precision_raw;
    group_query->add_option(
        "-p,--coverage-precision",
        coverage_precision_raw,
        "How precisely to compute coverage:\n"
        "'base'/'b' = per base pair (accurate, default).\n"
        "'node'/'n' = per whole node (faster, uses less memory).\n"
        "Use 'node' for very large or deep GAM files.\n"
        "Default: base."
    );

    // ============================================================
    // 3. PERFORMANCE / THREADS GROUP
    // ============================================================
    std::string worker_threads_arg = "auto";
    std::string decompression_threads_arg = "auto";

    app.add_option(
        "-t,--worker-threads",
        worker_threads_arg,
        "Number of threads used for computation (positive integer or 'auto').\n"
        "Default: half the machine threads."
    );

    app.add_option(
        "-T,--decompression-threads",
        decompression_threads_arg,
        "Number of threads used for decompression of the GAM file (positive integer or 'auto').\n "
        "Default: half the machine threads"
    );

    // ============================================================
    // 4. INSPECTION GROUP
    // ============================================================
    std::string inspect_val;
    // Allow flag to be called without a value (0 args) or with a component name (1 arg)
    auto *opt_inspect = app.add_option(
        "-i,--inspect",
        inspect_val,
        "Inspect CDX contents and exit.\n"
        "Without a value, display a summary of all components.\n"
        "With a component name/ID, display only the selected component."
    )->type_size(0, 1);

    // ============================================================
    // 5. OUTPUT CONFIGURATION GROUP
    // ============================================================
    auto *group_output = app.add_option_group("OUTPUT");

    group_output->add_option("-o,--output", args.output_directory, "Directory to save output files. Default: '.'");
    group_output->add_flag("--no-graph", args.no_graph, "Skip coverage graph generation.");
    group_output->add_flag("--no-stats", args.no_stats, "Skip writing summary statistics file.");
    group_output->add_flag("--no-table", args.no_table, "Skip writing per-node TSV table file.");

    // ============================================================
    // 6. GRAPH RENDERING OPTIONS
    // ============================================================
    auto *group_graph = app.add_option_group("GRAPH");

    int log_base_val = 10;
    // Optional value flag: Defaults to base 10 if flag is present without an explicit integer
    auto *opt_log = group_graph->add_option(
        "--log", log_base_val, "Use logarithmic coverage scale. Base defaults to 10 if omitted."
    )->type_size(0, 1)->check(CLI::Range(2, 10000))->default_str("10");

    group_graph->add_option(
                "--smoothing",
                args.smoothing,
                "Moving-average window fraction [0.0, 1.0].\n"
                "Default: 0.01.")
            ->check(CLI::Range(0.0, 1.0));

    group_graph->add_option(
        "--max-point,--max-points",
        args.max_plot_points,
        "Maximum points passed to plotting backend.\n"
        "Use 0 to disable downsampling and plot at full resolution.\n"
        "Default: 10000.");

    group_graph->add_option(
                "--dpi",
                args.dpi,
                "Output graph resolution in DPI.\n"
                "Default: 300.")
            ->check(CLI::PositiveNumber);

    std::string fig_size_raw;
    group_graph->add_option(
        "--fig-size",
        fig_size_raw,
        "Figure dimensions in inches (WIDTHxHEIGHT), e.g., '7x4.5'.");

    std::string line_color_raw;
    group_graph->add_option(
        "--color-line",
        line_color_raw,
        "Hexadecimal color for coverage line.\n"
        "Default: #1E3A8A.");

    std::string fill_color_raw;
    group_graph->add_option(
        "--color-filling",
        fill_color_raw,
        "Hexadecimal color for coverage fill.\n"
        "Default: #93C5FD.");

    // Parse raw CLI inputs; CLI11 automatically prints error/help/version text
    // via app.exit(). A zero exit code means --help/--version was requested,
    // which is not an error and still terminates the process immediately (no
    // caller ever needs to observe/test that path beyond its exit code). Any
    // other exit code represents a genuine parse failure, which is thrown
    // instead so it flows through the same std::exception handling as every
    // other validation error below (and so parse_args stays testable without
    // forking a subprocess).
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        const int exit_code = app.exit(e);
        if (exit_code == 0) {
            std::exit(exit_code);
        }
        throw std::runtime_error("Argument parsing failed (see usage above).");
    }

    // ============================================================
    // POST-PARSE VALIDATION AND DYNAMIC CONFIGURATION
    // ============================================================

    // 1. Process Inspection Mode Settings
    if (*opt_inspect) {
        args.inspect.enabled = true;
        if (!inspect_val.empty()) {
            args.inspect.component = inspect_val;
        }
    }

    // Enforce mandatory GAM file presence if NOT running in inspection mode
    if (!args.inspectMode() && args.gam_file.empty()) {
        throw std::runtime_error(
            "GAM file argument is required unless running in inspect mode (-i/--inspect)."
        );
    }

    // 2. Dynamic Thread Resolution based on Available Hardware
    const int machine_threads = get_hardware_threads();

    try {
        // Resolve worker thread count
        if (worker_threads_arg == "auto") {
            args.worker_threads = machine_threads;
        } else {
            const int requested = std::stoi(worker_threads_arg);
            if (requested <= 0) {
                throw std::invalid_argument("Worker threads must be positive or 'auto'.");
            }
            // Cap manually requested threads to available hardware maximum
            args.worker_threads = std::min(requested, machine_threads);
        }

        // Resolve decompression thread count (defaults to half of hardware threads)
        if (decompression_threads_arg == "auto") {
            args.decompression_threads = std::max(1, machine_threads / 2);
        } else {
            const int requested = std::stoi(decompression_threads_arg);
            if (requested <= 0) {
                throw std::invalid_argument("Decompression threads must be positive or 'auto'.");
            }
            // Cap manually requested threads to available hardware maximum
            args.decompression_threads = std::min(requested, machine_threads);
        }
    } catch (const std::exception &e) {
        throw std::runtime_error(std::string("Error parsing thread arguments: ") + e.what());
    }

    // 3. Parse Custom Data Formats (Query, Component Type, Hex Colors)
    if (!query_raw.empty()) {
        try {
            args.query = parse_query(query_raw);
        } catch (const std::exception &e) {
            throw std::runtime_error(std::string("Error parsing --query: ") + e.what());
        }
    }

    if (!comp_type_raw.empty()) {
        try {
            args.component_type = parse_component_type(comp_type_raw);
        } catch (const std::exception &e) {
            throw std::runtime_error(std::string("Error parsing --component-type: ") + e.what());
        }
    }

    if (!coverage_precision_raw.empty()) {
        try {
            args.coverage_precision = parse_coverage_precision(coverage_precision_raw);
        } catch (const std::exception &e) {
            throw std::runtime_error(std::string("Error parsing --coverage-precision: ") + e.what());
        }
    }

    if (*opt_log) {
        args.log_base = log_base_val;
    }

    if (!line_color_raw.empty()) {
        try {
            args.line_color = parse_hex_color(line_color_raw);
        } catch (const std::exception &e) {
            throw std::runtime_error(std::string("Error parsing --color-line: ") + e.what());
        }
    }

    if (!fill_color_raw.empty()) {
        try {
            args.fill_color = parse_hex_color(fill_color_raw);
        } catch (const std::exception &e) {
            throw std::runtime_error(std::string("Error parsing --color-filling: ") + e.what());
        }
    }

    // Ensure at least one output mode remains active
    if (args.no_graph && args.no_stats && args.no_table) {
        throw std::runtime_error("All outputs are disabled (--no-graph, --no-stats, --no-table).");
    }

    // 4. Determine Output Figure Dimensions
    if (!fig_size_raw.empty()) {
        // Parse explicit user dimensions (e.g., "7x4.5")
        try {
            args.custom_figure_size = parse_fig_size(fig_size_raw);
            args.figure_width = args.custom_figure_size->first;
            args.figure_height = args.custom_figure_size->second;
        } catch (const std::exception &e) {
            throw std::runtime_error(std::string("Error parsing --fig-size: ") + e.what());
        }
    } else {
        // Compute default canvas aspect ratios based on query presence and graph topology
        const bool is_circular = (args.component_type == ComponentType::Circular);
        const bool has_query = args.query.has_value();

        if (!has_query && !is_circular) {
            args.figure_width = 5.5; // Standard linear overview
            args.figure_height = 3.5;
        } else if (!has_query && is_circular) {
            args.figure_width = 4.0; // Square ratio for whole circular genome plots
            args.figure_height = 4.0;
        } else if (has_query && !is_circular) {
            args.figure_width = 7.0; // Wide aspect for detailed linear region focus
            args.figure_height = 4.5;
        } else {
            args.figure_width = 7.0; // Large square for detailed circular sub-region
            args.figure_height = 7.0;
        }
    }

    return args;
}
