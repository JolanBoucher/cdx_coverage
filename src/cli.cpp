#include "cli.hpp"
#include <CLI/CLI.hpp>

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

    // Custom parser pour --query (nom de composante ou ID + plage optionnelle)
    QuerySelection parse_query(const std::string &text) {
        static const std::regex query_regex(
            R"(^\s*([^\s]+)(?:\s+(-?\d+):(-?\d+))?\s*$)",
            std::regex::ECMAScript
        );

        std::smatch match;
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
        sel.component = match[1].str();

        if (match[2].matched && match[3].matched) {
            QueryRange r;
            r.start = std::stoll(match[2].str());
            r.end = std::stoll(match[3].str());
            sel.range = r;
        }

        return sel;
    }

    ComponentType parse_component_type(std::string value) {
        value = trim(value);

        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
                           return static_cast<char>(std::tolower(c));
                       }
        );

        if (value == "l" || value == "linear") return ComponentType::Linear;
        if (value == "c" || value == "circular") return ComponentType::Circular;
        throw std::invalid_argument("Component type must be 'linear'/'l' or 'circular'/'c'.");
    }

    std::pair<double, double> parse_fig_size(const std::string &value) {
        static const std::regex size_regex(
            R"(\s*(\d+(?:\.\d+)?)\s*[xX\xC3\x97]\s*(\d+(?:\.\d+)?)\s*)",
            std::regex::ECMAScript
        );

        std::smatch match;
        if (!std::regex_match(value, match, size_regex)) {
            throw std::invalid_argument("Figure size must use WIDTHxHEIGHT format, e.g., '7x4.5'.");
        }

        const double w = std::stod(match[1].str());
        const double h = std::stod(match[2].str());

        if (w <= 0.0 || h <= 0.0) {
            throw std::invalid_argument("Figure width and height must be strictly positive.");
        }

        return {w, h};
    }

    std::string parse_hex_color(std::string value) {
        value = trim(value);
        static const std::regex color_regex(
            R"(#(?:[0-9A-Fa-f]{3,4}|[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8}))",
            std::regex::ECMAScript
        );

        if (!std::regex_match(value, color_regex)) {
            throw std::invalid_argument("Color must be hexadecimal, for example #1E3A8A.");
        }

        std::transform(
            value.begin(),
            value.end(),
            value.begin(),
            [](unsigned char c) {
                return static_cast<char>(std::toupper(c));
            }
        );
        return value;
    }
} // namespace anonyme

CliArgs parse_args(const int argc, char **argv) {
    CliArgs args;

    CLI::App app{"Calculate per-node coverage from a GAM file using a CDX index."};
    app.usage("cdx_coverage <CDX> [GAM] [OPTIONS]");

    // ============================================================
    // 1. POSITIONAL INPUT FILES
    // ============================================================
    app.add_option("cdx", args.cdx_file, "Path to the binary CDX graph index.")
            ->required()
            ->check(CLI::ExistingFile);

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
        "Formats: COMPONENT or COMPONENT:START-END (accepts name or CID)\n"
        "Examples: -q chr1, -q 0, -q chr1 1000-5000"
    );

    std::string comp_type_raw;
    group_query->add_option(
        "-t,--component-type",
        comp_type_raw,
        "Graph coordinate mapping structure: 'linear'/'l' or 'circular'/'c'. Default: linear."
    );

    // ============================================================
    // 3. INSPECTION GROUP
    // ============================================================
    std::string inspect_val;
    auto *opt_inspect = app.add_option(
        "-i,--inspect",
        inspect_val,
        "Inspect CDX contents and exit.\n"
        "Without a value, display a summary of all components.\n"
        "With a component name/ID, display only the selected component."
    )->type_size(0, 1);

    // ============================================================
    // 4. OUTPUT CONFIGURATION GROUP
    // ============================================================
    auto *group_output = app.add_option_group("OUTPUT");

    group_output->add_option("-o,--output", args.output_directory, "Directory to save output files. Default: '.'");
    group_output->add_flag("--no-graph", args.no_graph, "Skip coverage graph generation.");
    group_output->add_flag("--no-stats", args.no_stats, "Skip writing summary statistics file.");
    group_output->add_flag("--no-table", args.no_table, "Skip writing per-node TSV table file.");

    // ============================================================
    // 5. GRAPH RENDERING OPTIONS
    // ============================================================
    auto *group_graph = app.add_option_group("GRAPH");

    int log_base_val = 10;
    auto *opt_log = group_graph->add_option(
        "--log", log_base_val, "Use logarithmic coverage scale. Base defaults to 10 if omitted."
    )->type_size(0, 1)->check(CLI::Range(2, 10000))->default_str("10");

    group_graph->add_option("--smoothing", args.smoothing, "Moving-average window fraction [0.0, 1.0]. Default: 0.01.")
            ->check(CLI::Range(0.0, 1.0));

    group_graph->add_option("--max-point,--max-points", args.max_plot_points,
                            "Maximum points passed to plotting backend (0 to disable). Default: 10000.");
    group_graph->add_option("--dpi", args.dpi, "Output graph resolution in DPI. Default: 300.")
            ->check(CLI::PositiveNumber);

    std::string fig_size_raw;
    group_graph->add_option("--fig-size", fig_size_raw, "Figure dimensions in inches (WIDTHxHEIGHT), e.g., '7x4.5'.");

    // Buffers temporaires pour éviter l'écrasement par chaîne vide de CLI11
    std::string line_color_raw;
    group_graph->add_option("--color-line", line_color_raw, "Hexadecimal color for coverage line. Default: #1E3A8A.");

    std::string fill_color_raw;
    group_graph->add_option("--color-filling", fill_color_raw,
                            "Hexadecimal color for coverage fill. Default: #93C5FD.");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        std::exit(app.exit(e));
    }

    // --- Post-parse validation logic ---

    if (*opt_inspect) {
        args.inspect.enabled = true;
        if (!inspect_val.empty()) {
            args.inspect.component = inspect_val;
        }
    }

    if (!args.inspectMode() && args.gam_file.empty()) {
        std::cerr << "Error: GAM file argument is required unless running in inspect mode (-i/--inspect).\n";
        std::exit(1);
    }

    if (!query_raw.empty()) {
        try {
            args.query = parse_query(query_raw);
        } catch (const std::exception &e) {
            std::cerr << "Error parsing --query: " << e.what() << "\n";
            std::exit(1);
        }
    }

    if (!comp_type_raw.empty()) {
        try {
            args.component_type = parse_component_type(comp_type_raw);
        } catch (const std::exception &e) {
            std::cerr << "Error parsing --component-type: " << e.what() << "\n";
            std::exit(1);
        }
    }

    if (*opt_log) {
        args.log_base = log_base_val;
    }

    // Assignation conditionnelle des couleurs
    if (!line_color_raw.empty()) {
        try {
            args.line_color = parse_hex_color(line_color_raw);
        } catch (const std::exception &e) {
            std::cerr << "Error parsing --color-line: " << e.what() << "\n";
            std::exit(1);
        }
    }

    if (!fill_color_raw.empty()) {
        try {
            args.fill_color = parse_hex_color(fill_color_raw);
        } catch (const std::exception &e) {
            std::cerr << "Error parsing --color-filling: " << e.what() << "\n";
            std::exit(1);
        }
    }

    if (args.no_graph && args.no_stats && args.no_table) {
        std::cerr << "Error: All outputs are disabled (--no-graph, --no-stats, --no-table).\n";
        std::exit(1);
    }

    if (!fig_size_raw.empty()) {
        try {
            args.custom_figure_size = parse_fig_size(fig_size_raw);
            args.figure_width = args.custom_figure_size->first;
            args.figure_height = args.custom_figure_size->second;
        } catch (const std::exception &e) {
            std::cerr << "Error parsing --fig-size: " << e.what() << "\n";
            std::exit(1);
        }
    } else {
        const bool is_circular = (args.component_type == ComponentType::Circular);
        const bool has_query = args.query.has_value();

        if (!has_query && !is_circular) {
            args.figure_width = 5.5;
            args.figure_height = 3.5;
        } else if (!has_query && is_circular) {
            args.figure_width = 4.0;
            args.figure_height = 4.0;
        } else if (has_query && !is_circular) {
            args.figure_width = 7.0;
            args.figure_height = 4.5;
        } else {
            args.figure_width = 7.0;
            args.figure_height = 7.0;
        }
    }

    return args;
}
