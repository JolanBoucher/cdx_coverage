#ifndef CDX_COVERAGE_CLI_HPP
#define CDX_COVERAGE_CLI_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>

struct QueryRange {
    std::int64_t start;
    std::int64_t end;
};

struct QuerySelection {
    std::string component; // Contient soit le nom ("chr1"), soit l'ID numérique ("0")
    std::optional<QueryRange> range;
};

enum class ComponentType {
    Linear,
    Circular
};

struct InspectOptions {
    bool enabled = false;
    std::optional<std::string> component; // Supporte aussi le nom ou l'ID ("chr1" ou "0")
};

struct CliArgs {
    // Positional input files
    std::string cdx_file;
    std::string gam_file;

    // Component Query Selection
    std::optional<QuerySelection> query;
    ComponentType component_type = ComponentType::Linear;

    // CDX Inspection Mode
    InspectOptions inspect;

    // threads configuration
    int worker_threads;
    int decompression_threads;

    // Output Configuration
    std::string output_directory = ".";
    bool no_graph = false;
    bool no_stats = false;
    bool no_table = false;

    // Graph Rendering Options
    std::optional<int> log_base;
    double smoothing = 0.01;
    std::size_t max_plot_points = 10000;
    int dpi = 300;
    std::optional<std::pair<double, double>> custom_figure_size;

    // Computed figure dimensions
    double figure_width = 5.5;
    double figure_height = 3.5;

    std::string line_color = "#1E3A8A";
    std::string fill_color = "#93C5FD";

    // --- Helpers ---

    [[nodiscard]] bool inspectMode() const noexcept {
        return inspect.enabled;
    }

    [[nodiscard]] bool generateGraph() const noexcept {
        return !no_graph;
    }

    [[nodiscard]] bool generateStats() const noexcept {
        return !no_stats;
    }

    [[nodiscard]] bool generateTable() const noexcept {
        return !no_table;
    }
};

/**
 * @brief Parses and validates command-line arguments for the coverage tool.
 *
 * Configures CLI11 options, positional arguments, option groups, and validation
 * rules for calculating GAM coverage using a CDX index. Handles thread auto-scaling,
 * component query parsing, logging parameters, and dynamic figure dimensions.
 *
 * @param argc Count of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return CliArgs Populated and validated argument structure ready for program execution.
 *
 * @throws std::runtime_error If argument parsing or post-parse validation fails for any reason
 *         (malformed/missing/out-of-range option values, missing required GAM file outside
 *         inspect mode, all outputs disabled, etc.). The exception message is ready to display
 *         to the user as-is.
 *
 * @note `--help`/`--version` are handled by calling std::exit(0) directly after printing their
 *       text (CLI11's own behavior) - they are not reported as errors and do not throw.
 */
CliArgs parse_args(int argc, char** argv);

#endif // CDX_COVERAGE_CLI_HPP