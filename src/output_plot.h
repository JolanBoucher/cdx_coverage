/**
 * @file output_plot.h
 * @brief Coverage-plot preparation and rendering utilities.
 *
 * This module provides the data structures and algorithms used to generate
 * coverage visualizations. It includes routines for smoothing and
 * downsampling coverage data, computing plot scales and tick marks,
 * preparing render-ready plotting packages, and generating both linear and
 * circular coverage plots.
 *
 * Plot preparation is performed in C++, producing compact data structures
 * that can be consumed directly by rendering backends. The module supports
 * optional logarithmic scaling, component-level and global visualizations,
 * and handling of circular genomic coordinates.
 */

#ifndef CDX_COVERAGE_OUTPUT_PLOT_H
#define CDX_COVERAGE_OUTPUT_PLOT_H

#include "cdx_types.h"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cdx {
    /**
     * @brief Render-ready data package for a linear coverage plot.
     *
     * Contains all precomputed information required to render a single linear
     * coverage track, including plotting coordinates, axis scaling parameters,
     * and tick metadata. Coverage preprocessing (smoothing, downsampling,
     * scaling, etc.) has already been performed before this package is
     * consumed by the rendering backend.
     */
    struct LinearPlotPackageBin {
        // Component and query-region metadata.
        std::string component_name;
        std::size_t query_start{0};
        std::size_t query_end{0};

        // Y-axis scaling configuration.
        bool logarithmic{false};
        int log_base{0};
        double y_upper_limit{1.0};

        // Regularly spaced genomic X coordinates. Position i is located at x_start + i * x_step.
        double x_start{0.0};
        double x_step{1.0};

        // Coverage values corresponding to the plotting grid.
        std::vector<double> y;

        // Y-axis tick positions and labels.
        std::vector<double> tick_positions;
        std::vector<std::string> tick_labels;
    };
} // namespace cdx

namespace output {

    // Small class to choose between two graph config
    enum class Topology {
        Linear,
        Circular
    };

    /**
      * @brief Configuration parameters controlling coverage-plot generation.
      *
      * Defines preprocessing, scaling, rendering, and styling options shared
      * by the plotting backends. These parameters influence smoothing,
      * downsampling, axis scaling, figure dimensions, and plot appearance.
      */
    struct PlotConfig {
        // Coverage preprocessing parameters.
        double smoothing = 0.01;
        std::size_t max_plot_points = 10000;

        // Rendering resolution and optional logarithmic scaling.
        int dpi = 300;
        std::optional<int> log_base;

        // Figure dimensions in inches.
        double figure_width = 7.0;
        double figure_height = 4.5;

        // Plot styling colors.
        std::string line_color = "#1E3A8A";
        std::string fill_color = "#93C5FD";
    };

    /**
     * @brief Coverage data prepared for plotting.
     *
     * Stores plotting coordinates, corresponding coverage values, and the
     * smoothing/downsampling window used to generate them.
     */
    struct PlotData {
        std::vector<std::size_t> x;
        std::vector<double> y;
        std::size_t window_size = 1;
    };

    /**
     * @brief Grid dimensions used to arrange multiple plots.
     *
     * Stores the number of rows and columns in a multi-panel plot layout.
     */
    struct GridLayout {
        std::size_t rows = 1;
        std::size_t columns = 1;
    };

    /**
     * @brief Tick-mark positions and axis limits for a linear coverage scale.
     *
     * Stores the coverage values used as axis ticks and the maximum axis limit
     * selected for plotting.
     */
    struct CoverageTicks {
        std::vector<double> values;
        double upper_limit = 1.0;
    };

    /**
     * @brief Tick-mark metadata for logarithmically scaled coverage plots.
     *
     * Stores coverage tick values in both their original (raw) coverage space
     * and their transformed display space, along with the upper limit of the
     * displayed logarithmic axis.
     */
    struct LogCoverageTicks {
        std::vector<double> raw_values;
        std::vector<double> display_values;
        double display_upper_limit = 1.0;
    };

     //----------- Shared numerical functions. ------------//

    /**
     * @brief Computes the optimal grid layout (rows and columns) for a given number of global graph components.
     *
     * For component counts between 1 and 9, hand-tuned layouts are used for ideal visual
     * aspect ratios. For 10 to 30 components, it dynamically calculates a balanced grid
     * based on the square root of the count.
     *
     * @param component_count The total number of graph components to arrange (must be between 1 and 30).
     * @return GridLayout A structure containing the calculated row and column dimensions.
     *
     * @throws std::invalid_argument If @p component_count is less than 1 or greater than 30.
     */
    [[nodiscard]]
    GridLayout chooseGlobalGraphGrid(
        std::size_t component_count
    );

    /**
     * @brief Overload that converts raw integer coverage data to double precision and delegates processing.
     *
     * This function casts each element of the `cdx::Coverage` vector into a `double` and forwards the
     * dataset to `prepareCoverageForPlotImpl`, which performs the primary downsampling and
     * moving-average smoothing calculations.
     */
    [[nodiscard]]
    PlotData prepareCoverageForPlot(
        const std::vector<cdx::Coverage> &coverage,
        double smoothing,
        std::size_t max_plot_points,
        Topology topology
    );

    /**
     * @brief Overload for pre-transformed double values (e.g., NaN-masked circular backends),
     *        delegating directly to `prepareCoverageForPlotImpl`.
     *
     * @param coverage Vector of double-precision coverage values (may contain NaN/sentinels).
     * @param smoothing Moving-average window size expressed as a fraction of total points [0.0, 1.0].
     * @param max_plot_points Maximum number of points allowed in the downsampled output (0 = full resolution).
     * @param topology Graph topology mode (Topology::Linear or Topology::Circular).
     * @return PlotData Processed structure containing downsampled X coordinates, smoothed Y values, and window size.
     */
    [[nodiscard]]
    PlotData prepareCoverageForPlot(
        const std::vector<double> &coverage,
        double smoothing,
        std::size_t max_plot_points,
        Topology topology
    );


    /**
     * @brief Computes a nice human-readable axis tick positions and an upper limit for linear coverage plots.
     *
     * This function calculates the best linear axis tick intervals based on the maximum coverage value
     * and a target tick count. It applies standard round number scaling algorithms (using multipliers
     * of 1, 2, 2.5, 3, 4, 5, and 10) to determine clean step sizes and upper boundary limits for plotting.
     *
     * @param maximum_coverage The peak coverage value encountered in the dataset.
     * @param target_tick_count The desired approximate number of tick intervals (must be at least 2).
     * @return CoverageTicks Structure containing calculated tick values and the upper axis limit.
     *
     * @throws std::invalid_argument If @p maximum_coverage is non-finite or @p target_tick_count is less than 2.
     */
    [[nodiscard]]
    CoverageTicks calculateCoverageTicks(
        double maximum_coverage,
        std::size_t target_tick_count = 4
    );

    /**
     * @brief Transforms a vector of coverage values onto a custom logarithmic scale.
     *
     * This function applies a modified log transformation to coverage data. Values less than
     * or equal to 1.0 are kept linear (preserving a smooth transition around zero and one),
     * while values greater than 1.0 are scaled logarithmically using the specified base.
     *
     * @param values Vector of raw double-precision coverage values to transform.
     * @param log_base The base of the logarithm used for transformation (must be strictly greater than 1).
     * @return std::vector<double> Vector containing the transformed coverage values.
     *
     * @throws std::invalid_argument If @p log_base is less than or equal to 1, or if any coverage value is negative.
     */
    [[nodiscard]]
    std::vector<double> logTransformCoverage(
        const std::vector<double> &values,
        int log_base
    );

    /**
     * @brief Computes raw and log-transformed axis tick positions for logarithmic coverage plots.
     *
     * This function calculates appropriate tick marks along a logarithmic scale based on the maximum
     * coverage value and the chosen log base. It generates raw power-of-base values, transforms them
     * for plotting, and determines the upper limit bound for display scaling.
     *
     * @param maximum_coverage The peak coverage value encountered in the dataset.
     * @param log_base The base used for logarithmic scaling (must be strictly greater than 1).
     * @return LogCoverageTicks Structure containing raw ticks, log-transformed display ticks, and the upper display limit.
     *
     * @throws std::invalid_argument If @p log_base is less than or equal to 1.
     */
    [[nodiscard]]
    LogCoverageTicks logCoverageTicks(
        double maximum_coverage,
        int log_base
    );

    /**
     * @brief Prepares and formats a linear coverage plot package for rendering and visualization.
     *
     * This function filters out-of-query sentinel values, applies downsampling and moving-average
     * smoothing via `prepareCoverageForPlot`, configures X-axis coordinate scaling steps, and handles
     * both linear and logarithmic Y-axis transformations (including tick placement and formatting labels).
     *
     * @param coverage Vector of raw coverage values for the component or region.
     * @param component_name Name or identifier of the graph component being processed.
     * @param offset Starting coordinate offset of the query region.
     * @param config Configuration parameters governing plotting options (smoothing, scaling, downsampling, log base).
     * @return cdx::LinearPlotPackageBin Populated package containing coordinate metadata, Y values, ticks, and labels.
     */
    [[nodiscard]]
    cdx::LinearPlotPackageBin prepareLinearPlotPackage(
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name,
        std::size_t offset,
        const PlotConfig &config
    );

     //---------- Linear and circular backend --------//

    /**
     * @brief Renders a linear plot for a single queried region within a component and saves it to a PNG file.
     *
     * @param output_png Target PNG output path.
     * @param coverage Vector of base-pair coverage data for the region.
     * @param component_name Name of the component sequence.
     * @param offset Start base-pair coordinate offset.
     * @param config Plot configuration options.
     */
    void writeLinearPlotQuery(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name,
        std::size_t offset,
        const PlotConfig &config
    );

    /**
     * @brief Generate a circular coverage plot for a single component or query.
     *
     * Prepares a render-ready CircularPlotPackage, serializes it into the
     * binary request format consumed by circular_plot.py, invokes the Python
     * renderer, and writes the resulting PNG image.
     *
     * A temporary working directory is created alongside the output image to
     * hold intermediate files required by the rendering pipeline and is removed
     * after successful completion.
     *
     * @param output_png Output PNG image path.
     * @param component_coverage Coverage values for the entire component.
     * @param component_name Component name displayed in the plot.
     * @param compo_length Component length in base pairs.
     * @param query_bound Inclusive query interval within the component.
     * @param config Plot rendering configuration.
     *
     * @throws std::runtime_error If request generation, script resolution,
     *         Python execution, or image creation fails.
     * @throws std::invalid_argument If the plotting inputs are inconsistent.
     * @throws std::out_of_range If the query bounds fall outside the component.
     */
    void writeCircularPlotQuery(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &component_coverage,
        const std::string &component_name,
        cdx::PosBp compo_length,
        std::pair<cdx::PosBp, cdx::PosBp> query_bound,
        const PlotConfig &config
    );

    /**
     * @brief Generates a multi-panel grid plot for global coverage across multiple components using multithreading.
     *
     * @param output_png Path where PNG output file will be written.
     * @param flat_bp_coverage Flattened coverage vector containing data for all components.
     * @param bp_component_offsets Array of index boundaries delineating components inside flat_bp_coverage.
     * @param component_names Array of component identifiers matching the sequence boundaries.
     * @param config Plot layout and configuration settings.
     *
     * @throws std::invalid_argument If offset array size or boundary ranges are inconsistent.
     */
    void writeLinearPlotGlobal(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &flat_bp_coverage,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names,
        const PlotConfig &config
    );

    /**
     * @brief Generate a multi-panel circular coverage plot for all components.
     *
     * Splits a flattened coverage table into component-specific coverage ranges,
     * prepares a CircularPlotPackage for each component, serializes all packages
     * into a single binary request, and invokes the Python rendering backend to
     * produce a grid of circular coverage plots.
     *
     * Numeric preprocessing for individual components is performed in parallel,
     * while request serialization and rendering are performed sequentially.
     * This minimizes rendering overhead by generating the complete figure with
     * a single Python process.
     *
     * @param output_png Output PNG image path.
     * @param flat_bp_coverage Flattened per-base coverage table containing all components.
     * @param bp_component_offsets Component boundaries within the flattened coverage table.
     * @param component_names Names associated with each component.
     * @param config Plot rendering configuration.
     *
     * @throws std::invalid_argument If component offsets are inconsistent with the coverage table.
     * @throws std::runtime_error If request generation or rendering fails.
     */
    void writeCircularPlotGlobal(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &flat_bp_coverage,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names,
        const PlotConfig &config
    );
} // namespace output

#endif // CDX_COVERAGE_OUTPUT_PLOT_H
