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
     * @brief Paquet de données prêt-à-dessiner pour un graphe de couverture
     *        linéaire (une composante ou une région).
     *
     * @note Malgré son nom historique (hérité de l'ancien format d'échange
     *       binaire avec le script Python), cette structure ne sert plus
     *       qu'en mémoire : le renderer Cairo la consomme directement.
     *       Anciennement dans bin_plot.h (fichier retiré : cette structure
     *       en était l'unique contenu, déplacée ici avec les autres types
     *       du module de rendu).
     */
    struct LinearPlotPackageBin {
        std::string component_name;
        std::size_t query_start{0};
        std::size_t query_end{0};
        bool logarithmic{false};
        int log_base{0};
        double y_upper_limit{1.0};

        // Grille régulière pour les coordonnées génomiques (X)
        double x_start{0.0};
        double x_step{1.0};

        std::vector<double> y;
        std::vector<double> tick_positions;
        std::vector<std::string> tick_labels;
    };
} // namespace cdx

namespace output {
    enum class Topology {
        Linear,
        Circular
    };

    struct PlotConfig {
        double smoothing = 0.01;
        std::size_t max_plot_points = 10000;
        int dpi = 300;
        std::optional<int> log_base;

        double figure_width = 7.0;
        double figure_height = 4.5;

        std::string line_color = "#1E3A8A";
        std::string fill_color = "#93C5FD";
    };

    struct PlotData {
        std::vector<std::size_t> x;
        std::vector<double> y;
        std::size_t window_size = 1;
    };

    struct GridLayout {
        std::size_t rows = 1;
        std::size_t columns = 1;
    };

    struct CoverageTicks {
        std::vector<double> values;
        double upper_limit = 1.0;
    };

    struct LogCoverageTicks {
        std::vector<double> raw_values;
        std::vector<double> display_values;
        double display_upper_limit = 1.0;
    };

    /*
     * Fonctions numériques communes.
     */

    [[nodiscard]]
    GridLayout chooseGlobalGraphGrid(
        std::size_t component_count
    );

    [[nodiscard]]
    PlotData prepareCoverageForPlot(
        const std::vector<cdx::Coverage> &coverage,
        double smoothing,
        std::size_t max_plot_points,
        Topology topology
    );

    /**
     * @brief Overload accepting pre-transformed double values (may contain
     *        NaN, treated as visual gaps: propagates the same way through
     *        the cumulative-sum smoothing here as it does through NumPy's
     *        cumsum in the Python reference - not "fixed", intentionally
     *        kept consistent with the established prototype behavior).
     *        Used by the circular backend, where sentinel/out-of-query
     *        positions are masked to NaN rather than replaced by 0 (unlike
     *        the linear backend, cf. prepareLinearPlotPackage).
     */
    [[nodiscard]]
    PlotData prepareCoverageForPlot(
        const std::vector<double> &coverage,
        double smoothing,
        std::size_t max_plot_points,
        Topology topology
    );

    [[nodiscard]]
    CoverageTicks calculateCoverageTicks(
        double maximum_coverage,
        std::size_t target_tick_count = 4
    );

    [[nodiscard]]
    std::vector<double> logTransformCoverage(
        const std::vector<double> &values,
        int log_base
    );

    [[nodiscard]]
    LogCoverageTicks logCoverageTicks(
        double maximum_coverage,
        int log_base
    );

    [[nodiscard]]
    cdx::LinearPlotPackageBin prepareLinearPlotPackage(
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name,
        std::size_t offset,
        const PlotConfig &config
    );

    /*
     * Backend linéaire et circulaire.
     */

    void writeLinearPlotQuery(
        const std::filesystem::path& output_png,
        const std::vector<cdx::Coverage>& coverage,
        const std::string& component_name,
        std::size_t offset,
        const PlotConfig& config
    );

    /**
     * @brief Rend un graphe de couverture circulaire pour une composante (ou une
     *        plage de celle-ci) via un sous-processus Python (pycirclize).
     *
     * @param component_coverage Couverture sur toute la longueur de la composante
     *        (les positions hors requête sont masquées avec des sentinelles, cf.
     *        trimCoverageToQuery) : contrairement au backend linéaire, on ne
     *        découpe PAS ce tableau avant l'appel, car le lissage circulaire et
     *        les requêtes traversant l'origine ont besoin du contexte complet.
     * @param compo_length Longueur totale de la composante en bp.
     * @param query_bound Bornes inclusives (query_start, query_end) de la requête.
     *        Vaut (0, compo_length - 1) pour représenter la composante entière.
     *        query_start > query_end signale une requête traversant l'origine.
     */
    void writeCircularPlotQuery(
        const std::filesystem::path& output_png,
        const std::vector<cdx::Coverage>& component_coverage,
        const std::string& component_name,
        cdx::PosBp compo_length,
        std::pair<cdx::PosBp, cdx::PosBp> query_bound,
        const PlotConfig& config
    );

    void writeLinearPlotGlobal(
        const std::filesystem::path& output_png,
        const std::vector<cdx::Coverage>& flat_bp_coverage,
        const std::vector<cdx::PosBp>& bp_component_offsets,
        const std::vector<std::string>& component_names,
        const PlotConfig& config
    );

    /**
     * @brief Rend un graphe de couverture circulaire multi-composantes (une
     *        piste par composante, chacune son propre idéogramme) via un
     *        sous-processus Python (pycirclize) : une seule invocation
     *        Python pour la grille entière (cf. python_script/circular_plot.py:
     *        _render_global), pas un appel par composante.
     */
    void writeCircularPlotGlobal(
        const std::filesystem::path& output_png,
        const std::vector<cdx::Coverage>& flat_bp_coverage,
        const std::vector<cdx::PosBp>& bp_component_offsets,
        const std::vector<std::string>& component_names,
        const PlotConfig& config
    );

} // namespace output

#endif // CDX_COVERAGE_OUTPUT_PLOT_H