/**
 * @file cairo_plot.cpp
 * @brief In-process linear coverage plot rendering via libcairo (no subprocess).
 *
 * Consomme les paquets produits par prepareLinearPlotPackage() (linear_plot.cpp)
 * et les dessine directement dans un cairo_surface_t rasterisé en PNG. Aucun
 * processus externe n'est lancé : tout reste dans le même processus C++, ce
 * qui élimine le coût de démarrage d'un interpréteur (Python ou gnuplot) et
 * garde le pipeline sous la seconde.
 */

#include "output_plot.h"
#include "config.h"

#include <cairo.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace output {
    namespace {
        struct Rgb {
            double r = 0.0;
            double g = 0.0;
            double b = 0.0;
        };

        struct Rect {
            double x = 0.0;
            double y = 0.0;
            double w = 0.0;
            double h = 0.0;
        };

        /**
         * @brief Description d'un panneau prêt à dessiner (un graphe = un panneau).
         */
        struct PanelSpec {
            std::string title;
            std::string x_label;
            std::string y_label;

            double x_min = 0.0;
            double x_max = 1.0;
            double y_min = 0.0;
            double y_max = 1.0;

            std::vector<std::pair<double, double> > points; // (x, y) en espace de données
            std::vector<double> x_ticks; // en espace de données X
            std::vector<double> y_ticks; // en espace d'affichage Y (déjà log-transformé si applicable)
            std::vector<std::string> y_tick_labels; // parallèle à y_ticks; vide => formatage automatique

            Rgb line_color{0.117, 0.227, 0.541};
            Rgb fill_color{0.576, 0.773, 0.992};
        };

        Rgb parseHexColor(const std::string &hex) {
            auto hex_digit = [](const char c) -> int {
                if (c >= '0' && c <= '9') return c - '0';
                if (c >= 'a' && c <= 'f') return c - 'a' + 10;
                if (c >= 'A' && c <= 'F') return c - 'A' + 10;
                return 0;
            };

            if (hex.size() < 7 || hex[0] != '#') {
                return {0.117, 0.227, 0.541};
            }

            const int r = hex_digit(hex[1]) * 16 + hex_digit(hex[2]);
            const int g = hex_digit(hex[3]) * 16 + hex_digit(hex[4]);
            const int b = hex_digit(hex[5]) * 16 + hex_digit(hex[6]);

            return {
                static_cast<double>(r) / 255.0,
                static_cast<double>(g) / 255.0,
                static_cast<double>(b) / 255.0
            };
        }

        /**
         * @brief Calcule des positions de graduations "rondes" sur [data_min, data_max].
         *
         * Même logique que output::calculateCoverageTicks (pas 1/2/2.5/3/4/5/10 x 10^n),
         * mais bornée à l'intervalle demandé plutôt que de partir de zéro : utilisable
         * pour un axe X quelconque, à l'image de matplotlib.ticker.MaxNLocator.
         */
        std::vector<double> niceAxisTicks(
            const double data_min,
            const double data_max,
            const std::size_t target_count
        ) {
            if (!(data_max > data_min) || target_count == 0) {
                return {data_min};
            }

            const double raw_step = (data_max - data_min) / static_cast<double>(target_count);
            const double magnitude = std::pow(10.0, std::floor(std::log10(raw_step)));
            const double normalized = raw_step / magnitude;

            static constexpr std::array<double, 7> multipliers{1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0};
            double multiplier = 10.0;
            for (const double m: multipliers) {
                if (m >= normalized) {
                    multiplier = m;
                    break;
                }
            }

            const double step = multiplier * magnitude;
            if (!(step > 0.0)) {
                return {data_min};
            }

            const double start = std::ceil(data_min / step) * step;

            std::vector<double> ticks;
            for (double v = start; v <= data_max + step * 1e-9; v += step) {
                ticks.push_back(v);
            }
            if (ticks.empty()) {
                ticks.push_back(data_min);
            }
            return ticks;
        }

        std::string formatPlainInteger(const double value) {
            return std::to_string(static_cast<long long>(std::llround(value)));
        }

        void setDashedStroke(cairo_t *cr, const double px_per_pt) {
            const std::array<double, 2> dashes{4.0 * px_per_pt, 3.0 * px_per_pt};
            cairo_set_dash(cr, dashes.data(), static_cast<int>(dashes.size()), 0.0);
        }

        void clearDash(cairo_t *cr) {
            cairo_set_dash(cr, nullptr, 0, 0.0);
        }

        /**
         * @brief Dessine un panneau complet (titre, grille, courbe, axes, ticks) dans `panel`.
         */
        void drawPanel(cairo_t *cr, const Rect &panel, const PanelSpec &spec, const double dpi) {
            const double px_per_pt = dpi / 72.0;

            const double title_size = 15.0 * px_per_pt;
            const double label_size = 12.0 * px_per_pt;
            const double tick_size = 10.5 * px_per_pt;

            // ---- Fond du panneau -------------------------------------------------
            cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
            cairo_rectangle(cr, panel.x, panel.y, panel.w, panel.h);
            cairo_fill(cr);

            const double margin_top = title_size * 2.3;
            const double margin_left = label_size * 4.8;
            const double margin_bottom = label_size * 3.4;
            const double margin_right = label_size * 1.1;

            const Rect axes{
                panel.x + margin_left,
                panel.y + margin_top,
                std::max(1.0, panel.w - margin_left - margin_right),
                std::max(1.0, panel.h - margin_top - margin_bottom)
            };

            const double x_range = (spec.x_max > spec.x_min) ? (spec.x_max - spec.x_min) : 1.0;
            const double y_range = (spec.y_max > spec.y_min) ? (spec.y_max - spec.y_min) : 1.0;

            auto mapX = [&](const double x) {
                return axes.x + (x - spec.x_min) / x_range * axes.w;
            };
            auto mapY = [&](const double y) {
                return axes.y + axes.h - (y - spec.y_min) / y_range * axes.h;
            };

            // ---- Titre (gras, aligné à gauche) ------------------------------------
            cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
            cairo_set_font_size(cr, title_size);
            cairo_set_source_rgb(cr, 0.06, 0.09, 0.16);
            cairo_move_to(cr, axes.x, panel.y + margin_top * 0.58);
            cairo_show_text(cr, spec.title.c_str());

            // ---- Grille pointillée (X et Y) ----------------------------------------
            cairo_save(cr);
            cairo_rectangle(cr, axes.x, axes.y, axes.w, axes.h);
            cairo_clip(cr);

            cairo_set_line_width(cr, 0.8 * px_per_pt);
            cairo_set_source_rgba(cr, 0.71, 0.76, 0.83, 0.55);
            setDashedStroke(cr, px_per_pt);

            for (const double tick: spec.y_ticks) {
                const double y_px = mapY(tick);
                cairo_move_to(cr, axes.x, y_px);
                cairo_line_to(cr, axes.x + axes.w, y_px);
                cairo_stroke(cr);
            }
            for (const double tick: spec.x_ticks) {
                const double x_px = mapX(tick);
                cairo_move_to(cr, x_px, axes.y);
                cairo_line_to(cr, x_px, axes.y + axes.h);
                cairo_stroke(cr);
            }
            clearDash(cr);

            // ---- Aire remplie + courbe de couverture --------------------------------
            if (!spec.points.empty()) {
                const double baseline_px = mapY(spec.y_min);

                cairo_move_to(cr, mapX(spec.points.front().first), baseline_px);
                for (const auto &[x, y]: spec.points) {
                    cairo_line_to(cr, mapX(x), mapY(y));
                }
                cairo_line_to(cr, mapX(spec.points.back().first), baseline_px);
                cairo_close_path(cr);

                cairo_set_source_rgba(cr, spec.fill_color.r, spec.fill_color.g, spec.fill_color.b, 0.18);
                cairo_fill(cr);

                cairo_move_to(cr, mapX(spec.points.front().first), mapY(spec.points.front().second));
                for (const auto &[x, y]: spec.points) {
                    cairo_line_to(cr, mapX(x), mapY(y));
                }
                cairo_set_source_rgb(cr, spec.line_color.r, spec.line_color.g, spec.line_color.b);
                cairo_set_line_width(cr, 1.3 * px_per_pt);
                cairo_stroke(cr);
            }

            cairo_restore(cr); // fin du clip sur la zone de tracé

            // ---- Bordures (spines gauche/bas seulement) -----------------------------
            cairo_set_source_rgb(cr, 0.58, 0.64, 0.72);
            cairo_set_line_width(cr, 1.0 * px_per_pt);
            cairo_move_to(cr, axes.x, axes.y);
            cairo_line_to(cr, axes.x, axes.y + axes.h);
            cairo_line_to(cr, axes.x + axes.w, axes.y + axes.h);
            cairo_stroke(cr);

            // ---- Graduations + labels Y ----------------------------------------------
            cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
            cairo_set_font_size(cr, tick_size);
            cairo_set_source_rgb(cr, 0.2, 0.25, 0.31);

            for (std::size_t i = 0; i < spec.y_ticks.size(); ++i) {
                const double y_px = mapY(spec.y_ticks[i]);

                cairo_move_to(cr, axes.x - 4.0 * px_per_pt, y_px);
                cairo_line_to(cr, axes.x, y_px);
                cairo_set_line_width(cr, 1.0 * px_per_pt);
                cairo_stroke(cr);

                const std::string label = (i < spec.y_tick_labels.size())
                                               ? spec.y_tick_labels[i]
                                               : formatPlainInteger(spec.y_ticks[i]);

                cairo_text_extents_t extents;
                cairo_text_extents(cr, label.c_str(), &extents);

                cairo_move_to(
                    cr,
                    axes.x - 8.0 * px_per_pt - extents.width - extents.x_bearing,
                    y_px - extents.y_bearing - extents.height / 2.0
                );
                cairo_show_text(cr, label.c_str());
            }

            // ---- Graduations + labels X -----------------------------------------------
            for (const double tick: spec.x_ticks) {
                const double x_px = mapX(tick);

                cairo_move_to(cr, x_px, axes.y + axes.h);
                cairo_line_to(cr, x_px, axes.y + axes.h + 4.0 * px_per_pt);
                cairo_stroke(cr);

                const std::string label = formatPlainInteger(tick);
                cairo_text_extents_t extents;
                cairo_text_extents(cr, label.c_str(), &extents);

                cairo_move_to(
                    cr,
                    x_px - extents.width / 2.0 - extents.x_bearing,
                    axes.y + axes.h + 8.0 * px_per_pt - extents.y_bearing
                );
                cairo_show_text(cr, label.c_str());
            }

            // ---- Label d'axe X ----------------------------------------------------------
            cairo_set_font_size(cr, label_size);
            {
                cairo_text_extents_t extents;
                cairo_text_extents(cr, spec.x_label.c_str(), &extents);
                cairo_move_to(
                    cr,
                    axes.x + axes.w / 2.0 - extents.width / 2.0 - extents.x_bearing,
                    panel.y + panel.h - label_size * 0.7
                );
                cairo_show_text(cr, spec.x_label.c_str());
            }

            // ---- Label d'axe Y (pivoté à 90°) --------------------------------------------
            {
                cairo_text_extents_t extents;
                cairo_text_extents(cr, spec.y_label.c_str(), &extents);

                cairo_save(cr);
                cairo_translate(cr, panel.x + label_size * 1.1, axes.y + axes.h / 2.0);
                cairo_rotate(cr, -M_PI / 2.0);
                cairo_move_to(cr, -extents.width / 2.0 - extents.x_bearing, -extents.y_bearing - extents.height / 2.0);
                cairo_show_text(cr, spec.y_label.c_str());
                cairo_restore(cr);
            }
        }

        /**
         * @brief Construit un PanelSpec à partir d'un paquet de données déjà préparé.
         */
        PanelSpec buildPanelSpec(
            const cdx::LinearPlotPackageBin &package,
            const bool query_mode,
            const PlotConfig &config
        ) {
            PanelSpec spec;

            if (query_mode) {
                std::ostringstream title;
                title << "Component " << package.component_name
                        << " (Region " << cfg::formatInteger(package.query_start)
                        << " - " << cfg::formatInteger(package.query_end) << ")";
                spec.title = title.str();
            } else {
                spec.title = "Component " + package.component_name;
            }

            spec.x_label = "Genomic Coordinate (bp)";
            spec.y_label = package.logarithmic
                               ? ("Depth of Coverage (log base " + std::to_string(package.log_base) + ")")
                               : "Depth of Coverage";

            spec.x_min = static_cast<double>(package.query_start);
            spec.x_max = (package.query_end > package.query_start)
                             ? static_cast<double>(package.query_end)
                             : static_cast<double>(package.query_start) + 1.0;

            spec.y_min = 0.0;
            spec.y_max = std::max(package.y_upper_limit, 1.0);

            spec.points.reserve(package.y.size());
            for (std::size_t i = 0; i < package.y.size(); ++i) {
                const double x = package.x_start + static_cast<double>(i) * package.x_step;
                spec.points.emplace_back(x, package.y[i]);
            }

            spec.x_ticks = niceAxisTicks(spec.x_min, spec.x_max, 8);

            if (package.logarithmic) {
                spec.y_ticks = package.tick_positions;
                spec.y_tick_labels = package.tick_labels;
            } else {
                spec.y_ticks = niceAxisTicks(spec.y_min, spec.y_max, 6);
                spec.y_tick_labels.reserve(spec.y_ticks.size());
                for (const double tick: spec.y_ticks) {
                    spec.y_tick_labels.push_back(cfg::formatInteger(static_cast<long long>(std::llround(tick))) + "x");
                }
            }

            spec.line_color = parseHexColor(config.line_color);
            spec.fill_color = parseHexColor(config.fill_color);

            return spec;
        }

        void writePngOrThrow(cairo_surface_t *surface, const std::filesystem::path &output_png) {
            cairo_surface_flush(surface);
            const cairo_status_t status = cairo_surface_write_to_png(surface, output_png.string().c_str());

            if (status != CAIRO_STATUS_SUCCESS) {
                throw std::runtime_error(
                    "Failed to write coverage plot PNG (" + std::string(cairo_status_to_string(status)) +
                    "): " + output_png.string()
                );
            }
        }
    } // namespace

    void writeLinearPlotQuery(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &coverage,
        const std::string &component_name,
        std::size_t offset,
        const PlotConfig &config
    ) {
        const cdx::LinearPlotPackageBin package = prepareLinearPlotPackage(coverage, component_name, offset, config);
        const PanelSpec spec = buildPanelSpec(package, /*query_mode=*/true, config);

        const auto parent_dir = output_png.parent_path();
        if (!parent_dir.empty()) {
            std::filesystem::create_directories(parent_dir);
        }

        const auto width_px = std::max(1, static_cast<int>(std::llround(config.figure_width * config.dpi)));
        const auto height_px = std::max(1, static_cast<int>(std::llround(config.figure_height * config.dpi)));

        cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width_px, height_px);
        cairo_t *cr = cairo_create(surface);

        drawPanel(
            cr,
            Rect{0.0, 0.0, static_cast<double>(width_px), static_cast<double>(height_px)},
            spec,
            static_cast<double>(config.dpi)
        );

        try {
            writePngOrThrow(surface, output_png);
        } catch (...) {
            cairo_destroy(cr);
            cairo_surface_destroy(surface);
            throw;
        }

        cairo_destroy(cr);
        cairo_surface_destroy(surface);
    }

    void writeLinearPlotGlobal(
        const std::filesystem::path &output_png,
        const std::vector<cdx::Coverage> &flat_bp_coverage,
        const std::vector<cdx::PosBp> &bp_component_offsets,
        const std::vector<std::string> &component_names,
        const PlotConfig &config
    ) {
        if (bp_component_offsets.size() < 2) {
            throw std::invalid_argument("bp_component_offsets must contain at least two offsets.");
        }

        const std::size_t component_count = bp_component_offsets.size() - 1;
        const GridLayout grid = chooseGlobalGraphGrid(component_count);

        const auto parent_dir = output_png.parent_path();
        if (!parent_dir.empty()) {
            std::filesystem::create_directories(parent_dir);
        }

        const auto panel_w = std::max(1, static_cast<int>(std::llround(config.figure_width * config.dpi)));
        const auto panel_h = std::max(1, static_cast<int>(std::llround(config.figure_height * config.dpi)));

        const int width_px = panel_w * static_cast<int>(grid.columns);
        const int height_px = panel_h * static_cast<int>(grid.rows);

        /*
         * Chaque panneau (préparation des données + tracé Cairo) est
         * indépendant des autres : on les rend en parallèle, chacun dans sa
         * propre cairo_surface_t privée (un cairo_t ne se partage pas entre
         * threads). La composition finale dans le canevas commun reste
         * séquentielle mais ne coûte qu'un simple blit par panneau.
         */
        std::vector<cairo_surface_t *> panel_surfaces(component_count, nullptr);

        std::exception_ptr worker_exception;
        std::mutex exception_mutex;

        auto render_component = [&](const std::size_t component_id) {
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

                const std::string name = (component_id < component_names.size())
                                              ? component_names[component_id]
                                              : std::to_string(component_id);

                const cdx::LinearPlotPackageBin package = prepareLinearPlotPackage(
                    component_coverage, name, 0, config
                );
                const PanelSpec spec = buildPanelSpec(package, /*query_mode=*/false, config);

                cairo_surface_t *panel_surface = cairo_image_surface_create(
                    CAIRO_FORMAT_ARGB32, panel_w, panel_h
                );
                cairo_t *panel_cr = cairo_create(panel_surface);

                drawPanel(
                    panel_cr,
                    Rect{0.0, 0.0, static_cast<double>(panel_w), static_cast<double>(panel_h)},
                    spec,
                    static_cast<double>(config.dpi)
                );

                cairo_destroy(panel_cr);
                panel_surfaces[component_id] = panel_surface;
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
                        render_component(idx);
                    }
                });
            }

            for (auto &worker: workers) {
                worker.join();
            }
        }

        if (worker_exception) {
            for (cairo_surface_t *panel_surface: panel_surfaces) {
                if (panel_surface != nullptr) {
                    cairo_surface_destroy(panel_surface);
                }
            }
            std::rethrow_exception(worker_exception);
        }

        cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width_px, height_px);
        cairo_t *cr = cairo_create(surface);

        cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
        cairo_paint(cr);

        try {
            for (std::size_t component_id = 0; component_id < component_count; ++component_id) {
                const std::size_t row = component_id / grid.columns;
                const std::size_t col = component_id % grid.columns;

                cairo_set_source_surface(
                    cr,
                    panel_surfaces[component_id],
                    static_cast<double>(col) * panel_w,
                    static_cast<double>(row) * panel_h
                );
                cairo_paint(cr);
                cairo_surface_destroy(panel_surfaces[component_id]);
            }

            writePngOrThrow(surface, output_png);
        } catch (...) {
            cairo_destroy(cr);
            cairo_surface_destroy(surface);
            throw;
        }

        cairo_destroy(cr);
        cairo_surface_destroy(surface);
    }
} // namespace output
