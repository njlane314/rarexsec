#pragma once

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TStyle.h"

#include <ROOT/RDataFrame.hxx>

namespace rarexsec {
namespace plot {

class Plotter;

class EventDisplay {
    friend class Plotter;

  public:
    enum class Mode { Detector,
                      Semantic };

    static Mode parse_mode(const std::string& s) {
        if (s == "semantic" || s == "Semantic")
            return Mode::Semantic;
        return Mode::Detector;
    }

    struct Spec {
        std::string id;
        std::string title;
        Mode mode{Mode::Detector};
        int grid_w{0};
        int grid_h{0};
    };

    struct Options {
        std::string out_dir = "plots";
        int canvas_size = 1400;
        double margin = 0.10;
        bool use_log_z = true;

        double det_threshold = 4.0;
        double det_min = 1.0;
        double det_max = 1000.0;

        bool show_legend = true;
        int legend_cols = 5;
    };

    using DetectorData = std::vector<float>;
    using SemanticData = std::vector<int>;

    void draw(TCanvas& canvas);

    void draw_and_save(const std::string& image_format = "png");
    void draw_and_save(const std::string& image_format, const std::string& file_override);

    struct BatchOptions {
        std::string selection_expr;
        unsigned long long n_events{1};

        std::string out_dir{"./plots/event_displays"};
        std::string image_format{"png"};
        std::string combined_pdf;
        std::string manifest_path;

        std::vector<std::string> planes{"U", "V", "W"};

        struct Columns {
            std::string run = "run";
            std::string sub = "sub";
            std::string evt = "evt";
            std::string det_u = "event_detector_image_u";
            std::string det_v = "event_detector_image_v";
            std::string det_w = "event_detector_image_w";
            std::string sem_u = "semantic_image_u";
            std::string sem_v = "semantic_image_v";
            std::string sem_w = "semantic_image_w";
        } cols;

        std::string file_pattern{"{plane}_{run}_{sub}_{evt}"};

        Mode mode{Mode::Detector};
        Options display;
    };

    static void render_from_rdf(ROOT::RDF::RNode df, const BatchOptions& opt);

  private:
    EventDisplay(Spec spec, Options opt, DetectorData data);
    EventDisplay(Spec spec, Options opt, SemanticData data);

    void setup_canvas(TCanvas& c) const;
    void build_histogram();

    void draw_detector(TCanvas& c);
    void draw_semantic(TCanvas& c);
    void draw_semantic_legend();

    static std::pair<int, int> deduce_grid(int requested_w,
                                           int requested_h,
                                           std::size_t flat_size);

  private:
    Spec spec_;
    Options opt_;

    std::variant<DetectorData, SemanticData> data_;

    std::unique_ptr<TH2F> hist_;
    std::unique_ptr<TLegend> legend_;
    std::vector<std::unique_ptr<TH1F>> legend_entries_;

    std::string plot_name_;
    std::string output_directory_;
};

}
}
