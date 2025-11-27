#include "rarexsec/Hub.h"
#include "rarexsec/Processor.h"
#include "rarexsec/plot/EventDisplay.h"

#include <ROOT/RDataFrame.hxx>
#include <TError.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static std::vector<std::string> split_periods(const char* env_val)
{
    std::vector<std::string> out;
    if (!env_val)
        return out;

    std::string s(env_val);
    std::istringstream iss(s);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        if (!tok.empty())
            out.push_back(tok);
    }
    return out;
}

static ROOT::RDF::RNode apply_mc_slice(ROOT::RDF::RNode node,
                                       const rarexsec::Entry& rec)
{
    using rarexsec::Slice;
    using rarexsec::Source;

    if (rec.source == Source::MC) {
        if (rec.slice == Slice::StrangenessInclusive)
            return node.Filter([](bool s) { return s; }, {"is_strange"});
        if (rec.slice == Slice::BeamInclusive)
            return node.Filter([](bool s) { return !s; }, {"is_strange"});
    } else if (rec.slice != Slice::None) {
        std::cerr << "[event_display] WARNING: slice requested for non-MC sample at "
                  << rec.beamline << "/" << rec.period << std::endl;
    }
    return node;
}

static rarexsec::Entry const* pick_entry(rarexsec::Hub& hub,
                                         const std::string& beamline,
                                         const std::vector<std::string>& periods,
                                         bool want_mc)
{
    if (want_mc) {
        auto sim_entries = hub.simulation_entries(beamline, periods);
        if (!sim_entries.empty())
            return sim_entries.front();
        std::cerr << "[event_display] No simulation entries found, falling back to data.\n";
    }

    auto data_entries = hub.data_entries(beamline, periods);
    if (!data_entries.empty())
        return data_entries.front();

    if (!want_mc) {
        auto sim_entries = hub.simulation_entries(beamline, periods);
        if (!sim_entries.empty()) {
            std::cerr << "[event_display] No data entries found, using MC instead.\n";
            return sim_entries.front();
        }
    }

    std::cerr << "[event_display] No suitable entries found for beamline=" << beamline << "\n";
    return nullptr;
}

static void run_event_display(bool use_semantic)
{
    gErrorIgnoreLevel = kWarning;

    const char* cfg_c     = std::getenv("RAREXSEC_CFG");
    const char* beam_c    = std::getenv("RAREXSEC_BEAMLINE");
    const char* periods_c = std::getenv("RAREXSEC_PERIODS");

    const std::string cfg      = cfg_c  ? cfg_c  : "data/samples.json";
    const std::string beamline = beam_c ? beam_c : "numi-fhc";

    std::vector<std::string> periods = split_periods(periods_c);
    if (periods.empty())
        periods.push_back("run1");

    std::cout << "[event_display] cfg=" << cfg
              << " beamline=" << beamline
              << " periods=";
    for (auto& p : periods) std::cout << " " << p;
    std::cout << std::endl;

    rarexsec::Hub hub(cfg);

    const bool want_mc = use_semantic;
    const rarexsec::Entry* rec_ptr = pick_entry(hub, beamline, periods, want_mc);
    if (!rec_ptr) {
        std::cerr << "[event_display] No entry to process, exiting.\n";
        return;
    }
    const rarexsec::Entry& rec = *rec_ptr;

    std::cout << "[event_display] Using "
              << (rec.source == rarexsec::Source::MC ? "MC" :
                  rec.source == rarexsec::Source::Data ? "DATA" : "EXT")
              << " sample, file=" << rec.file << std::endl;

    static const std::string tree_name = "nuselection/EventSelectionFilter";

    ROOT::RDataFrame df(tree_name, rec.files);

    ROOT::RDF::RNode node = rarexsec::processor().run(df, rec);
    node = apply_mc_slice(node, rec);

    using rarexsec::plot::EventDisplay;

    EventDisplay::BatchOptions opt;

    opt.out_dir       = use_semantic ? "plots/event_display_semantic"
                                     : "plots/event_display_detector";
    opt.image_format  = "png";
    opt.combined_pdf  = "";
    opt.manifest_path = opt.out_dir + "/manifest.json";

    // How many events to render per sample. Default is 10, but allow an
    // override via the environment for quick testing:
    //
    //   export RAREXSEC_N_EVENTS=25
    //
    int n_events = 10;
    if (const char* nev_c = std::getenv("RAREXSEC_N_EVENTS")) {
        try {
            n_events = std::max(1, std::stoi(nev_c));
        } catch (...) {
            std::cerr << "[event_display] Invalid RAREXSEC_N_EVENTS=" << nev_c
                      << ", using default n_events=" << n_events << '\n';
        }
    }
    opt.n_events = n_events;

    std::cout << "[event_display] Will render up to " << opt.n_events
              << " events per sample" << std::endl;

    std::string base_sel = "";
    if (rec.source == rarexsec::Source::MC && use_semantic) {
    }
    opt.selection_expr = base_sel;

    opt.planes = {"U", "V", "W"};

    opt.file_pattern = use_semantic
        ? "evd_sem_{plane}_run{run}_sub{sub}_evt{evt}"
        : "evd_det_{plane}_run{run}_sub{sub}_evt{evt}";

    opt.cols.run = "run";
    opt.cols.sub = "sub";
    opt.cols.evt = "evt";

    opt.cols.det_u = "detector_image_u";
    opt.cols.det_v = "detector_image_v";
    opt.cols.det_w = "detector_image_w";

    opt.cols.sem_u = "semantic_image_u";
    opt.cols.sem_v = "semantic_image_v";
    opt.cols.sem_w = "semantic_image_w";

    opt.mode = use_semantic ? EventDisplay::Mode::Semantic
                            : EventDisplay::Mode::Detector;

    opt.display.canvas_size   = 900;
    opt.display.margin        = 0.10;

    // For the detector view we now plot the raw ADC values and derive the
    // z-range per image inside EventDisplay, so these are just defaults.
    opt.display.det_threshold = 0.0f;   // raw filling, only used for log safety
    opt.display.det_min       = 0.0f;   // will be set per image
    opt.display.det_max       = 0.0f;   // will be set per image
    opt.display.use_log_z     = true;   // enable log-scale detector display

    opt.display.show_legend   = true;
    opt.display.legend_cols   = 4;

    std::cout << "[event_display] Rendering "
              << (use_semantic ? "semantic" : "detector")
              << " event displays..." << std::endl;

    EventDisplay::render_from_rdf(node, opt);

    std::cout << "[event_display] Done. Check " << opt.out_dir << std::endl;
}

void event_display_detector()
{
    run_event_display(false);
}

void event_display_semantic()
{
    run_event_display(true);
}

void event_display()
{
    event_display_detector();
}
