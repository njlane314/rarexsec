#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <cctype>
#include <exception>
#include <iostream>
#include <vector>

#include <rarexsec/Hub.h>
#include <rarexsec/proc/Env.h>
#include <rarexsec/plot/Channels.h>
#include <rarexsec/plot/Descriptors.h>
#include <rarexsec/plot/Plotter.h>

static std::vector<std::string> get_beamlines(const rarexsec::Env& env) {
    std::vector<std::string> out;

    if (const char* bls = gSystem->Getenv("RAREXSEC_BEAMLINES")) {
        std::string s{bls};
        std::string tok;

        auto flush = [&]() {
            if (!tok.empty()) {
                out.push_back(tok);
                tok.clear();
            }
        };

        for (char ch : s) {
            if (ch == ',' || std::isspace(static_cast<unsigned char>(ch))) {
                flush();
            } else {
                tok.push_back(ch);
            }
        }
        flush();
    }

    if (out.empty()) {
        out.push_back(env.beamline);
    }

    return out;
}

void plot_topology_variables() {
    try {
        ROOT::EnableImplicitMT();
        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        auto beamlines = get_beamlines(env);
        using SamplePtr = decltype(hub.simulation_entries(env.beamline, env.periods))::value_type;
        std::vector<SamplePtr> mc_samples;
        for (const auto& bl : beamlines) {
            auto sub = hub.simulation_entries(bl, env.periods);
            mc_samples.insert(mc_samples.end(), sub.begin(), sub.end());
        }

        std::cout << "Loaded beamlines";
        for (const auto& bl : beamlines) {
            std::cout << ' ' << bl;
        }
        std::cout << " for";
        for (const auto& period : env.periods) {
            std::cout << ' ' << period;
        }
        std::cout << " with " << mc_samples.size() << " simulation samples." << std::endl;

        rarexsec::plot::Options opt;
        opt.out_dir = "plots/selection";
        opt.use_log_y = true;
        opt.overlay_signal = true;
        opt.annotate_numbers = true;
        opt.image_format = "pdf";
        opt.legend_on_top = true;
        opt.beamline = env.beamline;
        opt.periods = env.periods;
        opt.analysis_region_label = "Empty Selection";
        opt.signal_channels = rarexsec::plot::Channels::signal_keys();

        rarexsec::plot::Plotter plotter(opt);

        rarexsec::plot::TH1DModel beam_pe;
        beam_pe.id = "optical_filter_pe_beam";
        beam_pe.title = ";Beamline PMT PE;Events";
        beam_pe.nbins = 100;
        beam_pe.xmin = 0.0;
        beam_pe.xmax = 200.0;
        beam_pe.sel = rarexsec::selection::Preset::Empty;

        rarexsec::plot::TH1DModel veto_pe = beam_pe;
        veto_pe.id = "optical_filter_pe_veto";
        veto_pe.title = ";Veto PMT PE;Events";
        veto_pe.xmax = 100.0;

        rarexsec::plot::TH1DModel software_trigger = beam_pe;
        software_trigger.id = "software_trigger";
        software_trigger.title = ";Software Trigger Decision;Events";
        software_trigger.nbins = 3;
        software_trigger.xmin = -0.5;
        software_trigger.xmax = 2.5;

        rarexsec::plot::TH1DModel num_slices = beam_pe;
        num_slices.id = "num_slices";
        num_slices.title = ";Number of Slices;Events";
        num_slices.nbins = 10;
        num_slices.xmin = -0.5;
        num_slices.xmax = 9.5;

        rarexsec::plot::TH1DModel topology_score = beam_pe;
        topology_score.id = "topological_score";
        topology_score.title = ";Topological Score;Events";
        topology_score.nbins = 100;
        topology_score.xmin = 0.0;
        topology_score.xmax = 1.0;

        rarexsec::plot::TH1DModel fiducial = beam_pe;
        fiducial.id = "in_reco_fiducial";
        fiducial.title = ";In Reconstructed Fiducial Volume;Events";
        fiducial.nbins = 2;
        fiducial.xmin = -0.5;
        fiducial.xmax = 1.5;

        rarexsec::plot::TH1DModel contained = topology_score;
        contained.id = "contained_fraction";
        contained.title = ";Contained Fraction;Events";

        rarexsec::plot::TH1DModel cluster = contained;
        cluster.id = "slice_cluster_fraction";
        cluster.title = ";Slice Cluster Fraction;Events";

        rarexsec::plot::TH1DModel muon_track_score = topology_score;
        muon_track_score.id = "muon_track_shower_score";
        muon_track_score.title = ";Muon Candidate Track Shower Score;Events";
        muon_track_score.expr =
            "([](const ROOT::RVec<float>& scores) -> float { return scores.empty() ? 0.f : ROOT::VecOps::Max(scores); })"
            "(track_shower_scores)";

        rarexsec::plot::TH1DModel muon_llr = topology_score;
        muon_llr.id = "muon_trk_llr_pid_v";
        muon_llr.title = ";Muon Candidate Track LLR PID;Events";
        muon_llr.xmin = -1.0;
        muon_llr.xmax = 1.0;
        muon_llr.expr =
            "([](const ROOT::RVec<float>& scores, const ROOT::RVec<float>& llrs) -> float {"
            " if (scores.empty()) { return -1.f; }"
            " auto idx = ROOT::VecOps::ArgMax(scores);"
            " return idx < llrs.size() ? llrs[idx] : -1.f; })"
            "(track_shower_scores, trk_llr_pid_v)";

        rarexsec::plot::TH1DModel muon_length = topology_score;
        muon_length.id = "muon_track_length";
        muon_length.title = ";Muon Candidate Track Length [cm];Events";
        muon_length.xmax = 200.0;
        muon_length.expr =
            "([](const ROOT::RVec<float>& scores, const ROOT::RVec<float>& lengths) -> float {"
            " if (scores.empty()) { return 0.f; }"
            " auto idx = ROOT::VecOps::ArgMax(scores);"
            " return idx < lengths.size() ? lengths[idx] : 0.f; })"
            "(track_shower_scores, track_length)";

        rarexsec::plot::TH1DModel muon_distance = topology_score;
        muon_distance.id = "muon_track_distance_to_vertex";
        muon_distance.title = ";Muon Candidate Track Distance to Vertex [cm];Events";
        muon_distance.xmax = 10.0;
        muon_distance.expr =
            "([](const ROOT::RVec<float>& scores, const ROOT::RVec<float>& distances) -> float {"
            " if (scores.empty()) { return 0.f; }"
            " auto idx = ROOT::VecOps::ArgMax(scores);"
            " return idx < distances.size() ? distances[idx] : 0.f; })"
            "(track_shower_scores, track_distance_to_vertex)";

        rarexsec::plot::TH1DModel muon_generation = topology_score;
        muon_generation.id = "muon_pfp_generation";
        muon_generation.title = ";Muon Candidate PFParticle Generation;Events";
        muon_generation.nbins = 9;
        muon_generation.xmin = -1.5;
        muon_generation.xmax = 7.5;
        muon_generation.expr =
            "([](const ROOT::RVec<float>& scores, const ROOT::RVec<unsigned>& gens) -> float {"
            " if (scores.empty()) { return -1.f; }"
            " auto idx = ROOT::VecOps::ArgMax(scores);"
            " return idx < gens.size() ? static_cast<float>(gens[idx]) : -1.f; })"
            "(track_shower_scores, pfp_generations)";

        const std::vector<rarexsec::plot::TH1DModel> specs = {
            beam_pe,
            veto_pe,
            software_trigger,
            num_slices,
            topology_score,
            fiducial,
            contained,
            cluster,
            muon_track_score,
            muon_llr,
            muon_length,
            muon_distance,
            muon_generation,
        };
        for (const auto& spec : specs) {
            plotter.draw_stack_by_channel(spec, mc_samples);
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
