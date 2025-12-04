#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <rarexsec/proc/Env.h>
#include <rarexsec/Hub.h>
#include <rarexsec/plot/Plotter.h>

#include <cctype>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Same helper as in your inspect_simulation_samples macro
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

void plot_inf_score_first() {
    try {
        ROOT::EnableImplicitMT();

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        auto beamlines = get_beamlines(env);

        using SamplePtr =
            decltype(hub.simulation_entries(env.beamline, env.periods))::value_type;

        std::vector<SamplePtr> samples;
        for (const auto& bl : beamlines) {
            auto sub = hub.simulation_entries(bl, env.periods);
            samples.insert(samples.end(), sub.begin(), sub.end());
        }

        std::cout << "Loaded beamlines";
        for (const auto& bl : beamlines) {
            std::cout << ' ' << bl;
        }
        std::cout << " for";
        for (const auto& period : env.periods) {
            std::cout << ' ' << period;
        }
        std::cout << " with " << samples.size() << " simulation samples." << std::endl;

        double total_pot_nom = 0.0;
        double total_pot_eqv = 0.0;
        double total_trig_nom = 0.0;
        double total_trig_eqv = 0.0;

        std::vector<const rarexsec::Entry*> mc_entries;

        auto origin_to_string = [](rarexsec::sample::origin kind) {
            switch (kind) {
            case rarexsec::sample::origin::data:
                return "data";
            case rarexsec::sample::origin::beam:
                return "beam";
            case rarexsec::sample::origin::strangeness:
                return "strangeness";
            case rarexsec::sample::origin::ext:
                return "ext";
            case rarexsec::sample::origin::dirt:
                return "dirt";
            case rarexsec::sample::origin::unknown:
            default:
                return "unknown";
            }
        };

        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }

            std::cout << "Sample kind '" << origin_to_string(entry->kind)
                      << "' from file " << entry->file << std::endl;

            auto final_count = entry->rnode().Count();
            auto eval_final_count = final_count.GetValue();
            std::cout << "  Final selection entries: " << eval_final_count << std::endl;

            for (const auto& detvar : entry->detvars) {
                if (!detvar.second.node) {
                    continue;
                }
                auto detvar_count = detvar.second.rnode().Count().GetValue();
                std::cout << "  Detector variation '" << detvar.first
                          << "' entries: " << detvar_count << std::endl;
            }

            total_pot_nom += entry->pot_nom;
            total_pot_eqv += (entry->pot_eqv > 0.0) ? entry->pot_eqv : entry->pot_nom;
            total_trig_nom += entry->trig_nom;
            total_trig_eqv += (entry->trig_eqv > 0.0) ? entry->trig_eqv : entry->trig_nom;

            mc_entries.push_back(entry);
        }

        std::cout << "Total POT (nominal): " << total_pot_nom << std::endl;
        std::cout << "Total POT (equivalent): " << total_pot_eqv << std::endl;
        std::cout << "Total triggers (nominal): " << total_trig_nom << std::endl;
        std::cout << "Total triggers (equivalent): " << total_trig_eqv << std::endl;

        rarexsec::plot::Options opt;
        opt.out_dir = "plots_inf_score";
        opt.image_format = "png";
        opt.show_ratio = false;
        opt.show_ratio_band = false;
        opt.y_title = "Events";
        opt.x_title = "First inference score";
        opt.beamline = env.beamline;
        opt.periods = env.periods;

        rarexsec::plot::TH1DModel spec;
        spec.id = "inf_score_first";
        spec.name = "inf_score_first";
        spec.title = ";First inference score;Events";
        spec.expr = "inf_scores.empty() ? -1.0f : inf_scores[0]";
        spec.weight = "w_nominal";
        spec.nbins = 50;
        spec.xmin = 0.0;
        spec.xmax = 1.0;
        spec.sel = rarexsec::selection::Preset::InclusiveMuCC;

        rarexsec::plot::Plotter plotter(opt);
        plotter.draw_stack_by_channel(spec, mc_entries);

    } catch (const std::exception& ex) {
        std::cerr << "Error in plot_inf_score_first: " << ex.what() << std::endl;
    }
}
