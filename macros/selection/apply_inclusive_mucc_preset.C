#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <rarexsec/Hub.h>
#include <rarexsec/proc/DataModel.h>
#include <rarexsec/proc/Env.h>
#include <rarexsec/proc/Selection.h>

#include <cctype>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

void apply_inclusive_mucc_preset() {
    try {
        ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT();

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        auto beamlines = get_beamlines(env);

        const auto preset = rarexsec::selection::Preset::InclusiveMuCC;
        auto preset_to_string = [](rarexsec::selection::Preset value) {
            switch (value) {
            case rarexsec::selection::Preset::Empty:
                return "Empty";
            case rarexsec::selection::Preset::Trigger:
                return "Trigger";
            case rarexsec::selection::Preset::Slice:
                return "Slice";
            case rarexsec::selection::Preset::Fiducial:
                return "Fiducial";
            case rarexsec::selection::Preset::Topology:
                return "Topology";
            case rarexsec::selection::Preset::Muon:
                return "Muon";
            case rarexsec::selection::Preset::InclusiveMuCC:
            default:
                return "InclusiveMuCC";
            }
        };

        std::cout << "Using preset: " << preset_to_string(preset) << "\n";

        using SamplePtr = decltype(hub.simulation_entries(env.beamline, env.periods))::value_type;
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

        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }

            auto node = rarexsec::selection::apply(entry->rnode(), preset, *entry);
            const auto selected = node.Count().GetValue();
            std::cout << "Sample '" << entry->file << "' selected entries: " << selected << std::endl;
        }

        const auto eval = rarexsec::selection::evaluate(
            samples,
            [](int ch) {
                const auto channel = static_cast<rarexsec::Channel>(ch);
                switch (channel) {
                case rarexsec::Channel::MuCC0pi_ge1p:
                case rarexsec::Channel::MuCC1pi:
                case rarexsec::Channel::MuCCPi0OrGamma:
                case rarexsec::Channel::MuCCNpi:
                case rarexsec::Channel::MuCCOther:
                    return true;
                default:
                    return false;
                }
            },
            preset);

        std::cout << "Selection evaluation:\n"
                  << "  Denominator (signal truth): " << eval.denom << '\n'
                  << "  Selected (all): " << eval.selected << '\n'
                  << "  Selected signal: " << eval.numer << '\n'
                  << "  Efficiency: " << eval.efficiency() << '\n'
                  << "  Purity: " << eval.purity() << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
