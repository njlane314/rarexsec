#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <rarexsec/Hub.h>
#include <rarexsec/proc/DataModel.h>
#include <rarexsec/proc/Env.h>
#include <rarexsec/proc/Selection.h>

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

void evaluate_cutflow() {
    const auto env = rarexsec::Env::from_env();
    auto hub = env.make_hub();
    auto beamlines = get_beamlines(env);
    using SamplePtr = decltype(hub.simulation_entries(env.beamline, env.periods))::value_type;
    std::vector<SamplePtr> mc;
    for (const auto& bl : beamlines) {
        auto sub = hub.simulation_entries(bl, env.periods);
        mc.insert(mc.end(), sub.begin(), sub.end());
    }

    auto is_signal = [](int ch_int) {
        const auto ch = static_cast<rarexsec::Channel>(ch_int);
        switch (ch) {
        case rarexsec::Channel::MuCC0pi_ge1p:
        case rarexsec::Channel::MuCC1pi:
        case rarexsec::Channel::MuCCPi0OrGamma:
        case rarexsec::Channel::MuCCNpi:
        case rarexsec::Channel::MuCCOther:
            return true;
        default:
            return false;
        }
    };

    auto sumw = [](ROOT::RDF::RNode n) {
        return static_cast<double>(n.Sum<float>("w_nominal").GetValue());
    };

    double denom = 0.0;
    for (const auto* rec : mc) {
        auto base = rec->nominal.rnode();
        denom += sumw(base.Filter([&](int ch){ return is_signal(ch); }, {"analysis_channels"}));
    }

    using Preset = rarexsec::selection::Preset;
    const std::vector<std::pair<std::string, Preset>> atoms = {
        {"Trigger",  Preset::Trigger},
        {"Slice",    Preset::Slice},
        {"Fiducial", Preset::Fiducial},
        {"Topology", Preset::Topology},
        {"Muon",     Preset::Muon}
    };

    std::string label;
    std::cout.setf(std::ios::fixed);
    std::cout.precision(6);
    constexpr auto stage_width = 36;
    constexpr auto value_width = 16;

    std::cout << std::left << std::setw(stage_width) << "Stage"
              << std::right << std::setw(value_width) << "Denom(signal)"
              << std::setw(value_width) << "Selected(all)"
              << std::setw(value_width) << "Selected(signal)"
              << std::setw(value_width) << "Efficiency"
              << std::setw(value_width) << "Purity" << '\n';

    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (!label.empty()) label += "+";
        label += atoms[i].first;

        double sel_all = 0.0;
        double sel_sig = 0.0;

        for (const auto* rec : mc) {
            auto node = rec->nominal.rnode();
            for (std::size_t j = 0; j <= i; ++j)
                node = rarexsec::selection::apply(node, atoms[j].second, *rec);

            sel_all += sumw(node);
            sel_sig += sumw(node.Filter([&](int ch){ return is_signal(ch); }, {"analysis_channels"}));
        }

        const double eff = denom > 0.0 ? sel_sig / denom : 0.0;
        const double pur = sel_all > 0.0 ? sel_sig / sel_all : 0.0;

        std::cout << std::left << std::setw(stage_width) << label
                  << std::right << std::setw(value_width) << denom
                  << std::setw(value_width) << sel_all
                  << std::setw(value_width) << sel_sig
                  << std::setw(value_width) << eff
                  << std::setw(value_width) << pur << '\n';
    }
}
