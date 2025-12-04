#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>

#include <cctype>
#include <TGraph.h>
#include <TStyle.h>
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

void evaluate_inclusive_mucc_cutflow_plot() {
    if (gSystem->Load("librarexsec") < 0) {
        std::cerr << "Failed to load librexsec\n";
        return;
    }

    //gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

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

    // Denominator: total signal truth (no selection)
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

    std::vector<std::string> labels;
    std::vector<double> effs, purs;

    std::string label;
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

        labels.push_back(label);
        effs.push_back(eff);
        purs.push_back(pur);
    }

    // --- Plot 1: Efficiency vs Stage (with bin labels) ---
    const int n = static_cast<int>(effs.size());
    TH1F h("h_eff_vs_stage",
           "Inclusive #mu#nu CC: Efficiency vs Stage;Stage;Efficiency",
           n, 0.0, static_cast<double>(n));
    for (int i = 0; i < n; ++i) {
        h.SetBinContent(i + 1, effs[i]);
        h.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
    h.SetMinimum(0.0);
    h.SetMaximum(1.0);

    TCanvas c1("c1", "eff_vs_stage", 900, 500);
    c1.SetBottomMargin(0.25);
    h.Draw("HIST");
    c1.SaveAs("inclusive_mucc_eff_vs_stage.pdf");

    // --- Plot 2: Efficiency vs Purity (scatter) ---
    TGraph gr(n, purs.data(), effs.data());
    gr.SetTitle("Inclusive #mu#nu CC: Efficiency vs Purity;Purity;Efficiency");
    gr.SetMarkerStyle(20);

    TCanvas c2("c2", "eff_vs_purity", 700, 600);
    gr.Draw("AP");
    c2.SaveAs("inclusive_mucc_eff_vs_purity.pdf");

    // Minimal confirmation to stdout
    std::cout << "Wrote: inclusive_mucc_eff_vs_stage.pdf, inclusive_mucc_eff_vs_purity.pdf\n";
}
