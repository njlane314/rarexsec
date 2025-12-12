#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TSystem.h>

#include <rarexsec/Hub.h>
#include <rarexsec/proc/Snapshot.h>
#include <rarexsec/proc/Env.h>
#include <rarexsec/proc/Selection.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
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

    if (out.empty())
        out.push_back(env.beamline);

    return out;
}

static bool has_column(ROOT::RDF::RNode node, const std::string& name) {
    const auto cols = node.GetColumnNames();
    return std::find(cols.begin(), cols.end(), name) != cols.end();
}

void snapshot_numu_selection() {
    try {
        ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT();

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        const auto beamlines = get_beamlines(env);

        std::vector<const rarexsec::Entry*> samples;
        for (const auto& bl : beamlines) {
            auto sub = hub.simulation_entries(bl, env.periods);
            samples.insert(samples.end(), sub.begin(), sub.end());
        }

        if (samples.empty()) {
            std::cout << "[snapshot] no simulation samples found for the requested configuration.\n";
            return;
        }

        rarexsec::snapshot::Options opt;
        opt.outdir = "snapshots";
        opt.tree = "analysis";
        opt.outfile = "numu_selection_train";
        if (!beamlines.empty()) {
            opt.outfile += "_" + beamlines.front();
            for (std::size_t i = 1; i < beamlines.size(); ++i)
                opt.outfile += "-" + beamlines[i];
        }
        for (const auto& period : env.periods)
            opt.outfile += "_" + period;
        opt.outfile += ".root";

        opt.columns = {
            "run",
            "sub",
            "evt",
            "w_nominal",
            "is_signal",
            "analysis_channels",
            "detector_image_u",
            "detector_image_v",
            "detector_image_w",
        };

        std::filesystem::create_directories(opt.outdir);
        const std::string outFile = (std::filesystem::path(opt.outdir) / opt.outfile).string();
        if (std::filesystem::exists(outFile))
            std::filesystem::remove(outFile);

        const auto preset = rarexsec::selection::Preset::InclusiveMuCC;

        bool fileExists = false;

        for (const auto* entry : samples) {
            if (!entry)
                continue;

            auto node = entry->rnode();
            if (!has_column(node, "is_training"))
                throw std::runtime_error("missing required column: is_training");
            node = node.Filter([](bool t) { return t; }, {"is_training"});
            node = rarexsec::selection::apply(node, preset, *entry);

            const auto cols = rarexsec::snapshot::intersect_cols(node, opt.columns);
            if (cols.empty())
                continue;

            ROOT::RDF::RSnapshotOptions sopt;
            sopt.fMode = fileExists ? "UPDATE" : "RECREATE";
            sopt.fOverwriteIfExists = true;
            sopt.fLazy = false;

            node.Snapshot(opt.tree, outFile, cols, sopt).GetValue();
            fileExists = true;
        }

        if (fileExists) {
            std::cout << "[snapshot] wrote selection snapshots to " << outFile << "\n";
        } else {
            std::cout << "[snapshot] no snapshots were written.\n";
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
    }
}
