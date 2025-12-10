#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TSystem.h>

#include <rarexsec/Hub.h>
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

    if (out.empty()) {
        out.push_back(env.beamline);
    }

    return out;
}

static std::string sample_label(const rarexsec::Entry& e) {
    using rarexsec::sample::origin;
    switch (e.kind) {
    case origin::data:
        return "data";
    case origin::ext:
        return "ext";
    case origin::dirt:
        return "dirt";
    case origin::strangeness:
        return "strangeness";
    case origin::beam:
    case origin::unknown:
    default:
        return "beam";
    }
}

static std::string sanitise(std::string s) {
    for (char& c : s) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_' || c == '.'))
            c = '_';
    }
    return s;
}

static std::string make_tree_name(const rarexsec::Entry& e, const std::string& detvar) {
    std::string name = sanitise(e.beamline) + "_" + sanitise(e.period) + "_" + sanitise(sample_label(e));
    if (!detvar.empty()) {
        name += "__" + sanitise(detvar);
    }
    return name;
}

static std::vector<std::string> intersect_cols(ROOT::RDF::RNode node, const std::vector<std::string>& wanted) {
    auto have = node.GetColumnNames();
    std::unordered_set<std::string> avail(have.begin(), have.end());

    std::vector<std::string> out;
    out.reserve(wanted.size());
    for (const auto& c : wanted) {
        if (avail.count(c)) {
            out.push_back(c);
        }
    }
    return out;
}

void snapshot_numu_selection() {
    try {
        ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT();

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        auto beamlines = get_beamlines(env);

        using SamplePtr = decltype(hub.simulation_entries(env.beamline, env.periods))::value_type;
        std::vector<SamplePtr> samples;
        for (const auto& bl : beamlines) {
            auto sub = hub.simulation_entries(bl, env.periods);
            samples.insert(samples.end(), sub.begin(), sub.end());
        }

        if (samples.empty()) {
            std::cout << "[snapshot] no simulation samples found for the requested configuration.\n";
            return;
        }

        std::filesystem::create_directories("snapshots");
        std::string outfile = "snapshots/numu_selection";
        if (!beamlines.empty()) {
            outfile += "_" + beamlines.front();
            for (std::size_t i = 1; i < beamlines.size(); ++i) {
                outfile += "-" + beamlines[i];
            }
        }
        for (const auto& period : env.periods) {
            outfile += "_" + period;
        }
        outfile += ".root";

        ROOT::RDF::RSnapshotOptions sopt;
        sopt.fOverwriteIfExists = true;
        sopt.fLazy = false;

        const std::vector<std::string> columns{
            "run",
            "subrun",
            "event",
            "w_nominal",
            "analysis_channels",
            "event_detector_image_u",
            "event_detector_image_v",
            "event_detector_image_w",
        };

        bool fileExists = std::filesystem::exists(outfile);
        auto snapshot_once = [&](ROOT::RDF::RNode node, const std::string& treeName) mutable {
            const auto cols = intersect_cols(node, columns);
            if (cols.empty()) {
                std::cout << "[snapshot] skipping tree '" << treeName
                          << "' because none of the requested columns are available.\n";
                return;
            }

            sopt.fMode = fileExists ? "UPDATE" : "RECREATE";
            node.Snapshot("analysis", outfile, cols, sopt).GetValue();
            fileExists = true;
        };

        const auto preset = rarexsec::selection::Preset::InclusiveMuCC;
        std::cout << "[snapshot] applying numu selection preset to " << samples.size()
                  << " simulation sample(s).\n";

        for (const auto* entry : samples) {
            if (!entry)
                continue;

            auto selected = rarexsec::selection::apply(entry->rnode(), preset, *entry);
            snapshot_once(selected, make_tree_name(*entry, ""));

            for (const auto& kv : entry->detvars) {
                const auto& tag = kv.first;
                const auto& dv = kv.second;
                auto dv_selected = rarexsec::selection::apply(dv.rnode(), preset, *entry);
                snapshot_once(dv_selected, make_tree_name(*entry, tag));
            }
        }

        if (fileExists) {
            std::cout << "[snapshot] wrote selection snapshots to " << outfile << "\n";
        } else {
            std::cout << "[snapshot] no snapshots were written.\n";
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
    }
}

