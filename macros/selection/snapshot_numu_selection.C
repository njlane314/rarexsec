#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TSystem.h>

#include <rarexsec/Hub.h>
#include <rarexsec/proc/Env.h>
#include <rarexsec/proc/Selection.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
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

static constexpr std::size_t kTrainingNSignal = 50000;
static constexpr std::size_t kTrainingNBackground = 50000;
static constexpr std::uint64_t kTrainingSeed = 12345u;

struct TrainCandidate {
    std::uint64_t event_key;
    double key;
};

static std::uint64_t make_event_key(int run, int sub, int evt) {
    const std::uint64_t r = static_cast<std::uint64_t>(static_cast<std::uint32_t>(run));
    const std::uint64_t s = static_cast<std::uint64_t>(static_cast<std::uint32_t>(sub));
    const std::uint64_t e = static_cast<std::uint64_t>(static_cast<std::uint32_t>(evt));
    return (r << 42) ^ (s << 21) ^ e;
}

static double stable_uniform(int run, int sub, int evt) {
    std::uint64_t key = make_event_key(run, sub, evt) ^ kTrainingSeed;
    std::uint64_t h = std::hash<std::uint64_t>{}(key);
    const double inv = 1.0 / static_cast<double>(std::numeric_limits<std::uint64_t>::max());
    return (static_cast<double>(h) + 0.5) * inv;
}

static std::unordered_set<std::uint64_t> select_top_ids(std::vector<TrainCandidate>& cands, std::size_t n) {
    std::unordered_set<std::uint64_t> out;
    if (n == 0 || cands.empty())
        return out;
    if (cands.size() <= n) {
        out.reserve(cands.size());
        for (const auto& c : cands)
            out.insert(c.event_key);
        return out;
    }
    auto nth = cands.end() - static_cast<std::vector<TrainCandidate>::difference_type>(n);
    std::nth_element(cands.begin(), nth, cands.end(),
                     [](const TrainCandidate& a, const TrainCandidate& b) {
                         return a.key < b.key;
                     });
    double cutoff = nth->key;
    out.reserve(n);
    std::size_t count = 0;
    for (const auto& c : cands) {
        if (c.key > cutoff) {
            if (out.insert(c.event_key).second)
                ++count;
        }
    }
    if (count < n) {
        for (const auto& c : cands) {
            if (c.key == cutoff) {
                if (out.insert(c.event_key).second)
                    ++count;
                if (count == n)
                    break;
            }
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

        if (std::filesystem::exists(outfile))
            std::filesystem::remove(outfile);

        ROOT::RDF::RSnapshotOptions sopt;
        sopt.fOverwriteIfExists = true;
        sopt.fLazy = false;

        const std::vector<std::string> columns{
            "run",
            "sub",
            "evt",
            "w_nominal",
            "analysis_channels",
            "detector_image_u",
            "detector_image_v",
            "detector_image_w",
        };

        bool fileExists = false;
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

        std::vector<TrainCandidate> signal_candidates;
        std::vector<TrainCandidate> background_candidates;

        for (const auto* entry : samples) {
            if (!entry)
                continue;

            auto base = rarexsec::selection::apply(entry->rnode(), preset, *entry)
                            .Define("event_key", make_event_key, {"run", "sub", "evt"})
                            .Define("valid_weight",
                                    [](float w) { return std::isfinite(w) && w > 0.0f; },
                                    {"w_nominal"})
                            .Filter("valid_weight")
                            .Define("u_rand", stable_uniform, {"run", "sub", "evt"})
                            .Define("es_key",
                                    [](double u, float w) { return std::log(u) / static_cast<double>(w); },
                                    {"u_rand", "w_nominal"});

            auto sig_df = base.Filter("is_signal");
            auto bkg_df = base.Filter("!is_signal");

            auto sig_keys = sig_df.Take<std::uint64_t>("event_key");
            auto sig_es = sig_df.Take<double>("es_key");

            auto bkg_keys = bkg_df.Take<std::uint64_t>("event_key");
            auto bkg_es = bkg_df.Take<double>("es_key");

            signal_candidates.reserve(signal_candidates.size() + sig_keys->size());
            for (std::size_t i = 0; i < sig_keys->size(); ++i) {
                signal_candidates.push_back(TrainCandidate{(*sig_keys)[i], (*sig_es)[i]});
            }

            background_candidates.reserve(background_candidates.size() + bkg_keys->size());
            for (std::size_t i = 0; i < bkg_keys->size(); ++i) {
                background_candidates.push_back(TrainCandidate{(*bkg_keys)[i], (*bkg_es)[i]});
            }
        }

        auto signal_ids = select_top_ids(signal_candidates, kTrainingNSignal);
        auto background_ids = select_top_ids(background_candidates, kTrainingNBackground);

        auto make_training_key = [](int run, int sub, int evt) {
            return make_event_key(run, sub, evt);
        };

        auto filter_training = [&](std::uint64_t key, bool is_signal) {
            if (is_signal)
                return signal_ids.find(key) != signal_ids.end();
            return background_ids.find(key) != background_ids.end();
        };

        std::size_t sample_index = 0;
        for (const auto* entry : samples) {
            ++sample_index;
            if (!entry)
                continue;

            const auto tree_name = make_tree_name(*entry, "");
            std::cout << "[snapshot] [" << sample_index << "/" << samples.size() << "] processing sample '"
                      << tree_name << "' with " << entry->detvars.size() << " detvar variation(s).\n";

            auto selected = rarexsec::selection::apply(entry->rnode(), preset, *entry)
                                .Define("training_event_key", make_training_key, {"run", "sub", "evt"})
                                .Filter(filter_training, {"training_event_key", "is_signal"});
            snapshot_once(selected, tree_name);

            for (const auto& kv : entry->detvars) {
                const auto& tag = kv.first;
                const auto& dv = kv.second;

                std::cout << "[snapshot]      detvar '" << tag << "'\n";

                auto dv_selected = rarexsec::selection::apply(dv.rnode(), preset, *entry)
                                       .Define("training_event_key", make_training_key, {"run", "sub", "evt"})
                                       .Filter(filter_training, {"training_event_key", "is_signal"});
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

