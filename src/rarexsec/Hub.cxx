#include "rarexsec/Hub.h"
#include "rarexsec/Processor.h"
#include "rarexsec/proc/Volume.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <nlohmann/json.hpp>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

using json = nlohmann::json;

namespace {

constexpr std::size_t kTrainingNSignal = 50000;
constexpr std::size_t kTrainingNBackground = 50000;
constexpr std::uint64_t kTrainingSeed = 12345u;

struct TrainCandidate {
    std::uint64_t event_key = 0;
    double key = 0.0;
    double weight = 0.0;
};

struct TrainSelection {
    std::unordered_set<std::uint64_t> ids;
    double sum_weight = 0.0;
};

std::uint64_t make_event_key(int run, int sub, int evt)
{
    const std::uint64_t r = static_cast<std::uint64_t>(static_cast<std::uint32_t>(run));
    const std::uint64_t s = static_cast<std::uint64_t>(static_cast<std::uint32_t>(sub));
    const std::uint64_t e = static_cast<std::uint64_t>(static_cast<std::uint32_t>(evt));
    return (r << 42) ^ (s << 21) ^ e;
}

double stable_uniform(std::uint64_t key)
{
    const std::uint64_t h = std::hash<std::uint64_t>{}(key ^ kTrainingSeed);
    const double inv = 1.0 / static_cast<double>(std::numeric_limits<std::uint64_t>::max());
    return (static_cast<double>(h) + 0.5) * inv;
}

TrainSelection select_top(std::vector<TrainCandidate>& cands, std::size_t n)
{
    TrainSelection out;
    if (cands.empty() || n == 0)
        return out;

    if (cands.size() <= n) {
        out.ids.reserve(cands.size());
        for (const auto& c : cands) {
            out.ids.insert(c.event_key);
            out.sum_weight += c.weight;
        }
        return out;
    }

    auto nth = cands.end() - static_cast<std::vector<TrainCandidate>::difference_type>(n);
    std::nth_element(cands.begin(), nth, cands.end(),
                     [](const TrainCandidate& a, const TrainCandidate& b) {
                         return a.key < b.key;
                     });

    const double cutoff = nth->key;
    out.ids.reserve(n);

    std::size_t count = 0;
    for (const auto& c : cands) {
        if (c.key > cutoff) {
            if (out.ids.insert(c.event_key).second) {
                out.sum_weight += c.weight;
                ++count;
            }
        }
    }

    if (count < n) {
        for (const auto& c : cands) {
            if (c.key == cutoff) {
                if (out.ids.insert(c.event_key).second) {
                    out.sum_weight += c.weight;
                    ++count;
                    if (count == n)
                        break;
                }
            }
        }
    }

    return out;
}

}

//____________________________________________________________________________
static std::string to_lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}
//____________________________________________________________________________
static rarexsec::Slice parse_slice_opt(const json& j)
{
    if (!j.contains("slice"))
        return rarexsec::Slice::None;
    const auto s = to_lower(j.at("slice").get<std::string>());
    if (s == "beam" || s == "beaminclusive")
        return rarexsec::Slice::BeamInclusive;
    if (s == "strange" || s == "strangeness" || s == "strangenessinclusive")
        return rarexsec::Slice::StrangenessInclusive;
    throw std::runtime_error("unknown slice: " + s);
}
//____________________________________________________________________________
static std::pair<rarexsec::Source, rarexsec::Slice>
parse_kind_slice(const std::string& kind, const json& s)
{
    if (kind == "data")
        return {rarexsec::Source::Data, rarexsec::Slice::None};
    if (kind == "ext" || kind == "external")
        return {rarexsec::Source::Ext, rarexsec::Slice::None};
    if (kind == "mc")
        return {rarexsec::Source::MC, parse_slice_opt(s)};
    if (kind == "beam")
        return {rarexsec::Source::MC, rarexsec::Slice::BeamInclusive};
    if (kind == "strangeness")
        return {rarexsec::Source::MC, rarexsec::Slice::StrangenessInclusive};
    if (kind == "dirt")
        return {rarexsec::Source::MC, rarexsec::Slice::None};
    throw std::runtime_error("unknown kind: " + kind);
}
//____________________________________________________________________________
rarexsec::Frame rarexsec::Hub::sample(const Entry& rec) const
{
    static const std::string tree = "nuselection/EventSelectionFilter";
    auto df_ptr = std::make_shared<ROOT::RDataFrame>(tree, rec.files);
    ROOT::RDF::RNode node = *df_ptr;

    node = processor().run(node, rec);
    node = apply_slice(node, rec);

    using rarexsec::Source;

    if (rec.source == Source::MC) {
        const auto* sig_ids = &training_signal_ids_;
        const auto* bkg_ids = &training_background_ids_;
        const double sig_scale = signal_analysis_scale_;
        const double bkg_scale = background_analysis_scale_;

        node = node.Define(
            "is_training",
            [sig_ids, bkg_ids](int run, int sub, int evt, bool is_signal) {
                const auto key = make_event_key(run, sub, evt);
                if (is_signal)
                    return sig_ids->find(key) != sig_ids->end();
                return bkg_ids->find(key) != bkg_ids->end();
            },
            {"run", "sub", "evt", "is_signal"});

        node = node.Define(
            "w_analysis",
            [sig_scale, bkg_scale](float w_nominal, bool is_signal, bool is_training) {
                if (!std::isfinite(w_nominal) || w_nominal <= 0.0f)
                    return 0.0f;
                if (is_training)
                    return 0.0f;
                const double scale = is_signal ? sig_scale : bkg_scale;
                const double w = static_cast<double>(w_nominal) * scale;
                if (!std::isfinite(w) || w <= 0.0)
                    return 0.0f;
                return static_cast<float>(w);
            },
            {"w_nominal", "is_signal", "is_training"});
    } else {
        node = node.Define("is_training", [] { return false; });
        node = node.Define(
            "w_analysis",
            [](float w_nominal) { return w_nominal; },
            {"w_nominal"});
    }

    return Frame{df_ptr, std::move(node)};
}
//____________________________________________________________________________
rarexsec::Hub::Hub(const std::string& path)
{
    std::ifstream cfg(path);
    if (!cfg)
        throw std::runtime_error("cannot open " + path);
    json j;
    cfg >> j;

    const auto& bl = j.at("beamlines");
    for (auto it_bl = bl.begin(); it_bl != bl.end(); ++it_bl) {
        const std::string beamline = it_bl.key();
        const auto& runs = it_bl.value();

        for (auto it_r = runs.begin(); it_r != runs.end(); ++it_r) {
            const std::string period = it_r.key();
            const auto& arr = it_r.value().at("samples");
            auto& bucket = db_[beamline][period];

            for (const auto& s : arr) {
                Entry rec;
                rec.beamline = beamline;
                rec.period = period;

                const auto kind_str = to_lower(s.at("kind").get<std::string>());
                std::tie(rec.source, rec.slice) = parse_kind_slice(kind_str, s);
                rec.kind = (kind_str == "dirt") ? sample::origin::dirt
                                                : sample::from_source_slice(rec.source, rec.slice);

                if (s.contains("files")) {
                    rec.files = s.at("files").get<std::vector<std::string>>();
                } else if (s.contains("file")) {
                    const auto f = s.at("file").get<std::string>();
                    rec.files = f.empty() ? std::vector<std::string>{} : std::vector<std::string>{f};
                } else {
                    throw std::runtime_error("sample missing 'file' or 'files'");
                }
                if (rec.files.empty())
                    throw std::runtime_error("empty 'files' for sample in " + beamline + "/" + period);
                rec.file = rec.files.front();

                if (rec.source == Source::Ext) {
                    rec.trig_nom = s.value("trig", 0.0);
                    rec.trig_eqv = s.value("trig_eff", 0.0);
                } else if (rec.source == Source::MC) {
                    rec.pot_nom = s.value("pot", 0.0);
                    rec.pot_eqv = s.value("pot_eff", 0.0);
                }

                rec.nominal = sample(rec);

                if (s.contains("detvars")) {
                    const auto& dvs = s.at("detvars");
                    for (auto it_dv = dvs.begin(); it_dv != dvs.end(); ++it_dv) {
                        const std::string tag = it_dv.key();
                        const auto& desc = it_dv.value();
                        std::vector<std::string> dv_files;
                        if (desc.contains("files"))
                            dv_files = desc.at("files").get<std::vector<std::string>>();
                        else if (desc.contains("file"))
                            dv_files = {desc.at("file").get<std::string>()};

                        if (!dv_files.empty()) {
                            Entry dv = rec;
                            dv.files = std::move(dv_files);
                            dv.file = dv.files.front();
                            rec.detvars.emplace(tag, sample(dv));
                        }
                    }
                }

                bucket.push_back(std::move(rec));
            }
        }
    }

    build_training_sets();
}
//____________________________________________________________________________
ROOT::RDF::RNode rarexsec::Hub::apply_slice(ROOT::RDF::RNode node, const Entry& rec)
{
    using rarexsec::Slice;
    using rarexsec::Source;

    if (rec.source == Source::MC) {
        if (rec.slice == Slice::StrangenessInclusive)
            return node.Filter([](bool s) { return s; }, {"is_strange"});
        if (rec.slice == Slice::BeamInclusive)
            return node.Filter([](bool s) { return !s; }, {"is_strange"});
        return node;
    }
    if (rec.slice != Slice::None) {
        throw std::runtime_error("Slice requested for non-MC sample at " + rec.beamline + "/" + rec.period);
    }
    return node;
}
//____________________________________________________________________________
std::vector<const rarexsec::Entry*>
rarexsec::Hub::simulation_entries(const std::string& beamline,
                                  const std::vector<std::string>& periods) const
{
    std::vector<const Entry*> out;
    auto it_bl = db_.find(beamline);
    if (it_bl == db_.end())
        return out;
    for (const auto& per : periods) {
        auto it_p = it_bl->second.find(per);
        if (it_p == it_bl->second.end())
            continue;
        for (const auto& rec : it_p->second) {
            if (rec.source != Source::Data)
                out.push_back(&rec);
        }
    }
    return out;
}
//____________________________________________________________________________
std::vector<const rarexsec::Entry*>
rarexsec::Hub::data_entries(const std::string& beamline,
                            const std::vector<std::string>& periods) const
{
    std::vector<const Entry*> out;
    auto it_bl = db_.find(beamline);
    if (it_bl == db_.end())
        return out;
    for (const auto& per : periods) {
        auto it_p = it_bl->second.find(per);
        if (it_p == it_bl->second.end())
            continue;
        for (const auto& rec : it_p->second)
            if (rec.source == Source::Data)
                out.push_back(&rec);
    }
    return out;
}
//____________________________________________________________________________
void rarexsec::Hub::build_training_sets()
{
    using rarexsec::Source;

    std::vector<TrainCandidate> sig_candidates;
    std::vector<TrainCandidate> bkg_candidates;
    double total_sig_weight = 0.0;
    double total_bkg_weight = 0.0;

    for (const auto& bl_pair : db_) {
        const auto& period_db = bl_pair.second;
        for (const auto& per_pair : period_db) {
            const auto& entries = per_pair.second;
            for (const auto& rec : entries) {
                if (rec.source != Source::MC)
                    continue;

                static const std::string tree = "nuselection/EventSelectionFilter";
                ROOT::RDataFrame df(tree, rec.files);
                ROOT::RDF::RNode node = df;

                node = processor().run(node, rec);
                node = apply_slice(node, rec);

                std::mutex mtx;
                node.Foreach(
                    [&](int run, int sub, int evt, bool is_signal, float w_nominal) {
                        if (!std::isfinite(w_nominal) || w_nominal <= 0.0f)
                            return;

                        const auto ev_key = make_event_key(run, sub, evt);
                        const double u = stable_uniform(ev_key);
                        if (!(u > 0.0 && u < 1.0))
                            return;

                        const double w = static_cast<double>(w_nominal);
                        const double key = std::log(u) / w;

                        TrainCandidate cand{ev_key, key, w};

                        std::lock_guard<std::mutex> lock(mtx);
                        if (is_signal) {
                            sig_candidates.emplace_back(cand);
                            total_sig_weight += w;
                        } else {
                            bkg_candidates.emplace_back(cand);
                            total_bkg_weight += w;
                        }
                    },
                    {"run", "sub", "evt", "is_signal", "w_nominal"});
            }
        }
    }

    const auto sig_sel = select_top(sig_candidates, kTrainingNSignal);
    const auto bkg_sel = select_top(bkg_candidates, kTrainingNBackground);

    training_signal_ids_ = sig_sel.ids;
    training_background_ids_ = bkg_sel.ids;

    const double sig_analysis_weight = total_sig_weight - sig_sel.sum_weight;
    const double bkg_analysis_weight = total_bkg_weight - bkg_sel.sum_weight;

    if (sig_analysis_weight > 0.0)
        signal_analysis_scale_ = total_sig_weight / sig_analysis_weight;
    else
        signal_analysis_scale_ = 1.0;

    if (bkg_analysis_weight > 0.0)
        background_analysis_scale_ = total_bkg_weight / bkg_analysis_weight;
    else
        background_analysis_scale_ = 1.0;
}
//____________________________________________________________________________
