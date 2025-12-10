#include "rarexsec/Hub.h"
#include "rarexsec/Processor.h"
#include "rarexsec/proc/Volume.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

using json = nlohmann::json;

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
