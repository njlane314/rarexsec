#pragma once
#include "rarexsec/Hub.h"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace rarexsec {
namespace snapshot {

struct Options {
    std::string outdir = "snapshots";
    std::string outfile = "all_samples.root";
    std::string tree = "analysis";
    std::vector<std::string> columns;
};

inline std::string source_to_string(Source s) {
    switch (s) {
    case Source::Data:
        return "data";
    case Source::Ext:
        return "ext";
    case Source::MC:
        return "mc";
    }
    return "unknown";
}

inline std::string slice_to_string(Slice s) {
    switch (s) {
    case Slice::None:
        return "none";
    case Slice::BeamInclusive:
        return "beam";
    case Slice::StrangenessInclusive:
        return "strangeness";
    }
    return "unknown";
}

inline std::string sample_label(const Entry& e) {
    if (e.kind == sample::origin::dirt)
        return "dirt";
    if (e.source == Source::MC) {
        const auto slice = slice_to_string(e.slice);
        if (slice == "none")
            return "mc";
        return slice;
    }
    return source_to_string(e.source);
}

inline std::string sanitise(std::string s) {
    for (char& c : s) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_' || c == '.'))
            c = '_';
    }
    return s;
}

inline const std::vector<std::string>& default_columns() {
    static const std::vector<std::string> cols{
        "run",
        "subrun",
        "event",
        "w_nominal",
        "analysis_channels",
    };
    return cols;
}

inline std::vector<std::string> intersect_cols(ROOT::RDF::RNode node, const std::vector<std::string>& wanted) {
    auto have = node.GetColumnNames();
    std::unordered_set<std::string> avail(have.begin(), have.end());
    const auto& req = wanted.empty() ? default_columns() : wanted;

    std::vector<std::string> out;
    out.reserve(req.size());
    for (const auto& c : req) {
        if (avail.count(c))
            out.push_back(c);
    }
    return out;
}

inline std::string make_out_path(const Options& opt, const Entry& e, const std::string& detvar) {
    const auto base = e.files.empty() ? std::string{}
                                      : std::filesystem::path(e.files.front()).filename().string();
    std::string name = sanitise(e.beamline) + "_" +
                       sanitise(e.period) + "_" +
                       sanitise(sample_label(e));
    if (!detvar.empty())
        name += "__" + sanitise(detvar);
    if (!base.empty())
        name += "__" + sanitise(base);

    std::filesystem::create_directories(opt.outdir);
    return (std::filesystem::path(opt.outdir) / name).string();
}

inline std::string make_out_file(const Options& opt) {
    std::filesystem::create_directories(opt.outdir);
    return (std::filesystem::path(opt.outdir) / opt.outfile).string();
}

inline std::string make_tree_name(const Options& opt, const Entry& e, const std::string& detvar) {
    std::string name = sanitise(opt.tree) + "_" + sanitise(sample_label(e));
    if (!detvar.empty())
        name += "__" + sanitise(detvar);
    return name;
}

inline std::vector<std::string> write(const std::vector<const Entry*>& samples,
                                      const Options& opt = {}) {
    std::vector<std::string> outputs;
    outputs.reserve(1);

    const std::string outFile = make_out_file(opt);
    bool fileExists = std::filesystem::exists(outFile);

    auto snapshot_once = [&](ROOT::RDF::RNode node,
                             const std::string& treeName,
                             const std::vector<std::string>& cols) {
        ROOT::RDF::RSnapshotOptions sopt;
        sopt.fMode = fileExists ? "UPDATE" : "RECREATE";
        sopt.fOverwriteIfExists = true;
        node.Snapshot(treeName, outFile, cols, sopt).GetValue();
        fileExists = true;
    };

    for (const Entry* e : samples) {
        if (!e)
            continue;

        const auto cols = intersect_cols(e->rnode(), opt.columns);
        const auto treeName = make_tree_name(opt, *e, "");
        snapshot_once(e->rnode(), treeName, cols);
    }

    if (!samples.empty())
        outputs.push_back(outFile);
    return outputs;
}

}
}
