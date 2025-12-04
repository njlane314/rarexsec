#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>

#include <cctype>
#include <rarexsec/proc/Env.h>
#include <rarexsec/Hub.h>
#include <rarexsec/proc/Snapshot.h>

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

void write_simulation_snapshots() {
    try {
        ROOT::EnableImplicitMT();

        if (gSystem->Load("librarexsec") < 0) {
            throw std::runtime_error("Failed to load librexsec");
        }

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        auto beamlines = get_beamlines(env);
        using SamplePtr = decltype(hub.simulation_entries(env.beamline, env.periods))::value_type;
        std::vector<SamplePtr> samples;
        for (const auto& bl : beamlines) {
            auto sub = hub.simulation_entries(bl, env.periods);
            samples.insert(samples.end(), sub.begin(), sub.end());
        }

        rarexsec::snapshot::Options opt;
        opt.outdir = "snapshots";
        opt.tree = env.tree;
        std::string outfile;
        for (std::size_t i = 0; i < beamlines.size(); ++i) {
            if (i != 0) {
                outfile += "-";
            }
            outfile += beamlines[i];
        }
        for (const auto& period : env.periods) {
            outfile += "_" + period;
        }
        outfile += ".root";
        opt.outfile = outfile;

        auto outputs = rarexsec::snapshot::write(samples, opt);

        if (outputs.empty()) {
            std::cout << "[snapshot] no files were written (no matching samples?).\n";
        } else {
            std::cout << "[snapshot] wrote " << outputs.size() << " file(s):\n";
            for (const auto& f : outputs) {
                std::cout << "  " << f << "\n";
            }
        }

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
    }
}
