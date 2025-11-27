#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TSystem.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <rarexsec/Hub.h>
#include <rarexsec/proc/Env.h>
#include <rarexsec/proc/DataModel.h>

namespace {
std::string origin_to_string(rarexsec::sample::origin kind) {
    switch (kind) {
    case rarexsec::sample::origin::data:
        return "data";
    case rarexsec::sample::origin::beam:
        return "beam";
    case rarexsec::sample::origin::strangeness:
        return "strangeness";
    case rarexsec::sample::origin::ext:
        return "ext";
    case rarexsec::sample::origin::dirt:
        return "dirt";
    case rarexsec::sample::origin::unknown:
    default:
        return "unknown";
    }
}

void summarize_entries(const std::vector<const rarexsec::Entry*>& entries,
                       const std::string& label,
                       bool include_exposure) {
    std::cout << label << " samples: " << entries.size() << std::endl;

    long long total_events = 0;
    double total_pot_nom = 0.0;
    double total_pot_eqv = 0.0;
    double total_trig_nom = 0.0;
    double total_trig_eqv = 0.0;

    for (const auto* entry : entries) {
        if (!entry) {
            continue;
        }

        std::cout << "Sample kind '" << origin_to_string(entry->kind) << "' from file " << entry->file << std::endl;

        auto nominal_count = entry->rnode().Count().GetValue();
        total_events += nominal_count;
        std::cout << "  Nominal entries: " << nominal_count << std::endl;

        for (const auto& detvar : entry->detvars) {
            if (!detvar.second.node) {
                continue;
            }
            auto detvar_count = detvar.second.rnode().Count().GetValue();
            std::cout << "  Detector variation '" << detvar.first << "' entries: " << detvar_count << std::endl;
        }

        if (include_exposure) {
            total_pot_nom += entry->pot_nom;
            total_pot_eqv += (entry->pot_eqv > 0.0) ? entry->pot_eqv : entry->pot_nom;
            total_trig_nom += entry->trig_nom;
            total_trig_eqv += (entry->trig_eqv > 0.0) ? entry->trig_eqv : entry->trig_nom;
        }
    }

    std::cout << "Total " << label << " events: " << total_events << std::endl;
    if (include_exposure) {
        std::cout << "Total POT (nominal): " << total_pot_nom << std::endl;
        std::cout << "Total POT (equivalent): " << total_pot_eqv << std::endl;
        std::cout << "Total triggers (nominal): " << total_trig_nom << std::endl;
        std::cout << "Total triggers (equivalent): " << total_trig_eqv << std::endl;
    }
}
}

void print_event_counts() {
    try {
        ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT();

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();

        const auto data_entries = hub.data_entries(env.beamline, env.periods);
        const auto sim_entries = hub.simulation_entries(env.beamline, env.periods);

        std::cout << "Loaded beamline " << env.beamline << " for";
        for (const auto& period : env.periods) {
            std::cout << ' ' << period;
        }
        std::cout << std::endl;

        summarize_entries(data_entries, "Data", /*include_exposure=*/false);
        summarize_entries(sim_entries, "Simulation", /*include_exposure=*/true);

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}

