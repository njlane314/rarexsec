#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TSystem.h>
#include <rarexsec/proc/Env.h>
#include <rarexsec/Hub.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void inspect_simulation_samples() {
    try {
        //ROOT::EnableThreadSafety();
        ROOT::EnableImplicitMT();

        const auto env = rarexsec::Env::from_env();
        auto hub = env.make_hub();
        const auto samples = hub.simulation_entries(env.beamline, env.periods);

        std::cout << "Loaded beamline " << env.beamline << " for";
        for (const auto& period : env.periods) {
            std::cout << ' ' << period;
        }
        std::cout << " with " << samples.size() << " simulation samples." << std::endl;

        auto origin_to_string = [](rarexsec::sample::origin kind) {
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
        };

        double total_pot_nom = 0.0;
        double total_pot_eqv = 0.0;
        double total_trig_nom = 0.0;
        double total_trig_eqv = 0.0;

        for (const auto* entry : samples) {
            if (!entry) {
                continue;
            }

            std::cout << "Sample kind '" << origin_to_string(entry->kind) << "' from file " << entry->file << std::endl;

            auto final_count = entry->rnode().Count();
            auto eval_final_count = final_count.GetValue();
            std::cout << "  Final selection entries: " << eval_final_count << std::endl;

            for (const auto& detvar : entry->detvars) {
                if (!detvar.second.node) {
                    continue;
                }
                auto detvar_count = detvar.second.rnode().Count().GetValue();
                std::cout << "  Detector variation '" << detvar.first << "' entries: " << detvar_count << std::endl;
            }

            total_pot_nom += entry->pot_nom;
            total_pot_eqv += (entry->pot_eqv > 0.0) ? entry->pot_eqv : entry->pot_nom;
            total_trig_nom += entry->trig_nom;
            total_trig_eqv += (entry->trig_eqv > 0.0) ? entry->trig_eqv : entry->trig_nom;
        }

        std::cout << "Total POT (nominal): " << total_pot_nom << std::endl;
        std::cout << "Total POT (equivalent): " << total_pot_eqv << std::endl;
        std::cout << "Total triggers (nominal): " << total_trig_nom << std::endl;
        std::cout << "Total triggers (equivalent): " << total_trig_eqv << std::endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}
