#pragma once

#include "rarexsec/proc/DataModel.h"
#include "rarexsec/Processor.h"
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace rarexsec {

class Hub {
  public:
    explicit Hub(const std::string& path);

    Frame sample(const Entry& rec) const;

    std::vector<const Entry*> simulation_entries(const std::string& beamline,
                                                 const std::vector<std::string>& periods) const;
    std::vector<const Entry*> data_entries(const std::string& beamline,
                                           const std::vector<std::string>& periods) const;

  private:
    void build_training_sets();

    static ROOT::RDF::RNode apply_slice(ROOT::RDF::RNode node, const Entry& rec);

    using PeriodDB = std::unordered_map<std::string, std::vector<Entry>>;
    std::unordered_map<std::string, PeriodDB> db_;

    std::unordered_set<std::uint64_t> training_signal_ids_;
    std::unordered_set<std::uint64_t> training_background_ids_;
    double signal_analysis_scale_ = 1.0;
    double background_analysis_scale_ = 1.0;
};

}
