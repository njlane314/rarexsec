#pragma once

#include "rarexsec/proc/DataModel.h"
#include "rarexsec/Processor.h"
#include <string>
#include <unordered_map>
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
    static ROOT::RDF::RNode apply_slice(ROOT::RDF::RNode node, const Entry& rec);

    using PeriodDB = std::unordered_map<std::string, std::vector<Entry>>;
    std::unordered_map<std::string, PeriodDB> db_;
};

}
