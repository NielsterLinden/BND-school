#include "TRExFitter/ReparametrizationManager.h"
#include "TRExFitter/Logger.h"

#include <fstream>

void ReparametrizationManager::ReadReparametrizationFile(const std::string& path) {

    fParametrizationMap.clear();

    std::fstream in(path, std::ios_base::in);
    if (in.bad() || !in.is_open()) {
        LOG(ERROR) << "Cannot open the file at: " << path << "\n";
        exit(EXIT_FAILURE);
    }

    std::string reg;
    std::string smp;
    int bin;
    std::string formula;
    std::string dependency;

    while(in >> reg >> smp >> bin >> formula >> dependency) {
        const std::string key = reg + "_" + smp;
        auto itr = fParametrizationMap.find(key);
        if (itr == fParametrizationMap.end()) {
            std::vector<std::pair<std::string, std::string> > value;
            value.emplace_back(std::make_pair(formula, dependency));
            fParametrizationMap.insert(std::make_pair(key, value));
        } else {
            itr->second.emplace_back(std::make_pair(formula, dependency));
        }
    }

    in.close();
}

std::vector<std::pair<std::string, std::string> > ReparametrizationManager::GetParametrization(const std::string& region, const std::string& sample) const {
    const std::string key = region + "_" + sample;
    auto itr = fParametrizationMap.find(key);
    if (itr == fParametrizationMap.end()) {
        LOG(ERROR) << "Cannot find region: " << region << " and sample: " << sample << " in the reparemetrization map\n";
        exit(EXIT_FAILURE);
    }

    return itr->second;
}