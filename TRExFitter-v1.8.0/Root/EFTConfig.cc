#include "TRExFitter/EFTConfig.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"

#include <algorithm>

EFTConfig::EFTConfig() :
  fEFTOrder(EFTConfig::EFTOrder::ALL),
  fSplitSamplesPerBin(false),
  fEFTPlotRatio(true)
{
}

void EFTConfig::AddParamMinMax(const std::string& param, const double min, const double max) {
    auto itr = fParamMinMax.find(param);
    if (itr != fParamMinMax.end()) {
        LOG(WARNING) << "Parameter: " << param << " already found in the map, ignoring\n";
        return;
    }

    if (min > max) {
        LOG(ERROR) << "Parameter: " << param << " min value is grater than max value\n";
        exit(EXIT_FAILURE);
    }

    if (min > 0 || max < 0) {
        LOG(WARNING) << "Parameter: " << param << " min value > 0 or max value < 0. This is weird for an EFT fit.\n";
    }

    fParamMinMax.insert({param, std::make_pair(min, max)});
}

std::pair<double, double> EFTConfig::GetParMinMax(const std::string& param) const {
    auto itr = fParamMinMax.find(param);
    if (itr != fParamMinMax.end()) {
        return itr->second;
    }

    LOG(DEBUG) << "Parameter: " << param << ", no value found in the map, setting the min max to: -50,50\n";
    return std::make_pair(-50, 50);
}

void EFTConfig::AddExcludedRegionSample(const std::string& param,
                                        const std::string& region,
                                        const std::string& sample) {

    const std::vector<std::string> matchedRegions = Common::MatchingElememnts(fRegions, region);
    const std::vector<std::string> matchedSamples = Common::MatchingElememnts(fReferenceSamples, sample);

    auto itrPar = fExcludedRegionSample.find(param);
    if (itrPar == fExcludedRegionSample.end()) {
        std::map<std::string, std::vector<std::string> > tmp;
        for (const auto& ireg : matchedRegions) {
            tmp.insert(std::make_pair(ireg, matchedSamples));
        }

        fExcludedRegionSample.insert(std::make_pair(param, std::move(tmp)));
    } else {
        for (const auto& ireg : matchedRegions) {
            auto itrReg = itrPar->second.find(ireg);
            if (itrReg == itrPar->second.end()) {
                itrPar->second.insert(std::make_pair(ireg, matchedSamples));
            } else {
                for (const auto& isample : matchedSamples) {
                    auto itrSample = std::find(itrReg->second.begin(), itrReg->second.end(), isample);
                    if (itrSample == itrReg->second.end()) {
                        itrReg->second.emplace_back(isample);
                    }
                }
            }
        }
    }
}

bool EFTConfig::IsFullyExcluded(const std::string& region, const std::string& sample) const {

    auto itr = fOperatorsForSample.find(sample);
    if (itr == fOperatorsForSample.end()) {
        LOG(ERROR) << "Cannot find sample: " << sample << " in the list of known SM reference samples\n";
        exit(EXIT_FAILURE);
    }

    std::size_t excludedParams(0);

    for (const auto& ipar : itr->second) {
        if (EFTConfig::ParamIsExcluded(ipar, region, sample)) ++excludedParams;
    }

    return excludedParams == itr->second.size();
}

bool EFTConfig::ParamIsExcluded(const std::string& param,
                                const std::string& region,
                                const std::string& sample) const {

    auto itrPar = fExcludedRegionSample.find(param);
    if (itrPar == fExcludedRegionSample.end()) return false;

    auto itrReg = itrPar->second.find(region);
    if (itrReg == itrPar->second.end()) return false;

    auto itr = std::find(itrReg->second.begin(), itrReg->second.end(), sample);
    if (itr == itrReg->second.end()) return false;

    return true;
}

std::vector<std::pair<std::string, std::string> > EFTConfig::GetUpdatedParametrization(const std::vector<std::pair<std::string, std::string> >& originalPars,
                                                                                       const std::string& region,
                                                                                       const std::string& sample) const {
    std::vector<std::pair<std::string, std::string> > result;
    for (const auto& iformulae : originalPars) {
        std::string formula = iformulae.first;
        for (const auto& ipar : fParams) {
            if (iformulae.first.find(ipar) == std::string::npos || !ParamIsExcluded(ipar, region, sample)) continue;
            // need to to update the values
            formula = Common::ReplaceString(formula, ipar, "0");
        }

        // update the dependancy
        const auto dep = Common::processString(iformulae.second);
        std::string dependancy("");
        for (const auto& idep : dep) {
            const std::string par = idep.first;
            if (ParamIsExcluded(par, region, sample)) continue;
            std::string range = par+"[";
            for (const auto ivalue : idep.second) {
                range += std::to_string(ivalue) + ",";
            }
            range.resize(range.size() - 1);
            range += "],";
            dependancy += range;
        }

        if (dependancy.empty()) {
            LOG(ERROR) << "Something went wrong with the parsing of the strings for the EFT exclusion\n";
            exit(EXIT_FAILURE);
        }

        dependancy.resize(dependancy.size() - 1);

        result.emplace_back(formula, dependancy);
    }

    return result;
}