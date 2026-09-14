#include "TRExFitter/Morphing.h"

#include "TRExFitter/Logger.h"

Morphing::Morphing() :
    fFitFunction("QUADRATIC"),
    fRangeScale(0.5),
    fChi2WarningLimit(1.5),
    fMake1DProjections(false)
{
}

void Morphing::AddStartValues(const std::string& target, const std::string& region, const int bin, const int index, const double val) {
    // bin == -999 indicates to set starting values across all bins
    auto itrTarget = fStartValues.find(target);
    if (itrTarget == fStartValues.end()) {
        std::map<std::string, std::map<int, std::map<int, double> > > tmp;
        std::map<int, std::map<int, double> > tmp2;
        std::map<int, double> tmp3;
        tmp3.insert({index, val});
        tmp2.insert({bin, tmp3});
        tmp.insert({region, tmp2});
        fStartValues.insert({target, tmp});
    } else {
        auto itr = itrTarget->second.find(region);
        if (itr == itrTarget->second.end()) {
            std::map<int, double> tmp;
            tmp.insert({index, val});
            std::map<int, std::map<int, double> > tmp2;
            tmp2.insert({bin, tmp});
            itrTarget->second.insert(std::make_pair(region, tmp2));
        }  else {
            auto itr2 = itr->second.find(bin);
            if (itr2 == itr->second.end()) {
                std::map<int, double> tmp;
                tmp.insert({index, val});
                itr->second.insert(std::make_pair(bin, tmp));
            } else {
                itr2->second.insert(std::make_pair(index, val));
            }
        }
    }
}

void Morphing::AddParamLimits(const std::string& target, const std::string& region, const int bin, const int index, const double min, const double max) {
    // bin == -999 indicates to set starting values across all bins
    auto itrTarget = fParamLimits.find(target);
    if (itrTarget == fParamLimits.end()) {
        std::map<std::string, std::map<int, std::map<int, std::pair<double, double> > > > tmp;
        std::map<int, std::map<int, std::pair<double, double> > > tmp2;
        std::map<int, std::pair<double, double> > tmp3;
        tmp3.insert({index, std::make_pair(min, max)});
        tmp2.insert({bin, tmp3});
        tmp.insert({region, tmp2});
        fParamLimits.insert({target, tmp});
    } else {
        auto itr = itrTarget->second.find(region);
        if (itr == itrTarget->second.end()) {
            std::map<int, std::pair<double,double> > tmp;
            tmp.insert({index, std::make_pair(min,max)});
            std::map<int, std::map<int, std::pair<double, double> > > tmp2;
            tmp2.insert({bin, tmp});
            itrTarget->second.insert(std::make_pair(region, tmp2));
        }  else {
            auto itr2 = itr->second.find(bin);
            if (itr2 == itr->second.end()) {
                std::map<int, std::pair<double, double> > tmp;
                tmp.insert({index, std::make_pair(min, max)});
                itr->second.insert(std::make_pair(bin, tmp));
            } else {
                itr2->second.insert(std::make_pair(index, std::make_pair(min, max)));
            }
        }
    }
}

bool Morphing::HasStartValue(const std::string& target, const std::string& region, const int bin, const int index) const {
    auto itrTarget = fStartValues.find(target);
    if (itrTarget == fStartValues.end()) return false;

    auto itrReg = itrTarget->second.find(region);

    if (itrReg == itrTarget->second.end()) return false;

    auto itr = itrReg->second.find(bin);
    if (itr == itrReg->second.end()) return false;

    auto itr2  = itr->second.find(index);
    return itr2 != itr->second.end();
}

bool Morphing::HasParamLimits(const std::string& target, const std::string& region, const int bin, const int index) const {

    auto itrTarget = fParamLimits.find(target);
    if (itrTarget == fParamLimits.end()) return false;

    auto itrReg = itrTarget->second.find(region);

    if (itrReg == itrTarget->second.end()) return false;

    auto itr = itrReg->second.find(bin);
    if (itr == itrReg->second.end()) return false;

    auto itr2 = itr->second.find(index);
    return itr2 != itr->second.end();
}

double Morphing::ParamStartValue(const std::string& target, const std::string& region, const int bin, const int index) const {
    if (!HasStartValue(target, region, bin, index)) {
        LOG(WARNING) << "Region: " << region << ", target: " << target << " bin: " << bin << " param index: " << index << " not found, returning -999\n";
        return -999;
    }

    return fStartValues.find(target)->second.find(region)->second.find(bin)->second.find(index)->second;
}

std::pair<double, double> Morphing::ParamLimits(const std::string& target, const std::string& region, const int bin, const int index) const {
    if (!HasParamLimits(target, region, bin, index)) {
        LOG(WARNING) << "Region: " << region << ", target: " << target << " bin: " << bin << " param index: " << index << " not found, returning -999\n";
        return {-999,-999};
    }

    return fParamLimits.find(target)->second.find(region)->second.find(bin)->second.find(index)->second;
}


void Morphing::AddExpression(const std::string& param, const std::pair<std::string, std::string>& expression) {
    auto itr = fExpressions.find(param);
    if (itr != fExpressions.end()) {
        LOG(WARNING) << "Parameter " << param << " is already registered, will not overwrite it!\n";
        return;
    }

    fExpressions.insert(std::make_pair(param, expression));
}
