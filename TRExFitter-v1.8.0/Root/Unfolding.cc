#include "TRExFitter/Unfolding.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/NormFactor.h"

#include <algorithm>
#include <limits>

Unfolding::Unfolding() :
    fName("Unfolding"),
    fMatrixOrientation(FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS),
    fTruthDistributionPath(""),
    fTruthDistributionFile(""),
    fTruthDistributionName(""),
    fNumberUnfoldingTruthBins(0),
    fUnfoldingResultMin(0),
    fUnfoldingResultMax(2),
    fUnfoldingTitleX("X axis"),
    fUnfoldingTitleY("Y axis"),
    fUnfoldingMinRangeY(0.0001),
    fUnfoldingScaleRangeY(-1),
    fUnfoldingTitleOffsetY(1.6),
    fUnfoldingTitleOffsetX(1),
    fUnfoldingRatioYmin(0.5),
    fUnfoldingRatioYmax(1.5),
    fUnfoldingLogX(false),
    fUnfoldingLogY(false),
    fMigrationLogX(false),
    fMigrationLogY(false),
    fMigrationTitleOffsetX(1),
    fMigrationTitleOffsetY(1.5),
    fMigrationZmin(0),
    fMigrationZmax(1),
    fResponseZmin(0),
    fResponseZmax(1),
    fPlotSystematicMigrations(false),
    fMigrationText(false),
    fNominalTruthSample("SetMe"),
    fAlternativeAsimovTruthSample(""),
    fUnfoldingDivideByBinWidth(false),
    fUnfoldingDivideByLumi(-1),
    fUnfoldingChi2Type(Unfolding::UnfoldingChi2Type::EMPTY),
    fUnfoldNormXSec(false),
    fUnfoldNormXSecBinN(-1),
    fErrorBandColor(-1),
    fStatErrorBandColor(-1),
    fRatioLineColor(kRed),
    fRatioLineStyle(2),
    fRatioLineWidth(2),
    fLegendXposition(0.5),
    fLegendTextSize(-1),
    fBreakdownRange(std::numeric_limits<double>::max())
{
}

//___________________________________________________________
//
std::string Unfolding::GetTitleRegionX(const std::string& regionName) const {
    return this->GetTitleRegion(regionName, true);
}

//___________________________________________________________
//
std::string Unfolding::GetTitleRegionY(const std::string& regionName) const {
    return this->GetTitleRegion(regionName, false);
}

//___________________________________________________________
//
std::string Unfolding::GetTitleRegion(const std::string& regionName, const bool isX) const {
    const auto& axis = isX ? fMigrationTitleX : fMigrationTitleY;

    if (axis.empty()) return "NOT SET";

    if (axis.size() == 1 && axis.at(0).first == "all") {
        return axis.at(0).second;
    }

    auto itr = std::find_if(axis.begin(), axis.end(), [&regionName](const auto& element){return element.first == regionName;});
    if (itr == axis.end()) {
        return "NOT FOUND";
    }

    return itr->second;
}

//___________________________________________________________
//
bool Unfolding::HasReparametrisation(const std::vector<std::shared_ptr<NormFactor> >& normFactors) const {

    for (int ibin = 0; ibin < fNumberUnfoldingTruthBins; ++ibin) {
        const std::string name = fName + "_Bin_" + Common::IntToFixLenStr(ibin+1) + "_mu";
        auto itr = std::find_if(normFactors.begin(), normFactors.end(), [&name](const auto& element){return name == element->fName;});
        if (itr != normFactors.end()) {
            if (!(*itr)->fExpression.first.empty()) return true;
        }
    }

    return false;
}