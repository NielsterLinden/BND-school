// Class include
#include "TRExFitter/PruningUtil.h"
#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"

// C++ includes
#include <memory>

// -------------------------------------------------------------------------------------------------
// class PruningUtil

//__________________________________________________________________________________
//
PruningUtil::PruningUtil() :
    fStrategy(0),
    fShapeOption(PruningUtil::SHAPEOPTION::MAXBIN),
    fThresholdNorm(-1),
    fThresholdShape(-1),
    fThresholdIsLarge(-1),
    fRemoveSystOnEmptySample(false)
{
}

//__________________________________________________________________________________
//
void PruningUtil::SetStrategy(const int strat) {
    fStrategy = strat;
}

//__________________________________________________________________________________
//
void PruningUtil::SetShapeOption(const PruningUtil::SHAPEOPTION opt) {
    fShapeOption = opt;
}

//__________________________________________________________________________________
//
void PruningUtil::SetThresholdNorm(const double thres) {
    fThresholdNorm = thres;
}

//__________________________________________________________________________________
//
void PruningUtil::SetThresholdShape(const double thres) {
    fThresholdShape = thres;
}

//__________________________________________________________________________________
//
void PruningUtil::SetThresholdIsLarge(const double thres) {
    fThresholdIsLarge = thres;
}

//__________________________________________________________________________________
//
void PruningUtil::SetRemoveSystOnEmptySample(const bool flag) {
    fRemoveSystOnEmptySample = flag;
}

//__________________________________________________________________________________
//
std::tuple<bool, bool, bool ,bool> PruningUtil::CheckSystPruning(const TH1* const hUp,
                                                                 const TH1* const hDown,
                                                                 const TH1* const hNom,
                                                                 double& magnitude,
                                                                 const TH1* hTot) {

    if(fStrategy!=0 && hTot==nullptr){
        LOG(ERROR) << "hTot set to 0 while asking for relative pruning... Reverting to sample-by-sample pruning.\n";
        fStrategy = 0;
    }
    const TH1* hRef(nullptr);
    if(fStrategy==0) hRef = hNom;
    else hRef = hTot;

    // create shape-only syst variations
    std::unique_ptr<TH1> hShapeUp        = nullptr;
    if(hUp) hShapeUp     = std::unique_ptr<TH1>(static_cast<TH1*>(hUp  ->Clone(Form("%s_shape",hUp  ->GetName()))));
    if(hShapeUp) hShapeUp->Scale( Common::EffIntegral(hNom)/Common::EffIntegral(hShapeUp.get()) );
    std::unique_ptr<TH1> hShapeDown      = nullptr;
    if(hDown) hShapeDown = std::unique_ptr<TH1>(static_cast<TH1*>(hDown->Clone(Form("%s_shape",hDown->GetName()))));
    if(hShapeDown) hShapeDown->Scale( Common::EffIntegral(hNom)/Common::EffIntegral(hShapeDown.get()) );
    if (hShapeUp) hShapeUp->SetDirectory(nullptr);
    if (hShapeDown) hShapeDown->SetDirectory(nullptr);

    // get norm effects
    const double normUp   = std::abs((Common::EffIntegral(hUp  )-Common::EffIntegral(hNom))/Common::EffIntegral(hRef));
    const double normDown = std::abs((Common::EffIntegral(hDown)-Common::EffIntegral(hNom))/Common::EffIntegral(hRef));
    const double normNom = Common::EffIntegral(hNom);

    magnitude = std::max(normUp,normDown);

    double shapeMax{0};
    // check if systematic has no shape --> 1
    bool hasShape(true);
    if(fThresholdShape>=0) {
        if (fShapeOption == PruningUtil::SHAPEOPTION::MAXBIN) {
            hasShape = HasShapeRelative(hNom,hShapeUp.get(),hShapeDown.get(),hRef,fThresholdShape,shapeMax);
        }
        if (fShapeOption == PruningUtil::SHAPEOPTION::KSTEST) {
            hasShape = HasShapeKS(hNom,hShapeUp.get(),hShapeDown.get(),fThresholdShape);
        }
    }

    if (shapeMax > magnitude) {
        magnitude = shapeMax;
    }

    // check if systematic norm effect is under threshold
    bool hasNorm = true;
    if(fThresholdNorm>=0) hasNorm = ((normUp >= fThresholdNorm) || (normDown >= fThresholdNorm));

    // now check for crazy systematics
    bool hasGoodShape = true;
    bool hasGoodNorm = true;
    if(fThresholdIsLarge>=0) {
        if (fShapeOption == PruningUtil::SHAPEOPTION::MAXBIN) {
            hasGoodShape = !HasShapeRelative(hNom,hShapeUp.get(),hShapeDown.get(),hRef,fThresholdIsLarge, shapeMax);
        }
        if (fShapeOption == PruningUtil::SHAPEOPTION::KSTEST) {
            hasGoodShape = !HasShapeKS(hNom,hShapeUp.get(),hShapeDown.get(),fThresholdIsLarge);
        }
        hasGoodNorm = ((normUp <= fThresholdIsLarge) && (normDown <= fThresholdIsLarge));
    }

    if (fRemoveSystOnEmptySample) {
        if (normNom < 1e-4) {
            hasShape = false;
            hasNorm = false;
        }
    }

    return std::make_tuple(!hasNorm, !hasShape, !hasGoodNorm, !hasGoodShape);
}

//_________________________________________________________________________
//
bool PruningUtil::HasShapeRelative(const TH1* const hNom,
                                   const TH1* const hUp,
                                   const TH1* const hDown,
                                   const TH1* const combined,
                                   const double threshold,
                                   double& magnitude) const {
    if (!hNom || !hUp || !hDown || !combined) return false;

    if (hUp->GetNbinsX() == 1) return false;

    const double& integralUp = hUp->Integral();
    const double& integralDown = hDown->Integral();
    const double& integralCombined = combined->Integral();

    if (std::isnan(integralUp) || integralUp == 0) return false;
    if (std::isnan(integralDown) || integralDown == 0) return false;
    if (std::isnan(integralCombined) || integralCombined == 0) return false;

    bool hasShape = false;

    double max{-1};
    for (int ibin = 1; ibin <= hUp->GetNbinsX(); ++ibin){
        const double& nominal  = hNom->GetBinContent(ibin);
        if(nominal<0) continue;
        const double& comb     = combined->GetBinContent(ibin);
        const double& up       = hUp->GetBinContent(ibin);
        const double& down     = hDown->GetBinContent(ibin);
        const double& up_err   = std::abs((up-nominal)/comb);
        const double& down_err = std::abs((down-nominal)/comb);
        if(up_err>=threshold || down_err>=threshold){
            hasShape = true;
        }
        double tmp = std::max(up_err,down_err);
        if (tmp > max) {
            max = tmp;
        }
    }

    magnitude = max;

    return hasShape;
}

//_________________________________________________________________________
//
bool PruningUtil::HasShapeKS(const TH1* const hNom,
                             const TH1* const hUp,
                             const TH1* const hDown,
                             const double threshold) const {
    if (!hNom || !hUp || !hDown) return false;

    if (hUp->GetNbinsX() == 1) return false;

    if (std::abs(hUp->Integral()) < 1e-6) return false;
    if (std::abs(hDown->Integral()) < 1e-6) return false;

    // check if the histograms are not identical as KS test in the toys setup fails in this case
    bool isIdentical(true);
    for (int ibin = 1; ibin <= hNom->GetNbinsX(); ++ibin) {
        const double nom  = hNom->GetBinContent(ibin);
        const double down = hDown->GetBinContent(ibin);
        const double up   = hUp->GetBinContent(ibin);
        if (std::abs(down-nom) > 1e-6 || std::abs(up-nom) > 1e-6) {
            isIdentical = false;
            break;
        }
    }

    if (isIdentical) return false;

    const double probThreshold = 1 - threshold;

    // first try up histogram
    const double upProb   = hUp->KolmogorovTest(hNom, "X");
    if (upProb <= probThreshold) {
        return true;
    }

    // up is not significant, try down
    const double& downProb = hDown->KolmogorovTest(hNom, "X");

    if (downProb <= probThreshold) {
        return true;
    }

    return false;
}
