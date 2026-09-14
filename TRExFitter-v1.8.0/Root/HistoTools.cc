/*

 HistoTools
 =========

 Contains all functions needed to handle properly the histograms:
    -> Variable Binning handling for fit/limits
    -> Symmetrisation of systematics
    -> Smoothing of systematics

 Call of the functions:
    -> #include "TRExFitter/HistoTools.C"
    -> Call of the function with HistoTools:: (Helps readability)

 */

// Class include
#include "TRExFitter/HistoTools.h"

// framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/SystematicHist.h"

#include "SystematicSmoothingTool/SmoothSystematics/SmoothHist.h"

// ROOT includes
#include "TH1.h"

// c++ includes
#include <algorithm>
#include <iostream>
#include <memory>

using namespace std;

//_________________________________________________________________________
//
std::unique_ptr<TH1D> HistoTools::TransformHistogramBinning(const TH1* originalHist,
                                                            const std::vector<int>& droppedBins,
                                                            const std::vector<double>& scale,
                                                            const bool useRegularBinning) {

    const int bins = originalHist->GetNbinsX() - droppedBins.size();

    std::unique_ptr<TH1D> hFinal(nullptr);
    if (useRegularBinning) {
        hFinal = std::make_unique<TH1D>(originalHist->GetName()+static_cast<TString>("_regBin"),originalHist->GetTitle(),bins,0,1);
    } else {
        if (droppedBins.empty()) {
            hFinal.reset(static_cast<TH1D*>(originalHist->Clone()));
            hFinal->SetName(originalHist->GetName()+static_cast<TString>("_regBin"));

        } else {
            hFinal = std::make_unique<TH1D>(originalHist->GetName()+static_cast<TString>("_dropBin"),originalHist->GetTitle(),bins,0,1);
        }
    }

    hFinal->SetDirectory(nullptr);
    int iBinNew = 1;
    for(int iBin = 1; iBin <= originalHist->GetNbinsX(); ++iBin){
        if (std::find(droppedBins.begin(), droppedBins.end(), iBin) != droppedBins.end()) continue;
        int currentBin = (useRegularBinning || !droppedBins.empty()) ? iBinNew : iBin;
        if (originalHist->GetBinContent(iBin) < 0) {
            hFinal->SetBinContent(currentBin,1e-6);
            hFinal->SetBinError(currentBin,1e-6);
        } else {
            hFinal->SetBinContent(currentBin,originalHist->GetBinContent(iBin));
            hFinal->SetBinError(currentBin,originalHist->GetBinError(iBin));

        }
        iBinNew++;
    }

    if (!scale.empty()) {
        if (static_cast<int>(scale.size()) != hFinal->GetNbinsX()) {
            LOG(WARNING) << "Scales are not empty, but the size does not match the bin size!\n";
        } else {
            for (int ibin = 1; ibin <= hFinal->GetNbinsX(); ++ibin) {
                hFinal->SetBinContent(ibin, hFinal->GetBinContent(ibin) * scale.at(ibin-1));
                hFinal->SetBinError  (ibin, hFinal->GetBinError(ibin)   * scale.at(ibin-1));
            }
        }
    }
    return hFinal;
}

//_________________________________________________________________________
//
void HistoTools::ManageHistograms(const int smoothingLevel,
                                  const SymmetrizationType& symType,
                                  const TH1* hNom,
                                  const TH1* originUp,
                                  const TH1* originDown,
                                  std::unique_ptr<TH1> &modifiedUp,
                                  std::unique_ptr<TH1> &modifiedDown,
                                  double scaleUp,
                                  double scaleDown,
                                  const SmoothOption &smoothOpt) {
    //
    // Only function called directly to handle operations on the histograms (symmetrisation and smoothing)
    //

    auto CreateNew = [&](bool isUp) {
        std::unique_ptr<TH1> result(nullptr);
        if (isUp) {
            if (modifiedUp) {
                result.reset(static_cast<TH1*>(modifiedUp->Clone()));
            } else {
                result.reset(static_cast<TH1*>(originUp->Clone()));
            }
        } else {
            if (modifiedDown) {
                result.reset(static_cast<TH1*>(modifiedDown->Clone()));
            } else {
                result.reset(static_cast<TH1*>(originDown->Clone()));
            }
        }
        result->SetDirectory(nullptr);
        return result;
    };

    // if one-sided & symmetrization asked, do smoothing first and symmetrization after
    if( symType == SymmetrizationType::SYMMETRIZEONESIDED ){
        SmoothHistograms(smoothingLevel, hNom, originUp, originDown, modifiedUp, modifiedDown, smoothOpt);
        // Here it gets tricky. We just allocated new memory to the pointer passes as reference.
        // This is our new "original" and we also need to delete the memeory
        std::unique_ptr<TH1> newUp   = CreateNew(true);
        std::unique_ptr<TH1> newDown = CreateNew(false);
        SymmetrizeHistograms(symType, hNom, newUp.get(), newDown.get(), modifiedUp, modifiedDown, scaleUp, scaleDown);
        if (!modifiedUp || !modifiedDown) {
            LOG(ERROR) << "Something went wring with the smoothing!\n";
            exit(EXIT_FAILURE);
        }
    }
    // otherwise, first symmetrization and then smoothing
    else{
        SymmetrizeHistograms(symType, hNom, originUp, originDown, modifiedUp, modifiedDown, scaleUp, scaleDown);
        std::unique_ptr<TH1> newUp   = CreateNew(true);
        std::unique_ptr<TH1> newDown = CreateNew(false);
        SmoothHistograms(smoothingLevel, hNom, newUp.get(), newDown.get(), modifiedUp, modifiedDown, smoothOpt);
        if (!modifiedUp) {
            modifiedUp.reset(static_cast<TH1*>(newUp->Clone()));
        }
        if (!modifiedDown) {
            modifiedDown.reset(static_cast<TH1*>(newDown->Clone()));
        }
    }
}

//_________________________________________________________________________
//
void HistoTools::SymmetrizeHistograms(const SymmetrizationType& symType,
                                      const TH1* hNom,
                                      const TH1* originUp,
                                      const TH1* originDown,
                                      std::unique_ptr<TH1> &modifiedUp,
                                      std::unique_ptr<TH1> &modifiedDown,
                                      double scaleUp,
                                      double scaleDown) {
    //##################################################
    //
    // FIRST STEP: SYMMETRISATION
    //
    //##################################################

    //Just to be sure, set sumw2() on histograms
    if(hNom == nullptr) {
        LOG(WARNING) << "Nominal Histogram is nullptr.\n";
        LOG(WARNING) << "Will not symmetrize.\n";
        return;
    }

    // Want to do the scaling before the Symmetrisation
    // to avoid improperly scaling bins that are limited at -100%

    std::unique_ptr<TH1> localOriginUp  (static_cast<TH1*>(originUp   -> Clone()));
    std::unique_ptr<TH1> localOriginDown(static_cast<TH1*>(originDown -> Clone()));
    std::unique_ptr<TH1> localHNom      (static_cast<TH1*>(hNom       -> Clone()));

    std::string nameUp, nameDown;

    if(localOriginUp != nullptr) {
        Scale(localOriginUp.get(),   localHNom.get(), scaleUp);
        nameUp = localOriginUp->GetName();
    }
    if(localOriginDown != nullptr) {
        Scale(localOriginDown.get(), localHNom.get(), scaleDown);
        nameDown = localOriginDown->GetName();
    }

    if( symType == SymmetrizationType::SYMMETRIZEONESIDED ) {
        bool isUp = true; //is the provided uncertainty the up or down variation (based on yield)
        if     (localOriginUp==nullptr && localOriginDown!=nullptr) isUp = false;
        else if(localOriginUp!=nullptr && localOriginDown==nullptr) isUp = true;
        else if(localOriginUp==nullptr && localOriginDown==nullptr){
            LOG(WARNING) << "Both up and down variations are empty.\n";
            LOG(WARNING) << "Will not symmetrize.\n";
            return;
        }
        // if both are non-empty, check the differences with the nominal
        else{
            const double separationUp = Separation(localHNom.get(),localOriginUp.get());
            const double separationDown = Separation(localHNom.get(),localOriginDown.get());
            if( separationUp > separationDown ) isUp = true;
            if( separationUp < separationDown ) isUp = false;
        }

        std::unique_ptr<TH1D> temp = nullptr;

        if(isUp){
            temp = std::unique_ptr<TH1D>(SymmetrizeOneSided(localHNom.get(), localOriginUp.get(), isUp));
            modifiedUp.reset(static_cast<TH1*>(localOriginUp -> Clone()));
            modifiedDown.reset(static_cast<TH1*>(temp -> Clone()));
        } else {
            temp = std::unique_ptr<TH1D>(SymmetrizeOneSided(localHNom.get(), localOriginDown.get(), isUp));
            modifiedUp.reset(static_cast<TH1*>(temp -> Clone()));
            modifiedDown.reset(static_cast<TH1*>(localOriginDown -> Clone()));
        }
    } else if ( symType == SymmetrizationType::SYMMETRIZETWOSIDED ) {
        modifiedUp = SymmetrizeTwoSided(localOriginUp.get(), localOriginDown.get(), localHNom.get());
        modifiedDown = InvertShift(modifiedUp.get(),localHNom.get());
    } else if ( symType == SymmetrizationType::SYMMETRIZEABSMEAN ) {
        modifiedUp = SymmetrizeAbsMean(localOriginUp.get(), localOriginDown.get(), localHNom.get());
        modifiedDown = InvertShift(modifiedUp.get(),localHNom.get());
    } else if ( symType == SymmetrizationType::SYMMETRIZEMAXIMUM ) {
        modifiedUp = SymmetrizeMaximum(localOriginUp.get(), localOriginDown.get(), localHNom.get());
        modifiedDown = InvertShift(modifiedUp.get(),localHNom.get());
    } else if ( symType == SymmetrizationType::SYMMETRIZEMAXIMUMKEEPSIGN ) {
        modifiedUp = SymmetrizeMaximumKeepSign(localOriginUp.get(), localOriginDown.get(), localHNom.get());
        modifiedDown = InvertShift(modifiedUp.get(),localHNom.get());
    } else {
        modifiedUp.reset(localOriginUp.release());
        modifiedDown.reset(localOriginDown.release());
    }

    modifiedDown -> SetName(nameDown.c_str());
    modifiedUp   -> SetName(nameUp.c_str());
}

//_________________________________________________________________________
//
void HistoTools::SmoothHistograms(int smoothingLevel,
                                  const TH1* hNom,
                                  const TH1* originUp,
                                  const TH1* originDown,
                                  std::unique_ptr<TH1> &modifiedUp,
                                  std::unique_ptr<TH1> &modifiedDown,
                                  const SmoothOption &smoothOpt) {
    //##################################################
    //
    // SECOND STEP: SMOOTHING
    //
    //##################################################

    // No smoothing
    if (smoothingLevel == 0){
        return;
    }

    std::unique_ptr<TH1D> nom_tmp(static_cast<TH1D*>(hNom->Clone()));
    std::unique_ptr<TH1D> up_tmp(static_cast<TH1D*>(originUp->Clone()));
    std::unique_ptr<TH1D> down_tmp(static_cast<TH1D*>(originDown->Clone()));

    // Initialize common smoothing tool
    SmoothHist smoothTool;

    if(hNom->GetNbinsX()==1){
        std::string temp =  hNom->GetName();
        LOG(DEBUG) << "Skipping smoothing for systematics on \"" << temp << "\" since just 1 bin.\n";
        return;
    }
    if (smoothOpt == TTBARRESONANCE) {
        if( ( smoothingLevel >= SMOOTHDEPENDENT ) && ( smoothingLevel < SMOOTHINDEPENDENT ) ){
            modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothTtresDependent")->Clone()));
            modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothTtresDependent")->Clone()));
        } else if( ( smoothingLevel >= SMOOTHINDEPENDENT ) && (smoothingLevel < UNKNOWN) ){
            modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothTtresIndependent")->Clone()));
            modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothTtresIndependent")->Clone()));
        } else {
            LOG(WARNING) <<  "Unknown smoothing level!\n";
            return;
        }
    } else if (smoothOpt == MAXVARIATION){
        smoothTool.setTRExTolerance(0.08); // This was also default before
        const int lvl = smoothingLevel / 10;
        smoothTool.setTRExNbins(lvl);
        modifiedUp   .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothTRExDefault")->Clone()));
        modifiedDown .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothTRExDefault")->Clone()));
    } else if (smoothOpt == COMMONTOOLSMOOTHMONOTONIC){
        modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothRebinMonotonic")->Clone()));
        modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothRebinMonotonic")->Clone()));
    } else if (smoothOpt == COMMONTOOLSMOOTHPARABOLIC){
        modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothRebinParabolic")->Clone()));
        modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothRebinParabolic")->Clone()));
    }else if (smoothOpt == TCHANNEL){
        modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothTchannel")->Clone()));
        modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothTchannel")->Clone()));
    } else if (smoothOpt == KERNELRATIOUNIFORM){
        modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothRatioUniformKernel")->Clone()));
        modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothRatioUniformKernel")->Clone()));
    } else if (smoothOpt == KERNELDELTAGAUSS){
        modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothDeltaGaussKernel")->Clone()));
        modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothDeltaGaussKernel")->Clone()));
    } else if (smoothOpt == KERNELRATIOGAUSS){
        modifiedUp      .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), up_tmp.get(),   "smoothRatioGaussKernel")->Clone()));
        modifiedDown    .reset(static_cast<TH1*>(smoothTool.Smooth(nom_tmp.get(), down_tmp.get(), "smoothRatioGaussKernel")->Clone()));
    } else {
        LOG(WARNING) << "Unknown smoothing option. Please check the config file.\n";
        LOG(WARNING) << "No smoothing will be applied.\n";
        return;
    }
}

//_________________________________________________________________________
//
std::unique_ptr<TH1D> HistoTools::SymmetrizeOneSided(const TH1* const h_nominal, const TH1* const h_syst, bool &isUp ){

    const double& yield_nominal     = h_nominal->Integral();
    const double& yield_syst        = h_syst->Integral();

    //Convention: one sided systematic leading to larger yield: "up" variation
    if(yield_syst>yield_nominal) isUp = true;

    if(Separation(h_nominal,h_syst)<1e-05){
        const std::string& temp1 =  h_nominal -> GetName();
        const std::string& temp2 =  h_syst -> GetName();
        LOG(WARNING) << "The two histograms for one-sided symmetrisation are the same (difference is < 1e-5)\n";
        LOG(WARNING) << "Will not symmetrize\n";
        LOG(WARNING) << "      --> Nominal : " << temp1 << "\n";
        LOG(WARNING) << "      --> Syst    : " << temp2 << "\n";
    }

    return InvertShift(h_syst,h_nominal);
}

//_________________________________________________________________________
//
std::unique_ptr<TH1D> HistoTools::InvertShift(const TH1* const h_syst, const TH1* const h_nominal){

    //Sanity check
    if(!h_syst || !h_nominal){
        return nullptr;
    }

    //Compute the symmetric histogram
    std::unique_ptr<TH1D> result(static_cast<TH1D*>(h_nominal->Clone()));
    result -> SetDirectory(0);
    result -> Add(h_syst,-1);
    result -> Add(h_nominal,1);
    for( int iBin = 1; iBin <= result -> GetNbinsX(); ++iBin ){
        result -> SetBinError(iBin, h_syst->GetBinError(iBin));
    }

    //Another sanity check: search for negative bins
    for( int iBin = 1; iBin <= result -> GetNbinsX(); ++iBin ){
        const double& content = result -> GetBinContent(iBin);
        if(content < 0){
            result -> SetBinContent(iBin, 0.);
        }
    }
    return result;
}

//_________________________________________________________________________
//
double HistoTools::Separation(const TH1* const h1, const TH1* const h2){
    double sep = 0.;
    for(int i_bin=1;i_bin<=h1->GetNbinsX();i_bin++){
        sep += std::abs( h1->GetBinContent(i_bin) - h2->GetBinContent(i_bin) );
    }
    return sep;
}

//_________________________________________________________________________
//
std::unique_ptr<TH1D> HistoTools::SymmetrizeTwoSided(const TH1* const var1, const TH1* const var2, const TH1* const hnom) {

    bool isLarge = false;
    if (std::abs(var1->Integral()/hnom->Integral()-1) > 0.005) isLarge = true;
    if (std::abs(var2->Integral()/hnom->Integral()-1) > 0.005) isLarge = true;
    //
    // Symmetrize a variation that is already two sided to smooth out possible fluctuations
    //
    //Nominal
    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(hnom->Clone()));

    //Up variation
    std::unique_ptr<TH1D> tmp1 (static_cast<TH1D*>(var1->Clone()));
    tmp1->Divide(nom.get());
    if(!tmp1->GetSumw2())tmp1->Sumw2();

    //Down variation
    std::unique_ptr<TH1D> tmp2(static_cast<TH1D*> (var2->Clone()));
    tmp2->Divide(nom.get());
    if(!tmp2->GetSumw2())tmp2->Sumw2();

    if (isLarge){
        CheckSameShift(var1, var2, hnom, tmp1.get(), tmp2.get());
    }

    //Flat (content = 0) histogram to substract
    std::unique_ptr<TH1D> unit (static_cast<TH1D*> (nom->Clone()));
    if(!unit->GetSumw2())unit->Sumw2();
    for (int bin=1; bin<= unit->GetNbinsX(); bin++){
        unit->SetBinContent(bin,1);
        unit->SetBinError(bin,0.0);
    }

    //Computes Var/Nom - 1
    tmp1->Add(unit.get(),-1);
    tmp2->Add(unit.get(),-1);

    //Computes Corrected = (DeltaUp-DeltaDown)/2 + 1
    tmp1->Add(tmp2.get(),-1);
    tmp1->Scale(0.5);
    tmp1->Add(unit.get());

    //Computed the new histogram Corrected*Nominal
    tmp1->Multiply(nom.get());

    //Protection against negative bin
    for (int bin=1; bin<= unit->GetNbinsX(); bin++){
        const double& content = tmp1->GetBinContent(bin);
        if(content<0){
            tmp1->SetBinContent(bin,0.);
        }
    }

    return tmp1;
}

//_________________________________________________________________________
//
std::unique_ptr<TH1D> HistoTools::SymmetrizeAbsMean(const TH1* const var1, const TH1* const var2, const TH1* const hnom) {
    //Nominal
    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(hnom->Clone()));

    //Up variation
    std::unique_ptr<TH1D> tmp1 (static_cast<TH1D*>(var1->Clone()));
    tmp1->Divide(nom.get());
    if(!tmp1->GetSumw2())tmp1->Sumw2();

    //Down variation
    std::unique_ptr<TH1D> tmp2(static_cast<TH1D*> (var2->Clone()));
    tmp2->Divide(nom.get());
    if(!tmp2->GetSumw2())tmp2->Sumw2();

    //Flat (content = 0) histogram to substract
    std::unique_ptr<TH1D> unit (static_cast<TH1D*> (nom->Clone()));
    if(!unit->GetSumw2())unit->Sumw2();
    for (int bin=1; bin<= unit->GetNbinsX(); bin++){
        unit->SetBinContent(bin,1);
        unit->SetBinError(bin,0.0);
    }

    //Computes Var/Nom - 1
    tmp1->Add(unit.get(),-1);
    tmp2->Add(unit.get(),-1);

    // change DeltaUp to abs(DeltaUp) and DeltaDown to abs(DeltaDown)
    for (int ibin = 1; ibin <= tmp1->GetNbinsX(); ++ibin){
        tmp1->SetBinContent(ibin, std::abs(tmp1->GetBinContent(ibin)));
        tmp2->SetBinContent(ibin, std::abs(tmp2->GetBinContent(ibin)));
    }

    //Computes Corrected = (DeltaUp + DeltaDown)/2 + 1
    tmp1->Add(tmp2.get());
    tmp1->Scale(0.5);
    tmp1->Add(unit.get());

    //Computed the new histogram Corrected*Nominal
    tmp1->Multiply(nom.get());

    //Protection against negative bin
    for (int bin=1; bin<= unit->GetNbinsX(); bin++){
        const double& content = tmp1->GetBinContent(bin);
        if(content<0){
            tmp1->SetBinContent(bin,0.);
        }
    }

    return tmp1;
}

//_________________________________________________________________________
//
std::unique_ptr<TH1D> HistoTools::SymmetrizeMaximum(const TH1* const var1, const TH1* const var2, const TH1* const hnom) {
    //Nominal
    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(hnom->Clone()));

    //Up variation
    std::unique_ptr<TH1D> tmp1 (static_cast<TH1D*>(var1->Clone()));
    tmp1->Divide(nom.get());
    if(!tmp1->GetSumw2())tmp1->Sumw2();

    //Down variation
    std::unique_ptr<TH1D> tmp2(static_cast<TH1D*> (var2->Clone()));
    tmp2->Divide(nom.get());
    if(!tmp2->GetSumw2())tmp2->Sumw2();

    // check which variation is larger and use that one
    for (int ibin = 1; ibin <= tmp1->GetNbinsX(); ++ibin){
        if (std::abs(tmp1->GetBinContent(ibin) - 1) < std::abs(tmp2->GetBinContent(ibin) - 1)){
            tmp1->SetBinContent(ibin, tmp2->GetBinContent(ibin));
        }
    }

    //Computed the new histogram Corrected*Nominal
    tmp1->Multiply(nom.get());

    //Protection against negative bin
    for (int bin=1; bin<= tmp1->GetNbinsX(); bin++){
        const double& content = tmp1->GetBinContent(bin);
        if(content<0){
            tmp1->SetBinContent(bin,0.);
        }
    }

    return tmp1;
}

//_________________________________________________________________________
//
std::unique_ptr<TH1D> HistoTools::SymmetrizeMaximumKeepSign(const TH1* const var1, const TH1* const var2, const TH1* const hnom) {
    //Nominal
    std::unique_ptr<TH1D> nom(static_cast<TH1D*>(hnom->Clone()));

    //Up variation
    std::unique_ptr<TH1D> tmp1 (static_cast<TH1D*>(var1->Clone()));
    tmp1->Divide(nom.get());
    if(!tmp1->GetSumw2())tmp1->Sumw2();

    //Down variation
    std::unique_ptr<TH1D> tmp2(static_cast<TH1D*> (var2->Clone()));
    tmp2->Divide(nom.get());
    if(!tmp2->GetSumw2())tmp2->Sumw2();

    // if up varaition is larger, use it, if not, use nominal - down
    // but it is in ratios so we comapre wrt "1"
    for (int ibin = 1; ibin <= tmp1->GetNbinsX(); ++ibin){
        double variation = tmp1->GetBinContent(ibin);
        if (std::abs(variation - 1) < std::abs(tmp2->GetBinContent(ibin) - 1)) {
            // 1 - (down - 1)
            variation = 2 - tmp2->GetBinContent(ibin);
        }
        tmp1->SetBinContent(ibin, variation);
    }

    //Computed the new histogram Corrected*Nominal
    tmp1->Multiply(nom.get());

    //Protection against negative bin
    for (int bin=1; bin<= tmp1->GetNbinsX(); bin++){
        const double& content = tmp1->GetBinContent(bin);
        if(content<0){
            tmp1->SetBinContent(bin,0.);
        }
    }

    return tmp1;
}
//_________________________________________________________________________
//
void HistoTools::Scale(TH1* h_syst, TH1* h_nominal, double factor){

    //Sanity check
    if(!h_syst || !h_nominal) return;

    //Just to be sure, set sumw2() on histograms
    if(!h_nominal->GetSumw2())h_nominal -> Sumw2();
    if(!h_syst->GetSumw2())h_syst -> Sumw2();

    // clone nominal histogram
    std::unique_ptr<TH1> h_nominal_tmp(static_cast<TH1*>(h_nominal->Clone()));
    // set errors to zero
    for(int i_bin=0;i_bin<=h_nominal_tmp->GetNbinsX()+1;i_bin++){
        h_nominal_tmp->SetBinError(i_bin,0.);
    }

    //scale difference to nominal
    h_syst -> Add(h_nominal_tmp.get(),-1);
    h_syst -> Scale(factor);
    h_syst -> Add(h_nominal_tmp.get(),1);

    //Another sanity check: search for negative bins
    for( int iBin = 1; iBin <= h_syst -> GetNbinsX(); ++iBin ){
        double content = h_syst -> GetBinContent(iBin);
        if(content < 0){
            h_syst -> SetBinContent(iBin, 0.);
        }
    }
}

//_________________________________________________________________________
//
bool HistoTools::CheckHistograms(TH1* nom, SystematicHist* sh, bool checkNullContent, bool causeCrash){

    //======================================================
    // This code performs sanity checks on the provided
    // histograms.
    //======================================================

    //
    // Two kinds of detection
    //    -> ERRORS: fundamental problems (non sense of fit, missing inputs, ...)
    //    -> WARNINGS: potential issues (possibly affecting the fits or the plots)
    //

    bool isGood = true;

    //
    // 1) Checks that the histograms exist
    //
    if(nom == nullptr){
        if (causeCrash){
            LOG(ERROR) << "The nominal histogram doesn't seem to exist!\n";
            exit(EXIT_FAILURE);
        } else {
            LOG(WARNING) << "The nominal histogram doesn't seem to exist!\n";
        }
        return false;
    }
    if(sh != nullptr){
        if(sh->fHistUp == nullptr){
            if (causeCrash){
                LOG(ERROR) << "The up variation histogram doesn't seem to exist!\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "The up variation histogram doesn't seem to exist!\n";
            }
            return false;
        }
        if(sh->fHistDown == nullptr){
            if (causeCrash) {
                LOG(ERROR) << "The down variation histogram doesn't seem to exist!\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "The down variation histogram doesn't seem to exist!\n";
            }
            return false;
        }
    }
    else return false;

    //
    // 2) Checks the binning is the same for the three histograms
    //
    const int NbinsNom  = nom->GetNbinsX();
    const int NbinsUp   = sh->fHistUp->GetNbinsX();
    const int NbinsDown = sh->fHistDown->GetNbinsX();
    if((NbinsNom != NbinsUp) || (NbinsNom != NbinsDown)){
        if (causeCrash) {
            LOG(ERROR) << "The number of bins is found inconsistent! Please check!\n";
            exit(EXIT_FAILURE);
        } else {
            LOG(WARNING) << "The number of bins is found inconsistent! Please check!\n";
        }
        return false;
    }

    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double lowEdgeNom   = nom->GetBinLowEdge(iBin);
        double lowEdgeUp    = sh->fHistUp->GetBinLowEdge(iBin);
        double lowEdgeDown  = sh->fHistDown->GetBinLowEdge(iBin);

        if( std::abs(lowEdgeNom-lowEdgeUp)>1e-05 || std::abs(lowEdgeNom-lowEdgeDown)>1e-05 || std::abs(lowEdgeDown-lowEdgeUp)>1e-05 ){
            if (causeCrash) {
                LOG(ERROR) << "The bin low edges are not consistent! Please check!\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "The bin low edges are not consistent! Please check!\n";
            }
            return false;
        }

        double binWidthNom   = nom->GetBinWidth(iBin);
        double binWidthUp    = sh->fHistUp->GetBinWidth(iBin);
        double binWidthDown  = sh->fHistDown->GetBinWidth(iBin);

        if( std::abs(binWidthNom-binWidthUp)>1e-05 || std::abs(binWidthNom-binWidthDown)>1e-05 || std::abs(binWidthDown-binWidthUp)>1e-05 ){
            if (causeCrash) {
                LOG(ERROR) << "The bin widths are not consistent! Please check!\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "The bin widths are not consistent! Please check!\n";
            }
            return false;
        }
    }

    //
    // 3) Checks the absence on bins with 0 content (for nominal)
    //
    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double content = nom->GetBinContent(iBin);
        if( ( checkNullContent && content<=0 ) || ( !checkNullContent && content<0 ) ){
            std::string temp = nom->GetName();
            isGood = false;
            if (causeCrash) {
                LOG(ERROR) << "In histo \"" << temp << "\", bin " << iBin << " has 0 content ! Please check\n";
                LOG(ERROR) << "Nominal: " << content << "\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "In histo \"" << temp << "\", bin " << iBin << " has 0 content ! Please check\n";
                LOG(WARNING) << "Nominal: " << content << "\n";
                //Corrects the nominal
                LOG(WARNING) << "Will set the bin content to 1e-06 pm 1e-07! Please check!\n";
                nom -> SetBinContent(iBin,1e-06);
                nom -> SetBinError(iBin, 1e-07);
            }
        }

        //Now, for those histograms, checks if a systematics also has 0 content to
        //avoid 100% down systematics
        if (std::abs(nom->GetBinContent(iBin) - 1e-06) < 1e-10) {
            //This bin has most likely been changed to the default non-zero value
            sh->fHistUp->SetBinContent(iBin,1e-06);
            sh->fHistDown->SetBinContent(iBin,1e-06);
        }

    }//loop over the bins


    //
    // 4) Check the presence of abnormal bins (content-based)
    //
    for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
        double contentNom   = nom->GetBinContent(iBin);
        double contentUp    = sh->fHistUp->GetBinContent(iBin);
        double contentDown  = sh->fHistDown->GetBinContent(iBin);

        //
        // 4.a) Checks the presence of negative bins for systematic
        //

        if(contentUp<0){
            std::string temp_string = sh->fHistUp->GetName();
            isGood = false;
            if(causeCrash){
                LOG(ERROR) << "In histo \"" << temp_string << "\", bin " << iBin << " has negative content! Please check!\n";
                LOG(ERROR) << "  Nominal: " << contentNom << "\n";
                LOG(ERROR) << "  Up: " << contentUp << "\n";
                LOG(ERROR) << "  Down: " << contentDown << "\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "In histo \"" << temp_string << "\", bin " << iBin << " has negative content! Please check!\n";
                LOG(WARNING) << "  Nominal: " << contentNom << "\n";
                LOG(WARNING) << "  Up: " << contentUp << "\n";
                LOG(WARNING) << "  Down: " << contentDown << "\n";
                LOG(WARNING) << "  => Setting Up to 1e-06";
                sh->fHistUp->SetBinContent(iBin,1e-06);
            }
        }
        if(contentDown<0){
            std::string temp_string = sh->fHistDown->GetName();
            isGood = false;
            if(causeCrash) {
                LOG(ERROR) << "In histo \"" << temp_string << "\", bin " << iBin << " has negative content! Please check!\n";
                LOG(ERROR) << "  Nominal: " << contentNom << "\n";
                LOG(ERROR) << "  Up: " << contentUp << "\n";
                LOG(ERROR) << "  Down: " << contentDown << "\n";
                exit(EXIT_FAILURE);
            }
            else{
                LOG(WARNING) << "In histo \"" << temp_string << "\", bin " << iBin << " has negative content! Please check!\n";
                LOG(WARNING) << "  Nominal: " << contentNom << "\n";
                LOG(WARNING) << "  Up: " << contentUp << "\n";
                LOG(WARNING) << "  Down: " << contentDown << "\n";
                LOG(WARNING) << "  => Setting Down to 1e-06";
                sh->fHistDown->SetBinContent(iBin,1e-06);
            }
        }

        //
        // 4.b) Checks that the systematics are not crazy (too large ratio, nan returned, ...)
        //
        contentNom   = nom->GetBinContent(iBin);
        contentUp    = sh->fHistUp->GetBinContent(iBin);
        contentDown  = sh->fHistDown->GetBinContent(iBin);
        //
        double ratioUp   = 0.;
        double ratioDown = 0.;
        if(contentNom != 0 ){
            ratioUp   = contentUp  /contentNom;
            ratioDown = contentDown/contentNom;
        } else if( std::abs(contentUp)>0 || std::abs(contentDown)>0 ) {
            std::string temp_string = sh->fHistUp->GetName();
            LOG(WARNING) << "In histo " << temp_string << ", bin " << iBin << " has null nominal content but systematics are >0! Hope you know that\n";
            continue;
        }

        // * first check is for nan
        // * second check is for negative values (normally impossible after step 3 and 4.a)
        // * third check is for abnormal systematic/nominal values (ratio > 100)
        if( (ratioUp!=ratioUp) || (ratioUp < 0) || (std::abs(ratioUp-1.) >= 100) ){
            std::string temp_string = sh->fHistUp->GetName();
            isGood = false;
            if(causeCrash) {
                LOG(ERROR) << "In histo \"" << temp_string << "\", bin " << iBin << " has weird content! Please check!\n";
                LOG(ERROR) << "  Nominal: " << contentNom << "\n";
                LOG(ERROR) << "  Up: " << contentUp << "\n";
                LOG(ERROR) << "  Down: " << contentDown << "\n";
                exit(EXIT_FAILURE);
            }
            LOG(WARNING) << "In histo \"" << temp_string << "\", bin " << iBin << " has weird content! Please check!\n";
            LOG(WARNING) << "  Nominal: " << contentNom << "\n";
            LOG(WARNING) << "  Up: " << contentUp << "\n";
            LOG(WARNING) << "  Down: " << contentDown << "\n";
            // try to fix it, if not aborting
            if(ratioUp!=ratioUp) sh->fHistUp->SetBinContent(iBin,contentNom);
            else return isGood;
        }
        if( (ratioDown!=ratioDown) || (ratioDown < 0) || (std::abs(ratioDown-1.) >= 100) ){
            std::string temp_string = sh->fHistDown->GetName();
            isGood = false;
            if(causeCrash) {
                LOG(ERROR) << "In histo \"" << temp_string << "\", bin " << iBin << " has weird content! Please check!\n";
                LOG(ERROR) << "  Nominal: " << contentNom << "\n";
                LOG(ERROR) << "  Up: " << contentUp << "\n";
                LOG(ERROR) << "  Down: " << contentDown << "\n";
                exit(EXIT_FAILURE);
            }
            LOG(WARNING) << "In histo \"" << temp_string << "\", bin " << iBin << " has weird content! Please check!\n";
            LOG(WARNING) << "  Nominal: " << contentNom << "\n";
            LOG(WARNING) << "  Up: " << contentUp << "\n";
            LOG(WARNING) << "  Down: " << contentDown << "\n";
            // try to fix it, if not aborting
            if(ratioDown!=ratioDown) sh->fHistDown->SetBinContent(iBin,contentNom);
            else return isGood;
        }
    }
    return isGood;
}

//_________________________________________________________________________
//
void HistoTools::CheckSameShift(const TH1* const var1, const TH1* const var2, const TH1* const hnom,
    const TH1D* const tmp1, const TH1D* const tmp2){
    std::vector<int> binSameShift{};

    std::unique_ptr<TH1D> nom_div(static_cast<TH1D*>(hnom->Clone()));
    nom_div->Divide(hnom);

    // check if the variations are not in the same direction
    for (int ibin = 1; ibin <= hnom->GetNbinsX(); ibin++){
        const double& nom_error = nom_div->GetBinError(ibin);
        const double& up_content = tmp1->GetBinContent(ibin);
        const double& up_error = tmp1->GetBinError(ibin);
        const double& down_content = tmp2->GetBinContent(ibin);
        const double& down_error = tmp2->GetBinError(ibin);

        const double total_uncert_up = std::sqrt(up_error*up_error + nom_error*nom_error);
        const double total_uncert_down = std::sqrt(down_error*down_error + nom_error*nom_error);

        // check if the variations are larger than stat uncertainty
        if ((std::abs(up_content - 1) > total_uncert_up) && (std::abs(down_content - 1) > total_uncert_down)){
            // check if the variations have the same shift
            // have same sign = same shift
            if ((up_content - 1) * (down_content - 1) > 0){
                binSameShift.emplace_back(ibin);
            }
        }
    }
    if (binSameShift.size() > 0){
        const std::string& name_up = var1->GetName();
        const std::string& name_down = var2->GetName();
        std::string tmp = "";
        for (const auto& ibin : binSameShift){
            tmp+= std::to_string(ibin) + " ";
        }
        LOG(WARNING) << "Histograms " << name_up << " and " << name_down << " in bins " << tmp << " have the same shift for up and down variation\n";
        LOG(WARNING) << "You should check this\n";
    }
}

//_________________________________________________________________________
//
void HistoTools::ForceShape(TH1* syst, const TH1* nominal, const HistoTools::FORCESHAPETYPE type) {

    if (type == HistoTools::FORCESHAPETYPE::NOSHAPE) return;

    if (type == HistoTools::FORCESHAPETYPE::LINEAR) {
        ForceShapeLinear(syst, nominal);
    } else if (type == HistoTools::FORCESHAPETYPE::TRIANGULAR) {
        ForceShapeTriangular(syst, nominal);
    } else {
        LOG(WARNING) << "Unknown type for ForceShape, ignoring\n";
    }
}

//_________________________________________________________________________
//
void HistoTools::ForceShapeLinear(TH1* syst, const TH1* nominal) {
    const int nbins = syst->GetNbinsX();
    if (nbins < 2) return;

    for (int ibin = 1; ibin <= nbins; ++ibin) {
        const double correction = static_cast<double>(1. - 2.*(ibin-1.)/(nbins-1.));
        const double content = (syst->GetBinContent(ibin) - nominal->GetBinContent(ibin)) * correction + nominal->GetBinContent(ibin);
        syst->SetBinContent(ibin, content);
    }
}

//_________________________________________________________________________
//
void HistoTools::ForceShapeTriangular(TH1* syst, const TH1* nominal) {
    const int nbins = syst->GetNbinsX();
    const int nbinsHalf = nbins / 2;
    if (nbins < 3) return;

    for (int ibin = 1; ibin <= nbins; ++ibin) {
        double correction;
        if (nbins % 2 == 0){
            correction = ibin <= nbinsHalf ?
                static_cast<double>((ibin - 1.)/nbinsHalf) :
                static_cast<double>(1. - 1.0*(ibin - nbinsHalf)/nbinsHalf);
        } else {
            if (ibin <= nbinsHalf) {
                correction = static_cast<double>((ibin - 1.)/nbinsHalf);
            } else if (ibin == (nbinsHalf+1)) {
                correction = 1.0;
            } else {
                correction = static_cast<double>(1. - 1.0*(ibin - nbinsHalf)/nbinsHalf);
            }
        }

        const double content = (syst->GetBinContent(ibin) - nominal->GetBinContent(ibin)) * correction + nominal->GetBinContent(ibin);
        syst->SetBinContent(ibin, content);
    }
}

//_________________________________________________________________________
//
void HistoTools::ScaleBinByBin(TH1* h, const std::vector<double>& scales) {
    const int nBins = h->GetNbinsX();

    if (nBins != static_cast<int>(scales.size())) {
        LOG(ERROR) << "Number of the bins and the scales do not match\n";
        exit(EXIT_FAILURE);
    }

    for (int ibin = 1; ibin <= nBins; ++ibin) {
        h->SetBinContent(ibin, scales.at(ibin-1)*h->GetBinContent(ibin));
        h->SetBinError(ibin, scales.at(ibin-1)*h->GetBinError(ibin));
    }
}
