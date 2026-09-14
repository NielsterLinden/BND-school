// Class include
#include "TRExFitter/SystematicHist.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/Systematic.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"

// -------------------------------------------------------------------------------------------------
// SystematicHist

//_____________________________________________________________________________
//
SystematicHist::SystematicHist(const std::string& name) :
    fName(name),
    fSystematic(nullptr),
    fIsShapeAccDecorr(false),
    fIsRegionDecorr(false),
    fIsSampleDecorr(false),
    fIsOverall(false),
    fHasShape(false),
    fSmoothType(0),
    fSymmetrisationType(HistoTools::NOSYMMETRIZATION),
    fShapePruned(false),
    fNormPruned(false),
    fBadShape(false),
    fBadNorm(false),
    fHistUp(nullptr),
    fHistUp_orig(nullptr),
    fHistUp_preSmooth(nullptr),
    fHistShapeUp(nullptr),
    fNormUp(0),
    fFileNameUp(""),
    fHistoNameUp(""),
    fFileNameShapeUp(""),
    fHistoNameShapeUp(""),
    fHistDown(nullptr),
    fHistDown_orig(nullptr),
    fHistDown_preSmooth(nullptr),
    fHistShapeDown(nullptr),
    fNormDown(0),
    fFileNameDown(""),
    fHistoNameDown(""),
    fFileNameShapeDown(""),
    fHistoNameShapeDown(""),
    fScaleUp(1.),
    fScaleDown(1.),
    fMagnitude(-1.) {
}

//_____________________________________________________________________________
//
SystematicHist::~SystematicHist(){
}

//_____________________________________________________________________________
//
void SystematicHist::WriteToFile(const std::vector<int>& droppedBins,
                                 const std::vector<double>& scales,
                                 const Common::FolderStructure& structure,
                                 const bool useRegularBinning,
                                 TFile* f,
                                 bool reWriteOrig) const{
    if (!f) {
        Common::WriteHistToFile(fHistUp.get(),fFileNameUp, structure);
        Common::WriteHistToFile(fHistDown.get(),fFileNameDown, structure);
        if (reWriteOrig) {
            Common::WriteHistToFile(fHistUp_orig.get(),fFileNameUp, structure);
            Common::WriteHistToFile(fHistDown_orig.get(),fFileNameDown, structure);
        }
        if (fHasShape) {
            Common::WriteHistToFile(fHistShapeUp.get(),fFileNameShapeUp, structure);
            Common::WriteHistToFile(fHistShapeDown.get(),fFileNameShapeDown, structure);
            std::unique_ptr<TH1D> shapeUp = HistoTools::TransformHistogramBinning(fHistShapeUp.get(), droppedBins, scales, useRegularBinning);
            std::unique_ptr<TH1D> shapeDown = HistoTools::TransformHistogramBinning(fHistShapeDown.get(), droppedBins, scales, useRegularBinning);
            Common::WriteHistToFile(shapeUp.get(),fFileNameShapeUp, structure);
            Common::WriteHistToFile(shapeDown.get(),fFileNameShapeDown, structure);
        }
        if (fSystematic->fType==Systematic::SHAPE) {
            std::unique_ptr<TH1D> up = HistoTools::TransformHistogramBinning(fHistUp.get(), droppedBins, scales, useRegularBinning);
            std::unique_ptr<TH1D> down = HistoTools::TransformHistogramBinning(fHistDown.get(), droppedBins, scales, useRegularBinning);
            Common::WriteHistToFile(up.get(),fFileNameUp, structure);
            Common::WriteHistToFile(down.get(),fFileNameDown, structure);
        }
    } else {
        Common::WriteHistToFile(fHistUp.get(),f, structure);
        Common::WriteHistToFile(fHistDown.get(),f, structure);
        if(reWriteOrig) {
            Common::WriteHistToFile(fHistUp_orig.get(),f, structure);
            Common::WriteHistToFile(fHistDown_orig.get(),f, structure);
        }
        if (fHasShape) {
            Common::WriteHistToFile(fHistShapeUp.get(),f, structure);
            Common::WriteHistToFile(fHistShapeDown.get(),f, structure);
            std::unique_ptr<TH1D> shapeUp = HistoTools::TransformHistogramBinning(fHistShapeUp.get(), droppedBins, scales, useRegularBinning);
            std::unique_ptr<TH1D> shapeDown = HistoTools::TransformHistogramBinning(fHistShapeDown.get(), droppedBins, scales, useRegularBinning);
            Common::WriteHistToFile(shapeUp.get(),f, structure);
            Common::WriteHistToFile(shapeDown.get(),f, structure);
        }
        if (fSystematic->fType==Systematic::SHAPE) {
            std::unique_ptr<TH1D> up = HistoTools::TransformHistogramBinning(fHistUp.get(), droppedBins, scales, useRegularBinning);
            std::unique_ptr<TH1D> down = HistoTools::TransformHistogramBinning(fHistDown.get(), droppedBins, scales, useRegularBinning);
            Common::WriteHistToFile(up.get(),f, structure);
            Common::WriteHistToFile(down.get(),f, structure);
        }
    }
}

//_____________________________________________________________________________
//
void SystematicHist::ReadFromFile(){
    fHistUp      = Common::HistFromFile(fFileNameUp,fHistoNameUp);
    fHistUp_orig = Common::HistFromFile(fFileNameUp,fHistoNameUp+"_orig");
    if(fHistUp_orig==nullptr) fHistUp_orig = std::unique_ptr<TH1>(static_cast<TH1D*>(fHistUp->Clone()));
    fHistShapeUp = Common::HistFromFile(fFileNameShapeUp,fHistoNameShapeUp);
    fHistDown      = Common::HistFromFile(fFileNameDown,fHistoNameDown);
    fHistDown_orig = Common::HistFromFile(fFileNameDown,fHistoNameDown+"_orig");
    if(fHistDown_orig==nullptr) fHistDown_orig = std::unique_ptr<TH1>(static_cast<TH1*>(fHistDown->Clone()));
    fHistShapeDown = Common::HistFromFile(fFileNameShapeDown,fHistoNameShapeDown);
}

//_____________________________________________________________________________
//
bool SystematicHist::IsShape() const{
    if(fHistUp!=nullptr || fHistDown!=nullptr) return true;
    return false;
}

//_____________________________________________________________________________
//
void SystematicHist::Print() const{
    std::string temp = "        Systematic: " + fName;
    if(fHistShapeUp==nullptr && fHistShapeDown==nullptr && fHistUp==nullptr && fHistDown==nullptr) temp += Form("\toverall (%.3f,%.3f)",fNormUp,fNormDown);
    LOG(INFO) << temp << "\n";
}

//_____________________________________________________________________________
//
void SystematicHist::Divide(TH1 *h){
    fHistUp->Divide(h);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Divide(h);
    fHistDown->Divide(h);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Divide(h);
}

//_____________________________________________________________________________
//
void SystematicHist::Divide(SystematicHist *syh){
    fHistUp->Divide(syh->fHistUp.get());
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Divide(  syh->fHistShapeUp.get());
    fHistDown->Divide(syh->fHistDown.get());
    if(fHistShapeDown!=nullptr) fHistShapeDown->Divide(syh->fHistShapeDown.get());
}

//_____________________________________________________________________________
//
void SystematicHist::Multiply(TH1 *h){
    fHistUp->Multiply(h);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Multiply(h);
    fHistDown->Multiply(h);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Multiply(h);
}

//_____________________________________________________________________________
//
void SystematicHist::Multiply(SystematicHist *syh){
    fHistUp->Multiply(syh->fHistUp.get());
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Multiply(  syh->fHistShapeUp.get());
    fHistDown->Multiply(syh->fHistDown.get());
    if(fHistShapeDown!=nullptr) fHistShapeDown->Multiply(syh->fHistShapeDown.get());
}

//_____________________________________________________________________________
//
void SystematicHist::Add(TH1 *h,double scale){
    fHistUp->Add(h,scale);
    if(fHistShapeUp!=nullptr)   fHistShapeUp->Add(h,scale);
    fHistDown->Add(h,scale);
    if(fHistShapeDown!=nullptr) fHistShapeDown->Add(h,scale);
}

//_____________________________________________________________________________
//
void SystematicHist::Add(SystematicHist *syh,double scale){
    fHistUp->Add(syh->fHistUp.get(),scale);
    if(fHistShapeUp!=nullptr)   {
        if (syh->fHistShapeUp != nullptr) fHistShapeUp->Add(syh->fHistShapeUp.get(),scale);
        else if (syh->fHistUp != nullptr){
          // the other syst is overall, while this is not: get by hand its dummy shape syst
          std::unique_ptr<TH1> htemp (static_cast<TH1*>(syh->fHistUp->Clone("hDummyShapeUp")));
          htemp->Scale(1.0/(1.0+syh->fNormUp));
          fHistShapeUp->Add(htemp.get(),scale);
        }
    }
    fHistDown->Add(syh->fHistDown.get(),scale);
    if(fHistShapeDown!=nullptr)   {
        if (syh->fHistShapeDown != nullptr) fHistShapeDown->Add(syh->fHistShapeDown.get(),scale);
        else if (syh->fHistDown != nullptr){// the other syst is overall, while this is not
          std::unique_ptr<TH1> htemp (static_cast<TH1*>(syh->fHistDown->Clone("hDummyShapeDown")));
          htemp->Scale(1.0/(1.0+syh->fNormDown));
          fHistShapeDown->Add(htemp.get(),scale);
        }
    }
}
