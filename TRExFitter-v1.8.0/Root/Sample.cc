// Class include
#include "TRExFitter/Sample.h"

// Framework includes
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"

// ROOT includes
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TChain.h"


// -------------------------------------------------------------------------------------------------
// Sample

//__________________________________________________________________________________
//
Sample::Sample(const std::string& name, Sample::SampleType type) :
    fName(name),
    fUseFriend(false),
    fType(type),
    fFitName(name),
    fTitle(name),
    fTexTitle(""),
    fGroup(""),
    fFillColor(10),
    fLineColor(kBlack),
    fNormalizedByTheory(true),
    fIgnoreSelection("FALSE"),
    fIgnoreWeight("FALSE"),
    fUseMCStat(true),
    fUseSystematics(true),
    fDivideBy(""),
    fMultiplyBy(""),
    fNormToSample(""),
    fSmooth(false),
    fSelection("1"),
    fMCweight("1"),
    fGridDID(""),
    fSeparateGammas(false),
    fMCstatScale(1.),
    fCorrelateGammasWithSample(""),
    fSystFromSample(""),
    fIsFolded(false),
    fUseGaussianShapeSysConstraint(false),
    fEFTSMReference("NONE"),
    fEFTParam(""),
    fEFTTitle(""),
    fEFTValue(""),
    fUnfoldingName(""),
    fUnfoldingSampleName(""),
    fNoPruning(false),
    fTemplateTarget("") {
}

//__________________________________________________________________________________
//
Sample::~Sample(){
}

//__________________________________________________________________________________
// cosmetics
void Sample::SetTitle(const std::string& title){
    fTitle = title;
}

//__________________________________________________________________________________
//
void Sample::SetFillColor(int color){
    fFillColor = color;
}

//__________________________________________________________________________________
//
void Sample::SetLineColor(int color){
    fLineColor = color;
}

//__________________________________________________________________________________
//
void Sample::NormalizedByTheory(const bool norm ){
    fNormalizedByTheory = norm;
}

//__________________________________________________________________________________
// read from ntuples
void Sample::SetMCweight(const std::string& weight){
    fMCweight = weight;
}

//__________________________________________________________________________________
//
void Sample::SetSelection(const std::string& selection){
    fSelection = selection;
}

//__________________________________________________________________________________
//
void Sample::AddNtuplePath(const std::string& path){
    fNtuplePaths.push_back(path);
}

//__________________________________________________________________________________
//
void Sample::AddNtupleFile(const std::string& file){
    fNtupleFiles.push_back(file);
}

//__________________________________________________________________________________
//
void Sample::AddNtupleName(const std::string& name){
    fNtupleNames.push_back(name);
}

//__________________________________________________________________________________
// read from histograms
void Sample::AddHistoPath(const std::string& path){
    fHistoPaths.push_back(path);
}

//__________________________________________________________________________________
//
void Sample::AddHistoFile(const std::string& file){
    fHistoFiles.push_back(file);
}

//__________________________________________________________________________________
//
void Sample::AddHistoName(const std::string& name){
    fHistoNames.push_back(name);
}

//__________________________________________________________________________________
// norm factors and systs
void Sample::AddNormFactor(std::shared_ptr<NormFactor> normFactor){
    fNormFactors.emplace_back(normFactor);
}

//__________________________________________________________________________________
//
void Sample::AddShapeFactor(std::shared_ptr<ShapeFactor> shapeFactor){
    fShapeFactors.emplace_back(shapeFactor);
}

//__________________________________________________________________________________
//
void Sample::AddSystematic(std::shared_ptr<Systematic> syst){
    fSystematics.emplace_back(syst);
}

//__________________________________________________________________________________
//
bool Sample::HasSystematic(const std::string& name) const{
    for(const auto& isyst : fSystematics) {
        if(isyst->fName==name) return true;
    }
    return false;
}

//__________________________________________________________________________________
//
bool Sample::HasNuisanceParameter(const std::string& name) const{
    for(const auto& isyst : fSystematics) {
        if(isyst->fNuisanceParameter==name) return true;
    }
    return false;
}

//__________________________________________________________________________________
//
bool Sample::HasNormFactor(const std::string& name) const{
    for(const auto& inorm : fNormFactors) {
        if(inorm->fName==name) return true;
    }
    return false;
}

//__________________________________________________________________________________
//
std::shared_ptr<NormFactor> Sample::AddNormFactor(const std::string& name,double nominal,double min,double max,bool isConst){
    fNormFactors.emplace_back(new NormFactor(name,nominal,min,max,isConst));
    return fNormFactors.back();
}

//__________________________________________________________________________________
//
std::shared_ptr<ShapeFactor> Sample::AddShapeFactor(const std::string& name,double nominal,double min,double max,bool isConst){
    fShapeFactors.emplace_back(new ShapeFactor(name,nominal,min,max,isConst));
    return fShapeFactors.back();
}

//__________________________________________________________________________________
//
std::shared_ptr<Systematic> Sample::AddSystematic(const std::string& name,int type,double up,double down){
    fSystematics.emplace_back(new Systematic(name,type,up,down));
    return fSystematics.back();
}
