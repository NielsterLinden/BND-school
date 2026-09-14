#include "TRExFitter/TruthSample.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Unfolding.h"

#include "TH1.h"

TruthSample::TruthSample(const std::string& name) :
    fUnfoldingName(""),
    fName(name),
    fTitle(""),
    fLineStyle(0),
    fLineColor(1),
    fLineWidth(2),
    fUseForPlotting(true),
    fTruthDistributionPath(""),
    fTruthDistributionFile(""),
    fTruthDistributionName("")
{
}

std::unique_ptr<TH1> TruthSample::GetHisto(const Unfolding* unfolding) const {

    const std::string path = fTruthDistributionPath == "" ? unfolding->fTruthDistributionPath : fTruthDistributionPath;
    const std::string file = fTruthDistributionFile == "" ? unfolding->fTruthDistributionFile : fTruthDistributionFile;
    const std::string name = fTruthDistributionName == "" ? unfolding->fTruthDistributionName : fTruthDistributionName;

    const std::string full = path+"/"+file+".root/"+name;
    
    return Common::HistFromFile(full);
}
