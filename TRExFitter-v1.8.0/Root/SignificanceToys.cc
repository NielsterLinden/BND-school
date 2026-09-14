#include "TRExFitter/SignificanceToys.h"
#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"

#include "RooRandom.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"

#include "TCanvas.h"
#include "TFile.h"

SignificanceToys::SignificanceToys() :
    fNtoysSplusB(0),
    fNtoysB(1000),
    fToysSeed(1234),
    fPlot(true),
    fOutputPath(""),
    fSignificanceSuffix("")
{
}

void SignificanceToys::RunToys(RooAbsData* data,
                               RooStats::ModelConfig* mcSplusb,
                               RooStats::ModelConfig* mcB) const {

    if (!data || !mcSplusb || !mcB) {
        LOG(ERROR) << "Data or one of the ModelConfigs is nullptr!\n";
        return;
    }

    LOG(INFO) << "Running toys for significance (B-only toys: " << fNtoysB << ", S+B toys" << fNtoysSplusB << ")\n";

    mcSplusb->SetName("SplusBModel");

    RooRandom::randomGenerator()->SetSeed(fToysSeed);

    RooRealVar* poi = static_cast<RooRealVar*>(mcB->GetParametersOfInterest()->first());
    mcB->SetSnapshot(*poi);
    RooRealVar* poisb = static_cast<RooRealVar*>(mcSplusb->GetParametersOfInterest()->first());
    mcSplusb->SetSnapshot(*poisb);

    RooStats::FrequentistCalculator freqCalc(*data, *mcSplusb, *mcB);

    auto plr = std::make_unique<RooStats::ProfileLikelihoodTestStat>(*mcB->GetPdf());
    plr->SetOneSidedDiscovery(true);

    const RooArgSet* glbObs = mcSplusb->GetGlobalObservables();
    RooStats::ToyMCSampler *toymcs = static_cast<RooStats::ToyMCSampler*>(freqCalc.GetTestStatSampler());
    if (glbObs) {
        toymcs->SetGlobalObservables(*glbObs);
    }
    toymcs->SetTestStatistic(plr.get());

    if (!mcSplusb->GetPdf()->canBeExtended()) {
        toymcs->SetNEventsPerToy(1);
    }

    freqCalc.SetToys(fNtoysB, fNtoysSplusB);

    RooStats::HypoTestResult* r = freqCalc.GetHypoTest();
    r->SetPValueIsRightTail(true);
    r->SetBackgroundAsAlt(false);

    LOG(INFO) << "----------------------------------------------------\n";
    LOG(INFO) << "Results of the toys for significance\n";
    LOG(INFO) << "----------------------------------------------------\n";
    LOG(INFO) << " Significance " << r->Significance() << " +- " << r->SignificanceError() << "\n";
    LOG(INFO) << "----------------------------------------------------\n";

    if (fPlot) {
        LOG(INFO) << "Printing output plot to " << fOutputPath << "\n";
        TCanvas c{};
        auto plot = std::make_unique<RooStats::HypoTestPlot>(*r);
        plot->SetLogYaxis(true);
        plot->Draw();
        c.Draw();
        Common::SaveCanvasAs(c, fOutputPath+"/SignificanceToysResult");
    }

    this->RootOutput(r);
}

void SignificanceToys::RootOutput(RooStats::HypoTestResult* result) const {
    LOG(INFO) << "Storing the toys results in: " << fOutputPath << "/Toys" << fSignificanceSuffix << ".root\n";
    std::unique_ptr<TFile> file(TFile::Open((fOutputPath+"/Toys"+fSignificanceSuffix+".root").c_str(), "RECREATE"));
    if (!file) {
        LOG(ERROR) << "Cannot open file at: " << fOutputPath << "/Toys" << fSignificanceSuffix << ".root\n";
        return;
    }
    result->Write();

    file->Write();
    file->Close();
}
