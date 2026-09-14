#include "TRExFitter/FitToys.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/LikelihoodScanManager.h"
#include "TRExFitter/Logger.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooMinimizer.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"

#include "StyleUtils/TRExLabels.h"
#include "StyleUtils/TRExUtils.h"

#include <iomanip>

FitToys::FitToys(const std::string& name, const std::string& suffix, const std::string& inputName) :
    fName(name),
    fSuffix(suffix),
    fInputName(inputName),
    fMinos({}),
    fIsRegularizedUnfolding(false),
    fNormFactors({}),
    fNToys(100),
    fRegularizationType(-1),
    fUseAutoDiff(false),
    fFitType(0),
    fHistoNbins(-1),
    fPOIs({}),
    fFitIsBlind(false),
    fSeed(123456),
    fStatOnlyFluctuation(false),
    fPlotLabel(""),
    fRunLHScanForAll(false)
{
}

void FitToys::RunToys(RooWorkspace* ws, FitResults* fr) {
    std::vector<std::pair<double,double> > binLimits;
    // get RooStats stuff
    RooStats::ModelConfig mc = *(static_cast<RooStats::ModelConfig*>(ws -> obj("ModelConfig")));
    RooSimultaneous simPdf = *(static_cast<RooSimultaneous*>((mc.GetPdf())));
    RooAbsPdf *pdf = mc.GetPdf();
    RooArgSet obsSet = *(mc.GetObservables());
    std::vector<RooRealVar*> nfs;
    for (const auto& nf: fNormFactors) {
        const std::string& iname = nf->fName;
        if (iname.find("morph_") != std::string::npos) continue;
        if(nf->fExpression.first!="") continue;
        if (iname.find("Expression") != std::string::npos) continue;
        nfs.emplace_back(static_cast<RooRealVar*>((& ws->allVars()[iname.c_str()])));
        binLimits.emplace_back(std::make_pair(nf->GetMin(), nf->GetMax()));
    }
    std::vector<RooRealVar*> nfs_and_nps_and_pois;
    if (mc.GetNuisanceParameters()) {
        for(auto ivar: *mc.GetNuisanceParameters()){
            RooRealVar* ivarReal = static_cast<RooRealVar*>(ivar);
            nfs_and_nps_and_pois.emplace_back(ivarReal);
        }
    }
    for(auto ivar: *mc.GetParametersOfInterest()){
        const std::string& iname = ivar->GetName();
        if (iname.find("morph_") != std::string::npos) continue;
        nfs_and_nps_and_pois.emplace_back(static_cast<RooRealVar*>(ivar));
    }

    //For the loop over NPs
    std::string varname{};

    //Create NLL only once
    std::unique_ptr<RooDataSet> dummy(pdf->generate(obsSet, RooFit::Extended()));

    // apply external constraints
    const RooArgSet* externalConstraints = FitUtils::GetExternalConstraints(ws, fNormFactors, fRegularizationType);

    const RooArgSet* glbObs = mc.GetGlobalObservables();
    const RooArgSet tmpDummy{};
    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (fUseAutoDiff) {
        nll = std::unique_ptr<RooAbsReal> (simPdf.createNLL(*dummy,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::Offset(1),
                                           NumCPU(1, RooFit::Hybrid),
                                           RooFit::ExternalConstraints(externalConstraints ? *externalConstraints : tmpDummy),
                                           RooFit::EvalBackend("codegen")));
    } else {
        nll = std::unique_ptr<RooAbsReal> (simPdf.createNLL(*dummy,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Offset(1),
                                           NumCPU(1, RooFit::Hybrid),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::ExternalConstraints(externalConstraints ? *externalConstraints : tmpDummy)));
    }

    auto isUnfolding = [&](std::size_t index) {
        if (fNormFactors.at(index)->fName.find("Bin_") != std::string::npos) {
            return true;
        }
        return false;
    };

    std::vector<Double_t> vals(nfs_and_nps_and_pois.size());
    std::vector<Double_t> val_errs(nfs_and_nps_and_pois.size());
    Double_t gcc;
    Int_t fitStatus;
    TTree *toys_tree = new TTree("toys","toys");
    for (size_t i=0; i<nfs_and_nps_and_pois.size(); ++i){
        toys_tree->Branch(nfs_and_nps_and_pois.at(i)->GetName(),&vals[i]);
        toys_tree->Branch((std::string(nfs_and_nps_and_pois.at(i)->GetName())+"_error").c_str(),&val_errs[i]);
    }
    if (fIsRegularizedUnfolding) {
        toys_tree->Branch("gcc",&gcc);
    }
    toys_tree->Branch("fit_status", &fitStatus);

    std::vector<TH1D> h_toys;
    std::vector<TH1D> h_pulls;
    for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
        h_toys.emplace_back(("h_toys_nf_"+std::to_string(inf)).c_str(),("h_toys_nf_"+std::to_string(inf)).c_str(),fHistoNbins,binLimits.at(inf).first,binLimits.at(inf).second);
        h_toys.back().SetDirectory(nullptr);
        if (fFitType == 3 && isUnfolding(inf)) {
            h_pulls.emplace_back(("h_pulls_nf_"+std::to_string(inf)).c_str(),("h_pulls_nf_"+std::to_string(inf)).c_str(),fHistoNbins,-3,3);
            h_pulls.back().SetDirectory(nullptr);
        }
    }

        // create non-const list of GlobalObservables
    RooArgSet toy_gobs;
    if (glbObs) {
        for(const auto np : *glbObs){
            toy_gobs.add(*np);
        }
    }

    // randomize GlobalObservables using ToyMCSampler
    RooStats::ProfileLikelihoodTestStat ts(*mc.GetPdf());
    RooStats::ToyMCSampler sampler(ts,fNToys);
    sampler.SetPdf(*mc.GetPdf());
    sampler.SetObservables(*mc.GetObservables());
    sampler.SetGlobalObservables(*mc.GetGlobalObservables());
    sampler.SetParametersForTestStat(*mc.GetParametersOfInterest());
    sampler.SetUseMultiGen(false);
    RooStats::ToyMCSampler::SetAlwaysUseMultiGen(false);

    RooArgSet nuisance;
    nuisance.add(FitUtils::GetNPsWithoutNFs(&mc, fNormFactors, fPOIs));
    RooArgSet* nullParams = static_cast<RooArgSet*>(nuisance.snapshot());
    ws->saveSnapshot("paramsToFitPE",nuisance);

    std::vector<double> gccVec;

    std::vector<std::string> toysLHscans;

    for(int i_toy = 0; i_toy < fNToys; ++i_toy) {
        ws->loadSnapshot("paramsToFitPE");

        auto FindNF = [this](const std::string& name) {
            for (auto&& inf : fNormFactors) {
                if (inf->fName == name) return inf;
            }
            LOG(ERROR) << "This should not happen\n";
            exit(EXIT_FAILURE);
        };

        // setting POIs and NFs to constant, not to allow it to fluctuate in toy creation
        for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
            nfs.at(inf)->setConstant(1);
            const std::string name = nfs.at(inf)->GetName();
            if (fFitIsBlind) {
                if (Common::FindInStringVector(fPOIs,name)>=0) {
                    auto itr = fFitPOIAsimov.find(name);
                    if (itr == fFitPOIAsimov.end()) {
                        auto&& nf = FindNF(name);
                        nfs.at(inf)->setVal(nf->GetNominal());
                    } else {
                        nfs.at(inf)->setVal(itr->second);
                    }
                } else {
                    auto&& nf = FindNF(name);
                    nfs.at(inf)->setVal(nf->GetNominal());
                }
            } else {
                const double postFitValue = fr->GetNuisParValue(name);
                LOG(DEBUG) << "Setting NF: " << name << " to: " << postFitValue << " for the toys generation\n";
                nfs.at(inf)->setVal(postFitValue);
            }
        }

        LOG(INFO) << "Generating toy n. " << (i_toy+1) << " out of " << fNToys << " toys\n";
        // set seed for toy dataset generation
        auto&& rnd = RooRandom::randomGenerator();
        rnd->SetSeed(fSeed+i_toy);
        // set GlobalObservables to randomized values and generate toy data
        std::unique_ptr<RooDataSet> toyData(nullptr);
        if (fStatOnlyFluctuation) {
            toyData = std::unique_ptr<RooDataSet>(pdf->generate(obsSet, RooFit::Extended()));
        } else {
            toyData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(sampler.GenerateToyData(*nullParams)));
        }

        // let the NF free float again
        for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
            if (fFitType==2 && Common::FindInStringVector(fPOIs,fNormFactors.at(inf)->fName)>=0) {
                nfs.at(inf)->setVal(0);
                nfs.at(inf)->setConstant(true);
            } else {
                auto&& nf = FindNF(nfs.at(inf)->GetName());
                const std::string name = nfs.at(inf)->GetName();
                if (nf->fConst) {
                    nfs.at(inf)->setConstant(true);
                } else {
                    nfs.at(inf)->setConstant(false);
                }
                nfs.at(inf)->setVal(nf->GetNominal());
                auto itr = fRandomPOIStartingValues.find(name);
                if (itr != fRandomPOIStartingValues.end()) {
                    const double initialValue = gRandom->Uniform(itr->second.first, itr->second.second);
                    nfs.at(inf)->setVal(initialValue);
                    LOG(DEBUG) << "Setting NormFactor " << name << " to a random initial value: " << initialValue << "\n";
                }
            }
        }

        // extract POI from fit result and fill histogram
        LOG(INFO) << "Fitting toy n. " << (i_toy+1) << "\n";
        //Set new dataset for NLL
        nll->setData(*toyData);
        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
        ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
        ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
        const double tol = ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //oAsymptoticCalculator enforces not less than 1 on this
        const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
        RooMinimizer m(*nll); // get MINUIT interface of fit
        m.optimizeConst(2);
        m.setPrintLevel(-1);
        m.setStrategy(1);
        m.setMinimizerType(minimType.Data());
        m.setEps(tol);
        m.migrad();
        std::unique_ptr<RooFitResult> r(m.save()); // save fit result
        fitStatus = r->status();

        std::size_t unfIndex(0);
        for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
            h_toys.at(inf).Fill(nfs.at(inf)->getVal());
            LOG(INFO) << "Toy n. " << (i_toy+1) << ", fitted value of NF: " << nfs.at(inf)->GetName() << ": " << nfs.at(inf)->getVal() << " +/- " << nfs.at(inf)->getError() << "\n";

            if (fFitType == 3 && isUnfolding(inf)) {
                const double value = nfs.at(inf)->getError() > 1e-6 ? (nfs.at(inf)->getVal() - 1.) / nfs.at(inf)->getError() : -9999;
                h_pulls.at(unfIndex).Fill(value);
                ++unfIndex;
            }
        }
        // fill pulls and constraints of all NPs, NFs, POIs into tree
        for (std::size_t inf = 0; inf < nfs_and_nps_and_pois.size(); ++inf) {
            vals.at(inf) = nfs_and_nps_and_pois.at(inf)->getVal();
            val_errs.at(inf) = nfs_and_nps_and_pois.at(inf)->getError();
        }
        if (fIsRegularizedUnfolding) {
            std::vector<std::string> poiNames;
            for (const auto& inf : fNormFactors) {
                if (inf->fName.find("Bin_") != std::string::npos) {
                    poiNames.emplace_back(inf->fName);
                }
            }
            gcc = Common::GCC(r.get(), poiNames);
            gccVec.emplace_back(gcc);
            LOG(INFO) << "Global Correlation Coefficient is " << gcc << "\n";
        }
        toys_tree->Fill();

        // run LH scan if requested
        for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
            const std::string name = nfs.at(inf)->GetName();
            auto itr = fLHScanCondition.find(name);
            if (itr == fLHScanCondition.end() && !fRunLHScanForAll) continue;
            if (fRunLHScanForAll) {
                auto itrPOI = std::find(fPOIs.begin(), fPOIs.end(), name);
                if (itrPOI == fPOIs.end()) continue;
            }

            const double min = itr->second.first;
            const double max = itr->second.second;
            std::string extension("");
            if (fRunLHScanForAll) {
                extension = "_toys_" + name + "_toyN_" + std::to_string(i_toy);
            } else {
                extension = "_toys_" + name + "_" + std::to_string(min) + "_" + std::to_string(max);
            }

            auto itrScan = std::find(toysLHscans.begin(), toysLHscans.end(), extension);
            bool runLHscan(false);
            if (fRunLHScanForAll) {
                runLHscan = true;
            } else {
                runLHscan = (itrScan == toysLHscans.end()) && (nfs.at(inf)->getVal() > min) && (nfs.at(inf)->getVal() < max);
            }

            if (runLHscan) {
                if (!fRunLHScanForAll) {
                    LOG(INFO) << "Will run LH scan for parameter: " << name << " that is between: " << min << " and " << max << "\n";
                }

                const double minValue = nfs.at(inf)->getMin();
                const double maxValue = nfs.at(inf)->getMax();

                this->LikelihoodScan(ws, name, toyData.get(), extension, minValue, maxValue);

                toysLHscans.emplace_back(extension);
            }
        }
    }

    if (fIsRegularizedUnfolding) {
        double gccAverage = std::accumulate(gccVec.begin(), gccVec.end(), 0.);
        gccAverage /= gccVec.size();
        double accum(0.);
        std::for_each(gccVec.begin(), gccVec.end(), [&accum,&gccAverage](const double d) {
            accum += (d - gccAverage) * (d - gccAverage);
        });
        const double stdev = std::sqrt(accum / (gccVec.size() - 1));
        LOG(INFO) << "Average Global Correlation Coefficient estimated from toys is " << gccAverage << " +- " << stdev << "\n";
    }

    // save individual pulls and constraints of all NPs into a root file
    std::unique_ptr<TFile> toys_out (TFile::Open((fName+"/Toys/Toys_NP_pulls"+fSuffix+"_Seed"+std::to_string(fSeed)+".root").c_str(), "RECREATE"));
    toys_tree->SetDirectory(toys_out.get());
    toys_out->Write();

    toys_out->Close();

    std::unique_ptr<TFile> out (TFile::Open((fName+"/Toys/Toys"+fSuffix+"_Seed"+std::to_string(fSeed)+".root").c_str(), "RECREATE"));
    std::vector<double> errors;
    for (std::size_t inf = 0; inf < nfs.size(); ++inf) {
        out->cd();
        TCanvas c("c","c",600,600);
        h_toys.at(inf).Draw("E");
        TF1 g("g","gaus",binLimits.at(inf).first,binLimits.at(inf).second);
        g.SetLineColor(kRed);
        h_toys.at(inf).Fit(&g,"RQ");
        g.Draw("same");
        h_toys.at(inf).GetXaxis()->SetTitle(nfs.at(inf)->GetName());
        h_toys.at(inf).GetYaxis()->SetTitle("Pseudo-experiements");
        myText(0.60,0.90,1,Form("Mean  = %.2f #pm %.2f",g.GetParameter(1),g.GetParError(1)));
        myText(0.60,0.85,1,Form("Sigma = %.2f #pm %.2f",g.GetParameter(2),g.GetParError(2)));
        errors.emplace_back(g.GetParameter(2));
        myText(0.60,0.80,1,Form("#chi^{2}/ndf = %.2f / %d",g.GetChisquare(),g.GetNDF()));
        Common::SaveCanvasAs(c, fName+"/Toys/ToysPlot_"+nfs.at(inf)->GetName()+"_Seed"+std::to_string(fSeed));

        out->cd();
        // Also create a ROOT file
        h_toys.at(inf).Write();
    }

    if (fFitType == 3) {
        DrawToyPullPlot(h_pulls, out.get());
    }

    for (auto& ihist : h_toys) {
        ihist.SetDirectory(nullptr);
    }
    for (auto& ihist : h_pulls) {
        ihist.SetDirectory(nullptr);
    }

    h_pulls.clear();
    h_toys.clear();

    out->Close();
}

void FitToys::DrawToyPullPlot(const std::vector<TH1D>& hist, TFile* out) const {
    std::vector<double> mean;
    std::vector<double> sigma;
    std::vector<double> meanError;
    std::vector<double> sigmaError;

    for (auto& ihist : hist) {
        mean.emplace_back(ihist.GetMean());
        sigma.emplace_back(ihist.GetRMS());
        meanError.emplace_back(ihist.GetMeanError());
        sigmaError.emplace_back(ihist.GetRMSError());

        out->cd();
        ihist.Write();
    }

    TH1D h_mean("","",hist.size(), 0, hist.size());
    TH1D h_sigma("","",hist.size(), 0, hist.size());

    for (std::size_t i = 0; i < hist.size(); ++i) {
        h_mean.SetBinContent(i+1, mean.at(i));
        h_mean.SetBinError(i+1, meanError.at(i));
        h_sigma.SetBinContent(i+1, sigma.at(i));
        h_sigma.SetBinError(i+1, sigmaError.at(i));

        h_mean.GetXaxis()->SetBinLabel(i+1, ("bin " + std::to_string(i+1)).c_str());
    }

    h_mean.GetYaxis()->SetRangeUser(-1, 3);
    h_mean.GetYaxis()->SetTitle("#frac{#mu_{fitted}-1}{#sigma_{fitted}}");
    h_mean.LabelsOption("v");
    h_mean.SetLineColor(kBlue);
    h_mean.SetMarkerColor(kBlue);
    h_sigma.SetLineColor(kRed);
    h_sigma.SetMarkerColor(kRed);

    TCanvas c("c","c",800,600);
    h_mean.Draw("PE");
    h_sigma.Draw("PE SAME");

    TLine l_mean(0, 0, hist.size(), 0);
    TLine l_sigma(0, 1, hist.size(), 1);

    l_mean.SetLineColor(kBlue);
    l_mean.SetLineWidth(2);
    l_mean.SetLineStyle(2);
    l_sigma.SetLineColor(kRed);
    l_sigma.SetLineWidth(2);
    l_sigma.SetLineStyle(2);

    l_mean.Draw("same");
    l_sigma.Draw("same");

    TLegend leg(0.7, 0.8, 0.9, 0.9);
    leg.SetFillStyle(0.);
    leg.SetBorderSize(0.);
    leg.AddEntry(&h_mean, "mean", "p");
    leg.AddEntry(&h_sigma, "sigma", "p");
    leg.Draw("same");

    if (fPlotLabel != "none") TRExLabel(0.3,0.85, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str());

    gPad->RedrawAxis();

    out->cd();
    h_mean.Write("pull_mean");
    h_sigma.Write("pull_sigma");

    Common::SaveCanvasAs(c, fName+"/Toys/PullPlot");
}

void FitToys::LikelihoodScan(RooWorkspace *ws,
                             const std::string& varName,
                             RooDataSet* data,
                             const std::string& toysExtension,
                             const double min,
                             const double max) {

    LOG(INFO) << "Running likelihood scan for the parameter = " << varName << "\n";

    const RooArgSet* externalConstraints = FitUtils::GetExternalConstraints(ws, fNormFactors, fRegularizationType);

    LikelihoodScanManager manager{};

    manager.SetScanParamsX(min, max, 50, -1);
    manager.SetNCPU(1);
    manager.SetOffSet(true);
    manager.SetUseNll(true);
    manager.SetExternalConstraints(externalConstraints);
    manager.SetToleranceScale(1.0);
    manager.SetUseAutoDiff(fUseAutoDiff);
    manager.SetNLLOffset("initial");

    FitUtils::FittingOptions options;
    options.noGammas = false;
    options.noSystematics = false;
    options.randomize = false;
    std::vector<std::string> constNPs;
    std::vector<double> constNPvalues;
    for (const auto& iparam : fFitFixedNPs) {
        constNPs.emplace_back(iparam.first);
        constNPvalues.emplace_back(iparam.second);
    }
    options.constNPs = std::move(constNPs);
    options.constNPvalues = std::move(constNPvalues);
    std::vector<std::string> npNames;
    std::vector<double> npValues;
    for(const auto& inf : fNormFactors) {
        if (Common::FindInStringVector(fPOIs,inf->fName)>=0) {
            continue;
        }
        npNames. emplace_back(inf->fName);
        npValues.emplace_back(inf->GetNominal());
    }
    options.initialNPs = std::move(npNames);
    options.initialNPvalues = std::move(npValues);
    options.constPOI = fFitType == 2;
    std::vector<std::pair<std::string, double> > poiVals;
    std::vector<std::string> constNFs;
    for (const auto& inf : fNormFactors) {
        poiVals.emplace_back(std::make_pair(inf->fName, inf->GetNominal()));
        if (inf->fConst) {
            constNFs.emplace_back(inf->fName);
        }
    }
    options.poiVals = std::move(poiVals);
    options.constNFs = std::move(constNFs);
    manager.SetFittingOptions(options);

    auto scanResult = manager.Run1DScan(ws, varName, data);

    const std::vector<double> x = scanResult.first;
    const std::vector<double> y = scanResult.second;

    if (x.empty() || y.empty() ) {
        LOG(WARNING) << "Skipping LHscan\n";
        return;
    }

    gSystem->mkdir((fName+"/LHoodPlots").c_str());

    TCanvas can("NLLscan");

    TGraph graph(x.size(), x.data(), y.data());
    graph.Draw("ALP");
    graph.GetXaxis()->SetRangeUser(x.at(0),x.back());

    // y axis
    graph.GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");
    if(TRExFitter::SYSTMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::SYSTMAP[varName].c_str());
    else if(TRExFitter::NPMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::NPMAP[varName].c_str());

    TString cname="";
    cname.Append("NLLscan_");
    cname.Append(varName);

    can.SetTitle(cname);
    can.SetName(cname);
    can.cd();

    TLatex tex{};
    tex.SetTextColor(kGray+2);

    TLine l1s(x.at(0),0.5,x.back(),0.5);
    l1s.SetLineStyle(kDashed);
    l1s.SetLineColor(kGray);
    l1s.SetLineWidth(2);
    if(graph.GetMaximum()>2){
        l1s.Draw();
        tex.DrawLatex(x.back(),0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
    }

    if(graph.GetMaximum()>2){
        TLine l2s(x.at(0),2,x.back(),2);
        l2s.SetLineStyle(kDashed);
        l2s.SetLineColor(kGray);
        l2s.SetLineWidth(2);
        l2s.Draw();
        tex.DrawLatex(x.back(),2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
    }

    if(graph.GetMaximum()>4.5){
        TLine l3s(x.at(0),4.5,x.back(),4.5);
        l3s.SetLineStyle(kDashed);
        l3s.SetLineColor(kGray);
        l3s.SetLineWidth(2);
        l3s.Draw();
        tex.DrawLatex(x.back(),4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
    }
    //
    TLine lv0(0,graph.GetMinimum(),0,graph.GetMaximum());
    lv0.Draw();
    //
    TLine lh0(x.at(0),0,x.back(),0);
    lh0.Draw();
    can.RedrawAxis();

    if (!fRunLHScanForAll) {
        Common::SaveCanvasAs(can, fName+"/LHoodPlots/NLLscan_"+varName + toysExtension +fSuffix);
    }

    // write it to a ROOT file as well
    std::unique_ptr<TFile> f(TFile::Open((fName+"/LHoodPlots/NLLscan_"+varName+fSuffix+ "_toys_curve.root").c_str(),"UPDATE"));
    if (!f) {
        LOG(WARNING) << "Cannot open ROOT file for likelihood scan!\n";
        return;
    }
    f->cd();
    graph.Write(("LHscan"+toysExtension).c_str(),TObject::kOverwrite);
    f->Close();

}
