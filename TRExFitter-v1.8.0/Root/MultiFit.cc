// Class include
#include "TRExFitter/MultiFit.h"

// Framework inncludes
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReader.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitToys.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/LikelihoodScanManager.h"
#include "TRExFitter/LimitToys.h"
#include "TRExFitter/LimitEvent.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/RankingManager.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SignificanceToys.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/UncertaintyPlotter.h"
#include "TRExFitter/Unfolding.h"
#include "TRExFitter/YamlConverter.h"

// Roofit includes
#include "RooDataSet.h"
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"

// HistFactory includes
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HypoTestInverterResult.h"

// Style stuff
#include "StyleUtils/TRExStyle.h"
#include "StyleUtils/TRExLabels.h"
#include "StyleUtils/TRExUtils.h"

// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"
#include "xRooFit/xRooNode.h"

#include <algorithm>
#include <filesystem>

using namespace std;
using namespace RooFit;
using namespace RooStats;

// -------------------------------------------------------------------------------------------------
// class MultiFit

//__________________________________________________________________________________
//
MultiFit::MultiFit(const string& name) :
    fCombine(false),
    fCompare(false),
    fStatOnly(false),
    fIncludeStatOnly(false),
    fCompareLimits(true),
    fComparePOI(true),
    fCompareEFT(true),
    fComparePulls(true),
    fPlotCombCorrMatrix(false),
    fName(name),
    fDir(""),
    fOutDir(""),
    fLabel(""),
    fShowObserved(false),
    fLimitTitle("95% CL limit on XXX"),
    fRankingOnly("all"),
    fGroupedImpactCategory("all"),
    fLimitMax(0),
    fUseRnd(false),
    fRndRange(0.1),
    fRndSeed(-999),
    fLumiLabel(""),
    fCmeLabel(""),
    fCombiLabel("Combined"),
    fConfig(std::unique_ptr<ConfigParser>(new ConfigParser())),
    fSaveSuf(""),
    fDataName("obsData"),
    fFitType(1), // 1: S+B, 2: B-only
    fFastFit(false),
    fFastFitForRanking(true),
    fNuisParListFile(""),
    fPlotSoverB(false),
    fSignalTitle("signal"),
    fFitResultsRootFile(""),
    fLimitsFile(""),
    fBonlySuffix(""),
    fShowSystForPOI(false),
    fLHscanMin(999999),
    fLHscanMax(-999999),
    fLHscanSteps(30),
    fLHscanStep(-1),
    fLHscanMinY(999999),
    fLHscanMaxY(-999999),
    fLHscanStepsY(30),
    fLHscanStepY(-1),
    fDoGroupedSystImpactTable(false),
    fLimitIsBlind(false),
    fSignalInjection(false),
    fSignalInjectionValue(0),
    fLimitParamName("parameter"),
    fLimitParamValue(0),
    fLimitOutputPrefixName("myLimit"),
    fLimitsConfidence(0.95),
    fLimitUseAutoDiff(false),
    fSignificanceIsBlind(false),
    fSignificanceDoInj(false),
    fSignificancePOIAsimov(0),
    fSignificanceParamName("parameter"),
    fSignificanceParamValue(0),
    fSignificanceOutputPrefixName("mySignificance"),
    fSignificanceUseAutoDiff(false),
    fShowTotalOnly(false),
    fHEPDataFormat(false),
    fUheppFormat(false),
    fFitStrategy(-1),
    fCPU(1),
    fBinnedLikelihood(false),
    fUsePOISinRanking(false),
    fUseHesseBeforeMigrad(false),
    fUseNllInLHscan(true),
    fLimitToysStepsSplusB(100),
    fLimitToysStepsB(100),
    fLimitToysScanSteps(21),
    fLimitToysScanMin(0.),
    fLimitToysScanMax(10.),
    fLimitToysSeed(1234),
    fLimitPlot(true),
    fLimitFile(true),
    fLimitType(TRExFit::LimitType::ASYMPTOTIC),
    fLimitFitStrategy(2),
    fLimitToysUsexRooFit(TRExFit::ToysUsexRooFit::FALSE),
    fSignificanceType(TRExFit::SignificanceType::ASYMPTOTIC),
    fSignificanceFitStrategy(2),
    fSignificanceToysStepsSplusB(0),
    fSignificanceToysStepsB(100),
    fSignificancePlot(true),
    fSignificanceToysSeed(1234),
    fSignificanceToysUsexRooFit(false),
    fSpeedUpFit(false),
    fNPCutOff(0.3),
    fPlotUnfoldedData(false),
    fUnfoldingShowStat(false),
    fUnfoldingShowUncertaintyBreakdown(false),
    fUnfoldingUncertaintyBreakdownTotal(true),
    fPlotGammas(false),
    fHasShapeFactorReparametrisation(false),
    fUseHesse(true),
    fMaximumNumberFCNcalls(-1),
    fToleranceScale(1.),
    fShiftGlobalObservablesInRanking(false),
    fUseAutoDiff(false),
    fNLLOffset("initial"),
    fFitToys(-1),
    fToysHistoNbins(50),
    fToysSeed(1234),
    fToysSetStatOnlyFluctuation(false),
    fToysLHScanForAll(false)
{
    fNPCategories.emplace_back("");
}

//__________________________________________________________________________________
//
void MultiFit::AddFitFromConfig(const std::string& configFile,
                                const std::string& opt,
                                const std::string& options,
                                const std::string& label,
                                const std::string& loadSuf,
                                const std::string& wsFile,
                                const bool useInFit,
                                const bool useInComparison){


    // check if the config is not already processed (but it might be intended, if comparing different fits from same config)
    if (std::find(fConfigPaths.begin(), fConfigPaths.end(), configFile) != fConfigPaths.end()){
        LOG(WARNING) << "Config " << configFile << " is added twice. Make sure you know what you are doing.\n";
    }

    fConfigPaths.emplace_back(configFile);

    // keep debug level
    const int debug = TRExFitter::DEBUGLEVEL;

    fFitList.emplace_back(std::make_unique<TRExFit>());

    fFitList.back()->fUseInFit = useInFit;

    fFitList.back()->fUseInComparison = useInComparison;

    // initialize config reader
    ConfigReader reader(fFitList.back().get());

    if (reader.ReadFullConfig(configFile,opt,options) != 0){
        LOG(ERROR) << "Failed to read the config file.\n";
        exit(EXIT_FAILURE);
    }

    fFitLabels.push_back(label);
    fFitSuffs.push_back(loadSuf);
    fWsFiles.push_back(wsFile);

    TRExFitter::DEBUGLEVEL = debug;
}

//__________________________________________________________________________________
//
std::unique_ptr<RooWorkspace> MultiFit::CombineWS() const{

    this->CheckSameRegionNames();

    gSystem->mkdir((fOutDir+"/RooStats").c_str());

    std::unique_ptr<RooStats::HistFactory::Measurement> measurement = nullptr;

    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        if (!fFitList.at(i_fit)->fUseInFit) continue;
        const std::string& fitName = fFitList[i_fit]->fInputName;
        const std::string& fitDir = fFitList[i_fit]->fName;
        LOG(DEBUG) << "Adding Fit: " << fitName << ", " << fFitLabels[i_fit] << ", " << fFitSuffs[i_fit] << fitDir << "\n";

        std::string fileName = fitDir + "/RooStats/" + fitName + "_combined_" + fitName + fFitSuffs[i_fit] + "_model.root";
        if(fWsFiles[i_fit]!="") fileName = fWsFiles[i_fit];
        LOG(DEBUG) << "Opening file " << fileName << "\n";
        TFile* rootFile = TFile::Open(fileName.c_str(),"read");
        if (!rootFile) {
            LOG(ERROR) << "Cannot open file: " << fileName << "\n";
            exit(EXIT_FAILURE);
        }
        LOG(DEBUG) << "Getting " << fitName << fFitSuffs[i_fit] << "\n";
        std::unique_ptr<RooStats::HistFactory::Measurement> meas(dynamic_cast<RooStats::HistFactory::Measurement*>(rootFile -> Get( (fitName+fFitSuffs[i_fit]).c_str())));
        if (!meas) {
            LOG(ERROR) << "Measurement is nullptr\n";
            exit(EXIT_FAILURE);
        }
        //
        // import measurement if not there yet
        if(!measurement){
            measurement = std::move(meas);
        } else {
            auto channels = meas->GetChannels();
            for (auto& channel : channels) {
                if (!measurement->HasChannel(channel.GetName())) measurement->AddChannel(channel);
            }
            // Add expressions
            auto functionobjects = meas->GetFunctionObjects();
            for(auto & f : functionobjects){
              measurement->AddFunctionObject(f);
            }
        }
    }

    for (const auto& ifit : fFitList) {
        for (const auto& inf : ifit->fNormFactors) {
            if (!inf->fConst) continue;
            LOG(INFO) << "Setting NF: " << inf->fName << " to constant\n";
            measurement->AddConstantParam(inf->fName);
        }
    }

    //
    // Create the HistoToWorkspaceFactoryFast object to perform safely the combination
    //
    if(!measurement){
        LOG(ERROR) << "The measurement object has not been retrieved! Please check.\n";
        exit(EXIT_FAILURE);
    }

    // propagate the POIs
    std::vector<std::string> pois;
    for (const auto& ifit : fFitList) {
        for (const auto& ipoi : ifit->fPOIs) {
            if (std::find(pois.begin(), pois.end(), ipoi) == pois.end()) {
                pois.emplace_back(ipoi);
                measurement->AddPOI(ipoi.c_str());
            }
        }
    }

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    RooStats::HistFactory::HistoToWorkspaceFactoryFast factory(*measurement);

    // Creating the combined model
    std::unique_ptr<RooWorkspace> ws(factory.MakeCombinedModel(*measurement));

    if (fBinnedLikelihood) {
        FitUtils::SetBinnedLikelihoodOptimisation(ws.get());
    }


    LOG(INFO) << "....................................\n";

    // Configure the workspace
    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement( "simPdf", ws.get(), *measurement );

    measurement->PrintXML((fOutDir + "/RooStats/").c_str());
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    return ws;
}

//__________________________________________________________________________________
//
void MultiFit::SaveCombinedWS() const{
    //
    // Creating the rootfile
    //
    std::unique_ptr<TFile> f(TFile::Open( (fOutDir+"/ws_combined"+fSaveSuf+".root").c_str(), "recreate"));
    if (!f) {
        LOG(ERROR) << "Cannot open file!\n";
        exit(EXIT_FAILURE);
    }
    //
    // Creating the workspace
    //
    std::unique_ptr<RooWorkspace> ws = CombineWS();

    if (fHasShapeFactorReparametrisation) {
        auto regions = this->GetFitRegions();
        std::vector<std::string> fitRegions;
        for (const auto& ireg : regions) {
            fitRegions.emplace_back(ireg->fName);
        }

        std::vector<std::shared_ptr<ShapeFactor> > shapeFactors;
        for (const auto& ifit : fFitList) {
            shapeFactors.insert(shapeFactors.end(), ifit->fShapeFactors.begin(), ifit->fShapeFactors.end());
        }
        FitUtils::ReparametrizeShapeFactors(ws.get(), shapeFactors, this->GetFitRegions(), fitRegions);
    }
    //
    // Save the workspace
    //
    f->cd();
    ws->Write("combWS");
    f->Close();
}

//__________________________________________________________________________________
//
std::map < std::string, double > MultiFit::FitCombinedWS(int fitType, const std::string& inputData, bool doLHscanOnly) const {
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    std::unique_ptr<TFile> f(TFile::Open((fOutDir+"/ws_combined"+fSaveSuf+".root").c_str()));
    if (!f) {
        LOG(ERROR) << "Cannot open file!\n";
        exit(EXIT_FAILURE);
    }
    RooWorkspace *ws = dynamic_cast<RooWorkspace*>(f->Get("combWS"));

    std::map < std::string, double > result;

    /////////////////////////////////
    //
    // Function performing a fit in a given configuration.
    //
    /////////////////////////////////

    //
    // Fit configuration (1: SPLUSB or 2: BONLY)
    //

    FittingTool fitTool{};
    fitTool.SetUseHesse(fUseHesse);
    fitTool.SetUseHesseBeforeMigrad(fUseHesseBeforeMigrad);
    fitTool.SetStrategy(fFitStrategy);
    fitTool.SetNCPU(fCPU);
    fitTool.SetSpeedUpFit(fSpeedUpFit);
    fitTool.SetMaxFCNcalls(fMaximumNumberFCNcalls);
    fitTool.SetToleranceScale(fToleranceScale);
    fitTool.SetUseAutoDiff(fUseAutoDiff);
    fitTool.SetNLLOffset(fNLLOffset);
    if (fSpeedUpFit) {
        // get the NPoder
        const std::vector<std::pair<std::string,double> >& smallNPs = this->GetSmallNPs();
        fitTool.SetSmallNPs(smallNPs);
    }
    if(fitType==2){
        for (const auto& iFit : fFitList) {
            for (const auto& inf : iFit->fNormFactors) {
                fitTool.AddValPOI(inf->fName, 0.);
            }
        }
        fitTool.ConstPOI(true);
    } else {
        for (const auto& iFit : fFitList) {
            for (const auto& inf : iFit->fNormFactors) {
                fitTool.AddValPOI(inf->fName, inf->GetNominal());
            }
        }
        for (const auto& i : fPOIInitials) {
            fitTool.ReplacePOIVal(i.first, i.second);
        }
        fitTool.ConstPOI(false);
    }
    if(fUseRnd) fitTool.SetRandomNP(fRndRange, fUseRnd, fRndSeed);

    auto nfs = this->GetFitNormFactors();
    std::vector<std::string> constNFs;
    for (const auto& inf : nfs) {
        if (inf->fConst) {
            constNFs.emplace_back(inf->fName);
        }
    }
    fitTool.SetConstantNFs(constNFs);

    //
    // Fit starting from custom point
    if(fFitResultsRootFile!="") {
        fFitList[0]->ReadFitResults(fFitResultsRootFile);
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp : fFitList[0]->fFitResults->GetNuisanceParameters()) {
            npNames.push_back(inp.first);
            npValues.push_back(inp.second->fFitValue);
        }
        fitTool.SetNPs(npNames,npValues);
    }
    // Fix NPs that are specified in the individual configs
    for (const auto& ifit : fFitList){
        if (!ifit->fUseInFit) continue;
        if(ifit->fFitFixedNPs.size()>0){
            std::vector<std::string> npNames;
            std::vector<double> npValues;
            for(const auto& nuisParToFix : ifit->fFitFixedNPs){
                npNames.push_back( nuisParToFix.first );
                npValues.push_back( nuisParToFix.second );
            }
            fitTool.FixNPs(npNames,npValues);
        }
    }

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
    std::vector<std::string> vVarNameMinos = this->GetMinos();
    if(vVarNameMinos.size()>0){
        LOG(DEBUG) << "Setting the variables to use MINOS with:\n";
        for(const auto& iminos : vVarNameMinos) {
            LOG(DEBUG) << "  " << iminos << "\n";
        }
        fitTool.UseMinos(vVarNameMinos);
    }

    //
    // Gets needed objects for the fit
    //
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
    RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>((mc->GetPdf()));
    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    //
    // Creates the data object
    //
    RooDataSet* data = nullptr;
    if(inputData=="asimovData"){
        RooArgSet empty;

        auto* pois = mc->GetParametersOfInterest();
        if (!pois) {
            LOG(ERROR) << "No POIs in the workspace, exiting...\n";
            exit(EXIT_FAILURE);
        }
        for (auto ipoi : *pois) {
            auto itr = std::find(fPOIs.begin(), fPOIs.end(), ipoi->GetName());
            if (itr == fPOIs.end()) continue;
            const std::size_t index = std::distance(fPOIs.begin(), itr);
            const double value = fPOIAsimov.at(index);
            static_cast<RooRealVar*>(ipoi)->setVal(value);
            LOG(INFO) << "Setting POI: " << fPOIs.at(index) << " to value: " << value << "\n";
        }

        data = static_cast<RooDataSet*>(RooStats::AsymptoticCalculator::MakeAsimovData( (*mc), RooArgSet(ws->allVars()), static_cast<RooArgSet&>(empty)));

    }
    else if(inputData!=""){
        data = static_cast<RooDataSet*>(ws->data(inputData.c_str()));
    } else {
        LOG(WARNING) << "You didn't specify inputData => will try with observed data!\n";
        data = static_cast<RooDataSet*>(ws->data("obsData"));
        if(!data){
            LOG(WARNING) << "Observed data not present => will use with asimov data!\n";
            data = static_cast<RooDataSet*>(ws->data("asimovData"));
        }
    }

    if (!data){
        LOG(ERROR) << "Data returns null ptr, probably wrong name in DataName?\n";
        exit(EXIT_FAILURE);
    }

    FitUtils::ApplyExternalConstraints(ws, &fitTool, simPdf, GetFitNormFactors());

    // Performs the fit
    gSystem -> mkdir((fOutDir+"/Fits/").c_str(),true);

    // save snapshot before fit
    ws->saveSnapshot("snapshot_BeforeFit_POI", *(mc->GetParametersOfInterest()) );
    if (mc->GetNuisanceParameters()) ws->saveSnapshot("snapshot_BeforeFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_BeforeFit_GO" , *(mc->GetGlobalObservables())    );

    // save snapshot after fit
    ws->saveSnapshot("snapshot_AfterFit_POI", *(mc->GetParametersOfInterest()) );
    if (mc->GetNuisanceParameters()) ws->saveSnapshot("snapshot_AfterFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(mc->GetGlobalObservables())    );

    std::map<std::string, std::string> mergedMap = this->GetSubcategoryMap();

    if (!doLHscanOnly){
        if (fBlindedParams.size() > 0) std::cout.setstate(std::ios_base::failbit);
        fitTool.FitPDF( mc, simPdf, data, fFastFit );
        if (fBlindedParams.size() > 0) std::cout.clear();

        auto fitResults = fitTool.GetFitResult();

        // error breakdown from the covariance matrix
        bool errorDec(false);
        for (const auto& ifit : fFitList) {
            if (!ifit) continue;
            if (ifit->fErrorDecomposition) {
                errorDec = true;
                break;
            }
        }
        const bool runBreakdown = errorDec && !fDoGroupedSystImpactTable && (fFitType != TRExFit::FitType::BONLY);
        if (runBreakdown) {
            fitTool.SetSystMap( mergedMap );
            LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
            LOG(INFO) << "-----------------------     ERROR DECOMPOSITION    -----------------------\n";
            LOG(INFO) << "  POI :  +/-tot.err  (+/-stat  +/-syst)\n";
            LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
            for (const auto& poi : fPOIs) {
                std::map<std::string,double> impactMap = fitTool.CalculateErrorDecomposition(mc,poi);
                if (impactMap.size() > 0) {
                    LOG(INFO) << Form("  %s :  +/-%g  (+/-%g  +/-%g)",poi.c_str(),impactMap["TOT_ERROR"],impactMap["STAT_ERROR"],impactMap["TOTSYST_ERROR"]) << "\n";
                }
            }
            LOG(INFO) << "  Full breakdown will be written to file in the Fits/ folder.\n";
            LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
            for (const auto& poi : fPOIs) {
                fitTool.ExportErrorDecompositionInTextFile(poi, fOutDir+"/Fits/"+fName+fSaveSuf+"_errDecomp_"+poi+".txt");
                fitTool.ExportErrorDecompositionGroupInTextFile(poi, fOutDir+"/Fits/"+fName+fSaveSuf, fLumiLabel, fCmeLabel, fOutDir, fHEPDataFormat);
            }
            this->PlotRankingFromCovMatrix(fitTool);
        }

        if (fFitType == TRExFit::FitType::UNFOLDING) {
            std::vector<std::string> categories;
            for (const auto& icategory : mergedMap) {
                categories.emplace_back(icategory.second);
            }
            this->PlotUnfoldedErrors(fitTool, categories);
        }

        xRooNode rooNode(*ws);
        rooNode.SetFitResult(*fitResults);
        auto pdfNode = rooNode["simPdf"];
        if (!pdfNode) {
            LOG(ERROR) << "Cannot read simPdf from the WS\n";
            exit(EXIT_FAILURE);
        }
        const double prob = FitUtils::SaturatedModelGoF(pdfNode, inputData);
        LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
        LOG(INFO) << "----------------------- GOODNESS OF FIT EVALUATION -----------------------\n";
        LOG(INFO) << "  probability = " << prob << "\n";
        LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
        fitTool.ExportFitResultInTextFile(fOutDir+"/Fits/"+fName+fSaveSuf+".txt", fBlindedParams);
        result = fitTool.ExportFitResultInMap();
    }

    //
    // grouped systematics impact
    if(fDoGroupedSystImpactTable && !doLHscanOnly){
        std::string outNameGroupedImpact = fOutDir+"/Fits/GroupedImpact"+fSaveSuf;
        if(fGroupedImpactCategory!="all") outNameGroupedImpact += "_"+fGroupedImpactCategory;
        outNameGroupedImpact += ".txt";
        // need to create a merged list of fSubCategoryImpactMap from the include Fits
        fitTool.SetSystMap( mergedMap );
        fitTool.GetGroupedImpact( mc, simPdf, data, ws, fGroupedImpactCategory, outNameGroupedImpact, fOutDir, fLumiLabel, fCmeLabel, fHEPDataFormat);
    }

    //
    // Calls the  function to create LH scan with respect to a parameter
    //
    const std::vector<std::string> parameters = FitUtils::GetAllParameters(mc);
    if(fVarNameLH.size()>0 && !doLHscanOnly && (fLHscanStep == -1 && fLHscanStepY == -1)){ // Skip 1Dscan when paralelizing 2D
        if (fVarNameLH[0]=="all"){
            for(const auto& iparam : parameters) {
                GetLikelihoodScan( ws, iparam, data, true);
            }
        } else{
            for(const auto& iLH : fVarNameLH){
                GetLikelihoodScan( ws, iLH, data, true);
            }
        }
    }
    if (doLHscanOnly && (fLHscanStep == -1 && fLHscanStepY == -1)){ // Skip 1Dscan when paralelizing 2D
        if (fVarNameLH.size() == 0){
            LOG(ERROR) << "Did not provide any LH scan parameter and running LH scan only. This is not correct.\n";
            exit(EXIT_FAILURE);
        }
        if (fVarNameLH[0]=="all"){
            LOG(WARNING) << "You are running LHscan only option but running it for all parameters. Will not parallelize!.\n";
            for(const auto& iparam : parameters) {
                GetLikelihoodScan( ws, iparam, data, true);
            }
        } else {
            GetLikelihoodScan( ws, fVarNameLH[0], data, true);
        }
    }
    // run 2D likelihood scan
    if(fVarName2DLH.size()>0){
        for (const auto & ipair : fVarName2DLH) {
            Get2DLikelihoodScan( ws, ipair, data);
        }
    }

    // Stat-only fit:
    // - read fit resutls
    // - fix all NP to fitted ones before fitting
    if(fIncludeStatOnly && !doLHscanOnly){
        LOG(INFO) << "Fitting stat-only: reading fit results from full fit from file:\n";
        LOG(INFO) << "  " + (fOutDir+"/Fits/"+fName+fSaveSuf+".root") << "\n";
        if (fBlindedParams.size() > 0) std::cout.setstate(std::ios_base::failbit);
        fFitList[0]->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".root");
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp : fFitList[0]->fFitResults->GetNuisanceParameters()) {
            bool isNF = false;
            for(const auto& ifit : fFitList) {
                if (!ifit->fUseInFit) continue;
                if(!ifit->fFixNPforStatOnlyFit &&
                    Common::FindInStringVector(ifit->fNormFactorNames,inp.second->fName)>=0){
                    isNF = true;
                    break;
                }
            }
            if(isNF) continue;
            npNames.push_back(inp.second->fName);
            npValues.push_back(inp.second->fFitValue);
        }
        fitTool.ResetFixedNP();
        fitTool.FixNPs(npNames,npValues);
        fitTool.FitPDF( mc, simPdf, data );
        fitTool.ExportFitResultInTextFile(fOutDir+"/Fits/"+fName+fSaveSuf+"_statOnly.txt", fBlindedParams);
        if (fBlindedParams.size() > 0) std::cout.clear();
    }

    if (fFitToys > 0) {
        this->RunToys();
    }

    return result;
}

//__________________________________________________________________________________
//
void MultiFit::RunLimitScan(xRooNLLVar::xRooHypoSpace& hs) const {
    if (fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
         hs.scan("cls");
    } else if (fLimitType == TRExFit::LimitType::TOYS) {
         if (fLimitToysUsexRooFit == TRExFit::ToysUsexRooFit::AUTO) {
             hs.scan("cls toys");
         } else if (fLimitToysUsexRooFit == TRExFit::ToysUsexRooFit::TRUE) {
            TString sToysScanArg = Common::xRooFitToysScanArg(fLimitToysStepsSplusB, fLimitToysStepsB);
            if (sToysScanArg == "") {
                LOG(ERROR) << "Invalid scan arg\n";
                exit(EXIT_FAILURE);
            }
            TString scanArgs = TString::Format("cls toys=%s", sToysScanArg.Data());
            hs.scan(scanArgs,
                    fLimitToysScanSteps,
                    fLimitToysScanMin,
                    fLimitToysScanMax);
         } else {
             // This should never happen because ConfigReader does not allow this
             LOG(ERROR) << "LimitType is TOYS but ToysUsexRooFit is False!\n";
         }
    } else {
        LOG(ERROR) << "Unkown LimitType\n";
    }
}


//__________________________________________________________________________________
//
void MultiFit::GetCombinedLimit(const std::string& inputData){

    if(fPOIforLimit==""){
        if(fPOIs.size() == 1){
            fPOIforLimit = fPOIs.at(0);
        }
        else{
            LOG(ERROR) << "No POI specified (in 'Limit').\n";
            exit(EXIT_FAILURE);
        }
    }

    const std::string wsFileName = fOutDir+"/ws_combined"+fSaveSuf+".root";

    std::unique_ptr<TFile> f(TFile::Open(wsFileName.c_str()));
    RooWorkspace* ws = f->Get<RooWorkspace>("combWS");
    if (!ws) {
        LOG(ERROR) << "Workspace is nullptr!\n";
        exit(EXIT_FAILURE);
    }

    RooDataSet* data = static_cast<RooDataSet*>(ws->data(inputData.c_str()));
    if (!data) {
        LOG(ERROR) << "Data is nullptr!\n";
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(WARNING) << "Cannot read ModelConfig\n";
        return;
    }

    mc->SetParametersOfInterest(fPOIforLimit.c_str());

    if (fLimitType == TRExFit::LimitType::ASYMPTOTIC ||
        (fLimitType == TRExFit::LimitType::TOYS &&
         fLimitToysUsexRooFit != TRExFit::ToysUsexRooFit::FALSE)) {
        gSystem->mkdir((fOutDir+"/Limits").c_str());
        std::string limitDir = fOutDir + "/Limits/" + (fLimitType == TRExFit::LimitType::ASYMPTOTIC? "Asymptotics" : "Toys");
        gSystem->mkdir(limitDir.c_str());
        const std::string limitFile = limitDir+"/"+fLimitOutputPrefixName+fSaveSuf+".root";
        LOG(INFO) << "Running limit estimate...\n";
        xRooNode model(*ws);

        LOG(DEBUG) << "Setting Hesse strategy to ... " << fLimitFitStrategy << "\n";
        const auto* constraints = FitUtils::GetExternalConstraints(ws, this->GetFitNormFactors(), this->GetRegularizationType());
        auto nll = model.nll(inputData.c_str());
        nll.SetOption(RooFit::Optimize(kTRUE));
        nll.SetOption(RooFit::NumCPU(fCPU, RooFit::Hybrid));
        if (constraints) {
            nll.SetOption(RooFit::ExternalConstraints(*constraints));
        }
        if (fSignificanceUseAutoDiff) {
            LOG(INFO) << "Will use automatic differentiation for the limit estimate\n";
            nll.SetOption(RooFit::EvalBackend("codegen"));
        }
        nll.fitConfigOptions()->SetValue("HesseStrategy",fLimitFitStrategy);
        if (TRExFitter::DEBUGLEVEL >= 2) {
            nll.reinitialize();
            nll.Print();
        }
        auto hs = nll.hypoSpace(fPOIforLimit.c_str());

        const std::string outWSFile = limitDir + "/PostLimitWS" + fSaveSuf + ".root";
        TFile outWS(outWSFile.c_str(), "UPDATE");
        RunLimitScan(hs);
        std::unique_ptr<RooStats::HypoTestInverterResult> result(hs.result());
        result->SetName("clsLimits");
        model.ws()->import(*result);

        std::unique_ptr<xRooNLLVar::xRooHypoSpace> hsInjected(nullptr);
        std::unique_ptr<RooStats::HypoTestInverterResult> injectedResult(nullptr);
        if (fSignalInjection && fSignalInjectionValue > 0) {
            LOG(INFO) << "Injecting signal into the limit estimate: " << fSignalInjectionValue << "\n";
            model.poi().at(fPOIforLimit.c_str())->SetContent(fSignalInjectionValue);

            hsInjected = std::make_unique<xRooNLLVar::xRooHypoSpace>(model.nll().hypoSpace(fPOIforLimit.c_str()));
            RunLimitScan(*hsInjected); // get the observed limit
            injectedResult.reset(hsInjected->result());
            injectedResult->SetName("CLSLimits_injected");
            model.ws()->import(*injectedResult);
        }

        Common::ProcessLimitOutput(hs, hsInjected.get(), "from the combined WS", fLimitIsBlind, limitFile, fLimitParamName, fLimitParamValue, fLimitType == TRExFit::LimitType::ASYMPTOTIC);
        outWS.Close();
        model.SaveAs(outWSFile.c_str(), "UPDATE");

        if (fLimitPlot && fLimitType == TRExFit::LimitType::TOYS) {
            TCanvas c;
            hs.Draw();
            Common::SaveCanvasAs(c, fName + "/Limits/Toys/LimitToysResults" + fLimitOutputPrefixName + fSaveSuf);
        }
    } else if (fLimitType == TRExFit::LimitType::TOYS) {
        RunLimitToys(data, ws);
    }
    f->Close();
}
//__________________________________________________________________________________
//
void MultiFit::GetCombinedSignificance(const std::string& inputData) {
    LOG(INFO) << "Running significance estiamte...\n";

    if(fPOIforSig==""){
        if(fPOIs.size() == 1){
            fPOIforSig = fPOIs.at(0);
        }
        else{
            LOG(ERROR) << "No POI specified (in 'Significance' block).\n";
            exit(EXIT_FAILURE);
        }
    }

    std::string wsFileName = fOutDir+"/ws_combined"+fSaveSuf+".root";

    std::unique_ptr<TFile> f(TFile::Open(wsFileName.c_str()));
    RooWorkspace* ws = f->Get<RooWorkspace>("combWS");
    if (!ws) {
        LOG(ERROR) << "Workspace is nullptr!\n";
        exit(EXIT_FAILURE);
    }
    RooDataSet* data = static_cast<RooDataSet*>(ws->data(inputData.c_str()));
    if (!data) {
        LOG(ERROR) << "Data is nullptr!\n";
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(WARNING) << "Cannot read ModelConfig\n";
        return;
    }

    mc->SetParametersOfInterest(fPOIforSig.c_str());

    if (fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC ||
        (fSignificanceType == TRExFit::SignificanceType::TOYS &&
         fSignificanceToysUsexRooFit)) {
        gSystem->mkdir((fOutDir+"/Significance").c_str());
        std::string significanceDir = fOutDir + "/Significance/" +
            (fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC? "Asymptotics" : "Toys");
        gSystem->mkdir(significanceDir.c_str());
        const std::string significanceFile = significanceDir + fSignificanceParamName + fSaveSuf + ".root";
        const std::string significanceWSfile = significanceDir + "/PostSignificanceWS" + fSaveSuf + ".root";

        const auto* constraints = FitUtils::GetExternalConstraints(ws, this->GetFitNormFactors(), this->GetRegularizationType());
        xRooNode model(*ws);
        auto nll = model.nll(inputData.c_str());
        nll.SetOption(RooFit::Optimize(kTRUE));
        nll.SetOption(RooFit::NumCPU(fCPU, RooFit::Hybrid));
        if (constraints) {
            nll.SetOption(RooFit::ExternalConstraints(*constraints));
        }
        if (fSignificanceUseAutoDiff) {
            LOG(INFO) << "Will use automatic differentiation for the significance estimate\n";
            nll.SetOption(RooFit::EvalBackend("codegen"));
        }
        nll.fitConfigOptions()->SetValue("HesseStrategy",fSignificanceFitStrategy);
        if (TRExFitter::DEBUGLEVEL >= 2) {
            nll.reinitialize();
            nll.Print();
        }
        auto hs = nll.hypoSpace(fPOIforSig.c_str(), xRooFit::Asymptotics::Uncapped, 1);

        TFile outWS(significanceWSfile.c_str(), "UPDATE");
        RunSignificanceScan(hs);

        std::unique_ptr<RooStats::HypoTestInverterResult> result(hs.result());
        const std::pair<double, double> pvalueWithErrorExpected = hs[0].pNull_asymp(0);
        const std::pair<double, double> pvalueWithErrorObserved  = hs[0].pNull_asymp();
        std::pair<double, double> pvalueWithErrorInjected = {-99, -99};
        model.ws()->import(*result);

        std::unique_ptr<xRooNLLVar::xRooHypoSpace> hsInjected(nullptr);
        std::unique_ptr<RooStats::HypoTestInverterResult> injectedResult(nullptr);
        if (fSignificanceDoInj && fSignificancePOIAsimov > 0) {
            LOG(INFO) << "Injecting signal into the significance estimate: " << fSignificancePOIAsimov << "\n";
            model.poi().at(fPOIforSig.c_str())->SetContent(fSignificancePOIAsimov);

            hsInjected = std::make_unique<xRooNLLVar::xRooHypoSpace>(model.nll().hypoSpace(fPOIforSig.c_str(), xRooFit::Asymptotics::Uncapped, fSignificancePOIAsimov));
            RunSignificanceScan(*hsInjected);
            pvalueWithErrorInjected = (*hsInjected)[0].pNull_asymp();
            injectedResult.reset(hsInjected->result());
            injectedResult->SetName(TString::Format("%s_%s=%g",injectedResult->GetName(),fPOIforSig.c_str(),fSignificancePOIAsimov));
            model.ws()->import(*injectedResult);
        }
        if (Common::ProcessSignificanceOutput(pvalueWithErrorObserved, pvalueWithErrorExpected, pvalueWithErrorInjected, "for combined fit", fSignificanceIsBlind, significanceFile)) {
            LOG(WARNING) << "NaN value discovered, printing more information\n";
            hs.Print();
            if (hsInjected) {
                LOG(WARNING) << "Injected significance information\n";
                hsInjected->Print();
            }
        }
        model.SaveAs(significanceWSfile.c_str());
    } else {
        RunSignificanceToys(data, ws);

        f->Close();
    }
}

//__________________________________________________________________________________
//
void MultiFit::RunSignificanceScan(xRooNLLVar::xRooHypoSpace& hs) const {
    if (fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
        hs.scan("pnull", 1, 0, 0);
    } else if (fSignificanceType == TRExFit::SignificanceType::TOYS) {
        if (fSignificanceToysUsexRooFit) {
            TString sToysScanArg = Common::xRooFitToysScanArg(fSignificanceToysStepsSplusB, fSignificanceToysStepsB);
            if (sToysScanArg == "") {
                LOG(ERROR) << "Invalid scan arg\n";
                exit(EXIT_FAILURE);
            }
            TString scanArgs = TString::Format("pnull toys=%s", sToysScanArg.Data());
            hs.scan(scanArgs, 1, 0, 0);
        } else {
            // This should never happen because ConfigReader does not allow this
            LOG(ERROR) << "SignificanceType is TOYS but ToysUsexRooFit is False!\n";
            exit(EXIT_FAILURE);
        }
    } else {
        LOG(ERROR) << "Unkown SignificanceType\n";
        exit(EXIT_FAILURE);
    }
}

//__________________________________________________________________________________
//
void MultiFit::ComparePOI(const string& POI, const std::size_t index) const {
    if (std::find(fBlindedParams.begin(), fBlindedParams.end(), POI) != fBlindedParams.end()) {
        LOG(INFO) << "Parameter " << POI << " is blinded. Will not produce the POI plot.\n";
        return;
    }
    double xmax(0);
    double xmin(0);

    if (fPOIMax.empty() || fPOIMin.empty()) {
        xmax = 20;
        xmin = 0;
    } else {
        xmax = fPOIMax.at(index) + (fPOIMax.at(index)-fPOIMin.at(index));
        xmin = fPOIMin.at(index);
    }

    string process = fLabel;

    // Fit titles
    vector<string> names;
    vector<string> dirs;
    vector<string> suffs;
    vector<string> titles;
    vector<string> pois;
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        LOG(INFO) << "Adding Fit: " << fFitList[i_fit]->fInputName << ", " << fFitLabels[i_fit] << ", " << fFitSuffs[i_fit] << "\n";
        names.push_back( fFitList[i_fit]->fInputName );
        dirs.push_back( fFitList[i_fit]->fName );
        titles.push_back( fFitLabels[i_fit] );
        suffs.push_back( fFitSuffs[i_fit] );
        if (fFitList[i_fit]->fPOIs.empty()) pois.push_back( POI );
        else {
            const int tmp = Common::FindInStringVector(fFitList[i_fit]->fPOIs,POI);
            if(tmp >= 0){
                pois.push_back(fFitList[i_fit]->fPOIs[tmp]);
            } else {
                pois.push_back(POI);
            }
        }
    }
    if(fCombine){
        LOG(INFO) << "Adding Combined Fit\n";
        names.push_back( fName );
        dirs.push_back( fOutDir );
        titles.push_back( fCombiLabel );
        suffs.push_back( "" );
        pois.push_back( POI );
    }

    int N = names.size();

    double ymin = -0.5;
    double ymax = N+1-0.5;
    double yaxis_extension = 0.5;

    TCanvas c("c","c",700,500);
    gStyle->SetEndErrorSize(6.);

    TGraph g_central(N);
    TGraphAsymmErrors g_stat(N);
    TGraphAsymmErrors g_tot(N);

    int Ndiv = N+1;

    NuisParameter *par;

    // get values
    TRExFit *fit = nullptr;
    for(int i=0;i<N;i++){
        if(!(fCombine && i==N-1)){
            if (!fFitList.at(i)->fUseInComparison) continue;
        }
        const bool isComb = (fCombine && i==N-1) ? true : false;
        //
        if(!isComb){
            fit = fFitList[i].get();
            if(fit->fFitResultsRootFile=="") fit->ReadFitResults(dirs[i]+"/Fits/"+names[i]+suffs[i]+".root");
            else                             fit->ReadFitResults(fit->fFitResultsRootFile);
        }
        else{
            if(fFitResultsRootFile=="") fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".root");
            else                        fit->ReadFitResults(fFitResultsRootFile);
        }
        bool found(false);
        for(const auto& inp : fit->fFitResults->GetNuisanceParameters()){
            par = inp.second.get();
            if( pois[i] == par->fName ){
                g_central.SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_stat   .SetPoint(N-i-1,par->fFitValue,N-i-1);
                g_tot    .SetPoint(N-i-1,par->fFitValue,N-i-1);
                //
                // temporary put the full uncertainty
                g_stat.SetPointEXhigh(N-i-1,par->fPostFitUp);
                g_stat.SetPointEXlow(N-i-1,-par->fPostFitDown);
                g_stat.SetPointEYhigh(N-i-1,0);
                g_stat.SetPointEYlow(N-i-1,0);
                //
                g_tot.SetPointEXhigh(N-i-1,par->fPostFitUp);
                g_tot.SetPointEXlow(N-i-1,-par->fPostFitDown);
                g_tot.SetPointEYhigh(N-i-1,0);
                g_tot.SetPointEYlow(N-i-1,0);
                //
                found = true;
                break;
            }
        }
        if(!found){
            g_central.SetPoint(N-i-1,-10,N-i-1);
            g_stat   .SetPoint(N-i-1,-10,N-i-1);
            g_tot    .SetPoint(N-i-1,-10,N-i-1);
            g_stat   .SetPointError(N-i-1,0,0,0,0);
            g_tot    .SetPointError(N-i-1,0,0,0,0);
        }
    }
    // stat error
    if (!fShowTotalOnly){
        for(int i=0;i<N;i++){
            const bool isComb = (fCombine && i==N-1);
            if(!isComb){
                fit = fFitList.at(i).get();
                fit->ReadStatOnlyFromErrorDecomposition(fit->fPOIs, dirs.at(i)+"/Fits/"+names.at(i)+suffs.at(i)+"_errDecomp");
            }
            else{
                fit->ReadStatOnlyFromErrorDecomposition(fPOIs, fOutDir+"/Fits/"+fName+fSaveSuf+"_errDecomp");
            }
            double statUp = 0;
            double statDown = 0;

            if (isComb) {
                statUp   = fit->GetStatOnlyErrorFromDecompositionUp(pois.at(i));
                statDown = -fit->GetStatOnlyErrorFromDecompositionDown(pois.at(i));
            } else {
                auto itr = std::find(fit->fPOIs.begin(), fit->fPOIs.end(), pois.at(i));
                if (itr == fit->fPOIs.end()) {
                    statUp = -1;
                    statDown = -1;
                } else {
                    statUp   = fit->GetStatOnlyErrorFromDecompositionUp(pois.at(i));
                    statDown = -fit->GetStatOnlyErrorFromDecompositionDown(pois.at(i));
                }
            }

            g_stat.SetPointEXhigh(N-i-1, statUp);
            g_stat.SetPointEXlow(N-i-1, statDown);
            g_stat.SetPointEYhigh(N-i-1,0);
            g_stat.SetPointEYlow(N-i-1,0);
        }
    }

    g_tot .SetLineWidth(3);
    g_stat.SetLineWidth(3);
    g_stat.SetMarkerStyle(kOpenCircle);
    g_tot .SetLineWidth(3);
    g_tot .SetMarkerStyle(kOpenCircle);
    g_stat.SetLineColor(kGreen-8);
    g_stat.SetMarkerSize(0);
    g_tot .SetLineColor(kBlack);
    g_central.SetMarkerStyle(kFullCircle);
    g_central.SetMarkerColor(kRed);
    g_central.SetMarkerSize(1.5);
    g_tot.SetMarkerSize(0);

    TH1D h_dummy("h_dummy","h_dummy",1,xmin,xmax);
    h_dummy.Draw();
    h_dummy.SetMinimum(ymin);
    h_dummy.SetMaximum(ymax+yaxis_extension);
    h_dummy.SetLineColor(kWhite);
    h_dummy.GetYaxis()->Set(N+1,ymin,ymax);
    h_dummy.GetYaxis()->SetNdivisions(Ndiv);

    TLatex tex{};

    const std::string precision = fPOIPrecision.empty() ? "1" : fPOIPrecision.at(index);

    for(int i=0;i<N;i++){
        if(!(fCombine && i==N-1)){
            if (!fFitList.at(i)->fUseInComparison) continue;
        }
        h_dummy.GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
        if(fShowSystForPOI){
            tex.SetTextSize(gStyle->GetTextSize()*1.2);
            tex.DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form(("#font[62]{%." + precision + "f}").c_str(),g_central.GetX()[N-i-1]));
            tex.DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{^{#plus%." + precision + "f}}").c_str(),g_tot.GetErrorXhigh(N-i-1)));
            tex.DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{_{#minus%." + precision + "f}}").c_str(),g_tot.GetErrorXlow(N-i-1)));
            tex.DrawLatex(xmin+0.69*(xmax-xmin),N-i-1,"(");
            if (!fShowTotalOnly){
                tex.DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + precision + "f}}").c_str(),g_stat.GetErrorXhigh(N-i-1)));
                tex.DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + precision + "f}}").c_str(),g_stat.GetErrorXlow(N-i-1)));
                tex.DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + precision + "f}}").c_str(),
                    std::sqrt( g_tot.GetErrorXhigh(N-i-1) * g_tot.GetErrorXhigh(N-i-1) - g_stat.GetErrorXhigh(N-i-1)* g_stat.GetErrorXhigh(N-i-1) ) ) );
                tex.DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + precision + "f}}").c_str(),
                    std::sqrt( g_tot.GetErrorXlow(N-i-1)*g_tot.GetErrorXlow(N-i-1) - g_stat.GetErrorXlow(N-i-1)*g_stat.GetErrorXlow(N-i-1) ) ) );
                tex.DrawLatex(xmin+0.94*(xmax-xmin),N-i-1,")");
            }
        }
        else{
            tex.DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form((fPOIName.at(index)+" = %." + precision + "f").c_str(),g_central.GetX()[N-i-1]));
            tex.DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("^{#plus%." + precision + "f}").c_str(),g_tot.GetErrorXhigh(N-i-1)));
            tex.DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("_{#minus%." + precision + "f}").c_str(),g_tot.GetErrorXlow(N-i-1)));
            if (!fShowTotalOnly){
                tex.DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("^{#plus%." + precision + "f}").c_str(),g_stat.GetErrorXhigh(N-i-1)));
                tex.DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("_{#minus%." + precision + "f}").c_str(),g_stat.GetErrorXlow(N-i-1)));
            }
        }
    }

    TLine l_SM(fPOINominal.at(index),-0.5,fPOINominal.at(index),N-0.5);
    l_SM.SetLineWidth(2);
    l_SM.SetLineColor(kGray);
    l_SM.Draw("same");

    TLine l_h(xmin,0.5,xmax,0.5);
    if(fCombine){
        l_h.SetLineWidth(2);
        l_h.SetLineColor(kBlack);
        l_h.SetLineStyle(kDashed);
        l_h.Draw("same");
    }

    g_tot.Draw("E same");
    if (!fShowTotalOnly){
        g_stat.Draw("E same");
    }
    g_central.Draw("P same");

    gPad->SetLeftMargin( 2*gPad->GetLeftMargin() );
    gPad->SetBottomMargin(0.184);
    gPad->SetTopMargin(0.09);
    h_dummy.GetXaxis()->SetTitle(fPOITitle.at(index).c_str());
    h_dummy.GetYaxis()->SetTickSize(0);

    c.RedrawAxis();

    // Placement of Labels is in NDC
    gPad->Update();
    double yLabelNDC = (((ymax+yaxis_extension)*1.01) - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());

    if (fFitList[0]->fPlotLabel != "none") TRExLabelNew(0.32, yLabelNDC, TRExFitter::EXPERIMENT_LABEL.c_str(), fFitList[0]->fPlotLabel.c_str(),kBlack,gStyle->GetTextSize(), 0.10);
    myText(0.68, yLabelNDC, kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    if(process!="") {
        double yProcessLabelNDC = ((ymax*0.98) - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
        myText(0.94, yProcessLabelNDC, kBlack,Form("#kern[-1]{%s}", process.c_str()));
    }

    double yLegendNDC = ((ymax*1.04) - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
    TLegend leg(0.35, yLegendNDC-0.125, 0.7, yLegendNDC);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(&g_tot,"tot.","l");
    if (!fShowTotalOnly){
        leg.AddEntry(&g_stat,"stat.","l");
    }
    leg.Draw();

    if(fShowSystForPOI){
        tex.DrawLatex(xmin+0.6*(xmax-xmin),N-0.4,"#font[62]{tot}");
        tex.DrawLatex(xmin+0.69*(xmax-xmin),N-0.4,"(");
        if (!fShowTotalOnly){
            tex.DrawLatex(xmin+0.72*(xmax-xmin),N-0.4,"stat");
        }
        tex.DrawLatex(xmin+0.83*(xmax-xmin),N-0.4,"syst");
        tex.DrawLatex(xmin+0.94*(xmax-xmin),N-0.4,")");
    }
    else{
        tex.DrawLatex(xmin+(0.7-0.02)*(xmax-xmin),N-0.4,"( tot )");
        if (!fShowTotalOnly){
            tex.DrawLatex(xmin+(0.85-0.02)*(xmax-xmin),N-0.4,"( stat )");
        }
    }

    Common::SaveCanvasAs(c, fOutDir+"/POI_"+fPOIs.at(index)+"_"+fSaveSuf);
}

//__________________________________________________________________________________
//
void MultiFit::CompareEFT() const {
    bool isEFT(false);
    for (const auto& ifit : fFitList) {
        if (ifit->fFitType == TRExFit::FitType::EFT) {
            isEFT = true;
            break;
        }
    }
    if (!isEFT) return;
    double xmax(0);
    double xmin(0);

    if (fPOIMax.empty() || fPOIMin.empty()) {
        xmax = 20;
        xmin = 0;
    } else {
        xmax = fPOIMax.at(0) + (fPOIMax.at(0)-fPOIMin.at(0));
        xmin = fPOIMin.at(0);
    }

    string process = fLabel;

    // Fit titles
    vector<string> names;
    vector<string> dirs;
    vector<string> suffs;
    vector<string> titles;
    vector<string> pois;
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        LOG(INFO) << "Adding Fit: " << fFitList[i_fit]->fInputName << ", " << fFitLabels[i_fit] << ", " << fFitSuffs[i_fit] << "\n";
        names.push_back( fFitList[i_fit]->fInputName );
        dirs.push_back( fFitList[i_fit]->fName );
        titles.push_back( fFitLabels[i_fit] );
        suffs.push_back( fFitSuffs[i_fit] );
        //const int imatch = Common::FindInStringVector(pois,par->fName);
        pois.push_back(fFitList[i_fit]->fPOIs[0]);
    }

    int N = names.size();

    double ymin = -0.5;
    double ymax = N+1-0.5;
    double yaxis_extension = 0.5;

    TCanvas c("c","c",700,500);
    gStyle->SetEndErrorSize(6.);

    //std::vector< TGraphAsymmErrors> g;
    TGraph g_central(N*4);
    TGraphAsymmErrors g_stat(N*4);
    TGraphAsymmErrors g_tot(N*4);

    int Ndiv = N+1;

    NuisParameter *par;
    // get values
    TRExFit *fit = nullptr;
    for(int i=0;i<N;i++){
        if(!(fCombine && i==N-1)){
            if (!fFitList.at(i)->fUseInComparison) continue;
        }
        const bool isComb = (fCombine && i==N-1) ? true : false;
        //
        if(!isComb){
            fit = fFitList.at(i).get();
            if(fit->fFitResultsRootFile=="") fit->ReadFitResults(dirs[i]+"/Fits/"+names[i]+suffs[i]+".root");
            else                             fit->ReadFitResults(fit->fFitResultsRootFile);
        }
        else{
            if(fFitResultsRootFile=="") fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".root");
            else                        fit->ReadFitResults(fFitResultsRootFile);
        }
        bool found(false);
        float ypos=N-i-1;
        int nPoint=ypos;
        for(const auto& inp : fit->fFitResults->GetNuisanceParameters()) {
            par = inp.second.get();
            if( par->fCategory == "EFT" ){
                g_central.SetPoint(nPoint,par->fFitValue,ypos);
                g_stat   .SetPoint(nPoint,par->fFitValue,ypos);
                g_tot    .SetPoint(nPoint,par->fFitValue,ypos);
                //
                // temporary put the full uncertainty
                g_stat.SetPointEXhigh(nPoint,par->fPostFitUp);
                g_stat.SetPointEXlow(nPoint,-par->fPostFitDown);
                g_stat.SetPointEYhigh(nPoint,0);
                g_stat.SetPointEYlow(nPoint,0);
                //
                g_tot.SetPointEXhigh(nPoint,par->fPostFitUp);
                g_tot.SetPointEXlow(nPoint,-par->fPostFitDown);
                g_tot.SetPointEYhigh(nPoint,0);
                g_tot.SetPointEYlow(nPoint,0);
                //
                found = true;
            }
        }
        if(!found){
            g_central.SetPoint(nPoint,-10,ypos);
            g_stat   .SetPoint(nPoint,-10,ypos);
            g_tot    .SetPoint(nPoint,-10,ypos);
            g_stat   .SetPointError(nPoint,0,0,0,0);
            g_tot    .SetPointError(nPoint,0,0,0,0);
        }
    }

    g_tot .SetLineWidth(3);
    g_stat.SetLineWidth(3);
    g_stat.SetMarkerStyle(kOpenCircle);
    g_tot .SetLineWidth(3);
    g_tot .SetMarkerStyle(kOpenCircle);
    g_stat.SetLineColor(kGreen-8);
    g_stat.SetMarkerSize(0);
    g_tot .SetLineColor(kBlack);
    g_central.SetMarkerStyle(kFullCircle);
    g_central.SetMarkerColor(kRed);
    g_central.SetMarkerSize(1.5);
    g_tot.SetMarkerSize(0);

    TH1D h_dummy("h_dummy","h_dummy",1,xmin,xmax);
    h_dummy.Draw();
    h_dummy.SetMinimum(ymin);
    h_dummy.SetMaximum(ymax+yaxis_extension);
    h_dummy.SetLineColor(kWhite);
    h_dummy.GetYaxis()->Set(N+1,ymin,ymax);
    h_dummy.GetYaxis()->SetNdivisions(Ndiv);

    TLatex tex{};

    const std::string precision = fPOIPrecision.empty() ? "1" : fPOIPrecision.at(0);

    for(int i=0;i<N;i++){
        if(!(fCombine && i==N-1)){
            if (!fFitList.at(i)->fUseInComparison) continue;
        }
        h_dummy.GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
        if(fShowSystForPOI){
            tex.SetTextSize(gStyle->GetTextSize()*1.2);
            tex.DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form(("#font[62]{%." + precision + "f}").c_str(),g_central.GetX()[N-i-1]));
            tex.DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{^{#plus%." + precision + "f}}").c_str(),g_tot.GetErrorXhigh(N-i-1)));
            tex.DrawLatex(xmin+0.6*(xmax-xmin),N-i-1,Form(("#font[62]{_{#minus%." + precision + "f}}").c_str(),g_tot.GetErrorXlow(N-i-1)));
            tex.DrawLatex(xmin+0.69*(xmax-xmin),N-i-1,"(");
            if (!fShowTotalOnly){
                tex.DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + precision + "f}}").c_str(),g_stat.GetErrorXhigh(N-i-1)));
                tex.DrawLatex(xmin+0.73*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + precision + "f}}").c_str(),g_stat.GetErrorXlow(N-i-1)));
                tex.DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{^{#plus%." + precision + "f}}").c_str(),
                    std::sqrt( g_tot.GetErrorXhigh(N-i-1) * g_tot.GetErrorXhigh(N-i-1) - g_stat.GetErrorXhigh(N-i-1)* g_stat.GetErrorXhigh(N-i-1) ) ) );
                tex.DrawLatex(xmin+0.84*(xmax-xmin),N-i-1,Form(("#font[42]{_{#minus%." + precision + "f}}").c_str(),
                    std::sqrt( g_tot.GetErrorXlow(N-i-1)*g_tot.GetErrorXlow(N-i-1) - g_stat.GetErrorXlow(N-i-1)*g_stat.GetErrorXlow(N-i-1) ) ) );
                tex.DrawLatex(xmin+0.94*(xmax-xmin),N-i-1,")");
            }
        }
        else{
            tex.DrawLatex(xmin+0.5*(xmax-xmin),N-i-1,Form(("%." + precision + "f").c_str(),g_central.GetX()[N-i-1]));
            tex.DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("^{#plus%." + precision + "f}").c_str(),g_tot.GetErrorXhigh(N-i-1)));
            tex.DrawLatex(xmin+0.7*(xmax-xmin),N-i-1,Form(("_{#minus%." + precision + "f}").c_str(),g_tot.GetErrorXlow(N-i-1)));
            if (!fShowTotalOnly){
                tex.DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("^{#plus%." + precision + "f}").c_str(),g_stat.GetErrorXhigh(N-i-1)));
                tex.DrawLatex(xmin+0.85*(xmax-xmin),N-i-1,Form(("_{#minus%." + precision + "f}").c_str(),g_stat.GetErrorXlow(N-i-1)));
            }
        }
    }

    TLine l_SM(fPOINominal.at(0),-0.5,fPOINominal.at(0),N-0.5);
    l_SM.SetLineWidth(2);
    l_SM.SetLineColor(kGray);
    l_SM.Draw("same");

    TLine l_h(xmin,0.5,xmax,0.5);
    if(fCombine){
        l_h.SetLineWidth(2);
        l_h.SetLineColor(kBlack);
        l_h.SetLineStyle(kDashed);
        l_h.Draw("same");
    }

    g_tot.Draw("E same");
    if (!fShowTotalOnly){
        g_stat.Draw("E same");
    }
    g_central.Draw("P same");

    gPad->SetLeftMargin( 2*gPad->GetLeftMargin() );
    gPad->SetBottomMargin(0.184);
    gPad->SetTopMargin(0.09);
    //h_dummy.GetXaxis()->SetTitle(fPOITitle.at(0).c_str());
    h_dummy.GetYaxis()->SetTickSize(0);

    c.RedrawAxis();

    // Placement of Labels is in NDC
    gPad->Update();
    double yLabelNDC = (((ymax+yaxis_extension)*1.01) - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());

    if (fFitList[0]->fPlotLabel != "none") TRExLabelNew(0.32, yLabelNDC, TRExFitter::EXPERIMENT_LABEL.c_str(), fFitList[0]->fPlotLabel.c_str(),kBlack,gStyle->GetTextSize(), 0.10);
    myText(0.68, yLabelNDC, kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    if(process!="") {
        double yProcessLabelNDC = ((ymax*0.98) - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
        myText(0.94, yProcessLabelNDC, kBlack,Form("#kern[-1]{%s}", process.c_str()));
    }

    double yLegendNDC = ((ymax*1.04) - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
    TLegend leg(0.35, yLegendNDC-0.125, 0.7, yLegendNDC);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(&g_tot,"tot.","l");
    if (!fShowTotalOnly){
        leg.AddEntry(&g_stat,"stat.","l");
    }
    leg.Draw();

    if(fShowSystForPOI){
        tex.DrawLatex(xmin+0.6*(xmax-xmin),N-0.4,"#font[62]{tot}");
        tex.DrawLatex(xmin+0.69*(xmax-xmin),N-0.4,"(");
        if (!fShowTotalOnly){
            tex.DrawLatex(xmin+0.72*(xmax-xmin),N-0.4,"stat");
        }
        tex.DrawLatex(xmin+0.83*(xmax-xmin),N-0.4,"syst");
        tex.DrawLatex(xmin+0.94*(xmax-xmin),N-0.4,")");
    }
    else{
        tex.DrawLatex(xmin+(0.7-0.02)*(xmax-xmin),N-0.4,"( tot )");
        if (!fShowTotalOnly){
            tex.DrawLatex(xmin+(0.85-0.02)*(xmax-xmin),N-0.4,"( stat )");
        }
    }

    Common::SaveCanvasAs(c, fOutDir+"/EFT_"+fPOIs.at(0)+"_"+fSaveSuf);
}

//__________________________________________________________________________________
//
void MultiFit::CompareLimit(){
    float xmax = 2.;
    const std::string process = fLabel;
    gStyle->SetEndErrorSize(0.);

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    for(unsigned int i_fit= 0; i_fit < fFitList.size(); ++i_fit) {
        LOG(INFO) << "Adding Fit: " << fFitList[i_fit]->fInputName << ", " << fFitLabels[i_fit] << ", " << fFitSuffs[i_fit] << "\n";
        // Every fit is an instance of TRExFit extracted from a config file
        // with fit file potentally overwritten by FitResultsFile option
        dirs.emplace_back(fFitList[i_fit]->fName);
        names.emplace_back(fFitList[i_fit]->fInputName);
        titles.emplace_back(fFitLabels[i_fit]);
        suffs.emplace_back(fFitSuffs[i_fit]);
    }
    if(fCombine){
        LOG(INFO) << "Adding combined limit\n";
        dirs.emplace_back(fOutDir);
        names.emplace_back(fName);
        titles.emplace_back(fCombiLabel);
        suffs.emplace_back(fSaveSuf);
    }

    // ---

    const std::size_t N = names.size();

    const double ymin = -0.5;
    const double ymax = N-0.5;

    TCanvas c("c","c",700,500);

    TGraphErrors g_obs(N);
    TGraphErrors g_exp(N);
    TGraphErrors g_inj(N);
    TGraphAsymmErrors g_1s(N);
    TGraphAsymmErrors g_2s(N);

    const std::size_t Ndiv = N+1;

    std::unique_ptr<TFile> f = nullptr;

    // get values
    for(unsigned int i = 0; i < N; ++i) {
        // fLimits file is set during MultiFit config reading to "" if no limits file is provided
        if(i > fLimitsFiles.size()-1){
            if(fLimitsFile != "") fLimitsFiles.emplace_back(fLimitsFile);
            else                  fLimitsFiles.emplace_back("");
        }

        std::string limitFilePath = fLimitsFiles.at(i);
        if(limitFilePath == ""){
            std::string name("");
            if (i < fFitList.size()) {
                name = fFitList.at(i)->fLimitOutputPrefixName;
            } else {
                name = fLimitOutputPrefixName;
            }
            limitFilePath = dirs.at(i)+"/Limits/Asymptotics/"+name+suffs.at(i)+".root";
            LOG(INFO) << "Reading file " << limitFilePath << "\n";
            f.reset(TFile::Open((limitFilePath).c_str()));
        }
        else{
            LOG(INFO) << "Reading file " << fLimitsFiles[i] << "\n";
            f.reset(TFile::Open(limitFilePath.c_str()));
        }
        if (!f) {
            LOG(WARNING) << "Cannot open the file " << limitFilePath << "\n";
            continue;
        }
        TTree* tree(f->Get<TTree>("stats"));
        if (!tree) {
            LOG(WARNING) << "Cannot read tree from the limit file\n";
            continue;
        }

        if (tree->GetEntries() == 0) {
            LOG(WARNING) << "The tree has 0 entries\n";
            continue;
        }

        LimitEvent event(tree);

        // read the first event = read all info
        event.GetEntry(0);

        if(fShowObserved) g_obs.SetPoint(N-i-1,event.obs_upperlimit,N-i-1);
        else g_obs.SetPoint(N-i-1,-1,N-i-1);
        g_exp.SetPoint(N-i-1,event.exp_upperlimit,N-i-1);
        if(fSignalInjection) g_inj.SetPoint(N-i-1,event.inj_upperlimit,N-i-1);
        g_1s.SetPoint(N-i-1,event.exp_upperlimit,N-i-1);
        g_2s.SetPoint(N-i-1,event.exp_upperlimit,N-i-1);
        g_obs.SetPointError(N-i-1,0,0.5);
        g_exp.SetPointError(N-i-1,0,0.5);
        g_inj.SetPointError(N-i-1,0,0.5);
        g_1s.SetPointError(N-i-1,event.exp_upperlimit-event.exp_upperlimit_minus1,event.exp_upperlimit_plus1-event.exp_upperlimit,0.5,0.5);
        g_2s.SetPointError(N-i-1,event.exp_upperlimit-event.exp_upperlimit_minus2,event.exp_upperlimit_plus2-event.exp_upperlimit,0.5,0.5);

        const std::vector<float> tmp = {event.obs_upperlimit, event.exp_upperlimit,
                                        event.exp_upperlimit_plus1, event.exp_upperlimit_minus1,
                                        event.exp_upperlimit_plus2, event.exp_upperlimit_minus2, xmax};

        xmax = *std::max_element(tmp.begin(), tmp.end());

        f->Close();
    }

    xmax*= 1.1;

    g_obs.SetLineWidth(3);
    g_exp.SetLineWidth(3);
    g_exp.SetLineStyle(2);
    g_inj.SetLineWidth(3);
    g_inj.SetLineStyle(2);
    g_inj.SetLineColor(kRed);
    g_1s.SetFillColor(kGreen);
    g_1s.SetLineWidth(3);
    g_1s.SetLineStyle(2);
    g_2s.SetFillColor(kYellow);
    g_2s.SetLineWidth(3);
    g_2s.SetLineStyle(2);

    g_2s.SetMarkerSize(0);
    g_1s.SetMarkerSize(0);
    g_exp.SetMarkerSize(0);
    g_obs.SetMarkerSize(0);
    g_inj.SetMarkerSize(0);

    if(fLimitMax!=0) xmax = fLimitMax;

    TH1D h_dummy("h_dummy","h_dummy",1,0,xmax);
    h_dummy.Draw();
    h_dummy.SetMinimum(ymin);
    h_dummy.SetMaximum(ymax);
    h_dummy.SetLineColor(kWhite);
    h_dummy.GetYaxis()->Set(N,ymin,ymax);
    h_dummy.GetYaxis()->SetNdivisions(Ndiv);
    for(unsigned int i = 0; i < N; ++i){
        h_dummy.GetYaxis()->SetBinLabel(N-i,titles[i].c_str());
    }

    g_2s.Draw("E2 same");
    g_1s.Draw("E2 same");
    g_exp.Draw("E same");
    if(fShowObserved) g_obs.Draw("E same");
    if(fSignalInjection) g_inj.Draw("E same");

    TLine l_SM(fPOINominal.at(0),-0.5,fPOINominal.at(0),N-0.5);
    l_SM.SetLineWidth(2);
    l_SM.SetLineColor(kGray);
    l_SM.Draw("same");

    c.RedrawAxis();

    gPad->SetLeftMargin(2*gPad->GetLeftMargin());
    gPad->SetBottomMargin(1.15*gPad->GetBottomMargin());
    gPad->SetTopMargin(1.8*gPad->GetTopMargin());
    h_dummy.GetXaxis()->SetTitle(fLimitTitle.c_str());

    if (fFitList[0]->fPlotLabel != "none") TRExLabel(0.32,0.93, TRExFitter::EXPERIMENT_LABEL.c_str(), fFitList[0]->fPlotLabel.c_str(),kBlack);
    myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(process!="") myText(0.94,0.85,kBlack,Form("#kern[-1]{%s}",process.c_str()));

    std::unique_ptr<TLegend> leg = nullptr;
    if(fShowObserved) leg = std::make_unique<TLegend>(0.65,0.2,0.95,0.40);
    else              leg = std::make_unique<TLegend>(0.65,0.2,0.95,0.35);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(&g_1s,"Expected #pm 1#sigma","lf");
    leg->AddEntry(&g_2s,"Expected #pm 2#sigma","lf");
    if(fShowObserved)    leg->AddEntry(&g_obs,"Observed","l");
    if(fSignalInjection) leg->AddEntry(&g_inj,("Expected ("+fPOIName.at(0)+"=1)").c_str(),"l");
    leg->Draw();

    Common::SaveCanvasAs(c, fOutDir+"/Limits" + fSaveSuf );
}

//__________________________________________________________________________________
//
void MultiFit::ComparePulls(const std::string& category) const{
    double ydist = 0.2;

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    vector<double>  yshift;
    static const std::vector<int> color = {kBlack,kRed,kBlue,kViolet,kOrange-7, kGreen+2};

    static const std::vector<int> style = {kFullCircle,kOpenCircle,kFullTriangleUp,kOpenTriangleDown,kOpenDiamond};

    unsigned int N = fFitList.size();
    if(fCombine) N++;

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1){
            LOG(INFO) << "Adding Combined Fit\n";
            dirs.push_back( fOutDir );
            names.push_back( fName );
            titles.push_back( fCombiLabel );
            suffs.push_back( "" );
        }
        else{
            dirs.push_back( fFitList[i_fit]->fName );
            names.push_back( fFitList[i_fit]->fInputName );
            titles.push_back( fFitLabels[i_fit] );
            suffs.push_back( fFitSuffs[i_fit] );
        }
        yshift.push_back( 0. - ydist*N/2. + ydist*i_fit );
    }

    double xmin = -2.9;
    double xmax = 2.9;
    double max = 0.;

    // create a list of Systematics
    std::vector< string > Names;
    std::vector< string > Titles;
    std::vector< string > Categories;
    std::string systName;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        for(const auto& isyst : fFitList[i_fit]->fSystematics) {
            systName = isyst->fNuisanceParameter;
            if(systName == "") systName = isyst->fName;
            if(Common::FindInStringVector(Names,systName)<0){
                Names.push_back(systName);
                Titles.push_back(isyst->fTitle);
                if(Common::FindInStringVector(fFitList[i_fit]->fDecorrSysts,isyst->fName)>=0){
                    Titles[Titles.size()-1] += fFitList[i_fit]->fDecorrSuff;
                }
                if (isyst->fSubCategory != "") {
                  Categories.push_back(isyst->fSubCategory);
                } else {
                  Categories.push_back(isyst->fCategory);
                }
            }
        }
    }
    unsigned int Nsyst = Names.size();

    // read fit resutls
    const NuisParameter* par(nullptr);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        if(fFitList[i_fit]->fFitResultsRootFile!="") fFitList[i_fit]->ReadFitResults(fFitList[i_fit]->fFitResultsRootFile);
        else                                         fFitList[i_fit]->ReadFitResults(dirs[i_fit]+"/Fits/"+names[i_fit]+suffs[i_fit]+".root");
    }

    // exclude unused systematics
    std::vector<string> NamesNew;
    std::vector<string> TitlesNew;
    std::vector<string> CategoriesNew;
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(fCombine && i_fit==N-1) break;
            const FitResults* fitRes = fFitList[i_fit]->fFitResults.get();
            for(const auto& inp : fitRes->GetNuisanceParameters()) {
                par = inp.second.get();
                systName = par->fName;
                if(systName==Names[i_syst]){
                    found = true;
                    break;
                }
            }
            if(found) break;
        }
        if(found){
            if(category=="" || category==Categories[i_syst]){
                NamesNew.push_back(Names[i_syst]);
                TitlesNew.push_back(Titles[i_syst]);
                CategoriesNew.push_back(Categories[i_syst]);
            }
        }
    }
    //
    Nsyst = NamesNew.size();
    Names.clear();
    Titles.clear();
    Categories.clear();
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        Names.push_back(NamesNew[i_syst]);
        Titles.push_back(TitlesNew[i_syst]);
        Categories.push_back(CategoriesNew[i_syst]);
    }
    if(fNuisParListFile!=""){
        //
        // reorder NPs
        Names.clear();
        Titles.clear();
        Categories.clear();
        ifstream in;
        in.open(fNuisParListFile.c_str());
        while(true){
            in >> systName;
            if(!in.good()) break;
            LOG(DEBUG) << "Looking for " << systName << "\n";
            std::string temp_string = "";
            for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
                if(NamesNew[i_syst]==systName){
                    temp_string+=  "found";
                    temp_string += ", title = " + TitlesNew[i_syst];
                    Names.push_back(NamesNew[i_syst]);
                    Titles.push_back(TitlesNew[i_syst]);
                    Categories.push_back(CategoriesNew[i_syst]);
                    break;
                }
            }
            LOG(DEBUG) << temp_string << "\n";
        }
        in.close();
    }
    Nsyst = Names.size();

    // fill stuff
    std::vector< TGraphAsymmErrors> g;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        // create maps for NP's
        std::map<string,double> centralMap;
        std::map<string,double> errUpMap;
        std::map<string,double> errDownMap;
        FitResults *fitRes;
        if(fCombine && i_fit==N-1){
            fitRes = new FitResults();
            if(fFitResultsRootFile!="") fitRes->ReadFromRootFile(fFitResultsRootFile);
            else                        fitRes->ReadFromRootFile(fOutDir+"/Fits/"+fName+fSaveSuf+".root");
        }
        else{
            fitRes = fFitList[i_fit]->fFitResults.get();
        }
        for (const auto& inp : fitRes->GetNuisanceParameters()) {
            par = inp.second.get();
            systName = par->fName;
            centralMap[systName] = par->fFitValue;
            errUpMap[systName]   = par->fPostFitUp;
            errDownMap[systName] = par->fPostFitDown;
        }
        //
        // create the graphs
        g.emplace_back(Nsyst);
        for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
            systName = Names[i_syst];
            if(centralMap[systName]!=0 || (errUpMap[systName]!=0 || errDownMap[systName]!=0)){
                g[i_fit].SetPoint(i_syst,centralMap[systName],(Nsyst-i_syst-1)+0.5+yshift[i_fit]);
                g[i_fit].SetPointEXhigh(i_syst,  errUpMap[systName]);
                g[i_fit].SetPointEXlow( i_syst, -errDownMap[systName]);
            }
            else{
                g[i_fit].SetPoint(i_syst,-10,-10);
                g[i_fit].SetPointEXhigh(i_syst, 0);
                g[i_fit].SetPointEXlow( i_syst, 0);
            }
        }
        if (fCombine && i_fit==N-1) {
            delete fitRes;
        }
    }

    max = Nsyst;

    int lineHeight = 20;
    int offsetUp = 50;
    int offsetDown = 60;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l0(0,0,0,max);
    l0.SetLineStyle(7);
    l0.SetLineColor(kBlack);
    TBox b1(-1,0,1,max);
    TBox b2(-2,0,2,max);
    b1.SetFillColor(kGreen);
    b2.SetFillColor(kYellow);
    b2.Draw("same");
    b1.Draw("same");
    l0.Draw("same");

    int i_style = 0;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if (!(fCombine && i_fit==N-1)){
            if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        g.at(i_fit).SetLineColor(color.at(i_style % color.size()));
        g.at(i_fit).SetMarkerColor(color.at(i_style % color.size()));
        g.at(i_fit).SetMarkerStyle(style.at(i_style % style.size()));
        g.at(i_fit).Draw("P same");
        i_style ++;
    }

    TLatex systs{};
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(unsigned int i_syst=0;i_syst<Nsyst;i_syst++){
        systs.DrawLatex(3.,(Nsyst-i_syst-1)+0.25,Titles[i_syst].c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    h_dummy.GetXaxis()->CenterTitle();
    h_dummy.GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    h_dummy.GetXaxis()->SetTitleOffset(1.2);

    float topMargin = c.GetTopMargin();
    TLegend leg(0.30, 1 - 1.1 * topMargin, 0.85, 0.99);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetNColumns(N);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if (!(fCombine && i_fit==N-1)){
            if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        leg.AddEntry(&g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg.Draw();

    if (fFitList[0]->fPlotLabel != "none") {
        TRExLabelNew(0.05, 1 - (4. / 5) * topMargin, TRExFitter::EXPERIMENT_LABEL.c_str(), fFitList[0]->fPlotLabel.c_str(), kBlack, gStyle->GetTextSize(), 0.10);
    }

    gPad->RedrawAxis();

    if(category=="") Common::SaveCanvasAs(c, fOutDir+"/NuisPar_comp"+fSaveSuf);
    else             Common::SaveCanvasAs(c, fOutDir+"/NuisPar_comp"+fSaveSuf+"_"+category);
}

//__________________________________________________________________________________
//
void MultiFit::CompareNormFactors(const std::string& category) const{
    double ydist = 0.2;

    // Fit titles
    vector<string> dirs;
    vector<string> names;
    vector<string> suffs;
    vector<string> titles;
    vector<double>  yshift;

    static const std::vector<int> color = {kBlack,kRed,kBlue,kViolet,kOrange};
    static const std::vector<int> style = {kFullCircle,kOpenCircle,kFullTriangleUp,kOpenTriangleDown,kOpenDiamond};

    unsigned int N = fFitList.size();
    if(fCombine) N++;

    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1){
            LOG(INFO) << "Adding Combined Fit\n";
            dirs.push_back( fOutDir );
            names.push_back( fName );
            titles.push_back(fCombiLabel);
            suffs.push_back( "" );
        }
        else{
            dirs.push_back( fFitList[i_fit]->fName );
            names.push_back( fFitList[i_fit]->fInputName );
            titles.push_back( fFitLabels[i_fit] );
            suffs.push_back( fFitSuffs[i_fit] );
        }
        yshift.push_back( 0. - ydist*N/2. + ydist*(N-i_fit-1) );
    }

    double xmin = -1.;
    double xmax = 10.;
    double max = 0.;

    // create a list of Norm Factors
    std::vector< string > Names;
    std::vector< string > Titles;
    std::vector< string > Categories;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        for(const auto& inorm : fFitList[i_fit]->fNormFactors) {
            const std::string normName = inorm->fName;
            if (Common::FindInStringVector(fPOIs,normName)>=0) continue;
            if (std::find(Names.begin(), Names.end(), normName) != Names.end()) continue;
            Names.push_back(normName);
            Titles.push_back(inorm->fTitle);
            Categories.push_back(inorm->fCategory);
        }
    }
    unsigned int Nnorm = Names.size();

    // read fit resutls
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(fCombine && i_fit==N-1) break;
        fFitList[i_fit]->ReadFitResults(dirs[i_fit]+"/Fits/"+names[i_fit]+suffs[i_fit]+".root");
    }

    // exclude norm factors
    std::vector<string> NamesNew;
    std::vector<string> TitlesNew;
    std::vector<string> CategoriesNew;
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        if (std::find(fBlindedParams.begin(), fBlindedParams.end(), Names[i_norm]) != fBlindedParams.end()) {
            LOG(INFO) << "Parameter " << Names[i_norm] << " is blinded, ignoring it.\n";
            break;
        }
        bool found = false;
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(fCombine && i_fit==N-1) break;
            const FitResults* fitRes = fFitList[i_fit]->fFitResults.get();
            for(const auto& inp : fitRes->GetNuisanceParameters()) {
                const NuisParameter* par = inp.second.get();
                const std::string& normName = par->fName;
                if(normName==Names[i_norm]){
                    found = true;
                    break;
                }
            }
            if(found) break;
        }
        if(found){
            if(category=="" || category==Categories[i_norm]){
                NamesNew.push_back(Names[i_norm]);
                TitlesNew.push_back(Titles[i_norm]);
                CategoriesNew.push_back(Categories[i_norm]);
            }
        }
    }
    //
    Nnorm = NamesNew.size();
    Names.clear();
    Titles.clear();
    Categories.clear();
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        Names.push_back(NamesNew[i_norm]);
        Titles.push_back(TitlesNew[i_norm]);
        Categories.push_back(CategoriesNew[i_norm]);
    }
    if(fNuisParListFile!=""){
        //
        // reorder NPs
        Names.clear();
        Titles.clear();
        Categories.clear();
        ifstream in;
        in.open(fNuisParListFile.c_str());
        std::string normName;
        while(true){
            in >> normName;
            if(!in.good()) break;
            LOG(DEBUG) << "Looking for " << normName << "... \n";
            std::string temp_string = "";
            for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
                if(NamesNew[i_norm]==normName){
                    temp_string+= "found";
                    temp_string+= ", title = " + TitlesNew[i_norm];
                    Names.push_back(NamesNew[i_norm]);
                    Titles.push_back(TitlesNew[i_norm]);
                    Categories.push_back(CategoriesNew[i_norm]);
                    break;
                }
            }
            LOG(DEBUG) << temp_string << "\n";
        }
        in.close();
    }
    Nnorm = Names.size();

    // fill stuff
    std::vector< TGraphAsymmErrors > g;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        // create maps for NP's
        std::map<string,double> centralMap;
        std::map<string,double> errUpMap;
        std::map<string,double> errDownMap;
        FitResults *fitRes;
        if(fCombine && i_fit==N-1){
            fitRes = new FitResults();
            if(fFitResultsRootFile!="") fitRes->ReadFromRootFile(fFitResultsRootFile);
            else                        fitRes->ReadFromRootFile(fOutDir+"/Fits/"+fName+fSaveSuf+".root");
        }
        else{
            fitRes = fFitList[i_fit]->fFitResults.get();
        }
        for(const auto& inp : fitRes->GetNuisanceParameters()) {
            const NuisParameter* par = inp.second.get();
            const std::string& normName = par->fName;
            centralMap[normName] = par->fFitValue;
            errUpMap[normName]   = par->fPostFitUp;
            errDownMap[normName] = par->fPostFitDown;
        }
        //
        // create the graphs
        g.emplace_back(Nnorm);
        for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
            const std::string& normName = Names[i_norm];
            if(centralMap[normName]!=0 || (errUpMap[normName]!=0 || errDownMap[normName]!=0)){
                g[i_fit].SetPoint(i_norm,centralMap[normName],(Nnorm-i_norm-1)+0.5+yshift[i_fit]);
                g[i_fit].SetPointEXhigh(i_norm,  errUpMap[normName]);
                g[i_fit].SetPointEXlow( i_norm, -errDownMap[normName]);
            }
            else{
                g[i_fit].SetPoint(i_norm,-10,-10);
                g[i_fit].SetPointEXhigh(i_norm, 0);
                g[i_fit].SetPointEXlow( i_norm, 0);
            }
        }
        if(fCombine && i_fit==N-1){
            delete fitRes;
        }
    }

    max = Nnorm;

    int lineHeight = 50;
    int offsetUp = 50;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.2/(8./6.));
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l1(1,0,1,max);
    l1.SetLineColor(kGray);
    l1.Draw("same");

    int i_style = 0;
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(!(fCombine && i_fit==N-1)){
            if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        g[i_fit].SetLineColor(color[i_style]);
        g[i_fit].SetMarkerColor(color[i_style]);
        g[i_fit].SetMarkerStyle(style[i_style]);
        g[i_fit].Draw("P same");
        i_style++;
    }

    TLatex norms{};
    norms.SetTextSize( norms.GetTextSize()*0.8 );
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        norms.DrawLatex(xmin-0.15*(xmax-xmin),(Nnorm-i_norm-1)+0.25,Titles[i_norm].c_str());
    }
    TLatex values{};
    values.SetTextSize( values.GetTextSize()*0.8 );
    for(unsigned int i_norm=0;i_norm<Nnorm;i_norm++){
        for(unsigned int i_fit=0;i_fit<N;i_fit++){
            if(!(fCombine && i_fit==N-1)){
                if (!fFitList.at(i_fit)->fUseInComparison) continue;
            }
            values.DrawLatex((xmin+(xmax-xmin)*(0.45+i_fit*0.55/N)),(Nnorm-i_norm-1)+0.25,
                             Form("%.2f^{+%.2f}_{-%.2f}",g[i_fit].GetX()[i_norm],g[i_fit].GetErrorXhigh(i_norm),g[i_fit].GetErrorXlow(i_norm)));
        }
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );

    TLegend leg(0.45,1.-0.02*(30./max),0.98,0.99);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetNColumns(N);
    for(unsigned int i_fit=0;i_fit<N;i_fit++){
        if(!(fCombine && i_fit==N-1)){
         if (!fFitList.at(i_fit)->fUseInComparison) continue;
        }
        leg.AddEntry(&g[i_fit],titles[i_fit].c_str(),"lp");
    }
    leg.Draw();

    if (fFitList[0]->fPlotLabel != "none") {
        float topMargin = c.GetTopMargin();
        TRExLabelNew(0.15, 1 - (2. / 3) * topMargin, TRExFitter::EXPERIMENT_LABEL.c_str(), fFitList[0]->fPlotLabel.c_str(), kBlack, gStyle->GetTextSize(), 0.10);
    }

    gPad->RedrawAxis();

    if(category=="") Common::SaveCanvasAs(c, fOutDir+"/NormFactors_comp"+fSaveSuf);
    else             Common::SaveCanvasAs(c, fOutDir+"/NormFactors_comp"+fSaveSuf+"_"+category);
}

//__________________________________________________________________________________
//
void MultiFit::PlotCombinedCorrelationMatrix() const{
    TRExFit *fit = fFitList.at(0).get();
    if(fit->fStatOnly){
        LOG(INFO) << "Stat only fit => No Correlation Matrix generated.\n";
        return;
    }
    //plot the correlation matrix (considering only correlations larger than TRExFitter::CORRELATIONTHRESHOLD)
    fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".root");
    if(fit->fFitResults){
        fit->fFitResults->SetOutputFolder(fOutDir);
        const std::string path = fOutDir+"/CorrMatrix_comb"+fSaveSuf;
        fit->fFitResults->DrawCorrelationMatrix(path,
                                                fHEPDataFormat,
                                                TRExFitter::CORRELATIONTHRESHOLD, false,
                                                fPOIs);
    }
}

//____________________________________________________________________________________
//
void MultiFit::ProduceNPRanking(const std::string& NPnames) const {
    LOG(INFO) << "....................................\n";
    LOG(INFO) << "Producing Ranking...\n";

    if(fFitType==2){
        LOG(ERROR) << "For ranking plots, the SPLUSB FitType is needed.\n";
        exit(EXIT_FAILURE);
    }

    const std::string& inputData = fDataName;

    // create a list of norm factors
    std::vector<std::shared_ptr<NormFactor> > vNormFactors = GetFitNormFactors();

    const std::string frFile = fOutDir + "/Fits/"+fName+fSaveSuf+".root";

    xRooNode rooNode = Common::ReadWSandFitResults(fOutDir+"/ws_combined"+fSaveSuf+".root",
                                                   frFile);

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf for postfit\n";
        exit(EXIT_FAILURE);
    }

    xRooNode rooNodePrefit = Common::ReadWSandFitResults(fOutDir+"/ws_combined"+fSaveSuf+".root",
                                                         "");

    auto pdfPrefit = rooNodePrefit["simPdf"];
    if (!pdfPrefit) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf for prefit\n";
        exit(EXIT_FAILURE);
    }

    //
    // List of systematics to check
    //
    RankingManager manager{};
    manager.SetUseHesseBeforeMigrad(fUseHesseBeforeMigrad);
    manager.SetRandomNPs(fRndRange, fUseRnd, fRndSeed);
    manager.SetMaximumFCNcalls(fMaximumNumberFCNcalls);
    manager.SetToleranceScale(fToleranceScale);
    manager.SetShiftGlobalObservables(fShiftGlobalObservablesInRanking);

    for (const auto& iparameter : pdf->floats()) {
        const std::string name = iparameter->GetName();
        if (NPnames != "all" && name.find(NPnames) == std::string::npos) continue;
                // check if it is NF
        auto itr = std::find_if(vNormFactors.begin(), vNormFactors.end(), [&name](const auto& nf){return name == nf->fName;});
        if (itr != vNormFactors.end()) {
            if (std::find(fPOIs.begin(), fPOIs.end(), name) != fPOIs.end()) {
                if (fUsePOISinRanking) {
                    manager.AddNuisPar(name, true);
                }
            }
            manager.AddNuisPar(name, true);
        } else {
            manager.AddNuisPar(name, false);
        }
    }

    //
    // Text files containing information necessary for drawing of ranking plot
    //
    std::string outName = fOutDir+"/Fits/NPRanking"+fSaveSuf;
    if(NPnames!="all") outName += "_"+NPnames;

    manager.SetOutputPath(outName);

    //
    // Get the combined model
    //
    std::unique_ptr<TFile> f(TFile::Open((fOutDir+"/ws_combined"+fSaveSuf+".root").c_str()));
    if (!f) {
        LOG(ERROR) << "Cannot open file!\n";
        exit(EXIT_FAILURE);
    }
    RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f->Get("combWS"));
    if (!ws) {
        LOG(ERROR) << "Cannot read workspace!\n";
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

    //
    // Creates the data object
    //
    RooDataSet* data = nullptr;
    if(inputData=="asimovData"){
        RooArgSet empty;// = RooArgSet();
        data = static_cast<RooDataSet*>(RooStats::AsymptoticCalculator::MakeAsimovData( (*mc), RooArgSet(ws->allVars()), static_cast<RooArgSet&>(empty)));
    }
    else if(inputData!=""){
        data = static_cast<RooDataSet*>(ws->data(inputData.c_str()));
    } else {
        LOG(WARNING) << "You didn't specify inputData!\n";
        exit(EXIT_FAILURE);
    }

    //
    // Create snapshot to keep inital values
    //
    ws -> saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));

    // find the last config that is used in a fit
    int validConfig(999999);
    for (int i = fFitList.size() - 1; i >= 0; --i) {
        if (fFitList.at(i)->fUseInFit) {
            validConfig = i;
            break;
        }
    }

    if (validConfig == 999999) {
        LOG(ERROR) << "Did not find a single config that should be used for a fit\n";
        exit(EXIT_FAILURE);
    }

    TRExFit *fit = fFitList.at(validConfig).get();
    fit->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".root");

    manager.SetInjectGlobalObservables(fit->fInjectGlobalObservables);
    manager.SetNPValues(fit->fFitNPValues);
    manager.SetFixedNPs(GetFixedNPs());
    manager.SetFitStrategy(fFitStrategy);
    if(fPOIs.size() > 0){
        manager.SetPOINames(fPOIs);
    }
    else{
        LOG(WARNING) << "No POI set. Not able to produce ranking.\n";
        return;
    }
    manager.SetNCPU(fCPU);
    manager.SetStatOnly(fStatOnly);
    manager.SetUsePOISinRanking(fUsePOISinRanking);
    manager.RunRanking(pdf, pdfPrefit, ws, data, GetFitNormFactors());

    f->Close();
}

//____________________________________________________________________________________
//
void MultiFit::PlotNPRankingManager() const{
  if(fFitList[0]->fRankingPlot=="MERGE"  || fFitList[0]->fRankingPlot=="ALL") PlotNPRanking(true,true);
  if(fFitList[0]->fRankingPlot=="SYSTS"  || fFitList[0]->fRankingPlot=="ALL") PlotNPRanking(true,false);
  if(fFitList[0]->fRankingPlot=="GAMMAS" || fFitList[0]->fRankingPlot=="ALL") PlotNPRanking(false,true);
}

//____________________________________________________________________________________
//
void MultiFit::PlotNPRanking(bool flagSysts, bool flagGammas) const {
    LOG(INFO) << "..................................../n";
    LOG(INFO) << "Plotting Ranking...\n";
    //
    for (const auto& poi : fPOIs) {
        const std::string fileToRead = fOutDir+"/Fits/NPRanking"+fSaveSuf+"_"+poi+".txt";
        std::ifstream in(fileToRead.c_str());
        if (!in.good()) { // file doesnt exist
            const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fOutDir+"/Fits/","NPRanking" + fSaveSuf + "_", poi);
            Common::MergeTxTFiles(inPaths, fileToRead);
        }
    }

    RankingManager manager{};
    manager.SetUseHesseBeforeMigrad(fUseHesseBeforeMigrad);
    manager.SetPlotLabel(fFitList[0]->fPlotLabel);
    manager.SetLumiLabel(fFitList[0]->fLumiLabel);
    manager.SetCmeLabel(fFitList[0]->fCmeLabel);
    manager.SetUseHEPDataFormat(fHEPDataFormat);
    manager.SetName(fOutDir);
    manager.SetSuffix(fSaveSuf);
    manager.SetMaxNPPlot(fFitList[0]->fRankingMaxNP);
    manager.SetRankingCanvasSize(fFitList[0]->fNPRankingCanvasSize);
    manager.SetShiftGlobalObservables(fShiftGlobalObservablesInRanking);

    std::string extra("");
    if (!flagSysts || !flagGammas) {
        extra = !flagSysts ? "_gammas" : "_systs";
    }

    std::vector<std::string> nfs;
    auto tmp = this->GetFitNormFactors();
    for (const auto& inf : tmp) {
        nfs.emplace_back(inf->fName);
    }

    std::vector<std::string> sfs;
    auto tmpSF = this->GetFitShapeFactors();
    for (const auto& isf : tmpSF) {
        sfs.emplace_back(isf->fName);
    }

    for (std::size_t ipoi = 0; ipoi < fPOIs.size(); ++ipoi) {
        const std::string fileToRead = fOutDir+"/Fits/NPRanking"+fSaveSuf+"_"+fPOIs.at(ipoi)+".txt";
        manager.SetOutputPath(fileToRead);
        manager.SetSuffix(fSaveSuf+"_"+fPOIs.at(ipoi)+extra);
        manager.SetRankingPOIName(fRankingPOIName.at(ipoi));
        manager.SetUpperAxisNdivisions(fRankingUpperAxisNdivision.at(ipoi));
        manager.SetRankingPOIAxisScale(fRankingPOIAxisScale.at(ipoi));
        manager.ReadRankingResults(flagSysts, flagGammas);
        manager.PlotRanking(GetFitRegions(), nfs, sfs, flagSysts, flagGammas);
        if (fShiftGlobalObservablesInRanking) {
            manager.SetSuffix(fSaveSuf+"_"+fPOIs.at(ipoi));
            manager.ProduceGroupedImpact(this->GetSubcategoryMap());
        }
    }
}

//__________________________________________________________________________________
//
void MultiFit::GetLikelihoodScan(const RooWorkspace *ws, const std::string& varName, RooDataSet* data,bool recreate) const{
    LOG(INFO) << "Running likelihood scan for the parameter = " << varName << "\n";
    double minVal;
    double maxVal;
    gSystem->mkdir((fOutDir+"/LHoodPlots").c_str());

    std::unique_ptr<TGraph> graph;
    if(recreate){
        LikelihoodScanManager manager{};
        manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps, fLHscanStep);
        manager.SetNCPU(fCPU);
        manager.SetOffSet(true);
        manager.SetUseNll(fUseNllInLHscan);
        manager.SetToleranceScale(fToleranceScale);
        manager.SetUseAutoDiff(fUseAutoDiff);
        manager.SetNLLOffset(fNLLOffset);

        FitUtils::FittingOptions options;
        options.noGammas = fStatOnly;
        options.noSystematics = fStatOnly;
        options.randomize = fUseRnd;
        options.randomValue = fRndRange;
        auto fixedNPs = this->GetFixedNPs();
        std::vector<std::string> constNPs;
        std::vector<double> constNPvalues;
        for (const auto& iparam : fixedNPs) {
            constNPs.emplace_back(iparam.first);
            constNPvalues.emplace_back(iparam.second);
        }
        options.constNPs = std::move(constNPs);
        options.constNPvalues = std::move(constNPvalues);
        auto nfs = this->GetFitNormFactors();
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inf : nfs) {
            if (Common::FindInStringVector(fPOIs,inf->fName)>=0) {
                continue;
            }
            npNames. emplace_back(inf->fName);
            npValues.emplace_back(inf->GetNominal());
        }
        options.initialNPs = std::move(npNames);
        options.initialNPvalues = std::move(npValues);
        options.constPOI = fFitType == TRExFit::FitType::BONLY;
        std::vector<std::pair<std::string, double> > poiVals;
        std::vector<std::string> constNFs;
        for (const auto& inf : nfs) {
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

        if (x.empty() || y.empty()) {
            LOG(WARNING) << "Skipping LHscan\n";
            return;
        }

        graph = std::make_unique<TGraph>(x.size(), &x[0], &y[0]);
        minVal = x[0];
        maxVal = x.back();

        YamlConverter converter{};
        converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
        converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
        converter.WriteLikelihoodScan(scanResult, fOutDir+"/LHoodPlots/NLLscan_"+varName+".yaml");
        if (fHEPDataFormat) {
            converter.WriteLikelihoodScanHEPData(scanResult, fOutDir, varName);
        }

    }
    else{
        TFile* f = TFile::Open((fOutDir+"/LHoodPlots/NLLscan_"+varName+"_curve.root").c_str(), "READ");
        graph.reset(dynamic_cast<TGraph*>(f->Get("LHscan")));
        graph->SetLineColor(kRed);
        graph->SetLineWidth(3);
        minVal = graph->GetXaxis()->GetXmin();
        maxVal = graph->GetXaxis()->GetXmax();
    }

    TCanvas can("NLLscan");
    can.SetTopMargin(0.1);
    TString cname="";
    cname.Append("NLLscan_");
    cname.Append(varName);

    can.SetTitle(cname);
    can.SetName(cname);
    can.cd();
    graph->Draw("ALP");

    std::size_t nValid(0);
    for (const auto& ifit : fFitList) {
        if (ifit->fUseInFit) ++nValid;
    }

    // take the LH curves also for other fits
    // cppcheck-suppress variableScope
    std::vector<TGraph*> curve_fit;
    TLegend leg(0.5,0.85-0.06*(nValid+1),0.75,0.85);
    leg.SetFillColor(kWhite);
    leg.SetBorderSize(0);
    leg.SetTextSize(gStyle->GetTextSize());
    leg.SetTextFont(gStyle->GetTextFont());
    if(fCompare){
        for(auto& fit : fFitList){
            if (!fit->fUseInFit) continue;
            TFile* f = nullptr;
            if(fit->fFitResultsRootFile!=""){
                std::vector<std::string> v = Common::Vectorize(fit->fFitResultsRootFile,'/');
                f = TFile::Open((v[0]+"/LHoodPlots/NLLscan_"+varName+"_curve.root").c_str());
            }
            else{
                f = TFile::Open((fit->fName+"/LHoodPlots/NLLscan_"+varName+"_curve.root").c_str());
            }
            if(f!=nullptr) {
                curve_fit.push_back(dynamic_cast<TGraph*>(f->Get("LHscan")));
            }
            else {
                curve_fit.push_back(nullptr);
            }
        }
        //
        int idx = 0;
        for(auto crv : curve_fit){
            if (!crv) continue;
            if(idx==0) crv->SetLineColor(kBlue);
            if(idx==1) crv->SetLineColor(kGreen);
            crv->Draw("LP same");
            leg.AddEntry(crv,fFitList[idx]->fLabel.c_str(),"l");
            idx++;
        }
        leg.AddEntry(graph.get(),"Combined","l");
    }

    graph->GetXaxis()->SetRangeUser(minVal,maxVal);

    // y axis
    graph->GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");
    if(TRExFitter::SYSTMAP[varName]!="") graph->GetXaxis()->SetTitle(TRExFitter::SYSTMAP[varName].c_str());
    else if(TRExFitter::NPMAP[varName]!="") graph->GetXaxis()->SetTitle(TRExFitter::NPMAP[varName].c_str());


    TLatex tex{};
    tex.SetTextColor(kGray+2);

    TLine l1s(minVal,0.5,maxVal,0.5);
    l1s.SetLineStyle(kDashed);
    l1s.SetLineColor(kGray);
    l1s.SetLineWidth(2);
    if(graph->GetMaximum()>2){
        l1s.Draw();
        tex.DrawLatex(maxVal,0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
    }

    if(graph->GetMaximum()>2){
        TLine l2s(minVal,2,maxVal,2);
        l2s.SetLineStyle(kDashed);
        l2s.SetLineColor(kGray);
        l2s.SetLineWidth(2);
        l2s.Draw();
        tex.DrawLatex(maxVal,2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
    }
    //
    if(graph->GetMaximum()>4.5){
        TLine l3s(minVal,4.5,maxVal,4.5);
        l3s.SetLineStyle(kDashed);
        l3s.SetLineColor(kGray);
        l3s.SetLineWidth(2);
        l3s.Draw();
        tex.DrawLatex(maxVal,4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
    }
    //
    TLine lv0(0,graph->GetMinimum(),0,graph->GetMaximum());
    lv0.Draw();
    //
    TLine lh0(minVal,0,maxVal,0);
    lh0.Draw();

    if(fCompare){
        leg.Draw();
        if (fFitList[0]->fPlotLabel != "none") TRExLabel(0.15,0.93, TRExFitter::EXPERIMENT_LABEL.c_str(), fFitList[0]->fPlotLabel.c_str(),kBlack);
        myText(0.68,0.93,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
        if(fLabel!="") myText(0.2,0.85,kBlack,Form("#kern[-1]{%s}",fLabel.c_str()));
    }

    can.RedrawAxis();

    Common::SaveCanvasAs(can, fOutDir+"/LHoodPlots/NLLscan_"+varName);

    if(recreate){
        // write it to a ROOT file as well
        TString scanStep="";
        if (fLHscanStep != -1) scanStep = "Step" + TString::Itoa(fLHscanStep, 10);
        std::unique_ptr<TFile> f(TFile::Open(fOutDir+"/LHoodPlots/NLLscan_"+varName+"_"+scanStep+"_curve.root","UPDATE"));
        f->cd();
        graph->Write("LHscan",TObject::kOverwrite);
        f->Close();
    }
}

//____________________________________________________________________________________
//
void MultiFit::Get2DLikelihoodScan(const RooWorkspace *ws, const std::vector<std::string>& varNames, RooDataSet* data) const{
    if (varNames.size() != 2){
        LOG(ERROR) << "Wrong number of parameters provided for 2D likelihood scan, returning\n";
        return;
    }

    LikelihoodScanManager manager{};
    manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps, fLHscanStep);
    manager.SetScanParamsY(fLHscanMinY, fLHscanMaxY, fLHscanStepsY, fLHscanStepY);
    manager.SetNCPU(fCPU);
    manager.SetOffSet( (fLHscanStep == -1 && fLHscanStepY == -1) );
    manager.SetUseNll(fUseNllInLHscan);
    manager.SetToleranceScale(fToleranceScale);
    manager.SetUseAutoDiff(fUseAutoDiff);
    manager.SetNLLOffset(fNLLOffset);
    FitUtils::FittingOptions options;
    options.noGammas = fStatOnly;
    options.noSystematics = fStatOnly;
    options.randomize = fUseRnd;
    options.randomValue = fRndRange;
    auto fixedNPs = this->GetFixedNPs();
    std::vector<std::string> constNPs;
    std::vector<double> constNPvalues;
    for (const auto& iparam : fixedNPs) {
        constNPs.emplace_back(iparam.first);
        constNPvalues.emplace_back(iparam.second);
    }
    options.constNPs = std::move(constNPs);
    options.constNPvalues = std::move(constNPvalues);
    auto nfs = this->GetFitNormFactors();
    std::vector<std::string> npNames;
    std::vector<double> npValues;
    for(const auto& inf : nfs) {
        if (Common::FindInStringVector(fPOIs,inf->fName)>=0) {
            continue;
        }
        npNames. emplace_back(inf->fName);
        npValues.emplace_back(inf->GetNominal());
    }
    options.initialNPs = std::move(npNames);
    options.initialNPvalues = std::move(npValues);
    options.constPOI = fFitType == TRExFit::FitType::BONLY;
    std::vector<std::pair<std::string, double> > poiVals;
    std::vector<std::string> constNFs;
    for (const auto& inf : nfs) {
        poiVals.emplace_back(std::make_pair(inf->fName, inf->GetNominal()));
        if (inf->fConst) {
            constNFs.emplace_back(inf->fName);
        }
    }
    options.poiVals = std::move(poiVals);
    options.constNFs = std::move(constNFs);
    manager.SetFittingOptions(options);

    const std::pair<std::string, std::string> names = std::make_pair(varNames.at(0), varNames.at(1));
    auto scanResult = manager.Run2DScan(ws, names, data);

    const std::vector<double> x = scanResult.x;
    const std::vector<double> y = scanResult.y;
    std::vector<std::vector<double> > z = scanResult.z;

    if (x.empty() || y.empty() || z.empty()) {
        LOG(WARNING) << "Skipping LHscan\n";
        return;
    }

    // find minimium in z
    auto findMin = [](const std::vector<std::vector<double> >& vec) {
        double min(9999999);
        for (const auto& i : vec) {
            for (const auto& j : i) {
                if (j < min) min = j;
            }
        }
        return min;
    };

    const double zmin = findMin(z);

    gSystem->mkdir((fOutDir+"/LHoodPlots").c_str());

    // make plots
    TCanvas can("2D_NLLscan");
    can.cd();

    TGraph2D graph(fLHscanSteps * fLHscanStepsY);
    TH2D h_nll("NLL", "NLL", fLHscanSteps, x[0], x.back(), fLHscanStepsY, y[0], y.back());

    // Only run minimisation if not in parallel
    for (int ipoint = 0; ipoint < fLHscanSteps; ++ipoint) {
        for (int jpoint = 0; jpoint < fLHscanStepsY; ++jpoint) {
              if ((fLHscanStep == -1 && fLHscanStepY == -1)) {
                  // shift the likelihood values to zero (unless single point)
                  z[ipoint][jpoint] -= zmin;
              }
              h_nll.SetBinContent(ipoint+1, jpoint+1, z[ipoint][jpoint]);
              unsigned int i = ipoint * fLHscanStepsY + jpoint;
              graph.SetPoint(i,x[ipoint],y[jpoint],z[ipoint][jpoint]);
        }
    }

    // Create output plot only if not running parallel processing
    if ( (fLHscanStep == -1 && fLHscanStepY == -1) ) {
        // Only draw and save graph when not running parallel
        gStyle->SetPalette(57); // Reset Palette to default (Pruning or Correlation matrinx changes this)
        graph.Draw("colz");
        graph.GetXaxis()->SetRangeUser(x[0],x.back());
        graph.GetYaxis()->SetRangeUser(y[0],y.back());

        // y axis
        graph.GetXaxis()->SetTitle(varNames.at(0).c_str());
        graph.GetYaxis()->SetTitle(varNames.at(1).c_str());

        // Print the canvas
        Common::SaveCanvasAs(can, fOutDir+"/LHoodPlots/NLLscan_"+varNames.at(0)+"_"+varNames.at(1));

        // write it to a ROOT file as well
        std::unique_ptr<TFile> f(TFile::Open((fOutDir+"/LHoodPlots/NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+"_curve.root").c_str(),"UPDATE"));
        f->cd();
        graph.Write(("LHscan_2D_"+varNames.at(0)+"_"+varNames.at(1)).c_str(),TObject::kOverwrite);
        f->Close();
    }

    // Write histogram to Root file as well
    TString scanStep="";
    if (fLHscanStep  != -1) {
        scanStep += "_StepX" + TString::Itoa(fLHscanStep, 10);
    }
    if (fLHscanStepY != -1) {
        scanStep += "_StepY" + TString::Itoa(fLHscanStepY, 10);
    }

    YamlConverter converter{};
    converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
    converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
    converter.WriteLikelihoodScan2D(scanResult, fOutDir+"/LHoodPlots/NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+".yaml");
    if (fHEPDataFormat) {
        converter.WriteLikelihoodScan2DHEPData(scanResult, fOutDir, varNames.at(0) + "_" + varNames.at(1));
    }

    std::unique_ptr<TFile> f2(TFile::Open(fOutDir+"/LHoodPlots/"+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+scanStep+"_histo.root","UPDATE"));
    h_nll.Write("NLL",TObject::kOverwrite);
    f2->Close();
}

//____________________________________________________________________________________
//
void MultiFit::PlotSummarySoverB() const {
    LOG(INFO) << "....................................\n";
    LOG(INFO) << "Producing S/B plot...\n";

    bool includeBonly = false;
    if(fBonlySuffix!="") includeBonly = true;

    const std::string fitFileName = fOutDir+"/Fits/"+fName+fSaveSuf+".root";
    double muFit(0.);
    if (std::filesystem::is_regular_file(fitFileName)) {
        fFitList[0]->ReadFitResults(fitFileName);
        muFit = fFitList[0]->fFitResults->GetNuisParValue(fPOIs.at(0));
    } else {
        muFit = 1;
    }
    const std::string name = fLimitOutputPrefixName;
    double muLimit(0);
    bool hasLimit(false);
    const std::string limitFileName = fOutDir+"/Limits/Asymptotics/"+name+fSaveSuf+".root";
    if (std::filesystem::is_regular_file(limitFileName)) {
        LOG(INFO) << "Reading file " << limitFileName << "\n";
        std::unique_ptr<TFile> f(TFile::Open((limitFileName).c_str()));
        TTree *tree = dynamic_cast<TTree*>(f->Get("stats"));
        if (!tree) {
            LOG(WARNING) << "Cannot read tree from the limit file\n";
            return;
        }
        if (tree->GetEntries() == 0) {
            LOG(WARNING) << "The tree has 0 entries\n";
            return;
        }

        LimitEvent event(tree);

        event.GetEntry(0);
        muLimit = event.obs_upperlimit;
        f->Close();
        hasLimit = true;
    }

    if (hasLimit) {
        LOG(INFO) << "Limit found, will add it\n";
    } else {
        LOG(INFO) << "Limit not found, will proceed without it\n";
    }

    std::vector<string> fileNames;
    std::vector<string> fileNamesBonly;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(unsigned int i_reg=0;i_reg<ifit->fRegions.size();i_reg++){
            if(ifit->fRegions[i_reg]->fRegionType==Region::VALIDATION) continue;
            if (!fOnlyRegions.empty() && Common::FindInStringVector(fOnlyRegions, ifit->fRegions[i_reg]->fName) < 0) continue;
            fileNames.emplace_back(ifit->fName+"/Histograms/"+ifit->fRegions[i_reg]->fName+"_postFit.root");
            if(includeBonly)
                fileNamesBonly.emplace_back(ifit->fName+"/Histograms/"+ifit->fRegions[i_reg]->fName+fBonlySuffix+"_postFit.root");
        }
    }
    const unsigned int Nhist = fileNames.size();

    //
    // create a list of all the samples
    std::vector<string> sigList;
    std::vector<string> bkgList;
    std::vector<string> dataList;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(unsigned int i_smp=0;i_smp<ifit->fSamples.size();i_smp++){
            if(ifit->fSamples[i_smp]->fType==Sample::SampleType::SIGNAL && Common::FindInStringVector(sigList,ifit->fSamples[i_smp]->fName)<0)
                sigList.push_back(ifit->fSamples[i_smp]->fName);
            if(ifit->fSamples[i_smp]->fType==Sample::SampleType::BACKGROUND && Common::FindInStringVector(bkgList,ifit->fSamples[i_smp]->fName)<0)
                bkgList.push_back(ifit->fSamples[i_smp]->fName);
            if(ifit->fSamples[i_smp]->fType==Sample::SampleType::DATA && Common::FindInStringVector(dataList,ifit->fSamples[i_smp]->fName)<0)
                dataList.push_back(ifit->fSamples[i_smp]->fName);
        }
    }
    //
    // create a list of all the systematics
    std::vector<string> systList;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        // actual systematics
        for(unsigned int i_syst=0;i_syst<ifit->fSystematics.size();i_syst++){
            if(Common::FindInStringVector(systList,ifit->fSystematics[i_syst]->fName)<0)
                systList.push_back(ifit->fSystematics[i_syst]->fName);
        }
        // norm factors
        for(unsigned int i_norm=0;i_norm<ifit->fNormFactors.size();i_norm++){
            if(Common::FindInStringVector(systList,ifit->fNormFactors[i_norm]->fName)<0)
                systList.push_back(ifit->fNormFactors[i_norm]->fName);
        }
    }

    unsigned int Nsyst = systList.size();

    std::vector<std::unique_ptr<TFile> > file;
    std::vector<std::unique_ptr<TFile> > fileBonly;
    std::vector<std::unique_ptr<TH1D> > h_sig;
    std::vector<std::unique_ptr<TH1D> > h_bkg;
    std::vector<std::unique_ptr<TH1D> > h_bkgBonly;
    std::vector<std::unique_ptr<TH1D> > h_tot_bkg_prefit;
    std::vector<std::unique_ptr<TH1D> > h_data;
    std::vector<std::vector<std::unique_ptr<TH1D> > > h_syst_up  (Nsyst);
    std::vector<std::vector<std::unique_ptr<TH1D> > > h_syst_down(Nsyst);
    for (auto& ihist : h_syst_up) {
        ihist.resize(Nhist);
    }
    for (auto& ihist : h_syst_down) {
        ihist.resize(Nhist);
    }

    // get histos
    for(unsigned int i_hist = 0; i_hist < Nhist; ++i_hist) {
        LOG(DEBUG) << "Opening file " << fileNames[i_hist] << "\n";
        file.emplace_back(std::move(TFile::Open(fileNames[i_hist].c_str())));
        if(includeBonly){
            LOG(DEBUG) << "Opening file " << fileNamesBonly[i_hist] << "\n";
            fileBonly.emplace_back(std::move(TFile::Open(fileNamesBonly[i_hist].c_str())));
        }
        //
        // initialize null histogram pointers
        h_sig.emplace_back(nullptr);
        h_bkg.emplace_back(nullptr);
        h_bkgBonly.emplace_back(nullptr);
        h_tot_bkg_prefit.emplace_back(nullptr);
        h_data.emplace_back(nullptr);
        for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
            h_syst_up  [i_syst][i_hist] = nullptr;
            h_syst_down[i_syst][i_hist] = nullptr;
        }
        //
        for(unsigned int i_sig=0;i_sig<sigList.size();i_sig++){
            LOG(DEBUG) << "  Getting histogram h_" << sigList[i_sig] << "_postFit\n";
            std::unique_ptr<TH1D> h_tmp(dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+sigList[i_sig]+"_postFit").c_str())));
            if(h_tmp!=nullptr){
                h_tmp->SetDirectory(nullptr);
                LOG(DEBUG) << " ... FOUND\n";
                if(h_sig[i_hist]==nullptr) {
                    h_sig[i_hist].reset(static_cast<TH1D*>(h_tmp->Clone()));
                    h_sig[i_hist]->SetDirectory(nullptr);
                } else {
                    h_sig[i_hist]->Add(h_tmp.get());
                }
            }
        }
        for(unsigned int i_bkg=0;i_bkg<bkgList.size();i_bkg++){
            LOG(DEBUG) << "  Getting histogram h_" << bkgList[i_bkg] << "_postFit\n";
            std::unique_ptr<TH1D> h_tmp(dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_postFit").c_str())));
            std::unique_ptr<TH1D> h_tmpBonly(nullptr);
            if(includeBonly) {
                h_tmpBonly.reset(dynamic_cast<TH1D*>(fileBonly[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_postFit").c_str())));
                h_tmpBonly->SetDirectory(nullptr);
            }
            if(h_tmp!=nullptr){
                h_tmp->SetDirectory(nullptr);
                LOG(DEBUG) << " ... FOUND\n";
                if(h_bkg[i_hist]==nullptr) {
                    h_bkg[i_hist].reset(static_cast<TH1D*>(h_tmp->Clone()));
                    h_bkg[i_hist]->SetDirectory(nullptr);
                } else {
                    h_bkg[i_hist]->Add(h_tmp.get());
                }
            }
            if(h_tmpBonly!=nullptr){
                h_tmpBonly->SetDirectory(nullptr);
                LOG(DEBUG) << " ... B-only FOUND\n";
                if(h_bkgBonly[i_hist]==nullptr) {
                    h_bkgBonly[i_hist].reset(static_cast<TH1D*>(h_tmpBonly->Clone()));
                    h_bkgBonly[i_hist]->SetDirectory(nullptr);
                } else {
                    h_bkgBonly[i_hist]->Add(h_tmpBonly.get());
                }
            }
            // syst variations
            for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
                // up
                std::unique_ptr<TH1D> h_tmpup(dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_"+systList[i_syst]+"_Up_postFit").c_str())));
                if(h_tmpup!=nullptr){
                    h_tmpup->SetDirectory(nullptr);
                    if(h_syst_up[i_syst][i_hist]==nullptr) {
                        h_syst_up[i_syst][i_hist].reset(static_cast<TH1D*>(h_tmpup->Clone()));
                        h_syst_up[i_syst][i_hist]->SetDirectory(nullptr);
                    } else {
                        h_syst_up[i_syst][i_hist]->Add(h_tmpup.get());
                    }
                }
                // down
                std::unique_ptr<TH1D> h_tmpdown(dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+bkgList[i_bkg]+"_"+systList[i_syst]+"_Down_postFit").c_str())));
                if(h_tmpup!=nullptr){
                    h_tmpdown->SetDirectory(nullptr);
                    if(h_syst_down[i_syst][i_hist]==nullptr) {
                        h_syst_down[i_syst][i_hist].reset(static_cast<TH1D*>(h_tmpdown->Clone()));
                        h_syst_down[i_syst][i_hist]->SetDirectory(nullptr);
                    } else {
                        h_syst_down[i_syst][i_hist]->Add(h_tmpdown.get());
                    }
                }
            }
        }
        for(unsigned int i_data=0;i_data<dataList.size();i_data++){
            std::unique_ptr<TH1D> h_tmp(dynamic_cast<TH1D*>(file[i_hist]->Get( ("h_"+dataList[i_data]).c_str())));
            if(h_tmp!=nullptr){
                h_tmp->SetDirectory(nullptr);
                if(h_data[i_hist]==nullptr) {
                    h_data[i_hist].reset(static_cast<TH1D*>(h_tmp->Clone()));
                    h_data[i_hist]->SetDirectory(nullptr);
                } else {
                    h_data[i_hist]->Add(h_tmp.get());
                }
            }
        }

        if(TRExFitter::PREFITONPOSTFIT) {
            h_tot_bkg_prefit[i_hist].reset(dynamic_cast<TH1D*>(file[i_hist]->Get("h_tot_bkg_prefit")));
            h_tot_bkg_prefit[i_hist]->SetDirectory(nullptr);
        }

        //
        // Fix eventually empty histograms
        if(h_sig[i_hist] ==nullptr){
            h_sig[i_hist].reset(static_cast<TH1D*>(h_bkg[i_hist]->Clone(Form("h_sig[%u]", i_hist))));
            h_sig[i_hist]->Scale(0.);
            h_sig[i_hist]->SetDirectory(nullptr);
        }
        if(h_data[i_hist]==nullptr){
            h_data[i_hist].reset(static_cast<TH1D*>(h_bkg[i_hist]->Clone(Form("h_data[%u]",i_hist))));
            h_data[i_hist]->Scale(0.);
            h_data[i_hist]->SetDirectory(nullptr);
        }
        for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
            // up
            if(h_syst_up[i_syst][i_hist]==nullptr){
                h_syst_up[i_syst][i_hist].reset(static_cast<TH1D*>(h_bkg[i_hist]->Clone(Form("h_syst_up[%d][%u]", i_syst,i_hist))));
                h_syst_up[i_syst][i_hist]->Scale(0.);
                h_syst_up[i_syst][i_hist]->SetDirectory(nullptr);
            }
            else{
                h_syst_up[i_syst][i_hist]->Add(h_bkg[i_hist].get(),-1);
            }
            // down
            if(h_syst_down[i_syst][i_hist]==nullptr){
                h_syst_down[i_syst][i_hist].reset(static_cast<TH1D*>(h_bkg[i_hist]->Clone(Form("h_syst_down[%d][%u]", i_syst,i_hist))));
                h_syst_down[i_syst][i_hist]->Scale(0.);
                h_syst_down[i_syst][i_hist]->SetDirectory(nullptr);
            }
            else{
                h_syst_down[i_syst][i_hist]->Add(h_bkg[i_hist].get(),-1);
            }
        }
    }

    // create combined histogram
    std::unique_ptr<TH1D> h_bkg_comb  = Combine(h_bkg);
    std::unique_ptr<TH1D> h_bkgBonly_comb = includeBonly ? Combine(h_bkgBonly) : nullptr;
    std::unique_ptr<TH1D> h_sig_comb  = Combine(h_sig);
    std::unique_ptr<TH1D> h_data_comb = Combine(h_data);
    std::unique_ptr<TH1D> h_tot_bkg_prefit_comb = TRExFitter::PREFITONPOSTFIT ? Combine(h_tot_bkg_prefit) : nullptr;

    std::vector<std::unique_ptr<TH1D> > h_syst_up_comb  (Nsyst);
    std::vector<std::unique_ptr<TH1D> > h_syst_down_comb(Nsyst);
    for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
        h_syst_up_comb  [i_syst] = Combine(h_syst_up  [i_syst]);
        h_syst_down_comb[i_syst] = Combine(h_syst_down[i_syst]);
    }

    std::vector<double> SoverSqrtB;

    for(int i_bin=1;i_bin<=h_bkg_comb->GetNbinsX();i_bin++){
        const double sig = h_sig_comb->GetBinContent(i_bin);
        const double bkg = h_bkg_comb->GetBinContent(i_bin);
        if (bkg > 1e-9) {
            SoverSqrtB.emplace_back(sig/bkg);
        } else {
            LOG(WARNING) << "Bkg < 1e-9, possible division by zero. Returning S/B = 9999.\n";
            SoverSqrtB.emplace_back(9999.);
        }
    }

    std::unique_ptr<TH1D> h_bkg_ord = Rebin(h_bkg_comb.get(),SoverSqrtB,false);
    std::unique_ptr<TH1D> h_bkgBonly_ord(nullptr);
    if(includeBonly) h_bkgBonly_ord = Rebin(h_bkgBonly_comb.get(),SoverSqrtB,false);
    std::unique_ptr<TH1D> h_sig_ord = Rebin(h_sig_comb.get(),SoverSqrtB,false);
    std::unique_ptr<TH1D> h_data_ord = Rebin(h_data_comb.get(),SoverSqrtB);
    std::unique_ptr<TH1D> h_tot_bkg_prefit_ord(nullptr);
    if(TRExFitter::PREFITONPOSTFIT) h_tot_bkg_prefit_ord = Rebin(h_tot_bkg_prefit_comb.get(),SoverSqrtB,false);

    std::vector<std::unique_ptr<TH1D> > h_syst_up_ord  (Nsyst);
    std::vector<std::unique_ptr<TH1D> > h_syst_down_ord(Nsyst);
    for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
        h_syst_up_ord  [i_syst] = Rebin(h_syst_up_comb  [i_syst].get(),SoverSqrtB,false);
        h_syst_down_ord[i_syst] = Rebin(h_syst_down_comb[i_syst].get(),SoverSqrtB,false);
    }

    for(int i_bin=0;i_bin<h_bkg_ord->GetNbinsX()+2;i_bin++){
        double err_tot = h_bkg_ord->GetBinError(i_bin); // this should be the stat unc
        double errUp(0);
        double errDown(0);
        double err(0);
        for(unsigned int i_syst=0;i_syst<systList.size();i_syst++){
            for(unsigned int j_syst=0;j_syst<systList.size();j_syst++){
                const double corr = fFitList[0]->fFitResults->GetCorrelationMatrix()->GetCorrelation( systList[i_syst],systList[j_syst] );
                errUp   += corr * h_syst_up_ord  [i_syst]->GetBinContent(i_bin) * h_syst_up_ord  [j_syst]->GetBinContent(i_bin);
                errDown += corr * h_syst_down_ord[i_syst]->GetBinContent(i_bin) * h_syst_down_ord[j_syst]->GetBinContent(i_bin);
            }
        }
        errUp   = std::sqrt(errUp);
        errDown = std::sqrt(errDown);
        err = 0.5*std::abs(errUp+errDown);
        err_tot = std::hypot(err_tot, err);
        h_bkg_ord->SetBinError(i_bin,err_tot);
    }

    TCanvas c("c","c",600,600);

    TPad pad0("pad0","pad0",0,0.28,1,1,0,0,0);
    pad0.SetTicks(1,1);
    pad0.SetTopMargin(0.05);
    pad0.SetBottomMargin(0);
    pad0.SetLeftMargin(0.14);
    pad0.SetRightMargin(0.05);
    pad0.SetFrameBorderMode(0);
    //
    TPad pad1("pad1","pad1",0,0,1,0.28,0,0,0);
    pad1.SetTicks(1,1);
    pad1.SetTopMargin(0.0);
    pad1.SetBottomMargin(0.37);
    pad1.SetLeftMargin(0.14);
    pad1.SetRightMargin(0.05);
    pad1.SetFrameBorderMode(0);

    pad1.Draw();
    pad0.Draw();
    pad0.cd();

    h_sig_ord->SetLineColor(kRed);
    h_sig_ord->SetFillColor(kRed);

    std::unique_ptr<TH1D> h_sig_ord_lim(static_cast<TH1D*>(h_sig_ord->Clone("h_sig_ord_lim")));
    h_sig_ord_lim->SetDirectory(nullptr);
    h_sig_ord_lim->Scale(muLimit/muFit);
    h_sig_ord_lim->SetFillColor(kOrange);
    h_sig_ord_lim->SetLineColor(kOrange);
    std::unique_ptr<TH1D> h_sig_ord_lim_diff(static_cast<TH1D*>(h_sig_ord_lim->Clone("h_sig_ord_lim_diff")));
    h_sig_ord_lim_diff->SetDirectory(nullptr);
    h_sig_ord_lim_diff->Add(h_sig_ord.get(),-1);

    THStack h_s{};
    h_s.Add(h_bkg_ord.get());
    h_s.Add(h_sig_ord.get());
    if (hasLimit) h_s.Add(h_sig_ord_lim_diff.get());
    h_data_ord->Draw("EX0");

    h_s.Draw("HISTsame");
    std::unique_ptr<TH1D> h_err(static_cast<TH1D*>(h_bkg_ord->Clone("h_err")));
    h_err->SetDirectory(nullptr);
    h_err->SetMarkerSize(0);
    h_err->SetFillColor(kBlack);
    h_err->SetFillStyle(3454);
    h_err->SetLineWidth(0);
    h_err->SetLineColor(kWhite);
    h_err->Draw("E2same");
    h_data_ord->Draw("EX0same");
    h_data_ord->SetMaximum(20*h_data_ord->GetMaximum());
    h_data_ord->SetMinimum(50);
    h_data_ord->SetLineWidth(2);
    h_data_ord->GetXaxis()->SetTitle("log_{10}(S/B)");
    h_data_ord->GetYaxis()->SetTitle("Events / 0.2");
    h_data_ord->GetYaxis()->SetTitleOffset(2.);
    h_data_ord->GetXaxis()->SetLabelSize(0);
    h_data_ord->GetXaxis()->SetTitleSize(0);

    if(includeBonly){
        h_bkgBonly_ord->SetLineColor(kBlack);
        h_bkgBonly_ord->SetLineStyle(kDashed);
        h_bkgBonly_ord->Draw("HISTsame");
    }

    if(TRExFitter::PREFITONPOSTFIT) {
      h_tot_bkg_prefit_ord->SetLineColor(kBlue);
      h_tot_bkg_prefit_ord->SetLineStyle(kDashed);
      h_tot_bkg_prefit_ord->Draw("HISTsame");
    }

    const std::string prefitLabel = fFitList.at(0)->fPreFitLabel;
    const std::string postfitLabel = fFitList.at(0)->fPostFitLabel;

    std::unique_ptr<TLegend> leg(nullptr);
    if(includeBonly) leg = std::make_unique<TLegend>(0.55,0.50,0.85,0.92);
    else             leg = std::make_unique<TLegend>(0.55,0.57,0.85,0.92);
    leg->SetFillStyle(0);
    leg->SetMargin(0.2);
    leg->SetBorderSize(0);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->AddEntry(h_data_ord.get(),"Data","lep");
    if (hasLimit) leg->AddEntry(h_sig_ord_lim.get(),Form(("%s ("+fPOIName.at(0)+"_{95%% excl.}=%.1f)").c_str(),fSignalTitle.c_str(),muLimit),"f");
    leg->AddEntry(h_sig_ord.get(),    Form(("%s ("+fPOIName.at(0)+"_{fit}=%.1f)"       ).c_str(),fSignalTitle.c_str(),muFit),"f");
    leg->AddEntry(h_bkg_ord.get(),"Background","f");
    leg->AddEntry(h_err.get(),"Bkgd. Unc.","f");
    if(includeBonly) leg->AddEntry(h_bkgBonly_ord.get(),("Bkgd. ("+fPOIName.at(0)+"=0)").c_str(),"l");
    if(TRExFitter::PREFITONPOSTFIT) leg->AddEntry(h_tot_bkg_prefit_comb.get(), (prefitLabel + " Bkgd.").c_str(),"l");
    leg->Draw();

    if (fFitList[0]->fPlotLabel != "none") TRExLabelNew(0.17,0.87, TRExFitter::EXPERIMENT_LABEL.c_str(), fFitList[0]->fPlotLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(0.17,0.80,kBlack,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()) );
    if(fLabel!="") myText(0.17,0.70,kBlack,Form("%s Combined",fLabel.c_str()) );
    else           myText(0.17,0.70,kBlack,"Combined");
    std::string channels = "";
    for(unsigned int i_fit=0;i_fit<fFitList.size();i_fit++){
        if (!fFitList.at(i_fit)->fUseInFit) continue;
        if(i_fit!=0){
            if(i_fit==fFitList.size()-1) channels += " and ";
            else                         channels += ", ";
        }
        channels += fFitList[i_fit]->fLabel;
    }
    myText(0.17,0.65,kBlack,channels.c_str());
    myText(0.17,0.57,kBlack,postfitLabel.c_str());

    pad0.RedrawAxis();
    pad0.SetLogy();

    pad1.cd();
    pad1.GetFrame()->SetY1(2);
    std::unique_ptr<TH1D> h_ratio(static_cast<TH1D*>(h_data_ord->Clone("h_ratio")));
    h_ratio->SetDirectory(nullptr);
    std::unique_ptr<TH1D> h_den(static_cast<TH1D*>(h_bkg_ord->Clone("h_den")));
    h_den->SetDirectory(nullptr);
    for(int i_bin=0;i_bin<h_den->GetNbinsX()+2;i_bin++){
        h_den->SetBinError(i_bin,0);
    }

    std::unique_ptr<TH1D> h_ratioBonly(nullptr);
    if(includeBonly){
        h_ratioBonly.reset(static_cast<TH1D*>(h_bkgBonly_ord->Clone("h_ratioBonly")));
        h_ratioBonly->SetDirectory(nullptr);
        h_ratioBonly->Divide(h_den.get());
        h_ratioBonly->SetLineStyle(kDashed);
        h_ratioBonly->SetLineColor(kBlack);
    }

    std::unique_ptr<TH1D> h_stackSig(static_cast<TH1D*>(h_sig_ord ->Clone("h_sig_ratio")));
    h_stackSig->SetDirectory(nullptr);
    h_stackSig->Add(h_bkg_ord.get());
    h_stackSig->Divide(h_den.get());
    h_stackSig->SetFillColor(0);
    h_stackSig->SetFillStyle(0);
    h_stackSig->SetLineColor(kRed);

    std::unique_ptr<TH1D> h_stackSigLim(static_cast<TH1D*>(h_sig_ord_lim ->Clone("h_sig_lim_ratio")));
    h_stackSigLim->SetDirectory(nullptr);
    h_stackSigLim->Add(h_bkg_ord.get());
    h_stackSigLim->Divide(h_den.get());
    h_stackSigLim->SetFillColor(0);
    h_stackSigLim->SetFillStyle(0);
    h_stackSigLim->SetLineStyle(kDashed);
    h_stackSigLim->SetLineColor(kOrange+1);

    std::unique_ptr<TH1D> h_ratio2(static_cast<TH1D*>(h_err->Clone("h_ratio2")));
    h_ratio2->SetDirectory(nullptr);
    h_ratio2->SetMarkerSize(0);
    h_ratio->SetTitle("Data/MC");
    h_ratio->GetYaxis()->SetTitle("Data / Bkgd.");
    h_ratio->GetYaxis()->SetTitleSize(20);
    h_ratio->GetYaxis()->SetTitleOffset(2.);
    h_ratio->GetYaxis()->SetLabelSize(20); // 0.04
    h_ratio ->Divide(h_den.get());
    h_ratio2->Divide(h_den.get());
    h_ratio->SetMarkerSize(1.2);
    h_ratio->SetLineWidth(2);
    gStyle->SetEndErrorSize(0.); // 4.
    h_ratio->GetYaxis()->CenterTitle();
    h_ratio->GetYaxis()->SetNdivisions(406);
    h_ratio->SetMinimum(0.6);
    h_ratio->SetMaximum(1.75);
    h_ratio->GetXaxis()->SetTitle(h_data_ord->GetXaxis()->GetTitle());
    h_ratio->GetXaxis()->SetTitleSize(20);
    h_ratio->GetXaxis()->SetTitleOffset(4.);
    h_ratio->GetXaxis()->SetLabelSize(20);
    TLine hline(h_ratio->GetXaxis()->GetXmin(),1,h_ratio->GetXaxis()->GetXmax(),1);
    h_ratio->Draw("E1");
    h_ratio2->Draw("same E2");
    hline.Draw();
    h_stackSig->Draw("same HIST");
    if (hasLimit) h_stackSigLim->Draw("same HIST");
    h_ratio->Draw("same E1");

    if(includeBonly) h_ratioBonly->Draw("same HIST");

    TLegend leg2(0.17,0.64,0.75,0.98);
    leg2.SetFillStyle(0);
    leg2.SetMargin(0.1);
    leg2.SetBorderSize(0);
    leg2.SetTextSize(gStyle->GetTextSize());
    if (hasLimit) leg2.AddEntry(h_stackSigLim.get(),Form(("%s ("+fPOIName.at(0)+"_{95%% excl.}=%.1f) + Bkgd.").c_str(),fSignalTitle.c_str(),muLimit),"l");
    leg2.AddEntry(h_stackSig.get(),   Form(("%s ("+fPOIName.at(0)+"_{fit}=%.1f) + Bkgd."       ).c_str(),fSignalTitle.c_str(),muFit)  ,"l");
    leg2.Draw();

    std::unique_ptr<TLegend> leg3;
    if(includeBonly){
        leg3 = std::make_unique<TLegend>(0.17,0.4,0.75,0.5);
        leg3->SetFillStyle(0);
        leg3->SetMargin(0.1);
        leg3->SetBorderSize(0);
        leg3->SetTextSize(gStyle->GetTextSize());
        leg3->AddEntry(h_ratioBonly.get(),"Bkgd. (from Bkgd-only fit)","l");
        leg3->Draw();
    }


    pad0.RedrawAxis();
    pad1.RedrawAxis();

    Common::SaveCanvasAs(c, fOutDir+"/SoverB_postFit");
    for (auto& ifile : file) {
        ifile->Close();
    }
    for (auto& ifile : fileBonly) {
        ifile->Close();
    }

}

//____________________________________________________________________________________
//
std::unique_ptr<TH1D> MultiFit::Combine(const std::vector<std::unique_ptr<TH1D> >& h) const{
    if (h.empty()) return nullptr;
    int Nbins = 0;
    for(std::size_t i_hist=0;i_hist < h.size();i_hist++){
        if(h[i_hist]==nullptr) LOG(WARNING) << "Empty histogram " << i_hist << "\n";
        else Nbins += h[i_hist]->GetNbinsX();
    }
    auto h_new = std::make_unique<TH1D>(Form("%s_comb",h[0]->GetName()),Form("%s_comb",h[0]->GetTitle()),Nbins,0,Nbins);
    int bin = 0;
    for(const auto& ihist : h) {
        if (!ihist) continue;
        for(int i_bin=1;i_bin <= ihist->GetNbinsX();i_bin++){
            bin++;
            h_new->SetBinContent(bin,ihist->GetBinContent(i_bin));
            h_new->SetBinError(bin,ihist->GetBinError(i_bin));
        }
    }
    h_new->SetDirectory(nullptr);
    return h_new;
}

//____________________________________________________________________________________
// merge bins in bins of SoverSqrtB
std::unique_ptr<TH1D> MultiFit::Rebin(TH1D* h, const vector<double>& vec, bool isData) const{
    std::unique_ptr<TH1D> h_new = std::make_unique<TH1D>(Form("%s_rebin",h->GetName()),Form("%s_rebin",h->GetTitle()),17,-3.8,-0.5);
    h_new->SetDirectory(nullptr);
    h_new->Sumw2();
    // new way
    for(int j_bin=1;j_bin<=h->GetNbinsX();j_bin++){
        double value=std::log10(vec[j_bin-1]);
        if ( value<h_new->GetXaxis()->GetXmin() ) value=0.9999*h_new->GetXaxis()->GetXmin();
        if ( value>h_new->GetXaxis()->GetXmax() ) {
            double tmpvalue=1.0001*h_new->GetXaxis()->GetXmax();
            LOG(DEBUG) << "Turning: " << value << " in: " << tmpvalue << "\n";
            value=tmpvalue;
        }
        int i_bin=h_new->FindBin(value);
        h_new->SetBinContent(i_bin,h_new->GetBinContent(i_bin)+h->GetBinContent(j_bin));
        if (!isData) {
            h_new->SetBinError(i_bin,std::hypot(h_new->GetBinError(i_bin), h->GetBinError(j_bin)));
        }
    }
    if (isData) {
        for(int j_bin=1;j_bin<=h_new->GetNbinsX();j_bin++){
            h_new->SetBinError(j_bin, std::sqrt(h_new->GetBinContent(j_bin) ) );
        }
    }
    h_new->SetMinimum(1);
    return h_new;
}

//____________________________________________________________________________________
// combine individual results from grouped impact evaluation into one table
void MultiFit::BuildGroupedImpactTable() const {
    LOG(INFO) << "Merging grouped impact evaluations\n";
    for (const auto& ipoi : fPOIs) {
        const std::string targetName = fOutDir+"/Fits/GroupedImpact_"+ipoi+fSaveSuf+".txt";

        if(std::ifstream(targetName).good()){
            LOG(WARNING) << "File " << targetName << " already exists, will not overwrite\n";
        } else {
            const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fOutDir+"/Fits/","GroupedImpact_", ipoi+fSaveSuf+".txt");
            Common::MergeTxTFiles(inPaths, targetName);
        }
    }

}

//__________________________________________________________________________________
//
std::vector<std::shared_ptr<NormFactor> > MultiFit::GetFitNormFactors() const {
    std::vector<std::shared_ptr<NormFactor> > result;

    std::vector<std::string> names;
    for (const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(const auto& nf : ifit->fNormFactors) {
            // only add the unique NFs
            if (std::find(names.begin(), names.end(), nf->fName) != names.end()) continue;
            names.emplace_back(nf->fName);

            result.emplace_back(nf);
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<std::shared_ptr<ShapeFactor> > MultiFit::GetFitShapeFactors() const {
    std::vector<std::shared_ptr<ShapeFactor> > result;

    std::vector<std::string> names;
    for (const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(const auto& sf : ifit->fShapeFactors) {
            // only add the unique NFs
            if (std::find(names.begin(), names.end(), sf->fName) != names.end()) continue;
            names.emplace_back(sf->fName);

            result.emplace_back(sf);
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::map<std::string, double> MultiFit::GetFixedNPs() const {

    std::map<std::string, double> result;
    for (const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for (const auto& ifixed : ifit->fFitFixedNPs) {
            auto it = result.find(ifixed.first);
            if (it != result.end()) continue;

            result[ifixed.first] = ifixed.second;
        }
    }
    return result;
}

//__________________________________________________________________________________
//
std::vector<std::shared_ptr<Systematic> > MultiFit::GetFitSystematics() const {

    std::vector<std::shared_ptr<Systematic> > result;

    std::vector< std::string > names;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(const auto& isyst : ifit->fSystematics) {
            const std::string systName = isyst->fNuisanceParameter;
            if(Common::FindInStringVector(names,systName) < 0) {
                names.push_back(systName);
                result.emplace_back(isyst);
            }
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<Region* > MultiFit::GetFitRegions() const {
    std::vector<Region*> result;

    std::vector< std::string > names;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for (auto& ireg : ifit->fRegions) {
            const std::string regName = ireg->fName;
            if (!fOnlyRegions.empty() && Common::FindInStringVector(fOnlyRegions, regName) < 0) continue;
            if(Common::FindInStringVector(names,regName) < 0) {
                names.emplace_back(regName);
                result.push_back(ireg.get());
            }
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<std::unique_ptr<Unfolding> > MultiFit::GetUnfoldings() const {
    std::vector<std::unique_ptr<Unfolding> > result;

    std::vector<std::string> names;
    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for (auto& iunf : ifit->fUnfolding) {
            const std::string unfName = iunf->fName;
            if(Common::FindInStringVector(names,unfName) < 0) {
                names.emplace_back(unfName);
                result.emplace_back(std::make_unique<Unfolding>(*iunf));
            }
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void MultiFit::RunLimitToys(RooAbsData* data, RooWorkspace* ws) const {
    FitUtils::SetBinnedLikelihoodOptimisation(ws);
    gSystem->mkdir((fOutDir + "/Limits").c_str());

    LimitToys toys{};
    toys.SetPlot(fLimitPlot);
    toys.SetFile(fLimitFile);
    toys.SetNToys(fLimitToysStepsSplusB, fLimitToysStepsB);
    toys.SetLimit(fLimitsConfidence);
    toys.SetScan(fLimitToysScanSteps, fLimitToysScanMin, fLimitToysScanMax);
    toys.SetOutputPath(fOutDir + "/Limits");
    toys.SetSeed(fLimitToysSeed);
    toys.SetLimitSuffix("_"+std::to_string(fLimitToysSeed));

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read ModelConfig\n";
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<RooStats::ModelConfig> mcBonly(static_cast<RooStats::ModelConfig*>(mc->Clone("BonlyModel")));
    RooRealVar* poi = static_cast<RooRealVar*>(mcBonly->GetParametersOfInterest()->first());
    poi->setVal(0);

    toys.RunToys(data, mc, mcBonly.get());
}

//__________________________________________________________________________________
//
void MultiFit::RunSignificanceToys(RooAbsData* data, RooWorkspace* ws) const {
    FitUtils::SetBinnedLikelihoodOptimisation(ws);
    gSystem->mkdir((fOutDir + "/Significance").c_str());

    SignificanceToys toys{};
    toys.SetPlot(fSignificancePlot);
    toys.SetNtoys(fSignificanceToysStepsSplusB, fSignificanceToysStepsB);
    toys.SetOutputPath(fOutDir + "/Significance");
    toys.SetSignificanceSuffix("_"+std::to_string(fSignificanceToysSeed));
    toys.SetToysSeed(fSignificanceToysSeed);

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read ModelConfig\n";
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<RooStats::ModelConfig> mcBonly(static_cast<RooStats::ModelConfig*>(mc->Clone("BonlyModel")));
    RooRealVar* poi = static_cast<RooRealVar*>(mcBonly->GetParametersOfInterest()->first());
    poi->setVal(0);

    toys.RunToys(data, mc, mcBonly.get());
}

//__________________________________________________________________________________
//
void MultiFit::AddPOI(const std::string& name) {
    if(Common::FindInStringVector(fPOIs,name)<0){
        fPOIs.emplace_back(name);
    }
    else{
        LOG(WARNING) << "POI " << name << " already set. Skipping.\n";
    }
}

//__________________________________________________________________________________
//
std::vector<std::pair<std::string,double> > MultiFit::GetSmallNPs() const {
    std::vector<std::pair<std::string,double> > result{};
    for(const auto& ifit : fFitList) {
        const std::string path = ifit->fName+"/NPorder.txt";
        std::ifstream file;
        file.open(path);
        if (!file.is_open() || !file.good()) {
            LOG(WARNING) << "Cannot open file at " << path << "\n";
            continue;
        }
        std::string syst;
        double value;
        while(file >> syst >> value) {
            auto itr = std::find_if(result.begin(), result.end(), [&syst](const std::pair<std::string,double>& element){return syst == element.first;});
            if (itr == result.end()) {
                result.emplace_back(syst, value);
            } else {
                if (itr->second < value) {
                    itr->second = value;
                }
            }
        }
    }

    std::sort(result.begin(), result.end(), [](const std::pair<std::string, double>& lhs, const std::pair<std::string, double>& rhs){return lhs.second < rhs.second;});
    result.resize(static_cast<std::size_t>(fNPCutOff * result.size()));
    return result;
}

//__________________________________________________________________________________
//
void MultiFit::PlotUnfoldedData() {

    std::vector<std::string> unfoldings;
    std::vector<const Unfolding*> combinedUnfoldings;
    std::vector<std::string> truthFileNames;

    for (auto& ifit : fFitList) {
        if (fUnfoldingShowStat) {
            ifit->fUnfoldingShowStat = fUnfoldingShowStat;
        }
        bool hasUnfolding(false);
        for (auto& iunfolding : ifit->fUnfolding) {
            const std::string name = iunfolding->fName;
            if (std::find(unfoldings.begin(), unfoldings.end(), name) != unfoldings.end()) continue;
            LOG(INFO) << "Plotting combined unfolding for " << name << "\n";
            unfoldings.emplace_back(name);
            combinedUnfoldings.emplace_back(iunfolding.get());
            hasUnfolding = true;
            ifit->PlotUnfoldedData(fOutDir+"/Fits/"+fName+fSaveSuf,
                                   fOutDir,
                                   fOutDir+"/ws_combined.root",
                                   "combWS");
        }
        if (hasUnfolding) {
            truthFileNames.emplace_back(ifit->fName);
        }
    }

    if (combinedUnfoldings.size() > 1) {
        fFitList.at(0)->PlotMultipleUnfoldingCovariance(combinedUnfoldings,
                                                        truthFileNames,
                                                        fOutDir+"/Fits/"+fName+fSaveSuf,
                                                        fOutDir+"/ws_combined.root",
                                                        fOutDir,
                                                        "combWS",
                                                        this->GetFitNormFactors());
    }

}

//__________________________________________________________________________________
//
void MultiFit::DrawGammas() {
    LOG(INFO) << "Drawing the gamma plot\n";
    if (fFitList.empty()) {
        LOG(ERROR) << "Empty fit list!\n";
        return;
    }

    auto sfs = this->GetFitShapeFactors();

    fFitList.at(0)->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".root");
    fFitList.at(0)->fFitResults->DrawGammaShapePulls(fOutDir+"/Gammas", {}, sfs, false);
    fFitList.at(0)->fFitResults->DrawGammaShapePulls(fOutDir+"/ShapeFactors", {}, sfs, true);
}

//__________________________________________________________________________________
//
void MultiFit::AdjustBlindedParams() {
    for (const auto& ifit : fFitList) {
        const std::vector<std::string>& blinded = ifit->fBlindedParameters;
        for (const auto& iparam : blinded) {
            auto itr = std::find(fBlindedParams.begin(), fBlindedParams.end(), iparam);
            if (itr == fBlindedParams.end()) {
                fBlindedParams.emplace_back(iparam);
            }
        }
    }
}

//__________________________________________________________________________________
//
void MultiFit::ProcessShapeFactorReparametrisation() {
    for (auto& ifit : fFitList) {
        ifit->ReadRegionBinning();
        ifit->ProcessTemplateMorphingShapeFactors();
        ifit->ProcessEFTShapeFactors();
        ifit->ProcessShapeFactorReparametrization();
        for (const auto& isf : ifit->fShapeFactors) {
            if (!isf->fExpression.empty()) {
                fHasShapeFactorReparametrisation = true;
            }
        }
    }
}

//__________________________________________________________________________________
//
void MultiFit::CheckSameRegionNames() const {
    if (fFitList.empty()) return;

    std::vector<std::string> regions;
    std::size_t startfit(0);
    for (; startfit < fFitList.size(); ++startfit) {
        if (!fFitList.at(startfit)->fUseInFit) continue;
        for (const auto& ireg : fFitList.at(startfit)->fRegions) {
            if (ireg->fRegionType == Region::RegionType::VALIDATION) continue;
            regions.emplace_back(ireg->fName);
        }
        break;
    }

    for (std::size_t ifit = startfit+1; ifit < fFitList.size(); ++ifit) {
        if (!fFitList.at(ifit)->fUseInFit) continue;
        for (const auto& ireg : fFitList.at(ifit)->fRegions) {
            if (ireg->fRegionType == Region::RegionType::VALIDATION) continue;
            auto itr = std::find(regions.begin(), regions.end(), ireg->fName);
            if (itr == regions.end()) {
                regions.emplace_back(ireg->fName);
            } else {
                LOG(WARNING) << "Region: \"" << ireg->fName << "\" appears multiple times in the config files\n";
                LOG(WARNING) << "Regions with the same name will NOT be combined in the fit\n";
            }
        }
    }
}

//__________________________________________________________________________________
//
void MultiFit::PlotRankingFromCovMatrix(const FittingTool& fitTool) const {
    if (fFitType == TRExFit::FitType::BONLY) return;
    RankingManager manager{};
    manager.SetPlotLabel(fFitList.at(0)->fPlotLabel);
    manager.SetLumiLabel(fFitList.at(0)->fLumiLabel);
    manager.SetCmeLabel(fFitList.at(0)->fCmeLabel);
    manager.SetUseHEPDataFormat(fHEPDataFormat);
    manager.SetName(fOutDir);
    manager.SetMaxNPPlot(fFitList.at(0)->fRankingMaxNP);
    manager.SetRankingCanvasSize(fFitList.at(0)->fNPRankingCanvasSize);
    manager.SetShiftGlobalObservables(true);

    const auto nfs = this->GetFitNormFactors();
    const auto sfs = this->GetFitShapeFactors();

    std::vector<std::string> nfNames;
    std::vector<std::string> sfNames;

    std::vector<std::string> exclude = fPOIs;
    for (const auto& inf : nfs) {
        exclude.emplace_back(inf->fName);
        nfNames.emplace_back(inf->fName);
    }
    for (const auto& isf : sfs) {
        exclude.emplace_back(isf->fName);
        sfNames.emplace_back(isf->fName);
    }

    for (std::size_t ipoi = 0; ipoi < fPOIs.size(); ++ipoi) {
        manager.SetRankingPOIName(fRankingPOIName.at(ipoi));
        manager.SetUpperAxisNdivisions(fRankingUpperAxisNdivision.at(ipoi));
        manager.SetRankingPOIAxisScale(fRankingPOIAxisScale.at(ipoi));
        manager.SetIsCovarianceBreakdown(true);
        if (fFitList.at(0)->fRankingPlot == "MERGE" || fFitList.at(0)->fRankingPlot == "ALL") {
            const auto container = fitTool.GetRankingContainer(fPOIs.at(ipoi), exclude, true, true);
            manager.SetRankingContainer(container);
            manager.SetSuffix(fSaveSuf+"_"+fPOIs.at(ipoi)+"_Breakdown");
            manager.PlotRanking(this->GetFitRegions(), nfNames, sfNames, true, true);
        }
        if (fFitList.at(0)->fRankingPlot == "SYSTS" || fFitList.at(0)->fRankingPlot == "ALL") {
            const auto container = fitTool.GetRankingContainer(fPOIs.at(ipoi), exclude, true, false);
            manager.SetRankingContainer(container);
            manager.SetSuffix(fSaveSuf+"_"+fPOIs.at(ipoi)+"_Breakdown_syst");
            manager.PlotRanking(this->GetFitRegions(), nfNames, sfNames, true, false);
        }
        if (fFitList.at(0)->fRankingPlot == "GAMMAS" || fFitList.at(0)->fRankingPlot == "ALL") {
            const auto container = fitTool.GetRankingContainer(fPOIs.at(ipoi), exclude, false, true);
            manager.SetRankingContainer(container);
            manager.SetSuffix(fSaveSuf+"_"+fPOIs.at(ipoi)+"_Breakdown_gammas");
            manager.PlotRanking(this->GetFitRegions(), nfNames, sfNames, false, true);
        }
    }
}

//__________________________________________________________________________________
//
std::map<std::string, std::string> MultiFit::GetSubcategoryMap() const {
    std::map<std::string, std::string> mergedMap;
    for(auto& fit : fFitList){
        if (!fit->fUseInFit) continue;
        fit->ProduceSystSubCategoryMap();
        for(const auto& m : fit->fSubCategoryImpactMap) {
            if(mergedMap[m.first] == "") {
                mergedMap[m.first] = m.second;
            } else if(mergedMap[m.first] != m.second){
                LOG(WARNING) << "Systematics for NP: " << m.first << " assigned to different SubCategory in the different included configs.\n";
                LOG(WARNING) << "Combined category: " << mergedMap[m.first] << ", current category for config " << fit->fName << ": " << m.second << "\n";
                LOG(WARNING) << "Keeping the current combined category (combined map)\n";
            }
        }
    }

    return mergedMap;
}

//__________________________________________________________________________________
//
void MultiFit::TranslateWSToHS3() const {

    if (fFitList.empty()) {
        LOG(ERROR) << "Fit list is empty, this should not happen\n";
        exit(EXIT_FAILURE);
    }

    const auto* fit = fFitList.at(0).get();

    if (!fit) {
        LOG(ERROR) << "Fit pointer is null, this should not happen\n";
        exit(EXIT_FAILURE);
    }

    const std::string wsFilePath = fOutDir+"/ws_combined"+fSaveSuf+".root";
    const std::string wsPath = "combWS";
    const std::string outPath = fOutDir+"/Workspace_HS3.json";

    fit->TranslateWSToHS3(wsFilePath, wsPath, outPath);
}

//__________________________________________________________________________________
//
void MultiFit::PlotUnfoldedErrors(const FittingTool& fitTool, const std::vector<std::string>& categories) const {

    if (!fUnfoldingShowUncertaintyBreakdown) return;

    auto isNorm =[](const std::vector<std::unique_ptr<Unfolding> >& unfolding) {
        for (const auto& iunf : unfolding) {
            if (iunf->fUnfoldNormXSec) return true;
        }
        return false;
    };

    gSystem->mkdir((fOutDir+"/UnfoldingPlots").c_str());
    UncertaintyPlotter plotter{};
    plotter.SetOutputFolder(fOutDir+"/UnfoldingPlots");
    plotter.SetUseTotalUncertainty(fUnfoldingUncertaintyBreakdownTotal);
    plotter.SetCMELabel(fCmeLabel);
    plotter.SetPlotLabel(fFitList[0]->fPlotLabel);
    plotter.SetLumiLabel(fLumiLabel);
    auto unfoldings = this->GetUnfoldings();
    plotter.SetIsNormalized(isNorm(unfoldings));
    plotter.SetCategories(categories);
    plotter.SetSubCategoryMap(this->GetSubcategoryMap());
    for (const auto& iunfolding : unfoldings) {
        std::vector<std::pair<std::string, int> > pois;
        std::vector<std::string> otherPOIs;
        for (int ibin = 0; ibin < iunfolding->fNumberUnfoldingTruthBins; ++ibin) {
            const std::string name = iunfolding->fName + "_Bin_" + Common::IntToFixLenStr(ibin+1) + "_mu";
            pois.emplace_back(std::make_pair(name, ibin));
        }

        for (const auto& ipoi : fPOIs) {
            if (Common::IsUnfoldingPOI(ipoi, unfoldings)) continue;
            auto itr = std::find_if(pois.begin(), pois.end(), [&ipoi](const auto& element){return element.first == ipoi;});
            if (itr == pois.end()) {
                otherPOIs.emplace_back(ipoi);
            }
        }

        plotter.SetPOIs(pois);
        plotter.SetOtherPOIs(otherPOIs);
        plotter.SetXaxisLabel(iunfolding->fUnfoldingTitleX);

        std::unique_ptr<TFile> input(TFile::Open((fFitList[0]->fName + "/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
        if (!input) {
            LOG(ERROR) << "Cannot read file from " << fFitList[0]->fName << "/UnfoldingHistograms/FoldedHistograms.root\n";
            exit(EXIT_FAILURE);
        }

        std::unique_ptr<TH1> truth(dynamic_cast<TH1*>(input->Get((iunfolding->fName+"_truth_distribution").c_str())));
        if (!truth) {
            LOG(ERROR) << "Cannot read the truth distribution\n";
            exit(EXIT_FAILURE);
        }
        truth->SetDirectory(nullptr);
        plotter.PlotUncertainties(truth.get(),
                                  *iunfolding,
                                  fitTool,
                                  fOutDir+"/ws_combined.root",
                                  "combWS",
                                  fOutDir+"/Fits/"+fName+fSaveSuf+".root",
                                  this->GetFitNormFactors());
        input->Close();
    }
}

//__________________________________________________________________________________
//
void MultiFit::RunToys() const {

    gSystem->mkdir((fName+"/Toys").c_str());

    LOG(INFO) << "Generating and fitting toys\n";
    LOG(DEBUG) << "Reading combined WS\n";
    auto ws = this->CombineWS();

    fFitList.at(0)->ReadFitResults(fOutDir+"/Fits/"+fName+fSaveSuf+".root");

    auto nfs = this->GetFitNormFactors();

    bool isRegularizedUnfolding(false);
    if (fFitType == TRExFit::FitType::UNFOLDING) {
        for (const auto& inf : nfs) {
            if (inf->fTau > 0) {
                isRegularizedUnfolding = true;
                break;
            }
        }
    }

    std::map<std::string, double> poiAsimov;
    if (fPOIAsimov.size() != fPOIName.size()) {
        LOG(ERROR) << "Size of POIAsimov and POINames do not match\n";
        return;
    }
    for (std::size_t i = 0; i < fPOIAsimov.size(); ++i) {
        poiAsimov.insert({fPOIName.at(i), fPOIAsimov.at(i)});
    }

    FitToys fitToys(fName, fSaveSuf, "");
    fitToys.SetMinos(this->GetMinos());
    fitToys.SetIsRegularizedUnfolding(isRegularizedUnfolding);
    fitToys.SetNormFactors(nfs);
    fitToys.SetNToys(fFitToys);
    fitToys.SetRegularizationType(this->GetRegularizationType());
    fitToys.SetUseAutoDiff(fUseAutoDiff);
    fitToys.SetFitType(fFitType);
    fitToys.SetHistoNbins(fToysHistoNbins);
    fitToys.SetPOIs(fPOIs);
    fitToys.SetFitIsBlind(fDataName == "asimovData");
    fitToys.SetPOIAsimov(poiAsimov);
    fitToys.SetSeed(fToysSeed);
    fitToys.SetStatOnlyFluctuation(fToysSetStatOnlyFluctuation);
    fitToys.SetRandomPOIStartingValue(fToysRandomPOIStartingValues);
    fitToys.SetPlotLabel(fFitList.at(0)->fPlotLabel);
    fitToys.SetFitFixedNPs(this->GetFixedNPs());
    fitToys.SetRunLHScanForAll(fToysLHScanForAll);
    fitToys.SetLHScanCondition(fToysLHScanCondition);

    fitToys.RunToys(ws.get(), fFitList.at(0)->fFitResults.get());
}

//__________________________________________________________________________________
//
std::vector<std::string> MultiFit::GetMinos() const {

    std::vector<std::string> result;

    for(const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;
        for(const auto& iminos : ifit->fVarNameMinos) {
            if(Common::FindInStringVector(result,iminos) < 0){
                result.push_back(iminos);
            }
        }
    }

    return result;
}

//__________________________________________________________________________________
//
int MultiFit::GetRegularizationType() const {
    int result(-999);

    for (const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;

        if (result > -10) {
            if (result != ifit->fRegularizationType) {
                LOG(ERROR) << "Inconsistent regularization types for the individual configs. This is not supported for multifit\n";
                exit(EXIT_FAILURE);
            }
        } else {
            result = ifit->fRegularizationType;
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void MultiFit::CheckDatasetConsistency() const {

    int nFitConfigs(0);
    int nDataInConfigs(0);

    for (const auto& ifit : fFitList) {
        if (!ifit->fUseInFit) continue;

        if (ifit->fHasDataSetInConfig) ++nDataInConfigs;
        ++nFitConfigs;
    }

    // no data is fine
    if (nDataInConfigs == 0) return;

    if (nFitConfigs != nDataInConfigs) {
        LOG(ERROR) << "Data samples are inconsistent in the fit configs\n";
        LOG(ERROR) << "You need to either provide all WSs with the Data samples or none, otherwise HistFactory would crash\n";
        LOG(ERROR) << "Please, fix this and rerun the \"w\" step in the single configs\n";
        exit(EXIT_FAILURE);
    }
}
