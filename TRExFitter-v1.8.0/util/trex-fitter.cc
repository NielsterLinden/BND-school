// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReader.h"
#include "TRExFitter/ConfigReaderMulti.h"
#include "TRExFitter/HistoReader.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/MultiFit.h"
#include "TRExFitter/NtupleReader.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/TRExFit.h"
#include "TRExFitter/TRExPlot.h"
#include "TRExFitter/UnfoldingSample.h"
#include "TRExFitter/UnfoldingSystematic.h"
#include "TRExFitter/YamlConverter.h"

// RooStatsIncludes
#include "RooStats/RooStatsUtils.h"

// Style stuff
#include "StyleUtils/TRExStyle.h"

// ROOT includes
#include "TH1.h"
#include "TSystem.h"
#include "TROOT.h"

// c++ includes
#include <memory>
#include <string>
#include <filesystem>
#include <fstream>

#define _STRINGIZE(x) #x
#define STRINGIZE(x) _STRINGIZE(x)

// -------------------------------------------------------
// -------------------------------------------------------

// trick to suppress the RooFit banner
// int doBanner(){ return 0; }

void FitExample(const std::string& opt="h", const std::string& configFile="config/myFit.config", const std::string& options=""){

    // pre-read the config just to extract the debug level
    const std::string debugStr = ReadValueFromConfig(configFile,"DebugLevel");
    if(debugStr=="") LOG(WARNING) << "Not able to pre-read the DebugLevel => keep the default (1).\n";
    else TRExFitter::SetDebugLevel(std::atoi(debugStr.c_str()));

    const std::string loggingFormat = ReadValueFromConfig(configFile,"LoggingFormat");
    if (loggingFormat != "") TRExFitter::SetLoggingFormat(loggingFormat);

    if (TRExFitter::DEBUGLEVEL < 2){
        gErrorIgnoreLevel = kError;
        RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }
    // pre-read the logo option
    std::string logoStr = ReadValueFromConfig(configFile,"Logo");
    std::transform(logoStr.begin(), logoStr.end(), logoStr.begin(), ::toupper);
    if(logoStr=="TRUE"){
        TString homeArea(gSystem->Getenv("$TREXFITTER_HOME") ? "$TREXFITTER_HOME" : "$DATAPATH");
#ifdef TREXFITTER_HOME
        homeArea = std::string(STRINGIZE(TREXFITTER_HOME));
#endif
        TString file = "logo.txt";
        gSystem->FindFile(gSystem->ExpandPathName(homeArea+":./"),file);
        if (file == "") {
            LOG(ERROR) << "logo.txt not found. Are you running from the main working directory?\n";
            exit(EXIT_FAILURE);
        }
        std::ifstream logoFile(file);
        std::string str;
        std::string logo = "";
        while(getline(logoFile,str)){
            if(!logoFile.good()) break;
            logo+=str;
            logo+="\n";
        }
        std::cout << logo << std::endl;
    }

    if(TRExFitter::DEBUGLEVEL<=0)      gErrorIgnoreLevel = kError;
    else if(TRExFitter::DEBUGLEVEL<=1) gErrorIgnoreLevel = kWarning;

    // now can set plot style
    SetTRExStyle();

    RooStats::UseNLLOffset(true);

    // interpret opt
    const bool readHistograms     = opt.find("h") != std::string::npos;
    const bool readNtuples        = opt.find("n") != std::string::npos;
    const bool rebinAndSmooth     = opt.find("b") != std::string::npos;
    const bool doEFTInputs        = opt.find("e") != std::string::npos;
    const bool createWorkspace    = opt.find("w") != std::string::npos;
    const bool doFit              = opt.find("f") != std::string::npos;
    const bool doRanking          = opt.find("r") != std::string::npos;
    const bool doLimit            = opt.find("l") != std::string::npos;
    const bool doSignificance     = opt.find("s") != std::string::npos;
    bool drawPreFit         = opt.find("d") != std::string::npos;
    const bool drawPostFit        = opt.find("p") != std::string::npos;
    const bool drawSeparation     = opt.find("a") != std::string::npos;
    const bool groupedImpact      = opt.find("i") != std::string::npos;
    const bool doLHscan           = opt.find("x") != std::string::npos;
    const bool prepareUnfolding   = opt.find("u") != std::string::npos;
    const bool plotFitResults     = opt.find("q") != std::string::npos;
    const bool prepareTemplate    = opt.find("t") != std::string::npos;
    const bool translateHS3       = opt.find("j") != std::string::npos;

    if(!readNtuples && !rebinAndSmooth){
        TH1::AddDirectory(kFALSE); // FIXME: it would be nice to have a solution which works always
    }

    // multi-fit
    const bool isMultiFit      = opt.find("m")!=std::string::npos;
    if(isMultiFit){
        std::unique_ptr<MultiFit> myMultiFit = std::make_unique<MultiFit>("MyMultiFit");
        ConfigReaderMulti confReaderMulti(myMultiFit.get());
        int sc = confReaderMulti.ReadFullConfig(configFile,opt,options) ;

        if (sc != 0){
            LOG(ERROR) << "Failed to read the config file for multifit.\n";
            LOG(ERROR) << "Found " << sc << " errors in the config.\n";
            LOG(ERROR) << "Please check all reported errors.\n";
            exit(EXIT_FAILURE);
        }
        myMultiFit->AdjustBlindedParams();
        myMultiFit->ProcessShapeFactorReparametrisation();
        if (myMultiFit->fHEPDataFormat) {
            for (const auto& ifit: myMultiFit->fFitList) {
                ifit->fHEPDataFormat = true;
            }
        }

        if (plotFitResults) {
            myMultiFit->PlotUnfoldedData();
        }

        if(myMultiFit->fCombine){
            if(createWorkspace){
                myMultiFit->CheckDatasetConsistency();
                LOG(INFO) << "Combining workspaces..." << std::endl;
                myMultiFit->SaveCombinedWS();
            }
            if(doFit){
                LOG(INFO) << "Fitting combined workspace..." << std::endl;
                myMultiFit->FitCombinedWS( myMultiFit->fFitType, myMultiFit->fDataName, false );
                if (myMultiFit->fPlotUnfoldedData) {
                    myMultiFit->PlotUnfoldedData();
                }
                if (myMultiFit->fPlotGammas) {
                    myMultiFit->DrawGammas();
                }
            }
            if(doLHscan){
                LOG(INFO) << "Doing likelihood scan for the combined workspace..." << std::endl;
                myMultiFit->FitCombinedWS( myMultiFit->fFitType, myMultiFit->fDataName, true );
            }
            if(doLimit){
                LOG(INFO) << "Getting combined limit..." << std::endl;
                myMultiFit->GetCombinedLimit( myMultiFit->fDataName );
            }
            if(doSignificance){
                LOG(INFO) << "Getting combined significance..." << std::endl;
                myMultiFit->GetCombinedSignificance( myMultiFit->fDataName );
            }
            if(doRanking){
                LOG(INFO) << "Getting combined ranking..." << std::endl;
                if (myMultiFit->fRankingOnly != "plot") myMultiFit->ProduceNPRanking(myMultiFit->fRankingOnly);
                if (myMultiFit->fRankingOnly=="all" || myMultiFit->fRankingOnly=="plot")  myMultiFit->PlotNPRankingManager();
            }
            if(groupedImpact){
                LOG(INFO) << "Getting combined grouped systematic impact..." << std::endl;
                myMultiFit->fDoGroupedSystImpactTable = true;
                if(myMultiFit->fGroupedImpactCategory!="combine") myMultiFit->FitCombinedWS( myMultiFit->fFitType, myMultiFit->fDataName, false );
                else                                              myMultiFit->BuildGroupedImpactTable();
            }
        }
        //
        if(myMultiFit->fCompare){
            LOG(INFO) << "Comparing fits..." << std::endl;
            if(myMultiFit->fComparePulls){
                for(const auto& icat : myMultiFit->fNPCategories) {
                    myMultiFit->ComparePulls(icat);
                }
                myMultiFit->CompareNormFactors("");
            }
            if(myMultiFit->fPlotCombCorrMatrix) myMultiFit->PlotCombinedCorrelationMatrix();
            if(myMultiFit->fComparePOI) {
                for (std::size_t i = 0; i < myMultiFit->fPOIs.size(); ++i) {
                    myMultiFit->ComparePOI(myMultiFit->fPOIs.at(i), i);
                }
            }
            if(myMultiFit->fCompareEFT) myMultiFit->CompareEFT();
            if(myMultiFit->fCompareLimits) myMultiFit->CompareLimit();
        }
        //
        if(myMultiFit->fPlotSoverB){
            LOG(INFO) << "Running S over B ..." << std::endl;
            myMultiFit->PlotSummarySoverB();
        }

        if (translateHS3) {
            myMultiFit->TranslateWSToHS3();
        }

        return;
    }

    // proceed if not multi-fit
    std::unique_ptr<TRExFit> myFit = std::make_unique<TRExFit>();
    {
        // initialize config reader
        ConfigReader reader(myFit.get());

        // read the actual config
        int sc = reader.ReadFullConfig(configFile,opt,options);
        if(sc!=0){
            LOG(ERROR) << "Failed to read the config file.\n";
            LOG(ERROR) << "Found " << sc << " errors in the config.\n";
            LOG(ERROR) << "Please check all reported errors.\n";
            exit(EXIT_FAILURE);
        }

        LOG(INFO) << "Successfully read config file.\n";
    }

    std::vector<std::string> binningsToCreate;
    if (readHistograms || readNtuples || rebinAndSmooth) {
        gSystem->mkdir((myFit->fName).c_str());
        for (const auto& ireg : myFit->fRegions) {
            const bool fileExists = std::filesystem::is_regular_file(myFit->fBinningsPath+ireg->fName+".txt");
            if (!fileExists) {
                binningsToCreate.emplace_back(ireg->fName);
            } else {
                if (myFit->fRecreateBinningFiles) {
                    LOG(INFO) << "Binnings file: " << myFit->fBinningsPath << ireg->fName << ".txt exists but you asked to overwrite it\n";
                    std::remove((myFit->fBinningsPath + ireg->fName + ".txt").c_str());
                    binningsToCreate.emplace_back(ireg->fName);
                } else {
                    LOG(INFO) << "Binnings file: " << myFit->fBinningsPath << ireg->fName << ".txt exists and will not overwrite it\n";
                    LOG(INFO) << "If you want to recreate it set \"RecreateBinningFiles: TRUE\" in the Job block or via the command line\n";
                }
            }
        }
        myFit->fBinningFilesToCreate = binningsToCreate;
    }

    if (myFit->fHEPDataFormat) {
        YamlConverter::SubmissionContainer container;
        container.useTables = myFit->fDoTables;
        std::vector<std::string> unfoldingNames;
        for (const auto& iunfold : myFit->fUnfolding) {
            unfoldingNames.emplace_back(iunfold->fName);
        }
        container.unfoldingNames = unfoldingNames;
        container.folder = myFit->fName;
        if (!myFit->fVarNameLH.empty() && myFit->fVarNameLH.at(0) != "all") {
            container.useLikelihoodScan = myFit->fVarNameLH;
        }
        for (const auto& ireg : myFit->fRegions) {
            container.regionNames.emplace_back(ireg->fName);
        }

        if (!myFit->fVarName2DLH.empty()) {
            std::vector<std::string> lhNames;
            for (const auto& i : myFit->fVarName2DLH) {
                lhNames.emplace_back(i.at(0)+"_"+i.at(1));
            }
            container.use2DLikelihoodScan = lhNames;
        }

        YamlConverter converter{};
        converter.WriteHEPDataSubmission(container, myFit->fPOIs);
    }

    // check compatibility between run option and config file
    if(readHistograms && myFit->fInputType!=TRExFit::HIST){
        LOG(ERROR) << "Option \"h\" requested but no HISTO InputType specified in the configuration file. Aborting.\n";
        return;
    }
    if(readNtuples && myFit->fInputType!=TRExFit::NTUP){
        LOG(ERROR) << "Option \"n\" requested but no NTUP InputType specified in the configuration file. Aborting.\n";
        return;
    }

    // -------------------------------------------------------
    myFit->PrintConfigSummary();

    myFit->fHasValidationRegions = myFit->HasValidationRegions();
    myFit->fHasDropBinRegions    = myFit->HasDropBinRegions();

    if (prepareUnfolding) {
        LOG(INFO) << "Running preparation step for unfolding..." << std::endl;
        if (myFit->fFitType != TRExFit::UNFOLDING) {
            LOG(ERROR) << "You want to run unfolding step but the fit type is not set to \"UNFOLDING\". Fix this please.\n";
            return;
        }
        myFit->PrepareUnfolding();
        myFit->CloseInputFiles();
    }

    // Free the memeory
    myFit->fUnfoldingSamples.clear();
    myFit->fUnfoldingSystematics.clear();

    // Fix expression to last bin for normalized-cross-section unfolding, taking bin contents of the truth histogram
    // We need to do it before the h or n steps since it produces some plots where this is needed
    myFit->FixUnfoldingExpressions();

    if(readHistograms){
        LOG(INFO) << "Reading histograms..." << std::endl;
        myFit->CreateRootFiles();
        {
            HistoReader histoReader(myFit.get());

            histoReader.ReadHistograms(binningsToCreate);
        }
        myFit->CreateEFTNominalSampleMaps();

        myFit->CorrectHistograms();
        myFit->MergeSystematics();
        myFit->CreateCustomAsimov();
        myFit->UnfoldingAlternativeAsimov();
        myFit->WriteHistos();
        if(TRExFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TRExFitter::SYSTDATAPLOT)     myFit->DrawSystPlotsSumSamples();
        myFit->CloseInputFiles();
    }
    else if(readNtuples){
        LOG(INFO) << "Reading ntuples..." << std::endl;
        myFit->CreateRootFiles();

        {
            NtupleReader reader(myFit.get());
            reader.ReadNtuples();
        }
        myFit->CreateEFTNominalSampleMaps();
        myFit->CorrectHistograms();
        myFit->MergeSystematics();
        myFit->CreateCustomAsimov();
        myFit->UnfoldingAlternativeAsimov();
        myFit->WriteHistos();
        if(TRExFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TRExFitter::SYSTDATAPLOT)     myFit->DrawSystPlotsSumSamples();
        myFit->CloseInputFiles();
    }
    else{
        if(drawPreFit || drawPostFit || createWorkspace || drawSeparation || rebinAndSmooth || doEFTInputs || prepareTemplate || myFit->fDoNonProfileFit) {
            HistoReader histoReader(myFit.get());

            bool checkHistogram = true;
            if (rebinAndSmooth || doEFTInputs || prepareTemplate) checkHistogram = false;

            histoReader.ReadTRExProducedHistograms(myFit->fFitType == TRExFit::FitType::UNFOLDING, checkHistogram);
            myFit->CreateEFTNominalSampleMaps();
            myFit->CloseInputFiles();
        }
    }


    if(rebinAndSmooth){
        LOG(INFO) << "Rebinning and smoothing..." << std::endl;
        const bool update = myFit->fUpdate;
        myFit->fUpdate = true;
        myFit->CreateRootFiles();
        myFit->fUpdate = update;
        myFit->CorrectHistograms();
        myFit->MergeSystematics();
        myFit->CombineSpecialSystematics();
        myFit->CreateCustomAsimov();
        myFit->UnfoldingAlternativeAsimov();
        myFit->WriteHistos(false);
        if(TRExFitter::SYSTCONTROLPLOTS) myFit->DrawSystPlots();
        if(TRExFitter::SYSTDATAPLOT)     myFit->DrawSystPlotsSumSamples();

        myFit->CloseInputFiles();

        if(drawPreFit || drawPostFit || createWorkspace || drawSeparation || doEFTInputs || prepareTemplate) {
            HistoReader histoReader(myFit.get());

            bool checkHistogram = true;
            if (doEFTInputs || prepareTemplate) checkHistogram = false;

            histoReader.ReadTRExProducedHistograms(myFit->fFitType == TRExFit::FitType::UNFOLDING, checkHistogram);
            myFit->CreateEFTNominalSampleMaps();
            myFit->CloseInputFiles();
        }
    }

    if (prepareUnfolding) {
        if (readHistograms) {
            myFit->ReadRegionBinning();
        }
    } else if (myFit->fWorkspaceFileName.empty()){
        myFit->ReadRegionBinning();
    }

    if (prepareTemplate) {
        myFit->PrepareTemplateMorphing();
    }

    if(drawPreFit || drawPostFit || createWorkspace || doFit || groupedImpact || doRanking || doLimit || doSignificance || doLHscan) {
        myFit->ProcessTemplateMorphingShapeFactors();
        myFit->ProcessEFTShapeFactors();
        myFit->ProcessShapeFactorReparametrization();
    }

    if(createWorkspace){
        LOG(INFO) << "Applying systematics pruning..." << std::endl;
        myFit->SystPruning();
    }

    if (doEFTInputs && !myFit->fEFTConfig.GetSplitSamplesPerBin()) {
        myFit->ProcessEFTInputs(true);
    }

    if (!readHistograms && !readNtuples && !rebinAndSmooth && myFit->fFitType == TRExFit::EFT && myFit->fEFTConfig.GetSplitSamplesPerBin()){
        // From here onwards we can go back to a SPLUSB fit
        myFit->fFitType = TRExFit::SPLUSB;
        if(doEFTInputs){
            // Even if existing fit results are present - force rerun
            LOG(INFO) << "Processing EFT inputs and fitting to get parametrisation..." << std::endl;
            myFit->ProcessEFTInputs(true);
        } else {
            LOG(INFO) << "Processing EFT inputs..." << std::endl;
            myFit->ProcessEFTInputs(false);
        }
    }

    if(createWorkspace){
        LOG(INFO) << "Creating workspace..." << std::endl;
        if(myFit->fFitType == TRExFit::EFT && myFit->fEFTConfig.GetSplitSamplesPerBin()){
            myFit->ProcessEFTInputs(false);
        }
        bool ForFitAndPlot = myFit->fWorkspaceCreationType == TRExFit::WorkspaceCreationType::BOTH;
        if (ForFitAndPlot || myFit->fWorkspaceCreationType == TRExFit::WorkspaceCreationType::FORFIT) {

            /* Create a workspace with all Bins if:

            - The workspace being created is for fit only or for fit+plot but there are no validation regions to plot
            AND
            - DropBins/AutomaticDropBins are used
            */

            myFit->CloseInputFiles();
            myFit->ToRooStats(false, false);

            if ( myFit->fHasDropBinRegions && (!ForFitAndPlot || (!myFit->fHasValidationRegions))){
                myFit->ToRooStats(false, true);
            }
        }
        if (ForFitAndPlot || myFit->fWorkspaceCreationType == TRExFit::WorkspaceCreationType::FORPLOTS) {
            if (myFit->fHasValidationRegions) {
                myFit->ToRooStats(true, false);
            } else {
                if (ForFitAndPlot){
                    LOG(WARNING) << "No VALIDATION regions found, will only produce workspace for fit regions\n";
                    LOG(WARNING) << "To add VALIDATION regions, rerun the \"w\" step with the option WorkspaceCreationType=FORPLOTS\n";
                }
                else{
                    LOG(ERROR) << "No VALIDATION regions found, while asking for WorkspaceCreationType=FORPLOTS\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
        myFit->CloseInputFiles();
    }

    if(doFit){
        LOG(INFO) << "Fitting..." << std::endl;
        myFit->Fit(false);
        myFit->PlotFittedNP();
        myFit->PlotCorrelationMatrix();
        myFit->PlotUnfoldedData(myFit->fName+"/Fits/"+myFit->fInputName+myFit->fSuffix,
                                myFit->fName,
                                myFit->fName+"/RooStats/"+myFit->fInputName+"_combined_"+myFit->fInputName+myFit->fSuffix+"_model.root",
                                "combined");
        myFit->CloseInputFiles();
    }
    if (plotFitResults) {
        myFit->PlotUnfoldedData(myFit->fName+"/Fits/"+myFit->fInputName+myFit->fSuffix,
                                myFit->fName,
                                myFit->fName+"/RooStats/"+myFit->fInputName+"_combined_"+myFit->fInputName+myFit->fSuffix+"_model.root",
                                "combined");
        myFit->CloseInputFiles();
    }
    if (doLHscan){
        LOG(INFO) << "Running LH scan only..." << std::endl;
        myFit->Fit(true);
    }
    if(doRanking){
        LOG(INFO) << "Doing ranking..." << std::endl;
        if (myFit->fRankingOnly!="plot") myFit->ProduceNPRanking(myFit->fRankingOnly);
        if (myFit->fRankingOnly=="all" || myFit->fRankingOnly=="plot") myFit->PlotNPRankingManager();
    }

    if(doLimit){
        LOG(INFO) << "Extracting limit..." << std::endl;
        myFit->GetLimit();
    }

    if(doSignificance){
        LOG(INFO) << "Extracting significance..." << std::endl;
        myFit->GetSignificance();
    }

    if(groupedImpact){
        LOG(INFO) << "Doing grouped systematics impact table..." << std::endl;
        myFit->fDoGroupedSystImpactTable = true;
        if(myFit->fGroupedImpactCategory!="combine") myFit->Fit(false);
        else                                         myFit->BuildGroupedImpactTable();
    }

    std::shared_ptr<TRExPlot> prefit_plot = nullptr;
    std::shared_ptr<TRExPlot> prefit_plot_valid = nullptr;
    if( drawPostFit && TRExFitter::PREFITONPOSTFIT ) {
        drawPreFit = true;
    }

    const std::string log = myFit->fSummaryLogY ?  "log " : "";


    //this is a trick to reduce the time to run destructors of different objects
    auto files = gROOT->GetListOfCleanups()->FindObject("Files");
    gROOT->GetListOfCleanups()->Remove(files);
    if (drawPreFit || drawPostFit){
        std::string resultFileAddition("");
        if (myFit->fHasValidationRegions){
            resultFileAddition =  "_allRegions";
        }
        else if (myFit->fHasDropBinRegions){
            resultFileAddition = "_allBinsFitRegions";
        }

        const std::string wsPath = myFit->fName+"/RooStats/"+myFit->fInputName+resultFileAddition +"_combined_"+myFit->fInputName+myFit->fSuffix+"_model.root";
        if (myFit->fHasValidationRegions){
            // If file does not exist
            if (gSystem->AccessPathName(wsPath.c_str()) == kTRUE){
                LOG(ERROR) << "Cannot find RooStats file containing VALIDATION regions: " << wsPath << "\n";
                LOG(ERROR) << "Did you creeate workspace with  WorkspaceCreationType=FORPLOTS or WorkspaceCreationType=BOTH ?\n";
                exit(EXIT_FAILURE);
            }
        }

        LOG(INFO) << "Drawing plots from the workspace: " << wsPath << "\n";

        // this is needed for blinded bins calculation
        myFit->SetBlindingInRegions(wsPath, drawPostFit);

        if(drawPreFit) {
            LOG(INFO) << "Drawing pre-fit plots..." << std::endl;
            if (TRExFitter::OPTION["PrefitRatioMax"]==2) {
                myFit->DrawAndSaveAll("prefit", wsPath);
                if (myFit->fDoMergedPlot) {
                    if (myFit->fRegionGroups.empty()) {
                        myFit->DrawMergedPlot("prefit");
                    }
                    for (const auto& iregGroup : myFit->fRegionGroups) {
                        myFit->DrawMergedPlot("prefit", iregGroup);
                    }
                }
                if (myFit->fDoSummaryPlot) {
                    prefit_plot       = myFit->DrawSummary(log + "prefit", nullptr, wsPath);
                    prefit_plot_valid = myFit->DrawSummary(log + "valid prefit", nullptr, wsPath);
                }
            } else {
                myFit->DrawAndSaveAll("", wsPath);
                if (myFit->fDoMergedPlot){
                    if (myFit->fRegionGroups.empty()) {
                        myFit->DrawMergedPlot("");
                    }
                    for (const auto& iregGroup : myFit->fRegionGroups) {
                        myFit->DrawMergedPlot("", iregGroup);
                    }
                }
                if (myFit->fDoSummaryPlot) {
                    prefit_plot       = myFit->DrawSummary(log, nullptr, wsPath);
                    prefit_plot_valid = myFit->DrawSummary(log + "valid", nullptr, wsPath);
                }
            }
            if (myFit->fDoTables) {
                myFit->BuildYieldTable("","",wsPath);
                for (const auto& iregGroup : myFit->fRegionGroups) {
                    myFit->BuildYieldTable("", iregGroup, wsPath);
                }
                myFit->PrintSystTables("", wsPath);
            }
            std::size_t nCols = 2;
            std::size_t nRows = 2;
            if (myFit->fRegions.size() > 4) {
                nCols = static_cast<int>(std::sqrt(myFit->fRegions.size()));
                if (std::sqrt(myFit->fRegions.size()) > nCols) nCols++;
                nRows = static_cast<int>(std::sqrt(myFit->fRegions.size()));
                if (nCols*nRows < myFit->fRegions.size()) nRows++;
            }
            if(myFit->fDoSignalRegionsPlot) myFit->DrawSignalRegionsPlot(nCols, nRows);
            if(myFit->fDoPieChartPlot)      myFit->DrawPieChartPlot("pre", nCols, nRows);
        } // end prefit plotting

        if (drawPostFit) {
            LOG(INFO) << "Drawing post-fit plots..." << std::endl;
            myFit->DrawAndSaveAll("post", wsPath);
            if(myFit->fDoMergedPlot) {
                if (myFit->fRegionGroups.empty()) {
                    myFit->DrawMergedPlot("post");
                }
                for (const auto& iregGroup : myFit->fRegionGroups) {
                    myFit->DrawMergedPlot("post", iregGroup);
                }
            }
            if (myFit->fDoSummaryPlot) {
                myFit->DrawSummary(log + "post",      prefit_plot, wsPath);
                myFit->DrawSummary(log + "post valid",prefit_plot_valid, wsPath);
            }
            if (myFit->fDoTables) {
                myFit->BuildYieldTable("post", "", wsPath);
                for (const auto& iregGroup : myFit->fRegionGroups) {
                    myFit->BuildYieldTable("post", iregGroup, wsPath);
                }
                myFit->PrintSystTables("post", wsPath);
            }
            std::size_t nCols = 2;
            std::size_t nRows = 2;
            if (myFit->fRegions.size() > 4) {
                nCols = static_cast<int>(std::sqrt(myFit->fRegions.size()));
                if (std::sqrt(myFit->fRegions.size()) > nCols) nCols++;
                nRows = static_cast<int>(std::sqrt(myFit->fRegions.size()));
                if (nCols*nRows < myFit->fRegions.size()) nRows++;
            }
            if (myFit->fDoPieChartPlot) myFit->DrawPieChartPlot("post", nCols, nRows);
        } // end postfit plotting
    }
    gROOT->GetListOfCleanups()->Add(files);

    if (drawSeparation) {
        LOG(INFO) << "Drawing separation plots..." << std::endl;
        myFit->DrawAndSaveSeparationPlots();
    }

    if (translateHS3) {
        const std::string wsFilePath = myFit->fName+"/RooStats/"+myFit->fInputName+"_combined_"+myFit->fInputName+myFit->fSuffix+"_model.root";
        const std::string wsPath = "combined";
        const std::string outPath = myFit->fName+"/Workspace_HS3.json";
        myFit->TranslateWSToHS3(wsFilePath, wsPath, outPath);
    }

    if (drawPreFit || drawPostFit || createWorkspace || drawSeparation || rebinAndSmooth || groupedImpact || prepareTemplate) myFit->CloseInputFiles();
}

// -------------------------------------------------------
// -------------------------------------------------------
// main function
// -------------------------------------------------------
// -------------------------------------------------------

int main(int argc, char **argv){
    #if TREXFITTER_AVX_OPTIMISATION == false
    RooFit::setBatchCompute("generic");
    #else
    RooFit::setBatchCompute("auto");
    #endif

    std::ifstream in;
    std::string version;
    TString homeArea(gSystem->Getenv("$TREXFITTER_HOME") ? "$TREXFITTER_HOME" : "$DATAPATH");
    #ifdef TREXFITTER_HOME
        homeArea = std::string(STRINGIZE(TREXFITTER_HOME));
    #endif
    TString file = "version.txt";
    auto* expandPath = gSystem->ExpandPathName(homeArea+":./");
    gSystem->FindFile(expandPath,file);
    if (file == "") {
        LOG(ERROR) << "version.txt not found. Are you running from the main working directory?\n";
        exit(EXIT_FAILURE);
    }
    delete [] expandPath;
    in.open(file);
    std::getline(in,version);
    in.close();
    std::cout << "\033[1mTRExFitter version -- " << version << "\n\033[0m";

    if (argc == 2){
        LOG(ERROR) << "You provided 1 parameter, but only 0, 2 or 3 parameters are allowed. Please fix this.'n";
        return 1;
    }

    if (argc > 4) {
        LOG(WARNING) << "You provided " << (argc-1) << " parameters, but only 0, 2 or 3 parameters are allowed.\n";
        LOG(WARNING) << "Ignoring parameters after third parameter.\n";
    }

    std::string opt="h";
    std::string config="config/myFit.config";
    std::string options="";

    if(argc>1) opt     = argv[1];
    if(argc>2) config  = argv[2];
    if(argc>3) options = argv[3];

    // call the function
    FitExample(opt,config,options);

    return 0;
}
