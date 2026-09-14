// Class include

#include "TRExFitter/ConfigReader.h"

// Framework inclused
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/Common.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/TRExFit.h"
#include "TRExFitter/TruthSample.h"
#include "TRExFitter/Unfolding.h"
#include "TRExFitter/UnfoldingSample.h"
#include "TRExFitter/UnfoldingSystematic.h"

// ROOT includes
#include "TColor.h"
#include "TSystem.h"

// c++ includes
#include <algorithm>
#include <iostream>

#define _STRINGIZE(x) #x
#define STRINGIZE(x) _STRINGIZE(x)

//__________________________________________________________________________________
//
ConfigReader::ConfigReader(TRExFit *fitter) :
    fFitter(fitter),
    fParser(new ConfigParser()),
    fAllowWrongRegionSample(false),
    fNonGhostIsSet(false),
    fRunLHscan(true),
    fLHscanVariable(""),
    fLHscanStep(-1),
    fLHscanStepY(-1),
    fHasAtLeastOneValidRegion(false),
    fHasAtLeastOneValidSample(false),
    fHasTemplateMorphing(false),
    fTemplateMorphingDimensions(-1),
    fAddNFsToPOIs(false) {
    LOG(INFO) << "Started reading the config\n";
}

//__________________________________________________________________________________
// Read the full config file
int ConfigReader::ReadFullConfig(const std::string& fileName, const std::string& opt, const std::string& option){
    // initialize ConfigParser for the actual config
    fParser->ReadFile(fileName);

    // initialize checker COnfigParser to cross check the input
    ConfigParser refConfig;
    TString homeArea(gSystem->Getenv("TREXFITTER_HOME") ? "$TREXFITTER_HOME" : "$DATAPATH");
#ifdef TREXFITTER_HOME
    homeArea = std::string(STRINGIZE(TREXFITTER_HOME));
#endif
    TString file = "jobSchema.config";
    auto* expandPath = gSystem->ExpandPathName(homeArea+":./");
    gSystem->FindFile(expandPath,file);
    if (file == "") {
        LOG(ERROR) << "jobSchema.config not found. Are you running from the main working directory?\n";
        exit(EXIT_FAILURE);
    }
    delete [] expandPath;
    refConfig.ReadFile(file.Data());
    int sc = fParser->CheckSyntax(&refConfig);

    if (sc != 0) return sc;

    // syntax of the config is ok
    // read different types of settings
    if (option != ""){
        sc+= ReadCommandLineOptions(option);
    }

    sc+= ReadJobOptions();

    sc+= ReadGeneralOptions();

    sc+= ReadEFTConfig();

    sc+= ReadFitOptions();

    sc+= ReadUnfoldingOptions();

    sc+= ReadLimitOptions();

    sc+= ReadSignificanceOptions();

    sc+= ReadRegionOptions(opt);

    sc+= ReadSampleOptions();

    sc+= ReadNormFactorOptions();

    sc+= ReadShapeFactorOptions();

    sc+= ReadSystOptions();

    sc+= ReadTruthSamples();

    sc+= ReadMorphing();

    sc+= ReadUnfoldingSampleOptions();

    sc+= ReadUnfoldingSystematicOptions();

    sc+= UnfoldingCorrections();

    if (option != ""){
        sc+= ReadCommandLineOptions(option);
    }

    sc+= PostConfig(opt);

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadCommandLineOptions(const std::string& option){
    std::vector< std::string > optVec = Common::Vectorize(option,':');
    std::map< std::string,std::string > optMap;

    int sc(0);

    for(const std::string& iopt : optVec){
        std::vector< std::string > optPair;
        optPair = Common::Vectorize(iopt,'=');
        if (optPair.size() < 2){
            LOG(ERROR) << "Cannot read your command line option, please check this!\n";
            ++sc;
        } else {
            optMap[optPair[0]] = optPair[1];
        }
    }

    // Iterate over map
    for (auto& optMapItr : optMap){
        if (optMapItr.first == "Regions") {
            fOnlyRegions = Common::Vectorize(optMapItr.second,',');
        }
        else if(optMapItr.first == "Samples") {
            fOnlySamples = Common::Vectorize(optMapItr.second,',');
        }
        else if(optMapItr.first == "Systematics") {
            fOnlySystematics = Common::Vectorize(optMapItr.second,',');
        }
        else if(optMapItr.first == "Exclude") {
            fToExclude = Common::Vectorize(optMapItr.second,',');
        }
        else if(optMapItr.first == "Suffix") {
            fFitter->fSuffix = optMapItr.second; // used for input & output  plots, txt files & workspaces - NOT for histograms file
        }
        else if(optMapItr.first == "CustomFunctionsExecutes") {
            fFitter->fCustomFunctionsExecutes = Common::Vectorize(optMapItr.second,';',false);
        }
        else if(optMapItr.first == "SaveSuffix") {
            fFitter->fSaveSuffix = optMapItr.second; // ... and this one for histograms file
        }
        else if(optMapItr.first == "Update" && optMapItr.second!="FALSE") {
            fFitter->fUpdate = true;
        }
        else if(optMapItr.first == "StatOnly" && optMapItr.second!="FALSE") {
            fFitter->fStatOnly = true;
        }
        else if(optMapItr.first == "StatOnlyFit" && optMapItr.second!="FALSE") {
            fFitter->fStatOnlyFit = true;
            fFitter->fUnfoldingShowStat = true;
            LOG(WARNING) << "StatOnlyFit option is deprecated for obtaining stat-only uncertainty. Please, use the covariance matrix breakdown\n";
        }
        else if(optMapItr.first == "Ranking") {
            fFitter->fRankingOnly = optMapItr.second;
        }
        else if(optMapItr.first == "Signal") {
            fOnlySignals.clear();
            fOnlySignals.emplace_back(optMapItr.second);
        }
        else if(optMapItr.first == "Signals") {
            fOnlySignals = Common::Vectorize(optMapItr.second,',');
        }
        else if(optMapItr.first == "NPValuesFromFitResultsFile") {
            fFitter->fFitNPValuesFromFitResultsFile = optMapItr.second;
        }
        else if(optMapItr.first == "FitResultsRootFile") {
            fFitter->fFitResultsRootFile = optMapItr.second;
        }
        else if(optMapItr.first == "FitType"){
            if(optMapItr.second=="SPLUSB") fFitter->SetFitType(TRExFit::SPLUSB);
            else if(optMapItr.second=="BONLY")  fFitter->SetFitType(TRExFit::BONLY);
            else if(optMapItr.second=="UNFOLDING")  fFitter->SetFitType(TRExFit::UNFOLDING);
            else if(optMapItr.second=="EFT")  fFitter->SetFitType(TRExFit::EFT);
            else {
                LOG(ERROR) << "Unrecognised value for command line option FitType, please check!\n";
                ++sc;
            }
        }
        else if(optMapItr.first == "LumiScale"){
            fFitter->fLumiScale = Common::convertStoNum<float>(optMapItr.second);
        }
        else if(optMapItr.first == "BootstrapIdx") {
            fFitter->fBootstrapIdx = Common::convertStoNum<int>(optMapItr.second);
        }
        else if(optMapItr.first == "BootstrapSyst") {
            fFitter->fBootstrapSyst = optMapItr.second;
        }
        else if(optMapItr.first == "BootstrapSample") {
            fFitter->fBootstrapSample = optMapItr.second;
        }
        else if(optMapItr.first == "BootstrapNomHistos") {
            fFitter->fBootstrapNomHistos = optMapItr.second;
        }
        else if(optMapItr.first == "GroupedImpact") {
            fFitter->fGroupedImpactCategory = optMapItr.second;
        }
        else if(optMapItr.first == "OutputDir") {
            fFitter->fDir = Common::RemoveQuotes(optMapItr.second);
            if(fFitter->fDir.back() != '/') fFitter->fDir += '/';
            gSystem->mkdir(fFitter->fDir.c_str());
        }
        else if(optMapItr.first == "Job") {
            if (fFitter->fDir.empty()) {
                fFitter->fName = Common::RemoveQuotes(optMapItr.second);
            } else {
                fFitter->fName = fFitter->fDir + Common::RemoveQuotes(optMapItr.second);
            }
        }
        else if(optMapItr.first == "LimitParamValue") {
            fFitter->fLimitParamValue = atof(optMapItr.second.c_str());
        }
        // By default, do not create or use the LH scan workspace cache
        else if (optMapItr.first == "CreateCache") {
            if (optMapItr.second != "FALSE"){
                fFitter->fCreateCache = true;
            }
        }
        else if (optMapItr.first == "UseCache") {
            if (optMapItr.second != "FALSE"){
                fFitter->fUseCache = true;
            }
        }
        // By default, LHScan should always be run
        else if (optMapItr.first == "RunLHscan" && optMapItr.second == "FALSE") {
            fRunLHscan = false;
        }
        else if(optMapItr.first == "LHscanVariable") {
            fLHscanVariable = optMapItr.second;
        }
        else if(optMapItr.first == "LHscanStep") {
            fLHscanStep = Common::convertStoNum<int>(optMapItr.second);
        }
        else if(optMapItr.first == "LHscanStepY") {
            fFitter->fLHscanStepY = Common::convertStoNum<int>(optMapItr.second);
        }
        else if(optMapItr.first == "FitBlind" && optMapItr.second!="FALSE") {
            fFitter->fFitIsBlind = true;
        }
        else if(optMapItr.first == "BlindedParameters") {
            fFitter->fBlindedParameters = Common::Vectorize(optMapItr.second,',');
        }
        else if(optMapItr.first == "POIAsimovOverride") {
            fFitter->fPOIAsimovOverride = optMapItr.second;
        }
        else if(optMapItr.first == "CustomAsimovOverride") {
            fFitter->fCustomAsimovOverride = optMapItr.second;
        }
        else if(optMapItr.first == "ConstNFs") {
            fConstNFs = Common::Vectorize(optMapItr.second, ',');
        }
        else if(optMapItr.first == "RecreateBinningFiles") {
            fFitter->fRecreateBinningFiles = Common::StringToBoolean(optMapItr.second);
        }
        else if(optMapItr.first == "ToysNpValuesFile") {
            fFitter->fToysNpValuesFile = optMapItr.second;
        }
        else if(optMapItr.first == "NumCPU") {
            try {
                fFitter->fCPU = Common::convertStoNum<int>(optMapItr.second);
                if (fFitter->fCPU <= 0) {
                    LOG(ERROR) << "NumCPU must be a positive integer!\n";
                    ++sc;
                }
            }
            catch (const std::invalid_argument & e) {
                LOG(ERROR) << "NumCPU must be an integer!\n";
                ++sc;
            }
        }
        else if(optMapItr.first == "LimitToysSuffix"){
            fFitter->fLimitToysSuffix = optMapItr.second; // suffix for outputs of toy limit setting - for parellelising
        }
        else if(optMapItr.first == "LimitToysSeed"){
            try {
                fFitter->fLimitToysSeed = Common::convertStoNum<int>(optMapItr.second);
            }
            catch (const std::invalid_argument & e) {
                LOG(ERROR) << "LimitToysSeed must be an int!\n";
                ++sc;
            }
        } else if(optMapItr.first == "SignificanceToysSeed"){
            try {
                fFitter->fSignificanceToysSeed = Common::convertStoNum<int>(optMapItr.second);
            }
            catch (const std::invalid_argument & e) {
                LOG(ERROR) << "SignificanceToysSeed must be an int!\n";
                ++sc;
            }
        } else if(optMapItr.first == "WorkspaceCreationType"){
            std::string param = optMapItr.second;
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "FORFIT") {
                fFitter->fWorkspaceCreationType = TRExFit::WorkspaceCreationType::FORFIT;
            } else if (param == "FORPLOTS") {
                fFitter->fWorkspaceCreationType = TRExFit::WorkspaceCreationType::FORPLOTS;
            } else if (param == "BOTH") {
                fFitter->fWorkspaceCreationType = TRExFit::WorkspaceCreationType::BOTH;
            } else {
                LOG(ERROR) << "Unsupported option for \"WorkspaceCreationType\": " << param << "\n";
                ++sc;
            }
        } else if(optMapItr.first == "ToysSeed"){
            fFitter->fToysSeed = Common::convertStoNum<int>(optMapItr.second);
        } else if(optMapItr.first == "AddNFsToPOIs"){
            fAddNFsToPOIs = Common::StringToBoolean(optMapItr.second);
        } else {
            LOG(ERROR) << "Unrecognised command line option " << optMapItr.first << ", please check this!\n";
            ++sc;
        }
    }

    if( fFitter->fBootstrapSyst!="" && fFitter->fBootstrapSample!=""  ){
        LOG(ERROR) << "Cannot do bootstrap on both Sample and Syst!\n";
        ++sc;
    }

    //
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Running options: \n";
    if(fOnlyRegions.size()>0){
        LOG(INFO) << "  Only these Regions: \n";
        for(const std::string& ireg : fOnlyRegions){
            LOG(INFO) << "    " << ireg << "\n";
        }
    }
    if(fOnlySamples.size()>0){
        LOG(INFO) << "  Only these Samples: \n";
        for(const std::string& isamp : fOnlySamples){
            LOG(INFO) << "    " << isamp << "\n";
        }
    }
    if(fOnlySystematics.size()>0){
        LOG(INFO) << "  Only these Systematics: \n";
        for(const std::string& isyst : fOnlySystematics){
            LOG(INFO) << "    " << isyst << "\n";
        }
    }
    if(fToExclude.size()>0){
        LOG(INFO) << "  Exclude: \n";
        for(const std::string& iexcl : fToExclude){
            LOG(INFO) << "    " << iexcl << "\n";
        }
    }
    if(fOnlySignals.size()>0){
        LOG(INFO) << "  Only Signals: \n";
        std::string toPrint = "";
        for(const auto& s : fOnlySignals) toPrint += " " + s;
        LOG(INFO) << "    " << toPrint << "\n";
    }
    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadJobOptions() {
    int sc(0);

    const ConfigSet *confSet = fParser->GetConfigSet("Job");
    if (confSet == nullptr){
        LOG(ERROR) << "You need to provide JOB settings!\n";
        return 1;
    }

    if (fFitter->fDir == "") {
        // default
        if (fFitter->fName == "MyMeasurement") fFitter->fName = Common::CheckName(confSet->GetValue());
    } else {
        if (fFitter->fName == "MyMeasurement") fFitter->fName = fFitter->fDir + Common::CheckName(confSet->GetValue());
        else fFitter->fName = fFitter->fDir + fFitter->fName;
        gSystem->mkdir(fFitter->fName.c_str());
    }
    fFitter->fInputName = Common::CheckName(confSet->GetValue());

    std::string param = confSet->Get("AllowWrongRegionSample");
    if( param != "") {
        fAllowWrongRegionSample = Common::StringToBoolean(param);
    }

    // Set outputDir
    param = confSet->Get("OutputDir");
    if(param != ""){
        if (fFitter->fDir.size() == 0){
            fFitter->fDir = Common::RemoveQuotes(param);
            if(fFitter->fDir.back() != '/') fFitter->fDir += '/';
            fFitter->fName = fFitter->fDir + fFitter->fName;
            gSystem->mkdir((fFitter->fName).c_str(), true);
        }
    }

    // Set Label
    param = confSet->Get("Label");
    if(param == "") fFitter->fLabel = fFitter->fName;
    else {
        if(param == "NONE") fFitter->fLabel = "";
        else fFitter->fLabel = Common::RemoveQuotes(param);
    }

    // Set POI and unit
    param = confSet->Get("POI");
    if(param!=""){
        std::vector<std::string> pois = Common::Vectorize(param, ',');
        std::vector<std::string> units;
        std::string unit = confSet->Get("POIUnit");
        if(unit!=""){
            units = Common::Vectorize(unit, ',');
            if(pois.size()!=units.size()){
                LOG(WARNING) << "Number of POI units set different from number of actual POI. Ignoring units.\n";
                units.clear();
            }
        }
        for(unsigned int i_poi=0;i_poi<pois.size();i_poi++){
            if(units.size()>i_poi){
                fFitter->AddPOI(Common::CheckName(pois[i_poi]),Common::RemoveQuotes(units[i_poi]));
            }
            else{
                fFitter->AddPOI(Common::CheckName(pois[i_poi]));
            }
        }
    }

    // Set reading option
    param = confSet->Get("ReadFrom");
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);
    if(      param=="HIST" || param=="HISTOGRAMS")  fFitter->fInputType = 0;
    else if( param=="NTUP" || param=="NTUPLES" )    fFitter->fInputType = 1;
    else{
        LOG(ERROR) << "Invalid \"ReadFrom\" argument. Options: \"HIST\", \"NTUP\"\n";
        ++sc;
    }

    // set default MERGEUNDEROVERFLOW
    if (fFitter->fFitType == TRExFit::FitType::UNFOLDING) {
        TRExFitter::MERGEUNDEROVERFLOW = false;
    } else {
        if(fFitter->fInputType==0)      TRExFitter::MERGEUNDEROVERFLOW = false;
        else if(fFitter->fInputType==1) TRExFitter::MERGEUNDEROVERFLOW = true;
    }

    // Set MergeUnderOverFlow from config
    param = confSet->Get("MergeUnderOverFlow");
    if(param!=""){
        TRExFitter::MERGEUNDEROVERFLOW = Common::StringToBoolean(param);
    }

    // Set paths for the unfolding preprocessing
    param = confSet->Get("ResponseMatrixName");
    if (param != "") {
        fFitter->fResponseMatrixNames.clear();
        fFitter->fResponseMatrixNames.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("ResponseMatrixNames");
    if (param != "") {
        fFitter->fResponseMatrixNames = Common::Vectorize(param, ',');
    }

    param = confSet->Get("ResponseMatrixFile");
    if (param != "") {
        fFitter->fResponseMatrixFiles.clear();
        fFitter->fResponseMatrixFiles.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("ResponseMatrixFiles");
    if (param != "") {
        fFitter->fResponseMatrixFiles = Common::Vectorize(param, ',');
    }

    param = confSet->Get("ResponseMatrixPath");
    if (param != "") {
        fFitter->fResponseMatrixPaths.clear();
        fFitter->fResponseMatrixPaths.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("ResponseMatrixPaths");
    if (param != "") {
        fFitter->fResponseMatrixPaths = Common::Vectorize(param, ',');
    }

    param = confSet->Get("ResponseMatrixNameNominal");
    if(param!=""){
      fFitter->fResponseMatrixNamesNominal.clear();
      fFitter->fResponseMatrixNamesNominal.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("AcceptanceName");
    if (param != "") {
        fFitter->fAcceptanceNames.clear();
        fFitter->fAcceptanceNames.emplace_back(Common::RemoveQuotes(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceNames");
    if (param != "") {
        fFitter->fAcceptanceNames = Common::Vectorize(param, ',');
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceFile");
    if (param != "") {
        fFitter->fAcceptanceFiles.clear();
        fFitter->fAcceptanceFiles.emplace_back(Common::RemoveQuotes(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceFiles");
    if (param != "") {
        fFitter->fAcceptanceFiles = Common::Vectorize(param, ',');
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptancePath");
    if (param != "") {
        fFitter->fAcceptancePaths.clear();
        fFitter->fAcceptancePaths.emplace_back(Common::RemoveQuotes(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptancePaths");
    if (param != "") {
        fFitter->fAcceptancePaths = Common::Vectorize(param, ',');
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("AcceptanceNameNominal");
    if(param!=""){
        fFitter->fAcceptanceNamesNominal.clear();
        fFitter->fAcceptanceNamesNominal.emplace_back(Common::RemoveQuotes(param));
        fFitter->fHasAcceptance = true;
    }

    param = confSet->Get("SelectionEffName");
    if (param != "") {
        fFitter->fSelectionEffNames.clear();
        fFitter->fSelectionEffNames.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("SelectionEffNames");
    if (param != "") {
        fFitter->fSelectionEffNames = Common::Vectorize(param, ',');
    }

    param = confSet->Get("SelectionEffFile");
    if (param != "") {
        fFitter->fSelectionEffFiles.clear();
        fFitter->fSelectionEffFiles.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("SelectionEffFiles");
    if (param != "") {
        fFitter->fSelectionEffFiles = Common::Vectorize(param, ',');
    }

    param = confSet->Get("SelectionEffPath");
    if (param != "") {
        fFitter->fSelectionEffPaths.clear();
        fFitter->fSelectionEffPaths.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("SelectionEffPaths");
    if (param != "") {
        fFitter->fSelectionEffPaths = Common::Vectorize(param, ',');
    }

    param = confSet->Get("SelectionEffNameNominal");
    if(param!=""){
      fFitter->fSelectionEffNamesNominal.clear();
      fFitter->fSelectionEffNamesNominal.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("MigrationName");
    if (param != "") {
        fFitter->fMigrationNames.clear();
        fFitter->fMigrationNames.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("MigrationNames");
    if (param != "") {
        fFitter->fMigrationNames = Common::Vectorize(param, ',');
    }

    param = confSet->Get("MigrationFile");
    if (param != "") {
        fFitter->fMigrationFiles.clear();
        fFitter->fMigrationFiles.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("MigrationFiles");
    if (param != "") {
        fFitter->fMigrationFiles = Common::Vectorize(param, ',');
    }

    param = confSet->Get("MigrationPath");
    if (param != "") {
        fFitter->fMigrationPaths.clear();
        fFitter->fMigrationPaths.emplace_back(Common::RemoveQuotes(param));
    }

    param = confSet->Get("MigrationPaths");
    if (param != "") {
        fFitter->fMigrationPaths = Common::Vectorize(param, ',');
    }

    param = confSet->Get("MigrationNameNominal");
    if(param!=""){
      fFitter->fMigrationNamesNominal.clear();
      fFitter->fMigrationNamesNominal.emplace_back(Common::RemoveQuotes(param));
    }

    // Set paths
    // HIST option only
    if(fFitter->fInputType==0){
        if (confSet->Get("MCweight") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'MCweight:' option, ignoring.\n";
        }
        if (confSet->Get("Selection") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'Selection:' option, ignoring.\n";
        }
        if (confSet->Get("NtupleName") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'NtupleName:' option, ignoring.\n";
        }
        if (confSet->Get("NtupleFile") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'NtupleFile:' option, ignoring.\n";
        }
        if (confSet->Get("NtuplePath") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'NtuplePath:' option, ignoring.\n";
        }
        if (confSet->Get("NtuplePaths") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'NtuplePaths:' option, ignoring.\n";
        }
        if (confSet->Get("FriendPath") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'FriendPath:' option, ignoring.\n";
        }
        if (confSet->Get("FriendPaths") != ""){
            LOG(WARNING) << "You specified HIST option but you provided 'FriendPaths:' option, ignoring.\n";
        }
        //
        param = confSet->Get("HistoName");
        if(param!=""){
            fFitter->fHistoNames.clear();
            fFitter->fHistoNames.emplace_back( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("HistoNames");
        if(param!=""){
            fFitter->fHistoNames = Common::Vectorize( param,',' );
        }
        param = confSet->Get("HistoNameNominal");
        if(param!=""){
          fFitter->fHistoNamesNominal.clear();
          fFitter->fHistoNamesNominal.emplace_back( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("HistoFile");
        if(param!=""){
            fFitter->fHistoFiles.clear();
            fFitter->fHistoFiles.emplace_back( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("HistoFiles");
        if(param!=""){
            fFitter->fHistoFiles = Common::Vectorize( param,',' );
        }
        param = confSet->Get("HistoPath");
        if(param!=""){
            fFitter->AddHistoPath( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("HistoPaths");
        if(param!=""){
            fFitter->fHistoPaths = Common::Vectorize( param,',' );
        }
    }
    // Setting for NTUP only
    if(fFitter->fInputType==1){
        if (confSet->Get("HistoName") != ""){
            LOG(WARNING) << "You specified NTUP option but you provided 'HistoName:' option, ignoring.\n";
        }
        if (confSet->Get("HistoFile") != ""){
            LOG(WARNING) << "You specified NTUP option but you provided 'HistoFile:' option, ignoring.\n";
        }
        if (confSet->Get("HistoPath") != ""){
            LOG(WARNING) << "You specified NTUP option but you provided 'HistoPath:' option, ignoring.\n";
        }
        if (confSet->Get("HistoPaths") != ""){
            LOG(WARNING) << "You specified NTUP option but you provided 'HistoPaths:' option, ignoring.\n";
        }
        //
        param = confSet->Get("NtupleName");
        if(param!=""){
            Common::CheckTreeName(param);
            fFitter->SetNtupleName( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("NtupleNames");
        if(param!=""){
            Common::CheckTreeName(param);
            fFitter->fNtupleNames = Common::Vectorize( param,',' );
        }
        param = confSet->Get("NtupleFile");
        if( param != "" ){
            fFitter->SetNtupleFile( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("NtupleFiles");
        if(param!=""){
            fFitter->fNtupleFiles = Common::Vectorize( param,',' );
        }
        param = confSet->Get("NtuplePath");
        if( param != "" ) {
            fFitter->AddNtuplePath( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("NtuplePaths");
        if( param != "" ){
            std::vector<std::string> paths = Common::Vectorize( param,',' );
            for(const std::string& ipath : paths){
                fFitter->AddNtuplePath( ipath );
            }
        }
        param = confSet->Get("FriendPath");
        if( param != "" ) {
            fFitter->AddFriendPath( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("FriendPaths");
        if( param != "" ){
            std::vector<std::string> fpaths = Common::Vectorize( param,',' );
            for(const std::string& ipath : fpaths){
                fFitter->AddFriendPath( ipath );
            }
        }
        param = confSet->Get("MCweight");
        if(param!=""){
            fFitter->SetMCweight( Common::RemoveQuotes(param) );
        }
        param = confSet->Get("Selection");
        if(param!=""){
            fFitter->SetSelection( Common::RemoveQuotes(param) );
        }
    }

    // Set lumi
    param = confSet->Get("Lumi");
    if( param != "" ) fFitter->SetLumi( atof(param.c_str()) );

    // Set LumiScale
    param = confSet->Get("LumiScale");
    if( param != "" ){
        LOG(WARNING) << "\"LumiScale\" is only done for quick tests since it is inefficient.\n";
        LOG(WARNING) << "To normalize all the samples to the luminosity, use \"Lumi\" instead.\n";
        fFitter->fLumiScale = atof(param.c_str());
    }

    // Set Smoothing option
    param = confSet->Get("SmoothingOption");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "MAXVARIATION" ) fFitter->fSmoothOption = HistoTools::SmoothOption::MAXVARIATION;
        else if (param == "TTBARRESONANCE") fFitter->fSmoothOption = HistoTools::SmoothOption::TTBARRESONANCE;
        else if (param == "COMMONTOOLSMOOTHMONOTONIC") fFitter->fSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHMONOTONIC;
        else if (param == "COMMONTOOLSMOOTHPARABOLIC") fFitter->fSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHPARABOLIC;
        else if (param == "TCHANNEL") fFitter->fSmoothOption = HistoTools::SmoothOption::TCHANNEL;
        else if (param == "KERNELRATIOUNIFORM") fFitter->fSmoothOption = HistoTools::SmoothOption::KERNELRATIOUNIFORM;
        else if (param == "KERNELDELTAGAUSS") fFitter->fSmoothOption = HistoTools::SmoothOption::KERNELDELTAGAUSS;
        else if (param == "KERNELRATIOGAUSS") fFitter->fSmoothOption = HistoTools::SmoothOption::KERNELRATIOGAUSS;
        else {
            LOG(WARNING) << "You specified 'SmoothingOption' option but you did not provide valid input. Using default (MAXVARIATION)\n";
            fFitter->fSmoothOption = HistoTools::SmoothOption::MAXVARIATION;
        }
    }

    // Set SystPruningShape
    param = confSet->Get("SystPruningShape");
    if( param != "") fFitter->fThresholdSystPruning_Shape = atof(param.c_str());

    // Set SystPruningNorm
    param = confSet->Get("SystPruningNorm");
    if( param != "")  fFitter->fThresholdSystPruning_Normalisation = atof(param.c_str());

    // Set SystLarge
    param = confSet->Get("SystLarge");
    if( param != "")  fFitter->fThresholdSystLarge = atof(param.c_str());

    // Set IntCode
    param = confSet->Get("IntCode");
    if( param != "")  fFitter->fIntCode = std::atoi(param.c_str());

    // Set MCstatThreshold
    param = confSet->Get("MCstatThreshold");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param=="NONE")  fFitter->SetStatErrorConfig( false, 0. );
        else{
            fFitter->SetStatErrorConfig( true, atof(param.c_str()));
        }
    }
    else{
        fFitter->SetStatErrorConfig( true, 0. );
    }

    //Set MCstatConstraint
    param = confSet->Get("MCstatConstraint");
    if( param != "") fFitter->fStatErrCons = Common::RemoveQuotes(param);

    // Set UseGammaPulls
    param = confSet->Get("UseGammaPulls");
    if( param != ""){
        fFitter->fUseGammaPulls = Common::StringToBoolean(param);
    }

    // plotting options are in special function
    sc+= SetJobPlot(confSet);

    // Set TableOptions
    param = confSet->Get("TableOptions");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        fFitter->fTableOptions = Common::RemoveQuotes(param);
    }

    // Set CorrelationThreshold
    param = confSet->Get("CorrelationThreshold");
    if( param != ""){
        // set it if it was not set already
        if (TRExFitter::CORRELATIONTHRESHOLD < 0){
            TRExFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
        }
    }

    // Set HistoChecks
    param = confSet->Get("HistoChecks");
    if(param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "NOCRASH" ){
            TRExFitter::HISTOCHECKCRASH = false;
        } else {
            LOG(WARNING) << "You specified 'HistoChecks' option but you did not set it to NOCRASH.\n";
        }
    }

    // Set LumiLabel
    param = confSet->Get("LumiLabel");
    if( param != "") fFitter->fLumiLabel = Common::RemoveQuotes(param);

    // Set CmeLabel
    param = confSet->Get("CmeLabel");
    if( param != "") fFitter->fCmeLabel = Common::RemoveQuotes(param);

    // Set SplitHistoFiles
    param = confSet->Get("SplitHistoFiles");
    if( param != ""){
        TRExFitter::SPLITHISTOFILES = Common::StringToBoolean(param);
    }

    // Set BlindingThreshold"
    param = confSet->Get("BlindingThreshold");
    if( param != ""){
        fFitter->fBlindingThreshold = atof(param.c_str());
    }

    // Set BlindingType
    param = confSet->Get("BlindingType");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param=="SOVERB" )                fFitter->fBlindingType = Common::SOVERB;
        else if ( param=="SOVERSPLUSB" )     fFitter->fBlindingType = Common::SOVERSPLUSB;
        else if ( param=="SOVERSQRTB" )      fFitter->fBlindingType = Common::SOVERSQRTB;
        else if ( param=="SOVERSQRTSPLUSB" ) fFitter->fBlindingType = Common::SOVERSQRTSPLUSB;
        else {
            LOG(WARNING) << "You specified 'BlindingType' option but did not provide a valid setting. Using default (SOVERB)\n";
            fFitter->fBlindingType = Common::SOVERB;
        }
    }

    // Set RankingMaxNP
    param = confSet->Get("RankingMaxNP");
    if( param != ""){
        fFitter->fRankingMaxNP = atoi(param.c_str());
    }

    // Set ImageFormat
    param = confSet->Get("ImageFormat");
    if( param != ""){
        std::vector<std::string> tmp = Common::Vectorize(param,',');
        if (tmp.size() > 0) fFitter->fImageFormat = tmp.at(0);
        else {
            LOG(ERROR) << "You specified 'ImageFormat' option but we cannot split the setting. Please check\n";
            ++sc;
        }
        TRExFitter::IMAGEFORMAT = tmp;
    }

    // Set StatOnly
    param = confSet->Get("StatOnly");
    if( param != "" ){
        fFitter->fStatOnly = Common::StringToBoolean(param);
    }

    // Set FixNPforStatOnly
    param = confSet->Get("FixNPforStatOnly");
    if( param != "" ){
        fFitter->fFixNPforStatOnlyFit = Common::StringToBoolean(param);
        LOG(WARNING) << "StatOnlyFit option is deprecated for obtaining stat-only uncertainty. Please, use the covariance matrix breakdown\n";
    }

    // Set InputFolder
    param = confSet->Get("InputFolder");
    if( param != "" ){
        fFitter->fInputFolder = Common::RemoveQuotes(param);
        gSystem->mkdir((fFitter->fInputFolder).c_str());
        fFitter->fBinningsPath = fFitter->fInputFolder + "/Binnings/";
    }
    else{
        fFitter->fBinningsPath = fFitter->fName + "/Histograms/Binnings/";
    }
    gSystem->mkdir((fFitter->fBinningsPath).c_str(), true);


    // Set InputName
    param = confSet->Get("InputName");
    if( param != "" ){
        fFitter->fInputName = Common::RemoveQuotes(param);
    }

    // Set WorkspaceFileName
    param = confSet->Get("WorkspaceFileName");
    if( param != "" ){
        fFitter->fWorkspaceFileName = Common::RemoveQuotes(param);
    }

    // Set KeepPruning
    param = confSet->Get("KeepPruning");
    if( param != "" ){
        fFitter->fKeepPruning = Common::StringToBoolean(param);
    }

    // Set PlotLabel
    param = confSet->Get("PlotLabel");
    if( param != "" ){
        std::string tmp = param;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
        if (tmp == "EMPTY") {
            fFitter->fPlotLabel = "";
        } else {
            fFitter->fPlotLabel = Common::RemoveQuotes(param);
        }
    }

    // Set CleanTables
    param = confSet->Get("CleanTables");
    if( param != "" ){
        fFitter->fCleanTables = Common::StringToBoolean(param);
    }

    // Set SystCategoryTables
    param = confSet->Get("SystCategoryTables");
    if( param != "" ){
        fFitter->fSystCategoryTables = Common::StringToBoolean(param);
    }

    // Set Suffix
    param = confSet->Get("Suffix");
    if( param != "" ){
        fFitter->fSuffix = Common::RemoveQuotes(param);
    }

    // Set SaveSuffix
    param = confSet->Get("SaveSuffix");
    if( param != "" ){
        fFitter->fSaveSuffix = Common::RemoveQuotes(param);
    }

    // Set SummaryPrefix
    param = confSet->Get("SummaryPrefix");
    if( param != "" ){
        fFitter->fSummaryPrefix = Common::RemoveQuotes(param);
    }

    // Set HideNP
    param = confSet->Get("HideNP");
    if( param != "" ){
        fFitter->fVarNameHide = Common::Vectorize(param,',');
    }

    // Set RegionGroups
    param = confSet->Get("RegionGroups");
    if( param != "" ) {
        std::vector<std::string> groups = Common::Vectorize(param,',');
        for(const std::string& igroup : groups) fFitter->fRegionGroups.emplace_back(igroup);
    }

    // Set KeepPrefitBlindedBins
    param = confSet->Get("KeepPrefitBlindedBins");
    if( param != "" ){
        fFitter->fKeepPrefitBlindedBins = Common::StringToBoolean(param);
    }

    // Set CustomAsimov
    param = confSet->Get("CustomAsimov");
    if ( fFitter->fCustomAsimovOverride != "" ) {
        param = fFitter->fCustomAsimovOverride;
        LOG(INFO) << "CustomAsimov values overriden by command line values\n";
    }
    if( param != "" ){
        fFitter->fCustomAsimov = Common::RemoveQuotes(param);
        param = confSet->Get("WriteCustomAsimovToWS");
        if( param != "" ){
          fFitter->fWriteCustomAsimovToWS = Common::StringToBoolean(param);
        }
    }

    // Set GetChi2
    param = confSet->Get("GetChi2");
    if( param != "" ){ // can be TRUE, SYST+STAT, STAT-ONLY... (if it contains STAT and no SYST => stat-only, ptherwise stat+syst)
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "TRUE" ){
            fFitter->fGetChi2 = 2;
        }
        else if( param.find("SYST")!=std::string::npos ){
            fFitter->fGetChi2 = 2;
        }
        else if( param.find("STAT")!=std::string::npos ){
            fFitter->fGetChi2 = 1;
        } else {
            LOG(ERROR) << "You specified 'GetChi2' option but you did not provide valid option. Check this\n";
            ++sc;
        }
    }

    // Set DoTables
    param = confSet->Get("DoTables");
    if( param != "" ){
        fFitter->fDoTables = Common::StringToBoolean(param);
    }

    // Set CustomFunctions
    param = confSet->Get("CustomFunctions");
    if( param != "" ) {
        fFitter->fCustomFunctions = Common::Vectorize(param,',');
    }

    // Set CustomIncludePaths
    param = confSet->Get("CustomIncludePaths");
    if( param != "" ) {
        fFitter->fCustomIncludePaths = Common::Vectorize(param,',');
    }

    // Set CustomFunctions Executes
    param = confSet->Get("CustomFunctionsExecutes");
    if( param != "" ) {
      fFitter->fCustomFunctionsExecutes = Common::Vectorize(Common::RemoveQuotes(param),';');
    }

    // Add Aliases
    param = confSet->Get("AddAliases");
    if( param != "" ) {
        fFitter->fAddAliases = Common::Vectorize(Common::RemoveQuotes(param),';');
    }

    // Set Bootstrap
    param = confSet->Get("Bootstrap");
    if( param != "" ){
        fFitter->fBootstrap = Common::RemoveQuotes(param);
    }

    // Set Bootstrap systematic
    param = confSet->Get("BootstrapSyst");
    if( param != "" ){
        fFitter->fBootstrapSyst = Common::RemoveQuotes(param);
    }

    // Set Bootstrap sample
    param = confSet->Get("BootstrapSample");
    if( param != "" ){
        fFitter->fBootstrapSample = Common::RemoveQuotes(param);
    }

    param = confSet->Get("BootstrapNomHistos");
    if( param != "" ){
        fFitter->fBootstrapNomHistos = Common::RemoveQuotes(param);
    }

    // Set DecorrSuff
    param = confSet->Get("DecorrSuff");
    if( param != ""){
        fFitter->fDecorrSuff = Common::RemoveQuotes(param);
    }

    // Set DecorrSysts
    param = confSet->Get("DecorrSysts");
    if( param != ""){
        fFitter->fDecorrSysts = Common::Vectorize(param,',');
    }

    // Set SmoothMorphingTemplates
    param = confSet->Get("SmoothMorphingTemplates");
    if ( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        fFitter->fSmoothMorphingTemplates = param;
    }

    // Set POIPrecision
    param = confSet->Get("POIPrecision");
    if( param != ""){
        fFitter->fPOIPrecision = stoi(param);
        if (fFitter->fPOIPrecision < 1 || fFitter->fPOIPrecision > 5){
            LOG(WARNING) << "Parameter POIPrecision has value smaller than 1 or larger than 5. Using default (2).\n";
            fFitter->fPOIPrecision = 2;
        }
    }

    // Set PropagateSystsForMorphing
    param = confSet->Get("PropagateSystsForMorphing");
    if( param != ""){
        fFitter->fPropagateSystsForMorphing = Common::StringToBoolean(param);
    }

    // Set UsePDGRounding
    param = confSet->Get("UsePDGRounding");
    if ( param != ""){
        fFitter->fUsePDGRoundingTxt = Common::StringToBoolean(param);
        fFitter->fUsePDGRoundingTex = Common::StringToBoolean(param);
    }
    param = confSet->Get("UsePDGRoundingTxt");
    if ( param != ""){
        fFitter->fUsePDGRoundingTxt = Common::StringToBoolean(param);
    }
    param = confSet->Get("UsePDGRoundingTex");
    if ( param != ""){
        fFitter->fUsePDGRoundingTex = Common::StringToBoolean(param);
    }

    // Set RankingPOIName
    param = confSet->Get("RankingPOIName");
    if( param != ""){
        fFitter->fRankingPOIName = Common::Vectorize(param, ',');
    }

    // Set RankingUpperAxisNdivision
    param = confSet->Get("RankingUpperAxisNdivision");
    if( param != ""){
        const auto tmp = Common::Vectorize(param, ',');
        std::vector<int> v;
        for (const auto& i : tmp) {
            v.emplace_back(Common::convertStoNum<int>(i));
        }
        fFitter->fRankingUpperAxisNdivision = v;
    }

    // Set RankingPOIAxisScale
    param = confSet->Get("RankingPOIAxisScale");
    if( param != ""){
        const auto tmp = Common::Vectorize(param, ',');
        std::vector<double> v;
        for (const auto& i : tmp) {
            v.emplace_back(Common::convertStoNum<float>(i));
        }
        fFitter->fRankingPOIAxisScale = v;
    }

    // Set PruningType
    param = confSet->Get("PruningType");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "SEPARATESAMPLE") {
            fFitter->fPruningType = TRExFit::SEPARATESAMPLE;
        } else if (param == "BACKGROUNDREFERENCE") {
            fFitter->fPruningType = TRExFit::BACKGROUNDREFERENCE;
        } else if (param == "COMBINEDREFERENCE") {
            fFitter->fPruningType = TRExFit::COMBINEDREFERENCE;
        } else if (param == "COMBINEDSIGNAL") {
            fFitter->fPruningType = TRExFit::COMBINEDSIGNAL;
        } else {
            LOG(WARNING) << "You specified PruningType option but did not provide valid parameter. Using default (SEPARATESAMPLE)\n";
            fFitter->fPruningType = TRExFit::SEPARATESAMPLE;
        }
    }

    param = confSet->Get("DoSystNormalizationPlots");
    if (param != "") {
        fFitter->fDoSystNormalizationPlots = Common::StringToBoolean(param);
    }

    // Set PrePostFitCanvasSize
    param = confSet->Get("PrePostFitCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Common::Vectorize(param, ',');
        if (tmp.size() != 2){
            LOG(ERROR) << "You specified PrePostFitCanvasSize option but did not provide 2 parameters. \n";
            ++sc;
        } else {
            const int x = Common::convertStoNum<int>(tmp.at(0));
            const int y = Common::convertStoNum<int>(tmp.at(1));

            if (x <= 100 || y <= 100){
                LOG(WARNING) << "You specified PrePostFitCanvasSize option but at least one parameter is <= 100.\n";
            }
            fFitter->fPrePostFitCanvasSize.emplace_back(x);
            fFitter->fPrePostFitCanvasSize.emplace_back(y);
        }
    }

    // Set SummaryCanvasSize
    param = confSet->Get("SummaryCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Common::Vectorize(param, ',');
        if (tmp.size() != 2){
            LOG(ERROR) << "You specified SummaryCanvasSize option but did not provide 2 parameters. Ignoring.\n";
            ++sc;
        } else {
            const int x = Common::convertStoNum<int>(tmp.at(0));
            const int y = Common::convertStoNum<int>(tmp.at(1));

            if (x <= 100 || y <= 100){
                LOG(WARNING) << "You specified SummaryCanvasSize option but at least one parameter is <= 100.\n";
            }
            fFitter->fSummaryCanvasSize.emplace_back(x);
            fFitter->fSummaryCanvasSize.emplace_back(y);
        }
    }

    // Set MergeCanvasSize
    param = confSet->Get("MergeCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Common::Vectorize(param, ',');
        if (tmp.size() != 2){
            LOG(ERROR) << "You specified MergeCanvasSize option but did not provide 2 parameters.\n";
            ++sc;
        } else {
            const int x = Common::convertStoNum<int>(tmp.at(0));
            const int y = Common::convertStoNum<int>(tmp.at(1));

            if (x <= 100 || y <= 100){
                LOG(WARNING) << "You specified MergeCanvasSize option but at least one parameter is <= 100.\n";
            }
            fFitter->fMergeCanvasSize.emplace_back(x);
            fFitter->fMergeCanvasSize.emplace_back(y);
        }
    }

    // Set PieChartCanvasSize
    param = confSet->Get("PieChartCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Common::Vectorize(param, ',');
        if (tmp.size() != 2){
            LOG(ERROR) << "You specified PieChartCanvasSize option but did not provide 2 parameters.\n";
            ++sc;
        } else {
            const int x = Common::convertStoNum<int>(tmp.at(0));
            const int y = Common::convertStoNum<int>(tmp.at(1));

            if (x <= 100 || y <= 100){
                LOG(WARNING) << "You specified PieChartCanvasSize option but at least one parameter is <= 100.\n";
            }
            fFitter->fPieChartCanvasSize.emplace_back(x);
            fFitter->fPieChartCanvasSize.emplace_back(y);
        }
    }

    // Set NPRankingCanvasSize
    param = confSet->Get("NPRankingCanvasSize");
    if( param != ""){
        std::vector<std::string> tmp = Common::Vectorize(param, ',');
        if (tmp.size() != 2){
            LOG(ERROR) << "You specified NPRankingCanvasSize option but did not provide 2 parameters.\n";
            ++sc;
        } else {
            const int x = Common::convertStoNum<int>(tmp.at(0));
            const int y = Common::convertStoNum<int>(tmp.at(1));

            if (x <= 100 || y <= 100){
                LOG(WARNING) << "You specified NPRankingCanvasSize option but at least one parameter is <= 100.\n";
            }
            fFitter->fNPRankingCanvasSize.emplace_back(x);
            fFitter->fNPRankingCanvasSize.emplace_back(y);
        }

    }
    // Set SeparationPlot
    param = confSet->Get("SeparationPlot");
    if( param != ""){
      fFitter->fSeparationPlot = Common::Vectorize(param, ',');
    }

    // Set ReorderNPs
    param = confSet->Get("ReorderNPs");
    if( param != "" ){
        fFitter->fReorderNPs = Common::StringToBoolean(param);
    }

    // Set BlindSRs
    param = confSet->Get("BlindSRs");
    if( param != "" ){
        fFitter->fBlindSRs = Common::StringToBoolean(param);
    }

    // Set HEPDataFormat
    param = confSet->Get("HEPDataFormat");
    if( param != "" ){
        fFitter->fHEPDataFormat = Common::StringToBoolean(param);
    }

    // Set UheppFormat
    param = confSet->Get("UheppFormat");
    if( param != "" ){
        fFitter->fUheppFormat = Common::StringToBoolean(param);
    }

    // Set CombinerFormat
    param = confSet->Get("CombinerFormat");
    if( param != "" ){
        fFitter->fCombinerFormat = Common::StringToBoolean(param);
    }

    // Set AlternativeShapeHistFactory
    param = confSet->Get("AlternativeShapeHistFactory");
    if( param != "" ){
        fFitter->fAlternativeShapeHistFactory = Common::StringToBoolean(param);
    }

    // Set RemoveLargeSyst
    param = confSet->Get("RemoveLargeSyst");
    if (param != "") {
        fFitter->fRemoveLargeSyst = Common::StringToBoolean(param);
    }

    // Set RemoveSystOnEmptySample
    param = confSet->Get("RemoveSystOnEmptySample");
    if( param != "" ){
        fFitter-> fRemoveSystOnEmptySample = Common::StringToBoolean(param);
    }

    // Set ShowValidationPruning
    param = confSet->Get("ShowValidationPruning");
    if( param != "" ){
        fFitter->fValidationPruning = Common::StringToBoolean(param);
    }

    // Set DoExtendedCovariances
    param = confSet->Get("DoExtendedCovariances");
    if(param != ""){
        fFitter->fDoExtendedCovariances = Common::StringToBoolean(param);
    }

    param = confSet->Get("UseFolderStructure");
    if(param != ""){
        TRExFitter::USEFOLDERSTRUCTURE = Common::StringToBoolean(param);
    }

    param = confSet->Get("UseHEPDataRounding");
    if(param != ""){
        TRExFitter::USEHEPDATAROUNDING = Common::StringToBoolean(param);
    }

    param = confSet->Get("UnfoldingShowStatOnlyError");
    if(param != ""){
        fFitter->fUnfoldingShowStat = Common::StringToBoolean(param);
    }

    param = confSet->Get("UnfoldingShowUncertaintyBreakdown");
    if(param != ""){
        fFitter->fUnfoldingShowUncertaintyBreakdown = Common::StringToBoolean(param);
    }

    param = confSet->Get("UnfoldingUncertaintyBreakdownTotal");
    if(param != ""){
        fFitter->fUnfoldingUncertaintyBreakdownTotal = Common::StringToBoolean(param);
    }

    param = confSet->Get("ApplyGammaCorrection");
    if (param != ""){
        fFitter->fApplyGammaCorrection = Common::StringToBoolean(param);
    }

    param = confSet->Get("RecreateBinningFiles");
    if (param != ""){
        fFitter->fRecreateBinningFiles = Common::StringToBoolean(param);
    }

    param = confSet->Get("NoPrePostFitLabel");
    if (param != ""){
        fFitter->fNoPrePostFitLabel = Common::StringToBoolean(param);
    }

    param = confSet->Get("WorkspaceCreationType");
    if (param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "FORFIT") {
            fFitter->fWorkspaceCreationType = TRExFit::WorkspaceCreationType::FORFIT;
        } else if (param == "FORPLOTS") {
            fFitter->fWorkspaceCreationType = TRExFit::WorkspaceCreationType::FORPLOTS;
        } else if (param == "BOTH") {
            fFitter->fWorkspaceCreationType = TRExFit::WorkspaceCreationType::BOTH;
        }
    }

    param = confSet->Get("ExperimentLabel");
    if (param != ""){
        TRExFitter::EXPERIMENT_LABEL = Common::RemoveQuotes(param);
    }

    param = confSet->Get("ProducePerRegionWS");
    if (param != ""){
        fFitter->fProducePerRegionWS = Common::StringToBoolean(param);
    }

    param = confSet->Get("DrawPruningPlot");
    if (param != ""){
        fFitter->fDrawPruningPlot = Common::StringToBoolean(param);
    }

    param = confSet->Get("PreFitLabel");
    if (param != ""){
        fFitter->fPreFitLabel = Common::RemoveQuotes(param);
    }

    param = confSet->Get("PostFitLabel");
    if (param != ""){
        fFitter->fPostFitLabel = Common::RemoveQuotes(param);
    }

    param = confSet->Get("ToysForErrorBand");
    if (param != ""){
        fFitter->fToysForErrorBand = Common::convertStoNum<int>(param);
    }

    // success
    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::SetJobPlot(const ConfigSet *confSet){
    int sc(0);

    // Plot option
    std::string param = confSet->Get("PlotOptions");
    std::vector<std::string> vec;
    if( param != ""){
        vec = Common::Vectorize(Common::RemoveQuotes(param),',');
        if( std::find(vec.begin(), vec.end(), "YIELDS") !=vec.end() )  TRExFitter::SHOWYIELDS     = true;
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TRExFitter::SHOWSTACKSIG   = false;
        if( std::find(vec.begin(), vec.end(), "SKIPSIG")!=vec.end() )  TRExFitter::ADDSTACKSIG    = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TRExFitter::SHOWNORMSIG    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TRExFitter::SHOWOVERLAYSIG = true;
        if( std::find(vec.begin(), vec.end(), "LEFT")   !=vec.end() )  TRExFitter::LEGENDLEFT     = true; // forces leg entry to be left-aligned when adding yields to legend
        if( std::find(vec.begin(), vec.end(), "RIGHT")  !=vec.end() )  TRExFitter::LEGENDRIGHT    = true; // forces leg entry to be right-aligned even when no yields in legend
        if( std::find(vec.begin(), vec.end(), "CHI2")   !=vec.end() )  TRExFitter::SHOWCHI2       = true;
        if( std::find(vec.begin(), vec.end(), "PREFITONPOSTFIT")   !=vec.end() )  TRExFitter::PREFITONPOSTFIT= true;
        if( std::find(vec.begin(), vec.end(), "XERR") !=vec.end() )    TRExFitter::REMOVEXERRORS  = false;
        if( std::find(vec.begin(), vec.end(), "OPRATIO") !=vec.end() ) TRExFitter::OPRATIO        = true;
        if( std::find(vec.begin(), vec.end(), "NORATIO") !=vec.end() ) TRExFitter::NORATIO        = true;
        if( std::find(vec.begin(), vec.end(), "KEEPDATAERRORS") !=vec.end() ) TRExFitter::KEEPDATAERRORS = true;
    }

    param = confSet->Get("OVERSIGSCALE");
    if( param != "" && TRExFitter::SHOWOVERLAYSIG == true){
       TRExFitter::SHOWOVERLAYSIG_CUSTOMSCALE = atof(param.c_str());
    }


    // Set PlotOptionsSummary
    param = confSet->Get("PlotOptionsSummary");
    if( param != ""){
        vec = Common::Vectorize(Common::RemoveQuotes(param),',');
        if( std::find(vec.begin(), vec.end(), "NOSIG")  !=vec.end() )  TRExFitter::SHOWSTACKSIG_SUMMARY   = false;
        if( std::find(vec.begin(), vec.end(), "NORMSIG")!=vec.end() )  TRExFitter::SHOWNORMSIG_SUMMARY    = true;
        if( std::find(vec.begin(), vec.end(), "OVERSIG")!=vec.end() )  TRExFitter::SHOWOVERLAYSIG_SUMMARY = true;
    }
    else{
        LOG(DEBUG) << "PlotOptionsSummary not specified setting Summary values to 'PlotOptions'\n";
        TRExFitter::SHOWSTACKSIG_SUMMARY   = TRExFitter::SHOWSTACKSIG    ;
        TRExFitter::SHOWNORMSIG_SUMMARY    = TRExFitter::SHOWNORMSIG     ;
        TRExFitter::SHOWOVERLAYSIG_SUMMARY = TRExFitter::SHOWOVERLAYSIG  ;
    }

    // Set SystControlPlots
    param = confSet->Get("SystControlPlots");
    if( param != ""){
        TRExFitter::SYSTCONTROLPLOTS = Common::StringToBoolean(param);
    }

    // Set SystDataPlots
    param = confSet->Get("SystDataPlots");
    if( param != "" ){
        TRExFitter::SYSTDATAPLOT = Common::StringToBoolean(param);
        fFitter->fSystDataPlot_upFrame = Common::StringToBoolean(param);
    }

    // Set SystErrorBars
    param = confSet->Get("SystErrorBars");
    if( param != ""){
        TRExFitter::SYSTERRORBARS = Common::StringToBoolean(param);
    }

    // Set GuessMCStatEmptyBins
    param = confSet->Get("GuessMCStatEmptyBins");
    if( param != ""){
        TRExFitter::GUESSMCSTATERROR = Common::StringToBoolean(param);
    }

    // Set CorrectNormForNegativeIntegral
    param = confSet->Get("CorrectNormForNegativeIntegral");
    if( param != ""){
        TRExFitter::CORRECTNORMFORNEGATIVEINTEGRAL = Common::StringToBoolean(param);
    }

    // Set SuppressNegativeBinWarnings
    param = confSet->Get("SuppressNegativeBinWarnings");
     if( param != ""){
        fFitter->fSuppressNegativeBinWarnings = Common::StringToBoolean(param);
    }

    // Set SignalRegionsPlot
    param = confSet->Get("SignalRegionsPlot");
    if(param != ""){
        fFitter->fRegionsToPlot = Common::Vectorize(param,',');
    }

    // Set SummaryPlotRegions
    param = confSet->Get("SummaryPlotRegions");
    if(param != ""){
        fFitter->fSummaryPlotRegions = Common::Vectorize(param,',');
    }

    // Set SummaryPlotLabels
    param = confSet->Get("SummaryPlotLabels");
    if(param != ""){
        fFitter->fSummaryPlotLabels = Common::Vectorize(param,',');
    }

    // Set SummaryPlotValidationRegions
    param = confSet->Get("SummaryPlotValidationRegions");
    if(param != ""){
        fFitter->fSummaryPlotValidationRegions = Common::Vectorize(param,',');
    }

    // Set SummaryPlotValidationLabels
    param = confSet->Get("SummaryPlotValidationLabels");
    if(param != ""){
        fFitter->fSummaryPlotValidationLabels = Common::Vectorize(param,',');
    }

    // Set SummaryPlotYmin
    param = confSet->Get("SummaryPlotYmin");
    if(param != "") fFitter->fYmin = atof(param.c_str());

    // Set SummaryPlotYmax
    param = confSet->Get("SummaryPlotYmax");
    if(param != "") fFitter->fYmax = atof(param.c_str());

    // Set SummaryLogY
    param = confSet->Get("SummaryLogY");
    if(param != "") {
        fFitter->fSummaryLogY = Common::StringToBoolean(param);
    }

    // Set RatioYmin
    param = confSet->Get("RatioYmin");
    if(param != "") {
        fFitter->fRatioYmin = atof(param.c_str());
        fFitter->fRatioYminPostFit = fFitter->fRatioYmin;
    }

    // Set RatioYmax
    param = confSet->Get("RatioYmax");
    if(param != ""){
        fFitter->fRatioYmax = atof(param.c_str());
        fFitter->fRatioYmaxPostFit = fFitter->fRatioYmax;
    }

    // Set RatioYminPostFit
    param = confSet->Get("RatioYminPostFit");
    if(param != "") fFitter->fRatioYminPostFit = atof(param.c_str());

    // Set RatioYmaxPostFit
    param = confSet->Get("RatioYmaxPostFit");
    if(param != "") fFitter->fRatioYmaxPostFit = atof(param.c_str());

    // Set RatioTitle
    param = confSet->Get("RatioYtitle");
    if(param != "") fFitter->fRatioYtitle = Common::RemoveQuotes(param);

    // Set RatioTitle
    param = confSet->Get("RatioType");
    if(param != "") {
        param = Common::RemoveQuotes(param);
        if (param == "DATAOVERMC") fFitter->fRatioType = TRExPlot::RATIOTYPE::DATAOVERMC;
        else if (param == "DATAOVERBKG") fFitter->fRatioType = TRExPlot::RATIOTYPE::DATAOVERB;
        else if (param == "SOVERB") fFitter->fRatioType = TRExPlot::RATIOTYPE::SOVERB;
        else if (param == "SOVERSQRT(B)") fFitter->fRatioType = TRExPlot::RATIOTYPE::SOVERSQRTB;
        else if (param == "SOVERSQRT(S+B)") fFitter->fRatioType = TRExPlot::RATIOTYPE::SOVERSQRTSPLUSB;
        else if (param == "DATAOVERBWITHS") fFitter->fRatioType = TRExPlot::RATIOTYPE::DATAOVERBWITHS;
    }

    // Set Label and Legend position
    param = confSet->Get("LabelX");
    if(param != "") fFitter->fLabelX = atof(param.c_str());
    param = confSet->Get("LabelY");
    if(param != "") fFitter->fLabelY = atof(param.c_str());
    param = confSet->Get("LegendX1");
    if(param != "") fFitter->fLegendX1 = atof(param.c_str());
    param = confSet->Get("LegendX2");
    if(param != "") fFitter->fLegendX2 = atof(param.c_str());
    param = confSet->Get("LegendY");
    if(param != "") fFitter->fLegendY = atof(param.c_str());

    // Set Label and Legend position for Summary
    param = confSet->Get("LabelXSummary");
    if(param != "") fFitter->fLabelXSummary = atof(param.c_str());
    param = confSet->Get("LabelYSummary");
    if(param != "") fFitter->fLabelYSummary = atof(param.c_str());
    param = confSet->Get("LegendX1Summary");
    if(param != "") fFitter->fLegendX1Summary = atof(param.c_str());
    param = confSet->Get("LegendX2Summary");
    if(param != "") fFitter->fLegendX2Summary = atof(param.c_str());
    param = confSet->Get("LegendYSummary");
    if(param != "") fFitter->fLegendYSummary = atof(param.c_str());

    // Set Label and Legend position for Merge
    param = confSet->Get("LabelXMerge");
    if(param != "") fFitter->fLabelXMerge = atof(param.c_str());
    param = confSet->Get("LabelYMerge");
    if(param != "") fFitter->fLabelYMerge = atof(param.c_str());
    param = confSet->Get("LegendX1Merge");
    if(param != "") fFitter->fLegendX1Merge = atof(param.c_str());
    param = confSet->Get("LegendX2Merge");
    if(param != "") fFitter->fLegendX2Merge = atof(param.c_str());
    param = confSet->Get("LegendYMerge");
    if(param != "") fFitter->fLegendYMerge = atof(param.c_str());

    // Set LegendNColumns
    param = confSet->Get("LegendNColumns");
    if(param != "") fFitter->fLegendNColumns = atoi(param.c_str());
    param = confSet->Get("LegendNColumnsMerge");
    if(param != "") fFitter->fLegendNColumnsMerge = atoi(param.c_str());
    param = confSet->Get("LegendNColumnsSummary");
    if(param != "") fFitter->fLegendNColumnsSummary = atoi(param.c_str());

    // Set DoSummaryPlot
    param = confSet->Get("DoSummaryPlot");
    if( param != "" ){
        fFitter->fDoSummaryPlot = Common::StringToBoolean(param);
    }

    // Set DoMergedPlot
    param = confSet->Get("DoMergedPlot");
    if( param != "" ){
        fFitter->fDoMergedPlot = Common::StringToBoolean(param);
    }

    // Set DoSignalRegionsPlot
    param = confSet->Get("DoSignalRegionsPlot");
    if( param != "" ){
        fFitter->fDoSignalRegionsPlot = Common::StringToBoolean(param);
    }
    param = confSet->Get("DoPieChartPlot");
    if( param != "" ){
        fFitter->fDoPieChartPlot = Common::StringToBoolean(param);
    }

    // Set RankingPlot
    param = confSet->Get("RankingPlot");
    if( param != ""){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        fFitter->fRankingPlot = Common::RemoveQuotes(param);
    }

    // exclude a template from morphing
    param = confSet->Get("ExcludeFromMorphing");
    if( param != ""){
        fFitter->fExcludeFromMorphing = Common::CheckName(param);
    }

    param = confSet->Get("MaxNtupleEvents");
    if(param != "") fFitter->fDebugNev = atoi(param.c_str());

    param = confSet->Get("PruningShapeOption");
    if(param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "MAXBIN") {
            fFitter->fPruningShapeOption = PruningUtil::SHAPEOPTION::MAXBIN;
        } else if (param == "KSTEST") {
            fFitter->fPruningShapeOption = PruningUtil::SHAPEOPTION::KSTEST;
        } else {
            LOG(WARNING) << "You specified 'PruningShapeOption' option but did not provide valid parameter. Using default (MAXBIN)\n";
            fFitter->fPruningShapeOption = PruningUtil::SHAPEOPTION::MAXBIN;
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadGeneralOptions(){
    int sc(0);
    const ConfigSet* confSet = fParser->GetConfigSet("Options");
    if (confSet != nullptr){
        for(int i=0; i < confSet->GetN(); i++){
            if(confSet->GetConfigValue(i) != ""){
                TRExFitter::OPTION[confSet->GetConfigName(i)] = atof(confSet->GetConfigValue(i).c_str());
            }
        }
    } else {
        LOG(DEBUG) << "You do not have 'Options' option in the config. It is ok, we just want to let you know.\n";
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadFitOptions() {
    int sc(0);

    const ConfigSet *confSet = fParser->GetConfigSet("Fit");
    if (confSet == nullptr){
        LOG(INFO) << "You do not have Fit option in the config. It is ok, we just want to let you know.\n";
        return 0; // it is ok to not have Fit set up
    }

    // Set FitType
    std::string param = confSet->Get("FitType");
    if( param != "" && fFitter->fFitType == TRExFit::UNDEFINED ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "SPLUSB" ){
            fFitter->SetFitType(TRExFit::SPLUSB);
        }
        else if( param == "BONLY" ){
            fFitter->SetFitType(TRExFit::BONLY);
        }
        else if( param == "UNFOLDING" ){
            fFitter->SetFitType(TRExFit::UNFOLDING);
        }
        else if( param == "EFT" ){
            fFitter->SetFitType(TRExFit::EFT);
        }
        else{
            LOG(ERROR) << "Unknown FitType argument : " << confSet->Get("FitType") << "\n";
            ++sc;
        }
    }
    else if( fFitter->fFitType == TRExFit::UNDEFINED ){
        LOG(WARNING) << "Setting default fit Type SPLUSB\n";
        fFitter->SetFitType(TRExFit::SPLUSB);
    }

    // Set FitRegion
    param = confSet->Get("FitRegion");
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);
    if( param != "" ){
        if( param == "CRONLY" ){
            fFitter->SetFitRegion(TRExFit::CRONLY);
        }
        else if( param == "CRSR" ){
            fFitter->SetFitRegion(TRExFit::CRSR);
        }
        else{
            LOG(ERROR) << "Unknown FitRegion argument : " << confSet->Get("FitRegion") << "\n";
            ++sc;
        }
    }

    // Set FitBlind
    param = confSet->Get("FitBlind");
    if( param != "" ){
        fFitter->fFitIsBlind = Common::StringToBoolean(param);
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if ( fFitter->fPOIAsimovOverride != "" ) {
        param = fFitter->fPOIAsimovOverride;
        LOG(INFO) << "POIAsimov Values overriden by command line values\n";
    }
    if( param != "" ){
        if (fFitter->fPOIs.empty()) {
            LOG(ERROR) << "POI is not set, you cannot ask for POIAsimov\n";
            ++sc;
        }
        std::vector<std::string> temp_vec = Common::Vectorize(param,',',false);
        for(const std::string& iPOI : temp_vec){
            std::vector<std::string> poi_value = Common::Vectorize(iPOI,'@');
            // single POIAsimov defined
            if(poi_value.size()==1 && temp_vec.size()==1){
                fFitter->fFitPOIAsimov[fFitter->fPOIs[0]] = Common::convertStoNum<float>(poi_value[0].c_str());
            }
            else if(poi_value.size()==1 && temp_vec.size()!=1){
                LOG(ERROR) << "Incorrect setting of 'POIAsimov'... changed recently by setting value with the @ symbol instead of the : symbol\n";
                LOG(ERROR) << "example in the config file: POIAsimov:POI1@2.0,POI2@-5,POI3@0.5\n";
                LOG(ERROR) << "example using command line: trex-fitter wf <config> POIAsimovOverride=POI1@2.0,POI2@-5,POI3@0.5\n";
                ++sc;
            }
            else if(poi_value.size()==2){
                fFitter->fFitPOIAsimov.insert(std::pair <std::string, double>(poi_value[0], Common::convertStoNum<float>(poi_value[1])));
            } else {
                LOG(WARNING) << "You specified 'POIAsimov' option but did not provide 2 parameters for each POI which is expected. Ignoring\n";
                LOG(WARNING) << "example in the config file: POIAsimov:POI1@2.0,POI2@-5,POI3@0.5\n";
                LOG(WARNING) << "example using command line: trex-fitter wf <config> POIAsimovOverride=POI1@2.0,POI2@-5,POI3@0.5\n";
            }
        }
    }


    // Set NPValues
    param = confSet->Get("NPValues");
    if( param != "" ){
        std::vector < std::string > temp_vec = Common::Vectorize(param,',',false);
        for(const std::string& iNP : temp_vec){
            std::vector < std::string > np_value = Common::Vectorize(iNP,':');
            if(np_value.size()==2){
                fFitter->fFitNPValues.insert( std::pair < std::string, double >( np_value[0], atof(np_value[1].c_str()) ) );
            } else {
                LOG(WARNING) << "You specified 'NPValues' option but did not provide 2 parameters for each NP which is expected. Ignoring\n";
            }
        }
    }

    // Get NPValues from fit reults (ROOT file)
    param = confSet->Get("NPValuesFromFitResultsFile");
    if( param != "" ){
        fFitter->fFitNPValuesFromFitResultsFile = Common::RemoveQuotes(param);
    }

    // Set InjectGlobalObservables
    param = confSet->Get("InjectGlobalObservables");
    if( param != ""){
        fFitter->fInjectGlobalObservables = Common::StringToBoolean(param);
    }

    // Set FixNPs
    param = confSet->Get("FixNPs");
    if( param != "" ){
        std::vector < std::string > temp_fixedNPs = Common::Vectorize(param,',',false);
        for(const std::string& iNP : temp_fixedNPs){
            std::vector < std::string > fixed_nps = Common::Vectorize(iNP,':');
            if(fixed_nps.size()==2){
                fFitter->fFitFixedNPs.insert( std::pair < std::string, double >( fixed_nps[0], std::atof(fixed_nps[1].c_str()) ) );
            } else {
                LOG(WARNING) << "You specified 'FixNPs' option but did not provide 2 parameters for each NP which is expected. Ignoring\n";
            }
        }
    }

    // Set doLHscan
    param = confSet->Get("doLHscan");
    if( fRunLHscan && param != "" ){
        if (fLHscanVariable==""){
            fFitter->fVarNameLH = Common::Vectorize(param,',');
        } else {
            fFitter->fVarNameLH.emplace_back(fLHscanVariable);
        }
    }

    // Set do2DLHscan
    param = confSet->Get("do2DLHscan");
    if( fRunLHscan && param != "" ){
        const std::vector<std::string> tmp = Common::Vectorize(param,':');
        for (const auto& ivec : tmp){
            const std::vector<std::string> v = Common::Vectorize(ivec,',');
            if (v.size() == 2){
                fFitter->fVarName2DLH.emplace_back(v);
            } else {
                LOG(WARNING) << "You specified 'do2DLHscan' option but did not provide correct input. Ignoring\n";
            }
        }
    }

    // Set LHscanMin
    param = confSet->Get("LHscanMin");
    if ( fRunLHscan && param != "" ) {
        if (fFitter->fVarNameLH.size() == 0 && fFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMin' option but did not set doLHscan. Ignoring\n";
        } else {
            fFitter->fLHscanMin = Common::convertStoNum<float>(param);
        }
    }

    // Set LHscanMax
    param = confSet->Get("LHscanMax");
    if ( fRunLHscan && param != "" ) {
        if (fFitter->fVarNameLH.size() == 0 && fFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMax' option but did not set doLHscan. Ignoring\n";
        } else {
            fFitter->fLHscanMax = Common::convertStoNum<float>(param);
        }
    }

    // Set LHscanSteps
    param = confSet->Get("LHscanSteps");
    if ( param != "" ) {
        if (fFitter->fVarNameLH.size() == 0 && fFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanSteps' option but did not set doLHscan. Ignoring\n";
        } else {
            fFitter->fLHscanSteps = Common::convertStoNum<int>(param);
            if(fFitter->fLHscanSteps < 3 || fFitter->fLHscanSteps > 500){
                LOG(WARNING) << "LHscanSteps is smaller than 3 or larger than 500, setting to default (30)\n";
                fFitter->fLHscanSteps = 30;
            }
        }
    }

    // Set LHscanStep (-1 = process all points, otherwise process the single point only)
    if (fLHscanStep != -1) {
        if (fLHscanStep < -1 || fLHscanStep >= fFitter->fLHscanSteps) {
            LOG(ERROR) << "Invalid command line option for LHscanStep, please check!\n";
            ++sc;
        }
        fFitter->fLHscanStep = fLHscanStep;
    }

    // Set LHscanMin for second variable
    param = confSet->Get("LHscanMinY");
    if ( fRunLHscan && param != "" ) {
        if (fFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMinY' option but did not set do2DLHscan. Ignoring\n";
        } else {
            fFitter->fLHscanMinY = Common::convertStoNum<float>(param);
        }
    }

    // Set LHscanMax for second variable
    param = confSet->Get("LHscanMaxY");
    if ( fRunLHscan && param != "" ) {
        if (fFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMaxY' option but did not set do2DLHscan.\n";
        } else {
            fFitter->fLHscanMaxY = Common::convertStoNum<float>(param);
        }
    }

    // Set LHscanSteps for second variable
    param = confSet->Get("LHscanStepsY");
    if ( param != "" ) {
        if (fFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanStepsY' option but did not set do2DLHscan.\n";
        } else {
            fFitter->fLHscanStepsY = Common::convertStoNum<int>(param);
            if(fFitter->fLHscanStepsY < 3 || fFitter->fLHscanStepsY > 100){
                LOG(WARNING) << "LHscanStepsY is smaller than 3 or larger than 100, setting to LHscanSteps\n";
                fFitter->fLHscanStepsY = fFitter->fLHscanSteps;
            }
        }
    }
    else {
        fFitter->fLHscanStepsY = fFitter->fLHscanSteps;
    }

    // Set LHscanStep (-1 = process all points, otherwise process the single point only)
    if (fLHscanStepY != -1) {
        if (fLHscanStepY < -1 || fLHscanStepY >= fFitter->fLHscanStepsY) {
            ++sc;
            LOG(ERROR) << "Invalid option for LHscanStepY, please check!\n";
        }
        fFitter->fLHscanStepY = fLHscanStepY;
    }

    // Set UseMinos
    param = confSet->Get("UseMinos");
    if( param != "" ){
        fFitter->fVarNameMinos = Common::Vectorize(param,',');
    }

    // Set SetRandomInitialNPval
    param = confSet->Get("SetRandomInitialNPval");
    if( param != ""){
        fFitter->fUseRnd = true;
        fFitter->fRndRange = std::atof(param.c_str());
    }

    // Set SetRandomInitialNPvalSeed
    param = confSet->Get("SetRandomInitialNPvalSeed");
    if( param != ""){
        fFitter->fRndSeed = std::atol(param.c_str());
    }

    // Set NumCPU
    param = confSet->Get("NumCPU");
    if( param != "" ){
        fFitter->fCPU = std::atoi( param.c_str());
        if (fFitter->fCPU <= 0) {
            LOG(ERROR) << "NumCPU must be a positive integer!\n";
            ++sc;
        }
    }

    // Set StatOnlyFit
    param = confSet->Get("StatOnlyFit");
    if( param != "" ){
        fFitter->fStatOnlyFit = Common::StringToBoolean(param);
        LOG(WARNING) << "StatOnlyFit option is deprecated for obtaining stat-only uncertainty. Please, use the covariance matrix breakdown\n";
    }

    // Set DoNonProfileFit
    param = confSet->Get("DoNonProfileFit");
    if( param != "" ){
        fFitter->fDoNonProfileFit = Common::StringToBoolean(param);
    }

    // Set NonProfileFitSystThreshold
    param = confSet->Get("NonProfileFitSystThreshold");
    if( param != "" ){
        fFitter->fNonProfileFitSystThreshold = std::atof(param.c_str());
    }

    // Set FitToys
    param = confSet->Get("FitToys");
    if( param != "" ){
        fFitter->fFitToys = std::atoi( param.c_str());
    }

    // Set ToysSeed
    param = confSet->Get("ToysSeed");
    if( param != "" ){
        fFitter->fToysSeed = std::atoi( param.c_str());
    }

    // Set ToysHistoNbins
    param = confSet->Get("ToysHistoNbins");
    if( param != "" ){
        fFitter->fToysHistoNbins = std::atoi( param.c_str());
        if (fFitter->fToysHistoNbins < 2){
            LOG(ERROR) << "Number of bins for toys is < 2\n";
            ++sc;
        }
    }

    param = confSet->Get("ToysStatOutput");
    if (param != ""){
        fFitter->fToysStatOutput = Common::StringToBoolean(param);
    }

    param = confSet->Get("ToysNpValuesFile");
    if (param != "" && fFitter->fToysNpValuesFile != "") {
        fFitter->fToysNpValuesFile = param;
    }

    // Set TemplateInterpolationOption
    param = confSet->Get("TemplateInterpolationOption");
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "LINEAR"){
            fFitter->fTemplateInterpolationOption = TRExFit::LINEAR;
        } else if (param == "SMOOTHLINEAR"){
            fFitter->fTemplateInterpolationOption = TRExFit::SMOOTHLINEAR;
        } else if (param == "SQUAREROOT"){
            fFitter->fTemplateInterpolationOption = TRExFit::SQUAREROOT;
        } else {
            LOG(WARNING) << "You specified 'TemplateInterpolationOption' option but did not provide valid parameter. Using default (LINEAR)\n";
            fFitter->fTemplateInterpolationOption = TRExFit::LINEAR;
        }
    }

    // Set BlindedParameters
    param = confSet->Get("BlindedParameters");
    if( param != "" ){
        const std::vector<std::string> tmp = Common::Vectorize( param,',');
        fFitter->fBlindedParameters = tmp;
    }

    // Set Minuit2 fit strategy
    param = confSet->Get("FitStrategy");
    if (param != "") {
        fFitter->fFitStrategy = Common::convertStoNum<int>(param);
        if (fFitter->fFitStrategy > 3) {
            LOG(WARNING) << "Fit strategy value is > 3, setting to default instead (-1)\n";
            fFitter->fFitStrategy = -1;
        }
    }

    param = confSet->Get("MaximumNumberFCNcalls");
    if (param != "") {
        const int max = Common::convertStoNum<int>(param);
        if (max < 0) {
            LOG(WARNING) << "\"MaximumNumberFCNcalls\" is < 0. Ignoring the option\n";
        } else {
            fFitter->fMaximumNumberFCNcalls = max;
        }
    }

    // Set BinnedLikelihoodOptimization
    param = confSet->Get("BinnedLikelihoodOptimization");
    if (param != "") {
        fFitter->fBinnedLikelihood = Common::StringToBoolean(param);
    }

    param = confSet->Get("UsePOISinRanking");
    if (param != "") {
        fFitter->fUsePOISinRanking = Common::StringToBoolean(param);
    }

    param = confSet->Get("ShiftGlobalObservablesInRanking");
    if (param != ""){
        fFitter->fShiftGlobalObservablesInRanking = Common::StringToBoolean(param);
    }

    param = confSet->Get("UseHesseBeforeMigrad");
    if (param != "") {
        fFitter->fUseHesseBeforeMigrad = Common::StringToBoolean(param);
    }

    param = confSet->Get("UseNLLwithoutOffsetInLHscan");
    if (param != "") {
        fFitter->fUseNllInLHscan = Common::StringToBoolean(param);
    }

    param = confSet->Get("DataWeighted");
    if (param != "") {
        fFitter->fDataWeighted = Common::StringToBoolean(param);
    }

    param = confSet->Get("RegularizationType");
    if (param != "") {
        fFitter->fRegularizationType = Common::convertStoNum<int>(param);
    }

    param = confSet->Get("MultiStageFit");
    if (param != "") {
        fFitter->fSpeedUpFit = Common::StringToBoolean(param);
    }

    param = confSet->Get("NPCutOff");
    if (param != "") {
        fFitter->fNPCutOff = Common::convertStoNum<float>(param);
        if (fFitter->fNPCutOff < 0. || fFitter->fNPCutOff > 1.0) {
            LOG(WARNING) << "NPCutOff is < 0 or > 1, this is not a valid option, setting to default (0.3)\n";
            fFitter->fNPCutOff = 0.3;
        }
    }

    param = confSet->Get("UseRegularBinning");
    if (param != "") {
        fFitter->fUseRebinned = Common::StringToBoolean(param);
    }

    param = confSet->Get("ErrorSigma");
    if (param != "") {
        fFitter->fErrorSigma = Common::convertStoNum<float>(param);
    }

    param = confSet->Get("FitResultsRootFile");
    if (param != "") {
        fFitter->fFitResultsRootFile = param;
    }

    param = confSet->Get("UseHesse");
    if (param != "") {
        fFitter->fUseHesse = Common::StringToBoolean(param);
    }

    param = confSet->Get("ToleranceScale");
    if (param != "") {
        fFitter->fToleranceScale = Common::convertStoNum<double>(param);
    }

    param = confSet->Get("ErrorDecomposition");
    if (param != "") {
        fFitter->fErrorDecomposition = Common::StringToBoolean(param);
    }

    param = confSet->Get("UseAutoDiff");
    if (param != "") {
        fFitter->fUseAutoDiff = Common::StringToBoolean(param);
    }

    param = confSet->Get("AddNFsToPOIs");
    if (param != "") {
        fAddNFsToPOIs = Common::StringToBoolean(param);
    }

    param = confSet->Get("FitToysRandomInitialNormFactorValues");
    if (param != "") {
        const std::vector<std::string> singleParam = Common::Vectorize(param, ',');
        for (const auto& iparam : singleParam) {
            const std::vector<std::string> values = Common::Vectorize(iparam, ':');
            if (values.size() != 3) {
                LOG(ERROR) << "Wrong input for \"FitToysRandomInitialNormFactorValues\"\n";
                return ++sc;
            }
            const std::string& nfName = values.at(0);
            const double min = Common::convertStoNum<double>(values.at(1));
            const double max = Common::convertStoNum<double>(values.at(2));

            auto itr = fFitter->fToyRandomPOIStartingValues.find(nfName);
            if (itr != fFitter->fToyRandomPOIStartingValues.end()) {
                LOG(WARNING) << nfName << " already set in the config for \"FitToysRandomInitialNormFactorValues\", ignoring\n";
            } else {
                fFitter->fToyRandomPOIStartingValues.insert({nfName, std::make_pair(min, max)});
            }
        }
    }

    param = confSet->Get("FitToysLHscanCondition");
    if (param != "") {
        const std::vector<std::string> singleParam = Common::Vectorize(param, ',');
        for (const auto& iparam : singleParam) {
            const std::vector<std::string> values = Common::Vectorize(iparam, ':');
            if (values.size() != 3) {
                LOG(ERROR) << "Wrong input for \"FitToysLHscanCondition\"\n";
                return ++sc;
            }
            const std::string& nfName = values.at(0);
            const double min = Common::convertStoNum<double>(values.at(1));
            const double max = Common::convertStoNum<double>(values.at(2));

            fFitter->fToysLHScanCondition.insert({nfName, std::make_pair(min, max)});
        }
    }

    param = confSet->Get("ToysLHScanForAll");
    if (param != "") {
        fFitter->fToysLHScanForAll = Common::StringToBoolean(param);
    }

    param = confSet->Get("ToysOnlyStatFluctuation");
    if (param != "") {
        fFitter->fToysOnlyStatFluctuation = Common::StringToBoolean(param);
    }

    param = confSet->Get("NLLOffset");
    if (param != "") {
        // Ensure lowercase to match RooFit documentation
        std::transform(param.begin(), param.end(), param.begin(), ::tolower);
        fFitter->fNLLOffset = param;
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadLimitOptions(){

    int sc(0);

    const ConfigSet* confSet = fParser->GetConfigSet("Limit");
    if (confSet == nullptr){
        LOG(DEBUG) << "You do not have Limit option in the config. It is ok, we just want to let you know.\n";
        return 0; // it is ok to not have Fit set up
    }

    // Set POI to be used for limit
    std::string param = confSet->Get("POI");
    if( param != "" ){
        fFitter->fPOIforLimit = Common::CheckName(param);
        // add this to the list of POIs
        fFitter->AddPOI(fFitter->fPOIforLimit);
    }

    // Set LimitType
    param = confSet->Get("LimitType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "ASYMPTOTIC" ){
            fFitter->SetLimitType(TRExFit::LimitType::ASYMPTOTIC);
        }
        else if( param == "TOYS" ){
            fFitter->SetLimitType(TRExFit::LimitType::TOYS);
        }
        else{
            LOG(ERROR) << "Unknown LimitType argument : " << confSet->Get("LimitType") << "\n";
            ++sc;
        }
    }

    // Set LimitBlind
    param = confSet->Get("LimitBlind");
    if( param != "" ){
        fFitter->fLimitIsBlind = Common::StringToBoolean(param);
    }

    // Set SignalInjection
    param = confSet->Get("SignalInjection");
    if( param != "" ){
        fFitter->fSignalInjection = Common::StringToBoolean(param);
    }

    // Set SignalInjectionValue
    param = confSet->Get("SignalInjectionValue");
    if( param != "" ){
        fFitter->fSignalInjectionValue = Common::convertStoNum<float>(param);
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fFitter->fLimitParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fFitter->fLimitParamValue = Common::convertStoNum<float>(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fFitter->fLimitOutputPrefixName = param;
    }

    param = confSet->Get("UseAutoDiff");
    if( param != "" ){
        fFitter->fLimitUseAutoDiff = Common::StringToBoolean(param);
    }

    param = confSet->Get("ConfidenceLevel");
    if( param != "" ){
        double conf = Common::convertStoNum<double>(param);
        if (conf <= 0. || conf >= 1.){
            LOG(WARNING) << "Confidence level is <= 0 or >=1. Setting to default 0.95\n";
            conf = 0.95;
        }
        fFitter->fLimitsConfidence = conf;
    }

    param = confSet->Get("SplusBToys");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "SplusBToys is available only when TOYS are used\n";
        } else {
            fFitter->fLimitToysStepsSplusB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("BonlyToys");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "BonlyToys is available only when TOYS are used\n";
        } else {
            fFitter->fLimitToysStepsB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("ScanSteps");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ScanSteps is available only when TOYS are used\n";
        } else {
            fFitter->fLimitToysScanSteps = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("ScanMin");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ScanMin is available only when TOYS are used\n";
        } else {
            fFitter->fLimitToysScanMin = Common::convertStoNum<float>(param);
        }
    }

    param = confSet->Get("ScanMax");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ScanMax is available only when TOYS are used\n";
        } else {
            fFitter->fLimitToysScanMax = Common::convertStoNum<float>(param);
        }
    }

    param = confSet->Get("ToysSeed");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ToysSeed is available only when TOYS are used\n";
        } else {
            fFitter->fLimitToysSeed = Common::convertStoNum<float>(param);
        }
    }


    param = confSet->Get("LimitPlot");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "LimitPlot is available only when TOYS are used\n";
        } else {
            fFitter->fLimitPlot = Common::StringToBoolean(param);
        }
    }

    param = confSet->Get("LimitFile");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "LimitFile is available only when TOYS are used\n";
        } else {
            fFitter->fLimitFile = Common::StringToBoolean(param);
        }
    }

    param = confSet->Get("FitStrategy");
    if( param != "" ) {
        if (fFitter->fLimitType != TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "FitStrategy is available only when ASYMPTOTIC limits are used\n";
        } else {
            fFitter->fLimitFitStrategy = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("ToysUsexRooFit");
    if( param != "" ) {
        if (fFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ToysUsexRooFit is available only when TOYS are used\n";
        } else {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "FALSE") {
                fFitter->fLimitToysUsexRooFit = TRExFit::ToysUsexRooFit::FALSE;
            } else if (param == "TRUE") {
                fFitter->fLimitToysUsexRooFit = TRExFit::ToysUsexRooFit::TRUE;
            } else if (param == "AUTO") {
                fFitter->fLimitToysUsexRooFit = TRExFit::ToysUsexRooFit::AUTO;
            } else {
                LOG(ERROR) << "Unknown ToysUsexRooFit argument : " << confSet->Get("ToysUsexRooFit") << "\n";
                ++sc;
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadSignificanceOptions() {
    int sc(0);

    const ConfigSet* confSet = fParser->GetConfigSet("Significance");
    if (confSet == nullptr){
        LOG(DEBUG) << "You do not have Significance option in the config. It is ok, we just want to let you know.\n";
        return 0; // it is ok to not have Fit set up
    }

    // Set POI to be used for sig
    std::string param = confSet->Get("POI");
    if( param != "" ){
        fFitter->fPOIforSig = Common::CheckName(param);
        // add this to the list of POIs
        fFitter->AddPOI(fFitter->fPOIforSig);
    }

    // Set LimitBlind
    param = confSet->Get("SignificanceBlind");
    if( param != "" ){
        fFitter->fSignificanceIsBlind = Common::StringToBoolean(param);
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if( param != "" ){
        fFitter->fSignificancePOIAsimov = atof(param.c_str());
        fFitter->fSignificanceDoInjection = true;
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fFitter->fSignificanceParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fFitter->fSignificanceParamValue = Common::convertStoNum<float>(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fFitter->fSignificanceOutputPrefixName = param;
    }

    param = confSet->Get("UseAutoDiff");
    if( param != "" ){
        fFitter->fSignificanceUseAutoDiff = Common::StringToBoolean(param);
    }

    // Set SignificanceType
    param = confSet->Get("SignificanceType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "ASYMPTOTIC" ){
            fFitter->SetSignificanceType(TRExFit::SignificanceType::ASYMPTOTIC);
        }
        else if( param == "TOYS" ){
            fFitter->SetSignificanceType(TRExFit::SignificanceType::TOYS);
        }
        else{
            LOG(ERROR) << "Unknown SignificanceType argument : " << confSet->Get("SignificanceType") << "\n";
            ++sc;
        }
    }

    param = confSet->Get("SplusBToys");
    if(param != "") {
        if (fFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "SplusBToys is available only when TOYS are used\n";
        } else {
            fFitter->fSignificanceToysStepsSplusB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("BonlyToys");
    if( param != "" ) {
        if (fFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "BonlyToys is available only when TOYS are used\n";
        } else {
            fFitter->fSignificanceToysStepsB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("ToysSeed");
    if( param != "" ) {
        if (fFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "ToysSeed is available only when TOYS are used\n";
        } else {
            fFitter->fSignificanceToysSeed = Common::convertStoNum<float>(param);
        }
    }

    param = confSet->Get("SignificancePlot");
    if( param != "" ) {
        if (fFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "SignificancePlot is available only when TOYS are used\n";
        } else {
            fFitter->fSignificancePlot = Common::StringToBoolean(param);
        }
    }

    param = confSet->Get("ToysUsexRooFit");
    if( param != "") {
        if (fFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "ToysUsexRooFit is available only when TOYS are used\n";
        } else {
            fFitter->fSignificanceToysUsexRooFit = Common::StringToBoolean(param);
        }
    }

    param = confSet->Get("FitStrategy");
    if( param != "" ) {
        if (fFitter->fSignificanceType != TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "FitStrategy is available only when ASYMPTOTIC limits are used\n";
        } else {
            fFitter->fSignificanceFitStrategy = Common::convertStoNum<int>(param);
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadRegionOptions(const std::string& opt) {

    int sc(0);

    fAvailableRegions = GetAvailableRegions();

    // Check of the regions from commands like exist
    if (fOnlyRegions.size() > 0){
        if (!CheckPresence(fOnlyRegions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "You set regions that do not exist in your command line options\n";
            } else {
                LOG(ERROR) << "You set regions that do not exist in your command line options\n";
                ++sc;
            }
        }
    }

    const auto* region_confsets = fParser->GetConfigSets("Region");
    if(!region_confsets) {
        LOG(ERROR) << "No Region sections found in the configuration file.\n";
        return ++sc;
    }

    fHasAtLeastOneValidRegion = false;
    for (const auto& confSetPtr : *region_confsets) {
        const ConfigSet *confSet = confSetPtr.get();

        if(fOnlyRegions.size()>0 && Common::FindInStringVector(fOnlyRegions,Common::RemoveQuotes(confSet->GetValue()))<0){
            fFitter->SkippedRegion(Common::CheckName(confSet->GetValue()));
            continue;
        }
        else if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,Common::RemoveQuotes(confSet->GetValue()))>=0) {
            fFitter->SkippedRegion(Common::CheckName(confSet->GetValue()));
            continue;
        };


        fRegNames.emplace_back( Common::CheckName(confSet->GetValue()) );
        fRegions.emplace_back( Common::CheckName(confSet->GetValue()) );

        Region *reg;
        reg = fFitter->NewRegion(Common::CheckName(confSet->GetValue()));
        reg->fJobName = fFitter->fName;
        reg->fJobBinningsPath = fFitter->fBinningsPath;
        reg->fGetChi2 = fFitter->fGetChi2;
        reg->SetVariableTitle(Common::RemoveQuotes(confSet->Get("VariableTitle")));
        reg->SetLabel(Common::RemoveQuotes(confSet->Get("Label")),Common::RemoveQuotes(confSet->Get("ShortLabel")));

        // Set Type
        std::string param = confSet->Get("Type");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param=="CONTROL" )     reg -> SetRegionType(Region::CONTROL);
            else if( param=="VALIDATION" )  reg -> SetRegionType(Region::VALIDATION);
            else if( param=="SIGNAL" )      reg -> SetRegionType(Region::SIGNAL);
            else {
                LOG(ERROR) << "You specified 'Type' option in region but did not provide valid parameter. Please check this!\n";
                ++sc;
            }
        }
        if(reg -> fRegionType != Region::VALIDATION){
            reg->fUseGammaPulls = fFitter->fUseGammaPulls;
            fHasAtLeastOneValidRegion = true;
        }

        // Set PlotSubdir
        param = confSet->Get("PlotSubdir");
        if( param != "" ) reg->fPlotSubdir = Common::RemoveQuotes(param);

        // Set axisTitle
        param = confSet->Get("YaxisTitle");
        if( param != "") reg->fYTitle = Common::RemoveQuotes(param);

        // Set YmaxScale
        param = confSet->Get("YmaxScale");
        if(param!="") reg->fYmaxScale = atof(param.c_str());

        // Set Ymin
        param = confSet->Get("Ymin");
        if(param!="") reg->fYmin = atof(param.c_str());

        // Set Ymax
        param = confSet->Get("Ymax");
        if(param!="") reg->fYmax = atof(param.c_str());

        // Set RatioYmin
        param = confSet->Get("RatioYmin");
        if(param!="") {
            reg->fRatioYmin = atof(param.c_str());
            reg->fRatioYminPostFit = reg->fRatioYmin;
        }

        // Set RatioYmax
        param = confSet->Get("RatioYmax");
        if(param!="") {
            reg->fRatioYmax = atof(param.c_str());
            reg->fRatioYmaxPostFit = reg->fRatioYmax;
        }

        // Set RatioYminPostFit
        param = confSet->Get("RatioYminPostFit");
        if(param!="") reg->fRatioYminPostFit = atof(param.c_str());

        // Set RatioYmaxPostFit
        param = confSet->Get("RatioYmaxPostFit");
        if(param!="") reg->fRatioYmaxPostFit = atof(param.c_str());

        // Set TexLabel
        param = confSet->Get("TexLabel");
        if( param != "") reg->fTexLabel = Common::RemoveQuotes(param);

        // Set LumiLabel
        param = confSet->Get("LumiLabel");
        if( param != "") reg->fLumiLabel = Common::RemoveQuotes(param);
        else reg->fLumiLabel = fFitter->fLumiLabel;

        // Set CmeLabel
        param = confSet->Get("CmeLabel");
        if( param != "") reg->fCmeLabel = Common::RemoveQuotes(param);
        else reg->fCmeLabel = fFitter->fCmeLabel;

        // Set LogScale
        param = confSet->Get("LogScale");
        if( param != "" ){
            reg->fLogScale = Common::StringToBoolean(param);
        }

        param = confSet->Get("LogScaleX");
        if( param != "" ){
            reg->fLogScaleX = Common::StringToBoolean(param);
        }

        // Set Group
        param = confSet->Get("Group");
        if( param != "") reg->fGroup = Common::RemoveQuotes(param);

        // Paths for the unfolding code
        param = confSet->Get("ResponseMatrixFile");
        if (param != "") {
            reg->fResponseMatrixFiles.clear();
            reg->fResponseMatrixFiles.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixFiles");
        if (param != "") {
            reg->fResponseMatrixFiles = Common::Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("ResponseMatrixName");
        if (param != "") {
            reg->fResponseMatrixNames.clear();
            reg->fResponseMatrixNames.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixNames");
        if (param != "") {
            reg->fResponseMatrixNames = Common::Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("ResponseMatrixPath");
        if (param != "") {
            reg->fResponseMatrixPaths.clear();
            reg->fResponseMatrixPaths.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixPaths");
        if (param != "") {
            reg->fResponseMatrixPaths = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixFileSuff");
        if (param != "") {
            reg->fResponseMatrixFileSuffs.clear();
            reg->fResponseMatrixFileSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixFileSuffs");
        if (param != "") {
            reg->fResponseMatrixFileSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixNameSuff");
        if (param != "") {
            reg->fResponseMatrixNameSuffs.clear();
            reg->fResponseMatrixNameSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixNameSuffs");
        if (param != "") {
            reg->fResponseMatrixNameSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixPathSuff");
        if (param != "") {
            reg->fResponseMatrixPathSuffs.clear();
            reg->fResponseMatrixPathSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixPathSuffs");
        if (param != "") {
            reg->fResponseMatrixPathSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("AcceptanceName");
        if (param != "") {
            reg->fAcceptanceNames.clear();
            reg->fAcceptanceNames.emplace_back(Common::RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceNames");
        if (param != "") {
            reg->fAcceptanceNames = Common::Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceFile");
        if (param != "") {
            reg->fAcceptanceFiles.clear();
            reg->fAcceptanceFiles.emplace_back(Common::RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceFiles");
        if (param != "") {
            reg->fAcceptanceFiles = Common::Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        // Paths for the unfolding code
        param = confSet->Get("AcceptancePath");
        if (param != "") {
            reg->fAcceptancePaths.clear();
            reg->fAcceptancePaths.emplace_back(Common::RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptancePaths");
        if (param != "") {
            reg->fAcceptancePaths = Common::Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceFileSuff");
        if (param != "") {
            reg->fAcceptanceFileSuffs.clear();
            reg->fAcceptanceFileSuffs.emplace_back(Common::RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceFileSuffs");
        if (param != "") {
            reg->fAcceptanceFileSuffs = Common::Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceNameSuff");
        if (param != "") {
            reg->fAcceptanceNameSuffs.clear();
            reg->fAcceptanceNameSuffs.emplace_back(Common::RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptanceNameSuffs");
        if (param != "") {
            reg->fAcceptanceNameSuffs = Common::Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptancePathSuff");
        if (param != "") {
            reg->fAcceptancePathSuffs.clear();
            reg->fAcceptancePathSuffs.emplace_back(Common::RemoveQuotes(param));
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("AcceptancePathSuffs");
        if (param != "") {
            reg->fAcceptancePathSuffs = Common::Vectorize(param, ',');
            reg->fHasAcceptance = true;
        }

        param = confSet->Get("SelectionEffName");
        if (param != "") {
            reg->fSelectionEffNames.clear();
            reg->fSelectionEffNames.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNames");
        if (param != "") {
            reg->fSelectionEffNames = Common::Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("SelectionEffPath");
        if (param != "") {
            reg->fSelectionEffPaths.clear();
            reg->fSelectionEffPaths.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPaths");
        if (param != "") {
            reg->fSelectionEffPaths = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFile");
        if (param != "") {
            reg->fSelectionEffFiles.clear();
            reg->fSelectionEffFiles.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffFiles");
        if (param != "") {
            reg->fSelectionEffFiles = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFileSuff");
        if (param != "") {
            reg->fSelectionEffFileSuffs.clear();
            reg->fSelectionEffFileSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffFileSuffs");
        if (param != "") {
            reg->fSelectionEffFileSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffNameSuff");
        if (param != "") {
            reg->fSelectionEffNameSuffs.clear();
            reg->fSelectionEffNameSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNameSuffs");
        if (param != "") {
            reg->fSelectionEffNameSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffPathSuff");
        if (param != "") {
            reg->fSelectionEffPathSuffs.clear();
            reg->fSelectionEffPathSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPathSuffs");
        if (param != "") {
            reg->fSelectionEffPathSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationName");
        if (param != "") {
            reg->fMigrationNames.clear();
            reg->fMigrationNames.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNames");
        if (param != "") {
            reg->fMigrationNames = Common::Vectorize(param, ',');
        }

        // Paths for the unfolding code
        param = confSet->Get("MigrationPath");
        if (param != "") {
            reg->fMigrationPaths.clear();
            reg->fMigrationPaths.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPaths");
        if (param != "") {
            reg->fMigrationPaths = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFile");
        if (param != "") {
            reg->fMigrationFiles.clear();
            reg->fMigrationFiles.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationFiles");
        if (param != "") {
            reg->fMigrationFiles = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFileSuff");
        if (param != "") {
            reg->fMigrationFileSuffs.clear();
            reg->fMigrationFileSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationFileSuffs");
        if (param != "") {
            reg->fMigrationFileSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationNameSuff");
        if (param != "") {
            reg->fMigrationNameSuffs.clear();
            reg->fMigrationNameSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNameSuffs");
        if (param != "") {
            reg->fMigrationNameSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationPathSuff");
        if (param != "") {
            reg->fMigrationPathSuffs.clear();
            reg->fMigrationPathSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPathSuffs");
        if (param != "") {
            reg->fMigrationPathSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("NumberOfRecoBins");
        if (param != "") {
            const int bins = Common::convertStoNum<int>(param);
            if (bins < 2 || bins > 100) {
                LOG(WARNING) << "Number of reco bins is < 2 or > 100. This does not seem correct\n";
            }
            reg->fNumberUnfoldingRecoBins = bins;
        }

        param = confSet->Get("NormalizeMigrationMatrix");
        if (param != "") {
            reg->fNormalizeMigrationMatrix = Common::StringToBoolean(param);
        }

        // Setting based on input type
        if (fFitter->fInputType == 0){
            sc+= SetRegionHIST(reg, confSet);
        } else if (fFitter->fInputType == 1){
            sc+= SetRegionNTUP(reg, confSet) != 0;
        } else {
            LOG(ERROR) << "Unknown input type: " << fFitter->fInputType << "\n";
            ++sc;
        }

        // Set Rebin
        param = confSet->Get("Rebin");
        if(param != "") {
            reg->Rebin(atoi(param.c_str()));
        }

        // Set post-process rebin ("Rebinning")
        param = confSet->Get("Rebinning");
        if(param != ""){
            std::vector < std::string > vec_bins = Common::Vectorize(param, ',');
            if (vec_bins.size() == 0){
                LOG(ERROR) << "You specified `Rebinning` option, but you did not provide any reasonable option. Check this!\n";
                ++sc;
            }
            // eventually add auto-binning
            const unsigned int nBounds = vec_bins.size();
            std::vector<double> bins(nBounds);
            for (unsigned int iBound = 0; iBound < nBounds; ++iBound){
                bins[iBound] = atof(vec_bins[iBound].c_str());
            }
            reg -> SetRebinning(nBounds-1,bins);
        }

        param = confSet->Get("TransfoErr");
        if(param != "")
        {
          reg->fTransfoErr=std::atof(param.c_str());
        }

        // Set Binning
        param = confSet->Get("Binning");
        if(param != "" && param !="-"){

            std::vector < std::string > vec_bins = Common::Vectorize(param, ',');
            if (vec_bins.size() == 0){
                LOG(ERROR) << "You specified `Binning` option, but you did not provide any reasonable option. Check this!\n";
                return sc;
            }
            if(vec_bins[0]=="AutoBin"){
                if (vec_bins.size() < 2){
                    LOG(ERROR) << "You specified `Binning` option with Autobin, but you did not provide any reasonable option. Check this!\n";
                    return sc;
                }
                reg -> fBinTransfo = vec_bins[1];
                if(vec_bins[1]=="TransfoD"){
                    if (vec_bins.size() < 4){
                        LOG(ERROR) << "You specified `Binning` option with TransfoD, but you did not provide any reasonable option. Check this!\n";
                        return sc;
                    }
                    reg -> fTransfoDzSig=Common::convertStoNum<double>(vec_bins[2]);
                    reg -> fTransfoDzBkg=Common::convertStoNum<double>(vec_bins[3]);
                    if(vec_bins.size()>4){
                        for(const std::string& ibkg : vec_bins){
                            reg -> fAutoBinBkgsInSig.emplace_back(ibkg);
                        }
                    }
                }
                else if(vec_bins[1]=="TransfoF"){
                    if (vec_bins.size() < 4){
                        LOG(ERROR) << "You specified `Binning` option with TransfoF, but you did not provide any reasonable option. Check this!\n";
                        return sc;
                    }
                    if ( confSet->Get("TransfoErr") == "" )
                      reg -> fTransfoErr = 0.10;

                    reg -> fTransfoFzSig=Common::convertStoNum<double>(vec_bins[2]);
                    reg -> fTransfoFzBkg=Common::convertStoNum<double>(vec_bins[3]);
                    if(vec_bins.size()>4){
                        for(const std::string& ibkg : vec_bins){
                            reg -> fAutoBinBkgsInSig.emplace_back(ibkg);
                        }
                    }
                }
                else if(vec_bins[1]=="TransfoJ"){
                    if(vec_bins.size() > 2) reg -> fTransfoJpar1=Common::convertStoNum<double>(vec_bins[2]);
                    else reg -> fTransfoJpar1 = 5.;
                    if(vec_bins.size() > 3) reg -> fTransfoJpar2=Common::convertStoNum<double>(vec_bins[3]);
                    else reg -> fTransfoJpar2 = 1.;
                    if(vec_bins.size() > 4) reg -> fTransfoJpar3=Common::convertStoNum<double>(vec_bins[4]);
                    else reg -> fTransfoJpar3 = 5.;
                    if(vec_bins.size()>5){
                        for(const std::string& ibkg : vec_bins){
                            reg -> fAutoBinBkgsInSig.emplace_back(ibkg);
                        }
                    }
                }
                else{
                    LOG(ERROR) << "Unknown transformation: " << vec_bins[1] << ", try again\n";
                    ++sc;
                }
            }
            else{
                const unsigned int nBounds = vec_bins.size();
                std::vector<double> bins(nBounds);
                for (unsigned int iBound = 0; iBound < nBounds; ++iBound){
                    bins[iBound] = atof(vec_bins[iBound].c_str());
                }
                reg -> SetBinning(nBounds-1,bins);
            }
        }

        // Set BinWidth
        param = confSet->Get("BinWidth");
        if(param != "") reg->fBinWidth = atof(param.c_str());

        // Set DataType
        param = confSet->Get("DataType");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param=="DATA" )     reg -> SetRegionDataType(Region::REALDATA);
            else if( param=="ASIMOV" )  reg -> SetRegionDataType(Region::ASIMOVDATA);
            else{
                LOG(ERROR) << "DataType is not recognised: " << param << "\n";
                ++sc;
            }
        }

        // Set SkipSmoothing
        param = confSet->Get("SkipSmoothing");
        if( param != "" ){
            reg->fSkipSmoothing = Common::StringToBoolean(param);
        }

        // Set AutomaticDropBins
        param = confSet->Get("AutomaticDropBins");
        if (param != "") {
            bool isOK(true);
            if (reg->fDropBins.size() > 0) {
                LOG(WARNING) << "You specified 'AutomaticDropBins' option but you previously set DropBins for region " << reg->fName << " ignoring the automatic option\n";
                isOK = false;
            }
            if (isOK) {
                reg->SetAutomaticDropBins(Common::StringToBoolean(param));
            }
        }

        // Set DropBins
        param = confSet->Get("DropBins");
        if( param != "" ){
            if (reg->GetAutomaticDropBins()) {
                LOG(WARNING) << "You specified set `AutomaticDropBins` to TRUE, but using DropBins will disable it for region " << reg->fName << "!\n";
            }
            reg->SetAutomaticDropBins(false);
            reg->fDropBins.clear();
            const std::vector<std::string>& droppedBins = Common::Vectorize( param,',' );
            for(const std::string& binNum : droppedBins){
                reg->fDropBins.emplace_back(atoi(binNum.c_str()));
            }
        }

        // Set BlindBins
        param = confSet->Get("BlindBins");
        if( param != "" ){
            reg->fManualBlindBins.clear();
            const std::vector<std::string>& blindedBins = Common::Vectorize( param,',' );
            for(const std::string& binNum : blindedBins){
                reg->fManualBlindBins.emplace_back(atoi(binNum.c_str()));
            }
        }

        // Set BinLabels
        param = confSet->Get("BinLabels");
        if( param != "" ){
            std::vector<std::string> vec_string = Common::Vectorize( param,',' );
            reg->fBinLabels = vec_string;
        }

        // Set BinDividers
        param = confSet->Get("BinDividers");
        if( param != "" ){
            std::vector<std::string> vec_string = Common::Vectorize( param,',' );
            for (const std::string& s : vec_string) {
                reg->fBinDividers.emplace_back(atoi(s.c_str()));
            }
        }

        // Set RangeLabels
        param = confSet->Get("RangeLabels");
        if( param != "" ){
            std::vector<std::string> vec_string = Common::Vectorize( param,',' );
            if ( reg->fBinDividers.size()==0 ) {
                LOG(WARNING) << "'RangeLabels' cannot be used if 'BinDividers' not defined. Ignoring.\n";
            }
            else if ( vec_string.size() != reg->fBinDividers.size()+1) {
                LOG(WARNING) << "'RangeLabels' must have size 1 larger than 'BinDividers'. Ignoring.\n";
            }
            else {
                reg->fRangeLabels = vec_string;
            }
        }

        // Set BinTickMarks
        param = confSet->Get("BinTickMarks");
        if( param != "" ){
            const std::vector<std::string> s = Common::Vectorize( param,',' );
            std::vector<int> bins;
            for(const std::string& is : s){
                bins.emplace_back(std::atoi(is.c_str()));
            }
            if (!std::is_sorted(bins.begin(),bins.end())) {
                LOG(WARNING) << "Bin numbers specified by 'BinTickMarks' is not sorted. Ignoring.\n";
            }
            else reg->fBinTickMarks = bins;
        }

        // Set BinTickMarksLabels
        param = confSet->Get("BinTickMarksLabels");
        if( param != "" ){
            const std::vector<std::string> vec_string = Common::Vectorize( param,',' );
            if (vec_string.size() != reg->fBinTickMarks.size()) {
                LOG(WARNING) << "'BinTickMarksLabels' doesn't make sense if it doesn't have the same size as 'BinTickMarks'. Ignoring.\n";
            }
            else {
                reg->fBinTickMarksLabels = vec_string;
            }
        }

        // Set XaxisRange
        param = confSet->Get("XaxisRange");
        if( param != "" ){
            std::vector<std::string> vec_string = Common::Vectorize( param,',' );
            if (vec_string.size() != 2){
                LOG(ERROR) << "Setting 'XaxisRange' needs exactly two parameters (floats).\n";
                ++sc;
            } else {
                const double min = Common::convertStoNum<double>(vec_string.at(0));
                const double max = Common::convertStoNum<double>(vec_string.at(1));

                std::vector<double> range{};
                if (min < max){
                    range.emplace_back(min);
                    range.emplace_back(max);
                } else {
                    LOG(WARNING) << "Setting 'XaxisRange' needs the first parameter to be smaller than the second parameter. Ignoring\n";
                }
                reg->fXaxisRange = range;
            }
        }

        // Inter-region smoothing
        param = confSet->Get("IsBinOfRegion");
        if( param != "" ){
            std::vector<std::string> vec_string = Common::Vectorize( param,':' );
            if (vec_string.size() != 2){
                LOG(WARNING) << "Setting 'IsBinOfRegion' needs exactly two parameters (in the form string:int). Ignoring.\n";
            }
            else{
                reg->fIsBinOfRegion[vec_string.at(0)] = Common::convertStoNum<float>(vec_string.at(1));
            }
        }

        // Inter-region smoothing
        param = confSet->Get("DrawDataMinusBkg");
        if( param != "" ){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "NONE") {
                reg->fDrawDataMinusBkg = false;
            } else if (param == "ERRONSIGNAL") {
                reg->fDrawDataMinusBkg = true;
                reg->fDrawDataMinusBkgOnSignal = true;
            } else if (param == "ERRONBKG") {
                reg->fDrawDataMinusBkg = true;
                reg->fDrawDataMinusBkgOnSignal = false;
            }
        }

        // Random scaling of the NFs for pre/post-fit plots
        param = confSet->Get("RandomNFscale");
        if (param != "") {
            const std::vector<std::string> nfs = Common::Vectorize(param, ',');
            for (const auto& inf : nfs) {
                const std::vector<std::string> oneNF = Common::Vectorize(inf, ':');
                if (oneNF.size() != 3) {
                    LOG(ERROR) << "RandomNFscale option has wrong format\n";
                    ++sc;
                } else {
                    const double min = Common::convertStoNum<double>(oneNF.at(1));
                    const double max = Common::convertStoNum<double>(oneNF.at(2));

                    if (min > max) {
                        LOG(ERROR) << "RandomNFscale option for " << oneNF.at(0) << " has min > max\n";
                        ++sc;
                    } else {
                        reg->m_randomNFmap.insert({oneNF.at(0),std::make_pair(min, max)});
                    }
                }
            }
        }

        // Set PruningPlotSamples
        param = confSet->Get("PruningPlotSamples");
        if( param != "" ){
            reg->fPruningPlotSamples = Common::Vectorize(param,',');
        }

    }

    if (!fHasAtLeastOneValidRegion && Common::OptionRunsFit(opt)){
        LOG(ERROR) << "You need to provide at least one region that is not Validation otherwise the fit will crash.\n";
        ++sc;
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::SetRegionHIST(Region* reg, const ConfigSet *confSet){
    int sc(0);

    // Set HistoFile
    std::string param = confSet->Get("HistoFile");
    if(param!=""){
        reg->fHistoFiles.clear();
        reg->fHistoFiles.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set HistoFiles
    param = confSet->Get("HistoFiles");
    if(param!="") reg->fHistoFiles = Common::Vectorize( param,',' );

    // Set HistoName
    param = confSet->Get("HistoName");
    if(param!=""){
        reg->fHistoNames.clear();
        reg->fHistoNames.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set HistoNames
    param = confSet->Get("HistoNames");
    if(param!="") reg->fHistoNames = Common::Vectorize( param,',' );

    // Set HistoPath
    param = confSet->Get("HistoPath");
    if(param!=""){
        reg->fHistoPaths.clear();
        reg->fHistoPaths.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set HistoFiles
    param = confSet->Get("HistoPaths");
    if(param!="") reg->fHistoPaths = Common::Vectorize( param,',' );

    // Set HistoFileSuff
    param = confSet->Get("HistoFileSuff");
    if(param !=""){
        reg->fHistoFileSuffs.clear();
        reg->fHistoFileSuffs.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set HistoPathSuffs
    param = confSet->Get("HistoFileSuffs");
    if(param!=""){
        reg->fHistoFileSuffs.clear();
        std::vector<std::string> paths = Common::Vectorize( param,',' );
        for(const std::string& ipath : paths){
            reg->fHistoFileSuffs.emplace_back( Common::RemoveQuotes(ipath) );
        }
    }

    // Set HistoNameSuff
    param = confSet->Get("HistoNameSuff");
    if(param !=""){
        reg->fHistoNameSuffs.clear();
        reg->fHistoNameSuffs.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set HistoNameSuffs
    param = confSet->Get("HistoNameSuffs");
    if(param!=""){
        reg->fHistoNameSuffs.clear();
        std::vector<std::string> paths = Common::Vectorize( param,',' );
        for(const std::string& ipath : paths){
            reg->fHistoNameSuffs.emplace_back( Common::RemoveQuotes(ipath) );
        }
    }

    // Set HistoPathSuff
    param = confSet->Get("HistoPathSuff");
    if(param !=""){
        reg->fHistoPathSuffs.clear();
        reg->fHistoPathSuffs.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set HistoPathSuffs
    param = confSet->Get("HistoPathSuffs");
    if(param!=""){
        reg->fHistoPathSuffs.clear();
        std::vector<std::string> paths = Common::Vectorize( param,',' );
        for(const std::string& ipath : paths){
            reg->fHistoPathSuffs.emplace_back( Common::RemoveQuotes(ipath) );
        }
    }

    // Check for NTUP inputs
    if (ConfigHasNTUP(confSet)){
        LOG(WARNING) << "Found some NTUP settings in Region, while input option is HIST. Ignoring them.\n";
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::SetRegionNTUP(Region* reg, const ConfigSet *confSet){
    int sc(0);

    // Set Variable
    std::vector<std::string> variable = Common::Vectorize(confSet->Get("Variable"),',');
    if (variable.size() == 0){
        LOG(ERROR) << "Variable option is required but not present. Please check it!\n";
        ++sc;
    } else {
        // fix variable vector if special functions are used
        if(variable[0].find("Alt$")!=std::string::npos || variable[0].find("MaxIf$")!=std::string::npos ||variable[0].find("MinIf$")!=std::string::npos ){
            if (variable.size() > 1){
                do{
                    variable[0]+=","+variable[1];
                    variable.erase(variable.begin()+1);
                }
                while(variable[1].find("Alt$")!=std::string::npos);
                    variable[0]+=","+variable[1];
                    variable.erase(variable.begin()+1);
            } else {
                LOG(ERROR) << "Variable option has weird input. Please check it!\n";
                ++sc;
            }
        }
    }

    std::vector<std::string> corrVar  = Common::Vectorize(variable[0],'|');
    if(corrVar.size()==2){
        if (variable.size() < 4){
            LOG(ERROR) << "Setting \"Variable\". It seems you are trying to process a correlated variable but provided less than 4 parameters for each variable. Check this!\n";
            LOG(ERROR) << "The correct format is \"Expression1,nbins1, min1, max1|Expression2,nbins2, min2, max2\".\n";
            ++sc;
        } else if (variable.size() == 4){
            LOG(DEBUG) << "Have a correlation variable in reg " << fRegNames.back() << " : \n";
            LOG(DEBUG) << corrVar[0] << " and " << corrVar[1] << "\n";
            reg->SetVariable(  "corr_"+corrVar[0]+"_"+corrVar[1], std::atoi(variable[1].c_str()), std::atof(variable[2].c_str()), std::atof(variable[3].c_str()), corrVar[0].c_str(), corrVar[1].c_str() );
        } else {
            LOG(ERROR) << "Unexpected error\n";
            ++sc;
        }
    }
    else {
        if (variable.size() < 4){
            LOG(ERROR) << "Setting \"Variable\". You provided less than 4 parameters that are needed. Check this!\n";
            LOG(ERROR) << "The correct format is \"Expression,nbins,min,max\".\n";
            ++sc;
        } else {
            LOG(DEBUG) << "Have a usual variable in reg " << fRegNames.back() << " : \n";
            LOG(DEBUG) << variable[0] << " and size of corrVar=" << corrVar.size() << "\n";
            reg->SetVariable(  variable[0], atoi(variable[1].c_str()), atof(variable[2].c_str()), atof(variable[3].c_str()) );
        }
    }

    // Set VariableForSample
    std::string param = confSet->Get("VariableForSample");
    if( param != "" ){
        std::vector < std::string > temp_samplesAndVars = Common::Vectorize(param,',',false);
        for(const std::string& ivar : temp_samplesAndVars){
          std::vector < std::string > vars = Common::Vectorize(ivar,':');
            if(vars.size()==2){
                reg->SetAlternativeVariable(vars[1], vars[0]);
            }
        }
    }

    // Set Selection
    param = confSet->Get("Selection");
    if(param != "") reg->AddSelection( Common::RemoveQuotes(param) );

    // Set SelectionForSample
    param = confSet->Get("SelectionForSample");
    if( param != "" ){
        std::vector < std::string > temp_samplesAndSels = Common::Vectorize(param,',',false);
        for(const std::string& ivar : temp_samplesAndSels){
            std::vector < std::string > vars = Common::Vectorize(ivar,':');
            if(vars.size()==2){
                reg->SetAlternativeSelection(vars[1], vars[0]);
            }
        }
    }

    // Set MCweight
    param = confSet->Get("MCweight");
    if (param != "") reg->fMCweight = Common::RemoveQuotes(param); // this will override the global MCweight, if any

    // Set NtupleFile
    param = confSet->Get("NtupleFile");
    if(param!=""){
        reg->fNtupleFiles.clear();
        reg->fNtupleFiles.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set NtupleFiles
    param = confSet->Get("NtupleFiles");
    if(param!="") reg->fNtupleFiles = Common::Vectorize( param,',' );

    // Set NtupleFileSuff
    param = confSet->Get("NtupleFileSuff");
    if(param!="") {
        reg->fNtupleFileSuffs.clear();
        reg->fNtupleFileSuffs.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set NtupleFileSuffs
    param = confSet->Get("NtupleFileSuffs");
    if( param != "" ){
        std::vector<std::string> paths = Common::Vectorize( param,',' );
        reg->fNtupleFileSuffs = paths;
    }

    // Set NtupleName
    param = confSet->Get("NtupleName");
    if(param!="") {
        Common::CheckTreeName(param);
        reg->fNtupleNames.clear();
        reg->fNtupleNames.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set NtupleNames
    param = confSet->Get("NtupleNames");
    if (param != "") {
        Common::CheckTreeName(param);
        reg->fNtupleNames = Common::Vectorize( param,',' );
    }

    // Set NtupleNameSuff
    param = confSet->Get("NtupleNameSuff");
    if(param!="") {
        Common::CheckTreeName(param);
        reg->fNtupleNameSuffs.clear();
        reg->fNtupleNameSuffs.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set NtupleNameSuffs
    param = confSet->Get("NtupleNameSuffs");
    if( param != "" ){
        Common::CheckTreeName(param);
        std::vector<std::string> paths = Common::Vectorize( param,',' );
        reg->fNtupleNameSuffs = paths;
    }

    // Set NtuplePath
    param = confSet->Get("NtuplePath");
    if(param!="") {
        reg->fNtuplePaths.clear();
        reg->fNtuplePaths.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set NtuplePaths
    param = confSet->Get("NtuplePaths");
    if(param!="") reg->fNtuplePaths = Common::Vectorize( param,',' );

    // Set NtuplePathSuff
    param = confSet->Get("NtuplePathSuff");
    if(param != "") {
        reg->fNtuplePathSuffs.clear();
        reg->fNtuplePathSuffs.emplace_back( Common::RemoveQuotes(param) );
    }
    // Set NtuplePathSuffs
    param = confSet->Get("NtuplePathSuffs");
    if( param != "" ){
        std::vector<std::string> paths = Common::Vectorize( param,',' );
        reg->fNtuplePathSuffs = paths;
    }
    // Set FriendPath
    param = confSet->Get("FriendPath");
    if(param!="") {
        reg->fFriendPaths.clear();
        reg->fFriendPaths.push_back( Common::RemoveQuotes(param) );
    }
    // Set FriendPaths
    param = confSet->Get("FriendPaths");
    if(param!="") reg->fFriendPaths = Common::Vectorize( param,',' );

    // Set FriendPathSuff
    param = confSet->Get("FriendPathSuff");
    if(param != "") {
        reg->fFriendPathSuffs.clear();
        reg->fFriendPathSuffs.push_back( Common::RemoveQuotes(param) );
    }
    // Set FriendPathSuffs
    param = confSet->Get("FriendPathSuffs");
    if( param != "" ){
        std::vector<std::string> paths = Common::Vectorize( param,',' );
        reg->fFriendPathSuffs = paths;
    }
    // Check for NTUP inputs
    if (ConfigHasHIST(confSet)){
        LOG(WARNING) << "Found some HIST settings in Region, while input option is NTUP. Ignoring them.\n";
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadSampleOptions() {
    int sc(0);

    fAvailableSamples = GetAvailableSamples();

    // Check of the Samples from commands like exist
    if (fOnlyRegions.size() > 0){
        if (!CheckPresence(fOnlySamples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "You set samples that do not exist in your command line options\n";
            } else {
                LOG(ERROR) << "You set samples that do not exist in your command line options\n";
                ++sc;
            }
        }
    }

    fHasAtLeastOneValidSample = false;

    std::vector< std::unique_ptr<ConfigSet> > EFTConfSets;
    std::map<std::string,std::string> newEFTSampleLists;

    const auto* sample_confsets = fParser->GetConfigSets("Sample");
    if(!sample_confsets){
        LOG(WARNING) << "No Sample config sets found!\n";
        return sc;
    }
    for (std::size_t isample = 0; isample < sample_confsets->size(); ++isample) {
        ConfigSet *confSet = const_cast<ConfigSet*>((*sample_confsets)[isample].get());

        std::shared_ptr<Sample> sample = nullptr;
        std::shared_ptr<NormFactor> nfactor = nullptr;
        std::shared_ptr<ShapeFactor> sfactor = nullptr;
        Sample::SampleType type;
        std::string param = "";

        if(fOnlySamples.size()>0 && Common::FindInStringVector(fOnlySamples,Common::RemoveQuotes(confSet->GetValue()))<0) continue;
        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,Common::RemoveQuotes(confSet->GetValue()))>=0) continue;
        type = Sample::SampleType::BACKGROUND;

        // Set Type
        param = confSet->Get("Type");
        if (param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "SIGNAL"){
                type = Sample::SampleType::SIGNAL;
                fNonGhostIsSet = true;
                fHasAtLeastOneValidSample = true;
            }
            else if(param == "DATA"){
                type = Sample::SampleType::DATA;
                fNonGhostIsSet = true;
                fFitter->fHasDataSetInConfig = true;
            }
            else if(param == "EFT"){
                //@TODO: Make sure non-SMReference samples are set first
                type = Sample::SampleType::EFT;
                fNonGhostIsSet = true;
                fEFTSamples.emplace_back(Common::RemoveQuotes(confSet->GetValue()));
            }
            else if(param == "GHOST"){
                if (fNonGhostIsSet){
                    LOG(ERROR) << "Please define GHOST samples first and then other samples\n";
                    ++sc;
                }
                type = Sample::SampleType::GHOST;
                fGhostSamples.emplace_back(Common::RemoveQuotes(confSet->GetValue()));
            }
            else if(param == "BACKGROUND"){
                type = Sample::SampleType::BACKGROUND;
                fNonGhostIsSet = true;
                fHasAtLeastOneValidSample = true;
            }
            else {
                LOG(WARNING) << "You specified 'Type' option in sample but did not provide valid parameter. Using default (BACKGROUND)\n";
            }
            if(fOnlySignals.size()>0){
                if(type==Sample::SampleType::SIGNAL){
                    bool skip = true;
                    for(const auto& s : fOnlySignals){
                        if(Common::CheckName(confSet->GetValue())==s) skip = false;
                    }
                    if(skip) continue;
                }
                else if(type==Sample::SampleType::GHOST){
                    for(const auto& s : fOnlySignals){
                        if(Common::CheckName(confSet->GetValue())==s) type = Sample::SampleType::SIGNAL;
                    }
                }
            }
        }

        sample = fFitter->NewSample(Common::CheckName(confSet->GetValue()),type);
        fSamples.emplace_back(Common::CheckName(confSet->GetValue()));

        // Set Title
        param = confSet->Get("Title");
        if (param != "") sample->SetTitle(Common::RemoveQuotes(param));

        // Set TexTitle
        param = confSet->Get("TexTitle");
        if(param!="") sample->fTexTitle = Common::RemoveQuotes(param);

        // Set Group
        param = confSet->Get("Group");
        if(param!="") sample->fGroup = Common::RemoveQuotes(param);

        // HIST input
        if (fFitter->fInputType == 0){
            // Set HistoFile
            param = confSet->Get("HistoFile");
            if(param!="") sample->fHistoFiles.emplace_back( Common::RemoveQuotes(param) );

            // Set HistoFiles
            param = confSet->Get("HistoFiles");
            if(param!="") sample->fHistoFiles = Common::Vectorize( param, ',' );

            // Set HistoFileSuff
            param = confSet->Get("HistoFileSuff");
            if(param!="") sample->fHistoFileSuffs.emplace_back( Common::RemoveQuotes(param) );

            // Set HistoFileSuffs
            param = confSet->Get("HistoFileSuffs");
            if(param!="") sample->fHistoFileSuffs = Common::Vectorize( param, ',' );

            // Set HistoName
            param = confSet->Get("HistoName");
            if(param!="") sample->fHistoNames.emplace_back( Common::RemoveQuotes(param) );

            // Set HistoNames
            param = confSet->Get("HistoNames");
            if(param!="") sample->fHistoNames = Common::Vectorize( param, ',' );

            // Set HistoNameSuff
            param = confSet->Get("HistoNameSuff");
            if(param!="") sample->fHistoNameSuffs.emplace_back( Common::RemoveQuotes(param) );

            // Set HistoNameSuffs
            param = confSet->Get("HistoNameSuffs");
            if(param!="") sample->fHistoNameSuffs = Common::Vectorize( param, ',' );

            // Set HistoPath
            param = confSet->Get("HistoPath");
            if(param!="") sample->fHistoPaths.emplace_back( Common::RemoveQuotes(param) );

            // Set HistoPaths
            param = confSet->Get("HistoPaths");
            if(param!="") sample->fHistoPaths = Common::Vectorize( param, ',' );

            // Set HistoPathSuff
            param = confSet->Get("HistoPathSuff");
            if(param!="") sample->fHistoPathSuffs.emplace_back( Common::RemoveQuotes(param) );

            // Set HistoPathSuffs
            param = confSet->Get("HistoPathSuffs");
            if(param!="") sample->fHistoPathSuffs = Common::Vectorize( param, ',' );

            if (ConfigHasNTUP(confSet)){
                LOG(WARNING) << "You provided some NTUP options but your input type is HIST. Options will be ignored\n";
            }
        } else if (fFitter->fInputType == 1){ // NTUP input
            // Set NtupleFile
            param = confSet->Get("NtupleFile");
            if(param!="") sample->fNtupleFiles.emplace_back( Common::RemoveQuotes(param) );

            // Set NtupleFiles
            param = confSet->Get("NtupleFiles");
            if(param!="") sample->fNtupleFiles = Common::Vectorize( param ,',' );

            param = confSet->Get("NtupleFileSuff");
            if(param!="") sample->fNtupleFileSuffs.emplace_back( Common::RemoveQuotes(param) );

            // Set NtupleFileSuffs
            param = confSet->Get("NtupleFileSuffs");
            if(param!="") sample->fNtupleFileSuffs = Common::Vectorize( param ,',' );

            // Set NtupleName
            param = confSet->Get("NtupleName");
            if(param!="") {
                Common::CheckTreeName(param);
                sample->fNtupleNames.emplace_back( Common::RemoveQuotes(param) );
            }

            // Set NtupleNames
            param = confSet->Get("NtupleNames");
            if(param!="") {
                Common::CheckTreeName(param);
                sample->fNtupleNames = Common::Vectorize( param ,',' );
            }

            // Set NtupleNameSuff
            param = confSet->Get("NtupleNameSuff");
            if(param!="") {
                Common::CheckTreeName(param);
                sample->fNtupleNameSuffs.emplace_back( Common::RemoveQuotes(param) );
            }

            // Set NtupleNameSuffs
            param = confSet->Get("NtupleNameSuffs");
            if( param != "" ) {
                Common::CheckTreeName(param);
                sample->fNtupleNameSuffs = Common::Vectorize( param,',' );
            }

            // Set NtuplePath
            param = confSet->Get("NtuplePath");
            if(param!="") sample->fNtuplePaths.emplace_back( Common::RemoveQuotes(param) );

            // Set NtuplePaths
            param = confSet->Get("NtuplePaths");
            if(param != "") sample->fNtuplePaths = Common::Vectorize( param ,',' );

            // Set NtuplePathSuff
            param = confSet->Get("NtuplePathSuff");
            if(param!="") sample->fNtuplePathSuffs.emplace_back( Common::RemoveQuotes(param) );

            // Set NtuplePathSuffs
            param = confSet->Get("NtuplePathSuffs");
            if(param != "") sample->fNtuplePathSuffs = Common::Vectorize( param ,',' );

            // Set FriendPath
            param = confSet->Get("FriendPath");
            if(param!="") sample->fFriendPaths.push_back( Common::RemoveQuotes(param) );

            // Set FriendPaths
            param = confSet->Get("FriendPaths");
            if(param != "") sample->fFriendPaths = Common::Vectorize( param ,',' );

            // Set FriendPathSuff
            param = confSet->Get("FriendPathSuff");
            if(param!="") sample->fFriendPathSuffs.push_back( Common::RemoveQuotes(param) );

            // Set FriendPathsSuffs
            param = confSet->Get("FriendPathSuffs");
            if(param != "") sample->fFriendPathSuffs = Common::Vectorize( param ,',' );

            // Set GridDID
            // indicates Sample is delivered by ServiceX
            param = confSet->Get("GridDID");
            if(param != "") sample->fGridDID = Common::RemoveQuotes(param);

            if (ConfigHasHIST(confSet)){
                LOG(WARNING) << "You provided some HIST options but your input type is NTUP. Options will be ignored\n";
            }
        } else {
            LOG(ERROR) << "No valid input type provided. Please check this!\n";
            ++sc;
        }

        // Set UseFriend
        param = confSet->Get("UseFriend");
        if(param !=""){
            sample->fUseFriend = Common::StringToBoolean(param);
        }

        // Set UseGaussianShapeSysConstraint
        param = confSet->Get("UseGaussianShapeSysConstraint");
        if(param != ""){
            sample->fUseGaussianShapeSysConstraint = Common::StringToBoolean(param);
        }

        // Set FillColor
        param = confSet->Get("FillColor");
        if(param != "") sample->SetFillColor(atoi(param.c_str()));

        // Set LineColor
        param = confSet->Get("LineColor");
        if(param != "") sample->SetLineColor(atoi(param.c_str()));

        // Set FillColor from RGB values if given
        param = confSet->Get("FillColorRGB");
        if (param != "") {
            const std::vector<std::string> rgb_strings = Common::Vectorize(param, ',');
            sample->fFillColor = Common::ColorFromRGB(rgb_strings);
        }

        // Set LineColor from RGB values if given
        param = confSet->Get("LineColorRGB");
        if (param != "") {
            const std::vector<std::string> rgb_strings = Common::Vectorize(param, ',');
            sample->fLineColor = Common::ColorFromRGB(rgb_strings);
        }

        // Set FillColor from HEX values if given
        std::string hexColor = confSet->Get("FillLineColorHEX");
        if (hexColor != "") {
            hexColor = Common::RemoveQuotes(hexColor);
            sample->SetFillColor(TColor::GetColor(hexColor.c_str()));
            sample->SetLineColor(TColor::GetColor(hexColor.c_str()));
        }

        // Set NormalizedByTheory
        param = confSet->Get("NormalizedByTheory");
        if(param != ""){
            sample->NormalizedByTheory(Common::StringToBoolean(param));
        }

        // Set MCweight and Selection
        if(fFitter->fInputType==1){
            param = confSet->Get("MCweight");
            if(param != "") sample->SetMCweight( Common::RemoveQuotes(param) );

            param = confSet->Get("Selection");
            if(param!="") sample->SetSelection( Common::RemoveQuotes(param) );
        }

        // to specify only certain regions
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        std::vector<std::string> regions = Common::Vectorize(regions_str,',');
        std::vector<std::string> exclude = Common::Vectorize(exclude_str,',');
        sample->fRegions.clear();

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
            } else {
                LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
                ++sc;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has regions to exclude set up that do not exist\n";
            } else {
                LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has regions to exclude set up that do not exist\n";
                ++sc;
            }
        }

        for(const auto& ireg : fFitter->fRegions) {
            const std::string regName = ireg->fName;
            if( (regions_str=="" || regions_str=="all" || Common::FindInStringVector(regions,regName)>=0)
                && Common::FindInStringVector(exclude,regName)<0 ){
                sample->fRegions.emplace_back(ireg->fName);
            }
        }

        // Set NormFactor
        param = confSet->Get("NormFactor");
        if(param!=""){
            // check if the normfactor is called just with the name or with full definition
            const unsigned int sz = Common::Vectorize(param,',').size();
            if (sz != 1 && sz != 4 && sz != 5){
                LOG(ERROR) << "No valid input for 'NormFactor' provided. Please check this!\n";
                ++sc;
            }
            if( sz > 1 ){
                bool isConst = false;
                if( Common::Vectorize(param,',').size()>4){
                    std::string tmp = Common::Vectorize(param,',')[4];
                    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
                    if (tmp == "TRUE"){
                        isConst = true;
                    }
                }
                if (sz > 3)
                    nfactor = sample->AddNormFactor(
                    Common::Vectorize(param,',')[0],
                    atof(Common::Vectorize(param,',')[1].c_str()),
                    atof(Common::Vectorize(param,',')[2].c_str()),
                    atof(Common::Vectorize(param,',')[3].c_str()),
                    isConst
                );
            }
            else{
                nfactor = sample->AddNormFactor( Common::Vectorize(param,',')[0] );
            }
            nfactor->fRegions = sample->fRegions;
            TRExFitter::SYSTMAP[nfactor->fName] = nfactor->fName;
            if( Common::FindInStringVector(fFitter->fNormFactorNames,nfactor->fName)<0 ){
                fFitter->fNormFactors.emplace_back( nfactor );
                fFitter->fNormFactorNames.emplace_back( nfactor->fName );
            }

            TRExFitter::NPMAP[nfactor->fName] = nfactor->fName;
        }

        // Set ShapeFactor
        param = confSet->Get("ShapeFactor");
        if(param!=""){
            // check if the normfactor is called just with the name or with full definition
            const unsigned int sz = Common::Vectorize(param,',').size();
            if (sz != 1 && sz != 4 && sz != 5){
                LOG(ERROR) << "No valid input for 'ShapeFactor' provided. Please check this!\n";
                ++sc;
            }
            if( sz > 1 ){
                bool isConst = false;
                if( Common::Vectorize(param,',').size()>4){
                    std::string tmp = Common::Vectorize(param,',')[4];
                    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
                    if (tmp == "TRUE"){
                        isConst = true;
                    }
                }
                if (sz > 3) {
                    sfactor = sample->AddShapeFactor(
                    Common::Vectorize(param,',')[0],
                    atof(Common::Vectorize(param,',')[1].c_str()),
                    atof(Common::Vectorize(param,',')[2].c_str()),
                    atof(Common::Vectorize(param,',')[3].c_str()),
                    isConst);
                    sfactor->fSamples = Common::ToVec(sample->fName);
                }
            }
            else{
                sfactor = sample->AddShapeFactor( Common::Vectorize(param,',')[0] );
                sfactor->fSamples = Common::ToVec(sample->fName);
            }
            sfactor->fRegions = sample->fRegions;
            if( Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName)<0 ){
                fFitter->fShapeFactors.emplace_back( sfactor );
                fFitter->fShapeFactorNames.emplace_back( sfactor->fName );
            }
        }

        // Set LumiScale
        param = confSet->Get("LumiScale");
        if(param!="") sample->fLumiScales.emplace_back( atof(param.c_str()) );

        // Set LumiScales
        param = confSet->Get("LumiScales");
        if(param!=""){
            std::vector<std::string> lumiScales_str = Common::Vectorize( param ,',' );
            for(const std::string& ilumiscale : lumiScales_str){
                sample->fLumiScales.emplace_back( atof(ilumiscale.c_str()) );
            }
        }

        // Set IgnoreSelection
        // to skip global & region selection for this sample
        param = confSet->Get("IgnoreSelection");
        if(param != ""){
            sample->fIgnoreSelection = Common::RemoveQuotes(param);
        }

        // Set IgnoreWeights
        // to skip global & region weights for this sample
        param = confSet->Get("IgnoreWeight");
        if(param!=""){
            sample->fIgnoreWeight = Common::RemoveQuotes(param);
        }

        // Set UseMCstat
        // to skip MC stat uncertainty for this sample
        param = confSet->Get("UseMCstat");
        if(param != ""){
            sample->fUseMCStat = Common::StringToBoolean(param);
        }

        // Set UseSystematics
        // to skip MC systematics for this sample
        param = confSet->Get("UseSystematics");
        // set it to false for ghost samples and data and true for other samples
        if(type == Sample::SampleType::GHOST || type == Sample::SampleType::DATA || type == Sample::SampleType::EFT) {
            sample->fUseSystematics = false;
        } else {
            sample->fUseSystematics = true;
        }
        if(param != ""){
            sample->fUseSystematics = Common::StringToBoolean(param);
        }

        // Set DivideBy
        param = confSet->Get("DivideBy");;
        if (param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), Common::RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for DivideBy that do not exist\n";
                } else {
                    LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for DivideBy that do not exist\n";
                    ++sc;
                }
            }
            sample->fDivideBy = Common::RemoveQuotes(param);
        }

        // Set MultiplyBy
        param = confSet->Get("MultiplyBy");
        if (param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), Common::RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for MultiplyBy that do not exist\n";
                } else {
                    LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for MultiplyBy that do not exist\n";
                    ++sc;
                }
            }
            sample->fMultiplyBy = Common::RemoveQuotes(param);
        }

        // Set SubtractSample
        param = confSet->Get("SubtractSample");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), Common::RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for SubtractSample that do not exist\n";
                } else {
                    LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for SubtractSample that do not exist\n";
                    ++sc;
                }
            }
            sample->fSubtractSamples.emplace_back( Common::RemoveQuotes(param) );
        }

        // Set SubtractSamples
        param = confSet->Get("SubtractSamples");
        if(param != ""){
            std::vector<std::string> tmp = Common::Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for SubtractSamples that do not exist\n";
                } else {
                    LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for SubtractSamples that do not exist\n";
                    ++sc;
                }
            }
            sample->fSubtractSamples = Common::Vectorize(param,',');
        }

        // Set AddSample
        param = confSet->Get("AddSample");
        if(param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), Common::RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for AddSample that do not exist\n";
                } else {
                    LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for AddSample that do not exist\n";
                    ++sc;
                }
            }
            sample->fAddSamples.emplace_back( Common::RemoveQuotes(param) );
        }

        // Set AddSamples
        param = confSet->Get("AddSamples");
        if(param!=""){
            std::vector<std::string> tmp = Common::Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for AddSamples that do not exist\n";
                } else {
                    LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for AddSamples that do not exist\n";
                    ++sc;
                }
            }
            sample->fAddSamples = Common::Vectorize(param,',');
        }

        // Set NormToSample
        param = confSet->Get("NormToSample");
        if(param != ""){
            if (std::find(fSamples.begin(), fSamples.end(), Common::RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for NormToSample that do not exist\n";
                } else {
                    LOG(ERROR) << "Sample: " << Common::CheckName(confSet->GetValue()) << " has samples set up for NormToSample that do not exist\n";
                    ++sc;
                }
            }
            sample->fNormToSample = Common::RemoveQuotes(param);
        }

        // Set Smooth
        // allow smoothing of nominal histogram?
        param = confSet->Get("Smooth");
        if(param!=""){
            sample->fSmooth = Common::StringToBoolean(param);
        }

        // Set AsimovReplacementFor
        // AsimovReplacementFor
        param = confSet->Get("AsimovReplacementFor");
        if( param!=""){
           std::vector<std::string> asimovreplacements = Common::Vectorize(param,',');
           for(const auto& asimovreplacementstring : asimovreplacements){
               if(Common::RemoveQuotes(asimovreplacementstring) != ""){
                   if(Common::Vectorize(asimovreplacementstring,'=').size() == 2){
                        std::vector<std::string> tmp = Common::Vectorize(asimovreplacementstring,'=');
                        std::pair<std::string,std::string> asimovreplacement = std::make_pair(tmp.at(0),tmp.at(1));
                        sample->fAsimovReplacementFor.push_back(asimovreplacement);
                   } else {
                       LOG(ERROR) << "You specified 'AsimovReplacementFor' option but did not provide 2 parameters, separated by an = sign. Please check this\n";
                       LOG(ERROR) << "example: AsimovReplacementFor: replacementname1=sampletoreplace1,replacementname2=sampletoreplace2\n";
                       ++sc;
                   }
               }
           }
        }

        // Set SeparateGammas
        // separate gammas
        param = confSet->Get("SeparateGammas");
        if(param != ""){
            if (!fFitter->fUseStatErr) {
                LOG(WARNING) << "MCstatThreshold set to NONE, skipping SeparateGammas\n";
            } else {
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if(param == "TRUE"){
                    sample->fSeparateGammas = true;
                    if(confSet->Get("UseMCstat") == "") sample->fUseMCStat = false; // remove the usual gammas for this sample (only if no UseMCstat is specified!!)
                } else sample->fSeparateGammas = false;
            }
        }

        // Scale up/down MC stat size
        param = confSet->Get("MCstatScale");
        if(param != ""){
            sample->fMCstatScale = atof(param.c_str());
        }

        // Set CorrelateGammasInRegions
        // in the form    CorrelateGammasInRegions: SR1:SR2,CR1:CR2:CR3
        param = confSet->Get("CorrelateGammasInRegions");
        if(param != ""){
            std::vector<std::string> sets = Common::Vectorize(param,',',false);
            for(const std::string& set : sets){
                std::vector<std::string> regions_corr = Common::Vectorize(set,':');
                LOG(DEBUG) << "Correlating gammas for this sample in regions " << set << "\n";
                sample->fCorrelateGammasInRegions.emplace_back(regions_corr);
            }
        }

        // Set CorrelateGammasWithSample
        param = confSet->Get("CorrelateGammasWithSample");
        if(param != ""){
            sample->fCorrelateGammasWithSample = Common::RemoveQuotes(param);
        }

        // Set SystFromSample
        param = confSet->Get("SystFromSample");
        if(param != ""){
            sample->fSystFromSample = Common::RemoveQuotes(param);
        }

        // Set Morphing
        param = confSet->Get("Morphing");
        if(param != ""){
            if(fFitter->fExcludeFromMorphing!=sample->fName){
                std::vector<std::string> morph_par = Common::Vectorize(param,',');
                if (morph_par.size() != 2){
                    LOG(ERROR) << "Morphing requires exactly 2 parameters, but " << morph_par.size() << " provided\n";
                    ++sc;
                }
                std::string name      = morph_par.at(0);
                double value = Common::convertStoNum<double>(morph_par.at(1));
                LOG(DEBUG) << "Morphing: Adding " << name << ", with value: " << value << "\n";
                if (!fFitter->MorphIsAlreadyPresent(name, value)) fFitter->AddTemplateWeight(name, value);
                // set proper normalization
                std::string morphName = "morph_"+name+"_"+Common::ReplaceString(std::to_string(value),"-","m");
                std::shared_ptr<NormFactor> nf = sample->AddNormFactor(morphName, 1, 0, 10, false);
                fFitter->fNormFactors.emplace_back( nf );
                fFitter->fNormFactorNames.emplace_back( nf->fName );
                sample->fIsMorph[name] = true;
                sample->fMorphValue[name] = value;
                if(Common::FindInStringVector(fFitter->fMorphParams,name)<0) fFitter->fMorphParams.emplace_back( name );
            }
        }

        // Set EFT Standard Model reference sample
        param = confSet->Get("EFTSMReference");
        if(param != ""){
            std::string tmpparam=param;
            std::transform(tmpparam.begin(), tmpparam.end(), tmpparam.begin(), ::toupper);
            if(tmpparam=="NONE") {
                sample->fEFTSMReference = "NONE";
            } else {
                sample->fEFTSMReference = Common::RemoveQuotes(param);
            }
        }

        // check that the EFT reference sample is listed after all eft samples
        if (type == Sample::SampleType::EFT && sample->fEFTSMReference != "NONE") {
            auto itr = std::find(fSamples.begin(), fSamples.end(), sample->fEFTSMReference);
            if (itr != fSamples.end()) {
                LOG(ERROR) << "The EFT reference sample is not defined after the EFT samples\n";
                LOG(ERROR) << "The current code implementation does not allow that\n";
                LOG(ERROR) << "Please, reorder the samples so that the reference sample appears after the EFT samples in the config\n";
                sc++;
                return sc;
            }
        }

        // Set EFT parameter name
        param = confSet->Get("EFTParam");
        if(param != ""){
            sample->fEFTParam = Common::RemoveQuotes(param);
        }

        // Set EFT parameter title
        param = confSet->Get("EFTTitle");
        if(param != ""){
            sample->fEFTTitle = Common::RemoveQuotes(param);
        }

        // Set EFT parameter value
        param = confSet->Get("EFTValue");
        if(param != ""){
            sample->fEFTValue = Common::RemoveQuotes(param);
            const std::vector<std::string> tmpPar = Common::Vectorize(sample->fEFTValue, ',');
            for (const auto& par : tmpPar) {
                const std::vector<std::string> tmp = Common::Vectorize(par, '=');
                if (tmp.size() != 2) {
                    LOG(ERROR) << "Wrong \"EFTValue\" setting for sample: " << sample->fName << "\n";
                    ++sc;
                } else {
                    auto itr = std::find(fEFTParams.begin(), fEFTParams.end(), tmp.at(0));
                    if (itr == fEFTParams.end()) {
                        fEFTParams.emplace_back(tmp.at(0));
                    }

                    auto itrMap = fEFTOperatorsPerSample.find(sample->fEFTSMReference);
                    if (itrMap == fEFTOperatorsPerSample.end()) {
                        fEFTOperatorsPerSample.insert(std::make_pair(sample->fEFTSMReference, std::vector<std::string>{tmp.at(0)}));
                    } else {
                        auto itrVec = std::find(itrMap->second.begin(), itrMap->second.end(), tmp.at(0));
                        if (itrVec == itrMap->second.end()) {
                            itrMap->second.emplace_back(tmp.at(0));
                        }
                    }
                }
            }
        }

        // Set if pruning should be used
        param = confSet->Get("NoPruning");
        if(param != ""){
            sample->fNoPruning = Common::StringToBoolean(param);
        }

        param = confSet->Get("HistoFolderName");
        if (param != "") {
            sample->fHistoFolderName = Common::RemoveQuotes(param);
        }

        param = confSet->Get("Template");
        if (param != "") {
            fHasTemplateMorphing = true;
            fFitter->fHasTemplateMorphing = true;
            sample->fAlternativeRefInTemplateMorphing = confSet->Get("AlternativeRefInTemplateMorphing");
            const auto initVec = Common::Vectorize(param, '@');
            bool hasSampleIdentifier(false);
            if (initVec.size() == 2) {
                hasSampleIdentifier = true;
            } else if (initVec.size() > 2) {
                LOG(ERROR) << "Invalid \"Template\" argument\n";
                ++sc;
                return sc;
            }

            if (hasSampleIdentifier) {
                sample->fTemplateTarget = initVec.at(0);
                auto itr = std::find(fTemplateTargets.begin(), fTemplateTargets.end(), sample->fTemplateTarget);
                if (itr == fTemplateTargets.end()) {
                    fTemplateTargets.emplace_back(sample->fTemplateTarget);
                }
            }

            const auto vec = Common::Vectorize(hasSampleIdentifier ? initVec.at(1) : param, ',');
            if (vec.size() > 2) {
                LOG(ERROR) << "Maximum 2D template morphing is supported\n";
                ++sc;
            }

            if (fTemplateMorphingDimensions < 0) {
                fTemplateMorphingDimensions = static_cast<int>(vec.size());
            } else {
                if (static_cast<int>(vec.size()) != fTemplateMorphingDimensions) {
                    LOG(ERROR) << "Inconsistent template dimensionality\n";
                    ++sc;
                }
            }

            for (const auto& i : vec) {
                const auto& tmp = Common::Vectorize(i, ':');
                if (tmp.size() != 2) {
                    LOG(ERROR) << "Wrong format of the \"Template\" config option for sample" << sample->fName << "\n";
                    return sc;
                }

                double value(0);
                try {
                    value = Common::convertStoNum<float>(tmp.at(1));
                } catch (...) {
                    LOG(ERROR) << "Cannot convert the second argument of \"Template\" to a number for sample " << sample->fName << "\n";
                    return sc;
                }
                sample->fTemplateMorphing.emplace_back(std::make_pair(tmp.at(0), value));
            }
        }

        if (fFitter->fEFTConfig.GetSplitSamplesPerBin()) {
            // EFT Sample Copying for SM References
            if(type==Sample::SampleType::EFT && sample->fEFTSMReference == "NONE"){
                LOG(DEBUG) << "Found EFT SM Reference sample " << sample->fName << "\n";
                for(const auto& ireg : fFitter->fRegions) {
                    const std::string regName = ireg->fName;
                    if( (regions_str!="" && regions_str!="all" && Common::FindInStringVector(regions,regName)<0) ) continue;
                    if( Common::FindInStringVector(exclude,regName)>=0 )continue;

                    // Creating new SIGNAL Samples based on bins in regions
                    for(int i=0; i < ireg->GetNbins(); ++i){

                        std::shared_ptr<Sample> EFTSM_s = std::make_shared<Sample>(*sample);
                        EFTSM_s->fName="SM_"+sample->fName+"_"+ireg->fName+"_bin"+i;
                        EFTSM_s->fUseSystematics = true;
                        EFTSM_s->fType=Sample::SampleType::SIGNAL;
                        EFTSM_s->fRegions=Common::ToVec(ireg->fName);
                        EFTSM_s->fEFTSMReference=sample->fName;
                        fSamples.emplace_back(EFTSM_s->fName);
                        fFitter->fSamples.emplace_back(EFTSM_s);
                        fAvailableSamples.emplace_back(EFTSM_s->fName);
                        LOG(DEBUG) << Form("   Created new Sample %s",EFTSM_s->fName.c_str()) << "\n";

                        // Make string list of new EFT samples
                        std::map<std::string,std::string>::iterator it;
                        it = newEFTSampleLists.find(sample->fName);
                        if (it != newEFTSampleLists.end()){
                            newEFTSampleLists[sample->fName]=newEFTSampleLists[sample->fName]+","+EFTSM_s->fName;
                        }else{
                            newEFTSampleLists[sample->fName]=EFTSM_s->fName;
                        }

                        EFTConfSets.emplace_back(std::make_unique<ConfigSet>());
                        EFTConfSets[EFTConfSets.size()-1]->Set("Sample",EFTSM_s->fName);
                        for(auto& iconf : confSet->fConfig){
                            if(iconf.fName=="UseSystematics"){iconf.fValue = "TRUE";}
                            else if(iconf.fName=="Type"){iconf.fValue = "SIGNAL";}
                            else if(iconf.fName=="Regions"){iconf.fValue=ireg->fName;}
                            else if(iconf.fName=="EFTSMReference"){iconf.fValue=sample->fName;}
                            EFTConfSets[EFTConfSets.size()-1]->SetConfig(iconf.fName,iconf.fValue,iconf.fFileName,iconf.fLineNumber);
                        }

                        // Adding NF for each EFT param
                        for(const auto& ismp : fFitter->fSamples) {
                            // We don't know what EFT params there are at Config reading level so just have to check each time and add them if they don't aready exist...
                            if(ismp->fType==Sample::SampleType::EFT && ismp->fEFTSMReference != "NONE"){
                                const std::string NFname=TString::Format("Expression_muEFT_%s_%s_bin%d",sample->fName.c_str(),ireg->fName.c_str(),i).Data();
                                if( Common::FindInStringVector(fFitter->fNormFactorNames,NFname)>=0 ) continue;
                                std::shared_ptr<NormFactor> munfactor = std::make_shared<NormFactor>(NFname);
                                LOG(DEBUG) << "    Adding EFT mu NormFactor: " << munfactor->fName << "\n";
                                TRExFitter::SYSTMAP[munfactor->fName] = munfactor->fName;
                                munfactor->SetNominalMinMax(1.0, -10., 10);
                                munfactor->fRegions = Common::ToVec(ireg->fName);

                                fFitter->fNormFactors.emplace_back( munfactor );
                                fFitter->fNormFactorNames.emplace_back( munfactor->fName );

                                EFTSM_s->AddNormFactor(munfactor);
                            }
                        }
                          }
                }
            }
        }
    }

    if (fFitter->fEFTConfig.GetSplitSamplesPerBin()) {
        for(const auto& eftconfset: EFTConfSets){
            auto& confsets = fParser->fConfSets.try_emplace("Sample").first->second;
            confsets.emplace_back(std::make_unique<ConfigSet>());
            auto *cs = confsets.back().get();
            cs->Set("Sample", eftconfset->fValue);

            for(const auto& iconf : eftconfset->fConfig){
                cs->SetConfig(iconf.fName,iconf.fValue,iconf.fFileName,iconf.fLineNumber);
            }
        }


        // Modify config sets to include new samples and replace all references to SM Reference with new samples
        for(const auto& [setName, confsets] : fParser->fConfSets) {
            for(auto& csPtr: confsets){
                ConfigSet& cs = *csPtr;
                for(auto& iconf : cs.fConfig){
                    if (!(iconf.fName == "Sample" || iconf.fName == "Samples" || iconf.fName == "AddSample" || iconf.fName == "AddSamples"  || iconf.fName == "Exclude"  || iconf.fName == "ExcludeRegionSample")) continue;
                    const std::vector<std::string> tmp_samples = Common::Vectorize(iconf.fValue,',');
                    std::string newvalue="";
                    for(const auto& tmp_sample : tmp_samples){
                        std::string newsample=tmp_sample;
                        std::map<std::string,std::string>::iterator it;
                        it = newEFTSampleLists.find(tmp_sample);
                        if (it != newEFTSampleLists.end()){
                            newsample=newEFTSampleLists[tmp_sample];
                        }

                        if (newvalue==""){
                            newvalue=newsample;
                        } else {
                            newvalue=newvalue+","+newsample;
                        }
                    }

                    if(iconf.fName == "Sample")iconf.fName = "Samples";
                    iconf.fValue=newvalue;
                }
            }
        }
    }

    fAvailableSamples = GetAvailableSamples();
    for(std::size_t i_smp = 0; i_smp < fFitter->fSamples.size(); ++i_smp) {
        for( auto asimovreplacements: fFitter->fSamples[i_smp]->fAsimovReplacementFor){
            if(asimovreplacements.first!=""){
                if (std::find(fSamples.begin(), fSamples.end(),asimovreplacements.second ) == fSamples.end()){
                    if (fAllowWrongRegionSample){
                        LOG(WARNING) << "Sample: " << fSamples[i_smp] << " has sample set up for AsimovReplacementFor that does not exist\n";
                    } else {
                        LOG(ERROR) << "Sample: " << fSamples[i_smp] << " has sample set up for AsimovReplacementFor that does not exist\n";
                        ++sc;
                    }
                }
                LOG(INFO) << "Creating sample " << asimovreplacements.first << "\n";
                std::shared_ptr<Sample> ca = fFitter->NewSample("customAsimov_"+asimovreplacements.first,Sample::SampleType::GHOST);
                ca->SetTitle("Pseudo-Data ("+asimovreplacements.first+")");
                ca->fUseSystematics = false;
            }
        }
    }
    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadNormFactorOptions(){
    int sc(0);

    const auto* normfactor_confsets = fParser->GetConfigSets("NormFactor");
    if(!normfactor_confsets){
        return sc;
    }
    for(const auto& confSetPtr : *normfactor_confsets){
        const ConfigSet *confSet = confSetPtr.get();

        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,Common::CheckName(confSet->GetValue()))>=0) continue;

        std::string param;
        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Common::Vectorize(samples_str,',');
        std::vector<std::string> regions = Common::Vectorize(regions_str,',');
        std::vector<std::string> exclude = Common::Vectorize(exclude_str,',');

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "NormFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
            } else {
                LOG(ERROR) << "NormFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
                ++sc;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "NormFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up for excluding that do not exist\n";
            } else {
                LOG(ERROR) << "NormFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up for excluding that do not exist\n";
                ++sc;
            }
        }

        if (samples.size() > 0 && !CheckPresence(samples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "NormFactor: " << Common::CheckName(confSet->GetValue()) << " has samples set up that do not exist\n";
            } else {
                LOG(ERROR) << "NormFactor: " << Common::CheckName(confSet->GetValue()) << " has samples set up that do not exist\n";
                ++sc;
            }
        }

        std::shared_ptr<NormFactor> nfactor = std::make_shared<NormFactor>(Common::CheckName(confSet->GetValue()));

        TRExFitter::SYSTMAP[nfactor->fName] = nfactor->fName;
        if( Common::FindInStringVector(fFitter->fNormFactorNames,nfactor->fName)<0 ){
            fFitter->fNormFactors.emplace_back( nfactor );
            fFitter->fNormFactorNames.emplace_back( nfactor->fName );
        }
        else{
            nfactor = fFitter->fNormFactors[ Common::FindInStringVector(fFitter->fNormFactorNames,nfactor->fName) ];
        }

        // Set NuisanceParameter = Name
        nfactor->fNuisanceParameter = Common::CheckName(nfactor->fName);
        TRExFitter::NPMAP[nfactor->fName] = Common::CheckName(nfactor->fName);

        if (SystHasProblematicName(nfactor->fNuisanceParameter)) {
            LOG(ERROR) << "NormFactor: + " << Common::CheckName(confSet->GetValue()) << " has a problematic NuisanceParameterName\n";
            ++sc;
        }

        // Set Constant
        param = confSet->Get("Constant");
        if(param != ""){
            nfactor->fConst = Common::StringToBoolean(param);
        }

        // Set Category
        param = confSet->Get("Category");
        if(param != "") nfactor->fCategory = Common::RemoveQuotes(param);

        // Set SubCategory
        param = confSet->Get("SubCategory");
        if(param != "") nfactor->fSubCategory = Common::RemoveQuotes(param);

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            nfactor->fTitle = Common::RemoveQuotes(param);
            TRExFitter::SYSTMAP[nfactor->fName] = nfactor->fTitle;
        }

        // Set TexTitle
        param = confSet->Get("TexTitle");
        if(param != "") TRExFitter::SYSTTEX[nfactor->fName] = Common::RemoveQuotes(param);

        {
            // set defaults for nominal value and bounds
            double nominal{1};
            double min{0};
            double max{10};
            // Set Min
            param = confSet->Get("Min");
            if(param!="") min = atof(param.c_str());

            // Set Max
            param = confSet->Get("Max");
            if(param!="") max = atof(param.c_str());

            // Set Nominal
            param = confSet->Get("Nominal");
            if(param!="") nominal = atof(param.c_str());
            nfactor->SetNominalMinMax(nominal, min, max);
        }

        // Set Expression
        param = confSet->Get("Expression");
        if(param!=""){
            std::vector<std::string> v = Common::Vectorize(param,':');
            if (v.size() < 2){
                LOG(ERROR) << "You specified 'Expression' option but did not provide 2 parameters. Please check this\n";
                ++sc;
            } else {
                nfactor->fExpression = std::make_pair(v[0],v[1]);
                // title will contain the expression FIXME
                nfactor->fTitle = "Expression_" + v[0];
                TRExFitter::SYSTMAP[nfactor->fName] = "Expression_" + v[0];
                // nuis-par will contain the nuis-par of the norm factor the expression depends on FIXME
                nfactor->fNuisanceParameter = "Expression_" + v[1];
                TRExFitter::NPMAP[nfactor->fName] = "Expression_" + v[1];
                // set nominal, min and max according to the norm factor the expression depends on FIXME
                for(const auto& nf : fFitter->fNormFactors){
                    if(nf->fNuisanceParameter == "Expression_" + v[1]) {
                        nfactor->SetNominalMinMax(nf->GetNominal(), nf->GetMin(), nf->GetMax());
                    }
                }
            }
        }

        // Set tau
        param = confSet->Get("Tau");
        if (param != "") {
            nfactor->fTau = Common::convertStoNum<float>(param);
        }

        // Set CorrelateUnfolding
        param = confSet->Get("CorrelateUnfolding");
        if (param != "") {
            if (fFitter->fFitType != TRExFit::FitType::UNFOLDING) {
                LOG(WARNING) << "\"CorrelateUnfolding\" has impact only on Unfolding type of a fit\n";
            } else {
                param = Common::RemoveQuotes(param);
                auto itr = std::find_if(fFitter->fUnfolding.begin(), fFitter->fUnfolding.end(), [&param](const std::unique_ptr<Unfolding>& element){return element->fName == param;});
                if (itr == fFitter->fUnfolding.end()) {
                    LOG(ERROR) << "\"CorrelateUnfolding\": " << param << " doesnt correspond to any known Unfolding object\n";
                    ++sc;
                } else {
                    nfactor->fCorrelateUnfolding = param;
                }
            }
        }

        // save list of
        if (regions.size() == 0 || exclude.size() == 0){
            LOG(ERROR) << "Region or exclude region size is equal to zero. Please check this\n";
            ++sc;
            return sc;
        }
        if(regions[0] == "all") {
            nfactor->fRegions = GetAvailableRegions();
        } else {
            nfactor->fRegions = regions;
        }
        if(exclude[0] != "") nfactor->fExclude = exclude;

        // this is needed only for the unfolding in case of custom NFs
        nfactor->fSamples = samples;
        // attach the syst to the proper samples
        for(auto& isample : fFitter->fSamples) {
            if(isample->fType == Sample::SampleType::DATA) continue;
            if((samples[0]=="all" || Common::FindInStringVector(samples, isample->fName)>=0 )
                && (exclude[0]==""    || Common::FindInStringVector(exclude, isample->fName)<0 ) ){
                isample->AddNormFactor(nfactor);
            }
        }
    }
    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadShapeFactorOptions() {
    int sc(0);

    const auto* shapefactor_confsets = fParser->GetConfigSets("ShapeFactor");
    if(!shapefactor_confsets){
        return sc;
    }
    for(const auto& confSetPtr : *shapefactor_confsets){
        const ConfigSet *confSet = confSetPtr.get();

        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,Common::CheckName(confSet->GetValue()))>=0) continue;
        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Common::Vectorize(samples_str,',');
        std::vector<std::string> regions = Common::Vectorize(regions_str,',');
        std::vector<std::string> exclude = Common::Vectorize(exclude_str,',');

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "ShapeFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
            } else {
                LOG(ERROR) << "ShapeFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
                ++sc;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "ShapeFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up for excluding that do not exist\n";
            } else {
                LOG(ERROR) << "ShapeFactor: " << Common::CheckName(confSet->GetValue()) << " has regions set up for excluding that do not exist\n";
                ++sc;
            }
        }

        if (samples.size() > 0 && !CheckPresence(samples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "ShapeFactor: " << Common::CheckName(confSet->GetValue()) << " has samples set up that do not exist\n";
            } else {
                LOG(ERROR) << "ShapeFactor: " << Common::CheckName(confSet->GetValue()) << " has samples set up that do not exist\n";
                ++sc;
            }
        }

        std::shared_ptr<ShapeFactor> sfactor = std::make_shared<ShapeFactor>(Common::CheckName(confSet->GetValue()));
        if( Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName)<0 ){
            fFitter->fShapeFactors.emplace_back( sfactor );
            fFitter->fShapeFactorNames.emplace_back( sfactor->fName );
        }
        else{
            sfactor = fFitter->fShapeFactors[ Common::FindInStringVector(fFitter->fShapeFactorNames,sfactor->fName) ];
        }

        // Set NuisanceParameter = Name
        sfactor->fNuisanceParameter = Common::CheckName(sfactor->fName);
        TRExFitter::NPMAP[sfactor->fName] = Common::CheckName(sfactor->fName);

        if (SystHasProblematicName(sfactor->fNuisanceParameter)) {
            LOG(ERROR) << "ShapeFactor: + " << Common::CheckName(confSet->GetValue()) << " has a problematic NuisanceParameterName\n";
            ++sc;
        }

        // Set Constant
        std::string param = confSet->Get("Constant");
        if(param != ""){
            sfactor->fConst = Common::StringToBoolean(param);
        }

        // Set Category
        param = confSet->Get("Category");
        if(param!="") sfactor->fCategory = Common::RemoveQuotes(param);

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            sfactor->fTitle = Common::RemoveQuotes(param);
            TRExFitter::SYSTMAP[sfactor->fName] = sfactor->fTitle;
        }
        param = confSet->Get("TexTitle");
        if(param != "") TRExFitter::SYSTTEX[sfactor->fName] = Common::RemoveQuotes(param);

        // Set Min
        param = confSet->Get("Min");
        if(param!="") sfactor->fMin = atof(param.c_str());

        // Set Max
        param = confSet->Get("Max");
        if(param!="") sfactor->fMax = atof(param.c_str());

        // Set Nominal
        param = confSet->Get("Nominal");
        if(param!="") sfactor->fNominal = atof(param.c_str());

        if (regions.size() == 0 || exclude.size() == 0){
            LOG(ERROR) << "Region or exclude region size is equal to zero. Please check this\n";
            ++sc;
            return sc;
        }
        // save list of
        if(regions[0] == "all") {
            sfactor->fRegions = GetAvailableRegions();
        } else {
            sfactor->fRegions = regions;
        }
        if(exclude[0]!="")    sfactor->fExclude = exclude;
        // attach the syst to the proper samples
        for(auto& isample : fFitter->fSamples) {
            if(isample->fType == Sample::SampleType::DATA) continue;
            if(   (samples[0]=="all" || Common::FindInStringVector(samples, isample->fName)>=0 )
               && (exclude[0]==""    || Common::FindInStringVector(exclude, isample->fName)<0 ) ){
                isample->AddShapeFactor(sfactor);
                sfactor->fSamples.emplace_back(isample->fName);
            }
        }

        param = confSet->Get("Expression");
        if (param != "") {
            if (sfactor->fSamples.size() > 1) {
                LOG(ERROR) << "You need to have only one sample if the \"Expression\" option is used\n";
                ++sc;
                return sc;
            }
            if (sfactor->fRegions.size() > 1) {
                LOG(ERROR) << "You need to have only one regiond if the \"Expression\" option is used\n";
                ++sc;
                return sc;
            }
            if ((sfactor->fRegions.at(0) == "all") && (this->GetAvailableRegions().size() > 1)) {
                LOG(ERROR) << "You need to have only one region if the \"Expression\" option is used\n";
                ++sc;
                return sc;
            }
            if ((sfactor->fSamples.at(0) == "all") && (this->GetAvailableSamples().size() > 1)) {
                LOG(ERROR) << "You need to have only one sample if the \"Expression\" option is used\n";
                ++sc;
                return sc;
            }
            const std::vector<std::string> oneParam = Common::Vectorize(param, '@');
            for (const auto& iparam : oneParam) {
                const std::vector<std::string> tmp = Common::Vectorize(iparam, ':');
                if (tmp.size() != 2) {
                    LOG(ERROR) << "Wrong \"Expression\" option syntax\n";
                    ++sc;
                    break;
                }
                sfactor->fExpression.emplace_back(std::make_pair(tmp.at(0), tmp.at(1)));
            }

        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadSystOptions(){
    int sc(0);

    if (fOnlySystematics.size() > 0){
        std::vector<std::string> availableSysts = GetAvailableSysts();
        if (!CheckPresence(fOnlySystematics, availableSysts)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "You set systematics that do not exist in your command line options\n";
            } else {
                LOG(ERROR) << "You set systematics that do not exist in your command line options\n";
                ++sc;
            }
        }
    }

    if (fFitter->fStatOnly) {
        int typed = Systematic::OVERALL;
        std::shared_ptr<Systematic> sysd = std::make_shared<Systematic>("Dummy",typed);
        sysd->fOverallUp   = 0.;
        sysd->fOverallDown = -0.;
        sysd->fScaleUp   = 1.;
        sysd->fScaleDown   = 1.;
        fFitter->fSystematics.emplace_back( sysd );
        TRExFitter::SYSTMAP[sysd->fName] = "Dummy";
        for(auto& isample : fFitter->fSamples) {
            if(isample->fType == Sample::SampleType::SIGNAL ) {
                isample->AddSystematic(sysd);
            }
        }
    }

    const auto* syst_confsets = fParser->GetConfigSets("Systematic");
    if(!syst_confsets){
        return sc;
    }
    for(const auto& confSetPtr : *syst_confsets){
        const ConfigSet *confSet = confSetPtr.get();

        std::string param = "";
        if(fOnlySystematics.size()>0 && Common::FindInStringVector(fOnlySystematics,Common::CheckName(confSet->GetValue()))<0) continue;
        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,Common::CheckName(confSet->GetValue()))>=0) continue;
        std::string samples_str = confSet->Get("Samples");
        std::string regions_str = confSet->Get("Regions");
        std::string exclude_str = confSet->Get("Exclude");
        std::string excludeRegionSample_str = confSet->Get("ExcludeRegionSample");
        if(samples_str=="") samples_str = "all";
        if(regions_str=="") regions_str = "all";
        std::vector<std::string> samples = Common::Vectorize(samples_str,',');
        std::vector<std::string> regions = Common::Vectorize(regions_str,',');
        std::vector<std::string> exclude = Common::Vectorize(exclude_str,',');

        //
        bool irrelevantSyst = false;
        if (fOnlySamples.size()>0 && samples_str != "all"){
            // if command line Samples are specified
            // loop over samples relevant to this syst, and check they are specified in the command line
            // systematic is irrelevant if all its samples are not in command line
            irrelevantSyst = std::all_of( samples.begin(), samples.end(),
                                          [&](const std::string& systSample) {
                                            return Common::FindInStringVector(fOnlySamples, Common::CheckName(systSample))<0;
                                            }
                                        );

        }
        if (irrelevantSyst) continue;
        if (fOnlyRegions.size()>0 && regions_str != "all" ){
            // if command line Regions are specified
            // loop over regions relevant to this syst, and check they are specified in the command line
            // systematic is irrelevant if all its regions are not in command line
            irrelevantSyst = std::all_of( regions.begin(), regions.end(),
                                          [&](const std::string& systRegion) {
                                            return Common::FindInStringVector(fOnlyRegions, Common::CheckName(systRegion))<0;
                                            }
                                        );
        }
        if (irrelevantSyst) continue;
        if (fToExclude.size()>0 && (regions_str != "all" ||  samples_str != "all")){
            // if command line Exclude is specified
            // loop over samples and regions relevant to this syst, and check they are being excluded
            if (regions_str != "all"){
                // systematic is irrelevant if all its regions are being excluded
                irrelevantSyst = std::all_of( regions.begin(), regions.end(),
                                            [&](const std::string& systRegion) {
                                                return Common::FindInStringVector(fToExclude, Common::CheckName(systRegion))>0;
                                                }
                                            );
            }
            if (irrelevantSyst) continue;
            if (samples_str != "all"){
                // systematic is irrelevant if all its samples are being excluded
                irrelevantSyst = std::all_of( samples.begin(), samples.end(),
                                [&](const std::string& systSample) {
                                return Common::FindInStringVector(fToExclude, Common::CheckName(systSample))>0;
                                }
                            );
            }
            if (irrelevantSyst) continue;
        }

        if(samples_str!="all" && confSet->Get("DummyForSamples")!=""){
            std::vector<std::string> addSamples = Common::Vectorize(confSet->Get("DummyForSamples"),',');
            for(const auto& addSmp : addSamples){
                if(Common::FindInStringVector(samples,addSmp)<0){
                    samples.emplace_back(addSmp);
                }
            }
        }

        if (regions.size() > 0 && !CheckPresence(regions, fAvailableRegions)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
            } else {
                LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up that do not exist\n";
                ++sc;
            }
        }

        if (exclude.size() > 0 && !CheckPresence(exclude, fAvailableRegions, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples/regions set up for excluding that do not exist\n";
            } else {
                LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples/regions set up for excluding that do not exist\n";
                ++sc;
            }
        }

        if (samples.size() > 0 && !CheckPresence(samples, fAvailableSamples)){
            if (fAllowWrongRegionSample){
                LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up that do not exist\n";
            } else {
                LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up that do not exist\n";
                ++sc;
            }
        }

        fExcludeRegionSample = Common::Vectorize(excludeRegionSample_str,',');
        int type = Systematic::HISTO;

        // Set type
        param = confSet->Get("Type");
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);

        if(param == "OVERALL") type = Systematic::OVERALL;
        else if(param == "SHAPE") type = Systematic::SHAPE;
        else if (param == "STAT") type = Systematic::STAT;

        std::string decorrelate = confSet->Get("Decorrelate");

        std::shared_ptr<Systematic> sys = std::make_shared<Systematic>(Common::CheckName(confSet->GetValue()),type);
        TRExFitter::SYSTMAP[sys->fName] = sys->fTitle;
        if(param == "OVERALL") sys->fIsNormOnly=true;

        sys->fSamples = samples;

        // SetCategory
        param = confSet->Get("Category");
        if(param != ""){
            sys->fCategory = Common::RemoveQuotes(param);
            sys->fSubCategory = Common::RemoveQuotes(param); //SubCategory defaults to the Category setting, if the Category is explicitly set
        }

        // CombineName
        param = confSet->Get("CombineName");
        if(param != ""){
            sys->fCombineName = Common::RemoveQuotes(param);
        }

        // CombineType
        param = confSet->Get("CombineType");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            const std::string tmp = Common::RemoveQuotes(param);

            if (tmp == "STANDARDDEVIATION") {
                sys->fCombineType = Systematic::COMBINATIONTYPE::STANDARDDEVIATION;
            } else if (tmp == "ENVELOPE") {
                sys->fCombineType = Systematic::COMBINATIONTYPE::ENVELOPE;
            } else if (tmp == "STANDARDDEVIATIONNODDOF") {
                sys->fCombineType = Systematic::COMBINATIONTYPE::STANDARDDEVIATIONNODDOF;
            } else if (tmp == "HESSIAN") {
                sys->fCombineType = Systematic::COMBINATIONTYPE::HESSIAN;
            } else if (tmp == "SUMINSQUARES") {
                sys->fCombineType = Systematic::COMBINATIONTYPE::SUMINSQUARES;
            } else {
                LOG(WARNING) << "You specified 'CombineType' option but did not provide valid parameter. Using default (ENVELOPE)\n";
                sys->fCombineType = Systematic::COMBINATIONTYPE::ENVELOPE;
            }
        }

        // SetSubCategory
        param = confSet->Get("SubCategory");
        if (param != ""){
            sys->fSubCategory = Common::RemoveQuotes(param); // note this needs to happen after Category was set, in order to overwrite the default if required
        }

        // Set IsFreeParameter
        // Experimental
        param = confSet->Get("IsFreeParameter");
        if(param != ""){
            sys->fIsFreeParameter = Common::StringToBoolean(param);
        }

        // Set IsCorrelated
        param = confSet->Get("IsCorrelated");
        if(param != ""){
            sys->fIsCorrelated = Common::StringToBoolean(param);
        }

        // Set StoredName
        // New: name to use when writing / reading the Histograms file
        param = confSet->Get("StoredName");
        if(param != "") sys->fStoredName = Common::RemoveQuotes(param);

        if(type==Systematic::HISTO || type==Systematic::SHAPE){
            bool hasUp(false);
            bool hasDown(false);

            if(fFitter->fInputType==0){ // HIST input
                param = confSet->Get("HistoPathUp");
                if(param!=""){
                    sys->fHistoPathsUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoPathDown");
                if(param!=""){
                    sys->fHistoPathsDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoPathSufUp");
                if(param!=""){
                    sys->fHistoPathSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoPathSufDown");
                if(param!=""){
                    sys->fHistoPathSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoFileUp");
                if(param!=""){
                    sys->fHistoFilesUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoFileDown");
                if(param!=""){
                    sys->fHistoFilesDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoFilesUp");
                if(param!=""){
                    sys->fHistoFilesUp = Common::Vectorize(param, ',');
                    hasUp   = true;
                }
                param = confSet->Get("HistoFilesDown");
                if(param!=""){
                    sys->fHistoFilesDown = Common::Vectorize(param, ',');
                    hasDown = true;
                }
                param = confSet->Get("HistoFileSufUp");
                if(param!=""){
                    sys->fHistoFileSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoFileSufDown");
                if(param!=""){
                    sys->fHistoFileSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoNameUp");
                if(param!=""){
                    sys->fHistoNamesUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoNameDown");
                if(param!=""){
                    sys->fHistoNamesDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoFolderNameUp");
                if(param!=""){
                    sys->fHistoFolderNamesUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoFolderNameDown");
                if(param!=""){
                    sys->fHistoFolderNamesDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoFolderSubtractNameUp");
                if(param!=""){
                    sys->fHistoFolderSubtractNamesUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("HistoFolderSubtractNameDown");
                if(param!=""){
                    sys->fHistoFolderSubtractNamesDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("HistoNameSufUp");
                if(param!=""){
                    sys->fHistoNameSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoNameSufDown");
                if(param!=""){
                    sys->fHistoNameSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoPathsUpSubtractSample");
                if(param!=""){
                    sys->fHistoPathsUpSubtractSample = Common::Vectorize(param,',');
                    hasUp = true;
                }
                param = confSet->Get("HistoFilesUpSubtractSample");
                if(param!=""){
                    sys->fHistoFilesUpSubtractSample = Common::Vectorize(param,',');
                    hasUp = true;
                }
                param = confSet->Get("HistoFileSufUpSubtractSample");
                if(param!=""){
                    sys->fHistoFileSufUpSubtractSample =  Common::RemoveQuotes(param);
                    hasUp = true;
                }
                param = confSet->Get("HistoNamesUpSubtractSample");
                if(param!=""){
                    sys->fHistoNamesUpSubtractSample = Common::Vectorize(param,',');
                    hasUp = true;
                }
                param = confSet->Get("HistoPathsDownSubtractSample");
                if(param!=""){
                    sys->fHistoPathsDownSubtractSample = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("HistoFilesDownSubtractSample");
                if(param!=""){
                    sys->fHistoFilesDownSubtractSample = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("HistoFileSufDownSubtractSample");
                if(param!=""){
                    sys->fHistoFileSufDownSubtractSample =  Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("HistoNamesDownSubtractSample");
                if(param!=""){
                    sys->fHistoNamesDownSubtractSample = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("HistoPathUpSubtractSample");
                if(param!=""){
                    sys->fHistoPathsUpSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasUp = true;
                }
                param = confSet->Get("HistoFileUpSubtractSample");
                if(param!=""){
                    sys->fHistoFilesUpSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasUp = true;
                }
                param = confSet->Get("HistoNameUpSubtractSample");
                if(param!=""){
                    sys->fHistoNamesUpSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasUp = true;
                }
                param = confSet->Get("HistoPathDownSubtractSample");
                if(param!=""){
                    sys->fHistoPathsDownSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("HistoFileDownSubtractSample");
                if(param!=""){
                    sys->fHistoFilesDownSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("HistoNameDownSubtractSample");
                if(param!=""){
                    sys->fHistoNamesDownSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("HistoNameSufUpSubtractSample");
                if(param!=""){
                    sys->fHistoNameSufUpSubtractSample = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("HistoNameSufDownSubtractSample");
                if(param!=""){
                    sys->fHistoNameSufDownSubtractSample = Common::RemoveQuotes(param);
                    hasDown = true;
                }
            }
            else if(fFitter->fInputType==1){ // NTUP option
                param = confSet->Get("NtuplePathUp");
                if(param!=""){
                    sys->fNtuplePathsUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathDown");
                if(param!=""){
                    sys->fNtuplePathsDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("NtuplePathsUp");
                if(param!=""){
                    sys->fNtuplePathsUp = Common::Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathsDown");
                if(param!=""){
                    sys->fNtuplePathsDown = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtuplePathSufUp");
                if(param!=""){
                    sys->fNtuplePathSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtuplePathSufDown");
                if(param!=""){
                    sys->fNtuplePathSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("FriendPathUp");
                if(param!=""){
                    sys->fFriendPathsUp.push_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("FriendPathDown");
                if(param!=""){
                    sys->fFriendPathsDown.push_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("FriendPathsUp");
                if(param!=""){
                    sys->fFriendPathsUp = Common::Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("FriendPathsDown");
                if(param!=""){
                    sys->fFriendPathsDown = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("FriendPathSufUp");
                if(param!=""){
                    sys->fFriendPathSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("FriendPathSufDown");
                if(param!=""){
                    sys->fFriendPathSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("FriendPathUpSubtractSample");
                if(param!=""){
                    sys->fFriendPathsUpSubtractSample.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("FriendPathDownSubtractSample");
                if(param!=""){
                    sys->fFriendPathsDownSubtractSample.emplace_back(Common::RemoveQuotes(param));
                    hasDown   = true;
                }
                param = confSet->Get("FriendPathsUpSubtractSample");
                if(param!=""){
                    sys->fFriendPathsUpSubtractSample = Common::Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("FriendPathsDownSubtractSample");
                if(param!=""){
                    sys->fFriendPathsDownSubtractSample = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileUp");
                if(param!=""){
                    sys->fNtupleFilesUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFileDown");
                if(param!=""){
                    sys->fNtupleFilesDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("NtupleFilesUp");
                if(param!=""){
                    sys->fNtupleFilesUp = Common::Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFilesDown");
                if(param!=""){
                    sys->fNtupleFilesDown = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileSufUp");
                if(param!=""){
                    sys->fNtupleFileSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtupleFileSufDown");
                if(param!=""){
                    sys->fNtupleFileSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameUp");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesUp.emplace_back(Common::RemoveQuotes(param));
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNameDown");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesDown.emplace_back(Common::RemoveQuotes(param));
                    hasDown = true;
                }
                param = confSet->Get("NtupleNamesUp");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesUp = Common::Vectorize(param,',');
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNamesDown");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesDown = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameSufUp");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNameSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("NtupleNameSufDown");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNameSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("WeightUp");
                if(param!=""){
                    sys->fWeightUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("WeightDown");
                if(param!=""){
                    sys->fWeightDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("WeightSufUp");
                if(param!=""){
                    sys->fWeightSufUp = Common::RemoveQuotes(param);
                    hasUp   = true;
                }
                param = confSet->Get("WeightSufDown");
                if(param!=""){
                    sys->fWeightSufDown = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("IgnoreWeight");
                if(param!=""){
                    sys->fIgnoreWeight = Common::RemoveQuotes(param);
                }
                param = confSet->Get("NtuplePathsUpSubtractSample");
                if(param!=""){
                    sys->fNtuplePathsUpSubtractSample = Common::Vectorize(param,',');
                    hasUp = true;
                }
                param = confSet->Get("NtupleFilesUpSubtractSample");
                if(param!=""){
                    sys->fNtupleFilesUpSubtractSample = Common::Vectorize(param,',');
                    hasUp = true;
                }
                param = confSet->Get("NtupleNamesUpSubtractSample");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesUpSubtractSample = Common::Vectorize(param,',');
                    hasUp = true;
                }
                param = confSet->Get("NtuplePathsDownSubtractSample");
                if(param!=""){
                    sys->fNtuplePathsDownSubtractSample = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleFilesDownSubtractSample");
                if(param!=""){
                    sys->fNtupleFilesDownSubtractSample = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtupleNamesDownSubtractSample");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesDownSubtractSample = Common::Vectorize(param,',');
                    hasDown = true;
                }
                param = confSet->Get("NtuplePathUpSubtractSample");
                if(param!=""){
                    sys->fNtuplePathsUpSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasUp = true;
                }
                param = confSet->Get("NtupleFileUpSubtractSample");
                if(param!=""){
                    sys->fNtupleFilesUpSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasUp = true;
                }
                param = confSet->Get("NtupleNameUpSubtractSample");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesUpSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasUp = true;
                }
                param = confSet->Get("NtuplePathDownSubtractSample");
                if(param!=""){
                    sys->fNtuplePathsDownSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileDownSubtractSample");
                if(param!=""){
                    sys->fNtupleFilesDownSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameDownSubtractSample");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNamesDownSubtractSample.emplace_back( Common::RemoveQuotes(param) );
                    hasDown = true;
                }
                param = confSet->Get("NtupleFileSufUpSubtractSample");
                if(param!=""){
                    sys->fNtupleFileSufUpSubtractSample = Common::RemoveQuotes(param);
                        hasUp   = true;
                }
                param = confSet->Get("NtupleFileSufDownSubtractSample");
                if(param!=""){
                    sys->fNtupleFileSufDownSubtractSample = Common::RemoveQuotes(param);
                    hasDown = true;
                }
                param = confSet->Get("NtupleNameSufUpSubtractSample");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNameSufUpSubtractSample = Common::RemoveQuotes(param);
                        hasUp   = true;
                }
                param = confSet->Get("NtupleNameSufDownSubtractSample");
                if(param!=""){
                    Common::CheckTreeName(param);
                    sys->fNtupleNameSufDownSubtractSample = Common::RemoveQuotes(param);
                    hasDown = true;
                }
            }
            sys->fHasUpVariation   = hasUp  ;
            sys->fHasDownVariation = hasDown;

            // Set Symmetrisation
            param = confSet->Get("Symmetrisation");
            if(param != ""){
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if(param == "ONESIDED"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEONESIDED;
                }
                else if(param == "TWOSIDED"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZETWOSIDED;
                }
                else if(param == "ABSMEAN"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEABSMEAN;
                }
                else if(param == "MAXIMUM"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEMAXIMUM;
                }
                else if(param == "MAXIMUMKEEPSIGN"){
                    sys->fSymmetrisationType = HistoTools::SYMMETRIZEMAXIMUMKEEPSIGN;
                }
                else {
                    LOG(ERROR) << "Symmetrisation scheme is not recognized ... \n";
                    ++sc;
                }
            }

            // Set Smoothing
            param = confSet->Get("Smoothing");
            if (param != "") {
                sys->fSmoothType = atoi(param.c_str());
            }

            param = confSet->Get("PreSmoothing");
            if (param != ""){
                sys->fPreSmoothing = Common::StringToBoolean(param);
            }

            param = confSet->Get("SmoothingOption");
            if(param != ""){
                sys->fSampleSmoothing = true;
                std::transform(param.begin(), param.end(), param.begin(), ::toupper);
                if( param == "MAXVARIATION" ) sys->fSampleSmoothOption = HistoTools::SmoothOption::MAXVARIATION;
                else if (param == "TTBARRESONANCE") sys->fSampleSmoothOption = HistoTools::SmoothOption::TTBARRESONANCE;
                else if (param == "COMMONTOOLSMOOTHMONOTONIC") sys->fSampleSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHMONOTONIC;
                else if (param == "COMMONTOOLSMOOTHPARABOLIC") sys->fSampleSmoothOption = HistoTools::SmoothOption::COMMONTOOLSMOOTHPARABOLIC;
                else if (param == "TCHANNEL") sys->fSampleSmoothOption = HistoTools::SmoothOption::TCHANNEL;
                else if (param == "KERNELRATIOUNIFORM") sys->fSampleSmoothOption = HistoTools::SmoothOption::KERNELRATIOUNIFORM;
                else if (param == "KERNELDELTAGAUSS") sys->fSampleSmoothOption = HistoTools::SmoothOption::KERNELDELTAGAUSS;
                else if (param == "KERNELRATIOGAUSS") sys->fSampleSmoothOption = HistoTools::SmoothOption::KERNELRATIOGAUSS;
                else {
                    LOG(WARNING) << "You specified 'SmoothingOption' option but you didn't provide valid input. Using default from job block\n";
                    sys->fSampleSmoothing = false;
                }
            }
        } // end of if(type==Systematic::HISTO || type==Systematic::SHAPE)

        // Set OverallUp
        param = confSet->Get("OverallUp");
        if (param != "") {
            if (type != Systematic::SystType::OVERALL) {
                LOG(ERROR) << "Systematic: " << confSet->GetValue() << " has OverallUp set, but the type is not OVERALL!\n";
                ++sc;
            } else {
                sys->fOverallUp = atof(param.c_str());
            }
        }

        // Set OverallUDown
        param = confSet->Get("OverallDown");
        if (param != "") {
            if (type != Systematic::SystType::OVERALL) {
                LOG(ERROR) << "Systematic: " << confSet->GetValue() << " has OverallDown set, but the type is not OVERALL!\n";
                ++sc;
            } else {
                sys->fOverallDown = atof(param.c_str());
            }
        }

        // Set ScaleUp
        param = confSet->Get("ScaleUp");
        if(param!=""){
            std::vector < std::string > temp_vec = Common::Vectorize(param,',',false);
            if(temp_vec.size()==1 && Common::Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleUp = atof(param.c_str());
            }
            else{
                for (const std::string& ireg : temp_vec){
                    std::vector < std::string > reg_value = Common::Vectorize(ireg,':');
                    if(reg_value.size()==2){
                        sys->fScaleUpRegions.insert( std::pair < std::string, double >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }

        // Set ScaleDown
        param = confSet->Get("ScaleDown");
        if(param!=""){
            std::vector < std::string > temp_vec = Common::Vectorize(param,',',false);
            if(temp_vec.size()==1 && Common::Vectorize(temp_vec[0],':').size()==1){
                sys->fScaleDown = atof(param.c_str());
            }
            else{
                for (const std::string& ireg : temp_vec){
                    std::vector < std::string > reg_value = Common::Vectorize(ireg,':');
                    if(reg_value.size()==2){
                        sys->fScaleDownRegions.insert( std::pair < std::string, double >( reg_value[0], atof(reg_value[1].c_str()) ) );
                    }
                }
            }
        }

        // Set SampleUp
        // a new way of defining systematics, just comparing directly with GHOST samples
        // --> this can be used only if this systematic is applied to a single sample
        param = confSet->Get("SampleUp");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), Common::RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in SampleUp that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in SampleUp that do not exist\n";
                    ++sc;
                }
            }
            sys->fSampleUp = Common::RemoveQuotes(param);
        }

        // Set SampleDown
        param = confSet->Get("SampleDown");
        if(param!=""){
            if (std::find(fSamples.begin(), fSamples.end(), Common::RemoveQuotes(param)) == fSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in SampleDown that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in SampleDown that do not exist\n";
                    ++sc;
                }
            }
            sys->fSampleDown = Common::RemoveQuotes(param);
        }

        // Set ReferenceSample
        // this to obtain syst variation relatively to given sample
        param = confSet->Get("ReferenceSample");
        if(param!=""){
            if (std::find(fAvailableSamples.begin(), fAvailableSamples.end(), Common::RemoveQuotes(param)) == fAvailableSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in ReferenceSample that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in ReferenceSample that do not exist\n";
                    ++sc;
                }
            }
            sys->fReferenceSample = Common::RemoveQuotes(param);
        }

        // Set ReferenceSmoothing
        // this is to obtain syst variation relatively to given sample
        param = confSet->Get("ReferenceSmoothing");
        if(param!=""){
            if (std::find(fAvailableSamples.begin(), fAvailableSamples.end(), Common::RemoveQuotes(param)) == fAvailableSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in ReferenceSmoothing that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in ReferenceSmoothing that do not exist\n";
                    ++sc;
                }
            }
            if (std::find(samples.begin(), samples.end(), Common::RemoveQuotes(param)) == samples.end()){
                LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " requires that the ReferenceSample appears in Samples for this systematic\n";
                ++sc;
            }
            sys->fReferenceSmoothing = Common::RemoveQuotes(param);
        }

        // Set ReferencePruning
        // this allows to prune wrt to a sample and then apply the result to all samples
        param = confSet->Get("ReferencePruning");
        if(param!=""){
            if (std::find(fAvailableSamples.begin(), fAvailableSamples.end(), Common::RemoveQuotes(param)) == fAvailableSamples.end()){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in ReferencePruning that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in ReferencePruning that do not exist\n";
                    ++sc;
                }
            }
            sys->fReferencePruning = Common::RemoveQuotes(param);
        }

        // Set DropShapeIn
        param = confSet->Get("DropShapeIn");
        if(param!="") {
            std::vector<std::string> tmp = Common::Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions)){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up in DropShapeIn that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up in DropShapeIn that do not exist\n";
                    ++sc;
                }
            }
            sys->fDropShapeIn = tmp;
        }

        // Set DropNorm
        param = confSet->Get("NoPruning");
        if (param!="") {
            sys->fNoPruning = Common::StringToBoolean(param);
        }

        // Set DropNorm
        param = confSet->Get("DropNorm");
        std::string tmp_param = param;
        if(param!=""){
            std::transform(tmp_param.begin(), tmp_param.end(), tmp_param.begin(), ::tolower);
            if (tmp_param == "all") param = tmp_param;
            std::vector<std::string> tmp = Common::Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up in DropNorm that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up in DropNorm that do not exist\n";
                    ++sc;
                }
            }
            sys->fDropNormIn = tmp;
        }

        // Set DropNormSpecial
        param = confSet->Get("DropNormSpecial");
        if(param!=""){
            std::vector<std::string> tmp = Common::Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up in DropNorm that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has regions set up in DropNorm that do not exist\n";
                    ++sc;
                }
            }
            sys->fDropNormSpecialIn = tmp;
        }

        // Set KeepNormForSamples
        param = confSet->Get("KeepNormForSamples");
        if(param!="") {
            std::vector<std::string> tmp = Common::Vectorize(param,',');
            if (tmp.size() > 0 && !CheckPresence(tmp, fAvailableRegions, fAvailableSamples)){
                if (fAllowWrongRegionSample){
                    LOG(WARNING) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in KeepNormForSamples that do not exist\n";
                } else {
                    LOG(ERROR) << "Systematic: " << Common::CheckName(confSet->GetValue()) << " has samples set up in KeepNormForSamples that do not exist\n";
                    ++sc;
                }
            }
            sys->fKeepNormForSamples = tmp;
        }

        if (regions.size() == 0 || exclude.size() == 0){
            LOG(ERROR) << "Region or exclude region size is equal to zero. Please check this\n";
            ++sc;
            return sc;
        }

        // Set DummyForSamples
        param = confSet->Get("DummyForSamples");
        if(param!="") {
            std::vector<std::string> tmp = Common::Vectorize(param,',');
            sys->fDummyForSamples = tmp;
        }

        // Set ForceShape
        param = confSet->Get("ForceShape");
        if(param!="") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "NOSHAPE") {
                sys->fForceShape = HistoTools::FORCESHAPETYPE::NOSHAPE;
            } else if (param == "LINEAR") {
                sys->fForceShape = HistoTools::FORCESHAPETYPE::LINEAR;
            } else if (param == "TRIANGULAR") {
                sys->fForceShape = HistoTools::FORCESHAPETYPE::TRIANGULAR;
            } else {
                LOG(WARNING) << "You specified 'ForceShape' option but did not provide valid parameter. Using default (NOSHAPE)\n";
                sys->fForceShape = HistoTools::FORCESHAPETYPE::NOSHAPE;
            }
        }

        if(regions[0]!="all") sys->fRegions = regions;
        if(exclude[0]!="")    sys->fExclude = exclude;

        for (const std::string& ier : fExcludeRegionSample){
            std::vector<std::string> pair_ERS = Common::Vectorize(ier,':');
            if (pair_ERS.size()==2){
                sys->fExcludeRegionSample.emplace_back(pair_ERS);
            }
        }
        if ( decorrelate == "" && type != Systematic::STAT) {
            sc += SetSystNoDecorelate(confSet, sys, samples, exclude);
        }
        else if (decorrelate == "REGION" || type == Systematic::STAT)  {
            sc+= SetSystRegionDecorelate(confSet, sys, samples, exclude, regions, type);
        }
        else if (decorrelate == "SAMPLE")  {
            sc+= SetSystSampleDecorelate(confSet, sys, samples, exclude);
        }
        else if (decorrelate == "SHAPEACC")  {
            sc+= SetSystShapeDecorelate(confSet, sys, samples, exclude);
        }
        else {
            LOG(ERROR) << "Decorrelate option: " << decorrelate  << "  not supported ...\n";
            LOG(ERROR) << "       PLEASE USE ONLY: REGION, SAMPLE, SHAPEACC\n";
            ++sc;
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystNoDecorelate(const ConfigSet *confSet, std::shared_ptr<Systematic> sys, const std::vector<std::string>& samples, const std::vector<std::string>& exclude){
    int sc(0);

    fFitter->fSystematics.emplace_back( sys );

    // Set NuisanceParameter
    std::string param = confSet->Get("NuisanceParameter");
    if(param!=""){
        sys->fNuisanceParameter = Common::CheckName(param);
        TRExFitter::NPMAP[sys->fName] = sys->fNuisanceParameter;
    }
    else{
        sys->fNuisanceParameter = Common::CheckName(sys->fName);
        TRExFitter::NPMAP[sys->fName] = sys->fName;
    }

    if (SystHasProblematicName(sys->fNuisanceParameter)) {
        LOG(ERROR) << "Systematic has a problematic nuisanceparemeter name: " << sys->fNuisanceParameter << "\n";
        ++sc;
    }

    // Set Title
    param = confSet->Get("Title");
    if(param != ""){
        sys->fTitle = Common::RemoveQuotes(param);
        TRExFitter::SYSTMAP[sys->fName] = sys->fTitle;
    }

    // Set TexTitle
    param = confSet->Get("TexTitle");
    if(param!="") TRExFitter::SYSTTEX[sys->fName] = Common::RemoveQuotes(param);

    // attach the syst to the proper samples
    for(auto& isample : fFitter->fSamples) {
        if(isample->fType == Sample::SampleType::DATA) continue;
        if(!isample->fUseSystematics) continue;
        if((samples[0]=="all" || Common::FindInStringVector(samples, isample->fName)>=0 )
           && (exclude[0]=="" || Common::FindInStringVector(exclude, isample->fName)<0 ) ){
            isample->AddSystematic(sys);
        }
    }
    if(Common::FindInStringVector(fFitter->fDecorrSysts,sys->fNuisanceParameter)>=0){
        LOG(INFO) << "Decorrelating systematic with NP = " << sys->fNuisanceParameter << "\n";
        sys->fNuisanceParameter += fFitter->fDecorrSuff;
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystRegionDecorelate(const ConfigSet *confSet,
                                          std::shared_ptr<Systematic> sys,
                                          const std::vector<std::string>& samples,
                                          const std::vector<std::string>& exclude,
                                          const std::vector<std::string>& regions,
                                          int type) {

    int sc(0);
    std::string param = "";

    for(const auto& ireg : fRegNames) {
        bool keepReg=false;
        if ( regions[0]=="all" ) keepReg=true;
        else {
            for ( const std::string& iGoodReg : regions) {
                if ( ireg == iGoodReg ) keepReg=true;
            }
        }
        for ( const std::string& iBadReg : exclude) {
            if ( iBadReg == ireg ) keepReg=false;
        }

        if (!keepReg) {
            LOG(INFO) << "IGNORING REGION: " << ireg << "\n";
            continue;
        }
        LOG(INFO) << "--> KEEPING IT!!! " << ireg << "\n";

        Region* reg = fFitter->GetRegion(ireg);

        if (type == Systematic::STAT) {
            unsigned int nbins = reg->fHistoNBinsRebin>0 ? reg->fHistoNBinsRebin : reg->GetNbins();
            LOG(INFO) << ireg << " " << nbins << "\n";
            // decorrelate by bin
            for (unsigned int i_bin = 0; i_bin < nbins; i_bin++) {
                std::shared_ptr<Systematic> mySys = std::make_shared<Systematic>(*sys);
                mySys->fName=(mySys->fName)+"_"+ireg+"_bin"+std::to_string(i_bin);
                std::vector<std::string> tmpReg;
                tmpReg.emplace_back( ireg );
                mySys->fRegions = tmpReg;
                std::vector<int> tmpBin;
                tmpBin.emplace_back( i_bin );
                mySys->fBins = tmpBin;
                fFitter->fSystematics.emplace_back( mySys );
                mySys->fIsRegionDecorr=true;

                // Set NuisanceParameter
                param = confSet->Get("NuisanceParameter");
                if(param != ""){
                    mySys->fNuisanceParameter = Common::CheckName(param)+"_"+ireg+"_bin"+std::to_string(i_bin);
                    TRExFitter::NPMAP[mySys->fName] = mySys->fNuisanceParameter;
                } else {
                    mySys->fNuisanceParameter = mySys->fName;
                    TRExFitter::NPMAP[mySys->fName] = mySys->fName;
                }

                if (SystHasProblematicName(mySys->fNuisanceParameter)) {
                    LOG(ERROR) << "Systematic has a problematic nuisanceparemeter name: " << mySys->fNuisanceParameter << "\n";
                    ++sc;
                }

                // Set Title
                param = confSet->Get("Title");
                if(param != ""){
                    mySys->fTitle = Common::RemoveQuotes(param)+" ("+reg->fLabel+", bin "+std::to_string(i_bin)+")";
                    TRExFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
                }
                for (auto& isample : fFitter->fSamples) {
                    if(isample->fType == Sample::SampleType::DATA) continue;
                    if(!isample->fUseSystematics) continue;
                    if(   (samples[0]=="all" || Common::FindInStringVector(samples, isample->fName)>=0 )
                       && (exclude[0]==""    || Common::FindInStringVector(exclude, isample->fName)<0 ) ){
                        isample->AddSystematic(mySys);
                    }
                }

            }
        } else {
            //
            // cloning the sys for each region
            std::shared_ptr<Systematic> mySys = std::make_shared<Systematic>(*sys);
            mySys->fName=(mySys->fName)+"_"+ireg;
            std::vector<std::string> tmpReg;
            tmpReg.emplace_back( ireg );
            mySys->fRegions = tmpReg;
            fFitter->fSystematics.emplace_back( mySys );
            mySys->fIsRegionDecorr=true;

            // Set NuisanceParameter
            param = confSet->Get("NuisanceParameter");
            if(param != ""){
                mySys->fNuisanceParameter = Common::CheckName(param)+"_"+ireg;
                TRExFitter::NPMAP[mySys->fName] = mySys->fNuisanceParameter;
            }
            else{
                mySys->fNuisanceParameter = Common::CheckName(mySys->fName);
                TRExFitter::NPMAP[mySys->fName] = Common::CheckName(mySys->fName);
            }

            if (SystHasProblematicName(mySys->fNuisanceParameter)) {
                LOG(ERROR) << "Systematic has a problematic nuisanceparemeter name: " << mySys->fNuisanceParameter << "\n";
                ++sc;
            }

            // Set Title
            param = confSet->Get("Title");
            if(param != ""){
                mySys->fTitle = Common::RemoveQuotes(param)+" ("+reg->fLabel+")";
                TRExFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
            }
            //
            for(auto& isample : fFitter->fSamples) {
                if(isample->fType == Sample::SampleType::DATA) continue;
                if(!isample->fUseSystematics) continue;
                if(   (samples[0]=="all" || Common::FindInStringVector(samples, isample->fName)>=0 )
                    && (exclude[0]==""    || Common::FindInStringVector(exclude, isample->fName)<0 ) ){
                    isample->AddSystematic(mySys);
                }
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystSampleDecorelate(const ConfigSet *confSet, std::shared_ptr<Systematic> sys, const std::vector<std::string> &samples, const std::vector<std::string> &exclude){

    int sc(0);

    std::string param = "";

    // (this is really messy)
    for(auto& isam : fFitter->fSamples) {
        if(isam->fType == Sample::SampleType::DATA) continue;
        if(isam->fType == Sample::SampleType::GHOST || isam->fType == Sample::SampleType::EFT) continue;
        bool keepSam=false;
        if ( samples[0]=="all" ) keepSam=true;
        else {
            if ( std::find(samples.begin(), samples.end(), isam->fName)!=samples.end() ) keepSam=true;
        }
        if ( find(exclude.begin(), exclude.end(), isam->fName)!=exclude.end() ) keepSam=false;
        if (!keepSam) {
            LOG(INFO) << "IGNORING SAMPLE: " << isam->fName << "\n";
            continue;
        }
        LOG(INFO) << "--> KEEPING SAMPLE: " << isam->fName << "\n";
        //
        // cloning the sys for each sample
        std::shared_ptr<Systematic> mySys = std::make_shared<Systematic>(*sys);
        mySys->fName=(mySys->fName)+"_"+isam->fName;
        fFitter->fSystematics.emplace_back( mySys );
        mySys->fIsSampleDecorr=true;

        // Set NuisanceParameter
        param = confSet->Get("NuisanceParameter");
        if(param != ""){
            mySys->fNuisanceParameter = Common::CheckName(param)+"_"+isam->fName;
            TRExFitter::NPMAP[mySys->fName] = mySys->fNuisanceParameter;
        }
        else{
            mySys->fNuisanceParameter = Common::CheckName(mySys->fName);
            TRExFitter::NPMAP[mySys->fName] = Common::CheckName(mySys->fName);
        }

        if (SystHasProblematicName(mySys->fNuisanceParameter)) {
            LOG(ERROR) << "Systematic has a problematic nuisanceparemeter name: " << mySys->fNuisanceParameter << "\n";
            ++sc;
        }

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            mySys->fTitle = Common::RemoveQuotes(param)+" "+isam->fTitle;
            TRExFitter::SYSTMAP[mySys->fName] = mySys->fTitle;
        }

        // Add sample/syst cross-reference
        isam->AddSystematic(mySys);
        mySys->fSamples = { isam->fName };

    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::SetSystShapeDecorelate(const ConfigSet *confSet, std::shared_ptr<Systematic> sys, const std::vector<std::string> &samples, const std::vector<std::string> &exclude){
    int sc(0);

    // cloning the sys
    std::shared_ptr<Systematic> mySys1 = std::make_shared<Systematic>(*sys);
    mySys1->fName=(mySys1->fName)+"_Acc";
    fFitter->fSystematics.emplace_back( mySys1 );
    mySys1->fIsNormOnly = true;
    mySys1->fIsShapeOnly = false;
    mySys1->fIsShapeAccDecorr=true;

    // Set NuisanceParameter
    std::string param = confSet->Get("NuisanceParameter");
    if(param != ""){
        mySys1->fNuisanceParameter = Common::CheckName(param)+"_Acc";
        TRExFitter::NPMAP[mySys1->fName] = mySys1->fNuisanceParameter;
    }
    else{
        mySys1->fNuisanceParameter = Common::CheckName(mySys1->fName);
        TRExFitter::NPMAP[mySys1->fName] = Common::CheckName(mySys1->fName);
    }

    if (SystHasProblematicName(mySys1->fNuisanceParameter)) {
        LOG(ERROR) << "Systematic has a problematic nuisanceparemeter name: " << mySys1->fNuisanceParameter << "\n";
        ++sc;
    }

    // Set Title
    param = confSet->Get("Title");
    if(param != ""){
        mySys1->fTitle = Common::RemoveQuotes(param)+" Acc.";
        TRExFitter::SYSTMAP[mySys1->fName] = mySys1->fTitle;
    }

    for(auto& isam : fFitter->fSamples) {
        if(isam->fType == Sample::SampleType::DATA) continue;
        if(!isam->fUseSystematics) continue;
        if(   (samples[0]=="all" || Common::FindInStringVector(samples, isam->fName)>=0 )
           && (exclude[0]==""    || Common::FindInStringVector(exclude, isam->fName)<0 ) ){
            isam->AddSystematic(mySys1);
        }
    }

    if ( sys->fType!=Systematic::OVERALL ) {
        // cloning the sys
        std::shared_ptr<Systematic> mySys2 = std::make_shared<Systematic>(*sys);
        mySys2->fName=(mySys2->fName)+"_Shape";
        mySys2->fIsNormOnly=false;
        mySys2->fIsShapeOnly=true;
        mySys2->fIsShapeAccDecorr=true;

        fFitter->fSystematics.emplace_back( mySys2 );

        //Set NuisanceParameter
        param = confSet->Get("NuisanceParameter");
        if(param != ""){
            mySys2->fNuisanceParameter = Common::CheckName(param)+"_Shape";
            TRExFitter::NPMAP[mySys2->fName] = mySys2->fNuisanceParameter;
        }
        else{
            mySys2->fNuisanceParameter = Common::CheckName(mySys2->fName);
            TRExFitter::NPMAP[mySys2->fName] = Common::CheckName(mySys2->fName);
        }

        if (SystHasProblematicName(mySys2->fNuisanceParameter)) {
            LOG(ERROR) << "Systematic has a problematic nuisanceparemeter name: " << mySys2->fNuisanceParameter << "\n";
            ++sc;
        }

        // Set Title
        param = confSet->Get("Title");
        if(param != ""){
            mySys2->fTitle = Common::RemoveQuotes(param)+" Shape";
            TRExFitter::SYSTMAP[mySys2->fName] = mySys2->fTitle;
        }

        for(auto& isam : fFitter->fSamples) {
            if(isam->fType == Sample::SampleType::DATA) continue;
            if(!isam->fUseSystematics) continue;
            if(   (samples[0]=="all" || Common::FindInStringVector(samples, isam->fName)>=0 )
               && (exclude[0]==""    || Common::FindInStringVector(exclude, isam->fName)<0 ) ){
                isam->AddSystematic(mySys2);
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::PostConfig(const std::string& opt) {

    int sc(0);

    if (!fHasAtLeastOneValidSample && Common::OptionRunsFit(opt)){
        LOG(ERROR) << "You need to provide at least one sample that is either SIGNAL or BACKGROUND, or provide an UnfoldingSample otherwise the fit will crash.\n";
        ++sc;
    }

    // if StatOnly, also sets to OFF the MC stat
    if(fFitter->fStatOnly){
        LOG(INFO) << "StatOnly option is setting to OFF the MC-stat (gammas) as well.\n";
        LOG(INFO) << "To keep them on use the command line option 'Systematics=NONE'\n";
        LOG(INFO) << "or comment out all Systematics in config file.\n";
        fFitter->SetStatErrorConfig( false, 0. );
    }

    // add nuisance parameter - systematic title correspondence
    for(const auto& syst : fFitter->fSystematics){
        if(syst->fNuisanceParameter!=syst->fName) TRExFitter::SYSTMAP[syst->fNuisanceParameter] = syst->fTitle;
    }
    // add nuisance parameter - norm-factor title correspondence & fix nuisance parameter
    for(const auto& norm : fFitter->fNormFactors){
        if(TRExFitter::NPMAP[norm->fName]=="") TRExFitter::NPMAP[norm->fName] = norm->fName;
        if(norm->fNuisanceParameter!=norm->fName) TRExFitter::SYSTMAP[norm->fNuisanceParameter] = norm->fTitle;
    }

    for (const auto& isample : fFitter->fSamples) {
        if (!SampleIsOk(isample.get())) {
            LOG(ERROR) << "Sample: " << isample->fName << " has the same name as one of the ghost samples or it includes itsef in some of the options\n";
            ++sc;
        }
    }

    sc+= CheckPOIs();

    if (fFitter->fRankingPOIName.empty()) {
        for (std::size_t i = 0; i < fFitter->fPOIs.size(); ++i) {
            fFitter->fRankingPOIName.emplace_back("#mu");
        }
    } else if (fFitter->fRankingPOIName.size() != fFitter->fPOIs.size()) {
        LOG(ERROR) << "Sizes of RankingPOIName and POIs do not match\n";
        ++sc;
    }

    if (fFitter->fRankingUpperAxisNdivision.empty()) {
        for (std::size_t i = 0; i < fFitter->fPOIs.size(); ++i) {
            fFitter->fRankingUpperAxisNdivision.emplace_back(510);
        }
    } else if (fFitter->fRankingUpperAxisNdivision.size() != fFitter->fPOIs.size()) {
        LOG(ERROR) << "Sizes of RankingUpperAxisNdivision and POIs do not match\n";
        ++sc;
    }

    if (fFitter->fRankingPOIAxisScale.empty()) {
        for (std::size_t i = 0; i < fFitter->fPOIs.size(); ++i) {
            fFitter->fRankingPOIAxisScale.emplace_back(1.0);
        }
    } else if (fFitter->fRankingPOIAxisScale.size() != fFitter->fPOIs.size()) {
        LOG(ERROR) << "Sizes of RankingPOIAxisScale and POIs do not match\n";
        ++sc;
    }

    // morphing
    if (fFitter->fMorphParams.size()!=0){
        // template fitting stuff
        fFitter->fTemplateWeightVec = fFitter->GetTemplateWeightVec(fFitter->fTemplateInterpolationOption);
        for(const TRExFit::TemplateWeight& itemp : fFitter->fTemplateWeightVec){
            std::string normName = "morph_"+itemp.name+"_"+Common::ReplaceString(std::to_string(itemp.value),"-","m");
            TRExFitter::SYSTMAP[normName] = itemp.function;
            TRExFitter::NPMAP[normName]   = itemp.name;
            // get the norm factor corresponding to each template
            for(auto norm : fFitter->fNormFactors){
                if(norm->fName == normName){
                    // find a norm factor in the config corresponding to the morphing parameter
                    // NB: it should be there in the config, otherwise an error message is shown and the code crashe
                    bool found = false;
                    for(const auto& norm2 : fFitter->fNormFactors){
                        if(norm2->fName == itemp.name){
                            norm->SetNominalMinMax(norm2->GetNominal(), norm2->GetMin(), norm2->GetMax());
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        LOG(ERROR) << "No NormFactor with name " << itemp.name << " found (needed for morphing\n";
                        LOG(ERROR) << "Please add to the config something like:\n";
                        LOG(ERROR) << "  NormFactor: " << itemp.name << "\n";
                        LOG(ERROR) << "    Min: <the min value for which you have provided template>\n";
                        LOG(ERROR) << "    Max: <the min value for which you have provided template>\n";
                        LOG(ERROR) << "    Samples: none\n";
                        ++sc;
                    }
                }
            }
        }
    }

    // Check and fix shape factor number of bins
    for(const auto& sfactor : fFitter->fShapeFactors){
        LOG(INFO) << "Checking consistency of shape factor " << sfactor->fName << " across regions...\n";
        int nbins = 0;
        for(const auto& regName : sfactor->fRegions){
            const Region *reg = fFitter->GetRegion(regName);
            if(reg==nullptr) continue; // skip non-existing regions
            if(nbins==0) nbins = reg->GetNbins();
            else{
                if(nbins!=reg->GetNbins()){
                    LOG(ERROR) << "Shape factor " << sfactor->fName << " assigned to regions with different number of bins. Please check.\n";
                    ++sc;
                }
            }
        }
        sfactor->fNbins = nbins;
    }

    // set fFitIsBlind, fSignificanceIsBlind and fLimitIsBlind if no sample is DATA
    bool hasData = false;
    for (const auto& smp : fFitter->fSamples) {
        if (smp->fType == Sample::SampleType::DATA) {
            hasData = true;
            break;
        }
    }
    if (!hasData) {
        fFitter->fFitIsBlind = true;
        fFitter->fLimitIsBlind = true;
        fFitter->fSignificanceIsBlind = true;
    }

    // if BlindSRs was set, force to use ASIMOV data in all SRs
    if (fFitter->fBlindSRs) {
        for (auto& reg : fFitter->fRegions) {
            if(reg->fRegionType == Region::SIGNAL) reg->SetRegionDataType(Region::ASIMOVDATA);
        }
    }

    // check all the Expressions, and eventually fix/add range
    for (auto& nf : fFitter->fNormFactors) {
        if (nf->fExpression.first != "") {
            // check the second argument of Expression (which should be the list of dependencies)
            std::string dep = nf->fExpression.second;
            std::string depNew = "";
            // if it contains [ and ], proceed as usual, but check if corresponding NormFactors have cinflicting [nominal,min,max] values
            if (dep.find("[")!=std::string::npos && dep.find("]")!=std::string::npos) {
                // split into arguments according to ']'
                std::vector<std::string> depVec = Common::Vectorize(dep,']');
                for (const auto& s : depVec) {
                    // check that we have a corresponding '['
                    if (Common::Vectorize(s,'[').size()!=2) {
                        LOG(ERROR) << "Problem in Expression for NormFactor " << nf->fName << ": " << nf->fExpression.first << "," << nf->fExpression.second << " (number of [ and ] not matching)\n";
                        ++sc;
                        continue;
                    }
                    // split the parameter part and the range part
                    std::string par = Common::Vectorize(s,'[')[0];
                    par = Common::ReplaceString(par,",",""); // eventually remove comma in front of the parameter name
                    std::string range = Common::Vectorize(s,'[')[1];
                    // check the range part
                    if (Common::Vectorize(range,',').size()!=3) {
                        LOG(ERROR) << "Problem in Expression for NormFactor " << nf->fName << ": " << nf->fExpression.first << "," << nf->fExpression.second << " (must specify nominal, min and max inside [])\n";
                        ++sc;
                        continue;
                    }
                    double nominal = atof(Common::Vectorize(range,',')[0].c_str());
                    double min = atof(Common::Vectorize(range,',')[1].c_str());
                    double max = atof(Common::Vectorize(range,',')[2].c_str());
                    // check if can find corresponding NormFactor and eventually fix these values accordingly
                    for (const auto& nf2 : fFitter->fNormFactors) {
                        if (nf2->fName == par) {
                            if (nf2->GetNominal()!=nominal) {
                                LOG(INFO) << "Fixing nominal value of NormFactor in Expression dependency: " << nominal << " -> " << nf2->GetNominal() << "\n";
                                nominal = nf2->GetNominal();
                            }
                            if (nf2->GetMin()!=min) {
                                LOG(INFO) << "Fixing min value of NormFactor in Expression dependency: " << min << " -> " << nf2->GetMin() << "\n";
                                min = nf2->GetMin();
                            }
                            if (nf2->GetMax()!=max) {
                                LOG(INFO) << "Fixing max value of NormFactor in Expression dependency: " << max << " -> " << nf2->GetMax() << "\n";
                                max = nf2->GetMax();
                            }
                        }
                    }
                    if (depNew!="") depNew += ",";
                    depNew += par + "[" + std::to_string(nominal) + "," + std::to_string(min) + "," + std::to_string(max) + "]";
                }
            }
            // if not, take nominal, min and max from the norm-factors
            else {
                // split into arguments according to ','
                std::vector<std::string> depVec = Common::Vectorize(dep,',');
                for (const auto& par : depVec) {
                    // look for corresponding NormFactor
                    double nominal = 0;
                    double min = 0;
                    double max = 0;
                    bool found = false;
                    for (const auto& nf2 : fFitter->fNormFactors) {
                        if (nf2->fName == par) {
                            nominal = nf2->GetNominal();
                            min = nf2->GetMin();
                            max = nf2->GetMax();
                            found = true;
                        }
                    }
                    if (!found) {
                        LOG(ERROR) << "Problem in Expression for NormFactor " << nf->fName << ": " << nf->fExpression.first << "," << nf->fExpression.second << " (range not specified, and not possible to get it from existing NormFactors)\n";
                        ++sc;
                    }
                    if (depNew!="") depNew += ",";
                    depNew += par + "[" + std::to_string(nominal) + "," + std::to_string(min) + "," + std::to_string(max) + "]";
                }
            }
            // replace with fixed dependencies string
            nf->fExpression.second = depNew;

            // remove the NF from POIs if they exist
            if (std::find(fFitter->fPOIs.begin(), fFitter->fPOIs.end(), nf->fName) != fFitter->fPOIs.end()) {
                fFitter->fPOIs.erase(std::remove(fFitter->fPOIs.begin(), fFitter->fPOIs.end(), nf->fName));
            }
        }
    }

    for (const auto& isf : fFitter->fShapeFactors) {
        if (!isf->fExpression.empty()) {
            fFitter->fShapeFactorReparametrisation = true;
            break;
        }
    }

    if (fHasTemplateMorphing) {
        fFitter->fShapeFactorReparametrisation = true;
        this->ProcessTemplateMorphing();
        this->CheckTemplates();
        if (fTemplateMorphingDimensions == 1) {
            fFitter->fTemplateMorphingFitDimensionality = TemplateMorpher::FIT_DIMENSIONALITY::ONE_DIMENSION;
        } else if (fTemplateMorphingDimensions == 2) {
            fFitter->fTemplateMorphingFitDimensionality = TemplateMorpher::FIT_DIMENSIONALITY::TWO_DIMENSIONS;
        }
    }

    if (fFitter->fFitType == TRExFit::FitType::EFT && !fFitter->fEFTConfig.GetSplitSamplesPerBin()) {
        this->ProcessEFTShapeFactors();
    }

    // propagate fixed NFs
    for (const auto& iConstNF : fConstNFs) {
        auto itr = std::find_if(fFitter->fNormFactors.begin(), fFitter->fNormFactors.end(),
        [&iConstNF](const auto& nf){return iConstNF == nf->fName;});

        if (itr == fFitter->fNormFactors.end()) {
            LOG(ERROR) << "NormFactor " << iConstNF << " asked to be constant does not exist\n";
            ++sc;
            return sc;
        }

        LOG(INFO) << "Fixing NormFactor " << iConstNF << "\n";
        (*itr)->fConst = true;
    }

    this->CheckNormFactors();

    if (fAddNFsToPOIs) {
        for (const auto& inf : fFitter->fNormFactors) {
            const std::string name = inf->fName;
            auto itr = std::find(fFitter->fPOIs.begin(), fFitter->fPOIs.end(), name);
            if (itr == fFitter->fPOIs.end()) {
                LOG(INFO) << "Adding NF: " << name << " to parameters of interest\n";
                fFitter->fPOIs.emplace_back(name);
                fFitter->fRankingPOIName.emplace_back(inf->fTitle);
                fFitter->fRankingUpperAxisNdivision.emplace_back(510);
                fFitter->fRankingPOIAxisScale.emplace_back(1.);
                fFitter->fFitPOIAsimov.insert({name, inf->GetNominal()});
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadUnfoldingOptions() {

    int sc(0);

    const auto* unfolding_confsets = fParser->GetConfigSets("Unfolding");
    if (!unfolding_confsets) {
        if (fFitter->fFitType == TRExFit::UNFOLDING) {
            LOG(ERROR) << "You set FitType == UNFOLDING, but didnt provide Unfolding block!\n";
            ++sc;
        }
        return sc;
    }

    for(const auto& confSetPtr : *unfolding_confsets){
        const ConfigSet *confSet = confSetPtr.get();

        auto unfItr = std::find_if(fFitter->fUnfolding.begin(), fFitter->fUnfolding.end(), [confSet](const std::unique_ptr<Unfolding>& element){return element->fName == Common::RemoveQuotes(confSet->GetValue());});
        if (unfItr != fFitter->fUnfolding.end()) {
            LOG(ERROR) << "Duplicate name found for Unfolding: " << Common::RemoveQuotes(confSet->GetValue()) << "\n";
            ++sc;
            continue;
        }
        fFitter->fUnfolding.emplace_back(std::make_unique<Unfolding>());
        auto& unfolding = fFitter->fUnfolding.back();
        unfolding->fName = Common::RemoveQuotes(confSet->GetValue());

        std::string param = confSet->Get("MatrixOrientation");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "TRUTHONHORIZONTAL") {
                unfolding->fMatrixOrientation = FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS;
            } else if (param == "TRUTHONVERTICAL") {
                unfolding->fMatrixOrientation = FoldingManager::MATRIXORIENTATION::TRUTHONVERTICALAXIS;
            } else {
                LOG(WARNING) << "You specified 'MatrixOrientation' option, but you didnt provide a valid config. Setting to TRUTHONHORIZONTAL.\n";
                unfolding->fMatrixOrientation = FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS;
            }
        }

        param = confSet->Get("TruthDistributionPath");
        if (param != "") {
            unfolding->fTruthDistributionPath = Common::RemoveQuotes(param);
        }

        param = confSet->Get("TruthDistributionFile");
        if (param != "") {
            unfolding->fTruthDistributionFile = Common::RemoveQuotes(param);
        }

        param = confSet->Get("TruthDistributionName");
        if (param != "") {
            unfolding->fTruthDistributionName = Common::RemoveQuotes(param);
        }

        param = confSet->Get("NumberOfTruthBins");
        if (param == "") {
            LOG(ERROR) << "You need to set the number of truth bins!\n";
            ++sc;
        } else {
            const int bins = Common::convertStoNum<int>(param);
            if (bins < 2 || bins > 100) {
                LOG(WARNING) << "You set the number of truth bins which is < 2 or > 100, that looks wrong\n";
            }
            unfolding->fNumberUnfoldingTruthBins = bins;
        }

        param = confSet->Get("Tau");
        if (param != "") {
            const std::vector<std::string>& tmp = Common::Vectorize(param, ',');
            if (tmp.size()==1) {
                if (tmp.at(0).find(":") == std::string::npos) {
                    std::vector<std::pair<int, double> > pair = {std::make_pair(-1, Common::convertStoNum<float>(param))};
                    fTaus.insert({unfolding->fName, pair});
                } else {
                    const std::vector<std::string> oneTau = Common::Vectorize(tmp.at(0), ':');
                    if (oneTau.size() != 2) {
                        LOG(ERROR) << "Wrong format for Tau!\n";
                        ++sc;
                    } else {
                        std::vector<std::pair<int, double> > pair = {std::make_pair(Common::convertStoNum<int>(oneTau.at(0)), Common::convertStoNum<float>(oneTau.at(1)))};
                        fTaus.insert({unfolding->fName, pair});
                    }
                }
            } else {
                for (std::size_t i = 0 ; i < tmp.size(); ++i) {
                    const std::vector<std::string> oneTau = Common::Vectorize(tmp.at(i), ':');
                    if (oneTau.size() != 2) {
                        LOG(ERROR) << "Wrong format for Tau!\n";
                        ++sc;
                    } else {
                        if (i == 0) {
                            std::vector<std::pair<int, double> > pair = {std::make_pair(Common::convertStoNum<int>(oneTau.at(0)), Common::convertStoNum<float>(oneTau.at(1)))};
                            fTaus.insert({unfolding->fName, pair});
                        } else {
                            auto itr = fTaus.find(unfolding->fName);
                            if (itr == fTaus.end()) {
                                LOG(ERROR) << "Something went wrong with the fTau map\n";
                                ++sc;
                            } else {
                                const int bin = Common::convertStoNum<int>(oneTau.at(0));
                                const double value = Common::convertStoNum<float>(oneTau.at(1));
                                itr->second.emplace_back(bin, value);
                            }
                        }
                    }
                }
            }
        }

        param = confSet->Get("UnfoldingResultMin");
        if (param != "") {
            unfolding->fUnfoldingResultMin = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("UnfoldingResultMax");
        if (param != "") {
            unfolding->fUnfoldingResultMax = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("TitleX");
        if (param != "") {
            unfolding->fUnfoldingTitleX = Common::RemoveQuotes(param);
        }

        param = confSet->Get("TitleY");
        if (param != "") {
            unfolding->fUnfoldingTitleY = Common::RemoveQuotes(param);
        }

        param = confSet->Get("ScaleRangeY");
        if (param != "") {
            unfolding->fUnfoldingScaleRangeY = Common::convertStoNum<float>(param);
        }

        param = confSet->Get("MinRangeY");
        if (param != "") {
            unfolding->fUnfoldingMinRangeY = Common::convertStoNum<float>(param);
        }

        param = confSet->Get("TitleOffsetX");
        if (param != "") {
            unfolding->fUnfoldingTitleOffsetX = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("TitleOffsetY");
        if (param != "") {
            unfolding->fUnfoldingTitleOffsetY = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("RatioYmax");
        if (param != "") {
            unfolding->fUnfoldingRatioYmax = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("RatioYmin");
        if (param != "") {
            unfolding->fUnfoldingRatioYmin = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("ErrorBandColor");
        if (param != "") {
            unfolding->fErrorBandColor = Common::convertStoNum<int>(param);
        }

        param = confSet->Get("StatErrorBandColor");
        if (param != "") {
            unfolding->fStatErrorBandColor = Common::convertStoNum<int>(param);
        }

        param = confSet->Get("RatioLineColor");
        if (param != "") {
            unfolding->fRatioLineColor = Common::convertStoNum<int>(param);
        }

        param = confSet->Get("RatioLineStyle");
        if (param != "") {
            unfolding->fRatioLineStyle = Common::convertStoNum<int>(param);
        }

        param = confSet->Get("RatioLineWidth");
        if (param != "") {
            unfolding->fRatioLineWidth = Common::convertStoNum<int>(param);
        }

        param = confSet->Get("LegendXposition");
        if (param != "") {
            unfolding->fLegendXposition = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("LegendTextSize");
        if (param != "") {
            unfolding->fLegendTextSize = Common::convertStoNum<int>(param);
        }

        param = confSet->Get("LogX");
        if (param != "") {
            unfolding->fUnfoldingLogX = Common::StringToBoolean(param);
        }

        param = confSet->Get("LogY");
        if (param != "") {
            unfolding->fUnfoldingLogY = Common::StringToBoolean(param);
        }

        param = confSet->Get("MigrationTitleX");
        if (param != "") {
            const std::vector<std::string> split = Common::Vectorize(param, ',', true, false);
            for (const auto& isplit : split) {
                const std::vector<std::string> oneParam = Common::Vectorize(isplit, ':', true, false);
                if (oneParam.size() == 1 && split.size() == 1) {
                    unfolding->fMigrationTitleX.emplace_back(std::make_pair("all", oneParam.at(0)));
                } else if (oneParam.size() == 2) {
                    unfolding->fMigrationTitleX.emplace_back(std::make_pair(oneParam.at(0), oneParam.at(1)));
                } else {
                    LOG(ERROR) << "Invalid option for \"MigrationTitleX\"\n";
                    ++sc;
                }
            }
        }

        param = confSet->Get("MigrationTitleY");
        if (param != "") {
            const std::vector<std::string> split = Common::Vectorize(param, ',', true, false);
            for (const auto& isplit : split) {
                const std::vector<std::string> oneParam = Common::Vectorize(isplit, ':', true, false);
                if (oneParam.size() == 1 && split.size() == 1) {
                    unfolding->fMigrationTitleY.emplace_back(std::make_pair("all", oneParam.at(0)));
                } else if (oneParam.size() == 2) {
                    unfolding->fMigrationTitleY.emplace_back(std::make_pair(oneParam.at(0), oneParam.at(1)));
                } else {
                    LOG(ERROR) << "Invalid option for \"MigrationTitleY\"\n";
                    ++sc;
                }
            }
        }

        param = confSet->Get("MigrationLogX");
        if (param != "") {
            unfolding->fMigrationLogX = Common::StringToBoolean(param);
        }

        param = confSet->Get("MigrationLogY");
        if (param != "") {
            unfolding->fMigrationLogY = Common::StringToBoolean(param);
        }

        param = confSet->Get("MigrationTitleOffsetX");
        if (param != "") {
            unfolding->fMigrationTitleOffsetX = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("MigrationTitleOffsetY");
        if (param != "") {
            unfolding->fMigrationTitleOffsetY = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("MigrationZmin");
        if (param != "") {
            unfolding->fMigrationZmin = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("MigrationZmax");
        if (param != "") {
            unfolding->fMigrationZmax = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("ResponseZmin");
        if (param != "") {
            unfolding->fResponseZmin = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("ResponseZmax");
        if (param != "") {
            unfolding->fResponseZmax = Common::convertStoNum<double>(param);
        }

        param = confSet->Get("PlotSystematicMigrations");
        if (param != "") {
            unfolding->fPlotSystematicMigrations = Common::StringToBoolean(param);
        }

        param = confSet->Get("MigrationText");
        if (param != "") {
            unfolding->fMigrationText = Common::StringToBoolean(param);
        }

        param = confSet->Get("NominalTruthSample");
        if (param == "") {
            LOG(ERROR) << "You need to set NominalTruthSample option!\n";
            ++sc;
        } else {
            unfolding->fNominalTruthSample = Common::RemoveQuotes(param);
        }

        param = confSet->Get("AlternativeAsimovTruthSample");
        if (param != "") {
            unfolding->fAlternativeAsimovTruthSample = Common::RemoveQuotes(param);
            fFitter->fFitIsBlind = false;
        }

        param = confSet->Get("UnfoldNormXSec");
        if (param != "") {
            unfolding->fUnfoldNormXSec = Common::StringToBoolean(param);
        }

        param = confSet->Get("UnfoldNormXSecBinN");
        if (param != "") {
            unfolding->fUnfoldNormXSecBinN = Common::convertStoNum<int>(param);
        }

        param = confSet->Get("DivideByBinWidth");
        if (param != "") {
            unfolding->fUnfoldingDivideByBinWidth = Common::StringToBoolean(param);
        }

        param = confSet->Get("UncertaintyBreakdownRange");
        if (param != "") {
            unfolding->fBreakdownRange = std::stod(param);
        }

        param = confSet->Get("DivideByLumi");
        if (param != "") {
            const double value = Common::convertStoNum<float>(param);
            if (value < 0) {
                LOG(WARNING) << "DivideByLumi < 0, did you set it properly?\n";
            }
            unfolding->fUnfoldingDivideByLumi = value;
        }

        param = confSet->Get("Expressions");
        if (param != "") {
            const std::vector<std::string>& tmp = Common::Vectorize(param, ',');
            for (std::size_t i = 0; i < tmp.size(); ++i) {
                const std::vector<std::string>& oneExpr = Common::Vectorize(tmp.at(i), '=');
                if (oneExpr.size() != 2) {
                    LOG(ERROR) << "Wrong format for Expressions!\n";
                    ++sc;
                }
                else {
                    if (i == 0) {
                        std::vector<std::pair<std::string, std::string> > vec = {std::make_pair(oneExpr.at(0),oneExpr.at(1))};
                        fExpressions.insert({unfolding->fName, vec});
                    } else {
                        auto itr = fExpressions.find(unfolding->fName);
                        if (itr == fExpressions.end()) {
                            LOG(ERROR) << "Something went wrong with the expressions map\n";
                            ++sc;
                        } else {
                            itr->second.emplace_back(oneExpr.at(0),oneExpr.at(1));
                        }
                    }
                }
            }
        }

        param = confSet->Get("UnfoldingChi2Type");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "EMPTY") {
                unfolding->fUnfoldingChi2Type = Unfolding::UnfoldingChi2Type::EMPTY;
            } else if (param == "CHI2") {
                unfolding->fUnfoldingChi2Type = Unfolding::UnfoldingChi2Type::CHI2;
            } else if (param == "CHI2NDF") {
                unfolding->fUnfoldingChi2Type = Unfolding::UnfoldingChi2Type::CHI2NDF;
            } else if (param == "PROB") {
                unfolding->fUnfoldingChi2Type = Unfolding::UnfoldingChi2Type::PROB;
            }
        }

    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadTruthSamples() {

    int sc(0);

    std::vector<bool> found(fFitter->fUnfolding.size(), false);

    std::vector<std::string> names;
    const auto* truthSample_confsets = fParser->GetConfigSets("TruthSample");
    if (!truthSample_confsets) return sc;

    for(const auto& confSetPtr : *truthSample_confsets){
        const ConfigSet *confSet = confSetPtr.get();
        auto sample = std::make_unique<TruthSample>(Common::RemoveQuotes(confSet->GetValue()));
        for (std::size_t i = 0; i < fFitter->fUnfolding.size(); ++i) {
            if (sample->GetName() == fFitter->fUnfolding.at(i)->fNominalTruthSample) {
                found.at(i) = true;
            }
        }

        if (std::find(names.begin(), names.end(), Common::RemoveQuotes(confSet->GetValue())) == names.end()) {
            names.emplace_back(Common::RemoveQuotes(confSet->GetValue()));
        } else {
            LOG(ERROR) << "Multiply defined TruthSample: " << Common::RemoveQuotes(confSet->GetValue()) << "\n";
            ++sc;
        }

        std::string param = confSet->Get("Title");
        if (param != "") {
            sample->SetTitle(Common::RemoveQuotes(param));
        }

        param = confSet->Get("LineStyle");
        if (param != "") {
            sample->SetLineStyle(Common::convertStoNum<int>(param));
        }

        param = confSet->Get("LineColor");
        if (param != "") {
            sample->SetLineColor(Common::convertStoNum<int>(param));
        }

        param = confSet->Get("LineColorRGB");
        if (param != "") {
            const std::vector<std::string> colours = Common::Vectorize(param, ',');
            sample->SetLineColor(Common::ColorFromRGB(colours));
        }

        param = confSet->Get("LineWidth");
        if (param != "") {
            const int width = Common::convertStoNum<int>(param);
            if (width < 1) {
                LOG(WARNING) << "\"LineWidth\" is set to be < 1, setting to 2\n";
                sample->SetLineWidth(2);
            } else {
                sample->SetLineWidth(width);
            }
        }

        param = confSet->Get("UseForPlotting");
        if (param != "") {
            sample->SetUseForPlotting(Common::StringToBoolean(param));
        }

        param = confSet->Get("TruthDistributionPath");
        if (param != "") {
            sample->SetTruthDistributionPath(Common::RemoveQuotes(param));
        }

        param = confSet->Get("TruthDistributionFile");
        if (param != "") {
            sample->SetTruthDistributionFile(Common::RemoveQuotes(param));
        }

        param = confSet->Get("TruthDistributionName");
        if (param != "") {
            sample->SetTruthDistributionName(Common::RemoveQuotes(param));
        }

        if (fFitter->fUnfolding.size() == 1) {
            sample->fUnfoldingName = fFitter->fUnfolding.at(0)->fName;
        } else {
            param = confSet->Get("UnfoldingName");
            if (param == "") {
                LOG(ERROR) << "Multiple unfoldings but no UnfoldingName provided for a TruthSample\n";
                ++sc;
            } else {
                param = Common::RemoveQuotes(param);
                sample->fUnfoldingName = param;
                if (std::find_if(fFitter->fUnfolding.begin(), fFitter->fUnfolding.end(), [&param](const std::unique_ptr<Unfolding>& unf){return unf->fName == param;})
                    == fFitter->fUnfolding.end()) {
                    LOG(ERROR) << "UnfoldingName: " << param << " does not match any Unfolding names\n";
                    ++sc;
                }
            }
        }

        fFitter->fTruthSamples.emplace_back(std::move(sample));
    }

    if (!fFitter->fTruthSamples.empty()) {
        for (std::size_t i = 0; i < found.size(); ++i) {
            if (!found.at(i)) {
                LOG(ERROR) << "The NominalTruthSample for Unfolding " << fFitter->fUnfolding.at(i)->fName << " not found in any of the TruthSample\n";
                ++sc;
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadUnfoldingSampleOptions() {

    int sc(0);
    const auto* unfoldingSample_confsets = fParser->GetConfigSets("UnfoldingSample");
    if (!unfoldingSample_confsets) return sc;

    for(const auto& confSetPtr : *unfoldingSample_confsets){
        const ConfigSet *confSet = confSetPtr.get();
        if (fFitter->fUnfolding.empty()) {
            LOG(ERROR) << "Need to define Unfolding blocks before the UnfoldingSample blocks\n";
            ++sc;
            return sc;
        }

        fHasAtLeastOneValidSample = true;

        if(fOnlySamples.size()>0 && Common::FindInStringVector(fOnlySamples,Common::RemoveQuotes(confSet->GetValue()))<0) continue;
        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,Common::RemoveQuotes(confSet->GetValue()))>=0) continue;

        auto sample = std::make_unique<UnfoldingSample>();
        sample->SetName(Common::RemoveQuotes(confSet->GetValue()));

        std::string param = confSet->Get("Title");
        if (param != "") {
            sample->SetTitle(Common::RemoveQuotes(param));
        }

        param = confSet->Get("Type");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "STANDARD") {
                sample->SetType(UnfoldingSample::TYPE::STANDARD);
            } else if (param == "GHOST") {
                sample->SetType(UnfoldingSample::TYPE::GHOST);
            }
        }

        param = confSet->Get("ResponseMatrixFile");
        if (param != "") {
            sample->fResponseMatrixFiles.clear();
            sample->fResponseMatrixFiles.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixFiles");
        if (param != "") {
            sample->fResponseMatrixFiles = Common::Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixName");
        if (param != "") {
            sample->fResponseMatrixNames.clear();
            sample->fResponseMatrixNames.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixNames");
        if (param != "") {
            sample->fResponseMatrixNames = Common::Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPath");
        if (param != "") {
            sample->fResponseMatrixPaths.clear();
            sample->fResponseMatrixPaths.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPaths");
        if (param != "") {
            sample->fResponseMatrixPaths = Common::Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixFileSuff");
        if (param != "") {
            sample->fResponseMatrixFileSuffs.clear();
            sample->fResponseMatrixFileSuffs.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixFileSuffs");
        if (param != "") {
            sample->fResponseMatrixFileSuffs = Common::Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixNameSuff");
        if (param != "") {
            sample->fResponseMatrixNameSuffs.clear();
            sample->fResponseMatrixNameSuffs.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixNameSuffs");
        if (param != "") {
            sample->fResponseMatrixNameSuffs = Common::Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPathSuff");
        if (param != "") {
            sample->fResponseMatrixPathSuffs.clear();
            sample->fResponseMatrixPathSuffs.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasResponse(true);
        }

        param = confSet->Get("ResponseMatrixPathSuffs");
        if (param != "") {
            sample->fResponseMatrixPathSuffs = Common::Vectorize(param, ',');
            sample->SetHasResponse(true);
        }

        param = confSet->Get("AcceptanceFile");
        if (param != "") {
            sample->fAcceptanceFiles.clear();
            sample->fAcceptanceFiles.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFiles");
        if (param != "") {
            sample->fAcceptanceFiles = Common::Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceName");
        if (param != "") {
            sample->fAcceptanceNames.clear();
            sample->fAcceptanceNames.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNames");
        if (param != "") {
            sample->fAcceptanceNames = Common::Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePath");
        if (param != "") {
            sample->fAcceptancePaths.clear();
            sample->fAcceptancePaths.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePaths");
        if (param != "") {
            sample->fAcceptancePaths = Common::Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuff");
        if (param != "") {
            sample->fAcceptanceFileSuffs.clear();
            sample->fAcceptanceFileSuffs.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffs");
        if (param != "") {
            sample->fAcceptanceFileSuffs = Common::Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuff");
        if (param != "") {
            sample->fAcceptanceNameSuffs.clear();
            sample->fAcceptanceNameSuffs.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffs");
        if (param != "") {
            sample->fAcceptanceNameSuffs = Common::Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuff");
        if (param != "") {
            sample->fAcceptancePathSuffs.clear();
            sample->fAcceptancePathSuffs.emplace_back(Common::RemoveQuotes(param));
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffs");
        if (param != "") {
            sample->fAcceptancePathSuffs = Common::Vectorize(param, ',');
            sample->SetHasAcceptance(true);
        }

        param = confSet->Get("SelectionEffFile");
        if (param != "") {
            sample->fSelectionEffFiles.clear();
            sample->fSelectionEffFiles.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffFiles");
        if (param != "") {
            sample->fSelectionEffFiles = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffName");
        if (param != "") {
            sample->fSelectionEffNames.clear();
            sample->fSelectionEffNames.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNames");
        if (param != "") {
            sample->fSelectionEffNames = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffPath");
        if (param != "") {
            sample->fSelectionEffPaths.clear();
            sample->fSelectionEffPaths.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPaths");
        if (param != "") {
            sample->fSelectionEffPaths = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFileSuff");
        if (param != "") {
            sample->fSelectionEffFileSuffs.clear();
            sample->fSelectionEffFileSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffFileSuffs");
        if (param != "") {
            sample->fSelectionEffFileSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffNameSuff");
        if (param != "") {
            sample->fSelectionEffNameSuffs.clear();
            sample->fSelectionEffNameSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNameSuffs");
        if (param != "") {
            sample->fSelectionEffNameSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffPathSuff");
        if (param != "") {
            sample->fSelectionEffPathSuffs.clear();
            sample->fSelectionEffPathSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPathSuffs");
        if (param != "") {
            sample->fSelectionEffPathSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFile");
        if (param != "") {
            sample->fMigrationFiles.clear();
            sample->fMigrationFiles.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationFiles");
        if (param != "") {
            sample->fMigrationFiles = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationName");
        if (param != "") {
            sample->fMigrationNames.clear();
            sample->fMigrationNames.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNames");
        if (param != "") {
            sample->fMigrationNames = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationPath");
        if (param != "") {
            sample->fMigrationPaths.clear();
            sample->fMigrationPaths.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPaths");
        if (param != "") {
            sample->fMigrationPaths = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFileSuff");
        if (param != "") {
            sample->fMigrationFileSuffs.clear();
            sample->fMigrationFileSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationFileSuffs");
        if (param != "") {
            sample->fMigrationFileSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationNameSuff");
        if (param != "") {
            sample->fMigrationNameSuffs.clear();
            sample->fMigrationNameSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNameSuffs");
        if (param != "") {
            sample->fMigrationNameSuffs = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationPathSuff");
        if (param != "") {
            sample->fMigrationPathSuffs.clear();
            sample->fMigrationPathSuffs.emplace_back(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPathSuffs");
        if (param != "") {
            sample->fMigrationPathSuffs = Common::Vectorize(param, ',');
        }

        // Set FillColor
        param = confSet->Get("FillColor");
        if(param != "") sample->SetFillColor(std::atoi(param.c_str()));

        param = confSet->Get("FillColorRGB");
        if(param != "") {
            const std::vector<std::string> colours = Common::Vectorize(param, ',');
            sample->SetFillColor(Common::ColorFromRGB(colours));
        }

        // Set LineColor
        param = confSet->Get("LineColor");
        if(param != "") sample->SetLineColor(std::atoi(param.c_str()));

        param = confSet->Get("LineColorRGB");
        if(param != "") {
            const std::vector<std::string> colours = Common::Vectorize(param, ',');
            sample->SetLineColor(Common::ColorFromRGB(colours));
        }

        // Set Regions
        param = confSet->Get("Regions");
        if(param == "") {
            sample->fRegions.emplace_back("all");
        } else {
            sample->fRegions = Common::Vectorize(param, ',');
        }

        // Set Regions for individual sub-samples (truth bins)
        param = confSet->Get("SubSampleRegions");
        if(param != "") {
            std::vector<std::string> tmpVec = Common::Vectorize(param, ',');
            for (const auto& s : tmpVec) {
                int bin = atoi(Common::Vectorize(s, ':')[0].c_str());
                std::string region = Common::Vectorize(s, ':')[1];
                sample->fSubSampleRegions[bin].emplace_back(region);
            }
        }

        // Set Gammas
        param = confSet->Get("Gammas");
        if(param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "SEPARATED") {
                sample->SetGammas(UnfoldingSample::GAMMAS::SEPARATED);
            } else if (param == "DISABLED") {
                sample->SetGammas(UnfoldingSample::GAMMAS::DISABLED);
            }
        }

        if (fFitter->fUnfolding.size() == 1) {
            sample->fUnfoldingName = fFitter->fUnfolding.at(0)->fName;
        } else {
            param = confSet->Get("UnfoldingName");
            if (param == "") {
                LOG(ERROR) << "Multiple unfoldings but no UnfoldingName provided for an UnfoldingSample\n";
                ++sc;
            } else {
                param = Common::RemoveQuotes(param);
                sample->fUnfoldingName = param;
                if (std::find_if(fFitter->fUnfolding.begin(), fFitter->fUnfolding.end(), [&param](const std::unique_ptr<Unfolding>& unf){return unf->fName == param;})
                    == fFitter->fUnfolding.end()) {
                    LOG(ERROR) << "UnfoldingName: " << param << " does not match any Unfolding names\n";
                    ++sc;
                }
            }
        }

        // Set NoPruning
        param = confSet->Get("NoPruning");
        if(param != ""){
            sample->SetNoPruning(Common::StringToBoolean(param));
        }

        fFitter->fUnfoldingSamples.emplace_back(std::move(sample));

    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadUnfoldingSystematicOptions() {

    int sc(0);

    const auto* unfoldingsyst_confsets = fParser->GetConfigSets("UnfoldingSystematic");
    if (unfoldingsyst_confsets == nullptr) return sc;

    for(const auto& confSetPtr : *unfoldingsyst_confsets) {
        const ConfigSet* confSet = confSetPtr.get();
        if(fOnlySystematics.size()>0 && Common::FindInStringVector(fOnlySystematics,Common::CheckName(confSet->GetValue()))<0) continue;
        if(fToExclude.size()>0 && Common::FindInStringVector(fToExclude,Common::CheckName(confSet->GetValue()))>=0) continue;

        auto syst = std::make_unique<UnfoldingSystematic>();
        syst->SetName(Common::RemoveQuotes(confSet->GetValue()));

        // Set NuisanceParameter
        std::string param = confSet->Get("NuisanceParameter");
        if(param != ""){
            syst->fNuisanceParameter = Common::CheckName(param);
            TRExFitter::NPMAP[syst->GetName()] = syst->fNuisanceParameter;
        } else {
            syst->fNuisanceParameter = Common::CheckName(syst->GetName());
            TRExFitter::NPMAP[syst->GetName()] = syst->GetName();
        }

        param = confSet->Get("Type");
        if (param == "") {
            syst->SetType(Systematic::HISTO);
        } else {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);

            if(param == "OVERALL") syst->SetType(Systematic::OVERALL);
            else if(param == "SHAPE") syst->SetType(Systematic::SHAPE);
            else if (param == "STAT") syst->SetType(Systematic::STAT);
            else if (param == "HISTO") syst->SetType(Systematic::HISTO);
        }

        param = confSet->Get("OverallUp");
        if (param != "") {
            if (syst->GetType() != Systematic::SystType::OVERALL) {
                LOG(ERROR) << "Systematic: " << confSet->GetValue() << " has OverallUp set, but the type is not OVERALL!\n";
                ++sc;
            } else {
                syst->SetOverallUp(Common::convertStoNum<float>(param));
            }
        }
        param = confSet->Get("OverallDown");
        if (param != "") {
            if (syst->GetType() != Systematic::SystType::OVERALL) {
                LOG(ERROR) << "Systematic: " << confSet->GetValue() << " has OverallDown set, but the type is not OVERALL!\n";
                ++sc;
            } else {
                syst->SetOverallDown(Common::convertStoNum<float>(param));
            }
        }

        param = confSet->Get("Samples");
        if (param == "") {
            syst->fSamples.emplace_back("all");
        } else {
            syst->fSamples = Common::Vectorize(param, ',');
        }

        param = confSet->Get("Regions");
        if (param == "") {
            syst->fRegions.emplace_back("all");
        } else {
            syst->fRegions = Common::Vectorize(param, ',');
        }

        param = confSet->Get("Title");
        if (param == "") {
            syst->SetTitle(Common::RemoveQuotes(confSet->GetValue()));
        } else {
            syst->SetTitle(Common::RemoveQuotes(param));
        }

        // SetCategory
        param = confSet->Get("Category");
        if(param != ""){
            syst->fCategory = Common::RemoveQuotes(param);
            syst->fSubCategory = Common::RemoveQuotes(param); //SubCategory defaults to the Category setting, if the Category is explicitly set
        }

        // SetSubCategory
        param = confSet->Get("SubCategory");
        if (param != ""){
            syst->fSubCategory = Common::RemoveQuotes(param); // note this needs to happen after Category was set, in order to overwrite the default if required
        }

        // CombineName
        param = confSet->Get("CombineName");
        if(param != ""){
            syst->fCombineName = Common::RemoveQuotes(param);
        }

        // CombineType
        param = confSet->Get("CombineType");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            const std::string tmp = Common::RemoveQuotes(param);

            if (tmp == "STANDARDDEVIATION") {
                syst->fCombineType = Systematic::COMBINATIONTYPE::STANDARDDEVIATION;
            } else if (tmp == "ENVELOPE") {
                syst->fCombineType = Systematic::COMBINATIONTYPE::ENVELOPE;
            } else if (tmp == "STANDARDDEVIATIONNODDOF") {
                syst->fCombineType = Systematic::COMBINATIONTYPE::STANDARDDEVIATIONNODDOF;
            } else if (tmp == "HESSIAN") {
                syst->fCombineType = Systematic::COMBINATIONTYPE::HESSIAN;
            } else if (tmp == "SUMINSQUARES") {
                syst->fCombineType = Systematic::COMBINATIONTYPE::SUMINSQUARES;
            } else {
                LOG(WARNING) << "You specified 'CombineType' option but did not provide valid parameter. Using default (ENVELOPE)\n";
                syst->fCombineType = Systematic::COMBINATIONTYPE::ENVELOPE;
            }
        }

        bool hasUp(false);
        bool hasDown(false);

        param = confSet->Get("ResponseMatrixPathUp");
        if (param != "") {
            syst->fResponseMatrixPathsUp.clear();
            syst->fResponseMatrixPathsUp.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathsUp");
        if (param != "") {
            syst->fResponseMatrixPathsUp = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathDown");
        if (param != "") {
            syst->fResponseMatrixPathsDown.clear();
            syst->fResponseMatrixPathsDown.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixPathsDown");
        if (param != "") {
            syst->fResponseMatrixPathsDown = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffUp");
        if (param != "") {
            syst->fResponseMatrixPathSuffsUp.clear();
            syst->fResponseMatrixPathSuffsUp.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffsUp");
        if (param != "") {
            syst->fResponseMatrixPathSuffsUp = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffDown");
        if (param != "") {
            syst->fResponseMatrixPathSuffsDown.clear();
            syst->fResponseMatrixPathSuffsDown.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixPathSuffsDown");
        if (param != "") {
            syst->fResponseMatrixPathSuffsDown = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNameUp");
        if (param != "") {
            syst->fResponseMatrixNamesUp.clear();
            syst->fResponseMatrixNamesUp.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNamesUp");
        if (param != "") {
            syst->fResponseMatrixNamesUp = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFolderNameUp");
        if (param != "") {
            syst->fResponseMatrixFolderNamesUp.clear();
            syst->fResponseMatrixFolderNamesUp.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNameDown");
        if (param != "") {
            syst->fResponseMatrixNamesDown.clear();
            syst->fResponseMatrixNamesDown.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNamesDown");
        if (param != "") {
            syst->fResponseMatrixNamesDown = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFolderNameDown");
        if (param != "") {
            syst->fResponseMatrixFolderNamesDown.clear();
            syst->fResponseMatrixFolderNamesDown.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffUp");
        if (param != "") {
            syst->fResponseMatrixNameSuffsUp.clear();
            syst->fResponseMatrixNameSuffsUp.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffsUp");
        if (param != "") {
            syst->fResponseMatrixNameSuffsUp = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffDown");
        if (param != "") {
            syst->fResponseMatrixNameSuffsDown.clear();
            syst->fResponseMatrixNameSuffsDown.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixNameSuffsDown");
        if (param != "") {
            syst->fResponseMatrixNameSuffsDown = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFileUp");
        if (param != "") {
            syst->fResponseMatrixFilesUp.clear();
            syst->fResponseMatrixFilesUp.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFilesUp");
        if (param != "") {
            syst->fResponseMatrixFilesUp = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFileDown");
        if (param != "") {
            syst->fResponseMatrixFilesDown.clear();
            syst->fResponseMatrixFilesDown.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFilesDown");
        if (param != "") {
            syst->fResponseMatrixFilesDown = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffUp");
        if (param != "") {
            syst->fResponseMatrixFileSuffsUp.clear();
            syst->fResponseMatrixFileSuffsUp.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffsUp");
        if (param != "") {
            syst->fResponseMatrixFileSuffsUp = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasUp = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffDown");
        if (param != "") {
            syst->fResponseMatrixFileSuffsDown.clear();
            syst->fResponseMatrixFileSuffsDown.emplace_back(Common::RemoveQuotes(param));
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixFileSuffsDown");
        if (param != "") {
            syst->fResponseMatrixFileSuffsDown = Common::Vectorize(param,',');
            syst->SetHasResponse(true);
            hasDown = true;
        }

        param = confSet->Get("AcceptancePathUp");
        if (param != "") {
            syst->fAcceptancePathsUp.clear();
            syst->fAcceptancePathsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathsUp");
        if (param != "") {
            syst->fAcceptancePathsUp = Common::Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathDown");
        if (param != "") {
            syst->fAcceptancePathsDown.clear();
            syst->fAcceptancePathsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathsDown");
        if (param != "") {
            syst->fAcceptancePathsDown = Common::Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffUp");
        if (param != "") {
            syst->fAcceptancePathSuffsUp.clear();
            syst->fAcceptancePathSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffsUp");
        if (param != "") {
            syst->fAcceptancePathSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffDown");
        if (param != "") {
            syst->fAcceptancePathSuffsDown.clear();
            syst->fAcceptancePathSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptancePathSuffsDown");
        if (param != "") {
            syst->fAcceptancePathSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameUp");
        if (param != "") {
            syst->fAcceptanceNamesUp.clear();
            syst->fAcceptanceNamesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNamesUp");
        if (param != "") {
            syst->fAcceptanceNamesUp = Common::Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFolderNameUp");
        if (param != "") {
            syst->fAcceptanceFolderNamesUp.clear();
            syst->fAcceptanceFolderNamesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameDown");
        if (param != "") {
            syst->fAcceptanceNamesDown.clear();
            syst->fAcceptanceNamesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNamesDown");
        if (param != "") {
            syst->fAcceptanceNamesDown = Common::Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFolderNameDown");
        if (param != "") {
            syst->fAcceptanceFolderNamesDown.clear();
            syst->fAcceptanceFolderNamesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffUp");
        if (param != "") {
            syst->fAcceptanceNameSuffsUp.clear();
            syst->fAcceptanceNameSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffsUp");
        if (param != "") {
            syst->fAcceptanceNameSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffDown");
        if (param != "") {
            syst->fAcceptanceNameSuffsDown.clear();
            syst->fAcceptanceNameSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceNameSuffsDown");
        if (param != "") {
            syst->fAcceptanceNameSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileUp");
        if (param != "") {
            syst->fAcceptanceFilesUp.clear();
            syst->fAcceptanceFilesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFilesUp");
        if (param != "") {
            syst->fAcceptanceFilesUp = Common::Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileDown");
        if (param != "") {
            syst->fAcceptanceFilesDown.clear();
            syst->fAcceptanceFilesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFilesDown");
        if (param != "") {
            syst->fAcceptanceFilesDown = Common::Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffUp");
        if (param != "") {
            syst->fAcceptanceFileSuffsUp.clear();
            syst->fAcceptanceFileSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffsUp");
        if (param != "") {
            syst->fAcceptanceFileSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffDown");
        if (param != "") {
            syst->fAcceptanceFileSuffsDown.clear();
            syst->fAcceptanceFileSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("AcceptanceFileSuffsDown");
        if (param != "") {
            syst->fAcceptanceFileSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
            syst->SetHasAcceptance(true);
        }
        param = confSet->Get("SelectionEffPathUp");
        if (param != "") {
            syst->fSelectionEffPathsUp.clear();
            syst->fSelectionEffPathsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("SelectionEffPathsUp");
        if (param != "") {
            syst->fSelectionEffPathsUp = Common::Vectorize(param,',');
            hasUp = true;
            syst->SetHasAcceptance(true);
        }

        param = confSet->Get("SelectionEffPathDown");
        if (param != "") {
            syst->fSelectionEffPathsDown.clear();
            syst->fSelectionEffPathsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffPathsDown");
        if (param != "") {
            syst->fSelectionEffPathsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffPathSuffUp");
        if (param != "") {
            syst->fSelectionEffPathSuffsUp.clear();
            syst->fSelectionEffPathSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffPathSuffsUp");
        if (param != "") {
            syst->fSelectionEffPathSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffPathSuffDown");
        if (param != "") {
            syst->fSelectionEffPathSuffsDown.clear();
            syst->fSelectionEffPathSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffPathSuffsDown");
        if (param != "") {
            syst->fSelectionEffPathSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNameUp");
        if (param != "") {
            syst->fSelectionEffNamesUp.clear();
            syst->fSelectionEffNamesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNamesUp");
        if (param != "") {
            syst->fSelectionEffNamesUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFolderNameUp");
        if (param != "") {
            syst->fSelectionEffFolderNamesUp.clear();
            syst->fSelectionEffFolderNamesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNameDown");
        if (param != "") {
            syst->fSelectionEffNamesDown.clear();
            syst->fSelectionEffNamesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNamesDown");
        if (param != "") {
            syst->fSelectionEffNamesDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFolderNameDown");
        if (param != "") {
            syst->fSelectionEffFolderNamesDown.clear();
            syst->fSelectionEffFolderNamesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNameSuffUp");
        if (param != "") {
            syst->fSelectionEffNameSuffsUp.clear();
            syst->fSelectionEffNameSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNameSuffsUp");
        if (param != "") {
            syst->fSelectionEffNameSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffNameSuffDown");
        if (param != "") {
            syst->fSelectionEffNameSuffsDown.clear();
            syst->fSelectionEffNameSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffNameSuffsDown");
        if (param != "") {
            syst->fSelectionEffNameSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFileUp");
        if (param != "") {
            syst->fSelectionEffFilesUp.clear();
            syst->fSelectionEffFilesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFilesUp");
        if (param != "") {
            syst->fSelectionEffFilesUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFileDown");
        if (param != "") {
            syst->fSelectionEffFilesDown.clear();
            syst->fSelectionEffFilesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFilesDown");
        if (param != "") {
            syst->fSelectionEffFilesDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFileSuffUp");
        if (param != "") {
            syst->fSelectionEffFileSuffsUp.clear();
            syst->fSelectionEffFileSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFileSuffsUp");
        if (param != "") {
            syst->fSelectionEffFileSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("SelectionEffFileSuffDown");
        if (param != "") {
            syst->fSelectionEffFileSuffsDown.clear();
            syst->fSelectionEffFileSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("SelectionEffFileSuffsDown");
        if (param != "") {
            syst->fSelectionEffFileSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationPathUp");
        if (param != "") {
            syst->fMigrationPathsUp.clear();
            syst->fMigrationPathsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationPathsUp");
        if (param != "") {
            syst->fMigrationPathsUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationPathDown");
        if (param != "") {
            syst->fMigrationPathsDown.clear();
            syst->fMigrationPathsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationPathsDown");
        if (param != "") {
            syst->fMigrationPathsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationPathSuffUp");
        if (param != "") {
            syst->fMigrationPathSuffsUp.clear();
            syst->fMigrationPathSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationPathSuffsUp");
        if (param != "") {
            syst->fMigrationPathSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationPathSuffDown");
        if (param != "") {
            syst->fMigrationPathSuffsDown.clear();
            syst->fMigrationPathSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationPathSuffsDown");
        if (param != "") {
            syst->fMigrationPathSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationNameUp");
        if (param != "") {
            syst->fMigrationNamesUp.clear();
            syst->fMigrationNamesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationNamesUp");
        if (param != "") {
            syst->fMigrationNamesUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationFolderNameUp");
        if (param != "") {
            syst->fMigrationFolderNamesUp.clear();
            syst->fMigrationFolderNamesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationNameDown");
        if (param != "") {
            syst->fMigrationNamesDown.clear();
            syst->fMigrationNamesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationNamesDown");
        if (param != "") {
            syst->fMigrationNamesDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationFolderNameDown");
        if (param != "") {
            syst->fMigrationFolderNamesDown.clear();
            syst->fMigrationFolderNamesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationNameSuffUp");
        if (param != "") {
            syst->fMigrationNameSuffsUp.clear();
            syst->fMigrationNameSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationNameSuffsUp");
        if (param != "") {
            syst->fMigrationNameSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationNameSuffDown");
        if (param != "") {
            syst->fMigrationNameSuffsDown.clear();
            syst->fMigrationNameSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationNameSuffsDown");
        if (param != "") {
            syst->fMigrationNameSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationFileUp");
        if (param != "") {
            syst->fMigrationFilesUp.clear();
            syst->fMigrationFilesUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationFilesUp");
        if (param != "") {
            syst->fMigrationFilesUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationFileDown");
        if (param != "") {
            syst->fMigrationFilesDown.clear();
            syst->fMigrationFilesDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationFilesDown");
        if (param != "") {
            syst->fMigrationFilesDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("MigrationFileSuffUp");
        if (param != "") {
            syst->fMigrationFileSuffsUp.clear();
            syst->fMigrationFileSuffsUp.emplace_back(Common::RemoveQuotes(param));
            hasUp = true;
        }

        param = confSet->Get("MigrationFileSuffsUp");
        if (param != "") {
            syst->fMigrationFileSuffsUp = Common::Vectorize(param,',');
            hasUp = true;
        }

        param = confSet->Get("MigrationFileSuffDown");
        if (param != "") {
            syst->fMigrationFileSuffsDown.clear();
            syst->fMigrationFileSuffsDown.emplace_back(Common::RemoveQuotes(param));
            hasDown = true;
        }

        param = confSet->Get("MigrationFileSuffsDown");
        if (param != "") {
            syst->fMigrationFileSuffsDown = Common::Vectorize(param,',');
            hasDown = true;
        }

        param = confSet->Get("ResponseMatrixPathUpSubtractSample");
        if (param != "") {
            syst->fResponseMatrixPathsUpSubtractSample.clear();
            syst->fResponseMatrixPathsUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixPathsUpSubtractSample");
        if (param != "") {
            syst->fResponseMatrixPathsUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixPathDownSubtractSample");
        if (param != "") {
            syst->fResponseMatrixPathsDownSubtractSample.clear();
            syst->fResponseMatrixPathsDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixPathsDownSubtractSample");
        if (param != "") {
            syst->fResponseMatrixPathsDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixFileUpSubtractSample");
        if (param != "") {
            syst->fResponseMatrixFilesUpSubtractSample.clear();
            syst->fResponseMatrixFilesUpSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("ResponseMatrixFilesUpSubtractSample");
        if (param != "") {
            syst->fResponseMatrixFilesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixFileDownSubtractSample");
        if (param != "") {
            syst->fResponseMatrixFilesDownSubtractSample.clear();
            syst->fResponseMatrixFilesDownSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("ResponseMatrixFilesDownSubtractSample");
        if (param != "") {
            syst->fResponseMatrixFilesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixNameUpSubtractSample");
        if (param != "") {
            syst->fResponseMatrixNamesUpSubtractSample.clear();
            syst->fResponseMatrixNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixNamesUpSubtractSample");
        if (param != "") {
            syst->fResponseMatrixNamesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixFolderNameUpSubtractSample");
        if (param != "") {
            syst->fResponseMatrixFolderNamesUpSubtractSample.clear();
            syst->fResponseMatrixFolderNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixNameDownSubtractSample");
        if (param != "") {
            syst->fResponseMatrixNamesDownSubtractSample.clear();
            syst->fResponseMatrixNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("ResponseMatrixNamesDownSubtractSample");
        if (param != "") {
            syst->fResponseMatrixNamesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("ResponseMatrixFolderNameDownSubtractSample");
        if (param != "") {
            syst->fResponseMatrixFolderNamesDownSubtractSample.clear();
            syst->fResponseMatrixFolderNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPathUpSubtractSample");
        if (param != "") {
            syst->fMigrationPathsUpSubtractSample.clear();
            syst->fMigrationPathsUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPathsUpSubtractSample");
        if (param != "") {
            syst->fMigrationPathsUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationPathDownSubtractSample");
        if (param != "") {
            syst->fMigrationPathsDownSubtractSample.clear();
            syst->fMigrationPathsDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationPathsDownSubtractSample");
        if (param != "") {
            syst->fMigrationPathsDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFileUpSubtractSample");
        if (param != "") {
            syst->fMigrationFilesUpSubtractSample.clear();
            syst->fMigrationFilesUpSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("MigrationFilesUpSubtractSample");
        if (param != "") {
            syst->fMigrationFilesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFileDownSubtractSample");
        if (param != "") {
            syst->fMigrationFilesDownSubtractSample.clear();
            syst->fMigrationFilesDownSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("MigrationFilesDownSubtractSample");
        if (param != "") {
            syst->fMigrationFilesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationNameUpSubtractSample");
        if (param != "") {
            syst->fMigrationNamesUpSubtractSample.clear();
            syst->fMigrationNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNamesUpSubtractSample");
        if (param != "") {
            syst->fMigrationNamesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFolderNameUpSubtractSample");
        if (param != "") {
            syst->fMigrationFolderNamesUpSubtractSample.clear();
            syst->fMigrationFolderNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNameDownSubtractSample");
        if (param != "") {
            syst->fMigrationNamesDownSubtractSample.clear();
            syst->fMigrationNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("MigrationNamesDownSubtractSample");
        if (param != "") {
            syst->fMigrationNamesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("MigrationFolderNameDownSubtractSample");
        if (param != "") {
            syst->fMigrationFolderNamesDownSubtractSample.clear();
            syst->fMigrationFolderNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPathUpSubtractSample");
        if (param != "") {
            syst->fSelectionEffPathsUpSubtractSample.clear();
            syst->fSelectionEffPathsUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPathsUpSubtractSample");
        if (param != "") {
            syst->fSelectionEffPathsUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffPathDownSubtractSample");
        if (param != "") {
            syst->fSelectionEffPathsDownSubtractSample.clear();
            syst->fSelectionEffPathsDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffPathsDownSubtractSample");
        if (param != "") {
            syst->fSelectionEffPathsDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFileUpSubtractSample");
        if (param != "") {
            syst->fSelectionEffFilesUpSubtractSample.clear();
            syst->fSelectionEffFilesUpSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("SelectionEffFilesUpSubtractSample");
        if (param != "") {
            syst->fSelectionEffFilesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFileDownSubtractSample");
        if (param != "") {
            syst->fSelectionEffFilesDownSubtractSample.clear();
            syst->fSelectionEffFilesDownSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("SelectionEffFilesDownSubtractSample");
        if (param != "") {
            syst->fSelectionEffFilesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffNameUpSubtractSample");
        if (param != "") {
            syst->fSelectionEffNamesUpSubtractSample.clear();
            syst->fSelectionEffNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNamesUpSubtractSample");
        if (param != "") {
            syst->fSelectionEffNamesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFolderNameUpSubtractSample");
        if (param != "") {
            syst->fSelectionEffFolderNamesUpSubtractSample.clear();
            syst->fSelectionEffFolderNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNameDownSubtractSample");
        if (param != "") {
            syst->fSelectionEffNamesDownSubtractSample.clear();
            syst->fSelectionEffNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("SelectionEffNamesDownSubtractSample");
        if (param != "") {
            syst->fSelectionEffNamesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("SelectionEffFolderNameDownSubtractSample");
        if (param != "") {
            syst->fSelectionEffFolderNamesDownSubtractSample.clear();
            syst->fSelectionEffFolderNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("AcceptancePathUpSubtractSample");
        if (param != "") {
            syst->fAcceptancePathsUpSubtractSample.clear();
            syst->fAcceptancePathsUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("AcceptancePathsUpSubtractSample");
        if (param != "") {
            syst->fAcceptancePathsUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("AcceptancePathDownSubtractSample");
        if (param != "") {
            syst->fAcceptancePathsDownSubtractSample.clear();
            syst->fAcceptancePathsDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("AcceptancePathsDownSubtractSample");
        if (param != "") {
            syst->fAcceptancePathsDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("AcceptanceFileUpSubtractSample");
        if (param != "") {
            syst->fAcceptanceFilesUpSubtractSample.clear();
            syst->fAcceptanceFilesUpSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("AcceptanceFilesUpSubtractSample");
        if (param != "") {
            syst->fAcceptanceFilesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("AcceptanceFileDownSubtractSample");
        if (param != "") {
            syst->fAcceptanceFilesDownSubtractSample.clear();
            syst->fAcceptanceFilesDownSubtractSample = Common::ToVec(param);
        }

        param = confSet->Get("AcceptanceFilesDownSubtractSample");
        if (param != "") {
            syst->fAcceptanceFilesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("AcceptanceNameUpSubtractSample");
        if (param != "") {
            syst->fAcceptanceNamesUpSubtractSample.clear();
            syst->fAcceptanceNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("AcceptanceNamesUpSubtractSample");
        if (param != "") {
            syst->fAcceptanceNamesUpSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("AcceptanceFolderNameUpSubtractSample");
        if (param != "") {
            syst->fAcceptanceFolderNamesUpSubtractSample.clear();
            syst->fAcceptanceFolderNamesUpSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("AcceptanceNameDownSubtractSample");
        if (param != "") {
            syst->fAcceptanceNamesDownSubtractSample.clear();
            syst->fAcceptanceNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        param = confSet->Get("AcceptanceNamesDownSubtractSample");
        if (param != "") {
            syst->fAcceptanceNamesDownSubtractSample = Common::Vectorize(param, ',');
        }

        param = confSet->Get("AcceptanceFolderNameDownSubtractSample");
        if (param != "") {
            syst->fAcceptanceFolderNamesDownSubtractSample.clear();
            syst->fAcceptanceFolderNamesDownSubtractSample = Common::ToVec(Common::RemoveQuotes(param));
        }

        syst->fHasUpVariation = hasUp;
        syst->fHasDownVariation = hasDown;

        param = confSet->Get("ReferenceSample");
        if (param != "") {
            syst->SetReferenceSample(Common::RemoveQuotes(param));
        }

        param = confSet->Get("DropShapeIn");
        if (param != "") {
            syst->fDropShapeIn = Common::Vectorize(param, ',');
        }

        param = confSet->Get("DropNorm");
        if (param != "") {
            syst->fDropNorm = Common::Vectorize(param, ',');
        }

        // Set Symmetrisation
        param = confSet->Get("Symmetrisation");
        if(param != ""){
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param == "ONESIDED"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZEONESIDED);
            }
            else if(param == "TWOSIDED"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZETWOSIDED);
            }
            else if(param == "ABSMEAN"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZEABSMEAN);
            }
            else if(param == "MAXIMUM"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZEMAXIMUM);
            }
            else if(param == "MAXIMUMKEEPSIGN"){
                syst->SetSymmetrisationType(HistoTools::SYMMETRIZEMAXIMUMKEEPSIGN);
            }
            else {
                LOG(ERROR) << "Symetrisation scheme is not recognized ... \n";
                ++sc;
            }
        }

        param = confSet->Get("SmoothingOption");
        if(param != ""){
            syst->fSampleSmoothing = true;
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if( param == "MAXVARIATION" ) syst->SetSmoothOption(HistoTools::SmoothOption::MAXVARIATION);
            else if (param == "TTBARRESONANCE") syst->SetSmoothOption(HistoTools::SmoothOption::TTBARRESONANCE);
            else if (param == "COMMONTOOLSMOOTHMONOTONIC") syst->SetSmoothOption(HistoTools::SmoothOption::COMMONTOOLSMOOTHMONOTONIC);
            else if (param == "COMMONTOOLSMOOTHPARABOLIC") syst->SetSmoothOption(HistoTools::SmoothOption::COMMONTOOLSMOOTHPARABOLIC);
            else if (param == "TCHANNEL") syst->SetSmoothOption(HistoTools::SmoothOption::TCHANNEL);
            else if (param == "KERNELRATIOUNIFORM") syst->SetSmoothOption(HistoTools::SmoothOption::KERNELRATIOUNIFORM);
            else if (param == "KERNELDELTAGAUSS") syst->SetSmoothOption(HistoTools::SmoothOption::KERNELDELTAGAUSS);
            else if (param == "KERNELRATIOGAUSS") syst->SetSmoothOption(HistoTools::SmoothOption::KERNELRATIOGAUSS);
            else {
                LOG(WARNING) << "You specified 'SmoothingOption' option but you didn't provide valid input. Using default from job block\n";
                syst->fSampleSmoothing = false;
            }
        }

        param = confSet->Get("Smoothing");
        if (param != "") {
            syst->SetSmoothingType(Common::convertStoNum<int>(param));
        }

        if (fFitter->fUnfolding.size() == 1) {
            syst->fUnfoldingName = fFitter->fUnfolding.at(0)->fName;
        } else {
            param = confSet->Get("UnfoldingName");
            if (param == "") {
                LOG(ERROR) << "Multiple unfoldings but no UnfoldingName provided for an UnfoldingSystematic\n";
                ++sc;
            } else {
                param = Common::RemoveQuotes(param);
                syst->fUnfoldingName = param;
                if (std::find_if(fFitter->fUnfolding.begin(), fFitter->fUnfolding.end(), [&param](const std::unique_ptr<Unfolding>& unf){return unf->fName == param;})
                    == fFitter->fUnfolding.end()) {
                    LOG(ERROR) << "UnfoldingName: " << param << " does not match any Unfolding names\n";
                    ++sc;
                }
            }
        }

        // Set ScaleUp
        param = confSet->Get("ScaleUp");
        if (param != "") {
            std::vector<std::string> temp_vec = Common::Vectorize(param,',',false);
            if (temp_vec.size() == 1 && Common::Vectorize(temp_vec.at(0), ':').size() == 1) {
                syst->fScaleUp = Common::convertStoNum<float>(param);
            } else{
                for (const std::string& ireg : temp_vec){
                    std::vector<std::string> reg_value = Common::Vectorize(ireg,':');
                    if (reg_value.size() == 2) {
                        syst->fScaleUpRegions.insert({reg_value.at(0), Common::convertStoNum<float>(reg_value.at(1).c_str())});
                    }
                }
            }
        }

        // Set ScaleDown
        param = confSet->Get("ScaleDown");
        if (param != "") {
            std::vector<std::string> temp_vec = Common::Vectorize(param,',',false);
            if (temp_vec.size() == 1 && Common::Vectorize(temp_vec.at(0), ':').size() == 1) {
                syst->fScaleDown = Common::convertStoNum<float>(param);
            } else{
                for (const std::string& ireg : temp_vec){
                    std::vector<std::string> reg_value = Common::Vectorize(ireg,':');
                    if (reg_value.size() == 2) {
                        syst->fScaleDownRegions.insert({reg_value.at(0), Common::convertStoNum<float>(reg_value.at(1).c_str())});
                    }
                }
            }
        }

        fFitter->fUnfoldingSystematics.emplace_back(std::move(syst));
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::UnfoldingCorrections() {
    if (fFitter->fFitType != TRExFit::FitType::UNFOLDING) return 0;

    int sc(0);

    // First process Samples
    sc += ProcessUnfoldingSamples();

    // Then process Systematics
    sc += ProcessUnfoldingSystematics();

    // Add norm factors
    sc += AddUnfoldingNormFactors();

    return sc;
}

//__________________________________________________________________________________
//
bool ConfigReader::ConfigHasNTUP(const ConfigSet* confSet){
    if (confSet->Get("Variable") != "" || confSet->Get("VariableForSample") != "" || confSet->Get("Selection") != "" || confSet->Get("NtupleName") != "" || confSet->Get("NtupleNameSuff") != "" || confSet->Get("MCweight") != "" || confSet->Get("NtuplePathSuff") != "" || confSet->Get("NtupleFile") != "" || confSet->Get("NtupleFiles") != "" || confSet->Get("NtupleNames") != "" || confSet->Get("NtuplePath") != "" || confSet->Get("NtuplePaths") != "" || confSet->Get("NtuplePathSuffs") != "" || confSet->Get("FriendPath") != "" || confSet->Get("FriendPaths") != "" || confSet->Get("FriendPathSuffs") != "") return true;
    else return false;
}

//__________________________________________________________________________________
//
bool ConfigReader::ConfigHasHIST(const ConfigSet* confSet){
    if (confSet->Get("HistoFile") != "" || confSet->Get("HistoName") != "" || confSet->Get("HistoPathSuff") != "" || confSet->Get("HistoPathSuffs") != "" || confSet->Get("HistoPath") != "" ) return true;
    else return false;
}


//__________________________________________________________________________________
//
bool ConfigReader::CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2){
    for (const auto& i : v1){
        if (i == "") continue;
        std::string s = i;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s == "NONE") continue;
        if (s == "ALL") continue;
        if (Common::FindInStringVector(v2, i) < 0){
            return false;
        }
    }

    return true;
}

//__________________________________________________________________________________
//
bool ConfigReader::CheckPresence(const std::vector<std::string> &v1, const std::vector<std::string> &v2, const std::vector<std::string> &v3){
    for (const auto& i : v1){
        if (i == "") continue;
        std::string s = i;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s == "NONE") continue;
        if (s == "ALL") continue;
        if (Common::FindInStringVector(v2, i) < 0){
            if (Common::FindInStringVector(v3, i) < 0){
                return false;
            }
        }
    }

    return true;
}


//__________________________________________________________________________________
//
std::vector<std::string> ConfigReader::GetAvailableRegions(){
    std::vector<std::string> availableRegions;
    const auto* region_confsets = fParser->GetConfigSets("Region");
    if(!region_confsets) return availableRegions;
    for(const auto& confSet : *region_confsets){
        std::string tmp = Common::RemoveQuotes(confSet->GetValue());
        availableRegions.emplace_back(Common::CheckName(tmp));
    }
    return availableRegions;
}

//__________________________________________________________________________________
//
std::vector<std::string> ConfigReader::GetAvailableSamples(){
    std::vector<std::string> availableSamples;
    const auto* sample_confsets = fParser->GetConfigSets("Sample");
    if(!sample_confsets) return availableSamples;
    for(const auto& confset : *sample_confsets){
        std::string tmp = Common::RemoveQuotes(confset->GetValue());
        availableSamples.emplace_back(Common::CheckName(tmp));
    }
    return availableSamples;
}

//__________________________________________________________________________________
//
std::vector<std::string> ConfigReader::GetAvailableSysts(){
    std::vector<std::string> availableSysts;
    const auto* syst_confsets = fParser->GetConfigSets("Systematic");
    if(!syst_confsets) return availableSysts;
    for(const auto& confset : *syst_confsets){
        std::string tmp = Common::RemoveQuotes(confset->GetValue());
        availableSysts.emplace_back(Common::CheckName(tmp));
    }
    return availableSysts;
}

//__________________________________________________________________________________
//
bool ConfigReader::SampleIsOk(const Sample* sample) const {

    if (sample->fType == Sample::SampleType::GHOST || sample->fType == Sample::SampleType::EFT) return true;

    if (std::find(fGhostSamples.begin(), fGhostSamples.end(), sample->fName) != fGhostSamples.end()) return false;
    if (std::find(fEFTSamples.begin(), fEFTSamples.end(), sample->fName) != fEFTSamples.end()) return false;
    if (sample->fDivideBy == sample->fName) return false;
    if (sample->fMultiplyBy == sample->fName) return false;
    if (sample->fNormToSample == sample->fName) return false;

    for (const auto& i : sample->fSubtractSamples) {
        if (i == sample->fName) return false;
    }

    for (const auto& i : sample->fAddSamples) {
        if (i == sample->fName) return false;
    }

    return true;
}

//__________________________________________________________________________________
//
bool ConfigReader::SystHasProblematicName(const std::string& name){
    if ((name.find("gamma") != std::string::npos)
            || (name.find("alpha") != std::string::npos)
            || (name.find("shape_") != std::string::npos)) {
        LOG(ERROR) << "NP " << name << " has a problematic name, please change it.\n";
        LOG(ERROR) << "You should not be using names with: \"gamma\", \"alpha\" or \"shape_\" as these are used internally and can cause problems.\n";
        return true;
    }

    return false;
}

//__________________________________________________________________________________
//
int ConfigReader::ProcessUnfoldingSamples() {

    int sc(0);

    for (const auto& ireg : fFitter->fRegions) {
        if (ireg->fRegionType != Region::RegionType::SIGNAL) continue;

        for (const auto& isample : fFitter->fUnfoldingSamples) {
            if(isample->fRegions[0] != "all" &&
                Common::FindInStringVector(isample->fRegions, ireg->fName) < 0) continue;
            if (isample->GetType() == UnfoldingSample::TYPE::GHOST) continue;

            const Unfolding* unfolding = fFitter->GetCorrespondingUnfolding(isample.get());
            if (!unfolding) {
                LOG(ERROR) << "The unfolding sample is nullptr!\n";
                exit(EXIT_FAILURE);
            }

            // Convert the UnfoldingSample to sample and adjust the paths
            const std::vector<std::shared_ptr<Sample> > samples = isample->ConvertToSample(ireg.get(), unfolding, fFitter->fName);
            // set regions to newly created sub-samples
            if (isample->fSubSampleRegions.size()>0) {
                for (unsigned int i_bin=1;i_bin<samples.size()+1;i_bin++) {
                    std::vector<std::string> regionVec;
                    // if Regions was "all", then set the new region list to what was set as SubSampleRegions
                    if (samples[i_bin-1]->fRegions[0] == "all") {
                        regionVec = isample->fSubSampleRegions[i_bin];
                    }
                    // otherwise, loop on all inidicated Regions and find which ones to keep
                    else {
                        for (const auto& region : samples[i_bin-1]->fRegions) {
                            if (Common::FindInStringVector(isample->fSubSampleRegions[i_bin],region) >= 0) {
                                regionVec.emplace_back(region);
                            }
                        }
                    }
                    // if regionVec is empy, set regions to "none"
                    if (regionVec.size()==0) {
                        regionVec.emplace_back("none");
                    }
                    samples[i_bin-1]->fRegions = regionVec;
                }
            }
            fFitter->fSamples.insert(fFitter->fSamples.end(), samples.begin(), samples.end());
        }

        for (const auto& unfolding : fFitter->fUnfolding) {
            // add custom asimov sample
            if (unfolding->fAlternativeAsimovTruthSample != "") {
                Sample* sample = new Sample("AlternativeSignal_"+unfolding->fName+"_"+ireg->fName, Sample::SampleType::GHOST);
                sample->fRegions = Common::ToVec(ireg->fName);
                sample->fHistoPaths = Common::ToVec(fFitter->fName+"/UnfoldingHistograms");
                sample->fHistoFiles = Common::ToVec("FoldedHistograms");
                const std::string histoName = unfolding->fName + "_" + ireg->fName + "_AlternativeAsimov";
                sample->fHistoNames = Common::ToVec(histoName);
                fFitter->fSamples.emplace_back(std::move(sample));
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ProcessUnfoldingSystematics() {

    int sc(0);

    for (const auto& ireg : fFitter->fRegions) {
        if (ireg->fRegionType != Region::RegionType::SIGNAL) continue;

        for (const auto& isyst : fFitter->fUnfoldingSystematics) {
            if(isyst->fRegions[0] != "all" &&
                Common::FindInStringVector(isyst->fRegions, ireg->fName) < 0) continue;

            const Unfolding* unfolding = fFitter->GetCorrespondingUnfolding(isyst.get());
            if (!unfolding) {
                LOG(ERROR) << "The unfolding sample is nullptr!\n";
                exit(EXIT_FAILURE);
            }

            if (fFitter->fUnfoldingSamples.size() == 0) {
                LOG(ERROR) << "No UnfoldingSamples set!\n";
                ++sc;
            }
            std::string unfoldingSampleName("");
            for (const auto& isample : fFitter->fUnfoldingSamples) {
                if(isample->fRegions[0] != "all" &&
                    Common::FindInStringVector(isample->fRegions, ireg->fName) < 0) continue;
                if (isample->GetType() == UnfoldingSample::TYPE::GHOST) continue;
                if (isample->fUnfoldingName != unfolding->fName) continue;

                unfoldingSampleName = isample->GetName();
            }

            if (unfoldingSampleName == "") {
                LOG(ERROR) << "Sample name not set!\n";
                ++sc;
            }

            const std::vector<std::shared_ptr<Systematic> > systs = isyst->ConvertToSystematic(ireg.get(),
                                                                                               unfolding,
                                                                                               fFitter->fName,
                                                                                               unfoldingSampleName,
                                                                                               fFitter->fSamples);

            fFitter->fSystematics.insert(fFitter->fSystematics.end(), systs.begin(), systs.end());
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::AddUnfoldingNormFactors() {

    int sc(0);

    // Add overall norm factor, for normalized xsec:
    for (auto& iunfolding : fFitter->fUnfolding) {
        if(iunfolding->fUnfoldNormXSec){
            std::shared_ptr<NormFactor> nfTot = std::make_shared<NormFactor>("TotalXsecOverTheory_"+iunfolding->fName);
            LOG(INFO) << "Adding overall norm-factor to all truth bins to get normalized cross-section\n";
            TRExFitter::SYSTMAP[nfTot->fName] = nfTot->fName;
            if(Common::FindInStringVector(fFitter->fNormFactorNames, nfTot->fName) < 0) {
               fFitter->fNormFactors.emplace_back(nfTot);
               fFitter->fNormFactorNames.emplace_back(nfTot->fName);
            }
            nfTot->fNuisanceParameter = Common::CheckName(nfTot->fName);
            TRExFitter::NPMAP[nfTot->fName] = nfTot->fName;
            nfTot->fCategory = "TotalXsec";
            nfTot->fTitle = "#sigma_{tot} / #sigma_{tot}^{theory}";
            nfTot->SetNominalMinMax(1., iunfolding->fUnfoldingResultMin, iunfolding->fUnfoldingResultMax);
            nfTot->fRegions = GetAvailableRegions();
            //
            for(auto& isample : fFitter->fSamples) {
                if (isample->fName.find("_Truth_bin_")==std::string::npos) continue;
                if (isample->fUnfoldingName != iunfolding->fName) continue;
                LOG(INFO) << "Adding to sample " << isample->fName << "\n";
                isample->AddNormFactor(nfTot);
            }
        }

        // in case no UnfoldNormXSecBinN was specified, set it to the last truth bin
        if(iunfolding->fUnfoldNormXSec && iunfolding->fUnfoldNormXSecBinN<=0) iunfolding->fUnfoldNormXSecBinN = iunfolding->fNumberUnfoldingTruthBins;

        for (int i = 0; i < iunfolding->fNumberUnfoldingTruthBins; ++i) {
            const std::string name = iunfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu"; // add leading zeros in order to have correct NF ordering in plots
            std::shared_ptr<NormFactor> nf = std::make_shared<NormFactor>(name);

            if(Common::FindInStringVector(fFitter->fNormFactorNames, nf->fName) < 0) {
               fFitter->fNormFactors.emplace_back(nf);
               fFitter->fNormFactorNames.emplace_back(nf->fName);
            }

            nf->fNuisanceParameter = Common::CheckName(nf->fName);
            TRExFitter::NPMAP[nf->fName] = nf->fName;
            nf->fCategory = "TruthBins";
            nf->fTitle = iunfolding->fName + " Unfolded Truth Bin "+std::to_string(i+1);
            TRExFitter::SYSTMAP[nf->fName] = nf->fTitle;
            nf->SetNominalMinMax(1., iunfolding->fUnfoldingResultMin, iunfolding->fUnfoldingResultMax);
            nf->fRegions = GetAvailableRegions();
            std::vector<std::string> sampleNames;
            for (const auto& ireg : fFitter->fRegions) {
                if (ireg->fRegionType == Region::RegionType::SIGNAL) {
                    const std::string sampleName = iunfolding->fName + "_" + ireg->fName + "_Truth_bin_" + std::to_string(i+1);
                    sampleNames.emplace_back(sampleName);
                }
            }
            for(auto& isample : fFitter->fSamples) {
                if (isample->fUnfoldingName != iunfolding->fName) continue;
                if (std::find(sampleNames.begin(), sampleNames.end(), isample->fName) == sampleNames.end()) continue;

                // add the unfolding NormFactor
                isample->AddNormFactor(nf);

                // now also add the custom NFs
                for (const auto& inf : fFitter->fNormFactors) {
                    // check if the sample name (or the sample's "UnfoldingSample" name) is in the sample-name list for this NF
                    if ( Common::FindInStringVector( inf->fSamples, isample->fName ) < 0
                        &&
                         Common::FindInStringVector( inf->fSamples, isample->fUnfoldingSampleName ) < 0
                       ) continue;

                    isample->AddNormFactor(inf);
                    LOG(INFO) << "Adding NF " << inf->fName << " to sample: " << isample->fName << "\n";
                }
            }

            // Add expression for last bin in case of norm xsec
            if(iunfolding->fUnfoldNormXSec && i==iunfolding->fUnfoldNormXSecBinN-1){
                LOG(INFO) << "Adding expression to turn cross section into normalized cross-section for Unfolding: " << iunfolding->fName << "\n";
                //
                // build expression:
                std::string v0 = "";
                std::string v1 = "";
                std::string num = "1.";
                std::string den = "N" + std::to_string(iunfolding->fUnfoldNormXSecBinN) + "/N";
                for (int j = 0; j < iunfolding->fNumberUnfoldingTruthBins; ++j){
                    if(j==iunfolding->fUnfoldNormXSecBinN-1) continue;
                    std::string NFname = iunfolding->fName +  "_Bin_" + Common::IntToFixLenStr(j+1) + "_mu";
                    num += "-" + NFname + "*(N" + std::to_string(j+1) + "/N)"; // "Nj/N" will be a place-holder here (replaced in TRExFit::
                    if(j!=0) v1 += ",";
                    v1 += NFname + "[1," + std::to_string(iunfolding->fUnfoldingResultMin) + "," + std::to_string(iunfolding->fUnfoldingResultMax) + "]";
                }
                v0 = "(" + num + ")/(" + den + ")";
                //
                nf->fExpression = std::make_pair(v0,v1);
                // title will contain the expression FIXME
                nf->fTitle = v0;
                TRExFitter::SYSTMAP[nf->fName] = v0;
                // nuis-par will contain the nuis-par of the norm factor the expression depends on FIXME
                nf->fNuisanceParameter = v1;
                TRExFitter::NPMAP[nf->fName] = v1;
            }
            else{
                // set this NF as POI
                fFitter->AddPOI(Common::CheckName(name));
                fFitter->fFitPOIAsimov[name] = 1;
            }
            // check if the bin is in the list specified taus
            auto itr = fTaus.find(iunfolding->fName);
            if (itr != fTaus.end()) {
                bool hasTau = false;
                double tauValue = 0.;
                if (itr->second.size()==1) {
                    if (itr->second[0].first==-1){
                        hasTau = true;
                        tauValue = itr->second[0].second;
                    } else if (i+1 == itr->second.at(0).first) {
                        hasTau = true;
                        tauValue = itr->second[0].second;
                    }
                }
                else {
                    auto it = std::find_if(itr->second.begin(), itr->second.end(),
                    [&i](const std::pair<int, double>& element){ return element.first == i+1;});
                    if (it != itr->second.end()) {
                        hasTau = true;
                        tauValue = it->second;
                    }
                }
                if (hasTau) {
                    LOG(INFO) << "Setting truth bin: " << (i+1) << " to use tau = " << tauValue << "\n";
                    nf->fTau = tauValue;
                }
            }
            // assign expressions to norm factors
            for (auto expressions : fExpressions) {
                if (expressions.first != iunfolding->fName) continue;
                for (auto expr : expressions.second) {
                    const std::string n = expr.first;
                    const std::vector<std::string>& v = Common::Vectorize( expr.second,':');
                    if (n!=nf->fName) continue;
                    if (v.size() != 2) {
                        LOG(ERROR) << "Wrong format for Expression for bin " << nf->fName << "\n";
                        ++sc;
                    }
                    else {
                        nf->fExpression = std::make_pair(v[0],v[1]);
                        // title will contain the expression FIXME
                        nf->fTitle = "Expression_" + v[0];
                        TRExFitter::SYSTMAP[nf->fName] = "Expression_" + v[0];
                        // nuis-par will contain the nuis-par of the norm factor the expression depends on FIXME
                        nf->fNuisanceParameter = "Expression_" + v[1];
                        TRExFitter::NPMAP[nf->fName] = "Expression_" + v[1];
                    }
                }
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::CheckPOIs() const {
    int sc(0);
    for (const auto& ipoi : fFitter->fPOIs) {
        auto it = std::find_if(fFitter->fNormFactors.begin(), fFitter->fNormFactors.end(),
            [&ipoi](const std::shared_ptr<NormFactor>& nf) {
                const std::string name = nf->fName;
                return ipoi == name;
            });

        if (it != fFitter->fNormFactors.end()) {
            continue;
        }

        // check if the is the reparametrised SF or Template
        bool found(false);
        for (const auto& isample : fFitter->fSamples) {
            for (const auto& i : isample->fTemplateMorphing) {
                if (i.first == ipoi) {
                    found = true;
                    break;
                }
            }
        }

        if (found) continue;

        // check the EFT reparametrisation
        if (fFitter->fFitType == TRExFit::FitType::EFT && !fFitter->fEFTConfig.GetSplitSamplesPerBin()) {
            auto itr = std::find(fEFTParams.begin(), fEFTParams.end(), ipoi);
            if (itr != fEFTParams.end()) {
                found = true;
                continue;
            }
        }

        for (const auto& isf : fFitter->fShapeFactors) {
            for (const auto& i : isf->fExpression) {
                if (i.first == ipoi) {
                    found = true;
                    break;
                }
            }
        }

        if (found) continue;

        // check if it is NP
        auto itNP = std::find_if(fFitter->fSystematics.begin(), fFitter->fSystematics.end(),
            [&ipoi](const std::shared_ptr<Systematic>& syst) {
                const std::string name = syst->fNuisanceParameter;
                return (name == ipoi || ipoi == "alpha_"+name);
            });

        if (itNP == fFitter->fSystematics.end()) {
            LOG(ERROR) << "POI: " << ipoi << "  is not NF nor NP!\n";
            ++sc;
        } else {
            LOG(WARNING) << "POI: " << ipoi << " is NP\n";
            LOG(WARNING) << "If this is intended, please ignore this warning\n";
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
void ConfigReader::CheckTemplates() const {

    /// region | temlate target | template name | number of occurences
    std::map< std::string, std::map<std::string, std::map<std::vector<std::string>, int> > > map;
    for (const auto& ireg : fFitter->fRegions) {
        std::map<std::string, std::map<std::vector<std::string>, int> > tmp;
        map.insert({ireg->fName, tmp});
    }

    for (const auto& isample : fFitter->fSamples) {
        if (isample->fType == Sample::SampleType::GHOST) continue;
        if (isample->fTemplateMorphing.empty()) continue;
        std::map<std::vector<std::string>, int> tmp;
        const std::string target = isample->fTemplateTarget;
        for (auto& imap : map) {
            auto itr = imap.second.find(target);
            if (itr == imap.second.end()) {
                imap.second.insert({target, tmp});
            }

        }
    }

    auto GetParams = [](const std::vector<std::pair<std::string, double> >& in){
        std::vector<std::string> result;
        for (const auto& i : in) {
            result.emplace_back(i.first);
        }

        return result;
    };

    for (const auto& isample : fFitter->fSamples) {
        if (isample->fType == Sample::SampleType::GHOST) continue;
        if (isample->fTemplateMorphing.empty()) continue;
        const std::vector<std::string> params = GetParams(isample->fTemplateMorphing);
        const std::string target = isample->fTemplateTarget;
        if (isample->fRegions.size() == 1 && isample->fRegions.at(0) == "all") {
            for (const auto& ireg : fFitter->fRegions) {
                auto itr = map.find(ireg->fName);
                auto itr2 = itr->second.find(target);
                auto itr3 = itr2->second.find(params);
                if (itr3 == itr2->second.end()) {
                    itr2->second.insert({params, 1});
                } else {
                    itr3->second++;
                }
            }
        } else {
            for (const auto& ireg : isample->fRegions) {
                auto itr = map.find(ireg);
                auto itr2 = itr->second.find(target);
                auto itr3 = itr2->second.find(params);
                if (itr3 == itr2->second.end()) {
                    itr2->second.insert({params, 1});
                } else {
                    itr3->second++;
                }
            }
        }
    }

    for (const auto& imap : map) {
        for (const auto& i : imap.second) {
            for (const auto& j : i.second) {
                if (j.second != 1) {
                    LOG(ERROR) << "Region: " << imap.first << " has templates applied to target: " << i.first << " " << j.second << " times to non-GHOST samples, while only one is allowed\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    for (const auto& isample : fFitter->fSamples) {
        if (isample->fTemplateMorphing.empty()) continue;
        if (!isample->fTemplateTarget.empty()) {
            const std::string sampleName = isample->fName;
            const std::string target = isample->fTemplateTarget;
            auto itr = std::find_if(fFitter->fSamples.begin(), fFitter->fSamples.end(), [&target](const auto& element){ return target == element->fName;});
            if (itr == fFitter->fSamples.end()) {
                LOG(ERROR) << "Template target " << target << " for sample " << sampleName << " does not match any sample\n";
                exit(EXIT_FAILURE);
            }
            if ((*itr)->fType != Sample::SampleType::SIGNAL) {
                LOG(ERROR) << "Template target " << target << " for sample " << sampleName << " does not match a SIGNAL sample\n";
                exit(EXIT_FAILURE);
            }
            if ((*itr)->fTemplateMorphing.empty()) {
                LOG(ERROR) << "Template target " << target << " for sample " << sampleName << " does not match a template SIGNAL sample\n";
                exit(EXIT_FAILURE);
            }
        }
    }
}

//__________________________________________________________________________________
//
void ConfigReader::ProcessTemplateMorphing() const {
    for (const auto& isample : fFitter->fSamples) {
        if (isample->fTemplateMorphing.empty()) continue;

        if (isample->fTemplateTarget.empty()) {
            // find the only sample with the template that is not a ghost
            for (const auto& isampleAll : fFitter->fSamples) {
                if (isampleAll->fType != Sample::SampleType::SIGNAL) continue;
                if (isampleAll->fTemplateMorphing.empty()) continue;
                isample->fTemplateTarget = isampleAll->fName;
                break;
            }
        }
        if (isample->fType == Sample::SampleType::GHOST) continue;
        if (isample->fRegions.at(0) == "all") {
            for (const auto& ireg : fFitter->fRegions) {
                auto sf = std::make_shared<ShapeFactor>("templateSF_"+ireg->fName+"_"+isample->fName);
                sf->fRegions = Common::ToVec(ireg->fName);
                sf->fSamples = Common::ToVec(isample->fName);
                isample->AddShapeFactor(sf);
                fFitter->fShapeFactorNames.emplace_back("templateSF_"+ireg->fName+"_"+isample->fName);
                fFitter->fShapeFactors.emplace_back(sf);
            }
        } else {
            for (const auto& ireg : isample->fRegions) {
                auto sf = std::make_shared<ShapeFactor>("templateSF_"+ireg+"_"+isample->fName);
                sf->fRegions = Common::ToVec(ireg);
                sf->fSamples = Common::ToVec(isample->fName);
                isample->AddShapeFactor(sf);
                fFitter->fShapeFactorNames.emplace_back("templateSF_"+ireg+"_"+isample->fName);
                fFitter->fShapeFactors.emplace_back(sf);
            }
        }
    }
}

//__________________________________________________________________________________
//
int ConfigReader::ReadMorphing() {

    int sc(0);

    const ConfigSet *confSet = fParser->GetConfigSet("Morphing");
    if (confSet == nullptr) return sc;

    std::string param = confSet->Get("FitFunction");
    if (param != "") {
        fFitter->fMorphingSetting.SetMorphingFunction(param);
    }

    param = confSet->Get("ParamStartValue");
    if (param != "") {
        const auto params = Common::Vectorize(param, ',');
        for (const auto& iparam : params) {
            const auto oneParam = Common::Vectorize(iparam, '@');
            if (oneParam.size() != 3) {
                LOG(ERROR) << "Problem parsing \"ParamStartValue\"\n";
                ++sc;
                break;
            }

            const std::string target = oneParam.at(0);
            auto itr = std::find(fTemplateTargets.begin(), fTemplateTargets.end(), target);
            if (itr == fTemplateTargets.end()) {
                if (target == "ALL") {
                    LOG(INFO) << "ALL morphing targets specified for ParamStartValues\n";
                }
                else {
                    LOG(ERROR) << "Cannot find target: " << target << " in the list of the morphing templates\n";
                    ++sc;
                    break;
                }
            }
            const std::string region = oneParam.at(1);
            if (region.find("SIGNAL") != std::string::npos || region.find("CONTROL") != std::string::npos || region.find("VALIDATION") != std::string::npos ) {
                LOG(INFO) << "Setting ParamStartValues across regions \n";
            }
            const auto tmp = Common::Vectorize(oneParam.at(2), ':');
            if (tmp.size() != 3) {
                LOG(ERROR) << "Problem parsing \"ParamStartValue\"\n";
                ++sc;
                break;
            }

            int index(0);
            int bin(0);
            double value(0);
            try {
                // ALL bins in config indicates use across bins, set to -999 behind scenes
                if (tmp.at(0)=="ALL") bin = -999;
                else bin = Common::convertStoNum<int>(tmp.at(0));
                index = Common::convertStoNum<int>(tmp.at(1));
                value = Common::convertStoNum<double>(tmp.at(2));
            } catch (...) {
                LOG(ERROR) << "Problem parsing \"ParamStartValue\"\n";
                ++sc;
                break;
            }
            if (target != "ALL" && region.find("SIGNAL") == std::string::npos && region.find("CONTROL") == std::string::npos && region.find("VALIDATION") == std::string::npos) {
                fFitter->fMorphingSetting.AddStartValues(target, region, bin, index, value);
            }
            else {
                std::vector<std::string> target_tmp;
                std::vector<std::string> regions_tmp;
                if (target == "ALL") {
                    target_tmp = fTemplateTargets;
                }
                else {
                    target_tmp.emplace_back(target);
                }
                if (region.find("SIGNAL") == std::string::npos && region.find("CONTROL") == std::string::npos && region.find("VALIDATION") == std::string::npos) {
                    regions_tmp.emplace_back(region);
                }
                else {
                    for (const auto& rgns : GetAvailableRegions()) {
                        const Region *regObj = fFitter->GetRegion(rgns);
                        if (region.find("SIGNAL") != std::string::npos && regObj->fRegionType == Region::RegionType::SIGNAL) {
                            regions_tmp.emplace_back(rgns);
                        }
                        if (region.find("CONTROL") != std::string::npos  && regObj->fRegionType == Region::RegionType::CONTROL) {
                            regions_tmp.emplace_back(rgns);
                        }
                        if (region.find("VALIDATION") != std::string::npos && regObj->fRegionType == Region::RegionType::VALIDATION) {
                            regions_tmp.emplace_back(rgns);
                        }
                    }
                }
                for (const auto& trgt : target_tmp) {
                    for (const auto& rgn : regions_tmp) {
                        fFitter->fMorphingSetting.AddStartValues(trgt, rgn, bin, index, value);
                    }
                }
            }
        }
    }

    param = confSet->Get("ParamLimits");
    if (param != "") {
        const auto params = Common::Vectorize(param, ',');
        for (const auto& iparam : params) {
            const auto oneParam = Common::Vectorize(iparam, '@');
            if (oneParam.size() != 3) {
                LOG(ERROR) << "Problem parsing \"ParamLimits\"\n";
                ++sc;
                break;
            }

            const std::string target = oneParam.at(0);
            auto itr = std::find(fTemplateTargets.begin(), fTemplateTargets.end(), target);
            if (itr == fTemplateTargets.end()) {
                if (target == "ALL") {
                    LOG(INFO) << "ALL morphing targets specified for ParamStartValues\n";
                }
                else {
                    LOG(ERROR) << "Cannot find target: " << target << " in the list of the morphing templates\n";
                    ++sc;
                    break;
                }
            }
            const std::string region = oneParam.at(1);
            if (region.find("SIGNAL") != std::string::npos || region.find("CONTROL") != std::string::npos || region.find("VALIDATION") != std::string::npos ) {
                LOG(INFO) << "Setting ParamStartValues across regions \n";
            }
            const auto tmp = Common::Vectorize(oneParam.at(2), ':');
            if (tmp.size() != 4) {
                LOG(ERROR) << "Problem parsing \"ParamLimits\"\n";
                ++sc;
                break;
            }

            int index(0);
            int bin(0);
            double min(0);
            double max(0);
            try {
                // ALL bins in config indicates use across bins, set to -999 behind scenes
                if (tmp.at(0)=="ALL") bin = -999;
                else bin = Common::convertStoNum<int>(tmp.at(0));
                index = Common::convertStoNum<int>(tmp.at(1));
                min = Common::convertStoNum<double>(tmp.at(2));
                max = Common::convertStoNum<double>(tmp.at(3));
            } catch (...) {
                LOG(ERROR) << "Problem parsing \"ParamLimits\"\n";
                ++sc;
                break;
            }
            if (target != "ALL" && region.find("SIGNAL") == std::string::npos && region.find("CONTROL") == std::string::npos && region.find("VALIDATION") == std::string::npos) {
                fFitter->fMorphingSetting.AddParamLimits(target, region, bin, index, min, max);
            }
            else {
                std::vector<std::string> target_tmp;
                std::vector<std::string> regions_tmp;
                if (target == "ALL") {
                    target_tmp = fTemplateTargets;
                }
                else {
                    target_tmp.emplace_back(target);
                }
                if (region.find("SIGNAL") == std::string::npos && region.find("CONTROL") == std::string::npos && region.find("VALIDATION") == std::string::npos) {
                    regions_tmp.emplace_back(region);
                }
                else {
                    for (const auto& rgns : GetAvailableRegions()) {
                        const Region *regObj = fFitter->GetRegion(rgns);
                        if (region.find("SIGNAL") != std::string::npos && regObj->fRegionType == Region::RegionType::SIGNAL) {
                            regions_tmp.emplace_back(rgns);
                        }
                        if (region.find("CONTROL") != std::string::npos  && regObj->fRegionType == Region::RegionType::CONTROL) {
                            regions_tmp.emplace_back(rgns);
                        }
                        if (region.find("VALIDATION") != std::string::npos && regObj->fRegionType == Region::RegionType::VALIDATION) {
                            regions_tmp.emplace_back(rgns);
                        }
                    }
                }
                for (const auto& trgt : target_tmp) {
                    for (const auto& rgn : regions_tmp) {
                        fFitter->fMorphingSetting.AddParamLimits(trgt, rgn, bin, index, min, max);
                    }
                }
            }
        }
    }

    param = confSet->Get("ParamRangeScale");
    if (param != "") {
        fFitter->fMorphingSetting.SetRangeScale(Common::convertStoNum<double>(param));
    }

    param = confSet->Get("Chi2WarningLimit");
    if (param != "") {
        fFitter->fMorphingSetting.SetChi2WarningLimit(Common::convertStoNum<double>(param));
    }

    param = confSet->Get("Make1DProjections");
    if (param != ""){
        fFitter->fMorphingSetting.SetMake1DProjections(Common::StringToBoolean(param));
    }

    param = confSet->Get("Expressions");
    if (param != ""){
        const std::vector<std::string> expressions = Common::Vectorize(param, ',');
        if (expressions.size() > 2) {
            LOG(ERROR) << "\"Expressions\" option has > 2 parameters reparametrised, this is not supported\n";
            ++sc;
        } else {
            for (const auto& iexpression : expressions) {
                const std::vector<std::string> oneExpression = Common::Vectorize(iexpression, '=');
                if (oneExpression.size() != 2) {
                    LOG(ERROR) << "\"Expressions\" option has a wrong format (=)\n";
                    ++sc;
                } else {
                    const std::string parameter = oneExpression.at(0);
                    const std::vector<std::string> formula = Common::Vectorize(oneExpression.at(1), ':');
                    if (formula.size() != 2) {
                       LOG(ERROR) << "\"Expressions\" option has a wrong format (:)\n";
                       ++sc;
                    } else {
                        fFitter->fMorphingSetting.AddExpression(parameter, std::make_pair(formula.at(0), formula.at(1)));
                    }
                }
            }
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReader::ReadEFTConfig() {

    int sc(0);

    const ConfigSet *confSet = fParser->GetConfigSet("EFTConfig");
    if (confSet == nullptr) return sc;

    std::string param = confSet->Get("ParameterMinMax");
    if (param != "") {
        const std::vector<std::string> par = Common::Vectorize(param, ',');
        for (const auto& ipar : par) {
            const std::vector<std::string> parMinMax = Common::Vectorize(ipar, ':');
            if (parMinMax.size() != 3) {
                LOG(ERROR) << "Wrong input for \"ParameterMinMax\"\n";
                ++sc;
            } else {
                const std::string& parName = parMinMax.at(0);
                const double min = Common::convertStoNum<double>(parMinMax.at(1));
                const double max = Common::convertStoNum<double>(parMinMax.at(2));

                fFitter->fEFTConfig.AddParamMinMax(parName, min, max);
            }
        }
    }

    param = confSet->Get("SplitSamplePerBin");
    if (param != ""){
        bool flag = Common::StringToBoolean(param);
        if (flag && fFitter->fInputType == 0) {
            LOG(ERROR) << "Using histograms input and setting \"EFTSplitSamplePerBin: TRUE\" not supported, switching to FALSE\n";
            flag = false;
        }

        fFitter->fEFTConfig.SetSplitSamplesPerBin(flag);
    }

    // Set EFTOrder
    param = confSet->Get("EFTOrder");
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);
    if (param != ""){
        if (param == "ALL" ) {
            fFitter->fEFTConfig.SetEFTOrder(EFTConfig::ALL);
        } else if (param == "LIN" || param == "LINEAR" ){
            fFitter->fEFTConfig.SetEFTOrder(EFTConfig::LIN);
        } else {
            LOG(ERROR) << "Unknown EFTOrder argument : " << param << "\n";
            ++sc;
        }
    } else {
        LOG(VERBOSE) << "Setting default EFT order ALL\n";
        fFitter->fEFTConfig.SetEFTOrder(EFTConfig::ALL);
    }

    param = confSet->Get("ExcludedParameterRegionSample");
    if (param != ""){
        const std::vector<std::string> exclude = Common::Vectorize(param, ',');
        for (const auto& iexclude : exclude) {
            const std::vector<std::string> tmp = Common::Vectorize(iexclude, ':');
            if (tmp.size() != 3) {
                LOG(ERROR) << "Wrong setting for option \"ExcludedParameterRegionSample\"\n";
                ++sc;
            } else {
                fEFTExclude.emplace_back(std::move(tmp));
            }
        }
    }

    param = confSet->Get("PlotRatio");
    if (param != ""){
        fFitter->fEFTConfig.SetEFTPlotRatio(Common::StringToBoolean(param));
    }

    return sc;
}

//__________________________________________________________________________________
//
void ConfigReader::CheckNormFactors() const {
    for (const auto& isample : fFitter->fSamples) {
        for (const auto& ireg : isample->fRegions) {
            std::vector<std::shared_ptr<NormFactor> > nfs;
            for (const auto& inf : isample->fNormFactors) {
                const auto& regions = inf->fRegions;
                // skip region/NF combination
                if (regions.empty()) continue;
                if (!(regions.size() == 1 && regions.at(0) == "all")) {
                    if (Common::FindInStringVector(regions, ireg) < 0) continue;
                }
                auto itr = std::find_if(nfs.begin(), nfs.end(), [&inf](const auto& element){return element->fName == inf->fName;});
                if (itr != nfs.end()) {
                    LOG(ERROR) << "Duplicate Norm Factor: " << inf->fName << ", for Sample: " << isample->fName << ", region: " << ireg << "\n";
                    exit(EXIT_FAILURE);
                } else {
                    nfs.emplace_back(inf);
                }
            }
        }
    }
}

//__________________________________________________________________________________
//
void ConfigReader::ProcessEFTShapeFactors() const {

    std::vector<std::string> referenceSamples;
    for (const auto& isample : fFitter->fSamples) {
        if (isample->fType != Sample::SampleType::EFT) continue;
        const std::string ref = isample->fEFTSMReference;
        auto itr = std::find (referenceSamples.begin(), referenceSamples.end(), ref);
        if (itr == referenceSamples.end()) {
            referenceSamples.emplace_back(ref);
        }
    }

    std::vector<std::string> regions;
    for (const auto& ireg : fFitter->fRegions) {
        if (ireg->fRegionType == Region::RegionType::VALIDATION) continue;
        regions.emplace_back(ireg->fName);
    }

    fFitter->fEFTConfig.SetReferenceSample(referenceSamples);
    fFitter->fEFTConfig.SetRegions(regions);
    fFitter->fEFTConfig.SetParameters(fEFTParams);
    fFitter->fEFTConfig.SetOperatorsPerSample(fEFTOperatorsPerSample);

    for (const auto& i : fEFTExclude) {
        fFitter->fEFTConfig.AddExcludedRegionSample(i.at(0), i.at(1), i.at(2));
    }

    for (const auto& isample : fFitter->fSamples) {
        if (isample->fType == Sample::SampleType::EFT) continue;
        auto itr = std::find(referenceSamples.begin(), referenceSamples.end(), isample->fName);
        if (itr == referenceSamples.end()) continue;
        if (isample->fRegions.at(0) == "all") {
            for (const auto& ireg : fFitter->fRegions) {
                if (fFitter->fEFTConfig.IsFullyExcluded(ireg->fName, isample->fName)) {
                    LOG(DEBUG) << "Region: " << ireg->fName << ", sample: " << isample->fName << ", is fully excluded for EFT\n";
                    continue;
                }
                LOG(DEBUG) << "Attaching ShapeFactors for EFT to to sample: " << isample->fName << ", region: " << ireg->fName << "\n";
                auto sf = std::make_shared<ShapeFactor>("EFTSF_"+ireg->fName+"_"+isample->fName);
                sf->fRegions = Common::ToVec(ireg->fName);
                sf->fSamples = Common::ToVec(isample->fName);
                isample->AddShapeFactor(sf);
                fFitter->fShapeFactorNames.emplace_back("EFTSF_"+ireg->fName+"_"+isample->fName);
                fFitter->fShapeFactors.emplace_back(sf);
            }
        } else {
            for (const auto& ireg : isample->fRegions) {
                if (fFitter->fEFTConfig.IsFullyExcluded(ireg, isample->fName)) {
                    LOG(DEBUG) << "Region: " << ireg << ", sample: " << isample->fName << ", is fully excluded for EFT\n";
                    continue;
                }
                LOG(DEBUG) << "Attaching ShapeFactors for EFT to to sample: " << isample->fName << ", region: " << ireg << "\n";
                auto sf = std::make_shared<ShapeFactor>("EFTSF_"+ireg+"_"+isample->fName);
                sf->fRegions = Common::ToVec(ireg);
                sf->fSamples = Common::ToVec(isample->fName);
                isample->AddShapeFactor(sf);
                fFitter->fShapeFactorNames.emplace_back("EFTSF_"+ireg+"_"+isample->fName);
                fFitter->fShapeFactors.emplace_back(sf);
            }
        }

        referenceSamples.erase(itr);
    }

    if (!referenceSamples.empty()) {
        LOG(ERROR) << "Some reference samples are not found.\n";
        for (const auto& iref : referenceSamples) {
            LOG(ERROR) << "\t missing reference sample: " << iref << "\n";
        }

        exit(EXIT_FAILURE);
    }

    fFitter->fShapeFactorReparametrisation = true;
}
