// Class include
#include "TRExFitter/MultiFit.h"

// Framework includes
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReaderMulti.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/TRExFit.h"

// ROOT includes
#include "TSystem.h"

// c++ includes
#include <algorithm>

#define _STRINGIZE(x) #x
#define STRINGIZE(x) _STRINGIZE(x)

//_______________________________________________________________________________________
//
ConfigReaderMulti::ConfigReaderMulti(MultiFit *multiFitter) :
    fMultiFitter(multiFitter),
    fGlobalSuffix(""),
    fRunLHscan(true),
    fLHscanVariable(""),
    fLHscanStep(-1),
    fLHscanStepY(-1){
    LOG(INFO) << "Started reading the config for multifit\n";
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadFullConfig(const std::string& fileName, const std::string& opt, const std::string& option){
    // initialize ConfigParser for the actual config
    fParser.ReadFile(fileName);

    // initialize checker COnfigParser to cross check the input
    ConfigParser refConfig;
    TString homeArea(gSystem->Getenv("TREXFITTER_HOME") ? "$TREXFITTER_HOME" : "$DATAPATH");
#ifdef TREXFITTER_HOME
    homeArea = std::string(STRINGIZE(TREXFITTER_HOME));
#endif
    TString file = "multiFitSchema.config";
    auto* expandPath = gSystem->ExpandPathName(homeArea+":./");
    gSystem->FindFile(expandPath,file);
    if (file == "") {
        LOG(ERROR) << "jobSchema.config not found. Are you running from the main working directory?\n";
        exit(EXIT_FAILURE);
    }
    delete [] expandPath;
    refConfig.ReadFile(static_cast<std::string>(file.Data()));
    int sc = fParser.CheckSyntax(&refConfig);

    if (sc != 0) return sc;

    // syntax of the config is ok
    // read different types of settings
    if (option != ""){
        sc+= ReadCommandLineOptions(option);
    }

    sc+= ReadJobOptions();

    sc+= ReadLimitOptions();

    sc+= ReadSignificanceOptions();

    sc+= ReadFitOptions(opt, option);

    if (option != ""){
        sc+= ReadCommandLineOptions(option);
    }

    sc += PostConfig();

    // make directory
    gSystem->mkdir(fMultiFitter->fOutDir.c_str());

    return sc;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadCommandLineOptions(const std::string &option){
    // Read options (to skip stuff, or include only some regions, samples, systs...)
    // Syntax: .. .. Regions=ge4jge2b:Exclude=singleTop,wjets
    std::map< std::string,std::string > optMap;

    int sc(0);

    std::vector< std::string > optVec = Common::Vectorize(option,':');
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

    if(optMap["Ranking"]!=""){
        fMultiFitter->fRankingOnly = optMap["Ranking"];
    }
    if(optMap["GroupedImpact"]!=""){
        fMultiFitter->fGroupedImpactCategory = optMap["GroupedImpact"];
    }

    if(optMap["Suffix"]!=""){
        fGlobalSuffix = optMap["Suffix"];
    }

    if(optMap["RunLHscan"]!="" && optMap["RunLHscan"]=="FALSE"){
        fRunLHscan = false;
    }
    if(optMap["LHscanVariable"]!=""){
        fLHscanVariable = optMap["LHscanVariable"];
    }
    if(optMap["LHscanStep"] != ""){
        fLHscanStep = atoi(optMap["LHscanStep"].c_str());
    }
    if(optMap["LHscanStepY"] != ""){
        fLHscanStepY = atoi(optMap["LHscanStepY"].c_str());
    }

    if(optMap["Regions"]!=""){
        fMultiFitter->fOnlyRegions = Common::Vectorize(optMap["Regions"],',');
    }
    if(optMap["BlindedParameters"]!=""){
        fMultiFitter->fBlindedParams = Common::Vectorize(optMap["BlindedParameters"],',');
    }
    if(optMap["NumCPU"]!=""){
        try {
            fMultiFitter->fCPU = Common::convertStoNum<int>(optMap["NumCPU"]);
            if (fMultiFitter->fCPU <= 0) {
                LOG(ERROR) << "NumCPU must be a positive integer!\n";
                ++sc;
            }
        }
        catch (const std::invalid_argument & e) {
            LOG(ERROR) << "NumCPU must be an integer!\n";
            ++sc;
        }
    } else if (optMap["SignificanceToysSeed"] != "") {
        try {
            fMultiFitter->fSignificanceToysSeed = Common::convertStoNum<int>(optMap["SignificanceToysSeed"]);
        }
        catch (const std::invalid_argument & e) {
            LOG(ERROR) << "SignificanceToysSeed must be an int!\n";
            ++sc;
        }
    } else if (optMap["LimitToysSeed"] != "") {
        try {
            fMultiFitter->fLimitToysSeed = Common::convertStoNum<int>(optMap["LimitToysSeed"]);
        }
        catch (const std::invalid_argument & e) {
            LOG(ERROR) << "LimitToysSeed must be an int!\n";
            ++sc;
        }
    } else if(optMap["ToysSeed"] != ""){
         fMultiFitter->fToysSeed = Common::convertStoNum<int>(optMap["ToysSeed"]);
    }

    return sc;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadJobOptions() {

    int sc(0);

    const ConfigSet *confSet = fParser.GetConfigSet("MultiFit");
    if (confSet == nullptr){
        LOG(ERROR) << "Cannot find 'MultiFit' in your config which is required. Please check this!\n";
        ++sc;
        return sc;
    }

    fMultiFitter->fName = Common::CheckName(confSet->GetValue());

    // Set POIName
    std::string param = confSet->Get("POIName");
    if( param != "" ) {
        const auto tmp =  Common::Vectorize(param, ',');
        for (const auto& i : tmp) {
            fMultiFitter->AddPOI(i);
        }
    }

    // Set OutputDir
    param = confSet->Get("OutputDir");
    if(param !=""){
        fMultiFitter->fDir = Common::CheckName(param);
        if(fMultiFitter->fDir.back() != '/') fMultiFitter->fDir += '/';
        fMultiFitter->fOutDir = fMultiFitter->fDir + fMultiFitter->fName;
        gSystem->mkdir((fMultiFitter->fOutDir).c_str(), true);
    }
    else{
        fMultiFitter->fOutDir = "./" + fMultiFitter->fName;
    }

    // Set Label
    param = confSet->Get("Label");
    if( param != "") fMultiFitter->fLabel = Common::RemoveQuotes(param);

    // Set LumiLabel
    param = confSet->Get("LumiLabel");
    if( param != "") fMultiFitter->fLumiLabel = Common::RemoveQuotes(param);

    // Set CmeLabel
    param = confSet->Get("CmeLabel");
    if( param != "") fMultiFitter->fCmeLabel = Common::RemoveQuotes(param);

    // Set Label of combination
    param = confSet->Get("CombiLabel");
    if( param != "") fMultiFitter->fCombiLabel = Common::RemoveQuotes(param);

    // Set SaveSuf
    param = confSet->Get("SaveSuf");
    if( param != "") fMultiFitter->fSaveSuf = Common::RemoveQuotes(param);
    else fMultiFitter->fSaveSuf             = fGlobalSuffix;

    // Set ShowObserved
    param = confSet->Get("ShowObserved");
    if( param != ""){
        fMultiFitter->fShowObserved = Common::StringToBoolean(param);
    }

    // Set LimitTitle
    param = confSet->Get("LimitTitle");
    if( param != "") fMultiFitter->fLimitTitle = Common::RemoveQuotes(param);
    // --- these lines should be removed at some point, since now we protect the text inside ""
    if(fMultiFitter->fLimitTitle.find("95CL")!=std::string::npos){
         fMultiFitter->fLimitTitle.replace(fMultiFitter->fLimitTitle.find("95CL"),4,"95% CL");
    }
    // ---

    // Ser POITitle
    param = confSet->Get("POITitle");
    if (param != "") {
        const auto vec = Common::Vectorize(param, ',');
        if (vec.size () != fMultiFitter->fPOIs.size()) {
            LOG(ERROR) << "You specified 'POITitle' option but you didn't pass the same number of parameters as the number of POIs. Please check this!\n";
            ++sc;
        }
        fMultiFitter->fPOITitle = vec;
     } else {
        for (std::size_t i = 0; i < fMultiFitter->fPOIs.size(); ++i) {
            fMultiFitter->fPOITitle.emplace_back("#mu");
        }
     }

    // Set CompareLimits
    param = confSet->Get("CompareLimits");
    if( param != ""){
        fMultiFitter->fCompareLimits = Common::StringToBoolean(param);
    }

    // Ser ComparePOI
    param = confSet->Get("ComparePOI");
    if( param != ""){
        fMultiFitter->fComparePOI = Common::StringToBoolean(param);
    }

    // Ser CompareEFT
    param = confSet->Get("CompareEFT");
    if( param != ""){
        fMultiFitter->fCompareEFT = Common::StringToBoolean(param);
    }

    // Set ComparePulls
    param = confSet->Get("ComparePulls");
    if( param != "" ){
        fMultiFitter->fComparePulls = Common::StringToBoolean(param);
    }

    // Set PlotCombCorrMatrix
    param = confSet->Get("PlotCombCorrMatrix");
    if (param != ""){
        fMultiFitter->fPlotCombCorrMatrix = Common::StringToBoolean(param);
    }

    // Set CorrelationThreshold
    param = confSet->Get("CorrelationThreshold");
    if( param != ""){
        TRExFitter::CORRELATIONTHRESHOLD = atof(param.c_str());
    }

    // Set Combine
    param = confSet->Get("Combine");
    if( param != ""){
        fMultiFitter->fCombine = Common::StringToBoolean(param);
    }

    // Set Compare
    param = confSet->Get("Compare");
    if( param != ""){
        fMultiFitter->fCompare = Common::StringToBoolean(param);
    }

    // Set StatOnly
    param = confSet->Get("StatOnly");
    if( param != "" ){
        fMultiFitter->fStatOnly = Common::StringToBoolean(param);
        LOG(WARNING) << "StatOnlyFit option is deprecated for obtaining stat-only uncertainty. Please, use the covariance matrix breakdown\n";
    }

    // Set IncludeStatOnly
    param = confSet->Get("IncludeStatOnly");
    if( param != ""){
        fMultiFitter->fIncludeStatOnly = Common::StringToBoolean(param);
        if (fMultiFitter->fIncludeStatOnly) {
            fMultiFitter->fUnfoldingShowStat = true;
        }
    }

    // Set POILabel
    param = confSet->Get("POILabel");
    if (param != "") {
        const auto vec = Common::Vectorize(param, ',');
        if (vec.size () != fMultiFitter->fPOIs.size()) {
            LOG(ERROR) << "You specified 'POILabel' option but you didn't pass the same number of parameters as the number of POIs. Please check this!\n";
            ++sc;
        }
        fMultiFitter->fPOIName = vec;
     } else {
        for (std::size_t i = 0; i < fMultiFitter->fPOIs.size(); ++i) {
            fMultiFitter->fPOIName.emplace_back("#mu");
        }
     }

    // Set POINominal
    param = confSet->Get("POINominal");
    if (param != "") {
        const auto vec = Common::Vectorize(param, ',');
        if (vec.size () != fMultiFitter->fPOIs.size()) {
            LOG(ERROR) << "You specified 'POINominal' option but you didn't pass the same number of parameters as the number of POIs. Please check this!\n";
            ++sc;
        }
        for (const auto& ivec : vec) {
            fMultiFitter->fPOINominal.emplace_back(Common::convertStoNum<float>(ivec));
        }
     } else {
        for (std::size_t i = 0; i < fMultiFitter->fPOIs.size(); ++i) {
            fMultiFitter->fPOINominal.emplace_back(1);
        }
     }

    // Set POIRange
    param = confSet->Get("POIRange");
    if( param != ""){
        const auto tmp = Common::Vectorize(param, ',');
        if (tmp.size() != fMultiFitter->fPOIs.size()) {
            LOG(ERROR) << "You specified 'POIRange' option but you didn't pass the same number of parametrs as the number of POIs. Please check this!\n";
            ++sc;
        }
        for (const auto& irange : tmp) {
            const auto vec = Common::Vectorize(irange,':');
            if (vec.size()==2 ) {
                fMultiFitter->fPOIMin.emplace_back(atof( vec[0].c_str() ));
                fMultiFitter->fPOIMax.emplace_back(atof( vec[1].c_str() ));
            } else {
                LOG(ERROR) << "You specified 'POIRange' option but you didn't provide valid setting. Please check this!\n";
                ++sc;
            }
        }
    }

    // Set LimitMax
    param = confSet->Get("LimitMax");
    if( param != "" ) {
        fMultiFitter->fLimitMax = atof( param.c_str() );
    }

    // Set POIPrecision
    param = confSet->Get("POIPrecision");
    if( param != "" ) {
        const auto tmp = Common::Vectorize(param, ',');
        if (tmp.size() != fMultiFitter->fPOIs.size()) {
            LOG(ERROR) << "You specified 'POIPrecision' option but you didn't pass the same number of parametrs as the number of POIs. Please check this!\n";
            ++sc;
        }
        for (const auto& i : tmp) {
            fMultiFitter->fPOIPrecision.emplace_back(Common::RemoveQuotes(i).c_str());
        }
    }

    //Set DataName
    param = confSet->Get("DataName");
    if(param == "") {
        LOG(ERROR) << "You did not set \"DataName\" in the config. This option is needed to decide which data is used in the combination\n";
        LOG(ERROR) << "You can set it to \"asimovData\" to use pure asimov dataset or \"obsData\" to use real data\n";
        LOG(ERROR) << "In case you want to use the mixed dataset, you need to use \"mixedAsimovData\", for custom asimov \"customAsimov\"\n";
        LOG(ERROR) << "In both of these special cases, you need to set the \"Workspace\" option in the \"Fit\" block to point to the right workspaces with the specal datasets\n";
        LOG(ERROR) << "For the mixed dataset, the WS ROOT file in the \"RooStats/\" folder has \"_realisticAsimov\" suffix, and for the custom asimov it has the \"_customAsimov\" suffix\n";
        ++sc;
    } else {
        fMultiFitter->fDataName = Common::RemoveQuotes(param);
    }

    // Set FitType
    param = confSet->Get("FitType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param=="SPLUSB")      fMultiFitter->fFitType = 1;
        else if(param=="BONLY")  fMultiFitter->fFitType = 2;
        else if(param=="UNFOLDING")  fMultiFitter->fFitType = 3;
        else {
            LOG(WARNING) << "You specified 'FitType' option but you didn't provide valid setting. Using default (SPLUSB)\n";
            fMultiFitter->fFitType = 1;
        }
    }

    // Set NPCategories
    param = confSet->Get("NPCategories");
    if( param != "" ) {
        std::vector<std::string> categ = Common::Vectorize(param,',');
        for(const std::string &icat : categ)
            fMultiFitter->fNPCategories.push_back(icat);
    }

    // Set SetRandomInitialNPval
    param = confSet->Get("SetRandomInitialNPval");
    if( param != ""){
        fMultiFitter->fUseRnd = true;
        fMultiFitter->fRndRange = atof(param.c_str());
    }

    // Set SetRandomInitialNPvalSeed
    param = confSet->Get("SetRandomInitialNPvalSeed");
    if( param != ""){
        fMultiFitter->fRndSeed = atol(param.c_str());
    }

    // Set NumCPU
    param = confSet->Get("NumCPU");
    if( param != "" ){
        fMultiFitter->fCPU = atoi( param.c_str());
        if (fMultiFitter->fCPU <= 0) {
            LOG(ERROR) << "NumCPU must be a positive integer!\n";
            ++sc;
        }
    }

    // Set FastFit
    param = confSet->Get("FastFit");
    if (param != ""){
        fMultiFitter->fFastFit = Common::StringToBoolean(param);
    }

    // Set FastFitForRanking
    param = confSet->Get("FastFitForRanking");
    if (param != ""){
        fMultiFitter->fFastFitForRanking = Common::StringToBoolean(param);
    }

    // Set NuisParListFile
    param = confSet->Get("NuisParListFile");
    if( param != "" ) fMultiFitter->fNuisParListFile = Common::RemoveQuotes(param);

    // Set PlotSoverB
    param = confSet->Get("PlotSoverB");
    if (param != ""){
        fMultiFitter->fPlotSoverB = Common::StringToBoolean(param);
    }

    // Set SignalTitle
    param = confSet->Get("SignalTitle");
    if( param != "" ) fMultiFitter->fSignalTitle = Common::RemoveQuotes(param);

    // Set FitResultsFile
    param = confSet->Get("FitResultsRootFile");
    if( param != "" ) fMultiFitter->fFitResultsRootFile = Common::RemoveQuotes(param);

    // Set LimitsFile
    param = confSet->Get("LimitsFile");
    if( param != "" ) fMultiFitter->fLimitsFile = Common::RemoveQuotes(param);

    // Set BonlySuffix
    param = confSet->Get("BonlySuffix");
    if( param != "" ) fMultiFitter->fBonlySuffix = Common::RemoveQuotes(param);

    // Set ShowSystForPOI
    param = confSet->Get("ShowSystForPOI");
    if( param != ""){
        fMultiFitter->fShowSystForPOI = Common::StringToBoolean(param);
    }

    // Set doLHscan
    param = confSet->Get("doLHscan");
    if( fRunLHscan && param != "" ){
        if (fLHscanVariable==""){
            fMultiFitter->fVarNameLH = Common::Vectorize(param,',');
        } else {
            fMultiFitter->fVarNameLH.emplace_back(fLHscanVariable);
        }
    }

    // Set do2DLHscan
    param = confSet->Get("do2DLHscan");
    if( param != "" ){
        const std::vector<std::string> tmp = Common::Vectorize(param,':');
        for (const auto& ivec : tmp){
            const std::vector<std::string> v = Common::Vectorize(ivec,',');
            if (v.size() == 2){
                fMultiFitter->fVarName2DLH.emplace_back(v);
            } else {
                LOG(WARNING) << "You specified 'do2DLHscan' option but did not provide correct input. Ignoring\n";
            }
        }
    }

    // Set LHscanMin
    param = confSet->Get("LHscanMin");
    if ( fRunLHscan && param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMin' option but did not set doLHscan. Ignoring\n";
        } else {
            fMultiFitter->fLHscanMin = Common::convertStoNum<float>(param);
        }
    }

    // Set LHscanMax
    param = confSet->Get("LHscanMax");
    if ( fRunLHscan && param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMax' option but did not set doLHscan. Ignoring\n";
        } else {
            fMultiFitter->fLHscanMax = Common::convertStoNum<float>(param);
        }
    }
    // Set LHscanMin for second variable
    param = confSet->Get("LHscanMinY");
    if ( param != "" ) {
        if (fMultiFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMinY' option but did not set do2DLHscan. Ignoring\n";
        } else {
            fMultiFitter->fLHscanMinY = Common::convertStoNum<float>(param);
        }
    }

    // Set LHscanMax for second variable
    param = confSet->Get("LHscanMaxY");
    if ( param != "" ) {
        if (fMultiFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanMaxY' option but did not set do2DLHscan. Ignoring\n";
        } else {
            fMultiFitter->fLHscanMaxY = Common::convertStoNum<float>(param);
        }
    }


    // Set LHscanSteps
    param = confSet->Get("LHscanSteps");
    if ( param != "" ) {
        if (fMultiFitter->fVarNameLH.size() == 0 && fMultiFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanSteps' option but did not set doLHscan. Ignoring\n";
        } else {
            fMultiFitter->fLHscanSteps = Common::convertStoNum<int>(param);
            if(fMultiFitter->fLHscanSteps < 3 || fMultiFitter->fLHscanSteps > 100){
                LOG(WARNING) << "LHscanSteps is smaller than 3 or larger than 100, setting to defaut (30)\n";
                fMultiFitter->fLHscanSteps = 30;
            }
        }
    }

    // Set LHscanSteps for second variable
    param = confSet->Get("LHscanStepsY");
    if ( param != "" ) {
        if (fMultiFitter->fVarName2DLH.size() == 0){
            LOG(WARNING) << "You specified 'LHscanStepsY' option but did not set do2DLHscan. Ignoring\n";
        } else {
            fMultiFitter->fLHscanStepsY = Common::convertStoNum<int>(param);
            if(fMultiFitter->fLHscanStepsY < 3 || fMultiFitter->fLHscanStepsY > 100){
                LOG(WARNING) << "LHscanStepsY is smaller than 3 or larger than 100, setting to LHscanSteps\n";
                fMultiFitter->fLHscanStepsY = fMultiFitter->fLHscanSteps;
            }
        }
    }
    else {
        fMultiFitter->fLHscanStepsY = fMultiFitter->fLHscanSteps;
    }



    // Set LHscanStep (-1 = process all points, otherwise process the single point only)
    if (fLHscanStep != -1) {
        if (fLHscanStep < -1 || fLHscanStep >= fMultiFitter->fLHscanSteps) {
            ++sc;
            LOG(ERROR) << "Invalid option for LHscanStep, please check!\n";
        }
        fMultiFitter->fLHscanStep = fLHscanStep;
    }

    if (fLHscanStepY != -1) {
        if (fLHscanStepY < -1 || fLHscanStepY >= fMultiFitter->fLHscanStepsY) {
            ++sc;
            LOG(ERROR) << "Invalid command line option for LHscanStepY, please check!\n";
        }
        fMultiFitter->fLHscanStepY = fLHscanStepY;
    }

    // Set ShowTotalOnly
    param = confSet->Get("ShowTotalOnly");
    if ( param != "" ) {
        fMultiFitter->fShowTotalOnly = Common::StringToBoolean(param);
    }

    // Set PlotOptions
    param = confSet->Get("PlotOptions");
    if( param != ""){
        std::vector<std::string> vec = Common::Vectorize(param,',');
        if( std::find(vec.begin(), vec.end(), "PREFITONPOSTFIT")   !=vec.end() )  TRExFitter::PREFITONPOSTFIT= true;
    }

    // Set POIAsimovl
    param = confSet->Get("POIAsimov");
    if( param != "" ) {
        const auto tmp = Common::Vectorize(param, ',');
        if (tmp.size() != fMultiFitter->fPOIs.size()) {
            LOG(ERROR) << "You specified 'POIAsimov' option but you didn't pass the same number of parametrs as the number of POIs. Please check this!\n";
            ++sc;
        }
        for (const auto& i : tmp) {
            fMultiFitter->fPOIAsimov.emplace_back(Common::convertStoNum<float>(i));
        }
    } else {
        for (std::size_t i = 0; i < fMultiFitter->fPOIs.size(); ++i) {
            fMultiFitter->fPOIAsimov.emplace_back(1);
        }
    }

    // Set POIInitial
    param = confSet->Get("POIInitial");
    if( param != ""){
        const std::vector<std::string> tmp = Common::Vectorize(param, ',');
        if (tmp.size() != fMultiFitter->fPOIs.size()) {
            LOG(ERROR) << "Number of POI initial values does not match the number of POIs\n";
            return 1;
        }

        for (std::size_t i = 0; i < tmp.size(); ++i) {
            double value(0);
            try {
                value = Common::convertStoNum<float>(tmp.at(i));
            } catch (...) {
                LOG(ERROR) << "Cannot convert POIInitial element to float\n";
                return 1;
            }
            fMultiFitter->fPOIInitials.emplace_back(fMultiFitter->fPOIs.at(i), value);
        }
    }

    // Set HEPDataFormat
    param = confSet->Get("HEPDataFormat");
    if( param != ""){
        fMultiFitter->fHEPDataFormat = Common::StringToBoolean(param);
    }

    // Set UheppFormat
    param = confSet->Get("UheppFormat");
    if( param != ""){
        fMultiFitter->fUheppFormat = Common::StringToBoolean(param);
    }

    // Set FitStrategy
    param = confSet->Get("FitStrategy");
    if (param != "") {
        fMultiFitter->fFitStrategy = Common::convertStoNum<int>(param);
        if (fMultiFitter->fFitStrategy > 3) {
            LOG(WARNING) << "FitStrategy > 3, setting to default (-1)\n";
            fMultiFitter->fFitStrategy = -1;
        }
    }

    // Set BinnedLikelihoodOptimization
    param = confSet->Get("BinnedLikelihoodOptimization");
    if (param != "") {
        fMultiFitter->fBinnedLikelihood = Common::StringToBoolean(param);
    }

    // Set UseHesseBeforeMigrad
    param = confSet->Get("UseHesseBeforeMigrad");
    if (param != "") {
        fMultiFitter->fUseHesseBeforeMigrad = Common::StringToBoolean(param);
    }

    // Set UsePOISinRanking
    param = confSet->Get("UsePOISinRanking");
    if (param != "") {
        fMultiFitter->fUsePOISinRanking = Common::StringToBoolean(param);
    }

    // Set Regions
    param = confSet->Get("Regions");
    if (param != "") {
        fMultiFitter->fOnlyRegions = Common::Vectorize(param,',');
    }

    param = confSet->Get("MultiStageFit");
    if (param != "") {
        fMultiFitter->fSpeedUpFit = Common::StringToBoolean(param);
    }

    param = confSet->Get("NPCutOff");
    if (param != "") {
        fMultiFitter->fNPCutOff = Common::convertStoNum<float>(param);
        if (fMultiFitter->fNPCutOff < 0. || fMultiFitter->fNPCutOff > 1.0) {
            LOG(WARNING) << "NPCutOff is < 0 or > 1, this is not a valid option, setting to default (0.3)\n";
            fMultiFitter->fNPCutOff = 0.3;
        }
    }

    param = confSet->Get("PlotUnfoldedData");
    if (param != "") {
        fMultiFitter->fPlotUnfoldedData = Common::StringToBoolean(param);
        if (fMultiFitter->fPlotUnfoldedData && fMultiFitter->fFitType != 3) {
            LOG(WARNING) << "PlotUnfoldedData set to TRUE, but the FitType is not UNFOLDING, ignoring\n";
            fMultiFitter->fPlotUnfoldedData = false;
        }
    }

    param = confSet->Get("UnfoldingShowStatOnlyError");
    if (param != "") {
        fMultiFitter->fUnfoldingShowStat = Common::StringToBoolean(param);
    }

    param = confSet->Get("UnfoldingShowUncertaintyBreakdown");
    if (param != "") {
        fMultiFitter->fUnfoldingShowUncertaintyBreakdown = Common::StringToBoolean(param);
    }

    param = confSet->Get("UnfoldingUncertaintyBreakdownTotal");
    if (param != "") {
        fMultiFitter->fUnfoldingUncertaintyBreakdownTotal = Common::StringToBoolean(param);
    }

    param = confSet->Get("PlotGammas");
    if (param != "") {
        fMultiFitter->fPlotGammas = Common::StringToBoolean(param);
    }

    param = confSet->Get("BlindedParameters");
    if( param != "" ){
        const std::vector<std::string> tmp = Common::Vectorize( param,',');
        fMultiFitter->fBlindedParams = std::move(tmp);
    }

    param = confSet->Get("RankingPOIName");
    if( param != ""){
        fMultiFitter->fRankingPOIName = Common::Vectorize(param, ',');
    }

    param = confSet->Get("ToleranceScale");
    if( param != ""){
        fMultiFitter->fToleranceScale = Common::convertStoNum<double>(param);
    }

    param = confSet->Get("ShiftGlobalObservablesInRanking");
    if (param != ""){
        fMultiFitter->fShiftGlobalObservablesInRanking = Common::StringToBoolean(param);
    }

    // Set RankingUpperAxisNdivision
    param = confSet->Get("RankingUpperAxisNdivision");
    if( param != ""){
        const auto tmp = Common::Vectorize(param, ',');
        std::vector<int> v;
        for (const auto& i : tmp) {
            v.emplace_back(Common::convertStoNum<int>(i));
        }
        fMultiFitter->fRankingUpperAxisNdivision = v;
    }

    // Set RankingPOIAxisScale
    param = confSet->Get("RankingPOIAxisScale");
    if( param != ""){
        const auto tmp = Common::Vectorize(param, ',');
        std::vector<double> v;
        for (const auto& i : tmp) {
            v.emplace_back(Common::convertStoNum<float>(i));
        }
        fMultiFitter->fRankingPOIAxisScale = v;
    }

    param = confSet->Get("MaximumNumberFCNcalls");
    if (param != "") {
        const int max = Common::convertStoNum<int>(param);
        if (max < 0) {
            LOG(WARNING) << "\"MaximumNumberFCNcalls\" is < 0. Ignoring the option\n";
        } else {
            fMultiFitter->fMaximumNumberFCNcalls = max;
        }
    }

    param = confSet->Get("UseAutoDiff");
    if (param != ""){
        fMultiFitter->fUseAutoDiff = Common::StringToBoolean(param);
    }

    param = confSet->Get("NLLOffset");
    if (param != "") {
        // Ensure lowercase to match RooFit documentation
        std::transform(param.begin(), param.end(), param.begin(), ::tolower);
        fMultiFitter->fNLLOffset = param;
    }

    param = confSet->Get("FitToys");
    if (param != ""){
        fMultiFitter->fFitToys = Common::convertStoNum<int>(param);
    }

    param = confSet->Get("ToysSeed");
    if (param != ""){
        fMultiFitter->fToysSeed = Common::convertStoNum<int>(param);
    }

    param = confSet->Get("ToysHistoNbins");
    if (param != ""){
        fMultiFitter->fToysHistoNbins = Common::convertStoNum<int>(param);
    }

    param = confSet->Get("ToysOnlyStatFluctuation");
    if (param != ""){
        fMultiFitter->fToysSetStatOnlyFluctuation = Common::StringToBoolean(param);
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

            auto itr = fMultiFitter->fToysRandomPOIStartingValues.find(nfName);
            if (itr != fMultiFitter->fToysRandomPOIStartingValues.end()) {
                LOG(WARNING) << nfName << " already set in the config for \"FitToysRandomInitialNormFactorValues\", ignoring\n";
            } else {
                fMultiFitter->fToysRandomPOIStartingValues.insert({nfName, std::make_pair(min, max)});
            }
        }
    }

    param = confSet->Get("ToysLHScanForAll");
    if (param != ""){
        fMultiFitter->fToysLHScanForAll = Common::StringToBoolean(param);
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

            fMultiFitter->fToysLHScanCondition.insert({nfName, std::make_pair(min, max)});
        }
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReaderMulti::ReadLimitOptions() {

    int sc(0);

    const ConfigSet* confSet = fParser.GetConfigSet("Limit");
    if (confSet == nullptr){
        LOG(DEBUG) << "You do not have Limit option in the config. It is ok, we just want to let you know.\n";
        return 0; // it is ok to not have Fit set up
    }

    // Set POI
    std::string param = confSet->Get("POI");
    if( param != "" ){
        fMultiFitter->fPOIforLimit = Common::RemoveQuotes(param);
        fMultiFitter->AddPOI(Common::RemoveQuotes(param));
    }

    // Set LimitType
    param = confSet->Get("LimitType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if( param == "ASYMPTOTIC" ){
            fMultiFitter->SetLimitType(TRExFit::LimitType::ASYMPTOTIC);
        }
        else if( param == "TOYS" ){
            fMultiFitter->SetLimitType(TRExFit::LimitType::TOYS);
        }
        else{
            LOG(ERROR) << "Unknown LimitType argument : " << confSet->Get("LimitType") << "\n";
            ++sc;
        }
    }

    // Set LimitBlind
    param = confSet->Get("LimitBlind");
    if( param != "" ){
        fMultiFitter->fLimitIsBlind = Common::StringToBoolean(param);
    }

    // Set SignalInjection
    param = confSet->Get("SignalInjection");
    if( param != "" ){
        fMultiFitter->fSignalInjection = Common::StringToBoolean(param);
    }

    // Set SignalInjectionValue
    param = confSet->Get("SignalInjectionValue");
    if( param != "" ){
        fMultiFitter->fSignalInjectionValue = Common::convertStoNum<float>(param);
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fMultiFitter->fLimitParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fMultiFitter->fLimitParamValue = Common::convertStoNum<float>(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fMultiFitter->fLimitOutputPrefixName = param;
    }

    param = confSet->Get("ConfidenceLevel");
    if( param != "" ){
        double conf = Common::convertStoNum<double>(param);
        if (conf <= 0. || conf >= 1.){
            LOG(WARNING) << "Confidence level is <= 0 or >=1. Setting to default 0.95\n";
            conf = 0.95;
        }
        fMultiFitter->fLimitsConfidence = conf;
    }

    param = confSet->Get("SplusBToys");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "SplusBToys is available only when TOYS are used\n";
        } else {
            fMultiFitter->fLimitToysStepsSplusB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("BonlyToys");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "BonlyToys is available only when TOYS are used\n";
        } else {
            fMultiFitter->fLimitToysStepsB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("ScanSteps");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ScanSteps is available only when TOYS are used\n";
        } else {
            fMultiFitter->fLimitToysScanSteps = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("ScanMin");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ScanMin is available only when TOYS are used\n";
        } else {
            fMultiFitter->fLimitToysScanMin = Common::convertStoNum<float>(param);
        }
    }

    param = confSet->Get("ScanMax");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ScanMax is available only when TOYS are used\n";
        } else {
            fMultiFitter->fLimitToysScanMax = Common::convertStoNum<float>(param);
        }
    }

    param = confSet->Get("LimitPlot");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "LimitPlot is available only when TOYS are used\n";
        } else {
            fMultiFitter->fLimitPlot = Common::StringToBoolean(param);
        }
    }

    param = confSet->Get("LimitFile");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "LimitFile is available only when TOYS are used\n";
        } else {
            fMultiFitter->fLimitFile = Common::StringToBoolean(param);
        }
    }

    param = confSet->Get("FitStrategy");
    if( param != "" ) {
        if (fMultiFitter->fLimitType != TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "FitStrategy is available only when ASYMPTOTIC limits are used\n";
        } else {
            fMultiFitter->fLimitFitStrategy = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("ToysSeed");
    if (param != ""){
        fMultiFitter->fLimitToysSeed = Common::convertStoNum<int>(param);
    }

    param = confSet->Get("UseNLLwithoutOffsetInLHscan");
    if (param != "") {
        fMultiFitter->fUseNllInLHscan = Common::StringToBoolean(param);
    }

    param = confSet->Get("UseHesse");
    if (param != "") {
        fMultiFitter->fUseHesse = Common::StringToBoolean(param);
    }

    param = confSet->Get("ToysUsexRooFit");
    if( param != "" ) {
        if (fMultiFitter->fLimitType == TRExFit::LimitType::ASYMPTOTIC) {
            LOG(WARNING) << "ToysUsexRooFit is available only when TOYS are used\n";
        } else {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if (param == "FALSE") {
                fMultiFitter->fLimitToysUsexRooFit = TRExFit::ToysUsexRooFit::FALSE;
            } else if (param == "TRUE") {
                fMultiFitter->fLimitToysUsexRooFit = TRExFit::ToysUsexRooFit::TRUE;
            } else if (param == "AUTO") {
                fMultiFitter->fLimitToysUsexRooFit = TRExFit::ToysUsexRooFit::AUTO;
            } else {
                LOG(ERROR) << "Unknown ToysUsexRooFit argument : " << confSet->Get("ToysUsexRooFit") << "\n";
                ++sc;
            }
        }
    }

    param = confSet->Get("UseAutoDiff");
    if( param != "" ){
        fMultiFitter->fLimitUseAutoDiff = Common::StringToBoolean(param);
    }

    return sc;
}

//__________________________________________________________________________________
//
int ConfigReaderMulti::ReadSignificanceOptions() {

    int sc(0);

    const ConfigSet* confSet = fParser.GetConfigSet("Significance");
    if (confSet == nullptr){
        LOG(DEBUG) << "You do not have Significance option in the config. It is ok, we just want to let you know.\n";
        return 0; // it is ok to not have Fit set up
    }

    // Set POI
    std::string param = confSet->Get("POI");
    if( param != "" ){
        fMultiFitter->fPOIforSig = Common::RemoveQuotes(param);
        fMultiFitter->AddPOI(Common::RemoveQuotes(param));
    }

    // Set SignificanceType
    param = confSet->Get("SignificanceType");
    if( param != "" ){
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if(param == "ASYMPTOTIC"){
            fMultiFitter->SetSignificanceType(TRExFit::SignificanceType::ASYMPTOTIC);
        }
        else if(param == "TOYS"){
            fMultiFitter->SetSignificanceType(TRExFit::SignificanceType::TOYS);
        }
        else{
            LOG(ERROR) << "Unknown SignificanceType argument : " << confSet->Get("SignificanceType") << "\n";
            ++sc;
        }
    }

    param = confSet->Get("SplusBToys");
    if( param != "" ) {
        if (fMultiFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "SplusBToys is available only when TOYS are used\n";
        } else {
            fMultiFitter->fSignificanceToysStepsSplusB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("BonlyToys");
    if( param != "" ) {
        if (fMultiFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "BonlyToys is available only when TOYS are used\n";
        } else {
            fMultiFitter->fSignificanceToysStepsB = Common::convertStoNum<int>(param);
        }
    }

    param = confSet->Get("SignificancePlot");
    if( param != "" ) {
        if (fMultiFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "SignificancePlot is available only when TOYS are used\n";
        } else {
            fMultiFitter->fSignificancePlot = Common::StringToBoolean(param);
        }
    }

    // Set SignificanceBlind
    param = confSet->Get("SignificanceBlind");
    if( param != "" ){
        fMultiFitter->fSignificanceIsBlind = Common::StringToBoolean(param);
    }

    // Set POIAsimov
    param = confSet->Get("POIAsimov");
    if( param != "" ){
        fMultiFitter->fSignificancePOIAsimov = atof(param.c_str());
        fMultiFitter->fSignificanceDoInj = true;
    }

    param = confSet->Get("ParamName");
    if( param != "" ){
        fMultiFitter->fSignificanceParamName = param;
    }

    param = confSet->Get("ParamValue");
    if( param != "" ){
        fMultiFitter->fSignificanceParamValue = Common::convertStoNum<float>(param);
    }

    param = confSet->Get("OutputPrefixName");
    if( param != "" ){
        fMultiFitter->fSignificanceOutputPrefixName = param;
    }

    param = confSet->Get("ToysSeed");
    if (param != ""){
        fMultiFitter->fSignificanceToysSeed = Common::convertStoNum<int>(param);
    }

    param = confSet->Get("ToysUsexRooFit");
    if( param != "") {
        if (fMultiFitter->fSignificanceType == TRExFit::SignificanceType::ASYMPTOTIC) {
            LOG(WARNING) << "ToysUsexRooFit is available only when TOYS are used\n";
        } else {
            fMultiFitter->fSignificanceToysUsexRooFit = Common::StringToBoolean(param);
        }
    }

    param = confSet->Get("UseAutoDiff");
    if( param != "" ){
        fMultiFitter->fSignificanceUseAutoDiff = Common::StringToBoolean(param);
    }

    param = confSet->Get("FitStrategy");
    if( param != "" ){
        fMultiFitter->fSignificanceFitStrategy = Common::convertStoNum<int>(param);
    }

    return sc;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::ReadFitOptions(const std::string& opt, const std::string& options){

    int sc(0);

    const auto* fit_confsets = fParser.GetConfigSets("Fit");
    if (fit_confsets == nullptr) {
        LOG(ERROR) << "You need to provide at least one 'Fit' option. Please check this!\n";
        return ++sc;
    }

    for(const auto& confSetPtr : *fit_confsets) {
        const ConfigSet *confSet = confSetPtr.get();
        if (confSet == nullptr) break;

        // Set Options
        std::string fullOptions;
        std::string param = confSet->Get("Options");
        if(param!="" && options!="") fullOptions = options+";"+Common::RemoveQuotes(param);
        else if(param!="") fullOptions = Common::RemoveQuotes(param);
        else fullOptions = options;

        // name
        fMultiFitter->fFitNames.push_back(Common::CheckName(confSet->GetValue()));

        // Set Label
        param = confSet->Get("Label");
        std::string label = Common::CheckName(confSet->GetValue());
        if(param!="") label = Common::RemoveQuotes(param);

        // Set suf
        param = confSet->Get("LoadSuf");
        const std::string loadSuf = param != "" ? Common::RemoveQuotes(param) : fGlobalSuffix;

        // config file
        std::string confFile = "";
        param = confSet->Get("ConfigFile");
        if(param!="") confFile = Common::RemoveQuotes(param);

        // workspace
        std::string wsFile = "";
        param = confSet->Get("Workspace");
        if(param!="") wsFile = Common::RemoveQuotes(param);

        bool useInFit(true);
        param = confSet->Get("UseInFit");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") {
                useInFit = false;
            }
        }

        bool useInComparison(true);
        param = confSet->Get("UseInComparison");
        if (param != "") {
            std::transform(param.begin(), param.end(), param.begin(), ::toupper);
            if(param=="FALSE") {
                useInComparison = false;
            }
        }


        fMultiFitter->AddFitFromConfig(confFile, opt, fullOptions, label, loadSuf, wsFile, useInFit, useInComparison);

        // Set FitResultsFile
        param = confSet->Get("FitResultsRootFile");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fFitResultsRootFile = Common::RemoveQuotes(param);
        fMultiFitter->fLimitsFiles.push_back("");

        // Set LimitsFile
        param = confSet->Get("LimitsFile");
        if( param != "" ) fMultiFitter->fLimitsFiles[fMultiFitter->fFitList.size()-1] = Common::RemoveQuotes(param);

        // Set POIName
        param = confSet->Get("POIName");
        if (param != "") {
            const std::vector<std::string> vec = Common::Vectorize(param, ',');
            for (const auto& i : vec) {
                fMultiFitter->AddPOI(i);
            }
        }

        // Set Directory
        param = confSet->Get("Directory");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fName = Common::RemoveQuotes(param);

        // Set InputName
        param = confSet->Get("InputName");
        if( param != "" ) fMultiFitter->fFitList[fMultiFitter->fFitList.size()-1]->fInputName = Common::RemoveQuotes(param);
    }

    return sc;
}

//_______________________________________________________________________________________
//
int ConfigReaderMulti::PostConfig() {

    int sc(0);

    if (fMultiFitter->fRankingPOIName.empty()) {
        for (std::size_t i = 0; i < fMultiFitter->fPOIs.size(); ++i) {
            fMultiFitter->fRankingPOIName.emplace_back("#mu");
        }
    } else if (fMultiFitter->fRankingPOIName.size() != fMultiFitter->fPOIs.size()) {
        LOG(ERROR) << "Sizes of RankingPOIName and POIs do not match\n";
        ++sc;
    }

    if (fMultiFitter->fRankingUpperAxisNdivision.empty()) {
        for (std::size_t i = 0; i < fMultiFitter->fPOIs.size(); ++i) {
            fMultiFitter->fRankingUpperAxisNdivision.emplace_back(510);
        }
    } else if (fMultiFitter->fRankingUpperAxisNdivision.size() != fMultiFitter->fPOIs.size()) {
        LOG(ERROR) << "Sizes of RankingUpperAxisNdivision and POIs do not match\n";
        ++sc;
    }

    if (fMultiFitter->fRankingPOIAxisScale.empty()) {
        for (std::size_t i = 0; i < fMultiFitter->fPOIs.size(); ++i) {
            fMultiFitter->fRankingPOIAxisScale.emplace_back(1.0);
        }
    } else if (fMultiFitter->fRankingPOIAxisScale.size() != fMultiFitter->fPOIs.size()) {
        LOG(ERROR) << "Sizes of RankingPOIAxisScale and POIs do not match\n";
        ++sc;
    }

    return sc;
}
