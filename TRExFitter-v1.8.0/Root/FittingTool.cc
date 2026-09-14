// Class include
#include "TRExFitter/FittingTool.h"

//Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/YamlConverter.h"

// Blinder
#include "Blinder/Blinder/Blinder.h"

//ROOR includes
#include "TCanvas.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFile.h"

//Roostats includes
#include "Math/MinimizerOptions.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ModelConfig.h"

//Roofit includes
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooArgSet.h"

//c++ includes
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

//________________________________________________________________________
//
FittingTool::FittingTool():
    m_CPU(1),
    m_useMinos(false),
    m_constPOI(false),
    m_fitResult(nullptr),
    m_noGammas(false),
    m_noSystematics(false),
    m_randomize(false),
    m_randomNP(0.1),
    m_randSeed(-999),
    m_minNll(999999),
    m_externalConstraints(nullptr),
    m_strategy(-1),
    m_useHesse(true),
    m_hesseBeforeMigrad(false),
    m_speedUpFit(false),
    m_errorSigma(1.),
    m_maxFCNcalls(-1),
    m_toleranceScale(1.),
    m_useAutoDiff(false),
    m_nllOffset("initial")
{
}

//________________________________________________________________________
//
void FittingTool::SetSubCategories() {
    LOG(DEBUG) << "Finding unique SubCategories\n";
    // loop over m_subCategoryMap to find all unique SubCategories, save in m_subCategories set
    for (const auto& it : m_subCategoryMap) {
        m_subCategories.insert(it.second);
    }
}

//________________________________________________________________________
//
void FittingTool::AddValPOI(const std::string& name, const double value) {
    auto it = std::find_if(m_valPOIs.begin(), m_valPOIs.end(),
    [&name](const std::pair<std::string, double>& element){ return element.first == name;});

    // take only unique ones
    if (it != m_valPOIs.end()) return;

    m_valPOIs.emplace_back(std::make_pair(name, value));
}

//____________________________________________________________________________________
//
void FittingTool::ReplacePOIVal(const std::string& name, const double value) {
    auto it = std::find_if(m_valPOIs.begin(), m_valPOIs.end(),
        [&](const std::pair<std::string, double>& element){return name == element.first;});

    if (it == m_valPOIs.end()) return;

    it->second = value;
}

//________________________________________________________________________
//
void FittingTool::SetUseAutoDiff(const bool flag) {
    m_useAutoDiff = flag;
    if (flag) {
        LOG(INFO) << "Will use automatic differentiation for fitting\n";
    }
}

//________________________________________________________________________
//
void FittingTool::SetNLLOffset(const std::string& val) {
    m_nllOffset = val;
    LOG(INFO) << "Likelihood offset set to: " << m_nllOffset << "\n";
}

//________________________________________________________________________
//
double FittingTool::FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, bool fastFit, bool noFit) {

    LOG(DEBUG) << "-> Entering in FitPDF function\n";

    //
    // Printing the whole model for information
    //
    if(TRExFitter::DEBUGLEVEL >= 2) model->Print();

    //
    // Get the global observables (nominal values)
    //
    const RooArgSet* glbObs = model->GetGlobalObservables();

    //
    // Create the likelihood based on fitpdf, fitData and the parameters
    //
    const RooArgSet dummy{};
    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (m_useAutoDiff) {
        nll = std::unique_ptr<RooAbsReal>(fitpdf->createNLL(*fitdata,
                                          RooFit::GlobalObservables(*glbObs),
                                          RooFit::Optimize(kTRUE),
                                          RooFit::ExternalConstraints(m_externalConstraints ? *m_externalConstraints : dummy),
                                          RooFit::EvalBackend("codegen")
                                          ));
    } else {
        nll = std::unique_ptr<RooAbsReal>(fitpdf->createNLL(*fitdata,
                                          RooFit::GlobalObservables(*glbObs),
					  RooFit::Offset(m_nllOffset),
                                          RooFit::NumCPU(m_CPU, RooFit::Hybrid),
                                          RooFit::Optimize(kTRUE),
                                          RooFit::ExternalConstraints(m_externalConstraints ? *m_externalConstraints : dummy)
                                          ));

    }

    //
    // Needed for Ranking plot, but also to set random initial values for the NPs
    //
    if(m_randSeed == -999){
        gRandom->SetSeed(time(nullptr));
    }
    else{
        gRandom->SetSeed(m_randSeed);
    }

    //
    // Getting the POIs
    //
    FitUtils::FittingOptions options;
    options.noGammas = m_noGammas;
    options.noSystematics = m_noSystematics;
    options.randomize = m_randomize;
    options.randomValue = m_randomNP;
    options.constNPs = m_constNP;
    options.constNPvalues = m_constNPvalue;
    options.initialNPs = m_initialNP;
    options.initialNPvalues = m_initialNPvalue;
    options.constPOI = m_constPOI;
    options.poiVals = m_valPOIs;
    options.constNFs = m_constNFs;

    FitUtils::SetStartingAndConstParams(model, options);

    double nllval = nll->getVal();

    LOG(DEBUG) << "   -> Initial value of the NLL = " << nllval << "\n";

    //
    // return here if specified not to perform the fit
    if(noFit) {
        m_minNll = nllval;
        return nllval;
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // run the fit to spped up thing if requested
    if (m_speedUpFit) {
        SpeedUpFit(nll.get(), model);
    }

    // Safe fit loop
    int strat = ::ROOT::Math::MinimizerOptions::DefaultStrategy();
    if(m_strategy >= 0) {
        strat = m_strategy;
        LOG(INFO) << "Manually setting strategy to " << strat << "\n";
    }

    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this

    const std::size_t nparams = FittingTool::GetNumberOfNPsPlusPOIs(model);

    int printLevel = TRExFitter::DEBUGLEVEL - 1;
    if (TRExFitter::DEBUGLEVEL == 2) printLevel = 0;

    RooMinimizer minim(*nll);
    minim.optimizeConst(2);
    minim.setStrategy(strat);
    minim.setMinimizerType(minimType.Data());
    minim.setPrintLevel(printLevel);
    minim.setEps(m_toleranceScale*tol);
    minim.setErrorLevel(0.5*m_errorSigma*m_errorSigma);
    if (m_maxFCNcalls < 0) {
        if (nparams > 100) {
            minim.setMaxFunctionCalls(5*nparams*nparams);
        }
    } else {
        minim.setMaxFunctionCalls(m_maxFCNcalls);
    }
    // it doesn't make sense to try more than 2 additional strategies
    const int maxRetries = 3 - strat;

    // fast fit - e.g. for ranking
    if(fastFit){
        minim.setStrategy(0);  // to be the same as ttH comb
        minim.setPrintLevel(0);
    }

    // always run one fit
    LOG(INFO) << "\n";
    LOG(INFO) << "\n";
    LOG(INFO) << "\n";
    LOG(INFO) << "Fit try no. 1\n";
    LOG(INFO) << "======================\n";
    LOG(INFO) << "\n";

    int status = minim.minimize(minimType.Data(),algorithm.Data());
    if (m_useHesse) {
        if (status == 0 || m_hesseBeforeMigrad) {
            minim.hesse();
        }
    }
    std::unique_ptr<RooFitResult> r(minim.save());
    if (!r) {
        LOG(ERROR) << "Cannot get minimisation status, exiting\n";
        exit(EXIT_FAILURE);
    }
    double edm = r->edm();
    status = r->status();

    // check if the fit converged
    bool fitIsNotGood = (status > 1) || (status < 0) || (edm > 0.0001*m_toleranceScale);

    // if not, loop maximum 2 times with increased strategy
    int nrItr = 0;
    while (nrItr<maxRetries && fitIsNotGood){
        LOG(WARNING) << "\n";
        LOG(WARNING) << "   *******************************\n";
        LOG(WARNING) << "   * Increasing Minuit strategy (was " << strat << ")\n";
        ++strat;
        LOG(WARNING) << "   * Fit failed with : \n";
        LOG(WARNING) << "      - minuit status " << (status % 100) << "\n";
        LOG(WARNING) << "      - hess status " << (status / 100) << "\n";
        LOG(WARNING) << "      - Edm = " << edm << "\n";
        LOG(WARNING) << "   * Retrying with strategy " << strat << "\n";
        LOG(WARNING) << "   ********************************\n";
        LOG(WARNING) << "\n";
        PrintMinuitHelp();
        minim.setStrategy(strat);
        status = minim.minimize(minimType.Data(),algorithm.Data());
        if (m_useHesse) {
            if (status <= 1 || m_hesseBeforeMigrad) {
                minim.hesse();
            }
        }
        r = std::unique_ptr<RooFitResult>(minim.save());
        if (!r) {
            LOG(ERROR) << "Cannot get minimisation status, exiting\n";
            exit(EXIT_FAILURE);
        }
        edm = r->edm();
        status = r->status();

        fitIsNotGood = (status > 1) || (status < 0) || (edm > 0.0001*m_toleranceScale);
        nrItr++;
    }

    // if the fit is not good even after retries print an error message
    if (fitIsNotGood) {
        LOG(ERROR) << "\n";
        LOG(ERROR) << "\n";
        LOG(ERROR) << "\n";
        LOG(ERROR) << "\n";
        LOG(ERROR) << "***********************************************************\n";
        LOG(ERROR) << "Fit failure unresolved with status " << status << "\n";
        LOG(ERROR) << "   Please investigate your workspace\n";
        LOG(ERROR) << "   Find a wall : you will need it to crash your head on it\n";
        LOG(ERROR) << "***********************************************************\n";
        LOG(ERROR) << "\n";
        LOG(ERROR) << "\n";
        LOG(ERROR) << "\n";
        PrintMinuitHelp();
        m_fitResult = nullptr;

        return 0;
    }

    int minosStatus = -1000;
    if(m_useMinos){
        // check if minos is run on all parameters
        if (!m_varMinos.empty() && m_varMinos.at(0) == "all") {
            minosStatus = minim.minos();
        } else {
            RooArgSet minosSet;
            // check POIs
            for (const auto ipoi : *model->GetParametersOfInterest()) {
                const std::string name = ipoi->GetName();
                auto itr = std::find(m_varMinos.begin(), m_varMinos.end(), name);
                if (itr != m_varMinos.end()) {
                    minosSet.add(*ipoi);
                }
            }
            // check NPs if they exist
            if (model->GetNuisanceParameters()) {
                for (const auto inp : *model->GetNuisanceParameters()) {
                    const std::string name = inp->GetName();

                    // name starts with "alpha_" for systematics so we need to check the substring match
                    auto itr = std::find_if(m_varMinos.begin(), m_varMinos.end(), [&name](const auto& element){return name.find(element) != std::string::npos;});
                    if (itr != m_varMinos.end()) {
                        minosSet.add(*inp);
                    }
                }
            }
            minosStatus = minim.minos(minosSet);
        }
    }

    r = std::unique_ptr<RooFitResult>(minim.save());
    LOG(INFO) << "\n";
    LOG(INFO) << "\n";
    LOG(INFO) << "\n";
    LOG(INFO) << "***********************************************************\n";
    LOG(INFO) << "         FIT FINALIZED SUCCESSFULLY : \n";
    LOG(INFO) << "            - minuit status " << (status % 100) << "\n";
    LOG(INFO) << "            - hess status " << (status / 100) << "\n";
    LOG(INFO) << "            - Edm = " << edm << "\n";
    if (minosStatus != -1000) LOG(INFO) << "            - minos status " << minosStatus << "\n";
    if (status == 1) LOG(WARNING) << "            Covariance was made pos-def!\n";
    LOG(INFO) << "***********************************************************\n";
    LOG(INFO) << "\n";
    LOG(INFO) << "\n";
    LOG(INFO) << "\n";

    if(r!=nullptr) m_fitResult = std::unique_ptr<RooFitResult>(static_cast<RooFitResult*>(r->Clone()));

    //
    // clean stuff

    nllval = nll->getVal();

    if(TRExFitter::DEBUGLEVEL >= 1) {
        std::streamsize ss = std::cout.precision();
        std::cout << std::fixed << std::setprecision(20);

        LOG(INFO) << "***********************************************************\n";
        LOG(INFO) << "  Final value of the NLL = " << nllval << "\n";
        LOG(INFO) << "***********************************************************\n";

        std::cout << resetiosflags( ios::fixed | ios::showpoint );
        std::cout << std::setprecision(ss);
    }
    m_minNll = nllval;
    return nllval;
}

//____________________________________________________________________________________
//
void FittingTool::SaveFitResult( const std::string &fileName )
{
    std::unique_ptr<TFile> f(TFile::Open(fileName.c_str(),"RECREATE"));
    if (!f) {
        LOG(WARNING) << "Cannot open file: " << fileName << " cannot write the fit result\n";
        return;
    }
    m_fitResult->Write("",TObject::kOverwrite);
    f->Close();
}

//____________________________________________________________________________________
//
void FittingTool::ExportFitResultInTextFile( const std::string &fileName, const std::vector<std::string>& blinded )
{
    if(!m_fitResult){
        LOG(ERROR) << "The FitResultObject seems not to be defined.\n";
        return;
    }

    //
    // Also saves fit result in root file with same name as txt file
    //
    TString fName = fileName;
    fName.ReplaceAll(".txt",".root");
    SaveFitResult(fName.Data());

    //
    // Printing the nuisance parameters post-fit values
    //
    ofstream nuisParAndCorr(fileName);
    nuisParAndCorr << "NUISANCE_PARAMETERS\n";

    for (auto var_tmp : m_fitResult->floatParsFinal()) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);

        // Not consider nuisance parameter being not associated to syst (yet)
        TString vname=var->GetName();
        vname.ReplaceAll("alpha_","");

        const double pull  = var->getVal(); // GetValue() return value in unit of sigma
        const double errorHi = var->getErrorHi();
        const double errorLo = var->getErrorLo();

        if (blinded.size() == 0){
            FittingTool::CheckUnderconstraint(var);
            nuisParAndCorr << vname << "  " << pull << " +" << std::abs(errorHi) << " -" << std::abs(errorLo)  << "\n";
        } else {
            std::string vname_s = vname.Data();
            if (std::find(blinded.begin(), blinded.end(), vname_s) == blinded.end()){
                FittingTool::CheckUnderconstraint(var);
                nuisParAndCorr << vname << "  " << pull << " +" << std::abs(errorHi) << " -" << std::abs(errorLo)  << "\n";
            } else {
                const std::string& hex = Blinder::DoubleToPseudoHex(pull);
                nuisParAndCorr << vname << "  " << hex << " +" << std::abs(errorHi) << " -" << std::abs(errorLo)  << "\n";
            }
        }
    }

    //
    // Correlation matrix
    //
    std::unique_ptr<TH2> h2Dcorrelation(m_fitResult->correlationHist());
    nuisParAndCorr << "\n\nCORRELATION_MATRIX\n";
    nuisParAndCorr << h2Dcorrelation->GetNbinsX() << "   " << h2Dcorrelation->GetNbinsY() << "\n";
    for(int kk=1; kk < h2Dcorrelation->GetNbinsX()+1; kk++) {
        for(int ll=1; ll < h2Dcorrelation->GetNbinsY()+1; ll++) {
            nuisParAndCorr << h2Dcorrelation->GetBinContent(kk,ll) << "   ";
        }
        nuisParAndCorr << "\n";
    }

    //
    // NLL value
    //
    nuisParAndCorr << "\n\nNLL\n";
    nuisParAndCorr << std::fixed << std::setprecision(5);
    nuisParAndCorr << m_minNll << "\n";

    //
    // Closing the output file
    //
    nuisParAndCorr << "\n";
    nuisParAndCorr.close();
}

//____________________________________________________________________________________
//
std::map <std::string, double> FittingTool::ExportFitResultInMap(){
    std::map <std::string, double> result;
    if(!m_fitResult){
        LOG(ERROR) << "The FitResultObject seems not to be defined.\n";
        return result;
    }
    for (auto var_tmp : m_fitResult->floatParsFinal()) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
        // Not consider nuisance parameter being not associated to syst
        std::string varname = var->GetName();
        const double pull  = var->getVal();
        result.insert( std::pair < std::string, double >(varname, pull) );
    }
    return result;
}

//____________________________________________________________________________________
//
void FittingTool::GetGroupedImpact(RooStats::ModelConfig* model,
                                   RooAbsPdf* fitpdf,
                                   RooAbsData* fitdata,
                                   RooWorkspace* ws,
                                   const std::string& categoryOfInterest,
                                   const std::string& outFileName,
                                   const std::string& fileName,
                                   const std::string& lumiLabel,
                                   const std::string& cmeLabel,
                                   const bool useHEPData) const {

    std::vector<RooRealVar*> pois = FitUtils::GetVectorPOI(model);
    std::vector<std::ofstream> outFiles(pois.size());

    if (pois.empty()) {
        LOG(ERROR) << "No POIs found!\n";
        return;
    }

    if(m_randSeed == -999){
        gRandom->SetSeed(time(nullptr));
    }
    else{
        gRandom->SetSeed(m_randSeed);
    }

    std::size_t iPOI(0);
    for (auto poi : pois) {
        poi -> setConstant(false);
        double value(0);
        const std::string name = poi->GetName();
        auto it = std::find_if(m_valPOIs.begin(), m_valPOIs.end(),
            [&name](const std::pair<std::string, double>& element){ return element.first == name;});

        if (it != m_valPOIs.end()) {
            value = it->second;
        }
        poi->setVal(value);
        if(!m_constPOI && m_randomize){
            poi->setVal(value + m_randomNP*(gRandom->Uniform(2)-1.) );
        }
        outFiles.at(iPOI).open((outFileName+"_"+name+".txt").c_str());
        ++iPOI;
    }

    // save snapshot of original workspace
    ws->saveSnapshot("snapshot_AfterFit_POI", *(model->GetParametersOfInterest()) );
    if (model->GetNuisanceParameters()) ws->saveSnapshot("snapshot_AfterFit_NP" , *(model->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(model->GetGlobalObservables())    );

    std::vector<std::string> associatedParams; // parameters associated to a SubCategory

    // repeat the nominal fit - done so that the initial randomization is the exact same as for the following fit(s)
    // this should help avoid issues with fits ending up in different local minima for groups with very small impact on the POI
    FitExcludingGroup(false, false, fitdata, fitpdf, model, ws, "Nominal", associatedParams);  // nothing held constant -> "snapshot_AfterFit_POI_Nominal"

    // loop over unique SubCategories
    for (const auto& cat : m_subCategories) {
        if(categoryOfInterest!="all" && cat != categoryOfInterest) continue; // if a category was specified via command line, only process that one

        LOG(INFO) << "Performing grouped systematics impact evaluation for: " << cat << "\n";

        // find all associated parameters per SubCategory
        associatedParams.clear();
        for(const auto& itSysts : m_subCategoryMap) {
            if (itSysts.second == cat) {
                associatedParams.push_back(itSysts.first);
            }
        }

        // special case for gammas
        if(cat == "Gammas") {
            FitExcludingGroup(true,  false, fitdata, fitpdf, model, ws, cat, associatedParams);
        }

        // special case for stat-only fit
        else if(cat == "FullSyst") {
            FitExcludingGroup(true,  true,  fitdata, fitpdf, model, ws, cat, associatedParams);
        }

        // default: perform a fit where parameters in SubCategory are held constant
        else {
            FitExcludingGroup(false, false, fitdata, fitpdf, model, ws, cat, associatedParams);
        }
    }

    // load original workspace again
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    if (model->GetNuisanceParameters()) {
        ws->loadSnapshot("snapshot_AfterFit_NP");
    }

    LOG(INFO) << "-----------------------------------------------------\n";

    for (std::size_t ipoi(0); ipoi < pois.size(); ++ipoi) {
        const std::string& name = pois.at(ipoi)->GetName();
        // report replication of nominal fit
        ws->loadSnapshot("snapshot_AfterFit_POI_Nominal");
        const double NomUp2=(pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi());
        const double NomLo2=(pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo());
        const double Nom   = 0.5*std::abs(pois.at(ipoi)->getErrorHi() - pois.at(ipoi)->getErrorLo());
        const double Nom2  =Nom*Nom;
        LOG(INFO) << "replicated nominal fit\n";
        LOG(INFO) << "POI: " << name << "\n";
        LOG(INFO) << "     " << pois.at(ipoi)->getVal() << " +/- " << Nom  <<
                                                         "    ( +" << pois.at(ipoi)->getErrorHi() << ", " << pois.at(ipoi)->getErrorLo() << " )\n";
        std::vector<YamlConverter::ImpactContainer> container;

        // report impact calculations, impact is obtained by quadrature subtraction from replicated nominal fit
        for (const auto& cat : m_subCategories) {
            if(categoryOfInterest!="all" && cat != categoryOfInterest) continue; // if a category was specified via command line, only process that one

            ws->loadSnapshot(("snapshot_AfterFit_POI_" + cat).c_str());
            const double errUp = pois.at(ipoi)->getErrorHi();
            const double errDown = pois.at(ipoi)->getErrorLo();
            const double errorAv = 0.5 * std::abs(errUp - errDown);
            LOG(INFO) << "-----------------------------------------------------\n";
            LOG(INFO) <<  "category: " << cat + " (fixed to best-fit values for fit)\n";
            LOG(INFO) <<  "POI is:   " << pois.at(ipoi)->getVal() << " +/- " << errorAv << "    ( +" << pois.at(ipoi)->getErrorHi() << ", " << pois.at(ipoi)->getErrorLo() << " )\n";
            if (cat == "FullSyst") LOG(DEBUG) << "  (corresponds to a stat-only fit)\n";
            const double impact = -(errorAv*errorAv) + Nom2;
            if (impact < 0) {
                LOG(WARNING) << "NaN impact encountered\n";
                LOG(WARNING) << "  - Nominal error squared: " << Nom2 << "\n";
                LOG(WARNING) << "  - Category fixed error squared: " << (errorAv*errorAv) << "\n";
                LOG(WARNING) << "  - Nominal error squared minus category error squared: " << impact << "\n";
            }
            LOG(INFO) << "           --> impact: " << std::sqrt(impact) << "    ( +" << std::sqrt(- (pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi()) + NomUp2) << ", -" << std::sqrt(-(pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo()) + NomLo2) << " )\n";
            // write results to file
            outFiles.at(ipoi) << cat << "    " << std::sqrt(impact) << "  ( +" << std::sqrt( - (pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi()) + NomUp2 ) << ", -" << std::sqrt( - (pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo()) + NomLo2 ) << " )\n";

            YamlConverter::ImpactContainer imp;
            imp.name    = cat;
            imp.error   = std::sqrt(impact);
            imp.errorHi = std::sqrt(-(pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi()) + NomUp2);
            imp.errorLo = std::sqrt(-(pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo()) + NomLo2);
            container.emplace_back(std::move(imp));
        }

        LOG(INFO) << "-----------------------------------------------------\n";

        YamlConverter converter{};
        converter.WriteImpact(container, outFileName+"_"+name+".yaml");
        if (useHEPData) {
            converter.SetLumi(Common::ReplaceString(lumiLabel, " fb^{-1}", ""));
            converter.SetCME(Common::ReplaceString(cmeLabel, " TeV", "000"));
            converter.WriteImpactHEPData(container, fileName, "_"+name);
        }

        outFiles.at(ipoi).close();
    }

    // load original workspace again
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    if (model->GetNuisanceParameters()) {
        ws->loadSnapshot("snapshot_AfterFit_NP");
    }
}

//____________________________________________________________________________________
//
// perform a fit where all parameters in "affectedParams" (usually coming from SubGroup "category") are set to constant, optionally also gammas or all parameters
void FittingTool::FitExcludingGroup(bool excludeGammas,
                                    bool statOnly,
                                    RooAbsData*& fitdata,
                                    RooAbsPdf*& fitpdf,
                                    RooStats::ModelConfig* mc,
                                    RooWorkspace* ws,
                                    const std::string& category,
                                    const std::vector<std::string>& affectedParams) const {

    if (!mc->GetNuisanceParameters()) return;
    std::vector<RooRealVar*> pois = FitUtils::GetVectorPOI(mc);
    // (VD): use this to fix nuisance parameter before the fit
    const RooArgSet* glbObs = mc->GetGlobalObservables();
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    ws->loadSnapshot("snapshot_AfterFit_NP");

    LOG(INFO) << "-----------------------------------------------------\n";
    LOG(INFO) << "           breakdown for " << category << "\n";
    LOG(INFO) << "-----------------------------------------------------\n";

    for (auto var_tmp : *mc->GetNuisanceParameters()) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
        std::string varname = var->GetName();

        // default: set everything non-constant
        var->setConstant(0);

        // if excludeGammas==true, set gammas to constant
        if (excludeGammas) {
            if (varname.find("gamma_stat")!=string::npos) var->setConstant(1);
        }

        // set all affectedParams constant
        if (std::find(affectedParams.begin(), affectedParams.end(), varname) != affectedParams.end()) {
            var->setConstant(1);
        } else {
            if (m_randomize && !statOnly && varname.find("alpha_") != std::string::npos) {
                var->setVal(m_randomNP*(gRandom->Uniform(2) - 1.));
            }
        }

        // for stat-only fits, set everything constant
        if (statOnly) {
            var->setConstant(1);
        }
    }

    // repeat the fit here ....
    const RooArgSet dummy{};
    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (m_useAutoDiff) {
        nll = std::unique_ptr<RooAbsReal>(fitpdf->createNLL(*fitdata,
                                          RooFit::GlobalObservables(*glbObs),
                                          RooFit::Optimize(kTRUE),
                                          RooFit::Offset(1),
                                          RooFit::NumCPU(m_CPU, RooFit::Hybrid),
                                          RooFit::ExternalConstraints(m_externalConstraints ? *m_externalConstraints : dummy),
                                          RooFit::EvalBackend("codegen")
                                          ));
    } else {
        nll = std::unique_ptr<RooAbsReal>(fitpdf->createNLL(*fitdata,
                                          RooFit::GlobalObservables(*glbObs),
                                          RooFit::Offset(1),
                                          RooFit::NumCPU(m_CPU, RooFit::Hybrid),
                                          RooFit::Optimize(kTRUE),
                                          RooFit::ExternalConstraints(m_externalConstraints ? *m_externalConstraints : dummy)
                                          ));
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this
    const std::size_t nparams = FittingTool::GetNumberOfNPsPlusPOIs(mc);

    RooMinimizer minim2(*nll);
    minim2.optimizeConst(2);
    minim2.setStrategy(m_strategy);
    minim2.setMinimizerType(minimType.Data());
    minim2.setPrintLevel(0);
    minim2.setEps(tol*m_toleranceScale);
    if (m_maxFCNcalls < 0) {
        if (nparams > 100) {
            minim2.setMaxFunctionCalls(5*nparams*nparams);
        }
    } else {
        minim2.setMaxFunctionCalls(m_maxFCNcalls);
    }
    const int status = minim2.minimize(minimType.Data(),algorithm.Data());
    const bool hessStatus = minim2.hesse();
    if (status != 0 || hessStatus) {
        LOG(ERROR) << "Unable to perform fit correctly!\n";
        LOG(ERROR) << "\tMigrad status: " << status << "\n";
        LOG(ERROR) << "\tHesse status: " << hessStatus << "\n";
    }

    if (m_useMinos) {
        RooArgSet minosSet;
        for (const auto ipoi : pois) {
            minosSet.add(*ipoi);
        }
        minim2.minos(minosSet);
    }

    bool first(true);
    for (auto poi : pois) {

        const double newPOIerrU=poi->getErrorHi();
        const double newPOIerrD=poi->getErrorLo();
        const double newPOIerr = 0.5*(std::abs(newPOIerrU - newPOIerrD));

        const std::string snapshotName = "snapshot_AfterFit_POI_" + category;
        if (first) ws->saveSnapshot(snapshotName.c_str(), *mc->GetParametersOfInterest() );

        ws->loadSnapshot("snapshot_AfterFit_POI_Nominal");
        const double oldPOIerrU=poi->getErrorHi();
        const double oldPOIerrD=poi->getErrorLo();
        const double oldPOIerr = 0.5*(oldPOIerrU - oldPOIerrD);

        // check if uncertainties have increased compared to nominal fit, with 0.5% tolerance
        if ( (std::abs(newPOIerrU)>std::abs(oldPOIerrU)*1.005) || (std::abs(newPOIerrD)>std::abs(oldPOIerrD)*1.005) ) {
            const std::string& name = poi->GetName();
            LOG(WARNING) << "POI: " << name << "\n";
            LOG(WARNING) << "uncertainty has increased for " << category << "! please check the fit\n";
            LOG(WARNING) << "old: " << oldPOIerr <<  " (+" << oldPOIerrU <<  ", "  << oldPOIerrD <<  ")\n";
            LOG(WARNING) << "new: " << newPOIerr <<  " (+" << newPOIerrU <<  ", "  << newPOIerrD <<  ")\n";
        }
        first = false;
    }
}

//____________________________________________________________________________________
//
// Check for underconstraints
void FittingTool::CheckUnderconstraint(const RooRealVar* const var) const {
    const std::string name = var->GetName();
    const double errorHi = var->getErrorHi();
    const double errorLo = var->getErrorLo();

    // dont check gamma parameters
    if (name.find("alpha_") == std::string::npos) return;

    if (errorHi > 1.001*m_errorSigma || errorLo < -1.001*m_errorSigma){
        LOG(WARNING) << "NuisanceParameter: " << name << " is underconstrained! This may indicate fit convergence problems!\n";
    }
}

//__________________________________________________________________________________
//
void FittingTool::PrintMinuitHelp() const {
    if (TRExFitter::DEBUGLEVEL < 2) return;
    LOG(INFO) << "\n";
    LOG(INFO) << "------------------------------------------\n";
    LOG(INFO) << "         Status code cheat-sheet          \n";
    LOG(INFO) << "------------------------------------------\n";
    LOG(INFO) << " status = MIGRAD status + 100*Hesse status\n";
    LOG(INFO) << " MIGRAD status\n";
    LOG(INFO) << "    - 0 = No problems\n";
    LOG(INFO) << "    - 1 = Covariance was made pos defined\n";
    LOG(INFO) << "    - 2 = Hesse is invalid\n";
    LOG(INFO) << "    - 3 = Edm is above max\n";
    LOG(INFO) << "    - 4 = Reached call limit\n";
    LOG(INFO) << "    - 5 = Any other failure\n";
    LOG(INFO) << " Hesse status\n";
    LOG(INFO) << "    - 0 = No problems\n";
    LOG(INFO) << "    - 1 = Hesse failed\n";
    LOG(INFO) << "    - 2 = Matrix inversion failed\n";
    LOG(INFO) << "    - 3 = Matrix is not pos defined\n";
    LOG(INFO) << "\n";
}

//__________________________________________________________________________________
//
void FittingTool::SpeedUpFit(RooAbsReal* nll, RooStats::ModelConfig* mc) const {
    if (mc->GetNuisanceParameters() == nullptr) return;

    LOG(INFO) << "Running the speed up fit\n";

    // disable gammas
    for (auto ivar : *mc->GetNuisanceParameters()) {
        RooRealVar* var = static_cast<RooRealVar*>(ivar);
        const std::string name = var->GetName();
        if (name.find("gamma_") != std::string::npos) {
            var->setVal(1.);
            var->setConstant(true);
        }
    }

    for (auto ivar : *mc->GetNuisanceParameters()) {
        RooRealVar* var = static_cast<RooRealVar*>(ivar);
        const std::string name = var->GetName();
        for (const auto& ismall : m_smallNPs) {
            if (name == ismall.first || name == "alpha_"+ismall.first) {
                var->setConstant(true);
            }
        }
    }

    // run the fit
    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance();

    int printLevel = TRExFitter::DEBUGLEVEL - 1;
    if (TRExFitter::DEBUGLEVEL == 2) printLevel = 0;

    RooMinimizer minim(*nll);
    minim.optimizeConst(2);
    minim.setStrategy(0);
    minim.setMinimizerType(minimType.Data());
    minim.setPrintLevel(printLevel);
    minim.setEps(tol*m_toleranceScale);

    minim.minimize(minimType.Data(),algorithm.Data());
    std::unique_ptr<RooFitResult> rfr(minim.save());
    const double edm = rfr->edm();
    const int status = rfr->status();

    if (edm > 0.0001*m_toleranceScale || status > 1) {
        LOG(WARNING) << "Speed up fit failed. Setting ALL NPs to 0 and all gammas to 1 and let them free float\n";
        // will just set everything back and bail out
        for (auto ivar : *mc->GetNuisanceParameters()) {
            RooRealVar* var = static_cast<RooRealVar*>(ivar);
            const std::string name = var->GetName();
            if (name.find("gamma_") != std::string::npos) {
                var->setVal(1.);
                var->setConstant(false);
            } else if (name.find("alpha_") != std::string::npos){
                var->setVal(0);
                var->setConstant(false);
            } else {
                var->setVal(1.);
                var->setConstant(false);
            }
        }
        return;
    }

    // enable gammas
    for (auto ivar : *mc->GetNuisanceParameters()) {
        RooRealVar* var = static_cast<RooRealVar*>(ivar);
        const std::string name = var->GetName();
        if (name.find("gamma_") != std::string::npos) {
            var->setVal(1.);
            var->setConstant(false);
       }
    }

    // propagate the fitted values now
    auto args = *mc->GetNuisanceParameters();
    auto args2 = *mc->GetParametersOfInterest();
    for (const auto& arg : args2) {
        args.add(*arg);
    }

    LOG(DEBUG) << "Will propagate these values\n";
    for (auto ivar : rfr->floatParsFinal()) {
        RooRealVar* var = static_cast<RooRealVar*>(ivar);
        const TString name = var->GetName();
        const double pull = var->getVal();
        RooRealVar* arg = static_cast<RooRealVar*>(&args[name]);
        arg->setVal(pull);
        LOG(DEBUG) << "\t" << name << ": value: " << pull << "\n";
    }
}


//__________________________________________________________________________________
//
std::size_t FittingTool::GetNumberOfNPsPlusPOIs(RooStats::ModelConfig* mc) {
    const std::size_t pois = mc->GetParametersOfInterest()->size();
    std::size_t nps(0);
    if (mc->GetNuisanceParameters()) {
        nps += mc->GetNuisanceParameters()->size();
    }

    return nps + pois;
}

//__________________________________________________________________________________
//
std::map<std::string,double> FittingTool::CalculateErrorDecomposition(RooStats::ModelConfig* mc, const std::string& poiName){
    // get fit result
    auto fr = this->GetFitResult();
    // check if fit-result is there
    if (!fr) {
        LOG(WARNING) << "No fit results found. Skipping error decomposition.\n";
        return m_errDecompMap[poiName];
    }
    //
    double totError = 0.;
    double systErrSq = 0.;
    double gammaErrSq = 0.;
    double nfErrSq = 0.;
    double shapeFactorSq = 0;
    //
    double totErrorUp = 0.;
    double systErrSqUp = 0.;
    double gammaErrSqUp = 0.;
    double nfErrSqUp = 0.;
    double shapeFactorSqUp = 0;
    //
    double totErrorDown = 0.;
    double systErrSqDown = 0.;
    double gammaErrSqDown = 0.;
    double nfErrSqDown = 0.;
    double shapeFactorSqDown = 0;
    //
    // store list of POIs
    std::vector<std::string> poiList;
    for (auto poi : *mc->GetParametersOfInterest()) {
        poiList.emplace_back(poi->GetName());
    }
    // check if specified POI is in the list: otherwise, throw a WARNING and skip error decomposition
    if (Common::FindInStringVector(poiList,poiName)<0) {
        LOG(WARNING) << "Specified POI not found in ModelConfig. Skipping error decomposition.\n";
        return m_errDecompMap[poiName];
    }
    //
    // store total error on POI
    for (auto arg : fr->floatParsFinal()) {
        RooRealVar *par = static_cast<RooRealVar*>(arg);
        std::string parName = par->GetName();
        if (parName == poiName) {
            totError = par->getError();
            totErrorUp = par->getErrorHi();
            totErrorDown = par->getErrorLo();
            break;
        }
    }
    bool hasSF(false);
    //
    // loop on all parameters in the fit-result
    for (auto arg : fr->floatParsFinal()) {
        RooRealVar *par = static_cast<RooRealVar*>(arg);
        const std::string parName = par->GetName();
        // skip POIs
        if (Common::FindInStringVector(poiList,parName)>=0) continue;
        // for other paramaters, evaluate impact
        const double corr = fr->correlation(parName.c_str(),poiName.c_str());
        const double err = par->getError();
        const double errUp = par->getErrorHi();
        const double errDown = par->getErrorLo();
        double impact = corr*err*totError;
        double impactUp = corr*errUp*totErrorUp;
        double impactDown = -corr*errDown*totErrorDown;
        if (Common::StartsWith(parName,"gamma_stat") || Common::StartsWith(parName,"gamma_shape")) {
            const double prefit = this->GetPrefitUncertainty(fr, parName, 0);
            const double prefitUp = this->GetPrefitUncertainty(fr, parName, 1);
            const double prefitDown = this->GetPrefitUncertainty(fr, parName, -1);
            impact /= prefit;
            impactUp /= prefitUp;
            impactDown /= prefitDown;
        }
        m_errDecompMap[poiName][parName] = impact;
        m_errDecompMapUp[poiName][parName] = impactUp;
        m_errDecompMapDown[poiName][parName] = impactDown;
        const double impactSq = impact*impact;
        // for asymmetric impacts, sum impacts to "up" or "down" depending on sign
        double impactSqUp = 0.;
        double impactSqDown = 0.;
        if (impactUp>0) impactSqUp += impactUp*impactUp;
        if (impactDown>0) impactSqUp += impactDown*impactDown;
        if (impactUp<0) impactSqDown += impactUp*impactUp;
        if (impactDown<0) impactSqDown += impactDown*impactDown;
        // if it's an "alpha" NP -> add impact to total systematic
        if (Common::StartsWith(parName,"alpha_") ) {
            systErrSq += impactSq;
            systErrSqUp += impactSqUp;
            systErrSqDown+= impactSqDown;
        }
        // if gamma_stat or gamma_shape (separate gammas) -> add to total gamma
        else if (Common::StartsWith(parName,"gamma_stat") || Common::StartsWith(parName,"gamma_shape")) {
            gammaErrSq += impactSq;
            gammaErrSqUp += impactSqUp;
            gammaErrSqDown += impactSqDown;
        }
        // if still starts with "gamma_" it is a shape factor
        else if (Common::StartsWith(parName,"gamma_")) {
            shapeFactorSq += impactSq;
            shapeFactorSqUp += impactSqUp;
            shapeFactorSqDown += impactSqDown;
            hasSF = true;
        }
        // otherwise it's a free parameter -> add impact to total norm-factor error
        else {
            nfErrSq += impactSq;
            nfErrSqUp += impactSqUp;
            nfErrSqDown += impactSqDown;
        }
    }
    //
    const double totAlphaError = std::sqrt(systErrSq);
    const double totGammaError = std::sqrt(gammaErrSq);
    const double totSystError = std::sqrt(systErrSq+gammaErrSq);
    const double totNfError = std::sqrt(nfErrSq);
    const double statError = std::sqrt(totError*totError - systErrSq - gammaErrSq);
    const double statErrorNoNF = std::sqrt(statError*statError - nfErrSq);
    const double statErrorNoSF = std::sqrt(statError*statError - shapeFactorSq);
    const double statErrorNoNFNoSF = std::sqrt(statError*statError - nfErrSq - shapeFactorSq);
    const double shapeFactorError = std::sqrt(shapeFactorSq);
    m_errDecompMap[poiName]["TOT_ERROR"] = totError;
    m_errDecompMap[poiName]["STAT_ERROR"] = statError;
    m_errDecompMap[poiName]["SYST_ERROR"] = totAlphaError;
    m_errDecompMap[poiName]["MCSTAT_ERROR"] = totGammaError;
    m_errDecompMap[poiName]["TOTSYST_ERROR"] = totSystError;
    m_errDecompMap[poiName]["STAT_ERR_NO_NF"] = statErrorNoNF;
    m_errDecompMap[poiName]["NF_ERROR"] = totNfError;
    if (hasSF) {
        m_errDecompMap[poiName]["SHAPE_FACTOR_ERROR"] = shapeFactorError;
        m_errDecompMap[poiName]["STAT_ERR_NO_SF"] = statErrorNoSF;
        m_errDecompMap[poiName]["STAT_ERR_NO_NF_NO_SF"] = statErrorNoNFNoSF;
    }
    //
    const double totAlphaErrorUp = std::sqrt(systErrSqUp);
    const double totGammaErrorUp = std::sqrt(gammaErrSqUp);
    const double totSystErrorUp = std::sqrt(systErrSqUp+gammaErrSqUp);
    const double totNfErrorUp = std::sqrt(nfErrSqUp);
    const double statErrorUp = std::sqrt(totErrorUp*totErrorUp - systErrSqUp - gammaErrSqUp);
    const double statErrorNoNFUp = std::sqrt(statErrorUp*statErrorUp - nfErrSqUp);
    const double statErrorNoSFUp = std::sqrt(statErrorUp*statErrorUp - shapeFactorSqUp);
    const double statErrorNoNFNoSFUp = std::sqrt(statErrorUp*statErrorUp - nfErrSqUp - shapeFactorSqUp);
    const double shapeFactorErrorUp = std::sqrt(shapeFactorSqUp);
    m_errDecompMapUp[poiName]["TOT_ERROR"] = totErrorUp;
    m_errDecompMapUp[poiName]["STAT_ERROR"] = statErrorUp;
    m_errDecompMapUp[poiName]["SYST_ERROR"] = totAlphaErrorUp;
    m_errDecompMapUp[poiName]["MCSTAT_ERROR"] = totGammaErrorUp;
    m_errDecompMapUp[poiName]["TOTSYST_ERROR"] = totSystErrorUp;
    m_errDecompMapUp[poiName]["STAT_ERR_NO_NF"] = statErrorNoNFUp;
    m_errDecompMapUp[poiName]["NF_ERROR"] = totNfErrorUp;
    if (hasSF) {
        m_errDecompMapUp[poiName]["SHAPE_FACTOR_ERROR"] = shapeFactorErrorUp;
        m_errDecompMapUp[poiName]["STAT_ERR_NO_SF"] = statErrorNoSFUp;
        m_errDecompMapUp[poiName]["STAT_ERR_NO_NF_NO_SF"] = statErrorNoNFNoSFUp;
    }
    //
    const double totAlphaErrorDown = std::sqrt(systErrSqDown);
    const double totGammaErrorDown = std::sqrt(gammaErrSqDown);
    const double totSystErrorDown = std::sqrt(systErrSqDown+gammaErrSqDown);
    const double totNfErrorDown = std::sqrt(nfErrSqDown);
    const double statErrorDown = std::sqrt(totErrorDown*totErrorDown - systErrSqDown - gammaErrSqDown);
    const double statErrorNoNFDown = std::sqrt(statErrorDown*statErrorDown - nfErrSqDown);
    const double statErrorNoSFDown = std::sqrt(statErrorDown*statErrorDown - shapeFactorSqDown);
    const double statErrorNoNFNoSFDown = std::sqrt(statErrorDown*statErrorDown - nfErrSqDown - shapeFactorSqDown);
    const double shapeFactorErrorDown = std::sqrt(shapeFactorSqDown);
    m_errDecompMapDown[poiName]["TOT_ERROR"] = totErrorDown;
    m_errDecompMapDown[poiName]["STAT_ERROR"] = -statErrorDown;
    m_errDecompMapDown[poiName]["SYST_ERROR"] = -totAlphaErrorDown;
    m_errDecompMapDown[poiName]["MCSTAT_ERROR"] = -totGammaErrorDown;
    m_errDecompMapDown[poiName]["TOTSYST_ERROR"] = -totSystErrorDown;
    m_errDecompMapDown[poiName]["STAT_ERR_NO_NF"] = -statErrorNoNFDown;
    m_errDecompMapDown[poiName]["NF_ERROR"] = -totNfErrorDown;
    if (hasSF) {
        m_errDecompMapDown[poiName]["SHAPE_FACTOR_ERROR"] = -shapeFactorErrorDown;
        m_errDecompMapDown[poiName]["STAT_ERR_NO_SF"] = -statErrorNoSFDown;
        m_errDecompMapDown[poiName]["STAT_ERR_NO_NF_NO_SF"] = -statErrorNoNFNoSFDown;
    }
    //
    return m_errDecompMap[poiName];
}

//__________________________________________________________________________________
//
void FittingTool::ExportErrorDecompositionInTextFile(const std::string& poiName, const std::string& fileName){
    // export breakdown to txt file
    if(fileName == "") return;

    std::ofstream out(fileName);
    if (!out.good() || !out.is_open()) {
        LOG(ERROR) << "Cannot open file: " << fileName << "\n";
        return;
    }

    out << "TOT_ERROR\t"      << m_errDecompMap[poiName]["TOT_ERROR"]      << "\t" << m_errDecompMapUp[poiName]["TOT_ERROR"]      << "\t" << m_errDecompMapDown[poiName]["TOT_ERROR"]      << std::endl;
    out << "STAT_ERROR\t"     << m_errDecompMap[poiName]["STAT_ERROR"]     << "\t" << m_errDecompMapUp[poiName]["STAT_ERROR"]     << "\t" << m_errDecompMapDown[poiName]["STAT_ERROR"]     << std::endl;
    out << "SYST_ERROR\t"     << m_errDecompMap[poiName]["SYST_ERROR"]     << "\t" << m_errDecompMapUp[poiName]["SYST_ERROR"]     << "\t" << m_errDecompMapDown[poiName]["SYST_ERROR"]     << std::endl;
    out << "MCSTAT_ERROR\t"   << m_errDecompMap[poiName]["MCSTAT_ERROR"]   << "\t" << m_errDecompMapUp[poiName]["MCSTAT_ERROR"]   << "\t" << m_errDecompMapDown[poiName]["MCSTAT_ERROR"]   << std::endl;
    out << "TOTSYST_ERROR\t"  << m_errDecompMap[poiName]["TOTSYST_ERROR"]  << "\t" << m_errDecompMapUp[poiName]["TOTSYST_ERROR"]  << "\t" << m_errDecompMapDown[poiName]["TOTSYST_ERROR"]  << std::endl;
    out << "STAT_ERR_NO_NF\t" << m_errDecompMap[poiName]["STAT_ERR_NO_NF"] << "\t" << m_errDecompMapUp[poiName]["STAT_ERR_NO_NF"] << "\t" << m_errDecompMapDown[poiName]["STAT_ERR_NO_NF"] << std::endl;
    out << "NF_ERROR\t"       << m_errDecompMap[poiName]["NF_ERROR"]       << "\t" << m_errDecompMapUp[poiName]["NF_ERROR"]       << "\t" << m_errDecompMapDown[poiName]["NF_ERROR"]       << std::endl;
    auto itr = m_errDecompMap[poiName].find("SHAPE_FACTOR_ERROR");
    if (itr != m_errDecompMap[poiName].end()) {
        out << "SHAPE_FACTOR_ERROR\t" << m_errDecompMap[poiName]["SHAPE_FACTOR_ERROR"] << "\t" << m_errDecompMapUp[poiName]["SHAPE_FACTOR_ERROR"] << "\t" << m_errDecompMapDown[poiName]["SHAPE_FACTOR_ERROR"] << std::endl;
        out << "STAT_ERR_NO_SF\t" << m_errDecompMap[poiName]["STAT_ERR_NO_SF"] << "\t" << m_errDecompMapUp[poiName]["STAT_ERR_NO_SF"] << "\t" << m_errDecompMapDown[poiName]["STAT_ERR_NO_SF"] << std::endl;
        out << "STAT_ERR_NO_NF_NO_SF\t" << m_errDecompMap[poiName]["STAT_ERR_NO_NF_NO_SF"] << "\t" << m_errDecompMapUp[poiName]["STAT_ERR_NO_NF_NO_SF"] << "\t" << m_errDecompMapDown[poiName]["STAT_ERR_NO_NF_NO_SF"] << std::endl;
    }
    for (auto np : m_errDecompMap[poiName]) {
        if (np.first == "TOT_ERROR" || np.first == "STAT_ERROR" || np.first == "SYST_ERROR") continue;
        if (np.first == "MCSTAT_ERROR" || np.first == "TOTSYST_ERROR") continue;
        if (np.first == "STAT_ERR_NO_NF" || np.first == "NF_ERROR") continue;
        if (np.first == "SHAPE_FACTOR_ERROR") continue;
        if (np.first == "STAT_ERR_NO_SF") continue;
        if (np.first == "STAT_ERR_NO_NF_NO_SF") continue;
        out << np.first << "\t" << np.second;
        out << "\t" << m_errDecompMapUp[poiName][np.first] << "\t" << m_errDecompMapDown[poiName][np.first];
        out << std::endl;
    }
    out << "\n";
    out.close();
}

//__________________________________________________________________________________
//
double FittingTool::GetPrefitUncertainty(const RooFitResult* fr, const std::string& param, int index) const {
    if (index < -1 || index > 1) {
        LOG(ERROR) << "Only -1, 0 and 1 indices are allowed\n";
        return -1;
    }

    for (auto arg : fr->floatParsInit()) {
        RooRealVar *par = static_cast<RooRealVar*>(arg);
        const std::string name = par->GetName();
        if (name != param) continue;

        if (index == -1) {
            return std::abs(par->getErrorLo());
        }
        if (index == 0) {
            return par->getError();
        }
        if (index == 1) {
            return par->getErrorHi();
        }
    }

    LOG(ERROR) << "Did not find parameter: " << param << "\n";
    // should not be reached
    return -1;
}

//__________________________________________________________________________________
//
double FittingTool::GetPrefitValue(const RooFitResult* fr, const std::string& param) const {
    for (auto arg : fr->floatParsInit()) {
        RooRealVar *par = static_cast<RooRealVar*>(arg);
        const std::string name = par->GetName();
        if (name != param) continue;

        return par->getVal();
    }

    LOG(ERROR) << "Did not find parameter: " << param << "\n";
    // should not be reached
    return -1;
}

//__________________________________________________________________________________
//
void FittingTool::ExportErrorDecompositionGroupInTextFile(const std::string& poiName,
                                                          const std::string& fileName,
                                                          const std::string& lumiLabel,
                                                          const std::string& cmeLabel,
                                                          const std::string& hepDataFolder,
                                                          const bool storeHEPData) {

    // export breakdown to txt file
    if(fileName == "") return;

    const std::string txtFileName = fileName + "_group_errDecomp_" + poiName + ".txt";

    LOG(INFO) << "Producing group uncertainty breakdown for POI: " << poiName << "\n";
    LOG(INFO) << "The txt file will be stored here: " << txtFileName << "\n";

    auto itrPOI     = m_errDecompMap.find(poiName);
    auto itrPOIUp   = m_errDecompMapUp.find(poiName);
    auto itrPOIDown = m_errDecompMapDown.find(poiName);
    if (itrPOI == m_errDecompMap.end() || itrPOIUp == m_errDecompMapUp.end() || itrPOIDown == m_errDecompMapDown.end()) {
        LOG(ERROR) << "Cannot find POI: " << poiName << " in the map\n";
        return;
    }

    std::ofstream out(txtFileName);
    if (!out.good() || !out.is_open()) {
        LOG(ERROR) << "Cannot open file: " << txtFileName << "\n";
        return;
    }

    const auto& map     = itrPOI->second;
    const auto& mapUp   = itrPOIUp->second;
    const auto& mapDown = itrPOIDown->second;

    double syst(0);
    double systUp(0);
    double systDown(0);

    std::vector<YamlConverter::ImpactContainer> container;

    std::map<std::string, double> tmp;
    for (const auto& icat : m_subCategories) {
        if (icat.empty()) continue;
        if (icat == "FullSyst") continue;
        if (icat == "NormFactors") continue;
        if (icat == "Gammas") continue;
        if (Common::StartsWith(icat, "gamma_")) continue;
        double unc(0.);
        double uncUp(0.);
        double uncDown(0.);
        for (const auto& iparam : m_subCategoryMap) {
            if (iparam.second != icat) continue;
            if (iparam.first == "DUMMY_STATONLY") continue;
            if (iparam.first == "DUMMY_GAMMAS") continue;
            if (Common::StartsWith(iparam.first, "gamma_")) continue;

            auto itrParam     = map.find(iparam.first);
            auto itrParamUp   = mapUp.find(iparam.first);
            auto itrParamDown = mapDown.find(iparam.first);
            if (itrParam == map.end() || itrParamUp == mapUp.end() || itrParamDown == mapDown.end()) {
                // the parameters do not exist in the fit results
                // which means they were pruned or fixed
                continue;
            }

            unc     += itrParam->second * itrParam->second;
            const double up   = itrParamUp->second;
            const double down = itrParamDown->second;

            if (up > 0) uncUp   += up*up;
            if (down > 0) uncUp += down*down;
            if (up < 0) uncDown   += up*up;
            if (down < 0) uncDown += down*down;
        }

        syst     += unc;
        systUp   += uncUp;
        systDown += uncDown;

        unc     = std::sqrt(unc);
        uncUp   = std::sqrt(uncUp);
        uncDown = std::sqrt(uncDown);

        out << icat << "    " << unc << "  ( +" << uncUp << ", -" << uncDown << " )\n";
        tmp.insert(std::make_pair(icat, unc));

        YamlConverter::ImpactContainer imp;
        imp.name = icat;
        imp.error = unc;
        imp.errorHi = uncUp;
        imp.errorLo = uncDown;

        container.emplace_back(imp);
    }

    // add gammas
    const double gamma     = map.at("MCSTAT_ERROR");
    const double gammaUp   = mapUp.at("MCSTAT_ERROR");
    const double gammaDown = mapDown.at("MCSTAT_ERROR");

    out << "Gammas" << "    " << gamma << "  ( +" << gammaUp << ", " << gammaDown << " )\n";

    syst     += gamma*gamma;
    systUp   += gammaUp*gammaUp;
    systDown += gammaDown*gammaDown;

    syst     = std::sqrt(syst);
    systUp   = std::sqrt(systUp);
    systDown = std::sqrt(systDown);

    {
        YamlConverter::ImpactContainer imp;
        imp.name = "MC stat.";
        imp.error = gamma;
        imp.errorHi = gammaUp;
        imp.errorLo = gammaDown;
        container.emplace_back(imp);
    }

    out << "FullSyst" << "    " << syst << "  ( +" << systUp << ", -" << systDown << " )\n";
    tmp.insert(std::make_pair("MC STAT", gamma));
    tmp.insert(std::make_pair("TOT_ERROR", map.at("TOT_ERROR")));
    tmp.insert(std::make_pair("TOTSYST_ERROR", map.at("TOTSYST_ERROR")));
    tmp.insert(std::make_pair("STAT_ERROR", map.at("STAT_ERROR")));
    m_categoryDecompositionMap.insert(std::make_pair(poiName, tmp));

    {
        YamlConverter::ImpactContainer imp;
        imp.name = "Total syst.";
        imp.error = map.at("TOTSYST_ERROR");
        imp.errorHi = mapUp.at("TOTSYST_ERROR");
        imp.errorLo = mapDown.at("TOTSYST_ERROR");
        container.emplace_back(imp);
    }

    {
        YamlConverter::ImpactContainer imp;
        imp.name = "Stat unc.";
        imp.error = map.at("STAT_ERROR");
        imp.errorHi = mapUp.at("STAT_ERROR");
        imp.errorLo = mapDown.at("STAT_ERROR");
        container.emplace_back(imp);
    }

    const std::string yamlFileName = fileName + "_group_errDecomp_" + poiName + ".yaml";
    YamlConverter converter{};
    converter.WriteImpact(container, yamlFileName);
    if (storeHEPData) {
        converter.SetLumi(Common::ReplaceString(lumiLabel, " fb^{-1}", ""));
        converter.SetCME(Common::ReplaceString(cmeLabel, " TeV", "000"));
        converter.WriteImpactHEPData(container, hepDataFolder, "_"+poiName);
    }

    out << "\n";
    out.close();
}

std::vector<YamlConverter::RankingContainer> FittingTool::GetRankingContainer(const std::string& poiName,
                                                                              const std::vector<std::string>& exclude,
                                                                              const bool flagSysts,
                                                                              const bool flagGammas) const {

    auto fr = this->GetFitResult();
    // check if fit-result is there
    if (!fr) {
        LOG(WARNING) << "No fit result, will return empty vector\n";
        return {};
    }

    auto itrPOIUp   = m_errDecompMapUp.find(poiName);
    auto itrPOIDown = m_errDecompMapDown.find(poiName);

    if (itrPOIUp == m_errDecompMapUp.end() || itrPOIDown == m_errDecompMapDown.end()) {
        LOG(ERROR) << "Cannot find POI: " << poiName << "\n";
        return {};
    }

    const auto& mapUp   = itrPOIUp->second;
    const auto& mapDown = itrPOIDown->second;

    std::vector<YamlConverter::RankingContainer> result;

    for (auto arg : fr->floatParsFinal()) {
        RooRealVar *par = static_cast<RooRealVar*>(arg);
        const std::string parName = par->GetName();
        if (poiName == parName) continue;

        auto itr = std::find_if(exclude.begin(), exclude.end(), [&parName](const auto& element){return parName.find(element) != std::string::npos;});
        if (itr != exclude.end()) continue;

        const bool isNP = Common::StartsWith(parName, "alpha_");
        const bool isGamma = Common::StartsWith(parName, "gamma_");

        if (isGamma && !flagGammas) continue;
        if (isNP    && !flagSysts) continue;

        auto itrParamUp   = mapUp.find(parName);
        auto itrParamDown = mapDown.find(parName);
        if (itrParamUp == mapUp.end() || itrParamDown == mapDown.end()) {
            LOG(ERROR) << "Cannot find parameter: " << parName << "\n";
            continue;
        }

        const double up   = itrParamUp->second;
        const double down = itrParamDown->second;

        YamlConverter::RankingContainer container{};
        if (isNP) {
            container.name = parName.substr(6, parName.size() - 6);;
        } else {
            container.name = parName;
        }
        if (isGamma) {
            container.nphat = (par->getVal() - 1.)/( par->getVal()>=1. ? std::abs(this->GetPrefitUncertainty(fr, parName, 1)) : std::abs(this->GetPrefitUncertainty(fr, parName, -1)) );
            container.nperrhi = par->getErrorHi() / this->GetPrefitUncertainty(fr, parName, 1);
            container.nperrlo = par->getErrorLo() / this->GetPrefitUncertainty(fr, parName, -1);
        } else {
            container.nphat = par->getVal();
            container.nperrhi = par->getErrorHi();
            container.nperrlo = par->getErrorLo();
        }
        container.poihi = up;
        container.poilo = down;
        container.poiprehi = 0;
        container.poiprelo = 0;

        result.emplace_back(std::move(container));
    }

    return result;
}
