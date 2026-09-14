#include "TRExFitter/LikelihoodScanManager.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/Logger.h"

#include "TRandom3.h"

#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"

#include <algorithm>

LikelihoodScanManager::LikelihoodScanManager() :
    fScanMinX(999999),
    fScanMinY(-999999),
    fStepsX(30),
    fStepX(-1),
    fScanMaxX(999999),
    fScanMaxY(-999999),
    fStepsY(30),
    fStepY(-1),
    fUseOffset(true),
    fCPU(1),
    fUseNllInLHscan(true),
    fMinValX(99999),
    fMinValY(99999),
    fMaxValX(-99999),
    fMaxValY(-99999),
    fExternalConstraints(nullptr),
    m_fitOptions({}),
    fToleranceScale(1.),
    fUseAutoDiff(false),
    fNLLOffset("initial")
{
    // shut-up RooFit!
    if(TRExFitter::DEBUGLEVEL<=1){
        if(TRExFitter::DEBUGLEVEL<=0) gErrorIgnoreLevel = kError;
        else gErrorIgnoreLevel = kWarning;
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Generation);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::LinkStateMgmt);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Optimization);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Tracing);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::Contents);
        RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
        RooMsgService::instance().setStreamStatus(1,false);
    }

}

//__________________________________________________________________________________
//
void LikelihoodScanManager::SetUseAutoDiff(const bool flag) {
    fUseAutoDiff = flag;
    if (flag) {
        LOG(INFO) << "Will use automatic differentiation for the likelihood scan\n";
    }
}

//________________________________________________________________________
//
void LikelihoodScanManager::SetNLLOffset(const std::string& val) {
    if (fUseOffset) {
        fNLLOffset = val;
        LOG(INFO) << "Likelihood offset set to: " << fNLLOffset << "\n";
    } else if (fNLLOffset != "none") {
        LOG(WARNING) << "SetNLLOffset call ignored since fUseOffset is false\n";
    }
}

//__________________________________________________________________________________
//
LikelihoodScanManager::ImpactScanResult1D LikelihoodScanManager::Run1DImpactScan(const RooWorkspace* ws,
                                                                                 const std::string& varNameX,
                                                                                 const std::string& varNameY,
                                                                                 RooDataSet* data) const {
    if (!ws || !data) {
        LOG(ERROR) << "Passed nullptr for ws\n";
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Passed nullptr for mc\n";
        exit(EXIT_FAILURE);
    }
    FitUtils::SetStartingAndConstParams(mc, m_fitOptions);
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    double min(-3);
    double max(3);
    LikelihoodScanManager::ImpactScanResult1D impactResult;

    auto FindVar = [&min,&max,&mc](RooRealVar*& var, const std::string& name, const bool isPOI) {
        if (!isPOI && !mc->GetNuisanceParameters()) return false;
        auto set = isPOI ? *mc->GetParametersOfInterest() : *mc->GetNuisanceParameters();
        for (auto var_tmp : set) {
            var = static_cast<RooRealVar*>(var_tmp);
            const std::string vname = var->GetName();
            if (vname == name) {
                LOG(INFO) << "GetImpactScan on parameter = " << vname << "\n";
                min = var->getMin();
                max = var->getMax();
                return true;
            }
        }
        return false;
    };

    RooRealVar* varY(nullptr);
    bool found = FindVar(varY, varNameY, true);

    if(!found) {
        LOG(WARNING) << "Cannot find POI " << varNameY << " to evaluate impact on.\n";
        return impactResult;
    }
    RooRealVar* varX(nullptr);
    found = FindVar(varX, varNameX, false);
    if (!found) {
        LOG(WARNING) << "Cannot find NP " << varNameX << " to evaluate impact of.\n";
        return impactResult;
    }

    if (!varX || !varY) {
        LOG(ERROR) << "Variables are nullpts!\n";
        return impactResult;
    }

    if (fScanMinX < 99999) { // is actually set
        min = fScanMinX;
    }

    if (fScanMaxX > -99999) { // is actually set
        max = fScanMaxX;
    }

    const auto offset = fUseOffset ? fNLLOffset : "none";
    const RooArgSet* glbObs = mc->GetGlobalObservables();

    const RooArgSet tmpDummy{};
    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (fUseAutoDiff) {
        nll = std::unique_ptr<RooAbsReal> (simPdf->createNLL(*data,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Offset(offset),
                                           RooFit::NumCPU(fCPU, RooFit::Hybrid),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::ExternalConstraints(fExternalConstraints ? *fExternalConstraints : tmpDummy),
                                           RooFit::EvalBackend("codegen")));
    } else {
        nll = std::unique_ptr<RooAbsReal> (simPdf->createNLL(*data,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Offset(offset),
                                           RooFit::NumCPU(fCPU, RooFit::Hybrid),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::ExternalConstraints(fExternalConstraints ? *fExternalConstraints : tmpDummy)));
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this

    varX->setConstant(kTRUE); // make POI constant in the fit

    const std::size_t nFree = FitUtils::NumberOfFreeParameters(mc);
    std::vector<double> steps;
    std::vector<double> strengths;
    steps.resize(fStepsX);
    strengths.resize(fStepsX);

    for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {

        RooMinimizer m(*nll);
        m.setPrintLevel(-1);
        m.setStrategy(1);
        m.optimizeConst(2);
        m.setMinimizerType(minimType.Data());
        m.setEps(tol*fToleranceScale);

        LOG(INFO) << "Running ImpactScan for point " << (ipoint+1) << " out of " << fStepsX << " points\n";
        steps[ipoint] = min+ipoint*(max-min)/(fStepsX - 1);
        *varX = steps[ipoint]; // set POI
        std::unique_ptr<RooFitResult> r(nullptr);
        if (nFree != 0) {
            m.migrad(); // minimize again with new posSigXsecOverSM value
            r = std::unique_ptr<RooFitResult>(m.save());
        }
        strengths[ipoint] = varY->getValV();
        LOG(DEBUG) << "Point: " << steps[ipoint] << ", impact: " << strengths[ipoint] << "\n";
    }
    varX->setConstant(kFALSE);

    TRandom3 rand(1234567);
    if (std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNameX) != fBlindedParameters.end()) {
        const double rndNumber = rand.Uniform(5);
        for (auto& ix : steps) {
            ix += rndNumber;
        }
    }

    impactResult.steps = steps;
    impactResult.strengths = strengths;
    return impactResult;
}



//__________________________________________________________________________________
//
LikelihoodScanManager::scanResult1D LikelihoodScanManager::Run1DScan(const RooWorkspace* ws,
                                                                     const std::string& varName,
                                                                     RooDataSet* data) const {

    if (!ws || !data) {
        LOG(ERROR) << "Passed nullptr for ws\n";
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Passed nullptr for mc\n";
        exit(EXIT_FAILURE);
    }

    FitUtils::SetStartingAndConstParams(mc, m_fitOptions);

    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    double min(-3);
    double max(3);

    bool found(false);
    RooRealVar* var(nullptr);
    {
        for (auto var_tmp : *mc->GetParametersOfInterest()) {
            var = static_cast<RooRealVar*>(var_tmp);
            const std::string vname = var->GetName();
            if (vname == varName) {
                LOG(INFO) << "GetLikelihoodScan for NP = " << vname << "\n";
                found=true;
                min = var->getMin();
                max = var->getMax();
                break;
            }
        }
    }

    if (!found) {
        if (mc->GetNuisanceParameters()) {
            for (auto var_tmp : *mc->GetNuisanceParameters()) {
                var = static_cast<RooRealVar*>(var_tmp);
                const std::string vname = var->GetName();
                if (vname == varName || vname == "alpha_"+varName) {
                    LOG(INFO) << "GetLikelihoodScan for NP = " << vname << "\n";
                    found=true;

                    if ((vname.find("gamma_") != std::string::npos)) {
                        min = 0.5;
                        max = 1.5;
                    } else {
                        min = var->getMin();
                        max = var->getMax();
                    }

                    break;
                }
            }
        }
    }

    if (fScanMinX < 99999) { // is actually set
        min = fScanMinX;
    }

    if (fScanMaxX > -99999) { // is actually set
        max = fScanMaxX;
    }

    LikelihoodScanManager::scanResult1D result;
    if (!found) {
        LOG(WARNING) << "Cannot find NP: " << varName << "\n";
        return result;
    }

    result.first.resize(fStepsX);
    result.second.resize(fStepsX);

    const auto offset = fUseOffset ? kTRUE : kFALSE;
    const RooArgSet* glbObs = mc->GetGlobalObservables();

    const RooArgSet tmpDummy{};
    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (fUseAutoDiff) {
        nll = std::unique_ptr<RooAbsReal> (simPdf->createNLL(*data,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::ExternalConstraints(fExternalConstraints ? *fExternalConstraints : tmpDummy),
                                           RooFit::EvalBackend("codegen")));
    } else {
        nll = std::unique_ptr<RooAbsReal> (simPdf->createNLL(*data,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Offset(offset),
                                           RooFit::NumCPU(fCPU, RooFit::Hybrid),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::ExternalConstraints(fExternalConstraints ? *fExternalConstraints : tmpDummy)));
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this

    var->setConstant(kTRUE); // make POI constant in the fit

    const std::size_t nFree = FitUtils::NumberOfFreeParameters(mc);

    for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {

        if (fStepX != -1 && fStepX != ipoint) {
            continue;
        }

        RooMinimizer m(*nll);
        m.setPrintLevel(-1);
        m.setStrategy(1);
        m.optimizeConst(2);
        m.setMinimizerType(minimType.Data());
        m.setEps(tol*fToleranceScale);

        LOG(INFO) << "Running LHscan for point " << (ipoint+1) << " out of " << fStepsX << " points\n";
        result.first[ipoint] = min+ipoint*(max-min)/(fStepsX - 1);
        *var = result.first[ipoint]; // set POI
        std::unique_ptr<RooFitResult> r(nullptr);
        if (nFree != 0) {
            m.migrad(); // minimize again with new posSigXsecOverSM value
            r = std::unique_ptr<RooFitResult>(m.save());
        }
        const double nllval = nll->getVal();
        if (fUseNllInLHscan || !r) {
            result.second[ipoint] = nllval;
        } else {
            result.second[ipoint] = r->minNll();
        }
        LOG(DEBUG) << "Point: " << result.first[ipoint] << ", nll: " << result.second[ipoint] << "\n";
    }
    var->setConstant(kFALSE);

    const double mnll = *std::min_element(result.second.begin(), result.second.end());

    if (fStepX == -1) {
        // Reduce minimum in default case only
        for (auto & iY : result.second) {
            iY = iY - mnll;
        }
    }

    TRandom3 rand(1234567);
    if (std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varName) != fBlindedParameters.end()) {
        const double rndNumber = rand.Uniform(5);
        for (auto& ix : result.first) {
            ix += rndNumber;
        }
    }

    // If parallel processing points, keep only the one point
    if (fStepX != -1) {
        LikelihoodScanManager::scanResult1D tempresult;
        tempresult.first.resize(1);
        tempresult.second.resize(1);
        tempresult.first[0] = result.first[fStepX];
        tempresult.second[0] = result.second[fStepX];
        return tempresult;
    }


    return result;
}

//__________________________________________________________________________________
//
LikelihoodScanManager::Result2D LikelihoodScanManager::Run2DScan(const RooWorkspace* ws,
                                                                 const std::pair<std::string, std::string>& varNames,
                                                                 RooDataSet* data) {

    if (!ws || !data) {
        LOG(ERROR) << "Passed nullptr for ws\n";
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Passed nullptr for mc\n";
        exit(EXIT_FAILURE);
    }
    FitUtils::SetStartingAndConstParams(mc, m_fitOptions);
    RooSimultaneous* simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    int count = 0;
    RooRealVar* varX = nullptr;
    RooRealVar* varY = nullptr;
    if (mc->GetNuisanceParameters()) {
        for (auto var_tmp : *mc->GetNuisanceParameters()) {
            RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
            const std::string vname = var->GetName();
            if (vname == varNames.first || vname == "alpha_"+varNames.first){
                varX = var;
                ++count;
            }
            if (vname == varNames.second || vname == "alpha_"+varNames.second){
                varY = var;
                ++count;
            }
            if (count == 2) break;
        }
    }

    if (count != 2) {
        for (auto var_tmp : *mc->GetParametersOfInterest()) {
            RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
            const std::string vname = var->GetName();
            if (vname == varNames.first || vname == "alpha_"+varNames.first){
                varX = var;
                ++count;
            }
            if (vname == varNames.second || vname == "alpha_"+varNames.second){
                varY = var;
                ++count;
            }
            if (count == 2) break;
        }
    }

    LikelihoodScanManager::Result2D result;

    if (count != 2) {
        LOG(WARNING) << "Did not find the two parameters you want to use in the 2D likelihood scan\n";
        return result;
    }

    double minValX = varX->getMin();
    double maxValX = varX->getMax();
    double minValY = varY->getMin();
    double maxValY = varY->getMax();

    if (fScanMinX < 99999) { // is actually set
        minValX = fScanMinX;
    }
    if (fScanMinY < 99999) { // is actually set
        minValY = fScanMinY;
    }
    if (fScanMaxX > -99999) { // is actually set
        maxValX = fScanMaxX;
    }
    if (fScanMaxY > -99999) { // is actually set
        maxValY = fScanMaxY;
    }

    const auto offset = fUseOffset ? kTRUE : kFALSE;
    const RooArgSet* glbObs = mc->GetGlobalObservables();

    const RooArgSet tmpDummy{};
    std::unique_ptr<RooAbsReal> nll(nullptr);
    if (fUseAutoDiff) {
        nll = std::unique_ptr<RooAbsReal> (simPdf->createNLL(*data,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::ExternalConstraints(fExternalConstraints ? *fExternalConstraints : tmpDummy),
                                           RooFit::EvalBackend("codegen")));
    } else {
        nll = std::unique_ptr<RooAbsReal> (simPdf->createNLL(*data,
                                           RooFit::GlobalObservables(*glbObs),
                                           RooFit::Offset(offset),
                                           RooFit::NumCPU(fCPU, RooFit::Hybrid),
                                           RooFit::Optimize(kTRUE),
                                           RooFit::ExternalConstraints(fExternalConstraints ? *fExternalConstraints : tmpDummy)));
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance();

    varX->setConstant(kTRUE); // make POI constant in the fit
    varY->setConstant(kTRUE); // make POI constant in the fit

    const std::size_t nFree = FitUtils::NumberOfFreeParameters(mc);

    result.x.resize(fStepsX);
    result.y.resize(fStepsY);
    result.z.resize(fStepsX);
    for (auto& iz : result.z) {
        iz.resize(fStepsY);
    }

    //values for parameter1, parameter2 and the NLL value
    for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {
        if (fStepX != -1 && fStepX != ipoint) {
            continue;
        }

        //Set both POIs to constant
        LOG(INFO) << "Running LHscan for point " << (ipoint+1) << " out of " << fStepsX << " points\n";
        result.x[ipoint] = minValX + ipoint * (maxValX - minValX) / (fStepsX - 1);
        *varX = result.x[ipoint]; // set POI
        for (int jpoint = 0; jpoint < fStepsY; ++jpoint) {
            if (fStepY != -1 && fStepY != jpoint) {
                continue;
            }
            RooMinimizer m(*nll);
            m.optimizeConst(2);
            m.setErrorLevel(-1);
            m.setPrintLevel(-1);
            m.setStrategy(1); // set precision to high
            m.setMinimizerType(minimType.Data());
            m.setEps(tol*fToleranceScale);
            LOG(INFO) << "Running LHscan for subpoint " << (jpoint+1) << " out of " << fStepsY <<  " points\n";
            result.y[jpoint] = minValY + jpoint * (maxValY - minValY) / (fStepsY - 1);
            *varY = result.y[jpoint]; // set POI
            std::unique_ptr<RooFitResult> r(nullptr);
            if (nFree != 0) {
                m.migrad(); // minimize again with new posSigXsecOverSM value
                r = std::unique_ptr<RooFitResult>(m.save());
            }
            const double nllval = nll->getVal();
            double z_tmp(0);
            if (fUseNllInLHscan || !r) {
                z_tmp = nllval;
            } else {
                z_tmp = r->minNll();
            }
            result.z[ipoint][jpoint] = z_tmp;
            LOG(DEBUG) << "Point x: " << result.x[ipoint] << ", y: " << result.y[jpoint] << ", nll: " << result.z[ipoint][jpoint] << "\n";
        }
    }
    varX->setConstant(kFALSE);
    varY->setConstant(kFALSE);

    TRandom3 rand(1234567);
    const bool blindX = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.first) != fBlindedParameters.end();
    const bool blindY = std::find(fBlindedParameters.begin(), fBlindedParameters.end(), varNames.second) != fBlindedParameters.end();
    const double rndNumber = rand.Uniform(5);
    if (blindX) {
        minValX += rndNumber;
        maxValX += rndNumber;
        for (auto& ix : result.x) {
            ix += rndNumber;
        }
    }

    if (blindY) {
        minValY += rndNumber;
        maxValY += rndNumber;
        for (auto& iy : result.y) {
            iy += rndNumber;
        }
    }

    fMinValX = minValX;
    fMinValY = minValY;
    fMaxValX = maxValX;
    fMaxValY = maxValY;


    // If parallel processing points, return a struct that only contains the filled points
    if (fStepX != -1 || fStepY != -1) {
        for (int ipoint = 0; ipoint < fStepsX; ++ipoint) {
	    for (int jpoint = 0; jpoint < fStepsY; ++jpoint) {
	      if (ipoint == fStepX && jpoint== fStepY) continue;
	      if (fStepY == -1 && ipoint == fStepX) continue;
	      if (fStepX == -1 && jpoint == fStepY) continue;
	      result.z[ipoint][jpoint] = 0.0;
	    }
	}
    }

    return result;
}

void LikelihoodScanManager::SetFittingOptions(const FitUtils::FittingOptions& options) {
    m_fitOptions.constNPs = options.constNPs;
    m_fitOptions.constNPvalues = options.constNPvalues;
    m_fitOptions.initialNPs = options.initialNPs;
    m_fitOptions.initialNPvalues = options.initialNPvalues;
    m_fitOptions.noGammas = options.noGammas;
    m_fitOptions.noSystematics = options.noSystematics;
    m_fitOptions.randomize = options.randomize;
    m_fitOptions.randomValue = options.randomValue;
    m_fitOptions.constPOI = options.constPOI;
    m_fitOptions.poiVals = options.poiVals;
}
