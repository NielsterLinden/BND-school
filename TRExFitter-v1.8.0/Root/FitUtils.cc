#include "TRExFitter/FitUtils.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/ShapeFactor.h"

#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooMultiVarGaussian.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooTFnBinding.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "xRooFit/xRooNode.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TVectorD.h"

#include <map>
#include <sstream>

//__________________________________________________________________________________
//
void FitUtils::ApplyExternalConstraints(RooWorkspace* ws,
                                        FittingTool* fitTool,
                                        RooSimultaneous* simPdf,
                                        const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                                        int type
                                       ) {

    const RooStats::ModelConfig* mc = dynamic_cast<const RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read the ModelConfig!\n";
        exit(EXIT_FAILURE);
    }
    const std::vector<std::string> params = FitUtils::GetAllParameters(mc);

    // Tikhonov regularization (for unfolding)
    // only add the unique NFs
    std::vector<std::string> names;
    std::vector<std::shared_ptr<NormFactor> > uniqueNF;
    for(const auto& nf : normFactors) {
        if (std::find(names.begin(), names.end(), nf->fName) != names.end()) continue;
        names.emplace_back(nf->fName);
        uniqueNF.emplace_back(nf);
    }

    std::vector<std::string> regularizationNames;

    // Simple form (individual constraint terms on signal strengths)
    if(type==0){
        std::vector<double> tauVec;
        RooArgList l;
        std::vector<double> nomVec;
        for(const auto& nf : uniqueNF) {
            if(nf->fTau == 0) continue;

            if(ws->var(nf->fName.c_str())) {
                l.add(*ws->var(nf->fName.c_str()));
            } else if (ws->function(nf->fName.c_str())) {
                l.add(*ws->function(nf->fName.c_str()));
            } else {
                LOG(WARNING) << "Cannot apply tau to norm factor " << nf->fName << "\n";
                continue;
            }

            nomVec.push_back( nf->GetNominal() );
            tauVec.push_back( nf->fTau );
        }

        if(tauVec.empty()) return;

        TVectorD nominal(nomVec.size());
        TMatrixDSym cov(tauVec.size());
        for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++){
            nominal(i_tau) = nomVec[i_tau];
            cov(i_tau,i_tau) = (1./tauVec[i_tau]) * (1./tauVec[i_tau]);
        }
        RooMultiVarGaussian r("regularization","regularization",l,nominal,cov);
        ws->import(r);
        regularizationNames.emplace_back("regularization");
    }
    // Discretized second derivative
    else {
        RooArgList l_prime;
        std::map<std::string, std::vector<std::shared_ptr<NormFactor> > > unfoldingNF;
        for(const auto& nf : uniqueNF) {
            if ((nf->fName).find("_Bin_") == std::string::npos) continue;
            // in case of multiple unfoldings, we need to match the correct ones
            const std::string unfoldingName = nf->fName.substr(0, nf->fName.find("_Bin_"));
            auto itr = unfoldingNF.find(unfoldingName);
            if (itr == unfoldingNF.end()) {
                std::vector<std::shared_ptr<NormFactor> > tmp = {nf};
                unfoldingNF.insert({unfoldingName, tmp});
            } else {
                itr->second.emplace_back(nf);
            }
        }


        for (auto& iunfolding : unfoldingNF) {
            std::vector<std::size_t> skip = {0, iunfolding.second.size() - 1};
            if (type > 1) {
                skip.emplace_back(1);
                skip.emplace_back(iunfolding.second.size() - 2);
            }
            if (type > 2) {
                skip.emplace_back(2);
                skip.emplace_back(iunfolding.second.size() - 3);
            }
            if (type > 3) {
                skip.emplace_back(3);
                skip.emplace_back(iunfolding.second.size() - 4);
            }
            std::vector<double> tauVec;
            for(std::size_t i_nf =0; i_nf < iunfolding.second.size(); ++i_nf) {
                if (std::find(skip.begin(), skip.end(), i_nf) != skip.end()) continue;
                if(iunfolding.second[i_nf]->fTau == 0) continue;
                tauVec.push_back(iunfolding.second[i_nf]->fTau);

                if (type == 1) {
                    RooFormulaVar f(Form("f_bin%ld",i_nf+1), "@0 - 2*@1 + @2", RooArgList(
                        *ws->function(iunfolding.second[i_nf-1]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+1]->fName.c_str())
                    ));

                    f.SetName(Form("mu_prime_%ld",i_nf+1));
                    ws->import(f);
                    l_prime.add(*ws->function(Form("mu_prime_%ld",i_nf+1)));
                } else if (type == 2) {
                    RooFormulaVar f(Form("f_bin%ld",i_nf+1), "@0 - 4*@1 + 6*@2 - 4*@3 + @4", RooArgList(
                        *ws->function(iunfolding.second[i_nf-2]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf-1]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+1]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+2]->fName.c_str())
                    ));

                    f.SetName(Form("mu_prime_%ld",i_nf+1));
                    ws->import(f);
                    l_prime.add(*ws->function(Form("mu_prime_%ld",i_nf+1)));
                } else if (type == 3) {
                    RooFormulaVar f(Form("f_bin%ld",i_nf+1), "@0 - 6*@1 + 15*@2 - 20*@3 + 15*@4 - 6*@5 + @6", RooArgList(
                        *ws->function(iunfolding.second[i_nf-3]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf-2]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf-1]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+1]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+2]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+3]->fName.c_str())
                    ));

                    f.SetName(Form("mu_prime_%ld",i_nf+1));
                    ws->import(f);
                    l_prime.add(*ws->function(Form("mu_prime_%ld",i_nf+1)));
                } else if (type == 4) {
                    RooFormulaVar f(Form("f_bin%ld",i_nf+1), "@0 - 8*@1 + 28*@2 - 56*@3 + 70*@4 - 56*@5 + 28*@6 - 8*@7 + @8", RooArgList(
                        *ws->function(iunfolding.second[i_nf-4]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf-3]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf-2]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf-1]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+1]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+2]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+3]->fName.c_str()),
                        *ws->function(iunfolding.second[i_nf+4]->fName.c_str())
                    ));

                    f.SetName(Form("mu_prime_%ld",i_nf+1));
                    ws->import(f);
                    l_prime.add(*ws->function(Form("mu_prime_%ld",i_nf+1)));
                } else {
                    LOG(ERROR) << "Unsupported regularisation option " << type << "\n";
                    exit(EXIT_FAILURE);
                }
            }

            TVectorD nominal(tauVec.size());
            TMatrixDSym cov(tauVec.size());
            for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++){
                nominal(i_tau) = 0;
                cov(i_tau,i_tau) = (1./tauVec[i_tau]) * (1./tauVec[i_tau]);
            }
            RooMultiVarGaussian r(("regularization_"+iunfolding.first).c_str(),("regularization_"+iunfolding.first).c_str(),l_prime,nominal,cov);
            ws->import(r);
            regularizationNames.emplace_back("regularization_"+iunfolding.first);
        }
    }

    // combine regularizationNames into a single, comma separated string
    std::string reg("");
    for (std::size_t i = 0; i < regularizationNames.size(); ++i) {
        reg += regularizationNames.at(i);
        if (i != regularizationNames.size() -1) {
            reg += ",";
        }
    }

    ws->defineSet("myConstraints",reg.c_str());
    simPdf->setStringAttribute("externalConstraints","myConstraints");

    if(simPdf->getStringAttribute("externalConstraints")){
        LOG(INFO) << Form("Building NLL with external constraints %s",simPdf->getStringAttribute("externalConstraints")) << "\n";
        const RooArgSet* externalConstraints = ws->set(simPdf->getStringAttribute("externalConstraints"));
        if (fitTool) fitTool->SetExternalConstraints( externalConstraints );
    }
}

//__________________________________________________________________________________
//
void FitUtils::SetBinnedLikelihoodOptimisation(RooWorkspace* ws) {
    for (auto arg : ws->components()) {
        if (arg->IsA() == RooRealSumPdf::Class()) {
            arg->setAttribute("BinnedLikelihood");
            const std::string temp_string = arg->GetName();
            LOG(DEBUG) << "Activating binned likelihood attribute for " << temp_string << "\n";
        }
    }
}

//__________________________________________________________________________________
//
void FitUtils::InjectGlobalObservables(RooWorkspace* ws, const FitResults* fr) {
    std::map<std::string, double> npValues;
    for (const auto& inp : fr->GetNuisanceParameters()) {
        const std::string name = inp.first;
        const double value = fr->GetNuisParValue(name);
        npValues.insert({name, value});
    }
    FitUtils::InjectGlobalObservables(ws, npValues);
}

//__________________________________________________________________________________
//
void FitUtils::InjectGlobalObservables(RooWorkspace* ws,
                                       const std::map< std::string, double >& npValues) {
    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read ModelCOnfig\n";
        exit(EXIT_FAILURE);
    }
    RooArgSet mc_globs = *mc->GetGlobalObservables();

    LOG(INFO) << "Injecting the following NP values to global observables\n";
    for(const auto& np_value : npValues) {
        const std::string this_name = np_value.first;
        double this_value = np_value.second;

        std::ostringstream tmp;

        // find the corresponding glob
        const std::string glob_name = "nom_" + this_name;
        const std::string glob_name_alpha = "nom_alpha_" + this_name;
        const std::string glob_name_gamma = "nom_gamma_" + this_name;
        RooRealVar* this_glob = nullptr;
        for ( auto glob_tmp : mc_globs) {
            RooRealVar* glob = static_cast<RooRealVar*>(glob_tmp);
            if(glob->GetName() == glob_name || glob->GetName() == glob_name_alpha || glob->GetName() == glob_name_gamma) {
                this_glob = glob;
                break;
            }
        }

        if(!this_glob) {
            LOG(WARNING) << "Could not find global observable " << glob_name << "\n";
            continue;
        }

        // set gamma values to gamma*nom_gamma
        if(glob_name.find("nom_gamma_") == 0) {
            this_value = this_value * this_glob->getVal();
        }

        tmp << this_glob->GetName() << ": " << this_value;
        LOG(INFO) << tmp.str() << "\n";
        this_glob->setVal(this_value);
    }
}

//__________________________________________________________________________________
//
void FitUtils::SetPOIinFile(const std::string& path, const std::string& wsName, const std::string& poi) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str()));
    if (!f) {
        LOG(ERROR) << "Cannot open input file at: " << path << "\n";
        return;
    }

   RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f->Get(wsName.c_str()));
    if (!ws) {
        LOG(ERROR) << "Cannot read workspace\n";
        return;
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read ModelConfig\n";
        return;
    }

    mc->SetParametersOfInterest(poi.c_str());

    f->Close();
}

//__________________________________________________________________________________
//
std::vector<std::string> FitUtils::GetAllParameters(const RooStats::ModelConfig* mc) {

    if (!mc) {
        LOG(ERROR) << "ModelConfig is nullptr\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> params;
    for (const auto& i : *mc->GetParametersOfInterest()) {
        params.emplace_back(i->GetName());
    }

    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            params.emplace_back(i->GetName());
        }
    }

    return params;
}


//__________________________________________________________________________________
//
std::size_t FitUtils::NumberOfFreeParameters(const RooStats::ModelConfig* mc) {
    if (!mc) {
        LOG(ERROR) << "ModelConfig is nullptr\n";
        exit(EXIT_FAILURE);
    }

    std::size_t result(0);

    for (const auto& i : *mc->GetParametersOfInterest()) {
        if (!i->isConstant()) ++result;
    }

    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            if (!i->isConstant()) ++result;
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void FitUtils::FixAllParameters(RooStats::ModelConfig* mc, const std::vector<std::string>& exceptions) {
    if (!mc) {
        LOG(ERROR) << "ModelConfig is nullptr\n";
        exit(EXIT_FAILURE);
    }

    for (const auto& i : *mc->GetParametersOfInterest()) {
        auto* tmp = static_cast<RooRealVar*>(i);
        tmp->setConstant(1);
    }

    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            auto* tmp = static_cast<RooRealVar*>(i);
            const std::string name = tmp->GetName();

            auto it = std::find_if(exceptions.begin(), exceptions.end(),
                [&name](const std::string& element){return name.find(element) != std::string::npos;});

            if (it != exceptions.end()) continue;
            tmp->setConstant(1);
        }
    }
}

//__________________________________________________________________________________
//
void FitUtils::FloatAllParameters(RooStats::ModelConfig* mc) {
    if (!mc) {
        LOG(ERROR) << "ModelCOnfig is nullptr\n";
        exit(EXIT_FAILURE);
    }

    for (const auto& i : *mc->GetParametersOfInterest()) {
        auto* tmp = static_cast<RooRealVar*>(i);
        tmp->setConstant(0);
    }

    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            auto* tmp = static_cast<RooRealVar*>(i);
            const std::string name = tmp->GetName();
            tmp->setConstant(0);
            if (name.find("gamma_") != std::string::npos) {
              tmp->setVal(1.);
            } else {
              tmp->setVal(0.);
            }
        }
    }
}

//__________________________________________________________________________________
//
std::vector<double> FitUtils::CalculateExpressionRoofit(RooWorkspace* ws,
                                                        TFile* fitResultFile,
                                                        const std::string& name,
                                                        const bool statOnly,
                                                        const std::vector<std::string>& freeParams) {

    if (!ws || !fitResultFile) {
        LOG(ERROR) << "Passed nullptr\n";
        exit(EXIT_FAILURE);
    }

    std::vector<double> result;

    std::unique_ptr<RooFitResult> fr(nullptr);
    for(auto *key : *fitResultFile->GetListOfKeys()){
        const std::string& keyName = key->GetName();
        const std::string nllName = "nll_simPdf";
        const std::string evaluatorName = "RooEvaluatorWrapper";
        auto itrNll = keyName.find(nllName);
        if (itrNll != std::string::npos) {
            fr.reset(fitResultFile->Get<RooFitResult>(keyName.c_str()));
            break;
        } else {
            if (keyName.find(evaluatorName) != std::string::npos) {
                fr.reset(fitResultFile->Get<RooFitResult>(keyName.c_str()));
                break;
            }
        }
    }

    if (!fr) {
        LOG(ERROR) << "FitResult is nullptr\n";
        exit(EXIT_FAILURE);
    }

    RooFormulaVar* muN(nullptr);
    for(auto *arg : ws->allFunctions()) {
        const std::string argName = arg->GetName();
        if(argName == name) {
            muN = &dynamic_cast<RooFormulaVar&>(*arg);
            break;
        }
    }

    if (!muN) {
        LOG(ERROR) << "Did not find the correct name: " << name << "\n";
        return result;
    }

    for(auto *par : fr->floatParsFinal()){
        ws->var( par->GetName() )->setVal( (static_cast<RooRealVar*>(par))->getVal() );
        ws->var( par->GetName() )->setError( (static_cast<RooRealVar*>(par))->getError() );
    }

    const double mean  = muN->getVal();

    double error = muN->getPropagatedError(*fr);

    // need to subtract the propagated impacts of the NPs
    if (statOnly) {
        double systImpactSq = 0;
        for (auto *par : fr->floatParsFinal()) {
            const std::string parName = par->GetName();
            if (parName == name) continue;
            auto itr = std::find(freeParams.begin(), freeParams.end(), parName);
            if (itr != freeParams.end()) continue;

            const double paramError = FitUtils::GetPropagatedCovariance(*muN, *static_cast<RooRealVar*>(par), *fr, {});
            systImpactSq += paramError*paramError;
        }

        error = std::sqrt(error*error - systImpactSq);
    }

    const double up    = mean + error;
    const double down  = mean - error;

    result.emplace_back(mean);
    result.emplace_back(up);
    result.emplace_back(down);
    result.emplace_back(error);

    return result;
}

//__________________________________________________________________________________
//
RooArgSet FitUtils::GetPOIsWithout(const RooStats::ModelConfig* mc, const std::vector<std::string>& names) {
    RooArgSet result;
    for (const auto i : *mc->GetParametersOfInterest()) {
        auto* tmp = static_cast<RooRealVar*>(i);
        const std::string name = tmp->GetName();
        auto it = std::find(names.begin(), names.end(), name);
        if (it != names.end()) continue;

        result.add(*i);
    }
    if (mc->GetNuisanceParameters()) {
        for (const auto i : *mc->GetNuisanceParameters()) {
            auto* tmp = static_cast<RooRealVar*>(i);
            const std::string name = tmp->GetName();
            auto it = std::find(names.begin(), names.end(), name);
            if (it != names.end()) continue;

            result.add(*i);
        }
     }

    return result;
}

//__________________________________________________________________________________
//
const RooArgSet* FitUtils::GetExternalConstraints(RooWorkspace* ws,
                                                  const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                                                  const int regularizationType) {

    RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());
    const RooArgSet* externalConstraints(nullptr);
    FitUtils::ApplyExternalConstraints(ws, nullptr, simPdf, normFactors, regularizationType);
    if (simPdf->getStringAttribute("externalConstraints")) {
        externalConstraints = ws->set(simPdf->getStringAttribute("externalConstraints"));
    }
    return externalConstraints;
}

//__________________________________________________________________________________
//
void FitUtils::FixParameters(RooStats::ModelConfig* mc,
                             const std::vector<std::pair<std::string, double> >& params) {

    auto* nps = mc->GetNuisanceParameters();
    if (nps) {
        for (auto ipar : *nps) {
            const std::string name = ipar->GetName();
            auto it = std::find_if(params.begin(), params.end(), [&name](const std::pair<std::string, double>& element){return name.find(element.first) != std::string::npos;});
            if (it != params.end()) {
                auto par = static_cast<RooRealVar*>(ipar);
                LOG(INFO) << "Fixing param: " << name << " to " << it->second << "\n";
                par->setVal(it->second);
                par->setConstant(true);
            }
        }
    }

    for (auto inf : *mc->GetParametersOfInterest()) {
        const std::string name = inf->GetName();
        auto it = std::find_if(params.begin(), params.end(), [&name](const auto& element){return element.first == name;});
        if (it != params.end()) {
            auto par = static_cast<RooRealVar*>(inf);
            LOG(INFO) << "Fixing param: " << name << " to " << it->second << "\n";
            par->setVal(it->second);
            par->setConstant(true);
        }
    }
}

//__________________________________________________________________________________
//
void FitUtils::ChangeInterpolationCode(RooWorkspace* ws, const int intCode) {
    if (!ws) {
        LOG(ERROR) << "Passed empty workspace\n";
        exit(EXIT_FAILURE);
    }

    RooArgSet funcs = ws->allFunctions();
    for (auto iobj : funcs) {
        RooStats::HistFactory::FlexibleInterpVar* flex = dynamic_cast<RooStats::HistFactory::FlexibleInterpVar*>(iobj);
        if (!flex) continue;
        flex->setAllInterpCodes(intCode);
    }
}

//__________________________________________________________________________________
//
void FitUtils::ReparametrizeShapeFactors(RooWorkspace* ws,
                                         const std::vector<std::shared_ptr<ShapeFactor> >& shapeFactors,
                                         const std::vector<Region*>& regions,
                                         const std::vector<std::string>& regionsUsedInFit) {
    LOG(INFO) << "Reparametrizing the shape factors in the WS\n";

    if (!ws) {
        LOG(ERROR) << "Cannot read the WS from the file\n";
        exit(EXIT_FAILURE);
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read the model config\n";
        exit(EXIT_FAILURE);
    }

    xRooNode(*ws).sterilize();

    for (const auto& sf : shapeFactors) {
        if (sf->fExpression.empty()) continue;
        const int nBins = sf->fExpression.size();
        const std::string regName = sf->fRegions.at(0);
        // skip the regions not used in the fit
        if (std::find(regionsUsedInFit.begin(), regionsUsedInFit.end(), regName) == regionsUsedInFit.end()) continue;
        auto itr = std::find_if(regions.begin(), regions.end(), [&regName](const auto& reg){return reg->fName == regName;});
        if (itr == regions.end()) {
            LOG(ERROR) << "Unknown region: " << regName << "\n";
            exit(EXIT_FAILURE);
        }

        if (nBins != (*itr)->GetNbins()) {
            LOG(ERROR) << "Number of reparametrized bins and the number of bins for region " << regName << " do not match\n";
            LOG(ERROR) << "  - Number of bins: " << (*itr)->GetNbins() << "\n";
            LOG(ERROR) << "  - Number of reparametrized bins: " << nBins << "\n";
            exit(EXIT_FAILURE);
        }

        /// we need to check every bin for its dependencies and find the unique ones that we pass to the WS
        for (int ibin = 0; ibin < nBins; ++ibin) {
            std::string dependencies("");
            const auto dep = Common::processString(sf->fExpression.at(ibin).second);

            for (const auto& i : dep) {
                if (!ws->factory(i.first.c_str())) {
                    // get the [nominal, min, max] part
                    std::string range("[");
                    for (const auto ivalue : i.second) {
                        range += std::to_string(ivalue) + ",";
                    }
                    range.resize(range.size() - 1);
                    range += "]";
                    ws->factory((i.first + range).c_str());
                    mc->SetParametersOfInterest(RooArgSet(*mc->GetParametersOfInterest(), *ws->var(i.first.c_str())));
                }
                dependencies += i.first + ",";
            }
            if (dependencies.size() == 0) {
                LOG(ERROR) << "Error in the dependencies parsing\n";
                exit(EXIT_FAILURE);
            }

            // remove the trailing ","
            dependencies.resize(dependencies.size() - 1);

            LOG(INFO) << "Adding the extra dependencies for scale factor " << sf->fName << " in bin " << ibin << " : " << dependencies << "\n";
            ws->defineSet((sf->fName + "_dependencies_"+std::to_string(ibin)).c_str(), dependencies.c_str());

            LOG(INFO) << "Using the following formula: " << sf->fExpression.at(ibin).first << "\n";
            RooFormulaVar sfFunc(("sfFunc_" + sf->fName + "_bin_" + std::to_string(ibin)).c_str(), (sf->fExpression.at(ibin).first).c_str(), RooArgList(*ws->set((sf->fName + "_dependencies_"+std::to_string(ibin)).c_str())));

            auto tmp = ws->var(("gamma_" + sf->fName + "_bin_" + std::to_string(ibin)).c_str());
            if (!tmp) {
                LOG(ERROR) << "Cannot find parameter: gamma_" << sf->fName << "_bin_" << ibin << " in the WS\n";
                exit(EXIT_FAILURE);
            }

            xRooNode node(*tmp);
            node.Replace(sfFunc);
        }
    }
}

//____________________________________________________________________________________
//
void FitUtils::SetStartingAndConstParams(RooStats::ModelConfig* mc, const FitUtils::FittingOptions& options) {

    std::vector<RooRealVar*> pois = FitUtils::GetVectorPOI(mc);

    if (pois.empty()) {
        LOG(WARNING) << "No POIs found!\n";
    }

    for (auto poi : pois) {
        poi->setConstant(options.constPOI);

        const std::string name = poi->GetName();
        double value(0);
        auto it = std::find_if(options.poiVals.begin(), options.poiVals.end(),
            [&name](const std::pair<std::string, double>& element){ return element.first == name;});

        if (it != options.poiVals.end()) {
            value = it->second;
        }
        poi->setVal(value);

        // randomize the POI
        if(!options.constPOI && options.randomize) {
            poi->setVal( value + options.randomValue*(gRandom->Uniform(2)-1.));
        }

        // check if it is not fixed
        auto itr = std::find(options.constNPs.begin(), options.constNPs.end(), name);
        if (itr != options.constNPs.end()) {
            const std::size_t index = std::distance(options.constNPs.begin(), itr);
            const double constVal = options.constNPvalues.at(index);
            poi->setVal(constVal);
            poi->setConstant(1);
        }

        LOG(INFO) << "POI: " << name << "\n";
        LOG(INFO) << "   -> Constant POI : " << poi->isConstant() << "\n";
        LOG(INFO) << "   -> Initial value of POI : " << poi->getVal() << "\n";
    }

    if (mc->GetNuisanceParameters()) {
        for(auto var_tmp : *mc->GetNuisanceParameters()) {
            RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
            const std::string np = var->GetName();
            bool found = false;
            //
            // first check if all systs, norm and gammas should be set to constant
            if((np.find("gamma_stat") != std::string::npos || np.find("gamma_shape_stat") != std::string::npos) && options.noGammas) {
                LOG(DEBUG) << "setting to constant : " << np << " at value " << var->getVal() << "\n";
                var->setConstant(1);
                var->setVal(1);
                found = true;
            }
            else if((np.find("alpha_") != std::string::npos || (np.find("gamma_shape") != std::string::npos && np.find("gamma_shape_stat") == std::string::npos)) && options.noSystematics) {
                var->setConstant(1);
                var->setVal(0);
                LOG(DEBUG) << "setting to constant : "<< np << " at value " << var->getVal() << "\n";
                found = true;
            }
            else if(np.find("alpha_") == std::string::npos && np.find("gamma_") == std::string::npos) {
                // is NFs, check if it is constant
                auto itr = std::find(options.constNFs.begin(), options.constNFs.end(), np);
                if (itr != options.constNFs.end()) {
                    LOG(DEBUG) << "setting to constant : " << np << " at value " << var->getVal() << "\n";
                    var->setConstant(1);
                    found = true;
                }
            }
            if(found) continue;
            //
            // loop on the NP specified to have custom starting value
            for(std::size_t inp = 0; inp < options.initialNPs.size(); ++inp){
                if(np == ("alpha_"+options.initialNPs.at(inp)) || np == options.initialNPs.at(inp)) {
                    var->setVal(options.initialNPvalues.at(inp));
                    LOG(INFO) << " ---> Setting " << options.initialNPs.at(inp) << " to "  << options.initialNPvalues.at(inp) << "\n";
                    found = true;
                    break;
                }
            }
            //
            // loop on the NP specified to be constant - This should be after setting to initial value
            for(std::size_t inp = 0; inp < options.constNPs.size(); ++inp) {
                if( np == ("alpha_"+options.constNPs.at(inp)) || np == options.constNPs.at(inp)
                    || np == ("gamma_"+options.constNPs.at(inp))
                ){
                    LOG(INFO) << "setting to constant : " << np << " at value " << options.constNPvalues.at(inp) << "\n";
                    var->setVal(options.constNPvalues.at(inp));
                    var->setConstant(1);
                    found = true;
                    break;
                }
            }
            if(!found){
                if(np.find("alpha_") != std::string::npos) {   // for syst NP
                    if(options.randomize) var->setVal(options.randomValue*(gRandom->Uniform(2)-1.));
                    else                  var->setVal(0);
                    var->setConstant(0);
                } else if (np.find("gamma_") != std::string::npos) {  // gammas, dont randomize
                    var->setVal(1);
                } else { // norm factors
                    if(options.randomize) var->setVal(1 + options.randomValue*(gRandom->Uniform(2)-1.));
                    else                  var->setVal(1);
                }
            }
        }
    }
}

//____________________________________________________________________________________
//
std::vector<RooRealVar*> FitUtils::GetVectorPOI(const RooStats::ModelConfig* model) {
    std::vector<RooRealVar*> result;
    for (auto var_tmp : *model->GetParametersOfInterest()) {
        auto tmp = dynamic_cast<RooRealVar*>(var_tmp);
        if (!tmp) {
            LOG(ERROR) << "Cannot find the parameter of interest!\n";
            result.clear();
            return result;
        }

        const std::string name = tmp->GetName();
        auto it = std::find_if(result.begin(), result.end(),
            [&name](const RooRealVar* var){return name == var->GetName();});

        if (it == result.end()) {
            result.emplace_back(tmp);
        }
    }
    return result;
}

//____________________________________________________________________________________
//
double FitUtils::SaturatedModelGoF(std::shared_ptr<xRooNode> node, const std::string& dataName) {
    if (!node) {
        LOG(ERROR) << "Passed node is nullptr!\n";
        return -9999;
    }

    auto data = node->nll(dataName.c_str());

    const double saturatedNLL = data.saturatedMainTermVal() + data.saturatedConstraintTermVal();
    const double nll = data.get()->getVal();
    const int dof = data.mainTermNdof();
    double twoTimesDiff = 2.*(nll - saturatedNLL);
    // this can happen for asimov due to numerics
    if (twoTimesDiff < 0) twoTimesDiff = 0.;
    const double prob = TMath::Prob(twoTimesDiff, dof);

    LOG(DEBUG) << "Goodness-of-fit evaluation\n";
    LOG(DEBUG) << "\tBest fit NLL: " << nll << "\n";
    LOG(DEBUG) << "\tSaturated model NLL: " << saturatedNLL << "\n";
    LOG(DEBUG) << "\tNDOF: " << dof << "\n";

    return prob;
}

//____________________________________________________________________________________
//
xRooNode FitUtils::GetReducedOrSingleSample(std::shared_ptr<xRooNode> node, const std::string& samples) {
    if (!node) {
        LOG(ERROR) << "Passed node is nullptr\n";
        exit(EXIT_FAILURE);
    }

    const bool multipleSample = samples.find(',') != std::string::npos;
    if (multipleSample) {
        return node->at("samples")->reduced(samples);
    }

    auto sample = node->at("samples")->find(samples);
    if (!sample) {
        LOG(ERROR) << "Sample: " << samples << " not found in the WS\n";
        exit(EXIT_FAILURE);
    }

    return *sample;
}

//____________________________________________________________________________________
//
RooArgSet FitUtils::GetNPsWithoutNFs(const RooStats::ModelConfig* mc,
                                     const std::vector<std::shared_ptr<NormFactor> >& nfs,
                                     const std::vector<std::string>& pois) {
    RooArgSet result;

    if (!mc->GetNuisanceParameters()) return result;

    for (const auto& inp : *mc->GetNuisanceParameters()) {
        const std::string name = inp->GetName();
        auto* tmp = static_cast<RooRealVar*>(inp);
        auto itr = std::find_if(nfs.begin(), nfs.end(), [&name](const auto& element){return element->fName == name;});
        auto itrPOI = std::find(pois.begin(), pois.end(), name);

        if (itr != nfs.end()) continue;
        if (itrPOI != pois.end()) continue;

        result.add(*tmp);
    }

    return result;
}

//____________________________________________________________________________________
//
double FitUtils::GetPropagatedCovariance(const RooAbsReal &f, const RooRealVar &p, const RooFitResult &fr, const RooArgSet &nset) {
    RooArgSet allParamsInAbsReal;
    f.getParameters(&nset, allParamsInAbsReal);

    const bool isGamma = Common::StartsWith(p.GetName(), "gamma_");

    const auto* paramPrefit = static_cast<RooRealVar*>(fr.floatParsInit().find(p.GetName()));
    if (!paramPrefit) {
        LOG(ERROR) << "Cannot find parameter: " << p.GetName() << " in the prefit parameter list\n";
        return -1;
    }

    const double prefit = paramPrefit->getError();

    // this will be the list of parameters the function depends on
    RooArgList paramList;
    for(auto * rrvFitRes : static_range_cast<RooRealVar*>(fr.floatParsFinal())) {

       auto rrvInAbsReal = static_cast<RooRealVar const*>(allParamsInAbsReal.find(*rrvFitRes));

       // If this RooAbsReal is a RooRealVar in the fit result, we don't need to
       // propagate anything and can just return the error in the fit result
       double pError = p.getError();
       if (isGamma) pError /= prefit;
       if(rrvFitRes->namePtr() == f.namePtr()) return fr.correlation(f,p) * pError*rrvFitRes->getError();;

       // Strip out parameters with zero error
       if (rrvFitRes->getError() <= std::abs(rrvFitRes->getVal()) * std::numeric_limits<double>::epsilon()) continue;

       // Ignore parameters in the fit result that this RooAbsReal doesn't depend on
       if(!rrvInAbsReal) continue;

       // Checking for float equality is a bad. We check if the values are
       // negligibly far away from each other, relative to the uncertainty.
       if(std::abs(rrvInAbsReal->getVal() - rrvFitRes->getVal()) > 0.01 * rrvFitRes->getError()) {
          LOG(ERROR) << "The parameters of the RooAbsReal don't have the same values as in the fit result! The logic of GetPropagatedCovariance is broken in this case\n";
          return -1.;
       }

       paramList.add(*rrvInAbsReal);
    }

    std::vector<double> plusVar;
    std::vector<double> minusVar;
    plusVar.reserve(paramList.size());
    minusVar.reserve(paramList.size());

    // Create std::vector of plus,minus variations for each parameter
    TMatrixDSym V(paramList.size() == fr.floatParsFinal().size() ?
        fr.covarianceMatrix() :
        fr.reducedCovarianceMatrix(paramList)) ;

    for (std::size_t ivar=0 ; ivar<paramList.size() ; ivar++) {

      auto& rrv = static_cast<RooRealVar&>(paramList[ivar]);

      double cenVal = rrv.getVal() ;
      double errVal = std::sqrt(V(ivar,ivar)) ;

      // Make Plus variation
      rrv.setVal(cenVal+errVal) ;
      plusVar.push_back(f.getVal(nset)) ;

      // Make Minus variation
      rrv.setVal(cenVal-errVal) ;
      minusVar.push_back(f.getVal(nset)) ;

      rrv.setVal(cenVal) ;
    }

    // Re-evaluate this RooAbsReal with the central parameters just to be
    // extra-safe that a call to `getPropagatedCovariance()` doesn't change any state.
    // It should not be necessary because thanks to the dirty flag propagation
    // the RooAbsReal is re-evaluated anyway the next time getVal() is called.
    // Still there are imaginable corner cases where it would not be triggered,
    // for example if the user changes the RooFit operation more after the error
    // propagation.
    f.getVal(nset);

    // Make std::vector of variations
    TVectorD F(plusVar.size()) ;
    for (std::size_t j=0 ; j<plusVar.size() ; j++) {
      F[j] = (plusVar[j]-minusVar[j]) * 0.5;
    }

    // create the sum we want as result:
    //   Cov(f,p_i) = sum_j F_j*rho_ij*err(p_i)
    //
    double pError = p.getError();
    if (isGamma) pError /= prefit;
    double sum = 0.;
    for (std::size_t i=0 ; i<paramList.size() ; i++) {
      sum += F[i] * fr.correlation(paramList[i].GetName(),p.GetName()) * pError;
    }

    return sum;
}

//____________________________________________________________________________________
//
void FitUtils::FixParametersInXRooFit(const xRooNode* node, const std::map<std::string, double>& toFix) {
    if (toFix.empty()) return;

    auto pars = node->pars();

    for (const auto& ifix : toFix) {
        auto par = pars[ifix.first]->get<RooAbsArg>();
        if (!par) {
            LOG(WARNING) << "Cannot find parameter: " << ifix.first << " will NOT fix it\n";
            continue;
        }

        LOG(INFO) << "Fixing parameter " << ifix.first << " to " << ifix.second << "\n";
        static_cast<RooRealVar*>(par)->setVal(ifix.second);
        par->setAttribute("Constant");
    }
}
