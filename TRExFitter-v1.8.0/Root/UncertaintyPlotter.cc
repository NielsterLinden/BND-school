#include "TRExFitter/UncertaintyPlotter.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"

#include "StyleUtils/TRExLabels.h"
#include "StyleUtils/TRExUtils.h"

#include "RooFormulaVar.h"
#include "RooRealVar.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"

#include <fstream>
#include <memory>

UncertaintyPlotter::UncertaintyPlotter() noexcept :
    m_outFolder("/set/me"),
    m_useTotalUncertainty(true),
    m_plotLabel("Internal"),
    m_CME("13 TeV"),
    m_Lumi("139 fb^{-1}"),
    m_isNormalized(false),
    m_xAxisLabel("x-axis label"),
    m_colours({kBlack,kViolet, kBlue, kCyan -3, kGreen+3, kRed, kMagenta+2, kCyan+1}),
    m_styles({1,2,3,5,8}),
    m_useCustomRange(false),
    m_customRange(1)
{
}

//__________________________________________________________________________________
//
void UncertaintyPlotter::PlotUncertainties(const TH1* hist,
                                           const Unfolding& unfolding,
                                           const FittingTool& fitTool,
                                           const std::string& wsName,
                                           const std::string& wsPath,
                                           const std::string& frName,
                                           const std::vector<std::shared_ptr<NormFactor> >& nfs) const {
    if (!hist) {
        LOG(ERROR) << "Histogram is nullptr\n";
        return;
    }

    // make a copy of the histograms to get he binning etc and set values to zero
    // but set error to total or total syst
    std::vector<std::string> freeParams;
    std::vector<double> total = m_useTotalUncertainty ? this->GetUncertainty(fitTool, "TOT_ERROR") : this->GetUncertainty(fitTool, "TOTSYST_ERROR");
    if (total.empty()) {
        LOG(ERROR) << "Total uncertainty is empty\n";
        return;
    }
    std::unique_ptr<TH1> copy(static_cast<TH1*>(hist->Clone()));
    {
        int ibin = 0;
        for (const auto& ipoi : m_pois) {
            copy->SetBinContent(ipoi.second + 1, 0.);
            copy->SetBinError(ipoi.second + 1, total.at(ibin));
            ++ibin;
            freeParams.emplace_back(ipoi.first);
        }
    }


    std::vector<double> reparametrisedStat(0);

    if (!m_reparametrisedBins.empty()) {
        for (const auto& inf : nfs) {
            const std::string name = inf->fName;
            auto itr = std::find(freeParams.begin(), freeParams.end(), name);
            if (itr != freeParams.end()) continue;

            freeParams.emplace_back(name);
        }
        std::unique_ptr<TFile> wsFile(TFile::Open(wsName.c_str(),"read"));
        if (!wsFile) {
            LOG(ERROR) << "Cannot open file with WS\n";
            exit(EXIT_FAILURE);
        }

        // - get fit-result
        std::unique_ptr<TFile> frFile(TFile::Open(frName.c_str(),"read"));
        if (!frFile) {
            LOG(ERROR) << "Cannot read the fit results file\n";
            return;
        }

        RooWorkspace* ws(wsFile->Get<RooWorkspace>(wsPath.c_str()));
        if (!ws) {
            LOG(ERROR) << "Cannot find the ws!\n";
            return;
        }

        std::unique_ptr<TFile> wsFileStat(TFile::Open(wsName.c_str(),"read"));
        RooWorkspace* wsStat(wsFileStat->Get<RooWorkspace>(wsPath.c_str()));
        for (const auto ibin : m_reparametrisedBins) {
            const std::string name = unfolding.fName + "_Bin_" + Common::IntToFixLenStr(ibin) + "_mu";
            const std::vector<double> values = FitUtils::CalculateExpressionRoofit(ws, frFile.get(), name, false, freeParams);
            if (values.size() != 4) {
                LOG(ERROR) << "Something went wrong with the propagation to the reparametrised bin\n";
                return;
            }

            copy->SetBinContent(ibin, 0);
            copy->SetBinError(ibin, 100*values.at(3)/values.at(0));

            const std::vector<double> valuesStat = FitUtils::CalculateExpressionRoofit(wsStat, frFile.get(), name, true, freeParams);
            if (valuesStat.size() != 4) {
                LOG(ERROR) << "Something went wrong with the propagation to the reparametrised bin for stat only\n";
                return;
            }
            reparametrisedStat.emplace_back(100*valuesStat.at(3)/valuesStat.at(0));
        }
        // need to close the file and open it again since somehting weird is happening with the WS...
        wsFile->Close();
    }

    copy->SetFillColor(kOrange-2);
    copy->SetLineColor(kOrange-2);
    copy->GetYaxis()->SetTitle("Fractional uncertainty [%]");
    copy->GetYaxis()->SetTitleSize(1.3*copy->GetYaxis()->GetTitleSize());
    copy->GetXaxis()->SetTitleSize(1.3*copy->GetXaxis()->GetTitleSize());
    copy->GetXaxis()->SetTitle(m_xAxisLabel.c_str());
    copy->SetMarkerSize(0);
    const double max = *std::max_element(total.begin(),total.end());
    if (m_useCustomRange) {
        copy->GetYaxis()->SetRangeUser(-m_customRange, m_customRange);
    } else {
        copy->GetYaxis()->SetRangeUser(-1.6*max,1.6*max);
    }

    // get stat only
    const std::vector<double> statValues = this->GetUncertainty(fitTool, "STAT_ERROR");
    std::unique_ptr<TH1> statHist(static_cast<TH1*>(hist->Clone()));
    {
        int ibin = 0;
        for (const auto& ipoi : m_pois) {
            statHist->SetBinContent(ipoi.second + 1, 0.);
            statHist->SetBinError(ipoi.second + 1, statValues.at(ibin));
            ++ibin;
        }
        for (std::size_t i = 0; i < m_reparametrisedBins.size(); ++i) {
            statHist->SetBinContent(m_reparametrisedBins.at(i), 0);
            statHist->SetBinError(m_reparametrisedBins.at(i), reparametrisedStat.at(i));
        }
    }
    statHist->SetFillColor(kOrange-3);
    statHist->SetLineColor(kOrange-3);
    statHist->SetMarkerSize(0);

    // for each category we need the relative impact in this bin
    const std::vector<std::vector<std::pair<std::string, double> > > impact = this->GetCategoryImpact(fitTool);

    std::vector<std::vector<double> > impactPropagated;

    for (const auto ibin : m_reparametrisedBins) {
        const std::string param = unfolding.fName + "_Bin_" + Common::IntToFixLenStr(ibin) + "_mu";
        impactPropagated.emplace_back(this->GetCategoryImpactReparametrised(fitTool, param, wsName, wsPath));
    }

    // regular category histos
    std::vector<std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D> > > categoryHistos;
    for (std::size_t icategory = 0; icategory < impact.at(0).size(); ++icategory) {
        categoryHistos.emplace_back(std::make_pair(static_cast<TH1D*>(copy->Clone()), static_cast<TH1D*>(copy->Clone())));
        if ((categoryHistos.back().first->GetNbinsX() != static_cast<int>(impact.size())) && (!m_isNormalized) && (m_reparametrisedBins.empty())) {
            LOG(ERROR) << "Mismatch between the bin numbers\n";
            exit(EXIT_FAILURE);
        }
    }

    int colourStyle(0);
    for (std::size_t icategory = 0; icategory < categoryHistos.size(); ++icategory) {
        int ibin = 0;
        for (const auto& ipoi : m_pois) {
            categoryHistos.at(icategory).first->SetBinError(ipoi.second+1, 0.);
            categoryHistos.at(icategory).first->SetBinContent(ipoi.second+1, impact.at(ibin).at(icategory).second);
            categoryHistos.at(icategory).second->SetBinError(ipoi.second+1, 0.);
            categoryHistos.at(icategory).second->SetBinContent(ipoi.second+1, -impact.at(ibin).at(icategory).second);
            ++ibin;
        }
        for (std::size_t bin = 0; bin < m_reparametrisedBins.size(); ++bin) {
            if (impactPropagated.at(bin).size() != impact.at(0).size()) {
                LOG(ERROR) << "Mismatch between the impact categories between the standard and reparametrised fits\n";
                LOG(ERROR) << "Size of impact: " << impact.at(0).size() << ", reparametrised: " << impactPropagated.at(bin).size() << "\n";
                return;
            }
            categoryHistos.at(icategory).first->SetBinError(m_reparametrisedBins.at(bin), 0.);
            categoryHistos.at(icategory).first->SetBinContent(m_reparametrisedBins.at(bin), impactPropagated.at(bin).at(icategory));
            categoryHistos.at(icategory).second->SetBinError(m_reparametrisedBins.at(bin), 0.);
            categoryHistos.at(icategory).second->SetBinContent(m_reparametrisedBins.at(bin), -impactPropagated.at(bin).at(icategory));
        }
        const std::size_t indexCol = colourStyle % m_colours.size();
        const std::size_t indexSty = colourStyle % m_styles.size();
        categoryHistos.at(icategory).first->SetLineColor(m_colours.at(indexCol));
        categoryHistos.at(icategory).first->SetLineStyle(m_styles.at(indexSty));
        categoryHistos.at(icategory).first->SetLineWidth(3);
        categoryHistos.at(icategory).first->SetFillColor(0);
        categoryHistos.at(icategory).second->SetLineColor(m_colours.at(indexCol));
        categoryHistos.at(icategory).second->SetLineStyle(m_styles.at(indexSty));
        categoryHistos.at(icategory).second->SetLineWidth(3);
        categoryHistos.at(icategory).second->SetFillColor(0);
        ++colourStyle;
    }

    // now do the actual plotting
    TCanvas canvas("","",1000,800);
    canvas.cd();

    copy->Draw("E2");
    if (m_useTotalUncertainty) {
        statHist->Draw("E2 same");
    }
    gPad->RedrawAxis();

    for (auto& icategory : categoryHistos) {
        icategory.first->Draw("HIST same");
        icategory.second->Draw("HIST same");
    }

    TLegend leg(0.2, 0.2, 0.8, 0.3);
    leg.SetNColumns(3);
    leg.AddEntry(copy.get(), m_useTotalUncertainty ? "Total unc." : "Systematic unc.", "f");
    if (m_useTotalUncertainty) {
        leg.AddEntry(statHist.get(), "Stat. unc.", "f");
    }
    for (std::size_t icat = 0; icat < categoryHistos.size(); ++icat) {
        leg.AddEntry(categoryHistos.at(icat).first.get(), impact.at(0).at(icat).first.c_str(), "l");
    }
    leg.Draw("same");

    std::unique_ptr<TLatex> latex;
    std::unique_ptr<TLatex> latex1;
    if (m_plotLabel != "none") {
        latex = std::make_unique<TLatex>();
        latex->SetNDC();
        latex->SetTextFont(73);
        latex->DrawLatex(0.2, 0.90, TRExFitter::EXPERIMENT_LABEL.c_str());
        latex1 = std::make_unique<TLatex>();
        latex1->SetNDC();
        latex1->SetTextFont(43);
        latex1->DrawLatex(0.28, 0.90, m_plotLabel.c_str());
    }
    myText(0.2,0.85,1,("#sqrt{s} = " + m_CME + ", " + m_Lumi).c_str());
    myText(0.2,0.82,1,(m_isNormalized ? "Normalised" : "Absolute"));

    Common::SaveCanvasAs(canvas, m_outFolder+"/Uncertainty_"+unfolding.fName);
}

//__________________________________________________________________________________
//
std::vector<std::vector<std::pair<std::string, double> > > UncertaintyPlotter::GetCategoryImpact(const FittingTool& fitTool) const {
    std::vector<std::vector<std::pair<std::string, double> > > result;

    std::vector<std::string> allPOIs;
    for (const auto& ipoi : m_pois) {
        allPOIs.emplace_back(ipoi.first);
    }
    allPOIs.insert(allPOIs.end(), m_otherPOIs.begin(), m_otherPOIs.end());

    // first need to get the reference
    for (const auto& ipoi : allPOIs) {

        auto itr = fitTool.GetErrDecompCategoryMap().find(ipoi);
        if (itr == fitTool.GetErrDecompCategoryMap().end()) {
            LOG(ERROR) << "Cannot find POI: " << ipoi << "\n";
            return result;
        }


        double nominal = static_cast<RooRealVar*>(fitTool.GetFitResult()->floatParsFinal().find(ipoi.c_str()))->getVal();
        if (nominal <= 0) {
            LOG(WARNING) << "POI: " << ipoi << " value <= 0, setting to 1e-3 for the fractional uncertainty plotting\n";
            nominal = 1e-3;
        }
        std::vector<std::pair<std::string, double> > oneNF;
        for (const auto& icategory : m_categories) {
            if (icategory == "FullSyst") continue;
            if (icategory == "NormFactors") continue;
            if (icategory == "") continue;
            if (Common::StartsWith(icategory, "gamma_")) continue;
            if (icategory == "Uncategorised") {
                LOG(WARNING) << "Uncategorised category found, ignoring.\n";
                continue;
            }
            std::string cat = icategory;
            if (cat == "Gammas") cat = "MC STAT";
            auto itrCat = itr->second.find(cat);
            if (itrCat == itr->second.end()) {
                LOG(ERROR) << "Cannot find category " << cat << "\n";
                exit(EXIT_FAILURE);
            }
            oneNF.emplace_back(std::make_pair(cat, 100*itrCat->second/nominal));
        }
        result.emplace_back(std::move(oneNF));
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<double> UncertaintyPlotter::GetCategoryImpactReparametrised(const FittingTool& fitTool,
                                                                        const std::string& param,
                                                                        const std::string& wsFile,
                                                                        const std::string& wsName) const {

    std::vector<double> result;

    std::unique_ptr<TFile> wsF(TFile::Open(wsFile.c_str(),"read"));
    if (!wsF) {
        LOG(ERROR) << "Cannot open file with WS\n";
        exit(EXIT_FAILURE);
    }

    RooWorkspace* ws = static_cast<RooWorkspace*>(wsF->Get(wsName.c_str()));
    if (!ws) return result;

    RooFormulaVar* muN(nullptr);
    for(auto *arg : ws->allFunctions()) {
        const std::string argName = arg->GetName();
        if(argName == param) {
            muN = &dynamic_cast<RooFormulaVar&>(*arg);
            break;
        }
    }

    if (!muN) {
        LOG(ERROR) << "Did not find the correct name: " << param << "\n";
        return result;
    }

    for(auto *par : fitTool.GetFitResult()->floatParsFinal()){
        ws->var( par->GetName() )->setVal( (static_cast<RooRealVar*>(par))->getVal() );
        ws->var( par->GetName() )->setError( (static_cast<RooRealVar*>(par))->getError() );
    }

    const double nominal = muN->getVal();

    for (const auto& icategory : m_categories) {
        if (icategory == "FullSyst") continue;
        if (icategory == "NormFactors") continue;
        if (icategory == "") continue;
        if (Common::StartsWith(icategory, "gamma_")) continue;
        if (icategory == "Uncategorised") continue;
        std::string cat = icategory;
        if (cat == "Gammas") cat = "MC STAT";

        double errorSq(0.);
        std::vector<std::string> variations;
        for (const auto& ivar : m_subCategoryMap) {
            if (icategory == "Gammas") {
                if (Common::StartsWith(ivar.first, "gamma_")) {
                    variations.emplace_back(ivar.first);
                }
            } else {
                if (ivar.second == icategory) {
                    variations.emplace_back(ivar.first);
                }
            }
        }

        for (const auto& ivar : variations) {

            // need to skip parameters that were pruned
            auto* par = ws->var(ivar.c_str());
            if (!par) continue;

            const auto *tmp = fitTool.GetFitResult()->floatParsFinal().find(ivar.c_str());
            if (!tmp) continue;

            const double err = FitUtils::GetPropagatedCovariance(*muN, *par, *fitTool.GetFitResult(), {});
            errorSq += err*err;
        }

        const double impact = 100.*std::sqrt(errorSq)/nominal;
        result.emplace_back(impact);
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<double> UncertaintyPlotter::GetUncertainty(const FittingTool& fitTool, const std::string& unc) const {
    std::vector<double> result;

    std::vector<std::string> allPOIs;
    for (const auto& ipoi : m_pois) {
        allPOIs.emplace_back(ipoi.first);
    }
    allPOIs.insert(allPOIs.end(), m_otherPOIs.begin(), m_otherPOIs.end());

    for (const auto& ipoi : allPOIs) {
        auto itr = fitTool.GetErrDecompCategoryMap().find(ipoi);
        if (itr == fitTool.GetErrDecompCategoryMap().end()) {
            LOG(ERROR) << "Cannot find POI: " << ipoi << "\n";
            return result;
        }

        auto* fitResult = fitTool.GetFitResult();
        if (!fitResult) {
            LOG(ERROR) << "Fit results is nullptr!\n";
            return result;
        }
        double nominal = static_cast<RooRealVar*>(fitTool.GetFitResult()->floatParsFinal().find(ipoi.c_str()))->getVal();
        if (nominal <= 0) {
            LOG(ERROR) << "POI: " << ipoi << " value <= 0, setting to 1e-3 for the fractional uncertainty plotting\n";
            nominal = 1e-3;
        }
        const double value = itr->second.find(unc)->second;
        result.emplace_back(100*value/nominal);
    }

    return result;
}
