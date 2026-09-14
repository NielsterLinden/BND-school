// Class include
#include "TRExFitter/EFTProcessor.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/EFTConfig.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "TFitResult.h"
#include "Math/WrappedMultiTF1.h"
#include "TPRegexp.h"

// c++ includes
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <memory>

//_____________________________________________________________________________
//
EFTProcessor::EFTProcessor() :
    fName(""),
    fTitle(""),
    fRefMap{},
    fExtrapMap{},
    fOperators{},
    fParams{},
    fTitles{},
    fPowers{},
    fOrder{},
    fEFTVariations{},
    fProduceShapeFactorParametrization(false)
{
}

EFTProcessor::EFTProcessor(const std::string &param, const std::string &title) :
    fName(param),
    fTitle(title),
    fRefMap{},
    fExtrapMap{},
    fOperators{},
    fParams{},
    fTitles{},
    fPowers{},
    fOrder{},
    fEFTVariations{},
    fProduceShapeFactorParametrization(false)
{
}

//__________________________________________________________________________________
//
void EFTProcessor::SetEFTOrder(const EFTConfig::EFTOrder& order) {
    fOrder = order;
}

//_____________________________________________________________________________
//
void EFTProcessor::Print() const {
    for (const auto &imap : fRefMap) {
        LOG(DEBUG) << "  - " + imap.first << "\n";
        for (const auto &eft_s : imap.second) {
            LOG(DEBUG) << "    - " << eft_s << "\n";
        }
    }
}

//_____________________________________________________________________________
//
std::pair<std::string, std::vector<int>> EFTProcessor::GetParametrisationAndTerms(const int op_idx_match) const {
    std::string eqn = "[0]";

    std::vector<int> term_vec(1, 0);
    int term_counter = 1;
    for (int term_idx = 0; term_idx < (int)fPowers.size(); term_idx++) {
        if (op_idx_match >= 0) {
            bool skip = false;
            for (int op_idx = 0; op_idx < (int)fPowers[term_idx].size(); op_idx++) {
                if (op_idx == op_idx_match) {
                    continue;
                }
                if (fPowers.at(term_idx).at(op_idx) > 0) {
                    skip = true;
                    break;
                }
            }
            if (skip) continue;
        }
        eqn = eqn + " + [" + (term_counter) + "]";

        term_counter += 1;
        term_vec.emplace_back(term_idx + 1);

        for (std::size_t op_idx = 0; op_idx < fPowers.at(term_idx).size(); op_idx++) {
            for (int pow_idx = 0; pow_idx < fPowers.at(term_idx).at(op_idx); pow_idx++) {
                if (op_idx_match < 0) {
                    eqn = eqn + "*x[" + op_idx + "]";
                } else {
                    eqn = eqn + "*x[0]";
                }
            }
        }
    }
    return std::make_pair(eqn, term_vec);
}

//_____________________________________________________________________________
//
std::string EFTProcessor::GetParametrisation() const
{
    std::string eqn = "[0]";
    for (std::size_t term_idx = 0; term_idx < fPowers.size(); term_idx++) {
        eqn = eqn + " + [" + (term_idx + 1) + "]";
        for (std::size_t op_idx = 0; op_idx < fPowers.at(term_idx).size(); op_idx++) {
            for (int pow_idx = 0; pow_idx < fPowers.at(term_idx).at(op_idx); pow_idx++) {
                eqn = eqn + "*x[" + op_idx + "]";
            }
        }
    }
    return eqn;
}

//_____________________________________________________________________________
//
void EFTProcessor::PrintCoeffMap() const {

    LOG(INFO) << "------------------------------------------------------------------------\n";
    LOG(INFO) << "f(x[]) = " << this->GetParametrisation() << "\n";
    LOG(INFO) << "\n";
    std::size_t i(0);
    for (const auto& iop : fOperators) {
        LOG(INFO) << "x[" << i << "] : " << iop << "\n";
        ++i;
    }
    LOG(INFO) << "\n";

    i = 0;
    LOG(INFO) << "[0] = 0\n";
    for (const auto& ipowerVec : fPowers) {
        LOG(INFO) << "[" << (i+1) << "] coefficient for:\n";
        std::size_t j = 0;
        for (const auto ipow : ipowerVec) {
    	    if (ipow == 0) {
	        ++j;
	        continue;
	    }
	    const auto tmp = fOperators.at(j) + "^" + std::to_string(ipow);
	    LOG(INFO) << "      " << tmp << "\n";
	    ++j;
        }
        ++i;
    }

    LOG(INFO) << "------------------------------------------------------------------------\n";
}

//_____________________________________________________________________________
//
void EFTProcessor::PrintValuesMap() const {
    LOG(INFO) << "------------------------------------------------------------------------\n";

    std::string params("");
    std::string operators("");

    int i(0);
    for (const auto& iop : fOperators) {
        operators += iop + "\t";
        params += "x["+std::to_string(i)+"]\t";
        ++i;
    }

    LOG(INFO) << "Parameters  : " << params << "\n";
    LOG(INFO) << "Coefficients: " + operators << "\n";
    for (const auto &eftVar : fEFTVariations) {
        std::string line("");
        for (const auto ivalue : eftVar.fValues) {
            line += std::to_string(ivalue) + "\t";
        }
        line += eftVar.fName;
        LOG(INFO) << "Parameter   : " << line << "\n";
    }

    std::vector<std::string> vec;
    for (const auto &EFTVar : fEFTVariations) {
        vec.emplace_back(EFTVar.fName);
    }
    LOG(INFO) << "------------------------------------------------------------------------\n";
}

//_____________________________________________________________________________
//
void EFTProcessor::Resize() {
    for (auto& ipower : fPowers) {
        ipower.resize(fOperators.size(), 0);
    }
    for (auto& ivariation : fEFTVariations) {
        ivariation.fValues.resize(fOperators.size(), 0.);
    }
}

//_____________________________________________________________________________
//
std::string EFTProcessor::MakeTitle(const std::string &instr) const {

    TString title = instr;
    for (const auto& param : fTitleMap)
    {
        title = title.ReplaceAll(param.first, param.second);
    }

    // Fix powers to look pretty in Latex
    TPRegexp r("(?<=\\^)\\d+");
    r.Substitute(title, "{$&}", "g");

    return title.Data();
}

//_____________________________________________________________________________
// this draws the control plots for each EFT param-sample combination
void EFTProcessor::DrawEFTInputs(const std::vector<std::unique_ptr<Region> >& Regions,
                                 const EFTConfig& eftConfig) const {

    LOG(DEBUG) << "SM Ref: " << fName << "\n";

    for (std::size_t par_idx = 0; par_idx < fParams.size(); par_idx++) {
        TString paramName = fParams.at(par_idx);
        paramName = paramName.ReplaceAll(" ", "");
        paramName = paramName.ReplaceAll("+", "_");
        paramName = paramName.ReplaceAll("^", "");
        paramName = paramName.ReplaceAll("*", "");

        LOG(DEBUG) << "  Parameters: " << fParams.at(par_idx) << "(" << paramName.Data() << ")\n";

        for (const auto &ireg : Regions) {
            if (ireg->fRegionType == Region::VALIDATION) continue;
            std::vector<std::unique_ptr<TH1> > EFT_hist_vec;
            const std::string regName = ireg->fName;
            LOG(DEBUG) << "   Region: " << regName << "\n";

            ////////////////////////
            // Get and draw SMRef and EFT histograms

            // Get SM Reference sample
            std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(fName);
            if (!SM_sh) {
                LOG(DEBUG) << "     - Skipping region\n";
                continue;
            }

            // Make output dir
            gSystem->mkdir((SM_sh->fFitName + "/EFT").c_str());
            gSystem->mkdir((SM_sh->fFitName + "/EFT/" + regName).c_str());
            gSystem->mkdir((SM_sh->fFitName + "/EFT/" + regName + "/" + fName).c_str());

            std::unique_ptr<TH1> SM_hist(static_cast<TH1*>(SM_sh->fHist->Clone()));
            SM_hist->SetDirectory(nullptr);
            const std::string SMTitle = SM_sh->fSample->fTitle;

            // Make SM hist EFT coeff=0 sample if exists
            for (const auto &EFTVar : fEFTVariations) {
                // Only keep relevant variations
                if (EFTVar.fParam != fParams.at(par_idx)) continue;
                std::shared_ptr<SampleHist> EFT_sh = ireg->GetSampleHist(EFTVar.fName);

                std::unique_ptr<TH1> EFT_hist(static_cast<TH1*>(EFT_sh->fHist->Clone()));
                EFT_hist->SetDirectory(nullptr);

                bool is_SM = true;

                for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
                    is_SM = (EFTVar.fValues.at(op_idx) == 0.);
                    if (!is_SM) break;
                }
                if (is_SM) {
                    LOG(DEBUG) << "   Found coeff==0 sample: " << EFTVar.fName << "\n";
                    for (int b = 0; b < SM_hist->GetNbinsX() + 1; b++) {
                        SM_hist->SetBinContent(b, EFT_hist->GetBinContent(b));
                        SM_hist->SetBinError(b, EFT_hist->GetBinError(b));
                    }
                    break;
                }
            }

            std::vector<Int_t> cols;
            cols.emplace_back(kRed + 1);
            cols.emplace_back(kRed - 4);
            cols.emplace_back(kGreen + 1);
            cols.emplace_back(kGreen - 4);
            cols.emplace_back(kBlue + 1);
            cols.emplace_back(kBlue - 4);
            cols.emplace_back(kMagenta + 1);
            cols.emplace_back(kMagenta - 4);
            //@TODO be more clever about this...

            // Plotting
            TCanvas c1("","",800,600);
            std::unique_ptr<TPad> p0(nullptr);
            std::unique_ptr<TPad> p1(nullptr);
            if (eftConfig.GetEFTPlotRatio()) {
                p0 = std::make_unique<TPad>("p0","p0",0,0.35,1,1);
                p1 = std::make_unique<TPad>("p1","p1",0,0,1,0.35);
                p0->SetBottomMargin(0);
                p1->SetTopMargin(0);
                p1->SetBottomMargin(0.3);
                p0->Draw();
                p1->Draw();
                p0->cd();
            }

            // Plot SM hist
            SM_hist->SetLineWidth(3);
            SM_hist->SetFillColor(0);
            SM_hist->SetLineColor(1);
            SM_hist->GetXaxis()->SetTitle(SM_sh->fVariableTitle.c_str());
            SM_hist->GetYaxis()->SetTitle("Number of events");
            SM_hist->GetYaxis()->SetTitleOffset(2.);

            TLegend leg(0.6, 0.6, 0.9, 0.87);
            leg.SetBorderSize(0);
            leg.SetNColumns(2);

            Int_t n = 0;
            double ymax = -999;
            double ymin = 99999999;

            //@TODO: make order coherent for legend
            // Cycle through all EFT variations and plot - each corresponds to a differrent sample type
            for (const auto &EFTVar : fEFTVariations) {
                // Only keep relevant variations
                if (EFTVar.fParam != fParams.at(par_idx)) continue;

                LOG(DEBUG) << "    Working on: " << EFTVar.fName << "\n";

                std::shared_ptr<SampleHist> EFT_sh = ireg->GetSampleHist(EFTVar.fName);

                EFT_hist_vec.emplace_back(std::unique_ptr<TH1>(static_cast<TH1*>(EFT_sh->fHist->Clone())));
                EFT_hist_vec.back()->SetDirectory(nullptr);
                EFT_hist_vec.back()->SetLineColor(cols.at(n % cols.size()));
                EFT_hist_vec.back()->SetLineWidth(3);
                EFT_hist_vec.back()->SetFillStyle(0);

                // Make line style depend on sign
                bool foundPos = false;
                for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
                    if (EFTVar.fValues.at(op_idx) > 0.) {
                        EFT_hist_vec.back()->SetLineStyle(kDashed);
                        foundPos = true;
                        break;
                    }
                }
                if (!foundPos) EFT_hist_vec.back()->SetLineStyle(kDotted);

                // Get good y scale
                if (EFT_hist_vec.back()->GetMaximum() > ymax) ymax = EFT_hist_vec.back()->GetMaximum();
                if (EFT_hist_vec.back()->GetMinimum() < ymin) ymin = EFT_hist_vec.back()->GetMinimum();

                leg.AddEntry(EFT_hist_vec.back().get(), MakeTitle(EFTVar.fValueString).c_str(), "l");
                n++;
            }

            SM_hist->GetYaxis()->SetRangeUser(std::max(ymin * 0.75, 0.001), ymax * 1.7);
            SM_hist->Draw("HIST");
            for (const auto& ihist : EFT_hist_vec) {
                ihist->Draw("HIST SAME");
            }
            leg.AddEntry(SM_hist.get(), "SM prediction", "l");
            leg.Draw("same");

            TLatex tex{};
            tex.SetNDC();
            tex.DrawLatex(0.2, 0.79, TString::Format("%s, %s", MakeTitle(fParams.at(par_idx)).c_str(), SMTitle.c_str()));
            tex.DrawLatex(0.2, 0.72, regName.c_str());

            std::vector<std::unique_ptr<TH1> > ratios;
            std::vector<double> maxRatios;
            std::vector<double> minRatios;
            std::unique_ptr<TLine> line;
            maxRatios.emplace_back(SM_hist->GetMaximum());
            minRatios.emplace_back(SM_hist->GetMinimum());
            if (eftConfig.GetEFTPlotRatio()) {
                p0->RedrawAxis();

                p1->cd();
                for (const auto& ihist : EFT_hist_vec) {
                    ratios.emplace_back(std::unique_ptr<TH1>(static_cast<TH1*>(ihist->Clone())));
                    ratios.back()->Divide(SM_hist.get());
                    maxRatios.emplace_back(ratios.back()->GetMaximum());
                    minRatios.emplace_back(ratios.back()->GetMinimum());
                }

                const double max = *std::max_element(maxRatios.begin(), maxRatios.end());
                const double min = *std::min_element(minRatios.begin(), minRatios.end());

                double lineMin(0);
                double lineMax(0);
                if (!ratios.empty()) {
                    ratios.at(0)->GetYaxis()->SetRangeUser(0.8*min, 1.2*max);
                    lineMin = ratios.at(0)->GetBinLowEdge(1);
                    lineMax = ratios.at(0)->GetBinLowEdge(ratios.at(0)->GetNbinsX()+1);
                }

                for (std::size_t iratio = 0; iratio < ratios.size(); ++iratio) {
                    if (iratio == 0) {
                        ratios.at(iratio)->GetYaxis()->SetNdivisions(504);
                        ratios.at(iratio)->GetYaxis()->SetTitle("#frac{EFT}{SM}");
                        ratios.at(iratio)->GetYaxis()->CenterTitle();
                        ratios.at(iratio)->GetYaxis()->SetTitleOffset(2.5);
                        ratios.at(iratio)->Draw("HIST");
                    } else {
                        ratios.at(iratio)->Draw("HIST SAME");
                    }
                }

                line = std::make_unique<TLine>(lineMin, 1, lineMax, 1);
                line->SetLineColor(kBlack);
                line->SetLineWidth(2);

                line->Draw("same");

                p1->RedrawAxis();
            }
            c1.RedrawAxis();
            const std::string outName = SM_sh->fFitName + "/EFT/" + regName + "/" + fName + "/" + static_cast<std::string>(paramName);

            Common::SaveCanvasAs(c1, outName);
        }
    }
}

//_____________________________________________________________________________
// this fits the EFT parametrisation bin-by-bin
void EFTProcessor::FitEFTInputs(const std::vector<std::unique_ptr<Region> >& Regions,
                                const std::string& folder,
                                const std::string& sampleName,
                                const EFTConfig& eftConfig) {

    // Open file for quadratic fit result output - NormFactor usage
    const std::string outFileName = folder + "/EFT_Fit_results_" + sampleName + ".txt";
    std::ofstream NFTXTFile(outFileName.c_str());

    const std::string sfAllOutFileName = folder + "/SF_ALL_EFT_Fit_results.txt";
    const std::string sfLinearOutFileName = folder + "/SF_LINEAR_EFT_Fit_results.txt";

    std::unique_ptr<std::ofstream> SHFTXTFile(nullptr);
    std::unique_ptr<std::ofstream> SHFLinTXTFile(nullptr);

    if (fProduceShapeFactorParametrization) {
        SHFTXTFile = std::make_unique<std::ofstream>(sfAllOutFileName, std::ios::app);
        SHFLinTXTFile = std::make_unique<std::ofstream> (sfLinearOutFileName, std::ios::app);
    }

    // Find functional form of EFT parametrisation for header
    TString headereqn = "0";
    TString opEqnSHF = "";
    std::vector<std::string> opStringVec; // this vector is to put into the ShapeFactor expression without re-doing any logic
    for (std::size_t term_idx = 0; term_idx < fPowers.size(); term_idx++) {
        opEqnSHF = "";
        headereqn = headereqn + " ";
        for (std::size_t op_idx = 0; op_idx < fPowers.at(term_idx).size(); op_idx++) {
            for (int pow_idx = 0; pow_idx < fPowers.at(term_idx).at(op_idx); pow_idx++) {
                if(pow_idx==0){
                  headereqn = headereqn + fOperators.at(op_idx);
                }
                else{
                  headereqn = headereqn + "*" + fOperators.at(op_idx);
                }
                opEqnSHF = opEqnSHF + "*" + fOperators.at(op_idx);
            }
        }
        opStringVec.emplace_back(opEqnSHF);
    }

    LOG(DEBUG) << "Found Header Equation: " << headereqn.Data() << "\n";

    // this begins the header line
    NFTXTFile << "NORMFACTORS "<< headereqn << "\n";
    LOG(DEBUG) << "SM Ref: " << fName << "\n";

    for (const auto &ireg : Regions) {
        if (ireg->fRegionType == Region::VALIDATION) continue;

        const std::string rname = ireg->fName;
        LOG(DEBUG) << "  Region: " << rname << "\n";

        // Get SM Reference sample
        std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(fName);
        if (!SM_sh) {
            LOG(DEBUG) << "  - Skipping region\n";
            continue;
        }

        // Make output dir
        gSystem->mkdir((SM_sh->fFitName + "/EFT").c_str());

        std::unique_ptr<TH1> SM_hist(static_cast<TH1*>(SM_sh->fHist->Clone()));
        const std::string SMRef = fName;
        const std::string SMTitle = SM_sh->fSample->fTitle;

        // Make SM hist EFT coeff=0 sample if exists
        for (const auto &EFTVar : fEFTVariations) {
            std::shared_ptr<SampleHist> EFT_sh = ireg->GetSampleHist(EFTVar.fName);

            std::unique_ptr<TH1> EFT_hist(static_cast<TH1 *>(EFT_sh->fHist->Clone()));
            EFT_hist->SetDirectory(nullptr);

            // TODO: what is happening here?
            bool is_SM = true;
            for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
                is_SM = (EFTVar.fValues.at(op_idx) == 0.);
                if (!is_SM) break;
            }

            if (is_SM) {
                LOG(DEBUG) << "   Found coeff==0 sample: " << EFTVar.fName << "\n";
                for (int b = 0; b < SM_hist->GetNbinsX() + 1; b++) {
                    SM_hist->SetBinContent(b, EFT_hist->GetBinContent(b));
                    SM_hist->SetBinError(b, EFT_hist->GetBinError(b));
                }
                break;
            }
        }

        ////////////////////////
        // Fitting for each bin

        // Setup canvas for each 1D operator plot
        const double canv_x = std::min((SM_hist->GetNbinsX() + 1) * 200, 2000);
        std::vector<std::unique_ptr<TCanvas> > canvs;

        for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
            std::unique_ptr<TCanvas> tmp_c = std::make_unique<TCanvas>(TString::Format("canv%ld", op_idx), TString::Format("canv%ld", op_idx), canv_x, 400);
            canvs.emplace_back(std::move(tmp_c));
            canvs.at(op_idx)->cd();

            canvs.at(op_idx)->Divide(SM_hist->GetNbinsX() + 1, 1, 0, 0);

            canvs.at(op_idx)->cd(1);
            // Global EFT param, Region and Sample labels
            TLatex btex{};
            btex.SetNDC();
            btex.SetTextSize(20);
            btex.DrawLatex(0.17, 0.9, fOperators.at(op_idx).c_str());
            btex.DrawLatex(0.17, 0.79, SMTitle.c_str());
            btex.DrawLatex(0.17, 0.72, rname.c_str());
        }
        double bottom_margin = 0.3;

        double max_rel_eft = 1.;
        double min_rel_eft = 1.;

        bool empty_bin_flag = false;

        std::vector<std::unique_ptr<TGraphErrors> > g;
        std::vector<std::unique_ptr<TGraphErrors> > g_other;
        std::vector<std::unique_ptr<TF1> > fit;

        // Build dataset and fit for each bin
        for (int b = 0; b < SM_hist->GetNbinsX(); b++) {
            LOG(DEBUG) << Form("    Bin: %d", b) << "\n";

            const std::string regionBinName = Form("%s_bin%d", rname.c_str(), b);
            const std::string nfName = Form("Expression_muEFT_%s_%s", fName.c_str(), regionBinName.c_str());

            // Make entry for this RegionBin
            std::vector<double> tmp_values;
            auto it = fExtrapMap.insert({regionBinName, tmp_values});

            // Set initial fit  parameter values
            it.first->second.emplace_back(1.);
            for (std::size_t pow_idx = 0; pow_idx < fPowers.size(); pow_idx++) {
                it.first->second.emplace_back(0.);
            }

            // Skip empty bins
            if (SM_hist->GetBinContent(b + 1) <= 1e-6) {
                empty_bin_flag = true;
            }

            // Define dataset for fit
            LOG(DEBUG) << Form("N Variations: %ld N Operators: %ld", fEFTVariations.size() + 1, fOperators.size()) << "\n";
            ROOT::Fit::BinData EFTdata(fEFTVariations.size() + 1, fOperators.size());

            // Fill SM point
            std::vector<double> tmp_SM;
            for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) tmp_SM.emplace_back(0.);
            EFTdata.Add(tmp_SM.data(), 1., SM_hist->GetBinError(b + 1) / SM_hist->GetBinContent(b + 1));

            // Cycle through all EFT variations to get yields and extract relative changes
            for (const auto &EFTVar : fEFTVariations) {
                LOG(DEBUG) << "      Reading in EFT variation: " << EFTVar.fName << "\n";
                std::shared_ptr<SampleHist> EFT_sh = ireg->GetSampleHist(EFTVar.fName);

                double rel_eft = EFT_sh->fHist->GetBinContent(b + 1) / SM_hist->GetBinContent(b + 1);
                double rel_eft_err = rel_eft * std::sqrt(((EFT_sh->fHist->GetBinError(b + 1) / EFT_sh->fHist->GetBinContent(b + 1)) * (EFT_sh->fHist->GetBinError(b + 1) / EFT_sh->fHist->GetBinContent(b + 1))) +
                                                         ((SM_hist->GetBinError(b + 1) / SM_hist->GetBinContent(b + 1)) * (SM_hist->GetBinError(b + 1) / SM_hist->GetBinContent(b + 1))));

                EFTdata.Add(EFTVar.fValues.data(), rel_eft, rel_eft_err);
            }

            // Define function
            LOG(DEBUG) << "      - Fitting using Function: " << GetParametrisation() << "\n";

            std::string func = GetParametrisation();

            // suppress error message from wrong dimensions in TF1 that is harmless
            const auto currentIgnoreLever = gErrorIgnoreLevel;
            gErrorIgnoreLevel = 6001;
            TF1 f("f", func.c_str());
            gErrorIgnoreLevel = currentIgnoreLever;

            LOG(DEBUG) << Form("N par: %d", f.GetNpar()) << "\n";
            LOG(DEBUG) << Form("N dim: %d", f.GetNdim()) << "\n";

            ROOT::Math::WrappedMultiTF1 fitFunction(f, f.GetNdim());
            ROOT::Fit::Fitter fitter;
            fitter.SetFunction(fitFunction, false);

            LOG(DEBUG) << "      - Running fit...\n";
            fitter.LeastSquareFit(EFTdata);
            TFitResult fitResult = fitter.Result();
            LOG(DEBUG) << "      - Fit results:\n";
            if (TRExFitter::DEBUGLEVEL > 1) {
                fitResult.Print();
            }

            // Save results
            // This is the start of the individual lines for the NormFactor expressions (but logic is utilised for ShapeFactor expressions, so do not remove)
            NFTXTFile << nfName;
            std::vector<float> coeffVec; // this vector is to put into the ShapeFactor expression without re-doing any logic
            for (std::size_t pow_idx = 0; pow_idx < fPowers.size() + 1; pow_idx++) {
                float val = -999;
                if (empty_bin_flag) {
                    LOG(DEBUG) << "  - Empty bin identified - setting parameterisation to SM only. Check output text file.\n";
                    if (pow_idx == 0) {
                        val = 1.0;
                    } else  {
                        val = 0.0;
                    }
                } else {
                    val = fitResult.Parameter(pow_idx);
                }
                coeffVec.emplace_back(val);
                it.first->second.at(pow_idx) = val;
                // Write results to file
                NFTXTFile << "  " << val;
            }
            NFTXTFile << "\n";
            // This is the end of the lines for the NormFactor expressions

            // This is the start of the individual lines for the ShapeFactor expressions
            if (fProduceShapeFactorParametrization) {
	        *SHFTXTFile << rname << " " << sampleName << " " << b;
	        *SHFLinTXTFile << rname << " " << sampleName << " " << b;
                for (std::size_t term_idx = 0; term_idx < coeffVec.size(); term_idx++){

		    if (term_idx == 0){ // coefficient to SM part of expression (should be as close to 1.0 as possible)
		        *SHFTXTFile << " " << std::scientific << std::setprecision(6) << coeffVec[term_idx];
		        *SHFLinTXTFile << " " << std::scientific << std::setprecision(6) << coeffVec[term_idx];
		    }
		    else{ // coefficients to EFT parts of expression
		        // all EFT terms
		        *SHFTXTFile << "+" << std::scientific << std::setprecision(6) << coeffVec[term_idx] << opStringVec[term_idx-1];

                        // linear EFT terms only
                        if (std::count(opStringVec[term_idx-1].begin(), opStringVec[term_idx-1].end(), '*') < 2){ // less than two asterisks, so there must only be one EFT parameter in this term
			    *SHFLinTXTFile << "+" << std::scientific << std::setprecision(6) << coeffVec[term_idx] << opStringVec[term_idx-1];
                        }
		    }
		}
                *SHFTXTFile << " ";
                *SHFLinTXTFile << " ";
                // Here we put the ranges on the end
                for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
                    const std::pair<double, double> minMax = eftConfig.GetParMinMax(fOperators.at(op_idx));
                    *SHFTXTFile << "" << fOperators.at(op_idx) << std::defaultfloat << "[0.0," << minMax.first << "," << minMax.second <<"]" << std::scientific << std::setprecision(6);
                    *SHFLinTXTFile << "" << fOperators.at(op_idx) << std::defaultfloat << "[0.0," << minMax.first << "," << minMax.second <<"]" << std::scientific << std::setprecision(6);
                    if(op_idx != fOperators.size()-1){
                        *SHFTXTFile << ",";
                        *SHFLinTXTFile << ",";
                    }
                }
                *SHFTXTFile << "\n";
                *SHFLinTXTFile << "\n";
            }
            // this is the end of the lines for the ShapeFactor expressions


            // Plotting 1D projections

            // Get data points and split into operator of interest and others
            for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
                canvs.at(op_idx)->cd();
                canvs.at(op_idx)->cd(b + 2);

                // Values for TGraphs, one for pure single operator points and one for the rest
                std::vector<double> x;
                std::vector<double> x_err;
                std::vector<double> y;
                std::vector<double> y_err;

                std::vector<double> x_other;
                std::vector<double> x_other_err;
                std::vector<double> y_other;
                std::vector<double> y_other_err;

                double min_c = 999;
                double max_c = -999;

                for (std::size_t e = 0; e < fEFTVariations.size() + 1; e++) {
                    double tmp_y = 0;
                    double tmp_y_err = 0;
                    const double *coords = EFTdata.GetPoint(e, tmp_y, tmp_y_err);
                    if (coords[op_idx] > max_c) max_c = coords[op_idx];
                    if (coords[op_idx] < min_c) min_c = coords[op_idx];
                    if (tmp_y > max_rel_eft) max_rel_eft = tmp_y;
                    if (tmp_y < min_rel_eft) min_rel_eft = tmp_y;

                    bool is_pure_op = true;
                    for (std::size_t op_idx2 = 0; op_idx2 < fOperators.size(); op_idx2++) {
                        if (op_idx == op_idx2) continue;
                        if (coords[op_idx2] != 0) {
                            is_pure_op = false;
                            break;
                        }
                    }
                    if (is_pure_op) {
                        y.emplace_back(tmp_y);
                        y_err.emplace_back(1. / tmp_y_err);
                        x.emplace_back(coords[op_idx]);
                        x_err.emplace_back(0.);
                    } else {
                        y_other.emplace_back(tmp_y);
                        y_other_err.emplace_back(1. / tmp_y_err);
                        x_other.emplace_back(coords[op_idx]);
                        x_other_err.emplace_back(0.);
                    }
                }

                // Build graphs from points
                g.emplace_back(std::make_unique<TGraphErrors>(y.size(), x.data(), y.data(), x_err.data(), y_err.data()));
                g_other.emplace_back(std::make_unique<TGraphErrors>(y_other.size(), x_other.data(), y_other.data(), x_other_err.data(), y_other_err.data()));
                g.back()->SetTitle("");
                g.back()->SetMarkerStyle(8);
                g.back()->SetMarkerSize(1);

                g_other.back()->SetMarkerStyle(8);
                g_other.back()->SetMarkerSize(1);
                g_other.back()->SetMarkerColor(kGray);
                g_other.back()->SetLineColor(kGray);

                // Create 1D fit function
                std::pair<std::string, std::vector<int>> eqn_params = GetParametrisationAndTerms(op_idx);
                LOG(DEBUG) << "        - 1D fit results:\n";
                LOG(DEBUG) << "          " + eqn_params.first << "\n";

                fit.emplace_back(std::make_unique<TF1>("fit1d", eqn_params.first.c_str(), min_c, max_c));

                for (std::size_t par = 0; par < eqn_params.second.size(); par++) {
                    LOG(DEBUG) << Form("           [%ld] (=[%d]): %.2f", par, eqn_params.second.at(par), fitResult.Parameter(eqn_params.second.at(par))) << "\n";
                    fit.back()->FixParameter(par, fitResult.Parameter(eqn_params.second.at(par)));
                }

                fit.back()->SetLineColor(kRed);
                fit.back()->SetLineWidth(2);

                // Plot the data and fit
                TH1 *frame = gPad->DrawFrame(min_c - (0.2 * std::abs(min_c)), 0.8 * min_rel_eft, max_c + (0.2 * std::abs(max_c)), 1.1 * max_rel_eft, TString::Format("frame%d", b + 2));
                frame->GetXaxis()->SetTitle(MakeTitle(fOperators.at(op_idx)).c_str());
                bottom_margin = gPad->GetBottomMargin();

                g.back()->Draw("P");
                fit.back()->Draw("SAME");
                g_other.back()->Draw("P SAME");
                g.back()->Draw("P SAME");

                // Also draw the fitted equation
                TLatex gtex{};
                gtex.SetNDC();
                gtex.SetTextSize(14);
                gtex.DrawLatex(0.17, 0.90, ireg->fVariableTitle.c_str());
                gtex.DrawLatex(0.17, 0.85, TString::Format("Bin%d", b));
                gtex.SetTextSize(13);

                double ypos = 0.8;
                TString eqn = eqn_params.first;
                eqn = eqn.ReplaceAll(" ", "").ReplaceAll("x[0]", fOperators.at(op_idx));
                std::vector<std::string> eqn_terms = Common::Vectorize(MakeTitle(eqn.Data()), '+', true, false);
                LOG(DEBUG) << Form("N terms: %ld eqn %s", eqn_terms.size(), MakeTitle(eqn.Data()).c_str()) << "\n";
                for (std::size_t par = 0; par < eqn_params.second.size(); par++) {
                    std::string termstr = "";
                    TString tmpterm = eqn_terms.at(par);

                    if (par == 0) {
                        tmpterm = tmpterm.ReplaceAll(Form("[%ld]", par), "");
                        termstr = TString::Format("y = %.2f%s", fitResult.Parameter(eqn_params.second.at(par)), tmpterm.Data()).Data();
                    } else {
                        tmpterm = tmpterm.ReplaceAll(Form("[%ld]*", par), "");
                        termstr = TString::Format("  + %.2f*%s", fitResult.Parameter(eqn_params.second.at(par)), tmpterm.Data()).Data();
                    }
                    gtex.DrawLatex(0.17, ypos, termstr.c_str());
                    ypos -= 0.04;
                }
            }
        }

        for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
            // Need this hack to get y-axis shown in a nice way
            canvs.at(op_idx)->cd(1);
            gPad->Range(-10, -1, 10, 1);
            TGaxis yaxis(10, (-1 + 2 * bottom_margin), 10, 1., 0.8 * min_rel_eft, 1.1 * max_rel_eft, 512, "");
            yaxis.SetTitle(TString::Format("#sigma(%s)/#sigma(SM)", fOperators.at(op_idx).c_str()));
            yaxis.SetTitleSize(0.08);
            yaxis.SetLabelSize(0.07);
            yaxis.Draw();

            // Draw legend (need some dummy objects)
            TH1D hleg("", "", 1, 0, 1);
            hleg.SetLineColor(kRed);
            hleg.SetLineWidth(2);
            TGraph gleg{};
            gleg.SetTitle("");
            gleg.SetMarkerStyle(8);
            gleg.SetMarkerSize(1);

            TLegend bleg(0.1, 0.2, 0.9, 0.5);
            bleg.SetBorderSize(0);
            bleg.SetTextSize(18);
            bleg.AddEntry(&gleg, TString::Format("#frac{#sigma(%s)}{#sigma(SM)}", fOperators.at(op_idx).c_str()), "p");
            bleg.AddEntry(&hleg, "Fit", "l");
            bleg.Draw();

            // Go back and fix all y-axis ranges to same value
            for (int ip = 2; ip <= SM_hist->GetNbinsX() + 1; ip++) {
                canvs.at(op_idx)->cd(ip);
                TH1 *tmp = static_cast<TH1 *>(gPad->GetPrimitive("hframe"));
                if (tmp) {
                    tmp->SetMaximum(1.1 * max_rel_eft);
                    tmp->SetMinimum(0.8 * min_rel_eft);
                } else {
                    gPad->GetListOfPrimitives()->Print();
                }

                gPad->RedrawAxis();
            }
	    gSystem->mkdir(Form("%s/EFT/Fits/%s/%s/%s", SM_sh->fFitName.c_str(), rname.c_str(), SMRef.c_str(), fOperators.at(op_idx).c_str()), 1);
            Common::SaveCanvasAs(*canvs.at(op_idx), Form("%s/EFT/Fits/%s/%s/%s", SM_sh->fFitName.c_str(), rname.c_str(), SMRef.c_str(), fOperators.at(op_idx).c_str()));
        }

    } // region loop

    NFTXTFile << "\n";
    NFTXTFile.close();
    if (fProduceShapeFactorParametrization) {
        SHFTXTFile->close();
        SHFLinTXTFile->close();
    }
}

//_____________________________________________________________________________
// Apply mu NormFactors to each bin based on quadratic fit output.
void EFTProcessor::ApplyMuFactExpressions(const std::vector<std::unique_ptr<Region> >& regions,
                                          std::vector<std::shared_ptr<NormFactor> >& normFactors) const {

    for (const auto& ireg : regions) {
        if (ireg->fRegionType == Region::VALIDATION)             continue;

        const std::string rname = ireg->fName;

        //  Get SM Reference sample
        std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(fName);
        if (!SM_sh) continue;

        for (int ib = 0; ib < ireg->GetNbins(); ib++) {

            // Form unique name for this region-bin-sample combination
            const std::string regionBinName = Form("%s_bin%d", rname.c_str(), ib);
            const std::string nfName = Form("Expression_muEFT_%s_%s", fName.c_str(), regionBinName.c_str());

            auto itr = fExtrapMap.find(regionBinName);
            if (itr == fExtrapMap.end()) {
                LOG(ERROR) << "Cannot find " << regionBinName << " in the map\n";
                exit(EXIT_FAILURE);
            }

            // Add the expression based on the result of the fit
            TString eqn = GetParametrisation();

            std::string dependancies = "";
            for (std::size_t op_idx = 0; op_idx < fOperators.size(); op_idx++) {
                eqn = eqn.ReplaceAll(Form("x[%ld]", op_idx), fOperators.at(op_idx));

                // need to find the min and max value for operators
                auto itrNF = std::find_if(normFactors.begin(), normFactors.end(), [this, op_idx](const auto& nf){return nf->fName == this->fOperators.at(op_idx);});
                if (itrNF == normFactors.end()) {
                    LOG(ERROR) << "NormFactor: " << fOperators.at(op_idx) << " not found. Did you provide it in the config?\n";
                    exit(EXIT_FAILURE);
                }

                const double min = (*itrNF)->GetMin();
                const double max = (*itrNF)->GetMax();

                if (op_idx == 0) dependancies = TString::Format("%s[%f,%f,%f]", fOperators.at(op_idx).c_str(), 0.0, min, max);
                else dependancies = dependancies + "," + TString::Format("%s[%f,%f,%f]", fOperators.at(op_idx).c_str(), 0.0, min, max);
            }

            for (std::size_t pow_idx = 0; pow_idx < fPowers.size() + 1; pow_idx++) {
                eqn = eqn.ReplaceAll(Form("[%ld]", pow_idx), Form("%f", itr->second.at(pow_idx)));
            }

            LOG(DEBUG) << "    Defining EFT mu NormFactor Expression for: " << nfName << "\n";
            LOG(DEBUG) << "             " << eqn << ":" << dependancies << "\n";

            // Add to matching NF
            for (auto& norm : normFactors) {
                if (norm->fName != nfName) continue;
                LOG(DEBUG) << "      --> Updating NF " << nfName << "\n";
                norm->fExpression = std::make_pair(eqn, dependancies);
                // Need to name like this because formula is used later on from fTitle...
                const std::string newname = ModifyEFTParamsOrder(Form("Expression_%s", eqn.Data()));
                norm->fTitle = newname;
                auto itrNF = TRExFitter::SYSTMAP.find(norm->fName);
                if (itrNF == TRExFitter::SYSTMAP.end()) {
                    TRExFitter::SYSTMAP.insert({norm->fName, newname});
                } else {
                    itrNF->second = newname;
                }
                // nuis-par will contain the nuis-par of the norm factor the expression depends on FIXME
                const std::string newnpname = Form("Expression_%s", dependancies.c_str());
                norm->fNuisanceParameter = newnpname;
                auto itrNP = TRExFitter::NPMAP.find(norm->fName);
                if (itrNP == TRExFitter::NPMAP.end()) {
                    TRExFitter::NPMAP.insert({norm->fName, newnpname});
                } else {
                    itrNP->second = newnpname;
                }

                LOG(DEBUG) << "         --> Now NF " << norm->fTitle << "\n";
                break;
            }
        }
    }
}

//_____________________________________________________________________________
// Read quadratic fit results from file and pass them to relevant EFTProcessor objects
void EFTProcessor::ReadEFTFitResults(const std::vector<std::unique_ptr<Region> >& regions, const std::string& fileName) {

    LOG(DEBUG) << " Opening file \"" << fileName << "\"\n";

    std::ifstream NFTXTFile;
    NFTXTFile.open(fileName.c_str());

    if (!NFTXTFile.is_open()) {
        LOG(ERROR) << "Could not open the file \"" << fileName << "\"\n";
        return;
    }

    std::string input;
    std::string line;
    bool readingNF = false;
    std::string name;

    //
    // read file line by line
    while (std::getline(NFTXTFile, line)) {
        if (line == "") continue;
        if (line.find("NORMFACTORS") != std::string::npos){
            LOG(DEBUG) << "--------------------\n";
            LOG(DEBUG) << "Reading Norm Factors...\n";
            LOG(DEBUG) << "--------------------\n";
            readingNF = true;
            continue;
        }
        std::istringstream iss(line);
        if (readingNF) {
            iss >> input;
            if (input == "") { // leaving NF block
                readingNF = false;
            }
            while (input.find("\\") != std::string::npos) input = input.replace(input.find("\\"), 1, "");
            name = input;

            // Find appropriate SM Ref
            bool foundEntry = false;
            for (const auto &ireg : regions) {
                if (ireg->fRegionType == Region::VALIDATION) continue;
                const std::string rname = ireg->fName;

                // Get SM Reference sample
                std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(fName);
                if (!SM_sh) continue;

                for (int ib = 0; ib < ireg->GetNbins(); ib++) {
                    const std::string regionBinName = Form("%s_bin%d", rname.c_str(), ib);
                    const std::string nfName = Form("Expression_muEFT_%s_%s", fName.c_str(), regionBinName.c_str());
                    if (name != nfName) continue;
                    auto itr = fExtrapMap.find(regionBinName);
                    if (itr == fExtrapMap.end()) {
                        itr = (fExtrapMap.insert({regionBinName, {}})).first;
                    }

                    LOG(VERBOSE) << "  -> Found matching ExtrapMap\n";

                    double p;
                    int termCounter = 0;
                    while (iss >> p) {
                        if (fOrder == EFTConfig::EFTOrder::LIN) {
                            if (termCounter != 0) { // 0 term is SM
                                LOG(DEBUG) << "termCounter:    \t" << termCounter << "\n";
                                // loop over all EFT coefficients sum the number of powers
                                int powercounter = 0;
                                for (const auto& coefficientpower : fPowers.at(termCounter - 1)) {
                                    LOG(DEBUG) << "coefficientpower: " << coefficientpower << "\n";
                                    powercounter += coefficientpower;
                                }
                                LOG(DEBUG) << "powercounter: " << powercounter << "\n";
                                if (powercounter > 1) {
                                    p = 0.0; // removal of effects of non-linear terms
                                }
                            }
                        }
                        LOG(VERBOSE) << "Putting " << p << " into fExtrapMap[" << regionBinName << "]\n";
                        itr->second.push_back(p);
                        termCounter++;
                    }

                    foundEntry = true;
                    break;
                }
                if (foundEntry) break;
            }
        }
    }

    NFTXTFile.close();
}

// Modify the order of the EFT params just before applying the mufactorepxressions
std::string EFTProcessor::ModifyEFTParamsOrder(const std::string& tmpname) const {

    if (fOrder == EFTConfig::EFTOrder::ALL) {
        return tmpname;
    } // all orders, no modifications

    std::stringstream ss(tmpname);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    std::string newname = vstrings.at(0);

    if (fOrder == EFTConfig::EFTOrder::LIN) { // linear terms only
        for (unsigned int i = 1; i < vstrings.size(); i++) {
            if (vstrings.at(i).find("+") != std::string::npos) {
                continue;
            }
            std::string::difference_type n = std::count(vstrings.at(i).begin(), vstrings.at(i).end(), '*');
            if (static_cast<int>(n) < 2) { // only append the terms which are not quadratic or worse (i.e. without muliple "*" symbols)
                newname += " + " + vstrings.at(i);
            }
        }

        return newname;
    }

    // if the ordering doesn't make sense, just return the original string
    return tmpname;
}
