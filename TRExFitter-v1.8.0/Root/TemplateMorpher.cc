#include "TRExFitter/TemplateMorpher.h"
#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/Morphing.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/SampleHist.h"

# include <regex>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TF12.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TSystem.h"

TemplateMorpher::TemplateMorpher(const std::string& name,
                                 const std::string& suffix) :
    fName(name),
    fSuffix(suffix),
    fFitDimensionality(TemplateMorpher::FIT_DIMENSIONALITY::ONE_DIMENSION),
    fFitFunction("QUADRATIC"),
    fMorphingSetting(nullptr)
{
    gSystem->mkdir((fName+"/Templates").c_str());
}

void TemplateMorpher::SetMorphingSetting(const Morphing* m) {
    fMorphingSetting = m;
    fFitFunction = m->GetMorphingFunction();
}

void TemplateMorpher::ProcessTemplatesAndPlot(const std::vector<std::unique_ptr<Region> >& regions) {
    LOG(INFO) << "Preparing histograms and formulae for the template morphing\n";
    std::ofstream out;
    out.open(fName+"/Templates/Parametrization.txt",std::ios::out);
    if (!out.is_open() || !out.good()) {
        LOG(ERROR) << "Cannot open file for reparametrisation output\n";
        exit(EXIT_FAILURE);
    }

    auto GetParamNames = [](const std::vector<std::pair<std::string, double> >& in) {
        std::vector<std::string> result;
        for (const auto& i : in) {
            result.emplace_back(i.first);
        }

        return result;
    };

    for (const auto& ireg : regions) {
        // get the list of targets for morphing
        std::vector<std::string> targets;
        for (const auto& isample : ireg->fSampleHists) {
            const auto& tmp = isample->fSample->fTemplateMorphing;
            if (tmp.empty()) continue;
            const std::string target = isample->fSample->fTemplateTarget;
            auto itr = std::find(targets.begin(), targets.end(), target);
            if (itr == targets.end()) {
                targets.emplace_back(target);
            }
        }

        std::vector<std::vector<std::string> > templateTypes;
        for (const auto& isample : ireg->fSampleHists) {
            const auto& tmp = isample->fSample->fTemplateMorphing;
            if (tmp.empty()) continue;
            const auto params = GetParamNames(tmp);
            auto itr = std::find(templateTypes.begin(), templateTypes.end(), params);
            if (itr == templateTypes.end()) {
                templateTypes.emplace_back(params);
            }
        }
        gSystem->mkdir((fName + "/Templates/" + ireg->fName).c_str());

        for (const auto& itarget : targets) {

            for (const auto& itype : templateTypes) {
                std::vector<std::shared_ptr<SampleHist> > morphSamples;
                std::string sampleName(""),sampleNameTemplate("");
                for (const auto& isample : ireg->fSampleHists) {
                    const auto& tmp = isample->fSample->fTemplateMorphing;
                    if (tmp.empty()) continue;
                    const auto params = GetParamNames(tmp);
                    if (params != itype) continue;
                    if (isample->fSample->fTemplateTarget != itarget && isample->fSample->fName != itarget) continue;
                    if (isample->fSample->fType != Sample::SampleType::GHOST) {
                        if (isample->fSample->fAlternativeRefInTemplateMorphing != ""){
                            sampleNameTemplate = isample->fSample->fAlternativeRefInTemplateMorphing;
                            sampleName = isample->fSample->fName;
                            continue;
                        } else {
                            sampleNameTemplate = isample->fSample->fName;
                            sampleName = isample->fSample->fName;
                        }
                    }
                    morphSamples.emplace_back(isample);
                }

                const auto params = this->ProcessAndPlotOneRegion(morphSamples, ireg->fName, sampleNameTemplate);

                /// @brief sampleName | bin | parameters
                const auto sampleParams = std::make_pair(sampleName, params);
                this->ParamsToString(ireg, sampleParams, &out);
            }
        }
    }
    out.close();
}

std::vector<std::vector<double> > TemplateMorpher::ProcessAndPlotOneRegion(const std::vector<std::shared_ptr<SampleHist> >& samples,
                                                                           const std::string& regName,
                                                                           const std::string& targetName) {

    if (samples.empty()) {
        LOG(ERROR) << "Samples passed are empty\n";
        exit(EXIT_FAILURE);
    }

    int nBins(0);

    std::size_t position(0);
    std::size_t nominalPosition(0);
    std::string param1("");
    std::string param2("");
    std::vector<std::pair<const TH1*, std::vector<std::pair<std::string, double> > > > histos;
    for (const auto& isample : samples) {
        const auto* tmp = isample->fHist.get();
        if (!tmp) {
            LOG(ERROR) << "Null ptr passed for sample: " << isample->fName << "\n";
            exit(EXIT_FAILURE);
        }
        nBins = tmp->GetNbinsX();
        if (isample->fSample->fType != Sample::SampleType::GHOST) nominalPosition = position;
        param1 = isample->fSample->fTemplateMorphing.at(0).first;
        if (fFitDimensionality == TemplateMorpher::FIT_DIMENSIONALITY::TWO_DIMENSIONS) {
            param2 = isample->fSample->fTemplateMorphing.at(1).first;
        }

        std::vector<std::pair<std::string, double> > templates;
        for (const auto& i : isample->fSample->fTemplateMorphing) {
            templates.emplace_back(i);
        }

        histos.emplace_back(std::make_pair(tmp, templates));
        ++position;
    }
    std::vector<std::vector<double> > params;

    for (int ibin = 1; ibin <= nBins; ++ibin) {
        /// @brief template histogram | templates | name of the template | value | nominal yield | yield error
        std::vector<std::tuple<std::vector<std::pair<std::string, double> >, double, double> > tuple;
        for (const auto& ihist : histos) {
            const double value = ihist.first->GetBinContent(ibin);
            const double error = ihist.first->GetBinError(ibin);
            tuple.emplace_back(std::make_tuple(ihist.second, value, error));
        }
        const double nominal = histos.at(nominalPosition).first->GetBinContent(ibin);

        if (fFitDimensionality == TemplateMorpher::FIT_DIMENSIONALITY::ONE_DIMENSION) {
            params.emplace_back(this->ProcessAndPlotOneBin1D(tuple, nominal, regName, ibin, param1, targetName));
        } else {
            params.emplace_back(this->ProcessAndPlotOneBin2D(tuple, nominal, regName, ibin, param1, param2, targetName));
        }
    }

    return params;
}

std::vector<double> TemplateMorpher::ProcessAndPlotOneBin1D(const std::vector<std::tuple<std::vector<std::pair<std::string, double> >, double, double> >& tuple,
                                                            const double nominal,
                                                            const std::string& regName,
                                                            const int bin,
                                                            const std::string& param1,
                                                            const std::string& targetName) const {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ex;
    std::vector<double> ey;

    for (const auto& itemplate : tuple) {
        x.emplace_back(std::get<0>(itemplate).at(0).second);
        ex.emplace_back(0.);
        y.emplace_back(std::get<1>(itemplate)/nominal);
        ey.emplace_back(std::get<2>(itemplate)/nominal);
    }

    const double min = *std::min_element(y.begin(), y.end());
    const double max = *std::max_element(y.begin(), y.end());

    TCanvas canvas("","",800,600);
    TGraphErrors graph(x.size(), x.data(), y.data(), ex.data(), ey.data());
    graph.GetXaxis()->SetTitle(param1.c_str());
    graph.GetYaxis()->SetTitle("Yield/Nominal template");
    graph.GetYaxis()->SetTitleOffset(1.5*graph.GetYaxis()->GetTitleOffset());
    graph.GetYaxis()->SetRangeUser(0.98*min, 1.03*max);
    graph.Draw("AP");
    const double minX = *std::min_element(x.begin(), x.end());
    const double maxX = *std::max_element(x.begin(), x.end());
    TF1 func = this->FitFunc(minX, maxX);
    const int nPars = func.GetNpar();
    for (int i = 0; i < nPars; ++i) {
        // priority goes to a specific starting value
        if (fMorphingSetting->HasStartValue(targetName, regName, bin-1, i)) {
            const double value = fMorphingSetting->ParamStartValue(targetName, regName, bin-1, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << i << " start value to " << value << "\n";
            func.SetParameter(i, value);
        }
        // the HasStartValue map bin of -999 indicates use across all bins
        else if (fMorphingSetting->HasStartValue(targetName, regName, -999, i)) {
            const double value = fMorphingSetting->ParamStartValue(targetName, regName, -999, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << i << " start value to " << value << "\n";
            func.SetParameter(i, value);
        }
        // priority goes to a specific starting value
        if (fMorphingSetting->HasParamLimits(targetName, regName, bin-1, i)) {
            const auto& tmp = fMorphingSetting->ParamLimits(targetName, regName, bin-1, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << i << " limits to " << tmp.first << ", " << tmp.second << "\n";
            func.SetParLimits(i, tmp.first, tmp.second);
        }
        // the HasParamLimits map bin of -999 indicates use across all bins
        else if (fMorphingSetting->HasParamLimits(targetName, regName, -999, i)) {
            const auto& tmp = fMorphingSetting->ParamLimits(targetName, regName, -999, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << i << " limits to " << tmp.first << ", " << tmp.second << "\n";
            func.SetParLimits(i, tmp.first, tmp.second);
        }
    }
    func.SetLineColor(kRed);
    graph.Fit(&func, "RQ");

    const double chi2 = func.GetChisquare();
    const int ndof = func.GetNDF();

    if (ndof > 0 && (chi2/ndof > fMorphingSetting->GetChi2WarningLimit())) {
        LOG(WARNING) << "Chi^2/NDF > " << fMorphingSetting->GetChi2WarningLimit() << " for region: " << regName << " in bin: " << bin << "\n";
    }

    TLegend leg(0.7,0.7,0.9,0.9);
    leg.AddEntry(&graph, "Template yields", "p");
    leg.AddEntry(&func, "Fit function", "l");
    leg.Draw("same");

    TLatex l;
    l.SetNDC();
    l.SetTextAlign(12);
    l.DrawLatex(0.3, 0.8,("Region: " + regName + ", bin: " + std::to_string(bin)).c_str());
    l.DrawLatex(0.3, 0.7,("Target: " + targetName).c_str());
    l.DrawLatex(0.3, 0.65,Form("#chi^{2}/NDF = %.2f/%i", chi2, ndof));

    const std::string printPath = fName + "/Templates/" + regName + "/" + targetName + "_Bin_" + std::to_string(bin);

    Common::SaveCanvasAs(canvas, printPath);

    std::vector<double> params(nPars);
    for (int i = 0; i < nPars; ++i) {
        params.at(i) = func.GetParameter(i);
    }

    return params;
}

std::vector<double> TemplateMorpher::ProcessAndPlotOneBin2D(const std::vector<std::tuple<std::vector<std::pair<std::string, double> >, double, double> >& tuple,
                                                            const double nominal,
                                                            const std::string& regName,
                                                            const int bin,
                                                            const std::string& param1,
                                                            const std::string& param2,
                                                            const std::string& targetName) const {

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> ex;
    std::vector<double> ey;
    std::vector<double> ez;

    for (const auto& itemplate : tuple) {
        x.emplace_back(std::get<0>(itemplate).at(0).second);
        ex.emplace_back(0.);
        y.emplace_back(std::get<0>(itemplate).at(1).second);
        ey.emplace_back(0.);
        z.emplace_back(std::get<1>(itemplate)/nominal);
        ez.emplace_back(std::get<2>(itemplate)/nominal);
    }

    const double min = *std::min_element(z.begin(), z.end());
    const double max = *std::max_element(z.begin(), z.end());
    TCanvas canvas("","",800,600);
    canvas.SetTopMargin(0.08);
    TGraph2DErrors graph(x.size(), x.data(), y.data(), z.data(), ex.data(), ey.data(), ez.data());
    graph.GetXaxis()->SetTitle(param1.c_str());
    graph.GetYaxis()->SetTitle(param2.c_str());
    graph.GetZaxis()->SetTitle("Yield/Nominal template");
    graph.GetZaxis()->SetRangeUser(0.98*min, 1.03*max);

    const double minX = *std::min_element(x.begin(), x.end());
    const double maxX = *std::max_element(x.begin(), x.end());
    const double minY = *std::min_element(y.begin(), y.end());
    const double maxY = *std::max_element(y.begin(), y.end());

    TF2 func = this->FitFunc(minX, maxX, minY, maxY);
    const int nPars = func.GetNpar();
    for (int i = 0; i < nPars; ++i) {
        // priority goes to a specific starting value        
        if (fMorphingSetting->HasStartValue(targetName, regName, bin-1, i)) {
            const double value = fMorphingSetting->ParamStartValue(targetName, regName, bin-1, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << (i) << " start value to " << value << "\n";
            func.SetParameter(i, value);
        }
        // the HasStartValue map bin of -999 indicates use across all bins
        else if (fMorphingSetting->HasStartValue(targetName, regName, -999, i)) {
            const double value = fMorphingSetting->ParamStartValue(targetName, regName, -999, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << (i) << " start value to " << value << "\n";
            func.SetParameter(i, value);
        }
        // priority goes to a specific starting value        
        if (fMorphingSetting->HasParamLimits(targetName, regName, bin-1, i)) {
            const auto& tmp = fMorphingSetting->ParamLimits(targetName, regName, bin-1, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << i << " limits to " << tmp.first << ", " << tmp.second << "\n";
            func.SetParLimits(i, tmp.first, tmp.second);
        }
        // the HasParamLimits map bin of -999 indicates use across all bins
        else if (fMorphingSetting->HasParamLimits(targetName, regName, -999, i)) {
            const auto& tmp = fMorphingSetting->ParamLimits(targetName, regName, -999, i);
            LOG(DEBUG) << "Region: " << regName << ", target: " << targetName << " bin: " << (bin-1) << " setting parameter: " << i << " limits to " << tmp.first << ", " << tmp.second << "\n";
            func.SetParLimits(i, tmp.first, tmp.second);
        }
    }
    func.SetLineColor(kRed);
    graph.Fit(&func, "RQ");
    func.Draw("surf1");
    graph.Draw("ap same");

    const double chi2 = func.GetChisquare();
    const int ndof = func.GetNDF();

    if (ndof > 0 && (chi2/ndof > fMorphingSetting->GetChi2WarningLimit())) {
        LOG(WARNING) << "Chi^2/NDF > " << fMorphingSetting->GetChi2WarningLimit() << " for region: " << regName << " in bin: " << bin << "\n";
    }

    TLegend leg(0.85,0.85,0.99,0.99);
    leg.AddEntry(&graph, "Template yields", "p");
    leg.AddEntry(&func, "Fit function", "l");
    leg.Draw("same");

    TLatex l;
    l.SetNDC();
    l.SetTextAlign(12);
    l.DrawLatex(0.1, 0.95,("Region: " + regName + ", bin: " + std::to_string(bin)).c_str());
    l.DrawLatex(0.1, 0.90,("Target: " + targetName).c_str());
    l.DrawLatex(0.1, 0.85,Form("#chi^{2}/NDF = %.2f/%i", chi2, ndof));

    const std::string printPath = fName + "/Templates/" + regName + "/" + targetName + "_Bin_" + std::to_string(bin);

    Common::SaveCanvasAs(canvas, printPath);

    if(fMorphingSetting->GetMake1DProjections()){
        //create and print 1D projections of 2D parametrisation
        for(int iparam=0; iparam < 2; iparam++){
            //only do it once per x-value!
            std::vector<double> processed_params;
            std::vector<double> curr_params;
            if(iparam == 0) {curr_params = x;}
            else {curr_params = y;}

            for(const auto & curr : curr_params){
                bool already_processed = false;
                for(const auto & old_x : processed_params){
                    if( std::abs(old_x - curr) < 1.e-4 ){
                        already_processed = true;
                    }
                }
                if(already_processed) continue;
                else {processed_params.push_back(curr); }
                std::unique_ptr<TF12> f12;
                if(iparam == 0) f12 = std::make_unique<TF12>("f12",&func,curr,"y");
                else f12 = std::make_unique<TF12>("f12",&func,curr,"x");
                //create a TGraphErrors object for the points with x = curr
                std::vector<double> proj_x, proj_y, proj_ex, proj_ey;
                for(int ip=0; ip < graph.GetN(); ip++){
                    double ipx,ipy,ipz;
                    if(iparam == 0) { graph.GetPoint(ip,ipx,ipy,ipz); }
                    else { graph.GetPoint(ip,ipy,ipx,ipz); }

                    if( std::abs(ipx - curr) < 1.e-4){
                        proj_x.push_back(ipy);
                        proj_y.push_back(ipz);
                        proj_ex.push_back(0.);
                        proj_ey.push_back(graph.GetErrorZ(ip));
                    }
                }

                double projx_min = *std::min_element(proj_x.begin(), proj_x.end());
                double projx_max = *std::max_element(proj_x.begin(), proj_x.end());
                double min1d = *std::min_element(proj_y.begin(), proj_y.end());
                double max1d = *std::max_element(proj_y.begin(), proj_y.end());
                if( f12->GetMinimum(projx_min,projx_max) < min1d) {min1d = f12->GetMinimum(projx_min,projx_max); }
                if( f12->GetMaximum(projx_min,projx_max) > max1d) {max1d = f12->GetMaximum(projx_min,projx_max); }
                TGraphErrors graph1d (proj_x.size(),proj_x.data(),proj_y.data(),proj_ex.data(),proj_ey.data());
                graph1d.GetYaxis()->SetRangeUser(0.98*min1d,1.03*max1d);
                graph1d.GetXaxis()->SetRangeUser(0.98*projx_min,1.03*projx_max);
                f12->GetYaxis()->SetRangeUser(0.98*min1d,1.03*max1d);
                f12->GetXaxis()->SetRangeUser(0.98*projx_min,1.03*projx_max);
                f12->SetLineColor(kRed);
                if(iparam == 0){f12->GetXaxis()->SetTitle(param2.c_str());}
                else {f12->GetXaxis()->SetTitle(param1.c_str());}
                f12->GetYaxis()->SetTitle("Template yields");
                f12->Draw();
                graph1d.Draw("P SAME");

                TLegend leg1d(0.7,0.7,0.9,0.9);
                leg1d.AddEntry(&graph1d, "Template yields", "p");
                leg1d.AddEntry(f12.get(), ("1D projection ("+param1+"="+std::to_string(curr)+")").c_str(), "l");
                leg1d.Draw("same");

                TLatex l1d;
                l1d.SetNDC();
                l1d.SetTextAlign(12);
                l1d.DrawLatex(0.3, 0.8,("Region: " + regName + ", bin: " + std::to_string(bin)).c_str());
                std::string strCurr = std::to_string(curr);
                std::replace(strCurr.begin(),strCurr.end(),'.','p');

                if(bin == 1) {
                    if(iparam == 0) gSystem->mkdir(( fName + "/Templates/" + regName + "/" + targetName + "_Proj1D_"+ param1  + strCurr ).c_str());
                    else gSystem->mkdir(( fName + "/Templates/" + regName + "/" + targetName + "_Proj1D_"+ param2  + strCurr ).c_str());
                }

                std::string printPath1d;
                if(iparam == 0) printPath1d = fName + "/Templates/" + regName + "/" + targetName + "_Proj1D_"+ param1  + strCurr + "/Bin_" + std::to_string(bin)+".";
                else printPath1d = fName + "/Templates/" + regName + "/" + targetName + "_Proj1D_"+ param2  + strCurr + "/Bin_" + std::to_string(bin);

                Common::SaveCanvasAs(canvas, printPath1d);
            }
        }
    }

    std::vector<double> params(nPars);
    for (int i = 0; i < nPars; ++i) {
        params.at(i) = func.GetParameter(i);
    }

    return params;
}

std::string TemplateMorpher::FitFunctionString() const {
    std::string upper = fFitFunction;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    if (fFitDimensionality == TemplateMorpher::FIT_DIMENSIONALITY::ONE_DIMENSION) {
        if (upper == "LINEAR") {
            return "[0]*x + [1]";
        } else if (upper == "QUADRATIC") {
            return "[0]*x*x + [1]*x + [2]";
        } else if (upper == "CUBIC") {
            return "[0]*x*x*x + [1]*x*x + [2]*x + [3]";
        } else if (upper == "QUARTIC") {
            return "[0]*x*x*x*x + [1]*x*x*x + [2]*x*x + [3]*x + [4]";
        } else if (upper == "COSINE") {
            return "[0]*cos([1] * x + [2]) + [3]";
        }
    } else {
        if (upper == "LINEAR") {
            return "[0] * x + [1] * y + [2]";
        } else if (upper == "QUADRATIC") {
            return "[0]*x*x + [1]*y*y + [2]*x*y + [3]*x + [4]*y + [5]";
        } else if (upper == "CUBIC") {
            return "[0]*x*x*x + [1]*x*x*y + [2]*x*y*y + [3]*y*y*y"
                   "+ [4]*x*x + [5]*x*y + [6]*y*y"
                   "+ [7]*x + [8]*y "
                   "+ [9]";
        } else if (upper == "QUARTIC") {
            return "[0]*x*x*x*x + [1]*x*x*x*y + [2]*x*x*y*y + [3]*x*y*y*y + [4]*y*y*y*y"
                   "+ [5]*x*x*x + [6]*x*x*y + [7]*x*y*y + [8]*y*y*y + [9]*x*x + [10]*x*y"
                   "+ [11]*y*y + [12]*x + [13]*y + [14]";
        }
	else if (upper == "COSINE") {
	    LOG(ERROR) << "`COSINE` fit function is not implemented for two-dimensional morphing.\n";
	}
    }

    return fFitFunction;
}

TF1 TemplateMorpher::FitFunc(const double minX, const double maxX) const {
    return TF1("fitFunc", this->FitFunctionString().c_str(), minX, maxX);
}

TF2 TemplateMorpher::FitFunc(const double minX, const double maxX, const double minY, const double maxY) const {
    return TF2("fitFunc", this->FitFunctionString().c_str(), minX, maxX, minY, maxY);
}

bool TemplateMorpher::ValidateFitFunction() const {
    // we need to make sure that the fit function contains N parameters from [0] to [N-1] without skipping any
    std::regex re(R"(\[(\d+)\])");
    const std::string fitfunction = this->FitFunctionString();
    std::sregex_iterator it(fitfunction.begin(), fitfunction.end(), re);
    std::sregex_iterator end;
    std::set<int> indices;
    while (it != end) {
        const int index = std::stoi((*it)[1].str());
        indices.insert(index);
        ++it;
    }
    if(indices.empty()) return true;

    const int min = *indices.begin();
    const int max = *indices.rbegin();
    if (min != 0) {
        return false;
    }
    if((int)indices.size() != (max + 1)) {
        return false;
    }
    // finally, check that the function actually compiles in ROOT
    TFormula formula("test", fitfunction.c_str());
    if(!formula.IsValid()) {
        return false;
    }
    const int ndim = formula.GetNdim();
    switch (fFitDimensionality) {
        case TemplateMorpher::FIT_DIMENSIONALITY::ONE_DIMENSION:
            return ndim == 1;
        case TemplateMorpher::FIT_DIMENSIONALITY::TWO_DIMENSIONS:
            return ndim == 2;
    }
    // this should never happen
    return false;
}

void TemplateMorpher::ParamsToString(const std::unique_ptr<Region>& region,
                                     const std::pair<std::string, std::vector<std::vector<double> > >& params,
                                     std::ofstream* out) const {

    std::string param1("");
    std::string param2("");
    std::vector<double> values1;
    std::vector<double> values2;
    double nominal1(0.);
    double nominal2(0.);
    for (const auto& isample : region->fSampleHists) {
        if (isample->fSample->fTemplateMorphing.empty()) continue;
        param1 = isample->fSample->fTemplateMorphing.at(0).first;
        values1.emplace_back(isample->fSample->fTemplateMorphing.at(0).second);
        if (isample->fSample->fType != Sample::SampleType::GHOST) nominal1 = isample->fSample->fTemplateMorphing.at(0).second;
        if (fFitDimensionality == TemplateMorpher::FIT_DIMENSIONALITY::TWO_DIMENSIONS) {
            param2 = isample->fSample->fTemplateMorphing.at(1).first;
            values2.emplace_back(isample->fSample->fTemplateMorphing.at(1).second);
            if (isample->fSample->fType != Sample::SampleType::GHOST) nominal2 = isample->fSample->fTemplateMorphing.at(1).second;
        }
    }

    const double min1 = *std::min_element(values1.begin(), values1.end());
    const double max1 = *std::max_element(values1.begin(), values1.end());
    const double extra1 = fMorphingSetting->GetRangeScale()*std::abs(max1 - min1);

    double min2(0);
    double max2(0);
    double extra2(0);
    if (fFitDimensionality == TemplateMorpher::FIT_DIMENSIONALITY::TWO_DIMENSIONS) {
        min2   = *std::min_element(values2.begin(), values2.end());
        max2   = *std::max_element(values2.begin(), values2.end());
        extra2 = fMorphingSetting->GetRangeScale()*std::abs(max2 - min2);
    }

    const std::string& regName = region->fName;
    std::size_t bin(0);

    const std::map<std::string, std::pair<std::string, std::string> > expressions = fMorphingSetting->GetExpressions();

    auto itrExpressions1 = expressions.find(param1);
    auto itrExpressions2 = expressions.find(param2);

    for (const auto& ibin : params.second) {
        std::string formula = this->FitFunctionString();
        for(unsigned i=0; i < ibin.size(); ++i){
            std::stringstream ss;
            ss << std::scientific << ibin.at(i);
            std::string placeholder = "[" + std::to_string(i) + "]";
            formula = Common::ReplaceString(formula, placeholder, ss.str());
	}
	// replace "x" and "y" by their respective variable names
	// only replace if they are not part of a word (e.g. exp, poly, ...)
	std::regex patternx(R"(\bx\b)");
    	formula = std::regex_replace(formula, patternx, param1);
	std::regex patterny(R"(\by\b)");
    	formula = std::regex_replace(formula, patterny, param2);
        formula = Common::ReplaceString(formula, " ", "");

        std::string parameters("");

        if (itrExpressions1 != expressions.end()) {
            formula = Common::ReplaceString(formula, param1, itrExpressions1->second.first);
            parameters = itrExpressions1->second.second;
        } else {
            parameters = param1 + "[" + std::to_string(nominal1) + "," + std::to_string(min1-extra1) + "," + std::to_string(max1+extra1) + "]";
        }
        if (param2 != "") {
            if (itrExpressions2 != expressions.end()) {
                formula = Common::ReplaceString(formula, param2, itrExpressions2->second.first);
                parameters += "," + itrExpressions2->second.second;
            } else {
                parameters += "," + param2 + "[" + std::to_string(nominal2) + "," + std::to_string(min2-extra2) + "," + std::to_string(max2+extra2) + "]";
            }
        }

        *out << regName << " " << params.first << " " << bin << " " << formula << " " << parameters << "\n";
        LOG(DEBUG) << "Region " << regName << " " << "sample: " << params.first << " bin: " << bin << ", formula: " << formula << "\n";
        ++bin;
    }
}
