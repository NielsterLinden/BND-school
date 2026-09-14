// Class include
#include "TRExFitter/FitResults.h"

// framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/ShapeFactor.h"

// Blinder
#include "Blinder/Blinder/Blinder.h"

// ROOT includes
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

#include "RooFitResult.h"
#include "RooRealVar.h"

// Style stuff
#include "StyleUtils/TRExLabels.h"

//c++ includes
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

//__________________________________________________________________________________
//
FitResults::FitResults() :
    fCorrMatrix(nullptr),
    fOutFolder(""),
    fPOIPrecision(2) {
}

//__________________________________________________________________________________
//
FitResults::~FitResults(){
}

//__________________________________________________________________________________
//
void FitResults::AddNuisPar(NuisParameter *par){
    const std::string p = par->fName;
    auto itr = fNuisPar.find(p);
    if (itr != fNuisPar.end()) {
        LOG(WARNING) << "Parameter " << p << " already in the list\n";
        return;
    }
    fNuisPar.insert({p,std::shared_ptr<NuisParameter>(par)});
}

//__________________________________________________________________________________
//
double FitResults::GetNuisParValue(const std::string& p) const {
    return this->GetNuisanceParameter(p)->fFitValue;
}

//__________________________________________________________________________________
//
double FitResults::GetNuisParErrUp(const std::string& p) const {
    return this->GetNuisanceParameter(p)->fPostFitUp;
}

//__________________________________________________________________________________
//
double FitResults::GetNuisParErrDown(const std::string& p) const {
    return this->GetNuisanceParameter(p)->fPostFitDown;
}

//__________________________________________________________________________________
//
bool FitResults::ReadFromRootFile(const std::string& fileName) {

    TFile* file = Common::GetFile(fileName);
    if (!file) {
        LOG(ERROR) << "Cannot open ROOT file fit RooFitResults at " << fileName << "\n";
        return false;
    }

    std::unique_ptr<CorrelationMatrix> matrix = std::make_unique<CorrelationMatrix>();

    std::unique_ptr<RooFitResult> fr = Common::RooFitResultFromFile(file);

    std::vector<std::string> matrixNPs;

    for (auto var_tmp : fr->floatParsFinal()) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);

        // Not consider nuisance parameter being not associated to syst (yet)
        std::string name = var->GetName();
        name = Common::ReplaceString(name,"alpha_","");
        name = Common::ReplaceString(name,"gamma_","");

        const double value  = var->getVal(); // GetValue() return value in unit of sigma
        const double errorHi = var->getErrorHi();
        const double errorLo = var->getErrorLo();

        this->AddNuisPar(new NuisParameter(name));
        NuisParameter *np = this->GetNuisanceParameter(name);
        np->fFitValue = value;
        np->fPostFitUp = errorHi;
        np->fPostFitDown = errorLo;
        LOG(VERBOSE) << "Reading parameter: " << name << ", central value: " << value << ", up error: " << errorHi << ", down error: " << errorLo << "\n";

        matrix->AddNuisPar(name);

        matrixNPs.emplace_back(name);
    }

    std::unique_ptr<TH2> h2Dcorrelation(fr->correlationHist());
    h2Dcorrelation->SetDirectory(nullptr);
    matrix->Resize(h2Dcorrelation->GetNbinsX());
    for(int i = 0; i < h2Dcorrelation->GetNbinsX(); i++) {
        for(int j = 0; j < h2Dcorrelation->GetNbinsY(); j++) {
            matrix->SetCorrelation(matrixNPs.at(matrixNPs.size() - j - 1), matrixNPs.at(i), h2Dcorrelation->GetBinContent(i+1, j+1));
        }
    }

    fCorrMatrix = std::move(matrix);

    fRoofitFitResult = std::move(fr);

    return true;
}

//__________________________________________________________________________________
//
void FitResults::DrawNormFactors(const std::string &path,
                                 const std::vector < std::shared_ptr<NormFactor> > &normFactors, const std::vector<std::string>& blinded) const {
    double xmin = 1000.;
    double xmax = -1000.;
    double max = 0.;

    TGraphAsymmErrors g{};

    std::vector<std::unique_ptr<NuisParameter> > selected_norm_factors;

    for (const auto& inp : fNuisPar){
        NuisParameter* par = inp.second.get();

        // skip the blinded NPs
        if (std::find(blinded.begin(), blinded.end(), par->fName) != blinded.end()) continue;

        // skip hidden NPs
        if (Common::FindInStringVector(fNuisParToHide,par->fName)>=0) continue;

        bool isNorm = false;
        for (const auto& norm : normFactors) {
            if(norm->fName==par->fName) {
                isNorm = true;
                break;
            }
        }
        if (!isNorm) continue;
        g.SetPoint(selected_norm_factors.size(),par->fFitValue,2*selected_norm_factors.size()+1);
        g.SetPointEXhigh(selected_norm_factors.size(), par->fPostFitUp);
        g.SetPointEXlow( selected_norm_factors.size(),-par->fPostFitDown);

        if( par->fFitValue+par->fPostFitUp > xmax ) xmax = par->fFitValue+par->fPostFitUp;
        if( par->fFitValue+par->fPostFitDown < xmin ) xmin = par->fFitValue+par->fPostFitDown;

        std::unique_ptr<NuisParameter> nuis(new NuisParameter(par->fName));
        nuis->fFitValue =    par -> fFitValue;
        nuis->fPostFitUp =   par -> fPostFitUp;
        nuis->fPostFitDown = par -> fPostFitDown;
        nuis->fTitle =       par -> fTitle;
        selected_norm_factors.emplace_back(nuis.release());
        if(2*selected_norm_factors.size() > max)  max = 2*selected_norm_factors.size();
    }
    xmax *= (xmax<0 ? 0.5 : 1.5);
    xmin *= (xmin>0 ? 0.5 : 1.5);
    if(xmin>0) xmin = 0.;
    xmax += (xmax-xmin)*0.25;

    int lineHeight = 40;
    int offsetUp = 60;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy( "h_dummy_norm","h_dummy_norm",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l0;
    TBox b1, b2;
    if(((TString)path.c_str()).Contains("EFTParams"))l0 = TLine(0,0,0,max);
    else l0 = TLine(1,0,1,max);
    l0.SetLineStyle(7);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    g.Draw("psame");

    TLatex systs{};
    if (fPlotLabel != "none") {
        TRExLabelNew(0.04, 1. - 40./newHeight, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str(), kBlack, gStyle->GetTextSize(), 0.09);
    }
    systs.SetTextSize( systs.GetTextSize() );
    for(unsigned int i=0;i<selected_norm_factors.size();i++){
        systs.DrawLatex(xmax+(xmax-xmin)*0.05,2*i+0.75,(selected_norm_factors[i]->fTitle).c_str());
        systs.DrawLatex(xmax-(xmax-xmin)*0.30,2*i+0.75,
            Form(("%."+std::to_string(fPOIPrecision)+"f ^{%."+std::to_string(fPOIPrecision)+"f}_{%."+std::to_string(fPOIPrecision)+"f}").c_str(),selected_norm_factors[i]->fFitValue, selected_norm_factors[i]->fPostFitUp, selected_norm_factors[i]->fPostFitDown ) );
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    gPad->RedrawAxis();

    Common::SaveCanvasAs(c, path);
}

//__________________________________________________________________________________
//
void FitResults::DrawGammaShapePulls(const std::string& path,
                                     const std::vector<std::string>& blinded,
                                     const std::vector<std::shared_ptr<ShapeFactor> >& shapeFactors,
                                     const bool isShape) const {

    if (isShape && shapeFactors.empty()) return;

    double xmin =  1000.;
    double xmax = -1000.;
    double max  =  0.;

    TGraphAsymmErrors g{};

    // get a copy of the NPs (I want an actual copy of the NP objects that I can
    // manipulate without changing the originals.
    std::vector<NuisParameter> myNPs;

    // make the copies, dropping non-gamma NPs and blinded NPs
    for (const auto& inp : fNuisPar) {
        std::string name = inp.first;
        name = Common::ReplaceString(name,"gamma_","");
        if (std::find(blinded.begin(), blinded.end(), name) != blinded.end()) continue;

        if (isShape) {
            auto itr = std::find_if(shapeFactors.begin(), shapeFactors.end(), [&name](const auto& element){return name.find(element->fName) != std::string::npos;});
            if (itr == shapeFactors.end()) continue;
        } else {
            if (inp.second->fName.find("stat_") == std::string::npos && inp.second->fName.find("shape_") == std::string::npos) continue;
        }

        // make a copy, clean the name and save it
        NuisParameter theNP(*inp.second);

        std::string clean_name = theNP.fTitle;
        clean_name = Common::ReplaceString(clean_name, "stat_", "#gamma ");
        clean_name = Common::ReplaceString(clean_name, "shape_", "#gamma ");
        clean_name = Common::ReplaceString(clean_name, "#gamma #gamma ", "#gamma ");
        clean_name = Common::ReplaceString(clean_name, "_", " ");


        theNP.fTitle = Common::pad_trail(clean_name);

        myNPs.push_back(theNP);
    }

    if (myNPs.empty()) return;

    // now sort myNPs by the cleaned and padded names
    // do a to_lower in the sort_func
    // sort using a custom struct
    auto my_sort_func = [](const NuisParameter& a, const NuisParameter& b){return a.fTitle < b.fTitle;};

    std::sort( myNPs.begin(), myNPs.end(), my_sort_func );

    // start making the plot
    for(unsigned int i = 0; i < myNPs.size(); ++i )
    {
        g.SetPoint(i, myNPs[i].fFitValue, i+0.5);
        g.SetPointEXhigh(i,  myNPs[i].fPostFitUp);
        g.SetPointEXlow (i, -myNPs[i].fPostFitDown);
        if((myNPs[i].fFitValue + myNPs[i].fPostFitUp) > xmax) xmax = myNPs[i].fFitValue + myNPs[i].fPostFitUp;
        if((myNPs[i].fFitValue + myNPs[i].fPostFitDown) < xmin) xmin = myNPs[i].fFitValue + myNPs[i].fPostFitDown;
    }

    max   = myNPs.size();
    xmax *= (1.2-(xmax<0));
    xmin *= (0.8+(xmin<0));

    int lineHeight = 20;
    int offsetUp = 60;
    int offsetDown = 40;
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    gPad->SetLeftMargin(0.05/(8./6.));
    gPad->SetRightMargin(0.5);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy( "h_dummy_gamma","h_dummy_gamma",10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);

    TLine l0;
    TBox b1, b2;
    l0 = TLine(1,0,1,max);
    l0.SetLineStyle(7);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    g.Draw("psame");

    TLatex systs{};
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    if (fPlotLabel != "none") {
        TRExLabelNew(0.04, 1. - 40./newHeight, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str(), kBlack, gStyle->GetTextSize(), 0.09);
    }

    for(std::size_t i=0; i< myNPs.size(); ++i) {
        systs.DrawLatex(xmax*1.05, i+0.25, myNPs[i].fTitle.c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize(h_dummy.GetXaxis()->GetLabelSize()*0.9);
    gPad->RedrawAxis();

    Common::SaveCanvasAs(c, path);
}

//__________________________________________________________________________________
//
void FitResults::DrawNPPulls(const std::string &path,
                             const std::string &category,
                             const std::vector<std::shared_ptr<NormFactor> > &normFactors,
                             const std::vector<std::shared_ptr<ShapeFactor> > &shapeFactors,
                             const std::vector<std::string>& blinded,
                             const PullSigType& pullSigType) const {
    double xmin = -2.9;

    if (pullSigType != PullSigType::NoPullSig) xmin = 0;

    double xmax = 2.9;
    double max = 0.;
    static const std::vector<std::string> npToExclude = {"gamma_", "stat_", "shape_"};

    // reorder the NPs
    std::vector<std::shared_ptr<NuisParameter> > nuisPar;
    if(fNuisParList.size()>0){
        for (const auto& npName : fNuisParList){
            for (auto& np : fNuisPar){
                if (np.second->fName == npName) {
                    nuisPar.emplace_back(np.second);
                }

            }
        }
    } else {
        for (const auto& np : fNuisPar){
            nuisPar.emplace_back(np.second);
        }
    }

    TGraphAsymmErrors g{};

    NuisParameter *par = nullptr;
    int idx = 0;
    std::vector<std::pair<double, std::string> > pullQuantityNames;

    for(unsigned int i = 0; i<nuisPar.size(); ++i){
        par = nuisPar.at(i).get();

        std::string name = par->fName;
        name = Common::ReplaceString(name,"alpha_","");

        if (std::find(blinded.begin(), blinded.end(), name) != blinded.end()) continue;

        if (category != "all" && category != par->fCategory) continue;
        if (Common::FindInStringVector(fNuisParToHide,par->fName) >= 0) continue;

        bool skip = false;
        for (const std::string& ii : npToExclude) {
            if (par->fName.find(ii) != std::string::npos) {
                skip = true;
                break;
            }
        }
        for (const auto& norm : normFactors) {
            if (norm->fName == par->fName) {
                skip = true;
                break;
            }
        }
        for (const auto& shape : shapeFactors) {
            if (par->fName.find(shape->fName + "_bin_") != std::string::npos) {
                skip = true;
                break;
            }
        }
        if(skip) continue;

        double sig(0);
        if (pullSigType != PullSigType::NoPullSig) {
            double postFitErrAvg = (std::abs(par->fPostFitUp) + std::abs(par->fPostFitDown))/2;
            if (postFitErrAvg == 1. && par->fFitValue == 0) {
                postFitErrAvg += 0.001;
            } else if (postFitErrAvg >= 1.) {
                LOG(DEBUG) << "Parameter " << par->fTitle << " has a post-fit error that is larger than 1 (under-constrained). This shouldn't happen. Setting it to 0.999 in pull significacne calculation.\n";
                postFitErrAvg = 0.999;
            }
            sig = std::abs(par->fFitValue)/std::sqrt(1- postFitErrAvg*postFitErrAvg);
        }

        if (pullSigType == PullSigType::Regular) {
            g.SetPoint(idx, 0.,idx+0.5);
            g.SetPointEXhigh(idx, sig);
            g.SetPointEXlow( idx, 0.);
            g.SetPointEYhigh(idx, 0.4);
            g.SetPointEYlow( idx, 0.4);

            pullQuantityNames.push_back(std::make_pair(sig,par->fTitle));
        } else if (pullSigType == PullSigType::Rank) {
            pullQuantityNames.push_back(std::make_pair(sig,par->fTitle));
        } else {
            g.SetPoint(idx,par->fFitValue,idx+0.5);
            g.SetPointEXhigh(idx, par->fPostFitUp);
            g.SetPointEXlow( idx,-par->fPostFitDown);
            pullQuantityNames.push_back(std::make_pair(par->fFitValue,par->fTitle));
        }

        idx++;
        if(idx > max) max = idx;
    }

    if (pullSigType == PullSigType::Rank && idx != 0) {
        std::sort(pullQuantityNames.begin(), pullQuantityNames.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

        if (pullQuantityNames.at(0).first > xmax){ xmax = pullQuantityNames.at(0).first + 0.1; }
        if (max > 20){
            max = 20;
            pullQuantityNames.resize(20);
        }

        std::reverse(pullQuantityNames.begin(), pullQuantityNames.end());

        for(unsigned int i = 0; i < max; ++i) {
            g.SetPoint(i, 0.,i+0.5);
            g.SetPointEXhigh(i, pullQuantityNames.at(i).first);
            g.SetPointEXlow( i, 0.);
            g.SetPointEYhigh(i, 0.4);
            g.SetPointEYlow( i, 0.4);
        }
    }

    int lineHeight = 20;
    int offsetUp = 60;
    int offsetDown = 60;
    if (max < 10){
        offsetDown = 65;
    }
    int offset = offsetUp + offsetDown;
    int newHeight = offset + max*lineHeight;
    TCanvas c("c","c",800,newHeight);
    c.SetTicks(1,0);
    if (pullSigType == PullSigType::Rank) {
        gPad->SetRightMargin(0.05/(8./6.));
        gPad->SetLeftMargin(0.4);
    } else {
        gPad->SetRightMargin(0.5);
        gPad->SetLeftMargin(0.05/(8./6.));
    }
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    if (pullSigType != PullSigType::NoPullSig) gPad->SetGrid(true);

    TH1D h_dummy( ("h_dummy"+category).c_str(),("h_dummy"+category).c_str(),10,xmin,xmax);
    h_dummy.SetMaximum(max);
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.SetMinimum(0.);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);
    std::unique_ptr<TBox> b1(nullptr);
    std::unique_ptr<TBox> b2(nullptr);
    std::unique_ptr<TLine> l0(nullptr);

    if (pullSigType != PullSigType::NoPullSig) {
        for(int i_bin = 0; i_bin < h_dummy.GetNbinsX()+1; i_bin++) {
            h_dummy.SetBinContent(i_bin,-10);
        }

        g.SetFillColor(kAzure-4);
        g.SetLineColor(g.GetFillColor());
        g.Draw("2 gon same");
    } else {
        l0 = std::make_unique<TLine>(0,0,0,max);
        l0->SetLineStyle(7);
        l0->SetLineColor(kBlack);
        b1 = std::make_unique<TBox>(-1,0,1,max);
        b2 = std::make_unique<TBox>(-2,0,2,max);
        b1->SetFillColor(kGreen);
        b2->SetFillColor(kYellow);
        b2->Draw("same");
        b1->Draw("same");
        l0->Draw("same");

        g.Draw("psame");
    }

    TLatex systs{};
    double yshift = 0.25;
    double xpos = 3;
    if (pullSigType == PullSigType::Rank) {
        systs.SetTextAlign(32);
        yshift = 0.5;
        xpos = xmin-0.1;
    }
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs.DrawLatex(xpos,i+yshift, pullQuantityNames.at(i).second.c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize(h_dummy.GetXaxis()->GetLabelSize()*0.9);
    h_dummy.GetXaxis()->CenterTitle();
    if (pullSigType != PullSigType::NoPullSig) {
        h_dummy.GetXaxis()->SetTitle("Pull significance");
    } else {
        h_dummy.GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    }

    if (max < 10){
        h_dummy.GetXaxis()->SetTitleOffset(0.9);
    } else {
        h_dummy.GetXaxis()->SetTitleOffset(1.15);
    }

    if (fPlotLabel != "none") {
        double label_xpos = 0.04;
        if (pullSigType == PullSigType::Rank) {
            label_xpos = 0.4;
        }
        TRExLabelNew(label_xpos, 1. - 40./newHeight, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str(), kBlack, gStyle->GetTextSize(), 0.09);
    }

    gPad->RedrawAxis();

    if(category!="all"){
        TLatex cat_legend{};
        if (pullSigType == PullSigType::Rank) {
            cat_legend.DrawLatexNDC(0.7,1-0.8*offsetUp/newHeight,category.c_str());
        } else if (pullSigType == PullSigType::Regular) {
            cat_legend.DrawLatexNDC(0.3,1-0.8*offsetUp/newHeight,category.c_str());
        } else {
            cat_legend.DrawLatexNDC(0.5,1-0.8*offsetUp/newHeight,category.c_str());
        }
    }

    Common::SaveCanvasAs(c, path);
}


//__________________________________________________________________________________
//
void FitResults::DrawCorrelationMatrix(const std::string& path, const bool useHEPDataFormat, const double corrMin, bool EFTonly, const std::vector<std::string>& POIs) {
    if(fCorrMatrix){
        fCorrMatrix->fOutFolder = fOutFolder;
        fCorrMatrix->fNuisParToHide = fNuisParToHide;
        fCorrMatrix->fNuisParList = fNuisParList;
	    if (EFTonly) fCorrMatrix->fEFTParList=fEFTNFs;

        fCorrMatrix->SetPlotLabel(fPlotLabel);
	    fCorrMatrix->Draw(path, useHEPDataFormat, corrMin, POIs);
    }
}

//__________________________________________________________________________________
//
NuisParameter* FitResults::GetNuisanceParameter(const std::string& param) const {
    auto itr = fNuisPar.find(param);
    if (itr == fNuisPar.end()) {
        LOG(ERROR) << "Cannot find nuisance parameter " << param << "\n";
        exit(EXIT_FAILURE);
    }

    return itr->second.get();
}
