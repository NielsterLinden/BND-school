// Class include
#include "TRExFitter/TRExPlot.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/TRExFit.h"

// Style stuff
#include "StyleUtils/TRExStyle.h"
#include "StyleUtils/TRExLabels.h"
#include "StyleUtils/TRExUtils.h"

// ROOT includes
#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGaxis.h"

// c++ includes
#include <algorithm>
#include <iostream>

using namespace std;

//_____________________________________________________________________________
//
TRExPlot::TRExPlot(const std::string& name,int canvasWidth,int canvasHeight,bool hideRatioPad) :
    fName(name),
    h_data(nullptr),
    g_data(nullptr),
    h_stack(new THStack("","")),
    h_tot(nullptr),
    g_tot(nullptr),
    h_tot_bkg_prefit(nullptr),
    h_dummy(nullptr),
    c(std::make_unique<TCanvas>("","",canvasWidth,canvasHeight)),
    leg(nullptr),
    leg1(nullptr),
    pad0(nullptr),
    pad1(nullptr),
    xtitle("Variable [GeV]"),
    ytitle("Events"),
    fDataName("Data"),
    fLumi("XXX fb^{-1}"),
    fCME("13 TeV"),
    fPlotLabel("none"),
    yMaxScale(2.),
    NDF(-1),
    Chi2val(-1),
    fChi2prob(-1),
    KSprob(-1),
    fYmax(0),
    fYmin(0),
    fRatioYmax(2.),
    fRatioYmin(0.),
    fBinWidth(-1),
    fIsNjet(false),
    fShowYields(false),
    fLumiScale(1.),
    fLegendNColumns(2),
    fRatioYtitle(""),
    fRatioType(TRExPlot::RATIOTYPE::DATAOVERMC),
    fLabelX(-1),
    fLabelY(-1),
    fLegendX1(-1),
    fLegendX2(-1),
    fLegendY(-1),
    fPreFitLabel("Pre-fit"),
    fBlinding(nullptr) {

    if(hideRatioPad){
        pad0 = new TPad("","",0,0,1,1,0,0,0);
    }
    else{
        pad0 = new TPad("","",0,0.20,1,1,0,0,0);
    }
    pad0->SetTicks(1,1);
    pad0->SetTopMargin(0.05*(700./canvasHeight));
    if(hideRatioPad){
        pad0->SetBottomMargin(0.14*(600./canvasHeight));
    }
    else{
        pad0->SetBottomMargin(0.1);
    }
    pad0->SetLeftMargin(0.14*(600./canvasWidth));
    pad0->SetRightMargin(0.05*(600./canvasWidth));
    pad0->SetFrameBorderMode(0);
    pad0->SetFillStyle(0);
    //
    if(hideRatioPad){
        pad1 = nullptr;
    }
    else{
        pad1 = new TPad("","",0,0,1,0.28,0,0,0);
        pad1->SetTicks(1,1);
        pad1->SetTopMargin(0.0);
        pad1->SetBottomMargin(0.37*(700./canvasHeight));
        pad1->SetLeftMargin(0.14*(600./canvasWidth));
        pad1->SetRightMargin(0.05*(600./canvasWidth));
        pad1->SetFrameBorderMode(0);
        pad1->SetFillStyle(0);
    }
    if(pad1!=nullptr) pad1->Draw();
    pad0->Draw();
    pad0->cd();

    fTickMarksFunction = nullptr;
}

//_____________________________________________________________________________
//
TRExPlot::~TRExPlot(){
    delete h_stack;
    delete h_tot_bkg_prefit;
    delete h_dummy;
    delete leg;
    delete leg1;
}

//_____________________________________________________________________________
//
void TRExPlot::AddLabel(const std::string& name){
    fLabels.push_back(name);
}

//_____________________________________________________________________________
//
void TRExPlot::SetLumi(const std::string& name){
    fLumi = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetLumiScale(double scale){
    fLumiScale = scale;
}

//_____________________________________________________________________________
//
void TRExPlot::SetCME(const std::string& name){
    fCME = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetXaxis(const std::string& name,bool isNjet){
    xtitle = name;
    fIsNjet = isNjet;
}

//_____________________________________________________________________________
//
void TRExPlot::SetYaxis(const std::string& name){
    ytitle = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetYmaxScale(double scale){
    yMaxScale = scale;
}

//_____________________________________________________________________________
//
void TRExPlot::ResizeBinLabel(const int n) {
    fBinLabel.resize(n);
}

//_____________________________________________________________________________
//
void TRExPlot::SetBinLabel(int bin, const std::string& name){
    fBinLabel[bin] = name;
}

//_____________________________________________________________________________
//
void TRExPlot::SetBinWidth(double width){
    fBinWidth = width;
}

//_____________________________________________________________________________
//
void TRExPlot::SetData(TH1* h,std::string name){
    h_data.reset(static_cast<TH1*>(h->Clone()));
    h_data->SetDirectory(nullptr);
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    fDataName = name;
}

//_____________________________________________________________________________
//
void TRExPlot::AddSignal(TH1* h,std::string name){
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    unsigned int idx = std::find(fSigNames.begin(),fSigNames.end(),name) - fSigNames.begin();
    if(idx<fSigNames.size()){
        h_signal[idx]->Add(h,fLumiScale);
    }
    else{
        h_signal.emplace_back(static_cast<TH1*>(h->Clone()));
        h_signal.back()->SetDirectory(nullptr);
        h_signal.back()->Scale(fLumiScale);
        fSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::AddNormSignal(TH1* h,std::string name){
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    unsigned int idx = std::find(fNormSigNames.begin(),fNormSigNames.end(),name) - fNormSigNames.begin();
    if(idx<fNormSigNames.size()){
        h_normsig[idx]->Add(h,fLumiScale);
    }
    else{
        h_normsig.emplace_back(static_cast<TH1*>(h->Clone()));
        h_normsig.back()->SetDirectory(nullptr);
        h_normsig.back()->Scale(fLumiScale);
        fNormSigNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::AddOverSignal(TH1* h, std::string name, double custom_scale){
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    unsigned int idx = std::find(fOverSigNames.begin(),fOverSigNames.end(),name) - fOverSigNames.begin();
    if(idx<fOverSigNames.size()){
        h_oversig[idx]->Add(h,fLumiScale*custom_scale);
    }
    else{
        h_oversig.emplace_back(static_cast<TH1*>(h->Clone()));
        h_oversig.back()->SetDirectory(nullptr);
        h_oversig.back()->Scale(fLumiScale*custom_scale);
        fOverSigNames.push_back(name);
    }
    if (custom_scale != 1.) {
        LOG(INFO) << "Scaling OVERSIG histogram: " << name << " by " << custom_scale << "\n";
    }
}

//_____________________________________________________________________________
//
void TRExPlot::AddBackground(TH1* h,std::string name){
    if(h_tot==nullptr) {
        h_tot.reset(static_cast<TH1*>(h->Clone()));
        h_tot->SetDirectory(nullptr);
    }
    else h_tot->Add(h);
    // if no name is given, take the histogram title
    if(name=="") name = h->GetTitle();
    //
    unsigned int idx = std::find(fBkgNames.begin(),fBkgNames.end(),name) - fBkgNames.begin();
    if(idx<fBkgNames.size()){
        h_bkg[idx]->Add(h,fLumiScale);
    }
    else{
        h_bkg.emplace_back(static_cast<TH1*>(h->Clone()));
        h_bkg.back()->SetDirectory(nullptr);
        h_bkg.back()->Scale(fLumiScale);
        fBkgNames.push_back(name);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SetTot(TH1* h){
    h_tot.reset(static_cast<TH1*>(h->Clone()));
    h_tot->Scale(fLumiScale);
    h_tot->SetDirectory(nullptr);
    g_tot = std::make_unique<TGraphAsymmErrors>(h);
    for(int i=0;i<g_tot->GetN();i++){
        g_tot->GetY()[i]      *= fLumiScale;
        g_tot->GetEYlow()[i]  *= fLumiScale;
        g_tot->GetEYhigh()[i] *= fLumiScale;
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SetTotAsym(TGraphAsymmErrors* g){
    g_tot.reset(static_cast<TGraphAsymmErrors*>(g->Clone()));
    for(int i=0;i<g_tot->GetN();i++){
        g_tot->GetY()[i] *= fLumiScale;
        g_tot->GetEYlow()[i]  *= fLumiScale;
        g_tot->GetEYhigh()[i] *= fLumiScale;
    }
    for(int i=1;i<h_tot->GetNbinsX()+1;i++){
        h_tot->SetBinContent(i,g_tot->GetY()[i-1]);
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SetChi2KS(double chi2prob, double ksprob, double chi2val, int ndf) {
    fChi2prob = chi2prob;
    KSprob = ksprob;
    Chi2val = chi2val;
    NDF = ndf;
}

//_____________________________________________________________________________
//
void TRExPlot::BlindData(){
    //
    // Eventually blind bins

    if(h_data) {
        fBlinding = Common::BlindDataHisto(h_data.get(), fBlindedBins);
    } else {
        LOG(WARNING) << "Data histogram not defined. Blinding not possible. Skipped.\n";
    }
}

//_____________________________________________________________________________
//
TH1* TRExPlot::GetTotBkg() const{
    TH1* h = static_cast<TH1*>(h_tot->Clone("h_tot_bkg"));
    if(TRExFitter::ADDSTACKSIG){
        for (unsigned int i=0; i<fSigNames.size(); ++i) {
            h->Add( h_signal[i].get(), -1);
        }
    }
    h->SetDirectory(nullptr);
    return h;
}

//_____________________________________________________________________________
//
void TRExPlot::Draw(const std::string& options){

    fArrows.clear();

    /////////////////////////
    //
    // Main function of the class
    // ==========================
    //   It takes the data, background, signal to perform the full comparison (stack, ratio plot, ...)
    //
    /////////////////////////

    //
    // Draws an empty histogram to reserve the upper pad and set style
    //
    gStyle->SetEndErrorSize(0);
    pad0->cd();
    h_dummy = static_cast<TH1*>(h_tot->Clone("h_dummy"));
    h_dummy->SetDirectory(nullptr);
    h_dummy->Scale(0);
    if(pad0->GetWw() > pad0->GetWh()){
        h_dummy->GetYaxis()->SetTickLength(0.01);
        if(gStyle->GetTickLength("x")==0) h_dummy->SetNdivisions(0);
        else h_dummy->GetXaxis()->SetTickLength(0.02);
    }
    if (fXaxisRange.size() > 1){
        h_dummy->GetXaxis()->SetRangeUser(fXaxisRange.at(0), fXaxisRange.at(1));
    }
    h_dummy->Draw("HIST");
    if(options.find("log")!=std::string::npos) pad0->SetLogy();
    if(options.find("LOGX")!=std::string::npos) {
        pad0->SetLogx();
        pad1->SetLogx();
    }
    if(TRExFitter::OPTION["LogXSignalRegionPlot"]){
        pad0->SetLogx();
        pad1->SetLogx();
    }

    if(g_tot==nullptr) g_tot = std::make_unique<TGraphAsymmErrors>(h_tot.get());

    //
    // Determines if the data is real (and computes the poisson uncertainty) or not
    //
    bool hasData = true;
    if(h_data){
        h_data->SetMarkerSize(1.4);
        h_data->SetLineWidth(2);
        // build asym data
        if (TRExFitter::KEEPDATAERRORS) {
            g_data = histToGraph(h_data.get());
        } else {
            g_data = poissonize(h_data.get());
        }
    }
    else{
        hasData = false;
        h_data.reset(static_cast<TH1D*>(h_tot->Clone("dummyData")));//tajes data = total
        h_data->SetTitle("Asimov Data");
        h_data->SetDirectory(nullptr);
        g_data = std::make_unique<TGraphAsymmErrors>(h_data.get());
    }

    //
    // Add Bkg's to the stack
    //
    for(int i_smp=fBkgNames.size()-1;i_smp>=0;i_smp--){
        h_bkg[i_smp]->SetLineWidth(1);
        h_stack->Add(h_bkg[i_smp].get());
    }

    //
    // Eventually add Signal(s)
    //
    for(int i_smp=fSigNames.size()-1;i_smp>=0;i_smp--){
        h_signal[i_smp]->SetLineWidth(1);
        h_stack->Add(h_signal[i_smp].get());
    }

    // protection against bad settings for BinTickMarks
    bool doTickMarks = (fBinTickMarks.size() !=0);
    if (doTickMarks && (fBinTickMarks.back() > h_tot->GetNbinsX()+1 || fBinTickMarks.front() < 1)) {
        LOG(WARNING) << "Inconsistent `BinTickMarks`: first value cannot be lower than 1, last value cannot be larger than number of bins +1. Ignoring.\n";
        doTickMarks = false;
    }

    //
    // Draw
    //
    h_stack->Draw("HIST same");

    if( TRExFitter::PREFITONPOSTFIT && h_tot_bkg_prefit ) {
      h_tot_bkg_prefit->SetFillColor(0);
      h_tot_bkg_prefit->SetLineStyle(kDashed);
      h_tot_bkg_prefit->SetLineColor(kBlue);
      h_tot_bkg_prefit->SetLineWidth(2);
      h_tot_bkg_prefit->Draw("HIST same");
      if(doTickMarks) h_tot_bkg_prefit->SetNdivisions(0);
    }

    //
    // Total error bands style setting
    //
    g_tot->SetFillStyle(3354);
    if(TRExFitter::OPTION["SystFillStyle"]!=0) g_tot->SetFillStyle(TRExFitter::OPTION["SystFillStyle"]);
    g_tot->SetFillColor(kBlue-7);
    g_tot->SetLineColor(kWhite);
    g_tot->SetLineWidth(0);
    g_tot->SetMarkerSize(0);
    g_tot->Draw("sameE2");

    //
    // Draw a normalized signal distribution
    //
    std::vector<double> signalScale;
    signalScale.resize(fNormSigNames.size());
    for(int i_smp=fNormSigNames.size()-1;i_smp>=0;i_smp--){
        if (std::abs(h_normsig[i_smp]->Integral()) < 1e-10) {
            // division by zero
            LOG(WARNING) << "--- Signal " << fNormSigNames[i_smp] << " has integral equal to zero - cannot scale, returning scale = 1\n";
            signalScale[i_smp] = 1;
        } else {
            signalScale[i_smp] = h_tot->Integral()/h_normsig[i_smp]->Integral();
            LOG(INFO) << "--- Signal " << fNormSigNames[i_smp] << " scaled by " << signalScale[i_smp] << "\n";
        }
        h_normsig[i_smp]->Scale(signalScale[i_smp]);
        h_normsig[i_smp]->SetLineColor(h_normsig[i_smp]->GetFillColor());
        h_normsig[i_smp]->SetFillColor(0);
        h_normsig[i_smp]->SetFillStyle(0);
        h_normsig[i_smp]->SetLineStyle(2);
        h_normsig[i_smp]->SetLineWidth(2);
        h_normsig[i_smp]->Draw("HISTsame");
        if(doTickMarks) h_normsig[i_smp]->SetNdivisions(0);
    }

    //
    // Draw a overlayed signal distribution
    //
    for(int i_smp=fOverSigNames.size()-1;i_smp>=0;i_smp--){
        h_oversig[i_smp]->SetLineColor(h_oversig[i_smp]->GetFillColor());
        h_oversig[i_smp]->SetFillColor(0);
        h_oversig[i_smp]->SetFillStyle(0);
        h_oversig[i_smp]->SetLineStyle(2);
        h_oversig[i_smp]->SetLineWidth(2);
        h_oversig[i_smp]->Draw("HISTsame");
        if(doTickMarks) h_oversig[i_smp]->SetNdivisions(0);
    }

    //
    // Draw data (if it is real data of course)
    //
    if(hasData) g_data->Draw("Ep1 same");

    //
    // Draw blinding markers
    //
    if(fBlinding) {
        h_blind.reset(static_cast<TH1D*>(fBlinding->Clone()));
        h_blind->SetDirectory(nullptr);
        h_blind->SetLineWidth(0);
        h_blind->SetLineColor(kGray);
        h_blind->SetFillColor(kGray);
        h_blind->SetFillStyle(3345);
        h_blind->Draw("same HIST");
        if(doTickMarks) h_blind->SetNdivisions(0);
    }

    //
    // Axes labelling and style
    //
    h_dummy->GetXaxis()->SetTitle(xtitle.c_str());
    h_dummy->GetYaxis()->SetTitle(ytitle.c_str());
    if(fIsNjet){
        for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
            int nj = (int)h_dummy->GetXaxis()->GetBinCenter(i_bin);
            if(i_bin<h_dummy->GetNbinsX()) h_dummy->GetXaxis()->SetBinLabel( i_bin,Form("%d",nj) );
            else                           h_dummy->GetXaxis()->SetBinLabel( i_bin,Form("#geq%d",nj) );
        }
    }
    else if ((int)fBinLabel.size() > h_dummy->GetNbinsX()) {
        for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
            if(fBinLabel[i_bin]!="") h_dummy->GetXaxis()->SetBinLabel( i_bin, fBinLabel[i_bin].c_str());
        }
    }
    if(fBinLabel.size() > 1 && fBinLabel[1]!="") h_dummy->GetXaxis()->LabelsOption("d");
    double offset = 1.9*pad0->GetWh()/pad0->GetWw();
    h_dummy->GetYaxis()->SetTitleOffset( offset );

    double textHeight = 0.05*(672./pad0->GetWh());
    if(pad1==nullptr) textHeight *= 0.8;

    //
    // labels
    //
    double labelX = 0.18*(600./pad0->GetWw());
    if(fLabelX>=0){
        labelX = fLabelX;
    }
    // was 0.84-textHeight+0.04
    double labelY = 1-0.08*(700./c->GetWh());
    if(fLabelY>=0){
        labelY = fLabelY;
    }
    // scale it down to give space to text (in this way one can set the same value for LabelY and LegendY and have them vertically alligned)
    labelY -= textHeight - 0.015;

    if(fPlotLabel!="none") TRExLabel(labelX,labelY, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str());
    myText(labelX,labelY-textHeight,1,Form("#sqrt{s} = %s, %s",fCME.c_str(),fLumi.c_str()));//,0.045);
    for(unsigned int i_lab=0;i_lab<fLabels.size();i_lab++){
        myText(labelX,labelY-(i_lab+2)*textHeight,1,Form("%s",fLabels[i_lab].c_str()));//,0.045);
    }

    double legX1 = 1-fLegendNColumns*0.2*(600./pad0->GetWw())-0.1*(600./pad0->GetWw());
    if(fLegendX1>=0){
        legX1 = fLegendX1;
    }
    double legX2 = 1-0.1*(600./pad0->GetWw());
    if(fLegendX2>=0){
        legX2 = fLegendX2;
    }
    double legXmid = legX1+0.5*(legX2-legX1);

    double legY = 1-0.08*(700./c->GetWh());
    if(fLegendY>=0){
        legY = fLegendY;
    }

    if(fShowYields){
        legXmid = legX1+0.6*(legX2-legX1);
        leg  = new TLegend(legX1,legY-(fBkgNames.size()+fSigNames.size()+2+(hasData))*textHeight, legXmid,legY);
        leg1 = new TLegend(legXmid,leg->GetY1(), legX2,leg->GetY2());
        //
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        if(!TRExFitter::LEGENDLEFT) leg->SetTextAlign(32);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetMargin(0.22);
        leg1->SetFillStyle(0);
        leg1->SetBorderSize(0);
        leg1->SetTextAlign(32);
        leg1->SetTextFont(gStyle->GetTextFont());
        leg1->SetTextSize(gStyle->GetTextSize());
        leg1->SetMargin(0.);

        // Add data in the legend if real data are here
        if(hasData){
            leg->AddEntry(h_data.get(),fDataName.c_str(),"lep");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_data->Integral()),"");
        }

        // Signal and background legends
        for(unsigned int i_smp=0;i_smp<fSigNames.size();i_smp++){
            leg->AddEntry(h_signal[i_smp].get(), fSigNames[i_smp].c_str(),"f");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_signal[i_smp]->Integral()),"");
        }
        for(unsigned int i_smp=0;i_smp<fBkgNames.size();i_smp++){
            leg->AddEntry(h_bkg[i_smp].get(), fBkgNames[i_smp].c_str(),"f");
            leg1->AddEntry((TObject*)0,Form("%.1f",h_bkg[i_smp]->Integral()),"");
        }
        leg->AddEntry((TObject*)0,"Total","");
        leg1->AddEntry((TObject*)0,Form("%.1f",h_tot->Integral()),"");
        leg->AddEntry(g_tot.get(),"Uncertainty","f");
        leg1->AddEntry((TObject*)0," ","");

        if(TRExFitter::PREFITONPOSTFIT && h_tot_bkg_prefit) {
            leg->AddEntry(h_tot_bkg_prefit, (fPreFitLabel + " Bkgd.").c_str(),"l");
            leg1->AddEntry((TObject*)0," ","");
        }

        leg->Draw();
        leg1->Draw();
    }
    else{
        int Nrows = fBkgNames.size()+fSigNames.size()+fNormSigNames.size()+fOverSigNames.size();
        if(hasData) Nrows ++;
        Nrows ++; // for "Uncertainty"
        double legHeight = ((Nrows+fLegendNColumns-1)/fLegendNColumns)*textHeight;
        leg  = new TLegend(legX1,legY-legHeight, legX2,legY);
        leg->SetNColumns(fLegendNColumns);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        if(TRExFitter::LEGENDRIGHT) leg->SetTextAlign(32);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetMargin(fLegendNColumns*((0.85*pad0->GetWh())/(1.*pad0->GetWw()))*textHeight/(legX2-legX1));
        //
        // Draws data in the legend only is real data
        if(hasData){
            if(TRExFitter::REMOVEXERRORS) leg->AddEntry(h_data.get(),fDataName.c_str(),"ep");
            else                          leg->AddEntry(h_data.get(),fDataName.c_str(),"lep");
        }
        //
        std::vector<std::string> legNames;
        // Signal and background legend
        for(std::size_t i_smp = 0; i_smp < fSigNames.size(); ++i_smp) {
            auto it = std::find(legNames.begin(), legNames.end(), fSigNames[i_smp]);
            if (it == legNames.end()) {
                leg->AddEntry(h_signal[i_smp].get(), fSigNames[i_smp].c_str(),"f");
                legNames.emplace_back(fSigNames[i_smp]);
            }
        }
        for(unsigned int i_smp=0;i_smp<fNormSigNames.size();i_smp++) leg->AddEntry(h_normsig[i_smp].get(), Form("%.1f x %s *",signalScale.at(i_smp), fNormSigNames[i_smp].c_str()),"l");
        for(unsigned int i_smp=0;i_smp<fOverSigNames.size();i_smp++) leg->AddEntry(h_oversig[i_smp].get(), fOverSigNames[i_smp].c_str(),"l");
        for(std::size_t i_smp = 0; i_smp < fBkgNames.size(); ++i_smp) {
            auto it = std::find(legNames.begin(), legNames.end(), fBkgNames[i_smp]);
            if (it == legNames.end()) {
                leg->AddEntry(h_bkg[i_smp].get(), fBkgNames[i_smp].c_str(),"f");
                legNames.emplace_back(fBkgNames[i_smp]);
            }
        }
        leg->AddEntry(g_tot.get(),"Uncertainty","f");
        //
        if(TRExFitter::PREFITONPOSTFIT && h_tot_bkg_prefit) leg->AddEntry(h_tot_bkg_prefit, (fPreFitLabel + " Bkgd.").c_str(),"l");
        //
        leg->Draw();
        //
        if(fNormSigNames.size()>0){
            const float corr = (fLegendNColumns == 1) ?  -0.12 : 0;
            myText(legX1 + corr,0.96,  1,"*: normalised to total prediction");
        }
    }

    //
    // Ratio pad: drawing dummy histogram
    //
    if(pad1!=nullptr){
        pad1->cd();
        pad1->GetFrame()->SetY1(2);
        h_dummy2.reset(static_cast<TH1*>(h_tot->Clone("h_dummy2")));
        h_dummy2->SetDirectory(nullptr);
        h_dummy2->Scale(0);
        if(pad0->GetWw() > pad0->GetWh()) h_dummy2->GetYaxis()->SetTickLength(0.01);
        if(gStyle->GetTickLength("x")==0) h_dummy2->SetNdivisions(0);
        h_dummy2->Draw("HIST");
        h_dummy2->GetYaxis()->SetTitleOffset(1.*h_dummy->GetYaxis()->GetTitleOffset());
        if (fXaxisRange.size() > 1){
            h_dummy2->GetXaxis()->SetRangeUser(fXaxisRange.at(0),fXaxisRange.at(1));
        }

        //
        // Initialising the ratios
        //    h_ratio: is the real Data/MC ratio
        //    h_ratio2: is a MC/MC ratio to plot the uncertainty band
        //
        if(fRatioType==TRExPlot::RATIOTYPE::SOVERB || fRatioType == TRExPlot::RATIOTYPE::SOVERSQRTB || fRatioType==TRExPlot::RATIOTYPE::SOVERSQRTSPLUSB){
            if(fSigNames.size()>0){
                h_ratio.reset(static_cast<TH1*>(h_signal[0] ->Clone()));
                h_ratio->SetDirectory(nullptr);
                h_ratio->SetLineColor(h_ratio->GetFillColor());
                h_ratio->SetFillStyle(0);
            }
            else if(fNormSigNames.size()>0){
                h_ratio.reset(static_cast<TH1*>(h_normsig[0]->Clone()));
                h_ratio->SetDirectory(nullptr);
                h_ratio->Scale(1./signalScale[0]);
            }
            else if(fOverSigNames.size()>0) {
                h_ratio.reset(static_cast<TH1*>(h_oversig[0]->Clone()));
                h_ratio->SetDirectory(nullptr);
            } else {
                h_ratio.reset(static_cast<TH1*>(h_tot->Clone()));
                h_ratio->SetDirectory(nullptr);
                h_ratio->Scale(0);
            }
        }
        else{
            h_ratio.reset(static_cast<TH1*>(h_data->Clone()));
            h_ratio->SetDirectory(nullptr);
        }

        // in case of S/B,.. and several signal samples, build other ratios
        std::vector<TH1*> h_addRatioVec;
        if(fRatioType==TRExPlot::RATIOTYPE::SOVERB || fRatioType == TRExPlot::RATIOTYPE::SOVERSQRTB || fRatioType==TRExPlot::RATIOTYPE::SOVERSQRTSPLUSB){
            if(fSigNames.size()>1){
                for(unsigned int i_sig=1;i_sig<fSigNames.size();i_sig++){
                    h_addRatioVec.push_back(static_cast<TH1*>(h_signal[i_sig] ->Clone()));
                    h_addRatioVec[i_sig-1]->SetLineColor(h_addRatioVec[i_sig-1]->GetFillColor());
                    h_addRatioVec[i_sig-1]->SetFillStyle(0);
                }
            }
            else if(fNormSigNames.size()>1){
                for(unsigned int i_sig=1;i_sig<fNormSigNames.size();i_sig++){
                    h_addRatioVec.push_back(static_cast<TH1*>(h_normsig[i_sig]->Clone(Form("h_ratio_%d",i_sig))));
                    h_addRatioVec[i_sig-1]->Scale(1./signalScale[i_sig]);
                }
            }
            else if(fOverSigNames.size()>1){
                for(unsigned int i_sig=1;i_sig<fOverSigNames.size();i_sig++){
                    h_addRatioVec.push_back(static_cast<TH1*>(h_oversig[i_sig]->Clone(Form("h_ratio_%d",i_sig))));
                }
            }
        }
        // do somehting similar also when plotting data/bkg and S+B/bkg
        else if(fRatioType==TRExPlot::RATIOTYPE::DATAOVERBWITHS){
            for(unsigned int i_sig=0;i_sig<fSigNames.size();i_sig++){
                h_addRatioVec.push_back(static_cast<TH1*>(h_signal[i_sig] ->Clone()));
                h_addRatioVec[i_sig]->SetLineColor(h_addRatioVec[i_sig]->GetFillColor());
                h_addRatioVec[i_sig]->SetFillStyle(0);
            }
            for(unsigned int i_sig=0;i_sig<fNormSigNames.size();i_sig++){
                h_addRatioVec.push_back(static_cast<TH1*>(h_normsig[i_sig]->Clone(Form("h_ratio_%d",i_sig))));
                h_addRatioVec[i_sig]->Scale(1./signalScale[i_sig]);
            }
            for(unsigned int i_sig=0;i_sig<fOverSigNames.size();i_sig++){
                h_addRatioVec.push_back(static_cast<TH1*>(h_oversig[i_sig]->Clone(Form("h_ratio_%d",i_sig))));
            }
        }

        h_tot_nosyst.reset(static_cast<TH1*>(h_tot->Clone("h_tot_nosyst")));
        h_tot_nosyst->SetDirectory(nullptr);
        for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){
            h_tot_nosyst->SetBinError(i_bin,0);
        }
        g_ratio2.reset(static_cast<TGraphAsymmErrors*>(g_tot->Clone("g_ratio2")));

        //
        // Plots style
        //
        std::string ratioTitle = "";
        if(fRatioYtitle == ""){
            if(fRatioType == TRExPlot::RATIOTYPE::DATAOVERMC)      ratioTitle = "Data / Pred.";
            if(fRatioType == TRExPlot::RATIOTYPE::DATAOVERB)       ratioTitle = "Data / Bkg.";
            if(fRatioType == TRExPlot::RATIOTYPE::SOVERB)          ratioTitle = "S / B";
            if(fRatioType == TRExPlot::RATIOTYPE::SOVERSQRTB)      ratioTitle = "S / #sqrt{B}";
            if(fRatioType == TRExPlot::RATIOTYPE::SOVERSQRTSPLUSB) ratioTitle = "S / #sqrt{S+B}";
            if(fRatioType == TRExPlot::RATIOTYPE::DATAOVERBWITHS)  ratioTitle = "Data / Bkg.";
        }
        else{
            ratioTitle = fRatioYtitle;
        }
        if(ratioTitle == "-") ratioTitle = "";
        h_dummy2->GetYaxis()->SetTitle(ratioTitle.c_str());
        h_dummy2->GetYaxis()->SetNdivisions(504,false);
        gStyle->SetEndErrorSize(0);


        //
        // Compute Data/MC ratio
        //
        if(fRatioType == TRExPlot::RATIOTYPE::SOVERSQRTB){
            for(int i_bin=1;i_bin<=h_ratio->GetNbinsX();i_bin++){
                double bkgYield(0);
                for (const auto& ibkg : h_bkg) {
                    bkgYield += ibkg->GetBinContent(i_bin);
                }
                h_ratio->SetBinContent(i_bin,h_ratio->GetBinContent(i_bin)/sqrt(bkgYield));
                for(auto h_tmp : h_addRatioVec) h_tmp->SetBinContent(i_bin,h_tmp->GetBinContent(i_bin)/sqrt(bkgYield));
            }
        }
        else if(fRatioType == TRExPlot::RATIOTYPE::SOVERSQRTSPLUSB){
            for(int i_bin=1;i_bin<=h_ratio->GetNbinsX();i_bin++){
                double bkgYield(0);
                for (const auto& ibkg : h_bkg) {
                    bkgYield += ibkg->GetBinContent(i_bin);
                }
                h_ratio->SetBinContent(i_bin,h_ratio->GetBinContent(i_bin)/sqrt(h_ratio->GetBinContent(i_bin)+bkgYield));
                for(auto h_tmp : h_addRatioVec) h_tmp->SetBinContent(i_bin,h_tmp->GetBinContent(i_bin)/sqrt(h_tmp->GetBinContent(i_bin)+bkgYield));
            }
        }
        else if (fRatioType==TRExPlot::RATIOTYPE::DATAOVERB) {
            std::unique_ptr<TH1> tmp(TRExPlot::GetTotBkg());
            for (int ibin = 1; ibin <= tmp->GetNbinsX(); ++ibin) {
                tmp->SetBinError(ibin, 0);
            }
            h_ratio->Divide(tmp.get());
        }
        else if (fRatioType==TRExPlot::RATIOTYPE::DATAOVERBWITHS) {
            std::unique_ptr<TH1> h_tot_bkg(TRExPlot::GetTotBkg());
            for (int ibin = 1; ibin <= h_tot_bkg->GetNbinsX(); ++ibin) {
                h_tot_bkg->SetBinError(ibin, 0);
            }
            h_ratio->Divide(h_tot_bkg.get());
            for (auto h_tmp : h_addRatioVec) {
                h_tmp->Add(h_tot_bkg.get());
                h_tmp->Divide(h_tot_bkg.get());
            }
        }
        else{
            h_ratio->Divide(h_tot_nosyst.get());
            for(auto h_tmp : h_addRatioVec) h_tmp->Divide(h_tot_nosyst.get());
        }
        if(TRExFitter::OPRATIO) h_ratio->SetMarkerStyle(24);
        else                    h_ratio->SetMarkerStyle(h_data->GetMarkerStyle());
        h_ratio->SetMarkerSize(1.4);
        h_ratio->SetMarkerColor(kBlack);
        h_ratio->SetLineWidth(2);
        g_ratio = histToGraph(h_ratio.get());
        std::unique_ptr<TH1> h_tot_bkg(TRExPlot::GetTotBkg());
        TH1* ref = (fRatioType==TRExPlot::RATIOTYPE::DATAOVERB || fRatioType==TRExPlot::RATIOTYPE::DATAOVERBWITHS) ? h_tot_bkg.get() : h_tot_nosyst.get();
        for(int i_bin=1;i_bin<=h_ratio->GetNbinsX();i_bin++){
            //For the ratio plot, the error is just to illustrate the "poisson uncertainty on the data"
            if(TRExFitter::REMOVEXERRORS){
                g_ratio->SetPointEXhigh(i_bin-1, 0.);
                g_ratio->SetPointEXlow( i_bin-1, 0.);
            }
            g_ratio->SetPointEYhigh(i_bin-1,g_data->GetErrorYhigh(i_bin-1)/ref->GetBinContent(i_bin));
            g_ratio->SetPointEYlow( i_bin-1,g_data->GetErrorYlow(i_bin-1) /ref->GetBinContent(i_bin));
        }

        //
        // Compute the MC/MC ratio (for uncertainty band in the bottom pad)
        //
        for(int i_bin=1;i_bin<h_tot_nosyst->GetNbinsX()+1;i_bin++){
            g_ratio2->SetPoint(i_bin-1,g_ratio2->GetX()[i_bin-1],g_ratio2->GetY()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
            g_ratio2->SetPointEXlow(i_bin-1,g_ratio2->GetEXlow()[i_bin-1]);
            g_ratio2->SetPointEXhigh(i_bin-1,g_ratio2->GetEXhigh()[i_bin-1]);
            if(h_tot_nosyst->GetBinContent(i_bin)>1e-4){
                g_ratio2->SetPointEYlow(i_bin-1,g_ratio2->GetEYlow()[i_bin-1]  /h_tot_nosyst->GetBinContent(i_bin));
                g_ratio2->SetPointEYhigh(i_bin-1,g_ratio2->GetEYhigh()[i_bin-1]/h_tot_nosyst->GetBinContent(i_bin));
            }
            else{
                g_ratio2->SetPointEYlow(i_bin-1,0.);
                g_ratio2->SetPointEYhigh(i_bin-1,0.);
            }
        }

        //
        // Now draws everything
        //
        if(fRatioType== TRExPlot::RATIOTYPE::DATAOVERMC || fRatioType== TRExPlot::RATIOTYPE::DATAOVERB || fRatioType== TRExPlot::RATIOTYPE::DATAOVERBWITHS){
            hline = std::make_unique<TLine>(h_dummy2->GetXaxis()->GetXmin(),1,h_dummy2->GetXaxis()->GetXmax(),1);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->SetLineStyle(2);
            hline->Draw();
        }
        //
        h_dummy2->SetMinimum(fRatioYmin);
        h_dummy2->SetMaximum(fRatioYmax);
        //
        h_dummy2->GetXaxis()->SetTitle(h_dummy->GetXaxis()->GetTitle());
        //
        h_dummy->GetXaxis()->SetTitle("");
        h_dummy->GetXaxis()->SetLabelSize(0);

        if(fRatioType== TRExPlot::RATIOTYPE::DATAOVERMC || fRatioType== TRExPlot::RATIOTYPE::DATAOVERB || fRatioType== TRExPlot::RATIOTYPE::DATAOVERBWITHS){
            g_ratio2->Draw("sameE2");
        }

        if(TRExFitter::PREFITONPOSTFIT && h_tot_bkg_prefit) {
            h_ratio_pre.reset(static_cast<TH1D*>(h_data->Clone()));
            h_ratio_pre->SetDirectory(nullptr);
            h_ratio_pre->Divide(h_tot_bkg_prefit);
            h_ratio_pre->SetFillColor(0);
            h_ratio_pre->SetLineStyle(kDashed);
            h_ratio_pre->SetLineColor(kBlue);
            h_ratio_pre->SetLineWidth(2);
            h_ratio_pre->Draw("HIST same");
        }

        bool customLabels = false;
        for(int i_bin=1;i_bin<h_dummy->GetNbinsX()+1;i_bin++){
            if(((std::string)h_dummy->GetXaxis()->GetBinLabel(i_bin))!=""){
                h_dummy2->GetXaxis()->SetBinLabel( i_bin, h_dummy->GetXaxis()->GetBinLabel(i_bin));
                customLabels = true;
            }
        }

        // if ratio style is S/B or similar, draw it and all the additional ratios
        if(fRatioType==TRExPlot::RATIOTYPE::SOVERB || fRatioType == TRExPlot::RATIOTYPE::SOVERSQRTB || fRatioType==TRExPlot::RATIOTYPE::SOVERSQRTSPLUSB){
            h_ratio->Draw("HIST same");
            if(doTickMarks) h_ratio->SetNdivisions(0);
            for(auto h_tmp : h_addRatioVec){
                h_tmp->Draw("HIST same");
                if(doTickMarks) h_tmp->SetNdivisions(0);
            }
        }
        // otherwise, draw only if there is data
        else if(hasData){
            g_ratio->Draw("pe0");
        }
        // regardless of data being there, in case of data/B + (S+B)/B, draw the additional ratios (showing (S+B)/B)
        if(fRatioType==TRExPlot::RATIOTYPE::DATAOVERBWITHS){
            for(auto h_tmp : h_addRatioVec){
                h_tmp->Draw("HIST same");
                if(doTickMarks) h_tmp->SetNdivisions(0);
            }
        }

        //
        // Mark blinded bins in ratio pad as  well
        //
        if(h_blind!=nullptr){
            h_blindratio.reset(static_cast<TH1D*>(h_blind->Clone()));
            h_blindratio->SetDirectory(nullptr);
            h_blindratio->Scale(2.);
            h_blindratio->Draw("HIST same");
            if(doTickMarks) h_blindratio->SetNdivisions(0);
        }

        if(fBinLabel.size() > 1 && fBinLabel[1]!="") h_dummy2->GetXaxis()->LabelsOption("d");
        h_dummy2->GetXaxis()->SetLabelOffset( h_dummy2->GetXaxis()->GetLabelOffset()+0.02 );
        if(customLabels && h_dummy->GetNbinsX()>10) h_dummy2->GetXaxis()->SetLabelSize(0.66*h_dummy2->GetXaxis()->GetLabelSize() );
        if(customLabels) h_dummy2->GetXaxis()->SetLabelOffset( h_dummy2->GetXaxis()->GetLabelOffset()+0.02 );
        if(fRangeLabels.size()) h_dummy2->GetXaxis()->SetLabelSize( h_dummy2->GetXaxis()->GetLabelSize()*0.9);

        h_dummy2->GetYaxis()->ChangeLabel(-1,-1,-1,-1,-1,-1," ");

        // Make custom tick marks for specific bins, if option BinTickMarks is used
        if(doTickMarks) CustomiseTickMarks(h_dummy2.get());

        gPad->RedrawAxis();

        //
        // Add arrows when the ratio is beyond the limits of the ratio plot
        //
        if(fRatioType== TRExPlot::RATIOTYPE::DATAOVERMC || fRatioType== TRExPlot::RATIOTYPE::DATAOVERB || fRatioType== TRExPlot::RATIOTYPE::DATAOVERBWITHS){
            for(int i_bin=0;i_bin<h_tot_nosyst->GetNbinsX()+2;i_bin++){

                if (i_bin==0 || i_bin>h_tot_nosyst->GetNbinsX()) continue; //skip under/overflow bins

                double val = h_ratio->GetBinContent(i_bin);

                double maxRange = h_dummy2->GetMaximum();
                double minRange = h_dummy2->GetMinimum();

                int isUp=0; //1==up, 0==nothing, -1==down
                if ( val<minRange ) isUp=-1;
                else if (val>maxRange ) isUp=1;
                if (val==0) isUp=0;

                if (isUp!=0) {
                    if (isUp==1) {
                        fArrows.emplace_back(std::make_unique<TArrow>(h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmax-0.05*(fRatioYmax-fRatioYmin), h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmax,0.030/(pad0->GetWw()/596.),"|>"));
                    } else {
                        fArrows.emplace_back(std::make_unique<TArrow>(h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmin+0.05*(fRatioYmax-fRatioYmin), h_ratio->GetXaxis()->GetBinCenter(i_bin),fRatioYmin,0.030/(pad0->GetWw()/596.),"|>"));
                    }
                    fArrows.back()->SetFillColor(10);
                    fArrows.back()->SetFillStyle(1001);
                    fArrows.back()->SetLineColor(kBlue-7);
                    fArrows.back()->SetLineWidth(2);
                    fArrows.back()->SetAngle(40);
                    fArrows.back()->Draw();
                }
            }
        }

        //
        // Draw bin dividers (ratio pad)
        //
        if(fBinDividers.size()>0) {
            std::vector<int> CheckedDividers;
            std::copy_if(fBinDividers.begin(), fBinDividers.end(), std::back_inserter(CheckedDividers),
                [&](int i_bin){
                    if(i_bin < 1 || i_bin > h_dummy->GetNbinsX()) {
                        LOG(WARNING) << "--- BinDivider at bin " << i_bin << " out of range\n";
                        return false;
                    }
                    return true;
                });
            fBinDividers = CheckedDividers;
            vline = std::make_unique<TLine>();
            vline->SetLineColor(kBlack);
            vline->SetLineStyle(2);
            vline->SetLineWidth(1);
            for(int i_bin : fBinDividers) {
                vline->DrawLine(h_dummy2->GetBinLowEdge(i_bin), fRatioYmin,
                                h_dummy2->GetBinLowEdge(i_bin), fRatioYmax);
            }
        }
        // ---

        //
        // Draw range labels
        //
        if(fRangeLabels.size()>0) {
            if(fRangeLabels.size()!=fBinDividers.size()+1) {
                LOG(WARNING) << "--- RangeLabels (length " << fRangeLabels.size() << ") incompatible with BinDividers (length " <<
                                    fBinDividers.size() << ". Will not draw RangeLabels\n";
            }
            else {
                auto rlab = std::make_unique<TLatex>();
                rlab->SetNDC(0);
                rlab->SetTextFont(42);
                rlab->SetTextSize(0.1);
                std::vector<int> RangeLabelLowBins = fBinDividers;
                RangeLabelLowBins.insert(RangeLabelLowBins.begin(), 1);
                for(size_t i=0; i<fRangeLabels.size(); i++) {
                    int i_bin = RangeLabelLowBins.at(i);
                    rlab->DrawLatex(h_dummy2->GetBinLowEdge(i_bin), fRatioYmin - 0.38*(fRatioYmax-fRatioYmin),
                                   fRangeLabels.at(i).c_str());
                }
            }
        }
        // ---

        pad1->cd();
        KSlab = std::make_unique<TLatex>();
        KSlab->SetNDC(1);
        KSlab->SetTextFont(42);
        KSlab->SetTextSize(0.1);
        std::string kslab = "";
        if(Chi2val >= 0)   kslab += Form("   #chi^{2}/ndf = %.1f",Chi2val);
        if(NDF >= 0)       kslab += Form(" / %d",NDF);
        if(fChi2prob >= 0) kslab += Form("  #chi^{2}prob = %.2f",fChi2prob);
        if(KSprob >= 0)    kslab += Form("  KS prob = %.2f",KSprob);
        KSlab->DrawLatex(0.15,0.9,kslab.c_str());
        //
        pad0->cd();
    }

    //
    // Set bin width and eventually divide larger bins by this bin width
    if(fBinWidth>0){
        for(unsigned int i_smp=0;i_smp<fSigNames.size();i_smp++)      SetHistBinWidth(h_signal[i_smp].get(), fBinWidth);
        for(unsigned int i_smp=0;i_smp<fNormSigNames.size();i_smp++)  SetHistBinWidth(h_normsig[i_smp].get(),fBinWidth);
        for(unsigned int i_smp=0;i_smp<fOverSigNames.size();i_smp++)  SetHistBinWidth(h_oversig[i_smp].get(),fBinWidth);
        for(unsigned int i_smp=0;i_smp<fBkgNames.size();i_smp++)      SetHistBinWidth(h_bkg[i_smp].get(),    fBinWidth);
        //
        if(h_tot) SetHistBinWidth(h_tot.get(),fBinWidth);
        if(g_tot) SetGraphBinWidth(g_tot.get(),fBinWidth);
        if(h_tot_bkg_prefit) SetHistBinWidth(h_tot_bkg_prefit, fBinWidth);
        if(h_data) SetHistBinWidth(h_data.get(),fBinWidth);
        if(g_data) SetGraphBinWidth(g_data.get(),fBinWidth);
        // try to guess y axis label...
        if(ytitle=="Events"){
            // extract unit from square brackets if present
            std::string unit = "";
            std::size_t bracketLeft = xtitle.find('[');
            std::size_t bracketRight = xtitle.find(']');
            if(bracketLeft != std::string::npos && bracketRight != std::string::npos && bracketRight > bracketLeft){
                unit = xtitle.substr(bracketLeft + 1, bracketRight - bracketLeft - 1);
            }

            if(!unit.empty()){
                if((int)fBinWidth==fBinWidth) ytitle = Form("Events / %.0f %s",fBinWidth,unit.c_str());
                else if((int)(fBinWidth*10)==(fBinWidth*10)) ytitle = Form("Events / %.1f %s",fBinWidth,unit.c_str());
                else if((int)(fBinWidth*100)==(fBinWidth*100)) ytitle = Form("Events / %.2f %s",fBinWidth,unit.c_str());
                else if((int)(fBinWidth*1000)==(fBinWidth*1000)) ytitle = Form("Events / %.3f %s",fBinWidth,unit.c_str());
                // ...
            }
            else{
                if((int)fBinWidth==fBinWidth) ytitle = Form("Events / %.0f",fBinWidth);
                else if((int)(fBinWidth*10)==(fBinWidth*10)) ytitle = Form("Events / %.1f",fBinWidth);
                else if((int)(fBinWidth*100)==(fBinWidth*100)) ytitle = Form("Events / %.2f",fBinWidth);
                else if((int)(fBinWidth*1000)==(fBinWidth*1000)) ytitle = Form("Events / %.3f",fBinWidth);
                // ...
            }
            h_dummy->GetYaxis()->SetTitle(ytitle.c_str());
        }
    }

    // turn off x-error bars
    if(TRExFitter::REMOVEXERRORS){
        for (int i=0; i < g_data->GetN(); i++) {
            g_data->SetPointEXlow(i,0);
            g_data->SetPointEXhigh(i,0);
        }
    }

    // Fix y max
    //
    double yMax = 0.;
    // take into account also total prediction uncertainty
    for(int i_bin=1;i_bin<h_tot->GetNbinsX()+1;i_bin++){
        double y = h_tot->GetBinContent(i_bin);
        if(y>yMax) yMax = y;
        if(hasData && h_data!=nullptr && g_data!=nullptr){
            if(h_data->Integral()>0 && h_data->GetBinContent(i_bin)>0 && g_data->GetY()[i_bin-1]>0 && g_data->GetEYhigh()[i_bin-1]>0){
                y = h_data->GetBinContent(i_bin)+g_data->GetEYhigh()[i_bin-1];
                if(y>yMax) yMax = y;
            }
        }
    }
    //
    if(options.find("log")==std::string::npos){
        if(fYmax!=0) h_dummy->SetMaximum(fYmax);
        else         h_dummy->SetMaximum(yMaxScale*yMax);
        if(fYmin>0)  h_dummy->SetMinimum(fYmin);
        else         h_dummy->SetMinimum(0.);
    }
    else{
        if(fYmax!=0) h_dummy->SetMaximum(fYmax);
        else         h_dummy->SetMaximum(yMax*pow(10,yMaxScale));
        if(fYmin>0)  h_dummy->SetMinimum(fYmin);
        else         h_dummy->SetMinimum(1.);
    }

    // Make custom tick marks for specific bins, if option BinTickMarks is used
    if(doTickMarks) CustomiseTickMarks(h_dummy);

    //
    // Fix / redraw axis
    //
    pad0->RedrawAxis();

    if(h_blind!=nullptr){
        h_blind->Scale(h_dummy->GetMaximum());
    }

    //
    // Draw bin dividers (top pad)
    //
    if(fBinDividers.size()>0 && vline) {
        for(int i_bin : fBinDividers) {
            if (i_bin >= 1 && i_bin <= h_dummy->GetNbinsX()) {
                vline->DrawLine(h_dummy->GetBinLowEdge(i_bin), h_dummy->GetMinimum(),
                                h_dummy->GetBinLowEdge(i_bin), yMax);
            }
        }
    }
}

//_____________________________________________________________________________
//
void TRExPlot::SaveAs(const std::string& name) const{
	Common::SaveCanvasAs(*c, name);
}

//_____________________________________________________________________________
//
void TRExPlot::SaveAsBkgOnly(const std::string& name) const{
    if (fBkgOnlyCanvas) Common::SaveCanvasAs(*fBkgOnlyCanvas, name);
}

//_____________________________________________________________________________
//
void TRExPlot::WriteToFile(const std::string& name) const{
    TDirectory *here = gDirectory;
    std::unique_ptr<TFile> f (TFile::Open(name.c_str(),"RECREATE"));
    if (!f) {
        LOG(WARNING) << "Cannot open file: " << name << ". Not writing into a file!\n";
        return;
    }
    f->cd();
    if(h_data) h_data->Write(Form("h_%s",fDataName.c_str()),TObject::kOverwrite);
    h_tot->Write("h_totErr",TObject::kOverwrite);
    if(g_tot) g_tot->Write("g_totErr",TObject::kOverwrite);
    for(int i_smp=fBkgNames.size()-1;i_smp>=0;i_smp--){
        h_bkg[i_smp]->Write(Form("h_%s",fBkgNames[i_smp].c_str()),TObject::kOverwrite);
    }
    for(int i_smp=fSigNames.size()-1;i_smp>=0;i_smp--){
        h_signal[i_smp]->Write(Form("h_%s",fSigNames[i_smp].c_str()),TObject::kOverwrite);
        if(h_normsig[i_smp]) h_normsig[i_smp]->Write(Form("h_%s_norm",fSigNames[i_smp].c_str()),TObject::kOverwrite);
    }
    here->cd();
    f->Close();
}

//_____________________________________________________________________________
//
void TRExPlot::SetBinBlinding(const std::vector<int>& bins){
    fBlindedBins = bins;
}

//_____________________________________________________________________________
//
void TRExPlot::DrawDataMinusBkg(const bool errOnSignal) {
    if (!h_data) return;
    fBkgOnlyCanvas = std::make_unique<TCanvas>("","",c->GetWw(),c->GetWh());
    fBkgOnlyCanvas->SetTicks(1,1);
    fBkgOnlyCanvas->SetBottomMargin(0.1);

    // get data - bkg
    std::unique_ptr<TH1> bkgStack(nullptr);
    fDataMinusBkg.reset(static_cast<TH1*>(h_data->Clone()));
    for (const auto& ibkg : h_bkg) {
        if (!bkgStack) {
            bkgStack.reset(static_cast<TH1*>(ibkg->Clone()));
        } else {
            bkgStack->Add(ibkg.get());
        }
    }
    for (const auto& isig : h_signal) {
        if (!fSigStack) {
            fSigStack.reset(static_cast<TH1*>(isig->Clone()));
        } else {
            fSigStack->Add(isig.get());
        }
    }

    fDataMinusBkg->Add(bkgStack.get(), -1);
    if(errOnSignal) {
        Common::SetGraphToHist(fBkgOnlyErr.get(), fSigStack.get());
    } else {
        Common::SetGraphToZero(fBkgOnlyErr.get());
    }

    fBkgOnlyErr->SetFillStyle(3354);
    if(TRExFitter::OPTION["SystFillStyle"]!=0) fBkgOnlyErr->SetFillStyle(TRExFitter::OPTION["SystFillStyle"]);
    fBkgOnlyErr->SetFillColor(kBlue-7);
    fBkgOnlyErr->SetLineColor(kWhite);
    fBkgOnlyErr->SetLineWidth(0);
    fBkgOnlyErr->SetMarkerSize(0);

    double min(9999.);
    double max(-9999.);
    for (int i = 0; i < fBkgOnlyErr->GetN(); ++i) {
        double x,y;
        fBkgOnlyErr->GetPoint(i, x, y);
        const double errDown = fBkgOnlyErr->GetErrorYlow(i);
        const double errUp = fBkgOnlyErr->GetErrorYhigh(i);
        const double totDown = y - errDown;
        const double totUp = y + errUp;
        if (totDown < min) min = totDown;
        if (totUp   > max) max = totUp;
    }
    max = std::max(fSigStack->GetMaximum(), max);
    min = std::min(fSigStack->GetMinimum(), min);

    double min_DataMinusBkg(9999.);
    double max_DataMinusBkg(-9999.);
    for (int i = 1; i <= fDataMinusBkg->GetNbinsX(); ++i) {
        double bincont = fDataMinusBkg->GetBinContent(i);
        const double err = fDataMinusBkg->GetBinError(i);
        const double totDown = bincont - err;
        const double totUp = bincont + err;
        if (totDown < min_DataMinusBkg) min_DataMinusBkg = totDown;
        if (totUp   > max_DataMinusBkg) max_DataMinusBkg = totUp;
    }

    min = std::min(min_DataMinusBkg, min);
    max = std::min(max_DataMinusBkg, max);

    min = min > 0 ? 0 : min*1.3;
    max = max > 0 ? max*1.6 : max*0.5;

    fSigStack->GetYaxis()->SetRangeUser(min, max);
    fSigStack->GetXaxis()->SetTitle(xtitle.c_str());
    fSigStack->GetYaxis()->SetTitle(ytitle.c_str());
    fSigStack->GetYaxis()->SetTitleOffset(fSigStack->GetYaxis()->GetTitleOffset()*1.5);

    fSigStack->Draw("hist");
    fBkgOnlyErr->Draw("E2 same");
    if (TRExFitter::REMOVEXERRORS) {
        fDataMinusBkg->Draw("P0EX0same");
    } else {
        fDataMinusBkg->Draw("P0Esame");
    }

    fBkgOnlyLegend = std::make_unique<TLegend>(0.19, 0.78, 0.37, 0.93);
    fBkgOnlyLegend->AddEntry(fDataMinusBkg.get(), "Data - background", "p");
    fBkgOnlyLegend->AddEntry(fSigStack.get(), "Signal", "f");
    fBkgOnlyLegend->AddEntry(fBkgOnlyErr.get(), "Post-fit bkg. uncertainty", "f");
    fBkgOnlyLegend->SetFillStyle(0);
    fBkgOnlyLegend->SetBorderSize(0);
    fBkgOnlyLegend->SetTextFont(gStyle->GetTextFont());
    fBkgOnlyLegend->SetTextSize(gStyle->GetTextSize());
    fBkgOnlyLegend->Draw("same");

    if(fPlotLabel!="none") TRExLabel(0.61,0.89, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str());
    myText(0.61,0.85,1,Form("#sqrt{s} = %s, %s",fCME.c_str(),fLumi.c_str()));
    for(unsigned int i_lab=0;i_lab<fLabels.size();i_lab++){
        myText(0.61,0.90-(i_lab+2)*0.04,1,Form("%s",fLabels[i_lab].c_str()));//,0.045);
    }

    // protection against bad settings for BinTickMarks
    bool doTickMarks = !fBinTickMarks.empty();
    if (doTickMarks && (fBinTickMarks.back() > fSigStack->GetNbinsX()+1 || fBinTickMarks.front() < 1)) {
        LOG(WARNING) << "Inconsistent `BinTickMarks`: first value cannot be lower than 1, last value cannot be larger than number of bins +1. Ignoring.\n";
        doTickMarks = false;
    }

    // Make custom tick marks for specific bins, if option BinTickMarks is used
    if(doTickMarks) CustomiseTickMarks(fSigStack.get());
}

// Customise the tick marks of the X-axis of the provided histogram, using fBinTickMarks and fBinTickMarksLabels
void TRExPlot::CustomiseTickMarks(TH1 *h){
    // first store necessary properties
    auto labelSize   = h->GetXaxis()->GetLabelSize();
    auto labelOffset = h->GetXaxis()->GetLabelOffset();
    auto labelFont   = h->GetXaxis()->GetLabelFont();
    auto titleOffset = h->GetXaxis()->GetTitleOffset();
    auto titleFont   = h->GetXaxis()->GetTitleFont();
    auto titleSize   = h->GetXaxis()->GetTitleSize();
    auto title       = h->GetXaxis()->GetTitle();
    auto ymin = static_cast<TH1*>(h->GetYaxis()->GetParent())->GetMinimum();
    auto ymax = static_cast<TH1*>(h->GetYaxis()->GetParent())->GetMaximum();
    auto xmin = h->GetXaxis()->GetBinLowEdge(1);
    auto xmax = h->GetXaxis()->GetBinUpEdge(h->GetNbinsX());

    // if no fBinTickMarksLabels given, make it using the bins low edges
    if (fBinTickMarksLabels.size() == 0) {
        for (const auto& ibin : fBinTickMarks) {
            std::ostringstream bin_label;
            bin_label << h->GetBinLowEdge(ibin);
            fBinTickMarksLabels.push_back(bin_label.str());
        }
    }

    // add the first and last ticks if missing, but make sure we won't show them
    if (fBinTickMarks.front()!=1) {
        fBinTickMarks.insert(fBinTickMarks.begin(),1);
        fBinTickMarksLabels.insert(fBinTickMarksLabels.begin(),std::string(""));// empty label to be used later
    }
    if (fBinTickMarks.back()!=h->GetNbinsX()+1) {
        fBinTickMarks.push_back(h->GetNbinsX()+1);
        fBinTickMarksLabels.push_back(std::string(""));// empty label to be used later
    }

    // first remove automatic tick marks
    h->SetNdivisions(0);

    // define picewise linear function, mapping the low bin edges values to the new tick mark indexes
    std::ostringstream faxis_formula;
    for (unsigned int tick = 0; tick < fBinTickMarks.size()-1; tick++) {
        faxis_formula << "(x<=" << tick+1 << ")?("
        << h->GetBinLowEdge(fBinTickMarks.at(tick)) <<"+"<<
        (h->GetBinLowEdge(fBinTickMarks.at(tick+1))-h->GetBinLowEdge(fBinTickMarks.at(tick))) << "*(x-" << tick << ")):";
    }
    faxis_formula << "(" << h->GetBinLowEdge(fBinTickMarks.back()) <<"+(x-"<<fBinTickMarks.size()-1<<"))";

    // define TF1 based on this formula, to be used in the TGaxis
    // no need to do it again if it is already defined
    if (fTickMarksFunction==nullptr) fTickMarksFunction = std::shared_ptr<TF1>(new TF1("fTickMarksFunction",faxis_formula.str().c_str(),0,fBinTickMarks.size()-1));

    // define the new axis and Draw it
    TGaxis * g = new TGaxis(xmin,ymin,xmax,ymin,"fTickMarksFunction",fBinTickMarks.size()-1,"");
    g->SetLabelSize(labelSize);

    // change the labels if they are specified by BinTickMarksLabels
    for (unsigned int itick=0; itick<fBinTickMarksLabels.size(); itick++) {
        auto s=fBinTickMarksLabels.at(itick);
        // the first tick mark start at 1
        // size=0 (third argument) means removing the label
        if (s == "") g->ChangeLabel(itick+1,-1,0,-1,-1,-1);
        else g->ChangeLabel(itick+1,-1,-1,-1,-1,-1,s.c_str());
    }

    // apply all the saved axis properties, and draw it
    g->SetLabelSize(labelSize);
    // increase offset to avoid overlap with Y-axis, when adding customised labels
    if (fBinTickMarksLabels.size() != 0) g->SetLabelOffset(labelOffset*1.35);
    else g->SetLabelOffset(labelOffset);
    g->SetLabelFont(labelFont);
    g->SetTitleOffset(titleOffset);
    g->SetTitleFont(titleFont);
    g->SetTitleSize(titleSize);
    g->SetTitle(title);
    g->Draw();

    // the same new axis, but on the top of the plot
    TGaxis * gg = new TGaxis(xmin,ymax,xmax,ymax,"fTickMarksFunction",fBinTickMarks.size()-1,"-");
    gg->SetLabelSize(0);
    gg->Draw();
}

//_____________________________________________________________________________
// function to get asymmetric error bars for hists (Used in WZ observation)
double GC_up(double data) {
    if (data == 0 ) return 0;
    return 0.5*TMath::ChisquareQuantile(1.-0.1586555,2.*(data+1))-data;
}

//_____________________________________________________________________________
//
double GC_down(double data) {
    if (data == 0 ) return 0;
    return data-0.5*TMath::ChisquareQuantile(0.1586555,2.*data);
}

//_____________________________________________________________________________
//
std::unique_ptr<TGraphAsymmErrors> poissonize(const TH1 *h) {
    std::unique_ptr<TGraphAsymmErrors> gr = std::make_unique<TGraphAsymmErrors>(h);
    for (int i = 0; i < gr->GetN(); i++) {
        double content = gr->GetErrorYhigh(i) * gr->GetErrorYhigh(i); // this to fix the case of the merged plots, where histograms (even data) are scaled; so the actual content is the square of the stat. error (right?)
        gr->SetPointError(i,0.499*h->GetBinWidth(i+1),0.5*h->GetBinWidth(i+1),GC_down(content),GC_up(content));
        if(h->GetBinContent(i+1)==0){
            gr->SetPoint(i,gr->GetX()[i],-1);
            gr->SetPointError(i,0,0,0,0);
        }
    }
    gr->SetMarkerSize(h->GetMarkerSize());
    gr->SetMarkerColor(h->GetMarkerColor());
    gr->SetMarkerStyle(h->GetMarkerStyle());
    gr->SetLineWidth(h->GetLineWidth());
    gr->SetLineColor(h->GetLineColor());
    gr->SetLineStyle(h->GetLineStyle());
    return gr;
}

//_____________________________________________________________________________
//
std::unique_ptr<TGraphAsymmErrors> histToGraph(const TH1* h){
    std::unique_ptr<TGraphAsymmErrors> gr = std::make_unique<TGraphAsymmErrors>(h);
    for (int i = 0; i < gr->GetN(); i++) {
        gr->SetPointEXlow(i,0.499*h->GetBinWidth(i+1));
        gr->SetPointEXhigh(i,0.5*h->GetBinWidth(i+1));
        if(h->GetBinContent(i+1)==0){
            gr->SetPoint(i,gr->GetX()[i],-1);
            gr->SetPointError(i,0,0,0,0);
        }
    }
    gr->SetMarkerStyle(h->GetMarkerStyle());
    gr->SetMarkerSize(h->GetMarkerSize());
    gr->SetMarkerColor(h->GetMarkerColor());
    gr->SetLineWidth(h->GetLineWidth());
    gr->SetLineColor(h->GetLineColor());
    gr->SetLineStyle(h->GetLineStyle());
    return gr;
}

//_____________________________________________________________________________
//
void SetHistBinWidth(TH1* h,double width){
    static const double epsilon = 0.00000001;
    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
        if(std::abs(h->GetBinWidth(i_bin)-width)>epsilon){
            h->SetBinContent(i_bin,h->GetBinContent(i_bin)*width/h->GetBinWidth(i_bin));
            h->SetBinError(  i_bin,h->GetBinError(i_bin)  *width/h->GetBinWidth(i_bin));
        }
    }
}

//_____________________________________________________________________________
//
void SetGraphBinWidth(TGraphAsymmErrors* g,double width){
    static const double epsilon = 0.00000001;
    for(int i_bin=0;i_bin<g->GetN();i_bin++){
        const double w = g->GetErrorXhigh(i_bin)+g->GetErrorXlow(i_bin);
        if(std::abs(w-width)>epsilon){
            g->SetPoint(      i_bin,g->GetX()[i_bin], g->GetY()[i_bin]*width/w);
            g->SetPointEYhigh(i_bin,g->GetErrorYhigh(i_bin)*width/w);
            g->SetPointEYlow( i_bin,g->GetErrorYlow(i_bin) *width/w);
        }
    }
}
