// Class include
#include "TRExFitter/CorrelationMatrix.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/YamlConverter.h"

// ROOT includes
#include "TCanvas.h"
#include "TH2D.h"
#include "TPad.h"
#include "TStyle.h"

// Style stuff
#include "StyleUtils/TRExLabels.h"

// c++ stuff
#include <algorithm>

//__________________________________________________________________________________
//
CorrelationMatrix::CorrelationMatrix() :
    fOutFolder(""),
    fPlotLabel("") {
}

//__________________________________________________________________________________
//
void CorrelationMatrix::AddNuisPar(const std::string& p){
    if (std::find(fNuisParNames.begin(), fNuisParNames.end(), p) != fNuisParNames.end()) {
        LOG(WARNING) << "Parameter " << p << " already in the list\n";
        return;
    }
    fNuisParNames.push_back(p);
}

//__________________________________________________________________________________
//
void CorrelationMatrix::Resize(const int size) {
    fMatrix.resize(size);
    for (auto& i : fMatrix) {
        i.resize(size);
    }
}

//__________________________________________________________________________________
//
void CorrelationMatrix::SetCorrelation(const std::string& p0, const std::string& p1, double corr){
    auto itr0 = std::find(fNuisParNames.begin(), fNuisParNames.end(), p0);
    auto itr1 = std::find(fNuisParNames.begin(), fNuisParNames.end(), p1);
    if (itr0 == fNuisParNames.end() || itr1 == fNuisParNames.end()) {
        LOG(ERROR) << "Cannot find parameter " << p0 << " or parameter " << p1 << "\n";
        return;
    }
    const std::size_t idx0 = std::distance(fNuisParNames.begin(), itr0);
    const std::size_t idx1 = std::distance(fNuisParNames.begin(), itr1);
    fMatrix[idx0][idx1] = corr;
}

//__________________________________________________________________________________
//
double CorrelationMatrix::GetCorrelation(const std::string& p0, const std::string& p1) const {
    bool isMorph_p0 = false;
    bool isMorph_p1 = false;
    if (p0.find("morph_") != std::string::npos) isMorph_p0 = true;
    if (p1.find("morph_") != std::string::npos) isMorph_p1 = true;
    auto itr0 = std::find(fNuisParNames.begin(), fNuisParNames.end(), p0);
    auto itr1 = std::find(fNuisParNames.begin(), fNuisParNames.end(), p1);
    // if one of the two is missing, return 1 or 0 (if name1==name2 ==> 1, not zero!)
    if(itr0 == fNuisParNames.end()) {
        if(!isMorph_p0) LOG(VERBOSE) << "NP " << p0 << " not found in correlation matrix. Returning correlation = " << (1.*(p0==p1)) << "\n";
        else            LOG(VERBOSE) << "NP " << p0 << " not found in correlation matrix. The NP is for morphing. Returning correlation = " << (1.*(p0==p1)) << "\n";
        return (p0 == p1) ? 1. : 0.;
    }
    if(itr1 == fNuisParNames.end()){
        if(!isMorph_p1) LOG(VERBOSE) << "NP " << p1 << " not found in correlation matrix. Returning correlation = " << (1.*(p0==p1)) << "\n";
        else            LOG(VERBOSE) << "NP " << p0 << " not found in correlation matrix. The NP is for morphing. Returning correlation = " << (1.*(p0==p1)) << "\n";
        return (p0 == p1) ? 1. : 0.;
    }
    const std::size_t idx0 = std::distance(fNuisParNames.begin(), itr0);
    const std::size_t idx1 = std::distance(fNuisParNames.begin(), itr1);
    return fMatrix[idx0][idx1];
}


//__________________________________________________________________________________
//
void CorrelationMatrix::Draw(const std::string& path, const bool useHEPDataFormat, const double minCorr, const std::vector<std::string>& POIs){
    //
    // 0) Determines the number of lines/columns
    //

    std::vector<std::vector<double> > correlations(fNuisParNames.size(), std::vector<double>(fNuisParNames.size()));

    std::vector <std::string> vec_NP;
    for(unsigned int iNP = 0; iNP < fNuisParNames.size(); ++iNP){
        const std::string iSystName = fNuisParNames.at(iNP);
        for(unsigned int jNP = 0; jNP < fNuisParNames.size(); ++jNP){
            const std::string jSystName = fNuisParNames.at(jNP);
            const double corr = GetCorrelation(iSystName, jSystName);

            correlations.at(iNP).at(jNP) = corr;

            if(jNP == iNP) continue;
            // if EFT mode just skip non-EFT items
            if(!fEFTParList.empty()){
              auto it_isys = std::find(fEFTParList.begin(), fEFTParList.end(), iSystName);
              if (it_isys == fEFTParList.end()) continue;
              auto it_jsys = std::find(fEFTParList.begin(), fEFTParList.end(), jSystName);
              if (it_jsys == fEFTParList.end()) continue;
            }
            // print NPs that have correlations above threshold
            if (minCorr > 0) {
                if(std::abs(corr)>=minCorr){
                    auto itr = std::find(vec_NP.begin(), vec_NP.end(), iSystName);
                    if (itr != vec_NP.end()) continue;
                    LOG(VERBOSE) << iSystName << " " << minCorr << "    " << corr << " (" << jSystName << ")\n";
                    vec_NP.push_back(iSystName);
                }
            } else {
                auto itr = std::find(vec_NP.begin(), vec_NP.end(), iSystName);
                if (itr != vec_NP.end()) continue;
                vec_NP.push_back(iSystName);
            }
        }
    }
    // also print the correlation of all POIs with the selected NPs that are above threshold
    for (const auto& poi: POIs) {
        if (std::find(vec_NP.begin(), vec_NP.end(), poi) == vec_NP.end()) vec_NP.push_back(poi);
    }
    int N = vec_NP.size();

    // pass everything to yaml converter
    YamlConverter converter{};
    if(!fEFTParList.empty()) converter.WriteCorrelation(fNuisParNames, correlations, fOutFolder, false,  "EFT_CorrelationMatrix");
    else converter.WriteCorrelation(fNuisParNames, correlations, fOutFolder, false);

    if (useHEPDataFormat) {
        if(!fEFTParList.empty()) converter.WriteCorrelationHEPData(fNuisParNames, correlations, fOutFolder, false, "EFT_Correlation");
        else converter.WriteCorrelationHEPData(fNuisParNames, correlations, fOutFolder, false);
    }

    //
    // 0.5) Skip some NPs
    //
    std::vector <std::string> vec_NP_old = vec_NP;
    vec_NP.clear();
    for(unsigned int iNP = 0; iNP < vec_NP_old.size(); ++iNP) {
        const std::string iSystName = vec_NP_old[iNP];
        bool skip(false);
        if (iSystName.find("Expression_") != std::string::npos) {
            skip = true;
        }
        if(skip) continue;
        if(Common::FindInStringVector(fNuisParToHide,iSystName) >=0) continue;
        vec_NP.push_back(iSystName);
    }
    N = vec_NP.size();

    //
    // 0.75) Reorder NPs
    if(fNuisParList.size()>0){
        vec_NP_old = vec_NP;
        vec_NP.clear();
        for(const auto& npName : fNuisParList){
            for(unsigned int iNP = 0; iNP < vec_NP_old.size(); ++iNP){
                if(vec_NP_old[iNP]==npName) vec_NP.emplace_back(vec_NP_old[iNP]);
            }
        }
        N = vec_NP.size();
    }

    //
    // 1) Performs the plot
    //
    TH2D h_corr("h_corr","",N,0,N,N,0,N);
    h_corr.SetDirectory(nullptr);

    for(unsigned int iNP = 0; iNP < vec_NP.size(); ++iNP){//line number
        const std::string iSystName = vec_NP[iNP];

        auto itr = TRExFitter::SYSTMAP.find(iSystName);
        if(itr != TRExFitter::SYSTMAP.end()){
            h_corr.GetXaxis()->SetBinLabel(iNP+1,itr->second.c_str());
            h_corr.GetYaxis()->SetBinLabel(N-iNP,itr->second.c_str());
        }
        else{
            h_corr.GetXaxis()->SetBinLabel(iNP+1,iSystName.c_str());
            h_corr.GetYaxis()->SetBinLabel(N-iNP,iSystName.c_str());
        }

        for(unsigned int jNP = 0; jNP < vec_NP.size(); ++jNP){//column number
            const std::string jSystName = vec_NP[jNP];

            h_corr.SetBinContent(iNP+1,N-jNP,100.*GetCorrelation(iSystName, jSystName));

        }
    }
    h_corr.SetMinimum(-100.);
    h_corr.SetMaximum(100.);

    int size = 200;
    if(vec_NP.size()>4){
        size = vec_NP.size()*50;
    }

    double markerSize = vec_NP.size() > 10 ? 0.5 : 1;
    if (vec_NP.size() > 50) markerSize = 0.2;

    //
    // 2) Style settings
    //
    TCanvas c1("","",0.,0.,size+300,size+300);
    gStyle->SetPalette(87);
    // cppcheck-suppress syntaxError
    #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
    h_corr.SetMarkerSize(0.75*1000);
    #else
    h_corr.SetMarkerSize(markerSize);
    #endif
    gStyle->SetPaintTextFormat(".1f");
    gPad->SetLeftMargin(240./(size+300));
    gPad->SetBottomMargin(240./(size+300));
    gPad->SetRightMargin(60./(size+300));
    gPad->SetTopMargin(60./(size+300));

    h_corr.GetXaxis()->LabelsOption("v");
    h_corr.GetXaxis()->SetLabelSize( h_corr.GetXaxis()->GetLabelSize()*0.75 );
    h_corr.GetYaxis()->SetLabelSize( h_corr.GetYaxis()->GetLabelSize()*0.75 );
    c1.SetTickx(0);
    c1.SetTicky(0);
    h_corr.GetYaxis()->SetTickLength(0);
    h_corr.GetXaxis()->SetTickLength(0);
    c1.SetGrid();
    h_corr.Draw("col TEXT");

    if (fPlotLabel != "none") {
        TRExLabelNew( 30./(size+300) , (size+240.+15.)/(size+300.), TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str(), kBlack, gStyle->GetTextSize(), (80. / (size+300)));
    }

    c1.RedrawAxis("g");

    if (!fEFTParList.empty()) {
        Common::SaveCanvasAs(c1, std::string(((TString)path).TString::ReplaceAll("CorrMatrix","EFT_CorrMatrix")));
    } else {
        Common::SaveCanvasAs(c1, path);
    }
}

//__________________________________________________________________________________
//
bool CorrelationMatrix::NPIsPresent(const std::string& np) const {
    auto itr = std::find(fNuisParNames.begin(), fNuisParNames.end(), np);
    return itr != fNuisParNames.end();
}
