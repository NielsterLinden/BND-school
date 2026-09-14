#include "TRExFitter/RankingManager.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Region.h"

#include "StyleUtils/TRExLabels.h"
#include "StyleUtils/TRExUtils.h"

#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooSimultaneous.h"

#include "xRooFit/xRooNode.h"

#include <algorithm>
#include <fstream>

//__________________________________________________________________________________
//
RankingManager::RankingManager() :
    fOutputPath(""),
    fInjectGlobalObservables(false),
    fFitStrategy(1),
    fRegularizationType(0),
    fCPU(1),
    fStatOnly(false),
    fPlotLabel("Internal"),
    fLumiLabel("139 fb^{-1}"),
    fCmeLabel("13 TeV"),
    fHEPDataFormat(false),
    fName("MyFit"),
    fSuffix(""),
    fRankingMaxNP(9999),
    fRankingPOIName(fName),
    fRankingPOIAxisScale(1),
    fUsePOISinRanking(false),
    fUseHesseBeforeMigrad(false),
    fSetNdivisions(510),
    fRandomNP(0.1),
    fRandomize(false),
    fRandomSeed(1234567),
    fMaximumNumberOfFCNcalls(-1),
    fToleranceScale(1.),
    fShiftGlobalObservables(false),
    fIsCovarianceBreakdown(false)
{
}

//__________________________________________________________________________________
//
void RankingManager::AddNuisPar(const std::string& name, const bool isNF) {

    // check if the NP already exists by checking the names
    auto it = std::find_if(fNuisPars.begin(), fNuisPars.end(),
                [&name](const std::pair<std::string, bool>& element){return element.first == name;});

    if (it != fNuisPars.end()) {
        LOG(WARNING) << "NP " << name << " already exists in the list, not adding it\n";
        return;
    }
    fNuisPars.emplace_back(name, isNF);
}

//__________________________________________________________________________________
//
void RankingManager::RunRanking(const std::shared_ptr<xRooNode>& pdf,
                                const std::shared_ptr<xRooNode>& pdfPrefit,
                                RooWorkspace* ws,
                                RooDataSet* data,
                                const std::vector<std::shared_ptr<NormFactor> >& nfs) const {

    if (fOutputPath == "") {
        LOG(ERROR) << "OutputPath not set, plese set it via SetOutputPath()\n";
        exit(EXIT_FAILURE);
    }

    if (!ws) {
        LOG(ERROR) << "Workspace is nullptr\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::ofstream> outFiles;
    for (const auto& poi : fPOINames) {
        outFiles.emplace_back((fOutputPath+"_"+poi+".txt").c_str());

        if (!outFiles.back().good() || !outFiles.back().is_open()) {
            LOG(ERROR) << "Cannot open file at " << fOutputPath << "_" << poi << ".txt\n";
            exit(EXIT_FAILURE);
        }
    }

    RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc){
        LOG(ERROR) << "ModelConfig is missing\n";
        exit(EXIT_FAILURE);
    }
    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());
    if (!simPdf || !data){
        LOG(ERROR) << "RooSimultaneous or data is missing\n";
        exit(EXIT_FAILURE);
    }

    if (fInjectGlobalObservables && !fFitValues.empty()) {
        FitUtils::InjectGlobalObservables(ws, fFitValues);
    }

    std::unique_ptr<RooArgSet> params(mc->GetPdf()->getParameters(data));

    ws->saveSnapshot("tmp_snapshot", *params);

    FittingTool fitTool{};
    fitTool.SetUseHesse(false);
    fitTool.SetUseHesseBeforeMigrad(fUseHesseBeforeMigrad);
    fitTool.SetStrategy(fFitStrategy);
    fitTool.SetRandomNP(fRandomNP, fRandomize, fRandomSeed);
    fitTool.SetMaxFCNcalls(fMaximumNumberOfFCNcalls);
    fitTool.SetToleranceScale(fToleranceScale);

    std::vector<std::string> constNFs;
    for(const auto& inf : nfs) {
        fitTool.AddValPOI(inf->fName, inf->GetNominal());
        if (inf->fConst) {
            constNFs.emplace_back(inf->fName);
        }
    }
    fitTool.SetNCPU(fCPU);
    fitTool.SetConstantNFs(constNFs);
    fitTool.ConstPOI(false);
    if(fStatOnly){
        fitTool.NoGammas();
        fitTool.NoSystematics();
    }

    std::vector<std::string> npNames;
    std::vector<double> npValues;
    for(const auto& inf : nfs) {
        if (!fUsePOISinRanking && Common::FindInStringVector(fPOINames, inf->fName) >= 0) continue;
        npNames. emplace_back(inf->fName);
        npValues.emplace_back(inf->GetNominal());
    }
    fitTool.SetNPs(npNames,npValues);

    FitUtils::ApplyExternalConstraints(ws, &fitTool, simPdf, nfs, fRegularizationType);

    std::vector<double> muhats;
    for (const auto& poi : fPOINames) {
        auto par = pdf->pars().find(poi);
        if (!par) {
            LOG(ERROR) << "Cannot find parameter " << poi << "\n";
            exit(EXIT_FAILURE);
        }
        muhats.emplace_back(par->get<RooRealVar>()->getVal());
    }
    for (const auto& iNP : fNuisPars){

        const std::string& npName = iNP.first;
        if (!fUsePOISinRanking && Common::FindInStringVector(fPOINames, npName) >= 0)  continue;
        auto par = pdf->pars().find(npName);
        if (!par) {
            LOG(ERROR) << "Cannot find parameter " << npName << "\n";
            exit(EXIT_FAILURE);
        }
        const auto parVar = par->get<RooRealVar>();

        //
        // Getting the postfit values of the nuisance parameter
        const double central = parVar->getVal();
        const double up      = parVar->getErrorHi();
        const double down    = parVar->getErrorLo();

        auto tmp = pdfPrefit->pars().find(iNP.first);
        if (!tmp) {
            LOG(ERROR) << "Cannot find parameter " << iNP.first << " for prefits\n";
            exit(EXIT_FAILURE);
        }
        const double centralPrefit = tmp->get<RooRealVar>()->getVal();

        RankingManager::RankingValues values;
        values.central = central;
        values.up = up;
        values.down = down;
        if (iNP.first.find("gamma_") != std::string::npos) {
            values.upPrefit = tmp->get<RooRealVar>()->getErrorHi();
            values.downPrefit = tmp->get<RooRealVar>()->getErrorLo();

        } else {
            values.upPrefit = 1.;
            values.downPrefit = -1.;
        }

        for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
            if (iNP.first == fPOINames.at(ipoi)) continue;

            if (iNP.first.find("gamma_") != std::string::npos) {
                outFiles.at(ipoi) << iNP.first << "   " << (central - 1.)/( central>=1. ? std::abs(values.upPrefit) : std::abs(values.downPrefit) ) << " +" << std::abs(up/values.upPrefit) << " -" << std::abs(down/values.downPrefit)<< "  ";
            } else {
                outFiles.at(ipoi) << Common::ReplaceString(iNP.first, "alpha_", "") << "   " << central - centralPrefit << " +" << std::abs(up) << " -" << std::abs(down)<< "  ";
            }
        }

        std::vector<double> dMuUp   = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, true, false, values, muhats);
        std::vector<double> dMuDown = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, false, false, values, muhats);

        for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
            if (iNP.first == fPOINames.at(ipoi)) continue;

            outFiles.at(ipoi) << dMuUp.at(ipoi) << "   " << dMuDown.at(ipoi) << "  ";
        }

        dMuUp   = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, true, true, values, muhats);
        dMuDown = RunSingleFit(&fitTool, ws, mc, simPdf, data, iNP, false, true, values, muhats);

        for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
            if (iNP.first == fPOINames.at(ipoi)) continue;

            outFiles.at(ipoi) << dMuUp.at(ipoi) << "   " << dMuDown.at(ipoi) << " " << std::endl;
        }

    }

    ws->loadSnapshot("tmp_snapshot");
    for (auto& ifile : outFiles) {
        ifile.close();
    }
}

//__________________________________________________________________________________
//
std::vector<double> RankingManager::RunSingleFit(FittingTool* fitTool,
                                                 RooWorkspace* ws,
                                                 RooStats::ModelConfig *mc,
                                                 RooSimultaneous *simPdf,
                                                 RooDataSet* data,
                                                 const std::pair<std::string, bool>& np,
                                                 const bool isUp,
                                                 const bool isPrefit,
                                                 const RankingManager::RankingValues& values,
                                                 const std::vector<double>& muhats) const {

    if (isPrefit && np.second) {
        std::vector<double> tmp(muhats.size(), 0);
        return tmp;
    }

    if (isPrefit && fShiftGlobalObservables) {
        std::vector<double> tmp(muhats.size(), 0);
        return tmp;
    }

    if (np.second && fShiftGlobalObservables) {
        std::vector<double> tmp(muhats.size(), 0);
        return tmp;
    }

    ws->loadSnapshot("tmp_snapshot");
    fitTool->ResetFixedNP();
    if(fFitFixedNPs.size()>0){
        for(const auto& nuisParToFix : fFitFixedNPs){
            fitTool->FixNP(nuisParToFix.first,nuisParToFix.second);
        }
    }

    double shift = 1.0;
    if (isPrefit) {
       shift = isUp ? values.upPrefit : values.downPrefit;
    } else {
       shift = isUp ? std::abs(values.up) : -std::abs(values.down);
    }

    // Experimental: reduce the range of ranking
    if(TRExFitter::OPTION["ReduceRanking"]!=0){
        shift *= TRExFitter::OPTION["ReduceRanking"];
    }

    if (fShiftGlobalObservables) {
        // shift global observables
        // will be restored by the snapshot call
        auto go = ws->var(("nom_"+np.first).c_str());
        if (!go) {
            LOG(ERROR) << "Cannot find global observable for NP: " << np.first << "\n";
            return std::vector<double>(muhats.size(), 0);
        }
        const double goVal = go->getVal();
        double goShift(isUp ? 1. : -1.);
        if (np.first.find("gamma_") != std::string::npos) {
            goShift = (goVal + (isUp ? 1./values.upPrefit : 1./values.downPrefit));
        }
        LOG(INFO) << "Shifting GlobalObservable for NP: " << np.first << " from: " << goVal << " to: " << goShift << "\n";
        go->setVal(goShift);
    } else {
        fitTool->FixNP(np.first, values.central + shift);
    }

    fitTool->FitPDF(mc, simPdf, data);

    std::vector<double> result(fPOINames.size());
    for (std::size_t ipoi = 0; ipoi < fPOINames.size(); ++ipoi) {
        result.at(ipoi) = fitTool->ExportFitResultInMap()[fPOINames.at(ipoi)] - muhats.at(ipoi);
    }

    if(TRExFitter::OPTION["ReduceRanking"]!=0){
        for (auto& i : result) {
            i /= TRExFitter::OPTION["ReduceRanking"];
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void RankingManager::PlotRanking(const std::vector<std::unique_ptr<Region> >& regions,
                                 const std::vector<std::string>& nfs,
                                 const std::vector<std::string>& sfs,
                                 const bool flagSysts,
                                 const bool flagGammas) {

    std::vector<Region*> reg;
    for (const auto& ireg : regions) {
        reg.push_back(ireg.get());
    }

    this->PlotRanking(reg, nfs, sfs, flagSysts, flagGammas);

}

//__________________________________________________________________________________
//
void RankingManager::PlotRanking(const std::vector<Region* >& regions,
                                 const std::vector<std::string>& nfs,
                                 const std::vector<std::string>& sfs,
                                 const bool flagSysts,
                                 const bool flagGammas) {

    std::size_t maxNP = fRankingMaxNP;
    std::vector<std::string> parname;
    std::vector<double> nuhat;
    std::vector<double> nuerrhi;
    std::vector<double> nuerrlo;
    std::vector<double> poiup;
    std::vector<double> poidown;
    std::vector<double> poinomup;
    std::vector<double> poinomdown;

    gSystem->mkdir((fName+"/Rankings/").c_str());

    if (fRankingContainer.empty()) {
        LOG(DEBUG) << "No ranking container set, will not produce ranking\n";
        return;
    }

    for (const auto& icontainer : fRankingContainer) {
        parname.emplace_back(icontainer.name);
        nuhat.emplace_back(icontainer.nphat);
        nuerrhi.emplace_back(icontainer.nperrhi);
        nuerrlo.emplace_back(icontainer.nperrlo);
        poiup.emplace_back(icontainer.poihi);
        poidown.emplace_back(icontainer.poilo);
        poinomup.emplace_back(icontainer.poiprehi);
        poinomdown.emplace_back(icontainer.poiprelo);
    }

    {

        std::sort(fRankingContainer.begin(), fRankingContainer.end(),
            [](const YamlConverter::RankingContainer& lhs, const YamlConverter::RankingContainer& rhs)
              {return (std::max(std::abs(lhs.poihi),std::abs(lhs.poilo))) > (std::max(std::abs(rhs.poihi),std::abs(rhs.poilo)));});

        fRankingContainer.resize(std::min(maxNP, fRankingContainer.size()));

        YamlConverter converter{};
        converter.WriteRanking(fRankingContainer, fName+"/Rankings/Ranking"+fSuffix+".yaml", fShiftGlobalObservables);
        if (fHEPDataFormat) {
            converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
            converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
            converter.WriteRankingHEPData(fRankingContainer, fName, fSuffix, fShiftGlobalObservables);
        }
    }

    unsigned int SIZE = parname.size();
    LOG(DEBUG) << "NP ordering...\n";
    std::vector<double> number;
    number.push_back(0.5);
    for (unsigned int i=1;i<SIZE;i++){
        number.push_back(i+0.5);
        double sumi = 0.0;
        int index=-1;
        sumi += std::max( std::abs(poiup[i]),std::abs(poidown[i]) );
        for (unsigned int j=1;j<=i;j++){
            double sumii = 0.0;
            sumii += std::max(std::abs(poiup[i-j]),std::abs(poidown[i-j]) );
            if (sumi<sumii){
                if (index==-1){
                    std::swap(poiup[i],poiup[i-j]);
                    std::swap(poidown[i],poidown[i-j]);
                    std::swap(poinomup[i],poinomup[i-j]);
                    std::swap(poinomdown[i],poinomdown[i-j]);
                    std::swap(nuhat[i],nuhat[i-j]);
                    std::swap(nuerrhi[i],nuerrhi[i-j]);
                    std::swap(nuerrlo[i],nuerrlo[i-j]);
                    std::swap(parname[i],parname[i-j]);
                    index=i-j;
                }
                else{
                    std::swap(poiup[index],poiup[i-j]);
                    std::swap(poidown[index],poidown[i-j]);
                    std::swap(poinomup[index],poinomup[i-j]);
                    std::swap(poinomdown[index],poinomdown[i-j]);
                    std::swap(nuhat[index],nuhat[i-j]);
                    std::swap(nuerrhi[index],nuerrhi[i-j]);
                    std::swap(nuerrlo[index],nuerrlo[i-j]);
                    std::swap(parname[index],parname[i-j]);
                    index=i-j;
                }
            }
            else{
                break;
            }
        }
    }
    number.push_back(parname.size()-0.5);

    double poimax = 0;
    for (unsigned int i=0;i<SIZE;i++) {
        poimax = std::max(poimax,std::max(std::abs(poiup[i]),std::abs(poidown[i]) ));
        poimax = std::max(poimax,std::max(std::abs(poinomup[i]),std::abs(poinomdown[i]) ));
        nuerrlo[i] = std::abs(nuerrlo[i]);
    }
    poimax *= 1.2;

    for (unsigned int i=0;i<SIZE;i++) {
        poiup[i]     *= (2./poimax);
        poidown[i]   *= (2./poimax);
        poinomup[i]  *= (2./poimax);
        poinomdown[i]*= (2./poimax);
    }

    poimax = poimax/fRankingPOIAxisScale;

    // Restrict to the first N
    if(SIZE>maxNP) SIZE = maxNP;

    // Graphical part - rewritten taking DrawPulls in TRExFitter
    double lineHeight  =  30.;
    double offsetUp    =  60.; // external
    double offsetDown  =  60.;
    double offsetUp1   = 100.; // internal
    double offsetDown1 =  15.;
    int offset = offsetUp + offsetDown + offsetUp1 + offsetDown1;
    int newHeight = offset + SIZE*lineHeight;

    double xmin = -2.;
    double xmax =  2.;
    double max  =  0.;

    TGraphAsymmErrors g{};
    TGraphAsymmErrors g1{};
    TGraphAsymmErrors g2{};
    TGraphAsymmErrors g1a{};
    TGraphAsymmErrors g2a{};

    int idx = 0;
    std::vector< std::string > Names;
    std::string parTitle;
    std::vector<std::unique_ptr<TMarker> > markers;

    for(unsigned int i = parname.size()-SIZE; i<parname.size(); ++i){
        g.SetPoint(idx, nuhat[i],  idx+0.5);
        g.SetPointEXhigh(      idx, nuerrhi[i]);
        g.SetPointEXlow(       idx, nuerrlo[i]);

        g1.SetPoint(      idx, 0.,idx+0.5);
        g1.SetPointEXhigh(idx, poiup[i]);
        g1.SetPointEXlow( idx, 0.);
        g1.SetPointEYhigh(idx, 0.4);
        g1.SetPointEYlow( idx, 0.4);

        g2.SetPoint(      idx, 0.,idx+0.5);
        g2.SetPointEXhigh(idx, poidown[i]);
        g2.SetPointEXlow( idx, 0.);
        g2.SetPointEYhigh(idx, 0.4);
        g2.SetPointEYlow( idx, 0.4);

        g1a.SetPoint(      idx, 0.,idx+0.5);
        g1a.SetPointEXhigh(idx, poinomup[i]);
        g1a.SetPointEXlow( idx, 0.);
        g1a.SetPointEYhigh(idx, 0.4);
        g1a.SetPointEYlow( idx, 0.4);

        g2a.SetPoint(      idx, 0.,idx+0.5);
        g2a.SetPointEXhigh(idx, poinomdown[i]);
        g2a.SetPointEXlow( idx, 0.);
        g2a.SetPointEYhigh(idx, 0.4);
        g2a.SetPointEYlow( idx, 0.4);

        double x,y;
        g.GetPoint(idx, x, y);
        // see if it is a NF
        auto itr = std::find(nfs.begin(), nfs.end(), parname.at(i));
        auto itrSF = std::find_if(sfs.begin(), sfs.end(), [&parname,&i](const auto& sf){return parname.at(i).find(sf) != std::string::npos;});
        if (itr != nfs.end()) {
            markers.emplace_back(std::make_unique<TMarker>(x,y,20));
            if (fShiftGlobalObservables) {
                markers.back()->SetMarkerColor(kViolet+2);
            } else {
                markers.back()->SetMarkerColor(kRed+1);
            }
            parTitle = TRExFitter::SYSTMAP[parname.at(i)];
        } else if (itrSF != sfs.end()) {
            markers.emplace_back(std::make_unique<TMarker>(x,y,20));
            markers.back()->SetMarkerColor(kRed+1);
            std::vector<std::string> tmpVec = Common::Vectorize(parname.at(i),'_');
            const std::string binName = tmpVec.back();
            parTitle = "ShapeFactor (bin " + binName + ")";
        } else if (parname[i].find("gamma")!=std::string::npos){
            markers.emplace_back(std::make_unique<TMarker>(x,y,24));
            markers.back()->SetMarkerColor(kBlack);
            // get name of the region
            std::vector<std::string> tmpVec = Common::Vectorize(parname.at(i),'_');
            std::size_t nWords = tmpVec.size();
            std::string regName = tmpVec.at(2);
            for (std::size_t i_word = 3; i_word < nWords-2; ++i_word) {
                regName += tmpVec.at(i_word);
            }
            // find the short label of this region
            std::string regTitle = regName;
            for(const auto* ireg : regions) {
                if (ireg->fName == regName){
                    regTitle = ireg->fShortLabel;
                    break;
                }
            }
            // build the title of the nuis par
            parTitle = "#gamma (" + regTitle + " bin " + tmpVec.at(nWords-1) + ")";
        } else {
            auto itrTitle = TRExFitter::SYSTMAP.find(parname.at(i));
            if (itrTitle == TRExFitter::SYSTMAP.end()) {
                parTitle = parname.at(i);
            } else {
                parTitle = itrTitle->second;
            }
            markers.emplace_back(std::make_unique<TMarker>(x,y,20));
            markers.back()->SetMarkerColor(kBlack);
        }

        Names.push_back(parTitle);

        idx ++;
        if(idx > max)  max = idx;
    }
    int newWidth = 600;
    if (fNPRankingCanvasSize.size() != 0){
        newWidth = fNPRankingCanvasSize.at(0);
        newHeight = fNPRankingCanvasSize.at(1);
    }
    TCanvas c("c","c",newWidth,newHeight);
    c.SetTicks(0,0);
    gPad->SetLeftMargin(0.4);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(1.*offsetUp/newHeight);
    gPad->SetBottomMargin(1.*offsetDown/newHeight);

    TH1D h_dummy("h_dummy","h_dummy",10,xmin,xmax);
    h_dummy.SetMaximum( SIZE + offsetUp1/lineHeight   );
    h_dummy.SetMinimum(      - offsetDown1/lineHeight );
    h_dummy.SetLineWidth(0);
    h_dummy.SetFillStyle(0);
    h_dummy.SetLineColor(kWhite);
    h_dummy.SetFillColor(kWhite);
    h_dummy.GetYaxis()->SetLabelSize(0);
    h_dummy.Draw();
    h_dummy.GetYaxis()->SetNdivisions(0);
    for(int i_bin=0;i_bin<h_dummy.GetNbinsX()+1;i_bin++){
        h_dummy.SetBinContent(i_bin,-10);
    }

    if (fShiftGlobalObservables) {
        if (fIsCovarianceBreakdown) {
            g1.SetFillColor(kGreen+2);
            g2.SetFillColor(kGreen-4);
        } else {
            g1.SetFillColor(kOrange+2);
            g2.SetFillColor(kOrange-2);
        }
    } else {
        g1.SetFillColor(kAzure-4);
        g2.SetFillColor(kCyan);
    }
    g1.SetLineColor(g1.GetFillColor());
    g2.SetLineColor(g2.GetFillColor());

    g1a.SetFillColor(kWhite);
    g2a.SetFillColor(kWhite);
    g1a.SetLineColor(kAzure-4);
    g2a.SetLineColor(kCyan);
    g1a.SetFillStyle(0);
    g2a.SetFillStyle(0);
    g1a.SetLineWidth(1);
    g2a.SetLineWidth(1);

    g.SetLineWidth(2);

    g.SetMarkerSize(0);

    if (!fShiftGlobalObservables) {
        g1a.Draw("5 same");
        g2a.Draw("5 same");
    }
    g1.Draw("2 same");
    g2.Draw("2 same");
    g.Draw("p same");

    TLatex systs{};
    systs.SetTextAlign(32);
    systs.SetTextSize( systs.GetTextSize()*0.8 );
    for(int i=0;i<max;i++){
        systs.DrawLatex(xmin-0.1,i+0.5,Names[i].c_str());
    }
    h_dummy.GetXaxis()->SetLabelSize( h_dummy.GetXaxis()->GetLabelSize()*0.9 );
    h_dummy.GetXaxis()->CenterTitle();
    h_dummy.GetXaxis()->SetTitle("(#hat{#theta}-#theta_{0})/#Delta#theta");
    h_dummy.GetXaxis()->SetTitleOffset(1.2);

    int POIScalePower = std::round(std::log10(fRankingPOIAxisScale));
    std::string POIScaleSuf = "";
    if(POIScalePower != 0){
        POIScaleSuf = Form(" #times 10^{%i}", POIScalePower);
    }

    TGaxis axis_up( -2, SIZE + (offsetUp1)/lineHeight, 2, SIZE + (offsetUp1)/lineHeight, -poimax,poimax, fSetNdivisions, "-" );
    axis_up.SetNdivisions(fSetNdivisions);
    axis_up.SetLabelOffset( 0.01 );
    axis_up.SetLabelSize(   h_dummy.GetXaxis()->GetLabelSize() );
    axis_up.SetLabelFont(   gStyle->GetTextFont() );
    axis_up.Draw();
    axis_up.CenterTitle();
    axis_up.SetTitle(("#Delta"+fRankingPOIName+POIScaleSuf).c_str());
    if((SIZE >= 15) || (parname.size() >= 15)) axis_up.SetTitleOffset(1.5);
    axis_up.SetTitleSize(   h_dummy.GetXaxis()->GetLabelSize() );
    axis_up.SetTitleFont(   gStyle->GetTextFont() );

    TPad pad1("p1","Pad High",0,(newHeight-offsetUp-offsetUp1)/newHeight,0.4,1);
    pad1.Draw();

    pad1.cd();
    TLegend leg1(0.02,0.7,1,1.0,("Pre-fit impact on "+fRankingPOIName+":").c_str());
    leg1.SetFillStyle(0);
    leg1.SetBorderSize(0);
    leg1.SetMargin(0.25);
    leg1.SetNColumns(2);
    leg1.SetTextFont(gStyle->GetTextFont());
    leg1.SetTextSize(gStyle->GetTextSize());
    leg1.AddEntry(&g1a,"#theta = #hat{#theta}+#Delta#theta","f");
    leg1.AddEntry(&g2a,"#theta = #hat{#theta}-#Delta#theta","f");
    if (!fShiftGlobalObservables) {
        leg1.Draw();
    }

    std::unique_ptr<TLatex> globalObs;
    if (fShiftGlobalObservables) {
        if (fIsCovarianceBreakdown) {
            globalObs = std::make_unique<TLatex>();
            globalObs->DrawLatex(0.03, 0.40,"#scale[0.9]{w/ #rho = linear correlation coef.}");
        } else {
            globalObs = std::make_unique<TLatex>();
            globalObs->DrawLatex(0.03, 0.40,"#scale[0.9]{w/ constraint term G(#font[12]{a}|#theta,#Delta#theta)}");
        }
    }
    std::unique_ptr<TLegend> leg2(nullptr);

    if (fShiftGlobalObservables) {
        leg2 = std::make_unique<TLegend>(0.02,0.52,1,0.9,("Unc. component on "+fRankingPOIName+":").c_str());
    } else {
        leg2 = std::make_unique<TLegend>(0.02,0.32,1,0.62,("Post-fit impact on "+fRankingPOIName+":").c_str());
    }
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.25);
    leg2->SetNColumns(2);
    leg2->SetTextFont(gStyle->GetTextFont());
    leg2->SetTextSize(gStyle->GetTextSize());
    if (fShiftGlobalObservables) {
        if (fIsCovarianceBreakdown) {
            leg2->AddEntry(&g1,"#sigma_{POI}^{up}#rho#sigma_{NP}^{up}","f");
            leg2->AddEntry(&g2,"-#sigma_{POI}^{dn}#rho#sigma_{NP}^{dn}","f");
        } else {
            leg2->AddEntry(&g1,"#font[12]{a}=#theta_{0}+#Delta#theta","f");
            leg2->AddEntry(&g2,"#font[12]{a}=#theta_{0}-#Delta#theta","f");
        }
    } else {
        leg2->AddEntry(&g1,"#theta = #hat{#theta}+#Delta#hat{#theta}","f");
        leg2->AddEntry(&g2,"#theta = #hat{#theta}-#Delta#hat{#theta}","f");
    }
    leg2->Draw();

    TLegend leg0(0.02,0.1,1,0.25);
    leg0.SetFillStyle(0);
    leg0.SetBorderSize(0);
    leg0.SetMargin(0.2);
    leg0.SetTextFont(gStyle->GetTextFont());
    leg0.SetTextSize(gStyle->GetTextSize());
    leg0.AddEntry(&g,"Nuis. Param. Pull","lp");
    leg0.Draw();

    c.cd();

    TLine l0(0,- offsetDown1/lineHeight,0,SIZE+0.5);// + offsetUp1/lineHeight);
    l0.SetLineStyle(kDashed);
    l0.SetLineColor(kBlack);
    l0.Draw("same");
    TLine l1 (-1,- offsetDown1/lineHeight,-1,SIZE+0.5);// + offsetUp1/lineHeight);
    l1.SetLineStyle(kDashed);
    l1.SetLineColor(kBlack);
    l1.Draw("same");
    TLine l2(1,- offsetDown1/lineHeight,1,SIZE+0.5);// + offsetUp1/lineHeight);
    l2.SetLineStyle(kDashed);
    l2.SetLineColor(kBlack);
    l2.Draw("same");

    for (auto& imarker : markers) {
        imarker->SetMarkerSize(1.2*imarker->GetMarkerSize());
        imarker->Draw("same");
    }

    if (fPlotLabel!= "none") TRExLabelNew(0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.6*offsetUp1)/newHeight), TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str(), kBlack, gStyle->GetTextSize());
    myText(0.42,(1.*(offsetDown+offsetDown1+SIZE*lineHeight+0.3*offsetUp1)/newHeight), 1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    gPad->RedrawAxis();

    if(flagGammas && flagSysts){
        Common::SaveCanvasAs(c, fName+"/Rankings/Ranking"+fSuffix);
    } else if(flagGammas){
        Common::SaveCanvasAs(c, fName+"/Rankings/RankingGammas"+fSuffix);
    } else if(flagSysts){
        Common::SaveCanvasAs(c, fName+"/Rankings/RankingSysts"+fSuffix);
    } else{
        LOG(WARNING) << "Your ranking plot fell in unknown category:\n";
        Common::SaveCanvasAs(c, fName+"/Rankings/RankingUnknown"+fSuffix);
    }
}

//__________________________________________________________________________________
//
void RankingManager::DumpCovariances(const std::vector<std::pair<std::string,double> >& pois,
                                     const std::vector<std::string>& nfs,
                                     const std::string& folder,
                                     const FitResults* fr,
                                     const std::map<std::string, double>& statOnly,
                                     const bool isUnfolding) const {

    gSystem->mkdir((fName+"/Covariances").c_str());

    // first, read the ranking txt files
    std::vector<std::string> params;
    std::vector<std::map<std::string, YamlConverter::RankingContainer> > rankingVector = this->GetRankingImpacts(folder, pois, nfs, params);

    std::vector<std::string> names;
    if (isUnfolding) {
        for (const auto& ipoi : pois) {
            const std::string betterName = Common::ReplaceString(ipoi.first, "_mu", "_yield");
            names.emplace_back(betterName);
        }
    } else {
        for (const auto& ipoi : pois) {
            names.emplace_back(ipoi.first);
        }
    }

    // for the total covariance calculation
    TMatrixD totalPre(pois.size(), pois.size());
    TMatrixD totalPost(pois.size(), pois.size());

    // initialize to zero
    for (std::size_t i = 0; i < pois.size(); ++i) {
        for (std::size_t j = 0; j < pois.size(); ++j) {
            totalPre(i,j) = 0.;
            totalPost(i,j) = 0.;
        }
    }

    // loop over parameters and fill the covariances
    for (const auto& ipar : params) {
        TMatrixD covPost(pois.size(), pois.size());
        TMatrixD covPre(pois.size(), pois.size());

        for (std::size_t i = 0; i < rankingVector.size(); ++i) {
            auto itr_i = rankingVector.at(i).find(ipar);
            if (itr_i == rankingVector.at(i).end()) {
                LOG(ERROR) << "Cannot find parameter " << ipar << " in the ranking txt file number " << i << "\n";
                return;
            }
            const auto& container_i = itr_i->second;
            const double deltaPost_i = this->Symmetrize(container_i.name, container_i.poihi, container_i.poilo);
            const double deltaPre_i = this->Symmetrize(container_i.name, container_i.poiprehi, container_i.poiprelo);
            for (std::size_t j = 0; j < rankingVector.size(); ++j) {
                auto itr_j = rankingVector.at(j).find(ipar);
                if (itr_j == rankingVector.at(j).end()) {
                    LOG(ERROR) << "Cannot find parameter " << ipar << " in the ranking txt file number " << j << "\n";
                    return;
                }
                const auto& container_j = itr_j->second;
                const double deltaPost_j = this->Symmetrize(container_j.name, container_j.poihi, container_j.poilo);
                const double deltaPre_j = this->Symmetrize(container_j.name, container_j.poiprehi, container_j.poiprelo);
                covPost(i,j) = deltaPost_i * deltaPost_j;
                covPre(i,j)  = deltaPre_i * deltaPre_j;
            }
        }


        // print the covariances
        YamlConverter converter;
        converter.WriteCorrelation(names, covPost, fName, true, "Covariances/Covariance_postfit_"+ipar);
        converter.WriteCorrelation(names, covPre, fName, true, "Covariances/Covariance_prefit_"+ipar);
        if (fHEPDataFormat) {
            converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
            converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
            converter.WriteCorrelationHEPData(names, covPost, fName, true, "Covariance_postfit_"+ipar);
            converter.WriteCorrelationHEPData(names, covPre, fName, true, "Covariance_prefit_"+ipar);
        }

        totalPre  += covPre;
        totalPost += covPost;
    }

    // add the stat only if it exists
    if (fr) {
        TMatrixD stat(pois.size(), pois.size());
        for (std::size_t ipoi = 0; ipoi < pois.size(); ++ipoi) {
            auto itr_i = statOnly.find(pois.at(ipoi).first);
            if (itr_i == statOnly.end()) {
                LOG(ERROR) << "Cannot find parameter: " << pois.at(ipoi).first << "\n";
                continue;
            }
            const double err_i = itr_i->second;
            for (std::size_t jpoi = 0; jpoi < pois.size(); ++jpoi) {
                auto itr_j = statOnly.find(pois.at(jpoi).first);
                if (itr_j == statOnly.end()) {
                    LOG(ERROR) << "Cannot find parameter: " << pois.at(jpoi).first << "\n";
                    continue;
                }
                const double err_j = itr_j->second;
                const double corr = (ipoi == jpoi) ? 1.0 : fr->GetCorrelationMatrix()->GetCorrelation(pois.at(ipoi).first, pois.at(jpoi).first);
                const double cov = err_i * corr * err_j;
                stat(ipoi,jpoi) = cov;
            }
        }

        totalPre  += stat;
        totalPost += stat;

        // print the covariances
        YamlConverter converter;
        converter.WriteCorrelation(names, stat, fName, true, "Covariances/Covariance_statOnly");
        if (fHEPDataFormat) {
            converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
            converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
            converter.WriteCorrelationHEPData(names, stat, fName, true, "Covariance_statOnly");
        }
    }

    this->PlotTotalCovariance(totalPre, names, true);
    this->PlotTotalCovariance(totalPost, names, false);

    YamlConverter converter;
    converter.WriteCorrelation(names, totalPre,  fName, true, "Covariances/Combined_Covariance_prefit");
    converter.WriteCorrelation(names, totalPost, fName, true, "Covariances/Combined_Covariance_postfit");
    if (fHEPDataFormat) {
        converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
        converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));

        converter.WriteCorrelationHEPData(names, totalPost, fName, true, "Combined_Covariance_postfit");
        converter.WriteCorrelationHEPData(names, totalPre, fName, true, "Combined_Covariance_prefit");
        // append the submission script
        converter.AppendSubmissionWithCovariances(params, fName + "/HEPData/submission.yaml");
    }
}

//__________________________________________________________________________________
//
void RankingManager::PlotTotalCovariance(const TMatrixD& matrix, const std::vector<std::string>& names, const bool isPrefit) const {

    const int n = matrix.GetNrows();

    TH2D h("", "", n, 0.5, n+0.5, n, 0.5, n+0.5);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            h.SetBinContent(i+1, n-j, matrix(i,j));
        }
        h.GetXaxis()->SetBinLabel(i+1, names.at(i).c_str());
        h.GetYaxis()->SetBinLabel(n-i, names.at(i).c_str());
    }

    TCanvas c1("","",0.,0.,1200,1200);
    gStyle->SetPalette(87);
    h.SetMarkerSize(0.75*1000);
    gStyle->SetPaintTextFormat("1.1e");
    c1.SetBottomMargin(0.2);
    c1.SetLeftMargin(0.2);

    h.GetXaxis()->LabelsOption("v");
    h.GetXaxis()->SetLabelSize( h.GetXaxis()->GetLabelSize()*0.75 );
    h.GetYaxis()->SetLabelSize( h.GetYaxis()->GetLabelSize()*0.75 );
    c1.SetTickx(0);
    c1.SetTicky(0);
    h.GetYaxis()->SetTickLength(0);
    h.GetXaxis()->SetTickLength(0);
    h.Draw("col TEXT");

    Common::SaveCanvasAs(c1, fName + "/Covariances/Combined_CovarianceMatrix_" + (isPrefit ? "preFit" : "postfit"));
}

//__________________________________________________________________________________
//
double RankingManager::Symmetrize(const std::string& param, const double up, const double down) const {
    // has the same sign
    if (up * down > 0) {
        LOG(WARNING) << "Values have the same shift for parameter: " << param << ", up: " << up << " down: " << down << "\n";
        return up;
    }

    return 0.5*(up - down);
}

//__________________________________________________________________________________
//
std::vector<std::map<std::string, YamlConverter::RankingContainer> > RankingManager::GetRankingImpacts(const std::string& folder,
                                                                                                       const std::vector<std::pair<std::string, double> >& pois,
                                                                                                       const std::vector<std::string>& nfs,
                                                                                                       std::vector<std::string>& params) const {

    std::vector<std::map<std::string, YamlConverter::RankingContainer> > result;
    for (const auto& ipoi : pois) {
        std::string param;
        double nuiphat;
        double nuiperrhi;
        double nuiperrlo;
        double PoiUp;
        double PoiDown;
        double PoiNomUp;
        double PoiNomDown;
        std::ifstream file(folder + "/Fits/NPRanking_"+ipoi.first+".txt");
        if (!file.is_open() || !file.good()) {
            LOG(ERROR) << "Cannot open file " << folder << "/Fits/NPRanking_" << ipoi.first << ".txt\n";
            return result;
        }

        std::map<std::string, YamlConverter::RankingContainer> rankingMap;
        YamlConverter::RankingContainer container;
        while (file >> param >> nuiphat >> nuiperrhi >> nuiperrlo >> PoiUp >> PoiDown >> PoiNomUp >> PoiNomDown) {
            if (std::find(nfs.begin(), nfs.end(), param) != nfs.end()) continue;
            container.name = param;
            container.nphat = nuiphat;
            container.nperrhi = nuiperrhi;
            container.nperrlo = nuiperrlo;
            container.poihi = PoiUp * ipoi.second;
            container.poilo = PoiDown * ipoi.second;
            container.poiprehi = PoiNomUp * ipoi.second;
            container.poiprelo = PoiNomDown * ipoi.second;

            auto itr = rankingMap.find(param);
            if (itr != rankingMap.end()) {
                LOG(WARNING) << "Parameter " << param << " already in the list\n";
                continue;
            }

            rankingMap.insert(std::make_pair(param, container));

            if (result.empty()) {
                if (std::find(params.begin(), params.end(), param) == params.end()) {
                    params.emplace_back(param);
                }
            }
        }
        result.emplace_back(std::move(rankingMap));
        file.close();
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<double> RankingManager::ImpactFromCorrelations(const CorrelationMatrix* corr,
                                                           const std::vector<std::pair<std::string, double> >& totalUncertainty,
                                                           const std::pair<std::string, bool>& np,
                                                           const RankingManager::RankingValues& values,
                                                           const bool isPrefit) const {

    if (isPrefit && np.second) {
        std::vector<double> tmp(totalUncertainty.size(), 0);
        return tmp;
    }

    std::vector<double> result(totalUncertainty.size());

    for (std::size_t ipoi = 0; ipoi < totalUncertainty.size(); ++ipoi) {
        std::string npNameSanitised = Common::ReplaceString(np.first, "alpha_", "");
        npNameSanitised = Common::ReplaceString(npNameSanitised, "gamma_", "");
        const double correlation = corr->GetCorrelation(npNameSanitised, totalUncertainty.at(ipoi).first);
        double uncertainty(0);
        if (isPrefit) {
            if (np.first.find("gamma_") != std::string::npos) {
                uncertainty = 0.5*(values.upPrefit - values.downPrefit);
            } else {
                uncertainty = 1.;
            }
        } else {
            if (np.first.find("gamma_") != std::string::npos) {
                uncertainty = (values.up - values.down)/(values.upPrefit - values.downPrefit);
            } else {
                if (np.second) {
                    uncertainty = 1;
                } else {
                    uncertainty = 0.5*(values.up - values.down);
                }
            }
        }
        result.at(ipoi) = correlation * uncertainty * totalUncertainty.at(ipoi).second;
        LOG(DEBUG) << "NP: " << np.first << ", correlation: " << correlation << ", uncertainty: " << uncertainty << ", totalUncertainty: " << totalUncertainty.at(ipoi).second << ", impact: " << result.at(ipoi) << "\n";
    }

    return result;
}

//__________________________________________________________________________________
//
void RankingManager::ReadRankingResults(const bool flagSysts, const bool flagGammas) {

    fRankingContainer.clear();

    std::string paramname{};
    double nphat;
    double nperrhi;
    double nperrlo;
    double poihi;
    double poilo;
    double poiprefithi;
    double poiprefitlo;

    std::ifstream fin(fOutputPath.c_str());
    const std::string temp_string = "Systematic called \"Luminosity\" found. This creates issues for the ranking plot. Skipping. Suggestion: rename this systematic as \"Lumi\" or \"luminosity\"";
    while (fin >> paramname >> nphat >> nperrhi >> nperrlo >> poihi >> poilo >> poiprefithi >> poiprefitlo){
        if (paramname == "Luminosity"){
            LOG(ERROR) << temp_string << "\n";
            return;
        }
        if(paramname.find("gamma")!=std::string::npos && !flagGammas){
            fin >> paramname >> nphat >> nperrhi >> nperrlo >> poihi >> poilo >> poiprefithi >> poiprefitlo;
            continue;
        }
        if(paramname.find("gamma")==std::string::npos && !flagSysts){
            fin >> paramname >> nphat >> nperrhi >> nperrlo >> poihi >> poilo >> poiprefithi >> poiprefitlo;
            continue;
        }

        YamlConverter::RankingContainer container;
        container.name = paramname;
        container.nphat = nphat;
        container.nperrhi = nperrhi;
        container.nperrlo = nperrlo;
        container.poihi = poihi;
        container.poilo = poilo;
        container.poiprehi = poiprefithi;
        container.poiprelo = poiprefitlo;

        fRankingContainer.emplace_back(std::move(container));
    }
}

//__________________________________________________________________________________
//
void RankingManager::ProduceGroupedImpact(const std::map<std::string, std::string>& categoryMap) const {

    const std::string outName = fName+"/GroupedImpactFromRanking" + fSuffix + ".txt";
    std::ofstream output(outName);
    if (!output.good() || !output.is_open()) {
        LOG(ERROR) << "Cannot open file: " << outName << "\n";
        return;
    }

    LOG(INFO) << "Creating grouped impact file at: " << outName << "\n";
    this->ProcessSingleGroupedImpact(&output, categoryMap);

    output.close();
}

//__________________________________________________________________________________
//
void RankingManager::ProcessSingleGroupedImpact(std::ofstream* file,
                                                const std::map<std::string, std::string>& categoryMap) const {

    std::map<std::string, std::tuple<double, double, double> > uncertainty;

    std::map<std::string, std::string> map = categoryMap;
    for (const auto& i : fRankingContainer) {
        if (Common::StartsWith(i.name, "gamma_")) {
            map.insert({i.name, "Gammas"});
        }
    }

    double total2(0);
    double totalUp2(0);
    double totalDown2(0);

    for (const auto& category : map) {
        if (category.first == "DUMMY_STATONLY") continue;
        if (category.first == "DUMMY_GAMMAS") continue;
        auto itrParam = std::find_if(fRankingContainer.begin(), fRankingContainer.end(), [&category](const auto& element) {return category.first.find(element.name) != std::string::npos;});
        if (itrParam == fRankingContainer.end()) {
            LOG(DEBUG) << "Cannot find parameter: " << category.first << "\n";
            continue;
        }

        const auto& ranking = (*itrParam);

        const double impactUp2   = ranking.poihi*ranking.poihi;
        const double impactDown2 = ranking.poilo*ranking.poilo;
        double impact2 = 0.5*(std::abs(ranking.poihi - ranking.poilo));
        impact2 = impact2 * impact2;

        total2     += impact2;
        totalUp2   += impactUp2;
        totalDown2 += impactDown2;

        auto itrCategory = uncertainty.find(category.first);
        if (itrCategory == uncertainty.end()) {
            std::tuple<double, double, double> tmp = std::make_tuple(impact2, impactUp2, impactDown2);
            uncertainty.insert({category.second, tmp});
        } else {
            double currentImpact2     = std::get<0>(itrCategory->second);
            double currentImpactUp2   = std::get<1>(itrCategory->second);
            double currentImpactDown2 = std::get<2>(itrCategory->second);
            std::tuple<double, double, double> tmp = std::make_tuple(impact2+currentImpact2, impactUp2 + currentImpactUp2, impactDown2 + currentImpactDown2);
            itrCategory->second = tmp;
        }
    }

    for (const auto& unc : uncertainty) {
        *file << unc.first << "    " << std::sqrt(std::get<0>(unc.second)) << "  ( +" << std::sqrt(std::get<1>(unc.second)) << ", -" << std::sqrt(std::get<2>(unc.second)) << " )\n";
    }
    *file << "FullSyst" << "    " << std::sqrt(total2) << "  ( +" << std::sqrt(totalUp2) << ", -" << std::sqrt(totalDown2) << " )\n";
}
