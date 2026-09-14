// Class include
#include "TRExFitter/Region.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/YamlConverter.h"

// ROOT includes
#include "Math/DistFunc.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TMatrixD.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "xRooFit/xRooFit.h"

// c++ includes
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

// -------------------------------------------------------------------------------------------------
// class Region

//__________________________________________________________________________________
//
Region::Region(const string& name) :
    fName(name),
    fUseFriend(false),
    fVariableTitle(""),
    fYTitle(""),
    fLabel(name),
    fShortLabel(name),
    fTexLabel(""),
    fFitName(""),
    fPlotSubdir(""),
    fRegionType(CONTROL),
    fRegionDataType(REALDATA),
    fHasData(false),
    fData(nullptr),
    fHasSig(false),
    fYmaxScale(0),
    fYmin(0),
    fYmax(0),
    fRatioYmin(0),
    fRatioYmax(2.),
    fRatioYminPostFit(0.5),
    fRatioYmaxPostFit(1.5),
    fRatioYtitle(""),
    fRatioType(TRExPlot::RATIOTYPE::DATAOVERMC),
    fTot(nullptr),
    fErr(nullptr),
    fTot_postFit(nullptr),
    fBkgOnly_postFit(nullptr),
    fErr_postFit(nullptr),
    fErrBkgOnly_postFit(nullptr),
    fBinTransfo(""),
    fTransfoDzBkg(0.),
    fTransfoDzSig(0.),
    fTransfoErr(10),
    fTransfoFzBkg(0.),
    fTransfoFzSig(0.),
    fTransfoJpar1(0.),
    fTransfoJpar2(0.),
    fTransfoJpar3(0.),
    fVariable(""),
    fCorrVar1(""),
    fCorrVar2(""),
    fXmin(0),
    fXmax(0),
    fSelection("1"),
    fMCweight("1"),
    fHistoNBinsRebin(-1),
    fHistoNBinsRebinPost(-1),
    fUseStatErr(false),
    fIntCode_overall(4),
    fIntCode_shape(0),
    fFitType(TRExFit::SPLUSB),
    fFitLabel(""),
    fLumiLabel(""),
    fCmeLabel(""),
    fLumiScale(1.),
    fLogScale(false),
    fLogScaleX(false),
    fBinWidth(0),
    fBlindingThreshold(-1),
    fBlindingType(Common::SOVERB),
    fSkipSmoothing(false),
    fPlotLabel("Internal"),
    fSuffix(""),
    fGroup(""),
    fKeepPrefitBlindedBins(false),
    fGetChi2(0),
    fChi2val(-1),
    fNDF(-1),
    fChi2prob(-1),
    fUseGammaPulls(false),
    fLabelX(-1),
    fLabelY(-1),
    fLegendX1(-1),
    fLegendX2(-1),
    fLegendY(-1),
    fLegendNColumns(2),
    fNumberUnfoldingRecoBins(-1),
    fNormalizeMigrationMatrix(true),
    fHasAcceptance(false),
    fFolder(""),
    fHEPDataFormat(false),
    fUheppFormat(false),
    fDrawDataMinusBkg(false),
    fDrawDataMinusBkgOnSignal(true),
    fJobName(""),
    fJobBinningsPath(""),
    fNoPrePostFitLabel(false),
    fPreFitLabel("Pre-fit"),
    fPostFitLabel("Post-fit"),
    fToysForErrorBand(-1),
    fAutomaticDropBins(false),
    fNbins(0)
{
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> Region::SetSampleHist(Sample *sample,
                                                  const std::string& histoName,
                                                  const std::string& fileName,
                                                  const bool& dropBinHistoNeeded,
                                                  const bool& checkHistogram,
                                                  const Common::FolderStructure* structure) {
    fSampleHists.emplace_back(new SampleHist( sample, histoName, fileName, dropBinHistoNeeded, checkHistogram, structure));
    if(sample->fType==Sample::SampleType::DATA){
        fHasData = true;
        fData = fSampleHists.back();
    }
    else if(sample->fType==Sample::SampleType::SIGNAL){
        fHasSig = true;
        fSig.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::SampleType::BACKGROUND){
        fBkg.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::SampleType::GHOST){
        LOG(DEBUG) << "Adding GHOST sample.\n";
    }
    else if(sample->fType==Sample::SampleType::EFT){
        LOG(DEBUG) << "Adding EFT sample.\n";
    }
    else{
        LOG(ERROR) << "SampleType not supported.\n";
    }
    fSampleHists.back()->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists.back()->fRegionName = fName;
    fSampleHists.back()->fRegionLabel = fLabel;
    fSampleHists.back()->fFitName = fFitName;
    fSampleHists.back()->fVariableTitle = fVariableTitle;
    fSampleHists.back()->fPlotLabel = fPlotLabel;
    return fSampleHists.back();
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> Region::SetSampleHist(Sample *sample, TH1* hist ){
    fSampleHists.emplace_back(new SampleHist( sample, hist ));
    if(sample->fType==Sample::SampleType::DATA){
        fHasData = true;
        fData = fSampleHists.back();
    }
    else if(sample->fType==Sample::SampleType::SIGNAL){
        fHasSig = true;
        fSig.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::SampleType::BACKGROUND){
        fBkg.emplace_back(fSampleHists.back());
    }
    else if(sample->fType==Sample::SampleType::GHOST){
        LOG(DEBUG) << "Adding GHOST sample.\n";
    }
    else if(sample->fType==Sample::SampleType::EFT){
        LOG(DEBUG) << "Adding EFT sample.\n";
    }
    else{
        LOG(ERROR) << "SampleType not supported.\n";
    }
    fSampleHists.back()->fHist->SetName(Form("%s_%s",fName.c_str(),sample->fName.c_str()));
    fSampleHists.back()->fRegionName = fName;
    fSampleHists.back()->fRegionLabel = fLabel;
    fSampleHists.back()->fFitName = fFitName;
    fSampleHists.back()->fVariableTitle = fVariableTitle;
    fSampleHists.back()->fPlotLabel = fPlotLabel;
    return fSampleHists.back();
}

//__________________________________________________________________________________
//
void Region::AddSample(Sample* sample){
    fSamples.emplace_back(std::move(sample));
}

//__________________________________________________________________________________
//
void Region::SetBinning(int N, const std::vector<double>& bins){
    fNbins = fHistoNBinsRebin = N;
    fHistoBins = bins;
    this->AddBinsToFile();
}

//__________________________________________________________________________________
//
void Region::Rebin(int N){
    fNbins = N;
    fHistoNBinsRebin = N;
    this->AddBinsToFile();
}

//__________________________________________________________________________________
//
void Region::SetRebinning(int N, const std::vector<double>& bins){
    fHistoNBinsRebinPost = N;
    fHistoBinsPost = bins;
}

//__________________________________________________________________________________
//
void Region::SetRegionType( RegionType type ){
    fRegionType = type;
}

//__________________________________________________________________________________
//
void Region::SetRegionDataType( DataType type ){
    fRegionDataType = type;
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> Region::GetSampleHist(const std::string &sampleName) const{
    for(const auto& isample : fSampleHists) {
        if(isample->fName == sampleName) return isample;
    }
    return nullptr;
}

//__________________________________________________________________________________
//
std::shared_ptr<TRExPlot> Region::DrawPrePostFit(const std::shared_ptr<xRooNode>& pdf,
                                                 const std::shared_ptr<xRooNode>& pdfPrefit,
                                                 const FitResults* fitRes,
                                                 const std::vector<int>& canvasSize,
                                                 const bool isPostFit,
                                                 const bool isBonly,
                                                 std::string opt) {

    int canvasWidth = 600;
    int canvasHeight = 700;
    const std::string cName = "c_"+fName;
    if(TRExFitter::OPTION["CanvasWidth"]!=0)  canvasWidth  = TRExFitter::OPTION["CanvasWidth"];
    if(TRExFitter::OPTION["CanvasHeight"]!=0) canvasHeight = TRExFitter::OPTION["CanvasHeight"];
    if (!canvasSize.empty()) {
        canvasWidth = canvasSize.at(0);
        canvasHeight = canvasSize.at(1);
    }

    std::shared_ptr<TRExPlot> plot = std::make_shared<TRExPlot>(cName,canvasWidth,canvasHeight,TRExFitter::NORATIO);

    if(isPostFit && TRExFitter::PREFITONPOSTFIT){
        if (!pdfPrefit) {
            LOG(ERROR) << "No prefit pdf provided for PREFITONPOSTFIT\n";
            exit(EXIT_FAILURE);
        }

        auto channel = pdfPrefit->find(fName);
        if (!channel) {
            LOG(ERROR) << "Cannot read channel: " << fName << "\n";
            exit(EXIT_FAILURE);
        }
        const std::string bkgSamples = Common::NonSignalSampleList(fSampleHists, fName);
        xRooNode bkgOnlySample = FitUtils::GetReducedOrSingleSample(channel, bkgSamples);

        auto histTotal = Common::HistoFromxRooNode(bkgOnlySample, fSampleHists.at(0)->fHist.get());

        plot->h_tot_bkg_prefit = static_cast<TH1*>(histTotal->Clone("h_tot_bkg_prefit"));
        plot->h_tot_bkg_prefit->SetDirectory(nullptr);
    }

    plot->fShowYields = TRExFitter::SHOWYIELDS;

    plot->SetXaxisRange(fXaxisRange);
    if(fYmaxScale==0) plot->SetYmaxScale(1.8);
    else              plot->SetYmaxScale(fYmaxScale);
    if(fYmax!=0) plot->fYmax = fYmax;
    if(fYmin!=0) plot->fYmin = fYmin;
    plot->fRatioYmax = isPostFit ? fRatioYmaxPostFit : fRatioYmax;
    plot->fRatioYmin = isPostFit ? fRatioYminPostFit : fRatioYmin;
    plot->SetXaxis(fVariableTitle,fVariableTitle.find("Number")!=string::npos);
    if(fYTitle!="") plot->SetYaxis(fYTitle);
    if (fFitLabel != "none") plot->AddLabel(fFitLabel);
    if (fLabel != "none") plot->AddLabel(fLabel);
    if (isPostFit) {
        if (!fNoPrePostFitLabel) {
            plot->AddLabel(fPostFitLabel);
        }
    } else {
        if (!fNoPrePostFitLabel) {
            plot->AddLabel(fPreFitLabel);
        }
    }
    plot->SetLumi(fLumiLabel);
    plot->SetCME(fCmeLabel);
    plot->SetLumiScale(fLumiScale);
    plot->fLegendNColumns = fLegendNColumns;

    if (!fBinLabels.empty()) {
        plot->ResizeBinLabel(fBinLabels.size() + 1);
        for(std::size_t i_bin=0; i_bin<fBinLabels.size(); i_bin++) {
            plot->SetBinLabel(i_bin+1,fBinLabels.at(i_bin));
        }
    }

    plot->fBinDividers = fBinDividers;
    plot->fRangeLabels = fRangeLabels;
    plot->fPreFitLabel = fPreFitLabel;

    YamlConverter::PlotContainer container;
    container.region = fName;
    container.xAxis = fVariableTitle;
    container.yAxis = fYTitle;

    auto channel = pdf->find(fName);
    if (!channel) {
        LOG(ERROR) << "Cannot read channel: " << fName << "\n";
        exit(EXIT_FAILURE);
    }

    // randomize NFs if set
    if (!isPostFit && !m_randomNFmap.empty()) {
        TRandom3 rng(12345);
        for (const auto& iparam : pdf->pars()) {
            auto itr = m_randomNFmap.find(iparam->GetName());
            if (itr != m_randomNFmap.end()) {
                const double value = rng.Uniform(itr->second.first, itr->second.second);
                channel->pars()[iparam->GetName()]->get<RooRealVar>()->setVal(value);
            }
        }
    }

    std::vector<std::shared_ptr<TH1> > allHistos;
    std::vector<std::shared_ptr<TH1> > bOnlyNormSigHistos;

    auto samples = channel->find("samples");
    if (!samples) {
        LOG(ERROR) << "Cannot find samples\n";
        exit(EXIT_FAILURE);
    }
    for (const auto& isample : fSampleHists) {
        const bool skip = isample->fSample->fType==Sample::SampleType::DATA || isample->fSample->fType==Sample::SampleType::GHOST
                    || isample->fSample->fType==Sample::SampleType::EFT;
        if (skip) {
            allHistos.emplace_back(static_cast<TH1*>(isample->fHist->Clone()));
            allHistos.back()->SetDirectory(nullptr);
            bOnlyNormSigHistos.emplace_back(static_cast<TH1*>(isample->fHist->Clone()));
            bOnlyNormSigHistos.back()->SetDirectory(nullptr);
            if (isPostFit) {
                isample->fHist_postFit.reset(static_cast<TH1*>(isample->fHist->Clone()));
                isample->fHist_postFit->SetDirectory(nullptr);
            } else {
                isample->fHist.reset(static_cast<TH1*>(isample->fHist->Clone()));
                isample->fHist->SetDirectory(nullptr);
            }
            continue;
        }
        auto sample = samples->find(isample->fName + "_" + fName + "_shapes");
        if (!sample) {
            LOG(ERROR) << "Cannot find sample: " << isample->fName << "_" << fName << "_shapes\n";
            exit(EXIT_FAILURE);
        }

        const bool bOnlyNorm = isample->fSample->fType == Sample::SampleType::SIGNAL && isPostFit && TRExFitter::SHOWNORMSIG && isBonly;
        std::shared_ptr<TH1> bOnlyNormSig(nullptr);
        // for signal in case of BONLY, postfit, NORMSIG we need to set its value by hand to 1 to be able to normalise it
        if (bOnlyNorm) {
            auto params = channel->pars();
            std::vector<std::string> poiNames;
            for (const auto& inf : isample->fSample->fNormFactors) {
                auto itr = std::find(fPOIs.begin(), fPOIs.end(), inf->fName);
                if (itr != fPOIs.end()) {
                    poiNames.emplace_back(inf->fName);
                }
            }
            for (const auto& ipoi : poiNames) {
                auto poi = params.find(ipoi);
                if (!poi) {
                    LOG(ERROR) << "Cannot find parameter: " << ipoi << ", in region: " << fName << "\n";
                    exit(EXIT_FAILURE);
                }
                auto poiParam = poi->get<RooRealVar>();
                poiParam->setVal(1);
            }
            bOnlyNormSig = Common::HistoFromxRooNode(sample, isample->fHist.get());

            // reset back pois to 0
            for (const auto& ipoi : poiNames) {
                auto poi = params.find(ipoi);
                auto poiParam = poi->get<RooRealVar>();
                poiParam->setVal(0);
            }
        }
        auto h = Common::HistoFromxRooNode(sample, isample->fHist.get());
        if (isPostFit) {
            isample->fHist_postFit.reset(static_cast<TH1*>(h->Clone()));
            isample->fHist_postFit->SetDirectory(nullptr);
        } else {
            isample->fHist.reset(static_cast<TH1*>(h->Clone()));
            isample->fHist->SetDirectory(nullptr);
        }

        if (isPostFit) {
            fPostfitYields.insert(std::make_pair(isample->fName, isample->fHist_postFit->Integral()));
        }
        if (bOnlyNorm) {
            fPostfitYieldsBonlyNorm.insert(std::make_pair(isample->fName, bOnlyNormSig->Integral()));
        }
        if (!isPostFit || TRExFitter::PREFITONPOSTFIT) {
            fPrefitYields.insert(std::make_pair(isample->fName, isample->fHist->Integral()));
        }
        allHistos.emplace_back(h);
        bOnlyNormSigHistos.emplace_back(bOnlyNormSig);
    }

    auto StoreYields = [this](std::vector<std::vector<double> >& containerYields, const std::vector<std::shared_ptr<TH1> >& hist, const Sample::SampleType type) {
        std::vector<std::string> unique;
        for(std::size_t i = 0; i < this->fSampleHists.size(); ++i){
            if (this->fSampleHists.at(i)->fSample->fType != type) continue;
            if (!hist.at(i)) continue;
            std::string title = this->fSampleHists.at(i)->fSample->fTitle;
            if(this->fSampleHists.at(i)->fSample->fGroup != "") title = this->fSampleHists.at(i)->fSample->fGroup;
            std::vector<double> tmp;
            for (int ibin = 1; ibin <= hist.at(i)->GetNbinsX(); ++ibin) {
                tmp.emplace_back(hist.at(i)->GetBinContent(ibin));
            }
            auto itr = std::find(unique.begin(), unique.end(), title);
            if (itr == unique.end()) {
                containerYields.emplace_back(std::move(tmp));
                unique.emplace_back(title);
            } else {
                const std::size_t pos = std::distance(unique.begin(), itr);
                for (std::size_t ii = 0; ii < containerYields.at(pos).size(); ++ii) {
                    containerYields.at(pos).at(ii) += tmp.at(ii);
                }
            }
        }
    };

    StoreYields(container.signalYields, allHistos, Sample::SampleType::SIGNAL);
    StoreYields(container.backgroundYields, allHistos, Sample::SampleType::BACKGROUND);

    auto AddSamplesToContainer = [this](std::vector<std::string>& containerSamples, const Sample::SampleType type) {
        for (const auto& isample : this->fSampleHists) {
            const auto sampleType = isample->fSample->fType;
            if (sampleType != type) continue;
            std::string title = isample->fSample->fTitle;
            if (isample->fSample->fGroup != "") title = isample->fSample->fGroup;
            auto itr = std::find(containerSamples.begin(), containerSamples.end(), title);
            if (itr != containerSamples.end()) continue;
            containerSamples.emplace_back(title);
        }
    };

    // add first the signals and then background
    AddSamplesToContainer(container.samples, Sample::SampleType::SIGNAL);
    AddSamplesToContainer(container.samples, Sample::SampleType::BACKGROUND);

    for (std::size_t i = 0; i < fSampleHists.size(); ++i) {
        std::string title = fSampleHists.at(i)->fSample->fTitle;
        const bool bOnlyNorm = fSampleHists.at(i)->fSample->fType == Sample::SampleType::SIGNAL && isPostFit && TRExFitter::SHOWNORMSIG && isBonly;
        if (fSampleHists.at(i)->fSample->fGroup != "") title = fSampleHists.at(i)->fSample->fGroup;
        if (fSampleHists.at(i)->fSample->fType == Sample::SampleType::SIGNAL) {
            if(TRExFitter::SHOWSTACKSIG) plot->AddSignal(allHistos.at(i).get(),title);
            if(TRExFitter::SHOWNORMSIG){
                if( (TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL)
                 || !TRExFitter::OPTION["NormSigSRonly"]) {
                    if (bOnlyNorm) {
                        plot->AddNormSignal(bOnlyNormSigHistos.at(i).get(),title);
                    } else {
                        plot->AddNormSignal(allHistos.at(i).get(),title);
                    }
                }
            } else {
                if (TRExFitter::OPTION["NormSigSRonly"] && fRegionType==SIGNAL) {
                    plot->AddNormSignal(allHistos.at(i).get(),title);
                }
            }
            if (TRExFitter::SHOWOVERLAYSIG) {
                if (TRExFitter::SHOWOVERLAYSIG_CUSTOMSCALE < 0) { // if no custom scale is set (default: -1), scale whatever is prediction or postfit
                    plot->AddOverSignal(allHistos.at(i).get(),title);
                } else {
                    if (isPostFit) {
                        double normSum = 0.;
                        for(const auto& nf: fSampleHists.at(i)->fSample->fNormFactors){
                            normSum += fitRes->GetNuisParValue(nf->fName);
                        }
                        double poiScale(0);
                        if (normSum != 0.) {
                            poiScale = 1./normSum * TRExFitter::SHOWOVERLAYSIG_CUSTOMSCALE;
                        } else {
                            poiScale = TRExFitter::SHOWOVERLAYSIG_CUSTOMSCALE;
                        }
                        plot->AddOverSignal((allHistos.at(i).get()), title, poiScale);
                    } else {
                        plot->AddOverSignal((allHistos.at(i).get()), title, TRExFitter::SHOWOVERLAYSIG_CUSTOMSCALE);
                    }
                }
            }
        } else if (fSampleHists.at(i)->fSample->fType == Sample::SampleType::BACKGROUND) {
            plot->AddBackground(allHistos.at(i).get(),title);
        }
    }

    if (fHasData && opt.find("blind") == string::npos) {
        for (int ibin = 1; ibin <= fData->fHist->GetNbinsX(); ++ibin) {
            container.data.emplace_back(fData->fHist->GetBinContent(ibin));
        }
        plot->SetData(fData->fHist.get(),fData->fSample->fTitle);
    }

    for(std::size_t i = 0; i < fSampleHists.size(); ++i) {
        if(fSampleHists.at(i)->fSample->fType==Sample::SampleType::DATA) continue;
        if(fSampleHists.at(i)->fSample->fType==Sample::SampleType::GHOST) continue;
        if(fSampleHists.at(i)->fSample->fType==Sample::SampleType::EFT) continue;
        if(fSampleHists.at(i)->fSample->fType==Sample::SampleType::SIGNAL && !(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG)) continue;
        if (isPostFit) {
            if(!fTot_postFit) {
                fTot_postFit.reset(static_cast<TH1*>(allHistos.at(i)->Clone("h_tot_postFit")));
            } else {
                fTot_postFit->Add(allHistos.at(i).get());
            }
            // add BkgOnly
            if (fSampleHists.at(i)->fSample->fType!=Sample::SampleType::SIGNAL) {
                if(!fBkgOnly_postFit) {
                    fBkgOnly_postFit.reset(static_cast<TH1*>(allHistos.at(i)->Clone()));
                    fBkgOnly_postFit->SetDirectory(nullptr);
                } else {
                    fBkgOnly_postFit->Add(allHistos.at(i).get());
                }
            }
        } else {
            if(!fTot) {
                fTot.reset(static_cast<TH1*>(allHistos.at(i)->Clone("h_tot_postFit")));
            } else {
                fTot->Add(allHistos.at(i).get());
            }
        }
    }

    if (isPostFit) {
        fTot_postFit->SetDirectory(nullptr);
        if (fBkgOnly_postFit) fBkgOnly_postFit->SetDirectory(nullptr);
    } else {
        fTot->SetDirectory(nullptr);
    }

    //
    // MCstat set to 0 if disabled
    //
    if(!fUseStatErr || fUseGammaPulls){
        if (isPostFit) {
            for(int ibin = 1; ibin <= fTot_postFit->GetNbinsX(); ibin++){
                fTot_postFit->SetBinError(ibin, 0);
                if (fBkgOnly_postFit) fBkgOnly_postFit->SetBinError(ibin, 0);
            }
        } else {
            for(int ibin = 1; ibin <= fTot->GetNbinsX(); ibin++){
                fTot->SetBinError(ibin, 0);
            }
        }
    }

    plot->SetTot(isPostFit ? fTot_postFit.get() : fTot.get());
    if(fBinWidth>0) plot->SetBinWidth(fBinWidth);

    if(fBlindingThreshold>=0 || fManualBlindBins.size()>0){
        // fBlindedBinsPostFit and fBlindedBins are set in TRExFit::DrawAndSaveAll
        container.blindedBins = (isPostFit && !fKeepPrefitBlindedBins) ? fBlindedBinsPostFit : fBlindedBins;
        plot->SetBinBlinding( (isPostFit && !fKeepPrefitBlindedBins) ? fBlindedBinsPostFit : fBlindedBins);
    }
    plot->BlindData();

    //
    // Build error band
    //
    this->BuildPrePostFitErrorHistxRooFit(channel, fSampleHists.back()->fHist.get(), fitRes, isPostFit);

    //
    // Print chi2 info
    //
    if(fGetChi2 != 0 && TRExFitter::SHOWCHI2 && fRegionDataType == REALDATA) {
        if (fBlindedBins.size() == 0){
            plot->SetChi2KS(fChi2prob, -1, fChi2val, fNDF);
        }
        else{
            LOG(DEBUG) << "You specified CHI2 in PlotOptions, but region" << fName << " has blinded bins -- not showing Chi2\n";
        }
    }

    if (isPostFit) {
        container.errors = fErr_postFit.get();
        plot->SetTotAsym(fErr_postFit.get());
    } else {
        container.errors = fErr.get();
        plot->SetTotAsym(fErr.get());
    }
    plot->fPlotLabel = fPlotLabel;
    plot->fRatioYtitle = fRatioYtitle;
    plot->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType==TRExPlot::RATIOTYPE::DATAOVERMC){
        plot->fRatioType = TRExPlot::RATIOTYPE::DATAOVERB;
    }
    plot->fLabelX = fLabelX;
    plot->fLabelY = fLabelY;
    plot->fLegendX1 = fLegendX1;
    plot->fLegendX2 = fLegendX2;
    plot->fLegendY = fLegendY;
    if(fLogScale) opt += " log";
    if(fLogScaleX) opt += " LOGX";
    plot->fBinTickMarks = fBinTickMarks;
    plot->fBinTickMarksLabels = fBinTickMarksLabels;

    plot->Draw(opt);
    if (fDrawDataMinusBkg) {
    if (isPostFit) {
        plot->fBkgOnlyErr = fErrBkgOnly_postFit;
            plot->DrawDataMinusBkg(fDrawDataMinusBkgOnSignal);
        }
    }

    // Yaml
    YamlConverter converter{};
    converter.WritePlot(container, fFolder, fPlotSubdir, isPostFit);

    if (fHEPDataFormat) {
        converter.WritePlotHEPData(container, fFolder, isPostFit);
    }

    if (fUheppFormat) {
        converter.WritePlotUhepp(container, fFolder, isPostFit);
    }

    //
    // Print bin content and errors
    //
    LOG(DEBUG) << "--------------------\n";
    LOG(DEBUG) << "Final bin contents\n";
    LOG(DEBUG) << "--------------------\n";
    if (isPostFit) {
        for (int ibin = 1; ibin <= fTot_postFit->GetNbinsX(); ibin++) {
            LOG(DEBUG) << ibin << ":\t" << fTot_postFit->GetBinContent(ibin) << " +" << fErr_postFit->GetErrorYhigh(ibin-1) << " -" << fErr_postFit->GetErrorYlow(ibin-1) << "\n";
        }
    } else {
        for (int ibin = 1; ibin <= fTot->GetNbinsX(); ibin++) {
            LOG(DEBUG) << ibin << ":\t" << fTot->GetBinContent(ibin) << " +" << fErr->GetErrorYhigh(ibin-1) << " -" << fErr->GetErrorYlow(ibin-1) << "\n";
        }
    }

    //
    // Save in a root file...
    //
    gSystem->mkdir((fFitName+"/Histograms").c_str());
    std::string fileName = fFitName+"/Histograms/"+fName+fSuffix;
    if (isPostFit) {
        fileName += "_postFit";
    }
    fileName += ".root";
    LOG(INFO) << "Writing file " << fileName << "\n";
    std::unique_ptr<TFile> f(TFile::Open(fileName.c_str(),"RECREATE"));
    f->cd();
    if (isPostFit) {
        fErr_postFit->Write("",TObject::kOverwrite);
        fTot_postFit->Write("",TObject::kOverwrite);
        if(plot->h_tot_bkg_prefit) {
           plot->h_tot_bkg_prefit->Write("",TObject::kOverwrite);
        }
    } else {
        fErr->Write("",TObject::kOverwrite);
        fTot->Write("",TObject::kOverwrite);
    }
    f->Close();
    return plot;
}

//__________________________________________________________________________________
//
void Region::AddSelection(const std::string& selection){
    if(selection=="") return;
    if(fSelection=="1" || fSelection=="") fSelection = selection;
    else fSelection += " && "+selection;
}

//__________________________________________________________________________________
//
void Region::AddMCweight(const std::string& weight){
    if(weight=="") return;
    if(fMCweight=="1" || fMCweight=="") fMCweight = weight;
    else fMCweight += " * "+weight;
}

//__________________________________________________________________________________
//
void Region::SetVariable(const std::string& variable,int nbin,double xmin,double xmax, const std::string& corrVar1, const std::string& corrVar2){
    fVariable = variable;
    fCorrVar1 = corrVar1;
    fCorrVar2 = corrVar2;
    fNbins = nbin;
    fXmin = xmin;
    fXmax = xmax;
}

//__________________________________________________________________________________
//
void Region::SetAlternativeVariable(const std::string& variable, const std::string& sample){
    fAlternativeVariables[sample] = variable;
}

//__________________________________________________________________________________
//
void Region::SetAlternativeSelection(const std::string& selection, const std::string& sample){
    fAlternativeSelections[sample] = selection;
}

//__________________________________________________________________________________
//
bool Region::UseAlternativeVariable(const std::string& sample){
    std::vector<std::string> tmpVec;
    for(const auto& tmp : fAlternativeVariables){
        tmpVec.push_back(tmp.first);
    }
    if (Common::FindInStringVector(tmpVec,sample)<0){
        return false;
    }
    return true;
}

//__________________________________________________________________________________
//
bool Region::UseAlternativeSelection(const std::string& sample){
    std::vector<std::string> tmpVec;
    for(const auto& tmp : fAlternativeSelections){
        tmpVec.push_back(tmp.first);
    }
    if (Common::FindInStringVector(tmpVec,sample)<0){
        return false;
    }
    return true;
}

//__________________________________________________________________________________
//
std::string Region::GetAlternativeVariable(const std::string& sample) const{
    std::vector<std::string> tmpVec;
    std::vector<std::string> tmpVec2;
    for(const auto& tmp : fAlternativeVariables){
        tmpVec.push_back(tmp.first);
        tmpVec2.push_back(tmp.second);
    }
    const int idx = Common::FindInStringVector(tmpVec,sample);
    if(idx<0){
        return "";
    }
    return tmpVec2[idx];
}

//__________________________________________________________________________________
//
std::string Region::GetAlternativeSelection(const std::string& sample) const{
    std::vector<std::string> tmpVec;
    std::vector<std::string> tmpVec2;
    for(const auto& tmp : fAlternativeSelections){
        tmpVec.push_back(tmp.first);
        tmpVec2.push_back(tmp.second);
    }
    const int idx = Common::FindInStringVector(tmpVec,sample);
    if(idx<0){
        return "";
    }
    return tmpVec2[idx];
}

//__________________________________________________________________________________
//
void Region::SetVariableTitle(const std::string& name){
    fVariableTitle = name;
}

//__________________________________________________________________________________
//
void Region::SetLabel(const std::string& label, const std::string& shortLabel){
    if(label=="NONE"){fLabel = "";}
    else {fLabel = label;}
    if(shortLabel=="") fShortLabel = label;
    else fShortLabel = shortLabel;
}

//__________________________________________________________________________________
//
void Region::Print() const{
    LOG(INFO) << "    Region: " << fName << "\n";
    for(const auto& isample : fSampleHists) {
        isample->Print();
    }
}

//__________________________________________________________________________________
//
void Region::PrintSystTable(std::shared_ptr<xRooNode> pdf, const std::string& opt) const {
    bool isPostFit  = false; if(opt.find("post")!=string::npos)     isPostFit  = true;
    bool doClean    = false; if(opt.find("clean")!=string::npos)    doClean    = true;
    bool doCategory = false; if(opt.find("category")!=string::npos) doCategory = true;
    bool standalone = false; if(opt.find("standalone")!=string::npos)standalone= true;
    bool landscape  = false; if(opt.find("landscape")!=string::npos) landscape = true;
    bool footnotesize=false; if(opt.find("footnotesize")!=string::npos)footnotesize=true;;

    ofstream out;
    ofstream texout;
    ofstream out_cat;
    ofstream texout_cat;
    gSystem->mkdir(fFitName.c_str());
    gSystem->mkdir((fFitName+"/Tables").c_str());
    if(isPostFit){
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_postFit.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category_postFit.tex").c_str());
        }
    }
    else{
        out.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.txt").c_str());
        texout.open((fFitName+"/Tables/"+fName+fSuffix+"_syst.tex").c_str());
        if(doCategory){
            out_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.txt").c_str());
            texout_cat.open((fFitName+"/Tables/"+fName+fSuffix+"_syst_category.tex").c_str());
        }
    }

    out << " | ";
    if (standalone) {
        texout << "\\documentclass[10pt]{article}" << endl;
        texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
        texout << "\\begin{document}" << endl;
    }
    if (landscape) texout << "\\begin{landscape}" << endl;
    texout << "\\begin{table}[htbp]" << endl;
    texout << "\\begin{center}" << endl;
    if (footnotesize) texout << "\\footnotesize" << endl;
    texout << "\\begin{tabular}{|c" ;

    if(doCategory){
        out_cat << " | ";
        if (standalone) {
            texout_cat << "\\documentclass[10pt]{article}" << endl;
            texout_cat << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << endl;
            texout_cat << "\\begin{document}" << endl;
        }
        texout_cat << "\\begin{table}[htbp]" << endl;
        texout_cat << "\\begin{center}" << endl;
        texout_cat << "\\begin{tabular}{|c" ;
    }


    for(std::size_t ismp = 0; ismp < fSampleHists.size(); ismp++){
        const SampleHist* sh = fSampleHists.at(ismp).get();
        const Sample* s = sh->fSample;
        if(s->fType==Sample::SampleType::DATA) continue;
        if(s->fType==Sample::SampleType::GHOST) continue;
        if(s->fType==Sample::SampleType::EFT) continue;
        texout << "|c";
        if(doCategory) texout_cat << "|c";
    }
    texout << "|}" << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        texout_cat << "|}" << endl;
        texout_cat << "\\hline " << endl;
    }

    for(std::size_t ismp = 0; ismp < fSampleHists.size(); ismp++){
        const SampleHist* sh = fSampleHists.at(ismp).get();
        const Sample* s = sh->fSample;
        if(s->fType==Sample::SampleType::DATA) continue;
        if(s->fType==Sample::SampleType::GHOST) continue;
        if(s->fType==Sample::SampleType::EFT) continue;
        std::string title = s->fTitle;
        if(s->fTexTitle!="") title = s->fTexTitle;
        out << "      | " << s->fTitle;
        texout << "      & " << title;
        if(doCategory){
            out_cat << "      | " << s->fTitle;
            texout_cat << "      & " << title;
        }
    }
    out << " |" << endl;
    texout << " \\\\ " << endl;
    texout << "\\hline " << endl;
    if(doCategory){
        out_cat << " |" << endl;
        texout_cat << " \\\\ " << endl;
        texout_cat << "\\hline " << endl;
    }

    auto channel = pdf->find(fName);
    if (!channel) {
        LOG(ERROR) << "Cannot read channel: " << fName << "\n";
        exit(EXIT_FAILURE);
    }

    auto pars = channel->pars();

    std::vector<RooRealVar*> params;
    for (const auto& iparam : channel->floats()) {
        const std::string name = iparam->GetName();
        if (!isPostFit) {
            auto itr = std::find_if(fNormFactors.begin(), fNormFactors.end(), [&name](const auto& nf){return nf->fName == name;});
            if (itr != fNormFactors.end()) continue;
        }
        params.emplace_back(iparam->get<RooRealVar>());
    }

    // cache values to speed up the code
    std::vector<std::shared_ptr<xRooNode> > samples;
    for (const auto& isample : fSampleHists) {
        if (isample->fSample->fType == Sample::SampleType::DATA) continue;
        if (isample->fSample->fType == Sample::SampleType::GHOST) continue;
        if (isample->fSample->fType == Sample::SampleType::EFT) continue;
        samples.emplace_back(channel->at("samples")->find(isample->fName + "_" + fName + "_shapes"));
    }

    std::vector<double> sampleNominal;
    for (const auto& isample : samples) {
        sampleNominal.emplace_back(isample->GetContent());
    }

    for (const auto& iparam : params) {
        auto itrNames = std::find_if(TRExFitter::SYSTMAP.begin(), TRExFitter::SYSTMAP.end(),
                        [&iparam](const auto& element){return static_cast<std::string>(iparam->GetName()).find(element.first) != std::string::npos;});
        if (itrNames == TRExFitter::SYSTMAP.end()) {
            out << " | " << iparam->GetName();
        } else {
            out << " | " << itrNames->second;
        }
        auto itrTex = std::find_if(TRExFitter::SYSTTEX.begin(), TRExFitter::SYSTTEX.end(),
                      [&iparam](const auto& element){return static_cast<std::string>(iparam->GetName()).find(element.first) != std::string::npos;});
        if (itrTex == TRExFitter::SYSTTEX.end()) {
            if (itrNames == TRExFitter::SYSTMAP.end()) {
                texout << "  " << iparam->GetName();
            } else {
                std::string fixedTitle = itrNames->second;
                fixedTitle = Common::ReplaceString(fixedTitle,"#geq","$\\geq$");
                texout << "  " << fixedTitle;
            }
        } else {
            texout << "  " << itrTex->second;
        }

        const double value = iparam->getVal();
        const double up    = iparam->getErrorHi();
        const double down  = iparam->getErrorLo();

        std::size_t sampleIndex(0);
        for (const auto& isample : samples) {
            if (!isample) {
                out << " |    NA   ";
                texout << " &    NA   ";
                continue;
            }

            const double effectNominal = sampleNominal.at(sampleIndex);
            if (effectNominal == 0.) {
                out << " |    NA   ";
                texout << " &    NA   ";
                ++sampleIndex;
                continue;
            }

            // set to +1 sigma
            iparam->setVal(value+up);
            const double effectUp = isample->GetContent();

            // set to -1 sigma
            iparam->setVal(value+down);
            const double effectDown = isample->GetContent();

            // restore
            iparam->setVal(value);

            const double normUp = effectUp/effectNominal - 1.;
            const double normDown = effectDown/effectNominal -1.;

            out << " | " << normUp;
            texout << setprecision(3) << " & " << normUp;
            out << " / " << normDown;
            texout << setprecision(3) << " / " << normDown;
            ++sampleIndex;
        }
        out << " |\n";
        texout << " \\\\ \n";
    }

    if(doCategory){
        std::map<std::string, std::vector<std::string> > categories;
        for (const auto& isample : fSampleHists) {
            if (isample->fSample->fType == Sample::SampleType::DATA) continue;
            if (isample->fSample->fType == Sample::SampleType::GHOST) continue;
            if (isample->fSample->fType == Sample::SampleType::EFT) continue;
            for (const auto& isyst : isample->fSample->fSystematics) {
                auto itr = categories.find(isyst->fCategory);
                if (itr == categories.end()) {
                    categories.insert({isyst->fCategory, {isyst->fName}});
                } else {
                    auto itr2 = std::find(itr->second.begin(), itr->second.end(), isyst->fName);
                    if (itr2 == itr->second.end()) {
                        itr->second.emplace_back(isyst->fName);
                    }
                }
            }
        }

        for (const auto& icategory : categories) {
            out_cat << " | " << icategory.first;
            texout_cat << " " << icategory.first;

            std::vector<double> nominal;
            for (const auto& isample : samples) {
                nominal.emplace_back(isample->GetContent());
            }
            std::vector<std::pair<std::string, double> > nomValues;
            for (const auto& iparam : params) {
                auto itr = std::find_if(icategory.second.begin(), icategory.second.end(), [&iparam](const auto& element){return static_cast<std::string>(iparam->GetName()).find(element) != std::string::npos;});
                if (itr == icategory.second.end()) continue;
                const double value = iparam->getVal();
                const double up    = iparam->getErrorHi();

                nomValues.emplace_back(std::make_pair(iparam->GetName(), value));
                iparam->setVal(value+up);
            }

            std::vector<double> shifted;
            for (const auto& isample : samples) {
                shifted.emplace_back(isample->GetContent());
            }

            // restore
            for (const auto& nom : nomValues) {
                pars[nom.first]->get<RooRealVar>()->setVal(nom.second);
            }

            if (nominal.size() != shifted.size()) {
                LOG(ERROR) << "Sizes of the vectors do not match\n";
                return;
            }

            for (std::size_t i = 0; i < nominal.size(); ++i) {
                const double err = nominal.at(i) != 0. ? shifted.at(i)/nominal.at(i) - 1. : -999;
                out_cat << setprecision(3) << " | " << err;
                texout_cat << setprecision(3) << " & " << err;
            }

            out_cat << " |\n";
            texout_cat << " \\\\ \n";
        }
    }

    texout << "\\hline " << endl;
    texout << "\\end{tabular} " << endl;
    texout << "\\caption{Relative effect of each systematic on the yields.} " << endl;
    texout << "\\end{center} " << endl;
    texout << "\\end{table} " << endl;

    if(doCategory){
        texout_cat << "\\hline " << endl;
        texout_cat << "\\end{tabular} " << endl;
        texout_cat << "\\caption{Realtive effect of each group of systematics on the yields.} " << endl;
        texout_cat << "\\end{center} " << endl;
        texout_cat << "\\end{table} " << endl;
    }
    if (landscape) texout << "\\end{landscape}" << endl;
    if (standalone) {
        texout << "\\end{document}" << endl;
    }

    if (doCategory && standalone) {
        texout_cat << "\\end{document}" << endl;
    }

    if(doClean){
        std::string shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
        if(isPostFit) shellcommand += "_postFit";
        shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
        shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst_clean.tex";
        gSystem->Exec(shellcommand.c_str());
        if(doCategory){
            shellcommand = "cat "+fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += ".tex|sed -e \"s/\\#/ /g\" > ";
            shellcommand += fFitName+"/Tables/"+fName+fSuffix+"_syst";
            shellcommand += "_category";
            if(isPostFit) shellcommand += "_postFit";
            shellcommand += "_clean.tex";
            gSystem->Exec(shellcommand.c_str());
        }
    }
}

//___________________________________________________________
// function to get pre/post-fit agreement
std::pair<double,int> Region::GetChi2Test(const bool isPostFit,
                                          const TH1* h_data,
                                          const TH1* h_nominal,
                                          const std::vector< std::shared_ptr<TH1> >& h_up,
                                          const std::vector< string >& systNames,
                                          const CorrelationMatrix *matrix) const{
    const unsigned int nbins = h_nominal->GetNbinsX();
    int ndf = 0;
    for(unsigned int i=0;i<nbins;++i){
        if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), i+1) != fBlindedBins.end()) continue;
        if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), i+1) != fBlindedBinsPostFit.end()) continue;
        if (std::find(fDropBins.begin(), fDropBins.end(), i+1) != fDropBins.end()) continue;
        if (h_data->GetBinContent(i+1) == 0. && h_nominal->GetBinContent(i+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos (e.g. because X-axis range contains non-physical values)
        ndf++;
    }
    //
    //Speed Up: remove irrelevant systematics (which would give in any case 0 correlation)
    std::vector< string > EffectiveSystNames;
    std::vector< unsigned int > EffectiveSystIndex;
    for(unsigned int n=0;n<systNames.size();++n){
        if(matrix!=nullptr){
            if (matrix->NPIsPresent(systNames[n])) {
               EffectiveSystNames.push_back(systNames[n]);
               EffectiveSystIndex.push_back(n);
            }
            else LOG(DEBUG) << "Will skip syst. " << systNames[n] << "\n";
        }
        else {
            EffectiveSystNames.push_back(systNames[n]);
            EffectiveSystIndex.push_back(n);
        }
    }
    const unsigned int nsyst=EffectiveSystNames.size();

    // precalculate correlations
    std::vector<std::vector<double> > correlation(nsyst, std::vector<double>(nsyst));
    if (matrix) {
        for(unsigned int n = 0; n < nsyst; ++n){ //n!=m, run only across correlated systs
            for(unsigned int m = 0; m < nsyst; ++m){ //n!=m, run only across correlated systs
                if (n == m) continue;;
                correlation[n][m] = matrix->GetCorrelation(EffectiveSystNames[n],EffectiveSystNames[m]);
            }
        }
    }

    //
    TMatrixD C(ndf,ndf);
    //
    int ibin = 0;
    for(unsigned int i = 0; i < nbins; ++i) {
        if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), i+1) != fBlindedBins.end()) continue;
        if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), i+1) != fBlindedBinsPostFit.end()) continue;
        if (std::find(fDropBins.begin(), fDropBins.end(), i+1) != fDropBins.end()) continue;
        if (h_data->GetBinContent(i+1) == 0. && h_nominal->GetBinContent(i+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
        const double ynom_i = h_nominal->GetBinContent(i+1);
        int jbin = ibin;
        for(unsigned int j = i; j < nbins; ++j) {
            double sum = 0.;
            if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), j+1) != fBlindedBins.end()) continue;
            if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), j+1) != fBlindedBinsPostFit.end()) continue;
            if (std::find(fDropBins.begin(), fDropBins.end(), j+1) != fDropBins.end()) continue;
            if (h_data->GetBinContent(j+1) == 0. && h_nominal->GetBinContent(j+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
            const double ynom_j = h_nominal->GetBinContent(j+1);
            for(unsigned int n = 0; n < nsyst; ++n) { //n!=m, run only across correlated systs
                if (!matrix) continue;
                const double ysyst_i_n = h_up[EffectiveSystIndex[n]]->GetBinContent(i+1);
                for(unsigned int m = 0; m < nsyst; ++m) {
                    if (n == m) continue;
                    const double ysyst_j_m = h_up[EffectiveSystIndex[m]]->GetBinContent(j+1);
                    const double corr = correlation[n][m];
                    sum += (ysyst_i_n-ynom_i) * corr * (ysyst_j_m-ynom_j); // exploiting the symmetry
                }
            }
            for(unsigned int n=0;n<systNames.size();++n){ // n==m, all systs, corr=1
                const double ysyst_i_n = h_up[n]->GetBinContent(i+1), ysyst_j_n = h_up[n]->GetBinContent(j+1);
                sum += (ysyst_i_n-ynom_i) /* * 1.0 */ * (ysyst_j_n-ynom_j);
            }
            if(i==j && ynom_i>0) {
                sum += ynom_i;  // add stat uncertainty to diagonal
                // add MC stat to prefit or when gamma pulls are not used
                if (!isPostFit || !fUseGammaPulls) sum += h_nominal->GetBinError(i+1) * h_nominal->GetBinError(i+1);
            }
            C[ibin][jbin] = sum;
            if (ibin != jbin) {
                C[jbin][ibin] = sum;
            }
            ++jbin;
        }
        ++ibin;
    }
    //
    if(TRExFitter::DEBUGLEVEL > 1) C.Print();
    //
    // Invert the matrix
    C.Invert();
    if(TRExFitter::DEBUGLEVEL > 1) C.Print();
    //
    double chi2 = 0.;
    ibin = 0;
    for(unsigned int i=0;i<nbins;++i){
        if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), i+1) != fBlindedBins.end()) continue;
        if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), i+1) != fBlindedBinsPostFit.end()) continue;
        if (std::find(fDropBins.begin(), fDropBins.end(), i+1) != fDropBins.end()) continue;
        if (h_data->GetBinContent(i+1) == 0. && h_nominal->GetBinContent(i+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
        const double ynom_i = h_nominal->GetBinContent(i+1);
        const double ydata_i = h_data->GetBinContent(i+1);
        int jbin = 0;
        for(unsigned int j=0;j<nbins;j++){
            if ((!isPostFit || fKeepPrefitBlindedBins) && std::find(fBlindedBins.begin(), fBlindedBins.end(), j+1) != fBlindedBins.end()) continue;
            if ((isPostFit && !fKeepPrefitBlindedBins) && std::find(fBlindedBinsPostFit.begin(), fBlindedBinsPostFit.end(), j+1) != fBlindedBinsPostFit.end()) continue;
            if (std::find(fDropBins.begin(), fDropBins.end(), j+1) != fDropBins.end()) continue;
            if (h_data->GetBinContent(j+1) == 0. && h_nominal->GetBinContent(j+1) <= fSampleHists.size()*1e-6) continue; // skip bins when such bins are not filled in both data and nominal histos
            const double ynom_j = h_nominal->GetBinContent(j+1);
            const double ydata_j = h_data->GetBinContent(j+1);
            chi2 += (ydata_i - ynom_i)*C[ibin][jbin]*(ydata_j-ynom_j);
            ++jbin;
        }
        ++ibin;
    }
    //
    return std::make_pair(chi2,ndf);
}

//___________________________________________________________
//
std::vector<std::pair<std::string,double> > Region::SystPruning(PruningUtil *pu) {
    std::vector<std::pair<std::string, double> > result;
    std::unique_ptr<TH1> hTot(nullptr);
    if(pu->GetStrategy() == 1){
        hTot = GetTotHist(false); // don't include signal
    }
    else if(pu->GetStrategy() == 2 || pu->GetStrategy() == 3){
        hTot = GetTotHist(true); // include signal
    }
    for(auto& sh : fSampleHists){
        if (!(sh->fSample->fNoPruning)) {
            sh->SystPruning(pu, hTot.get());
            //
            // flag overall systematics as no shape also for pruning purposes
            for(const auto& syh : sh->fSyst){
                if(!syh) continue;
                if(!syh->fSystematic) continue;
                if(syh->fSystematic->fNoPruning) {
                    syh->fShapePruned = false;
                    syh->fNormPruned  = false;
                    syh->fBadShape    = false;
                    syh->fBadNorm     = false;
                    continue;
                }
                if(syh->fSystematic->fType==Systematic::OVERALL){
                    syh->fShapePruned = true;
                }
            }
        }
        //
        // add by-hand dropping of shape or norm
        for(const auto& syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            if(syh->fSystematic->fNoPruning) continue;
            if( Common::FindInStringVector(syh->fSystematic->fDropShapeIn,fName)>=0  ||
                Common::FindInStringVector(syh->fSystematic->fDropShapeIn,fName)>=0 ||
                Common::FindInStringVector(syh->fSystematic->fDropShapeIn, "all")>=0
            ){
                syh->fShapePruned = true;
            }
            if( Common::FindInStringVector(syh->fSystematic->fDropNormIn,fName)>=0 ||
                Common::FindInStringVector(syh->fSystematic->fDropNormIn,fName)>=0 ||
                Common::FindInStringVector(syh->fSystematic->fDropNormIn, "all")>=0
            ){
                if(sh->fName.rfind("SM_", 0) == 0 && sh->fName.find("_bin") != std::string::npos){
                    LOG(WARNING) << "Although you specified drop norm for pruning this will not work on an EFT sample for now\n";
                    LOG(WARNING) << "will not drop norm for sample " <<  sh->fName << "\n";
                    LOG(WARNING) << "check " << static_cast<int>(sh->fSample->fType) << "\n";
                }
                else{
                    syh->fNormPruned = true;
                }
            }
        }

        // get the heuristic impact
        for(const auto& syh : sh->fSyst) {
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            const std::string& name = syh->fSystematic->fNuisanceParameter;
            const double impact = syh->fMagnitude;
            result.emplace_back(name, impact);
        }
    }
    //
    // reference pruning
    for(const auto& sh : fSampleHists){
        for(const auto& syh : sh->fSyst){
            if(!syh) continue;
            if(!syh->fSystematic) continue;
            std::shared_ptr<Systematic> syst = syh->fSystematic;
            if(syst->fReferencePruning == "") continue;
            std::shared_ptr<SampleHist> refSmpH = GetSampleHist(syst->fReferencePruning);
            if(refSmpH==nullptr){
                LOG(WARNING) << "Cannot find reference pruning sample: " << syst->fReferencePruning << " in region: " << fName << "\n";
                continue;
            }
            std::shared_ptr<SystematicHist> refSysH = refSmpH->GetSystematic(syst->fName);
            if(refSysH==nullptr){
                LOG(WARNING) << "Cannot find systematic " << syst->fName << " for reference pruning sample " << syst->fReferencePruning << "\n";
                continue;
            }
            syh->fNormPruned = refSysH->fNormPruned;
            syh->fShapePruned = refSysH->fShapePruned;
            syh->fBadShape = refSysH->fBadShape;
            syh->fBadNorm = refSysH->fBadNorm;
        }
    }
    return result;
}

//___________________________________________________________
//
std::unique_ptr<TH1> Region::GetTotHist(bool includeSignal) {
    std::unique_ptr<TH1> hTot(nullptr);
    for(const auto& sh : fSampleHists){
        if(!sh->fSample) continue;
        if(sh->fSample->fType==Sample::SampleType::GHOST) continue;
        if(sh->fSample->fType==Sample::SampleType::DATA) continue;
        if(sh->fSample->fType==Sample::SampleType::EFT) continue;
        if(!includeSignal && sh->fSample->fType==Sample::SampleType::SIGNAL) continue;
        if(!sh->fHist) continue;
        std::unique_ptr<TH1> hTmp(static_cast<TH1*>(sh->fHist->Clone()));
        if(!hTot) {
            hTot = std::move(hTmp);
        } else {
            hTot->Add(hTmp.get());
        }
    }
    if (hTot) hTot->SetDirectory(nullptr);
    return hTot;
}

//___________________________________________________________
//
void Region::SavePreFitUncertaintyAndTotalMCObjects() {
    auto originalDir = gDirectory;
    TString uncBandFileName = TString::Format("%s/Histograms/%s%s_preFit.root", fFitName.c_str(), fName.c_str(), fSuffix.c_str());
    gSystem->mkdir((fFitName + "/Histograms").c_str());
    LOG(DEBUG) << TString::Format("Writing file %s", uncBandFileName.Data()) << "\n";
    std::unique_ptr<TFile> f(TFile::Open(uncBandFileName, "RECREATE"));
    f->cd();
    fErr->Write("", TObject::kOverwrite);
    fTot->Write("", TObject::kOverwrite);
    originalDir->cd();
    f->Close();
}

//___________________________________________________________
//
void Region::AddBinsToFile() const {
    std::fstream file;
    file.open(fJobBinningsPath + fName + ".txt", std::fstream::app);

    if (!file.is_open() || !file.good()) {
        LOG(ERROR) << "Cannot open the " << fJobBinningsPath << fName << ".txt file\n";
        return;
    }

    file << fNbins << "\n";

    file.close();
}

//___________________________________________________________
//
void Region::BuildPrePostFitErrorHistxRooFit(const std::shared_ptr<xRooNode>& channel,
                                             const TH1* nominal,
                                             const FitResults* fitRes,
                                             const bool isPostFit) {

    if (isPostFit && fToysForErrorBand > 0) {
        LOG(INFO) << "Using " << fToysForErrorBand << " toys to build the error band\n";
    }

    if (isPostFit) {
        if (TRExFitter::SHOWSTACKSIG) {
            fErr_postFit = Common::ErrorFromxRooNode(*channel, nominal, fitRes, fToysForErrorBand);
        } else {
            const std::string bkgSamples = Common::NonSignalSampleList(fSampleHists, fName);
            xRooNode bkgOnlySample = FitUtils::GetReducedOrSingleSample(channel, bkgSamples);
            fErr_postFit = Common::ErrorFromxRooNode(bkgOnlySample, nominal, fitRes, fToysForErrorBand);
        }
    } else {
        if (TRExFitter::SHOWSTACKSIG) {
            fErr = Common::ErrorFromxRooNode(*channel, nominal, nullptr, -1);
        } else {
            const std::string bkgSamples = Common::NonSignalSampleList(fSampleHists, fName);
            xRooNode bkgOnlySample = FitUtils::GetReducedOrSingleSample(channel, bkgSamples);
            fErr = Common::ErrorFromxRooNode(bkgOnlySample, nominal, nullptr, -1);
        }
    }

    if (isPostFit && fDrawDataMinusBkg) {
        // generate the bkg only error band by selecting samples that are not signal samples
        const std::string bkgSamples = Common::NonSignalSampleList(fSampleHists, fName);
        xRooNode bkgOnlySample = FitUtils::GetReducedOrSingleSample(channel, bkgSamples);

        fErrBkgOnly_postFit = Common::ErrorFromxRooNode(bkgOnlySample, nominal, fitRes, fToysForErrorBand);
    }

    //build the per NP variations for chi^2
    std::vector<std::string> params;
    std::vector<std::shared_ptr<TH1> > systHistos;

    if (fGetChi2 == 2 && fHasData && fRegionDataType == REALDATA) {
        for (const auto& iparam : channel->floats()) {
            const std::string name = iparam->GetName();
            // skip NFs for prefit as they are not defined
            if (!isPostFit) {
                auto itr = std::find_if(fNormFactors.begin(), fNormFactors.end(), [&name](const auto& element){return element->fName == name;});
                if (itr != fNormFactors.end()) continue;
            }
            params.emplace_back(name);
        }

        auto pars = channel->pars();

        // fill the histograms needed for chi^2
        for (const auto& iparam : params) {
            auto par = pars[iparam]->get<RooRealVar>();
            const double value = par->getVal();
            const double up    = par->getErrorHi();

            // set to +1 sigma
            par->setVal(value+up);
            systHistos.emplace_back(Common::HistoFromxRooNode(channel, nominal));

            // restore
            par->setVal(value);
        }
    }

    // Goodness of post-fit
    if(fGetChi2 != 0 && fRegionDataType == REALDATA) {
        if (!fHasData){
            LOG(WARNING) << "Data histogram is nullptr, cannot calculate Chi2 agreement.\n";
            LOG(WARNING) << "Maybe you do not have data sample defined?\n";
            return ;
        }
        // remove blinded bins
        std::unique_ptr<TH1> h_data(static_cast<TH1*>(fData->fHist->Clone()));
        h_data->SetDirectory(nullptr);
        std::pair<double, double> chi2;
        if (isPostFit) {
            chi2 = GetChi2Test(true, h_data.get(), fTot_postFit.get(), systHistos, params, fitRes->GetCorrelationMatrix());
        } else {
            chi2 = GetChi2Test(false, h_data.get(), fTot.get(), systHistos, params, nullptr);
        }
        fChi2val = chi2.first;
        fNDF = chi2.second;

        // the chi^2/NDF is only chi^2 distributed IF:
        // a) it is for a region not used in the fit (validation region)
        // b) AND it does not have large statistical correlation with the observables used in the fit
        // OR if it is prefit
        fChi2prob = -1;
        if (fRegionType == Region::RegionType::VALIDATION || !isPostFit) {
            fChi2prob = ROOT::Math::chisquared_cdf_c(chi2.first, chi2.second);
        }

        // Do not print chi2 results if any bins are blinded in this region
        if (fBlindedBins.size() == 0 ){
            LOG(INFO) << "----------------------- ---------------------------- -----------------------\n";
            if (isPostFit) {
                LOG(INFO) << "----------------------- POST-FIT AGREEMENT EVALUATION -----------------------\n";
            } else {
                LOG(INFO) << "----------------------- PRE-FIT AGREEMENT EVALUATION -----------------------\n";
            }
            if(fGetChi2==1) {
                LOG(INFO) << "----------------------- -------- STAT-ONLY --------- -----------------------\n";
            }
            LOG(INFO) << "--- REGION " << fName << ":\n";
            LOG(INFO) << "  chi2        = " << fChi2val << "\n";
            LOG(INFO) << "  ndof        = " << fNDF << "\n";
            if (fRegionType == Region::RegionType::VALIDATION) {
                LOG(INFO) << "  probability = " << fChi2prob << "\n";
            }
            LOG(INFO) << "----------------------- ---------------------------- -----------------------\n";
            LOG(INFO) << "----------------------- ---------------------------- -----------------------\n";
        }
    }
}

//___________________________________________________________
//
void Region::SetAndAddNbins(const int n, const bool addToFile) {
    fNbins = n;
    if (addToFile) this->AddBinsToFile();
}

//___________________________________________________________
//
double Region::GetIntegral(const std::string& sampleName, const bool isPostfit) const {
    const auto& map = isPostfit ? fPostfitYields : fPrefitYields;
    auto itr = map.find(sampleName);

    if (itr == map.end()) {
        LOG(WARNING) << "Cannot find sample: " << sampleName << " returning -1\n";
        return -1;
    }

    return itr->second;
}

//___________________________________________________________
//
double Region::GetIntegralBonlyNorm(const std::string& sampleName) const {
    auto itr = fPostfitYieldsBonlyNorm.find(sampleName);

    if (itr == fPostfitYieldsBonlyNorm.end()) {
        LOG(WARNING) << "Cannot find sample: " << sampleName << " returning -1\n";
        return -1;
    }

    return itr->second;
}
