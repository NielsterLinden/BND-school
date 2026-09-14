// Class include
#include "TRExFitter/TRExFit.h"

// Framework includes
#include "TRExFitter/ConfigParser.h"
#include "TRExFitter/ConfigReader.h"
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/EFTProcessor.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FittingTool.h"
#include "TRExFitter/FitToys.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/LikelihoodScanManager.h"
#include "TRExFitter/LimitToys.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/NuisParameter.h"
#include "TRExFitter/PruningUtil.h"
#include "TRExFitter/RankingManager.h"
#include "TRExFitter/ReparametrizationManager.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/SignificanceToys.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/TRExPlot.h"
#include "TRExFitter/TruthSample.h"
#include "TRExFitter/UncertaintyPlotter.h"
#include "TRExFitter/Unfolding.h"
#include "TRExFitter/UnfoldingSample.h"
#include "TRExFitter/UnfoldingSystematic.h"
#include "TRExFitter/YamlConverter.h"

// UnfoldingCode includes
#include "UnfoldingCode/UnfoldingCode/UnfoldingTools.h"
#include "UnfoldingCode/UnfoldingCode/UnfoldingResult.h"

//Roofit headers
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFitHS3/RooJSONFactoryWSTool.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"

// RooStats includes
#include "RooStats/HypoTestInverterResult.h"

//HistFactory headers
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"

// Style stuff
#include "StyleUtils/TRExStyle.h"
#include "StyleUtils/TRExLabels.h"
#include "StyleUtils/TRExUtils.h"

// ROOT includes
#include "TCanvas.h"
#include "TF1.h"
#include "TF3.h"
#include "TFile.h"
#include "v5/TFormula.h"  // v5 to increase operator limit
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPie.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "xRooFit/xRooFit.h"
#include "xRooFit/xRooHypoSpace.h"
#include "xRooFit/xRooNode.h"

// c++ includes
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <filesystem>
#include <fstream>

using namespace RooFit;

// -------------------------------------------------------------------------------------------------
// class TRExFit

//__________________________________________________________________________________
//
TRExFit::TRExFit(const std::string& name) :
    fName(name),
    fUseFriend(false),
    fDir(""),
    fLabel(""),
    fInputFolder(""),
    fInputName(name),
    fBinningsPath(name+ "/Binnings/"),
    fUseStatErr(false),
    fStatErrThres(0.05),
    fUseGammaPulls(true),
    fLumi(1.),
    fLumiScale(1.),
    fThresholdSystPruning_Normalisation(-1),
    fThresholdSystPruning_Shape(-1),
    fThresholdSystLarge(-1),
    fMCweight("1"),
    fSelection("1"),
    fFitResults(nullptr),
    fWithPullTables(false),
    fIntCode(4),
    fInputType(HIST),
    fSystDataPlot_upFrame(false),
    fStatOnly(false),
    fGammasInStatOnly(false),
    fStatOnlyFit(false),
    fFixNPforStatOnlyFit(false),
    fYmin(0),
    fYmax(0),
    fRatioYmin(0.5),
    fRatioYmax(1.5),
    fRatioYminPostFit(0.5),
    fRatioYmaxPostFit(1.5),
    fRatioYtitle(""),
    fRatioType(TRExPlot::RATIOTYPE::DATAOVERMC),
    fLumiLabel("XX fb^{-1}"),
    fCmeLabel("13 TeV"),
    fSuffix(""),
    fSaveSuffix(""),
    fUpdate(false),
    fKeepPruning(false),
    fBlindingThreshold(-1),
    fBlindingType(Common::SOVERB),
    fRankingMaxNP(10),
    fRankingOnly("all"),
    fRankingPlot("MERGE"),
    fImageFormat("png"),
    fPlotLabel("Internal"),
    fDoSummaryPlot(true),
    fDoMergedPlot(false),
    fDoTables(true),
    fDoSignalRegionsPlot(true),
    fDoPieChartPlot(true),
    fGroupedImpactCategory("all"),
    fSummaryPrefix(""),
    fFitType(UNDEFINED),
    fFitRegion(CRSR),
    fFitNPValuesFromFitResultsFile(""),
    fInjectGlobalObservables(false),
    fPOIAsimovOverride(""),
    fFitIsBlind(false),
    fUseRnd(false),
    fRndRange(0.1),
    fRndSeed(-999),
    fCreateCache(false),
    fUseCache(false),
    fLHscanMin(999999),
    fLHscanMax(-999999),
    fLHscanSteps(30),
    fLHscanStep(-1),
    fLHscanMinY(999999),
    fLHscanMaxY(-999999),
    fLHscanStepsY(30),
    fLHscanStepY(-1),
    fWorkspaceFileName(""),
    fDoGroupedSystImpactTable(false),
    fLimitType(LimitType::ASYMPTOTIC),
    fLimitIsBlind(false),
    fSignalInjection(false),
    fSignalInjectionValue(0),
    fLimitParamName("parameter"),
    fLimitParamValue(0),
    fLimitOutputPrefixName("myLimit"),
    fLimitsConfidence(0.95),
    fLimitUseAutoDiff(false),
    fSignificanceType(SignificanceType::ASYMPTOTIC),
    fSignificanceIsBlind(false),
    fSignificanceDoInjection(false),
    fSignificancePOIAsimov(0),
    fSignificanceParamName("parameter"),
    fSignificanceParamValue(0),
    fSignificanceOutputPrefixName("mySignificance"),
    fSignificanceUseAutoDiff(false),
    fCleanTables(false),
    fSystCategoryTables(false),
    fKeepPrefitBlindedBins(false),
    fCustomAsimov(""),
    fWriteCustomAsimovToWS(false),
    fTableOptions("STANDALONE"),
    fGetChi2(0), // 0: no, 1: stat-only, 2: with syst
    fSmoothOption(HistoTools::SmoothOption::MAXVARIATION),
    fSuppressNegativeBinWarnings(false),
    fTemplateInterpolationOption(TRExFit::LINEAR),
    fBootstrap(""),
    fBootstrapSyst(""),
    fBootstrapSample(""),
    fBootstrapNomHistos(""),
    fBootstrapIdx(-1),
    fDecorrSuff("_decor"),
    fDoNonProfileFit(false),
    fNonProfileFitSystThreshold(0),
    fFitToys(0),
    fToysHistoNbins(50),
    fSmoothMorphingTemplates(""),
    fPOIPrecision(2),
    fUsePDGRounding(false),
    fUsePDGRoundingTxt(false),
    fUsePDGRoundingTex(false),
    fPropagateSystsForMorphing(false),
    fPruningType(SEPARATESAMPLE),
    fLabelX(-1),
    fLabelY(-1),
    fLegendX1(-1),
    fLegendX2(-1),
    fLegendY(-1),
    fLabelXSummary(-1),
    fLabelYSummary(-1),
    fLegendX1Summary(-1),
    fLegendX2Summary(-1),
    fLegendYSummary(-1),
    fLabelXMerge(-1),
    fLabelYMerge(-1),
    fLegendX1Merge(-1),
    fLegendX2Merge(-1),
    fLegendYMerge(-1),
    fLegendNColumns(2),
    fLegendNColumnsSummary(3),
    fLegendNColumnsMerge(3),
    fShowRatioPad(true),
    fShowRatioPadSummary(true),
    fShowRatioPadMerge(true),
    fExcludeFromMorphing(""),
    fDoSystNormalizationPlots(true),
    fDebugNev(-1),
    fCPU(1),
    fHasAcceptance(false),
    fPruningShapeOption(PruningUtil::SHAPEOPTION::MAXBIN),
    fSummaryLogY(true),
    fUseInFit(true),
    fUseInComparison(true),
    fReorderNPs(false),
    fBlindSRs(false),
    fHEPDataFormat(false),
    fUheppFormat(false),
    fAlternativeShapeHistFactory(false),
    fFitStrategy(-1),
    fBinnedLikelihood(false),
    fRemoveLargeSyst(false),
    fRemoveSystOnEmptySample(false),
    fValidationPruning(false),
    fUsePOISinRanking(false),
    fUseHesseBeforeMigrad(false),
    fUseNllInLHscan(true),
    fLimitToysStepsSplusB(100),
    fLimitToysStepsB(100),
    fLimitToysScanSteps(21),
    fLimitToysScanMin(0.),
    fLimitToysScanMax(10.),
    fToysSeed(1234),
    fLimitToysSeed(1234),
    fLimitPlot(true),
    fLimitFile(true),
    fLimitToysSuffix(""),
    fLimitFitStrategy(2),
    fLimitToysUsexRooFit(ToysUsexRooFit::FALSE),
    fSignificanceToysStepsSplusB(0),
    fSignificanceToysStepsB(100),
    fSignificanceToysSeed(1234),
    fSignificancePlot(true),
    fSignificanceFitStrategy(2),
    fSignificanceToysUsexRooFit(ToysUsexRooFit::FALSE),
    fDataWeighted(false),
    fRegularizationType(0),
    fSpeedUpFit(false),
    fNPCutOff(0.3),
    fDoExtendedCovariances(false),
    fUnfoldingShowStat(false),
    fUnfoldingShowUncertaintyBreakdown(false),
    fUnfoldingUncertaintyBreakdownTotal(true),
    fCombinerFormat(false),
    fToysStatOutput(false),
    fToysNpValuesFile(""),
    fUseRebinned(false),
    fApplyGammaCorrection(true),
    fErrorSigma(1.0),
    fShapeFactorReparametrisation(false),
    fHasTemplateMorphing(false),
    fTemplateMorphingFitDimensionality(TemplateMorpher::FIT_DIMENSIONALITY::ONE_DIMENSION),
    fMorphingSetting(Morphing()),
    fFitResultsRootFile(""),
    fRecreateBinningFiles(false),
    fUseHesse(true),
    fWorkspaceCreationType(TRExFit::WorkspaceCreationType::BOTH),
    fHasValidationRegions(true),
    fHasDropBinRegions(false),
    fErrorDecomposition(true),
    fMaximumNumberFCNcalls(-1),
    fToleranceScale(1.),
    fShiftGlobalObservablesInRanking(false),
    fUseAutoDiff(false),
    fNoPrePostFitLabel(false),
    fToysLHScanForAll(false),
    fToysOnlyStatFluctuation(false),
    fNLLOffset("initial"),
    fEFTConfig(EFTConfig()),
    fProducePerRegionWS(false),
    fDrawPruningPlot(true),
    fPreFitLabel("Pre-fit"),
    fPostFitLabel("Post-fit"),
    fToysForErrorBand(-1),
    fHasDataSetInConfig(false)
{
    TRExFitter::IMAGEFORMAT.emplace_back("png");
    // increase operator limit to be able to handle very long expressions (cuts/weights)
    ROOT::v5::TFormula::SetMaxima(100000,1000,1000000);
    // cppcheck-suppress syntaxError
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,38,0)
    RooRealVar::enableSilentClipping();
    #endif
}

//__________________________________________________________________________________
//
TRExFit::~TRExFit() {
}

//__________________________________________________________________________________
//
void TRExFit::SetPOI(const std::string& name, const std::string& unit){
    fPOIs.clear();
    fPOIs.emplace_back(name);
    fPOIunit[name] = unit;
}

//__________________________________________________________________________________
//
void TRExFit::AddPOI(const std::string& name, const std::string& unit){
    if(Common::FindInStringVector(fPOIs,name)<0){
        fPOIs.emplace_back(name);
        fPOIunit[name] = unit;
    }
    else{
        LOG(WARNING) << "POI " << name << " already set. Skipping.\n";
    }
}

//__________________________________________________________________________________
//
void TRExFit::SetStatErrorConfig(bool useIt, double thres, const std::string& cons){
    fUseStatErr = useIt;
    fStatErrThres = thres;
    fStatErrCons = cons;
}

//__________________________________________________________________________________
//
void TRExFit::SetLumi(const double lumi){
    fLumi = lumi;
}

//__________________________________________________________________________________
//
void TRExFit::SetFitType(FitType type){
    fFitType = type;
}

//__________________________________________________________________________________
//
void TRExFit::SetLimitType(LimitType type){
    fLimitType = type;
}

//__________________________________________________________________________________
//
void TRExFit::SetSignificanceType(SignificanceType type){
    fSignificanceType = type;
}

//__________________________________________________________________________________
//
void TRExFit::SetFitRegion(FitRegion region){
    fFitRegion = region;
}

//__________________________________________________________________________________
//
std::shared_ptr<Sample> TRExFit::NewSample(const std::string& name, Sample::SampleType type){
    fSamples.emplace_back(new Sample(name,type));
    //
    return fSamples.back();
}

//__________________________________________________________________________________
//
std::shared_ptr<Systematic> TRExFit::NewSystematic(const std::string& name){
    fSystematics.emplace_back(new Systematic(name));
    return fSystematics.back();
}

//__________________________________________________________________________________
//
Region* TRExFit::NewRegion(const std::string& name){
    fRegions.emplace_back(std::make_unique<Region>(name));
    //
    fRegions.back()->fFitName = fName;
    fRegions.back()->fSuffix = fSuffix;
    fRegions.back()->fFitLabel = fLabel;
    fRegions.back()->fFitType = fFitType;
    fRegions.back()->fPOIs = fPOIs;
    fRegions.back()->fIntCode_overall = fIntCode;
    fRegions.back()->fIntCode_shape   = (fIntCode == 4 ? 0 : fIntCode);
    fRegions.back()->fLumiScale = fLumiScale;
    fRegions.back()->fBlindingThreshold = fBlindingThreshold;
    fRegions.back()->fBlindingType = fBlindingType;
    fRegions.back()->fKeepPrefitBlindedBins = fKeepPrefitBlindedBins;
    fRegions.back()->fRatioYmax = fRatioYmax;
    fRegions.back()->fRatioYmin = fRatioYmin;
    fRegions.back()->fRatioYmaxPostFit = fRatioYmaxPostFit;
    fRegions.back()->fRatioYminPostFit = fRatioYminPostFit;
    fRegions.back()->fRatioYtitle = fRatioYtitle;
    fRegions.back()->fRatioType = fRatioType;
    fRegions.back()->fLabelX = fLabelX;
    fRegions.back()->fLabelY = fLabelY;
    fRegions.back()->fLegendX1 = fLegendX1;
    fRegions.back()->fLegendX2 = fLegendX2;
    fRegions.back()->fLegendY = fLegendY;
    fRegions.back()->fLegendNColumns = fLegendNColumns;
    fRegions.back()->fNoPrePostFitLabel = fNoPrePostFitLabel;
    //
    return fRegions.back().get();
}

//__________________________________________________________________________________
//
void TRExFit::SkippedRegion(const std::string& name){
    fSkippedRegions.emplace_back(name);
}

//__________________________________________________________________________________
//
void TRExFit::AddNtuplePath(const std::string& path){
    fNtuplePaths.push_back(path);
}

//__________________________________________________________________________________
//
void TRExFit::AddFriendPath(const std::string& path){
    fFriendPaths.push_back(path);
}

//__________________________________________________________________________________
//
void TRExFit::SetMCweight(const std::string &weight){
    fMCweight = weight;
}

//__________________________________________________________________________________
//
void TRExFit::SetSelection(const std::string& selection){
    fSelection = selection;
}

//__________________________________________________________________________________
//
void TRExFit::SetNtupleName(const std::string& name){
    fNtupleNames.clear();
    fNtupleNames.push_back(name);
}

//__________________________________________________________________________________
//
void TRExFit::SetNtupleFile(const std::string& name){
    fNtupleFiles.clear();
    fNtupleFiles.push_back(name);
}

//__________________________________________________________________________________
//
void TRExFit::AddHistoPath(const std::string& path){
    fHistoPaths.push_back(path);
}

//__________________________________________________________________________________
// apply smoothing to systematics
void TRExFit::SmoothSystematics(std::string syst){
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Smoothing and/or Symmetrising Systematic Variations ...\n";

    for(std::size_t i_ch = 0; i_ch< fRegions.size(); ++i_ch){
        //
        // Scale systematics according to smoothing of another region (inter-region smoothing)
        // Loop on previous regions
        for(std::size_t j_ch=0; j_ch < i_ch; ++j_ch){
            if(fRegions[i_ch]->fIsBinOfRegion[fRegions[j_ch]->fName] <=0) continue; // NB: these bins have to start with 1, not 0 !! 0 means not filled (due to map implementation)
            const int binIdx = fRegions[i_ch]->fIsBinOfRegion[fRegions[j_ch]->fName];
            for(const auto& sh : fRegions[i_ch]->fSampleHists){
                for(const auto& syh : sh->fSyst){
                    float scaleUp = 1.;
                    float scaleDown = 1.;
                    // get scale factors to apply according to reference bin of reference region
                    if(fRegions[j_ch]->GetSampleHist(sh->fSample->fName)!=nullptr){
                        std::shared_ptr<SampleHist> sh_ref = fRegions[j_ch]->GetSampleHist(sh->fSample->fName);
                        std::shared_ptr<SystematicHist> syh_ref = sh_ref->GetSystematic(syh->fSystematic->fName);
                        if(syh_ref!=nullptr){
                            float systVarUp_orig = syh_ref->fHistUp_orig->GetBinContent(binIdx);
                            float systVarDown_orig = syh_ref->fHistDown_orig->GetBinContent(binIdx);
                            float systVarUp = syh_ref->fHistUp->GetBinContent(binIdx);
                            float systVarDown = syh_ref->fHistDown->GetBinContent(binIdx);
                            if(systVarUp_orig!=0) scaleUp = systVarUp/systVarUp_orig;
                            else LOG(WARNING) << "In inter-region smoothing attempting to divide by zero. Skipping scaling.\n";
                            if(systVarDown_orig!=0) scaleDown = systVarDown/systVarDown_orig;
                            else LOG(WARNING) << "In inter-region smoothing attempting to divide by zero. Skipping scaling.\n";
                            if(syh->fSystematic->fSymmetrisationType==HistoTools::SYMMETRIZEONESIDED){
                                bool isUp = HistoTools::Separation(sh->fHist.get(),syh->fHistUp.get()) >= HistoTools::Separation(sh->fHist.get(),syh->fHistDown.get());
                                if(isUp) scaleDown = 1.;
                                else     scaleUp   = 1.;
                            }
                        }
                    }
                    else{
                        LOG(WARNING) << "Sample not found in region indicated as reference for inter-region smoothing.\n";
                    }
                    // scale
                    syh->fHistUp->Scale(scaleUp);
                    syh->fHistDown->Scale(scaleDown);
                }
            }
        }

        // collect information which systematics contain reference smoothing samples
        std::vector<std::string> referenceSmoothSysts{};
        for (const auto& isyst : fSystematics){
            if (std::find(isyst->fRegions.begin(), isyst->fRegions.end(), fRegions[i_ch]->fName) == isyst->fRegions.end()) continue;
            if (isyst->fReferenceSmoothing != ""){
                referenceSmoothSysts.emplace_back(isyst->fName);
            }
        }

        // if there are no reference smoothing samples, proceed as usual
        if (referenceSmoothSysts.size() == 0){
            for(auto& isample : fRegions[i_ch]->fSampleHists) {
                isample->SmoothSyst(fSmoothOption, fAlternativeShapeHistFactory, syst, false);
            }
        } else {
            std::vector<std::size_t> usedSysts{};
            for (auto& isample : fRegions[i_ch]->fSampleHists) {
                for (std::size_t i_syst = 0; i_syst < fSystematics.size(); ++i_syst){
                    if (fSystematics.at(i_syst) == nullptr) continue;
                    // check only systematics for the samples that are specified
                    if (std::find(fSystematics.at(i_syst)->fSamples.begin(), fSystematics.at(i_syst)->fSamples.end(), isample->GetSample()->fName) == fSystematics.at(i_syst)->fSamples.end()) continue;
                    // take only systematics that belong to this region
                    if (std::find(fSystematics.at(i_syst)->fRegions.begin(), fSystematics.at(i_syst)->fRegions.end(), fRegions[i_ch]->fName) == fSystematics.at(i_syst)->fRegions.end()) continue;
                    if (fSystematics.at(i_syst)->fReferenceSmoothing == "") {
                        // the systemtic is not using special smoothing
                        isample->SmoothSyst(fSmoothOption, fAlternativeShapeHistFactory, fSystematics.at(i_syst)->fName, true);
                    } else {
                        // check if the syst has been smoothed already
                        if (std::find(usedSysts.begin(), usedSysts.end(), i_syst) != usedSysts.end()) continue;
                        // Need to apply special smoothing
                        // smooth the reference sample
                        std::shared_ptr<SampleHist> sh = GetSampleHistFromName(fRegions[i_ch].get(), fSystematics.at(i_syst)->fReferenceSmoothing);
                            if (sh == nullptr){
                            LOG(ERROR) << "Cannot find ReferenceSmoothing in the list of samples!\n";
                            exit(EXIT_FAILURE);
                        }

                        std::unique_ptr<TH1> nominal_cpy = nullptr;
                        std::unique_ptr<TH1> up_cpy = nullptr;
                        std::unique_ptr<TH1> down_cpy = nullptr;

                        int systIndex = -1;
                        // smooth on the sample that is specified in ReferenceSmoothing
                        for (auto& jsample : fRegions[i_ch]->fSampleHists) {
                            if (jsample->GetSample()->fName == fSystematics.at(i_syst)->fReferenceSmoothing){
                                sh->SmoothSyst(fSmoothOption, fAlternativeShapeHistFactory, fSystematics.at(i_syst)->fName, true);

                                // save the smoothed histograms
                                nominal_cpy = std::unique_ptr<TH1>(static_cast<TH1*>(jsample->fHist->Clone()));
                                systIndex = GetSystIndex(jsample.get(), fSystematics.at(i_syst)->fName);
                                if (systIndex < 0){
                                    LOG(WARNING) << "Cannot find systematic in the list wont smooth!\n";
                                    return;
                                }
                                up_cpy = std::unique_ptr<TH1>(static_cast<TH1*>(jsample->fSyst[systIndex]->fHistUp->Clone()));
                                down_cpy = std::unique_ptr<TH1>(static_cast<TH1*>(jsample->fSyst[systIndex]->fHistDown->Clone()));
                                break;
                            }
                        }

                        // finally, apply the same smoothing to all other samples, bin-by-bin
                        for (auto& jsample : fRegions[i_ch]->fSampleHists) {
                            // skip samples that do not belong to this systematics
                            if (std::find(fSystematics.at(i_syst)->fSamples.begin(), fSystematics.at(i_syst)->fSamples.end(), jsample->GetSample()->fName) ==
                                fSystematics.at(i_syst)->fSamples.end()) continue;
                            // skip the one that has already been smoothed, the ReferenceSmoothing
                            if (jsample->GetSample()->fName == fSystematics.at(i_syst)->fReferenceSmoothing) continue;

                            if (systIndex < 0){
                                LOG(WARNING) << "Cannot find systematic in the list wont smooth!\n";
                                return;
                            }
                            jsample->fSyst[systIndex]->fHistUp.reset(CopySmoothedHisto(jsample.get(),nominal_cpy.get(),up_cpy.get(),down_cpy.get(),true));
                            jsample->fSyst[systIndex]->fHistDown.reset(CopySmoothedHisto(jsample.get(),nominal_cpy.get(),up_cpy.get(),down_cpy.get(),false));
                        }

                        usedSysts.emplace_back(i_syst);
                    }
                } // loop over systs
            } // loop over samples
        }
    }
}

//
// Try to split root file creation and histogram wiriting
//__________________________________________________________________________________
// create new root file(s)
void TRExFit::CreateRootFiles(){
    bool recreate = !fUpdate;
    gSystem->mkdir( fName.c_str());
    gSystem->mkdir( (fName + "/Histograms/").c_str() );
    std::string fileName;
    bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    //
    if(singleOutputFile){
        if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_histos" + fSaveSuffix + ".root";
        else                 fileName = fName + "/Histograms/" + fInputName + "_histos" + fSaveSuffix + ".root";
        // Bootstrap
        if(fBootstrap!="" && fBootstrapIdx>=0){
            if (fBootstrapSyst!=""){
                fileName = Common::ReplaceString(fileName,Form("_histos%s.root", fSaveSuffix.c_str()),Form("_histos%s_%s_%d.root",fSaveSuffix.c_str(),fBootstrapSyst.c_str(),fBootstrapIdx));
            } else{
                fileName = Common::ReplaceString(fileName,Form("_histos%s.root", fSaveSuffix.c_str()),Form("_histos%s_%s_%d.root",fSaveSuffix.c_str(),fBootstrapSample.c_str(),fBootstrapIdx));
            }
        }
        LOG(INFO) << "-------------------------------------------\n";
        LOG(INFO) << "Creating/updating file " << fileName << " ...\n";
        if(recreate) fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"RECREATE")));
        else         fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"UPDATE")));
        // we need to replace the pointer since if the file was opened with a READ option, we would not
        // be able to write into that file
        auto itr = TRExFitter::TFILEMAP.find(fileName);
        if (itr != TRExFitter::TFILEMAP.end()) {
            itr->second.reset(std::move(fFiles.back()));
        } else {
            TRExFitter::TFILEMAP.insert(std::make_pair(fileName,std::move(fFiles.back())));
        }
    }
    else{
        for(const auto& ireg : fRegions) {
            if(fInputFolder!="") fileName = fInputFolder           + fInputName + "_" + ireg->fName + "_histos" + fSaveSuffix + ".root";
            else                 fileName = fName + "/Histograms/" + fInputName + "_" + ireg->fName + "_histos" + fSaveSuffix + ".root";
            // Bootstrap
            if(fBootstrap!="" && fBootstrapIdx>=0){
                if (fBootstrapSyst!=""){
                    fileName = Common::ReplaceString(fileName,Form("_histos%s.root", fSaveSuffix.c_str()),Form("_histos%s_%s_%d.root",fSaveSuffix.c_str(),fBootstrapSyst.c_str(),fBootstrapIdx));
                } else{
                    fileName = Common::ReplaceString(fileName,Form("_histos%s.root", fSaveSuffix.c_str()),Form("_histos%s_%s_%d.root",fSaveSuffix.c_str(),fBootstrapSample.c_str(),fBootstrapIdx));
                }
            }
            LOG(INFO) << "-------------------------------------------\n";
            LOG(INFO) << "Creating/updating file " << fileName << " ...\n";
            if(recreate) fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"RECREATE")));
            else         fFiles.emplace_back(std::move(TFile::Open(fileName.c_str(),"UPDATE")));
            // we need to replace the pointer since if the file was opened with a READ option, we would not
            // be able to write into that file
            auto itr = TRExFitter::TFILEMAP.find(fileName);
            if (itr != TRExFitter::TFILEMAP.end()) {
                itr->second.reset(std::move(fFiles.back()));
            } else {
                TRExFitter::TFILEMAP.insert(std::make_pair(fileName,std::move(fFiles.back())));
            }
        }
    }
}

//__________________________________________________________________________________
// fill files with all the histograms
void TRExFit::WriteHistos(bool reWriteOrig) const{
    bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    std::string fileName;
    for(std::size_t i_ch = 0; i_ch < fRegions.size(); ++i_ch) {
        //
        if(singleOutputFile) fileName = fFiles[0]->GetName();
        else                 fileName = fFiles[i_ch]->GetName();

        LOG(INFO) << "-------------------------------------------\n";
        LOG(INFO) << "Writing histograms to file " << fileName << " for region: " << fRegions.at(i_ch)->fName << " ...\n";

        // find data and get scales per bin if data are weighted
        std::vector<double> binScales;
        if (fDataWeighted) {
            for (const auto& isample : fSamples) {
                if (isample->fType != Sample::SampleType::DATA) continue;
                const std::shared_ptr<SampleHist>& sh = fRegions[i_ch]->GetSampleHist(isample->fName);
                binScales = sh->GetDataScales();
                break;
            }
        }

        std::vector<double> binPruningReference{};

        if (fPruningType == TRExFit::BACKGROUNDREFERENCE || fPruningType == TRExFit::COMBINEDREFERENCE || fPruningType == TRExFit::COMBINEDSIGNAL) {
            bool isFirst(true);
            for (const auto& isample : fSamples) {
                if (isample->fType == Sample::SampleType::DATA) continue;
                if (isample->fType == Sample::SampleType::GHOST) continue;
                if (isample->fType == Sample::SampleType::EFT) continue;
                if ((fPruningType == TRExFit::BACKGROUNDREFERENCE) && (isample->fType == Sample::SampleType::SIGNAL)) continue;
                if ((fPruningType == TRExFit::COMBINEDSIGNAL) && (isample->fType == Sample::SampleType::BACKGROUND)) continue;
                std::shared_ptr<SampleHist> sh = fRegions[i_ch]->GetSampleHist(isample->fName);
                if (!sh) {
                    LOG(DEBUG) << "Skipping pruning for not-found sample " << isample->fName << " in region " << fRegions[i_ch]->fName << "\n";
                    continue;
                }
                for (int ibin = 1; ibin <= sh->fHist->GetNbinsX(); ++ibin) {
                    if (isFirst) binPruningReference.emplace_back(sh->fHist->GetBinContent(ibin));
                    else         binPruningReference.at(ibin-1) += sh->fHist->GetBinContent(ibin);
                }
                isFirst = false;
            }
        }

        // Get the groups of samples that share gammas
        std::vector<std::vector<std::string> > gammaCorrelationSampleNames;
        for (const auto& isample : fSamples) {
            if (isample->fCorrelateGammasWithSample == "") continue;
            //check if the given sample and the sample it is correlated with are in the list
            bool found(false);
            for (std::size_t i = 0; i < gammaCorrelationSampleNames.size(); ++i) {
                auto it1 = std::find(gammaCorrelationSampleNames.at(i).begin(), gammaCorrelationSampleNames.at(i).end(), isample->fCorrelateGammasWithSample);
                auto it2 = std::find(gammaCorrelationSampleNames.at(i).begin(), gammaCorrelationSampleNames.at(i).end(), isample->fName);
                if ((it1 != gammaCorrelationSampleNames.at(i).end()) || (it2 != gammaCorrelationSampleNames.at(i).end())) {
                    found = true;
                    if (it1 == gammaCorrelationSampleNames.at(i).end()) {
                        gammaCorrelationSampleNames.at(i).emplace_back(isample->fCorrelateGammasWithSample);
                    }
                    if (it2 == gammaCorrelationSampleNames.at(i).end()) {
                        gammaCorrelationSampleNames.at(i).emplace_back(isample->fName);
                    }
                    break;
                }
            }
            if (!found) {
                std::vector<std::string> tmp;
                tmp.emplace_back(isample->fCorrelateGammasWithSample);
                tmp.emplace_back(isample->fName);
                gammaCorrelationSampleNames.emplace_back(tmp);
            }
        }

        std::vector<std::vector<const TH1*> > corrHistos(gammaCorrelationSampleNames.size());

        // get the histograms for each group
        for (const auto& isample : fSamples) {
            const std::string name = isample->fName;
            std::size_t index(0);
            for (const auto& icor : gammaCorrelationSampleNames) {
                auto itr = std::find(icor.begin(), icor.end(), name);
                if (itr != icor.end()) {
                    // doesnt exist for this region
                    if ((Common::FindInStringVector(isample->fRegions, fRegions[i_ch]->fName) < 0) && (isample->fRegions.at(0) != "all")) {
                        corrHistos.at(index).emplace_back(nullptr);
                        break;
                    }
                    const auto& sh = fRegions[i_ch]->GetSampleHist(name);
                    if (!sh) {
                        LOG(ERROR) << "Cannot find SampleHist: " << name << " in region " << fRegions[i_ch]->fName << "\n";
                        exit(EXIT_FAILURE);
                    }
                    const TH1* h = sh->fHist.get();
                    corrHistos.at(index).emplace_back(h);
                    break;
                }
                ++index;
            }
        }

        // for each correlation group, get the correction per sample per bin
        std::vector<std::vector<std::vector<double> > > errorCor;
        for (const auto& igroup : corrHistos) {
            errorCor.emplace_back(SampleHist::GetGammaCorrelationUncertainties(igroup));
        }

        for (std::size_t i_smp = 0; i_smp < fSamples.size(); ++i_smp) {
            std::shared_ptr<SampleHist> sh = fRegions[i_ch]->GetSampleHist(fSamples[i_smp]->fName);
            if(!sh){
                LOG(DEBUG) << "SampleHist[" << i_smp << "] for sample " << fSamples[i_smp]->fName << " not there.\n";
                continue;
            }
            std::vector<double> corr;

            // for correlations, get the proper corrections
            if (!gammaCorrelationSampleNames.empty() && fApplyGammaCorrection) {
                std::size_t groupIndex(0);
                std::size_t vectorIndex(0);
                bool found(false);
                for (const auto& icor : gammaCorrelationSampleNames) {
                    auto itr = std::find(icor.begin(), icor.end(), fSamples[i_smp]->fName);
                    if (itr != icor.end()) {
                        vectorIndex = std::distance(icor.begin(), itr);
                        found = true;
                        break;
                    }
                    ++groupIndex;
                }
                if (found && !errorCor.at(groupIndex).empty()) corr = errorCor.at(groupIndex).at(vectorIndex);
            }
            // set file and histo names for nominal
            sh->fHistoName = sh->fHist->GetName();
            sh->fFileName = fileName;
            // set file and histo names for systematics
            for (const auto& isyst : sh->fSyst) {
                if(isyst->fHistUp  ==nullptr) continue;
                if(isyst->fHistDown==nullptr) continue;
                isyst->fFileNameUp    = fileName;
                isyst->fHistoNameUp   = isyst->fHistUp->GetName();
                isyst->fFileNameDown  = fileName;
                isyst->fHistoNameDown = isyst->fHistDown->GetName();
                if (isyst->fHasShape) {
                    isyst->fFileNameShapeUp    = fileName;
                    isyst->fHistoNameShapeUp   = isyst->fHistShapeUp->GetName();
                    isyst->fFileNameShapeDown  = fileName;
                    isyst->fHistoNameShapeDown = isyst->fHistShapeDown->GetName();
                }
            }
            const std::vector<int>& droppedBins = fRegions[i_ch]->GetAutomaticDropBins() ?
                                                  fRegions[i_ch]->fComputedBlindedBins : fRegions[i_ch]->fDropBins;
            const double threshold = fUseStatErr ? fStatErrThres : -1;

            if(singleOutputFile) sh->WriteToFile(droppedBins, binScales, threshold, binPruningReference, fRegions[i_ch]->fName, fUseRebinned, fFiles[0]   , reWriteOrig, corr);
            else                 sh->WriteToFile(droppedBins, binScales, threshold, binPruningReference, fRegions[i_ch]->fName, fUseRebinned, fFiles[i_ch], reWriteOrig, corr);
        }
    }
    LOG(INFO) << "-------------------------------------------\n";
}

//__________________________________________________________________________________
// Draw morphing plots
void TRExFit::DrawMorphingPlots(const std::string& name) const{
    for(const auto& reg : fRegions){
        TCanvas c("c","c",600,600);
        TPad p0("p0","p0",0,0.35,1,1);
        TPad p1("p1","p1",0,0,1,0.35);
        p0.SetBottomMargin(0);
        p1.SetTopMargin(0);
        p1.SetBottomMargin(0.3);
        p0.Draw();
        p1.Draw();
        p0.cd();
        int nTemp = 0;
        std::vector<std::unique_ptr<TH1> > hVec;
        // cppcheck-suppress variableScope
        std::vector<std::unique_ptr<TH1> > hVecRatio;
        for(const auto& sh : reg->fSampleHists){
            Sample* smp = sh->fSample;
            // if the sample has morphing
            if(smp->fIsMorph[name]){
                std::unique_ptr<TH1> h(static_cast<TH1*>(sh->fHist->Clone(("h_temp_"+smp->fName).c_str())));
                if(h->GetFillColor()!=0) h->SetLineColor(h->GetFillColor());
                h->SetFillStyle(0);
                h->SetLineWidth(2);
                h->Scale(1./h->Integral());
                if(nTemp==0) h->Draw("HIST");
                else         h->Draw("HIST same");
                hVec.push_back(std::move(h));
                nTemp++;
            }
        }
        if(hVec.size()>0){
            hVec[0]->SetMaximum(1.25*hVec[0]->GetMaximum());
            hVec[0]->GetYaxis()->SetTitle("Fraction of events");
            hVec[0]->GetYaxis()->SetTitleOffset(1.75);
            // ratio
            p1.cd();
            for(const auto& hh : hVec){
                hVecRatio.push_back(std::unique_ptr<TH1>(static_cast<TH1*>(hh->Clone())));
                hVecRatio[hVecRatio.size()-1]->Divide(hVec[0].get());
                if(hVecRatio.size()-1==0) hVecRatio[hVecRatio.size()-1]->Draw("HIST");
                else                      hVecRatio[hVecRatio.size()-1]->Draw("HIST same");
            }
            hVecRatio[0]->SetMinimum(0.91);
            hVecRatio[0]->SetMaximum(1.09);
            hVecRatio[0]->GetXaxis()->SetTitle(reg->fVariableTitle.c_str());
            hVecRatio[0]->GetYaxis()->SetTitle("Ratio");
            hVecRatio[0]->GetYaxis()->SetTitleOffset(1.75);
            hVecRatio[0]->GetXaxis()->SetTitleOffset(3);
            Common::SaveCanvasAs(c, fName+"/Morphing/Templates_"+name+"_"+reg->fName);
        }
    }
}

//__________________________________________________________________________________
// Draw syst plots
void TRExFit::DrawSystPlots() const{
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Systematics").c_str());
    for(const auto& ireg : fRegions) {
        LOG(INFO) << "Drawing red-blue plots for region: " << ireg->fName << "\n";
        for(const auto& isample : ireg->fSampleHists) {
            isample->DrawSystPlot(nullptr, false, false, fFitType == TRExFit::FitType::UNFOLDING);
        }
    }
}

//__________________________________________________________________________________
// Draw syst plots for combined samples
void TRExFit::DrawSystPlotsSumSamples() const{
    std::unique_ptr<TH1> h_dataCopy(nullptr);
    for(const auto& reg : fRegions){
        LOG(INFO) << "-------------------------------------------\n";
        LOG(INFO) << "Drawing combined plots of syst effects on data for region: " << reg->fName << "...\n";
        SampleHist hist{};
        bool empty = true;
        std::set<std::string> systNames;
        for(const auto& isample : reg->fSampleHists) {
            for(std::size_t i_smSyst=0; i_smSyst < isample->fSyst.size(); i_smSyst++){
                systNames.insert(isample->fSyst[i_smSyst]->fName);
            }
        }
        for(const auto& isample : reg->fSampleHists){
            if(isample->fSample->fType==Sample::SampleType::DATA) h_dataCopy=std::unique_ptr<TH1>(static_cast<TH1*>(isample->fHist->Clone()));
            else if(isample->fSample->fType==Sample::SampleType::GHOST) continue;
            else if(isample->fSample->fType==Sample::SampleType::EFT) continue;
            else {
                if(empty){
                    hist.CloneSampleHist(isample.get(),systNames, 1);
                    hist.fName = reg->fName + "_Combined";
                    empty=false;
                } else {
                    hist.SampleHistAdd(isample.get(), 1);
                }
            }
        }
        // Get blinded bins
        hist.DrawSystPlot(h_dataCopy.get(), true, fSystDataPlot_upFrame, false, reg->fBlindedBins);
    }
}

//__________________________________________________________________________________
// this method creates a map for all the EFT SM samples per region
void TRExFit::CreateEFTNominalSampleMaps(){
    if (fFitType != TRExFit::FitType::EFT) return;
    if (!fEFTConfig.GetSplitSamplesPerBin()) return;
    for (const auto& reg : fRegions){
        std::map< std::string, std::vector<std::shared_ptr<SampleHist> > > regionSampleMaps;
        if (fEFTNominalSampleMaps.find(reg->fName) != fEFTNominalSampleMaps.end()) {
            regionSampleMaps = fEFTNominalSampleMaps[reg->fName];
        }
        for(const auto& smp : reg->fSampleHists){
            LOG(DEBUG) << "  hist " << smp->fName << " in region " << reg->fName << "\n";
            if(smp->fName.rfind("SM_", 0) == 0 && smp->fName.find("_bin") != std::string::npos){
                LOG(DEBUG) << "  found hist " << smp->fName << " in region " << reg->fName << "\n";

                size_t found = smp->fName.find("_bin");
                //int binnumber = Common::convertStoNum<int>(smp->fName.substr(found+4));
                std::string sample_name = smp->fName.substr(0,found);
                if (regionSampleMaps.find(sample_name) == regionSampleMaps.end()) {
                    regionSampleMaps[sample_name] = {smp};
                }
                else{
                    regionSampleMaps[sample_name].emplace_back(smp);
                }
            }
        }
        fEFTNominalSampleMaps[reg->fName] = regionSampleMaps;
    }
}

//__________________________________________________________________________________
// this method takes care of rebinning, smoothing, fixing
void TRExFit::CorrectHistograms(){
    //
    // loop on regions, and then perform a set of operations for each of them
    for(auto& reg : fRegions){
        int nbins(0);
        //
        // 1. Reset histograms to the ones save as "_orig" (both for nominal and systematics
        for(const auto& smp : fSamples){
            //
            // eventually skip sample / region combination
            if(Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            int fillcolor = sh->fHist->GetFillColor();
            int linecolor = sh->fHist->GetLineColor();
            TH1* h_orig = sh->fHist_orig.get();
            TH1* h = nullptr;
            if(h_orig!=nullptr) h = static_cast<TH1*>(h_orig->Clone(sh->fHist->GetName()));
            sh->fHist = std::unique_ptr<TH1>(h);
            if(sh->fHist==nullptr) continue;
            sh->fHist->SetLineColor(linecolor);
            sh->fHist->SetFillColor(fillcolor);

            // Scale MC stat
            Common::ScaleMCstatInHist(sh->fHist.get(), smp->fMCstatScale);

            // loop on systematics
            for(const auto& syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                // skip also separate gamma systs
                if(syst->fName.find("stat_")!=std::string::npos) continue;
                //
                // get the original syst histograms & reset the syst histograms
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic( syst->fName );
                //
                if(syh==nullptr) continue;
                TH1* hUp_orig   = syh->fHistUp_orig.get();
                TH1* hDown_orig = syh->fHistDown_orig.get();
                //
                // if Overall only => fill SystematicHist
                if(syst->fType==Systematic::OVERALL){
                    for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
                        if(hUp_orig!=nullptr)   hUp_orig  ->SetBinContent(i_bin,h_orig->GetBinContent(i_bin)*(1.+syst->fOverallUp));
                        if(hDown_orig!=nullptr) hDown_orig->SetBinContent(i_bin,h_orig->GetBinContent(i_bin)*(1.+syst->fOverallDown));
                    }
                }
                // cppcheck-suppress constVariablePointer
                TH1* hUp   = nullptr;
                // cppcheck-suppress constVariablePointer
                TH1* hDown = nullptr;
                if(hUp_orig!=nullptr)   hUp   = static_cast<TH1*>(hUp_orig->Clone( syh->fHistUp->GetName()));
                if(hDown_orig!=nullptr) hDown = static_cast<TH1*>(hDown_orig->Clone(syh->fHistDown->GetName()));
                syh->fHistUp.reset(hUp);
                syh->fHistDown.reset(hDown);
            }
        }
        //
        // 2. Rebin
        for(const auto& smp : fSamples){
            //
            // eventually skip sample / region combination
            if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            //
            // Rebinning (FIXME: better to introduce a method Region::Rebin() ?)
            if(reg->fHistoNBinsRebinPost>0){
                LOG(DEBUG) << "Rebinning " << smp->fName << " to " << reg->fHistoNBinsRebinPost << " bins.\n";
                sh->fHist = std::unique_ptr<TH1>(sh->fHist->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]));
                if (TRExFitter::MERGEUNDEROVERFLOW && sh->fHist) Common::MergeUnderOverFlow(sh->fHist.get());
                for(auto& syh : sh->fSyst){
                    LOG(DEBUG) << "  systematic " << syh->fName << " to " << reg->fHistoNBinsRebinPost << " bins.\n";
                    if(syh==nullptr) continue;
                    if(syh->fSystematic->fSampleUp==""   && syh->fSystematic->fHasUpVariation   && syh->fHistUp!=nullptr) {
                        syh->fHistUp.reset(syh->fHistUp  ->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]));
                    } else {
                        const std::string name = Form("h_%s_%s_%sUp",reg->fName.c_str(), smp->fName.c_str(), syh->fSystematic->fStoredName.c_str());
                        syh->fHistUp.reset(static_cast<TH1*>(sh->fHist->Clone(name.c_str())));
                    }
                    if(syh->fSystematic->fSampleDown=="" && syh->fSystematic->fHasDownVariation && syh->fHistDown!=nullptr) {
                        syh->fHistDown.reset(syh->fHistDown->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]));
                    } else {
                        const std::string name = Form("h_%s_%s_%sDown",reg->fName.c_str(), smp->fName.c_str(), syh->fSystematic->fStoredName.c_str());
                        syh->fHistDown.reset(static_cast<TH1*>(sh->fHist->Clone(name.c_str())));
                    }
                    if (TRExFitter::MERGEUNDEROVERFLOW && syh->fHistUp)   Common::MergeUnderOverFlow(syh->fHistUp.get());
                    if (TRExFitter::MERGEUNDEROVERFLOW && syh->fHistDown) Common::MergeUnderOverFlow(syh->fHistDown.get());
                }
                //
                // rebin also separate-gamma hists!
                if(smp->fSeparateGammas){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic( "stat_"+smp->fName );
                    if(syh==nullptr) continue;
                    if(syh->fHistUp!=nullptr)   syh->fHistUp  ->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]);
                    if(syh->fHistDown!=nullptr) syh->fHistDown->Rebin(reg->fHistoNBinsRebinPost,"",&reg->fHistoBinsPost[0]);
                    if (TRExFitter::MERGEUNDEROVERFLOW && syh->fHistUp)   Common::MergeUnderOverFlow(syh->fHistUp.get());
                    if (TRExFitter::MERGEUNDEROVERFLOW && syh->fHistDown) Common::MergeUnderOverFlow(syh->fHistDown.get());
                }
            }
            nbins = sh->fHist->GetNbinsX();
        }
        // Randomize MC (before add/multiply/scale)
        if(TRExFitter::OPTION["RandomizeMC"]!=0){
            gRandom->SetSeed(TRExFitter::OPTION["RandomizeMC"]);
            for(const auto& smp : fSamples){
                //
                // eventually skip sample / region combination
                if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                //
                std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                if(sh==nullptr) continue;
                if(sh->fHist==nullptr) continue;
                //
                if(smp->fUseMCStat){
                    TH1* hTmp = sh->fHist.get();
                    for(int i_bin=1;i_bin<=hTmp->GetNbinsX();i_bin++){
                        hTmp->SetBinContent(i_bin,gRandom->Poisson( hTmp->GetBinContent(i_bin) ));
                    }
                    for(const auto& syh : sh->fSyst){
                        for(int i_ud=0;i_ud<2;i_ud++){
                            if(i_ud==0) hTmp = syh->fHistUp.get();
                            else        hTmp = syh->fHistDown.get();
                            for(int i_bin=1;i_bin<=hTmp->GetNbinsX();i_bin++){
                                hTmp->SetBinContent(i_bin,gRandom->Poisson( hTmp->GetBinContent(i_bin) ));
                            }
                        }
                    }
                }
            }
        }

        // 3. Add/Multiply/Scale
        for(auto smp : fSamples){
            //
            // eventually skip sample / region combination
            if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
            //
            std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
            if(sh==nullptr) continue;
            if(sh->fHist==nullptr) continue;
            int fillcolor = sh->fHist->GetFillColor();
            int linecolor = sh->fHist->GetLineColor();
            //
            // Subtraction / Addition of sample
            for(const auto& sample : smp->fSubtractSamples){
                LOG(DEBUG) << "Subtracting sample " << sample << " from sample " << smp->fName << "\n";
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(sample);
                if(smph0!=nullptr) sh->Add(smph0.get(),-1);
                else LOG(WARNING) << "Sample Hist of sample " << sample << " not found ...\n";
            }
            for(const auto& sample : smp->fAddSamples){
                LOG(DEBUG) << "Adding sample " << sample << " to sample " << smp->fName << "\n";
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(sample);
                if(smph0!=nullptr) sh->Add(smph0.get());
                else LOG(WARNING) << "Sample Hist of sample " << sample << " not found ...\n";
            }
            // Division & Multiplication by other samples
            if(smp->fMultiplyBy!=""){
                LOG(DEBUG) << "Multiplying " << smp->fName << " by sample " << smp->fMultiplyBy << "\n";
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(smp->fMultiplyBy);
                if(smph0!=nullptr) sh->Multiply(smph0.get());
                else LOG(WARNING) << "Sample Hist of sample " << smp->fMultiplyBy << " not found ...\n";
            }
            if(smp->fDivideBy!=""){
                LOG(DEBUG) << "Dividing " << smp->fName << " by sample " << smp->fDivideBy << " from sample " << smp->fName << "\n";
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(smp->fDivideBy);
                if(smph0!=nullptr) sh->Divide(smph0.get());
                else LOG(WARNING) << "Sample Hist of sample " << smp->fDivideBy << " not found ...\n";
            }
            // Norm to sample
            if(smp->fNormToSample!=""){
                LOG(DEBUG) << "Normalizing " << smp->fName << " to sample " << smp->fNormToSample << "\n";
                std::shared_ptr<SampleHist> smph0 = reg->GetSampleHist(smp->fNormToSample);
                if(smph0!=nullptr) sh->Scale(smph0->fHist->Integral()/sh->fHist->Integral());
                else LOG(WARNING) << "Sample Hist of sample " << smp->fNormToSample << " not found ...\n";
            }

            //
            // For SampleUp / SampleDown
            for(auto& syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                // get the original syst histograms & reset the syst histograms
                std::shared_ptr<SystematicHist> syh;
                //
                // if syst defined with SampleUp / SampleDown
                if( syst->fSampleUp != "" || syst->fSampleDown != "" ){
                    LOG(DEBUG) << "SampleUp/SampleDown set for systematic " << syst->fName << ".\n";
                    bool isDummy = ( syst->fDummyForSamples.size()>0 && Common::FindInStringVector(syst->fDummyForSamples,smp->fName)>=0 );
                    std::unique_ptr<TH1> h_up = nullptr;
                    if(syst->fSampleUp   !="" && !isDummy){
                        if(reg->GetSampleHist(syst->fSampleUp  )){
                            h_up.reset(static_cast<TH1*>(reg->GetSampleHist(syst->fSampleUp  )->fHist->Clone("h_tmp_up")));
                        }
                    }
                    else{
                        h_up.reset(static_cast<TH1*>(sh->fHist->Clone("h_tmp_up")));
                    }
                    std::unique_ptr<TH1> h_down = nullptr;
                    if(syst->fSampleDown !="" && !isDummy){
                        if(reg->GetSampleHist(syst->fSampleDown)){
                            h_down.reset(static_cast<TH1*>(reg->GetSampleHist(syst->fSampleDown)->fHist.get()->Clone("h_tmp_down")));
                        }
                    }
                    else{
                        h_down.reset(static_cast<TH1*>(sh->fHist.get()->Clone("h_tmp_down")));
                    }
                    //
                    // if systematic also uses ReferenceSample, produce syst variations according to the refefence sample instead of nominal
                    if(syst->fReferenceSample!=""){
                        LOG(DEBUG) << "ReferenceSample set for a systematic with SampleUp/SampleDown. Building proper systematic variation.\n";
                        std::shared_ptr<SampleHist> refSh = reg->GetSampleHist(syst->fReferenceSample);
                        if(refSh!=nullptr){
                            if(syst->fSampleUp != ""){
                                h_up->Divide(refSh->fHist.get());
                                h_up->Multiply(sh->fHist.get());
                            }
                            if(syst->fSampleDown != ""){
                                h_down->Divide(refSh->fHist.get());
                                h_down->Multiply(sh->fHist.get());
                            }
                        }
                        else{
                            LOG(WARNING) << "ReferenceSample for systematc " << syst->fName << " set but no corresponding sample found. Ignoring.\n";
                        }
                    }
                    syh = sh->AddHistoSyst(syst->fName,syst->fStoredName,h_up.get(),h_down.get());
                    syh->fSystematic = syst;
                }
            }

            //
            // Save to _preSmooth histograms (to be shown in syst plots) at this point
            sh->fHist_preSmooth.reset(static_cast<TH1*>(sh->fHist->Clone(Form("%s_preSmooth",sh->fHist->GetName()))));
            sh->fHist_preSmooth->SetDirectory(nullptr);
            for(const auto& syh : sh->fSyst){
                if(syh!=nullptr){
                    if(syh->fHistUp!=nullptr)   syh->fHistUp_preSmooth.reset(static_cast<TH1*>(syh->fHistUp->Clone(Form("%s_preSmooth",syh->fHistUp->GetName()))));
                    else                        syh->fHistUp_preSmooth.reset(static_cast<TH1*>(sh->fHist_preSmooth->Clone()));
                    syh->fHistUp_preSmooth->SetDirectory(nullptr);
                    if(syh->fHistDown!=nullptr) syh->fHistDown_preSmooth.reset(static_cast<TH1*>(syh->fHistDown->Clone(Form("%s_preSmooth",syh->fHistDown->GetName()))));
                    else                        syh->fHistDown_preSmooth.reset(static_cast<TH1*>(sh->fHist_preSmooth->Clone()));
                    syh->fHistDown_preSmooth->SetDirectory(nullptr);
                }
            }

            //
            // Fix empty bins
            if(smp->fType!=Sample::SampleType::DATA && smp->fType!=Sample::SampleType::SIGNAL){
                sh->FixEmptyBins(fSuppressNegativeBinWarnings);
            }

            //
            // Eventually smooth nominal histogram  (use with caution...)
            TH1* h_correction = nullptr;
            bool isFlat = false;
            if(smp->fSmooth && !reg->fSkipSmoothing){
                h_correction = static_cast<TH1*>(sh->fHist->Clone( Form("%s_corr",sh->fHist->GetName())));
                TH1* h0 = static_cast<TH1*>(sh->fHist->Clone(Form("%s_orig0",sh->fHist->GetName())));
                if (fSmoothOption == HistoTools::TTBARRESONANCE) {
                    isFlat = false;
                    Common::SmoothHistogramTtres(sh->fHist.get());
                } else {
                    isFlat = Common::SmoothHistogram(sh->fHist.get());
                }
                h_correction->Divide(h0);
            }

            //
            // Systematics
            for(const auto& syst : smp->fSystematics){
                if(syst==nullptr) continue;
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                //
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic( syst->fName );
                if(syh==nullptr) continue;
                TH1* hUp   = syh->fHistUp.get();
                TH1* hDown = syh->fHistDown.get();
                //
                // if Overall only, re-create it if smoothing was applied
                if(syst->fType==Systematic::OVERALL){
                    if(h_correction!=nullptr && smp->fSmooth){
                        for(int i_bin=1;i_bin<=sh->fHist->GetNbinsX();i_bin++){
                            hUp  ->SetBinContent(i_bin,sh->fHist->GetBinContent(i_bin)*(1.+syst->fOverallUp));
                            hDown->SetBinContent(i_bin,sh->fHist->GetBinContent(i_bin)*(1.+syst->fOverallDown));
                        }
                    }
                    continue;
                }
                //
                // correct according to the sample nominal smoothing
                if(h_correction!=nullptr && smp->fSmooth){
                    if(hUp!=nullptr  ) Common::SmoothHistogram( hUp  , isFlat );
                    if(hDown!=nullptr) Common::SmoothHistogram( hDown, isFlat );
                }

                //
                // Histogram smoothing, Symmetrisation, Massaging...
                if(!reg->fSkipSmoothing) syh -> fSmoothType = syst -> fSmoothType;
                else                                syh -> fSmoothType = 0;
                syh -> fSymmetrisationType = syst -> fSymmetrisationType;

            }  // end syst loop
            //
            // Histograms checking
            for(const auto& syst : smp->fSystematics){
                //
                // eventually skip systematic / region combination
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions,reg->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude,reg->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,reg->fName, smp->fName)>=0 ) continue;
                if( sh->GetSystematic( syst->fName )==nullptr ) continue;
                //
                HistoTools::CheckHistograms( reg->GetSampleHist(smp->fName)->fHist.get() /*nominal*/,
                                            sh->GetSystematic( syst->fName ).get() /*systematic*/,
                                            smp->fType!=Sample::SampleType::SIGNAL/*check bins with content=0*/,
                                            TRExFitter::HISTOCHECKCRASH /*cause crash if problem*/);
            }

            // Check on negative content in nominal histograms

            TH1* nom = reg->GetSampleHist(smp->fName)->fHist.get();
            const int NbinsNom  = nom->GetNbinsX();

            for( int iBin = 1; iBin <= NbinsNom; ++iBin ){
                double content = nom->GetBinContent(iBin);

                if ( content<0 ){
                    std::string name = nom->GetName();
                    if(TRExFitter::HISTOCHECKCRASH) {
                        LOG(ERROR) << "In histo \"" << name << "\", bin " << iBin << " has negative content ! Please check\n";
                        LOG(ERROR) << "Nominal: " << content << "\n";
                        exit(EXIT_FAILURE);
                    } else {
                        LOG(WARNING) << "In histo \"" << name << "\", bin " << iBin << " has negative content ! Please check\n";
                        LOG(WARNING) << "Nominal: " << content << "\n";
                        // Corrects the nominal
                        LOG(WARNING) << "Will set the bin content to 1e-06 pm 1e-07! Please check!\n";
                        nom -> SetBinContent(iBin,1e-06);
                        nom -> SetBinError(iBin, 1e-07);
                    }
                }

            }//loop over the bins

            // set the fill color
            sh->fHist->SetFillColor(fillcolor);
            sh->fHist->SetLineColor(linecolor);
        } // end sample loop
        auto itr = std::find(fBinningFilesToCreate.begin(), fBinningFilesToCreate.end(), reg->fName);
        reg->SetAndAddNbins(nbins, itr != fBinningFilesToCreate.end());
        //
    } // end region loop

    //
    // Morph smoothing
    if(fSmoothMorphingTemplates!=""){
        for(const auto& par : fMorphParams){
            LOG(INFO) << "Smoothing morphing templates for parameter " << par << "\n";
            if(fSmoothMorphingTemplates=="TRUE") SmoothMorphTemplates(par);
            else SmoothMorphTemplates(par,fSmoothMorphingTemplates);
            // to add: possibility to set initial values of parameters
        }
    }

    //
    // Plot Morphing templates
    if (fMorphParams.size() > 0){
        gSystem->mkdir((fName+"/Morphing").c_str());
        for(const auto& par : fMorphParams){
            DrawMorphingPlots(par);
        }
    }

    // drop normalisation part of systematic according to fDropNormIn
    for(const auto& reg : fRegions){
        std::map<std::string, std::map<std::string, std::vector<double> > > integrals_EFT;
        for(auto& sh : reg->fSampleHists){
            if(sh->fHist==nullptr) continue;
            for(auto syst : fSystematics){
                if(  Common::FindInStringVector(syst->fDropNormIn, reg->fName)>=0
                  || Common::FindInStringVector(syst->fDropNormIn, sh->fSample->fName)>=0
                  || Common::FindInStringVector(syst->fDropNormIn, "all")>=0
                  ){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    if(sh->fHist->Integral()!=0){
                        LOG(DEBUG) << "  Normalising syst " << syst->fName << " for sample " << sh->fSample->fName << "\n";
                        if(sh->fName.rfind("SM_", 0) == 0 && sh->fName.find("_bin") != std::string::npos && fFitType == TRExFit::FitType::EFT && fEFTConfig.GetSplitSamplesPerBin()) {
                            LOG(DEBUG) << "  in rfind(SM) if\n";
                            LOG(DEBUG) << "  Normalising syst " << syst->fName << " for sample " << sh->fSample->fName << "\n";

                             //have to find the full integral somehow now
                            std::string _smp = sh->fSample->fName;
                            size_t found = _smp.find("_bin");
                            std::string samplename = _smp.substr(0,found);

                             if (integrals_EFT.find(syst->fName) == integrals_EFT.end()) {
                                if (integrals_EFT[syst->fName].find(samplename) == integrals_EFT[syst->fName].end()) {
                                    std::map<std::string, std::vector<double> > integral_map;
                                    std::vector<double> integrals = {0,0,0};
                                    for(auto& efthist: fEFTNominalSampleMaps[reg->fName][samplename]){
                                        double err(0);
                                        const double _intNom = Common::CorrectIntegral(efthist->fHist.get(), &err);
                                        std::shared_ptr<SystematicHist> syh_backup = efthist->GetSystematic(syst->fName);
                                        TH1* _hUp = syh_backup->fHistUp.get();
                                        TH1* _hDown = syh_backup->fHistDown.get();

                                        const double _intUp = Common::CorrectIntegral(_hUp, &err);
                                        const double _intDown = Common::CorrectIntegral(_hDown, &err);

                                        integrals.at(0)+=_intNom;
                                        integrals.at(1)+=_intUp;
                                        integrals.at(2)+=_intDown;
                                    }
                                    integral_map[samplename] = integrals;
                                    integrals_EFT[syst->fName] = integral_map;
                                }
                             }
                            syh->fHistUp.get()->Scale(integrals_EFT[syst->fName][samplename].at(0)/integrals_EFT[syst->fName][samplename].at(1));
                            syh->fHistDown.get()->Scale(integrals_EFT[syst->fName][samplename].at(0)/integrals_EFT[syst->fName][samplename].at(2));
                            //Common::DropNorm(syh->fHistUp.get(), syh->fHistDown.get(), intNom);
                        }
                        else{
                            LOG(DEBUG) << "  in else\n";
                            Common::DropNorm(syh->fHistUp.get(), syh->fHistDown.get(), sh->fHist.get());
                        }
                    }
                }
            }
        }
    }

    // drop shape part of systematic according to fDropShapeIn
    for(const auto& reg : fRegions){
        for(const auto& sh : reg->fSampleHists){
            if(sh->fHist==nullptr) continue;
            for(const auto& syst : fSystematics){
                if(  Common::FindInStringVector(syst->fDropShapeIn, reg->fName)>=0
                  || Common::FindInStringVector(syst->fDropShapeIn, sh->fSample->fName)>=0
                  || Common::FindInStringVector(syst->fDropShapeIn, "all")>=0
                  ){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    LOG(DEBUG) << "  Removing shape component of syst " << syst->fName << " for sample " << sh->fSample->fName << "\n";
                    Common::DropShape(syh->fHistUp.get(), syh->fHistDown.get(), sh->fHist.get());
                }
            }
        }
    }

    //
    // Smooth systematics
    SmoothSystematics("all");

    //
    // Artifificially set all systematics not to affect overall normalisation for sample or set of samples
    // (the form should be KeepNormForSamples: ttlight+ttc+ttb,wjets
    //
    for(const auto& reg : fRegions){
        for(const auto& syst : fSystematics){
            if(syst->fKeepNormForSamples.size()==0) continue;
            for(unsigned int ii=0;ii<syst->fKeepNormForSamples.size();ii++){
                std::vector<std::string> subSamples = Common::Vectorize(syst->fKeepNormForSamples[ii],'+');
                // get nominal yield and syst yields for this sum of samples
                double yieldNominal = 0.;
                double yieldUp = 0.;
                double yieldDown = 0.;
                for(const auto& smp : fSamples){
                    if(Common::FindInStringVector(subSamples,smp->fName)<0) continue;
                    std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                    if(sh==nullptr) continue;
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    yieldNominal += sh ->fHist    ->Integral();
                    yieldUp      += syh->fHistUp  ->Integral();
                    yieldDown    += syh->fHistDown->Integral();
                }
                // scale each syst variation
                for(const auto& smp : fSamples){
                    if(Common::FindInStringVector(subSamples,smp->fName)<0) continue;
                    std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                    if(sh==nullptr) continue;
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(syst->fName);
                    if(syh==nullptr) continue;
                    LOG(DEBUG) << "Normalising syst " << syst->fName << " for sample " << smp->fName << "\n";
                    LOG(DEBUG) << "Scaling by " << (yieldNominal/yieldUp) << " (up), " << (yieldNominal/yieldDown) << " (down)\n";
                    syh->fHistUp  ->Scale(yieldNominal/yieldUp);
                    syh->fHistDown->Scale(yieldNominal/yieldDown);
                }
            }
        }
    }

    // Systematics for morphing samples inherited from nominal sample
    if (fPropagateSystsForMorphing){
        for(const auto& par : fMorphParams){
            for(const auto& reg : fRegions){
                // find nominal morphing sample Hist
                double nominalValue = 0.;
                for(const auto& norm : fNormFactors){
                    if(norm->fName==par) nominalValue = norm->GetNominal();
                }
                std::shared_ptr<SampleHist> shNominal = nullptr;
                for(const auto& sh : reg->fSampleHists){
                    if(!sh->fSample->fIsMorph[par]) continue;
                    if(sh->fSample->fMorphValue[par]==nominalValue){ // FIXME: eventually add something to flag a sample as nominal for morphing
                        shNominal = sh;
                        break;
                    }
                }
                // loop on all other samples
                for(auto& sh : reg->fSampleHists){
                    if(!sh->fSample->fIsMorph[par]) continue;
                    if(sh != shNominal){
                        for(const auto& syh : shNominal->fSyst){
                            std::shared_ptr<Systematic> syst = syh->fSystematic;
                            if(syst->fIsNormOnly){
                                std::shared_ptr<SystematicHist> syhNew = sh->AddOverallSyst(syst->fName,syst->fStoredName,syst->fOverallUp,syst->fOverallDown);
                                syhNew->fSystematic = syst;
                            }
                            else{
                                TH1* hUpNew   = static_cast<TH1*>(syh->fHistUp->Clone());
                                TH1* hDownNew = static_cast<TH1*>(syh->fHistDown->Clone());
                                hUpNew->Divide(shNominal->fHist.get());
                                hUpNew->Multiply(sh->fHist.get());
                                hDownNew->Divide(shNominal->fHist.get());
                                hDownNew->Multiply(sh->fHist.get());
                                std::shared_ptr<SystematicHist> syhNew = sh->AddHistoSyst(syst->fName,syst->fStoredName,hUpNew,hDownNew);
                                syhNew->fSystematic = syst;
                                sh->fSample->fUseSystematics = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Propagate all systematics from another sample
    for(const auto& reg : fRegions){
        for(const auto& smp : fSamples){
            if(smp->fSystFromSample != ""){
                // eventually skip sample / region combination
                if( Common::FindInStringVector(smp->fRegions,reg->fName)<0 ) continue;
                std::shared_ptr<SampleHist> sh = reg->GetSampleHist(smp->fName);
                if(sh==nullptr) continue;
                sh->fSample->fUseSystematics = true;
                //
                std::shared_ptr<SampleHist> shReference = reg->GetSampleHist(smp->fSystFromSample);
                for(const auto& syh : shReference->fSyst){
                    std::shared_ptr<Systematic> syst = syh->fSystematic;
                    if(syst->fIsNormOnly){
                        std::shared_ptr<SystematicHist> syhNew = sh->AddOverallSyst(syst->fName,syst->fStoredName,syst->fOverallUp,syst->fOverallDown);
                        syhNew->fSystematic = syst;
                    }
                    else{
                        TH1* hUpNew   = static_cast<TH1*>(syh->fHistUp->Clone());
                        TH1* hDownNew = static_cast<TH1*>(syh->fHistDown->Clone());
                        hUpNew->Divide(shReference->fHist.get());
                        hUpNew->Multiply(sh->fHist.get());
                        hDownNew->Divide(shReference->fHist.get());
                        hDownNew->Multiply(sh->fHist.get());
                        std::shared_ptr<SystematicHist> syhNew = sh->AddHistoSyst(syst->fName,syst->fStoredName,hUpNew,hDownNew);
                        syhNew->fSystematic = syst;
                    }
                }
            }
        }
    }

    //
    // set the hasData flag
    bool hasData = false;
    for(const auto& smp : fSamples){
        if(smp->fType==Sample::SampleType::DATA){
            hasData = true;
            break;
        }
    }

    //
    // Poissonize data
    if(hasData && TRExFitter::OPTION["PoissonizeData"]!=0){
        for(const auto& reg : fRegions){
            if(reg->fData!=nullptr){
                if(reg->fData->fHist!=nullptr){
                    TH1 *hdata = reg->fData->fHist.get();
                    if(TRExFitter::OPTION["PoissonizeData"]>0) gRandom->SetSeed(TRExFitter::OPTION["PoissonizeData"]);
                    for(int i_bin=1;i_bin<=hdata->GetNbinsX();i_bin++){
                        hdata->SetBinContent(i_bin,gRandom->Poisson( hdata->GetBinContent(i_bin) ));
                        hdata->SetBinError(i_bin,std::sqrt(hdata->GetBinContent(i_bin)));
                    }
                }
            }
        }
    }

    // Manually change the shape of systematics
    RunForceShape();
}

//__________________________________________________________________________________
//
void TRExFit::CloseInputFiles(){
    //
    // Close all input files
    for(auto& it : TRExFitter::TFILEMAP){
        TDirectory *dir = gDirectory;
        TFile* f = it.second.get();
        if(f) {
            dir->cd();
            f->Close();
        }
    }
    TRExFitter::TFILEMAP.clear();
}

//__________________________________________________________________________________
//
void TRExFit::DrawAndSaveAll(std::string opt, const std::string& wsPath){
    bool isPostFit = opt.find("post")!=std::string::npos;
    //
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    opt += " poissonize";
    if(isPostFit) {
        if(fFitResultsRootFile!=""){
            ReadFitResults(fFitResultsRootFile);
        } else {
            if (fBootstrap!="" && fBootstrapIdx>=0){
                ReadFitResults(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".root");
            } else {
                ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
            }
        }
    }

    std::string frFile("");
    if (isPostFit) {
        if (fFitResultsRootFile != "") {
            frFile = fFitResultsRootFile;
        } else {
            frFile = fName + "/Fits/"+fInputName+fSuffix+".root";
        }
    }


    xRooNode rooNode = Common::ReadWSandFitResults(wsPath, frFile);

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf\n";
        exit(EXIT_FAILURE);
    }

    std::shared_ptr<xRooNode> pdfPrefit(nullptr);
    if (TRExFitter::PREFITONPOSTFIT && isPostFit) {

        xRooNode rooNodePrefit = Common::ReadWSandFitResults(wsPath, "");
        pdfPrefit = rooNodePrefit["simPdf"];
        if (!pdfPrefit) {
            LOG(ERROR) << "Cannot read ModelConfig: simPdf for prefit\n";
            exit(EXIT_FAILURE);
        }
    }

    for(auto& ireg : fRegions) {
        std::shared_ptr<TRExPlot> p(nullptr);
        ireg->fUseStatErr = fUseStatErr;
        ireg->fPlotLabel = fPlotLabel;
        ireg->fPreFitLabel = fPreFitLabel;
        ireg->fPostFitLabel = fPostFitLabel;
        ireg->fToysForErrorBand = fToysForErrorBand;

        if(fCustomAsimov!=""){
            std::string name = "customAsimov_"+fCustomAsimov;
            std::shared_ptr<SampleHist> cash = ireg->GetSampleHist(name);
            if(cash==nullptr){
                LOG(WARNING) << "No Custom Asimov " << fCustomAsimov << " available. Taking regular Asimov.\n";
            }
            else{
                std::string s = cash->fHist->GetName();
                LOG(DEBUG) << "  Adding Custom-Asimov Data: " << s << "\n";
                ireg->fData = cash;
            }
        }
        for (const auto& iregion : fRegions) {
            iregion->fFolder = fName;
            iregion->fHEPDataFormat = fHEPDataFormat;
            iregion->fUheppFormat = fUheppFormat;
        }
        //

        std::string outputDir("");
        if (ireg->fPlotSubdir != ""){
            outputDir = fName + "/Plots/" + ireg->fPlotSubdir + "/";
            gSystem->mkdir(outputDir.c_str());
        } else {
            outputDir = fName + "/Plots/";
        }
        ireg->fNormFactors = fNormFactors;
        if(isPostFit){
            gSystem->mkdir( (fName + "/Histograms/").c_str() );
            if(ireg->fRegionDataType==Region::ASIMOVDATA) p = ireg->DrawPrePostFit(pdf, pdfPrefit, fFitResults.get(),fPrePostFitCanvasSize,true, fFitType == FitType::BONLY, opt+" blind");
            else                                          p = ireg->DrawPrePostFit(pdf, pdfPrefit, fFitResults.get(),fPrePostFitCanvasSize,true, fFitType == FitType::BONLY, opt);
            p->SaveAs(outputDir+ireg->fName+"_postFit"+fSuffix);
            p->SaveAsBkgOnly(outputDir+ireg->fName+"_bkgOnly_postFit"+fSuffix);
        }
        else{
            if(ireg->fRegionDataType==Region::ASIMOVDATA) p = ireg->DrawPrePostFit(pdf, nullptr, nullptr, fPrePostFitCanvasSize, false, fFitType == FitType::BONLY, opt+" blind");
            else                                          p = ireg->DrawPrePostFit(pdf, nullptr, nullptr, fPrePostFitCanvasSize, false, fFitType == FitType::BONLY, opt);
            // this line to fix the y-axis maximum getting doubled in some cases (FIXME)
            if((ireg->fYmin==0) && (ireg->fYmax==0) && (ireg->fYmaxScale==0)){
                if(!ireg->fLogScale) p->h_dummy->GetYaxis()->SetRangeUser(p->h_dummy->GetYaxis()->GetXmin(),p->h_dummy->GetMaximum());
                else                 p->h_dummy->GetYaxis()->SetRangeUser(1                                ,p->h_dummy->GetMaximum());
            }
            p->SaveAs(outputDir+ireg->fName+fSuffix);
        }
    }
}

//__________________________________________________________________________________
//
std::shared_ptr<TRExPlot> TRExFit::DrawSummary(std::string opt, std::shared_ptr<TRExPlot> prefit_plot, const std::string& wsPath) {
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Building Summary Plot...\n";
    gSystem->mkdir(fName.c_str(),true);
    const bool isPostFit = opt.find("post")!=std::string::npos;
    const bool checkVR = opt.find("valid")!=std::string::npos;
    opt += " poissonize";
    // build one bin per region
    std::unique_ptr<TH1D> data = nullptr;

    // Vectors of histos in the form { Histogram(sample1, allRegions), Histogram(sample2, allRegions), ... }
    std::vector<std::shared_ptr<TH1D> > sigPostfit;
    std::vector<std::shared_ptr<TH1D> > sigPostfitBonlyNorm;
    std::vector<std::shared_ptr<TH1D> > bkgPostfit;
    std::vector<std::shared_ptr<TH1D> > sigPrefit;
    std::vector<std::shared_ptr<TH1D> > bkgPrefit;

    std::string frFile("");
    if (isPostFit) {
        if (fFitResultsRootFile != "") {
            frFile = fFitResultsRootFile;
        } else {
            frFile = fName + "/Fits/"+fInputName+fSuffix+".root";
        }
    }

    xRooNode rooNode = Common::ReadWSandFitResults(wsPath, frFile);

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> regions;
    std::vector<int> divisionVec;
    for (const auto& ireg : fRegions) {
        if (checkVR) {
            if (!fSummaryPlotValidationRegions.empty()) {
                if (Common::FindInStringVector(fSummaryPlotValidationRegions, ireg->fName) < 0) continue;
            }
        } else {
            if (ireg->fRegionType == Region::RegionType::VALIDATION) continue;
            if (!fSummaryPlotRegions.empty()) {
                if (Common::FindInStringVector(fSummaryPlotRegions, ireg->fName) < 0) continue;
            }
        }

        regions.emplace_back(ireg->fName);
    }

    if(regions.empty()) return nullptr;
    const int nbins = regions.size();

    for(const auto& isample : fSamples) {
        if(isample->fType == Sample::SampleType::GHOST) continue;
        if(isample->fType == Sample::SampleType::EFT) continue;
        const std::string name = isample->fName;
        std::string title = isample->fTitle;
        if(isample->fGroup != "") title = isample->fGroup;
        // look for the first SampleHist defined for this sample
        const SampleHist* sh(nullptr);
        for (const auto& ireg : fRegions) {
            sh = ireg->GetSampleHist(name).get();
            if (sh) break;
        }
        // skip sample if no SampleHist found
        if(!sh) continue;
        if(!sh->fHist) continue;
        //
        const int lineColor = sh->fHist->GetLineColor();
        const int fillColor = sh->fHist->GetFillColor();
        const int lineWidth = sh->fHist->GetLineWidth();

        auto hPre = std::make_shared<TH1D>("", title.c_str(), nbins, 0, nbins);
        auto hPost = std::make_shared<TH1D>("", title.c_str(), nbins, 0, nbins);
        auto hPostBonlyNorm = std::make_shared<TH1D>("", title.c_str(), nbins, 0, nbins);

        if(isample->fType != Sample::SampleType::DATA) {
            int ibin(1);
            for (const auto& ireg : regions) {
                // skip samples not in this region
                if ((isample->fRegions.at(0) != "all") && (Common::FindInStringVector(isample->fRegions, ireg) < 0)) {
                    ++ibin;
                    continue;
                }
                auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&ireg](const auto& region){return ireg == region->fName;});
                if (itr == fRegions.end()) {
                    LOG(ERROR) << "Cannot find region: " << ireg << "\n";
                    exit(EXIT_FAILURE);
                }
                double content       = (*itr)->GetIntegral(name, isPostFit);
                double contentPrefit = (!isPostFit || TRExFitter::PREFITONPOSTFIT) ? (*itr)->GetIntegral(name, false) : content;

                hPre->SetBinContent(ibin, contentPrefit);
                hPre->SetBinError(ibin, 0);
                hPost->SetBinContent(ibin, content);
                hPost->SetBinError(ibin, 0);
                if (fFitType == FitType::BONLY && isample->fType == Sample::SampleType::SIGNAL && TRExFitter::SHOWNORMSIG && isPostFit) {
                    const double c = (*itr)->GetIntegralBonlyNorm(name);
                    hPostBonlyNorm->SetBinContent(ibin, c);
                    hPostBonlyNorm->SetBinError(ibin, 0);
                }

                ++ibin;
            }
        }
        hPre->SetLineColor(lineColor);
        hPre->SetFillColor(fillColor);
        hPre->SetLineWidth(lineWidth);
        hPost->SetLineColor(lineColor);
        hPost->SetFillColor(fillColor);
        hPost->SetLineWidth(lineWidth);
        hPostBonlyNorm->SetLineColor(lineColor);
        hPostBonlyNorm->SetFillColor(fillColor);
        hPostBonlyNorm->SetLineWidth(lineWidth);
        //
        if(isample->fType == Sample::SampleType::SIGNAL){
            sigPrefit.emplace_back(std::move(hPre));
            sigPostfit.emplace_back(std::move(hPost));
            sigPostfitBonlyNorm.emplace_back(std::move(hPostBonlyNorm));
        } else if(isample->fType == Sample::SampleType::BACKGROUND){
            bkgPrefit.emplace_back(std::move(hPre));
            bkgPostfit.emplace_back(std::move(hPost));
        } else if (isample->fType == Sample::SampleType::DATA) {
            data = std::make_unique<TH1D>("",title.c_str(), nbins, 0, nbins);
            data->SetDirectory(nullptr);
            int ibin = 1;
            for (const auto& ireg : regions) {
                auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&ireg](const auto& region){return ireg == region->fName;});
                if (itr == fRegions.end()) {
                    LOG(ERROR) << "Cannot find region: " << ireg << "\n";
                    exit(EXIT_FAILURE);
                }
                if ((*itr)->fRegionDataType == Region::ASIMOVDATA) {
                    data->SetBinContent(ibin, 0);
                } else if ((*itr)->fData && (*itr)->fData->fHist) {
                    data->SetBinContent(ibin, (*itr)->fData->fHist->Integral());
                } else {
                    LOG(WARNING) << "Data not provided for Region: " << (*itr)->fName << " setting the content to 0\n";
                    data->SetBinContent(ibin, 0);
                }

                ++ibin;
            }
        }
    }
    //
    std::shared_ptr<TRExPlot> p;
    //
    if (fSummaryCanvasSize.size() == 0){
        p = std::make_shared<TRExPlot>(fInputName+"_summary",900,700,TRExFitter::NORATIO);
    } else {
        p = std::make_shared<TRExPlot>(fInputName+"_summary",fSummaryCanvasSize.at(0),fSummaryCanvasSize.at(1),TRExFitter::NORATIO);
    }
    if(fYmin!=0) p->fYmin = fYmin;
    else         p->fYmin = 1;
    if(fYmax!=0) p->fYmax = fYmax;
    else         p->SetYmaxScale(2);
    p->SetXaxis("",false);
    p->AddLabel(fLabel);
    if (!fNoPrePostFitLabel) {
        if(isPostFit) p->AddLabel(fPostFitLabel);
        else          p->AddLabel(fPreFitLabel);
    }
    //
    if(isPostFit) p->fRatioYmax = fRatioYmaxPostFit;
    else          p->fRatioYmax = fRatioYmax;
    if(isPostFit) p->fRatioYmin = fRatioYminPostFit;
    else          p->fRatioYmin = fRatioYmin;
    //
    // propagate settings from Job to plot
    p->fRatioYtitle = fRatioYtitle;
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType==TRExPlot::RATIOTYPE::DATAOVERMC){
        p->fRatioType = TRExPlot::RATIOTYPE::DATAOVERB;
    }
    p->fPlotLabel = fPlotLabel;
    p->SetLumi(fLumiLabel);
    p->SetCME(fCmeLabel);
    p->SetLumiScale(fLumiScale);
    p->fLabelX = fLabelXSummary;
    p->fLabelY = fLabelYSummary;
    p->fLegendX1 = fLegendX1Summary;
    p->fLegendX2 = fLegendX2Summary;
    p->fLegendY = fLegendYSummary;
    p->fLegendNColumns = fLegendNColumnsSummary;


    std::vector<int> blindedBinsFromRegionTotals;
    if(fBlindingThreshold >= 0) {

        std::unique_ptr<TH1> signal(nullptr);
        std::unique_ptr<TH1> bkg(nullptr);

        if (isPostFit && !fKeepPrefitBlindedBins) {
            if (sigPostfit.size() > 0) {
                signal = Common::CombineHistosFromHistosVec(sigPostfit);
            }
            bkg = Common::CombineHistosFromHistosVec(bkgPostfit);
        }
        else{
            if (sigPrefit.size() > 0) {
                signal = Common::CombineHistosFromHistosVec(sigPrefit);
            }
            bkg = Common::CombineHistosFromHistosVec(bkgPrefit);
        }

        blindedBinsFromRegionTotals = Common::ComputeBlindedBins(signal.get(),
                                                                 bkg.get(),
                                                                 fBlindingType,
                                                                 fBlindingThreshold);
    }
    //
    if (data) p->SetData(data.get(), data->GetTitle());
    for (std::size_t i = 0; i < sigPostfit.size(); i++) {
        if (isPostFit) {
            if(TRExFitter::SHOWSTACKSIG_SUMMARY) p->AddSignal(sigPostfit.at(i).get(),sigPostfit.at(i)->GetTitle());
            if(TRExFitter::SHOWNORMSIG_SUMMARY) {
                if (fFitType == FitType::BONLY && TRExFitter::SHOWNORMSIG) {
                    p->AddNormSignal(sigPostfitBonlyNorm.at(i).get(),sigPostfit.at(i)->GetTitle());
                } else {
                    p->AddNormSignal(sigPostfit.at(i).get(),sigPostfit.at(i)->GetTitle());
                }
            }
            if(TRExFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(sigPostfit.at(i).get(),sigPostfit.at(i)->GetTitle());
        } else {
            if(TRExFitter::SHOWSTACKSIG_SUMMARY)   p->AddSignal(    sigPrefit.at(i).get(),sigPrefit.at(i)->GetTitle());
            if(TRExFitter::SHOWNORMSIG_SUMMARY)    p->AddNormSignal(sigPrefit.at(i).get(),sigPrefit.at(i)->GetTitle());
            if(TRExFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(sigPrefit.at(i).get(),sigPrefit.at(i)->GetTitle());
        }
    }
    for (std::size_t i = 0; i < bkgPostfit.size(); i++) {
        if (isPostFit) {
            p->AddBackground(bkgPostfit.at(i).get(),bkgPostfit.at(i)->GetTitle());
        } else {
            p->AddBackground(bkgPrefit.at(i).get(),bkgPrefit.at(i)->GetTitle());
        }
    }

    if( TRExFitter::PREFITONPOSTFIT && isPostFit) {
        p->h_tot_bkg_prefit = static_cast<TH1*>(prefit_plot->GetTotBkg()->Clone("h_tot_bkg_prefit"));
    }

    //
    // Build tot
    //
    std::unique_ptr<TH1D> tot = std::make_unique<TH1D>("h_Tot_summary","h_Tot_summary", nbins,0,nbins);
    std::unique_ptr<TGraphAsymmErrors> error = std::make_unique<TGraphAsymmErrors>(nbins);
    int ibin = 1;

    std::vector<int> blindedRegionBinNum;
    for (const auto& ireg : regions) {
        auto channel = pdf->find(ireg);

        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&ireg](const auto& region){return ireg == region->fName;});
        std::vector<int> blindedBins = Common::CombineVectors((*itr)->fDropBins, (*itr)->fBlindedBins,true);
        if (blindedBins.size() > 0) {
            blindedRegionBinNum.push_back(ibin);
        }

        tot->SetBinContent(ibin, channel->GetContent());
        tot->SetBinError(ibin, 0.);
        error->SetPoint(ibin-1, ibin-0.5, channel->GetContent());
        error->SetPointError(ibin-1, 0.5, 0.5, channel->GetError(), channel->GetError());
        ++ibin;
    }

    p->SetTot(tot.get());
    p->SetTotAsym(error.get());
    // Blind if one bin in region is blinded/dropped, entire region is blinded based on integral
    blindedRegionBinNum = Common::CombineVectors(blindedRegionBinNum, blindedBinsFromRegionTotals, true);
    p->SetBinBlinding(blindedRegionBinNum);
    p->BlindData();

    p->ResizeBinLabel(nbins+1);
    ibin = 1;
    for (const auto& ireg : regions) {
        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&ireg](const auto& region){return ireg == region->fName;});
        p->SetBinLabel(ibin, (*itr)->fShortLabel.c_str());
        ++ibin;
    }

    p->Draw(opt);

    if(divisionVec.size()>0){
        p->pad0->cd();
        TLine line(0,1,0,1);
        line.SetNDC(0);
        line.SetLineStyle(7);
        line.SetLineColor(kBlack);
        line.SetLineWidth(2);
        line.DrawLine(divisionVec.at(0),(static_cast<TH1D*>(p->pad0->GetPrimitive("h_dummy")))->GetMinimum(),
                      divisionVec.at(0),std::pow((static_cast<TH1D*>(p->pad0->GetPrimitive("h_dummy")))->GetMaximum(),0.73) );
        p->pad1->cd();
        line.DrawLine(divisionVec.at(0),(static_cast<TH1D*>(p->pad1->GetPrimitive("h_dummy2")))->GetMinimum(),
                      divisionVec.at(0),(static_cast<TH1D*>(p->pad1->GetPrimitive("h_dummy2")))->GetMaximum());
    }
    //
    p->pad0->cd();
    if(!checkVR){
        if(fSummaryPlotLabels.size()>0){
            TLatex tex;
            tex.SetNDC(0);
            tex.SetTextAlign(20);
            //
            for(unsigned int ii=0;ii<=divisionVec.size();ii++){
                if(fSummaryPlotLabels.size()<ii+1) break;
                double xmax = nbins;
                double xmin = 0.;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                double xpos = xmin + 0.5*(xmax - xmin);
                double ypos = std::pow((static_cast<TH1D*>(p->pad0->GetPrimitive("h_dummy")))->GetMaximum(), 0.61);
                tex.DrawLatex(xpos,ypos,fSummaryPlotLabels[ii].c_str());
            }
        }
    }
    else{
        if(fSummaryPlotValidationLabels.size()>0){
            TLatex tex;
            tex.SetNDC(0);
            tex.SetTextAlign(20);
            //
            for(unsigned int ii=0;ii<=divisionVec.size();ii++){
                if(fSummaryPlotValidationLabels.size()<ii+1) break;
                double xmax = nbins;
                double xmin = 0.;
                if(divisionVec.size()>ii) xmax = divisionVec[ii];
                if(ii>0) xmin = divisionVec[ii-1];
                double xpos = xmin + 0.5*(xmax - xmin);
                double ypos = std::pow((static_cast<TH1D*>(p->pad0->GetPrimitive("h_dummy")))->GetMaximum(), 0.61 );
                tex.DrawLatex(xpos,ypos,fSummaryPlotValidationLabels[ii].c_str());
            }
        }
    }
    //
    for(int i_bin = 1; i_bin <= nbins; i_bin++){
        LOG(DEBUG) << i_bin << ":\t" << tot->GetBinContent(i_bin) << "\t+" << error->GetErrorYhigh(i_bin-1) << "\t-" << error->GetErrorYlow(i_bin-1) << "\n";
    }
    //
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    if(fSummaryPrefix!=""){
        if(isPostFit)  p->SaveAs(fName+"/Plots/"+fSummaryPrefix+"_Summary_postFit"+(checkVR?"_VR":"")+fSuffix);
        else           p->SaveAs(fName+"/Plots/"+fSummaryPrefix+"_Summary"        +(checkVR?"_VR":"")+fSuffix);
    }
    else{
        if(isPostFit)  p->SaveAs(fName+"/Plots/Summary_postFit"+(checkVR?"_VR":"")+fSuffix);
        else           p->SaveAs(fName+"/Plots/Summary"        +(checkVR?"_VR":"")+fSuffix);
    }
    return p;
}

//__________________________________________________________________________________
//
void TRExFit::DrawMergedPlot(std::string opt, const std::string& group) const{
    std::vector<Region*> regions;
    if(group=="") {
        for (auto& ireg : fRegions) {
            regions.push_back(ireg.get());
        }
    }
    else{
        for(auto& region : fRegions){
            if(region->fGroup == group) regions.push_back(region.get());
        }
    }
    bool isPostFit = false;
    if(opt.find("post")!=std::string::npos) isPostFit = true;
    opt += " poissonize";
    // start with total prediction, which should be always there
    // build a vector of histograms
    int i_ch = 0;
    std::vector<std::unique_ptr<TH1> > hTotVec;
    std::vector<double> edges;
    std::vector<std::unique_ptr<TGaxis> > xaxis;
    std::vector<std::unique_ptr<TGaxis> > yaxis;
    //
    double ymax0 = -1.; // ymax0 is the max y of the first region
    double ymax  = -1.;
    double ymin  = 0.;
    if(TRExFitter::OPTION["MergeLogY"]>0) ymin = 1.;
    if(TRExFitter::OPTION["MergeLogYmin"]>0) ymin = TRExFitter::OPTION["MergeLogYmin"];
    // retrieve from config file whether last bin width needs to be fixed for merge plot
    bool fixLastBinWidth = TRExFitter::OPTION["MergeFixLastBinWidth"]>0;
    for(auto region : regions){
        std::unique_ptr<TH1> h_tmp  = nullptr;
        if(isPostFit) h_tmp = std::unique_ptr<TH1>(static_cast<TH1*>(region->fTot_postFit->Clone()));
        else          h_tmp = std::unique_ptr<TH1>(static_cast<TH1*>(region->fTot->Clone()));
        std::unique_ptr<TH1> h_data = nullptr;
        if(region->fData!=nullptr) h_data = std::unique_ptr<TH1>(static_cast<TH1*>(region->fData->fHist->Clone()));
        //
        for(int i_bin=1;i_bin<=h_tmp->GetNbinsX();i_bin++){
            if(isPostFit) h_tmp->SetBinError( i_bin,region->fErr_postFit->GetErrorY(i_bin-1) );
            else          h_tmp->SetBinError( i_bin,region->fErr->GetErrorY(i_bin-1) );
        }
        // find max y (don't rely on GetMaximum, since it could have been modified
        double ymaxTmp = 0;
        for(int i_bin=1;i_bin<=h_tmp->GetNbinsX();i_bin++){
            float k = 1.; // factor != 1 in case of y axis normalized by bin width
            if(regions[0]->fBinWidth>0){
                k = regions[0]->fBinWidth/h_tmp->GetXaxis()->GetBinWidth(i_bin);
            }
            if(k*h_tmp->GetBinContent(i_bin)>ymaxTmp) ymaxTmp = k*h_tmp->GetBinContent(i_bin);
            if(h_data!=nullptr)
                if(k*(h_data->GetBinContent(i_bin)+h_data->GetBinError(i_bin))>ymaxTmp) ymaxTmp = k*(h_data->GetBinContent(i_bin)+h_data->GetBinError(i_bin));
        }
        if(i_ch>0){
            if(TRExFitter::OPTION["MergeLogY"]>0) ymaxTmp = pow(ymaxTmp,1.1);
            else ymaxTmp *= 1.1;
        }
        h_tmp->SetMaximum(ymaxTmp);
        // set max for first hist
        if(ymax0<0){
            ymax0 = ymaxTmp;
            if(TRExFitter::OPTION["MergeYfrac"]==0) TRExFitter::OPTION["MergeYfrac"] = 0.7;
            ymax  = TRExFitter::OPTION["MergeYfrac"]*ymax0;
        }
        double binUpEdge = h_tmp->GetXaxis()->GetBinUpEdge(h_tmp->GetNbinsX());
        if(fixLastBinWidth){
            binUpEdge = h_tmp->GetXaxis()->GetBinLowEdge(h_tmp->GetNbinsX()) + h_tmp->GetXaxis()->GetBinUpEdge(h_tmp->GetNbinsX()-1) - h_tmp->GetXaxis()->GetBinLowEdge(h_tmp->GetNbinsX()-1);
        }
        if(i_ch>0) edges.push_back(binUpEdge - h_tmp->GetXaxis()->GetBinLowEdge(1) + edges.at(i_ch-1));
        else       edges.push_back(binUpEdge);
        double binUpTmp = h_tmp->GetXaxis()->GetBinLowEdge(h_tmp->GetNbinsX()) + h_tmp->GetXaxis()->GetBinWidth(h_tmp->GetNbinsX()-1);
        if(i_ch>0) xaxis.emplace_back(std::make_unique<TGaxis>(edges.at(i_ch-1),                   0,edges.at(i_ch),0, h_tmp->GetXaxis()->GetBinLowEdge(1),binUpTmp ,510,"+S"));
        else       xaxis.emplace_back(std::make_unique<TGaxis>(h_tmp->GetXaxis()->GetBinLowEdge(1),0,edges.at(i_ch),0, h_tmp->GetXaxis()->GetBinLowEdge(1),binUpTmp ,510,"+S"));
        if(TRExFitter::OPTION["MergeUniqueY"]==0){
            // get yaxes
            if(i_ch>0) yaxis.emplace_back(std::make_unique<TGaxis>(edges.at(i_ch-1),                   ymin,edges.at(i_ch-1),                   ymax, ymin,ymaxTmp, 510,"-S"));
            else       yaxis.emplace_back(std::make_unique<TGaxis>(h_tmp->GetXaxis()->GetBinLowEdge(1),ymin,h_tmp->GetXaxis()->GetBinLowEdge(1),ymax, ymin,ymaxTmp, 510,"-S"));
        }
        i_ch ++;
        hTotVec.emplace_back(std::move(h_tmp));
    }
    // then proceed with data, singnal and bkg
    std::vector<TH1*> hDataVec;
    std::vector<std::vector<TH1*>> hSignalVec;
    std::vector<std::vector<TH1*>> hBackgroundVec;
    for(const auto& sample : fSamples){
        if(sample->fType==Sample::SampleType::GHOST) continue;
        if(sample->fType==Sample::SampleType::EFT) continue;
        std::vector<TH1*> tmpVec;
        i_ch = 0;
        for(const auto* const region : regions){
            TH1* h_tmp = nullptr;
            for(const auto& sampleHist : region->fSampleHists){
                if(sampleHist->fSample->fName == sample->fName){
                    if(isPostFit && sample->fType!=Sample::SampleType::DATA){
                        if(sampleHist->fHist_postFit!=nullptr){
                            h_tmp = static_cast<TH1*>(sampleHist->fHist_postFit->Clone());
                        }
                    }
                    else{
                        if(sampleHist->fHist!=nullptr){
                            h_tmp = static_cast<TH1*>(sampleHist->fHist->Clone());
                            if(!sampleHist->fSample->fUseMCStat && !sampleHist->fSample->fSeparateGammas){
                                for(int i_bin=0;i_bin<h_tmp->GetNbinsX()+2;i_bin++) h_tmp->SetBinError(i_bin,0.);
                            }
                        }
                    }
                    break;
                }
            }
            // if the sample was not in the region...
            if(h_tmp==nullptr){
                h_tmp = static_cast<TH1*>(hTotVec.at(i_ch)->Clone());
                h_tmp->Scale(0);
            }
            if(sample->fGroup!="") h_tmp->SetTitle(sample->fGroup.c_str());
            else                   h_tmp->SetTitle(sample->fTitle.c_str());
            tmpVec.push_back(h_tmp);
            //
            i_ch ++;
        }
        if(sample->fType==Sample::SampleType::DATA)            hDataVec = tmpVec;
        else if(sample->fType==Sample::SampleType::SIGNAL)     hSignalVec.push_back(tmpVec);
        else if(sample->fType==Sample::SampleType::BACKGROUND) hBackgroundVec.push_back(tmpVec);
    }
    //
    // scale them (but the first region)
    if(TRExFitter::OPTION["MergeUniqueY"]==0){
        for(unsigned int i_channel=1;i_channel<regions.size();i_channel++){
            double scale = ymax/hTotVec.at(i_channel)->GetMaximum();
            hTotVec.at(i_channel)->Scale( scale );
            if(hDataVec.size()>0){
                for(int i_bin=1;i_bin<=hDataVec.at(i_channel)->GetNbinsX();i_bin++){
                    hDataVec.at(i_channel)->SetBinError(i_bin,std::sqrt(hDataVec.at(i_channel)->GetBinContent(i_bin)));
                }
                hDataVec.at(i_channel)->Scale( scale );
            }
            for(const auto& hVec : hSignalVec)     hVec.at(i_channel)->Scale( scale );
            for(const auto& hVec : hBackgroundVec) hVec.at(i_channel)->Scale( scale );
        }
    }
    //
    // merge and plot them
    std::shared_ptr<TRExPlot> p;
    int cWidthMerge = 1200;
    int cHeightMerge = 600;
    if(fMergeCanvasSize.size()>1){
        cWidthMerge = fMergeCanvasSize.at(0);
        cHeightMerge = fMergeCanvasSize.at(1);
    }
    p = std::make_shared<TRExPlot>(fInputName+"_merge",cWidthMerge,cHeightMerge,TRExFitter::NORATIO);
    //

    if(hDataVec.size()>0){
        std::unique_ptr<TH1> merged(Common::MergeHistograms(hDataVec,fixLastBinWidth));
        p->SetData(merged.get(),"");
    }
    for(unsigned int i_sig=0;i_sig<hSignalVec.size();i_sig++){
        std::unique_ptr<TH1> tmp(Common::MergeHistograms(hSignalVec.at(i_sig),fixLastBinWidth));
        if(TRExFitter::SHOWSTACKSIG_SUMMARY)   p->AddSignal(    tmp.get(),"");
        if(TRExFitter::SHOWNORMSIG_SUMMARY)    p->AddNormSignal(tmp.get(),"");
        if(TRExFitter::SHOWOVERLAYSIG_SUMMARY) p->AddOverSignal(tmp.get(),"");

    }
    for(unsigned int i_bkg=0;i_bkg<hBackgroundVec.size();i_bkg++) {
        std::unique_ptr<TH1> tmp(Common::MergeHistograms(hBackgroundVec.at(i_bkg),fixLastBinWidth));
        p->AddBackground(tmp.get(),"");
    }
    std::unique_ptr<TH1> total(Common::MergeHistograms(hTotVec,fixLastBinWidth));
    p->SetTot(total.get());
    //
    p->SetCME(fCmeLabel);
    p->SetLumi(fLumiLabel);
    p->fPlotLabel = fPlotLabel;
    p->fLabelX = fLabelXMerge;
    p->fLabelY = fLabelYMerge;
    p->fLegendX1 = fLegendX1Merge;
    p->fLegendX2 = fLegendX2Merge;
    p->fLegendY = fLegendYMerge;
    p->fLegendNColumns = fLegendNColumnsMerge;
    p->SetXaxis(regions.at(0)->fVariableTitle);
    if(regions[0]->fBinWidth>0) p->SetBinWidth(regions[0]->fBinWidth); // Take bin width from first region
    p->fRatioType = fRatioType;
    if(!(TRExFitter::SHOWSTACKSIG && TRExFitter::ADDSTACKSIG) && fRatioType==TRExPlot::RATIOTYPE::DATAOVERMC){
        p->fRatioType = TRExPlot::RATIOTYPE::DATAOVERB;
    }
    if(fBlindingThreshold >= 0) {
        std::unique_ptr<TH1> signal(nullptr);
        if (hSignalVec.size() > 0) {
            signal.reset(static_cast<TH1*>(Common::MergeHistograms(hSignalVec.at(0),fixLastBinWidth)->Clone()));
        }
        std::unique_ptr<TH1> bkg(nullptr);
        for (std::size_t i = 0; i < hBackgroundVec.size(); ++i) {
            std::unique_ptr<TH1> tmp(Common::MergeHistograms(hBackgroundVec.at(i),fixLastBinWidth));
            if (!tmp) continue;
            if (!bkg) bkg.reset(static_cast<TH1*>(tmp->Clone()));
            else      bkg->Add(tmp.get());
        }
        const std::vector<int>& blindedBins = Common::ComputeBlindedBins(signal.get(),
                                                                         bkg.get(),
                                                                         fBlindingType,
                                                                         fBlindingThreshold);
        p->SetBinBlinding(blindedBins);
    }

    if(TRExFitter::OPTION["MergeYmaxScale"]==0) TRExFitter::OPTION["MergeYmaxScale"] = 1.3;
    if(TRExFitter::OPTION["MergeLogY"]>0){
        p->fYmax = pow(ymax0,TRExFitter::OPTION["MergeYmaxScale"])/ymin;
    }
    else{
        p->fYmax = TRExFitter::OPTION["MergeYmaxScale"]*ymax0;
    }
    if(fRatioYmax>0) p->fRatioYmax = fRatioYmax;
    if(fRatioYmin>0) p->fRatioYmin = fRatioYmin;
    if(isPostFit && fRatioYmaxPostFit>0) p->fRatioYmax = fRatioYmaxPostFit;
    if(isPostFit && fRatioYminPostFit>0) p->fRatioYmin = fRatioYminPostFit;
    gStyle->SetTickLength(0.0,"x"); // this will switch off the automatic tick marks on the x axis - new custom axes will be added on top of the x axis
    p->Draw(opt);
    //
    // manipulate canvas / pad
    //
    // dahsed line in ratio
    p->pad1->cd();
    std::vector<std::unique_ptr<TLine> > l;
    for(auto edge : edges){
        std::unique_ptr<TLine> l_tmp = std::make_unique<TLine>(edge,(static_cast<TH1*>(p->pad1->GetPrimitive("h_dummy2")))->GetMinimum(),
                                                               edge,(static_cast<TH1*>(p->pad1->GetPrimitive("h_dummy2")))->GetMaximum());
        l_tmp->SetLineStyle(kDashed);
        l_tmp->Draw("same");
        l.emplace_back(std::move(l_tmp));
    }
    //
    // (dahsed) line in main pad
    p->pad0->cd();
    std::vector<std::unique_ptr<TLine> > l1;
    for(auto edge : edges){
        std::unique_ptr<TLine> l_tmp = std::make_unique<TLine>(edge,ymin,edge,TRExFitter::OPTION["MergeLogY"]>0 ? pow(ymax,1.25)/ymin : 1.25*ymax);
        if(TRExFitter::OPTION["MergeUniqueY"]!=0) l_tmp->SetLineStyle(kDashed);
        l_tmp->Draw("same");
        l1.emplace_back(std::move(l_tmp));
    }
    //
    TH1* h_dummy = static_cast<TH1*>(p->pad0->GetPrimitive("h_dummy"));
    //
    // y-axis
    p->pad0->cd();
    int i_yaxis = 0;
    for(auto& a : yaxis){
        if(i_yaxis>0){
            if(TRExFitter::OPTION["MergeLogY"]>0) a->SetOption("G");
            a->SetTickLength((ymax0/ymax)*h_dummy->GetYaxis()->GetTickLength());
            a->SetLabelFont(gStyle->GetTextFont());
            a->SetLabelSize(gStyle->GetTextSize());
            a->SetNdivisions(805);
            a->ChangeLabel(1,-1,-1,-1,-1,-1," ");
            a->Draw();
        }
        i_yaxis ++;
    }
    //
    // x-axis
    p->pad0->cd();
    h_dummy->GetXaxis()->SetLabelSize(0);
    double totXaxisLenght = 0.;
    for(auto& a : xaxis){
        totXaxisLenght += a->GetX2()-a->GetX1();
    }
    for(auto& a : xaxis){
        a->SetLabelFont(gStyle->GetTextFont());
        a->SetLabelSize(gStyle->GetTextSize());
        a->SetNdivisions(805);
        a->SetTickLength( 0.02 * totXaxisLenght / (a->GetX2()-a->GetX1()) );
        (static_cast<TGaxis*>((a->DrawClone())))->SetLabelSize(0);
    }
    h_dummy->GetYaxis()->SetTitleOffset(2.0);
    // ratio
    p->pad1->cd();
    p->pad1->SetTickx(0);
    TH1* h_dummy2 = static_cast<TH1*>(p->pad1->GetPrimitive("h_dummy2"));
    h_dummy2->GetXaxis()->SetLabelSize(0);
    h_dummy2->GetYaxis()->SetTitleOffset(2.0);
    unsigned int i_reg = 0;
    for(auto& a : xaxis){
        a->SetLabelFont(gStyle->GetTextFont());
        a->SetLabelSize(gStyle->GetTextSize());
        a->SetNdivisions(805);
        a->SetTickLength( 0.06 * totXaxisLenght / (a->GetX2()-a->GetX1()) );
        TGaxis *ga = static_cast<TGaxis*>(a->DrawClone());
        if(fRatioYmin!=0) { ga->SetY1(fRatioYmin); ga->SetY2(fRatioYmin); }
        if(isPostFit && fRatioYminPostFit!=0) { ga->SetY1(fRatioYminPostFit); ga->SetY2(fRatioYminPostFit); }
        ga->SetTitle(regions.at(i_reg)->fVariableTitle.c_str());
        ga->SetTitleOffset(h_dummy2->GetXaxis()->GetTitleOffset()*0.45*(1200./600.)*(TRExFitter::OPTION["CanvasHeight"]/TRExFitter::OPTION["CanvasWidthMerge"]));
        ga->SetTitleSize(gStyle->GetTextSize());
        ga->SetTitleFont(gStyle->GetTextFont());
        if(i_reg<regions.size()-1){
            ga->ChangeLabel(-1,-1,-1,-1,-1,-1," "); // shut up the last label
        }
        i_reg ++;
    }
    //
    // additional labels
    p->pad0->cd();
    TLatex tex;
    tex.SetTextSize(gStyle->GetTextSize());
    tex.SetTextFont(gStyle->GetTextFont());
    for(unsigned int i_channel=0;i_channel<regions.size();i_channel++){
        tex.SetNDC(0);
        float y = TRExFitter::OPTION["MergeLogY"]>0 ? pow(ymax,1.05)/ymin : 1.15*ymax;
        // MergeRegionLabelAlign: 1=left, 2=center, 3=right - if not set, will use center
        float x = 0.;
        std::string label = "";
        if(TRExFitter::OPTION["MergeRegionLabelAlign"]==0){
            TRExFitter::OPTION["MergeRegionLabelAlign"] = 2;
        }
        if(TRExFitter::OPTION["MergeRegionLabelAlign"]==1){
            tex.SetTextAlign(11);
            x = h_dummy->GetXaxis()->GetXmin();
            if(i_channel>0) x = edges.at(i_channel-1);
            label = "  "+regions.at(i_channel)->fLabel; // add leading space (easiest fix for the label position)
        }
        else if(TRExFitter::OPTION["MergeRegionLabelAlign"]==2){
            tex.SetTextAlign(21);
            float x1 = h_dummy->GetXaxis()->GetXmin();
            float x2 = h_dummy->GetXaxis()->GetXmax();
            // if it's not the first region
            if(i_channel>0) x1 = edges.at(i_channel-1);
            // if it's not the last one
            if(i_channel<regions.size()-1) x2 = edges.at(i_channel);
            x = 0.5*(x1+x2);
            label = regions.at(i_channel)->fLabel;
        }
        else if(TRExFitter::OPTION["MergeRegionLabelAlign"]==3){
            tex.SetTextAlign(31);
            x = h_dummy->GetXaxis()->GetXmax();
            if(i_channel<regions.size()-1) x = edges.at(i_channel);
            label = regions.at(i_channel)->fLabel + " "; // add trailing space (easiest fix for the label position)
        }
        tex.DrawLatex(x,y,label.c_str());
    }
    //
    tex.SetNDC(1);
    double textHeight = 0.05*(672./p->pad0->GetWh());
    double labelY = 1-0.08*(700./p->c->GetWh());
    if(p->fLabelY>=0) labelY = p->fLabelY;
    labelY -= textHeight - 0.015;
    tex.DrawLatex(0.33,labelY,fLabel.c_str());
    if(isPostFit) tex.DrawLatex(0.33,labelY-textHeight,fPostFitLabel.c_str());
    else          tex.DrawLatex(0.33,labelY-textHeight,fPreFitLabel.c_str());
    //
    // final fixes
    if(TRExFitter::OPTION["MergeLogY"]!=0){
        h_dummy->SetMinimum(ymin);
        p->pad0->SetLogy();
    }
    //
    // save image
    p->pad0->RedrawAxis();
    p->pad1->RedrawAxis();
    std::string saveName = fName+"/Plots/";
    if(fSummaryPrefix!="") saveName += fSummaryPrefix+"_";
    saveName += "Merge";
    if(group!="") saveName += "_"+group;
    if(isPostFit) saveName += "_postFit";
    saveName += fSuffix;
    p->SaveAs(saveName);
}

//__________________________________________________________________________________
//
void TRExFit::BuildYieldTable(const std::string& opt, const std::string& group, const std::string& wsPath) const{
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Building Yields Table...\n";
    bool isPostFit = opt.find("post")!=std::string::npos;
    std::ofstream out;
    std::ofstream texout;
    gSystem->mkdir(fName.c_str(),true);
    gSystem->mkdir((fName+"/Tables").c_str());
    std::string suffix = "";
    if(group!="") suffix += "_"+group;
    suffix += fSuffix;
    if(!isPostFit){
        out.open(   (fName+"/Tables/Yields"+suffix+".txt").c_str());
        texout.open((fName+"/Tables/Yields"+suffix+".tex").c_str());
    }
    else{
        out.open(   (fName+"/Tables/Yields_postFit"+suffix+".txt").c_str());
        texout.open((fName+"/Tables/Yields_postFit"+suffix+".tex").c_str());
    }

    std::string frFile("");
    if (isPostFit) {
        if (fFitResultsRootFile != "") {
            frFile = fFitResultsRootFile;
        } else {
            frFile = fName + "/Fits/"+fInputName+fSuffix+".root";
        }
    }

    xRooNode rooNode = Common::ReadWSandFitResults(wsPath, frFile);

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf\n";
        exit(EXIT_FAILURE);
    }

    YamlConverter::TableContainer container;

    std::vector<std::string> regions;
    for (const auto& ireg : fRegions) {
        if(group != "" && ireg->fGroup != group) continue;
        regions.emplace_back(ireg->fName);
    }
    if(regions.empty()) return;

    out << " |       | ";
    for (const auto& ireg : regions) {
        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&ireg](const auto& region){return ireg == region->fName;});
        out << (*itr)->fLabel << " | ";
        container.regionNames.emplace_back((*itr)->fLabel);
    }
    out << std::endl;
    if(fTableOptions.find("STANDALONE")!=std::string::npos){
        texout << "\\documentclass[10pt]{article}" << std::endl;
        texout << "\\usepackage{siunitx}" << std::endl;
        texout << "\\sisetup{separate-uncertainty,table-format=6.3(6)}  % hint: modify table-format to best fit your tables" << std::endl;
        texout << "\\usepackage[margin=0.1in,landscape,papersize={210mm,350mm}]{geometry}" << std::endl;
        texout << "\\begin{document}" << std::endl;
    }
    // if not STANDALONE, add a comment in the tex saying that one needs to include siunitx
    else{
        texout << "% NB: add to main document: " << std::endl;
        texout << "% \\usepackage{siunitx} " << std::endl;
        texout << "% \\sisetup{separate-uncertainty,table-format=6.3(6)}  % hint: modify table-format to best fit your tables" << std::endl;
    }
    if(fTableOptions.find("LANDSCAPE")!=std::string::npos){
        texout << "\\begin{landscape}" << std::endl;
    }
    texout << "\\begin{table}[htbp]" << std::endl;
    texout << "\\begin{center}" << std::endl;
    if(fTableOptions.find("FOOTNOTESIZE")!=std::string::npos){
        texout << "\\footnotesize" << std::endl;
    }
    texout << "\\begin{tabular}{|l" ;
    for(std::size_t ibin = 1; ibin <= regions.size(); ++ibin) {
        texout << "|S";
    }
    texout << "|}" << std::endl;
    texout << "\\hline " << std::endl;
    for (const auto& ireg : regions) {
        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&ireg](const auto& region){return ireg == region->fName;});
        if((*itr)->fTexLabel!="") texout << " & {" << (*itr)->fTexLabel << "}";
        else                      texout << " & {" << (*itr)->fLabel    << "}";
    }
    texout << "\\\\" << std::endl;
    texout << "\\hline " << std::endl;

    // for each unique (unique title), for each region, build a list of samples to consider
    std::vector<std::pair<std::string, std::map<std::string, std::vector<std::string> > > > samplesToConsider;

    for (const auto& isample : fSamples) {
        if (isample->fType==Sample::SampleType::GHOST) continue;
        if (isample->fType==Sample::SampleType::EFT && isample->fEFTSMReference != "") continue;
        if (isample->fType==Sample::SampleType::DATA) continue;
        if (isample->fRegions.empty()) continue;

        auto itr = std::find_if(samplesToConsider.begin(), samplesToConsider.end(), [&isample](const auto& element){return isample->fTitle == element.first;});
        if (itr == samplesToConsider.end()) {
            std::map<std::string, std::vector<std::string> > toConsider;
            for (const auto& ireg : regions) {
                if (isample->fRegions.at(0) == "all" || Common::FindInStringVector(isample->fRegions, ireg) >= 0) {
                    std::vector<std::string> tmp;
                    tmp.emplace_back(isample->fName);
                    toConsider.insert({ireg, tmp});
                }
            }
            samplesToConsider.emplace_back(std::make_pair(isample->fTitle, toConsider));
        } else {
            for (const auto& ireg : regions) {
                auto itrReg = itr->second.find(ireg);
                if (isample->fRegions.at(0) == "all" || Common::FindInStringVector(isample->fRegions, ireg) >= 0) {
                    if (itrReg == itr->second.end()) {
                        std::vector<std::string> tmp;
                        tmp.emplace_back(isample->fName);
                        itr->second.insert({ireg, tmp});
                    } else {
                        itrReg->second.emplace_back(isample->fName);
                    }
                }
            }
        }
    }

    // add tot uncertainty on each sample
    for(const auto& isample : samplesToConsider) {
        //
        // print values
        out << " | " << isample.first << " | ";
        container.sampleNames.emplace_back(isample.first);
        std::vector<double> yamlTmpYields;
        std::vector<double> yamlTmpErrors;

        auto itr = std::find_if(fSamples.begin(), fSamples.end(), [&isample](const auto& element){return isample.first == element->fTitle;});
        if (itr == fSamples.end()) {
            LOG(ERROR) << "Cannot find sample: " << isample.first << "\n";
            exit(EXIT_FAILURE);
        }

        if((*itr)->fTexTitle!="") texout << "  " << (*itr)->fTexTitle << "  ";
        else                      texout << "  " << isample.first << "  ";
        for (const auto& ireg : regions) {
            std::string samples("");
            bool sampleAvailableInRegion(true);
            auto itrReg = isample.second.find(ireg);
            if (itrReg == isample.second.end()) sampleAvailableInRegion = false;
            if (sampleAvailableInRegion) {
                for (const auto& s : itrReg->second) {
                    const std::string name = s + "_" + ireg + "_shapes";
                    samples += samples.empty() ? name : "," + name;
                }
            }

            auto channel = pdf->find(ireg);
            if (!channel) {
                LOG(ERROR) << "Cannot read node: " << ireg << "\n";
                exit(EXIT_FAILURE);
            }

            if (sampleAvailableInRegion) {
                xRooNode sample = FitUtils::GetReducedOrSingleSample(channel, samples);
                double mean = sample.GetContent();
                double uncertainty = sample.GetError();
                double mean_rounded = mean;
                double uncertainty_rounded = uncertainty;
                yamlTmpYields.emplace_back(mean);
                yamlTmpErrors.emplace_back(uncertainty);
                int n = -1; // this will contain the number of decimal places
                if (fUsePDGRoundingTxt || fUsePDGRoundingTex){
                    n = Common::ApplyPDGrounding(mean_rounded, uncertainty_rounded);
                }
                if(fUsePDGRoundingTxt){
                    out << mean_rounded << " pm " << uncertainty_rounded << " | ";
                }
                else{
                    out << mean << " pm " << uncertainty << " | ";
                }
                if (fUsePDGRoundingTex){
                    texout << " & ";
                    if(n<0) texout << mean_rounded;
                    else    texout << Form(("%."+std::to_string(n)+"f").c_str(),mean_rounded);
                    if(uncertainty==0){ // to fix Latex siunitx issue
                        if(n<0) texout << " (" << uncertainty_rounded << ")";
                        else    texout << " (" << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded) << ")";
                    }
                    else{
                        if(n<0) texout << " \\pm " << uncertainty_rounded;
                        else    texout << " \\pm " << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded);
                    }
                }
                else{
                    texout << " & " << mean << " \\pm " << uncertainty;
                }
            } else {
                yamlTmpYields.emplace_back(0);
                yamlTmpErrors.emplace_back(0);
                out << " 0 pm 0 | ";
                texout << " & 0 \\pm 0";
            }
        }
        out << std::endl;
        texout << " \\\\ ";
        texout << std::endl;
        container.mcYields.emplace_back(yamlTmpYields);
        container.mcErrors.emplace_back(yamlTmpErrors);
    }

    out << " | Total | ";
    texout << "\\hline " << std::endl;
    texout << "  Total ";
    container.sampleNames.emplace_back("Total");
    std::vector<double> yamlTmpYields;
    std::vector<double> yamlTmpErrors;
    for (const auto& ireg : regions) {
        auto channel = pdf->find(ireg);

        double mean = channel->GetContent();
        double uncertainty = channel->GetError();
        double mean_rounded = mean;
        double uncertainty_rounded = uncertainty;
        yamlTmpYields.emplace_back(mean);
        yamlTmpErrors.emplace_back(uncertainty);
        int n = -1; // this will contain the number of decimal places
        if (fUsePDGRoundingTxt || fUsePDGRoundingTex){
            n = Common::ApplyPDGrounding(mean_rounded, uncertainty_rounded);
        }
        if (fUsePDGRoundingTxt){
            out << mean_rounded << " pm " << uncertainty_rounded << " | ";
        }
        else{
            out << mean << " pm " << uncertainty << " | ";
        }
        if(fUsePDGRoundingTex){
            texout << " & ";
            if(n<0) texout << mean_rounded;
            else    texout << Form(("%."+std::to_string(n)+"f").c_str(),mean_rounded);
            if(uncertainty==0){ // to fix Latex siunitx issue
                if(n<0) texout << " (" << uncertainty_rounded << ")";
                else    texout << " (" << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded) << ")";
            }
            else{
                if(n<0) texout << " \\pm " << uncertainty_rounded;
                else    texout << " \\pm " << Form(("%."+std::to_string(n)+"f").c_str(),uncertainty_rounded);
            }
        }
        else{
            texout << " & " << mean << " \\pm " << uncertainty;
        }
    }
    container.mcYields.emplace_back(yamlTmpYields);
    container.mcErrors.emplace_back(yamlTmpErrors);
    out << std::endl;
    texout << " \\\\ ";
    texout << std::endl;

    //
    // Print data
    if (!fFitIsBlind) {
        texout << "\\hline " << std::endl;
        std::string title("");
        std::string texTitle("");
        for (const auto& isample : fSamples) {
            if (isample->fType!=Sample::SampleType::DATA) continue;
            title = isample->fTitle;
            if (isample->fTexTitle!="") texTitle = isample->fTexTitle;
            else                        texTitle = isample->fTitle;
        }
        std::vector<double> yamlTempDataYields;
        out << " | " << title << " | ";
        texout << "  " << texTitle << "  ";
        container.dataNames.emplace_back(title);
        for (const auto& ireg : regions) {
            auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&ireg](const auto& region){return ireg == region->fName;});

            std::vector<int> blindedBins = Common::CombineVectors((*itr)->fDropBins, (*itr)->fBlindedBins, true);
            if (blindedBins.empty()) {
                const double data = (*itr)->fData->fHist->Integral();
                yamlTempDataYields.emplace_back(data);
                texout << " & ";
                out << data;
                texout << Form("%.0f",data);
                out << " | ";
            }
            else{
                const double data = -1.0;
                yamlTempDataYields.emplace_back(data);
                texout << " & ";
                out << data;
                texout << Form("%.0f",data);
                out << " | ";
            }
        }
        out << std::endl;
        texout << " \\\\ ";
        texout << std::endl;
        container.dataYields.emplace_back(yamlTempDataYields);
    }
    //
    texout << "\\hline " << std::endl;
    texout << "\\end{tabular} " << std::endl;
    texout << "\\caption{Yields of the analysis} " << std::endl;
    texout << "\\end{center} " << std::endl;
    texout << "\\end{table} " << std::endl;
    if(fTableOptions.find("LANDSCAPE")!=std::string::npos){
        texout << "\\end{landscape}" << std::endl;
    }
    if(fTableOptions.find("STANDALONE")!=std::string::npos){
        texout << "\\end{document}" << std::endl;
    }
    if(fCleanTables){
        std::string prepost_suffix = suffix;
        if (isPostFit) prepost_suffix = suffix + "_postFit";
        std::string shellcommand = "cat "+fName+"/Tables/Yields"+prepost_suffix+".tex|sed -e 's:#:\\\\:g' > "+fName+"/Tables/Yields"+prepost_suffix;
        shellcommand += "_clean.tex";
        gSystem->Exec(shellcommand.c_str());
    }

    // Write YAML
    YamlConverter converter{};
    converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
    converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
    converter.WriteTables(container, fName, isPostFit);

    if (fHEPDataFormat) {
        converter.WriteTablesHEPData(container, fName, isPostFit);
    }
}

//__________________________________________________________________________________
//
void TRExFit::DrawSignalRegionsPlot(int nCols,int nRows) const{
    std::vector<std::string> regions;
    if(fRegionsToPlot.size() > 0) {
        nCols = 1;
        nRows = 1;
        // first loop
        int nRegInRow = 0;
        for(std::size_t i = 0; i < fRegionsToPlot.size(); ++i) {
            LOG(DEBUG) << "Regions to Plot: " << fRegionsToPlot[i] << "\n";
            if(fRegionsToPlot.at(i).find("ENDL") != std::string::npos){
                nRows++;
                if(nRegInRow > nCols) nCols = nRegInRow;
                nRegInRow = 0;
            } else if(fRegionsToPlot.at(i).find("EMPTY") != std::string::npos){
                regions.emplace_back("");
                nRegInRow++;
            } else{
                regions.emplace_back(fRegionsToPlot.at(i));
                nRegInRow++;
            }
        }
    } else {
        for (const auto& ireg : fRegions) {
            regions.push_back(ireg->fName);
        }
    }
    DrawSignalRegionsPlot(nCols,nRows,regions);
}

//__________________________________________________________________________________
//
void TRExFit::DrawSignalRegionsPlot(int nCols,int nRows, const std::vector<std::string>& regions) const{
    gSystem->mkdir(fName.c_str(), true);
    double Hp = 250.; // height of one mini-plot, in pixels
    double Wp = 200.; // width of one mini-plot, in pixels
    double H0 = 100.; // height of the top label pad
    if(TRExFitter::OPTION["SignalRegionSize"]!=0){
        Hp = TRExFitter::OPTION["SignalRegionSize"];
        Wp = (200./250.)*TRExFitter::OPTION["SignalRegionSize"];
    }
    double H = H0 + nRows*Hp; // tot height of the canvas
    double W = nCols*Wp; // tot width of the canvas
    W += 0.1; // to fix eps format (why is this needed?)

    TCanvas c("c","c",W,H);
    TPad pTop("c0","c0",0,1-H0/H,1,1);
    pTop.Draw();
    pTop.cd();
    if (fPlotLabel != "none") TRExLabel(0.1/(W/200.),1.-0.3, TRExFitter::EXPERIMENT_LABEL.c_str(), ("Simulation " + fPlotLabel).c_str());
    myText(    0.1/(W/200.),1.-0.6,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));
    if(fLabel!="-") myText(    0.1/(W/200.),1.-0.9,1,Form("%s",fLabel.c_str()));

    std::unique_ptr<TLegend> leg = nullptr;

    c.cd();

    TPad pLeft("c1","c1",0,0,0+(W-nCols*Wp)/W,1-H0/H);
    pLeft.Draw();
    pLeft.cd();
    TLatex tex0{};
    tex0.SetNDC();
    tex0.SetTextAngle(90);
    tex0.SetTextAlign(23);
    tex0.DrawLatex(0.4,0.5,"S / #sqrt{ B }");

    c.cd();

    TPad pBottom("c1","c1",0+(W-nCols*Wp)/W,0,1,1-H0/H);
    pBottom.Draw();
    pBottom.cd();

    pBottom.Divide(nCols,nRows);
    std::size_t Nreg = static_cast<unsigned int> (nRows*nCols);
    if(Nreg>regions.size()) Nreg = regions.size();
    std::vector<double> S(Nreg);
    std::vector<double> B(Nreg);
    std::vector<double> xbins = {0.0,0.1,0.9,1.0};
    TLatex tex{};
    tex.SetNDC();
    tex.SetTextSize(gStyle->GetTextSize());
    pBottom.cd(1);

    //
    // Get the values
    //
    std::string resultFileAddition("");
    if (fHasValidationRegions){
        resultFileAddition =  "_allRegions";
    }
    else if (fHasDropBinRegions){
        resultFileAddition = "_allBinsFitRegions";
    }

    const std::string wsPath = fName+"/RooStats/"+fInputName+resultFileAddition +"_combined_"+fInputName+fSuffix+"_model.root";
    xRooNode rooNode = Common::ReadWSandFitResults(wsPath, "");

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf\n";
        exit(EXIT_FAILURE);
    }

    for(std::size_t i = 0; i < Nreg; ++i) {
        S.at(i) = 0.;
        B.at(i) = 0.;
        const std::string regName = regions.at(i);
        if (regName.empty()) continue;
        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&regName](const auto& element){return regName == element->fName;});
        if (itr == fRegions.end()) {
            LOG(ERROR) << "Cannot find region: " << regName << "\n";
            exit(EXIT_FAILURE);
        }

        for(const auto& isig : (*itr)->fSig) {
            if(isig == nullptr) continue;
            auto sample = Common::XRooNodeSampleFromNode(pdf, (*itr)->fName, isig->fName);
            S.at(i) += sample->GetContent();
        }
        for(const auto& ibkg : (*itr)->fBkg) {
            if(ibkg == nullptr) continue;
            auto sample = Common::XRooNodeSampleFromNode(pdf, (*itr)->fName, ibkg->fName);
            B.at(i) += sample->GetContent();
        }
        // to avoid nan or inf...
        if(B.at(i)==0) B.at(i) = 1e-10;
        // scale up for projections
        if(fLumiScale!=1){
            S.at(i)*=fLumiScale;
            B.at(i)*=fLumiScale;
        }
    }
    //
    double yMax = 0;
    //
    bool hasSR = false;
    bool hasCR = false;
    bool hasVR = false;

    std::vector<TH1D> h;
    h.reserve(Nreg);

    for(std::size_t i = 0; i < Nreg; ++i) {
        const std::string regName = regions.at(i);
        if (regName.empty()) continue;
        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&regName](const auto& element){return regName == element->fName;});
        if (itr == fRegions.end()) {
            LOG(ERROR) << "Cannot find region: " << regName << "\n";
            exit(EXIT_FAILURE);
        }
        pBottom.cd(i+1);
        if(TRExFitter::OPTION["LogSignalRegionPlot"]) gPad->SetLogy();
        if(TRExFitter::OPTION["LogXSignalRegionPlot"])gPad->SetLogx();
        std::string label = (*itr)->fShortLabel;
        h.emplace_back(Form("h[%ld]",i),label.c_str(),3,xbins.data());
        h.back().SetBinContent(2,S.at(i)/std::sqrt(B.at(i)));
        h.back().GetYaxis()->CenterTitle();
        h.back().GetYaxis()->SetLabelSize( h.back().GetYaxis()->GetLabelSize() * (Wp/200.) );
        h.back().GetXaxis()->SetTickLength(0);
        if(TRExFitter::OPTION["LogSignalRegionPlot"]==0) h.back().GetYaxis()->SetNdivisions(3);
        else TGaxis::SetMaxDigits(5);
        yMax = std::max(yMax,h.back().GetMaximum());
        h.back().GetXaxis()->SetLabelSize(0);
        h.back().SetLineWidth(1);
        h.back().SetLineColor(kBlack);
        if((*itr)->fRegionType==Region::SIGNAL)          h.back().SetFillColor(kRed+1);
        else if((*itr)->fRegionType==Region::VALIDATION) h.back().SetFillColor(kGray);
        else                                             h.back().SetFillColor(kAzure-4);
        if(leg!=nullptr){
            if((*itr)->fRegionType==Region::CONTROL && !hasCR)    {
                leg->AddEntry(&h.back(),"Control Regions","f");
                hasCR = true;
            }
            if((*itr)->fRegionType==Region::VALIDATION && !hasVR) {
                leg->AddEntry(&h.back(),"Validation Regions","f");
                hasVR = true;
            }
            if((*itr)->fRegionType==Region::SIGNAL && !hasSR)     {
                leg->AddEntry(&h.back(),"Signal Regions","f");
                hasSR = true;
            }
        }
        h.back().Draw();
        gPad->SetLeftMargin( gPad->GetLeftMargin()*2.4 );
        gPad->SetRightMargin(gPad->GetRightMargin()*0.1);
        gPad->SetTicky(0);
        gPad->RedrawAxis();
        tex.DrawLatex(0.27,0.85,label.c_str());
        const double SoB = S.at(i)/B.at(i);
        std::string SB = Form("%.1f%%",(100.*SoB));
        SB = "#scale[0.75]{S/B} = "+SB;
        tex.DrawLatex(0.27,0.72,SB.c_str());
    }
    //
    for(std::size_t i = 0; i < Nreg; ++i) {
        if (regions.at(i).empty()) continue;
        if ((h.size() - 1)  <= i) break;
        if(TRExFitter::OPTION["LogSignalRegionPlot"]!=0){
            h.at(i).SetMaximum(yMax*200);
            h.at(i).SetMinimum(2e-4);
        }
        else{
            h.at(i).SetMaximum(yMax*1.5);
            h.at(i).SetMinimum(0.);
        }
    }
    //
    Common::SaveCanvasAs(c, fName+"/SignalRegions"+fSuffix);

}

//__________________________________________________________________________________
//
void TRExFit::DrawPieChartPlot(const std::string &opt, int nCols, int nRows) const{

    std::vector<std::string> regions;
    if(fRegionsToPlot.size()>0){
        nCols = 1;
        nRows = 1;
        // first loop
        int nRegInRow = 0;
        for(std::size_t i = 0; i < fRegionsToPlot.size(); ++i) {
            LOG(DEBUG) << "Regions to plot: " << fRegionsToPlot[i] << "\n";
            if (fRegionsToPlot[i].find("ENDL") != std::string::npos) {
                nRows++;
                if(nRegInRow>nCols) nCols = nRegInRow;
                nRegInRow = 0;
            } else if(fRegionsToPlot.at(i).find("EMPTY") != std::string::npos){
                regions.emplace_back("");
                nRegInRow++;
            } else {
                regions.emplace_back(fRegionsToPlot.at(i));
                nRegInRow++;
            }
        }
    } else {
        for (const auto& ireg : fRegions) {
            regions.push_back(ireg->fName);
        }
    }
    DrawPieChartPlot(opt, nCols,nRows,regions);

}


//__________________________________________________________________________________
//
void TRExFit::DrawPieChartPlot(const std::string &opt, int nCols,int nRows, const std::vector<std::string> &regions) const{

    double Hp = 250.; // height of one mini-plot, in pixels
    double Wp = 250.; // width of one mini-plot, in pixels
    double H0 = 100.; // height of the top label pad

    if(TRExFitter::OPTION["PieChartSize"]!=0){
        Hp = TRExFitter::OPTION["PieChartSize"];
        Wp = TRExFitter::OPTION["PieChartSize"];
    }

    double H = H0 + nRows*Hp; // tot height of the canvas
    double W = nCols*Wp; // tot width of the canvas

    bool isPostFit = opt.find("post")!=std::string::npos;

    //
    // Create the canvas
    //
    if (fPieChartCanvasSize.size() != 0){
        W = fPieChartCanvasSize.at(0);
        H = fPieChartCanvasSize.at(1);
    }
    TCanvas c("c","c",W,H);
    TPad pTop("c0","c0",0,1-H0/H,1,1);
    pTop.Draw();
    pTop.cd();

    if (fPlotLabel != "none") TRExLabel(0.05 / (W/200),0.7, TRExFitter::EXPERIMENT_LABEL.c_str(), ("Simulation " + fPlotLabel).c_str());
    myText(    0.05 / (W/200),0.4,1,Form("#sqrt{s} = %s",fCmeLabel.c_str()));
    if(fLabel!="-") myText(    0.05 / (W/200),0.1,1,Form("%s",fLabel.c_str()));

    c.cd();
    TPad pBottom("c1","c1",0,0,1,1-H0/H);
    pBottom.Draw();
    pBottom.cd();
    pBottom.Divide(nCols,nRows);
    int Nreg = nRows*nCols;
    if(Nreg>(int)regions.size()) Nreg = regions.size();
    TLatex tex{};
    tex.SetNDC();
    tex.SetTextSize(gStyle->GetTextSize());
    pBottom.cd(1);

    //
    // Create the map to store all the needed information
    //
    std::map < std::string, int > map_for_legend;
    std::vector < std::map < std::string, double > > results;
    std::vector < std::map < std::string, int > > results_color;

    std::string resultFileAddition("");
    if (fHasValidationRegions){
        resultFileAddition =  "_allRegions";
    }
    else if (fHasDropBinRegions){
        resultFileAddition = "_allBinsFitRegions";
    }

    std::string frFile("");
    if (fFitResultsRootFile != "") {
        frFile = fFitResultsRootFile;
    } else {
        frFile = fName + "/Fits/"+fInputName+fSuffix+".root";
    }

    const std::string wsPath = fName+"/RooStats/"+fInputName+resultFileAddition +"_combined_"+fInputName+fSuffix+"_model.root";
    xRooNode rooNode = Common::ReadWSandFitResults(wsPath, isPostFit ? frFile : "");

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf\n";
        exit(EXIT_FAILURE);
    }

    //
    // Get the values
    //
    for(int i=0;i<Nreg;i++){
        std::map < std::string, double > temp_map_for_region;
        std::map < std::string, int > temp_map_for_region_color;
        const std::string regName = regions.at(i);
        if (regName.empty()) continue;
        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&regName](const auto& element){return regName == element->fName;});
        if (itr == fRegions.end()) {
            LOG(ERROR) << "Cannot find region: " << regName << "\n";
            exit(EXIT_FAILURE);
        }

        for(int i_bkg = (*itr)->fBkg.size()-1; i_bkg >= 0; --i_bkg) {
            if (!(*itr)->fBkg.at(i_bkg)) continue;
            std::string title = (*itr)->fBkg.at(i_bkg)->fSample->fTitle;
            if((*itr)->fBkg.at(i_bkg)->fSample->fGroup != "") title = (*itr)->fBkg.at(i_bkg)->fSample->fGroup.c_str();

            auto sample = Common::XRooNodeSampleFromNode(pdf, (*itr)->fName, (*itr)->fBkg.at(i_bkg)->fName);
            double integral = sample->GetContent() * fLumiScale;

            if(temp_map_for_region.find(title)!=temp_map_for_region.end()){
                temp_map_for_region[title] += integral;
            } else {
                temp_map_for_region.insert( std::pair < std::string, double > (title,integral) );
                temp_map_for_region_color.insert( std::pair < std::string, int > (title, (*itr)->fBkg.at(i_bkg)->fSample->fFillColor) );
            }
            map_for_legend[title] = (*itr)->fBkg.at(i_bkg)->fSample->fFillColor;
        }
        results.push_back(temp_map_for_region);
        results_color.push_back(temp_map_for_region_color);
    }

    //
    // Finally writting the pie chart
    //
    std::vector<std::unique_ptr<TPie> > pie;
    int validIndex(0);
    for(int i=0;i<Nreg;i++){
        const std::string regName = regions.at(i);
        if (regName.empty()) continue;
        auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&regName](const auto& element){return regName == element->fName;});
        if (itr == fRegions.end()) {
            LOG(ERROR) << "Cannot find region: " << regName << "\n";
            exit(EXIT_FAILURE);
        }

        pBottom.cd(i+1);
        std::string label = (*itr)->fShortLabel;

        const unsigned int back_n = results.at(validIndex).size();
        std::vector<double> values(back_n);
        std::vector<int> colors(back_n);
        for( unsigned int iTemp = 0; iTemp < back_n; ++iTemp ){
            values.at(iTemp) = 0.;
            colors.at(iTemp) = 0;
        }

        int count = 0;
        for (const auto& temp_pair : results.at(validIndex) ){
            values.at(count) = temp_pair.second;
            colors.at(count) = results_color.at(validIndex)[temp_pair.first];
            count++;
        }

        pie.emplace_back(std::make_unique<TPie>(("pie_"+label).c_str()," ",back_n, values.data(), colors.data()));
        pie.back()->SetRadius( pie.back()->GetRadius() * 0.8 );
        for(int iEntry = 0; iEntry < pie.back()->GetEntries(); ++iEntry) {
            pie.back()->SetEntryLabel(iEntry,"");
        }
        pie.back()->Draw();
        tex.DrawLatex(0.1,0.85,label.c_str());
        validIndex++;
    }

    c.cd();

    //
    // Adding the legend in the top panel
    //
    pTop.cd();
    std::unique_ptr<TLegend> leg = std::make_unique<TLegend>(0.7,0.1,0.95,0.90);

    leg->SetLineStyle(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetBorderSize(0);
    leg->SetTextFont( gStyle->GetTextFont() );
    leg->SetTextSize( gStyle->GetTextSize() );

    std::vector<std::string> legVec;
    for ( const std::pair < std::string, int > legend_entry : map_for_legend ) {
        legVec.push_back(legend_entry.first);
    }
    std::vector<std::unique_ptr<TH1D> > dummy;
    for(int i_leg=legVec.size()-1;i_leg>=0;i_leg--){
        dummy.emplace_back(std::make_unique<TH1D>(("legend_entry_" + legVec[i_leg]).c_str(), "",1,0,1));
        dummy.back()->SetFillColor(map_for_legend[legVec.at(i_leg)]);
        dummy.back()->SetLineColor(kBlack);
        dummy.back()->SetLineWidth(1);
        leg->AddEntry(dummy.back().get(),legVec.at(i_leg).c_str(),"f");
    }
    leg->Draw();

    //
    // Stores the pie chart in the desired format
    //
    Common::SaveCanvasAs(c, fName+"/PieChart" + fSuffix + (isPostFit ? "_postFit" : ""));
}

//__________________________________________________________________________________
// called before w in case of CustomAsimov
void TRExFit::CreateCustomAsimov() const {
    LOG(DEBUG) << "Running CreateCustomAsimov\n";
    // get a list of all CustomAsimov to create
    std::vector<std::string> customAsimovList;
    for(const auto& isample : fSamples) {
        for( auto asimovreplacements: isample->fAsimovReplacementFor){
            if(asimovreplacements.first!="" && Common::FindInStringVector(customAsimovList,asimovreplacements.first)<0) {
                customAsimovList.push_back(asimovreplacements.first);
            }
        }
    }
    //
    // fill a different CustomAsimov data-set for each element in the list
    for(const auto& customAsimov : customAsimovList){
        LOG(DEBUG) << "CustomAsimov: " << customAsimov << "\n";
        std::shared_ptr<Sample> ca = GetSample("customAsimov_"+customAsimov);
        // create a new data sample taking the nominal S and B
        for(const auto& ireg : fRegions) {
            // Now we need to clone a histogram, but need to find one that is valid
            std::shared_ptr<SampleHist> sample_hist = ireg->fData;
            if (!sample_hist) {
                // try to clone signal
                for (const auto& isig : ireg->fSig) {
                    if (isig != nullptr) {
                        sample_hist = isig;
                        break;
                    }
                }
            }
            if (!sample_hist) {
                // try to clone background
                for (const auto& ibkg : ireg->fBkg) {
                    if (ibkg != nullptr) {
                        sample_hist = ibkg;
                        break;
                    }
                }
            }

            if (!sample_hist) {
                LOG(ERROR) << "Cannot copy a valid sample hist!\n";
                exit(EXIT_FAILURE);
            }
            std::shared_ptr<SampleHist> cash = ireg->SetSampleHist(ca.get(),static_cast<TH1*>(sample_hist->fHist->Clone()));

            cash->fHist_orig->SetName( Form("%s_orig",cash->fHist->GetName()) ); // fix the name
            cash->fHist->Scale(0.);
            //
            std::vector<std::string> smpToExclude;

            for(const auto& isample : fSamples) {
                std::shared_ptr<SampleHist> h = ireg->GetSampleHist(isample->fName);
                if(!h) continue;
                if(h->fSample->fType==Sample::SampleType::DATA) continue;
                if(h->fSample->fType==Sample::SampleType::GHOST || h->fSample->fType==Sample::SampleType::EFT ) {
                    bool continue_sample = true;
                    for( auto asimovreplacements: h->fSample->fAsimovReplacementFor){
                        if(asimovreplacements.first!=customAsimov){
                            continue;
                        }
                        //can probably get rid of the if statement now
                        if(asimovreplacements.second!="" ){
                            smpToExclude.push_back(asimovreplacements.second);
                            continue_sample = false;
                        }
                    }
                    if(continue_sample){
                        continue;
                    }
                }
                if( Common::FindInStringVector(smpToExclude,isample->fName) >= 0 ) continue;
                if (fFitType == TRExFit::FitType::EFT){
                    if(h->fSample->fName.rfind("SM_", 0) == 0 && h->fSample->fName.find("_bin") != std::string::npos){
                        std::string eft_samplename = h->fSample->fName;
                        eft_samplename = eft_samplename.substr(3);
                        const size_t found = eft_samplename.find(ireg->fName);
                        if (found == std::string::npos || found == 0) {
                            LOG(ERROR) << "ProblemParsing EFT name: " << eft_samplename << "\n";
                            return;
                        }
                        eft_samplename.resize(found-1);
                        if( Common::FindInStringVector(smpToExclude,eft_samplename) >= 0 ) continue;
                    }
                }
                //
                // bug-fix: change normalisation factors to nominal value!
                double factor = 1.;
                for(const auto& norm : isample->fNormFactors) {
                    if (std::find(norm->fRegions.begin(), norm->fRegions.end(), ireg->fName) != norm->fRegions.end()) {
                        LOG(DEBUG) << "Setting norm factor to " << norm->GetNominal() << "\n";
                        factor *= norm->GetNominal();
                    }
                }
                LOG(INFO) << "Creating asimov dataset for " << customAsimov << ", adding sample " << h->fSample->fName << "\n";
                //
                cash->fHist->Add(h->fHist.get(),factor);
            }
            cash->fHist->Sumw2(false);
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::UnfoldingAlternativeAsimov() {
    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        UnfoldingAlternativeAsimov(fUnfolding.at(i).get());
    }
}
//__________________________________________________________________________________
//
void TRExFit::UnfoldingAlternativeAsimov(const Unfolding* unfolding) {
    if (fFitType != TRExFit::FitType::UNFOLDING) return;
    if (unfolding->fAlternativeAsimovTruthSample == "") return;

    LOG(INFO) << "Replacing data with alternative asimov\n";

    // loop over regions
    for (const auto& ireg : fRegions) {
        std::shared_ptr<SampleHist> sh = ireg->fData;

        // Try getting signal
        if (!sh) {
            for (const auto& isig : ireg->fSig) {
                if (isig) {
                    sh = isig;
                    break;
                }
            }
        }

        if (!sh) {
            LOG(ERROR) << "No data or signal found, this should not happen!\n";
            exit(EXIT_FAILURE);
        }

        // now replace fData
        TH1* hist = static_cast<TH1*>(sh->fHist->Clone());
        hist->Reset();

        // Add new signal
        std::shared_ptr<Sample> newAsimov = GetSample("AlternativeSignal_"+unfolding->fName+"_"+ireg->fName);
        if (!newAsimov) {
            LOG(ERROR) << "Cannot read the new asimov sample\n";
            exit(EXIT_FAILURE);
        }
        std::shared_ptr<SampleHist> newsh = ireg->GetSampleHist(newAsimov->fName);
        if (!newsh) {
            LOG(ERROR) << "Cannot read the new asimov SampleHist\n";
            exit(EXIT_FAILURE);
        }
        hist->Add(newsh->fHist.get());

        // add bkgs bkg
        for (const auto& isample : fSamples) {
            std::shared_ptr<SampleHist> sampleHist = ireg->GetSampleHist(isample->fName);
            if (!sampleHist) continue;
            if (sampleHist->fSample->fType != Sample::SampleType::BACKGROUND) continue;
            if(Common::FindInStringVector(isample->fRegions,ireg->fName) <0) continue;
            hist->Add(sampleHist->fHist.get());
        }

        hist->Sumw2(false);
        ireg->fData->fHist.reset(static_cast<TH1*>(hist->Clone()));
    }
}

//__________________________________________________________________________________
// turn to RooStats::HistFactory
void TRExFit::ToRooStats(const bool allRegions, const bool allBinsFitRegions) const {

    // need to check if we have some valid regions for the actual fit
    if (!allRegions) {
        bool hasNonValidation(false);
        for (const auto& ireg : fRegions) {
            if (ireg->fRegionType != Region::RegionType::VALIDATION) {
                hasNonValidation = true;
                break;
            }
        }
        if (!hasNonValidation) {
            LOG(WARNING) << "Only validation regions are provided, skipping the generation of the WS for fitting\n";
            return;
        }
    }

    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Exporting to RooStats...\n";

    if (fPOIs.empty()){
        LOG(ERROR) << "The configuration file is missing the declaration of a POI name in the JOB block. This breaks RooStats.\n";
        exit(EXIT_FAILURE);
    }

    RooStats::HistFactory::Measurement meas((fInputName+fSuffix).c_str(), (fInputName+fSuffix).c_str());

    std::string wsDir = fName+"/RooStats/";
    std::string fileName = wsDir+fInputName;
    if(fBootstrap!="" && fBootstrapIdx>=0) {
        wsDir += fBootstrapSyst+fBootstrapSample+"_BSId"+Form("%d",fBootstrapIdx) + "/";
        gSystem->mkdir(wsDir.c_str(), true);
        fileName = wsDir+fInputName;
    }


    if (allRegions) {
        LOG(INFO) << "Creating workspace using all regions (including validation regions) and all bins\n";
        LOG(DEBUG) << "Output workspace path will have suffix _allRegions...\n";
        meas.SetOutputFilePrefix((fileName+"_allRegions").c_str());
    }
    else if (allBinsFitRegions){
        LOG(INFO) << "Creating workspace only with fit regions and all bins (no validation regions)\n";
        LOG(DEBUG) << "Output workspace path will have suffix _allBinsFitRegions...\n";
        meas.SetOutputFilePrefix((fileName+"_allBinsFitRegions").c_str());
    }
    else {
        LOG(INFO) << "Creating workspace only with fit regions and dropped bins (no validation regions)\n";
        LOG(DEBUG) << "Output workspace path will have no special suffix ...\n";
        meas.SetOutputFilePrefix(fileName.c_str());
    }


    LOG(INFO) << "-------------------------------------------\n";

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);

    #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
    meas.SetExportOnly(true);
    #endif

    // insert POIs (in reverse order, since the POI is added "at beginning of vector of PoIs")
    for(int i_poi=fPOIs.size()-1;i_poi>=0;i_poi--){
        meas.SetPOI(fPOIs[i_poi].c_str());
    }
    meas.SetLumi(fLumiScale);
    // TRExFitter does not use the HistFactory built-in Lumi parameter with uncertainty LumiRelErr,
    // instead luminosity uncertainties can be implemented as OVERALL type systematic uncertainties
    meas.AddConstantParam("Lumi");
    meas.SetLumiRelErr(0.0);

    std::vector<std::string> regionsUsedInFit;
    for(std::size_t ich = 0; ich < fRegions.size(); ++ich) {

        if(!allRegions && fRegions.at(ich)->fRegionType==Region::VALIDATION) continue;

        LOG(DEBUG) << "Adding Channel: " << fRegions.at(ich)->fName << "\n";
        regionsUsedInFit.emplace_back(fRegions.at(ich)->fName);
        const RooStats::HistFactory::Channel chan = OneChannelToRooStats(&meas, ich, allRegions, allBinsFitRegions);

        meas.AddChannel(chan);
    }
    // Experimental: turn off constraints for given systematics
    for(const auto& isyst : fSystematics) {
        if(isyst->fIsFreeParameter) meas.AddUniformSyst(isyst->fNuisanceParameter.c_str());
    }

    // morphing
    for(const TRExFit::TemplateWeight& itemp : fTemplateWeightVec){
        const std::string normName = "morph_"+itemp.name+"_"+Common::ReplaceString(std::to_string(itemp.value),"-","m");
        LOG(DEBUG) << "Morphing: normName: " << normName << "\n";
        meas.AddPreprocessFunction(normName, itemp.function, itemp.range);
    }
    for(const auto& nf : fNormFactors){
        if(nf->fExpression.first!=""){
            meas.AddPreprocessFunction(nf->fName,nf->fExpression.first,nf->fExpression.second);
        }
    }

    if (!allRegions) meas.PrintXML((wsDir).c_str());

    meas.CollectHistograms();
    if (TRExFitter::DEBUGLEVEL > 1) {
        meas.PrintTree();
    }

    RooStats::HistFactory::HistoToWorkspaceFactoryFast::Configuration cfg{};
    cfg.createPerRegionWorkspaces = fProducePerRegionWS;
    cfg.createWorkspaceFile = false;
    std::unique_ptr<RooWorkspace> work(RooStats::HistFactory::MakeModelAndMeasurementFast(meas, cfg));

    const RooArgSet* externalConstraints = FitUtils::GetExternalConstraints(work.get(), fNormFactors, fRegularizationType);
    if (externalConstraints) {
        work->import(*externalConstraints);
        RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(work->obj("ModelConfig"));
        if (!mc) {
            LOG(ERROR) << "Cannot retrieve the Model Config\n";
            exit(EXIT_FAILURE);
        }
        mc->SetExternalConstraints(*externalConstraints);
    }

    if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();

    if (allRegions) fileName += "_allRegions";
    else if (allBinsFitRegions) fileName += "_allBinsFitRegions";
    else fileName += "";
    fileName += "_combined_"+fInputName+fSuffix+"_model.root";

    this->AddWSMetadata(work.get(), allRegions);

    if (fShapeFactorReparametrisation) {
        FitUtils::ReparametrizeShapeFactors(work.get(), fShapeFactors, this->GetRawRegionVector(), regionsUsedInFit);
    }

    work->writeToFile(fileName.c_str());

    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "UPDATE"));
    if (!file) {
        LOG(ERROR) << "Cannot open file to add the measurement object\n";
        exit(EXIT_FAILURE);
    }

    file->cd();
    meas.writeToFile(file.get());

    file->Close();
}

//__________________________________________________________________________________
//
RooStats::HistFactory::Channel TRExFit::OneChannelToRooStats(RooStats::HistFactory::Measurement* meas,
                                                             const int i_ch,
                                                             const bool allRegions,
                                                             const bool allBinsFitRegions) const {
    RooStats::HistFactory::Channel chan(fRegions[i_ch]->fName.c_str());


    //Suffix used for the regular bin transformed histogram
    std::string suffix_regularBinning = (allRegions||allBinsFitRegions) ? "" : "_regBin";
    if (! (allRegions||allBinsFitRegions)){
        const std::vector<int>& droppedBins = fRegions[i_ch]->GetAutomaticDropBins() ?
                                              fRegions[i_ch]->fComputedBlindedBins : fRegions[i_ch]->fDropBins;
        if (!droppedBins.empty()){
            suffix_regularBinning = "_dropBin";
        }
    }

    //Checks if a data sample exists
    bool hasData = false;
    for(const auto& isample : fSamples) {
        if(isample->fType == Sample::SampleType::DATA){
            hasData = true;
            break;
        }
    }

    if(fCustomAsimov != "") {
        const std::string name = "customAsimov_"+fCustomAsimov;
        std::shared_ptr<SampleHist> cash = fRegions[i_ch]->GetSampleHist(name);
        if(cash==nullptr){
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.clear();
            LOG(WARNING) << "No Custom Asimov " << fCustomAsimov << " available. Taking regular Asimov.\n";
            if (TRExFitter::DEBUGLEVEL < 2) std::cout.setstate(std::ios_base::failbit);
        } else{
            const std::string folder = TRExFitter::USEFOLDERSTRUCTURE ? fRegions[i_ch]->fName + "/" + cash->fSample->fName + "/nominal/" : "" ;
            const std::string temp_string = cash->fHist->GetName();
            LOG(DEBUG) << "  Adding Custom-Asimov Data: " << temp_string << "\n";
            chan.SetData(cash->fHistoName+suffix_regularBinning, cash->fFileName, folder);
        }
    } else if(hasData) {
        if (!fRegions[i_ch]->fData || !fRegions[i_ch]->fData->fSample) {
            LOG(ERROR) << "Data sample not provided for Region: " << fRegions[i_ch]->fName << ". Provide Data sample for all regions, or remove it completely\n";
            exit(EXIT_FAILURE);
        }
        const std::string folder = TRExFitter::USEFOLDERSTRUCTURE ? fRegions[i_ch]->fName + "/" + fRegions[i_ch]->fData->fSample->fName + "/nominal/" : "" ;
        if (!fRegions[i_ch]->fData->fHist) {
            LOG(ERROR) << "Data sample not provided for Region: " << fRegions[i_ch]->fName << ". Provide Data sample for all regions, or remove it completely\n";
            exit(EXIT_FAILURE);
        } else {
            const std::string temp_string = fRegions[i_ch]->fData->fHist->GetName();
            LOG(DEBUG) << "  Adding Data: " << temp_string << "\n";
            chan.SetData(fRegions[i_ch]->fData->fHistoName+suffix_regularBinning, fRegions[i_ch]->fData->fFileName, folder);
        }
    } else {
        chan.SetData("", "");
    }

    // fStatErrCons is upper case after config reading if the MCstatThreshold option is used, otherwise it defaults to "Poisson"
    // HistFactory expects the constraint not in all uppercase, but in form "Poisson"/"Gaussian" instead
    if(fStatErrCons=="Poisson" || fStatErrCons=="POISSON") chan.SetStatErrorConfig(fStatErrThres, "Poisson");
    else if(fStatErrCons=="GAUSSIAN")                      chan.SetStatErrorConfig(fStatErrThres, "Gaussian");

    // Here we can check if the total nominal Histogram has negative bins:
    std::unique_ptr<TH1> hTot =  fRegions.at(i_ch)->GetTotHist(true);
    if (!hTot) {
        LOG(ERROR) << "Cannot read the total histogram for region: " << fRegions.at(i_ch)->fName << "\n";
        exit(EXIT_FAILURE);
    }
    for (int ibin = 1; ibin <= hTot->GetNbinsX(); ++ibin) {
        if (hTot->GetBinContent(ibin) < 0) {
            LOG(ERROR) << "Total nominal histogram in the region" << fRegions.at(i_ch)->fName << " has negative bins. This is not supported by HistFactory - Consider rebinning your histogram.\n";
            exit(EXIT_FAILURE);
        }
    }

    for(std::size_t i_smp = 0; i_smp < fSamples.size(); ++i_smp) {
        std::shared_ptr<SampleHist> h = fRegions.at(i_ch)->GetSampleHist(fSamples.at(i_smp)->fName);
        if (!h) continue;
        if (h->fSample->fType == Sample::SampleType::DATA) continue;
        if (h->fSample->fType == Sample::SampleType::GHOST) continue;
        if (h->fSample->fType == Sample::SampleType::EFT) continue;

        LOG(DEBUG) << "  Adding Sample: " << fSamples.at(i_smp)->fName << "\n";

        const RooStats::HistFactory::Sample sample = OneSampleToRooStats(meas, h.get(), i_ch, i_smp, allRegions, allBinsFitRegions);
        chan.AddSample(sample);
    }

    return chan;
}

//__________________________________________________________________________________
//
RooStats::HistFactory::Sample TRExFit::OneSampleToRooStats(RooStats::HistFactory::Measurement* meas,
                                                           const SampleHist* h,
                                                           const int i_ch,
                                                           const int i_smp,
                                                           const bool allRegions,
                                                           const bool allBinsFitRegions) const {
    RooStats::HistFactory::Sample sample(fSamples[i_smp]->fName.c_str());

    std::string suffix_regularBinning = (allRegions||allBinsFitRegions) ? "" : "_regBin";
    if (! (allRegions||allBinsFitRegions)){
        const std::vector<int>& droppedBins = fRegions[i_ch]->GetAutomaticDropBins() ?
                                               fRegions[i_ch]->fComputedBlindedBins : fRegions[i_ch]->fDropBins;

        if (!droppedBins.empty()){
            suffix_regularBinning = "_dropBin";
        }
    }
    const std::string folder = TRExFitter::USEFOLDERSTRUCTURE ? fRegions[i_ch]->fName + "/" + h->fSample->fName + "/nominal/" : "" ;

    if(fUseStatErr && fSamples[i_smp]->fUseMCStat) sample.ActivateStatError();
    sample.SetHistoName(h->fHistoName+suffix_regularBinning);
    sample.SetHistoPath(folder);
    sample.SetInputFile(h->fFileName);
    sample.SetNormalizeByTheory(fSamples[i_smp]->fNormalizedByTheory);
    // norm factors
    for(const auto& inorm : h->fSample->fNormFactors) {

        if (Common::FindInStringVector(inorm->fExclude, fRegions[i_ch]->fName) >= 0) continue;
        if ((inorm->fRegions.size() > 0) && Common::FindInStringVector(inorm->fRegions, fRegions[i_ch]->fName) < 0) continue;

        LOG(DEBUG) << "    Adding NormFactor: " << inorm->fName << ", " << inorm->GetNominal() << "\n";
        sample.AddNormFactor(inorm->fName,
                             inorm->GetNominal(),
                             inorm->GetMin(),
                             inorm->GetMax());
        if (inorm->fConst) meas->AddConstantParam(inorm->fName);
        if (fStatOnly && fFixNPforStatOnlyFit && Common::FindInStringVector(fPOIs,inorm->fName)<0) {
            meas->AddConstantParam(inorm->fName);
        }
    }

    // shape factors
    for(const auto& ishape : h->fSample->fShapeFactors) {
        if (Common::FindInStringVector(ishape->fExclude, fRegions[i_ch]->fName) >= 0) continue;
        if ((ishape->fRegions.size() > 0) && Common::FindInStringVector(ishape->fRegions, fRegions[i_ch]->fName) < 0) continue;

        LOG(DEBUG) << "    Adding ShapeFactor: " << ishape->fName << ", " << ishape->fNominal << "\n";
        sample.AddShapeFactor(ishape->fName);
        if (ishape->fConst
            || (fStatOnly && fFixNPforStatOnlyFit && Common::FindInStringVector(fPOIs,ishape->fName)<0) ) {
            for(int i_bin=0; i_bin < ishape->fNbins; ++i_bin) {
                meas->AddConstantParam( "gamma_" + ishape->fName + "_bin_" + std::to_string(i_bin) );
            }
        }
    }

    if (fStatOnly) {
        sample.AddOverallSys( "Dummy",1,1 );
        return sample;
    }

    const bool useGaussianShapeSysConstraint = h->fSample->fUseGaussianShapeSysConstraint;
    // systematics
    for(std::size_t i_syst = 0; i_syst < h->fSyst.size(); ++i_syst) {
        std::string temp_systname = h->fSyst[i_syst]->fName;
        if(h->fSyst[i_syst]->fSystematic->fIsShapeAccDecorr){
            if(h->fSyst[i_syst]->fSystematic->fIsNormOnly)          temp_systname = Common::ReplaceString(h->fSyst[i_syst]->fName, "_Acc", "");
            else if(h->fSyst[i_syst]->fSystematic->fIsShapeOnly)    temp_systname = Common::ReplaceString(h->fSyst[i_syst]->fName, "_Shape", "");
        }

        if(h->fSyst[i_syst]->fSystematic->fIsRegionDecorr){
            temp_systname = Common::ReplaceString(h->fSyst[i_syst]->fName, "_" + fRegions[i_ch]->fName, "");
        }

        if(h->fSyst[i_syst]->fSystematic->fIsSampleDecorr){
            temp_systname = Common::ReplaceString(h->fSyst[i_syst]->fName, "_" + fSamples[i_smp]->fName, "");
        }

        const std::string prefix = TRExFitter::USEFOLDERSTRUCTURE ? fRegions[i_ch]->fName + "/" + h->fSample->fName + "/" + temp_systname + "/" : "" ;
        // add normalization part
        LOG(DEBUG) << "    Adding Systematic: " << h->fSyst[i_syst]->fName << "\n";
        if ( h->fSyst[i_syst]->fSystematic->fType==Systematic::SHAPE){
            std::string npName = "shape_" + h->fSyst[i_syst]->fSystematic->fNuisanceParameter+"_";
            std::string regionName = fRegions[i_ch]->fName;
            if(h->fSyst[i_syst]->fSystematic->fNuisanceParameter.find("stat_")!=std::string::npos) {
                // see if there are regions to correlate with others
                for(const auto& set : h->fSample->fCorrelateGammasInRegions){
                    for(std::size_t i_reg = 0; i_reg < set.size(); ++i_reg){
                        if(i_reg != 0 && regionName == set[i_reg]) {
                            regionName = set[0];
                            break;
                        }
                    }
                }
                // eventually correlate MC stat with other samples
                if(h->fSample->fCorrelateGammasWithSample != ""){
                    npName = "shape_stat_"+h->fSample->fCorrelateGammasWithSample+"_";
                }
            }
            npName += regionName;
            if (useGaussianShapeSysConstraint) {
                sample.AddShapeSys( npName,
                                    RooStats::HistFactory::Constraint::Gaussian,
                                    h->fSyst[i_syst]->fHistoNameUp+"_Var"+suffix_regularBinning,
                                    h->fSyst[i_syst]->fFileNameUp,
                                    prefix);
            } else {
                sample.AddShapeSys( npName,
                                    RooStats::HistFactory::Constraint::Poisson,
                                    h->fSyst[i_syst]->fHistoNameUp+"_Var"+suffix_regularBinning,
                                    h->fSyst[i_syst]->fFileNameUp,
                                    prefix);

            }
        } else {
            if ( !h->fSyst[i_syst]->fSystematic->fIsShapeOnly &&
                !h->fSyst[i_syst]->fNormPruned
              ) {
                sample.AddOverallSys( h->fSyst[i_syst]->fSystematic->fNuisanceParameter,
                                      1+h->fSyst[i_syst]->fNormDown,
                                      1+h->fSyst[i_syst]->fNormUp);
            }
            // eventually add shape part
            if ( h->fSyst[i_syst]->fHasShape &&
                !h->fSyst[i_syst]->fSystematic->fIsNormOnly &&
                !h->fSyst[i_syst]->fShapePruned
              ){
                sample.AddHistoSys( h->fSyst[i_syst]->fSystematic->fNuisanceParameter,
                                    h->fSyst[i_syst]->fHistoNameShapeDown+suffix_regularBinning,
                                    h->fSyst[i_syst]->fFileNameShapeDown,
                                    prefix,
                                    h->fSyst[i_syst]->fHistoNameShapeUp+suffix_regularBinning,
                                    h->fSyst[i_syst]->fFileNameShapeUp,
                                    prefix);
            }
        }
    }

    return sample;
}

//__________________________________________________________________________________
//
void TRExFit::SystPruning() {
    LOG(INFO) << "------------------------------------------------------\n";
    LOG(INFO) << "Apply Systematics Pruning ...\n";
    if (fPruningShapeOption == PruningUtil::SHAPEOPTION::KSTEST) {
        LOG(INFO) << "Will run KS test to determine the shape pruning. This is slow compared to the default option (MAXBIN). Patience young padawan.\n";
    }
    if(fSystematics.size()==0 || fStatOnly){
        LOG(INFO) << "No systematics => No Pruning applied.\n";
        return;
    }

    // if COMBINEDSIGNAL is used (works on separate gammas) revert to sample by sample for other systematics
    int strategy = static_cast<int>(fPruningType);
    if (fPruningType == TRExFit::COMBINEDSIGNAL) strategy = 0;

    PruningUtil pu{};
    pu.SetShapeOption(fPruningShapeOption);
    pu.SetStrategy(strategy);
    pu.SetThresholdNorm(fThresholdSystPruning_Normalisation);
    pu.SetThresholdShape(fThresholdSystPruning_Shape);
    pu.SetThresholdIsLarge(fThresholdSystLarge);
    pu.SetRemoveSystOnEmptySample(fRemoveSystOnEmptySample);

    std::vector<std::pair<std::string, double> > impact;
    for(auto& reg : fRegions){
        // if want to skip validation regions from pruning, add a condition here
        const std::vector<std::pair<std::string, double> >& tmp = reg->SystPruning(&pu);
        for (const auto& iimpact : tmp) {
            auto itr = std::find_if(impact.begin(), impact.end(), [&iimpact](const std::pair<std::string, double>& pair){return iimpact.first == pair.first;});
            if (itr == impact.end()) {
                // doesnt exist, add it
                impact.emplace_back(std::move(iimpact));
            } else {
                // exists, check if it is larger in this region
                if (itr->second < iimpact.second) {
                    itr->second = iimpact.second;
                }
            }
        }
    }

    // sort the impact
    std::sort(impact.begin(), impact.end(), [](const std::pair<std::string,double>& lhs, const std::pair<std::string,double>& rhs){return lhs.second < rhs.second;});

    // Draw plot with normalisation effect of each systematic
    if (fDoSystNormalizationPlots) {
        DrawSystematicNormalisationSummary();
    }

    // it also writes the txt file actually
    DrawPruningPlot();

    // drop the bad uncertianties if this is requested
    if (fRemoveLargeSyst) {
        for (const auto& ireg : fRegions) {
            for (const auto& isample : ireg->fSampleHists) {
                for (const auto& isyst : isample->fSyst) {
                    if (isyst->fBadNorm)  isyst->fNormPruned = true;
                    if (isyst->fBadShape) isyst->fShapePruned = true;
                }
            }
        }
    }

    // dumpe the impact table t oa text file
    std::ofstream file(fName+"/NPorder.txt");
    if (!file.good() || !file.is_open()) {
        LOG(ERROR) << "Cannot open file: " << fName << "/NPorder.txt\n";
        return;
    }

    for (const auto& i : impact) {
        file << i.first << " " << i.second << "\n";
    }
    file.close();
}

//__________________________________________________________________________________
//
void TRExFit::DrawSystematicNormalisationSummary() const {


    const std::vector<std::string> uniqueSysts = GetUniqueSystNamesWithoutGamma();
    const std::vector<Region*> regions = GetNonValidationRegions();
    const std::vector<std::shared_ptr<Sample> > samples = GetNonDataNonGhostSamples();

    const int xBins = regions.size()*samples.size();

    TH2D histo("", "", xBins, 0, xBins, 2*uniqueSysts.size(), 0, 2*uniqueSysts.size());
    histo.GetXaxis()->SetTickSize(0);
    histo.GetYaxis()->SetTickSize(0);
    for (std::size_t ireg = 0; ireg < regions.size(); ++ireg) {
        for (std::size_t ibin = 0; ibin < samples.size(); ++ibin) {
            histo.GetXaxis()->SetBinLabel(ireg*samples.size()+ibin+1,
                                          (regions.at(ireg)->fName+"_"+samples.at(ibin)->fName).c_str());
        }
    }
    for (std::size_t isyst = 0; isyst < uniqueSysts.size(); ++isyst) {
        histo.GetYaxis()->SetBinLabel(2*isyst+1, (uniqueSysts.at(isyst)+"_Up").c_str());
        histo.GetYaxis()->SetBinLabel(2*isyst+2, (uniqueSysts.at(isyst)+"_Dn").c_str());
    }
    histo.GetXaxis()->LabelsOption("v");

    // Fill the histograms
    for (std::size_t ireg = 0; ireg < regions.size(); ++ireg) {
        for (std::size_t isample = 0; isample < samples.size(); ++isample) {
            std::shared_ptr<SampleHist> sh = regions.at(ireg)->GetSampleHist(samples.at(isample)->fName);
            for (std::size_t isyst = 0; isyst < uniqueSysts.size(); ++isyst) {
                if (!sh) {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                    continue;
                }

                // actually calcualte the normalisation effect
                const TH1* nominal = sh->fHist.get();
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic(uniqueSysts.at(isyst));
                if (!syh) {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                    continue;
                }

                if (!nominal) {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                    continue;
                }

                const TH1* up   = syh->fHistUp.get();
                const TH1* down = syh->fHistDown.get();

                if (up) {
                    const double norm = 100*(Common::EffIntegral(up) - Common::EffIntegral(nominal))/Common::EffIntegral(nominal);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        norm);
                } else {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+1,
                                        0);
                }

                if (down) {
                    const double norm = 100*(Common::EffIntegral(down) - Common::EffIntegral(nominal))/Common::EffIntegral(nominal);
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        norm);
                } else {
                    histo.SetBinContent(ireg*samples.size()+isample+1,
                                        2*isyst+2,
                                        0);
                }
            }
        }
    }

    static const int upSize = 50;
    static const int loSize = 150;
    static const int leftSize = 250;
    static const int separation = 10;
    const int regionSize = 20*samples.size();
    const int mainHeight = 2*uniqueSysts.size()*20;
    const int mainWidth = regions.size()*(regionSize+separation);

    TCanvas c("","", leftSize+mainWidth, upSize+mainHeight+loSize);
    c.SetTopMargin(0.05);
    c.SetBottomMargin(0.4);
    c.SetLeftMargin(0.4);
    c.SetRightMargin(0.0);
    gStyle->SetPalette(87);
    c.SetGrid();

    #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
    histo.SetMarkerSize(450);
    #else
    histo.SetMarkerSize(1.0);
    #endif
    histo.GetXaxis()->SetLabelOffset(0.3*histo.GetXaxis()->GetLabelOffset());
    gStyle->SetPaintTextFormat(".1f");
    histo.Draw("col TEXT");
    c.RedrawAxis("g");

    Common::SaveCanvasAs(c, fName+"/NormalisationPlot"+fSuffix);
}

//__________________________________________________________________________________
//
void TRExFit::DrawPruningPlot() const{
    //
    std::ofstream out;
    out.open((fName+"/PruningText.txt").c_str());
    out << "-------///////                 ///////-------" << std::endl ;
    out << "-------/////// IN PRUNING PLOT ///////-------" << std::endl ;
    out << "-------///////                 ///////-------" << std::endl ;
    //
    std::vector< std::unique_ptr<TH2F> > histPrun;
    std::vector< std::unique_ptr<TH2F> > histPrun_toSave;
    // make a list of non-data, non-ghost samples
    std::vector<std::vector<std::shared_ptr<Sample> > > samplesVec;
    int iMorph = 0;
    for (const auto& ireg : fRegions) {
        std::vector<std::shared_ptr<Sample> > samples;
        for(const auto& isample : fSamples) {
            if(isample->fType==Sample::SampleType::DATA) continue;
            if(isample->fType==Sample::SampleType::GHOST) continue;
            if(isample->fType==Sample::SampleType::EFT) continue;
            if (iMorph > 0 && !isample->fMorphValue.empty()) continue;
            if (!ireg->fPruningPlotSamples.empty() &&
                (std::find(ireg->fPruningPlotSamples.begin(), ireg->fPruningPlotSamples.end(), isample->fName) == ireg->fPruningPlotSamples.end())) continue;
            if (!isample->fMorphValue.empty()) ++iMorph;
            samples.emplace_back(isample);
        }
        samplesVec.emplace_back(std::move(samples));
    }
    // make a list of non-gamma systematics only
    std::vector<std::shared_ptr<Systematic> > nonGammaSystematics;
    const std::vector<std::string>& uniqueSyst = GetUniqueSystNamesWithoutGamma();
    for(const auto& isyst : fSystematics) {
        if(isyst->fType == Systematic::SHAPE) continue;

        nonGammaSystematics.push_back(isyst);
    }
    const size_t NnonGammaSyst = nonGammaSystematics.size();
    if(NnonGammaSyst==0){
        LOG(INFO) << "No non-gamma systematics found => No Pruning plot generated.\n";
        return;
    }
    //
    const std::string suffix = fRemoveLargeSyst ? " (dropped)" : "";
    int iReg = 0;
    for(const auto& ireg : fRegions) {
        if(!fValidationPruning && ireg->fRegionType==Region::VALIDATION) continue;

        const std::size_t size = samplesVec.at(iReg).size();
        out << "In Region : " << ireg->fName << std::endl ;
        histPrun.emplace_back(std::make_unique<TH2F>(Form("h_prun_%s", ireg->fName.c_str()  ),ireg->fShortLabel.c_str(),size,0,size, uniqueSyst.size(),0,uniqueSyst.size()));
        histPrun.back()->SetDirectory(nullptr);

        for(std::size_t i_smp=0;i_smp<size;i_smp++){
            out << " -> In Sample : " << samplesVec.at(iReg).at(i_smp)->fName << std::endl;

            for(std::size_t uniqueIndex = 0; uniqueIndex < uniqueSyst.size(); ++uniqueIndex){
               histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -1 );
            }

            std::shared_ptr<SampleHist> sh = ireg->GetSampleHist(samplesVec.at(iReg).at(i_smp)->fName);
            if (sh == nullptr) continue;

            for(size_t i_syst=0;i_syst<NnonGammaSyst;i_syst++){
                // find the corresponding index of unique syst
                auto it = std::find(uniqueSyst.begin(), uniqueSyst.end(), nonGammaSystematics.at(i_syst)->fName);
                const std::size_t uniqueIndex = std::distance(uniqueSyst.begin(), it);
                out << " --->>  " << nonGammaSystematics[i_syst]->fName << "     " ;
                if( (Common::FindInStringVector(nonGammaSystematics[i_syst]->fSamples,samplesVec.at(iReg).at(i_smp)->fName)>=0 || nonGammaSystematics[i_syst]->fSamples[0] == "all")
                    && sh->HasSyst(nonGammaSystematics[i_syst]->fName)
                ){
                    std::shared_ptr<SystematicHist> syh = sh->GetSystematic(nonGammaSystematics[i_syst]->fName);
                    histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 0 );
                    const bool forgeDropShape = nonGammaSystematics[i_syst]->fIsNormOnly;
                    const bool forgeDropNorm  = nonGammaSystematics[i_syst]->fIsShapeOnly;
                    //
                    if(syh->fShapePruned && syh->fNormPruned) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 3 );
                    else if(syh->fShapePruned || forgeDropShape) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 1 );
                    else if(syh->fNormPruned || forgeDropNorm) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), 2 );
                    //
                    if(syh->fBadNorm) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -2 );
                    if(syh->fBadShape) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -3 );
                    if(syh->fBadShape && syh->fBadNorm) histPrun[iReg]->SetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex), -4 );
                    //
                }
                if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -1 ) out << " is not present" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 0 ) out << " is kept" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 1 ) out << " is norm only" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 2 ) out << " is shape only" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== 3 ) out << " is dropped" << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -2 ) out << " has bad norm" << suffix << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -3 ) out << " has bad shape" << suffix << std::endl;
                else if( histPrun[iReg]->GetBinContent( histPrun[iReg]->FindBin(i_smp,uniqueIndex) )== -4 ) out << " is bad" << suffix << std::endl;
            }
        }
        //
        histPrun_toSave.emplace_back(std::unique_ptr<TH2F>(static_cast<TH2F*>(histPrun[iReg]->Clone(Form("%s_toSave",histPrun[iReg]->GetName())))));
        histPrun_toSave[iReg]->SetDirectory(0);
        //
        iReg++;
    }
    if (fDrawPruningPlot) {
        int totalSize(0);
        for (const auto& ihist : histPrun) {
            totalSize += ihist->GetNbinsX();
        }
        if (totalSize > 1000) {
            LOG(WARNING) << "Too many (>1000) samples x regions for the pruning plot. This might cause running time and memory issues\n";
            LOG(WARNING) << "Consider skipping the plotting part for the pruning plot with \"DrawPruningPlot: FALSE\"\n";
        }
        //
        // draw the histograms
        int upSize = 50;
        int loSize = 150;
        int mainHeight = uniqueSyst.size()*20;
        int leftSize = 250;
        int regionSize = 20*totalSize/iReg;
        int separation = 10;
        int mainWidth = iReg*(regionSize+2*separation);
        //
        TCanvas c("c_pruning","Canvas - Pruning",leftSize+mainWidth,upSize+mainHeight+loSize);
        std::vector<Int_t> colors = {kBlack,6,kBlue, kGray, 8, kYellow, kOrange-3, kRed}; // #colors >= #levels - 1
        gStyle->SetPalette(colors.size(), &colors[0]);
        TPad pUp("pUp","Pad High",0,(1.*loSize+mainHeight)/(upSize+mainHeight+loSize),1,1);
        pUp.Draw();
        c.cd();
        std::vector<std::unique_ptr<TPad> > pReg(100);
        for(std::size_t i_reg=0;i_reg<histPrun.size();i_reg++){
            c.cd();
            if(i_reg==0){
                pReg[i_reg] = std::make_unique<TPad>(Form("pReg[%zu]",i_reg),"Pad Region",
                                      0,   0,
                                      (leftSize+1.*i_reg*(regionSize+separation)+regionSize)/(leftSize+mainWidth),   (1.*loSize+mainHeight)/(upSize+mainHeight+loSize));
                pReg[i_reg]->SetLeftMargin( (1.*leftSize) / (1.*leftSize+regionSize));
            }
            else{
                pReg[i_reg] = std::make_unique<TPad>(Form("pReg[%zu]",i_reg),"Pad Region",
                                      (leftSize+1.*i_reg*(regionSize+separation))           /(leftSize+mainWidth),   0,
                                      (leftSize+1.*i_reg*(regionSize+separation)+regionSize)/(leftSize+mainWidth),   (1.*loSize+mainHeight)/(upSize+mainHeight+loSize));
                pReg[i_reg]->SetLeftMargin(0);
            }
            pReg[i_reg]->SetBottomMargin( (1.*loSize) / (1.*loSize+mainHeight) );
            pReg[i_reg]->Draw();
            pReg[i_reg]->cd();
            gPad->SetGridy();
            for(int i_bin=1;i_bin<=histPrun[i_reg]->GetNbinsX();i_bin++){
                histPrun[i_reg]       ->GetXaxis()->SetBinLabel(i_bin,samplesVec.at(i_reg).at(i_bin-1)->fTitle.c_str());
                histPrun_toSave[i_reg]->GetXaxis()->SetBinLabel(i_bin,samplesVec.at(i_reg).at(i_bin-1)->fName.c_str());
            }
            for(int i_bin=1;i_bin<=histPrun[i_reg]->GetNbinsY();i_bin++){
                if(i_reg==0) {
                    histPrun[i_reg]->GetYaxis()->SetBinLabel(i_bin,TRExFitter::SYSTMAP[uniqueSyst[i_bin-1]].c_str());
                }
                else {
                    histPrun[i_reg]->GetYaxis()->SetBinLabel(i_bin,"");
                }
                histPrun_toSave[i_reg]->GetYaxis()->SetBinLabel(i_bin,uniqueSyst[i_bin-1].c_str());
            }
            histPrun[i_reg]->Draw("COL");
            histPrun[i_reg]->GetYaxis()->SetLabelOffset(0.03);
            gPad->SetTopMargin(0);
            gPad->SetRightMargin(0);
            histPrun[i_reg]->GetXaxis()->LabelsOption("v");
            histPrun[i_reg]->GetXaxis()->SetLabelSize( histPrun[i_reg]->GetXaxis()->GetLabelSize()*0.75 );
            histPrun[i_reg]->GetYaxis()->SetLabelSize( histPrun[i_reg]->GetYaxis()->GetLabelSize()*0.75 );
            gPad->SetTickx(0);
            gPad->SetTicky(0);
            histPrun[i_reg]->SetMinimum(-4);
            histPrun[i_reg]->SetMaximum( 3.1);
            histPrun[i_reg]->GetYaxis()->SetTickLength(0);
            histPrun[i_reg]->GetXaxis()->SetTickLength(0);
            gPad->SetGrid();
            //
            pUp.cd();
            myText((leftSize+1.*i_reg*(regionSize+separation))/(leftSize+mainWidth),0.1 ,1,histPrun[i_reg]->GetTitle());
        }
        c.cd();
        TPad pLo("pLo","Pad Low",0,0,(1.*leftSize)/(leftSize+mainWidth),(1.*loSize)/(upSize+mainHeight+loSize));
        pLo.Draw();
        //
        c.cd();
        pUp.cd();
        myText(0.01,0.5,1,fLabel.c_str());
        //
        pLo.cd();
        TLegend leg(0.005,0,0.95,0.95);
        TH1D hGray   ("hGray"  ,"hGray"  ,1,0,1);    hGray.SetFillColor(kGray);         hGray.SetLineWidth(0);
        TH1D hYellow ("hYellow","hYellow",1,0,1);    hYellow.SetFillColor(kYellow);     hYellow.SetLineWidth(0);
        TH1D hOrange ("hOrange","hOrange",1,0,1);    hOrange.SetFillColor(kOrange-3);   hOrange.SetLineWidth(0);
        TH1D hRed    ("hRed"   ,"hRed"   ,1,0,1);    hRed.SetFillColor(kRed);           hRed.SetLineWidth(0);
        TH1D hGreen  ("hGreen" ,"hGree"  ,1,0,1);    hGreen.SetFillColor(8);            hGreen.SetLineWidth(0);
        TH1D hBlue   ("hBlue"  ,"hBlue"  ,1,0,1);    hBlue.SetFillColor(kBlue);         hBlue.SetLineWidth(0);
        TH1D hPurple ("hPurple","hPurple",1,0,1);    hPurple.SetFillColor(6);           hPurple.SetLineWidth(0);
        TH1D hBlack  ("hBlack" ,"hBlack" ,1,0,1);    hBlack.SetFillColor(kBlack);       hBlack.SetLineWidth(0);
        std::string sysLarg="Norm >" + std::to_string((int)(fThresholdSystLarge*100))+"%"+suffix;
        leg.SetBorderSize(0);
        leg.SetMargin(0.1);
        leg.SetFillStyle(0);
        leg.AddEntry(&hGray,"Not present","f");
        leg.AddEntry(&hGreen,"Kept","f");
        leg.AddEntry(&hYellow, "Shape dropped","f");
        leg.AddEntry(&hOrange, "Norm. dropped","f");
        leg.AddEntry(&hRed, "Dropped","f");
        if (fThresholdSystLarge > -1) {
            leg.AddEntry(&hBlue  , sysLarg.c_str() ,"f");
            leg.AddEntry(&hPurple, ("Bad shape" + suffix).c_str() ,"f");
            leg.AddEntry(&hBlack , ("Bad shape & norm." + suffix).c_str() ,"f");
        }
        #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
        leg.SetTextSize(0.85*gStyle->GetTextSize());
        #endif
        leg.Draw();
        //
        Common::SaveCanvasAs(c, fName+"/Pruning"+fSuffix);
    }

    //
    // Save prunign hist for future usage
    std::unique_ptr<TFile> filePrun = nullptr;
    // - checking if Pruning.root exists
    // if yes
    if(!gSystem->AccessPathName( (fName+"/Pruning.root").c_str() )){
        // ...
        filePrun.reset(TFile::Open( (fName+"/Pruning.root").c_str() ));
    }
    else{
        filePrun.reset(TFile::Open( (fName+"/Pruning.root").c_str(),"RECREATE" ));
        for(std::size_t i_reg=0;i_reg<histPrun.size();i_reg++){
            histPrun_toSave[i_reg]->Write("",TObject::kOverwrite);
        }
    }
    if (filePrun != nullptr){
        filePrun->Close();
    }
}

//__________________________________________________________________________________
//
void TRExFit::Fit(bool isLHscanOnly){

    std::unique_ptr<RooDataSet> data(nullptr);
    std::unique_ptr<RooWorkspace> ws(nullptr);

    //
    // Read NPvalues from fit-result file
    //
    if (fFitNPValuesFromFitResultsFile!=""){
        LOG(INFO) << "Setting NP values for Asimov data-set creation from fit results stored in file " << fFitNPValuesFromFitResultsFile << "...\n";
        fFitNPValues = NPValuesFromFitResultsFile(fFitNPValuesFromFitResultsFile);
    }

    //
    // If fDoNonProfileFit => set stat-only
    //
    if (fDoNonProfileFit){
        LOG(INFO) << "In non-profile mode => Setting to stat-only.\n";
        fStatOnly = true;
    }

    std::pair<std::unique_ptr<RooWorkspace>, std::unique_ptr<RooDataSet> > wsAndData;
    //
    // If there's a workspace specified, go on with simple fit, without looking for separate workspaces per region
    //
    std::unique_ptr<TFile> rootFile(nullptr);

    // Create cache and use cache must be mutually exclusive operations
    if (fCreateCache && fUseCache) fUseCache = false;

    // If user requests to use a cache, make sure it exists then read the workspace and data
    // If cache does not exist, abort
    if (fUseCache) {
        std::string cacheFileName = "FitWorkspace_POI";
        for (const auto& val : fFitPOIAsimov) {
            cacheFileName += "_" + val.first + "_" + std::to_string(val.second);
        }
        const std::string cacheFilePath = fName+"/RooStats/" + cacheFileName + ".root";

        if (gSystem->AccessPathName(cacheFilePath.c_str())) {
            LOG(ERROR) << "Request to use cache, but cache does not exist! Aborting!\n";
            return;
        }

        rootFile.reset(TFile::Open(cacheFilePath.c_str(), "read"));
        ws = std::unique_ptr<RooWorkspace>(dynamic_cast<RooWorkspace*>(rootFile->Get("combined")));
        if(!ws){
            LOG(ERROR) << "The workspace (\"combined\") cannot be found in file " << cacheFilePath << ". Please check!\n";
            return;
        }

        data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data("newasimovData")));
        if(!data){
            LOG(ERROR) << "Cannot read the data from custom WS. Please check!\n";
            return;
        }
        LOG(INFO) << "Request to use cache, data and workspace successfully read.\n";
    }
    else if(fWorkspaceFileName!=""){
        LOG(INFO) << "\n";
        LOG(INFO) << "-------------------------------------------\n";
        LOG(INFO) << "Performing nominal fit on pre-specified workspace...\n";
        rootFile.reset(TFile::Open(fWorkspaceFileName.c_str(), "read"));
        if (!rootFile) {
            LOG(ERROR) << "Cannot open custom WS ROOT file at: " << fWorkspaceFileName << "\n";
            exit(EXIT_FAILURE);
        }
        ws = std::unique_ptr<RooWorkspace>(dynamic_cast<RooWorkspace*>(rootFile->Get("combined")));
        if(!ws){
            LOG(ERROR) << "The workspace (\"combined\") cannot be found in file " << fWorkspaceFileName << ". Please check!\n";
            exit(EXIT_FAILURE);
        }
        if(!fFitIsBlind){
            data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data("obsData")));
        }
        else{
            if(fCustomAsimov!=""){
                //not sure if this works
                std::string custom_asimov_name = "customAsimov_"+fCustomAsimov;
                data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data(custom_asimov_name.c_str())));
            }
            else{
                    data = std::unique_ptr<RooDataSet>(static_cast<RooDataSet*>(ws->data("asimovData")));
            }
        }

        if(!data){
            LOG(ERROR) << "Cannot read the data from custom WS. Please check!\n";
            exit(EXIT_FAILURE);
        }
    }
    //
    // Otherwise go on with normal fit
    //
    else{
        wsAndData = this->PrepareMixedDataset(TRExFit::WorkspaceType::FIT);
        ws = std::move(wsAndData.first);
        data = std::move(wsAndData.second);

        if (DoingMixedFitting()) {
            std::unique_ptr<TFile> rootFileCombined(nullptr);

            if(fBootstrap!="" && fBootstrapIdx>=0) {
                rootFileCombined.reset(TFile::Open( (fName+"/RooStats/"+fBootstrapSyst+fBootstrapSample+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read"));
            } else {
                rootFileCombined.reset(TFile::Open( (fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read"));
            }

            if (!rootFileCombined) {
                LOG(ERROR) << "Cannot read the WS\n";
                exit(EXIT_FAILURE);
            }

            auto measurement = std::unique_ptr<RooStats::HistFactory::Measurement>(static_cast<RooStats::HistFactory::Measurement*>(rootFileCombined -> Get( (fInputName+fSuffix).c_str())));
            // store data in the WS
            const auto name = data->GetName();
            data->SetName("mixedAsimovData");
            ws->import(*(data.get()));
            ws->writeToFile((fName+"/RooStats/"+fInputName+"_realisticAsimov.root").c_str());
            data->SetName(name);

            std::unique_ptr<TFile> out(TFile::Open((fName+"/RooStats/"+fInputName+"_realisticAsimov.root").c_str(), "UPdate"));
            if (!out) {
                LOG(ERROR) << "Cannot open the out file to write the measurement object\n";
                exit(EXIT_FAILURE);
            }
            measurement->writeToFile(out.get());

            out->Close();
            rootFileCombined->Close();
        }
        if(fWriteCustomAsimovToWS) {
            // store data in the WS
            const auto name = data->GetName();
            data->SetName("customAsimov");
            ws->import(*(data.get()));
            ws->writeToFile((fName+"/RooStats/"+fInputName+"_customAsimov.root").c_str());
            data->SetName(name);
            LOG(INFO) << "Wrote out RooStats workspace containing custom Asimov.\n";
        }

        // "data" is a new dataset not saved in the workspace
        // import it then save to a root file for later reuse
        if (fCreateCache) {
            std::string cacheFileName = "FitWorkspace_POI";
            for (const auto& val : fFitPOIAsimov) {
                cacheFileName += "_" + val.first + "_" + std::to_string(val.second);
            }
            ws->import(*(data.get()));
            ws->writeToFile((fName+"/RooStats/" + cacheFileName + ".root").c_str(), true);

            // When the scan step is rerun, it will use the saved workspace. Terminate now
            if (isLHscanOnly && fLHscanStep != -1) {
                LOG(INFO) << "You are parallel-processing the Likelihood scan.\n";
                LOG(INFO) << "The workspace has been saved and will be reloaded\n";
                LOG(INFO) << "Please rerun this scan step\n";
                return;
            }
        }
    }


    if (!isLHscanOnly){
        LOG(INFO) << "\n";
        LOG(INFO) << "-------------------------------------------\n";
        LOG(INFO) << "Performing nominal fit...\n";

        //
        // Calls the PerformFit() function to actually do the fit
        //
        PerformFit(ws.get(), data.get(), fFitType, true);
    }

    //
    // Toys
    //
    if(fFitToys>0 && !isLHscanOnly){
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
        RunToys();
    }

    //
    // Fit result on Asimov with shifted systematics
    //
    if(fDoNonProfileFit && !isLHscanOnly){
        this->NonProfiledFit();
    }

    //
    // Calls the  function to create LH scan with respect to a parameter
    //
    // get list of all parameters
    const RooStats::ModelConfig* mc = dynamic_cast<const RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Passed nullptr for MC\n";
        exit(EXIT_FAILURE);
    }
    const std::vector<std::string> parameters = FitUtils::GetAllParameters(mc);

    bool running2Dparallel = ( (fLHscanStepY != -1 || fLHscanStep != -1) && fVarName2DLH.size() != 0);

    if(fVarNameLH.size()>0 && !isLHscanOnly && !running2Dparallel){
        //
        // Don't do it if you did a non-profile fit (FIXME)
        if(fDoNonProfileFit){
            LOG(WARNING) << "Better not to perform LH scan if you did non-profile fit with scan on systematics. Skipping LH scan.\n";
        }
        else{
            if (fVarNameLH[0]=="all") {
                for(const auto& iparam : parameters) {
                    GetLikelihoodScan( ws.get(), iparam, data.get());
                }
            }
            else{
                for(const auto& iparam : fVarNameLH) {
                    GetLikelihoodScan( ws.get(), iparam, data.get());
                }
            }
            if (fVarName2DLH.size() > 0) {
                for (const auto & ipair : fVarName2DLH) {
                    Get2DLikelihoodScan( ws.get(), ipair, data.get());
                }
            }
        }
    }
    if (isLHscanOnly && !running2Dparallel){
        if (fVarNameLH.size() == 0 && fVarName2DLH.size() == 0){
            LOG(ERROR) << "Did not provide any LH scan parameter and running LH scan only. This is not correct.\n";
            exit(EXIT_FAILURE);
        }
        if (fVarNameLH.size() > 0) {
            if (fVarNameLH[0]=="all"){
                LOG(WARNING) << "You are running LHscan only option but running it for all parameters. Will not parallelize!\n";
                for(const auto& iparam : parameters) {
                    GetLikelihoodScan( ws.get(), iparam, data.get());
                }
            } else {
                GetLikelihoodScan( ws.get(), fVarNameLH[0], data.get());
            }
        }
        if (fVarName2DLH.size() > 0) {
            for (const auto & ipair : fVarName2DLH) {
                Get2DLikelihoodScan( ws.get(), ipair, data.get());
            }
        }
    }

    // If we read a workspace from file, close the TFile to avoid memory leak.
    // Workspace "ws" owns the "data" object in this case. To avoid a double delete
    // of "data", release data without deleting it. It will be deleted when ws deletes.
    if (rootFile){
        // cppcheck-suppress ignoredReturnValue
        data.release();
        rootFile->Close();
    }

    // we need to relaase the memory i nthis case otherwise root would result in a double delete
    if(fWorkspaceFileName!=""){
        // cppcheck-suppress ignoredReturnValue
        data.release();
        // cppcheck-suppress ignoredReturnValue
        ws.release();
    }
}


//__________________________________________________________________________________
//
RooDataSet* TRExFit::DumpData( RooWorkspace *ws, const std::map <std::string, int>& regionDataType, const std::map <std::string, double>& npValues, const std::map <std::string, double>& poiValues ){
    //
    // This function dumps a RooDataSet object using the input informations provided by the user
    //    |-> Used when testing Fit response (inject one NP in data and check fit result)
    //    |-> Used when using fit results in some regions to generate Asimov data in blinded regions
    //
    LOG(DEBUG) << "Dumping data with the following parameters\n";
    LOG(DEBUG) << "    * Regions data type\n";
    for(const auto& dataType : regionDataType ){
        LOG(DEBUG) << "       - Region: " << dataType.first << "       DataType: " << dataType.second << "\n";
    }
    if(npValues.size()){
        LOG(DEBUG) << "    * Injected NP values \n";
        for (const auto& npValue : npValues ){
            LOG(DEBUG) << "       - NP: " << npValue.first << "       Value: " << npValue.second << "\n";
        }
    }
    else {
        LOG(DEBUG) << "    * No NP values injected\n";
    }
    if(poiValues.size()){
        LOG(DEBUG) << "    * Injected POI values\n";
        for (const auto& poiValue : poiValues ){
            LOG(DEBUG) << "       - POI: " << poiValue.first << "       Value: " << poiValue.second << "\n";
        }
    }
    else {
        LOG(DEBUG) << "    * No POI values injected\n";
    }

    RooStats::ModelConfig *mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));

    //Save the initial values of the NP
    if (!ws->getSnapshot("InitialStateModelGlob")) {
        ws->saveSnapshot("InitialStateModelGlob",   *mc->GetGlobalObservables());
    }
    if (!(fStatOnly && fFitIsBlind)){
        if (mc->GetNuisanceParameters()) ws->saveSnapshot("InitialStateModelNuis",   *mc->GetNuisanceParameters());
    }

    //Be sure to take the initial values of the NP
    ws->loadSnapshot("InitialStateModelGlob");
    if (!fStatOnly){
        ws->loadSnapshot("InitialStateModelNuis");
    }

    if (fBinnedLikelihood) {
        FitUtils::SetBinnedLikelihoodOptimisation(ws);
    }

    if (fIntCode != 4) {
        FitUtils::ChangeInterpolationCode(ws, fIntCode);
    }

    //Creating a set
    static const std::string weightName="weightVar";
    const auto* observables = mc->GetObservables();
    RooArgSet obsAndWeight;
    obsAndWeight.add(*observables);

    RooRealVar* weightVar = nullptr;
    if ( !(weightVar = ws->var(weightName.c_str())) ){
        RooRealVar tmp(weightName.c_str(), weightName.c_str(), 1);
        ws->import(tmp);
        weightVar = ws->var(weightName.c_str());
    }
    obsAndWeight.add(*ws->var(weightName.c_str()));
    ws->defineSet("obsAndWeight",obsAndWeight);

    //
    // Getting observed data (in case some regions are unblinded)
    //
    RooDataSet* realData = static_cast<RooDataSet*>(ws->data("obsData"));

    //
    // Set some parameters for the Asimov production
    //     |-> Values of NPs
    //     |-> Values of POI
    //

    //-- POIs
    for (auto poi_tmp : *mc->GetParametersOfInterest()) {
        RooRealVar* poi = static_cast<RooRealVar*>(poi_tmp);
        auto it_poiValue = poiValues.find( poi -> GetName() );
        if( it_poiValue != poiValues.end() ){
            poi -> setVal(it_poiValue -> second);
        }
    }

    //-- Nuisance parameters
    if (mc->GetNuisanceParameters()) {
        for (auto var_tmp : *mc->GetNuisanceParameters()) {
            RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
            auto it_npValue = npValues.find( var -> GetName() );
            if( it_npValue != npValues.end() ){
                var -> setVal(it_npValue -> second);
            }
        }
    }

    //Looping over regions
    std::stack<std::unique_ptr<RooAbsData> > ownedAsimovData;
    std::map<std::string, RooDataSet*> asimovDataMap;
    RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
    RooCategory* channelCat = const_cast<RooCategory*>(static_cast<const RooCategory*>(&simPdf->indexCat()));
    for (std::size_t icat = 0; icat < channelCat->size(); ++icat) {
        channelCat->setIndex(icat);
        //Check the type of data to store for this region !
        int dataType = Region::ASIMOVDATA;//default is AsimovData
        auto it_dataType = regionDataType.find(channelCat->getLabel());
        if(it_dataType == regionDataType.end()) {
            const std::string temp_string = channelCat->getLabel();
            LOG(WARNING) << "The following region is not specified in the inputs to the function (" << temp_string << "): use Asimov\n";
            LOG(WARNING) << "   This can happen for CRONLY fits. Please check if everything is fine!\n";
            continue;
        } else {
            dataType = it_dataType->second;
        }

        //A protection: if there is no real observed data, use only ASIMOV (but print a warning)
        if(dataType==Region::REALDATA && !realData){
            std::string temp_string = channelCat->getLabel();
            LOG(WARNING) << "You want real data for channel " << temp_string << " but none is available in the workspace. Using Asimov instead.\n";
            dataType = Region::ASIMOVDATA;
        }

        if(dataType==Region::ASIMOVDATA){
            if(fCustomAsimov!=""){
                ownedAsimovData.emplace(realData->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),channelCat->getLabel())));
                asimovDataMap[channelCat->getLabel()] = static_cast<RooDataSet*>(ownedAsimovData.top().get());
            }
            else{
                // Get pdf associated with state from simpdf
                RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel());

                // Generate observables defined by the pdf associated with this state
                std::unique_ptr<RooArgSet> obstmp(pdftmp->getObservables(*observables));

                ownedAsimovData.emplace(std::make_unique<RooDataSet>(Form("combAsimovData%ld",icat),Form("combAsimovData%ld",icat),RooArgSet(obsAndWeight,*channelCat),RooFit::WeightVar(*weightVar)));
                RooRealVar* thisObs = static_cast<RooRealVar*>(obstmp->first());
                const double expectedEvents = pdftmp->expectedEvents(*obstmp);

                for(int jj=0; jj<thisObs->numBins(); ++jj){
                    thisObs->setBin(jj);
                    const double thisNorm=pdftmp->getVal(obstmp.get())*thisObs->getBinWidth(jj);
                    if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < 1e18) {
                        ownedAsimovData.top()->add(*observables, thisNorm*expectedEvents);
                    }
                }
                if (TRExFitter::DEBUGLEVEL >= 2) {
                    ownedAsimovData.top()->Print();
                }
                if (std::isnan(ownedAsimovData.top()->sumEntries())) {
                    LOG(ERROR) << "Sum entries is NaN\n";
                    exit(EXIT_FAILURE);
                }
                asimovDataMap[channelCat->getLabel()] = static_cast<RooDataSet*>(ownedAsimovData.top().get());
            }

        } else if(dataType==Region::REALDATA) {
            ownedAsimovData.emplace(realData->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),channelCat->getLabel())));
            asimovDataMap[channelCat->getLabel()] = static_cast<RooDataSet*>(ownedAsimovData.top().get());
        }
    }

    RooDataSet *asimovData = new RooDataSet("newasimovData",
                                            "newasimovData",
                                            RooArgSet(obsAndWeight,*channelCat),
                                            Index(*channelCat),
                                            Import(asimovDataMap),
                                            WeightVar(*weightVar));

    ws->loadSnapshot("InitialStateModelGlob");
    if (!fStatOnly){
        ws->loadSnapshot("InitialStateModelNuis");
    }

    return asimovData;
}

//__________________________________________________________________________________
//
std::map < std::string, double > TRExFit::PerformFit( RooWorkspace *ws, RooDataSet* inputData, FitType fitType, bool save){

    std::map < std::string, double > result;

    /////////////////////////////////
    //
    // Function performing a fit in a given configuration.
    //
    /////////////////////////////////

    std::vector<std::pair<std::string,double> > smallNPs;
    if (fSpeedUpFit) {
        // get the list of ordered NPs
        std::ifstream file;
        file.open(fName+"/NPorder.txt");
        if (!file.is_open() || !file.good()) {
            LOG(WARNING) << "Cannot open the txt file needed to identify small NPs! Will ignore it\n";
        } else {
            std::string syst;
            double value;
            while (file >> syst >> value) {
                smallNPs.emplace_back(syst,value);
            }
        }
        smallNPs.resize(static_cast<std::size_t>(fNPCutOff * smallNPs.size()));
    }

    // prepare vectors for starting point of normfactors
    std::vector<std::string> NPnames;
    std::vector<double> NPvalues;
    for(const auto& inf : fNormFactors) {
        if (Common::FindInStringVector(fPOIs,inf->fName)>=0) {
            continue;
        }
        NPnames. emplace_back(inf->fName);
        NPvalues.emplace_back(inf->GetNominal());
    }
    //
    // Fit configuration (SPLUSB or BONLY)
    //
    FittingTool fitTool{};
    fitTool.SetUseHesse(fUseHesse);
    fitTool.SetUseHesseBeforeMigrad(fUseHesseBeforeMigrad);
    fitTool.SetStrategy(fFitStrategy);
    fitTool.SetNCPU(fCPU);
    fitTool.SetSpeedUpFit(fSpeedUpFit);
    fitTool.SetSmallNPs(smallNPs);
    fitTool.SetErrorSigma(fErrorSigma);
    fitTool.SetMaxFCNcalls(fMaximumNumberFCNcalls);
    fitTool.SetToleranceScale(fToleranceScale);
    fitTool.SetUseAutoDiff(fUseAutoDiff);
    fitTool.SetNLLOffset(fNLLOffset);
    if(fitType == BONLY){
        for (const auto& inf : fNormFactors) {
            fitTool.AddValPOI(inf->fName, 0.);
        }
        fitTool.ConstPOI(true);
    } else if(fitType == SPLUSB || fitType == UNFOLDING || fitType == EFT){
        for (const auto& inf : fNormFactors) {
            fitTool.AddValPOI(inf->fName, inf->GetNominal());
        }
        fitTool.ConstPOI(false);
    }
    fitTool.SetNPs( NPnames,NPvalues );
    fitTool.SetRandomNP(fRndRange, fUseRnd, fRndSeed);
    if(fStatOnly){
        if(!fGammasInStatOnly) fitTool.NoGammas();
        fitTool.NoSystematics();
    }
    std::vector<std::string> constNFs;
    for (const auto& inf : fNormFactors) {
        if (inf->fConst) {
            constNFs.emplace_back(inf->fName);
        }
    }
    fitTool.SetConstantNFs(constNFs);

    //
    // Fit starting from custom point
    if(fFitResultsRootFile!=""){
        ReadFitResults(fFitResultsRootFile);
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp : fFitResults->GetNuisanceParameters()) {
            npNames.push_back( inp.second->fName);
            npValues.push_back(inp.second->fFitValue);
        }
        fitTool.SetNPs(npNames,npValues);
    }

    //
    // Set Minos
    if(fVarNameMinos.size()>0){
        LOG(DEBUG) << "Setting the variables to use MINOS with\n";
        fitTool.UseMinos(fVarNameMinos);
    }

    //
    // Gets needed objects for the fit
    //
    RooStats::ModelConfig* mc = static_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    RooSimultaneous *simPdf = static_cast<RooSimultaneous*>(mc->GetPdf());

    //
    // Creates the data object
    //
    RooDataSet* data = nullptr;
    bool usingObsOrAsimovData(false);
    if(inputData){
        data = inputData;
    } else {
        LOG(WARNING) << "You didn't provide inputData => will use the observed data!\n";
        data = static_cast<RooDataSet*>(ws->data("obsData"));
        if(data==nullptr){
            LOG(WARNING) << "No observedData found => will use the Asimov data!\n";
            data = static_cast<RooDataSet*>(ws->data("asimovData"));
        }
        inputData = data;
        usingObsOrAsimovData = true;
    }
    if (!usingObsOrAsimovData) {
        ws->import(*data);
    }

    //
    // For stat-only fit on data:
    // - read fit resutls
    // - fix all NP to fitted ones before fitting
    if(fStatOnlyFit){
        LOG(DEBUG) << "Fitting stat-only: reading fit results from full fit from file:\n";
        LOG(DEBUG) << "  " << fName+"/Fits/"+fInputName+fSuffix+".root\n";
        const bool exist = ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
        if (!exist) {
            LOG(ERROR) << "No fit results file available, please run the full fit before running the StatOnlyFit fit\n";
            exit(EXIT_FAILURE);
        }
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& inp : fFitResults->GetNuisanceParameters()) {
            if(!fFixNPforStatOnlyFit && Common::FindInStringVector(fNormFactorNames,inp.second->fName)>=0) continue;
            auto itr = std::find_if(fShapeFactorNames.begin(), fShapeFactorNames.end(), [&inp](const auto& element){return inp.second->fName.find(element+"_bin") != std::string::npos;});
            if(!fFixNPforStatOnlyFit && itr != fShapeFactorNames.end()) continue;
            npNames.push_back( inp.second->fName);
            npValues.push_back(inp.second->fFitValue);
        }
        fitTool.FixNPs(npNames,npValues);
    }

    // FixNP
    if(fFitFixedNPs.size()>0){
        std::vector<std::string> npNames;
        std::vector<double> npValues;
        for(const auto& nuisParToFix : fFitFixedNPs){
            npNames.push_back( nuisParToFix.first );
            npValues.push_back( nuisParToFix.second );
        }
        fitTool.FixNPs(npNames,npValues);
    }

    FitUtils::ApplyExternalConstraints(ws, &fitTool, simPdf, fNormFactors, fRegularizationType);

    // save snapshot before fit
    if (!ws->getSnapshot("snapshot_BeforeFit_POI")) {
        ws->saveSnapshot("snapshot_BeforeFit_POI", *(mc->GetParametersOfInterest()) );
        if (mc->GetNuisanceParameters()) ws->saveSnapshot("snapshot_BeforeFit_NP" , *(mc->GetNuisanceParameters())   );
        ws->saveSnapshot("snapshot_BeforeFit_GO" , *(mc->GetGlobalObservables())    );
    }

    //
    // Get initial ikelihood value from Asimov
    if (fBlindedParameters.size() > 0) std::cout.setstate(std::ios_base::failbit);

    // save snapshot after fit
    ws->saveSnapshot("snapshot_AfterFit_POI", *(mc->GetParametersOfInterest()) );
    if (mc->GetNuisanceParameters()) ws->saveSnapshot("snapshot_AfterFit_NP" , *(mc->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(mc->GetGlobalObservables())    );

    // Performs the fit
    fitTool.FitPDF( mc, simPdf, data );
    auto fitResults = fitTool.GetFitResult();
    if (fBlindedParameters.size() == 0) std::cout.clear();

    // Error decomposition
    if (fErrorDecomposition && !fDoGroupedSystImpactTable && !fStatOnlyFit && (fFitType != FitType::BONLY) && (fErrorSigma < 2)) {
        LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
        LOG(INFO) << "-----------------------     ERROR DECOMPOSITION    -----------------------\n";
        LOG(INFO) << "  POI :  +/-tot.err  (+/-stat  +/-syst)\n";
        LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
        for (const auto& poi : fPOIs) {
            std::map<std::string,double> impactMap = fitTool.CalculateErrorDecomposition(mc,poi);
            if (impactMap.size() > 0) {
                LOG(INFO) << Form("  %s :  +/-%g  (+/-%g  +/-%g)",poi.c_str(),impactMap["TOT_ERROR"],impactMap["STAT_ERROR"],impactMap["TOTSYST_ERROR"]) << "\n";
            }
        }
        if(save) LOG(INFO) << "  Full breakdown will be written to file in the Fits/ folder.\n";
        LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
    }

    // Goodness of fit
    if(!fDoGroupedSystImpactTable && fitResults){
        xRooNode rooNode(*ws);
        rooNode.SetFitResult(*fitResults);
        auto pdfNode = rooNode["simPdf"];
        if (!pdfNode) {
            LOG(ERROR) << "Cannot read simPdf from the WS\n";
            exit(EXIT_FAILURE);
        }
        const double prob = FitUtils::SaturatedModelGoF(pdfNode, data->GetName());
        LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
        LOG(INFO) << "----------------------- GOODNESS OF FIT EVALUATION -----------------------\n";
        LOG(INFO) << "  probability = " << prob << "\n";
        LOG(INFO) << "----------------------- -------------------------- -----------------------\n";
    }
    if(save){
        ProduceSystSubCategoryMap();
        fitTool.SetSystMap( fSubCategoryImpactMap );
        if(fBootstrap!="" && fBootstrapIdx>=0){
            gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
            if(fStatOnlyFit) fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+"_statOnly.txt", fBlindedParameters);
            else             fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".txt", fBlindedParameters);
            if(fErrorDecomposition && !fStatOnlyFit && !fDoGroupedSystImpactTable && (fFitType != FitType::BONLY) && (fErrorSigma < 2)){
                for (const auto& poi : fPOIs) {
                    fitTool.ExportErrorDecompositionInTextFile(poi,fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+"_errDecomp_"+poi+".txt");
                    fitTool.ExportErrorDecompositionGroupInTextFile(poi,fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix, fLumiLabel, fCmeLabel, fName, fHEPDataFormat);
                }
                this->PlotRankingFromCovMatrix(fitTool);
            }
        }
        else{
            gSystem -> mkdir((fName+"/Fits/").c_str(),true);
            if(fStatOnlyFit) fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+"_statOnly.txt", fBlindedParameters);
            else             fitTool.ExportFitResultInTextFile(fName+"/Fits/"+fInputName+fSuffix+".txt", fBlindedParameters);
            if(fErrorDecomposition && !fStatOnlyFit && !fDoGroupedSystImpactTable && (fFitType != FitType::BONLY) && (fErrorSigma < 2)){
                for (const auto& poi : fPOIs) {
                    fitTool.ExportErrorDecompositionInTextFile(poi,fName+"/Fits/"+fInputName+fSuffix+"_errDecomp_"+poi+".txt");
                    fitTool.ExportErrorDecompositionGroupInTextFile(poi,fName+"/Fits/"+fInputName+fSuffix, fLumiLabel, fCmeLabel, fName, fHEPDataFormat);
                }
                this->PlotRankingFromCovMatrix(fitTool);
            }
        }
    }
    result = fitTool.ExportFitResultInMap();
    if (fBlindedParameters.size() > 0) std::cout.clear();

    if (fFitType == FitType::UNFOLDING && save) {
        ProduceSystSubCategoryMap();

        std::vector<std::string> categories;
        for (const auto& icategory : fSubCategoryImpactMap) {
            // only take unique categories
            if (std::find(categories.begin(), categories.end(), icategory.second) != categories.end()) continue;
            categories.emplace_back(icategory.second);
        }
        this->PlotUnfoldingErrors(fitTool, categories);
    }

    // grouped systematics impact
    if(fDoGroupedSystImpactTable){
        // name of file to write results to
        std::string outNameGroupedImpact = fName+"/Fits/GroupedImpact"+fSuffix;
        if(fBootstrap!="" && fBootstrapIdx>=0){
            gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
            outNameGroupedImpact = fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+"GroupedImpact"+fSuffix;
        }
        if(fGroupedImpactCategory!="all") outNameGroupedImpact += "_"+fGroupedImpactCategory;

        ProduceSystSubCategoryMap();                        // fill fSubCategoryImpactMap first
        fitTool.SetSystMap( fSubCategoryImpactMap );     // hand over the map to the FittingTool
        fitTool.GetGroupedImpact( mc, simPdf, data, ws, fGroupedImpactCategory, outNameGroupedImpact, fName, fLumiLabel, fCmeLabel, fHEPDataFormat);
    }

    return result;
}

//__________________________________________________________________________________
//
std::unique_ptr<RooWorkspace> TRExFit::PerformWorkspaceCombinationxRooFit(const std::vector <std::string>& regionsToFit) const {

    if (regionsToFit.empty()) {
        LOG(ERROR) << "List of regions to fit is empty\n";
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<TFile> rootFileCombined(nullptr);

    if(fBootstrap!="" && fBootstrapIdx>=0) {
        rootFileCombined.reset(TFile::Open( (fName+"/RooStats/"+fBootstrapSyst+fBootstrapSample+"_BSId"+Form("%d",fBootstrapIdx)+"/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read"));
    } else {
        rootFileCombined.reset(TFile::Open( (fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root").c_str(),"read"));
    }

    if (!rootFileCombined) {
        LOG(ERROR) << "Cannot read the WS\n";
        exit(EXIT_FAILURE);
    }

    auto measurement = std::unique_ptr<RooStats::HistFactory::Measurement>(static_cast<RooStats::HistFactory::Measurement*>(rootFileCombined -> Get( (fInputName+fSuffix).c_str())));

    if (!measurement) {
        LOG(ERROR) << "Cannot read the measurement object from the root file\n";
        exit(EXIT_FAILURE);
    }

    auto inputWS = std::unique_ptr<RooWorkspace>(rootFileCombined->Get<RooWorkspace>("combined"));
    if (!inputWS) {
        LOG(ERROR) << "Cannot read the workspace from the root file\n";
        exit(EXIT_FAILURE);
    }

    xRooNode node(*inputWS);

    std::string regions("");
    for (const auto& ireg : regionsToFit) {
        regions += ireg;
        regions += ",";
    }

    regions.resize(regions.size() - 1);

    auto reduced = node["pdfs/simPdf"]->reduced(regions);

    std::unique_ptr<RooWorkspace> result = std::make_unique<RooWorkspace>("combined","combined");
    xRooNode ws(*result);
    auto addedPdf = ws.Add(reduced);
    addedPdf.at(0)->coords();
    addedPdf.get<TNamed>()->SetName("simPdf");
    for(auto ds : reduced.datasets()) {
       ws.Add(*ds);
    }
    RooStats::ModelConfig mc("ModelConfig","ModelConfig",result.get());
    mc.SetPdf(addedPdf->GetName());
    mc.SetObservables(*addedPdf.robs().get<RooArgList>());
    mc.SetGlobalObservables(*addedPdf.globs().get<RooArgList>());
    mc.SetNuisanceParameters(*addedPdf.np().get<RooArgList>());
    mc.SetParametersOfInterest(*addedPdf.poi().get<RooArgList>());
    ws.Add(mc);

    RooStats::HistFactory::HistoToWorkspaceFactoryFast::ConfigureWorkspaceForMeasurement("simPdf", result.get(), *measurement );

    return result;
}

//__________________________________________________________________________________
//
void TRExFit::PlotFittedNP(){
    if (fBootstrap!="" && fBootstrapIdx>=0){
        ReadFitResults(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".root");
    } else {
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
    }
    if (fFitResults) {
        fFitResults->SetParsToHide(fVarNameHide);
        std::set < std::string > npCategories;

        for(const auto& syst : fSystematics) {
            if (syst->fCategory.find("gamma_") != std::string::npos) continue;
            npCategories.insert(syst->fCategory);
        }

        // Create directories for organisation
        gSystem->mkdir((fName+"/Pulls/All").c_str(), true);
        gSystem->mkdir((fName+"/Pulls/PullSig").c_str(), true);
        gSystem->mkdir((fName+"/Pulls/RankedPullSig").c_str(), true);
        // Add a category enclosing all NPs
        npCategories.insert("all");

        if (fStatOnly) {
            LOG(INFO) << "Stat only fit => No NP Pull plots generated.\n";
        } else {
            for(const std::string& cat : npCategories){
                std::string cat_for_name = cat;
                if (cat != "all"){
                    cat_for_name.insert(0, "_");
                    std::replace(cat_for_name.begin(), cat_for_name.end(), ' ', '_');
                    std::replace(cat_for_name.begin(), cat_for_name.end(), '#', '_');
                    std::replace(cat_for_name.begin(), cat_for_name.end(), '{', '_');
                    std::replace(cat_for_name.begin(), cat_for_name.end(), '}', '_');
                } else {
                    cat_for_name = "";
                }

                // Plot pulls
                fFitResults->DrawNPPulls(fName+"/Pulls/All/NuisPar"+cat_for_name+fSuffix, cat, fNormFactors, fShapeFactors, fBlindedParameters, FitResults::PullSigType::NoPullSig);
                // Plot pull significances
                fFitResults->DrawNPPulls(fName+"/Pulls/PullSig/NuisParSig"+cat_for_name+fSuffix, cat, fNormFactors, fShapeFactors, fBlindedParameters, FitResults::PullSigType::Regular);
                // Plot Ranked pull significances
                fFitResults->DrawNPPulls(fName+"/Pulls/RankedPullSig/NuisParSigRank"+cat_for_name+fSuffix, cat, fNormFactors, fShapeFactors, fBlindedParameters, FitResults::PullSigType::Rank);
            }

            fFitResults->DrawGammaShapePulls(fName+"/Gammas"+fSuffix, fBlindedParameters, fShapeFactors, false);
            fFitResults->DrawGammaShapePulls(fName+"/ShapeFactors"+fSuffix, fBlindedParameters, fShapeFactors, true);
        }
        // Check for EFT params
        std::vector<std::shared_ptr<NormFactor> > eftNFs;
        std::vector<std::shared_ptr<NormFactor> > nonEFTNFs;
        for(const auto& norm : fNormFactors){
            if(norm->fCategory == "EFT" ) {
                eftNFs.emplace_back(norm);
            } else {
                nonEFTNFs.emplace_back(norm);
            }
        }
        if(!eftNFs.empty()) {
            fFitResults->DrawNormFactors(fName+"/EFTParams"+fSuffix, eftNFs, fBlindedParameters);
            fFitResults->DrawNormFactors(fName+"/NormFactors"+fSuffix, nonEFTNFs, fBlindedParameters);
        } else {
            fFitResults->DrawNormFactors(fName+"/NormFactors"+fSuffix, fNormFactors, fBlindedParameters);
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::PlotCorrelationMatrix(){
    if(fFitType != FitType::UNFOLDING && (fStatOnly || fStatOnlyFit)) {
        LOG(INFO) << "Stat only fit => No Correlation Matrix generated.\n";
        return;
    }
    //plot the correlation matrix (considering only correlations larger than TRExFitter::CORRELATIONTHRESHOLD)
    if (fBootstrap!="" && fBootstrapIdx>=0){
        ReadFitResults(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".root");
    } else{
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
    }

    if(fFitResults){
        fFitResults->SetParsToHide(fVarNameHide);
        fFitResults->SetOutputFolder(fName);
        const std::string path = fName+"/CorrMatrix"+fSuffix;
        fFitResults->DrawCorrelationMatrix(path,
                                           fHEPDataFormat,
                                           TRExFitter::CORRELATIONTHRESHOLD, false,
                                           fPOIs);


        std::vector < std::string > EFTNFs;
        for(auto norm : fNormFactors){
          if(norm->fCategory == "EFT" )EFTNFs.emplace_back(norm->fName);
        }
        if(EFTNFs.size()){
            fFitResults->SetEFTNFs(EFTNFs);
            fFitResults->DrawCorrelationMatrix(path,
                                               fHEPDataFormat,
                                               TRExFitter::CORRELATIONTHRESHOLD, true,
                                               EFTNFs);
            }
    }
}

//__________________________________________________________________________________
//
void TRExFit::PlotUnfoldedData(const std::string& frPath,
                               const std::string& outputPath,
                               const std::string& wsFileName,
                               const std::string& wsName) {
    if (fFitType != TRExFit::FitType::UNFOLDING) return;
    // Read the correlation matrix
    ReadFitResults(frPath+".root");

    bool hasTau(false);
    for (const auto& inf : fNormFactors) {
        if (inf->fTau > 0) {
            hasTau = true;
            break;
        }
    }

    if (hasTau) {
        const std::vector<double> globalCorrelation = this->CalculateGlobalCorrelation();
        if (globalCorrelation.size() != fUnfolding.size()) {
            LOG(WARNING) << "Mismatch between size global correlations and unfoldings\n";
        } else {
            for (std::size_t iunf = 0; iunf < globalCorrelation.size(); ++iunf) {
                LOG(INFO) << "Global correlation coefficient for unfolding: " << fUnfolding.at(iunf)->fName << " is: " << globalCorrelation.at(iunf) << "\n";
            }
        }
    }

    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        PlotUnfoldedData(fUnfolding.at(i).get(), frPath, wsFileName, outputPath, wsName);
    }

    if (fUnfolding.size() > 1) {
        std::vector<const Unfolding*> unfoldings;
        std::vector<std::string> names;
        for (const auto& iunf: fUnfolding) {
            names.emplace_back(fName);
            unfoldings.emplace_back(iunf.get());
        }
        this->PlotMultipleUnfoldingCovariance(unfoldings, names, frPath, wsFileName, fName, wsName, fNormFactors);
    }

    // add the cross correlation terms
    if (fUnfolding.size() == 1) return;
    std::vector<std::string> bins;
    for (const auto& iunf : fUnfolding) {
        for (int ibin = 0; ibin < iunf->fNumberUnfoldingTruthBins; ++ibin) {
            const std::string binName = iunf->fName + "_Bin_" + Common::IntToFixLenStr(ibin+1) + "_mu";
            bins.emplace_back(binName);
        }
    }

    auto text = std::make_unique<std::ofstream>();
    text->open(outputPath +"/Fits/Correlations_UnfoldedResults.txt");
    if (!text->is_open() || !text->good()) {
        LOG(ERROR) << "Cannot open text file in: " << fName + "/Fits/Correlations_UnfoldedResults.txt\n";
        exit(EXIT_FAILURE);
    }

    *text << "List of Parameters\n";
    for (const auto& i : bins) {
        *text << i << "\n";
    }
    *text << "\n\nCorrelations\n";
    *text << std::fixed << std::setprecision(4);
    for (const auto& i : bins) {
        for (const auto& j : bins) {
            const double corr = (i == j) ? 1. : fFitResults->GetCorrelationMatrix()->GetCorrelation(i, j);
            *text << corr << " ";
        }
        *text << "\n";
    }
    text->close();

}
//__________________________________________________________________________________
//
void TRExFit::PlotUnfoldedData(const Unfolding* unfolding,
                               const std::string& frPath,
                               const std::string& wsFileName,
                               const std::string& outputPath,
                               const std::string& wsName) {

    LOG(INFO) << "Producing unfolded plots...\n";

    std::unique_ptr<TFile> input(TFile::Open((fName + "/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
    if (!input) {
        LOG(ERROR) << "Cannot read file from " + fName + "/UnfoldingHistograms/FoldedHistograms.root\n";
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<TH1> truth(dynamic_cast<TH1*>(input->Get((unfolding->fName+"_truth_distribution").c_str())));
    if (!truth) {
        LOG(ERROR) << "Cannot read the truth distribution\n";
        exit(EXIT_FAILURE);
    }
    truth->SetDirectory(nullptr);
    truth->Scale(fLumiScale);

    UnfoldingResult unfolded;
    unfolded.SetTruthDistribution(truth.get());
    unfolded.SetNormXSec(unfolding->fUnfoldNormXSec);

    UnfoldingResult unfoldedStat;
    unfoldedStat.SetTruthDistribution(truth.get());
    unfoldedStat.SetNormXSec(unfolding->fUnfoldNormXSec);

    std::unique_ptr<TGraphAsymmErrors> statError(nullptr);
    // Get all the relevant values
    Common::GetUnfoldingResult(unfolded, unfolding, truth.get(), fFitResults.get(), wsFileName, wsName, frPath+".root", fNormFactors, {}, false);
    if (fUnfoldingShowStat) {
        this->ReadStatOnlyFromErrorDecomposition(fPOIs, frPath + "_errDecomp");
        Common::GetUnfoldingResult(unfoldedStat, unfolding, truth.get(), fFitResults.get(), wsFileName, wsName, frPath+".root", fNormFactors, fStatOnlyErrorDecomposition, true);
        statError = unfoldedStat.GetUnfoldedResultErrorBand();
    }

    std::unique_ptr<TH1D> data               = unfolded.GetUnfoldedResult();
    std::unique_ptr<TGraphAsymmErrors> error = unfolded.GetUnfoldedResultErrorBand();

    YamlConverter converter{};
    converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
    converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
    converter.WriteUnfolding(error.get(), outputPath, unfolding->fName);

    PlotUnfold(data.get(), error.get(), statError.get(), unfolding, outputPath, converter);

    // Dump results into a text file
    auto text = std::make_unique<std::ofstream>();
    text->open(outputPath+"/Fits/"+unfolding->fName+"_UnfoldedResults.txt");
    if (!text->is_open() || !text->good()) {
        LOG(ERROR) << "Cannot open text file in: " << outputPath+unfolding->fName+"_UnfoldedResults.txt\n";
        exit(EXIT_FAILURE);
    }

    if (unfolding->fUnfoldingDivideByLumi > 0 && !unfolding->fUnfoldNormXSec) {
        unfolded.SetDivideByLumi(unfolding->fUnfoldingDivideByLumi);
    }
    unfolded.SetDivideByBinWidth(unfolding->fUnfoldingDivideByBinWidth);
    unfolded.DumpResults(text.get());
    *text << std::fixed << std::setprecision(4);
    *text << "\n\n Correlations \n";
    const int n = unfolding->fNumberUnfoldingTruthBins;
    for (int i = 0; i < n; ++i) {
        const std::string name_i = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
        for (int j = 0; j < n; ++j) {
            const std::string name_j = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(j+1) + "_mu";
            const double corr = (i == j) ? 1. : fFitResults->GetCorrelationMatrix()->GetCorrelation(name_i, name_j);
            *text << corr << " ";
        }
        *text << "\n";
    }

    this->PlotCovarianceMatrix(unfolding, unfolded, outputPath);

    text->close();

    input->Close();
}

//__________________________________________________________________________________
//
void TRExFit::GetLimit() {
    //
    // Make sure a proper POI is present
    //
    if(fPOIforLimit=="") {
        if(fPOIs.size() == 1) {
            fPOIforLimit = fPOIs.at(0);
        } else {
            LOG(ERROR) << "No POI specified (in 'Limit' block).\n";
            return;
        }
    }

    if (fFitType == TRExFit::FitType::BONLY) {
        LOG(WARNING) << "Running limit estimate with \"FitType: BONLY\". This setting will be ignored for the limit estimate.\n";
    }

    //
    // Read NPvalues from fit-result file
    //
    if (fFitNPValuesFromFitResultsFile != "") {
        LOG(INFO) << "Setting NP values for Asimov data-set creation from fit results stored in file " << fFitNPValuesFromFitResultsFile << "...\n";
        fFitNPValues = NPValuesFromFitResultsFile(fFitNPValuesFromFitResultsFile);
    }

    //
    // If a workspace file name is specified, do simple limit
    //
    gSystem->mkdir((fName+"/Limits").c_str());
    std::string limitDir = fName + "/Limits/" + (fLimitType == LimitType::ASYMPTOTIC? "Asymptotics" : "Toys");
    gSystem->mkdir(limitDir.c_str());
    std::string dataName = "obsData";
    if(fLimitIsBlind) dataName = "asimovData";
    const std::string limitFile = limitDir + "/" + fLimitOutputPrefixName+fSuffix +".root";
    const std::string limitWSfile = limitDir + "/PostLimitWS" + fSuffix + ".root";
    if(fWorkspaceFileName!=""){
        std::unique_ptr<TFile> f(TFile::Open(fWorkspaceFileName.c_str()));
        RooAbsData* data = f->Get<RooAbsData>(dataName.c_str());
        RooWorkspace* ws = f->Get<RooWorkspace>("combined");
        if (fLimitType == LimitType::ASYMPTOTIC ||
            (fLimitType == LimitType::TOYS &&
             fLimitToysUsexRooFit != ToysUsexRooFit::FALSE)) {
            RunLimit(ws, dataName, limitFile, limitWSfile);
        } else {
            RunLimitToys(data, ws);
        }
        f->Close();
    } else {
        auto wsAndData = this->PrepareMixedDataset(TRExFit::WorkspaceType::LIMIT);

        if (this->DoingMixedFitting()) {
            LOG(INFO) << "Using mixed Data/Asimov dataset in the limit estimate\n";
            const std::string newDataName = "mixedAsimov";
            dataName = newDataName;
            wsAndData.second->SetName(newDataName.c_str());
            wsAndData.first->import(*wsAndData.second);
            fLimitIsBlind = true;
        }

        if (fLimitType == LimitType::ASYMPTOTIC ||
            (fLimitType == LimitType::TOYS &&
             fLimitToysUsexRooFit != ToysUsexRooFit::FALSE)) {
            RunLimit(wsAndData.first.get(), dataName, limitFile, limitWSfile);
        } else {
            RunLimitToys(wsAndData.second.get(), wsAndData.first.get());
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::GetSignificance() {
    //
    // Make sure a proper POI is present
    //
    if (fPOIforSig == "") {
        if (fPOIs.size() == 1) {
            fPOIforSig = fPOIs.at(0);
        } else {
            LOG(ERROR) << "No POI specified (in 'Significance' block).\n";
            return;
        }
    }

    if (fFitType == TRExFit::FitType::BONLY) {
        LOG(WARNING) << "Running significance estimate with \"FitType: BONLY\". This setting will be ignored for the significance estimate.\n";
    }

    //
    // Read NPvalues from fit-result file
    //
    if (fFitNPValuesFromFitResultsFile != "") {
        LOG(INFO) << "Setting NP values for Asimov data-set creation from fit results stored in file " << fFitNPValuesFromFitResultsFile << "...\n";
        fFitNPValues = NPValuesFromFitResultsFile(fFitNPValuesFromFitResultsFile);
    }

    //
    // If a workspace file name is specified, do simple significance
    //
    std::string dataName = "obsData";
    if(fSignificanceIsBlind) dataName = "asimovData";
    gSystem->mkdir((fName+"/Significance").c_str());
    std::string significanceDir = fName + "/Significance/" +
        (fSignificanceType == SignificanceType::ASYMPTOTIC? "Asymptotics" : "Toys");
    gSystem->mkdir(significanceDir.c_str());
    const std::string significanceFile = significanceDir + "/" + fSignificanceOutputPrefixName+fSuffix +".root";
    const std::string significanceWSfile = significanceDir + "/PostSignificanceWS" + fSuffix + ".root";
    if(fWorkspaceFileName!=""){
        std::unique_ptr<TFile> f(TFile::Open(fWorkspaceFileName.c_str()));
        RooAbsData* data = f->Get<RooAbsData>(dataName.c_str());
        RooWorkspace* ws = f->Get<RooWorkspace>("combined");
        if (fSignificanceType == SignificanceType::ASYMPTOTIC ||
            (fSignificanceType == SignificanceType::TOYS &&
             fSignificanceToysUsexRooFit)) {
            RunSignificance(ws, dataName, significanceFile, significanceWSfile);
        } else {
            RunSignificanceToys(data, ws);
        }
        f->Close();
    } else {
        auto wsAndData = this->PrepareMixedDataset(TRExFit::WorkspaceType::SIGNIFICANCE);

        if (this->DoingMixedFitting()) {
            LOG(INFO) << "Using mixed Data/Asimov dataset in the significance estimate\n";
            const std::string newDataName = "mixedAsimov";
            dataName = newDataName;
            wsAndData.second->SetName(newDataName.c_str());
            wsAndData.first->import(*wsAndData.second);
            fSignificanceIsBlind = true;
        }

        if (fSignificanceType == SignificanceType::ASYMPTOTIC ||
            (fSignificanceType == SignificanceType::TOYS &&
             fSignificanceToysUsexRooFit != ToysUsexRooFit::FALSE)) {
            RunSignificance(wsAndData.first.get(), dataName, significanceFile, significanceWSfile);
        } else {
            RunSignificanceToys(wsAndData.second.get(), wsAndData.first.get());
        }
    }
}

//__________________________________________________________________________________
//
bool TRExFit::ReadFitResults(const std::string& fileName) {
    LOG(INFO) << "------------------------------------------------------\n";
    LOG(INFO) << "Reading fit results from file\n";
    fFitResults.reset(new FitResults());
    fFitResults->SetPOIPrecision(fPOIPrecision);
    fFitResults->SetPlotLabel(fPlotLabel);

    bool exist(false);

    if (fileName.find(".root")!=std::string::npos) {
        exist = fFitResults->ReadFromRootFile(fileName);
    }
    // make a list of systematics from all samples...
    // ...
    // assign to each NP in the FitResults a title, and a category according to the syst in the fitter

    // note: some NPs are assigned to multiple systematics (those which are correlated)
    // we will just keep overwriting, so the title and catagory
    // will be from the last systematic with that NP name

    for(const auto& inp : fFitResults->GetNuisanceParameters()) {

        for(unsigned int j_sys = 0; j_sys < fSystematics.size(); ++j_sys) {
            // the systematic fName doesn't necessarily equate to the NP fName
            // compare the fSystematics[j_sys]->fNuisanceParameter instead!

            if (fSystematics[j_sys]->fNuisanceParameter == inp.second->fName) {
                inp.second->fTitle = fSystematics[j_sys]->fTitle;
                inp.second->fCategory = fSystematics[j_sys]->fCategory;
            }
        }
        for (unsigned int j = 0; j < fNormFactors.size(); j++) {
            if (fNormFactors[j]->fName == inp.second->fName) {
                inp.second->fTitle = fNormFactors[j]->fTitle;
                inp.second->fCategory = fNormFactors[j]->fCategory;
            }
        }
        // FIXME SF probably there are several NPs associated to it
        for (unsigned int j = 0; j < fShapeFactors.size(); j++) {
            if (fShapeFactors[j]->fName == inp.second->fName) {
                inp.second->fTitle = fShapeFactors[j]->fTitle;
                inp.second->fCategory = fShapeFactors[j]->fCategory;
            }
        }
    }

    // create an ordered list of NPs
    if (fReorderNPs) {
        fFitResults->ClearNuisParList();
        for(auto norm : fNormFactors){
            if (norm->fExpression.first=="") {
                if(Common::FindInStringVector(fFitResults->GetNuisParList(),norm->fNuisanceParameter)<0) fFitResults->AddToNuisParList(norm->fNuisanceParameter);
            }
        }
        for(auto syst : fSystematics){
            // if non-SHAPE
            if(syst->fType!=Systematic::SHAPE){
                if(Common::FindInStringVector(fFitResults->GetNuisParList(),syst->fNuisanceParameter)<0) fFitResults->AddToNuisParList(syst->fNuisanceParameter);
            }
            // FIXME: something to add for actual SHAPE systematics (not "separate gammas")
        }
        // FIXME: something to add for shape factors
    }
    return exist;
}

//__________________________________________________________________________________
//
void TRExFit::PrintConfigSummary() const{
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Job name: " << fName << "\n";
    LOG(INFO) << "Reading the following regions:\n";
    for(const auto& ireg : fRegions) {
        ireg->Print();
    }
    LOG(INFO) << "Reading the following samples:\n";
    for(const auto& isample : fSamples) {
        LOG(INFO) << "     " << isample->fName << "\n";
    }
    LOG(INFO) << "Reading the following systematics:\n";
    std::vector<std::string> tmp{};
    for(const auto& isyst : fSystematics) {
        if (std::find(tmp.begin(), tmp.end(), isyst->fName) == tmp.end()){
            LOG(INFO) << " " << isyst->fName << "\n";
            tmp.emplace_back(isyst->fName);
        }
    }
    LOG(INFO) << "-------------------------------------------\n";
}

//__________________________________________________________________________________
//
std::shared_ptr<Sample> TRExFit::GetSample(const std::string& name) const{
    for(unsigned int i=0;i<fSamples.size();i++){
        if(fSamples[i]->fName == name) return fSamples[i];
    }
    return nullptr;
}

//__________________________________________________________________________________
//
std::size_t TRExFit::GetSampleIndex(const std::string& name) const{
    for(std::size_t i=0; i<fSamples.size(); ++i){
        if(fSamples[i]->fName == name) return i;
    }
    return 99999;
}

//__________________________________________________________________________________
//
void TRExFit::DrawAndSaveSeparationPlots() const{

    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/Plots").c_str());
    gSystem->mkdir((fName+"/Plots/Separation").c_str());


    // loop over regions
    for(unsigned int i_ch=0; i_ch < fRegions.size(); i_ch++){
        // begin plotting
        TCanvas dummy3 ("dummy3", "dummy3", 600,600);
        dummy3.cd();

        if(fSeparationPlot.size()==0) {
            if(fRegions[i_ch]->fSig.size() ==  0){
                LOG(ERROR) << "No Signal found in region " << fRegions[i_ch]->fName << "\n";
                continue;
            }
        }

        TLegend legend3(0.55,0.77-0.05*std::max(0,int(fSeparationPlot.size()-4)),0.94,0.87);
        legend3.SetTextFont(gStyle->GetTextFont());
        legend3.SetTextSize(gStyle->GetTextSize());

        if(fSeparationPlot.size()>0) {
            std::vector<std::unique_ptr<TH1D> > templ;
            float Ymaximum = 0;
            int linestyle = 0;
            for(const auto& isample : fRegions[i_ch]->fSampleHists) {
                std::unique_ptr<TH1D> tmp_smp(static_cast<TH1D*>(isample->fHist->Clone()));
                tmp_smp->SetDirectory(nullptr);
                if (Common::FindInStringVector(fSeparationPlot, isample->fSample->fName)<0) continue;

                int linecolor = tmp_smp->GetFillColor();
                tmp_smp->SetLineColor( linecolor );
                tmp_smp->SetMarkerColor( linecolor );
                tmp_smp->SetLineWidth( 3 );
                tmp_smp->SetFillStyle( 0 );
                tmp_smp->SetLineStyle( ++linestyle );
                tmp_smp->SetMarkerSize(0);
                tmp_smp->Scale(1./tmp_smp->Integral());

                legend3.AddEntry(tmp_smp.get(), isample->fSample->fTitle.c_str() , "l");
                Ymaximum = (tmp_smp->GetMaximum() > Ymaximum) ? tmp_smp->GetMaximum() : Ymaximum;

                templ.emplace_back(std::move(tmp_smp));
            }
            if(templ.size()==0)
                continue;

            legend3.SetFillStyle(0) ;
            legend3.SetBorderSize(0);

            std::string xaxis = fRegions[i_ch]->fVariableTitle;

            templ[0]->GetYaxis()->SetTitle("Arbitrary units");
            templ[0]->GetXaxis()->SetTitle(xaxis.c_str());
            templ[0]->GetYaxis()->SetTitleOffset(1.6);
            templ[0]->GetYaxis()->SetNdivisions(506);
            templ[0]->GetYaxis()->SetRangeUser(0.,Ymaximum*1.5);
            templ[0]->Draw("hist e");
            for(unsigned int i_smp=0; i_smp< templ.size(); i_smp++){
                templ[i_smp]->Draw("histsame e");
            }
            legend3.Draw("same");

            myText(0.20,0.78,1,fLabel.c_str());
            myText(0.20,0.73,1,fRegions[i_ch]->fLabel.c_str());

            std::string cme = fRegions[i_ch]->fCmeLabel;
            std::string lumi = fRegions[i_ch]->fLumiLabel;

            myText(0.20,0.83,1,Form("#sqrt{s} = %s, %s", cme.c_str(), lumi.c_str()));
            if(fPlotLabel!="none") TRExLabelNew(0.20,0.84+0.04, TRExFitter::EXPERIMENT_LABEL.c_str(), (fPlotLabel+"  Simulation").c_str(), kBlack, gStyle->GetTextSize());

            std::ostringstream SEP;
            SEP.precision(3);
            if(templ.size()>1) SEP << "Separation: " << Common::GetSeparation(templ[0].get(),templ[1].get())*100 << "%";
            myText(0.55,0.73-0.05*std::max(0,int(fSeparationPlot.size()-4)),1,SEP.str().c_str());

            LOG(INFO) << SEP.str() << " in region " << fRegions[i_ch]->fName << "\n";

            Common::SaveCanvasAs(dummy3, fName+"/Plots/Separation/"+fRegions[i_ch]->fName+fSuffix);
        }
        else {
            std::unique_ptr<TH1D> sig(static_cast<TH1D*>(fRegions[i_ch]->fSig[0]->fHist->Clone()));

            std::unique_ptr<TH1D> bkg (static_cast<TH1D*>(fRegions[i_ch]->fBkg[0]->fHist->Clone())); // clone the first bkg
            for(std::size_t i_bkg=1; i_bkg< fRegions[i_ch] -> fBkg.size(); i_bkg++){
                bkg->Add(fRegions[i_ch]->fBkg[i_bkg]->fHist.get()); // add the rest
            }

            sig->SetLineColor( 2 );
            sig->SetLineWidth( 3 );
            sig->SetFillStyle( 0 );
            sig->SetLineStyle( 2 );

            bkg->SetLineColor( kBlue );
            bkg->SetLineWidth( 3 );
            bkg->SetFillStyle( 0 );
            bkg->SetLineStyle( 1 );

            legend3.AddEntry(bkg.get(), "Total background" , "l");
            legend3.AddEntry(sig.get(), fRegions[i_ch]->fSig[0]->fSample->fTitle.c_str() , "l");
            legend3.SetFillStyle(0) ;
            legend3.SetBorderSize(0);

            std::string xaxis = fRegions[i_ch]->fVariableTitle;

            sig->GetYaxis()->SetTitle("Arbitrary units");
            sig->GetXaxis()->SetTitle(xaxis.c_str());

            sig->GetYaxis()->SetTitleOffset(1.6);

            bkg->GetYaxis()->SetTitle("Arbitrary units");
            bkg->GetXaxis()->SetTitle(xaxis.c_str());

            bkg->GetYaxis()->SetTitleOffset(1.6);

            sig->GetYaxis()->SetNdivisions(506);
            bkg->GetYaxis()->SetNdivisions(506);

            sig->Scale(1./sig->Integral());
            bkg->Scale(1./bkg->Integral());


            if(bkg->GetMaximum() > sig->GetMaximum()){
                bkg->GetYaxis()->SetRangeUser(0.,bkg->GetMaximum()*1.5);
                bkg->Draw("hist");
                sig->Draw("histsame");
            }
            else {
                sig->GetYaxis()->SetRangeUser(0.,sig->GetMaximum()*1.5);
                sig->Draw("hist");
                bkg->Draw("histsame");
                sig->Draw("histsame");
            }

            legend3.Draw("same");

            if (fLabel != "none") myText(0.20,0.78,1,fLabel.c_str());
            myText(0.20,0.73,1,fRegions[i_ch]->fLabel.c_str());

            std::string cme = fRegions[i_ch]->fCmeLabel;
            std::string lumi = fRegions[i_ch]->fLumiLabel;

            myText(0.20,0.83,1,Form("#sqrt{s} = %s, %s", cme.c_str(), lumi.c_str()));

            if (fPlotLabel!="none") TRExLabelNew(0.20,0.84+0.04, TRExFitter::EXPERIMENT_LABEL.c_str(), (fPlotLabel+"  Simulation").c_str(), kBlack, gStyle->GetTextSize());

            std::ostringstream sep;
            sep.precision(3);
            sep << "Separation: " << Common::GetSeparation(sig.get(),bkg.get())*100 << "%";
            myText(0.55,0.73,1,sep.str().c_str());

            LOG(INFO) << sep.str() << " in region " << fRegions[i_ch]->fName << "\n";

            Common::SaveCanvasAs(dummy3, fName+"/Plots/Separation/"+fRegions[i_ch]->fName+fSuffix);
        }
    }// regions

   return;
}

//____________________________________________________________________________________
//
void TRExFit::ProduceNPRanking(const std::string& NPnames) {

    if(fFitType==BONLY){
        LOG(ERROR) << "For ranking plots, the SPLUSB FitType is needed.\n";
        exit(EXIT_FAILURE);
    }

    if (fBootstrap!="" && fBootstrapIdx>=0){
        ReadFitResults(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".root");
    } else{
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
    }

    std::string frFile("");
    if (fFitResultsRootFile != "") {
        frFile = fFitResultsRootFile;
    } else {
        frFile = fName + "/Fits/"+fInputName+fSuffix+".root";
    }

    xRooNode rooNode = Common::ReadWSandFitResults(fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root",
                                                   frFile);

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf for postfit\n";
        exit(EXIT_FAILURE);
    }

    xRooNode rooNodePrefit = Common::ReadWSandFitResults(fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root",
                                                         "");

    auto pdfPrefit = rooNodePrefit["simPdf"];
    if (!pdfPrefit) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf for prefit\n";;
        exit(EXIT_FAILURE);
    }

    //
    // List of systematics to check
    //
    RankingManager manager{};
    manager.SetUseHesseBeforeMigrad(fUseHesseBeforeMigrad);
    manager.SetMaximumFCNcalls(fMaximumNumberFCNcalls);
    manager.SetRandomNPs(fRndRange, fUseRnd, fRndSeed);
    manager.SetToleranceScale(fToleranceScale);
    manager.SetShiftGlobalObservables(fShiftGlobalObservablesInRanking);
    for (const auto& iparameter : pdf->floats()) {
        const std::string name = iparameter->GetName();
        if (NPnames != "all" && name.find(NPnames) == std::string::npos) continue;

        // check if it is NF
        auto itr = std::find_if(fNormFactors.begin(), fNormFactors.end(), [&name](const auto& nf){return name == nf->fName;});
        if (itr != fNormFactors.end()) {
            if (fShiftGlobalObservablesInRanking) continue;
            if (std::find(fPOIs.begin(), fPOIs.end(), name) != fPOIs.end()) {
                if (fUsePOISinRanking) {
                    manager.AddNuisPar(name, true);
                }
            }
            manager.AddNuisPar(name, true);
        } else {
            //check if it is ShapeFactor
            auto itrSF = std::find_if(fShapeFactors.begin(), fShapeFactors.end(), [&name](const auto& sf){return name.find(sf->fName) != std::string::npos;});
            if (itrSF != fShapeFactors.end()) {
                if (fShiftGlobalObservablesInRanking) continue;
                manager.AddNuisPar(name, true);
            } else {
                manager.AddNuisPar(name, false);
            }
        }
    }

    //
    // Text files containing information necessary for drawing of ranking plot
    //
    std::string outName = fName+"/Fits/NPRanking"+fSuffix;
    if(fBootstrap!="" && fBootstrapIdx>=0){
        gSystem -> mkdir((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)).c_str(),true);
        outName = fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+"NPRanking"+fSuffix;
    }
    if(NPnames!="all") outName += "_"+NPnames;
    manager.SetOutputPath(outName);

    //
    // Creating the combined model
    //
    std::unique_ptr<RooDataSet> data(nullptr);
    std::unique_ptr<RooWorkspace> ws(nullptr);
    std::unique_ptr<TFile> customWSfile(nullptr);

    //
    // Read NPvalues from fit-result file
    //
    if (fFitNPValuesFromFitResultsFile != ""){
        LOG(INFO) << "Setting NP values for Asimov data-set creation from fit results stored in file " << fFitNPValuesFromFitResultsFile << "...\n";
        fFitNPValues = NPValuesFromFitResultsFile(fFitNPValuesFromFitResultsFile);
    }

    //
    // If there's a workspace specified, go on with simple fit, without looking for separate workspaces per region
    //
    if(fWorkspaceFileName!=""){
        customWSfile.reset(TFile::Open(fWorkspaceFileName.c_str(),"read"));
        ws = std::unique_ptr<RooWorkspace>(dynamic_cast<RooWorkspace*>(customWSfile->Get("combined")));
        if (!ws) {
            LOG(ERROR) << "Cannot read the custom WS!\n";
            return;
        }
        if(!fFitIsBlind) data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ws->data("obsData")));
        else             data = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(ws->data("asimovData")));
        if (!data) {
            LOG(ERROR) << "Cannot read the custom data from WS!\n";
            return;
        }
    }
    else{
        auto wsAndData = this->PrepareMixedDataset(TRExFit::WorkspaceType::FIT);
        ws = std::move(wsAndData.first);
        data = std::move(wsAndData.second);
    }

    const RooStats::ModelConfig *mc = dynamic_cast<const RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Model config is nullptr\n";
        exit(EXIT_FAILURE);
    }

    if (manager.GetNPsize() == 0) {
        LOG(WARNING) << "No parameter provided for ranking, maybe pruned?\n";
        LOG(WARNING) << "Will not produce any ranking output file for this fit\n";
        return;
    }

    manager.SetInjectGlobalObservables(fInjectGlobalObservables);
    manager.SetNPValues(fFitNPValues);
    manager.SetFixedNPs(fFitFixedNPs);
    manager.SetFitStrategy(fFitStrategy);
    manager.SetRegularizationType(fRegularizationType);

    if (fPOIs.size() > 0) {
        manager.SetPOINames(fPOIs);
    } else {
        LOG(WARNING) << "No POI set. Not able to produce ranking.\n";
        return;
    }
    manager.SetNCPU(fCPU);
    manager.SetStatOnly(fStatOnly);
    manager.SetUsePOISinRanking(fUsePOISinRanking);
    manager.RunRanking(pdf, pdfPrefit, ws.get(), data.get(), fNormFactors);

    if(customWSfile!=nullptr) {
        ws.reset(nullptr);
        data.reset(nullptr);
        customWSfile->Close();
    }
}

//____________________________________________________________________________________
//
void TRExFit::PlotNPRankingManager() {
    if(fRankingPlot=="MERGE"  || fRankingPlot=="ALL") PlotNPRanking(true,true);
    if(fRankingPlot=="SYSTS"  || fRankingPlot=="ALL") PlotNPRanking(true,false);
    if(fRankingPlot=="GAMMAS" || fRankingPlot=="ALL") PlotNPRanking(false,true);
}

//____________________________________________________________________________________
//
void TRExFit::PlotNPRanking(const bool flagSysts, const bool flagGammas) {

    for (const auto& poi : fPOIs) {
        const std::string fileToRead = fName+"/Fits/NPRanking"+fSuffix+"_"+poi+".txt";
        std::ifstream in(fileToRead.c_str());
        if (!in.good()) { // file doesnt exist
            const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fName+"/Fits/","NPRanking"+fSuffix+"_", poi);
            Common::MergeTxTFiles(inPaths, fileToRead);
        }
    }

    RankingManager manager{};
    manager.SetUseHesseBeforeMigrad(fUseHesseBeforeMigrad);
    manager.SetPlotLabel(fPlotLabel);
    manager.SetLumiLabel(fLumiLabel);
    manager.SetCmeLabel(fCmeLabel);
    manager.SetUseHEPDataFormat(fHEPDataFormat);
    manager.SetName(fName);
    manager.SetMaxNPPlot(fRankingMaxNP);
    manager.SetRankingCanvasSize(fNPRankingCanvasSize);
    manager.SetShiftGlobalObservables(fShiftGlobalObservablesInRanking);

    std::string extra("");
    if (!flagSysts || !flagGammas) {
        extra = !flagSysts ? "_gammas" : "_systs";
    }

    for (std::size_t ipoi = 0; ipoi < fPOIs.size(); ++ipoi) {
        const std::string fileToRead = fName+"/Fits/NPRanking"+fSuffix+"_"+fPOIs.at(ipoi)+".txt";
        manager.SetOutputPath(fileToRead);
        manager.SetSuffix(fSuffix+"_"+fPOIs.at(ipoi)+extra);
        manager.SetRankingPOIName(fRankingPOIName.at(ipoi));
        manager.SetUpperAxisNdivisions(fRankingUpperAxisNdivision.at(ipoi));
        manager.SetRankingPOIAxisScale(fRankingPOIAxisScale.at(ipoi));
        manager.ReadRankingResults(flagSysts, flagGammas);
        manager.PlotRanking(fRegions, fNormFactorNames, fShapeFactorNames, flagSysts, flagGammas);
        if (fShiftGlobalObservablesInRanking) {
            this->ProduceSystSubCategoryMap();
            manager.SetSuffix(fSuffix+"_"+fPOIs.at(ipoi));
            manager.ProduceGroupedImpact(fSubCategoryImpactMap);
        }
    }

    std::vector<std::pair<std::string,double> > parameters;
    if (fFitType == FitType::UNFOLDING) {
        for (const auto& iunfold : fUnfolding) {
            // read the truth file
            std::unique_ptr<TFile> file(TFile::Open((fName + "/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
            if (!file) {
                LOG(ERROR) << "Cannot open the truth file\n";
                exit(EXIT_FAILURE);
            }

            std::unique_ptr<TH1> histo(dynamic_cast<TH1*>(file->Get((iunfold->fName + "_truth_distribution").c_str())));
            if (!histo) {
                LOG(ERROR) << "Cannot find the truth histogram\n";
                exit(EXIT_FAILURE);
            }

            if (iunfold->fUnfoldNormXSec) {
                histo->Scale(1./histo->Integral());
            }

            if (iunfold->fUnfoldingDivideByBinWidth) {
                Common::ScaleByBinWidth(histo.get());
            }

            if ((iunfold->fUnfoldingDivideByLumi > 0) && !(iunfold->fUnfoldNormXSec)) {
                histo->Scale(1./iunfold->fUnfoldingDivideByLumi);
            }

            for (int i = 0; i < iunfold->fNumberUnfoldingTruthBins; ++i) {
                const std::string name = iunfold->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
                auto itr = std::find(fPOIs.begin(), fPOIs.end(), name);
                if (itr != fPOIs.end()) {
                    const double yield = histo->GetBinContent(i+1);
                    parameters.emplace_back(std::make_pair(name, yield));
                }
            }
        }

        // now add the other POIs
        for (const auto& ipoi : fPOIs) {
            auto itr = std::find_if(parameters.begin(), parameters.end(), [&ipoi](const auto& element){return element.first == ipoi;});
            if (itr != parameters.end()) continue;
            parameters.emplace_back(std::make_pair(ipoi, 1.));
        }
    } else {
        for (const auto& ipoi : fPOIs) {
            parameters.emplace_back(std::make_pair(ipoi, 1.));
        }
    }

    if (fDoExtendedCovariances) {
        // check if the statOnly fit result is there
        std::string tmp = fName;
        auto pos =  tmp.find(fDir);
        if (pos != std::string::npos) {
            tmp.erase(pos, fDir.length());
        }
        this->ReadStatOnlyFromErrorDecomposition(fPOIs, fName + "/Fits/" + tmp + "_errDecomp");
        manager.DumpCovariances(parameters, fNormFactorNames, fName, fFitResults.get(), fStatOnlyErrorDecomposition, fFitType == FitType::UNFOLDING);
    }

    if (fCombinerFormat) {
        if (fBootstrap!="" && fBootstrapIdx>=0){
            this->ReadFitResults(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".root");
        } else{
            this->ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
        }
        this->ReadStatOnlyFromErrorDecomposition(fPOIs, fName+"/Fits/"+fInputName+fSuffix + "_errDecomp");
        std::vector<std::string> params{};
        auto rankingVec = manager.GetRankingImpacts(fName, parameters, fNormFactorNames, params);
        YamlConverter::Measurement measurement;
        measurement.pois = fPOIs;
        std::vector<double> values;
        for (const auto& ipoi : fPOIs) {
            values.emplace_back(fFitResults->GetNuisParValue(ipoi));
        }
        measurement.values = values;
        for (const auto& inp : params) {
            auto itrNF = std::find(fNormFactorNames.begin(), fNormFactorNames.end(), inp);
            if (itrNF != fNormFactorNames.end()) {
                LOG(ERROR) << "NP: " << inp << " is NormFactor. Please fix.\n";
                exit(EXIT_FAILURE);
            }
            YamlConverter::Systematic syst;
            std::string name = Common::ReplaceString(inp,"alpha_","");
            name = Common::ReplaceString(name,"gamma_","");
            syst.name = name;
            if (inp.find("gamma_") != std::string::npos) {
                syst.pull = 0.;
                syst.constraint = 1.;
            } else {
                syst.pull = fFitResults->GetNuisParValue(name);
                const double up = fFitResults->GetNuisParErrUp(name);
                const double down = fFitResults->GetNuisParErrDown(name);
                syst.constraint = 0.5*(std::abs(up) + std::abs(down));
            }
            std::vector<double> impact;
            for (const auto& ipoi : rankingVec) {
                auto itr = ipoi.find(inp);
                if (itr == ipoi.end()) {
                    LOG(ERROR) << "Cannot find NP: " << name << "\n";
                    exit(EXIT_FAILURE);
                }
                if (inp.find("gamma_") != std::string::npos) {
                  impact.emplace_back(0.5*(itr->second.poihi - itr->second.poilo));
                } else {
                  impact.emplace_back(0.5*(itr->second.poiprehi - itr->second.poiprelo));
                }
            }
            syst.impact = impact;
            measurement.systematics.emplace_back(std::move(syst));
        }
        TMatrixDSym cor(params.size());
        for (std::size_t i = 0; i < params.size(); ++i) {
            for (std::size_t j = i; j < params.size(); ++j) {
                std::string param_i = Common::ReplaceString(params.at(i),"alpha_","");
                param_i = Common::ReplaceString(param_i,"gamma_","");
                std::string param_j = Common::ReplaceString(params.at(j),"alpha_","");
                param_j = Common::ReplaceString(param_j,"gamma_","");
                double corr = fFitResults->GetCorrelationMatrix()->GetCorrelation(param_i, param_j);
                cor(i,j) = corr;
                if (i != j) {
                    cor(j,i) = corr;
                }
            }
        }
        measurement.npCorrelation.ResizeTo(cor);
        measurement.npCorrelation = cor;

        // stat only covariance matrices
        std::string tmp = fName;
        auto pos =  tmp.find(fDir);
        if (pos != std::string::npos) {
            tmp.erase(pos, fDir.length());
        }
        TMatrixDSym cov(fPOIs.size());
        for (std::size_t i = 0; i < fPOIs.size(); ++i) {
            for (std::size_t j = 0; j < fPOIs.size(); ++j) {
                const double corr = fFitResults->GetCorrelationMatrix()->GetCorrelation(fPOIs.at(i), fPOIs.at(j));
                const double sigma_i = this->GetStatOnlyErrorFromDecomposition(fPOIs.at(i));
                const double sigma_j = this->GetStatOnlyErrorFromDecomposition(fPOIs.at(j));
                cov(i,j) = sigma_i * corr * sigma_j;
            }
        }
        measurement.statCov.ResizeTo(cov);
        measurement.statCov = cov;

        YamlConverter converter{};
        converter.WriteCombiner(fName, measurement, fName + "/CombinerOutput_" + fName + ".yml");
    }
}

//____________________________________________________________________________________
//
void TRExFit::PrintSystTables(std::string opt, const std::string& wsPath) const{
    LOG(INFO) << "Printing syst tables\n";
    if(fCleanTables) opt += "clean";
    if(fSystCategoryTables) opt += "category";
    if(fTableOptions.find("STANDALONE")!=std::string::npos) opt += "standalone";
    if(fTableOptions.find("LANDSCAPE")!=std::string::npos) opt +="landscape";
    if(fTableOptions.find("FOOTNOTESIZE")!=std::string::npos) opt +="footnotesize";

    const bool isPostFit = opt.find("post") != std::string::npos;

    std::string frFile("");
    if (isPostFit) {
        if (fFitResultsRootFile != "") {
            frFile = fFitResultsRootFile;
        } else {
            frFile = fName + "/Fits/"+fInputName+fSuffix+".root";
        }
    }

    xRooNode rooNode = Common::ReadWSandFitResults(wsPath, frFile);

    auto pdf = rooNode["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read ModelConfig: simPdf\n";
        exit(EXIT_FAILURE);
    }
    for(const auto& ireg : fRegions) {
        ireg->PrintSystTable(pdf, opt);
    }
}

//____________________________________________________________________________________
// this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter
void TRExFit::MergeSystematics(){
    // loop on systematics, see if any of them has name != nuisance
    for(const auto& syst : fSystematics){
        if(syst->fName == syst->fNuisanceParameter) continue;
        // if so, loop on other systematics to find one with name = nuisance parameter
        for(const auto& syst1 : fSystematics){
            if(syst->fName==syst1->fName) continue;
            if(!(syst->fNuisanceParameter==syst1->fName && syst1->fNuisanceParameter==syst1->fName)) continue;
            // now merge all SystematicHist in all regions
            LOG(DEBUG) << "Found NP(syst) " << syst->fNuisanceParameter << "(" << syst->fName << ") = to syst name " << syst1->fName << "\n";
            for(const auto& reg : fRegions){
                LOG(DEBUG) << "Region: " << reg->fName << "\n";
                for(const auto& sh : reg->fSampleHists){
                    std::shared_ptr<SystematicHist> syh  = sh->GetSystematic(syst ->fName);
                    std::shared_ptr<SystematicHist> syh1 = sh->GetSystematic(syst1->fName);
                    if(syh==nullptr || syh1 ==nullptr) continue;
                    // FIXME...
                    // the issue here is that to combine uncertainties one has to act differently depending on the fact that the different sources come from a multiplication/division or not...
                    syh1 ->Add(syh.get());
                    syh1 ->Add(sh->fHist.get(),-1);
                    LOG(DEBUG) << "Adding syst of " << syh->fName << " to " << syh1->fName << "\n";
                    LOG(DEBUG) << "Setting to 0 all Up/Down of " << syh->fName << "\n";
                    //
                    // set to zero the other syst
                    syh->fHistUp.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistUp->GetName())));
                    syh->fHistDown.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistDown->GetName())));
                    syh->fHistShapeUp.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistShapeUp->GetName())));
                    syh->fHistShapeDown.reset(static_cast<TH1*>(sh->fHist->Clone(syh->fHistShapeDown->GetName())));
                    syh->fNormUp   = 0.;
                    syh->fNormDown = 0.;
                    syh->fNormPruned  = true;
                    syh->fShapePruned = true;
                }
            }
        }
    }
}

//____________________________________________________________________________________
// this will combine special systematics into a single systematic (e.g. envelope)
void TRExFit::CombineSpecialSystematics() {
    LOG(INFO) << "Combining special systematics\n";
    std::vector<std::string> combineNames;

    // First we need to know how many special systematics there are
    for (const auto& syst : fSystematics) {
        if (syst->fCombineName == "") continue;
        if (std::find(combineNames.begin(), combineNames.end(), syst->fCombineName) == combineNames.end()) {
            combineNames.emplace_back(syst->fCombineName);
        }
    }

    if (combineNames.size() == 0) return;

    for (const auto& icombine : combineNames) {
        Systematic::COMBINATIONTYPE type(Systematic::COMBINATIONTYPE::ENVELOPE);
        // loop over regions
        std::vector<std::string> dropNorm;
        for (std::size_t iReg = 0; iReg < fRegions.size(); ++iReg) {
            // loop over systematics and get the list of systematics to combine
            std::vector<std::string> names;
            for (const auto& isyst : fSystematics) {
                if (isyst->fCombineName != icombine) continue;
                if (std::find(names.begin(), names.end(), isyst->fName) != names.end()) continue;
                names.emplace_back(isyst->fName);
                type = isyst->fCombineType;

                if (dropNorm.size() == 0) {
                    if (isyst->fDropNormSpecialIn.size() != 0) {
                        dropNorm = isyst->fDropNormSpecialIn;
                    }
                }
            }

            const bool dropFromAll = std::find(dropNorm.begin(), dropNorm.end(), "all") != dropNorm.end() ? true : false;

            // now loop over the systematics to combine
            std::vector<std::vector<std::shared_ptr<SystematicHist> > > sh_vec;

            for (const auto& sh : fRegions[iReg]->fSampleHists) {
                std::vector<std::shared_ptr<SystematicHist> > tmp;
                for (const auto& ispecial : names) {
                    auto sampleHist = sh->GetSystematic(ispecial);
                    if (!sampleHist) continue;
                    tmp.emplace_back(sh->GetSystematic(ispecial));
                }
                sh_vec.emplace_back(tmp);
            }

            if (sh_vec.size() == 0) {
                LOG(WARNING) << "No systematics to combine for " << icombine << " as no samples are affected\n";
                return;
            }

            // Now we have a list of systematicHist to combine and we change the first one, and then remove the rest
            for (std::size_t ish = 0; ish < fRegions[iReg]->fSampleHists.size(); ++ish) {
                auto newSystHist = fRegions[iReg]->fSampleHists.at(ish)->GetSystematic(names.at(0));
                if (!newSystHist) continue;
                newSystHist = CombineSpecialHistos(newSystHist, sh_vec.at(ish), type, fRegions[iReg]->fSampleHists.at(ish).get());

                // apply the drop normalisation
                if (dropFromAll || std::find(dropNorm.begin(), dropNorm.end(), fRegions[iReg]->fName) != dropNorm.end()) {
                    const SampleHist* sh = fRegions[iReg]->fSampleHists.at(ish).get();
                    if (!sh) continue;
                    const TH1* nominal = sh->fHist.get();
                    if (!nominal) continue;

                    newSystHist->fHistUp  ->Scale(nominal->Integral()/newSystHist->fHistUp  ->Integral());
                    newSystHist->fHistDown->Scale(nominal->Integral()/newSystHist->fHistDown->Integral());
                    newSystHist->fNormPruned = true;
                    newSystHist->fShapePruned = false;
                    newSystHist->fNormUp = 0.;
                    newSystHist->fNormDown = 0.;
                }
            }

            // And set the remaining systematicHist to zero
            for (std::size_t ispecial = 1; ispecial < names.size(); ++ispecial) {
                for (const auto& sh : fRegions[iReg]->fSampleHists) {
                    auto sampleHist = sh->GetSystematic(names.at(ispecial));
                    if (!sampleHist) continue;
                    sampleHist->fHistUp.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistUp->GetName())));
                    sampleHist->fHistDown.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistDown->GetName())));
                    sampleHist->fHistShapeUp.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistShapeUp->GetName())));
                    sampleHist->fHistShapeDown.reset(static_cast<TH1*>(sh->fHist->Clone(sampleHist->fHistShapeDown->GetName())));
                    sampleHist->fNormUp   = 0.;
                    sampleHist->fNormDown = 0.;
                    sampleHist->fNormPruned  = true;
                    sampleHist->fShapePruned = true;
                }
            }
        }
    }
}

//____________________________________________________________________________________
//
void TRExFit::ComputeBinning(int regIter){
    //
    //Creating histograms to rebin
    std::unique_ptr<TH1D> hsig(nullptr);
    std::unique_ptr<TH1D> hbkg(nullptr);
    bool nDefSig=true;
    bool nDefBkg=true;
    std::string fullSelection;
    std::string fullMCweight;
    std::vector<std::string> fullPaths;
    std::vector<std::string> friendPaths;
    bool bkgReg=false;
    bool flatBkg=false;
    if(fRegions[regIter]->fRegionType==Region::CONTROL) bkgReg=true;
    if(bkgReg && fRegions[regIter]->fTransfoDzSig<1e-3) flatBkg=true;
    //
    LOG(DEBUG) << "Will compute binning with the following options:\n";

    if((fRegions[regIter]->fTransfoDzSig>1e-3 || fRegions[regIter]->fTransfoDzBkg>1e-3) ) {
        LOG(DEBUG) << "TransfoD - zSig=" << fRegions[regIter]->fTransfoDzSig << " - zBkg=" << fRegions[regIter]->fTransfoDzBkg << "\n";
    }
    if((fRegions[regIter]->fTransfoFzSig>1e-3 || fRegions[regIter]->fTransfoFzBkg>1e-3) ) {
        LOG(DEBUG) << "TransfoF - zSig=" << fRegions[regIter]->fTransfoFzSig << " - zBkg=" << fRegions[regIter]->fTransfoFzBkg << "\n";
    }
    if((fRegions[regIter]->fTransfoJpar1>1e-3 || fRegions[regIter]->fTransfoJpar2>1e-3 || fRegions[regIter]->fTransfoJpar3>1e-3) ) {
        LOG(DEBUG) << "TransfoJ - z1=" << fRegions[regIter]->fTransfoJpar1 << " - z2=" << fRegions[regIter]->fTransfoJpar2 << " - z3=" << fRegions[regIter]->fTransfoJpar3 << "\n";
    }

    if(bkgReg) LOG(DEBUG) << " - bkg reg\n";
    else LOG(DEBUG) << " - sig reg\n";

    for(const auto& isample : fSamples) {
        //
        // using NTuples
        if(fInputType==1){
            if(isample->fType == Sample::SampleType::DATA) continue;
            if(isample->fType == Sample::SampleType::GHOST) continue;
            if(isample->fType == Sample::SampleType::EFT) continue;
            if(Common::FindInStringVector(isample->fRegions,fRegions[regIter]->fName) < 0) continue;
            //
            fullSelection = FullSelection(  fRegions[regIter].get(),isample.get());
            fullMCweight  = FullWeight(     fRegions[regIter].get(),isample.get());
            fullPaths     = FullNtuplePaths(fRegions[regIter].get(),isample.get());
            friendPaths = {}; //set it to an empty vector

            if (isample->fUseFriend)
                friendPaths   = FullNtuplePaths(fRegions[regIter].get(),isample.get(), nullptr, true, false, true);
            if(friendPaths.size() == 0) {
                for (unsigned int i=0; i < fullPaths.size(); i++){
                    friendPaths.push_back("");
                }
            }
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                const int tmp_debugLevel=TRExFitter::DEBUGLEVEL;
                TRExFitter::SetDebugLevel(0);
                std::unique_ptr<TH1D> htmp(Common::HistFromNtuple( fullPaths.at(i_path),
                                           friendPaths.at(i_path),
                                           fRegions.at(regIter)->fVariable,
                                           10000,
                                           fRegions.at(regIter)->fXmin,
                                           fRegions.at(regIter)->fXmax,
                                           fullSelection,
                                           fullMCweight,
                                           fAddAliases,
                                           fDebugNev));
                TRExFitter::SetDebugLevel(tmp_debugLevel);
                //
                // Pre-processing of histograms (rebinning, lumi scaling)
                if(isample->fType != Sample::SampleType::DATA && isample->fNormalizedByTheory) htmp -> Scale(fLumi);
                //
                if(isample->fLumiScales.size()>i_path) htmp -> Scale(isample->fLumiScales[i_path]);
                else if(isample->fLumiScales.size()==1) htmp -> Scale(isample->fLumiScales[0]);
                //
                // Importing the histogram in TRExFitter
                if(isample->fType == Sample::SampleType::SIGNAL){
                    if(nDefSig){
                        hsig.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                        nDefSig=false;
                    }
                    else hsig->Add(htmp.get());
                }
                else{
                    if(bkgReg && !flatBkg){
                        bool usedInSig=false;
                        for(unsigned int i_bkgs=0; i_bkgs<fRegions[regIter]->fAutoBinBkgsInSig.size(); ++i_bkgs){
                            if(isample->fName == fRegions[regIter]->fAutoBinBkgsInSig[i_bkgs]){
                                usedInSig=true;
                                break;
                            }
                        }
                        if(usedInSig){
                            LOG(DEBUG) << "Using " << isample->fName << " as signal\n";
                            if(nDefSig){
                                hsig.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                                nDefSig=false;
                            }
                            else hsig->Add(htmp.get());
                        }
                        else{
                            if(nDefBkg){
                                hbkg.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                                nDefBkg=false;
                            }
                            else hbkg->Add(htmp.get());
                        }
                    }
                    else{
                        if(nDefBkg){
                            hbkg.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                            nDefBkg=false;
                        }
                        else hbkg->Add(htmp.get());
                    }
                }
            }
        }
        //
        // Input with hists
        else if(fInputType == 0){
            if(isample->fType==Sample::SampleType::DATA) continue;
            if(isample->fType==Sample::SampleType::GHOST) continue;
            if(isample->fType==Sample::SampleType::EFT) continue;
            if(Common::FindInStringVector(isample->fRegions,fRegions[regIter]->fName) < 0) continue;
            //
            fullPaths     = FullHistogramPaths(fRegions[regIter].get(),isample.get());
            for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
                const int tmp_debugLevel=TRExFitter::DEBUGLEVEL;
                TRExFitter::SetDebugLevel(0);
                std::unique_ptr<TH1> htmp = Common::HistFromFile( fullPaths[i_path] );
                if (!htmp) {
                    if (!TRExFitter::HISTOCHECKCRASH) continue;
                    LOG(ERROR) << "Histo pointer is empty cannot continue running the code\n";
                    exit(EXIT_FAILURE);
                }
                TRExFitter::SetDebugLevel(tmp_debugLevel);
                //
                // Pre-processing of histograms (rebinning, lumi scaling)
                if(fRegions[regIter]->fHistoBins.size() > 0){
                    const std::string hname = htmp->GetName();
                    std::unique_ptr<TH1> tmp_copy(static_cast<TH1*>(htmp->Rebin(fRegions[regIter]->fHistoNBinsRebin, "tmp_copy", &fRegions[regIter]->fHistoBins[0])));
                    htmp.reset(tmp_copy.release());
                    htmp->SetName(hname.c_str());
                    if(TRExFitter::MERGEUNDEROVERFLOW) Common::MergeUnderOverFlow(htmp.get());
                }
                else if(fRegions[regIter]->fHistoNBinsRebin != -1) {
                    htmp->Rebin(fRegions[regIter]->fHistoNBinsRebin);
                }
                //
                if(isample->fType!=Sample::SampleType::DATA && isample->fNormalizedByTheory) htmp -> Scale(fLumi);
                //
                if(isample->fLumiScales.size()>i_path) htmp -> Scale(isample->fLumiScales[i_path]);
                else if(isample->fLumiScales.size()==1) htmp -> Scale(isample->fLumiScales[0]);
                //
                // apply histogram to signal or background
                if(isample->fType==Sample::SampleType::SIGNAL){
                    if(nDefSig){
                        hsig.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                        nDefSig=false;
                    }
                    else hsig->Add(htmp.get());
                }
                else{
                if(bkgReg && !flatBkg){
                    bool usedInSig=false;
                    for(unsigned int i_bkgs=0; i_bkgs<fRegions[regIter]->fAutoBinBkgsInSig.size(); ++i_bkgs){
                        if(isample->fName==fRegions[regIter]->fAutoBinBkgsInSig[i_bkgs]){
                            usedInSig=true;
                            break;
                        }
                    }
                    if(usedInSig){
                        LOG(DEBUG) << "Using " << isample->fName << " as signal\n";
                        if(nDefSig){
                            hsig.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                            nDefSig=false;
                        }
                        else hsig->Add(htmp.get());
                    }
                    else{
                        if(nDefBkg){
                            hbkg.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                            nDefBkg=false;
                        }
                        else hbkg->Add(htmp.get());
                    }
                }
                    else{
                        if(nDefBkg){
                            hbkg.reset(static_cast<TH1D*>(htmp->Clone(Form("h_%s_%s",fRegions[regIter]->fName.c_str(),isample->fName.c_str()))));
                            nDefBkg=false;
                        }
                        else hbkg->Add(htmp.get());
                    }
                }
            }
        }
    }
    //
    //computing new bins
    //
    // get a vector of bins where to rebin to get an uncertainty <= 5% per bin.
    // starting from highest bin!
    // the numbers give the lowest bin included in the new bin
    // overflowbin+1 and underflow bins are returned as the first and last element in the vector, respectively.
    std::vector<int> bins_vec;
    //
    if (!hbkg || !hsig) {
        LOG(ERROR) << "Please provide signal and background histograms!\n";
        gSystem -> Exit(1);
    }
    int nBins_Vec = hbkg -> GetNbinsX();
    int iBin = nBins_Vec; // skip overflow bin
    bins_vec.push_back(nBins_Vec + 1);
    double nBkg = hbkg -> Integral(1, nBins_Vec );
    double nSig = hsig -> Integral(1, nBins_Vec );
    bool jBreak = false;
    double jTarget = 1e5;
    double jGoal = fRegions[regIter]->fTransfoJpar1;
    //
    while (iBin > 0) {
        double sumBkg = 0;
        double sumSig = 0;
        double sumSigL = 0;
        double err2Bkg = 0;
        bool pass = false;
        double dist = 1e10;
        double distPrev = 1e10;
        //
        while (!pass && iBin > 0) {
            double nBkgBin = hbkg -> GetBinContent(iBin);
            double nSigBin = hsig -> GetBinContent(iBin);
            sumBkg += nBkgBin;
            sumSig += nSigBin;
            if (nBkgBin > 0 && nSigBin > 0) {
                sumSigL += nSigBin * std::log1p(nSigBin / nBkgBin);
            }
            err2Bkg += (hbkg -> GetBinError(iBin))*(hbkg -> GetBinError(iBin));
            //
            double err2RelBkg = 1;
            if (sumBkg != 0) {
                err2RelBkg = err2Bkg / (sumBkg*sumBkg);
            }
            //
            double err2Rel = 1.;
            if(fRegions[regIter]->fBinTransfo == "TransfoD"){
                // "trafo D"
                if (sumBkg != 0 && sumSig != 0)
                  err2Rel = 1. / (sumBkg / (nBkg / fRegions[regIter]->fTransfoDzBkg) + sumSig / (nSig / fRegions[regIter]->fTransfoDzSig));
                else if (sumBkg != 0)
                  err2Rel = (nBkg / fRegions[regIter]->fTransfoDzBkg) / sumBkg;
                else if (sumSig != 0)
                  err2Rel = (nSig / fRegions[regIter]->fTransfoDzSig) / sumSig;

                pass = std::sqrt(err2Rel) < 1 && std::sqrt(err2RelBkg) < fRegions[regIter]->fTransfoErr ;

                if ( err2RelBkg > 1 && pass) LOG(WARNING) << "There was a relative error larger than 1 in region " << fRegions[regIter]->fName << ".\n";
                // distance
                dist = std::abs(err2Rel - 1);
            }
            else if(fRegions[regIter]->fBinTransfo == "TransfoF"){
                // "trafo F" with 5% bkg stat unc
                if (sumBkg != 0 && sumSigL != 0)
                  err2Rel = 1 / (std::sqrt(sumBkg / (nBkg / fRegions[regIter]->fTransfoFzBkg)) + std::sqrt(sumSigL / (1 / fRegions[regIter]->fTransfoFzSig)));
                else if (sumBkg != 0)
                  err2Rel = std::sqrt((nBkg / fRegions[regIter]->fTransfoFzBkg) / sumBkg);
                else if (sumSigL != 0)
                  err2Rel = std::sqrt((1 / fRegions[regIter]->fTransfoFzSig) / sumSigL);
                pass = std::sqrt(err2Rel) < 1 && std::sqrt(err2RelBkg) < fRegions[regIter]->fTransfoErr ;
            }
            else if(fRegions[regIter]->fBinTransfo == "TransfoJ"){
                if (!jBreak) pass = (sumBkg >  jGoal);
                else pass = (sumBkg > jTarget);
                if( pass && !jBreak ){
                    if( (sumSig/sumBkg) <  fRegions[regIter]->fTransfoJpar2*(nSig/nBkg) ){
                        jBreak = true;
                        jTarget = hbkg->Integral(0,iBin)/ fRegions[regIter]->fTransfoJpar3;
                    }
                    else{
                        jGoal = jGoal+1;
                    }
                }
            }
            else{
                LOG(ERROR) << "Transformation method '" << fRegions[regIter]->fBinTransfo << "' unknown, try again!\n";
                exit(EXIT_FAILURE);
            }
            if (!(pass && dist > distPrev)) {
                iBin--;
            } // else use previous bin
            distPrev = dist;
        }
        // remove last bin
        if (iBin == 0 && bins_vec.size() > 1) {
            if (fRegions[regIter]->fBinTransfo == "TransfoF") {
                bins_vec.pop_back();
            }
            else if (fRegions[regIter]->fBinTransfo == "TransfoD" && bins_vec.size() > fRegions[regIter]->fTransfoDzSig + fRegions[regIter]->fTransfoDzBkg + 0.01) {
                // remove last bin if Nbin > Zsig + Zbkg
                // (1% threshold to capture rounding issues)
                bins_vec.pop_back();
            }
        }
        bins_vec.push_back(iBin + 1);
    }
    //
    //transform bin numbers in histo edges
    int nBins = bins_vec.size();
    std::vector<double> bins(nBins);
    bins[0] = hbkg->GetBinLowEdge(1);
    for(std::size_t i = 1; i < bins_vec.size()-1; ++i){
        bins[i] = hbkg->GetBinLowEdge(bins_vec[nBins-i-1]);
    }
    bins[nBins-1] = hbkg->GetBinLowEdge(hbkg->GetNbinsX() + 1);
    LOG(INFO) << "Your final binning from automatic binning function is:\n";
    std::string temp_string = "";
    for(std::size_t ibins = 0; ibins < bins_vec.size(); ++ibins) {
      temp_string+= std::to_string(bins[ibins]) + " - ";
    }
    LOG(INFO) << "  " << temp_string << "\n";
    //
    fRegions[regIter]->SetBinning(nBins-1, bins);
}

//__________________________________________________________________________________
//
void TRExFit::GetLikelihoodScan(RooWorkspace *ws,
                                const std::string& varName,
                                RooDataSet* data,
                                const std::string& toysExtension,
                                const bool plot) const {

    LOG(INFO) << "Running likelihood scan for the parameter = " << varName << "\n";

    const RooArgSet* externalConstraints = FitUtils::GetExternalConstraints(ws, fNormFactors, fRegularizationType);

    LikelihoodScanManager manager{};
    manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps, fLHscanStep);
    manager.SetNCPU(fCPU);
    manager.SetOffSet(true);
    manager.SetBlindedParameters(fBlindedParameters);
    manager.SetUseNll(fUseNllInLHscan);
    manager.SetExternalConstraints(externalConstraints);
    manager.SetToleranceScale(fToleranceScale);
    manager.SetUseAutoDiff(fUseAutoDiff);
    manager.SetNLLOffset(fNLLOffset);

    FitUtils::FittingOptions options;
    options.noGammas = fStatOnly && fGammasInStatOnly;
    options.noSystematics = fStatOnly;
    options.randomize = fUseRnd;
    options.randomValue = fRndRange;
    std::vector<std::string> constNPs;
    std::vector<double> constNPvalues;
    for (const auto& iparam : fFitFixedNPs) {
        constNPs.emplace_back(iparam.first);
        constNPvalues.emplace_back(iparam.second);
    }
    options.constNPs = std::move(constNPs);
    options.constNPvalues = std::move(constNPvalues);
    std::vector<std::string> npNames;
    std::vector<double> npValues;
    for(const auto& inf : fNormFactors) {
        if (Common::FindInStringVector(fPOIs,inf->fName)>=0) {
            continue;
        }
        npNames. emplace_back(inf->fName);
        npValues.emplace_back(inf->GetNominal());
    }
    options.initialNPs = std::move(npNames);
    options.initialNPvalues = std::move(npValues);
    options.constPOI = fFitType == TRExFit::FitType::BONLY;
    std::vector<std::pair<std::string, double> > poiVals;
    std::vector<std::string> constNFs;
    for (const auto& inf : fNormFactors) {
        poiVals.emplace_back(std::make_pair(inf->fName, inf->GetNominal()));
        if (inf->fConst) {
            constNFs.emplace_back(inf->fName);
        }
    }
    options.poiVals = std::move(poiVals);
    options.constNFs = std::move(constNFs);
    manager.SetFittingOptions(options);

    auto scanResult = manager.Run1DScan(ws, varName, data);

    const std::vector<double> x = scanResult.first;
    const std::vector<double> y = scanResult.second;

    if (x.empty() || y.empty() ) {
        LOG(WARNING) << "Skipping LHscan\n";
        return;
    }

    gSystem->mkdir((fName+"/LHoodPlots").c_str());

    std::string scanStep="";
    if (fLHscanStep != -1) scanStep = "Step" + std::to_string(fLHscanStep);

    if (toysExtension.empty()) {
        YamlConverter converter{};
        converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
        converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
        converter.WriteLikelihoodScan(scanResult, fName+"/LHoodPlots/NLLscan_"+varName+fSuffix+scanStep+".yaml");
        if (fHEPDataFormat) {
            converter.WriteLikelihoodScanHEPData(scanResult, fName, varName+fSuffix);
        }
    }

    TCanvas can("NLLscan");

    TGraph graph(x.size(), x.data(), y.data());
    graph.Draw("ALP");
    graph.GetXaxis()->SetRangeUser(x.at(0),x.back());

    // y axis
    graph.GetYaxis()->SetTitle("-#Delta #kern[-0.1]{ln(#it{L})}");
    if(TRExFitter::SYSTMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::SYSTMAP[varName].c_str());
    else if(TRExFitter::NPMAP[varName]!="") graph.GetXaxis()->SetTitle(TRExFitter::NPMAP[varName].c_str());

    TString cname="";
    cname.Append("NLLscan_");
    cname.Append(varName);

    can.SetTitle(cname);
    can.SetName(cname);
    can.cd();

    TLatex tex{};
    tex.SetTextColor(kGray+2);

    TLine l1s(x.at(0),0.5,x.back(),0.5);
    l1s.SetLineStyle(kDashed);
    l1s.SetLineColor(kGray);
    l1s.SetLineWidth(2);
    if(graph.GetMaximum()>2){
        l1s.Draw();
        tex.DrawLatex(x.back(),0.5,"#lower[-0.1]{#kern[-1]{1 #it{#sigma}   }}");
    }

    if(graph.GetMaximum()>2){
        TLine l2s(x.at(0),2,x.back(),2);
        l2s.SetLineStyle(kDashed);
        l2s.SetLineColor(kGray);
        l2s.SetLineWidth(2);
        l2s.Draw();
        tex.DrawLatex(x.back(),2,"#lower[-0.1]{#kern[-1]{2 #it{#sigma}   }}");
    }
    //
    if(graph.GetMaximum()>4.5){
        TLine l3s(x.at(0),4.5,x.back(),4.5);
        l3s.SetLineStyle(kDashed);
        l3s.SetLineColor(kGray);
        l3s.SetLineWidth(2);
        l3s.Draw();
        tex.DrawLatex(x.back(),4.5,"#lower[-0.1]{#kern[-1]{3 #it{#sigma}   }}");
    }
    //
    TLine lv0(0,graph.GetMinimum(),0,graph.GetMaximum());
    lv0.Draw();
    //
    TLine lh0(x.at(0),0,x.back(),0);
    lh0.Draw();
    can.RedrawAxis();

    if (plot) {
        Common::SaveCanvasAs(can, fName+"/LHoodPlots/NLLscan_"+varName + toysExtension +fSuffix);
    }

    // write it to a ROOT file as well
    std::unique_ptr<TFile> f(TFile::Open(fName+"/LHoodPlots/NLLscan_"+varName+fSuffix+TString(scanStep)+ (toysExtension.empty() ? "" : "_toys")+"_curve.root","UPDATE"));
    if (!f) {
        LOG(WARNING) << "Cannot open ROOT file for likelihood scan!\n";
        return;
    }
    f->cd();
    graph.Write(("LHscan"+toysExtension).c_str(),TObject::kOverwrite);
    f->Close();
}

//____________________________________________________________________________________
//
void TRExFit::Get2DLikelihoodScan( RooWorkspace *ws, const std::vector<std::string>& varNames, RooDataSet* data) const{
    if (varNames.size() != 2){
        LOG(ERROR) << "Wrong number of parameters provided for 2D likelihood scan, returning\n";
        return;
    }

    const RooArgSet* externalConstraints = FitUtils::GetExternalConstraints(ws, fNormFactors, fRegularizationType);

    LikelihoodScanManager manager{};
    manager.SetScanParamsX(fLHscanMin, fLHscanMax, fLHscanSteps, fLHscanStep);
    manager.SetScanParamsY(fLHscanMinY, fLHscanMaxY, fLHscanStepsY, fLHscanStepY);
    manager.SetNCPU(fCPU);
    manager.SetOffSet( (fLHscanStep == -1 && fLHscanStepY == -1) );
    manager.SetBlindedParameters(fBlindedParameters);
    manager.SetUseNll(fUseNllInLHscan);
    manager.SetExternalConstraints(externalConstraints);
    manager.SetToleranceScale(fToleranceScale);
    manager.SetUseAutoDiff(fUseAutoDiff);
    manager.SetNLLOffset(fNLLOffset);

    FitUtils::FittingOptions options;
    options.noGammas = fStatOnly && fGammasInStatOnly;
    options.noSystematics = fStatOnly;
    options.randomize = fUseRnd;
    options.randomValue = fRndRange;
    std::vector<std::string> constNPs;
    std::vector<double> constNPvalues;
    for (const auto& iparam : fFitFixedNPs) {
        constNPs.emplace_back(iparam.first);
        constNPvalues.emplace_back(iparam.second);
    }
    options.constNPs = std::move(constNPs);
    options.constNPvalues = std::move(constNPvalues);
    std::vector<std::string> npNames;
    std::vector<double> npValues;
    for(const auto& inf : fNormFactors) {
        if (Common::FindInStringVector(fPOIs,inf->fName)>=0) {
            continue;
        }
        npNames. emplace_back(inf->fName);
        npValues.emplace_back(inf->GetNominal());
    }
    options.initialNPs = std::move(npNames);
    options.initialNPvalues = std::move(npValues);
    options.constPOI = fFitType == TRExFit::FitType::BONLY;
    std::vector<std::pair<std::string, double> > poiVals;
    std::vector<std::string> constNFs;
    for (const auto& inf : fNormFactors) {
        poiVals.emplace_back(std::make_pair(inf->fName, inf->GetNominal()));
        if (inf->fConst) {
            constNFs.emplace_back(inf->fName);
        }
    }
    options.poiVals = std::move(poiVals);
    options.constNFs = std::move(constNFs);
    manager.SetFittingOptions(options);

    const std::pair<std::string, std::string> names = std::make_pair(varNames.at(0), varNames.at(1));
    auto scanResult = manager.Run2DScan(ws, names, data);

    const std::vector<double> x = scanResult.x;
    const std::vector<double> y = scanResult.y;
    std::vector<std::vector<double> > z = scanResult.z;

    if (x.empty() || y.empty() || z.empty()) {
        LOG(WARNING) << "Skipping LHscan\n";
        return;
    }

    // find minimium in z
    auto findMin = [](const std::vector<std::vector<double> >& vec) {
        double min(9999999);
        for (const auto& i : vec) {
            for (const auto& j : i) {
                if (j < min) min = j;
            }
        }
        return min;
    };

    const double zmin = findMin(z);

    TString LHDir("LHoodPlots/");
    if (system(TString("mkdir -vp ")+fName+"/"+LHDir) != 0) {
        LOG(ERROR) << "Cannot open the LHdir folder\n";
        return;
    }

    // make plots
    TCanvas can("NLLscan_2D_");
    can.cd();

    TGraph2D graph(fLHscanSteps * fLHscanStepsY);

    TH2D h_nll("NLL", "NLL", fLHscanSteps, manager.GetMinValX(), manager.GetMaxValX(), fLHscanStepsY, manager.GetMinValY(), manager.GetMaxValY());

    // Transfer result into 2D histogram
    for (int ipoint = 0; ipoint < fLHscanSteps; ++ipoint) {
        for (int jpoint = 0; jpoint < fLHscanStepsY; ++jpoint) {
            // shift the likelihood values to zero (unless single point)
            if ( (fLHscanStep == -1 && fLHscanStepY == -1) ) {
                z[ipoint][jpoint] -= zmin;
            }
            h_nll.SetBinContent(ipoint+1, jpoint+1, z[ipoint][jpoint]);
            unsigned int i = ipoint * fLHscanStepsY + jpoint;
            graph.SetPoint(i,x[ipoint],y[jpoint],z[ipoint][jpoint]);
        }
    }

    // Create output plot only if not running parallel processing
    if ( (fLHscanStep == -1 && fLHscanStepY == -1) ) {
        gStyle->SetPalette(57); // Reset Palette to default (Pruning or Correlation matrinx changes this)
        graph.Draw("colz");
        graph.GetXaxis()->SetRangeUser(x[0],x.back());
        graph.GetYaxis()->SetRangeUser(y[0],y.back());

        // y axis
        graph.GetXaxis()->SetTitle(varNames.at(0).c_str());
        graph.GetYaxis()->SetTitle(varNames.at(1).c_str());

        // Print the canvas
        Common::SaveCanvasAs(can, std::string(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix));

        // write it to a ROOT file as well
        std::unique_ptr<TFile> f(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+"_curve.root","UPDATE"));
        f->cd();
        graph.Write(("LHscan_2D_"+varNames.at(0)+"_"+varNames.at(1)).c_str(),TObject::kOverwrite);
        f->Close();
    }

    // store yaml format
    YamlConverter converter{};
    converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
    converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
    converter.WriteLikelihoodScan2D(scanResult, fName+"/LHoodPlots/NLLscan_"+varNames.at(0)+"_" + varNames.at(1)+fSuffix+".yaml");
    if (fHEPDataFormat) {
        converter.WriteLikelihoodScan2DHEPData(scanResult, fName, varNames.at(0)+ "_" + varNames.at(1)+fSuffix);
    }

    // Write histogram to Root file as well
    TString scanStep="";
    if (fLHscanStep  != -1) {
        scanStep += "_StepX" + TString::Itoa(fLHscanStep, 10);
    }
    if (fLHscanStepY != -1) {
        scanStep += "_StepY" + TString::Itoa(fLHscanStepY, 10);
    }

    std::unique_ptr<TFile> f2(TFile::Open(fName+"/"+LHDir+"NLLscan_"+varNames.at(0)+"_"+varNames.at(1)+fSuffix+scanStep+"_histo.root","UPDATE"));
    h_nll.Write("NLL",TObject::kOverwrite);
    f2->Close();
}

//__________________________________________________________________________________
//
void TRExFit::AddTemplateWeight(const std::string& name, double value){
    std::pair<double, std::string> temp = std::make_pair(value, name);
    fTemplatePair.push_back(temp);
}

//__________________________________________________________________________________
//
std::vector<TRExFit::TemplateWeight> TRExFit::GetTemplateWeightVec(const TRExFit::TemplateInterpolationOption& opt){
    std::vector<TRExFit::TemplateWeight> vec;
    for(const auto& name : fMorphParams){
        // create map only for values of the specified parameter
        std::vector<std::pair<double,std::string> > templatePair;
        for(auto tp : fTemplatePair){
            if(tp.second==name) templatePair.push_back(tp);
        }
        // first sort vector of inputs for templates
        if (templatePair.size() < 2){
            LOG(ERROR) << "You need to provide at least 2 templates for template fit to work, but you provided: " << fTemplatePair.size() << "\n";
        }
        std::sort(templatePair.begin(), templatePair.end());
        // find min and max for range
        double min = templatePair.at(0).first;
        double max = templatePair.at(templatePair.size() -1).first;

        // find the norm factor with the name of the template and set the nominal value based onm that
        auto itr = std::find_if(fNormFactors.begin(), fNormFactors.end(), [&templatePair](const auto& nf){return nf->fName == templatePair.at(0).second;});
        if (itr == fNormFactors.end()) {
            LOG(ERROR) << "Cannot find norm factor: " << templatePair.at(0).second << "\n";
            exit(EXIT_FAILURE);
        }

        const double nominal = (*itr)->GetNominal();

        for (unsigned int itemp = 0; itemp < (templatePair.size() ); itemp++){
            LOG(DEBUG) << "Morphing: Template " << itemp << "\n";
            TRExFit::TemplateWeight tmp;
            tmp.name = templatePair.at(itemp).second;
            tmp.value = templatePair.at(itemp).first;
            LOG(DEBUG) << "Morphing:   " << tmp.name << " = " << tmp.value << "\n";
            tmp.range = tmp.name+"["+std::to_string(nominal)+","+std::to_string(min)+","+std::to_string(max)+"]";
            // calculate the actual function
            tmp.function = TRExFit::GetWeightFunction(templatePair, itemp, opt);
            vec.push_back(tmp);
        }
    }
    return vec;
}

//__________________________________________________________________________________
//
std::string TRExFit::GetWeightFunction(const std::vector<std::pair<double,std::string> >& templatePair, unsigned int itemp, const TRExFit::TemplateInterpolationOption& opt) const{
    std::string fun = "";
    double x_i;
    std::string name;
    if (itemp < templatePair.size()){
        x_i = templatePair.at(itemp).first;
        name = templatePair.at(itemp).second;
    }
    else return fun;
    //
    if (opt == TRExFit::LINEAR){
        double deltaXp(-1.); // |x(i+1)-x(i)|
        double deltaXm(-1.); // |x(i-1)-x(i)|
        if ((itemp+1) < templatePair.size() ){
            deltaXp = std::abs(templatePair.at(itemp+1).first - templatePair.at(itemp).first);
        }
        if (((int)itemp-1) >=0 ){
            deltaXm = std::abs(templatePair.at(itemp-1).first - templatePair.at(itemp).first);
        }
        if(deltaXp<0 && deltaXm<0){
            LOG(ERROR) << "Morphing: delta X = " << deltaXp << ", " << deltaXm << "\n";
            return fun;
        }
        fun  = "(";
        if(deltaXm>0) fun += "((("+name+"-"+std::to_string(x_i)+")< 0)&&(abs("+name+"-"+std::to_string(x_i)+")<"+std::to_string(deltaXm)+"))*(1.-(abs("+name+"-"+std::to_string(x_i)+"))/"+std::to_string(deltaXm)+")";
        else fun += "0.";
        fun += "+";
        if(deltaXp>0) fun += "((("+name+"-"+std::to_string(x_i)+")>=0)&&(abs("+name+"-"+std::to_string(x_i)+")<"+std::to_string(deltaXp)+"))*(1.-(abs("+name+"-"+std::to_string(x_i)+"))/"+std::to_string(deltaXp)+")";
        else fun += "0.";
        fun += ")";
    } else if (opt == TRExFit::SMOOTHLINEAR) {
        // this will return a string that represents integral of hyperbolic tangent function that
        // approximates absolute value
        fun = GetSmoothLinearInterpolation(itemp);
    } else if (opt == TRExFit::SQUAREROOT) {
        fun = GetSquareRootLinearInterpolation(itemp);
    }
    // ...
    fun = Common::ReplaceString(fun,"--","+");
    LOG(DEBUG) << "Morphing:   weight function = " << fun << "\n";
    return fun;
}

//__________________________________________________________________________________
//
std::string TRExFit::GetSmoothLinearInterpolation(unsigned int itemp) const {
    // The idea is simple: use integral of hyperbolic tangent
    //
    // -2/width*(1/k)*corr*(log(e^(k*(x-x_mean)+e^(-(k*(x-x_mean))))) - ln(2)) + 1
    //
    // but the function is split into two parts to allow also non-equal steps in
    // the templates used for the template fit
    // "width" represents the x-axis size that is used for that particular function
    // "corr" represents correction to the function so that for x_min/x_max the functional value
    // is 0, without this correction it won't exactly be zero, because it is an anproximation

    if (itemp >= fTemplatePair.size()){
        return "";
    }

    // parameter that controls how close to a linear function we want to be
    static const double k_init(80.);

    double x_left = -99999.;
    double x_right = -99999.;
    double corr_left = 1.;
    double corr_right = 1.;

    if (itemp == 0) { // first template
        x_left = 2*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp+1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    } else if (itemp == (fTemplatePair.size()-1)) { // last template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = 2*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp-1).first;
    } else { // general template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    }

    double x_mean = fTemplatePair.at(itemp).first;
    double width_left = 2*std::abs(x_mean - x_left);
    double width_right = 2*std::abs(x_mean - x_right);

    // apply correction to the k parameter depending on the range of the x axis
    const double k_left = k_init/width_left;
    const double k_right = k_init/width_right;

    // calculate correction to the function to get y= 0 at x_min and x_max
    // use iterative process to find something which is close enough
    corr_left= 1+GetCorrection(k_left, width_left, x_mean, x_left);
    corr_left+= GetCorrection(k_left, width_left, x_mean, x_left, corr_left);
    corr_left+= GetCorrection(k_left, width_left, x_mean, x_left, corr_left);
    corr_right= 1+GetCorrection(k_right, width_right, x_mean, x_right);
    corr_right+= GetCorrection(k_right, width_right, x_mean, x_right, corr_right);
    corr_right+= GetCorrection(k_right, width_right, x_mean, x_right, corr_right);

    // prepare the actual string as "function" + "step function"
    std::string name = fTemplatePair.at(itemp).second;
    std::string step_left = "";
    std::string step_right = "";

    // the step function
    if (itemp == 0) {
        step_left = "(("+name+"-"+std::to_string(x_mean)+"<0)&&("+name+"-"+std::to_string(x_mean)+">0))";
        step_right = "(("+name+"-"+std::to_string(x_mean)+">=0) && ("+name+"<"+std::to_string(x_right)+"))";
    } else if (itemp == (fTemplatePair.size()-1)) {
        step_left = "(("+name+">="+std::to_string(x_left)+")&&("+name+"<"+std::to_string(x_mean)+"))";
        step_right = "(("+name+"-"+std::to_string(x_mean)+"<0)&&("+name+">"+std::to_string(x_right)+"))";
    } else {
        step_left = "((("+name+"-"+std::to_string(x_mean)+")<=0)&&("+name+">"+std::to_string(x_left)+"))";
        step_right = "((("+name+"-"+std::to_string(x_mean)+")>0)&&("+name+"<"+std::to_string(x_right)+"))";
    }

    // the functions
    std::string fun_left = "(-2/("+std::to_string(width_left*k_left)+")*"+std::to_string(corr_left)+"*(log(exp("+std::to_string(k_left)+"*("+name+"-"+std::to_string(x_mean)+")) + exp(-"+std::to_string(k_left)+"*("+name+"-"+std::to_string(x_mean)+")"+")) - log(2)) +1) * " + step_left;

    std::string fun_right = "(-2/("+std::to_string(width_right*k_right)+")*"+std::to_string(corr_right)+"*(log(exp("+std::to_string(k_right)+"*("+name+"-"+std::to_string(x_mean)+")) + exp(-"+std::to_string(k_right)+"*("+name+"-"+std::to_string(x_mean)+")"+")) - log(2)) +1) * " + step_right;

    return ("(("+fun_left+") + (" +fun_right+"))");
}

//__________________________________________________________________________________
//
double TRExFit::GetCorrection(double k, double width, double x_mean, double x_left, double init) const {
    double logterm = 0;
    double corr = 0;

    // since we are calculating logarithm of potentially wery large numbers (e^20)
    // we need to help the code in case we get overflow, when this happens we simply discard
    // the smaller contribution and set log(e^(x)) = x manually;

    //check if the number is inf
    if (std::isinf(std::exp(k*(x_left-x_mean)))) {
        logterm = k*(x_left-x_mean);
    } else if (std::isinf(std::exp(-k*(x_left-x_mean)))) {
        logterm = -k*(x_left-x_mean);
    } else {
        logterm = std::log(exp(k*(x_left-x_mean)) + std::exp(-k*(x_left-x_mean)));
    }
    corr = ((-2/(k*width))*init*(logterm - std::log(2))+1);

    return corr;
}

//__________________________________________________________________________________
//
std::string TRExFit::GetSquareRootLinearInterpolation(unsigned int itemp) const {
    double epsilon = 0.0000001;

    double x_i = fTemplatePair.at(itemp).first;
    double x_left = -99999.;
    double x_right = -99999.;

    if (itemp == 0) { // first template
        x_left = 2.*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp+1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    } else if (itemp == (fTemplatePair.size()-1)) { // last template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = 2.*fTemplatePair.at(itemp).first - fTemplatePair.at(itemp-1).first;
    } else { // general template
        x_left = fTemplatePair.at(itemp-1).first;
        x_right = fTemplatePair.at(itemp+1).first;
    }

    //apply correction
    double a_left = 0;
    double b_left = 0;
    double a_right = 0;
    double b_right = 0;

    GetSquareCorrection(&a_left, &b_left, x_i, x_left, epsilon);
    GetSquareCorrection(&a_right, &b_right, x_i, x_right, epsilon);

    // prepare the actual string as "function" + "step function"
    std::string name = fTemplatePair.at(itemp).second;
    std::string step_left = "";
    std::string step_right = "";

    // the step function
    if (itemp == 0) {
        step_left = "(("+name+"-"+std::to_string(x_i)+"<0)&&("+name+"-"+std::to_string(x_i)+">0))";
        step_right = "(("+name+"-"+std::to_string(x_i)+">=0) && ("+name+"<"+std::to_string(x_right)+"))";
    } else if (itemp == (fTemplatePair.size()-1)) {
        step_left = "(("+name+">="+std::to_string(x_left)+")&&("+name+"<"+std::to_string(x_i)+"))";
        step_right = "(("+name+"-"+std::to_string(x_i)+"<0)&&("+name+">"+std::to_string(x_right)+"))";
    } else {
        step_left = "((("+name+"-"+std::to_string(x_i)+")<=0)&&("+name+">"+std::to_string(x_left)+"))";
        step_right = "((("+name+"-"+std::to_string(x_i)+")>0)&&("+name+"<"+std::to_string(x_right)+"))";
    }

    std::string fun_left = "("+step_left+")*(-"+std::to_string(a_left)+"*sqrt(("+name+"-"+std::to_string(x_i)+")*("+name+"-"+std::to_string(x_i)+")+"+std::to_string(epsilon)+")+"+std::to_string(b_left)+")";
    std::string fun_right = "("+step_right+")*(-"+std::to_string(a_right)+"*sqrt(("+name+"-"+std::to_string(x_i)+")*("+name+"-"+std::to_string(x_i)+")+"+std::to_string(epsilon)+")+"+std::to_string(b_right)+")";

    return ("("+fun_left+"+"+fun_right+")");
}

//__________________________________________________________________________________
//
void TRExFit::GetSquareCorrection(double *a, double *b, double x_i, double x_left, double epsilon) const {
    if (x_left == 0) {
        x_left = 2*x_i;
    }

    // this can be analytically calculated
    double k = std::sqrt(((x_i - x_left)*(x_i - x_left)/epsilon) + 1);
    *b = k/(k-1);
    *a = (*b ) / std::sqrt((x_i-x_left)*(x_i-x_left) + epsilon);
}

//__________________________________________________________________________________
//
void TRExFit::SmoothMorphTemplates(const std::string& name,const std::string& formula,double *p) const{
    TCanvas c("c","c",600,600);
    // find NF associated to this morph param
    std::shared_ptr<NormFactor> nf = nullptr;
    for(const auto& norm : fNormFactors){
        if(norm->fName == name) nf = norm;
    }
    // get one histogram per bin (per region)
    for(const auto& reg : fRegions){
        std::map<double,TH1*> hMap; // map (paramater-value,histogram)
        TH1* h_tmp = nullptr;
        int nTemplates = 0;
        double min = -999.;
        double max = -999.;
        for(const auto& sh : reg->fSampleHists){
            Sample* smp = sh->fSample;
            // if the sample has morphing
            if(smp->fIsMorph[name]){
                hMap[smp->fMorphValue[name]] = sh->fHist.get();
                if(smp->fMorphValue[name]<min) min = smp->fMorphValue[name];
                if(smp->fMorphValue[name]>max) max = smp->fMorphValue[name];
                nTemplates++;
                h_tmp = sh->fHist.get();
            }
        }
        if(h_tmp==nullptr) return;
        for(int i_bin=1;i_bin<=h_tmp->GetNbinsX();i_bin++){
            TGraphErrors g_bin(nTemplates);
            int i_pt = 0;
            for(auto vh : hMap){
                g_bin.SetPoint(i_pt,vh.first,vh.second->GetBinContent(i_bin));
                g_bin.SetPointError(i_pt,0,vh.second->GetBinError(i_bin));
                // if it's the nominal sample, set error to very small value => forced not to change nominal!
                if(nf!=nullptr) if(nf->GetNominal()==vh.first) g_bin.SetPointError(i_pt,0,vh.second->GetBinError(i_bin)*0.001);
                i_pt++;
            }
            c.cd();
            g_bin.Draw("epa");
            TF1 l("l",formula.c_str(),min,max);
            if(p!=0x0) l.SetParameters(p);
            g_bin.Fit("l","RQN");
            l.SetLineColor(kRed);
            l.Draw("same");
            gSystem->mkdir((fName+"/Morphing/").c_str());
            Common::SaveCanvasAs(c, fName+"/Morphing/g_"+name+"_"+reg->fName+"_bin"+std::to_string(i_bin));
            for(auto vh : hMap){
                vh.second->SetBinContent(i_bin,l.Eval(vh.first));
            }
        }
    }
}

//____________________________________________________________________________________
//
bool TRExFit::MorphIsAlreadyPresent(const std::string& name, const double value) const {
    for (const std::pair<double, std::string>& itemp : fTemplatePair){
        if ((itemp.second == name) && (itemp.first == value)){
            return true;
        }
    }
    return false;
}

//____________________________________________________________________________________
// create a map associating parameters to their SubCategory
void TRExFit::ProduceSystSubCategoryMap(){
    LOG(DEBUG) << "Filling SubCategory map\n";

    fSubCategoryImpactMap.clear();

    // special treatment needed for two cases:
    // 1) stat-only fit where all parameters are fixed, see FittingTool::GetGroupedImpact()
    // 2) fit with all Gammas fixed, see FittingTool::GetGroupedImpact()
    fSubCategoryImpactMap.insert(std::make_pair("DUMMY_STATONLY", "FullSyst"));
    fSubCategoryImpactMap.insert(std::make_pair("DUMMY_GAMMAS", "Gammas"));

    // add all systematics, here an "alpha_" prefix is needed
    for(const auto& isyst : fSystematics) {
        if(isyst->fSubCategory=="Gammas" || isyst->fSubCategory=="FullSyst" || isyst->fSubCategory=="combine")
             LOG(WARNING) << "Use of \"Gammas\", \"FullSyst\" or \"combine\" as SubCategory names is not supported, you will likely run into issues\n";
        if(isyst->fType!=Systematic::SHAPE){
            fSubCategoryImpactMap.insert(std::make_pair(("alpha_" + isyst->fNuisanceParameter).c_str(), isyst->fSubCategory));
        }
        else{
            // treat SHAPE systematics separately, since they are not prefixed with "alpha_", but "gamma_shape_" instead
            // need one per bin per region
            for(const auto& reg : fRegions){
                if(reg->GetNbins() == 0){
                    LOG(ERROR) << "Cannot determine binning (no samples assigned to region?), exiting\n";
                    exit(EXIT_FAILURE);
                }

                // determine the amount of bins in this region, requires that ReadHistos() was used
                const int nRegionBins = reg->GetNbins();

                for(int i_bin=1; i_bin < nRegionBins+1; i_bin++){
                    fSubCategoryImpactMap.insert(std::make_pair(Form("gamma_shape_%s_%s_bin_%d",(isyst->fNuisanceParameter).c_str(),(reg->fName).c_str(),i_bin-1), isyst->fSubCategory));
                }
            }
        }
    }

    // need to add separate gammas
    for (const auto& isample : fSamples) {
        if (!isample->fSeparateGammas) continue;
        const std::string systName = "gamma_"+isample->fName;
        if (isample->fRegions.at(0) == "all") {
            for (const auto& ireg : fRegions) {
                if (ireg->fRegionType == Region::RegionType::VALIDATION) continue;
                for(int ibin = 1; ibin < ireg->GetNbins()+1; ++ibin) {
                    const std::string paramName = Form("gamma_shape_stat_%s_%s_bin_%d",isample->fName.c_str(),(ireg->fName).c_str(),ibin-1);
                    fSubCategoryImpactMap.insert(std::make_pair(paramName, systName));
                }
            }
        } else {
            for (const auto& reg : isample->fRegions) {
                auto ireg = std::find_if(fRegions.begin(), fRegions.end(), [&reg](const auto& element){return reg == element->fName;});
                if (ireg == fRegions.end()) {
                    LOG(ERROR) << "Cannot find region " << reg << "\n";
                    exit(EXIT_FAILURE);
                }
                if ((*ireg)->fRegionType == Region::RegionType::VALIDATION) continue;
                for(int ibin = 1; ibin < (*ireg)->GetNbins()+1; ++ibin) {
                    const std::string paramName = Form("gamma_shape_stat_%s_%s_bin_%d",isample->fName.c_str(),((*ireg)->fName).c_str(),ibin-1);
                    fSubCategoryImpactMap.insert(std::make_pair(paramName, systName));
                }
            }
        }
    }

    // also add norm factors, no "alpha_" needed
    for(const auto& inorm : fNormFactors) {
        if(inorm->fSubCategory=="Gammas" || inorm->fSubCategory=="FullSyst" || inorm->fSubCategory=="combine")
             LOG(WARNING) << "Use of \"Gammas\", \"FullSyst\" or \"combine\" as SubCategory names is not supported, you will likely run into issues\n";
        if (Common::FindInStringVector(fPOIs,inorm->fName)<0) {
            fSubCategoryImpactMap.insert(std::make_pair(inorm->fNuisanceParameter, inorm->fSubCategory));
        }
    }
}

//____________________________________________________________________________________
// combine individual results from grouped impact evaluation into one table
void TRExFit::BuildGroupedImpactTable() {
    LOG(INFO) << "Merging grouped impact evaluations\n";
    for (const auto& ipoi : fPOIs) {
        const std::string targetName = fName+"/Fits/GroupedImpact_"+ipoi+fSuffix+".txt";

        if(std::ifstream(targetName).good()){
            LOG(WARNING) << "File " << targetName << " already exists, will not overwrite\n";
        } else {
            const std::vector<std::string>& inPaths = Common::GetFilesMatchingString(fName+"/Fits/","GroupedImpact_", ipoi+fSuffix+".txt");
            Common::MergeTxTFiles(inPaths, targetName);
        }
    }
}

//____________________________________________________________________________________
//
void TRExFit::RunToys() {
    gSystem->mkdir( (fName+"/Toys").c_str());
    LOG(INFO) << "\n";
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Generating and fitting toys...\n";
    LOG(INFO) << "-------------------------------------------\n";

    if (fCPU > 1) {
        LOG(WARNING) << "Toys are not supported in the multi-processor mode, setting fCPU = 1\n";
        fCPU = 1;
    }

    std::vector < std:: string > regionsToFit;
    for(const auto& ireg : fRegions) {
        if (fFitRegion == CRONLY && ireg->fRegionType == Region::CONTROL) {
            regionsToFit.emplace_back(ireg->fName);
        } else if (fFitRegion == CRSR && (ireg->fRegionType == Region::CONTROL || ireg->fRegionType == Region::SIGNAL)) {
            regionsToFit.emplace_back(ireg->fName);
        }
    }
    std::unique_ptr<RooWorkspace> ws = PerformWorkspaceCombinationxRooFit(regionsToFit);
    if (!ws){
        LOG(ERROR) << "Cannot retrieve the workspace, exiting!\n";
        exit(EXIT_FAILURE);
    }
    if (fBinnedLikelihood) {
        FitUtils::SetBinnedLikelihoodOptimisation(ws.get());
    }
    if (fIntCode != 4) {
        FitUtils::ChangeInterpolationCode(ws.get(), fIntCode);
    }

    bool isRegularisedUnfolding(false);
    if (fFitType == TRExFit::FitType::UNFOLDING) {
        for (const auto& inf : fNormFactors) {
            if (inf->fTau > 0) {
                isRegularisedUnfolding = true;
                break;
            }
        }
    }

    FitToys fitToys(fName, fSuffix, fInputName);
    fitToys.SetMinos(fVarNameMinos);
    fitToys.SetIsRegularizedUnfolding(isRegularisedUnfolding);
    fitToys.SetNormFactors(fNormFactors);
    fitToys.SetNToys(fFitToys);
    fitToys.SetRegularizationType(fRegularizationType);
    fitToys.SetUseAutoDiff(fUseAutoDiff);
    fitToys.SetFitType(fFitType);
    fitToys.SetHistoNbins(fToysHistoNbins);
    fitToys.SetPOIs(fPOIs);
    fitToys.SetFitIsBlind(fFitIsBlind);
    fitToys.SetPOIAsimov(fFitPOIAsimov);
    fitToys.SetSeed(fToysSeed);
    fitToys.SetStatOnlyFluctuation(fToysOnlyStatFluctuation);
    fitToys.SetRandomPOIStartingValue(fToyRandomPOIStartingValues);
    fitToys.SetPlotLabel(fPlotLabel);
    fitToys.SetRunLHScanForAll(fToysLHScanForAll);
    fitToys.SetLHScanCondition(fToysLHScanCondition);
    fitToys.SetFitFixedNPs(fFitFixedNPs);

    if (fToysNpValuesFile != "") {
        LOG(INFO) << "Inject the NP values from " << fToysNpValuesFile << " file and shifting the global observables to these values\n";
        ReadFitResults(fToysNpValuesFile);
        FitUtils::InjectGlobalObservables(ws.get(), fFitResults.get());
    }

    fitToys.RunToys(ws.get(), fFitResults.get());
}

//__________________________________________________________________________________
// Computes the variable string to be used when reading ntuples, for a given region, sample combination
std::string TRExFit::Variable(Region *reg, const Sample *smp){
    // protection against nullptr
    if(reg==nullptr){
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    //
    std::string variable = "";
    // from Region
    if(reg->UseAlternativeVariable(smp->fName)){
        variable = reg->GetAlternativeVariable(smp->fName);
    }
    else{
        variable = reg->fVariable;
    }
    // check the final expression
    if(!Common::CheckExpression(variable)){
        LOG(ERROR) << "Variable expression not valid. Please check: " << variable << "\n";
        exit(EXIT_FAILURE);
    }
    //
    return variable;
}

//__________________________________________________________________________________
// Computes the full selection string to be used when reading ntuples, for a given region, sample combination
std::string TRExFit::FullSelection(Region *reg, const Sample *smp){
    // protection against nullptr
    if(reg==nullptr){
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    //
    std::string selection = "";
    // from Job
    // smp->fGridDID=="" to ignore selections for sample delivered by ServiceX
    if(fSelection!="" && fSelection!="1" && smp->fGridDID==""){
        selection = "("+fSelection+")";
    }
    // from Region
    if(reg->UseAlternativeSelection(smp->fName)){
        if(selection!="") selection += " && ";
        selection += "("+reg->GetAlternativeSelection(smp->fName)+")";
    }
    else if(reg->fSelection!="" && reg->fSelection!="1"){
        if(selection!="") selection += " && ";
        selection += "("+reg->fSelection+")";
    }
    // eventually apply IgnoreSelection Sample-option
    if(smp->fIgnoreSelection=="TRUE"){
        selection = "";
    }
    else if(smp->fIgnoreSelection!="FALSE" && smp->fIgnoreSelection!=""){
        selection = Common::ReplaceString(selection,smp->fIgnoreSelection,"1");
    }
    // from Sample
    // smp->fGridDID=="" to ignore selections for sample delivered by ServiceX
    if(smp->fSelection!="" && smp->fSelection!="1" && smp->fGridDID==""){
        if(selection!="") selection += " && ";
        selection += "("+smp->fSelection+")";
    }
    // check the final expression
    if(!Common::CheckExpression(selection)){
        LOG(ERROR) << "Full selection expression not valid. Please check: " << selection << "\n";
        exit(EXIT_FAILURE);
    }
    //
    if(selection=="") selection = "1";
    return selection;
}

//__________________________________________________________________________________
// Computes the full weight string to be used when reading ntuples, for a given region, sample and systematic combination
std::string TRExFit::FullWeight(const Region *reg, const Sample *smp, const Systematic *syst,bool isUp){
    // protection against nullptr
    if(reg==nullptr){
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    // if it's Data, just return "1"
    if(smp->fType==Sample::SampleType::DATA) return "1";
    //
    std::string weight = "";
    // from Job (only for fNormalizedByTheory samples)
    if(fMCweight!="" && fMCweight!="1" && smp->fNormalizedByTheory){
        weight += "("+fMCweight+")";
    }
    // from Region
    if(reg->fMCweight!="" && reg->fMCweight!="1"){
        if(weight!="") weight += " * ";
        weight += "("+reg->fMCweight+")";
    }
    // eventually apply IgnoreWeight Sample-option
    if(smp->fIgnoreWeight=="TRUE"){
        weight = "";
    }
    else if(smp->fIgnoreWeight!="FALSE" && smp->fIgnoreWeight!=""){
        weight = Common::ReplaceString(weight,smp->fIgnoreWeight,"1");
    }
    // from Sample (nominal...
    std::string sampleWeight = "";
    if(smp->fMCweight!="" && smp->fMCweight!="1"){
        sampleWeight = "("+smp->fMCweight+")";
    }
    // ... and systematics)
    if(syst!=nullptr){
        if(syst->fIgnoreWeight!=""){
            weight = Common::ReplaceString(weight,syst->fIgnoreWeight,"1");
            sampleWeight = Common::ReplaceString(sampleWeight,syst->fIgnoreWeight,"1");
        }
        if(isUp){
            if(syst->fWeightUp!=""){
                sampleWeight = syst->fWeightUp;
            }
            else if(syst->fWeightSufUp!=""){
                if(sampleWeight!="") sampleWeight += " * ";
                sampleWeight += "("+syst->fWeightSufUp+")";
            }
        }
        else{
            if(syst->fWeightDown!=""){
                sampleWeight = syst->fWeightDown;
            }
            else if(syst->fWeightSufDown!=""){
                if(sampleWeight!="") sampleWeight += " * ";
                sampleWeight += "("+syst->fWeightSufDown+")";
            }
        }
    }
    if(sampleWeight!=""){
        if(weight!="") weight += " * ";
        weight += "("+sampleWeight+")";
    }
    // add Bootstrap weights
    // fBootstrapSyst means only systematic variation is bootstrapped, not nominal
    // WiP fBootstrapSample means bootstrap on sample and all the correlated systs (so no fntupleFiles -> not sure it captures everything)
    if(fBootstrap!="" && fBootstrapIdx>=0
      && !(fBootstrapSyst!="" && syst==nullptr)
      && (fBootstrapSample=="" || ( smp->fName==fBootstrapSample && ( !syst || syst->fIsCorrelated ) ) ) ){
        if(weight!="") weight += " * ";
        weight += "("+Common::ReplaceString(fBootstrap,"BootstrapIdx",Form("%d",fBootstrapIdx))+")";
        gRandom->SetSeed(fBootstrapIdx);
    }
    // check the final expression
    LOG(DEBUG) << "Full weight expression : " << weight << "\n";
    if(!Common::CheckExpression(weight)){
        LOG(ERROR) << "Full weight expression not valid. Please check: " << weight << "\n";
        exit(EXIT_FAILURE);
    }
    //
    if(weight=="") weight = "1";
    return weight;
}

//__________________________________________________________________________________
// Computes the full list of path + file-name + ntuple-name string to be used when reading ntuples, for a given region, sample and systematic combination
std::vector<std::string> TRExFit::FullNtuplePaths(Region *reg,
                                                  Sample *smp,
                                                  Systematic *syst,
                                                  bool isUp,
                                                  bool isSubtract,
                                                  bool isFriend){
    // protection against nullptr
    if(reg==nullptr){
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    bool isData = (smp->fType==Sample::SampleType::DATA);
    if(syst && !isData) {
        if (isSubtract) { // "special" systematics like JER
            if (isUp) {
                if (!isFriend) {
                    paths = syst->fNtuplePathsUpSubtractSample;
                } else {
                    paths = syst->fFriendPathsUpSubtractSample;
                }
                files = syst->fNtupleFilesUpSubtractSample;
                names = syst->fNtupleNamesUpSubtractSample;
            } else {
                if (!isFriend) {
                    paths = syst->fNtuplePathsDownSubtractSample;
                } else {
                    paths = syst->fFriendPathsDownSubtractSample;
                }
                files = syst->fNtupleFilesDownSubtractSample;
                names = syst->fNtupleNamesDownSubtractSample;
            }
        } else { // "normal" systematics
            if(isUp){
                if (!isFriend) {
                    if(syst->fNtuplePathsUp.size()  >0) paths = syst->fNtuplePathsUp;
                } else {
                    if(syst->fFriendPathsUp.size()  >0) paths = syst->fFriendPathsUp;
                }
                if(syst->fNtupleFilesUp.size()  >0) files = syst->fNtupleFilesUp;
                if(syst->fNtupleNamesUp.size()  >0) names = syst->fNtupleNamesUp;
            } else{
                if (!isFriend) {
                    if(syst->fNtuplePathsDown.size()>0) paths = syst->fNtuplePathsDown;
                } else {
                    if(syst->fFriendPathsDown.size()>0) paths = syst->fFriendPathsDown;
                }
                if(syst->fNtupleFilesDown.size()>0) files = syst->fNtupleFilesDown;
                if(syst->fNtupleNamesDown.size()>0) names = syst->fNtupleNamesDown;
            }
        }
    }
    if (!isFriend) {
        if(paths.empty() && smp->fNtuplePaths.size()>0) paths = smp->fNtuplePaths;
    } else {
        if(paths.empty() && smp->fFriendPaths.size()>0) paths = smp->fFriendPaths;
    }
    if(files.empty() && smp->fNtupleFiles.size()>0) files = smp->fNtupleFiles;
    if(names.empty() && smp->fNtupleNames.size()>0) names = smp->fNtupleNames;
    //
    if (!isFriend) {
        if(paths.empty() && reg->fNtuplePaths.size()>0) paths = reg->fNtuplePaths;
    } else {
        if(paths.empty() && reg->fFriendPaths.size()>0) paths = reg->fFriendPaths;
    }
    if(files.empty() && reg->fNtupleFiles.size()>0) files = reg->fNtupleFiles;
    if(names.empty() && reg->fNtupleNames.size()>0) names = reg->fNtupleNames;
    //
    if (!isFriend) {
        if(paths.empty() && fNtuplePaths.size()>0) paths = fNtuplePaths;
    } else {
        if(paths.empty() && fFriendPaths.size()>0) paths = fFriendPaths;
    }
    if(files.empty() && fNtupleFiles.size()>0) files = fNtupleFiles;
    if(names.empty() && fNtupleNames.size()>0) names = fNtupleNames;
    //
    // now combining suffs, with all the combinations instead of giving priority
    // (same order as above used for suffix order - OK? FIXME)
    if(syst && !isData) {
        if(isUp) {
            if (!isFriend) {
                pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), Common::ToVec(syst->fNtuplePathSufUp) );
            } else {
                pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fFriendPathSuffs, smp->fFriendPathSuffs ), Common::ToVec(syst->fFriendPathSufUp) );
            }
        } else{
            if (!isFriend) {
                pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs ), Common::ToVec(syst->fNtuplePathSufDown) );
            } else{
                pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fFriendPathSuffs, smp->fFriendPathSuffs ), Common::ToVec(syst->fFriendPathSufDown) );
            }
        }
    } else {
        if (!isFriend) {
            pathSuffs = Common::CombinePathSufs( reg->fNtuplePathSuffs, smp->fNtuplePathSuffs );
        } else{
            pathSuffs = Common::CombinePathSufs( reg->fFriendPathSuffs, smp->fFriendPathSuffs );
        }
    }
    //
    if(syst && !isData) {
        if(isSubtract) {
           if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufUpSubtractSample) );
           else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufDownSubtractSample) );
       }
       else{
           if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufUp) );
           else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs ), Common::ToVec(syst->fNtupleFileSufDown) );
       }
    } else {
        fileSuffs = Common::CombinePathSufs( reg->fNtupleFileSuffs, smp->fNtupleFileSuffs );
    }
    //
    if(syst && !isData) {
        if(isSubtract) {
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufUpSubtractSample) );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufDownSubtractSample) );
        }
        else {
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufUp) );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs ), Common::ToVec(syst->fNtupleNameSufDown) );
        }
    } else {
        nameSuffs = Common::CombinePathSufs( reg->fNtupleNameSuffs, smp->fNtupleNameSuffs );
    }
    //
    // And finally put everything together
    fullPaths = Common::CreatePathsList( paths,pathSuffs, files,fileSuffs, names,nameSuffs );
    return fullPaths;
}

//__________________________________________________________________________________
// Computes the full list of path + file-name + histogram-name string to be used when reading ntuples, for a given region, sample and systematic combination
std::vector<std::string> TRExFit::FullHistogramPaths(Region *reg,
                                                     Sample *smp,
                                                     Systematic *syst,
                                                     bool isUp,
                                                     const bool isFolded,
                                                     const bool isSubtract) {
    // protection against nullptr
    if(reg==nullptr){
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(smp==nullptr){
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    bool isData = (smp->fType==Sample::SampleType::DATA);
    if(syst && !isData){
        if(isUp){
            if(syst->fHistoPathsUp.size()  >0) paths = syst->fHistoPathsUp;
            if(syst->fHistoFilesUp.size()  >0) files = syst->fHistoFilesUp;
            if(syst->fHistoNamesUp.size()  >0) names = syst->fHistoNamesUp;
        }
        else{
            if(syst->fHistoPathsDown.size()>0) paths = syst->fHistoPathsDown;
            if(syst->fHistoFilesDown.size()>0) files = syst->fHistoFilesDown;
            if(syst->fHistoNamesDown.size()>0) names = syst->fHistoNamesDown;
        }
        if (isSubtract) { // "special" systematics like JER
            if (isUp) {
                if (!syst->fHistoPathsUpSubtractSample.empty()) paths = syst->fHistoPathsUpSubtractSample;
                if (!syst->fHistoFilesUpSubtractSample.empty()) files = syst->fHistoFilesUpSubtractSample;
                if (!syst->fHistoNamesUpSubtractSample.empty()) names = syst->fHistoNamesUpSubtractSample;
            } else {
                if (!syst->fHistoPathsDownSubtractSample.empty()) paths = syst->fHistoPathsDownSubtractSample;
                if (!syst->fHistoFilesDownSubtractSample.empty()) files = syst->fHistoFilesDownSubtractSample;
                if (!syst->fHistoNamesDownSubtractSample.empty()) names = syst->fHistoNamesDownSubtractSample;
            }
        }
    }
    if(paths.size()==0 && smp->fHistoPaths.size()>0) paths = smp->fHistoPaths;
    if(files.size()==0 && smp->fHistoFiles.size()>0) files = smp->fHistoFiles;
    if(names.size()==0 && smp->fHistoNames.size()>0) names = smp->fHistoNames;
    //
    if(paths.size()==0 && reg->fHistoPaths.size()>0) paths = reg->fHistoPaths;
    if(files.size()==0 && reg->fHistoFiles.size()>0) files = reg->fHistoFiles;
    if(names.size()==0 && reg->fHistoNames.size()>0) names = reg->fHistoNames;
    //
    if(paths.size()==0 && fHistoPaths.size()>0) paths = fHistoPaths;
    if(files.size()==0 && fHistoFiles.size()>0) files = fHistoFiles;
    if(names.size()==0 && fHistoNames.size()>0) names = fHistoNames;

    if (!smp->fHistoFolderName.empty()) {
        for (auto& iname : names) {
            iname = Common::ReplaceFolderName(iname, smp->fHistoFolderName);
        }
    }

    // now we need to propagate the folder name updates
    if (syst && !isData) {

      auto broadcast = [&](auto& updates, std::size_t idx){
          if (updates.size() == names.size())
              return updates.at(idx);
          if (updates.size() == 1)
              return updates.at(0);
          return names.at(idx);
      };

        if (isSubtract) {
            if (isUp) {
                if (syst->fHistoFolderSubtractNamesUp.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fHistoFolderSubtractNamesUp, i));
                    }
                }
            } else {
                if (syst->fHistoFolderSubtractNamesDown.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fHistoFolderSubtractNamesDown, i));
                    }
                }
            }
        } else {
            if (isUp) {
                if (syst->fHistoFolderNamesUp.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fHistoFolderNamesUp, i));
                    }
                }
            } else {
                if (syst->fHistoFolderNamesDown.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fHistoFolderNamesDown, i));
                    }
                }
            }
        }

        //
        // now combining suffs, with all the combinations instead of giving priority
        // (same order as above used for suffix order - OK? FIXME)
        if(isUp) pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs, isFolded ), Common::ToVec(syst->fHistoPathSufUp), isFolded );
        else     pathSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs, isFolded ), Common::ToVec(syst->fHistoPathSufDown), isFolded );
    } else {
        pathSuffs = Common::CombinePathSufs( reg->fHistoPathSuffs, smp->fHistoPathSuffs, isFolded );
    }
    //
    if(syst && !isData){
        if (isSubtract) { // "special" systematics like JER
            if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs, isFolded ), Common::ToVec(syst->fHistoFileSufUpSubtractSample), isFolded );
            else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs, isFolded ), Common::ToVec(syst->fHistoFileSufDownSubtractSample), isFolded );
        } else {
            if(isUp) fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs, isFolded ), Common::ToVec(syst->fHistoFileSufUp), isFolded );
            else     fileSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs, isFolded ), Common::ToVec(syst->fHistoFileSufDown), isFolded );
        }
    }
    else{
        fileSuffs = Common::CombinePathSufs( reg->fHistoFileSuffs, smp->fHistoFileSuffs, isFolded );
    }
    //
    if(syst && !isData){
        if(isSubtract) {
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs, isFolded ), Common::ToVec(syst->fHistoNameSufUpSubtractSample), isFolded );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs, isFolded ), Common::ToVec(syst->fHistoNameSufDownSubtractSample), isFolded );
        }
        else {
            if(isUp) nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs, isFolded ), Common::ToVec(syst->fHistoNameSufUp), isFolded );
            else     nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs, isFolded ), Common::ToVec(syst->fHistoNameSufDown), isFolded );
        }
    }
    else{
      nameSuffs = Common::CombinePathSufs( Common::CombinePathSufs( reg->fHistoNameSuffs, smp->fHistoNameSuffs, isFolded), fHistoNamesNominal, isFolded );
    }
    //
    // And finally put everything together
    fullPaths = Common::CreatePathsList( paths,pathSuffs, files,fileSuffs, names,nameSuffs );
    return fullPaths;
}

//__________________________________________________________________________________
//
std::shared_ptr<SampleHist> TRExFit::GetSampleHistFromName(const Region* const reg, const std::string& name) const{
    for (const auto& isample : reg->fSampleHists) {
       if (isample->fName == name){
            return isample;
        }
    }
    return nullptr;
}

//__________________________________________________________________________________
//
TH1* TRExFit::CopySmoothedHisto(const SampleHist* const sh, const TH1* const nominal, const TH1* const up, const TH1* const down, const bool isUp) const{
    TH1* currentNominal = sh->fHist.get();

    if (currentNominal->GetNbinsX() != nominal->GetNbinsX()){
        LOG(ERROR) << "Histograms to smooth have different binning (nominal)\n";
        exit(EXIT_FAILURE);
    }
    if (currentNominal->GetNbinsX() != up->GetNbinsX()){
        LOG(ERROR) << "Histograms to smooth have different binning (nominal vs up)\n";
        exit(EXIT_FAILURE);
    }
    if (currentNominal->GetNbinsX() != down->GetNbinsX()){
        LOG(ERROR) << "Histograms to smooth have different binning (nominal vs down)\n";
        exit(EXIT_FAILURE);
    }

    TH1* result = static_cast<TH1*>(currentNominal->Clone());

    for (int ibin = 1; ibin <= currentNominal->GetNbinsX(); ++ibin){
        double ratio = 1;
        if (nominal->GetBinContent(ibin) != 0){
            if (isUp){
                ratio = up->GetBinContent(ibin)/nominal->GetBinContent(ibin);
            } else {
                ratio = down->GetBinContent(ibin)/nominal->GetBinContent(ibin);
            }
        }
        ratio*= currentNominal->GetBinContent(ibin);

        result->SetBinContent(ibin, ratio);
    }

    return result;
}

//__________________________________________________________________________________
//
int TRExFit::GetSystIndex(const SampleHist* const sh, const std::string& name) const{
    for (std::size_t i = 0; i < sh->fSyst.size(); ++i){
        if (sh->fSyst[i]->fName == name){
            return i;
        }
    }

    return -1;
}

//__________________________________________________________________________________
//
std::shared_ptr<SystematicHist> TRExFit::CombineSpecialHistos(std::shared_ptr<SystematicHist> orig,
                                                              const std::vector<std::shared_ptr<SystematicHist> >& vec,
                                                              Systematic::COMBINATIONTYPE type,
                                                              const SampleHist* sh) const {
    if (vec.size() == 0) return nullptr;
    if (!orig) return nullptr;
    if (!orig->fHistUp) return nullptr;
    if (!orig->fHistDown) return nullptr;
    if (!sh) return nullptr;
    if (!sh->fHist) return nullptr;

    const int nbins = sh->fHist->GetNbinsX();
    if (type == Systematic::COMBINATIONTYPE::ENVELOPE) {
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            std::vector<double> hist_max;
            std::vector<double> hist_min;
            for (const auto& isyst : vec) {
                if (!isyst) continue;
                if (!isyst->fHistUp) continue;
                if (!isyst->fHistDown) continue;
                const double diff_up = isyst->fHistUp->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                const double diff_down = isyst->fHistDown->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                if (diff_up > 0) {
                    hist_max.emplace_back(diff_up);
                } else {
                    hist_min.emplace_back(diff_up);
                }
                if (diff_down > 0) {
                    hist_max.emplace_back(diff_down);
                } else {
                    hist_min.emplace_back(diff_down);
                }
            }

            double max(0);
            double min(0);

            if (hist_max.size() > 0) {
                max = *std::max_element(hist_max.begin(), hist_max.end());
            }
            if (hist_min.size() > 0) {
                min = *std::min_element(hist_min.begin(), hist_min.end());
            }

            if (max >= 0) {
                orig->fHistUp->SetBinContent(ibin, max + sh->fHist->GetBinContent(ibin));
            }
            if (min <= 0) {
                orig->fHistDown->SetBinContent(ibin, min + sh->fHist->GetBinContent(ibin));
            }
        }
    } else if (type == Systematic::COMBINATIONTYPE::STANDARDDEVIATION
                || type == Systematic::COMBINATIONTYPE::STANDARDDEVIATIONNODDOF
                || type == Systematic::COMBINATIONTYPE::HESSIAN) {
        for (int ibin = 1; ibin <= nbins; ++ ibin) {
            std::vector<double> content;
            for (const auto& isyst : vec) {
                if (!isyst) continue;
                if (!isyst->fHistUp) continue;
                if (!isyst->fHistDown) continue;
                const double up = isyst->fHistUp->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                const double down = isyst->fHistDown->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                // take the larger (signed) variation
                const double max = std::max(std::abs(up),std::abs(down));
                content.emplace_back(max);
            }

            if (content.size() < 2) continue;

            // calculate the standard deviation in each bin
            double sum(0);
            double sigma(0);

            for (const auto& value : content) {
                sum+= value;
            }

            // get mean
            const double mean = static_cast<double>(sum / content.size());

            // get std
            for (const auto& value : content) {
                sigma+= (value-mean)*(value-mean);
            }

            double denominator(0.);

            if (type == Systematic::COMBINATIONTYPE::STANDARDDEVIATIONNODDOF) {
                denominator = content.size();
            } else if (type == Systematic::COMBINATIONTYPE::HESSIAN) {
                denominator = 1.;
            } else { // normal standard deviation
                denominator = content.size() - 1;
            }

            sigma = static_cast<double>(sigma / denominator);
            sigma = std::sqrt(sigma);

            // now set the bin content
            orig->fHistUp->SetBinContent(ibin,    sigma + sh->fHist->GetBinContent(ibin));
            orig->fHistDown->SetBinContent(ibin, -sigma + sh->fHist->GetBinContent(ibin));
        }
    } else if (type == Systematic::COMBINATIONTYPE::SUMINSQUARES) {
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            double up(0);
            double down(0);
            for (const auto& isyst : vec) {
                if (!isyst) continue;
                if (!isyst->fHistUp) continue;
                if (!isyst->fHistDown) continue;
                const double diff_up = isyst->fHistUp->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                const double diff_down = isyst->fHistDown->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);

                if (diff_up > 0) {
                   up = std::hypot(diff_up, up);
                } else {
                   down = -std::hypot(diff_up, down);
                }
                if (diff_down > 0) {
                   up = std::hypot(diff_down, up);
                } else {
                   down = -std::hypot(diff_down, down);
                }
            }

            orig->fHistUp->  SetBinContent(ibin, sh->fHist->GetBinContent(ibin) + up);
            orig->fHistDown->SetBinContent(ibin, sh->fHist->GetBinContent(ibin) + down);
        }
    } else {
        LOG(ERROR) << "Unknown combination type, this should not happen!\n";
        exit(EXIT_FAILURE);
    }

    // modify the shape histos as well
    orig->fHistShapeUp.reset(static_cast<TH1*>(orig->fHistUp->Clone()));
    orig->fHistShapeUp->Scale(Common::EffIntegral(sh->fHist.get())/Common::EffIntegral(orig->fHistUp.get()));
    orig->fHistShapeDown.reset(static_cast<TH1*>(orig->fHistDown->Clone()));
    orig->fHistShapeDown->Scale(Common::EffIntegral(sh->fHist.get())/Common::EffIntegral(orig->fHistDown.get()));
    orig->fHistShapeUp->SetDirectory(nullptr);
    orig->fHistShapeDown->SetDirectory(nullptr);
    orig->fHistShapeUp->SetName(orig->fHistoNameShapeUp.c_str());
    orig->fHistShapeDown->SetName(orig->fHistoNameShapeDown.c_str());
    orig->fHasShape = true;

    return orig;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::GetUniqueSystNamesWithoutGamma() const {
    std::vector<std::string> result;
    for(const auto& isyst : fSystematics) {
        if (isyst->fType == Systematic::SHAPE) continue;
        if (std::find(result.begin(), result.end(), isyst->fName) == result.end()) {
            result.emplace_back(isyst->fName);
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<Region*> TRExFit::GetNonValidationRegions() const {
    std::vector<Region*> result;
    for (const auto& ireg : fRegions) {
        if (ireg->fRegionType == Region::VALIDATION) continue;
        result.emplace_back(ireg.get());
    }

    return result;
}

//__________________________________________________________________________________
//
std::vector<std::shared_ptr<Sample> > TRExFit::GetNonDataNonGhostSamples() const {
    std::vector<std::shared_ptr<Sample> > result;

    for (const auto& isample : fSamples) {
        if (isample->fType == Sample::SampleType::DATA) continue;
        if (isample->fType == Sample::SampleType::GHOST) continue;
        if (isample->fType == Sample::SampleType::EFT) continue;

        result.emplace_back(isample);
    }

    return result;
}

//__________________________________________________________________________________
//
void TRExFit::PrepareUnfolding() {
    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        PrepareUnfolding(fUnfolding.at(i).get(), i);
    }
}
//__________________________________________________________________________________
//
void TRExFit::PrepareUnfolding(const Unfolding* unfolding, const int index) {
    // Prepare the folder structure
    gSystem->mkdir(fName.c_str());
    gSystem->mkdir((fName+"/UnfoldingHistograms").c_str());

    // Open the otuput ROOT file
    std::unique_ptr<TFile> outputFile(nullptr);
    if (index == 0) {
        outputFile.reset(TFile::Open((fName+"/UnfoldingHistograms/FoldedHistograms.root").c_str(), "RECREATE"));
    } else {
        outputFile.reset(TFile::Open((fName+"/UnfoldingHistograms/FoldedHistograms.root").c_str(), "UPDATE"));
    }
    if (!outputFile) {
        LOG(ERROR) << "Cannot open the output file at: " << fName << "/UnfoldingHistograms/FoldedHistograms.root\n";
        exit(EXIT_FAILURE);
    }

    FoldingManager manager{};
    manager.SetMatrixOrientation(unfolding->fMatrixOrientation);

    const bool horizontal = (unfolding->fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);

    std::unique_ptr<TH1> alternativeTruth(nullptr);

    {
        std::unique_ptr<TH1> truth(nullptr);
        for (const auto& itruth : fTruthSamples) {
            if (itruth->GetName() == unfolding->fNominalTruthSample) {
                // check if it is assigned to the correct Unfolding object
                if (itruth->fUnfoldingName != unfolding->fName) {
                    LOG(ERROR) << "The assigned truth sample: " << itruth->GetName() << " does not have UnfoldingName set properly!\n";
                    exit(EXIT_FAILURE);
                }
                truth = itruth->GetHisto(unfolding);
                break;
            }
        }
        if (!truth) {
            LOG(ERROR) << "The truth histogram is nullptr\n";
            exit(EXIT_FAILURE);
        }
        if (truth->GetNbinsX() != unfolding->fNumberUnfoldingTruthBins) {
            LOG(ERROR) << "The number of truth bins doesnt match the value from the config\n";
            exit(EXIT_FAILURE);
        }
        manager.SetTruthDistribution(truth.get());
        manager.WriteTruthToHisto(outputFile.get(), "", unfolding->fName + "_truth_distribution");

        // read the alternative truth sample
        if (unfolding->fAlternativeAsimovTruthSample != "") {
            for (const auto& itruth : fTruthSamples) {
                if (itruth->GetName() == unfolding->fAlternativeAsimovTruthSample) {
                    if (itruth->fUnfoldingName != unfolding->fName) {
                        LOG(ERROR) << "The assigned truth sample: " << itruth->GetName() << " does not have UnfoldingName set properly!\n";
                        exit(EXIT_FAILURE);
                    }
                    alternativeTruth = itruth->GetHisto(unfolding);
                    break;
                }
            }
        }
    }

    // loop over regions
    for (const auto& ireg : fRegions) {

        // only signal regions are processed at this step
        if (ireg->fRegionType != Region::RegionType::SIGNAL) continue;

        // loop over all samples
        for (const auto& isample : fUnfoldingSamples) {
            if (isample->fUnfoldingName != unfolding->fName) continue;

            // skip samples not associated to the region
            if(isample->fRegions[0] != "all" &&
                Common::FindInStringVector(isample->fRegions, ireg->fName) < 0) continue;

            // first process nominal
            if (isample->GetHasResponse()) {
                const std::vector<std::string>& fullResponsePaths = FullResponseMatrixPaths(ireg.get(), isample.get());

                std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(fullResponsePaths);
                if (!matrix) {
                    LOG(ERROR) << "Cannot read the response matrix!\n";
                    exit(EXIT_FAILURE);
                }
                const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
                const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
                if (ireg->fNumberUnfoldingRecoBins > 0 && nRecoBins != ireg->fNumberUnfoldingRecoBins) {
                    LOG(ERROR) << "Number of reco bins do not match the number of reco bins for the response matrix in region: " << ireg->fName << "\n";
                    exit(EXIT_FAILURE);
                }
                if (nTruthBins != unfolding->fNumberUnfoldingTruthBins) {
                    LOG(ERROR) << "Number of truth bins do not match the number of truth bins for the response matrix in regoin: " << ireg->fName << "\n";
                    exit(EXIT_FAILURE);
                }

                if (isample->GetType() != UnfoldingSample::TYPE::GHOST) PlotMigrationResponse(matrix.get(), false, ireg->fName, "", unfolding);

                manager.SetResponseMatrix(matrix.get());
                outputFile->cd();
                matrix->Write((unfolding->fName + "_" + ireg->fName + "_" + isample->GetName() + "_response").c_str());
            } else {
                // need to add acceptance, selection and migration
                {
                    const std::vector<std::string>& fullMigrationMatrixPaths = FullMigrationMatrixPaths(ireg.get(), isample.get());
                    std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(fullMigrationMatrixPaths);
                    if (!matrix) {
                        exit(EXIT_FAILURE);
                    }
                    const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
                    const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
                    if (ireg->fNumberUnfoldingRecoBins > 0 && nRecoBins != ireg->fNumberUnfoldingRecoBins) {
                        LOG(ERROR) << "Number of reco bins do not match the number of reco bins for the migration matrix in region: " << ireg->fName << "\n";
                        exit(EXIT_FAILURE);
                    }
                    if (nTruthBins != unfolding->fNumberUnfoldingTruthBins) {
                        LOG(ERROR) << "Number of truth bins do not match the number of truth bins for the migration matrix in region: " << ireg->fName << "\n";
                        exit(EXIT_FAILURE);
                    }

                    UnfoldingTools::NormalizeMatrix(matrix.get(), !horizontal);

                    if (isample->GetType() != UnfoldingSample::TYPE::GHOST) PlotMigrationResponse(matrix.get(), true, ireg->fName, "", unfolding);

                    // pass the migration to the tool
                    manager.SetMigrationMatrix(matrix.get(), false);
                    outputFile->cd();
                    matrix->Write((unfolding->fName + "_" + ireg->fName + "_" + isample->GetName() + "_migration").c_str());
                }

                // add selection eff
                {
                    const std::vector<std::string>& fullSelectionEffPaths = FullSelectionEffPaths(ireg.get(), isample.get());
                    std::unique_ptr<TH1> eff = Common::CombineHistosFromFullPaths(fullSelectionEffPaths);
                    if (!eff) {
                        exit(EXIT_FAILURE);
                    }
                    const int nbins = eff->GetNbinsX();
                    if (nbins != unfolding->fNumberUnfoldingTruthBins) {
                        LOG(ERROR) << "Number of efficiency selection bins doesnt match the number of truth bins\n";
                        exit(EXIT_FAILURE);
                    }

                    manager.SetSelectionEfficiency(eff.get());
                }

                // add acceptance
                if (fHasAcceptance || isample->GetHasAcceptance() || ireg->fHasAcceptance) {
                    const std::vector<std::string>& fullAcceptancePaths = FullAcceptancePaths(ireg.get(), isample.get());
                    std::unique_ptr<TH1> acc = Common::CombineHistosFromFullPaths(fullAcceptancePaths);
                    if (!acc) {
                        exit(EXIT_FAILURE);
                    }
                    const int nbins = acc->GetNbinsX();
                    if (ireg->fNumberUnfoldingRecoBins > 0 && nbins != ireg->fNumberUnfoldingRecoBins) {
                        LOG(ERROR) << "Number of acceptance bins doesnt match the number of reco bins in region " << ireg->fName << "\n";
                        exit(EXIT_FAILURE);
                    }

                    manager.SetAcceptance(acc.get());
                }

                manager.CalculateResponseMatrix(true);
                if (isample->GetType() != UnfoldingSample::TYPE::GHOST) PlotMigrationResponse(manager.GetResponseMatrix(), false, ireg->fName, "", unfolding);
                outputFile->cd();
                manager.GetResponseMatrix()->Write((unfolding->fName + "_" + ireg->fName + "_" + isample->GetName() + "_response").c_str());
            }

            std::unique_ptr<TH2> nominal(static_cast<TH2*>(manager.GetResponseMatrix()->Clone()));

            manager.FoldTruth();

            // create the folder structure
            const TDirectory* dir = dynamic_cast<const TDirectory*>(outputFile->Get("nominal"));
            if (!dir) {
                outputFile->cd();
                outputFile->mkdir("nominal");
            }

            const std::string histoName = unfolding->fName + "_" + ireg->fName + "_" + isample->GetName();
            manager.WriteFoldedToHisto(outputFile.get(), "nominal", histoName);

            // fold and store the altenative truth sample
            if (unfolding->fAlternativeAsimovTruthSample != "") {
                std::unique_ptr<TH1> tmp = manager.TotalFold(alternativeTruth.get());
                outputFile->cd();
                tmp->Write((unfolding->fName + "_" + ireg->fName + "_AlternativeAsimov").c_str());
            }

            if (isample->GetType() == UnfoldingSample::TYPE::GHOST) continue;

            // Process systematics
            for (const auto& isyst : fUnfoldingSystematics) {
                if (!isyst) continue;
                if (isyst->GetName() == "Dummy") continue;
                if (isyst->fUnfoldingName != unfolding->fName) continue;

                if(isyst->fRegions.at(0) != "all" &&
                     Common::FindInStringVector(isyst->fRegions, ireg->fName) < 0) continue;
                if(isyst->fSamples.at(0) != "all" &&
                     Common::FindInStringVector(isyst->fSamples, isample->GetName()) < 0) continue;

                ProcessUnfoldingSystematics(unfolding,
                                            &manager,
                                            outputFile.get(),
                                            ireg.get(),
                                            isample.get(),
                                            isyst.get(),
                                            nominal.get());
            }
        }
    }

    outputFile->Close();
}

//__________________________________________________________________________________
//
void TRExFit::ProcessUnfoldingSystematics(const Unfolding* unfolding,
                                          FoldingManager* manager,
                                          TFile* file,
                                          const Region* reg,
                                          const UnfoldingSample* sample,
                                          const UnfoldingSystematic* syst,
                                          const TH2* nominal) const {

    const bool horizontal = (unfolding->fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);

    // lambda to propagate reference sample
    auto PropagateRefSample = [&](TH2* response) {
        auto it = std::find_if(fUnfoldingSamples.begin(), fUnfoldingSamples.end(), [&syst](const std::unique_ptr<UnfoldingSample>& smp) {
            return syst->GetReferenceSample() == smp->GetName();
        });

        if (it == fUnfoldingSamples.end()) {
            LOG(ERROR) << "Cannot find the reference sample\n";
            exit(EXIT_FAILURE);
        }


        std::unique_ptr<TH2> referenceMatrix(nullptr);;
        if (syst->GetHasResponse()) {
            const std::vector<std::string>& paths = FullResponseMatrixPaths(reg, it->get());
            referenceMatrix = Common::CombineHistos2DFromFullPaths(paths);
        } else {
            FoldingManager mgr{};
            mgr.SetMatrixOrientation(unfolding->fMatrixOrientation);
            {
                const std::vector<std::string>& paths = FullMigrationMatrixPaths(reg, it->get());
                std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(paths);
                if (!matrix) {
                    exit(EXIT_FAILURE);
                }
                UnfoldingTools::NormalizeMatrix(matrix.get(), !horizontal);
                mgr.SetMigrationMatrix(matrix.get(), false);
            }
            {
                const std::vector<std::string>& paths = FullSelectionEffPaths(reg, it->get());
                std::unique_ptr<TH1> eff = Common::CombineHistosFromFullPaths(paths);
                if (!eff) {
                    exit(EXIT_FAILURE);
                }

                mgr.SetSelectionEfficiency(eff.get());
            }

            if (fHasAcceptance || syst->GetHasAcceptance() || reg->fHasAcceptance) {
                const std::vector<std::string>& paths = FullAcceptancePaths(reg, it->get());
                std::unique_ptr<TH1> acc = Common::CombineHistosFromFullPaths(paths);
                if (!acc) {
                    exit(EXIT_FAILURE);
                }

                mgr.SetAcceptance(acc.get());
            }
            mgr.CalculateResponseMatrix(true);
            referenceMatrix.reset(static_cast<TH2*>(mgr.GetResponseMatrix()->Clone()));
        }

        if (!referenceMatrix) {
            LOG(ERROR) << "ReferenceResponse is nullptr\n";
            exit(EXIT_FAILURE);
        }
        response->Add(referenceMatrix.get(), -1.);
        response->Add(nominal);
    };

    // lambda for processing one (up or down) variation
    auto ProcessOneVariation = [&](const bool isUp, const bool isSubtract) {
        if (isSubtract && !syst->HasSubtract()) return;
        if (syst->GetHasResponse()) {
            const std::vector<std::string>& paths = FullResponseMatrixPaths(reg, sample, syst, isUp, isSubtract);
            std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(paths);
            if (!matrix) {
                exit(EXIT_FAILURE);
            }
            const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
            const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
            if (reg->fNumberUnfoldingRecoBins > 0 && nRecoBins != reg->fNumberUnfoldingRecoBins) {
                LOG(ERROR) << "Number of reco bins do not match the number of reco bins for the response matrix in region: " << reg->fName << "\n";
                exit(EXIT_FAILURE);
            }
            if (nTruthBins != unfolding->fNumberUnfoldingTruthBins) {
                LOG(ERROR) << "Number of truth bins do not match the number of truth bins for the response matrix: " << reg->fName << "\n";
                exit(EXIT_FAILURE);
            }
            if (unfolding->fPlotSystematicMigrations) {
                PlotMigrationResponse(matrix.get(), false, reg->fName, syst->GetName(), unfolding);
            }

            if (syst->GetReferenceSample() != "") {
                PropagateRefSample(matrix.get());
            }

            manager->SetResponseMatrix(matrix.get());

            if (unfolding->fPlotSystematicMigrations) {
                file->cd();
                matrix->Write((unfolding->fName + "_" + reg->fName + "_" + syst->GetName() + "_response").c_str());
            }
        } else {
            {
                /// migration first
                const std::vector<std::string>& paths = FullMigrationMatrixPaths(reg, sample, syst, isUp, isSubtract);
                std::unique_ptr<TH2> matrix = Common::CombineHistos2DFromFullPaths(paths);
                if (!matrix) {
                    exit(EXIT_FAILURE);
                }

                const int nRecoBins  = horizontal ? matrix->GetNbinsY() : matrix->GetNbinsX();
                const int nTruthBins = horizontal ? matrix->GetNbinsX() : matrix->GetNbinsY();
                if (reg->fNumberUnfoldingRecoBins > 0 && nRecoBins != reg->fNumberUnfoldingRecoBins) {
                    LOG(ERROR) << "Number of reco bins do not match the number of reco bins for the migration matrix in region: " << reg->fName << "\n";
                    exit(EXIT_FAILURE);
                }
                if (nTruthBins != unfolding->fNumberUnfoldingTruthBins) {
                    LOG(ERROR) << "Number of truth bins do not match the number of truth bins for the migration matrix: " << reg->fName << "\n";
                    exit(EXIT_FAILURE);
                }
                UnfoldingTools::NormalizeMatrix(matrix.get(), !horizontal);

                if (unfolding->fPlotSystematicMigrations) {
                    PlotMigrationResponse(matrix.get(), true, reg->fName, syst->GetName(), unfolding);
                }

                manager->SetMigrationMatrix(matrix.get(), false);
                if (unfolding->fPlotSystematicMigrations) {
                    file->cd();
                    matrix->Write((unfolding->fName + "_" + reg->fName + "_" + syst->GetName() + "_migration").c_str());
                }
            }

            // selectio eff now
            {
                const std::vector<std::string>& paths = FullSelectionEffPaths(reg, sample, syst, isUp, isSubtract);
                std::unique_ptr<TH1> eff = Common::CombineHistosFromFullPaths(paths);
                if (!eff) {
                    exit(EXIT_FAILURE);
                }

                const int nbins = eff->GetNbinsX();
                if (nbins != unfolding->fNumberUnfoldingTruthBins) {
                    LOG(ERROR) << "Number of efficiency selection bins doesnt match the number of truth bins\n";
                    exit(EXIT_FAILURE);
                }

                manager->SetSelectionEfficiency(eff.get());
            }

            if (fHasAcceptance || syst->GetHasAcceptance() || reg->fHasAcceptance) {
                const std::vector<std::string>& paths = FullAcceptancePaths(reg, sample, syst, isUp, isSubtract);
                std::unique_ptr<TH1> acc = Common::CombineHistosFromFullPaths(paths);
                if (!acc) {
                    exit(EXIT_FAILURE);
                }

                const int nbins = acc->GetNbinsX();
                if (reg-> fNumberUnfoldingRecoBins > 0 && nbins != reg->fNumberUnfoldingRecoBins) {
                    LOG(ERROR) << "Number of acceptance bins doesnt match the number of reco bins in region " << reg->fName << "\n";
                    exit(EXIT_FAILURE);
                }

                manager->SetAcceptance(acc.get());
            }
            manager->CalculateResponseMatrix(true);

            std::unique_ptr<TH2> matrix(static_cast<TH2*>(manager->GetResponseMatrix()->Clone()));
            if (syst->GetReferenceSample() != "") {
                PropagateRefSample(matrix.get());
                manager->SetResponseMatrix(matrix.get());
            }

            if (unfolding->fPlotSystematicMigrations) {
                PlotMigrationResponse(matrix.get(), false, reg->fName, syst->GetName(), unfolding);
                file->cd();
                const std::string finalName = unfolding->fName + "_" + reg->fName + "_" + (isSubtract ? syst->GetName() + "_subtract" : syst->GetName()) + "_response";
                matrix->Write(finalName.c_str());
            }
        }
        manager->FoldTruth();

        // create the folder structure
        const std::string folderName = isUp ? (isSubtract ? syst->GetName() + "_subtract_Up" : syst->GetName() + "_Up")
                                            : (isSubtract ? syst->GetName() + "_subtract_Down" : syst->GetName() + "_Down");
        const TDirectory* dir = dynamic_cast<const TDirectory*>(file->Get(folderName.c_str()));
        if (!dir) {
            file->cd();
            file->mkdir(folderName.c_str());
        }

        const std::string histoName = unfolding->fName + "_" + reg->fName + "_" + sample->GetName();
        manager->WriteFoldedToHisto(file, folderName, histoName);
    };

    if (syst->fHasUpVariation) {
        ProcessOneVariation(true, false);
        ProcessOneVariation(true, true);
    }
    if (syst->fHasDownVariation) {
        ProcessOneVariation(false, false);
        ProcessOneVariation(false, true);
    }
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullResponseMatrixPaths(const Region* reg,
                                                          const UnfoldingSample* smp,
                                                          const UnfoldingSystematic* syst,
                                                          const bool isUp,
                                                          const bool isSubtract) const {
    // protection against nullptr
    if(!reg) {
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fResponseMatrixPathsUp.size()  >0) paths = syst->fResponseMatrixPathsUp;
            if(syst->fResponseMatrixFilesUp.size()  >0) files = syst->fResponseMatrixFilesUp;
            if(syst->fResponseMatrixNamesUp.size()  >0) names = syst->fResponseMatrixNamesUp;
        } else {
            if(syst->fResponseMatrixPathsDown.size()>0) paths = syst->fResponseMatrixPathsDown;
            if(syst->fResponseMatrixFilesDown.size()>0) files = syst->fResponseMatrixFilesDown;
            if(syst->fResponseMatrixNamesDown.size()>0) names = syst->fResponseMatrixNamesDown;
        }
        if (isSubtract) {
            if (isUp) {
                if (!syst->fResponseMatrixPathsUpSubtractSample.empty()) paths = syst->fResponseMatrixPathsUpSubtractSample;
                if (!syst->fResponseMatrixFilesUpSubtractSample.empty()) files = syst->fResponseMatrixFilesUpSubtractSample;
                if (!syst->fResponseMatrixNamesUpSubtractSample.empty()) names = syst->fResponseMatrixNamesUpSubtractSample;
            } else {
                if (!syst->fResponseMatrixPathsDownSubtractSample.empty()) paths = syst->fResponseMatrixPathsDownSubtractSample;
                if (!syst->fResponseMatrixFilesDownSubtractSample.empty()) files = syst->fResponseMatrixFilesDownSubtractSample;
                if (!syst->fResponseMatrixNamesDownSubtractSample.empty()) names = syst->fResponseMatrixNamesDownSubtractSample;
            }
        }
    }
    if(paths.size()==0 && smp->fResponseMatrixPaths.size()>0) paths = smp->fResponseMatrixPaths;
    if(files.size()==0 && smp->fResponseMatrixFiles.size()>0) files = smp->fResponseMatrixFiles;
    if(names.size()==0 && smp->fResponseMatrixNames.size()>0) names = smp->fResponseMatrixNames;

    if(paths.size()==0 && reg->fResponseMatrixPaths.size()>0) paths = reg->fResponseMatrixPaths;
    if(files.size()==0 && reg->fResponseMatrixFiles.size()>0) files = reg->fResponseMatrixFiles;
    if(names.size()==0 && reg->fResponseMatrixNames.size()>0) names = reg->fResponseMatrixNames;

    if(paths.size()==0 && fResponseMatrixPaths.size()>0) paths = fResponseMatrixPaths;
    if(files.size()==0 && fResponseMatrixFiles.size()>0) files = fResponseMatrixFiles;
    if(names.size()==0 && fResponseMatrixNames.size()>0) names = fResponseMatrixNames;

    // now we need to propagate the folder name updates
    if (syst) {

        auto broadcast = [&](auto& updates, std::size_t idx){
            if (updates.size() == names.size())
                return updates.at(idx);
            if (updates.size() == 1)
                return updates.at(0);
            return names.at(idx);
        };

        if (isSubtract) {
            if (isUp) {
                if (syst->fResponseMatrixFolderNamesUpSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fResponseMatrixFolderNamesUpSubtractSample, i));
                    }
                }
            } else {
                if (syst->fResponseMatrixFolderNamesDownSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fResponseMatrixFolderNamesDownSubtractSample, i));
                    }
                }
            }
        } else {
            if (isUp) {
                if (syst->fResponseMatrixFolderNamesUp.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fResponseMatrixFolderNamesUp, i));
                    }
                }
            } else {
                if (syst->fResponseMatrixFolderNamesDown.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fResponseMatrixFolderNamesDown, i));
                    }
                }
            }
        }
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixPathSuffs,smp->fResponseMatrixPathSuffs), syst->fResponseMatrixPathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixPathSuffs,smp->fResponseMatrixPathSuffs), syst->fResponseMatrixPathSuffsDown);
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixFileSuffs, smp->fResponseMatrixFileSuffs), syst->fResponseMatrixFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixFileSuffs, smp->fResponseMatrixFileSuffs), syst->fResponseMatrixFileSuffsDown);
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixNameSuffs, smp->fResponseMatrixNameSuffs), syst->fResponseMatrixNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixNameSuffs, smp->fResponseMatrixNameSuffs), syst->fResponseMatrixNameSuffsDown);
    } else {
        pathSuffs = Common::CombinePathSufs(reg->fResponseMatrixPathSuffs, smp->fResponseMatrixPathSuffs);
        fileSuffs = Common::CombinePathSufs(reg->fResponseMatrixFileSuffs, smp->fResponseMatrixFileSuffs);
        nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fResponseMatrixNameSuffs, smp->fResponseMatrixNameSuffs), fResponseMatrixNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullMigrationMatrixPaths(const Region* reg,
                                                           const UnfoldingSample* smp,
                                                           const UnfoldingSystematic* syst,
                                                           const bool isUp,
                                                           const bool isSubtract) const {
    // protection against nullptr
    if(!reg) {
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fMigrationPathsUp.size()  >0) paths = syst->fMigrationPathsUp;
            if(syst->fMigrationFilesUp.size()  >0) files = syst->fMigrationFilesUp;
            if(syst->fMigrationNamesUp.size()  >0) names = syst->fMigrationNamesUp;
        } else {
            if(syst->fMigrationPathsDown.size()>0) paths = syst->fMigrationPathsDown;
            if(syst->fMigrationFilesDown.size()>0) files = syst->fMigrationFilesDown;
            if(syst->fMigrationNamesDown.size()>0) names = syst->fMigrationNamesDown;
        }
        if (isSubtract) {
            if (isUp) {
                if (!syst->fMigrationPathsUpSubtractSample.empty()) paths = syst->fMigrationPathsUpSubtractSample;
                if (!syst->fMigrationFilesUpSubtractSample.empty()) files = syst->fMigrationFilesUpSubtractSample;
                if (!syst->fMigrationNamesUpSubtractSample.empty()) names = syst->fMigrationNamesUpSubtractSample;
            } else {
                if (!syst->fMigrationPathsDownSubtractSample.empty()) paths = syst->fMigrationPathsDownSubtractSample;
                if (!syst->fMigrationFilesDownSubtractSample.empty()) files = syst->fMigrationFilesDownSubtractSample;
                if (!syst->fMigrationNamesDownSubtractSample.empty()) names = syst->fMigrationNamesDownSubtractSample;
            }
        }
    }
    if(paths.size()==0 && smp->fMigrationPaths.size()>0) paths = smp->fMigrationPaths;
    if(files.size()==0 && smp->fMigrationFiles.size()>0) files = smp->fMigrationFiles;
    if(names.size()==0 && smp->fMigrationNames.size()>0) names = smp->fMigrationNames;

    if(paths.size()==0 && reg->fMigrationPaths.size()>0) paths = reg->fMigrationPaths;
    if(files.size()==0 && reg->fMigrationFiles.size()>0) files = reg->fMigrationFiles;
    if(names.size()==0 && reg->fMigrationNames.size()>0) names = reg->fMigrationNames;

    if(paths.size()==0 && fMigrationPaths.size()>0) paths = fMigrationPaths;
    if(files.size()==0 && fMigrationFiles.size()>0) files = fMigrationFiles;
    if(names.size()==0 && fMigrationNames.size()>0) names = fMigrationNames;

    // now we need to propagate the folder name updates
    if (syst) {

        auto broadcast = [&](auto& updates, std::size_t idx){
            if (updates.size() == names.size())
                return updates.at(idx);
            if (updates.size() == 1)
                return updates.at(0);
            return names.at(idx);
        };

        if (isSubtract) {
            if (isUp) {
                if (syst->fMigrationFolderNamesUpSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fMigrationFolderNamesUpSubtractSample, i));
                    }
                }
            } else {
                if (syst->fMigrationFolderNamesDownSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fMigrationFolderNamesDownSubtractSample, i));
                    }
                }
            }
        } else {
            if (isUp) {
                if (syst->fMigrationFolderNamesUp.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fMigrationFolderNamesUp, i));
                    }
                }
            } else {
                if (syst->fMigrationFolderNamesDown.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fMigrationFolderNamesDown, i));
                    }
                }
            }
        }
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationPathSuffs,smp->fMigrationPathSuffs), syst->fMigrationPathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationPathSuffs,smp->fMigrationPathSuffs), syst->fMigrationPathSuffsDown);
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationFileSuffs, smp->fMigrationFileSuffs), syst->fMigrationFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationFileSuffs, smp->fMigrationFileSuffs), syst->fMigrationFileSuffsDown);
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationNameSuffs, smp->fMigrationNameSuffs), syst->fMigrationNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationNameSuffs, smp->fMigrationNameSuffs), syst->fMigrationNameSuffsDown);
    } else {
        pathSuffs = Common::CombinePathSufs(reg->fMigrationPathSuffs, smp->fMigrationPathSuffs);
        fileSuffs = Common::CombinePathSufs(reg->fMigrationFileSuffs, smp->fMigrationFileSuffs);
        nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fMigrationNameSuffs, smp->fMigrationNameSuffs), fMigrationNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullAcceptancePaths(const Region* reg,
                                                      const UnfoldingSample* smp,
                                                      const UnfoldingSystematic* syst,
                                                      const bool isUp,
                                                      const bool isSubtract) const {
    // protection against nullptr
    if(!reg) {
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fAcceptancePathsUp.size()  >0) paths = syst->fAcceptancePathsUp;
            if(syst->fAcceptanceFilesUp.size()  >0) files = syst->fAcceptanceFilesUp;
            if(syst->fAcceptanceNamesUp.size()  >0) names = syst->fAcceptanceNamesUp;
        } else {
            if(syst->fAcceptancePathsDown.size()>0) paths = syst->fAcceptancePathsDown;
            if(syst->fAcceptanceFilesDown.size()>0) files = syst->fAcceptanceFilesDown;
            if(syst->fAcceptanceNamesDown.size()>0) names = syst->fAcceptanceNamesDown;
        }
        if (isSubtract) {
            if (isUp) {
                if (!syst->fAcceptancePathsUpSubtractSample.empty()) paths = syst->fAcceptancePathsUpSubtractSample;
                if (!syst->fAcceptanceFilesUpSubtractSample.empty()) files = syst->fAcceptanceFilesUpSubtractSample;
                if (!syst->fAcceptanceNamesUpSubtractSample.empty()) names = syst->fAcceptanceNamesUpSubtractSample;
            } else {
                if (!syst->fAcceptancePathsDownSubtractSample.empty()) paths = syst->fAcceptancePathsDownSubtractSample;
                if (!syst->fAcceptanceFilesDownSubtractSample.empty()) files = syst->fAcceptanceFilesDownSubtractSample;
                if (!syst->fAcceptanceNamesDownSubtractSample.empty()) names = syst->fAcceptanceNamesDownSubtractSample;
            }
        }
    }
    if(paths.size()==0 && smp->fAcceptancePaths.size()>0) paths = smp->fAcceptancePaths;
    if(files.size()==0 && smp->fAcceptanceFiles.size()>0) files = smp->fAcceptanceFiles;
    if(names.size()==0 && smp->fAcceptanceNames.size()>0) names = smp->fAcceptanceNames;

    if(paths.size()==0 && reg->fAcceptancePaths.size()>0) paths = reg->fAcceptancePaths;
    if(files.size()==0 && reg->fAcceptanceFiles.size()>0) files = reg->fAcceptanceFiles;
    if(names.size()==0 && reg->fAcceptanceNames.size()>0) names = reg->fAcceptanceNames;

    if(paths.size()==0 && fAcceptancePaths.size()>0) paths = fAcceptancePaths;
    if(files.size()==0 && fAcceptanceFiles.size()>0) files = fAcceptanceFiles;
    if(names.size()==0 && fAcceptanceNames.size()>0) names = fAcceptanceNames;

    // now we need to propagate the folder name updates
    if (syst) {

        auto broadcast = [&](auto& updates, std::size_t idx){
            if (updates.size() == names.size())
                return updates.at(idx);
            if (updates.size() == 1)
                return updates.at(0);
            return names.at(idx);
        };

        if (isSubtract) {
            if (isUp) {
                if (syst->fAcceptanceFolderNamesUpSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fAcceptanceFolderNamesUpSubtractSample, i));
                    }
                }
            } else {
                if (syst->fAcceptanceFolderNamesDownSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fAcceptanceFolderNamesDownSubtractSample, i));
                    }
                }
            }
        } else {
            if (isUp) {
                if (syst->fAcceptanceFolderNamesUp.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fAcceptanceFolderNamesUp, i));
                    }
                }
            } else {
                if (syst->fAcceptanceFolderNamesDown.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fAcceptanceFolderNamesDown, i));
                    }
                }
            }
        }
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptancePathSuffs,smp->fAcceptancePathSuffs), syst->fAcceptancePathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptancePathSuffs,smp->fAcceptancePathSuffs), syst->fAcceptancePathSuffsDown);
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceFileSuffs, smp->fAcceptanceFileSuffs), syst->fAcceptanceFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceFileSuffs, smp->fAcceptanceFileSuffs), syst->fAcceptanceFileSuffsDown);
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceNameSuffs, smp->fAcceptanceNameSuffs), syst->fAcceptanceNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceNameSuffs, smp->fAcceptanceNameSuffs), syst->fAcceptanceNameSuffsDown);
    } else {
        pathSuffs = Common::CombinePathSufs(reg->fAcceptancePathSuffs, smp->fAcceptancePathSuffs);
        fileSuffs = Common::CombinePathSufs(reg->fAcceptanceFileSuffs, smp->fAcceptanceFileSuffs);
        nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fAcceptanceNameSuffs, smp->fAcceptanceNameSuffs), fAcceptanceNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
std::vector<std::string> TRExFit::FullSelectionEffPaths(const Region* reg,
                                                        const UnfoldingSample* smp,
                                                        const UnfoldingSystematic* syst,
                                                        const bool isUp,
                                                        const bool isSubtract) const {
    // protection against nullptr
    if(!reg) {
        LOG(ERROR) << "Null pointer for Region.\n";
        exit(EXIT_FAILURE);
    }
    if(!smp) {
        LOG(ERROR) << "Null pointer for Sample.\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> fullPaths;
    std::vector<std::string> paths;
    std::vector<std::string> pathSuffs;
    std::vector<std::string> files;
    std::vector<std::string> fileSuffs;
    std::vector<std::string> names;
    std::vector<std::string> nameSuffs;
    // precendence:
    // 1. Systematic
    // 2. Sample
    // 3. Region
    // 4. Job
    if(syst){
        if(isUp) {
            if(syst->fSelectionEffPathsUp.size()  >0) paths = syst->fSelectionEffPathsUp;
            if(syst->fSelectionEffFilesUp.size()  >0) files = syst->fSelectionEffFilesUp;
            if(syst->fSelectionEffNamesUp.size()  >0) names = syst->fSelectionEffNamesUp;
        } else {
            if(syst->fSelectionEffPathsDown.size()>0) paths = syst->fSelectionEffPathsDown;
            if(syst->fSelectionEffFilesDown.size()>0) files = syst->fSelectionEffFilesDown;
            if(syst->fSelectionEffNamesDown.size()>0) names = syst->fSelectionEffNamesDown;
        }
        if (isSubtract) {
            if (isUp) {
                if (!syst->fSelectionEffPathsUpSubtractSample.empty()) paths = syst->fSelectionEffPathsUpSubtractSample;
                if (!syst->fSelectionEffFilesUpSubtractSample.empty()) files = syst->fSelectionEffFilesUpSubtractSample;
                if (!syst->fSelectionEffNamesUpSubtractSample.empty()) names = syst->fSelectionEffNamesUpSubtractSample;
            } else {
                if (!syst->fSelectionEffPathsDownSubtractSample.empty()) paths = syst->fSelectionEffPathsDownSubtractSample;
                if (!syst->fSelectionEffFilesDownSubtractSample.empty()) files = syst->fSelectionEffFilesDownSubtractSample;
                if (!syst->fSelectionEffNamesDownSubtractSample.empty()) names = syst->fSelectionEffNamesDownSubtractSample;
            }
        }
    }
    if(paths.size()==0 && smp->fSelectionEffPaths.size()>0) paths = smp->fSelectionEffPaths;
    if(files.size()==0 && smp->fSelectionEffFiles.size()>0) files = smp->fSelectionEffFiles;
    if(names.size()==0 && smp->fSelectionEffNames.size()>0) names = smp->fSelectionEffNames;

    if(paths.size()==0 && reg->fSelectionEffPaths.size()>0) paths = reg->fSelectionEffPaths;
    if(files.size()==0 && reg->fSelectionEffFiles.size()>0) files = reg->fSelectionEffFiles;
    if(names.size()==0 && reg->fSelectionEffNames.size()>0) names = reg->fSelectionEffNames;

    if(paths.size()==0 && fSelectionEffPaths.size()>0) paths = fSelectionEffPaths;
    if(files.size()==0 && fSelectionEffFiles.size()>0) files = fSelectionEffFiles;
    if(names.size()==0 && fSelectionEffNames.size()>0) names = fSelectionEffNames;

    // now we need to propagate the folder name updates
    if (syst) {

        auto broadcast = [&](auto& updates, std::size_t idx){
            if (updates.size() == names.size())
                return updates.at(idx);
            if (updates.size() == 1)
                return updates.at(0);
            return names.at(idx);
        };

        if (isSubtract) {
            if (isUp) {
                if (syst->fSelectionEffFolderNamesUpSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fSelectionEffFolderNamesUpSubtractSample, i));
                    }
                }
            } else {
                if (syst->fSelectionEffFolderNamesDownSubtractSample.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fSelectionEffFolderNamesDownSubtractSample, i));
                    }
                }
            }
        } else {
            if (isUp) {
                if (syst->fSelectionEffFolderNamesUp.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fSelectionEffFolderNamesUp, i));
                    }
                }
            } else {
                if (syst->fSelectionEffFolderNamesDown.size() > 0) {
                    for (std::size_t i = 0; i < names.size(); ++i) {
                        names.at(i) = Common::ReplaceFolderName(names.at(i), broadcast(syst->fSelectionEffFolderNamesDown, i));
                    }
                }
            }
        }
        if(isUp) pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffPathSuffs,smp->fSelectionEffPathSuffs), syst->fSelectionEffPathSuffsUp);
        else     pathSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffPathSuffs,smp->fSelectionEffPathSuffs), syst->fSelectionEffPathSuffsDown);
        if(isUp) fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffFileSuffs, smp->fSelectionEffFileSuffs), syst->fSelectionEffFileSuffsUp);
        else     fileSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffFileSuffs, smp->fSelectionEffFileSuffs), syst->fSelectionEffFileSuffsDown);
        if(isUp) nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffNameSuffs, smp->fSelectionEffNameSuffs), syst->fSelectionEffNameSuffsUp);
        else     nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffNameSuffs, smp->fSelectionEffNameSuffs), syst->fSelectionEffNameSuffsDown);
    } else {
        pathSuffs = Common::CombinePathSufs(reg->fSelectionEffPathSuffs, smp->fSelectionEffPathSuffs);
        fileSuffs = Common::CombinePathSufs(reg->fSelectionEffFileSuffs, smp->fSelectionEffFileSuffs);
        nameSuffs = Common::CombinePathSufs(Common::CombinePathSufs(reg->fSelectionEffNameSuffs, smp->fSelectionEffNameSuffs), fSelectionEffNamesNominal);
    }

    // And finally put everything together
    fullPaths = Common::CreatePathsList(paths, pathSuffs, files, fileSuffs, names, nameSuffs);
    return fullPaths;
}

//__________________________________________________________________________________
//
void TRExFit::PlotUnfold(TH1D* data,
                         TGraphAsymmErrors* total,
                         TGraphAsymmErrors* statOnly,
                         const Unfolding* unfolding,
                         const std::string& outputPath,
                         const YamlConverter& converter) {

    std::vector<std::unique_ptr<TH1> > mc;
    std::vector<std::string> legendNames;
    for (const auto& isample : fTruthSamples) {
        if (!isample->GetUseForPlotting()) continue;
        if (isample->fUnfoldingName != unfolding->fName) continue;
        auto hist = isample->GetHisto(unfolding);
        if (!hist) {
            LOG(ERROR) << "Histogram for the truth sample is nullptr\n";
            return;
        }
        hist->SetDirectory(nullptr);
        hist->Scale(fLumiScale);
        mc.emplace_back(std::move(hist));
        mc.back()->SetLineColor(isample->GetLineColor());
        mc.back()->SetLineStyle(isample->GetLineStyle());
        mc.back()->SetLineWidth(isample->GetLineWidth());
        legendNames.emplace_back(isample->GetTitle());
    }

    if (mc.empty()) {
        LOG(WARNING) << "No MC samples set for plotting. Will not create the final plots\n";
        return;
    }

    if(unfolding->fUnfoldNormXSec){
        for (auto& itruth : mc) {
            itruth->Scale(1./itruth->Integral());
        }
    }

    if (unfolding->fUnfoldingDivideByBinWidth) {
        Common::ScaleByBinWidth(data);
        Common::ScaleByBinWidth(total);
        if (statOnly) {
            Common::ScaleByBinWidth(statOnly);
        }
        for (auto& imc : mc) {
            Common::ScaleByBinWidth(imc.get());
        }
    }

    // divide by lumi if set, but only for absolute differential xsection
    if ((unfolding->fUnfoldingDivideByLumi > 0) && !(unfolding->fUnfoldNormXSec)) {
        data->Scale(1./unfolding->fUnfoldingDivideByLumi);
        Common::ScaleByConst(total, 1./unfolding->fUnfoldingDivideByLumi);
        if (statOnly) {
            Common::ScaleByConst(statOnly, 1./unfolding->fUnfoldingDivideByLumi);
        }
        for (auto& imc : mc) {
            imc->Scale(1./unfolding->fUnfoldingDivideByLumi);
        }
    }

    TCanvas c("","",600,600);
    TPad pad1("pad1","pad1",0.0, 0.3, 1.0, 1.00);
    TPad pad2("pad2","pad2", 0.0, 0.010, 1.0, 0.3);
    pad1.SetBottomMargin(0.001);
    pad1.SetBorderMode(0);
    pad2.SetBottomMargin(0.5);
    pad1.SetTicks(1,1);
    pad2.SetTicks(1,1);
    pad1.Draw();
    pad2.Draw();

    if (unfolding->fUnfoldingLogX) {
        pad1.SetLogx();
        pad2.SetLogx();
    }
    if (unfolding->fUnfoldingLogY) {
        pad1.SetLogy();
    }

    pad1.cd();
    total->SetMarkerStyle(20);
    total->SetMarkerSize(1.2);
    total->SetLineColor(kBlack);
    total->SetLineStyle(1);

    total->GetYaxis()->SetLabelSize(0.05);
    total->GetYaxis()->SetLabelFont(42);
    total->GetYaxis()->SetTitleFont(42);
    total->GetYaxis()->SetTitleSize(0.07);
    total->GetYaxis()->SetTitleOffset(1.1);

    std::unique_ptr<TH1> h_dummy(static_cast<TH1*>(mc[0]->Clone()));
    const double corr = unfolding->fUnfoldingScaleRangeY > 0 ? unfolding->fUnfoldingScaleRangeY : (unfolding->fUnfoldingLogY ? 1e6 : 1.5);
    h_dummy->GetYaxis()->SetRangeUser(unfolding->fUnfoldingMinRangeY, corr*h_dummy->GetMaximum());
    h_dummy->GetYaxis()->SetTitle(unfolding->fUnfoldingTitleY.c_str());
    h_dummy->GetYaxis()->SetTitleOffset(unfolding->fUnfoldingTitleOffsetY*h_dummy->GetYaxis()->GetTitleOffset());
    h_dummy->SetLineWidth(0);
    h_dummy->SetLineColor(kWhite);
    h_dummy->DrawClone("HIST");

    std::unique_ptr<TGraphAsymmErrors> error(static_cast<TGraphAsymmErrors*>(total->Clone()));
    error->SetFillStyle(1001);
    if (unfolding->fErrorBandColor > 0) {
        error->SetFillColor(unfolding->fErrorBandColor);
        error->SetLineColor(unfolding->fErrorBandColor);
    } else {
        error->SetFillColor(kGray+1);
        error->SetLineColor(kGray+1);
    }
    error->Draw("E2 same");

    std::unique_ptr<TGraphAsymmErrors> statError(nullptr);
    if (statOnly) {
        statError.reset(static_cast<TGraphAsymmErrors*>(statOnly->Clone()));
        statError->SetFillStyle(1001);
        if (unfolding->fStatErrorBandColor > 0) {
            statError->SetFillColor(unfolding->fStatErrorBandColor);
            statError->SetLineColor(unfolding->fStatErrorBandColor);
        } else {
            statError->SetFillColor(kGray);
            statError->SetLineColor(kGray);
        }
        statError->Draw("E2 same");
    }

    std::vector<double> chi2Vec{};
    std::vector<double> chi2NDFVec{};
    std::vector<double> probVec{};

    for (std::size_t i = 0; i < mc.size(); ++i) {
        mc.at(i)->Draw("HIST same");

        const double chi2 = GetUnfoldedChi2(error.get(), mc.at(i).get(), fFitResults->GetCorrelationMatrix(), unfolding);
        const int ndf     = unfolding->fUnfoldNormXSec ? (error->GetN() -1) : error->GetN();
        const double prob = ROOT::Math::chisquared_cdf_c(chi2, ndf);
        chi2Vec.emplace_back(chi2);
        chi2NDFVec.emplace_back(static_cast<double>(chi2/ndf));
        probVec.emplace_back(prob);

        LOG(INFO) << "MC: " << legendNames.at(i) << "\n";
        LOG(INFO) << "\tchi2: " << chi2 << "\n";
        LOG(INFO) << "\tndf: " << ndf << "\n";
        LOG(INFO) << "\tprob: " << prob << "\n";
    }

    total->Draw("pX SAME");

    float legH = (mc.size()+2)*0.07;
    TLegend leg(unfolding->fLegendXposition, 0.9-legH, unfolding->fLegendXposition+0.2, 0.9);
    if (unfolding->fAlternativeAsimovTruthSample == "") {
        leg.AddEntry(total, "Data", "p");
    } else {
        leg.AddEntry(total, "Pseudo-data", "p");
    }
    for (std::size_t imc = 0; imc < mc.size(); ++imc) {
        if (unfolding->fUnfoldingChi2Type == Unfolding::UnfoldingChi2Type::EMPTY) {
            leg.AddEntry(mc.at(imc).get(), legendNames.at(imc).c_str(), "l");
        } else if (unfolding->fUnfoldingChi2Type == Unfolding::UnfoldingChi2Type::CHI2) {
            leg.AddEntry(mc.at(imc).get(), Form("%s(#chi^{2}: %.1f)",legendNames.at(imc).c_str(), chi2Vec.at(imc)), "l");
        } else if (unfolding->fUnfoldingChi2Type == Unfolding::UnfoldingChi2Type::CHI2NDF) {
            leg.AddEntry(mc.at(imc).get(), Form("%s(#chi^{2}/NDF: %.1f)",legendNames.at(imc).c_str(), chi2NDFVec.at(imc)), "l");
        } else if (unfolding->fUnfoldingChi2Type == Unfolding::UnfoldingChi2Type::PROB) {
            leg.AddEntry(mc.at(imc).get(), Form("%s(prob: %.2f)",legendNames.at(imc).c_str(), probVec.at(imc)), "l");
        }
    }
    if (statError) {
        leg.AddEntry(statError.get(), "Stat. uncertainty","f");
    }
    leg.AddEntry(error.get(), "Total uncertainty","f");

    leg.SetFillColor(0);
    leg.SetLineColor(0);
    leg.SetBorderSize(0);
    leg.SetTextFont(gStyle->GetTextFont());
    if (unfolding->fLegendTextSize > 0) {
        leg.SetTextSize(unfolding->fLegendTextSize);
    } else {
        leg.SetTextSize(gStyle->GetTextSize());
    }
    leg.Draw("SAME");

    if (fPlotLabel != "none") TRExLabel(0.2,0.87, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str());
    myText(0.2,0.8,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    pad1.RedrawAxis();

    // Plot ratio
    pad2.cd();
    h_dummy->GetXaxis()->SetTitleOffset(unfolding->fUnfoldingTitleOffsetX*h_dummy->GetXaxis()->GetTitleOffset());
    h_dummy->GetYaxis()->SetTitleOffset(unfolding->fUnfoldingTitleOffsetY*0.55*h_dummy->GetYaxis()->GetTitleOffset());
    h_dummy->GetYaxis()->SetRangeUser(unfolding->fUnfoldingRatioYmin,unfolding->fUnfoldingRatioYmax);
    h_dummy->GetYaxis()->SetTitle("#frac{Prediction}{Data}");
    h_dummy->GetYaxis()->SetNdivisions(505);
    h_dummy->GetXaxis()->SetTitle(unfolding->fUnfoldingTitleX.c_str());
    h_dummy->Draw("HIST");

    // Plot the error band
    std::unique_ptr<TGraphAsymmErrors> band = Common::GetRatioBand(total, data);
    band->SetFillStyle(1001);
    if (unfolding->fErrorBandColor > 0) {
        band->SetFillColor(unfolding->fErrorBandColor);
        band->SetLineColor(unfolding->fErrorBandColor);
    } else {
        band->SetFillColor(kGray+1);
        band->SetLineColor(kGray+1);
    }
    band->SetMarkerStyle(0);
    band->SetLineWidth(2);
    band->Draw("E2 SAME");

    std::unique_ptr<TGraphAsymmErrors> bandStat;
    if (statOnly) {
        bandStat = Common::GetRatioBand(statOnly, data);
        bandStat->SetFillStyle(1001);
        if (unfolding->fStatErrorBandColor > 0) {
            bandStat->SetFillColor(unfolding->fStatErrorBandColor);
            bandStat->SetLineColor(unfolding->fStatErrorBandColor);
        } else {
            bandStat->SetFillColor(kGray);
            bandStat->SetLineColor(kGray);
        }
        bandStat->SetMarkerStyle(0);
        bandStat->SetLineWidth(2);
        bandStat->Draw("E2 SAME");
    }

    // Plot the theory predictions
    std::vector<std::unique_ptr<TH1D> > ratios;
    for (const auto& itruth : mc) {
        ratios.emplace_back(static_cast<TH1D*>(itruth->Clone()));
        ratios.back()->Divide(data);
        ratios.back()->Draw("HIST same");
    }

    const double min = mc.at(0)->GetXaxis()->GetXmin();
    const double max = mc.at(0)->GetXaxis()->GetXmax();

    TLine line(min, 1., max, 1.);
    line.SetLineStyle(unfolding->fRatioLineStyle);
    line.SetLineColor(unfolding->fRatioLineColor);
    line.SetLineWidth(unfolding->fRatioLineWidth);
    line.Draw("same");

    pad2.RedrawAxis();

    Common::SaveCanvasAs(c, outputPath+"/"+unfolding->fName+"_UnfoldedData");
    if (fHEPDataFormat) {
        converter.WriteUnfoldingHEPData(total, unfolding->fUnfoldingTitleX, outputPath, unfolding->fName, mc, legendNames);
    }

}

//__________________________________________________________________________________
//
void TRExFit::PlotMigrationResponse(const TH2* matrix,
                                    const bool isMigration,
                                    const std::string& regionName,
                                    const std::string& systematicName,
                                    const Unfolding* unfolding) const {

    gSystem->mkdir((fName+"/UnfoldingPlots/").c_str());

    std::unique_ptr<TH2D> m(static_cast<TH2D*>(matrix->Clone()));

    gStyle->SetPalette(87);
    if (isMigration) gStyle->SetPaintTextFormat("1.2f");
    else gStyle->SetPaintTextFormat("2.1f");
    TCanvas c("","",600,600);
    c.cd();

    for (int ibin = 1; ibin <= matrix->GetNbinsX(); ++ibin) {
        for (int jbin = 1; jbin <= matrix->GetNbinsY(); ++jbin) {
            if (m->GetBinContent(ibin, jbin) <= 0) {
                m->SetBinContent(ibin, jbin, 1e-6);
            }
        }
    }

    TPad pad("","",0,0.0,1,1);
    pad.SetRightMargin(0.15);
    pad.SetLeftMargin(0.15);
    pad.SetTopMargin(0.15);
    pad.SetBottomMargin(0.15);
    pad.Draw();
    pad.cd();

    if (unfolding->fMigrationLogX) {
        pad.SetLogx();
    }
    if (unfolding->fMigrationLogY) {
        pad.SetLogy();
    }

    double markerSize = matrix->GetNbinsX() < 10 ? 1. : 0.5;

    #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
    m->SetMarkerSize(0.75*1000);
    #else
    m->SetMarkerSize(markerSize);
    #endif
    m->SetMarkerColor(kBlack);
    m->GetXaxis()->SetTitle(unfolding->GetTitleRegionX(regionName).c_str());
    m->GetXaxis()->SetTitleOffset(unfolding->fMigrationTitleOffsetX * m->GetXaxis()->GetTitleOffset());
    m->GetYaxis()->SetTitleOffset(unfolding->fMigrationTitleOffsetY * m->GetYaxis()->GetTitleOffset());
    m->GetYaxis()->SetTitle(unfolding->GetTitleRegionY(regionName).c_str());
    m->GetXaxis()->SetNdivisions(505);
    m->GetZaxis()->SetTitleOffset(1.3);
    if (isMigration) {
        m->GetZaxis()->SetTitle("Migration");
        m->GetZaxis()->SetRangeUser(unfolding->fMigrationZmin,unfolding->fMigrationZmax);
    } else {
        m->GetZaxis()->SetTitle("Response");
        m->GetZaxis()->SetRangeUser(unfolding->fResponseZmin,unfolding->fResponseZmax);
    }

    c.SetGrid();
    if(unfolding->fMigrationText) m->Draw("COLZ TEXT");
    else m->Draw("COLZ");

    if (fPlotLabel != "none") TRExLabel(0.03,0.92, TRExFitter::EXPERIMENT_LABEL.c_str(), fPlotLabel.c_str());
    myText(0.68,0.92,1,Form("#sqrt{s} = %s, %s",fCmeLabel.c_str(),fLumiLabel.c_str()));

    c.RedrawAxis("g");

    if (systematicName != "") {
        gSystem->mkdir((fName+"/Systematics").c_str());
        gSystem->mkdir((fName+"/Systematics/"+systematicName).c_str());
    }
    const std::string tmp = isMigration ? "migration" : "response";
    const std::string name = systematicName == "" ? tmp + "_" + unfolding->fName + "_" + regionName : "/Systematics/" + systematicName+"/" + unfolding->fName + "_" + tmp + "_" + regionName;

    Common::SaveCanvasAs(c, fName+"/UnfoldingPlots/" + name);
    if (fHEPDataFormat) {
        // store hepdata
        const std::string hepDataName = unfolding->fName + "_" + regionName + "_" + systematicName;
        const bool trueHorizontal = (unfolding->fMatrixOrientation == FoldingManager::MATRIXORIENTATION::TRUTHONHORIZONTALAXIS);
        YamlConverter converter{};
        converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
        converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
        converter.WriteMigrationResponseHEPData(matrix, fName, !isMigration, hepDataName, regionName, trueHorizontal);
    }
}

//__________________________________________________________________________________
//
void TRExFit::RunForceShape() {
    for (const auto& ireg : fRegions) {
        for (const auto& ismp : fSamples) {
            if(Common::FindInStringVector(ismp->fRegions, ireg->fName) < 0) continue;
            std::shared_ptr<SampleHist> sh = ireg->GetSampleHist(ismp->fName);
            if(!sh) continue;
            for (const auto& isyst : ismp->fSystematics) {
                if (isyst->fForceShape == HistoTools::FORCESHAPETYPE::NOSHAPE) continue;
                if(isyst->fRegions.size()>0 && Common::FindInStringVector(isyst->fRegions,ireg->fName)<0  ) continue;
                if(isyst->fExclude.size()>0 && Common::FindInStringVector(isyst->fExclude,ireg->fName)>=0 ) continue;
                if(isyst->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(isyst->fExcludeRegionSample, ireg->fName, ismp->fName)>=0 ) continue;
                std::shared_ptr<SystematicHist> syh = sh->GetSystematic(isyst->fName);
                if(!syh) continue;

                HistoTools::ForceShape(syh->fHistUp.get(), sh->fHist.get(), isyst->fForceShape);
                // Symmetrise
                for (int ibin = 1; ibin <= sh->fHist->GetNbinsX(); ++ibin) {
                    const double diff = syh->fHistUp->GetBinContent(ibin) - sh->fHist->GetBinContent(ibin);
                    syh->fHistDown->SetBinContent(ibin, sh->fHist->GetBinContent(ibin) - diff);
                }
            }
        }
    }
}

//__________________________________________________________________________________
//
bool TRExFit::DoingMixedFitting() const {
    if (fFitNPValues.size() > 0) return false;
    if (fFitNPValuesFromFitResultsFile != "") return false;
    int dataType = -1;
    // loop on regions
    for (const auto& reg : fRegions) {
        if (reg->fRegionType == Region::VALIDATION) continue;
        if (dataType >= 0 && reg->fRegionDataType != dataType) return true;
        dataType = reg->fRegionDataType;
    }
    return false;
}

//__________________________________________________________________________________
//
std::vector < std:: string > TRExFit::ListRegionsToFit(const bool useFitRegions, int dataType) const {
    std::vector < std:: string > list;
    for (const auto& reg : fRegions) {
        if (reg->fRegionType == Region::VALIDATION) continue;
        if (dataType >= 0 && reg->fRegionDataType != dataType) continue;
        if (useFitRegions && fFitRegion == CRONLY && reg->fRegionType != Region::CONTROL) continue;
        list.emplace_back( reg->fName );
    }
    return list;
}

//__________________________________________________________________________________
//
std::map < std::string, int > TRExFit::MapRegionDataTypes(const std::vector<std::string>& regionList,bool isBlind) const {
    std::map < std::string, int > dataTypes;
    for (const auto& name : regionList) {
        if (isBlind) {
            dataTypes[name] = Region::ASIMOVDATA;
            continue;
        }
        const Region* reg = GetRegion(name);
        if (!reg) {
            LOG(ERROR) << "Trying to access a non-exisiting region: " << name << "...\n";
            exit(EXIT_FAILURE);
        }
        dataTypes[name] = reg->fRegionDataType;
    }
    return dataTypes;
}

//__________________________________________________________________________________
//
std::map < std::string, double > TRExFit::NPValuesFromFitResultsFile(const std::string& fitResultsFile) {
    std::map < std::string, double > npValues;
    ReadFitResults(fitResultsFile);
    for(const auto& inp : fFitResults->GetNuisanceParameters()) {
        // now we need to append some prefixes as we stripped the "gamma_" prefix when reading the txt
        // file in FitResults::ReadFromTXT, while the "alpha_" prefix is not even stored in the txt file
        auto nfIt = std::find_if(fNormFactors.begin(), fNormFactors.end(), [&inp](const std::shared_ptr<NormFactor>& nf){return nf->fName == inp.second->fName;});
        // dont modify NFs
        if (nfIt != fNormFactors.end()) {
            npValues[inp.second->fName] = inp.second->fFitValue;
            continue;
        }

        auto npIt = std::find_if(fSystematics.begin(), fSystematics.end(), [&inp](const std::shared_ptr<Systematic>& np){return np->fNuisanceParameter == inp.second->fName;});

        // if it is a NP we need to add "alpha_", if it is not, it is a gamma NP so we add "gamma_"
        const std::string name = ((npIt == fSystematics.end()) ? "gamma_" : "alpha_") + inp.second->fName;

        npValues[name] = inp.second->fFitValue;
    }

    return npValues;
}

//__________________________________________________________________________________
//
void TRExFit::FixUnfoldingExpressions() {
    if (fFitType != TRExFit::FitType::UNFOLDING) return;

    // Get truth histogram
    std::unique_ptr<TFile> input(TFile::Open((fName + "/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
    if (!input) {
        LOG(ERROR) << "Cannot read file from " << fName + "/UnfoldingHistograms/FoldedHistograms.root\n";
        exit(EXIT_FAILURE);
    }
    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        FixUnfoldingExpressions(fUnfolding.at(i).get(), input.get());
    }

    // process CorrelateUnfolding
    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        ProcessCorrelateUnfolding(fUnfolding.at(i).get(), input.get());
    }

    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        ProcessReplacementUnfolding(fUnfolding.at(i).get(), input.get());
    }

    input->Close();
}

//__________________________________________________________________________________
//
void TRExFit::FixUnfoldingExpressions(const Unfolding* unfolding, TFile* input) {
    if (!unfolding->fUnfoldNormXSec) return;

    std::unique_ptr<TH1> truth(dynamic_cast<TH1*>(input->Get((unfolding->fName + "_truth_distribution").c_str())));
    if (!truth) {
        LOG(ERROR) << "Cannot read the truth distribution\n";
        exit(EXIT_FAILURE);
    }
    truth->SetDirectory(nullptr);
    const double N = truth->Integral();
    // Get expression of last bin norm-factor
    const std::string NFname = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(unfolding->fUnfoldNormXSecBinN) + "_mu";
    for (const auto& inorm : fNormFactors) {
        if (inorm->fName==NFname) {
            TString expression(inorm->fExpression.first);
            for (int i = 0; i < unfolding->fNumberUnfoldingTruthBins; ++i) {
                // replace all "Nj/N" with truth->bincontent/truth->integral
                expression.ReplaceAll("N"+std::to_string(i+1)+"/N",std::to_string(truth->GetBinContent(i+1)/N));
            }
            inorm->fExpression.first = expression.Data();
            // title will contain the expression FIXME
            inorm->fTitle = inorm->fExpression.first;
            TRExFitter::SYSTMAP[inorm->fName] = inorm->fExpression.first;
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::ProcessCorrelateUnfolding(const Unfolding* unfolding, TFile* input) {
    // Get truth histogram
    std::unique_ptr<TH1> truth(dynamic_cast<TH1*>(input->Get((unfolding->fName + "_truth_distribution").c_str())));
    if (!truth) {
        LOG(ERROR) << "Cannot read the truth distribution\n";
        exit(EXIT_FAILURE);
    }
    truth->SetDirectory(nullptr);
    const double N = truth->Integral();

    for (const auto& inorm : fNormFactors) {
        if (inorm->fCorrelateUnfolding != unfolding->fName) continue;
        if (unfolding->fUnfoldNormXSec) {
            const std::string name = "TotalXsecOverTheory_" + unfolding->fName;
            const std::string expressionParam = name+"[1,"+std::to_string(unfolding->fUnfoldingResultMin)+","+std::to_string(unfolding->fUnfoldingResultMax)+"]";
            inorm->fExpression = std::make_pair(name, expressionParam);
        } else {
            std::string formula("");
            std::string params("");
            // formula is /sum_i \mu_i * N_i/N
            for (int i = 0; i < unfolding->fNumberUnfoldingTruthBins; ++i) {
                const double N_i = truth->GetBinContent(i+1);
                const std::string NFname = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
                formula += "(" + NFname + "*" + std::to_string(N_i) + "/" + std::to_string(N) + ")";
                params  += NFname+"[1,"+std::to_string(unfolding->fUnfoldingResultMin)+","+std::to_string(unfolding->fUnfoldingResultMax)+"]";
                if (i != unfolding->fNumberUnfoldingTruthBins - 1) {
                    formula += "+";
                    params += ",";
                }
            }
            inorm->fExpression = std::make_pair(formula, params);
        }
        LOG(DEBUG) << "Adding Expression: " << inorm->fExpression.first << inorm->fExpression.second << " to Unfolding: " << inorm->fName << "\n";
        TRExFitter::SYSTMAP[inorm->fName] = "Expression_" + inorm->fExpression.first;
        inorm->fNuisanceParameter = "Expression_" + inorm->fExpression.second;
        inorm->fTitle = "Expression_" + inorm->fExpression.second;
    }
}


//__________________________________________________________________________________
//
void TRExFit::ProcessReplacementUnfolding(const Unfolding* unfolding, TFile* input) {
    // Get truth histogram
    std::unique_ptr<TH1> truth(dynamic_cast<TH1*>(input->Get((unfolding->fName + "_truth_distribution").c_str())));
    if (!truth) {
        LOG(ERROR) << "Cannot read the truth distribution\n";
        exit(EXIT_FAILURE);
    }
    truth->SetDirectory(nullptr);
    static const std::vector<std::string> matches = {"@Yield", "@Width", "@Edge", "@Center", "@Integral"};

    auto Replacer = [&truth](const std::string& match, const int i) {
        if (i < 1 || i > truth->GetNbinsX()) {
            LOG(ERROR) << "Index for replacement is smaller than 1 or larger than NbinsX\n";
            exit(EXIT_FAILURE);
        }
        double result(0.);
        if (match == "@Yield") {
            result = truth->GetBinContent(i);
        } else if (match == "@Width") {
            result = truth->GetBinWidth(i);
        } else if (match == "@Edge") {
            result = truth->GetBinLowEdge(i);
        } else if (match == "@Center") {
            result = truth->GetBinCenter(i);
        } else if (match == "@Integral") {
            result = truth->Integral();
        }
        return std::to_string(result);
    };

    for (const auto& inorm : fNormFactors) {
        if (inorm->fExpression.first.empty()) continue;
        LOG(INFO) << "NF: " << inorm->fName << ", original Expression: " << inorm->fExpression.first << ", " << inorm->fExpression.second << "\n";
        for (const auto& imatch : matches) {
            size_t start_pos = 0;
            while((start_pos = inorm->fExpression.first.find(imatch, start_pos)) != std::string::npos) {
                std::string comb{""};
                size_t p = start_pos + imatch.length();
                while (p != std::string::npos && isdigit(inorm->fExpression.first[p])) {
                    comb += inorm->fExpression.first[p];
                    p++;
                }

                if (comb.empty() && imatch != "@Integral") {
                    LOG(ERROR) << "Matcher: " << imatch << " needs integer\n";
                    exit(EXIT_FAILURE);
                }
                const int fixed = comb.empty() ? 1 : stoi(comb);
                const std::string to = Replacer(imatch, fixed);
                inorm->fExpression.first.replace(start_pos, imatch.length() + comb.length(), to);
                start_pos += to.length();
            }
            size_t start_pos_second = 0;
            while((start_pos_second = inorm->fExpression.second.find(imatch, start_pos_second)) != std::string::npos) {
                std::string comb{""};
                size_t p = start_pos_second + imatch.length();
                while (p != std::string::npos && isdigit(inorm->fExpression.second[p])) {
                    comb += inorm->fExpression.second[p];
                    p++;
                }

                if (comb.empty() && imatch != "@Integral") {
                    LOG(ERROR) << "Matcher: " << imatch << " needs integer\n";
                    exit(EXIT_FAILURE);
                }
                const int fixed = comb.empty() ? 1 : stoi(comb);
                const std::string to = Replacer(imatch, fixed);
                inorm->fExpression.second.replace(start_pos_second, imatch.length() + comb.length(), to);
                start_pos_second += to.length();
            }

            // now we need to change the strings of type "NF1[10+10,-1,1],NF2[4+4.2,-1,1]" to
            // NF1[20,-1,1],NF2[8.2,-1,1], i.e. we need to evaluate the expressions
            std::vector<std::string> singleFactor = Common::Vectorize(inorm->fExpression.second,']');
            for (auto& istring : singleFactor) {
                istring += "]";
                if (istring.at(0) == ',') istring.erase(0,1);
                // now find the "[" and "]"
                std::size_t first = istring.find("[");
                std::size_t last = istring.find("]");
                std::string sub = istring.substr(first+1, last-first-1);
                // inner is a string of form "4+4,1,1"
                std::vector<std::string> singleFormula = Common::Vectorize(sub, ',');
                // now we need to replace each formula with the proper value
                for (auto& iformula : singleFormula) {
                    const std::string finalFormula = iformula += "+x";
                    TFormula tf("", finalFormula.c_str());
                    const double result = tf.Eval(0);
                    iformula = std::to_string(result);
                }
                // now we need to put it all back together
                std::string updated{""};
                for (const auto& iformula : singleFormula) {
                    updated += iformula + ",";
                }
                updated.pop_back(); // remove tha last ','

                // now replace it
                for( std::size_t pos = 0; ; pos += updated.length() ) {
                    // Locate the substring to replace
                    pos = istring.find( sub, pos );
                    if( pos == std::string::npos ) break;
                    // Replace by erasing and inserting
                    istring.erase( pos, sub.length() );
                    istring.insert( pos, updated );
                }
            }

            // now put everything together
            std::string res{""};
            for (const auto& i : singleFactor) {
                res += i + ",";
            }
            res.pop_back();
            inorm->fExpression.second = res;

        }
         // title will contain the expression FIXME
        inorm->fTitle = inorm->fExpression.first;
        TRExFitter::SYSTMAP[inorm->fName] = inorm->fExpression.first;
        TRExFitter::NPMAP[inorm->fName] = inorm->fExpression.second;
        LOG(INFO) << "\tReplaced Expression: " << inorm->fExpression.first << ", " << inorm->fExpression.second << "\n";
    }
}

//__________________________________________________________________________________
//
void TRExFit::RunLimitScan(xRooNLLVar::xRooHypoSpace& hs) const {
    if (fLimitType == LimitType::ASYMPTOTIC) {
        hs.scan("cls");
    } else if (fLimitType == LimitType::TOYS) {
        if (fLimitToysUsexRooFit == ToysUsexRooFit::AUTO) {
            hs.scan("cls toys");
        } else if (fLimitToysUsexRooFit == ToysUsexRooFit::TRUE) {
            TString sToysScanArg = Common::xRooFitToysScanArg(fLimitToysStepsSplusB, fLimitToysStepsB);
            if (sToysScanArg == "") {
                LOG(ERROR) << "Invalid scan arg\n";
                exit(EXIT_FAILURE);
            }
            TString scanArgs = TString::Format("cls toys=%s", sToysScanArg.Data());
            LOG(DEBUG) << "Calling xRooNLLVar::xRooHypoSpace::scan() with " << scanArgs << "\n";
            hs.scan(scanArgs,
                    fLimitToysScanSteps,
                    fLimitToysScanMin,
                    fLimitToysScanMax);
        } else {
            // This should never happen because ConfigReader does not allow this
            LOG(ERROR) << "LimitType is TOYS but ToysUsexRooFit is False!\n";
            exit(EXIT_FAILURE);
        }
    } else {
        LOG(ERROR) << "Unkown LimitType\n";
        exit(EXIT_FAILURE);
    }
}

//__________________________________________________________________________________
//
void TRExFit::RunLimit(RooWorkspace* ws, const std::string& dataName, const std::string& limitFile, const std::string& limitWSfile) const {
    const auto* constraints = FitUtils::GetExternalConstraints(ws, fNormFactors, fRegularizationType);
    xRooNode model(*ws);
    FitUtils::FixParametersInXRooFit(&model, fFitFixedNPs);

    LOG(DEBUG) << "Setting Hesse strategy to: " << fLimitFitStrategy << "\n";

    auto nll = model.nll(dataName.c_str());
    nll.SetOption(RooFit::Optimize(kTRUE));
    nll.SetOption(RooFit::NumCPU(fCPU, RooFit::Hybrid));
    if (constraints) {
        nll.SetOption(RooFit::ExternalConstraints(*constraints));
    }
    if (fLimitUseAutoDiff) {
        LOG(INFO) << "Will use automatic differentiation for the limit estimate\n";
        nll.SetOption(RooFit::EvalBackend("codegen"));
    }
    nll.fitConfigOptions()->SetValue("HesseStrategy",fLimitFitStrategy);
    if (TRExFitter::DEBUGLEVEL >= 2) {
        nll.reinitialize();
        nll.Print();
    }

    auto hs = nll.hypoSpace(fPOIforLimit.c_str());
    TFile outWS(limitWSfile.c_str(), "RECREATE");

    RunLimitScan(hs);

    std::unique_ptr<RooStats::HypoTestInverterResult> result(hs.result());
    result->SetName("clsLimits");
    model.ws()->import(*result);

    std::unique_ptr<xRooNLLVar::xRooHypoSpace> hsInjected(nullptr);
    std::unique_ptr<RooStats::HypoTestInverterResult> injectedResult(nullptr);

    if (fSignalInjection && fSignalInjectionValue > 0) {
       LOG(INFO) << "Injecting signal into the limit estimate: " << fSignalInjectionValue << "\n";
       model.poi().at(fPOIforLimit.c_str())->SetContent(fSignalInjectionValue);

       hsInjected = std::make_unique<xRooNLLVar::xRooHypoSpace>(model.nll().hypoSpace(fPOIforLimit.c_str()));
       RunLimitScan(*hsInjected); // get the observed limit
       injectedResult.reset(hsInjected->result());
       injectedResult->SetName("CLSLimits_injected");
       model.ws()->import(*injectedResult);
    }
    Common::ProcessLimitOutput(hs, hsInjected.get(), (fWorkspaceFileName != "")? "from a custom file" : "", fLimitIsBlind, limitFile, fLimitParamName, fLimitParamValue, fLimitType == LimitType::ASYMPTOTIC);
    outWS.Close();
    model.SaveAs(limitWSfile.c_str(), "UPDATE");

    if (fLimitPlot && fLimitType == TRExFit::LimitType::TOYS) {
       TCanvas c;
       hs.Draw();
       Common::SaveCanvasAs(c, fName + "/Limits/Toys/LimitToysResults" + fLimitToysSuffix);
    }

}

//__________________________________________________________________________________
//
void TRExFit::RunLimitToys(RooAbsData* data, RooWorkspace* ws) const {

    // Set the likelihood optimisation as it takes half-infinite time anyway
    FitUtils::SetBinnedLikelihoodOptimisation(ws);
    if (fIntCode != 4) {
        FitUtils::ChangeInterpolationCode(ws, fIntCode);
    }

    gSystem->mkdir((fName + "/Limits").c_str());

    LimitToys toys{};
    toys.SetNToys(fLimitToysStepsSplusB, fLimitToysStepsB);
    toys.SetLimit(fLimitsConfidence);
    toys.SetScan(fLimitToysScanSteps, fLimitToysScanMin, fLimitToysScanMax);
    toys.SetPlot(fLimitPlot);
    toys.SetFile(fLimitFile);
    toys.SetSeed(fLimitToysSeed);
    toys.SetOutputPath(fName + "/Limits");
    toys.SetLimitSuffix(fLimitToysSuffix);

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read ModelConfig\n";
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<RooStats::ModelConfig> mcBonly(static_cast<RooStats::ModelConfig*>(mc->Clone("BonlyModel")));
    RooRealVar* poi = static_cast<RooRealVar*>(mcBonly->GetParametersOfInterest()->first());
    poi->setVal(0);

    toys.RunToys(data, mc, mcBonly.get());
}

//__________________________________________________________________________________
//
void TRExFit::RunSignificance(RooWorkspace* ws, const std::string& dataName, const std::string& significanceFile, const std::string& significanceWSfile) const {

    const auto* constraints = FitUtils::GetExternalConstraints(ws, fNormFactors, fRegularizationType);

    xRooNode model(*ws);
    FitUtils::FixParametersInXRooFit(&model, fFitFixedNPs);

    auto nll = model.nll(dataName.c_str());
    nll.SetOption(RooFit::Optimize(kTRUE));
    nll.SetOption(RooFit::NumCPU(fCPU, RooFit::Hybrid));
    if (constraints) {
        nll.SetOption(RooFit::ExternalConstraints(*constraints));
    }
    if (fSignificanceUseAutoDiff) {
        LOG(INFO) << "Will use automatic differentiation for the significance estimate\n";
        nll.SetOption(RooFit::EvalBackend("codegen"));
    }
    nll.fitConfigOptions()->SetValue("HesseStrategy",fSignificanceFitStrategy);
    if (TRExFitter::DEBUGLEVEL >= 2) {
        nll.reinitialize();
        nll.Print();
    }
    auto hs = nll.hypoSpace(fPOIforSig.c_str(), xRooFit::Asymptotics::Uncapped, 1);
    TFile outWS(significanceWSfile.c_str(), "RECREATE");

    RunSignificanceScan(hs);

    std::unique_ptr<RooStats::HypoTestInverterResult> result(hs.result());
    const std::pair<double, double> pvalueWithErrorExpected = hs[0].pNull_asymp(0);
    const std::pair<double, double> pvalueWithErrorObserved  = hs[0].pNull_asymp();
    std::pair<double, double> pvalueWithErrorInjected = {-99, -99};
    model.ws()->import(*result);

    std::unique_ptr<xRooNLLVar::xRooHypoSpace> hsInjected(nullptr);
    std::unique_ptr<RooStats::HypoTestInverterResult> injectedResult(nullptr);
    if (fSignificanceDoInjection && fSignificancePOIAsimov > 0) {
        LOG(INFO) << "Injecting signal into the significance estimate: " << fSignificancePOIAsimov << "\n";
        model.poi().at(fPOIforSig.c_str())->SetContent(fSignificancePOIAsimov);

        hsInjected = std::make_unique<xRooNLLVar::xRooHypoSpace>(
            model.nll().hypoSpace(fPOIforSig.c_str(), xRooFit::Asymptotics::Uncapped, fSignificancePOIAsimov));
        RunSignificanceScan(*hsInjected);
        pvalueWithErrorInjected = (*hsInjected)[0].pNull_asymp();
        injectedResult.reset(hsInjected->result());
        injectedResult->SetName(
            TString::Format("%s_%s=%g",injectedResult->GetName(),fPOIforSig.c_str(),fSignificancePOIAsimov));
        model.ws()->import(*injectedResult);
    }

    if (Common::ProcessSignificanceOutput(pvalueWithErrorObserved, pvalueWithErrorExpected, pvalueWithErrorInjected,
                                          (fWorkspaceFileName != "")? "from a custom file" : "",
                                          fSignificanceIsBlind, significanceFile)) {
        LOG(WARNING) << "NaN value discovered, printing more information\n";
        hs.Print();
        if (hsInjected) {
            LOG(WARNING) << "Injected significance information\n";
            hsInjected->Print();
        }
    }

    outWS.Close();
    model.SaveAs(significanceWSfile.c_str(), "UPDATE");

    if (fSignificancePlot && fSignificanceType == TRExFit::SignificanceType::TOYS) {
        TCanvas c;
        hs.Draw();
        Common::SaveCanvasAs(c, fName + "/Significances/Toys/SignificanceToysResults");
     }

}

//__________________________________________________________________________________
//
void TRExFit::RunSignificanceScan(xRooNLLVar::xRooHypoSpace& hs) const {
    if (fSignificanceType == SignificanceType::ASYMPTOTIC) {
        hs.scan("pnull", 1, 0, 0);
    } else if (fSignificanceType == SignificanceType::TOYS) {
        if (fSignificanceToysUsexRooFit) {
            TString sToysScanArg = Common::xRooFitToysScanArg(fSignificanceToysStepsSplusB, fSignificanceToysStepsB);
            if (sToysScanArg == "") {
                LOG(ERROR) << "Invalid scan arg\n";
                exit(EXIT_FAILURE);
            }
            TString scanArgs = TString::Format("pnull toys=%s", sToysScanArg.Data());
            hs.scan(scanArgs, 1, 0, 0);
        } else {
            // This should never happen because ConfigReader does not allow this
            LOG(ERROR) << "SignificanceType is TOYS but ToysUsexRooFit is False!\n";
            exit(EXIT_FAILURE);
        }
    } else {
        LOG(ERROR) << "Unkown SignificanceType\n";
        exit(EXIT_FAILURE);
    }
}


//__________________________________________________________________________________
//
void TRExFit::RunSignificanceToys(RooAbsData* data, RooWorkspace* ws) const {

    // Set the likelihood optimisation as it takes half-infinite time anyway
    FitUtils::SetBinnedLikelihoodOptimisation(ws);
    if (fIntCode != 4) {
        FitUtils::ChangeInterpolationCode(ws, fIntCode);
    }

    gSystem->mkdir((fName + "/Significance").c_str());

    SignificanceToys toys{};
    toys.SetNtoys(fSignificanceToysStepsSplusB, fSignificanceToysStepsB);
    toys.SetPlot(fSignificancePlot);
    toys.SetToysSeed(fSignificanceToysSeed);
    toys.SetOutputPath(fName + "/Significance");
    toys.SetSignificanceSuffix("_"+std::to_string(fSignificanceToysSeed));

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot read ModelConfig\n";
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<RooStats::ModelConfig> mcBonly(static_cast<RooStats::ModelConfig*>(mc->Clone("BonlyModel")));
    RooRealVar* poi = static_cast<RooRealVar*>(mcBonly->GetParametersOfInterest()->first());
    poi->setVal(0);

    toys.RunToys(data, mc, mcBonly.get());
}

//__________________________________________________________________________________
//
void TRExFit::ProcessEFTInputs(bool overwrite) {


    std::map<std::string,std::unique_ptr<EFTProcessor> > EFTProcs;

    // Cycle through all EFT samples to get map of global EFT parameters, SM References and associated EFT Samples
    for (const auto& isample : fSamples) {

        // This is not a SM Reference sample
        if(isample->fType == Sample::SampleType::EFT && isample->fEFTSMReference != "NONE"){

            LOG(DEBUG) << "--------------------------------------------------------------------\n";
            LOG(DEBUG) << "EFT things happening here!\n";
            const std::string EFTParam = isample->fEFTParam;
            const std::string EFTName = isample->fName;
            const std::string EFTValue = isample->fEFTValue;
            const std::string EFTTitle = isample->fEFTTitle;
            const std::string SMRef = isample->fEFTSMReference;

            //@TODO: Add map of coeff values and file names!

            // Add global EFT parameter to the map if not already there
            auto itr = EFTProcs.find(SMRef);
            if(itr == EFTProcs.end()){
                LOG(DEBUG) << "Found new EFT sample: " << SMRef << "\n";
                std::unique_ptr<EFTProcessor> tmp_EFTProc = std::make_unique<EFTProcessor>(SMRef,"");
                itr = (EFTProcs.insert(std::make_pair(SMRef,std::move(tmp_EFTProc)))).first;
            }

            itr->second->fProduceShapeFactorParametrization = !fEFTConfig.GetSplitSamplesPerBin();

            LOG(DEBUG) << "Param: " << EFTParam << "\n";
            LOG(DEBUG) << "Title: " << EFTTitle << "\n";
            LOG(DEBUG) << "Name:  " << EFTName << "\n";
            LOG(DEBUG) << "Value: " << EFTValue << "\n";
            LOG(DEBUG) << "SMRef: " << SMRef << "\n";

            // Add new parameter combinations and their titles to EFTProcessor
            if (std::find(itr->second->fParams.begin(), itr->second->fParams.end(), EFTParam) == itr->second->fParams.end()) {
                itr->second->fParams.emplace_back(EFTParam);
            }

            // Split values to individual coeffs
            std::vector<std::string> coefftitles = Common::Vectorize(EFTTitle,',',true,false);
            for (const auto& coefftitle : coefftitles) {
                // Split into coeff and value
                std::vector<std::string> coeff_title = Common::Vectorize(coefftitle,'=',true,false);
                if (coeff_title.size() != 2) {
                    LOG(ERROR) << "Invalid EFTTitle setting: " << coefftitle << "\n";
                    exit(EXIT_FAILURE);
                }
                std::string coeff=coeff_title[0];
                std::string title=coeff_title[1];

                if(itr->second->fTitleMap.find(coeff) == itr->second->fTitleMap.end()){
                    itr->second->fTitleMap[coeff] = title;
                }
            }

            // Creat EFT Variation object
            EFTProcessor::EFTVariation tmp_EFTVar;
            tmp_EFTVar.fName = EFTName;
            tmp_EFTVar.fParam = EFTParam;
            tmp_EFTVar.fTitle = EFTTitle;
            tmp_EFTVar.fValueString = EFTValue;

            // Split parametrisation into individual terms
            std::vector<std::string> terms = Common::Vectorize(isample->fEFTParam,'+');
            for(const auto& term : terms) {
                // Split the powers in to coeffient and powers
                std::vector<std::string> coeff_power = Common::Vectorize(term,'^');
                std::vector<std::string> coeffs = Common::Vectorize(coeff_power.at(0),'*');
                std::string power = "1";
                if (coeff_power.size() > 1) power = coeff_power.at(1);

                std::vector<int> tmp_powers;
                tmp_powers.resize(itr->second->fOperators.size());

                for(const auto& coeff : coeffs) {
                    // If coeff/operator value not already present add it, otherwise just get index
                    auto op_it = std::find(itr->second->fOperators.begin(), itr->second->fOperators.end(), coeff);
                    int op_index = -1;
                    if(op_it == itr->second->fOperators.end()){
                        LOG(DEBUG) << "Found new operator: " << coeff << "\n";
                        itr->second->fOperators.emplace_back(coeff);
                        op_index = itr->second->fOperators.size() - 1;

                        // Resize all existing vectors to allow space for the new operator powers
                        itr->second->Resize();

                        tmp_powers.resize(itr->second->fOperators.size());
                    } else {
                        op_index = op_it - itr->second->fOperators.begin();
                    }
                    LOG(DEBUG) << TString::Format("coeff: %s (%i) power: %s",coeff.c_str(),op_index,power.c_str()).Data() << "\n";

                    tmp_powers[op_index] += Common::convertStoNum<int>(power);
                }

                // Add new term if not already present
                if (!std::count(itr->second->fPowers.begin(), itr->second->fPowers.end(), tmp_powers)) {
                    itr->second->fPowers.emplace_back(tmp_powers);
                    LOG(DEBUG) << "Found new term: " << term << "\n";
                }
            }

            itr->second->Resize();

            // Recording coefficient--value points from all samples

            // Split values to individual coeffs
            std::vector<std::string> coeffvalues = Common::Vectorize(isample->fEFTValue,',');
            for(const auto& coeffvalue : coeffvalues) {
                // Split into coeff and value
                std::vector<std::string> coeff_value = Common::Vectorize(coeffvalue,'=');
                if (coeff_value.size() != 2) {
                    LOG(ERROR) << "Invalid EFTValue setting: " << coeffvalue << "\n";
                    exit(EXIT_FAILURE);
                }
                std::string coeff = coeff_value.at(0);
                std::string value = coeff_value.at(1);
                LOG(DEBUG) << "Found new coeff value setting: " << coeff << " = " << value << "\n";

                // If coeff found store the Sample name and value
                auto op_it = std::find(itr->second->fOperators.begin(), itr->second->fOperators.end(), coeff);
                if (op_it == itr->second->fOperators.end()) {
                    LOG(ERROR) << "Couldn't find matching EFT coefficient for setting " << coeffvalue << "\n";
                    exit(EXIT_FAILURE);
                } else {
                    // Add EFT sample to EFTValues map if not already there
                    tmp_EFTVar.fValues.resize(itr->second->fOperators.size(),0.);
                    tmp_EFTVar.fValues[std::distance(itr->second->fOperators.begin(), op_it)] = Common::convertStoNum<double>(value);
                }
            }

            itr->second->fEFTVariations.emplace_back(tmp_EFTVar);
            itr->second->Resize();
            LOG(DEBUG) << "--------------------------------------------------------------------\n";
        }
    }

    for (const auto& EFTProc : EFTProcs) {
        EFTProc.second->PrintCoeffMap();
        EFTProc.second->PrintValuesMap();
    }

    if (!fEFTConfig.GetSplitSamplesPerBin()) {
        const std::string allFileName = fName + "/EFT/SF_ALL_EFT_Fit_results.txt";
        const std::string linearFileName = fName + "/EFT/SF_LINEAR_EFT_Fit_results.txt";
        std::remove(allFileName.c_str());
        std::remove(linearFileName.c_str());
    }

    // Go back through all global EFT params found to perform fits and add NFs
    for (const auto& imap : EFTProcs) {
        LOG(DEBUG) << " ==== " << imap.first << " ====\n";

        // Extract EFT inputs, plot them and fit
        const std::string fileName = fName + "/EFT/EFT_Fit_results_" + imap.first + ".txt";
        LOG(DEBUG) << " Trying to open " << fileName << "...\n";
        std::ifstream in(fileName.c_str());
        bool found = in.good();
        in.close();

        if (!fEFTConfig.GetSplitSamplesPerBin()) {
            imap.second->DrawEFTInputs(fRegions, fEFTConfig);
            imap.second->FitEFTInputs(fRegions, fName + "/EFT", imap.first.c_str(), fEFTConfig);
            continue;
        }

        // this is needed only for the old style of processing
        if (!found || overwrite) { // Rerun plotting and fitting from scratch
            if (found){ // file doesnt exist
                LOG(DEBUG) << " Running EFT plotting and fitting for first time (" << fileName << " not found)\n";
            } else {
                LOG(DEBUG) << " Overwriting existing EFT plotting and fitting (" << fileName << " already exists)\n";
            }
            imap.second->DrawEFTInputs(fRegions, fEFTConfig);
            imap.second->FitEFTInputs(fRegions, fName + "/EFT", imap.first.c_str(), fEFTConfig);
        } else {
            LOG(DEBUG) << " Reading EFT fit results from file (" << fileName << " not found)\n";
            imap.second->SetEFTOrder(fEFTConfig.GetEFTOrder());
            imap.second->ReadEFTFitResults(fRegions, fileName);
        }
        // Take fit results and apply to mu NormFactors
        imap.second->SetEFTOrder(fEFTConfig.GetEFTOrder()); // again for case of not yet done
        imap.second->ApplyMuFactExpressions(fRegions,fNormFactors);
    }

}
//__________________________________________________________________________________
//
double TRExFit::GetUnfoldedChi2(const TGraphAsymmErrors* error, const TH1* mc, const CorrelationMatrix* matrix, const Unfolding* unfolding) const {

    if (!matrix) {
        LOG(WARNING) << "Correlation matrix is nullptr, returning -1\n";
        return -1;
    }

    const int n = unfolding->fUnfoldNormXSec ? error->GetN() - 1 : error->GetN();
    if (n != (unfolding->fUnfoldNormXSec ? mc->GetNbinsX() - 1 : mc->GetNbinsX())) {
        LOG(WARNING) << "Number of bins do not match, returning -1\n";
        return -1;
    }

    // Fill the covariance
    TMatrixD cov(n,n);
    const int loopRange = unfolding->fUnfoldNormXSec ? (n + 1) : n;
    int indexI = 0;
    for (int i = 0; i < loopRange; ++i) {
        if (unfolding->fUnfoldNormXSec && unfolding->fUnfoldNormXSecBinN == i+1) continue;
        const double sigma_i = error->GetErrorY(i);
        const std::string name_i = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
        int indexJ = 0;
        for (int j = 0; j < loopRange; ++j) {
            if (unfolding->fUnfoldNormXSec && unfolding->fUnfoldNormXSecBinN == j+1) continue;
            const double sigma_j = error->GetErrorY(j);
            const std::string name_j = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(j+1) + "_mu";
            const double corr = matrix->GetCorrelation(name_i,name_j);
            const double value = (indexI == indexJ) ? sigma_i * sigma_i : sigma_i * corr * sigma_j;
            cov[indexI][indexJ] = value;
            ++indexJ;
        }
        ++indexI;
    }

    LOG(DEBUG) << "Covariance matrix\n";
    if(TRExFitter::DEBUGLEVEL > 1) cov.Print();
    cov.Invert();
    LOG(DEBUG) << "Inverted covariance matrix\n";
    if(TRExFitter::DEBUGLEVEL > 1) cov.Print();

    // calculate the chi2
    double chi2(0.);
    indexI = 0;
    for (int i = 0; i < loopRange; ++i) {
        if (unfolding->fUnfoldNormXSec && unfolding->fUnfoldNormXSecBinN == i+1) continue;
        double x,y;
        error->GetPoint(i, x, y);
        const double data_i = y;
        const double mc_i = mc->GetBinContent(i+1);
        int indexJ = 0;
        for (int j = 0; j < loopRange; ++j) {
            if (unfolding->fUnfoldNormXSec && unfolding->fUnfoldNormXSecBinN == j+1) continue;
            error->GetPoint(j, x, y);
            const double data_j = y;
            const double mc_j = mc->GetBinContent(j+1);

            chi2 += (mc_i - data_i) * cov[indexI][indexJ] * (mc_j - data_j);
            ++indexJ;
        }
        ++indexI;
    }

    return chi2;
}

//__________________________________________________________________________________
//
const Unfolding* TRExFit::GetCorrespondingUnfolding(const UnfoldingSample* sample) const {
    if (!sample) return nullptr;

    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        if (fUnfolding.at(i)->fName == sample->fUnfoldingName) return fUnfolding.at(i).get();
    }

    return nullptr;
}

//__________________________________________________________________________________
//
const Unfolding* TRExFit::GetCorrespondingUnfolding(const UnfoldingSystematic* syst) const {
    if (!syst) return nullptr;

    for (std::size_t i = 0; i < fUnfolding.size(); ++i) {
        if (fUnfolding.at(i)->fName == syst->fUnfoldingName) return fUnfolding.at(i).get();
    }

    return nullptr;
}

//__________________________________________________________________________________
//
void TRExFit::PlotCovarianceMatrix(const Unfolding* unfolding,
                                   const UnfoldingResult& unfolded,
                                   const std::string& outputPath) const {

    gSystem->mkdir((outputPath + "/UnfoldingPlots/").c_str());

    const int n = unfolding->fNumberUnfoldingTruthBins;
    const std::vector<UnfoldingResult::FitValue>& results = unfolded.GetFittedResults();
    const std::string unfoldingName = unfolding->fName;

    YamlConverter converter{};
    converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
    converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
    std::vector<std::vector<double> > covarianceMatrix(n, std::vector<double>(n));
    std::vector<std::string> names;
    TH2D histo("","", n, 0.5, n+0.5, n, 0.5, n+0.5);

    for (int i = 0; i < n; ++i) {
        const std::string name_i = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
        for (int j = 0; j < n; ++j) {
            const std::string name_j = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(j+1) + "_mu";
            const double corr = (i == j) ? 1. : fFitResults->GetCorrelationMatrix()->GetCorrelation(name_i, name_j);
            const double up_i   = results.at(i).up;
            const double down_i = results.at(i).down;
            const double up_j   = results.at(j).up;
            const double down_j = results.at(j).down;
            const double sigma_i = 0.5 * (std::abs(up_i) + std::abs(down_i));
            const double sigma_j = 0.5 * (std::abs(up_j) + std::abs(down_j));

            // this is covariance i,j
            const double cov = sigma_i * corr * sigma_j;
            histo.SetBinContent(i+1,n-j, cov);
            covarianceMatrix[i][j] = cov;
        }
        names.emplace_back(unfolding->fName + " bin " + std::to_string(i+1));
        histo.GetXaxis()->SetBinLabel(i+1, (unfolding->fName + " bin " + std::to_string(i+1)).c_str());
        histo.GetYaxis()->SetBinLabel(n-i, (unfolding->fName + " bin " + std::to_string(i+1)).c_str());
    }

    converter.WriteCorrelation(names, covarianceMatrix, outputPath + "/UnfoldingPlots/", true, "CovarianceMatrix_"+unfoldingName);
    if (fHEPDataFormat) {
        converter.WriteCorrelationHEPData(names, covarianceMatrix, outputPath, true, "CovarianceMatrix_"+unfoldingName);
    }

    TCanvas c1("","",0.,0.,1200,1200);
    gStyle->SetPalette(87);
    #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
    histo.SetMarkerSize(0.75*1000);
    #else
    histo.SetMarkerSize(1.0);
    #endif
    gStyle->SetPaintTextFormat("1.1e");
    c1.SetBottomMargin(0.2);
    c1.SetLeftMargin(0.2);

    histo.GetXaxis()->LabelsOption("v");
    histo.GetXaxis()->SetLabelSize( histo.GetXaxis()->GetLabelSize()*0.75 );
    histo.GetYaxis()->SetLabelSize( histo.GetYaxis()->GetLabelSize()*0.75 );
    c1.SetTickx(0);
    c1.SetTicky(0);
    histo.GetYaxis()->SetTickLength(0);
    histo.GetXaxis()->SetTickLength(0);
    histo.Draw("col TEXT");

    Common::SaveCanvasAs(c1, outputPath + "/UnfoldingPlots/CovarianceMatrix_"+unfoldingName);
}

//__________________________________________________________________________________
//
std::pair<std::unique_ptr<RooWorkspace>, std::unique_ptr<RooDataSet> > TRExFit::PrepareMixedDataset(const TRExFit::WorkspaceType type) {
    std::unique_ptr<RooDataSet> data(nullptr);
    //
    // Fills a vector of regions to consider for fit
    //
    bool blind(false);
    switch (type) {
        case (TRExFit::WorkspaceType::FIT):          blind = fFitIsBlind; break;
        case (TRExFit::WorkspaceType::SIGNIFICANCE): blind = fSignificanceIsBlind; break;
        case (TRExFit::WorkspaceType::LIMIT):        blind = fLimitIsBlind; break;
        default:
            LOG(ERROR) << "Unknown workspace type\n";
            exit(EXIT_FAILURE);
    }
    const bool cronlyFlag = fFitRegion == TRExFit::FitRegion::CRONLY;
    if (type == TRExFit::WorkspaceType::SIGNIFICANCE || type == TRExFit::WorkspaceType::LIMIT) {
        if (cronlyFlag) {
            LOG(WARNING) << "FitRegion set to CRONLY and running limit or significance. Only CRs will be used. Check that this is what you want\n";
        }
    }
    const std::vector <std::string>& regions = ListRegionsToFit(cronlyFlag);
    std::map <std::string, int> regionsDataType = MapRegionDataTypes(regions,blind);
    std::map <std::string, double> npValues;
    //
    // flag if mixed Data / Asimov required
    // if mixed limit, perform a first fit on the regions with data only
    const bool isMixedFit = DoingMixedFitting();
    if(isMixedFit && !blind) {
        const std::vector<std::string>& regionsForFit = ListRegionsToFit(cronlyFlag, Region::REALDATA);
        //
        // Creates a combined workspace with the regions to be used *in the fit*
        //
        LOG(INFO) << "Creating ws for regions with real data only...\n";
        std::unique_ptr<RooWorkspace> ws = PerformWorkspaceCombinationxRooFit(regionsForFit);
        if (!ws) {
            LOG(ERROR) << "Cannot retrieve the workspace, exiting!\n";
            exit(EXIT_FAILURE);
        }
        //
        // Calls the PerformFit() function to actually do the fit
        //
        LOG(INFO) << "Performing a fit in regions with real data only...\n";

        const bool currentErrorDecompositionFlag = fErrorDecomposition;
        // temporarily disable error decompositon and minos
        fErrorDecomposition = false;
        const auto minosCopy = fVarNameMinos;
        fVarNameMinos.clear();

        npValues = PerformFit(ws.get(), data.get(), FitType::BONLY, false);

        fErrorDecomposition = currentErrorDecompositionFlag;
        fVarNameMinos = minosCopy;
        LOG(INFO) << "Now will use the fit results to create the Asimov in the regions without real data!\n";
    }
    else{
        npValues = fFitNPValues;
    }

    //
    // Create the final asimov dataset
    //
    std::unique_ptr<RooWorkspace> ws = PerformWorkspaceCombinationxRooFit(regions);
    if (!ws) {
        LOG(ERROR) << "Cannot retrieve the workspace, exiting!\n";
        exit(EXIT_FAILURE);
    }
    std::map<std::string,double> poiValues;
    std::string poi{};
    double poiVal{};
    if (type == TRExFit::WorkspaceType::FIT) {
        for(const auto& ipoi : fPOIs){
            // if POI found in npValues, use that value
            if(npValues.find(ipoi)!=npValues.end()){
                poiValues[ipoi] = npValues[ipoi];
            }
            // else, if found in POIasimov, use that
            else if(fFitPOIAsimov.find(ipoi)!=fFitPOIAsimov.end()){
                poiValues[ipoi] = fFitPOIAsimov.find(ipoi)->second;
            }
        }
    } else {
        poi    = (type == TRExFit::WorkspaceType::SIGNIFICANCE) ? fPOIforSig : fPOIforLimit;
        poiVal = (type == TRExFit::WorkspaceType::SIGNIFICANCE) ? fSignificancePOIAsimov : 0.;
    }

    if (type != TRExFit::WorkspaceType::FIT) {
        if(npValues.find(poi) != npValues.end()){
            poiValues[poi] = npValues[poi];
        }
        else{
            poiValues[poi] = poiVal;
        }
    }
    data = std::unique_ptr<RooDataSet>(DumpData(ws.get(), regionsDataType, npValues, poiValues ));

    //
    // Set the global observables to the specified NPValues if the corresponding flag is TRUE
    //
    if(fInjectGlobalObservables && !fFitNPValues.empty()) {
        FitUtils::InjectGlobalObservables(ws.get(), fFitNPValues);
    }

    return std::make_pair(std::move(ws), std::move(data));
}

//__________________________________________________________________________________
//
void TRExFit::PrepareTemplateMorphing() const {
    TemplateMorpher morpher(fName, fSuffix);
    morpher.SetFitDimensionality(fTemplateMorphingFitDimensionality);
    morpher.SetMorphingSetting(&fMorphingSetting);
    morpher.ProcessTemplatesAndPlot(fRegions);
}

//__________________________________________________________________________________
//
void TRExFit::ProcessShapeFactorReparametrization() {
    if (!fShapeFactorReparametrisation) return;
    for (const auto& isf : fShapeFactors) {
        if (isf->fExpression.empty()) continue;
        for (const auto& iexpr : isf->fExpression) {
            const auto idep = Common::processString(iexpr.second);
            for (const auto& i : idep) {
                const std::string& param = i.first;
                const std::vector<double>& range = i.second;

                // add the NF and POIs
                auto itr = std::find(fNormFactorNames.begin(), fNormFactorNames.end(), param);
                if (itr != fNormFactorNames.end()) continue;
                LOG(INFO) << "Adding NormFactor: " << param << " used in the ShapeFactor reparametrization\n";
                std::shared_ptr<NormFactor> nf = std::make_shared<NormFactor>(param);
                TRExFitter::SYSTMAP[nf->fName] = nf->fName;
                TRExFitter::NPMAP[nf->fName] = nf->fName;
                fNormFactorNames.emplace_back(param);
                nf->fNuisanceParameter = param;
                nf->fTitle = param;
                nf->fRegions = isf->fRegions;
                nf->SetNominalMinMax(range.at(0), range.at(1), range.at(2));

                fNormFactors.emplace_back(nf);
            }
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::ProcessTemplateMorphingShapeFactors() {
    if (!fHasTemplateMorphing) return;
    const std::string filePath = fName + "/Templates/Parametrization.txt";
    ReparametrizationManager manager{};
    manager.ReadReparametrizationFile(filePath);
    for (const auto& isample : fSamples) {
        if (isample->fTemplateMorphing.empty()) continue;
        for (const auto& isf : isample->fShapeFactors) {
            const auto reparam = manager.GetParametrization(isf->fRegions.at(0), isample->fName);
            isf->fExpression = reparam;
        }
    }
}

//__________________________________________________________________________________
//
Region* TRExFit::GetRegion(const std::string& regName) const {
    auto itr = std::find_if(fRegions.begin(), fRegions.end(), [&regName](const auto& reg){return reg->fName == regName;});

    if (itr == fRegions.end()) {
        return nullptr;
    }

    return itr->get();
}

//__________________________________________________________________________________
//
std::vector<double> TRExFit::CalculateGlobalCorrelation() {
    if (fBootstrap!="" && fBootstrapIdx>=0){
        ReadFitResults(fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fInputName+fSuffix+".root");
    } else {
        ReadFitResults(fName+"/Fits/"+fInputName+fSuffix+".root");
    }

    std::vector<double> results;

    for (const auto& iunfolding : fUnfolding) {
        std::vector<std::string> bins;
        for (int ibin = 0; ibin < iunfolding->fNumberUnfoldingTruthBins; ++ibin) {
            if (iunfolding->fUnfoldNormXSec && (iunfolding->fUnfoldNormXSecBinN == ibin + 1)) continue;
            const std::string binName = iunfolding->fName + "_Bin_" + Common::IntToFixLenStr(ibin+1) + "_mu";
            bins.emplace_back(binName);
        }

        TMatrixDSym cov(bins.size());
        for (int i = 0; i < cov.GetNrows(); ++i) {
            for (int j = i; j < cov.GetNrows(); ++j) {
                const double cor = fFitResults->GetCorrelationMatrix()->GetCorrelation(bins.at(i), bins.at(j));
                double sigma_i = 0.5 * std::abs(fFitResults->GetNuisParErrUp(bins.at(i)) - fFitResults->GetNuisParErrDown(bins.at(i)));
                double sigma_j = 0.5 * std::abs(fFitResults->GetNuisParErrUp(bins.at(j)) - fFitResults->GetNuisParErrDown(bins.at(j)));
                cov(i,j) = sigma_i * cor * sigma_j;
                if (i != j) {
                    cov(j,i) = cov(i,j);
                }
            }
        }

        TMatrixDSym covInverted(cov);
        covInverted.Invert();

        //calculate the coefficient as <sqrt(1 - (V_ii V_ii^-1)^-1)>
        double coef(0.);
        for (std::size_t i = 0; i < bins.size(); ++i) {
            coef += std::sqrt(1 - 1./(cov(i,i)* covInverted(i,i)));
        }

        results.emplace_back(coef/bins.size());
    }

    return results;
}

//__________________________________________________________________________________
//
void TRExFit::ReadRegionBinning() {
    for (const auto& ireg : fRegions) {
        std::ifstream file;
        file.open(fBinningsPath + ireg->fName+".txt");

        if (!file.is_open() || !file.good()) {
            LOG(ERROR) << "Cannot read " << fBinningsPath + ireg->fName+".txt file\n";
            exit(EXIT_FAILURE);
        }

        int nbins;

        while (file >> nbins) {
            auto skippeditr = std::find(fSkippedRegions.begin(), fSkippedRegions.end(), ireg->fName);
            if (skippeditr != fSkippedRegions.end()) continue;
            ireg->SetNbins(nbins);
        }
        file.close();
    }
}

//__________________________________________________________________________________
//
std::vector<Region*> TRExFit::GetRawRegionVector() const {
    std::vector<Region*> result;
    for (const auto& reg : fRegions) {
        result.emplace_back(reg.get());
    }

    return result;
}

//__________________________________________________________________________________
//
void TRExFit::NonProfiledFit() {
    if(fPOIs.size() != 1) {
        LOG(ERROR) << "Non-profiled fit available only for single POI. Number of POI defined: " << fPOIs.size() << "\n";
        exit(EXIT_FAILURE);
    }

    if (fFitPOIAsimov.empty()) {
        auto itr = std::find_if(fNormFactors.begin(), fNormFactors.end(), [this](const auto& element){return element->fName == this->fPOIs.at(0);});
        if (itr == fNormFactors.end()) {
            LOG(ERROR) << "Cannot find NormFactor: " << fPOIs.at(0) << "\n";
            exit(EXIT_FAILURE);
        }

        fFitPOIAsimov.insert({fPOIs.at(0), (*itr)->GetNominal()});
    }

    LOG(INFO) << "\n";
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Scan of systematics for non-profile fit...\n";
    std::ofstream out;
    std::ofstream tex;
    std::ofstream out2;
    std::ofstream tex2;
    std::vector<std::string> regionsToFit;
    if(fBootstrap != "" && fBootstrapIdx >= 0) {
        out.open((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts.txt").c_str());
        tex.open((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts.tex").c_str());
        tex2.open((fName+"/Fits/"+fBootstrapSyst+fBootstrapSample+Form("_BSId%d/",fBootstrapIdx)+fName+fSuffix+"_nonProfiledSysts_grouped.tex").c_str());
    } else {
        out.open((fName+"/Fits/"+fInputName+fSuffix+"_nonProfiledSysts.txt").c_str());
        tex.open((fName+"/Fits/"+fInputName+fSuffix+"_nonProfiledSysts.tex").c_str());
        tex2.open((fName+"/Fits/"+fInputName+fSuffix+"_nonProfiledSysts_grouped.tex").c_str());
    }

    if (!out.good() || !out.is_open()) {
        LOG(ERROR) << "Cannot open the txt file with the results\n";
        exit(EXIT_FAILURE);
    }

    if (!tex.good() || !tex.is_open()) {
        LOG(ERROR) << "Cannot open the tex file with the results\n";
        exit(EXIT_FAILURE);
    }

    if (!tex2.good() || !tex2.is_open()) {
        LOG(ERROR) << "Cannot open the tex2 file with the results\n";
        exit(EXIT_FAILURE);
    }

    std::map<std::string, int> regionDataType;
    for(const auto& reg : fRegions) {
        regionDataType.insert(std::make_pair(reg->fName, Region::DataType::ASIMOVDATA));
    }
    for(const auto& ireg : fRegions) {
        if (fFitRegion == CRONLY && ireg->fRegionType == Region::RegionType::CONTROL) {
            regionsToFit.emplace_back(ireg->fName);
        } else if (fFitRegion == CRSR && (ireg->fRegionType == Region::RegionType::CONTROL || ireg->fRegionType == Region::RegionType::SIGNAL)) {
            regionsToFit.emplace_back(ireg->fName);
        }
    }
    auto ws = PerformWorkspaceCombinationxRooFit(regionsToFit);
    if (!ws) {
        LOG(ERROR) << "Cannot retrieve the workspace, exiting!\n";
        exit(EXIT_FAILURE);
    }
    // nominal fit on Asimov
    RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        LOG(ERROR) << "Cannot retrieve the model config, exiting!\n";
        exit(EXIT_FAILURE);
    }
    ws->saveSnapshot("InitialStateModelGlob", *mc->GetGlobalObservables());
    if (mc->GetNuisanceParameters()) {
        ws->saveSnapshot("InitialStateModelNuis", *mc->GetNuisanceParameters());
    }
    LOG(INFO) << "Fitting nominal Asimov...\n";
    auto data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, fFitNPValues, fFitPOIAsimov ));
    std::map<std::string, double> npValues = PerformFit( ws.get(), data.get(), fFitType, false);
    auto itr = npValues.find(fPOIs.at(0));
    if (itr == npValues.end()) {
        LOG(ERROR) << "Cannot find parameter: " << fPOIs.at(0) << ", exiting!\n";
        exit(EXIT_FAILURE);
    }
    const double nominalPOIval = itr->second;
    const RooRealVar* poiVar = static_cast<RooRealVar*>(& ws->allVars()[fPOIs.at(0).c_str()]);
    const double statUp = poiVar->getErrorHi();
    const double statDo = poiVar->getErrorLo();

    // temporary switch off minos (not needed)
    std::vector<std::string> varMinosTmp = fVarNameMinos;
    fVarNameMinos.clear();

    // MC stat
    double MCstatUp = 0.;
    double MCstatDo = 0.;
    if(fUseStatErr) {
        fGammasInStatOnly = true;
        LOG(INFO) << "MC stat...\n";
        for(const auto& reg : fRegions){
            std::unique_ptr<TH1> hTot = nullptr;
            for(auto& sh : reg->fSampleHists){
                if(sh->fSample->fSeparateGammas) continue;
                if(sh->fSample->fType==Sample::SampleType::GHOST) continue;
                if(sh->fSample->fType==Sample::SampleType::EFT) continue;
                if(sh->fSample->fType==Sample::SampleType::DATA) continue;
                if(!sh->fSample->fUseMCStat) continue; // need to fix something for separate gammas
                bool skip = false;
                for(const auto& morphPar : fMorphParams) {
                    // find nominal value of this morph parameter
                    double nfVal = 0.;
                    for(const auto& nf : fNormFactors) {
                        if(nf->fName == morphPar){
                            nfVal = nf->GetNominal();
                            break;
                        }
                    }
                    if(sh->fSample->fIsMorph[morphPar] && sh->fSample->fMorphValue[morphPar] != nfVal) skip = true;
                }
                if(skip) continue;
                LOG(INFO) << "  Including sample " << sh->fSample->fName << "\n";
                if(!hTot && sh->fHist) hTot.reset(static_cast<TH1*>(sh->fHist->Clone("h_tot")));
                else if(sh->fHist) hTot->Add(sh->fHist.get());
            }
            if(!hTot) continue;
            for(int ibin = 1; ibin <= hTot->GetNbinsX(); ++ibin) {
                const double statErr =  hTot->GetBinError(ibin)/hTot->GetBinContent(ibin);
                const std::string gammaName = "gamma_stat_"+reg->fName+"_bin_"+std::to_string(ibin-1);
                // up
                std::map<std::string, double> npVal = fFitNPValues;
                {
                    auto itrGamma = npVal.find(gammaName);
                    if (itrGamma == npVal.end()) {
                        npVal.insert(std::make_pair(gammaName, 1+statErr));
                    } else {
                        itrGamma->second = 1+statErr;
                    }
                }
                LOG(DEBUG) << "Setting " << gammaName << " to " << (1+statErr) << "\n";
                data = std::unique_ptr<RooDataSet>(DumpData(ws.get(), regionDataType, npVal, fFitPOIAsimov));
                npVal.at(gammaName) = 1;
                ws->loadSnapshot("InitialStateModelGlob");
                if (mc->GetNuisanceParameters()) {
                    ws->loadSnapshot("InitialStateModelNuis");
                }
                npValues = PerformFit( ws.get(), data.get(), fFitType, false);
                MCstatUp = std::hypot(MCstatUp, (npValues.at(fPOIs.at(0))-nominalPOIval));
                // down
                npVal = fFitNPValues;
                {
                    auto itrGamma = npVal.find(gammaName);
                    if (itrGamma == npVal.end()) {
                        npVal.insert(std::make_pair(gammaName, 1+statErr));
                    } else {
                        itrGamma->second = 1+statErr;
                    }
                }
                LOG(DEBUG) << "Setting " << gammaName << " to " << (1+statErr) << "\n";
                npVal.at(gammaName) = 1-statErr;
                LOG(DEBUG) << "Setting " << gammaName << " to " << (1-statErr) << "\n";
                data = std::unique_ptr<RooDataSet>(DumpData(ws.get(), regionDataType, npVal, fFitPOIAsimov));
                npVal.at(gammaName) = 1;
                ws->loadSnapshot("InitialStateModelGlob");
                if (mc->GetNuisanceParameters()) {
                    ws->loadSnapshot("InitialStateModelNuis");
                }
                npValues = PerformFit(ws.get(), data.get(), fFitType, false);
                MCstatDo = -std::hypot(MCstatDo,(npValues.at(fPOIs.at(0))-nominalPOIval));
            }
        }
        fGammasInStatOnly = false;
    }
    // MC-stat for specific samples (those using separate gammas)
    std::map<std::string,double> MCstatUpSample;
    std::map<std::string,double> MCstatDoSample;
    std::map<std::string,std::string> smpTexTitle;
    if(fUseStatErr) {
        std::map <std::string, double> npVal;
        for(const auto& reg : fRegions) {
            for(auto& sh : reg->fSampleHists) {
                std::unique_ptr<TH1> hTot = nullptr;
                if(sh->fSample->fType==Sample::SampleType::GHOST) continue;
                if(sh->fSample->fType==Sample::SampleType::EFT) continue;
                if(sh->fSample->fType==Sample::SampleType::DATA) continue;
                if(!sh->fSample->fSeparateGammas) continue;
                bool skip = false;
                for(const auto& morphPar : fMorphParams){
                    // find nominal value of this morph parameter
                    double nfVal = 0.;
                    for(const auto& nf : fNormFactors){
                        if(nf->fName==morphPar){
                            nfVal = nf->GetNominal();
                            break;
                        }
                    }
                    if(sh->fSample->fIsMorph[morphPar] && sh->fSample->fMorphValue[morphPar]!=nfVal) skip = true;
                }
                if (skip) continue;
                LOG(INFO) << "MC stat for sample " << sh->fSample->fName << "...\n";
                hTot.reset(static_cast<TH1*>(sh->fHist->Clone("h_tot")));
                if (!hTot) continue;
                for (int ibin = 1; ibin <= hTot->GetNbinsX(); ++ibin) {
                    const double statErr = hTot->GetBinError(ibin)/hTot->GetBinContent(ibin);
                    const std::string gammaName = "gamma_shape_stat_"+sh->fSample->fName+"_"+reg->fName+"_bin_"+std::to_string(ibin-1);
                    npVal = fFitNPValues;
                    auto itrGamma = npVal.find(gammaName);
                    if (itrGamma == npVal.end()) {
                        npVal.insert(std::make_pair(gammaName, 1+statErr));
                    } else {
                        itrGamma->second = 1+statErr;
                    }
                    LOG(DEBUG) << "Setting " << gammaName << " to " << (1+statErr) << "\n";
                    data = std::unique_ptr<RooDataSet>(DumpData(ws.get(), regionDataType, npVal, fFitPOIAsimov));
                    npVal.at(gammaName) = 1;
                    ws->loadSnapshot("InitialStateModelGlob");
                    if (mc->GetNuisanceParameters()) {
                        ws->loadSnapshot("InitialStateModelNuis");
                    }
                    npValues = PerformFit(ws.get(), data.get(), fFitType, false);
                    MCstatUpSample.insert({sh->fSample->fName, std::hypot(MCstatUpSample[sh->fSample->fName], (npValues.at(fPOIs.at(0))-nominalPOIval))});
                    npVal = fFitNPValues;
                    npVal.at(gammaName) = 1-statErr;
                    LOG(DEBUG) << "Setting " << gammaName << " to " << (1-statErr) << "\n";
                    data = std::unique_ptr<RooDataSet>(DumpData(ws.get(), regionDataType, npVal, fFitPOIAsimov));
                    npVal.at(gammaName) = 1;
                    ws->loadSnapshot("InitialStateModelGlob");
                    if (mc->GetNuisanceParameters()) {
                        ws->loadSnapshot("InitialStateModelNuis");
                    }
                    npValues = PerformFit(ws.get(), data.get(), fFitType, false);
                    MCstatDoSample.insert({sh->fSample->fName, -std::hypot(MCstatDoSample[sh->fSample->fName],(npValues.at(fPOIs.at(0))-nominalPOIval))});
                    smpTexTitle.insert({sh->fSample->fName, sh->fSample->fTexTitle});
                }
            }
        }
    }
    //
    std::map<std::string, double> newPOIvalUp;
    std::map<std::string, double> newPOIvalDo;
    std::map<std::string, double> systGroups;
    std::vector<std::string> systGroupNames;
    std::vector<std::string> npList;
    for(const auto& syst : fSystematics){
        if(Common::FindInStringVector(npList,syst->fNuisanceParameter)<0){
            npList.emplace_back(syst->fNuisanceParameter);
            for(int ud = 0; ud < 2; ++ud){
                //Be sure to take the initial values of the NP
                ws->loadSnapshot("InitialStateModelGlob");
                if (mc->GetNuisanceParameters()) {
                    ws->loadSnapshot("InitialStateModelNuis");
                }
                // - create Asimov with that NP fixed to +/-1sigma
                std::map<std::string, double> npVal;
                npVal = fFitNPValues;
                if(ud==0) npVal["alpha_"+syst->fNuisanceParameter] =  1;
                if(ud==1) npVal["alpha_"+syst->fNuisanceParameter] = -1;
                LOG(INFO) << "Systematic " << syst->fNuisanceParameter << "...\n";
                data = std::unique_ptr<RooDataSet>(DumpData( ws.get(), regionDataType, npVal, fFitPOIAsimov ));
                // - again a stat-only fit to that Asimov
                fFitFixedNPs[syst->fNuisanceParameter] = 0;
                fFitFixedNPs["alpha_"+syst->fNuisanceParameter] = 0;
                npValues = PerformFit( ws.get(), data.get(), fFitType, false);
                double newPOIval = npValues.at(fPOIs.at(0));
                if(ud==0) newPOIvalUp[syst->fNuisanceParameter] = newPOIval;
                if(ud==1) newPOIvalDo[syst->fNuisanceParameter] = newPOIval;
            }
            std::string category = syst->fCategory;
            if(syst->fSubCategory!="") {
                category = syst->fSubCategory;
            }
            if(Common::FindInStringVector(systGroupNames,category)<0) {
                systGroupNames.push_back(category);
            }
            if(std::abs(newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval) > fNonProfileFitSystThreshold
            || std::abs(newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval) > fNonProfileFitSystThreshold)
                systGroups[category] = std::hypot(systGroups[category], ((std::abs(newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval)+std::abs(newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval))/2));
        }
    }
    // print:
    LOG(INFO) << "Results of non-profile fit:\n";
    LOG(INFO) << "-----------------------------------\n";
    LOG(INFO) << "StatisticalError\t" << statUp << "\t" << statDo << "\n";
    LOG(INFO) << "-----------------------------------\n";
    out       << "StatisticalError\t" << statUp << "\t" << statDo << std::endl;
    tex       << "\\begin{tabular}{lr}" << std::endl;
    tex       << "\\hline" << std::endl;
    tex       << "\\hline" << std::endl;
    tex       << "  Source & Shift up / down";
    if(fPOIunit.at(fPOIs.at(0))!="") tex       << " [" << fPOIunit.at(fPOIs.at(0)) << "]";
    tex       << "\\\\" << std::endl;
    tex       << "\\hline" << std::endl;
    tex       << "  Statistical & $+" << Form("%.2f",statUp) << "$ / $" << Form("%.2f",statDo) << "$ \\\\" << std::endl;
    tex       << "\\hline" << std::endl;
    //
    double totUp = 0.;
    double totDo = 0.;
    npList.clear();
    // MC stat
    LOG(INFO) << "Stat.MC\t" << MCstatUp << "\t" << MCstatDo << "\n";
    out       << "Stat.MC\t" << MCstatUp << "\t" << MCstatDo << std::endl;
    tex       << "  MC-stat & $+" << Form("%.2f",MCstatUp) << "$ / $" << Form("%.2f",MCstatDo) << "$ \\\\" << std::endl;
    tex       << "\\hline" << std::endl;
    totUp = std::hypot(totUp,MCstatUp);
    totDo = std::hypot(totDo,MCstatDo);
    // MC stat for separate gamma samples
    for(const auto& sepGammaPair : MCstatUpSample){
        std::string smpName = sepGammaPair.first;
        LOG(INFO) << "Stat." << smpName << "\t" << MCstatUpSample[smpName] << "\t" << MCstatDoSample[smpName] << "\n";
        out       << "Stat." << smpName << "\t" << MCstatUpSample[smpName] << "\t" << MCstatDoSample[smpName] << std::endl;
        tex       << "  Stat (" << smpTexTitle[smpName] << ") & $+" << Form("%.2f",MCstatUpSample[smpName]) << "$ / $" << Form("%.2f",MCstatDoSample[smpName]) << "$ \\\\" << std::endl;
        tex       << "\\hline" << std::endl;
        totUp = std::hypot(totUp,MCstatUpSample[smpName]);
        totDo = std::hypot(totDo,MCstatDoSample[smpName]);
    }
    LOG(INFO) << "-----------------------------------\n";
    // systematics
    for(const auto& syst : fSystematics) {
        if(syst->fName.find("stat_")!=std::string::npos && syst->fType==Systematic::SHAPE) continue;
        if(Common::FindInStringVector(npList,syst->fNuisanceParameter)<0){
            npList.push_back(syst->fNuisanceParameter);
            LOG(INFO) << syst->fNuisanceParameter << "\n";
            out       << syst->fNuisanceParameter;
            if(TRExFitter::SYSTTEX[syst->fNuisanceParameter]!="") tex << "  " << TRExFitter::SYSTTEX[syst->fNuisanceParameter];
            else                                                  tex << "  " << TRExFitter::SYSTMAP[syst->fNuisanceParameter];
            // - up and down
            for(int ud=0;ud<2;ud++){
                double valUp = newPOIvalUp[syst->fNuisanceParameter]-nominalPOIval;
                double valDo = newPOIvalDo[syst->fNuisanceParameter]-nominalPOIval;
                if (ud == 0) LOG(INFO) << "\t" << valUp << "\n";
                if (ud == 1) LOG(INFO) << "\t" << valDo << "\n";
                if(ud==0) out       << "\t" << valUp;
                if(ud==1) out       << "\t" << valDo;
                if(ud==0) tex       << " & " << Form("$%s%.2f$",(valUp>=0 ? "+" : "-"),std::abs(valUp));
                if(ud==1) tex       << " / " << Form("$%s%.2f$",(valDo>=0 ? "+" : "-"),std::abs(valDo));
                if(std::abs(valUp)>fNonProfileFitSystThreshold || std::abs(valDo)>fNonProfileFitSystThreshold) {
                    if(ud==0 && valUp>0) totUp = std::hypot(totUp,valUp);
                    if(ud==0 && valUp<0) totDo = std::hypot(totDo,valUp);
                    if(ud==1 && valDo>0) totUp = std::hypot(totUp,valDo);
                    if(ud==1 && valDo<0) totDo = std::hypot(totDo,valDo);
                }
            }
            LOG(INFO) << "\n";
            out       << std::endl;
            tex       << "\\\\" << std::endl;
        }
    }
    totUp = std::hypot(totUp,fNonProfileFitSystThreshold);
    totDo = std::hypot(totDo,fNonProfileFitSystThreshold);
    LOG(INFO) << "-----------------------------------\n";
    LOG(INFO) << "TotalSystematic\t" << totUp << "\t" << totDo << "\n";
    out       << "TotalSystematic\t" << totUp << "\t-" << totDo << std::endl;
    LOG(INFO) << "-----------------------------------\n";
    LOG(INFO) << "TotalStat+Syst\t" << std::hypot(totUp, statUp) << "\t" << std::hypot(totDo, statDo) << "\n";
    out       << "TotalStat+Syst\t" << std::hypot(totUp,statUp) << "\t-" << std::hypot(totDo,statDo) << std::endl;
    LOG(INFO) << "-----------------------------------\n";
    out.close();
    tex << "\\hline" << std::endl;
    tex << "  Total systematics & $+" << Form("%.2f",totUp) << "$ / $-" << Form("%.2f",totDo) << "$ \\\\" << std::endl;
    tex << "\\hline" << std::endl;
    tex << "  Total stat+syst & $+"   << Form("%.2f",std::hypot(totUp,statUp))     << "$ / $-" << Form("%.2f",std::hypot(totDo,statDo))     << "$ \\\\" << std::endl;
    tex << "\\hline" << std::endl;
    tex << "\\hline" << std::endl;
    tex << "\\end{tabular}" << std::endl;
    tex.close();
    fVarNameMinos = varMinosTmp; // retore Minos settings
    //
    // Systematics merged according to syst groups
    LOG(INFO) << "-----------------------------------\n";
    LOG(INFO) << "- Systematic impact per category  -\n";
    LOG(INFO) << "-----------------------------------\n";
    tex2      << "\\begin{tabular}{lr}" << std::endl;
    tex2      << "\\hline" << std::endl;
    tex2      << "\\hline" << std::endl;
    LOG(INFO) << "Data statistics\t" << ((std::abs(statUp)+std::abs(statDo))/2) << "\n";
    out2      << "Data statistics" << "\t" << (std::abs(statUp)+std::abs(statDo))/2. << std::endl;
    tex2      << "Data statistics" << " & " << Form("$%.2f$",(std::abs(statUp)+std::abs(statDo))/2.) << " \\\\" << std::endl;
    LOG(INFO) << "-----------------------------------\n";
    tex2      << "\\hline" << std::endl;
    LOG(INFO) << "MC background stat.\t" << ((std::abs(MCstatUp)+std::abs(MCstatDo))/2) << "\n";
    out2      << "MC background stat." << "\t" << ((std::abs(MCstatUp)+std::abs(MCstatDo))/2.) << std::endl;
    tex2      << "MC background stat." << " & " << Form("$%.2f$",(std::abs(MCstatUp)+std::abs(MCstatDo))/2.) << " \\\\" << std::endl;
    for(auto sepGammaPair : MCstatUpSample){
        const std::string smpName = sepGammaPair.first;
        LOG(INFO) << smpTexTitle[smpName] << " stat.\t" << (std::abs(MCstatUpSample[smpName])+std::abs(MCstatDoSample[smpName])/2) << "\n";
        out2      << smpTexTitle[smpName] << " stat." << "\t" << (std::abs(MCstatUpSample[smpName])+std::abs(MCstatDoSample[smpName]))/2. << std::endl;
        tex2      << smpTexTitle[smpName] << " stat." << " & " << Form("$%.2f$",(std::abs(MCstatUpSample[smpName])+std::abs(MCstatDoSample[smpName]))/2.) << " \\\\" << std::endl;
    }
    LOG(INFO) << "-----------------------------------\n";
    tex2      << "\\hline" << std::endl;
    for(const auto& systGroupName : systGroupNames){
        if(systGroupName=="") continue;
        LOG(INFO) << systGroupName << "\t" << systGroups[systGroupName] << "\n";
        out2      << systGroupName << "\t" << systGroups[systGroupName] << std::endl;
        tex2      << systGroupName << " & " << Form("$%.2f$",systGroups[systGroupName]) << " \\\\" << std::endl;
    }
    LOG(INFO) << "-----------------------------------\n";
    LOG(INFO) << "TotalSystematic\t" << ((std::abs(totUp)+std::abs(totDo))/2) << "\n";
    LOG(INFO) << "TotalStat+Syst\t" << ((std::hypot(totUp,statUp)+std::hypot(totDo,statDo))/2.) << "\n";
    LOG(INFO) << "-----------------------------------\n";
    tex2 << "\\hline" << std::endl;
    tex2 << "Total systematic uncertainty & " << Form("$%.2f$",(std::abs(totUp)+std::abs(totDo))/2.) << " \\\\" << std::endl;
    tex2 << "\\hline" << std::endl;
    tex2 << "Total & "   << Form("$%.2f$",(std::hypot(totUp,statUp)+std::hypot(totDo,statDo))/2.) << " \\\\" << std::endl;
    tex2 << "\\hline" << std::endl;
    tex2 << "\\hline" << std::endl;
    tex2 << "\\end{tabular}" << std::endl;
    tex2.close();
}

//__________________________________________________________________________________
//
void TRExFit::PlotMultipleUnfoldingCovariance(const std::vector<const Unfolding*>& unfoldings,
                                              const std::vector<std::string>& truthFileNames,
                                              const std::string& frPath,
                                              const std::string& wsFileName,
                                              const std::string& outputPath,
                                              const std::string& wsName,
                                              const std::vector<std::shared_ptr<NormFactor> >& nfs) {

    if (unfoldings.size() < 2) {
        LOG(ERROR) << "Number of unfoldings is < 2\n";
        return;
    }

    if (unfoldings.size() != truthFileNames.size()) {
        LOG(ERROR) << "Sizes of unfoldings and truth file names do not match\n";
        return;
    }

    for (const auto* iunf : unfoldings) {
        if (iunf->HasReparametrisation(nfs)) {
            LOG(INFO) << "Unfolding: " << iunf->fName << " has reparametrised NormFactors, will not produce the combined covariance matrix\n";
            return;
        }
    }

    gSystem->mkdir((outputPath + "/UnfoldingPlots/").c_str());

    this->ReadFitResults(frPath+".root");

    std::vector<UnfoldingResult> unfoldedResults(unfoldings.size());

    int size(0);
    std::vector<std::string> names;
    std::vector<std::string> plotNames;
    std::vector<UnfoldingResult::FitValue> results;
    std::vector<std::vector<std::string> > binNames;
    std::vector<std::string> unfoldingNames;
    std::size_t index(0);
    for (const auto* iunf : unfoldings) {
        LOG(INFO) << "Adding unfolding " << iunf->fName << "\n";
        std::unique_ptr<TFile> file(TFile::Open((truthFileNames.at(index)+"/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
        if (!file) {
            LOG(ERROR) << "Cannot open the input file for the truth histogram\n";
            return;
        }

        std::unique_ptr<TH1> truth(file->Get<TH1>((iunf->fName+"_truth_distribution").c_str()));
        if (!truth) {
            LOG(ERROR) << "Cannot read the input histogram: " << iunf->fName+"_truth_distribution\n";
            return;
        }
        truth->SetDirectory(nullptr);

        unfoldedResults.at(index).SetTruthDistribution(truth.get());
        unfoldedResults.at(index).SetNormXSec(iunf->fUnfoldNormXSec);
        unfoldedResults.at(index).SetDivideByBinWidth(iunf->fUnfoldingDivideByBinWidth);
        if (iunf->fUnfoldingDivideByLumi && !iunf->fUnfoldNormXSec) {
            unfoldedResults.at(index).SetDivideByLumi(iunf->fUnfoldingDivideByLumi);
        }
        Common::GetUnfoldingResult(unfoldedResults.at(index), iunf, truth.get(), fFitResults.get(), wsFileName, wsName, frPath+".root", nfs, {}, false);

        file->Close();
        size += iunf->fNumberUnfoldingTruthBins;
        std::vector<std::string> tmpVec;
        for (int i = 0; i < iunf->fNumberUnfoldingTruthBins; ++i) {
            const std::string name = iunf->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
            names.emplace_back(std::move(name));
            plotNames.emplace_back(iunf->fName + " Bin " + Common::IntToFixLenStr(i+1));
            tmpVec.emplace_back(iunf->fName + " Bin " + Common::IntToFixLenStr(i+1));
        }
        binNames.emplace_back(tmpVec);
        const auto& tmp = unfoldedResults.at(index).GetFittedResults();
        results.insert(results.end(), tmp.begin(), tmp.end());
        unfoldingNames.emplace_back(iunf->fName);
        index++;
    }

    std::vector<std::vector<double> > cov(size, std::vector<double>(size));
    std::vector<std::vector<double> > cor(size, std::vector<double>(size));
    TH2D histo("","", size, 0.5, size+0.5, size, 0.5, size+0.5);
    TH2D histoCorr("","", size, 0.5, size+0.5, size, 0.5, size+0.5);
    for (int i = 0; i < size; ++i) {
        const std::string& name_i = names.at(i);
        for (int j = 0; j < size; ++j) {
            const std::string& name_j = names.at(j);
            const double corr = (i == j) ? 1. : fFitResults->GetCorrelationMatrix()->GetCorrelation(name_i, name_j);
            const double up_i   = results.at(i).up;
            const double down_i = results.at(i).down;
            const double up_j   = results.at(j).up;
            const double down_j = results.at(j).down;
            const double sigma_i = 0.5 * (std::abs(up_i) + std::abs(down_i));
            const double sigma_j = 0.5 * (std::abs(up_j) + std::abs(down_j));

            // this is covariance i,j
            const double covariance = sigma_i * corr * sigma_j;
            histo.SetBinContent(i+1,size-j, covariance);
            histoCorr.SetBinContent(i+1,size-j, corr);
            cov.at(i).at(j) = covariance;
            cor.at(i).at(j) = corr;
        }
        histo.GetXaxis()->SetBinLabel(i+1, (plotNames.at(i)).c_str());
        histo.GetYaxis()->SetBinLabel(size-i, (plotNames.at(i)).c_str());
        histoCorr.GetXaxis()->SetBinLabel(i+1, (plotNames.at(i)).c_str());
        histoCorr.GetYaxis()->SetBinLabel(size-i, (plotNames.at(i)).c_str());
    }

    YamlConverter converter{};
    converter.SetLumi(Common::ReplaceString(fLumiLabel, " fb^{-1}", ""));
    converter.SetCME(Common::ReplaceString(fCmeLabel, " TeV", "000"));
    converter.WriteCorrelation(plotNames, cov, outputPath + "/UnfoldingPlots", true, "CovarianceMatrix_combined");
    converter.WriteCorrelation(plotNames, cor, outputPath + "/UnfoldingPlots", false, "CorrelationMatrix_combined");
    if (fHEPDataFormat) {
        converter.WriteCorrelationHEPData(plotNames, cov, outputPath, true, "CovarianceMatrix_combined");
        converter.WriteCorrelationHEPData(plotNames, cor, outputPath, false, "CorrelationMatrix_combined");
    }

    TCanvas c1("","",0.,0.,1200,1200);
    gStyle->SetPalette(87);
    #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
    histo.SetMarkerSize(0.75*1000);
    #else
    histo.SetMarkerSize(1.0);
    #endif
    gStyle->SetPaintTextFormat("1.1e");
    c1.SetBottomMargin(0.2);
    c1.SetLeftMargin(0.2);

    histo.GetXaxis()->LabelsOption("v");
    histo.GetXaxis()->SetLabelSize(histo.GetXaxis()->GetLabelSize()*0.75);
    histo.GetYaxis()->SetLabelSize(histo.GetYaxis()->GetLabelSize()*0.75);
    c1.SetTickx(0);
    c1.SetTicky(0);
    histo.GetYaxis()->SetTickLength(0);
    histo.GetXaxis()->SetTickLength(0);
    histo.Draw("col TEXT");

    Common::SaveCanvasAs(c1, outputPath + "/UnfoldingPlots/CovarianceMatrix_combined_unfolding");

    TCanvas c2("","",0.,0.,1200,1200);
    #if ROOT_VERSION_CODE < ROOT_VERSION(6,35,0)
    histo.SetMarkerSize(0.75*1000);
    #else
    histo.SetMarkerSize(1.0);
    #endif
    c2.SetBottomMargin(0.2);
    c2.SetLeftMargin(0.2);

    histoCorr.GetXaxis()->LabelsOption("v");
    histoCorr.GetXaxis()->SetLabelSize(histoCorr.GetXaxis()->GetLabelSize()*0.75);
    histoCorr.GetYaxis()->SetLabelSize(histoCorr.GetYaxis()->GetLabelSize()*0.75);
    c2.SetTickx(0);
    c2.SetTicky(0);
    histoCorr.GetYaxis()->SetTickLength(0);
    histoCorr.GetXaxis()->SetTickLength(0);
    histoCorr.Draw("col TEXT");

    Common::SaveCanvasAs(c2, outputPath + "/UnfoldingPlots/CorrelationMatrix_combined_unfolding");

    this->SplitUnfoldingCorrelationMatrix(cor, binNames, unfoldingNames, outputPath);
}

//__________________________________________________________________________________
//
void TRExFit::SplitUnfoldingCorrelationMatrix(const std::vector<std::vector<double> >& cor,
                                              const std::vector<std::vector<std::string> >& binNames,
                                              const std::vector<std::string>& unfoldingNames,
                                              const std::string& folder) const {

    // some checks first
    if (cor.empty()) {
        LOG(ERROR) << "Correlation matrix is empty\n";
        return;
    }

    if (binNames.empty()) {
        LOG(ERROR) << "Names of the bins are empty\n";
        return;
    }

    if (binNames.size() != unfoldingNames.size()) {
        LOG(ERROR) << "Names of the unfolded measurements and the bin Names do not match\n";
        return;
    }

    auto SumBefore = [](const std::vector<std::vector<std::string> >& vec, const std::size_t index) {
        std::size_t sum(0);
        for (std::size_t i = 0; i < vec.size(); ++i) {
            if (i < index) {
                sum += vec.at(i).size();
            }
        }

        return sum;
    };

    for (std::size_t imeas = 0; imeas < binNames.size(); ++imeas) {
        for (std::size_t jmeas = imeas+1; jmeas < binNames.size(); ++jmeas) {
            const std::size_t offsetX = SumBefore(binNames, jmeas);
            const std::size_t offsetY = SumBefore(binNames, imeas);
            std::vector<std::vector<double> > correlations;
            for (std::size_t ibin = 0; ibin < binNames.at(jmeas).size(); ++ibin) {
                std::vector<double> tmp;
                for (std::size_t jbin = 0; jbin < binNames.at(imeas).size(); ++jbin) {
                    const double corr = cor.at(offsetX + ibin).at(offsetY + jbin);
                    tmp.emplace_back(corr);
                }
                correlations.emplace_back(std::move(tmp));
            }

            // store in the output
            std::ofstream out(folder+"/Fits/Correlations_unfolding_"+unfoldingNames.at(imeas) + "_" + unfoldingNames.at(jmeas) + ".txt");
            if (!out.is_open() || !out.good()) {
                LOG(ERROR) << "Cannot open file at: " << folder+"/Correlations_unfolding_"+unfoldingNames.at(imeas) + "_" + unfoldingNames.at(jmeas) + ".txt\n";
                return;
            }

            out << "Correlations\n";
            out << "bin | ";
            for (const auto& i : binNames.at(imeas)) {
                out << i << " | ";
            }
            out << "\n";
            for (std::size_t i = 0; i < correlations.size(); ++i) {
                out << binNames.at(jmeas).at(i) << " | ";
                for (const auto& j : correlations.at(i)) {
                    out << j << " ";
                }
                out << "\n";
            }

            out.close();
        }
    }
}

//__________________________________________________________________________________
//
bool TRExFit::HasValidationRegions() const {
    auto itr = std::find_if(fRegions.begin(), fRegions.end(), [](const auto& ireg){return ireg->fRegionType == Region::RegionType::VALIDATION;});

    return itr != fRegions.end();
}

//__________________________________________________________________________________
//
bool TRExFit::HasDropBinRegions() const {
    auto itr = std::find_if(fRegions.begin(), fRegions.end(), [](const auto& ireg){return (ireg->fDropBins.size() > 0 || ireg->GetAutomaticDropBins());});
    return itr != fRegions.end();
}

//__________________________________________________________________________________
//
void TRExFit::AddWSMetadata(RooWorkspace* ws, const bool allRegions) const {

    LOG(INFO) << "Attaching metadata (titles, colours) to the WS\n";
    ws->import(*gROOT->GetListOfColors(),true);
    xRooNode node(*ws);

    this->AddParameterTitlesToMetadata(node);

    auto pdf = node["simPdf"];
    if (!pdf) {
        LOG(ERROR) << "Cannot read \"pdf\"\n";
        exit(EXIT_FAILURE);
    }
    node["asimovData"]->SetTitle("Asimov Data");
    auto data = node["obsData"];
    if (data) {
        data->SetTitle("Data");
    }

    for(const auto& ireg : fRegions) {
        if(!allRegions && ireg->fRegionType==Region::VALIDATION) continue;
        auto regNode = pdf->find(ireg->fName);
        if (!regNode) {
            LOG(ERROR) << "Cannot find region node: " << ireg->fName << " from WS\n";
            exit(EXIT_FAILURE);
        }
        regNode->SetTitle(ireg->fLabel.c_str());
        auto samples = regNode->find("samples");
        if (!samples) {
            LOG(ERROR) << "Cannot find samples in region node: " << ireg->fName << " from WS\n";
            exit(EXIT_FAILURE);
        }
        for (const auto& isample : ireg->fSampleHists) {
            if(isample->fSample->fType==Sample::SampleType::DATA) continue;
            if(isample->fSample->fType==Sample::SampleType::GHOST) continue;
            if(isample->fSample->fType==Sample::SampleType::EFT) continue;

            auto sampleNode = samples->find((isample->fSample->fName+"_"+ireg->fName+"_shapes").c_str());
            if (!sampleNode) {
                LOG(ERROR) << "Cannot find sample: " << isample->fSample->fName << ", in region: " << ireg->fName << "\n";
                exit(EXIT_FAILURE);
            }
            const std::string& title = isample->fSample->fTitle;
            sampleNode->SetTitle(title.c_str());
            sampleNode->styles().get<TStyle>()->SetFillColor(isample->fSample->fFillColor);
            sampleNode->styles().get<TStyle>()->SetLineColor(isample->fSample->fLineColor);
        }
    }
}

void TRExFit::AddParameterTitlesToMetadata(const xRooNode& node) const {

    for (auto& iparam : node.vars()) {
        const std::string parName = iparam->GetName();

        // test NPs
        auto itrSyst = std::find_if(fSystematics.begin(), fSystematics.end(),
            [parName](const auto& element){return parName.find(element->fNuisanceParameter) != std::string::npos;});

        if (itrSyst != fSystematics.end()) {
            iparam->SetTitle((*itrSyst)->fTitle.c_str());
        } else {
            // test NFs
            auto itrNF = std::find_if(fNormFactors.begin(), fNormFactors.end(),
                [parName](const auto& element){return parName == element->fNuisanceParameter;});

            if (itrNF != fNormFactors.end()) {
                iparam->SetTitle((*itrNF)->fTitle.c_str());
            }
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::PlotRankingFromCovMatrix(const FittingTool& fitTool) const {

    RankingManager manager{};
    manager.SetPlotLabel(fPlotLabel);
    manager.SetLumiLabel(fLumiLabel);
    manager.SetCmeLabel(fCmeLabel);
    manager.SetUseHEPDataFormat(fHEPDataFormat);
    manager.SetName(fName);
    manager.SetMaxNPPlot(fRankingMaxNP);
    manager.SetRankingCanvasSize(fNPRankingCanvasSize);
    manager.SetShiftGlobalObservables(true);

    std::vector<std::string> exclude = fPOIs;
    for (const auto& inf : fNormFactorNames) {
        exclude.emplace_back(inf);
    }
    for (const auto& isf : fShapeFactorNames) {
        exclude.emplace_back(isf);
    }

    for (std::size_t ipoi = 0; ipoi < fPOIs.size(); ++ipoi) {
        if (fFitType == FitType::BONLY) continue;
        manager.SetRankingPOIName(fRankingPOIName.at(ipoi));
        manager.SetUpperAxisNdivisions(fRankingUpperAxisNdivision.at(ipoi));
        manager.SetRankingPOIAxisScale(fRankingPOIAxisScale.at(ipoi));
        manager.SetIsCovarianceBreakdown(true);
        if (fRankingPlot=="MERGE"  || fRankingPlot=="ALL") {
            const auto container = fitTool.GetRankingContainer(fPOIs.at(ipoi), exclude, true, true);
            manager.SetRankingContainer(container);
            manager.SetSuffix(fSuffix+"_"+fPOIs.at(ipoi)+"_Breakdown");
            manager.PlotRanking(fRegions, fNormFactorNames, fShapeFactorNames, true, true);
        }
        if (fRankingPlot=="SYSTS"  || fRankingPlot=="ALL") {
            const auto container = fitTool.GetRankingContainer(fPOIs.at(ipoi), exclude, true, false);
            manager.SetRankingContainer(container);
            manager.SetSuffix(fSuffix+"_"+fPOIs.at(ipoi)+"_Breakdown_syst");
            manager.PlotRanking(fRegions, fNormFactorNames, fShapeFactorNames, true, false);
        }
        if (fRankingPlot=="GAMMAS" || fRankingPlot=="ALL") {
            const auto container = fitTool.GetRankingContainer(fPOIs.at(ipoi), exclude, false, true);
            manager.SetRankingContainer(container);
            manager.SetSuffix(fSuffix+"_"+fPOIs.at(ipoi)+"_Breakdown_gammas");
            manager.PlotRanking(fRegions, fNormFactorNames, fShapeFactorNames, false, true);
        }
    }

}

//__________________________________________________________________________________
//
void TRExFit::SetBlindingInRegions(const std::string& wsPath, const bool hasPostfit) {

    std::string frFile("");
    if (hasPostfit) {
        if (fFitResultsRootFile != "") {
            frFile = fFitResultsRootFile;
        } else {
            frFile = fName + "/Fits/"+fInputName+fSuffix+".root";
        }
    }

    xRooNode rooNodePrefit = Common::ReadWSandFitResults(wsPath, "");

    auto pdfPrefit = rooNodePrefit["simPdf"];
    if (!pdfPrefit) {
        LOG(ERROR) << "Cannot read simPdf for prefit\n";
        exit(EXIT_FAILURE);
    }

    xRooNode rooNodePostfit{};
    std::shared_ptr<xRooNode> pdfPostfit(nullptr);
    if (hasPostfit) {
        rooNodePostfit = Common::ReadWSandFitResults(wsPath, frFile);
        pdfPostfit = rooNodePostfit["simPdf"];
        if (!pdfPostfit) {
            LOG(ERROR) << "Cannot read simPdf for postfit\n";
            exit(EXIT_FAILURE);
        }
    }

    /*
        We can set which bins are blinded for each region once
        and for all to be used in all plots and tables
    */
    for (const auto& ireg : fRegions) {
        // Blinded bins according to blinding threshold
        auto regionPrefit = pdfPrefit->find(ireg->fName);
        if (!regionPrefit) {
            LOG(ERROR) << "Cannot find channel: " << ireg->fName << " in the prefit WS\n";
            exit(EXIT_FAILURE);
        }
        std::vector<int> computedBlindedBins = Common::GetBlindedBins(ireg.get(), regionPrefit, fBlindingType, fBlindingThreshold);
        ireg->fBlindedBins = Common::CombineVectors(ireg->fManualBlindBins, computedBlindedBins, true);
        ireg->fComputedBlindedBins = computedBlindedBins;
        if (hasPostfit) {
            auto regionPostfit = pdfPostfit->find(ireg->fName);
            if (!regionPostfit) {
                LOG(ERROR) << "Cannot find channel: " << ireg->fName << " in the postfit WS\n";
                exit(EXIT_FAILURE);
            }
            // Blinded bins post-fit, from threshold and those manually blinded
            computedBlindedBins = Common::GetBlindedBins(ireg.get(), regionPostfit, fBlindingType, fBlindingThreshold);
            ireg->fBlindedBinsPostFit = Common::CombineVectors(ireg->fManualBlindBins, computedBlindedBins, true);
        }
    }
}

//__________________________________________________________________________________
//
void TRExFit::TranslateWSToHS3(const std::string& wsFilePath,
                               const std::string& wsPath,
                               const std::string& outPath) const {

    LOG(INFO) << "Will produce the workspace in HS3 format\n";

    std::unique_ptr<TFile> wsFile(TFile::Open(wsFilePath.c_str(), "READ"));
    if (!wsFile) {
        LOG(ERROR) << "Cannot open WS file: " << wsFilePath << ". Did you run the \"w\" step before?\n";
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<RooWorkspace> ws(wsFile->Get<RooWorkspace>(wsPath.c_str()));
    if (!ws) {
        LOG(ERROR) << "Cannot read ws from file: " << wsPath << "\n";
        exit(EXIT_FAILURE);
    }

    bool isRegularisedUnfolding(false);
    if (fFitType == TRExFit::FitType::UNFOLDING) {
        for (const auto& inf : fNormFactors) {
            if (inf->fTau > 0) {
                isRegularisedUnfolding = true;
                break;
            }
        }
    }

    if (isRegularisedUnfolding) {
        ws->set("externalConstraints");
    }

    RooJSONFactoryWSTool tool(*ws);

    LOG(INFO) << "Writing the workspace into HS3 JSON format: " << outPath << "\n";
    tool.exportJSON(outPath);

    wsFile->Close();
}

//__________________________________________________________________________________
//
void TRExFit::ReadStatOnlyFromErrorDecomposition(const std::string& poi, const std::string& path) {
    LOG(DEBUG) << "Reading error decomposition from: " << path << "\n";
    std::ifstream in(path.c_str());
    if (!in.is_open() || !in.good()) {
        LOG(ERROR) << "Cannot open file: " << path << "\n";
        return;
    }

    auto itr     = fStatOnlyErrorDecomposition.find(poi);
    auto itrUp   = fStatOnlyErrorDecompositionUp.find(poi);
    auto itrDown = fStatOnlyErrorDecompositionDown.find(poi);

    std::string name;
    double error;
    double errorUp;
    double errorDown;

    while (in >> name >> error >> errorUp >> errorDown) {
        if (name != "STAT_ERROR") continue;

        if (itr == fStatOnlyErrorDecomposition.end()) {
            fStatOnlyErrorDecomposition.insert(std::make_pair(poi, error));
        } else {
            itr->second = error;
        }
        if (itrUp == fStatOnlyErrorDecompositionUp.end()) {
            fStatOnlyErrorDecompositionUp.insert(std::make_pair(poi, errorUp));
        } else {
            itrUp->second = errorUp;
        }
        if (itrDown == fStatOnlyErrorDecompositionDown.end()) {
            fStatOnlyErrorDecompositionDown.insert(std::make_pair(poi, errorDown));
        } else {
            itrDown->second = errorDown;
        }

        in.close();
        return;
    }

    // should not happen
    LOG(ERROR) << "Cannot find POI: " << poi << "\n";

    in.close();
}

//__________________________________________________________________________________
//
void TRExFit::ReadStatOnlyFromErrorDecomposition(const std::vector<std::string>& pois, const std::string& path) {
    for (const auto& ipoi : pois) {
        const std::string fullPath = path + "_" + ipoi + ".txt";
        this->ReadStatOnlyFromErrorDecomposition(ipoi, fullPath);
    }
}

//__________________________________________________________________________________
//
double TRExFit::GetStatOnlyErrorFromDecomposition(const std::string& poi) const {
    auto itr = fStatOnlyErrorDecomposition.find(poi);
    if (itr == fStatOnlyErrorDecomposition.end()) {
        LOG(ERROR) << "Cannot find POI: " << poi << "\n";
        return -1;
    }

    return itr->second;
}

//__________________________________________________________________________________
//
double TRExFit::GetStatOnlyErrorFromDecompositionUp(const std::string& poi) const {
    auto itr = fStatOnlyErrorDecompositionUp.find(poi);
    if (itr == fStatOnlyErrorDecompositionUp.end()) {
        LOG(ERROR) << "Cannot find POI: " << poi << "\n";
        return -1;
    }

    return itr->second;
}

//__________________________________________________________________________________
//
double TRExFit::GetStatOnlyErrorFromDecompositionDown(const std::string& poi) const {
    auto itr = fStatOnlyErrorDecompositionDown.find(poi);
    if (itr == fStatOnlyErrorDecompositionDown.end()) {
        LOG(ERROR) << "Cannot find POI: " << poi << "\n";
        return -1;
    }

    return itr->second;
}

//__________________________________________________________________________________
//
void TRExFit::PlotUnfoldingErrors(const FittingTool& fitTool, const std::vector<std::string>& categories) {
    auto isNorm =[](const std::vector<std::unique_ptr<Unfolding> >& unfolding) {
        for (const auto& iunf : unfolding) {
            if (iunf->fUnfoldNormXSec) return true;
        }
        return false;
    };

    this->ProduceSystSubCategoryMap();

    gSystem->mkdir((fName+"/UnfoldingPlots").c_str());
    UncertaintyPlotter plotter{};
    plotter.SetOutputFolder(fName+"/UnfoldingPlots");
    plotter.SetUseTotalUncertainty(fUnfoldingUncertaintyBreakdownTotal);
    plotter.SetCMELabel(fCmeLabel);
    plotter.SetPlotLabel(fPlotLabel);
    plotter.SetLumiLabel(fLumiLabel);
    plotter.SetIsNormalized(isNorm(fUnfolding));
    plotter.SetCategories(categories);
    plotter.SetSubCategoryMap(fSubCategoryImpactMap);
    for (const auto& iunfolding : fUnfolding) {
        if (iunfolding->fBreakdownRange < std::numeric_limits<double>::max()) { // is set
            plotter.SetCustomRange(iunfolding->fBreakdownRange);
        } else {
            plotter.SetUseCustomRange(false);
        }
        std::vector<int> reparametrisedBins;
        std::vector<std::pair<std::string, int> > pois;
        std::vector<std::string> otherPOIs;
        for (int ibin = 0; ibin < iunfolding->fNumberUnfoldingTruthBins; ++ibin) {
            if (iunfolding->fUnfoldNormXSec && iunfolding->fUnfoldNormXSecBinN == ibin+1) {
                reparametrisedBins.emplace_back(ibin+1);
                continue;
            }
            const std::string name = iunfolding->fName + "_Bin_" + Common::IntToFixLenStr(ibin+1) + "_mu";
            auto itr = std::find_if(fNormFactors.begin(), fNormFactors.end(), [&name](const auto& element){return element->fName == name;});
            if (itr == fNormFactors.end()) {
                LOG(ERROR) << "Cannot find NormFactor: " << name << "\n";
                return;
            }
            if (!((*itr)->fExpression.first.empty())) {
                reparametrisedBins.emplace_back(ibin+1);
                continue;
            }
            pois.emplace_back(std::make_pair(name, ibin));
        }
        for (const auto& ipoi : fPOIs) {
            if (Common::IsUnfoldingPOI(ipoi, fUnfolding)) continue;
            auto itr = std::find_if(pois.begin(), pois.end(), [&ipoi](const auto& element){return element.first == ipoi;});
            if (itr == pois.end()) {
                otherPOIs.emplace_back(ipoi);
            }
        }
        plotter.SetReparametrisedBins(reparametrisedBins);
        plotter.SetPOIs(pois);
        plotter.SetOtherPOIs(otherPOIs);
        plotter.SetXaxisLabel(iunfolding->fUnfoldingTitleX);
        std::unique_ptr<TFile> input(TFile::Open((fName + "/UnfoldingHistograms/FoldedHistograms.root").c_str(), "READ"));
        if (!input) {
            LOG(ERROR) << "Cannot read file from " << fName << "/UnfoldingHistograms/FoldedHistograms.root\n";
            exit(EXIT_FAILURE);
        }

        std::unique_ptr<TH1> truth(dynamic_cast<TH1*>(input->Get((iunfolding->fName+"_truth_distribution").c_str())));
        if (!truth) {
            LOG(ERROR) << "Cannot read the truth distribution\n";
            exit(EXIT_FAILURE);
        }
        truth->SetDirectory(nullptr);

        plotter.PlotUncertainties(truth.get(),
                                  *iunfolding,
                                  fitTool,
                                  fName+"/RooStats/"+fInputName+"_combined_"+fInputName+fSuffix+"_model.root",
                                  "combined",
                                  fName+"/Fits/"+fInputName+fSuffix+".root", fNormFactors);
        input->Close();
    }
}

//__________________________________________________________________________________
//
void TRExFit::ProcessEFTShapeFactors() {
    if (fFitType != TRExFit::FitType::EFT) return;
    if (fEFTConfig.GetSplitSamplesPerBin()) return;

    std::vector<std::string> referenceSamples;
    for (const auto& isample : fSamples) {
        if (isample->fType != Sample::SampleType::EFT) continue;
        const std::string ref = isample->fEFTSMReference;
        auto itr = std::find(referenceSamples.begin(), referenceSamples.end(), ref);
        if (itr == referenceSamples.end()) {
            referenceSamples.emplace_back(ref);
        }
    }

    std::string filePath = fName;
    filePath += "/EFT/SF_";
    filePath += (fEFTConfig.GetEFTOrder() == EFTConfig::EFTOrder::LIN ? "LINEAR" : "ALL");
    filePath += "_EFT_Fit_results.txt";
    LOG(INFO) << "Processing EFT shape factors from file: " << filePath << "\n";
    ReparametrizationManager manager{};
    manager.ReadReparametrizationFile(filePath);
    for (const auto& isample : fSamples) {
        if (isample->fType == Sample::SampleType::EFT) continue;
        auto itr = std::find(referenceSamples.begin(), referenceSamples.end(), isample->fName);
        if (itr == referenceSamples.end()) continue;
        LOG(DEBUG) << "Reparametrizing sample: " << isample->fName << "\n";
        for (const auto& isf : isample->fShapeFactors) {
            const auto reparam = manager.GetParametrization(isf->fRegions.at(0), isample->fName);

            const std::vector<std::pair<std::string, std::string> > updatedReparametrization = fEFTConfig.GetUpdatedParametrization(reparam, isf->fRegions.at(0), isample->fName);

            isf->fExpression = updatedReparametrization;
        }
    }
}
