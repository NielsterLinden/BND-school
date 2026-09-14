#ifndef MULTIFIT_H
#define MULTIFIT_H

/// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/TRExFit.h"

/// c++ includes
#include <map>
#include <memory>
#include <vector>

/// Forwards declaration
class ConfigParser;
class FittingTool;
class RooAbsData;
class RooDataSet;
class RooSimultaneous;
class RooWorkspace;
class TH1D;

class MultiFit {
public:

    explicit MultiFit(const std::string& name);

    ~MultiFit() = default;
    MultiFit(const MultiFit& m) = delete;
    MultiFit(MultiFit&& m) = delete;
    MultiFit& operator=(const MultiFit& m) = delete;
    MultiFit& operator=(MultiFit&& m) = delete;

    void AddFitFromConfig(const std::string& configFile,
                          const std::string& opt,
                          const std::string& options,
                          const std::string& label,
                          const std::string& loadSuf,
                          const std::string& wsFile,
                          const bool useInFit,
                          const bool useInComparison);

    void AddPOI(const std::string& name);

    std::unique_ptr<RooWorkspace> CombineWS() const;
    void SaveCombinedWS() const;
    std::map < std::string, double > FitCombinedWS( int fitType, const std::string& inputData, bool doLHscanOnly ) const;
    void RunLimitScan(xRooNLLVar::xRooHypoSpace& hs) const;
    void RunSignificanceScan(xRooNLLVar::xRooHypoSpace& hs) const;
    void GetCombinedLimit(const std::string& inputData="obsData"); // or asimovData
    void GetCombinedSignificance(const std::string& inputData="obsData"); // or asimovData

    void ComparePOI(const std::string& POI, const std::size_t index) const;
    void CompareEFT() const;
    void CompareLimit();
    void ComparePulls(const std::string& category="") const;
    void CompareNormFactors(const std::string& category="") const;
    void PlotCombinedCorrelationMatrix() const;

    /**
     * @brief Produce the NP ranking txt file
     *
     * @param NPnames flag to tell which NPs to consider
     */
    void ProduceNPRanking(const std::string& NPnames) const;

    /**
     * @brief Code that decides whether the ranking or plotting or both is run
     *
     */
    void PlotNPRankingManager() const;

    /**
     * @brief Plot the NP ranking
     *
     * @param flagSysts which systematics?
     * @param flagGammas what about gammas?
     */
    void PlotNPRanking(bool flagSysts=true, bool flagGammas=false) const;
    void PlotSummarySoverB() const;
    void GetLikelihoodScan(const RooWorkspace *ws, const std::string& varName, RooDataSet* data,bool recreate=true) const;
    void Get2DLikelihoodScan(const RooWorkspace *ws, const std::vector<std::string>& varName, RooDataSet* data) const;
    void BuildGroupedImpactTable() const;

    /**
      * A helper function to get vector of unique normfactors used in a fit
      * @return the vector of norm factors
      */
    std::vector<std::shared_ptr<NormFactor> > GetFitNormFactors() const;

    /**
      * A helper function to get vector of unique shapefactors used in a fit
      * @return the vector of shape factors
      */
    std::vector<std::shared_ptr<ShapeFactor> > GetFitShapeFactors() const;

    /**
      * A helper function to get vector of unique systematics used in a fit
      * @return the vector of systematics
      */
    std::vector<std::shared_ptr<Systematic> > GetFitSystematics() const;

    /**
      * A helper function to get map of fixed NPs and their values from the individual configs
      * @return the map
      */
    std::map<std::string, double> GetFixedNPs() const;

    /**
      * A helper function to get vector of regions
      * @return the vector
      */
    std::vector<Region* > GetFitRegions() const;

    /**
     * @brief A helper function to get the vector of unfoldings
     *
     * @return std::vector<Unfolding*>
     */
    std::vector<std::unique_ptr<Unfolding> > GetUnfoldings() const;

    /**
      * @brief Function to run Limit estimation using toys
      * @param Data
      * @param Workspace
      */
    void RunLimitToys(RooAbsData* data, RooWorkspace* ws) const;

    /**
     * @brief Run the significance from toys
     *
     * @param data
     * @param ws
     */
    void RunSignificanceToys(RooAbsData* data, RooWorkspace* ws) const;

    /**
     * @brief Plots the unfolded data
     *
     */
    void PlotUnfoldedData();

    /**
     * @brief Draw gamma plots
     *
     */
    void DrawGammas();

    inline void SetLimitType(const TRExFit::LimitType type){fLimitType = type;}

    inline void SetSignificanceType(const TRExFit::SignificanceType type){fSignificanceType = type;}

    std::unique_ptr<TH1D> Combine(const std::vector<std::unique_ptr<TH1D> >& hists) const;
    std::unique_ptr<TH1D> Rebin(TH1D* h, const std::vector<double>& vec, bool isData=true) const;
    std::vector<std::pair<std::string,double> > GetSmallNPs() const;

    /**
     * @brief Add the blinded params form the individual configs
     *
     */
    void AdjustBlindedParams();

    /**
     * @brief Pre0process tep needed for reparametrised shape factors and template morphing
     *
     */
    void ProcessShapeFactorReparametrisation();

    /**
     * @brief Check that regions have the same name
     *
     */
    void CheckSameRegionNames() const;

    /**
     * @brief Plot ranking from covariance matrix
     *
     * @param fitTool
     */
    void PlotRankingFromCovMatrix(const FittingTool& fitTool) const;

    /**
     * @brief Get the Subcategory Map
     *
     * @return std::map<std::string, std::string>
     */
    std::map<std::string, std::string> GetSubcategoryMap() const;

    /**
     * @brief Take the WS and store it as HS3 format
     *
     */
    void TranslateWSToHS3() const;

    /**
     * @brief Plot the uncertainties in each bin
     *
     * @param fitTool
     * @param categories
     */
    void PlotUnfoldedErrors(const FittingTool& fitTool, const std::vector<std::string>& categories) const;

    /**
     * @brief Run pseudoexperiments
     *
     */
    void RunToys() const;

    /**
     * @brief Get the combined list of varaibles to run MINOS on
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> GetMinos() const;

    /**
     * @brief Get the Regularization Type object
     *
     * @return int
     */
    int GetRegularizationType() const;

    /**
     * @brief Check if all the datasets are provided in all WS (or none) to prevent a crash in HistFactory
     *
     */
    void CheckDatasetConsistency() const;

    std::vector<std::string> fFitNames;
    std::vector<std::unique_ptr<TRExFit> > fFitList;
    std::vector<std::string> fFitLabels;
    std::vector<std::string> fFitSuffs;
    std::vector<std::string> fWsFiles;
    std::vector<std::string> fDirectory;
    std::vector<std::string> fInputName;

    std::vector<std::string> fNPCategories;

    bool fCombine;
    bool fCompare;
    bool fStatOnly;
    bool fIncludeStatOnly;

    bool fCompareLimits;
    bool fComparePOI;
    bool fCompareEFT;
    bool fComparePulls;
    bool fPlotCombCorrMatrix;

    std::string fName;
    std::string fDir;
    std::string fOutDir;
    std::string fLabel;
    bool fShowObserved;
    std::string fLimitTitle;
    std::vector<std::string> fPOITitle;
    std::string fRankingOnly;
    std::string fGroupedImpactCategory;

    std::vector<std::string> fPOIs;
    std::string fPOIforLimit;
    std::string fPOIforSig;
    std::vector<double> fPOIMin;
    std::vector<double> fPOIMax;
    std::vector<std::string> fPOIPrecision;
    double fLimitMax;

    bool fUseRnd;
    double fRndRange;
    long int fRndSeed;

    std::string fLumiLabel;
    std::string fCmeLabel;
    std::string fCombiLabel;

    std::unique_ptr<ConfigParser> fConfig;

    std::string fSaveSuf;

    std::string fDataName;
    int fFitType;

    bool fFastFit;
    bool fFastFitForRanking;
    std::string fNuisParListFile;

    bool fPlotSoverB;
    std::string fSignalTitle;

    std::string fFitResultsRootFile;
    std::string fLimitsFile;
    std::vector<std::string> fLimitsFiles;
    std::string fBonlySuffix;

    bool fShowSystForPOI;

    std::vector<std::string> fVarNameLH;
    std::vector<std::vector<std::string> > fVarName2DLH;
    double fLHscanMin;
    double fLHscanMax;
    int fLHscanSteps;
    int fLHscanStep;
    double fLHscanMinY;
    double fLHscanMaxY;
    int fLHscanStepsY;
    int fLHscanStepY;
    bool fDoGroupedSystImpactTable;

    std::vector<std::string> fPOIName;
    std::vector<double> fPOINominal;
    std::vector<double> fPOIAsimov;

    //
    // Limit parameters
    //
    bool fLimitIsBlind;
    bool fSignalInjection;
    double fSignalInjectionValue;
    std::string fLimitParamName;
    double fLimitParamValue;
    std::string fLimitOutputPrefixName;
    double fLimitsConfidence;
    bool fLimitUseAutoDiff;

    //
    // Significance parameters
    //
    bool fSignificanceIsBlind;
    bool fSignificanceDoInj;
    double fSignificancePOIAsimov;
    std::string fSignificanceParamName;
    double fSignificanceParamValue;
    std::string fSignificanceOutputPrefixName;
    bool fSignificanceUseAutoDiff;

    bool fShowTotalOnly;

    std::vector<std::pair<std::string, double> > fPOIInitials;
    bool fHEPDataFormat;
    bool fUheppFormat;
    std::vector<std::string> fConfigPaths;
    int fFitStrategy;
    int fCPU;
    bool fBinnedLikelihood;
    bool fUsePOISinRanking;
    bool fUseHesseBeforeMigrad;
    bool fUseNllInLHscan;
    int fLimitToysStepsSplusB;
    int fLimitToysStepsB;
    int fLimitToysScanSteps;
    double fLimitToysScanMin;
    double fLimitToysScanMax;
    int fLimitToysSeed;
    bool fLimitPlot;
    bool fLimitFile;
    TRExFit::LimitType fLimitType;
    int fLimitFitStrategy;
    TRExFit::ToysUsexRooFit fLimitToysUsexRooFit;
    TRExFit::SignificanceType fSignificanceType;
    int fSignificanceFitStrategy;
    int fSignificanceToysStepsSplusB;
    int fSignificanceToysStepsB;
    bool fSignificancePlot;
    int fSignificanceToysSeed;
    bool fSignificanceToysUsexRooFit;
    std::vector<std::string> fOnlyRegions;
    bool fSpeedUpFit;
    double fNPCutOff;
    bool fPlotUnfoldedData;
    bool fUnfoldingShowStat;
    bool fUnfoldingShowUncertaintyBreakdown;
    bool fUnfoldingUncertaintyBreakdownTotal;
    bool fPlotGammas;
    std::vector<std::string> fBlindedParams;
    std::vector<std::string> fRankingPOIName;
    std::vector<int> fRankingUpperAxisNdivision;
    std::vector<double> fRankingPOIAxisScale;
    bool fHasShapeFactorReparametrisation;
    bool fUseHesse;
    int fMaximumNumberFCNcalls;
    double fToleranceScale;
    bool fShiftGlobalObservablesInRanking;
    bool fUseAutoDiff;
    std::string fNLLOffset;

    // toys
    int fFitToys;
    int fToysHistoNbins;
    int fToysSeed;
    bool fToysSetStatOnlyFluctuation;
    std::map<std::string, std::pair<double, double> > fToysRandomPOIStartingValues;
    bool fToysLHScanForAll;
    std::map<std::string, std::pair<double, double> > fToysLHScanCondition;

};

#endif
