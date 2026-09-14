#ifndef REGION_H
#define REGION_H

/// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/PruningUtil.h"
#include "TRExFitter/TRExFit.h"
#include "TRExFitter/TRExPlot.h"

/// c++ includes
#include <map>
#include <string>
#include <vector>

/// Forwards class declaration
class CorrelationMatrix;
class FitResults;
class TFile;
class TGraphAsymmErrors;
class TH1;
class THStack;
class NormFactor;
class Sample;
class SampleHist;
class ShapeFactor;
class Systematic;
class TRExFit;
class TRExPlot;

class Region {
public:

    enum RegionType {
        CONTROL = 1,
        VALIDATION = 2,
        SIGNAL = 3
    };

    enum DataType {
        REALDATA = 0,
        ASIMOVDATA = 1
    };

   explicit Region(const std::string& name);
    ~Region() = default;
    Region(const Region& r) = delete;
    Region(Region&& r) = delete;
    Region& operator=(const Region& r) = delete;
    Region& operator=(Region&& r) = delete;

    // -------
    // Methods
    // -------

    /**
     * @brief Set the sampleHist - calls the constructor of SampleHist and that reads the histograms
     *
     * @param sample Sample object
     * @param histoName name of the histogram
     * @param fileName file path
     * @param structure if nullptr, no folder structure will be used
     * @return std::shared_ptr<SampleHist>
     */
    std::shared_ptr<SampleHist> SetSampleHist(Sample *sample,
                                              const std::string& histoName,
                                              const std::string& fileName,
                                              const bool& dropBinHistoNeeded,
                                              const bool& checkHistogram,
                                              const Common::FolderStructure* structure = nullptr);
    std::shared_ptr<SampleHist> SetSampleHist(Sample *sample, TH1* hist );
    std::shared_ptr<SampleHist> GetSampleHist(const std::string &sampleName) const;

    void SavePreFitUncertaintyAndTotalMCObjects();

    /**
     * @brief Build the error histogram
     *
     * @param channel node for the channel
     * @param nominal histogram to ge the binning
     * @param fitRes fit results for the chi^2 calculation
     * @param isPostFit flag to tell the code if it is pre-fit or post-fit
     */
    void BuildPrePostFitErrorHistxRooFit(const std::shared_ptr<xRooNode>& channel,
                                         const TH1* nominal,
                                         const FitResults* fitRes,
                                         const bool isPostFit);

    /**
     * @brief Draw pre/post fit plots
     *
     * @param pdf the simPdf object from the WS
     * @param pdfPrefit neede for PREFITONPOSTFIT
     * @param fitRes fit results
     * @param canvasSize size of the canvas
     * @param isPostFit is postfit plot
     * @param isBonly is BONLY fit
     * @param opt more options
     * @return std::shared_ptr<TRExPlot>
     */
    std::shared_ptr<TRExPlot> DrawPrePostFit(const std::shared_ptr<xRooNode>& pdf,
                                             const std::shared_ptr<xRooNode>& prefPrefit,
                                             const FitResults* fitRes,
                                             const std::vector<int>& canvasSize,
                                             const bool isPostFit,
                                             const bool isBonly,
                                             std::string opt="");

    void SetBinning(int N, const std::vector<double>& bins);
    void Rebin(int N);
    void SetRebinning(int N, const std::vector<double>& bins);
    void SetRegionType(RegionType type);
    void SetRegionDataType( DataType type );
    void AddSample(Sample *sample);

    void AddSelection(const std::string& selection);
    void AddMCweight(const std::string& weight);
    void SetVariable(const std::string& variable,int nbin,double xmin,double xmax,const std::string& corrVar1="", const std::string& corrVar2="");
    void SetAlternativeVariable(const std::string& variable, const std::string& sample);
    bool UseAlternativeVariable(const std::string& sample);
    std::string GetAlternativeVariable(const std::string& sample) const;
    void SetAlternativeSelection(const std::string& selection, const std::string& sample);
    bool UseAlternativeSelection(const std::string& sample);
    std::string GetAlternativeSelection(const std::string& sample) const;

    void AddSystematic(Systematic *syst);

    // cosmetics
    void SetVariableTitle(const std::string& name);
    void SetLabel(const std::string& label, const std::string& shortLabel="");

    // log
    void Print() const;

    /**
     * @brief Print effect of normalisation from each systematic source
     *
     * @param pdf the simPdf object from the WS
     * @param opt other options
     */
    void PrintSystTable(std::shared_ptr<xRooNode> pdf, const std::string& opt="") const;

    /**
     * Function that calls systematics pruning through the PruningUtil class
     * @param pointer to PruningUtil instance
     */
    std::vector<std::pair<std::string, double> > SystPruning(PruningUtil *pu);

    /**
      * Helper function to get a "total prediction" histogram
      * @param bool specifying whether signal sample have to be included in the sum or not (true by default)
      * @return combined histogram
      */
    std::unique_ptr<TH1> GetTotHist(bool includeSignal);

    /**
     * @brief Set the Automatic Drop Bins object
     *
     * @return true
     * @return false
     */
    void SetAutomaticDropBins(const bool flag) {fAutomaticDropBins = flag;}

    /**
     * @brief Get the Automatic Drop Bins object
     *
     * @return true
     * @return false
     */
    bool GetAutomaticDropBins() const {return fAutomaticDropBins;}

    /**
     * @brief Get the Nbins object
     *
     * @return int
     */
    int GetNbins() const {return fNbins;}

    /**
     * @brief Store the bin numbers in a file
     *
     */
    void AddBinsToFile() const;

    /**
     * @brief Set the Nbins object
     *
     * @param n
     */
    void SetNbins(const int n) {fNbins = n;}

    /**
     * @brief Set the And Add Nbins object
     *
     * @param n
     * @param addToFile
     */
    void SetAndAddNbins(const int n, const bool addToFile);

    /**
     * @brief Get the integral for a given sample
     *
     * @param sampleName
     * @param isPostfit
     * @return double
     */
    double GetIntegral(const std::string& sampleName, const bool isPostfit) const;

    /**
     * @brief Get the Integral for bOnly for normalised
     *
     * @param sampleName
     * @return double
     */
    double GetIntegralBonlyNorm(const std::string& sampleName) const;

    // -------
    // Members
    // -------

    std::string fName;
    bool fUseFriend;
    std::string fVariableTitle;
    std::string fYTitle;
    std::string fLabel; // something like "e/mu + 6 j, >=4 b b"
    std::string fShortLabel; // something like "6j,3b"
    std::string fTexLabel;
    std::string fFitName;
    std::string fPlotSubdir;
    RegionType fRegionType;
    DataType fRegionDataType;
    bool fHasData;
    std::shared_ptr<SampleHist> fData;
    bool fHasSig;
    std::vector<std::shared_ptr<SampleHist> > fSig;
    std::vector<std::shared_ptr<SampleHist> > fBkg;
    std::vector < std::shared_ptr<SampleHist> > fSampleHists;
    std::vector < std::shared_ptr<Sample> > fSamples;
    double fYmaxScale;
    double fYmin;
    double fYmax;
    double fRatioYmin;
    double fRatioYmax;
    double fRatioYminPostFit;
    double fRatioYmaxPostFit;
    std::string fRatioYtitle;
    TRExPlot::RATIOTYPE fRatioType;

    // to draw
    std::unique_ptr<TH1> fTot;
    std::unique_ptr<TGraphAsymmErrors> fErr;

    // post fit
    std::unique_ptr<TH1> fTot_postFit;
    std::unique_ptr<TH1> fBkgOnly_postFit;
    std::unique_ptr<TGraphAsymmErrors> fErr_postFit;
    std::shared_ptr<TGraphAsymmErrors> fErrBkgOnly_postFit;

    // ntuple stuff
    std::string fBinTransfo;
    double fTransfoDzBkg;
    double fTransfoDzSig;
    double fTransfoErr;
    double fTransfoFzBkg;
    double fTransfoFzSig;
    double fTransfoJpar1;
    double fTransfoJpar2;
    double fTransfoJpar3;
    std::vector<std::string> fAutoBinBkgsInSig;
    std::string fVariable;
    std::map<std::string, std::string> fAlternativeVariables;
    std::map<std::string, std::string> fAlternativeSelections;
    std::string fCorrVar1;
    std::string fCorrVar2;
    double fXmin;
    double fXmax;
    std::string fSelection;
    std::string fMCweight;
    std::vector<std::string> fNtuplePaths;
    std::vector<std::string> fNtuplePathSuffs;
    std::vector<std::string> fFriendPaths;
    std::vector<std::string> fFriendPathSuffs;
    std::vector<std::string> fNtupleFiles;
    std::vector<std::string> fNtupleFileSuffs;
    std::vector<std::string> fNtupleNames;
    std::vector<std::string> fNtupleNameSuffs;

    // histogram stuff
    std::vector<double> fHistoBins;
    int fHistoNBinsRebin;
    std::vector<double> fHistoBinsPost;
    int fHistoNBinsRebinPost;
    std::vector<std::string> fResponseMatrixPaths;
    std::vector<std::string> fResponseMatrixPathSuffs;
    std::vector<std::string> fResponseMatrixFiles;
    std::vector<std::string> fResponseMatrixFileSuffs;
    std::vector<std::string> fResponseMatrixNames;
    std::vector<std::string> fResponseMatrixNameSuffs;
    std::vector<std::string> fAcceptancePaths;
    std::vector<std::string> fAcceptancePathSuffs;
    std::vector<std::string> fAcceptanceFiles;
    std::vector<std::string> fAcceptanceFileSuffs;
    std::vector<std::string> fAcceptanceNames;
    std::vector<std::string> fAcceptanceNameSuffs;
    std::vector<std::string> fSelectionEffPaths;
    std::vector<std::string> fSelectionEffPathSuffs;
    std::vector<std::string> fSelectionEffFiles;
    std::vector<std::string> fSelectionEffFileSuffs;
    std::vector<std::string> fSelectionEffNames;
    std::vector<std::string> fSelectionEffNameSuffs;
    std::vector<std::string> fMigrationPaths;
    std::vector<std::string> fMigrationPathSuffs;
    std::vector<std::string> fMigrationFiles;
    std::vector<std::string> fMigrationFileSuffs;
    std::vector<std::string> fMigrationNames;
    std::vector<std::string> fMigrationNameSuffs;
    std::vector<std::string> fHistoPaths;
    std::vector<std::string> fHistoPathSuffs;
    std::vector<std::string> fHistoFiles;
    std::vector<std::string> fHistoFileSuffs;
    std::vector<std::string> fHistoNames;
    std::vector<std::string> fHistoNameSuffs;

    bool fUseStatErr;

    int fIntCode_overall;
    int fIntCode_shape;

    std::vector< std::string > fSystNames;
    std::vector< std::string > fNpNames;

    TRExFit::FitType fFitType;
    std::vector<std::string> fPOIs;
    std::string fFitLabel;

    std::string fLumiLabel;
    std::string fCmeLabel;

    double fLumiScale;

    bool fLogScale;
    bool fLogScaleX;

    double fBinWidth;

    double fBlindingThreshold;
    Common::BlindingType fBlindingType;

    bool fSkipSmoothing;

    std::string fPlotLabel;
    std::string fSuffix;

    std::string fGroup; // used to split yield tables

    bool fKeepPrefitBlindedBins;
    int fGetChi2;

    std::vector<int> fDropBins;
    std::vector<int> fManualBlindBins;
    std::vector<int> fComputedBlindedBins;
    std::vector<int> fBlindedBins;
    std::vector<int> fBlindedBinsPostFit;

    std::vector<std::string> fBinLabels;
    std::vector<int> fBinDividers;
    std::vector<std::string> fRangeLabels;

    // two parameters to chose put tick marks at specific bin edges and how to label them
    std::vector<int> fBinTickMarks;
    std::vector<std::string> fBinTickMarksLabels;

    double fChi2val;
    int fNDF;
    double fChi2prob;

    bool fUseGammaPulls;

    std::vector<double> fXaxisRange;

    double fLabelX;
    double fLabelY;
    double fLegendX1;
    double fLegendX2;
    double fLegendY;

    int fLegendNColumns;

    std::map<std::string,int> fIsBinOfRegion;

    int fNumberUnfoldingRecoBins;
    bool fNormalizeMigrationMatrix;
    bool fHasAcceptance;

    std::string fFolder;
    bool fHEPDataFormat;
    bool fUheppFormat;
    std::unordered_map<std::string, std::pair<double,double> > m_randomNFmap;
    bool fDrawDataMinusBkg;
    bool fDrawDataMinusBkgOnSignal;
    std::vector<std::string> fPruningPlotSamples;
    std::vector<std::shared_ptr<NormFactor> > fNormFactors;
    std::string fJobName;
    std::string fJobBinningsPath;
    std::map<std::string, double> fPrefitYields;
    std::map<std::string, double> fPostfitYields;
    std::map<std::string, double> fPostfitYieldsBonlyNorm;
    bool fNoPrePostFitLabel;
    std::string fPreFitLabel;
    std::string fPostFitLabel;
    int fToysForErrorBand;

private:

    bool fAutomaticDropBins;
    int fNbins;

    /**
     * @brief Get the Chi 2 agreement between histograms
     *
     * @param isPostFit Flag to separate pre- and post-fit
     * @param h_data Data histogram
     * @param h_nominal Total prediction
     * @param h_up Shifted total histogram per NP
     * @param systNames Names of the NPs
     * @param matrix Post-fit correlation matrix
     * @return std::pair<double,int>
     */
    std::pair<double,int> GetChi2Test(const bool isPostFit,
                                      const TH1* h_data,
                                      const TH1* h_nominal,
                                      const std::vector< std::shared_ptr<TH1> >& h_up,
                                      const std::vector< std::string >& systNames,
                                      const CorrelationMatrix *matrix=nullptr) const;

};

#endif
