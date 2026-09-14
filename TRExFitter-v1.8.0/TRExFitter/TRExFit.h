#ifndef TRExFIT_H
#define TRExFIT_H

/// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/EFTConfig.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Morphing.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/TemplateMorpher.h"
#include "TRExFitter/TRExPlot.h"
#include "TRExFitter/YamlConverter.h"

// Unfolding includes
#include "UnfoldingCode/UnfoldingCode/FoldingManager.h"

// RooFit
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"
#include "RooMultiVarGaussian.h"

// RooStats
#include "RooStats/HistFactory/Measurement.h"

/// c++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

/// Forwards class declaration
class ConfigParser;
class CorrelationMatrix;
class FitResults;
class FittingTool;
class NormFactor;
class RooDataSet;
class RooWorkspace;
class Region;
class Sample;
class SampleHist;
class ShapeFactor;
class TColor;
class TGraphAsymmErrors;
class TruthSample;
class TFile;
class Unfolding;
class UnfoldingSample;
class UnfoldingSystematic;
class UnfoldingResult;

class TRExFit {
public:

    enum FitType {
        UNDEFINED = 0,
        SPLUSB = 1,
        BONLY = 2 ,
        UNFOLDING = 3,
        EFT = 4
    };

    enum FitRegion {
        CRONLY = 1,
        CRSR = 2
    };

    enum InputType {
        HIST = 0,
        NTUP = 1
    };

    enum class LimitType {
        ASYMPTOTIC = 0,
        TOYS = 1
    };

    enum class SignificanceType {
        ASYMPTOTIC = 0,
        TOYS = 1
    };

    enum TemplateInterpolationOption{
        LINEAR = 0,
        SMOOTHLINEAR = 1,
        SQUAREROOT = 2
    };

    enum PruningType{
        SEPARATESAMPLE = 0,
        BACKGROUNDREFERENCE = 1,
        COMBINEDREFERENCE = 2,
        COMBINEDSIGNAL = 3
    };

    enum WorkspaceType {
        FIT = 0,
        SIGNIFICANCE = 1,
        LIMIT = 2
    };

    enum class WorkspaceCreationType {
        FORFIT,
        FORPLOTS,
        BOTH
    };

    struct TemplateWeight{
        std::string function;
        std::string range;
        std::string name;
        double value;
    };

    enum ToysUsexRooFit {
        FALSE,
        TRUE,
        AUTO,
    };

    explicit TRExFit(const std::string& name="MyMeasurement");
    ~TRExFit();

    TRExFit(const TRExFit& t) = delete;
    TRExFit(TRExFit&& t) = delete;
    TRExFit& operator=(const TRExFit& t) = delete;
    TRExFit& operator=(TRExFit&& t) = delete;

    void SetPOI(const std::string& name="SigXsecOverSM", const std::string& unit="");
    void AddPOI(const std::string& name="SigXsecOverSM", const std::string& unit="");
    void SetStatErrorConfig(bool useIt=true, double thres=0.05, const std::string& cons="Poisson");
    void SetLumi(const double lumi);
    void SetFitType(FitType type);
    void SetLimitType( LimitType type );
    void SetSignificanceType( SignificanceType type );
    void SetFitRegion(FitRegion region);

    std::shared_ptr<Sample> NewSample(const std::string& name, Sample::SampleType type);
    std::shared_ptr<Systematic> NewSystematic(const std::string& name);
    Region* NewRegion(const std::string& name);
    void SkippedRegion(const std::string& name);

    // ntuple stuff
    void AddNtuplePath(const std::string& path);
    void AddFriendPath(const std::string& path);
    void SetMCweight(const std::string& weight);
    void SetSelection(const std::string& selection);
    void SetNtupleName(const std::string& name);
    void SetNtupleFile(const std::string& name);
    void ComputeBinning(int regIter);
    void DefineVariable(int regIter);

    // histogram stuff
    void AddHistoPath(const std::string& path);
    void SetHistoName(const std::string& name);

    void SmoothSystematics(std::string syst="all");

    // create new root file with all the histograms
    void CreateRootFiles();
    void WriteHistos(bool reWriteOrig=true) const;

    void DrawSystPlots() const;
    void DrawSystPlotsSumSamples() const;

    // read from ..
    void CloseInputFiles();
    void CreateEFTNominalSampleMaps();

    void CorrectHistograms();

    void DrawAndSaveAll(std::string opt="", const std::string& wsPath="");

    // separation plots
    void DrawAndSaveSeparationPlots() const;

    std::shared_ptr<TRExPlot> DrawSummary(std::string opt="", std::shared_ptr<TRExPlot> = nullptr, const std::string& wsPath="");
    void DrawMergedPlot(std::string opt="", const std::string& group="") const;
    void BuildYieldTable(const std::string& opt="", const std::string& group="", const std::string& wsPath="") const;

    // regions examples:
    // ...
    void DrawSignalRegionsPlot(int nCols=0,int nRows=0) const;
    void DrawSignalRegionsPlot(int nCols,int nRows, const std::vector<std::string>& regions) const;
    void DrawPieChartPlot(const std::string &opt="", int nCols=0,int nRows=0) const;
    void DrawPieChartPlot(const std::string &opt, int nCols,int nRows, const std::vector<std::string> &regions) const;

    void CreateCustomAsimov() const;

    /**
      * Runs code that replaces asimov data with custom asimov for unfolding
      */
    void UnfoldingAlternativeAsimov();
    void UnfoldingAlternativeAsimov(const Unfolding* unfolding);

    /**
     * @brief Creates the WS used for all purposes
     *
     * @param allRegions Consider all regions
     */
    void ToRooStats(const bool allRegions, const bool allBinsFitRegions) const;

    /**
     * @brief WS input for one channel
     *
     * @param meas measurement
     * @param ichan index of the channel
     * @param allRegions flag to tell if all regions are used
     * @return RooStats::HistFactory::Channel
     */
    RooStats::HistFactory::Channel OneChannelToRooStats(RooStats::HistFactory::Measurement* meas,
                                                        const int ichan,
                                                        const bool allRegions,
                                                        const bool allBinsFitRegions) const;

    /**
     * @brief WS input for one sample
     *
     * @param meas mesurement
     * @param h SampleHist for the sample
     * @param i_ch channel index
     * @param i_smp sample index
     * @param allRegions flag to tell if all regions are used
     * @return RooStats::HistFactory::Sample
     */
    RooStats::HistFactory::Sample OneSampleToRooStats(RooStats::HistFactory::Measurement* meas,
                                                      const SampleHist* h,
                                                      const int i_ch,
                                                      const int i_smp,
                                                      const bool allRegions,
                                                      const bool allBinsFitRegions) const;

    void SystPruning();
    void DrawPruningPlot() const;

    /**
      * A helper function to draw the plots with normalisation for each systematic
      */
    void DrawSystematicNormalisationSummary() const;

    // Extras for EFT
    void ProcessEFTInputs(bool overwrite=false);

    // fit etc...
    void Fit(bool isLHscanOnly);
    RooDataSet* DumpData( RooWorkspace *ws, const std::map <std::string, int>& regionDataType, const std::map <std::string, double>& npValues, const std::map <std::string, double>& poiValues);
    std::map < std::string, double > PerformFit( RooWorkspace *ws, RooDataSet* inputData, FitType fitType=SPLUSB, bool save=false);

    /**
     * @brief Get the reduced set of Regions from a WS
     *
     * @param regionsToFit
     * @return std::unique_ptr<RooWorkspace>
     */
    std::unique_ptr<RooWorkspace> PerformWorkspaceCombinationxRooFit(const std::vector <std::string>& regionsToFit) const;

    void PlotFittedNP();
    /**
     * @brief Plots the correlation matrix of the parameters from the fit
     *
     */
    void PlotCorrelationMatrix();

    /**
     * @brief Plots the unfolded data and the truth MC, will loop over the Unfolding objects
     *
     * @param frPath path to the fit results (withtout the file extension)
     * @param outputPath path to the output folder
     * @param wsFileName path to the workspace (with the extension)
     * @param wsName name of the WS in the ROOT file
     */
    void PlotUnfoldedData(const std::string& frPath,
                          const std::string& outputPath,
                          const std::string& wsFileName,
                          const std::string& wsName);

    /**
     * @brief Plots the unfolded data per unfolding object
     *
     * @param unfolding the unfolding object
     * @param frPath path to the fit results (withtout the file extension)
     * @param outputPath path to the output folder
     * @param wsFileName path to the workspace (with the extension)
     * @param wsName name of the WS in the ROOT file
     */
    void PlotUnfoldedData(const Unfolding* unfolding,
                          const std::string& frPath,
                          const std::string& wsFileName,
                          const std::string& outputPath,
                          const std::string& wsName);


    /**
     * @brief Plot covariance matrix of multiple unfoldings
     *
     * @param unfoldings vector of unfoldings
     * @param truthFileNames prefix for the names of the files where the histogram is read from (basically fName)
     * @param frPath path to the fit results
     * @param wsFileName path to the ROOT file with the WS
     * @param outputPath output pat
     * @param wsName name of the WS in the root file
     * @param nfs NormFactors
     */
    void PlotMultipleUnfoldingCovariance(const std::vector<const Unfolding*>& unfoldings,
                                         const std::vector<std::string>& truthFileNames,
                                         const std::string& frPath,
                                         const std::string& wsFileName,
                                         const std::string& outputPath,
                                         const std::string& wsName,
                                         const std::vector<std::shared_ptr<NormFactor> >& nfs);

    void GetLimit();
    void GetSignificance();

    /**
     * @brief Run the 1D LH scan
     *
     * @param ws workspace
     * @param varName name of the parameter to scan
     * @param data data
     * @param toysExtensions extra string extension that can be set for toys
     * @param plot flat to plot the curve
     */
    void GetLikelihoodScan(RooWorkspace *ws, const std::string& varName, RooDataSet* data, const std::string& toysExtensions = "", const bool plot = true) const;
    void Get2DLikelihoodScan(RooWorkspace *ws, const std::vector<std::string>& varName, RooDataSet* data) const;

    /**
     * @brief Read the fit results from the txt file
     *
     * @param fileName path tothe file
     * @return true file exists
     * @return false file does not exist
     */
    bool ReadFitResults(const std::string& fileName);

    void PrintConfigSummary() const;

    /**
     * @brief Returns the regions based on the the name
     *
     * @param name
     * @return Region*
     */
    Region* GetRegion(const std::string& name) const;
    std::shared_ptr<Sample> GetSample(const std::string& name) const;
    std::size_t GetSampleIndex(const std::string& name) const;

    /**
     * @brief Produce NP ranking file
     *
     * @param NPnames flag to tell which NPs to check
     */
    void ProduceNPRanking(const std::string& NPnames);

    /**
     * @brief Plot NP ranking
     *
     * @param flagSysts which systematics?
     * @param flagGammas what about gammas?
     */
    void PlotNPRanking(const bool flagSysts, const bool flagGammas);

    /**
     * @brief Decide whether the code is producing the ranking or plotting or both
     *
     */
    void PlotNPRankingManager();

    void PrintSystTables(std::string opt="", const std::string& wsPath="") const;

    void MergeSystematics(); // this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter
    void CombineSpecialSystematics(); // this will merge into single SystematicHist all the SystematicHist from systematics with same nuisance parameter

    // for template fitting
    void AddTemplateWeight(const std::string& name, double);

    std::vector<TemplateWeight> GetTemplateWeightVec(const TemplateInterpolationOption& opt);

    std::string GetWeightFunction(const std::vector<std::pair<double,std::string> >& templatePair, unsigned int itemp, const TemplateInterpolationOption& opt) const;

    /**
     * Function that returns string that represents smoothed abs value function
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSmoothLinearInterpolation(unsigned int itemp) const;

    /**
     * Helper function to calualte numerical correction to the smoothed linear function
     * @param parameter in the argument of the hyperbolic tangent function
     * @param size of the x axis interval
     * @param central position of the function
     * @param left position of the function on x axis
     * @param right position of the function on x axis
     * @param parameter of iteration, set to 1 for the first iteration
     * @return correction
     */
    double GetCorrection(double k, double width, double x_mean, double x_left, double init = 1.) const;

    /**
     * Helper function to approximate absolute value by sqrt(x^2+e)
     * @param index of the template
     * @return function in the string form
     */
    std::string GetSquareRootLinearInterpolation(unsigned int itemp) const;

    /**
     * Helper function to apply correction to square root aproximation
     * @param will return value for a from -a*sqrt(x^2+epsilon) +b
     * @param will return value for b from -a*sqrt(x^2+epsilon) +b
     * @param central position of the function
     * @param left position of the function on x axis
     * @param epsilon = precision of the approximation
     */
    void GetSquareCorrection(double *a, double *b, double x_i, double x_left, double epsilon) const;

    /**
     * Helper function to smooth morphing templates according to a given functional form (bin-by-bin)
     * @param name of the morphing parameter
     * @param functional form
     * @param pointer to an array of parameter values
     */
    void SmoothMorphTemplates(const std::string& name,const std::string& formula="pol1",double *p=0x0) const;

    /**
     * Helper function that draws plots with morphijng templates
     * @param name of the morphing parameter
     */
    void DrawMorphingPlots(const std::string& name) const;

    bool MorphIsAlreadyPresent(const std::string& name, const double value) const;

    // for grouped impact evaluation
    void ProduceSystSubCategoryMap();

    void BuildGroupedImpactTable();

    /**
     * Helper function that runs toys experiments
     */
    void RunToys();

    /**
     * Helper function to compute the variable string to be used when reading ntuples, for a given region, sample combination
     * @param pointer to the Region
     * @param pointer to the Sample
     */
    std::string Variable(Region *reg, const Sample *smp);

    /**
     * Helper function to compute the selection string to be used when reading ntuples, for a given region, sample combination
     * @param pointer to the Region
     * @param pointer to the Sample
     */
    std::string FullSelection(Region *reg, const Sample *smp);

    /**
     * Helper function to compute the weight string to be used when reading ntuples, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     */
    std::string FullWeight(const Region *reg, const Sample *smp, const Systematic *syst=nullptr,bool isUp=true);

    /**
     * Helper function to compute the full paths to be used when reading ntuples, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @param bool to specify if we are reading a sample that is used only for internal subtraction e.g. for JER
     */
    std::vector<std::string> FullNtuplePaths(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true, bool isSubtract=false, bool isFriend=false);

    /**
     * Helper function to compute the full paths to be used when reading histograms, for a given region, sample and systematic combination
     * @param pointer to the Region
     * @param pointer to the Sample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @param bool to specify if subtraction histograms need to be read
     */
    std::vector<std::string> FullHistogramPaths(Region *reg,Sample *smp,Systematic *syst=nullptr,bool isUp=true, const bool isFolded = false, const bool isSubtract=false);

    /**
     * A helper function to compute the fgull paths for a response matrix
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @param isSubtract to specify if the paths are the subtract paths (for JER)
     * @return the full path
     */
    std::vector<std::string> FullResponseMatrixPaths(const Region* reg,
                                                     const UnfoldingSample* smp,
                                                     const UnfoldingSystematic* syst = nullptr,
                                                     const bool isUp = true,
                                                     const bool isSubtract = false) const;

    /**
     * A helper function to compute the fgull paths for a migration matrix
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @param isSubtract to specify if the paths are the subtract paths (for JER)
     * @return the full path
     */
    std::vector<std::string> FullMigrationMatrixPaths(const Region* reg,
                                                      const UnfoldingSample* smp,
                                                      const UnfoldingSystematic* syst = nullptr,
                                                      const bool isUp = true,
                                                      const bool isSubtract = false) const;

    /**
     * A helper function to compute the fgull paths for acceptance
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @param isSubtract to specify if the paths are the subtract paths (for JER)
     * @return the full path
     */
    std::vector<std::string> FullAcceptancePaths(const Region* reg,
                                                 const UnfoldingSample* smp,
                                                 const UnfoldingSystematic* syst = nullptr,
                                                 const bool isUp = true,
                                                 const bool isSubtract = false) const;

    /**
     * A helper function to compute the fgull paths for selection efficiency
     * @param pointer to the Region
     * @param pointer to the UnfoldingSample
     * @param pointer to the Systematic (default = nullptr)
     * @param bool to specify up (true) or down (false) syst variation
     * @param isSubtract to specify if the paths are the subtract paths (for JER)
     * @return the full path
     */
    std::vector<std::string> FullSelectionEffPaths(const Region* reg,
                                                   const UnfoldingSample* smp,
                                                   const UnfoldingSystematic* syst = nullptr,
                                                   const bool isUp = true,
                                                   const bool isSubtract = false) const;

    /**
      * A helper function to combine paths for the truth distributions
      * @return a vector of the paths
      */
    std::vector<std::string> FullTruthPaths() const;

    /**
    * A helper function to get SampleHisto from a region that matches a name of the sample
    * @param Region
    * @@param name
    * @return SampleHist
    */
    std::shared_ptr<SampleHist> GetSampleHistFromName(const Region* const reg, const std::string& name) const;

    /**
     * A helper function to Copy a smoothing from a reference histogram to other histograms bin by bin
     * @param SampleHist
     * @param nominal histogram
     * @param up variation of histogram
     * @param down variation
     * @param flag if we want up or down variation returned
     * @return Up/Down variation for the new histogram
     */
    TH1* CopySmoothedHisto(const SampleHist* const sh, const TH1* const nominal, const TH1* const up, const TH1* const down, const bool isUp) const;

    /**
     * A helper function to get an index of systemaati variation that matches some name
     * @param SampleHist
     * @param name of the systematics
     * @return index
     */
    int GetSystIndex(const SampleHist* const sh, const std::string& name) const;

    std::shared_ptr<SystematicHist> CombineSpecialHistos(std::shared_ptr<SystematicHist> orig,
                                                         const std::vector<std::shared_ptr<SystematicHist> >& vec,
                                                         Systematic::COMBINATIONTYPE type,
                                                         const SampleHist* sh) const;

    /**
      *  A helper function to get the list of unique names of non-gamma systematics
      *  @return the list of unique non-gamma systematics
      */
    std::vector<std::string> GetUniqueSystNamesWithoutGamma() const;

    /**
      * A helper function to get the vector of non-validation regions
      * @return the vector on non-validation regions
      */
    std::vector<Region*> GetNonValidationRegions() const;

    /**
      * A helper function to get the vector of non-data, non-ghost samples
      * @return the vector of non-data, non-ghost samples
      */
    std::vector<std::shared_ptr<Sample> > GetNonDataNonGhostSamples() const;

    /**
      * A function that prepares signal inputs for unfolding.
      * Folded distributions are created
      */
    void PrepareUnfolding();

    /**
      * A function that prepares signal inputs for ONE unfolding object.
      * Folded distributions are created
      */
    void PrepareUnfolding(const Unfolding* unfolding, const int index);

    /**
      * A helper function to fold systematic distributions needed for unfolding
      * @param Unfolding
      * @param Folding manager
      * @param output file
      * @param Region
      * @param UnfoldingSample
      * @param Current UnfoldingSystystematics
      * @param Nominal migration matrix
      */
    void ProcessUnfoldingSystematics(const Unfolding* unfolding,
                                     FoldingManager* manager,
                                     TFile* file,
                                     const Region* reg,
                                     const UnfoldingSample* sample,
                                     const UnfoldingSystematic* syst,
                                     const TH2* nominal) const;

    /** A helper function that does the actual plotting of unfolded data
      * @param unfoded data
      * @param error band
      * @param error statonly error band
      * @param unfolding object
      * @param converter object
      */
    void PlotUnfold(TH1D* data,
                    TGraphAsymmErrors* band,
                    TGraphAsymmErrors* statOnly,
                    const Unfolding* unfolding,
                    const std::string& outputPath,
                    const YamlConverter& converter);

    /**
      * A helper function to plot migration or reposne matrix
      * @param matrix
      * @param flag if sample is migration
      * @param name of the region
      * @param name of the systematic
      * @param Unfolding object
      */
    void PlotMigrationResponse(const TH2* matrix,
                               const bool isMigration,
                               const std::string& regionName,
                               const std::string& systematicName,
                               const Unfolding* unfolding) const;

    /**
      * A helper function to force shape on some systematics
      */
    void RunForceShape();

    /**
     * A helper function to check if the fit/limit/significance is done on a mixture of real-data and Asimov-data
     */
    bool DoingMixedFitting() const;

    /**
     * A helper function to create a list of all regions to fit
     * @param Flag to consider FitRegion setting
     * @param If set, only regions with this DataType are included
     */
    std::vector<std::string> ListRegionsToFit(const bool useFitRegion, int dataType=-1) const;

    /**
     * A helper function to create a map with specified regions and corresponding DataType
     * @param List of regions to consider
     * @param Set to true to force to use ASIMOV everywhere
     */
    std::map<std::string, int> MapRegionDataTypes(const std::vector<std::string>& regionList,bool isBlind=false) const;

    /**
     * A helper function to read a fit result root file and store the NP values into a map
     * @param Text file containing fit results
     */
    std::map<std::string, double> NPValuesFromFitResultsFile(const std::string& fitResultsFile);

    /**
     * A helper function to replace the norm-factor expression of the last bin of a truth distribution when performing normalized-cross-section unfolding
     */
    void FixUnfoldingExpressions();
    void FixUnfoldingExpressions(const Unfolding* unfolding, TFile* input);
    void ProcessCorrelateUnfolding(const Unfolding* unfolding, TFile* input);
    void ProcessReplacementUnfolding(const Unfolding* unfolding, TFile* input);

    /**
     * @brief Function to run asymptotic or toy-based limits using xRooFit
     * @param ws the RooWorkspace
     * @param dataName the data to use for nll creation
     * @limitFile limitFile path to file where limit results are saved
     * @limitWSFile limitWSFile path to file where limit workspace is saved
     */
    void RunLimit(RooWorkspace* ws, const std::string& dataName, const std::string& limitFile, const std::string& limitWSfile) const;

    /**
     * @brief Function to run asymptotic or toy-based limit scan using xRooFit
     * @param hs the xRooHypoSpace to scan
     */
    void RunLimitScan(xRooNLLVar::xRooHypoSpace& hs) const;

    /**
      * @brief Function to run Limit estimation using toys
      * @param Data
      * @param Workspace
      */
    void RunLimitToys(RooAbsData* data, RooWorkspace* ws) const;

    /**
     * @brief Function to run asymptotic or toy-based significances using xRooFit
     * @param ws the RooWorkspace
     * @param dataName the data to use for nll creation
     * @significanceFile significanceFile path to file where significance results are saved
     * @significanceWSFile significanceWSFile path to file where significance workspace is saved
     */
    void RunSignificance(RooWorkspace* ws, const std::string& dataName, const std::string& significanceFile, const std::string& significanceWSfile) const;

    /**
     * @brief Function to run asymptotic or toy-based significance scan using xRooFit
     * @param hs the xRooHypoSpace to scan
     */
    void RunSignificanceScan(xRooNLLVar::xRooHypoSpace& hs) const;

    /**
     * @brief Run significance estimation from limits
     *
     * @param data
     * @param ws
     */
    void RunSignificanceToys(RooAbsData* data, RooWorkspace* ws) const;

    /**
      * Function to calculate chi^2 between unfolded data and a MC prediction
      * @param Errors
      * @param MC prediction
      * @param Correlation matrix
      * @param Unfolding object
      */
    double GetUnfoldedChi2(const TGraphAsymmErrors* error, const TH1* mc, const CorrelationMatrix* matrix, const Unfolding* unfolding) const;

    /**
      * Function to Get the Unfolding object corresponding to a given UnfoldingSample
      * @param UnfoldingSample
      */
    const Unfolding* GetCorrespondingUnfolding(const UnfoldingSample* sample) const;

    /**
      * Function to Get the Unfolding object corresponding to a given UnfoldingSystematic
      * @param UnfoldingSystematic
      */
    const Unfolding* GetCorrespondingUnfolding(const UnfoldingSystematic* syst) const;

    /**
     * @brief Plot the full covariance matrix and produce the yaml formats
     *
     * @param unfolding Pointer to the unfolding object
     * @param unfolded Pointer to the object containing the results
     * @param outputPath path to the output folder
     */
    void PlotCovarianceMatrix(const Unfolding* unfolding, const UnfoldingResult& unfolded, const std::string& outputPath) const;

    /**
     * @brief Method to run a fit to get the mixed data-MC workspace and data
     *
     * @param type type of the fit: FIT, LIMIT or SIGNIFICANCE
     * @return std::pair<std::unique_ptr<RooWorkspace>, std::unique_ptr<RooDataSet> >
     */
    std::pair<std::unique_ptr<RooWorkspace>, std::unique_ptr<RooDataSet> > PrepareMixedDataset(const WorkspaceType type);

    /**
     * @brief Produce the shape factor parametrisation and the validation plots for template morping
     *
     */
    void PrepareTemplateMorphing() const;

    /**
     * @brief Add the NFs that are used in the shape factor reparametrisation
     *
     */
    void ProcessShapeFactorReparametrization();

    /**
     * @brief Add the shape factor reparametrisation for every region based on the fitted values
     *
     */
    void ProcessTemplateMorphingShapeFactors();

    /**
     * @brief Add the shape factor reparametrisation for the EFT SM reference sanples in all regions
     *
     */
    void ProcessEFTShapeFactors();

    /**
     * @brief Function to calculate the global correlation parameter per unfolding
     * < sqrt( 1 - (V_ii V_ii^-1)^-1) >
     *
     * @return vector<double>
     */
    std::vector<double> CalculateGlobalCorrelation();

    /**
     * @brief Read the binning from the txt files
     *
     */
    void ReadRegionBinning();

    /**
     * @brief Get the Raw Region Vector object
     *
     * @return std::vector<Region*>
     */
    std::vector<Region*> GetRawRegionVector() const;

    /**
     * @brief Run the non-profiled fit
     *
     */
    void NonProfiledFit();

    /**
     * @brief Store the cross-unfolding correlations - this is useful e.g. for EFTFitter
     *
     * @param cor Covariance matrix
     * @param binNames Names of the bins for each observable
     * @param unfoldingNames Names of the individual measurements
     * @param folder Output folder path
     */
    void SplitUnfoldingCorrelationMatrix(const std::vector<std::vector<double> >& cor,
                                         const std::vector<std::vector<std::string> >& binNames,
                                         const std::vector<std::string>& unfoldingNames,
                                         const std::string& folder) const;


    /**
     * @brief Check if validation regions are set
     *
     * @return true
     * @return false
     */
    bool HasValidationRegions() const;

    /**
     * @brief Check if the region got a drop bins
     *
     * @param reg
     * @return true
     * @return false
     */
    bool HasDropBinRegions() const;

    /**
     * @brief Add metadata to the WS
     *
     * @param ws the workspace
     * @param allRegions include validation regions?
     */
    void AddWSMetadata(RooWorkspace* ws, const bool allRegions) const;

    /**
     * @brief Add titles for NPs and NormFactors to the WS metadata
     *
     * @param node
     */
    void AddParameterTitlesToMetadata(const xRooNode& node) const;

     /**
      * @brief Plot ranking from the covariance matrix breakdown
      *
      * @param fitTool Fitting tool
      */
    void PlotRankingFromCovMatrix(const FittingTool& fitTool) const;

     /**
      * @brief set the blinding properties in all regions
      *
      * @param wsPath
      * @param hasPostfit
      */
    void SetBlindingInRegions(const std::string& wsPath, const bool hasPostfit);

    /**
     * @brief Take the WS and store it as HS3 format
     *
     * @param wsFilePath
     * @param wsPath
     * @param outputPath
     */
    void TranslateWSToHS3(const std::string& wsFilePath,
                          const std::string& wsPath,
                          const std::string& outputPath) const;

    /**
     * @brief Read Stat only result from the error decomposition txt file
     *
     * @param poi
     * @param path
     */
    void ReadStatOnlyFromErrorDecomposition(const std::string& poi, const std::string& path);

    /**
     * @brief Read Stat only results from the error decomposition for all POIs
     *
     * @param pois
     * @param path
     */
    void ReadStatOnlyFromErrorDecomposition(const std::vector<std::string>& pois, const std::string& path);

    /**
     * @brief Get the Stat Only Error From Decomposition object
     *
     * @param poi
     * @return double
     */
    double GetStatOnlyErrorFromDecomposition(const std::string& poi) const;

    /**
     * @brief Get the Stat Only Error From Decomposition Up object
     *
     * @param poi
     * @return double
     */
    double GetStatOnlyErrorFromDecompositionUp(const std::string& poi) const;

    /**
     * @brief Get the Stat Only Error From Decomposition Down object
     *
     * @param poi
     * @return double
     */
    double GetStatOnlyErrorFromDecompositionDown(const std::string& poi) const;

    /**
     * @brief Plot the uncertainty breakdown
     *
     * @param fitTool
     * @param categories
     */
    void PlotUnfoldingErrors(const FittingTool& fitTool, const std::vector<std::string>& categories);

    // -------------------------

    std::string fName;
    bool fUseFriend;
    std::string fDir;
    std::string fLabel;
    std::string fInputFolder;
    std::string fInputName;
    std::string fBinningsPath;

    std::vector < TFile* > fFiles;

    std::vector < std::unique_ptr<Region> > fRegions;
    std::vector < std::string > fSkippedRegions;
    std::vector < std::shared_ptr<Sample> > fSamples;
    std::vector < std::shared_ptr<Systematic> > fSystematics;
    std::vector < std::shared_ptr<NormFactor> >fNormFactors;
    std::vector < std::shared_ptr<ShapeFactor> > fShapeFactors;
    std::vector < std::string > fSystematicNames;
    std::vector < std::string > fNormFactorNames;
    std::vector < std::string > fShapeFactorNames;

    std::map<std::string, std::map< std::string, std::vector<std::shared_ptr<SampleHist> > > > fEFTNominalSampleMaps;

    std::vector<std::string> fPOIs;
    std::map<std::string,std::string> fPOIunit;
    std::string fPOIforLimit;
    std::string fPOIforSig;
    bool fUseStatErr;
    double fStatErrThres;
    std::string fStatErrCons;
    bool fUseGammaPulls;

    double fLumi;
    double fLumiScale;

    double fThresholdSystPruning_Normalisation;
    double fThresholdSystPruning_Shape;
    double fThresholdSystLarge;
    std::vector<std::string> fNtuplePaths;
    std::vector<std::string> fNtupleFiles;
    std::vector<std::string> fNtupleNames;
    std::vector<std::string> fFriendPaths;
    std::vector<std::string> fFriendFiles;
    std::vector<std::string> fFriendNames;
    std::string fMCweight;
    std::string fSelection;

    std::vector<std::string> fResponseMatrixNames;
    std::vector<std::string> fResponseMatrixFiles;
    std::vector<std::string> fResponseMatrixPaths;
    std::vector<std::string> fResponseMatrixNamesNominal;
    std::vector<std::string> fAcceptanceNames;
    std::vector<std::string> fAcceptanceFiles;
    std::vector<std::string> fAcceptancePaths;
    std::vector<std::string> fAcceptanceNamesNominal;
    std::vector<std::string> fSelectionEffNames;
    std::vector<std::string> fSelectionEffFiles;
    std::vector<std::string> fSelectionEffPaths;
    std::vector<std::string> fSelectionEffNamesNominal;
    std::vector<std::string> fMigrationNames;
    std::vector<std::string> fMigrationFiles;
    std::vector<std::string> fMigrationPaths;
    std::vector<std::string> fMigrationNamesNominal;

    std::vector<std::string> fHistoPaths;
    std::vector<std::string> fHistoFiles;
    std::vector<std::string> fHistoNames;
    std::vector<std::string> fHistoNamesNominal;

    std::unique_ptr<FitResults> fFitResults;

    bool fWithPullTables;

    int fIntCode;

    int fInputType; // 0: histo, 1: ntup

    bool fSystDataPlot_upFrame;
    bool fStatOnly;
    bool fGammasInStatOnly;
    bool fStatOnlyFit;
    bool fFixNPforStatOnlyFit;

    std::vector<std::string> fRegionsToPlot;
    std::vector<std::string> fSummaryPlotRegions;
    std::vector<std::string> fSummaryPlotLabels;
    std::vector<std::string> fSummaryPlotValidationRegions;
    std::vector<std::string> fSummaryPlotValidationLabels;

    double fYmin;
    double fYmax;
    double fRatioYmin;
    double fRatioYmax;
    double fRatioYminPostFit;
    double fRatioYmaxPostFit;
    std::string fRatioYtitle;
    TRExPlot::RATIOTYPE fRatioType;

    std::string fLumiLabel;
    std::string fCmeLabel;

    std::string fSuffix;
    std::string fSaveSuffix;

    bool fUpdate;
    bool fKeepPruning;

    double fBlindingThreshold;
    Common::BlindingType fBlindingType;

    int fRankingMaxNP;
    std::string fRankingOnly;
    std::string fRankingPlot;
    std::string fImageFormat;
    std::string fPlotLabel;

    bool fDoSummaryPlot;
    bool fDoMergedPlot;
    bool fDoTables;
    bool fDoSignalRegionsPlot;
    bool fDoPieChartPlot;

    std::string fGroupedImpactCategory;

    std::string fSummaryPrefix;

    //
    // Fit characteristics
    //
    FitType fFitType;
    FitRegion fFitRegion;
    std::map< std::string, double > fFitNPValues;
    std::map< std::string, double > fFitFixedNPs;
    std::string fFitNPValuesFromFitResultsFile;
    bool fInjectGlobalObservables;
    std::string fPOIAsimovOverride;
    std::string fCustomAsimovOverride;
    std::map< std::string, double > fFitPOIAsimov;
    bool fFitIsBlind;
    bool fUseRnd;
    double fRndRange;
    long int fRndSeed;
    std::vector<std::string> fVarNameLH;
    std::vector<std::vector<std::string> > fVarName2DLH;
    bool fCreateCache;
    bool fUseCache;
    double fLHscanMin;
    double fLHscanMax;
    int fLHscanSteps;
    int fLHscanStep;
    double fLHscanMinY;
    double fLHscanMaxY;
    int fLHscanStepsY;
    int fLHscanStepY;
    std::vector<std::string> fVarNameMinos;
    std::vector<std::string> fVarNameHide;
    std::string fWorkspaceFileName;
    bool fDoGroupedSystImpactTable;
    std::map<std::string, std::string> fSubCategoryImpactMap;

    //
    // Limit parameters
    //
    LimitType fLimitType;
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
    SignificanceType fSignificanceType;
    bool fSignificanceIsBlind;
    bool fSignificanceDoInjection;
    double fSignificancePOIAsimov;
    std::string fSignificanceParamName;
    double fSignificanceParamValue;
    std::string fSignificanceOutputPrefixName;
    bool fSignificanceUseAutoDiff;

    bool fCleanTables;
    bool fSystCategoryTables;

    std::vector< std::string > fRegionGroups;

    bool fKeepPrefitBlindedBins;
    std::unique_ptr<TH1D> fBlindedBins;

    std::string fCustomAsimov;

    /// flag to control if custom asimov dataset should be written into a Workspace file
    bool fWriteCustomAsimovToWS;

    std::string fTableOptions;

    int fGetChi2;

    HistoTools::SmoothOption fSmoothOption;

    bool fSuppressNegativeBinWarnings;

    std::vector<std::string> fAddAliases;

    std::vector<std::string> fCustomFunctions;
    std::vector<std::string> fCustomIncludePaths;
    std::vector<std::string> fCustomFunctionsExecutes;

    std::vector<std::string> fMorphParams;
    std::vector<std::pair<double,std::string> > fTemplatePair;
    std::vector<TRExFit::TemplateWeight> fTemplateWeightVec;
    TemplateInterpolationOption fTemplateInterpolationOption;

    std::string fBootstrap;
    std::string fBootstrapSyst;
    std::string fBootstrapSample;
    std::string fBootstrapNomHistos;
    int fBootstrapIdx;

    std::vector<std::string> fDecorrSysts;
    std::string fDecorrSuff;

    bool fDoNonProfileFit;
    double fNonProfileFitSystThreshold;
    int fFitToys;
    int fToysHistoNbins;
    std::string fSmoothMorphingTemplates;
    int fPOIPrecision;

    std::vector<std::string> fRankingPOIName;
    std::vector<int> fRankingUpperAxisNdivision;
    std::vector<double> fRankingPOIAxisScale;
    bool fUsePDGRounding;
    bool fUsePDGRoundingTxt;
    bool fUsePDGRoundingTex;
    bool fPropagateSystsForMorphing;
    PruningType fPruningType;

    std::vector<int> fPrePostFitCanvasSize;
    std::vector<int> fSummaryCanvasSize;
    std::vector<int> fMergeCanvasSize;
    std::vector<int> fPieChartCanvasSize;
    std::vector<int> fNPRankingCanvasSize;

    std::vector<std::string> fBlindedParameters;

    double fLabelX;
    double fLabelY;
    double fLegendX1;
    double fLegendX2;
    double fLegendY;

    double fLabelXSummary;
    double fLabelYSummary;
    double fLegendX1Summary;
    double fLegendX2Summary;
    double fLegendYSummary;

    double fLabelXMerge;
    double fLabelYMerge;
    double fLegendX1Merge;
    double fLegendX2Merge;
    double fLegendYMerge;

    int fLegendNColumns;
    int fLegendNColumnsSummary;
    int fLegendNColumnsMerge;

    bool fShowRatioPad;
    bool fShowRatioPadSummary;
    bool fShowRatioPadMerge;

    std::string fExcludeFromMorphing;

    bool fDoSystNormalizationPlots;

    int fDebugNev;

    int fCPU;

    std::vector< std::string > fSeparationPlot;

    std::vector<std::unique_ptr<Unfolding> > fUnfolding;

    std::vector<std::unique_ptr<UnfoldingSample> > fUnfoldingSamples;
    std::vector<std::unique_ptr<UnfoldingSystematic> > fUnfoldingSystematics;
    bool fHasAcceptance;
    std::vector<std::unique_ptr<TruthSample> > fTruthSamples;
    PruningUtil::SHAPEOPTION fPruningShapeOption;
    bool fSummaryLogY;
    /// This variable is needed only for multifit
    bool fUseInFit;
    bool fUseInComparison;
    bool fReorderNPs;
    bool fBlindSRs;
    bool fHEPDataFormat;
    bool fUheppFormat;
    bool fAlternativeShapeHistFactory;
    int fFitStrategy;
    bool fBinnedLikelihood;
    bool fRemoveLargeSyst;
    bool fRemoveSystOnEmptySample;
    bool fValidationPruning;
    bool fUsePOISinRanking;
    bool fUseHesseBeforeMigrad;
    bool fUseNllInLHscan;
    int fLimitToysStepsSplusB;
    int fLimitToysStepsB;
    int fLimitToysScanSteps;
    double fLimitToysScanMin;
    double fLimitToysScanMax;
    int fToysSeed;
    int fLimitToysSeed;
    bool fLimitPlot;
    bool fLimitFile;
    std::string fLimitToysSuffix;
    int fLimitFitStrategy;
    ToysUsexRooFit fLimitToysUsexRooFit;

    int fSignificanceToysStepsSplusB;
    int fSignificanceToysStepsB;
    int fSignificanceToysSeed;
    bool fSignificancePlot;
    int fSignificanceFitStrategy;
    bool fSignificanceToysUsexRooFit;

    bool fDataWeighted;
    int fRegularizationType;
    bool fSpeedUpFit;
    double fNPCutOff;
    bool fDoExtendedCovariances;
    bool fUnfoldingShowStat;
    bool fUnfoldingShowUncertaintyBreakdown;
    bool fUnfoldingUncertaintyBreakdownTotal;
    bool fCombinerFormat;
    bool fToysStatOutput;
    std::string fToysNpValuesFile;
    bool fUseRebinned;
    bool fApplyGammaCorrection;
    double fErrorSigma;
    bool fShapeFactorReparametrisation;
    bool fHasTemplateMorphing;
    TemplateMorpher::FIT_DIMENSIONALITY fTemplateMorphingFitDimensionality;
    Morphing fMorphingSetting;
    std::string fFitResultsRootFile;
    bool fRecreateBinningFiles;
    std::vector<std::string> fBinningFilesToCreate;
    bool fUseHesse;
    WorkspaceCreationType fWorkspaceCreationType;
    bool fHasValidationRegions;
    bool fHasDropBinRegions;
    bool fErrorDecomposition;
    int fMaximumNumberFCNcalls;
    double fToleranceScale;
    bool fShiftGlobalObservablesInRanking;
    bool fUseAutoDiff;
    bool fNoPrePostFitLabel;
    std::map<std::string, std::pair<double, double> > fToyRandomPOIStartingValues;
    std::map<std::string, std::pair<double, double> > fToysLHScanCondition;
    bool fToysLHScanForAll;
    bool fToysOnlyStatFluctuation;
    std::string fNLLOffset;

    std::map<std::string, double> fStatOnlyErrorDecomposition;
    std::map<std::string, double> fStatOnlyErrorDecompositionUp;
    std::map<std::string, double> fStatOnlyErrorDecompositionDown;

    EFTConfig fEFTConfig;
    bool fProducePerRegionWS;
    bool fDrawPruningPlot;
    std::string fPreFitLabel;
    std::string fPostFitLabel;
    int fToysForErrorBand;
    bool fHasDataSetInConfig;
};

#endif
