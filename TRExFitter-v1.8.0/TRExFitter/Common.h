#ifndef COMMON_H
#define COMMON_H

/// TRExFitter stuff
#include "TRExFitter/Logger.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/NormFactor.h"

#include "UnfoldingCode/UnfoldingCode/UnfoldingResult.h"

// ROOT stuff
#include "TF1.h"
#include "TChain.h"
#include "TCanvas.h"

#include "xRooFit/xRooNode.h"

/// c++ stuff
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class FitResults;
class Region;
class RooFitResult;
class ShapeFactor;
class TDirectory;
class TFile;
class TGraphAsymmErrors;
class TH1;
class TH2;
class TH1D;
class TPad;
class Unfolding;

namespace TRExFitter{
    extern int DEBUGLEVEL;
    void SetDebugLevel(int level=0);
    void SetLoggingFormat(std::string param);
    extern bool SHOWYIELDS; // flag to show or not yields in plots
    extern bool SHOWSTACKSIG;  // flag to show signal or not
    extern bool ADDSTACKSIG;  // flag to add signal to total or not
    extern bool SHOWNORMSIG;  // flag to show normalized signal or not
    extern bool SHOWOVERLAYSIG;  // flag to show overlayed signal or not
    extern double SHOWOVERLAYSIG_CUSTOMSCALE; //if SHOWOVERLAYSIG is set, signal can be scaled by custom constant scale (post-fit POI scale factor will not be applied)
    extern bool SHOWCHI2;
    extern bool SHOWSTACKSIG_SUMMARY;  // flag to show signal or not in Summary Plot
    extern bool SHOWNORMSIG_SUMMARY;  // flag to show normalized signal or not in Summary Plot
    extern bool SHOWOVERLAYSIG_SUMMARY;  // flag to show overlayed signal or not in Summary Plot
    extern bool LEGENDLEFT;  // flag to show sample names on left aligned in the legend
    extern bool LEGENDRIGHT;  // flag to show sample names on right aligned in the legend
    extern bool PREFITONPOSTFIT;  // flag to show prefit background as dashed line on postfit plots
    extern bool SYSTCONTROLPLOTS;
    extern bool SYSTDATAPLOT;
    extern bool SYSTERRORBARS;
    extern bool SPLITHISTOFILES;
    extern bool HISTOCHECKCRASH;
    extern bool REMOVEXERRORS;
    extern bool OPRATIO;
    extern bool NORATIO; // flag to hide ratio pad
    extern double CORRELATIONTHRESHOLD;
    extern bool MERGEUNDEROVERFLOW;
    extern std::map< std::string,std::string > SYSTMAP;
    extern std::map< std::string,std::string > SYSTTEX;
    extern std::map< std::string,std::string > NPMAP;
    extern std::vector< std::string > IMAGEFORMAT;
    //
    extern std::map< std::string, double > OPTION;
    extern std::map<std::string, std::unique_ptr<TFile> > TFILEMAP;
    extern bool GUESSMCSTATERROR;
    extern bool CORRECTNORMFORNEGATIVEINTEGRAL;
    extern bool USEFOLDERSTRUCTURE;
    extern bool USEHEPDATAROUNDING;
    extern bool KEEPDATAERRORS;
    extern std::string EXPERIMENT_LABEL;
}

namespace Common {

enum BlindingType {
    SOVERB = 1,
    SOVERSPLUSB = 2,
    SOVERSQRTB = 3,
    SOVERSQRTSPLUSB = 4
};

/**
 * @brief A "container" that holds the names of the Region, Sample and Systematics that will be used to
 * build the internal directory structure in a ROOT file
 *
 */
struct FolderStructure {
    std::string region;
    std::string sample;
    std::string systematics;
};

/**
    * @brief Combine a vector of histograms into 1
    * @param vec std::vector<std::shared_ptr<TH1> >
    * @return std::unique_ptr<TH1>
**/

std::unique_ptr<TH1D> CombineHistosFromHistosVec(const std::vector<std::shared_ptr<TH1D> >& vec);

/**
    * @brief Combine two vectors with the same type into one
    *
    * @param v1 std::vector<T>
    * @param v2 std::vector<T>
    * @return std::vector<T>
**/

template <class T>
std::vector<T> CombineVectors(std::vector<T>& a, std::vector<T>& b, bool unique = false) {
    std::vector<T> c;
    if (unique) {
        std::set<T> s;
        for (const auto& i : a) s.insert(i);
        for (const auto& i : b) s.insert(i);
        for (const auto& i : s) c.push_back(i);
    } else {
        c = a;
        c.insert(c.end(), b.begin(), b.end());
    }
    return c;
}

/** Convert a string to a number, error out if there's any issue
 *
 * @brief Safely convert a string to a number
 * @param toConvert the string to convert
 * @return the converted number
**/
template <typename T>
T convertStoNum(const std::string& toConvert){
    T converted;
    std::string::size_type pos;
    try{
        using U = std::decay_t<T>;
        if constexpr(std::is_same_v<U, float>) {
            converted = std::stof(toConvert, &pos);
        } else  if constexpr(std::is_same_v<U, double>) {
            converted = std::stod(toConvert, &pos);
        } else  if constexpr(std::is_same_v<U, int>) {
            converted = std::stoi(toConvert, &pos);
        } else {
            static_assert(false, "Invalid type passed to convertStoNum");
        }
        if(pos != toConvert.size()){
            LOG(ERROR) << "Partially converted string: " << toConvert << "\n";
            exit(EXIT_FAILURE);
        }
    }
    catch(const std::exception& err){
        LOG(ERROR) << "Exception caught: " << toConvert << " " << err.what() << "\n";
        exit(EXIT_FAILURE);
    }
    return converted;
}


/**
 * @brief Get the ROOT file from a path. If the file is already open, return the pointer. If it is not open, open it and add it to the map of open files
 *
 * @param fileName path
 * @return TFile*
 */
void Glob(const std::string& pattern, const std::string& tree, TChain& t);
TFile* GetFile(const std::string& fileName);
TH1D* HistFromNtuple(const std::string& ntuple, const std::string& frnd, const std::string& variable, int nbin, double xmin, double xmax, const std::string& selection, const std::string& weight, const std::vector<std::string>& aliases={}, int Nev=-1);
TH1D* HistFromNtupleBinArr(const std::string& ntuple, const std::string& frnd, const std::string& variable, int nbin, double *bins, const std::string& selection, const std::string& weight, const std::vector<std::string>& aliases={}, int Nev=-1);

/**
 * @brief Get a histogram from a file
 *
 * @param fullName name of the file and the name of the histogram as one string
 * @param structure folder structure
 * @return std::unique_ptr<TH1>
 */
std::unique_ptr<TH1> HistFromFile(const std::string& fullName, const bool& dropBinHistoNeeded = false, const bool& checkHistogram = false, const FolderStructure* structure = nullptr );

/**
 * @brief Gets Histogram from a file and the histo path
 *
 * @param fileName path of the ROOT file
 * @param histoName name of the histogram
 * @param structure structure for the folder path, if nullptr, no folder structure is used
 * @return std::unique_ptr<TH1>
 */
std::unique_ptr<TH1> HistFromFile(const std::string& fileName, const std::string& histoName, const bool& dropBinHistoNeeded = false, const bool& checkHistogram = false, const FolderStructure* structure = nullptr);

/**
 * @brief Get 2D histogram from a file
 *
 * @param fullName name of the file and the name of the histogram as one string
 * @return std::unique_ptr<TH2>
 */
std::unique_ptr<TH2> Hist2DFromFile(const std::string& fullName);

/**
 * @brief Get 2D histogram from a filepath and a histo path
 *
 * @param fileName path of the ROOT file
 * @param histoName histo path
 * @return std::unique_ptr<TH2>
 */
std::unique_ptr<TH2> Hist2DFromFile(const std::string& fileName, const std::string& histoName);

/**
 * @brief Write histogram to a file that will be opened
 *
 * @param h The histogram
 * @param fileName Path to the file
 * @param structure the folder structure
 * @param option TFile::Open parameter
 */
void WriteHistToFile(TH1* h, const std::string& fileName, const FolderStructure& structure, const std::string& option = "UPDATE");

/**
 * @brief Write histogram to a file
 *
 * @param h the histogram
 * @param f the file
 * @param structure the folder structure
 */
void WriteHistToFile(TH1* h, TFile* f, const FolderStructure& structure);


void MergeUnderOverFlow(TH1* h);
std::vector<std::string> CreatePathsList(std::vector<std::string> paths, std::vector<std::string> pathSufs,
                                         std::vector<std::string> files, std::vector<std::string> fileSufs,
                                         std::vector<std::string> names, std::vector<std::string> nameSufs);
std::vector<std::string> CombinePathSufs(std::vector<std::string> pathSufs, std::vector<std::string> newPathSufs, const bool isFolded = false);
std::vector<std::string> ToVec(const std::string& s);
std::string ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace);
std::vector< std::pair < std::string,std::vector<double> > > processString(std::string target);

bool StringsMatch(const std::string& s1, const std::string& s2);
int wildcmp(const char *wild, const char *string);

int FindInStringVector(const std::vector<std::string>& v, const std::string& s);
int FindInStringVectorOfVectors(const std::vector<std::vector<std::string> >& v, const std::string& s, const std::string& ss);
double GetSeparation( TH1D* S1, TH1D* B1 );

/**
  * Function to blind data and retrieve the blinding histogram
  * @param data histogram
  * @param indices ob blinded bins
  * @return Histogram with non-zero bins on positions to be blinded
  */
std::unique_ptr<TH1D> BlindDataHisto(TH1* h_data, const std::vector<int>& blindedBins);

bool SmoothHistogram( TH1* h, double nsigma=2. ); // forceFlat: 0 force no flat, 1 force flat, -1 keep it free
void SmoothHistogramTtres( TH1* h);

double CorrectIntegral(TH1* h, double *err=0);

TH1D* MergeHistograms(const std::vector<std::unique_ptr<TH1> >& hVec,bool fixLastBinWidth=false);
TH1D* MergeHistograms(const std::vector<TH1*>& hVec,bool fixLastBinWidth=false);

/**
  * A function to apply PDG rounding rules to values
  * @param A reference to mean value
  * @param A reference to uncertainty
  */
int ApplyPDGrounding(double& mean, double& error);

/**
  * A helper function to round error according to PDG rules
  * @param The value of error that will be rounded
  * @return number of iterations of multiplication/division by 10 needed to reach the same precision for nominal value
  */
int ApplyErrorRounding(double& error, int& sig);

/**
  * A helper function to round value to n decimal places
  * @param A value that needs to be rounded
  * @param Number of multiplications/divisions by 10 needed to get the value that can be rounded
  */
void RoundToSig(double& value, const int& n);

std::string KeepSignificantDigits(double value, const int n);

TH1* CloneNoError(TH1* h,const char* name="");

unsigned int NCharactersInString(const std::string& s,const char c);

bool CheckExpression(const std::string& s);

/**
 * Helper function to parsee the string to identify if the chosen option needs to run the fit
 * @return true if needs to run the fit
 */
bool OptionRunsFit(const std::string& opt);

/**
 * Helper function to make a copy of histogram with no errors in bins
 * This is useful when doing some scaling operations like Add/Divide
 * without modifying the original uncertainty in bins
 * @param histogram to be copied
 * @return histogram with no errors
 */
std::unique_ptr<TH1> GetHistCopyNoError(const TH1* const hist);

void ScaleMCstatInHist(TH1* hist, const double scale);

/// BW added functions to help pad bin numbers in gamma NP plot

std::vector<std::string> mysplit(const std::string & s, char delimiter);
std::string addpad( const std::string & input, const char filler, const unsigned width );
std::string pad_trail( const std::string & input );

// Helper functions to drop norm or shape part from systematic variations
void DropNorm(TH1* hUp,TH1* hDown,TH1* hNom);
void DropNorm(TH1* hUp,TH1* hDown, const double intNom);
void DropShape(TH1* hUp,TH1* hDown,TH1* hNom);

void SetHistoBinsFromOtherHist(TH1* toSet, const TH1* other);

/**
 * Helper function to get the integral of a histogram only considering positive bins
 * @param pinter to histogram
 * @return the effective integral
 */
double EffIntegral(const TH1* const h);

/**
  * A helper function that gets the indices of the blinded bins in a region
  * @param reg the given Region
  * @param node xRooNode for the region
  * @param type blinding type
  * @param threshold blinding threshold
  * @return indices of the bins (ROOT index convention)
  */
std::vector<int> GetBlindedBins(const Region* reg,
                                std::shared_ptr<xRooNode> node,
                                const BlindingType type,
                                const double threshold);

/**
  * A helper function to retrieve the blinded bins from histograms
  * @param signal histogram
  * @param background histogram
  * @param blinding type
  * @param blinding threshold
  * @return blinded bins
  */
std::vector<int> ComputeBlindedBins(const TH1* signal,
                                    const TH1* bkg,
                                    const BlindingType type,
                                    const double threshold);

/**
  * A helper function to combine histograms from a vector of full paths
  * @param A vector where each element represents the full path
  * @return a combined histogram
  */
std::unique_ptr<TH1> CombineHistosFromFullPaths(const std::vector<std::string>& paths);

/**
  * A helper function to combine 2D histograms from a vector of full paths
  * @param A vector where each element represents the full path
  * @return a combined 2D histogram
  */
std::unique_ptr<TH2> CombineHistos2DFromFullPaths(const std::vector<std::string>& paths);

/**
  * A helper function to calculate error band on ratio
  * @param Graph with total uncertainties
  * @param data
  * @return ratio graph
  */
std::unique_ptr<TGraphAsymmErrors> GetRatioBand(const TGraphAsymmErrors* total, const TH1D* data);

/**
 * A helper function to divide bin content by bin width
 * @param histogram
 */
void ScaleByBinWidth(TH1* h);

/**
 * A helper function to scale tgraph by constant value
 * @param histogram
 * @param the value
 */
void ScaleByConst(TGraphAsymmErrors* g, const double scale);

/**
 * A helper function to divide bin content by bin width
 * @param graph
 */
void ScaleByBinWidth(TGraphAsymmErrors* g);

/**
  * A helper function to transform an integer into a string, with leading zeros
  * @param input int
  * @param number of places (default = 3)
  * @return string
  */
std::string IntToFixLenStr(int i,int n=3);

/**
  * Case insensitive string to boolean
  * @param string
  * @return conversion result
  */
bool StringToBoolean(std::string param);

/**
  * Get paths to files containing a key from a folder
  * @param folder path
  * @param key
  * @param key2
  * @return vector of paths
  */
std::vector<std::string> GetFilesMatchingString(const std::string& folder, const std::string& key, const std::string& key2);

/**
  * Merge txt files into another one
  * @param paths to the input files
  * @param path to the output file
  * @return outputfile path
  */
void MergeTxTFiles(const std::vector<std::string>& input, const std::string& out);

/**
  * A helper function to check the validity of a string
  * @param
  * return name
  */
std::string CheckName(const std::string& name);

/**
 * @brief Check the name of the TTree to make sure it does not contain any "weird" characters
 *
 * @param name
 * @return std::string
 */
void CheckTreeName(const std::string& name);

/**
  * A helper function to remove quotes
  * @param
  * return name
  */
std::string RemoveQuotes(const std::string& name);

/**
  * A helper function to remove spaces
  * @param
  * return name
  */
std::string RemoveSpaces(const std::string& name);

/**
  * A helper function to remove comments
  * @param
  * return name
  */
std::string RemoveComments(const std::string& s);

/**
 * @brief A helper function to split a string
 *
 * @param s original string
 * @param c splitter
 * @param removeQuotes remove quotes?
 * @param removeComments remove comments?
 * @return std::vector<std::string>
 */
std::vector<std::string> Vectorize(const std::string& s,char c,bool removeQuotes=true,bool removeComments=true);

/**
  * A helper function to calculate powers of number with integer exponent
  * @param value
  * @param exponent
  * return value^exponent
  */
double IntPow(const double value, const int exp);

/**
 * @brief Recalculate the impact of the shape factor reparametrisation
 *
 * @param sf ShapeFactor
 * @param nfs NormFactors
 * @param ibin Bin index
 * @param isPostFit Flag to tell if we should use postfir results
 * @param fitResult Stores the postfit results
 * @return std::vector<double> nominal,min,max
 */
std::vector<double> CalculateShapeFactorReparametrization(const std::shared_ptr<ShapeFactor>& reg,
                                                          const std::vector<std::shared_ptr<NormFactor> >& nfs,
                                                          const int ibin,
                                                          const bool isPostFit,
                                                          const FitResults* fitResult);

/**
 * @brief Set mean values of the graph to the histo
 *
 * @param graph the graph to be updated
 * @param hist the histogram that the graph will be set to
 */
void SetGraphToHist(TGraphAsymmErrors* graph, TH1* hist);

/**
 * @brief Set the graph mean values on y axis to zero
 *
 * @param graph
 */
void SetGraphToZero(TGraphAsymmErrors* graph);

/**
 * @brief Prepare the folder structure in a ROOT file if ti does not exist
 *
 * @param file ROOT file
 * @param structure the folder structure
 * @return TDirectory* pointer to the final directory
 */
TDirectory* SetFolderStructure(TFile* file,
                               const FolderStructure& structure);

/**
 * @brief Propagate uncertainties to the Unfolding Result object
 *
 * @param result
 * @param unfolding
 * @param truth
 * @param fitResults
 * @param wsFileName
 * @param wsName
 * @param frFileName
 * @param normFactors
 * @param statErrorMap
 * @param statOnly
 */
void GetUnfoldingResult(UnfoldingResult& result,
                        const Unfolding* unfolding,
                        const TH1* truth,
                        const FitResults* fitResults,
                        const std::string& wsFileName,
                        const std::string& wsName,
                        const std::string& frFileName,
                        const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                        const std::map<std::string, double>& statErrorMap,
                        const bool statOnly);

std::string ReplaceAllCharacters(std::string string, const char from, const char to);

/**
 * @brief Replaces the string representing a path to a histogram with a string
 * that has a new folder (before "/")
 *
 * @param name string to be replaced
 * @param newFolder new folder name
 * @return std::string
 */
std::string ReplaceFolderName(std::string name, const std::string& newFolder);

/**
 * @brief The same thing as separator.join(strings_to_join) in python
 */
std::string join_strings(const std::vector<std::string> &strings_to_join, const std::string& separator);

/**
 * @brief The same thing as input_string.split(separator) in python
 */
std::vector<std::string> SplitString(std::string input_string, const std::string& separator);

/**
 * @brief Strip characters in chars_to_remove from beginning and end of the input_string
 */
void StripString(std::string *input_string, const std::string &chars_to_remove = " \n\t\r");

/**
 * @brief Split the string based on separator and then strip white chars
 */
std::vector<std::string> SplitAndStripString(const std::string &input_string, const std::string& separator);

/**
 * @brief Return true if main_string starts with prefix
 */
bool StartsWith(const std::string &main_string, const std::string &prefix);

/**
 * @brief Get a histogram from a node
 *
 * @param node
 * @return std::shared_ptr<TH1>
 */
std::shared_ptr<TH1> HistoFromxRooNode(const std::shared_ptr<xRooNode>& node, const TH1* hist);

/**
 * @brief Get a histogram from a node
 *
 * @param node
 * @return std::shared_ptr<TH1>
 */
std::shared_ptr<TH1> HistoFromxRooNode(const xRooNode& node, const TH1* hist);

/**
 * @brief Get the error band
 *
 * @param node
 * @param hist
 * @param fr
 * @param nToys
 * @return std::shared_ptr<TH1>
 */
std::unique_ptr<TGraphAsymmErrors> ErrorFromxRooNode(const xRooNode& node,
                                                     const TH1* hist,
                                                     const FitResults* fr,
                                                     const int nToys);

/**
 * @brief Get the xRooNode from a ROOT file, if Fit Results path is provided also attach it
 * xRooNode keeps the file open as long the node is not out of scope
 *
 * @param wsPath Path to the workspace
 * @param frPath Path to the Fit Results file, empty string if prefit
 * @return xRooNode
 */
xRooNode ReadWSandFitResults(const std::string& wsPath, const std::string& frPath);

/**
 * @brief Get xRooNode for a sample in a given region
 *
 * @param node
 * @param region
 * @param sample
 * @return std::shared_ptr<xRooNode>
 */
std::shared_ptr<xRooNode> XRooNodeSampleFromNode(std::shared_ptr<xRooNode> node, const std::string& region, const std::string& sample);

/**
 * @brief Get the list of non-signal samples
 *
 * @param samples
 * @param region name of the region
 * @return std::string
 */
std::string NonSignalSampleList(const std::vector<std::shared_ptr<SampleHist> >& samples, const std::string& region);

/**
 * @brief Get the Normalisation Component for systematic uncertainty
 *
 * @param syst
 * @param nominal
 * @return double
 */
double GetNormalisationComponent(const TH1* syst, const TH1* nominal);

/**
 * @brief Get the global correlation coefficient from RooFit result
 *
 * @param fr Results form the fit
 * @param pois THe list of the POIs
 * @return double
 */
double GCC(const RooFitResult* fr, const std::vector<std::string>& pois);

/**
 * @brief Print the results of the limits
 *
 * @param hs Hypospace
 * @param injectedHs Hypospace for the injected limit
 * @param name name for some printing
 * @param isBlind fit is blind?
 * @param outputFilePath path to the ROOT output file
 * @param paramName name of the parameter
 * @param paramValue the value of the parameter to be stored in the output file
 */
void ProcessLimitOutput(const xRooNLLVar::xRooHypoSpace& hs,
                        const xRooNLLVar::xRooHypoSpace* injectedHs,
                        const std::string& name,
                        const bool isBlind,
                        const std::string& outputFilePath,
                        const std::string& paramName,
                        float paramValue,
                        bool asymptotic = true);

/**
 * @brief Print the results of the significance calculation
 *
 * @param significanceObserved significance and uncertainty
 * @param significanceExpected significance and uncertainty
 * @param significanceInjects significance and uncertainty for the injected significance
 * @param name mame for some printing
 * @param isBlind significance is blind?
 * @param outputFilePath path to the ROOT file with the output file
 *
 * @return true/false
 */
bool ProcessSignificanceOutput(const std::pair<double, double>& significanceObserved,
                               const std::pair<double, double>& significanceExpected,
                               const std::pair<double, double>& significanceInjected,
                               const std::string& name,
                               const bool isBlind,
                               const std::string& outputFilePath);
/**
 * @brief Save canvas in file at given path
 *
 * @param c canvas to save
 * param path path to file to save canvas in
 */
void SaveCanvasAs(const TCanvas& c, const std::string& path);

/**
 * @brief Get the list of elements matching a string with a wildcard
 *
 * @param vector
 * @param match
 * @return std::vector<std::string>
 */
std::vector<std::string> MatchingElememnts(const std::vector<std::string>& vector, const std::string& match);

/**
 * @brief Get colour from RGB
 *
 * @param colours
 * @return int
 */
int ColorFromRGB(const std::vector<std::string>& colours);

/**
 * @brief Read RooFitResult from a ROOT file
 *
 * @param file
 * @return std::unique_ptr<RooFitResult>
 */
std::unique_ptr<RooFitResult> RooFitResultFromFile(TFile* file);


/**
 * @brief Format the first argument to xRooNLLVar::xRooHypoSpace::scan()
 *
 * @param stepsSplusB Number of S+B toys
 * @param stepsB Number of B-only toys
 * @return the formatted string, or an empty string in case of errors
 */
TString xRooFitToysScanArg(int stepsSplusB, int stepsB);

/**
 * @brief Remove objects from cleanup to speed up the destructors
 *
 * @param pad
 */
void CleanTPad(TPad* pad);

/**
 * @brief Is the given POI an unfolding POI?
 *
 * @param poi
 * @param unfoldings
 * @return true
 * @return false
 */
bool IsUnfoldingPOI(const std::string& poi, const std::vector<std::unique_ptr<Unfolding> >& unfoldings);

}

#endif
