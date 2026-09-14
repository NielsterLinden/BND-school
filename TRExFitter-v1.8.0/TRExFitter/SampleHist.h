#ifndef SAMPLEHIST_H
#define SAMPLEHIST_H

/// Framework includes
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/PruningUtil.h"

/// ROOT includes
#include "Rtypes.h"

/// c++ includes
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class TFile;
class TH1;
class TPad;
class Sample;
class NormFactor;
class ShapeFactor;
class SystematicHist;

namespace Common {
    struct FolderStructure;
}

class SampleHist {
public:
    explicit SampleHist();
    explicit SampleHist(Sample *sample,TH1 *hist);

    /**
     * @brief Construct a new Sample Hist object
     *
     * @param sample Pointer to the sample
     * @param histoName name of the histogram
     * @param fileName path of the ROOT file
     * @param structure if se to nullptr, will not use any folder structure
     */
    explicit SampleHist(Sample *sample,
                        const std::string& histoName,
                        const std::string& fileName,
                        const bool& dropBinHistoNeeded,
                        const bool& checkHistogram,
                        const Common::FolderStructure* structure = nullptr
                        );
    ~SampleHist() = default;
    SampleHist(SampleHist&& s) = delete;
    SampleHist& operator=(SampleHist&& s) = delete;

    SampleHist(const SampleHist& s) = delete;  // copy constructor
    SampleHist& operator=(const SampleHist& s) = delete;

    TH1* GetHist() const;
    const Sample* GetSample() const {return fSample;}
    std::shared_ptr<SystematicHist> AddOverallSyst(const std::string& name,const std::string& storedName,double up,double down);
    std::shared_ptr<SystematicHist> AddStatSyst(const std::string& name,const std::string& storedName,int i_bin);
    std::shared_ptr<SystematicHist> AddHistoSyst(const std::string& name,const std::string& storedName,TH1* h_up,TH1* h_down);

    /**
     * @brief Adds Histo systematic histogram
     *
     * @param name name of the histogram
     * @param storedName name used for storing
     * @param histoName_up name of the up variation histogram
     * @param fileName_up path to the up ROOT file
     * @param histoName_down name of the down variation
     * @param fileName_down path to the down ROOT file
     * @param pruned flag for pruning
     * @param structure if set to nullptr, no folder structure will be used
     * @return std::shared_ptr<SystematicHist>
     */
    std::shared_ptr<SystematicHist> AddHistoSyst(const std::string& name,
                                                 const std::string& storedName,
                                                 const std::string& histoName_up,
                                                 const std::string& fileName_up,
                                                 const std::string& histoName_down,
                                                 const std::string& fileName_down,
                                                 const bool& dropBinHistoNeeded=false,
                                                 const bool& checkHistogram=false,
                                                 int pruned=0,
                                                 const Common::FolderStructure* structure = nullptr);
    std::shared_ptr<SystematicHist> GetSystematic(const std::string& systName) const;
    std::shared_ptr<SystematicHist> GetSystFromNP(const std::string& NuisParName) const;
    std::shared_ptr<NormFactor> AddNormFactor(const std::string& name,double nominal, double min, double max);
    std::shared_ptr<NormFactor> AddNormFactor(std::shared_ptr<NormFactor> normFactor);
    std::shared_ptr<NormFactor> GetNormFactor(const std::string& name) const;
    std::shared_ptr<ShapeFactor> AddShapeFactor(const std::string& name,double nominal, double min, double max);
    std::shared_ptr<ShapeFactor> AddShapeFactor(std::shared_ptr<ShapeFactor> shapeFactor);
    std::shared_ptr<ShapeFactor> GetShapeFactor(const std::string& name) const;

    bool HasSyst(const std::string& name) const;
    bool HasNorm(const std::string& name) const;
    bool HasShapeFactor(const std::string& name) const;

    /**
     * @brief White the histogram to a file
     *
     * @param blindedBins list of blinded bins
     * @param scales scaling of the histograms
     * @param gammaThreshold threshold used for MC stat pruning
     * @param pruningReference reference for pruning - bin yields
     * @param regName region name
     * @param f file
     * @param reWriteOrig flag to set if the original histo should be rewritten
     * @param histErrors Updated errors for histograms if correlate gammas for samples is used
     */
    void WriteToFile(const std::vector<int>& blindedBins,
                     const std::vector<double>& scales,
                     const double gammaThreshold,
                     const std::vector<double>& pruningReference,
                     const std::string& regName,
                     const bool useRegularBinning,
                     TFile* f=nullptr,
                     bool reWriteOrig=true,
                     const std::vector<double>& histErrors = {});
    void ReadFromFile();

    void FixEmptyBins(const bool suppress);
    void NegativeTotalYieldWarning(TH1* hist, double yield) const;

    void Print() const;

    void Rebin(int ngroup = 2, const Double_t* xbins = 0);
    void DrawSystPlot(TH1* h_data,
                      bool SumAndData,
                      bool bothPanels,
                      const bool isUnfolding,
                      const std::vector<int>& blindedBins={} ) const;
    void SmoothSyst(const HistoTools::SmoothOption &opt, const bool useAlternativeShapeHistFactory, const std::string& syst="all", bool force=false);

    void Divide(  SampleHist* sh);
    void Multiply(SampleHist* sh);
    void Add(     SampleHist* sh,double scale=1.);
    void Scale(double scale);

    void SampleHistAdd(SampleHist* h, double scale = 1.);
    void SampleHistAddNominal(SampleHist* h, double scale);
    void CloneSampleHist(SampleHist* h, const std::set<std::string>& names, double scale = 1.);

    /**
     * @brief Function to apply the systematic pruning
     *
     * @param pu PruningUtil pointer
     * @param hTot total histogram - needed for some pruning types
     */
    void SystPruning(PruningUtil *pu, const TH1* hTot=nullptr);

    void DrawSystPlotUpper(TPad* pad0,
                           TH1* nominal,
                           TH1* nominal_orig,
                           TH1* syst_up,
                           TH1* syst_up_orig,
                           TH1* syst_down,
                           TH1* syst_down_orig,
                           TH1* data,
                           TH1* tmp,
                           bool SumAndData,
                           bool bothPanels) const;

    void DrawSystPlotRatio(TPad* pad1,
                           TH1* nominal,
                           TH1* nominal_orig,
                           TH1* syst_up,
                           TH1* syst_up_orig,
                           TH1* syst_down,
                           TH1* syst_down_orig,
                           TH1* data,
                           TH1* tmp,
                           bool SumAndData,
                           const std::vector<int>& blindedBins) const;

    std::vector<double> GetDataScales() const;

    static std::vector<std::vector<double> > GetGammaCorrelationUncertainties(const std::vector<const TH1*>& histos);

    std::string fName;
    Sample *fSample;
    std::unique_ptr<TH1> fHist;
    std::unique_ptr<TH1> fHist_orig;
    std::unique_ptr<TH1> fHist_regBin;
    std::unique_ptr<TH1> fHist_preSmooth; // new - to use only for syst plots
    std::shared_ptr<TH1> fHist_postFit;
    std::string fFileName;
    std::string fHistoName;
    bool fIsData;
    bool fIsSig;
    std::map<std::string,bool> fIsMorph;

    std::vector < std::shared_ptr<SystematicHist> > fSyst;

    std::vector <  std::string > fNormFactorNames;
    std::vector <  std::string > fShapeFactorNames;

    // other useful info
    std::string fFitName;
    std::string fRegionName;
    std::string fRegionLabel;
    std::string fVariableTitle;
    bool fSystSmoothed;

    std::string fPlotLabel;
};

#endif

