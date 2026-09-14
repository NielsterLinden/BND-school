#ifndef RANKINGMANAGER_H
#define RANKINGMANAGER_H

#include "TRExFitter/YamlConverter.h"

#include "TMatrixD.h"

#include <string>
#include <map>
#include <memory>
#include <vector>

class CorrelationMatrix;
class FitResults;
class FittingTool;
class RooDataSet;
class RooWorkspace;
class RooSimultaneous;
class NormFactor;
class Region;
class xRooNode;

namespace RooStats {
    class ModelConfig;
}

class RankingManager {
public:

    struct RankingValues {
        double central;
        double up;
        double down;
        double upPrefit;
        double downPrefit;
    };

    explicit RankingManager();
    ~RankingManager() = default;
    RankingManager(const RankingManager& r) = delete;
    RankingManager(RankingManager&& r) = delete;
    RankingManager& operator=(const RankingManager& r) = delete;
    RankingManager& operator=(RankingManager&& r) = delete;

    inline std::size_t GetNPsize() const {return fNuisPars.size();}

    inline void SetOutputPath(const std::string& path){fOutputPath = path;}
    inline void SetInjectGlobalObservables(const bool flag){fInjectGlobalObservables = flag;}
    inline void SetNPValues(const std::map<std::string, double>& m){fFitValues = m;}
    inline void SetFixedNPs(const std::map<std::string, double>& m){fFitFixedNPs = m;}
    inline void SetFitStrategy(const int s){fFitStrategy = s;}
    inline void SetPOINames(const std::vector<std::string>& name){fPOINames = name;}
    inline void SetNCPU(const int n){fCPU = n;}
    inline void SetStatOnly(const bool flag){fStatOnly = flag;}
    inline void SetPlotLabel(const std::string& l){fPlotLabel = l;}
    inline void SetLumiLabel(const std::string& l){fLumiLabel = l;}
    inline void SetCmeLabel(const std::string& l){fCmeLabel = l;}
    inline void SetUseHEPDataFormat(const bool flag){fHEPDataFormat = flag;}
    inline void SetName(const std::string& l){fName = l;}
    inline void SetSuffix(const std::string& l){fSuffix = l;}
    inline void SetMaxNPPlot(const std::size_t n){fRankingMaxNP = n;}
    inline void SetRankingPOIName(const std::string& l){fRankingPOIName = l;}
    inline void SetRankingPOIAxisScale(const double scale){fRankingPOIAxisScale = scale;}
    inline void SetRankingCanvasSize(const std::vector<int>& s){fNPRankingCanvasSize = s;}
    inline void SetUsePOISinRanking(const bool flag){fUsePOISinRanking = flag;}
    inline void SetUseHesseBeforeMigrad(const bool flag){fUseHesseBeforeMigrad = flag;}
    inline void SetUpperAxisNdivisions(const int value){fSetNdivisions = value;}
    inline void SetRegularizationType(const int value){fRegularizationType = value;}
    inline void SetRandomNPs(const double value, const bool randomize, const long int seed) {
        fRandomNP = value;
        fRandomize = randomize;
        fRandomSeed = seed;
    }

    inline void SetToleranceScale(const double value) {fToleranceScale = value;}

    void AddNuisPar(const std::string& name, const bool isNF);

    void RunRanking(const std::shared_ptr<xRooNode>& pdf,
                    const std::shared_ptr<xRooNode>& pdfPrefit,
                    RooWorkspace* ws,
                    RooDataSet* data,
                    const std::vector<std::shared_ptr<NormFactor> >& nfs) const;

    /**
     * @brief Plot the NP ranking
     *
     * @param reg list of regions
     * @param nfs names of the NFs
     * @param sfs names of the SFs
     * @param flagSysts show systs?
     * @param flagGammas show gammas?
     */
    void PlotRanking(const std::vector<Region* >& reg,
                     const std::vector<std::string>& nfs,
                     const std::vector<std::string>& sfs,
                     const bool flagSysts,
                     const bool flagGammas);

    /**
     * @brief Plot the NP ranking
     *
     * @param reg list of regions
     * @param nfs names of the NFs
     * @param sfs names of the SFs
     * @param flagSysts show systs?
     * @param flagGammas show gammas?
     */
    void PlotRanking(const std::vector<std::unique_ptr<Region> >& reg,
                     const std::vector<std::string>& nfs,
                     const std::vector<std::string>& sfs,
                     const bool flagSysts,
                     const bool flagGammas);

    /**
     * @brief Dump individual syst + stat covariance matrices
     *
     * @param pois
     * @param nfs
     * @param folder
     * @param fr
     * @param statOnly
     * @param isUnfolding
     */
    void DumpCovariances(const std::vector<std::pair<std::string,double> >& pois,
                         const std::vector<std::string>& nfs,
                         const std::string& folder,
                         const FitResults* fr,
                         const std::map<std::string, double>& statOnly,
                         const bool isUnfolding) const;

    std::vector<std::map<std::string, YamlConverter::RankingContainer> > GetRankingImpacts(const std::string& folder,
                                                                                           const std::vector<std::pair<std::string, double> >& pois,
                                                                                           const std::vector<std::string>& nfs,
                                                                                           std::vector<std::string>& params) const;

    /**
     * @brief Set the Maximum FCN calls
     *
     * @param max
     */
    inline void SetMaximumFCNcalls(const int max) {fMaximumNumberOfFCNcalls = max;}

    /**
     * @brief Set the Shift Global Observables
     *
     * @param flag
     */
    inline void SetShiftGlobalObservables(const bool flag) {fShiftGlobalObservables = flag;}

    /**
     * @brief Set the Is Covariance Breakdown
     *
     * @param flag
     */
    inline void SetIsCovarianceBreakdown(const bool flag) {fIsCovarianceBreakdown = flag;}

    /**
     * @brief Set the Ranking Container object
     *
     * @param container
     */
    inline void SetRankingContainer(const std::vector<YamlConverter::RankingContainer>& container) {fRankingContainer = container;}

    /**
     * @brief Read results from the ranking txt file and set the RankingContainer
     *
     * @param flagSysts
     * @param flagGammas
     *
     */
    void ReadRankingResults(const bool flagSysts, const bool flagGammas);

    /**
     * @brief Produce grouped impact from the results
     *
     * @param categoryMap
     */
    void ProduceGroupedImpact(const std::map<std::string, std::string>& categoryMap) const;

    /**
     * @brief Produce grouped impact for a single POI
     *
     * @param file output file
     * @param categoryMap caregory map
     */
    void ProcessSingleGroupedImpact(std::ofstream* file, const std::map<std::string, std::string>& categoryMap) const;

private:

    std::string fOutputPath;
    std::vector<std::pair<std::string,bool> > fNuisPars;
    bool fInjectGlobalObservables;
    std::map<std::string, double> fFitValues;
    std::map<std::string, double> fFitFixedNPs;
    int fFitStrategy;
    int fRegularizationType;

    std::vector<std::string> fPOINames;
    int fCPU;
    bool fStatOnly;
    std::string fPlotLabel;
    std::string fLumiLabel;
    std::string fCmeLabel;
    bool fHEPDataFormat;
    std::string fName;
    std::string fSuffix;
    std::size_t fRankingMaxNP;
    std::string fRankingPOIName;
    double fRankingPOIAxisScale;
    std::vector<int> fNPRankingCanvasSize;
    bool fUsePOISinRanking;
    bool fUseHesseBeforeMigrad;
    int fSetNdivisions;
    double fRandomNP;
    bool fRandomize;
    long int fRandomSeed;
    int fMaximumNumberOfFCNcalls;
    double fToleranceScale;
    bool fShiftGlobalObservables;
    bool fIsCovarianceBreakdown;
    std::vector<YamlConverter::RankingContainer> fRankingContainer;

    std::vector<double> RunSingleFit(FittingTool* fitTool,
                                     RooWorkspace* ws,
                                     RooStats::ModelConfig *mc,
                                     RooSimultaneous *simPdf,
                                     RooDataSet* data,
                                     const std::pair<std::string, bool>& np,
                                     const bool isUp,
                                     const bool isPrefit,
                                     const RankingValues& values,
                                     const std::vector<double>& muhat) const;

    void PlotTotalCovariance(const TMatrixD& matrix,
                             const std::vector<std::string>& names,
                             const bool isPrefit) const;

    double Symmetrize(const std::string& param, const double up, const double down) const;

    /**
     * @brief Get impact on all POIs from ranking by taking "NP uncertainty * correlation * total uncertainty"
     *
     * @param corr Correlation matrix
     * @param totalUncertainty total uncertainty per POI
     * @param np nuisance parameter
     * @param values NP central value, uncertainties
     * @param isPrefit flag
     * @return std::vector<double> impact on each POI
     */
    std::vector<double> ImpactFromCorrelations(const CorrelationMatrix* corr,
                                               const std::vector<std::pair<std::string, double> >& totalUncertainty,
                                               const std::pair<std::string, bool>& np,
                                               const RankingValues& values,
                                               const bool isPrefit) const;
};

#endif
