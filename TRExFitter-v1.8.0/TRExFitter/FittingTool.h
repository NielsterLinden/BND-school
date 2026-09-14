//
// Fit.h
// ==========
// Allows to implement all necessary informations and functions to perform likelihood fits.
//

#ifndef FITTINGTOOL_H
#define FITTINGTOOL_H

#include "TRExFitter/YamlConverter.h"

/// RooStats include
#include "RooStats/ModelConfig.h"

// c++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

/// Forward declaration
class RooFitResult;
class TString;
class RooAbsPdf;
class RooAbsData;

class FittingTool {

public:

    //
    // Standard C++ functions
    //
    explicit FittingTool();
    ~FittingTool() = default;
    FittingTool(const FittingTool& rhs) = delete;
    FittingTool& operator=(const FittingTool& rhs) = delete;
    FittingTool(FittingTool&& rhs) = delete;
    FittingTool& operator=(FittingTool&& rhs) = delete;

    //
    // Gettters and setters
    //
    inline void SetNCPU (const int cpu){ m_CPU = cpu; }

    /**
     * @brief Get the Fit Result object
     *
     * @return const RooFitResult*
     */
    inline const RooFitResult* GetFitResult() const {return m_fitResult.get();}

    void AddValPOI(const std::string& name, const double value);
    void ReplacePOIVal(const std::string& name, const double value);

    inline void ConstPOI(const bool constant) { m_constPOI = constant; }

    inline void NoGammas()      { m_noGammas=true;      }
    inline void NoSystematics() { m_noSystematics=true; }

    inline void SetRandomNP(const double rndNP,
                            const bool rndize,
                            const long int rndSeed=-999) { m_randomNP = rndNP; m_randomize = rndize; m_randSeed = rndSeed; }

    void SetSubCategories();
    void SetSystMap(const std::map<std::string, std::string>& subCategoryMap) { m_subCategoryMap = subCategoryMap; SetSubCategories(); } // fills both m_subCategoryMap and m_subCategories

    inline void ResetFixedNP() { m_constNP.clear(); m_constNPvalue.clear(); };
    inline void FixNP(const std::string& np, const double value) { m_constNP.emplace_back(np); m_constNPvalue.emplace_back(value); }
    inline void FixNPs(const std::vector<std::string>& nps, const std::vector<double>& values) { m_constNP = nps; m_constNPvalue = values; }
    inline void SetNPs(const std::vector<std::string>& nps, const std::vector<double>& values) { m_initialNP = nps; m_initialNPvalue = values; }

    inline void UseMinos(const std::vector<std::string>& minosvar){ m_useMinos = true; m_varMinos = minosvar; }

    inline void SetExternalConstraints(const RooArgSet* externalConstraints = 0){ m_externalConstraints = externalConstraints; }

    inline void SetStrategy(const int strategy){m_strategy = strategy;}

    inline void SetUseHesse(const bool flag){m_useHesse = flag;}

    inline void SetUseHesseBeforeMigrad(const bool flag){m_hesseBeforeMigrad = flag;}

    inline void SetSpeedUpFit(const bool flag){m_speedUpFit = flag;}

    inline void SetSmallNPs(const std::vector<std::pair<std::string,double> >& v) {m_smallNPs = v;}

    inline void SetErrorSigma(const float sigma) {m_errorSigma = sigma;}

    inline void SetConstantNFs(const std::vector<std::string>& constNFs) {m_constNFs = constNFs;}

    inline void SetToleranceScale(const double value) {m_toleranceScale = value;};

    //
    // Specific functions
    //
    double FitPDF(RooStats::ModelConfig* model,
                  RooAbsPdf* fitpdf,
                  RooAbsData* fitdata,
                  bool fastFit = false,
                  bool noFit = false);

    void SaveFitResult(const std::string& fileName);

    void ExportFitResultInTextFile(const std::string& finalName,
                                   const std::vector<std::string>& blinded);

    std::map < std::string, double > ExportFitResultInMap();

    void GetGroupedImpact(RooStats::ModelConfig* model,
                          RooAbsPdf* fitpdf,
                          RooAbsData* fitdata,
                          RooWorkspace* ws,
                          const std::string& categoryOfInterest,
                          const std::string& outFileName,
                          const std::string& fileName,
                          const std::string& lumiLabel,
                          const std::string& cmeLabel,
                          const bool useHEPData ) const;

    void FitExcludingGroup(bool excludeGammas,
                           bool statOnly,
                           RooAbsData*& fitdata,
                           RooAbsPdf*& fitpdf,
                           RooStats::ModelConfig* mc,
                           RooWorkspace* ws,
                           const std::string& category,
                           const std::vector<std::string>& affectedParams) const;

    /**
     * @brief Get the Number of all parameters (NPs + POIs)
     *
     * @param mc Model config
     * @return std::size_t
     */
    static std::size_t GetNumberOfNPsPlusPOIs(RooStats::ModelConfig* mc);

    /**
     * @brief Set the Max FCN calls
     *
     * @param max
     */
    inline void SetMaxFCNcalls(const int max) {m_maxFCNcalls = max;}

    /**
     * @brief Calculate uncertainty decomposition from covariance matrix, for the given POI, and print-out the splitting
     *
     * @param mc Model config
     * @param poiName Name of the POI to run the decomposition for
     * @return std::map<std::string,double>
     */
    std::map<std::string,double> CalculateErrorDecomposition(RooStats::ModelConfig* mc, const std::string& poiName);

    /**
     * @brief Export error decomposition from covariance matrix
     *
     * @param poiName Name of the POI to export the decomposition for
     * @param fileName Output file name
     */
    void ExportErrorDecompositionInTextFile(const std::string& poiName, const std::string& fileName);

    /**
     * @brief Export error decomposition per category from covariance matrix
     *
     * @param poiName Name of the POI to export the decomposition for
     * @param fileName Output file name
     * @param lumiLabel Lumi label
     * @param cmeLabel CME label
     * @param hepDataFolder HEPData folder
     * @param storeHEPData store HEPData format
     */
    void ExportErrorDecompositionGroupInTextFile(const std::string& poiName,
                                                 const std::string& fileName,
                                                 const std::string& lumiLabel,
                                                 const std::string& cmeLabel,
                                                 const std::string& hepDataFolder,
                                                 const bool storeHEPData);

    /**
     * @brief Set the AD option
     *
     * @param flag
     */
    void SetUseAutoDiff(const bool flag);

    /**
     * @brief Set the likelihood offset method
     *
     * @param flag (valid: "none", "initial", "bin")
     */
    void SetNLLOffset(const std::string& val);

    /**
     * @brief Get the Prefit uncertainty for a parameter
     *
     * @param fr Fit results
     * @param param Parameter name
     * @param index -1 for down uncertainty, 0 for summetrised, +1 for up uncertainty
     * @return double
     */
    double GetPrefitUncertainty(const RooFitResult* fr, const std::string& param, int index) const;

    /**
     * @brief Get the Prefit Value
     *
     * @param fr
     * @param param
     * @return double
     */
    double GetPrefitValue(const RooFitResult* fr, const std::string& param) const;

    /**
     * @brief Get the Ranking Container object
     *
     * @param poiName
     * @param exclude
     * @param flagSysts
     * @param flagGammas
     * @return std::vector<YamlConverter::RankingContainer>
     */
    std::vector<YamlConverter::RankingContainer> GetRankingContainer(const std::string& poiName,
                                                                     const std::vector<std::string>& exclude,
                                                                     const bool flagSysts,
                                                                     const bool flagGammas) const;

    /**
     * @brief Get the Err Decomp Map object
     *
     * @return const std::map<std::string, std::map<std::string, double> >&
     */
    const std::map<std::string, std::map<std::string, double> >& GetErrDecompCategoryMap() const {return m_categoryDecompositionMap;}

private:


    void CheckUnderconstraint(const RooRealVar* const var) const;

    void PrintMinuitHelp() const;

    void SpeedUpFit(RooAbsReal* nll, RooStats::ModelConfig* mc) const;

    int m_CPU;
    std::vector<std::pair<std::string, double> > m_valPOIs;
    bool m_useMinos;
    std::vector<std::string> m_varMinos;
    bool m_constPOI;
    std::unique_ptr<RooFitResult> m_fitResult;
    bool m_noGammas;
    bool m_noSystematics;
    bool m_randomize;
    double m_randomNP;
    long int m_randSeed;
    double m_minNll;

    std::vector<std::string> m_constNP;
    std::vector<double> m_constNPvalue;
    std::vector<std::string> m_initialNP;
    std::vector<double> m_initialNPvalue;
    std::vector<std::string> m_constNFs;

    std::map<std::string, std::string> m_subCategoryMap;
    std::set<std::string> m_subCategories;

    const RooArgSet* m_externalConstraints;
    int m_strategy;
    bool m_useHesse;
    bool m_hesseBeforeMigrad;
    bool m_speedUpFit;
    std::vector<std::pair<std::string,double> > m_smallNPs;
    double m_errorSigma;
    int m_maxFCNcalls;
    double m_toleranceScale;
    bool m_useAutoDiff;
    std::string m_nllOffset;

    std::map<std::string,std::map<std::string,double> > m_errDecompMap;
    std::map<std::string,std::map<std::string,double> > m_errDecompMapUp;
    std::map<std::string,std::map<std::string,double> > m_errDecompMapDown;
    std::map<std::string, std::map<std::string, double > > m_categoryDecompositionMap;
};


#endif //FitTools
