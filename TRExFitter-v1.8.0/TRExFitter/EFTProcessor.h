#ifndef EFTPROCESSOR_H
#define EFTPROCESSOR_H

/// Framework includes
#include "TRExFitter/EFTConfig.h"
#include "TRExFitter/TRExFit.h"

/// ROOT includes
#include "Rtypes.h"

/// c++ includes
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

/// Forward class declaration
class EFTConfig;
class NormFactor;
class Region;

class EFTProcessor {
public:

    explicit EFTProcessor();
    explicit EFTProcessor(const std::string& param, const std::string& title);
    ~EFTProcessor() = default;

    class EFTVariation {
    public:
        std::string fParam;
        std::string fTitle;
        std::string fName;
        std::string fValueString;
        std::vector<double> fValues;
    };

    void SetEFTOrder(const EFTConfig::EFTOrder& order);

    void Print() const;

    void PrintCoeffMap() const;

    std::string GetParametrisation() const;

    std::pair<std::string, std::vector<int> > GetParametrisationAndTerms(const int op_idx_match = -1) const;

    void PrintValuesMap() const;

    std::string MakeTitle(const std::string& instr) const;

    void Resize();

    void DrawEFTInputs(const std::vector<std::unique_ptr<Region> >& Regions,
                       const EFTConfig& eftConfig) const;

    void FitEFTInputs(const std::vector<std::unique_ptr<Region> >& Regions,
                      const std::string& folder,
                      const std::string& sampleName,
                      const EFTConfig& eftConfig);

    void ApplyMuFactExpressions(const std::vector<std::unique_ptr<Region> >& regions,
    	                        std::vector<std::shared_ptr<NormFactor> >& normFactors) const;

    void ReadEFTFitResults(const std::vector<std::unique_ptr<Region> >& negions, const std::string& fileName);

    std::string ModifyEFTParamsOrder(const std::string& name) const;



    std::string fName;
    std::string fTitle;
    std::map<std::string,std::vector<std::string> > fRefMap;
    std::map<std::string,std::vector<double> > fExtrapMap;

    /**
     * @brief List of operators
     *
     */
    std::vector<std::string> fOperators;
    std::vector<std::string> fParams;
    std::vector<std::string> fTitles;
    std::map<std::string,std::string> fTitleMap;

    /**
     * @brief Stores which powers each coefficient has for each operator
     *
     */
    std::vector<std::vector<int> > fPowers;

    /**
     * @brief Linear vs quadratic EFT
     *
     */
    EFTConfig::EFTOrder fOrder;
    std::vector<EFTVariation> fEFTVariations;
    bool fProduceShapeFactorParametrization;
};

#endif


