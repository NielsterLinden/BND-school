#ifndef LIKELIHOOHSVANMANAGER_H
#define LIKELIHOOHSVANMANAGER_H

#include "TRExFitter/FitUtils.h"

#include <memory>
#include <string>
#include <vector>

class RooArgSet;
class RooDataSet;
class RooWorkspace;

class LikelihoodScanManager {
public:

    struct Result2D {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<std::vector<double> > z;
    };

    explicit LikelihoodScanManager();
    ~LikelihoodScanManager() = default;

    LikelihoodScanManager(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager(LikelihoodScanManager&& m) = delete;
    LikelihoodScanManager& operator=(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager& operator=(LikelihoodScanManager&& m) = delete;

    inline void SetScanParamsX(const double min, const double max, const int steps, const int step) {
        fScanMinX = min;
        fScanMaxX = max;
        fStepsX = steps;
        fStepX = step;
    }

    inline void SetScanParamsY(const double min, const double max, const int steps, const int step) {
        fScanMinY = min;
        fScanMaxY = max;
        fStepsY = steps;
        fStepY = step;
    }

    inline void SetNCPU(const int n){fCPU = n;}

    inline void SetOffSet(const bool flag){fUseOffset = flag;}

    inline void SetBlindedParameters(const std::vector<std::string>& par){fBlindedParameters = par;}

    inline void SetUseNll(const bool flag) {fUseNllInLHscan = flag;}

    inline void SetExternalConstraints(const RooArgSet* const external){fExternalConstraints = external;}

    void SetFittingOptions(const FitUtils::FittingOptions& options);

    inline void SetToleranceScale(const double value) {fToleranceScale = value;}

    void SetUseAutoDiff(const bool flag);

    void SetNLLOffset(const std::string& val);

    inline float GetMinValX() const {return fMinValX;}
    inline float GetMinValY() const {return fMinValY;}
    inline float GetMaxValX() const {return fMaxValX;}
    inline float GetMaxValY() const {return fMaxValY;}

    typedef std::pair<std::vector<double>, std::vector<double> > scanResult1D;
    struct ImpactScanResult1D {
        std::vector<double> steps;
        std::vector<double> strengths;
    };

    scanResult1D Run1DScan(const RooWorkspace* ws,
                           const std::string& varName,
                           RooDataSet* data) const;

    ImpactScanResult1D Run1DImpactScan(const RooWorkspace* ws,
                                       const std::string& varNameX,
                                       const std::string& varNameY,
                                       RooDataSet* data) const;


    Result2D Run2DScan(const RooWorkspace* ws,
                       const std::pair<std::string,std::string>& varNames,
                       RooDataSet* data);

private:

    double fScanMinX;
    double fScanMinY;
    int fStepsX;
    int fStepX;
    double fScanMaxX;
    double fScanMaxY;
    int fStepsY;
    int fStepY;
    bool fUseOffset;
    int fCPU;
    bool fUseNllInLHscan;
    std::vector<std::string> fBlindedParameters;
    float fMinValX;
    float fMinValY;
    float fMaxValX;
    float fMaxValY;
    const RooArgSet* fExternalConstraints;
    FitUtils::FittingOptions m_fitOptions;
    double fToleranceScale;
    bool fUseAutoDiff;
    std::string fNLLOffset;

};

#endif
