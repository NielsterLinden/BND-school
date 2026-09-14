#ifndef LIMITTOYS_H
#define LIMITTOYS_H

#include <string>

class RooAbsData;
namespace RooStats{
    class ModelConfig;
    class HypoTestInverterResult;
}

class LimitToys {
public:
    explicit LimitToys();
    ~LimitToys() = default;

    LimitToys(const LimitToys& l) = delete;
    LimitToys(LimitToys&& l) = delete;
    LimitToys& operator=(const LimitToys& l) = delete;
    LimitToys& operator=(LimitToys&& l) = delete;

    inline void SetNToys(const int splusb, const int b) {
        fNtoysSplusB = splusb;
        fNtoysB = b;
    }

    inline void SetLimit(const float limit){fLimit = limit;}

    inline void SetScan(const int steps, const float min, const float max) {
        fScanSteps = steps;
        fScanMin = min;
        fScanMax = max;
    }

    inline void SetSeed(const int seed){ fToysSeed = seed;}

    inline void SetPlot(const bool flag){fPlot = flag;}

    inline void SetFile(const bool flag){fFile = flag;}

    inline void SetOutputPath(const std::string& path){fOutputPath = path;}

    inline void SetLimitSuffix(const std::string& suffix){fLimitSuffix = suffix;}

    void RunToys(RooAbsData* data,
                 RooStats::ModelConfig* mcSplusb,
                 RooStats::ModelConfig* mcB) const;

private:
    int fNtoysSplusB;
    int fNtoysB;
    float fLimit;
    int fScanSteps;
    float fScanMin;
    float fScanMax;
    int fToysSeed;
    bool fPlot;
    bool fFile;
    std::string fOutputPath;
    std::string fLimitSuffix;

    void RootOutput(RooStats::HypoTestInverterResult* result) const;
};

#endif
