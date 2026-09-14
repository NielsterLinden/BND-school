#ifndef PruningUtil_H
#define PruningUtil_H

#include "TH1F.h"
#include <iostream>

#include <tuple>

class PruningUtil {
public:

    enum class SHAPEOPTION {
        MAXBIN = 0,
        KSTEST = 1
    };

    explicit PruningUtil();
    ~PruningUtil() = default;

    PruningUtil(const PruningUtil& p) = delete;
    PruningUtil(PruningUtil&& p) = delete;
    PruningUtil& operator=(const PruningUtil& p) = delete;
    PruningUtil& operator=(PruningUtil&& p) = delete;

    // 0 = sample-by-sample, 1 = relative to tot background, 2 = relative to tot S+B
    void SetStrategy(const int strat);
    inline int GetStrategy() const {return fStrategy;}
    void SetShapeOption(const SHAPEOPTION opt);
    void SetThresholdNorm(const double thres);
    void SetThresholdShape(const double thres);
    void SetThresholdIsLarge(const double thres);
    void SetRemoveSystOnEmptySample(const bool flag);

    // if hTot is not set, just prun w.r.t. only this nominal,
    // otherwise use hTot and do relative pruning
    std::tuple<bool,bool,bool,bool> CheckSystPruning(const TH1* const hUp,
                                    const TH1* const hDown,
                                    const TH1* const hNom,
                                    double& magnitude,
                                    const TH1* hTot = nullptr);

    bool HasShapeRelative(const TH1* const hNom,
                          const TH1* const hUp,
                          const TH1* const hDown,
                          const TH1* const combined,
                          const double threshold,
                          double& magnitude) const;

    bool HasShapeKS(const TH1* const hNom,
                    const TH1* const hUp,
                    const TH1* const hDown,
                    const double threshold) const;

 private:
    int fStrategy;
    SHAPEOPTION fShapeOption;
    double fThresholdNorm;
    double fThresholdShape;
    double fThresholdIsLarge;
    bool fRemoveSystOnEmptySample;
};

#endif
