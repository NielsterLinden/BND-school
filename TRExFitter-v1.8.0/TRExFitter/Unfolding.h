#ifndef UNFOLDING_H
#define UNFOLDING_H

// Unfolding includes
#include "UnfoldingCode/UnfoldingCode/FoldingManager.h"

#include <string>

class NormFactor;

class Unfolding {
public:

    enum UnfoldingChi2Type {
        EMPTY = 0,
        CHI2 = 1,
        CHI2NDF = 2,
        PROB = 3
    };

    explicit Unfolding();
    Unfolding(const Unfolding& unf) = default;
    Unfolding& operator=(const Unfolding& unf) = default;
    Unfolding(Unfolding&& unf) = default;
    Unfolding& operator=(Unfolding&& unf) = default;

    ~Unfolding() = default;

    std::string GetTitleRegionX(const std::string& regionName) const;
    std::string GetTitleRegionY(const std::string& regionName) const;

    std::string fName;
    FoldingManager::MATRIXORIENTATION fMatrixOrientation;
    std::string fTruthDistributionPath;
    std::string fTruthDistributionFile;
    std::string fTruthDistributionName;
    int fNumberUnfoldingTruthBins;
    double fUnfoldingResultMin;
    double fUnfoldingResultMax;
    std::string fUnfoldingTitleX;
    std::string fUnfoldingTitleY;
    double fUnfoldingMinRangeY;
    double fUnfoldingScaleRangeY;
    double fUnfoldingTitleOffsetY;
    double fUnfoldingTitleOffsetX;
    double fUnfoldingRatioYmin;
    double fUnfoldingRatioYmax;
    bool fUnfoldingLogX;
    bool fUnfoldingLogY;
    std::vector<std::pair<std::string, std::string> > fMigrationTitleX;
    std::vector<std::pair<std::string, std::string> > fMigrationTitleY;
    bool fMigrationLogX;
    bool fMigrationLogY;
    double fMigrationTitleOffsetX;
    double fMigrationTitleOffsetY;
    double fMigrationZmin;
    double fMigrationZmax;
    double fResponseZmin;
    double fResponseZmax;
    bool fPlotSystematicMigrations;
    bool fMigrationText;
    std::string fNominalTruthSample;
    std::string fAlternativeAsimovTruthSample;
    bool fUnfoldingDivideByBinWidth;
    double fUnfoldingDivideByLumi;
    UnfoldingChi2Type fUnfoldingChi2Type;
    bool fUnfoldNormXSec;
    int fUnfoldNormXSecBinN;
    int fErrorBandColor;
    int fStatErrorBandColor;
    int fRatioLineColor;
    int fRatioLineStyle;
    int fRatioLineWidth;
    double fLegendXposition;
    int fLegendTextSize;
    double fBreakdownRange;

    /**
     * @brief Does it have norm factor that are reparametrised?
     *
     * @param norrmFactors
     * @return true
     * @return false
     */
    bool HasReparametrisation(const std::vector<std::shared_ptr<NormFactor> >& normFactors) const;

private:
    std::string GetTitleRegion(const std::string& regionName, const bool isX) const;
};

#endif
