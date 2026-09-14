#ifndef TRUTHSAMPLE_H_
#define TRUTHSAMPLE_H_

#include <memory>
#include <string>

class TH1;
class Unfolding;

class TruthSample {

public:

    explicit TruthSample(const std::string& name);
    ~TruthSample() = default;

    TruthSample(const TruthSample& s) = default;
    TruthSample& operator=(const TruthSample& s) = default;
    TruthSample(TruthSample&& s) = default;
    TruthSample& operator=(TruthSample&& s) = default;

    inline const std::string& GetName() const {return fName;};
    inline void SetTitle(const std::string& s) {fTitle = s;}
    inline const std::string& GetTitle() const {return fTitle;}
    inline void SetLineStyle(const int s) {fLineStyle = s;}
    inline int GetLineStyle() const {return fLineStyle;}
    inline void SetLineColor(const int c) {fLineColor = c;}
    inline int GetLineColor() const {return fLineColor;}
    inline void SetLineWidth(const int w) {fLineWidth = w;}
    inline int GetLineWidth() const {return fLineWidth;}
    inline void SetUseForPlotting(const bool b) {fUseForPlotting = b;}
    inline int GetUseForPlotting() const {return fUseForPlotting;}
    inline void SetTruthDistributionPath(const std::string& s) {fTruthDistributionPath = s;}
    inline void SetTruthDistributionFile(const std::string& s) {fTruthDistributionFile = s;}
    inline void SetTruthDistributionName(const std::string& s) {fTruthDistributionName = s;}

    std::unique_ptr<TH1> GetHisto(const Unfolding* unfolding) const;
    std::string fUnfoldingName;

private:

    std::string fName;
    std::string fTitle;
    int fLineStyle;
    int fLineColor;
    int fLineWidth;
    bool fUseForPlotting;
    std::string fTruthDistributionPath;
    std::string fTruthDistributionFile;
    std::string fTruthDistributionName;
};

#endif
