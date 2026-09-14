#ifndef TRExPLOT_H
#define TRExPLOT_H


/// c++ includes
#include <memory>
#include <string>
#include <vector>

/// Forwards class declaration
class TArrow;
class TCanvas;
class TGraphAsymmErrors;
class TH1;
class TH1D;
class THStack;
class TLatex;
class TLegend;
class TLine;
class TPad;
class TF1;

class TRExPlot {
  public:

    enum class RATIOTYPE {
        DATAOVERMC = 1,
        DATAOVERB = 2,
        SOVERB = 3,
        SOVERSQRTB = 4,
        SOVERSQRTSPLUSB = 5,
        DATAOVERBWITHS = 6
    };

    explicit TRExPlot(const std::string& name="c",int canvasWidth=600,int canvasHeight=700,bool hideRatioPad=false);

    ~TRExPlot();

    TRExPlot(const TRExPlot& t) = delete;
    TRExPlot(TRExPlot&& t) = delete;
    TRExPlot& operator=(const TRExPlot& t) = delete;
    TRExPlot& operator=(TRExPlot&& t) = delete;

    void AddLabel(const std::string& name);
    void SetLumi(const std::string& name);
    void SetLumiScale(double scale);
    void SetCME(const std::string& name);
    void SetXaxis(const std::string& name,bool isNjet=false);
    void SetYaxis(const std::string& name);
    void SetYmaxScale(double scale);
    void ResizeBinLabel(const int n);
    void SetBinLabel(int bin, const std::string& name);
    void SetBinWidth(double width);

    void SetData(TH1* h,std::string name="Data");
    void AddSignal(TH1* h,std::string name="Signal");
    void AddNormSignal(TH1* h,std::string name="Signal");
    void AddOverSignal(TH1* h,std::string name="Signal", double custom_scale=1.0);
    void AddBackground(TH1* h,std::string name="MC");
    void SetTotAsym(TGraphAsymmErrors* g);
    void SetTot(TH1* h);

    void SetChi2KS(double chi2prob, double ksprob=-1., double chi2val=-1., int ndf=-1);
    void BlindData();

    void SetXaxisRange(const std::vector<double>& vec) {fXaxisRange = vec;}

    void Draw(const std::string& options="");
    void SaveAs(const std::string& name) const;
    void SaveAsBkgOnly(const std::string& name) const;
    void WriteToFile(const std::string& name) const;

    std::string fName;
    std::unique_ptr<TH1> h_data;
    std::unique_ptr<TGraphAsymmErrors> g_data;
    std::vector<std::unique_ptr<TH1> > h_bkg;
    std::vector<std::unique_ptr<TH1> > h_signal;
    std::vector<std::unique_ptr<TH1> > h_normsig;
    std::vector<std::unique_ptr<TH1> > h_oversig;
    THStack* h_stack;
    std::unique_ptr<TH1> h_tot;
    std::unique_ptr<TGraphAsymmErrors> g_tot;
    TH1* h_tot_bkg_prefit;
    TH1* h_dummy;
    std::shared_ptr<TH1> h_ratio;
    std::shared_ptr<TH1> h_ratio_pre;
    std::shared_ptr<TH1> h_tot_nosyst;
    std::shared_ptr<TGraphAsymmErrors> g_ratio;
    std::shared_ptr<TGraphAsymmErrors> g_ratio2;
    std::shared_ptr<TH1D> h_blind;
    std::shared_ptr<TH1D> h_blindratio;
    std::shared_ptr<TH1> h_dummy2;

    std::unique_ptr<TCanvas> c;
    TLegend* leg;
    TLegend* leg1;
    TPad* pad0;
    TPad* pad1;
    std::unique_ptr<TLatex> KSlab;
    std::unique_ptr<TLine> hline, vline;

    std::unique_ptr<TCanvas> fBkgOnlyCanvas;
    std::unique_ptr<TLegend> fBkgOnlyLegend;
    std::unique_ptr<TH1> fDataMinusBkg;
    std::unique_ptr<TH1> fSigStack;
    std::shared_ptr<TGraphAsymmErrors> fBkgOnlyErr;

    std::string xtitle;
    std::string ytitle;
    std::string fDataName;
    std::vector< std::string > fBkgNames;
    std::vector< std::string > fSigNames;
    std::vector< std::string > fNormSigNames;
    std::vector< std::string > fOverSigNames;
    std::vector< std::string > fLabels;
    std::string fLumi;
    std::string fCME;
    std::string fPlotLabel;
    double yMaxScale;
    int NDF;
    double Chi2val;
    double fChi2prob;
    double KSprob;

    double fYmax;
    double fYmin;
    double fRatioYmax;
    double fRatioYmin;
    double fBinWidth;
    bool fIsNjet;
    bool fShowYields;
    std::vector<std::string> fBinLabel;
    std::vector<int> fBinDividers;
    std::vector<std::string> fRangeLabels;
    double fLumiScale;
    int fLegendNColumns;
    std::vector<double> fXaxisRange;
    std::string fRatioYtitle;

    RATIOTYPE fRatioType;

    double fLabelX;
    double fLabelY;
    double fLegendX1;
    double fLegendX2;
    double fLegendY;

    std::vector<int> fBlindedBins;

    // two parameters to chose put tick marks at specific bin edges and how to label them
    std::vector<int> fBinTickMarks;
    std::vector<std::string> fBinTickMarksLabels;
    std::string fPreFitLabel;

public:

    /**
     * @brief Return the total background histogram
     *
     * @return TH1*
     */
    TH1* GetTotBkg() const;

    void SetBinBlinding(const std::vector<int>& bins);

    void DrawDataMinusBkg(const bool errorOnSignal);

private:
    std::unique_ptr<TH1D> fBlinding;

    // Customise the tick marks of the X-axis of the provided histogram, using fBinTickMarks and fBinTickMarksLabels
    // takes as input pointer h to the histogram used for the plotting in the first place
    void CustomiseTickMarks(TH1 *h);
    // when the feature is used, this function needs to exist until the plots are saved
    std::shared_ptr<TF1> fTickMarksFunction;
    std::vector<std::unique_ptr<TArrow> > fArrows;


};

// function to get asymmetric error bars for hists
double GC_up(double data);
double GC_down(double data);

/**
 * @brief Turn histogram into a graph and use poisson uncertainty per bin
 *
 * @param h
 * @return std::unique_ptr<TGraphAsymmErrors>
 */
std::unique_ptr<TGraphAsymmErrors> poissonize(const TH1 *h);

/**
 * @brief Tun Histogram into a graph and keep the errors as set in the histogram
 *
 * @param h
 * @return std::unique_ptr<TGraphAsymmErrors>
 */
std::unique_ptr<TGraphAsymmErrors> histToGraph(const TH1* h);
void SetHistBinWidth(TH1* h,double width);
void SetGraphBinWidth(TGraphAsymmErrors* g,double width);

#endif
