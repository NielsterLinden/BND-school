#ifndef HISTOTOOLS_H
#define HISTOTOOLS_H

#include "TH1D.h"

/// c++ includes
#include <memory>
#include <string>
#include <vector>

/// Foward class declaration
class TH1;
class SystematicHist;

namespace HistoTools {

    /// Enum that stores different symmetrization options
    enum SymmetrizationType {
        NOSYMMETRIZATION = 0, /// no symmetrization applied
        SYMMETRIZEONESIDED = 1, ///symmetrize one-sided systematics (e.g. JER)
        SYMMETRIZETWOSIDED = 2, /// symmetrize two-sided systematics (protects from statistical fluctuations)
        SYMMETRIZEABSMEAN = 3, /// symmetrize two-sided systematics by taking mean of the absolute values from up and down shift (protects from statistical fluctuations)
        SYMMETRIZEMAXIMUM = 4, /// symmetrize two-sided systematics by taking the larger variation from up and down and symmetrizing (protects from statistical fluctuations)
        SYMMETRIZEMAXIMUMKEEPSIGN = 5 /// symmetrize two-sided systematic by taking scaling the variatio to the max of the up and down but keeping the sign
    };

    /// Enum that stores different smoothing level options
    enum SmoothLevel {
        SMOOTHDEPENDENT = 10,
        SMOOTHINDEPENDENT = 100,
        UNKNOWN = 1000
    };

    /// Enum that stores different smoothing options
    enum SmoothOption {
        MAXVARIATION = 0,
        TTBARRESONANCE = 1,
        COMMONTOOLSMOOTHMONOTONIC = 2,
        COMMONTOOLSMOOTHPARABOLIC = 3,
        KERNELRATIOUNIFORM = 4,
        KERNELDELTAGAUSS = 5,
        KERNELRATIOGAUSS = 6,
        TCHANNEL = 7
    };

    enum class FORCESHAPETYPE {
        NOSHAPE = 0,
        LINEAR = 1,
        TRIANGULAR
    };

    /**
    * In RooStats, input histogram variable binning is not supported => convert to a constant binning
    * by creating an histogram with the same number of bins but with constant binning between 0 and 1
    * - now in case some bins are < 0 (due to bin drop functionality), they are ignored for the regBin histos
    * @param originalHist Original histogram
    * @param blindedBins List of blinded bins
    * @param scale per-bin scale
    * @param useRegularBinning flag to tell code to use regular binning
    * @return Trasformed histogram
    */
    std::unique_ptr<TH1D> TransformHistogramBinning(const TH1* originalHist,
                                                    const std::vector<int>& blindedBins,
                                                    const std::vector<double>& scale,
                                                    const bool useRegularBinning);

    /**
     * A helper function to smooth and symmetrize histograms
     * @param level of smoothing
     * @param type of symmetrization
     * @param nominal histogram
     * @param up variation histogram
     * @param down variation histogram
     * @param symmetrized/smoothed up variation
     * @param symmetrized/smoothed down variation
     * @param scale for the up variation
     * @param scale for the down varaition
     * @param smoothing option
     * @param apply ttbar resonance smoothing
     */
    void ManageHistograms(const int smoothingLevel,
                          const SymmetrizationType& symType,
                          const TH1* hNom,
                          const TH1* originUp,
                          const TH1* originDown,
                          std::unique_ptr<TH1> &modifiedUp,
                          std::unique_ptr<TH1> &modifiedDown,
                          double scaleUp,
                          double scaleDown,
                          const SmoothOption &smoothOpt);

    /**
     * A helper function to symmetrize histograms
     * @param symmetrization type
     * @param nominal histogram
     * @param up variation histogram
     * @param down variation histogram
     * @param symmetrized/smoothed up variation
     * @param symmetrized/smoothed down variation
     * @param scale for the up variation
     * @param scale for the down varaition
     */
    void SymmetrizeHistograms(const SymmetrizationType& symType,
                              const TH1* hNom,
                              const TH1* originUp,
                              const TH1* originDown,
                              std::unique_ptr<TH1> &modifiedUp,
                              std::unique_ptr<TH1> &modifiedDown,
                              double scaleUp,
                              double scaleDown);

    /**
     * A helper function to smooth histograms
     * @param level of smoothing
     * @param nominal histogram
     * @param up variation histogram
     * @param down variation histogram
     * @param symmetrized/smoothed up variation
     * @param symmetrized/smoothed down variation
     * @param smooth option
     * @param apply ttbar resonance smoothing
     */
    void SmoothHistograms(const int smoothingLevels,
                          const TH1* hNom,
                          const TH1* originUp,
                          const TH1* originDown,
                          std::unique_ptr<TH1> &modifiedUp,
                          std::unique_ptr<TH1> &modifiedDown,
                          const SmoothOption &smoothOpt);

    /**
     * A helper function to Symmetrize systematics using one sided method
     * @param Nominal histogram
     * @param Systematic histogram
     * @param Flag for up/down variation
     * @return Modified histogram
     */
    std::unique_ptr<TH1D> SymmetrizeOneSided(const TH1* const h_nominal, const TH1* const h_syst, bool &isUp );

    /**
     * A helper function to invert one sided shift
     * @param Systematic histogram
     * @param Nominal histogram
     * @return Modified histogram
     */
    std::unique_ptr<TH1D> InvertShift(const TH1* const h_syst, const TH1* const h_nominal);

    /**
     * A helper function to calculate separation between histograms
     * @param First histogram
     * @param Second histogram
     * @return Separation value
     */
    double Separation(const TH1* const h1, const TH1* const h2);

    /**
     * A helper function to calculate symmetrization using two sided variation
     * @param Original variation up histogram
     * @param Original variation down histogram
     * @param Nominal histogram
     * @return Up variation that is calculated, will be ionverted for down  variation later
     */
    std::unique_ptr<TH1D> SymmetrizeTwoSided(const TH1* const var1, const TH1* const var2, const TH1* const hnom);

    /**
     * A helper function to calculate symmetrization using mean of the absolute values
     * @param Original variation up histogram
     * @param Original variation down histogram
     * @param Nominal histogram
     * @return Up variation that is calculated, will be ionverted for down  variation later
     */
    std::unique_ptr<TH1D> SymmetrizeAbsMean(const TH1* const var1, const TH1* const var2, const TH1* const hnom);

    /**
     * A helper function to calculate symmetrization using maximum variation
     * @param Original variation up histogram
     * @param Original variation down histogram
     * @param Nominal histogram
     * @return Up variation that is calculated, will be ionverted for down  variation later
     */
    std::unique_ptr<TH1D> SymmetrizeMaximum(const TH1* const var1, const TH1* const var2, const TH1* const hnom);
    
    /**
     * A helper function to calculate symmetrization using maximum variation but preserving the sign
     * @param Original variation up histogram
     * @param Original variation down histogram
     * @param Nominal histogram
     * @return Up variation that is calculated, will be ionverted for down  variation later
     */
    std::unique_ptr<TH1D> SymmetrizeMaximumKeepSign(const TH1* const var1, const TH1* const var2, const TH1* const hnom);

    /**
     * A helper function to print warning if the two variations have the same shift and TwoSided
     * symmetrization is used.
     * @param Original variation up
     * @param Original variation down
     * @param Nominal histogram
     * @param Modified syst variation up
     * @param Modifiec syst variation down
     */
    void CheckSameShift(const TH1* const var1, const TH1* const var2, const TH1* const hnom,
    const TH1D* const tmp1, const TH1D* const tmp2);

    /**
     * A helper function to scale systeamtic histograms
     * @param Systematic histogram that will be scaled
     * @param Npminal histogram
     * @param scale factor
     */
    void Scale(TH1* h_syst, TH1* h_nominal, double factor);

    /**
     * A helper function to check for various weird features
     * @param Nominal histogram
     * @param Systematic histogram
     * @param Check for nullpointer
     * @param Crash when histo is nullptr?
     * @return Histo is OK
     */
    bool CheckHistograms(TH1* nom, SystematicHist* sh, bool checkNull = true, bool causeCrash = false);

    /**
      * A helper function to manually modify the shape of the systematic variation
      * @param Systematic variation
      * @param Nominal
      * @param type of the shape
      */ 
    void ForceShape(TH1* syst, const TH1* nominal, const FORCESHAPETYPE type); 

    /**
      * A helper function to force linear shape
      * @param Systematic variation
      * @param Nominal
      */ 
    void ForceShapeLinear(TH1* syst, const TH1* nominal);

    /**
      * A helper function to force triangular shape
      * @param Systematic variation
      * @param Nominal
      */ 
    void ForceShapeTriangular(TH1* syst, const TH1* nominal);

    /**
     * @brief scale a histograms bin by bin 
     * 
     * @param h histogram to be updated 
     * @param scales scale per bin
     */
    void ScaleBinByBin(TH1* h, const std::vector<double>& scales);

}
#endif
