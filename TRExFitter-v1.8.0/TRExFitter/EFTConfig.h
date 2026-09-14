#pragma once

#include <map>
#include <string>
#include <vector>

/**
 * @brief Class holding EFT configuration
 *
 */
class EFTConfig {

public:

    /**
     * @brief Linear vs quadratic parametrisation
     *
     */
    enum EFTOrder {
      ALL = 0,
      LIN = 1
    };

    /**
     * @brief Construct a new EFTConfig object
     *
     */
    EFTConfig();

    /**
     * @brief Destroy the EFTConfig object
     *
     */
    ~EFTConfig() = default;

    /**
     * @brief Add minimum and maximum for a parameter
     *
     * @param param
     * @param min
     * @param max
     */
    void AddParamMinMax(const std::string& param, const double min, const double max);

    /**
     * @brief Get the Par Min Max object
     *
     * @param param
     * @return std::pair<double, double>
     */
    std::pair<double, double> GetParMinMax(const std::string& param) const;

    /**
     * @brief Set the parametristion order
     *
     * @param order
     */
    inline void SetEFTOrder(const EFTOrder order) {fEFTOrder = order;}

    /**
     * @brief Get the parametrisation order
     *
     * @return EFTOrder
     */
    inline EFTOrder GetEFTOrder() const {return fEFTOrder;}

    /**
     * @brief Set the Split Samples Per Bin object. True = using NormFactor reparametrisation, False = using ShapeFactor reparametrisation
     *
     * @param flag
     */
    inline void SetSplitSamplesPerBin(const bool flag) {fSplitSamplesPerBin = flag;}

    /**
     * @brief Get the Split Samples Per Bin object
     *
     * @return true
     * @return false
     */
    inline bool GetSplitSamplesPerBin() const {return fSplitSamplesPerBin;}

    /**
     * @brief Set the Reference Sample object
     *
     * @param samples
     */
    inline void SetReferenceSample(const std::vector<std::string>& samples) {fReferenceSamples = samples;}

    /**
     * @brief Set the Parameters object
     *
     * @param pars
     */
    inline void SetParameters(const std::vector<std::string>& pars) {fParams = pars;}

    /**
     * @brief Set the Regions object
     *
     * @param regions
     */
    inline void SetRegions(const std::vector<std::string>& regions) {fRegions = regions;}

    /**
     * @brief Add excluded parameter, region, sample combination. Wildcards for the region and sample are supported
     *
     * @param param
     * @param region
     * @param sample
     */
    void AddExcludedRegionSample(const std::string& param,
                                 const std::string& region,
                                 const std::string& sample);

    /**
     * @brief Are all EFT parameters excluded in a given region and sample
     *
     * @param region
     * @param sample
     * @return true
     * @return false
     */
    bool IsFullyExcluded(const std::string& region, const std::string& sample) const;

    /**
     * @brief Is a given parameter excluded for a given region and sample
     *
     * @param param
     * @param region
     * @param sample
     * @return true
     * @return false
     */
    bool ParamIsExcluded(const std::string& param,
                         const std::string& region,
                         const std::string& sample) const;

    /**
     * @brief Get the formula and dependancy after the excluded param/sample/region settings
     *
     * @param originalPars
     * @param region
     * @param sample
     * @return std::vector<std::pair<std::string, std::string> >
     */
    std::vector<std::pair<std::string, std::string> > GetUpdatedParametrization(const std::vector<std::pair<std::string, std::string> >& originalPars,
                                                                                const std::string& region,
                                                                                const std::string& sample) const;

    /**
     * @brief Set the Operators Per Sample object
     *
     * @param map
     */
    inline void SetOperatorsPerSample(const std::map<std::string, std::vector<std::string> >& map) {fOperatorsForSample = map;}

    /**
     * @brief Set the flag for the ratio plot
     *
     * @param flag
     */
    inline void SetEFTPlotRatio(const bool flag) {fEFTPlotRatio = flag;}

    /**
     * @brief
     *
     * @return true
     * @return false
     */
    bool GetEFTPlotRatio() const {return fEFTPlotRatio;}

private:
    std::map<std::string, std::pair<double, double> > fParamMinMax;
    EFTOrder fEFTOrder;
    bool fSplitSamplesPerBin;
    std::vector<std::string> fReferenceSamples;
    std::vector<std::string> fParams;
    std::vector<std::string> fRegions;
    std::map<std::string,std::vector<std::string> > fOperatorsForSample;
    std::map<std::string, std::map<std::string, std::vector<std::string> > > fExcludedRegionSample;
    bool fEFTPlotRatio;
};
