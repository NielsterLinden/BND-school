#ifndef NORMFACTOR_H
#define NORMFACTOR_H

/// c++ includes
#include <string>
#include <vector>

class NormFactor{
public:
    explicit NormFactor();

    explicit NormFactor(const std::string& name, double nominal=1., double min=0., double max=10., bool isConst=false, const std::string& subCategory="NormFactors");

    ~NormFactor() = default;
    NormFactor(const NormFactor& n) = delete;
    NormFactor(NormFactor&& n) = delete;
    NormFactor& operator=(const NormFactor& n) = delete;
    NormFactor& operator=(NormFactor&& n) = delete;

    void Print() const;

    /**
     * @brief Get the nominal value
     *
     * @return double
     */
    inline double GetNominal() const {return fNominal;}

    /**
     * @brief Get the Min value
     *
     * @return double
     */
    inline double GetMin() const {return fMin;}

    /**
     * @brief Get the Max value
     *
     * @return double
     */
    inline double GetMax() const {return fMax;}

    /**
     * @brief Set the nominal, min and max values
     *
     * @param nominal
     * @param min
     * @param max
     */
    void SetNominalMinMax(const double nominal, const double min, const double max);

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;
    std::string fSubCategory;

    bool fConst;

    std::vector<std::string> fRegions;
    std::vector<std::string> fExclude;

    std::pair<std::string,std::string> fExpression;

    double fTau; // for Tikhonov regularization
    std::string fCorrelateUnfolding;
    std::vector<std::string> fSamples;

private:
    double fNominal;
    double fMin;
    double fMax;

};

#endif
