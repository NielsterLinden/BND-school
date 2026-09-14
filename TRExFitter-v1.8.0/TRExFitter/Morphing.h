#ifndef MORPHING_H
#define MORPHING_H

#include <map>
#include <string>

/**
 * @brief Class that hold several settings for the template morphing
 *
 */
class Morphing {
public:

    /**
     * @brief Construct a new Morphing object
     *
     */
    explicit Morphing();

    /**
     * @brief Destroy the Morphing object
     *
     */
    ~Morphing() = default;

    /**
     * @brief Set the fit function
     *
     * @param f Fit function
     */
    inline void SetMorphingFunction(const std::string f) {fFitFunction = f;}

    /**
     * @brief Get the fit function
     *
     * @return std::string
     */
    std::string GetMorphingFunction() const {return fFitFunction;};

    /**
     * @brief Add a starting parameter value
     *
     * @param target Template target
     * @param region Region name
     * @param bin Bin index starting from 0
     * @param index Parameter index
     * @param value starting value
     */
    void AddStartValues(const std::string& target, const std::string& region, const int bin, const int index, const double value);

    /**
     * @brief Add parameter limits
     *
     * @param target Template target
     * @param region Region name
     * @param bin Bin index starting from 0
     * @param index Parameter index
     * @param mim minimum
     * @param max maximum
     */
    void AddParamLimits(const std::string& target, const std::string& region, const int bin, const int index, const double mim, const double max);

    /**
     * @brief Check if a starting value vas set
     *
     * @param target Template target
     * @param region Region name
     * @param bin Bin index starting from 0
     * @param index Parameter index
     * @return true
     * @return false
     */
    bool HasStartValue(const std::string& target, const std::string& region, const int bin, const int index) const;

    /**
     * @brief Check if parameter limits were set
     *
     * @param target Template target
     * @param region Region name
     * @param bin Bin index starting from 0
     * @param index Parameter index
     * @return true
     * @return false
     */
    bool HasParamLimits(const std::string& target, const std::string& region, const int bin, const int index) const;

    /**
     * @brief Get parameter starting value
     *
     * @param target Template target
     * @param region Region name
     * @param bin Bin index starting from 0
     * @param index Parameter index
     * @return double
     */
    double ParamStartValue(const std::string& target, const std::string& region, const int bin, const int index) const;

    /**
     * @brief Get parameter limits
     *
     * @param target Template target
     * @param region Region name
     * @param bin Bin index starting from 0
     * @param index Parameter index
     * @return std::pair<double, double> min | max
     */
    std::pair<double, double> ParamLimits(const std::string& target, const std::string& region, const int bin, const int index) const;

    /**
     * @brief Set the scale of the range used in the parameter limits
     *
     * @param range
     */
    inline void SetRangeScale(const double range) {fRangeScale = range;}

    /**
     * @brief Get the scale of the range used in the parameter limits
     *
     * @return double
     */
    inline double GetRangeScale() const {return fRangeScale;}

    /**
     * @brief Set the Chi 2 Warning Limit
     *
     * @param limit
     */
    inline void SetChi2WarningLimit(const double limit) {fChi2WarningLimit = limit;}

    /**
     * @brief Get the Chi 2 Warning Limit
     *
     * @return double
     */
    inline double GetChi2WarningLimit() const {return fChi2WarningLimit;}

    /**
     * @brief Set bool if 1D projections of 2D parametrisation should be made
     *
     * @param make1DProjections
     */
    inline void SetMake1DProjections(const bool make1DProjections) {fMake1DProjections= make1DProjections;}

    /**
     * @brief Get bool if 1D projections of 2D parametrisation should be made
     *
     * @return bool
     */
    inline bool GetMake1DProjections() const {return fMake1DProjections;}

    /**
     * @brief Add an expression to the map
     *
     * @param param
     * @param expression
     */
    void AddExpression(const std::string& param, const std::pair<std::string, std::string>& expression);

    /**
     * @brief Get the Expressions object
     *
     * @return const std::map<std::string, std::pair<std::string, std::string> >
     */
    const std::map<std::string, std::pair<std::string, std::string> >& GetExpressions() const {return fExpressions;}

private:

    /**
     * @brief Fit function
     *
     */
    std::string fFitFunction;

    /**
     * @brief Starting values. target | region | bin index | parameter index | value
     *
     */
    std::map<std::string, std::map<std::string, std::map<int, std::map<int, double> > > > fStartValues;

    /**
     * @brief Parameter limits. target | region | bin index | paramter index | min | max
     *
     */
    std::map<std::string, std::map<std::string, std::map<int, std::map<int, std::pair<double,double> > > > > fParamLimits;

    /**
     * @brief Range scale
     *
     */
    double fRangeScale;

    /**
     * @brief Chi2 warning limit
     *
     */
    double fChi2WarningLimit;

    /**
     * @brief Make 1D projection plots for 2D parametrisation
     *
     */
    bool fMake1DProjections;

    /**
     * @brief Map of expressions for reparametrisation
     *
     */
    std::map<std::string, std::pair<std::string, std::string> > fExpressions;
};
#endif
