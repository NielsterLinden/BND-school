#ifndef TEMPLATEMORPHER_H_
#define TEMPLATEMORPHER_H_

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "TF1.h"
#include "TF2.h"

class Morphing;
class Region;
class SampleHist;

/**
 * @brief Class used for doing the per-bin fitting for the template mophing
 *
 */
class TemplateMorpher {

public:

  /**
   * @brief 1D vs 2D fits
   *
   */
  enum FIT_DIMENSIONALITY {
    ONE_DIMENSION = 0,
    TWO_DIMENSIONS = 1
  };

  /**
   * @brief Construct a new Template Morpher object
   *
   * @param name TRExFitter fName
   * @param suffix TRExFitter fSuffix
   */
  TemplateMorpher(const std::string& name, const std::string& suffix);

  /**
   * @brief Destroy the Template Morpher object
   *
   */
  ~TemplateMorpher() = default;

  /**
   * @brief Main method that does the fitting and plotting
   *
   * @param regions list of regions
   */
  void ProcessTemplatesAndPlot(const std::vector<std::unique_ptr<Region> >& regions);

  /**
   * @brief Set the Fit Dimensionality object
   *
   * @param dim
   */
  inline void SetFitDimensionality(const FIT_DIMENSIONALITY dim) {fFitDimensionality = dim;}

  /**
   * @brief Set the Morphing Setting object
   *
   * @param setting Morphing settings
   */
  void SetMorphingSetting(const Morphing* setting);

private:

  /**
   * @brief Fit and plot one region
   *
   * @param samples List of template samples
   * @param regName Name of the region
   * @param targetName Name of the sample to apply the templates to
   * @return std::vector<std::vector<double> > fitted parameters per bin
   */
  std::vector<std::vector<double> > ProcessAndPlotOneRegion(const std::vector<std::shared_ptr<SampleHist> >& samples,
                                                            const std::string& regName,
                                                            const std::string& targetName);

  /**
   * @brief Plot and fit one bin for the 1D fit
   *
   * @param tuple For each template, the list of parameter names and values, yield, yield error
   * @param nominal Yield of the nominal template
   * @param regName Region name
   * @param bin Index of the bin
   * @param param1 Name of the parameter
   * @param targetName Name of the sample to apply the templates to
   * @return std::vector<double> fitted parameters
   */
  std::vector<double> ProcessAndPlotOneBin1D(const std::vector<std::tuple<std::vector<std::pair<std::string, double> >, double, double> >& tuple,
                                             const double nominal,
                                             const std::string& regName,
                                             const int bin,
                                             const std::string& param1,
                                             const std::string& targetName) const;

  /**
   * @brief Plot and fit one bin for the 2D fit
   *
   * @param tuple For each template, the list of parameter names and values, yield, yield error
   * @param nominal Yield of the nominal template
   * @param regName Region name
   * @param bin Index of the bin
   * @param param1 Name of the first parameter
   * @param param2 Name of the second parameter
   * @param targetName Name of the sample to apply the templates to
   * @return std::vector<double> fitted parameters
   */
  std::vector<double> ProcessAndPlotOneBin2D(const std::vector<std::tuple<std::vector<std::pair<std::string, double> >, double, double> >& tuple,
                                             const double nominal,
                                             const std::string& regName,
                                             const int bin,
                                             const std::string& param1,
                                             const std::string& param2,
                                             const std::string& targetName) const;

  /**
   * @brief The string representing the formula
   *
   * @return std::string
   */
  std::string FitFunctionString() const;

  /**
   * @brief Get the Fit Function object for the 1D fit
   *
   * @param minX minimum x value
   * @param maxX maximum x value
   * @return TF1
   */
  TF1 FitFunc(const double minX, const double maxX) const;

  /**
   * @brief Get the Fit Function object for the 2D fit
   *
   * @param minX minimum x value
   * @param maxX maximum x value
   * @param minY minimum y value
   * @param maxY maximum y value
   * @return TF2
   */
  TF2 FitFunc(const double minX, const double maxX, const double minY, const double maxY) const;

  /**
   * @brief Validate that the fit function contains correct number of fit parameters and variables
   *
   * @return bool
   */
  bool ValidateFitFunction() const;

  /**
   * @brief Turns the fitted parameter values into a formula as a string
   *
   * @param reg Region
   * @param params parameter values per bin
   * @param ostream output stream
   */
  void ParamsToString(const std::unique_ptr<Region>& reg,
                      const std::pair<std::string, std::vector<std::vector<double> > >& params,
                      std::ofstream* ostream) const;

  /**
   * @brief Name from TRExFitter fName
   *
   */
  std::string fName;

  /**
   * @brief Suffix from TRExFitter fSuffix
   *
   */
  std::string fSuffix;

  /**
   * @brief Fit dimensionality
   *
   */
  FIT_DIMENSIONALITY fFitDimensionality;

  /**
   * @brief Fit function
   *
   */
  std::string fFitFunction;

  /**
   * @brief Different Morphing settings
   *
   */
  const Morphing* fMorphingSetting;
};

#endif
