#pragma once

#include "TRExFitter/FittingTool.h"
#include "TRExFitter/Unfolding.h"

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

/// forward class declarations
class TH1;
class FitResults;
class NormFactor;

/**
 * @brief Class to make nice plots from the grouped impact
 *
 */
class UncertaintyPlotter {

public:

  /**
   * @brief Construct a new Uncertainty Plotter object
   *
   */
  explicit UncertaintyPlotter() noexcept;

  /**
   * @brief Destroy the Uncertainty Plotter object
   *
   */
  ~UncertaintyPlotter() = default;

  /**
   * @brief move/copy/assign
   */
  UncertaintyPlotter(const UncertaintyPlotter& rhs) = delete;
  UncertaintyPlotter(UncertaintyPlotter&& rhs) = delete;
  UncertaintyPlotter operator=(const UncertaintyPlotter& rhs) = delete;
  UncertaintyPlotter operator=(UncertaintyPlotter&& rhs) = delete;

  /**
   * @brief Set the Output Folder path
   *
   * @param folder
   */
  inline void SetOutputFolder(const std::string& folder) {m_outFolder = folder;}

  /**
   * @brief Tells the code if total or syst only uncertainty is used
   *
   * @param flag
   */
  inline void SetUseTotalUncertainty(const bool flag) {m_useTotalUncertainty = flag;}

  /**
   * @brief Pass the names of the POIs, Must match the order of the grouped impact
   *
   * @param pois
   */
  inline void SetPOIs(const std::vector<std::pair<std::string, int> >& pois) {m_pois = pois;}

  /**
   * @brief Set the Plot label string, if set to "none" will not show it at all
   *
   * @param label
   */
  inline void SetPlotLabel(const std::string& label) {m_plotLabel = label;}

  /**
   * @brief Set CME label
   *
   * @param label
   */
  inline void SetCMELabel(const std::string& label) {m_CME = label;}

  /**
   * @brief Set the lumi Label
   *
   * @param label
   */
  inline void SetLumiLabel(const std::string& label) {m_Lumi = label;}

  /**
   * @brief Tell the code if the plot is normalized
   *
   * @param flag
   */
  inline void SetIsNormalized(const bool flag) {m_isNormalized = flag;}

  /**
   * @brief Set the x-axis label
   *
   * @param label
   */
  inline void SetXaxisLabel(const std::string& label) {m_xAxisLabel = label;}

  /**
   * @brief  Plot uncertainty breakdown
   *
   * @param histo
   * @param unfolding
   * @param fitTool
   * @param wsName
   * @param wsPath
   * @param frName
   * @param ns
   */
  void PlotUncertainties(const TH1* histo,
                         const Unfolding& unfolding,
                         const FittingTool& fitTool,
                         const std::string& wsName,
                         const std::string& wsPath,
                         const std::string& frName,
                         const std::vector<std::shared_ptr<NormFactor> >& ns) const;

  /**
   * @brief Set the Categories object
   *
   * @param categories
   */
  void SetCategories(const std::vector<std::string>& categories) {m_categories = categories;}

  /**
   * @brief Get the Uncertainty object
   *
   * @param fitTool
   * @param unc
   * @return std::vector<double>
   */
  std::vector<double> GetUncertainty(const FittingTool& fitTool, const std::string& unc) const;

  /**
   * @brief Set the Sub Caterory Map object
   *
   * @param map
   */
  void SetSubCategoryMap(const std::map<std::string, std::string>& map) {m_subCategoryMap = map;}

  /**
   * @brief Set the Reparametrised Bins
   *
   * @param bins
   */
  void SetReparametrisedBins(const std::vector<int>& bins) {m_reparametrisedBins = bins;}

  /**
   * @brief Set the Custom Range object
   *
   * @param range
   */
  inline void SetCustomRange(const double range) {m_useCustomRange = true; m_customRange = range;}

  /**
   * @brief Set the Use Custom Range object
   *
   * @param flag
   */
  inline void SetUseCustomRange(const bool flag) {m_useCustomRange = flag;}

  /**
   * @brief Set other POIs (not form the bins)
   *
   * @param pois
   */
  inline void SetOtherPOIs(const std::vector<std::string>& pois) {m_otherPOIs = pois;}

private:

  /**
   * @brief A helper function to get the impact of each category. The first dimsion is the POI, second is the category
   *
   * @param fitTool
   * @return std::vector< std::vector<std::pair<std::string, double> > >
   */
  std::vector< std::vector<std::pair<std::string, double> > > GetCategoryImpact(const FittingTool& fitTool) const;

  /**
   * @brief Get the Category Impact for the reparametrised object
   *
   * @param fitTool
   * @param param
   * @param wsFile
   * @param wsName
   * @return std::vector<double>
   */
  std::vector<double > GetCategoryImpactReparametrised(const FittingTool& fitTool,
                                                       const std::string& param,
                                                       const std::string& wsFile,
                                                       const std::string& wsName) const;

  /**
   * @brief Output folder
   *
   */
  std::string m_outFolder;

  /**
   * @brief Flag to use the total uncertainty
   *
   */
  bool m_useTotalUncertainty;

  /**
   * @brief Names of the POIs
   *
   */
  std::vector<std::pair<std::string, int> > m_pois;

  /**
   * @brief Plot label
   *
   */
  std::string m_plotLabel;

  /**
   * @brief CME label
   *
   */
  std::string m_CME;

  /**
   * @brief Luminosity label
   *
   */
  std::string m_Lumi;

  /**
   * @brief Flag to tell the code if the plot is normalized
   *
   */
  bool m_isNormalized;

  /**
   * @brief x axis label
   *
   */
  std::string m_xAxisLabel;

  /**
   * @brief List of colours
   *
   */
  std::vector<int> m_colours;

  /**
   * @brief List of styles
   *
   */
  std::vector<int> m_styles;

  std::vector<int> m_reparametrisedBins;

  std::vector<std::string> m_categories;

  std::map<std::string, std::string> m_subCategoryMap;

  bool m_useCustomRange;

  double m_customRange;

  std::vector<std::string> m_otherPOIs;

};

