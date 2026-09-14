#pragma once

#include "TRExFitter/NormFactor.h"

#include "TH1D.h"

#include "RooWorkspace.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

class FitResults;
class RooDataSet;
class RooWorkspace;

class FitToys {
public:

  explicit FitToys(const std::string& name, const std::string& suffix, const std::string& inputName);
  ~FitToys() = default;

  inline void SetMinos(const std::vector<std::string>& minos) {fMinos = minos;}

  inline void SetIsRegularizedUnfolding(const bool flag) {fIsRegularizedUnfolding = flag;}

  inline void SetNormFactors(const std::vector<std::shared_ptr<NormFactor> >& nfs) {fNormFactors = nfs;}

  inline void SetNToys(const int toys) {fNToys = toys;}

  inline void SetRegularizationType(const int type) {fRegularizationType = type;}

  inline void SetUseAutoDiff(const bool flag) {fUseAutoDiff = flag;}

  inline void SetFitType(const int type) {fFitType = type;}

  inline void SetHistoNbins(const int bins) {fHistoNbins = bins;}

  inline void SetPOIs(const std::vector<std::string>& pois) {fPOIs = pois;}

  inline void SetFitIsBlind(const bool flag) {fFitIsBlind = flag;}

  inline void SetPOIAsimov(const std::map<std::string, double>& map) {fFitPOIAsimov = map;}

  inline void SetSeed(const int seed) {fSeed = seed;}

  inline void SetStatOnlyFluctuation(const bool flag) {fStatOnlyFluctuation = flag;}

  inline void SetRandomPOIStartingValue(const std::map<std::string, std::pair<double, double> >& map) {fRandomPOIStartingValues = map;}

  inline void SetPlotLabel(const std::string& label) {fPlotLabel = label;}

  inline void SetRunLHScanForAll(const bool flag) {fRunLHScanForAll = flag;}

  inline void SetLHScanCondition(const std::map<std::string, std::pair<double, double> >& map) {fLHScanCondition = map;}

  inline void SetFitFixedNPs(const std::map<std::string, double>& map) {fFitFixedNPs = map;}

  void RunToys(RooWorkspace* ws, FitResults* fr);

private:

  void DrawToyPullPlot(const std::vector<TH1D>& hist, TFile* out) const;

  /**
   * @brief Write uncertianty from toys to the standard txt file
   *
   * @param errors
   * @param fr
   */
  void WriteToysStat(const std::vector<double>& errors, const FitResults* fr);

  void LikelihoodScan(RooWorkspace *ws,
                      const std::string& varName,
                      RooDataSet* data,
                      const std::string& toysExtension,
                      const double min,
                      const double max);

  std::string fName;
  std::string fSuffix;
  std::string fInputName;
  std::vector<std::string> fMinos;
  bool fIsRegularizedUnfolding;
  std::vector<std::shared_ptr<NormFactor> > fNormFactors;
  int fNToys;
  int fRegularizationType;
  bool fUseAutoDiff;
  int fFitType;
  int fHistoNbins;
  std::vector<std::string> fPOIs;
  bool fFitIsBlind;
  std::map<std::string, double> fFitPOIAsimov;
  int fSeed;
  bool fStatOnlyFluctuation;
  std::map<std::string, std::pair<double, double> > fRandomPOIStartingValues;
  std::string fPlotLabel;
  bool fRunLHScanForAll;
  std::map<std::string, std::pair<double, double> > fLHScanCondition;
  std::map<std::string, double> fFitFixedNPs;
};
