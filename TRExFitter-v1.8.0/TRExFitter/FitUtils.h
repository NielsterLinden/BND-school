#ifndef FITUTILS_H
#define FITUTILS_H

#include "xRooFit/xRooNode.h"

#include "RooDataSet.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

class FittingTool;
class FitResults;
class NormFactor;
class Region;
class RooArgSet;
class RooRealVar;
class RooSimultaneous;
class RooWorkspace;
class ShapeFactor;
class TDirectory;
class TFile;

namespace RooStats {
    class ModelConfig;
}

namespace FitUtils {

struct FittingOptions {
    bool noGammas;
    bool noSystematics;
    std::vector<std::string> initialNPs;
    std::vector<double> initialNPvalues;
    std::vector<std::string> constNPs;
    std::vector<double> constNPvalues;
    bool randomize;
    double randomValue;
    bool constPOI;
    std::vector<std::pair<std::string, double> > poiVals;
    std::vector<std::string> constNFs;
};

/**
 * Implementing external constraints on the workspace
 * @param workspace
 * @param Fitting tool
 * @param the PDF
 * @param vector of norma factors
 */
void ApplyExternalConstraints(RooWorkspace* ws,
                              FittingTool* fitTool,
                              RooSimultaneous* simPdf,
                              const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                              int type = 0);

/**
 * A helper function to set BinnedLikelihoodOptimisation
 * @param workspace
 */
void SetBinnedLikelihoodOptimisation(RooWorkspace* ws);

/**
 * A helper function to injects values of NPs
 * @param workspace
 * @param map of NP values
 */
void InjectGlobalObservables(RooWorkspace* ws, const std::map< std::string, double >& npValues);

/**
 * @brief A helper function to injects values of NPs
 *
 * @param ws
 * @param fr
 */
void InjectGlobalObservables(RooWorkspace* ws, const FitResults* fr);

/**
 * A helper function to set POI in a file
 * @param path to the file
 * @param WS name
 * @param POI
 */
void SetPOIinFile(const std::string& path, const std::string& wsName, const std::string& poi);

/**
 * A helper function to set POI in a file
 * @param ModelConfig
 * @return parameters
 */
std::vector<std::string> GetAllParameters(const RooStats::ModelConfig* mc);

/**
 * A helper function to get number of free parameters
 * @param ModelConfig
 * @return number of free parameter
 */
std::size_t NumberOfFreeParameters(const RooStats::ModelConfig* mc);

/**
 * A helper function to set all parameters to const except the specified ones
 * @param ModelConfig
 * @param names of the exceptions
 */
void FixAllParameters(RooStats::ModelConfig* mc, const std::vector<std::string>& except);

/**
 * A helper function to set all parameters to floating
 * @param ModelConfig
 */
void FloatAllParameters(RooStats::ModelConfig* mc);

/**
 * A helper function to propagate Expression using Roofit
 * @param Workspace
 * @param File with roofit results
 * @param name of the parameter to be calculated
 * @param statOnly flag if the uncertainties are stat only
 * @param freeParams normFactors + unfolding params
 * @return Returns a vector of doubles with a following convention:
 *         - element 0: nominal
 *         - element 1: up
 *         - element 2: down
 *         - element 3: error
 * Return 0 size vector if parameter cannot be found
 */
std::vector<double> CalculateExpressionRoofit(RooWorkspace* ws,
                                              TFile* fitResultFile,
                                              const std::string& name,
                                              const bool statOnly,
                                              const std::vector<std::string>& freeParams);

/**
 * A helper function to get RooArgSet without some parameter
 * @param MC
 * @param names of the morph parameters
 */
RooArgSet GetPOIsWithout(const RooStats::ModelConfig* mc, const std::vector<std::string>& names);

/**
 * @brief Get the External Constraints object
 *
 * @param ws
 * @param normFactors
 * @param regularizationType
 * @return const RooArgSet*
 */
const RooArgSet* GetExternalConstraints(RooWorkspace* ws,
                                        const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                                        const int regularizationType);

/**
 * @brief Fix parameters to predefined values
 *
 * @param mc
 * @param params
 */
void FixParameters(RooStats::ModelConfig* mc,
                   const std::vector<std::pair<std::string, double> >& params);

/**
 * @brief Code to change the interpolation setting for a workspace
 *
 * @param ws workspace
 * @param intCode interpolation code
 */
void ChangeInterpolationCode(RooWorkspace* ws, const int intCode);

/**
 * @brief Reparametrise the shape factors based on the custom formulae
 *
 * @param ws
 * @param shapeFactors
 * @param regions
 * @param regionUsedInFit
 */
void ReparametrizeShapeFactors(RooWorkspace* ws,
                               const std::vector<std::shared_ptr<ShapeFactor> >& shapeFactors,
                               const std::vector<Region* >& regions,
                               const std::vector<std::string>& regionsUsedInFit);

/**
 * @brief Set starting values and fixes parameters based on the config settings
 *
 * @param mc Model Config
 * @param options options
 */
void SetStartingAndConstParams(RooStats::ModelConfig* mc, const FittingOptions& options);

/**
 * @brief Get the Vector of POIs
 *
 * @param model
 * @return std::vector<RooRealVar*>
 */
std::vector<RooRealVar*> GetVectorPOI(const RooStats::ModelConfig* model);

/**
 * @brief Get the goodness of fit value from the saturated model
 *
 * @param node
 * @param dataName
 * @return double
 */
double SaturatedModelGoF(std::shared_ptr<xRooNode> node, const std::string& dataName);

/**
 * @brief Get the reduced set of samples or a single sample
 *
 * @param node
 * @param samples
 * @return xRooNode
 */
xRooNode GetReducedOrSingleSample(std::shared_ptr<xRooNode> node, const std::string& samples);

/**
 * @brief Get NPs without NFs and POIs
 *
 * @param mc
 * @param nfs
 * @param pois
 * @return RooArgSet
 */
RooArgSet GetNPsWithoutNFs(const RooStats::ModelConfig* mc,
                           const std::vector<std::shared_ptr<NormFactor> >& nfs,
                           const std::vector<std::string>& pois);

/**
 * @brief Get the Propagated Covariance object
 *
 * @param f The expression to be evaluated, i.e. the RooFormulaVar object
 * @param p
 * @param fr Fit results object
 * @param nset
 * @return double
 */
double GetPropagatedCovariance(const RooAbsReal &f, const RooRealVar &p, const RooFitResult &fr, const RooArgSet &nset);

/**
 * @brief Fix parameters in xRooFit node
 *
 * @param node
 * @param toFix
 */
void FixParametersInXRooFit(const xRooNode* node, const std::map<std::string, double>& toFix);

}

#endif
