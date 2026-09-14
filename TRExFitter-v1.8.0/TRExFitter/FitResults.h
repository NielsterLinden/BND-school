#ifndef FITRESULTS_H
#define FITRESULTS_H

/// c++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

/// Forward class declaration
class CorrelationMatrix;
class NormFactor;
class NuisParameter;
class ShapeFactor;
class RooFitResult;

class FitResults {
public:
    explicit FitResults();

    ~FitResults();

    FitResults(const FitResults& f) = delete;
    FitResults(FitResults&& f) = delete;
    FitResults& operator=(const FitResults& f) = delete;
    FitResults& operator=(FitResults&& f) = delete;

    //
    // Functions
    //
    void AddNuisPar(NuisParameter *par);

    /**
     * @brief Get the NP value
     *
     * @param p name of the NP
     * @return double
     */
    double GetNuisParValue(const std::string& p) const;

    /**
     * @brief Get the NP uncertainty up
     *
     * @param p name of the NP
     * @return double
     */
    double GetNuisParErrUp(const std::string& p) const;

    /**
     * @brief Get the NP uncertainty down
     *
     * @param p name of the NP
     * @return double
     */
    double GetNuisParErrDown(const std::string& p) const;

    /**
     * @brief Read RooFitResults from a ROOT file
     *
     * @param fileName
     * @return true
     * @return false
     */
    bool ReadFromRootFile(const std::string& fileName);

    /**
     * @brief  Enum to hold the pull significance combos
     *
     */
    enum class PullSigType {Regular, Rank, NoPullSig};

    /**
     * @brief
     *
     * @param path
     * @param category
     * @param normFactors
     * @param shapeFactors
     * @param blinded
     * @param pullSigType
     */
    void DrawNPPulls(const std::string &path,
                     const std::string &category,
                     const std::vector<std::shared_ptr<NormFactor> > &normFactors,
                     const std::vector<std::shared_ptr<ShapeFactor> > &shapeFactors,
                     const std::vector<std::string>& blinded,
                     const PullSigType& pullSigType) const;

    void DrawNormFactors(const std::string &path, const std::vector<std::shared_ptr<NormFactor> > &normFactor, const std::vector<std::string>& blinded ) const;

    /**
     * @brief Draw gammas or shape factors
     *
     * @param path
     * @param blinded
     * @param shapeFactors
     * @param isShape
     */
    void DrawGammaShapePulls(const std::string &path,
                             const std::vector<std::string>& blinded,
                             const std::vector<std::shared_ptr<ShapeFactor> >& shapeFactors,
                             const bool isShape) const;

    /**
     * @brief Draw the correlation matrix
     *
     * @param path
     * @param useHEPDataFormat
     * @param corrMin
     * @param EFTonly
     * @param POIs
     */
    void DrawCorrelationMatrix(const std::string& path, const bool useHEPDataFormat, const double corrMin = -1., bool EFTonly=false, const std::vector<std::string>& POIs={} );

    /**
     * @brief Set POI precision
     *
     * @param precision
     */
    void SetPOIPrecision(const int& precision){fPOIPrecision = precision;}

    /**
     * @brief Set the Plot Label object
     *
     * @param l
     */
    void SetPlotLabel(const std::string& l){fPlotLabel = l;}

    /**
     * @brief Get the Correlation Matrix object
     *
     * @return CorrelationMatrix*
     */
    CorrelationMatrix* GetCorrelationMatrix() const {return fCorrMatrix.get();}

    /**
     * @brief Get the Nuisance Parameters object
     *
     * @return const std::map<std::string, std::shared_ptr<NuisParameter> >&
     */
    const std::map<std::string, std::shared_ptr<NuisParameter> >& GetNuisanceParameters() const {return fNuisPar;}

    /**
     * @brief Get the Output Folder object
     *
     * @return const std::string&
     */
    const std::string& GetOutputFolder() const {return fOutFolder;}

    /**
     * @brief Set the Output Folder object
     *
     * @param folder
     */
    void SetOutputFolder(const std::string& folder) {fOutFolder = folder;}

    /**
     * @brief Set the Pars To Hide object
     *
     * @param params
     */
    void SetParsToHide(const std::vector<std::string>& params) {fNuisParToHide = params;}

    /**
     * @brief Set EFT NFs
     *
     * @param params
     */
    void SetEFTNFs(const std::vector<std::string>& params) {fEFTNFs = params;}

    /**
     * @brief  Clear NP list
     *
     */
    void ClearNuisParList() {fNuisParList.clear();}

    /**
     * @brief Get the Nuis Par List object
     *
     * @return const std::vector<std::string>&
     */
    const std::vector<std::string>& GetNuisParList() const {return fNuisParList;}

    /**
     * @brief  Add to the nuis par list
     *
     * @param param
     */
    void AddToNuisParList(const std::string& param) {fNuisParList.emplace_back(param);}

    /**
     * @brief Get the RooFitResult object
     *
     * @return const RooFitResult* const
     */
    const RooFitResult* GetRooFitResult() const {return fRoofitFitResult.get();}

private:

    NuisParameter* GetNuisanceParameter(const std::string& param) const;

    //
    // Data members
    //
    std::vector<std::string> fNuisParToHide; // NPs to hide
    std::vector<std::string> fNuisParList; // NPs to show, ordered
    std::vector<std::string> fEFTNFs; // Only EFT NFs

    std::map<std::string, std::shared_ptr<NuisParameter> > fNuisPar;
    std::unique_ptr<CorrelationMatrix> fCorrMatrix;
    std::unique_ptr<RooFitResult> fRoofitFitResult;
    std::string fOutFolder;

    int fPOIPrecision;
    std::string fPlotLabel;
};

#endif
