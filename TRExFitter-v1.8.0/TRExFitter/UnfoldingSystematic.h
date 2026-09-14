#ifndef UNFOLDINGSYSTEMATIC_H_
#define UNFOLDINGSYSTEMATIC_H_

#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/Unfolding.h"

#include <memory>
#include <string>
#include <vector>

class Region;
class Sample;

class UnfoldingSystematic {
public:

    explicit UnfoldingSystematic();
    ~UnfoldingSystematic() = default;

    UnfoldingSystematic(const UnfoldingSystematic& s) = default;
    UnfoldingSystematic& operator=(const UnfoldingSystematic& s) = default;
    UnfoldingSystematic(UnfoldingSystematic&& s) = default;
    UnfoldingSystematic& operator=(UnfoldingSystematic&& s) = default;

    inline void SetName(const std::string& name) {fName = name;}
    inline const std::string& GetName() const {return fName;}
    inline void SetTitle(const std::string& title) {fTitle = title;}
    inline const std::string& GetTitle() const {return fTitle;}
    inline void SetType(const int type) {fType = type;}
    inline int GetType() const {return fType;}
    inline void SetOverallUp(const double value) {fOverallUp = value;}
    inline void SetOverallDown(const double value) {fOverallDown = value;}
    inline void SetSymmetrisationType(const HistoTools::SymmetrizationType type) {fSymmetrisationType = type;}
    inline void SetSmoothOption(const HistoTools::SmoothOption opt) {fSampleSmoothingOption = opt;}
    inline void SetSmoothingType(const int type) {fSmoothingType = type;}
    inline void SetHasResponse(const bool r) {fHasResponse = r;}
    inline bool GetHasResponse() const {return fHasResponse;}
    inline void SetHasAcceptance(const bool a) {fHasAcceptance = a;}
    inline bool GetHasAcceptance() const {return fHasAcceptance;}
    inline void SetReferenceSample(const std::string& sample) {fReferenceSample = sample;}
    inline const std::string& GetReferenceSample() const {return fReferenceSample;}

    std::vector<std::shared_ptr<Systematic> > ConvertToSystematic(const Region* reg,
                                                                  const Unfolding* unfolding,
                                                                  const std::string& name,
                                                                  const std::string& unfoldingSampleName,
                                                                  const std::vector<std::shared_ptr<Sample> >& samples) const;

    bool HasSubtract() const;

    std::vector<std::string> fRegions;
    std::vector<std::string> fSamples;
    std::string fCategory;
    std::string fSubCategory;
    std::vector<std::string> fResponseMatrixPathsUp;
    std::vector<std::string> fResponseMatrixPathsDown;
    std::vector<std::string> fResponseMatrixNamesUp;
    std::vector<std::string> fResponseMatrixNamesDown;
    std::vector<std::string> fResponseMatrixFolderNamesUp;
    std::vector<std::string> fResponseMatrixFolderNamesDown;
    std::vector<std::string> fResponseMatrixFilesUp;
    std::vector<std::string> fResponseMatrixFilesDown;
    std::vector<std::string> fResponseMatrixPathSuffsUp;
    std::vector<std::string> fResponseMatrixPathSuffsDown;
    std::vector<std::string> fResponseMatrixNameSuffsUp;
    std::vector<std::string> fResponseMatrixNameSuffsDown;
    std::vector<std::string> fResponseMatrixFileSuffsUp;
    std::vector<std::string> fResponseMatrixFileSuffsDown;
    std::vector<std::string> fAcceptancePathsUp;
    std::vector<std::string> fAcceptancePathsDown;
    std::vector<std::string> fAcceptanceNamesUp;
    std::vector<std::string> fAcceptanceNamesDown;
    std::vector<std::string> fAcceptanceFolderNamesUp;
    std::vector<std::string> fAcceptanceFolderNamesDown;
    std::vector<std::string> fAcceptanceFilesUp;
    std::vector<std::string> fAcceptanceFilesDown;
    std::vector<std::string> fAcceptancePathSuffsUp;
    std::vector<std::string> fAcceptancePathSuffsDown;
    std::vector<std::string> fAcceptanceNameSuffsUp;
    std::vector<std::string> fAcceptanceNameSuffsDown;
    std::vector<std::string> fAcceptanceFileSuffsUp;
    std::vector<std::string> fAcceptanceFileSuffsDown;
    std::vector<std::string> fSelectionEffPathsUp;
    std::vector<std::string> fSelectionEffPathsDown;
    std::vector<std::string> fSelectionEffNamesUp;
    std::vector<std::string> fSelectionEffNamesDown;
    std::vector<std::string> fSelectionEffFolderNamesUp;
    std::vector<std::string> fSelectionEffFolderNamesDown;
    std::vector<std::string> fSelectionEffFilesUp;
    std::vector<std::string> fSelectionEffFilesDown;
    std::vector<std::string> fSelectionEffPathSuffsUp;
    std::vector<std::string> fSelectionEffPathSuffsDown;
    std::vector<std::string> fSelectionEffNameSuffsUp;
    std::vector<std::string> fSelectionEffNameSuffsDown;
    std::vector<std::string> fSelectionEffFileSuffsUp;
    std::vector<std::string> fSelectionEffFileSuffsDown;
    std::vector<std::string> fMigrationPathsUp;
    std::vector<std::string> fMigrationPathsDown;
    std::vector<std::string> fMigrationNamesUp;
    std::vector<std::string> fMigrationNamesDown;
    std::vector<std::string> fMigrationFolderNamesUp;
    std::vector<std::string> fMigrationFolderNamesDown;
    std::vector<std::string> fMigrationFilesUp;
    std::vector<std::string> fMigrationFilesDown;
    std::vector<std::string> fMigrationPathSuffsUp;
    std::vector<std::string> fMigrationPathSuffsDown;
    std::vector<std::string> fMigrationNameSuffsUp;
    std::vector<std::string> fMigrationNameSuffsDown;
    std::vector<std::string> fMigrationFileSuffsUp;
    std::vector<std::string> fMigrationFileSuffsDown;

    bool fHasUpVariation;
    bool fHasDownVariation;
    bool fSampleSmoothing;
    std::string fNuisanceParameter;

    std::string fCombineName;
    Systematic::COMBINATIONTYPE fCombineType;
    std::string fUnfoldingName;

    std::vector<std::string> fResponseMatrixPathsUpSubtractSample;
    std::vector<std::string> fResponseMatrixPathsDownSubtractSample;
    std::vector<std::string> fResponseMatrixFilesUpSubtractSample;
    std::vector<std::string> fResponseMatrixFilesDownSubtractSample;
    std::vector<std::string> fResponseMatrixNamesUpSubtractSample;
    std::vector<std::string> fResponseMatrixNamesDownSubtractSample;
    std::vector<std::string> fResponseMatrixFolderNamesUpSubtractSample;
    std::vector<std::string> fResponseMatrixFolderNamesDownSubtractSample;
    std::vector<std::string> fResponseMatrixPathSuffsUpSubtractSample;
    std::vector<std::string> fResponseMatrixPathSuffsDownSubtractSample;
    std::vector<std::string> fResponseMatrixFilesSuffsUpSubtractSample;
    std::vector<std::string> fResponseMatrixFilesSuffsDownSubtractSample;
    std::vector<std::string> fResponseMatrixNamesSuffsUpSubtractSample;
    std::vector<std::string> fResponseMatrixNamesSuffsDownSubtractSample;
    std::vector<std::string> fMigrationPathsUpSubtractSample;
    std::vector<std::string> fMigrationPathsDownSubtractSample;
    std::vector<std::string> fMigrationFilesUpSubtractSample;
    std::vector<std::string> fMigrationFilesDownSubtractSample;
    std::vector<std::string> fMigrationNamesUpSubtractSample;
    std::vector<std::string> fMigrationNamesDownSubtractSample;
    std::vector<std::string> fMigrationFolderNamesUpSubtractSample;
    std::vector<std::string> fMigrationFolderNamesDownSubtractSample;
    std::vector<std::string> fMigrationPathSuffsUpSubtractSample;
    std::vector<std::string> fMigrationPathSuffsDownSubtractSample;
    std::vector<std::string> fMigrationFilesSuffsUpSubtractSample;
    std::vector<std::string> fMigrationFilesSuffsDownSubtractSample;
    std::vector<std::string> fMigrationNamesSuffsUpSubtractSample;
    std::vector<std::string> fMigrationNamesSuffsDownSubtractSample;
    std::vector<std::string> fSelectionEffPathsUpSubtractSample;
    std::vector<std::string> fSelectionEffPathsDownSubtractSample;
    std::vector<std::string> fSelectionEffFilesUpSubtractSample;
    std::vector<std::string> fSelectionEffFilesDownSubtractSample;
    std::vector<std::string> fSelectionEffNamesUpSubtractSample;
    std::vector<std::string> fSelectionEffNamesDownSubtractSample;
    std::vector<std::string> fSelectionEffFolderNamesUpSubtractSample;
    std::vector<std::string> fSelectionEffFolderNamesDownSubtractSample;
    std::vector<std::string> fSelectionEffPathSuffsUpSubtractSample;
    std::vector<std::string> fSelectionEffPathSuffsDownSubtractSample;
    std::vector<std::string> fSelectionEffFilesSuffsUpSubtractSample;
    std::vector<std::string> fSelectionEffFilesSuffsDownSubtractSample;
    std::vector<std::string> fSelectionEffNamesSuffsUpSubtractSample;
    std::vector<std::string> fSelectionEffNamesSuffsDownSubtractSample;
    std::vector<std::string> fAcceptancePathsUpSubtractSample;
    std::vector<std::string> fAcceptancePathsDownSubtractSample;
    std::vector<std::string> fAcceptanceFilesUpSubtractSample;
    std::vector<std::string> fAcceptanceFilesDownSubtractSample;
    std::vector<std::string> fAcceptanceNamesUpSubtractSample;
    std::vector<std::string> fAcceptanceNamesDownSubtractSample;
    std::vector<std::string> fAcceptanceFolderNamesUpSubtractSample;
    std::vector<std::string> fAcceptanceFolderNamesDownSubtractSample;
    std::vector<std::string> fAcceptancePathSuffsUpSubtractSample;
    std::vector<std::string> fAcceptancePathSuffsDownSubtractSample;
    std::vector<std::string> fAcceptanceFilesSuffsUpSubtractSample;
    std::vector<std::string> fAcceptanceFilesSuffsDownSubtractSample;
    std::vector<std::string> fAcceptanceNamesSuffsUpSubtractSample;
    std::vector<std::string> fAcceptanceNamesSuffsDownSubtractSample;
    std::vector<std::string> fDropShapeIn;
    std::vector<std::string> fDropNorm;

    double fScaleUp;
    double fScaleDown;

    std::map<std::string,double> fScaleUpRegions;
    std::map<std::string,double> fScaleDownRegions;

private:
    std::string fName;
    std::string fTitle;
    int fType;
    HistoTools::SymmetrizationType fSymmetrisationType;
    HistoTools::SmoothOption fSampleSmoothingOption;
    int fSmoothingType;
    bool fHasResponse;
    bool fHasAcceptance;
    std::string fReferenceSample;
    double fOverallUp;
    double fOverallDown;

};

#endif
