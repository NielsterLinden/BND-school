#include "TRExFitter/UnfoldingSystematic.h"

#include "TRExFitter/Logger.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/Unfolding.h"

UnfoldingSystematic::UnfoldingSystematic() :
    fCategory(""),
    fSubCategory(""),
    fHasUpVariation(false),
    fHasDownVariation(false),
    fSampleSmoothing(false),
    fNuisanceParameter(""),
    fCombineName(""),
    fCombineType(Systematic::COMBINATIONTYPE::ENVELOPE),
    fUnfoldingName(""),
    fScaleUp(1.0),
    fScaleDown(1.0),
    fName(""),
    fTitle(""),
    fType(0),
    fSymmetrisationType(HistoTools::SymmetrizationType::NOSYMMETRIZATION),
    fSampleSmoothingOption(HistoTools::SmoothOption::MAXVARIATION),
    fSmoothingType(0),
    fHasResponse(false),
    fHasAcceptance(false),
    fReferenceSample(""),
    fOverallUp(0.),
    fOverallDown(0.)
{
}

std::vector<std::shared_ptr<Systematic> > UnfoldingSystematic::ConvertToSystematic(const Region* reg,
                                                                                   const Unfolding* unfolding,
                                                                                   const std::string& name,
                                                                                   const std::string& unfoldingSampleName,
                                                                                   const std::vector<std::shared_ptr<Sample> >& samples) const {

    if (samples.size() < 1) {
        LOG(ERROR) << "Samples size < 1\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::shared_ptr<Systematic> > result;
    for (int ibin = 0; ibin < unfolding->fNumberUnfoldingTruthBins; ++ibin) {
        const std::string sampleName = unfolding->fName + "_" + reg->fName + "_Truth_bin_" + std::to_string(ibin+1);

        std::shared_ptr<Systematic> syst = std::make_shared<Systematic>(fName, fType);
        if(fNuisanceParameter != ""){
            syst->fNuisanceParameter = fNuisanceParameter;
            TRExFitter::NPMAP[syst->fName] = syst->fNuisanceParameter;
        }
        else{
            syst->fNuisanceParameter = syst->fName;
            TRExFitter::NPMAP[syst->fName] = syst->fName;
        }
        syst->fStoredName = fName;
        syst->fTitle = fTitle;
        syst->fHasUpVariation = fHasUpVariation;
        syst->fHasDownVariation = fHasDownVariation;
        syst->fRegions = Common::ToVec(reg->fName);
        syst->fSamples = Common::ToVec(sampleName);
        syst->fCategory = fCategory;
        syst->fSubCategory = fSubCategory;
        syst->fSampleSmoothing = fSampleSmoothing;
        syst->fSymmetrisationType = fSymmetrisationType;
        syst->fSampleSmoothOption = fSampleSmoothingOption;
        syst->fSmoothType = fSmoothingType;
        syst->fOverallUp = fOverallUp;
        syst->fOverallDown = fOverallDown;
        syst->fCombineName = fCombineName;
        syst->fCombineType = fCombineType;
        syst->fDropShapeIn = fDropShapeIn;
        syst->fDropNormIn = fDropNorm;
        syst->fScaleUp = fScaleUp;
        syst->fScaleDown = fScaleDown;
        syst->fScaleUpRegions = fScaleUpRegions;
        syst->fScaleDownRegions = fScaleDownRegions;

        TRExFitter::SYSTMAP[syst->fName] = syst->fTitle;

        // Paths
        if (fHasUpVariation) {
            syst->fHistoPathsUp = Common::ToVec(name + "/UnfoldingHistograms");
            syst->fHistoFilesUp = Common::ToVec("FoldedHistograms");
            const std::string histoName = fName + "_Up/" + unfolding->fName + "_" + reg->fName + "_" + unfoldingSampleName + "_bin_" + std::to_string(ibin);
            syst->fHistoNamesUp = Common::ToVec(histoName);
            if (this->HasSubtract()) {
                const std::string histoNameSub = fName + "_subtract_Up/" + unfolding->fName + "_" + reg->fName + "_" + unfoldingSampleName + "_bin_" + std::to_string(ibin);
                syst->fHistoNamesUpSubtractSample = Common::ToVec(histoNameSub);
            }
        }
        if (fHasDownVariation) {
            syst->fHistoPathsDown = Common::ToVec(name + "/UnfoldingHistograms");
            syst->fHistoFilesDown = Common::ToVec("FoldedHistograms");
            const std::string histoName = fName + "_Down/" + unfolding->fName + "_" + reg->fName + "_" + unfoldingSampleName + "_bin_" + std::to_string(ibin);
            syst->fHistoNamesDown = Common::ToVec(histoName);
            if (this->HasSubtract()) {
                const std::string histoNameSub = fName + "_subtract_Down/" + unfolding->fName + "_" + reg->fName + "_" + unfoldingSampleName + "_bin_" + std::to_string(ibin);
                syst->fHistoNamesDownSubtractSample = Common::ToVec(histoNameSub);
            }
        }

        for(auto isample : samples) {
            if(!isample->fUseSystematics) continue;
            if(Common::FindInStringVector(syst->fSamples, isample->fName) < 0) continue;

            isample->AddSystematic(syst);
        }

        result.emplace_back(syst);
    }

    return result;
}

bool UnfoldingSystematic::HasSubtract() const {
    return !(fResponseMatrixPathsUpSubtractSample.empty() && fResponseMatrixPathsDownSubtractSample.empty() &&
            fResponseMatrixFilesUpSubtractSample.empty() && fResponseMatrixFilesDownSubtractSample.empty() &&
            fResponseMatrixNamesUpSubtractSample.empty() && fResponseMatrixNamesDownSubtractSample.empty() &&
            fMigrationPathsUpSubtractSample.empty() && fMigrationPathsDownSubtractSample.empty() &&
            fMigrationFilesUpSubtractSample.empty() && fMigrationFilesDownSubtractSample.empty() &&
            fMigrationNamesUpSubtractSample.empty() && fMigrationNamesDownSubtractSample.empty() &&
            fSelectionEffPathsUpSubtractSample.empty() && fSelectionEffPathsDownSubtractSample.empty() &&
            fSelectionEffFilesUpSubtractSample.empty() && fSelectionEffFilesDownSubtractSample.empty() &&
            fSelectionEffNamesUpSubtractSample.empty() && fSelectionEffNamesDownSubtractSample.empty() &&
            fAcceptancePathsUpSubtractSample.empty() && fAcceptancePathsDownSubtractSample.empty() &&
            fAcceptanceFilesUpSubtractSample.empty() && fAcceptanceFilesDownSubtractSample.empty() &&
            fAcceptanceNamesUpSubtractSample.empty() && fAcceptanceNamesDownSubtractSample.empty() &&
            fResponseMatrixFolderNamesUpSubtractSample.empty() && fResponseMatrixFolderNamesDownSubtractSample.empty() &&
            fMigrationFolderNamesUpSubtractSample.empty() && fMigrationFolderNamesDownSubtractSample.empty() &&
            fSelectionEffFolderNamesUpSubtractSample.empty() && fSelectionEffFolderNamesDownSubtractSample.empty() &&
            fAcceptanceFolderNamesUpSubtractSample.empty() && fAcceptanceFolderNamesDownSubtractSample.empty()
            );
}