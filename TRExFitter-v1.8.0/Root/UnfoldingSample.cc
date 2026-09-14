#include "TRExFitter/UnfoldingSample.h"

#include "TRExFitter/Logger.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/Unfolding.h"

UnfoldingSample::UnfoldingSample() :
    fUnfoldingName(""),
    fName(""),
    fTitle(""),
    fFillColor(0),
    fLineColor(0),
    fHasResponse(false),
    fHasAcceptance(false),
    fType(UnfoldingSample::TYPE::STANDARD),
    fGammas(UnfoldingSample::GAMMAS::SEPARATED),
    fNoPruning(false)
{
}

std::vector<std::shared_ptr<Sample> > UnfoldingSample::ConvertToSample(const Region* reg,
                                                                       const Unfolding* unfolding,
                                                                       const std::string& name) const {

    std::vector<std::shared_ptr<Sample> > result;
    for (int ibin = 0; ibin < unfolding->fNumberUnfoldingTruthBins; ++ibin) {
        const std::string sampleName = unfolding->fName + "_" + reg->fName + "_Truth_bin_" + std::to_string(ibin+1);

        std::shared_ptr<Sample> sample = std::make_shared<Sample>(sampleName, Sample::SampleType::SIGNAL);
        sample->SetTitle(fTitle);
        sample->SetFillColor(fFillColor);
        sample->SetLineColor(fLineColor);
        sample->fRegions = Common::ToVec(reg->fName);
        sample->fNoPruning = fNoPruning;

        // paths
        sample->fHistoPaths = Common::ToVec(name + "/UnfoldingHistograms");
        sample->fHistoFiles = Common::ToVec("FoldedHistograms");
        const std::string histoName = "nominal/" + unfolding->fName + "_" + reg->fName + "_" + fName + "_bin_" + std::to_string(ibin);
        sample->fHistoNames = Common::ToVec(histoName);
        sample->fIsFolded = true;
        sample->fUnfoldingName = unfolding->fName;
        sample->fUnfoldingSampleName = fName;

        if (fGammas == UnfoldingSample::GAMMAS::SEPARATED) {
            sample->fSeparateGammas = true;
            sample->fUseMCStat = false;
        } else if (fGammas == UnfoldingSample::GAMMAS::DISABLED) {
            sample->fUseMCStat = false;
        } else {
            LOG(ERROR) << "Unknown GAMMAS type\n";
            exit(EXIT_FAILURE);
        }

        result.emplace_back(sample);
    }

    return result;
}
