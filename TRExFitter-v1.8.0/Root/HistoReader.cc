#include "TRExFitter/HistoReader.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/SampleHist.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/Systematic.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/TRExFit.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <algorithm>

HistoReader::HistoReader(TRExFit* fitter) :
    fFitter(fitter)
{
}

HistoReader::~HistoReader() {
}

void HistoReader::ReadHistograms(const std::vector<std::string>& binningsToCreate) {
    //
    // Loop on regions and samples
    //
    for(std::size_t i_ch = 0; i_ch < fFitter->fRegions.size(); ++i_ch) {
        LOG(INFO) << "  Region " << fFitter->fRegions[i_ch]->fName << " ...\n";
        //
        if(TRExFitter::SPLITHISTOFILES) fFitter->fFiles[i_ch]->cd();
        //
        if(fFitter->fRegions[i_ch]->fBinTransfo != "") fFitter->ComputeBinning(i_ch);
        auto itr = std::find(binningsToCreate.begin(), binningsToCreate.end(), fFitter->fRegions.at(i_ch)->fName);
        // first we must read the DATA samples
        ReadOneRegion(i_ch, true, false);

        // then we can read the other samples
        ReadOneRegion(i_ch, false, itr != binningsToCreate.end());
    }
}

std::unique_ptr<TH1> HistoReader::ReadSingleHistogram(const std::vector<std::string>& fullPaths,
                                                      const std::vector<std::string>& subtractPaths,
                                                      Systematic* syst,
                                                      int i_ch,
                                                      int i_smp,
                                                      bool isUp,
                                                      bool isMC) {

    std::unique_ptr<TH1> hsub(nullptr);
    // Need to get the hitogram that is used for subtracting
    if (!subtractPaths.empty()) {
        std::size_t i(0);
        for (const auto& ipath : subtractPaths) {
            std::unique_ptr<TH1> tmp = Common::HistFromFile(ipath);
            if (!tmp) {
                LOG(ERROR) << "Histo pointer is nullptr, cannot continue running the code\n";
                exit(EXIT_FAILURE);
            }
            //Pre-processing of histograms (rebinning, lumi scaling)
            if(fFitter->fRegions[i_ch]->fHistoBins.size() > 0){
                const std::string hname = tmp->GetName();
                std::unique_ptr<TH1> htmp_copy(static_cast<TH1*>(tmp->Rebin(fFitter->fRegions[i_ch]->fHistoNBinsRebin, "tmp_copy", &(fFitter->fRegions[i_ch]->fHistoBins[0]))));
                tmp.reset(htmp_copy.release());
                tmp->SetName(hname.c_str());
                if(TRExFitter::MERGEUNDEROVERFLOW) Common::MergeUnderOverFlow(tmp.get());
            }
            else if(fFitter->fRegions[i_ch]->fHistoNBinsRebin != -1) {
                tmp->Rebin(fFitter->fRegions[i_ch]->fHistoNBinsRebin);
            }

            if (isMC){
                if(fFitter->fSamples[i_smp]->fNormalizedByTheory){
                    tmp->Scale(fFitter->fLumi);
                }
            }

            if(fFitter->fSamples[i_smp]->fLumiScales.size() > i){
                 tmp->Scale(fFitter->fSamples[i_smp]->fLumiScales[i]);
            }
            else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1){
                tmp->Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
            }
            if (i == 0) {
                hsub.reset(static_cast<TH1*>(tmp->Clone()));
            } else {
                hsub->Add(tmp.get());
            }
            ++i;
        }

        if(syst->fReferenceSample != "") {
            this->PropagateReference(hsub.get(), syst, i_ch, i_smp);
        }
    }

    std::unique_ptr<TH1> result(nullptr);
    for(unsigned int i_path = 0; i_path < fullPaths.size(); ++i_path){
        std::unique_ptr<TH1> htmp = Common::HistFromFile( fullPaths.at(i_path) );
        if (!htmp) {
            if (!TRExFitter::HISTOCHECKCRASH) continue;
            LOG(ERROR) << "Histo pointer is nullptr, cannot continue running the code\n";
            exit(EXIT_FAILURE);
        }
        //Pre-processing of histograms (rebinning, lumi scaling)
        if(fFitter->fRegions[i_ch]->fHistoBins.size() > 0){
            const std::string hname = htmp->GetName();
            std::unique_ptr<TH1> tmp_copy(static_cast<TH1*>(htmp->Rebin(fFitter->fRegions[i_ch]->fHistoNBinsRebin, "tmp_copy", &(fFitter->fRegions[i_ch]->fHistoBins[0]))));
            htmp.reset(tmp_copy.release());
            htmp->SetName(hname.c_str());
            if(TRExFitter::MERGEUNDEROVERFLOW) Common::MergeUnderOverFlow(htmp.get());
        }
        else if(fFitter->fRegions[i_ch]->fHistoNBinsRebin != -1) {
            htmp->Rebin(fFitter->fRegions[i_ch]->fHistoNBinsRebin);
        }

        if (isMC){
            if(fFitter->fSamples[i_smp]->fNormalizedByTheory){
                htmp->Scale(fFitter->fLumi);
            }
        }

        if(fFitter->fSamples[i_smp]->fLumiScales.size()>i_path){
             htmp->Scale(fFitter->fSamples[i_smp]->fLumiScales[i_path]);
        }
        else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1){
            htmp->Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
        }

        if (isMC && syst != nullptr) {
            // obtain relative variation and apply it to proper sample
            // & try to keep also the same total relative variation
            if(syst->fReferenceSample != "") {
                this->PropagateReference(htmp.get(), syst, i_ch, i_smp);
            }
        }

        if(i_path == 0 || !result){
            if (!syst) { // is nominal
                result.reset(static_cast<TH1*>(htmp->Clone(Form("h_%s_%s",fFitter->fRegions[i_ch]->fName.c_str(),
                    fFitter->fSamples[i_smp]->fName.c_str()))));
            } else { // is syst
                if (isUp){ // up variation
                    result.reset(static_cast<TH1*>(htmp->Clone(Form("h_%s_%s_%sUp",fFitter->fRegions[i_ch]->fName.c_str(),
                        fFitter->fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()))));
                } else { // down variation
                    result.reset(static_cast<TH1*>(htmp->Clone(Form("h_%s_%s_%sDown",fFitter->fRegions[i_ch]->fName.c_str(),
                        fFitter->fSamples[i_smp]->fName.c_str(),syst->fStoredName.c_str()))));
                }
            }
        }
        else{
            result->Add(htmp.get());
        }
    }

    if (result) {
        if (hsub) {
            // final plot: systematic variation - variation on pseudodata + nominal
            const TH1* hnom = fFitter->fRegions[i_ch]->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();
            result->Add(hsub.get(), -1);
            result->Add(hnom);
        }
        result->SetDirectory(nullptr);
    }
    return result;
}

void HistoReader::ReadTRExProducedHistograms(const bool isUnfolding, const bool& checkHistogram) {
    std::string fileName("");
    std::string fileNameNom("");

    const bool singleOutputFile = !TRExFitter::SPLITHISTOFILES;
    std::string fileNameEnd = "_histos"+ fFitter->fSaveSuffix +".root";
    if(singleOutputFile){
        if(fFitter->fInputFolder!="") fileName = fFitter->fInputFolder           + fFitter->fInputName + fileNameEnd;
        else                          fileName = fFitter->fName + "/Histograms/" + fFitter->fInputName + fileNameEnd;

        // Bootstrap
        fileNameNom = fileName;
        if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0){

            if (fFitter->fBootstrapSyst!="") {
                fileName = Common::ReplaceString(fileName, fileNameEnd,Form("_histos%s_%s_%d.root", fFitter->fSaveSuffix.c_str(), fFitter->fBootstrapSyst.c_str(),fFitter->fBootstrapIdx));
                fileNameNom = fileName;
            } else {
                fileName = Common::ReplaceString(fileName, fileNameEnd,Form("_histos%s_%s_%d.root",fFitter->fSaveSuffix.c_str(), fFitter->fBootstrapSample.c_str(),fFitter->fBootstrapIdx));
                fileNameNom = fileName;
            }
            if(fFitter->fBootstrapNomHistos != "") fileNameNom = fFitter->fBootstrapNomHistos + fFitter->fInputName + fileNameEnd;
        }
        LOG(INFO) << "-----------------------------\n";
        LOG(INFO) << "Reading histograms from file " << fileName << " ...\n";
        if (fFitter->fBootstrapSyst!="" && fFitter->fBootstrapIdx>=0) {
            LOG(INFO) << "-----------------------------\n";
            LOG(INFO) << "Reading non-bootstrapped histograms from file " << fileNameNom << " ...\n";
        }
    }
    //
    std::vector< std::unique_ptr<TH2> > histPrun;
    std::unique_ptr<TFile> filePrun(nullptr);
    if( fFitter->fKeepPruning ){
        filePrun.reset(TFile::Open((fFitter->fName+"/Pruning.root").c_str()));
        if(!filePrun) fFitter->fKeepPruning = false;
    }
    //
    // when we multply/divide by or subtract/add other samples, need to add systematics on the other samples
    for(auto& isample : fFitter->fSamples) {
        if(!isample->fUseSystematics) continue;
        if(isample->fDivideBy!=""){
            std::shared_ptr<Sample> smp = fFitter->GetSample(isample->fDivideBy);
            for(const auto& isyst : smp->fSystematics) {
                const std::string systNPName = isyst->fNuisanceParameter;
                if(!isample->HasNuisanceParameter(systNPName)){
                    LOG(DEBUG) << " The sample " << isample->fName << " doesn't have natively NP " << systNPName << "\n";
                    LOG(DEBUG) << "                Inheriting it from " << smp->fName << "\n";
                    std::shared_ptr<Systematic> tmpsyst = std::make_shared<Systematic>(*isyst);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    isample->AddSystematic(tmpsyst);
                    for(const auto& jsyst : fFitter->fSystematics) {
                       if(jsyst->fName==tmpsyst->fName) {
                          if( Common::FindInStringVector(jsyst->fSamples,
                                                 isample->fName)<0 ) jsyst->fSamples.push_back(isample->fName);
                       }
                    }
                }
            }
        }
        if(isample->fMultiplyBy!=""){
            std::shared_ptr<Sample> smp = fFitter->GetSample(isample->fMultiplyBy);
            for(const auto& isyst : smp->fSystematics) {
                const std::string systNPName = isyst->fNuisanceParameter;
                if(!isample->HasNuisanceParameter(systNPName)){
                    std::shared_ptr<Systematic> tmpsyst = std::make_shared<Systematic>(*isyst);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    LOG(DEBUG) << " The sample " << isample->fName << " doesn't have natively NP " << systNPName << "\n";
                    LOG(DEBUG) << "                Inheriting it from " << smp->fName << "\n";
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    isample->AddSystematic(tmpsyst);
                    for(const auto& jsyst : fFitter->fSystematics) {
                       if(jsyst->fName==tmpsyst->fName) {
                          if( Common::FindInStringVector(jsyst->fSamples,
                                                 isample->fName)<0 ) jsyst->fSamples.push_back(isample->fName);
                       }
                    }
                }
            }
        }
        for(const auto& sample : isample->fSubtractSamples){
            std::shared_ptr<Sample> smp = fFitter->GetSample(sample);
            for(const auto& isyst : smp->fSystematics) {
                const std::string systNPName = isyst->fNuisanceParameter;
                if(!isample->HasNuisanceParameter(systNPName)){
                    LOG(DEBUG) << " The sample " << isample->fName << " doesn't have natively NP " << systNPName << "\n";
                    LOG(DEBUG) << "                Inheriting it from " << smp->fName << "\n";
                    std::shared_ptr<Systematic> tmpsyst = std::make_shared<Systematic>(*isyst);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    isample->AddSystematic(tmpsyst);
                    for(const auto& jsyst : fFitter->fSystematics) {
                       if(jsyst->fName==tmpsyst->fName) {
                          if(Common::FindInStringVector(jsyst->fSamples,
                                                isample->fName)<0 ) jsyst->fSamples.push_back(isample->fName);
                       }
                    }
                }
            }
        }
        for(const auto& sample : isample->fAddSamples){
            std::shared_ptr<Sample> smp = fFitter->GetSample(sample);
            for(const auto& isyst : smp->fSystematics) {
                const std::string systNPName = isyst->fNuisanceParameter;
                if(!isample->HasNuisanceParameter(systNPName)){
                    LOG(DEBUG) << " The sample " << isample->fName << " doesn't have natively NP " << systNPName << "\n";
                    LOG(DEBUG) << "                Inheriting it from " << smp->fName << "\n";
                    std::shared_ptr<Systematic> tmpsyst = std::make_shared<Systematic>(*isyst);
                    tmpsyst->fName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    tmpsyst->fStoredName = systNPName; // want to inherit the triggering systematic, not the derived ones
                    if (tmpsyst->fType == Systematic::OVERALL ) {
                          tmpsyst->fType = Systematic::HISTO; // even if it was overall for "inheritors", that's not guaranteed for the "inheritand"
                          tmpsyst->fIsNormOnly = false;
                    }
                    isample->AddSystematic(tmpsyst);
                    for(const auto& jsyst : fFitter->fSystematics) {
                       if(jsyst->fName==tmpsyst->fName) {
                          if(Common::FindInStringVector(jsyst->fSamples,
                                                isample->fName)<0 ) jsyst->fSamples.push_back(isample->fName);
                       }
                    }
                }
            }
        }
    }
    //
    // Syst for morphing samples inherited from nominal sample
    if (fFitter->fPropagateSystsForMorphing){
        for(const auto& par : fFitter->fMorphParams){
            double nominalValue = 0.;
            for(const auto& norm : fFitter->fNormFactors){
                if(norm->fName == par) nominalValue = norm->GetNominal();
            }
            std::shared_ptr<Sample> smpNominal = nullptr;
            for(const auto& smp : fFitter->fSamples){
                if(!smp->fIsMorph[par]) continue;
                if(smp->fMorphValue[par] == nominalValue){ // FIXME: eventually add something to flag a sample as nominal for morphing
                    smpNominal = smp;
                    break;
                }
            }
            for(const auto& smp : fFitter->fSamples){
                if(!smp->fIsMorph[par]) continue;
                if(smp == smpNominal) continue;
                for(const auto& syst : smpNominal->fSystematics){
                    smp->AddSystematic(syst);
                }
            }
        }
    }
    //
    // fSystFromSample
    for(const auto& smp : fFitter->fSamples){
        if(smp->fSystFromSample!=""){
            std::shared_ptr<Sample> smpReference = nullptr;
            for(const auto& smp2 : fFitter->fSamples){
                if(smp2->fName == smp->fName) continue;
                if(smp2->fName == smp->fSystFromSample){
                    smpReference = smp2;
                    break;
                }
            }
            if(smpReference!=nullptr){
                for(const auto& syst : smpReference->fSystematics){
                    smp->AddSystematic(syst);
                }
            }
        }
    }
    //
    for(std::size_t i_ch = 0; i_ch<fFitter->fRegions.size(); ++i_ch) {
        if(fFitter->fKeepPruning){
            histPrun.emplace_back(static_cast<TH2*>(filePrun->Get(Form("h_prun_%s_toSave", fFitter->fRegions[i_ch]->fName.c_str()))));
        }
        const std::string regionName = fFitter->fRegions[i_ch]->fName;
        LOG(DEBUG) << "  Reading region " << regionName << "\n";

        const std::vector<int>& droppedBins = fFitter->fRegions[i_ch]->GetAutomaticDropBins() ?
                                              fFitter->fRegions[i_ch]->fComputedBlindedBins : fFitter->fRegions[i_ch]->fDropBins;

        bool dropBinHistoNeeded = droppedBins.size() > 0;

        //
        if(!singleOutputFile){
            if(fFitter->fInputFolder!="") fileName = fFitter->fInputFolder           + fFitter->fInputName + "_" + regionName + fileNameEnd;
            else                          fileName = fFitter->fName + "/Histograms/" + fFitter->fInputName + "_" + regionName + fileNameEnd;

            // Bootstrap
            fileNameNom = fileName;
            if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0){
                if (fFitter->fBootstrapSyst!="") {
                    fileName = Common::ReplaceString(fileName,fileNameEnd,Form("_histos%s_%s_%d.root", fFitter->fSaveSuffix.c_str(), fFitter->fBootstrapSyst.c_str(),fFitter->fBootstrapIdx));
                    fileNameNom = fileName;
                } else {
                    fileName = Common::ReplaceString(fileName,fileNameEnd,Form("_histos%s_%s_%d.root", fFitter->fSaveSuffix.c_str(), fFitter->fBootstrapSample.c_str(),fFitter->fBootstrapIdx));
                    fileNameNom = fileName;
                }
                if(fFitter->fBootstrapNomHistos != "") fileNameNom = fFitter->fBootstrapNomHistos + fFitter->fInputName + "_" + regionName  + fileNameEnd;
            }
            LOG(INFO) << "-----------------------------\n";
            LOG(INFO) << "Reading histograms from file " << fileName << " ...\n";
            if (fFitter->fBootstrapSyst!="" && fFitter->fBootstrapIdx>=0) {
                LOG(INFO) << "-----------------------------\n";
                LOG(INFO) << "Reading non-bootstrapped histograms from file " << fileNameNom << " ...\n";
            }
        }
        //
        for(const auto& isample : fFitter->fSamples) {
            //
            // eventually skip sample / region combination
            //
            if(Common::FindInStringVector(isample->fRegions,regionName)<0 && isample->fName.find("customAsimov_")==std::string::npos ) continue;
            //
            const std::string sampleName = isample->fName;
            Common::FolderStructure structure;
            structure.region = regionName;
            structure.sample = sampleName;
            structure.systematics = "nominal";
            LOG(DEBUG) << "    Reading sample " << sampleName << "\n";

            if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0 && fFitter->fBootstrapSample == sampleName ){
                fFitter->fRegions[i_ch]->SetSampleHist(isample.get(),regionName+"_"+sampleName,fileName, dropBinHistoNeeded, checkHistogram, &structure);
            }
            else{
                fFitter->fRegions[i_ch]->SetSampleHist(isample.get(),regionName+"_"+sampleName,fileNameNom, dropBinHistoNeeded, checkHistogram, &structure);
            }



            std::shared_ptr<SampleHist> sh = fFitter->fRegions[i_ch]->GetSampleHist(sampleName);
            if(!sh) continue;

            // separate gammas -> Add systematic
            if(isample->fSeparateGammas){
                std::string systName = "stat_"+isample->fName;
                std::string systStoredName = systName;
                structure.systematics = systName;
                LOG(DEBUG) << "Adding separate gammas as SHAPE systematic " << systName << "\n";
                std::shared_ptr<SystematicHist> syh_tmp = sh->AddHistoSyst(systName,
                                                                           systStoredName,
                                                                           Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                                           Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                                                           dropBinHistoNeeded,
                                                                           checkHistogram,
                                                                           0,
                                                                           &structure
                                                                          );
                if(!syh_tmp){
                    LOG(WARNING) << "No histogram found for separate gamma, but may be you will create it right now.\n";
                }
                else{
                    std::shared_ptr<Systematic> gamma = nullptr;
                    if(Common::FindInStringVector(fFitter->fSystematicNames,systName)>=0) gamma = fFitter->fSystematics[Common::FindInStringVector(fFitter->fSystematicNames,systName)];  //GetSystematic(systName);
                    if(gamma==nullptr) gamma = fFitter->NewSystematic(systName);
                    LOG(DEBUG) << "Adding separate gammas as SHAPE systematic " << systName << "\n";
                    gamma->fType = Systematic::SHAPE;
                    gamma->fRegions.clear();
                    gamma->fRegions.push_back(fFitter->fRegions[i_ch]->fName);
                    syh_tmp->fSystematic = gamma;
                    gamma->fNuisanceParameter = gamma->fName;
                    TRExFitter::NPMAP[gamma->fName] = gamma->fNuisanceParameter;
                    gamma->fCategory = (isUnfolding && isample->fType == Sample::SampleType::SIGNAL) ? "gamma_response" : "gamma_" + isample->fName;
                    gamma->fSubCategory = (isUnfolding && isample->fType == Sample::SampleType::SIGNAL) ? "gamma_response" : "gamma_" + isample->fName;
                }
            }
            //
            // norm factors
            for(const auto& inorm : isample->fNormFactors) {
                //
                // eventually skip norm factor / region combination
                if(inorm->fRegions.size()>0 && Common::FindInStringVector(inorm->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if(inorm->fExclude.size()>0 && Common::FindInStringVector(inorm->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                const std::string normName = inorm->fName;
                LOG(DEBUG) << "      Reading norm " << normName << "\n";
                // norm only
                sh->AddNormFactor(inorm);
            }
            //
            // shape factors
            for(const auto& ishape : isample->fShapeFactors) {
                //
                // eventually skip shape factor / region combination
                if(ishape->fRegions.size()>0 && Common::FindInStringVector(ishape->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if(ishape->fExclude.size()>0 && Common::FindInStringVector(ishape->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                //
                const std::string shapeName = ishape->fName;
                LOG(DEBUG) << "      Reading shape " << shapeName << "\n";
                // shape only
                sh->AddShapeFactor(ishape);
            }
            //
            // systematics
            for(std::size_t i_syst = 0; i_syst< isample->fSystematics.size(); ++i_syst) {
                Common::FolderStructure systStructure;
                systStructure.region = regionName;
                systStructure.sample = sampleName;
                // here we check if systematic has been a clone from a shape-acc systematic
                // and if this is the case, remove the _shape and _acc from the name before looking for it in hist files
                std::string temp_systname = isample->fSystematics[i_syst]->fName;
                if(isample->fSystematics[i_syst]->fIsShapeAccDecorr){

                    if(isample->fSystematics[i_syst]->fIsNormOnly)          temp_systname = Common::ReplaceString(isample->fSystematics[i_syst]->fName, "_Acc", "");
                    else if(isample->fSystematics[i_syst]->fIsShapeOnly)    temp_systname = Common::ReplaceString(isample->fSystematics[i_syst]->fName, "_Shape", "");
                }
                if(isample->fSystematics[i_syst]->fIsRegionDecorr){
                    temp_systname = Common::ReplaceString(isample->fSystematics[i_syst]->fName, "_" + fFitter->fRegions[i_ch]->fName, "");
                }

                if(isample->fSystematics[i_syst]->fIsSampleDecorr){
                    temp_systname = Common::ReplaceString(isample->fSystematics[i_syst]->fName, "_" + isample->fName, "");
                }
                systStructure.systematics = temp_systname;
                //
                // eventually skip systematic / region combination
                if(isample->fSystematics[i_syst]->fRegions.size()>0 && Common::FindInStringVector(isample->fSystematics[i_syst]->fRegions,fFitter->fRegions[i_ch]->fName)<0  ) continue;
                if(isample->fSystematics[i_syst]->fExclude.size()>0 && Common::FindInStringVector(isample->fSystematics[i_syst]->fExclude,fFitter->fRegions[i_ch]->fName)>=0 ) continue;
                if(isample->fSystematics[i_syst]->fExcludeRegionSample.size()>0 && Common::FindInStringVectorOfVectors(isample->fSystematics[i_syst]->fExcludeRegionSample,fFitter->fRegions[i_ch]->fName, isample->fName)>=0 ) continue;
                //
                const std::string systName       = isample->fSystematics[i_syst]->fName;
                const std::string systStoredName = isample->fSystematics[i_syst]->fStoredName; // if no StoredName specified in the config, this should be == fName
                //
                // eventually skip systematics if pruned
                int binContent(0);
                if( fFitter->fKeepPruning && histPrun[i_ch] != nullptr ){
                    const int xbin = histPrun[i_ch]->GetXaxis()->FindBin( sampleName.c_str() ); // sample
                    const int ybin = histPrun[i_ch]->GetYaxis()->FindBin( systName.c_str() ); // syst
                    const int bin = histPrun[i_ch]->GetBin(xbin,ybin);
                    binContent = histPrun[i_ch]->GetBinContent(bin);
                    if( binContent <= -4 || binContent == -1 || binContent >= 3 ){
                        LOG(DEBUG) << "SKIPPING systematic " << systName << "\n";
                        continue;
                    }
                }
                LOG(DEBUG) << "      Reading syst " << systName << "\n";
                // norm only
                std::shared_ptr<SystematicHist> syh(nullptr);
                if(isample->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    if( fFitter->fKeepPruning ){
                        if( binContent == -2 || binContent == 2 ) continue;
                    }
                    syh = sh->AddOverallSyst(systName,
                                             systStoredName,
                                             isample->fSystematics[i_syst]->fOverallUp,
                                             isample->fSystematics[i_syst]->fOverallDown);
                }
                // histo syst
                else{
                    int pruned = 0;
                    if(fFitter->fKeepPruning){
                        if(binContent==1 || binContent==-2) pruned = 1;
                        if(binContent==2 || binContent==-3) pruned = 2;
                    }

                    if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0 && ( fFitter->fBootstrapSyst == isample->fSystematics[i_syst]->fNuisanceParameter ||  fFitter->fBootstrapSyst == systName ) ){
                        syh = sh->AddHistoSyst(systName,
                                               systStoredName,
                                               Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                               Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileName,
                                               dropBinHistoNeeded,
                                               checkHistogram,
                                               pruned,
                                               &systStructure
                                              );
                    } else {
                        syh = sh->AddHistoSyst(systName,
                                                systStoredName,
                                                Form("%s_%s_%s_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileNameNom,
                                                Form("%s_%s_%s_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str()), fileNameNom,
                                                dropBinHistoNeeded,
                                                checkHistogram,
                                                pruned,
                                                &systStructure
                                                );
                    }

                    if(!syh) {
                        if (!pruned) LOG(WARNING) << "No syst histo found for syst " << systName << ", sample " << sampleName << ", region " << regionName << "\n";
                        continue;
                    }
                }
                // for both
                syh->fSystematic = isample->fSystematics[i_syst];
                syh->fHistoNameShapeUp   = Form("%s_%s_%s_Shape_Up",  regionName.c_str(),sampleName.c_str(),systStoredName.c_str());
                syh->fHistoNameShapeDown = Form("%s_%s_%s_Shape_Down",regionName.c_str(),sampleName.c_str(),systStoredName.c_str());

                if(fFitter->fBootstrap!="" && fFitter->fBootstrapIdx>=0 && ( fFitter->fBootstrapSyst == isample->fSystematics[i_syst]->fNuisanceParameter ||  fFitter->fBootstrapSyst == systName ) ){
                    syh->fFileNameShapeUp    = fileName;
                    syh->fFileNameShapeDown  = fileName;
                } else {
                    syh->fFileNameShapeUp    = fileNameNom;
                    syh->fFileNameShapeDown  = fileNameNom;
                }



                syh->fScaleUp = isample->fSystematics[i_syst]->fScaleUp;
                if(isample->fSystematics[i_syst]->fScaleUpRegions.size()!=0){
                    if(isample->fSystematics[i_syst]->fScaleUpRegions[regionName]!=0){
                        syh->fScaleUp *= isample->fSystematics[i_syst]->fScaleUpRegions[regionName];
                    }
                }
                syh->fScaleDown = isample->fSystematics[i_syst]->fScaleDown;
                if(isample->fSystematics[i_syst]->fScaleDownRegions.size()!=0){
                    if(isample->fSystematics[i_syst]->fScaleDownRegions[regionName]!=0){
                        syh->fScaleDown *= isample->fSystematics[i_syst]->fScaleDownRegions[regionName];
                    }
                }
                //
                if(isample->fSystematics[i_syst]->fType == Systematic::OVERALL){
                    syh->fNormUp   *= syh->fScaleUp;
                    syh->fNormDown *= syh->fScaleDown;
                }
            }
        }
    }

    if (filePrun) {
        filePrun->Close();
    }
}

void HistoReader::ReadOneRegion(const int i_ch, const bool is_data, const bool recreateBinnings) {
    int nbins(0);
    for(const auto& ismp : fFitter->fSamples) {
        if (is_data && ismp->fType!=Sample::SampleType::DATA) continue;
        if (!is_data && ismp->fType==Sample::SampleType::DATA) continue;

        auto itr = std::find(fFitter->fSamples.begin(), fFitter->fSamples.end(), ismp);
        const std::size_t i_smp = std::distance(fFitter->fSamples.begin(), itr);
        LOG(DEBUG) << "  Reading " << ismp->fName << "\n";
        //
        // eventually skip sample / region combination
        //
        if(Common::FindInStringVector(ismp->fRegions,fFitter->fRegions[i_ch]->fName)<0) continue;
        //
        // read nominal
        //
        std::vector<std::string> fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch].get(),ismp.get(), nullptr, true, ismp->fIsFolded);

        std::unique_ptr<TH1> h = ReadSingleHistogram(fullPaths, {}, nullptr, i_ch, i_smp, true, !is_data); // is MC

        if (!h) {
            if (!TRExFitter::HISTOCHECKCRASH) continue;
            LOG(ERROR) << "Histo pointer is nullptr, cannot continue running the code\n";
            exit(EXIT_FAILURE);
        }

        // Save the original histogram
        TH1* h_orig = static_cast<TH1*>(h->Clone( Form("%s_orig",h->GetName())));

        // Importing the histogram in TRExFitter
        std::shared_ptr<SampleHist> sh = fFitter->fRegions[i_ch]->SetSampleHist(ismp.get(), h.get());
        sh->fHist_orig.reset(h_orig);
        sh->fHist_orig->SetName(Form("%s_orig",sh->fHist->GetName())); // fix the name
        sh->fHist_orig->SetDirectory(nullptr);


        // end here no systematics allowed (e.g. generally for GHOST samples)
        if (!(ismp->fUseSystematics) && !is_data) continue;

        if (!is_data) {
            ReadNormShape(sh.get(), i_ch, ismp.get());
        }

        nbins = sh->fHist_orig->GetNbinsX();

        for(const auto& isyst : ismp->fSystematics) {

            Systematic *syst = isyst.get();
            if (is_data) continue;

            // eventually skip systematic / region combination
            if(syst->fRegions.size()>0 &&
               Common::FindInStringVector(syst->fRegions,fFitter->fRegions[i_ch]->fName)<0 ) continue;
            if(syst->fExclude.size()>0 &&
               Common::FindInStringVector(syst->fExclude,fFitter->fRegions[i_ch]->fName)>=0) continue;
            if(syst->fExcludeRegionSample.size()>0 &&
                Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample,
                                                    fFitter->fRegions[i_ch]->fName,
                                                    ismp->fName)>=0) continue;

            LOG(DEBUG) << "Adding syst " << syst->fName << "\n";

            if (!is_data) {
                if (SetSystematics(i_ch, ismp.get(), isyst)) continue;
            }
            std::unique_ptr<TH1> hUp = GetSystHisto(i_ch, i_smp, syst, is_data, true);
            std::unique_ptr<TH1> hDown = GetSystHisto(i_ch, i_smp, syst, is_data, false);
            std::shared_ptr<SystematicHist> syh = sh->AddHistoSyst(isyst->fName,
                                                                   isyst->fStoredName,
                                                                   hUp.get(),
                                                                   hDown.get());
            syh->fSystematic = isyst;
            syh->fScaleUp = isyst->fScaleUp;
            if(isyst->fScaleUpRegions.size()!=0) {
                if(isyst->fScaleUpRegions[fFitter->fRegions[i_ch]->fName]!=0) {
                    syh->fScaleUp *= isyst->fScaleUpRegions[fFitter->fRegions[i_ch]->fName];
                }
            }
            syh->fScaleDown = isyst->fScaleDown;
            if(isyst->fScaleDownRegions.size()!=0) {
                if(isyst->fScaleDownRegions[fFitter->fRegions[i_ch]->fName]!=0) {
                    syh->fScaleDown *= isyst->fScaleDownRegions[fFitter->fRegions[i_ch]->fName];
                }
            }
        }
    }
    if (!is_data) fFitter->fRegions.at(i_ch)->SetAndAddNbins(nbins, recreateBinnings);
}

void HistoReader::ReadNormShape(SampleHist* sh,
                                const int i_ch,
                                const Sample* smp) {
    // read norm factors
    for(const auto& inorm : smp->fNormFactors) {
        // eventually skip systematic / region combination
        if(inorm->fRegions.size()>0 && Common::FindInStringVector(inorm->fRegions,fFitter->fRegions[i_ch]->fName)<0) continue;
        if(inorm->fExclude.size()>0 && Common::FindInStringVector(inorm->fExclude,fFitter->fRegions[i_ch]->fName)>=0) continue;

        LOG(DEBUG) << "Adding norm " << inorm->fName << "\n";

        sh->AddNormFactor(inorm);
    }

    // read shape factors
    for(const auto& ishape : smp->fShapeFactors) {

        // eventually skip systematic / region combination
        if(ishape->fRegions.size()>0 && Common::FindInStringVector(ishape->fRegions,fFitter->fRegions[i_ch]->fName)<0) continue;
        if(ishape->fExclude.size()>0 && Common::FindInStringVector(ishape->fExclude,fFitter->fRegions[i_ch]->fName)>=0) continue;

        LOG(DEBUG) << "Adding shape " << ishape->fName << "\n";

        sh->AddShapeFactor(ishape);
    }
}

bool HistoReader::SetSystematics(const int i_ch,
                                 const Sample* ismp,
                                 std::shared_ptr<Systematic> syst) {

    const Region* reg = fFitter->fRegions[i_ch].get();

    // if Overall only ...
    if(syst->fType==Systematic::OVERALL) {
        std::shared_ptr<SystematicHist> syh = reg->GetSampleHist(ismp->fName)->AddOverallSyst(syst->fName,
                                                                                              syst->fStoredName,
                                                                                              syst->fOverallUp,
                                                                                              syst->fOverallDown);
        syh->fSystematic = syst;
        syh->fScaleUp = syst->fScaleUp;
        if(syst->fScaleUpRegions.size()!=0) {
            if(syst->fScaleUpRegions[reg->fName]!=0) {
                syh->fScaleUp *= syst->fScaleUpRegions[reg->fName];
            }
        }
        syh->fScaleDown = syst->fScaleDown;
        if(syst->fScaleDownRegions.size()!=0) {
            if(syst->fScaleDownRegions[reg->fName]!=0) {
                syh->fScaleDown *= syst->fScaleDownRegions[reg->fName];
            }
        }
        return true;
    }
    return false;
}

std::unique_ptr<TH1> HistoReader::GetSystHisto(const int i_ch,
                                               const std::size_t i_smp,
                                               Systematic* syst,
                                               bool is_data,
                                               bool is_up) {

    std::unique_ptr<TH1> result(nullptr);

    if ((syst->fHasDownVariation && !is_up) || (syst->fHasUpVariation && is_up)) {
        const auto& fullPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch].get(),
                                                            fFitter->fSamples[i_smp].get(),
                                                            syst,
                                                            is_up,
                                                            fFitter->fSamples[i_smp]->fIsFolded);

        bool hasSubtract(false);
        if (is_up) {
            hasSubtract = !syst->fHistoPathsUpSubtractSample.empty()   || !syst->fHistoFilesUpSubtractSample.empty()   || !syst->fHistoNamesUpSubtractSample.empty() || syst->fHistoNameSufUpSubtractSample!="" || !syst->fHistoFolderSubtractNamesUp.empty() || !syst->fHistoFileSufUpSubtractSample.empty();
        } else {
            hasSubtract = !syst->fHistoPathsDownSubtractSample.empty() || !syst->fHistoFilesDownSubtractSample.empty() || !syst->fHistoNamesDownSubtractSample.empty() || syst->fHistoNameSufDownSubtractSample!="" || !syst->fHistoFolderSubtractNamesDown.empty() || !syst->fHistoFileSufDownSubtractSample.empty();
        }

        std::vector<std::string> subtractPaths{};
        if (hasSubtract) {
            subtractPaths = fFitter->FullHistogramPaths(fFitter->fRegions[i_ch].get(),
                                                        fFitter->fSamples[i_smp].get(),
                                                        syst,
                                                        is_up,
                                                        fFitter->fSamples[i_smp]->fIsFolded,
                                                        true);
        }

        result = ReadSingleHistogram(fullPaths,
                                     subtractPaths,
                                     syst,
                                     i_ch,
                                     i_smp,
                                     is_up,
                                     !is_data);
    }

    if(!result) {
        result.reset(static_cast<TH1*>(fFitter->fRegions[i_ch]->GetSampleHist(fFitter->fSamples[i_smp]->fName)->fHist->Clone()));
    }

    result->SetDirectory(nullptr);
    return result;
}


void HistoReader::PropagateReference(TH1* h, Systematic* syst, const int i_ch, const std::size_t i_smp) const {
    if (fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample) == nullptr){
        LOG(ERROR) << "Reference sample: " << syst->fReferenceSample << " does not exist for region: " << fFitter->fRegions[i_ch]->fName << ". Please check this!\n";
        exit(EXIT_FAILURE);
    }
    TH1* href = fFitter->fRegions[i_ch]->GetSampleHist(syst->fReferenceSample)->fHist.get();
    const TH1* hnom = fFitter->fRegions[i_ch]->GetSampleHist(fFitter->fSamples[i_smp]->fName )->fHist.get();
    // Protection added: fix empty bins before starting to divide and multiply
    for (int i_bin = 0; i_bin < href->GetNbinsX() + 2; i_bin++) {
        if (href->GetBinContent(i_bin) <= 1e-6) {
            href->SetBinContent(i_bin,1e-6);
        }
    }
    for (int i_bin=0; i_bin < h->GetNbinsX() + 2; i_bin++) {
        if (h->GetBinContent(i_bin) <= 1e-6) {
            h->SetBinContent(i_bin,1e-6);
        }
    }
    auto hrefTmp = Common::GetHistCopyNoError(href);
    auto hnomTmp = Common::GetHistCopyNoError(hnom);
    h->Divide(  hrefTmp.get());
    h->Multiply(hnomTmp.get());
}
