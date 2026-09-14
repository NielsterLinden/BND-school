#include "TRExFitter/NtupleReader.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/TRExFit.h"

#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"

NtupleReader::NtupleReader(TRExFit* fitter) :
    fFitter(fitter)
{
}

NtupleReader::~NtupleReader() {
}

TH1D* NtupleReader::GetHistogram(Region* region, Sample* sample, Systematic* syst, const bool isUpSyst, const bool isSubtract){


    std::string variable                   = fFitter->Variable(       region, sample);
    std::string fullSelection              = fFitter->FullSelection(  region, sample);
    std::string fullMCweight               = fFitter->FullWeight(     region, sample, syst, isUpSyst);
    std::vector<std::string> fullPaths     = fFitter->FullNtuplePaths(region, sample, syst, isUpSyst);

    /*
        Handle reading of friend trees/ ntuples
    */
    std::vector<std::string> friendPaths;
    if (sample->fUseFriend){
        friendPaths   = fFitter->FullNtuplePaths(region, sample, syst, isUpSyst, false,true); // isSubtract=False, isFriend=True
    }

    /*
        Handle reading of Subtraction histograms (replace sample's nominal paths with subtracted paths)
    */
    if (isSubtract){
        fullPaths =  fFitter->FullNtuplePaths(region,sample,syst,isUpSyst,true,false); // isSubtract=False, isFriend=False
        if (sample->fUseFriend){
            friendPaths = fFitter->FullNtuplePaths(region,sample,syst,isUpSyst,true,true); // isSubtract=True, isFriend=True
        }
    }

    /*
        Safety, hen there is no firendPaths ensure friendPaths vector
        as long as the FullPaths vector, but filled with empties
    */
    if(friendPaths.size() == 0) {
        for (unsigned int i=0; i < fullPaths.size(); i++){
            friendPaths.push_back("");
        }
    }

    std::string forWhat("for:");

    forWhat += " Region: " + region->fName;
    forWhat += " Sample: " + sample->fName;
    if (syst != nullptr){
        forWhat += " Systematic: " + syst->fName;
        if (isUpSyst) forWhat += " Up variation";
    }
    if (isSubtract){
        forWhat += " in Subtraction paths";
    }
    if (fullPaths.size() != friendPaths.size()) {
        LOG(ERROR) << "Sizes of paths and friend paths do not match for " << forWhat << "\n";
        exit(EXIT_FAILURE);
    }

    /*
        Read the Histograms from all paths and combine them
    */
    TH1D* histogram = nullptr;
    std::string histoName;
    if (syst == nullptr)    histoName = "h_"+region->fName+"_"+sample->fName;
    else{
        if (isUpSyst)   histoName = "h_"+region->fName+"_"+sample->fName+"_"+syst->fName+"_Up";
        else    histoName = "h_"+region->fName+"_"+sample->fName+"_"+syst->fName+"_Down";
    }

    for(unsigned int i_path = 0; i_path < fullPaths.size(); i_path++){
        TH1D* hpath = nullptr;

        if(region->fHistoBins.size() > 0){
            hpath = Common::HistFromNtupleBinArr(   fullPaths[i_path],
                                                    friendPaths[i_path],
                                                    variable,
                                                    region->fHistoNBinsRebin,
                                                    &(region->fHistoBins[0]),
                                                    fullSelection,
                                                    fullMCweight,
                                                    fFitter->fAddAliases,
                                                    fFitter->fDebugNev);
        }
        else{
            hpath = Common::HistFromNtuple( fullPaths[i_path],
                                            friendPaths[i_path],
                                            variable,
                                            region->GetNbins(),
                                            region->fXmin,
                                            region->fXmax,
                                            fullSelection,
                                            fullMCweight,
                                            fFitter->fAddAliases,
                                            fFitter->fDebugNev);
            //Pre-processing of histograms (rebinning, lumi scaling)
            if(region->fHistoNBinsRebin != -1){
                hpath->Rebin(region->fHistoNBinsRebin);
            }
        }

        /* Scaling of MC if NormalizedByTheory */
        if( sample->fNormalizedByTheory &&
            sample->fType!=Sample::SampleType::DATA
            ){
                hpath -> Scale(fFitter->fLumi);
            }

        /* Apply luminosity scalings */
        if (sample->fLumiScales.size() > i_path){
            hpath -> Scale(sample->fLumiScales[i_path]);
        }
        else if (sample->fLumiScales.size()==1){
            hpath -> Scale(sample->fLumiScales[0]);
        }

        /* Add Up! */
        if(i_path==0)   histogram = static_cast<TH1D*>(hpath->Clone(histoName.c_str()));
        else    histogram->Add(hpath);

        /* Clean up histogram from the individual path */
        delete hpath;
    } /* End of loop over paths */

    if (fFitter->fEFTConfig.GetSplitSamplesPerBin()) {

        // Leave here if histogram is null, no need to do EFT manipulation
        if (histogram == nullptr){
            LOG(WARNING) << "Report this. Could not build histogram for " << forWhat << "\n";
            if (syst == nullptr)    return histogram;
            else    histogram = static_cast<TH1D*>(region->GetSampleHist(sample->fName )->fHist.get());
        }

        // EFT Manipulation - Need to find EFT samples and zero out all bins except the one of interest
        if( sample->fType==Sample::SampleType::SIGNAL &&
            sample->fEFTSMReference != "NONE" &&
            sample->fEFTSMReference != ""){

            LOG(DEBUG) << "Found EFT SM Reference " << sample->fName << "\n";

            for(int ib=0; ib<region->GetNbins(); ++ib){
                std::string ending = TString::Format("_%s_bin%d", region->fName.c_str(),ib).Data();

                LOG(DEBUG) << Form("Checking %s",ending.c_str()) << "\n";

                if (    sample->fName.length() >= ending.length() &&
                        sample->fName.compare(sample->fName.length() - ending.length(), ending.length(), ending) == 0
                    ){
                   LOG(DEBUG) << Form("Sample %s - Keep this bin", sample->fName.c_str()) << "\n";
                }
                else {
                    LOG(DEBUG) << Form("Sample %s - Delete this bin", sample->fName.c_str()) << "\n";
                    histogram->SetBinContent(ib+1,0.);
                }
            }
        }
    }
    return histogram;
}

void NtupleReader::ReadNtuples(){
    LOG(INFO) << "-------------------------------------------\n";
    LOG(INFO) << "Reading ntuples...\n";
    LOG(INFO) << "-------------------------------------------\n";
    //
    // Import custom functions from .C files
    //
    for(const auto& path : fFitter->fCustomIncludePaths){
        LOG(INFO) << "  Adding include path " << path << " ...\n";
        gROOT->ProcessLineSync((".include "+path).c_str());
    }
    for(const auto& file : fFitter->fCustomFunctions){
        LOG(INFO) << "  Loading function from " << file << " ...\n";
        gROOT->ProcessLineSync((".L "+file+"+").c_str());
    }
    for(const auto& line : fFitter->fCustomFunctionsExecutes){
        LOG(INFO) << "  Executing function line: " << line << " ...\n";
        gROOT->ProcessLine(line.c_str());
    }
    //
    // Loop on regions
    //
    int ChannelIdx = 0;
    for(auto& region: fFitter->fRegions) {

        LOG(INFO) << "  Reading region " << region->fName << " ...\n";

        if (TRExFitter::SPLITHISTOFILES) {
            // Check if the input directory has been found before writing to it, avoid a segfault
            if(!fFitter->fFiles.at(ChannelIdx)) {
                LOG(ERROR) << "Input Histograms folder is not found. If you specify InputFolder in Job block, check the directory exist\n";
                exit(EXIT_FAILURE);
            }
            fFitter->fFiles.at(ChannelIdx)->cd();
        } else {
            if (fFitter->fFiles.empty() || !fFitter->fFiles.at(0)) {
                LOG(ERROR) << "Input Histograms folder is not found. If you specify InputFolder in Job block, check the directory exist\n";
                exit(EXIT_FAILURE);
            }
            fFitter->fFiles.at(0)->cd();
        }

        //
        if(region->fBinTransfo != "") fFitter->ComputeBinning(ChannelIdx);
        if(region->fCorrVar1 != ""){
            if(region->fCorrVar2 == ""){
                LOG(WARNING) << "Only first correlation variable defined, do not read region : " << region->fName << "\n";
                continue;
            }
            LOG(DEBUG) << "Calling the function 'DefineVariable(ChannelIdx)'\n";
            DefineVariable(ChannelIdx);
        }
        else if(region->fCorrVar2 != ""){
            LOG(WARNING) << "Only second correlation variable defined, do not read region : " << region->fName << "\n";
            continue;
        }
        //
        // Loop on samples
        //
        for(auto& sample: fFitter->fSamples) {

            //
            // skip sample if it shouldn't be in current region
            //
            if( Common::FindInStringVector(sample->fRegions,region->fName)<0 ) continue;
            //
            LOG(INFO) << "    Reading sample " << sample->fName << "\n";
            //
            /*  Read nominal histogram */
            TH1D* hNom = GetHistogram(region.get(), sample.get(), nullptr, true, false); // last arg irrelevant if syst = nullptr
            TH1* h_orig = static_cast<TH1*>(hNom->Clone( Form("%s_orig",hNom->GetName()) )); // _orig before smoothing/symmetrising

            // Save histogram into TREx object (SampleHist)
            std::shared_ptr<SampleHist> sh = region->SetSampleHist( sample.get(), hNom);
            // From now on sh->fHist can be used to get nominal histogram
            sh->fHist_orig.reset(h_orig);
            sh->fHist_orig->SetName( Form("%s_orig", sh->fHist->GetName()) ); // fix the name


            /* Add NormFactors */
            for(const auto& inorm : sample->fNormFactors) {
                //
                // skip NF if it shouldn't be in current region
                //
                if(inorm->fRegions.size()>0 && Common::FindInStringVector(inorm->fRegions, region->fName)<0  ) continue;
                if(inorm->fExclude.size()>0 && Common::FindInStringVector(inorm->fExclude, region->fName)>=0 ) continue;
                //
                LOG(DEBUG) << "Adding NormFactor " << inorm->fName << "\n";
                //
                sh->AddNormFactor(inorm);
            }

            /* Add ShapeFactors */
            for(const auto& ishape : sample->fShapeFactors) {
                //
                // skip ShapeFactor if it shouldn't be in current region
                //
                if(ishape->fRegions.size()>0 && Common::FindInStringVector(ishape->fRegions, region->fName)<0  ) continue;
                if(ishape->fExclude.size()>0 && Common::FindInStringVector(ishape->fExclude, region->fName)>=0 ) continue;
                //
                LOG(DEBUG) << "Adding ShapeFactor " << ishape->fName << "\n";
                //
                sh->AddShapeFactor(ishape);
            }

            /*
                Build systematic histograms
            */


            for(const auto& syst_shared_ptr :sample->fSystematics) {

                Systematic* syst = syst_shared_ptr.get();

                //
                // skip Systemtic if it shouldn't be in current region/sample
                //
                if( syst->fRegions.size()>0 && Common::FindInStringVector(syst->fRegions, region->fName)<0  ) continue;
                if( syst->fExclude.size()>0 && Common::FindInStringVector(syst->fExclude, region->fName)>=0 ) continue;
                if( syst->fExcludeRegionSample.size()>0 &&
                    Common::FindInStringVectorOfVectors(syst->fExcludeRegionSample, region->fName, sample->fName)>=0 ) continue;

                //
                LOG(DEBUG) << "    Reading systematic " << syst->fName << "\n";
                //


                /* Handle OVERALL systematics */
                if(syst->fType==Systematic::OVERALL){

                    std::shared_ptr<SystematicHist> SystHisto = region->GetSampleHist(sample->fName)->AddOverallSyst(syst->fName,syst->fStoredName,syst->fOverallUp,syst->fOverallDown);
                    SystHisto->fSystematic = syst_shared_ptr;

                    /* Apply scaling of systematic (make them larger/smaller) */
                    SystHisto->fScaleUp    = syst->fScaleUp;
                    if(syst->fScaleUpRegions.size()!=0)
                        if(syst->fScaleUpRegions[region->fName]!=0)
                            SystHisto->fScaleUp *= syst->fScaleUpRegions[region->fName];

                    SystHisto->fScaleDown = syst->fScaleDown;
                    if(syst->fScaleDownRegions.size()!=0)
                        if(syst->fScaleDownRegions[region->fName]!=0)
                            SystHisto->fScaleDown *= syst->fScaleDownRegions[region->fName];
                    // Move on to next syst, no more processing
                    continue;
                }

                /* Handle MCStats (gammas) */
                if(syst->fType == Systematic::STAT){
                    std::shared_ptr<SystematicHist> SystHisto = region->GetSampleHist(sample->fName)->AddStatSyst(syst->fName,syst->fStoredName,syst->fBins[0]);
                    SystHisto->fSystematic = syst_shared_ptr;
                    // Move on to next syst, no more processing
                    continue;
                }

                /*
                    Now we process other types of systematics that need to read new samples/trees/branches
                    i.e. Histogram systematics, they have an up and a down variation defined
                */

                // Do a check on existence of reference sample in config for a syst if needed
                if(syst->fReferenceSample!="") {
                    // If reference sample not parsed from config, break
                    if( fFitter->GetSample(syst->fReferenceSample) == nullptr ){
                        LOG(ERROR) << "Reference sample: " << syst->fReferenceSample << " is not defined for region: " << region->fName << ".\n";
                        LOG(ERROR) << "This probably means that you run over a specific sample, you need to run over the reference sample as well.\n";
                        LOG(ERROR) << "Ignoring SeparateSample setting.\n";
                        exit(EXIT_FAILURE);
                    }
                    // If reference sample parsed but its histogram has not been built, break (ordering issue)
                    else if ( region->GetSampleHist(syst->fReferenceSample) ==nullptr ){
                        LOG(ERROR) << "Histogram of reference sample: " << syst->fReferenceSample << " is not built for region: " << region->fName << ".\n";
                        LOG(ERROR) << "Make sure you define non-GHOST reference samples BEFORE the sample using them for systematics.\n";
                        exit(EXIT_FAILURE);
                    }
                    else {
                        //sample = fFitter->GetSample(syst->fReferenceSample);
                    }
                }


                TH1D* hUp = nullptr;
                TH1D* hDown = nullptr;

                /* If systematic is a DUMMY */
                if(Common::FindInStringVector(syst->fDummyForSamples, sample->fName)>=0){
                    //
                    LOG(INFO) << "Systematic " << syst->fName << " set as dummy for sample " << sample->fName << " in region " << region->fName << "\n";
                    //

                    hUp   = static_cast<TH1D*>(sh->fHist->Clone());
                    hDown = static_cast<TH1D*>(sh->fHist->Clone());
                    std::shared_ptr<SystematicHist> SystHisto = sh->AddHistoSyst(syst->fName, syst->fStoredName, hUp, hDown);
                    SystHisto->fSystematic = syst_shared_ptr;

                    /* Apply scaling of systematic (make them larger/smaller) */
                    SystHisto->fScaleUp = syst->fScaleUp;
                    if(syst->fScaleUpRegions.size()!=0)
                        if(syst->fScaleUpRegions[region->fName]!=0)
                            SystHisto->fScaleUp *= syst->fScaleUpRegions[region->fName];

                    SystHisto->fScaleDown = syst->fScaleDown;
                    if(syst->fScaleDownRegions.size()!=0)
                        if(syst->fScaleDownRegions[region->fName]!=0)
                            SystHisto->fScaleDown *= syst->fScaleDownRegions[region->fName];
                    // Move on to next syst, no more processing
                    continue;
                }

                /* If systematic is a not DUMMY, process up and down variations */

                std::vector<std::string> variations = {"Up", "Down"};
                std::vector<TH1D*> VarHistos; // vector to store up and down histograms

                for (const auto& variation: variations){
                    bool VarIsUp   = (variation == "Up");
                    bool hasVar;
                    if (VarIsUp)    hasVar = syst->fHasUpVariation;
                    else            hasVar = syst->fHasDownVariation;

                    TH1D* hVar(nullptr);

                    if (!hasVar) {
                        hVar   = static_cast<TH1D*>(region->GetSampleHist(sample->fName )->fHist.get());
                        LOG(DEBUG) << "  Syst " << syst->fName << " has no " << variation << " variation defined. Will use nominal for now.\n";
                        VarHistos.push_back(hVar);
                        continue;
                    }

                    hVar = GetHistogram(region.get(), sample.get(), syst, VarIsUp, false);

                    /*
                        Handle variation subtraction e.g. for JER
                    */
                    std::vector<std::string> NtuplePathsVarSubtractSample;
                    std::vector<std::string> NtupleFilesVarSubtractSample;
                    std::vector<std::string> NtupleNamesVarSubtractSample;
                    std::string NtupleFileSufVarSubtractSample;
                    std::string NtupleNameSufVarSubtractSample;


                    if (variation == "Up"){
                        NtuplePathsVarSubtractSample = syst->fNtuplePathsUpSubtractSample;
                        NtupleFilesVarSubtractSample = syst->fNtupleFilesUpSubtractSample;
                        NtupleNamesVarSubtractSample = syst->fNtupleNamesUpSubtractSample;
                        NtupleFileSufVarSubtractSample = syst->fNtupleFileSufUpSubtractSample;
                        NtupleNameSufVarSubtractSample = syst->fNtupleNameSufUpSubtractSample;
                    }
                    else{
                        NtuplePathsVarSubtractSample = syst->fNtuplePathsDownSubtractSample;
                        NtupleFilesVarSubtractSample = syst->fNtupleFilesDownSubtractSample;
                        NtupleNamesVarSubtractSample = syst->fNtupleNamesDownSubtractSample;
                        NtupleFileSufVarSubtractSample = syst->fNtupleFileSufDownSubtractSample;
                        NtupleNameSufVarSubtractSample = syst->fNtupleNameSufDownSubtractSample;
                    }

                    const bool VarSubtract   = ( !NtuplePathsVarSubtractSample.empty() ||
                                                 !NtupleFilesVarSubtractSample.empty() ||
                                                 !NtupleNamesVarSubtractSample.empty() ||
                                                 NtupleFileSufVarSubtractSample!=""    ||
                                                 NtupleNameSufVarSubtractSample!="");


                    if (VarSubtract) {
                        TH1D* hSubtract = GetHistogram(region.get(), sample.get(), syst, VarIsUp, VarSubtract);
                        // Systematic = variation - subtracted + nominal (later replaced with Reference if available)
                        hVar->Add(hSubtract, -1);
                        hVar->Add(sh->fHist.get());
                    }

                    // Obtain relative variation with respect to a reference sample if available
                    if(syst->fReferenceSample!="" && region->GetSampleHist(syst->fReferenceSample)!=nullptr){
                        // EFT block
                        if(sh->fName.rfind("SM_", 0) == 0 && sh->fName.find("_bin") != std::string::npos && sample->fType == Sample::SampleType::EFT){
                            size_t found = sh->fName.find("_bin");
                            int BinNum = Common::convertStoNum<int>(sh->fName.substr(found+4));

                            LOG(DEBUG) << "Sample being varied is an EFT sample at a specific bin\n";
                            LOG(DEBUG) << "Bin number:\t" << BinNum << "\n";
                            LOG(DEBUG) << variation << " Variation BEFORE reference sample:\n";
                            LOG(DEBUG) << Form("%f",hVar->GetBinContent(BinNum+1)) << "\n";

                            const TH1* href = region->GetSampleHist(syst->fReferenceSample)->fHist.get();
                            const TH1* hnom = region->GetSampleHist(sample->fName )->fHist.get();

                            // get copies with no error on bins
                            auto hrefTmp = Common::GetHistCopyNoError(href);
                            auto hnomTmp = Common::GetHistCopyNoError(hnom);

                            double divider    = hrefTmp.get()->GetBinContent(BinNum+1);
                            double multiplier = hnomTmp.get()->GetBinContent(BinNum+1);


                            double newcontent = hVar->GetBinContent(BinNum+1)*multiplier/divider;
                            hVar->SetBinContent(BinNum+1, newcontent);
                            LOG(DEBUG) << variation << " Variation AFTER reference sample:\n";
                            LOG(DEBUG) << Form("%f",hVar->GetBinContent(BinNum+1)) << "\n";
                        }
                        else{
                            TH1* href = region->GetSampleHist(syst->fReferenceSample)->fHist.get();
                            TH1* hnom = region->GetSampleHist(sample->fName )->fHist.get();
                            // Protection added: fix empty bins before starting to divide and multiply
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++)  if(href->GetBinContent(i_bin)<=1e-6) href->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<hVar->GetNbinsX()+2;i_bin++) if(hVar ->GetBinContent(i_bin)<=1e-6) hVar ->SetBinContent(i_bin,1e-6);
                            for(int i_bin=0;i_bin<href->GetNbinsX()+2;i_bin++)  if(href->GetBinContent(i_bin)<=1e-6) hVar ->SetBinContent(i_bin,1e-6); // avoid multiplying by 1e6
                            //
                            if (VarSubtract){
                                hVar->Add( hnom, -1 );
                                hVar->Add( href );
                            }
                            hVar->Divide(  href );
                            hVar->Multiply( hnom );
                        }
                    } // End reference sample block
                    //if (hVar==nullptr)  hVar   = static_cast<TH1D*>(region->GetSampleHist(sample->fName )->fHist.get());
                    VarHistos.push_back(hVar);

                } // end loop over up down

                hUp = static_cast<TH1D*>(VarHistos.at(0)->Clone(Form("h_%s_%s_%sUp",region->fName.c_str(),sample->fName.c_str(),syst->fStoredName.c_str())));
                hDown = static_cast<TH1D*>(VarHistos.at(1)->Clone(Form("h_%s_%s_%sDown",region->fName.c_str(),sample->fName.c_str(),syst->fStoredName.c_str())));

                /*
                    Add systematic (both variations) now to sample trex object
                */
                std::shared_ptr<SystematicHist> SystHisto = sh->AddHistoSyst(syst->fName,
                                                                             syst->fStoredName,
                                                                             hUp,
                                                                             hDown);
                SystHisto->fSystematic = syst_shared_ptr;

                SystHisto->fScaleUp = syst->fScaleUp;
                if(syst->fScaleUpRegions.size()!=0)
                    if(syst->fScaleUpRegions[region->fName]!=0)
                        SystHisto->fScaleUp *= syst->fScaleUpRegions[region->fName];

                SystHisto->fScaleDown = syst->fScaleDown;
                if(syst->fScaleDownRegions.size()!=0)
                    if(syst->fScaleDownRegions[region->fName]!=0)
                        SystHisto->fScaleDown *= syst->fScaleDownRegions[region->fName];
            } // end syst loop
        }// end sample loop
        ChannelIdx++;
    }   // end region loop
}

void NtupleReader::DefineVariable(int regIter){
    TH1::StatOverflows(true);  //////  What is the defaut in root for this ???
    LOG(DEBUG) << "//////// --------\n";
    LOG(DEBUG) << "// DEBUG CORR VAR\n";
    TH1* h1 = new TH1D("h1","h1",1,-2000.,1000.);
    TH1* h2 = new TH1D("h2","h2",1,-2000.,1000.);
    std::string fullSelection;
    std::string fullMCweight;
    std::vector<std::string> fullPaths;

    // copy of NtupleReading function.
    for(std::size_t i_smp = 0; i_smp < fFitter->fSamples.size(); ++i_smp) {
        LOG(DEBUG) << "Processing sample : " << fFitter->fSamples[i_smp]->fName << "\n";
        if(fFitter->fSamples[i_smp]->fType==Sample::SampleType::DATA) continue;
        if(Common::FindInStringVector(fFitter->fSamples[i_smp]->fRegions,fFitter->fRegions[regIter]->fName)<0 ) continue;
        LOG(DEBUG) << " -> is used in the considered region\n";
        //
        // set selection, weight and paths (no variables)
        fullSelection = fFitter->FullSelection(  fFitter->fRegions[regIter].get(),fFitter->fSamples[i_smp].get());
        fullMCweight  = fFitter->FullWeight(     fFitter->fRegions[regIter].get(),fFitter->fSamples[i_smp].get());
        fullPaths     = fFitter->FullNtuplePaths(fFitter->fRegions[regIter].get(),fFitter->fSamples[i_smp].get());
        //
        for(unsigned int i_path=0;i_path<fullPaths.size();i_path++){
            LOG(DEBUG) << " -> Retrieving : " << fFitter->fRegions[regIter]->fCorrVar1 << " w/ weight " << fullMCweight << "*" << fullSelection <<
                                                        " from " <<  fullPaths[i_path] << "\n";
            TH1* htmp1 = new TH1D("htmp1","htmp1",1,-2000.,1000.);
            TH1* htmp2 = new TH1D("htmp2","htmp2",1,-2000.,1000.);
            TChain *t = new TChain();
            t->Add(fullPaths[i_path].c_str());
            t->Draw(Form("%s>>htmp1",fFitter->fRegions[regIter]->fCorrVar1.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
            t->Draw( Form("%s>>htmp2",fFitter->fRegions[regIter]->fCorrVar2.c_str()), Form("(%s)*(%s)",fullMCweight.c_str(),fullSelection.c_str()), "goff");
            delete t;
            //
            if(fFitter->fSamples[i_smp]->fType!=Sample::SampleType::DATA && fFitter->fSamples[i_smp]->fNormalizedByTheory) htmp1 -> Scale(fFitter->fLumi);
            if(fFitter->fSamples[i_smp]->fLumiScales.size()>i_path)  htmp1 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[i_path]);
            else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1) htmp1 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
            //
            if(fFitter->fSamples[i_smp]->fType!=Sample::SampleType::DATA && fFitter->fSamples[i_smp]->fNormalizedByTheory) htmp2 -> Scale(fFitter->fLumi);
            if(fFitter->fSamples[i_smp]->fLumiScales.size()>i_path)  htmp2 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[i_path]);
            else if(fFitter->fSamples[i_smp]->fLumiScales.size()==1) htmp2 -> Scale(fFitter->fSamples[i_smp]->fLumiScales[0]);
            //
            h1->Add(htmp1);
            h2->Add(htmp2);
            delete htmp1;
            delete htmp2;
        }
    }
    double mean1 = h1->GetMean();
    double rms1 = h1->GetRMS();
    double mean2 = h2->GetMean();
    double rms2 = h2->GetRMS();
    LOG(DEBUG) << "The new variable : ( ( (" << fFitter->fRegions[regIter]->fCorrVar1 << ") - " << mean1 << " )*( (" << fFitter->fRegions[regIter]->fCorrVar2 << ")-" << mean2 << " ) )/( " << rms1 << " * " << rms2 << ")\n";
    fFitter->fRegions[regIter]->fVariable = Form("( ( (%s)-%f )*( (%s)-%f ) )/( %f * %f )",fFitter->fRegions[regIter]->fCorrVar1.c_str(),mean1,fFitter->fRegions[regIter]->fCorrVar2.c_str(),mean2,rms1,rms2);
    TH1::StatOverflows(false);  //////  What is the defaut in root for this ???
    delete h1;
    delete h2;
}
