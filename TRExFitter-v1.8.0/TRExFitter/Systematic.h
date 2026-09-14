#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

/// Framework includes
#include "TRExFitter/HistoTools.h"

/// c++ includes
#include <map>
#include <string>
#include <vector>

class Systematic {
public:

    enum COMBINATIONTYPE {
        STANDARDDEVIATION = 0,
        ENVELOPE = 1,
        STANDARDDEVIATIONNODDOF = 2,
        HESSIAN = 3,
        SUMINSQUARES = 4
    };

    enum SystType{
         OVERALL, // 0
         SHAPE, // 1
         HISTO, // 2
         STAT // 3
    };


    explicit Systematic(const std::string& name,int type=0,double up=0.,double down=0.);
    
    Systematic(const Systematic& sys) = default;  // copy constructor
    Systematic& operator=(const Systematic& sys) = default;
    Systematic(Systematic&& sys) = delete;
    Systematic& operator=(Systematic&& sys) = delete;

    ~Systematic() = default;

    // -------
    // Members
    // -------

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;
    std::string fSubCategory;
    std::string fStoredName;
    int fType;
    int fSmoothType;
    bool fPreSmoothing;
    bool fSampleSmoothing;
    HistoTools::SmoothOption fSampleSmoothOption;
    HistoTools::SymmetrizationType fSymmetrisationType;
    std::string fReferenceSample;
    std::string fReferenceSmoothing;
    std::string fReferencePruning;

    double fOverallUp;
    double fOverallDown;

    double fScaleUp;
    double fScaleDown;

    std::map<std::string,double> fScaleUpRegions;
    std::map<std::string,double> fScaleDownRegions;

    bool fHasUpVariation;
    bool fHasDownVariation;

    bool fIsFreeParameter;
    bool fIsCorrelated;
    bool fIsShapeAccDecorr;
    bool fIsRegionDecorr;
    bool fIsSampleDecorr;
    bool fIsShapeOnly;
    bool fIsNormOnly;
    bool fNoPruning;

    std::vector<std::string> fRegions;
    std::vector<std::string> fExclude;
    std::vector<std::vector<std::string> > fExcludeRegionSample;
    std::vector<std::string> fDropShapeIn;
    std::vector<std::string> fDropNormIn;
    std::vector<std::string> fDropNormSpecialIn;
    std::vector<std::string> fKeepNormForSamples;
    std::vector<std::string> fDummyForSamples;
    std::vector<int> fBins;

    // from ntuples - up
    std::string fWeightUp;
    std::string fWeightSufUp;
    std::vector<std::string> fNtuplePathsUp;
    std::string fNtuplePathSufUp;
    std::vector<std::string> fFriendPathsUp;
    std::string fFriendPathSufUp;
    std::vector<std::string> fNtupleFilesUp;
    std::string fNtupleFileSufUp;
    std::vector<std::string> fNtupleNamesUp;
    std::string fNtupleNameSufUp;

    std::vector<std::string> fNtuplePathsUpSubtractSample;
    std::vector<std::string> fFriendPathsUpSubtractSample;
    std::vector<std::string> fNtupleFilesUpSubtractSample;
    std::vector<std::string> fNtupleNamesUpSubtractSample;
    std::string fNtupleFileSufUpSubtractSample;
    std::string fNtupleNameSufUpSubtractSample;
    
    std::vector<std::string> fHistoPathsUpSubtractSample;
    std::vector<std::string> fHistoFilesUpSubtractSample;
    std::vector<std::string> fHistoNamesUpSubtractSample;
    std::string fHistoNameSufUpSubtractSample;

    // from ntuples - down
    std::string fWeightDown;
    std::string fWeightSufDown;
    std::vector<std::string> fNtuplePathsDown;
    std::string fNtuplePathSufDown;
    std::vector<std::string> fNtupleFilesDown;
    std::string fNtupleFileSufDown;
    std::vector<std::string> fNtupleNamesDown;
    std::string fNtupleNameSufDown;
    std::vector<std::string> fFriendPathsDown;
    std::string fFriendPathSufDown;

    std::vector<std::string> fNtuplePathsDownSubtractSample;
    std::vector<std::string> fFriendPathsDownSubtractSample;
    std::vector<std::string> fNtupleFilesDownSubtractSample;
    std::vector<std::string> fNtupleNamesDownSubtractSample;
    std::string fNtupleFileSufDownSubtractSample;
    std::string fNtupleNameSufDownSubtractSample;
    std::string fHistoFileSufUpSubtractSample;

    std::vector<std::string> fHistoPathsDownSubtractSample;
    std::vector<std::string> fHistoFilesDownSubtractSample;
    std::vector<std::string> fHistoNamesDownSubtractSample;
    std::string fHistoNameSufDownSubtractSample;
    std::string fHistoFileSufDownSubtractSample;

    std::string fIgnoreWeight;

    // from histos - up
    std::vector<std::string> fHistoPathsUp;
    std::string fHistoPathSufUp;
    std::vector<std::string> fHistoFilesUp;
    std::string fHistoFileSufUp;
    std::vector<std::string> fHistoNamesUp;
    std::string fHistoNameSufUp;

    // from histos - down
    std::vector<std::string> fHistoPathsDown;
    std::string fHistoPathSufDown;
    std::vector<std::string> fHistoFilesDown;
    std::string fHistoFileSufDown;
    std::vector<std::string> fHistoNamesDown;
    std::string fHistoNameSufDown;
    std::vector<std::string> fHistoFolderNamesUp;
    std::vector<std::string> fHistoFolderNamesDown;
    std::vector<std::string> fHistoFolderSubtractNamesUp;
    std::vector<std::string> fHistoFolderSubtractNamesDown;

    //
    std::string fSampleUp;
    std::string fSampleDown;

    std::vector<std::string> fSamples;

    std::string fCombineName;
    COMBINATIONTYPE fCombineType;
    HistoTools::FORCESHAPETYPE fForceShape;

};

#endif
