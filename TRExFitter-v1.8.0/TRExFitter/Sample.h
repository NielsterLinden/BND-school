#ifndef SAMPLE_H
#define SAMPLE_H

/// c++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

/// Forward class declaration
class NormFactor;
class ShapeFactor;
class Systematic;

class Sample {
public:

    enum class SampleType{
        BACKGROUND, // 0
        SIGNAL, // 1
        DATA, // 2
        GHOST, // 3
        EFT // 4
    };

    explicit Sample(const std::string& name, SampleType type);
    ~Sample();
    //Sample(const Sample& s) = delete;
    Sample(const Sample& s) = default;  // copy constructor
    Sample(Sample&& s) = delete;
    //Sample& operator=(const Sample& s) = delete;
    Sample& operator=(const Sample& s) = default;
    Sample& operator=(Sample&& s) = delete;




    // -------
    // Methods
    // -------

    // comsetics
    void SetTitle(const std::string& title);
    void SetFillColor(int color);
    void SetLineColor(int color);
    void NormalizedByTheory(const bool norm);

    // read from ntupes
    void AddNtuplePath(const std::string& path);
    void AddNtupleFile(const std::string& file);
    void AddNtupleName(const std::string& name);
    void SetMCweight(const std::string& weight);
    void SetSelection(const std::string& selection);

    // read from histos
    void AddHistoPath(const std::string& path);
    void AddHistoFile(const std::string& file);
    void AddHistoName(const std::string& name);

    // norm factors and systs
    void AddNormFactor(std::shared_ptr<NormFactor> factor);
    void AddShapeFactor(std::shared_ptr<ShapeFactor> factor);
    void AddSystematic(std::shared_ptr<Systematic> syst);
    std::shared_ptr<NormFactor> AddNormFactor(const std::string& name,double nominal=1.,double min=0.,double max=10.,bool isConst=false);
    std::shared_ptr<ShapeFactor> AddShapeFactor(const std::string& name,double nominal=1.,double min=0.,double max=10.,bool isConst=false);
    std::shared_ptr<Systematic> AddSystematic(const std::string& name,int type=0,double up=0.,double down=0.);
    bool HasNormFactor(const std::string& name) const;
    bool HasSystematic(const std::string& name) const;
    bool HasNuisanceParameter(const std::string& name) const;

    // -------
    // Members
    // -------

    std::string fName;
    bool fUseFriend;
    SampleType fType;
    std::string fFitName;
    std::string fTitle;
    std::string fTexTitle;
    std::string fGroup;
    int fFillColor;
    int fLineColor;
    bool fNormalizedByTheory;
    std::vector<std::string> fRegions;
    std::vector<double> fLumiScales;
    std::string fIgnoreSelection;
    std::string fIgnoreWeight;
    bool fUseMCStat;
    bool fUseSystematics;
    std::string fDivideBy;
    std::string fMultiplyBy;
    std::vector<std::string> fSubtractSamples;
    std::vector<std::string> fAddSamples;
    std::string fNormToSample;
    bool fSmooth;
    std::map<std::string,bool> fIsMorph;
    std::map<std::string,double> fMorphValue;

    // to read from ntuples
    std::string fSelection;
    std::string fMCweight;
    std::vector<std::string> fNtuplePaths;
    std::vector<std::string> fNtuplePathSuffs;
    std::vector<std::string> fFriendPaths;
    std::vector<std::string> fFriendPathSuffs;
    std::vector<std::string> fNtupleFiles;
    std::vector<std::string> fNtupleFileSuffs;
    std::vector<std::string> fNtupleNames;
    std::vector<std::string> fNtupleNameSuffs;
    std::string fGridDID;

    // to read from histograms
    // <path>/<file>.root/<name>
    std::vector<std::string> fHistoPaths;
    std::vector<std::string> fHistoPathSuffs;
    std::vector<std::string> fHistoFiles;
    std::vector<std::string> fHistoFileSuffs;
    std::vector<std::string> fHistoNames;
    std::vector<std::string> fHistoNameSuffs;
    std::string fHistoFolderName;

    // systematics & norm.factors
    std::vector < std::shared_ptr<Systematic> > fSystematics;
    std::vector < std::shared_ptr<NormFactor> > fNormFactors;
    std::vector < std::shared_ptr<ShapeFactor> > fShapeFactors;

    std::vector<std::pair<std::string,std::string>> fAsimovReplacementFor = {};

    bool fSeparateGammas;
    double fMCstatScale;
    std::vector<std::vector<std::string>> fCorrelateGammasInRegions;
    std::string fCorrelateGammasWithSample;

    std::string fSystFromSample;
    bool fIsFolded;
    bool fUseGaussianShapeSysConstraint;

    std::string fEFTSMReference;
    std::string fEFTParam;
    std::string fEFTTitle;
    std::string fEFTValue;
    std::string fUnfoldingName;
    std::string fUnfoldingSampleName;
    bool fNoPruning;
    std::vector<std::pair<std::string, double> > fTemplateMorphing;
    std::string fAlternativeRefInTemplateMorphing;
    std::string fTemplateTarget;
};

#endif

