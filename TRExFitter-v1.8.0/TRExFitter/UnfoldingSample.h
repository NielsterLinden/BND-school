#ifndef UNFOLDINGSAMPLE_H_
#define UNFOLDINGSAMPLE_H_

#include <memory>
#include <string>
#include <vector>
#include <map>

class Region;
class Sample;
class Unfolding;

class UnfoldingSample {

public:

    enum class TYPE {
        STANDARD,
        GHOST
    };

    enum class GAMMAS {
        SEPARATED,
        DISABLED
    };

    explicit UnfoldingSample();
    ~UnfoldingSample() = default;

    UnfoldingSample(const UnfoldingSample& s) = default;
    UnfoldingSample(UnfoldingSample&& s) = default;
    UnfoldingSample& operator=(const UnfoldingSample& s) = default;
    UnfoldingSample& operator=(UnfoldingSample&& s) = default;

    inline void SetName(const std::string& s) {fName = s;}
    inline const std::string& GetName() const {return fName;}
    inline void SetTitle(const std::string& s) {fTitle = s;}
    inline const std::string& GetTitle() const {return fTitle;}
    inline void SetFillColor(const int c) {fFillColor = c;}
    inline int GetFillColor() const {return fFillColor;}
    inline void SetLineColor(const int c) {fLineColor = c;}
    inline int GetLineColor() const {return fLineColor;}
    inline void SetHasResponse(const bool r) {fHasResponse = r;}
    inline bool GetHasResponse() const {return fHasResponse;}
    inline void SetHasAcceptance(const bool a) {fHasAcceptance = a;}
    inline bool GetHasAcceptance() const {return fHasAcceptance;}
    inline void SetType(const TYPE type) {fType = type;}
    inline TYPE GetType() const {return fType;}
    inline void SetGammas(const GAMMAS gammas){fGammas = gammas;}
    inline GAMMAS GetGammas() const {return fGammas;}
    inline void SetNoPruning(const bool noPruning){fNoPruning = noPruning;}
    inline bool GetNoPruning() const {return fNoPruning;}

    std::vector<std::shared_ptr<Sample> > ConvertToSample(const Region* reg,
                                                          const Unfolding* unfolding,
                                                          const std::string& name) const;
    
    std::vector<std::string> fResponseMatrixFiles;
    std::vector<std::string> fResponseMatrixNames;
    std::vector<std::string> fResponseMatrixPaths;
    std::vector<std::string> fResponseMatrixFileSuffs;
    std::vector<std::string> fResponseMatrixNameSuffs;
    std::vector<std::string> fResponseMatrixPathSuffs;
    
    std::vector<std::string> fAcceptanceFiles;
    std::vector<std::string> fAcceptanceNames;
    std::vector<std::string> fAcceptancePaths;
    std::vector<std::string> fAcceptanceFileSuffs;
    std::vector<std::string> fAcceptanceNameSuffs;
    std::vector<std::string> fAcceptancePathSuffs;
    
    std::vector<std::string> fSelectionEffFiles;
    std::vector<std::string> fSelectionEffNames;
    std::vector<std::string> fSelectionEffPaths;
    std::vector<std::string> fSelectionEffFileSuffs;
    std::vector<std::string> fSelectionEffNameSuffs;
    std::vector<std::string> fSelectionEffPathSuffs;
    
    std::vector<std::string> fMigrationFiles;
    std::vector<std::string> fMigrationNames;
    std::vector<std::string> fMigrationPaths;
    std::vector<std::string> fMigrationFileSuffs;
    std::vector<std::string> fMigrationNameSuffs;
    std::vector<std::string> fMigrationPathSuffs;
    
    std::vector<std::string> fRegions;
    std::map<int,std::vector<std::string>> fSubSampleRegions;
    std::string fUnfoldingName;

private:
    std::string fName;
    std::string fTitle;
    int fFillColor;
    int fLineColor;
    bool fHasResponse;
    bool fHasAcceptance;
    TYPE fType;
    GAMMAS fGammas;
    bool fNoPruning;
};

#endif
