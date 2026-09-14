// Header include
#include "TRExFitter/Common.h"

// Framework includes
#include "TRExFitter/CorrelationMatrix.h"
#include "TRExFitter/FitResults.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/Logger.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/ShapeFactor.h"
#include "TRExFitter/SystematicHist.h"
#include "TRExFitter/Unfolding.h"

// ROOT stuff
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TPad.h"
#include "RooFitResult.h"
#include "TString.h"
#include "TSystem.h"
#include "xRooFit/xRooFit.h"
#include "xRooFit/xRooHypoSpace.h"

#include "RooStats/RooStatsUtils.h"

// c++ stuff
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <regex>

namespace fs = std::filesystem;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// VARIABLES
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int TRExFitter::DEBUGLEVEL = 1;
bool TRExFitter::SHOWYIELDS = false;
bool TRExFitter::SHOWSTACKSIG = true;
bool TRExFitter::ADDSTACKSIG = true;
bool TRExFitter::SHOWNORMSIG = false;
bool TRExFitter::SHOWOVERLAYSIG = false;
double TRExFitter::SHOWOVERLAYSIG_CUSTOMSCALE = -1.0;
bool TRExFitter::SHOWCHI2 = false;
bool TRExFitter::SHOWSTACKSIG_SUMMARY = true;
bool TRExFitter::SHOWNORMSIG_SUMMARY = false;
bool TRExFitter::SHOWOVERLAYSIG_SUMMARY = false;
bool TRExFitter::LEGENDLEFT = false;
bool TRExFitter::LEGENDRIGHT = false;
bool TRExFitter::PREFITONPOSTFIT = false;
bool TRExFitter::SYSTCONTROLPLOTS = true;
bool TRExFitter::SYSTERRORBARS = true;
bool TRExFitter::SYSTDATAPLOT = false;
bool TRExFitter::SPLITHISTOFILES = false;
bool TRExFitter::HISTOCHECKCRASH = true;
bool TRExFitter::GUESSMCSTATERROR = true;
bool TRExFitter::CORRECTNORMFORNEGATIVEINTEGRAL = false;
bool TRExFitter::REMOVEXERRORS = true;
double TRExFitter::CORRELATIONTHRESHOLD = -1.;
bool TRExFitter::MERGEUNDEROVERFLOW = false;
bool TRExFitter::OPRATIO = false;
bool TRExFitter::NORATIO = false;
bool TRExFitter::USEFOLDERSTRUCTURE = true;
bool TRExFitter::USEHEPDATAROUNDING = false;
bool TRExFitter::KEEPDATAERRORS = false;
std::map <std::string,std::string> TRExFitter::SYSTMAP;
std::map <std::string,std::string> TRExFitter::SYSTTEX;
std::map <std::string,std::string> TRExFitter::NPMAP;
std::vector <std::string> TRExFitter::IMAGEFORMAT;
std::string TRExFitter::EXPERIMENT_LABEL = "ATLAS";
//
std::map<std::string,double> TRExFitter::OPTION;
std::map<std::string, std::unique_ptr<TFile> > TRExFitter::TFILEMAP;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
void TRExFitter::SetDebugLevel(int level){
    DEBUGLEVEL = level;
    if (level == 0) {
        Logger::get().setLogLevel(LoggingLevel::WARNING);
    } else if (level == 1) {
        Logger::get().setLogLevel(LoggingLevel::INFO);
    } else if (level == 2) {
        Logger::get().setLogLevel(LoggingLevel::DEBUG);
    } else if (level >= 3) {
        Logger::get().setLogLevel(LoggingLevel::VERBOSE);
    }
}

void TRExFitter::SetLoggingFormat(std::string param) {
    if (param != "") {
        std::transform(param.begin(), param.end(), param.begin(), ::toupper);
        if (param == "NONE") {
            Logger::get().setIncludeTime(false);
            Logger::get().setIncludeDate(false);
            Logger::get().setIncludeLines(false);
            Logger::get().setIncludeFunction(false);
        } else {
            Logger::get().setIncludeTime(false);
            Logger::get().setIncludeDate(false);
            Logger::get().setIncludeLines(false);
            Logger::get().setIncludeFunction(false);
            const std::vector split = Common::Vectorize(param, ',');
            for (const auto& p : split) {
                if (p == "TIME") {
                    Logger::get().setIncludeTime(true);
                } else if (p == "DATE") {
                    Logger::get().setIncludeDate(true);
                } else if (p == "FILELINES") {
                    Logger::get().setIncludeLines(true);
                } else if (p == "FUNCTION") {
                    Logger::get().setIncludeFunction(true);
                }
            }
        }
    }
}

std::unique_ptr<TH1D>  Common::CombineHistosFromHistosVec(const std::vector<std::shared_ptr<TH1D> >& vec){
    std::unique_ptr<TH1D> result(nullptr);
    for (const auto& ihist : vec) {
        if (!ihist) continue;
        if (result) {
            result->Add(ihist.get());
        } else {
            result.reset(static_cast<TH1D*>(ihist->Clone()));
            result->SetDirectory(nullptr);
        }
    }
    return result;
}

#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__))
#include <glob.h>
void Common::Glob(const std::string& pattern, const std::string& tree, TChain& t) {
    glob_t g;
    glob(pattern.c_str(), GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
    std::vector<std::string> filelist;
    filelist.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        filelist.emplace_back(g.gl_pathv[i]);

        std::string path = g.gl_pathv[i];
        std::string fullpath = path + "/" + tree;
        t.Add(fullpath.c_str());
    }
    globfree(&g);
    if (filelist.size() == 0){
        LOG(WARNING) << "You used wildcards, but added zero files from " << pattern << "\n";
        LOG(WARNING) << "Note that the globbing functionality will not work with the xrootd protocol\n";
    }
}
#else
void Common::Glob(const std::string& pattern, const std::string& tree, TChain& t) {
    LOG(WARNING) << "Detected you are in a Windows environment ... \n";
    LOG(WARNING) << "UNIX-style globbing for NTuple paths is not supported ... \n";
    std::string fullpath = pattern+'/'+tree;
    if (t.Add(fullpath.c_str()) == 0){
        LOG(WARNING) << "You used wildcards, but added zero files from\n";
    }
}
#endif

//__________________________________________________________________________________
//
TH1D* Common::HistFromNtuple(const std::string& ntuple,
                             const std::string& frnd,
                             const std::string& variable,
                             int nbin,
                             double xmin,
                             double xmax,
                             const std::string& selection,
                             const std::string& weight,
                             const std::vector<std::string>& aliases,
                             int Nev) {
    TH1D* h = new TH1D("h","h",nbin,xmin,xmax);
    LOG(VERBOSE) << "    Extracting histogram " << variable << " from  " << ntuple << "  ...\n";
    LOG(VERBOSE) << "        with weight  (" << weight << ")*(" << selection << ")  ...\n";

    bool hasWildcard = false;
    // check whether file actually exists, AccessPathName() returns FALSE if file can be accessed
    // see https://root.cern.ch/root/html602/TSystem.html#TSystem:AccessPathName

    // the string is now "path/to/the/file.root/path/to/the/tree"
    const std::string fileName = ntuple.substr(0,ntuple.find_last_of(".")+5); // find .root and take what is left of it to get the path
    const std::string treeName = ntuple.substr(ntuple.find_last_of(".")+6); // find .root/ and take what is right of it

    if (fileName.find('*') != std::string::npos) hasWildcard = true;
    if (gSystem->AccessPathName(fileName.c_str()) == kTRUE && !hasWildcard ){
        if (TRExFitter::HISTOCHECKCRASH) {
            LOG(ERROR) << "Cannot find input file in: " << fileName << "\n";
            exit(EXIT_FAILURE);
        } else {
            LOG(WARNING) << "Cannot find input file in: " << fileName << "\n";
        }
    }

    TChain t{};
    if (fileName.find('*') == std::string::npos) {
        // does not contain globbing
        t.Add((fileName + "/" + treeName).c_str());
    } else {
        // uses globbing
        Glob(fileName, treeName, t);
    }

    if (frnd != ""){
        LOG(VERBOSE) << "    FriendTree: Extracting histogram " << variable << " from  " << frnd << "  ...\n";
        LOG(VERBOSE) << "                with weight  (" << weight << ")*(" << selection << ")  ...\n";
        const std::string friendName = frnd.substr(0,frnd.find_last_of("/"));
        const std::string friendTree = frnd.substr(frnd.find_last_of("/")+1);
        hasWildcard = false;
        if (friendName.find('*') != std::string::npos) hasWildcard = true;
        if (gSystem->AccessPathName(friendName.c_str()) == kTRUE && !hasWildcard ){
            if (TRExFitter::HISTOCHECKCRASH) {
                LOG(ERROR) << "Cannot find input friend trees in: " << friendName << "\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "Cannot find input friend trees in: " << friendName << "\n";
            }
        }

        if ((t.AddFriend(friendTree.c_str(), friendName.c_str()) == 0 ) && hasWildcard){
            LOG(WARNING) << "You used wildcards for friend trees, but added zero files from " << friendName << "\n";
        }
    }

    for (const auto& alias : aliases) {
        auto sub_str = alias.find_first_of(":");
        std::string str_alias = alias.substr(0,sub_str);
        std::string str_formula = alias.substr(sub_str+1);
        if (str_alias != str_formula){
            t.SetAlias(str_alias.c_str(),str_formula.c_str());
        }
        else {
            LOG(WARNING) << "Alias and formula must be splitted by colon, e.g. HT:pt1+pt2+pt3\n";
        }
    }
    h->Sumw2();
    TString drawVariable = Form("%s>>h",variable.c_str()), drawWeight = Form("(%s)*(%s)",weight.c_str(),selection.c_str());
    if(Nev>=0) t.Draw(drawVariable, drawWeight, "goff", Nev);
    else       t.Draw(drawVariable, drawWeight, "goff");
    if(TRExFitter::MERGEUNDEROVERFLOW) Common::MergeUnderOverFlow(h);
    return h;
}

//__________________________________________________________________________________
//
TH1D* Common::HistFromNtupleBinArr(const std::string& ntuple,
                                   const std::string& frnd,
                                   const std::string& variable,
                                   int nbin,
                                   double *bins,
                                   const std::string& selection,
                                   const std::string& weight,
                                   const std::vector<std::string>& aliases,
                                   int Nev) {
    TH1D* h = new TH1D("h","h",nbin,bins);
    LOG(VERBOSE) << "  Extracting histogram " << variable << " from  " << ntuple << "  ...\n";
    LOG(VERBOSE) << "      with weight  (" << weight << ")*(" << selection << ")  ...\n";


    bool hasWildcard = false;
    // check whether file actually exists, AccessPathName() returns FALSE if file can be accessed
    // see https://root.cern.ch/root/html602/TSystem.html#TSystem:AccessPathName
    // the string is now "path/to/the/file.root/path/to/the/tree"
    const std::string fileName = ntuple.substr(0,ntuple.find_last_of(".")+5); // find .root and take what is left of it to get the path
    const std::string treeName = ntuple.substr(ntuple.find_last_of(".")+6); // find .root/ and take what is right of it
    if (fileName.find('*') != std::string::npos) hasWildcard = true;
    if (gSystem->AccessPathName(fileName.c_str()) == kTRUE && !hasWildcard ){
        if (TRExFitter::HISTOCHECKCRASH) {
            LOG(WARNING) << "Cannot find input file in: " << fileName << "\n";
            exit(EXIT_FAILURE);
        } else {
            LOG(WARNING) << "Cannot find input file in: " << fileName << "\n";
        }
    }

    TChain t{};
    if (fileName.find('*') == std::string::npos) {
        // does not contain globbing
        t.Add((fileName + "/" + treeName).c_str());
    } else {
        // uses globbing
        Glob(fileName, treeName, t);
    }

    if (frnd != ""){
        LOG(VERBOSE) << "  FriendTree Extracting histogram " << variable << " from  " << frnd << "  ...\n";
        LOG(VERBOSE) << "             with weight  (" << weight << ")*(" << selection << ")  ...\n";
        const std::string friendName = frnd.substr(0,  frnd.find_last_of("/"));
        const std::string friendTree = frnd.substr(frnd.find_last_of("/")+1);
        hasWildcard = false;
        if (friendName.find('*') != std::string::npos) hasWildcard = true;
        if (gSystem->AccessPathName(friendName.c_str()) == kTRUE && !hasWildcard ){
            if (TRExFitter::HISTOCHECKCRASH) {
                LOG(ERROR) << "Cannot find input friend tree file in: " << friendName << "\n";
                exit(EXIT_FAILURE);
            } else {
                LOG(WARNING) << "Cannot find input friend tree file in: " << friendName << "\n";
            }
        }

        if ((t.AddFriend(friendTree.c_str(), friendName.c_str()) == 0 ) && hasWildcard){
            LOG(WARNING) << "You used wildcards, but added zero friend tree files from " << fileName << "\n";
        }
    }

    for (const auto& alias : aliases) {
        auto sub_str = alias.find_first_of(":");
        std::string str_alias = alias.substr(0,sub_str);
        std::string str_formula = alias.substr(sub_str+1);
        if (str_alias != str_formula){
            t.SetAlias(str_alias.c_str(),str_formula.c_str());
        }
        else {
            LOG(WARNING) << "Alias and formula must be splitted by colon, e.g. HT:pt1+pt2+pt3\n";
        }
    }
    h->Sumw2();
    TString drawVariable = Form("%s>>h",variable.c_str()), drawWeight = Form("(%s)*(%s)",weight.c_str(),selection.c_str());
    if(Nev>=0) t.Draw(drawVariable, drawWeight, "goff", Nev);
    else       t.Draw(drawVariable, drawWeight, "goff");
    if(TRExFitter::MERGEUNDEROVERFLOW) Common::MergeUnderOverFlow(h);
    return h;
}

//__________________________________________________________________________________
//
TFile* Common::GetFile(const std::string& fileName) {
    auto it = TRExFitter::TFILEMAP.find(fileName);
    if(it != TRExFitter::TFILEMAP.end()) return it->second.get();
    else {
        std::unique_ptr<TFile> f(TFile::Open(fileName.c_str()));
        TRExFitter::TFILEMAP.insert(std::make_pair(fileName, std::move(f)));
        return TRExFitter::TFILEMAP[fileName].get();
    }
}

//__________________________________________________________________________________
//
std::unique_ptr<TH1> Common::HistFromFile(const std::string& fullName,
                                          const bool& dropBinHistoNeeded,
                                          const bool& checkHistogram,
                                          const Common::FolderStructure* structure
                                          ) {
    const std::string fileName  = fullName.substr(0,fullName.find_last_of(".")+5);
    const std::string histoName = fullName.substr(fullName.find_last_of(".")+6,std::string::npos);
    return HistFromFile(fileName,histoName, dropBinHistoNeeded, checkHistogram, structure);
}

//__________________________________________________________________________________
//
std::unique_ptr<TH1> Common::HistFromFile(const std::string& fileName,
                                          const std::string& histoName,
                                          const bool& dropBinHistoNeeded,
                                          const bool& checkHistogram,
                                          const Common::FolderStructure* structure
                                          ) {
    if(fileName=="") return nullptr;
    if(histoName=="") return nullptr;
    bool hasCustomAsimov = false;
    if (fileName.find("customAsimov") != std::string::npos) hasCustomAsimov = true;
    LOG(VERBOSE) << "  Extracting histogram    " << histoName << "  from file    " << fileName << "    ...\n";
    TFile* f = Common::GetFile(fileName);
    if(!f){
        LOG(ERROR) << "Cannot find input file '" << fileName << "'\n";
        return nullptr;
    }
    const std::string folder = (structure && TRExFitter::USEFOLDERSTRUCTURE) ? structure->region+"/"+structure->sample+"/"+structure->systematics+"/" : "" ;

    auto isHistoAvailable = [f, folder, histoName](const std::string& suffix){

        const bool folderDirIsOkay = f->cd(folder.c_str());
        if (!folderDirIsOkay) {
            LOG(ERROR) << "Cannot find folder: " << folder << "\n";
            exit(EXIT_FAILURE);
        }

        if ( ((histoName.find("_Up") != std::string::npos) || (histoName.find("_Down") != std::string::npos))) {

            std::string histoNameShape;
            if (histoName.find("_Up") != std::string::npos) histoNameShape = ReplaceString(histoName, "_Up", "_Shape_Up");
            else                                            histoNameShape = ReplaceString(histoName, "_Down", "_Shape_Down");

            if (gDirectory->GetListOfKeys()->Contains((histoNameShape).c_str())) {
                if (!gDirectory->GetListOfKeys()->Contains((histoNameShape+suffix).c_str())) {
                    LOG(ERROR) << suffix << " histogram was not found for " << histoName << "\n";
                    LOG(ERROR) << "Make sure you run n/b/h step again with one of the following:\n";
                    LOG(ERROR) << "A) bin dropping activated in config if you are trying to access bin-dropped histogram (_dropBin error).\n";
                    LOG(ERROR) << "B) bin dropping de-activated in config if you are trying to access a complete histogram (_regBin error).\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
        else{
            if (!gDirectory->GetListOfKeys()->Contains((histoName+suffix).c_str())) {
                LOG(ERROR) << suffix << " histogram was not found for " << histoName << " even though bins are being dropped.\n";
                LOG(ERROR) << "Make sure you run n/b/h step again with one of the following:\n";
                LOG(ERROR) << "A) bin dropping activated in config if you are trying to access bin-dropped histogram (_dropBin error).\n";
                LOG(ERROR) << "B) bin dropping de-activated in config if you are trying to access a complete histogram (_regBin error).\n";
                exit(EXIT_FAILURE);
            }
        }
    };

    if (histoName.find("_orig") == std::string::npos && dropBinHistoNeeded && checkHistogram){
        isHistoAvailable(std::string("_dropBin"));
    }
    else if (histoName.find("_orig") == std::string::npos && !dropBinHistoNeeded && checkHistogram){
        isHistoAvailable("_regBin");
    }

    std::unique_ptr<TH1> h(f->Get<TH1>((folder+histoName).c_str()));
    if(!h){
        if (!hasCustomAsimov) {
            LOG(ERROR) << "Cannot find histogram '" << folder << histoName << "' from input file '" << fileName << "'\n";
            if (TRExFitter::USEFOLDERSTRUCTURE) {
                LOG(ERROR) << "If you get this crash reading the histograms produced by TRExFitter, maybe you need to use \"UseFolderStructure: FALSE\" in the Job block if you are reading histograms produced in an older version of TRExFitter.\n";
            }
        }
        else LOG(DEBUG) << "Cannot find histogram '" << folder << histoName << "' from input file '" << fileName << "', but its customAsimov histogram so this should not be a problem\n";
        return nullptr;
    }
    h->SetDirectory(nullptr);
    if(TRExFitter::MERGEUNDEROVERFLOW) Common::MergeUnderOverFlow(h.get());
    return h;
}

//__________________________________________________________________________________
//
std::unique_ptr<TH2> Common::Hist2DFromFile(const std::string& fullName) {
    const std::string fileName  = fullName.substr(0,fullName.find_last_of(".")+5);
    const std::string histoName = fullName.substr(fullName.find_last_of(".")+6,std::string::npos);
    return Hist2DFromFile(fileName,histoName);
}

//__________________________________________________________________________________
//
std::unique_ptr<TH2> Common::Hist2DFromFile(const std::string& fileName,
                                            const std::string& histoName) {
    if(fileName=="") return nullptr;
    if(histoName=="") return nullptr;
    LOG(VERBOSE) << "  Extracting histogram    " << histoName << "  from file    " << fileName + "    ...\n";
    TFile* f = Common::GetFile(fileName);
    if(!f){
        LOG(ERROR) << "Cannot find input file '" << fileName << "'\n";
        return nullptr;
    }
    std::unique_ptr<TH2> h(dynamic_cast<TH2*>(f->Get(histoName.c_str())));
    if(!h){
        LOG(ERROR) << "Cannot find histogram '" << histoName << "' from input file '" << fileName << "'\n";
        return nullptr;
    }
    h->SetDirectory(0);
    return h;
}

//__________________________________________________________________________________
//
void Common::WriteHistToFile(TH1* h,
                             const std::string& fileName,
                             const Common::FolderStructure& structure,
                             const std::string& option) {
    TDirectory *dir = gDirectory;
    std::unique_ptr<TFile> f(TFile::Open(fileName.c_str(),option.c_str()));
    if (!f) {
        LOG(ERROR) << "Cannot open file at: " << fileName << "\n";
        return;
    }

    if (TRExFitter::USEFOLDERSTRUCTURE) {
        TDirectory* currentDir = Common::SetFolderStructure(f.get(), structure);
        currentDir->cd();
    }

    h->Write("",TObject::kOverwrite);
    h->SetDirectory(nullptr);
    f->Close();
    dir->cd();
}

//__________________________________________________________________________________
//
void Common::WriteHistToFile(TH1* h,
                             TFile* f,
                             const Common::FolderStructure& structure) {
    TDirectory *dir = gDirectory;
    f->cd();

    if (TRExFitter::USEFOLDERSTRUCTURE) {
        TDirectory* currentDir = Common::SetFolderStructure(f, structure);
        currentDir->cd();
    }

    h->Write("",TObject::kOverwrite);
    h->SetDirectory(nullptr);
    dir->cd();
}

//__________________________________________________________________________________
//
void Common::MergeUnderOverFlow(TH1* h) {
    int nbins = h->GetNbinsX();
    h->AddBinContent( 1, h->GetBinContent(0) ); // merge first bin with underflow bin
    h->SetBinError(   1, std::hypot( h->GetBinError(1),h->GetBinError(0))); // increase the stat uncertainty as well
    h->AddBinContent( nbins, h->GetBinContent(nbins+1) ); // merge first bin with overflow bin
    h->SetBinError(   nbins, std::hypot(h->GetBinError(nbins),h->GetBinError(nbins+1))); // increase the stat uncertainty as well
    // set under/overflow bins and its errors to 0
    h->SetBinContent( 0, 0. );
    h->SetBinContent( nbins+1, 0. );
    h->SetBinError(   0, 0. );
    h->SetBinError(   nbins+1, 0. );
}

//__________________________________________________________________________________
//
std::vector<std::string> Common::CreatePathsList(std::vector<std::string> paths,
                                                 std::vector<std::string> pathSufs,
                                                 std::vector<std::string> files,
                                                 std::vector<std::string> fileSufs,
                                                 std::vector<std::string> names,
                                                 std::vector<std::string> nameSufs) {
    // turn the empty vectors into vectors containing one "" entry
    if(paths.size()==0) paths.push_back("");
    if(pathSufs.size()==0) pathSufs.push_back("");
    if(files.size()==0) files.push_back("");
    if(fileSufs.size()==0) fileSufs.push_back("");
    if(names.size()==0) names.push_back("");
    if(nameSufs.size()==0) nameSufs.push_back("");
    //
    std::vector<std::string> output;
    for (const auto& ipath : paths) {
        for(const auto& ipathSuf : pathSufs){
            for(const auto& ifile : files){
                for(const auto& ifileSuf : fileSufs){
                    for(const auto& iname : names){
                        for(const auto& inameSuf : nameSufs){
                            std::string fullPath = ipath;
                            fullPath += ipathSuf;
                            fullPath += "/";
                            fullPath += ifile;
                            fullPath += ifileSuf;
                            fullPath += ".root";
                            if(iname !="" || inameSuf!=""){
                                fullPath += "/";
                                fullPath += iname;
                                fullPath += inameSuf;
                            }
                            output.emplace_back( fullPath );
                        }
                    }
                }
            }
        }
    }
    return output;
}

//__________________________________________________________________________________
//
std::vector<std::string> Common::CombinePathSufs(std::vector<std::string> pathSufs,
                                                 std::vector<std::string> newPathSufs,
                                                 const bool isFolded) {
    std::vector<std::string> output;
    if (isFolded) {
        output.emplace_back("");
        return output;
    }
    if(pathSufs.size()==0) pathSufs.emplace_back("");
    if(newPathSufs.size()==0) newPathSufs.emplace_back("");
    for(const auto& i : pathSufs) {
        for(const auto& j : newPathSufs) {
            output.push_back(i+j);
        }
    }
    return output;
}

//__________________________________________________________________________________
//
std::vector<std::string> Common::ToVec(const std::string& s) {
    std::vector<std::string> output;
    output.emplace_back(s);
    return output;
}

//__________________________________________________________________________________
//
std::string Common::ReplaceString(std::string subject,
                                  const std::string& search,
                                  const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}

//__________________________________________________________________________________
//
std::vector< std::pair < std::string,std::vector<double> > > Common::processString(std::string target) {
    size_t pos = 0;
    std::vector<std::pair <std::string,std::vector<double> > > output;
    while((pos = target.find("[",pos)) !=std::string::npos) {
        std::pair <std::string, std::vector<double> > onePair;
        std::vector<double> values;
        double oneValue;
        int length = target.find("]",pos) - pos;
        std::stringstream ss(target.substr(pos+1,length-1));
        while (ss>>oneValue){
            values.push_back(oneValue);
            if (ss.peek() == ','){
                ss.ignore();
            }
        }
        onePair.first = target.substr(0,pos);
        onePair.second = values;
        output.push_back(onePair);
        target.erase(0,pos+length+2);
        pos = 0;
    }
    return output;
}

//__________________________________________________________________________________
// taking into account wildcards on both
bool Common::StringsMatch(const std::string& s1, const std::string& s2){
    if(Common::wildcmp(s1.c_str(),s2.c_str())>0 || Common::wildcmp(s2.c_str(),s1.c_str())>0) return true;
    return false;
}

//__________________________________________________________________________________
// taking into account wildcards on first argument
int Common::wildcmp(const char *wild, const char *string) {
    // Written by Jack Handy - <A href="mailto:jakkhandy@hotmail.com">jakkhandy@hotmail.com</A>
    const char *cp = nullptr;
    const char *mp = nullptr;
    while ((*string) && (*wild != '*')) {
        if ((*wild != *string) && (*wild != '?')) {
            return 0;
        }
        wild++;
        string++;
    }
    while (*string) {
        if (*wild == '*') {
          if (!*++wild) {
              return 1;
          }
          mp = wild;
          cp = string+1;
        } else if ((*wild == *string) || (*wild == '?')) {
            wild++;
            string++;
        } else {
            wild = mp;
            string = cp++;
        }
    }
    while (*wild == '*') {
        wild++;
    }
    return !*wild;
}

//__________________________________________________________________________________
//
int Common::FindInStringVector(const std::vector<std::string>& v,
                               const std::string& s) {
    int idx = -1;
    for(unsigned int i=0;i<v.size();i++){
        std::string s1 = v[i];
        if(Common::StringsMatch(s1,s)){
            idx = (int)i;
            break;
        }
    }
    return idx;
}

//__________________________________________________________________________________
//
int Common::FindInStringVectorOfVectors(const std::vector< std::vector<std::string> >& v,
                                        const std::string& s,
                                        const std::string& ss) {
    int idx = -1;
    for(unsigned int i=0;i<v.size();i++){
        std::string s1 = v[i][0];
        std::string s2 = v[i][1];
        if(Common::StringsMatch(s1,s) && Common::StringsMatch(s2,ss)){
            idx = (int)i;
            break;
        }
    }
    return idx;
}

//__________________________________________________________________________________
//
double Common::GetSeparation(TH1D* S1, TH1D* B1) {
    // taken from TMVA!!!
    std::unique_ptr<TH1> S(static_cast<TH1*>(S1->Clone()));
    std::unique_ptr<TH1> B(static_cast<TH1*>(B1->Clone()));
    Double_t separation = 0;
    if ((S->GetNbinsX() != B->GetNbinsX()) || (S->GetNbinsX() <= 0)) {
        LOG(ERROR) << "Signal and background histograms have different number of bins: " << S->GetNbinsX() << " : " << B->GetNbinsX() << "\n";
    }
    if (S->GetXaxis()->GetXmin() != B->GetXaxis()->GetXmin() ||
            S->GetXaxis()->GetXmax() != B->GetXaxis()->GetXmax() ||
            S->GetXaxis()->GetXmax() <= S->GetXaxis()->GetXmin()) {
        LOG(ERROR) << "Signal and background histograms have different or invalid dimensions:\n";
        LOG(ERROR) << "Signal Xmin: " << S->GetXaxis()->GetXmin() << ", background Xmin " << B->GetXaxis()->GetXmin() << "\n";
        LOG(ERROR) << "Signal Xmax: " << S->GetXaxis()->GetXmax() << ", background Xmax " << B->GetXaxis()->GetXmax() << "\n";
    }
    Int_t nstep     = S->GetNbinsX();
    Double_t intBin = (S->GetXaxis()->GetXmax() - S->GetXaxis()->GetXmin())/nstep;
    Double_t nS     = S->GetSumOfWeights()*intBin;
    Double_t nB     = B->GetSumOfWeights()*intBin;
    if (nS > 0 && nB > 0) {
        for (Int_t bin=0; bin <= nstep + 1; bin++) {
            Double_t s = S->GetBinContent( bin )/Double_t(nS);
            Double_t b = B->GetBinContent( bin )/Double_t(nB);
    if (s + b > 0) separation += 0.5*(s - b)*(s - b)/(s + b);
        }
        separation *= intBin;
    }
    else {
        LOG(ERROR) << "Histograms with zero entries: signal: " << nS << " : background : " << nB << " cannot compute separation\n";
        separation = 0;
    }
    return separation;
}

//__________________________________________________________________________________
//
std::unique_ptr<TH1D> Common::BlindDataHisto(TH1* h_data,
                                             const std::vector<int>& blindedBins) {

    std::unique_ptr<TH1D> h_blind(static_cast<TH1D*>(h_data->Clone("h_blind")));
    h_blind->SetDirectory(nullptr);
    for(int i_bin = 1; i_bin <= h_data->GetNbinsX(); ++i_bin) {
        if(std::find(blindedBins.begin(), blindedBins.end(), i_bin) != blindedBins.end()) {
            LOG(DEBUG) << "Blinding bin n." << i_bin << "\n";
            h_data->SetBinContent(i_bin,0.);
            h_data->SetBinError(i_bin,0.);
            h_blind->SetBinContent(i_bin,1.);
        }
        else{
            h_blind->SetBinContent(i_bin,0.);
        }
    }
    return h_blind;
}

//__________________________________________________________________________________
//
void Common::SmoothHistogramTtres(TH1* h) {
    double origIntegral = h->Integral();

    h->Smooth(2);

    if(h->Integral()!=0){
        h->Scale(origIntegral/h->Integral());
    }
}

//__________________________________________________________________________________
// to smooth a nominal histogram, taking into account the statistical uncertinaty on each bin (note: no empty bins, please!!)
bool Common::SmoothHistogram(TH1* h,
                             double nsigma){
    int nbinsx = h->GetNbinsX();
    double error = 0.;
    double integral = h->IntegralAndError(1,h->GetNbinsX(),error);
    //
    // if not flat, go on with the smoothing
    int Nmax = 5;
    for(int i=0;i<Nmax;i++){
        std::unique_ptr<TH1> h0(static_cast<TH1*>(h->Clone("h0")));
        h->Smooth();
        bool changesApplied = false;
        for(int i_bin=1;i_bin<=nbinsx;i_bin++){
            if( std::abs(h->GetBinContent(i_bin) - h0->GetBinContent(i_bin)) > nsigma*h0->GetBinError(i_bin) ){
                h->SetBinContent(i_bin,h0->GetBinContent(i_bin));
            }
            else{
                changesApplied = true;
            }
            // bring bins < 1e-6 to 1e-06
            if(h->GetBinContent(i_bin)<1e-06) h->SetBinContent(i_bin,1e-06);
        }
        if(!changesApplied) break;
    }

    //
    // try to see if it's consistent with being flat
    //
    // make sure you didn't change the integral
    if(h->Integral()>0){
        h->Scale(integral/h->Integral());
    }
    //
    // fix stat error so that the total stat error is unchanged, and it's distributed among all bins
    for(int i_bin=1;i_bin<=nbinsx;i_bin++){
        double N = integral;
        double E = error;
        double n = h->GetBinContent(i_bin);
        h->SetBinError(i_bin,E*std::sqrt(n)/std::sqrt(N));
    }
    //
    return false; // this is actual behaviour that was implemented previously
}

//__________________________________________________________________________________
//
double Common::CorrectIntegral(TH1* h, double* err) {
    double integral = 0.;
    double error = 0.;
    for( int i_bin=1; i_bin <= h->GetNbinsX(); i_bin++){
        if(h->GetBinContent(i_bin)<0) continue;
        integral += h->GetBinContent(i_bin);
        if(h->GetBinError(i_bin)<=0) continue;
        error += h->GetBinError(i_bin) * h->GetBinError(i_bin);
    }
    if(err!=0) *err = std::sqrt(error);
    return integral;
}

//__________________________________________________________________________________
//
TH1D* Common::MergeHistograms(const std::vector<std::unique_ptr<TH1> >& hVec,bool fixLastBinWidth){
    std::vector<TH1*> vec;
    for (const auto& i : hVec) {
        vec.emplace_back(i.get());
    }
    return Common::MergeHistograms(vec,fixLastBinWidth);
}

//__________________________________________________________________________________
//
TH1D* Common::MergeHistograms(const std::vector<TH1*>& hVec,bool fixLastBinWidth){
    if(hVec.size()==0) return nullptr;
    if(hVec[0]==nullptr) return nullptr;
    // build vector of bin edges
    std::vector<double> binVec;
    binVec.push_back( hVec[0]->GetXaxis()->GetBinLowEdge(1) );
    // define the offset, which will be increased by the last bin UpEdge of a histogram at the end of the loop on its bins
    double offset = 0;
    //
    for(unsigned int i_h=0;i_h<hVec.size();i_h++){
        TH1* h = hVec[i_h];
        for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
            double binUpEdge = h->GetXaxis()->GetBinUpEdge(i_bin);
            // if fixLastBinWidth is TRUE, set bin width of last bin to the next-to-last one
            if(fixLastBinWidth && i_bin==h->GetNbinsX() && i_bin>1){
                double binWidth = h->GetXaxis()->GetBinUpEdge(i_bin-1) - h->GetXaxis()->GetBinLowEdge(i_bin-1);
                binUpEdge = h->GetXaxis()->GetBinLowEdge(i_bin) + binWidth;
            }
            if(i_h==0) binVec.push_back( binUpEdge + offset );
            else       binVec.push_back( binUpEdge - h->GetXaxis()->GetBinLowEdge(1) + offset );
            if(i_bin==h->GetNbinsX()){
                if(i_h==0) offset += binUpEdge;
                else       offset += binUpEdge - h->GetXaxis()->GetBinLowEdge(1);
            }
        }
    }
    int Nbins = binVec.size()-1;
    // create the new histogram
    TH1D* hOut = new TH1D("h_merge","h_merge",Nbins,&binVec[0]);
    hOut->SetTitle(hVec[0]->GetTitle());
    hOut->SetLineColor(hVec[0]->GetLineColor());
    hOut->SetLineStyle(hVec[0]->GetLineStyle());
    hOut->SetLineWidth(hVec[0]->GetLineWidth());
    hOut->SetFillColor(hVec[0]->GetFillColor());
    hOut->SetFillStyle(hVec[0]->GetFillStyle());
    // fill it
    int k_bin = 1;
    for(const auto& h : hVec){
        for(int i_bin=1;i_bin<=h->GetNbinsX();i_bin++){
            hOut->SetBinContent(k_bin,h->GetBinContent(i_bin));
            hOut->SetBinError(k_bin,h->GetBinError(i_bin));
            k_bin ++;
        }
    }
    // return
    return hOut;
}

//___________________________________________________________
//
int Common::ApplyPDGrounding(double &mean,
                             double &error) {
    if (error < 0 ){
        LOG(WARNING) << "Error value is < 0. Not applying rounding.\n";
        return -1;
    }

    int sig = 0;
    const int iterations = Common::ApplyErrorRounding(error,sig);
    if (iterations > 100) { // something went wrong
        LOG(WARNING) << "Problem with applying PDG rounding rules to error.\n";
        return -1;
    }

    // now apply the correct rounding for nominal value
    Common::RoundToSig(mean, iterations);

    // return the number of decimal digits (for later printing avoiding exponent...)
    int decPlaces = iterations;
    if(iterations<0) decPlaces = 0;
    return decPlaces;
}

//___________________________________________________________
// FIXME : still to fix the 100
int Common::ApplyErrorRounding(double& error,
                               int& sig) {
    int iterations = 0;

    if (error == 0) {
        LOG(WARNING) << "Error is zero, you should have a look at this.\n";
        return 0;
    }

    while (error < 100) {
        error*= 10;
        iterations++;
        if (iterations > 15){
            LOG(WARNING) << "Too many iterations in determination of decimal places. Not applying rounding\n";
            return 999;
        }
    }

    while (error >= 1000) {
        error/= 10;
        iterations--;
        if (iterations < -15){
            LOG(WARNING) << "Too many iterations in determination of decimal places. Not applying rounding\n";
            return 999;
        }
    }

    // PDG rounding rules
    if (error >= 100 && error < 355){
        sig = 2;
    } else if (error >= 355 && error < 950) {
        sig = 1;
    } else if (error >= 950) {
        error = 1000;
        sig = 1;
    } else {
        LOG(WARNING) << "3 significant digit are < 100 or > 999. This should not happen.\n";
        return 999;
    }

    // have three significant digits, now round
    // according to the number of decimal places
    error/= IntPow(10, (3-sig));
    error = std::round(error);
    error*= IntPow(10, (3-sig));

    // now we need to get back to original value
    // this is not optimal but should be optimized by compiler
    if(iterations>0) error/= IntPow(10, std::abs(iterations));
    if(iterations<0) error*= IntPow(10, std::abs(iterations));

    // return number of iterations needed minus 2
    // this will be used to match precision of mean to precision
    // of rounded error
    return (iterations - 3 + sig);
}

//___________________________________________________________
//
void Common::RoundToSig(double& value,
                        const int& n){
    if (n == 0) {
        value = std::round(value);
        return;
    }

    if (n > 0) { // will multiply
        value*= IntPow(10,n);
        value = std::round(value);
        value/= IntPow(10,n);
    } else { // will divide
        value/= IntPow(10,std::abs(n));
        value = std::round(value);
        value*= IntPow(10,std::abs(n));
    }
}

//___________________________________________________________
//
std::string Common::KeepSignificantDigits(double value, const int n) {
    if (n < 1) {
        LOG(WARNING) << "Number of significant digits < 1\n";
        return "n/a";
    }
    int iterations(0);
    if (value > 0) {
        if (value > IntPow(10,n)) {
            while (value > IntPow(10,n)) {
                value/=10;
                ++iterations;
            }
            value = std::round(value);
            value*= IntPow(10,iterations);
            return Form("%.f", value);
        } else {
            while (value < IntPow(10,n-1)) {
                value*=10;
                ++iterations;
            }
            value = std::round(value);
            value/= IntPow(10,iterations);
            return Form(("%."+std::to_string(iterations)+"f").c_str(), value);
        }
    } else if (value < 0){
        if (value < -IntPow(10,n)) {
            while (value < -IntPow(10,n)) {
                value/=10;
                ++iterations;
            }
            value = std::round(value);
            value*= IntPow(10,iterations);
            return Form("%.f", value);
        } else {
            while (value > -IntPow(10,n-1)) {
                value*=10;
                ++iterations;
            }
            value = std::round(value);
            value/= IntPow(10,iterations);
            return Form(("%."+std::to_string(iterations)+"f").c_str(), value);
        }
    } else {
        return "0";
    }

    LOG(WARNING) << "This should never be reached\n";
    return "n/a";
}

//___________________________________________________________
//
unsigned int Common::NCharactersInString(const std::string& s,
                                         const char c){
    unsigned int N = 0;
    for(const auto& i_c : s) {
        if(i_c == c) N++;
    }
    return N;
}

//___________________________________________________________
// for the moment just checks the number of parenthesis, but can be expanded
bool Common::CheckExpression(const std::string& s) {
    if(s.find("Alt$")!=std::string::npos){
        return true;
    }
    int nParOpen = Common::NCharactersInString(s,'(');
    int nParClose = Common::NCharactersInString(s,')');
    if(nParOpen!=nParClose) return false;
    // ...
    return true;
}

//___________________________________________________________
//
bool Common::OptionRunsFit(const std::string& opt){
    if (opt.find("f")!=std::string::npos) return true;
    if (opt.find("l")!=std::string::npos) return true;
    if (opt.find("s")!=std::string::npos) return true;
    if (opt.find("r")!=std::string::npos) return true;
    if (opt.find("i")!=std::string::npos) return true;
    if (opt.find("x")!=std::string::npos) return true;
    return false;
}

//___________________________________________________________
//
std::unique_ptr<TH1> Common::GetHistCopyNoError(const TH1* const hist){
    if (hist == nullptr) return nullptr;
    std::unique_ptr<TH1> result(static_cast<TH1*>(hist->Clone()));

    for (int ibin = 0; ibin <= hist->GetNbinsX(); ++ibin){
        result->SetBinError(ibin, 0.);
    }
    result->SetDirectory(nullptr);

    return result;
}

// BW helper functions to pad bin numbers for gamma plots
// replaces them with zero padded versions.  "Gamma Bin 1" -> "Gamma Bin 0001"

std::vector<std::string> Common::mysplit(const std::string& s,
                                         const char delimiter) {
    std::vector<std::string> answer;
    std::string token;

    // this converts a single char into a string
    // the constructor std:string( n, char )
    // produces a string of n copies of char
    std::string localDelim( 1, delimiter );
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        //keep the delimiter for easy reconstruction
        answer.push_back(token+localDelim);
    }

    //remove the trailing delim from the last token
    std::string last = answer.back();
    if( !last.empty() ) last.pop_back();

    answer.pop_back();
    answer.push_back( last );

    return answer;
}

//___________________________________________________________
//
std::string Common::addpad(const std::string& input,
                           const char filler,
                           const unsigned width ) {
    std::stringstream mySS;

    mySS.fill(filler);
    mySS.width(width);

    mySS << input;

    return mySS.str();

}

//___________________________________________________________
//
std::string Common::pad_trail(const std::string& input) {

    std::vector<std::string> words = mysplit( input, ' ' );

    std::string paddedValue = addpad( words.back(), '0', 4 );

    words.pop_back();
    words.push_back( paddedValue );

    std::string answer = std::accumulate( words.begin(), words.end(), std::string("") );

    return answer;
}

//___________________________________________________________
// Helper functions to drop norm or shape part from systematic variations
void Common::DropNorm(TH1* hUp,
                      TH1* hDown,
                      TH1* hNom) {
    double err(0);
    const double intNom = Common::CorrectIntegral(hNom, &err);
    Common::DropNorm(hUp, hDown, intNom);
}

void Common::DropNorm(TH1* hUp,
                      TH1* hDown,
                      const double intNom) {
    double err(0);

    if(hUp!=nullptr){
        const double intUp = Common::CorrectIntegral(hUp, &err);
        if(std::abs(intUp) > 1e-6) hUp->Scale(intNom/intUp);
        else LOG(WARNING) << "Integral of up variation = 0. Cannot drop normalization.\n";
    }
    if(hDown!=nullptr){
        const double intDown = Common::CorrectIntegral(hDown, &err);

        if(std::abs(intDown) > 1e-6) hDown->Scale(intNom/intDown);
        else LOG(WARNING) << "Integral of down variation = 0. Cannot drop normalization.\n";
    }
}

//___________________________________________________________
//
void Common::DropShape(TH1* hUp,
                       TH1* hDown,
                       TH1* hNom){

    double err(0);
    const double intNom = Common::CorrectIntegral(hNom, &err);
    if(std::abs(intNom) < 1e-6) {
        LOG(WARNING) << "Integral of nominal histogram = 0. Cannot drop shape of syst variations.\n";
        return;
    }
    if(hUp!=nullptr){
        const double intUp = Common::CorrectIntegral(hUp, &err);
        if (std::abs(intUp) > 1e-6) {
            Common::SetHistoBinsFromOtherHist(hUp, hNom);
            hUp->Scale(intUp/intNom);
        } else {
            LOG(WARNING) << "Integral of up variation = 0. Cannot drop shape.\n";
        }
    }
    if(hDown!=nullptr){
        const double intDown = Common::CorrectIntegral(hDown, &err);
        if (std::abs(intDown) > 1e-6) {
            Common::SetHistoBinsFromOtherHist(hDown, hNom);
            hDown->Scale(intDown/intNom);
        } else {
            LOG(WARNING) << "Integral of down variation = 0. Cannot drop shape.\n";
        }
    }
}

//___________________________________________________________
//
void Common::ScaleMCstatInHist(TH1* hist,
                               const double scale) {
    if (std::abs(scale-1) < 1e-6) return; // basically scale == 1 but floating precision

    for (int ibin = 1; ibin <= hist->GetNbinsX(); ++ibin) {
        hist->SetBinError(ibin, scale * hist->GetBinError(ibin));
    }
}

//___________________________________________________________
//
void Common::SetHistoBinsFromOtherHist(TH1* toSet,
                                       const TH1* other) {
    if (!toSet) return;
    if (!other) return;

    const int nbins = toSet->GetNbinsX();
    if (other->GetNbinsX() != nbins) {
        LOG(WARNING) << "Bin sizes are different! Skipping\n";
        return;
    }

    for (int ibin = 1; ibin <= nbins; ++ibin) {
        toSet->SetBinContent(ibin, other->GetBinContent(ibin));
        toSet->SetBinError  (ibin, other->GetBinError(ibin));
    }
}

//___________________________________________________________
//
double Common::EffIntegral(const TH1* const h) {
    double integral = 0.;
    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin){
        if(h->GetBinContent(ibin)>=0) integral += h->GetBinContent(ibin);
    }
    return integral;
}

//__________________________________________________________________________________
//
std::vector<int> Common::GetBlindedBins(const Region* reg,
                                        std::shared_ptr<xRooNode> node,
                                        const Common::BlindingType type,
                                        const double threshold) {
    if (reg->fSampleHists.empty()) {
        LOG(WARNING) << "Region: " << reg->fName << " has empty file list, returning empty list of blinded bins\n";
        return {};
    }

    if (!reg->fSampleHists.at(0)) {
        LOG(WARNING) << "Region: " << reg->fName << ", first SampleHist is nullptr, returning empty list of blinded bins\n";
        return {};
    }

    std::string sigSampleList("");
    std::string bkgSampleList("");

    /// stack the samplehists
    for(const auto& isample : reg->fSampleHists) {
        if (isample->fSample->fType==Sample::SampleType::GHOST) continue;
        if (isample->fSample->fType==Sample::SampleType::EFT) continue;
        if (isample->fSample->fType==Sample::SampleType::DATA) continue;
        if (isample->fSample->fType==Sample::SampleType::SIGNAL) {
            sigSampleList += isample->fSample->fName + "_" + reg->fName + "_shapes" + ",";
        } else if (isample->fSample->fType==Sample::SampleType::BACKGROUND) {
            bkgSampleList += isample->fSample->fName + "_" + reg->fName + "_shapes" + ",";
        } else {
            LOG(ERROR) << "Unknown sample type!\n";
            exit(EXIT_FAILURE);
        }
    }

    std::shared_ptr<TH1> sigHist(nullptr);
    std::shared_ptr<TH1> bkgHist(nullptr);

    if (!sigSampleList.empty()) {
        // remove trailing ","
        sigSampleList.resize(sigSampleList.size()-1);
        auto sigSample = node->at("samples")->reduced(sigSampleList);
        sigHist = Common::HistoFromxRooNode(sigSample, reg->fSampleHists.at(0)->fHist.get());
    }
    if (!bkgSampleList.empty()) {
        // remove trailing ","
        bkgSampleList.resize(bkgSampleList.size()-1);
        auto bkgSample = node->at("samples")->reduced(bkgSampleList);
        bkgHist = Common::HistoFromxRooNode(bkgSample, reg->fSampleHists.at(0)->fHist.get());
    }

    return Common::ComputeBlindedBins(sigHist.get(), bkgHist.get(), type, threshold);
}

//__________________________________________________________________________________
//
std::vector<int> Common::ComputeBlindedBins(const TH1* signal,
                                            const TH1* bkg,
                                            const Common::BlindingType type,
                                            const double threshold) {

    std::vector<int> result;
    if (threshold < 0) return result;
    if (!signal) return result;
    std::unique_ptr<TH1> combined(static_cast<TH1*>(signal->Clone()));
    if (bkg) {
        combined->Add(bkg);
    }
    for (int ibin = 1; ibin <= signal->GetNbinsX(); ++ibin) {
        double soverb(-1);
        double soversplusb(-1);
        double soversqrtb(-1);
        double soversqrtsplusb(-1);
        if (bkg) {
            if (bkg->GetBinContent(ibin) > 1e-9) {
                soverb = signal->GetBinContent(ibin)/bkg->GetBinContent(ibin);
                soversqrtb = signal->GetBinContent(ibin)/std::sqrt(bkg->GetBinContent(ibin));
            } else {
                soverb = 99999;
                soversqrtb = 99999;
            }
        }
        if (combined) {
            if (combined->GetBinContent(ibin) > 1e-9) {
                soversplusb = signal->GetBinContent(ibin)/combined->GetBinContent(ibin);
                soversqrtsplusb = signal->GetBinContent(ibin)/std::sqrt(combined->GetBinContent(ibin));
            }
        }
        switch(type) {
            case Common::SOVERB:
                if (soverb > threshold) {
                    result.emplace_back(ibin);
                }
                break;
            case Common::SOVERSPLUSB:
                if (soversplusb > threshold) {
                    result.emplace_back(ibin);
                }
                break;
            case Common::SOVERSQRTB:
                if (soversqrtb > threshold) {
                    result.emplace_back(ibin);
                }
                break;
            case Common::SOVERSQRTSPLUSB:
                if (soversqrtsplusb > threshold) {
                    result.emplace_back(ibin);
                }
                break;
            default:
                LOG(ERROR) << "Unknown blinding type\n";
                exit(EXIT_FAILURE);
        }
    }

    return result;
}

//__________________________________________________________________________________
//
std::unique_ptr<TH1> Common::CombineHistosFromFullPaths(const std::vector<std::string>& paths) {
    std::unique_ptr<TH1> result(nullptr);

    for (const auto& ipath : paths) {
        if (!result) {
            result = Common::HistFromFile(ipath);
        } else {
            std::unique_ptr<TH1> tmp = Common::HistFromFile(ipath);
            if (!tmp) {
                LOG(WARNING) << "Cannot add histogram from: " << ipath << ", skipping\n";
                continue;
            }
            result->Add(tmp.get());
        }
    }
    result->SetDirectory(nullptr);

    return result;
}

//__________________________________________________________________________________
//
std::unique_ptr<TH2> Common::CombineHistos2DFromFullPaths(const std::vector<std::string>& paths) {
    std::unique_ptr<TH2> result(nullptr);

    for (const auto& ipath : paths) {
        if (!result) {
            result = Common::Hist2DFromFile(ipath);
            if (!result) {
                LOG(ERROR) << "Cannot read the input histogram\n";
                exit(EXIT_FAILURE);
            }
        } else {
            std::unique_ptr<TH2> tmp = Common::Hist2DFromFile(ipath);
            if (!tmp) {
                LOG(ERROR) << "Cannot add histogram from: " << ipath << ", skipping\n";
                continue;
            }
            result->Add(tmp.get());
        }
    }

    if (result) result->SetDirectory(nullptr);
    return result;
}

//__________________________________________________________________________________
//
std::unique_ptr<TGraphAsymmErrors> Common::GetRatioBand(const TGraphAsymmErrors* total,
                                                        const TH1D* data) {

    std::unique_ptr<TGraphAsymmErrors> result(static_cast<TGraphAsymmErrors*>(total->Clone()));

    std::unique_ptr<TH1D> up(static_cast<TH1D*>(data->Clone()));
    std::unique_ptr<TH1D> down(static_cast<TH1D*>(data->Clone()));

    for (int ibin = 0; ibin < result->GetN(); ++ibin) {
        up->SetBinContent(ibin+1, result->GetErrorYhigh(ibin));
        up->SetBinError(ibin+1, 0);
        down->SetBinContent(ibin+1, result->GetErrorYlow(ibin));
        down->SetBinError(ibin+1, 0);
    }

    std::unique_ptr<TH1D> ratio_up(static_cast<TH1D*>(up->Clone()));
    std::unique_ptr<TH1D> ratio_down(static_cast<TH1D*>(down->Clone()));
    ratio_up->Divide(data);
    ratio_down->Divide(data);

    for (int ibin = 0; ibin < result->GetN(); ++ibin) {
        double x,y;
        result->GetPoint(ibin, x, y);
        result->SetPoint(ibin, x, 1.);

        result->SetPointEYhigh(ibin, ratio_up->GetBinContent(ibin+1));
        result->SetPointEYlow (ibin, ratio_down->GetBinContent(ibin+1));
    }

    return result;
}

//__________________________________________________________________________________
//
void Common::ScaleByBinWidth(TH1* h) {
    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
        const double width = h->GetBinWidth(ibin);
        h->SetBinContent(ibin,h->GetBinContent(ibin)/width);
        h->SetBinError(  ibin,h->GetBinError(ibin)  /width);
    }
}

//__________________________________________________________________________________
//
void Common::ScaleByBinWidth(TGraphAsymmErrors* g) {
    for(int ibin = 0; ibin < g->GetN(); ++ibin) {
        const double width = g->GetErrorXhigh(ibin) + g->GetErrorXlow(ibin);
        g->SetPoint(      ibin,g->GetX()[ibin], g->GetY()[ibin]/width);
        g->SetPointEYhigh(ibin,g->GetErrorYhigh(ibin)/width);
        g->SetPointEYlow( ibin,g->GetErrorYlow(ibin) /width);
    }
}

//__________________________________________________________________________________
//
void Common::ScaleByConst(TGraphAsymmErrors* g, const double scale) {
    for(int ibin = 0; ibin < g->GetN(); ++ibin) {
        g->SetPoint(      ibin,g->GetX()[ibin], g->GetY()[ibin]*scale);
        g->SetPointEYhigh(ibin,g->GetErrorYhigh(ibin)*scale);
        g->SetPointEYlow( ibin,g->GetErrorYlow(ibin) *scale);
    }
}
//__________________________________________________________________________________
//
std::string Common::IntToFixLenStr(int i,int n){
    std::stringstream ss;
    ss << std::setw(n) << std::setfill('0') << i;
    std::string s = ss.str();
    return s;
}

//__________________________________________________________________________________
//
bool Common::StringToBoolean(std::string param) {
    std::transform(param.begin(), param.end(), param.begin(), ::toupper);

    if (param == "TRUE") return true;

    return false;
}

//__________________________________________________________________________________
//
std::vector<std::string> Common::GetFilesMatchingString(const std::string& folder,
                                                        const std::string& key,
                                                        const std::string& key2) {

    std::vector<std::string> result;

    for (const auto& ifile : fs::directory_iterator(folder)) {
        if (ifile.path().u8string().find(key) != std::string::npos &&
            ifile.path().u8string().find(key2) != std::string::npos) {
            result.emplace_back(ifile.path().u8string());
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void Common::MergeTxTFiles(const std::vector<std::string>& input, const std::string& out) {

    std::ofstream outFile;
    outFile.open(out.c_str(), std::ios::trunc);
    if (!outFile.is_open() || !outFile.good()) {
        LOG(WARNING) << "Cannot open output file: " << out << "\n";
        return;
    }

    for (const auto& ifile : input) {
        std::ifstream in(ifile.c_str());
        if (!in.is_open() || !in.good()) {
            LOG(WARNING) << "Cannot open file: " << ifile << "\n";
            continue;
        }

        std::string line;
        while(std::getline(in, line)) {
            outFile << line << "\n";
        }
        in.close();
    }

    outFile.close();
}
//__________________________________________________________________________________
//
std::string Common::CheckName(const std::string& name) {
    const std::string result = RemoveQuotes(name);
    const bool hasWhiteSpace = (std::find_if(result.begin(), result.end(), isspace) != result.end());
    if(((result.size() > 0) && (std::isdigit(result.at(0)))) || hasWhiteSpace) {
        LOG(ERROR) << "Failed to browse name: " << name << "\n";
        LOG(ERROR) << "           Either a number has been detected at the beginning or a whitespace character has been found\n";
        LOG(ERROR) << "           This can lead to unexpected behaviour in HistFactory. Please change the name.\n";
        exit(EXIT_FAILURE);
    } else {
        return result;
    }
}

//__________________________________________________________________________________
//
void Common::CheckTreeName(const std::string& name) {
    bool isOkay(true);

    if (name.find('.') != std::string::npos) {
        isOkay = false;
    }

    if (!isOkay) {
        LOG(ERROR) << "Tree name " << name << " used for ntuple processing contains characters that are not allowed: \".\"\n";
        exit(EXIT_FAILURE);
    }
}

//__________________________________________________________________________________
// Removes '"'
std::string Common::RemoveQuotes(const std::string& s){
    if(s=="") return "";
    std::string ss = s;
    replace(ss.begin(), ss.end(), '"', ' ');
    return Common::RemoveSpaces(ss);
}

//__________________________________________________________________________________
// Removes leading and trailing white spaces
std::string Common::RemoveSpaces(const std::string& s){
    if(s=="") return "";
    std::string ss = s;
    if (ss.find_first_not_of(' ')>=std::string::npos){
        ss = "";
    }
    else if (ss.find_first_not_of(' ')>0){
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' '));
    }
    else{
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' ')+1);
    }
    return ss;
}

//__________________________________________________________________________________
// Removes everything after '%' or '#', but only if not inside quotation marks!!
std::string Common::RemoveComments(const std::string& s){
    if(s=="") return "";
    std::string ss = "";
    bool insideQuotes = false;
    for(const auto& i : s) {
        if(i == '"'){
            if(!insideQuotes) insideQuotes = true;
            else              insideQuotes = false;
        }
        if((i == '%' || i == '#') && !insideQuotes) break;
        ss += i;
    }
    return Common::RemoveSpaces(ss);
}

//__________________________________________________________________________________
//
std::vector<std::string> Common::Vectorize(const std::string& s, char c, bool removeQuotes, bool removeComments) {
    std::vector<std::string> v;
    std::string ss = s;

    if(removeComments) ss=Common::RemoveComments(s);

    if(ss==""){
        v.emplace_back("");
        return v;
    }
    std::string t;
    bool insideQuotes = false;
    for(const auto& i : ss) {
        if(!insideQuotes && i == c){
            if(removeQuotes) v.emplace_back(Common::RemoveQuotes(t));
            else             v.emplace_back(Common::RemoveSpaces(t));
            t = "";
        }
        else if(!insideQuotes && i == '"'){
            insideQuotes = true;
            t += i;
        }
        else if(insideQuotes && i == '"'){
            insideQuotes = false;
            t += i;
        }
        else{
            if (insideQuotes) {
                t += i;
            } else {
                if (removeComments) {
                    if (!std::isspace(i)) {
                        t += i;
                    }
                } else {
                    t += i;
                }
            }
        }
    }
    if(Common::RemoveQuotes(t)!=""){
        if(removeQuotes) v.emplace_back(Common::RemoveQuotes(t));
        else             v.emplace_back(Common::RemoveSpaces(t));
    }
    return v;
}

//__________________________________________________________________________________
//
double Common::IntPow(const double value, const int exp) {
    if (exp == 0) return 1.;
    if (exp == 1) return value;
    if (exp == -1) return 1./value;
    if (value == 0) return 0.;

    const int absexp = std::abs(exp);
    double result = value;

    for (int i = 1; i < absexp; ++i) {
        result*= value;
    }

    return exp > 0 ? result : 1./result;
}

//___________________________________________________________
//
std::vector<double> Common::CalculateShapeFactorReparametrization(const std::shared_ptr<ShapeFactor>& sf,
                                                                  const std::vector<std::shared_ptr<NormFactor> >& nfs,
                                                                  const int ibin,
                                                                  const bool isPostFit,
                                                                  const FitResults* fitResult) {

    std::vector<double> result;

    std::string formula = sf->fExpression.at(ibin).first;
    LOG(VERBOSE) << "ShapeFactor: " << sf->fName << ", bin: " << ibin <<", formula: " << formula << "\n";
    const auto& processedString = Common::processString(sf->fExpression.at(ibin).second);
    for (std::size_t i = 0; i < processedString.size(); ++i) {
        formula = Common::ReplaceString(formula, processedString.at(i).first, "x["+std::to_string(i)+"]");
    }

    TFormula func("func",formula.c_str());
    std::vector<double> params;
    std::vector<double> paramsUp;
    std::vector<double> paramsDown;
    if (isPostFit) {
        if (!fitResult) {
            LOG(ERROR) << "Is post-fit but FitResults object is nullptr\n";
            exit(EXIT_FAILURE);
        }

        for (const auto& i : processedString) {
            params.emplace_back(fitResult->GetNuisParValue(i.first));
            const double up   = fitResult->GetNuisParErrUp(i.first);
            const double down = fitResult->GetNuisParErrDown(i.first);
            paramsUp.emplace_back(up);
            paramsDown.emplace_back(down);
        }

    } else {
        for (const auto& i : processedString) {
            auto itr = std::find_if(nfs.begin(), nfs.end(), [&i](const auto& nf){return i.first == nf->fName;});
            if (itr == nfs.end()) {
                LOG(ERROR) << "Cannot find NF:" << i.first << "\n";
                exit(EXIT_FAILURE);
            }
            params.emplace_back((*itr)->GetNominal());
        }
    }

    LOG(VERBOSE) << "Nominal values for the parameters:\n";
    for (const auto ival : params) {
        LOG(VERBOSE) << "   " << ival << "\n";
    }
    result.emplace_back(func.EvalPar(params.data(), nullptr));
    if (isPostFit) {
        const double upEffect   = func.EvalPar(paramsUp.data(), nullptr);
        const double downEffect = func.EvalPar(paramsDown.data(), nullptr);

        const double finalUp   = (upEffect > downEffect) ? upEffect : downEffect;
        const double finalDown = (downEffect < upEffect) ? downEffect : upEffect;

        result.emplace_back(finalUp);
        result.emplace_back(finalDown);
        LOG(VERBOSE) << "Nominal values for up parameters:\n";
        for (const auto ival : paramsUp) {
            LOG(VERBOSE) << "   " << ival << "\n";
        }
        LOG(VERBOSE) << "Nominal values for down parameters:\n";
        for (const auto ival : paramsDown) {
            LOG(VERBOSE) << "   " << ival << "\n";
        }
    }
    LOG(VERBOSE) << "Recalculated central value: " << result.at(0) << "\n";
    if (isPostFit) {
        LOG(VERBOSE) << "Recalculated up value: " << result.at(1) << "\n";
        LOG(VERBOSE) << "Recalculated down value: " << result.at(2) << "\n";
    }

    return result;
}

//___________________________________________________________
//
void Common::SetGraphToHist(TGraphAsymmErrors* graph, TH1* hist) {
    if (!graph || !hist) {
        LOG(ERROR) << "Nullptr passed\n";
        return;
    }

    const int bins = graph->GetN();
    if (bins != hist->GetNbinsX()) {
        LOG(ERROR) << "Incompatible bins between the histogram and graph\n";
        return;
    }

    for (int ibin = 0; ibin < bins; ++ibin) {
        double x, y;
        graph->GetPoint(ibin, x, y);
        y = hist->GetBinContent(ibin + 1);
        graph->SetPoint(ibin, x, y);
    }
}

//___________________________________________________________
//
void Common::SetGraphToZero(TGraphAsymmErrors* graph) {
    if (!graph) {
        LOG(ERROR) << "Nullptr passed\n";
        return;
    }

    const int bins = graph->GetN();

    for (int ibin = 0; ibin < bins; ++ibin) {
        double x, y;
        graph->GetPoint(ibin, x, y);
        y = 0.;
        graph->SetPoint(ibin, x, y);
    }
}

//___________________________________________________________
//
TDirectory* Common::SetFolderStructure(TFile* file,
                                       const Common::FolderStructure& structure) {

    if (!file) {
        LOG(ERROR) << "File is nullptr\n";
        exit(EXIT_FAILURE);
    }

    // add region folder
    if (!file->Get(structure.region.c_str())) {
        file->cd();
        gDirectory->mkdir(structure.region.c_str());
    }

    // add sample folder
    if (!file->Get((structure.region + "/" + structure.sample).c_str())) {
        file->cd(structure.region.c_str());
        gDirectory->mkdir(structure.sample.c_str());
    }

    // add systematics folder
    if (!file->Get((structure.region + "/" + structure.sample + "/" + structure.systematics).c_str())) {
        file->cd((structure.region + "/" + structure.sample).c_str());
        gDirectory->mkdir(structure.systematics.c_str());
    }

    return static_cast<TDirectory*>(file->Get((structure.region + "/" + structure.sample + "/" + structure.systematics).c_str()));
}

//_______________________________________________________________________________________
//
void Common::GetUnfoldingResult(UnfoldingResult& result,
                                const Unfolding* unfolding,
                                const TH1* truth,
                                const FitResults* fitResults,
                                const std::string& wsFileName,
                                const std::string& wsName,
                                const std::string& frFileName,
                                const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                                const std::map<std::string, double>& statErrorMap,
                                const bool statOnly) {

    // - get workspace
    std::unique_ptr<TFile> wsFile(TFile::Open(wsFileName.c_str(),"read"));
    if (!wsFile) {
        LOG(ERROR) << "Cannot open file with WS\n";
        exit(EXIT_FAILURE);
    }

    // - get fit-result
    std::unique_ptr<TFile> frFile(TFile::Open(frFileName.c_str(),"read"));
    if (!frFile) {
        LOG(ERROR) << "Cannot read the fit results file\n";
        return;
    }
    std::unique_ptr<RooFitResult> fr(nullptr);

    std::vector<std::string> freeParameters;
    for (int i = 0; i < unfolding->fNumberUnfoldingTruthBins; ++i) {
        const std::string name = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
        if(unfolding->fUnfoldNormXSec && i==unfolding->fUnfoldNormXSecBinN-1) continue;
        freeParameters.emplace_back(name);
    }

    for (const auto& inf : normFactors) {
        const std::string name = inf->fName;
        auto itr = std::find(freeParameters.begin(), freeParameters.end(), name);
        if (itr != freeParameters.end()) continue;
        freeParameters.emplace_back(name);
    }

    // pass the fit results to the tool
    for (int i = 0; i < unfolding->fNumberUnfoldingTruthBins; ++i) {
        const std::string name = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
        // if it's last bin and getting norm xsec:
        if(unfolding->fUnfoldNormXSec && i==unfolding->fUnfoldNormXSecBinN-1){
            // calculate propagated value and error in 2 ways:
            // 1. by hand:
            double num = 1.;
            double den = 1.;
            const double N = truth->Integral();
            const double Ni = truth->GetBinContent(i+1);
            for (int j = 0; j < unfolding->fNumberUnfoldingTruthBins; ++j){
                if(j==unfolding->fUnfoldNormXSecBinN-1) continue;
                const double Nj = truth->GetBinContent(j+1);
                num -= fitResults->GetNuisParValue(unfolding->fName + "_Bin_" + Common::IntToFixLenStr(j+1) + "_mu") * Nj/N;
                den -= Nj/N;
            }
            const double mean0 = num/den;
            double error0 = 0.;
            for (int j = 0; j < unfolding->fNumberUnfoldingTruthBins; ++j){
                if(j==unfolding->fUnfoldNormXSecBinN-1) continue;
                const std::string par_j = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(j+1) + "_mu";
                const double Nj = truth->GetBinContent(j+1);
                for (int k = 0; k < unfolding->fNumberUnfoldingTruthBins; ++k){
                    if(k==unfolding->fUnfoldNormXSecBinN-1) continue;
                    const std::string par_k = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(k+1) + "_mu";
                    const double Nk = truth->GetBinContent(k+1);
                    double rho = fitResults->GetCorrelationMatrix()->GetCorrelation(par_j,par_k);
                    double err_j(0);
                    double err_k(0);
                    if (statOnly) {
                        auto itr_j = statErrorMap.find(par_j);
                        auto itr_k = statErrorMap.find(par_k);
                        if (itr_j == statErrorMap.end() || itr_k == statErrorMap.end()) {
                            LOG(ERROR) << "Cannot read parameters from the stat only map\n";
                        } else {
                            err_j = std::abs(itr_j->second);
                            err_k = std::abs(itr_k->second);
                        }
                    } else {
                        err_j = 0.5*std::abs(fitResults->GetNuisParErrUp(par_j) - fitResults->GetNuisParErrDown(par_j));
                        err_k = 0.5*std::abs(fitResults->GetNuisParErrUp(par_k) - fitResults->GetNuisParErrDown(par_k));
                    }
                    error0 += rho*(Nj/Ni)*(Nk/Ni)*err_j*err_k;
                }
            }
            error0 = std::sqrt(error0);
            const double up0   = mean0 + error0;
            const double down0 = mean0 - error0;
            // 2. within RooFit:
            RooWorkspace* ws = dynamic_cast<RooWorkspace*>(wsFile->Get(wsName.c_str())); // to avoid a crash, had to re-get the ws from the TFile
            if (!ws) {
                LOG(ERROR) << "Cannot read the WS\n";
                exit(EXIT_FAILURE);
            }
            const std::vector<double> params = FitUtils::CalculateExpressionRoofit(ws, frFile.get(), name, statOnly, freeParameters);
            // Choose which one to take
            if(params.size() != 4){
                // if problems in getting roofit fit result, just use the by-hand one:
                LOG(WARNING) << "No suitable RooFitResult or norm factor in WS found: taking by-hand calculation.\n";
                result.AddFitValue(mean0, up0, down0);
            }
            else{
                // otherwise, check compatibility:
                // - for mean (1% threshold)
                if(std::abs(params.at(0)-mean0)/(0.5*(params.at(0)+mean0))>0.01) LOG(WARNING) << "Propagation of nominal unfolding result to non-independent bin giving incompatible results: " << mean0 << " vs. " << params.at(0) << ". Note that the latter will be used.\n";
                // - for error (1% threshold)
                if(std::abs(params.at(3)-error0)/(0.5*(params.at(3)+error0))>0.01) LOG(WARNING) << "Propagation of unfolding error to non-independent bin giving incompatible results: " << error0 << " vs. " << params.at(3) << ". Note that the latter will be used.\n";
                // in any case use the roofit one
                result.AddFitValue(params.at(0), params.at(1), params.at(2));
            }
        }
        else{
            auto it = std::find_if(normFactors.begin(), normFactors.end(), [&name](std::shared_ptr<NormFactor> nf){return nf->fName == name;});
            if (it == normFactors.end()) {
                LOG(ERROR) << "Unexpected error in NormFactor name\n";
                exit(EXIT_FAILURE);
            }

            // Does not use reparametrisation
            if ((*it)->fExpression.first == "") {
                const double mean = fitResults->GetNuisParValue(name);
                double up(0);
                double down(0);
                if (statOnly) {
                    auto itr = statErrorMap.find(name);
                    if (itr == statErrorMap.end()) {
                        LOG(ERROR) << "Cannot read parameters from the stat only map\n";
                    } else {
                        up   = mean + std::abs(itr->second);
                        down = mean - std::abs(itr->second);
                    }
                } else {
                    up   = mean + fitResults->GetNuisParErrUp(name);
                    down = mean + fitResults->GetNuisParErrDown(name);
                }
                result.AddFitValue(mean, up, down);
            } else {
                const std::string binName = unfolding->fName + "_Bin_" + Common::IntToFixLenStr(i+1) + "_mu";
                RooWorkspace* ws = dynamic_cast<RooWorkspace*>(wsFile->Get(wsName.c_str())); // to avoid a crash, had to re-get the ws from the TFile
                if (!ws) {
                    LOG(ERROR) << "Cannot read the WS\n";
                    exit(EXIT_FAILURE);
                }
                const std::vector<double> params = FitUtils::CalculateExpressionRoofit(ws, frFile.get(), binName, statOnly, freeParameters);
                if (params.size() != 4) {
                    LOG(ERROR) << "Params size is not 4\n";
                    exit(EXIT_FAILURE);
                }
                result.AddFitValue(params.at(0), params.at(1), params.at(2));
            }
        }
    }
    wsFile->Close();
    frFile->Close();
}

//___________________________________________________________
//
std::string Common::ReplaceAllCharacters(std::string string, const char from, const char to) {
    std::replace(string.begin(), string.end(), from, to);
    return string;
}

//___________________________________________________________
//
std::string Common::ReplaceFolderName(std::string name, const std::string& newFolder) {
    auto it = name.find("/");
    if (it == std::string::npos) {
        LOG(WARNING) << "String: " << name << " does not contain \"/\", new folder name will be prepended\n";
        return newFolder + "/" + name;
    } else {
        return name.replace(0, it, newFolder);
    }
}

//___________________________________________________________
//
std::string Common::join_strings(const std::vector<std::string> &strings_to_join, const std::string& separator) {
    std::string result = "";
    for (unsigned int i = 0; i < strings_to_join.size(); i++)   {
        result = result + strings_to_join[i];
        if (i+1 != strings_to_join.size())  {
           result = result + separator;
        }
    }
    return result;
}

//___________________________________________________________
//
std::vector<std::string> Common::SplitString(std::string input_string, const std::string& separator)    {
    std::vector<std::string> result;
    size_t pos = 0;
    while ((pos = input_string.find(separator)) != std::string::npos) {
        result.push_back(input_string.substr(0, pos));
        input_string.erase(0, pos + separator.length());
    }
    if (input_string.length() > 0) result.push_back(input_string);
    return result;
}

//___________________________________________________________
//
void Common::StripString(std::string *input_string, const std::string &chars_to_remove)    {
    input_string->erase(0,input_string->find_first_not_of(chars_to_remove));
    input_string->erase(input_string->find_last_not_of(chars_to_remove)+1);
}

//___________________________________________________________
//
std::vector<std::string> Common::SplitAndStripString(const std::string &input_string,const std::string& separator) {
    std::vector<std::string> result = SplitString(input_string, separator);
    for (std::string &x : result)    {
        StripString(&x, " \n\t\r");
    }
    return result;
}

//___________________________________________________________
//
bool Common::StartsWith(const std::string &main_string, const std::string &prefix)    {
    const unsigned int prefix_lenght = prefix.length();
    const unsigned int string_lenght = main_string.length();
    if (prefix_lenght > string_lenght) return false;
    return prefix == main_string.substr(0, prefix_lenght);
}

//___________________________________________________________
//
std::shared_ptr<TH1> Common::HistoFromxRooNode(const xRooNode& node, const TH1* hist) {
    std::shared_ptr<TH1> result(static_cast<TH1*>(hist->Clone()));

    const auto bins = node.GetBinContents();

    for (int ibin = 1; ibin <= hist->GetNbinsX(); ++ibin) {
        result->SetBinContent(ibin, bins.at(ibin-1));
        result->SetBinError(ibin, 0.);
    }

    result->SetDirectory(nullptr);

    return result;
}

//___________________________________________________________
//
std::shared_ptr<TH1> Common::HistoFromxRooNode(const std::shared_ptr<xRooNode>& node, const TH1* hist) {
    return Common::HistoFromxRooNode(*node, hist);
}

//___________________________________________________________
//
std::unique_ptr<TGraphAsymmErrors> Common::ErrorFromxRooNode(const xRooNode& node,
                                                             const TH1* hist,
                                                             const FitResults* fr,
                                                             const int nToys) {
    auto error = std::make_unique<TGraphAsymmErrors>(hist);

    const int nbins = hist->GetNbinsX();

    const std::vector<double> contents = node.GetBinContents();
    std::vector<double> errors;
    std::vector<double> errorsUp;
    std::vector<double> errorsDown;
    if (fr && nToys > 0) {
        errorsUp   = node.GetBinErrorsHi(1, 0, *(fr->GetRooFitResult()), nToys);
        errorsDown = node.GetBinErrorsLo(1, 0, *(fr->GetRooFitResult()), nToys);
    } else {
        errors = node.GetBinErrors();
    }

    for (int ibin = 0; ibin < nbins; ++ibin) {
        const double central = contents.at(ibin);
        error->SetPointY(ibin, central);
        const double ex = error->GetErrorX(ibin);
        if (fr && nToys > 0) {
            const double uncUp = errorsUp.at(ibin);
            const double uncDown = errorsDown.at(ibin);
            error->SetPointError(ibin, ex, ex, std::abs(uncDown), std::abs(uncUp));
        } else {
            const double unc = errors.at(ibin);
            error->SetPointError(ibin, ex, ex, unc, unc);
        }
    }
    return error;
}

//___________________________________________________________
//
xRooNode Common::ReadWSandFitResults(const std::string& wsPath, const std::string& frPath) {
    xRooNode node(wsPath.c_str());

    if (frPath != "") {
        std::unique_ptr<TFile> frFile(TFile::Open(frPath.c_str(), "READ"));
        if (!frFile) {
            LOG(ERROR) << "Cannot open the file needed for RooFitResult\n";
            exit(EXIT_FAILURE);
        }
        std::unique_ptr<RooFitResult> fr = Common::RooFitResultFromFile(frFile.get());
        node.SetFitResult(*fr);
    }

    return node;
}

//___________________________________________________________
//
std::string Common::NonSignalSampleList(const std::vector<std::shared_ptr<SampleHist> >& samples, const std::string& region) {
    std::string bkgSamples("");
    for (const auto& isample : samples) {
        if (isample->fSample->fType == Sample::SampleType::SIGNAL) continue;
        if (isample->fSample->fType == Sample::SampleType::DATA) continue;
        const std::string name = isample->fSample->fName + "_" + region + "_shapes";
        bkgSamples += bkgSamples.empty() ? name : "," + name;
    }

    return bkgSamples;
}

//___________________________________________________________
//
double Common::GetNormalisationComponent(const TH1* syst, const TH1* nominal) {
    return (Common::EffIntegral(syst) - Common::EffIntegral(nominal) / Common::EffIntegral(nominal));
}

//___________________________________________________________
//
double Common::GCC(const RooFitResult* fr, const std::vector<std::string>& pois) {
    if (pois.size() < 2) {
        LOG(WARNING) << "Less than 2 POIs provided. Returning 0.\n";
        return 0.;
    }

    if (!fr) {
        LOG(WARNING) << "FitResult is nullptr. Returning 0.\n";
        return 0.;
    }

    std::vector<std::size_t> poiIndices;

    std::size_t index(0);
    for (const auto& par : fr->floatParsFinal()) {
        const std::string name = par->GetName();
        if (std::find(pois.begin(), pois.end(), name) != pois.end()) {
            poiIndices.emplace_back(index);
        }
        ++index;
    }

    if (poiIndices.size() != pois.size()) {
        LOG(WARNING) << "Could not find all POIs in the results. Returning 0.\n";
        return 0.;
    }

    const auto covMatrix = fr->covarianceMatrix();

    TMatrixDSym cov(pois.size());
    for (std::size_t i = 0; i < poiIndices.size(); ++i) {
        for (std::size_t j = 0; j < poiIndices.size(); ++j) {
            const double element = covMatrix(poiIndices.at(i), poiIndices.at(j));
            cov(i,j) = element;
        }
    }

    // calculate GCC
    TMatrixDSym covInverted(cov);
    covInverted.Invert();

    double gcc(0.);
    for (std::size_t i = 0; i < pois.size(); ++i) {
        gcc += std::sqrt(1 - 1./(cov(i,i) * covInverted(i,i)));
    }

    return gcc/pois.size();
}

//___________________________________________________________
//
void Common::ProcessLimitOutput(const xRooNLLVar::xRooHypoSpace& hs,
                                const xRooNLLVar::xRooHypoSpace* injectedHs,
                                const std::string& name,
                                const bool isBlind,
                                const std::string& outputFilePath,
                                const std::string& paramName,
                                float paramValue,
                                bool asymptotic) {
    const auto expected      = hs.limit("cls",0);
    const auto expected1up   = hs.limit("cls",1);
    const auto expected1down = hs.limit("cls",-1);
    const auto expected2up   = hs.limit("cls",2);
    const auto expected2down = hs.limit("cls",-2);
    const auto observed      = hs.limit("cls");
    Float_t injectedLimit    = -1;

    bool isNaN(false);
    if (std::isnan(expected.value())) isNaN = true;
    if (std::isinf(expected.value())) isNaN = true;
    if (asymptotic)
        LOG(INFO) << "Printing results of the asymptotic limit estimate " << name << "\n";
    else
        LOG(INFO) << "Printing results of the toy-based limit estimate " << name << "\n";
    LOG(INFO) << "Expected limit (median): " << expected.value() << ", worst case error: " << expected.error() << "\n";
    LOG(INFO) << "Expected limit (-1 sig): " << expected1down.value() << ", worst case error: " << expected1down.error() << "\n";
    LOG(INFO) << "Expected limit ( 1 sig): " << expected1up.value() << ", worst case error: " << expected1up.error() << "\n";
    LOG(INFO) << "Expected limit (-2 sig): " << expected2down.value() << ", worst case error: " << expected2down.error() << "\n";
    LOG(INFO) << "Expected limit ( 2 sig): " << expected2up.value() << ", worst case error: " << expected2up.error() << "\n";
    if (injectedHs) {
        const auto injected = injectedHs->limit("cls");
        injectedLimit = injected.value();
        if (std::isnan(injected.value()) || std::isinf(injected.value())) {
            LOG(WARNING) << "NaN values identified in the injected hypospace, printing more information\n";
            injectedHs->Print();
        }
        LOG(INFO) << "Injected limit         : " << injected.value() << ", worst case error: " << injected.error() << "\n";
    }
    if (!isBlind) {
        LOG(INFO) << "Observed limit         : " << observed.value() << ", worst case error: " << observed.error() << "\n";
        if (std::isnan(observed.value()) || std::isinf(observed.value())) isNaN = true;
    }

    if (isNaN) {
        LOG(WARNING) << "NaN values identified, printing more information\n";
        hs.Print();
    }

    std::unique_ptr<TFile> out(TFile::Open(outputFilePath.c_str(), "RECREATE"));
    if (!out) {
        LOG(ERROR) << "Cannot open the output file\n";
        return;
    }

    TTree *tree = new TTree("stats","Limit results");
    Float_t upperLimit(0);
    Float_t expectedLimit(0);
    Float_t expectedLimit_plus1(0);
    Float_t expectedLimit_plus2(0);
    Float_t expectedLimit_minus1(9);
    Float_t expectedLimit_minus2(0);

    tree->Branch("obs_upperlimit", &upperLimit);
    tree->Branch("exp_upperlimit", &expectedLimit);
    tree->Branch("exp_upperlimit_plus1", &expectedLimit_plus1);
    tree->Branch("exp_upperlimit_plus2", &expectedLimit_plus2);
    tree->Branch("exp_upperlimit_minus1", &expectedLimit_minus1);
    tree->Branch("exp_upperlimit_minus2", &expectedLimit_minus2);
    tree->Branch("inj_upperlimit", &injectedLimit);
    tree->Branch(paramName.c_str(), &paramValue);

    upperLimit           = isBlind ? -1 : observed.value();
    expectedLimit        = expected.value();
    expectedLimit_plus1  = expected1up.value();
    expectedLimit_plus2  = expected2up.value();
    expectedLimit_minus1  = expected1down.value();
    expectedLimit_minus2  = expected2down.value();
    tree->Fill();

    tree->SetDirectory(out.get());
    out->Write();
    out->Close();
}

//___________________________________________________________
//
bool Common::ProcessSignificanceOutput(const std::pair<double, double>& sigObserved,
                                       const std::pair<double, double>& sigExpected,
                                       const std::pair<double, double>& sigInjected,
                                       const std::string& name,
                                       const bool isBlind,
                                       const std::string& outputFilePath) {

    const double observed    = RooStats::PValueToSignificance(sigObserved.first);
    const double observedErr = RooStats::PValueToSignificance(sigObserved.second + sigObserved.first) - observed;
    const double expected    = RooStats::PValueToSignificance(sigExpected.first);
    const double expectedErr = RooStats::PValueToSignificance(sigExpected.second + sigExpected.first) - expected;
    const double injected    = RooStats::PValueToSignificance(sigInjected.first);
    const double injectedErr = RooStats::PValueToSignificance(sigInjected.second + sigInjected.first) - injected;

    bool isNaN(false);
    if (std::isnan(observed)) isNaN = true;
    if (std::isnan(expected)) isNaN = true;
    if (std::isnan(injected)) isNaN = true;
    if (std::isinf(observed)) isNaN = true;
    if (std::isinf(expected)) isNaN = true;
    if (std::isinf(injected)) isNaN = true;

    LOG(INFO) << "Printing results of the asymptotic significance estimate " << name << "\n";
    LOG(INFO) << "Expected p-value mu = 1 (median)     : " << sigExpected.first << ", error: " << sigExpected.second << "\n";
    if (sigInjected.first > -90) {
        LOG(INFO) << "Injected p-value (median)            : " << sigInjected.first << ", error: " << sigInjected.second << "\n";
    }
    if (!isBlind) {
        LOG(INFO) << "Observed p-value (median)            : " << sigObserved.first << ", error: " << sigObserved.second << "\n";
    }
    LOG(INFO) << "\n";
    LOG(INFO) << "Expected significance mu = 1 (median): " << expected << ", error: " << expectedErr << "\n";
    if (sigInjected.first > -90) {
        LOG(INFO) << "Injected significance (median)       : " << injected << ", error: " << injectedErr << "\n";
    }
    if (!isBlind) {
        LOG(INFO) << "Observed significance (median)       : " << observed << ", error: " << observedErr << "\n";
    }

    std::unique_ptr<TFile> out(TFile::Open(outputFilePath.c_str(), "RECREATE"));
    if (!out) {
        LOG(ERROR) << "Cannot open the output file\n";
        return false;
    }

    TTree *tree = new TTree("stats","Significance results");
    Float_t observed_pvalue(0);
    Float_t expected_pvalue(0);
    Float_t injected_pvalue(0);
    Float_t observed_significance(0);
    Float_t expected_significance(0);
    Float_t injected_significance(0);

    tree->Branch("obs_pvalue", &observed_pvalue);
    tree->Branch("exp_pvalue", &expected_pvalue);
    tree->Branch("inj_pvalue", &injected_pvalue);
    tree->Branch("obs_significance", &observed_significance);
    tree->Branch("exp_significance", &expected_significance);
    tree->Branch("inj_significance", &injected_significance);

    observed_pvalue = isBlind ? -1 : sigObserved.first;
    expected_pvalue = sigExpected.first;
    if (sigInjected.first > -90) {
        injected_pvalue = sigInjected.first;
    }
    observed_significance = isBlind ? -1 : observed;
    expected_significance = expected;
    if (sigInjected.first > -90) {
        injected_significance = injected;
    }
    tree->Fill();

    tree->SetDirectory(out.get());
    out->Write();
    out->Close();

    return isNaN;
}

//___________________________________________________________
//
void Common::SaveCanvasAs(const TCanvas& c, const std::string& path) {
    for(const auto& format: TRExFitter::IMAGEFORMAT){
        c.SaveAs((path+"."+format).c_str());
    }
}

//___________________________________________________________
//
std::shared_ptr<xRooNode> Common::XRooNodeSampleFromNode(std::shared_ptr<xRooNode> node, const std::string& region, const std::string& sample) {
    auto channel = node->find(region);
    if (!channel) {
        LOG(ERROR) << "Cannot read channel: " << region << "\n";
        exit(EXIT_FAILURE);
    }

    auto samples = channel->find("samples");
    if (!samples) {
        LOG(ERROR) << "Cannot find samples\n";
        exit(EXIT_FAILURE);
    }

    return samples->find(sample + "_" + region + "_shapes");
}

//___________________________________________________________
//
std::vector<std::string> Common::MatchingElememnts(const std::vector<std::string>& vector, const std::string& match) {
    std::vector<std::string> result;
    for (const auto& ielement : vector) {
        if (Common::StringsMatch(ielement, match)) {
            result.emplace_back(ielement);
        }
    }

    return result;
}

//___________________________________________________________
//
int Common::ColorFromRGB(const std::vector<std::string>& colours) {
    if (colours.size() != 3) {
        LOG(ERROR) << "Need to provide three colours\n";
        throw std::invalid_argument("");
    }

    auto convert_str_to_RGB = [] (const std::string& str) {
      int num = Common::convertStoNum<int>(str);
      if (num < 0) throw std::invalid_argument("RGB value out of range [0, 255]");
      if (num > 255) throw std::invalid_argument("RGB value out of range [0, 255]");
      return num;
    };

    // Convert a vector of strings to a vector of RGB values
    auto create_RGB_array = [&convert_str_to_RGB] (const std::vector<std::string>& str_vec) {
      std::array<int, 3> int_arr{{-1}};
      for (auto itr = str_vec.begin(); itr != str_vec.end(); ++itr) {
        int_arr.at(itr - str_vec.begin()) = convert_str_to_RGB(*itr);
      }
      return int_arr;
    };

    auto col_arr = create_RGB_array(colours);
    return TColor::GetColor(col_arr[0], col_arr[1], col_arr[2]);
}

//___________________________________________________________
//
std::unique_ptr<RooFitResult> Common::RooFitResultFromFile(TFile* file) {

    std::unique_ptr<RooFitResult> fr(nullptr);
    for(auto *key : *file->GetListOfKeys()){
        const std::string& keyName = key->GetName();
        const std::string nllName = "nll_simPdf";
        const std::string evaluatorName = "RooEvaluatorWrapper";
        const std::string evaluatorName2 = "nll_func_wrapper";
        auto itrNll = keyName.find(nllName);
        if (itrNll != std::string::npos) {
            fr.reset(file->Get<RooFitResult>(keyName.c_str()));
            break;
        } else if (keyName.find(evaluatorName) != std::string::npos) {
            fr.reset(file->Get<RooFitResult>(keyName.c_str()));
            break;
        } else if (keyName.find(evaluatorName2) != std::string::npos) {
            fr.reset(file->Get<RooFitResult>(keyName.c_str()));
            break;
        }
    }

    if (!fr) {
        LOG(ERROR) << "Cannot read the fit results\n";
        exit(EXIT_FAILURE);
    }

    return fr;
}

//___________________________________________________________
//
TString Common::xRooFitToysScanArg(int stepsSplusB, int stepsB) {
    // Format for number of toys argument is "nullToys.altToysFraction"
    // where the number of alt toys is "0.altToysFraction * nullToys".
    // An altToysFraction of 100% is represented by omitting the .altToysFraction

    if (stepsB < 0 || stepsSplusB < 0) {
        LOG(WARNING) << "Number of S+B and/or B-only toys < 0\n";
        return "";
    }

    if (stepsB > stepsSplusB) {
        LOG(WARNING) << "Number of B-only toys > S+B toys, not supported by xRooFit\n";
        return "";
    }

    TString sAltToysFraction;
    if  (stepsB < stepsSplusB) {
        sAltToysFraction = TString::Format("%f", 1.0 * stepsB / stepsSplusB);
        LOG(DEBUG) << "sAltToysFraction=" << sAltToysFraction << "\n";
        assert(sAltToysFraction.Index("0.") == 0);
        // 1 instead of 2, to reuse the period ('.')
        sAltToysFraction = sAltToysFraction(1, sAltToysFraction.Length());
    }

    return TString::Format("%d%s", stepsSplusB, sAltToysFraction.Data());
}

//___________________________________________________________
//
void Common::CleanTPad(TPad* pad) {
    TObject *obj;

    TIter next(pad->GetListOfPrimitives());
    while ((obj = next())) {
        if (auto hist = dynamic_cast<TH1*>(obj)) {
            hist->SetBit(kMustCleanup, false);
        }
        if (auto tl = dynamic_cast<TLatex*>(obj)) {
            tl->SetBit(kMustCleanup, false);
        }
        if (auto tli = dynamic_cast<TLine*>(obj)) {
            tli->SetBit(kMustCleanup, false);
        }
    }
}

//___________________________________________________________
//
bool Common::IsUnfoldingPOI(const std::string& poi, const std::vector<std::unique_ptr<Unfolding> >& unfoldings) {
    for (const auto& iunf : unfoldings) {
        if (poi == "TotalXsecOverTheory_"+iunf->fName) return true;
        for (int ibin = 0; ibin < iunf->fNumberUnfoldingTruthBins; ++ibin) {
            const std::string name = iunf->fName + "_Bin_" + Common::IntToFixLenStr(ibin+1) + "_mu";
            if (poi == name) return true;
        }
    }

    return false;
}