// The scripts behave exactly as "hadd -f", with the notable difference 
// that it will not add identical TH1, keeping just one copy of them
// So if calling
//     hupdate targetfile file1 file2
// with
//     file1 containing two TH1:   h1 h2
//     file2 containing two TH1:   h2 h3
// targetfile will contain:        h1 h2 h3, 
// with h2 taken only from file1 (not added to its copy in file2)
// The script is not (yet) iterating over subdirectories, not needed for its purposes
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

struct Folders {
   std::vector<std::string> regions;
   std::vector<std::string> samples;
   std::vector<std::string> systematics;
};

void MergeRootfile(std::unique_ptr<TFile>& target, const std::vector<std::unique_ptr<TFile> >& sources, const Folders& folders);

Folders CreateStructure(std::unique_ptr<TFile>& target, std::vector<std::unique_ptr<TFile> >& sources);

// main function

int main(int argc, char **argv){
   if (argc<3) {
      std::cout << "Too few arguments ..." <<std::endl    ;
      return -1;
   }

   std::cout << "Creating target mergin histograms from source files. Histograms not summed if have the same name: only the first one is kept." << std::endl;

   std::cout << "Target file:\t" << argv[1] << std::endl;

   std::unique_ptr<TFile> target(TFile::Open( argv[1], "RECREATE" ));
   std::vector<std::unique_ptr<TFile> > sources;

   for (int i=2; i<argc; i++){
      std::cout << "Source file("<< i-1<<"):\t" << argv[i] <<std::endl;
      std::unique_ptr<TFile> file(TFile::Open(argv[i], "READ"));
      if (!file) {
         std::cerr << "Cannot open file at " << argv[i] << "\n";
         std::exit(EXIT_FAILURE); 
      }
      sources.emplace_back(std::move(file));
   }

   Folders folders = CreateStructure(target, sources);

   // call the function
   MergeRootfile(target, sources, folders);

   target->Close();
   for (auto& ifile : sources) {
      ifile->Close();
   }

   std::cout << "Done." << std::endl;   

   return 0;
}

Folders CreateStructure(std::unique_ptr<TFile>& target, std::vector<std::unique_ptr<TFile> >& sources) {
   std::vector<std::string> regions;
   std::vector<std::string> samples;
   std::vector<std::string> systematics;

   auto FillVector =[](std::vector<std::string>& vec, const TDirectory* dir) {
      TIter next(dir->GetListOfKeys());
      TKey *key;
      while ((key = static_cast<TKey*>(next()))) {
         const std::string name = key->GetName();
         auto itr_reg = std::find(vec.begin(), vec.end(), name);
         if (itr_reg == vec.end()) {
            vec.emplace_back(name);
         }
      }
   };

   for (auto& ifile : sources) {
      ifile->cd();
      FillVector(regions, gDirectory);
   }
   for (auto& ifile : sources) {
      for (const auto& ireg : regions) {
         if (!ifile->Get<TDirectory>(ireg.c_str())) continue;
         ifile->cd(ireg.c_str());
         FillVector(samples, gDirectory);
      }
   }
   for (auto& ifile : sources) {
      for (const auto& ireg : regions) {
         for (const auto& isample : samples) {
            if (!ifile->Get<TDirectory>((ireg + "/" + isample).c_str())) continue;
            ifile->cd((ireg+"/"+isample).c_str());
            FillVector(systematics, gDirectory);
         }
      }
   }

   // make the folder
   for (const auto& ireg : regions) {
      target->cd();
      gDirectory->mkdir(ireg.c_str());
      for (const auto& isample : samples) {
         target->cd(ireg.c_str());
         gDirectory->mkdir(isample.c_str());
         for (const auto& isystematic : systematics) {
            target->cd((ireg+"/"+isample).c_str());
            gDirectory->mkdir(isystematic.c_str());
         }
      }
   }

   Folders result;
   result.regions = regions;
   result.samples = samples;
   result.systematics = systematics;

   return result;
}

void MergeRootfile(std::unique_ptr<TFile>& target, const std::vector<std::unique_ptr<TFile> >& sources, const Folders& folders) {
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   for (const auto& ifile : sources) {
      for (const auto& iregion : folders.regions) {
         for (const auto& isample : folders.samples) {
            for (const auto& isyst : folders.systematics) {
               TDirectory* dir = ifile->Get<TDirectory>((iregion + "/" + isample + "/" + isyst).c_str());
               if (!dir) continue;

               // iterate over all entries
               dir->cd();
               TIter next(dir->GetListOfKeys());
               TKey *key;
               while ((key = static_cast<TKey*>(next()))) {
                  TObject *obj = key->ReadObj();
                  if (!obj->IsA()->InheritsFrom(TH1::Class())) continue;

                  //check if it already exists in the target
                  const std::string name = key->GetName();
                  const TH1* h = target->Get<TH1>(((iregion + "/" + isample + "/" + isyst + "/" + name)).c_str());  
                  // is already there
                  if (h) continue;

                  // is not there, add it
                  target->cd((iregion + "/" + isample + "/" + isyst).c_str());
                  obj->Write(name.c_str());
               }
            }
         }
      }
   }

   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
} 
