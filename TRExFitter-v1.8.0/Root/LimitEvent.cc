#include "TRExFitter/LimitEvent.h"

LimitEvent::LimitEvent(TTree *tree) : fChain(0) 
{
   Init(tree);
}

LimitEvent::~LimitEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}
Int_t LimitEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t LimitEvent::LoadTree(Long64_t entry) {
// Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
       fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}

void LimitEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("obs_upperlimit", &obs_upperlimit, &b_obs_upperlimit);
   fChain->SetBranchAddress("inj_upperlimit", &inj_upperlimit, &b_inj_upperlimit);
   fChain->SetBranchAddress("exp_upperlimit", &exp_upperlimit, &b_exp_upperlimit);
   fChain->SetBranchAddress("exp_upperlimit_plus1", &exp_upperlimit_plus1, &b_exp_upperlimit_plus1);
   fChain->SetBranchAddress("exp_upperlimit_plus2", &exp_upperlimit_plus2, &b_exp_upperlimit_plus2);
   fChain->SetBranchAddress("exp_upperlimit_minus1", &exp_upperlimit_minus1, &b_exp_upperlimit_minus1);
   fChain->SetBranchAddress("exp_upperlimit_minus2", &exp_upperlimit_minus2, &b_exp_upperlimit_minus2);
}
