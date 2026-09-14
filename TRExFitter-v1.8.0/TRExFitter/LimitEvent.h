#ifndef LimitEvent_h
#define LimitEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class LimitEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         obs_upperlimit;
   Float_t         inj_upperlimit;
   Float_t         exp_upperlimit;
   Float_t         exp_upperlimit_plus1;
   Float_t         exp_upperlimit_plus2;
   Float_t         exp_upperlimit_minus1;
   Float_t         exp_upperlimit_minus2;

   // List of branches
   TBranch        *b_obs_upperlimit;   //!
   TBranch        *b_inj_upperlimit;   //!
   TBranch        *b_exp_upperlimit;   //!
   TBranch        *b_exp_upperlimit_plus1;   //!
   TBranch        *b_exp_upperlimit_plus2;   //!
   TBranch        *b_exp_upperlimit_minus1;   //!
   TBranch        *b_exp_upperlimit_minus2;   //!
   
   explicit LimitEvent(TTree *tree=0);
   virtual ~LimitEvent();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
};

#endif
