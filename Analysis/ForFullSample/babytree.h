//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 23 12:09:46 2014 by ROOT version 5.34/18
// from TTree tree/A Baby Ntuple
// found on file: ../../baby/baby_DoubleMu_Run2012D-PromptReco2.root
//////////////////////////////////////////////////////////

#include <TTree.h>
#include <TChain.h>
#include <vector>
#include <string>

   // Declaration of leaf types
   int           run;
   int           lumiblock;
   int           event;
//   int           Nfatjet;
   int           Nskinnyjet;
   int           Npv;
   float         Npu;
   float         EventWeight;
   float         mll;
//   float         MJ;
   float         HT;
//   vector<float>   *mj;
//   vector<float>   *FatjetPt;
//   vector<float>   *FatjetEta;
//   vector<float>   *FatjetPhi;
   vector<float>   *RA4leptonPt;
   vector<float>   *RA4leptonEta;
   vector<float>   *RA4leptonPhi;
   vector<float>   *RA4leptonId;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_event;   //!
//   TBranch        *b_Nfatjet;   //!
   TBranch        *b_Nskinnyjet;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_Npu;   //!
   TBranch        *b_EventWeight;   //!
   TBranch        *b_mll;   //!
//   TBranch        *b_MJ;   //!
   TBranch        *b_HT;   //!
//   TBranch        *b_mj;   //!
//   TBranch        *b_FatjetPt;   //!
//   TBranch        *b_FatjetEta;   //!
//   TBranch        *b_FatjetPhi;   //!
   TBranch        *b_RA4leptonPt;   //!
   TBranch        *b_RA4leptonEta;   //!
   TBranch        *b_RA4leptonPhi;   //!
   TBranch        *b_RA4leptonId;   //!

void InitTree(TChain *fChain)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
//   mj = 0;
//   FatjetPt = 0;
//   FatjetEta = 0;
//   FatjetPhi = 0;
   RA4leptonPt = 0;
   RA4leptonEta = 0;
   RA4leptonPhi = 0;
   RA4leptonId = 0;
   // Set branch addresses and branch pointers
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("event", &event, &b_event);
//   fChain->SetBranchAddress("Nfatjet", &Nfatjet, &b_Nfatjet);
   fChain->SetBranchAddress("Nskinnyjet", &Nskinnyjet, &b_Nskinnyjet);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("Npu", &Npu, &b_Npu);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
   fChain->SetBranchAddress("mll", &mll, &b_mll);
//   fChain->SetBranchAddress("MJ", &MJ, &b_MJ);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
//   fChain->SetBranchAddress("mj", &mj, &b_mj);
//   fChain->SetBranchAddress("FatjetPt", &FatjetPt, &b_FatjetPt);
//   fChain->SetBranchAddress("FatjetEta", &FatjetEta, &b_FatjetEta);
//   fChain->SetBranchAddress("FatjetPhi", &FatjetPhi, &b_FatjetPhi);
   fChain->SetBranchAddress("RA4leptonPt", &RA4leptonPt, &b_RA4leptonPt);
   fChain->SetBranchAddress("RA4leptonEta", &RA4leptonEta, &b_RA4leptonEta);
   fChain->SetBranchAddress("RA4leptonPhi", &RA4leptonPhi, &b_RA4leptonPhi);
   fChain->SetBranchAddress("RA4leptonId", &RA4leptonId, &b_RA4leptonId);
}
