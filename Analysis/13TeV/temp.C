#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>

#include "MakeHists.C" 

using namespace std;
void temp()
{
    // Load macros 
    //int loadMakeHists   = gROOT->LoadMacro("MakeHists.C+");

    //cout << "Loading MakeHists.C    : " << loadMakeHists << endl;

    float Lumi = 40.; // pb-1

    // ----------------------------------------
    //  Define chains  
    // ----------------------------------------
    TChain *ch_ttbar_sl     = new TChain("tree", "TT_sl");
  
    //TString BabyDir = "/Volumes/Data/Rtuples/2015_07_22/skim_1lht400/";
    TString BabyDir = "/Users/jaehyeok/scratch/2015_05_25/skim/";
    
    // TT 
    ch_ttbar_sl->Add(BabyDir+"*TTJets*.root");
    

    cout << "ttbarl             : " << ch_ttbar_sl->GetEntries()    << endl;

    //
    // Loop over SR and CR : make sure that these regions exist in "PassSelection.h"
    //
    MakeHists(ch_ttbar_sl,  (char*)"Method1",	Lumi); 
}
