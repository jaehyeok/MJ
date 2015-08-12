#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>

#include "MakeHists.C" 
#include "Make1DPlots.C"
#include "Make2DPlots.C"
#include "MakeTables.C"
#include "MakeCards.C"

using namespace std;

//
// main 
//
void DoAnalysis(bool OnlyDraw=false) 
{
    // Style
    gROOT->ProcessLine(".L /Users/jaehyeok/macros/rootlogon.C");
/*    
    // Load macros 
    int loadMakeHists   = gROOT->LoadMacro("MakeHists.C+");
    int loadMake1DPlots = gROOT->LoadMacro("Make1DPlots.C+");
    int loadMake2DPlots =  gROOT->LoadMacro("Make2DPlots.C+");
    int loadMakeTables  =  gROOT->LoadMacro("MakeTables.C+");
    int loadMakeCards   =  gROOT->LoadMacro("MakeCards.C+");

    cout << "Loading MakeHists.C    : " << (loadMakeHists==0?"Loaded":"Not loaded")   << endl;
    cout << "Loading Make1DPlots.C  : " << (loadMake1DPlots==0?"Loaded":"Not loaded") << endl;
    cout << "Loading Make2DPlots.C  : " << (loadMake2DPlots==0?"Loaded":"Not loaded") << endl;
    cout << "Loading MakeTables.C   : " << (loadMakeTables==0?"Loaded":"Not loaded")  << endl;
    cout << "Loading MakeCards.C    : " << (loadMakeCards==0?"Loaded":"Not loaded")   << endl;
    cout << endl; 
*/
    float Lumi = 41.6; // pb-1

    // ----------------------------------------
    //  Define chains  
    // ----------------------------------------
    TChain *ch_data         = new TChain("tree", "DATA");
    TChain *ch_ttbar_sl     = new TChain("tree", "TT_sl");
    TChain *ch_ttbar_ll     = new TChain("tree", "TT_ll");
    TChain *ch_wjets        = new TChain("tree", "WJets");
    TChain *ch_dy           = new TChain("tree", "DY");
    TChain *ch_t            = new TChain("tree", "T");
    TChain *ch_ttv          = new TChain("tree", "TTV");
    TChain *ch_f1500_100    = new TChain("tree", "T1tttt_f1500_100");
    TChain *ch_f1200_800    = new TChain("tree", "T1tttt_f1200_800");
  

    TString BabyDir = "/Users/jaehyeok/scratch/2015_07_26/skim_1lht400/";
    
    // Data
    //ch_data->Add(BabyDir+"*JetHT*.root");                            
    ch_data->Add(BabyDir+"*HTMHT*.root");                            
    
    // TT 
    ch_ttbar_sl->Add(BabyDir+"*TTJets*25ns*.root");
    ch_ttbar_ll->Add(BabyDir+"*TTJets*25ns*.root");
    // WJets 
    ch_wjets->Add(BabyDir+"*WJetsToLNu*.root");
    // DY 
    ch_dy->Add(BabyDir+"*DYJetsToLL*.root");
    // Singla top 
    ch_t->Add(BabyDir+"*T*channel*.root");
    // TTV 
    ch_ttv->Add(BabyDir+"*TTW*.root");
    ch_ttv->Add(BabyDir+"*TTZ*.root");
    ch_ttv->Add(BabyDir+"*WH_HToBB*.root");

    // Signal
    ch_f1500_100->Add(BabyDir+"*_*mGl-1500_mLSP-100*.root");
    ch_f1200_800->Add(BabyDir+"*_*mGl-1200_mLSP-800*.root");
    
    // ----------------------------------------
    //  Get number of entries 
    // ----------------------------------------
    cout << "data               : " << ch_data->GetEntries()        << endl;
    cout << "ttbarl             : " << ch_ttbar_sl->GetEntries()    << endl;
    cout << "ttbarll            : " << ch_ttbar_ll->GetEntries()    << endl;
    cout << "wjets              : " << ch_wjets->GetEntries()       << endl;
    cout << "dy                 : " << ch_dy->GetEntries()          << endl;
    cout << "Single top         : " << ch_t->GetEntries()           << endl;
    cout << "TTV                : " << ch_ttv->GetEntries()         << endl;
    cout << "T1tttt(1500,100)   : " << ch_f1500_100->GetEntries()   << endl;
    cout << "T1tttt(1200,8000)  : " << ch_f1200_800->GetEntries()   << endl;

    //
    // Loop over SR and CR : make sure that these regions exist in "PassSelection.h"
    //
    char* Region[] = {(char*)"TEST"}; 
    int NRegion = sizeof(Region)/sizeof(Region[0]);

    for(int iregion=0; iregion<NRegion; iregion++)
    {
        cout << endl;
        cout << "[MJ Analysis] Analyzing " << Region[iregion] << endl;
        cout << endl;
        cout << "[MJ Analysis] Making directory for figures : Figures/" << Region[iregion] << endl;
        gSystem->mkdir(Form("Figures/%s",Region[iregion])); 
       
        if(!OnlyDraw) 
        {
            // ----------------------------------------
            //  Fill histrograms 
            // ----------------------------------------
            MakeHists(ch_data,	    Region[iregion],	Lumi); 
            MakeHists(ch_ttbar_sl,  Region[iregion],	Lumi); 
            MakeHists(ch_ttbar_ll,  Region[iregion],	Lumi); 
            MakeHists(ch_wjets,     Region[iregion],	Lumi); 
            MakeHists(ch_dy,	    Region[iregion],	Lumi); 
            MakeHists(ch_t,         Region[iregion],	Lumi); 
            MakeHists(ch_ttv,       Region[iregion],	Lumi); 
            MakeHists(ch_f1500_100, Region[iregion],	Lumi);  
            MakeHists(ch_f1200_800, Region[iregion],	Lumi); 
            
            // ----------------------------------------
            //  Make the final histogram file
            // ----------------------------------------
            cout << "[MJ Analysis] Merging result files" << endl;
            gSystem->Exec(Form("rm HistFiles/Hist_%s.root", Region[iregion]));
            gSystem->Exec(Form("hadd -f HistFiles/Hist_%s.root HistFiles/*_%s.root", Region[iregion], Region[iregion]));
            gSystem->Exec(Form("mv HistFiles/Hist_%s.root HistFiles/Hist_%s.root.tmp", Region[iregion], Region[iregion]));
            gSystem->Exec(Form("rm HistFiles/*_%s.root", Region[iregion]));
            gSystem->Exec(Form("mv HistFiles/Hist_%s.root.tmp HistFiles/Hist_%s.root", Region[iregion], Region[iregion]));

        }


        // ----------------------------------------
        //  Draw histograms 
        // ---------------------------------------- 

        // basic kinematic variables
        Make1DPlots("muspT",        Region[iregion],	1,	Lumi);
        Make1DPlots("musPhi",       Region[iregion],	1,	Lumi);
        Make1DPlots("musEta",       Region[iregion],	1,	Lumi);
        Make1DPlots("elspT",        Region[iregion],	1,	Lumi);
        Make1DPlots("elsPhi",       Region[iregion],	1,	Lumi);
        Make1DPlots("elsEta",       Region[iregion],	1,	Lumi);
        Make1DPlots("mT",           Region[iregion],	1,	Lumi);
        Make1DPlots("mj",           Region[iregion],	1,	Lumi);
        Make1DPlots("MJ",           Region[iregion],	2,	Lumi);
        Make1DPlots("HT",           Region[iregion],	2,	Lumi);
        Make1DPlots("Nfatjet",      Region[iregion],	1,	Lumi);
        Make1DPlots("Nskinnyjet",   Region[iregion],	1,	Lumi);
        Make1DPlots("Ncsvm",        Region[iregion],	1,	Lumi);
        Make1DPlots("MET",          Region[iregion],	1,	Lumi);
        Make1DPlots("mj1",          Region[iregion],	1,	Lumi);
        Make1DPlots("mj2",          Region[iregion],	1,	Lumi);
        Make1DPlots("mj3",          Region[iregion],	1,	Lumi);
        Make1DPlots("mj4",          Region[iregion],	1,	Lumi);
        
        // Corroborators
        //Make1DPlots("mj08_1",        Region[iregion],	1,	Lumi);
        //Make1DPlots("mj08_2",        Region[iregion],	1,	Lumi);
        //Make1DPlots("mj08_3",        Region[iregion],	1,	Lumi);
        //Make1DPlots("mj08_4",        Region[iregion],	1,	Lumi);
        //Make1DPlots("mindPhibb",     Region[iregion],	1,	Lumi);
        
        // ----------------------------------------
        //  Make table of yields 
        // ---------------------------------------- 
        MakeTables(0,   Region[iregion], false,	Lumi);
        //MakeTables(11,  Region[iregion], false,	Lumi);
        //MakeTables(13,  Region[iregion], false,	Lumi);
        
        // ----------------------------------------
        //  Make cards for combine/LandS 
        // ---------------------------------------- 
        MakeCards(0,   Region[iregion],	Lumi);
        //MakeCards(11,  Region[iregion],	Lumi);
        //MakeCards(13,  Region[iregion],	Lumi);

    } //for(int iregion=0; iregion<NRegion; iregion++)

}
