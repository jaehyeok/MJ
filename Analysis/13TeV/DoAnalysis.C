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
    
    float Lumi = 2100; // pb-1

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
    TChain *ch_qcd          = new TChain("tree", "QCD");
    TChain *ch_others       = new TChain("tree", "Others");
    
    TChain *ch_f1500_100    = new TChain("tree", "T1tttt_f1500_100");
    TChain *ch_f1200_800    = new TChain("tree", "T1tttt_f1200_800");
  

    TString BabyDirData = "/Users/jaehyeok/Research/cms/UCSB/MJ/Analysis/susy_cfa_babies/2015_11_20/data/singlelep/combined/skim_abcd/";
    TString BabyDirMC = "/Users/jaehyeok/Research/cms/UCSB/MJ/Analysis/susy_cfa_babies/2015_11_28/mc/skim_abcd/";
    
    // Data
    ch_data->Add(BabyDirData+"*Single*.root");                            
    
    // TT 
    ch_ttbar_sl->Add(BabyDirMC+"*TTJets*Lept*");
    ch_ttbar_sl->Add(BabyDirMC+"*TTJets_HT*");
    ch_ttbar_ll->Add(BabyDirMC+"*TTJets*Lept*");
    ch_ttbar_ll->Add(BabyDirMC+"*TTJets_HT*");
    // WJets 
    ch_wjets->Add(BabyDirMC+"*_WJetsToLNu*");
    // Singla top 
    ch_t->Add(BabyDirMC+"*_ST_*");
    // TTV 
    ch_ttv->Add(BabyDirMC+"*_TTWJets*");
    ch_ttv->Add(BabyDirMC+"*_TTZTo*");
    ch_ttv->Add(BabyDirMC+"*_TTG*");
    // QCD   
    ch_qcd->Add(BabyDirMC+"*_QCD_HT*");
    //ch_qcd->Add(BabyDirMC+"*_TTJets_TuneCUET*");
    // DY 
    ch_dy->Add(BabyDirMC+"*DYJetsToLL*");
    // Others 
    ch_others->Add(BabyDirMC+"*_ZJet*");
    ch_others->Add(BabyDirMC+"*_WWTo*");
    ch_others->Add(BabyDirMC+"*_WZTo*");
    ch_others->Add(BabyDirMC+"*ggZH_HToBB*");
    ch_others->Add(BabyDirMC+"*ttHJetTobb*");
    ch_others->Add(BabyDirMC+"*_TTTT*");

    // Signal
    ch_f1500_100->Add(BabyDirMC+"*1500*100*");
    ch_f1200_800->Add(BabyDirMC+"*1200*800*");
    
    if(!OnlyDraw) 
    {
    // ----------------------------------------
    //  Get number of entries 
    // ----------------------------------------
    cout << "data               : " << ch_data->GetEntries()        << endl;
    cout << "ttbar(1l)          : " << ch_ttbar_sl->GetEntries()    << endl;
    cout << "ttbar(2l)          : " << ch_ttbar_ll->GetEntries()    << endl;
    cout << "wjets              : " << ch_wjets->GetEntries()       << endl;
    cout << "dy                 : " << ch_dy->GetEntries()          << endl;
    cout << "single top         : " << ch_t->GetEntries()           << endl;
    cout << "ttb                : " << ch_ttv->GetEntries()         << endl;
    cout << "others             : " << ch_others->GetEntries()      << endl;
    cout << "T1tttt(1500,100)   : " << ch_f1500_100->GetEntries()   << endl;
    cout << "T1tttt(1200,8000)  : " << ch_f1200_800->GetEntries()   << endl;
    }

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
            MakeHists(ch_qcd,       Region[iregion],	Lumi); 
            MakeHists(ch_others,    Region[iregion],	Lumi); 
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
        Make1DPlots("lepspT",       Region[iregion],	1,	false,	Lumi);
        Make1DPlots("lepsPhi",      Region[iregion],	1,	false,	Lumi);
        Make1DPlots("lepsEta",      Region[iregion],	1,	false,	Lumi);
        Make1DPlots("mT",           Region[iregion],	1,	false,	Lumi);
        Make1DPlots("mj",           Region[iregion],	1,	false,	Lumi);
        Make1DPlots("MJ",           Region[iregion],	2,	false,	Lumi);
        Make1DPlots("HT",           Region[iregion],	2,	false,	Lumi);
        Make1DPlots("Nfatjet",      Region[iregion],	1,	false,	Lumi);
        Make1DPlots("Nskinnyjet",   Region[iregion],	1,	false,	Lumi);
        Make1DPlots("Ncsvm",        Region[iregion],	1,	false,	Lumi);
        Make1DPlots("MET",          Region[iregion],	1,	false,	Lumi);
        Make1DPlots("mj1",          Region[iregion],	1,	false,	Lumi);
        Make1DPlots("mj2",          Region[iregion],	1,	false,	Lumi);
        Make1DPlots("mj3",          Region[iregion],	1,	false,	Lumi);
        Make1DPlots("mj4",          Region[iregion],	1,	false,	Lumi);
        Make1DPlots("mbb",          Region[iregion],	1,	true,	Lumi);
        
        // Corroborators
        //Make1DPlots("mj08_1",        Region[iregion],	1,	false,	Lumi);
        //Make1DPlots("mj08_2",        Region[iregion],	1,	false,	Lumi);
        //Make1DPlots("mj08_3",        Region[iregion],	1,	false,	Lumi);
        //Make1DPlots("mj08_4",        Region[iregion],	1,	false,	Lumi);
        //Make1DPlots("mindPhibb",     Region[iregion],	1,	false,	Lumi);
        
        // ----------------------------------------
        //  Make table of yields 
        // ---------------------------------------- 
        MakeTables(0,   Region[iregion], false,	Lumi);
        //MakeTables(11,  Region[iregion], false,	Lumi);
        //MakeTables(13,  Region[iregion], false,	Lumi);
        
        MakeTablesAllRegions(0,   Region[iregion], false,	Lumi);
        
        // ----------------------------------------
        //  Make cards for combine/LandS 
        // ---------------------------------------- 
        MakeCards(0,   Region[iregion],	Lumi);
        //MakeCards(11,  Region[iregion],	Lumi);
        //MakeCards(13,  Region[iregion],	Lumi);

    } //for(int iregion=0; iregion<NRegion; iregion++)

}
