
//
// main 
//
void DoAnalysis(bool OnlyDraw=false) 
{
    // Style
    gROOT->ProcessLine(".L /Users/jaehyeok/macros/rootlogon.C");
    // Load macros 
    gROOT->LoadMacro("MakeHists.C+");
    gROOT->LoadMacro("Make1DPlots.C+");
    gROOT->LoadMacro("Make2DPlots.C+");
    gROOT->LoadMacro("MakeTables.C+");

    // ----------------------------------------
    //  Define chains  
    // ----------------------------------------
    TChain *ch_data         = new TChain("tree", "DATA");
    TChain *ch_ttbar_sl     = new TChain("tree", "TT_sl");
    TChain *ch_ttbar_ll     = new TChain("tree", "TT_ll");
    TChain *ch_wjets        = new TChain("tree", "WJets");
    TChain *ch_dy           = new TChain("tree", "DY");
    TChain *ch_t            = new TChain("tree", "T");
    TChain *ch_f1500_100    = new TChain("tree", "T1tttt_f1500_100");
    TChain *ch_f1200_800    = new TChain("tree", "T1tttt_f1200_800");
  
    //TString BabyDir = "/Users/jaehyeok/Research/Tools/fastjet-3.0.6/example/babies/13TeV/HT750MET250/";
    TString BabyDir = "/Users/jaehyeok/Research/Tools/fastjet-3.0.6/example/babies/13TeV/Phys14/HT750MET250/";
    
    // Data
    //ch_data->Add(BabyDir+"baby_MuHad_*.root");                            
    
    // TT 
    ch_ttbar_sl->Add(BabyDir+"baby_TTJets*.root");
    ch_ttbar_ll->Add(BabyDir+"baby_TTJets*.root");
    // WJets 
    //ch_wjets->Add(BabyDir+"baby_WJetsToLNu*f1To*.root");
    ch_wjets->Add(BabyDir+"baby_WJetsToLNu*.root");
    // DY 
    //ch_dy->Add(BabyDir+"baby_DYJetsToLL*f1To*.root");
    ch_dy->Add(BabyDir+"baby_DYJetsToLL*.root");
    // Singla top 
    ch_t->Add(BabyDir+"baby_*channel*_f*.root");

    // Signal
    ch_f1500_100->Add(BabyDir+"baby_*1500*.root");
    ch_f1200_800->Add(BabyDir+"baby_*1200*.root");
    
    
    // ----------------------------------------
    //  Get number of entries 
    // ----------------------------------------
    cout << "data               : " << ch_data->GetEntries()        << endl;
    cout << "ttbarl             : " << ch_ttbar_sl->GetEntries()    << endl;
    cout << "ttbarll            : " << ch_ttbar_ll->GetEntries()    << endl;
    cout << "wjets              : " << ch_wjets->GetEntries()       << endl;
    cout << "dy                 : " << ch_dy->GetEntries()          << endl;
    cout << "Single top         : " << ch_t->GetEntries()           << endl;
    cout << "T1tttt(1500,100)   : " << ch_f1500_100->GetEntries()    << endl;
    cout << "T1tttt(1200,8000)  : " << ch_f1200_800->GetEntries()  << endl;
   

    //
    // Loop over SR and CR : make sure that these regions exist in "PassSelection.h"
    //
    //char* Region[11] = {"Baseline","SR0", "SR1", "SR2", "SR3", "SR4", "SR5", "SR6", "SR7", "SR8", "SR9"}; 
    char* Region[1] = {"SR20p1"}; 

    for(int iregion=0; iregion<1; iregion++)
    {
        cout << endl;
        cout << "[MJ Analysis] Analyzing " << Region[iregion] << endl;
        cout << endl;
        
        if(!OnlyDraw) 
        {
            // ----------------------------------------
            //  Fill histrograms 
            // ----------------------------------------
            MakeHists(ch_data,	    Region[iregion]); 
            MakeHists(ch_ttbar_sl,  Region[iregion]); 
            MakeHists(ch_ttbar_ll,  Region[iregion]); 
            MakeHists(ch_wjets,     Region[iregion]); 
            MakeHists(ch_dy,	    Region[iregion]); 
            MakeHists(ch_t,         Region[iregion]); 
            MakeHists(ch_f1500_100, Region[iregion]);  
            MakeHists(ch_f1200_800, Region[iregion]); 

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
        Make1DPlots("muspT",        Region[iregion]);
        Make1DPlots("musPhi",       Region[iregion]);
        Make1DPlots("musEta",       Region[iregion]);
        Make1DPlots("elspT",        Region[iregion]);
        Make1DPlots("elsPhi",       Region[iregion]);
        Make1DPlots("elsEta",       Region[iregion]);
        Make1DPlots("mT",           Region[iregion]);
        Make1DPlots("mj",           Region[iregion]);
        Make1DPlots("MJ",           Region[iregion]);
        Make1DPlots("HT",           Region[iregion]);
        Make1DPlots("Nfatjet",      Region[iregion]);
        Make1DPlots("Nskinnyjet",   Region[iregion]);
        Make1DPlots("Ncsvm",        Region[iregion]);
        Make1DPlots("MET",          Region[iregion]);
        Make1DPlots("METPhi",       Region[iregion]);
        Make1DPlots("WpT",          Region[iregion]);
        Make1DPlots("FatjetPt1",    Region[iregion]);
        Make1DPlots("FatjetPt2",    Region[iregion]);
        Make1DPlots("FatjetPt3",    Region[iregion]);
        Make1DPlots("FatjetPt4",    Region[iregion]);
        Make1DPlots("FatjetPhi1",   Region[iregion]);
        Make1DPlots("FatjetPhi2",   Region[iregion]);
        Make1DPlots("FatjetPhi3",   Region[iregion]);
        Make1DPlots("FatjetPhi4",   Region[iregion]);
        Make1DPlots("FatjetEta1",   Region[iregion]);
        Make1DPlots("FatjetEta2",   Region[iregion]);
        Make1DPlots("FatjetEta3",   Region[iregion]);
        Make1DPlots("FatjetEta4",   Region[iregion]);
        Make1DPlots("mj1",          Region[iregion]);
        Make1DPlots("mj2",          Region[iregion]);
        Make1DPlots("mj3",          Region[iregion]);
        Make1DPlots("mj4",          Region[iregion]);
        Make1DPlots("mj1OverMJ",    Region[iregion]);
        Make1DPlots("mj2OverMJ",    Region[iregion]);
        Make1DPlots("mj3OverMJ",    Region[iregion]);
        Make1DPlots("mj4OverMJ",    Region[iregion]);
        Make1DPlots("N1",           Region[iregion]);
        Make1DPlots("N2",           Region[iregion]);
        Make1DPlots("N3",           Region[iregion]);
        Make1DPlots("N4",           Region[iregion]);
        Make1DPlots("mjOverPt1",    Region[iregion]);
        Make1DPlots("mjOverPt2",    Region[iregion]);
        Make1DPlots("mjOverPt3",    Region[iregion]);
        Make1DPlots("mjOverPt4",    Region[iregion]);
        Make1DPlots("mj3overmj2",   Region[iregion]);
        Make1DPlots("mj2overmj1",   Region[iregion]);

        // ----------------------------------------
        //  Make table of yields 
        // ---------------------------------------- 
        MakeTables(0,   Region[iregion], false);
        MakeTables(11,  Region[iregion], false);
        MakeTables(13,  Region[iregion], false);
        
    } //for(int iregion=0; iregion<2; iregion++)

}
