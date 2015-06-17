
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
    gROOT->LoadMacro("MakeCards.C+");

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
  
    //TString BabyDir = "/Users/jaehyeok/Research/Tools/fastjet-3.0.6/example/babies/13TeV/Phys14/HT750MET250/";
    //TString BabyDir = "/Users/jaehyeok/Research/Tools/fastjet-3.0.6/example/babies/13TeV/Phys14_JetPt20_16Mar2015/"/*_HT750MET250/"*/;
    TString BabyDir = "/Users/jaehyeok/scratch/2015_05_25/skim/"/*_HT750MET250/"*/;
    
    // Data
    //ch_data->Add(BabyDir+"baby_MuHad_*.root");                            
    
    // TT 
    ch_ttbar_sl->Add(BabyDir+"*TTJets*.root");
    ch_ttbar_ll->Add(BabyDir+"*TTJets*.root");
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
    cout << "TTV                : " << ch_ttv->GetEntries()       << endl;
    cout << "T1tttt(1500,100)   : " << ch_f1500_100->GetEntries()    << endl;
    cout << "T1tttt(1200,8000)  : " << ch_f1200_800->GetEntries()  << endl;
   

    //
    // Loop over SR and CR : make sure that these regions exist in "PassSelection.h"
    //
    //char* Region[] = {"Baseline","SR0", "SR1", "SR2", "SR3", "SR4", "SR5", "SR6", "SR7", "SR8", "SR9"}; 
    char* Region[] = {"TEST"}; 
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
            MakeHists(ch_data,	    Region[iregion]); 
            MakeHists(ch_ttbar_sl,  Region[iregion]); 
            MakeHists(ch_ttbar_ll,  Region[iregion]); 
            MakeHists(ch_wjets,     Region[iregion]); 
            MakeHists(ch_dy,	    Region[iregion]); 
            MakeHists(ch_t,         Region[iregion]); 
            MakeHists(ch_ttv,       Region[iregion]); 
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
/*
        Make1DPlots("dRlep",        Region[iregion]);
        Make1DPlots("dPhiMET",      Region[iregion]);
        Make1DPlots("dRbmin",       Region[iregion]);
        Make1DPlots("dPhiMETlep",   Region[iregion]);
        Make1DPlots("muspT",        Region[iregion]);
        Make1DPlots("musPhi",       Region[iregion]);
        Make1DPlots("musEta",       Region[iregion]);
        Make1DPlots("elspT",        Region[iregion]);
        Make1DPlots("elsPhi",       Region[iregion]);
        Make1DPlots("elsEta",       Region[iregion]);
        Make1DPlots("mT",           Region[iregion]);
        Make1DPlots("mj",           Region[iregion]);
        Make1DPlots("MJ",           Region[iregion]);
        Make1DPlots("MJ_ISR",       Region[iregion]);
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
*/
        Make1DPlots("mj08_1",        Region[iregion]);
        Make1DPlots("mj08_2",        Region[iregion]);
        Make1DPlots("mj08_3",        Region[iregion]);
        Make1DPlots("mj08_4",        Region[iregion]);
        Make1DPlots("mindPhibb",     Region[iregion]);
        // ----------------------------------------
        //  Make table of yields 
        // ---------------------------------------- 
        MakeTables(0,   Region[iregion], false);
        MakeTables(11,  Region[iregion], false);
        MakeTables(13,  Region[iregion], false);
        
        // ----------------------------------------
        //  Make cards for combine/LandS 
        // ---------------------------------------- 
        //MakeCards(0,   Region[iregion]);
        //MakeCards(11,  Region[iregion]);
        //MakeCards(13,  Region[iregion]);
    } //for(int iregion=0; iregion<2; iregion++)

}
