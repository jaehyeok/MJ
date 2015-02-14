#include "babytree.h"

using namespace std;

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
   
    //TString BabyDir = "HT500MET100/";  
    TString BabyDir = "HT750MET250/";  
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
    ch_f1500_100->Add(BabyDir+"baby_*f1500_100.root");
    ch_f1200_800->Add(BabyDir+"baby_*f1200_800.root");
    
    
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
    
    if(!OnlyDraw) 
    {
        // ----------------------------------------
        //  Fill histrograms 
        // ----------------------------------------
        MakeHists(ch_data	    ); 
        MakeHists(ch_ttbar_sl	); 
        MakeHists(ch_ttbar_ll	); 
        MakeHists(ch_wjets	    ); 
        MakeHists(ch_dy	        ); 
        MakeHists(ch_t	        ); 
        MakeHists(ch_f1500_100	);  
        MakeHists(ch_f1200_800	); 

        // ----------------------------------------
        //  Make the final histogram file
        // ----------------------------------------
        cout << "[MJ Analysis] Merging result files" << endl;
        gSystem->Exec("rm HistFiles/Hist.root");
        gSystem->Exec("hadd -f HistFiles/Hist.root HistFiles/*.root");
    }

    // ----------------------------------------
    //  Draw histograms 
    // ---------------------------------------- 
    
    Make1DPlots("muspT"    	    );
    Make1DPlots("musPhi"        );
    Make1DPlots("musEta"   	    );
    Make1DPlots("elspT"    	    );
    Make1DPlots("elsPhi"        );
    Make1DPlots("elsEta"   	    );
    Make1DPlots("mT"       	    );
    Make1DPlots("mj"       	    );
    Make1DPlots("MJ"       	    );
    Make1DPlots("HT"       	    );
    Make1DPlots("Nfatjet"       );
    Make1DPlots("Nskinnyjet"    );
    Make1DPlots("Ncsvm" 	    );
    Make1DPlots("MET"           );
    Make1DPlots("METPhi"        );
    Make1DPlots("WpT"           );
    Make1DPlots("FatjetPt1"     );
    Make1DPlots("FatjetPt2"     );
    Make1DPlots("FatjetPt3"     );
    Make1DPlots("FatjetPt4"     );
    Make1DPlots("FatjetPhi1"    );
    Make1DPlots("FatjetPhi2"    );
    Make1DPlots("FatjetPhi3"    );
    Make1DPlots("FatjetPhi4"    );
    Make1DPlots("FatjetEta1"    );
    Make1DPlots("FatjetEta2"    );
    Make1DPlots("FatjetEta3"    );
    Make1DPlots("FatjetEta4"    );
    Make1DPlots("mj1"   	  	);
    Make1DPlots("mj2"   	    );
    Make1DPlots("mj3"   	    );
    Make1DPlots("mj4"   	    );
    Make1DPlots("mj3overmj2"    );
    Make1DPlots("mj2overmj1"    );

// Yield table
    MakeTables(0, false);
    MakeTables(11, false);
    MakeTables(13, false);

}
