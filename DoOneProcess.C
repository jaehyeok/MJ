#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"

#include "Branch.h"
#include "ObjectSelector_Sync.h"

using namespace std;

//
// Get mj 
//
float Getmj(double px, double py, double pz, double E){

    float mj = TMath::Sqrt(E*E - px*px - py*py - pz*pz); 
    return mj;
}

//
// Get MJ 
//
float GetMJ(vector<float> Vectormj){
    
    if(Vectormj.size()==0) {
   
        return -999.;
    
    } else {

        float MJ = 0.;
        for(int imj=0; imj<(int)Vectormj.size(); imj++) MJ = MJ + Vectormj.at(imj);
        return MJ;
    
    }
}

//
// main function
//
void DoOneProcess(TString InputCfAName, TString ProcessName, float weight) { 
    
    // 
    // Get tree 
    // 
    TChain * chainA = new TChain("/configurableAnalysis/eventA");   
    TChain * chainB = new TChain("/configurableAnalysis/eventB");   
    chainA->Add(InputCfAName);
    chainB->Add(InputCfAName);
  
    InitializeA(chainA);
    InitializeB(chainB);

    // 
    // histograms
    // 
    TH1F *h1_result       = new TH1F("h1_result","h1_result", 1, 0, 1);
    
    //
    // main event loop
    //
    Int_t nentries = (Int_t)chainA->GetEntries();
    cout<<"The number of entries is: "<<nentries<<endl;
    // Progress tracking 
    int i_permille_old = 0; 
    TDatime DTStart;
    int StartDate = DTStart.GetDate(); 
    int StartTime = DTStart.GetTime(); 
    cout << "Start time : " << StartTime << endl; 
    
    for(int i = 0; i<nentries; i++) {

        // Progress tracking begin -------------------------------- 
        int i_permille = (int)floor(1000 * i / float(nentries));
        TDatime DTCurrent;
        int CurrentDate = DTCurrent.GetDate();
        int CurrentTime = DTCurrent.GetTime();
        int TimeLaps = (CurrentDate-StartDate)*1000000+(CurrentTime-StartTime);
        int TimeToRun = (float)nentries/(float)i*TimeLaps;
        if (i_permille != i_permille_old) {
            // xterm magic from L. Vacavant and A. Cerri
            if (isatty(1)) {
                printf("\015\033[32m Processed :: \033[1m\033[31m%4.1f %%" 
                       "\033[0m\033[32m   Expected processing time :: \033[1m\033[31m%i:%i:%i \033[0m\015",
                        i_permille/10., (TimeToRun/10000)%100<60 ? (TimeToRun/10000)%100 : (TimeToRun/10000)%100-40, 
                                        (TimeToRun/100)%100<60 ? (TimeToRun/100)%100 : (TimeToRun/100)%100-40, 
                                        (TimeToRun%100)<60 ? (TimeToRun)%100 : (TimeToRun)%100-40 );
                fflush(stdout);
            }
            i_permille_old = i_permille;
        } 
        // Progress tracking end ----------------------------------
        
        // get an entry of an event 
        chainA->GetEntry(i);
        chainB->GetEntry(i);
        
        //
        // Core analysis 
        //

        // Getting good muons
        vector<int> RA4MuonVeto; RA4MuonVeto.clear();
        vector<int> RA4Muon = GetRA4Muon(RA4MuonVeto);
        // Getting good electrons
        vector<int> RA4ElecVeto; RA4ElecVeto.clear();
        vector<int> RA4Elec = GetRA4Elec(RA4ElecVeto, "", 0, true);
        // Containers for B-tagged jets 
        vector<int> LooseBJet; 
        vector<int> MediumBJet; 
        // HT
        double HT=-999.; 
        // MJ 
        vector<float> Vector_mj;   // mj
        double MJ=-999.; 

        // Getting good skinny jets 
        vector<int> GoodJets_AK5PFclean = GetJets(RA4Muon,RA4Elec,RA4MuonVeto,RA4ElecVeto,
                                                  HT,LooseBJet,MediumBJet, 
                                                  2.4, 30, 0.3); 
        for(int i=0; i<GoodJets_AK5PFclean.size(); i++) cout << event << " :: " << GoodJets_AK5PFclean.at(i) << endl;
       

        h1_result->Fill(0.5);

        /* 
        // variables 
        vector<float> Vector_mj;   // mj
        float MJ=0;
        int Nfastjets=0;
        for(int ifastjet=0; ifastjet<(int)fastjets_AK5PF_px_->size(); ifastjet++) {
            
            // Number of jets 
            float pT = TMath::Sqrt( fastjets_AK5PF_px_->at(ifastjet)*fastjets_AK5PF_px_->at(ifastjet)
                                   +fastjets_AK5PF_py_->at(ifastjet)*fastjets_AK5PF_py_->at(ifastjet));

            float temp_mj = Getmj(fastjets_AK5PF_px_->at(ifastjet), fastjets_AK5PF_py_->at(ifastjet),
                                  fastjets_AK5PF_pz_->at(ifastjet), fastjets_AK5PF_energy_->at(ifastjet));

            Vector_mj.push_back(temp_mj); 
            if(Npv_ > 25) {
                h1_mj_fastjets_Npv26toInf->Fill(temp_mj);
            } else if(Npv_>15) {
                h1_mj_fastjets_Npv16to25->Fill(temp_mj);
            } else {
                h1_mj_fastjets_Npv0to15->Fill(temp_mj);
            }
            Nfastjets++;
        }

        // Get MJ        
        MJ = GetMJ(Vector_mj);

        // fill histogram for MJ and Njets
        if(Npv_ > 25) {
            h1_njets_fastjets_Npv26toInf->Fill( TMath::Min((Float_t)Nfastjets,(Float_t)29.499) );
            h1_MJ_fastjets_Npv26toInf->Fill(MJ);
        } else if(Npv_>15) {
            h1_njets_fastjets_Npv16to25->Fill( TMath::Min((Float_t)Nfastjets,(Float_t)29.499) );
            h1_MJ_fastjets_Npv16to25->Fill(MJ);
        } else {
            h1_njets_fastjets_Npv0to15->Fill( TMath::Min((Float_t)Nfastjets,(Float_t)29.499) );
            h1_MJ_fastjets_Npv0to15->Fill(MJ);
        }
       */
    
    } // event loop
    cout << endl;
    cout << "... Looping events done" << endl;

    //
    // Write the result 
    //
    TString ResultFileName = ProcessName+".root";
    //TString ResultFileName = InputCfAName;
    //ResultFileName.ReplaceAll(".root", "");
    //ResultFileName.ReplaceAll("/", "_");
    //ResultFileName = "Result_"+ResultFileName+".root";
    //ResultFileName.ReplaceAll("__", "_");
    cout << "Writing " << ResultFileName << endl;
    TFile *ResultFile = new TFile(ResultFileName, "RECREATE");
    gROOT->cd();
    ResultFile->cd();
    h1_result->SetDirectory(0);  h1_result->Write();
    ResultFile->Close();
    
    TDatime DTEnd;
    int EndTime = DTEnd.GetTime(); 
    cout << "End time : " << EndTime << endl; 

    // cleanup
    delete chainA;
    delete chainB;
    delete ResultFile;
}
