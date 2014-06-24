#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TMath.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TLegend.h"

using namespace std;

void h1cosmetic(TH1F* &h1, char* title, int linecolor=kRed, int linewidth=1, int fillcolor=0, TString var=""){

    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.1);

}

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


void DoAnalysis() { 

    // 
    // Get tree 
    // 
    TChain * chainB_QCD_R1p2usingR0p5AK5PF        = new TChain("/configurableAnalysis/eventB");   
    TChain * chainB_QCD_R1p2usingR0p5AK5PFclean   = new TChain("/configurableAnalysis/eventB");   
    chainB_QCD_R1p2usingR0p5AK5PF->Add("../cfA/NopT10Cut/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2048_v71_*_TestR1p2usingR0p5jetsAK5PF.root");
    chainB_QCD_R1p2usingR0p5AK5PFclean->Add("../cfA/NopT10Cut/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2048_v71_*_TestR1p2usingR0p5jetsAK5PFclean.root");

    // 
    // set address of variables to read 
    // 
    
    // 
    // chainB_QCD 
    // 
    // event info
    UInt_t   event_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("event", &event_);
    UInt_t   run_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("run", &run_);
    UInt_t   lumiblock_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("lumiblock", &lumiblock_);
    Int_t   Npv_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("Npv", &Npv_);

    // 
    // chainB_QCD_R1p2usingR0p5AK5PF
    // 
    vector<float>   *jets_AK5PF_R1p2_px_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PF->SetBranchAddress("jets_AK5PF_R1p2_px", &jets_AK5PF_R1p2_px_);
    vector<float>   *jets_AK5PF_R1p2_py_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PF->SetBranchAddress("jets_AK5PF_R1p2_py", &jets_AK5PF_R1p2_py_);
    vector<float>   *jets_AK5PF_R1p2_pz_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PF->SetBranchAddress("jets_AK5PF_R1p2_pz", &jets_AK5PF_R1p2_pz_);
    vector<float>   *jets_AK5PF_R1p2_energy_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PF->SetBranchAddress("jets_AK5PF_R1p2_energy", &jets_AK5PF_R1p2_energy_);
    vector<float>   *jets_AK5PF_R1p2_phi_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PF->SetBranchAddress("jets_AK5PF_R1p2_phi", &jets_AK5PF_R1p2_phi_);
    vector<float>   *jets_AK5PF_R1p2_eta_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PF->SetBranchAddress("jets_AK5PF_R1p2_eta", &jets_AK5PF_R1p2_eta_);

    // 
    // chainB_QCD_R1p2usingR0p5AK5PFclean
    // 
    vector<float>   *jets_AK5PFclean_R1p2_px_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("jets_AK5PF_R1p2_px", &jets_AK5PFclean_R1p2_px_);
    vector<float>   *jets_AK5PFclean_R1p2_py_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("jets_AK5PF_R1p2_py", &jets_AK5PFclean_R1p2_py_);
    vector<float>   *jets_AK5PFclean_R1p2_pz_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("jets_AK5PF_R1p2_pz", &jets_AK5PFclean_R1p2_pz_);
    vector<float>   *jets_AK5PFclean_R1p2_energy_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("jets_AK5PF_R1p2_energy", &jets_AK5PFclean_R1p2_energy_);
    vector<float>   *jets_AK5PFclean_R1p2_phi_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("jets_AK5PF_R1p2_phi", &jets_AK5PFclean_R1p2_phi_);
    vector<float>   *jets_AK5PFclean_R1p2_eta_ = 0;
    chainB_QCD_R1p2usingR0p5AK5PFclean->SetBranchAddress("jets_AK5PF_R1p2_eta", &jets_AK5PFclean_R1p2_eta_);

    // 
    //  Histograms
    //  A : AK5PF
    //  B : AK5PFclean
    // 
    TH1F *h1_pt_jets_A      = new TH1F("h1_pt_jets_A","h1_pt_jets_A",           100, 0, 1000);
    TH1F *h1_pt_jets_B      = new TH1F("h1_pt_jets_B","h1_pt_jets_B",           100, 0, 1000);
    TH1F *h1_eta_jets_A     = new TH1F("h1_eta_jets_A","h1_eta_jets_A",         50, -5, 5);
    TH1F *h1_eta_jets_B     = new TH1F("h1_eta_jets_B","h1_eta_jets_B",         50, -5, 5);
    TH1F *h1_phi_jets_A     = new TH1F("h1_phi_jets_A","h1_phi_jets_A",         50, -3.141592, 3.141592);
    TH1F *h1_phi_jets_B     = new TH1F("h1_phi_jets_B","h1_phi_jets_B",         50, -3.141592, 3.141592);
    TH1F *h1_mj_jets_A      = new TH1F("h1_mj_jets_A","h1_mj_jets_A",           50, 0, 1500);
    TH1F *h1_mj_jets_B      = new TH1F("h1_mj_jets_B","h1_mj_jets_B",           50, 0, 1500);
    TH1F *h1_MJ_jets_A      = new TH1F("h1_MJ_jets_A","h1_MJ_jets_A",           50, 0, 3000);
    TH1F *h1_MJ_jets_B      = new TH1F("h1_MJ_jets_B","h1_MJ_jets_B",           50, 0, 3000);
    TH1F *h1_njets_jets_A   = new TH1F("h1_njets_jets_A","h1_njets_jets_A",     30, -0.5, 29.5);
    TH1F *h1_njets_jets_B   = new TH1F("h1_njets_jets_B","h1_njets_jets_B",     30, -0.5, 29.5);

    //
    // pT threshold for Njets counting
    //
    const float JetpTthres=50;

    //
    // main event loop
    //
    Int_t nentries = (Int_t)chainB_QCD_R1p2usingR0p5AK5PF->GetEntries();
    cout<<"The number of entries is: "<<nentries<<endl;

    for(int ib = 0; ib<nentries; ib++) {
        //for(int ib = 0; ib<10; ib++) { // DEBUG

        if(ib%500==0) cout << "Entry : " << ib << endl; 

        // get entry of an event 
        chainB_QCD_R1p2usingR0p5AK5PF->GetEntry(ib);
        chainB_QCD_R1p2usingR0p5AK5PFclean->GetEntry(ib); 
        
        //
        // Loop over jets and fill the histograms
        //
        vector<float> Vector_mj_A, Vector_mj_B;   // mj
        float MJ_A=0, MJ_B=0;
        int Njets_A=0, Njets_B=0;
        for(int ijet=0; ijet<(int)jets_AK5PF_R1p2_px_->size(); ijet++) {
            float pT_A = TMath::Sqrt( jets_AK5PF_R1p2_px_->at(ijet)
                                   *jets_AK5PF_R1p2_px_->at(ijet)
                                   +jets_AK5PF_R1p2_py_->at(ijet)
                                   *jets_AK5PF_R1p2_py_->at(ijet));
            if(pT_A<JetpTthres) continue;  // pT cut for jets to calculate mj with
            Njets_A++; 
            h1_pt_jets_A->Fill(pT_A);
            h1_eta_jets_A->Fill(jets_AK5PF_R1p2_eta_->at(ijet));
            h1_phi_jets_A->Fill(jets_AK5PF_R1p2_phi_->at(ijet));
            
            float temp_mj = Getmj(jets_AK5PF_R1p2_px_->at(ijet), jets_AK5PF_R1p2_py_->at(ijet),
                                  jets_AK5PF_R1p2_pz_->at(ijet), jets_AK5PF_R1p2_energy_->at(ijet));
            h1_mj_jets_A->Fill(temp_mj);
            Vector_mj_A.push_back(temp_mj); 
    
        }

        for(int ijet=0; ijet<(int)jets_AK5PFclean_R1p2_px_->size(); ijet++) {
            float pT_B = TMath::Sqrt( jets_AK5PFclean_R1p2_px_->at(ijet)
                                            *jets_AK5PFclean_R1p2_px_->at(ijet)
                                            +jets_AK5PFclean_R1p2_py_->at(ijet)
                                            *jets_AK5PFclean_R1p2_py_->at(ijet));
            if(pT_B<JetpTthres) continue;  // pT cut for jets to calculate mj with
            Njets_B++;
            h1_pt_jets_B->Fill(pT_B);
            h1_eta_jets_B->Fill(jets_AK5PFclean_R1p2_eta_->at(ijet));
            h1_phi_jets_B->Fill(jets_AK5PFclean_R1p2_phi_->at(ijet));
            
            float temp_mj = Getmj(jets_AK5PFclean_R1p2_px_->at(ijet), jets_AK5PFclean_R1p2_py_->at(ijet),
                                  jets_AK5PFclean_R1p2_pz_->at(ijet), jets_AK5PFclean_R1p2_energy_->at(ijet));
            h1_mj_jets_B->Fill(temp_mj);
            Vector_mj_B.push_back(temp_mj); 
        }

        // Get MJ
        MJ_A = GetMJ(Vector_mj_A);
        MJ_B = GetMJ(Vector_mj_B);
        
        // Njets
        h1_njets_jets_A->Fill( TMath::Min((Float_t)Njets_A,(Float_t)29.499) );
        h1_MJ_jets_A->Fill(MJ_A);
        h1_njets_jets_B->Fill( TMath::Min((Float_t)Njets_B,(Float_t)29.499) );
        h1_MJ_jets_B->Fill(MJ_B);

    } // event loop

    //
    // cosmetics for histograms
    //
    //void h1cosmetic(TH1F* &h1, char* title, int linecolor=kRed, int linewidth=1, int fillcolor=0, TString var=""){
    h1cosmetic(h1_pt_jets_A,     Form("pT >%i GeV",(int)JetpTthres), kBlack, 2, 0, "p_{T} [GeV]");
    h1cosmetic(h1_pt_jets_B,     Form("pT >%i GeV",(int)JetpTthres), kBlue, 2, 0, "p_{T} [GeV]");
    h1cosmetic(h1_eta_jets_A,    "Eta", kBlack, 2, 0, "#eta");
    h1cosmetic(h1_eta_jets_B,    "Eta", kBlue, 2, 0, "#eta");
    h1cosmetic(h1_phi_jets_A,    "Phi", kBlack, 2, 0, "#phi");
    h1cosmetic(h1_phi_jets_B,    "Phi", kBlue, 2, 0, "#phi");
    h1cosmetic(h1_mj_jets_A,     "mj", kBlack, 2, 0, "mj");
    h1cosmetic(h1_mj_jets_B,     "mj", kBlue, 2, 0, "#mj");
    h1cosmetic(h1_MJ_jets_A,     "MJ", kBlack, 2, 0, "MJ");
    h1cosmetic(h1_MJ_jets_B,     "MJ", kBlue, 2, 0, "#Mj");
    h1cosmetic(h1_njets_jets_A,  "Njets", kBlack, 2, 0, "Njets");
    h1cosmetic(h1_njets_jets_B,  "Njets", kBlue, 2, 0, "Njets");


    // 
    // Legend 
    // 
    TLegend *l1 = new TLegend(0.15, 0.60, 0.85, 0.85);
    l1->SetFillColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1_njets_jets_A,   "Using R=0.5 AK5PF branch",          "l");
    l1->AddEntry(h1_njets_jets_B,   "Using R=0.5 AK5PFclean branch",     "l");

    // 
    // Canvas
    // 
    TCanvas *c = new TCanvas("c","c",1200,800);
    c->Divide(3,2);
    c->cd(1);
    h1_pt_jets_A->SetMaximum(h1_pt_jets_A->GetMaximum()*2);
    h1_pt_jets_A->Draw("HIST");
    h1_pt_jets_B->Draw("HIST SAME");
    l1->Draw();
    c->cd(2);
    h1_eta_jets_A->SetMaximum(h1_eta_jets_A->GetMaximum()*2);
    h1_eta_jets_A->Draw("HIST");
    h1_eta_jets_B->Draw("HIST SAME");
    c->cd(3);
    h1_phi_jets_A->SetMaximum(h1_phi_jets_A->GetMaximum()*2);
    h1_phi_jets_A->Draw("HIST");
    h1_phi_jets_B->Draw("HIST SAME");
    c->cd(4)->SetLogy(1);
    h1_mj_jets_A->SetMaximum(h1_mj_jets_A->GetMaximum()*2);
    h1_mj_jets_A->Draw("HIST");
    h1_mj_jets_B->Draw("HIST SAME");
    c->cd(5)->SetLogy(1);
    h1_MJ_jets_A->SetMaximum(h1_MJ_jets_A->GetMaximum()*2);
    h1_MJ_jets_A->Draw("HIST");
    h1_MJ_jets_B->Draw("HIST SAME");
    c->cd(6);
    h1_njets_jets_A->SetMaximum(h1_njets_jets_A->GetMaximum()*2);
    h1_njets_jets_A->Draw("HIST");
    h1_njets_jets_B->Draw("HIST SAME");
    c->SaveAs(Form("QCDMC_pt_eta_phi_mj_MJ_Njet_AK5PFVsAK5PFclean.pdf", (int)JetpTthres)); 

    // cleanup
    //delete f;
}
