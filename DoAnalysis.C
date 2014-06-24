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
    TChain * chainB_QCD_R1p2                    = new TChain("/configurableAnalysis/eventB");   
    TChain * chainB_QCD_R1p2usingR0p5   = new TChain("/configurableAnalysis/eventB");   
    chainB_QCD_R1p2->Add("../cfA/NopT10Cut/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2048_v71_*_TestR1p2.root");
    chainB_QCD_R1p2usingR0p5->Add("../cfA/NopT10Cut/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2048_v71_*_TestR1p2usingR0p5jetsAK5PF.root");

    // 
    // set address of variables to read 
    // 
    
    // 
    // chainB_QCD 
    // 
    // event info
    UInt_t   event_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("event", &event_);
    UInt_t   run_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("run", &run_);
    UInt_t   lumiblock_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("lumiblock", &lumiblock_);
    Int_t   Npv_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("Npv", &Npv_);

    // chainB_QCD 
    vector<float>   *fastjets_AK5PF_R1p2_px_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("fastjets_AK5PF_px", &fastjets_AK5PF_R1p2_px_);
    vector<float>   *fastjets_AK5PF_R1p2_py_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("fastjets_AK5PF_py", &fastjets_AK5PF_R1p2_py_);
    vector<float>   *fastjets_AK5PF_R1p2_pz_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("fastjets_AK5PF_pz", &fastjets_AK5PF_R1p2_pz_);
    vector<float>   *fastjets_AK5PF_R1p2_energy_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("fastjets_AK5PF_energy", &fastjets_AK5PF_R1p2_energy_);
    vector<float>   *fastjets_AK5PF_R1p2_phi_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("fastjets_AK5PF_phi", &fastjets_AK5PF_R1p2_phi_);
    vector<float>   *fastjets_AK5PF_R1p2_eta_ = 0;
    chainB_QCD_R1p2->SetBranchAddress("fastjets_AK5PF_eta", &fastjets_AK5PF_R1p2_eta_);
    
    // 
    // chainB_QCD_R1p2usingR0p5
    // 
    vector<float>   *jets_AK5PF_R1p2_px_ = 0;
    chainB_QCD_R1p2usingR0p5->SetBranchAddress("jets_AK5PF_R1p2_px", &jets_AK5PF_R1p2_px_);
    vector<float>   *jets_AK5PF_R1p2_py_ = 0;
    chainB_QCD_R1p2usingR0p5->SetBranchAddress("jets_AK5PF_R1p2_py", &jets_AK5PF_R1p2_py_);
    vector<float>   *jets_AK5PF_R1p2_pz_ = 0;
    chainB_QCD_R1p2usingR0p5->SetBranchAddress("jets_AK5PF_R1p2_pz", &jets_AK5PF_R1p2_pz_);
    vector<float>   *jets_AK5PF_R1p2_energy_ = 0;
    chainB_QCD_R1p2usingR0p5->SetBranchAddress("jets_AK5PF_R1p2_energy", &jets_AK5PF_R1p2_energy_);
    vector<float>   *jets_AK5PF_R1p2_phi_ = 0;
    chainB_QCD_R1p2usingR0p5->SetBranchAddress("jets_AK5PF_R1p2_phi", &jets_AK5PF_R1p2_phi_);
    vector<float>   *jets_AK5PF_R1p2_eta_ = 0;
    chainB_QCD_R1p2usingR0p5->SetBranchAddress("jets_AK5PF_R1p2_eta", &jets_AK5PF_R1p2_eta_);

    // 
    // histograms
    // 
    TH1F *h1_pt_fastjets_R1p2               = new TH1F("h1_pt_fastjets_R1p2","h1_pt_fastjets_R1p2",                             100, 0, 1000);
    TH1F *h1_pt_fastjets_R1p2usingR0p5      = new TH1F("h1_pt_fastjets_R1p2usingR0p5","h1_pt_fastjets_R1p2usingR0p5",           100, 0, 1000);
    TH1F *h1_eta_fastjets_R1p2              = new TH1F("h1_eta_fastjets_R1p2","h1_eta_fastjets_R1p2",                           25, -5, 5);
    TH1F *h1_eta_fastjets_R1p2usingR0p5     = new TH1F("h1_eta_fastjets_R1p2usingR0p5","h1_eta_fastjets_R1p2usingR0p5",         25, -5, 5);
    TH1F *h1_phi_fastjets_R1p2              = new TH1F("h1_phi_fastjets_R1p2","h1_phi_fastjets_R1p2",                           25, -3.141592, 3.141592);
    TH1F *h1_phi_fastjets_R1p2usingR0p5     = new TH1F("h1_phi_fastjets_R1p2usingR0p5","h1_phi_fastjets_R1p2usingR0p5",         25, -3.141592, 3.141592);
    TH1F *h1_mj_fastjets_R1p2               = new TH1F("h1_mj_fastjets_R1p2","h1_mj_fastjets_R1p2",                             25, 0, 1500);
    TH1F *h1_mj_fastjets_R1p2usingR0p5      = new TH1F("h1_mj_fastjets_R1p2usingR0p5","h1_mj_fastjets_R1p2usingR0p5",           25, 0, 1500);
    TH1F *h1_MJ_fastjets_R1p2               = new TH1F("h1_MJ_fastjets_R1p2","h1_MJ_fastjets_R1p2",                             25, 0, 3000);
    TH1F *h1_MJ_fastjets_R1p2usingR0p5      = new TH1F("h1_MJ_fastjets_R1p2usingR0p5","h1_MJ_fastjets_R1p2usingR0p5",           25, 0, 3000);
    TH1F *h1_njets_fastjets_R1p2            = new TH1F("h1_njets_fastjets_R1p2","h1_njets_fastjets_R1p2",                   30, -0.5, 29.5);
    TH1F *h1_njets_fastjets_R1p2usingR0p5   = new TH1F("h1_njets_fastjets_R1p2usingR0p5","h1_njets_fastjets_R1p2usingR0p5", 30, -0.5, 29.5);

    //
    // pT threshold for Njets counting
    //
    const float JetpTthres=50;

    //
    // main event loop
    //
    Int_t nentries = (Int_t)chainB_QCD_R1p2->GetEntries();
    cout<<"The number of entries is: "<<nentries<<endl;

    for(int ib = 0; ib<nentries; ib++) {
        //for(int ib = 0; ib<10; ib++) { // DEBUG

        if(ib%500==0) cout << "Entry : " << ib << endl; 

        // get entry of an event 
        chainB_QCD_R1p2->GetEntry(ib);
        chainB_QCD_R1p2usingR0p5->GetEntry(ib); 
        
        //
        // Loop over jets and fill the histograms
        //
        vector<float> Vector_mj, Vector_mj_usingR0p5;   // mj
        float MJ=0, MJ_usingR0p5=0;
        int Nfastjets=0, Nfastjets_usingR0p5=0;
        for(int ifastjet=0; ifastjet<(int)fastjets_AK5PF_R1p2_px_->size(); ifastjet++) {
            float pT = TMath::Sqrt( fastjets_AK5PF_R1p2_px_->at(ifastjet)
                                   *fastjets_AK5PF_R1p2_px_->at(ifastjet)
                                   +fastjets_AK5PF_R1p2_py_->at(ifastjet)
                                   *fastjets_AK5PF_R1p2_py_->at(ifastjet));
            if(pT<JetpTthres) continue;  // pT cut for jets to calculate mj with
            Nfastjets++; 
            h1_pt_fastjets_R1p2->Fill(pT);
            h1_eta_fastjets_R1p2->Fill(fastjets_AK5PF_R1p2_eta_->at(ifastjet));
            h1_phi_fastjets_R1p2->Fill(fastjets_AK5PF_R1p2_phi_->at(ifastjet));
            
            float temp_mj = Getmj(fastjets_AK5PF_R1p2_px_->at(ifastjet), fastjets_AK5PF_R1p2_py_->at(ifastjet),
                                  fastjets_AK5PF_R1p2_pz_->at(ifastjet), fastjets_AK5PF_R1p2_energy_->at(ifastjet));
            h1_mj_fastjets_R1p2->Fill(temp_mj);
            Vector_mj.push_back(temp_mj); 
    
        }

        for(int ifastjet=0; ifastjet<(int)jets_AK5PF_R1p2_px_->size(); ifastjet++) {
            float pTusingR0p5 = TMath::Sqrt( jets_AK5PF_R1p2_px_->at(ifastjet)
                                            *jets_AK5PF_R1p2_px_->at(ifastjet)
                                            +jets_AK5PF_R1p2_py_->at(ifastjet)
                                            *jets_AK5PF_R1p2_py_->at(ifastjet));
            if(pTusingR0p5<JetpTthres) continue;  // pT cut for jets to calculate mj with
            Nfastjets_usingR0p5++;
            h1_pt_fastjets_R1p2usingR0p5->Fill(pTusingR0p5);
            h1_eta_fastjets_R1p2usingR0p5->Fill(jets_AK5PF_R1p2_eta_->at(ifastjet));
            h1_phi_fastjets_R1p2usingR0p5->Fill(jets_AK5PF_R1p2_phi_->at(ifastjet));
            
            float temp_mj = Getmj(jets_AK5PF_R1p2_px_->at(ifastjet), jets_AK5PF_R1p2_py_->at(ifastjet),
                                  jets_AK5PF_R1p2_pz_->at(ifastjet), jets_AK5PF_R1p2_energy_->at(ifastjet));
            h1_mj_fastjets_R1p2usingR0p5->Fill(temp_mj);
            Vector_mj_usingR0p5.push_back(temp_mj); 
        }

        // Get MJ
        MJ = GetMJ(Vector_mj);
        MJ_usingR0p5 = GetMJ(Vector_mj_usingR0p5);
        
        // Njets
        h1_njets_fastjets_R1p2->Fill( TMath::Min((Float_t)Nfastjets,(Float_t)29.499) );
        h1_MJ_fastjets_R1p2->Fill(MJ);
        h1_njets_fastjets_R1p2usingR0p5->Fill( TMath::Min((Float_t)Nfastjets_usingR0p5,(Float_t)29.499) );
        h1_MJ_fastjets_R1p2usingR0p5->Fill(MJ_usingR0p5);

    } // event loop

    //
    // cosmetics for histograms
    //
    //void h1cosmetic(TH1F* &h1, char* title, int linecolor=kRed, int linewidth=1, int fillcolor=0, TString var=""){
    h1cosmetic(h1_pt_fastjets_R1p2,             "pT", kBlack, 2, 0, "p_{T} [GeV]");
    h1cosmetic(h1_pt_fastjets_R1p2usingR0p5,    "pT", kBlue, 2, 0, "p_{T} [GeV]");
    h1cosmetic(h1_eta_fastjets_R1p2,            "Eta", kBlack, 2, 0, "#eta");
    h1cosmetic(h1_eta_fastjets_R1p2usingR0p5,   "Eta", kBlue, 2, 0, "#eta");
    h1cosmetic(h1_phi_fastjets_R1p2,            "Phi", kBlack, 2, 0, "#phi");
    h1cosmetic(h1_phi_fastjets_R1p2usingR0p5,   "Phi", kBlue, 2, 0, "#phi");
    h1cosmetic(h1_mj_fastjets_R1p2,            "mj", kBlack, 2, 0, "mj");
    h1cosmetic(h1_mj_fastjets_R1p2usingR0p5,   "mj", kBlue, 2, 0, "#mj");
    h1cosmetic(h1_MJ_fastjets_R1p2,            "MJ", kBlack, 2, 0, "MJ");
    h1cosmetic(h1_MJ_fastjets_R1p2usingR0p5,   "MJ", kBlue, 2, 0, "#Mj");
    h1cosmetic(h1_njets_fastjets_R1p2,         "Njets", kBlack, 2, 0, "Njets");
    h1cosmetic(h1_njets_fastjets_R1p2usingR0p5,"Njets", kBlue, 2, 0, "Njets");


    // 
    // Legend 
    // 
    TLegend *l1 = new TLegend(0.5, 0.60, 0.85, 0.85);
    l1->SetFillColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1_njets_fastjets_R1p2,            "Using PFcands",        "l");
    l1->AddEntry(h1_njets_fastjets_R1p2usingR0p5,   "Using R=0.5 jets",     "l");

    // 
    // Canvas
    // 
    TCanvas *c = new TCanvas("c","c",1200,800);
    c->Divide(3,2);
    c->cd(1);
    h1_pt_fastjets_R1p2->SetMaximum(h1_pt_fastjets_R1p2->GetMaximum()*2);
    h1_pt_fastjets_R1p2->Draw("HIST");
    h1_pt_fastjets_R1p2usingR0p5->Draw("HIST SAME");
    l1->Draw();
    c->cd(2);
    h1_eta_fastjets_R1p2->SetMaximum(h1_eta_fastjets_R1p2->GetMaximum()*2);
    h1_eta_fastjets_R1p2->Draw("HIST");
    h1_eta_fastjets_R1p2usingR0p5->Draw("HIST SAME");
    c->cd(3);
    h1_phi_fastjets_R1p2->SetMaximum(h1_phi_fastjets_R1p2->GetMaximum()*2);
    h1_phi_fastjets_R1p2->Draw("HIST");
    h1_phi_fastjets_R1p2usingR0p5->Draw("HIST SAME");
    c->cd(4)->SetLogy(1);
    h1_mj_fastjets_R1p2->SetMaximum(h1_mj_fastjets_R1p2->GetMaximum()*2);
    h1_mj_fastjets_R1p2->Draw("HIST");
    h1_mj_fastjets_R1p2usingR0p5->Draw("HIST SAME");
    c->cd(5)->SetLogy(1);
    h1_MJ_fastjets_R1p2->SetMaximum(h1_MJ_fastjets_R1p2->GetMaximum()*2);
    h1_MJ_fastjets_R1p2->Draw("HIST");
    h1_MJ_fastjets_R1p2usingR0p5->Draw("HIST SAME");
    c->cd(6);
    h1_njets_fastjets_R1p2->SetMaximum(h1_njets_fastjets_R1p2->GetMaximum()*2);
    h1_njets_fastjets_R1p2->Draw("HIST");
    h1_njets_fastjets_R1p2usingR0p5->Draw("HIST SAME");
    c->SaveAs(Form("QCDMC_pt_eta_phi_mj_MJ_Njet_PFcandsVsAK5PF.pdf", (int)JetpTthres)); 

    // cleanup
    //delete f;
}
