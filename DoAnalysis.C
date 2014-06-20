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
    TChain * chainB_QCD = new TChain("/configurableAnalysis/eventB");   
    chainB_QCD->Add("../cfA/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2029_v71_*_TestR1p2.root");
    //chainB_QCD->Add("../cfA/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2029_v71_f2_1_71V_Test.root");

    // 
    // set address of variables to read 
    // 

    // event info
    UInt_t   event_ = 0;
    chainB_QCD->SetBranchAddress("event", &event_);
    UInt_t   run_ = 0;
    chainB_QCD->SetBranchAddress("run", &run_);
    UInt_t   lumiblock_ = 0;
    chainB_QCD->SetBranchAddress("lumiblock", &lumiblock_);
    Int_t   Npv_ = 0;
    chainB_QCD->SetBranchAddress("Npv", &Npv_);

    // fastjets
    vector<float>   *fastjets_AK5PF_px_ = 0;
    chainB_QCD->SetBranchAddress("fastjets_AK5PF_px", &fastjets_AK5PF_px_);
    vector<float>   *fastjets_AK5PF_py_ = 0;
    chainB_QCD->SetBranchAddress("fastjets_AK5PF_py", &fastjets_AK5PF_py_);
    vector<float>   *fastjets_AK5PF_pz_ = 0;
    chainB_QCD->SetBranchAddress("fastjets_AK5PF_pz", &fastjets_AK5PF_pz_);
    vector<float>   *fastjets_AK5PF_energy_ = 0;
    chainB_QCD->SetBranchAddress("fastjets_AK5PF_energy", &fastjets_AK5PF_energy_);
    vector<float>   *fastjets_AK5PF_phi_ = 0;
    chainB_QCD->SetBranchAddress("fastjets_AK5PF_phi", &fastjets_AK5PF_phi_);
    vector<float>   *fastjets_AK5PF_eta_ = 0;
    chainB_QCD->SetBranchAddress("fastjets_AK5PF_eta", &fastjets_AK5PF_eta_);

    // 
    // histograms
    // 
    TH1F *h1_mj_fastjets_Npv0to15       = new TH1F("h1_mj_fastjets_Npv0to15","h1_mj_fastjets_Npv0to15", 50, 0, 1500);
    TH1F *h1_MJ_fastjets_Npv0to15       = new TH1F("h1_MJ_fastjets_Npv0to15","h1_MJ_fastjets_Npv0to15", 50, 0, 3000);
    TH1F *h1_mj_fastjets_Npv16to25      = new TH1F("h1_mj_fastjets_Npv16to25","h1_mj_fastjets_Npv16to25", 50, 0, 1500);
    TH1F *h1_MJ_fastjets_Npv16to25      = new TH1F("h1_MJ_fastjets_Npv16to25","h1_MJ_fastjets_Npv16to25", 50, 0, 3000);
    TH1F *h1_mj_fastjets_Npv26toInf     = new TH1F("h1_mj_fastjets_Npv26toInf","h1_mj_fastjets_Npv26toInf", 50, 0, 1500);
    TH1F *h1_MJ_fastjets_Npv26toInf     = new TH1F("h1_MJ_fastjets_Npv26toInf","h1_MJ_fastjets_Npv26toInf", 50, 0, 3000);
    TH1F *h1_njets_fastjets_Npv0to15    = new TH1F("h1_njets_fastjets_Npv0to15","h1_njets_fastjets_Npv0to15", 30, -0.5, 29.5);
    TH1F *h1_njets_fastjets_Npv16to25   = new TH1F("h1_njets_fastjets_Npv16to25","h1_njets_fastjets_Npv16to25", 30, -0.5, 29.5);
    TH1F *h1_njets_fastjets_Npv26toInf  = new TH1F("h1_njets_fastjets_Npv26toInf","h1_njets_fastjets_Npv26toInf", 30, -0.5, 29.5);

    //
    // pT threshold for Njets counting
    //
    const float JetpTthres_mj=300;

    //
    // main event loop
    //
    Int_t nentries = (Int_t)chainB_QCD->GetEntries();
    cout<<"The number of entries is: "<<nentries<<endl;

    for(int ib = 0; ib<nentries; ib++) {
        //for(int ib = 0; ib<10; ib++) { // DEBUG

        if(ib%500==0) cout << "Entry : " << ib << endl; 

        // get entry of an event 
        chainB_QCD->GetEntry(ib);

        // variables 
        vector<float> Vector_mj;   // mj
        float MJ=0;
        int Nfastjets=0;
        for(int ifastjet=0; ifastjet<(int)fastjets_AK5PF_px_->size(); ifastjet++) {
            // Number of jets 
            float pT = TMath::Sqrt( fastjets_AK5PF_px_->at(ifastjet)*fastjets_AK5PF_px_->at(ifastjet)
                    +fastjets_AK5PF_py_->at(ifastjet)*fastjets_AK5PF_py_->at(ifastjet));

            // Get mj
            if(pT<JetpTthres_mj) continue;  // pT cut for jets to calculate mj with

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


    } // event loop

    //
    // cosmetics for histograms
    //
    //void h1cosmetic(TH1F* &h1, char* title, int linecolor=kRed, int linewidth=1, int fillcolor=0, TString var=""){
    // mj 
    h1cosmetic(h1_mj_fastjets_Npv0to15, Form("mj (pT>%i GeV)", (int)JetpTthres_mj), kBlack, 2, 0, "m_{j} [GeV]");
    h1cosmetic(h1_mj_fastjets_Npv16to25, Form("mj (pT>%i GeV)", (int)JetpTthres_mj), kRed, 2, 0, "m_{j} [GeV]");
    h1cosmetic(h1_mj_fastjets_Npv26toInf, Form("mj (pT>%i GeV)", (int)JetpTthres_mj), kBlue, 2, 0, "m_{j} [GeV]");
    // MJ 
    h1cosmetic(h1_MJ_fastjets_Npv0to15, Form("MJ (pT>%i GeV)", (int)JetpTthres_mj), kBlack, 2, 0, "M_{J} [GeV]");
    h1cosmetic(h1_MJ_fastjets_Npv16to25, Form("MJ (pT>%i GeV)", (int)JetpTthres_mj), kRed, 2, 0, "M_{J} [GeV]");
    h1cosmetic(h1_MJ_fastjets_Npv26toInf, Form("MJ (pT>%i GeV)", (int)JetpTthres_mj), kBlue, 2, 0, "M_{J} [GeV]");
    // Njets 
    h1cosmetic(h1_njets_fastjets_Npv0to15, Form("Njets (pT>%i GeV)", (int)JetpTthres_mj), kBlack, 2, 0, "Njets");
    h1cosmetic(h1_njets_fastjets_Npv16to25, Form("Njets (pT>%i GeV)", (int)JetpTthres_mj), kRed, 2, 0, "Njets");
    h1cosmetic(h1_njets_fastjets_Npv26toInf, Form("Njets (pT>%i GeV)", (int)JetpTthres_mj), kBlue, 2, 0, "Njets");

    // 
    // Legend 
    // 
    TLegend *l1 = new TLegend(0.5, 0.60, 0.85, 0.85);
    l1->SetFillColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1_mj_fastjets_Npv0to15,        "Npv = 0 - 15 ",   "l");
    l1->AddEntry(h1_mj_fastjets_Npv16to25,       "Npv = 16 - 25 ",  "l");
    l1->AddEntry(h1_mj_fastjets_Npv26toInf,      "Npv = 26 -  ",    "l");

    // 
    // Canvas
    // 
    TCanvas *c = new TCanvas("c","c",1200,400);
    c->Divide(3,1);
    c->cd(1);
    c->cd(1)->SetLogy(1);
    h1_mj_fastjets_Npv0to15->SetMaximum(h1_mj_fastjets_Npv0to15->GetMaximum()*2);
    h1_mj_fastjets_Npv0to15->DrawNormalized("HIST");
    h1_mj_fastjets_Npv16to25->DrawNormalized("HIST SAME");
    h1_mj_fastjets_Npv26toInf->DrawNormalized("HIST SAME");
    l1->Draw();
    c->cd(2);
    c->cd(2)->SetLogy(1);
    h1_MJ_fastjets_Npv0to15->SetMaximum(h1_MJ_fastjets_Npv0to15->GetMaximum()*2);
    h1_MJ_fastjets_Npv0to15->DrawNormalized("HIST");
    h1_MJ_fastjets_Npv16to25->DrawNormalized("HIST SAME");
    h1_MJ_fastjets_Npv26toInf->DrawNormalized("HIST SAME");
    c->cd(3);
    h1_njets_fastjets_Npv0to15->SetMaximum(h1_njets_fastjets_Npv0to15->GetMaximum()*2);
    h1_njets_fastjets_Npv0to15->DrawNormalized("HIST");
    h1_njets_fastjets_Npv16to25->DrawNormalized("HIST SAME");
    h1_njets_fastjets_Npv26toInf->DrawNormalized("HIST SAME");
    c->SaveAs(Form("QCDMC_R1p2_mj_MJ_Njet_pT%i.pdf", (int)JetpTthres_mj)); 

    // cleanup
    //delete f;
}
