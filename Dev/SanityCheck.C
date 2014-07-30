#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"

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

void SanityCheck() { 
    
    // 
    // Get tree 
    // 
    //TFile *f = TFile::Open("cfA/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2029_v71_f1_1_YWH_Test.root", "READ");
    //TFile *f = TFile::Open("cfA/cfA_JetHT_Run2012D-PromptReco-v1_AOD_UCSB2030_v71_f1_1_fip_Test.root", "READ");
    TFile *f = TFile::Open("cfA/NopT10Cut/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2048_v71_f1_1_DGh_Test.root", "READ");
    //TFile *f = TFile::Open("cfA/NopT10Cut/cfA_JetHT_Run2012D-PromptReco-v1_AOD_UCSB2049_v71_f1000_1_473_Test.root", "READ");
    TDirectory* dir = f->GetDirectory("configurableAnalysis");
    dir->cd();
    TTree *eventB = (TTree*)dir->Get("eventB");
    
    // 
    // set address of variables to read 
    // 

    // event info
    UInt_t   event_ = 0;
    eventB->SetBranchAddress("event", &event_);
    UInt_t   run_ = 0;
    eventB->SetBranchAddress("run", &run_);
    UInt_t   lumiblock_ = 0;
    eventB->SetBranchAddress("lumiblock", &lumiblock_);
    // jets
    vector<float>   *jets_AK5PFclean_px_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_px", &jets_AK5PFclean_px_);
    vector<float>   *jets_AK5PFclean_py_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_py", &jets_AK5PFclean_py_);
    vector<float>   *jets_AK5PFclean_pz_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_pz", &jets_AK5PFclean_pz_);
    vector<float>   *jets_AK5PFclean_energy_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_energy", &jets_AK5PFclean_energy_);
    vector<float>   *jets_AK5PFclean_phi_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_phi", &jets_AK5PFclean_phi_);
    vector<float>   *jets_AK5PFclean_eta_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_eta", &jets_AK5PFclean_eta_);
    vector<float>   *jets_AK5PFclean_rawPt_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_rawPt", &jets_AK5PFclean_rawPt_);
    vector<float>   *jets_AK5PF_rawPt_ = 0;
    eventB->SetBranchAddress("jets_AK5PF_rawPt", &jets_AK5PF_rawPt_);
    // fastjets
    vector<float>   *fastjets_AK5PF_px_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_px", &fastjets_AK5PF_px_);
    vector<float>   *fastjets_AK5PF_py_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_py", &fastjets_AK5PF_py_);
    vector<float>   *fastjets_AK5PF_pz_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_pz", &fastjets_AK5PF_pz_);
    vector<float>   *fastjets_AK5PF_energy_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_energy", &fastjets_AK5PF_energy_);
    vector<float>   *fastjets_AK5PF_phi_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_phi", &fastjets_AK5PF_phi_);
    vector<float>   *fastjets_AK5PF_eta_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_eta", &fastjets_AK5PF_eta_);
   
    // 
    //  loop over entries
    // 
    Int_t nentries = (Int_t)eventB->GetEntries();
    cout<<"The number of entries is: "<<nentries<<endl;

    TH1F *h1_pT_jets            = new TH1F("h1_pT_jets","h1_pT_jets", 100, 0, 1000);
    TH1F *h1_pT_fastjets        = new TH1F("h1_pT_fastjets","h1_pT_fastjets", 100, 0, 1000);
    TH1F *h1_eta_jets           = new TH1F("h1_eta_jets","h1_eta_jets", 100, -5, 5);
    TH1F *h1_eta_fastjets       = new TH1F("h1_eta_fastjets","h1_eta_fastjets", 100, -5, 5);
    TH1F *h1_phi_jets           = new TH1F("h1_phi_jets","h1_phi_jets", 100, -0.5, 0.5);
    TH1F *h1_phi_fastjets       = new TH1F("h1_phi_fastjets","h1_phi_fastjets", 100, -0.5, 0.5);
    TH1F *h1_njets_jets         = new TH1F("h1_njets_jets","h1_njets_jets", 15, -0.5, 14.5);
    TH1F *h1_njets_fastjets     = new TH1F("h1_njets_fastjets","h1_njets_fastjets", 15, -0.5, 14.5);
    
    TH1F *h1_pT_diff       = new TH1F("h1_pT_diff","h1_pT_diff", 100, -10, 10);
    //TH1F *h1_pT_diff       = new TH1F("h1_pT_diff","h1_pT_diff", 100, -0.5, 0.5);
    TH1F *h1_eta_diff      = new TH1F("h1_eta_diff","h1_eta_diff", 100, -0.1, 0.1);
    TH1F *h1_phi_diff      = new TH1F("h1_phi_diff","h1_phi_diff", 100, -0.1, 0.1);

    // pT threshold for Njets counting
    const float JetpTthres=30;

    // main event loop
    for(int ib = 0; ib<nentries; ib++) {
    //for(int ib = 0; ib<10; ib++) { // DEBUG
       
        if(ib%500==0) cout << "Entry : " << ib << endl; 

        // get entry of an event 
        eventB->GetEntry(ib);
        // 
        // Before matching 
        // 
        // Number of jets 
        int Njets=0;
        for(int ijet=0; ijet<(int)jets_AK5PFclean_rawPt_->size(); ijet++) {
            if(jets_AK5PFclean_rawPt_->at(ijet)>JetpTthres) Njets++;
        }
        int Nfastjets=0;
        for(int ifastjet=0; ifastjet<(int)fastjets_AK5PF_px_->size(); ifastjet++) {
            float pT = TMath::Sqrt( fastjets_AK5PF_px_->at(ifastjet)*fastjets_AK5PF_px_->at(ifastjet)
                                   +fastjets_AK5PF_py_->at(ifastjet)*fastjets_AK5PF_py_->at(ifastjet));
            if(pT>JetpTthres) Nfastjets++; 
        }
        
        // pT
        for(int ifastjet=0; ifastjet<(int)fastjets_AK5PF_px_->size(); ifastjet++) {
            float fastpT = TMath::Sqrt( fastjets_AK5PF_px_->at(ifastjet)*fastjets_AK5PF_px_->at(ifastjet)
                                       +fastjets_AK5PF_py_->at(ifastjet)*fastjets_AK5PF_py_->at(ifastjet));   
                h1_pT_fastjets->Fill( fastpT );
        }
        
        for(int ijet=0; ijet<(int)jets_AK5PFclean_rawPt_->size(); ijet++) {
                h1_pT_jets->Fill( jets_AK5PFclean_rawPt_->at(ijet) );
        }
        
        //
        // After matching 
        //
        for(int ifastjet=0; ifastjet<(int)fastjets_AK5PF_px_->size(); ifastjet++) {
            float fastpT = TMath::Sqrt( fastjets_AK5PF_px_->at(ifastjet)*fastjets_AK5PF_px_->at(ifastjet)
                                       +fastjets_AK5PF_py_->at(ifastjet)*fastjets_AK5PF_py_->at(ifastjet));   
            bool matched = false;
            int imatchedjet = 0;
            int imatchedetajet = 0;
            int imatchedphijet = 0;
            float detamin = 999.;
            float deta = 999.;
            float dphimin = 999.;
            float dphi = 999.;
            for(int ijet=0; ijet<(int)jets_AK5PFclean_rawPt_->size(); ijet++) {
                deta = TMath::Abs(jets_AK5PFclean_eta_->at(ijet) - fastjets_AK5PF_eta_->at(ifastjet));
                dphi = TMath::Abs(jets_AK5PFclean_phi_->at(ijet) - fastjets_AK5PF_phi_->at(ifastjet));

                if(deta<detamin){  
                    imatchedetajet = ijet;
                    detamin = deta;
                }
                if(dphi<dphimin){  
                    imatchedphijet = ijet;
                    dphimin = dphi;
                }

            }
            
            if( imatchedetajet==imatchedphijet //&& 
                //dphimin < 0.01 && detamin < 0.05
              ) { 
                matched = true; 
                imatchedjet=imatchedphijet;
            }

            if (matched) {  
                h1_eta_jets->Fill( jets_AK5PFclean_eta_->at(imatchedjet) );
                h1_eta_fastjets->Fill( fastjets_AK5PF_eta_->at(ifastjet) );
                h1_phi_jets->Fill( jets_AK5PFclean_phi_->at(imatchedjet) );
                h1_phi_fastjets->Fill( fastjets_AK5PF_phi_->at(ifastjet) );
                //h1_pT_diff->Fill( (fastpT-jets_AK5PFclean_rawPt_->at(imatchedjet))/jets_AK5PFclean_rawPt_->at(imatchedjet) );
                h1_pT_diff->Fill( (fastpT-jets_AK5PFclean_rawPt_->at(imatchedjet))/1.);
                h1_eta_diff->Fill( (fastjets_AK5PF_eta_->at(ifastjet)-jets_AK5PFclean_eta_->at(imatchedjet))/TMath::Abs(jets_AK5PFclean_eta_->at(imatchedjet)));
                h1_phi_diff->Fill( (fastjets_AK5PF_phi_->at(ifastjet)-jets_AK5PFclean_phi_->at(imatchedjet))/jets_AK5PFclean_phi_->at(imatchedjet) );

            } else {
                h1_eta_jets->Fill( -999. );
                h1_eta_fastjets->Fill( fastjets_AK5PF_eta_->at(ifastjet) );
                h1_phi_jets->Fill( -999. );
                h1_phi_fastjets->Fill( fastjets_AK5PF_phi_->at(ifastjet) );
                h1_pT_diff->Fill( -999. );
                h1_eta_diff->Fill( -999. );
                h1_phi_diff->Fill( -999.); 
            }
        } //if(matched)

        // njets needs to be outside of for loop
        h1_njets_jets->Fill( TMath::Min((Float_t)Njets,(Float_t)14.499) );
        h1_njets_fastjets->Fill( TMath::Min((Float_t)Nfastjets,(Float_t)14.999) );
    
    } // event loop

    h1cosmetic(h1_pT_jets, (char*)"pT", kBlack, 1, 0, "pT (GeV)");
    h1cosmetic(h1_pT_fastjets, (char*)"CMS fastjets", kRed, 1, 0, "pT (GeV)");
    h1cosmetic(h1_njets_jets, Form("Njets pT>%i GeV",(int)JetpTthres), kBlack, 1, 0, "Njets");
    h1cosmetic(h1_njets_fastjets, (char*)"Njets", kRed, 1, 0, "Njets");
    h1cosmetic(h1_pT_diff, (char*)"pT difference(Fastjet - cfA)", kBlack, 1, 0, "#Delta pT (GeV)");
    
    TCanvas *c = new TCanvas("c","c",1000,500);
    c->Divide(2,1);
    c->cd(1);
    c->cd(1)->SetLogy(1);
    h1_pT_jets->SetMaximum(h1_pT_jets->GetMaximum()*2);
    h1_pT_jets->Draw();
    h1_pT_fastjets->Draw("HIST SAME");
    c->cd(2);
    h1_pT_diff->SetMaximum(h1_pT_diff->GetMaximum()*1.2);
    h1_pT_diff->Draw();
    h1_pT_diff->Draw("HIST SAME");
    c->SaveAs(Form("SanityCheckpT.pdf")); 

    // njet distribution
    TCanvas *c_njets = new TCanvas();
    c_njets->cd(1);
    h1_njets_jets->SetMaximum(h1_njets_jets->GetMaximum()*1.4);
    h1_njets_jets->Draw();
    h1_njets_fastjets->Draw("HIST SAME");
    c_njets->SaveAs(Form("SanityCheckNjetsPT%i.pdf", (int)JetpTthres)); 
   
    // cleanup
    //delete f;
}
