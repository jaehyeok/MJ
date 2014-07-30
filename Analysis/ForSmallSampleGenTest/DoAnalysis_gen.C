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

int R0p5pTcut = 10;
int JetpTthres_mj = 0;
bool TEST = false;
char* jettype = Form("fastjets_AK5PFclean_gen_R1p2_R0p5pT%i",R0p5pTcut);

void h1cosmetic(TH1F* &h1, char* title, int linecolor=kRed, int linewidth=1, int fillcolor=0, TString var=""){

    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.1);

}

void h2cosmetic(TH2F* &h2, char* title, TString varX="x", TString varY="y"){

    h2->SetTitle(title);
    h2->SetXTitle(varX);
    h2->SetYTitle(varY);
    h2->SetStats(0);
    //h2->SetMinimum(0.1);

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
float GetMJscalar(vector<float> Vectormj){
    
    if(Vectormj.size()==0) {
   
        return -999.;
    
    } else {

        float MJ = 0.;
        for(int imj=0; imj<(int)Vectormj.size(); imj++) MJ = MJ + Vectormj.at(imj);
        return MJ;
    
    }
}

void DoAnalysis_gen_one(TString Input, 
                        TH1F* &h1_MJscalar, TH1F* &h1_MJvec, TH1F* &h1_Ratio, 
                        TH2F* &h2_MJscalarVsMJvec, TH2F* &h2_MJscalarVsMJRatio, 
                        float xsec) { 
    

    // 
    // Get tree 
    // 
    TChain * chainB = new TChain("/configurableAnalysis/eventB");   
    chainB->Add(Input);
    // 
    // set address of variables to read 
    // 

    // event info
    UInt_t   event_ = 0;
    chainB->SetBranchAddress("event", &event_);
    UInt_t   run_ = 0;
    chainB->SetBranchAddress("run", &run_);
    UInt_t   lumiblock_ = 0;
    chainB->SetBranchAddress("lumiblock", &lumiblock_);
    Int_t   Npv_ = 0;
    chainB->SetBranchAddress("Npv", &Npv_);

    // fastjets
    vector<float>   *fastjets_AK5PF_px_ = 0;
    chainB->SetBranchAddress(Form("%s_px",jettype), &fastjets_AK5PF_px_);
    vector<float>   *fastjets_AK5PF_py_ = 0;
    chainB->SetBranchAddress(Form("%s_py",jettype), &fastjets_AK5PF_py_);
    vector<float>   *fastjets_AK5PF_pz_ = 0;
    chainB->SetBranchAddress(Form("%s_pz",jettype), &fastjets_AK5PF_pz_);
    vector<float>   *fastjets_AK5PF_energy_ = 0;
    chainB->SetBranchAddress(Form("%s_energy",jettype), &fastjets_AK5PF_energy_);
    vector<float>   *fastjets_AK5PF_phi_ = 0;
    chainB->SetBranchAddress(Form("%s_phi",jettype), &fastjets_AK5PF_phi_);
    vector<float>   *fastjets_AK5PF_eta_ = 0;
    chainB->SetBranchAddress(Form("%s_eta",jettype), &fastjets_AK5PF_eta_);
    
    //
    // main event loop
    //
    Int_t nentries = (Int_t)chainB->GetEntries();
    cout<<"The number of entries is: "<<nentries<<endl;
    int i_permille_old = 0;

    float weight = xsec/nentries;

    for(int ib = 0; ib<nentries; ib++) {

        int i_permille = (int)floor(1000 * ib / float(nentries));
        if (i_permille != i_permille_old) {
            // xterm magic from L. Vacavant and A. Cerri
            if (isatty(1)) {
                printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                        "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                fflush(stdout);
            }
            i_permille_old = i_permille;
        }
    
        // get entry of an event 
        chainB->GetEntry(ib);
        
        // variables 
        vector<float> Vector_mj;   // mj
        float MJscalar=0;
        float MJvec=0;
        float px=0, py=0, pz=0, energy=0;
        for(int ifastjet=0; ifastjet<(int)fastjets_AK5PF_px_->size(); ifastjet++) {
            
            // Number of jets 
            float pT = TMath::Sqrt( fastjets_AK5PF_px_->at(ifastjet)*fastjets_AK5PF_px_->at(ifastjet)
                                   +fastjets_AK5PF_py_->at(ifastjet)*fastjets_AK5PF_py_->at(ifastjet));

            // Get mj
            if(pT<JetpTthres_mj) continue;  // pT cut for jets to calculate mj with

            float temp_mj = Getmj(fastjets_AK5PF_px_->at(ifastjet), fastjets_AK5PF_py_->at(ifastjet),
                                  fastjets_AK5PF_pz_->at(ifastjet), fastjets_AK5PF_energy_->at(ifastjet));
           
            Vector_mj.push_back(temp_mj);

            px = px + fastjets_AK5PF_px_->at(ifastjet);  
            py = py + fastjets_AK5PF_py_->at(ifastjet);  
            pz = pz + fastjets_AK5PF_pz_->at(ifastjet);  
            energy = energy + fastjets_AK5PF_energy_->at(ifastjet);  
        }

        // Get MJ        
        MJscalar = GetMJscalar(Vector_mj);
        h1_MJscalar->Fill(MJscalar, weight);
        MJvec = TMath::Sqrt(energy*energy - px*px - py*py - pz*pz);
        h1_MJvec->Fill(MJvec, weight);
        h1_Ratio->Fill(MJscalar/MJvec, weight);
        h2_MJscalarVsMJvec->Fill(MJscalar, MJvec, weight);
        h2_MJscalarVsMJRatio->Fill(MJscalar, MJscalar/MJvec, weight);

        //cout << MJscalar << " " << MJvec << endl; // DEBUG
    } // event loop
}

void DoAnalysis_gen() {
    
    // TH1F
    TH1F *h1_MJscalar_qcd1000       = new TH1F("h1_MJscalar_qcd1000","h1_MJscalar_qcd1000", 50, 0, 2000);
    TH1F *h1_MJvec_qcd1000          = new TH1F("h1_MJvec_qcd1000","h1_MJvec_qcd1000", 50, 0, 6000);
    TH1F *h1_MJRatio_qcd1000        = new TH1F("h1_MJRatio_qcd1000","h1_MJRatio_qcd1000", 50, 0, 1.5);
    TH1F *h1_MJscalar_qcd500        = new TH1F("h1_MJscalar_qcd500","h1_MJscalar_qcd500", 50, 0, 2000);
    TH1F *h1_MJvec_qcd500           = new TH1F("h1_MJvec_qcd500","h1_MJvec_qcd500", 50, 0, 6000);
    TH1F *h1_MJRatio_qcd500         = new TH1F("h1_MJRatio_qcd500","h1_MJRatio_qcd500", 50, 0, 1.5);
    TH1F *h1_MJscalar_qcd250        = new TH1F("h1_MJscalar_qcd250","h1_MJscalar_qcd250", 50, 0, 2000);
    TH1F *h1_MJvec_qcd250           = new TH1F("h1_MJvec_qcd250","h1_MJvec_qcd250", 50, 0, 6000);
    TH1F *h1_MJRatio_qcd250         = new TH1F("h1_MJRatio_qcd250","h1_MJRatio_qcd250", 50, 0, 1.5);
    TH1F *h1_MJscalar_qcd100        = new TH1F("h1_MJscalar_qcd100","h1_MJscalar_qcd100", 50, 0, 2000);
    TH1F *h1_MJvec_qcd100           = new TH1F("h1_MJvec_qcd100","h1_MJvec_qcd100", 50, 0, 6000);
    TH1F *h1_MJRatio_qcd100         = new TH1F("h1_MJRatio_qcd100","h1_MJRatio_qcd100", 50, 0, 1.5);
    TH1F *h1_MJscalar_qcd           = new TH1F("h1_MJscalar_qcd","h1_MJscalar_qcd", 50, 0, 2000);
    TH1F *h1_MJvec_qcd              = new TH1F("h1_MJvec_qcd","h1_MJvec_qcd", 50, 0, 6000);
    TH1F *h1_MJRatio_qcd            = new TH1F("h1_MJRatio_qcd","h1_MJRatio_qcd", 50, 0, 1.5);
    TH1F *h1_MJscalar_ttbar         = new TH1F("h1_MJscalar_ttbar","h1_MJscalar_ttbar", 50, 0, 2000);
    TH1F *h1_MJvec_ttbar            = new TH1F("h1_MJvec_ttbar","h1_MJvec_ttbar", 50, 0, 6000);
    TH1F *h1_MJRatio_ttbar          = new TH1F("h1_MJRatio_ttbar","h1_MJRatio_ttbar", 50, 0, 1.5);
    TH1F *h1_MJscalar_sig           = new TH1F("h1_MJscalar_sig","h1_MJscalar_sig", 50, 0, 2000);
    TH1F *h1_MJvec_sig              = new TH1F("h1_MJvec_sig","h1_MJvec_sig", 50, 0, 6000);
    TH1F *h1_MJRatio_sig            = new TH1F("h1_MJRatio_sig","h1_MJRatio_sig", 50, 0, 1.5);
    TH1F *h1_MJscalar_sig2          = new TH1F("h1_MJscalar_sig2","h1_MJscalar_sig2", 50, 0, 2000);
    TH1F *h1_MJvec_sig2             = new TH1F("h1_MJvec_sig2","h1_MJvec_sig2", 50, 0, 6000);
    TH1F *h1_MJRatio_sig2           = new TH1F("h1_MJRatio_sig2","h1_MJRatio_sig2", 50, 0, 1.5);
    TH1F *h1_MJscalar_sig3          = new TH1F("h1_MJscalar_sig3","h1_MJscalar_sig3", 50, 0, 2000);
    TH1F *h1_MJvec_sig3             = new TH1F("h1_MJvec_sig3","h1_MJvec_sig3", 50, 0, 6000);
    TH1F *h1_MJRatio_sig3           = new TH1F("h1_MJRatio_sig3","h1_MJRatio_sig3", 50, 0, 1.5);
    
    // TH2F
    TH2F *h2_MJscalarVsMJvec_qcd1000    = new TH2F("h2_MJscalarVsMJvec_qcd1000","h2_MJscalarVsMJvec_qcd1000", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_qcd1000  = new TH2F("h2_MJscalarVsMJRatio_qcd1000","h2_MJscalarVsMJRatio_qcd1000", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_qcd500     = new TH2F("h2_MJscalarVsMJvec_qcd500","h2_MJscalarVsMJvec_qcd500", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_qcd500   = new TH2F("h2_MJscalarVsMJRatio_qcd500","h2_MJscalarVsMJRatio_qcd500", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_qcd250     = new TH2F("h2_MJscalarVsMJvec_qcd250","h2_MJscalarVsMJvec_qcd250", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_qcd250   = new TH2F("h2_MJscalarVsMJRatio_qcd250","h2_MJscalarVsMJRatio_qcd250", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_qcd100     = new TH2F("h2_MJscalarVsMJvec_qcd100","h2_MJscalarVsMJvec_qcd100", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_qcd100   = new TH2F("h2_MJscalarVsMJRatio_qcd100","h2_MJscalarVsMJRatio_qcd100", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_qcd        = new TH2F("h2_MJscalarVsMJvec_qcd","h2_MJscalarVsMJvec_qcd", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_qcd      = new TH2F("h2_MJscalarVsMJRatio_qcd","h2_MJscalarVsMJRatio_qcd", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_ttbar      = new TH2F("h2_MJscalarVsMJvec_ttbar","h2_MJscalarVsMJvec_ttbar", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_ttbar    = new TH2F("h2_MJscalarVsMJRatio_ttbar","h2_MJscalarVsMJRatio_ttbar", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_sig        = new TH2F("h2_MJscalarVsMJvec_sig","h2_MJscalarVsMJvec_sig", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_sig      = new TH2F("h2_MJscalarVsMJRatio_sig","h2_MJscalarVsMJRatio_sig", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_sig2       = new TH2F("h2_MJscalarVsMJvec_sig2","h2_MJscalarVsMJvec_sig2", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_sig2     = new TH2F("h2_MJscalarVsMJRatio_sig2","h2_MJscalarVsMJRatio_sig2", 25, 0, 2000, 30, 0, 1.5);
    TH2F *h2_MJscalarVsMJvec_sig3       = new TH2F("h2_MJscalarVsMJvec_sig3","h2_MJscalarVsMJvec_sig3", 25, 0, 2000, 30, 0, 6000);
    TH2F *h2_MJscalarVsMJRatio_sig3     = new TH2F("h2_MJscalarVsMJRatio_sig3","h2_MJscalarVsMJRatio_sig3", 25, 0, 2000, 30, 0, 1.5);
                        
    
    //
    if(TEST) {
        // QCD xsec ref : CMS AN AN-12-350
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f1*.root", 
                h1_MJscalar_qcd1000, h1_MJvec_qcd1000, h1_MJRatio_qcd1000, h2_MJscalarVsMJvec_qcd1000, h2_MJscalarVsMJRatio_qcd1000, 204);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f1*.root", 
                h1_MJscalar_qcd500, h1_MJvec_qcd500, h1_MJRatio_qcd500, h2_MJscalarVsMJvec_qcd500, h2_MJscalarVsMJRatio_qcd500, 8426);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f1*.root", 
                h1_MJscalar_qcd250, h1_MJvec_qcd250, h1_MJRatio_qcd250, h2_MJscalarVsMJvec_qcd250, h2_MJscalarVsMJRatio_qcd250, 276000);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f1*.root", 
                h1_MJscalar_qcd100, h1_MJvec_qcd100, h1_MJRatio_qcd100, h2_MJscalarVsMJvec_qcd100, h2_MJscalarVsMJRatio_qcd100, 10360000);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_*_f1*.root", 
                h1_MJscalar_ttbar, h1_MJvec_ttbar, h1_MJRatio_ttbar, h2_MJscalarVsMJvec_ttbar, h2_MJscalarVsMJRatio_ttbar, 100);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_*_f1400_600_*1696*.root", // 2046
                h1_MJscalar_sig, h1_MJvec_sig, h1_MJRatio_sig, h2_MJscalarVsMJvec_sig, h2_MJscalarVsMJRatio_sig, 100);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-25to500_50GeVX50GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_*_f1400_175_*1070*.root", // 1732
                h1_MJscalar_sig2, h1_MJvec_sig2, h1_MJRatio_sig2, h2_MJscalarVsMJvec_sig2, h2_MJscalarVsMJRatio_sig2, 100);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-400to750_mLSP-25to550_50GeVX50GeV_Binning_Summer12-START52_V9_FSIM-v1_AODSIM_UCSB1733reshuf_v68_f700_175_0_717.root", // 1733 
                h1_MJscalar_sig3, h1_MJvec_sig3, h1_MJRatio_sig3, h2_MJscalarVsMJvec_sig3, h2_MJscalarVsMJRatio_sig3, 100);  
    } else {
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f*.root", 
                h1_MJscalar_qcd1000, h1_MJvec_qcd1000, h1_MJRatio_qcd1000, h2_MJscalarVsMJvec_qcd1000, h2_MJscalarVsMJRatio_qcd1000, 204);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f*.root", 
                h1_MJscalar_qcd500, h1_MJvec_qcd500, h1_MJRatio_qcd500, h2_MJscalarVsMJvec_qcd500, h2_MJscalarVsMJRatio_qcd500, 8426);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f*.root", 
                h1_MJscalar_qcd250, h1_MJvec_qcd250, h1_MJRatio_qcd250, h2_MJscalarVsMJvec_qcd250, h2_MJscalarVsMJRatio_qcd250, 276000);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_*_f*.root", 
                h1_MJscalar_qcd100, h1_MJvec_qcd100, h1_MJRatio_qcd100, h2_MJscalarVsMJvec_qcd100, h2_MJscalarVsMJRatio_qcd100, 10360000);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_*_f*.root", 
                h1_MJscalar_ttbar, h1_MJvec_ttbar, h1_MJRatio_ttbar, h2_MJscalarVsMJvec_ttbar, h2_MJscalarVsMJRatio_ttbar, 100);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_*_f1400_600_*.root", // 2046
                h1_MJscalar_sig, h1_MJvec_sig, h1_MJRatio_sig, h2_MJscalarVsMJvec_sig, h2_MJscalarVsMJRatio_sig, 100);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-25to500_50GeVX50GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_*_f1400_175_*.root", // 1732
                h1_MJscalar_sig2, h1_MJvec_sig2, h1_MJRatio_sig2, h2_MJscalarVsMJvec_sig2, h2_MJscalarVsMJRatio_sig2, 100);  
        DoAnalysis_gen_one("/Users/jaehyeok/scratch/genTest/cfA_SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-400to750_mLSP-25to550_50GeVX50GeV_Binning_Summer12-START52_V9_FSIM-v1_AODSIM_UCSB1733reshuf_v68_f700_175_*.root",  // 1733
                h1_MJscalar_sig3, h1_MJvec_sig3, h1_MJRatio_sig3, h2_MJscalarVsMJvec_sig3, h2_MJscalarVsMJRatio_sig3, 100);  
    }

    h1_MJscalar_qcd->Add(h1_MJscalar_qcd1000);
    h1_MJscalar_qcd->Add(h1_MJscalar_qcd500);
    h1_MJscalar_qcd->Add(h1_MJscalar_qcd250);
    h1_MJscalar_qcd->Add(h1_MJscalar_qcd100);
    h1_MJvec_qcd->Add(h1_MJvec_qcd1000);
    h1_MJvec_qcd->Add(h1_MJvec_qcd500);
    h1_MJvec_qcd->Add(h1_MJvec_qcd250);
    h1_MJvec_qcd->Add(h1_MJvec_qcd100);
    h1_MJRatio_qcd->Add(h1_MJRatio_qcd1000);
    h1_MJRatio_qcd->Add(h1_MJRatio_qcd500);
    h1_MJRatio_qcd->Add(h1_MJRatio_qcd250);
    h1_MJRatio_qcd->Add(h1_MJRatio_qcd100);
    h2_MJscalarVsMJvec_qcd->Add(h2_MJscalarVsMJvec_qcd1000);
    h2_MJscalarVsMJvec_qcd->Add(h2_MJscalarVsMJvec_qcd500);
    h2_MJscalarVsMJvec_qcd->Add(h2_MJscalarVsMJvec_qcd250);
    h2_MJscalarVsMJvec_qcd->Add(h2_MJscalarVsMJvec_qcd100);
    h2_MJscalarVsMJRatio_qcd->Add(h2_MJscalarVsMJRatio_qcd1000);
    h2_MJscalarVsMJRatio_qcd->Add(h2_MJscalarVsMJRatio_qcd500);
    h2_MJscalarVsMJRatio_qcd->Add(h2_MJscalarVsMJRatio_qcd250);
    h2_MJscalarVsMJRatio_qcd->Add(h2_MJscalarVsMJRatio_qcd100);

    //
    // cosmetics for histograms
    //
    h1cosmetic(h1_MJscalar_qcd1000, "", kBlue, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_qcd1000, "", kBlue, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_qcd1000, "", kBlue, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_qcd500, "", kBlue, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_qcd500, "", kBlue, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_qcd500, "", kBlue, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_qcd250, "", kBlue, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_qcd250, "", kBlue, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_qcd250, "", kBlue, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_qcd100, "", kBlue, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_qcd100, "", kBlue, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_qcd100, "", kBlue, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_qcd, "", kBlue, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_qcd, "", kBlue, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_qcd, "", kBlue, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_ttbar, "", kBlack, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_ttbar, "", kBlack, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_ttbar, "", kBlack, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_sig, "", kRed, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_sig, "", kRed, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_sig, "", kRed, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_sig2, "", kOrange, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_sig2, "", kOrange, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_sig2, "", kOrange, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    h1cosmetic(h1_MJscalar_sig3, "", kGreen, 2, 0, "MJ(scalar sum)");
    h1cosmetic(h1_MJvec_sig3, "", kGreen, 2, 0, "MJ(vector sum)");
    h1cosmetic(h1_MJRatio_sig3, "", kGreen, 2, 0, "MJ(scalar sum)/MJ(vector sum)");
    
    h2cosmetic(h2_MJscalarVsMJvec_qcd, "QCD", "MJ(scalar sum)", "MJ(vector sum)");
    h2cosmetic(h2_MJscalarVsMJRatio_qcd, "QCD", "MJ(scalar sum)", "R");
    h2cosmetic(h2_MJscalarVsMJvec_ttbar, "TT", "MJ(scalar sum)", "MJ(vector sum)");
    h2cosmetic(h2_MJscalarVsMJRatio_ttbar, "TT", "MJ(scalar sum)", "R");
    h2cosmetic(h2_MJscalarVsMJvec_sig, "m_{#tilde{g}}=1400 GeV m_{LSP}=575 GeV", "MJ(scalar sum)", "MJ(vector sum)");
    h2cosmetic(h2_MJscalarVsMJRatio_sig, "m_{#tilde{g}}=1400 GeV m_{LSP}=575 GeV", "MJ(scalar sum)", "R");
    h2cosmetic(h2_MJscalarVsMJvec_sig2, "m_{#tilde{g}}=1400 GeV m_{LSP}=175 GeV", "MJ(scalar sum)", "MJ(vector sum)");
    h2cosmetic(h2_MJscalarVsMJRatio_sig2, "m_{#tilde{g}}=1400 GeV m_{LSP}=175 GeV", "MJ(scalar sum)", "R");
    h2cosmetic(h2_MJscalarVsMJvec_sig3, "m_{#tilde{g}}=700 GeV m_{LSP}=175 GeV", "MJ(scalar sum)", "MJ(vector sum)");
    h2cosmetic(h2_MJscalarVsMJRatio_sig3, "m_{#tilde{g}}=700 GeV m_{LSP}=175 GeV", "MJ(scalar sum)", "R");
    
    // 
    // Legend 
    // 
    TLegend *l1 = new TLegend(0.55, 0.66, 0.85, 0.85);
    l1->SetFillColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1_MJscalar_qcd1000,        "QCD",   "l");
    l1->AddEntry(h1_MJscalar_ttbar,        "TT",   "l");
    l1->AddEntry(h1_MJscalar_sig3,       "m_{#tilde{g}}=700 GeV m_{LSP}=175 GeV",   "l");
    l1->AddEntry(h1_MJscalar_sig2,       "m_{#tilde{g}}=1400 GeV m_{LSP}=175 GeV",   "l");
    l1->AddEntry(h1_MJscalar_sig,        "m_{#tilde{g}}=1400 GeV m_{LSP}=575 GeV",   "l");
    
    // 
    // Canvas
    // 
    TCanvas *c = new TCanvas("c","c",800,800);
    c->Divide(2,2);
    c->cd(1);
    c->cd(1)->SetLogy(1);
    h1_MJscalar_qcd->SetMaximum(h1_MJscalar_qcd->GetMaximum()*5);
    h1_MJscalar_qcd->DrawNormalized("HIST");
    h1_MJscalar_ttbar->DrawNormalized("HIST SAME");
    h1_MJscalar_sig->DrawNormalized("HIST SAME");
    h1_MJscalar_sig2->DrawNormalized("HIST SAME");
    h1_MJscalar_sig3->DrawNormalized("HIST SAME");
    //h1_MJscalar_qcd1000->DrawNormalized("HIST SAME");
    //h1_MJscalar_qcd500->DrawNormalized("HIST SAME");
    //h1_MJscalar_qcd250->DrawNormalized("HIST SAME");
    l1->Draw();
    c->cd(2);
    c->cd(2)->SetLogy(1);
    h1_MJvec_qcd->SetMaximum(h1_MJvec_qcd->GetMaximum()*5);
    h1_MJvec_qcd->DrawNormalized("HIST");
    h1_MJvec_ttbar->DrawNormalized("HIST SAME");
    h1_MJvec_sig->DrawNormalized("HIST SAME");
    h1_MJvec_sig2->DrawNormalized("HIST SAME");
    h1_MJvec_sig3->DrawNormalized("HIST SAME");
    //h1_MJvec_qcd1000->DrawNormalized("HIST SAME");
    //h1_MJvec_qcd500->DrawNormalized("HIST SAME");
    //h1_MJvec_qcd250->DrawNormalized("HIST SAME");
    c->cd(3);
    c->cd(3)->SetLogy(1);
    h1_MJRatio_qcd->SetMaximum(h1_MJRatio_qcd->GetMaximum()*5);
    h1_MJRatio_qcd->DrawNormalized("HIST");
    h1_MJRatio_ttbar->DrawNormalized("HIST SAME");
    h1_MJRatio_sig->DrawNormalized("HIST SAME");
    h1_MJRatio_sig2->DrawNormalized("HIST SAME");
    h1_MJRatio_sig3->DrawNormalized("HIST SAME");
    
    c->SaveAs(Form("MJvec_MJscalar_R0p5pTcut%i_gen.pdf",R0p5pTcut)); 
    
    TCanvas *c2D = new TCanvas("c2D","c2D",1500,600);
    c2D->Divide(5,2);
    c2D->cd(1);
    h2_MJscalarVsMJvec_qcd->Draw("colz");
    c2D->cd(2);
    h2_MJscalarVsMJvec_ttbar->Draw("colz");
    c2D->cd(3);
    h2_MJscalarVsMJvec_sig->Draw("colz");
    c2D->cd(4);
    h2_MJscalarVsMJvec_sig2->Draw("colz");
    c2D->cd(5);
    h2_MJscalarVsMJvec_sig3->Draw("colz");
    c2D->cd(6);
    h2_MJscalarVsMJRatio_qcd->Draw("colz");
    c2D->cd(7);
    h2_MJscalarVsMJRatio_ttbar->Draw("colz");
    c2D->cd(8);
    h2_MJscalarVsMJRatio_sig->Draw("colz");
    c2D->cd(9);
    h2_MJscalarVsMJRatio_sig2->Draw("colz");
    c2D->cd(10);
    h2_MJscalarVsMJRatio_sig3->Draw("colz");
    
    c2D->SaveAs(Form("MJvec_MJscalar_2D_R0p5pTcut%i_gen.pdf",R0p5pTcut)); 
    
    // cleanup
    //delete f;
}
