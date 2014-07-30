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
#include "TSystem.h"
#include "TLorentzVector.h"

//#include "Branch.h"
#include "slimmedBranch.h"
#include "ObjectSelector_Sync.h"
#include "Utilities.h"
#include "inJSON2012.h"

using namespace std;

//
// Get mj 
//
float Getmj(double px, double py, double pz, double E)
{

    float mj = TMath::Sqrt(E*E - px*px - py*py - pz*pz); 
    return mj;
}

//
// Get MJ 
//
float GetMJ(vector<float> Vectormj)
{
    
    if(Vectormj.size()>0) 
    {
        float MJ = 0.;
        for(int imj=0; imj<(int)Vectormj.size(); imj++)
        {    
            MJ = MJ + Vectormj.at(imj);
        }
        return MJ;
    } else 
    {
        return -999.;
    }
}

//
// main function
//
void DoOneProcess(TString InputName, TString ProcessName, int ibegin, int iend, bool isData, float Lumi) 
{
    //
    cout << "[MJ Analysis] ------------------------------------------------------------------------------------"<<endl; 
    cout << "[MJ Analysis] Processing : " << ProcessName << endl;
    cout << "[MJ Analysis] Input dir  : " << InputName << endl;

    TFile *babyFile_ = new TFile(Form("baby_%s_f%iTo%i.root", ProcessName.Data(), ibegin, iend), "RECREATE");
    babyFile_->cd();
    TTree *babyTree_ = new TTree("tree", "A Baby Ntuple");
    
    // 
    // Get tree 
    // 
    TChain * chainA = new TChain("/configurableAnalysis/eventA");   
    TChain * chainB = new TChain("/configurableAnalysis/eventB");  
    TChain * chainATotal = new TChain("/configurableAnalysis/eventB");  
    
    for(int i=ibegin; i<=iend; i++)  
    {
        gSystem->Exec(Form("ls %s/*f%i_*.root", InputName.Data(), i));
        chainA->Add(Form("%s/*f%i_*.root", InputName.Data(), i));
        chainB->Add(Form("%s/*f%i_*.root", InputName.Data(), i));
    } 
    TList *l = (TList*)chainA->GetListOfFiles();
    l->Print();

    InitializeA(chainA);
    InitializeB(chainB);
    
    // Get total number of events of a sample 
    chainATotal->Add(Form("%s/*.root", InputName.Data()));
    int TotalNEntries = (int)chainATotal->GetEntries();

    //
    // Baby variables 
    //
    int run_; 
    int lumiblock_; 
    int event_; 
    int Nfatjet_pT10_; 
    int Nfatjet_pT20_; 
    int Nfatjet_pT30_; 
    int Nskinnyjet_; 
    int Npv_; 
    float Npu_; 
    float EventWeight_; 
    float mll_; 
    float MJ_pT10_; 
    float MJ_pT20_; 
    float MJ_pT30_; 
    float HT_;  
    vector<float> mj_pT10_;
    vector<float> mj_pT20_;
    vector<float> mj_pT30_;
    vector<float> FatjetPt_pT10_;
    vector<float> FatjetEta_pT10_;
    vector<float> FatjetPhi_pT10_;
    vector<float> FatjetPt_pT20_;
    vector<float> FatjetEta_pT20_;
    vector<float> FatjetPhi_pT20_;
    vector<float> FatjetPt_pT30_;
    vector<float> FatjetEta_pT30_;
    vector<float> FatjetPhi_pT30_;
    vector<float> RA4leptonPt_;
    vector<float> RA4leptonEta_;
    vector<float> RA4leptonPhi_;
    vector<float> RA4leptonId_;
    
    babyTree_->Branch("run",            &run_);    
    babyTree_->Branch("lumiblock",      &lumiblock_); 
    babyTree_->Branch("event",          &event_);     
    babyTree_->Branch("Nfatjet_pT10",   &Nfatjet_pT10_);   
    babyTree_->Branch("Nfatjet_pT20",   &Nfatjet_pT20_);   
    babyTree_->Branch("Nfatjet_pT30",   &Nfatjet_pT30_);   
    babyTree_->Branch("Nskinnyjet",     &Nskinnyjet_);
    babyTree_->Branch("Npv",            &Npv_);       
    babyTree_->Branch("Npu",            &Npu_);       
    babyTree_->Branch("EventWeight",    &EventWeight_);
    babyTree_->Branch("mll",            &mll_);       
    babyTree_->Branch("MJ_pT10",        &MJ_pT10_);        
    babyTree_->Branch("MJ_pT20",        &MJ_pT20_);        
    babyTree_->Branch("MJ_pT30",        &MJ_pT30_);        
    babyTree_->Branch("HT",             &HT_);        
    babyTree_->Branch("mj_pT10",        &mj_pT10_);     
    babyTree_->Branch("mj_pT20",        &mj_pT20_);     
    babyTree_->Branch("mj_pT30",        &mj_pT30_);     
    babyTree_->Branch("FatjetPt_pT10",  &FatjetPt_pT10_); 
    babyTree_->Branch("FatjetEta_pT10", &FatjetEta_pT10_);
    babyTree_->Branch("FatjetPhi_pT10", &FatjetPhi_pT10_);
    babyTree_->Branch("FatjetPt_pT20",  &FatjetPt_pT20_); 
    babyTree_->Branch("FatjetEta_pT20", &FatjetEta_pT20_);
    babyTree_->Branch("FatjetPhi_pT20", &FatjetPhi_pT20_);
    babyTree_->Branch("FatjetPt_pT30",  &FatjetPt_pT30_); 
    babyTree_->Branch("FatjetEta_pT30", &FatjetEta_pT30_);
    babyTree_->Branch("FatjetPhi_pT30", &FatjetPhi_pT30_);
    babyTree_->Branch("RA4leptonPt",    &RA4leptonPt_);
    babyTree_->Branch("RA4leptonEta",   &RA4leptonEta_);
    babyTree_->Branch("RA4leptonPhi",   &RA4leptonPhi_);
    babyTree_->Branch("RA4leptonId",    &RA4leptonId_);
    
    // 
    // Event weights
    //
    // PU for MC
    TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_19p5ifb.root"); 
    TH1F *h1PU = (TH1F*)(fPUFile->Get("puWeights"));
    
    // Get cross section for MC 
    float Xsec=1; 
    if(!isData) Xsec=GetXsec(ProcessName);

    // json for DATA 
    std::vector< std::vector<int> > VRunLumi; VRunLumi.clear();
    if(isData) 
    {
        if(ProcessName.Contains("PromptReco")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("13Jul2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt"<< endl;       
        } else if(ProcessName.Contains("06Aug2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("24Aug2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("11Dec2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("22Jan2013")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" << endl;       
        } else 
        {
            cout << "[Error] No proper choice of JSON files!!" << endl;
            return ;
        }
    } else 
    {
        cout << "[MJ Analysis] No JSON files applied because it is MC" << endl;
    }

    //
    // main event loop
    //
    Int_t nentries = (Int_t)chainA->GetEntries();
    cout<<"[MJ Analysis] Number of entries is: "<<nentries<<endl;
    cout<<"[MJ Analysis] Number of entries of total sample is: "<<TotalNEntries<<endl;
    // Progress tracking 
    int i_permille_old = 0; 
    TDatime DTStart;
    int StartDate = DTStart.GetDate(); 
    int StartTime = DTStart.GetTime(); 
    cout << "[MJ Analysis] Start time : " << (StartTime/10000)%100 << ":"
         << (StartTime/100)%100 << ":" << StartTime%100
         << endl; 
    
    for(int i = 0; i<nentries; i++) 
    {
        // Progress indicator begin -------------------------------- 
        int i_permille = (int)floor(1000 * i / float(nentries));
        TDatime DTCurrent;
        int CurrentDate = DTCurrent.GetDate();
        int CurrentTime = DTCurrent.GetTime();
        int TimeLaps = (CurrentDate-StartDate)*1000000+(CurrentTime-StartTime);
        int TimeToRun = (int)((float)nentries/(float)i)*TimeLaps;
        if (i_permille != i_permille_old) 
        {
            // xterm magic from L. Vacavant and A. Cerri
            if (isatty(1)) 
            {
//                printf("\015\033[32m Processed :: \033[1m\033[31m%4.1f %%" 
//                       "\033[0m\033[32m   Expected processing time :: \033[1m\033[31m%i:%i:%i \033[0m\015",
//                        i_permille/10., (TimeToRun/10000)%100<60 ? (TimeToRun/10000)%100 : (TimeToRun/10000)%100-40, 
//                                        (TimeToRun/100)%100<60 ? (TimeToRun/100)%100 : (TimeToRun/100)%100-40, 
//                                        (TimeToRun%100)<60 ? (TimeToRun)%100 : (TimeToRun)%100-40 );
//                fflush(stdout);
            }
            i_permille_old = i_permille;
        } 
        // Progress indicator end ----------------------------------
        
        // Access to the event  
        chainA->GetEntry(i);
        chainB->GetEntry(i);
      
        // initialize baby variables
        run_        =   -1;
        lumiblock_  =   -1;
        event_      =   -1;
        Nfatjet_pT10_    =   -1;
        Nfatjet_pT20_    =   -1;
        Nfatjet_pT30_    =   -1;
        Nskinnyjet_ =   -1;
        Npv_        =   -1;
        Npu_        =   -1;
        EventWeight_=   1.;
        mll_        =-999.;
        MJ_pT10_    =-999.;
        MJ_pT20_    =-999.;
        MJ_pT30_    =-999.;
        HT_         =-999.;
        mj_pT10_.clear();
        mj_pT20_.clear();
        mj_pT30_.clear();
        FatjetPt_pT10_.clear();
        FatjetEta_pT10_.clear();
        FatjetPhi_pT10_.clear();
        FatjetPt_pT20_.clear();
        FatjetEta_pT20_.clear();
        FatjetPhi_pT20_.clear();
        FatjetPt_pT30_.clear();
        FatjetEta_pT30_.clear();
        FatjetPhi_pT30_.clear();
        RA4leptonPt_.clear();
        RA4leptonEta_.clear();
        RA4leptonPhi_.clear();
        RA4leptonId_.clear();

        // Flag weather or not to recored this event 
        bool FillResult = false;

        //
        // Core analysis 
        //
        // Get event weight 
        float EventWeight = 1;
        if(!isData) 
        {
            EventWeight = Xsec/TotalNEntries*Lumi; // scale for 1 pb * Lumi
            // need more weights if needed 
        } else 
        {
            if(!inJSON(VRunLumi,run,lumiblock)) continue; // JSON
        }

        // Get good RA4 muons
        vector<int> RA4MuonVeto; RA4MuonVeto.clear();
        vector<int> RA4Muon = GetRA4Muon(RA4MuonVeto);
        // Get good RA4 electrons
        vector<int> RA4ElecVeto; RA4ElecVeto.clear();
        vector<int> RA4Elec = GetRA4Elec(RA4ElecVeto, "", 0, true);
        // Get good skinny jets, HT and B-tagged jets 
        double HT=-999.; 
        vector<int> LooseBJet; 
        vector<int> MediumBJet; 
        vector<int> GoodJets_AK5PFclean = GetJets(RA4Muon,RA4Elec,RA4MuonVeto,RA4ElecVeto,
                                                  HT,LooseBJet,MediumBJet, 
                                                  2.4, 30, 0.3); 
        // 
        // Get mj and MJ 
        // 

        //
        // pT(R=0.5) > 10 GeV
        //
        double MJ_pT10=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjet_pT10; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PFclean_R1p2_R0p5pT10_px->size(); ifatjet++) 
        {
            float temp_pT_pT10 = TMath::Sqrt(fastjets_AK5PFclean_R1p2_R0p5pT10_px->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT10_px->at(ifatjet)
                                        +fastjets_AK5PFclean_R1p2_R0p5pT10_py->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT10_py->at(ifatjet));
            if(temp_pT_pT10<50) continue;
            TLorentzVector temp_GoodFatjet_pT10( fastjets_AK5PFclean_R1p2_R0p5pT10_px->at(ifatjet), 
                                            fastjets_AK5PFclean_R1p2_R0p5pT10_py->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT10_pz->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT10_energy->at(ifatjet));
            Vector_GoodFatjet_pT10.push_back(temp_GoodFatjet_pT10);
        } 
       
        vector<float> Vector_mj_pT10;   // mj
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT10.size(); igoodfatjet++) 
        {
            float temp_mj_pT10 = Getmj(Vector_GoodFatjet_pT10.at(igoodfatjet).Px(), 
                                  Vector_GoodFatjet_pT10.at(igoodfatjet).Py(),
                                  Vector_GoodFatjet_pT10.at(igoodfatjet).Pz(),
                                  Vector_GoodFatjet_pT10.at(igoodfatjet).E());
            Vector_mj_pT10.push_back(temp_mj_pT10);
        }
        MJ_pT10 = GetMJ(Vector_mj_pT10);

        int Nfatjet_pT10 = Vector_GoodFatjet_pT10.size(); 
        
        //
        // pT(R=0.5) > 20 GeV
        //
        double MJ_pT20=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjet_pT20; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PFclean_R1p2_R0p5pT20_px->size(); ifatjet++) 
        {
            float temp_pT_pT20 = TMath::Sqrt(fastjets_AK5PFclean_R1p2_R0p5pT20_px->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT20_px->at(ifatjet)
                                        +fastjets_AK5PFclean_R1p2_R0p5pT20_py->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT20_py->at(ifatjet));
            if(temp_pT_pT20<50) continue;
            TLorentzVector temp_GoodFatjet_pT20( fastjets_AK5PFclean_R1p2_R0p5pT20_px->at(ifatjet), 
                                            fastjets_AK5PFclean_R1p2_R0p5pT20_py->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT20_pz->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT20_energy->at(ifatjet));
            Vector_GoodFatjet_pT20.push_back(temp_GoodFatjet_pT20);
        } 
       
        vector<float> Vector_mj_pT20;   // mj
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT20.size(); igoodfatjet++) 
        {
            float temp_mj_pT20 = Getmj(Vector_GoodFatjet_pT20.at(igoodfatjet).Px(), 
                                  Vector_GoodFatjet_pT20.at(igoodfatjet).Py(),
                                  Vector_GoodFatjet_pT20.at(igoodfatjet).Pz(),
                                  Vector_GoodFatjet_pT20.at(igoodfatjet).E());
            Vector_mj_pT20.push_back(temp_mj_pT20);
        }
        MJ_pT20 = GetMJ(Vector_mj_pT20);

        int Nfatjet_pT20 = Vector_GoodFatjet_pT20.size(); 
        
        // 
        // pT(R=0.5) > 30 GeV
        // 
        double MJ_pT30=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjet_pT30; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PFclean_R1p2_R0p5pT30_px->size(); ifatjet++) 
        {
            float temp_pT_pT30 = TMath::Sqrt(fastjets_AK5PFclean_R1p2_R0p5pT30_px->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT30_px->at(ifatjet)
                                        +fastjets_AK5PFclean_R1p2_R0p5pT30_py->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT30_py->at(ifatjet));
            if(temp_pT_pT30<50) continue;
            TLorentzVector temp_GoodFatjet_pT30( fastjets_AK5PFclean_R1p2_R0p5pT30_px->at(ifatjet), 
                                            fastjets_AK5PFclean_R1p2_R0p5pT30_py->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT30_pz->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT30_energy->at(ifatjet));
            Vector_GoodFatjet_pT30.push_back(temp_GoodFatjet_pT30);
        } 
       
        vector<float> Vector_mj_pT30;   // mj
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT30.size(); igoodfatjet++) 
        {
            float temp_mj_pT30 = Getmj(Vector_GoodFatjet_pT30.at(igoodfatjet).Px(), 
                                  Vector_GoodFatjet_pT30.at(igoodfatjet).Py(),
                                  Vector_GoodFatjet_pT30.at(igoodfatjet).Pz(),
                                  Vector_GoodFatjet_pT30.at(igoodfatjet).E());
            Vector_mj_pT30.push_back(temp_mj_pT30);
        }
        MJ_pT30 = GetMJ(Vector_mj_pT30);

        int Nfatjet_pT30 = Vector_GoodFatjet_pT30.size(); 
        
        
        //
        // Select DYmm events 
        //
        // pT>20/20 GeV 
        // only 2 good muons 
        // MET < 20 GeV : not applied b/c it can kill high jet mulipicity events 
        float mll=-999.;
        int NpT20Muons = 0;
        int imus1=-1, imus2=-1;
        for(unsigned int imus=0; imus<RA4Muon.size(); imus++) 
        {
            if(mus_pt->at(RA4Muon.at(imus))>20) NpT20Muons++; 
            if(NpT20Muons==1) imus1 = RA4Muon.at(imus);
            if(NpT20Muons==2) imus2 = RA4Muon.at(imus);
        }
        if( NpT20Muons==2 &&                                // exactly two muons
            mus_charge->at(imus1)*mus_charge->at(imus2)<0)   // opposite sign 
        {
            TLorentzVector mus1P4(mus_px->at(imus1), mus_py->at(imus1), mus_pz->at(imus1), mus_energy->at(imus1));
            TLorentzVector mus2P4(mus_px->at(imus2), mus_py->at(imus2), mus_pz->at(imus2), mus_energy->at(imus2));
            mll = (mus1P4+mus2P4).M();  
        }
        
        //
        // Fill the baby variables 
        //
        run_            =   run;
        lumiblock_      =   lumiblock;
        event_          =   event;
        Nfatjet_pT10_   =   Nfatjet_pT10;
        Nfatjet_pT20_   =   Nfatjet_pT20;
        Nfatjet_pT30_   =   Nfatjet_pT30;
        Nskinnyjet_     =   GoodJets_AK5PFclean.size();
        Npv_            =   Npv;
        if(!isData) Npu_=   PU_TrueNumInteractions->at(1);
        EventWeight_    =   EventWeight;
        mll_            =   mll;
        MJ_pT10_        =   MJ_pT10;
        MJ_pT20_        =   MJ_pT20;
        MJ_pT30_        =   MJ_pT30;
        HT_             =   HT;
        
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT10.size(); igoodfatjet++) 
        { 
            mj_pT10_.push_back(Vector_mj_pT10.at(igoodfatjet));
            FatjetPt_pT10_.push_back(Vector_GoodFatjet_pT10.at(igoodfatjet).Pt());
            FatjetEta_pT10_.push_back(Vector_GoodFatjet_pT10.at(igoodfatjet).Eta());
            FatjetPhi_pT10_.push_back(Vector_GoodFatjet_pT10.at(igoodfatjet).Phi());
        }
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT20.size(); igoodfatjet++) 
        { 
            mj_pT20_.push_back(Vector_mj_pT20.at(igoodfatjet));
            FatjetPt_pT20_.push_back(Vector_GoodFatjet_pT20.at(igoodfatjet).Pt());
            FatjetEta_pT20_.push_back(Vector_GoodFatjet_pT20.at(igoodfatjet).Eta());
            FatjetPhi_pT20_.push_back(Vector_GoodFatjet_pT20.at(igoodfatjet).Phi());
        }
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT30.size(); igoodfatjet++) 
        { 
            mj_pT30_.push_back(Vector_mj_pT30.at(igoodfatjet));
            FatjetPt_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Pt());
            FatjetEta_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Eta());
            FatjetPhi_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Phi());
        }

        for(unsigned int imus=0; imus<RA4Muon.size(); imus++) 
        {
            RA4leptonPt_.push_back(mus_pt->at(RA4Muon.at(imus)));
            RA4leptonEta_.push_back(mus_eta->at(RA4Muon.at(imus)));
            RA4leptonPhi_.push_back(mus_phi->at(RA4Muon.at(imus)));
            RA4leptonId_.push_back(-1*13*mus_charge->at(RA4Muon.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4Elec.size(); iels++) 
        {
            RA4leptonPt_.push_back(els_pt->at(RA4Elec.at(iels)));
            RA4leptonEta_.push_back(els_eta->at(RA4Elec.at(iels)));
            RA4leptonPhi_.push_back(els_phi->at(RA4Elec.at(iels)));
            RA4leptonId_.push_back(-1*11*els_charge->at(RA4Elec.at(iels)));
        }
        
        babyTree_->Fill(); // Fill all events

        //for(int i=0; i<GoodJets_AK5PFclean.size(); i++) cout << event << " :: " << HT << endl;
        
        // Clean fat jets? 
        // (1) identify if all skinny jets associated with a given fat jet 
        //     are good jets by comparing index of skinny jet in GoodJets_AK5PFclean
        //     and fastjet_.._R1p2pTxx_index
        // (2) if a skinny jet is not in the good jet, its four vector is subtracted 
        //     from the fat jet and mj is calculated again
        // (3) MJ is then calculated


    } // event loop
    cout << endl;
    cout << "[MJ Analysis] Looping over events has been done" << endl;
    
    // clean up  
    fPUFile->Close();

    //
    // Write the baby file 
    //
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();

    // 
    //  
    // 
    TDatime DTEnd;
    int EndTime = DTEnd.GetTime(); 
    cout << "[MJ Analysis] End time : " << (EndTime/10000)%100 << ":"
         << (EndTime/100)%100 << ":" << EndTime%100
         << endl; 
    cout << "[MJ Analysis] Done with " << ProcessName << endl; 
    cout << "[MJ Analysis] ------------------------------------------------------------------------------------"<<endl; 
   
    // cleanup
    delete chainA;
    delete chainB;
}
