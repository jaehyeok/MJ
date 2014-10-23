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

#include "Branch_Fatjet_v1.h"
#include "ObjectSelector_Sync.h"
#include "EventSelector.h"
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
void DoOneSignal(TString InputName, TString ProcessName, int mg, int mlsp) 
{
    //
    cout << "[MJ Analysis] ------------------------------------------------------------------------------------"<<endl; 
    cout << "[MJ Analysis] Processing : " << ProcessName << endl;
    cout << "[MJ Analysis] Input dir  : " << InputName << endl;

    TFile *babyFile_ = new TFile(Form("baby_%s_f%i_%i.root", ProcessName.Data(), mg, mlsp), "RECREATE");
    babyFile_->cd();
    TTree *babyTree_ = new TTree("tree", "A Baby Ntuple");
    
    // 
    // Get tree 
    // 
    TChain * chainA = new TChain("/configurableAnalysis/eventA");   
    TChain * chainB = new TChain("/configurableAnalysis/eventB");  
    
    gSystem->Exec(Form("ls %s/*_f%i_%i_*.root", InputName.Data(), mg, mlsp));
    chainA->Add(Form("%s/*_f%i_%i_*.root", InputName.Data(), mg, mlsp));
    chainB->Add(Form("%s/*_f%i_%i_*.root", InputName.Data(), mg, mlsp));
    
    TList *l = (TList*)chainA->GetListOfFiles();
    l->Print();

    InitializeA(chainA);
    InitializeB(chainB);
    
    // Get total number of events of a sample 
    int TotalNEntries=1;
    
    //
    // Baby variables 
    //
    int run_; 
    int lumiblock_; 
    int event_; 
    vector<bool> filter_;  
    bool TrigMuon_;  
    bool TrigSingleMuon_;  
    bool TrigHTMuon_;  
    int Nfatjet_pT30_; 
//    int NfatjetCHS_pT30_; 
    int Nskinnyjet_; 
    int NBtagCSVM_; 
//    int NskinnyjetCHS_; 
//    int NBtagCHSCSVM_; 
    int Npv_; 
    float Npu_; 
    float EventWeight_; 
    float MJ_pT30_; 
//    float MJCHS_pT30_; 
    float MET_;  
    float METPhi_;  
    float HT_;  
//    float HTCHS_;  
    vector<float> mj_pT30_;
//    vector<float> mjCHS_pT30_;
    vector<float> FatjetPt_pT30_;
    vector<float> FatjetEta_pT30_;
    vector<float> FatjetPhi_pT30_;
//    vector<float> FatjetCHSPt_pT30_;
//    vector<float> FatjetCHSEta_pT30_;
//    vector<float> FatjetCHSPhi_pT30_;
    vector<float> RA4ElsPt_;
    vector<float> RA4ElsEta_;
    vector<float> RA4ElsPhi_;
    vector<float> RA4MusPt_;
    vector<float> RA4MusEta_;
    vector<float> RA4MusPhi_;
    vector<float> RA4ElsVetoPt_;
    vector<float> RA4ElsVetoEta_;
    vector<float> RA4ElsVetoPhi_;
    vector<float> RA4MusVetoPt_;
    vector<float> RA4MusVetoEta_;
    vector<float> RA4MusVetoPhi_;
    vector<float> JetPt_;
    vector<float> JetEta_;
    vector<float> JetPhi_;
    vector<float> JetCSV_;
    vector<float> GenPt_;
    vector<float> GenEta_;
    vector<float> GenPhi_;
    vector<float> GenId_;
    vector<float> GenMId_;
    vector<float> GenGMId_;
    
    babyTree_->Branch("run",            	&run_);    
    babyTree_->Branch("lumiblock",      	&lumiblock_); 
    babyTree_->Branch("event",          	&event_);     
    babyTree_->Branch("filter",          	&filter_);     
    babyTree_->Branch("TrigMuon",          	&TrigMuon_);     
    babyTree_->Branch("TrigSingleMuon",    	&TrigMuon_);     
    babyTree_->Branch("TrigHTMuon",       	&TrigHTMuon_);     
    babyTree_->Branch("Nfatjet_pT30",   	&Nfatjet_pT30_);   
//    babyTree_->Branch("NatjetCHS_pT30", 	&NfatjetCHS_pT30_);   
    babyTree_->Branch("Nskinnyjet",     	&Nskinnyjet_);
    babyTree_->Branch("NBtagCSVM",     	    &NBtagCSVM_);
//    babyTree_->Branch("NskinnyjetCHS",  	&NskinnyjetCHS_);
//    babyTree_->Branch("NBtagCHSCSVM",     	&NBtagCHSCSVM_);
    babyTree_->Branch("Npv",            	&Npv_);       
    babyTree_->Branch("Npu",            	&Npu_);       
    babyTree_->Branch("EventWeight",    	&EventWeight_);
    babyTree_->Branch("MJ_pT30",        	&MJ_pT30_);        
//    babyTree_->Branch("MJCHS_pT30",        	&MJCHS_pT30_);        
    babyTree_->Branch("MET",            	&MET_);        
    babyTree_->Branch("METPhi",            	&METPhi_);        
    babyTree_->Branch("HT",             	&HT_);        
    babyTree_->Branch("mj_pT30",        	&mj_pT30_);     
//    babyTree_->Branch("mjCHS_pT30",        	&mjCHS_pT30_);     
    babyTree_->Branch("FatjetPt_pT30",  	&FatjetPt_pT30_); 
    babyTree_->Branch("FatjetEta_pT30", 	&FatjetEta_pT30_);
    babyTree_->Branch("FatjetPhi_pT30",     &FatjetPhi_pT30_);
//    babyTree_->Branch("FatjetCHSPt_pT30",   &FatjetCHSPt_pT30_); 
//    babyTree_->Branch("FatjetCHSEta_pT30",  &FatjetCHSEta_pT30_);
//    babyTree_->Branch("FatjetCHSPhi_pT30",  &FatjetCHSPhi_pT30_);
    babyTree_->Branch("RA4ElsPt",           &RA4ElsPt_);
    babyTree_->Branch("RA4ElsEta",          &RA4ElsEta_);
    babyTree_->Branch("RA4ElsPhi",          &RA4ElsPhi_);
    babyTree_->Branch("RA4MusPt",           &RA4MusPt_);
    babyTree_->Branch("RA4MusEta",          &RA4MusEta_);
    babyTree_->Branch("RA4MusPhi",          &RA4MusPhi_);
    babyTree_->Branch("RA4ElsVetoPt",       &RA4ElsVetoPt_);
    babyTree_->Branch("RA4ElsVetoEta",      &RA4ElsVetoEta_);
    babyTree_->Branch("RA4ElsVetoPhi",      &RA4ElsVetoPhi_);
    babyTree_->Branch("RA4MusVetoPt",       &RA4MusVetoPt_);
    babyTree_->Branch("RA4MusVetoEta",      &RA4MusVetoEta_);
    babyTree_->Branch("RA4MusVetoPhi",      &RA4MusVetoPhi_);
    babyTree_->Branch("JetPt",              &JetPt_);
    babyTree_->Branch("JetEta",             &JetEta_);
    babyTree_->Branch("JetPhi",             &JetPhi_);
    babyTree_->Branch("JetCSV",             &JetCSV_);
    babyTree_->Branch("GenPt",              &GenPt_);
    babyTree_->Branch("GenEta",             &GenEta_);
    babyTree_->Branch("GenPhi",             &GenPhi_);
    babyTree_->Branch("GenId",              &GenId_);
    babyTree_->Branch("GenMId",             &GenMId_);
    babyTree_->Branch("GenGMId",            &GenGMId_);
   
    // 
    // Event weights
    //
    // PU for MC
    TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_19p5ifb.root"); 
    TH1F *h1PU = (TH1F*)(fPUFile->Get("puWeights"));

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
        if(i%10000==0) cout << "Progress : " << i << "/" << nentries << endl;
        // Progress indicator end ----------------------------------
        
        // Access to the event  
        chainA->GetEntry(i);
        chainB->GetEntry(i);
      
        // initialize baby variables
        run_                =   -1;
        lumiblock_          =   -1;
        event_              =   -1;
        TrigMuon_           =   1;
        TrigSingleMuon_     =   1;
        TrigHTMuon_         =   1;
        Nfatjet_pT30_       =   -1;
//        NfatjetCHS_pT30_    =   -1;
        Nskinnyjet_         =   -1;
        NBtagCSVM_          =   -1;
//        NskinnyjetCHS_      =   -1;
//        NBtagCHSCSVM_       =   -1;
        Npv_                =   -1;
        Npu_                =   -1;
        EventWeight_        =   1.;
        MJ_pT30_            =-999.;
//        MJCHS_pT30_         =-999.;
        MET_                =-999.;
        METPhi_             =-999.;
        HT_                 =-999.;
//        HTCHS_              =-999.;
        filter_.clear(); 
        mj_pT30_.clear();
//        mjCHS_pT30_.clear();
        FatjetPt_pT30_.clear();
        FatjetEta_pT30_.clear();
        FatjetPhi_pT30_.clear();
//        FatjetCHSPt_pT30_.clear();
//        FatjetCHSEta_pT30_.clear();
//        FatjetCHSPhi_pT30_.clear();
        RA4ElsPt_.clear();
        RA4ElsEta_.clear();
        RA4ElsPhi_.clear();
        RA4MusPt_.clear();
        RA4MusEta_.clear();
        RA4MusPhi_.clear();
        RA4ElsVetoPt_.clear();
        RA4ElsVetoEta_.clear();
        RA4ElsVetoPhi_.clear();
        RA4MusVetoPt_.clear();
        RA4MusVetoEta_.clear();
        RA4MusVetoPhi_.clear();
        JetPt_.clear();
        JetEta_.clear();
        JetPhi_.clear();
        JetCSV_.clear();
        GenPt_.clear();
        GenEta_.clear();
        GenPhi_.clear();
        GenId_.clear();
        GenMId_.clear();
        GenGMId_.clear();

        //
        // Core analysis 
        //
        // Get event weight 
        float EventWeight = 1;

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
//        double HTCHS=-999.; 
//        vector<int> LooseBJetCHS; 
//        vector<int> MediumBJetCHS; 
//        vector<int> GoodJets_AK5PF = GetJetsCHS(RA4Muon,RA4Elec,RA4MuonVeto,RA4ElecVeto,
//                                                HTCHS,LooseBJetCHS,MediumBJetCHS, 
//                                                2.4, 30, 0.3); 
        // 
        // pT(R=0.5) > 30 GeV
        // 
        
        // Regular jets
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
       
        // CHS jets
//        double MJCHS_pT30=-999.; 
//        // first, select good fat jets 
//        vector<TLorentzVector> Vector_GoodFatjetCHS_pT30; 
//        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PF_R1p2_R0p5pT30_px->size(); ifatjet++) 
//        {
//            float temp_pT_pT30 = TMath::Sqrt(fastjets_AK5PF_R1p2_R0p5pT30_px->at(ifatjet)*fastjets_AK5PF_R1p2_R0p5pT30_px->at(ifatjet)
//                                        +fastjets_AK5PF_R1p2_R0p5pT30_py->at(ifatjet)*fastjets_AK5PF_R1p2_R0p5pT30_py->at(ifatjet));
//            if(temp_pT_pT30<50) continue;
//            TLorentzVector temp_GoodFatjetCHS_pT30( fastjets_AK5PF_R1p2_R0p5pT30_px->at(ifatjet), 
//                                            fastjets_AK5PF_R1p2_R0p5pT30_py->at(ifatjet),
//                                            fastjets_AK5PF_R1p2_R0p5pT30_pz->at(ifatjet),
//                                            fastjets_AK5PF_R1p2_R0p5pT30_energy->at(ifatjet));
//            Vector_GoodFatjetCHS_pT30.push_back(temp_GoodFatjetCHS_pT30);
//        } 
//       
//        vector<float> Vector_mjCHS_pT30;   // mj
//        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjetCHS_pT30.size(); igoodfatjet++) 
//        {
//            float temp_mj_pT30 = Getmj(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Px(), 
//                                       Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Py(),
//                                       Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Pz(),
//                                       Vector_GoodFatjetCHS_pT30.at(igoodfatjet).E());
//            Vector_mjCHS_pT30.push_back(temp_mj_pT30);
//        }
//        MJCHS_pT30 = GetMJ(Vector_mjCHS_pT30);
//
//        int NfatjetCHS_pT30 = Vector_GoodFatjetCHS_pT30.size(); 
//        

        //
        // Fill the baby variables 
        //
        run_                =   run;
        lumiblock_          =   lumiblock;
        event_              =   event; 
        TrigMuon_           =   PassMuonTrig(); 
        TrigHTMuon_         =   PassHTMuonTrig(); 
        TrigSingleMuon_     =   PassSingleMuonTrig(); 
        Nfatjet_pT30_       =   Nfatjet_pT30;
//        NfatjetCHS_pT30_    =   NfatjetCHS_pT30;
        Nskinnyjet_         =   GoodJets_AK5PFclean.size();
        NBtagCSVM_          =   MediumBJet.size();
//        NskinnyjetCHS_      =   GoodJets_AK5PF.size();
//        NBtagCHSCSVM_       =   MediumBJetCHS.size();
        Npv_                =   Npv;
        Npu_                =   PU_TrueNumInteractions->at(0);
        EventWeight_        =   EventWeight;
        MJ_pT30_            =   MJ_pT30;
//        MJCHS_pT30_         =   MJCHS_pT30;
        MET_                =   pfTypeImets_et->at(0);
        METPhi_             =   pfTypeImets_phi->at(0);
        HT_                 =   HT;
//        HTCHS_              =   HTCHS;
        
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT30.size(); igoodfatjet++) 
        { 
            mj_pT30_.push_back(Vector_mj_pT30.at(igoodfatjet));
            FatjetPt_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Pt());
            FatjetEta_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Eta());
            FatjetPhi_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Phi());
        }

//        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjetCHS_pT30.size(); igoodfatjet++) 
//        { 
//            mjCHS_pT30_.push_back(Vector_mjCHS_pT30.at(igoodfatjet));
//            FatjetCHSPt_pT30_.push_back(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Pt());
//            FatjetCHSEta_pT30_.push_back(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Eta());
//            FatjetCHSPhi_pT30_.push_back(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Phi());
//        }
        
        for(unsigned int imus=0; imus<RA4Muon.size(); imus++) 
        {
            RA4MusPt_.push_back(mus_pt->at(RA4Muon.at(imus)));
            RA4MusEta_.push_back(mus_eta->at(RA4Muon.at(imus)));
            RA4MusPhi_.push_back(mus_phi->at(RA4Muon.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4Elec.size(); iels++) 
        {
            RA4ElsPt_.push_back(els_pt->at(RA4Elec.at(iels)));
            RA4ElsEta_.push_back(els_eta->at(RA4Elec.at(iels)));
            RA4ElsPhi_.push_back(els_phi->at(RA4Elec.at(iels)));
        }
        for(unsigned int imus=0; imus<RA4MuonVeto.size(); imus++) 
        {
            RA4MusVetoPt_.push_back(mus_pt->at(RA4MuonVeto.at(imus)));
            RA4MusVetoEta_.push_back(mus_eta->at(RA4MuonVeto.at(imus)));
            RA4MusVetoPhi_.push_back(mus_phi->at(RA4MuonVeto.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4ElecVeto.size(); iels++) 
        {
            RA4ElsVetoPt_.push_back(els_pt->at(RA4ElecVeto.at(iels)));
            RA4ElsVetoEta_.push_back(els_eta->at(RA4ElecVeto.at(iels)));
            RA4ElsVetoPhi_.push_back(els_phi->at(RA4ElecVeto.at(iels)));
        }
        
        for(int igoodjet=0; igoodjet<(int)GoodJets_AK5PFclean.size(); igoodjet++) 
        {
          int ijet = GoodJets_AK5PFclean.at(igoodjet); 
          JetPt_.push_back(jets_AK5PFclean_pt->at(ijet)); 
          JetEta_.push_back(jets_AK5PFclean_eta->at(ijet)); 
          JetPhi_.push_back(jets_AK5PFclean_phi->at(ijet)); 
          JetCSV_.push_back(jets_AK5PFclean_btag_secVertexCombined->at(ijet)); 
        }
        for(int igen=0; igen<mc_doc_id->size(); igen++) 
        { 
            GenPt_.push_back(   mc_doc_pt->at(igen));  
            GenPhi_.push_back(  mc_doc_phi->at(igen));  
            GenEta_.push_back(  mc_doc_eta->at(igen));  
            GenId_.push_back(   mc_doc_id->at(igen));  
            GenMId_.push_back(  mc_doc_mother_id->at(igen));  
            GenGMId_.push_back( mc_doc_grandmother_id->at(igen));  
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
