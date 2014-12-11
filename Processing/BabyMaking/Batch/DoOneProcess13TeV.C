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

#include "Branch_v75.h" 
#include "ra4_objects_13TeV.h"
//#include "Branch_v77.h" 
//#include "r/a4_objects_13TeV_new.h"
#include "Utilities_13TeV.h"

// include necessary fastjet files
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

//
// Make fat jets 
//
vector<TLorentzVector> makeFatJet( vector<TLorentzVector> FatJetConstituent, 
                                   double Rparam=1.2, 
                                   int ConstituentpTcut=30, 
                                   float ConstituentEtacut=100) 
{
    vector<TLorentzVector> FatJets;

    // Loop over R=0.5 jets, form into PseudoJets vector
    vector<fastjet::PseudoJet> input_particles;
    double FatjetConstituent_px_tmp, FatjetConstituent_py_tmp, FatjetConstituent_pz_tmp, FatjetConstituent_energy_tmp;

    for(int ijet = 0; ijet<(int)FatJetConstituent.size(); ijet++) 
    { 

        FatjetConstituent_px_tmp        = FatJetConstituent.at(ijet).Px();
        FatjetConstituent_py_tmp        = FatJetConstituent.at(ijet).Py();
        FatjetConstituent_pz_tmp        = FatJetConstituent.at(ijet).Pz();
        FatjetConstituent_energy_tmp    = FatJetConstituent.at(ijet).E();	  
    
//        cout << ijet << " :: " 
//             << FatjetConstituent_px_tmp << " " 
//             << FatjetConstituent_py_tmp << " " 
//             << FatjetConstituent_pz_tmp << " " 
//             << FatjetConstituent_energy_tmp << " " 
//             << endl;

        if(TMath::Sqrt( FatjetConstituent_px_tmp*FatjetConstituent_px_tmp
                       +FatjetConstituent_py_tmp*FatjetConstituent_py_tmp)<ConstituentpTcut) continue;

        if(TMath::Abs(FatJetConstituent.at(ijet).Eta())>ConstituentEtacut) continue;

        input_particles.push_back(fastjet::PseudoJet( FatjetConstituent_px_tmp, FatjetConstituent_py_tmp,
                                                      FatjetConstituent_pz_tmp, FatjetConstituent_energy_tmp));
    }
    
    //
    // Run Fastjet to reconstuct jets 
    //

    // Create an object that represents your choice of jet algorithm and the associated parameters
    fastjet::Strategy strategy = fastjet::Best;
    fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, recomb_scheme, strategy);

    // run the jet clustering with the above jet definition
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);

    // 
    // Get p4 of the reconstructed jets  
    //
    double ptmin = 0.0; // could use 3.0 here, instead of applying later
    vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
    //Sort by pt
    vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(inclusive_jets);
    //fill fastjet output into vectors, continue as original code
    for(int isortjets = 0; isortjets< (int)sorted_jets.size(); isortjets++)
    {
        //store only if pt >3 GeV to match CMS jets
        if(TMath::Sqrt( sorted_jets[isortjets].px()*sorted_jets[isortjets].px()
                       +sorted_jets[isortjets].py()*sorted_jets[isortjets].py())>50) 
        {
            TLorentzVector FatJet_tmp( sorted_jets[isortjets].px(), sorted_jets[isortjets].py(), 
                                       sorted_jets[isortjets].pz(), sorted_jets[isortjets].E());
            FatJets.push_back(FatJet_tmp);

//            cout << isortjets << " :: "  
//                 << sorted_jets[isortjets].px() << " " 
//                 << sorted_jets[isortjets].py() << " " 
//                 << sorted_jets[isortjets].pz() << " " 
//                 << sorted_jets[isortjets].E() << " " 
//                 << endl;
        }
    }

    return FatJets;
}

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
void DoOneProcess13TeV(TString InputName, TString ProcessName, int ibegin, int iend, bool isData, float Lumi) 
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
    
    for(int i=ibegin; i<=iend; i++)  
    {
        gSystem->Exec(Form("ls %s/*_f%i.root", InputName.Data(), i));
        chainA->Add(Form("%s/*_f%i.root", InputName.Data(), i));
        chainB->Add(Form("%s/*_f%i.root", InputName.Data(), i));
        gSystem->Exec(Form("ls %s/*_f%i_*.root", InputName.Data(), i));
        chainA->Add(Form("%s/*_f%i_*.root", InputName.Data(), i));
        chainB->Add(Form("%s/*_f%i_*.root", InputName.Data(), i));
        gSystem->Exec(Form("ls %s/configurableAnalysis_%i.root", InputName.Data(), i)); // for non-published cfA samples 
        chainA->Add(Form("%s/configurableAnalysis_%i.root", InputName.Data(), i));
        chainB->Add(Form("%s/configurableAnalysis_%i.root", InputName.Data(), i));
    } 
    TList *l = (TList*)chainA->GetListOfFiles();
    l->Print();

    InitializeA(chainA);
    InitializeB(chainB);
    
    // Get total number of events of a sample 
    int TotalNEntries=1;
    
    if(!isData)
    { 
        TChain * chainATotal = new TChain("/configurableAnalysis/eventA");  
        chainATotal->Add(Form("%s/*.root", InputName.Data()));
        TotalNEntries = (int)chainATotal->GetEntries();
    }

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
    int Nskinnyjet_; 
    int NBtagCSVM_; 
    int NBtagCSVL_; 
    int Npv_; 
    float Npu_; 
    float EventWeight_; 
    float MJ_pT30_; 
    float MET_;  
    float METPhi_;  
    float HT_;  
    float top1pT_;  
    float top1Eta_;  
    float top1Phi_;  
    float top2pT_;  
    float top2Eta_;  
    float top2Phi_;  
    vector<float> mj_pT30_;
    vector<float> FatjetPt_pT30_;
    vector<float> FatjetEta_pT30_;
    vector<float> FatjetPhi_pT30_;
    vector<float> FatjetN_pT30_;
    vector<float> mj_R0p8_pT30_Eta5_;
    vector<float> FatjetPt_R0p8_pT30_Eta5_;
    vector<float> FatjetEta_R0p8_pT30_Eta5_;
    vector<float> FatjetPhi_R0p8_pT30_Eta5_;
    vector<float> mj_R0p9_pT30_Eta5_;
    vector<float> FatjetPt_R0p9_pT30_Eta5_;
    vector<float> FatjetEta_R0p9_pT30_Eta5_;
    vector<float> FatjetPhi_R0p9_pT30_Eta5_;
    vector<float> mj_R1p0_pT30_Eta5_;
    vector<float> FatjetPt_R1p0_pT30_Eta5_;
    vector<float> FatjetEta_R1p0_pT30_Eta5_;
    vector<float> FatjetPhi_R1p0_pT30_Eta5_;
    vector<float> mj_R1p1_pT30_Eta5_;
    vector<float> FatjetPt_R1p1_pT30_Eta5_;
    vector<float> FatjetEta_R1p1_pT30_Eta5_;
    vector<float> FatjetPhi_R1p1_pT30_Eta5_;
    vector<float> mj_R1p2_pT30_Eta5_;
    vector<float> FatjetPt_R1p2_pT30_Eta5_;
    vector<float> FatjetEta_R1p2_pT30_Eta5_;
    vector<float> FatjetPhi_R1p2_pT30_Eta5_;
    vector<vector<float> >   *Fatjet_sjmj_R1p2_pT30_Eta5_ = 0;
    vector<float> mj_R1p3_pT30_Eta5_;
    vector<float> FatjetPt_R1p3_pT30_Eta5_;
    vector<float> FatjetEta_R1p3_pT30_Eta5_;
    vector<float> FatjetPhi_R1p3_pT30_Eta5_;
    vector<float> mj_R1p4_pT30_Eta5_;
    vector<float> FatjetPt_R1p4_pT30_Eta5_;
    vector<float> FatjetEta_R1p4_pT30_Eta5_;
    vector<float> FatjetPhi_R1p4_pT30_Eta5_;
    vector<float> mj_R1p5_pT30_Eta5_;
    vector<float> FatjetPt_R1p5_pT30_Eta5_;
    vector<float> FatjetEta_R1p5_pT30_Eta5_;
    vector<float> FatjetPhi_R1p5_pT30_Eta5_;
    vector<float> mj_R0p8_pT30_Eta2p5_;
    vector<float> FatjetPt_R0p8_pT30_Eta2p5_;
    vector<float> FatjetEta_R0p8_pT30_Eta2p5_;
    vector<float> FatjetPhi_R0p8_pT30_Eta2p5_;
    vector<float> mj_R0p9_pT30_Eta2p5_;
    vector<float> FatjetPt_R0p9_pT30_Eta2p5_;
    vector<float> FatjetEta_R0p9_pT30_Eta2p5_;
    vector<float> FatjetPhi_R0p9_pT30_Eta2p5_;
    vector<float> mj_R1p0_pT30_Eta2p5_;
    vector<float> FatjetPt_R1p0_pT30_Eta2p5_;
    vector<float> FatjetEta_R1p0_pT30_Eta2p5_;
    vector<float> FatjetPhi_R1p0_pT30_Eta2p5_;
    vector<float> mj_R1p1_pT30_Eta2p5_;
    vector<float> FatjetPt_R1p1_pT30_Eta2p5_;
    vector<float> FatjetEta_R1p1_pT30_Eta2p5_;
    vector<float> FatjetPhi_R1p1_pT30_Eta2p5_;
    vector<float> mj_R1p2_pT30_Eta2p5_;
    vector<float> FatjetPt_R1p2_pT30_Eta2p5_;
    vector<float> FatjetEta_R1p2_pT30_Eta2p5_;
    vector<float> FatjetPhi_R1p2_pT30_Eta2p5_;
    vector<float> mj_R1p3_pT30_Eta2p5_;
    vector<float> FatjetPt_R1p3_pT30_Eta2p5_;
    vector<float> FatjetEta_R1p3_pT30_Eta2p5_;
    vector<float> FatjetPhi_R1p3_pT30_Eta2p5_;
    vector<float> mj_R1p4_pT30_Eta2p5_;
    vector<float> FatjetPt_R1p4_pT30_Eta2p5_;
    vector<float> FatjetEta_R1p4_pT30_Eta2p5_;
    vector<float> FatjetPhi_R1p4_pT30_Eta2p5_;
    vector<float> mj_R1p5_pT30_Eta2p5_;
    vector<float> FatjetPt_R1p5_pT30_Eta2p5_;
    vector<float> FatjetEta_R1p5_pT30_Eta2p5_;
    vector<float> FatjetPhi_R1p5_pT30_Eta2p5_;
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
    babyTree_->Branch("Nskinnyjet",     	&Nskinnyjet_);
    babyTree_->Branch("NBtagCSVM",     	    &NBtagCSVM_);
    babyTree_->Branch("NBtagCSVL",     	    &NBtagCSVL_);
    babyTree_->Branch("Npv",            	&Npv_);       
    babyTree_->Branch("Npu",            	&Npu_);       
    babyTree_->Branch("EventWeight",    	&EventWeight_);
    babyTree_->Branch("MJ_pT30",        	&MJ_pT30_);        
    babyTree_->Branch("MET",            	&MET_);        
    babyTree_->Branch("METPhi",            	&METPhi_);        
    babyTree_->Branch("HT",             	&HT_);        
    babyTree_->Branch("top1pT",          	&top1pT_);        
    babyTree_->Branch("top1Eta",          	&top1Eta_);        
    babyTree_->Branch("top1Phi",          	&top1Phi_);        
    babyTree_->Branch("top2pT",          	&top2pT_);        
    babyTree_->Branch("top2Eta",          	&top2Eta_);        
    babyTree_->Branch("top2Phi",          	&top2Phi_);        
    babyTree_->Branch("mj_pT30",        	&mj_pT30_);     
    babyTree_->Branch("FatjetPt_pT30",  	&FatjetPt_pT30_); 
    babyTree_->Branch("FatjetEta_pT30", 	&FatjetEta_pT30_);
    babyTree_->Branch("FatjetPhi_pT30",     &FatjetPhi_pT30_);
    babyTree_->Branch("FatjetN_pT30",       &FatjetN_pT30_);
    babyTree_->Branch("mj_R0p8_pT30_Eta5",        	&mj_R0p8_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R0p8_pT30_Eta5",  	&FatjetPt_R0p8_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R0p8_pT30_Eta5", 	&FatjetEta_R0p8_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R0p8_pT30_Eta5",     &FatjetPhi_R0p8_pT30_Eta5_);
    babyTree_->Branch("mj_R0p9_pT30_Eta5",        	&mj_R0p9_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R0p9_pT30_Eta5",  	&FatjetPt_R0p9_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R0p9_pT30_Eta5", 	&FatjetEta_R0p9_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R0p9_pT30_Eta5",     &FatjetPhi_R0p9_pT30_Eta5_);
    babyTree_->Branch("mj_R1p0_pT30_Eta5",        	&mj_R1p0_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p0_pT30_Eta5",  	&FatjetPt_R1p0_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p0_pT30_Eta5", 	&FatjetEta_R1p0_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p0_pT30_Eta5",     &FatjetPhi_R1p0_pT30_Eta5_);
    babyTree_->Branch("mj_R1p1_pT30_Eta5",        	&mj_R1p1_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p1_pT30_Eta5",  	&FatjetPt_R1p1_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p1_pT30_Eta5", 	&FatjetEta_R1p1_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p1_pT30_Eta5",     &FatjetPhi_R1p1_pT30_Eta5_);
    babyTree_->Branch("mj_R1p2_pT30_Eta5",        	&mj_R1p2_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p2_pT30_Eta5",  	&FatjetPt_R1p2_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p2_pT30_Eta5", 	&FatjetEta_R1p2_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p2_pT30_Eta5",     &FatjetPhi_R1p2_pT30_Eta5_);
    babyTree_->Branch("Fatjet_sjmj_R1p2_pT30_Eta5",     &Fatjet_sjmj_R1p2_pT30_Eta5_);
    babyTree_->Branch("mj_R1p3_pT30_Eta5",        	&mj_R1p3_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p3_pT30_Eta5",  	&FatjetPt_R1p3_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p3_pT30_Eta5", 	&FatjetEta_R1p3_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p3_pT30_Eta5",     &FatjetPhi_R1p3_pT30_Eta5_);
    babyTree_->Branch("mj_R1p4_pT30_Eta5",        	&mj_R1p4_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p4_pT30_Eta5",  	&FatjetPt_R1p4_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p4_pT30_Eta5", 	&FatjetEta_R1p4_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p4_pT30_Eta5",     &FatjetPhi_R1p4_pT30_Eta5_);
    babyTree_->Branch("mj_R1p5_pT30_Eta5",        	&mj_R1p5_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p5_pT30_Eta5",  	&FatjetPt_R1p5_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p5_pT30_Eta5", 	&FatjetEta_R1p5_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p5_pT30_Eta5",     &FatjetPhi_R1p5_pT30_Eta5_);
    babyTree_->Branch("mj_R0p8_pT30_Eta2p5",        	&mj_R0p8_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R0p8_pT30_Eta2p5",  	&FatjetPt_R0p8_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R0p8_pT30_Eta2p5", 	&FatjetEta_R0p8_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R0p8_pT30_Eta2p5",     &FatjetPhi_R0p8_pT30_Eta2p5_);
    babyTree_->Branch("mj_R0p9_pT30_Eta2p5",        	&mj_R0p9_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R0p9_pT30_Eta2p5",  	&FatjetPt_R0p9_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R0p9_pT30_Eta2p5", 	&FatjetEta_R0p9_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R0p9_pT30_Eta2p5",     &FatjetPhi_R0p9_pT30_Eta2p5_);
    babyTree_->Branch("mj_R1p0_pT30_Eta2p5",        	&mj_R1p0_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p0_pT30_Eta2p5",  	&FatjetPt_R1p0_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p0_pT30_Eta2p5", 	&FatjetEta_R1p0_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p0_pT30_Eta2p5",     &FatjetPhi_R1p0_pT30_Eta2p5_);
    babyTree_->Branch("mj_R1p1_pT30_Eta2p5",        	&mj_R1p1_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p1_pT30_Eta2p5",  	&FatjetPt_R1p1_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p1_pT30_Eta2p5", 	&FatjetEta_R1p1_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p1_pT30_Eta2p5",     &FatjetPhi_R1p1_pT30_Eta2p5_);
    babyTree_->Branch("mj_R1p2_pT30_Eta2p5",        	&mj_R1p2_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p2_pT30_Eta2p5",  	&FatjetPt_R1p2_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p2_pT30_Eta2p5", 	&FatjetEta_R1p2_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p2_pT30_Eta2p5",     &FatjetPhi_R1p2_pT30_Eta2p5_);
    babyTree_->Branch("mj_R1p3_pT30_Eta2p5",        	&mj_R1p3_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p3_pT30_Eta2p5",  	&FatjetPt_R1p3_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p3_pT30_Eta2p5", 	&FatjetEta_R1p3_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p3_pT30_Eta2p5",     &FatjetPhi_R1p3_pT30_Eta2p5_);
    babyTree_->Branch("mj_R1p4_pT30_Eta2p5",        	&mj_R1p4_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p4_pT30_Eta2p5",  	&FatjetPt_R1p4_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p4_pT30_Eta2p5", 	&FatjetEta_R1p4_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p4_pT30_Eta2p5",     &FatjetPhi_R1p4_pT30_Eta2p5_);
    babyTree_->Branch("mj_R1p5_pT30_Eta2p5",        	&mj_R1p5_pT30_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p5_pT30_Eta2p5",  	&FatjetPt_R1p5_pT30_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p5_pT30_Eta2p5", 	&FatjetEta_R1p5_pT30_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p5_pT30_Eta2p5",     &FatjetPhi_R1p5_pT30_Eta2p5_);
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
    
    // Get cross section for MC 
    float Xsec=1; 
    if(!isData) Xsec=GetXsec(ProcessName);

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
//        int TimeToRun = (int)((float)nentries/(float)i)*TimeLaps;
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
        Nskinnyjet_         =   -1;
        NBtagCSVM_          =   -1;
        NBtagCSVL_          =   -1;
        Npv_                =   -1;
        Npu_                =   -1;
        EventWeight_        =   1.;
        MJ_pT30_            =-999.;
        MET_                =-999.;
        METPhi_             =-999.;
        HT_                 =-999.;
        top1pT_             =-999.;
        top1Phi_            =-999.;
        top1Eta_            =-999.;
        top2pT_             =-999.;
        top2Phi_            =-999.;
        top2Eta_            =-999.;
        filter_.clear(); 
        mj_pT30_.clear();
        FatjetPt_pT30_.clear();
        FatjetEta_pT30_.clear();
        FatjetPhi_pT30_.clear();
        FatjetN_pT30_.clear();
        mj_R0p8_pT30_Eta5_.clear();
        FatjetPt_R0p8_pT30_Eta5_.clear();
        FatjetEta_R0p8_pT30_Eta5_.clear();
        FatjetPhi_R0p8_pT30_Eta5_.clear();
        mj_R0p9_pT30_Eta5_.clear();
        FatjetPt_R0p9_pT30_Eta5_.clear();
        FatjetEta_R0p9_pT30_Eta5_.clear();
        FatjetPhi_R0p9_pT30_Eta5_.clear();
        mj_R1p0_pT30_Eta5_.clear();
        FatjetPt_R1p0_pT30_Eta5_.clear();
        FatjetEta_R1p0_pT30_Eta5_.clear();
        FatjetPhi_R1p0_pT30_Eta5_.clear();
        mj_R1p1_pT30_Eta5_.clear();
        FatjetPt_R1p1_pT30_Eta5_.clear();
        FatjetEta_R1p1_pT30_Eta5_.clear();
        FatjetPhi_R1p1_pT30_Eta5_.clear();
        mj_R1p2_pT30_Eta5_.clear();
        FatjetPt_R1p2_pT30_Eta5_.clear();
        FatjetEta_R1p2_pT30_Eta5_.clear();
        FatjetPhi_R1p2_pT30_Eta5_.clear();
        Fatjet_sjmj_R1p2_pT30_Eta5_.clear(); 
        mj_R1p3_pT30_Eta5_.clear();
        FatjetPt_R1p3_pT30_Eta5_.clear();
        FatjetEta_R1p3_pT30_Eta5_.clear();
        FatjetPhi_R1p3_pT30_Eta5_.clear();
        mj_R1p4_pT30_Eta5_.clear();
        FatjetPt_R1p4_pT30_Eta5_.clear();
        FatjetEta_R1p4_pT30_Eta5_.clear();
        FatjetPhi_R1p4_pT30_Eta5_.clear();
        mj_R1p5_pT30_Eta5_.clear();
        FatjetPt_R1p5_pT30_Eta5_.clear();
        FatjetEta_R1p5_pT30_Eta5_.clear();
        FatjetPhi_R1p5_pT30_Eta5_.clear();
        mj_R0p8_pT30_Eta2p5_.clear();
        FatjetPt_R0p8_pT30_Eta2p5_.clear();
        FatjetEta_R0p8_pT30_Eta2p5_.clear();
        FatjetPhi_R0p8_pT30_Eta2p5_.clear();
        mj_R0p9_pT30_Eta2p5_.clear();
        FatjetPt_R0p9_pT30_Eta2p5_.clear();
        FatjetEta_R0p9_pT30_Eta2p5_.clear();
        FatjetPhi_R0p9_pT30_Eta2p5_.clear();
        mj_R1p0_pT30_Eta2p5_.clear();
        FatjetPt_R1p0_pT30_Eta2p5_.clear();
        FatjetEta_R1p0_pT30_Eta2p5_.clear();
        FatjetPhi_R1p0_pT30_Eta2p5_.clear();
        mj_R1p1_pT30_Eta2p5_.clear();
        FatjetPt_R1p1_pT30_Eta2p5_.clear();
        FatjetEta_R1p1_pT30_Eta2p5_.clear();
        FatjetPhi_R1p1_pT30_Eta2p5_.clear();
        mj_R1p2_pT30_Eta2p5_.clear();
        FatjetPt_R1p2_pT30_Eta2p5_.clear();
        FatjetEta_R1p2_pT30_Eta2p5_.clear();
        FatjetPhi_R1p2_pT30_Eta2p5_.clear();
        mj_R1p3_pT30_Eta2p5_.clear();
        FatjetPt_R1p3_pT30_Eta2p5_.clear();
        FatjetEta_R1p3_pT30_Eta2p5_.clear();
        FatjetPhi_R1p3_pT30_Eta2p5_.clear();
        mj_R1p4_pT30_Eta2p5_.clear();
        FatjetPt_R1p4_pT30_Eta2p5_.clear();
        FatjetEta_R1p4_pT30_Eta2p5_.clear();
        FatjetPhi_R1p4_pT30_Eta2p5_.clear();
        mj_R1p5_pT30_Eta2p5_.clear();
        FatjetPt_R1p5_pT30_Eta2p5_.clear();
        FatjetEta_R1p5_pT30_Eta2p5_.clear();
        FatjetPhi_R1p5_pT30_Eta2p5_.clear();
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
        if(!isData) 
        {
            EventWeight = Xsec/TotalNEntries*Lumi; // scale for 1 pb * Lumi
            // need more weights if needed 
        } 


        // Get good RA4 muons
        vector<int> RA4MuonVeto = GetMuons(false);
        vector<int> RA4Muon     = GetMuons(true);
        // Get good RA4 electrons
        vector<int> RA4ElecVeto = GetElectrons(false);
        vector<int> RA4Elec     = GetElectrons(true);
        // Get good skinny jets, HT and B-tagged jets 
        float HT=-999.; 
        vector<int> GoodJets_AK4 = GetJets(RA4Elec,RA4Muon,RA4ElecVeto,RA4MuonVeto, HT);


        // Skim : HT>750 MET>250
        if(HT<750) continue;
        if(mets_et->at(0)<250) continue;

        // Nbtag 
        int Ncsvm=0;
        int Ncsvl=0;
        for(int j=0; j<(int)GoodJets_AK4.size(); j++)
        { 
            if(jets_AK4_btag_secVertexCombined->at(GoodJets_AK4.at(j)) > /*0.814*/0.679) Ncsvm++; 
            if(jets_AK4_btag_secVertexCombined->at(GoodJets_AK4.at(j)) > /*0.423*/0.244) Ncsvl++; 
        }
        
        // 
        // pT(R=0.5) > 30 GeV
        // 

        // Regular jets
        double MJ_pT30=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjet_pT30; 
        vector<int> Vector_GoodFatjet_pT30_Index; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK4_R1p2_R0p5pT30_px->size(); ifatjet++) 
        {
            float temp_pT_pT30 = TMath::Sqrt(fastjets_AK4_R1p2_R0p5pT30_px->at(ifatjet)*fastjets_AK4_R1p2_R0p5pT30_px->at(ifatjet)
                                        +fastjets_AK4_R1p2_R0p5pT30_py->at(ifatjet)*fastjets_AK4_R1p2_R0p5pT30_py->at(ifatjet));
            if(temp_pT_pT30<50) continue;
            TLorentzVector temp_GoodFatjet_pT30( fastjets_AK4_R1p2_R0p5pT30_px->at(ifatjet), 
                                            fastjets_AK4_R1p2_R0p5pT30_py->at(ifatjet),
                                            fastjets_AK4_R1p2_R0p5pT30_pz->at(ifatjet),
                                            fastjets_AK4_R1p2_R0p5pT30_energy->at(ifatjet));
            Vector_GoodFatjet_pT30.push_back(temp_GoodFatjet_pT30);
            Vector_GoodFatjet_pT30_Index.push_back(ifatjet);
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
        // Fat jets on-the-fly vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv 
        // 
        vector<TLorentzVector> FatJetConstituent; 
        for(unsigned int ijet=0;ijet<jets_AK4_pt->size();ijet++)
        {
            TLorentzVector tmp(jets_AK4_px->at(ijet),jets_AK4_py->at(ijet),
                               jets_AK4_pz->at(ijet), jets_AK4_energy->at(ijet));
            FatJetConstituent.push_back(tmp);
        }

        //  R=0.8-1.5  pT(SJ)>30  |eta|<5
        vector<TLorentzVector>  FatJet_R0p8_pT30_Eta5 = makeFatJet(FatJetConstituent, 0.8, 30, 5);
        vector<TLorentzVector>  FatJet_R0p9_pT30_Eta5 = makeFatJet(FatJetConstituent, 0.9, 30, 5);
        vector<TLorentzVector>  FatJet_R1p0_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.0, 30, 5);
        vector<TLorentzVector>  FatJet_R1p1_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.1, 30, 5);
        vector<TLorentzVector>  FatJet_R1p2_pT30_Eta5 = makeFatJet(FatJetConstituent, FatJet_sjmj_R1p2_pT30_Eta5_, 1.2, 30, 5);
        vector<TLorentzVector>  FatJet_R1p3_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.3, 30, 5);
        vector<TLorentzVector>  FatJet_R1p4_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.4, 30, 5);
        vector<TLorentzVector>  FatJet_R1p5_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.5, 30, 5);
        for(unsigned int ifj=0; ifj<FatJet_R0p8_pT30_Eta5.size();ifj++)
        {
            mj_R0p8_pT30_Eta5_.push_back(FatJet_R0p8_pT30_Eta5.at(ifj).M());
            FatjetPt_R0p8_pT30_Eta5_.push_back(FatJet_R0p8_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R0p8_pT30_Eta5_.push_back(FatJet_R0p8_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R0p8_pT30_Eta5_.push_back(FatJet_R0p8_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R0p9_pT30_Eta5.size();ifj++)
        {
            mj_R0p9_pT30_Eta5_.push_back(FatJet_R0p9_pT30_Eta5.at(ifj).M());
            FatjetPt_R0p9_pT30_Eta5_.push_back(FatJet_R0p9_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R0p9_pT30_Eta5_.push_back(FatJet_R0p9_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R0p9_pT30_Eta5_.push_back(FatJet_R0p9_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p0_pT30_Eta5.size();ifj++)
        {
            mj_R1p0_pT30_Eta5_.push_back(FatJet_R1p0_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p0_pT30_Eta5_.push_back(FatJet_R1p0_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p0_pT30_Eta5_.push_back(FatJet_R1p0_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p0_pT30_Eta5_.push_back(FatJet_R1p0_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p1_pT30_Eta5.size();ifj++)
        {
            mj_R1p1_pT30_Eta5_.push_back(FatJet_R1p1_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p1_pT30_Eta5_.push_back(FatJet_R1p1_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p1_pT30_Eta5_.push_back(FatJet_R1p1_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p1_pT30_Eta5_.push_back(FatJet_R1p1_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT30_Eta5.size();ifj++)
        {
            mj_R1p2_pT30_Eta5_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p2_pT30_Eta5_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p2_pT30_Eta5_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p2_pT30_Eta5_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p3_pT30_Eta5.size();ifj++)
        {
            mj_R1p3_pT30_Eta5_.push_back(FatJet_R1p3_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p3_pT30_Eta5_.push_back(FatJet_R1p3_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p3_pT30_Eta5_.push_back(FatJet_R1p3_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p3_pT30_Eta5_.push_back(FatJet_R1p3_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p4_pT30_Eta5.size();ifj++)
        {
            mj_R1p4_pT30_Eta5_.push_back(FatJet_R1p4_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p4_pT30_Eta5_.push_back(FatJet_R1p4_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p4_pT30_Eta5_.push_back(FatJet_R1p4_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p4_pT30_Eta5_.push_back(FatJet_R1p4_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p5_pT30_Eta5.size();ifj++)
        {
            mj_R1p5_pT30_Eta5_.push_back(FatJet_R1p5_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p5_pT30_Eta5_.push_back(FatJet_R1p5_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p5_pT30_Eta5_.push_back(FatJet_R1p5_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p5_pT30_Eta5_.push_back(FatJet_R1p5_pT30_Eta5.at(ifj).Phi());
        }
        //  R=0.8-1.5  pT(SJ)>30  |eta|<2.5
        vector<TLorentzVector>  FatJet_R0p8_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 0.8, 30, 2.5);
        vector<TLorentzVector>  FatJet_R0p9_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 0.9, 30, 2.5);
        vector<TLorentzVector>  FatJet_R1p0_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 1.0, 30, 2.5);
        vector<TLorentzVector>  FatJet_R1p1_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 1.1, 30, 2.5);
        vector<TLorentzVector>  FatJet_R1p2_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 1.2, 30, 2.5);
        vector<TLorentzVector>  FatJet_R1p3_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 1.3, 30, 2.5);
        vector<TLorentzVector>  FatJet_R1p4_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 1.4, 30, 2.5);
        vector<TLorentzVector>  FatJet_R1p5_pT30_Eta2p5 = makeFatJet(FatJetConstituent, 1.5, 30, 2.5);
        for(unsigned int ifj=0; ifj<FatJet_R0p8_pT30_Eta2p5.size();ifj++)
        {
            mj_R0p8_pT30_Eta2p5_.push_back(FatJet_R0p8_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R0p8_pT30_Eta2p5_.push_back(FatJet_R0p8_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R0p8_pT30_Eta2p5_.push_back(FatJet_R0p8_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R0p8_pT30_Eta2p5_.push_back(FatJet_R0p8_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R0p9_pT30_Eta2p5.size();ifj++)
        {
            mj_R0p9_pT30_Eta2p5_.push_back(FatJet_R0p9_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R0p9_pT30_Eta2p5_.push_back(FatJet_R0p9_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R0p9_pT30_Eta2p5_.push_back(FatJet_R0p9_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R0p9_pT30_Eta2p5_.push_back(FatJet_R0p9_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p0_pT30_Eta2p5.size();ifj++)
        {
            mj_R1p0_pT30_Eta2p5_.push_back(FatJet_R1p0_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R1p0_pT30_Eta2p5_.push_back(FatJet_R1p0_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p0_pT30_Eta2p5_.push_back(FatJet_R1p0_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p0_pT30_Eta2p5_.push_back(FatJet_R1p0_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p1_pT30_Eta2p5.size();ifj++)
        {
            mj_R1p1_pT30_Eta2p5_.push_back(FatJet_R1p1_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R1p1_pT30_Eta2p5_.push_back(FatJet_R1p1_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p1_pT30_Eta2p5_.push_back(FatJet_R1p1_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p1_pT30_Eta2p5_.push_back(FatJet_R1p1_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT30_Eta2p5.size();ifj++)
        {
            mj_R1p2_pT30_Eta2p5_.push_back(FatJet_R1p2_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R1p2_pT30_Eta2p5_.push_back(FatJet_R1p2_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p2_pT30_Eta2p5_.push_back(FatJet_R1p2_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p2_pT30_Eta2p5_.push_back(FatJet_R1p2_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p3_pT30_Eta2p5.size();ifj++)
        {
            mj_R1p3_pT30_Eta2p5_.push_back(FatJet_R1p3_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R1p3_pT30_Eta2p5_.push_back(FatJet_R1p3_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p3_pT30_Eta2p5_.push_back(FatJet_R1p3_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p3_pT30_Eta2p5_.push_back(FatJet_R1p3_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p4_pT30_Eta2p5.size();ifj++)
        {
            mj_R1p4_pT30_Eta2p5_.push_back(FatJet_R1p4_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R1p4_pT30_Eta2p5_.push_back(FatJet_R1p4_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p4_pT30_Eta2p5_.push_back(FatJet_R1p4_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p4_pT30_Eta2p5_.push_back(FatJet_R1p4_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p5_pT30_Eta2p5.size();ifj++)
        {
            mj_R1p5_pT30_Eta2p5_.push_back(FatJet_R1p5_pT30_Eta2p5.at(ifj).M());
            FatjetPt_R1p5_pT30_Eta2p5_.push_back(FatJet_R1p5_pT30_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p5_pT30_Eta2p5_.push_back(FatJet_R1p5_pT30_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p5_pT30_Eta2p5_.push_back(FatJet_R1p5_pT30_Eta2p5.at(ifj).Phi());
        }
        // Fat jets on-the-fly ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
       
        //
        // Fill the baby variables 
        //
        run_                =   run;
        lumiblock_          =   lumiblock;
        event_              =   event; 
        //TrigMuon_           =   PassMuonTrig(); 
        //TrigHTMuon_         =   PassHTMuonTrig(); 
        //TrigSingleMuon_     =   PassSingleMuonTrig(); 
        Nfatjet_pT30_       =   Nfatjet_pT30;
        Nskinnyjet_         =   GoodJets_AK4.size();
        NBtagCSVM_          =   Ncsvm;
        NBtagCSVL_          =   Ncsvl;
        EventWeight_        =   EventWeight;
        MJ_pT30_            =   MJ_pT30;
        MET_                =   mets_et->at(0);
        METPhi_             =   mets_phi->at(0);
        HT_                 =   HT;
        Npv_                =   Npv;
        for(unsigned int bc(0); bc<PU_bunchCrossing->size(); ++bc)
        {
            if(PU_bunchCrossing->at(bc)==0)
            {
                Npu_ = PU_TrueNumInteractions->at(bc);
                break;
            }
        }
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT30.size(); igoodfatjet++) 
        { 
            mj_pT30_.push_back(Vector_mj_pT30.at(igoodfatjet));
            FatjetPt_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Pt());
            FatjetEta_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Eta());
            FatjetPhi_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Phi());
            FatjetN_pT30_.push_back(fastjets_AK4_R1p2_R0p5pT30_nconstituents->at(Vector_GoodFatjet_pT30_Index.at(igoodfatjet)));
        }

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
        
        for(int igoodjet=0; igoodjet<(int)GoodJets_AK4.size(); igoodjet++) 
        {
          int ijet = GoodJets_AK4.at(igoodjet); 
          JetPt_.push_back(jets_AK4_pt->at(ijet)); 
          JetEta_.push_back(jets_AK4_eta->at(ijet)); 
          JetPhi_.push_back(jets_AK4_phi->at(ijet)); 
          JetCSV_.push_back(jets_AK4_btag_secVertexCombined->at(ijet)); 
        }
        // gen top pT for TTbar samples  
        for(int igen=0; igen<(int)mc_doc_id->size(); igen++) 
        { 
            GenPt_.push_back(   mc_doc_pt->at(igen));  
            GenPhi_.push_back(  mc_doc_phi->at(igen));  
            GenEta_.push_back(  mc_doc_eta->at(igen));  
            GenId_.push_back(   mc_doc_id->at(igen));  
            GenMId_.push_back(  mc_doc_mother_id->at(igen));  
            GenGMId_.push_back( mc_doc_grandmother_id->at(igen));  

            if(ProcessName.Contains("TT")) 
            {
                if(mc_doc_id->at(igen)==6)  
                {
                    top1pT_ = mc_doc_pt->at(igen);  
                    top1Phi_ = mc_doc_phi->at(igen);  
                    top1Eta_ = mc_doc_eta->at(igen);  
                }
                if(mc_doc_id->at(igen)==-6)
                {
                    top2pT_ = mc_doc_pt->at(igen);  
                    top2Phi_ = mc_doc_phi->at(igen);  
                    top2Eta_ = mc_doc_eta->at(igen);  
                }
            }
        }

        babyTree_->Fill(); // Fill all events

        //for(int i=0; i<GoodJets_AK4.size(); i++) cout << event << " :: " << HT << endl;
        
        // Clean fat jets? 
        // (1) identify if all skinny jets associated with a given fat jet 
        //     are good jets by comparing index of skinny jet in GoodJets_AK4
        //     and fastjet_.._R1p2pTxx_index
        // (2) if a skinny jet is not in the good jet, its four vector is subtracted 
        //     from the fat jet and mj is calculated again
        // (3) MJ is then calculated

    } // event loop
    cout << endl;
    cout << "[MJ Analysis] Looping over events has been done" << endl;
    
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
    elete chainA;
    delete chainB;
}
