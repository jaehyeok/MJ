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

//#include "Branch_v75.h" 
//#include "ra4_objects_13TeV.h"
#include "Branch/Branch_v77.h" 
#include "core/ra4_objects_13TeV.h"
#include "core/Utilities_13TeV.h"

// include necessary fastjet files
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

bool isLO(TString ProcessName)
{
    // single top and TTH samples are generated with aMC@NLO
    // and require special treatment of negative weights
    if(ProcessName.Contains("amcatnlo")) return false;
    else if(ProcessName.Contains("aMCatNLO")) return false;
    else if(ProcessName.Contains("TToLeptons")) return false;
    else if(ProcessName.Contains("TBarToLeptons")) return false;
    else if(ProcessName.Contains("TTH")) return false;
    else return true;

}

// only care about the sign of the weight;
// the average of the weights is the cross section
int getWeights(TChain *chB)
{
    // gimmicky way to quickly get number of positive and negative weights
    // without having to call SetBranchStatus manually
    TH1F *weights = new TH1F("weights", "weights", 2, -1e9, 1e9);
    chB->Project("weights", "weight");

    float weightSum = weights->GetBinContent(2)-weights->GetBinContent(1);
    cout << weightSum << " events after taking into account negative weights" << endl;

    return (int)weightSum;
}

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
    TChain * chainA = new TChain("/cfA/eventA");   
    TChain * chainB = new TChain("/cfA/eventB");  
    
    for(int i=ibegin; i<=iend; i++)  
    {
        gSystem->Exec(Form("ls %s/*_f%i.root", InputName.Data(), i));
        chainA->Add(Form("%s/*_f%i.root", InputName.Data(), i));
        chainB->Add(Form("%s/*_f%i.root", InputName.Data(), i));
    } 
    TList *l = (TList*)chainA->GetListOfFiles();
    l->Print();

    InitializeA(chainA);
    InitializeB(chainB);
    
    // Get total number of events of a sample 
    int TotalNEntries=1;
    int TotalNEntriesNeg=1;
    
    if(!isData)
    { 
        TChain * chainATotal = new TChain("/cfA/eventA");  
        chainATotal->Add(Form("%s/*.root", InputName.Data()));
        TotalNEntries = (int)chainATotal->GetEntries();
        if(isLO(ProcessName)) TotalNEntriesNeg = (int)chainATotal->GetEntries();
        else TotalNEntriesNeg = getWeights(chainB);
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
    int Npuminusone_; 
    int Npuplusone_; 
    float Npu_; 
    float EventWeight_; 
    float EventWeightNeg_; 
    float MJ_pT30_; 
    float MET_;  
    float METPhi_;  
    float GenMET_;  
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
    /*
    vector<float> mj_R1p6_pT30_Eta5_;
    vector<float> FatjetPt_R1p6_pT30_Eta5_;
    vector<float> FatjetEta_R1p6_pT30_Eta5_;
    vector<float> FatjetPhi_R1p6_pT30_Eta5_;
    vector<float> mj_R1p7_pT30_Eta5_;
    vector<float> FatjetPt_R1p7_pT30_Eta5_;
    vector<float> FatjetEta_R1p7_pT30_Eta5_;
    vector<float> FatjetPhi_R1p7_pT30_Eta5_;
    vector<float> mj_R1p8_pT30_Eta5_;
    vector<float> FatjetPt_R1p8_pT30_Eta5_;
    vector<float> FatjetEta_R1p8_pT30_Eta5_;
    vector<float> FatjetPhi_R1p8_pT30_Eta5_;
    */
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
    vector<float> mj_R1p2_pT35_Eta2p5_;
    vector<float> FatjetPt_R1p2_pT35_Eta2p5_;
    vector<float> FatjetEta_R1p2_pT35_Eta2p5_;
    vector<float> FatjetPhi_R1p2_pT35_Eta2p5_; 
    vector<float> mj_R1p2_pT40_Eta2p5_;
    vector<float> FatjetPt_R1p2_pT40_Eta2p5_;
    vector<float> FatjetEta_R1p2_pT40_Eta2p5_;
    vector<float> FatjetPhi_R1p2_pT40_Eta2p5_; 
    vector<float> mj_R1p2_pT45_Eta2p5_;
    vector<float> FatjetPt_R1p2_pT45_Eta2p5_;
    vector<float> FatjetEta_R1p2_pT45_Eta2p5_;
    vector<float> FatjetPhi_R1p2_pT45_Eta2p5_; 
    vector<float> mj_R1p2_pT50_Eta2p5_;
    vector<float> FatjetPt_R1p2_pT50_Eta2p5_;
    vector<float> FatjetEta_R1p2_pT50_Eta2p5_;
    vector<float> FatjetPhi_R1p2_pT50_Eta2p5_; 
    vector<float> mj_R1p2_pT60_Eta2p5_;
    vector<float> FatjetPt_R1p2_pT60_Eta2p5_;
    vector<float> FatjetEta_R1p2_pT60_Eta2p5_;
    vector<float> FatjetPhi_R1p2_pT60_Eta2p5_; 
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
    vector<float> JetE_;
    vector<float> JetPt_;
    vector<float> JetEta_;
    vector<float> JetPhi_;
    vector<float> JetCSV_;
    vector<float> RA4ElsPt_mi_;
    vector<float> RA4ElsEta_mi_;
    vector<float> RA4ElsPhi_mi_;
    vector<float> RA4MusPt_mi_;
    vector<float> RA4MusEta_mi_;
    vector<float> RA4MusPhi_mi_;
    vector<float> RA4ElsVetoPt_mi_;
    vector<float> RA4ElsVetoEta_mi_;
    vector<float> RA4ElsVetoPhi_mi_;
    vector<float> RA4MusVetoPt_mi_;
    vector<float> RA4MusVetoEta_mi_;
    vector<float> RA4MusVetoPhi_mi_;
    vector<float> IsoTrkVetoElsPt_;
    vector<float> IsoTrkVetoElsEta_;
    vector<float> IsoTrkVetoElsPhi_;
    vector<float> IsoTrkVetoMusPt_;
    vector<float> IsoTrkVetoMusEta_;
    vector<float> IsoTrkVetoMusPhi_;
    vector<float> IsoTrkVetoHadPt_;
    vector<float> IsoTrkVetoHadEta_;
    vector<float> IsoTrkVetoHadPhi_;
    vector<float> JetE_mi_;
    vector<float> JetPt_mi_;
    vector<float> JetEta_mi_;
    vector<float> JetPhi_mi_;
    vector<float> JetCSV_mi_;
    vector<float> GenPt_;
    vector<float> GenEta_;
    vector<float> GenPhi_;
    vector<float> GenStatus_;
    vector<float> GenId_;
    vector<float> GenMId_;
    vector<float> GenGMId_;
    vector<float> GenJetPt_;
    vector<float> GenJetEta_;
    vector<float> GenJetPhi_;
    vector<float> MCJetE_;
    vector<float> MCJetPt_;
    vector<float> MCJetEta_;
    vector<float> MCJetPhi_;
    vector<float> Genmj_R0p8_pT30_Eta2p5_;
    vector<float> GenFatjetPt_R0p8_pT30_Eta2p5_;
    vector<float> GenFatjetEta_R0p8_pT30_Eta2p5_;
    vector<float> GenFatjetPhi_R0p8_pT30_Eta2p5_;
    vector<float> Genmj_R1p0_pT30_Eta2p5_;
    vector<float> GenFatjetPt_R1p0_pT30_Eta2p5_;
    vector<float> GenFatjetEta_R1p0_pT30_Eta2p5_;
    vector<float> GenFatjetPhi_R1p0_pT30_Eta2p5_;
    vector<float> Genmj_R1p2_pT30_Eta2p5_;
    vector<float> GenFatjetPt_R1p2_pT30_Eta2p5_;
    vector<float> GenFatjetEta_R1p2_pT30_Eta2p5_;
    vector<float> GenFatjetPhi_R1p2_pT30_Eta2p5_;
    vector<float> Genmj_R1p4_pT30_Eta2p5_;
    vector<float> GenFatjetPt_R1p4_pT30_Eta2p5_;
    vector<float> GenFatjetEta_R1p4_pT30_Eta2p5_;
    vector<float> GenFatjetPhi_R1p4_pT30_Eta2p5_;
    vector<float> Genmj_R1p2_pT35_Eta2p5_;
    vector<float> GenFatjetPt_R1p2_pT35_Eta2p5_;
    vector<float> GenFatjetEta_R1p2_pT35_Eta2p5_;
    vector<float> GenFatjetPhi_R1p2_pT35_Eta2p5_;
    vector<float> Genmj_R1p2_pT40_Eta2p5_;
    vector<float> GenFatjetPt_R1p2_pT40_Eta2p5_;
    vector<float> GenFatjetEta_R1p2_pT40_Eta2p5_;
    vector<float> GenFatjetPhi_R1p2_pT40_Eta2p5_;
    vector<float> Genmj_R1p2_pT45_Eta2p5_;
    vector<float> GenFatjetPt_R1p2_pT45_Eta2p5_;
    vector<float> GenFatjetEta_R1p2_pT45_Eta2p5_;
    vector<float> GenFatjetPhi_R1p2_pT45_Eta2p5_;
    vector<float> Genmj_R1p2_pT50_Eta2p5_;
    vector<float> GenFatjetPt_R1p2_pT50_Eta2p5_;
    vector<float> GenFatjetEta_R1p2_pT50_Eta2p5_;
    vector<float> GenFatjetPhi_R1p2_pT50_Eta2p5_;
    vector<float> Genmj_R1p2_pT60_Eta2p5_;
    vector<float> GenFatjetPt_R1p2_pT60_Eta2p5_;
    vector<float> GenFatjetEta_R1p2_pT60_Eta2p5_;
    vector<float> GenFatjetPhi_R1p2_pT60_Eta2p5_;
    
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
    babyTree_->Branch("Npuminusone",       	&Npuminusone_);       
    babyTree_->Branch("Npuplusone",        	&Npuplusone_);       
    babyTree_->Branch("Npu",            	&Npu_);       
    babyTree_->Branch("EventWeight",    	&EventWeight_);
    babyTree_->Branch("EventWeightNeg",    	&EventWeightNeg_);
    babyTree_->Branch("MJ_pT30",        	&MJ_pT30_);        
    babyTree_->Branch("MET",            	&MET_);        
    babyTree_->Branch("METPhi",            	&METPhi_);        
    babyTree_->Branch("GenMET",            	&GenMET_);        
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
    /*
    babyTree_->Branch("mj_R1p6_pT30_Eta5",        	&mj_R1p6_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p6_pT30_Eta5",  	&FatjetPt_R1p6_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p6_pT30_Eta5", 	&FatjetEta_R1p6_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p6_pT30_Eta5",     &FatjetPhi_R1p6_pT30_Eta5_);
    babyTree_->Branch("mj_R1p7_pT30_Eta5",        	&mj_R1p7_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p7_pT30_Eta5",  	&FatjetPt_R1p7_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p7_pT30_Eta5", 	&FatjetEta_R1p7_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p7_pT30_Eta5",     &FatjetPhi_R1p7_pT30_Eta5_);
    babyTree_->Branch("mj_R1p8_pT30_Eta5",        	&mj_R1p8_pT30_Eta5_);     
    babyTree_->Branch("FatjetPt_R1p8_pT30_Eta5",  	&FatjetPt_R1p8_pT30_Eta5_); 
    babyTree_->Branch("FatjetEta_R1p8_pT30_Eta5", 	&FatjetEta_R1p8_pT30_Eta5_);
    babyTree_->Branch("FatjetPhi_R1p8_pT30_Eta5",     &FatjetPhi_R1p8_pT30_Eta5_);
    */
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
    babyTree_->Branch("mj_R1p2_pT35_Eta2p5",        	&mj_R1p2_pT35_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p2_pT35_Eta2p5",  	&FatjetPt_R1p2_pT35_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p2_pT35_Eta2p5", 	&FatjetEta_R1p2_pT35_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p2_pT35_Eta2p5",     &FatjetPhi_R1p2_pT35_Eta2p5_);
    babyTree_->Branch("mj_R1p2_pT40_Eta2p5",        	&mj_R1p2_pT40_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p2_pT40_Eta2p5",  	&FatjetPt_R1p2_pT40_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p2_pT40_Eta2p5", 	&FatjetEta_R1p2_pT40_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p2_pT40_Eta2p5",     &FatjetPhi_R1p2_pT40_Eta2p5_);
    babyTree_->Branch("mj_R1p2_pT45_Eta2p5",        	&mj_R1p2_pT45_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p2_pT45_Eta2p5",  	&FatjetPt_R1p2_pT45_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p2_pT45_Eta2p5", 	&FatjetEta_R1p2_pT45_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p2_pT45_Eta2p5",     &FatjetPhi_R1p2_pT45_Eta2p5_);
    babyTree_->Branch("mj_R1p2_pT50_Eta2p5",        	&mj_R1p2_pT50_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p2_pT50_Eta2p5",  	&FatjetPt_R1p2_pT50_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p2_pT50_Eta2p5", 	&FatjetEta_R1p2_pT50_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p2_pT50_Eta2p5",     &FatjetPhi_R1p2_pT50_Eta2p5_);
    babyTree_->Branch("mj_R1p2_pT60_Eta2p5",        	&mj_R1p2_pT60_Eta2p5_);     
    babyTree_->Branch("FatjetPt_R1p2_pT60_Eta2p5",  	&FatjetPt_R1p2_pT60_Eta2p5_); 
    babyTree_->Branch("FatjetEta_R1p2_pT60_Eta2p5", 	&FatjetEta_R1p2_pT60_Eta2p5_);
    babyTree_->Branch("FatjetPhi_R1p2_pT60_Eta2p5",     &FatjetPhi_R1p2_pT60_Eta2p5_);
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
    babyTree_->Branch("JetE",              &JetE_);
    babyTree_->Branch("JetPt",              &JetPt_);
    babyTree_->Branch("JetEta",             &JetEta_);
    babyTree_->Branch("JetPhi",             &JetPhi_);
    babyTree_->Branch("JetCSV",             &JetCSV_);
    babyTree_->Branch("RA4ElsPt_mi",        &RA4ElsPt_mi_);
    babyTree_->Branch("RA4ElsEta_mi",       &RA4ElsEta_mi_);
    babyTree_->Branch("RA4ElsPhi_mi",       &RA4ElsPhi_mi_);
    babyTree_->Branch("RA4MusPt_mi",        &RA4MusPt_mi_);
    babyTree_->Branch("RA4MusEta_mi",       &RA4MusEta_mi_);
    babyTree_->Branch("RA4MusPhi_mi",       &RA4MusPhi_mi_);
    babyTree_->Branch("RA4ElsVetoPt_mi",    &RA4ElsVetoPt_mi_);
    babyTree_->Branch("RA4ElsVetoEta_mi",   &RA4ElsVetoEta_mi_);
    babyTree_->Branch("RA4ElsVetoPhi_mi",   &RA4ElsVetoPhi_mi_);
    babyTree_->Branch("RA4MusVetoPt_mi",    &RA4MusVetoPt_mi_);
    babyTree_->Branch("RA4MusVetoEta_mi",   &RA4MusVetoEta_mi_);
    babyTree_->Branch("RA4MusVetoPhi_mi",   &RA4MusVetoPhi_mi_);
    babyTree_->Branch("IsoTrkVetoElsPt",    &IsoTrkVetoElsPt_);
    babyTree_->Branch("IsoTrkVetoElsEta",   &IsoTrkVetoElsEta_);
    babyTree_->Branch("IsoTrkVetoElsPhi",   &IsoTrkVetoElsPhi_);
    babyTree_->Branch("IsoTrkVetoMusPt",    &IsoTrkVetoMusPt_);
    babyTree_->Branch("IsoTrkVetoMusEta",   &IsoTrkVetoMusEta_);
    babyTree_->Branch("IsoTrkVetoMusPhi",   &IsoTrkVetoMusPhi_);
    babyTree_->Branch("IsoTrkVetoHadPt",    &IsoTrkVetoHadPt_);
    babyTree_->Branch("IsoTrkVetoHadEta",   &IsoTrkVetoHadEta_);
    babyTree_->Branch("IsoTrkVetoHadPhi",   &IsoTrkVetoHadPhi_);
    babyTree_->Branch("JetE_mi",            &JetE_mi_);
    babyTree_->Branch("JetPt_mi",           &JetPt_mi_);
    babyTree_->Branch("JetEta_mi",          &JetEta_mi_);
    babyTree_->Branch("JetPhi_mi",          &JetPhi_mi_);
    babyTree_->Branch("JetCSV_mi",          &JetCSV_mi_);
    babyTree_->Branch("GenPt",              &GenPt_);
    babyTree_->Branch("GenEta",             &GenEta_);
    babyTree_->Branch("GenPhi",             &GenPhi_);
    babyTree_->Branch("GenStatus",          &GenStatus_);
    babyTree_->Branch("GenId",              &GenId_);
    babyTree_->Branch("GenMId",             &GenMId_);
    babyTree_->Branch("GenGMId",            &GenGMId_);
    babyTree_->Branch("GenJetPt",           &GenJetPt_);
    babyTree_->Branch("GenJetEta",          &GenJetEta_);
    babyTree_->Branch("GenJetPhi",          &GenJetPhi_);
    babyTree_->Branch("MCJetE",             &MCJetE_);
    babyTree_->Branch("MCJetPt",            &MCJetPt_);
    babyTree_->Branch("MCJetEta",           &MCJetEta_);
    babyTree_->Branch("MCJetPhi",           &MCJetPhi_);
    babyTree_->Branch("Genmj_R0p8_pT30_Eta2p5",        	&Genmj_R0p8_pT30_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R0p8_pT30_Eta2p5",  	&GenFatjetPt_R0p8_pT30_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R0p8_pT30_Eta2p5", 	&GenFatjetEta_R0p8_pT30_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R0p8_pT30_Eta2p5",  &GenFatjetPhi_R0p8_pT30_Eta2p5_);
    babyTree_->Branch("Genmj_R1p0_pT30_Eta2p5",        	&Genmj_R1p0_pT30_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p0_pT30_Eta2p5",  	&GenFatjetPt_R1p0_pT30_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p0_pT30_Eta2p5", 	&GenFatjetEta_R1p0_pT30_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p0_pT30_Eta2p5",  &GenFatjetPhi_R1p0_pT30_Eta2p5_);
    babyTree_->Branch("Genmj_R1p2_pT30_Eta2p5",        	&Genmj_R1p2_pT30_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p2_pT30_Eta2p5",  	&GenFatjetPt_R1p2_pT30_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p2_pT30_Eta2p5", 	&GenFatjetEta_R1p2_pT30_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p2_pT30_Eta2p5",  &GenFatjetPhi_R1p2_pT30_Eta2p5_);
    babyTree_->Branch("Genmj_R1p4_pT30_Eta2p5",        	&Genmj_R1p4_pT30_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p4_pT30_Eta2p5",  	&GenFatjetPt_R1p4_pT30_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p4_pT30_Eta2p5", 	&GenFatjetEta_R1p4_pT30_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p4_pT30_Eta2p5",  &GenFatjetPhi_R1p4_pT30_Eta2p5_);
    babyTree_->Branch("Genmj_R1p2_pT35_Eta2p5",        	&Genmj_R1p2_pT35_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p2_pT35_Eta2p5",  	&GenFatjetPt_R1p2_pT35_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p2_pT35_Eta2p5", 	&GenFatjetEta_R1p2_pT35_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p2_pT35_Eta2p5",  &GenFatjetPhi_R1p2_pT35_Eta2p5_);
    babyTree_->Branch("Genmj_R1p2_pT40_Eta2p5",        	&Genmj_R1p2_pT40_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p2_pT40_Eta2p5",  	&GenFatjetPt_R1p2_pT40_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p2_pT40_Eta2p5", 	&GenFatjetEta_R1p2_pT40_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p2_pT40_Eta2p5",  &GenFatjetPhi_R1p2_pT40_Eta2p5_);
    babyTree_->Branch("Genmj_R1p2_pT45_Eta2p5",        	&Genmj_R1p2_pT45_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p2_pT45_Eta2p5",  	&GenFatjetPt_R1p2_pT45_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p2_pT45_Eta2p5", 	&GenFatjetEta_R1p2_pT45_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p2_pT45_Eta2p5",  &GenFatjetPhi_R1p2_pT45_Eta2p5_);
    babyTree_->Branch("Genmj_R1p2_pT50_Eta2p5",        	&Genmj_R1p2_pT50_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p2_pT50_Eta2p5",  	&GenFatjetPt_R1p2_pT50_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p2_pT50_Eta2p5", 	&GenFatjetEta_R1p2_pT50_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p2_pT50_Eta2p5",  &GenFatjetPhi_R1p2_pT50_Eta2p5_);
    babyTree_->Branch("Genmj_R1p2_pT60_Eta2p5",        	&Genmj_R1p2_pT60_Eta2p5_);     
    babyTree_->Branch("GenFatjetPt_R1p2_pT60_Eta2p5",  	&GenFatjetPt_R1p2_pT60_Eta2p5_); 
    babyTree_->Branch("GenFatjetEta_R1p2_pT60_Eta2p5", 	&GenFatjetEta_R1p2_pT60_Eta2p5_);
    babyTree_->Branch("GenFatjetPhi_R1p2_pT60_Eta2p5",  &GenFatjetPhi_R1p2_pT60_Eta2p5_);
   
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
    cout<<"[MJ Analysis] Number of entries of total sample is(taking into account negative weights in aMC@NLO samples): "<<TotalNEntriesNeg<<endl;
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
        Npuminusone_        =   -1;
        Npuplusone_         =   -1;
        Npu_                =   -1;
        EventWeight_        =   1.;
        EventWeightNeg_     =   1.;
        MJ_pT30_            =-999.;
        MET_                =-999.;
        METPhi_             =-999.;
        GenMET_             =-999.;
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
        /*
        mj_R1p6_pT30_Eta5_.clear();
        FatjetPt_R1p6_pT30_Eta5_.clear();
        FatjetEta_R1p6_pT30_Eta5_.clear();
        FatjetPhi_R1p6_pT30_Eta5_.clear();
        mj_R1p7_pT30_Eta5_.clear();
        FatjetPt_R1p7_pT30_Eta5_.clear();
        FatjetEta_R1p7_pT30_Eta5_.clear();
        FatjetPhi_R1p7_pT30_Eta5_.clear();
        mj_R1p8_pT30_Eta5_.clear();
        FatjetPt_R1p8_pT30_Eta5_.clear();
        FatjetEta_R1p8_pT30_Eta5_.clear();
        FatjetPhi_R1p8_pT30_Eta5_.clear();
        */
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
        mj_R1p2_pT35_Eta2p5_.clear();
        FatjetPt_R1p2_pT35_Eta2p5_.clear();
        FatjetEta_R1p2_pT35_Eta2p5_.clear();
        FatjetPhi_R1p2_pT35_Eta2p5_.clear();
        mj_R1p2_pT40_Eta2p5_.clear();
        FatjetPt_R1p2_pT40_Eta2p5_.clear();
        FatjetEta_R1p2_pT40_Eta2p5_.clear();
        FatjetPhi_R1p2_pT40_Eta2p5_.clear();
        mj_R1p2_pT45_Eta2p5_.clear();
        FatjetPt_R1p2_pT45_Eta2p5_.clear();
        FatjetEta_R1p2_pT45_Eta2p5_.clear();
        FatjetPhi_R1p2_pT45_Eta2p5_.clear();
        mj_R1p2_pT50_Eta2p5_.clear();
        FatjetPt_R1p2_pT50_Eta2p5_.clear();
        FatjetEta_R1p2_pT50_Eta2p5_.clear();
        FatjetPhi_R1p2_pT50_Eta2p5_.clear();
        mj_R1p2_pT60_Eta2p5_.clear();
        FatjetPt_R1p2_pT60_Eta2p5_.clear();
        FatjetEta_R1p2_pT60_Eta2p5_.clear();
        FatjetPhi_R1p2_pT60_Eta2p5_.clear();
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
        JetE_.clear();
        JetPt_.clear();
        JetEta_.clear();
        JetPhi_.clear();
        JetCSV_.clear();
        RA4ElsPt_mi_.clear();
        RA4ElsEta_mi_.clear();
        RA4ElsPhi_mi_.clear();
        RA4MusPt_mi_.clear();
        RA4MusEta_mi_.clear();
        RA4MusPhi_mi_.clear();
        RA4ElsVetoPt_mi_.clear();
        RA4ElsVetoEta_mi_.clear();
        RA4ElsVetoPhi_mi_.clear();
        RA4MusVetoPt_mi_.clear();
        RA4MusVetoEta_mi_.clear();
        RA4MusVetoPhi_mi_.clear();
        IsoTrkVetoElsPt_.clear();
        IsoTrkVetoElsEta_.clear();
        IsoTrkVetoElsPhi_.clear();
        IsoTrkVetoMusPt_.clear();
        IsoTrkVetoMusEta_.clear();
        IsoTrkVetoMusPhi_.clear();
        IsoTrkVetoHadPt_.clear();
        IsoTrkVetoHadEta_.clear();
        IsoTrkVetoHadPhi_.clear();
        JetE_mi_.clear();
        JetPt_mi_.clear();
        JetEta_mi_.clear();
        JetPhi_mi_.clear();
        JetCSV_mi_.clear();
        GenPt_.clear();
        GenEta_.clear();
        GenPhi_.clear();
        GenStatus_.clear();
        GenId_.clear();
        GenMId_.clear();
        GenGMId_.clear();
        GenJetPt_.clear();
        GenJetEta_.clear();
        GenJetPhi_.clear();
        MCJetE_.clear();
        MCJetPt_.clear();
        MCJetEta_.clear();
        MCJetPhi_.clear();
        Genmj_R0p8_pT30_Eta2p5_.clear();
        GenFatjetPt_R0p8_pT30_Eta2p5_.clear();
        GenFatjetEta_R0p8_pT30_Eta2p5_.clear();
        GenFatjetPhi_R0p8_pT30_Eta2p5_.clear();
        Genmj_R1p0_pT30_Eta2p5_.clear();
        GenFatjetPt_R1p0_pT30_Eta2p5_.clear();
        GenFatjetEta_R1p0_pT30_Eta2p5_.clear();
        GenFatjetPhi_R1p0_pT30_Eta2p5_.clear();
        Genmj_R1p2_pT30_Eta2p5_.clear();
        GenFatjetPt_R1p2_pT30_Eta2p5_.clear();
        GenFatjetEta_R1p2_pT30_Eta2p5_.clear();
        GenFatjetPhi_R1p2_pT30_Eta2p5_.clear();
        Genmj_R1p4_pT30_Eta2p5_.clear();
        GenFatjetPt_R1p4_pT30_Eta2p5_.clear();
        GenFatjetEta_R1p4_pT30_Eta2p5_.clear();
        GenFatjetPhi_R1p4_pT30_Eta2p5_.clear();
        Genmj_R1p2_pT35_Eta2p5_.clear();
        GenFatjetPt_R1p2_pT35_Eta2p5_.clear();
        GenFatjetEta_R1p2_pT35_Eta2p5_.clear();
        GenFatjetPhi_R1p2_pT35_Eta2p5_.clear();
        Genmj_R1p2_pT40_Eta2p5_.clear();
        GenFatjetPt_R1p2_pT40_Eta2p5_.clear();
        GenFatjetEta_R1p2_pT40_Eta2p5_.clear();
        GenFatjetPhi_R1p2_pT40_Eta2p5_.clear();
        Genmj_R1p2_pT45_Eta2p5_.clear();
        GenFatjetPt_R1p2_pT45_Eta2p5_.clear();
        GenFatjetEta_R1p2_pT45_Eta2p5_.clear();
        GenFatjetPhi_R1p2_pT45_Eta2p5_.clear();
        Genmj_R1p2_pT50_Eta2p5_.clear();
        GenFatjetPt_R1p2_pT50_Eta2p5_.clear();
        GenFatjetEta_R1p2_pT50_Eta2p5_.clear();
        GenFatjetPhi_R1p2_pT50_Eta2p5_.clear();
        Genmj_R1p2_pT60_Eta2p5_.clear();
        GenFatjetPt_R1p2_pT60_Eta2p5_.clear();
        GenFatjetEta_R1p2_pT60_Eta2p5_.clear();
        GenFatjetPhi_R1p2_pT60_Eta2p5_.clear();

        //
        // Core analysis 
        //
        // Get event weight 
        float EventWeight = 1;
        float EventWeightNeg = 1;
        if(!isData) 
        {
            EventWeight         = Xsec/TotalNEntries*Lumi; // scale for 1 pb * Lumi
            EventWeightNeg      = Xsec/TotalNEntriesNeg*Lumi; // scale for 1 pb * Lumi
            // take into account negative weights in NLO samples
            if(weight<0) EventWeightNeg*=-1.0;
            // need more weights if needed 
        } 
    
        //
        // Reliso 
        //
        // Get good RA4 muons
        vector<int> RA4MuonVeto = GetMuons(false, false);
        vector<int> RA4Muon     = GetMuons(true, false);
        // Get good RA4 electrons
        vector<int> RA4ElecVeto = GetElectrons(false, false);
        vector<int> RA4Elec     = GetElectrons(true, false);
        // Get good skinny jets, HT and B-tagged jets 
        float HT=-999.; 
        vector<int> GoodJets_AK4 = GetJets(RA4Elec,RA4Muon,RA4ElecVeto,RA4MuonVeto,HT);
        // Nbtag 
        int Ncsvm=0;
        int Ncsvl=0;
        for(int j=0; j<(int)GoodJets_AK4.size(); j++)
        { 
            if(jets_AK4_btag_inc_secVertexCombined->at(GoodJets_AK4.at(j)) > 0.814/*0.679*/) Ncsvm++; 
            if(jets_AK4_btag_inc_secVertexCombined->at(GoodJets_AK4.at(j)) > 0.423/*0.244*/) Ncsvl++; 
        }
        
        //
        // MiniIso 
        //
        // Get good RA4 muons
        vector<int> RA4MuonVeto_mini = GetMuons(false,  true);
        vector<int> RA4Muon_mini     = GetMuons(true,   true);
        // Get good RA4 electrons
        vector<int> RA4ElecVeto_mini = GetElectrons(false,  true);
        vector<int> RA4Elec_mini     = GetElectrons(true,   true);
        // Get good skinny jets, HT and B-tagged jets 
        float HT_mini=-999.; 
        vector<int> GoodJets_AK4_mini = GetJets(RA4Elec_mini,RA4Muon_mini,RA4ElecVeto_mini,RA4MuonVeto_mini,HT_mini);
        // Nbtag 
        int Ncsvm_mini=0;
        int Ncsvl_mini=0;
        for(int jmini=0; jmini<(int)GoodJets_AK4_mini.size(); jmini++)
        { 
            if(jets_AK4_btag_inc_secVertexCombined->at(GoodJets_AK4_mini.at(jmini)) > 0.814/*0.679*/) Ncsvm_mini++; 
            if(jets_AK4_btag_inc_secVertexCombined->at(GoodJets_AK4_mini.at(jmini)) > 0.423/*0.244*/) Ncsvl_mini++; 
        }

        // Iso track veto 
        std::vector<std::pair<int,double> > eCands;
        std::vector<std::pair<int,double> > muCands;
        std::vector<std::pair<int,double> > hadCands;
        GetIsoTracks(eCands,muCands,hadCands);

        // Skim : HT>500 MET>250
        if( HT<500 ) continue;
        //if(pfType1mets_et->at(0)<200) continue;
        if( (RA4Muon_mini.size()+RA4Elec_mini.size())<2 && (RA4Muon.size()+RA4Elec.size())<2 ) continue;

        // 
        // pT(R=0.5) > 30 GeV
        // 

        // Regular jets
        double MJ_pT30=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjet_pT30; 
        vector<int> Vector_GoodFatjet_pT30_Index; 
        /*
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
        */

        vector<float> Vector_mj_pT30;   // mj
        for(int ifj=0; ifj<fjets30_m->size(); ifj++) 
        {
            Vector_mj_pT30.push_back(fjets30_m->at(ifj));
        }
        MJ_pT30 = GetMJ(Vector_mj_pT30);

        int Nfatjet_pT30 = fjets30_m->size(); 
        
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
        vector<TLorentzVector>  FatJet_R1p2_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.2, 30, 5);
        vector<TLorentzVector>  FatJet_R1p3_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.3, 30, 5);
        vector<TLorentzVector>  FatJet_R1p4_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.4, 30, 5);
        vector<TLorentzVector>  FatJet_R1p5_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.5, 30, 5);
/*
        vector<TLorentzVector>  FatJet_R1p6_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.6, 30, 5);
        vector<TLorentzVector>  FatJet_R1p7_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.7, 30, 5);
        vector<TLorentzVector>  FatJet_R1p8_pT30_Eta5 = makeFatJet(FatJetConstituent, 1.8, 30, 5);
*/
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
/*
        for(unsigned int ifj=0; ifj<FatJet_R1p6_pT30_Eta5.size();ifj++)
        {
            mj_R1p6_pT30_Eta5_.push_back(FatJet_R1p6_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p6_pT30_Eta5_.push_back(FatJet_R1p6_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p6_pT30_Eta5_.push_back(FatJet_R1p6_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p6_pT30_Eta5_.push_back(FatJet_R1p6_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p7_pT30_Eta5.size();ifj++)
        {
            mj_R1p7_pT30_Eta5_.push_back(FatJet_R1p7_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p7_pT30_Eta5_.push_back(FatJet_R1p7_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p7_pT30_Eta5_.push_back(FatJet_R1p7_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p7_pT30_Eta5_.push_back(FatJet_R1p7_pT30_Eta5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p8_pT30_Eta5.size();ifj++)
        {
            mj_R1p8_pT30_Eta5_.push_back(FatJet_R1p8_pT30_Eta5.at(ifj).M());
            FatjetPt_R1p8_pT30_Eta5_.push_back(FatJet_R1p8_pT30_Eta5.at(ifj).Pt());
            FatjetEta_R1p8_pT30_Eta5_.push_back(FatJet_R1p8_pT30_Eta5.at(ifj).Eta());
            FatjetPhi_R1p8_pT30_Eta5_.push_back(FatJet_R1p8_pT30_Eta5.at(ifj).Phi());
        }
*/
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
        
        //  R=1.2  pT(SJ)>35,40,45,50,50  |eta|<2.5
        vector<TLorentzVector>  FatJet_R1p2_pT35_Eta2p5 = makeFatJet(FatJetConstituent, 1.2, 35, 2.5);
        vector<TLorentzVector>  FatJet_R1p2_pT40_Eta2p5 = makeFatJet(FatJetConstituent, 1.2, 40, 2.5);
        vector<TLorentzVector>  FatJet_R1p2_pT45_Eta2p5 = makeFatJet(FatJetConstituent, 1.2, 45, 2.5);
        vector<TLorentzVector>  FatJet_R1p2_pT50_Eta2p5 = makeFatJet(FatJetConstituent, 1.2, 50, 2.5);
        vector<TLorentzVector>  FatJet_R1p2_pT60_Eta2p5 = makeFatJet(FatJetConstituent, 1.2, 60, 2.5);
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT35_Eta2p5.size();ifj++)
        {
            mj_R1p2_pT35_Eta2p5_.push_back(FatJet_R1p2_pT35_Eta2p5.at(ifj).M());
            FatjetPt_R1p2_pT35_Eta2p5_.push_back(FatJet_R1p2_pT35_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p2_pT35_Eta2p5_.push_back(FatJet_R1p2_pT35_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p2_pT35_Eta2p5_.push_back(FatJet_R1p2_pT35_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT40_Eta2p5.size();ifj++)
        {
            mj_R1p2_pT40_Eta2p5_.push_back(FatJet_R1p2_pT40_Eta2p5.at(ifj).M());
            FatjetPt_R1p2_pT40_Eta2p5_.push_back(FatJet_R1p2_pT40_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p2_pT40_Eta2p5_.push_back(FatJet_R1p2_pT40_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p2_pT40_Eta2p5_.push_back(FatJet_R1p2_pT40_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT45_Eta2p5.size();ifj++)
        {
            mj_R1p2_pT45_Eta2p5_.push_back(FatJet_R1p2_pT45_Eta2p5.at(ifj).M());
            FatjetPt_R1p2_pT45_Eta2p5_.push_back(FatJet_R1p2_pT45_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p2_pT45_Eta2p5_.push_back(FatJet_R1p2_pT45_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p2_pT45_Eta2p5_.push_back(FatJet_R1p2_pT45_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT50_Eta2p5.size();ifj++)
        {
            mj_R1p2_pT50_Eta2p5_.push_back(FatJet_R1p2_pT50_Eta2p5.at(ifj).M());
            FatjetPt_R1p2_pT50_Eta2p5_.push_back(FatJet_R1p2_pT50_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p2_pT50_Eta2p5_.push_back(FatJet_R1p2_pT50_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p2_pT50_Eta2p5_.push_back(FatJet_R1p2_pT50_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT60_Eta2p5.size();ifj++)
        {
            mj_R1p2_pT60_Eta2p5_.push_back(FatJet_R1p2_pT60_Eta2p5.at(ifj).M());
            FatjetPt_R1p2_pT60_Eta2p5_.push_back(FatJet_R1p2_pT60_Eta2p5.at(ifj).Pt());
            FatjetEta_R1p2_pT60_Eta2p5_.push_back(FatJet_R1p2_pT60_Eta2p5.at(ifj).Eta());
            FatjetPhi_R1p2_pT60_Eta2p5_.push_back(FatJet_R1p2_pT60_Eta2p5.at(ifj).Phi());
        }
     
        //
        // Cluster GenJet to make FJ 
        //
        vector<TLorentzVector> GenFatJetConstituent; 
        for(int imcjet=0; imcjet<mc_jets_pt->size(); imcjet++) 
        {
            float px_tmp = mc_jets_pt->at(imcjet)*TMath::Cos(mc_jets_phi->at(imcjet));; 
            float py_tmp = mc_jets_pt->at(imcjet)*TMath::Sin(mc_jets_phi->at(imcjet));; 
            float pz_tmp = mc_jets_pt->at(imcjet)/TMath::Tan(2*TMath::ATan(TMath::Exp(-mc_jets_eta->at(imcjet)))); 

            TLorentzVector tmp(px_tmp, py_tmp, pz_tmp, mc_jets_energy->at(imcjet));
            GenFatJetConstituent.push_back(tmp);
        }
        vector<TLorentzVector>  GenFatJet_R0p8_pT30_Eta2p5 = makeFatJet(GenFatJetConstituent, 0.8, 30, 2.5);
        vector<TLorentzVector>  GenFatJet_R1p0_pT30_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.0, 30, 2.5);
        vector<TLorentzVector>  GenFatJet_R1p2_pT30_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.2, 30, 2.5);
        vector<TLorentzVector>  GenFatJet_R1p4_pT30_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.4, 30, 2.5);
        
        vector<TLorentzVector>  GenFatJet_R1p2_pT35_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.2, 35, 2.5);
        vector<TLorentzVector>  GenFatJet_R1p2_pT40_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.2, 40, 2.5);
        vector<TLorentzVector>  GenFatJet_R1p2_pT45_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.2, 45, 2.5);
        vector<TLorentzVector>  GenFatJet_R1p2_pT50_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.2, 50, 2.5);
        vector<TLorentzVector>  GenFatJet_R1p2_pT60_Eta2p5 = makeFatJet(GenFatJetConstituent, 1.2, 60, 2.5);
        
        for(unsigned int ifj=0; ifj<GenFatJet_R0p8_pT30_Eta2p5.size();ifj++)
        {
            Genmj_R0p8_pT30_Eta2p5_.push_back(GenFatJet_R0p8_pT30_Eta2p5.at(ifj).M());
            GenFatjetPt_R0p8_pT30_Eta2p5_.push_back(GenFatJet_R0p8_pT30_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R0p8_pT30_Eta2p5_.push_back(GenFatJet_R0p8_pT30_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R0p8_pT30_Eta2p5_.push_back(GenFatJet_R0p8_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p0_pT30_Eta2p5.size();ifj++)
        {
            Genmj_R1p0_pT30_Eta2p5_.push_back(GenFatJet_R1p0_pT30_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p0_pT30_Eta2p5_.push_back(GenFatJet_R1p0_pT30_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p0_pT30_Eta2p5_.push_back(GenFatJet_R1p0_pT30_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p0_pT30_Eta2p5_.push_back(GenFatJet_R1p0_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p2_pT30_Eta2p5.size();ifj++)
        {
            Genmj_R1p2_pT30_Eta2p5_.push_back(GenFatJet_R1p2_pT30_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p2_pT30_Eta2p5_.push_back(GenFatJet_R1p2_pT30_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p2_pT30_Eta2p5_.push_back(GenFatJet_R1p2_pT30_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p2_pT30_Eta2p5_.push_back(GenFatJet_R1p2_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p4_pT30_Eta2p5.size();ifj++)
        {
            Genmj_R1p4_pT30_Eta2p5_.push_back(GenFatJet_R1p4_pT30_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p4_pT30_Eta2p5_.push_back(GenFatJet_R1p4_pT30_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p4_pT30_Eta2p5_.push_back(GenFatJet_R1p4_pT30_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p4_pT30_Eta2p5_.push_back(GenFatJet_R1p4_pT30_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p2_pT35_Eta2p5.size();ifj++)
        {
            Genmj_R1p2_pT35_Eta2p5_.push_back(GenFatJet_R1p2_pT35_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p2_pT35_Eta2p5_.push_back(GenFatJet_R1p2_pT35_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p2_pT35_Eta2p5_.push_back(GenFatJet_R1p2_pT35_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p2_pT35_Eta2p5_.push_back(GenFatJet_R1p2_pT35_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p2_pT40_Eta2p5.size();ifj++)
        {
            Genmj_R1p2_pT40_Eta2p5_.push_back(GenFatJet_R1p2_pT40_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p2_pT40_Eta2p5_.push_back(GenFatJet_R1p2_pT40_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p2_pT40_Eta2p5_.push_back(GenFatJet_R1p2_pT40_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p2_pT40_Eta2p5_.push_back(GenFatJet_R1p2_pT40_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p2_pT45_Eta2p5.size();ifj++)
        {
            Genmj_R1p2_pT45_Eta2p5_.push_back(GenFatJet_R1p2_pT45_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p2_pT45_Eta2p5_.push_back(GenFatJet_R1p2_pT45_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p2_pT45_Eta2p5_.push_back(GenFatJet_R1p2_pT45_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p2_pT45_Eta2p5_.push_back(GenFatJet_R1p2_pT45_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p2_pT50_Eta2p5.size();ifj++)
        {
            Genmj_R1p2_pT50_Eta2p5_.push_back(GenFatJet_R1p2_pT50_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p2_pT50_Eta2p5_.push_back(GenFatJet_R1p2_pT50_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p2_pT50_Eta2p5_.push_back(GenFatJet_R1p2_pT50_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p2_pT50_Eta2p5_.push_back(GenFatJet_R1p2_pT50_Eta2p5.at(ifj).Phi());
        }
        for(unsigned int ifj=0; ifj<GenFatJet_R1p2_pT60_Eta2p5.size();ifj++)
        {
            Genmj_R1p2_pT60_Eta2p5_.push_back(GenFatJet_R1p2_pT60_Eta2p5.at(ifj).M());
            GenFatjetPt_R1p2_pT60_Eta2p5_.push_back(GenFatJet_R1p2_pT60_Eta2p5.at(ifj).Pt());
            GenFatjetEta_R1p2_pT60_Eta2p5_.push_back(GenFatJet_R1p2_pT60_Eta2p5.at(ifj).Eta());
            GenFatjetPhi_R1p2_pT60_Eta2p5_.push_back(GenFatJet_R1p2_pT60_Eta2p5.at(ifj).Phi());
        }
        // Fat jets on-the-fly ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 


        // FIXME 
        //cout << NIsoTkVeto() << endl; // need to find correct branch names in cfa 


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
        EventWeightNeg_     =   EventWeightNeg;
        MJ_pT30_            =   MJ_pT30;
        MET_                =   pfType1mets_et->at(0);
        METPhi_             =   pfType1mets_phi->at(0);
        HT_                 =   HT;
        Npv_                =   Npv;
        for(unsigned int bc(0); bc<PU_bunchCrossing->size(); ++bc)
        {
            if(PU_bunchCrossing->at(bc)==-1) 
                Npuminusone_ = PU_NumInteractions->at(bc);
            if(PU_bunchCrossing->at(bc)==1) 
                Npuplusone_  = PU_NumInteractions->at(bc);
            if(PU_bunchCrossing->at(bc)==0)
                Npu_ = PU_NumInteractions->at(bc);
        }
        for(int ifj=0; ifj<(int)fjets30_m->size(); ifj++) 
        { 
            mj_pT30_.push_back(fjets30_m->at(ifj));
            FatjetPt_pT30_.push_back(fjets30_pt->at(ifj));
            FatjetEta_pT30_.push_back(fjets30_eta->at(ifj));
            FatjetPhi_pT30_.push_back(fjets30_phi->at(ifj));
            //FatjetN_pT30_.push_back(fastjets_AK4_R1p2_R0p5pT30_nconstituents->at(Vector_GoodFatjet_pT30_Index.at(igoodfatjet)));
        }
        // RelIso
        for(unsigned int imus=0; imus<RA4Muon.size(); imus++) 
        {
            //cout << "[DEBUG] [Muon] mini-isolation : " << imus << " :: " << GetIsolation(imus, 13, 0.2, false, true, true, true, false) << endl;
            RA4MusPt_.push_back(mus_pt->at(RA4Muon.at(imus)));
            RA4MusEta_.push_back(mus_eta->at(RA4Muon.at(imus)));
            RA4MusPhi_.push_back(mus_phi->at(RA4Muon.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4Elec.size(); iels++) 
        {
            //cout << "[DEBUG] [Electron] mini-isolation : " << iels << " :: " << GetIsolation(iels, 11, 0.2, false, true, true, true, false) << endl;
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
          JetE_.push_back(jets_AK4_energy->at(ijet)); 
          JetPt_.push_back(jets_AK4_pt->at(ijet)); 
          JetEta_.push_back(jets_AK4_eta->at(ijet)); 
          JetPhi_.push_back(jets_AK4_phi->at(ijet)); 
          JetCSV_.push_back(jets_AK4_btag_inc_secVertexCombined->at(ijet)); 
        }
        // MiniIso
        for(unsigned int imus=0; imus<RA4Muon_mini.size(); imus++) 
        {
            //cout << "[DEBUG] [Muon] mini-isolation : " << imus << " :: " << GetIsolation(imus, 13, 0.2, false, true, true, true, false) << endl;
            RA4MusPt_mi_.push_back(mus_pt->at(RA4Muon_mini.at(imus)));
            RA4MusEta_mi_.push_back(mus_eta->at(RA4Muon_mini.at(imus)));
            RA4MusPhi_mi_.push_back(mus_phi->at(RA4Muon_mini.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4Elec_mini.size(); iels++) 
        {
            //cout << "[DEBUG] [Electron] mini-isolation : " << iels << " :: " << GetIsolation(iels, 11, 0.2, false, true, true, true, false) << endl;
            RA4ElsPt_mi_.push_back(els_pt->at(RA4Elec_mini.at(iels)));
            RA4ElsEta_mi_.push_back(els_eta->at(RA4Elec_mini.at(iels)));
            RA4ElsPhi_mi_.push_back(els_phi->at(RA4Elec_mini.at(iels)));
        }
        for(unsigned int imus=0; imus<RA4MuonVeto_mini.size(); imus++) 
        {
            RA4MusVetoPt_mi_.push_back(mus_pt->at(RA4MuonVeto_mini.at(imus)));
            RA4MusVetoEta_mi_.push_back(mus_eta->at(RA4MuonVeto_mini.at(imus)));
            RA4MusVetoPhi_mi_.push_back(mus_phi->at(RA4MuonVeto_mini.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4ElecVeto_mini.size(); iels++) 
        {
            RA4ElsVetoPt_mi_.push_back(els_pt->at(RA4ElecVeto_mini.at(iels)));
            RA4ElsVetoEta_mi_.push_back(els_eta->at(RA4ElecVeto_mini.at(iels)));
            RA4ElsVetoPhi_mi_.push_back(els_phi->at(RA4ElecVeto_mini.at(iels)));
        }
        for(unsigned int iisotrkels=0; iisotrkels<eCands.size(); iisotrkels++) 
        {
            if(eCands.at(iisotrkels).second>0.2) continue;
            unsigned int ipfcand = eCands.at(iisotrkels).first;
            IsoTrkVetoElsPt_.push_back(pfcand_pt->at(ipfcand));
            IsoTrkVetoElsEta_.push_back(pfcand_eta->at(ipfcand));
            IsoTrkVetoElsPhi_.push_back(pfcand_phi->at(ipfcand));
        }
        for(unsigned int iisotrkmus=0; iisotrkmus<muCands.size(); iisotrkmus++) 
        {
            if(muCands.at(iisotrkmus).second>0.2) continue;
            unsigned int ipfcand = muCands.at(iisotrkmus).first;
            IsoTrkVetoMusPt_.push_back(pfcand_pt->at(ipfcand));
            IsoTrkVetoMusEta_.push_back(pfcand_eta->at(ipfcand));
            IsoTrkVetoMusPhi_.push_back(pfcand_phi->at(ipfcand));
        }
        for(unsigned int iisotrkhad=0; iisotrkhad<hadCands.size(); iisotrkhad++) 
        {
            if(hadCands.at(iisotrkhad).second>0.1) continue;
            unsigned int ipfcand = hadCands.at(iisotrkhad).first;
            IsoTrkVetoHadPt_.push_back(pfcand_pt->at(ipfcand));
            IsoTrkVetoHadEta_.push_back(pfcand_eta->at(ipfcand));
            IsoTrkVetoHadPhi_.push_back(pfcand_phi->at(ipfcand));
        }

        for(int igoodjet=0; igoodjet<(int)GoodJets_AK4_mini.size(); igoodjet++) 
        {
          int ijet = GoodJets_AK4_mini.at(igoodjet); 
          JetE_mi_.push_back(jets_AK4_energy->at(ijet)); 
          JetPt_mi_.push_back(jets_AK4_pt->at(ijet)); 
          JetEta_mi_.push_back(jets_AK4_eta->at(ijet)); 
          JetPhi_mi_.push_back(jets_AK4_phi->at(ijet)); 
          JetCSV_mi_.push_back(jets_AK4_btag_inc_secVertexCombined->at(ijet)); 
        }

        // gen top pT for TTbar samples  
        for(int igen=0; igen<(int)mc_doc_id->size(); igen++) 
        {   
            // ME ISR 
            //  status=23 && same parents as top with status=22 or proton(2212)
            // PS ISR 
            //  status!=23 && same parents as top with status=22 or proton(2212)
            // FSR 
            //  parent is top with status=22 
            GenPt_.push_back(   mc_doc_pt->at(igen));  
            GenEta_.push_back(  mc_doc_eta->at(igen));  
            GenPhi_.push_back(  mc_doc_phi->at(igen));  
            GenStatus_.push_back(  mc_doc_status->at(igen));  
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
        for(int igenjet=0; igenjet<jets_AK4_gen_pt->size(); igenjet++) 
        { 
            GenJetPt_.push_back(jets_AK4_gen_pt->at(igenjet)); 
            //GenJetEta_.push_back(jets_AK4_gen_eta->at(igenjet)); 
            //GenJetPhi_.push_back(jets_AK4_gen_phi->at(igenjet)); 
        }
        GenMET_                =   pfType1mets_gen_et->at(0);
        for(int imcjet=0; imcjet<mc_jets_pt->size(); imcjet++) 
        { 
            MCJetE_.push_back(mc_jets_energy->at(imcjet)); 
            MCJetPt_.push_back(mc_jets_pt->at(imcjet)); 
            MCJetEta_.push_back(mc_jets_eta->at(imcjet)); 
            MCJetPhi_.push_back(mc_jets_phi->at(imcjet)); 
        }

        babyTree_->Fill(); // Fill all events

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
    delete chainA;
    delete chainB;
}

