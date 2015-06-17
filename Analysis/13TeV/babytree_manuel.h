#include <vector>
#include <string>
#include "TChain.h"

using namespace std;

// bool status = true;
int             event_;
float           EventWeight_;
float           EventWeightNeg_;
int             Npu_;
int             Npv_;
int             Nfatjet_;
int             Nskinnyjet_;
int             NBtagCSVM_;
float           MJ_;
float           MET_;
float           HT_;
float           METPhi_;
vector<float>   *mj_;
vector<float>   *mj08_;
vector<float>   *FatjetPt_;
vector<float>   *FatjetEta_;
vector<float>   *FatjetPhi_;
vector<float>   *RA4MusPt_;
vector<float>   *RA4MusPhi_;
vector<float>   *RA4MusEta_;
vector<float>   *RA4MusQ_;
vector<float>   *RA4ElsPt_;
vector<float>   *RA4ElsPhi_;
vector<float>   *RA4ElsEta_;
vector<float>   *RA4ElsQ_;
vector<float>   *JetPt_;
vector<float>   *JetEta_;
vector<float>   *JetPhi_;
vector<float>   *JetCSV_;
vector<bool>    *JetIsLep_;
vector<float>   *RA4MusVetoPt_;
vector<float>   *RA4ElsVetoPt_;
vector<int>     *GenId_;
vector<float>   *GenStatus_;
vector<int>     *GenMId_;
vector<int>     *GenGMId_;
vector<float>   *GenPt_;
vector<float>   *GenEta_;
vector<float>   *GenPhi_;

int nels_;
int nmus_;
int nvels_;
int nvmus_;
float  mt_;
unsigned int  mc_type_;
int ntruleps_;
int ntruels_;
int ntrumus_;
int ntrutaush_;
int ntrutausl_;

//vector<float>   *GenE_;

bool        TrigMuon_;
bool        TrigSingleMuon_;
float        top1pT_;
float        top1Phi_;
float        top2pT_;
float        top2Phi_;

TBranch         *b_event;   //!
TBranch         *b_EventWeight;   //!
TBranch         *b_EventWeightNeg;   //!
TBranch         *b_Npu;   //!
TBranch         *b_Npv;   //!
TBranch         *b_Nfatjet;   //!
TBranch         *b_Nskinnyjet;   //!
TBranch         *b_NBtagCSVM;   //!
TBranch         *b_MJ;   //!
TBranch         *b_MET;   //!
TBranch         *b_HT;   //!
TBranch         *b_METPhi;   //!
TBranch         *b_mj;   //!
TBranch         *b_mj08;   //!
TBranch         *b_FatjetPt;   //!
TBranch         *b_FatjetEta;   //!
TBranch         *b_FatjetPhi;   //!
TBranch         *b_mus_pt;   //!
TBranch         *b_mus_phi;   //!
TBranch         *b_mus_eta;   //!
TBranch         *b_mus_charge;   //!
TBranch         *b_els_pt;   //!
TBranch         *b_els_phi;   //!
TBranch         *b_els_eta;   //!
TBranch         *b_els_charge;   //!
TBranch         *b_JetPt;   //!
TBranch         *b_JetEta;   //!
TBranch         *b_JetPhi;   //!
TBranch         *b_JetCSV;   //!
TBranch         *b_JetIsLep;   //!
TBranch         *b_RA4MusVetoPt;   //!
TBranch         *b_RA4ElsVetoPt;   //!
TBranch         *b_GenId;   //!
TBranch         *b_GenStatus;   //!
TBranch         *b_GenPt;   //!
TBranch         *b_GenPhi;   //!
TBranch         *b_GenEta;   //!
//TBranch         *b_GenE;   //!
TBranch         *b_GenMId;   //!
TBranch         *b_GenGMId;   //!
TBranch     *b_TrigMuon; 
TBranch     *b_TrigSingleMuon; 
TBranch     *b_top1pT; 
TBranch     *b_top1Phi; 
TBranch     *b_top2pT; 
TBranch     *b_top2Phi;
TBranch *b_nels;
TBranch *b_nmus;
TBranch *b_nvels;
TBranch *b_nvmus;
TBranch * b_mt;
TBranch * b_mc_type;
TBranch *b_ntruleps;
TBranch *b_ntrumus;
TBranch *b_ntruels;
TBranch *b_ntrutausl;
TBranch *b_ntrutaush;

void InitBaby(TChain *ch)
{
    mj_        = 0;
    mj08_      = 0;
    FatjetPt_  = 0;
    FatjetEta_ = 0;
    FatjetPhi_ = 0;
    RA4MusPt_  = 0;
    RA4MusPhi_ = 0;
    RA4MusEta_ = 0;
    RA4MusQ_  = 0;
    RA4ElsPt_  = 0;
    RA4ElsPhi_  = 0;
    RA4ElsEta_  = 0;
    RA4ElsQ_  = 0;
    JetPt_  = 0;
    JetEta_ = 0;
    JetPhi_ = 0;
    JetCSV_ = 0;
    JetIsLep_ = 0;
    RA4MusVetoPt_  = 0;
    RA4ElsVetoPt_  = 0;
    GenId_  = 0;
    GenPt_  = 0;
    GenPhi_  = 0;
    GenEta_  = 0;
    //GenE_  = 0;
    GenStatus_  = 0;
    GenMId_  = 0;
    GenGMId_  = 0;


    ch->SetBranchAddress("event",           &event_,         &b_event);
    ch->SetBranchAddress("weight",     &EventWeight_,   &b_EventWeight);
    ch->SetBranchAddress("ntrupv",             &Npu_,           &b_Npu);
    ch->SetBranchAddress("npv",             &Npv_,           &b_Npv);
    ch->SetBranchAddress("nfjets",    &Nfatjet_,       &b_Nfatjet);
    ch->SetBranchAddress("njets",      &Nskinnyjet_,    &b_Nskinnyjet);
    ch->SetBranchAddress("nbm",       &NBtagCSVM_,     &b_NBtagCSVM);
    ch->SetBranchAddress("mj",         &MJ_,            &b_MJ);
    ch->SetBranchAddress("met",             &MET_,           &b_MET);
    ch->SetBranchAddress("ht",              &HT_,            &b_HT);
    ch->SetBranchAddress("met_phi",          &METPhi_,        &b_METPhi);
    ch->SetBranchAddress("fjets_m",         &mj_,            &b_mj);
    ch->SetBranchAddress("fjets_pt",   &FatjetPt_,      &b_FatjetPt);
    ch->SetBranchAddress("fjets_eta",  &FatjetEta_,     &b_FatjetEta);
    ch->SetBranchAddress("fjets_phi",  &FatjetPhi_,     &b_FatjetPhi);
    ch->SetBranchAddress("fjets08_m",         &mj08_,            &b_mj08);
//    ch->SetBranchAddress("fjets08_pt",   &FatjetPt_,      &b_FatjetPt);
//    ch->SetBranchAddress("fjets08_eta",  &FatjetEta_,     &b_FatjetEta);
//    ch->SetBranchAddress("fjets08_phi",  &FatjetPhi_,     &b_FatjetPhi);

    ch->SetBranchAddress("nels",        &nels_,      &b_nels);
    ch->SetBranchAddress("nmus",        &nmus_,      &b_nmus);
    ch->SetBranchAddress("nvels",        &nvels_,      &b_nvels);
    ch->SetBranchAddress("nvmus",        &nvmus_,      &b_nvmus);
    ch->SetBranchAddress("mt",        &mt_,      &b_mt);
//    ch->SetBranchAddress("mc_type",        &mc_type_,      &b_mc_type);
    ch->SetBranchAddress("ntruleps",        &ntruleps_,      &b_ntruleps);
    ch->SetBranchAddress("ntrumus",        &ntrumus_,      &b_ntrumus);
    ch->SetBranchAddress("ntruels",        &ntruels_,      &b_ntruels);
    ch->SetBranchAddress("ntrutausl",        &ntrutausl_,      &b_ntrutausl);
    ch->SetBranchAddress("ntrutaush",        &ntrutaush_,      &b_ntrutaush);

    ch->SetBranchAddress("mus_pt",          &RA4MusPt_,         &b_mus_pt);
    ch->SetBranchAddress("mus_phi",         &RA4MusPhi_,        &b_mus_phi);
    ch->SetBranchAddress("mus_eta",         &RA4MusEta_,        &b_mus_eta);
    ch->SetBranchAddress("mus_charge",      &RA4MusQ_,          &b_mus_charge);
    ch->SetBranchAddress("els_pt",          &RA4ElsPt_,         &b_els_pt);
    ch->SetBranchAddress("els_phi",         &RA4ElsPhi_,        &b_els_phi);
    ch->SetBranchAddress("els_eta",         &RA4ElsEta_,        &b_els_eta);
    ch->SetBranchAddress("els_charge",      &RA4ElsQ_,          &b_els_charge);

    ch->SetBranchAddress("jets_pt",           &JetPt_,         &b_JetPt);
    ch->SetBranchAddress("jets_phi",          &JetPhi_,        &b_JetPhi);
    ch->SetBranchAddress("jets_eta",          &JetEta_,        &b_JetEta);
    ch->SetBranchAddress("jets_csv",          &JetCSV_,        &b_JetCSV);
    ch->SetBranchAddress("jets_islep",          &JetIsLep_,        &b_JetIsLep);
    /* ch->SetBranchAddress("RA4MusVetoPt_mi",    &RA4MusVetoPt_, &b_RA4MusVetoPt);
       ch->SetBranchAddress("RA4ElsVetoPt_mi",    &RA4ElsVetoPt_, &b_RA4ElsVetoPt);*/

    //ch->SetBranchAddress("mc_id",           &GenId_,         &b_GenId);
    //ch->SetBranchAddress("mc_pt",           &GenPt_,         &b_GenPt);
    //ch->SetBranchAddress("mc_phi",           &GenPhi_,         &b_GenPhi);
    //ch->SetBranchAddress("mc_eta",           &GenEta_,         &b_GenEta);
    // ch->SetBranchAddress("GenE",           &GenE_,         &b_GenE);
    //ch->SetBranchAddress("mc_mom",          &GenMId_,        &b_GenMId);
    //if(status)ch->SetBranchAddress("mc_status",          &GenStatus_,        &b_GenStatus);

    //ch->SetBranchAddress("mc_gmom",         &GenGMId_,       &b_GenGMId);
    /* ch->SetBranchAddress("top1pT",        &top1pT_, &b_top1pT);
       ch->SetBranchAddress("top1Phi",       &top1Phi_, &b_top1Phi);
       ch->SetBranchAddress("top2pT",        &top2pT_, &b_top2pT);
       ch->SetBranchAddress("top2Phi",       &top2Phi_, &b_top2Phi);*/
    //ch->SetBranchAddress("TrigMuon",      &TrigMuon_, &b_TrigMuon);
    //ch->SetBranchAddress("TrigSingleMuon", &TrigSingleMuon_, &b_TrigSingleMuon);
}
