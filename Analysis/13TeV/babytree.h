#include <vector>
#include <string>
#include "TChain.h"

    int             event;
    float           EventWeight;
    float           Npu;
    int             Npv;
    int             Nfatjet;
    int             Nskinnyjet;
    int             NBtagCSVM;
    float           MJ;
    float           MET;
    float           HT;
    float           METPhi;
    vector<float>   *mj;
    vector<float>   *FatjetPt;
    vector<float>   *FatjetEta;
    vector<float>   *FatjetPhi;
    vector<float>   *RA4MusPt;
    vector<float>   *RA4MusPhi;
    vector<float>   *RA4MusEta;
    vector<float>   *RA4ElsPt;
    vector<float>   *RA4ElsPhi;
    vector<float>   *RA4ElsEta;
    vector<float>   *JetPt;
    vector<float>   *JetEta;
    vector<float>   *JetPhi;
    vector<float>   *JetCSV;
    vector<float>   *RA4MusVetoPt;
    vector<float>   *RA4ElsVetoPt;
    vector<float>   *GenId;
    vector<float>   *GenMId;
    vector<float>   *GenGMId;
    bool        TrigMuon;
    bool        TrigSingleMuon;
    float        top1pT;
    float        top1Phi;
    float        top2pT;
    float        top2Phi;
    
    TBranch         *b_event;   //!
    TBranch         *b_EventWeight;   //!
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
    TBranch         *b_FatjetPt;   //!
    TBranch         *b_FatjetEta;   //!
    TBranch         *b_FatjetPhi;   //!
    TBranch         *b_RA4MusPt;   //!
    TBranch         *b_RA4MusPhi;   //!
    TBranch         *b_RA4MusEta;   //!
    TBranch         *b_RA4ElsPt;   //!
    TBranch         *b_RA4ElsPhi;   //!
    TBranch         *b_RA4ElsEta;   //!
    TBranch         *b_JetPt;   //!
    TBranch         *b_JetEta;   //!
    TBranch         *b_JetPhi;   //!
    TBranch         *b_JetCSV;   //!
    TBranch         *b_RA4MusVetoPt;   //!
    TBranch         *b_RA4ElsVetoPt;   //!
    TBranch         *b_GenId;   //!
    TBranch         *b_GenMId;   //!
    TBranch         *b_GenGMId;   //!
    TBranch     *b_TrigMuon; 
    TBranch     *b_TrigSingleMuon; 
    TBranch     *b_top1pT; 
    TBranch     *b_top1Phi; 
    TBranch     *b_top2pT; 
    TBranch     *b_top2Phi; 
    
void InitBaby(TChain *ch)
{
    mj        = 0;
    FatjetPt  = 0;
    FatjetEta = 0;
    FatjetPhi = 0;
    RA4MusPt  = 0;
    RA4MusPhi = 0;
    RA4MusEta = 0;
    RA4ElsPt  = 0;
    RA4ElsPhi  = 0;
    RA4ElsEta  = 0;
    JetPt  = 0;
    JetEta = 0;
    JetPhi = 0;
    JetCSV = 0;
    RA4MusVetoPt  = 0;
    RA4ElsVetoPt  = 0;
    GenId  = 0;
    GenMId  = 0;
    GenGMId  = 0;

    ch->SetBranchAddress("event",       &event,         &b_event);
    ch->SetBranchAddress("EventWeight", &EventWeight,   &b_EventWeight);
    ch->SetBranchAddress("Npu",         &Npu,           &b_Npu);
    ch->SetBranchAddress("Npv",         &Npv,           &b_Npv);
    ch->SetBranchAddress("Nfatjet_pT30",     &Nfatjet,       &b_Nfatjet);
    ch->SetBranchAddress("Nskinnyjet",  &Nskinnyjet,    &b_Nskinnyjet);
    ch->SetBranchAddress("NBtagCSVM",   &NBtagCSVM,     &b_NBtagCSVM);
    ch->SetBranchAddress("MJ_pT30",          &MJ,            &b_MJ);
    ch->SetBranchAddress("MET",         &MET,           &b_MET);
    ch->SetBranchAddress("HT",          &HT,            &b_HT);
    ch->SetBranchAddress("METPhi",      &METPhi,        &b_METPhi);
    ch->SetBranchAddress("mj_pT30",          &mj,            &b_mj);
    ch->SetBranchAddress("FatjetPt_pT30",    &FatjetPt,      &b_FatjetPt);
    ch->SetBranchAddress("FatjetEta_pT30",   &FatjetEta,     &b_FatjetEta);
    ch->SetBranchAddress("FatjetPhi_pT30",   &FatjetPhi,     &b_FatjetPhi);
    ch->SetBranchAddress("RA4MusPt",    &RA4MusPt,      &b_RA4MusPt);
    ch->SetBranchAddress("RA4MusPhi",   &RA4MusPhi,     &b_RA4MusPhi);
    ch->SetBranchAddress("RA4MusEta",   &RA4MusEta,     &b_RA4MusEta);
    ch->SetBranchAddress("RA4ElsPt",    &RA4ElsPt,      &b_RA4ElsPt);
    ch->SetBranchAddress("RA4ElsPhi",   &RA4ElsPhi,     &b_RA4ElsPhi);
    ch->SetBranchAddress("RA4ElsEta",   &RA4ElsEta,     &b_RA4ElsEta);
    ch->SetBranchAddress("JetPt",       &JetPt,         &b_JetPt);
    ch->SetBranchAddress("JetPhi",      &JetPhi,        &b_JetPhi);
    ch->SetBranchAddress("JetEta",      &JetEta,        &b_JetEta);
    ch->SetBranchAddress("JetCSV",      &JetCSV,        &b_JetCSV);
    ch->SetBranchAddress("RA4MusVetoPt", &RA4MusVetoPt, &b_RA4MusVetoPt);
    ch->SetBranchAddress("RA4ElsVetoPt", &RA4ElsVetoPt, &b_RA4ElsVetoPt);
    ch->SetBranchAddress("GenId",       &GenId,         &b_GenId);
    ch->SetBranchAddress("GenMId",      &GenMId,        &b_GenMId);
    ch->SetBranchAddress("GenGMId",     &GenGMId,       &b_GenGMId);
    //ch->SetBranchAddress("top1pT", &top1pT, &b_top1pT);
    //ch->SetBranchAddress("top1Phi", &top1Phi, &b_top1Phi);
    //ch->SetBranchAddress("top2pT", &top2pT, &b_top2pT);
    //ch->SetBranchAddress("top2Phi", &top2Phi, &b_top2Phi);
    //ch->SetBranchAddress("TrigMuon", &TrigMuon, &b_TrigMuon);
    //ch->SetBranchAddress("TrigSingleMuon", &TrigSingleMuon, &b_TrigSingleMuon);
}
