#include <vector>
#include <string>
#include "TChain.h"

using namespace std;

// event info  
Long64_t        event_;
int             run_;
int             lumi_;
// weights, filters, trigger 
bool            stitch_; 
float           weight_;
bool            pass_; // need to add more pass_xx 
bool            pass_hbhe_;
bool            pass_cschalo_;
bool            pass_eebadsc_;
bool            pass_goodv_; 
vector<bool>    *trig_;
// PU
int             ntrupv_;
float           ntrupv_mean_;
int             npv_;
// event variables 
float           ht_;
float           met_;
float           met_phi_;
float           mt_;
// fat jets 
float           mj_;
int             nfjets_;
vector<float>   *fjets_m_;
vector<float>   *fjets_pt_;
vector<float>   *fjets_phi_;
vector<float>   *fjets_eta_;
float           mj08_;
int             nfjets08_;
vector<float>   *fjets08_m_;
vector<float>   *fjets08_pt_;
vector<float>   *fjets08_phi_; 
vector<float>   *fjets08_eta_;
// leptons
int             nleps_; 
int             nvleps_; 
vector<float>   *leps_pt_;
vector<float>   *leps_phi_;
vector<float>   *leps_eta_;
vector<float>   *leps_id_;
int             nmus_;
int             nvmus_;
vector<float>   *mus_pt_;
vector<float>   *mus_phi_;
vector<float>   *mus_eta_;
vector<float>   *mus_charge_;
vector<float>   *mus_miniso_;
vector<bool>    *mus_sigid_;
int             nels_;
int             nvels_;
vector<float>   *els_pt_;
vector<float>   *els_phi_;
vector<float>   *els_eta_;
vector<float>   *els_charge_;
vector<float>   *els_miniso_;
vector<bool>    *els_sigid_;
// jets
int             njets_;
int             nbm_;
vector<float>   *jets_m_;
vector<float>   *jets_pt_;
vector<float>   *jets_phi_;
vector<float>   *jets_eta_;
vector<float>   *jets_csv_;
vector<bool>    *jets_islep_;
// mc truth 
int ntruleps_;
int ntruels_;
int ntrumus_;
int ntrutaush_;
int ntrutausl_;
vector<int>     *mc_id_;
vector<int>     *mc_mom_;
vector<float>   *mc_pt_;
vector<float>   *mc_phi_;
vector<float>   *mc_eta_;
vector<float>   *mc_mass_;

//////////////////////////////////////////////////////

// event info  
TBranch         *b_event;  
TBranch         *b_run;  
TBranch         *b_lumi;    
// weights, filters, trigger 
TBranch         *b_stitch;
TBranch         *b_weight;  
TBranch         *b_pass;
TBranch         *b_pass_hbhe;
TBranch         *b_pass_cschalo;
TBranch         *b_pass_eebadsc;
TBranch         *b_pass_goodv;
TBranch         *b_trig;   //!
// PU
TBranch         *b_ntrupv;  
TBranch         *b_ntrupv_mean;  
TBranch         *b_npv;   
// event variables 
TBranch         *b_ht;   
TBranch         *b_met;   
TBranch         *b_met_phi;   
TBranch         *b_mt;  
// fat jets
TBranch         *b_mj;  
TBranch         *b_nfjets;   
TBranch         *b_fjets_m;   //!
TBranch         *b_fjets_pt;   //!
TBranch         *b_fjets_phi;   //!
TBranch         *b_fjets_eta;   //!
TBranch         *b_mj08;  
TBranch         *b_nfjets08;   
TBranch         *b_fjets08_m;   //!
TBranch         *b_fjets08_pt;   //!
TBranch         *b_fjets08_phi;   //!
TBranch         *b_fjets08_eta;   //!
// leptons 
TBranch         *b_nleps;   //!
TBranch         *b_nvleps;   //!
TBranch         *b_leps_pt;   //!
TBranch         *b_leps_phi;   //!
TBranch         *b_leps_eta;   //!
TBranch         *b_leps_id;   //!
TBranch         *b_nmus;
TBranch         *b_nvmus;
TBranch         *b_mus_pt;   //!
TBranch         *b_mus_phi;   //!
TBranch         *b_mus_eta;   //!
TBranch         *b_mus_charge;   //!
TBranch         *b_mus_miniso;   //!
TBranch         *b_mus_sigid;   //!
TBranch         *b_nels;
TBranch         *b_nvels;
TBranch         *b_els_pt;   //!
TBranch         *b_els_phi;   //!
TBranch         *b_els_eta;   //!
TBranch         *b_els_charge;   //!
TBranch         *b_els_miniso;   //!
TBranch         *b_els_sigid;   //!
// jets
TBranch         *b_njets;   
TBranch         *b_nbm;   
TBranch         *b_jets_m;   //!
TBranch         *b_jets_pt;   //!
TBranch         *b_jets_phi;   //!
TBranch         *b_jets_eta;   //!
TBranch         *b_jets_csv;   //!
TBranch         *b_jets_islep;   //!
// mc truth
TBranch         *b_ntruleps;
TBranch         *b_ntrumus;
TBranch         *b_ntruels;
TBranch         *b_ntrutausl;
TBranch         *b_ntrutaush;
TBranch         *b_mc_id;   //!
TBranch         *b_mc_mom;   //!
TBranch         *b_mc_pt;   //!
TBranch         *b_mc_phi;   //!
TBranch         *b_mc_eta;   //!
TBranch         *b_mc_mass;   //!

void InitBaby(TChain *ch)
{
    trig_        = 0;
    mj_          = 0;
    fjets_m_     = 0;
    fjets_pt_    = 0;
    fjets_phi_   = 0;
    fjets_eta_   = 0;
    fjets08_m_   = 0;
    fjets08_pt_  = 0;
    fjets08_phi_ = 0;
    fjets08_eta_ = 0;
    leps_pt_     = 0;
    leps_phi_    = 0;
    leps_eta_    = 0;
    leps_id_     = 0;
    mus_pt_      = 0;
    mus_phi_     = 0;
    mus_eta_     = 0;
    mus_charge_  = 0;
    els_pt_      = 0;
    els_phi_     = 0;
    els_eta_     = 0;
    els_charge_  = 0;
    jets_m_      = 0;
    jets_pt_     = 0;
    jets_phi_    = 0;
    jets_eta_    = 0;
    jets_csv_    = 0;
    jets_islep_  = 0;
    mc_id_       = 0;
    mc_mom_      = 0;
    mc_pt_       = 0;
    mc_phi_      = 0;
    mc_eta_      = 0;
    mc_mass_     = 0;

    // event info  
    ch->SetBranchAddress("event",           &event_,        &b_event);
    ch->SetBranchAddress("run",             &run_,          &b_run);
    ch->SetBranchAddress("lumiblock",       &lumi_,         &b_lumi); 
    // weights, filters, trigger 
    ch->SetBranchAddress("weight",          &weight_,       &b_weight);
    ch->SetBranchAddress("stitch",          &stitch_,       &b_stitch);
    ch->SetBranchAddress("pass",            &pass_,         &b_pass);
    ch->SetBranchAddress("pass_hbhe",       &pass_hbhe_,    &b_pass_hbhe);
    ch->SetBranchAddress("pass_cschalo",    &pass_cschalo_, &b_pass_cschalo);
    ch->SetBranchAddress("pass_eebadsc",    &pass_eebadsc_, &b_pass_eebadsc);
    ch->SetBranchAddress("pass_goodv",      &pass_goodv_,   &b_pass_goodv);
    ch->SetBranchAddress("trig",            &trig_,         &b_trig);
    // PU 
    ch->SetBranchAddress("ntrupv",          &ntrupv_,       &b_ntrupv);
    ch->SetBranchAddress("ntrupv_mean",     &ntrupv_mean_,  &b_ntrupv_mean);
    ch->SetBranchAddress("npv",             &npv_,          &b_npv);
    // event variables 
    ch->SetBranchAddress("ht",              &ht_,           &b_ht);
    ch->SetBranchAddress("met",             &met_,          &b_met);
    ch->SetBranchAddress("met_phi",         &met_phi_,      &b_met_phi); 
    ch->SetBranchAddress("mt",              &mt_,           &b_mt);
    // fat jets
    ch->SetBranchAddress("mj",              &mj_,           &b_mj);
    ch->SetBranchAddress("nfjets",          &nfjets_,       &b_nfjets);
    ch->SetBranchAddress("fjets_m",         &fjets_m_,      &b_fjets_m);
    ch->SetBranchAddress("fjets_pt",        &fjets_pt_,     &b_fjets_pt);
    ch->SetBranchAddress("fjets_phi",       &fjets_phi_,    &b_fjets_phi);
    ch->SetBranchAddress("fjets_eta",       &fjets_eta_,    &b_fjets_eta);
    ch->SetBranchAddress("mj08",            &mj08_,         &b_mj08);
    ch->SetBranchAddress("nfjets08",        &nfjets08_,     &b_nfjets08);
    ch->SetBranchAddress("fjets08_m",       &fjets08_m_,    &b_fjets08_m);
    ch->SetBranchAddress("fjets08_pt",      &fjets08_pt_,   &b_fjets08_pt);
    ch->SetBranchAddress("fjets08_phi",     &fjets08_phi_,  &b_fjets08_phi);
    ch->SetBranchAddress("fjets08_eta",     &fjets08_eta_,  &b_fjets08_eta);
    // leptons 
    ch->SetBranchAddress("nleps",           &nleps_,        &b_nleps);
    ch->SetBranchAddress("nvleps",          &nvleps_,       &b_nvleps);
    ch->SetBranchAddress("leps_pt",         &leps_pt_,      &b_leps_pt);
    ch->SetBranchAddress("leps_phi",        &leps_phi_,     &b_leps_phi);
    ch->SetBranchAddress("leps_eta",        &leps_eta_,     &b_leps_eta);
    ch->SetBranchAddress("leps_id",         &leps_id_,      &b_leps_id);
    ch->SetBranchAddress("nmus",            &nmus_,         &b_nmus);
    ch->SetBranchAddress("nvmus",           &nvmus_,        &b_nvmus);
    ch->SetBranchAddress("mus_pt",          &mus_pt_,       &b_mus_pt);
    ch->SetBranchAddress("mus_phi",         &mus_phi_,      &b_mus_phi);
    ch->SetBranchAddress("mus_eta",         &mus_eta_,      &b_mus_eta);
    ch->SetBranchAddress("mus_charge",      &mus_charge_,   &b_mus_charge);
    ch->SetBranchAddress("mus_miniso",      &mus_miniso_,   &b_mus_miniso);
    ch->SetBranchAddress("mus_sigid",       &mus_sigid_,    &b_mus_sigid);
    ch->SetBranchAddress("nels",            &nels_,         &b_nels);
    ch->SetBranchAddress("nvels",           &nvels_,        &b_nvels);
    ch->SetBranchAddress("els_pt",          &els_pt_,       &b_els_pt);
    ch->SetBranchAddress("els_phi",         &els_phi_,      &b_els_phi);
    ch->SetBranchAddress("els_eta",         &els_eta_,      &b_els_eta);
    ch->SetBranchAddress("els_charge",      &els_charge_,   &b_els_charge);
    ch->SetBranchAddress("els_miniso",      &els_miniso_,   &b_els_miniso);
    ch->SetBranchAddress("els_sigid",       &els_sigid_,    &b_els_sigid);
    // jets
    ch->SetBranchAddress("njets",           &njets_,        &b_njets);
    ch->SetBranchAddress("nbm",             &nbm_,          &b_nbm);
    ch->SetBranchAddress("jets_m",          &jets_m_,       &b_jets_m);
    ch->SetBranchAddress("jets_pt",         &jets_pt_,      &b_jets_pt);
    ch->SetBranchAddress("jets_phi",        &jets_phi_,     &b_jets_phi);
    ch->SetBranchAddress("jets_eta",        &jets_eta_,     &b_jets_eta);
    ch->SetBranchAddress("jets_csv",        &jets_csv_,     &b_jets_csv);
    ch->SetBranchAddress("jets_islep",      &jets_islep_,   &b_jets_islep);
    // mc truth 
    ch->SetBranchAddress("ntruleps",        &ntruleps_,     &b_ntruleps);
    ch->SetBranchAddress("ntrumus",         &ntrumus_,      &b_ntrumus);
    ch->SetBranchAddress("ntruels",         &ntruels_,      &b_ntruels);
    ch->SetBranchAddress("ntrutausl",       &ntrutausl_,    &b_ntrutausl);
    ch->SetBranchAddress("ntrutaush",       &ntrutaush_,    &b_ntrutaush);
    ch->SetBranchAddress("mc_id",           &mc_id_,        &b_mc_id);
    ch->SetBranchAddress("mc_mom",          &mc_mom_,       &b_mc_mom);
    ch->SetBranchAddress("mc_pt",           &mc_pt_,        &b_mc_pt);
    ch->SetBranchAddress("mc_phi",          &mc_phi_,       &b_mc_phi);
    ch->SetBranchAddress("mc_eta",          &mc_eta_,       &b_mc_eta);
    ch->SetBranchAddress("mc_mass",         &mc_mass_,      &b_mc_mass);
   
}
