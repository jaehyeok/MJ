//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 23 11:24:43 2014 by ROOT version 5.34/18
// from TTree eventB/StringBasedNTupler tree
// found on file: slim_cfA_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66_f1000_1_aK6.root
//////////////////////////////////////////////////////////

#include <TTree.h>
#include <TChain.h>
#include <vector>
#include <string>

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *els_conversion_dist;
   vector<float>   *els_conversion_dcot;
   vector<float>   *els_PFchargedHadronIsoR03;
   vector<float>   *els_PFphotonIsoR03;
   vector<float>   *els_PFneutralHadronIsoR03;
   vector<bool>    *els_hasMatchedConversion;
   vector<float>   *PU_TrueNumInteractions;
   Float_t         rho_kt6PFJetsForIsolation2011;
   
   vector<float>   *beamSpot_x;
   vector<float>   *beamSpot_y;
   vector<float>   *els_energy;
   vector<float>   *els_et;
   vector<float>   *els_eta;
   vector<float>   *els_phi;
   vector<float>   *els_pt;
   vector<float>   *els_px;
   vector<float>   *els_py;
   vector<float>   *els_pz;
   vector<float>   *els_robustTightId;
   vector<float>   *els_simpleEleId95relIso;
   vector<float>   *els_simpleEleId90relIso;
   vector<float>   *els_simpleEleId85relIso;
   vector<float>   *els_simpleEleId80relIso;
   vector<float>   *els_simpleEleId70relIso;
   vector<float>   *els_simpleEleId95cIso;
   vector<float>   *els_simpleEleId90cIso;
   vector<float>   *els_simpleEleId85cIso;
   vector<float>   *els_simpleEleId80cIso;
   vector<float>   *els_simpleEleId70cIso;
   vector<float>   *els_cIso;
   vector<float>   *els_tIso;
   vector<float>   *els_ecalIso;
   vector<float>   *els_hcalIso;
   vector<float>   *els_charge;
   vector<float>   *els_caloEnergy;
   vector<float>   *els_hadOverEm;
   vector<float>   *els_eOverPIn;
   vector<float>   *els_sigmaIEtaIEta;
   vector<float>   *els_scEnergy;
   vector<float>   *els_scEta;
   vector<float>   *els_scE1x5;
   vector<float>   *els_scE2x5Max;
   vector<float>   *els_scE5x5;
   vector<float>   *els_isEB;
   vector<float>   *els_isEE;
   vector<float>   *els_dEtaIn;
   vector<float>   *els_dPhiIn;
   vector<float>   *els_dEtaOut;
   vector<float>   *els_dPhiOut;
   vector<float>   *els_numlosthits;
   vector<float>   *els_tk_pz;
   vector<float>   *els_tk_pt;
   vector<float>   *els_tk_phi;
   vector<float>   *els_d0dum;
   vector<float>   *els_vx;
   vector<float>   *els_vy;
   vector<float>   *els_vz;
   vector<float>   *els_ptError;
   vector<float>   *els_n_inner_layer;
   vector<float>   *els_dr03EcalRecHitSumEt;
   vector<float>   *els_dr03HcalTowerSumEt;
   vector<float>   *els_dr03HcalDepth1TowerSumEt;
   vector<float>   *els_dr03HcalDepth2TowerSumEt;
   vector<float>   *els_dr03TkSumPt;
   vector<float>   *jets_AK5PF_phi;
   vector<float>   *jets_AK5PF_pt;
   vector<float>   *jets_AK5PF_pz;
   vector<float>   *jets_AK5PF_px;
   vector<float>   *jets_AK5PF_py;
   vector<float>   *jets_AK5PF_eta;
   vector<float>   *jets_AK5PF_et;
   vector<float>   *jets_AK5PF_energy;
   vector<float>   *jets_AK5PF_parton_Id;
   vector<float>   *jets_AK5PF_parton_motherId;
   vector<float>   *jets_AK5PF_gen_pt;
   vector<float>   *jets_AK5PF_gen_phi;
   vector<float>   *jets_AK5PF_partonFlavour;
   vector<float>   *jets_AK5PF_btag_TC_highPur;
   vector<float>   *jets_AK5PF_btag_TC_highEff;
   vector<float>   *jets_AK5PF_btag_jetProb;
   vector<float>   *jets_AK5PF_btag_jetBProb;
   vector<float>   *jets_AK5PF_btag_secVertexHighPur;
   vector<float>   *jets_AK5PF_btag_secVertexHighEff;
   vector<float>   *jets_AK5PF_btag_secVertexCombined;
   vector<float>   *jets_AK5PF_jetCharge;
   vector<float>   *jets_AK5PF_chgEmE;
   vector<float>   *jets_AK5PF_chgHadE;
   vector<float>   *jets_AK5PF_photonEnergy;
   vector<float>   *jets_AK5PF_chg_Mult;
   vector<float>   *jets_AK5PF_neutralEmE;
   vector<float>   *jets_AK5PF_neutralHadE;
   vector<float>   *jets_AK5PF_neutral_Mult;
   vector<float>   *jets_AK5PF_mu_Mult;
   vector<float>   *jets_AK5PF_ehf;
   vector<float>   *jets_AK5PF_corrFactorRaw;
   vector<float>   *jets_AK5PFclean_phi;
   vector<float>   *jets_AK5PFclean_pt;
   vector<float>   *jets_AK5PFclean_pz;
   vector<float>   *jets_AK5PFclean_px;
   vector<float>   *jets_AK5PFclean_py;
   vector<float>   *jets_AK5PFclean_eta;
   vector<float>   *jets_AK5PFclean_et;
   vector<float>   *jets_AK5PFclean_energy;
   vector<float>   *jets_AK5PFclean_parton_Id;
   vector<float>   *jets_AK5PFclean_parton_motherId;
   vector<float>   *jets_AK5PFclean_gen_pt;
   vector<float>   *jets_AK5PFclean_gen_phi;
   vector<float>   *jets_AK5PFclean_partonFlavour;
   vector<float>   *jets_AK5PFclean_btag_TC_highPur;
   vector<float>   *jets_AK5PFclean_btag_TC_highEff;
   vector<float>   *jets_AK5PFclean_btag_jetProb;
   vector<float>   *jets_AK5PFclean_btag_jetBProb;
   vector<float>   *jets_AK5PFclean_btag_secVertexHighPur;
   vector<float>   *jets_AK5PFclean_btag_secVertexHighEff;
   vector<float>   *jets_AK5PFclean_btag_secVertexCombined;
   vector<float>   *jets_AK5PFclean_jetCharge;
   vector<float>   *jets_AK5PFclean_chgEmE;
   vector<float>   *jets_AK5PFclean_chgHadE;
   vector<float>   *jets_AK5PFclean_photonEnergy;
   vector<float>   *jets_AK5PFclean_chg_Mult;
   vector<float>   *jets_AK5PFclean_neutralEmE;
   vector<float>   *jets_AK5PFclean_neutralHadE;
   vector<float>   *jets_AK5PFclean_neutral_Mult;
   vector<float>   *jets_AK5PFclean_mu_Mult;
   vector<float>   *jets_AK5PFclean_ehf;
   vector<float>   *jets_AK5PFclean_corrFactorRaw;
   vector<float>   *jets_AK5PFclean_rawPt;
   vector<float>   *mus_energy;
   vector<float>   *mus_et;
   vector<float>   *mus_eta;
   vector<float>   *mus_phi;
   vector<float>   *mus_pt;
   vector<float>   *mus_px;
   vector<float>   *mus_py;
   vector<float>   *mus_pz;
   vector<float>   *mus_cIso;
   vector<float>   *mus_tIso;
   vector<float>   *mus_ecalIso;
   vector<float>   *mus_hcalIso;
   vector<float>   *mus_ecalvetoDep;
   vector<float>   *mus_hcalvetoDep;
   vector<float>   *mus_pfIsolationR03_sumChargedHadronPt;
   vector<float>   *mus_pfIsolationR03_sumChargedParticlePt;
   vector<float>   *mus_pfIsolationR03_sumNeutralHadronEt;
   vector<float>   *mus_pfIsolationR03_sumNeutralHadronEtHighThreshold;
   vector<float>   *mus_pfIsolationR03_sumPhotonEt;
   vector<float>   *mus_pfIsolationR03_sumPhotonEtHighThreshold;
   vector<float>   *mus_pfIsolationR03_sumPUPt;
   vector<float>   *mus_pfIsolationR04_sumChargedHadronPt;
   vector<float>   *mus_pfIsolationR04_sumChargedParticlePt;
   vector<float>   *mus_pfIsolationR04_sumNeutralHadronEt;
   vector<float>   *mus_pfIsolationR04_sumNeutralHadronEtHighThreshold;
   vector<float>   *mus_pfIsolationR04_sumPhotonEt;
   vector<float>   *mus_pfIsolationR04_sumPhotonEtHighThreshold;
   vector<float>   *mus_pfIsolationR04_sumPUPt;
   vector<float>   *mus_charge;
   vector<float>   *mus_cm_chi2;
   vector<float>   *mus_cm_ndof;
   vector<float>   *mus_cm_pt;
   vector<float>   *mus_cm_ptErr;
   vector<float>   *mus_tk_chi2;
   vector<float>   *mus_tk_ndof;
   vector<float>   *mus_tk_pt;
   vector<float>   *mus_tk_px;
   vector<float>   *mus_tk_py;
   vector<float>   *mus_tk_pz;
   vector<float>   *mus_tk_phi;
   vector<float>   *mus_tk_d0dum;
   vector<float>   *mus_tk_vx;
   vector<float>   *mus_tk_vy;
   vector<float>   *mus_tk_vz;
   vector<float>   *mus_tk_numvalhits;
   vector<float>   *mus_tk_ptErr;
   vector<float>   *mus_tk_numvalPixelhits;
   vector<float>   *mus_tk_numpixelWthMeasr;
   vector<float>   *mus_stamu_pt;
   vector<float>   *mus_stamu_ptErr;
   vector<float>   *mus_num_matches;
   vector<float>   *mus_isPFMuon;
   vector<float>   *mus_isTrackerMuon;
   vector<float>   *mus_isGlobalMuon;
   vector<float>   *mus_id_AllGlobalMuons;
   vector<float>   *mus_id_AllTrackerMuons;
   vector<float>   *mus_id_GlobalMuonPromptTight;
   vector<float>   *mus_tk_LayersWithMeasurement;
   vector<float>   *mus_dB;
   vector<float>   *mus_numberOfMatchedStations;
   vector<float>   *pfTypeImets_et;
   vector<float>   *pfTypeImets_phi;
   vector<float>   *pfTypeImets_ex;
   vector<float>   *pfTypeImets_ey;
   vector<float>   *pfTypeImets_gen_et;
   vector<float>   *pfTypeImets_gen_phi;
   vector<float>   *pfTypeImets_sumEt;
   vector<float>   *pf_els_energy;
   vector<float>   *pf_els_et;
   vector<float>   *pf_els_eta;
   vector<float>   *pf_els_phi;
   vector<float>   *pf_els_pt;
   vector<float>   *pf_els_px;
   vector<float>   *pf_els_py;
   vector<float>   *pf_els_pz;
   vector<float>   *pf_els_robustTightId;
   vector<float>   *pf_els_simpleEleId95relIso;
   vector<float>   *pf_els_simpleEleId90relIso;
   vector<float>   *pf_els_simpleEleId85relIso;
   vector<float>   *pf_els_simpleEleId80relIso;
   vector<float>   *pf_els_simpleEleId70relIso;
   vector<float>   *pf_els_simpleEleId95cIso;
   vector<float>   *pf_els_simpleEleId90cIso;
   vector<float>   *pf_els_simpleEleId85cIso;
   vector<float>   *pf_els_simpleEleId80cIso;
   vector<float>   *pf_els_simpleEleId70cIso;
   vector<float>   *pf_els_cIso;
   vector<float>   *pf_els_tIso;
   vector<float>   *pf_els_ecalIso;
   vector<float>   *pf_els_hcalIso;
   vector<float>   *pf_els_chargedHadronIso;
   vector<float>   *pf_els_photonIso;
   vector<float>   *pf_els_neutralHadronIso;
   vector<float>   *pf_els_charge;
   vector<float>   *pf_els_hadOverEm;
   vector<float>   *pf_els_eOverPIn;
   vector<float>   *pf_els_sigmaIEtaIEta;
   vector<float>   *pf_els_scEnergy;
   vector<float>   *pf_els_scEta;
   vector<float>   *pf_els_scE1x5;
   vector<float>   *pf_els_scE2x5Max;
   vector<float>   *pf_els_scE5x5;
   vector<float>   *pf_els_isEB;
   vector<float>   *pf_els_isEE;
   vector<float>   *pf_els_dEtaIn;
   vector<float>   *pf_els_dPhiIn;
   vector<float>   *pf_els_dEtaOut;
   vector<float>   *pf_els_dPhiOut;
   vector<float>   *pf_els_numlosthits;
   vector<float>   *pf_els_tk_phi;
   vector<float>   *pf_els_d0dum;
   vector<float>   *pf_els_vx;
   vector<float>   *pf_els_vy;
   vector<float>   *pf_els_vz;
   vector<float>   *pf_els_ptError;
   vector<float>   *pf_els_n_inner_layer;
   vector<float>   *pf_els_dr03EcalRecHitSumEt;
   vector<float>   *pf_els_dr03HcalTowerSumEt;
   vector<float>   *pf_els_dr03HcalDepth1TowerSumEt;
   vector<float>   *pf_els_dr03HcalDepth2TowerSumEt;
   vector<float>   *pf_els_dr03TkSumPt;
   vector<float>   *pf_mus_energy;
   vector<float>   *pf_mus_et;
   vector<float>   *pf_mus_eta;
   vector<float>   *pf_mus_phi;
   vector<float>   *pf_mus_pt;
   vector<float>   *pf_mus_px;
   vector<float>   *pf_mus_py;
   vector<float>   *pf_mus_pz;
   vector<float>   *pf_mus_cIso;
   vector<float>   *pf_mus_tIso;
   vector<float>   *pf_mus_ecalIso;
   vector<float>   *pf_mus_hcalIso;
   vector<float>   *pf_mus_neutralHadronIso;
   vector<float>   *pf_mus_chargedHadronIso;
   vector<float>   *pf_mus_photonIso;
   vector<float>   *pf_mus_charge;
   vector<float>   *pf_mus_cm_chi2;
   vector<float>   *pf_mus_cm_ndof;
   vector<float>   *pf_mus_cm_pt;
   vector<float>   *pf_mus_cm_ptErr;
   vector<float>   *pf_mus_tk_chi2;
   vector<float>   *pf_mus_tk_ndof;
   vector<float>   *pf_mus_tk_pt;
   vector<float>   *pf_mus_tk_phi;
   vector<float>   *pf_mus_tk_d0dum;
   vector<float>   *pf_mus_tk_vz;
   vector<float>   *pf_mus_tk_numvalhits;
   vector<float>   *pf_mus_tk_ptErr;
   vector<float>   *pf_mus_tk_numvalPixelhits;
   vector<float>   *pf_mus_tk_numpixelWthMeasr;
   vector<float>   *pf_mus_stamu_pt;
   vector<float>   *pf_mus_stamu_ptErr;
   vector<float>   *pf_mus_num_matches;
   vector<float>   *pf_mus_isTrackerMuon;
   vector<float>   *pf_mus_isGlobalMuon;
   vector<float>   *pf_mus_id_GlobalMuonPromptTight;
   vector<float>   *pfcand_pdgId;
   vector<float>   *pfcand_particleId;
   vector<float>   *pfcand_pt;
   vector<float>   *pfcand_pz;
   vector<float>   *pfcand_px;
   vector<float>   *pfcand_py;
   vector<float>   *pfcand_eta;
   vector<float>   *pfcand_phi;
   vector<float>   *pfcand_theta;
   vector<float>   *pfcand_energy;
   vector<float>   *pfcand_charge;
   vector<float>   *pfmets_et;
   vector<float>   *pfmets_phi;
   vector<float>   *pfmets_ex;
   vector<float>   *pfmets_ey;
   vector<float>   *pfmets_gen_et;
   vector<float>   *pfmets_gen_phi;
   vector<float>   *pfmets_sumEt;
   UInt_t          Npv;
   vector<float>   *pv_x;
   vector<float>   *pv_y;
   vector<float>   *pv_z;
   vector<float>   *pv_ndof;
   vector<float>   *pv_isFake;
   vector<float>   *pv_tracksSize;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          bunchCrossing;
   Float_t         weight;
   string          *model_params;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT10_px;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT10_py;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT10_pz;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT10_energy;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT10_phi;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT10_eta;
   vector<vector<int> > *fastjets_AK5PFclean_R1p2_R0p5pT10_index;
   vector<int>     *fastjets_AK5PFclean_R1p2_R0p5pT10_nconstituents;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT15_px;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT15_py;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT15_pz;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT15_energy;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT15_phi;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT15_eta;
   vector<vector<int> > *fastjets_AK5PFclean_R1p2_R0p5pT15_index;
   vector<int>     *fastjets_AK5PFclean_R1p2_R0p5pT15_nconstituents;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT20_px;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT20_py;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT20_pz;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT20_energy;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT20_phi;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT20_eta;
   vector<vector<int> > *fastjets_AK5PFclean_R1p2_R0p5pT20_index;
   vector<int>     *fastjets_AK5PFclean_R1p2_R0p5pT20_nconstituents;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT25_px;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT25_py;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT25_pz;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT25_energy;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT25_phi;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT25_eta;
   vector<vector<int> > *fastjets_AK5PFclean_R1p2_R0p5pT25_index;
   vector<int>     *fastjets_AK5PFclean_R1p2_R0p5pT25_nconstituents;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT30_px;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT30_py;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT30_pz;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT30_energy;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT30_phi;
   vector<float>   *fastjets_AK5PFclean_R1p2_R0p5pT30_eta;
   vector<vector<int> > *fastjets_AK5PFclean_R1p2_R0p5pT30_index;
   vector<int>     *fastjets_AK5PFclean_R1p2_R0p5pT30_nconstituents;

   // List of branches
   TBranch        *b_els_conversion_dist;   //!
   TBranch        *b_els_conversion_dcot;   //!
   TBranch        *b_els_PFchargedHadronIsoR03;   //!
   TBranch        *b_els_PFphotonIsoR03;   //!
   TBranch        *b_els_PFneutralHadronIsoR03;   //!
   TBranch        *b_els_hasMatchedConversion;   //!
   TBranch        *b_PU_TrueNumInteractions;   //!
   TBranch        *b_rho_kt6PFJetsForIsolation2011;   //!

   TBranch        *b_beamSpot_x;   //!
   TBranch        *b_beamSpot_y;   //!
   TBranch        *b_els_energy;   //!
   TBranch        *b_els_et;   //!
   TBranch        *b_els_eta;   //!
   TBranch        *b_els_phi;   //!
   TBranch        *b_els_pt;   //!
   TBranch        *b_els_px;   //!
   TBranch        *b_els_py;   //!
   TBranch        *b_els_pz;   //!
   TBranch        *b_els_robustTightId;   //!
   TBranch        *b_els_simpleEleId95relIso;   //!
   TBranch        *b_els_simpleEleId90relIso;   //!
   TBranch        *b_els_simpleEleId85relIso;   //!
   TBranch        *b_els_simpleEleId80relIso;   //!
   TBranch        *b_els_simpleEleId70relIso;   //!
   TBranch        *b_els_simpleEleId95cIso;   //!
   TBranch        *b_els_simpleEleId90cIso;   //!
   TBranch        *b_els_simpleEleId85cIso;   //!
   TBranch        *b_els_simpleEleId80cIso;   //!
   TBranch        *b_els_simpleEleId70cIso;   //!
   TBranch        *b_els_cIso;   //!
   TBranch        *b_els_tIso;   //!
   TBranch        *b_els_ecalIso;   //!
   TBranch        *b_els_hcalIso;   //!
   TBranch        *b_els_charge;   //!
   TBranch        *b_els_caloEnergy;   //!
   TBranch        *b_els_hadOverEm;   //!
   TBranch        *b_els_eOverPIn;   //!
   TBranch        *b_els_sigmaIEtaIEta;   //!
   TBranch        *b_els_scEnergy;   //!
   TBranch        *b_els_scEta;   //!
   TBranch        *b_els_scE1x5;   //!
   TBranch        *b_els_scE2x5Max;   //!
   TBranch        *b_els_scE5x5;   //!
   TBranch        *b_els_isEB;   //!
   TBranch        *b_els_isEE;   //!
   TBranch        *b_els_dEtaIn;   //!
   TBranch        *b_els_dPhiIn;   //!
   TBranch        *b_els_dEtaOut;   //!
   TBranch        *b_els_dPhiOut;   //!
   TBranch        *b_els_numlosthits;   //!
   TBranch        *b_els_tk_pz;   //!
   TBranch        *b_els_tk_pt;   //!
   TBranch        *b_els_tk_phi;   //!
   TBranch        *b_els_d0dum;   //!
   TBranch        *b_els_vx;   //!
   TBranch        *b_els_vy;   //!
   TBranch        *b_els_vz;   //!
   TBranch        *b_els_ptError;   //!
   TBranch        *b_els_n_inner_layer;   //!
   TBranch        *b_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_els_dr03TkSumPt;   //!
   TBranch        *b_jets_AK5PF_phi;   //!
   TBranch        *b_jets_AK5PF_pt;   //!
   TBranch        *b_jets_AK5PF_pz;   //!
   TBranch        *b_jets_AK5PF_px;   //!
   TBranch        *b_jets_AK5PF_py;   //!
   TBranch        *b_jets_AK5PF_eta;   //!
   TBranch        *b_jets_AK5PF_et;   //!
   TBranch        *b_jets_AK5PF_energy;   //!
   TBranch        *b_jets_AK5PF_parton_Id;   //!
   TBranch        *b_jets_AK5PF_parton_motherId;   //!
   TBranch        *b_jets_AK5PF_gen_pt;   //!
   TBranch        *b_jets_AK5PF_gen_phi;   //!
   TBranch        *b_jets_AK5PF_partonFlavour;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PF_btag_jetProb;   //!
   TBranch        *b_jets_AK5PF_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PF_jetCharge;   //!
   TBranch        *b_jets_AK5PF_chgEmE;   //!
   TBranch        *b_jets_AK5PF_chgHadE;   //!
   TBranch        *b_jets_AK5PF_photonEnergy;   //!
   TBranch        *b_jets_AK5PF_chg_Mult;   //!
   TBranch        *b_jets_AK5PF_neutralEmE;   //!
   TBranch        *b_jets_AK5PF_neutralHadE;   //!
   TBranch        *b_jets_AK5PF_neutral_Mult;   //!
   TBranch        *b_jets_AK5PF_mu_Mult;   //!
   TBranch        *b_jets_AK5PF_ehf;   //!
   TBranch        *b_jets_AK5PF_corrFactorRaw;   //!
   TBranch        *b_jets_AK5PFclean_phi;   //!
   TBranch        *b_jets_AK5PFclean_pt;   //!
   TBranch        *b_jets_AK5PFclean_pz;   //!
   TBranch        *b_jets_AK5PFclean_px;   //!
   TBranch        *b_jets_AK5PFclean_py;   //!
   TBranch        *b_jets_AK5PFclean_eta;   //!
   TBranch        *b_jets_AK5PFclean_et;   //!
   TBranch        *b_jets_AK5PFclean_energy;   //!
   TBranch        *b_jets_AK5PFclean_parton_Id;   //!
   TBranch        *b_jets_AK5PFclean_parton_motherId;   //!
   TBranch        *b_jets_AK5PFclean_gen_pt;   //!
   TBranch        *b_jets_AK5PFclean_gen_phi;   //!
   TBranch        *b_jets_AK5PFclean_partonFlavour;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PFclean_jetCharge;   //!
   TBranch        *b_jets_AK5PFclean_chgEmE;   //!
   TBranch        *b_jets_AK5PFclean_chgHadE;   //!
   TBranch        *b_jets_AK5PFclean_photonEnergy;   //!
   TBranch        *b_jets_AK5PFclean_chg_Mult;   //!
   TBranch        *b_jets_AK5PFclean_neutralEmE;   //!
   TBranch        *b_jets_AK5PFclean_neutralHadE;   //!
   TBranch        *b_jets_AK5PFclean_neutral_Mult;   //!
   TBranch        *b_jets_AK5PFclean_mu_Mult;   //!
   TBranch        *b_jets_AK5PFclean_ehf;   //!
   TBranch        *b_jets_AK5PFclean_corrFactorRaw;   //!
   TBranch        *b_jets_AK5PFclean_rawPt;   //!
   TBranch        *b_mus_energy;   //!
   TBranch        *b_mus_et;   //!
   TBranch        *b_mus_eta;   //!
   TBranch        *b_mus_phi;   //!
   TBranch        *b_mus_pt;   //!
   TBranch        *b_mus_px;   //!
   TBranch        *b_mus_py;   //!
   TBranch        *b_mus_pz;   //!
   TBranch        *b_mus_cIso;   //!
   TBranch        *b_mus_tIso;   //!
   TBranch        *b_mus_ecalIso;   //!
   TBranch        *b_mus_hcalIso;   //!
   TBranch        *b_mus_ecalvetoDep;   //!
   TBranch        *b_mus_hcalvetoDep;   //!
   TBranch        *b_mus_pfIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_mus_pfIsolationR03_sumChargedParticlePt;   //!
   TBranch        *b_mus_pfIsolationR03_sumNeutralHadronEt;   //!
   TBranch        *b_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR03_sumPhotonEt;   //!
   TBranch        *b_mus_pfIsolationR03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR03_sumPUPt;   //!
   TBranch        *b_mus_pfIsolationR04_sumChargedHadronPt;   //!
   TBranch        *b_mus_pfIsolationR04_sumChargedParticlePt;   //!
   TBranch        *b_mus_pfIsolationR04_sumNeutralHadronEt;   //!
   TBranch        *b_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR04_sumPhotonEt;   //!
   TBranch        *b_mus_pfIsolationR04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mus_pfIsolationR04_sumPUPt;   //!
   TBranch        *b_mus_charge;   //!
   TBranch        *b_mus_cm_chi2;   //!
   TBranch        *b_mus_cm_ndof;   //!
   TBranch        *b_mus_cm_pt;   //!
   TBranch        *b_mus_cm_ptErr;   //!
   TBranch        *b_mus_tk_chi2;   //!
   TBranch        *b_mus_tk_ndof;   //!
   TBranch        *b_mus_tk_pt;   //!
   TBranch        *b_mus_tk_px;   //!
   TBranch        *b_mus_tk_py;   //!
   TBranch        *b_mus_tk_pz;   //!
   TBranch        *b_mus_tk_phi;   //!
   TBranch        *b_mus_tk_d0dum;   //!
   TBranch        *b_mus_tk_vx;   //!
   TBranch        *b_mus_tk_vy;   //!
   TBranch        *b_mus_tk_vz;   //!
   TBranch        *b_mus_tk_numvalhits;   //!
   TBranch        *b_mus_tk_ptErr;   //!
   TBranch        *b_mus_tk_numvalPixelhits;   //!
   TBranch        *b_mus_tk_numpixelWthMeasr;   //!
   TBranch        *b_mus_stamu_pt;   //!
   TBranch        *b_mus_stamu_ptErr;   //!
   TBranch        *b_mus_num_matches;   //!
   TBranch        *b_mus_isPFMuon;   //!
   TBranch        *b_mus_isTrackerMuon;   //!
   TBranch        *b_mus_isGlobalMuon;   //!
   TBranch        *b_mus_id_AllGlobalMuons;   //!
   TBranch        *b_mus_id_AllTrackerMuons;   //!
   TBranch        *b_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_mus_tk_LayersWithMeasurement;   //!
   TBranch        *b_mus_dB;   //!
   TBranch        *b_mus_numberOfMatchedStations;   //!
   TBranch        *b_pfTypeImets_et;   //!
   TBranch        *b_pfTypeImets_phi;   //!
   TBranch        *b_pfTypeImets_ex;   //!
   TBranch        *b_pfTypeImets_ey;   //!
   TBranch        *b_pfTypeImets_gen_et;   //!
   TBranch        *b_pfTypeImets_gen_phi;   //!
   TBranch        *b_pfTypeImets_sumEt;   //!
   TBranch        *b_pf_els_energy;   //!
   TBranch        *b_pf_els_et;   //!
   TBranch        *b_pf_els_eta;   //!
   TBranch        *b_pf_els_phi;   //!
   TBranch        *b_pf_els_pt;   //!
   TBranch        *b_pf_els_px;   //!
   TBranch        *b_pf_els_py;   //!
   TBranch        *b_pf_els_pz;   //!
   TBranch        *b_pf_els_robustTightId;   //!
   TBranch        *b_pf_els_simpleEleId95relIso;   //!
   TBranch        *b_pf_els_simpleEleId90relIso;   //!
   TBranch        *b_pf_els_simpleEleId85relIso;   //!
   TBranch        *b_pf_els_simpleEleId80relIso;   //!
   TBranch        *b_pf_els_simpleEleId70relIso;   //!
   TBranch        *b_pf_els_simpleEleId95cIso;   //!
   TBranch        *b_pf_els_simpleEleId90cIso;   //!
   TBranch        *b_pf_els_simpleEleId85cIso;   //!
   TBranch        *b_pf_els_simpleEleId80cIso;   //!
   TBranch        *b_pf_els_simpleEleId70cIso;   //!
   TBranch        *b_pf_els_cIso;   //!
   TBranch        *b_pf_els_tIso;   //!
   TBranch        *b_pf_els_ecalIso;   //!
   TBranch        *b_pf_els_hcalIso;   //!
   TBranch        *b_pf_els_chargedHadronIso;   //!
   TBranch        *b_pf_els_photonIso;   //!
   TBranch        *b_pf_els_neutralHadronIso;   //!
   TBranch        *b_pf_els_charge;   //!
   TBranch        *b_pf_els_hadOverEm;   //!
   TBranch        *b_pf_els_eOverPIn;   //!
   TBranch        *b_pf_els_sigmaIEtaIEta;   //!
   TBranch        *b_pf_els_scEnergy;   //!
   TBranch        *b_pf_els_scEta;   //!
   TBranch        *b_pf_els_scE1x5;   //!
   TBranch        *b_pf_els_scE2x5Max;   //!
   TBranch        *b_pf_els_scE5x5;   //!
   TBranch        *b_pf_els_isEB;   //!
   TBranch        *b_pf_els_isEE;   //!
   TBranch        *b_pf_els_dEtaIn;   //!
   TBranch        *b_pf_els_dPhiIn;   //!
   TBranch        *b_pf_els_dEtaOut;   //!
   TBranch        *b_pf_els_dPhiOut;   //!
   TBranch        *b_pf_els_numlosthits;   //!
   TBranch        *b_pf_els_tk_phi;   //!
   TBranch        *b_pf_els_d0dum;   //!
   TBranch        *b_pf_els_vx;   //!
   TBranch        *b_pf_els_vy;   //!
   TBranch        *b_pf_els_vz;   //!
   TBranch        *b_pf_els_ptError;   //!
   TBranch        *b_pf_els_n_inner_layer;   //!
   TBranch        *b_pf_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_pf_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_pf_els_dr03TkSumPt;   //!
   TBranch        *b_pf_mus_energy;   //!
   TBranch        *b_pf_mus_et;   //!
   TBranch        *b_pf_mus_eta;   //!
   TBranch        *b_pf_mus_phi;   //!
   TBranch        *b_pf_mus_pt;   //!
   TBranch        *b_pf_mus_px;   //!
   TBranch        *b_pf_mus_py;   //!
   TBranch        *b_pf_mus_pz;   //!
   TBranch        *b_pf_mus_cIso;   //!
   TBranch        *b_pf_mus_tIso;   //!
   TBranch        *b_pf_mus_ecalIso;   //!
   TBranch        *b_pf_mus_hcalIso;   //!
   TBranch        *b_pf_mus_neutralHadronIso;   //!
   TBranch        *b_pf_mus_chargedHadronIso;   //!
   TBranch        *b_pf_mus_photonIso;   //!
   TBranch        *b_pf_mus_charge;   //!
   TBranch        *b_pf_mus_cm_chi2;   //!
   TBranch        *b_pf_mus_cm_ndof;   //!
   TBranch        *b_pf_mus_cm_pt;   //!
   TBranch        *b_pf_mus_cm_ptErr;   //!
   TBranch        *b_pf_mus_tk_chi2;   //!
   TBranch        *b_pf_mus_tk_ndof;   //!
   TBranch        *b_pf_mus_tk_pt;   //!
   TBranch        *b_pf_mus_tk_phi;   //!
   TBranch        *b_pf_mus_tk_d0dum;   //!
   TBranch        *b_pf_mus_tk_vz;   //!
   TBranch        *b_pf_mus_tk_numvalhits;   //!
   TBranch        *b_pf_mus_tk_ptErr;   //!
   TBranch        *b_pf_mus_tk_numvalPixelhits;   //!
   TBranch        *b_pf_mus_tk_numpixelWthMeasr;   //!
   TBranch        *b_pf_mus_stamu_pt;   //!
   TBranch        *b_pf_mus_stamu_ptErr;   //!
   TBranch        *b_pf_mus_num_matches;   //!
   TBranch        *b_pf_mus_isTrackerMuon;   //!
   TBranch        *b_pf_mus_isGlobalMuon;   //!
   TBranch        *b_pf_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_pfcand_pdgId;   //!
   TBranch        *b_pfcand_particleId;   //!
   TBranch        *b_pfcand_pt;   //!
   TBranch        *b_pfcand_pz;   //!
   TBranch        *b_pfcand_px;   //!
   TBranch        *b_pfcand_py;   //!
   TBranch        *b_pfcand_eta;   //!
   TBranch        *b_pfcand_phi;   //!
   TBranch        *b_pfcand_theta;   //!
   TBranch        *b_pfcand_energy;   //!
   TBranch        *b_pfcand_charge;   //!
   TBranch        *b_pfmets_et;   //!
   TBranch        *b_pfmets_phi;   //!
   TBranch        *b_pfmets_ex;   //!
   TBranch        *b_pfmets_ey;   //!
   TBranch        *b_pfmets_gen_et;   //!
   TBranch        *b_pfmets_gen_phi;   //!
   TBranch        *b_pfmets_sumEt;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_pv_tracksSize;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_model_params;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_px;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_py;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_pz;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_energy;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_phi;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_eta;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_index;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT10_nconstituents;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_px;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_py;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_pz;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_energy;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_phi;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_eta;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_index;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT15_nconstituents;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_px;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_py;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_pz;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_energy;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_phi;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_eta;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_index;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT20_nconstituents;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_px;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_py;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_pz;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_energy;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_phi;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_eta;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_index;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT25_nconstituents;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_px;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_py;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_pz;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_energy;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_phi;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_eta;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_index;   //!
   TBranch        *b_fastjets_AK5PFclean_R1p2_R0p5pT30_nconstituents;   //!

void InitializeB(TChain *fChain)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   beamSpot_x = 0;
   beamSpot_y = 0;
   els_energy = 0;
   els_et = 0;
   els_eta = 0;
   els_phi = 0;
   els_pt = 0;
   els_px = 0;
   els_py = 0;
   els_pz = 0;
   els_robustTightId = 0;
   els_simpleEleId95relIso = 0;
   els_simpleEleId90relIso = 0;
   els_simpleEleId85relIso = 0;
   els_simpleEleId80relIso = 0;
   els_simpleEleId70relIso = 0;
   els_simpleEleId95cIso = 0;
   els_simpleEleId90cIso = 0;
   els_simpleEleId85cIso = 0;
   els_simpleEleId80cIso = 0;
   els_simpleEleId70cIso = 0;
   els_cIso = 0;
   els_tIso = 0;
   els_ecalIso = 0;
   els_hcalIso = 0;
   els_charge = 0;
   els_caloEnergy = 0;
   els_hadOverEm = 0;
   els_eOverPIn = 0;
   els_sigmaIEtaIEta = 0;
   els_scEnergy = 0;
   els_scEta = 0;
   els_scE1x5 = 0;
   els_scE2x5Max = 0;
   els_scE5x5 = 0;
   els_isEB = 0;
   els_isEE = 0;
   els_dEtaIn = 0;
   els_dPhiIn = 0;
   els_dEtaOut = 0;
   els_dPhiOut = 0;
   els_numlosthits = 0;
   els_tk_pz = 0;
   els_tk_pt = 0;
   els_tk_phi = 0;
   els_d0dum = 0;
   els_vx = 0;
   els_vy = 0;
   els_vz = 0;
   els_ptError = 0;
   els_n_inner_layer = 0;
   els_dr03EcalRecHitSumEt = 0;
   els_dr03HcalTowerSumEt = 0;
   els_dr03HcalDepth1TowerSumEt = 0;
   els_dr03HcalDepth2TowerSumEt = 0;
   els_dr03TkSumPt = 0;
   jets_AK5PF_phi = 0;
   jets_AK5PF_pt = 0;
   jets_AK5PF_pz = 0;
   jets_AK5PF_px = 0;
   jets_AK5PF_py = 0;
   jets_AK5PF_eta = 0;
   jets_AK5PF_et = 0;
   jets_AK5PF_energy = 0;
   jets_AK5PF_parton_Id = 0;
   jets_AK5PF_parton_motherId = 0;
   jets_AK5PF_gen_pt = 0;
   jets_AK5PF_gen_phi = 0;
   jets_AK5PF_partonFlavour = 0;
   jets_AK5PF_btag_TC_highPur = 0;
   jets_AK5PF_btag_TC_highEff = 0;
   jets_AK5PF_btag_jetProb = 0;
   jets_AK5PF_btag_jetBProb = 0;
   jets_AK5PF_btag_secVertexHighPur = 0;
   jets_AK5PF_btag_secVertexHighEff = 0;
   jets_AK5PF_btag_secVertexCombined = 0;
   jets_AK5PF_jetCharge = 0;
   jets_AK5PF_chgEmE = 0;
   jets_AK5PF_chgHadE = 0;
   jets_AK5PF_photonEnergy = 0;
   jets_AK5PF_chg_Mult = 0;
   jets_AK5PF_neutralEmE = 0;
   jets_AK5PF_neutralHadE = 0;
   jets_AK5PF_neutral_Mult = 0;
   jets_AK5PF_mu_Mult = 0;
   jets_AK5PF_ehf = 0;
   jets_AK5PF_corrFactorRaw = 0;
   jets_AK5PFclean_phi = 0;
   jets_AK5PFclean_pt = 0;
   jets_AK5PFclean_pz = 0;
   jets_AK5PFclean_px = 0;
   jets_AK5PFclean_py = 0;
   jets_AK5PFclean_eta = 0;
   jets_AK5PFclean_et = 0;
   jets_AK5PFclean_energy = 0;
   jets_AK5PFclean_parton_Id = 0;
   jets_AK5PFclean_parton_motherId = 0;
   jets_AK5PFclean_gen_pt = 0;
   jets_AK5PFclean_gen_phi = 0;
   jets_AK5PFclean_partonFlavour = 0;
   jets_AK5PFclean_btag_TC_highPur = 0;
   jets_AK5PFclean_btag_TC_highEff = 0;
   jets_AK5PFclean_btag_jetProb = 0;
   jets_AK5PFclean_btag_jetBProb = 0;
   jets_AK5PFclean_btag_secVertexHighPur = 0;
   jets_AK5PFclean_btag_secVertexHighEff = 0;
   jets_AK5PFclean_btag_secVertexCombined = 0;
   jets_AK5PFclean_jetCharge = 0;
   jets_AK5PFclean_chgEmE = 0;
   jets_AK5PFclean_chgHadE = 0;
   jets_AK5PFclean_photonEnergy = 0;
   jets_AK5PFclean_chg_Mult = 0;
   jets_AK5PFclean_neutralEmE = 0;
   jets_AK5PFclean_neutralHadE = 0;
   jets_AK5PFclean_neutral_Mult = 0;
   jets_AK5PFclean_mu_Mult = 0;
   jets_AK5PFclean_ehf = 0;
   jets_AK5PFclean_corrFactorRaw = 0;
   jets_AK5PFclean_rawPt = 0;
   mus_energy = 0;
   mus_et = 0;
   mus_eta = 0;
   mus_phi = 0;
   mus_pt = 0;
   mus_px = 0;
   mus_py = 0;
   mus_pz = 0;
   mus_cIso = 0;
   mus_tIso = 0;
   mus_ecalIso = 0;
   mus_hcalIso = 0;
   mus_ecalvetoDep = 0;
   mus_hcalvetoDep = 0;
   mus_pfIsolationR03_sumChargedHadronPt = 0;
   mus_pfIsolationR03_sumChargedParticlePt = 0;
   mus_pfIsolationR03_sumNeutralHadronEt = 0;
   mus_pfIsolationR03_sumNeutralHadronEtHighThreshold = 0;
   mus_pfIsolationR03_sumPhotonEt = 0;
   mus_pfIsolationR03_sumPhotonEtHighThreshold = 0;
   mus_pfIsolationR03_sumPUPt = 0;
   mus_pfIsolationR04_sumChargedHadronPt = 0;
   mus_pfIsolationR04_sumChargedParticlePt = 0;
   mus_pfIsolationR04_sumNeutralHadronEt = 0;
   mus_pfIsolationR04_sumNeutralHadronEtHighThreshold = 0;
   mus_pfIsolationR04_sumPhotonEt = 0;
   mus_pfIsolationR04_sumPhotonEtHighThreshold = 0;
   mus_pfIsolationR04_sumPUPt = 0;
   mus_charge = 0;
   mus_cm_chi2 = 0;
   mus_cm_ndof = 0;
   mus_cm_pt = 0;
   mus_cm_ptErr = 0;
   mus_tk_chi2 = 0;
   mus_tk_ndof = 0;
   mus_tk_pt = 0;
   mus_tk_px = 0;
   mus_tk_py = 0;
   mus_tk_pz = 0;
   mus_tk_phi = 0;
   mus_tk_d0dum = 0;
   mus_tk_vx = 0;
   mus_tk_vy = 0;
   mus_tk_vz = 0;
   mus_tk_numvalhits = 0;
   mus_tk_ptErr = 0;
   mus_tk_numvalPixelhits = 0;
   mus_tk_numpixelWthMeasr = 0;
   mus_stamu_pt = 0;
   mus_stamu_ptErr = 0;
   mus_num_matches = 0;
   mus_isPFMuon = 0;
   mus_isTrackerMuon = 0;
   mus_isGlobalMuon = 0;
   mus_id_AllGlobalMuons = 0;
   mus_id_AllTrackerMuons = 0;
   mus_id_GlobalMuonPromptTight = 0;
   mus_tk_LayersWithMeasurement = 0;
   mus_dB = 0;
   mus_numberOfMatchedStations = 0;
   pfTypeImets_et = 0;
   pfTypeImets_phi = 0;
   pfTypeImets_ex = 0;
   pfTypeImets_ey = 0;
   pfTypeImets_gen_et = 0;
   pfTypeImets_gen_phi = 0;
   pfTypeImets_sumEt = 0;
   pf_els_energy = 0;
   pf_els_et = 0;
   pf_els_eta = 0;
   pf_els_phi = 0;
   pf_els_pt = 0;
   pf_els_px = 0;
   pf_els_py = 0;
   pf_els_pz = 0;
   pf_els_robustTightId = 0;
   pf_els_simpleEleId95relIso = 0;
   pf_els_simpleEleId90relIso = 0;
   pf_els_simpleEleId85relIso = 0;
   pf_els_simpleEleId80relIso = 0;
   pf_els_simpleEleId70relIso = 0;
   pf_els_simpleEleId95cIso = 0;
   pf_els_simpleEleId90cIso = 0;
   pf_els_simpleEleId85cIso = 0;
   pf_els_simpleEleId80cIso = 0;
   pf_els_simpleEleId70cIso = 0;
   pf_els_cIso = 0;
   pf_els_tIso = 0;
   pf_els_ecalIso = 0;
   pf_els_hcalIso = 0;
   pf_els_chargedHadronIso = 0;
   pf_els_photonIso = 0;
   pf_els_neutralHadronIso = 0;
   pf_els_charge = 0;
   pf_els_hadOverEm = 0;
   pf_els_eOverPIn = 0;
   pf_els_sigmaIEtaIEta = 0;
   pf_els_scEnergy = 0;
   pf_els_scEta = 0;
   pf_els_scE1x5 = 0;
   pf_els_scE2x5Max = 0;
   pf_els_scE5x5 = 0;
   pf_els_isEB = 0;
   pf_els_isEE = 0;
   pf_els_dEtaIn = 0;
   pf_els_dPhiIn = 0;
   pf_els_dEtaOut = 0;
   pf_els_dPhiOut = 0;
   pf_els_numlosthits = 0;
   pf_els_tk_phi = 0;
   pf_els_d0dum = 0;
   pf_els_vx = 0;
   pf_els_vy = 0;
   pf_els_vz = 0;
   pf_els_ptError = 0;
   pf_els_n_inner_layer = 0;
   pf_els_dr03EcalRecHitSumEt = 0;
   pf_els_dr03HcalTowerSumEt = 0;
   pf_els_dr03HcalDepth1TowerSumEt = 0;
   pf_els_dr03HcalDepth2TowerSumEt = 0;
   pf_els_dr03TkSumPt = 0;
   pf_mus_energy = 0;
   pf_mus_et = 0;
   pf_mus_eta = 0;
   pf_mus_phi = 0;
   pf_mus_pt = 0;
   pf_mus_px = 0;
   pf_mus_py = 0;
   pf_mus_pz = 0;
   pf_mus_cIso = 0;
   pf_mus_tIso = 0;
   pf_mus_ecalIso = 0;
   pf_mus_hcalIso = 0;
   pf_mus_neutralHadronIso = 0;
   pf_mus_chargedHadronIso = 0;
   pf_mus_photonIso = 0;
   pf_mus_charge = 0;
   pf_mus_cm_chi2 = 0;
   pf_mus_cm_ndof = 0;
   pf_mus_cm_pt = 0;
   pf_mus_cm_ptErr = 0;
   pf_mus_tk_chi2 = 0;
   pf_mus_tk_ndof = 0;
   pf_mus_tk_pt = 0;
   pf_mus_tk_phi = 0;
   pf_mus_tk_d0dum = 0;
   pf_mus_tk_vz = 0;
   pf_mus_tk_numvalhits = 0;
   pf_mus_tk_ptErr = 0;
   pf_mus_tk_numvalPixelhits = 0;
   pf_mus_tk_numpixelWthMeasr = 0;
   pf_mus_stamu_pt = 0;
   pf_mus_stamu_ptErr = 0;
   pf_mus_num_matches = 0;
   pf_mus_isTrackerMuon = 0;
   pf_mus_isGlobalMuon = 0;
   pf_mus_id_GlobalMuonPromptTight = 0;
   pfcand_pdgId = 0;
   pfcand_particleId = 0;
   pfcand_pt = 0;
   pfcand_pz = 0;
   pfcand_px = 0;
   pfcand_py = 0;
   pfcand_eta = 0;
   pfcand_phi = 0;
   pfcand_theta = 0;
   pfcand_energy = 0;
   pfcand_charge = 0;
   pfmets_et = 0;
   pfmets_phi = 0;
   pfmets_ex = 0;
   pfmets_ey = 0;
   pfmets_gen_et = 0;
   pfmets_gen_phi = 0;
   pfmets_sumEt = 0;
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_ndof = 0;
   pv_isFake = 0;
   pv_tracksSize = 0;
   model_params = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_px = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_py = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_pz = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_energy = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_phi = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_eta = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_index = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT10_nconstituents = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_px = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_py = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_pz = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_energy = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_phi = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_eta = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_index = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT15_nconstituents = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_px = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_py = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_pz = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_energy = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_phi = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_eta = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_index = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT20_nconstituents = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_px = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_py = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_pz = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_energy = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_phi = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_eta = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_index = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT25_nconstituents = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_px = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_py = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_pz = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_energy = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_phi = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_eta = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_index = 0;
   fastjets_AK5PFclean_R1p2_R0p5pT30_nconstituents = 0;

   fChain->SetBranchAddress("beamSpot_x", &beamSpot_x, &b_beamSpot_x);
   fChain->SetBranchAddress("beamSpot_y", &beamSpot_y, &b_beamSpot_y);
   fChain->SetBranchAddress("els_energy", &els_energy, &b_els_energy);
   fChain->SetBranchAddress("els_et", &els_et, &b_els_et);
   fChain->SetBranchAddress("els_eta", &els_eta, &b_els_eta);
   fChain->SetBranchAddress("els_phi", &els_phi, &b_els_phi);
   fChain->SetBranchAddress("els_pt", &els_pt, &b_els_pt);
   fChain->SetBranchAddress("els_px", &els_px, &b_els_px);
   fChain->SetBranchAddress("els_py", &els_py, &b_els_py);
   fChain->SetBranchAddress("els_pz", &els_pz, &b_els_pz);
   fChain->SetBranchAddress("els_robustTightId", &els_robustTightId, &b_els_robustTightId);
   fChain->SetBranchAddress("els_simpleEleId95relIso", &els_simpleEleId95relIso, &b_els_simpleEleId95relIso);
   fChain->SetBranchAddress("els_simpleEleId90relIso", &els_simpleEleId90relIso, &b_els_simpleEleId90relIso);
   fChain->SetBranchAddress("els_simpleEleId85relIso", &els_simpleEleId85relIso, &b_els_simpleEleId85relIso);
   fChain->SetBranchAddress("els_simpleEleId80relIso", &els_simpleEleId80relIso, &b_els_simpleEleId80relIso);
   fChain->SetBranchAddress("els_simpleEleId70relIso", &els_simpleEleId70relIso, &b_els_simpleEleId70relIso);
   fChain->SetBranchAddress("els_simpleEleId95cIso", &els_simpleEleId95cIso, &b_els_simpleEleId95cIso);
   fChain->SetBranchAddress("els_simpleEleId90cIso", &els_simpleEleId90cIso, &b_els_simpleEleId90cIso);
   fChain->SetBranchAddress("els_simpleEleId85cIso", &els_simpleEleId85cIso, &b_els_simpleEleId85cIso);
   fChain->SetBranchAddress("els_simpleEleId80cIso", &els_simpleEleId80cIso, &b_els_simpleEleId80cIso);
   fChain->SetBranchAddress("els_simpleEleId70cIso", &els_simpleEleId70cIso, &b_els_simpleEleId70cIso);
   fChain->SetBranchAddress("els_cIso", &els_cIso, &b_els_cIso);
   fChain->SetBranchAddress("els_tIso", &els_tIso, &b_els_tIso);
   fChain->SetBranchAddress("els_ecalIso", &els_ecalIso, &b_els_ecalIso);
   fChain->SetBranchAddress("els_hcalIso", &els_hcalIso, &b_els_hcalIso);
   fChain->SetBranchAddress("els_charge", &els_charge, &b_els_charge);
   fChain->SetBranchAddress("els_caloEnergy", &els_caloEnergy, &b_els_caloEnergy);
   fChain->SetBranchAddress("els_hadOverEm", &els_hadOverEm, &b_els_hadOverEm);
   fChain->SetBranchAddress("els_eOverPIn", &els_eOverPIn, &b_els_eOverPIn);
   fChain->SetBranchAddress("els_sigmaIEtaIEta", &els_sigmaIEtaIEta, &b_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("els_scEnergy", &els_scEnergy, &b_els_scEnergy);
   fChain->SetBranchAddress("els_scEta", &els_scEta, &b_els_scEta);
   fChain->SetBranchAddress("els_scE1x5", &els_scE1x5, &b_els_scE1x5);
   fChain->SetBranchAddress("els_scE2x5Max", &els_scE2x5Max, &b_els_scE2x5Max);
   fChain->SetBranchAddress("els_scE5x5", &els_scE5x5, &b_els_scE5x5);
   fChain->SetBranchAddress("els_isEB", &els_isEB, &b_els_isEB);
   fChain->SetBranchAddress("els_isEE", &els_isEE, &b_els_isEE);
   fChain->SetBranchAddress("els_dEtaIn", &els_dEtaIn, &b_els_dEtaIn);
   fChain->SetBranchAddress("els_dPhiIn", &els_dPhiIn, &b_els_dPhiIn);
   fChain->SetBranchAddress("els_dEtaOut", &els_dEtaOut, &b_els_dEtaOut);
   fChain->SetBranchAddress("els_dPhiOut", &els_dPhiOut, &b_els_dPhiOut);
   fChain->SetBranchAddress("els_numlosthits", &els_numlosthits, &b_els_numlosthits);
   fChain->SetBranchAddress("els_tk_pz", &els_tk_pz, &b_els_tk_pz);
   fChain->SetBranchAddress("els_tk_pt", &els_tk_pt, &b_els_tk_pt);
   fChain->SetBranchAddress("els_tk_phi", &els_tk_phi, &b_els_tk_phi);
   fChain->SetBranchAddress("els_d0dum", &els_d0dum, &b_els_d0dum);
   fChain->SetBranchAddress("els_vx", &els_vx, &b_els_vx);
   fChain->SetBranchAddress("els_vy", &els_vy, &b_els_vy);
   fChain->SetBranchAddress("els_vz", &els_vz, &b_els_vz);
   fChain->SetBranchAddress("els_ptError", &els_ptError, &b_els_ptError);
   fChain->SetBranchAddress("els_n_inner_layer", &els_n_inner_layer, &b_els_n_inner_layer);
   fChain->SetBranchAddress("els_dr03EcalRecHitSumEt", &els_dr03EcalRecHitSumEt, &b_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr03HcalTowerSumEt", &els_dr03HcalTowerSumEt, &b_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth1TowerSumEt", &els_dr03HcalDepth1TowerSumEt, &b_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth2TowerSumEt", &els_dr03HcalDepth2TowerSumEt, &b_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr03TkSumPt", &els_dr03TkSumPt, &b_els_dr03TkSumPt);
   fChain->SetBranchAddress("jets_AK5PF_phi", &jets_AK5PF_phi, &b_jets_AK5PF_phi);
   fChain->SetBranchAddress("jets_AK5PF_pt", &jets_AK5PF_pt, &b_jets_AK5PF_pt);
   fChain->SetBranchAddress("jets_AK5PF_pz", &jets_AK5PF_pz, &b_jets_AK5PF_pz);
   fChain->SetBranchAddress("jets_AK5PF_px", &jets_AK5PF_px, &b_jets_AK5PF_px);
   fChain->SetBranchAddress("jets_AK5PF_py", &jets_AK5PF_py, &b_jets_AK5PF_py);
   fChain->SetBranchAddress("jets_AK5PF_eta", &jets_AK5PF_eta, &b_jets_AK5PF_eta);
   fChain->SetBranchAddress("jets_AK5PF_et", &jets_AK5PF_et, &b_jets_AK5PF_et);
   fChain->SetBranchAddress("jets_AK5PF_energy", &jets_AK5PF_energy, &b_jets_AK5PF_energy);
   fChain->SetBranchAddress("jets_AK5PF_parton_Id", &jets_AK5PF_parton_Id, &b_jets_AK5PF_parton_Id);
   fChain->SetBranchAddress("jets_AK5PF_parton_motherId", &jets_AK5PF_parton_motherId, &b_jets_AK5PF_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PF_gen_pt", &jets_AK5PF_gen_pt, &b_jets_AK5PF_gen_pt);
   fChain->SetBranchAddress("jets_AK5PF_gen_phi", &jets_AK5PF_gen_phi, &b_jets_AK5PF_gen_phi);
   fChain->SetBranchAddress("jets_AK5PF_partonFlavour", &jets_AK5PF_partonFlavour, &b_jets_AK5PF_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highPur", &jets_AK5PF_btag_TC_highPur, &b_jets_AK5PF_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highEff", &jets_AK5PF_btag_TC_highEff, &b_jets_AK5PF_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetProb", &jets_AK5PF_btag_jetProb, &b_jets_AK5PF_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetBProb", &jets_AK5PF_btag_jetBProb, &b_jets_AK5PF_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighPur", &jets_AK5PF_btag_secVertexHighPur, &b_jets_AK5PF_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighEff", &jets_AK5PF_btag_secVertexHighEff, &b_jets_AK5PF_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexCombined", &jets_AK5PF_btag_secVertexCombined, &b_jets_AK5PF_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PF_jetCharge", &jets_AK5PF_jetCharge, &b_jets_AK5PF_jetCharge);
   fChain->SetBranchAddress("jets_AK5PF_chgEmE", &jets_AK5PF_chgEmE, &b_jets_AK5PF_chgEmE);
   fChain->SetBranchAddress("jets_AK5PF_chgHadE", &jets_AK5PF_chgHadE, &b_jets_AK5PF_chgHadE);
   fChain->SetBranchAddress("jets_AK5PF_photonEnergy", &jets_AK5PF_photonEnergy, &b_jets_AK5PF_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PF_chg_Mult", &jets_AK5PF_chg_Mult, &b_jets_AK5PF_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PF_neutralEmE", &jets_AK5PF_neutralEmE, &b_jets_AK5PF_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PF_neutralHadE", &jets_AK5PF_neutralHadE, &b_jets_AK5PF_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PF_neutral_Mult", &jets_AK5PF_neutral_Mult, &b_jets_AK5PF_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PF_mu_Mult", &jets_AK5PF_mu_Mult, &b_jets_AK5PF_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PF_ehf", &jets_AK5PF_ehf, &b_jets_AK5PF_ehf);
   fChain->SetBranchAddress("jets_AK5PF_corrFactorRaw", &jets_AK5PF_corrFactorRaw, &b_jets_AK5PF_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK5PFclean_phi", &jets_AK5PFclean_phi, &b_jets_AK5PFclean_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_pt", &jets_AK5PFclean_pt, &b_jets_AK5PFclean_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_pz", &jets_AK5PFclean_pz, &b_jets_AK5PFclean_pz);
   fChain->SetBranchAddress("jets_AK5PFclean_px", &jets_AK5PFclean_px, &b_jets_AK5PFclean_px);
   fChain->SetBranchAddress("jets_AK5PFclean_py", &jets_AK5PFclean_py, &b_jets_AK5PFclean_py);
   fChain->SetBranchAddress("jets_AK5PFclean_eta", &jets_AK5PFclean_eta, &b_jets_AK5PFclean_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_et", &jets_AK5PFclean_et, &b_jets_AK5PFclean_et);
   fChain->SetBranchAddress("jets_AK5PFclean_energy", &jets_AK5PFclean_energy, &b_jets_AK5PFclean_energy);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_Id", &jets_AK5PFclean_parton_Id, &b_jets_AK5PFclean_parton_Id);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_motherId", &jets_AK5PFclean_parton_motherId, &b_jets_AK5PFclean_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_pt", &jets_AK5PFclean_gen_pt, &b_jets_AK5PFclean_gen_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_phi", &jets_AK5PFclean_gen_phi, &b_jets_AK5PFclean_gen_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_partonFlavour", &jets_AK5PFclean_partonFlavour, &b_jets_AK5PFclean_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highPur", &jets_AK5PFclean_btag_TC_highPur, &b_jets_AK5PFclean_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highEff", &jets_AK5PFclean_btag_TC_highEff, &b_jets_AK5PFclean_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetProb", &jets_AK5PFclean_btag_jetProb, &b_jets_AK5PFclean_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetBProb", &jets_AK5PFclean_btag_jetBProb, &b_jets_AK5PFclean_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighPur", &jets_AK5PFclean_btag_secVertexHighPur, &b_jets_AK5PFclean_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighEff", &jets_AK5PFclean_btag_secVertexHighEff, &b_jets_AK5PFclean_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexCombined", &jets_AK5PFclean_btag_secVertexCombined, &b_jets_AK5PFclean_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PFclean_jetCharge", &jets_AK5PFclean_jetCharge, &b_jets_AK5PFclean_jetCharge);
   fChain->SetBranchAddress("jets_AK5PFclean_chgEmE", &jets_AK5PFclean_chgEmE, &b_jets_AK5PFclean_chgEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_chgHadE", &jets_AK5PFclean_chgHadE, &b_jets_AK5PFclean_chgHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_photonEnergy", &jets_AK5PFclean_photonEnergy, &b_jets_AK5PFclean_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PFclean_chg_Mult", &jets_AK5PFclean_chg_Mult, &b_jets_AK5PFclean_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralEmE", &jets_AK5PFclean_neutralEmE, &b_jets_AK5PFclean_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralHadE", &jets_AK5PFclean_neutralHadE, &b_jets_AK5PFclean_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutral_Mult", &jets_AK5PFclean_neutral_Mult, &b_jets_AK5PFclean_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_mu_Mult", &jets_AK5PFclean_mu_Mult, &b_jets_AK5PFclean_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_ehf", &jets_AK5PFclean_ehf, &b_jets_AK5PFclean_ehf);
   fChain->SetBranchAddress("jets_AK5PFclean_corrFactorRaw", &jets_AK5PFclean_corrFactorRaw, &b_jets_AK5PFclean_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK5PFclean_rawPt", &jets_AK5PFclean_rawPt, &b_jets_AK5PFclean_rawPt);
   fChain->SetBranchAddress("mus_energy", &mus_energy, &b_mus_energy);
   fChain->SetBranchAddress("mus_et", &mus_et, &b_mus_et);
   fChain->SetBranchAddress("mus_eta", &mus_eta, &b_mus_eta);
   fChain->SetBranchAddress("mus_phi", &mus_phi, &b_mus_phi);
   fChain->SetBranchAddress("mus_pt", &mus_pt, &b_mus_pt);
   fChain->SetBranchAddress("mus_px", &mus_px, &b_mus_px);
   fChain->SetBranchAddress("mus_py", &mus_py, &b_mus_py);
   fChain->SetBranchAddress("mus_pz", &mus_pz, &b_mus_pz);
   fChain->SetBranchAddress("mus_cIso", &mus_cIso, &b_mus_cIso);
   fChain->SetBranchAddress("mus_tIso", &mus_tIso, &b_mus_tIso);
   fChain->SetBranchAddress("mus_ecalIso", &mus_ecalIso, &b_mus_ecalIso);
   fChain->SetBranchAddress("mus_hcalIso", &mus_hcalIso, &b_mus_hcalIso);
   fChain->SetBranchAddress("mus_ecalvetoDep", &mus_ecalvetoDep, &b_mus_ecalvetoDep);
   fChain->SetBranchAddress("mus_hcalvetoDep", &mus_hcalvetoDep, &b_mus_hcalvetoDep);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumChargedHadronPt", &mus_pfIsolationR03_sumChargedHadronPt, &b_mus_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumChargedParticlePt", &mus_pfIsolationR03_sumChargedParticlePt, &b_mus_pfIsolationR03_sumChargedParticlePt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumNeutralHadronEt", &mus_pfIsolationR03_sumNeutralHadronEt, &b_mus_pfIsolationR03_sumNeutralHadronEt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumNeutralHadronEtHighThreshold", &mus_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumPhotonEt", &mus_pfIsolationR03_sumPhotonEt, &b_mus_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumPhotonEtHighThreshold", &mus_pfIsolationR03_sumPhotonEtHighThreshold, &b_mus_pfIsolationR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR03_sumPUPt", &mus_pfIsolationR03_sumPUPt, &b_mus_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumChargedHadronPt", &mus_pfIsolationR04_sumChargedHadronPt, &b_mus_pfIsolationR04_sumChargedHadronPt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumChargedParticlePt", &mus_pfIsolationR04_sumChargedParticlePt, &b_mus_pfIsolationR04_sumChargedParticlePt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumNeutralHadronEt", &mus_pfIsolationR04_sumNeutralHadronEt, &b_mus_pfIsolationR04_sumNeutralHadronEt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumNeutralHadronEtHighThreshold", &mus_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumPhotonEt", &mus_pfIsolationR04_sumPhotonEt, &b_mus_pfIsolationR04_sumPhotonEt);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumPhotonEtHighThreshold", &mus_pfIsolationR04_sumPhotonEtHighThreshold, &b_mus_pfIsolationR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mus_pfIsolationR04_sumPUPt", &mus_pfIsolationR04_sumPUPt, &b_mus_pfIsolationR04_sumPUPt);
   fChain->SetBranchAddress("mus_charge", &mus_charge, &b_mus_charge);
   fChain->SetBranchAddress("mus_cm_chi2", &mus_cm_chi2, &b_mus_cm_chi2);
   fChain->SetBranchAddress("mus_cm_ndof", &mus_cm_ndof, &b_mus_cm_ndof);
   fChain->SetBranchAddress("mus_cm_pt", &mus_cm_pt, &b_mus_cm_pt);
   fChain->SetBranchAddress("mus_cm_ptErr", &mus_cm_ptErr, &b_mus_cm_ptErr);
   fChain->SetBranchAddress("mus_tk_chi2", &mus_tk_chi2, &b_mus_tk_chi2);
   fChain->SetBranchAddress("mus_tk_ndof", &mus_tk_ndof, &b_mus_tk_ndof);
   fChain->SetBranchAddress("mus_tk_pt", &mus_tk_pt, &b_mus_tk_pt);
   fChain->SetBranchAddress("mus_tk_px", &mus_tk_px, &b_mus_tk_px);
   fChain->SetBranchAddress("mus_tk_py", &mus_tk_py, &b_mus_tk_py);
   fChain->SetBranchAddress("mus_tk_pz", &mus_tk_pz, &b_mus_tk_pz);
   fChain->SetBranchAddress("mus_tk_phi", &mus_tk_phi, &b_mus_tk_phi);
   fChain->SetBranchAddress("mus_tk_d0dum", &mus_tk_d0dum, &b_mus_tk_d0dum);
   fChain->SetBranchAddress("mus_tk_vx", &mus_tk_vx, &b_mus_tk_vx);
   fChain->SetBranchAddress("mus_tk_vy", &mus_tk_vy, &b_mus_tk_vy);
   fChain->SetBranchAddress("mus_tk_vz", &mus_tk_vz, &b_mus_tk_vz);
   fChain->SetBranchAddress("mus_tk_numvalhits", &mus_tk_numvalhits, &b_mus_tk_numvalhits);
   fChain->SetBranchAddress("mus_tk_ptErr", &mus_tk_ptErr, &b_mus_tk_ptErr);
   fChain->SetBranchAddress("mus_tk_numvalPixelhits", &mus_tk_numvalPixelhits, &b_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("mus_tk_numpixelWthMeasr", &mus_tk_numpixelWthMeasr, &b_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("mus_stamu_pt", &mus_stamu_pt, &b_mus_stamu_pt);
   fChain->SetBranchAddress("mus_stamu_ptErr", &mus_stamu_ptErr, &b_mus_stamu_ptErr);
   fChain->SetBranchAddress("mus_num_matches", &mus_num_matches, &b_mus_num_matches);
   fChain->SetBranchAddress("mus_isPFMuon", &mus_isPFMuon, &b_mus_isPFMuon);
   fChain->SetBranchAddress("mus_isTrackerMuon", &mus_isTrackerMuon, &b_mus_isTrackerMuon);
   fChain->SetBranchAddress("mus_isGlobalMuon", &mus_isGlobalMuon, &b_mus_isGlobalMuon);
   fChain->SetBranchAddress("mus_id_AllGlobalMuons", &mus_id_AllGlobalMuons, &b_mus_id_AllGlobalMuons);
   fChain->SetBranchAddress("mus_id_AllTrackerMuons", &mus_id_AllTrackerMuons, &b_mus_id_AllTrackerMuons);
   fChain->SetBranchAddress("mus_id_GlobalMuonPromptTight", &mus_id_GlobalMuonPromptTight, &b_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("mus_tk_LayersWithMeasurement", &mus_tk_LayersWithMeasurement, &b_mus_tk_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_dB", &mus_dB, &b_mus_dB);
   fChain->SetBranchAddress("mus_numberOfMatchedStations", &mus_numberOfMatchedStations, &b_mus_numberOfMatchedStations);
   fChain->SetBranchAddress("pfTypeImets_et", &pfTypeImets_et, &b_pfTypeImets_et);
   fChain->SetBranchAddress("pfTypeImets_phi", &pfTypeImets_phi, &b_pfTypeImets_phi);
   fChain->SetBranchAddress("pfTypeImets_ex", &pfTypeImets_ex, &b_pfTypeImets_ex);
   fChain->SetBranchAddress("pfTypeImets_ey", &pfTypeImets_ey, &b_pfTypeImets_ey);
   fChain->SetBranchAddress("pfTypeImets_gen_et", &pfTypeImets_gen_et, &b_pfTypeImets_gen_et);
   fChain->SetBranchAddress("pfTypeImets_gen_phi", &pfTypeImets_gen_phi, &b_pfTypeImets_gen_phi);
   fChain->SetBranchAddress("pfTypeImets_sumEt", &pfTypeImets_sumEt, &b_pfTypeImets_sumEt);
   fChain->SetBranchAddress("pf_els_energy", &pf_els_energy, &b_pf_els_energy);
   fChain->SetBranchAddress("pf_els_et", &pf_els_et, &b_pf_els_et);
   fChain->SetBranchAddress("pf_els_eta", &pf_els_eta, &b_pf_els_eta);
   fChain->SetBranchAddress("pf_els_phi", &pf_els_phi, &b_pf_els_phi);
   fChain->SetBranchAddress("pf_els_pt", &pf_els_pt, &b_pf_els_pt);
   fChain->SetBranchAddress("pf_els_px", &pf_els_px, &b_pf_els_px);
   fChain->SetBranchAddress("pf_els_py", &pf_els_py, &b_pf_els_py);
   fChain->SetBranchAddress("pf_els_pz", &pf_els_pz, &b_pf_els_pz);
   fChain->SetBranchAddress("pf_els_robustTightId", &pf_els_robustTightId, &b_pf_els_robustTightId);
   fChain->SetBranchAddress("pf_els_simpleEleId95relIso", &pf_els_simpleEleId95relIso, &b_pf_els_simpleEleId95relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90relIso", &pf_els_simpleEleId90relIso, &b_pf_els_simpleEleId90relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85relIso", &pf_els_simpleEleId85relIso, &b_pf_els_simpleEleId85relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80relIso", &pf_els_simpleEleId80relIso, &b_pf_els_simpleEleId80relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70relIso", &pf_els_simpleEleId70relIso, &b_pf_els_simpleEleId70relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId95cIso", &pf_els_simpleEleId95cIso, &b_pf_els_simpleEleId95cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90cIso", &pf_els_simpleEleId90cIso, &b_pf_els_simpleEleId90cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85cIso", &pf_els_simpleEleId85cIso, &b_pf_els_simpleEleId85cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80cIso", &pf_els_simpleEleId80cIso, &b_pf_els_simpleEleId80cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70cIso", &pf_els_simpleEleId70cIso, &b_pf_els_simpleEleId70cIso);
   fChain->SetBranchAddress("pf_els_cIso", &pf_els_cIso, &b_pf_els_cIso);
   fChain->SetBranchAddress("pf_els_tIso", &pf_els_tIso, &b_pf_els_tIso);
   fChain->SetBranchAddress("pf_els_ecalIso", &pf_els_ecalIso, &b_pf_els_ecalIso);
   fChain->SetBranchAddress("pf_els_hcalIso", &pf_els_hcalIso, &b_pf_els_hcalIso);
   fChain->SetBranchAddress("pf_els_chargedHadronIso", &pf_els_chargedHadronIso, &b_pf_els_chargedHadronIso);
   fChain->SetBranchAddress("pf_els_photonIso", &pf_els_photonIso, &b_pf_els_photonIso);
   fChain->SetBranchAddress("pf_els_neutralHadronIso", &pf_els_neutralHadronIso, &b_pf_els_neutralHadronIso);
   fChain->SetBranchAddress("pf_els_charge", &pf_els_charge, &b_pf_els_charge);
   fChain->SetBranchAddress("pf_els_hadOverEm", &pf_els_hadOverEm, &b_pf_els_hadOverEm);
   fChain->SetBranchAddress("pf_els_eOverPIn", &pf_els_eOverPIn, &b_pf_els_eOverPIn);
   fChain->SetBranchAddress("pf_els_sigmaIEtaIEta", &pf_els_sigmaIEtaIEta, &b_pf_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("pf_els_scEnergy", &pf_els_scEnergy, &b_pf_els_scEnergy);
   fChain->SetBranchAddress("pf_els_scEta", &pf_els_scEta, &b_pf_els_scEta);
   fChain->SetBranchAddress("pf_els_scE1x5", &pf_els_scE1x5, &b_pf_els_scE1x5);
   fChain->SetBranchAddress("pf_els_scE2x5Max", &pf_els_scE2x5Max, &b_pf_els_scE2x5Max);
   fChain->SetBranchAddress("pf_els_scE5x5", &pf_els_scE5x5, &b_pf_els_scE5x5);
   fChain->SetBranchAddress("pf_els_isEB", &pf_els_isEB, &b_pf_els_isEB);
   fChain->SetBranchAddress("pf_els_isEE", &pf_els_isEE, &b_pf_els_isEE);
   fChain->SetBranchAddress("pf_els_dEtaIn", &pf_els_dEtaIn, &b_pf_els_dEtaIn);
   fChain->SetBranchAddress("pf_els_dPhiIn", &pf_els_dPhiIn, &b_pf_els_dPhiIn);
   fChain->SetBranchAddress("pf_els_dEtaOut", &pf_els_dEtaOut, &b_pf_els_dEtaOut);
   fChain->SetBranchAddress("pf_els_dPhiOut", &pf_els_dPhiOut, &b_pf_els_dPhiOut);
   fChain->SetBranchAddress("pf_els_numlosthits", &pf_els_numlosthits, &b_pf_els_numlosthits);
   fChain->SetBranchAddress("pf_els_tk_phi", &pf_els_tk_phi, &b_pf_els_tk_phi);
   fChain->SetBranchAddress("pf_els_d0dum", &pf_els_d0dum, &b_pf_els_d0dum);
   fChain->SetBranchAddress("pf_els_vx", &pf_els_vx, &b_pf_els_vx);
   fChain->SetBranchAddress("pf_els_vy", &pf_els_vy, &b_pf_els_vy);
   fChain->SetBranchAddress("pf_els_vz", &pf_els_vz, &b_pf_els_vz);
   fChain->SetBranchAddress("pf_els_ptError", &pf_els_ptError, &b_pf_els_ptError);
   fChain->SetBranchAddress("pf_els_n_inner_layer", &pf_els_n_inner_layer, &b_pf_els_n_inner_layer);
   fChain->SetBranchAddress("pf_els_dr03EcalRecHitSumEt", &pf_els_dr03EcalRecHitSumEt, &b_pf_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalTowerSumEt", &pf_els_dr03HcalTowerSumEt, &b_pf_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth1TowerSumEt", &pf_els_dr03HcalDepth1TowerSumEt, &b_pf_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth2TowerSumEt", &pf_els_dr03HcalDepth2TowerSumEt, &b_pf_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03TkSumPt", &pf_els_dr03TkSumPt, &b_pf_els_dr03TkSumPt);
   fChain->SetBranchAddress("pf_mus_energy", &pf_mus_energy, &b_pf_mus_energy);
   fChain->SetBranchAddress("pf_mus_et", &pf_mus_et, &b_pf_mus_et);
   fChain->SetBranchAddress("pf_mus_eta", &pf_mus_eta, &b_pf_mus_eta);
   fChain->SetBranchAddress("pf_mus_phi", &pf_mus_phi, &b_pf_mus_phi);
   fChain->SetBranchAddress("pf_mus_pt", &pf_mus_pt, &b_pf_mus_pt);
   fChain->SetBranchAddress("pf_mus_px", &pf_mus_px, &b_pf_mus_px);
   fChain->SetBranchAddress("pf_mus_py", &pf_mus_py, &b_pf_mus_py);
   fChain->SetBranchAddress("pf_mus_pz", &pf_mus_pz, &b_pf_mus_pz);
   fChain->SetBranchAddress("pf_mus_cIso", &pf_mus_cIso, &b_pf_mus_cIso);
   fChain->SetBranchAddress("pf_mus_tIso", &pf_mus_tIso, &b_pf_mus_tIso);
   fChain->SetBranchAddress("pf_mus_ecalIso", &pf_mus_ecalIso, &b_pf_mus_ecalIso);
   fChain->SetBranchAddress("pf_mus_hcalIso", &pf_mus_hcalIso, &b_pf_mus_hcalIso);
   fChain->SetBranchAddress("pf_mus_neutralHadronIso", &pf_mus_neutralHadronIso, &b_pf_mus_neutralHadronIso);
   fChain->SetBranchAddress("pf_mus_chargedHadronIso", &pf_mus_chargedHadronIso, &b_pf_mus_chargedHadronIso);
   fChain->SetBranchAddress("pf_mus_photonIso", &pf_mus_photonIso, &b_pf_mus_photonIso);
   fChain->SetBranchAddress("pf_mus_charge", &pf_mus_charge, &b_pf_mus_charge);
   fChain->SetBranchAddress("pf_mus_cm_chi2", &pf_mus_cm_chi2, &b_pf_mus_cm_chi2);
   fChain->SetBranchAddress("pf_mus_cm_ndof", &pf_mus_cm_ndof, &b_pf_mus_cm_ndof);
   fChain->SetBranchAddress("pf_mus_cm_pt", &pf_mus_cm_pt, &b_pf_mus_cm_pt);
   fChain->SetBranchAddress("pf_mus_cm_ptErr", &pf_mus_cm_ptErr, &b_pf_mus_cm_ptErr);
   fChain->SetBranchAddress("pf_mus_tk_chi2", &pf_mus_tk_chi2, &b_pf_mus_tk_chi2);
   fChain->SetBranchAddress("pf_mus_tk_ndof", &pf_mus_tk_ndof, &b_pf_mus_tk_ndof);
   fChain->SetBranchAddress("pf_mus_tk_pt", &pf_mus_tk_pt, &b_pf_mus_tk_pt);
   fChain->SetBranchAddress("pf_mus_tk_phi", &pf_mus_tk_phi, &b_pf_mus_tk_phi);
   fChain->SetBranchAddress("pf_mus_tk_d0dum", &pf_mus_tk_d0dum, &b_pf_mus_tk_d0dum);
   fChain->SetBranchAddress("pf_mus_tk_vz", &pf_mus_tk_vz, &b_pf_mus_tk_vz);
   fChain->SetBranchAddress("pf_mus_tk_numvalhits", &pf_mus_tk_numvalhits, &b_pf_mus_tk_numvalhits);
   fChain->SetBranchAddress("pf_mus_tk_ptErr", &pf_mus_tk_ptErr, &b_pf_mus_tk_ptErr);
   fChain->SetBranchAddress("pf_mus_tk_numvalPixelhits", &pf_mus_tk_numvalPixelhits, &b_pf_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_tk_numpixelWthMeasr", &pf_mus_tk_numpixelWthMeasr, &b_pf_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("pf_mus_stamu_pt", &pf_mus_stamu_pt, &b_pf_mus_stamu_pt);
   fChain->SetBranchAddress("pf_mus_stamu_ptErr", &pf_mus_stamu_ptErr, &b_pf_mus_stamu_ptErr);
   fChain->SetBranchAddress("pf_mus_num_matches", &pf_mus_num_matches, &b_pf_mus_num_matches);
   fChain->SetBranchAddress("pf_mus_isTrackerMuon", &pf_mus_isTrackerMuon, &b_pf_mus_isTrackerMuon);
   fChain->SetBranchAddress("pf_mus_isGlobalMuon", &pf_mus_isGlobalMuon, &b_pf_mus_isGlobalMuon);
   fChain->SetBranchAddress("pf_mus_id_GlobalMuonPromptTight", &pf_mus_id_GlobalMuonPromptTight, &b_pf_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("pfcand_pdgId", &pfcand_pdgId, &b_pfcand_pdgId);
   fChain->SetBranchAddress("pfcand_particleId", &pfcand_particleId, &b_pfcand_particleId);
   fChain->SetBranchAddress("pfcand_pt", &pfcand_pt, &b_pfcand_pt);
   fChain->SetBranchAddress("pfcand_pz", &pfcand_pz, &b_pfcand_pz);
   fChain->SetBranchAddress("pfcand_px", &pfcand_px, &b_pfcand_px);
   fChain->SetBranchAddress("pfcand_py", &pfcand_py, &b_pfcand_py);
   fChain->SetBranchAddress("pfcand_eta", &pfcand_eta, &b_pfcand_eta);
   fChain->SetBranchAddress("pfcand_phi", &pfcand_phi, &b_pfcand_phi);
   fChain->SetBranchAddress("pfcand_theta", &pfcand_theta, &b_pfcand_theta);
   fChain->SetBranchAddress("pfcand_energy", &pfcand_energy, &b_pfcand_energy);
   fChain->SetBranchAddress("pfcand_charge", &pfcand_charge, &b_pfcand_charge);
   fChain->SetBranchAddress("pfmets_et", &pfmets_et, &b_pfmets_et);
   fChain->SetBranchAddress("pfmets_phi", &pfmets_phi, &b_pfmets_phi);
   fChain->SetBranchAddress("pfmets_ex", &pfmets_ex, &b_pfmets_ex);
   fChain->SetBranchAddress("pfmets_ey", &pfmets_ey, &b_pfmets_ey);
   fChain->SetBranchAddress("pfmets_gen_et", &pfmets_gen_et, &b_pfmets_gen_et);
   fChain->SetBranchAddress("pfmets_gen_phi", &pfmets_gen_phi, &b_pfmets_gen_phi);
   fChain->SetBranchAddress("pfmets_sumEt", &pfmets_sumEt, &b_pfmets_sumEt);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("pv_tracksSize", &pv_tracksSize, &b_pv_tracksSize);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("model_params", &model_params, &b_model_params);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_px", &fastjets_AK5PFclean_R1p2_R0p5pT10_px, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_px);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_py", &fastjets_AK5PFclean_R1p2_R0p5pT10_py, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_py);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_pz", &fastjets_AK5PFclean_R1p2_R0p5pT10_pz, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_pz);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_energy", &fastjets_AK5PFclean_R1p2_R0p5pT10_energy, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_energy);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_phi", &fastjets_AK5PFclean_R1p2_R0p5pT10_phi, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_phi);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_eta", &fastjets_AK5PFclean_R1p2_R0p5pT10_eta, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_eta);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_index", &fastjets_AK5PFclean_R1p2_R0p5pT10_index, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_index);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT10_nconstituents", &fastjets_AK5PFclean_R1p2_R0p5pT10_nconstituents, &b_fastjets_AK5PFclean_R1p2_R0p5pT10_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_px", &fastjets_AK5PFclean_R1p2_R0p5pT15_px, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_px);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_py", &fastjets_AK5PFclean_R1p2_R0p5pT15_py, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_py);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_pz", &fastjets_AK5PFclean_R1p2_R0p5pT15_pz, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_pz);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_energy", &fastjets_AK5PFclean_R1p2_R0p5pT15_energy, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_energy);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_phi", &fastjets_AK5PFclean_R1p2_R0p5pT15_phi, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_phi);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_eta", &fastjets_AK5PFclean_R1p2_R0p5pT15_eta, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_eta);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_index", &fastjets_AK5PFclean_R1p2_R0p5pT15_index, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_index);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT15_nconstituents", &fastjets_AK5PFclean_R1p2_R0p5pT15_nconstituents, &b_fastjets_AK5PFclean_R1p2_R0p5pT15_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_px", &fastjets_AK5PFclean_R1p2_R0p5pT20_px, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_px);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_py", &fastjets_AK5PFclean_R1p2_R0p5pT20_py, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_py);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_pz", &fastjets_AK5PFclean_R1p2_R0p5pT20_pz, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_pz);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_energy", &fastjets_AK5PFclean_R1p2_R0p5pT20_energy, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_energy);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_phi", &fastjets_AK5PFclean_R1p2_R0p5pT20_phi, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_phi);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_eta", &fastjets_AK5PFclean_R1p2_R0p5pT20_eta, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_eta);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_index", &fastjets_AK5PFclean_R1p2_R0p5pT20_index, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_index);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT20_nconstituents", &fastjets_AK5PFclean_R1p2_R0p5pT20_nconstituents, &b_fastjets_AK5PFclean_R1p2_R0p5pT20_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_px", &fastjets_AK5PFclean_R1p2_R0p5pT25_px, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_px);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_py", &fastjets_AK5PFclean_R1p2_R0p5pT25_py, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_py);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_pz", &fastjets_AK5PFclean_R1p2_R0p5pT25_pz, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_pz);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_energy", &fastjets_AK5PFclean_R1p2_R0p5pT25_energy, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_energy);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_phi", &fastjets_AK5PFclean_R1p2_R0p5pT25_phi, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_phi);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_eta", &fastjets_AK5PFclean_R1p2_R0p5pT25_eta, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_eta);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_index", &fastjets_AK5PFclean_R1p2_R0p5pT25_index, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_index);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT25_nconstituents", &fastjets_AK5PFclean_R1p2_R0p5pT25_nconstituents, &b_fastjets_AK5PFclean_R1p2_R0p5pT25_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_px", &fastjets_AK5PFclean_R1p2_R0p5pT30_px, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_px);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_py", &fastjets_AK5PFclean_R1p2_R0p5pT30_py, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_py);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_pz", &fastjets_AK5PFclean_R1p2_R0p5pT30_pz, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_pz);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_energy", &fastjets_AK5PFclean_R1p2_R0p5pT30_energy, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_energy);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_phi", &fastjets_AK5PFclean_R1p2_R0p5pT30_phi, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_phi);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_eta", &fastjets_AK5PFclean_R1p2_R0p5pT30_eta, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_eta);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_index", &fastjets_AK5PFclean_R1p2_R0p5pT30_index, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_index);
   fChain->SetBranchAddress("fastjets_AK5PFclean_R1p2_R0p5pT30_nconstituents", &fastjets_AK5PFclean_R1p2_R0p5pT30_nconstituents, &b_fastjets_AK5PFclean_R1p2_R0p5pT30_nconstituents);

}


void InitializeA(TChain *fChain)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   els_conversion_dist = 0;
   els_conversion_dcot = 0;
   els_PFchargedHadronIsoR03 = 0;
   els_PFphotonIsoR03 = 0;
   els_PFneutralHadronIsoR03 = 0;
   els_hasMatchedConversion = 0;
   PU_TrueNumInteractions = 0;

   fChain->SetBranchAddress("els_conversion_dist", &els_conversion_dist, &b_els_conversion_dist);
   fChain->SetBranchAddress("els_conversion_dcot", &els_conversion_dcot, &b_els_conversion_dcot);
   fChain->SetBranchAddress("els_PFchargedHadronIsoR03", &els_PFchargedHadronIsoR03, &b_els_PFchargedHadronIsoR03);
   fChain->SetBranchAddress("els_PFphotonIsoR03", &els_PFphotonIsoR03, &b_els_PFphotonIsoR03);
   fChain->SetBranchAddress("els_PFneutralHadronIsoR03", &els_PFneutralHadronIsoR03, &b_els_PFneutralHadronIsoR03);
   fChain->SetBranchAddress("els_hasMatchedConversion", &els_hasMatchedConversion, &b_els_hasMatchedConversion);
   fChain->SetBranchAddress("PU_TrueNumInteractions", &PU_TrueNumInteractions, &b_PU_TrueNumInteractions);
   fChain->SetBranchAddress("rho_kt6PFJetsForIsolation2011", &rho_kt6PFJetsForIsolation2011, &b_rho_kt6PFJetsForIsolation2011);

}
