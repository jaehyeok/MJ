#include <vector>
#include <string>
#include "TTree.h"
#include "TChain.h"

   // Declaration of leaf types
   vector<bool>    *trigger_decision;
   vector<string>  *trigger_name;
   vector<float>   *trigger_prescalevalue;
   vector<float>   *standalone_triggerobject_pt;
   vector<float>   *standalone_triggerobject_px;
   vector<float>   *standalone_triggerobject_py;
   vector<float>   *standalone_triggerobject_pz;
   vector<float>   *standalone_triggerobject_et;
   vector<float>   *standalone_triggerobject_energy;
   vector<float>   *standalone_triggerobject_phi;
   vector<float>   *standalone_triggerobject_eta;
   vector<string>  *standalone_triggerobject_collectionname;
   vector<vector<float> > *PU_zpositions;
   vector<vector<float> > *PU_sumpT_lowpT;
   vector<vector<float> > *PU_sumpT_highpT;
   vector<vector<int> > *PU_ntrks_lowpT;
   vector<vector<int> > *PU_ntrks_highpT;
   vector<int>     *PU_NumInteractions;
   vector<int>     *PU_bunchCrossing;
   vector<float>   *PU_TrueNumInteractions;
   Int_t           trackingfailurefilter_decision;
   Int_t           goodVerticesfilter_decision;
   Int_t           cschalofilter_decision;
   Int_t           trkPOGfilter_decision;
   Int_t           trkPOG_logErrorTooManyClustersfilter_decision;
   Int_t           EcalDeadCellTriggerPrimitivefilter_decision;
   Int_t           ecallaserfilter_decision;
   Int_t           trkPOG_manystripclus53Xfilter_decision;
   Int_t           eebadscfilter_decision;
   Int_t           METFiltersfilter_decision;
   Int_t           HBHENoisefilter_decision;
   Int_t           trkPOG_toomanystripclus53Xfilter_decision;
   Int_t           hcallaserfilter_decision;
   vector<bool>    *els_isPF;
   vector<bool>    *mus_isPF;
   vector<int>     *jets_AK4_maxpt_id;
   vector<int>     *jets_AK4_mu_ind;
   vector<int>     *jets_AK4_el_ind;
   vector<int>     *taus_el_ind;
   vector<int>     *taus_mu_ind;
   vector<int>     *els_jet_ind;
   vector<int>     *mus_jet_ind;

   // List of branches
   TBranch        *b_trigger_decision;   //!
   TBranch        *b_trigger_name;   //!
   TBranch        *b_trigger_prescalevalue;   //!
   TBranch        *b_standalone_triggerobject_pt;   //!
   TBranch        *b_standalone_triggerobject_px;   //!
   TBranch        *b_standalone_triggerobject_py;   //!
   TBranch        *b_standalone_triggerobject_pz;   //!
   TBranch        *b_standalone_triggerobject_et;   //!
   TBranch        *b_standalone_triggerobject_energy;   //!
   TBranch        *b_standalone_triggerobject_phi;   //!
   TBranch        *b_standalone_triggerobject_eta;   //!
   TBranch        *b_standalone_triggerobject_collectionname;   //!
   TBranch        *b_PU_zpositions;   //!
   TBranch        *b_PU_sumpT_lowpT;   //!
   TBranch        *b_PU_sumpT_highpT;   //!
   TBranch        *b_PU_ntrks_lowpT;   //!
   TBranch        *b_PU_ntrks_highpT;   //!
   TBranch        *b_PU_NumInteractions;   //!
   TBranch        *b_PU_bunchCrossing;   //!
   TBranch        *b_PU_TrueNumInteractions;   //!
   TBranch        *b_rackingfailurefilter_decision;   //!
   TBranch        *b_oodVerticesfilter_decision;   //!
   TBranch        *b_schalofilter_decision;   //!
   TBranch        *b_rkPOGfilter_decision;   //!
   TBranch        *b_rkPOG_logErrorTooManyClustersfilter_decision;   //!
   TBranch        *b_calDeadCellTriggerPrimitivefilter_decision;   //!
   TBranch        *b_ecallaserfilter_decision;   //!
   TBranch        *b_rkPOG_manystripclus53Xfilter_decision;   //!
   TBranch        *b_ebadscfilter_decision;   //!
   TBranch        *b_ETFiltersfilter_decision;   //!
   TBranch        *b_BHENoisefilter_decision;   //!
   TBranch        *b_rkPOG_toomanystripclus53Xfilter_decision;   //!
   TBranch        *b_hcallaserfilter_decision;   //!
   TBranch        *b_els_isPF;   //!
   TBranch        *b_mus_isPF;   //!
   TBranch        *b_jets_AK4_maxpt_id;   //!
   TBranch        *b_jets_AK4_mu_ind;   //!
   TBranch        *b_jets_AK4_el_ind;   //!
   TBranch        *b_taus_el_ind;   //!
   TBranch        *b_taus_mu_ind;   //!
   TBranch        *b_els_jet_ind;   //!
   TBranch        *b_mus_jet_ind;   //!

   // Declaration of leaf types
   UInt_t          NbeamSpot;
   vector<float>   *beamSpot_x;
   vector<float>   *beamSpot_y;
   vector<float>   *beamSpot_z;
   vector<float>   *beamSpot_x0Error;
   vector<float>   *beamSpot_y0Error;
   vector<float>   *beamSpot_z0Error;
   vector<float>   *beamSpot_sigmaZ;
   vector<float>   *beamSpot_sigmaZ0Error;
   vector<float>   *beamSpot_dxdz;
   vector<float>   *beamSpot_dxdzError;
   vector<float>   *beamSpot_dydz;
   vector<float>   *beamSpot_dydzError;
   vector<float>   *beamSpot_beamWidthX;
   vector<float>   *beamSpot_beamWidthY;
   vector<float>   *beamSpot_beamWidthXError;
   vector<float>   *beamSpot_beamWidthYError;
   UInt_t          Nels;
   vector<float>   *els_energy;
   vector<float>   *els_et;
   vector<float>   *els_eta;
   vector<float>   *els_phi;
   vector<float>   *els_pt;
   vector<float>   *els_px;
   vector<float>   *els_py;
   vector<float>   *els_pz;
   vector<float>   *els_status;
   vector<float>   *els_theta;
   vector<float>   *els_pfIsolationR03_sumChargedHadronPt;
   vector<float>   *els_pfIsolationR03_sumNeutralHadronEt;
   vector<float>   *els_pfIsolationR03_sumPhotonEt;
   vector<float>   *els_pfIsolationR03_sumPUPt;
   vector<float>   *els_full5x5_sigmaIetaIeta;
   vector<float>   *els_gen_id;
   vector<float>   *els_gen_phi;
   vector<float>   *els_gen_pt;
   vector<float>   *els_gen_pz;
   vector<float>   *els_gen_px;
   vector<float>   *els_gen_py;
   vector<float>   *els_gen_eta;
   vector<float>   *els_gen_theta;
   vector<float>   *els_gen_et;
   vector<float>   *els_gen_mother_id;
   vector<float>   *els_gen_mother_phi;
   vector<float>   *els_gen_mother_pt;
   vector<float>   *els_gen_mother_pz;
   vector<float>   *els_gen_mother_px;
   vector<float>   *els_gen_mother_py;
   vector<float>   *els_gen_mother_eta;
   vector<float>   *els_gen_mother_theta;
   vector<float>   *els_gen_mother_et;
   vector<float>   *els_tightId;
   vector<float>   *els_looseId;
   vector<float>   *els_robustTightId;
   vector<float>   *els_robustLooseId;
   vector<float>   *els_robustHighEnergyId;
   vector<float>   *els_simpleEleId95relIso;
   vector<float>   *els_simpleEleId90relIso;
   vector<float>   *els_simpleEleId85relIso;
   vector<float>   *els_simpleEleId80relIso;
   vector<float>   *els_simpleEleId70relIso;
   vector<float>   *els_simpleEleId60relIso;
   vector<float>   *els_simpleEleId95cIso;
   vector<float>   *els_simpleEleId90cIso;
   vector<float>   *els_simpleEleId85cIso;
   vector<float>   *els_simpleEleId80cIso;
   vector<float>   *els_simpleEleId70cIso;
   vector<float>   *els_simpleEleId60cIso;
   vector<float>   *els_cIso;
   vector<float>   *els_tIso;
   vector<float>   *els_ecalIso;
   vector<float>   *els_hcalIso;
   vector<float>   *els_chi2;
   vector<float>   *els_charge;
   vector<float>   *els_caloEnergy;
   vector<float>   *els_hadOverEm;
   vector<float>   *els_hcalOverEcalBc;
   vector<float>   *els_eOverPIn;
   vector<float>   *els_eSeedOverPOut;
   vector<float>   *els_sigmaEtaEta;
   vector<float>   *els_sigmaIEtaIEta;
   vector<float>   *els_scEnergy;
   vector<float>   *els_scRawEnergy;
   vector<float>   *els_scSeedEnergy;
   vector<float>   *els_scEta;
   vector<float>   *els_scPhi;
   vector<float>   *els_scEtaWidth;
   vector<float>   *els_scPhiWidth;
   vector<float>   *els_scE1x5;
   vector<float>   *els_scE2x5Max;
   vector<float>   *els_scE5x5;
   vector<float>   *els_isEB;
   vector<float>   *els_isEE;
   vector<float>   *els_dEtaIn;
   vector<float>   *els_dPhiIn;
   vector<float>   *els_dEtaOut;
   vector<float>   *els_dPhiOut;
   vector<float>   *els_numvalhits;
   vector<float>   *els_numlosthits;
   vector<float>   *els_basicClustersSize;
   vector<float>   *els_tk_pz;
   vector<float>   *els_tk_pt;
   vector<float>   *els_tk_phi;
   vector<float>   *els_tk_eta;
   vector<float>   *els_d0dum;
   vector<float>   *els_dz;
   vector<float>   *els_vx;
   vector<float>   *els_vy;
   vector<float>   *els_vz;
   vector<float>   *els_ndof;
   vector<float>   *els_ptError;
   vector<float>   *els_d0dumError;
   vector<float>   *els_dzError;
   vector<float>   *els_etaError;
   vector<float>   *els_phiError;
   vector<float>   *els_tk_charge;
   vector<float>   *els_core_ecalDrivenSeed;
   vector<float>   *els_n_inner_layer;
   vector<float>   *els_n_outer_layer;
   vector<float>   *els_ctf_tk_id;
   vector<float>   *els_ctf_tk_charge;
   vector<float>   *els_ctf_tk_eta;
   vector<float>   *els_ctf_tk_phi;
   vector<float>   *els_fbrem;
   vector<float>   *els_shFracInnerHits;
   vector<float>   *els_dr03EcalRecHitSumEt;
   vector<float>   *els_dr03HcalTowerSumEt;
   vector<float>   *els_dr03HcalDepth1TowerSumEt;
   vector<float>   *els_dr03HcalDepth2TowerSumEt;
   vector<float>   *els_dr03TkSumPt;
   vector<float>   *els_dr04EcalRecHitSumEt;
   vector<float>   *els_dr04HcalTowerSumEt;
   vector<float>   *els_dr04HcalDepth1TowerSumEt;
   vector<float>   *els_dr04HcalDepth2TowerSumEt;
   vector<float>   *els_dr04TkSumPt;
   vector<float>   *els_cpx;
   vector<float>   *els_cpy;
   vector<float>   *els_cpz;
   vector<float>   *els_vpx;
   vector<float>   *els_vpy;
   vector<float>   *els_vpz;
   vector<float>   *els_cx;
   vector<float>   *els_cy;
   vector<float>   *els_cz;
   vector<float>   *els_PATpassConversionVeto;
   UInt_t          Njets_AK4;
   vector<float>   *jets_AK4_status;
   vector<float>   *jets_AK4_phi;
   vector<float>   *jets_AK4_pt;
   vector<float>   *jets_AK4_pz;
   vector<float>   *jets_AK4_px;
   vector<float>   *jets_AK4_py;
   vector<float>   *jets_AK4_eta;
   vector<float>   *jets_AK4_theta;
   vector<float>   *jets_AK4_et;
   vector<float>   *jets_AK4_energy;
   vector<float>   *jets_AK4_parton_Id;
   vector<float>   *jets_AK4_parton_motherId;
   vector<float>   *jets_AK4_parton_pt;
   vector<float>   *jets_AK4_parton_phi;
   vector<float>   *jets_AK4_parton_eta;
   vector<float>   *jets_AK4_parton_Energy;
   vector<float>   *jets_AK4_parton_mass;
   vector<float>   *jets_AK4_gen_et;
   vector<float>   *jets_AK4_gen_pt;
   vector<float>   *jets_AK4_gen_eta;
   vector<float>   *jets_AK4_gen_phi;
   vector<float>   *jets_AK4_gen_mass;
   vector<float>   *jets_AK4_gen_Energy;
   vector<float>   *jets_AK4_gen_Id;
   vector<float>   *jets_AK4_gen_motherID;
   vector<float>   *jets_AK4_gen_threeCharge;
   vector<float>   *jets_AK4_partonFlavour;
   vector<float>   *jets_AK4_btag_TC_highPur;
   vector<float>   *jets_AK4_btag_TC_highEff;
   vector<float>   *jets_AK4_btag_jetProb;
   vector<float>   *jets_AK4_btag_jetBProb;
   vector<float>   *jets_AK4_btag_softEle;
   vector<float>   *jets_AK4_btag_softMuon;
   vector<float>   *jets_AK4_btag_secVertexHighPur;
   vector<float>   *jets_AK4_btag_secVertexHighEff;
   vector<float>   *jets_AK4_btag_secVertexCombined;
   vector<float>   *jets_AK4_jetCharge;
   vector<float>   *jets_AK4_chgEmE;
   vector<float>   *jets_AK4_chgHadE;
   vector<float>   *jets_AK4_photonEnergy;
   vector<float>   *jets_AK4_chgMuE;
   vector<float>   *jets_AK4_chg_Mult;
   vector<float>   *jets_AK4_neutralEmE;
   vector<float>   *jets_AK4_neutralHadE;
   vector<float>   *jets_AK4_neutral_Mult;
   vector<float>   *jets_AK4_mu_Mult;
   vector<float>   *jets_AK4_emf;
   vector<float>   *jets_AK4_ehf;
   vector<float>   *jets_AK4_n60;
   vector<float>   *jets_AK4_n90;
   vector<float>   *jets_AK4_etaetaMoment;
   vector<float>   *jets_AK4_etaphiMoment;
   vector<float>   *jets_AK4_phiphiMoment;
   vector<float>   *jets_AK4_n90Hits;
   vector<float>   *jets_AK4_fHPD;
   vector<float>   *jets_AK4_fRBX;
   vector<float>   *jets_AK4_hitsInN90;
   vector<float>   *jets_AK4_nECALTowers;
   vector<float>   *jets_AK4_nHCALTowers;
   vector<float>   *jets_AK4_fSubDetector1;
   vector<float>   *jets_AK4_fSubDetector2;
   vector<float>   *jets_AK4_fSubDetector3;
   vector<float>   *jets_AK4_fSubDetector4;
   vector<float>   *jets_AK4_area;
   vector<float>   *jets_AK4_corrFactorRaw;
   vector<float>   *jets_AK4_rawPt;
   vector<float>   *jets_AK4_mass;
   UInt_t          Nmc_doc;
   vector<float>   *mc_doc_id;
   vector<float>   *mc_doc_pt;
   vector<float>   *mc_doc_px;
   vector<float>   *mc_doc_py;
   vector<float>   *mc_doc_pz;
   vector<float>   *mc_doc_eta;
   vector<float>   *mc_doc_phi;
   vector<float>   *mc_doc_theta;
   vector<float>   *mc_doc_energy;
   vector<float>   *mc_doc_status;
   vector<float>   *mc_doc_charge;
   vector<float>   *mc_doc_mother_id;
   vector<float>   *mc_doc_grandmother_id;
   vector<float>   *mc_doc_ggrandmother_id;
   vector<float>   *mc_doc_mother_pt;
   vector<float>   *mc_doc_vertex_x;
   vector<float>   *mc_doc_vertex_y;
   vector<float>   *mc_doc_vertex_z;
   vector<float>   *mc_doc_mass;
   vector<float>   *mc_doc_numOfDaughters;
   vector<float>   *mc_doc_numOfMothers;
   UInt_t          Nmc_electrons;
   vector<float>   *mc_electrons_id;
   vector<float>   *mc_electrons_pt;
   vector<float>   *mc_electrons_px;
   vector<float>   *mc_electrons_py;
   vector<float>   *mc_electrons_pz;
   vector<float>   *mc_electrons_eta;
   vector<float>   *mc_electrons_phi;
   vector<float>   *mc_electrons_theta;
   vector<float>   *mc_electrons_status;
   vector<float>   *mc_electrons_energy;
   vector<float>   *mc_electrons_charge;
   vector<float>   *mc_electrons_mother_id;
   vector<float>   *mc_electrons_mother_pt;
   vector<float>   *mc_electrons_grandmother_id;
   vector<float>   *mc_electrons_ggrandmother_id;
   vector<float>   *mc_electrons_vertex_x;
   vector<float>   *mc_electrons_vertex_y;
   vector<float>   *mc_electrons_vertex_z;
   vector<float>   *mc_electrons_mass;
   vector<float>   *mc_electrons_numOfDaughters;
   UInt_t          Nmc_final;
   vector<float>   *mc_final_id;
   vector<float>   *mc_final_pt;
   vector<float>   *mc_final_px;
   vector<float>   *mc_final_py;
   vector<float>   *mc_final_pz;
   vector<float>   *mc_final_eta;
   vector<float>   *mc_final_phi;
   vector<float>   *mc_final_theta;
   vector<float>   *mc_final_energy;
   vector<float>   *mc_final_status;
   vector<float>   *mc_final_charge;
   vector<float>   *mc_final_mother_id;
   vector<float>   *mc_final_grandmother_id;
   vector<float>   *mc_final_ggrandmother_id;
   vector<float>   *mc_final_mother_pt;
   vector<float>   *mc_final_vertex_x;
   vector<float>   *mc_final_vertex_y;
   vector<float>   *mc_final_vertex_z;
   vector<float>   *mc_final_mass;
   vector<float>   *mc_final_numOfDaughters;
   vector<float>   *mc_final_numOfMothers;
   UInt_t          Nmc_jets;
   vector<float>   *mc_jets_phi;
   vector<float>   *mc_jets_pt;
   vector<float>   *mc_jets_pz;
   vector<float>   *mc_jets_px;
   vector<float>   *mc_jets_py;
   vector<float>   *mc_jets_eta;
   vector<float>   *mc_jets_theta;
   vector<float>   *mc_jets_et;
   vector<float>   *mc_jets_energy;
   vector<float>   *mc_jets_emEnergy;
   vector<float>   *mc_jets_hadEnergy;
   vector<float>   *mc_jets_invisibleEnergy;
   vector<float>   *mc_jets_auxiliaryEnergy;
   vector<float>   *mc_jets_etaetaMoment;
   vector<float>   *mc_jets_etaphiMoment;
   vector<float>   *mc_jets_phiphiMoment;
   vector<float>   *mc_jets_mass;
   UInt_t          Nmc_mus;
   vector<float>   *mc_mus_id;
   vector<float>   *mc_mus_pt;
   vector<float>   *mc_mus_px;
   vector<float>   *mc_mus_py;
   vector<float>   *mc_mus_pz;
   vector<float>   *mc_mus_eta;
   vector<float>   *mc_mus_phi;
   vector<float>   *mc_mus_theta;
   vector<float>   *mc_mus_status;
   vector<float>   *mc_mus_energy;
   vector<float>   *mc_mus_charge;
   vector<float>   *mc_mus_mother_id;
   vector<float>   *mc_mus_mother_pt;
   vector<float>   *mc_mus_grandmother_id;
   vector<float>   *mc_mus_ggrandmother_id;
   vector<float>   *mc_mus_vertex_x;
   vector<float>   *mc_mus_vertex_y;
   vector<float>   *mc_mus_vertex_z;
   vector<float>   *mc_mus_mass;
   vector<float>   *mc_mus_numOfDaughters;
   UInt_t          Nmc_nues;
   vector<float>   *mc_nues_id;
   vector<float>   *mc_nues_pt;
   vector<float>   *mc_nues_px;
   vector<float>   *mc_nues_py;
   vector<float>   *mc_nues_pz;
   vector<float>   *mc_nues_eta;
   vector<float>   *mc_nues_phi;
   vector<float>   *mc_nues_theta;
   vector<float>   *mc_nues_status;
   vector<float>   *mc_nues_energy;
   vector<float>   *mc_nues_charge;
   vector<float>   *mc_nues_mother_id;
   vector<float>   *mc_nues_mother_pt;
   vector<float>   *mc_nues_grandmother_id;
   vector<float>   *mc_nues_ggrandmother_id;
   vector<float>   *mc_nues_vertex_x;
   vector<float>   *mc_nues_vertex_y;
   vector<float>   *mc_nues_vertex_z;
   vector<float>   *mc_nues_mass;
   vector<float>   *mc_nues_numOfDaughters;
   UInt_t          Nmc_numus;
   vector<float>   *mc_numus_id;
   vector<float>   *mc_numus_pt;
   vector<float>   *mc_numus_px;
   vector<float>   *mc_numus_py;
   vector<float>   *mc_numus_pz;
   vector<float>   *mc_numus_eta;
   vector<float>   *mc_numus_phi;
   vector<float>   *mc_numus_theta;
   vector<float>   *mc_numus_status;
   vector<float>   *mc_numus_energy;
   vector<float>   *mc_numus_charge;
   vector<float>   *mc_numus_mother_id;
   vector<float>   *mc_numus_mother_pt;
   vector<float>   *mc_numus_grandmother_id;
   vector<float>   *mc_numus_ggrandmother_id;
   vector<float>   *mc_numus_vertex_x;
   vector<float>   *mc_numus_vertex_y;
   vector<float>   *mc_numus_vertex_z;
   vector<float>   *mc_numus_mass;
   vector<float>   *mc_numus_numOfDaughters;
   UInt_t          Nmc_nutaus;
   vector<float>   *mc_nutaus_id;
   vector<float>   *mc_nutaus_pt;
   vector<float>   *mc_nutaus_px;
   vector<float>   *mc_nutaus_py;
   vector<float>   *mc_nutaus_pz;
   vector<float>   *mc_nutaus_eta;
   vector<float>   *mc_nutaus_phi;
   vector<float>   *mc_nutaus_theta;
   vector<float>   *mc_nutaus_status;
   vector<float>   *mc_nutaus_energy;
   vector<float>   *mc_nutaus_charge;
   vector<float>   *mc_nutaus_mother_id;
   vector<float>   *mc_nutaus_mother_pt;
   vector<float>   *mc_nutaus_grandmother_id;
   vector<float>   *mc_nutaus_ggrandmother_id;
   vector<float>   *mc_nutaus_vertex_x;
   vector<float>   *mc_nutaus_vertex_y;
   vector<float>   *mc_nutaus_vertex_z;
   vector<float>   *mc_nutaus_mass;
   vector<float>   *mc_nutaus_numOfDaughters;
   UInt_t          Nmc_photons;
   vector<float>   *mc_photons_id;
   vector<float>   *mc_photons_pt;
   vector<float>   *mc_photons_px;
   vector<float>   *mc_photons_py;
   vector<float>   *mc_photons_pz;
   vector<float>   *mc_photons_eta;
   vector<float>   *mc_photons_phi;
   vector<float>   *mc_photons_theta;
   vector<float>   *mc_photons_status;
   vector<float>   *mc_photons_energy;
   vector<float>   *mc_photons_charge;
   vector<float>   *mc_photons_mother_id;
   vector<float>   *mc_photons_mother_pt;
   vector<float>   *mc_photons_grandmother_id;
   vector<float>   *mc_photons_ggrandmother_id;
   vector<float>   *mc_photons_vertex_x;
   vector<float>   *mc_photons_vertex_y;
   vector<float>   *mc_photons_vertex_z;
   vector<float>   *mc_photons_mass;
   vector<float>   *mc_photons_numOfDaughters;
   UInt_t          Nmc_taus;
   vector<float>   *mc_taus_id;
   vector<float>   *mc_taus_pt;
   vector<float>   *mc_taus_px;
   vector<float>   *mc_taus_py;
   vector<float>   *mc_taus_pz;
   vector<float>   *mc_taus_eta;
   vector<float>   *mc_taus_phi;
   vector<float>   *mc_taus_theta;
   vector<float>   *mc_taus_status;
   vector<float>   *mc_taus_energy;
   vector<float>   *mc_taus_charge;
   vector<float>   *mc_taus_mother_id;
   vector<float>   *mc_taus_mother_pt;
   vector<float>   *mc_taus_grandmother_id;
   vector<float>   *mc_taus_ggrandmother_id;
   vector<float>   *mc_taus_vertex_x;
   vector<float>   *mc_taus_vertex_y;
   vector<float>   *mc_taus_vertex_z;
   vector<float>   *mc_taus_mass;
   vector<float>   *mc_taus_numOfDaughters;
   UInt_t          Nmets;
   vector<float>   *mets_et;
   vector<float>   *mets_phi;
   vector<float>   *mets_ex;
   vector<float>   *mets_ey;
   vector<float>   *mets_gen_et;
   vector<float>   *mets_gen_phi;
   vector<float>   *mets_sign;
   vector<float>   *mets_sumEt;
   vector<float>   *mets_unCPhi;
   vector<float>   *mets_unCPt;
   UInt_t          Nmus;
   vector<float>   *mus_energy;
   vector<float>   *mus_et;
   vector<float>   *mus_eta;
   vector<float>   *mus_phi;
   vector<float>   *mus_pt;
   vector<float>   *mus_px;
   vector<float>   *mus_py;
   vector<float>   *mus_pz;
   vector<float>   *mus_status;
   vector<float>   *mus_theta;
   vector<float>   *mus_gen_id;
   vector<float>   *mus_gen_phi;
   vector<float>   *mus_gen_pt;
   vector<float>   *mus_gen_pz;
   vector<float>   *mus_gen_px;
   vector<float>   *mus_gen_py;
   vector<float>   *mus_gen_eta;
   vector<float>   *mus_gen_theta;
   vector<float>   *mus_gen_et;
   vector<float>   *mus_gen_mother_id;
   vector<float>   *mus_gen_mother_phi;
   vector<float>   *mus_gen_mother_pt;
   vector<float>   *mus_gen_mother_pz;
   vector<float>   *mus_gen_mother_px;
   vector<float>   *mus_gen_mother_py;
   vector<float>   *mus_gen_mother_eta;
   vector<float>   *mus_gen_mother_theta;
   vector<float>   *mus_gen_mother_et;
   vector<float>   *mus_tkHits;
   vector<float>   *mus_cIso;
   vector<float>   *mus_tIso;
   vector<float>   *mus_ecalIso;
   vector<float>   *mus_hcalIso;
   vector<float>   *mus_ecalvetoDep;
   vector<float>   *mus_hcalvetoDep;
   vector<float>   *mus_calEnergyEm;
   vector<float>   *mus_calEnergyHad;
   vector<float>   *mus_calEnergyHo;
   vector<float>   *mus_calEnergyEmS9;
   vector<float>   *mus_calEnergyHadS9;
   vector<float>   *mus_calEnergyHoS9;
   vector<float>   *mus_iso03_emVetoEt;
   vector<float>   *mus_iso03_hadVetoEt;
   vector<float>   *mus_iso03_sumPt;
   vector<float>   *mus_iso03_emEt;
   vector<float>   *mus_iso03_hadEt;
   vector<float>   *mus_iso03_hoEt;
   vector<float>   *mus_iso03_nTracks;
   vector<float>   *mus_iso05_sumPt;
   vector<float>   *mus_iso05_emEt;
   vector<float>   *mus_iso05_hadEt;
   vector<float>   *mus_iso05_hoEt;
   vector<float>   *mus_iso05_nTracks;
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
   vector<float>   *mus_cm_chg;
   vector<float>   *mus_cm_pt;
   vector<float>   *mus_cm_px;
   vector<float>   *mus_cm_py;
   vector<float>   *mus_cm_pz;
   vector<float>   *mus_cm_eta;
   vector<float>   *mus_cm_phi;
   vector<float>   *mus_cm_theta;
   vector<float>   *mus_cm_d0dum;
   vector<float>   *mus_cm_dz;
   vector<float>   *mus_cm_vx;
   vector<float>   *mus_cm_vy;
   vector<float>   *mus_cm_vz;
   vector<float>   *mus_cm_numvalhits;
   vector<float>   *mus_cm_numlosthits;
   vector<float>   *mus_cm_numvalMuonhits;
   vector<float>   *mus_cm_d0dumErr;
   vector<float>   *mus_cm_dzErr;
   vector<float>   *mus_cm_ptErr;
   vector<float>   *mus_cm_etaErr;
   vector<float>   *mus_cm_phiErr;
   vector<float>   *mus_tk_id;
   vector<float>   *mus_tk_chi2;
   vector<float>   *mus_tk_ndof;
   vector<float>   *mus_tk_chg;
   vector<float>   *mus_tk_pt;
   vector<float>   *mus_tk_px;
   vector<float>   *mus_tk_py;
   vector<float>   *mus_tk_pz;
   vector<float>   *mus_tk_eta;
   vector<float>   *mus_tk_phi;
   vector<float>   *mus_tk_theta;
   vector<float>   *mus_tk_d0dum;
   vector<float>   *mus_tk_dz;
   vector<float>   *mus_tk_vx;
   vector<float>   *mus_tk_vy;
   vector<float>   *mus_tk_vz;
   vector<float>   *mus_tk_numvalhits;
   vector<float>   *mus_tk_numlosthits;
   vector<float>   *mus_tk_d0dumErr;
   vector<float>   *mus_tk_dzErr;
   vector<float>   *mus_tk_ptErr;
   vector<float>   *mus_tk_etaErr;
   vector<float>   *mus_tk_phiErr;
   vector<float>   *mus_tk_numvalPixelhits;
   vector<float>   *mus_tk_numpixelWthMeasr;
   vector<float>   *mus_stamu_chi2;
   vector<float>   *mus_stamu_ndof;
   vector<float>   *mus_stamu_chg;
   vector<float>   *mus_stamu_pt;
   vector<float>   *mus_stamu_px;
   vector<float>   *mus_stamu_py;
   vector<float>   *mus_stamu_pz;
   vector<float>   *mus_stamu_eta;
   vector<float>   *mus_stamu_phi;
   vector<float>   *mus_stamu_theta;
   vector<float>   *mus_stamu_d0dum;
   vector<float>   *mus_stamu_dz;
   vector<float>   *mus_stamu_vx;
   vector<float>   *mus_stamu_vy;
   vector<float>   *mus_stamu_vz;
   vector<float>   *mus_stamu_numvalhits;
   vector<float>   *mus_stamu_numlosthits;
   vector<float>   *mus_stamu_d0dumErr;
   vector<float>   *mus_stamu_dzErr;
   vector<float>   *mus_stamu_ptErr;
   vector<float>   *mus_stamu_etaErr;
   vector<float>   *mus_stamu_phiErr;
   vector<float>   *mus_num_matches;
   vector<float>   *mus_isPFMuon;
   vector<float>   *mus_isTrackerMuon;
   vector<float>   *mus_isStandAloneMuon;
   vector<float>   *mus_isCaloMuon;
   vector<float>   *mus_isGlobalMuon;
   vector<float>   *mus_isElectron;
   vector<float>   *mus_isConvertedPhoton;
   vector<float>   *mus_isPhoton;
   vector<float>   *mus_id_All;
   vector<float>   *mus_id_AllGlobalMuons;
   vector<float>   *mus_id_AllStandAloneMuons;
   vector<float>   *mus_id_AllTrackerMuons;
   vector<float>   *mus_id_TrackerMuonArbitrated;
   vector<float>   *mus_id_AllArbitrated;
   vector<float>   *mus_id_GlobalMuonPromptTight;
   vector<float>   *mus_id_TMLastStationLoose;
   vector<float>   *mus_id_TMLastStationTight;
   vector<float>   *mus_id_TM2DCompatibilityLoose;
   vector<float>   *mus_id_TM2DCompatibilityTight;
   vector<float>   *mus_id_TMOneStationLoose;
   vector<float>   *mus_id_TMOneStationTight;
   vector<float>   *mus_id_TMLastStationOptimizedLowPtLoose;
   vector<float>   *mus_id_TMLastStationOptimizedLowPtTight;
   vector<float>   *mus_tk_LayersWithMeasurement;
   vector<float>   *mus_tk_PixelLayersWithMeasurement;
   vector<float>   *mus_tk_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *mus_tk_LayersWithoutMeasurement;
   vector<float>   *mus_tk_ExpectedHitsInner;
   vector<float>   *mus_tk_ExpectedHitsOuter;
   vector<float>   *mus_cm_LayersWithMeasurement;
   vector<float>   *mus_cm_PixelLayersWithMeasurement;
   vector<float>   *mus_cm_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *mus_cm_LayersWithoutMeasurement;
   vector<float>   *mus_cm_ExpectedHitsInner;
   vector<float>   *mus_cm_ExpectedHitsOuter;
   vector<float>   *mus_picky_LayersWithMeasurement;
   vector<float>   *mus_picky_PixelLayersWithMeasurement;
   vector<float>   *mus_picky_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *mus_picky_LayersWithoutMeasurement;
   vector<float>   *mus_picky_ExpectedHitsInner;
   vector<float>   *mus_picky_ExpectedHitsOuter;
   vector<float>   *mus_tpfms_LayersWithMeasurement;
   vector<float>   *mus_tpfms_PixelLayersWithMeasurement;
   vector<float>   *mus_tpfms_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *mus_tpfms_LayersWithoutMeasurement;
   vector<float>   *mus_tpfms_ExpectedHitsInner;
   vector<float>   *mus_tpfms_ExpectedHitsOuter;
   vector<float>   *mus_picky_id;
   vector<float>   *mus_picky_chi2;
   vector<float>   *mus_picky_ndof;
   vector<float>   *mus_picky_chg;
   vector<float>   *mus_picky_pt;
   vector<float>   *mus_picky_px;
   vector<float>   *mus_picky_py;
   vector<float>   *mus_picky_pz;
   vector<float>   *mus_picky_eta;
   vector<float>   *mus_picky_phi;
   vector<float>   *mus_picky_theta;
   vector<float>   *mus_picky_d0dum;
   vector<float>   *mus_picky_dz;
   vector<float>   *mus_picky_vx;
   vector<float>   *mus_picky_vy;
   vector<float>   *mus_picky_vz;
   vector<float>   *mus_picky_numvalhits;
   vector<float>   *mus_picky_numlosthits;
   vector<float>   *mus_picky_d0dumErr;
   vector<float>   *mus_picky_dzErr;
   vector<float>   *mus_picky_ptErr;
   vector<float>   *mus_picky_etaErr;
   vector<float>   *mus_picky_phiErr;
   vector<float>   *mus_picky_numvalPixelhits;
   vector<float>   *mus_tpfms_id;
   vector<float>   *mus_tpfms_chi2;
   vector<float>   *mus_tpfms_ndof;
   vector<float>   *mus_tpfms_chg;
   vector<float>   *mus_tpfms_pt;
   vector<float>   *mus_tpfms_px;
   vector<float>   *mus_tpfms_py;
   vector<float>   *mus_tpfms_pz;
   vector<float>   *mus_tpfms_eta;
   vector<float>   *mus_tpfms_phi;
   vector<float>   *mus_tpfms_theta;
   vector<float>   *mus_tpfms_d0dum;
   vector<float>   *mus_tpfms_dz;
   vector<float>   *mus_tpfms_vx;
   vector<float>   *mus_tpfms_vy;
   vector<float>   *mus_tpfms_vz;
   vector<float>   *mus_tpfms_numvalhits;
   vector<float>   *mus_tpfms_numlosthits;
   vector<float>   *mus_tpfms_d0dumErr;
   vector<float>   *mus_tpfms_dzErr;
   vector<float>   *mus_tpfms_ptErr;
   vector<float>   *mus_tpfms_etaErr;
   vector<float>   *mus_tpfms_phiErr;
   vector<float>   *mus_tpfms_numvalPixelhits;
   vector<float>   *mus_dB;
   vector<float>   *mus_numberOfMatchedStations;
   UInt_t          Npfcand;
   vector<float>   *pfcand_pdgId;
   vector<float>   *pfcand_pt;
   vector<float>   *pfcand_pz;
   vector<float>   *pfcand_px;
   vector<float>   *pfcand_py;
   vector<float>   *pfcand_eta;
   vector<float>   *pfcand_phi;
   vector<float>   *pfcand_theta;
   vector<float>   *pfcand_energy;
   vector<float>   *pfcand_charge;
   UInt_t          Nphotons;
   vector<float>   *photons_energy;
   vector<float>   *photons_et;
   vector<float>   *photons_eta;
   vector<float>   *photons_phi;
   vector<float>   *photons_pt;
   vector<float>   *photons_px;
   vector<float>   *photons_py;
   vector<float>   *photons_pz;
   vector<float>   *photons_status;
   vector<float>   *photons_theta;
   vector<float>   *photons_hadOverEM;
   vector<float>   *photons_hadTowOverEM;
   vector<float>   *photons_scEnergy;
   vector<float>   *photons_scRawEnergy;
   vector<float>   *photons_scEta;
   vector<float>   *photons_scPhi;
   vector<float>   *photons_scEtaWidth;
   vector<float>   *photons_scPhiWidth;
   vector<float>   *photons_tIso;
   vector<float>   *photons_ecalIso;
   vector<float>   *photons_hcalIso;
   vector<float>   *photons_isoEcalRecHitDR04;
   vector<float>   *photons_isoHcalRecHitDR04;
   vector<float>   *photons_isoSolidTrkConeDR04;
   vector<float>   *photons_isoHollowTrkConeDR04;
   vector<float>   *photons_nTrkSolidConeDR04;
   vector<float>   *photons_nTrkHollowConeDR04;
   vector<float>   *photons_isoEcalRecHitDR03;
   vector<float>   *photons_isoHcalRecHitDR03;
   vector<float>   *photons_isoSolidTrkConeDR03;
   vector<float>   *photons_isoHollowTrkConeDR03;
   vector<float>   *photons_nTrkSolidConeDR03;
   vector<float>   *photons_nTrkHollowConeDR03;
   vector<float>   *photons_isAlsoElectron;
   vector<float>   *photons_hasPixelSeed;
   vector<float>   *photons_isConverted;
   vector<float>   *photons_isEBGap;
   vector<float>   *photons_isEEGap;
   vector<float>   *photons_isEBEEGap;
   vector<float>   *photons_isEBPho;
   vector<float>   *photons_isEEPho;
   vector<float>   *photons_isLoosePhoton;
   vector<float>   *photons_isTightPhoton;
   vector<float>   *photons_maxEnergyXtal;
   vector<float>   *photons_e1x5;
   vector<float>   *photons_e2x5;
   vector<float>   *photons_e3x3;
   vector<float>   *photons_e5x5;
   vector<float>   *photons_sigmaEtaEta;
   vector<float>   *photons_sigmaIetaIeta;
   vector<float>   *photons_r9;
   vector<float>   *photons_gen_et;
   vector<float>   *photons_gen_eta;
   vector<float>   *photons_gen_phi;
   vector<float>   *photons_gen_id;
   UInt_t          Npv;
   vector<float>   *pv_x;
   vector<float>   *pv_y;
   vector<float>   *pv_z;
   vector<float>   *pv_xErr;
   vector<float>   *pv_yErr;
   vector<float>   *pv_zErr;
   vector<float>   *pv_chi2;
   vector<float>   *pv_ndof;
   vector<float>   *pv_isFake;
   vector<float>   *pv_isValid;
   vector<float>   *pv_tracksSize;
   UInt_t          Ntaus;
   vector<float>   *taus_status;
   vector<float>   *taus_phi;
   vector<float>   *taus_pt;
   vector<float>   *taus_pz;
   vector<float>   *taus_px;
   vector<float>   *taus_py;
   vector<float>   *taus_eta;
   vector<float>   *taus_theta;
   vector<float>   *taus_et;
   vector<float>   *taus_energy;
   vector<float>   *taus_charge;
   vector<float>   *taus_emf;
   vector<float>   *taus_hcalTotOverPLead;
   vector<float>   *taus_hcalMaxOverPLead;
   vector<float>   *taus_hcal3x3OverPLead;
   vector<float>   *taus_ecalStripSumEOverPLead;
   vector<float>   *taus_elecPreIdOutput;
   vector<float>   *taus_elecPreIdDecision;
   vector<float>   *taus_leadPFChargedHadrCand_pt;
   vector<float>   *taus_leadPFChargedHadrCand_charge;
   vector<float>   *taus_leadPFChargedHadrCand_eta;
   vector<float>   *taus_leadPFChargedHadrCand_ECAL_eta;
   vector<float>   *taus_leadPFChargedHadrCand_phi;
   vector<float>   *taus_isoPFGammaCandsEtSum;
   vector<float>   *taus_isoPFChargedHadrCandsPtSum;
   vector<float>   *taus_leadingTrackFinding;
   vector<float>   *taus_leadingTrackPtCut;
   vector<float>   *taus_trackIsolation;
   vector<float>   *taus_ecalIsolation;
   vector<float>   *taus_byIsolation;
   vector<float>   *taus_againstElectron;
   vector<float>   *taus_againstMuon;
   vector<float>   *taus_taNC_quarter;
   vector<float>   *taus_taNC_one;
   vector<float>   *taus_taNC_half;
   vector<float>   *taus_taNC_tenth;
   vector<float>   *taus_taNC;
   vector<float>   *taus_byIsoUsingLeadingPi;
   vector<float>   *taus_tkIsoUsingLeadingPi;
   vector<float>   *taus_ecalIsoUsingLeadingPi;
   vector<float>   *taus_againstElectronLoose;
   vector<float>   *taus_againstElectronMedium;
   vector<float>   *taus_againstElectronTight;
   vector<float>   *taus_againstElectronMVA;
   vector<float>   *taus_againstMuonLoose;
   vector<float>   *taus_againstMuonMedium;
   vector<float>   *taus_againstMuonTight;
   vector<float>   *taus_decayModeFinding;
   vector<float>   *taus_byVLooseIsolation;
   vector<float>   *taus_byLooseIsolation;
   vector<float>   *taus_byMediumIsolation;
   vector<float>   *taus_byTightIsolation;
   vector<float>   *taus_byVLooseIsolationDeltaBetaCorr;
   vector<float>   *taus_byLooseIsolationDeltaBetaCorr;
   vector<float>   *taus_byMediumIsolationDeltaBetaCorr;
   vector<float>   *taus_byTightIsolationDeltaBetaCorr;
   vector<float>   *taus_signalPFChargedHadrCandsSize;
   vector<float>   *taus_muDecision;
   vector<float>   *taus_Nprongs;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          experimentType;
   UInt_t          bunchCrossing;
   UInt_t          orbitNumber;
   Float_t         weight;
   string          *model_params;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT10_px;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT10_py;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT10_pz;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT10_energy;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT10_phi;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT10_eta;
   vector<vector<int> > *fastjets_AK4_R1p2_R0p5pT10_index;
   vector<int>     *fastjets_AK4_R1p2_R0p5pT10_nconstituents;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT15_px;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT15_py;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT15_pz;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT15_energy;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT15_phi;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT15_eta;
   vector<vector<int> > *fastjets_AK4_R1p2_R0p5pT15_index;
   vector<int>     *fastjets_AK4_R1p2_R0p5pT15_nconstituents;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT20_px;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT20_py;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT20_pz;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT20_energy;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT20_phi;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT20_eta;
   vector<vector<int> > *fastjets_AK4_R1p2_R0p5pT20_index;
   vector<int>     *fastjets_AK4_R1p2_R0p5pT20_nconstituents;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT25_px;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT25_py;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT25_pz;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT25_energy;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT25_phi;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT25_eta;
   vector<vector<int> > *fastjets_AK4_R1p2_R0p5pT25_index;
   vector<int>     *fastjets_AK4_R1p2_R0p5pT25_nconstituents;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT30_px;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT30_py;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT30_pz;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT30_energy;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT30_phi;
   vector<float>   *fastjets_AK4_R1p2_R0p5pT30_eta;
   vector<vector<int> > *fastjets_AK4_R1p2_R0p5pT30_index;
   vector<int>     *fastjets_AK4_R1p2_R0p5pT30_nconstituents;

   // List of branches
   TBranch        *b_NbeamSpot;   //!
   TBranch        *b_beamSpot_x;   //!
   TBranch        *b_beamSpot_y;   //!
   TBranch        *b_beamSpot_z;   //!
   TBranch        *b_beamSpot_x0Error;   //!
   TBranch        *b_beamSpot_y0Error;   //!
   TBranch        *b_beamSpot_z0Error;   //!
   TBranch        *b_beamSpot_sigmaZ;   //!
   TBranch        *b_beamSpot_sigmaZ0Error;   //!
   TBranch        *b_beamSpot_dxdz;   //!
   TBranch        *b_beamSpot_dxdzError;   //!
   TBranch        *b_beamSpot_dydz;   //!
   TBranch        *b_beamSpot_dydzError;   //!
   TBranch        *b_beamSpot_beamWidthX;   //!
   TBranch        *b_beamSpot_beamWidthY;   //!
   TBranch        *b_beamSpot_beamWidthXError;   //!
   TBranch        *b_beamSpot_beamWidthYError;   //!
   TBranch        *b_Nels;   //!
   TBranch        *b_els_energy;   //!
   TBranch        *b_els_et;   //!
   TBranch        *b_els_eta;   //!
   TBranch        *b_els_phi;   //!
   TBranch        *b_els_pt;   //!
   TBranch        *b_els_px;   //!
   TBranch        *b_els_py;   //!
   TBranch        *b_els_pz;   //!
   TBranch        *b_els_status;   //!
   TBranch        *b_els_theta;   //!
   TBranch        *b_els_pfIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_els_pfIsolationR03_sumNeutralHadronEt;   //!
   TBranch        *b_els_pfIsolationR03_sumPhotonEt;   //!
   TBranch        *b_els_pfIsolationR03_sumPUPt;   //!
   TBranch        *b_els_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_els_gen_id;   //!
   TBranch        *b_els_gen_phi;   //!
   TBranch        *b_els_gen_pt;   //!
   TBranch        *b_els_gen_pz;   //!
   TBranch        *b_els_gen_px;   //!
   TBranch        *b_els_gen_py;   //!
   TBranch        *b_els_gen_eta;   //!
   TBranch        *b_els_gen_theta;   //!
   TBranch        *b_els_gen_et;   //!
   TBranch        *b_els_gen_mother_id;   //!
   TBranch        *b_els_gen_mother_phi;   //!
   TBranch        *b_els_gen_mother_pt;   //!
   TBranch        *b_els_gen_mother_pz;   //!
   TBranch        *b_els_gen_mother_px;   //!
   TBranch        *b_els_gen_mother_py;   //!
   TBranch        *b_els_gen_mother_eta;   //!
   TBranch        *b_els_gen_mother_theta;   //!
   TBranch        *b_els_gen_mother_et;   //!
   TBranch        *b_els_tightId;   //!
   TBranch        *b_els_looseId;   //!
   TBranch        *b_els_robustTightId;   //!
   TBranch        *b_els_robustLooseId;   //!
   TBranch        *b_els_robustHighEnergyId;   //!
   TBranch        *b_els_simpleEleId95relIso;   //!
   TBranch        *b_els_simpleEleId90relIso;   //!
   TBranch        *b_els_simpleEleId85relIso;   //!
   TBranch        *b_els_simpleEleId80relIso;   //!
   TBranch        *b_els_simpleEleId70relIso;   //!
   TBranch        *b_els_simpleEleId60relIso;   //!
   TBranch        *b_els_simpleEleId95cIso;   //!
   TBranch        *b_els_simpleEleId90cIso;   //!
   TBranch        *b_els_simpleEleId85cIso;   //!
   TBranch        *b_els_simpleEleId80cIso;   //!
   TBranch        *b_els_simpleEleId70cIso;   //!
   TBranch        *b_els_simpleEleId60cIso;   //!
   TBranch        *b_els_cIso;   //!
   TBranch        *b_els_tIso;   //!
   TBranch        *b_els_ecalIso;   //!
   TBranch        *b_els_hcalIso;   //!
   TBranch        *b_els_chi2;   //!
   TBranch        *b_els_charge;   //!
   TBranch        *b_els_caloEnergy;   //!
   TBranch        *b_els_hadOverEm;   //!
   TBranch        *b_els_hcalOverEcalBc;   //!
   TBranch        *b_els_eOverPIn;   //!
   TBranch        *b_els_eSeedOverPOut;   //!
   TBranch        *b_els_sigmaEtaEta;   //!
   TBranch        *b_els_sigmaIEtaIEta;   //!
   TBranch        *b_els_scEnergy;   //!
   TBranch        *b_els_scRawEnergy;   //!
   TBranch        *b_els_scSeedEnergy;   //!
   TBranch        *b_els_scEta;   //!
   TBranch        *b_els_scPhi;   //!
   TBranch        *b_els_scEtaWidth;   //!
   TBranch        *b_els_scPhiWidth;   //!
   TBranch        *b_els_scE1x5;   //!
   TBranch        *b_els_scE2x5Max;   //!
   TBranch        *b_els_scE5x5;   //!
   TBranch        *b_els_isEB;   //!
   TBranch        *b_els_isEE;   //!
   TBranch        *b_els_dEtaIn;   //!
   TBranch        *b_els_dPhiIn;   //!
   TBranch        *b_els_dEtaOut;   //!
   TBranch        *b_els_dPhiOut;   //!
   TBranch        *b_els_numvalhits;   //!
   TBranch        *b_els_numlosthits;   //!
   TBranch        *b_els_basicClustersSize;   //!
   TBranch        *b_els_tk_pz;   //!
   TBranch        *b_els_tk_pt;   //!
   TBranch        *b_els_tk_phi;   //!
   TBranch        *b_els_tk_eta;   //!
   TBranch        *b_els_d0dum;   //!
   TBranch        *b_els_dz;   //!
   TBranch        *b_els_vx;   //!
   TBranch        *b_els_vy;   //!
   TBranch        *b_els_vz;   //!
   TBranch        *b_els_ndof;   //!
   TBranch        *b_els_ptError;   //!
   TBranch        *b_els_d0dumError;   //!
   TBranch        *b_els_dzError;   //!
   TBranch        *b_els_etaError;   //!
   TBranch        *b_els_phiError;   //!
   TBranch        *b_els_tk_charge;   //!
   TBranch        *b_els_core_ecalDrivenSeed;   //!
   TBranch        *b_els_n_inner_layer;   //!
   TBranch        *b_els_n_outer_layer;   //!
   TBranch        *b_els_ctf_tk_id;   //!
   TBranch        *b_els_ctf_tk_charge;   //!
   TBranch        *b_els_ctf_tk_eta;   //!
   TBranch        *b_els_ctf_tk_phi;   //!
   TBranch        *b_els_fbrem;   //!
   TBranch        *b_els_shFracInnerHits;   //!
   TBranch        *b_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_els_dr03TkSumPt;   //!
   TBranch        *b_els_dr04EcalRecHitSumEt;   //!
   TBranch        *b_els_dr04HcalTowerSumEt;   //!
   TBranch        *b_els_dr04HcalDepth1TowerSumEt;   //!
   TBranch        *b_els_dr04HcalDepth2TowerSumEt;   //!
   TBranch        *b_els_dr04TkSumPt;   //!
   TBranch        *b_els_cpx;   //!
   TBranch        *b_els_cpy;   //!
   TBranch        *b_els_cpz;   //!
   TBranch        *b_els_vpx;   //!
   TBranch        *b_els_vpy;   //!
   TBranch        *b_els_vpz;   //!
   TBranch        *b_els_cx;   //!
   TBranch        *b_els_cy;   //!
   TBranch        *b_els_cz;   //!
   TBranch        *b_els_PATpassConversionVeto;   //!
   TBranch        *b_Njets_AK4;   //!
   TBranch        *b_jets_AK4_status;   //!
   TBranch        *b_jets_AK4_phi;   //!
   TBranch        *b_jets_AK4_pt;   //!
   TBranch        *b_jets_AK4_pz;   //!
   TBranch        *b_jets_AK4_px;   //!
   TBranch        *b_jets_AK4_py;   //!
   TBranch        *b_jets_AK4_eta;   //!
   TBranch        *b_jets_AK4_theta;   //!
   TBranch        *b_jets_AK4_et;   //!
   TBranch        *b_jets_AK4_energy;   //!
   TBranch        *b_jets_AK4_parton_Id;   //!
   TBranch        *b_jets_AK4_parton_motherId;   //!
   TBranch        *b_jets_AK4_parton_pt;   //!
   TBranch        *b_jets_AK4_parton_phi;   //!
   TBranch        *b_jets_AK4_parton_eta;   //!
   TBranch        *b_jets_AK4_parton_Energy;   //!
   TBranch        *b_jets_AK4_parton_mass;   //!
   TBranch        *b_jets_AK4_gen_et;   //!
   TBranch        *b_jets_AK4_gen_pt;   //!
   TBranch        *b_jets_AK4_gen_eta;   //!
   TBranch        *b_jets_AK4_gen_phi;   //!
   TBranch        *b_jets_AK4_gen_mass;   //!
   TBranch        *b_jets_AK4_gen_Energy;   //!
   TBranch        *b_jets_AK4_gen_Id;   //!
   TBranch        *b_jets_AK4_gen_motherID;   //!
   TBranch        *b_jets_AK4_gen_threeCharge;   //!
   TBranch        *b_jets_AK4_partonFlavour;   //!
   TBranch        *b_jets_AK4_btag_TC_highPur;   //!
   TBranch        *b_jets_AK4_btag_TC_highEff;   //!
   TBranch        *b_jets_AK4_btag_jetProb;   //!
   TBranch        *b_jets_AK4_btag_jetBProb;   //!
   TBranch        *b_jets_AK4_btag_softEle;   //!
   TBranch        *b_jets_AK4_btag_softMuon;   //!
   TBranch        *b_jets_AK4_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK4_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK4_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK4_jetCharge;   //!
   TBranch        *b_jets_AK4_chgEmE;   //!
   TBranch        *b_jets_AK4_chgHadE;   //!
   TBranch        *b_jets_AK4_photonEnergy;   //!
   TBranch        *b_jets_AK4_chgMuE;   //!
   TBranch        *b_jets_AK4_chg_Mult;   //!
   TBranch        *b_jets_AK4_neutralEmE;   //!
   TBranch        *b_jets_AK4_neutralHadE;   //!
   TBranch        *b_jets_AK4_neutral_Mult;   //!
   TBranch        *b_jets_AK4_mu_Mult;   //!
   TBranch        *b_jets_AK4_emf;   //!
   TBranch        *b_jets_AK4_ehf;   //!
   TBranch        *b_jets_AK4_n60;   //!
   TBranch        *b_jets_AK4_n90;   //!
   TBranch        *b_jets_AK4_etaetaMoment;   //!
   TBranch        *b_jets_AK4_etaphiMoment;   //!
   TBranch        *b_jets_AK4_phiphiMoment;   //!
   TBranch        *b_jets_AK4_n90Hits;   //!
   TBranch        *b_jets_AK4_fHPD;   //!
   TBranch        *b_jets_AK4_fRBX;   //!
   TBranch        *b_jets_AK4_hitsInN90;   //!
   TBranch        *b_jets_AK4_nECALTowers;   //!
   TBranch        *b_jets_AK4_nHCALTowers;   //!
   TBranch        *b_jets_AK4_fSubDetector1;   //!
   TBranch        *b_jets_AK4_fSubDetector2;   //!
   TBranch        *b_jets_AK4_fSubDetector3;   //!
   TBranch        *b_jets_AK4_fSubDetector4;   //!
   TBranch        *b_jets_AK4_area;   //!
   TBranch        *b_jets_AK4_corrFactorRaw;   //!
   TBranch        *b_jets_AK4_rawPt;   //!
   TBranch        *b_jets_AK4_mass;   //!
   TBranch        *b_Nmc_doc;   //!
   TBranch        *b_mc_doc_id;   //!
   TBranch        *b_mc_doc_pt;   //!
   TBranch        *b_mc_doc_px;   //!
   TBranch        *b_mc_doc_py;   //!
   TBranch        *b_mc_doc_pz;   //!
   TBranch        *b_mc_doc_eta;   //!
   TBranch        *b_mc_doc_phi;   //!
   TBranch        *b_mc_doc_theta;   //!
   TBranch        *b_mc_doc_energy;   //!
   TBranch        *b_mc_doc_status;   //!
   TBranch        *b_mc_doc_charge;   //!
   TBranch        *b_mc_doc_mother_id;   //!
   TBranch        *b_mc_doc_grandmother_id;   //!
   TBranch        *b_mc_doc_ggrandmother_id;   //!
   TBranch        *b_mc_doc_mother_pt;   //!
   TBranch        *b_mc_doc_vertex_x;   //!
   TBranch        *b_mc_doc_vertex_y;   //!
   TBranch        *b_mc_doc_vertex_z;   //!
   TBranch        *b_mc_doc_mass;   //!
   TBranch        *b_mc_doc_numOfDaughters;   //!
   TBranch        *b_mc_doc_numOfMothers;   //!
   TBranch        *b_Nmc_electrons;   //!
   TBranch        *b_mc_electrons_id;   //!
   TBranch        *b_mc_electrons_pt;   //!
   TBranch        *b_mc_electrons_px;   //!
   TBranch        *b_mc_electrons_py;   //!
   TBranch        *b_mc_electrons_pz;   //!
   TBranch        *b_mc_electrons_eta;   //!
   TBranch        *b_mc_electrons_phi;   //!
   TBranch        *b_mc_electrons_theta;   //!
   TBranch        *b_mc_electrons_status;   //!
   TBranch        *b_mc_electrons_energy;   //!
   TBranch        *b_mc_electrons_charge;   //!
   TBranch        *b_mc_electrons_mother_id;   //!
   TBranch        *b_mc_electrons_mother_pt;   //!
   TBranch        *b_mc_electrons_grandmother_id;   //!
   TBranch        *b_mc_electrons_ggrandmother_id;   //!
   TBranch        *b_mc_electrons_vertex_x;   //!
   TBranch        *b_mc_electrons_vertex_y;   //!
   TBranch        *b_mc_electrons_vertex_z;   //!
   TBranch        *b_mc_electrons_mass;   //!
   TBranch        *b_mc_electrons_numOfDaughters;   //!
   TBranch        *b_Nmc_final;   //!
   TBranch        *b_mc_final_id;   //!
   TBranch        *b_mc_final_pt;   //!
   TBranch        *b_mc_final_px;   //!
   TBranch        *b_mc_final_py;   //!
   TBranch        *b_mc_final_pz;   //!
   TBranch        *b_mc_final_eta;   //!
   TBranch        *b_mc_final_phi;   //!
   TBranch        *b_mc_final_theta;   //!
   TBranch        *b_mc_final_energy;   //!
   TBranch        *b_mc_final_status;   //!
   TBranch        *b_mc_final_charge;   //!
   TBranch        *b_mc_final_mother_id;   //!
   TBranch        *b_mc_final_grandmother_id;   //!
   TBranch        *b_mc_final_ggrandmother_id;   //!
   TBranch        *b_mc_final_mother_pt;   //!
   TBranch        *b_mc_final_vertex_x;   //!
   TBranch        *b_mc_final_vertex_y;   //!
   TBranch        *b_mc_final_vertex_z;   //!
   TBranch        *b_mc_final_mass;   //!
   TBranch        *b_mc_final_numOfDaughters;   //!
   TBranch        *b_mc_final_numOfMothers;   //!
   TBranch        *b_Nmc_jets;   //!
   TBranch        *b_mc_jets_phi;   //!
   TBranch        *b_mc_jets_pt;   //!
   TBranch        *b_mc_jets_pz;   //!
   TBranch        *b_mc_jets_px;   //!
   TBranch        *b_mc_jets_py;   //!
   TBranch        *b_mc_jets_eta;   //!
   TBranch        *b_mc_jets_theta;   //!
   TBranch        *b_mc_jets_et;   //!
   TBranch        *b_mc_jets_energy;   //!
   TBranch        *b_mc_jets_emEnergy;   //!
   TBranch        *b_mc_jets_hadEnergy;   //!
   TBranch        *b_mc_jets_invisibleEnergy;   //!
   TBranch        *b_mc_jets_auxiliaryEnergy;   //!
   TBranch        *b_mc_jets_etaetaMoment;   //!
   TBranch        *b_mc_jets_etaphiMoment;   //!
   TBranch        *b_mc_jets_phiphiMoment;   //!
   TBranch        *b_mc_jets_mass;   //!
   TBranch        *b_Nmc_mus;   //!
   TBranch        *b_mc_mus_id;   //!
   TBranch        *b_mc_mus_pt;   //!
   TBranch        *b_mc_mus_px;   //!
   TBranch        *b_mc_mus_py;   //!
   TBranch        *b_mc_mus_pz;   //!
   TBranch        *b_mc_mus_eta;   //!
   TBranch        *b_mc_mus_phi;   //!
   TBranch        *b_mc_mus_theta;   //!
   TBranch        *b_mc_mus_status;   //!
   TBranch        *b_mc_mus_energy;   //!
   TBranch        *b_mc_mus_charge;   //!
   TBranch        *b_mc_mus_mother_id;   //!
   TBranch        *b_mc_mus_mother_pt;   //!
   TBranch        *b_mc_mus_grandmother_id;   //!
   TBranch        *b_mc_mus_ggrandmother_id;   //!
   TBranch        *b_mc_mus_vertex_x;   //!
   TBranch        *b_mc_mus_vertex_y;   //!
   TBranch        *b_mc_mus_vertex_z;   //!
   TBranch        *b_mc_mus_mass;   //!
   TBranch        *b_mc_mus_numOfDaughters;   //!
   TBranch        *b_Nmc_nues;   //!
   TBranch        *b_mc_nues_id;   //!
   TBranch        *b_mc_nues_pt;   //!
   TBranch        *b_mc_nues_px;   //!
   TBranch        *b_mc_nues_py;   //!
   TBranch        *b_mc_nues_pz;   //!
   TBranch        *b_mc_nues_eta;   //!
   TBranch        *b_mc_nues_phi;   //!
   TBranch        *b_mc_nues_theta;   //!
   TBranch        *b_mc_nues_status;   //!
   TBranch        *b_mc_nues_energy;   //!
   TBranch        *b_mc_nues_charge;   //!
   TBranch        *b_mc_nues_mother_id;   //!
   TBranch        *b_mc_nues_mother_pt;   //!
   TBranch        *b_mc_nues_grandmother_id;   //!
   TBranch        *b_mc_nues_ggrandmother_id;   //!
   TBranch        *b_mc_nues_vertex_x;   //!
   TBranch        *b_mc_nues_vertex_y;   //!
   TBranch        *b_mc_nues_vertex_z;   //!
   TBranch        *b_mc_nues_mass;   //!
   TBranch        *b_mc_nues_numOfDaughters;   //!
   TBranch        *b_Nmc_numus;   //!
   TBranch        *b_mc_numus_id;   //!
   TBranch        *b_mc_numus_pt;   //!
   TBranch        *b_mc_numus_px;   //!
   TBranch        *b_mc_numus_py;   //!
   TBranch        *b_mc_numus_pz;   //!
   TBranch        *b_mc_numus_eta;   //!
   TBranch        *b_mc_numus_phi;   //!
   TBranch        *b_mc_numus_theta;   //!
   TBranch        *b_mc_numus_status;   //!
   TBranch        *b_mc_numus_energy;   //!
   TBranch        *b_mc_numus_charge;   //!
   TBranch        *b_mc_numus_mother_id;   //!
   TBranch        *b_mc_numus_mother_pt;   //!
   TBranch        *b_mc_numus_grandmother_id;   //!
   TBranch        *b_mc_numus_ggrandmother_id;   //!
   TBranch        *b_mc_numus_vertex_x;   //!
   TBranch        *b_mc_numus_vertex_y;   //!
   TBranch        *b_mc_numus_vertex_z;   //!
   TBranch        *b_mc_numus_mass;   //!
   TBranch        *b_mc_numus_numOfDaughters;   //!
   TBranch        *b_Nmc_nutaus;   //!
   TBranch        *b_mc_nutaus_id;   //!
   TBranch        *b_mc_nutaus_pt;   //!
   TBranch        *b_mc_nutaus_px;   //!
   TBranch        *b_mc_nutaus_py;   //!
   TBranch        *b_mc_nutaus_pz;   //!
   TBranch        *b_mc_nutaus_eta;   //!
   TBranch        *b_mc_nutaus_phi;   //!
   TBranch        *b_mc_nutaus_theta;   //!
   TBranch        *b_mc_nutaus_status;   //!
   TBranch        *b_mc_nutaus_energy;   //!
   TBranch        *b_mc_nutaus_charge;   //!
   TBranch        *b_mc_nutaus_mother_id;   //!
   TBranch        *b_mc_nutaus_mother_pt;   //!
   TBranch        *b_mc_nutaus_grandmother_id;   //!
   TBranch        *b_mc_nutaus_ggrandmother_id;   //!
   TBranch        *b_mc_nutaus_vertex_x;   //!
   TBranch        *b_mc_nutaus_vertex_y;   //!
   TBranch        *b_mc_nutaus_vertex_z;   //!
   TBranch        *b_mc_nutaus_mass;   //!
   TBranch        *b_mc_nutaus_numOfDaughters;   //!
   TBranch        *b_Nmc_photons;   //!
   TBranch        *b_mc_photons_id;   //!
   TBranch        *b_mc_photons_pt;   //!
   TBranch        *b_mc_photons_px;   //!
   TBranch        *b_mc_photons_py;   //!
   TBranch        *b_mc_photons_pz;   //!
   TBranch        *b_mc_photons_eta;   //!
   TBranch        *b_mc_photons_phi;   //!
   TBranch        *b_mc_photons_theta;   //!
   TBranch        *b_mc_photons_status;   //!
   TBranch        *b_mc_photons_energy;   //!
   TBranch        *b_mc_photons_charge;   //!
   TBranch        *b_mc_photons_mother_id;   //!
   TBranch        *b_mc_photons_mother_pt;   //!
   TBranch        *b_mc_photons_grandmother_id;   //!
   TBranch        *b_mc_photons_ggrandmother_id;   //!
   TBranch        *b_mc_photons_vertex_x;   //!
   TBranch        *b_mc_photons_vertex_y;   //!
   TBranch        *b_mc_photons_vertex_z;   //!
   TBranch        *b_mc_photons_mass;   //!
   TBranch        *b_mc_photons_numOfDaughters;   //!
   TBranch        *b_Nmc_taus;   //!
   TBranch        *b_mc_taus_id;   //!
   TBranch        *b_mc_taus_pt;   //!
   TBranch        *b_mc_taus_px;   //!
   TBranch        *b_mc_taus_py;   //!
   TBranch        *b_mc_taus_pz;   //!
   TBranch        *b_mc_taus_eta;   //!
   TBranch        *b_mc_taus_phi;   //!
   TBranch        *b_mc_taus_theta;   //!
   TBranch        *b_mc_taus_status;   //!
   TBranch        *b_mc_taus_energy;   //!
   TBranch        *b_mc_taus_charge;   //!
   TBranch        *b_mc_taus_mother_id;   //!
   TBranch        *b_mc_taus_mother_pt;   //!
   TBranch        *b_mc_taus_grandmother_id;   //!
   TBranch        *b_mc_taus_ggrandmother_id;   //!
   TBranch        *b_mc_taus_vertex_x;   //!
   TBranch        *b_mc_taus_vertex_y;   //!
   TBranch        *b_mc_taus_vertex_z;   //!
   TBranch        *b_mc_taus_mass;   //!
   TBranch        *b_mc_taus_numOfDaughters;   //!
   TBranch        *b_Nmets;   //!
   TBranch        *b_mets_et;   //!
   TBranch        *b_mets_phi;   //!
   TBranch        *b_mets_ex;   //!
   TBranch        *b_mets_ey;   //!
   TBranch        *b_mets_gen_et;   //!
   TBranch        *b_mets_gen_phi;   //!
   TBranch        *b_mets_sign;   //!
   TBranch        *b_mets_sumEt;   //!
   TBranch        *b_mets_unCPhi;   //!
   TBranch        *b_mets_unCPt;   //!
   TBranch        *b_Nmus;   //!
   TBranch        *b_mus_energy;   //!
   TBranch        *b_mus_et;   //!
   TBranch        *b_mus_eta;   //!
   TBranch        *b_mus_phi;   //!
   TBranch        *b_mus_pt;   //!
   TBranch        *b_mus_px;   //!
   TBranch        *b_mus_py;   //!
   TBranch        *b_mus_pz;   //!
   TBranch        *b_mus_status;   //!
   TBranch        *b_mus_theta;   //!
   TBranch        *b_mus_gen_id;   //!
   TBranch        *b_mus_gen_phi;   //!
   TBranch        *b_mus_gen_pt;   //!
   TBranch        *b_mus_gen_pz;   //!
   TBranch        *b_mus_gen_px;   //!
   TBranch        *b_mus_gen_py;   //!
   TBranch        *b_mus_gen_eta;   //!
   TBranch        *b_mus_gen_theta;   //!
   TBranch        *b_mus_gen_et;   //!
   TBranch        *b_mus_gen_mother_id;   //!
   TBranch        *b_mus_gen_mother_phi;   //!
   TBranch        *b_mus_gen_mother_pt;   //!
   TBranch        *b_mus_gen_mother_pz;   //!
   TBranch        *b_mus_gen_mother_px;   //!
   TBranch        *b_mus_gen_mother_py;   //!
   TBranch        *b_mus_gen_mother_eta;   //!
   TBranch        *b_mus_gen_mother_theta;   //!
   TBranch        *b_mus_gen_mother_et;   //!
   TBranch        *b_mus_tkHits;   //!
   TBranch        *b_mus_cIso;   //!
   TBranch        *b_mus_tIso;   //!
   TBranch        *b_mus_ecalIso;   //!
   TBranch        *b_mus_hcalIso;   //!
   TBranch        *b_mus_ecalvetoDep;   //!
   TBranch        *b_mus_hcalvetoDep;   //!
   TBranch        *b_mus_calEnergyEm;   //!
   TBranch        *b_mus_calEnergyHad;   //!
   TBranch        *b_mus_calEnergyHo;   //!
   TBranch        *b_mus_calEnergyEmS9;   //!
   TBranch        *b_mus_calEnergyHadS9;   //!
   TBranch        *b_mus_calEnergyHoS9;   //!
   TBranch        *b_mus_iso03_emVetoEt;   //!
   TBranch        *b_mus_iso03_hadVetoEt;   //!
   TBranch        *b_mus_iso03_sumPt;   //!
   TBranch        *b_mus_iso03_emEt;   //!
   TBranch        *b_mus_iso03_hadEt;   //!
   TBranch        *b_mus_iso03_hoEt;   //!
   TBranch        *b_mus_iso03_nTracks;   //!
   TBranch        *b_mus_iso05_sumPt;   //!
   TBranch        *b_mus_iso05_emEt;   //!
   TBranch        *b_mus_iso05_hadEt;   //!
   TBranch        *b_mus_iso05_hoEt;   //!
   TBranch        *b_mus_iso05_nTracks;   //!
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
   TBranch        *b_mus_cm_chg;   //!
   TBranch        *b_mus_cm_pt;   //!
   TBranch        *b_mus_cm_px;   //!
   TBranch        *b_mus_cm_py;   //!
   TBranch        *b_mus_cm_pz;   //!
   TBranch        *b_mus_cm_eta;   //!
   TBranch        *b_mus_cm_phi;   //!
   TBranch        *b_mus_cm_theta;   //!
   TBranch        *b_mus_cm_d0dum;   //!
   TBranch        *b_mus_cm_dz;   //!
   TBranch        *b_mus_cm_vx;   //!
   TBranch        *b_mus_cm_vy;   //!
   TBranch        *b_mus_cm_vz;   //!
   TBranch        *b_mus_cm_numvalhits;   //!
   TBranch        *b_mus_cm_numlosthits;   //!
   TBranch        *b_mus_cm_numvalMuonhits;   //!
   TBranch        *b_mus_cm_d0dumErr;   //!
   TBranch        *b_mus_cm_dzErr;   //!
   TBranch        *b_mus_cm_ptErr;   //!
   TBranch        *b_mus_cm_etaErr;   //!
   TBranch        *b_mus_cm_phiErr;   //!
   TBranch        *b_mus_tk_id;   //!
   TBranch        *b_mus_tk_chi2;   //!
   TBranch        *b_mus_tk_ndof;   //!
   TBranch        *b_mus_tk_chg;   //!
   TBranch        *b_mus_tk_pt;   //!
   TBranch        *b_mus_tk_px;   //!
   TBranch        *b_mus_tk_py;   //!
   TBranch        *b_mus_tk_pz;   //!
   TBranch        *b_mus_tk_eta;   //!
   TBranch        *b_mus_tk_phi;   //!
   TBranch        *b_mus_tk_theta;   //!
   TBranch        *b_mus_tk_d0dum;   //!
   TBranch        *b_mus_tk_dz;   //!
   TBranch        *b_mus_tk_vx;   //!
   TBranch        *b_mus_tk_vy;   //!
   TBranch        *b_mus_tk_vz;   //!
   TBranch        *b_mus_tk_numvalhits;   //!
   TBranch        *b_mus_tk_numlosthits;   //!
   TBranch        *b_mus_tk_d0dumErr;   //!
   TBranch        *b_mus_tk_dzErr;   //!
   TBranch        *b_mus_tk_ptErr;   //!
   TBranch        *b_mus_tk_etaErr;   //!
   TBranch        *b_mus_tk_phiErr;   //!
   TBranch        *b_mus_tk_numvalPixelhits;   //!
   TBranch        *b_mus_tk_numpixelWthMeasr;   //!
   TBranch        *b_mus_stamu_chi2;   //!
   TBranch        *b_mus_stamu_ndof;   //!
   TBranch        *b_mus_stamu_chg;   //!
   TBranch        *b_mus_stamu_pt;   //!
   TBranch        *b_mus_stamu_px;   //!
   TBranch        *b_mus_stamu_py;   //!
   TBranch        *b_mus_stamu_pz;   //!
   TBranch        *b_mus_stamu_eta;   //!
   TBranch        *b_mus_stamu_phi;   //!
   TBranch        *b_mus_stamu_theta;   //!
   TBranch        *b_mus_stamu_d0dum;   //!
   TBranch        *b_mus_stamu_dz;   //!
   TBranch        *b_mus_stamu_vx;   //!
   TBranch        *b_mus_stamu_vy;   //!
   TBranch        *b_mus_stamu_vz;   //!
   TBranch        *b_mus_stamu_numvalhits;   //!
   TBranch        *b_mus_stamu_numlosthits;   //!
   TBranch        *b_mus_stamu_d0dumErr;   //!
   TBranch        *b_mus_stamu_dzErr;   //!
   TBranch        *b_mus_stamu_ptErr;   //!
   TBranch        *b_mus_stamu_etaErr;   //!
   TBranch        *b_mus_stamu_phiErr;   //!
   TBranch        *b_mus_num_matches;   //!
   TBranch        *b_mus_isPFMuon;   //!
   TBranch        *b_mus_isTrackerMuon;   //!
   TBranch        *b_mus_isStandAloneMuon;   //!
   TBranch        *b_mus_isCaloMuon;   //!
   TBranch        *b_mus_isGlobalMuon;   //!
   TBranch        *b_mus_isElectron;   //!
   TBranch        *b_mus_isConvertedPhoton;   //!
   TBranch        *b_mus_isPhoton;   //!
   TBranch        *b_mus_id_All;   //!
   TBranch        *b_mus_id_AllGlobalMuons;   //!
   TBranch        *b_mus_id_AllStandAloneMuons;   //!
   TBranch        *b_mus_id_AllTrackerMuons;   //!
   TBranch        *b_mus_id_TrackerMuonArbitrated;   //!
   TBranch        *b_mus_id_AllArbitrated;   //!
   TBranch        *b_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_mus_id_TMLastStationLoose;   //!
   TBranch        *b_mus_id_TMLastStationTight;   //!
   TBranch        *b_mus_id_TM2DCompatibilityLoose;   //!
   TBranch        *b_mus_id_TM2DCompatibilityTight;   //!
   TBranch        *b_mus_id_TMOneStationLoose;   //!
   TBranch        *b_mus_id_TMOneStationTight;   //!
   TBranch        *b_mus_id_TMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_mus_id_TMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_mus_tk_LayersWithMeasurement;   //!
   TBranch        *b_mus_tk_PixelLayersWithMeasurement;   //!
   TBranch        *b_mus_tk_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_mus_tk_LayersWithoutMeasurement;   //!
   TBranch        *b_mus_tk_ExpectedHitsInner;   //!
   TBranch        *b_mus_tk_ExpectedHitsOuter;   //!
   TBranch        *b_mus_cm_LayersWithMeasurement;   //!
   TBranch        *b_mus_cm_PixelLayersWithMeasurement;   //!
   TBranch        *b_mus_cm_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_mus_cm_LayersWithoutMeasurement;   //!
   TBranch        *b_mus_cm_ExpectedHitsInner;   //!
   TBranch        *b_mus_cm_ExpectedHitsOuter;   //!
   TBranch        *b_mus_picky_LayersWithMeasurement;   //!
   TBranch        *b_mus_picky_PixelLayersWithMeasurement;   //!
   TBranch        *b_mus_picky_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_mus_picky_LayersWithoutMeasurement;   //!
   TBranch        *b_mus_picky_ExpectedHitsInner;   //!
   TBranch        *b_mus_picky_ExpectedHitsOuter;   //!
   TBranch        *b_mus_tpfms_LayersWithMeasurement;   //!
   TBranch        *b_mus_tpfms_PixelLayersWithMeasurement;   //!
   TBranch        *b_mus_tpfms_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_mus_tpfms_LayersWithoutMeasurement;   //!
   TBranch        *b_mus_tpfms_ExpectedHitsInner;   //!
   TBranch        *b_mus_tpfms_ExpectedHitsOuter;   //!
   TBranch        *b_mus_picky_id;   //!
   TBranch        *b_mus_picky_chi2;   //!
   TBranch        *b_mus_picky_ndof;   //!
   TBranch        *b_mus_picky_chg;   //!
   TBranch        *b_mus_picky_pt;   //!
   TBranch        *b_mus_picky_px;   //!
   TBranch        *b_mus_picky_py;   //!
   TBranch        *b_mus_picky_pz;   //!
   TBranch        *b_mus_picky_eta;   //!
   TBranch        *b_mus_picky_phi;   //!
   TBranch        *b_mus_picky_theta;   //!
   TBranch        *b_mus_picky_d0dum;   //!
   TBranch        *b_mus_picky_dz;   //!
   TBranch        *b_mus_picky_vx;   //!
   TBranch        *b_mus_picky_vy;   //!
   TBranch        *b_mus_picky_vz;   //!
   TBranch        *b_mus_picky_numvalhits;   //!
   TBranch        *b_mus_picky_numlosthits;   //!
   TBranch        *b_mus_picky_d0dumErr;   //!
   TBranch        *b_mus_picky_dzErr;   //!
   TBranch        *b_mus_picky_ptErr;   //!
   TBranch        *b_mus_picky_etaErr;   //!
   TBranch        *b_mus_picky_phiErr;   //!
   TBranch        *b_mus_picky_numvalPixelhits;   //!
   TBranch        *b_mus_tpfms_id;   //!
   TBranch        *b_mus_tpfms_chi2;   //!
   TBranch        *b_mus_tpfms_ndof;   //!
   TBranch        *b_mus_tpfms_chg;   //!
   TBranch        *b_mus_tpfms_pt;   //!
   TBranch        *b_mus_tpfms_px;   //!
   TBranch        *b_mus_tpfms_py;   //!
   TBranch        *b_mus_tpfms_pz;   //!
   TBranch        *b_mus_tpfms_eta;   //!
   TBranch        *b_mus_tpfms_phi;   //!
   TBranch        *b_mus_tpfms_theta;   //!
   TBranch        *b_mus_tpfms_d0dum;   //!
   TBranch        *b_mus_tpfms_dz;   //!
   TBranch        *b_mus_tpfms_vx;   //!
   TBranch        *b_mus_tpfms_vy;   //!
   TBranch        *b_mus_tpfms_vz;   //!
   TBranch        *b_mus_tpfms_numvalhits;   //!
   TBranch        *b_mus_tpfms_numlosthits;   //!
   TBranch        *b_mus_tpfms_d0dumErr;   //!
   TBranch        *b_mus_tpfms_dzErr;   //!
   TBranch        *b_mus_tpfms_ptErr;   //!
   TBranch        *b_mus_tpfms_etaErr;   //!
   TBranch        *b_mus_tpfms_phiErr;   //!
   TBranch        *b_mus_tpfms_numvalPixelhits;   //!
   TBranch        *b_mus_dB;   //!
   TBranch        *b_mus_numberOfMatchedStations;   //!
   TBranch        *b_Npfcand;   //!
   TBranch        *b_pfcand_pdgId;   //!
   TBranch        *b_pfcand_pt;   //!
   TBranch        *b_pfcand_pz;   //!
   TBranch        *b_pfcand_px;   //!
   TBranch        *b_pfcand_py;   //!
   TBranch        *b_pfcand_eta;   //!
   TBranch        *b_pfcand_phi;   //!
   TBranch        *b_pfcand_theta;   //!
   TBranch        *b_pfcand_energy;   //!
   TBranch        *b_pfcand_charge;   //!
   TBranch        *b_Nphotons;   //!
   TBranch        *b_photons_energy;   //!
   TBranch        *b_photons_et;   //!
   TBranch        *b_photons_eta;   //!
   TBranch        *b_photons_phi;   //!
   TBranch        *b_photons_pt;   //!
   TBranch        *b_photons_px;   //!
   TBranch        *b_photons_py;   //!
   TBranch        *b_photons_pz;   //!
   TBranch        *b_photons_status;   //!
   TBranch        *b_photons_theta;   //!
   TBranch        *b_photons_hadOverEM;   //!
   TBranch        *b_photons_hadTowOverEM;   //!
   TBranch        *b_photons_scEnergy;   //!
   TBranch        *b_photons_scRawEnergy;   //!
   TBranch        *b_photons_scEta;   //!
   TBranch        *b_photons_scPhi;   //!
   TBranch        *b_photons_scEtaWidth;   //!
   TBranch        *b_photons_scPhiWidth;   //!
   TBranch        *b_photons_tIso;   //!
   TBranch        *b_photons_ecalIso;   //!
   TBranch        *b_photons_hcalIso;   //!
   TBranch        *b_photons_isoEcalRecHitDR04;   //!
   TBranch        *b_photons_isoHcalRecHitDR04;   //!
   TBranch        *b_photons_isoSolidTrkConeDR04;   //!
   TBranch        *b_photons_isoHollowTrkConeDR04;   //!
   TBranch        *b_photons_nTrkSolidConeDR04;   //!
   TBranch        *b_photons_nTrkHollowConeDR04;   //!
   TBranch        *b_photons_isoEcalRecHitDR03;   //!
   TBranch        *b_photons_isoHcalRecHitDR03;   //!
   TBranch        *b_photons_isoSolidTrkConeDR03;   //!
   TBranch        *b_photons_isoHollowTrkConeDR03;   //!
   TBranch        *b_photons_nTrkSolidConeDR03;   //!
   TBranch        *b_photons_nTrkHollowConeDR03;   //!
   TBranch        *b_photons_isAlsoElectron;   //!
   TBranch        *b_photons_hasPixelSeed;   //!
   TBranch        *b_photons_isConverted;   //!
   TBranch        *b_photons_isEBGap;   //!
   TBranch        *b_photons_isEEGap;   //!
   TBranch        *b_photons_isEBEEGap;   //!
   TBranch        *b_photons_isEBPho;   //!
   TBranch        *b_photons_isEEPho;   //!
   TBranch        *b_photons_isLoosePhoton;   //!
   TBranch        *b_photons_isTightPhoton;   //!
   TBranch        *b_photons_maxEnergyXtal;   //!
   TBranch        *b_photons_e1x5;   //!
   TBranch        *b_photons_e2x5;   //!
   TBranch        *b_photons_e3x3;   //!
   TBranch        *b_photons_e5x5;   //!
   TBranch        *b_photons_sigmaEtaEta;   //!
   TBranch        *b_photons_sigmaIetaIeta;   //!
   TBranch        *b_photons_r9;   //!
   TBranch        *b_photons_gen_et;   //!
   TBranch        *b_photons_gen_eta;   //!
   TBranch        *b_photons_gen_phi;   //!
   TBranch        *b_photons_gen_id;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_xErr;   //!
   TBranch        *b_pv_yErr;   //!
   TBranch        *b_pv_zErr;   //!
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_pv_isValid;   //!
   TBranch        *b_pv_tracksSize;   //!
   TBranch        *b_Ntaus;   //!
   TBranch        *b_taus_status;   //!
   TBranch        *b_taus_phi;   //!
   TBranch        *b_taus_pt;   //!
   TBranch        *b_taus_pz;   //!
   TBranch        *b_taus_px;   //!
   TBranch        *b_taus_py;   //!
   TBranch        *b_taus_eta;   //!
   TBranch        *b_taus_theta;   //!
   TBranch        *b_taus_et;   //!
   TBranch        *b_taus_energy;   //!
   TBranch        *b_taus_charge;   //!
   TBranch        *b_taus_emf;   //!
   TBranch        *b_taus_hcalTotOverPLead;   //!
   TBranch        *b_taus_hcalMaxOverPLead;   //!
   TBranch        *b_taus_hcal3x3OverPLead;   //!
   TBranch        *b_taus_ecalStripSumEOverPLead;   //!
   TBranch        *b_taus_elecPreIdOutput;   //!
   TBranch        *b_taus_elecPreIdDecision;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_pt;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_charge;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_eta;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_ECAL_eta;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_phi;   //!
   TBranch        *b_taus_isoPFGammaCandsEtSum;   //!
   TBranch        *b_taus_isoPFChargedHadrCandsPtSum;   //!
   TBranch        *b_taus_leadingTrackFinding;   //!
   TBranch        *b_taus_leadingTrackPtCut;   //!
   TBranch        *b_taus_trackIsolation;   //!
   TBranch        *b_taus_ecalIsolation;   //!
   TBranch        *b_taus_byIsolation;   //!
   TBranch        *b_taus_againstElectron;   //!
   TBranch        *b_taus_againstMuon;   //!
   TBranch        *b_taus_taNC_quarter;   //!
   TBranch        *b_taus_taNC_one;   //!
   TBranch        *b_taus_taNC_half;   //!
   TBranch        *b_taus_taNC_tenth;   //!
   TBranch        *b_taus_taNC;   //!
   TBranch        *b_taus_byIsoUsingLeadingPi;   //!
   TBranch        *b_taus_tkIsoUsingLeadingPi;   //!
   TBranch        *b_taus_ecalIsoUsingLeadingPi;   //!
   TBranch        *b_taus_againstElectronLoose;   //!
   TBranch        *b_taus_againstElectronMedium;   //!
   TBranch        *b_taus_againstElectronTight;   //!
   TBranch        *b_taus_againstElectronMVA;   //!
   TBranch        *b_taus_againstMuonLoose;   //!
   TBranch        *b_taus_againstMuonMedium;   //!
   TBranch        *b_taus_againstMuonTight;   //!
   TBranch        *b_taus_decayModeFinding;   //!
   TBranch        *b_taus_byVLooseIsolation;   //!
   TBranch        *b_taus_byLooseIsolation;   //!
   TBranch        *b_taus_byMediumIsolation;   //!
   TBranch        *b_taus_byTightIsolation;   //!
   TBranch        *b_taus_byVLooseIsolationDeltaBetaCorr;   //!
   TBranch        *b_taus_byLooseIsolationDeltaBetaCorr;   //!
   TBranch        *b_taus_byMediumIsolationDeltaBetaCorr;   //!
   TBranch        *b_taus_byTightIsolationDeltaBetaCorr;   //!
   TBranch        *b_taus_signalPFChargedHadrCandsSize;   //!
   TBranch        *b_taus_muDecision;   //!
   TBranch        *b_taus_Nprongs;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_experimentType;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_model_params;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_px;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_py;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_pz;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_energy;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_phi;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_eta;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_index;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT10_nconstituents;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_px;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_py;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_pz;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_energy;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_phi;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_eta;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_index;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT15_nconstituents;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_px;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_py;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_pz;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_energy;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_phi;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_eta;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_index;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT20_nconstituents;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_px;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_py;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_pz;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_energy;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_phi;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_eta;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_index;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT25_nconstituents;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_px;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_py;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_pz;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_energy;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_phi;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_eta;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_index;   //!
   TBranch        *b_fastjets_AK4_R1p2_R0p5pT30_nconstituents;   //!

void InitializeA(TChain *fChain)
{

   // Set object pointer
   trigger_decision = 0;
   trigger_name = 0;
   trigger_prescalevalue = 0;
   standalone_triggerobject_pt = 0;
   standalone_triggerobject_px = 0;
   standalone_triggerobject_py = 0;
   standalone_triggerobject_pz = 0;
   standalone_triggerobject_et = 0;
   standalone_triggerobject_energy = 0;
   standalone_triggerobject_phi = 0;
   standalone_triggerobject_eta = 0;
   standalone_triggerobject_collectionname = 0;
   PU_zpositions = 0;
   PU_sumpT_lowpT = 0;
   PU_sumpT_highpT = 0;
   PU_ntrks_lowpT = 0;
   PU_ntrks_highpT = 0;
   PU_NumInteractions = 0;
   PU_bunchCrossing = 0;
   PU_TrueNumInteractions = 0;
   els_isPF = 0;
   mus_isPF = 0;
   jets_AK4_maxpt_id = 0;
   jets_AK4_mu_ind = 0;
   jets_AK4_el_ind = 0;
   taus_el_ind = 0;
   taus_mu_ind = 0;
   els_jet_ind = 0;
   mus_jet_ind = 0;
   fChain->SetBranchAddress("trigger_decision", &trigger_decision, &b_trigger_decision);
   fChain->SetBranchAddress("trigger_name", &trigger_name, &b_trigger_name);
   fChain->SetBranchAddress("trigger_prescalevalue", &trigger_prescalevalue, &b_trigger_prescalevalue);
   fChain->SetBranchAddress("standalone_triggerobject_pt", &standalone_triggerobject_pt, &b_standalone_triggerobject_pt);
   fChain->SetBranchAddress("standalone_triggerobject_px", &standalone_triggerobject_px, &b_standalone_triggerobject_px);
   fChain->SetBranchAddress("standalone_triggerobject_py", &standalone_triggerobject_py, &b_standalone_triggerobject_py);
   fChain->SetBranchAddress("standalone_triggerobject_pz", &standalone_triggerobject_pz, &b_standalone_triggerobject_pz);
   fChain->SetBranchAddress("standalone_triggerobject_et", &standalone_triggerobject_et, &b_standalone_triggerobject_et);
   fChain->SetBranchAddress("standalone_triggerobject_energy", &standalone_triggerobject_energy, &b_standalone_triggerobject_energy);
   fChain->SetBranchAddress("standalone_triggerobject_phi", &standalone_triggerobject_phi, &b_standalone_triggerobject_phi);
   fChain->SetBranchAddress("standalone_triggerobject_eta", &standalone_triggerobject_eta, &b_standalone_triggerobject_eta);
   fChain->SetBranchAddress("standalone_triggerobject_collectionname", &standalone_triggerobject_collectionname, &b_standalone_triggerobject_collectionname);
   fChain->SetBranchAddress("PU_zpositions", &PU_zpositions, &b_PU_zpositions);
   fChain->SetBranchAddress("PU_sumpT_lowpT", &PU_sumpT_lowpT, &b_PU_sumpT_lowpT);
   fChain->SetBranchAddress("PU_sumpT_highpT", &PU_sumpT_highpT, &b_PU_sumpT_highpT);
   fChain->SetBranchAddress("PU_ntrks_lowpT", &PU_ntrks_lowpT, &b_PU_ntrks_lowpT);
   fChain->SetBranchAddress("PU_ntrks_highpT", &PU_ntrks_highpT, &b_PU_ntrks_highpT);
   fChain->SetBranchAddress("PU_NumInteractions", &PU_NumInteractions, &b_PU_NumInteractions);
   fChain->SetBranchAddress("PU_bunchCrossing", &PU_bunchCrossing, &b_PU_bunchCrossing);
   fChain->SetBranchAddress("PU_TrueNumInteractions", &PU_TrueNumInteractions, &b_PU_TrueNumInteractions);
   fChain->SetBranchAddress("trackingfailurefilter_decision", &trackingfailurefilter_decision, &b_rackingfailurefilter_decision);
   fChain->SetBranchAddress("goodVerticesfilter_decision", &goodVerticesfilter_decision, &b_oodVerticesfilter_decision);
   fChain->SetBranchAddress("cschalofilter_decision", &cschalofilter_decision, &b_schalofilter_decision);
   fChain->SetBranchAddress("trkPOGfilter_decision", &trkPOGfilter_decision, &b_rkPOGfilter_decision);
   fChain->SetBranchAddress("trkPOG_logErrorTooManyClustersfilter_decision", &trkPOG_logErrorTooManyClustersfilter_decision, &b_rkPOG_logErrorTooManyClustersfilter_decision);
   fChain->SetBranchAddress("EcalDeadCellTriggerPrimitivefilter_decision", &EcalDeadCellTriggerPrimitivefilter_decision, &b_calDeadCellTriggerPrimitivefilter_decision);
   fChain->SetBranchAddress("ecallaserfilter_decision", &ecallaserfilter_decision, &b_ecallaserfilter_decision);
   fChain->SetBranchAddress("trkPOG_manystripclus53Xfilter_decision", &trkPOG_manystripclus53Xfilter_decision, &b_rkPOG_manystripclus53Xfilter_decision);
   fChain->SetBranchAddress("eebadscfilter_decision", &eebadscfilter_decision, &b_ebadscfilter_decision);
   fChain->SetBranchAddress("METFiltersfilter_decision", &METFiltersfilter_decision, &b_ETFiltersfilter_decision);
   fChain->SetBranchAddress("HBHENoisefilter_decision", &HBHENoisefilter_decision, &b_BHENoisefilter_decision);
   fChain->SetBranchAddress("trkPOG_toomanystripclus53Xfilter_decision", &trkPOG_toomanystripclus53Xfilter_decision, &b_rkPOG_toomanystripclus53Xfilter_decision);
   fChain->SetBranchAddress("hcallaserfilter_decision", &hcallaserfilter_decision, &b_hcallaserfilter_decision);
   fChain->SetBranchAddress("els_isPF", &els_isPF, &b_els_isPF);
   fChain->SetBranchAddress("mus_isPF", &mus_isPF, &b_mus_isPF);
   fChain->SetBranchAddress("jets_AK4_maxpt_id", &jets_AK4_maxpt_id, &b_jets_AK4_maxpt_id);
   fChain->SetBranchAddress("jets_AK4_mu_ind", &jets_AK4_mu_ind, &b_jets_AK4_mu_ind);
   fChain->SetBranchAddress("jets_AK4_el_ind", &jets_AK4_el_ind, &b_jets_AK4_el_ind);
   fChain->SetBranchAddress("taus_el_ind", &taus_el_ind, &b_taus_el_ind);
   fChain->SetBranchAddress("taus_mu_ind", &taus_mu_ind, &b_taus_mu_ind);
   fChain->SetBranchAddress("els_jet_ind", &els_jet_ind, &b_els_jet_ind);
   fChain->SetBranchAddress("mus_jet_ind", &mus_jet_ind, &b_mus_jet_ind);

}

void InitializeB(TChain *fChain)
{

   // Set object pointer
   beamSpot_x = 0;
   beamSpot_y = 0;
   beamSpot_z = 0;
   beamSpot_x0Error = 0;
   beamSpot_y0Error = 0;
   beamSpot_z0Error = 0;
   beamSpot_sigmaZ = 0;
   beamSpot_sigmaZ0Error = 0;
   beamSpot_dxdz = 0;
   beamSpot_dxdzError = 0;
   beamSpot_dydz = 0;
   beamSpot_dydzError = 0;
   beamSpot_beamWidthX = 0;
   beamSpot_beamWidthY = 0;
   beamSpot_beamWidthXError = 0;
   beamSpot_beamWidthYError = 0;
   els_energy = 0;
   els_et = 0;
   els_eta = 0;
   els_phi = 0;
   els_pt = 0;
   els_px = 0;
   els_py = 0;
   els_pz = 0;
   els_status = 0;
   els_theta = 0;
   els_pfIsolationR03_sumChargedHadronPt = 0;
   els_pfIsolationR03_sumNeutralHadronEt = 0;
   els_pfIsolationR03_sumPhotonEt = 0;
   els_pfIsolationR03_sumPUPt = 0;
   els_full5x5_sigmaIetaIeta = 0;
   els_gen_id = 0;
   els_gen_phi = 0;
   els_gen_pt = 0;
   els_gen_pz = 0;
   els_gen_px = 0;
   els_gen_py = 0;
   els_gen_eta = 0;
   els_gen_theta = 0;
   els_gen_et = 0;
   els_gen_mother_id = 0;
   els_gen_mother_phi = 0;
   els_gen_mother_pt = 0;
   els_gen_mother_pz = 0;
   els_gen_mother_px = 0;
   els_gen_mother_py = 0;
   els_gen_mother_eta = 0;
   els_gen_mother_theta = 0;
   els_gen_mother_et = 0;
   els_tightId = 0;
   els_looseId = 0;
   els_robustTightId = 0;
   els_robustLooseId = 0;
   els_robustHighEnergyId = 0;
   els_simpleEleId95relIso = 0;
   els_simpleEleId90relIso = 0;
   els_simpleEleId85relIso = 0;
   els_simpleEleId80relIso = 0;
   els_simpleEleId70relIso = 0;
   els_simpleEleId60relIso = 0;
   els_simpleEleId95cIso = 0;
   els_simpleEleId90cIso = 0;
   els_simpleEleId85cIso = 0;
   els_simpleEleId80cIso = 0;
   els_simpleEleId70cIso = 0;
   els_simpleEleId60cIso = 0;
   els_cIso = 0;
   els_tIso = 0;
   els_ecalIso = 0;
   els_hcalIso = 0;
   els_chi2 = 0;
   els_charge = 0;
   els_caloEnergy = 0;
   els_hadOverEm = 0;
   els_hcalOverEcalBc = 0;
   els_eOverPIn = 0;
   els_eSeedOverPOut = 0;
   els_sigmaEtaEta = 0;
   els_sigmaIEtaIEta = 0;
   els_scEnergy = 0;
   els_scRawEnergy = 0;
   els_scSeedEnergy = 0;
   els_scEta = 0;
   els_scPhi = 0;
   els_scEtaWidth = 0;
   els_scPhiWidth = 0;
   els_scE1x5 = 0;
   els_scE2x5Max = 0;
   els_scE5x5 = 0;
   els_isEB = 0;
   els_isEE = 0;
   els_dEtaIn = 0;
   els_dPhiIn = 0;
   els_dEtaOut = 0;
   els_dPhiOut = 0;
   els_numvalhits = 0;
   els_numlosthits = 0;
   els_basicClustersSize = 0;
   els_tk_pz = 0;
   els_tk_pt = 0;
   els_tk_phi = 0;
   els_tk_eta = 0;
   els_d0dum = 0;
   els_dz = 0;
   els_vx = 0;
   els_vy = 0;
   els_vz = 0;
   els_ndof = 0;
   els_ptError = 0;
   els_d0dumError = 0;
   els_dzError = 0;
   els_etaError = 0;
   els_phiError = 0;
   els_tk_charge = 0;
   els_core_ecalDrivenSeed = 0;
   els_n_inner_layer = 0;
   els_n_outer_layer = 0;
   els_ctf_tk_id = 0;
   els_ctf_tk_charge = 0;
   els_ctf_tk_eta = 0;
   els_ctf_tk_phi = 0;
   els_fbrem = 0;
   els_shFracInnerHits = 0;
   els_dr03EcalRecHitSumEt = 0;
   els_dr03HcalTowerSumEt = 0;
   els_dr03HcalDepth1TowerSumEt = 0;
   els_dr03HcalDepth2TowerSumEt = 0;
   els_dr03TkSumPt = 0;
   els_dr04EcalRecHitSumEt = 0;
   els_dr04HcalTowerSumEt = 0;
   els_dr04HcalDepth1TowerSumEt = 0;
   els_dr04HcalDepth2TowerSumEt = 0;
   els_dr04TkSumPt = 0;
   els_cpx = 0;
   els_cpy = 0;
   els_cpz = 0;
   els_vpx = 0;
   els_vpy = 0;
   els_vpz = 0;
   els_cx = 0;
   els_cy = 0;
   els_cz = 0;
   els_PATpassConversionVeto = 0;
   jets_AK4_status = 0;
   jets_AK4_phi = 0;
   jets_AK4_pt = 0;
   jets_AK4_pz = 0;
   jets_AK4_px = 0;
   jets_AK4_py = 0;
   jets_AK4_eta = 0;
   jets_AK4_theta = 0;
   jets_AK4_et = 0;
   jets_AK4_energy = 0;
   jets_AK4_parton_Id = 0;
   jets_AK4_parton_motherId = 0;
   jets_AK4_parton_pt = 0;
   jets_AK4_parton_phi = 0;
   jets_AK4_parton_eta = 0;
   jets_AK4_parton_Energy = 0;
   jets_AK4_parton_mass = 0;
   jets_AK4_gen_et = 0;
   jets_AK4_gen_pt = 0;
   jets_AK4_gen_eta = 0;
   jets_AK4_gen_phi = 0;
   jets_AK4_gen_mass = 0;
   jets_AK4_gen_Energy = 0;
   jets_AK4_gen_Id = 0;
   jets_AK4_gen_motherID = 0;
   jets_AK4_gen_threeCharge = 0;
   jets_AK4_partonFlavour = 0;
   jets_AK4_btag_TC_highPur = 0;
   jets_AK4_btag_TC_highEff = 0;
   jets_AK4_btag_jetProb = 0;
   jets_AK4_btag_jetBProb = 0;
   jets_AK4_btag_softEle = 0;
   jets_AK4_btag_softMuon = 0;
   jets_AK4_btag_secVertexHighPur = 0;
   jets_AK4_btag_secVertexHighEff = 0;
   jets_AK4_btag_secVertexCombined = 0;
   jets_AK4_jetCharge = 0;
   jets_AK4_chgEmE = 0;
   jets_AK4_chgHadE = 0;
   jets_AK4_photonEnergy = 0;
   jets_AK4_chgMuE = 0;
   jets_AK4_chg_Mult = 0;
   jets_AK4_neutralEmE = 0;
   jets_AK4_neutralHadE = 0;
   jets_AK4_neutral_Mult = 0;
   jets_AK4_mu_Mult = 0;
   jets_AK4_emf = 0;
   jets_AK4_ehf = 0;
   jets_AK4_n60 = 0;
   jets_AK4_n90 = 0;
   jets_AK4_etaetaMoment = 0;
   jets_AK4_etaphiMoment = 0;
   jets_AK4_phiphiMoment = 0;
   jets_AK4_n90Hits = 0;
   jets_AK4_fHPD = 0;
   jets_AK4_fRBX = 0;
   jets_AK4_hitsInN90 = 0;
   jets_AK4_nECALTowers = 0;
   jets_AK4_nHCALTowers = 0;
   jets_AK4_fSubDetector1 = 0;
   jets_AK4_fSubDetector2 = 0;
   jets_AK4_fSubDetector3 = 0;
   jets_AK4_fSubDetector4 = 0;
   jets_AK4_area = 0;
   jets_AK4_corrFactorRaw = 0;
   jets_AK4_rawPt = 0;
   jets_AK4_mass = 0;
   mc_doc_id = 0;
   mc_doc_pt = 0;
   mc_doc_px = 0;
   mc_doc_py = 0;
   mc_doc_pz = 0;
   mc_doc_eta = 0;
   mc_doc_phi = 0;
   mc_doc_theta = 0;
   mc_doc_energy = 0;
   mc_doc_status = 0;
   mc_doc_charge = 0;
   mc_doc_mother_id = 0;
   mc_doc_grandmother_id = 0;
   mc_doc_ggrandmother_id = 0;
   mc_doc_mother_pt = 0;
   mc_doc_vertex_x = 0;
   mc_doc_vertex_y = 0;
   mc_doc_vertex_z = 0;
   mc_doc_mass = 0;
   mc_doc_numOfDaughters = 0;
   mc_doc_numOfMothers = 0;
   mc_electrons_id = 0;
   mc_electrons_pt = 0;
   mc_electrons_px = 0;
   mc_electrons_py = 0;
   mc_electrons_pz = 0;
   mc_electrons_eta = 0;
   mc_electrons_phi = 0;
   mc_electrons_theta = 0;
   mc_electrons_status = 0;
   mc_electrons_energy = 0;
   mc_electrons_charge = 0;
   mc_electrons_mother_id = 0;
   mc_electrons_mother_pt = 0;
   mc_electrons_grandmother_id = 0;
   mc_electrons_ggrandmother_id = 0;
   mc_electrons_vertex_x = 0;
   mc_electrons_vertex_y = 0;
   mc_electrons_vertex_z = 0;
   mc_electrons_mass = 0;
   mc_electrons_numOfDaughters = 0;
   mc_final_id = 0;
   mc_final_pt = 0;
   mc_final_px = 0;
   mc_final_py = 0;
   mc_final_pz = 0;
   mc_final_eta = 0;
   mc_final_phi = 0;
   mc_final_theta = 0;
   mc_final_energy = 0;
   mc_final_status = 0;
   mc_final_charge = 0;
   mc_final_mother_id = 0;
   mc_final_grandmother_id = 0;
   mc_final_ggrandmother_id = 0;
   mc_final_mother_pt = 0;
   mc_final_vertex_x = 0;
   mc_final_vertex_y = 0;
   mc_final_vertex_z = 0;
   mc_final_mass = 0;
   mc_final_numOfDaughters = 0;
   mc_final_numOfMothers = 0;
   mc_jets_phi = 0;
   mc_jets_pt = 0;
   mc_jets_pz = 0;
   mc_jets_px = 0;
   mc_jets_py = 0;
   mc_jets_eta = 0;
   mc_jets_theta = 0;
   mc_jets_et = 0;
   mc_jets_energy = 0;
   mc_jets_emEnergy = 0;
   mc_jets_hadEnergy = 0;
   mc_jets_invisibleEnergy = 0;
   mc_jets_auxiliaryEnergy = 0;
   mc_jets_etaetaMoment = 0;
   mc_jets_etaphiMoment = 0;
   mc_jets_phiphiMoment = 0;
   mc_jets_mass = 0;
   mc_mus_id = 0;
   mc_mus_pt = 0;
   mc_mus_px = 0;
   mc_mus_py = 0;
   mc_mus_pz = 0;
   mc_mus_eta = 0;
   mc_mus_phi = 0;
   mc_mus_theta = 0;
   mc_mus_status = 0;
   mc_mus_energy = 0;
   mc_mus_charge = 0;
   mc_mus_mother_id = 0;
   mc_mus_mother_pt = 0;
   mc_mus_grandmother_id = 0;
   mc_mus_ggrandmother_id = 0;
   mc_mus_vertex_x = 0;
   mc_mus_vertex_y = 0;
   mc_mus_vertex_z = 0;
   mc_mus_mass = 0;
   mc_mus_numOfDaughters = 0;
   mc_nues_id = 0;
   mc_nues_pt = 0;
   mc_nues_px = 0;
   mc_nues_py = 0;
   mc_nues_pz = 0;
   mc_nues_eta = 0;
   mc_nues_phi = 0;
   mc_nues_theta = 0;
   mc_nues_status = 0;
   mc_nues_energy = 0;
   mc_nues_charge = 0;
   mc_nues_mother_id = 0;
   mc_nues_mother_pt = 0;
   mc_nues_grandmother_id = 0;
   mc_nues_ggrandmother_id = 0;
   mc_nues_vertex_x = 0;
   mc_nues_vertex_y = 0;
   mc_nues_vertex_z = 0;
   mc_nues_mass = 0;
   mc_nues_numOfDaughters = 0;
   mc_numus_id = 0;
   mc_numus_pt = 0;
   mc_numus_px = 0;
   mc_numus_py = 0;
   mc_numus_pz = 0;
   mc_numus_eta = 0;
   mc_numus_phi = 0;
   mc_numus_theta = 0;
   mc_numus_status = 0;
   mc_numus_energy = 0;
   mc_numus_charge = 0;
   mc_numus_mother_id = 0;
   mc_numus_mother_pt = 0;
   mc_numus_grandmother_id = 0;
   mc_numus_ggrandmother_id = 0;
   mc_numus_vertex_x = 0;
   mc_numus_vertex_y = 0;
   mc_numus_vertex_z = 0;
   mc_numus_mass = 0;
   mc_numus_numOfDaughters = 0;
   mc_nutaus_id = 0;
   mc_nutaus_pt = 0;
   mc_nutaus_px = 0;
   mc_nutaus_py = 0;
   mc_nutaus_pz = 0;
   mc_nutaus_eta = 0;
   mc_nutaus_phi = 0;
   mc_nutaus_theta = 0;
   mc_nutaus_status = 0;
   mc_nutaus_energy = 0;
   mc_nutaus_charge = 0;
   mc_nutaus_mother_id = 0;
   mc_nutaus_mother_pt = 0;
   mc_nutaus_grandmother_id = 0;
   mc_nutaus_ggrandmother_id = 0;
   mc_nutaus_vertex_x = 0;
   mc_nutaus_vertex_y = 0;
   mc_nutaus_vertex_z = 0;
   mc_nutaus_mass = 0;
   mc_nutaus_numOfDaughters = 0;
   mc_photons_id = 0;
   mc_photons_pt = 0;
   mc_photons_px = 0;
   mc_photons_py = 0;
   mc_photons_pz = 0;
   mc_photons_eta = 0;
   mc_photons_phi = 0;
   mc_photons_theta = 0;
   mc_photons_status = 0;
   mc_photons_energy = 0;
   mc_photons_charge = 0;
   mc_photons_mother_id = 0;
   mc_photons_mother_pt = 0;
   mc_photons_grandmother_id = 0;
   mc_photons_ggrandmother_id = 0;
   mc_photons_vertex_x = 0;
   mc_photons_vertex_y = 0;
   mc_photons_vertex_z = 0;
   mc_photons_mass = 0;
   mc_photons_numOfDaughters = 0;
   mc_taus_id = 0;
   mc_taus_pt = 0;
   mc_taus_px = 0;
   mc_taus_py = 0;
   mc_taus_pz = 0;
   mc_taus_eta = 0;
   mc_taus_phi = 0;
   mc_taus_theta = 0;
   mc_taus_status = 0;
   mc_taus_energy = 0;
   mc_taus_charge = 0;
   mc_taus_mother_id = 0;
   mc_taus_mother_pt = 0;
   mc_taus_grandmother_id = 0;
   mc_taus_ggrandmother_id = 0;
   mc_taus_vertex_x = 0;
   mc_taus_vertex_y = 0;
   mc_taus_vertex_z = 0;
   mc_taus_mass = 0;
   mc_taus_numOfDaughters = 0;
   mets_et = 0;
   mets_phi = 0;
   mets_ex = 0;
   mets_ey = 0;
   mets_gen_et = 0;
   mets_gen_phi = 0;
   mets_sign = 0;
   mets_sumEt = 0;
   mets_unCPhi = 0;
   mets_unCPt = 0;
   mus_energy = 0;
   mus_et = 0;
   mus_eta = 0;
   mus_phi = 0;
   mus_pt = 0;
   mus_px = 0;
   mus_py = 0;
   mus_pz = 0;
   mus_status = 0;
   mus_theta = 0;
   mus_gen_id = 0;
   mus_gen_phi = 0;
   mus_gen_pt = 0;
   mus_gen_pz = 0;
   mus_gen_px = 0;
   mus_gen_py = 0;
   mus_gen_eta = 0;
   mus_gen_theta = 0;
   mus_gen_et = 0;
   mus_gen_mother_id = 0;
   mus_gen_mother_phi = 0;
   mus_gen_mother_pt = 0;
   mus_gen_mother_pz = 0;
   mus_gen_mother_px = 0;
   mus_gen_mother_py = 0;
   mus_gen_mother_eta = 0;
   mus_gen_mother_theta = 0;
   mus_gen_mother_et = 0;
   mus_tkHits = 0;
   mus_cIso = 0;
   mus_tIso = 0;
   mus_ecalIso = 0;
   mus_hcalIso = 0;
   mus_ecalvetoDep = 0;
   mus_hcalvetoDep = 0;
   mus_calEnergyEm = 0;
   mus_calEnergyHad = 0;
   mus_calEnergyHo = 0;
   mus_calEnergyEmS9 = 0;
   mus_calEnergyHadS9 = 0;
   mus_calEnergyHoS9 = 0;
   mus_iso03_emVetoEt = 0;
   mus_iso03_hadVetoEt = 0;
   mus_iso03_sumPt = 0;
   mus_iso03_emEt = 0;
   mus_iso03_hadEt = 0;
   mus_iso03_hoEt = 0;
   mus_iso03_nTracks = 0;
   mus_iso05_sumPt = 0;
   mus_iso05_emEt = 0;
   mus_iso05_hadEt = 0;
   mus_iso05_hoEt = 0;
   mus_iso05_nTracks = 0;
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
   mus_cm_chg = 0;
   mus_cm_pt = 0;
   mus_cm_px = 0;
   mus_cm_py = 0;
   mus_cm_pz = 0;
   mus_cm_eta = 0;
   mus_cm_phi = 0;
   mus_cm_theta = 0;
   mus_cm_d0dum = 0;
   mus_cm_dz = 0;
   mus_cm_vx = 0;
   mus_cm_vy = 0;
   mus_cm_vz = 0;
   mus_cm_numvalhits = 0;
   mus_cm_numlosthits = 0;
   mus_cm_numvalMuonhits = 0;
   mus_cm_d0dumErr = 0;
   mus_cm_dzErr = 0;
   mus_cm_ptErr = 0;
   mus_cm_etaErr = 0;
   mus_cm_phiErr = 0;
   mus_tk_id = 0;
   mus_tk_chi2 = 0;
   mus_tk_ndof = 0;
   mus_tk_chg = 0;
   mus_tk_pt = 0;
   mus_tk_px = 0;
   mus_tk_py = 0;
   mus_tk_pz = 0;
   mus_tk_eta = 0;
   mus_tk_phi = 0;
   mus_tk_theta = 0;
   mus_tk_d0dum = 0;
   mus_tk_dz = 0;
   mus_tk_vx = 0;
   mus_tk_vy = 0;
   mus_tk_vz = 0;
   mus_tk_numvalhits = 0;
   mus_tk_numlosthits = 0;
   mus_tk_d0dumErr = 0;
   mus_tk_dzErr = 0;
   mus_tk_ptErr = 0;
   mus_tk_etaErr = 0;
   mus_tk_phiErr = 0;
   mus_tk_numvalPixelhits = 0;
   mus_tk_numpixelWthMeasr = 0;
   mus_stamu_chi2 = 0;
   mus_stamu_ndof = 0;
   mus_stamu_chg = 0;
   mus_stamu_pt = 0;
   mus_stamu_px = 0;
   mus_stamu_py = 0;
   mus_stamu_pz = 0;
   mus_stamu_eta = 0;
   mus_stamu_phi = 0;
   mus_stamu_theta = 0;
   mus_stamu_d0dum = 0;
   mus_stamu_dz = 0;
   mus_stamu_vx = 0;
   mus_stamu_vy = 0;
   mus_stamu_vz = 0;
   mus_stamu_numvalhits = 0;
   mus_stamu_numlosthits = 0;
   mus_stamu_d0dumErr = 0;
   mus_stamu_dzErr = 0;
   mus_stamu_ptErr = 0;
   mus_stamu_etaErr = 0;
   mus_stamu_phiErr = 0;
   mus_num_matches = 0;
   mus_isPFMuon = 0;
   mus_isTrackerMuon = 0;
   mus_isStandAloneMuon = 0;
   mus_isCaloMuon = 0;
   mus_isGlobalMuon = 0;
   mus_isElectron = 0;
   mus_isConvertedPhoton = 0;
   mus_isPhoton = 0;
   mus_id_All = 0;
   mus_id_AllGlobalMuons = 0;
   mus_id_AllStandAloneMuons = 0;
   mus_id_AllTrackerMuons = 0;
   mus_id_TrackerMuonArbitrated = 0;
   mus_id_AllArbitrated = 0;
   mus_id_GlobalMuonPromptTight = 0;
   mus_id_TMLastStationLoose = 0;
   mus_id_TMLastStationTight = 0;
   mus_id_TM2DCompatibilityLoose = 0;
   mus_id_TM2DCompatibilityTight = 0;
   mus_id_TMOneStationLoose = 0;
   mus_id_TMOneStationTight = 0;
   mus_id_TMLastStationOptimizedLowPtLoose = 0;
   mus_id_TMLastStationOptimizedLowPtTight = 0;
   mus_tk_LayersWithMeasurement = 0;
   mus_tk_PixelLayersWithMeasurement = 0;
   mus_tk_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_tk_LayersWithoutMeasurement = 0;
   mus_tk_ExpectedHitsInner = 0;
   mus_tk_ExpectedHitsOuter = 0;
   mus_cm_LayersWithMeasurement = 0;
   mus_cm_PixelLayersWithMeasurement = 0;
   mus_cm_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_cm_LayersWithoutMeasurement = 0;
   mus_cm_ExpectedHitsInner = 0;
   mus_cm_ExpectedHitsOuter = 0;
   mus_picky_LayersWithMeasurement = 0;
   mus_picky_PixelLayersWithMeasurement = 0;
   mus_picky_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_picky_LayersWithoutMeasurement = 0;
   mus_picky_ExpectedHitsInner = 0;
   mus_picky_ExpectedHitsOuter = 0;
   mus_tpfms_LayersWithMeasurement = 0;
   mus_tpfms_PixelLayersWithMeasurement = 0;
   mus_tpfms_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_tpfms_LayersWithoutMeasurement = 0;
   mus_tpfms_ExpectedHitsInner = 0;
   mus_tpfms_ExpectedHitsOuter = 0;
   mus_picky_id = 0;
   mus_picky_chi2 = 0;
   mus_picky_ndof = 0;
   mus_picky_chg = 0;
   mus_picky_pt = 0;
   mus_picky_px = 0;
   mus_picky_py = 0;
   mus_picky_pz = 0;
   mus_picky_eta = 0;
   mus_picky_phi = 0;
   mus_picky_theta = 0;
   mus_picky_d0dum = 0;
   mus_picky_dz = 0;
   mus_picky_vx = 0;
   mus_picky_vy = 0;
   mus_picky_vz = 0;
   mus_picky_numvalhits = 0;
   mus_picky_numlosthits = 0;
   mus_picky_d0dumErr = 0;
   mus_picky_dzErr = 0;
   mus_picky_ptErr = 0;
   mus_picky_etaErr = 0;
   mus_picky_phiErr = 0;
   mus_picky_numvalPixelhits = 0;
   mus_tpfms_id = 0;
   mus_tpfms_chi2 = 0;
   mus_tpfms_ndof = 0;
   mus_tpfms_chg = 0;
   mus_tpfms_pt = 0;
   mus_tpfms_px = 0;
   mus_tpfms_py = 0;
   mus_tpfms_pz = 0;
   mus_tpfms_eta = 0;
   mus_tpfms_phi = 0;
   mus_tpfms_theta = 0;
   mus_tpfms_d0dum = 0;
   mus_tpfms_dz = 0;
   mus_tpfms_vx = 0;
   mus_tpfms_vy = 0;
   mus_tpfms_vz = 0;
   mus_tpfms_numvalhits = 0;
   mus_tpfms_numlosthits = 0;
   mus_tpfms_d0dumErr = 0;
   mus_tpfms_dzErr = 0;
   mus_tpfms_ptErr = 0;
   mus_tpfms_etaErr = 0;
   mus_tpfms_phiErr = 0;
   mus_tpfms_numvalPixelhits = 0;
   mus_dB = 0;
   mus_numberOfMatchedStations = 0;
   pfcand_pdgId = 0;
   pfcand_pt = 0;
   pfcand_pz = 0;
   pfcand_px = 0;
   pfcand_py = 0;
   pfcand_eta = 0;
   pfcand_phi = 0;
   pfcand_theta = 0;
   pfcand_energy = 0;
   pfcand_charge = 0;
   photons_energy = 0;
   photons_et = 0;
   photons_eta = 0;
   photons_phi = 0;
   photons_pt = 0;
   photons_px = 0;
   photons_py = 0;
   photons_pz = 0;
   photons_status = 0;
   photons_theta = 0;
   photons_hadOverEM = 0;
   photons_hadTowOverEM = 0;
   photons_scEnergy = 0;
   photons_scRawEnergy = 0;
   photons_scEta = 0;
   photons_scPhi = 0;
   photons_scEtaWidth = 0;
   photons_scPhiWidth = 0;
   photons_tIso = 0;
   photons_ecalIso = 0;
   photons_hcalIso = 0;
   photons_isoEcalRecHitDR04 = 0;
   photons_isoHcalRecHitDR04 = 0;
   photons_isoSolidTrkConeDR04 = 0;
   photons_isoHollowTrkConeDR04 = 0;
   photons_nTrkSolidConeDR04 = 0;
   photons_nTrkHollowConeDR04 = 0;
   photons_isoEcalRecHitDR03 = 0;
   photons_isoHcalRecHitDR03 = 0;
   photons_isoSolidTrkConeDR03 = 0;
   photons_isoHollowTrkConeDR03 = 0;
   photons_nTrkSolidConeDR03 = 0;
   photons_nTrkHollowConeDR03 = 0;
   photons_isAlsoElectron = 0;
   photons_hasPixelSeed = 0;
   photons_isConverted = 0;
   photons_isEBGap = 0;
   photons_isEEGap = 0;
   photons_isEBEEGap = 0;
   photons_isEBPho = 0;
   photons_isEEPho = 0;
   photons_isLoosePhoton = 0;
   photons_isTightPhoton = 0;
   photons_maxEnergyXtal = 0;
   photons_e1x5 = 0;
   photons_e2x5 = 0;
   photons_e3x3 = 0;
   photons_e5x5 = 0;
   photons_sigmaEtaEta = 0;
   photons_sigmaIetaIeta = 0;
   photons_r9 = 0;
   photons_gen_et = 0;
   photons_gen_eta = 0;
   photons_gen_phi = 0;
   photons_gen_id = 0;
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_xErr = 0;
   pv_yErr = 0;
   pv_zErr = 0;
   pv_chi2 = 0;
   pv_ndof = 0;
   pv_isFake = 0;
   pv_isValid = 0;
   pv_tracksSize = 0;
   taus_status = 0;
   taus_phi = 0;
   taus_pt = 0;
   taus_pz = 0;
   taus_px = 0;
   taus_py = 0;
   taus_eta = 0;
   taus_theta = 0;
   taus_et = 0;
   taus_energy = 0;
   taus_charge = 0;
   taus_emf = 0;
   taus_hcalTotOverPLead = 0;
   taus_hcalMaxOverPLead = 0;
   taus_hcal3x3OverPLead = 0;
   taus_ecalStripSumEOverPLead = 0;
   taus_elecPreIdOutput = 0;
   taus_elecPreIdDecision = 0;
   taus_leadPFChargedHadrCand_pt = 0;
   taus_leadPFChargedHadrCand_charge = 0;
   taus_leadPFChargedHadrCand_eta = 0;
   taus_leadPFChargedHadrCand_ECAL_eta = 0;
   taus_leadPFChargedHadrCand_phi = 0;
   taus_isoPFGammaCandsEtSum = 0;
   taus_isoPFChargedHadrCandsPtSum = 0;
   taus_leadingTrackFinding = 0;
   taus_leadingTrackPtCut = 0;
   taus_trackIsolation = 0;
   taus_ecalIsolation = 0;
   taus_byIsolation = 0;
   taus_againstElectron = 0;
   taus_againstMuon = 0;
   taus_taNC_quarter = 0;
   taus_taNC_one = 0;
   taus_taNC_half = 0;
   taus_taNC_tenth = 0;
   taus_taNC = 0;
   taus_byIsoUsingLeadingPi = 0;
   taus_tkIsoUsingLeadingPi = 0;
   taus_ecalIsoUsingLeadingPi = 0;
   taus_againstElectronLoose = 0;
   taus_againstElectronMedium = 0;
   taus_againstElectronTight = 0;
   taus_againstElectronMVA = 0;
   taus_againstMuonLoose = 0;
   taus_againstMuonMedium = 0;
   taus_againstMuonTight = 0;
   taus_decayModeFinding = 0;
   taus_byVLooseIsolation = 0;
   taus_byLooseIsolation = 0;
   taus_byMediumIsolation = 0;
   taus_byTightIsolation = 0;
   taus_byVLooseIsolationDeltaBetaCorr = 0;
   taus_byLooseIsolationDeltaBetaCorr = 0;
   taus_byMediumIsolationDeltaBetaCorr = 0;
   taus_byTightIsolationDeltaBetaCorr = 0;
   taus_signalPFChargedHadrCandsSize = 0;
   taus_muDecision = 0;
   taus_Nprongs = 0;
   model_params = 0;
   fastjets_AK4_R1p2_R0p5pT10_px = 0;
   fastjets_AK4_R1p2_R0p5pT10_py = 0;
   fastjets_AK4_R1p2_R0p5pT10_pz = 0;
   fastjets_AK4_R1p2_R0p5pT10_energy = 0;
   fastjets_AK4_R1p2_R0p5pT10_phi = 0;
   fastjets_AK4_R1p2_R0p5pT10_eta = 0;
   fastjets_AK4_R1p2_R0p5pT10_index = 0;
   fastjets_AK4_R1p2_R0p5pT10_nconstituents = 0;
   fastjets_AK4_R1p2_R0p5pT15_px = 0;
   fastjets_AK4_R1p2_R0p5pT15_py = 0;
   fastjets_AK4_R1p2_R0p5pT15_pz = 0;
   fastjets_AK4_R1p2_R0p5pT15_energy = 0;
   fastjets_AK4_R1p2_R0p5pT15_phi = 0;
   fastjets_AK4_R1p2_R0p5pT15_eta = 0;
   fastjets_AK4_R1p2_R0p5pT15_index = 0;
   fastjets_AK4_R1p2_R0p5pT15_nconstituents = 0;
   fastjets_AK4_R1p2_R0p5pT20_px = 0;
   fastjets_AK4_R1p2_R0p5pT20_py = 0;
   fastjets_AK4_R1p2_R0p5pT20_pz = 0;
   fastjets_AK4_R1p2_R0p5pT20_energy = 0;
   fastjets_AK4_R1p2_R0p5pT20_phi = 0;
   fastjets_AK4_R1p2_R0p5pT20_eta = 0;
   fastjets_AK4_R1p2_R0p5pT20_index = 0;
   fastjets_AK4_R1p2_R0p5pT20_nconstituents = 0;
   fastjets_AK4_R1p2_R0p5pT25_px = 0;
   fastjets_AK4_R1p2_R0p5pT25_py = 0;
   fastjets_AK4_R1p2_R0p5pT25_pz = 0;
   fastjets_AK4_R1p2_R0p5pT25_energy = 0;
   fastjets_AK4_R1p2_R0p5pT25_phi = 0;
   fastjets_AK4_R1p2_R0p5pT25_eta = 0;
   fastjets_AK4_R1p2_R0p5pT25_index = 0;
   fastjets_AK4_R1p2_R0p5pT25_nconstituents = 0;
   fastjets_AK4_R1p2_R0p5pT30_px = 0;
   fastjets_AK4_R1p2_R0p5pT30_py = 0;
   fastjets_AK4_R1p2_R0p5pT30_pz = 0;
   fastjets_AK4_R1p2_R0p5pT30_energy = 0;
   fastjets_AK4_R1p2_R0p5pT30_phi = 0;
   fastjets_AK4_R1p2_R0p5pT30_eta = 0;
   fastjets_AK4_R1p2_R0p5pT30_index = 0;
   fastjets_AK4_R1p2_R0p5pT30_nconstituents = 0;

   fChain->SetBranchAddress("NbeamSpot", &NbeamSpot, &b_NbeamSpot);
   fChain->SetBranchAddress("beamSpot_x", &beamSpot_x, &b_beamSpot_x);
   fChain->SetBranchAddress("beamSpot_y", &beamSpot_y, &b_beamSpot_y);
   fChain->SetBranchAddress("beamSpot_z", &beamSpot_z, &b_beamSpot_z);
   fChain->SetBranchAddress("beamSpot_x0Error", &beamSpot_x0Error, &b_beamSpot_x0Error);
   fChain->SetBranchAddress("beamSpot_y0Error", &beamSpot_y0Error, &b_beamSpot_y0Error);
   fChain->SetBranchAddress("beamSpot_z0Error", &beamSpot_z0Error, &b_beamSpot_z0Error);
   fChain->SetBranchAddress("beamSpot_sigmaZ", &beamSpot_sigmaZ, &b_beamSpot_sigmaZ);
   fChain->SetBranchAddress("beamSpot_sigmaZ0Error", &beamSpot_sigmaZ0Error, &b_beamSpot_sigmaZ0Error);
   fChain->SetBranchAddress("beamSpot_dxdz", &beamSpot_dxdz, &b_beamSpot_dxdz);
   fChain->SetBranchAddress("beamSpot_dxdzError", &beamSpot_dxdzError, &b_beamSpot_dxdzError);
   fChain->SetBranchAddress("beamSpot_dydz", &beamSpot_dydz, &b_beamSpot_dydz);
   fChain->SetBranchAddress("beamSpot_dydzError", &beamSpot_dydzError, &b_beamSpot_dydzError);
   fChain->SetBranchAddress("beamSpot_beamWidthX", &beamSpot_beamWidthX, &b_beamSpot_beamWidthX);
   fChain->SetBranchAddress("beamSpot_beamWidthY", &beamSpot_beamWidthY, &b_beamSpot_beamWidthY);
   fChain->SetBranchAddress("beamSpot_beamWidthXError", &beamSpot_beamWidthXError, &b_beamSpot_beamWidthXError);
   fChain->SetBranchAddress("beamSpot_beamWidthYError", &beamSpot_beamWidthYError, &b_beamSpot_beamWidthYError);
   fChain->SetBranchAddress("Nels", &Nels, &b_Nels);
   fChain->SetBranchAddress("els_energy", &els_energy, &b_els_energy);
   fChain->SetBranchAddress("els_et", &els_et, &b_els_et);
   fChain->SetBranchAddress("els_eta", &els_eta, &b_els_eta);
   fChain->SetBranchAddress("els_phi", &els_phi, &b_els_phi);
   fChain->SetBranchAddress("els_pt", &els_pt, &b_els_pt);
   fChain->SetBranchAddress("els_px", &els_px, &b_els_px);
   fChain->SetBranchAddress("els_py", &els_py, &b_els_py);
   fChain->SetBranchAddress("els_pz", &els_pz, &b_els_pz);
   fChain->SetBranchAddress("els_status", &els_status, &b_els_status);
   fChain->SetBranchAddress("els_theta", &els_theta, &b_els_theta);
   fChain->SetBranchAddress("els_pfIsolationR03_sumChargedHadronPt", &els_pfIsolationR03_sumChargedHadronPt, &b_els_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("els_pfIsolationR03_sumNeutralHadronEt", &els_pfIsolationR03_sumNeutralHadronEt, &b_els_pfIsolationR03_sumNeutralHadronEt);
   fChain->SetBranchAddress("els_pfIsolationR03_sumPhotonEt", &els_pfIsolationR03_sumPhotonEt, &b_els_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("els_pfIsolationR03_sumPUPt", &els_pfIsolationR03_sumPUPt, &b_els_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("els_full5x5_sigmaIetaIeta", &els_full5x5_sigmaIetaIeta, &b_els_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("els_gen_id", &els_gen_id, &b_els_gen_id);
   fChain->SetBranchAddress("els_gen_phi", &els_gen_phi, &b_els_gen_phi);
   fChain->SetBranchAddress("els_gen_pt", &els_gen_pt, &b_els_gen_pt);
   fChain->SetBranchAddress("els_gen_pz", &els_gen_pz, &b_els_gen_pz);
   fChain->SetBranchAddress("els_gen_px", &els_gen_px, &b_els_gen_px);
   fChain->SetBranchAddress("els_gen_py", &els_gen_py, &b_els_gen_py);
   fChain->SetBranchAddress("els_gen_eta", &els_gen_eta, &b_els_gen_eta);
   fChain->SetBranchAddress("els_gen_theta", &els_gen_theta, &b_els_gen_theta);
   fChain->SetBranchAddress("els_gen_et", &els_gen_et, &b_els_gen_et);
   fChain->SetBranchAddress("els_gen_mother_id", &els_gen_mother_id, &b_els_gen_mother_id);
   fChain->SetBranchAddress("els_gen_mother_phi", &els_gen_mother_phi, &b_els_gen_mother_phi);
   fChain->SetBranchAddress("els_gen_mother_pt", &els_gen_mother_pt, &b_els_gen_mother_pt);
   fChain->SetBranchAddress("els_gen_mother_pz", &els_gen_mother_pz, &b_els_gen_mother_pz);
   fChain->SetBranchAddress("els_gen_mother_px", &els_gen_mother_px, &b_els_gen_mother_px);
   fChain->SetBranchAddress("els_gen_mother_py", &els_gen_mother_py, &b_els_gen_mother_py);
   fChain->SetBranchAddress("els_gen_mother_eta", &els_gen_mother_eta, &b_els_gen_mother_eta);
   fChain->SetBranchAddress("els_gen_mother_theta", &els_gen_mother_theta, &b_els_gen_mother_theta);
   fChain->SetBranchAddress("els_gen_mother_et", &els_gen_mother_et, &b_els_gen_mother_et);
   fChain->SetBranchAddress("els_tightId", &els_tightId, &b_els_tightId);
   fChain->SetBranchAddress("els_looseId", &els_looseId, &b_els_looseId);
   fChain->SetBranchAddress("els_robustTightId", &els_robustTightId, &b_els_robustTightId);
   fChain->SetBranchAddress("els_robustLooseId", &els_robustLooseId, &b_els_robustLooseId);
   fChain->SetBranchAddress("els_robustHighEnergyId", &els_robustHighEnergyId, &b_els_robustHighEnergyId);
   fChain->SetBranchAddress("els_simpleEleId95relIso", &els_simpleEleId95relIso, &b_els_simpleEleId95relIso);
   fChain->SetBranchAddress("els_simpleEleId90relIso", &els_simpleEleId90relIso, &b_els_simpleEleId90relIso);
   fChain->SetBranchAddress("els_simpleEleId85relIso", &els_simpleEleId85relIso, &b_els_simpleEleId85relIso);
   fChain->SetBranchAddress("els_simpleEleId80relIso", &els_simpleEleId80relIso, &b_els_simpleEleId80relIso);
   fChain->SetBranchAddress("els_simpleEleId70relIso", &els_simpleEleId70relIso, &b_els_simpleEleId70relIso);
   fChain->SetBranchAddress("els_simpleEleId60relIso", &els_simpleEleId60relIso, &b_els_simpleEleId60relIso);
   fChain->SetBranchAddress("els_simpleEleId95cIso", &els_simpleEleId95cIso, &b_els_simpleEleId95cIso);
   fChain->SetBranchAddress("els_simpleEleId90cIso", &els_simpleEleId90cIso, &b_els_simpleEleId90cIso);
   fChain->SetBranchAddress("els_simpleEleId85cIso", &els_simpleEleId85cIso, &b_els_simpleEleId85cIso);
   fChain->SetBranchAddress("els_simpleEleId80cIso", &els_simpleEleId80cIso, &b_els_simpleEleId80cIso);
   fChain->SetBranchAddress("els_simpleEleId70cIso", &els_simpleEleId70cIso, &b_els_simpleEleId70cIso);
   fChain->SetBranchAddress("els_simpleEleId60cIso", &els_simpleEleId60cIso, &b_els_simpleEleId60cIso);
   fChain->SetBranchAddress("els_cIso", &els_cIso, &b_els_cIso);
   fChain->SetBranchAddress("els_tIso", &els_tIso, &b_els_tIso);
   fChain->SetBranchAddress("els_ecalIso", &els_ecalIso, &b_els_ecalIso);
   fChain->SetBranchAddress("els_hcalIso", &els_hcalIso, &b_els_hcalIso);
   fChain->SetBranchAddress("els_chi2", &els_chi2, &b_els_chi2);
   fChain->SetBranchAddress("els_charge", &els_charge, &b_els_charge);
   fChain->SetBranchAddress("els_caloEnergy", &els_caloEnergy, &b_els_caloEnergy);
   fChain->SetBranchAddress("els_hadOverEm", &els_hadOverEm, &b_els_hadOverEm);
   fChain->SetBranchAddress("els_hcalOverEcalBc", &els_hcalOverEcalBc, &b_els_hcalOverEcalBc);
   fChain->SetBranchAddress("els_eOverPIn", &els_eOverPIn, &b_els_eOverPIn);
   fChain->SetBranchAddress("els_eSeedOverPOut", &els_eSeedOverPOut, &b_els_eSeedOverPOut);
   fChain->SetBranchAddress("els_sigmaEtaEta", &els_sigmaEtaEta, &b_els_sigmaEtaEta);
   fChain->SetBranchAddress("els_sigmaIEtaIEta", &els_sigmaIEtaIEta, &b_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("els_scEnergy", &els_scEnergy, &b_els_scEnergy);
   fChain->SetBranchAddress("els_scRawEnergy", &els_scRawEnergy, &b_els_scRawEnergy);
   fChain->SetBranchAddress("els_scSeedEnergy", &els_scSeedEnergy, &b_els_scSeedEnergy);
   fChain->SetBranchAddress("els_scEta", &els_scEta, &b_els_scEta);
   fChain->SetBranchAddress("els_scPhi", &els_scPhi, &b_els_scPhi);
   fChain->SetBranchAddress("els_scEtaWidth", &els_scEtaWidth, &b_els_scEtaWidth);
   fChain->SetBranchAddress("els_scPhiWidth", &els_scPhiWidth, &b_els_scPhiWidth);
   fChain->SetBranchAddress("els_scE1x5", &els_scE1x5, &b_els_scE1x5);
   fChain->SetBranchAddress("els_scE2x5Max", &els_scE2x5Max, &b_els_scE2x5Max);
   fChain->SetBranchAddress("els_scE5x5", &els_scE5x5, &b_els_scE5x5);
   fChain->SetBranchAddress("els_isEB", &els_isEB, &b_els_isEB);
   fChain->SetBranchAddress("els_isEE", &els_isEE, &b_els_isEE);
   fChain->SetBranchAddress("els_dEtaIn", &els_dEtaIn, &b_els_dEtaIn);
   fChain->SetBranchAddress("els_dPhiIn", &els_dPhiIn, &b_els_dPhiIn);
   fChain->SetBranchAddress("els_dEtaOut", &els_dEtaOut, &b_els_dEtaOut);
   fChain->SetBranchAddress("els_dPhiOut", &els_dPhiOut, &b_els_dPhiOut);
   fChain->SetBranchAddress("els_numvalhits", &els_numvalhits, &b_els_numvalhits);
   fChain->SetBranchAddress("els_numlosthits", &els_numlosthits, &b_els_numlosthits);
   fChain->SetBranchAddress("els_basicClustersSize", &els_basicClustersSize, &b_els_basicClustersSize);
   fChain->SetBranchAddress("els_tk_pz", &els_tk_pz, &b_els_tk_pz);
   fChain->SetBranchAddress("els_tk_pt", &els_tk_pt, &b_els_tk_pt);
   fChain->SetBranchAddress("els_tk_phi", &els_tk_phi, &b_els_tk_phi);
   fChain->SetBranchAddress("els_tk_eta", &els_tk_eta, &b_els_tk_eta);
   fChain->SetBranchAddress("els_d0dum", &els_d0dum, &b_els_d0dum);
   fChain->SetBranchAddress("els_dz", &els_dz, &b_els_dz);
   fChain->SetBranchAddress("els_vx", &els_vx, &b_els_vx);
   fChain->SetBranchAddress("els_vy", &els_vy, &b_els_vy);
   fChain->SetBranchAddress("els_vz", &els_vz, &b_els_vz);
   fChain->SetBranchAddress("els_ndof", &els_ndof, &b_els_ndof);
   fChain->SetBranchAddress("els_ptError", &els_ptError, &b_els_ptError);
   fChain->SetBranchAddress("els_d0dumError", &els_d0dumError, &b_els_d0dumError);
   fChain->SetBranchAddress("els_dzError", &els_dzError, &b_els_dzError);
   fChain->SetBranchAddress("els_etaError", &els_etaError, &b_els_etaError);
   fChain->SetBranchAddress("els_phiError", &els_phiError, &b_els_phiError);
   fChain->SetBranchAddress("els_tk_charge", &els_tk_charge, &b_els_tk_charge);
   fChain->SetBranchAddress("els_core_ecalDrivenSeed", &els_core_ecalDrivenSeed, &b_els_core_ecalDrivenSeed);
   fChain->SetBranchAddress("els_n_inner_layer", &els_n_inner_layer, &b_els_n_inner_layer);
   fChain->SetBranchAddress("els_n_outer_layer", &els_n_outer_layer, &b_els_n_outer_layer);
   fChain->SetBranchAddress("els_ctf_tk_id", &els_ctf_tk_id, &b_els_ctf_tk_id);
   fChain->SetBranchAddress("els_ctf_tk_charge", &els_ctf_tk_charge, &b_els_ctf_tk_charge);
   fChain->SetBranchAddress("els_ctf_tk_eta", &els_ctf_tk_eta, &b_els_ctf_tk_eta);
   fChain->SetBranchAddress("els_ctf_tk_phi", &els_ctf_tk_phi, &b_els_ctf_tk_phi);
   fChain->SetBranchAddress("els_fbrem", &els_fbrem, &b_els_fbrem);
   fChain->SetBranchAddress("els_shFracInnerHits", &els_shFracInnerHits, &b_els_shFracInnerHits);
   fChain->SetBranchAddress("els_dr03EcalRecHitSumEt", &els_dr03EcalRecHitSumEt, &b_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr03HcalTowerSumEt", &els_dr03HcalTowerSumEt, &b_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth1TowerSumEt", &els_dr03HcalDepth1TowerSumEt, &b_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth2TowerSumEt", &els_dr03HcalDepth2TowerSumEt, &b_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr03TkSumPt", &els_dr03TkSumPt, &b_els_dr03TkSumPt);
   fChain->SetBranchAddress("els_dr04EcalRecHitSumEt", &els_dr04EcalRecHitSumEt, &b_els_dr04EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr04HcalTowerSumEt", &els_dr04HcalTowerSumEt, &b_els_dr04HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr04HcalDepth1TowerSumEt", &els_dr04HcalDepth1TowerSumEt, &b_els_dr04HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr04HcalDepth2TowerSumEt", &els_dr04HcalDepth2TowerSumEt, &b_els_dr04HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr04TkSumPt", &els_dr04TkSumPt, &b_els_dr04TkSumPt);
   fChain->SetBranchAddress("els_cpx", &els_cpx, &b_els_cpx);
   fChain->SetBranchAddress("els_cpy", &els_cpy, &b_els_cpy);
   fChain->SetBranchAddress("els_cpz", &els_cpz, &b_els_cpz);
   fChain->SetBranchAddress("els_vpx", &els_vpx, &b_els_vpx);
   fChain->SetBranchAddress("els_vpy", &els_vpy, &b_els_vpy);
   fChain->SetBranchAddress("els_vpz", &els_vpz, &b_els_vpz);
   fChain->SetBranchAddress("els_cx", &els_cx, &b_els_cx);
   fChain->SetBranchAddress("els_cy", &els_cy, &b_els_cy);
   fChain->SetBranchAddress("els_cz", &els_cz, &b_els_cz);
   fChain->SetBranchAddress("els_PATpassConversionVeto", &els_PATpassConversionVeto, &b_els_PATpassConversionVeto);
   fChain->SetBranchAddress("Njets_AK4", &Njets_AK4, &b_Njets_AK4);
   fChain->SetBranchAddress("jets_AK4_status", &jets_AK4_status, &b_jets_AK4_status);
   fChain->SetBranchAddress("jets_AK4_phi", &jets_AK4_phi, &b_jets_AK4_phi);
   fChain->SetBranchAddress("jets_AK4_pt", &jets_AK4_pt, &b_jets_AK4_pt);
   fChain->SetBranchAddress("jets_AK4_pz", &jets_AK4_pz, &b_jets_AK4_pz);
   fChain->SetBranchAddress("jets_AK4_px", &jets_AK4_px, &b_jets_AK4_px);
   fChain->SetBranchAddress("jets_AK4_py", &jets_AK4_py, &b_jets_AK4_py);
   fChain->SetBranchAddress("jets_AK4_eta", &jets_AK4_eta, &b_jets_AK4_eta);
   fChain->SetBranchAddress("jets_AK4_theta", &jets_AK4_theta, &b_jets_AK4_theta);
   fChain->SetBranchAddress("jets_AK4_et", &jets_AK4_et, &b_jets_AK4_et);
   fChain->SetBranchAddress("jets_AK4_energy", &jets_AK4_energy, &b_jets_AK4_energy);
   fChain->SetBranchAddress("jets_AK4_parton_Id", &jets_AK4_parton_Id, &b_jets_AK4_parton_Id);
   fChain->SetBranchAddress("jets_AK4_parton_motherId", &jets_AK4_parton_motherId, &b_jets_AK4_parton_motherId);
   fChain->SetBranchAddress("jets_AK4_parton_pt", &jets_AK4_parton_pt, &b_jets_AK4_parton_pt);
   fChain->SetBranchAddress("jets_AK4_parton_phi", &jets_AK4_parton_phi, &b_jets_AK4_parton_phi);
   fChain->SetBranchAddress("jets_AK4_parton_eta", &jets_AK4_parton_eta, &b_jets_AK4_parton_eta);
   fChain->SetBranchAddress("jets_AK4_parton_Energy", &jets_AK4_parton_Energy, &b_jets_AK4_parton_Energy);
   fChain->SetBranchAddress("jets_AK4_parton_mass", &jets_AK4_parton_mass, &b_jets_AK4_parton_mass);
   fChain->SetBranchAddress("jets_AK4_gen_et", &jets_AK4_gen_et, &b_jets_AK4_gen_et);
   fChain->SetBranchAddress("jets_AK4_gen_pt", &jets_AK4_gen_pt, &b_jets_AK4_gen_pt);
   fChain->SetBranchAddress("jets_AK4_gen_eta", &jets_AK4_gen_eta, &b_jets_AK4_gen_eta);
   fChain->SetBranchAddress("jets_AK4_gen_phi", &jets_AK4_gen_phi, &b_jets_AK4_gen_phi);
   fChain->SetBranchAddress("jets_AK4_gen_mass", &jets_AK4_gen_mass, &b_jets_AK4_gen_mass);
   fChain->SetBranchAddress("jets_AK4_gen_Energy", &jets_AK4_gen_Energy, &b_jets_AK4_gen_Energy);
   fChain->SetBranchAddress("jets_AK4_gen_Id", &jets_AK4_gen_Id, &b_jets_AK4_gen_Id);
   fChain->SetBranchAddress("jets_AK4_gen_motherID", &jets_AK4_gen_motherID, &b_jets_AK4_gen_motherID);
   fChain->SetBranchAddress("jets_AK4_gen_threeCharge", &jets_AK4_gen_threeCharge, &b_jets_AK4_gen_threeCharge);
   fChain->SetBranchAddress("jets_AK4_partonFlavour", &jets_AK4_partonFlavour, &b_jets_AK4_partonFlavour);
   fChain->SetBranchAddress("jets_AK4_btag_TC_highPur", &jets_AK4_btag_TC_highPur, &b_jets_AK4_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK4_btag_TC_highEff", &jets_AK4_btag_TC_highEff, &b_jets_AK4_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK4_btag_jetProb", &jets_AK4_btag_jetProb, &b_jets_AK4_btag_jetProb);
   fChain->SetBranchAddress("jets_AK4_btag_jetBProb", &jets_AK4_btag_jetBProb, &b_jets_AK4_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK4_btag_softEle", &jets_AK4_btag_softEle, &b_jets_AK4_btag_softEle);
   fChain->SetBranchAddress("jets_AK4_btag_softMuon", &jets_AK4_btag_softMuon, &b_jets_AK4_btag_softMuon);
   fChain->SetBranchAddress("jets_AK4_btag_secVertexHighPur", &jets_AK4_btag_secVertexHighPur, &b_jets_AK4_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK4_btag_secVertexHighEff", &jets_AK4_btag_secVertexHighEff, &b_jets_AK4_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK4_btag_secVertexCombined", &jets_AK4_btag_secVertexCombined, &b_jets_AK4_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK4_jetCharge", &jets_AK4_jetCharge, &b_jets_AK4_jetCharge);
   fChain->SetBranchAddress("jets_AK4_chgEmE", &jets_AK4_chgEmE, &b_jets_AK4_chgEmE);
   fChain->SetBranchAddress("jets_AK4_chgHadE", &jets_AK4_chgHadE, &b_jets_AK4_chgHadE);
   fChain->SetBranchAddress("jets_AK4_photonEnergy", &jets_AK4_photonEnergy, &b_jets_AK4_photonEnergy);
   fChain->SetBranchAddress("jets_AK4_chgMuE", &jets_AK4_chgMuE, &b_jets_AK4_chgMuE);
   fChain->SetBranchAddress("jets_AK4_chg_Mult", &jets_AK4_chg_Mult, &b_jets_AK4_chg_Mult);
   fChain->SetBranchAddress("jets_AK4_neutralEmE", &jets_AK4_neutralEmE, &b_jets_AK4_neutralEmE);
   fChain->SetBranchAddress("jets_AK4_neutralHadE", &jets_AK4_neutralHadE, &b_jets_AK4_neutralHadE);
   fChain->SetBranchAddress("jets_AK4_neutral_Mult", &jets_AK4_neutral_Mult, &b_jets_AK4_neutral_Mult);
   fChain->SetBranchAddress("jets_AK4_mu_Mult", &jets_AK4_mu_Mult, &b_jets_AK4_mu_Mult);
   fChain->SetBranchAddress("jets_AK4_emf", &jets_AK4_emf, &b_jets_AK4_emf);
   fChain->SetBranchAddress("jets_AK4_ehf", &jets_AK4_ehf, &b_jets_AK4_ehf);
   fChain->SetBranchAddress("jets_AK4_n60", &jets_AK4_n60, &b_jets_AK4_n60);
   fChain->SetBranchAddress("jets_AK4_n90", &jets_AK4_n90, &b_jets_AK4_n90);
   fChain->SetBranchAddress("jets_AK4_etaetaMoment", &jets_AK4_etaetaMoment, &b_jets_AK4_etaetaMoment);
   fChain->SetBranchAddress("jets_AK4_etaphiMoment", &jets_AK4_etaphiMoment, &b_jets_AK4_etaphiMoment);
   fChain->SetBranchAddress("jets_AK4_phiphiMoment", &jets_AK4_phiphiMoment, &b_jets_AK4_phiphiMoment);
   fChain->SetBranchAddress("jets_AK4_n90Hits", &jets_AK4_n90Hits, &b_jets_AK4_n90Hits);
   fChain->SetBranchAddress("jets_AK4_fHPD", &jets_AK4_fHPD, &b_jets_AK4_fHPD);
   fChain->SetBranchAddress("jets_AK4_fRBX", &jets_AK4_fRBX, &b_jets_AK4_fRBX);
   fChain->SetBranchAddress("jets_AK4_hitsInN90", &jets_AK4_hitsInN90, &b_jets_AK4_hitsInN90);
   fChain->SetBranchAddress("jets_AK4_nECALTowers", &jets_AK4_nECALTowers, &b_jets_AK4_nECALTowers);
   fChain->SetBranchAddress("jets_AK4_nHCALTowers", &jets_AK4_nHCALTowers, &b_jets_AK4_nHCALTowers);
   fChain->SetBranchAddress("jets_AK4_fSubDetector1", &jets_AK4_fSubDetector1, &b_jets_AK4_fSubDetector1);
   fChain->SetBranchAddress("jets_AK4_fSubDetector2", &jets_AK4_fSubDetector2, &b_jets_AK4_fSubDetector2);
   fChain->SetBranchAddress("jets_AK4_fSubDetector3", &jets_AK4_fSubDetector3, &b_jets_AK4_fSubDetector3);
   fChain->SetBranchAddress("jets_AK4_fSubDetector4", &jets_AK4_fSubDetector4, &b_jets_AK4_fSubDetector4);
   fChain->SetBranchAddress("jets_AK4_area", &jets_AK4_area, &b_jets_AK4_area);
   fChain->SetBranchAddress("jets_AK4_corrFactorRaw", &jets_AK4_corrFactorRaw, &b_jets_AK4_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK4_rawPt", &jets_AK4_rawPt, &b_jets_AK4_rawPt);
   fChain->SetBranchAddress("jets_AK4_mass", &jets_AK4_mass, &b_jets_AK4_mass);
   fChain->SetBranchAddress("Nmc_doc", &Nmc_doc, &b_Nmc_doc);
   fChain->SetBranchAddress("mc_doc_id", &mc_doc_id, &b_mc_doc_id);
   fChain->SetBranchAddress("mc_doc_pt", &mc_doc_pt, &b_mc_doc_pt);
   fChain->SetBranchAddress("mc_doc_px", &mc_doc_px, &b_mc_doc_px);
   fChain->SetBranchAddress("mc_doc_py", &mc_doc_py, &b_mc_doc_py);
   fChain->SetBranchAddress("mc_doc_pz", &mc_doc_pz, &b_mc_doc_pz);
   fChain->SetBranchAddress("mc_doc_eta", &mc_doc_eta, &b_mc_doc_eta);
   fChain->SetBranchAddress("mc_doc_phi", &mc_doc_phi, &b_mc_doc_phi);
   fChain->SetBranchAddress("mc_doc_theta", &mc_doc_theta, &b_mc_doc_theta);
   fChain->SetBranchAddress("mc_doc_energy", &mc_doc_energy, &b_mc_doc_energy);
   fChain->SetBranchAddress("mc_doc_status", &mc_doc_status, &b_mc_doc_status);
   fChain->SetBranchAddress("mc_doc_charge", &mc_doc_charge, &b_mc_doc_charge);
   fChain->SetBranchAddress("mc_doc_mother_id", &mc_doc_mother_id, &b_mc_doc_mother_id);
   fChain->SetBranchAddress("mc_doc_grandmother_id", &mc_doc_grandmother_id, &b_mc_doc_grandmother_id);
   fChain->SetBranchAddress("mc_doc_ggrandmother_id", &mc_doc_ggrandmother_id, &b_mc_doc_ggrandmother_id);
   fChain->SetBranchAddress("mc_doc_mother_pt", &mc_doc_mother_pt, &b_mc_doc_mother_pt);
   fChain->SetBranchAddress("mc_doc_vertex_x", &mc_doc_vertex_x, &b_mc_doc_vertex_x);
   fChain->SetBranchAddress("mc_doc_vertex_y", &mc_doc_vertex_y, &b_mc_doc_vertex_y);
   fChain->SetBranchAddress("mc_doc_vertex_z", &mc_doc_vertex_z, &b_mc_doc_vertex_z);
   fChain->SetBranchAddress("mc_doc_mass", &mc_doc_mass, &b_mc_doc_mass);
   fChain->SetBranchAddress("mc_doc_numOfDaughters", &mc_doc_numOfDaughters, &b_mc_doc_numOfDaughters);
   fChain->SetBranchAddress("mc_doc_numOfMothers", &mc_doc_numOfMothers, &b_mc_doc_numOfMothers);
   fChain->SetBranchAddress("Nmc_electrons", &Nmc_electrons, &b_Nmc_electrons);
   fChain->SetBranchAddress("mc_electrons_id", &mc_electrons_id, &b_mc_electrons_id);
   fChain->SetBranchAddress("mc_electrons_pt", &mc_electrons_pt, &b_mc_electrons_pt);
   fChain->SetBranchAddress("mc_electrons_px", &mc_electrons_px, &b_mc_electrons_px);
   fChain->SetBranchAddress("mc_electrons_py", &mc_electrons_py, &b_mc_electrons_py);
   fChain->SetBranchAddress("mc_electrons_pz", &mc_electrons_pz, &b_mc_electrons_pz);
   fChain->SetBranchAddress("mc_electrons_eta", &mc_electrons_eta, &b_mc_electrons_eta);
   fChain->SetBranchAddress("mc_electrons_phi", &mc_electrons_phi, &b_mc_electrons_phi);
   fChain->SetBranchAddress("mc_electrons_theta", &mc_electrons_theta, &b_mc_electrons_theta);
   fChain->SetBranchAddress("mc_electrons_status", &mc_electrons_status, &b_mc_electrons_status);
   fChain->SetBranchAddress("mc_electrons_energy", &mc_electrons_energy, &b_mc_electrons_energy);
   fChain->SetBranchAddress("mc_electrons_charge", &mc_electrons_charge, &b_mc_electrons_charge);
   fChain->SetBranchAddress("mc_electrons_mother_id", &mc_electrons_mother_id, &b_mc_electrons_mother_id);
   fChain->SetBranchAddress("mc_electrons_mother_pt", &mc_electrons_mother_pt, &b_mc_electrons_mother_pt);
   fChain->SetBranchAddress("mc_electrons_grandmother_id", &mc_electrons_grandmother_id, &b_mc_electrons_grandmother_id);
   fChain->SetBranchAddress("mc_electrons_ggrandmother_id", &mc_electrons_ggrandmother_id, &b_mc_electrons_ggrandmother_id);
   fChain->SetBranchAddress("mc_electrons_vertex_x", &mc_electrons_vertex_x, &b_mc_electrons_vertex_x);
   fChain->SetBranchAddress("mc_electrons_vertex_y", &mc_electrons_vertex_y, &b_mc_electrons_vertex_y);
   fChain->SetBranchAddress("mc_electrons_vertex_z", &mc_electrons_vertex_z, &b_mc_electrons_vertex_z);
   fChain->SetBranchAddress("mc_electrons_mass", &mc_electrons_mass, &b_mc_electrons_mass);
   fChain->SetBranchAddress("mc_electrons_numOfDaughters", &mc_electrons_numOfDaughters, &b_mc_electrons_numOfDaughters);
   fChain->SetBranchAddress("Nmc_final", &Nmc_final, &b_Nmc_final);
   fChain->SetBranchAddress("mc_final_id", &mc_final_id, &b_mc_final_id);
   fChain->SetBranchAddress("mc_final_pt", &mc_final_pt, &b_mc_final_pt);
   fChain->SetBranchAddress("mc_final_px", &mc_final_px, &b_mc_final_px);
   fChain->SetBranchAddress("mc_final_py", &mc_final_py, &b_mc_final_py);
   fChain->SetBranchAddress("mc_final_pz", &mc_final_pz, &b_mc_final_pz);
   fChain->SetBranchAddress("mc_final_eta", &mc_final_eta, &b_mc_final_eta);
   fChain->SetBranchAddress("mc_final_phi", &mc_final_phi, &b_mc_final_phi);
   fChain->SetBranchAddress("mc_final_theta", &mc_final_theta, &b_mc_final_theta);
   fChain->SetBranchAddress("mc_final_energy", &mc_final_energy, &b_mc_final_energy);
   fChain->SetBranchAddress("mc_final_status", &mc_final_status, &b_mc_final_status);
   fChain->SetBranchAddress("mc_final_charge", &mc_final_charge, &b_mc_final_charge);
   fChain->SetBranchAddress("mc_final_mother_id", &mc_final_mother_id, &b_mc_final_mother_id);
   fChain->SetBranchAddress("mc_final_grandmother_id", &mc_final_grandmother_id, &b_mc_final_grandmother_id);
   fChain->SetBranchAddress("mc_final_ggrandmother_id", &mc_final_ggrandmother_id, &b_mc_final_ggrandmother_id);
   fChain->SetBranchAddress("mc_final_mother_pt", &mc_final_mother_pt, &b_mc_final_mother_pt);
   fChain->SetBranchAddress("mc_final_vertex_x", &mc_final_vertex_x, &b_mc_final_vertex_x);
   fChain->SetBranchAddress("mc_final_vertex_y", &mc_final_vertex_y, &b_mc_final_vertex_y);
   fChain->SetBranchAddress("mc_final_vertex_z", &mc_final_vertex_z, &b_mc_final_vertex_z);
   fChain->SetBranchAddress("mc_final_mass", &mc_final_mass, &b_mc_final_mass);
   fChain->SetBranchAddress("mc_final_numOfDaughters", &mc_final_numOfDaughters, &b_mc_final_numOfDaughters);
   fChain->SetBranchAddress("mc_final_numOfMothers", &mc_final_numOfMothers, &b_mc_final_numOfMothers);
   fChain->SetBranchAddress("Nmc_jets", &Nmc_jets, &b_Nmc_jets);
   fChain->SetBranchAddress("mc_jets_phi", &mc_jets_phi, &b_mc_jets_phi);
   fChain->SetBranchAddress("mc_jets_pt", &mc_jets_pt, &b_mc_jets_pt);
   fChain->SetBranchAddress("mc_jets_pz", &mc_jets_pz, &b_mc_jets_pz);
   fChain->SetBranchAddress("mc_jets_px", &mc_jets_px, &b_mc_jets_px);
   fChain->SetBranchAddress("mc_jets_py", &mc_jets_py, &b_mc_jets_py);
   fChain->SetBranchAddress("mc_jets_eta", &mc_jets_eta, &b_mc_jets_eta);
   fChain->SetBranchAddress("mc_jets_theta", &mc_jets_theta, &b_mc_jets_theta);
   fChain->SetBranchAddress("mc_jets_et", &mc_jets_et, &b_mc_jets_et);
   fChain->SetBranchAddress("mc_jets_energy", &mc_jets_energy, &b_mc_jets_energy);
   fChain->SetBranchAddress("mc_jets_emEnergy", &mc_jets_emEnergy, &b_mc_jets_emEnergy);
   fChain->SetBranchAddress("mc_jets_hadEnergy", &mc_jets_hadEnergy, &b_mc_jets_hadEnergy);
   fChain->SetBranchAddress("mc_jets_invisibleEnergy", &mc_jets_invisibleEnergy, &b_mc_jets_invisibleEnergy);
   fChain->SetBranchAddress("mc_jets_auxiliaryEnergy", &mc_jets_auxiliaryEnergy, &b_mc_jets_auxiliaryEnergy);
   fChain->SetBranchAddress("mc_jets_etaetaMoment", &mc_jets_etaetaMoment, &b_mc_jets_etaetaMoment);
   fChain->SetBranchAddress("mc_jets_etaphiMoment", &mc_jets_etaphiMoment, &b_mc_jets_etaphiMoment);
   fChain->SetBranchAddress("mc_jets_phiphiMoment", &mc_jets_phiphiMoment, &b_mc_jets_phiphiMoment);
   fChain->SetBranchAddress("mc_jets_mass", &mc_jets_mass, &b_mc_jets_mass);
   fChain->SetBranchAddress("Nmc_mus", &Nmc_mus, &b_Nmc_mus);
   fChain->SetBranchAddress("mc_mus_id", &mc_mus_id, &b_mc_mus_id);
   fChain->SetBranchAddress("mc_mus_pt", &mc_mus_pt, &b_mc_mus_pt);
   fChain->SetBranchAddress("mc_mus_px", &mc_mus_px, &b_mc_mus_px);
   fChain->SetBranchAddress("mc_mus_py", &mc_mus_py, &b_mc_mus_py);
   fChain->SetBranchAddress("mc_mus_pz", &mc_mus_pz, &b_mc_mus_pz);
   fChain->SetBranchAddress("mc_mus_eta", &mc_mus_eta, &b_mc_mus_eta);
   fChain->SetBranchAddress("mc_mus_phi", &mc_mus_phi, &b_mc_mus_phi);
   fChain->SetBranchAddress("mc_mus_theta", &mc_mus_theta, &b_mc_mus_theta);
   fChain->SetBranchAddress("mc_mus_status", &mc_mus_status, &b_mc_mus_status);
   fChain->SetBranchAddress("mc_mus_energy", &mc_mus_energy, &b_mc_mus_energy);
   fChain->SetBranchAddress("mc_mus_charge", &mc_mus_charge, &b_mc_mus_charge);
   fChain->SetBranchAddress("mc_mus_mother_id", &mc_mus_mother_id, &b_mc_mus_mother_id);
   fChain->SetBranchAddress("mc_mus_mother_pt", &mc_mus_mother_pt, &b_mc_mus_mother_pt);
   fChain->SetBranchAddress("mc_mus_grandmother_id", &mc_mus_grandmother_id, &b_mc_mus_grandmother_id);
   fChain->SetBranchAddress("mc_mus_ggrandmother_id", &mc_mus_ggrandmother_id, &b_mc_mus_ggrandmother_id);
   fChain->SetBranchAddress("mc_mus_vertex_x", &mc_mus_vertex_x, &b_mc_mus_vertex_x);
   fChain->SetBranchAddress("mc_mus_vertex_y", &mc_mus_vertex_y, &b_mc_mus_vertex_y);
   fChain->SetBranchAddress("mc_mus_vertex_z", &mc_mus_vertex_z, &b_mc_mus_vertex_z);
   fChain->SetBranchAddress("mc_mus_mass", &mc_mus_mass, &b_mc_mus_mass);
   fChain->SetBranchAddress("mc_mus_numOfDaughters", &mc_mus_numOfDaughters, &b_mc_mus_numOfDaughters);
   fChain->SetBranchAddress("Nmc_nues", &Nmc_nues, &b_Nmc_nues);
   fChain->SetBranchAddress("mc_nues_id", &mc_nues_id, &b_mc_nues_id);
   fChain->SetBranchAddress("mc_nues_pt", &mc_nues_pt, &b_mc_nues_pt);
   fChain->SetBranchAddress("mc_nues_px", &mc_nues_px, &b_mc_nues_px);
   fChain->SetBranchAddress("mc_nues_py", &mc_nues_py, &b_mc_nues_py);
   fChain->SetBranchAddress("mc_nues_pz", &mc_nues_pz, &b_mc_nues_pz);
   fChain->SetBranchAddress("mc_nues_eta", &mc_nues_eta, &b_mc_nues_eta);
   fChain->SetBranchAddress("mc_nues_phi", &mc_nues_phi, &b_mc_nues_phi);
   fChain->SetBranchAddress("mc_nues_theta", &mc_nues_theta, &b_mc_nues_theta);
   fChain->SetBranchAddress("mc_nues_status", &mc_nues_status, &b_mc_nues_status);
   fChain->SetBranchAddress("mc_nues_energy", &mc_nues_energy, &b_mc_nues_energy);
   fChain->SetBranchAddress("mc_nues_charge", &mc_nues_charge, &b_mc_nues_charge);
   fChain->SetBranchAddress("mc_nues_mother_id", &mc_nues_mother_id, &b_mc_nues_mother_id);
   fChain->SetBranchAddress("mc_nues_mother_pt", &mc_nues_mother_pt, &b_mc_nues_mother_pt);
   fChain->SetBranchAddress("mc_nues_grandmother_id", &mc_nues_grandmother_id, &b_mc_nues_grandmother_id);
   fChain->SetBranchAddress("mc_nues_ggrandmother_id", &mc_nues_ggrandmother_id, &b_mc_nues_ggrandmother_id);
   fChain->SetBranchAddress("mc_nues_vertex_x", &mc_nues_vertex_x, &b_mc_nues_vertex_x);
   fChain->SetBranchAddress("mc_nues_vertex_y", &mc_nues_vertex_y, &b_mc_nues_vertex_y);
   fChain->SetBranchAddress("mc_nues_vertex_z", &mc_nues_vertex_z, &b_mc_nues_vertex_z);
   fChain->SetBranchAddress("mc_nues_mass", &mc_nues_mass, &b_mc_nues_mass);
   fChain->SetBranchAddress("mc_nues_numOfDaughters", &mc_nues_numOfDaughters, &b_mc_nues_numOfDaughters);
   fChain->SetBranchAddress("Nmc_numus", &Nmc_numus, &b_Nmc_numus);
   fChain->SetBranchAddress("mc_numus_id", &mc_numus_id, &b_mc_numus_id);
   fChain->SetBranchAddress("mc_numus_pt", &mc_numus_pt, &b_mc_numus_pt);
   fChain->SetBranchAddress("mc_numus_px", &mc_numus_px, &b_mc_numus_px);
   fChain->SetBranchAddress("mc_numus_py", &mc_numus_py, &b_mc_numus_py);
   fChain->SetBranchAddress("mc_numus_pz", &mc_numus_pz, &b_mc_numus_pz);
   fChain->SetBranchAddress("mc_numus_eta", &mc_numus_eta, &b_mc_numus_eta);
   fChain->SetBranchAddress("mc_numus_phi", &mc_numus_phi, &b_mc_numus_phi);
   fChain->SetBranchAddress("mc_numus_theta", &mc_numus_theta, &b_mc_numus_theta);
   fChain->SetBranchAddress("mc_numus_status", &mc_numus_status, &b_mc_numus_status);
   fChain->SetBranchAddress("mc_numus_energy", &mc_numus_energy, &b_mc_numus_energy);
   fChain->SetBranchAddress("mc_numus_charge", &mc_numus_charge, &b_mc_numus_charge);
   fChain->SetBranchAddress("mc_numus_mother_id", &mc_numus_mother_id, &b_mc_numus_mother_id);
   fChain->SetBranchAddress("mc_numus_mother_pt", &mc_numus_mother_pt, &b_mc_numus_mother_pt);
   fChain->SetBranchAddress("mc_numus_grandmother_id", &mc_numus_grandmother_id, &b_mc_numus_grandmother_id);
   fChain->SetBranchAddress("mc_numus_ggrandmother_id", &mc_numus_ggrandmother_id, &b_mc_numus_ggrandmother_id);
   fChain->SetBranchAddress("mc_numus_vertex_x", &mc_numus_vertex_x, &b_mc_numus_vertex_x);
   fChain->SetBranchAddress("mc_numus_vertex_y", &mc_numus_vertex_y, &b_mc_numus_vertex_y);
   fChain->SetBranchAddress("mc_numus_vertex_z", &mc_numus_vertex_z, &b_mc_numus_vertex_z);
   fChain->SetBranchAddress("mc_numus_mass", &mc_numus_mass, &b_mc_numus_mass);
   fChain->SetBranchAddress("mc_numus_numOfDaughters", &mc_numus_numOfDaughters, &b_mc_numus_numOfDaughters);
   fChain->SetBranchAddress("Nmc_nutaus", &Nmc_nutaus, &b_Nmc_nutaus);
   fChain->SetBranchAddress("mc_nutaus_id", &mc_nutaus_id, &b_mc_nutaus_id);
   fChain->SetBranchAddress("mc_nutaus_pt", &mc_nutaus_pt, &b_mc_nutaus_pt);
   fChain->SetBranchAddress("mc_nutaus_px", &mc_nutaus_px, &b_mc_nutaus_px);
   fChain->SetBranchAddress("mc_nutaus_py", &mc_nutaus_py, &b_mc_nutaus_py);
   fChain->SetBranchAddress("mc_nutaus_pz", &mc_nutaus_pz, &b_mc_nutaus_pz);
   fChain->SetBranchAddress("mc_nutaus_eta", &mc_nutaus_eta, &b_mc_nutaus_eta);
   fChain->SetBranchAddress("mc_nutaus_phi", &mc_nutaus_phi, &b_mc_nutaus_phi);
   fChain->SetBranchAddress("mc_nutaus_theta", &mc_nutaus_theta, &b_mc_nutaus_theta);
   fChain->SetBranchAddress("mc_nutaus_status", &mc_nutaus_status, &b_mc_nutaus_status);
   fChain->SetBranchAddress("mc_nutaus_energy", &mc_nutaus_energy, &b_mc_nutaus_energy);
   fChain->SetBranchAddress("mc_nutaus_charge", &mc_nutaus_charge, &b_mc_nutaus_charge);
   fChain->SetBranchAddress("mc_nutaus_mother_id", &mc_nutaus_mother_id, &b_mc_nutaus_mother_id);
   fChain->SetBranchAddress("mc_nutaus_mother_pt", &mc_nutaus_mother_pt, &b_mc_nutaus_mother_pt);
   fChain->SetBranchAddress("mc_nutaus_grandmother_id", &mc_nutaus_grandmother_id, &b_mc_nutaus_grandmother_id);
   fChain->SetBranchAddress("mc_nutaus_ggrandmother_id", &mc_nutaus_ggrandmother_id, &b_mc_nutaus_ggrandmother_id);
   fChain->SetBranchAddress("mc_nutaus_vertex_x", &mc_nutaus_vertex_x, &b_mc_nutaus_vertex_x);
   fChain->SetBranchAddress("mc_nutaus_vertex_y", &mc_nutaus_vertex_y, &b_mc_nutaus_vertex_y);
   fChain->SetBranchAddress("mc_nutaus_vertex_z", &mc_nutaus_vertex_z, &b_mc_nutaus_vertex_z);
   fChain->SetBranchAddress("mc_nutaus_mass", &mc_nutaus_mass, &b_mc_nutaus_mass);
   fChain->SetBranchAddress("mc_nutaus_numOfDaughters", &mc_nutaus_numOfDaughters, &b_mc_nutaus_numOfDaughters);
   fChain->SetBranchAddress("Nmc_photons", &Nmc_photons, &b_Nmc_photons);
   fChain->SetBranchAddress("mc_photons_id", &mc_photons_id, &b_mc_photons_id);
   fChain->SetBranchAddress("mc_photons_pt", &mc_photons_pt, &b_mc_photons_pt);
   fChain->SetBranchAddress("mc_photons_px", &mc_photons_px, &b_mc_photons_px);
   fChain->SetBranchAddress("mc_photons_py", &mc_photons_py, &b_mc_photons_py);
   fChain->SetBranchAddress("mc_photons_pz", &mc_photons_pz, &b_mc_photons_pz);
   fChain->SetBranchAddress("mc_photons_eta", &mc_photons_eta, &b_mc_photons_eta);
   fChain->SetBranchAddress("mc_photons_phi", &mc_photons_phi, &b_mc_photons_phi);
   fChain->SetBranchAddress("mc_photons_theta", &mc_photons_theta, &b_mc_photons_theta);
   fChain->SetBranchAddress("mc_photons_status", &mc_photons_status, &b_mc_photons_status);
   fChain->SetBranchAddress("mc_photons_energy", &mc_photons_energy, &b_mc_photons_energy);
   fChain->SetBranchAddress("mc_photons_charge", &mc_photons_charge, &b_mc_photons_charge);
   fChain->SetBranchAddress("mc_photons_mother_id", &mc_photons_mother_id, &b_mc_photons_mother_id);
   fChain->SetBranchAddress("mc_photons_mother_pt", &mc_photons_mother_pt, &b_mc_photons_mother_pt);
   fChain->SetBranchAddress("mc_photons_grandmother_id", &mc_photons_grandmother_id, &b_mc_photons_grandmother_id);
   fChain->SetBranchAddress("mc_photons_ggrandmother_id", &mc_photons_ggrandmother_id, &b_mc_photons_ggrandmother_id);
   fChain->SetBranchAddress("mc_photons_vertex_x", &mc_photons_vertex_x, &b_mc_photons_vertex_x);
   fChain->SetBranchAddress("mc_photons_vertex_y", &mc_photons_vertex_y, &b_mc_photons_vertex_y);
   fChain->SetBranchAddress("mc_photons_vertex_z", &mc_photons_vertex_z, &b_mc_photons_vertex_z);
   fChain->SetBranchAddress("mc_photons_mass", &mc_photons_mass, &b_mc_photons_mass);
   fChain->SetBranchAddress("mc_photons_numOfDaughters", &mc_photons_numOfDaughters, &b_mc_photons_numOfDaughters);
   fChain->SetBranchAddress("Nmc_taus", &Nmc_taus, &b_Nmc_taus);
   fChain->SetBranchAddress("mc_taus_id", &mc_taus_id, &b_mc_taus_id);
   fChain->SetBranchAddress("mc_taus_pt", &mc_taus_pt, &b_mc_taus_pt);
   fChain->SetBranchAddress("mc_taus_px", &mc_taus_px, &b_mc_taus_px);
   fChain->SetBranchAddress("mc_taus_py", &mc_taus_py, &b_mc_taus_py);
   fChain->SetBranchAddress("mc_taus_pz", &mc_taus_pz, &b_mc_taus_pz);
   fChain->SetBranchAddress("mc_taus_eta", &mc_taus_eta, &b_mc_taus_eta);
   fChain->SetBranchAddress("mc_taus_phi", &mc_taus_phi, &b_mc_taus_phi);
   fChain->SetBranchAddress("mc_taus_theta", &mc_taus_theta, &b_mc_taus_theta);
   fChain->SetBranchAddress("mc_taus_status", &mc_taus_status, &b_mc_taus_status);
   fChain->SetBranchAddress("mc_taus_energy", &mc_taus_energy, &b_mc_taus_energy);
   fChain->SetBranchAddress("mc_taus_charge", &mc_taus_charge, &b_mc_taus_charge);
   fChain->SetBranchAddress("mc_taus_mother_id", &mc_taus_mother_id, &b_mc_taus_mother_id);
   fChain->SetBranchAddress("mc_taus_mother_pt", &mc_taus_mother_pt, &b_mc_taus_mother_pt);
   fChain->SetBranchAddress("mc_taus_grandmother_id", &mc_taus_grandmother_id, &b_mc_taus_grandmother_id);
   fChain->SetBranchAddress("mc_taus_ggrandmother_id", &mc_taus_ggrandmother_id, &b_mc_taus_ggrandmother_id);
   fChain->SetBranchAddress("mc_taus_vertex_x", &mc_taus_vertex_x, &b_mc_taus_vertex_x);
   fChain->SetBranchAddress("mc_taus_vertex_y", &mc_taus_vertex_y, &b_mc_taus_vertex_y);
   fChain->SetBranchAddress("mc_taus_vertex_z", &mc_taus_vertex_z, &b_mc_taus_vertex_z);
   fChain->SetBranchAddress("mc_taus_mass", &mc_taus_mass, &b_mc_taus_mass);
   fChain->SetBranchAddress("mc_taus_numOfDaughters", &mc_taus_numOfDaughters, &b_mc_taus_numOfDaughters);
   fChain->SetBranchAddress("Nmets", &Nmets, &b_Nmets);
   fChain->SetBranchAddress("mets_et", &mets_et, &b_mets_et);
   fChain->SetBranchAddress("mets_phi", &mets_phi, &b_mets_phi);
   fChain->SetBranchAddress("mets_ex", &mets_ex, &b_mets_ex);
   fChain->SetBranchAddress("mets_ey", &mets_ey, &b_mets_ey);
   fChain->SetBranchAddress("mets_gen_et", &mets_gen_et, &b_mets_gen_et);
   fChain->SetBranchAddress("mets_gen_phi", &mets_gen_phi, &b_mets_gen_phi);
   fChain->SetBranchAddress("mets_sign", &mets_sign, &b_mets_sign);
   fChain->SetBranchAddress("mets_sumEt", &mets_sumEt, &b_mets_sumEt);
   fChain->SetBranchAddress("mets_unCPhi", &mets_unCPhi, &b_mets_unCPhi);
   fChain->SetBranchAddress("mets_unCPt", &mets_unCPt, &b_mets_unCPt);
   fChain->SetBranchAddress("Nmus", &Nmus, &b_Nmus);
   fChain->SetBranchAddress("mus_energy", &mus_energy, &b_mus_energy);
   fChain->SetBranchAddress("mus_et", &mus_et, &b_mus_et);
   fChain->SetBranchAddress("mus_eta", &mus_eta, &b_mus_eta);
   fChain->SetBranchAddress("mus_phi", &mus_phi, &b_mus_phi);
   fChain->SetBranchAddress("mus_pt", &mus_pt, &b_mus_pt);
   fChain->SetBranchAddress("mus_px", &mus_px, &b_mus_px);
   fChain->SetBranchAddress("mus_py", &mus_py, &b_mus_py);
   fChain->SetBranchAddress("mus_pz", &mus_pz, &b_mus_pz);
   fChain->SetBranchAddress("mus_status", &mus_status, &b_mus_status);
   fChain->SetBranchAddress("mus_theta", &mus_theta, &b_mus_theta);
   fChain->SetBranchAddress("mus_gen_id", &mus_gen_id, &b_mus_gen_id);
   fChain->SetBranchAddress("mus_gen_phi", &mus_gen_phi, &b_mus_gen_phi);
   fChain->SetBranchAddress("mus_gen_pt", &mus_gen_pt, &b_mus_gen_pt);
   fChain->SetBranchAddress("mus_gen_pz", &mus_gen_pz, &b_mus_gen_pz);
   fChain->SetBranchAddress("mus_gen_px", &mus_gen_px, &b_mus_gen_px);
   fChain->SetBranchAddress("mus_gen_py", &mus_gen_py, &b_mus_gen_py);
   fChain->SetBranchAddress("mus_gen_eta", &mus_gen_eta, &b_mus_gen_eta);
   fChain->SetBranchAddress("mus_gen_theta", &mus_gen_theta, &b_mus_gen_theta);
   fChain->SetBranchAddress("mus_gen_et", &mus_gen_et, &b_mus_gen_et);
   fChain->SetBranchAddress("mus_gen_mother_id", &mus_gen_mother_id, &b_mus_gen_mother_id);
   fChain->SetBranchAddress("mus_gen_mother_phi", &mus_gen_mother_phi, &b_mus_gen_mother_phi);
   fChain->SetBranchAddress("mus_gen_mother_pt", &mus_gen_mother_pt, &b_mus_gen_mother_pt);
   fChain->SetBranchAddress("mus_gen_mother_pz", &mus_gen_mother_pz, &b_mus_gen_mother_pz);
   fChain->SetBranchAddress("mus_gen_mother_px", &mus_gen_mother_px, &b_mus_gen_mother_px);
   fChain->SetBranchAddress("mus_gen_mother_py", &mus_gen_mother_py, &b_mus_gen_mother_py);
   fChain->SetBranchAddress("mus_gen_mother_eta", &mus_gen_mother_eta, &b_mus_gen_mother_eta);
   fChain->SetBranchAddress("mus_gen_mother_theta", &mus_gen_mother_theta, &b_mus_gen_mother_theta);
   fChain->SetBranchAddress("mus_gen_mother_et", &mus_gen_mother_et, &b_mus_gen_mother_et);
   fChain->SetBranchAddress("mus_tkHits", &mus_tkHits, &b_mus_tkHits);
   fChain->SetBranchAddress("mus_cIso", &mus_cIso, &b_mus_cIso);
   fChain->SetBranchAddress("mus_tIso", &mus_tIso, &b_mus_tIso);
   fChain->SetBranchAddress("mus_ecalIso", &mus_ecalIso, &b_mus_ecalIso);
   fChain->SetBranchAddress("mus_hcalIso", &mus_hcalIso, &b_mus_hcalIso);
   fChain->SetBranchAddress("mus_ecalvetoDep", &mus_ecalvetoDep, &b_mus_ecalvetoDep);
   fChain->SetBranchAddress("mus_hcalvetoDep", &mus_hcalvetoDep, &b_mus_hcalvetoDep);
   fChain->SetBranchAddress("mus_calEnergyEm", &mus_calEnergyEm, &b_mus_calEnergyEm);
   fChain->SetBranchAddress("mus_calEnergyHad", &mus_calEnergyHad, &b_mus_calEnergyHad);
   fChain->SetBranchAddress("mus_calEnergyHo", &mus_calEnergyHo, &b_mus_calEnergyHo);
   fChain->SetBranchAddress("mus_calEnergyEmS9", &mus_calEnergyEmS9, &b_mus_calEnergyEmS9);
   fChain->SetBranchAddress("mus_calEnergyHadS9", &mus_calEnergyHadS9, &b_mus_calEnergyHadS9);
   fChain->SetBranchAddress("mus_calEnergyHoS9", &mus_calEnergyHoS9, &b_mus_calEnergyHoS9);
   fChain->SetBranchAddress("mus_iso03_emVetoEt", &mus_iso03_emVetoEt, &b_mus_iso03_emVetoEt);
   fChain->SetBranchAddress("mus_iso03_hadVetoEt", &mus_iso03_hadVetoEt, &b_mus_iso03_hadVetoEt);
   fChain->SetBranchAddress("mus_iso03_sumPt", &mus_iso03_sumPt, &b_mus_iso03_sumPt);
   fChain->SetBranchAddress("mus_iso03_emEt", &mus_iso03_emEt, &b_mus_iso03_emEt);
   fChain->SetBranchAddress("mus_iso03_hadEt", &mus_iso03_hadEt, &b_mus_iso03_hadEt);
   fChain->SetBranchAddress("mus_iso03_hoEt", &mus_iso03_hoEt, &b_mus_iso03_hoEt);
   fChain->SetBranchAddress("mus_iso03_nTracks", &mus_iso03_nTracks, &b_mus_iso03_nTracks);
   fChain->SetBranchAddress("mus_iso05_sumPt", &mus_iso05_sumPt, &b_mus_iso05_sumPt);
   fChain->SetBranchAddress("mus_iso05_emEt", &mus_iso05_emEt, &b_mus_iso05_emEt);
   fChain->SetBranchAddress("mus_iso05_hadEt", &mus_iso05_hadEt, &b_mus_iso05_hadEt);
   fChain->SetBranchAddress("mus_iso05_hoEt", &mus_iso05_hoEt, &b_mus_iso05_hoEt);
   fChain->SetBranchAddress("mus_iso05_nTracks", &mus_iso05_nTracks, &b_mus_iso05_nTracks);
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
   fChain->SetBranchAddress("mus_cm_chg", &mus_cm_chg, &b_mus_cm_chg);
   fChain->SetBranchAddress("mus_cm_pt", &mus_cm_pt, &b_mus_cm_pt);
   fChain->SetBranchAddress("mus_cm_px", &mus_cm_px, &b_mus_cm_px);
   fChain->SetBranchAddress("mus_cm_py", &mus_cm_py, &b_mus_cm_py);
   fChain->SetBranchAddress("mus_cm_pz", &mus_cm_pz, &b_mus_cm_pz);
   fChain->SetBranchAddress("mus_cm_eta", &mus_cm_eta, &b_mus_cm_eta);
   fChain->SetBranchAddress("mus_cm_phi", &mus_cm_phi, &b_mus_cm_phi);
   fChain->SetBranchAddress("mus_cm_theta", &mus_cm_theta, &b_mus_cm_theta);
   fChain->SetBranchAddress("mus_cm_d0dum", &mus_cm_d0dum, &b_mus_cm_d0dum);
   fChain->SetBranchAddress("mus_cm_dz", &mus_cm_dz, &b_mus_cm_dz);
   fChain->SetBranchAddress("mus_cm_vx", &mus_cm_vx, &b_mus_cm_vx);
   fChain->SetBranchAddress("mus_cm_vy", &mus_cm_vy, &b_mus_cm_vy);
   fChain->SetBranchAddress("mus_cm_vz", &mus_cm_vz, &b_mus_cm_vz);
   fChain->SetBranchAddress("mus_cm_numvalhits", &mus_cm_numvalhits, &b_mus_cm_numvalhits);
   fChain->SetBranchAddress("mus_cm_numlosthits", &mus_cm_numlosthits, &b_mus_cm_numlosthits);
   fChain->SetBranchAddress("mus_cm_numvalMuonhits", &mus_cm_numvalMuonhits, &b_mus_cm_numvalMuonhits);
   fChain->SetBranchAddress("mus_cm_d0dumErr", &mus_cm_d0dumErr, &b_mus_cm_d0dumErr);
   fChain->SetBranchAddress("mus_cm_dzErr", &mus_cm_dzErr, &b_mus_cm_dzErr);
   fChain->SetBranchAddress("mus_cm_ptErr", &mus_cm_ptErr, &b_mus_cm_ptErr);
   fChain->SetBranchAddress("mus_cm_etaErr", &mus_cm_etaErr, &b_mus_cm_etaErr);
   fChain->SetBranchAddress("mus_cm_phiErr", &mus_cm_phiErr, &b_mus_cm_phiErr);
   fChain->SetBranchAddress("mus_tk_id", &mus_tk_id, &b_mus_tk_id);
   fChain->SetBranchAddress("mus_tk_chi2", &mus_tk_chi2, &b_mus_tk_chi2);
   fChain->SetBranchAddress("mus_tk_ndof", &mus_tk_ndof, &b_mus_tk_ndof);
   fChain->SetBranchAddress("mus_tk_chg", &mus_tk_chg, &b_mus_tk_chg);
   fChain->SetBranchAddress("mus_tk_pt", &mus_tk_pt, &b_mus_tk_pt);
   fChain->SetBranchAddress("mus_tk_px", &mus_tk_px, &b_mus_tk_px);
   fChain->SetBranchAddress("mus_tk_py", &mus_tk_py, &b_mus_tk_py);
   fChain->SetBranchAddress("mus_tk_pz", &mus_tk_pz, &b_mus_tk_pz);
   fChain->SetBranchAddress("mus_tk_eta", &mus_tk_eta, &b_mus_tk_eta);
   fChain->SetBranchAddress("mus_tk_phi", &mus_tk_phi, &b_mus_tk_phi);
   fChain->SetBranchAddress("mus_tk_theta", &mus_tk_theta, &b_mus_tk_theta);
   fChain->SetBranchAddress("mus_tk_d0dum", &mus_tk_d0dum, &b_mus_tk_d0dum);
   fChain->SetBranchAddress("mus_tk_dz", &mus_tk_dz, &b_mus_tk_dz);
   fChain->SetBranchAddress("mus_tk_vx", &mus_tk_vx, &b_mus_tk_vx);
   fChain->SetBranchAddress("mus_tk_vy", &mus_tk_vy, &b_mus_tk_vy);
   fChain->SetBranchAddress("mus_tk_vz", &mus_tk_vz, &b_mus_tk_vz);
   fChain->SetBranchAddress("mus_tk_numvalhits", &mus_tk_numvalhits, &b_mus_tk_numvalhits);
   fChain->SetBranchAddress("mus_tk_numlosthits", &mus_tk_numlosthits, &b_mus_tk_numlosthits);
   fChain->SetBranchAddress("mus_tk_d0dumErr", &mus_tk_d0dumErr, &b_mus_tk_d0dumErr);
   fChain->SetBranchAddress("mus_tk_dzErr", &mus_tk_dzErr, &b_mus_tk_dzErr);
   fChain->SetBranchAddress("mus_tk_ptErr", &mus_tk_ptErr, &b_mus_tk_ptErr);
   fChain->SetBranchAddress("mus_tk_etaErr", &mus_tk_etaErr, &b_mus_tk_etaErr);
   fChain->SetBranchAddress("mus_tk_phiErr", &mus_tk_phiErr, &b_mus_tk_phiErr);
   fChain->SetBranchAddress("mus_tk_numvalPixelhits", &mus_tk_numvalPixelhits, &b_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("mus_tk_numpixelWthMeasr", &mus_tk_numpixelWthMeasr, &b_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("mus_stamu_chi2", &mus_stamu_chi2, &b_mus_stamu_chi2);
   fChain->SetBranchAddress("mus_stamu_ndof", &mus_stamu_ndof, &b_mus_stamu_ndof);
   fChain->SetBranchAddress("mus_stamu_chg", &mus_stamu_chg, &b_mus_stamu_chg);
   fChain->SetBranchAddress("mus_stamu_pt", &mus_stamu_pt, &b_mus_stamu_pt);
   fChain->SetBranchAddress("mus_stamu_px", &mus_stamu_px, &b_mus_stamu_px);
   fChain->SetBranchAddress("mus_stamu_py", &mus_stamu_py, &b_mus_stamu_py);
   fChain->SetBranchAddress("mus_stamu_pz", &mus_stamu_pz, &b_mus_stamu_pz);
   fChain->SetBranchAddress("mus_stamu_eta", &mus_stamu_eta, &b_mus_stamu_eta);
   fChain->SetBranchAddress("mus_stamu_phi", &mus_stamu_phi, &b_mus_stamu_phi);
   fChain->SetBranchAddress("mus_stamu_theta", &mus_stamu_theta, &b_mus_stamu_theta);
   fChain->SetBranchAddress("mus_stamu_d0dum", &mus_stamu_d0dum, &b_mus_stamu_d0dum);
   fChain->SetBranchAddress("mus_stamu_dz", &mus_stamu_dz, &b_mus_stamu_dz);
   fChain->SetBranchAddress("mus_stamu_vx", &mus_stamu_vx, &b_mus_stamu_vx);
   fChain->SetBranchAddress("mus_stamu_vy", &mus_stamu_vy, &b_mus_stamu_vy);
   fChain->SetBranchAddress("mus_stamu_vz", &mus_stamu_vz, &b_mus_stamu_vz);
   fChain->SetBranchAddress("mus_stamu_numvalhits", &mus_stamu_numvalhits, &b_mus_stamu_numvalhits);
   fChain->SetBranchAddress("mus_stamu_numlosthits", &mus_stamu_numlosthits, &b_mus_stamu_numlosthits);
   fChain->SetBranchAddress("mus_stamu_d0dumErr", &mus_stamu_d0dumErr, &b_mus_stamu_d0dumErr);
   fChain->SetBranchAddress("mus_stamu_dzErr", &mus_stamu_dzErr, &b_mus_stamu_dzErr);
   fChain->SetBranchAddress("mus_stamu_ptErr", &mus_stamu_ptErr, &b_mus_stamu_ptErr);
   fChain->SetBranchAddress("mus_stamu_etaErr", &mus_stamu_etaErr, &b_mus_stamu_etaErr);
   fChain->SetBranchAddress("mus_stamu_phiErr", &mus_stamu_phiErr, &b_mus_stamu_phiErr);
   fChain->SetBranchAddress("mus_num_matches", &mus_num_matches, &b_mus_num_matches);
   fChain->SetBranchAddress("mus_isPFMuon", &mus_isPFMuon, &b_mus_isPFMuon);
   fChain->SetBranchAddress("mus_isTrackerMuon", &mus_isTrackerMuon, &b_mus_isTrackerMuon);
   fChain->SetBranchAddress("mus_isStandAloneMuon", &mus_isStandAloneMuon, &b_mus_isStandAloneMuon);
   fChain->SetBranchAddress("mus_isCaloMuon", &mus_isCaloMuon, &b_mus_isCaloMuon);
   fChain->SetBranchAddress("mus_isGlobalMuon", &mus_isGlobalMuon, &b_mus_isGlobalMuon);
   fChain->SetBranchAddress("mus_isElectron", &mus_isElectron, &b_mus_isElectron);
   fChain->SetBranchAddress("mus_isConvertedPhoton", &mus_isConvertedPhoton, &b_mus_isConvertedPhoton);
   fChain->SetBranchAddress("mus_isPhoton", &mus_isPhoton, &b_mus_isPhoton);
   fChain->SetBranchAddress("mus_id_All", &mus_id_All, &b_mus_id_All);
   fChain->SetBranchAddress("mus_id_AllGlobalMuons", &mus_id_AllGlobalMuons, &b_mus_id_AllGlobalMuons);
   fChain->SetBranchAddress("mus_id_AllStandAloneMuons", &mus_id_AllStandAloneMuons, &b_mus_id_AllStandAloneMuons);
   fChain->SetBranchAddress("mus_id_AllTrackerMuons", &mus_id_AllTrackerMuons, &b_mus_id_AllTrackerMuons);
   fChain->SetBranchAddress("mus_id_TrackerMuonArbitrated", &mus_id_TrackerMuonArbitrated, &b_mus_id_TrackerMuonArbitrated);
   fChain->SetBranchAddress("mus_id_AllArbitrated", &mus_id_AllArbitrated, &b_mus_id_AllArbitrated);
   fChain->SetBranchAddress("mus_id_GlobalMuonPromptTight", &mus_id_GlobalMuonPromptTight, &b_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("mus_id_TMLastStationLoose", &mus_id_TMLastStationLoose, &b_mus_id_TMLastStationLoose);
   fChain->SetBranchAddress("mus_id_TMLastStationTight", &mus_id_TMLastStationTight, &b_mus_id_TMLastStationTight);
   fChain->SetBranchAddress("mus_id_TM2DCompatibilityLoose", &mus_id_TM2DCompatibilityLoose, &b_mus_id_TM2DCompatibilityLoose);
   fChain->SetBranchAddress("mus_id_TM2DCompatibilityTight", &mus_id_TM2DCompatibilityTight, &b_mus_id_TM2DCompatibilityTight);
   fChain->SetBranchAddress("mus_id_TMOneStationLoose", &mus_id_TMOneStationLoose, &b_mus_id_TMOneStationLoose);
   fChain->SetBranchAddress("mus_id_TMOneStationTight", &mus_id_TMOneStationTight, &b_mus_id_TMOneStationTight);
   fChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtLoose", &mus_id_TMLastStationOptimizedLowPtLoose, &b_mus_id_TMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtTight", &mus_id_TMLastStationOptimizedLowPtTight, &b_mus_id_TMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("mus_tk_LayersWithMeasurement", &mus_tk_LayersWithMeasurement, &b_mus_tk_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_tk_PixelLayersWithMeasurement", &mus_tk_PixelLayersWithMeasurement, &b_mus_tk_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_tk_ValidStripLayersWithMonoAndStereoHit", &mus_tk_ValidStripLayersWithMonoAndStereoHit, &b_mus_tk_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_tk_LayersWithoutMeasurement", &mus_tk_LayersWithoutMeasurement, &b_mus_tk_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_tk_ExpectedHitsInner", &mus_tk_ExpectedHitsInner, &b_mus_tk_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_tk_ExpectedHitsOuter", &mus_tk_ExpectedHitsOuter, &b_mus_tk_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_cm_LayersWithMeasurement", &mus_cm_LayersWithMeasurement, &b_mus_cm_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_cm_PixelLayersWithMeasurement", &mus_cm_PixelLayersWithMeasurement, &b_mus_cm_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_cm_ValidStripLayersWithMonoAndStereoHit", &mus_cm_ValidStripLayersWithMonoAndStereoHit, &b_mus_cm_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_cm_LayersWithoutMeasurement", &mus_cm_LayersWithoutMeasurement, &b_mus_cm_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_cm_ExpectedHitsInner", &mus_cm_ExpectedHitsInner, &b_mus_cm_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_cm_ExpectedHitsOuter", &mus_cm_ExpectedHitsOuter, &b_mus_cm_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_picky_LayersWithMeasurement", &mus_picky_LayersWithMeasurement, &b_mus_picky_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_picky_PixelLayersWithMeasurement", &mus_picky_PixelLayersWithMeasurement, &b_mus_picky_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_picky_ValidStripLayersWithMonoAndStereoHit", &mus_picky_ValidStripLayersWithMonoAndStereoHit, &b_mus_picky_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_picky_LayersWithoutMeasurement", &mus_picky_LayersWithoutMeasurement, &b_mus_picky_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_picky_ExpectedHitsInner", &mus_picky_ExpectedHitsInner, &b_mus_picky_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_picky_ExpectedHitsOuter", &mus_picky_ExpectedHitsOuter, &b_mus_picky_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_tpfms_LayersWithMeasurement", &mus_tpfms_LayersWithMeasurement, &b_mus_tpfms_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_tpfms_PixelLayersWithMeasurement", &mus_tpfms_PixelLayersWithMeasurement, &b_mus_tpfms_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_tpfms_ValidStripLayersWithMonoAndStereoHit", &mus_tpfms_ValidStripLayersWithMonoAndStereoHit, &b_mus_tpfms_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_tpfms_LayersWithoutMeasurement", &mus_tpfms_LayersWithoutMeasurement, &b_mus_tpfms_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_tpfms_ExpectedHitsInner", &mus_tpfms_ExpectedHitsInner, &b_mus_tpfms_ExpectedHitsInner);
   fChain->SetBranchAddress("mus_tpfms_ExpectedHitsOuter", &mus_tpfms_ExpectedHitsOuter, &b_mus_tpfms_ExpectedHitsOuter);
   fChain->SetBranchAddress("mus_picky_id", &mus_picky_id, &b_mus_picky_id);
   fChain->SetBranchAddress("mus_picky_chi2", &mus_picky_chi2, &b_mus_picky_chi2);
   fChain->SetBranchAddress("mus_picky_ndof", &mus_picky_ndof, &b_mus_picky_ndof);
   fChain->SetBranchAddress("mus_picky_chg", &mus_picky_chg, &b_mus_picky_chg);
   fChain->SetBranchAddress("mus_picky_pt", &mus_picky_pt, &b_mus_picky_pt);
   fChain->SetBranchAddress("mus_picky_px", &mus_picky_px, &b_mus_picky_px);
   fChain->SetBranchAddress("mus_picky_py", &mus_picky_py, &b_mus_picky_py);
   fChain->SetBranchAddress("mus_picky_pz", &mus_picky_pz, &b_mus_picky_pz);
   fChain->SetBranchAddress("mus_picky_eta", &mus_picky_eta, &b_mus_picky_eta);
   fChain->SetBranchAddress("mus_picky_phi", &mus_picky_phi, &b_mus_picky_phi);
   fChain->SetBranchAddress("mus_picky_theta", &mus_picky_theta, &b_mus_picky_theta);
   fChain->SetBranchAddress("mus_picky_d0dum", &mus_picky_d0dum, &b_mus_picky_d0dum);
   fChain->SetBranchAddress("mus_picky_dz", &mus_picky_dz, &b_mus_picky_dz);
   fChain->SetBranchAddress("mus_picky_vx", &mus_picky_vx, &b_mus_picky_vx);
   fChain->SetBranchAddress("mus_picky_vy", &mus_picky_vy, &b_mus_picky_vy);
   fChain->SetBranchAddress("mus_picky_vz", &mus_picky_vz, &b_mus_picky_vz);
   fChain->SetBranchAddress("mus_picky_numvalhits", &mus_picky_numvalhits, &b_mus_picky_numvalhits);
   fChain->SetBranchAddress("mus_picky_numlosthits", &mus_picky_numlosthits, &b_mus_picky_numlosthits);
   fChain->SetBranchAddress("mus_picky_d0dumErr", &mus_picky_d0dumErr, &b_mus_picky_d0dumErr);
   fChain->SetBranchAddress("mus_picky_dzErr", &mus_picky_dzErr, &b_mus_picky_dzErr);
   fChain->SetBranchAddress("mus_picky_ptErr", &mus_picky_ptErr, &b_mus_picky_ptErr);
   fChain->SetBranchAddress("mus_picky_etaErr", &mus_picky_etaErr, &b_mus_picky_etaErr);
   fChain->SetBranchAddress("mus_picky_phiErr", &mus_picky_phiErr, &b_mus_picky_phiErr);
   fChain->SetBranchAddress("mus_picky_numvalPixelhits", &mus_picky_numvalPixelhits, &b_mus_picky_numvalPixelhits);
   fChain->SetBranchAddress("mus_tpfms_id", &mus_tpfms_id, &b_mus_tpfms_id);
   fChain->SetBranchAddress("mus_tpfms_chi2", &mus_tpfms_chi2, &b_mus_tpfms_chi2);
   fChain->SetBranchAddress("mus_tpfms_ndof", &mus_tpfms_ndof, &b_mus_tpfms_ndof);
   fChain->SetBranchAddress("mus_tpfms_chg", &mus_tpfms_chg, &b_mus_tpfms_chg);
   fChain->SetBranchAddress("mus_tpfms_pt", &mus_tpfms_pt, &b_mus_tpfms_pt);
   fChain->SetBranchAddress("mus_tpfms_px", &mus_tpfms_px, &b_mus_tpfms_px);
   fChain->SetBranchAddress("mus_tpfms_py", &mus_tpfms_py, &b_mus_tpfms_py);
   fChain->SetBranchAddress("mus_tpfms_pz", &mus_tpfms_pz, &b_mus_tpfms_pz);
   fChain->SetBranchAddress("mus_tpfms_eta", &mus_tpfms_eta, &b_mus_tpfms_eta);
   fChain->SetBranchAddress("mus_tpfms_phi", &mus_tpfms_phi, &b_mus_tpfms_phi);
   fChain->SetBranchAddress("mus_tpfms_theta", &mus_tpfms_theta, &b_mus_tpfms_theta);
   fChain->SetBranchAddress("mus_tpfms_d0dum", &mus_tpfms_d0dum, &b_mus_tpfms_d0dum);
   fChain->SetBranchAddress("mus_tpfms_dz", &mus_tpfms_dz, &b_mus_tpfms_dz);
   fChain->SetBranchAddress("mus_tpfms_vx", &mus_tpfms_vx, &b_mus_tpfms_vx);
   fChain->SetBranchAddress("mus_tpfms_vy", &mus_tpfms_vy, &b_mus_tpfms_vy);
   fChain->SetBranchAddress("mus_tpfms_vz", &mus_tpfms_vz, &b_mus_tpfms_vz);
   fChain->SetBranchAddress("mus_tpfms_numvalhits", &mus_tpfms_numvalhits, &b_mus_tpfms_numvalhits);
   fChain->SetBranchAddress("mus_tpfms_numlosthits", &mus_tpfms_numlosthits, &b_mus_tpfms_numlosthits);
   fChain->SetBranchAddress("mus_tpfms_d0dumErr", &mus_tpfms_d0dumErr, &b_mus_tpfms_d0dumErr);
   fChain->SetBranchAddress("mus_tpfms_dzErr", &mus_tpfms_dzErr, &b_mus_tpfms_dzErr);
   fChain->SetBranchAddress("mus_tpfms_ptErr", &mus_tpfms_ptErr, &b_mus_tpfms_ptErr);
   fChain->SetBranchAddress("mus_tpfms_etaErr", &mus_tpfms_etaErr, &b_mus_tpfms_etaErr);
   fChain->SetBranchAddress("mus_tpfms_phiErr", &mus_tpfms_phiErr, &b_mus_tpfms_phiErr);
   fChain->SetBranchAddress("mus_tpfms_numvalPixelhits", &mus_tpfms_numvalPixelhits, &b_mus_tpfms_numvalPixelhits);
   fChain->SetBranchAddress("mus_dB", &mus_dB, &b_mus_dB);
   fChain->SetBranchAddress("mus_numberOfMatchedStations", &mus_numberOfMatchedStations, &b_mus_numberOfMatchedStations);
   fChain->SetBranchAddress("Npfcand", &Npfcand, &b_Npfcand);
   fChain->SetBranchAddress("pfcand_pdgId", &pfcand_pdgId, &b_pfcand_pdgId);
   fChain->SetBranchAddress("pfcand_pt", &pfcand_pt, &b_pfcand_pt);
   fChain->SetBranchAddress("pfcand_pz", &pfcand_pz, &b_pfcand_pz);
   fChain->SetBranchAddress("pfcand_px", &pfcand_px, &b_pfcand_px);
   fChain->SetBranchAddress("pfcand_py", &pfcand_py, &b_pfcand_py);
   fChain->SetBranchAddress("pfcand_eta", &pfcand_eta, &b_pfcand_eta);
   fChain->SetBranchAddress("pfcand_phi", &pfcand_phi, &b_pfcand_phi);
   fChain->SetBranchAddress("pfcand_theta", &pfcand_theta, &b_pfcand_theta);
   fChain->SetBranchAddress("pfcand_energy", &pfcand_energy, &b_pfcand_energy);
   fChain->SetBranchAddress("pfcand_charge", &pfcand_charge, &b_pfcand_charge);
   fChain->SetBranchAddress("Nphotons", &Nphotons, &b_Nphotons);
   fChain->SetBranchAddress("photons_energy", &photons_energy, &b_photons_energy);
   fChain->SetBranchAddress("photons_et", &photons_et, &b_photons_et);
   fChain->SetBranchAddress("photons_eta", &photons_eta, &b_photons_eta);
   fChain->SetBranchAddress("photons_phi", &photons_phi, &b_photons_phi);
   fChain->SetBranchAddress("photons_pt", &photons_pt, &b_photons_pt);
   fChain->SetBranchAddress("photons_px", &photons_px, &b_photons_px);
   fChain->SetBranchAddress("photons_py", &photons_py, &b_photons_py);
   fChain->SetBranchAddress("photons_pz", &photons_pz, &b_photons_pz);
   fChain->SetBranchAddress("photons_status", &photons_status, &b_photons_status);
   fChain->SetBranchAddress("photons_theta", &photons_theta, &b_photons_theta);
   fChain->SetBranchAddress("photons_hadOverEM", &photons_hadOverEM, &b_photons_hadOverEM);
   fChain->SetBranchAddress("photons_hadTowOverEM", &photons_hadTowOverEM, &b_photons_hadTowOverEM);
   fChain->SetBranchAddress("photons_scEnergy", &photons_scEnergy, &b_photons_scEnergy);
   fChain->SetBranchAddress("photons_scRawEnergy", &photons_scRawEnergy, &b_photons_scRawEnergy);
   fChain->SetBranchAddress("photons_scEta", &photons_scEta, &b_photons_scEta);
   fChain->SetBranchAddress("photons_scPhi", &photons_scPhi, &b_photons_scPhi);
   fChain->SetBranchAddress("photons_scEtaWidth", &photons_scEtaWidth, &b_photons_scEtaWidth);
   fChain->SetBranchAddress("photons_scPhiWidth", &photons_scPhiWidth, &b_photons_scPhiWidth);
   fChain->SetBranchAddress("photons_tIso", &photons_tIso, &b_photons_tIso);
   fChain->SetBranchAddress("photons_ecalIso", &photons_ecalIso, &b_photons_ecalIso);
   fChain->SetBranchAddress("photons_hcalIso", &photons_hcalIso, &b_photons_hcalIso);
   fChain->SetBranchAddress("photons_isoEcalRecHitDR04", &photons_isoEcalRecHitDR04, &b_photons_isoEcalRecHitDR04);
   fChain->SetBranchAddress("photons_isoHcalRecHitDR04", &photons_isoHcalRecHitDR04, &b_photons_isoHcalRecHitDR04);
   fChain->SetBranchAddress("photons_isoSolidTrkConeDR04", &photons_isoSolidTrkConeDR04, &b_photons_isoSolidTrkConeDR04);
   fChain->SetBranchAddress("photons_isoHollowTrkConeDR04", &photons_isoHollowTrkConeDR04, &b_photons_isoHollowTrkConeDR04);
   fChain->SetBranchAddress("photons_nTrkSolidConeDR04", &photons_nTrkSolidConeDR04, &b_photons_nTrkSolidConeDR04);
   fChain->SetBranchAddress("photons_nTrkHollowConeDR04", &photons_nTrkHollowConeDR04, &b_photons_nTrkHollowConeDR04);
   fChain->SetBranchAddress("photons_isoEcalRecHitDR03", &photons_isoEcalRecHitDR03, &b_photons_isoEcalRecHitDR03);
   fChain->SetBranchAddress("photons_isoHcalRecHitDR03", &photons_isoHcalRecHitDR03, &b_photons_isoHcalRecHitDR03);
   fChain->SetBranchAddress("photons_isoSolidTrkConeDR03", &photons_isoSolidTrkConeDR03, &b_photons_isoSolidTrkConeDR03);
   fChain->SetBranchAddress("photons_isoHollowTrkConeDR03", &photons_isoHollowTrkConeDR03, &b_photons_isoHollowTrkConeDR03);
   fChain->SetBranchAddress("photons_nTrkSolidConeDR03", &photons_nTrkSolidConeDR03, &b_photons_nTrkSolidConeDR03);
   fChain->SetBranchAddress("photons_nTrkHollowConeDR03", &photons_nTrkHollowConeDR03, &b_photons_nTrkHollowConeDR03);
   fChain->SetBranchAddress("photons_isAlsoElectron", &photons_isAlsoElectron, &b_photons_isAlsoElectron);
   fChain->SetBranchAddress("photons_hasPixelSeed", &photons_hasPixelSeed, &b_photons_hasPixelSeed);
   fChain->SetBranchAddress("photons_isConverted", &photons_isConverted, &b_photons_isConverted);
   fChain->SetBranchAddress("photons_isEBGap", &photons_isEBGap, &b_photons_isEBGap);
   fChain->SetBranchAddress("photons_isEEGap", &photons_isEEGap, &b_photons_isEEGap);
   fChain->SetBranchAddress("photons_isEBEEGap", &photons_isEBEEGap, &b_photons_isEBEEGap);
   fChain->SetBranchAddress("photons_isEBPho", &photons_isEBPho, &b_photons_isEBPho);
   fChain->SetBranchAddress("photons_isEEPho", &photons_isEEPho, &b_photons_isEEPho);
   fChain->SetBranchAddress("photons_isLoosePhoton", &photons_isLoosePhoton, &b_photons_isLoosePhoton);
   fChain->SetBranchAddress("photons_isTightPhoton", &photons_isTightPhoton, &b_photons_isTightPhoton);
   fChain->SetBranchAddress("photons_maxEnergyXtal", &photons_maxEnergyXtal, &b_photons_maxEnergyXtal);
   fChain->SetBranchAddress("photons_e1x5", &photons_e1x5, &b_photons_e1x5);
   fChain->SetBranchAddress("photons_e2x5", &photons_e2x5, &b_photons_e2x5);
   fChain->SetBranchAddress("photons_e3x3", &photons_e3x3, &b_photons_e3x3);
   fChain->SetBranchAddress("photons_e5x5", &photons_e5x5, &b_photons_e5x5);
   fChain->SetBranchAddress("photons_sigmaEtaEta", &photons_sigmaEtaEta, &b_photons_sigmaEtaEta);
   fChain->SetBranchAddress("photons_sigmaIetaIeta", &photons_sigmaIetaIeta, &b_photons_sigmaIetaIeta);
   fChain->SetBranchAddress("photons_r9", &photons_r9, &b_photons_r9);
   fChain->SetBranchAddress("photons_gen_et", &photons_gen_et, &b_photons_gen_et);
   fChain->SetBranchAddress("photons_gen_eta", &photons_gen_eta, &b_photons_gen_eta);
   fChain->SetBranchAddress("photons_gen_phi", &photons_gen_phi, &b_photons_gen_phi);
   fChain->SetBranchAddress("photons_gen_id", &photons_gen_id, &b_photons_gen_id);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_xErr", &pv_xErr, &b_pv_xErr);
   fChain->SetBranchAddress("pv_yErr", &pv_yErr, &b_pv_yErr);
   fChain->SetBranchAddress("pv_zErr", &pv_zErr, &b_pv_zErr);
   fChain->SetBranchAddress("pv_chi2", &pv_chi2, &b_pv_chi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("pv_isValid", &pv_isValid, &b_pv_isValid);
   fChain->SetBranchAddress("pv_tracksSize", &pv_tracksSize, &b_pv_tracksSize);
   fChain->SetBranchAddress("Ntaus", &Ntaus, &b_Ntaus);
   fChain->SetBranchAddress("taus_status", &taus_status, &b_taus_status);
   fChain->SetBranchAddress("taus_phi", &taus_phi, &b_taus_phi);
   fChain->SetBranchAddress("taus_pt", &taus_pt, &b_taus_pt);
   fChain->SetBranchAddress("taus_pz", &taus_pz, &b_taus_pz);
   fChain->SetBranchAddress("taus_px", &taus_px, &b_taus_px);
   fChain->SetBranchAddress("taus_py", &taus_py, &b_taus_py);
   fChain->SetBranchAddress("taus_eta", &taus_eta, &b_taus_eta);
   fChain->SetBranchAddress("taus_theta", &taus_theta, &b_taus_theta);
   fChain->SetBranchAddress("taus_et", &taus_et, &b_taus_et);
   fChain->SetBranchAddress("taus_energy", &taus_energy, &b_taus_energy);
   fChain->SetBranchAddress("taus_charge", &taus_charge, &b_taus_charge);
   fChain->SetBranchAddress("taus_emf", &taus_emf, &b_taus_emf);
   fChain->SetBranchAddress("taus_hcalTotOverPLead", &taus_hcalTotOverPLead, &b_taus_hcalTotOverPLead);
   fChain->SetBranchAddress("taus_hcalMaxOverPLead", &taus_hcalMaxOverPLead, &b_taus_hcalMaxOverPLead);
   fChain->SetBranchAddress("taus_hcal3x3OverPLead", &taus_hcal3x3OverPLead, &b_taus_hcal3x3OverPLead);
   fChain->SetBranchAddress("taus_ecalStripSumEOverPLead", &taus_ecalStripSumEOverPLead, &b_taus_ecalStripSumEOverPLead);
   fChain->SetBranchAddress("taus_elecPreIdOutput", &taus_elecPreIdOutput, &b_taus_elecPreIdOutput);
   fChain->SetBranchAddress("taus_elecPreIdDecision", &taus_elecPreIdDecision, &b_taus_elecPreIdDecision);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_pt", &taus_leadPFChargedHadrCand_pt, &b_taus_leadPFChargedHadrCand_pt);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_charge", &taus_leadPFChargedHadrCand_charge, &b_taus_leadPFChargedHadrCand_charge);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_eta", &taus_leadPFChargedHadrCand_eta, &b_taus_leadPFChargedHadrCand_eta);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_ECAL_eta", &taus_leadPFChargedHadrCand_ECAL_eta, &b_taus_leadPFChargedHadrCand_ECAL_eta);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_phi", &taus_leadPFChargedHadrCand_phi, &b_taus_leadPFChargedHadrCand_phi);
   fChain->SetBranchAddress("taus_isoPFGammaCandsEtSum", &taus_isoPFGammaCandsEtSum, &b_taus_isoPFGammaCandsEtSum);
   fChain->SetBranchAddress("taus_isoPFChargedHadrCandsPtSum", &taus_isoPFChargedHadrCandsPtSum, &b_taus_isoPFChargedHadrCandsPtSum);
   fChain->SetBranchAddress("taus_leadingTrackFinding", &taus_leadingTrackFinding, &b_taus_leadingTrackFinding);
   fChain->SetBranchAddress("taus_leadingTrackPtCut", &taus_leadingTrackPtCut, &b_taus_leadingTrackPtCut);
   fChain->SetBranchAddress("taus_trackIsolation", &taus_trackIsolation, &b_taus_trackIsolation);
   fChain->SetBranchAddress("taus_ecalIsolation", &taus_ecalIsolation, &b_taus_ecalIsolation);
   fChain->SetBranchAddress("taus_byIsolation", &taus_byIsolation, &b_taus_byIsolation);
   fChain->SetBranchAddress("taus_againstElectron", &taus_againstElectron, &b_taus_againstElectron);
   fChain->SetBranchAddress("taus_againstMuon", &taus_againstMuon, &b_taus_againstMuon);
   fChain->SetBranchAddress("taus_taNC_quarter", &taus_taNC_quarter, &b_taus_taNC_quarter);
   fChain->SetBranchAddress("taus_taNC_one", &taus_taNC_one, &b_taus_taNC_one);
   fChain->SetBranchAddress("taus_taNC_half", &taus_taNC_half, &b_taus_taNC_half);
   fChain->SetBranchAddress("taus_taNC_tenth", &taus_taNC_tenth, &b_taus_taNC_tenth);
   fChain->SetBranchAddress("taus_taNC", &taus_taNC, &b_taus_taNC);
   fChain->SetBranchAddress("taus_byIsoUsingLeadingPi", &taus_byIsoUsingLeadingPi, &b_taus_byIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_tkIsoUsingLeadingPi", &taus_tkIsoUsingLeadingPi, &b_taus_tkIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_ecalIsoUsingLeadingPi", &taus_ecalIsoUsingLeadingPi, &b_taus_ecalIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_againstElectronLoose", &taus_againstElectronLoose, &b_taus_againstElectronLoose);
   fChain->SetBranchAddress("taus_againstElectronMedium", &taus_againstElectronMedium, &b_taus_againstElectronMedium);
   fChain->SetBranchAddress("taus_againstElectronTight", &taus_againstElectronTight, &b_taus_againstElectronTight);
   fChain->SetBranchAddress("taus_againstElectronMVA", &taus_againstElectronMVA, &b_taus_againstElectronMVA);
   fChain->SetBranchAddress("taus_againstMuonLoose", &taus_againstMuonLoose, &b_taus_againstMuonLoose);
   fChain->SetBranchAddress("taus_againstMuonMedium", &taus_againstMuonMedium, &b_taus_againstMuonMedium);
   fChain->SetBranchAddress("taus_againstMuonTight", &taus_againstMuonTight, &b_taus_againstMuonTight);
   fChain->SetBranchAddress("taus_decayModeFinding", &taus_decayModeFinding, &b_taus_decayModeFinding);
   fChain->SetBranchAddress("taus_byVLooseIsolation", &taus_byVLooseIsolation, &b_taus_byVLooseIsolation);
   fChain->SetBranchAddress("taus_byLooseIsolation", &taus_byLooseIsolation, &b_taus_byLooseIsolation);
   fChain->SetBranchAddress("taus_byMediumIsolation", &taus_byMediumIsolation, &b_taus_byMediumIsolation);
   fChain->SetBranchAddress("taus_byTightIsolation", &taus_byTightIsolation, &b_taus_byTightIsolation);
   fChain->SetBranchAddress("taus_byVLooseIsolationDeltaBetaCorr", &taus_byVLooseIsolationDeltaBetaCorr, &b_taus_byVLooseIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_byLooseIsolationDeltaBetaCorr", &taus_byLooseIsolationDeltaBetaCorr, &b_taus_byLooseIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_byMediumIsolationDeltaBetaCorr", &taus_byMediumIsolationDeltaBetaCorr, &b_taus_byMediumIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_byTightIsolationDeltaBetaCorr", &taus_byTightIsolationDeltaBetaCorr, &b_taus_byTightIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("taus_signalPFChargedHadrCandsSize", &taus_signalPFChargedHadrCandsSize, &b_taus_signalPFChargedHadrCandsSize);
   fChain->SetBranchAddress("taus_muDecision", &taus_muDecision, &b_taus_muDecision);
   fChain->SetBranchAddress("taus_Nprongs", &taus_Nprongs, &b_taus_Nprongs);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("experimentType", &experimentType, &b_experimentType);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("model_params", &model_params, &b_model_params);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_px", &fastjets_AK4_R1p2_R0p5pT10_px, &b_fastjets_AK4_R1p2_R0p5pT10_px);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_py", &fastjets_AK4_R1p2_R0p5pT10_py, &b_fastjets_AK4_R1p2_R0p5pT10_py);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_pz", &fastjets_AK4_R1p2_R0p5pT10_pz, &b_fastjets_AK4_R1p2_R0p5pT10_pz);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_energy", &fastjets_AK4_R1p2_R0p5pT10_energy, &b_fastjets_AK4_R1p2_R0p5pT10_energy);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_phi", &fastjets_AK4_R1p2_R0p5pT10_phi, &b_fastjets_AK4_R1p2_R0p5pT10_phi);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_eta", &fastjets_AK4_R1p2_R0p5pT10_eta, &b_fastjets_AK4_R1p2_R0p5pT10_eta);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_index", &fastjets_AK4_R1p2_R0p5pT10_index, &b_fastjets_AK4_R1p2_R0p5pT10_index);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT10_nconstituents", &fastjets_AK4_R1p2_R0p5pT10_nconstituents, &b_fastjets_AK4_R1p2_R0p5pT10_nconstituents);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_px", &fastjets_AK4_R1p2_R0p5pT15_px, &b_fastjets_AK4_R1p2_R0p5pT15_px);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_py", &fastjets_AK4_R1p2_R0p5pT15_py, &b_fastjets_AK4_R1p2_R0p5pT15_py);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_pz", &fastjets_AK4_R1p2_R0p5pT15_pz, &b_fastjets_AK4_R1p2_R0p5pT15_pz);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_energy", &fastjets_AK4_R1p2_R0p5pT15_energy, &b_fastjets_AK4_R1p2_R0p5pT15_energy);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_phi", &fastjets_AK4_R1p2_R0p5pT15_phi, &b_fastjets_AK4_R1p2_R0p5pT15_phi);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_eta", &fastjets_AK4_R1p2_R0p5pT15_eta, &b_fastjets_AK4_R1p2_R0p5pT15_eta);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_index", &fastjets_AK4_R1p2_R0p5pT15_index, &b_fastjets_AK4_R1p2_R0p5pT15_index);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT15_nconstituents", &fastjets_AK4_R1p2_R0p5pT15_nconstituents, &b_fastjets_AK4_R1p2_R0p5pT15_nconstituents);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_px", &fastjets_AK4_R1p2_R0p5pT20_px, &b_fastjets_AK4_R1p2_R0p5pT20_px);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_py", &fastjets_AK4_R1p2_R0p5pT20_py, &b_fastjets_AK4_R1p2_R0p5pT20_py);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_pz", &fastjets_AK4_R1p2_R0p5pT20_pz, &b_fastjets_AK4_R1p2_R0p5pT20_pz);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_energy", &fastjets_AK4_R1p2_R0p5pT20_energy, &b_fastjets_AK4_R1p2_R0p5pT20_energy);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_phi", &fastjets_AK4_R1p2_R0p5pT20_phi, &b_fastjets_AK4_R1p2_R0p5pT20_phi);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_eta", &fastjets_AK4_R1p2_R0p5pT20_eta, &b_fastjets_AK4_R1p2_R0p5pT20_eta);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_index", &fastjets_AK4_R1p2_R0p5pT20_index, &b_fastjets_AK4_R1p2_R0p5pT20_index);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT20_nconstituents", &fastjets_AK4_R1p2_R0p5pT20_nconstituents, &b_fastjets_AK4_R1p2_R0p5pT20_nconstituents);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_px", &fastjets_AK4_R1p2_R0p5pT25_px, &b_fastjets_AK4_R1p2_R0p5pT25_px);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_py", &fastjets_AK4_R1p2_R0p5pT25_py, &b_fastjets_AK4_R1p2_R0p5pT25_py);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_pz", &fastjets_AK4_R1p2_R0p5pT25_pz, &b_fastjets_AK4_R1p2_R0p5pT25_pz);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_energy", &fastjets_AK4_R1p2_R0p5pT25_energy, &b_fastjets_AK4_R1p2_R0p5pT25_energy);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_phi", &fastjets_AK4_R1p2_R0p5pT25_phi, &b_fastjets_AK4_R1p2_R0p5pT25_phi);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_eta", &fastjets_AK4_R1p2_R0p5pT25_eta, &b_fastjets_AK4_R1p2_R0p5pT25_eta);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_index", &fastjets_AK4_R1p2_R0p5pT25_index, &b_fastjets_AK4_R1p2_R0p5pT25_index);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT25_nconstituents", &fastjets_AK4_R1p2_R0p5pT25_nconstituents, &b_fastjets_AK4_R1p2_R0p5pT25_nconstituents);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_px", &fastjets_AK4_R1p2_R0p5pT30_px, &b_fastjets_AK4_R1p2_R0p5pT30_px);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_py", &fastjets_AK4_R1p2_R0p5pT30_py, &b_fastjets_AK4_R1p2_R0p5pT30_py);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_pz", &fastjets_AK4_R1p2_R0p5pT30_pz, &b_fastjets_AK4_R1p2_R0p5pT30_pz);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_energy", &fastjets_AK4_R1p2_R0p5pT30_energy, &b_fastjets_AK4_R1p2_R0p5pT30_energy);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_phi", &fastjets_AK4_R1p2_R0p5pT30_phi, &b_fastjets_AK4_R1p2_R0p5pT30_phi);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_eta", &fastjets_AK4_R1p2_R0p5pT30_eta, &b_fastjets_AK4_R1p2_R0p5pT30_eta);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_index", &fastjets_AK4_R1p2_R0p5pT30_index, &b_fastjets_AK4_R1p2_R0p5pT30_index);
   fChain->SetBranchAddress("fastjets_AK4_R1p2_R0p5pT30_nconstituents", &fastjets_AK4_R1p2_R0p5pT30_nconstituents, &b_fastjets_AK4_R1p2_R0p5pT30_nconstituents);

}

