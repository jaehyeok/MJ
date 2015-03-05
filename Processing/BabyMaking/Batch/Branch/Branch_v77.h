#include <vector>
#include <string>
#include "TTree.h"
//#include "loader.C"
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
   Float_t         genHT;
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
   vector<float>   *isotk_pt;
   vector<float>   *isotk_phi;
   vector<float>   *isotk_eta;
   vector<float>   *isotk_iso;
   vector<float>   *isotk_dzpv;
   vector<int>     *isotk_charge;
   vector<int>     *taus_n_pfcands;
   vector<int>     *taus_decayMode;
   vector<float>   *taus_CombinedIsolationDeltaBetaCorrRaw3Hits;
   vector<bool>    *taus_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *taus_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *taus_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<bool>    *taus_byDecayModeFinding;
   vector<bool>    *taus_byDecayModeFindingNewDMs;
   vector<float>   *taus_chargedIsoPtSum;
   vector<float>   *taus_neutralIsoPtSum;
   vector<float>   *taus_puCorrPtSum;
   vector<bool>    *taus_againstMuonLoose3;
   vector<bool>    *taus_againstElectronLooseMVA5;
   vector<float>   *fjets30_pt;
   vector<float>   *fjets30_eta;
   vector<float>   *fjets30_phi;
   vector<float>   *fjets30_energy;
   vector<float>   *fjets30_m;
   Float_t         pfType1mets_uncert_JetEnUp_dpx;
   Float_t         pfType1mets_uncert_JetEnUp_dpy;
   Float_t         pfType1mets_uncert_JetEnUp_sumEt;
   Float_t         pfType1mets_uncert_JetEnDown_dpx;
   Float_t         pfType1mets_uncert_JetEnDown_dpy;
   Float_t         pfType1mets_uncert_JetEnDown_sumEt;
   Float_t         pfType1mets_uncert_JetResUp_dpx;
   Float_t         pfType1mets_uncert_JetResUp_dpy;
   Float_t         pfType1mets_uncert_JetResUp_sumEt;
   Float_t         pfType1mets_uncert_JetResDown_dpx;
   Float_t         pfType1mets_uncert_JetResDown_dpy;
   Float_t         pfType1mets_uncert_JetResDown_sumEt;

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
   TBranch        *b_genHT;   //!
   TBranch        *b_trackingfailurefilter_decision;   //!
   TBranch        *b_goodVerticesfilter_decision;   //!
   TBranch        *b_cschalofilter_decision;   //!
   TBranch        *b_trkPOGfilter_decision;   //!
   TBranch        *b_trkPOG_logErrorTooManyClustersfilter_decision;   //!
   TBranch        *b_ecalDeadCellTriggerPrimitivefilter_decision;   //!
   TBranch        *b_ecallaserfilter_decision;   //!
   TBranch        *b_trkPOG_manystripclus53Xfilter_decision;   //!
   TBranch        *b_eebadscfilter_decision;   //!
   TBranch        *b_METFiltersfilter_decision;   //!
   TBranch        *b_HBHENoisefilter_decision;   //!
   TBranch        *b_trkPOG_toomanystripclus53Xfilter_decision;   //!
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
   TBranch        *b_isotk_pt;   //!
   TBranch        *b_isotk_phi;   //!
   TBranch        *b_isotk_eta;   //!
   TBranch        *b_isotk_iso;   //!
   TBranch        *b_isotk_dzpv;   //!
   TBranch        *b_isotk_charge;   //!
   TBranch        *b_taus_n_pfcands;   //!
   TBranch        *b_taus_decayMode;   //!
   TBranch        *b_taus_CombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_taus_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_taus_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_taus_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_taus_byDecayModeFinding;   //!
   TBranch        *b_taus_byDecayModeFindingNewDMs;   //!
   TBranch        *b_taus_chargedIsoPtSum;   //!
   TBranch        *b_taus_neutralIsoPtSum;   //!
   TBranch        *b_taus_puCorrPtSum;   //!
   TBranch        *b_taus_againstMuonLoose3;   //!
   TBranch        *b_taus_againstElectronLooseMVA5;   //!
   TBranch        *b_fjets30_pt;   //!
   TBranch        *b_fjets30_eta;   //!
   TBranch        *b_fjets30_phi;   //!
   TBranch        *b_fjets30_energy;   //!
   TBranch        *b_fjets30_m;   //!
   TBranch        *b_pfType1mets_uncert_JetEnUp_dpx;   //!
   TBranch        *b_pfType1mets_uncert_JetEnUp_dpy;   //!
   TBranch        *b_pfType1mets_uncert_JetEnUp_sumEt;   //!
   TBranch        *b_pfType1mets_uncert_JetEnDown_dpx;   //!
   TBranch        *b_pfType1mets_uncert_JetEnDown_dpy;   //!
   TBranch        *b_pfType1mets_uncert_JetEnDown_sumEt;   //!
   TBranch        *b_pfType1mets_uncert_JetResUp_dpx;   //!
   TBranch        *b_pfType1mets_uncert_JetResUp_dpy;   //!
   TBranch        *b_pfType1mets_uncert_JetResUp_sumEt;   //!
   TBranch        *b_pfType1mets_uncert_JetResDown_dpx;   //!
   TBranch        *b_pfType1mets_uncert_JetResDown_dpy;   //!
   TBranch        *b_pfType1mets_uncert_JetResDown_sumEt;   //!

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
   vector<float>   *els_expectedMissingInnerHits;
   vector<float>   *els_tightId;
   vector<float>   *els_looseId;
   vector<float>   *els_robustTightId;
   vector<float>   *els_robustLooseId;
   vector<float>   *els_robustHighEnergyId;
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
   vector<float>   *jets_AK4_energy;
   vector<float>   *jets_AK4_et;
   vector<float>   *jets_AK4_eta;
   vector<float>   *jets_AK4_phi;
   vector<float>   *jets_AK4_pt;
   vector<float>   *jets_AK4_px;
   vector<float>   *jets_AK4_py;
   vector<float>   *jets_AK4_pz;
   vector<float>   *jets_AK4_status;
   vector<float>   *jets_AK4_theta;
   vector<float>   *jets_AK4_btag_jetBProb;
   vector<float>   *jets_AK4_btag_jetProb;
   vector<float>   *jets_AK4_btag_TC_highPur;
   vector<float>   *jets_AK4_btag_TC_highEff;
   vector<float>   *jets_AK4_btag_secVertexHighEff;
   vector<float>   *jets_AK4_btag_secVertexHighPur;
   vector<float>   *jets_AK4_btag_inc_secVertexCombined;
   vector<float>   *jets_AK4_btag_pf_secVertexCombined;
   vector<float>   *jets_AK4_btag_MVA;
   vector<float>   *jets_AK4_pileupID_MVA;
   vector<float>   *jets_AK4_parton_Id;
   vector<float>   *jets_AK4_parton_motherId;
   vector<float>   *jets_AK4_parton_grandmotherID;
   vector<float>   *jets_AK4_parton_pt;
   vector<float>   *jets_AK4_parton_phi;
   vector<float>   *jets_AK4_parton_eta;
   vector<float>   *jets_AK4_parton_Energy;
   vector<float>   *jets_AK4_parton_mass;
   vector<float>   *jets_AK4_partonFlavour;
   vector<float>   *jets_AK4_gen_pt;
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
   UInt_t          Nmc_final;
   vector<float>   *mc_final_id;
   vector<float>   *mc_final_pt;
   vector<float>   *mc_final_eta;
   vector<float>   *mc_final_phi;
   vector<float>   *mc_final_energy;
   vector<float>   *mc_final_charge;
   vector<float>   *mc_final_mother_id;
   vector<float>   *mc_final_grandmother_id;
   vector<float>   *mc_final_ggrandmother_id;
   vector<float>   *mc_final_numOfMothers;
   UInt_t          Nmc_jets;
   vector<float>   *mc_jets_pt;
   vector<float>   *mc_jets_eta;
   vector<float>   *mc_jets_phi;
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
   vector<float>   *mus_globalTrack_normalizedChi2;
   vector<float>   *mus_trkPositionMatch;
   vector<float>   *mus_trkKink;
   vector<float>   *mus_tkHits;
   vector<float>   *mus_tkHitsFrac;
   vector<float>   *mus_segmentCompatibility;
   vector<float>   *mus_caloCompatibility;
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
   vector<float>   *mus_isGlobalMuon;
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
   vector<float>   *mus_cm_LayersWithMeasurement;
   vector<float>   *mus_cm_PixelLayersWithMeasurement;
   vector<float>   *mus_cm_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *mus_cm_LayersWithoutMeasurement;
   vector<float>   *mus_dB;
   vector<float>   *mus_numberOfMatchedStations;
   UInt_t          NpfType1mets;
   vector<float>   *pfType1mets_et;
   vector<float>   *pfType1mets_phi;
   vector<float>   *pfType1mets_ex;
   vector<float>   *pfType1mets_ey;
   vector<float>   *pfType1mets_NeutralEMFraction;
   vector<float>   *pfType1mets_NeutralHadEtFraction;
   vector<float>   *pfType1mets_ChargedEMEtFraction;
   vector<float>   *pfType1mets_ChargedHadEtFraction;
   vector<float>   *pfType1mets_MuonEtFraction;
   vector<float>   *pfType1mets_Type6EtFraction;
   vector<float>   *pfType1mets_Type7EtFraction;
   vector<float>   *pfType1mets_sumEt;
   vector<float>   *pfType1mets_unCPhi;
   vector<float>   *pfType1mets_unCPt;
   vector<float>   *pfType1mets_gen_et;
   vector<float>   *pfType1mets_gen_phi;
   UInt_t          Npfcand;
   vector<float>   *pfcand_pdgId;
   vector<float>   *pfcand_pt;
   vector<float>   *pfcand_eta;
   vector<float>   *pfcand_phi;
   vector<float>   *pfcand_energy;
   vector<float>   *pfcand_charge;
   vector<float>   *pfcand_dz;
   vector<float>   *pfcand_dxy;
   vector<float>   *pfcand_fromPV;
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
   vector<float>   *photons_passElectronVeto;
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
   UInt_t          Ntaus;
   vector<float>   *taus_energy;
   vector<float>   *taus_et;
   vector<float>   *taus_eta;
   vector<float>   *taus_phi;
   vector<float>   *taus_pt;
   vector<float>   *taus_px;
   vector<float>   *taus_py;
   vector<float>   *taus_pz;
   vector<float>   *taus_status;
   vector<float>   *taus_theta;
   vector<float>   *taus_charge;
   vector<float>   *taus_leadChargedHadrCand_pt;
   vector<float>   *taus_leadChargedHadrCand_charge;
   vector<float>   *taus_leadChargedHadrCand_eta;
   vector<float>   *taus_leadChargedHadrCand_phi;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          experimentType;
   UInt_t          bunchCrossing;
   UInt_t          orbitNumber;
   Float_t         weight;
   vector<int>     *weightIndex;
   vector<float>   *weightVector;
   string          *model_params;

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
   TBranch        *b_els_expectedMissingInnerHits;   //!
   TBranch        *b_els_tightId;   //!
   TBranch        *b_els_looseId;   //!
   TBranch        *b_els_robustTightId;   //!
   TBranch        *b_els_robustLooseId;   //!
   TBranch        *b_els_robustHighEnergyId;   //!
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
   TBranch        *b_jets_AK4_energy;   //!
   TBranch        *b_jets_AK4_et;   //!
   TBranch        *b_jets_AK4_eta;   //!
   TBranch        *b_jets_AK4_phi;   //!
   TBranch        *b_jets_AK4_pt;   //!
   TBranch        *b_jets_AK4_px;   //!
   TBranch        *b_jets_AK4_py;   //!
   TBranch        *b_jets_AK4_pz;   //!
   TBranch        *b_jets_AK4_status;   //!
   TBranch        *b_jets_AK4_theta;   //!
   TBranch        *b_jets_AK4_btag_jetBProb;   //!
   TBranch        *b_jets_AK4_btag_jetProb;   //!
   TBranch        *b_jets_AK4_btag_TC_highPur;   //!
   TBranch        *b_jets_AK4_btag_TC_highEff;   //!
   TBranch        *b_jets_AK4_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK4_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK4_btag_inc_secVertexCombined;   //!
   TBranch        *b_jets_AK4_btag_pf_secVertexCombined;   //!
   TBranch        *b_jets_AK4_btag_MVA;   //!
   TBranch        *b_jets_AK4_pileupID_MVA;   //!
   TBranch        *b_jets_AK4_parton_Id;   //!
   TBranch        *b_jets_AK4_parton_motherId;   //!
   TBranch        *b_jets_AK4_parton_grandmotherID;   //!
   TBranch        *b_jets_AK4_parton_pt;   //!
   TBranch        *b_jets_AK4_parton_phi;   //!
   TBranch        *b_jets_AK4_parton_eta;   //!
   TBranch        *b_jets_AK4_parton_Energy;   //!
   TBranch        *b_jets_AK4_parton_mass;   //!
   TBranch        *b_jets_AK4_partonFlavour;   //!
   TBranch        *b_jets_AK4_gen_pt;   //!
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
   TBranch        *b_Nmc_final;   //!
   TBranch        *b_mc_final_id;   //!
   TBranch        *b_mc_final_pt;   //!
   TBranch        *b_mc_final_eta;   //!
   TBranch        *b_mc_final_phi;   //!
   TBranch        *b_mc_final_energy;   //!
   TBranch        *b_mc_final_charge;   //!
   TBranch        *b_mc_final_mother_id;   //!
   TBranch        *b_mc_final_grandmother_id;   //!
   TBranch        *b_mc_final_ggrandmother_id;   //!
   TBranch        *b_mc_final_numOfMothers;   //!
   TBranch        *b_Nmc_jets;   //!
   TBranch        *b_mc_jets_pt;   //!
   TBranch        *b_mc_jets_eta;   //!
   TBranch        *b_mc_jets_phi;   //!
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
   TBranch        *b_mus_globalTrack_normalizedChi2;   //!
   TBranch        *b_mus_trkPositionMatch;   //!
   TBranch        *b_mus_trkKink;   //!
   TBranch        *b_mus_tkHits;   //!
   TBranch        *b_mus_tkHitsFrac;   //!
   TBranch        *b_mus_segmentCompatibility;   //!
   TBranch        *b_mus_caloCompatibility;   //!
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
   TBranch        *b_mus_isGlobalMuon;   //!
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
   TBranch        *b_mus_cm_LayersWithMeasurement;   //!
   TBranch        *b_mus_cm_PixelLayersWithMeasurement;   //!
   TBranch        *b_mus_cm_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_mus_cm_LayersWithoutMeasurement;   //!
   TBranch        *b_mus_dB;   //!
   TBranch        *b_mus_numberOfMatchedStations;   //!
   TBranch        *b_NpfType1mets;   //!
   TBranch        *b_pfType1mets_et;   //!
   TBranch        *b_pfType1mets_phi;   //!
   TBranch        *b_pfType1mets_ex;   //!
   TBranch        *b_pfType1mets_ey;   //!
   TBranch        *b_pfType1mets_NeutralEMFraction;   //!
   TBranch        *b_pfType1mets_NeutralHadEtFraction;   //!
   TBranch        *b_pfType1mets_ChargedEMEtFraction;   //!
   TBranch        *b_pfType1mets_ChargedHadEtFraction;   //!
   TBranch        *b_pfType1mets_MuonEtFraction;   //!
   TBranch        *b_pfType1mets_Type6EtFraction;   //!
   TBranch        *b_pfType1mets_Type7EtFraction;   //!
   TBranch        *b_pfType1mets_sumEt;   //!
   TBranch        *b_pfType1mets_unCPhi;   //!
   TBranch        *b_pfType1mets_unCPt;   //!
   TBranch        *b_pfType1mets_gen_et;   //!
   TBranch        *b_pfType1mets_gen_phi;   //!
   TBranch        *b_Npfcand;   //!
   TBranch        *b_pfcand_pdgId;   //!
   TBranch        *b_pfcand_pt;   //!
   TBranch        *b_pfcand_eta;   //!
   TBranch        *b_pfcand_phi;   //!
   TBranch        *b_pfcand_energy;   //!
   TBranch        *b_pfcand_charge;   //!
   TBranch        *b_pfcand_dz;   //!
   TBranch        *b_pfcand_dxy;   //!
   TBranch        *b_pfcand_fromPV;   //!
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
   TBranch        *b_photons_passElectronVeto;   //!
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
   TBranch        *b_Ntaus;   //!
   TBranch        *b_taus_energy;   //!
   TBranch        *b_taus_et;   //!
   TBranch        *b_taus_eta;   //!
   TBranch        *b_taus_phi;   //!
   TBranch        *b_taus_pt;   //!
   TBranch        *b_taus_px;   //!
   TBranch        *b_taus_py;   //!
   TBranch        *b_taus_pz;   //!
   TBranch        *b_taus_status;   //!
   TBranch        *b_taus_theta;   //!
   TBranch        *b_taus_charge;   //!
   TBranch        *b_taus_leadChargedHadrCand_pt;   //!
   TBranch        *b_taus_leadChargedHadrCand_charge;   //!
   TBranch        *b_taus_leadChargedHadrCand_eta;   //!
   TBranch        *b_taus_leadChargedHadrCand_phi;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_experimentType;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_weightIndex;   //!
   TBranch        *b_weightVector;   //!
   TBranch        *b_model_params;   //!

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
   isotk_pt = 0;
   isotk_phi = 0;
   isotk_eta = 0;
   isotk_iso = 0;
   isotk_dzpv = 0;
   isotk_charge = 0;
   taus_n_pfcands = 0;
   taus_decayMode = 0;
   taus_CombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   taus_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   taus_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   taus_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   taus_byDecayModeFinding = 0;
   taus_byDecayModeFindingNewDMs = 0;
   taus_chargedIsoPtSum = 0;
   taus_neutralIsoPtSum = 0;
   taus_puCorrPtSum = 0;
   taus_againstMuonLoose3 = 0;
   taus_againstElectronLooseMVA5 = 0;
   fjets30_pt = 0;
   fjets30_eta = 0;
   fjets30_phi = 0;
   fjets30_energy = 0;
   fjets30_m = 0;

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
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("trackingfailurefilter_decision", &trackingfailurefilter_decision, &b_trackingfailurefilter_decision);
   fChain->SetBranchAddress("goodVerticesfilter_decision", &goodVerticesfilter_decision, &b_goodVerticesfilter_decision);
   fChain->SetBranchAddress("cschalofilter_decision", &cschalofilter_decision, &b_cschalofilter_decision);
   fChain->SetBranchAddress("trkPOGfilter_decision", &trkPOGfilter_decision, &b_trkPOGfilter_decision);
   fChain->SetBranchAddress("trkPOG_logErrorTooManyClustersfilter_decision", &trkPOG_logErrorTooManyClustersfilter_decision, &b_trkPOG_logErrorTooManyClustersfilter_decision);
   fChain->SetBranchAddress("EcalDeadCellTriggerPrimitivefilter_decision", &EcalDeadCellTriggerPrimitivefilter_decision, &b_ecalDeadCellTriggerPrimitivefilter_decision);
   fChain->SetBranchAddress("ecallaserfilter_decision", &ecallaserfilter_decision, &b_ecallaserfilter_decision);
   fChain->SetBranchAddress("trkPOG_manystripclus53Xfilter_decision", &trkPOG_manystripclus53Xfilter_decision, &b_trkPOG_manystripclus53Xfilter_decision);
   fChain->SetBranchAddress("eebadscfilter_decision", &eebadscfilter_decision, &b_eebadscfilter_decision);
   fChain->SetBranchAddress("METFiltersfilter_decision", &METFiltersfilter_decision, &b_METFiltersfilter_decision);
   fChain->SetBranchAddress("HBHENoisefilter_decision", &HBHENoisefilter_decision, &b_HBHENoisefilter_decision);
   fChain->SetBranchAddress("trkPOG_toomanystripclus53Xfilter_decision", &trkPOG_toomanystripclus53Xfilter_decision, &b_trkPOG_toomanystripclus53Xfilter_decision);
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
   fChain->SetBranchAddress("isotk_pt", &isotk_pt, &b_isotk_pt);
   fChain->SetBranchAddress("isotk_phi", &isotk_phi, &b_isotk_phi);
   fChain->SetBranchAddress("isotk_eta", &isotk_eta, &b_isotk_eta);
   fChain->SetBranchAddress("isotk_iso", &isotk_iso, &b_isotk_iso);
   fChain->SetBranchAddress("isotk_dzpv", &isotk_dzpv, &b_isotk_dzpv);
   fChain->SetBranchAddress("isotk_charge", &isotk_charge, &b_isotk_charge);
   fChain->SetBranchAddress("taus_n_pfcands", &taus_n_pfcands, &b_taus_n_pfcands);
   fChain->SetBranchAddress("taus_decayMode", &taus_decayMode, &b_taus_decayMode);
   fChain->SetBranchAddress("taus_CombinedIsolationDeltaBetaCorrRaw3Hits", &taus_CombinedIsolationDeltaBetaCorrRaw3Hits, &b_taus_CombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("taus_byLooseCombinedIsolationDeltaBetaCorr3Hits", &taus_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_taus_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("taus_byMediumCombinedIsolationDeltaBetaCorr3Hits", &taus_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_taus_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("taus_byTightCombinedIsolationDeltaBetaCorr3Hits", &taus_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_taus_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("taus_byDecayModeFinding", &taus_byDecayModeFinding, &b_taus_byDecayModeFinding);
   fChain->SetBranchAddress("taus_byDecayModeFindingNewDMs", &taus_byDecayModeFindingNewDMs, &b_taus_byDecayModeFindingNewDMs);
   fChain->SetBranchAddress("taus_chargedIsoPtSum", &taus_chargedIsoPtSum, &b_taus_chargedIsoPtSum);
   fChain->SetBranchAddress("taus_neutralIsoPtSum", &taus_neutralIsoPtSum, &b_taus_neutralIsoPtSum);
   fChain->SetBranchAddress("taus_puCorrPtSum", &taus_puCorrPtSum, &b_taus_puCorrPtSum);
   fChain->SetBranchAddress("taus_againstMuonLoose3", &taus_againstMuonLoose3, &b_taus_againstMuonLoose3);
   fChain->SetBranchAddress("taus_againstElectronLooseMVA5", &taus_againstElectronLooseMVA5, &b_taus_againstElectronLooseMVA5);
   fChain->SetBranchAddress("fjets30_pt", &fjets30_pt, &b_fjets30_pt);
   fChain->SetBranchAddress("fjets30_eta", &fjets30_eta, &b_fjets30_eta);
   fChain->SetBranchAddress("fjets30_phi", &fjets30_phi, &b_fjets30_phi);
   fChain->SetBranchAddress("fjets30_energy", &fjets30_energy, &b_fjets30_energy);
   fChain->SetBranchAddress("fjets30_m", &fjets30_m, &b_fjets30_m);
   fChain->SetBranchAddress("pfType1mets_uncert_JetEnUp_dpx", &pfType1mets_uncert_JetEnUp_dpx, &b_pfType1mets_uncert_JetEnUp_dpx);
   fChain->SetBranchAddress("pfType1mets_uncert_JetEnUp_dpy", &pfType1mets_uncert_JetEnUp_dpy, &b_pfType1mets_uncert_JetEnUp_dpy);
   fChain->SetBranchAddress("pfType1mets_uncert_JetEnUp_sumEt", &pfType1mets_uncert_JetEnUp_sumEt, &b_pfType1mets_uncert_JetEnUp_sumEt);
   fChain->SetBranchAddress("pfType1mets_uncert_JetEnDown_dpx", &pfType1mets_uncert_JetEnDown_dpx, &b_pfType1mets_uncert_JetEnDown_dpx);
   fChain->SetBranchAddress("pfType1mets_uncert_JetEnDown_dpy", &pfType1mets_uncert_JetEnDown_dpy, &b_pfType1mets_uncert_JetEnDown_dpy);
   fChain->SetBranchAddress("pfType1mets_uncert_JetEnDown_sumEt", &pfType1mets_uncert_JetEnDown_sumEt, &b_pfType1mets_uncert_JetEnDown_sumEt);
   fChain->SetBranchAddress("pfType1mets_uncert_JetResUp_dpx", &pfType1mets_uncert_JetResUp_dpx, &b_pfType1mets_uncert_JetResUp_dpx);
   fChain->SetBranchAddress("pfType1mets_uncert_JetResUp_dpy", &pfType1mets_uncert_JetResUp_dpy, &b_pfType1mets_uncert_JetResUp_dpy);
   fChain->SetBranchAddress("pfType1mets_uncert_JetResUp_sumEt", &pfType1mets_uncert_JetResUp_sumEt, &b_pfType1mets_uncert_JetResUp_sumEt);
   fChain->SetBranchAddress("pfType1mets_uncert_JetResDown_dpx", &pfType1mets_uncert_JetResDown_dpx, &b_pfType1mets_uncert_JetResDown_dpx);
   fChain->SetBranchAddress("pfType1mets_uncert_JetResDown_dpy", &pfType1mets_uncert_JetResDown_dpy, &b_pfType1mets_uncert_JetResDown_dpy);
   fChain->SetBranchAddress("pfType1mets_uncert_JetResDown_sumEt", &pfType1mets_uncert_JetResDown_sumEt, &b_pfType1mets_uncert_JetResDown_sumEt);
   
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
   els_expectedMissingInnerHits = 0;
   els_tightId = 0;
   els_looseId = 0;
   els_robustTightId = 0;
   els_robustLooseId = 0;
   els_robustHighEnergyId = 0;
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
   jets_AK4_energy = 0;
   jets_AK4_et = 0;
   jets_AK4_eta = 0;
   jets_AK4_phi = 0;
   jets_AK4_pt = 0;
   jets_AK4_px = 0;
   jets_AK4_py = 0;
   jets_AK4_pz = 0;
   jets_AK4_status = 0;
   jets_AK4_theta = 0;
   jets_AK4_btag_jetBProb = 0;
   jets_AK4_btag_jetProb = 0;
   jets_AK4_btag_TC_highPur = 0;
   jets_AK4_btag_TC_highEff = 0;
   jets_AK4_btag_secVertexHighEff = 0;
   jets_AK4_btag_secVertexHighPur = 0;
   jets_AK4_btag_inc_secVertexCombined = 0;
   jets_AK4_btag_pf_secVertexCombined = 0;
   jets_AK4_btag_MVA = 0;
   jets_AK4_pileupID_MVA = 0;
   jets_AK4_parton_Id = 0;
   jets_AK4_parton_motherId = 0;
   jets_AK4_parton_grandmotherID = 0;
   jets_AK4_parton_pt = 0;
   jets_AK4_parton_phi = 0;
   jets_AK4_parton_eta = 0;
   jets_AK4_parton_Energy = 0;
   jets_AK4_parton_mass = 0;
   jets_AK4_partonFlavour = 0;
   jets_AK4_gen_pt = 0;
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
   mc_final_id = 0;
   mc_final_pt = 0;
   mc_final_eta = 0;
   mc_final_phi = 0;
   mc_final_energy = 0;
   mc_final_charge = 0;
   mc_final_mother_id = 0;
   mc_final_grandmother_id = 0;
   mc_final_ggrandmother_id = 0;
   mc_final_numOfMothers = 0;
   mc_jets_pt = 0;
   mc_jets_eta = 0;
   mc_jets_phi = 0;
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
   mus_globalTrack_normalizedChi2 = 0;
   mus_trkPositionMatch = 0;
   mus_trkKink = 0;
   mus_tkHits = 0;
   mus_tkHitsFrac = 0;
   mus_segmentCompatibility = 0;
   mus_caloCompatibility = 0;
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
   mus_isGlobalMuon = 0;
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
   mus_cm_LayersWithMeasurement = 0;
   mus_cm_PixelLayersWithMeasurement = 0;
   mus_cm_ValidStripLayersWithMonoAndStereoHit = 0;
   mus_cm_LayersWithoutMeasurement = 0;
   mus_dB = 0;
   mus_numberOfMatchedStations = 0;
   pfType1mets_et = 0;
   pfType1mets_phi = 0;
   pfType1mets_ex = 0;
   pfType1mets_ey = 0;
   pfType1mets_NeutralEMFraction = 0;
   pfType1mets_NeutralHadEtFraction = 0;
   pfType1mets_ChargedEMEtFraction = 0;
   pfType1mets_ChargedHadEtFraction = 0;
   pfType1mets_MuonEtFraction = 0;
   pfType1mets_Type6EtFraction = 0;
   pfType1mets_Type7EtFraction = 0;
   pfType1mets_sumEt = 0;
   pfType1mets_unCPhi = 0;
   pfType1mets_unCPt = 0;
   pfType1mets_gen_et = 0;
   pfType1mets_gen_phi = 0;
   pfcand_pdgId = 0;
   pfcand_pt = 0;
   pfcand_eta = 0;
   pfcand_phi = 0;
   pfcand_energy = 0;
   pfcand_charge = 0;
   pfcand_dz = 0;
   pfcand_dxy = 0;
   pfcand_fromPV = 0;
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
   photons_passElectronVeto = 0;
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
   taus_energy = 0;
   taus_et = 0;
   taus_eta = 0;
   taus_phi = 0;
   taus_pt = 0;
   taus_px = 0;
   taus_py = 0;
   taus_pz = 0;
   taus_status = 0;
   taus_theta = 0;
   taus_charge = 0;
   taus_leadChargedHadrCand_pt = 0;
   taus_leadChargedHadrCand_charge = 0;
   taus_leadChargedHadrCand_eta = 0;
   taus_leadChargedHadrCand_phi = 0;
   weightIndex = 0;
   weightVector = 0;
   model_params = 0;

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
   fChain->SetBranchAddress("els_expectedMissingInnerHits", &els_expectedMissingInnerHits, &b_els_expectedMissingInnerHits);
   fChain->SetBranchAddress("els_tightId", &els_tightId, &b_els_tightId);
   fChain->SetBranchAddress("els_looseId", &els_looseId, &b_els_looseId);
   fChain->SetBranchAddress("els_robustTightId", &els_robustTightId, &b_els_robustTightId);
   fChain->SetBranchAddress("els_robustLooseId", &els_robustLooseId, &b_els_robustLooseId);
   fChain->SetBranchAddress("els_robustHighEnergyId", &els_robustHighEnergyId, &b_els_robustHighEnergyId);
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
   fChain->SetBranchAddress("jets_AK4_energy", &jets_AK4_energy, &b_jets_AK4_energy);
   fChain->SetBranchAddress("jets_AK4_et", &jets_AK4_et, &b_jets_AK4_et);
   fChain->SetBranchAddress("jets_AK4_eta", &jets_AK4_eta, &b_jets_AK4_eta);
   fChain->SetBranchAddress("jets_AK4_phi", &jets_AK4_phi, &b_jets_AK4_phi);
   fChain->SetBranchAddress("jets_AK4_pt", &jets_AK4_pt, &b_jets_AK4_pt);
   fChain->SetBranchAddress("jets_AK4_px", &jets_AK4_px, &b_jets_AK4_px);
   fChain->SetBranchAddress("jets_AK4_py", &jets_AK4_py, &b_jets_AK4_py);
   fChain->SetBranchAddress("jets_AK4_pz", &jets_AK4_pz, &b_jets_AK4_pz);
   fChain->SetBranchAddress("jets_AK4_status", &jets_AK4_status, &b_jets_AK4_status);
   fChain->SetBranchAddress("jets_AK4_theta", &jets_AK4_theta, &b_jets_AK4_theta);
   fChain->SetBranchAddress("jets_AK4_btag_jetBProb", &jets_AK4_btag_jetBProb, &b_jets_AK4_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK4_btag_jetProb", &jets_AK4_btag_jetProb, &b_jets_AK4_btag_jetProb);
   fChain->SetBranchAddress("jets_AK4_btag_TC_highPur", &jets_AK4_btag_TC_highPur, &b_jets_AK4_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK4_btag_TC_highEff", &jets_AK4_btag_TC_highEff, &b_jets_AK4_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK4_btag_secVertexHighEff", &jets_AK4_btag_secVertexHighEff, &b_jets_AK4_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK4_btag_secVertexHighPur", &jets_AK4_btag_secVertexHighPur, &b_jets_AK4_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK4_btag_inc_secVertexCombined", &jets_AK4_btag_inc_secVertexCombined, &b_jets_AK4_btag_inc_secVertexCombined);
   fChain->SetBranchAddress("jets_AK4_btag_pf_secVertexCombined", &jets_AK4_btag_pf_secVertexCombined, &b_jets_AK4_btag_pf_secVertexCombined);
   fChain->SetBranchAddress("jets_AK4_btag_MVA", &jets_AK4_btag_MVA, &b_jets_AK4_btag_MVA);
   fChain->SetBranchAddress("jets_AK4_pileupID_MVA", &jets_AK4_pileupID_MVA, &b_jets_AK4_pileupID_MVA);
   fChain->SetBranchAddress("jets_AK4_parton_Id", &jets_AK4_parton_Id, &b_jets_AK4_parton_Id);
   fChain->SetBranchAddress("jets_AK4_parton_motherId", &jets_AK4_parton_motherId, &b_jets_AK4_parton_motherId);
   fChain->SetBranchAddress("jets_AK4_parton_grandmotherID", &jets_AK4_parton_grandmotherID, &b_jets_AK4_parton_grandmotherID);
   fChain->SetBranchAddress("jets_AK4_parton_pt", &jets_AK4_parton_pt, &b_jets_AK4_parton_pt);
   fChain->SetBranchAddress("jets_AK4_parton_phi", &jets_AK4_parton_phi, &b_jets_AK4_parton_phi);
   fChain->SetBranchAddress("jets_AK4_parton_eta", &jets_AK4_parton_eta, &b_jets_AK4_parton_eta);
   fChain->SetBranchAddress("jets_AK4_parton_Energy", &jets_AK4_parton_Energy, &b_jets_AK4_parton_Energy);
   fChain->SetBranchAddress("jets_AK4_parton_mass", &jets_AK4_parton_mass, &b_jets_AK4_parton_mass);
   fChain->SetBranchAddress("jets_AK4_partonFlavour", &jets_AK4_partonFlavour, &b_jets_AK4_partonFlavour);
   fChain->SetBranchAddress("jets_AK4_gen_pt", &jets_AK4_gen_pt, &b_jets_AK4_gen_pt);
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
   fChain->SetBranchAddress("Nmc_final", &Nmc_final, &b_Nmc_final);
   fChain->SetBranchAddress("mc_final_id", &mc_final_id, &b_mc_final_id);
   fChain->SetBranchAddress("mc_final_pt", &mc_final_pt, &b_mc_final_pt);
   fChain->SetBranchAddress("mc_final_eta", &mc_final_eta, &b_mc_final_eta);
   fChain->SetBranchAddress("mc_final_phi", &mc_final_phi, &b_mc_final_phi);
   fChain->SetBranchAddress("mc_final_energy", &mc_final_energy, &b_mc_final_energy);
   fChain->SetBranchAddress("mc_final_charge", &mc_final_charge, &b_mc_final_charge);
   fChain->SetBranchAddress("mc_final_mother_id", &mc_final_mother_id, &b_mc_final_mother_id);
   fChain->SetBranchAddress("mc_final_grandmother_id", &mc_final_grandmother_id, &b_mc_final_grandmother_id);
   fChain->SetBranchAddress("mc_final_ggrandmother_id", &mc_final_ggrandmother_id, &b_mc_final_ggrandmother_id);
   fChain->SetBranchAddress("mc_final_numOfMothers", &mc_final_numOfMothers, &b_mc_final_numOfMothers);
   fChain->SetBranchAddress("Nmc_jets", &Nmc_jets, &b_Nmc_jets);
   fChain->SetBranchAddress("mc_jets_pt", &mc_jets_pt, &b_mc_jets_pt);
   fChain->SetBranchAddress("mc_jets_eta", &mc_jets_eta, &b_mc_jets_eta);
   fChain->SetBranchAddress("mc_jets_phi", &mc_jets_phi, &b_mc_jets_phi);
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
   fChain->SetBranchAddress("mus_globalTrack_normalizedChi2", &mus_globalTrack_normalizedChi2, &b_mus_globalTrack_normalizedChi2);
   fChain->SetBranchAddress("mus_trkPositionMatch", &mus_trkPositionMatch, &b_mus_trkPositionMatch);
   fChain->SetBranchAddress("mus_trkKink", &mus_trkKink, &b_mus_trkKink);
   fChain->SetBranchAddress("mus_tkHits", &mus_tkHits, &b_mus_tkHits);
   fChain->SetBranchAddress("mus_tkHitsFrac", &mus_tkHitsFrac, &b_mus_tkHitsFrac);
   fChain->SetBranchAddress("mus_segmentCompatibility", &mus_segmentCompatibility, &b_mus_segmentCompatibility);
   fChain->SetBranchAddress("mus_caloCompatibility", &mus_caloCompatibility, &b_mus_caloCompatibility);
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
   fChain->SetBranchAddress("mus_isGlobalMuon", &mus_isGlobalMuon, &b_mus_isGlobalMuon);
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
   fChain->SetBranchAddress("mus_cm_LayersWithMeasurement", &mus_cm_LayersWithMeasurement, &b_mus_cm_LayersWithMeasurement);
   fChain->SetBranchAddress("mus_cm_PixelLayersWithMeasurement", &mus_cm_PixelLayersWithMeasurement, &b_mus_cm_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("mus_cm_ValidStripLayersWithMonoAndStereoHit", &mus_cm_ValidStripLayersWithMonoAndStereoHit, &b_mus_cm_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("mus_cm_LayersWithoutMeasurement", &mus_cm_LayersWithoutMeasurement, &b_mus_cm_LayersWithoutMeasurement);
   fChain->SetBranchAddress("mus_dB", &mus_dB, &b_mus_dB);
   fChain->SetBranchAddress("mus_numberOfMatchedStations", &mus_numberOfMatchedStations, &b_mus_numberOfMatchedStations);
   fChain->SetBranchAddress("NpfType1mets", &NpfType1mets, &b_NpfType1mets);
   fChain->SetBranchAddress("pfType1mets_et", &pfType1mets_et, &b_pfType1mets_et);
   fChain->SetBranchAddress("pfType1mets_phi", &pfType1mets_phi, &b_pfType1mets_phi);
   fChain->SetBranchAddress("pfType1mets_ex", &pfType1mets_ex, &b_pfType1mets_ex);
   fChain->SetBranchAddress("pfType1mets_ey", &pfType1mets_ey, &b_pfType1mets_ey);
   fChain->SetBranchAddress("pfType1mets_NeutralEMFraction", &pfType1mets_NeutralEMFraction, &b_pfType1mets_NeutralEMFraction);
   fChain->SetBranchAddress("pfType1mets_NeutralHadEtFraction", &pfType1mets_NeutralHadEtFraction, &b_pfType1mets_NeutralHadEtFraction);
   fChain->SetBranchAddress("pfType1mets_ChargedEMEtFraction", &pfType1mets_ChargedEMEtFraction, &b_pfType1mets_ChargedEMEtFraction);
   fChain->SetBranchAddress("pfType1mets_ChargedHadEtFraction", &pfType1mets_ChargedHadEtFraction, &b_pfType1mets_ChargedHadEtFraction);
   fChain->SetBranchAddress("pfType1mets_MuonEtFraction", &pfType1mets_MuonEtFraction, &b_pfType1mets_MuonEtFraction);
   fChain->SetBranchAddress("pfType1mets_Type6EtFraction", &pfType1mets_Type6EtFraction, &b_pfType1mets_Type6EtFraction);
   fChain->SetBranchAddress("pfType1mets_Type7EtFraction", &pfType1mets_Type7EtFraction, &b_pfType1mets_Type7EtFraction);
   fChain->SetBranchAddress("pfType1mets_sumEt", &pfType1mets_sumEt, &b_pfType1mets_sumEt);
   fChain->SetBranchAddress("pfType1mets_unCPhi", &pfType1mets_unCPhi, &b_pfType1mets_unCPhi);
   fChain->SetBranchAddress("pfType1mets_unCPt", &pfType1mets_unCPt, &b_pfType1mets_unCPt);
   fChain->SetBranchAddress("pfType1mets_gen_et", &pfType1mets_gen_et, &b_pfType1mets_gen_et);
   fChain->SetBranchAddress("pfType1mets_gen_phi", &pfType1mets_gen_phi, &b_pfType1mets_gen_phi);
   fChain->SetBranchAddress("Npfcand", &Npfcand, &b_Npfcand);
   fChain->SetBranchAddress("pfcand_pdgId", &pfcand_pdgId, &b_pfcand_pdgId);
   fChain->SetBranchAddress("pfcand_pt", &pfcand_pt, &b_pfcand_pt);
   fChain->SetBranchAddress("pfcand_eta", &pfcand_eta, &b_pfcand_eta);
   fChain->SetBranchAddress("pfcand_phi", &pfcand_phi, &b_pfcand_phi);
   fChain->SetBranchAddress("pfcand_energy", &pfcand_energy, &b_pfcand_energy);
   fChain->SetBranchAddress("pfcand_charge", &pfcand_charge, &b_pfcand_charge);
   fChain->SetBranchAddress("pfcand_dz", &pfcand_dz, &b_pfcand_dz);
   fChain->SetBranchAddress("pfcand_dxy", &pfcand_dxy, &b_pfcand_dxy);
   fChain->SetBranchAddress("pfcand_fromPV", &pfcand_fromPV, &b_pfcand_fromPV);
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
   fChain->SetBranchAddress("photons_passElectronVeto", &photons_passElectronVeto, &b_photons_passElectronVeto);
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
   fChain->SetBranchAddress("Ntaus", &Ntaus, &b_Ntaus);
   fChain->SetBranchAddress("taus_energy", &taus_energy, &b_taus_energy);
   fChain->SetBranchAddress("taus_et", &taus_et, &b_taus_et);
   fChain->SetBranchAddress("taus_eta", &taus_eta, &b_taus_eta);
   fChain->SetBranchAddress("taus_phi", &taus_phi, &b_taus_phi);
   fChain->SetBranchAddress("taus_pt", &taus_pt, &b_taus_pt);
   fChain->SetBranchAddress("taus_px", &taus_px, &b_taus_px);
   fChain->SetBranchAddress("taus_py", &taus_py, &b_taus_py);
   fChain->SetBranchAddress("taus_pz", &taus_pz, &b_taus_pz);
   fChain->SetBranchAddress("taus_status", &taus_status, &b_taus_status);
   fChain->SetBranchAddress("taus_theta", &taus_theta, &b_taus_theta);
   fChain->SetBranchAddress("taus_charge", &taus_charge, &b_taus_charge);
   fChain->SetBranchAddress("taus_leadChargedHadrCand_pt", &taus_leadChargedHadrCand_pt, &b_taus_leadChargedHadrCand_pt);
   fChain->SetBranchAddress("taus_leadChargedHadrCand_charge", &taus_leadChargedHadrCand_charge, &b_taus_leadChargedHadrCand_charge);
   fChain->SetBranchAddress("taus_leadChargedHadrCand_eta", &taus_leadChargedHadrCand_eta, &b_taus_leadChargedHadrCand_eta);
   fChain->SetBranchAddress("taus_leadChargedHadrCand_phi", &taus_leadChargedHadrCand_phi, &b_taus_leadChargedHadrCand_phi);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("experimentType", &experimentType, &b_experimentType);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("weightIndex", &weightIndex, &b_weightIndex);
   fChain->SetBranchAddress("weightVector", &weightVector, &b_weightVector);
   fChain->SetBranchAddress("model_params", &model_params, &b_model_params);

}

