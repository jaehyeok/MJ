#include <vector>
#include <string>
#include "TTree.h"
#include "TChain.h"

   // Declaration of leaf types
   vector<float>   *trigger_prescalevalue;
   vector<string>  *trigger_name;
   vector<float>   *trigger_decision;
   vector<string>  *trigger_lastfiltername;
   vector<vector<float> > *triggerobject_pt;
   vector<vector<float> > *triggerobject_px;
   vector<vector<float> > *triggerobject_py;
   vector<vector<float> > *triggerobject_pz;
   vector<vector<float> > *triggerobject_et;
   vector<vector<float> > *triggerobject_energy;
   vector<vector<float> > *triggerobject_phi;
   vector<vector<float> > *triggerobject_eta;
   vector<vector<string> > *triggerobject_collectionname;
   vector<float>   *standalone_triggerobject_pt;
   vector<float>   *standalone_triggerobject_px;
   vector<float>   *standalone_triggerobject_py;
   vector<float>   *standalone_triggerobject_pz;
   vector<float>   *standalone_triggerobject_et;
   vector<float>   *standalone_triggerobject_energy;
   vector<float>   *standalone_triggerobject_phi;
   vector<float>   *standalone_triggerobject_eta;
   vector<string>  *standalone_triggerobject_collectionname;
   vector<float>   *L1trigger_bit;
   vector<float>   *L1trigger_techTrigger;
   vector<float>   *L1trigger_prescalevalue;
   vector<string>  *L1trigger_name;
   vector<string>  *L1trigger_alias;
   vector<float>   *L1trigger_decision;
   vector<float>   *L1trigger_decision_nomask;
   vector<float>   *els_conversion_dist;
   vector<float>   *els_conversion_dcot;
   vector<float>   *els_PFchargedHadronIsoR03;
   vector<float>   *els_PFphotonIsoR03;
   vector<float>   *els_PFneutralHadronIsoR03;
   vector<bool>    *els_hasMatchedConversion;
   vector<float>   *pf_els_PFchargedHadronIsoR03;
   vector<float>   *pf_els_PFphotonIsoR03;
   vector<float>   *pf_els_PFneutralHadronIsoR03;
   vector<bool>    *pf_els_hasMatchedConversion;
   Int_t           trk_nTOBTEC;
   Float_t         trk_ratioAllTOBTEC;
   Float_t         trk_ratioJetTOBTEC;
   Int_t           hbhefilter_decision;
   Int_t           trackingfailurefilter_decision;
   Int_t           cschalofilter_decision;
   Int_t           ecalTPfilter_decision;
   Int_t           ecalBEfilter_decision;
   Int_t           scrapingVeto_decision;
   Int_t           greedymuonfilter_decision;
   Int_t           inconsistentPFmuonfilter_decision;
   Int_t           hcallaserfilter_decision;
   Int_t           ecallaserfilter_decision;
   Int_t           eenoisefilter_decision;
   Int_t           eebadscfilter_decision;
   Int_t           trackercoherentnoisefilter1_decision;
   Int_t           trackercoherentnoisefilter2_decision;
   Int_t           trackertoomanyclustersfilter_decision;
   Int_t           trackertoomanytripletsfilter_decision;
   Int_t           trackertoomanyseedsfilter_decision;
   Int_t           passprescalePFHT350filter_decision;
   Int_t           passprescaleHT250filter_decision;
   Int_t           passprescaleHT300filter_decision;
   Int_t           passprescaleHT350filter_decision;
   Int_t           passprescaleHT400filter_decision;
   Int_t           passprescaleHT450filter_decision;
   Int_t           passprescaleJet30MET80filter_decision;
   Float_t         MPT;
   Float_t         genHT;
   vector<float>   *jets_AK5PFclean_corrL2L3;
   vector<float>   *jets_AK5PFclean_corrL2L3Residual;
   vector<float>   *jets_AK5PFclean_corrL1FastL2L3;
   vector<float>   *jets_AK5PFclean_corrL1L2L3;
   vector<float>   *jets_AK5PFclean_corrL1FastL2L3Residual;
   vector<float>   *jets_AK5PFclean_corrL1L2L3Residual;
   vector<float>   *jets_AK5PFclean_Uncert;
   vector<vector<float> > *PU_zpositions;
   vector<vector<float> > *PU_sumpT_lowpT;
   vector<vector<float> > *PU_sumpT_highpT;
   vector<vector<int> > *PU_ntrks_lowpT;
   vector<vector<int> > *PU_ntrks_highpT;
   vector<int>     *PU_NumInteractions;
   vector<int>     *PU_bunchCrossing;
   vector<float>   *PU_TrueNumInteractions;
   Float_t         rho_kt6PFJetsForIsolation2011;
   Float_t         rho_kt6PFJetsForIsolation2012;
   Float_t         pfmets_fullSignif;
   Float_t         pfmets_fullSignifCov00;
   Float_t         pfmets_fullSignifCov10;
   Float_t         pfmets_fullSignifCov11;
   Float_t         softjetUp_dMEx;
   Float_t         softjetUp_dMEy;
   vector<float>   *pdfweights_cteq;
   vector<float>   *pdfweights_mstw;
   vector<float>   *pdfweights_nnpdf;
   vector<float>   *photon_chIsoValues;
   vector<float>   *photon_phIsoValues;
   vector<float>   *photon_nhIsoValues;
   vector<bool>    *photon_passElectronVeto;
   vector<vector<float> > *puJet_rejectionBeta;
   vector<vector<float> > *puJet_rejectionMVA;
   Float_t         pfmets_fullSignif_2012;
   Float_t         pfmets_fullSignifCov00_2012;
   Float_t         pfmets_fullSignifCov10_2012;
   Float_t         pfmets_fullSignifCov11_2012;
   Float_t         pfmets_fullSignif_2012_dataRes;
   Float_t         pfmets_fullSignifCov00_2012_dataRes;
   Float_t         pfmets_fullSignifCov10_2012_dataRes;
   Float_t         pfmets_fullSignifCov11_2012_dataRes;
   vector<float>   *isotk_pt;
   vector<float>   *isotk_phi;
   vector<float>   *isotk_eta;
   vector<float>   *isotk_iso;
   vector<float>   *isotk_dzpv;
   vector<int>     *isotk_charge;

   // List of branches
   TBranch        *b_trigger_prescalevalue;   //!
   TBranch        *b_trigger_name;   //!
   TBranch        *b_trigger_decision;   //!
   TBranch        *b_trigger_lastfiltername;   //!
   TBranch        *b_triggerobject_pt;   //!
   TBranch        *b_triggerobject_px;   //!
   TBranch        *b_triggerobject_py;   //!
   TBranch        *b_triggerobject_pz;   //!
   TBranch        *b_triggerobject_et;   //!
   TBranch        *b_triggerobject_energy;   //!
   TBranch        *b_triggerobject_phi;   //!
   TBranch        *b_triggerobject_eta;   //!
   TBranch        *b_triggerobject_collectionname;   //!
   TBranch        *b_standalone_triggerobject_pt;   //!
   TBranch        *b_standalone_triggerobject_px;   //!
   TBranch        *b_standalone_triggerobject_py;   //!
   TBranch        *b_standalone_triggerobject_pz;   //!
   TBranch        *b_standalone_triggerobject_et;   //!
   TBranch        *b_standalone_triggerobject_energy;   //!
   TBranch        *b_standalone_triggerobject_phi;   //!
   TBranch        *b_standalone_triggerobject_eta;   //!
   TBranch        *b_standalone_triggerobject_collectionname;   //!
   TBranch        *b_L1trigger_bit;   //!
   TBranch        *b_L1trigger_techTrigger;   //!
   TBranch        *b_L1trigger_prescalevalue;   //!
   TBranch        *b_L1trigger_name;   //!
   TBranch        *b_L1trigger_alias;   //!
   TBranch        *b_L1trigger_decision;   //!
   TBranch        *b_L1trigger_decision_nomask;   //!
   TBranch        *b_els_conversion_dist;   //!
   TBranch        *b_els_conversion_dcot;   //!
   TBranch        *b_els_PFchargedHadronIsoR03;   //!
   TBranch        *b_els_PFphotonIsoR03;   //!
   TBranch        *b_els_PFneutralHadronIsoR03;   //!
   TBranch        *b_els_hasMatchedConversion;   //!
   TBranch        *b_pf_els_PFchargedHadronIsoR03;   //!
   TBranch        *b_pf_els_PFphotonIsoR03;   //!
   TBranch        *b_pf_els_PFneutralHadronIsoR03;   //!
   TBranch        *b_pf_els_hasMatchedConversion;   //!
   TBranch        *b_Ctrk_nTOBTEC;   //!
   TBranch        *b_trk_ratioAllTOBTEC;   //!
   TBranch        *b_trk_ratioJetTOBTEC;   //!
   TBranch        *b_hbhefilter_decision;   //!
   TBranch        *b_trackingfailurefilter_decision;   //!
   TBranch        *b_cschalofilter_decision;   //!
   TBranch        *b_ecalTPfilter_decision;   //!
   TBranch        *b_ecalBEfilter_decision;   //!
   TBranch        *b_scrapingVeto_decision;   //!
   TBranch        *b_greedymuonfilter_decision;   //!
   TBranch        *b_inconsistentPFmuonfilter_decision;   //!
   TBranch        *b_hcallaserfilter_decision;   //!
   TBranch        *b_ecallaserfilter_decision;   //!
   TBranch        *b_eenoisefilter_decision;   //!
   TBranch        *b_eebadscfilter_decision;   //!
   TBranch        *b_trackercoherentnoisefilter1 ;   //!
   TBranch        *b_trackercoherentnoisefilter2 ;   //!
   TBranch        *b_trackertoomanyclustersfilter;   //!
   TBranch        *b_trackertoomanytripletsfilter;   //!
   TBranch        *b_trackertoomanyseedsfilter ;   //!
   TBranch        *b_passprescalePFHT350filter_decision;   //!
   TBranch        *b_passprescaleHT250filter_decision;   //!
   TBranch        *b_passprescaleHT300filter_decision;   //!
   TBranch        *b_passprescaleHT350filter_decision;   //!
   TBranch        *b_passprescaleHT400filter_decision;   //!
   TBranch        *b_passprescaleHT450filter_decision;   //!
   TBranch        *b_passprescaleJet30MET80filter_decision;   //!
   TBranch        *b_MPT;   //!
   TBranch        *b_genHT;   //!
   TBranch        *b_jets_AK5PFclean_corrL2L3;   //!
   TBranch        *b_jets_AK5PFclean_corrL2L3Residual;   //!
   TBranch        *b_jets_AK5PFclean_corrL1FastL2L3;   //!
   TBranch        *b_jets_AK5PFclean_corrL1L2L3;   //!
   TBranch        *b_jets_AK5PFclean_corrL1FastL2L3Residual;   //!
   TBranch        *b_jets_AK5PFclean_corrL1L2L3Residual;   //!
   TBranch        *b_jets_AK5PFclean_Uncert;   //!
   TBranch        *b_PU_zpositions;   //!
   TBranch        *b_PU_sumpT_lowpT;   //!
   TBranch        *b_PU_sumpT_highpT;   //!
   TBranch        *b_PU_ntrks_lowpT;   //!
   TBranch        *b_PU_ntrks_highpT;   //!
   TBranch        *b_PU_NumInteractions;   //!
   TBranch        *b_PU_bunchCrossing;   //!
   TBranch        *b_PU_TrueNumInteractions;   //!
   TBranch        *b_rho_kt6PFJetsForIsolation2011;   //!
   TBranch        *b_rho_kt6PFJetsForIsolation2012;   //!
   TBranch        *b_pfmets_fullSignif;   //!
   TBranch        *b_pfmets_fullSignifCov00;   //!
   TBranch        *b_pfmets_fullSignifCov10;   //!
   TBranch        *b_pfmets_fullSignifCov11;   //!
   TBranch        *b_softjetUp_dMEx;   //!
   TBranch        *b_softjetUp_dMEy;   //!
   TBranch        *b_pdfweights_cteq;   //!
   TBranch        *b_pdfweights_mstw;   //!
   TBranch        *b_pdfweights_nnpdf;   //!
   TBranch        *b_photon_chIsoValues;   //!
   TBranch        *b_photon_phIsoValues;   //!
   TBranch        *b_photon_nhIsoValues;   //!
   TBranch        *b_photon_passElectronVeto;   //!
   TBranch        *b_puJet_rejectionBeta;   //!
   TBranch        *b_puJet_rejectionMVA;   //!
   TBranch        *b_pfmets_fullSignif_2012;   //!
   TBranch        *b_pfmets_fullSignifCov00_2012;   //!
   TBranch        *b_pfmets_fullSignifCov10_2012;   //!
   TBranch        *b_pfmets_fullSignifCov11_2012;   //!
   TBranch        *b_pfmets_fullSignif_2012_dataRes;   //!
   TBranch        *b_pfmets_fullSignifCov00_2012_dataRes;   //!
   TBranch        *b_pfmets_fullSignifCov10_2012_dataRes;   //!
   TBranch        *b_pfmets_fullSignifCov11_2012_dataRes;   //!
   TBranch        *b_isotk_pt;   //!
   TBranch        *b_isotk_phi;   //!
   TBranch        *b_isotk_eta;   //!
   TBranch        *b_isotk_iso;   //!
   TBranch        *b_isotk_dzpv;   //!
   TBranch        *b_isotk_charge;   //!
   
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
   UInt_t          Njets_AK5PF;
   vector<float>   *jets_AK5PF_status;
   vector<float>   *jets_AK5PF_phi;
   vector<float>   *jets_AK5PF_pt;
   vector<float>   *jets_AK5PF_pz;
   vector<float>   *jets_AK5PF_px;
   vector<float>   *jets_AK5PF_py;
   vector<float>   *jets_AK5PF_eta;
   vector<float>   *jets_AK5PF_theta;
   vector<float>   *jets_AK5PF_et;
   vector<float>   *jets_AK5PF_energy;
   vector<float>   *jets_AK5PF_parton_Id;
   vector<float>   *jets_AK5PF_parton_motherId;
   vector<float>   *jets_AK5PF_parton_pt;
   vector<float>   *jets_AK5PF_parton_phi;
   vector<float>   *jets_AK5PF_parton_eta;
   vector<float>   *jets_AK5PF_parton_Energy;
   vector<float>   *jets_AK5PF_parton_mass;
   vector<float>   *jets_AK5PF_gen_et;
   vector<float>   *jets_AK5PF_gen_pt;
   vector<float>   *jets_AK5PF_gen_eta;
   vector<float>   *jets_AK5PF_gen_phi;
   vector<float>   *jets_AK5PF_gen_mass;
   vector<float>   *jets_AK5PF_gen_Energy;
   vector<float>   *jets_AK5PF_gen_Id;
   vector<float>   *jets_AK5PF_gen_motherID;
   vector<float>   *jets_AK5PF_gen_threeCharge;
   vector<float>   *jets_AK5PF_partonFlavour;
   vector<float>   *jets_AK5PF_btag_TC_highPur;
   vector<float>   *jets_AK5PF_btag_TC_highEff;
   vector<float>   *jets_AK5PF_btag_jetProb;
   vector<float>   *jets_AK5PF_btag_jetBProb;
   vector<float>   *jets_AK5PF_btag_softEle;
   vector<float>   *jets_AK5PF_btag_softMuon;
   vector<float>   *jets_AK5PF_btag_secVertexHighPur;
   vector<float>   *jets_AK5PF_btag_secVertexHighEff;
   vector<float>   *jets_AK5PF_btag_secVertexCombined;
   vector<float>   *jets_AK5PF_jetCharge;
   vector<float>   *jets_AK5PF_chgEmE;
   vector<float>   *jets_AK5PF_chgHadE;
   vector<float>   *jets_AK5PF_photonEnergy;
   vector<float>   *jets_AK5PF_chgMuE;
   vector<float>   *jets_AK5PF_chg_Mult;
   vector<float>   *jets_AK5PF_neutralEmE;
   vector<float>   *jets_AK5PF_neutralHadE;
   vector<float>   *jets_AK5PF_neutral_Mult;
   vector<float>   *jets_AK5PF_mu_Mult;
   vector<float>   *jets_AK5PF_emf;
   vector<float>   *jets_AK5PF_ehf;
   vector<float>   *jets_AK5PF_n60;
   vector<float>   *jets_AK5PF_n90;
   vector<float>   *jets_AK5PF_etaetaMoment;
   vector<float>   *jets_AK5PF_etaphiMoment;
   vector<float>   *jets_AK5PF_phiphiMoment;
   vector<float>   *jets_AK5PF_n90Hits;
   vector<float>   *jets_AK5PF_fHPD;
   vector<float>   *jets_AK5PF_fRBX;
   vector<float>   *jets_AK5PF_hitsInN90;
   vector<float>   *jets_AK5PF_nECALTowers;
   vector<float>   *jets_AK5PF_nHCALTowers;
   vector<float>   *jets_AK5PF_fSubDetector1;
   vector<float>   *jets_AK5PF_fSubDetector2;
   vector<float>   *jets_AK5PF_fSubDetector3;
   vector<float>   *jets_AK5PF_fSubDetector4;
   vector<float>   *jets_AK5PF_area;
   vector<float>   *jets_AK5PF_corrFactorRaw;
   vector<float>   *jets_AK5PF_rawPt;
   vector<float>   *jets_AK5PF_mass;
   UInt_t          Njets_AK5PFclean;
   vector<float>   *jets_AK5PFclean_status;
   vector<float>   *jets_AK5PFclean_phi;
   vector<float>   *jets_AK5PFclean_pt;
   vector<float>   *jets_AK5PFclean_pz;
   vector<float>   *jets_AK5PFclean_px;
   vector<float>   *jets_AK5PFclean_py;
   vector<float>   *jets_AK5PFclean_eta;
   vector<float>   *jets_AK5PFclean_theta;
   vector<float>   *jets_AK5PFclean_et;
   vector<float>   *jets_AK5PFclean_energy;
   vector<float>   *jets_AK5PFclean_parton_Id;
   vector<float>   *jets_AK5PFclean_parton_motherId;
   vector<float>   *jets_AK5PFclean_parton_pt;
   vector<float>   *jets_AK5PFclean_parton_phi;
   vector<float>   *jets_AK5PFclean_parton_eta;
   vector<float>   *jets_AK5PFclean_parton_Energy;
   vector<float>   *jets_AK5PFclean_parton_mass;
   vector<float>   *jets_AK5PFclean_gen_et;
   vector<float>   *jets_AK5PFclean_gen_pt;
   vector<float>   *jets_AK5PFclean_gen_eta;
   vector<float>   *jets_AK5PFclean_gen_phi;
   vector<float>   *jets_AK5PFclean_gen_mass;
   vector<float>   *jets_AK5PFclean_gen_Energy;
   vector<float>   *jets_AK5PFclean_gen_Id;
   vector<float>   *jets_AK5PFclean_partonFlavour;
   vector<float>   *jets_AK5PFclean_btag_TC_highPur;
   vector<float>   *jets_AK5PFclean_btag_TC_highEff;
   vector<float>   *jets_AK5PFclean_btag_jetProb;
   vector<float>   *jets_AK5PFclean_btag_jetBProb;
   vector<float>   *jets_AK5PFclean_btag_softEle;
   vector<float>   *jets_AK5PFclean_btag_softMuon;
   vector<float>   *jets_AK5PFclean_btag_secVertexHighPur;
   vector<float>   *jets_AK5PFclean_btag_secVertexHighEff;
   vector<float>   *jets_AK5PFclean_btag_secVertexCombined;
   vector<float>   *jets_AK5PFclean_jetCharge;
   vector<float>   *jets_AK5PFclean_chgEmE;
   vector<float>   *jets_AK5PFclean_chgHadE;
   vector<float>   *jets_AK5PFclean_photonEnergy;
   vector<float>   *jets_AK5PFclean_chgMuE;
   vector<float>   *jets_AK5PFclean_chg_Mult;
   vector<float>   *jets_AK5PFclean_neutralEmE;
   vector<float>   *jets_AK5PFclean_neutralHadE;
   vector<float>   *jets_AK5PFclean_neutral_Mult;
   vector<float>   *jets_AK5PFclean_mu_Mult;
   vector<float>   *jets_AK5PFclean_emf;
   vector<float>   *jets_AK5PFclean_ehf;
   vector<float>   *jets_AK5PFclean_n60;
   vector<float>   *jets_AK5PFclean_n90;
   vector<float>   *jets_AK5PFclean_etaetaMoment;
   vector<float>   *jets_AK5PFclean_etaphiMoment;
   vector<float>   *jets_AK5PFclean_phiphiMoment;
   vector<float>   *jets_AK5PFclean_n90Hits;
   vector<float>   *jets_AK5PFclean_fHPD;
   vector<float>   *jets_AK5PFclean_fRBX;
   vector<float>   *jets_AK5PFclean_hitsInN90;
   vector<float>   *jets_AK5PFclean_nECALTowers;
   vector<float>   *jets_AK5PFclean_nHCALTowers;
   vector<float>   *jets_AK5PFclean_fSubDetector1;
   vector<float>   *jets_AK5PFclean_fSubDetector2;
   vector<float>   *jets_AK5PFclean_fSubDetector3;
   vector<float>   *jets_AK5PFclean_fSubDetector4;
   vector<float>   *jets_AK5PFclean_area;
   vector<float>   *jets_AK5PFclean_corrFactorRaw;
   vector<float>   *jets_AK5PFclean_rawPt;
   vector<float>   *jets_AK5PFclean_mass;
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
   UInt_t          Nmc_pdf;
   vector<float>   *mc_pdf_x1;
   vector<float>   *mc_pdf_x2;
   vector<float>   *mc_pdf_q;
   vector<float>   *mc_pdf_id1;
   vector<float>   *mc_pdf_id2;
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
   UInt_t          NmetsHO;
   vector<float>   *metsHO_et;
   vector<float>   *metsHO_phi;
   vector<float>   *metsHO_ex;
   vector<float>   *metsHO_ey;
   vector<float>   *metsHO_sumEt;
   UInt_t          Nmets_AK5;
   vector<float>   *mets_AK5_et;
   vector<float>   *mets_AK5_phi;
   vector<float>   *mets_AK5_ex;
   vector<float>   *mets_AK5_ey;
   vector<float>   *mets_AK5_gen_et;
   vector<float>   *mets_AK5_gen_phi;
   vector<float>   *mets_AK5_sign;
   vector<float>   *mets_AK5_sumEt;
   vector<float>   *mets_AK5_unCPhi;
   vector<float>   *mets_AK5_unCPt;
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
   UInt_t          NpfTypeINoXYCorrmets;
   vector<float>   *pfTypeINoXYCorrmets_et;
   vector<float>   *pfTypeINoXYCorrmets_phi;
   vector<float>   *pfTypeINoXYCorrmets_ex;
   vector<float>   *pfTypeINoXYCorrmets_ey;
   vector<float>   *pfTypeINoXYCorrmets_gen_et;
   vector<float>   *pfTypeINoXYCorrmets_gen_phi;
   vector<float>   *pfTypeINoXYCorrmets_sign;
   vector<float>   *pfTypeINoXYCorrmets_sumEt;
   vector<float>   *pfTypeINoXYCorrmets_unCPhi;
   vector<float>   *pfTypeINoXYCorrmets_unCPt;
   UInt_t          NpfTypeIType0mets;
   vector<float>   *pfTypeIType0mets_et;
   vector<float>   *pfTypeIType0mets_phi;
   vector<float>   *pfTypeIType0mets_ex;
   vector<float>   *pfTypeIType0mets_ey;
   vector<float>   *pfTypeIType0mets_gen_et;
   vector<float>   *pfTypeIType0mets_gen_phi;
   vector<float>   *pfTypeIType0mets_sign;
   vector<float>   *pfTypeIType0mets_sumEt;
   vector<float>   *pfTypeIType0mets_unCPhi;
   vector<float>   *pfTypeIType0mets_unCPt;
   UInt_t          NpfTypeImets;
   vector<float>   *pfTypeImets_et;
   vector<float>   *pfTypeImets_phi;
   vector<float>   *pfTypeImets_ex;
   vector<float>   *pfTypeImets_ey;
   vector<float>   *pfTypeImets_gen_et;
   vector<float>   *pfTypeImets_gen_phi;
   vector<float>   *pfTypeImets_sign;
   vector<float>   *pfTypeImets_sumEt;
   vector<float>   *pfTypeImets_unCPhi;
   vector<float>   *pfTypeImets_unCPt;
   UInt_t          Npf_els;
   vector<float>   *pf_els_energy;
   vector<float>   *pf_els_et;
   vector<float>   *pf_els_eta;
   vector<float>   *pf_els_phi;
   vector<float>   *pf_els_pt;
   vector<float>   *pf_els_px;
   vector<float>   *pf_els_py;
   vector<float>   *pf_els_pz;
   vector<float>   *pf_els_status;
   vector<float>   *pf_els_theta;
   vector<float>   *pf_els_gen_id;
   vector<float>   *pf_els_gen_phi;
   vector<float>   *pf_els_gen_pt;
   vector<float>   *pf_els_gen_pz;
   vector<float>   *pf_els_gen_px;
   vector<float>   *pf_els_gen_py;
   vector<float>   *pf_els_gen_eta;
   vector<float>   *pf_els_gen_theta;
   vector<float>   *pf_els_gen_et;
   vector<float>   *pf_els_gen_mother_id;
   vector<float>   *pf_els_gen_mother_phi;
   vector<float>   *pf_els_gen_mother_pt;
   vector<float>   *pf_els_gen_mother_pz;
   vector<float>   *pf_els_gen_mother_px;
   vector<float>   *pf_els_gen_mother_py;
   vector<float>   *pf_els_gen_mother_eta;
   vector<float>   *pf_els_gen_mother_theta;
   vector<float>   *pf_els_gen_mother_et;
   vector<float>   *pf_els_tightId;
   vector<float>   *pf_els_looseId;
   vector<float>   *pf_els_robustTightId;
   vector<float>   *pf_els_robustLooseId;
   vector<float>   *pf_els_robustHighEnergyId;
   vector<float>   *pf_els_simpleEleId95relIso;
   vector<float>   *pf_els_simpleEleId90relIso;
   vector<float>   *pf_els_simpleEleId85relIso;
   vector<float>   *pf_els_simpleEleId80relIso;
   vector<float>   *pf_els_simpleEleId70relIso;
   vector<float>   *pf_els_simpleEleId60relIso;
   vector<float>   *pf_els_simpleEleId95cIso;
   vector<float>   *pf_els_simpleEleId90cIso;
   vector<float>   *pf_els_simpleEleId85cIso;
   vector<float>   *pf_els_simpleEleId80cIso;
   vector<float>   *pf_els_simpleEleId70cIso;
   vector<float>   *pf_els_simpleEleId60cIso;
   vector<float>   *pf_els_cIso;
   vector<float>   *pf_els_tIso;
   vector<float>   *pf_els_ecalIso;
   vector<float>   *pf_els_hcalIso;
   vector<float>   *pf_els_chargedHadronIso;
   vector<float>   *pf_els_photonIso;
   vector<float>   *pf_els_neutralHadronIso;
   vector<float>   *pf_els_chi2;
   vector<float>   *pf_els_charge;
   vector<float>   *pf_els_caloEnergy;
   vector<float>   *pf_els_hadOverEm;
   vector<float>   *pf_els_hcalOverEcalBc;
   vector<float>   *pf_els_eOverPIn;
   vector<float>   *pf_els_eSeedOverPOut;
   vector<float>   *pf_els_sigmaEtaEta;
   vector<float>   *pf_els_sigmaIEtaIEta;
   vector<float>   *pf_els_scEnergy;
   vector<float>   *pf_els_scRawEnergy;
   vector<float>   *pf_els_scSeedEnergy;
   vector<float>   *pf_els_scEta;
   vector<float>   *pf_els_scPhi;
   vector<float>   *pf_els_scEtaWidth;
   vector<float>   *pf_els_scPhiWidth;
   vector<float>   *pf_els_scE1x5;
   vector<float>   *pf_els_scE2x5Max;
   vector<float>   *pf_els_scE5x5;
   vector<float>   *pf_els_isEB;
   vector<float>   *pf_els_isEE;
   vector<float>   *pf_els_dEtaIn;
   vector<float>   *pf_els_dPhiIn;
   vector<float>   *pf_els_dEtaOut;
   vector<float>   *pf_els_dPhiOut;
   vector<float>   *pf_els_numvalhits;
   vector<float>   *pf_els_numlosthits;
   vector<float>   *pf_els_basicClustersSize;
   vector<float>   *pf_els_tk_pz;
   vector<float>   *pf_els_tk_pt;
   vector<float>   *pf_els_tk_phi;
   vector<float>   *pf_els_tk_eta;
   vector<float>   *pf_els_d0dum;
   vector<float>   *pf_els_dz;
   vector<float>   *pf_els_vx;
   vector<float>   *pf_els_vy;
   vector<float>   *pf_els_vz;
   vector<float>   *pf_els_ndof;
   vector<float>   *pf_els_ptError;
   vector<float>   *pf_els_d0dumError;
   vector<float>   *pf_els_dzError;
   vector<float>   *pf_els_etaError;
   vector<float>   *pf_els_phiError;
   vector<float>   *pf_els_tk_charge;
   vector<float>   *pf_els_core_ecalDrivenSeed;
   vector<float>   *pf_els_n_inner_layer;
   vector<float>   *pf_els_n_outer_layer;
   vector<float>   *pf_els_ctf_tk_id;
   vector<float>   *pf_els_ctf_tk_charge;
   vector<float>   *pf_els_ctf_tk_eta;
   vector<float>   *pf_els_ctf_tk_phi;
   vector<float>   *pf_els_fbrem;
   vector<float>   *pf_els_shFracInnerHits;
   vector<float>   *pf_els_dr03EcalRecHitSumEt;
   vector<float>   *pf_els_dr03HcalTowerSumEt;
   vector<float>   *pf_els_dr03HcalDepth1TowerSumEt;
   vector<float>   *pf_els_dr03HcalDepth2TowerSumEt;
   vector<float>   *pf_els_dr03TkSumPt;
   vector<float>   *pf_els_dr04EcalRecHitSumEt;
   vector<float>   *pf_els_dr04HcalTowerSumEt;
   vector<float>   *pf_els_dr04HcalDepth1TowerSumEt;
   vector<float>   *pf_els_dr04HcalDepth2TowerSumEt;
   vector<float>   *pf_els_dr04TkSumPt;
   vector<float>   *pf_els_cpx;
   vector<float>   *pf_els_cpy;
   vector<float>   *pf_els_cpz;
   vector<float>   *pf_els_vpx;
   vector<float>   *pf_els_vpy;
   vector<float>   *pf_els_vpz;
   vector<float>   *pf_els_cx;
   vector<float>   *pf_els_cy;
   vector<float>   *pf_els_cz;
   vector<float>   *pf_els_PATpassConversionVeto;
   UInt_t          Npf_mus;
   vector<float>   *pf_mus_energy;
   vector<float>   *pf_mus_et;
   vector<float>   *pf_mus_eta;
   vector<float>   *pf_mus_phi;
   vector<float>   *pf_mus_pt;
   vector<float>   *pf_mus_px;
   vector<float>   *pf_mus_py;
   vector<float>   *pf_mus_pz;
   vector<float>   *pf_mus_status;
   vector<float>   *pf_mus_theta;
   vector<float>   *pf_mus_gen_id;
   vector<float>   *pf_mus_gen_phi;
   vector<float>   *pf_mus_gen_pt;
   vector<float>   *pf_mus_gen_pz;
   vector<float>   *pf_mus_gen_px;
   vector<float>   *pf_mus_gen_py;
   vector<float>   *pf_mus_gen_eta;
   vector<float>   *pf_mus_gen_theta;
   vector<float>   *pf_mus_gen_et;
   vector<float>   *pf_mus_gen_mother_id;
   vector<float>   *pf_mus_gen_mother_phi;
   vector<float>   *pf_mus_gen_mother_pt;
   vector<float>   *pf_mus_gen_mother_pz;
   vector<float>   *pf_mus_gen_mother_px;
   vector<float>   *pf_mus_gen_mother_py;
   vector<float>   *pf_mus_gen_mother_eta;
   vector<float>   *pf_mus_gen_mother_theta;
   vector<float>   *pf_mus_gen_mother_et;
   vector<float>   *pf_mus_tkHits;
   vector<float>   *pf_mus_cIso;
   vector<float>   *pf_mus_tIso;
   vector<float>   *pf_mus_ecalIso;
   vector<float>   *pf_mus_hcalIso;
   vector<float>   *pf_mus_iso03_emVetoEt;
   vector<float>   *pf_mus_iso03_hadVetoEt;
   vector<float>   *pf_mus_calEnergyEm;
   vector<float>   *pf_mus_calEnergyHad;
   vector<float>   *pf_mus_calEnergyHo;
   vector<float>   *pf_mus_calEnergyEmS9;
   vector<float>   *pf_mus_calEnergyHadS9;
   vector<float>   *pf_mus_calEnergyHoS9;
   vector<float>   *pf_mus_iso03_sumPt;
   vector<float>   *pf_mus_iso03_emEt;
   vector<float>   *pf_mus_iso03_hadEt;
   vector<float>   *pf_mus_iso03_hoEt;
   vector<float>   *pf_mus_iso03_nTracks;
   vector<float>   *pf_mus_iso05_sumPt;
   vector<float>   *pf_mus_iso05_emEt;
   vector<float>   *pf_mus_iso05_hadEt;
   vector<float>   *pf_mus_iso05_hoEt;
   vector<float>   *pf_mus_iso05_nTracks;
   vector<float>   *pf_mus_neutralHadronIso;
   vector<float>   *pf_mus_chargedHadronIso;
   vector<float>   *pf_mus_photonIso;
   vector<float>   *pf_mus_charge;
   vector<float>   *pf_mus_cm_chi2;
   vector<float>   *pf_mus_cm_ndof;
   vector<float>   *pf_mus_cm_chg;
   vector<float>   *pf_mus_cm_pt;
   vector<float>   *pf_mus_cm_px;
   vector<float>   *pf_mus_cm_py;
   vector<float>   *pf_mus_cm_pz;
   vector<float>   *pf_mus_cm_eta;
   vector<float>   *pf_mus_cm_phi;
   vector<float>   *pf_mus_cm_theta;
   vector<float>   *pf_mus_cm_d0dum;
   vector<float>   *pf_mus_cm_dz;
   vector<float>   *pf_mus_cm_vx;
   vector<float>   *pf_mus_cm_vy;
   vector<float>   *pf_mus_cm_vz;
   vector<float>   *pf_mus_cm_numvalhits;
   vector<float>   *pf_mus_cm_numlosthits;
   vector<float>   *pf_mus_cm_numvalMuonhits;
   vector<float>   *pf_mus_cm_d0dumErr;
   vector<float>   *pf_mus_cm_dzErr;
   vector<float>   *pf_mus_cm_ptErr;
   vector<float>   *pf_mus_cm_etaErr;
   vector<float>   *pf_mus_cm_phiErr;
   vector<float>   *pf_mus_tk_id;
   vector<float>   *pf_mus_tk_chi2;
   vector<float>   *pf_mus_tk_ndof;
   vector<float>   *pf_mus_tk_chg;
   vector<float>   *pf_mus_tk_pt;
   vector<float>   *pf_mus_tk_px;
   vector<float>   *pf_mus_tk_py;
   vector<float>   *pf_mus_tk_pz;
   vector<float>   *pf_mus_tk_eta;
   vector<float>   *pf_mus_tk_phi;
   vector<float>   *pf_mus_tk_theta;
   vector<float>   *pf_mus_tk_d0dum;
   vector<float>   *pf_mus_tk_dz;
   vector<float>   *pf_mus_tk_vx;
   vector<float>   *pf_mus_tk_vy;
   vector<float>   *pf_mus_tk_vz;
   vector<float>   *pf_mus_tk_numvalhits;
   vector<float>   *pf_mus_tk_numlosthits;
   vector<float>   *pf_mus_tk_d0dumErr;
   vector<float>   *pf_mus_tk_dzErr;
   vector<float>   *pf_mus_tk_ptErr;
   vector<float>   *pf_mus_tk_etaErr;
   vector<float>   *pf_mus_tk_phiErr;
   vector<float>   *pf_mus_tk_numvalPixelhits;
   vector<float>   *pf_mus_tk_numpixelWthMeasr;
   vector<float>   *pf_mus_stamu_chi2;
   vector<float>   *pf_mus_stamu_ndof;
   vector<float>   *pf_mus_stamu_chg;
   vector<float>   *pf_mus_stamu_pt;
   vector<float>   *pf_mus_stamu_px;
   vector<float>   *pf_mus_stamu_py;
   vector<float>   *pf_mus_stamu_pz;
   vector<float>   *pf_mus_stamu_eta;
   vector<float>   *pf_mus_stamu_phi;
   vector<float>   *pf_mus_stamu_theta;
   vector<float>   *pf_mus_stamu_d0dum;
   vector<float>   *pf_mus_stamu_dz;
   vector<float>   *pf_mus_stamu_vx;
   vector<float>   *pf_mus_stamu_vy;
   vector<float>   *pf_mus_stamu_vz;
   vector<float>   *pf_mus_stamu_numvalhits;
   vector<float>   *pf_mus_stamu_numlosthits;
   vector<float>   *pf_mus_stamu_d0dumErr;
   vector<float>   *pf_mus_stamu_dzErr;
   vector<float>   *pf_mus_stamu_ptErr;
   vector<float>   *pf_mus_stamu_etaErr;
   vector<float>   *pf_mus_stamu_phiErr;
   vector<float>   *pf_mus_num_matches;
   vector<float>   *pf_mus_isTrackerMuon;
   vector<float>   *pf_mus_isStandAloneMuon;
   vector<float>   *pf_mus_isCaloMuon;
   vector<float>   *pf_mus_isGlobalMuon;
   vector<float>   *pf_mus_isElectron;
   vector<float>   *pf_mus_isConvertedPhoton;
   vector<float>   *pf_mus_isPhoton;
   vector<float>   *pf_mus_id_All;
   vector<float>   *pf_mus_id_AllGlobalMuons;
   vector<float>   *pf_mus_id_AllStandAloneMuons;
   vector<float>   *pf_mus_id_AllTrackerMuons;
   vector<float>   *pf_mus_id_TrackerMuonArbitrated;
   vector<float>   *pf_mus_id_AllArbitrated;
   vector<float>   *pf_mus_id_GlobalMuonPromptTight;
   vector<float>   *pf_mus_id_TMLastStationLoose;
   vector<float>   *pf_mus_id_TMLastStationTight;
   vector<float>   *pf_mus_id_TM2DCompatibilityLoose;
   vector<float>   *pf_mus_id_TM2DCompatibilityTight;
   vector<float>   *pf_mus_id_TMOneStationLoose;
   vector<float>   *pf_mus_id_TMOneStationTight;
   vector<float>   *pf_mus_id_TMLastStationOptimizedLowPtLoose;
   vector<float>   *pf_mus_id_TMLastStationOptimizedLowPtTight;
   vector<float>   *pf_mus_tk_LayersWithMeasurement;
   vector<float>   *pf_mus_tk_PixelLayersWithMeasurement;
   vector<float>   *pf_mus_tk_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *pf_mus_tk_LayersWithoutMeasurement;
   vector<float>   *pf_mus_tk_ExpectedHitsInner;
   vector<float>   *pf_mus_tk_ExpectedHitsOuter;
   vector<float>   *pf_mus_cm_LayersWithMeasurement;
   vector<float>   *pf_mus_cm_PixelLayersWithMeasurement;
   vector<float>   *pf_mus_cm_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *pf_mus_cm_LayersWithoutMeasurement;
   vector<float>   *pf_mus_cm_ExpectedHitsInner;
   vector<float>   *pf_mus_cm_ExpectedHitsOuter;
   vector<float>   *pf_mus_picky_LayersWithMeasurement;
   vector<float>   *pf_mus_picky_PixelLayersWithMeasurement;
   vector<float>   *pf_mus_picky_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *pf_mus_picky_LayersWithoutMeasurement;
   vector<float>   *pf_mus_picky_ExpectedHitsInner;
   vector<float>   *pf_mus_picky_ExpectedHitsOuter;
   vector<float>   *pf_mus_tpfms_LayersWithMeasurement;
   vector<float>   *pf_mus_tpfms_PixelLayersWithMeasurement;
   vector<float>   *pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit;
   vector<float>   *pf_mus_tpfms_LayersWithoutMeasurement;
   vector<float>   *pf_mus_tpfms_ExpectedHitsInner;
   vector<float>   *pf_mus_tpfms_ExpectedHitsOuter;
   vector<float>   *pf_mus_picky_id;
   vector<float>   *pf_mus_picky_chi2;
   vector<float>   *pf_mus_picky_ndof;
   vector<float>   *pf_mus_picky_chg;
   vector<float>   *pf_mus_picky_pt;
   vector<float>   *pf_mus_picky_px;
   vector<float>   *pf_mus_picky_py;
   vector<float>   *pf_mus_picky_pz;
   vector<float>   *pf_mus_picky_eta;
   vector<float>   *pf_mus_picky_phi;
   vector<float>   *pf_mus_picky_theta;
   vector<float>   *pf_mus_picky_d0dum;
   vector<float>   *pf_mus_picky_dz;
   vector<float>   *pf_mus_picky_vx;
   vector<float>   *pf_mus_picky_vy;
   vector<float>   *pf_mus_picky_vz;
   vector<float>   *pf_mus_picky_numvalhits;
   vector<float>   *pf_mus_picky_numlosthits;
   vector<float>   *pf_mus_picky_d0dumErr;
   vector<float>   *pf_mus_picky_dzErr;
   vector<float>   *pf_mus_picky_ptErr;
   vector<float>   *pf_mus_picky_etaErr;
   vector<float>   *pf_mus_picky_phiErr;
   vector<float>   *pf_mus_picky_numvalPixelhits;
   vector<float>   *pf_mus_tpfms_id;
   vector<float>   *pf_mus_tpfms_chi2;
   vector<float>   *pf_mus_tpfms_ndof;
   vector<float>   *pf_mus_tpfms_chg;
   vector<float>   *pf_mus_tpfms_pt;
   vector<float>   *pf_mus_tpfms_px;
   vector<float>   *pf_mus_tpfms_py;
   vector<float>   *pf_mus_tpfms_pz;
   vector<float>   *pf_mus_tpfms_eta;
   vector<float>   *pf_mus_tpfms_phi;
   vector<float>   *pf_mus_tpfms_theta;
   vector<float>   *pf_mus_tpfms_d0dum;
   vector<float>   *pf_mus_tpfms_dz;
   vector<float>   *pf_mus_tpfms_vx;
   vector<float>   *pf_mus_tpfms_vy;
   vector<float>   *pf_mus_tpfms_vz;
   vector<float>   *pf_mus_tpfms_numvalhits;
   vector<float>   *pf_mus_tpfms_numlosthits;
   vector<float>   *pf_mus_tpfms_d0dumErr;
   vector<float>   *pf_mus_tpfms_dzErr;
   vector<float>   *pf_mus_tpfms_ptErr;
   vector<float>   *pf_mus_tpfms_etaErr;
   vector<float>   *pf_mus_tpfms_phiErr;
   vector<float>   *pf_mus_tpfms_numvalPixelhits;
   vector<float>   *pf_mus_pfIsolationR03_sumChargedHadronPt;
   vector<float>   *pf_mus_pfIsolationR03_sumChargedParticlePt;
   vector<float>   *pf_mus_pfIsolationR03_sumNeutralHadronEt;
   vector<float>   *pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold;
   vector<float>   *pf_mus_pfIsolationR03_sumPhotonEt;
   vector<float>   *pf_mus_pfIsolationR03_sumPhotonEtHighThreshold;
   vector<float>   *pf_mus_pfIsolationR03_sumPUPt;
   vector<float>   *pf_mus_pfIsolationR04_sumChargedHadronPt;
   vector<float>   *pf_mus_pfIsolationR04_sumChargedParticlePt;
   vector<float>   *pf_mus_pfIsolationR04_sumNeutralHadronEt;
   vector<float>   *pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold;
   vector<float>   *pf_mus_pfIsolationR04_sumPhotonEt;
   vector<float>   *pf_mus_pfIsolationR04_sumPhotonEtHighThreshold;
   vector<float>   *pf_mus_pfIsolationR04_sumPUPt;
   vector<float>   *pf_mus_dB;
   vector<float>   *pf_mus_numberOfMatchedStations;
   vector<float>   *pf_mus_isPFMuon;
   UInt_t          Npf_photons;
   vector<float>   *pf_photons_energy;
   vector<float>   *pf_photons_et;
   vector<float>   *pf_photons_eta;
   vector<float>   *pf_photons_phi;
   vector<float>   *pf_photons_pt;
   vector<float>   *pf_photons_px;
   vector<float>   *pf_photons_py;
   vector<float>   *pf_photons_pz;
   vector<float>   *pf_photons_status;
   vector<float>   *pf_photons_theta;
   vector<float>   *pf_photons_hadOverEM;
   vector<float>   *pf_photons_hadTowOverEM;
   vector<float>   *pf_photons_scEnergy;
   vector<float>   *pf_photons_scRawEnergy;
   vector<float>   *pf_photons_scEta;
   vector<float>   *pf_photons_scPhi;
   vector<float>   *pf_photons_scEtaWidth;
   vector<float>   *pf_photons_scPhiWidth;
   vector<float>   *pf_photons_isAlsoElectron;
   vector<float>   *pf_photons_hasPixelSeed;
   vector<float>   *pf_photons_isConverted;
   vector<float>   *pf_photons_isEBGap;
   vector<float>   *pf_photons_isEEGap;
   vector<float>   *pf_photons_isEBEEGap;
   vector<float>   *pf_photons_isEBPho;
   vector<float>   *pf_photons_isEEPho;
   vector<float>   *pf_photons_maxEnergyXtal;
   vector<float>   *pf_photons_e1x5;
   vector<float>   *pf_photons_e2x5;
   vector<float>   *pf_photons_e3x3;
   vector<float>   *pf_photons_e5x5;
   vector<float>   *pf_photons_sigmaEtaEta;
   vector<float>   *pf_photons_sigmaIetaIeta;
   vector<float>   *pf_photons_r9;
   vector<float>   *pf_photons_chIso;
   vector<float>   *pf_photons_nhIso;
   vector<float>   *pf_photons_phIso;
   UInt_t          Npfcand;
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
   UInt_t          Npfmets;
   vector<float>   *pfmets_et;
   vector<float>   *pfmets_phi;
   vector<float>   *pfmets_ex;
   vector<float>   *pfmets_ey;
   vector<float>   *pfmets_gen_et;
   vector<float>   *pfmets_gen_phi;
   vector<float>   *pfmets_sign;
   vector<float>   *pfmets_sumEt;
   vector<float>   *pfmets_unCPhi;
   vector<float>   *pfmets_unCPt;
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
   UInt_t          Ntcmets;
   vector<float>   *tcmets_et;
   vector<float>   *tcmets_phi;
   vector<float>   *tcmets_ex;
   vector<float>   *tcmets_ey;
   vector<float>   *tcmets_sumEt;
   UInt_t          Ntracks;
   vector<float>   *tracks_chi2;
   vector<float>   *tracks_ndof;
   vector<float>   *tracks_chg;
   vector<float>   *tracks_pt;
   vector<float>   *tracks_px;
   vector<float>   *tracks_py;
   vector<float>   *tracks_pz;
   vector<float>   *tracks_eta;
   vector<float>   *tracks_phi;
   vector<float>   *tracks_d0dum;
   vector<float>   *tracks_dz;
   vector<float>   *tracks_vx;
   vector<float>   *tracks_vy;
   vector<float>   *tracks_vz;
   vector<float>   *tracks_numvalhits;
   vector<float>   *tracks_numlosthits;
   vector<float>   *tracks_d0dumErr;
   vector<float>   *tracks_dzErr;
   vector<float>   *tracks_ptErr;
   vector<float>   *tracks_etaErr;
   vector<float>   *tracks_phiErr;
   vector<float>   *tracks_highPurity;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          experimentType;
   UInt_t          bunchCrossing;
   UInt_t          orbitNumber;
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
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT10_px;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT10_py;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT10_pz;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT10_energy;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT10_phi;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT10_eta;
   vector<vector<int> > *fastjets_AK5PF_R1p2_R0p5pT10_index;
   vector<int>     *fastjets_AK5PF_R1p2_R0p5pT10_nconstituents;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT15_px;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT15_py;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT15_pz;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT15_energy;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT15_phi;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT15_eta;
   vector<vector<int> > *fastjets_AK5PF_R1p2_R0p5pT15_index;
   vector<int>     *fastjets_AK5PF_R1p2_R0p5pT15_nconstituents;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT20_px;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT20_py;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT20_pz;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT20_energy;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT20_phi;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT20_eta;
   vector<vector<int> > *fastjets_AK5PF_R1p2_R0p5pT20_index;
   vector<int>     *fastjets_AK5PF_R1p2_R0p5pT20_nconstituents;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT25_px;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT25_py;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT25_pz;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT25_energy;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT25_phi;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT25_eta;
   vector<vector<int> > *fastjets_AK5PF_R1p2_R0p5pT25_index;
   vector<int>     *fastjets_AK5PF_R1p2_R0p5pT25_nconstituents;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT30_px;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT30_py;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT30_pz;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT30_energy;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT30_phi;
   vector<float>   *fastjets_AK5PF_R1p2_R0p5pT30_eta;
   vector<vector<int> > *fastjets_AK5PF_R1p2_R0p5pT30_index;
   vector<int>     *fastjets_AK5PF_R1p2_R0p5pT30_nconstituents;

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
   TBranch        *b_Njets_AK5PF;   //!
   TBranch        *b_jets_AK5PF_status;   //!
   TBranch        *b_jets_AK5PF_phi;   //!
   TBranch        *b_jets_AK5PF_pt;   //!
   TBranch        *b_jets_AK5PF_pz;   //!
   TBranch        *b_jets_AK5PF_px;   //!
   TBranch        *b_jets_AK5PF_py;   //!
   TBranch        *b_jets_AK5PF_eta;   //!
   TBranch        *b_jets_AK5PF_theta;   //!
   TBranch        *b_jets_AK5PF_et;   //!
   TBranch        *b_jets_AK5PF_energy;   //!
   TBranch        *b_jets_AK5PF_parton_Id;   //!
   TBranch        *b_jets_AK5PF_parton_motherId;   //!
   TBranch        *b_jets_AK5PF_parton_pt;   //!
   TBranch        *b_jets_AK5PF_parton_phi;   //!
   TBranch        *b_jets_AK5PF_parton_eta;   //!
   TBranch        *b_jets_AK5PF_parton_Energy;   //!
   TBranch        *b_jets_AK5PF_parton_mass;   //!
   TBranch        *b_jets_AK5PF_gen_et;   //!
   TBranch        *b_jets_AK5PF_gen_pt;   //!
   TBranch        *b_jets_AK5PF_gen_eta;   //!
   TBranch        *b_jets_AK5PF_gen_phi;   //!
   TBranch        *b_jets_AK5PF_gen_mass;   //!
   TBranch        *b_jets_AK5PF_gen_Energy;   //!
   TBranch        *b_jets_AK5PF_gen_Id;   //!
   TBranch        *b_jets_AK5PF_gen_motherID;   //!
   TBranch        *b_jets_AK5PF_gen_threeCharge;   //!
   TBranch        *b_jets_AK5PF_partonFlavour;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PF_btag_jetProb;   //!
   TBranch        *b_jets_AK5PF_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PF_btag_softEle;   //!
   TBranch        *b_jets_AK5PF_btag_softMuon;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PF_jetCharge;   //!
   TBranch        *b_jets_AK5PF_chgEmE;   //!
   TBranch        *b_jets_AK5PF_chgHadE;   //!
   TBranch        *b_jets_AK5PF_photonEnergy;   //!
   TBranch        *b_jets_AK5PF_chgMuE;   //!
   TBranch        *b_jets_AK5PF_chg_Mult;   //!
   TBranch        *b_jets_AK5PF_neutralEmE;   //!
   TBranch        *b_jets_AK5PF_neutralHadE;   //!
   TBranch        *b_jets_AK5PF_neutral_Mult;   //!
   TBranch        *b_jets_AK5PF_mu_Mult;   //!
   TBranch        *b_jets_AK5PF_emf;   //!
   TBranch        *b_jets_AK5PF_ehf;   //!
   TBranch        *b_jets_AK5PF_n60;   //!
   TBranch        *b_jets_AK5PF_n90;   //!
   TBranch        *b_jets_AK5PF_etaetaMoment;   //!
   TBranch        *b_jets_AK5PF_etaphiMoment;   //!
   TBranch        *b_jets_AK5PF_phiphiMoment;   //!
   TBranch        *b_jets_AK5PF_n90Hits;   //!
   TBranch        *b_jets_AK5PF_fHPD;   //!
   TBranch        *b_jets_AK5PF_fRBX;   //!
   TBranch        *b_jets_AK5PF_hitsInN90;   //!
   TBranch        *b_jets_AK5PF_nECALTowers;   //!
   TBranch        *b_jets_AK5PF_nHCALTowers;   //!
   TBranch        *b_jets_AK5PF_fSubDetector1;   //!
   TBranch        *b_jets_AK5PF_fSubDetector2;   //!
   TBranch        *b_jets_AK5PF_fSubDetector3;   //!
   TBranch        *b_jets_AK5PF_fSubDetector4;   //!
   TBranch        *b_jets_AK5PF_area;   //!
   TBranch        *b_jets_AK5PF_corrFactorRaw;   //!
   TBranch        *b_jets_AK5PF_rawPt;   //!
   TBranch        *b_jets_AK5PF_mass;   //!
   TBranch        *b_Njets_AK5PFclean;   //!
   TBranch        *b_jets_AK5PFclean_status;   //!
   TBranch        *b_jets_AK5PFclean_phi;   //!
   TBranch        *b_jets_AK5PFclean_pt;   //!
   TBranch        *b_jets_AK5PFclean_pz;   //!
   TBranch        *b_jets_AK5PFclean_px;   //!
   TBranch        *b_jets_AK5PFclean_py;   //!
   TBranch        *b_jets_AK5PFclean_eta;   //!
   TBranch        *b_jets_AK5PFclean_theta;   //!
   TBranch        *b_jets_AK5PFclean_et;   //!
   TBranch        *b_jets_AK5PFclean_energy;   //!
   TBranch        *b_jets_AK5PFclean_parton_Id;   //!
   TBranch        *b_jets_AK5PFclean_parton_motherId;   //!
   TBranch        *b_jets_AK5PFclean_parton_pt;   //!
   TBranch        *b_jets_AK5PFclean_parton_phi;   //!
   TBranch        *b_jets_AK5PFclean_parton_eta;   //!
   TBranch        *b_jets_AK5PFclean_parton_Energy;   //!
   TBranch        *b_jets_AK5PFclean_parton_mass;   //!
   TBranch        *b_jets_AK5PFclean_gen_et;   //!
   TBranch        *b_jets_AK5PFclean_gen_pt;   //!
   TBranch        *b_jets_AK5PFclean_gen_eta;   //!
   TBranch        *b_jets_AK5PFclean_gen_phi;   //!
   TBranch        *b_jets_AK5PFclean_gen_mass;   //!
   TBranch        *b_jets_AK5PFclean_gen_Energy;   //!
   TBranch        *b_jets_AK5PFclean_gen_Id;   //!
   TBranch        *b_jets_AK5PFclean_partonFlavour;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PFclean_btag_softEle;   //!
   TBranch        *b_jets_AK5PFclean_btag_softMuon;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PFclean_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PFclean_jetCharge;   //!
   TBranch        *b_jets_AK5PFclean_chgEmE;   //!
   TBranch        *b_jets_AK5PFclean_chgHadE;   //!
   TBranch        *b_jets_AK5PFclean_photonEnergy;   //!
   TBranch        *b_jets_AK5PFclean_chgMuE;   //!
   TBranch        *b_jets_AK5PFclean_chg_Mult;   //!
   TBranch        *b_jets_AK5PFclean_neutralEmE;   //!
   TBranch        *b_jets_AK5PFclean_neutralHadE;   //!
   TBranch        *b_jets_AK5PFclean_neutral_Mult;   //!
   TBranch        *b_jets_AK5PFclean_mu_Mult;   //!
   TBranch        *b_jets_AK5PFclean_emf;   //!
   TBranch        *b_jets_AK5PFclean_ehf;   //!
   TBranch        *b_jets_AK5PFclean_n60;   //!
   TBranch        *b_jets_AK5PFclean_n90;   //!
   TBranch        *b_jets_AK5PFclean_etaetaMoment;   //!
   TBranch        *b_jets_AK5PFclean_etaphiMoment;   //!
   TBranch        *b_jets_AK5PFclean_phiphiMoment;   //!
   TBranch        *b_jets_AK5PFclean_n90Hits;   //!
   TBranch        *b_jets_AK5PFclean_fHPD;   //!
   TBranch        *b_jets_AK5PFclean_fRBX;   //!
   TBranch        *b_jets_AK5PFclean_hitsInN90;   //!
   TBranch        *b_jets_AK5PFclean_nECALTowers;   //!
   TBranch        *b_jets_AK5PFclean_nHCALTowers;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector1;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector2;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector3;   //!
   TBranch        *b_jets_AK5PFclean_fSubDetector4;   //!
   TBranch        *b_jets_AK5PFclean_area;   //!
   TBranch        *b_jets_AK5PFclean_corrFactorRaw;   //!
   TBranch        *b_jets_AK5PFclean_rawPt;   //!
   TBranch        *b_jets_AK5PFclean_mass;   //!
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
   TBranch        *b_Nmc_pdf;   //!
   TBranch        *b_mc_pdf_x1;   //!
   TBranch        *b_mc_pdf_x2;   //!
   TBranch        *b_mc_pdf_q;   //!
   TBranch        *b_mc_pdf_id1;   //!
   TBranch        *b_mc_pdf_id2;   //!
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
   TBranch        *b_NmetsHO;   //!
   TBranch        *b_metsHO_et;   //!
   TBranch        *b_metsHO_phi;   //!
   TBranch        *b_metsHO_ex;   //!
   TBranch        *b_metsHO_ey;   //!
   TBranch        *b_metsHO_sumEt;   //!
   TBranch        *b_Nmets_AK5;   //!
   TBranch        *b_mets_AK5_et;   //!
   TBranch        *b_mets_AK5_phi;   //!
   TBranch        *b_mets_AK5_ex;   //!
   TBranch        *b_mets_AK5_ey;   //!
   TBranch        *b_mets_AK5_gen_et;   //!
   TBranch        *b_mets_AK5_gen_phi;   //!
   TBranch        *b_mets_AK5_sign;   //!
   TBranch        *b_mets_AK5_sumEt;   //!
   TBranch        *b_mets_AK5_unCPhi;   //!
   TBranch        *b_mets_AK5_unCPt;   //!
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
   TBranch        *b_NpfTypeINoXYCorrmets;   //!
   TBranch        *b_pfTypeINoXYCorrmets_et;   //!
   TBranch        *b_pfTypeINoXYCorrmets_phi;   //!
   TBranch        *b_pfTypeINoXYCorrmets_ex;   //!
   TBranch        *b_pfTypeINoXYCorrmets_ey;   //!
   TBranch        *b_pfTypeINoXYCorrmets_gen_et;   //!
   TBranch        *b_pfTypeINoXYCorrmets_gen_phi;   //!
   TBranch        *b_pfTypeINoXYCorrmets_sign;   //!
   TBranch        *b_pfTypeINoXYCorrmets_sumEt;   //!
   TBranch        *b_pfTypeINoXYCorrmets_unCPhi;   //!
   TBranch        *b_pfTypeINoXYCorrmets_unCPt;   //!
   TBranch        *b_NpfTypeIType0mets;   //!
   TBranch        *b_pfTypeIType0mets_et;   //!
   TBranch        *b_pfTypeIType0mets_phi;   //!
   TBranch        *b_pfTypeIType0mets_ex;   //!
   TBranch        *b_pfTypeIType0mets_ey;   //!
   TBranch        *b_pfTypeIType0mets_gen_et;   //!
   TBranch        *b_pfTypeIType0mets_gen_phi;   //!
   TBranch        *b_pfTypeIType0mets_sign;   //!
   TBranch        *b_pfTypeIType0mets_sumEt;   //!
   TBranch        *b_pfTypeIType0mets_unCPhi;   //!
   TBranch        *b_pfTypeIType0mets_unCPt;   //!
   TBranch        *b_NpfTypeImets;   //!
   TBranch        *b_pfTypeImets_et;   //!
   TBranch        *b_pfTypeImets_phi;   //!
   TBranch        *b_pfTypeImets_ex;   //!
   TBranch        *b_pfTypeImets_ey;   //!
   TBranch        *b_pfTypeImets_gen_et;   //!
   TBranch        *b_pfTypeImets_gen_phi;   //!
   TBranch        *b_pfTypeImets_sign;   //!
   TBranch        *b_pfTypeImets_sumEt;   //!
   TBranch        *b_pfTypeImets_unCPhi;   //!
   TBranch        *b_pfTypeImets_unCPt;   //!
   TBranch        *b_Npf_els;   //!
   TBranch        *b_pf_els_energy;   //!
   TBranch        *b_pf_els_et;   //!
   TBranch        *b_pf_els_eta;   //!
   TBranch        *b_pf_els_phi;   //!
   TBranch        *b_pf_els_pt;   //!
   TBranch        *b_pf_els_px;   //!
   TBranch        *b_pf_els_py;   //!
   TBranch        *b_pf_els_pz;   //!
   TBranch        *b_pf_els_status;   //!
   TBranch        *b_pf_els_theta;   //!
   TBranch        *b_pf_els_gen_id;   //!
   TBranch        *b_pf_els_gen_phi;   //!
   TBranch        *b_pf_els_gen_pt;   //!
   TBranch        *b_pf_els_gen_pz;   //!
   TBranch        *b_pf_els_gen_px;   //!
   TBranch        *b_pf_els_gen_py;   //!
   TBranch        *b_pf_els_gen_eta;   //!
   TBranch        *b_pf_els_gen_theta;   //!
   TBranch        *b_pf_els_gen_et;   //!
   TBranch        *b_pf_els_gen_mother_id;   //!
   TBranch        *b_pf_els_gen_mother_phi;   //!
   TBranch        *b_pf_els_gen_mother_pt;   //!
   TBranch        *b_pf_els_gen_mother_pz;   //!
   TBranch        *b_pf_els_gen_mother_px;   //!
   TBranch        *b_pf_els_gen_mother_py;   //!
   TBranch        *b_pf_els_gen_mother_eta;   //!
   TBranch        *b_pf_els_gen_mother_theta;   //!
   TBranch        *b_pf_els_gen_mother_et;   //!
   TBranch        *b_pf_els_tightId;   //!
   TBranch        *b_pf_els_looseId;   //!
   TBranch        *b_pf_els_robustTightId;   //!
   TBranch        *b_pf_els_robustLooseId;   //!
   TBranch        *b_pf_els_robustHighEnergyId;   //!
   TBranch        *b_pf_els_simpleEleId95relIso;   //!
   TBranch        *b_pf_els_simpleEleId90relIso;   //!
   TBranch        *b_pf_els_simpleEleId85relIso;   //!
   TBranch        *b_pf_els_simpleEleId80relIso;   //!
   TBranch        *b_pf_els_simpleEleId70relIso;   //!
   TBranch        *b_pf_els_simpleEleId60relIso;   //!
   TBranch        *b_pf_els_simpleEleId95cIso;   //!
   TBranch        *b_pf_els_simpleEleId90cIso;   //!
   TBranch        *b_pf_els_simpleEleId85cIso;   //!
   TBranch        *b_pf_els_simpleEleId80cIso;   //!
   TBranch        *b_pf_els_simpleEleId70cIso;   //!
   TBranch        *b_pf_els_simpleEleId60cIso;   //!
   TBranch        *b_pf_els_cIso;   //!
   TBranch        *b_pf_els_tIso;   //!
   TBranch        *b_pf_els_ecalIso;   //!
   TBranch        *b_pf_els_hcalIso;   //!
   TBranch        *b_pf_els_chargedHadronIso;   //!
   TBranch        *b_pf_els_photonIso;   //!
   TBranch        *b_pf_els_neutralHadronIso;   //!
   TBranch        *b_pf_els_chi2;   //!
   TBranch        *b_pf_els_charge;   //!
   TBranch        *b_pf_els_caloEnergy;   //!
   TBranch        *b_pf_els_hadOverEm;   //!
   TBranch        *b_pf_els_hcalOverEcalBc;   //!
   TBranch        *b_pf_els_eOverPIn;   //!
   TBranch        *b_pf_els_eSeedOverPOut;   //!
   TBranch        *b_pf_els_sigmaEtaEta;   //!
   TBranch        *b_pf_els_sigmaIEtaIEta;   //!
   TBranch        *b_pf_els_scEnergy;   //!
   TBranch        *b_pf_els_scRawEnergy;   //!
   TBranch        *b_pf_els_scSeedEnergy;   //!
   TBranch        *b_pf_els_scEta;   //!
   TBranch        *b_pf_els_scPhi;   //!
   TBranch        *b_pf_els_scEtaWidth;   //!
   TBranch        *b_pf_els_scPhiWidth;   //!
   TBranch        *b_pf_els_scE1x5;   //!
   TBranch        *b_pf_els_scE2x5Max;   //!
   TBranch        *b_pf_els_scE5x5;   //!
   TBranch        *b_pf_els_isEB;   //!
   TBranch        *b_pf_els_isEE;   //!
   TBranch        *b_pf_els_dEtaIn;   //!
   TBranch        *b_pf_els_dPhiIn;   //!
   TBranch        *b_pf_els_dEtaOut;   //!
   TBranch        *b_pf_els_dPhiOut;   //!
   TBranch        *b_pf_els_numvalhits;   //!
   TBranch        *b_pf_els_numlosthits;   //!
   TBranch        *b_pf_els_basicClustersSize;   //!
   TBranch        *b_pf_els_tk_pz;   //!
   TBranch        *b_pf_els_tk_pt;   //!
   TBranch        *b_pf_els_tk_phi;   //!
   TBranch        *b_pf_els_tk_eta;   //!
   TBranch        *b_pf_els_d0dum;   //!
   TBranch        *b_pf_els_dz;   //!
   TBranch        *b_pf_els_vx;   //!
   TBranch        *b_pf_els_vy;   //!
   TBranch        *b_pf_els_vz;   //!
   TBranch        *b_pf_els_ndof;   //!
   TBranch        *b_pf_els_ptError;   //!
   TBranch        *b_pf_els_d0dumError;   //!
   TBranch        *b_pf_els_dzError;   //!
   TBranch        *b_pf_els_etaError;   //!
   TBranch        *b_pf_els_phiError;   //!
   TBranch        *b_pf_els_tk_charge;   //!
   TBranch        *b_pf_els_core_ecalDrivenSeed;   //!
   TBranch        *b_pf_els_n_inner_layer;   //!
   TBranch        *b_pf_els_n_outer_layer;   //!
   TBranch        *b_pf_els_ctf_tk_id;   //!
   TBranch        *b_pf_els_ctf_tk_charge;   //!
   TBranch        *b_pf_els_ctf_tk_eta;   //!
   TBranch        *b_pf_els_ctf_tk_phi;   //!
   TBranch        *b_pf_els_fbrem;   //!
   TBranch        *b_pf_els_shFracInnerHits;   //!
   TBranch        *b_pf_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_pf_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_pf_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_pf_els_dr03TkSumPt;   //!
   TBranch        *b_pf_els_dr04EcalRecHitSumEt;   //!
   TBranch        *b_pf_els_dr04HcalTowerSumEt;   //!
   TBranch        *b_pf_els_dr04HcalDepth1TowerSumEt;   //!
   TBranch        *b_pf_els_dr04HcalDepth2TowerSumEt;   //!
   TBranch        *b_pf_els_dr04TkSumPt;   //!
   TBranch        *b_pf_els_cpx;   //!
   TBranch        *b_pf_els_cpy;   //!
   TBranch        *b_pf_els_cpz;   //!
   TBranch        *b_pf_els_vpx;   //!
   TBranch        *b_pf_els_vpy;   //!
   TBranch        *b_pf_els_vpz;   //!
   TBranch        *b_pf_els_cx;   //!
   TBranch        *b_pf_els_cy;   //!
   TBranch        *b_pf_els_cz;   //!
   TBranch        *b_pf_els_PATpassConversionVeto;   //!
   TBranch        *b_Npf_mus;   //!
   TBranch        *b_pf_mus_energy;   //!
   TBranch        *b_pf_mus_et;   //!
   TBranch        *b_pf_mus_eta;   //!
   TBranch        *b_pf_mus_phi;   //!
   TBranch        *b_pf_mus_pt;   //!
   TBranch        *b_pf_mus_px;   //!
   TBranch        *b_pf_mus_py;   //!
   TBranch        *b_pf_mus_pz;   //!
   TBranch        *b_pf_mus_status;   //!
   TBranch        *b_pf_mus_theta;   //!
   TBranch        *b_pf_mus_gen_id;   //!
   TBranch        *b_pf_mus_gen_phi;   //!
   TBranch        *b_pf_mus_gen_pt;   //!
   TBranch        *b_pf_mus_gen_pz;   //!
   TBranch        *b_pf_mus_gen_px;   //!
   TBranch        *b_pf_mus_gen_py;   //!
   TBranch        *b_pf_mus_gen_eta;   //!
   TBranch        *b_pf_mus_gen_theta;   //!
   TBranch        *b_pf_mus_gen_et;   //!
   TBranch        *b_pf_mus_gen_mother_id;   //!
   TBranch        *b_pf_mus_gen_mother_phi;   //!
   TBranch        *b_pf_mus_gen_mother_pt;   //!
   TBranch        *b_pf_mus_gen_mother_pz;   //!
   TBranch        *b_pf_mus_gen_mother_px;   //!
   TBranch        *b_pf_mus_gen_mother_py;   //!
   TBranch        *b_pf_mus_gen_mother_eta;   //!
   TBranch        *b_pf_mus_gen_mother_theta;   //!
   TBranch        *b_pf_mus_gen_mother_et;   //!
   TBranch        *b_pf_mus_tkHits;   //!
   TBranch        *b_pf_mus_cIso;   //!
   TBranch        *b_pf_mus_tIso;   //!
   TBranch        *b_pf_mus_ecalIso;   //!
   TBranch        *b_pf_mus_hcalIso;   //!
   TBranch        *b_pf_mus_iso03_emVetoEt;   //!
   TBranch        *b_pf_mus_iso03_hadVetoEt;   //!
   TBranch        *b_pf_mus_calEnergyEm;   //!
   TBranch        *b_pf_mus_calEnergyHad;   //!
   TBranch        *b_pf_mus_calEnergyHo;   //!
   TBranch        *b_pf_mus_calEnergyEmS9;   //!
   TBranch        *b_pf_mus_calEnergyHadS9;   //!
   TBranch        *b_pf_mus_calEnergyHoS9;   //!
   TBranch        *b_pf_mus_iso03_sumPt;   //!
   TBranch        *b_pf_mus_iso03_emEt;   //!
   TBranch        *b_pf_mus_iso03_hadEt;   //!
   TBranch        *b_pf_mus_iso03_hoEt;   //!
   TBranch        *b_pf_mus_iso03_nTracks;   //!
   TBranch        *b_pf_mus_iso05_sumPt;   //!
   TBranch        *b_pf_mus_iso05_emEt;   //!
   TBranch        *b_pf_mus_iso05_hadEt;   //!
   TBranch        *b_pf_mus_iso05_hoEt;   //!
   TBranch        *b_pf_mus_iso05_nTracks;   //!
   TBranch        *b_pf_mus_neutralHadronIso;   //!
   TBranch        *b_pf_mus_chargedHadronIso;   //!
   TBranch        *b_pf_mus_photonIso;   //!
   TBranch        *b_pf_mus_charge;   //!
   TBranch        *b_pf_mus_cm_chi2;   //!
   TBranch        *b_pf_mus_cm_ndof;   //!
   TBranch        *b_pf_mus_cm_chg;   //!
   TBranch        *b_pf_mus_cm_pt;   //!
   TBranch        *b_pf_mus_cm_px;   //!
   TBranch        *b_pf_mus_cm_py;   //!
   TBranch        *b_pf_mus_cm_pz;   //!
   TBranch        *b_pf_mus_cm_eta;   //!
   TBranch        *b_pf_mus_cm_phi;   //!
   TBranch        *b_pf_mus_cm_theta;   //!
   TBranch        *b_pf_mus_cm_d0dum;   //!
   TBranch        *b_pf_mus_cm_dz;   //!
   TBranch        *b_pf_mus_cm_vx;   //!
   TBranch        *b_pf_mus_cm_vy;   //!
   TBranch        *b_pf_mus_cm_vz;   //!
   TBranch        *b_pf_mus_cm_numvalhits;   //!
   TBranch        *b_pf_mus_cm_numlosthits;   //!
   TBranch        *b_pf_mus_cm_numvalMuonhits;   //!
   TBranch        *b_pf_mus_cm_d0dumErr;   //!
   TBranch        *b_pf_mus_cm_dzErr;   //!
   TBranch        *b_pf_mus_cm_ptErr;   //!
   TBranch        *b_pf_mus_cm_etaErr;   //!
   TBranch        *b_pf_mus_cm_phiErr;   //!
   TBranch        *b_pf_mus_tk_id;   //!
   TBranch        *b_pf_mus_tk_chi2;   //!
   TBranch        *b_pf_mus_tk_ndof;   //!
   TBranch        *b_pf_mus_tk_chg;   //!
   TBranch        *b_pf_mus_tk_pt;   //!
   TBranch        *b_pf_mus_tk_px;   //!
   TBranch        *b_pf_mus_tk_py;   //!
   TBranch        *b_pf_mus_tk_pz;   //!
   TBranch        *b_pf_mus_tk_eta;   //!
   TBranch        *b_pf_mus_tk_phi;   //!
   TBranch        *b_pf_mus_tk_theta;   //!
   TBranch        *b_pf_mus_tk_d0dum;   //!
   TBranch        *b_pf_mus_tk_dz;   //!
   TBranch        *b_pf_mus_tk_vx;   //!
   TBranch        *b_pf_mus_tk_vy;   //!
   TBranch        *b_pf_mus_tk_vz;   //!
   TBranch        *b_pf_mus_tk_numvalhits;   //!
   TBranch        *b_pf_mus_tk_numlosthits;   //!
   TBranch        *b_pf_mus_tk_d0dumErr;   //!
   TBranch        *b_pf_mus_tk_dzErr;   //!
   TBranch        *b_pf_mus_tk_ptErr;   //!
   TBranch        *b_pf_mus_tk_etaErr;   //!
   TBranch        *b_pf_mus_tk_phiErr;   //!
   TBranch        *b_pf_mus_tk_numvalPixelhits;   //!
   TBranch        *b_pf_mus_tk_numpixelWthMeasr;   //!
   TBranch        *b_pf_mus_stamu_chi2;   //!
   TBranch        *b_pf_mus_stamu_ndof;   //!
   TBranch        *b_pf_mus_stamu_chg;   //!
   TBranch        *b_pf_mus_stamu_pt;   //!
   TBranch        *b_pf_mus_stamu_px;   //!
   TBranch        *b_pf_mus_stamu_py;   //!
   TBranch        *b_pf_mus_stamu_pz;   //!
   TBranch        *b_pf_mus_stamu_eta;   //!
   TBranch        *b_pf_mus_stamu_phi;   //!
   TBranch        *b_pf_mus_stamu_theta;   //!
   TBranch        *b_pf_mus_stamu_d0dum;   //!
   TBranch        *b_pf_mus_stamu_dz;   //!
   TBranch        *b_pf_mus_stamu_vx;   //!
   TBranch        *b_pf_mus_stamu_vy;   //!
   TBranch        *b_pf_mus_stamu_vz;   //!
   TBranch        *b_pf_mus_stamu_numvalhits;   //!
   TBranch        *b_pf_mus_stamu_numlosthits;   //!
   TBranch        *b_pf_mus_stamu_d0dumErr;   //!
   TBranch        *b_pf_mus_stamu_dzErr;   //!
   TBranch        *b_pf_mus_stamu_ptErr;   //!
   TBranch        *b_pf_mus_stamu_etaErr;   //!
   TBranch        *b_pf_mus_stamu_phiErr;   //!
   TBranch        *b_pf_mus_num_matches;   //!
   TBranch        *b_pf_mus_isTrackerMuon;   //!
   TBranch        *b_pf_mus_isStandAloneMuon;   //!
   TBranch        *b_pf_mus_isCaloMuon;   //!
   TBranch        *b_pf_mus_isGlobalMuon;   //!
   TBranch        *b_pf_mus_isElectron;   //!
   TBranch        *b_pf_mus_isConvertedPhoton;   //!
   TBranch        *b_pf_mus_isPhoton;   //!
   TBranch        *b_pf_mus_id_All;   //!
   TBranch        *b_pf_mus_id_AllGlobalMuons;   //!
   TBranch        *b_pf_mus_id_AllStandAloneMuons;   //!
   TBranch        *b_pf_mus_id_AllTrackerMuons;   //!
   TBranch        *b_pf_mus_id_TrackerMuonArbitrated;   //!
   TBranch        *b_pf_mus_id_AllArbitrated;   //!
   TBranch        *b_pf_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_pf_mus_id_TMLastStationLoose;   //!
   TBranch        *b_pf_mus_id_TMLastStationTight;   //!
   TBranch        *b_pf_mus_id_TM2DCompatibilityLoose;   //!
   TBranch        *b_pf_mus_id_TM2DCompatibilityTight;   //!
   TBranch        *b_pf_mus_id_TMOneStationLoose;   //!
   TBranch        *b_pf_mus_id_TMOneStationTight;   //!
   TBranch        *b_pf_mus_id_TMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_pf_mus_id_TMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_pf_mus_tk_LayersWithMeasurement;   //!
   TBranch        *b_pf_mus_tk_PixelLayersWithMeasurement;   //!
   TBranch        *b_pf_mus_tk_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_pf_mus_tk_LayersWithoutMeasurement;   //!
   TBranch        *b_pf_mus_tk_ExpectedHitsInner;   //!
   TBranch        *b_pf_mus_tk_ExpectedHitsOuter;   //!
   TBranch        *b_pf_mus_cm_LayersWithMeasurement;   //!
   TBranch        *b_pf_mus_cm_PixelLayersWithMeasurement;   //!
   TBranch        *b_pf_mus_cm_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_pf_mus_cm_LayersWithoutMeasurement;   //!
   TBranch        *b_pf_mus_cm_ExpectedHitsInner;   //!
   TBranch        *b_pf_mus_cm_ExpectedHitsOuter;   //!
   TBranch        *b_pf_mus_picky_LayersWithMeasurement;   //!
   TBranch        *b_pf_mus_picky_PixelLayersWithMeasurement;   //!
   TBranch        *b_pf_mus_picky_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_pf_mus_picky_LayersWithoutMeasurement;   //!
   TBranch        *b_pf_mus_picky_ExpectedHitsInner;   //!
   TBranch        *b_pf_mus_picky_ExpectedHitsOuter;   //!
   TBranch        *b_pf_mus_tpfms_LayersWithMeasurement;   //!
   TBranch        *b_pf_mus_tpfms_PixelLayersWithMeasurement;   //!
   TBranch        *b_pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit;   //!
   TBranch        *b_pf_mus_tpfms_LayersWithoutMeasurement;   //!
   TBranch        *b_pf_mus_tpfms_ExpectedHitsInner;   //!
   TBranch        *b_pf_mus_tpfms_ExpectedHitsOuter;   //!
   TBranch        *b_pf_mus_picky_id;   //!
   TBranch        *b_pf_mus_picky_chi2;   //!
   TBranch        *b_pf_mus_picky_ndof;   //!
   TBranch        *b_pf_mus_picky_chg;   //!
   TBranch        *b_pf_mus_picky_pt;   //!
   TBranch        *b_pf_mus_picky_px;   //!
   TBranch        *b_pf_mus_picky_py;   //!
   TBranch        *b_pf_mus_picky_pz;   //!
   TBranch        *b_pf_mus_picky_eta;   //!
   TBranch        *b_pf_mus_picky_phi;   //!
   TBranch        *b_pf_mus_picky_theta;   //!
   TBranch        *b_pf_mus_picky_d0dum;   //!
   TBranch        *b_pf_mus_picky_dz;   //!
   TBranch        *b_pf_mus_picky_vx;   //!
   TBranch        *b_pf_mus_picky_vy;   //!
   TBranch        *b_pf_mus_picky_vz;   //!
   TBranch        *b_pf_mus_picky_numvalhits;   //!
   TBranch        *b_pf_mus_picky_numlosthits;   //!
   TBranch        *b_pf_mus_picky_d0dumErr;   //!
   TBranch        *b_pf_mus_picky_dzErr;   //!
   TBranch        *b_pf_mus_picky_ptErr;   //!
   TBranch        *b_pf_mus_picky_etaErr;   //!
   TBranch        *b_pf_mus_picky_phiErr;   //!
   TBranch        *b_pf_mus_picky_numvalPixelhits;   //!
   TBranch        *b_pf_mus_tpfms_id;   //!
   TBranch        *b_pf_mus_tpfms_chi2;   //!
   TBranch        *b_pf_mus_tpfms_ndof;   //!
   TBranch        *b_pf_mus_tpfms_chg;   //!
   TBranch        *b_pf_mus_tpfms_pt;   //!
   TBranch        *b_pf_mus_tpfms_px;   //!
   TBranch        *b_pf_mus_tpfms_py;   //!
   TBranch        *b_pf_mus_tpfms_pz;   //!
   TBranch        *b_pf_mus_tpfms_eta;   //!
   TBranch        *b_pf_mus_tpfms_phi;   //!
   TBranch        *b_pf_mus_tpfms_theta;   //!
   TBranch        *b_pf_mus_tpfms_d0dum;   //!
   TBranch        *b_pf_mus_tpfms_dz;   //!
   TBranch        *b_pf_mus_tpfms_vx;   //!
   TBranch        *b_pf_mus_tpfms_vy;   //!
   TBranch        *b_pf_mus_tpfms_vz;   //!
   TBranch        *b_pf_mus_tpfms_numvalhits;   //!
   TBranch        *b_pf_mus_tpfms_numlosthits;   //!
   TBranch        *b_pf_mus_tpfms_d0dumErr;   //!
   TBranch        *b_pf_mus_tpfms_dzErr;   //!
   TBranch        *b_pf_mus_tpfms_ptErr;   //!
   TBranch        *b_pf_mus_tpfms_etaErr;   //!
   TBranch        *b_pf_mus_tpfms_phiErr;   //!
   TBranch        *b_pf_mus_tpfms_numvalPixelhits;   //!
   TBranch        *b_pf_mus_pfIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_pf_mus_pfIsolationR03_sumChargedParticlePt;   //!
   TBranch        *b_pf_mus_pfIsolationR03_sumNeutralHadronEt;   //!
   TBranch        *b_pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_pf_mus_pfIsolationR03_sumPhotonEt;   //!
   TBranch        *b_pf_mus_pfIsolationR03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_pf_mus_pfIsolationR03_sumPUPt;   //!
   TBranch        *b_pf_mus_pfIsolationR04_sumChargedHadronPt;   //!
   TBranch        *b_pf_mus_pfIsolationR04_sumChargedParticlePt;   //!
   TBranch        *b_pf_mus_pfIsolationR04_sumNeutralHadronEt;   //!
   TBranch        *b_pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_pf_mus_pfIsolationR04_sumPhotonEt;   //!
   TBranch        *b_pf_mus_pfIsolationR04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_pf_mus_pfIsolationR04_sumPUPt;   //!
   TBranch        *b_pf_mus_dB;   //!
   TBranch        *b_pf_mus_numberOfMatchedStations;   //!
   TBranch        *b_pf_mus_isPFMuon;   //!
   TBranch        *b_Npf_photons;   //!
   TBranch        *b_pf_photons_energy;   //!
   TBranch        *b_pf_photons_et;   //!
   TBranch        *b_pf_photons_eta;   //!
   TBranch        *b_pf_photons_phi;   //!
   TBranch        *b_pf_photons_pt;   //!
   TBranch        *b_pf_photons_px;   //!
   TBranch        *b_pf_photons_py;   //!
   TBranch        *b_pf_photons_pz;   //!
   TBranch        *b_pf_photons_status;   //!
   TBranch        *b_pf_photons_theta;   //!
   TBranch        *b_pf_photons_hadOverEM;   //!
   TBranch        *b_pf_photons_hadTowOverEM;   //!
   TBranch        *b_pf_photons_scEnergy;   //!
   TBranch        *b_pf_photons_scRawEnergy;   //!
   TBranch        *b_pf_photons_scEta;   //!
   TBranch        *b_pf_photons_scPhi;   //!
   TBranch        *b_pf_photons_scEtaWidth;   //!
   TBranch        *b_pf_photons_scPhiWidth;   //!
   TBranch        *b_pf_photons_isAlsoElectron;   //!
   TBranch        *b_pf_photons_hasPixelSeed;   //!
   TBranch        *b_pf_photons_isConverted;   //!
   TBranch        *b_pf_photons_isEBGap;   //!
   TBranch        *b_pf_photons_isEEGap;   //!
   TBranch        *b_pf_photons_isEBEEGap;   //!
   TBranch        *b_pf_photons_isEBPho;   //!
   TBranch        *b_pf_photons_isEEPho;   //!
   TBranch        *b_pf_photons_maxEnergyXtal;   //!
   TBranch        *b_pf_photons_e1x5;   //!
   TBranch        *b_pf_photons_e2x5;   //!
   TBranch        *b_pf_photons_e3x3;   //!
   TBranch        *b_pf_photons_e5x5;   //!
   TBranch        *b_pf_photons_sigmaEtaEta;   //!
   TBranch        *b_pf_photons_sigmaIetaIeta;   //!
   TBranch        *b_pf_photons_r9;   //!
   TBranch        *b_pf_photons_chIso;   //!
   TBranch        *b_pf_photons_nhIso;   //!
   TBranch        *b_pf_photons_phIso;   //!
   TBranch        *b_Npfcand;   //!
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
   TBranch        *b_Npfmets;   //!
   TBranch        *b_pfmets_et;   //!
   TBranch        *b_pfmets_phi;   //!
   TBranch        *b_pfmets_ex;   //!
   TBranch        *b_pfmets_ey;   //!
   TBranch        *b_pfmets_gen_et;   //!
   TBranch        *b_pfmets_gen_phi;   //!
   TBranch        *b_pfmets_sign;   //!
   TBranch        *b_pfmets_sumEt;   //!
   TBranch        *b_pfmets_unCPhi;   //!
   TBranch        *b_pfmets_unCPt;   //!
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
   TBranch        *b_Ntcmets;   //!
   TBranch        *b_tcmets_et;   //!
   TBranch        *b_tcmets_phi;   //!
   TBranch        *b_tcmets_ex;   //!
   TBranch        *b_tcmets_ey;   //!
   TBranch        *b_tcmets_sumEt;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_tracks_chi2;   //!
   TBranch        *b_tracks_ndof;   //!
   TBranch        *b_tracks_chg;   //!
   TBranch        *b_tracks_pt;   //!
   TBranch        *b_tracks_px;   //!
   TBranch        *b_tracks_py;   //!
   TBranch        *b_tracks_pz;   //!
   TBranch        *b_tracks_eta;   //!
   TBranch        *b_tracks_phi;   //!
   TBranch        *b_tracks_d0dum;   //!
   TBranch        *b_tracks_dz;   //!
   TBranch        *b_tracks_vx;   //!
   TBranch        *b_tracks_vy;   //!
   TBranch        *b_tracks_vz;   //!
   TBranch        *b_tracks_numvalhits;   //!
   TBranch        *b_tracks_numlosthits;   //!
   TBranch        *b_tracks_d0dumErr;   //!
   TBranch        *b_tracks_dzErr;   //!
   TBranch        *b_tracks_ptErr;   //!
   TBranch        *b_tracks_etaErr;   //!
   TBranch        *b_tracks_phiErr;   //!
   TBranch        *b_tracks_highPurity;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_experimentType;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
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
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_px;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_py;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_pz;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_energy;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_phi;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_eta;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_index;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT10_nconstituents;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_px;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_py;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_pz;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_energy;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_phi;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_eta;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_index;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT15_nconstituents;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_px;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_py;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_pz;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_energy;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_phi;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_eta;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_index;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT20_nconstituents;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_px;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_py;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_pz;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_energy;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_phi;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_eta;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_index;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT25_nconstituents;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_px;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_py;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_pz;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_energy;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_phi;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_eta;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_index;   //!
   TBranch        *b_fastjets_AK5PF_R1p2_R0p5pT30_nconstituents;   //!

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
   jets_AK5PF_status = 0;
   jets_AK5PF_phi = 0;
   jets_AK5PF_pt = 0;
   jets_AK5PF_pz = 0;
   jets_AK5PF_px = 0;
   jets_AK5PF_py = 0;
   jets_AK5PF_eta = 0;
   jets_AK5PF_theta = 0;
   jets_AK5PF_et = 0;
   jets_AK5PF_energy = 0;
   jets_AK5PF_parton_Id = 0;
   jets_AK5PF_parton_motherId = 0;
   jets_AK5PF_parton_pt = 0;
   jets_AK5PF_parton_phi = 0;
   jets_AK5PF_parton_eta = 0;
   jets_AK5PF_parton_Energy = 0;
   jets_AK5PF_parton_mass = 0;
   jets_AK5PF_gen_et = 0;
   jets_AK5PF_gen_pt = 0;
   jets_AK5PF_gen_eta = 0;
   jets_AK5PF_gen_phi = 0;
   jets_AK5PF_gen_mass = 0;
   jets_AK5PF_gen_Energy = 0;
   jets_AK5PF_gen_Id = 0;
   jets_AK5PF_gen_motherID = 0;
   jets_AK5PF_gen_threeCharge = 0;
   jets_AK5PF_partonFlavour = 0;
   jets_AK5PF_btag_TC_highPur = 0;
   jets_AK5PF_btag_TC_highEff = 0;
   jets_AK5PF_btag_jetProb = 0;
   jets_AK5PF_btag_jetBProb = 0;
   jets_AK5PF_btag_softEle = 0;
   jets_AK5PF_btag_softMuon = 0;
   jets_AK5PF_btag_secVertexHighPur = 0;
   jets_AK5PF_btag_secVertexHighEff = 0;
   jets_AK5PF_btag_secVertexCombined = 0;
   jets_AK5PF_jetCharge = 0;
   jets_AK5PF_chgEmE = 0;
   jets_AK5PF_chgHadE = 0;
   jets_AK5PF_photonEnergy = 0;
   jets_AK5PF_chgMuE = 0;
   jets_AK5PF_chg_Mult = 0;
   jets_AK5PF_neutralEmE = 0;
   jets_AK5PF_neutralHadE = 0;
   jets_AK5PF_neutral_Mult = 0;
   jets_AK5PF_mu_Mult = 0;
   jets_AK5PF_emf = 0;
   jets_AK5PF_ehf = 0;
   jets_AK5PF_n60 = 0;
   jets_AK5PF_n90 = 0;
   jets_AK5PF_etaetaMoment = 0;
   jets_AK5PF_etaphiMoment = 0;
   jets_AK5PF_phiphiMoment = 0;
   jets_AK5PF_n90Hits = 0;
   jets_AK5PF_fHPD = 0;
   jets_AK5PF_fRBX = 0;
   jets_AK5PF_hitsInN90 = 0;
   jets_AK5PF_nECALTowers = 0;
   jets_AK5PF_nHCALTowers = 0;
   jets_AK5PF_fSubDetector1 = 0;
   jets_AK5PF_fSubDetector2 = 0;
   jets_AK5PF_fSubDetector3 = 0;
   jets_AK5PF_fSubDetector4 = 0;
   jets_AK5PF_area = 0;
   jets_AK5PF_corrFactorRaw = 0;
   jets_AK5PF_rawPt = 0;
   jets_AK5PF_mass = 0;
   jets_AK5PFclean_status = 0;
   jets_AK5PFclean_phi = 0;
   jets_AK5PFclean_pt = 0;
   jets_AK5PFclean_pz = 0;
   jets_AK5PFclean_px = 0;
   jets_AK5PFclean_py = 0;
   jets_AK5PFclean_eta = 0;
   jets_AK5PFclean_theta = 0;
   jets_AK5PFclean_et = 0;
   jets_AK5PFclean_energy = 0;
   jets_AK5PFclean_parton_Id = 0;
   jets_AK5PFclean_parton_motherId = 0;
   jets_AK5PFclean_parton_pt = 0;
   jets_AK5PFclean_parton_phi = 0;
   jets_AK5PFclean_parton_eta = 0;
   jets_AK5PFclean_parton_Energy = 0;
   jets_AK5PFclean_parton_mass = 0;
   jets_AK5PFclean_gen_et = 0;
   jets_AK5PFclean_gen_pt = 0;
   jets_AK5PFclean_gen_eta = 0;
   jets_AK5PFclean_gen_phi = 0;
   jets_AK5PFclean_gen_mass = 0;
   jets_AK5PFclean_gen_Energy = 0;
   jets_AK5PFclean_gen_Id = 0;
   jets_AK5PFclean_partonFlavour = 0;
   jets_AK5PFclean_btag_TC_highPur = 0;
   jets_AK5PFclean_btag_TC_highEff = 0;
   jets_AK5PFclean_btag_jetProb = 0;
   jets_AK5PFclean_btag_jetBProb = 0;
   jets_AK5PFclean_btag_softEle = 0;
   jets_AK5PFclean_btag_softMuon = 0;
   jets_AK5PFclean_btag_secVertexHighPur = 0;
   jets_AK5PFclean_btag_secVertexHighEff = 0;
   jets_AK5PFclean_btag_secVertexCombined = 0;
   jets_AK5PFclean_jetCharge = 0;
   jets_AK5PFclean_chgEmE = 0;
   jets_AK5PFclean_chgHadE = 0;
   jets_AK5PFclean_photonEnergy = 0;
   jets_AK5PFclean_chgMuE = 0;
   jets_AK5PFclean_chg_Mult = 0;
   jets_AK5PFclean_neutralEmE = 0;
   jets_AK5PFclean_neutralHadE = 0;
   jets_AK5PFclean_neutral_Mult = 0;
   jets_AK5PFclean_mu_Mult = 0;
   jets_AK5PFclean_emf = 0;
   jets_AK5PFclean_ehf = 0;
   jets_AK5PFclean_n60 = 0;
   jets_AK5PFclean_n90 = 0;
   jets_AK5PFclean_etaetaMoment = 0;
   jets_AK5PFclean_etaphiMoment = 0;
   jets_AK5PFclean_phiphiMoment = 0;
   jets_AK5PFclean_n90Hits = 0;
   jets_AK5PFclean_fHPD = 0;
   jets_AK5PFclean_fRBX = 0;
   jets_AK5PFclean_hitsInN90 = 0;
   jets_AK5PFclean_nECALTowers = 0;
   jets_AK5PFclean_nHCALTowers = 0;
   jets_AK5PFclean_fSubDetector1 = 0;
   jets_AK5PFclean_fSubDetector2 = 0;
   jets_AK5PFclean_fSubDetector3 = 0;
   jets_AK5PFclean_fSubDetector4 = 0;
   jets_AK5PFclean_area = 0;
   jets_AK5PFclean_corrFactorRaw = 0;
   jets_AK5PFclean_rawPt = 0;
   jets_AK5PFclean_mass = 0;
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
   mc_pdf_x1 = 0;
   mc_pdf_x2 = 0;
   mc_pdf_q = 0;
   mc_pdf_id1 = 0;
   mc_pdf_id2 = 0;
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
   metsHO_et = 0;
   metsHO_phi = 0;
   metsHO_ex = 0;
   metsHO_ey = 0;
   metsHO_sumEt = 0;
   mets_AK5_et = 0;
   mets_AK5_phi = 0;
   mets_AK5_ex = 0;
   mets_AK5_ey = 0;
   mets_AK5_gen_et = 0;
   mets_AK5_gen_phi = 0;
   mets_AK5_sign = 0;
   mets_AK5_sumEt = 0;
   mets_AK5_unCPhi = 0;
   mets_AK5_unCPt = 0;
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
   pfTypeINoXYCorrmets_et = 0;
   pfTypeINoXYCorrmets_phi = 0;
   pfTypeINoXYCorrmets_ex = 0;
   pfTypeINoXYCorrmets_ey = 0;
   pfTypeINoXYCorrmets_gen_et = 0;
   pfTypeINoXYCorrmets_gen_phi = 0;
   pfTypeINoXYCorrmets_sign = 0;
   pfTypeINoXYCorrmets_sumEt = 0;
   pfTypeINoXYCorrmets_unCPhi = 0;
   pfTypeINoXYCorrmets_unCPt = 0;
   pfTypeIType0mets_et = 0;
   pfTypeIType0mets_phi = 0;
   pfTypeIType0mets_ex = 0;
   pfTypeIType0mets_ey = 0;
   pfTypeIType0mets_gen_et = 0;
   pfTypeIType0mets_gen_phi = 0;
   pfTypeIType0mets_sign = 0;
   pfTypeIType0mets_sumEt = 0;
   pfTypeIType0mets_unCPhi = 0;
   pfTypeIType0mets_unCPt = 0;
   pfTypeImets_et = 0;
   pfTypeImets_phi = 0;
   pfTypeImets_ex = 0;
   pfTypeImets_ey = 0;
   pfTypeImets_gen_et = 0;
   pfTypeImets_gen_phi = 0;
   pfTypeImets_sign = 0;
   pfTypeImets_sumEt = 0;
   pfTypeImets_unCPhi = 0;
   pfTypeImets_unCPt = 0;
   pf_els_energy = 0;
   pf_els_et = 0;
   pf_els_eta = 0;
   pf_els_phi = 0;
   pf_els_pt = 0;
   pf_els_px = 0;
   pf_els_py = 0;
   pf_els_pz = 0;
   pf_els_status = 0;
   pf_els_theta = 0;
   pf_els_gen_id = 0;
   pf_els_gen_phi = 0;
   pf_els_gen_pt = 0;
   pf_els_gen_pz = 0;
   pf_els_gen_px = 0;
   pf_els_gen_py = 0;
   pf_els_gen_eta = 0;
   pf_els_gen_theta = 0;
   pf_els_gen_et = 0;
   pf_els_gen_mother_id = 0;
   pf_els_gen_mother_phi = 0;
   pf_els_gen_mother_pt = 0;
   pf_els_gen_mother_pz = 0;
   pf_els_gen_mother_px = 0;
   pf_els_gen_mother_py = 0;
   pf_els_gen_mother_eta = 0;
   pf_els_gen_mother_theta = 0;
   pf_els_gen_mother_et = 0;
   pf_els_tightId = 0;
   pf_els_looseId = 0;
   pf_els_robustTightId = 0;
   pf_els_robustLooseId = 0;
   pf_els_robustHighEnergyId = 0;
   pf_els_simpleEleId95relIso = 0;
   pf_els_simpleEleId90relIso = 0;
   pf_els_simpleEleId85relIso = 0;
   pf_els_simpleEleId80relIso = 0;
   pf_els_simpleEleId70relIso = 0;
   pf_els_simpleEleId60relIso = 0;
   pf_els_simpleEleId95cIso = 0;
   pf_els_simpleEleId90cIso = 0;
   pf_els_simpleEleId85cIso = 0;
   pf_els_simpleEleId80cIso = 0;
   pf_els_simpleEleId70cIso = 0;
   pf_els_simpleEleId60cIso = 0;
   pf_els_cIso = 0;
   pf_els_tIso = 0;
   pf_els_ecalIso = 0;
   pf_els_hcalIso = 0;
   pf_els_chargedHadronIso = 0;
   pf_els_photonIso = 0;
   pf_els_neutralHadronIso = 0;
   pf_els_chi2 = 0;
   pf_els_charge = 0;
   pf_els_caloEnergy = 0;
   pf_els_hadOverEm = 0;
   pf_els_hcalOverEcalBc = 0;
   pf_els_eOverPIn = 0;
   pf_els_eSeedOverPOut = 0;
   pf_els_sigmaEtaEta = 0;
   pf_els_sigmaIEtaIEta = 0;
   pf_els_scEnergy = 0;
   pf_els_scRawEnergy = 0;
   pf_els_scSeedEnergy = 0;
   pf_els_scEta = 0;
   pf_els_scPhi = 0;
   pf_els_scEtaWidth = 0;
   pf_els_scPhiWidth = 0;
   pf_els_scE1x5 = 0;
   pf_els_scE2x5Max = 0;
   pf_els_scE5x5 = 0;
   pf_els_isEB = 0;
   pf_els_isEE = 0;
   pf_els_dEtaIn = 0;
   pf_els_dPhiIn = 0;
   pf_els_dEtaOut = 0;
   pf_els_dPhiOut = 0;
   pf_els_numvalhits = 0;
   pf_els_numlosthits = 0;
   pf_els_basicClustersSize = 0;
   pf_els_tk_pz = 0;
   pf_els_tk_pt = 0;
   pf_els_tk_phi = 0;
   pf_els_tk_eta = 0;
   pf_els_d0dum = 0;
   pf_els_dz = 0;
   pf_els_vx = 0;
   pf_els_vy = 0;
   pf_els_vz = 0;
   pf_els_ndof = 0;
   pf_els_ptError = 0;
   pf_els_d0dumError = 0;
   pf_els_dzError = 0;
   pf_els_etaError = 0;
   pf_els_phiError = 0;
   pf_els_tk_charge = 0;
   pf_els_core_ecalDrivenSeed = 0;
   pf_els_n_inner_layer = 0;
   pf_els_n_outer_layer = 0;
   pf_els_ctf_tk_id = 0;
   pf_els_ctf_tk_charge = 0;
   pf_els_ctf_tk_eta = 0;
   pf_els_ctf_tk_phi = 0;
   pf_els_fbrem = 0;
   pf_els_shFracInnerHits = 0;
   pf_els_dr03EcalRecHitSumEt = 0;
   pf_els_dr03HcalTowerSumEt = 0;
   pf_els_dr03HcalDepth1TowerSumEt = 0;
   pf_els_dr03HcalDepth2TowerSumEt = 0;
   pf_els_dr03TkSumPt = 0;
   pf_els_dr04EcalRecHitSumEt = 0;
   pf_els_dr04HcalTowerSumEt = 0;
   pf_els_dr04HcalDepth1TowerSumEt = 0;
   pf_els_dr04HcalDepth2TowerSumEt = 0;
   pf_els_dr04TkSumPt = 0;
   pf_els_cpx = 0;
   pf_els_cpy = 0;
   pf_els_cpz = 0;
   pf_els_vpx = 0;
   pf_els_vpy = 0;
   pf_els_vpz = 0;
   pf_els_cx = 0;
   pf_els_cy = 0;
   pf_els_cz = 0;
   pf_els_PATpassConversionVeto = 0;
   pf_mus_energy = 0;
   pf_mus_et = 0;
   pf_mus_eta = 0;
   pf_mus_phi = 0;
   pf_mus_pt = 0;
   pf_mus_px = 0;
   pf_mus_py = 0;
   pf_mus_pz = 0;
   pf_mus_status = 0;
   pf_mus_theta = 0;
   pf_mus_gen_id = 0;
   pf_mus_gen_phi = 0;
   pf_mus_gen_pt = 0;
   pf_mus_gen_pz = 0;
   pf_mus_gen_px = 0;
   pf_mus_gen_py = 0;
   pf_mus_gen_eta = 0;
   pf_mus_gen_theta = 0;
   pf_mus_gen_et = 0;
   pf_mus_gen_mother_id = 0;
   pf_mus_gen_mother_phi = 0;
   pf_mus_gen_mother_pt = 0;
   pf_mus_gen_mother_pz = 0;
   pf_mus_gen_mother_px = 0;
   pf_mus_gen_mother_py = 0;
   pf_mus_gen_mother_eta = 0;
   pf_mus_gen_mother_theta = 0;
   pf_mus_gen_mother_et = 0;
   pf_mus_tkHits = 0;
   pf_mus_cIso = 0;
   pf_mus_tIso = 0;
   pf_mus_ecalIso = 0;
   pf_mus_hcalIso = 0;
   pf_mus_iso03_emVetoEt = 0;
   pf_mus_iso03_hadVetoEt = 0;
   pf_mus_calEnergyEm = 0;
   pf_mus_calEnergyHad = 0;
   pf_mus_calEnergyHo = 0;
   pf_mus_calEnergyEmS9 = 0;
   pf_mus_calEnergyHadS9 = 0;
   pf_mus_calEnergyHoS9 = 0;
   pf_mus_iso03_sumPt = 0;
   pf_mus_iso03_emEt = 0;
   pf_mus_iso03_hadEt = 0;
   pf_mus_iso03_hoEt = 0;
   pf_mus_iso03_nTracks = 0;
   pf_mus_iso05_sumPt = 0;
   pf_mus_iso05_emEt = 0;
   pf_mus_iso05_hadEt = 0;
   pf_mus_iso05_hoEt = 0;
   pf_mus_iso05_nTracks = 0;
   pf_mus_neutralHadronIso = 0;
   pf_mus_chargedHadronIso = 0;
   pf_mus_photonIso = 0;
   pf_mus_charge = 0;
   pf_mus_cm_chi2 = 0;
   pf_mus_cm_ndof = 0;
   pf_mus_cm_chg = 0;
   pf_mus_cm_pt = 0;
   pf_mus_cm_px = 0;
   pf_mus_cm_py = 0;
   pf_mus_cm_pz = 0;
   pf_mus_cm_eta = 0;
   pf_mus_cm_phi = 0;
   pf_mus_cm_theta = 0;
   pf_mus_cm_d0dum = 0;
   pf_mus_cm_dz = 0;
   pf_mus_cm_vx = 0;
   pf_mus_cm_vy = 0;
   pf_mus_cm_vz = 0;
   pf_mus_cm_numvalhits = 0;
   pf_mus_cm_numlosthits = 0;
   pf_mus_cm_numvalMuonhits = 0;
   pf_mus_cm_d0dumErr = 0;
   pf_mus_cm_dzErr = 0;
   pf_mus_cm_ptErr = 0;
   pf_mus_cm_etaErr = 0;
   pf_mus_cm_phiErr = 0;
   pf_mus_tk_id = 0;
   pf_mus_tk_chi2 = 0;
   pf_mus_tk_ndof = 0;
   pf_mus_tk_chg = 0;
   pf_mus_tk_pt = 0;
   pf_mus_tk_px = 0;
   pf_mus_tk_py = 0;
   pf_mus_tk_pz = 0;
   pf_mus_tk_eta = 0;
   pf_mus_tk_phi = 0;
   pf_mus_tk_theta = 0;
   pf_mus_tk_d0dum = 0;
   pf_mus_tk_dz = 0;
   pf_mus_tk_vx = 0;
   pf_mus_tk_vy = 0;
   pf_mus_tk_vz = 0;
   pf_mus_tk_numvalhits = 0;
   pf_mus_tk_numlosthits = 0;
   pf_mus_tk_d0dumErr = 0;
   pf_mus_tk_dzErr = 0;
   pf_mus_tk_ptErr = 0;
   pf_mus_tk_etaErr = 0;
   pf_mus_tk_phiErr = 0;
   pf_mus_tk_numvalPixelhits = 0;
   pf_mus_tk_numpixelWthMeasr = 0;
   pf_mus_stamu_chi2 = 0;
   pf_mus_stamu_ndof = 0;
   pf_mus_stamu_chg = 0;
   pf_mus_stamu_pt = 0;
   pf_mus_stamu_px = 0;
   pf_mus_stamu_py = 0;
   pf_mus_stamu_pz = 0;
   pf_mus_stamu_eta = 0;
   pf_mus_stamu_phi = 0;
   pf_mus_stamu_theta = 0;
   pf_mus_stamu_d0dum = 0;
   pf_mus_stamu_dz = 0;
   pf_mus_stamu_vx = 0;
   pf_mus_stamu_vy = 0;
   pf_mus_stamu_vz = 0;
   pf_mus_stamu_numvalhits = 0;
   pf_mus_stamu_numlosthits = 0;
   pf_mus_stamu_d0dumErr = 0;
   pf_mus_stamu_dzErr = 0;
   pf_mus_stamu_ptErr = 0;
   pf_mus_stamu_etaErr = 0;
   pf_mus_stamu_phiErr = 0;
   pf_mus_num_matches = 0;
   pf_mus_isTrackerMuon = 0;
   pf_mus_isStandAloneMuon = 0;
   pf_mus_isCaloMuon = 0;
   pf_mus_isGlobalMuon = 0;
   pf_mus_isElectron = 0;
   pf_mus_isConvertedPhoton = 0;
   pf_mus_isPhoton = 0;
   pf_mus_id_All = 0;
   pf_mus_id_AllGlobalMuons = 0;
   pf_mus_id_AllStandAloneMuons = 0;
   pf_mus_id_AllTrackerMuons = 0;
   pf_mus_id_TrackerMuonArbitrated = 0;
   pf_mus_id_AllArbitrated = 0;
   pf_mus_id_GlobalMuonPromptTight = 0;
   pf_mus_id_TMLastStationLoose = 0;
   pf_mus_id_TMLastStationTight = 0;
   pf_mus_id_TM2DCompatibilityLoose = 0;
   pf_mus_id_TM2DCompatibilityTight = 0;
   pf_mus_id_TMOneStationLoose = 0;
   pf_mus_id_TMOneStationTight = 0;
   pf_mus_id_TMLastStationOptimizedLowPtLoose = 0;
   pf_mus_id_TMLastStationOptimizedLowPtTight = 0;
   pf_mus_tk_LayersWithMeasurement = 0;
   pf_mus_tk_PixelLayersWithMeasurement = 0;
   pf_mus_tk_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_tk_LayersWithoutMeasurement = 0;
   pf_mus_tk_ExpectedHitsInner = 0;
   pf_mus_tk_ExpectedHitsOuter = 0;
   pf_mus_cm_LayersWithMeasurement = 0;
   pf_mus_cm_PixelLayersWithMeasurement = 0;
   pf_mus_cm_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_cm_LayersWithoutMeasurement = 0;
   pf_mus_cm_ExpectedHitsInner = 0;
   pf_mus_cm_ExpectedHitsOuter = 0;
   pf_mus_picky_LayersWithMeasurement = 0;
   pf_mus_picky_PixelLayersWithMeasurement = 0;
   pf_mus_picky_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_picky_LayersWithoutMeasurement = 0;
   pf_mus_picky_ExpectedHitsInner = 0;
   pf_mus_picky_ExpectedHitsOuter = 0;
   pf_mus_tpfms_LayersWithMeasurement = 0;
   pf_mus_tpfms_PixelLayersWithMeasurement = 0;
   pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit = 0;
   pf_mus_tpfms_LayersWithoutMeasurement = 0;
   pf_mus_tpfms_ExpectedHitsInner = 0;
   pf_mus_tpfms_ExpectedHitsOuter = 0;
   pf_mus_picky_id = 0;
   pf_mus_picky_chi2 = 0;
   pf_mus_picky_ndof = 0;
   pf_mus_picky_chg = 0;
   pf_mus_picky_pt = 0;
   pf_mus_picky_px = 0;
   pf_mus_picky_py = 0;
   pf_mus_picky_pz = 0;
   pf_mus_picky_eta = 0;
   pf_mus_picky_phi = 0;
   pf_mus_picky_theta = 0;
   pf_mus_picky_d0dum = 0;
   pf_mus_picky_dz = 0;
   pf_mus_picky_vx = 0;
   pf_mus_picky_vy = 0;
   pf_mus_picky_vz = 0;
   pf_mus_picky_numvalhits = 0;
   pf_mus_picky_numlosthits = 0;
   pf_mus_picky_d0dumErr = 0;
   pf_mus_picky_dzErr = 0;
   pf_mus_picky_ptErr = 0;
   pf_mus_picky_etaErr = 0;
   pf_mus_picky_phiErr = 0;
   pf_mus_picky_numvalPixelhits = 0;
   pf_mus_tpfms_id = 0;
   pf_mus_tpfms_chi2 = 0;
   pf_mus_tpfms_ndof = 0;
   pf_mus_tpfms_chg = 0;
   pf_mus_tpfms_pt = 0;
   pf_mus_tpfms_px = 0;
   pf_mus_tpfms_py = 0;
   pf_mus_tpfms_pz = 0;
   pf_mus_tpfms_eta = 0;
   pf_mus_tpfms_phi = 0;
   pf_mus_tpfms_theta = 0;
   pf_mus_tpfms_d0dum = 0;
   pf_mus_tpfms_dz = 0;
   pf_mus_tpfms_vx = 0;
   pf_mus_tpfms_vy = 0;
   pf_mus_tpfms_vz = 0;
   pf_mus_tpfms_numvalhits = 0;
   pf_mus_tpfms_numlosthits = 0;
   pf_mus_tpfms_d0dumErr = 0;
   pf_mus_tpfms_dzErr = 0;
   pf_mus_tpfms_ptErr = 0;
   pf_mus_tpfms_etaErr = 0;
   pf_mus_tpfms_phiErr = 0;
   pf_mus_tpfms_numvalPixelhits = 0;
   pf_mus_pfIsolationR03_sumChargedHadronPt = 0;
   pf_mus_pfIsolationR03_sumChargedParticlePt = 0;
   pf_mus_pfIsolationR03_sumNeutralHadronEt = 0;
   pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold = 0;
   pf_mus_pfIsolationR03_sumPhotonEt = 0;
   pf_mus_pfIsolationR03_sumPhotonEtHighThreshold = 0;
   pf_mus_pfIsolationR03_sumPUPt = 0;
   pf_mus_pfIsolationR04_sumChargedHadronPt = 0;
   pf_mus_pfIsolationR04_sumChargedParticlePt = 0;
   pf_mus_pfIsolationR04_sumNeutralHadronEt = 0;
   pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold = 0;
   pf_mus_pfIsolationR04_sumPhotonEt = 0;
   pf_mus_pfIsolationR04_sumPhotonEtHighThreshold = 0;
   pf_mus_pfIsolationR04_sumPUPt = 0;
   pf_mus_dB = 0;
   pf_mus_numberOfMatchedStations = 0;
   pf_mus_isPFMuon = 0;
   pf_photons_energy = 0;
   pf_photons_et = 0;
   pf_photons_eta = 0;
   pf_photons_phi = 0;
   pf_photons_pt = 0;
   pf_photons_px = 0;
   pf_photons_py = 0;
   pf_photons_pz = 0;
   pf_photons_status = 0;
   pf_photons_theta = 0;
   pf_photons_hadOverEM = 0;
   pf_photons_hadTowOverEM = 0;
   pf_photons_scEnergy = 0;
   pf_photons_scRawEnergy = 0;
   pf_photons_scEta = 0;
   pf_photons_scPhi = 0;
   pf_photons_scEtaWidth = 0;
   pf_photons_scPhiWidth = 0;
   pf_photons_isAlsoElectron = 0;
   pf_photons_hasPixelSeed = 0;
   pf_photons_isConverted = 0;
   pf_photons_isEBGap = 0;
   pf_photons_isEEGap = 0;
   pf_photons_isEBEEGap = 0;
   pf_photons_isEBPho = 0;
   pf_photons_isEEPho = 0;
   pf_photons_maxEnergyXtal = 0;
   pf_photons_e1x5 = 0;
   pf_photons_e2x5 = 0;
   pf_photons_e3x3 = 0;
   pf_photons_e5x5 = 0;
   pf_photons_sigmaEtaEta = 0;
   pf_photons_sigmaIetaIeta = 0;
   pf_photons_r9 = 0;
   pf_photons_chIso = 0;
   pf_photons_nhIso = 0;
   pf_photons_phIso = 0;
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
   pfmets_sign = 0;
   pfmets_sumEt = 0;
   pfmets_unCPhi = 0;
   pfmets_unCPt = 0;
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
   tcmets_et = 0;
   tcmets_phi = 0;
   tcmets_ex = 0;
   tcmets_ey = 0;
   tcmets_sumEt = 0;
   tracks_chi2 = 0;
   tracks_ndof = 0;
   tracks_chg = 0;
   tracks_pt = 0;
   tracks_px = 0;
   tracks_py = 0;
   tracks_pz = 0;
   tracks_eta = 0;
   tracks_phi = 0;
   tracks_d0dum = 0;
   tracks_dz = 0;
   tracks_vx = 0;
   tracks_vy = 0;
   tracks_vz = 0;
   tracks_numvalhits = 0;
   tracks_numlosthits = 0;
   tracks_d0dumErr = 0;
   tracks_dzErr = 0;
   tracks_ptErr = 0;
   tracks_etaErr = 0;
   tracks_phiErr = 0;
   tracks_highPurity = 0;
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
   fastjets_AK5PF_R1p2_R0p5pT10_px = 0;
   fastjets_AK5PF_R1p2_R0p5pT10_py = 0;
   fastjets_AK5PF_R1p2_R0p5pT10_pz = 0;
   fastjets_AK5PF_R1p2_R0p5pT10_energy = 0;
   fastjets_AK5PF_R1p2_R0p5pT10_phi = 0;
   fastjets_AK5PF_R1p2_R0p5pT10_eta = 0;
   fastjets_AK5PF_R1p2_R0p5pT10_index = 0;
   fastjets_AK5PF_R1p2_R0p5pT10_nconstituents = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_px = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_py = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_pz = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_energy = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_phi = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_eta = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_index = 0;
   fastjets_AK5PF_R1p2_R0p5pT15_nconstituents = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_px = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_py = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_pz = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_energy = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_phi = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_eta = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_index = 0;
   fastjets_AK5PF_R1p2_R0p5pT20_nconstituents = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_px = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_py = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_pz = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_energy = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_phi = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_eta = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_index = 0;
   fastjets_AK5PF_R1p2_R0p5pT25_nconstituents = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_px = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_py = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_pz = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_energy = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_phi = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_eta = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_index = 0;
   fastjets_AK5PF_R1p2_R0p5pT30_nconstituents = 0;

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
   fChain->SetBranchAddress("Njets_AK5PF", &Njets_AK5PF, &b_Njets_AK5PF);
   fChain->SetBranchAddress("jets_AK5PF_status", &jets_AK5PF_status, &b_jets_AK5PF_status);
   fChain->SetBranchAddress("jets_AK5PF_phi", &jets_AK5PF_phi, &b_jets_AK5PF_phi);
   fChain->SetBranchAddress("jets_AK5PF_pt", &jets_AK5PF_pt, &b_jets_AK5PF_pt);
   fChain->SetBranchAddress("jets_AK5PF_pz", &jets_AK5PF_pz, &b_jets_AK5PF_pz);
   fChain->SetBranchAddress("jets_AK5PF_px", &jets_AK5PF_px, &b_jets_AK5PF_px);
   fChain->SetBranchAddress("jets_AK5PF_py", &jets_AK5PF_py, &b_jets_AK5PF_py);
   fChain->SetBranchAddress("jets_AK5PF_eta", &jets_AK5PF_eta, &b_jets_AK5PF_eta);
   fChain->SetBranchAddress("jets_AK5PF_theta", &jets_AK5PF_theta, &b_jets_AK5PF_theta);
   fChain->SetBranchAddress("jets_AK5PF_et", &jets_AK5PF_et, &b_jets_AK5PF_et);
   fChain->SetBranchAddress("jets_AK5PF_energy", &jets_AK5PF_energy, &b_jets_AK5PF_energy);
   fChain->SetBranchAddress("jets_AK5PF_parton_Id", &jets_AK5PF_parton_Id, &b_jets_AK5PF_parton_Id);
   fChain->SetBranchAddress("jets_AK5PF_parton_motherId", &jets_AK5PF_parton_motherId, &b_jets_AK5PF_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PF_parton_pt", &jets_AK5PF_parton_pt, &b_jets_AK5PF_parton_pt);
   fChain->SetBranchAddress("jets_AK5PF_parton_phi", &jets_AK5PF_parton_phi, &b_jets_AK5PF_parton_phi);
   fChain->SetBranchAddress("jets_AK5PF_parton_eta", &jets_AK5PF_parton_eta, &b_jets_AK5PF_parton_eta);
   fChain->SetBranchAddress("jets_AK5PF_parton_Energy", &jets_AK5PF_parton_Energy, &b_jets_AK5PF_parton_Energy);
   fChain->SetBranchAddress("jets_AK5PF_parton_mass", &jets_AK5PF_parton_mass, &b_jets_AK5PF_parton_mass);
   fChain->SetBranchAddress("jets_AK5PF_gen_et", &jets_AK5PF_gen_et, &b_jets_AK5PF_gen_et);
   fChain->SetBranchAddress("jets_AK5PF_gen_pt", &jets_AK5PF_gen_pt, &b_jets_AK5PF_gen_pt);
   fChain->SetBranchAddress("jets_AK5PF_gen_eta", &jets_AK5PF_gen_eta, &b_jets_AK5PF_gen_eta);
   fChain->SetBranchAddress("jets_AK5PF_gen_phi", &jets_AK5PF_gen_phi, &b_jets_AK5PF_gen_phi);
   fChain->SetBranchAddress("jets_AK5PF_gen_mass", &jets_AK5PF_gen_mass, &b_jets_AK5PF_gen_mass);
   fChain->SetBranchAddress("jets_AK5PF_gen_Energy", &jets_AK5PF_gen_Energy, &b_jets_AK5PF_gen_Energy);
   fChain->SetBranchAddress("jets_AK5PF_gen_Id", &jets_AK5PF_gen_Id, &b_jets_AK5PF_gen_Id);
   fChain->SetBranchAddress("jets_AK5PF_gen_motherID", &jets_AK5PF_gen_motherID, &b_jets_AK5PF_gen_motherID);
   fChain->SetBranchAddress("jets_AK5PF_gen_threeCharge", &jets_AK5PF_gen_threeCharge, &b_jets_AK5PF_gen_threeCharge);
   fChain->SetBranchAddress("jets_AK5PF_partonFlavour", &jets_AK5PF_partonFlavour, &b_jets_AK5PF_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highPur", &jets_AK5PF_btag_TC_highPur, &b_jets_AK5PF_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highEff", &jets_AK5PF_btag_TC_highEff, &b_jets_AK5PF_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetProb", &jets_AK5PF_btag_jetProb, &b_jets_AK5PF_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetBProb", &jets_AK5PF_btag_jetBProb, &b_jets_AK5PF_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_softEle", &jets_AK5PF_btag_softEle, &b_jets_AK5PF_btag_softEle);
   fChain->SetBranchAddress("jets_AK5PF_btag_softMuon", &jets_AK5PF_btag_softMuon, &b_jets_AK5PF_btag_softMuon);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighPur", &jets_AK5PF_btag_secVertexHighPur, &b_jets_AK5PF_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighEff", &jets_AK5PF_btag_secVertexHighEff, &b_jets_AK5PF_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexCombined", &jets_AK5PF_btag_secVertexCombined, &b_jets_AK5PF_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PF_jetCharge", &jets_AK5PF_jetCharge, &b_jets_AK5PF_jetCharge);
   fChain->SetBranchAddress("jets_AK5PF_chgEmE", &jets_AK5PF_chgEmE, &b_jets_AK5PF_chgEmE);
   fChain->SetBranchAddress("jets_AK5PF_chgHadE", &jets_AK5PF_chgHadE, &b_jets_AK5PF_chgHadE);
   fChain->SetBranchAddress("jets_AK5PF_photonEnergy", &jets_AK5PF_photonEnergy, &b_jets_AK5PF_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PF_chgMuE", &jets_AK5PF_chgMuE, &b_jets_AK5PF_chgMuE);
   fChain->SetBranchAddress("jets_AK5PF_chg_Mult", &jets_AK5PF_chg_Mult, &b_jets_AK5PF_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PF_neutralEmE", &jets_AK5PF_neutralEmE, &b_jets_AK5PF_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PF_neutralHadE", &jets_AK5PF_neutralHadE, &b_jets_AK5PF_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PF_neutral_Mult", &jets_AK5PF_neutral_Mult, &b_jets_AK5PF_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PF_mu_Mult", &jets_AK5PF_mu_Mult, &b_jets_AK5PF_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PF_emf", &jets_AK5PF_emf, &b_jets_AK5PF_emf);
   fChain->SetBranchAddress("jets_AK5PF_ehf", &jets_AK5PF_ehf, &b_jets_AK5PF_ehf);
   fChain->SetBranchAddress("jets_AK5PF_n60", &jets_AK5PF_n60, &b_jets_AK5PF_n60);
   fChain->SetBranchAddress("jets_AK5PF_n90", &jets_AK5PF_n90, &b_jets_AK5PF_n90);
   fChain->SetBranchAddress("jets_AK5PF_etaetaMoment", &jets_AK5PF_etaetaMoment, &b_jets_AK5PF_etaetaMoment);
   fChain->SetBranchAddress("jets_AK5PF_etaphiMoment", &jets_AK5PF_etaphiMoment, &b_jets_AK5PF_etaphiMoment);
   fChain->SetBranchAddress("jets_AK5PF_phiphiMoment", &jets_AK5PF_phiphiMoment, &b_jets_AK5PF_phiphiMoment);
   fChain->SetBranchAddress("jets_AK5PF_n90Hits", &jets_AK5PF_n90Hits, &b_jets_AK5PF_n90Hits);
   fChain->SetBranchAddress("jets_AK5PF_fHPD", &jets_AK5PF_fHPD, &b_jets_AK5PF_fHPD);
   fChain->SetBranchAddress("jets_AK5PF_fRBX", &jets_AK5PF_fRBX, &b_jets_AK5PF_fRBX);
   fChain->SetBranchAddress("jets_AK5PF_hitsInN90", &jets_AK5PF_hitsInN90, &b_jets_AK5PF_hitsInN90);
   fChain->SetBranchAddress("jets_AK5PF_nECALTowers", &jets_AK5PF_nECALTowers, &b_jets_AK5PF_nECALTowers);
   fChain->SetBranchAddress("jets_AK5PF_nHCALTowers", &jets_AK5PF_nHCALTowers, &b_jets_AK5PF_nHCALTowers);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector1", &jets_AK5PF_fSubDetector1, &b_jets_AK5PF_fSubDetector1);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector2", &jets_AK5PF_fSubDetector2, &b_jets_AK5PF_fSubDetector2);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector3", &jets_AK5PF_fSubDetector3, &b_jets_AK5PF_fSubDetector3);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector4", &jets_AK5PF_fSubDetector4, &b_jets_AK5PF_fSubDetector4);
   fChain->SetBranchAddress("jets_AK5PF_area", &jets_AK5PF_area, &b_jets_AK5PF_area);
   fChain->SetBranchAddress("jets_AK5PF_corrFactorRaw", &jets_AK5PF_corrFactorRaw, &b_jets_AK5PF_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK5PF_rawPt", &jets_AK5PF_rawPt, &b_jets_AK5PF_rawPt);
   fChain->SetBranchAddress("jets_AK5PF_mass", &jets_AK5PF_mass, &b_jets_AK5PF_mass);
   fChain->SetBranchAddress("Njets_AK5PFclean", &Njets_AK5PFclean, &b_Njets_AK5PFclean);
   fChain->SetBranchAddress("jets_AK5PFclean_status", &jets_AK5PFclean_status, &b_jets_AK5PFclean_status);
   fChain->SetBranchAddress("jets_AK5PFclean_phi", &jets_AK5PFclean_phi, &b_jets_AK5PFclean_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_pt", &jets_AK5PFclean_pt, &b_jets_AK5PFclean_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_pz", &jets_AK5PFclean_pz, &b_jets_AK5PFclean_pz);
   fChain->SetBranchAddress("jets_AK5PFclean_px", &jets_AK5PFclean_px, &b_jets_AK5PFclean_px);
   fChain->SetBranchAddress("jets_AK5PFclean_py", &jets_AK5PFclean_py, &b_jets_AK5PFclean_py);
   fChain->SetBranchAddress("jets_AK5PFclean_eta", &jets_AK5PFclean_eta, &b_jets_AK5PFclean_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_theta", &jets_AK5PFclean_theta, &b_jets_AK5PFclean_theta);
   fChain->SetBranchAddress("jets_AK5PFclean_et", &jets_AK5PFclean_et, &b_jets_AK5PFclean_et);
   fChain->SetBranchAddress("jets_AK5PFclean_energy", &jets_AK5PFclean_energy, &b_jets_AK5PFclean_energy);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_Id", &jets_AK5PFclean_parton_Id, &b_jets_AK5PFclean_parton_Id);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_motherId", &jets_AK5PFclean_parton_motherId, &b_jets_AK5PFclean_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_pt", &jets_AK5PFclean_parton_pt, &b_jets_AK5PFclean_parton_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_phi", &jets_AK5PFclean_parton_phi, &b_jets_AK5PFclean_parton_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_eta", &jets_AK5PFclean_parton_eta, &b_jets_AK5PFclean_parton_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_Energy", &jets_AK5PFclean_parton_Energy, &b_jets_AK5PFclean_parton_Energy);
   fChain->SetBranchAddress("jets_AK5PFclean_parton_mass", &jets_AK5PFclean_parton_mass, &b_jets_AK5PFclean_parton_mass);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_et", &jets_AK5PFclean_gen_et, &b_jets_AK5PFclean_gen_et);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_pt", &jets_AK5PFclean_gen_pt, &b_jets_AK5PFclean_gen_pt);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_eta", &jets_AK5PFclean_gen_eta, &b_jets_AK5PFclean_gen_eta);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_phi", &jets_AK5PFclean_gen_phi, &b_jets_AK5PFclean_gen_phi);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_mass", &jets_AK5PFclean_gen_mass, &b_jets_AK5PFclean_gen_mass);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_Energy", &jets_AK5PFclean_gen_Energy, &b_jets_AK5PFclean_gen_Energy);
   fChain->SetBranchAddress("jets_AK5PFclean_gen_Id", &jets_AK5PFclean_gen_Id, &b_jets_AK5PFclean_gen_Id);
   fChain->SetBranchAddress("jets_AK5PFclean_partonFlavour", &jets_AK5PFclean_partonFlavour, &b_jets_AK5PFclean_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highPur", &jets_AK5PFclean_btag_TC_highPur, &b_jets_AK5PFclean_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_TC_highEff", &jets_AK5PFclean_btag_TC_highEff, &b_jets_AK5PFclean_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetProb", &jets_AK5PFclean_btag_jetProb, &b_jets_AK5PFclean_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_jetBProb", &jets_AK5PFclean_btag_jetBProb, &b_jets_AK5PFclean_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_softEle", &jets_AK5PFclean_btag_softEle, &b_jets_AK5PFclean_btag_softEle);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_softMuon", &jets_AK5PFclean_btag_softMuon, &b_jets_AK5PFclean_btag_softMuon);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighPur", &jets_AK5PFclean_btag_secVertexHighPur, &b_jets_AK5PFclean_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexHighEff", &jets_AK5PFclean_btag_secVertexHighEff, &b_jets_AK5PFclean_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PFclean_btag_secVertexCombined", &jets_AK5PFclean_btag_secVertexCombined, &b_jets_AK5PFclean_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PFclean_jetCharge", &jets_AK5PFclean_jetCharge, &b_jets_AK5PFclean_jetCharge);
   fChain->SetBranchAddress("jets_AK5PFclean_chgEmE", &jets_AK5PFclean_chgEmE, &b_jets_AK5PFclean_chgEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_chgHadE", &jets_AK5PFclean_chgHadE, &b_jets_AK5PFclean_chgHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_photonEnergy", &jets_AK5PFclean_photonEnergy, &b_jets_AK5PFclean_photonEnergy);
   fChain->SetBranchAddress("jets_AK5PFclean_chgMuE", &jets_AK5PFclean_chgMuE, &b_jets_AK5PFclean_chgMuE);
   fChain->SetBranchAddress("jets_AK5PFclean_chg_Mult", &jets_AK5PFclean_chg_Mult, &b_jets_AK5PFclean_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralEmE", &jets_AK5PFclean_neutralEmE, &b_jets_AK5PFclean_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutralHadE", &jets_AK5PFclean_neutralHadE, &b_jets_AK5PFclean_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PFclean_neutral_Mult", &jets_AK5PFclean_neutral_Mult, &b_jets_AK5PFclean_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_mu_Mult", &jets_AK5PFclean_mu_Mult, &b_jets_AK5PFclean_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PFclean_emf", &jets_AK5PFclean_emf, &b_jets_AK5PFclean_emf);
   fChain->SetBranchAddress("jets_AK5PFclean_ehf", &jets_AK5PFclean_ehf, &b_jets_AK5PFclean_ehf);
   fChain->SetBranchAddress("jets_AK5PFclean_n60", &jets_AK5PFclean_n60, &b_jets_AK5PFclean_n60);
   fChain->SetBranchAddress("jets_AK5PFclean_n90", &jets_AK5PFclean_n90, &b_jets_AK5PFclean_n90);
   fChain->SetBranchAddress("jets_AK5PFclean_etaetaMoment", &jets_AK5PFclean_etaetaMoment, &b_jets_AK5PFclean_etaetaMoment);
   fChain->SetBranchAddress("jets_AK5PFclean_etaphiMoment", &jets_AK5PFclean_etaphiMoment, &b_jets_AK5PFclean_etaphiMoment);
   fChain->SetBranchAddress("jets_AK5PFclean_phiphiMoment", &jets_AK5PFclean_phiphiMoment, &b_jets_AK5PFclean_phiphiMoment);
   fChain->SetBranchAddress("jets_AK5PFclean_n90Hits", &jets_AK5PFclean_n90Hits, &b_jets_AK5PFclean_n90Hits);
   fChain->SetBranchAddress("jets_AK5PFclean_fHPD", &jets_AK5PFclean_fHPD, &b_jets_AK5PFclean_fHPD);
   fChain->SetBranchAddress("jets_AK5PFclean_fRBX", &jets_AK5PFclean_fRBX, &b_jets_AK5PFclean_fRBX);
   fChain->SetBranchAddress("jets_AK5PFclean_hitsInN90", &jets_AK5PFclean_hitsInN90, &b_jets_AK5PFclean_hitsInN90);
   fChain->SetBranchAddress("jets_AK5PFclean_nECALTowers", &jets_AK5PFclean_nECALTowers, &b_jets_AK5PFclean_nECALTowers);
   fChain->SetBranchAddress("jets_AK5PFclean_nHCALTowers", &jets_AK5PFclean_nHCALTowers, &b_jets_AK5PFclean_nHCALTowers);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector1", &jets_AK5PFclean_fSubDetector1, &b_jets_AK5PFclean_fSubDetector1);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector2", &jets_AK5PFclean_fSubDetector2, &b_jets_AK5PFclean_fSubDetector2);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector3", &jets_AK5PFclean_fSubDetector3, &b_jets_AK5PFclean_fSubDetector3);
   fChain->SetBranchAddress("jets_AK5PFclean_fSubDetector4", &jets_AK5PFclean_fSubDetector4, &b_jets_AK5PFclean_fSubDetector4);
   fChain->SetBranchAddress("jets_AK5PFclean_area", &jets_AK5PFclean_area, &b_jets_AK5PFclean_area);
   fChain->SetBranchAddress("jets_AK5PFclean_corrFactorRaw", &jets_AK5PFclean_corrFactorRaw, &b_jets_AK5PFclean_corrFactorRaw);
   fChain->SetBranchAddress("jets_AK5PFclean_rawPt", &jets_AK5PFclean_rawPt, &b_jets_AK5PFclean_rawPt);
   fChain->SetBranchAddress("jets_AK5PFclean_mass", &jets_AK5PFclean_mass, &b_jets_AK5PFclean_mass);
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
   fChain->SetBranchAddress("Nmc_pdf", &Nmc_pdf, &b_Nmc_pdf);
   fChain->SetBranchAddress("mc_pdf_x1", &mc_pdf_x1, &b_mc_pdf_x1);
   fChain->SetBranchAddress("mc_pdf_x2", &mc_pdf_x2, &b_mc_pdf_x2);
   fChain->SetBranchAddress("mc_pdf_q", &mc_pdf_q, &b_mc_pdf_q);
   fChain->SetBranchAddress("mc_pdf_id1", &mc_pdf_id1, &b_mc_pdf_id1);
   fChain->SetBranchAddress("mc_pdf_id2", &mc_pdf_id2, &b_mc_pdf_id2);
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
   fChain->SetBranchAddress("NmetsHO", &NmetsHO, &b_NmetsHO);
   fChain->SetBranchAddress("metsHO_et", &metsHO_et, &b_metsHO_et);
   fChain->SetBranchAddress("metsHO_phi", &metsHO_phi, &b_metsHO_phi);
   fChain->SetBranchAddress("metsHO_ex", &metsHO_ex, &b_metsHO_ex);
   fChain->SetBranchAddress("metsHO_ey", &metsHO_ey, &b_metsHO_ey);
   fChain->SetBranchAddress("metsHO_sumEt", &metsHO_sumEt, &b_metsHO_sumEt);
   fChain->SetBranchAddress("Nmets_AK5", &Nmets_AK5, &b_Nmets_AK5);
   fChain->SetBranchAddress("mets_AK5_et", &mets_AK5_et, &b_mets_AK5_et);
   fChain->SetBranchAddress("mets_AK5_phi", &mets_AK5_phi, &b_mets_AK5_phi);
   fChain->SetBranchAddress("mets_AK5_ex", &mets_AK5_ex, &b_mets_AK5_ex);
   fChain->SetBranchAddress("mets_AK5_ey", &mets_AK5_ey, &b_mets_AK5_ey);
   fChain->SetBranchAddress("mets_AK5_gen_et", &mets_AK5_gen_et, &b_mets_AK5_gen_et);
   fChain->SetBranchAddress("mets_AK5_gen_phi", &mets_AK5_gen_phi, &b_mets_AK5_gen_phi);
   fChain->SetBranchAddress("mets_AK5_sign", &mets_AK5_sign, &b_mets_AK5_sign);
   fChain->SetBranchAddress("mets_AK5_sumEt", &mets_AK5_sumEt, &b_mets_AK5_sumEt);
   fChain->SetBranchAddress("mets_AK5_unCPhi", &mets_AK5_unCPhi, &b_mets_AK5_unCPhi);
   fChain->SetBranchAddress("mets_AK5_unCPt", &mets_AK5_unCPt, &b_mets_AK5_unCPt);
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
   fChain->SetBranchAddress("NpfTypeINoXYCorrmets", &NpfTypeINoXYCorrmets, &b_NpfTypeINoXYCorrmets);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_et", &pfTypeINoXYCorrmets_et, &b_pfTypeINoXYCorrmets_et);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_phi", &pfTypeINoXYCorrmets_phi, &b_pfTypeINoXYCorrmets_phi);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_ex", &pfTypeINoXYCorrmets_ex, &b_pfTypeINoXYCorrmets_ex);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_ey", &pfTypeINoXYCorrmets_ey, &b_pfTypeINoXYCorrmets_ey);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_gen_et", &pfTypeINoXYCorrmets_gen_et, &b_pfTypeINoXYCorrmets_gen_et);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_gen_phi", &pfTypeINoXYCorrmets_gen_phi, &b_pfTypeINoXYCorrmets_gen_phi);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_sign", &pfTypeINoXYCorrmets_sign, &b_pfTypeINoXYCorrmets_sign);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_sumEt", &pfTypeINoXYCorrmets_sumEt, &b_pfTypeINoXYCorrmets_sumEt);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_unCPhi", &pfTypeINoXYCorrmets_unCPhi, &b_pfTypeINoXYCorrmets_unCPhi);
   fChain->SetBranchAddress("pfTypeINoXYCorrmets_unCPt", &pfTypeINoXYCorrmets_unCPt, &b_pfTypeINoXYCorrmets_unCPt);
   fChain->SetBranchAddress("NpfTypeIType0mets", &NpfTypeIType0mets, &b_NpfTypeIType0mets);
   fChain->SetBranchAddress("pfTypeIType0mets_et", &pfTypeIType0mets_et, &b_pfTypeIType0mets_et);
   fChain->SetBranchAddress("pfTypeIType0mets_phi", &pfTypeIType0mets_phi, &b_pfTypeIType0mets_phi);
   fChain->SetBranchAddress("pfTypeIType0mets_ex", &pfTypeIType0mets_ex, &b_pfTypeIType0mets_ex);
   fChain->SetBranchAddress("pfTypeIType0mets_ey", &pfTypeIType0mets_ey, &b_pfTypeIType0mets_ey);
   fChain->SetBranchAddress("pfTypeIType0mets_gen_et", &pfTypeIType0mets_gen_et, &b_pfTypeIType0mets_gen_et);
   fChain->SetBranchAddress("pfTypeIType0mets_gen_phi", &pfTypeIType0mets_gen_phi, &b_pfTypeIType0mets_gen_phi);
   fChain->SetBranchAddress("pfTypeIType0mets_sign", &pfTypeIType0mets_sign, &b_pfTypeIType0mets_sign);
   fChain->SetBranchAddress("pfTypeIType0mets_sumEt", &pfTypeIType0mets_sumEt, &b_pfTypeIType0mets_sumEt);
   fChain->SetBranchAddress("pfTypeIType0mets_unCPhi", &pfTypeIType0mets_unCPhi, &b_pfTypeIType0mets_unCPhi);
   fChain->SetBranchAddress("pfTypeIType0mets_unCPt", &pfTypeIType0mets_unCPt, &b_pfTypeIType0mets_unCPt);
   fChain->SetBranchAddress("NpfTypeImets", &NpfTypeImets, &b_NpfTypeImets);
   fChain->SetBranchAddress("pfTypeImets_et", &pfTypeImets_et, &b_pfTypeImets_et);
   fChain->SetBranchAddress("pfTypeImets_phi", &pfTypeImets_phi, &b_pfTypeImets_phi);
   fChain->SetBranchAddress("pfTypeImets_ex", &pfTypeImets_ex, &b_pfTypeImets_ex);
   fChain->SetBranchAddress("pfTypeImets_ey", &pfTypeImets_ey, &b_pfTypeImets_ey);
   fChain->SetBranchAddress("pfTypeImets_gen_et", &pfTypeImets_gen_et, &b_pfTypeImets_gen_et);
   fChain->SetBranchAddress("pfTypeImets_gen_phi", &pfTypeImets_gen_phi, &b_pfTypeImets_gen_phi);
   fChain->SetBranchAddress("pfTypeImets_sign", &pfTypeImets_sign, &b_pfTypeImets_sign);
   fChain->SetBranchAddress("pfTypeImets_sumEt", &pfTypeImets_sumEt, &b_pfTypeImets_sumEt);
   fChain->SetBranchAddress("pfTypeImets_unCPhi", &pfTypeImets_unCPhi, &b_pfTypeImets_unCPhi);
   fChain->SetBranchAddress("pfTypeImets_unCPt", &pfTypeImets_unCPt, &b_pfTypeImets_unCPt);
   fChain->SetBranchAddress("Npf_els", &Npf_els, &b_Npf_els);
   fChain->SetBranchAddress("pf_els_energy", &pf_els_energy, &b_pf_els_energy);
   fChain->SetBranchAddress("pf_els_et", &pf_els_et, &b_pf_els_et);
   fChain->SetBranchAddress("pf_els_eta", &pf_els_eta, &b_pf_els_eta);
   fChain->SetBranchAddress("pf_els_phi", &pf_els_phi, &b_pf_els_phi);
   fChain->SetBranchAddress("pf_els_pt", &pf_els_pt, &b_pf_els_pt);
   fChain->SetBranchAddress("pf_els_px", &pf_els_px, &b_pf_els_px);
   fChain->SetBranchAddress("pf_els_py", &pf_els_py, &b_pf_els_py);
   fChain->SetBranchAddress("pf_els_pz", &pf_els_pz, &b_pf_els_pz);
   fChain->SetBranchAddress("pf_els_status", &pf_els_status, &b_pf_els_status);
   fChain->SetBranchAddress("pf_els_theta", &pf_els_theta, &b_pf_els_theta);
   fChain->SetBranchAddress("pf_els_gen_id", &pf_els_gen_id, &b_pf_els_gen_id);
   fChain->SetBranchAddress("pf_els_gen_phi", &pf_els_gen_phi, &b_pf_els_gen_phi);
   fChain->SetBranchAddress("pf_els_gen_pt", &pf_els_gen_pt, &b_pf_els_gen_pt);
   fChain->SetBranchAddress("pf_els_gen_pz", &pf_els_gen_pz, &b_pf_els_gen_pz);
   fChain->SetBranchAddress("pf_els_gen_px", &pf_els_gen_px, &b_pf_els_gen_px);
   fChain->SetBranchAddress("pf_els_gen_py", &pf_els_gen_py, &b_pf_els_gen_py);
   fChain->SetBranchAddress("pf_els_gen_eta", &pf_els_gen_eta, &b_pf_els_gen_eta);
   fChain->SetBranchAddress("pf_els_gen_theta", &pf_els_gen_theta, &b_pf_els_gen_theta);
   fChain->SetBranchAddress("pf_els_gen_et", &pf_els_gen_et, &b_pf_els_gen_et);
   fChain->SetBranchAddress("pf_els_gen_mother_id", &pf_els_gen_mother_id, &b_pf_els_gen_mother_id);
   fChain->SetBranchAddress("pf_els_gen_mother_phi", &pf_els_gen_mother_phi, &b_pf_els_gen_mother_phi);
   fChain->SetBranchAddress("pf_els_gen_mother_pt", &pf_els_gen_mother_pt, &b_pf_els_gen_mother_pt);
   fChain->SetBranchAddress("pf_els_gen_mother_pz", &pf_els_gen_mother_pz, &b_pf_els_gen_mother_pz);
   fChain->SetBranchAddress("pf_els_gen_mother_px", &pf_els_gen_mother_px, &b_pf_els_gen_mother_px);
   fChain->SetBranchAddress("pf_els_gen_mother_py", &pf_els_gen_mother_py, &b_pf_els_gen_mother_py);
   fChain->SetBranchAddress("pf_els_gen_mother_eta", &pf_els_gen_mother_eta, &b_pf_els_gen_mother_eta);
   fChain->SetBranchAddress("pf_els_gen_mother_theta", &pf_els_gen_mother_theta, &b_pf_els_gen_mother_theta);
   fChain->SetBranchAddress("pf_els_gen_mother_et", &pf_els_gen_mother_et, &b_pf_els_gen_mother_et);
   fChain->SetBranchAddress("pf_els_tightId", &pf_els_tightId, &b_pf_els_tightId);
   fChain->SetBranchAddress("pf_els_looseId", &pf_els_looseId, &b_pf_els_looseId);
   fChain->SetBranchAddress("pf_els_robustTightId", &pf_els_robustTightId, &b_pf_els_robustTightId);
   fChain->SetBranchAddress("pf_els_robustLooseId", &pf_els_robustLooseId, &b_pf_els_robustLooseId);
   fChain->SetBranchAddress("pf_els_robustHighEnergyId", &pf_els_robustHighEnergyId, &b_pf_els_robustHighEnergyId);
   fChain->SetBranchAddress("pf_els_simpleEleId95relIso", &pf_els_simpleEleId95relIso, &b_pf_els_simpleEleId95relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90relIso", &pf_els_simpleEleId90relIso, &b_pf_els_simpleEleId90relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85relIso", &pf_els_simpleEleId85relIso, &b_pf_els_simpleEleId85relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80relIso", &pf_els_simpleEleId80relIso, &b_pf_els_simpleEleId80relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70relIso", &pf_els_simpleEleId70relIso, &b_pf_els_simpleEleId70relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId60relIso", &pf_els_simpleEleId60relIso, &b_pf_els_simpleEleId60relIso);
   fChain->SetBranchAddress("pf_els_simpleEleId95cIso", &pf_els_simpleEleId95cIso, &b_pf_els_simpleEleId95cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId90cIso", &pf_els_simpleEleId90cIso, &b_pf_els_simpleEleId90cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId85cIso", &pf_els_simpleEleId85cIso, &b_pf_els_simpleEleId85cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId80cIso", &pf_els_simpleEleId80cIso, &b_pf_els_simpleEleId80cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId70cIso", &pf_els_simpleEleId70cIso, &b_pf_els_simpleEleId70cIso);
   fChain->SetBranchAddress("pf_els_simpleEleId60cIso", &pf_els_simpleEleId60cIso, &b_pf_els_simpleEleId60cIso);
   fChain->SetBranchAddress("pf_els_cIso", &pf_els_cIso, &b_pf_els_cIso);
   fChain->SetBranchAddress("pf_els_tIso", &pf_els_tIso, &b_pf_els_tIso);
   fChain->SetBranchAddress("pf_els_ecalIso", &pf_els_ecalIso, &b_pf_els_ecalIso);
   fChain->SetBranchAddress("pf_els_hcalIso", &pf_els_hcalIso, &b_pf_els_hcalIso);
   fChain->SetBranchAddress("pf_els_chargedHadronIso", &pf_els_chargedHadronIso, &b_pf_els_chargedHadronIso);
   fChain->SetBranchAddress("pf_els_photonIso", &pf_els_photonIso, &b_pf_els_photonIso);
   fChain->SetBranchAddress("pf_els_neutralHadronIso", &pf_els_neutralHadronIso, &b_pf_els_neutralHadronIso);
   fChain->SetBranchAddress("pf_els_chi2", &pf_els_chi2, &b_pf_els_chi2);
   fChain->SetBranchAddress("pf_els_charge", &pf_els_charge, &b_pf_els_charge);
   fChain->SetBranchAddress("pf_els_caloEnergy", &pf_els_caloEnergy, &b_pf_els_caloEnergy);
   fChain->SetBranchAddress("pf_els_hadOverEm", &pf_els_hadOverEm, &b_pf_els_hadOverEm);
   fChain->SetBranchAddress("pf_els_hcalOverEcalBc", &pf_els_hcalOverEcalBc, &b_pf_els_hcalOverEcalBc);
   fChain->SetBranchAddress("pf_els_eOverPIn", &pf_els_eOverPIn, &b_pf_els_eOverPIn);
   fChain->SetBranchAddress("pf_els_eSeedOverPOut", &pf_els_eSeedOverPOut, &b_pf_els_eSeedOverPOut);
   fChain->SetBranchAddress("pf_els_sigmaEtaEta", &pf_els_sigmaEtaEta, &b_pf_els_sigmaEtaEta);
   fChain->SetBranchAddress("pf_els_sigmaIEtaIEta", &pf_els_sigmaIEtaIEta, &b_pf_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("pf_els_scEnergy", &pf_els_scEnergy, &b_pf_els_scEnergy);
   fChain->SetBranchAddress("pf_els_scRawEnergy", &pf_els_scRawEnergy, &b_pf_els_scRawEnergy);
   fChain->SetBranchAddress("pf_els_scSeedEnergy", &pf_els_scSeedEnergy, &b_pf_els_scSeedEnergy);
   fChain->SetBranchAddress("pf_els_scEta", &pf_els_scEta, &b_pf_els_scEta);
   fChain->SetBranchAddress("pf_els_scPhi", &pf_els_scPhi, &b_pf_els_scPhi);
   fChain->SetBranchAddress("pf_els_scEtaWidth", &pf_els_scEtaWidth, &b_pf_els_scEtaWidth);
   fChain->SetBranchAddress("pf_els_scPhiWidth", &pf_els_scPhiWidth, &b_pf_els_scPhiWidth);
   fChain->SetBranchAddress("pf_els_scE1x5", &pf_els_scE1x5, &b_pf_els_scE1x5);
   fChain->SetBranchAddress("pf_els_scE2x5Max", &pf_els_scE2x5Max, &b_pf_els_scE2x5Max);
   fChain->SetBranchAddress("pf_els_scE5x5", &pf_els_scE5x5, &b_pf_els_scE5x5);
   fChain->SetBranchAddress("pf_els_isEB", &pf_els_isEB, &b_pf_els_isEB);
   fChain->SetBranchAddress("pf_els_isEE", &pf_els_isEE, &b_pf_els_isEE);
   fChain->SetBranchAddress("pf_els_dEtaIn", &pf_els_dEtaIn, &b_pf_els_dEtaIn);
   fChain->SetBranchAddress("pf_els_dPhiIn", &pf_els_dPhiIn, &b_pf_els_dPhiIn);
   fChain->SetBranchAddress("pf_els_dEtaOut", &pf_els_dEtaOut, &b_pf_els_dEtaOut);
   fChain->SetBranchAddress("pf_els_dPhiOut", &pf_els_dPhiOut, &b_pf_els_dPhiOut);
   fChain->SetBranchAddress("pf_els_numvalhits", &pf_els_numvalhits, &b_pf_els_numvalhits);
   fChain->SetBranchAddress("pf_els_numlosthits", &pf_els_numlosthits, &b_pf_els_numlosthits);
   fChain->SetBranchAddress("pf_els_basicClustersSize", &pf_els_basicClustersSize, &b_pf_els_basicClustersSize);
   fChain->SetBranchAddress("pf_els_tk_pz", &pf_els_tk_pz, &b_pf_els_tk_pz);
   fChain->SetBranchAddress("pf_els_tk_pt", &pf_els_tk_pt, &b_pf_els_tk_pt);
   fChain->SetBranchAddress("pf_els_tk_phi", &pf_els_tk_phi, &b_pf_els_tk_phi);
   fChain->SetBranchAddress("pf_els_tk_eta", &pf_els_tk_eta, &b_pf_els_tk_eta);
   fChain->SetBranchAddress("pf_els_d0dum", &pf_els_d0dum, &b_pf_els_d0dum);
   fChain->SetBranchAddress("pf_els_dz", &pf_els_dz, &b_pf_els_dz);
   fChain->SetBranchAddress("pf_els_vx", &pf_els_vx, &b_pf_els_vx);
   fChain->SetBranchAddress("pf_els_vy", &pf_els_vy, &b_pf_els_vy);
   fChain->SetBranchAddress("pf_els_vz", &pf_els_vz, &b_pf_els_vz);
   fChain->SetBranchAddress("pf_els_ndof", &pf_els_ndof, &b_pf_els_ndof);
   fChain->SetBranchAddress("pf_els_ptError", &pf_els_ptError, &b_pf_els_ptError);
   fChain->SetBranchAddress("pf_els_d0dumError", &pf_els_d0dumError, &b_pf_els_d0dumError);
   fChain->SetBranchAddress("pf_els_dzError", &pf_els_dzError, &b_pf_els_dzError);
   fChain->SetBranchAddress("pf_els_etaError", &pf_els_etaError, &b_pf_els_etaError);
   fChain->SetBranchAddress("pf_els_phiError", &pf_els_phiError, &b_pf_els_phiError);
   fChain->SetBranchAddress("pf_els_tk_charge", &pf_els_tk_charge, &b_pf_els_tk_charge);
   fChain->SetBranchAddress("pf_els_core_ecalDrivenSeed", &pf_els_core_ecalDrivenSeed, &b_pf_els_core_ecalDrivenSeed);
   fChain->SetBranchAddress("pf_els_n_inner_layer", &pf_els_n_inner_layer, &b_pf_els_n_inner_layer);
   fChain->SetBranchAddress("pf_els_n_outer_layer", &pf_els_n_outer_layer, &b_pf_els_n_outer_layer);
   fChain->SetBranchAddress("pf_els_ctf_tk_id", &pf_els_ctf_tk_id, &b_pf_els_ctf_tk_id);
   fChain->SetBranchAddress("pf_els_ctf_tk_charge", &pf_els_ctf_tk_charge, &b_pf_els_ctf_tk_charge);
   fChain->SetBranchAddress("pf_els_ctf_tk_eta", &pf_els_ctf_tk_eta, &b_pf_els_ctf_tk_eta);
   fChain->SetBranchAddress("pf_els_ctf_tk_phi", &pf_els_ctf_tk_phi, &b_pf_els_ctf_tk_phi);
   fChain->SetBranchAddress("pf_els_fbrem", &pf_els_fbrem, &b_pf_els_fbrem);
   fChain->SetBranchAddress("pf_els_shFracInnerHits", &pf_els_shFracInnerHits, &b_pf_els_shFracInnerHits);
   fChain->SetBranchAddress("pf_els_dr03EcalRecHitSumEt", &pf_els_dr03EcalRecHitSumEt, &b_pf_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalTowerSumEt", &pf_els_dr03HcalTowerSumEt, &b_pf_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth1TowerSumEt", &pf_els_dr03HcalDepth1TowerSumEt, &b_pf_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03HcalDepth2TowerSumEt", &pf_els_dr03HcalDepth2TowerSumEt, &b_pf_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr03TkSumPt", &pf_els_dr03TkSumPt, &b_pf_els_dr03TkSumPt);
   fChain->SetBranchAddress("pf_els_dr04EcalRecHitSumEt", &pf_els_dr04EcalRecHitSumEt, &b_pf_els_dr04EcalRecHitSumEt);
   fChain->SetBranchAddress("pf_els_dr04HcalTowerSumEt", &pf_els_dr04HcalTowerSumEt, &b_pf_els_dr04HcalTowerSumEt);
   fChain->SetBranchAddress("pf_els_dr04HcalDepth1TowerSumEt", &pf_els_dr04HcalDepth1TowerSumEt, &b_pf_els_dr04HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr04HcalDepth2TowerSumEt", &pf_els_dr04HcalDepth2TowerSumEt, &b_pf_els_dr04HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("pf_els_dr04TkSumPt", &pf_els_dr04TkSumPt, &b_pf_els_dr04TkSumPt);
   fChain->SetBranchAddress("pf_els_cpx", &pf_els_cpx, &b_pf_els_cpx);
   fChain->SetBranchAddress("pf_els_cpy", &pf_els_cpy, &b_pf_els_cpy);
   fChain->SetBranchAddress("pf_els_cpz", &pf_els_cpz, &b_pf_els_cpz);
   fChain->SetBranchAddress("pf_els_vpx", &pf_els_vpx, &b_pf_els_vpx);
   fChain->SetBranchAddress("pf_els_vpy", &pf_els_vpy, &b_pf_els_vpy);
   fChain->SetBranchAddress("pf_els_vpz", &pf_els_vpz, &b_pf_els_vpz);
   fChain->SetBranchAddress("pf_els_cx", &pf_els_cx, &b_pf_els_cx);
   fChain->SetBranchAddress("pf_els_cy", &pf_els_cy, &b_pf_els_cy);
   fChain->SetBranchAddress("pf_els_cz", &pf_els_cz, &b_pf_els_cz);
   fChain->SetBranchAddress("pf_els_PATpassConversionVeto", &pf_els_PATpassConversionVeto, &b_pf_els_PATpassConversionVeto);
   fChain->SetBranchAddress("Npf_mus", &Npf_mus, &b_Npf_mus);
   fChain->SetBranchAddress("pf_mus_energy", &pf_mus_energy, &b_pf_mus_energy);
   fChain->SetBranchAddress("pf_mus_et", &pf_mus_et, &b_pf_mus_et);
   fChain->SetBranchAddress("pf_mus_eta", &pf_mus_eta, &b_pf_mus_eta);
   fChain->SetBranchAddress("pf_mus_phi", &pf_mus_phi, &b_pf_mus_phi);
   fChain->SetBranchAddress("pf_mus_pt", &pf_mus_pt, &b_pf_mus_pt);
   fChain->SetBranchAddress("pf_mus_px", &pf_mus_px, &b_pf_mus_px);
   fChain->SetBranchAddress("pf_mus_py", &pf_mus_py, &b_pf_mus_py);
   fChain->SetBranchAddress("pf_mus_pz", &pf_mus_pz, &b_pf_mus_pz);
   fChain->SetBranchAddress("pf_mus_status", &pf_mus_status, &b_pf_mus_status);
   fChain->SetBranchAddress("pf_mus_theta", &pf_mus_theta, &b_pf_mus_theta);
   fChain->SetBranchAddress("pf_mus_gen_id", &pf_mus_gen_id, &b_pf_mus_gen_id);
   fChain->SetBranchAddress("pf_mus_gen_phi", &pf_mus_gen_phi, &b_pf_mus_gen_phi);
   fChain->SetBranchAddress("pf_mus_gen_pt", &pf_mus_gen_pt, &b_pf_mus_gen_pt);
   fChain->SetBranchAddress("pf_mus_gen_pz", &pf_mus_gen_pz, &b_pf_mus_gen_pz);
   fChain->SetBranchAddress("pf_mus_gen_px", &pf_mus_gen_px, &b_pf_mus_gen_px);
   fChain->SetBranchAddress("pf_mus_gen_py", &pf_mus_gen_py, &b_pf_mus_gen_py);
   fChain->SetBranchAddress("pf_mus_gen_eta", &pf_mus_gen_eta, &b_pf_mus_gen_eta);
   fChain->SetBranchAddress("pf_mus_gen_theta", &pf_mus_gen_theta, &b_pf_mus_gen_theta);
   fChain->SetBranchAddress("pf_mus_gen_et", &pf_mus_gen_et, &b_pf_mus_gen_et);
   fChain->SetBranchAddress("pf_mus_gen_mother_id", &pf_mus_gen_mother_id, &b_pf_mus_gen_mother_id);
   fChain->SetBranchAddress("pf_mus_gen_mother_phi", &pf_mus_gen_mother_phi, &b_pf_mus_gen_mother_phi);
   fChain->SetBranchAddress("pf_mus_gen_mother_pt", &pf_mus_gen_mother_pt, &b_pf_mus_gen_mother_pt);
   fChain->SetBranchAddress("pf_mus_gen_mother_pz", &pf_mus_gen_mother_pz, &b_pf_mus_gen_mother_pz);
   fChain->SetBranchAddress("pf_mus_gen_mother_px", &pf_mus_gen_mother_px, &b_pf_mus_gen_mother_px);
   fChain->SetBranchAddress("pf_mus_gen_mother_py", &pf_mus_gen_mother_py, &b_pf_mus_gen_mother_py);
   fChain->SetBranchAddress("pf_mus_gen_mother_eta", &pf_mus_gen_mother_eta, &b_pf_mus_gen_mother_eta);
   fChain->SetBranchAddress("pf_mus_gen_mother_theta", &pf_mus_gen_mother_theta, &b_pf_mus_gen_mother_theta);
   fChain->SetBranchAddress("pf_mus_gen_mother_et", &pf_mus_gen_mother_et, &b_pf_mus_gen_mother_et);
   fChain->SetBranchAddress("pf_mus_tkHits", &pf_mus_tkHits, &b_pf_mus_tkHits);
   fChain->SetBranchAddress("pf_mus_cIso", &pf_mus_cIso, &b_pf_mus_cIso);
   fChain->SetBranchAddress("pf_mus_tIso", &pf_mus_tIso, &b_pf_mus_tIso);
   fChain->SetBranchAddress("pf_mus_ecalIso", &pf_mus_ecalIso, &b_pf_mus_ecalIso);
   fChain->SetBranchAddress("pf_mus_hcalIso", &pf_mus_hcalIso, &b_pf_mus_hcalIso);
   fChain->SetBranchAddress("pf_mus_iso03_emVetoEt", &pf_mus_iso03_emVetoEt, &b_pf_mus_iso03_emVetoEt);
   fChain->SetBranchAddress("pf_mus_iso03_hadVetoEt", &pf_mus_iso03_hadVetoEt, &b_pf_mus_iso03_hadVetoEt);
   fChain->SetBranchAddress("pf_mus_calEnergyEm", &pf_mus_calEnergyEm, &b_pf_mus_calEnergyEm);
   fChain->SetBranchAddress("pf_mus_calEnergyHad", &pf_mus_calEnergyHad, &b_pf_mus_calEnergyHad);
   fChain->SetBranchAddress("pf_mus_calEnergyHo", &pf_mus_calEnergyHo, &b_pf_mus_calEnergyHo);
   fChain->SetBranchAddress("pf_mus_calEnergyEmS9", &pf_mus_calEnergyEmS9, &b_pf_mus_calEnergyEmS9);
   fChain->SetBranchAddress("pf_mus_calEnergyHadS9", &pf_mus_calEnergyHadS9, &b_pf_mus_calEnergyHadS9);
   fChain->SetBranchAddress("pf_mus_calEnergyHoS9", &pf_mus_calEnergyHoS9, &b_pf_mus_calEnergyHoS9);
   fChain->SetBranchAddress("pf_mus_iso03_sumPt", &pf_mus_iso03_sumPt, &b_pf_mus_iso03_sumPt);
   fChain->SetBranchAddress("pf_mus_iso03_emEt", &pf_mus_iso03_emEt, &b_pf_mus_iso03_emEt);
   fChain->SetBranchAddress("pf_mus_iso03_hadEt", &pf_mus_iso03_hadEt, &b_pf_mus_iso03_hadEt);
   fChain->SetBranchAddress("pf_mus_iso03_hoEt", &pf_mus_iso03_hoEt, &b_pf_mus_iso03_hoEt);
   fChain->SetBranchAddress("pf_mus_iso03_nTracks", &pf_mus_iso03_nTracks, &b_pf_mus_iso03_nTracks);
   fChain->SetBranchAddress("pf_mus_iso05_sumPt", &pf_mus_iso05_sumPt, &b_pf_mus_iso05_sumPt);
   fChain->SetBranchAddress("pf_mus_iso05_emEt", &pf_mus_iso05_emEt, &b_pf_mus_iso05_emEt);
   fChain->SetBranchAddress("pf_mus_iso05_hadEt", &pf_mus_iso05_hadEt, &b_pf_mus_iso05_hadEt);
   fChain->SetBranchAddress("pf_mus_iso05_hoEt", &pf_mus_iso05_hoEt, &b_pf_mus_iso05_hoEt);
   fChain->SetBranchAddress("pf_mus_iso05_nTracks", &pf_mus_iso05_nTracks, &b_pf_mus_iso05_nTracks);
   fChain->SetBranchAddress("pf_mus_neutralHadronIso", &pf_mus_neutralHadronIso, &b_pf_mus_neutralHadronIso);
   fChain->SetBranchAddress("pf_mus_chargedHadronIso", &pf_mus_chargedHadronIso, &b_pf_mus_chargedHadronIso);
   fChain->SetBranchAddress("pf_mus_photonIso", &pf_mus_photonIso, &b_pf_mus_photonIso);
   fChain->SetBranchAddress("pf_mus_charge", &pf_mus_charge, &b_pf_mus_charge);
   fChain->SetBranchAddress("pf_mus_cm_chi2", &pf_mus_cm_chi2, &b_pf_mus_cm_chi2);
   fChain->SetBranchAddress("pf_mus_cm_ndof", &pf_mus_cm_ndof, &b_pf_mus_cm_ndof);
   fChain->SetBranchAddress("pf_mus_cm_chg", &pf_mus_cm_chg, &b_pf_mus_cm_chg);
   fChain->SetBranchAddress("pf_mus_cm_pt", &pf_mus_cm_pt, &b_pf_mus_cm_pt);
   fChain->SetBranchAddress("pf_mus_cm_px", &pf_mus_cm_px, &b_pf_mus_cm_px);
   fChain->SetBranchAddress("pf_mus_cm_py", &pf_mus_cm_py, &b_pf_mus_cm_py);
   fChain->SetBranchAddress("pf_mus_cm_pz", &pf_mus_cm_pz, &b_pf_mus_cm_pz);
   fChain->SetBranchAddress("pf_mus_cm_eta", &pf_mus_cm_eta, &b_pf_mus_cm_eta);
   fChain->SetBranchAddress("pf_mus_cm_phi", &pf_mus_cm_phi, &b_pf_mus_cm_phi);
   fChain->SetBranchAddress("pf_mus_cm_theta", &pf_mus_cm_theta, &b_pf_mus_cm_theta);
   fChain->SetBranchAddress("pf_mus_cm_d0dum", &pf_mus_cm_d0dum, &b_pf_mus_cm_d0dum);
   fChain->SetBranchAddress("pf_mus_cm_dz", &pf_mus_cm_dz, &b_pf_mus_cm_dz);
   fChain->SetBranchAddress("pf_mus_cm_vx", &pf_mus_cm_vx, &b_pf_mus_cm_vx);
   fChain->SetBranchAddress("pf_mus_cm_vy", &pf_mus_cm_vy, &b_pf_mus_cm_vy);
   fChain->SetBranchAddress("pf_mus_cm_vz", &pf_mus_cm_vz, &b_pf_mus_cm_vz);
   fChain->SetBranchAddress("pf_mus_cm_numvalhits", &pf_mus_cm_numvalhits, &b_pf_mus_cm_numvalhits);
   fChain->SetBranchAddress("pf_mus_cm_numlosthits", &pf_mus_cm_numlosthits, &b_pf_mus_cm_numlosthits);
   fChain->SetBranchAddress("pf_mus_cm_numvalMuonhits", &pf_mus_cm_numvalMuonhits, &b_pf_mus_cm_numvalMuonhits);
   fChain->SetBranchAddress("pf_mus_cm_d0dumErr", &pf_mus_cm_d0dumErr, &b_pf_mus_cm_d0dumErr);
   fChain->SetBranchAddress("pf_mus_cm_dzErr", &pf_mus_cm_dzErr, &b_pf_mus_cm_dzErr);
   fChain->SetBranchAddress("pf_mus_cm_ptErr", &pf_mus_cm_ptErr, &b_pf_mus_cm_ptErr);
   fChain->SetBranchAddress("pf_mus_cm_etaErr", &pf_mus_cm_etaErr, &b_pf_mus_cm_etaErr);
   fChain->SetBranchAddress("pf_mus_cm_phiErr", &pf_mus_cm_phiErr, &b_pf_mus_cm_phiErr);
   fChain->SetBranchAddress("pf_mus_tk_id", &pf_mus_tk_id, &b_pf_mus_tk_id);
   fChain->SetBranchAddress("pf_mus_tk_chi2", &pf_mus_tk_chi2, &b_pf_mus_tk_chi2);
   fChain->SetBranchAddress("pf_mus_tk_ndof", &pf_mus_tk_ndof, &b_pf_mus_tk_ndof);
   fChain->SetBranchAddress("pf_mus_tk_chg", &pf_mus_tk_chg, &b_pf_mus_tk_chg);
   fChain->SetBranchAddress("pf_mus_tk_pt", &pf_mus_tk_pt, &b_pf_mus_tk_pt);
   fChain->SetBranchAddress("pf_mus_tk_px", &pf_mus_tk_px, &b_pf_mus_tk_px);
   fChain->SetBranchAddress("pf_mus_tk_py", &pf_mus_tk_py, &b_pf_mus_tk_py);
   fChain->SetBranchAddress("pf_mus_tk_pz", &pf_mus_tk_pz, &b_pf_mus_tk_pz);
   fChain->SetBranchAddress("pf_mus_tk_eta", &pf_mus_tk_eta, &b_pf_mus_tk_eta);
   fChain->SetBranchAddress("pf_mus_tk_phi", &pf_mus_tk_phi, &b_pf_mus_tk_phi);
   fChain->SetBranchAddress("pf_mus_tk_theta", &pf_mus_tk_theta, &b_pf_mus_tk_theta);
   fChain->SetBranchAddress("pf_mus_tk_d0dum", &pf_mus_tk_d0dum, &b_pf_mus_tk_d0dum);
   fChain->SetBranchAddress("pf_mus_tk_dz", &pf_mus_tk_dz, &b_pf_mus_tk_dz);
   fChain->SetBranchAddress("pf_mus_tk_vx", &pf_mus_tk_vx, &b_pf_mus_tk_vx);
   fChain->SetBranchAddress("pf_mus_tk_vy", &pf_mus_tk_vy, &b_pf_mus_tk_vy);
   fChain->SetBranchAddress("pf_mus_tk_vz", &pf_mus_tk_vz, &b_pf_mus_tk_vz);
   fChain->SetBranchAddress("pf_mus_tk_numvalhits", &pf_mus_tk_numvalhits, &b_pf_mus_tk_numvalhits);
   fChain->SetBranchAddress("pf_mus_tk_numlosthits", &pf_mus_tk_numlosthits, &b_pf_mus_tk_numlosthits);
   fChain->SetBranchAddress("pf_mus_tk_d0dumErr", &pf_mus_tk_d0dumErr, &b_pf_mus_tk_d0dumErr);
   fChain->SetBranchAddress("pf_mus_tk_dzErr", &pf_mus_tk_dzErr, &b_pf_mus_tk_dzErr);
   fChain->SetBranchAddress("pf_mus_tk_ptErr", &pf_mus_tk_ptErr, &b_pf_mus_tk_ptErr);
   fChain->SetBranchAddress("pf_mus_tk_etaErr", &pf_mus_tk_etaErr, &b_pf_mus_tk_etaErr);
   fChain->SetBranchAddress("pf_mus_tk_phiErr", &pf_mus_tk_phiErr, &b_pf_mus_tk_phiErr);
   fChain->SetBranchAddress("pf_mus_tk_numvalPixelhits", &pf_mus_tk_numvalPixelhits, &b_pf_mus_tk_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_tk_numpixelWthMeasr", &pf_mus_tk_numpixelWthMeasr, &b_pf_mus_tk_numpixelWthMeasr);
   fChain->SetBranchAddress("pf_mus_stamu_chi2", &pf_mus_stamu_chi2, &b_pf_mus_stamu_chi2);
   fChain->SetBranchAddress("pf_mus_stamu_ndof", &pf_mus_stamu_ndof, &b_pf_mus_stamu_ndof);
   fChain->SetBranchAddress("pf_mus_stamu_chg", &pf_mus_stamu_chg, &b_pf_mus_stamu_chg);
   fChain->SetBranchAddress("pf_mus_stamu_pt", &pf_mus_stamu_pt, &b_pf_mus_stamu_pt);
   fChain->SetBranchAddress("pf_mus_stamu_px", &pf_mus_stamu_px, &b_pf_mus_stamu_px);
   fChain->SetBranchAddress("pf_mus_stamu_py", &pf_mus_stamu_py, &b_pf_mus_stamu_py);
   fChain->SetBranchAddress("pf_mus_stamu_pz", &pf_mus_stamu_pz, &b_pf_mus_stamu_pz);
   fChain->SetBranchAddress("pf_mus_stamu_eta", &pf_mus_stamu_eta, &b_pf_mus_stamu_eta);
   fChain->SetBranchAddress("pf_mus_stamu_phi", &pf_mus_stamu_phi, &b_pf_mus_stamu_phi);
   fChain->SetBranchAddress("pf_mus_stamu_theta", &pf_mus_stamu_theta, &b_pf_mus_stamu_theta);
   fChain->SetBranchAddress("pf_mus_stamu_d0dum", &pf_mus_stamu_d0dum, &b_pf_mus_stamu_d0dum);
   fChain->SetBranchAddress("pf_mus_stamu_dz", &pf_mus_stamu_dz, &b_pf_mus_stamu_dz);
   fChain->SetBranchAddress("pf_mus_stamu_vx", &pf_mus_stamu_vx, &b_pf_mus_stamu_vx);
   fChain->SetBranchAddress("pf_mus_stamu_vy", &pf_mus_stamu_vy, &b_pf_mus_stamu_vy);
   fChain->SetBranchAddress("pf_mus_stamu_vz", &pf_mus_stamu_vz, &b_pf_mus_stamu_vz);
   fChain->SetBranchAddress("pf_mus_stamu_numvalhits", &pf_mus_stamu_numvalhits, &b_pf_mus_stamu_numvalhits);
   fChain->SetBranchAddress("pf_mus_stamu_numlosthits", &pf_mus_stamu_numlosthits, &b_pf_mus_stamu_numlosthits);
   fChain->SetBranchAddress("pf_mus_stamu_d0dumErr", &pf_mus_stamu_d0dumErr, &b_pf_mus_stamu_d0dumErr);
   fChain->SetBranchAddress("pf_mus_stamu_dzErr", &pf_mus_stamu_dzErr, &b_pf_mus_stamu_dzErr);
   fChain->SetBranchAddress("pf_mus_stamu_ptErr", &pf_mus_stamu_ptErr, &b_pf_mus_stamu_ptErr);
   fChain->SetBranchAddress("pf_mus_stamu_etaErr", &pf_mus_stamu_etaErr, &b_pf_mus_stamu_etaErr);
   fChain->SetBranchAddress("pf_mus_stamu_phiErr", &pf_mus_stamu_phiErr, &b_pf_mus_stamu_phiErr);
   fChain->SetBranchAddress("pf_mus_num_matches", &pf_mus_num_matches, &b_pf_mus_num_matches);
   fChain->SetBranchAddress("pf_mus_isTrackerMuon", &pf_mus_isTrackerMuon, &b_pf_mus_isTrackerMuon);
   fChain->SetBranchAddress("pf_mus_isStandAloneMuon", &pf_mus_isStandAloneMuon, &b_pf_mus_isStandAloneMuon);
   fChain->SetBranchAddress("pf_mus_isCaloMuon", &pf_mus_isCaloMuon, &b_pf_mus_isCaloMuon);
   fChain->SetBranchAddress("pf_mus_isGlobalMuon", &pf_mus_isGlobalMuon, &b_pf_mus_isGlobalMuon);
   fChain->SetBranchAddress("pf_mus_isElectron", &pf_mus_isElectron, &b_pf_mus_isElectron);
   fChain->SetBranchAddress("pf_mus_isConvertedPhoton", &pf_mus_isConvertedPhoton, &b_pf_mus_isConvertedPhoton);
   fChain->SetBranchAddress("pf_mus_isPhoton", &pf_mus_isPhoton, &b_pf_mus_isPhoton);
   fChain->SetBranchAddress("pf_mus_id_All", &pf_mus_id_All, &b_pf_mus_id_All);
   fChain->SetBranchAddress("pf_mus_id_AllGlobalMuons", &pf_mus_id_AllGlobalMuons, &b_pf_mus_id_AllGlobalMuons);
   fChain->SetBranchAddress("pf_mus_id_AllStandAloneMuons", &pf_mus_id_AllStandAloneMuons, &b_pf_mus_id_AllStandAloneMuons);
   fChain->SetBranchAddress("pf_mus_id_AllTrackerMuons", &pf_mus_id_AllTrackerMuons, &b_pf_mus_id_AllTrackerMuons);
   fChain->SetBranchAddress("pf_mus_id_TrackerMuonArbitrated", &pf_mus_id_TrackerMuonArbitrated, &b_pf_mus_id_TrackerMuonArbitrated);
   fChain->SetBranchAddress("pf_mus_id_AllArbitrated", &pf_mus_id_AllArbitrated, &b_pf_mus_id_AllArbitrated);
   fChain->SetBranchAddress("pf_mus_id_GlobalMuonPromptTight", &pf_mus_id_GlobalMuonPromptTight, &b_pf_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationLoose", &pf_mus_id_TMLastStationLoose, &b_pf_mus_id_TMLastStationLoose);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationTight", &pf_mus_id_TMLastStationTight, &b_pf_mus_id_TMLastStationTight);
   fChain->SetBranchAddress("pf_mus_id_TM2DCompatibilityLoose", &pf_mus_id_TM2DCompatibilityLoose, &b_pf_mus_id_TM2DCompatibilityLoose);
   fChain->SetBranchAddress("pf_mus_id_TM2DCompatibilityTight", &pf_mus_id_TM2DCompatibilityTight, &b_pf_mus_id_TM2DCompatibilityTight);
   fChain->SetBranchAddress("pf_mus_id_TMOneStationLoose", &pf_mus_id_TMOneStationLoose, &b_pf_mus_id_TMOneStationLoose);
   fChain->SetBranchAddress("pf_mus_id_TMOneStationTight", &pf_mus_id_TMOneStationTight, &b_pf_mus_id_TMOneStationTight);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationOptimizedLowPtLoose", &pf_mus_id_TMLastStationOptimizedLowPtLoose, &b_pf_mus_id_TMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("pf_mus_id_TMLastStationOptimizedLowPtTight", &pf_mus_id_TMLastStationOptimizedLowPtTight, &b_pf_mus_id_TMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("pf_mus_tk_LayersWithMeasurement", &pf_mus_tk_LayersWithMeasurement, &b_pf_mus_tk_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tk_PixelLayersWithMeasurement", &pf_mus_tk_PixelLayersWithMeasurement, &b_pf_mus_tk_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tk_ValidStripLayersWithMonoAndStereoHit", &pf_mus_tk_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_tk_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_tk_LayersWithoutMeasurement", &pf_mus_tk_LayersWithoutMeasurement, &b_pf_mus_tk_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_tk_ExpectedHitsInner", &pf_mus_tk_ExpectedHitsInner, &b_pf_mus_tk_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_tk_ExpectedHitsOuter", &pf_mus_tk_ExpectedHitsOuter, &b_pf_mus_tk_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_cm_LayersWithMeasurement", &pf_mus_cm_LayersWithMeasurement, &b_pf_mus_cm_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_cm_PixelLayersWithMeasurement", &pf_mus_cm_PixelLayersWithMeasurement, &b_pf_mus_cm_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_cm_ValidStripLayersWithMonoAndStereoHit", &pf_mus_cm_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_cm_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_cm_LayersWithoutMeasurement", &pf_mus_cm_LayersWithoutMeasurement, &b_pf_mus_cm_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_cm_ExpectedHitsInner", &pf_mus_cm_ExpectedHitsInner, &b_pf_mus_cm_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_cm_ExpectedHitsOuter", &pf_mus_cm_ExpectedHitsOuter, &b_pf_mus_cm_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_picky_LayersWithMeasurement", &pf_mus_picky_LayersWithMeasurement, &b_pf_mus_picky_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_picky_PixelLayersWithMeasurement", &pf_mus_picky_PixelLayersWithMeasurement, &b_pf_mus_picky_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_picky_ValidStripLayersWithMonoAndStereoHit", &pf_mus_picky_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_picky_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_picky_LayersWithoutMeasurement", &pf_mus_picky_LayersWithoutMeasurement, &b_pf_mus_picky_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_picky_ExpectedHitsInner", &pf_mus_picky_ExpectedHitsInner, &b_pf_mus_picky_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_picky_ExpectedHitsOuter", &pf_mus_picky_ExpectedHitsOuter, &b_pf_mus_picky_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_tpfms_LayersWithMeasurement", &pf_mus_tpfms_LayersWithMeasurement, &b_pf_mus_tpfms_LayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tpfms_PixelLayersWithMeasurement", &pf_mus_tpfms_PixelLayersWithMeasurement, &b_pf_mus_tpfms_PixelLayersWithMeasurement);
   fChain->SetBranchAddress("pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit", &pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit, &b_pf_mus_tpfms_ValidStripLayersWithMonoAndStereoHit);
   fChain->SetBranchAddress("pf_mus_tpfms_LayersWithoutMeasurement", &pf_mus_tpfms_LayersWithoutMeasurement, &b_pf_mus_tpfms_LayersWithoutMeasurement);
   fChain->SetBranchAddress("pf_mus_tpfms_ExpectedHitsInner", &pf_mus_tpfms_ExpectedHitsInner, &b_pf_mus_tpfms_ExpectedHitsInner);
   fChain->SetBranchAddress("pf_mus_tpfms_ExpectedHitsOuter", &pf_mus_tpfms_ExpectedHitsOuter, &b_pf_mus_tpfms_ExpectedHitsOuter);
   fChain->SetBranchAddress("pf_mus_picky_id", &pf_mus_picky_id, &b_pf_mus_picky_id);
   fChain->SetBranchAddress("pf_mus_picky_chi2", &pf_mus_picky_chi2, &b_pf_mus_picky_chi2);
   fChain->SetBranchAddress("pf_mus_picky_ndof", &pf_mus_picky_ndof, &b_pf_mus_picky_ndof);
   fChain->SetBranchAddress("pf_mus_picky_chg", &pf_mus_picky_chg, &b_pf_mus_picky_chg);
   fChain->SetBranchAddress("pf_mus_picky_pt", &pf_mus_picky_pt, &b_pf_mus_picky_pt);
   fChain->SetBranchAddress("pf_mus_picky_px", &pf_mus_picky_px, &b_pf_mus_picky_px);
   fChain->SetBranchAddress("pf_mus_picky_py", &pf_mus_picky_py, &b_pf_mus_picky_py);
   fChain->SetBranchAddress("pf_mus_picky_pz", &pf_mus_picky_pz, &b_pf_mus_picky_pz);
   fChain->SetBranchAddress("pf_mus_picky_eta", &pf_mus_picky_eta, &b_pf_mus_picky_eta);
   fChain->SetBranchAddress("pf_mus_picky_phi", &pf_mus_picky_phi, &b_pf_mus_picky_phi);
   fChain->SetBranchAddress("pf_mus_picky_theta", &pf_mus_picky_theta, &b_pf_mus_picky_theta);
   fChain->SetBranchAddress("pf_mus_picky_d0dum", &pf_mus_picky_d0dum, &b_pf_mus_picky_d0dum);
   fChain->SetBranchAddress("pf_mus_picky_dz", &pf_mus_picky_dz, &b_pf_mus_picky_dz);
   fChain->SetBranchAddress("pf_mus_picky_vx", &pf_mus_picky_vx, &b_pf_mus_picky_vx);
   fChain->SetBranchAddress("pf_mus_picky_vy", &pf_mus_picky_vy, &b_pf_mus_picky_vy);
   fChain->SetBranchAddress("pf_mus_picky_vz", &pf_mus_picky_vz, &b_pf_mus_picky_vz);
   fChain->SetBranchAddress("pf_mus_picky_numvalhits", &pf_mus_picky_numvalhits, &b_pf_mus_picky_numvalhits);
   fChain->SetBranchAddress("pf_mus_picky_numlosthits", &pf_mus_picky_numlosthits, &b_pf_mus_picky_numlosthits);
   fChain->SetBranchAddress("pf_mus_picky_d0dumErr", &pf_mus_picky_d0dumErr, &b_pf_mus_picky_d0dumErr);
   fChain->SetBranchAddress("pf_mus_picky_dzErr", &pf_mus_picky_dzErr, &b_pf_mus_picky_dzErr);
   fChain->SetBranchAddress("pf_mus_picky_ptErr", &pf_mus_picky_ptErr, &b_pf_mus_picky_ptErr);
   fChain->SetBranchAddress("pf_mus_picky_etaErr", &pf_mus_picky_etaErr, &b_pf_mus_picky_etaErr);
   fChain->SetBranchAddress("pf_mus_picky_phiErr", &pf_mus_picky_phiErr, &b_pf_mus_picky_phiErr);
   fChain->SetBranchAddress("pf_mus_picky_numvalPixelhits", &pf_mus_picky_numvalPixelhits, &b_pf_mus_picky_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_tpfms_id", &pf_mus_tpfms_id, &b_pf_mus_tpfms_id);
   fChain->SetBranchAddress("pf_mus_tpfms_chi2", &pf_mus_tpfms_chi2, &b_pf_mus_tpfms_chi2);
   fChain->SetBranchAddress("pf_mus_tpfms_ndof", &pf_mus_tpfms_ndof, &b_pf_mus_tpfms_ndof);
   fChain->SetBranchAddress("pf_mus_tpfms_chg", &pf_mus_tpfms_chg, &b_pf_mus_tpfms_chg);
   fChain->SetBranchAddress("pf_mus_tpfms_pt", &pf_mus_tpfms_pt, &b_pf_mus_tpfms_pt);
   fChain->SetBranchAddress("pf_mus_tpfms_px", &pf_mus_tpfms_px, &b_pf_mus_tpfms_px);
   fChain->SetBranchAddress("pf_mus_tpfms_py", &pf_mus_tpfms_py, &b_pf_mus_tpfms_py);
   fChain->SetBranchAddress("pf_mus_tpfms_pz", &pf_mus_tpfms_pz, &b_pf_mus_tpfms_pz);
   fChain->SetBranchAddress("pf_mus_tpfms_eta", &pf_mus_tpfms_eta, &b_pf_mus_tpfms_eta);
   fChain->SetBranchAddress("pf_mus_tpfms_phi", &pf_mus_tpfms_phi, &b_pf_mus_tpfms_phi);
   fChain->SetBranchAddress("pf_mus_tpfms_theta", &pf_mus_tpfms_theta, &b_pf_mus_tpfms_theta);
   fChain->SetBranchAddress("pf_mus_tpfms_d0dum", &pf_mus_tpfms_d0dum, &b_pf_mus_tpfms_d0dum);
   fChain->SetBranchAddress("pf_mus_tpfms_dz", &pf_mus_tpfms_dz, &b_pf_mus_tpfms_dz);
   fChain->SetBranchAddress("pf_mus_tpfms_vx", &pf_mus_tpfms_vx, &b_pf_mus_tpfms_vx);
   fChain->SetBranchAddress("pf_mus_tpfms_vy", &pf_mus_tpfms_vy, &b_pf_mus_tpfms_vy);
   fChain->SetBranchAddress("pf_mus_tpfms_vz", &pf_mus_tpfms_vz, &b_pf_mus_tpfms_vz);
   fChain->SetBranchAddress("pf_mus_tpfms_numvalhits", &pf_mus_tpfms_numvalhits, &b_pf_mus_tpfms_numvalhits);
   fChain->SetBranchAddress("pf_mus_tpfms_numlosthits", &pf_mus_tpfms_numlosthits, &b_pf_mus_tpfms_numlosthits);
   fChain->SetBranchAddress("pf_mus_tpfms_d0dumErr", &pf_mus_tpfms_d0dumErr, &b_pf_mus_tpfms_d0dumErr);
   fChain->SetBranchAddress("pf_mus_tpfms_dzErr", &pf_mus_tpfms_dzErr, &b_pf_mus_tpfms_dzErr);
   fChain->SetBranchAddress("pf_mus_tpfms_ptErr", &pf_mus_tpfms_ptErr, &b_pf_mus_tpfms_ptErr);
   fChain->SetBranchAddress("pf_mus_tpfms_etaErr", &pf_mus_tpfms_etaErr, &b_pf_mus_tpfms_etaErr);
   fChain->SetBranchAddress("pf_mus_tpfms_phiErr", &pf_mus_tpfms_phiErr, &b_pf_mus_tpfms_phiErr);
   fChain->SetBranchAddress("pf_mus_tpfms_numvalPixelhits", &pf_mus_tpfms_numvalPixelhits, &b_pf_mus_tpfms_numvalPixelhits);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumChargedHadronPt", &pf_mus_pfIsolationR03_sumChargedHadronPt, &b_pf_mus_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumChargedParticlePt", &pf_mus_pfIsolationR03_sumChargedParticlePt, &b_pf_mus_pfIsolationR03_sumChargedParticlePt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumNeutralHadronEt", &pf_mus_pfIsolationR03_sumNeutralHadronEt, &b_pf_mus_pfIsolationR03_sumNeutralHadronEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold", &pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_pf_mus_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumPhotonEt", &pf_mus_pfIsolationR03_sumPhotonEt, &b_pf_mus_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumPhotonEtHighThreshold", &pf_mus_pfIsolationR03_sumPhotonEtHighThreshold, &b_pf_mus_pfIsolationR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR03_sumPUPt", &pf_mus_pfIsolationR03_sumPUPt, &b_pf_mus_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumChargedHadronPt", &pf_mus_pfIsolationR04_sumChargedHadronPt, &b_pf_mus_pfIsolationR04_sumChargedHadronPt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumChargedParticlePt", &pf_mus_pfIsolationR04_sumChargedParticlePt, &b_pf_mus_pfIsolationR04_sumChargedParticlePt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumNeutralHadronEt", &pf_mus_pfIsolationR04_sumNeutralHadronEt, &b_pf_mus_pfIsolationR04_sumNeutralHadronEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold", &pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_pf_mus_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumPhotonEt", &pf_mus_pfIsolationR04_sumPhotonEt, &b_pf_mus_pfIsolationR04_sumPhotonEt);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumPhotonEtHighThreshold", &pf_mus_pfIsolationR04_sumPhotonEtHighThreshold, &b_pf_mus_pfIsolationR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("pf_mus_pfIsolationR04_sumPUPt", &pf_mus_pfIsolationR04_sumPUPt, &b_pf_mus_pfIsolationR04_sumPUPt);
   fChain->SetBranchAddress("pf_mus_dB", &pf_mus_dB, &b_pf_mus_dB);
   fChain->SetBranchAddress("pf_mus_numberOfMatchedStations", &pf_mus_numberOfMatchedStations, &b_pf_mus_numberOfMatchedStations);
   fChain->SetBranchAddress("pf_mus_isPFMuon", &pf_mus_isPFMuon, &b_pf_mus_isPFMuon);
   fChain->SetBranchAddress("Npf_photons", &Npf_photons, &b_Npf_photons);
   fChain->SetBranchAddress("pf_photons_energy", &pf_photons_energy, &b_pf_photons_energy);
   fChain->SetBranchAddress("pf_photons_et", &pf_photons_et, &b_pf_photons_et);
   fChain->SetBranchAddress("pf_photons_eta", &pf_photons_eta, &b_pf_photons_eta);
   fChain->SetBranchAddress("pf_photons_phi", &pf_photons_phi, &b_pf_photons_phi);
   fChain->SetBranchAddress("pf_photons_pt", &pf_photons_pt, &b_pf_photons_pt);
   fChain->SetBranchAddress("pf_photons_px", &pf_photons_px, &b_pf_photons_px);
   fChain->SetBranchAddress("pf_photons_py", &pf_photons_py, &b_pf_photons_py);
   fChain->SetBranchAddress("pf_photons_pz", &pf_photons_pz, &b_pf_photons_pz);
   fChain->SetBranchAddress("pf_photons_status", &pf_photons_status, &b_pf_photons_status);
   fChain->SetBranchAddress("pf_photons_theta", &pf_photons_theta, &b_pf_photons_theta);
   fChain->SetBranchAddress("pf_photons_hadOverEM", &pf_photons_hadOverEM, &b_pf_photons_hadOverEM);
   fChain->SetBranchAddress("pf_photons_hadTowOverEM", &pf_photons_hadTowOverEM, &b_pf_photons_hadTowOverEM);
   fChain->SetBranchAddress("pf_photons_scEnergy", &pf_photons_scEnergy, &b_pf_photons_scEnergy);
   fChain->SetBranchAddress("pf_photons_scRawEnergy", &pf_photons_scRawEnergy, &b_pf_photons_scRawEnergy);
   fChain->SetBranchAddress("pf_photons_scEta", &pf_photons_scEta, &b_pf_photons_scEta);
   fChain->SetBranchAddress("pf_photons_scPhi", &pf_photons_scPhi, &b_pf_photons_scPhi);
   fChain->SetBranchAddress("pf_photons_scEtaWidth", &pf_photons_scEtaWidth, &b_pf_photons_scEtaWidth);
   fChain->SetBranchAddress("pf_photons_scPhiWidth", &pf_photons_scPhiWidth, &b_pf_photons_scPhiWidth);
   fChain->SetBranchAddress("pf_photons_isAlsoElectron", &pf_photons_isAlsoElectron, &b_pf_photons_isAlsoElectron);
   fChain->SetBranchAddress("pf_photons_hasPixelSeed", &pf_photons_hasPixelSeed, &b_pf_photons_hasPixelSeed);
   fChain->SetBranchAddress("pf_photons_isConverted", &pf_photons_isConverted, &b_pf_photons_isConverted);
   fChain->SetBranchAddress("pf_photons_isEBGap", &pf_photons_isEBGap, &b_pf_photons_isEBGap);
   fChain->SetBranchAddress("pf_photons_isEEGap", &pf_photons_isEEGap, &b_pf_photons_isEEGap);
   fChain->SetBranchAddress("pf_photons_isEBEEGap", &pf_photons_isEBEEGap, &b_pf_photons_isEBEEGap);
   fChain->SetBranchAddress("pf_photons_isEBPho", &pf_photons_isEBPho, &b_pf_photons_isEBPho);
   fChain->SetBranchAddress("pf_photons_isEEPho", &pf_photons_isEEPho, &b_pf_photons_isEEPho);
   fChain->SetBranchAddress("pf_photons_maxEnergyXtal", &pf_photons_maxEnergyXtal, &b_pf_photons_maxEnergyXtal);
   fChain->SetBranchAddress("pf_photons_e1x5", &pf_photons_e1x5, &b_pf_photons_e1x5);
   fChain->SetBranchAddress("pf_photons_e2x5", &pf_photons_e2x5, &b_pf_photons_e2x5);
   fChain->SetBranchAddress("pf_photons_e3x3", &pf_photons_e3x3, &b_pf_photons_e3x3);
   fChain->SetBranchAddress("pf_photons_e5x5", &pf_photons_e5x5, &b_pf_photons_e5x5);
   fChain->SetBranchAddress("pf_photons_sigmaEtaEta", &pf_photons_sigmaEtaEta, &b_pf_photons_sigmaEtaEta);
   fChain->SetBranchAddress("pf_photons_sigmaIetaIeta", &pf_photons_sigmaIetaIeta, &b_pf_photons_sigmaIetaIeta);
   fChain->SetBranchAddress("pf_photons_r9", &pf_photons_r9, &b_pf_photons_r9);
   fChain->SetBranchAddress("pf_photons_chIso", &pf_photons_chIso, &b_pf_photons_chIso);
   fChain->SetBranchAddress("pf_photons_nhIso", &pf_photons_nhIso, &b_pf_photons_nhIso);
   fChain->SetBranchAddress("pf_photons_phIso", &pf_photons_phIso, &b_pf_photons_phIso);
   fChain->SetBranchAddress("Npfcand", &Npfcand, &b_Npfcand);
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
   fChain->SetBranchAddress("Npfmets", &Npfmets, &b_Npfmets);
   fChain->SetBranchAddress("pfmets_et", &pfmets_et, &b_pfmets_et);
   fChain->SetBranchAddress("pfmets_phi", &pfmets_phi, &b_pfmets_phi);
   fChain->SetBranchAddress("pfmets_ex", &pfmets_ex, &b_pfmets_ex);
   fChain->SetBranchAddress("pfmets_ey", &pfmets_ey, &b_pfmets_ey);
   fChain->SetBranchAddress("pfmets_gen_et", &pfmets_gen_et, &b_pfmets_gen_et);
   fChain->SetBranchAddress("pfmets_gen_phi", &pfmets_gen_phi, &b_pfmets_gen_phi);
   fChain->SetBranchAddress("pfmets_sign", &pfmets_sign, &b_pfmets_sign);
   fChain->SetBranchAddress("pfmets_sumEt", &pfmets_sumEt, &b_pfmets_sumEt);
   fChain->SetBranchAddress("pfmets_unCPhi", &pfmets_unCPhi, &b_pfmets_unCPhi);
   fChain->SetBranchAddress("pfmets_unCPt", &pfmets_unCPt, &b_pfmets_unCPt);
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
   fChain->SetBranchAddress("Ntcmets", &Ntcmets, &b_Ntcmets);
   fChain->SetBranchAddress("tcmets_et", &tcmets_et, &b_tcmets_et);
   fChain->SetBranchAddress("tcmets_phi", &tcmets_phi, &b_tcmets_phi);
   fChain->SetBranchAddress("tcmets_ex", &tcmets_ex, &b_tcmets_ex);
   fChain->SetBranchAddress("tcmets_ey", &tcmets_ey, &b_tcmets_ey);
   fChain->SetBranchAddress("tcmets_sumEt", &tcmets_sumEt, &b_tcmets_sumEt);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("tracks_chi2", &tracks_chi2, &b_tracks_chi2);
   fChain->SetBranchAddress("tracks_ndof", &tracks_ndof, &b_tracks_ndof);
   fChain->SetBranchAddress("tracks_chg", &tracks_chg, &b_tracks_chg);
   fChain->SetBranchAddress("tracks_pt", &tracks_pt, &b_tracks_pt);
   fChain->SetBranchAddress("tracks_px", &tracks_px, &b_tracks_px);
   fChain->SetBranchAddress("tracks_py", &tracks_py, &b_tracks_py);
   fChain->SetBranchAddress("tracks_pz", &tracks_pz, &b_tracks_pz);
   fChain->SetBranchAddress("tracks_eta", &tracks_eta, &b_tracks_eta);
   fChain->SetBranchAddress("tracks_phi", &tracks_phi, &b_tracks_phi);
   fChain->SetBranchAddress("tracks_d0dum", &tracks_d0dum, &b_tracks_d0dum);
   fChain->SetBranchAddress("tracks_dz", &tracks_dz, &b_tracks_dz);
   fChain->SetBranchAddress("tracks_vx", &tracks_vx, &b_tracks_vx);
   fChain->SetBranchAddress("tracks_vy", &tracks_vy, &b_tracks_vy);
   fChain->SetBranchAddress("tracks_vz", &tracks_vz, &b_tracks_vz);
   fChain->SetBranchAddress("tracks_numvalhits", &tracks_numvalhits, &b_tracks_numvalhits);
   fChain->SetBranchAddress("tracks_numlosthits", &tracks_numlosthits, &b_tracks_numlosthits);
   fChain->SetBranchAddress("tracks_d0dumErr", &tracks_d0dumErr, &b_tracks_d0dumErr);
   fChain->SetBranchAddress("tracks_dzErr", &tracks_dzErr, &b_tracks_dzErr);
   fChain->SetBranchAddress("tracks_ptErr", &tracks_ptErr, &b_tracks_ptErr);
   fChain->SetBranchAddress("tracks_etaErr", &tracks_etaErr, &b_tracks_etaErr);
   fChain->SetBranchAddress("tracks_phiErr", &tracks_phiErr, &b_tracks_phiErr);
   fChain->SetBranchAddress("tracks_highPurity", &tracks_highPurity, &b_tracks_highPurity);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("experimentType", &experimentType, &b_experimentType);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
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
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_px", &fastjets_AK5PF_R1p2_R0p5pT10_px, &b_fastjets_AK5PF_R1p2_R0p5pT10_px);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_py", &fastjets_AK5PF_R1p2_R0p5pT10_py, &b_fastjets_AK5PF_R1p2_R0p5pT10_py);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_pz", &fastjets_AK5PF_R1p2_R0p5pT10_pz, &b_fastjets_AK5PF_R1p2_R0p5pT10_pz);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_energy", &fastjets_AK5PF_R1p2_R0p5pT10_energy, &b_fastjets_AK5PF_R1p2_R0p5pT10_energy);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_phi", &fastjets_AK5PF_R1p2_R0p5pT10_phi, &b_fastjets_AK5PF_R1p2_R0p5pT10_phi);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_eta", &fastjets_AK5PF_R1p2_R0p5pT10_eta, &b_fastjets_AK5PF_R1p2_R0p5pT10_eta);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_index", &fastjets_AK5PF_R1p2_R0p5pT10_index, &b_fastjets_AK5PF_R1p2_R0p5pT10_index);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT10_nconstituents", &fastjets_AK5PF_R1p2_R0p5pT10_nconstituents, &b_fastjets_AK5PF_R1p2_R0p5pT10_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_px", &fastjets_AK5PF_R1p2_R0p5pT15_px, &b_fastjets_AK5PF_R1p2_R0p5pT15_px);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_py", &fastjets_AK5PF_R1p2_R0p5pT15_py, &b_fastjets_AK5PF_R1p2_R0p5pT15_py);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_pz", &fastjets_AK5PF_R1p2_R0p5pT15_pz, &b_fastjets_AK5PF_R1p2_R0p5pT15_pz);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_energy", &fastjets_AK5PF_R1p2_R0p5pT15_energy, &b_fastjets_AK5PF_R1p2_R0p5pT15_energy);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_phi", &fastjets_AK5PF_R1p2_R0p5pT15_phi, &b_fastjets_AK5PF_R1p2_R0p5pT15_phi);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_eta", &fastjets_AK5PF_R1p2_R0p5pT15_eta, &b_fastjets_AK5PF_R1p2_R0p5pT15_eta);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_index", &fastjets_AK5PF_R1p2_R0p5pT15_index, &b_fastjets_AK5PF_R1p2_R0p5pT15_index);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT15_nconstituents", &fastjets_AK5PF_R1p2_R0p5pT15_nconstituents, &b_fastjets_AK5PF_R1p2_R0p5pT15_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_px", &fastjets_AK5PF_R1p2_R0p5pT20_px, &b_fastjets_AK5PF_R1p2_R0p5pT20_px);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_py", &fastjets_AK5PF_R1p2_R0p5pT20_py, &b_fastjets_AK5PF_R1p2_R0p5pT20_py);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_pz", &fastjets_AK5PF_R1p2_R0p5pT20_pz, &b_fastjets_AK5PF_R1p2_R0p5pT20_pz);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_energy", &fastjets_AK5PF_R1p2_R0p5pT20_energy, &b_fastjets_AK5PF_R1p2_R0p5pT20_energy);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_phi", &fastjets_AK5PF_R1p2_R0p5pT20_phi, &b_fastjets_AK5PF_R1p2_R0p5pT20_phi);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_eta", &fastjets_AK5PF_R1p2_R0p5pT20_eta, &b_fastjets_AK5PF_R1p2_R0p5pT20_eta);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_index", &fastjets_AK5PF_R1p2_R0p5pT20_index, &b_fastjets_AK5PF_R1p2_R0p5pT20_index);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT20_nconstituents", &fastjets_AK5PF_R1p2_R0p5pT20_nconstituents, &b_fastjets_AK5PF_R1p2_R0p5pT20_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_px", &fastjets_AK5PF_R1p2_R0p5pT25_px, &b_fastjets_AK5PF_R1p2_R0p5pT25_px);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_py", &fastjets_AK5PF_R1p2_R0p5pT25_py, &b_fastjets_AK5PF_R1p2_R0p5pT25_py);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_pz", &fastjets_AK5PF_R1p2_R0p5pT25_pz, &b_fastjets_AK5PF_R1p2_R0p5pT25_pz);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_energy", &fastjets_AK5PF_R1p2_R0p5pT25_energy, &b_fastjets_AK5PF_R1p2_R0p5pT25_energy);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_phi", &fastjets_AK5PF_R1p2_R0p5pT25_phi, &b_fastjets_AK5PF_R1p2_R0p5pT25_phi);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_eta", &fastjets_AK5PF_R1p2_R0p5pT25_eta, &b_fastjets_AK5PF_R1p2_R0p5pT25_eta);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_index", &fastjets_AK5PF_R1p2_R0p5pT25_index, &b_fastjets_AK5PF_R1p2_R0p5pT25_index);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT25_nconstituents", &fastjets_AK5PF_R1p2_R0p5pT25_nconstituents, &b_fastjets_AK5PF_R1p2_R0p5pT25_nconstituents);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_px", &fastjets_AK5PF_R1p2_R0p5pT30_px, &b_fastjets_AK5PF_R1p2_R0p5pT30_px);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_py", &fastjets_AK5PF_R1p2_R0p5pT30_py, &b_fastjets_AK5PF_R1p2_R0p5pT30_py);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_pz", &fastjets_AK5PF_R1p2_R0p5pT30_pz, &b_fastjets_AK5PF_R1p2_R0p5pT30_pz);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_energy", &fastjets_AK5PF_R1p2_R0p5pT30_energy, &b_fastjets_AK5PF_R1p2_R0p5pT30_energy);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_phi", &fastjets_AK5PF_R1p2_R0p5pT30_phi, &b_fastjets_AK5PF_R1p2_R0p5pT30_phi);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_eta", &fastjets_AK5PF_R1p2_R0p5pT30_eta, &b_fastjets_AK5PF_R1p2_R0p5pT30_eta);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_index", &fastjets_AK5PF_R1p2_R0p5pT30_index, &b_fastjets_AK5PF_R1p2_R0p5pT30_index);
   fChain->SetBranchAddress("fastjets_AK5PF_R1p2_R0p5pT30_nconstituents", &fastjets_AK5PF_R1p2_R0p5pT30_nconstituents, &b_fastjets_AK5PF_R1p2_R0p5pT30_nconstituents);

}


void InitializeA(TChain *fChain) 
{
   // Set object pointer
   trigger_prescalevalue = 0;
   trigger_name = 0;
   trigger_decision = 0;
   trigger_lastfiltername = 0;
   triggerobject_pt = 0;
   triggerobject_px = 0;
   triggerobject_py = 0;
   triggerobject_pz = 0;
   triggerobject_et = 0;
   triggerobject_energy = 0;
   triggerobject_phi = 0;
   triggerobject_eta = 0;
   triggerobject_collectionname = 0;
   standalone_triggerobject_pt = 0;
   standalone_triggerobject_px = 0;
   standalone_triggerobject_py = 0;
   standalone_triggerobject_pz = 0;
   standalone_triggerobject_et = 0;
   standalone_triggerobject_energy = 0;
   standalone_triggerobject_phi = 0;
   standalone_triggerobject_eta = 0;
   standalone_triggerobject_collectionname = 0;
   L1trigger_bit = 0;
   L1trigger_techTrigger = 0;
   L1trigger_prescalevalue = 0;
   L1trigger_name = 0;
   L1trigger_alias = 0;
   L1trigger_decision = 0;
   L1trigger_decision_nomask = 0;
   els_conversion_dist = 0;
   els_conversion_dcot = 0;
   els_PFchargedHadronIsoR03 = 0;
   els_PFphotonIsoR03 = 0;
   els_PFneutralHadronIsoR03 = 0;
   els_hasMatchedConversion = 0;
   pf_els_PFchargedHadronIsoR03 = 0;
   pf_els_PFphotonIsoR03 = 0;
   pf_els_PFneutralHadronIsoR03 = 0;
   pf_els_hasMatchedConversion = 0;
   jets_AK5PFclean_corrL2L3 = 0;
   jets_AK5PFclean_corrL2L3Residual = 0;
   jets_AK5PFclean_corrL1FastL2L3 = 0;
   jets_AK5PFclean_corrL1L2L3 = 0;
   jets_AK5PFclean_corrL1FastL2L3Residual = 0;
   jets_AK5PFclean_corrL1L2L3Residual = 0;
   jets_AK5PFclean_Uncert = 0;
   PU_zpositions = 0;
   PU_sumpT_lowpT = 0;
   PU_sumpT_highpT = 0;
   PU_ntrks_lowpT = 0;
   PU_ntrks_highpT = 0;
   PU_NumInteractions = 0;
   PU_bunchCrossing = 0;
   PU_TrueNumInteractions = 0;
   pdfweights_cteq = 0;
   pdfweights_mstw = 0;
   pdfweights_nnpdf = 0;
   photon_chIsoValues = 0;
   photon_phIsoValues = 0;
   photon_nhIsoValues = 0;
   photon_passElectronVeto = 0;
   puJet_rejectionBeta = 0;
   puJet_rejectionMVA = 0;
   isotk_pt = 0;
   isotk_phi = 0;
   isotk_eta = 0;
   isotk_iso = 0;
   isotk_dzpv = 0;
   isotk_charge = 0;

   fChain->SetBranchAddress("trigger_prescalevalue", &trigger_prescalevalue, &b_trigger_prescalevalue);
   fChain->SetBranchAddress("trigger_name", &trigger_name, &b_trigger_name);
   fChain->SetBranchAddress("trigger_decision", &trigger_decision, &b_trigger_decision);
   fChain->SetBranchAddress("trigger_lastfiltername", &trigger_lastfiltername, &b_trigger_lastfiltername);
   fChain->SetBranchAddress("triggerobject_pt", &triggerobject_pt, &b_triggerobject_pt);
   fChain->SetBranchAddress("triggerobject_px", &triggerobject_px, &b_triggerobject_px);
   fChain->SetBranchAddress("triggerobject_py", &triggerobject_py, &b_triggerobject_py);
   fChain->SetBranchAddress("triggerobject_pz", &triggerobject_pz, &b_triggerobject_pz);
   fChain->SetBranchAddress("triggerobject_et", &triggerobject_et, &b_triggerobject_et);
   fChain->SetBranchAddress("triggerobject_energy", &triggerobject_energy, &b_triggerobject_energy);
   fChain->SetBranchAddress("triggerobject_phi", &triggerobject_phi, &b_triggerobject_phi);
   fChain->SetBranchAddress("triggerobject_eta", &triggerobject_eta, &b_triggerobject_eta);
   fChain->SetBranchAddress("triggerobject_collectionname", &triggerobject_collectionname, &b_triggerobject_collectionname);
   fChain->SetBranchAddress("standalone_triggerobject_pt", &standalone_triggerobject_pt, &b_standalone_triggerobject_pt);
   fChain->SetBranchAddress("standalone_triggerobject_px", &standalone_triggerobject_px, &b_standalone_triggerobject_px);
   fChain->SetBranchAddress("standalone_triggerobject_py", &standalone_triggerobject_py, &b_standalone_triggerobject_py);
   fChain->SetBranchAddress("standalone_triggerobject_pz", &standalone_triggerobject_pz, &b_standalone_triggerobject_pz);
   fChain->SetBranchAddress("standalone_triggerobject_et", &standalone_triggerobject_et, &b_standalone_triggerobject_et);
   fChain->SetBranchAddress("standalone_triggerobject_energy", &standalone_triggerobject_energy, &b_standalone_triggerobject_energy);
   fChain->SetBranchAddress("standalone_triggerobject_phi", &standalone_triggerobject_phi, &b_standalone_triggerobject_phi);
   fChain->SetBranchAddress("standalone_triggerobject_eta", &standalone_triggerobject_eta, &b_standalone_triggerobject_eta);
   fChain->SetBranchAddress("standalone_triggerobject_collectionname", &standalone_triggerobject_collectionname, &b_standalone_triggerobject_collectionname);
   fChain->SetBranchAddress("L1trigger_bit", &L1trigger_bit, &b_L1trigger_bit);
   fChain->SetBranchAddress("L1trigger_techTrigger", &L1trigger_techTrigger, &b_L1trigger_techTrigger);
   fChain->SetBranchAddress("L1trigger_prescalevalue", &L1trigger_prescalevalue, &b_L1trigger_prescalevalue);
   fChain->SetBranchAddress("L1trigger_name", &L1trigger_name, &b_L1trigger_name);
   fChain->SetBranchAddress("L1trigger_alias", &L1trigger_alias, &b_L1trigger_alias);
   fChain->SetBranchAddress("L1trigger_decision", &L1trigger_decision, &b_L1trigger_decision);
   fChain->SetBranchAddress("L1trigger_decision_nomask", &L1trigger_decision_nomask, &b_L1trigger_decision_nomask);
   fChain->SetBranchAddress("els_conversion_dist", &els_conversion_dist, &b_els_conversion_dist);
   fChain->SetBranchAddress("els_conversion_dcot", &els_conversion_dcot, &b_els_conversion_dcot);
   fChain->SetBranchAddress("els_PFchargedHadronIsoR03", &els_PFchargedHadronIsoR03, &b_els_PFchargedHadronIsoR03);
   fChain->SetBranchAddress("els_PFphotonIsoR03", &els_PFphotonIsoR03, &b_els_PFphotonIsoR03);
   fChain->SetBranchAddress("els_PFneutralHadronIsoR03", &els_PFneutralHadronIsoR03, &b_els_PFneutralHadronIsoR03);
   fChain->SetBranchAddress("els_hasMatchedConversion", &els_hasMatchedConversion, &b_els_hasMatchedConversion);
   fChain->SetBranchAddress("pf_els_PFchargedHadronIsoR03", &pf_els_PFchargedHadronIsoR03, &b_pf_els_PFchargedHadronIsoR03);
   fChain->SetBranchAddress("pf_els_PFphotonIsoR03", &pf_els_PFphotonIsoR03, &b_pf_els_PFphotonIsoR03);
   fChain->SetBranchAddress("pf_els_PFneutralHadronIsoR03", &pf_els_PFneutralHadronIsoR03, &b_pf_els_PFneutralHadronIsoR03);
   fChain->SetBranchAddress("pf_els_hasMatchedConversion", &pf_els_hasMatchedConversion, &b_pf_els_hasMatchedConversion);
   fChain->SetBranchAddress("trk_nTOBTEC", &trk_nTOBTEC, &b_Ctrk_nTOBTEC);
   fChain->SetBranchAddress("trk_ratioAllTOBTEC", &trk_ratioAllTOBTEC, &b_trk_ratioAllTOBTEC);
   fChain->SetBranchAddress("trk_ratioJetTOBTEC", &trk_ratioJetTOBTEC, &b_trk_ratioJetTOBTEC);
   fChain->SetBranchAddress("hbhefilter_decision", &hbhefilter_decision, &b_hbhefilter_decision);
   fChain->SetBranchAddress("trackingfailurefilter_decision", &trackingfailurefilter_decision, &b_trackingfailurefilter_decision);
   fChain->SetBranchAddress("cschalofilter_decision", &cschalofilter_decision, &b_cschalofilter_decision);
   fChain->SetBranchAddress("ecalTPfilter_decision", &ecalTPfilter_decision, &b_ecalTPfilter_decision);
   fChain->SetBranchAddress("ecalBEfilter_decision", &ecalBEfilter_decision, &b_ecalBEfilter_decision);
   fChain->SetBranchAddress("scrapingVeto_decision", &scrapingVeto_decision, &b_scrapingVeto_decision);
   fChain->SetBranchAddress("greedymuonfilter_decision", &greedymuonfilter_decision, &b_greedymuonfilter_decision);
   fChain->SetBranchAddress("inconsistentPFmuonfilter_decision", &inconsistentPFmuonfilter_decision, &b_inconsistentPFmuonfilter_decision);
   fChain->SetBranchAddress("hcallaserfilter_decision", &hcallaserfilter_decision, &b_hcallaserfilter_decision);
   fChain->SetBranchAddress("ecallaserfilter_decision", &ecallaserfilter_decision, &b_ecallaserfilter_decision);
   fChain->SetBranchAddress("eenoisefilter_decision", &eenoisefilter_decision, &b_eenoisefilter_decision);
   fChain->SetBranchAddress("eebadscfilter_decision", &eebadscfilter_decision, &b_eebadscfilter_decision);
   fChain->SetBranchAddress("trackercoherentnoisefilter1_decision", &trackercoherentnoisefilter1_decision, &b_trackercoherentnoisefilter1 );
   fChain->SetBranchAddress("trackercoherentnoisefilter2_decision", &trackercoherentnoisefilter2_decision, &b_trackercoherentnoisefilter2 );
   fChain->SetBranchAddress("trackertoomanyclustersfilter_decision", &trackertoomanyclustersfilter_decision, &b_trackertoomanyclustersfilter);
   fChain->SetBranchAddress("trackertoomanytripletsfilter_decision", &trackertoomanytripletsfilter_decision, &b_trackertoomanytripletsfilter);
   fChain->SetBranchAddress("trackertoomanyseedsfilter_decision", &trackertoomanyseedsfilter_decision, &b_trackertoomanyseedsfilter );
   fChain->SetBranchAddress("passprescalePFHT350filter_decision", &passprescalePFHT350filter_decision, &b_passprescalePFHT350filter_decision);
   fChain->SetBranchAddress("passprescaleHT250filter_decision", &passprescaleHT250filter_decision, &b_passprescaleHT250filter_decision);
   fChain->SetBranchAddress("passprescaleHT300filter_decision", &passprescaleHT300filter_decision, &b_passprescaleHT300filter_decision);
   fChain->SetBranchAddress("passprescaleHT350filter_decision", &passprescaleHT350filter_decision, &b_passprescaleHT350filter_decision);
   fChain->SetBranchAddress("passprescaleHT400filter_decision", &passprescaleHT400filter_decision, &b_passprescaleHT400filter_decision);
   fChain->SetBranchAddress("passprescaleHT450filter_decision", &passprescaleHT450filter_decision, &b_passprescaleHT450filter_decision);
   fChain->SetBranchAddress("passprescaleJet30MET80filter_decision", &passprescaleJet30MET80filter_decision, &b_passprescaleJet30MET80filter_decision);
   fChain->SetBranchAddress("MPT", &MPT, &b_MPT);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL2L3", &jets_AK5PFclean_corrL2L3, &b_jets_AK5PFclean_corrL2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL2L3Residual", &jets_AK5PFclean_corrL2L3Residual, &b_jets_AK5PFclean_corrL2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1FastL2L3", &jets_AK5PFclean_corrL1FastL2L3, &b_jets_AK5PFclean_corrL1FastL2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1L2L3", &jets_AK5PFclean_corrL1L2L3, &b_jets_AK5PFclean_corrL1L2L3);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1FastL2L3Residual", &jets_AK5PFclean_corrL1FastL2L3Residual, &b_jets_AK5PFclean_corrL1FastL2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_corrL1L2L3Residual", &jets_AK5PFclean_corrL1L2L3Residual, &b_jets_AK5PFclean_corrL1L2L3Residual);
   fChain->SetBranchAddress("jets_AK5PFclean_Uncert", &jets_AK5PFclean_Uncert, &b_jets_AK5PFclean_Uncert);
   fChain->SetBranchAddress("PU_zpositions", &PU_zpositions, &b_PU_zpositions);
   fChain->SetBranchAddress("PU_sumpT_lowpT", &PU_sumpT_lowpT, &b_PU_sumpT_lowpT);
   fChain->SetBranchAddress("PU_sumpT_highpT", &PU_sumpT_highpT, &b_PU_sumpT_highpT);
   fChain->SetBranchAddress("PU_ntrks_lowpT", &PU_ntrks_lowpT, &b_PU_ntrks_lowpT);
   fChain->SetBranchAddress("PU_ntrks_highpT", &PU_ntrks_highpT, &b_PU_ntrks_highpT);
   fChain->SetBranchAddress("PU_NumInteractions", &PU_NumInteractions, &b_PU_NumInteractions);
   fChain->SetBranchAddress("PU_bunchCrossing", &PU_bunchCrossing, &b_PU_bunchCrossing);
   fChain->SetBranchAddress("PU_TrueNumInteractions", &PU_TrueNumInteractions, &b_PU_TrueNumInteractions);
   fChain->SetBranchAddress("rho_kt6PFJetsForIsolation2011", &rho_kt6PFJetsForIsolation2011, &b_rho_kt6PFJetsForIsolation2011);
   fChain->SetBranchAddress("rho_kt6PFJetsForIsolation2012", &rho_kt6PFJetsForIsolation2012, &b_rho_kt6PFJetsForIsolation2012);
   fChain->SetBranchAddress("pfmets_fullSignif", &pfmets_fullSignif, &b_pfmets_fullSignif);
   fChain->SetBranchAddress("pfmets_fullSignifCov00", &pfmets_fullSignifCov00, &b_pfmets_fullSignifCov00);
   fChain->SetBranchAddress("pfmets_fullSignifCov10", &pfmets_fullSignifCov10, &b_pfmets_fullSignifCov10);
   fChain->SetBranchAddress("pfmets_fullSignifCov11", &pfmets_fullSignifCov11, &b_pfmets_fullSignifCov11);
   fChain->SetBranchAddress("softjetUp_dMEx", &softjetUp_dMEx, &b_softjetUp_dMEx);
   fChain->SetBranchAddress("softjetUp_dMEy", &softjetUp_dMEy, &b_softjetUp_dMEy);
   fChain->SetBranchAddress("pdfweights_cteq", &pdfweights_cteq, &b_pdfweights_cteq);
   fChain->SetBranchAddress("pdfweights_mstw", &pdfweights_mstw, &b_pdfweights_mstw);
   fChain->SetBranchAddress("pdfweights_nnpdf", &pdfweights_nnpdf, &b_pdfweights_nnpdf);
   fChain->SetBranchAddress("photon_chIsoValues", &photon_chIsoValues, &b_photon_chIsoValues);
   fChain->SetBranchAddress("photon_phIsoValues", &photon_phIsoValues, &b_photon_phIsoValues);
   fChain->SetBranchAddress("photon_nhIsoValues", &photon_nhIsoValues, &b_photon_nhIsoValues);
   fChain->SetBranchAddress("photon_passElectronVeto", &photon_passElectronVeto, &b_photon_passElectronVeto);
   fChain->SetBranchAddress("puJet_rejectionBeta", &puJet_rejectionBeta, &b_puJet_rejectionBeta);
   fChain->SetBranchAddress("puJet_rejectionMVA", &puJet_rejectionMVA, &b_puJet_rejectionMVA);
   fChain->SetBranchAddress("pfmets_fullSignif_2012", &pfmets_fullSignif_2012, &b_pfmets_fullSignif_2012);
   fChain->SetBranchAddress("pfmets_fullSignifCov00_2012", &pfmets_fullSignifCov00_2012, &b_pfmets_fullSignifCov00_2012);
   fChain->SetBranchAddress("pfmets_fullSignifCov10_2012", &pfmets_fullSignifCov10_2012, &b_pfmets_fullSignifCov10_2012);
   fChain->SetBranchAddress("pfmets_fullSignifCov11_2012", &pfmets_fullSignifCov11_2012, &b_pfmets_fullSignifCov11_2012);
   fChain->SetBranchAddress("pfmets_fullSignif_2012_dataRes", &pfmets_fullSignif_2012_dataRes, &b_pfmets_fullSignif_2012_dataRes);
   fChain->SetBranchAddress("pfmets_fullSignifCov00_2012_dataRes", &pfmets_fullSignifCov00_2012_dataRes, &b_pfmets_fullSignifCov00_2012_dataRes);
   fChain->SetBranchAddress("pfmets_fullSignifCov10_2012_dataRes", &pfmets_fullSignifCov10_2012_dataRes, &b_pfmets_fullSignifCov10_2012_dataRes);
   fChain->SetBranchAddress("pfmets_fullSignifCov11_2012_dataRes", &pfmets_fullSignifCov11_2012_dataRes, &b_pfmets_fullSignifCov11_2012_dataRes);
   fChain->SetBranchAddress("isotk_pt", &isotk_pt, &b_isotk_pt);
   fChain->SetBranchAddress("isotk_phi", &isotk_phi, &b_isotk_phi);
   fChain->SetBranchAddress("isotk_eta", &isotk_eta, &b_isotk_eta);
   fChain->SetBranchAddress("isotk_iso", &isotk_iso, &b_isotk_iso);
   fChain->SetBranchAddress("isotk_dzpv", &isotk_dzpv, &b_isotk_dzpv);
   fChain->SetBranchAddress("isotk_charge", &isotk_charge, &b_isotk_charge);
   
}
