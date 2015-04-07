#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <iostream> 

#define MinSignalLeptonPt 20
#define MinVetoLeptonPt 10 /*15*/
#define MinJetPt 40 /*40*/ 

using namespace std;

//
// dPhi, dR, ... 
//
float getDPhi(float phi1, float phi2)
{
    float absdphi = abs(phi1-phi2);
    if(absdphi < TMath::Pi()) return absdphi;
    else return (2*TMath::Pi() - absdphi);
}
float getDR(float dphi, float deta)
{
    return TMath::Sqrt(dphi*dphi+deta*deta);
}
float getDR(float eta1, float eta2, float phi1, float phi2)
{
    return getDR(getDPhi(phi1, phi2), eta1-eta2);
}


double getDZ(double vx, double vy, double vz, double px, double py, double pz, int firstGoodVertex)
{
    return vz - pv_z->at(firstGoodVertex) 
        -((vx-pv_x->at(firstGoodVertex))*px+(vy-pv_y->at(firstGoodVertex))*py)*pz/(px*px+py*py); 
}

//
// Code for mini-isolation
// Taken from Jack's code : https://github.com/jbradmil/csa14/blob/master/src/event_handler.cpp#L3028
// 
double GetIsolation(const int ilep, const int ParticleType, const double rmax, const bool mini, const bool addCH, const bool addPH, const bool addNH, const bool usePFweight) 
{ 

    double ptThresh(0.5);
    double lep_pt(0.), lep_eta(0.), lep_phi(0.), deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);;
    if(ParticleType==11) {
        lep_pt = els_pt->at(ilep);
        lep_eta = els_eta->at(ilep);
        lep_phi = els_phi->at(ilep);
        ptThresh = 0;
        if (fabs(lep_eta)>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    }else if(ParticleType==13){
        lep_pt = mus_pt->at(ilep);
        lep_eta = mus_eta->at(ilep);
        lep_phi = mus_phi->at(ilep);
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else{
        lep_pt = pfcand_pt->at(ilep);
        lep_eta = pfcand_eta->at(ilep);
        lep_phi = pfcand_phi->at(ilep);
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // Using muon cones
    }

    bool need_pfweight = false;
    if(usePFweight) need_pfweight = true;

    double riso_max = std::max(0.4,10./lep_pt);
    // find the PF cands that matches the lepton
    double drmin = DBL_MAX;
    uint match_index = 9999999;
    for (unsigned int icand = 0; icand < pfcand_pt->size(); icand++) {
        if(isnan(pfcand_eta->at(icand))
                || isnan(pfcand_phi->at(icand))) continue;
        double dr = getDR(pfcand_eta->at(icand), lep_eta, pfcand_phi->at(icand), lep_phi);
        if (dr < drmin){
            drmin = dr;
            match_index = icand;
        }
    }

    // 11, 13, 22 for ele/mu/gamma, 211 for charged hadrons, 130 for neutral hadrons,
    // 1 and 2 for hadronic and em particles in HF
    double iso_nh(0.), iso_ph(0.), iso_ch(0.), iso_pu(0.);
    // get cone size--shrinking if mini, otherwise fixed
    double R(0);
    if (mini) R=std::max(0.05,std::min(rmax, 10./lep_pt));
    else R=rmax;

    for (unsigned int icand = 0; icand < pfcand_pt->size(); icand++) {
        if (icand==match_index) continue;
        uint pdgId = TMath::Nint(pfcand_pdgId->at(icand));
        if (pdgId<7) continue;
        if(isnan(pfcand_pt->at(icand))
                || isnan(pfcand_eta->at(icand))
                || isnan(pfcand_phi->at(icand))) continue;
        double dr = getDR(pfcand_eta->at(icand), lep_eta, pfcand_phi->at(icand), lep_phi);
        if (dr > riso_max) continue;
        ////////////////// NEUTRALS /////////////////////////
        if (pfcand_charge->at(icand)==0){
            if (pfcand_pt->at(icand)>ptThresh) {
                double wpv(0.), wpu(0.), wpf(1.);
                for (unsigned int jcand = 0; need_pfweight && jcand < pfcand_pt->size(); jcand++) {
                    if (pfcand_charge->at(icand)!=0 || icand==jcand) continue;
                    double jpt = pfcand_pt->at(jcand);
                    double jdr = getDR(pfcand_eta->at(icand), pfcand_eta->at(jcand),
                            pfcand_phi->at(icand), pfcand_phi->at(jcand));
                    if(jdr<=0) continue; // We can either not count it, or count it infinitely...
                    if (pfcand_fromPV->at(icand)>1) wpv += log(jpt/jdr);
                    else wpu += log(jpt/jdr);
                }
                /////////// PHOTONS ////////////
                if (pdgId==22) {
                    if(dr < deadcone_ph) continue;
                    wpf = (usePFweight)?(wpv/(wpv+wpu)):1.;
                    if (dr<R) iso_ph += wpf*pfcand_pt->at(icand);
                }
                /////////// NEUTRAL HADRONS ////////////
                else if (pdgId==130) {
                    if(dr < deadcone_nh) continue;
                    wpf = (usePFweight)?(wpv/(wpv+wpu)):1.;
                    if (dr<R) iso_nh += wpf*pfcand_pt->at(icand);
                }
            }
            ////////////////// CHARGED from PV /////////////////////////
        } else if (pfcand_fromPV->at(icand)>1){
            if (fabs(pdgId)==211) {
                if(dr < deadcone_ch) continue;
                if (dr<R) {
                    // if (ParticleType==13) cout << "Adding to cone.." << endl;
                    // if (ParticleType==13) printf("imu %d, ch pfcand %d: pt=%f, deadcone=%f, dr cut=%f, dr=%f\n",ilep, icand, pfcand_pt->at(icand), deadcone_ch, R, dr); 
                    iso_ch += pfcand_pt->at(icand);
                }
            }
            ////////////////// CHARGED from PU /////////////////////////
        } else {
            if (pfcand_pt->at(icand)>ptThresh){
                if(dr < deadcone_pu) continue;
                if (dr<R) iso_pu += pfcand_pt->at(icand);}
        }
    }

    // now add the components
    double isolation(0.);
    if(addPH) isolation += iso_ph;
    if(addNH) isolation += iso_nh;
    if(addPH && addNH && !usePFweight) isolation -= iso_pu/2.;
    if(isolation < 0) isolation = 0;
    if(addCH) isolation += iso_ch;

    // if (ParticleType==13){
    //   printf("Muon %d: pt=%3.3f, mini_iso=%3.3f, iso_ph=%3.3f, iso_nh=%3.3f, iso_ch=%3.3f, iso_pu=%3.3f\n", ilep, mus_pt->at(ilep), isolation, iso_ph, iso_nh, iso_ch, iso_pu); 
    // }

    return isolation/lep_pt;
}


/////////////////////////////////////////////////////////////////////////
////////////////////////////////  MUONS  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// Ref : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon
// The veto muon Id is the same as the tight but with lower pT cut. 
// The signal muon isolation is R= 0.4, iso < 0.12 
// The veto muon isolation is R= 0.4, iso < 0.2
float GetMuonIsolation(int imu)
{
    if(imu >= (int)mus_pt->size()) return -999;
    //double sumEt = mus_pfIsolationR03_sumNeutralHadronEt->at(imu) + mus_pfIsolationR03_sumPhotonEt->at(imu) 
    //    - 0.5*mus_pfIsolationR03_sumPUPt->at(imu);
    double sumEt = mus_pfIsolationR04_sumNeutralHadronEt->at(imu) + mus_pfIsolationR04_sumPhotonEt->at(imu) 
        - 0.5*mus_pfIsolationR04_sumPUPt->at(imu);
    if(sumEt<0.0) sumEt=0.0;
    //return (mus_pfIsolationR03_sumChargedHadronPt->at(imu) + sumEt)/mus_pt->at(imu);
    return (mus_pfIsolationR04_sumChargedHadronPt->at(imu) + sumEt)/mus_pt->at(imu);
}

bool IsBasicMuon(int imu)
{
    if(imu >= (int)mus_pt->size()) return false;

    float d0PV = mus_tk_d0dum->at(imu)-pv_x->at(0)*sin(mus_tk_phi->at(imu))+pv_y->at(0)*cos(mus_tk_phi->at(imu));

    return ( //mus_isGlobalMuon->at(imu) > 0
            mus_isPF->at(imu)
            // recoMu.globalTrack()->normalizedChi2() < 10.
            // recoMu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
            && mus_id_GlobalMuonPromptTight->at(imu)> 0  // this includes above two and isGlobalMuon 
            && mus_numberOfMatchedStations->at(imu) > 1
            //&& fabs(mus_dB->at(imu)) < 0.02
            && fabs(d0PV) < 0.2/*0.02*/
            //&& fabs(getDZ(mus_tk_vx->at(imu), mus_tk_vy->at(imu), mus_tk_vz->at(imu), mus_tk_px->at(imu), 
            //        mus_tk_py->at(imu), mus_tk_pz->at(imu), 0)) < 0.5
            && fabs(mus_tk_vz->at(imu)-pv_z->at(0))<0.5
            && mus_tk_numvalPixelhits->at(imu) > 0
            && mus_tk_LayersWithMeasurement->at(imu) > 5
            
            && mus_pt->at(imu) >= MinSignalLeptonPt
            && fabs(mus_eta->at(imu)) <= 2.4);
}

bool IsSignalMuon(int imu, bool doMiniIso)
{
    if(imu >= (int)mus_pt->size()) return false;

    bool passIso=false; 
    if(!doMiniIso)  passIso = GetMuonIsolation(imu)<0.12;  
    if(doMiniIso)   passIso = GetIsolation(imu, 13, 0.2, false, true, true, true, false)<0.2;  
    
    return (IsBasicMuon(imu) && passIso); 
}

bool IsVetoMuon(int imu, bool doMiniIso)
{
    if(imu >= (int)mus_pt->size()) return false;
    if(IsSignalMuon(imu, doMiniIso)) return false; // not signal muon

    float d0PV = mus_tk_d0dum->at(imu)-pv_x->at(0)*sin(mus_tk_phi->at(imu))+pv_y->at(0)*cos(mus_tk_phi->at(imu));
    
    bool passIso=false; 
    if(!doMiniIso)  passIso = GetMuonIsolation(imu)<0.2;  
    if(doMiniIso)   passIso = GetIsolation(imu, 13, 0.2, false, true, true, true, false)<0.2;  // FIXME : make sure the cut value here  
/*
    return ((mus_isGlobalMuon->at(imu) >0 || mus_isTrackerMuon->at(imu) >0)
            && mus_isPF->at(imu) 
            && fabs(getDZ(mus_tk_vx->at(imu), mus_tk_vy->at(imu), mus_tk_vz->at(imu), mus_tk_px->at(imu), 
                    mus_tk_py->at(imu), mus_tk_pz->at(imu), 0)) < 0.5 
            && mus_pt->at(imu) >= MinVetoLeptonPt
            && fabs(mus_eta->at(imu)) <= 2.5
*/
    return (//mus_isGlobalMuon->at(imu) > 0
            mus_isPF->at(imu)
            // recoMu.globalTrack()->normalizedChi2() < 10.
            // recoMu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
            && mus_id_GlobalMuonPromptTight->at(imu) > 0   // this includes above two and isGlobalMuon
            && mus_numberOfMatchedStations->at(imu) > 1
            //&& fabs(mus_dB->at(imu)) < 0.02
            && fabs(d0PV) < 0.2/*0.02*/
            //&& fabs(getDZ(mus_tk_vx->at(imu), mus_tk_vy->at(imu), mus_tk_vz->at(imu), mus_tk_px->at(imu), 
            //        mus_tk_py->at(imu), mus_tk_pz->at(imu), 0)) < 0.5
            && fabs(mus_tk_vz->at(imu)-pv_z->at(0))<0.5
            && mus_tk_numvalPixelhits->at(imu) > 0
            && mus_tk_LayersWithMeasurement->at(imu) > 5
            
            && mus_pt->at(imu) >= MinVetoLeptonPt
            && fabs(mus_eta->at(imu)) <= 2.4
            && passIso);
}


vector<int> GetMuons(bool doSignal, bool doMiniIso)
{
    vector<int> muons;
    for(int index=0; index<(int)mus_pt->size(); index++)
        if(doSignal)
        {
            if(IsSignalMuon(index, doMiniIso)) muons.push_back(index);
        } else 
        {
            if(IsVetoMuon(index, doMiniIso)) muons.push_back(index);
        }
    return muons;
}


/////////////////////////////////////////////////////////////////////////
//////////////////////////////  ELECTRONS  //////////////////////////////
/////////////////////////////////////////////////////////////////////////
// Ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#CSA14_selection_conditions_25ns

float GetElectronIsolation(int iel)
{
    float absiso = els_pfIsolationR03_sumChargedHadronPt->at(iel) + std::max(0.0 , els_pfIsolationR03_sumNeutralHadronEt->at(iel) + els_pfIsolationR03_sumPhotonEt->at(iel) - 0.5 * els_pfIsolationR03_sumPUPt->at(iel) );
    return absiso/els_pt->at(iel);
}

bool IsBasicElectron(int iel) // Medium working point 
{
    if(iel >= (int)els_pt->size()) return false;

    float d0PV = els_d0dum->at(iel)-pv_x->at(0)*sin(els_tk_phi->at(iel))+pv_y->at(0)*cos(els_tk_phi->at(iel));
    
    return (els_pt->at(iel) > MinSignalLeptonPt
            && fabs(els_scEta->at(iel)) < 2.5
            && ((els_isEB->at(iel) // Barrel selection
                 && fabs(els_dEtaIn->at(iel)) < 0.0106/*0.004*/
                 && fabs(els_dPhiIn->at(iel)) < 0.0323/*0.06*/
                 && els_full5x5_sigmaIetaIeta->at(iel) < 0.0107/*0.01*/
                 && els_hadOverEm->at(iel) < 0.067/*0.12*/
                 && fabs(d0PV) < 0.0131/*0.02*/
                 //&& fabs(getDZ(els_vx->at(iel), els_vy->at(iel), els_vz->at(iel), cos(els_tk_phi->at(iel))*els_tk_pt->at(iel), 
                 //        sin(els_tk_phi->at(iel))*els_tk_pt->at(iel), els_tk_pz->at(iel), 0)) < 0.2231/*0.1*/
                 && fabs(els_vz->at(iel) - pv_z->at(0))<0.2231
                 && fabs(1./els_caloEnergy->at(iel) - els_eOverPIn->at(iel)/els_caloEnergy->at(iel)) < 0.1043/*0.05*/
                && els_PATpassConversionVeto->at(iel)             // was !els_hasMatchedConversion->at(iel)
                && els_expectedMissingInnerHits->at(iel) <= 1     // was els_n_inner_layer->at(iel) <= 1 
                ) ||
                (els_isEE->at(iel)  // Endcap selection
                 && fabs(els_dEtaIn->at(iel)) < 0.0108/* 0.007*/
                 && fabs(els_dPhiIn->at(iel)) < 0.0455/*0.03*/
                 && els_full5x5_sigmaIetaIeta->at(iel) < 0.0318/*0.03*/
                 && els_hadOverEm->at(iel) < 0.097/*0.10*/ 
                 && fabs(d0PV) < 0.0845/*0.02*/ 
                 //&& fabs(getDZ(els_vx->at(iel), els_vy->at(iel), els_vz->at(iel), cos(els_tk_phi->at(iel))*els_tk_pt->at(iel), 
                 //        sin(els_tk_phi->at(iel))*els_tk_pt->at(iel), els_tk_pz->at(iel), 0)) < 0.7523/*0.1*/
                 && fabs(els_vz->at(iel) - pv_z->at(0))<0.7523
                 && fabs(1./els_caloEnergy->at(iel) - els_eOverPIn->at(iel)/els_caloEnergy->at(iel)) < 0.1201/*0.05*/ 
                 && els_PATpassConversionVeto->at(iel)            // was !els_hasMatchedConversion->at(iel)
                 && els_expectedMissingInnerHits->at(iel) <= 1    // was els_n_inner_layer->at(iel) <= 1 
                ))
           );
}

bool IsSignalElectron(int iel, bool doMiniIso)
{
    if(iel >= (int)els_pt->size()) return false;
   
    bool passIso=false;
    float isocut=0.2179;  // Medium working point 
    if(els_isEE->at(iel)) isocut=0.254;
    if(!doMiniIso)  passIso=GetElectronIsolation(iel)<isocut;
    if(doMiniIso)   passIso=GetIsolation(iel, 11, 0.2, false, true, true, true, false)<0.1;
    
    return (IsBasicElectron(iel) && passIso);
}

bool IsVetoElectron(int iel, bool doMiniIso)
{
    if(iel >= (int)els_pt->size()) return false;
    if(IsSignalElectron(iel, doMiniIso)) return false; // not signal electron 

    float d0PV = els_d0dum->at(iel)-pv_x->at(0)*sin(els_tk_phi->at(iel))+pv_y->at(0)*cos(els_tk_phi->at(iel));
    // 
    bool passIso=false;
    float isocut=0.3313;   
    if(els_isEE->at(iel)) isocut=0.3816;
    if(!doMiniIso)  passIso=GetElectronIsolation(iel)<isocut;
    if(doMiniIso)   passIso=GetIsolation(iel, 11, 0.2, false, true, true, true, false)<0.1;

    return (els_pt->at(iel) > MinVetoLeptonPt
            && fabs(els_scEta->at(iel)) < 2.5
            && ((els_isEB->at(iel) // Endcap selection
                    && fabs(els_dEtaIn->at(iel)) < 0.02/*0.007*/
                    && fabs(els_dPhiIn->at(iel)) < 0.2579/*0.8*/
                    && els_full5x5_sigmaIetaIeta->at(iel) < 0.0125/*0.01*/
                    && els_hadOverEm->at(iel) < 0.2564/*0.15*/
                    && fabs(d0PV) < 0.025/*0.04*/ 
                    //&& fabs(getDZ(els_vx->at(iel), els_vy->at(iel), els_vz->at(iel), cos(els_tk_phi->at(iel))*els_tk_pt->at(iel), 
                    //              sin(els_tk_phi->at(iel))*els_tk_pt->at(iel), els_tk_pz->at(iel), 0)) < 0.5863/*0.2*/
                    && fabs(els_vz->at(iel) - pv_z->at(0))<0.5863
                    && fabs(1./els_caloEnergy->at(iel) - els_eOverPIn->at(iel)/els_caloEnergy->at(iel)) < 0.1508 
                    && els_PATpassConversionVeto->at(iel)             // was !els_hasMatchedConversion->at(iel)
                    && els_expectedMissingInnerHits->at(iel) <= 2     // was els_n_inner_layer->at(iel) <= 2  
                    && passIso
                    ) ||
                (els_isEE->at(iel)  // Barrel selection
                    && fabs(els_dEtaIn->at(iel)) < 0.0141/*0.01*/
                    && fabs(els_dPhiIn->at(iel)) < 0.2591/*0.7*/
                    && els_full5x5_sigmaIetaIeta->at(iel) < 0.0371/*0.03*/
                    && els_hadOverEm->at(iel) < 0.1335/*0.15*/
                    && fabs(d0PV) < 0.2232/*0.04*/ 
                    //&& fabs(getDZ(els_vx->at(iel), els_vy->at(iel), els_vz->at(iel), cos(els_tk_phi->at(iel))*els_tk_pt->at(iel), 
                    //           sin(els_tk_phi->at(iel))*els_tk_pt->at(iel), els_tk_pz->at(iel), 0)) < 0.9513/*0.2*/
                    && fabs(els_vz->at(iel) - pv_z->at(0))<0.9513
                    && fabs(1./els_caloEnergy->at(iel) - els_eOverPIn->at(iel)/els_caloEnergy->at(iel)) < 0.1542 
                    && els_PATpassConversionVeto->at(iel)             // was !els_hasMatchedConversion->at(iel)
                    && els_expectedMissingInnerHits->at(iel) <= 3     // was els_n_inner_layer->at(iel) <= 3 
                    && passIso
                 ))
           );  
}

vector<int> GetElectrons(bool doSignal, bool doMiniIso)
{
    vector<int> electrons;
    for(int index=0; index<(int)els_pt->size(); index++)
        if(doSignal)
        {
            if(IsSignalElectron(index, doMiniIso)) electrons.push_back(index);
        }   
        else 
        {
            if(IsVetoElectron(index, doMiniIso)) electrons.push_back(index);
        }
    return electrons;
}


/////////////////////////////////////////////////////////////////////////
////////////////////////////  Track Veto  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////
// From Jack's code : https://github.com/jbradmil/csa14/blob/master/src/event_handler.cpp
bool PassIsoTrackBaseline(const uint itk)
{
    if (static_cast<int>(pfcand_charge->at(itk))==0) return false;
    if (pfcand_pt->at(itk)<3) return false;
    if (fabs(pfcand_eta->at(itk))>2.5) return false;
    if (pfcand_fromPV->at(itk)<=1) return false; // pileup suppression
    if (static_cast<int>(pfcand_pdgId->at(itk))==0) return false;
    return true;
}

double GetPFCandIsolation(const uint indexA)
{ // absolute, not relative -- charged tracks only
    double isoSum(0);
    for (uint other(0); other<pfcand_pdgId->size(); other++) {
        if (isnan(pfcand_charge->at(other)) || isnan(pfcand_pt->at(other))  || isnan(pfcand_fromPV->at(other)) || isnan(pfcand_phi->at(other)) || isnan(pfcand_eta->at(other)) ) continue;
        if (other==indexA) continue; // don't count track in it's own isolation sum
        if (static_cast<int>(pfcand_charge->at(other))==0) continue; // only consider charged tracks
        if (fabs(pfcand_fromPV->at(other))<=1) continue; // pileup suppression
        double deltaR = getDR(pfcand_eta->at(indexA), pfcand_eta->at(other), 
                              pfcand_phi->at(indexA), pfcand_phi->at(other)); // isolation cone
        if (deltaR>0.3) continue;
        isoSum+=pfcand_pt->at(other);
    }
    return isoSum;
}

void GetIsoTracks(std::vector<std::pair<int,double> > &eCands, 
                  std::vector<std::pair<int,double> > &muCands, 
                  std::vector<std::pair<int,double> > &hadCands) 
{
    eCands.clear();
    muCands.clear();
    hadCands.clear();
    // cout << "Found " << pfcand_pt->size() << " PFCands." << endl;
    for (uint itk(0); itk<pfcand_pdgId->size(); itk++) {
        if (isnan(pfcand_charge->at(itk)) || isnan(pfcand_pt->at(itk))  || isnan(pfcand_dz->at(itk)) || isnan(pfcand_phi->at(itk)) || isnan(pfcand_eta->at(itk)) ) continue;
        if (!PassIsoTrackBaseline(itk)) continue;
        double pt = pfcand_pt->at(itk);
        if (pt<5) continue;
        //  int type = static_cast<int>(pfcand_pdgId->at(itk));
        double iso = GetPFCandIsolation(itk);
        double relIso=iso/pfcand_pt->at(itk);
        // note: not cutting here on isolation!
        switch (abs(TMath::Nint(pfcand_pdgId->at(itk)))) {
            case 11:
                eCands.push_back(std::make_pair(itk, relIso));
                break;
            case 13:
                muCands.push_back(std::make_pair(itk, relIso));
                break;
            case 211:
                // if (pfcand_pt->at(itk)>10) printf("pfcand %d: pdgId=%d, pt=%f, ch_rel_iso=%f, mT=%f\n", itk, TMath::Nint(pfcand_pdgId->at(itk)), pfcand_pt->at(itk), relIso, GetMTW(pfcand_pt->at(itk),pfTypeImets_et->at(0),pfcand_phi->at(itk),pfTypeImets_phi->at(0)));
                if (pt>10) hadCands.push_back(std::make_pair(itk, relIso));
                break;
            default:
                continue;
        }
    }
    return;
}
/*
int NIsoTkVeto()
{
    int isotks = 0;
    for(size_t itk = 0; itk < tks_pt().size() && isotks < 2; ++itk){
        switch(abs(tks_id().at(itk))){
            case 11:
            case 13:
                if(tks_pt().at(itk)>5.
                        && tks_r03_ch().at(itk)<0.2) ++isotks;
                break;
            default:
                if(tks_pt().at(itk)>10
                        && tks_r03_ch().at(itk)<0.1) ++isotks;
                break;
        }
    }
    return isotks;
}
*/

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////  JETS  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////
bool IsBasicJet( unsigned int ijet) 
{
    double rawRatio =(jets_AK4_rawPt->at(ijet)/jets_AK4_pt->at(ijet)); // Same as jets_AK4_corrFactorRaw
    double jetenergy = jets_AK4_energy->at(ijet) * rawRatio;
    double NEF = -999., CEF = -999., NHF=-999., CHF=-999.;
    double chgMult=jets_AK4_chg_Mult->at(ijet);
    double numConst=jets_AK4_mu_Mult->at(ijet)+jets_AK4_neutral_Mult->at(ijet)+jets_AK4_chg_Mult->at(ijet);

    if(jetenergy > 0)
    {
        NEF = jets_AK4_neutralEmE->at(ijet)/jetenergy;
        CEF = jets_AK4_chgEmE->at(ijet)/jetenergy;
        NHF = jets_AK4_neutralHadE->at(ijet)/jetenergy;
        CHF = jets_AK4_chgHadE->at(ijet)/jetenergy;   
    }

    return (NEF < 0.99 && CEF < 0.99 && NHF < 0.99 && CHF > 0 &&
            chgMult > 0 && numConst > 1);
}

bool IsGoodJet( unsigned int ijet,  double ptThresh,  double etaThresh) 
{
    if(jets_AK4_pt->size()<=ijet) return false;
    if(!IsBasicJet(ijet)) return false;
    if(jets_AK4_pt->at(ijet)<ptThresh || fabs(jets_AK4_eta->at(ijet))>etaThresh) return false;
    return true;
}

vector<int> GetJets(vector<int> SigEl, vector<int> SigMu, vector<int> VetoEl, vector<int> VetoMu, float &HT)
{
    vector<int> jets;
    vector<bool> jet_is_lepton(jets_AK4_pt->size(), false);
    HT = SigEl.size()+VetoEl.size()+SigMu.size()+VetoMu.size(); // To avoid warnings
    HT = 0;
    // Finding jets that contain good leptons
    for(int index = 0; index < (int)SigEl.size(); index++) 
    {
        int ijet = els_jet_ind->at(SigEl[index]);
        if(ijet >= 0) 
        {
            jet_is_lepton[ijet] = true;
        }
    }
    for(int index = 0; index < (int)VetoEl.size(); index++) 
    {
        int ijet = els_jet_ind->at(VetoEl[index]);
        if(ijet >= 0) jet_is_lepton[ijet] = true;
    }

    for(int index = 0; index < (int)SigMu.size(); index++) 
    {
        int ijet = mus_jet_ind->at(SigMu[index]);
        if(ijet >= 0) 
        {
            jet_is_lepton[ijet] = true;
        }
    }
    for(int index = 0; index < (int)VetoMu.size(); index++) 
    {
        int ijet = mus_jet_ind->at(VetoMu[index]);
        if(ijet >= 0) jet_is_lepton[ijet] = true;
    }

    // Tau/photon cleaning, and calculation of HT
    for(int ijet = 0; ijet<(int)jets_AK4_pt->size(); ijet++) 
    {
        if(!IsGoodJet(ijet, MinJetPt, 2.4) || jet_is_lepton[ijet]) continue;

        // double tmpdR, partp, jetp = sqrt(pow(jets_AK4_px->at(ijet),2)+pow(jets_AK4_py->at(ijet),2)+pow(jets_AK4_pz->at(ijet),2));
        // bool useJet = true;
        // Tau cleaning: jet rejected if withing deltaR = 0.4 of tau, and momentum at least 60% from tau
        // for(int index = 0; index < taus_pt->size(); index++) 
        //{
        //   tmpdR = dR(jets_AK4_eta->at(ijet), taus_eta->at(index), jets_AK4_phi->at(ijet), taus_phi->at(index)); 
        //   partp = sqrt(pow(taus_px->at(index),2)+pow(taus_py->at(index),2)+pow(taus_pz->at(index),2));
        //   if(tmpdR < 0.4 && partp/jetp >= 0.6)
        //{
        //    useJet = false; break;}
        // }
        // if(!useJet) continue;

        // // Photon cleaning: jet rejected if withing deltaR = 0.4 of photon, and momentum at least 60% from photon
        // for(int index = 0; index < photons_pt->size(); index++) {
        //   tmpdR = dR(jets_AK4_eta->at(ijet), photons_eta->at(index), jets_AK4_phi->at(ijet), photons_phi->at(index));   
        //   partp = sqrt(pow(photons_px->at(index),2)+pow(photons_py->at(index),2)+pow(photons_pz->at(index),2));
        //   if(tmpdR < 0.4 && partp/jetp >= 0.6){useJet = false; break;}
        // }
        // if(!useJet) continue;

        if(jets_AK4_pt->at(ijet) > MinJetPt) 
        {
            HT += jets_AK4_pt->at(ijet);
            jets.push_back(ijet);
        }
    } // Loop over jets
    return jets;
}


/////////////////////////////////////////////////////////////////////////
////////////////////////////  TRUTH-MATCHING  ///////////////////////////
/////////////////////////////////////////////////////////////////////////
int GetTrueParticle(double RecEta, double RecPhi, double &closest_dR)
{
    int closest_imc = -1; 
    double dR = 9999.; closest_dR = 9999.;
    double MCEta, MCPhi;
    for(int imc=0; imc < (int)mc_doc_id->size(); imc++)
    {
        MCEta = mc_doc_eta->at(imc); MCPhi = mc_doc_phi->at(imc);
        dR = sqrt(pow(RecEta-MCEta,2) + pow(RecPhi-MCPhi,2));
        if(dR < closest_dR) 
        {
            closest_dR = dR;
            closest_imc = imc;
        }
    }
    return closest_imc;
}

/////////////////////////////////////////////////////////////////////////
////////////////////////////  EVENT CLEANING  ///////////////////////////
/////////////////////////////////////////////////////////////////////////
bool PassesPVCut()
{
    if(beamSpot_x->size()<1 || pv_x->size()<1) return false;
     double pv_rho(sqrt(pv_x->at(0)*pv_x->at(0) + pv_y->at(0)*pv_y->at(0)));
    if(pv_ndof->at(0)>4 && fabs(pv_z->at(0))<24. && pv_rho<2.0 && pv_isFake->at(0)==0) return true;
    return false;
}

bool PassesMETCleaningCut() 
{
    return cschalofilter_decision
        && hcallaserfilter_decision 
        && trackingfailurefilter_decision 
        && eebadscfilter_decision; 
}

