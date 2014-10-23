/*
Author: Wing To, HEP, UCSB
Org Date: 02-Sep-2010
Helper function for RA4Ana
Original location : /net/cms2/cms2r0/cfA/src/ObjectSelector_Sync.h
*/
#ifndef OBJECTSELECTOR_SYNC_H
#define OBJECTSELECTOR_SYNC_H

#include <iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Rotation3D.h"
#include "Math/EulerAngles.h"
#include "Math/AxisAngle.h"
#include "Math/Quaternion.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/LorentzRotation.h"
#include "Math/Boost.h"
#include "Math/BoostX.h"
#include "Math/BoostY.h"
#include "Math/BoostZ.h"
#include "Math/Transform3D.h"
#include "Math/Plane3D.h"
#include "Math/VectorUtil.h"

#define TIGHTLEPTONPT 20
#define LEPTONVETOPT 15
bool isFunnyJet()  
{
    return false;
    return event==975572;
    //  return false;
    /*   return ((event==41677 && lumiblock==12500598) || (event==41678 && lumiblock==12500611) || (event==41678 && lumiblock==12500809)  */
    /* 	  || (event==41679 && lumiblock==12500986) || (event==41679 && lumiblock==12501146)); */
}

bool isFunnyMuon()
{
    return false;
    return event==975572;
    //  return ((event==41677 && lumiblock==12500598) || (event==41678 && lumiblock==12500611) || (event==41678 && lumiblock==12500809) 
    //	  || (event==41679 && lumiblock==12500986) || (event==41679 && lumiblock==12501146));
}

bool isFunnyElectron()
{
    return isFunnyMuon();
    //  return false;
    //  return event==22307835;
    /*   return ((event==41677 && lumiblock==12500598) || (event==41678 && lumiblock==12500611) || (event==41678 && lumiblock==12500809)  */
    /* 	  || (event==41679 && lumiblock==12500986) || (event==41679 && lumiblock==12501146)); */

}

double GetMax ( double a, double b ) {
    return (b<a)?a:b;
}

int numGoodVertices()
{

    int nvtx = 0;

    for(unsigned int ai = 0;ai<pv_x->size();ai++)
    {      
        if  ( pv_isFake->at(ai)    <    1 && 
                fabs( pv_z->at(ai) ) <=  24 && 
                pv_ndof->at(ai)      >    4 && 
                sqrt( pow(pv_x->at(ai),2)+pow(pv_y->at(ai),2) ) <= 2 ) nvtx++;

    }

    return nvtx;
}

// deltaphi accounting for discontinuity at phi=pi.
double GetDeltaPhi(double phi1, double phi2)
{
    double result = phi1-phi2;
    while (result>TMath::Pi()) result -= 2*TMath::Pi();
    while (result<=-TMath::Pi()) result += 2*TMath::Pi();
    return fabs(result);
}

double GetDeltaR(double eta1,double eta2,double phi1,double phi2){
    return sqrt(pow(GetDeltaPhi(phi1,phi2),2)+pow(eta1-eta2,2));
}

float GetEffectiveArea(float SCEta, bool isMC=false);

    namespace particleId {
    enum leptonType { 
        X=0,
        electron=2,
        muon=3,
    };
}

bool hasPFMatch(int i, particleId::leptonType type, int &pfIdx)
{
    double deltaRVal = 999.;
    double deltaPT = 999.;
    double leptonEta=-999., leptonPhi=-999., leptonPt=-999.;
    if(type == particleId::muon ) {
        leptonEta = mus_eta->at(i);
        leptonPhi = mus_phi->at(i);
        leptonPt = mus_pt->at(i);
    }
    else if(type == particleId::electron) {
        leptonEta = els_eta->at(i);
        leptonPhi = els_phi->at(i);
        leptonPt = els_pt->at(i);
    }

    for(unsigned iCand=0; iCand<pfcand_pt->size(); iCand++) {
        if(pfcand_particleId->at(iCand)==type) {
            double tempDeltaR = GetDeltaR(leptonEta, pfcand_eta->at(iCand), leptonPhi, pfcand_phi->at(iCand));
            if(tempDeltaR < deltaRVal) {
                deltaRVal = tempDeltaR;
                //		      	deltaPT = fabs(leptonPt-pfcand_pt->at(iCand))/leptonPt;
                deltaPT = fabs(leptonPt-pfcand_pt->at(iCand));
                pfIdx=iCand;
            }
        }
    }

    if(isFunnyElectron()) std::cout << "event: " << event << " deltaPT: " << deltaPT << std::endl;

    if(type == particleId::electron) {
        return (deltaPT<10);
    } 
    else {
        return (deltaPT<5);
    }
}

double getDZ(double vx, double vy, double vz, double px, double py, double pz, int firstGoodVertex)
{
    return vz - pv_z->at(firstGoodVertex) -((vx-pv_x->at(firstGoodVertex))*px+(vy-pv_y->at(firstGoodVertex))*py)*pz/(px*px+py*py);   // = dsz/cos(lambda)
    //  return vz - pv_z->at(firstGoodVertex)  -((vx-beamSpot_x->at(0))*px+(vy-beamSpot_y->at(0))*py)*pz/(px*px+py*py);   // = dsz/cos(lambda)
}



vector<int> GetRA4Muon(vector<int> & RA4MuonVeto, std::string cutname="", int firstGoodVertex=0){
    vector<int> RA4Muon;
    bool noPT= (cutname=="PT")?true:false;
    bool noEta= (cutname=="Eta")?true:false;
    bool noIso= (cutname=="Iso")?true:false;
    bool no_isGlobalMuon= (cutname=="isGlobalMuon")?true:false;
    bool no_isPFMuon= (cutname=="isPFMuon")?true:false;
    bool no_GlobalMuonPromptTight= (cutname=="GlobalMuonPromptTight")?true:false;
    bool no_tk_LayersWithMeasurement= (cutname=="tk_LayersWithMeasurement")?true:false;
    bool no_tk_numvalPixelhits= (cutname=="tk_numvalPixelhits")?true:false;
    bool no_numberOfMatchedStations= (cutname=="numberOfMatchedStations")?true:false;
    bool no_d0= (cutname=="D0")?true:false;
    bool no_dz= (cutname=="DZ")?true:false;

    //       std::cout << "\n new event: " << std::endl;

    for(unsigned int i=0;i<mus_pt->size();i++){

        float d0PV, relIso;
        d0PV = mus_tk_d0dum->at(i)-pv_x->at(firstGoodVertex)*sin(mus_tk_phi->at(i))
            +pv_y->at(firstGoodVertex)*cos(mus_tk_phi->at(i));
        double max = mus_pfIsolationR03_sumNeutralHadronEt->at(i) + mus_pfIsolationR03_sumPhotonEt->at(i) - 0.5*mus_pfIsolationR03_sumPUPt->at(i);
        if(max<0.0) max=0.0;
        relIso = (mus_pfIsolationR03_sumChargedHadronPt->at(i) + max)/mus_pt->at(i);

        if(isFunnyMuon() || isFunnyElectron()) {
            std::cout << "Muon:  run " << run << " lumi " << lumiblock << " event " << event
                << " pt " << mus_pt->at(i) << " eta " << mus_eta->at(i) << " phi " << mus_phi->at(i)
                << " isGlobalMuon: " << mus_isGlobalMuon->at(i)
                << " dB: " << mus_dB->at(i)
                << " d0PV: " << d0PV
                << " relIso: " << relIso
                << " sumNeutralHadronEt_R03: " << mus_pfIsolationR03_sumNeutralHadronEt->at(i)
                << " sumPhotonEt_R03: " << mus_pfIsolationR03_sumPhotonEt->at(i)
                << " sumPUPt_R03: " << mus_pfIsolationR03_sumPUPt->at(i)
                << " sumChargedHadronPt_R03: " << mus_pfIsolationR03_sumChargedHadronPt->at(i);

        }
        int pfIdx=-1;
        // 2012 RA4 selection
        if( (mus_isGlobalMuon->at(i) || no_isGlobalMuon)
                && (mus_isPFMuon->at(i) || no_isPFMuon)
                && (mus_id_GlobalMuonPromptTight->at(i)> 0 || no_GlobalMuonPromptTight)
                // included in GlobalMuonPromptTight
                /*  	     && mus_tk_chi2->at(i)/mus_tk_ndof->at(i) < 10 */
                /* 	     && mus_cm_numvalMuonhits->at(i) > 0 */
                && (mus_tk_LayersWithMeasurement->at(i) > 5 || no_tk_LayersWithMeasurement)
                && (mus_tk_numvalPixelhits->at(i) > 0 || no_tk_numvalPixelhits)
                && (mus_numberOfMatchedStations->at(i) > 1 || no_numberOfMatchedStations)
                && (fabs(d0PV) < 0.02 || no_d0)
                //	     && (fabs(mus_dB->at(i)) < 0.02 || no_d0)
                && (fabs(getDZ(mus_tk_vx->at(i), mus_tk_vy->at(i), mus_tk_vz->at(i), mus_tk_px->at(i), mus_tk_py->at(i), mus_tk_pz->at(i), firstGoodVertex)) < 0.5 || no_dz)
                && (mus_pt->at(i) >= TIGHTLEPTONPT || noPT)
                && (fabs(mus_eta->at(i)) <= 2.4 || noEta)
                && (relIso < 0.12 || noIso)
                && hasPFMatch(i, particleId::muon, pfIdx)
          ){

            /* 	   Hist1D["hMuonEtaRECO"]->Fill(mus_eta->at(i)); */
            /* 	   if( hasPFMatch(i, particleId::muon, pfIdx)) { */
            /* 	     Hist1D["hMuonEtaRECOPF"]->Fill(mus_eta->at(i)); */
            /* 	     Hist2D["hMuonPFRecoPT"]->Fill(mus_pt->at(i), mus_pt->at(i)-pfcand_pt->at(pfIdx)); */
            /* 	     Hist2D["hMuonPFRecoPTRatio"]->Fill(mus_pt->at(i), (mus_pt->at(i)-pfcand_pt->at(pfIdx))/mus_pt->at(i)); */
            RA4Muon.push_back(i);
            if(isFunnyMuon() || isFunnyElectron()) {   std::cout << " pass: 1" << std::endl; }
            /* 	   } */
        }
        else if((mus_isGlobalMuon->at(i) > 0 || mus_isTrackerMuon->at(i) > 0) 
                && mus_isPFMuon->at(i) > 0
                && mus_pt->at(i) >= LEPTONVETOPT
                //&&  fabs(mus_eta->at(i)) < 2.4
                && fabs(mus_eta->at(i)) <= 2.5
                // 26JUL2012 commented the following lines out
                && relIso < 0.2
                //&& (fabs(mus_dB->at(i)) < 0.2 || no_d0)
                && (fabs(d0PV) < 0.2 || no_d0)
                && fabs(getDZ(mus_tk_vx->at(i), mus_tk_vy->at(i), mus_tk_vz->at(i), mus_tk_px->at(i), mus_tk_py->at(i), mus_tk_pz->at(i), firstGoodVertex)) < 0.5
               ){
            if(isFunnyMuon() || isFunnyElectron()) {   std::cout << " passveto: 1" << std::endl; }
            RA4MuonVeto.push_back(i);
        } // RA4 Veto selection
        else {
            if(isFunnyMuon() || isFunnyElectron()) {   std::cout << " pass: 0" << std::endl; }
        }
    }// loop over reco muons.

    return RA4Muon;
}

vector<int> GetRA4Muon2011(vector<int> & RA4MuonVeto){
    vector<int> RA4Muon;
    for(unsigned int i=0;i<mus_pt->size();i++){

        float d0BS = mus_tk_d0dum->at(i)-beamSpot_x->at(0)*sin(mus_tk_phi->at(i))
            +beamSpot_y->at(0)*cos(mus_tk_phi->at(i));
        float CRI03 = (mus_tIso->at(i)+mus_ecalIso->at(i)+mus_hcalIso->at(i))/mus_pt->at(i);

        // RA4 selection
        if(mus_num_matches->at(i) > 1
                && mus_id_GlobalMuonPromptTight->at(i)> 0
                && mus_id_AllTrackerMuons->at(i)> 0
                && mus_pt->at(i) >= TIGHTLEPTONPT
                && fabs(mus_eta->at(i)) < 2.1
                && CRI03 < 0.1
                && fabs(d0BS) < 0.02
                && mus_tk_numvalhits->at(i) >= 11
                && mus_num_matches->at(i) >= 2
                && mus_tk_numpixelWthMeasr->at(i) >= 1
                //getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), els_px->at(i), els_py->at(i), els_pz->at(i), firstGoodVertex)
                //		&& fabs(mus_tk_vz->at(i)-pv_z->at(0)) < 1
          ){
            RA4Muon.push_back(i);
        }
        else if(mus_id_GlobalMuonPromptTight->at(i) > 0 // applies globalchi2 < 10 and nMuon hit > 0
                && mus_pt->at(i) >= LEPTONVETOPT
                && fabs(mus_eta->at(i)) < 2.5
                && CRI03 < 0.15 // follow RA6's Isolation.
                && fabs(d0BS) < 0.04 // x2 of normal. is 0.02 too small?
                && fabs(mus_tk_vz->at(i)-pv_z->at(0)) < 1 // x2 of normal
               ){
            RA4MuonVeto.push_back(i);
        } // RA4 selection
    }// loop over reco muons.
    return RA4Muon;
}

vector<int> GetRA4Elec(vector<int> & ElecVeto, std::string cutname="", int firstGoodVertex=0, bool isMC=false){
    vector<int> Elec;
    bool noPT= (cutname=="PT")?true:false;
    bool noEta= (cutname=="Eta")?true:false;
    bool noIso= (cutname=="Iso")?true:false;
    bool nodEtaIn= (cutname=="dEtaIn")?true:false;
    bool nodPhiIn= (cutname=="dPhiIn")?true:false;
    bool nosigmaIEtaIEta= (cutname=="sigmaIEtaIEta")?true:false;
    bool nohadOverEm= (cutname=="hadOverEm")?true:false;
    bool nod0PV= (cutname=="d0PV")?true:false;
    bool nodzPV= (cutname=="dzPV")?true:false;
    bool noOneOverEminusOneOverP= (cutname=="OneOverEminusOneOverP")?true:false;
    bool nohasMatchedConversion= (cutname=="hasMatchedConversion")?true:false;
    bool non_inner_layer= (cutname=="n_inner_layer")?true:false;
    for(int i =  int(els_pt->size())-1; i>=0; i--){

        float d0PV = els_d0dum->at(i)-pv_x->at(firstGoodVertex)*sin(els_tk_phi->at(i))
            +pv_y->at(firstGoodVertex)*cos(els_tk_phi->at(i));
        double max = els_PFphotonIsoR03->at(i) + els_PFneutralHadronIsoR03->at(i) - rho_kt6PFJetsForIsolation2011 * GetEffectiveArea(els_scEta->at(i), isMC);
        if(max<0.0) max=0;
        double relIso = (els_PFchargedHadronIsoR03->at(i) + max)/els_pt->at(i);

        if(isFunnyElectron() || isFunnyMuon()) {
            std::cout << "Electron:  run " << run << " lumi " << lumiblock << " event " << event
                << " pt " << els_pt->at(i) << " eta " << els_scEta->at(i)
                << " dEtaIn: " << els_dEtaIn->at(i)
                << " dPhiIn: " << els_dPhiIn->at(i)
                << " sigmaIEtaIEta "
                << els_sigmaIEtaIEta->at(i)
                << " H/E:" << els_hadOverEm->at(i)
                << " d0: " << d0PV
                //<< " dz: " << fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), els_px->at(i), els_py->at(i), els_pz->at(i), firstGoodVertex))
                << " dz: " << fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), cos(els_tk_phi->at(i))*els_tk_pt->at(i), sin(els_tk_phi->at(i))*els_tk_pt->at(i), els_tk_pz->at(i), firstGoodVertex))
                << " firstGoodVertex: " << firstGoodVertex
                << " |1/E-1/p| " << fabs(1./els_caloEnergy->at(i) - els_eOverPIn->at(i)/els_caloEnergy->at(i))
                << " :  ChargedIso "
                << els_PFchargedHadronIsoR03->at(i) << " PhotonIso "
                << els_PFphotonIsoR03->at(i) << " NeutralHadron Iso "
                << els_PFneutralHadronIsoR03->at(i)
                << " matched conversion: " <<   els_hasMatchedConversion->at(i)
                << " conversion: missing hits: " <<  els_n_inner_layer->at(i);
            //		<< std::endl;
        }

        int pfIdx=-1;
        // RA4 Selection
        if((els_pt->at(i) >= TIGHTLEPTONPT || noPT)
                && (fabs(els_scEta->at(i)) <= 2.5 || noEta)
                && (
                    //EB selection
                    (els_isEB->at(i)
                     && (fabs(els_dEtaIn->at(i)) < 0.004 || nodEtaIn)
                     && (fabs(els_dPhiIn->at(i)) < 0.06 || nodPhiIn)
                     && (els_sigmaIEtaIEta->at(i) < 0.01 || nosigmaIEtaIEta)
                     && (els_hadOverEm->at(i) < 0.12 || nohadOverEm)
                     && (fabs(d0PV) < 0.02 || nod0PV)
                     //&& (fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), els_px->at(i), els_py->at(i), els_pz->at(i), firstGoodVertex)) < 0.1 || nodzPV)
                     && (fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), cos(els_tk_phi->at(i))*els_tk_pt->at(i), sin(els_tk_phi->at(i))*els_tk_pt->at(i), els_tk_pz->at(i), firstGoodVertex)) < 0.1 || nodzPV)
                     && (fabs(1./els_caloEnergy->at(i) - els_eOverPIn->at(i)/els_caloEnergy->at(i)) < 0.05  || noOneOverEminusOneOverP)
                     && (relIso < 0.15 || noIso)
                     && (!els_hasMatchedConversion->at(i)  || nohasMatchedConversion)
                     && (els_n_inner_layer->at(i) <= 1 || non_inner_layer)
                     && hasPFMatch(i, particleId::electron, pfIdx)
                    )
                    ||
                    //EE selection
                    (els_isEE->at(i)
                     && (fabs(els_dEtaIn->at(i)) < 0.007 || nodEtaIn)
                     && (fabs(els_dPhiIn->at(i)) < 0.03 || nodPhiIn)
                     && (els_sigmaIEtaIEta->at(i) < 0.03 || nosigmaIEtaIEta)
                     && (els_hadOverEm->at(i) < 0.10 || nohadOverEm)
                     && (fabs(d0PV) < 0.02 || nod0PV)
                     //&& (fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), els_px->at(i), els_py->at(i), els_pz->at(i), firstGoodVertex)) < 0.1 || nodzPV)
                     && (fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), cos(els_tk_phi->at(i))*els_tk_pt->at(i), sin(els_tk_phi->at(i))*els_tk_pt->at(i), els_tk_pz->at(i), firstGoodVertex)) < 0.1 || nodzPV)
                     && (fabs(1./els_caloEnergy->at(i) - els_eOverPIn->at(i)/els_caloEnergy->at(i)) < 0.05  || noOneOverEminusOneOverP)
                     && (relIso < 0.15 || noIso)
                     && (!els_hasMatchedConversion->at(i)  || nohasMatchedConversion)
                     && (els_n_inner_layer->at(i) <= 1 || non_inner_layer)
                     && hasPFMatch(i, particleId::electron, pfIdx)
                    )
                    )
                    ) {
                        if( hasPFMatch(i, particleId::electron, pfIdx)) {

                            Elec.push_back(i);
                            if(isFunnyElectron()|| isFunnyMuon()) {	std::cout << " pass: 1" << std::endl;     }
                        }
                    }// RA4 Selection
        else if( els_pt->at(i) >= LEPTONVETOPT
                && fabs(els_scEta->at(i)) <= 2.5
                && (
                    //EB selection
                    (els_isEB->at(i)
                     //(fabs(els_scEta->at(i)) < 1.4442
                    && fabs(els_dEtaIn->at(i)) < 0.007
                    && fabs(els_dPhiIn->at(i)) < 0.8
                    && els_sigmaIEtaIEta->at(i) < 0.01
                    && els_hadOverEm->at(i) < 0.15
                    && fabs(d0PV) < 0.04
                    //&& fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), els_px->at(i), els_py->at(i), els_pz->at(i), firstGoodVertex)) < 0.2
                    && fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), cos(els_tk_phi->at(i))*els_tk_pt->at(i), sin(els_tk_phi->at(i))*els_tk_pt->at(i), els_tk_pz->at(i), firstGoodVertex)) < 0.2 
                    && relIso < 0.15
                    )
                     ||
                     //EE selection
                     //(fabs(els_scEta->at(i))>1.566
                        (els_isEE->at(i)
                         && fabs(els_dEtaIn->at(i)) < 0.01
                         && fabs(els_dPhiIn->at(i)) < 0.7 && els_sigmaIEtaIEta->at(i) < 0.03
                         && fabs(d0PV) < 0.04
                         //&& fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), els_px->at(i), els_py->at(i), els_pz->at(i), firstGoodVertex)) < 0.2
                         && fabs(getDZ(els_vx->at(i), els_vy->at(i), els_vz->at(i), cos(els_tk_phi->at(i))*els_tk_pt->at(i), sin(els_tk_phi->at(i))*els_tk_pt->at(i), els_tk_pz->at(i), firstGoodVertex)) < 0.2 
                         && relIso < 0.15
                        )
                        )
                        ){
                            ElecVeto.push_back(i);
                            if(isFunnyElectron()|| isFunnyMuon()) {	std::cout << " passveto: 1" << std::endl;      }
                        }
        else {
            if(isFunnyElectron()|| isFunnyMuon()) {	std::cout << " pass: 0" << std::endl;     }
        }
    }
    return Elec;
}

float GetEffectiveArea(float SCEta, bool isMC)
{
    float EffectiveArea;

    if(isMC) {
        if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.110;
        if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.130;
        if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.089;
        if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.130;
        if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.150;
        if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.160;
        if (fabs(SCEta) >= 2.4) EffectiveArea = 0.190;
    }
    else {
        //kEleGammaAndNeutralHadronIso03 from 2011 data
        //obtained from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.3&view=markup
        if (fabs(SCEta) >= 0.0 && fabs(SCEta) < 1.0 ) EffectiveArea = 0.100;
        if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) EffectiveArea = 0.120;
        if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) EffectiveArea = 0.085;
        if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.2 ) EffectiveArea = 0.110;
        if (fabs(SCEta) >= 2.2 && fabs(SCEta) < 2.3 ) EffectiveArea = 0.120;
        if (fabs(SCEta) >= 2.3 && fabs(SCEta) < 2.4 ) EffectiveArea = 0.120;
        if (fabs(SCEta) >= 2.4) EffectiveArea = 0.130;
    }
    return EffectiveArea;
}

vector<int> GetNIsoElec(){
    vector<int> Elec;
    for(int i = 0; i < int(els_pt->size());i++){
        float d0BS = els_d0dum->at(i)-beamSpot_x->at(0)*sin(els_tk_phi->at(i))
            +beamSpot_y->at(0)*cos(els_tk_phi->at(i));
        float EIso;
        if(fabs(els_scEta->at(i)) < 1.4442){
            EIso=GetMax(0,els_dr03EcalRecHitSumEt->at(i)-1);
        }
        else{
            EIso=els_dr03EcalRecHitSumEt->at(i);
        }
        double CRIso = (els_dr03TkSumPt->at(i)+EIso
                +els_dr03HcalTowerSumEt->at(i))/els_pt->at(i);
        d0BS=d0BS;
        // RA4 Selection
        if(els_pt->at(i) > TIGHTLEPTONPT
                && (fabs(els_scEta->at(i)) < 1.4442
                    ||fabs(els_scEta->at(i)) > 1.566)
                && fabs(els_scEta->at(i)) < 2.4
                && (els_simpleEleId80cIso->at(i) ==7
                    || els_simpleEleId80cIso->at(i) == 5)
                && fabs(els_vz->at(i)-pv_z->at(0)) < 1.0
                //&& fabs(d0BS) < 0.02
                && CRIso < 1.5){
            Elec.push_back(i);
        }// RA4 Selection
    }
    return Elec;
}

vector<int> GetNonDupElec(vector<int> RA4Elec){
    vector<int> CleanElec;
    for(int i=0;i<int(RA4Elec.size());i++){
        bool ADup=false;
        for(int j=0;j<int(CleanElec.size());j++){
            double DeltaP =sqrt(pow(els_px->at(RA4Elec[i])-els_px->at(CleanElec[j]),2)+					
                    pow(els_py->at(RA4Elec[i])-els_py->at(CleanElec[j]),2)+
                    pow(els_pz->at(RA4Elec[i])-els_pz->at(CleanElec[j]),2));
            if(DeltaP < 0.1){
                ADup = true;
            }
        }

        if(!ADup) {
            CleanElec.push_back(RA4Elec[i]);			
        }
    }
    return CleanElec;
}

vector<int> GetNIsoMuon(){
    vector<int> NIsoMuon;
    for(unsigned int i=0;i<mus_pt->size();i++){
        float d0BS = mus_tk_d0dum->at(i)-beamSpot_x->at(0)*sin(mus_tk_phi->at(i))
            +beamSpot_y->at(0)*cos(mus_tk_phi->at(i));
        float CRI03 = (mus_tIso->at(i)+mus_ecalIso->at(i)+mus_hcalIso->at(i))/mus_pt->at(i);
        d0BS=d0BS;
        // NIso selection	
        if(mus_id_GlobalMuonPromptTight->at(i)> 0
                && mus_id_AllTrackerMuons->at(i)> 0
                && mus_pt->at(i) >= TIGHTLEPTONPT
                && fabs(mus_eta->at(i)) < 2.1
                && CRI03 < 1.5
                //&& fabs(d0BS) < 0.02
                && fabs(mus_tk_vz->at(i)-pv_z->at(0)) < 1
                && mus_tk_numvalhits->at(i) >= 11
                && mus_num_matches->at(i) >= 1
                && mus_tk_numpixelWthMeasr->at(i) >= 1
          ){
            NIsoMuon.push_back(i);
        }
    }// loop over reco muons.
    return NIsoMuon;
}


vector<int> GetJets(vector<int> RA4Muon, vector<int> RA4Elec,
        vector<int> RA4MuonVeto, vector<int> RA4ElecVeto,
        double& HT, vector<int>& LooseBJet, vector<int>& MediumBJet,
        double etaMax=2.4, double ptMin=40, double minDeltaR=0.3){

    HT=0.0;

    vector<int> Jets; RA4MuonVeto=RA4MuonVeto; RA4ElecVeto=RA4ElecVeto;
    for(unsigned int i=0;i<jets_AK5PFclean_pt->size();i++){
        double rawRatio =(jets_AK5PFclean_rawPt->at(i)/jets_AK5PFclean_pt->at(i)); 
        if(jets_AK5PFclean_pt->at(i) >= ptMin
                && fabs(jets_AK5PFclean_eta->at(i)) <= etaMax
                && jets_AK5PFclean_neutralEmE->at(i)/jets_AK5PFclean_energy->at(i)/rawRatio < 0.99
                && jets_AK5PFclean_chgEmE->at(i)/jets_AK5PFclean_energy->at(i)/rawRatio < 0.99
                && jets_AK5PFclean_neutralHadE->at(i)/jets_AK5PFclean_energy->at(i)/rawRatio < 0.99
                && jets_AK5PFclean_chgHadE->at(i)/jets_AK5PFclean_energy->at(i)/rawRatio > 0
                && jets_AK5PFclean_chg_Mult->at(i) > 0
                && jets_AK5PFclean_chg_Mult->at(i)
                +jets_AK5PFclean_neutral_Mult->at(i)
                +jets_AK5PFclean_mu_Mult->at(i) > 1){
            // Elec XCleaning
            float DeltaR=0;
            float MinElecDeltaR=9.99;
            for(int iElec=0;iElec<int(RA4Elec.size());iElec++){
                DeltaR= sqrt( pow( GetDeltaPhi(jets_AK5PFclean_phi->at(i),els_phi->at(RA4Elec[iElec])),2)
                        +pow(jets_AK5PFclean_eta->at(i)-els_eta->at(RA4Elec[iElec]),2));
                if(DeltaR < MinElecDeltaR) MinElecDeltaR = DeltaR;
            }
            float MinElecVetoDeltaR=9.99;
            for(int iElecVeto=0;iElecVeto<int(RA4ElecVeto.size());iElecVeto++){
                DeltaR= sqrt( pow( GetDeltaPhi(jets_AK5PFclean_phi->at(i),els_phi->at(RA4ElecVeto[iElecVeto])),2)
                        +pow(jets_AK5PFclean_eta->at(i)-els_eta->at(RA4ElecVeto[iElecVeto]),2));
                if(DeltaR < MinElecVetoDeltaR) MinElecVetoDeltaR = DeltaR;
            }

            // Muon XCleaning
            DeltaR=0;
            float MinMuonDeltaR=9.99;
            for(int iMuon=0;iMuon<int(RA4Muon.size());iMuon++){
                DeltaR = sqrt(pow( GetDeltaPhi(jets_AK5PFclean_phi->at(i),mus_phi->at(RA4Muon[iMuon])),2)
                        +pow(jets_AK5PFclean_eta->at(i)-mus_eta->at(RA4Muon[iMuon]),2));
                if(DeltaR < MinMuonDeltaR) MinMuonDeltaR = DeltaR;
            }
            float MinMuonVetoDeltaR=9.99;
            for(int iMuonVeto=0;iMuonVeto<int(RA4MuonVeto.size());iMuonVeto++){
                DeltaR = sqrt(pow( GetDeltaPhi(jets_AK5PFclean_phi->at(i),mus_phi->at(RA4MuonVeto[iMuonVeto])),2)
                        +pow(jets_AK5PFclean_eta->at(i)-mus_eta->at(RA4MuonVeto[iMuonVeto]),2));
                if(DeltaR < MinMuonVetoDeltaR) MinMuonVetoDeltaR = DeltaR;
            }
            if(MinElecDeltaR >= minDeltaR && MinMuonDeltaR >= minDeltaR && MinElecVetoDeltaR >= minDeltaR && MinMuonVetoDeltaR >= minDeltaR){ 		  
                //	if(event==556) {
                //	  std::cout << i << " " << HT << " " << jets_AK5PFclean_pt->at(i) << std::endl; 
                HT=HT+jets_AK5PFclean_pt->at(i);
                Jets.push_back(i);
                //	  std::cout << i << " " << HT << " " << jets_AK5PFclean_pt->at(i) << std::endl; 
                //	}
                // Combined Secondary Vertex Loose=0.244, Medium = 0.679, Tight=0.898
                if(jets_AK5PFclean_btag_secVertexCombined->at(i) > 0.244){
                    LooseBJet.push_back(i);
                }
                if(jets_AK5PFclean_btag_secVertexCombined->at(i) > 0.679){
                    MediumBJet.push_back(i);
                }
            }
        }
    }

    return Jets;
}

vector<int> GetJetsCHS(vector<int> RA4Muon, vector<int> RA4Elec,
                       vector<int> RA4MuonVeto, vector<int> RA4ElecVeto,
                       double& HT, vector<int>& LooseBJet, vector<int>& MediumBJet,
                       double etaMax=2.4, double ptMin=40, double minDeltaR=0.3){

    HT=0.0;

    vector<int> Jets; RA4MuonVeto=RA4MuonVeto; RA4ElecVeto=RA4ElecVeto;
    for(unsigned int i=0;i<jets_AK5PF_pt->size();i++){
        double rawRatio =(jets_AK5PF_rawPt->at(i)/jets_AK5PF_pt->at(i)); 
        if(jets_AK5PF_pt->at(i) >= ptMin
                && fabs(jets_AK5PF_eta->at(i)) <= etaMax
                && jets_AK5PF_neutralEmE->at(i)/jets_AK5PF_energy->at(i)/rawRatio < 0.99
                && jets_AK5PF_chgEmE->at(i)/jets_AK5PF_energy->at(i)/rawRatio < 0.99
                && jets_AK5PF_neutralHadE->at(i)/jets_AK5PF_energy->at(i)/rawRatio < 0.99
                && jets_AK5PF_chgHadE->at(i)/jets_AK5PF_energy->at(i)/rawRatio > 0
                && jets_AK5PF_chg_Mult->at(i) > 0
                && jets_AK5PF_chg_Mult->at(i)
                +jets_AK5PF_neutral_Mult->at(i)
                +jets_AK5PF_mu_Mult->at(i) > 1){
            // Elec XCleaning
            float DeltaR=0;
            float MinElecDeltaR=9.99;
            for(int iElec=0;iElec<int(RA4Elec.size());iElec++){
                DeltaR= sqrt( pow( GetDeltaPhi(jets_AK5PF_phi->at(i),els_phi->at(RA4Elec[iElec])),2)
                        +pow(jets_AK5PF_eta->at(i)-els_eta->at(RA4Elec[iElec]),2));
                if(DeltaR < MinElecDeltaR) MinElecDeltaR = DeltaR;
            }
            float MinElecVetoDeltaR=9.99;
            for(int iElecVeto=0;iElecVeto<int(RA4ElecVeto.size());iElecVeto++){
                DeltaR= sqrt( pow( GetDeltaPhi(jets_AK5PF_phi->at(i),els_phi->at(RA4ElecVeto[iElecVeto])),2)
                        +pow(jets_AK5PF_eta->at(i)-els_eta->at(RA4ElecVeto[iElecVeto]),2));
                if(DeltaR < MinElecVetoDeltaR) MinElecVetoDeltaR = DeltaR;
            }

            // Muon XCleaning
            DeltaR=0;
            float MinMuonDeltaR=9.99;
            for(int iMuon=0;iMuon<int(RA4Muon.size());iMuon++){
                DeltaR = sqrt(pow( GetDeltaPhi(jets_AK5PF_phi->at(i),mus_phi->at(RA4Muon[iMuon])),2)
                        +pow(jets_AK5PF_eta->at(i)-mus_eta->at(RA4Muon[iMuon]),2));
                if(DeltaR < MinMuonDeltaR) MinMuonDeltaR = DeltaR;
            }
            float MinMuonVetoDeltaR=9.99;
            for(int iMuonVeto=0;iMuonVeto<int(RA4MuonVeto.size());iMuonVeto++){
                DeltaR = sqrt(pow( GetDeltaPhi(jets_AK5PF_phi->at(i),mus_phi->at(RA4MuonVeto[iMuonVeto])),2)
                        +pow(jets_AK5PF_eta->at(i)-mus_eta->at(RA4MuonVeto[iMuonVeto]),2));
                if(DeltaR < MinMuonVetoDeltaR) MinMuonVetoDeltaR = DeltaR;
            }
            if(MinElecDeltaR >= minDeltaR && MinMuonDeltaR >= minDeltaR && MinElecVetoDeltaR >= minDeltaR && MinMuonVetoDeltaR >= minDeltaR){ 		  
                //	if(event==556) {
                //	  std::cout << i << " " << HT << " " << jets_AK5PF_pt->at(i) << std::endl; 
                HT=HT+jets_AK5PF_pt->at(i);
                Jets.push_back(i);
                //	  std::cout << i << " " << HT << " " << jets_AK5PF_pt->at(i) << std::endl; 
                //	}
                // Combined Secondary Vertex Loose=0.244, Medium = 0.679, Tight=0.898
                if(jets_AK5PF_btag_secVertexCombined->at(i) > 0.244){
                    LooseBJet.push_back(i);
                }
                if(jets_AK5PF_btag_secVertexCombined->at(i) > 0.679){
                    MediumBJet.push_back(i);
                }
            }
        }
    }

    return Jets;
}


#endif //OBJECTSELECTOR_SYNC_H
