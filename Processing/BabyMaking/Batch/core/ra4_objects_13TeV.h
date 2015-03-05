#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <iostream> 

#define MinSignalLeptonPt 20
#define MinVetoLeptonPt 15
#define MinJetPt 40

using namespace std;

double getDZ(double vx, double vy, double vz, double px, double py, double pz, int firstGoodVertex)
{
    return vz - pv_z->at(firstGoodVertex) 
        -((vx-pv_x->at(firstGoodVertex))*px+(vy-pv_y->at(firstGoodVertex))*py)*pz/(px*px+py*py); 
}

/////////////////////////////////////////////////////////////////////////
////////////////////////////////  MUONS  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////
float GetMuonIsolation(int imu)
{
    if(imu >= (int)mus_pt->size()) return -999;
    double sumEt = mus_pfIsolationR03_sumNeutralHadronEt->at(imu) + mus_pfIsolationR03_sumPhotonEt->at(imu) 
        - 0.5*mus_pfIsolationR03_sumPUPt->at(imu);
    if(sumEt<0.0) sumEt=0.0;
    return (mus_pfIsolationR03_sumChargedHadronPt->at(imu) + sumEt)/mus_pt->at(imu);
}

bool IsBasicMuon(int imu)
{
    if(imu >= (int)mus_pt->size()) return false;

    float d0PV = mus_tk_d0dum->at(imu)-pv_x->at(0)*sin(mus_tk_phi->at(imu))+pv_y->at(0)*cos(mus_tk_phi->at(imu));

    return (mus_isGlobalMuon->at(imu) > 0
            && mus_isPF->at(imu)
            && mus_id_GlobalMuonPromptTight->at(imu)> 0 
            && mus_tk_LayersWithMeasurement->at(imu) > 5
            && mus_tk_numvalPixelhits->at(imu) > 0
            && mus_numberOfMatchedStations->at(imu) > 1
            //&& fabs(mus_dB->at(imu)) < 0.02
            && fabs(d0PV) < 0.02
            && fabs(getDZ(mus_tk_vx->at(imu), mus_tk_vy->at(imu), mus_tk_vz->at(imu), mus_tk_px->at(imu), 
                    mus_tk_py->at(imu), mus_tk_pz->at(imu), 0)) < 0.5
            && mus_pt->at(imu) >= MinSignalLeptonPt
            && fabs(mus_eta->at(imu)) <= 2.4);
}

bool IsSignalMuon(int imu)
{
    if(imu >= (int)mus_pt->size()) return false;

    float relIso = GetMuonIsolation(imu);  
    return (IsBasicMuon(imu) && relIso < 0.12); 
}

bool IsVetoMuon(int imu)
{
    if(imu >= (int)mus_pt->size()) return false;
    if(IsSignalMuon(imu)) return false; // not signal muon

    float relIso = GetMuonIsolation(imu);

    return ((mus_isGlobalMuon->at(imu) >0 || mus_isTrackerMuon->at(imu) >0)
            && mus_isPF->at(imu) 
            && fabs(getDZ(mus_tk_vx->at(imu), mus_tk_vy->at(imu), mus_tk_vz->at(imu), mus_tk_px->at(imu), 
                    mus_tk_py->at(imu), mus_tk_pz->at(imu), 0)) < 0.5 
            && mus_pt->at(imu) >= MinVetoLeptonPt
            && fabs(mus_eta->at(imu)) <= 2.5
            && relIso < 0.2);
}


vector<int> GetMuons(bool doSignal)
{
    vector<int> muons;
    for(int index=0; index<(int)mus_pt->size(); index++)
        if(doSignal)
        {
            if(IsSignalMuon(index)) muons.push_back(index);
        } else 
        {
            if(IsVetoMuon(index)) muons.push_back(index);
        }
    return muons;
}


/////////////////////////////////////////////////////////////////////////
//////////////////////////////  ELECTRONS  //////////////////////////////
/////////////////////////////////////////////////////////////////////////
float GetElectronIsolation(int iel)
{
    float absiso = els_pfIsolationR03_sumChargedHadronPt->at(iel) + std::max(0.0 , els_pfIsolationR03_sumNeutralHadronEt->at(iel) + els_pfIsolationR03_sumPhotonEt->at(iel) - 0.5 * els_pfIsolationR03_sumPUPt->at(iel) );
    return absiso/els_pt->at(iel);
}

bool IsBasicElectron(int iel)
{
    if(iel >= (int)els_pt->size()) return false;

    float d0PV = els_d0dum->at(iel)-pv_x->at(0)*sin(els_tk_phi->at(iel))+pv_y->at(0)*cos(els_tk_phi->at(iel));

    return (els_pt->at(iel) > MinSignalLeptonPt
            && fabs(els_scEta->at(iel)) < 2.5
            && !els_hasMatchedConversion->at(iel)
            && els_n_inner_layer->at(iel) <= 1 // FIXME why this does not exist in PHYS14? 
            && fabs(getDZ(els_vx->at(iel), els_vy->at(iel), els_vz->at(iel), cos(els_tk_phi->at(iel))*els_tk_pt->at(iel), 
                    sin(els_tk_phi->at(iel))*els_tk_pt->at(iel), els_tk_pz->at(iel), 0)) < 0.1
            && fabs(1./els_caloEnergy->at(iel) - els_eOverPIn->at(iel)/els_caloEnergy->at(iel)) < 0.05 
            && fabs(d0PV) < 0.02 
            && ((els_isEB->at(iel) // Endcap selection
                    && fabs(els_dEtaIn->at(iel)) < 0.004
                    && fabs(els_dPhiIn->at(iel)) < 0.06
                    && els_sigmaIEtaIEta->at(iel) < 0.01
                    && els_hadOverEm->at(iel) < 0.12 ) ||
                (els_isEE->at(iel)  // Barrel selection
                 && fabs(els_dEtaIn->at(iel)) < 0.007
                 && fabs(els_dPhiIn->at(iel)) < 0.03
                 && els_sigmaIEtaIEta->at(iel) < 0.03
                 && els_hadOverEm->at(iel) < 0.10 ))
           );
}

bool IsSignalElectron(int iel)
{
    if(iel >= (int)els_pt->size()) return false;

    double relIso = GetElectronIsolation(iel);
    return (IsBasicElectron(iel) && relIso < 0.15);
}

bool IsVetoElectron(int iel)
{
    if(iel >= (int)els_pt->size()) return false;
    if(IsSignalElectron(iel)) return false; // not signal electron 

    float d0PV = els_d0dum->at(iel)-pv_x->at(0)*sin(els_tk_phi->at(iel))+pv_y->at(0)*cos(els_tk_phi->at(iel));
    double relIso = GetElectronIsolation(iel);

    return (els_pt->at(iel) > MinVetoLeptonPt
            && fabs(els_scEta->at(iel)) < 2.5
            && relIso < 0.15
            && fabs(getDZ(els_vx->at(iel), els_vy->at(iel), els_vz->at(iel), cos(els_tk_phi->at(iel))*els_tk_pt->at(iel), 
                    sin(els_tk_phi->at(iel))*els_tk_pt->at(iel), els_tk_pz->at(iel), 0)) < 0.2
            && fabs(d0PV) < 0.04 
            && ((els_isEB->at(iel) // Endcap selection
                    && fabs(els_dEtaIn->at(iel)) < 0.007
                    && fabs(els_dPhiIn->at(iel)) < 0.8
                    && els_sigmaIEtaIEta->at(iel) < 0.01
                    && els_hadOverEm->at(iel) < 0.15) ||
                (els_isEE->at(iel)  // Barrel selection
                 && fabs(els_dEtaIn->at(iel)) < 0.01
                 && fabs(els_dPhiIn->at(iel)) < 0.7
                 && els_sigmaIEtaIEta->at(iel) < 0.03))
           );  
}

vector<int> GetElectrons(bool doSignal)
{
    vector<int> electrons;
    for(int index=0; index<(int)els_pt->size(); index++)
        if(doSignal)
        {
            if(IsSignalElectron(index)) electrons.push_back(index);
        }   
        else 
        {
            if(IsVetoElectron(index)) electrons.push_back(index);
        }
    return electrons;
}

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

        if(jets_AK4_pt->at(ijet) > 40) 
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

/*
int GetTrueMuon(int index, int &momID, double &closest_dR)
{
    if(index < 0 || index >= static_cast<int>(mus_eta->size())) return -1;

    int closest_imc = -1, idLepton = 0; 
    double dR = 9999.; closest_dR = 9999.;
    double MCEta, MCPhi;
    double RecEta = mus_eta->at(index), RecPhi = mus_phi->at(index);
    for(int imc=0; imc < (int)mc_mus_id->size(); imc++)
    {
        MCEta = mc_mus_eta->at(imc); MCPhi = mc_mus_phi->at(imc);
        dR = sqrt(pow(RecEta-MCEta,2) + pow(RecPhi-MCPhi,2));
        if(dR < closest_dR) 
        {
            closest_dR = dR;
            closest_imc = imc;
        }
    }
    if(closest_imc >= 0)
    {
        idLepton = static_cast<int>(mc_mus_id->at(closest_imc));
        momID = static_cast<int>(mc_mus_mother_id->at(closest_imc));
        if(idLepton == momID) momID = static_cast<int>(mc_mus_ggrandmother_id->at(closest_imc));
    } 
    else 
    {
        closest_imc = GetTrueParticle(RecEta, RecPhi, closest_dR);
        if(closest_imc >= 0)
        {
            momID = static_cast<int>(mc_doc_mother_id->at(closest_imc));
            idLepton = static_cast<int>(mc_doc_id->at(closest_imc));
        } 
        else 
        {
            momID = 0;
            idLepton = 0;
        }
    }
    return idLepton;
}
*/
/*
int GetTrueElectron(int index, int &momID, double &closest_dR)
{
    if(index < 0 || index >= static_cast<int>(els_eta->size())) return -1;

    int closest_imc = -1, idLepton = 0; 
    double dR = 9999.; closest_dR = 9999.;
    double MCEta, MCPhi;
    double RecEta = els_eta->at(index), RecPhi = els_phi->at(index);
    for(int imc=0; imc < (int)mc_electrons_id->size(); imc++)
    {
        MCEta = mc_electrons_eta->at(imc); MCPhi = mc_electrons_phi->at(imc);
        dR = sqrt(pow(RecEta-MCEta,2) + pow(RecPhi-MCPhi,2));
        if(dR < closest_dR) 
        {
            closest_dR = dR;
            closest_imc = imc;
        }
    }
    if(closest_imc >= 0)
    {
        idLepton = static_cast<int>(mc_electrons_id->at(closest_imc));
        momID = static_cast<int>(mc_electrons_mother_id->at(closest_imc));
        if(idLepton == momID) momID = static_cast<int>(mc_electrons_ggrandmother_id->at(closest_imc));
    } 
    else 
    {
        closest_imc = GetTrueParticle(RecEta, RecPhi, closest_dR);
        if(closest_imc >= 0)
        {
            momID = static_cast<int>(mc_doc_mother_id->at(closest_imc));
            idLepton = static_cast<int>(mc_doc_id->at(closest_imc));
        } 
        else 
        {
            momID = 0;
            idLepton = 0;
        }
    }
    return idLepton;
}
*/

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

