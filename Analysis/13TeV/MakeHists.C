#include <iostream>

#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TMath.h"

#include "babytree_manuel.h"
#include "PassSelection.h"

using namespace std;
bool    DoLog           = 1;
int     FatjetpTthres   = 0;
int     mjthres         = 0;

//
//TH1F initialization
//
TH1F* InitTH1F(char* Name, char* Title, int Nbins, double XMin, double XMax)
{
    TH1F *h1 = new TH1F(Name, Title, Nbins, XMin, XMax);
    h1->Sumw2();
    return h1;
}
//
//TH2F initialization
//
TH2F* InitTH2F(char* Name, char* Title, int NXbins, double XMin, double XMax, int NYbins, double YMin, double YMax)
{
    TH2F *h2 = new TH2F(Name, Title, NXbins, XMin, XMax, NYbins, YMin, YMax);
    h2->Sumw2();
    return h2;
}

//
// Fill TH1F
//
void FillTH1F(TH1F* &h1, double var, double weight)
{
    if(var >= h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())) 
        var=h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())-0.00001;
    if(var < h1->GetXaxis()->GetBinLowEdge(1)) 
        var=h1->GetXaxis()->GetBinLowEdge(1)+0.00001;
    h1->Fill(var, weight);        
}
void FillTH1FAll(TH1F* h1[7], int NFJ, double var, double weight)
{
    FillTH1F(h1[NFJ], var,  weight);
    FillTH1F(h1[6], var,  weight);
}

//
// Fill TH2F
//
void FillTH2F(TH2F* &h2, double varX, double varY, double weight){
    if(varX >= h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins())) 
        varX=h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins())-0.00001;
    if(varY >= h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins())) 
        varY=h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins())-0.00001;
    if(varX < h2->GetXaxis()->GetBinLowEdge(1)) 
        varX=h2->GetXaxis()->GetBinLowEdge(1)+0.00001;
    if(varY < h2->GetYaxis()->GetBinLowEdge(1)) 
        varY=h2->GetYaxis()->GetBinLowEdge(1)+0.00001;
    h2->Fill(varX, varY, weight);        
}
void FillTH2FAll(TH2F* h2[7], int NFJ, double varX, double varY, double weight)
{
    FillTH2F(h2[NFJ], varX,  varY, weight);
    FillTH2F(h2[6], varX,  varY, weight);
}

//
// loading PU reweigting histogram
//
double nPUScaleFactor2012(TH1F* h1PU, float npu){
    // Get PU scale factor histogram
    double mynpu = TMath::Min(npu,(float)49.499);
    Int_t npuxbin = h1PU->GetXaxis()->FindBin(mynpu);
    return h1PU->GetBinContent(npuxbin);
}

float getDPhi(float phi1, float phi2) 
{ 
    float absdphi = abs(phi1-phi2);
    if(absdphi < TMath::Pi()) return absdphi;
    else return (2*TMath::Pi() - absdphi);
    // Ryan's : http://www.aidansean.com/cheatsheets/?p=105
    //return TMath::Abs(TMath::Abs(TMath::Abs(phi1 - phi2) - TMath::Pi())-TMath::Pi());  

}
float getDEta(float eta1, float eta2)
{
    return abs(eta1-eta2);
}
float getDR(float dphi, float deta)
{
    return TMath::Sqrt(dphi*dphi+deta*deta);
}
float getDR(float eta1, float eta2, float phi1, float phi2)
{
    return getDR(getDPhi(phi1, phi2), eta1-eta2);
}

//
float getISRSF(float ISRpT)
{
    // ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMST2ccMadgraph8TeV
    if( ISRpT>250) return 0.8;
    else if( ISRpT>150 ) return 0.9;
    else if( ISRpT>120 ) return 0.95;
    else return 1;
}

//
// main function
//
void MakeHists(TChain *ch, char* Region) 
{ 

    InitBaby(ch); 
    TString ChainName = ch->GetTitle();
    cout << "[MJ Analysis] " << ChainName << endl;  

    //
    // Define histograms
    //
    TH1F *h1_yields[7];
    TH1F *h1_dPhitop[7]; // gen quantities
    TH1F *h1_dRlep[7], *h1_dPhiMET[7], *h1_dRbmin[7], *h1_dPhiMETlep[7]; // temporary histograms 
    TH1F *h1_MJ[7], *h1_mT[7], *h1_Nskinnyjet[7], *h1_Ncsvm[7], *h1_MJ_ISR[7],
         *h1_muspT[7], *h1_muspTminusMET[7], *h1_musEta[7], *h1_musPhi[7], 
         *h1_elspT[7], *h1_elspTminusMET[7], *h1_elsEta[7], *h1_elsPhi[7], 
         *h1_mj[7], *h1_FatjetPt[7], *h1_FatjetEta[7], 
         *h1_FatjetPt1[7],  *h1_FatjetPt2[7],  *h1_FatjetPt3[7],  *h1_FatjetPt4[7],
         *h1_mj1[7],  *h1_mj2[7],  *h1_mj3[7],  *h1_mj4[7], 
         *h1_mj08_1[7],  *h1_mj08_2[7],  *h1_mj08_3[7],  *h1_mj08_4[7], 
         *h1_mj1OverMJ[7],  *h1_mj2OverMJ[7],  *h1_mj3OverMJ[7],  *h1_mj4OverMJ[7], 
         *h1_N1[7],  *h1_N2[7],  *h1_N3[7],  *h1_N4[7], 
         *h1_mjOverPt1[7],  *h1_mjOverPt2[7],  *h1_mjOverPt3[7],  *h1_mjOverPt4[7], 
         *h1_mj3overmj2[7], *h1_mj2overmj1[7],
         *h1_FatjetPhi1[7], *h1_FatjetPhi2[7], *h1_FatjetPhi3[7], *h1_FatjetPhi4[7],
         *h1_FatjetEta1[7], *h1_FatjetEta2[7], *h1_FatjetEta3[7], *h1_FatjetEta4[7],
         *h1_dRFJ[7], *h1_dPhiFJ[7], *h1_dEtaFJ[7], *h1_mindPhibb[7],
         *h1_HT[7], *h1_MET[7], *h1_METPhi[7], *h1_METx[7], *h1_METy[7], *h1_DPhi[7], *h1_Nfatjet[7], *h1_WpT[7];
    TH2F *h2_mj1vsmj2[7], *h2_mj2vsmj3[7], *h2_mj3vsmj4[7];
    TH2F *h2_mj1vspt1[7], *h2_mj2vspt2[7], *h2_mj3vspt3[7], *h2_mj4vspt4[7];
    TH2F *h2_HTMET[7], *h2_MJmT[7];
    TH2F *h2_HTmT[7], *h2_MJMET[7];
    TH2F *h2_HTMJ[7], *h2_METmT[7];
    TH2F *h2_dPhidpTtop[7], *h2_dPhidpTttbar[7];

    for(int i=0; i<7; i++) 
    {   
        //
        // yields : for table                   
        //
        h1_yields[i] = InitTH1F( Form("h1_%s_yields_%ifatjet", ch->GetTitle(), i), 
                                 Form("h1_%s_yields_%ifatjet", ch->GetTitle(), i), 
                                 2, 0, 2); // bin1 for e, bin2 for mu
        //
        // h1                    
        //
        h1_dPhitop[i] = InitTH1F( Form("h1_%s_dPhitop_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_dPhitop_%ifatjet", ch->GetTitle(), i), 
                                10, 0, TMath::Pi());
        h1_dRlep[i] = InitTH1F( Form("h1_%s_dRlep_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_dRlep_%ifatjet", ch->GetTitle(), i), 
                                20, 0, 4);
        h1_dPhiMET[i] = InitTH1F( Form("h1_%s_dPhiMET_%ifatjet", ch->GetTitle(), i), 
                                  Form("h1_%s_dPhiMET_%ifatjet", ch->GetTitle(), i), 
                                  20, 0, TMath::Pi());
        h1_dRbmin[i] = InitTH1F( Form("h1_%s_dRbmin_%ifatjet", ch->GetTitle(), i), 
                                 Form("h1_%s_dRbmin_%ifatjet", ch->GetTitle(), i), 
                                 20, 0, 4);
        h1_dPhiMETlep[i] = InitTH1F( Form("h1_%s_dPhiMETlep_%ifatjet", ch->GetTitle(), i), 
                                     Form("h1_%s_dPhiMETlep_%ifatjet", ch->GetTitle(), i), 
                                     20, 0, TMath::Pi());
        h1_MJ[i] = InitTH1F( Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 2000);
        h1_MJ_ISR[i] = InitTH1F( Form("h1_%s_MJ_ISR_%ifatjet", ch->GetTitle(), i), 
                                 Form("h1_%s_MJ_ISR_%ifatjet", ch->GetTitle(), i), 
                                 20, 0, 2000);
        h1_mj[i] = InitTH1F( Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             40, 0, 800);
        h1_Nfatjet[i] = InitTH1F( Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  11, -0.5, 10.5);
        h1_Nskinnyjet[i] = InitTH1F( Form("h1_%s_Nskinnyjet_%ifatjet", ch->GetTitle(), i), 
                                     Form("h1_%s_Nskinnyjet_%ifatjet", ch->GetTitle(), i), 
                                     16, -0.5, 15.5);
        h1_Ncsvm[i] = InitTH1F( Form("h1_%s_Ncsvm_%ifatjet", ch->GetTitle(), i), 
                                     Form("h1_%s_Ncsvm_%ifatjet", ch->GetTitle(), i), 
                                     7, -0.5, 6.5);
        h1_mT[i] = InitTH1F( Form("h1_%s_mT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mT_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 200);
                             16, 0, 800);
        h1_muspT[i] = InitTH1F( Form("h1_%s_muspT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_muspT_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 200);
        h1_muspTminusMET[i] = InitTH1F( Form("h1_%s_muspTminusMET_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_muspTminusMET_%ifatjet", ch->GetTitle(), i), 
                             //20, 100, 400);
                             20, -2, 8);
        h1_musEta[i] = InitTH1F( Form("h1_%s_musEta_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_musEta_%ifatjet", ch->GetTitle(), i), 
                                20, -2.5, 2.5);
                                //100, -2.1, 2.1);
        h1_musPhi[i] = InitTH1F( Form("h1_%s_musPhi_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_musPhi_%ifatjet", ch->GetTitle(), i), 
                                20, -TMath::Pi(), TMath::Pi());
        h1_elspT[i] = InitTH1F( Form("h1_%s_elspT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_elspT_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 200);
        h1_elspTminusMET[i] = InitTH1F( Form("h1_%s_elspTminusMET_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_elspTminusMET_%ifatjet", ch->GetTitle(), i), 
                             //20, 100, 400);
                             20, -2, 8);
        h1_elsEta[i] = InitTH1F( Form("h1_%s_elsEta_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_elsEta_%ifatjet", ch->GetTitle(), i), 
                                20, -2.5, 2.5);
                                //100, -2.1, 2.1);
        h1_elsPhi[i] = InitTH1F( Form("h1_%s_elsPhi_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_elsPhi_%ifatjet", ch->GetTitle(), i), 
                                20, -TMath::Pi(), TMath::Pi());
        h1_HT[i] = InitTH1F( Form("h1_%s_HT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_HT_%ifatjet", ch->GetTitle(), i), 
                             //20, 350, 1350);
                             28, 500, 4000);
        h1_MET[i] = InitTH1F( Form("h1_%s_MET_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MET_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             18, 100, 1000);
        h1_METPhi[i] = InitTH1F( Form("h1_%s_METPhi_%ifatjet", ch->GetTitle(), i), 
                                 Form("h1_%s_METPhi_%ifatjet", ch->GetTitle(), i), 
                                 20, -TMath::Pi(), TMath::Pi());
        h1_METx[i] = InitTH1F( Form("h1_%s_METx_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_METx_%ifatjet", ch->GetTitle(), i), 
                               100, -300, 300);
        h1_METy[i] = InitTH1F( Form("h1_%s_METy_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_METy_%ifatjet", ch->GetTitle(), i), 
                               100, -300, 300);
        h1_DPhi[i] = InitTH1F( Form("h1_%s_DPhi_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_DPhi_%ifatjet", ch->GetTitle(), i), 
                               32, 0, TMath::Pi());
        h1_dRFJ[i] = InitTH1F( Form("h1_%s_dRFJ_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dRFJ_%ifatjet", ch->GetTitle(), i), 
                               20, 1, 5);
        h1_dPhiFJ[i] = InitTH1F( Form("h1_%s_dPhiFJ_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dPhiFJ_%ifatjet", ch->GetTitle(), i), 
                               32, 0, TMath::Pi());
        h1_dEtaFJ[i] = InitTH1F( Form("h1_%s_dEtaFJ_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dEtaFJ_%ifatjet", ch->GetTitle(), i), 
                               20, 0, 4);
        h1_FatjetPt[i] = InitTH1F( Form("h1_%s_FatjetPt_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_FatjetPt_%ifatjet", ch->GetTitle(), i), 
                                   20, 0, 800);
        h1_FatjetPt1[i] = InitTH1F( Form("h1_%s_FatjetPt1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt1_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 2000);
        h1_FatjetPt2[i] = InitTH1F( Form("h1_%s_FatjetPt2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt2_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1200);
        h1_FatjetPt3[i] = InitTH1F( Form("h1_%s_FatjetPt3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt3_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 800);
        h1_FatjetPt4[i] = InitTH1F( Form("h1_%s_FatjetPt4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt4_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 500);
        h1_mj1[i] = InitTH1F( Form("h1_%s_mj1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj1_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 500);
        h1_mj2[i] = InitTH1F( Form("h1_%s_mj2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj2_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 250);
        h1_mj3[i] = InitTH1F( Form("h1_%s_mj3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj3_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 200);
        h1_mj4[i] = InitTH1F( Form("h1_%s_mj4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj4_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 100);
        h1_mj08_1[i] = InitTH1F( Form("h1_%s_mj08_1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj08_1_%ifatjet", ch->GetTitle(), i), 
                                    12, 0, 600);
        h1_mj08_2[i] = InitTH1F( Form("h1_%s_mj08_2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj08_2_%ifatjet", ch->GetTitle(), i), 
                                    12, 0, 300);
        h1_mj08_3[i] = InitTH1F( Form("h1_%s_mj08_3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj08_3_%ifatjet", ch->GetTitle(), i), 
                                    12, 0, 240);
        h1_mj08_4[i] = InitTH1F( Form("h1_%s_mj08_4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj08_4_%ifatjet", ch->GetTitle(), i), 
                                    12, 0, 120);
        h1_N1[i] = InitTH1F( Form("h1_%s_N1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_N1_%ifatjet", ch->GetTitle(), i), 
                                    5, 0.5, 5.5);
        h1_N2[i] = InitTH1F( Form("h1_%s_N2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_N2_%ifatjet", ch->GetTitle(), i), 
                                    5, 0.5, 5.5);
        h1_N3[i] = InitTH1F( Form("h1_%s_N3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_N3_%ifatjet", ch->GetTitle(), i), 
                                    5, 0.5, 5.5);
        h1_N4[i] = InitTH1F( Form("h1_%s_N4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_N4_%ifatjet", ch->GetTitle(), i), 
                                    5, 0.5, 5.5);
        h1_mjOverPt1[i] = InitTH1F( Form("h1_%s_mjOverPt1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mjOverPt1_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mjOverPt2[i] = InitTH1F( Form("h1_%s_mjOverPt2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mjOverPt2_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mjOverPt3[i] = InitTH1F( Form("h1_%s_mjOverPt3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mjOverPt3_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mjOverPt4[i] = InitTH1F( Form("h1_%s_mjOverPt4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mjOverPt4_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mj1OverMJ[i] = InitTH1F( Form("h1_%s_mj1OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj1OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mj2OverMJ[i] = InitTH1F( Form("h1_%s_mj2OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj2OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mj3OverMJ[i] = InitTH1F( Form("h1_%s_mj3OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj3OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mj4OverMJ[i] = InitTH1F( Form("h1_%s_mj4OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj4OverMJ_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mj3overmj2[i] = InitTH1F( Form("h1_%s_mj3overmj2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj3overmj2_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_mj2overmj1[i] = InitTH1F( Form("h1_%s_mj2overmj1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj2overmj1_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 1);
        h1_FatjetPhi1[i] = InitTH1F( Form("h1_%s_FatjetPhi1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi1_%ifatjet", ch->GetTitle(), i), 
                                    32, -TMath::Pi(), TMath::Pi());
        h1_FatjetPhi2[i] = InitTH1F( Form("h1_%s_FatjetPhi2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi2_%ifatjet", ch->GetTitle(), i), 
                                    32, -TMath::Pi(), TMath::Pi());
        h1_FatjetPhi3[i] = InitTH1F( Form("h1_%s_FatjetPhi3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi3_%ifatjet", ch->GetTitle(), i), 
                                    32, -TMath::Pi(), TMath::Pi());
        h1_FatjetPhi4[i] = InitTH1F( Form("h1_%s_FatjetPhi4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi4_%ifatjet", ch->GetTitle(), i), 
                                    32, -TMath::Pi(), TMath::Pi());
        h1_FatjetEta1[i] = InitTH1F( Form("h1_%s_FatjetEta1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta1_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_FatjetEta2[i] = InitTH1F( Form("h1_%s_FatjetEta2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta2_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_FatjetEta3[i] = InitTH1F( Form("h1_%s_FatjetEta3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta3_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_FatjetEta4[i] = InitTH1F( Form("h1_%s_FatjetEta4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta4_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_FatjetEta[i] = InitTH1F( Form("h1_%s_FatjetEta_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta_%ifatjet", ch->GetTitle(), i), 
                                    20, -2.5, 2.5);
        h1_WpT[i] = InitTH1F( Form("h1_%s_WpT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_WpT_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             20, 0, 1000);
        h1_mindPhibb[i] = InitTH1F( Form("h1_%s_mindPhibb_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_mindPhibb_%ifatjet", ch->GetTitle(), i), 
                                10, 0, TMath::Pi());
        //
        // h2                    
        //
        h2_dPhidpTtop[i] =   InitTH2F(Form("h2_%s_dPhidpTtop_%ifatjet", ch->GetTitle(), i),
                                      Form("h2_%s_dPhidpTtop_%ifatjet", ch->GetTitle(), i), 
                                      10, 0, TMath::Pi(), 10, 0, 500);
        h2_dPhidpTttbar[i] =   InitTH2F(Form("h2_%s_dPhidpTttbar_%ifatjet", ch->GetTitle(), i),
                                      Form("h2_%s_dPhidpTttbar_%ifatjet", ch->GetTitle(), i), 
                                      10, 0, TMath::Pi(), 10, 0, 1000);
        h2_HTMET[i]    =   InitTH2F(Form("h2_%s_HTMET_%ifatjet", ch->GetTitle(), i),
                                    Form("h2_%s_HTMET_%ifatjet", ch->GetTitle(), i), 
                                    25, 500, 1750, 10, 250, 750);
        h2_MJmT[i]    =   InitTH2F(Form("h2_%s_MJmT_%ifatjet", ch->GetTitle(), i),
                                   Form("h2_%s_MJmT_%ifatjet", ch->GetTitle(), i), 
                                   15, 0, 750, 10, 0, 500);
        h2_HTmT[i]    =   InitTH2F(Form("h2_%s_HTmT_%ifatjet", ch->GetTitle(), i),
                                   Form("h2_%s_HTmT_%ifatjet", ch->GetTitle(), i), 
                                   25, 500, 1750, 10, 0, 500);
        h2_MJMET[i]    =   InitTH2F(Form("h2_%s_MJMET_%ifatjet", ch->GetTitle(), i),
                                   Form("h2_%s_MJMET_%ifatjet", ch->GetTitle(), i), 
                                   15, 0, 750, 10, 250, 750);
        h2_HTMJ[i]    =   InitTH2F(Form("h2_%s_HTMJ_%ifatjet", ch->GetTitle(), i),
                                   Form("h2_%s_HTMJ_%ifatjet", ch->GetTitle(), i), 
                                   25, 500, 1750, 15, 0, 750);
        h2_METmT[i]    =   InitTH2F(Form("h2_%s_METmT_%ifatjet", ch->GetTitle(), i),
                                   Form("h2_%s_METmT_%ifatjet", ch->GetTitle(), i), 
                                   10, 250, 750, 10, 0, 500);
        h2_mj1vsmj2[i]    =   InitTH2F(Form("h2_%s_mj1vsmj2_%ifatjet", ch->GetTitle(), i),
                                       Form("h2_%s_mj1vsmj2_%ifatjet", ch->GetTitle(), i), 
                                       10, 0, 500, 10, 0, 500);
        h2_mj2vsmj3[i]    =   InitTH2F(Form("h2_%s_mj2vsmj3_%ifatjet", ch->GetTitle(), i),
                                       Form("h2_%s_mj2vsmj3_%ifatjet", ch->GetTitle(), i), 
                                       10, 0, 500, 10, 0, 500);
        h2_mj3vsmj4[i]    =   InitTH2F(Form("h2_%s_mj3vsmj4_%ifatjet", ch->GetTitle(), i),
                                       Form("h2_%s_mj4vsmj4_%ifatjet", ch->GetTitle(), i), 
                                       10, 0, 500, 10, 0, 500);
        h2_mj1vspt1[i]    =   InitTH2F(Form("h2_%s_mj1vspt1_%ifatjet", ch->GetTitle(), i),
                                       Form("h2_%s_mj1vspt1_%ifatjet", ch->GetTitle(), i), 
                                       10, 0, 500, 10, 0, 2000);
        h2_mj2vspt2[i]    =   InitTH2F(Form("h2_%s_mj2vspt2_%ifatjet", ch->GetTitle(), i),
                                       Form("h2_%s_mj2vspt2_%ifatjet", ch->GetTitle(), i), 
                                       10, 0, 250, 10, 0, 1200);
        h2_mj3vspt3[i]    =   InitTH2F(Form("h2_%s_mj3vspt3_%ifatjet", ch->GetTitle(), i),
                                       Form("h2_%s_mj3vspt3_%ifatjet", ch->GetTitle(), i), 
                                       10, 0, 150, 10, 0, 800);
        h2_mj4vspt4[i]    =   InitTH2F(Form("h2_%s_mj4vspt4_%ifatjet", ch->GetTitle(), i),
                                       Form("h2_%s_mj4vspt4_%ifatjet", ch->GetTitle(), i), 
                                       10, 0, 100, 10, 0, 500);
    }
    
    //
    // Progress indicator
    //
    int i_permille_old = 0; 
    TDatime DTStart;
    int StartDate = DTStart.GetDate();
    int StartTime = DTStart.GetTime();
    cout << "[MJ Analysis] Start time : " << (StartTime/10000)%100 << ":"
        << (StartTime/100)%100 << ":" << StartTime%100
        << endl;
 
    //
    // Loop over events
    //
    Int_t nentries = (Int_t)ch->GetEntries();
    for(int i = 0; i<nentries; i++)
    {
        ch->GetEntry(i); 

        // Progress indicator begin --------------------------------
        int i_permille = (int)floor(1000 * i / float(nentries));
        TDatime DTCurrent;
        int CurrentDate = DTCurrent.GetDate();
        int CurrentTime = DTCurrent.GetTime();
        int TimeLaps = (CurrentDate-StartDate)*1000000+(CurrentTime-StartTime);
        int TimeToRun = (int)((float)nentries/(float)i)*TimeLaps;
        if (i_permille != i_permille_old)
        {
            // xterm magic from L. Vacavant and A. Cerri
            if ( isatty(1) )
            {
                printf("\015\033[32m Processed :: \033[1m\033[31m%4.1f %%"
                        "\033[0m\033[32m   Expected processing time :: \033[1m\033[31m%i:%i:%i \033[0m\015",
                        i_permille/10., (TimeToRun/10000)%100<60 ? (TimeToRun/10000)%100 : (TimeToRun/10000)%100-40,
                        (TimeToRun/100)%100<60 ? (TimeToRun/100)%100 : (TimeToRun/100)%100-40,
                        (TimeToRun%100)<60 ? (TimeToRun)%100 : (TimeToRun)%100-40 );
                fflush(stdout);
            }
            i_permille_old = i_permille;
        }
        // Progress indicator end ----------------------------------
   
        // Dividing ttbar to semi-leptonic and full-leptonic 
        int Ngenlep = ntrumus_ + ntruels_ + ntrutaush_ + ntrutausl_;
        if(ChainName.Contains("TT_sl") && Ngenlep!=1) continue;  
        if(ChainName.Contains("TT_ll") && Ngenlep!=2) continue;  

        // 
        // weights 
        // 
        EventWeight_ *= 10.; // 10 fb-1
       
        //
        // Pileup 
        //
        // blah
        
        //
        // Filters
        //
        // blah

        // 
        // single-lepton events selection
        // 
        if(!PassNLep(1))  continue; 

        float HT_thres=HT_;
        int Nskinny_thres=Nskinnyjet_;
        int Nbcsvm_thres=NBtagCSVM_;
        float mT = mt_;

        // Nfatjet, MJ, mj sorting 
        int Nfatjet_thres = 0;
        double MJ_thres=0; 
        vector<double> mj_thres_sorted; 
        vector<int> mj_thres_sorted_index; 
        for(int ifj=0; ifj<(int)FatjetPt_->size(); ifj++)
        {   
            if(FatjetPt_->at(ifj)<FatjetpTthres) continue; 
            if(mj_->at(ifj)<mjthres) continue; 
            MJ_thres = MJ_thres + mj_->at(ifj);   

            //
            // do some matching
            //
            // ------------------------- matching begins -------------------------------
            // lepton
            float dRlep=999.;
            if( RA4MusPt_->size()==1 ) dRlep = getDR( RA4MusEta_->at(0), FatjetEta_->at(ifj), 
                                                    RA4MusPhi_->at(0), FatjetPhi_->at(ifj) );
            if( RA4ElsPt_->size()==1 ) dRlep = getDR( RA4ElsEta_->at(0), FatjetEta_->at(ifj), 
                                                    RA4ElsPhi_->at(0), FatjetPhi_->at(ifj) );

            // MET 
            float dPhiMET=getDPhi(FatjetPhi_->at(ifj), METPhi_);
             
            // B-jet 
            float dRbmin=999.;
            for(unsigned int ijet=0; ijet<JetPt_->size(); ijet++) 
            { 

                if(JetPt_->at(ijet)<30 || TMath::Abs(JetEta_->at(ijet))>2.5) continue;

                if(JetCSV_->at(ijet)<0.814) continue;

                float dRb_tmp = getDR( JetEta_->at(ijet), FatjetEta_->at(ifj), 
                                       JetPhi_->at(ijet), FatjetPhi_->at(ifj) );
                if( dRb_tmp<dRbmin ) dRbmin=dRb_tmp;                         
            }
            
            // dPhi(MET,lep)
            float dPhiMETlep=999; 
            if( RA4MusPt_->size()==1 ) dPhiMETlep=getDPhi(RA4MusPhi_->at(0), METPhi_);
            if( RA4ElsPt_->size()==1 ) dPhiMETlep=getDPhi(RA4ElsPhi_->at(0), METPhi_);
           
            // 
            FillTH1FAll(h1_dRlep,       3, dRlep,   EventWeight_);  
            FillTH1FAll(h1_dPhiMET,     3, dPhiMET, EventWeight_);  
            FillTH1FAll(h1_dRbmin,      3, dRbmin,  EventWeight_);  
            FillTH1FAll(h1_dPhiMETlep,  3, dPhiMETlep, EventWeight_);  
            
            // ------------------------- matching done -------------------------------
            
            Nfatjet_thres++;
            mj_thres_sorted.push_back(mj_->at(ifj));
        }
        sort(mj_thres_sorted.begin(), mj_thres_sorted.end());
        reverse(mj_thres_sorted.begin(), mj_thres_sorted.end());

        // Get indices of sorted mj 
        // I know this is not the best way of doing this, but ... 
        for(int imj_sorted=0; imj_sorted<(int)mj_thres_sorted.size(); imj_sorted++) 
        { 
            for(int imj=0; imj<(int)mj_->size(); imj++)
            { 
                if(FatjetPt_->at(imj)<FatjetpTthres) continue; 
                if(mj_->at(imj)<mjthres) continue; 
                if(mj_->at(imj) == mj_thres_sorted.at(imj_sorted) ) 
                {   
                    mj_thres_sorted_index.push_back(imj);
                    continue;
                }
            }
        }
        if(mj_thres_sorted.size() != mj_thres_sorted_index.size() ) 
        { 
            cout << "[MJ Analysis] !! Caution : Something is wrong with sorted mj index !!" << endl;
            continue;
        }
        
        // 
        // R=0.8 for corroborators
        // 
        vector<double> mj08_thres_sorted; 
        vector<int> mj08_thres_sorted_index; 
        for(int imj08=0; imj08<(int)mj08_->size(); imj08++)
        {
            mj08_thres_sorted.push_back(mj08_->at(imj08));
        }
        sort(mj08_thres_sorted.begin(), mj08_thres_sorted.end());
        reverse(mj08_thres_sorted.begin(), mj08_thres_sorted.end());
        
        for(int imj08_sorted=0; imj08_sorted<(int)mj08_thres_sorted.size(); imj08_sorted++) 
        { 
            for(int imj08=0; imj08<(int)mj08_->size(); imj08++)
            { 
                if(mj08_->at(imj08) == mj08_thres_sorted.at(imj08_sorted) ) 
                {   
                    mj08_thres_sorted_index.push_back(imj08);
                    continue;
                }
            }
        }
        if(mj08_thres_sorted.size() != mj08_thres_sorted_index.size() ) 
        { 
            cout << "[MJ Analysis] !! Caution : Something is wrong with sorted mj08 index !!" << endl;
            continue;
        }


        //
        // Apply selection   
        //
        // baseline selection
        if( !PassBaselineSelection(HT_thres, MET_, Nbcsvm_thres, Nskinny_thres) ) continue; 
        if( !PassSelection(Region, HT_thres, MET_, Nbcsvm_thres, Nskinny_thres, mT, MJ_thres)) continue;
        
        //
        // NFJbin
        //
        int NFJbin = -1;
        if(Nfatjet_thres>4) NFJbin=5;
        else if(Nfatjet_thres>3) NFJbin=4;
        else if(Nfatjet_thres>2) NFJbin=3;
        else if(Nfatjet_thres>1) NFJbin=2;
        else if(Nfatjet_thres>0) NFJbin=1;
        else NFJbin=0;

        if(NFJbin==-1) 
        {   
            cout << "[MJ Analysis] ERROR : NFJ is cannot be negative" << endl;
            continue;
        }
        
/*      
        // gen top info
        float top1Phi=-999;  
        float top2Phi=-999; 
        float top1pT=-999;  
        float top2pT=-999; 
        for(unsigned int igen=0; igen<GenId_->size(); igen++)
        { 
            if( GenId_->at(igen)==6 && GenMId_->at(igen)!=6)
            {       
                top1Phi = GenPhi_->at(igen);  
                top1pT  = GenPt_->at(igen);  
            }
            if( GenId_->at(igen)==-6 && GenMId_->at(igen)!=-6)
            {
                top2Phi = GenPhi_->at(igen);  
                top2pT  = GenPt_->at(igen);  
            }

        }
        float ISRweight=1.;
        if(ChainName.Contains("TT")) 
        {   
            FillTH1FAll(h1_dPhitop,     NFJbin, getDPhi(top1Phi, top2Phi), EventWeight_);  
            FillTH2FAll(h2_dPhidpTtop,  NFJbin, getDPhi(top1Phi, top2Phi), TMath::Abs(top1pT-top2pT), EventWeight_); 
            float ISRpx = top1pT*TMath::Cos(top1Phi) + top2pT*TMath::Cos(top2Phi);
            float ISRpy = top1pT*TMath::Sin(top1Phi) + top2pT*TMath::Sin(top2Phi);
            float ISRpT = TMath::Sqrt(ISRpx*ISRpx+ISRpy*ISRpy);
            ISRweight = getISRSF(ISRpT);
            FillTH2FAll(h2_dPhidpTttbar,  NFJbin, getDPhi(top1Phi, top2Phi), ISRpT, EventWeight_); 
        }
*/

        //
        // Fill histogams 
        //

        // yields
	    if(nels_==1) FillTH1FAll(h1_yields, NFJbin, 0.5, EventWeight_);   
	    if(nmus_==1) FillTH1FAll(h1_yields, NFJbin, 1.5, EventWeight_);

        if(NFJbin>1) FillTH2FAll(h2_mj1vsmj2, NFJbin, mj_->at(mj_thres_sorted_index.at(0)), mj_->at(mj_thres_sorted_index.at(1)), EventWeight_);           
        if(NFJbin>2) FillTH2FAll(h2_mj2vsmj3, NFJbin, mj_->at(mj_thres_sorted_index.at(1)), mj_->at(mj_thres_sorted_index.at(2)), EventWeight_);           
        if(NFJbin>3) FillTH2FAll(h2_mj3vsmj4, NFJbin, mj_->at(mj_thres_sorted_index.at(2)), mj_->at(mj_thres_sorted_index.at(3)), EventWeight_);           
        FillTH2FAll(h2_HTMET, NFJbin, HT_thres, MET_, EventWeight_);         
        FillTH2FAll(h2_MJmT, NFJbin, MJ_thres, mT, EventWeight_);       
        FillTH2FAll(h2_HTmT, NFJbin, HT_thres, mT, EventWeight_);            
        FillTH2FAll(h2_MJMET, NFJbin, MJ_thres, MET_, EventWeight_);    
        FillTH2FAll(h2_HTMJ, NFJbin, HT_thres, MJ_thres,EventWeight_);       
        FillTH2FAll(h2_METmT, NFJbin, MET_, mT, EventWeight_);           
        FillTH1FAll(h1_mT,  NFJbin, mT, EventWeight_);                 
        FillTH1FAll(h1_HT, NFJbin, HT_thres, EventWeight_);                
        FillTH1FAll(h1_MJ, NFJbin, MJ_thres, EventWeight_); 
        FillTH1FAll(h1_MET, NFJbin, MET_, EventWeight_);              
        FillTH1FAll(h1_METPhi, NFJbin, METPhi_, EventWeight_);        
        FillTH1FAll(h1_METx, NFJbin, MET_*TMath::Cos(METPhi_), EventWeight_);           
        FillTH1FAll(h1_METy, NFJbin, MET_*TMath::Sin(METPhi_), EventWeight_);           
        
        if(Nfatjet_thres>0) 
        {
            FillTH1FAll(h1_FatjetPt1,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(0)),       EventWeight_);   
            FillTH1FAll(h1_mj1,         NFJbin, mj_->at(mj_thres_sorted_index.at(0)),  EventWeight_);    
            FillTH1FAll(h1_mj1OverMJ,   NFJbin, mj_->at(mj_thres_sorted_index.at(0))/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt1,   NFJbin, mj_->at(mj_thres_sorted_index.at(0))/FatjetPt_->at(mj_thres_sorted_index.at(0)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi1,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(0)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta1,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(0)),      EventWeight_); 
            FillTH2FAll(h2_mj1vspt1, NFJbin, mj_->at(mj_thres_sorted_index.at(0)), FatjetPt_->at(mj_thres_sorted_index.at(0)), EventWeight_);           
        }
        if(Nfatjet_thres>1) 
        {
            FillTH1FAll(h1_FatjetPt2,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(1)),       EventWeight_);   
            FillTH1FAll(h1_mj2,         NFJbin, mj_->at(mj_thres_sorted_index.at(1)),  EventWeight_);    
            FillTH1FAll(h1_mj2OverMJ,   NFJbin, mj_->at(mj_thres_sorted_index.at(1))/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt2,   NFJbin, mj_->at(mj_thres_sorted_index.at(1))/FatjetPt_->at(mj_thres_sorted_index.at(1)),  EventWeight_);    
            FillTH1FAll(h1_mj2overmj1,  NFJbin, mj_->at(mj_thres_sorted_index.at(1))/mj_->at(mj_thres_sorted_index.at(0)), EventWeight_);   
            FillTH1FAll(h1_FatjetPhi2,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(1)),      EventWeight_);  
            FillTH1FAll(h1_FatjetEta2,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(1)),      EventWeight_);  
            FillTH2FAll(h2_mj2vspt2, NFJbin, mj_->at(mj_thres_sorted_index.at(1)), FatjetPt_->at(mj_thres_sorted_index.at(1)), EventWeight_);           
        }
        if(Nfatjet_thres>2) 
        {
            FillTH1FAll(h1_FatjetPt3,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(2)),       EventWeight_);   
            FillTH1FAll(h1_mj3,         NFJbin, mj_->at(mj_thres_sorted_index.at(2)),  EventWeight_);    
            FillTH1FAll(h1_mj3OverMJ,   NFJbin, mj_->at(mj_thres_sorted_index.at(2))/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt3,   NFJbin, mj_->at(mj_thres_sorted_index.at(2))/FatjetPt_->at(mj_thres_sorted_index.at(2)),  EventWeight_);    
            FillTH1FAll(h1_mj3overmj2,  NFJbin, mj_->at(mj_thres_sorted_index.at(2))/mj_->at(mj_thres_sorted_index.at(1)), EventWeight_);   
            FillTH1FAll(h1_FatjetPhi3,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(2)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta3,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(2)),      EventWeight_); 
            FillTH2FAll(h2_mj3vspt3, NFJbin, mj_->at(mj_thres_sorted_index.at(2)), FatjetPt_->at(mj_thres_sorted_index.at(2)), EventWeight_);           
        }
        if(Nfatjet_thres>3) 
        {
            FillTH1FAll(h1_FatjetPt4,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(3)),       EventWeight_);   
            FillTH1FAll(h1_mj4,         NFJbin, mj_->at(mj_thres_sorted_index.at(3)),  EventWeight_);    
            FillTH1FAll(h1_mj4OverMJ,   NFJbin, mj_->at(mj_thres_sorted_index.at(3))/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt4,   NFJbin, mj_->at(mj_thres_sorted_index.at(3))/FatjetPt_->at(mj_thres_sorted_index.at(3)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi4,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(3)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta4,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(3)),      EventWeight_); 
            FillTH2FAll(h2_mj4vspt4, NFJbin, mj_->at(mj_thres_sorted_index.at(3)), FatjetPt_->at(mj_thres_sorted_index.at(3)), EventWeight_);           
        }
        if(Nfatjet_thres>4) 
        {
            // add 5th jet to h1_FatjetPt4
            FillTH1FAll(h1_FatjetPt4,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(4)),       EventWeight_);   
            FillTH1FAll(h1_mj4,         NFJbin, mj_->at(mj_thres_sorted_index.at(4)),  EventWeight_);    
            FillTH1FAll(h1_mj4OverMJ,   NFJbin, mj_->at(mj_thres_sorted_index.at(4))/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt4,   NFJbin, mj_->at(mj_thres_sorted_index.at(4))/FatjetPt_->at(mj_thres_sorted_index.at(4)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi4,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(4)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta4,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(4)),      EventWeight_); 
            FillTH2FAll(h2_mj4vspt4, NFJbin, mj_->at(mj_thres_sorted_index.at(4)), FatjetPt_->at(mj_thres_sorted_index.at(4)), EventWeight_);           
        }
        if(Nfatjet_thres>5) 
        {
            // add 6th jet to h1_FatjetPt4
            FillTH1FAll(h1_FatjetPt4,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(5)),       EventWeight_);   
            FillTH1FAll(h1_mj4,         NFJbin, mj_->at(mj_thres_sorted_index.at(5)),  EventWeight_);    
            FillTH1FAll(h1_mj4OverMJ,   NFJbin, mj_->at(mj_thres_sorted_index.at(5))/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt4,   NFJbin, mj_->at(mj_thres_sorted_index.at(5))/FatjetPt_->at(mj_thres_sorted_index.at(5)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi4,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(5)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta4,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(5)),      EventWeight_); 
            FillTH2FAll(h2_mj4vspt4, NFJbin, mj_->at(mj_thres_sorted_index.at(5)), FatjetPt_->at(mj_thres_sorted_index.at(5)), EventWeight_);           
        }

        for(int imj=0; imj<(int)mj_->size(); imj++) 
        { 
            if(FatjetPt_->at(imj)<FatjetpTthres) continue; 
            if(mj_->at(imj)<mjthres) continue; 
            
            FillTH1FAll(h1_mj,          NFJbin, mj_->at(imj),           EventWeight_);                            
            FillTH1FAll(h1_FatjetPt,    NFJbin, FatjetPt_->at(imj),     EventWeight_);                
            FillTH1FAll(h1_FatjetEta,   NFJbin, FatjetEta_->at(imj),    EventWeight_);              
        }

        FillTH1FAll(h1_Nfatjet,     NFJbin, Nfatjet_thres,  EventWeight_);       
        FillTH1FAll(h1_Nskinnyjet,  NFJbin, Nskinny_thres,    EventWeight_);      
        FillTH1FAll(h1_Ncsvm,       NFJbin, Nbcsvm_thres,     EventWeight_);            
      
        //
        // corroborators
        //

        // mJ R=0.8 
        float MJ08=0;
        for(int imj08=0; imj08<(int)mj08_->size(); imj08++) MJ08 += mj08_->at(imj08);
        if(MJ08>0)
        {
            if(mj08_->size()>0) FillTH1FAll(h1_mj08_1,         NFJbin, mj08_->at(mj08_thres_sorted_index.at(0)),  EventWeight_);
            if(mj08_->size()>1) FillTH1FAll(h1_mj08_2,         NFJbin, mj08_->at(mj08_thres_sorted_index.at(1)),  EventWeight_);
            if(mj08_->size()>2) FillTH1FAll(h1_mj08_3,         NFJbin, mj08_->at(mj08_thres_sorted_index.at(2)),  EventWeight_);
            if(mj08_->size()>3) FillTH1FAll(h1_mj08_4,         NFJbin, mj08_->at(mj08_thres_sorted_index.at(3)),  EventWeight_);
        }

        // min dPhi bb 
        float   mindPhibb=999.;
        for(unsigned int ijet=0; ijet<JetPt_->size(); ijet++)
        {
            if( JetCSV_->at(ijet)<0.814 ) continue;
            for(unsigned int jjet=0; jjet<JetPt_->size(); jjet++) 
            { 
                if(ijet==jjet) continue;
                if( JetCSV_->at(jjet)<0.814 ) continue;
                float dPhitemp = getDPhi(JetPhi_->at(ijet), JetPhi_->at(jjet));
                if( dPhitemp<mindPhibb ) mindPhibb = dPhitemp;
            }
        }
        FillTH1FAll(h1_mindPhibb,         NFJbin, mindPhibb,  EventWeight_);

    } // for(int i = 0; i<nentries; i++)
   
    TString HistFileName = ch->GetTitle();
    HistFileName = Form("HistFiles/%s_%s.root", HistFileName.Data(), Region);
    cout << "[MJ Analysis] Writing " << HistFileName << endl;
    TFile *HistFile = new TFile(HistFileName, "RECREATE");
    gROOT->cd();
    HistFile->cd();
    // write histograms
    for(int i=0; i<7; i++)  
    {
        h1_yields[i]->SetDirectory(0);                      h1_yields[i]->Write();
        h1_dPhitop[i]->SetDirectory(0);                     h1_dPhitop[i]->Write();
        h1_dRlep[i]->SetDirectory(0);                       h1_dRlep[i]->Write();
        h1_dRbmin[i]->SetDirectory(0);                      h1_dRbmin[i]->Write();
        h1_dPhiMET[i]->SetDirectory(0);                     h1_dPhiMET[i]->Write();
        h1_dPhiMETlep[i]->SetDirectory(0);                  h1_dPhiMETlep[i]->Write();
        h2_dPhidpTtop[i]->SetDirectory(0);                  h2_dPhidpTtop[i]->Write();
        h2_dPhidpTttbar[i]->SetDirectory(0);                h2_dPhidpTttbar[i]->Write();
        h1_MJ[i]->SetDirectory(0);                          h1_MJ[i]->Write();
        h1_MJ_ISR[i]->SetDirectory(0);                      h1_MJ_ISR[i]->Write();
        h1_HT[i]->SetDirectory(0);                          h1_HT[i]->Write();
        h1_MET[i]->SetDirectory(0);                         h1_MET[i]->Write();
        h1_METPhi[i]->SetDirectory(0);                      h1_METPhi[i]->Write();
        h1_METx[i]->SetDirectory(0);                        h1_METx[i]->Write();
        h1_METy[i]->SetDirectory(0);                        h1_METy[i]->Write();
        h1_FatjetPt[i]->SetDirectory(0);                    h1_FatjetPt[i]->Write();
        h1_FatjetPt1[i]->SetDirectory(0);                   h1_FatjetPt1[i]->Write();
        h1_FatjetPt2[i]->SetDirectory(0);                   h1_FatjetPt2[i]->Write();
        h1_FatjetPt3[i]->SetDirectory(0);                   h1_FatjetPt3[i]->Write();
        h1_FatjetPt4[i]->SetDirectory(0);                   h1_FatjetPt4[i]->Write();
        h1_mj1[i]->SetDirectory(0);                         h1_mj1[i]->Write();
        h1_mj2[i]->SetDirectory(0);                         h1_mj2[i]->Write();
        h1_mj3[i]->SetDirectory(0);                         h1_mj3[i]->Write();
        h1_mj4[i]->SetDirectory(0);                         h1_mj4[i]->Write();
        h1_mj08_1[i]->SetDirectory(0);                      h1_mj08_1[i]->Write();
        h1_mj08_2[i]->SetDirectory(0);                      h1_mj08_2[i]->Write();
        h1_mj08_3[i]->SetDirectory(0);                      h1_mj08_3[i]->Write();
        h1_mj08_4[i]->SetDirectory(0);                      h1_mj08_4[i]->Write();
        h1_mj1OverMJ[i]->SetDirectory(0);                   h1_mj1OverMJ[i]->Write();
        h1_mj2OverMJ[i]->SetDirectory(0);                   h1_mj2OverMJ[i]->Write();
        h1_mj3OverMJ[i]->SetDirectory(0);                   h1_mj3OverMJ[i]->Write();
        h1_mj4OverMJ[i]->SetDirectory(0);                   h1_mj4OverMJ[i]->Write();
        h1_N1[i]->SetDirectory(0);                          h1_N1[i]->Write();
        h1_N1[i]->SetDirectory(0);                          h1_N1[i]->Write();
        h1_N2[i]->SetDirectory(0);                          h1_N2[i]->Write();
        h1_N3[i]->SetDirectory(0);                          h1_N3[i]->Write();
        h1_N4[i]->SetDirectory(0);                          h1_N4[i]->Write();
        h1_mjOverPt1[i]->SetDirectory(0);                   h1_mjOverPt1[i]->Write();
        h1_mjOverPt2[i]->SetDirectory(0);                   h1_mjOverPt2[i]->Write();
        h1_mjOverPt3[i]->SetDirectory(0);                   h1_mjOverPt3[i]->Write();
        h1_mjOverPt4[i]->SetDirectory(0);                   h1_mjOverPt4[i]->Write();
        h1_mj3overmj2[i]->SetDirectory(0);                  h1_mj3overmj2[i]->Write();
        h1_mj2overmj1[i]->SetDirectory(0);                  h1_mj2overmj1[i]->Write();
        h1_FatjetPhi1[i]->SetDirectory(0);                  h1_FatjetPhi1[i]->Write();
        h1_FatjetPhi2[i]->SetDirectory(0);                  h1_FatjetPhi2[i]->Write();
        h1_FatjetPhi3[i]->SetDirectory(0);                  h1_FatjetPhi3[i]->Write();
        h1_FatjetPhi4[i]->SetDirectory(0);                  h1_FatjetPhi4[i]->Write();
        h1_FatjetEta1[i]->SetDirectory(0);                  h1_FatjetEta1[i]->Write();
        h1_FatjetEta2[i]->SetDirectory(0);                  h1_FatjetEta2[i]->Write();
        h1_FatjetEta3[i]->SetDirectory(0);                  h1_FatjetEta3[i]->Write();
        h1_FatjetEta4[i]->SetDirectory(0);                  h1_FatjetEta4[i]->Write();
        h1_dRFJ[i]->SetDirectory(0);                        h1_dRFJ[i]->Write();
        h1_dPhiFJ[i]->SetDirectory(0);                      h1_dPhiFJ[i]->Write();
        h1_dEtaFJ[i]->SetDirectory(0);                      h1_dEtaFJ[i]->Write();
        h1_mindPhibb[i]->SetDirectory(0);                   h1_mindPhibb[i]->Write();
        h2_HTMET[i]->SetDirectory(0);                       h2_HTMET[i]->Write();
        h2_MJmT[i]->SetDirectory(0);                        h2_MJmT[i]->Write();
        h2_HTmT[i]->SetDirectory(0);                        h2_HTmT[i]->Write();
        h2_MJMET[i]->SetDirectory(0);                       h2_MJMET[i]->Write();
        h2_HTMJ[i]->SetDirectory(0);                        h2_HTMJ[i]->Write();
        h2_METmT[i]->SetDirectory(0);                       h2_METmT[i]->Write();
        h1_FatjetEta[i]->SetDirectory(0);                   h1_FatjetEta[i]->Write();
        h1_DPhi[i]->SetDirectory(0);                        h1_DPhi[i]->Write();
        h1_mj[i]->SetDirectory(0);                          h1_mj[i]->Write();
        h1_mT[i]->SetDirectory(0);                          h1_mT[i]->Write();
        h1_WpT[i]->SetDirectory(0);                         h1_WpT[i]->Write();
        h1_muspT[i]->SetDirectory(0);                       h1_muspT[i]->Write();
        h1_muspTminusMET[i]->SetDirectory(0);               h1_muspTminusMET[i]->Write();
        h1_musEta[i]->SetDirectory(0);                      h1_musEta[i]->Write();
        h1_musPhi[i]->SetDirectory(0);                      h1_musPhi[i]->Write();
        h1_elspT[i]->SetDirectory(0);                       h1_elspT[i]->Write();
        h1_elspTminusMET[i]->SetDirectory(0);               h1_elspTminusMET[i]->Write();
        h1_elsEta[i]->SetDirectory(0);                      h1_elsEta[i]->Write();
        h1_elsPhi[i]->SetDirectory(0);                      h1_elsPhi[i]->Write();
        h1_Nfatjet[i]->SetDirectory(0);                     h1_Nfatjet[i]->Write();
        h1_Nskinnyjet[i]->SetDirectory(0);                  h1_Nskinnyjet[i]->Write();
        h1_Ncsvm[i]->SetDirectory(0);                       h1_Ncsvm[i]->Write();
        h2_mj1vsmj2[i]->SetDirectory(0);                    h2_mj1vsmj2[i]->Write();
        h2_mj2vsmj3[i]->SetDirectory(0);                    h2_mj2vsmj3[i]->Write();
        h2_mj3vsmj4[i]->SetDirectory(0);                    h2_mj3vsmj4[i]->Write();
        h2_mj1vspt1[i]->SetDirectory(0);                    h2_mj1vspt1[i]->Write();
        h2_mj2vspt2[i]->SetDirectory(0);                    h2_mj2vspt2[i]->Write();
        h2_mj3vspt3[i]->SetDirectory(0);                    h2_mj3vspt3[i]->Write();
        h2_mj4vspt4[i]->SetDirectory(0);                    h2_mj4vspt4[i]->Write();
    }
    HistFile->Close();
}
