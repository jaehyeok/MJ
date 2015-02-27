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

#include "babytree.h"
#include "PassSelection.h"

using namespace std;
bool    DoLog           = 1;
int     FatjetpTthres   = 50;
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


float GetMuonEff(float pt, float eta)
{ 
   float eff = 1; 
 
   // For SingleMu dataset
   if(pt<25) eff = 0;
/* 
   
   if(pt>25 && pt<30)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.9;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.81;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.78;
   }
   
   if(pt>30 && pt<35)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.91;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.83;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.80;
   }
   
   if(pt>35 && pt<40)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.92;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.84;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.80;
   }
   
   if(pt>40 && pt<50)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.94;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.85;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.82;
   }
   
   if(pt>50)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.94;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.86;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.82;
   }

   if(TMath::Abs(eta)>2.1) eff = 0; 

   eff = 0.88;
*/

    if(TMath::Abs(eta)<0.9) eff = 0.98;
    else eff = 0.84;
   
    return eff;
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
// per process
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
    TH1F *h1_MJ[7], *h1_mT[7], *h1_Nskinnyjet[7], *h1_Ncsvm[7], 
         *h1_muspT[7], *h1_muspTminusMET[7], *h1_musEta[7], *h1_musPhi[7], 
         *h1_elspT[7], *h1_elspTminusMET[7], *h1_elsEta[7], *h1_elsPhi[7], 
         *h1_mj[7], *h1_FatjetPt[7], *h1_FatjetEta[7], 
         *h1_FatjetPt1[7],  *h1_FatjetPt2[7],  *h1_FatjetPt3[7],  *h1_FatjetPt4[7],
         *h1_mj1[7],  *h1_mj2[7],  *h1_mj3[7],  *h1_mj4[7], 
         *h1_mj1OverMJ[7],  *h1_mj2OverMJ[7],  *h1_mj3OverMJ[7],  *h1_mj4OverMJ[7], 
         *h1_N1[7],  *h1_N2[7],  *h1_N3[7],  *h1_N4[7], 
         *h1_mjOverPt1[7],  *h1_mjOverPt2[7],  *h1_mjOverPt3[7],  *h1_mjOverPt4[7], 
         *h1_mj3overmj2[7], *h1_mj2overmj1[7],
         *h1_FatjetPhi1[7], *h1_FatjetPhi2[7], *h1_FatjetPhi3[7], *h1_FatjetPhi4[7],
         *h1_FatjetEta1[7], *h1_FatjetEta2[7], *h1_FatjetEta3[7], *h1_FatjetEta4[7],
         *h1_dRFJ[7], *h1_dPhiFJ[7], *h1_dEtaFJ[7],
         *h1_HT[7], *h1_MET[7], *h1_METPhi[7], *h1_METx[7], *h1_METy[7], *h1_DPhi[7], *h1_Nfatjet[7], *h1_WpT[7];
    TH2F *h2_mj1vsmj2[7], *h2_mj2vsmj3[7], *h2_mj3vsmj4[7];
    TH2F *h2_HTMET[7], *h2_MJmT[7];
    TH2F *h2_HTmT[7], *h2_MJMET[7];
    TH2F *h2_HTMJ[7], *h2_METmT[7];

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
        h1_MJ[i] = InitTH1F( Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
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
                                    20, 0, 150);
        h1_mj4[i] = InitTH1F( Form("h1_%s_mj4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj4_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 100);
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
        //
        // h2                    
        //
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
    }
    
    //
    //
    //
    int i_permille_old = 0; 
    TDatime DTStart;
    int StartDate = DTStart.GetDate();
    int StartTime = DTStart.GetTime();
    cout << "[MJ Analysis] Start time : " << (StartTime/10000)%100 << ":"
        << (StartTime/100)%100 << ":" << StartTime%100
        << endl;
  
    // 
    // TT composition 
    // 
    float Nll=0.;
    float Nlt=0.;
    float Ntt=0.;
    float Nltl=0.;
    float Nlth=0.;
    //
    //
    //
    //InitTree(ch);
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
    
        // separate TT_sl and TT_ll
        int Ngenlep=0;
        for(unsigned int igen=0; igen<GenId_->size(); igen++)
        { 
            if( TMath::Abs(GenMId_->at(igen))!=24 ) continue;
            if( TMath::Abs(GenId_->at(igen))!=11 && TMath::Abs(GenId_->at(igen))!=13 && TMath::Abs(GenId_->at(igen))!=15 ) continue;
            Ngenlep++;
        }
        if(ChainName.Contains("TT_sl") && Ngenlep!=1) continue;  
        if(ChainName.Contains("TT_ll") && Ngenlep!=2) continue;  

        // 
        // weights 
        // 
        // Temp fixes for wrong event weights 
        // TTJets cross section  
        //if(ChainName.Contains("TT"))  
        //{ 
            //EventWeight_ = EventWeight_*806.1/832.; // I used 832 pb while M used 806.1 pb.
        //}
        // Pileup 

        // 
        // single-lepton events selection
        // 
        if(!PassNLep(1))  continue; // need this upfront because of mT calculation

        // Nfatjet, MJ, mj sorting 
        int Nfatjet_thres = 0;
        double MJ_thres=0; 
        vector<double> mj_thres_sorted; 
        vector<int> mj_thres_sorted_index; 
        for(int ifj=0; ifj<(int)FatjetPt_->size(); ifj++)
        {   
            if(FatjetPt_->at(ifj)<FatjetpTthres) continue; 
            if(mj_->at(ifj)<mjthres) continue; 
            Nfatjet_thres++;
            MJ_thres = MJ_thres + mj_->at(ifj);   
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
        // DEBUG  
        //cout << "_____________________________" << endl;
        for(int imj_sorted_index=0; imj_sorted_index<(int)mj_thres_sorted_index.size(); imj_sorted_index++) 
        { 
            //cout << imj_sorted_index << " : " 
            //     << mj_thres_sorted.at(imj_sorted_index) << " " 
            //     << mj_->at(mj_thres_sorted_index.at(imj_sorted_index)) << " "
            //     << mj_->at(imj_sorted_index) << endl;
            if( mj_thres_sorted.at(imj_sorted_index) != mj_->at(mj_thres_sorted_index.at(imj_sorted_index))) 
                    cout << " Wrong index matching" << endl;
        }

        //
        // Calculate variables 
        //
        double mT;
        double WpT;
        if(RA4MusPt_->size()==1) {
            mT  = TMath::Sqrt( 2*MET_*RA4MusPt_->at(0)*(1-TMath::Cos(METPhi_-RA4MusPhi_->at(0))) ); 
            WpT = TMath::Sqrt(  
                        (RA4MusPt_->at(0)*TMath::Cos(RA4MusPhi_->at(0)) + MET_*TMath::Cos(METPhi_))
                       *(RA4MusPt_->at(0)*TMath::Cos(RA4MusPhi_->at(0)) + MET_*TMath::Cos(METPhi_))
                       +(RA4MusPt_->at(0)*TMath::Sin(RA4MusPhi_->at(0)) + MET_*TMath::Sin(METPhi_))
                       *(RA4MusPt_->at(0)*TMath::Sin(RA4MusPhi_->at(0)) + MET_*TMath::Sin(METPhi_))  ); 
        }
        if(RA4ElsPt_->size()==1) {
            mT  = TMath::Sqrt( 2*MET_*RA4ElsPt_->at(0)*(1-TMath::Cos(METPhi_-RA4ElsPhi_->at(0))) ); 
            WpT = TMath::Sqrt(  
                        (RA4ElsPt_->at(0)*TMath::Cos(RA4ElsPhi_->at(0)) + MET_*TMath::Cos(METPhi_))
                       *(RA4ElsPt_->at(0)*TMath::Cos(RA4ElsPhi_->at(0)) + MET_*TMath::Cos(METPhi_))
                       +(RA4ElsPt_->at(0)*TMath::Sin(RA4ElsPhi_->at(0)) + MET_*TMath::Sin(METPhi_))
                       *(RA4ElsPt_->at(0)*TMath::Sin(RA4ElsPhi_->at(0)) + MET_*TMath::Sin(METPhi_))  ); 
        }

        //
        // Apply selection   
        //
        // baseline selection
        if( !PassBaselineSelection() ) continue; 
        if( !PassSelection(Region, HT_, MET_, NBtagCSVM_, Nskinnyjet_, mT, MJ_thres)) continue;
 
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
        
        //     
        // ttbar gen composition  
        //     
        int Ne=0;
        int Nm=0;
        int Nt=0;
        int Ntn=0; int Nen=0; int Nmn=0;
        for(unsigned int igen=0; igen<GenId_->size(); igen++)
        { 
            if( TMath::Abs(GenMId_->at(igen))==24 && TMath::Abs(GenId_->at(igen))==11 ) Ne++;
            if( TMath::Abs(GenMId_->at(igen))==24 && TMath::Abs(GenId_->at(igen))==13 ) Nm++;
            if( TMath::Abs(GenMId_->at(igen))==24 && TMath::Abs(GenId_->at(igen))==15 ) Nt++;
            if( (TMath::Abs(GenGMId_->at(igen))==24||TMath::Abs(GenGMId_->at(igen))==15) && TMath::Abs(GenMId_->at(igen))==15 && TMath::Abs(GenId_->at(igen))==12 ) Nen++;
            if( (TMath::Abs(GenGMId_->at(igen))==24||TMath::Abs(GenGMId_->at(igen))==15) && TMath::Abs(GenMId_->at(igen))==15 && TMath::Abs(GenId_->at(igen))==14 ) Nmn++;
            if( (TMath::Abs(GenGMId_->at(igen))==24||TMath::Abs(GenGMId_->at(igen))==15) && TMath::Abs(GenMId_->at(igen))==15 && TMath::Abs(GenId_->at(igen))==16 ) Ntn++;
        }
        //cout << Ne << " " << Nm << " " << Nt << " " << Nen << " " << Nmn << " " << Ntn << endl;
        if( Ne==2 || (Ne==1 && Nm==1) || Nm==2 )                                    Nll  = Nll + EventWeight_; 
        if( (Ne==1 || Nm==1) && Nt==1 )                                             Nlt  = Nlt + EventWeight_; 
        if( Nt==2 )                                                                 Ntt  = Ntt + EventWeight_; 
        if( (Ne==1 || Nm==1) && Nt==1 && (Ntn==1 && (Nen==1 || Nmn==1)) )           Nltl = Nltl + EventWeight_; 
        if( (Ne==1 || Nm==1) && Nt==1 && (Ntn==1 && (Nen==0 && Nmn==0)) )           Nlth = Nlth + EventWeight_; 

        //
        // Fill histogams 
        //

        // yields
        if(RA4ElsPt_->size()==1) FillTH1FAll(h1_yields, NFJbin, 0.5, EventWeight_);   
        if(RA4MusPt_->size()==1) FillTH1FAll(h1_yields, NFJbin, 1.5, EventWeight_);  
        
        // plots
        if(RA4MusPt_->size()==1) {
            FillTH1FAll(h1_muspT,   NFJbin,  RA4MusPt_->at(0),  EventWeight_); 
            FillTH1FAll(h1_musEta,  NFJbin, RA4MusEta_->at(0), EventWeight_); 
            FillTH1FAll(h1_musPhi,  NFJbin, RA4MusPhi_->at(0), EventWeight_); 
            if(getDPhi(RA4MusPhi_->at(0),METPhi_)<0.4) 
            { 
                FillTH1FAll(h1_muspTminusMET, NFJbin, (MET_-RA4MusPt_->at(0))/RA4MusPt_->at(0), EventWeight_);   
            }
        }
        if(RA4ElsPt_->size()==1) {
            FillTH1FAll(h1_elspT, NFJbin, RA4ElsPt_->at(0), EventWeight_);   
            FillTH1FAll(h1_elsEta, NFJbin, RA4ElsEta_->at(0), EventWeight_); 
            FillTH1FAll(h1_elsPhi, NFJbin, RA4ElsPhi_->at(0), EventWeight_); 
            if(getDPhi(RA4ElsPhi_->at(0),METPhi_)<0.4) 
            { 
                FillTH1FAll(h1_elspTminusMET, NFJbin, (MET_-RA4ElsPt_->at(0))/RA4ElsPt_->at(0), EventWeight_);   
            }
        }

        if(NFJbin>1) FillTH2FAll(h2_mj1vsmj2, NFJbin, mj_thres_sorted.at(0), mj_thres_sorted.at(1), EventWeight_);           
        if(NFJbin>2) FillTH2FAll(h2_mj2vsmj3, NFJbin, mj_thres_sorted.at(1), mj_thres_sorted.at(2), EventWeight_);           
        if(NFJbin>3) FillTH2FAll(h2_mj3vsmj4, NFJbin, mj_thres_sorted.at(2), mj_thres_sorted.at(3), EventWeight_);           
        FillTH2FAll(h2_HTMET, NFJbin, HT_, MET_, EventWeight_);         
        FillTH2FAll(h2_MJmT, NFJbin, MJ_thres, mT, EventWeight_);       
        FillTH2FAll(h2_HTmT, NFJbin, HT_, mT, EventWeight_);            
        FillTH2FAll(h2_MJMET, NFJbin, MJ_thres, MET_, EventWeight_);    
        FillTH2FAll(h2_HTMJ, NFJbin, HT_, MJ_thres,EventWeight_);       
        FillTH2FAll(h2_METmT, NFJbin, MET_, mT, EventWeight_);           
        FillTH1FAll(h1_mT,  NFJbin, mT, EventWeight_);                 
        FillTH1FAll(h1_WpT, NFJbin, WpT, EventWeight_);               
        FillTH1FAll(h1_HT, NFJbin, HT_, EventWeight_);                
        FillTH1FAll(h1_MJ, NFJbin, MJ_thres, EventWeight_); 
        FillTH1FAll(h1_MET, NFJbin, MET_, EventWeight_);              
        FillTH1FAll(h1_METPhi, NFJbin, METPhi_, EventWeight_);        
        FillTH1FAll(h1_METx, NFJbin, MET_*TMath::Cos(METPhi_), EventWeight_);           
        FillTH1FAll(h1_METy, NFJbin, MET_*TMath::Sin(METPhi_), EventWeight_);           
        if(RA4MusPt_->size()==1) FillTH1FAll(h1_DPhi, NFJbin, getDPhi(RA4MusPhi_->at(0),METPhi_), EventWeight_);   
        if(RA4ElsPt_->size()==1) FillTH1FAll(h1_DPhi, NFJbin, getDPhi(RA4ElsPhi_->at(0),METPhi_), EventWeight_);   
        
        if(Nfatjet_thres>0) 
        {
            FillTH1FAll(h1_FatjetPt1,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(0)),       EventWeight_);   
            FillTH1FAll(h1_mj1,         NFJbin, mj_thres_sorted.at(0),  EventWeight_);    
            FillTH1FAll(h1_mj1OverMJ,   NFJbin, mj_thres_sorted.at(0)/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt1,   NFJbin, mj_thres_sorted.at(0)/FatjetPt_->at(mj_thres_sorted_index.at(0)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi1,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(0)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta1,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(0)),      EventWeight_); 
//            FillTH1FAll(h1_N1,          NFJbin, FatjetN_->at(0),        EventWeight_); 
        }
        if(Nfatjet_thres>1) 
        {
            FillTH1FAll(h1_FatjetPt2,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(1)),       EventWeight_);   
            FillTH1FAll(h1_mj2,         NFJbin, mj_thres_sorted.at(1),  EventWeight_);    
            FillTH1FAll(h1_mj2OverMJ,   NFJbin, mj_thres_sorted.at(1)/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt2,   NFJbin, mj_thres_sorted.at(1)/FatjetPt_->at(mj_thres_sorted_index.at(1)),  EventWeight_);    
            FillTH1FAll(h1_mj2overmj1,  NFJbin, mj_thres_sorted.at(1)/mj_thres_sorted.at(0), EventWeight_);   
            FillTH1FAll(h1_FatjetPhi2,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(1)),      EventWeight_);  
            FillTH1FAll(h1_FatjetEta2,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(1)),      EventWeight_);  
//            FillTH1FAll(h1_N2,          NFJbin, FatjetN_->at(1),      EventWeight_); 
        }
        if(Nfatjet_thres>2) 
        {
            FillTH1FAll(h1_FatjetPt3,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(2)),       EventWeight_);   
            FillTH1FAll(h1_mj3,         NFJbin, mj_thres_sorted.at(2),  EventWeight_);    
            FillTH1FAll(h1_mj3OverMJ,   NFJbin, mj_thres_sorted.at(2)/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt3,   NFJbin, mj_thres_sorted.at(2)/FatjetPt_->at(mj_thres_sorted_index.at(2)),  EventWeight_);    
            FillTH1FAll(h1_mj3overmj2,  NFJbin, mj_thres_sorted.at(2)/mj_thres_sorted.at(1), EventWeight_);   
            FillTH1FAll(h1_FatjetPhi3,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(2)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta3,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(2)),      EventWeight_); 
//            FillTH1FAll(h1_N3,          NFJbin, FatjetN_->at(2),      EventWeight_); 
        }
        if(Nfatjet_thres>3) 
        {
            FillTH1FAll(h1_FatjetPt4,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(3)),       EventWeight_);   
            FillTH1FAll(h1_mj4,         NFJbin, mj_thres_sorted.at(3),  EventWeight_);    
            FillTH1FAll(h1_mj4OverMJ,   NFJbin, mj_thres_sorted.at(3)/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt4,   NFJbin, mj_thres_sorted.at(3)/FatjetPt_->at(mj_thres_sorted_index.at(3)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi4,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(3)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta4,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(3)),      EventWeight_); 
//            FillTH1FAll(h1_N4,          NFJbin, FatjetN_->at(3),      EventWeight_); 
        }
        if(Nfatjet_thres>4) 
        {
            // add 5th jet to h1_FatjetPt4
            FillTH1FAll(h1_FatjetPt4,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(4)),       EventWeight_);   
            FillTH1FAll(h1_mj4,         NFJbin, mj_thres_sorted.at(4),  EventWeight_);    
            FillTH1FAll(h1_mj4OverMJ,   NFJbin, mj_thres_sorted.at(4)/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt4,   NFJbin, mj_thres_sorted.at(4)/FatjetPt_->at(mj_thres_sorted_index.at(4)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi4,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(4)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta4,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(4)),      EventWeight_); 
//            FillTH1FAll(h1_N4,          NFJbin, FatjetN_->at(4),      EventWeight_); 
        }
        if(Nfatjet_thres>5) 
        {
            // add 6th jet to h1_FatjetPt4
            FillTH1FAll(h1_FatjetPt4,   NFJbin, FatjetPt_->at(mj_thres_sorted_index.at(5)),       EventWeight_);   
            FillTH1FAll(h1_mj4,         NFJbin, mj_thres_sorted.at(5),  EventWeight_);    
            FillTH1FAll(h1_mj4OverMJ,   NFJbin, mj_thres_sorted.at(5)/MJ_thres,  EventWeight_);    
            FillTH1FAll(h1_mjOverPt4,   NFJbin, mj_thres_sorted.at(5)/FatjetPt_->at(mj_thres_sorted_index.at(5)),  EventWeight_);    
            FillTH1FAll(h1_FatjetPhi4,  NFJbin, FatjetPhi_->at(mj_thres_sorted_index.at(5)),      EventWeight_); 
            FillTH1FAll(h1_FatjetEta4,  NFJbin, FatjetEta_->at(mj_thres_sorted_index.at(5)),      EventWeight_); 
//            FillTH1FAll(h1_N4,          NFJbin, FatjetN_->at(5),      EventWeight_); 
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
        FillTH1FAll(h1_Nskinnyjet,  NFJbin, Nskinnyjet_,    EventWeight_);      
        FillTH1FAll(h1_Ncsvm,       NFJbin, NBtagCSVM_,     EventWeight_);            
            
    } // for(int i = 0; i<nentries; i++)
   
    // 
    if(false && ChainName.Contains("TT_ll")) 
    { 
        cout << "---------------------------------" << endl;
        cout << "Nll  : " << Nll << endl;
        cout << "Nlt  : " << Nlt << endl;
        cout << "Ntt  : " << Ntt << endl;
        cout << "Nlth : " << Nlth << endl;
        cout << "Nltl : " << Nltl << endl;
        cout << "---------------------------------" << endl;
    }

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
        h1_MJ[i]->SetDirectory(0);                          h1_MJ[i]->Write();
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
    }
    HistFile->Close();
}
