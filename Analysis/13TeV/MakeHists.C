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
TH1F* InitTH1F(char* Name, char* Title, int Nbins, double XMin, double XMax){
    TH1F *h1 = new TH1F(Name, Title, Nbins, XMin, XMax);
    h1->Sumw2();
    return h1;
}
//
//TH2F initialization
//
TH2F* InitTH2F(char* Name, char* Title, int NXbins, double XMin, double XMax, int NYbins, double YMin, double YMax){
    TH2F *h2 = new TH2F(Name, Title, NXbins, XMin, XMax, NYbins, YMin, YMax);
    h2->Sumw2();
    return h2;
}
//
// Fill TH1F
//
TH1F* FillTH1F(TH1F* &h1, double var, double weight){
    if(var >= h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())) 
        var=h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())-0.00001;
    if(var < h1->GetXaxis()->GetBinLowEdge(1)) 
        var=h1->GetXaxis()->GetBinLowEdge(1)+0.00001;
    h1->Fill(var, weight);        
    return h1;
}
//
// Fill TH2F
//
TH2F* FillTH2F(TH2F* &h2, double varX, double varY, double weight){
    if(varX >= h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins())) 
        varX=h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins())-0.00001;
    if(varY >= h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins())) 
        varY=h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins())-0.00001;
    if(varX < h2->GetXaxis()->GetBinLowEdge(1)) 
        varX=h2->GetXaxis()->GetBinLowEdge(1)+0.00001;
    if(varY < h2->GetYaxis()->GetBinLowEdge(1)) 
        varY=h2->GetYaxis()->GetBinLowEdge(1)+0.00001;
    h2->Fill(varX, varY, weight);        
    return h2;
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
/* 
   // For SingleMu dataset
   if(pt<25) eff = 0;
   
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
void MakeHists(TChain *ch) 
{ 

    InitBaby(ch); 
    TString ChainName = ch->GetTitle();
    cout << "[MJ Analysis] " << ChainName << endl;  

    //
    //
    //
    TH1F *h1_yields[7];
    TH1F *h1_MJ[7], *h1_mT[7], *h1_Nskinnyjet[7], *h1_Ncsvm[7], 
         *h1_muspT[7], *h1_muspTminusMET[7], *h1_musEta[7], *h1_musPhi[7], 
         *h1_elspT[7], *h1_elspTminusMET[7], *h1_elsEta[7], *h1_elsPhi[7], 
         *h1_mj[7], *h1_FatjetPt[7], *h1_FatjetEta[7], 
         *h1_FatjetPt1[7],  *h1_FatjetPt2[7],  *h1_FatjetPt3[7],  *h1_FatjetPt4[7],
         *h1_mj1[7],  *h1_mj2[7],  *h1_mj3[7],  *h1_mj4[7], 
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
        // yields                    
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
                                    20, 0, 500);
        h1_mj3[i] = InitTH1F( Form("h1_%s_mj3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj3_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 500);
        h1_mj4[i] = InitTH1F( Form("h1_%s_mj4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_mj4_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 500);
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
        for(int i=0; i<GenId->size(); i++)
        { 
            if( TMath::Abs(GenMId->at(i))!=24 ) continue;
            if( TMath::Abs(GenId->at(i))!=11 && TMath::Abs(GenId->at(i))!=13 && TMath::Abs(GenId->at(i))!=15 ) continue;
            Ngenlep++;
        }
        if(ChainName.Contains("TT_sl") && Ngenlep!=1) continue;  
        if(ChainName.Contains("TT_ll") && Ngenlep!=2) continue;  

        // 
        // weights 
        // 
        // Temp fixes for wrong event weights 

        // TTJets cross section  
        if(ChainName.Contains("TT"))  
        { 
            EventWeight = EventWeight*806.1/832.; // I used 832 pb while M used 806.1 pb.
        }
        // Pileup 

        // 
        // single-lepton events selection
        // 
        if(!PassNLep(1))  continue;

        // Nfatjet counting with threshold 
        int Nfatjet_thres = 0;
        for(int ifatjet=0; ifatjet<(int)FatjetPt->size(); ifatjet++)
        {   
            float FatjetPt_cor = FatjetPt->at(ifatjet); 
            if(FatjetPt_cor>FatjetpTthres) Nfatjet_thres++;
        }

        //
        // Fill histograms 
        //
        float mT;
        float WpT;
        if(RA4MusPt->size()==1) {
            mT  = TMath::Sqrt( 2*MET*RA4MusPt->at(0)*(1-TMath::Cos(METPhi-RA4MusPhi->at(0))) ); 
            WpT = TMath::Sqrt(  
                        (RA4MusPt->at(0)*TMath::Cos(RA4MusPhi->at(0)) + MET*TMath::Cos(METPhi))
                       *(RA4MusPt->at(0)*TMath::Cos(RA4MusPhi->at(0)) + MET*TMath::Cos(METPhi))
                       +(RA4MusPt->at(0)*TMath::Sin(RA4MusPhi->at(0)) + MET*TMath::Sin(METPhi))
                       *(RA4MusPt->at(0)*TMath::Sin(RA4MusPhi->at(0)) + MET*TMath::Sin(METPhi))  ); 
        }
        if(RA4ElsPt->size()==1) {
            mT  = TMath::Sqrt( 2*MET*RA4ElsPt->at(0)*(1-TMath::Cos(METPhi-RA4ElsPhi->at(0))) ); 
            WpT = TMath::Sqrt(  
                        (RA4ElsPt->at(0)*TMath::Cos(RA4ElsPhi->at(0)) + MET*TMath::Cos(METPhi))
                       *(RA4ElsPt->at(0)*TMath::Cos(RA4ElsPhi->at(0)) + MET*TMath::Cos(METPhi))
                       +(RA4ElsPt->at(0)*TMath::Sin(RA4ElsPhi->at(0)) + MET*TMath::Sin(METPhi))
                       *(RA4ElsPt->at(0)*TMath::Sin(RA4ElsPhi->at(0)) + MET*TMath::Sin(METPhi))  ); 
        }

        // sort mj   
        vector<float> mj_sorted; 
        for(int i=0; i<(int)mj->size(); i++)  mj_sorted.push_back(mj->at(i));
        sort(mj_sorted.begin(), mj_sorted.end());
        reverse(mj_sorted.begin(), mj_sorted.end());

        ///// MJ 
        float MJ_tmp=0; 
        for(int i=0; i<(int)mj->size(); i++) 
        {
            if(FatjetPt->at(i)<FatjetpTthres) continue;
            MJ_tmp = MJ_tmp + mj->at(i);   
        }
        MJ = MJ_tmp;
        
        //
        // Selection   
        //
        
        // baseline selection
        if( !PassBaselineSelection() ) continue; 
        if( !PassSelection("SR0", HT, MET, NBtagCSVM, Nskinnyjet, mT, MJ)) continue;
 
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

        //if(NFJbin<2) continue; 

        // TT gen composition  
        int Ne=0;
        int Nm=0;
        int Nt=0;
        int Nth=0;
        int Ntl=0;
        int Ntn=0; int Nen=0; int Nmn=0;
        for(int i=0; i<GenId->size(); i++)
        { 
            if( TMath::Abs(GenMId->at(i))==24 && TMath::Abs(GenId->at(i))==11 ) Ne++;
            if( TMath::Abs(GenMId->at(i))==24 && TMath::Abs(GenId->at(i))==13 ) Nm++;
            if( TMath::Abs(GenMId->at(i))==24 && TMath::Abs(GenId->at(i))==15 ) Nt++;
            if( (TMath::Abs(GenGMId->at(i))==24||TMath::Abs(GenGMId->at(i))==15) && TMath::Abs(GenMId->at(i))==15 && TMath::Abs(GenId->at(i))==12 ) Nen++;
            if( (TMath::Abs(GenGMId->at(i))==24||TMath::Abs(GenGMId->at(i))==15) && TMath::Abs(GenMId->at(i))==15 && TMath::Abs(GenId->at(i))==14 ) Nmn++;
            if( (TMath::Abs(GenGMId->at(i))==24||TMath::Abs(GenGMId->at(i))==15) && TMath::Abs(GenMId->at(i))==15 && TMath::Abs(GenId->at(i))==16 ) Ntn++;
        }
        //cout << Ne << " " << Nm << " " << Nt << " " << Nen << " " << Nmn << " " << Ntn << endl;
        if( (Ne==2) || (Ne==1 && Nm==1) || (Nm==2) )                                    Nll  = Nll + EventWeight; 
        if( (Ne==1 || Nm==1) && (Nt==1) )                                               Nlt  = Nlt + EventWeight; 
        if( (Nt==2) )                                                                   Ntt  = Ntt + EventWeight; 
        if( (Ne==1 || Nm==1) && (Nt==1) && (Ntn==1 && (Nen==1 || Nmn==1)) )             Nltl = Nltl + EventWeight; 
        if( (Ne==1 || Nm==1) && (Nt==1) && (Ntn==1 && (Nen==0 && Nmn==0)) )             Nlth = Nlth + EventWeight; 

        //
        // Fill histogams 
        //
        
        // yields
        if(RA4ElsPt->size()==1) {FillTH1F(h1_yields[NFJbin], 0.5, EventWeight);   FillTH1F(h1_yields[6], 0.5, EventWeight);}
        if(RA4MusPt->size()==1) {FillTH1F(h1_yields[NFJbin], 1.5, EventWeight);   FillTH1F(h1_yields[6], 1.5, EventWeight);}
        
        // plots
        if(RA4MusPt->size()==1) {
            FillTH1F(h1_muspT[NFJbin],  RA4MusPt->at(0),  EventWeight); FillTH1F(h1_muspT[6],RA4MusPt->at(0), EventWeight);
            FillTH1F(h1_musEta[NFJbin], RA4MusEta->at(0), EventWeight); FillTH1F(h1_musEta[6], RA4MusEta->at(0), EventWeight);
            FillTH1F(h1_musPhi[NFJbin], RA4MusPhi->at(0), EventWeight); FillTH1F(h1_musPhi[6], RA4MusPhi->at(0), EventWeight);
            if(getDPhi(RA4MusPhi->at(0),METPhi)<0.4) 
            { 
                FillTH1F(h1_muspTminusMET[NFJbin], (MET-RA4MusPt->at(0))/RA4MusPt->at(0), EventWeight);   
                FillTH1F(h1_muspTminusMET[6],(MET-RA4MusPt->at(0))/RA4MusPt->at(0), EventWeight);
            }
        }
        if(RA4ElsPt->size()==1) {
            FillTH1F(h1_elspT[NFJbin], RA4ElsPt->at(0), EventWeight);   FillTH1F(h1_elspT[6],RA4ElsPt->at(0), EventWeight);
            FillTH1F(h1_elsEta[NFJbin], RA4ElsEta->at(0), EventWeight); FillTH1F(h1_elsEta[6], RA4ElsEta->at(0), EventWeight);
            FillTH1F(h1_elsPhi[NFJbin], RA4ElsPhi->at(0), EventWeight); FillTH1F(h1_elsPhi[6], RA4ElsPhi->at(0), EventWeight);
            if(getDPhi(RA4ElsPhi->at(0),METPhi)<0.4) 
            { 
                FillTH1F(h1_elspTminusMET[NFJbin], (MET-RA4ElsPt->at(0))/RA4ElsPt->at(0), EventWeight);   
                FillTH1F(h1_elspTminusMET[6],(MET-RA4ElsPt->at(0))/RA4ElsPt->at(0), EventWeight);
            }
        }

        if(NFJbin>1){ FillTH2F(h2_mj1vsmj2[NFJbin], mj_sorted.at(0), mj_sorted.at(1), EventWeight);           
                      FillTH2F(h2_mj1vsmj2[6], mj_sorted.at(0), mj_sorted.at(1), EventWeight); }           
        if(NFJbin>2){ FillTH2F(h2_mj2vsmj3[NFJbin], mj_sorted.at(1), mj_sorted.at(2), EventWeight);           
                      FillTH2F(h2_mj2vsmj3[6], mj_sorted.at(1), mj_sorted.at(2), EventWeight); }           
        if(NFJbin>3){ FillTH2F(h2_mj3vsmj4[NFJbin], mj_sorted.at(2), mj_sorted.at(3), EventWeight);           
                      FillTH2F(h2_mj3vsmj4[6], mj_sorted.at(2), mj_sorted.at(3), EventWeight); }           
        FillTH2F(h2_HTMET[NFJbin], HT, MET, EventWeight);           FillTH2F(h2_HTMET[6], HT, MET, EventWeight);
        FillTH2F(h2_MJmT[NFJbin], MJ, mT, EventWeight);             FillTH2F(h2_MJmT[6], MJ, mT, EventWeight);
        FillTH2F(h2_HTmT[NFJbin], HT, mT, EventWeight);             FillTH2F(h2_HTmT[6], HT, mT, EventWeight);
        FillTH2F(h2_MJMET[NFJbin], MJ, MET, EventWeight);           FillTH2F(h2_MJMET[6], MJ, MET, EventWeight);
        FillTH2F(h2_HTMJ[NFJbin], HT, MJ,EventWeight);              FillTH2F(h2_HTMJ[6], HT, MJ, EventWeight);
        FillTH2F(h2_METmT[NFJbin], MET, mT, EventWeight);           FillTH2F(h2_METmT[6], MET, mT, EventWeight);
        FillTH1F(h1_mT[NFJbin], mT, EventWeight);                   FillTH1F(h1_mT[6], mT, EventWeight);
        FillTH1F(h1_WpT[NFJbin], WpT, EventWeight);                 FillTH1F(h1_WpT[6], WpT, EventWeight);
        FillTH1F(h1_MJ[NFJbin], MJ, EventWeight);                   FillTH1F(h1_MJ[6], MJ, EventWeight);
        FillTH1F(h1_HT[NFJbin], HT, EventWeight);                   FillTH1F(h1_HT[6], HT, EventWeight);
        FillTH1F(h1_MET[NFJbin], MET, EventWeight);                 FillTH1F(h1_MET[6], MET, EventWeight);
        FillTH1F(h1_METPhi[NFJbin], METPhi, EventWeight);           FillTH1F(h1_METPhi[6], METPhi, EventWeight);
        FillTH1F(h1_METx[NFJbin], MET*TMath::Cos(METPhi), EventWeight);           FillTH1F(h1_METx[6], MET*TMath::Cos(METPhi), EventWeight);
        FillTH1F(h1_METy[NFJbin], MET*TMath::Sin(METPhi), EventWeight);           FillTH1F(h1_METy[6], MET*TMath::Sin(METPhi), EventWeight);
        if(RA4MusPt->size()==1) {FillTH1F(h1_DPhi[NFJbin], getDPhi(RA4MusPhi->at(0),METPhi), EventWeight);   FillTH1F(h1_DPhi[6], getDPhi(RA4MusPhi->at(0),METPhi), EventWeight);}
        if(RA4ElsPt->size()==1) {FillTH1F(h1_DPhi[NFJbin], getDPhi(RA4ElsPhi->at(0),METPhi), EventWeight);   FillTH1F(h1_DPhi[6], getDPhi(RA4ElsPhi->at(0),METPhi), EventWeight);}
        
        if(Nfatjet_thres>0) 
        {
            FillTH1F(h1_FatjetPt1[NFJbin], FatjetPt->at(0), EventWeight);   FillTH1F(h1_FatjetPt1[6], FatjetPt->at(0), EventWeight); 
            FillTH1F(h1_mj1[NFJbin], mj_sorted.at(0), EventWeight);   FillTH1F(h1_mj1[6], mj_sorted.at(0), EventWeight); 
            FillTH1F(h1_FatjetPhi1[NFJbin], FatjetPhi->at(0), EventWeight); FillTH1F(h1_FatjetPhi1[6], FatjetPhi->at(0), EventWeight); 
            FillTH1F(h1_FatjetEta1[NFJbin], FatjetEta->at(0), EventWeight); FillTH1F(h1_FatjetEta1[6], FatjetEta->at(0), EventWeight); 
        }
        if(Nfatjet_thres>1) 
        {
            FillTH1F(h1_FatjetPt2[NFJbin], FatjetPt->at(1), EventWeight);   FillTH1F(h1_FatjetPt2[6], FatjetPt->at(1), EventWeight); 
            FillTH1F(h1_mj2[NFJbin], mj_sorted.at(1), EventWeight);   FillTH1F(h1_mj2[6], mj_sorted.at(1), EventWeight); 
            FillTH1F(h1_mj2overmj1[NFJbin], mj_sorted.at(1)/mj_sorted.at(0), EventWeight);   FillTH1F(h1_mj2overmj1[6], mj_sorted.at(1)/mj_sorted.at(0), EventWeight); 
            FillTH1F(h1_FatjetPhi2[NFJbin], FatjetPhi->at(1), EventWeight); FillTH1F(h1_FatjetPhi2[6], FatjetPhi->at(1), EventWeight); 
            FillTH1F(h1_FatjetEta2[NFJbin], FatjetEta->at(1), EventWeight); FillTH1F(h1_FatjetEta2[6], FatjetEta->at(1), EventWeight); 
        }
        if(Nfatjet_thres>2) 
        {
            FillTH1F(h1_FatjetPt3[NFJbin], FatjetPt->at(2), EventWeight);   FillTH1F(h1_FatjetPt3[6], FatjetPt->at(2), EventWeight); 
            FillTH1F(h1_mj3[NFJbin], mj_sorted.at(2), EventWeight);   FillTH1F(h1_mj3[6], mj_sorted.at(2), EventWeight); 
            FillTH1F(h1_mj3overmj2[NFJbin], mj_sorted.at(2)/mj_sorted.at(1), EventWeight);   FillTH1F(h1_mj3overmj2[6], mj_sorted.at(2)/mj_sorted.at(1), EventWeight); 
            FillTH1F(h1_FatjetPhi3[NFJbin], FatjetPhi->at(2), EventWeight); FillTH1F(h1_FatjetPhi3[6], FatjetPhi->at(2), EventWeight); 
            FillTH1F(h1_FatjetEta3[NFJbin], FatjetEta->at(2), EventWeight); FillTH1F(h1_FatjetEta3[6], FatjetEta->at(2), EventWeight); 
        }
        if(Nfatjet_thres>3) 
        {
            FillTH1F(h1_FatjetPt4[NFJbin], FatjetPt->at(3), EventWeight);   FillTH1F(h1_FatjetPt4[6], FatjetPt->at(3), EventWeight); 
            FillTH1F(h1_mj4[NFJbin], mj_sorted.at(3), EventWeight);   FillTH1F(h1_mj4[6], mj_sorted.at(3), EventWeight); 
            FillTH1F(h1_FatjetPhi4[NFJbin], FatjetPhi->at(3), EventWeight); FillTH1F(h1_FatjetPhi4[6], FatjetPhi->at(3), EventWeight); 
            FillTH1F(h1_FatjetEta4[NFJbin], FatjetEta->at(3), EventWeight); FillTH1F(h1_FatjetEta4[6], FatjetEta->at(3), EventWeight); 
        }
        if(Nfatjet_thres>4) 
        {
            // add 5th jet to h1_FatjetPt4
            FillTH1F(h1_FatjetPt4[NFJbin], FatjetPt->at(4), EventWeight);   FillTH1F(h1_FatjetPt4[6], FatjetPt->at(4), EventWeight); 
            FillTH1F(h1_mj4[NFJbin], mj_sorted.at(4), EventWeight);   FillTH1F(h1_mj4[6], mj_sorted.at(4), EventWeight); 
            FillTH1F(h1_FatjetPhi4[NFJbin], FatjetPhi->at(4), EventWeight); FillTH1F(h1_FatjetPhi4[6], FatjetPhi->at(4), EventWeight); 
            FillTH1F(h1_FatjetEta4[NFJbin], FatjetEta->at(4), EventWeight); FillTH1F(h1_FatjetEta4[6], FatjetEta->at(4), EventWeight); 
        }
        if(Nfatjet_thres>5) 
        {
            // add 6th jet to h1_FatjetPt4
            FillTH1F(h1_FatjetPt4[NFJbin], FatjetPt->at(5), EventWeight);   FillTH1F(h1_FatjetPt4[6], FatjetPt->at(5), EventWeight); 
            FillTH1F(h1_mj4[NFJbin], mj_sorted.at(5), EventWeight);   FillTH1F(h1_mj4[6], mj_sorted.at(5), EventWeight); 
            FillTH1F(h1_FatjetPhi4[NFJbin], FatjetPhi->at(5), EventWeight); FillTH1F(h1_FatjetPhi4[6], FatjetPhi->at(5), EventWeight); 
            FillTH1F(h1_FatjetEta4[NFJbin], FatjetEta->at(5), EventWeight); FillTH1F(h1_FatjetEta4[6], FatjetEta->at(5), EventWeight); 
        }

        for(int i=0; i<(int)mj->size(); i++) 
        { 
            if(FatjetPt->at(i)<FatjetpTthres) continue;

            FillTH1F(h1_mj[NFJbin], mj->at(i), EventWeight);                            FillTH1F(h1_mj[6], mj->at(i), EventWeight);
            FillTH1F(h1_FatjetPt[NFJbin], FatjetPt->at(i), EventWeight);                FillTH1F(h1_FatjetPt[6], FatjetPt->at(i), EventWeight);
            FillTH1F(h1_FatjetEta[NFJbin], FatjetEta->at(i), EventWeight);              FillTH1F(h1_FatjetEta[6], FatjetEta->at(i), EventWeight);
        }

        FillTH1F(h1_Nfatjet[NFJbin], Nfatjet_thres, EventWeight);   FillTH1F(h1_Nfatjet[6], Nfatjet_thres, EventWeight);
        FillTH1F(h1_Nskinnyjet[NFJbin], Nskinnyjet, EventWeight);   FillTH1F(h1_Nskinnyjet[6], Nskinnyjet, EventWeight);
        FillTH1F(h1_Ncsvm[NFJbin], NBtagCSVM, EventWeight);   FillTH1F(h1_Ncsvm[6], NBtagCSVM, EventWeight);
            
    } // for(int i = 0; i<nentries; i++)
    if(ChainName.Contains("TT_ll")) 
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
    HistFileName = Form("HistFiles/%s.root", HistFileName.Data());
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
        h1_mj3overmj2[i]->SetDirectory(0);                   h1_mj3overmj2[i]->Write();
        h1_mj2overmj1[i]->SetDirectory(0);                   h1_mj2overmj1[i]->Write();
        h1_FatjetPhi1[i]->SetDirectory(0);                   h1_FatjetPhi1[i]->Write();
        h1_FatjetPhi2[i]->SetDirectory(0);                   h1_FatjetPhi2[i]->Write();
        h1_FatjetPhi3[i]->SetDirectory(0);                   h1_FatjetPhi3[i]->Write();
        h1_FatjetPhi4[i]->SetDirectory(0);                   h1_FatjetPhi4[i]->Write();
        h1_FatjetEta1[i]->SetDirectory(0);                   h1_FatjetEta1[i]->Write();
        h1_FatjetEta2[i]->SetDirectory(0);                   h1_FatjetEta2[i]->Write();
        h1_FatjetEta3[i]->SetDirectory(0);                   h1_FatjetEta3[i]->Write();
        h1_FatjetEta4[i]->SetDirectory(0);                   h1_FatjetEta4[i]->Write();
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
