#include <iostream>

#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLorentzVector.h"

//#include "babytree.h"

using namespace std;
bool DoLog          = 0;
bool OnlyDraw       = 1;
bool RemoveMuon     = 0;
bool Do3FatjetScale = 0;
int pTR0p5thres     = 30;
int FatjetpTthres   = 50;

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
// h1 cosmetics
//
void h1cosmetic(TH1F* &h1, char* title, int linecolor=kBlack, int linewidth=1, int fillcolor=0, TString var="")
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.1);
}

//
// h2 cosmetics
//
void h2cosmetic(TH2F* &h2, char* title, TString Xvar="", TString Yvar="", TString Zvar="")
{
    h2->SetTitle(title);
    h2->SetXTitle(Xvar);
    h2->SetYTitle(Yvar);
    h2->SetZTitle(Yvar);
    h2->SetStats(0);
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
    else 2*TMath::Pi() - absdphi;
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
void DoOneProcess(TChain *ch, int pTR0p5) 
{ 

    float           EventWeight;
    float           Npu;
    int             Npv;
    int             Nfatjet;
    int             Nskinnyjet;
    int             NBtagCSVM;
    float           MJ;
    float           MET;
    float           HT;
    float           METPhi;
    vector<float>   *mj;
    vector<float>   *FatjetPt;
    vector<float>   *FatjetEta;
    vector<float>   *FatjetPhi;
    vector<float>   *RA4MusPt;
    vector<float>   *RA4MusPhi;
    vector<float>   *RA4MusEta;
    vector<float>   *JetPt;
    vector<float>   *JetEta;
    vector<float>   *JetPhi;
    vector<float>   *JetCSV;
    vector<float>   *RA4MusVetoPt;
    vector<float>   *RA4ElsVetoPt;

    TBranch         *b_EventWeight;   //!
    TBranch         *b_Npu;   //!
    TBranch         *b_Npv;   //!
    TBranch         *b_Nfatjet;   //!
    TBranch         *b_Nskinnyjet;   //!
    TBranch         *b_NBtagCSVM;   //!
    TBranch         *b_MJ;   //!
    TBranch         *b_MET;   //!
    TBranch         *b_HT;   //!
    TBranch         *b_METPhi;   //!
    TBranch         *b_mj;   //!
    TBranch         *b_FatjetPt;   //!
    TBranch         *b_FatjetEta;   //!
    TBranch         *b_FatjetPhi;   //!
    TBranch         *b_RA4MusPt;   //!
    TBranch         *b_RA4MusPhi;   //!
    TBranch         *b_RA4MusEta;   //!
    TBranch         *b_JetPt;   //!
    TBranch         *b_JetEta;   //!
    TBranch         *b_JetPhi;   //!
    TBranch         *b_JetCSV;   //!
    TBranch         *b_RA4MusVetoPt;   //!
    TBranch         *b_RA4ElsVetoPt;   //!
    mj        = 0;
    FatjetPt  = 0;
    FatjetEta = 0;
    FatjetPhi = 0;
    RA4MusPt  = 0;
    RA4MusPhi = 0;
    RA4MusEta = 0;
    JetPt  = 0;
    JetEta = 0;
    JetPhi = 0;
    JetCSV = 0;
    RA4MusVetoPt  = 0;
    RA4ElsVetoPt  = 0;

    ch->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
    ch->SetBranchAddress("Npu", &Npu, &b_Npu);
    ch->SetBranchAddress("Npv", &Npv, &b_Npv);
    ch->SetBranchAddress(Form("Nfatjet_pT%i", pTR0p5), &Nfatjet, &b_Nfatjet);
    ch->SetBranchAddress("Nskinnyjet", &Nskinnyjet, &b_Nskinnyjet);
    ch->SetBranchAddress("NBtagCSVM", &NBtagCSVM, &b_NBtagCSVM);
    ch->SetBranchAddress(Form("MJ_pT%i", pTR0p5), &MJ, &b_MJ);
    ch->SetBranchAddress("MET", &MET, &b_MET);
    ch->SetBranchAddress("HT", &HT, &b_HT);
    ch->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
    ch->SetBranchAddress(Form("mj_pT%i", pTR0p5), &mj, &b_mj);
    ch->SetBranchAddress(Form("FatjetPt_pT%i", pTR0p5), &FatjetPt, &b_FatjetPt);
    ch->SetBranchAddress(Form("FatjetEta_pT%i", pTR0p5), &FatjetEta, &b_FatjetEta);
    ch->SetBranchAddress(Form("FatjetPhi_pT%i", pTR0p5), &FatjetPhi, &b_FatjetPhi);
    ch->SetBranchAddress("RA4MusPt", &RA4MusPt, &b_RA4MusPt);
    ch->SetBranchAddress("RA4MusPhi", &RA4MusPhi, &b_RA4MusPhi);
    ch->SetBranchAddress("RA4MusEta", &RA4MusEta, &b_RA4MusEta);
    ch->SetBranchAddress("JetPt",  &JetPt,  &b_JetPt);
    ch->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
    ch->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
    ch->SetBranchAddress("JetCSV", &JetCSV, &b_JetCSV);
    ch->SetBranchAddress("RA4MusVetoPt", &RA4MusVetoPt, &b_RA4MusVetoPt);
    ch->SetBranchAddress("RA4ElsVetoPt", &RA4ElsVetoPt, &b_RA4ElsVetoPt);
    
    TString ChainName = ch->GetTitle();
    cout << "[MJ Analysis] " << ChainName << endl;  

    bool        TrigMuon;
    bool        TrigSingleMuon;
    TBranch     *b_TrigMuon; 
    TBranch     *b_TrigSingleMuon; 
    if(ChainName.Contains("DATA"))  
    { 
        ch->SetBranchAddress("TrigMuon", &TrigMuon, &b_TrigMuon);
        ch->SetBranchAddress("TrigSingleMuon", &TrigSingleMuon, &b_TrigSingleMuon);
    }
    
    float        top1pT;
    float        top1Phi;
    float        top2pT;
    float        top2Phi;
    TBranch     *b_top1pT; 
    TBranch     *b_top1Phi; 
    TBranch     *b_top2pT; 
    TBranch     *b_top2Phi; 
    if(ChainName.Contains("TT") )  
    { 
        ch->SetBranchAddress("top1pT", &top1pT, &b_top1pT);
        ch->SetBranchAddress("top1Phi", &top1Phi, &b_top1Phi);
        ch->SetBranchAddress("top2pT", &top2pT, &b_top2pT);
        ch->SetBranchAddress("top2Phi", &top2Phi, &b_top2Phi);
    }

    //
    //
    //
    TH1F *h1_MJ[6], *h1_mT[6], *h1_Nskinnyjet[6], *h1_muspT[6], *h1_musEta[6], *h1_musPhi[6], 
         *h1_mj[6], *h1_FatjetPt[6], 
         *h1_mj_mumatch[6], *h1_FatjetPt_mumatch[6], *h1_mj_notmumatch[6], *h1_FatjetPt_notmumatch[6], 
         *h1_mj_bmatch[6], *h1_FatjetPt_bmatch[6], *h1_mj_notbmatch[6], *h1_FatjetPt_notbmatch[6], 
         *h1_mj_metmatch[6], *h1_FatjetPt_metmatch[6], *h1_mj_notmetmatch[6], *h1_FatjetPt_notmetmatch[6], 
         *h1_HT[6], *h1_MET[6], *h1_METPhi[6], *h1_DPhi[6], *h1_Nfatjet[6], *h1_WpT[6], *h1_dRmin[6],
         *h1_genpTttbar[6];
    TH2F *h2_mjvsFatjetPt[6];
    TH2F *h2_mjvsmj_2fatjet;

    for(int i=0; i<6; i++) 
    {
        h1_MJ[i] = InitTH1F( Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 1000);
        h1_mj[i] = InitTH1F( Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             40, 0, 500);
        h1_mj_mumatch[i] = InitTH1F( Form("h1_%s_mj_mumatch_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_mj_mumatch_%ifatjet", ch->GetTitle(), i), 
                                   40, 0, 500);
        h1_mj_notmumatch[i] = InitTH1F( Form("h1_%s_mj_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                      Form("h1_%s_mj_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                      40, 0, 500);
        h1_mj_bmatch[i] = InitTH1F( Form("h1_%s_mj_bmatch_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_mj_bmatch_%ifatjet", ch->GetTitle(), i), 
                                   40, 0, 500);
        h1_mj_notbmatch[i] = InitTH1F( Form("h1_%s_mj_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                      Form("h1_%s_mj_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                      40, 0, 500);
        h1_mj_metmatch[i] = InitTH1F( Form("h1_%s_mj_metmatch_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_mj_metmatch_%ifatjet", ch->GetTitle(), i), 
                                   40, 0, 500);
        h1_mj_notmetmatch[i] = InitTH1F( Form("h1_%s_mj_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                      Form("h1_%s_mj_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                      40, 0, 500);
        h1_Nfatjet[i] = InitTH1F( Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  11, -0.5, 10.5);
        h1_Nskinnyjet[i] = InitTH1F( Form("h1_%s_Nskinnyjet_%ifatjet", ch->GetTitle(), i), 
                                     Form("h1_%s_Nskinnyjet_%ifatjet", ch->GetTitle(), i), 
                                     20, -0.5, 19.5);
        h1_mT[i] = InitTH1F( Form("h1_%s_mT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mT_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 200);
                             20, 0, 600);
        h1_genpTttbar[i] = InitTH1F( Form("h1_%s_genpTttbar_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_genpTttbar_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 200);
                             20, 0, 600);
        h1_muspT[i] = InitTH1F( Form("h1_%s_muspT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_muspT_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 200);
        h1_musEta[i] = InitTH1F( Form("h1_%s_musEta_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_musEta_%ifatjet", ch->GetTitle(), i), 
                                20, -2.5, 2.5);
                                //100, -2.1, 2.1);
        h1_musPhi[i] = InitTH1F( Form("h1_%s_musPhi_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_musPhi_%ifatjet", ch->GetTitle(), i), 
                                20, -3.2, 3.2);
        h1_HT[i] = InitTH1F( Form("h1_%s_HT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_HT_%ifatjet", ch->GetTitle(), i), 
                             //20, 350, 1350);
                             20, 500, 2000);
        h1_MET[i] = InitTH1F( Form("h1_%s_MET_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MET_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             20, 0, 1000);
        h1_METPhi[i] = InitTH1F( Form("h1_%s_METPhi_%ifatjet", ch->GetTitle(), i), 
                                 Form("h1_%s_METPhi_%ifatjet", ch->GetTitle(), i), 
                                 20, -3.2, 3.2);
        h1_DPhi[i] = InitTH1F( Form("h1_%s_DPhi_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_DPhi_%ifatjet", ch->GetTitle(), i), 
                               20, -2, 2);
        h1_dRmin[i] = InitTH1F( Form("h1_%s_dRmin_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dRmin_%ifatjet", ch->GetTitle(), i), 
                               20, 0, 2);
        h1_FatjetPt[i] = InitTH1F( Form("h1_%s_FatjetPt_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_FatjetPt_%ifatjet", ch->GetTitle(), i), 
                                   20, 0, 800);
        h1_FatjetPt_mumatch[i] = InitTH1F( Form("h1_%s_FatjetPt_mumatch_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_FatjetPt_mumatch_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 800);
        h1_FatjetPt_notmumatch[i] = InitTH1F( Form("h1_%s_FatjetPt_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_FatjetPt_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                            20, 0, 800);
        h1_FatjetPt_bmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_bmatch_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_FatjetPt_bmatch_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 800);
        h1_FatjetPt_notbmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_FatjetPt_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                            20, 0, 800);
        h1_FatjetPt_metmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_metmatch_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_FatjetPt_metmatch_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 800);
        h1_FatjetPt_notmetmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_FatjetPt_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                            20, 0, 800);
        h1_WpT[i] = InitTH1F( Form("h1_%s_WpT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_WpT_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             20, 0, 1000);
        h2_mjvsFatjetPt[i]    =   InitTH2F(Form("h2_%s_mjvsFatjetPt_%ifatjet", ch->GetTitle(), i),
                                           Form("h2_%s_mjvsFatjetPt_%ifatjet", ch->GetTitle(), i), 
                                           20, 0, 500, 20, 0, 800);
    }
    h2_mjvsmj_2fatjet = InitTH2F(Form("h2_%s_mjvsmj_%ifatjet", ch->GetTitle(), 2),
                                 Form("h2_%s_mjvsmj_%ifatjet", ch->GetTitle(), 2),
                                 20, 0, 500, 20, 0, 500);;
    // Pile up reweighting hist
    //TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_1p317ipb.root");
    TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_19p5ifb.root");
    TH1F *h1PU = (TH1F*)(fPUFile->Get("puWeights"));
    
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
      
        bool cut_others = false;

        // 
        // weights 
        // 
        // Temp fixes for wrong event weights 

        // Data : trigger
        if(ChainName.Contains("DATA"))  
        { 
            EventWeight = EventWeight * (TrigMuon); 
        }
        // Pileup 
        if(!ChainName.Contains("DATA")) 
        {
            EventWeight = EventWeight/19500*19600;                      // luminostiry correction   
            EventWeight = EventWeight*nPUScaleFactor2012(h1PU, Npu);    // PU   
            EventWeight = EventWeight*0.94;                             // Btagging SF (https://twiki.cern.ch/twiki/pub/CMSPublic/TWikiBTV_Moriond2013/Fig17b.pdf) 
        }
        // TT reweighting
        if(ChainName.Contains("TT")) 
        { 
            float weight_top1pT = TMath::Exp(0.156-0.00137*top1pT);
            float weight_top2pT = TMath::Exp(0.156-0.00137*top2pT);
            EventWeight = EventWeight * TMath::Sqrt(weight_top1pT*weight_top2pT);
        }
        // Signal cross section
        // 0.000871201 in SUSY xsec twiki : https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVgluglu 
        // 0.000256511 in the cfA model param
        // f1400_25   : 137498 events 
        // f1400_1000 : 126811 events 
        if(ChainName.Contains("f1400_25")) 
        {
            EventWeight = EventWeight*19600*0.000871201/137498*1000;                      
        }
        if(ChainName.Contains("f1400_1000")) 
        {
            EventWeight = EventWeight*19600*0.000871201/126811*1000;                      
        }

        //
        // cuts applied  
        //

        // select only single-muon events
        if( RA4MusPt->size()!=1 ) continue;
        if( RA4MusVetoPt->size()>0 ) continue;
        if( RA4ElsVetoPt->size()>0 ) continue;

        // single lepton trigger efficienc from HWW 
        if(!ChainName.Contains("DATA")) EventWeight = EventWeight * GetMuonEff(RA4MusPt->at(0), RA4MusEta->at(0));

        // selection
        if( HT>500                              && 
            MET>100                             && 
            RA4MusPt->at(0)>20                  &&  
            NBtagCSVM>1        
           ) cut_others = true;

        // Nfatjet counting with threshold 
        int Nfatjet_thres = 0;
        for(int ifatjet=0; ifatjet<(int)FatjetPt->size(); ifatjet++)
        {   
            float FatjetPt_cor = FatjetPt->at(ifatjet); 
//            if(RemoveMuon) 
//            {
//                for(int ilep=0; ilep<(int)RA4MusPt->size(); ilep++)  
//                {
//                    float dEta = FatjetEta->at(ifatjet) - RA4MusEta->at(ilep);
//                    float dPhi = FatjetPhi->at(ifatjet) - RA4MusPhi->at(ilep);
//                    if(TMath::Sqrt(dEta*dEta+dPhi*dPhi)<1.0) FatjetPt_cor = FatjetPt_cor - RA4MusPt->at(ilep);
//                }
//            } 
            if(FatjetPt_cor>FatjetpTthres) Nfatjet_thres++;
        }

        //
        // Fill histograms 
        //
        ///// mT 
        float mT  = TMath::Sqrt( 2*MET*RA4MusPt->at(0)*(1-TMath::Cos(METPhi-RA4MusPhi->at(0))) ); 
        ///// WpT 
        float WpT =  TMath::Sqrt(  
                        (RA4MusPt->at(0)*TMath::Cos(RA4MusPhi->at(0)) + MET*TMath::Cos(METPhi))
                       *(RA4MusPt->at(0)*TMath::Cos(RA4MusPhi->at(0)) + MET*TMath::Cos(METPhi))
                       +(RA4MusPt->at(0)*TMath::Sin(RA4MusPhi->at(0)) + MET*TMath::Sin(METPhi))
                       *(RA4MusPt->at(0)*TMath::Sin(RA4MusPhi->at(0)) + MET*TMath::Sin(METPhi))  ); 
        ///// MJ 
        float MJ_tmp=0; 
        for(int i=0; i<(int)mj->size(); i++) 
        {
            if(FatjetPt->at(i)<FatjetpTthres) continue;
            MJ_tmp = MJ_tmp + mj->at(i);   
        }
        MJ = MJ_tmp;
        ///// gen ttbar pT
        //float genpTttbar = TMath::Sqrt( 
        //                       (top1pT*TMath::Sin(top1Phi)+top2pT*TMath::Sin(top2Phi))*(top1pT*TMath::Sin(top1Phi)+top2pT*TMath::Sin(top2Phi))
        //                      +(top1pT*TMath::Cos(top1Phi)+top2pT*TMath::Cos(top2Phi))*(top1pT*TMath::Cos(top1Phi)+top2pT*TMath::Cos(top2Phi)) ); 

        if(cut_others) {

            if(Nfatjet_thres>3) 
            {
                h1_muspT[4]->Fill(RA4MusPt->at(0), EventWeight);
                h1_musEta[4]->Fill(RA4MusEta->at(0), EventWeight);
                h1_musPhi[4]->Fill(RA4MusPhi->at(0), EventWeight);
                h1_mT[4]->Fill(mT, EventWeight);
                //h1_genpTttbar[4]->Fill(genpTttbar, EventWeight);
                h1_WpT[4]->Fill(WpT, EventWeight);
                h1_MJ[4]->Fill(MJ, EventWeight);
                h1_HT[4]->Fill(HT, EventWeight);
                h1_MET[4]->Fill(MET, EventWeight);
                h1_METPhi[4]->Fill(METPhi, EventWeight);
                h1_DPhi[4]->Fill(RA4MusPhi->at(0)-METPhi, EventWeight);
                for(int i=0; i<(int)mj->size(); i++) 
                { 
                    if(FatjetPt->at(i)<FatjetpTthres) continue;

                    h1_mj[4]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt[4]->Fill(FatjetPt->at(i), EventWeight);
                    h2_mjvsFatjetPt[4]->Fill(mj->at(i), FatjetPt->at(i), EventWeight);
                    // met matching 
                    if(getDPhi(FatjetPhi->at(i), METPhi)<1)
                    { 
                        h1_mj_metmatch[4]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_metmatch[4]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    else 
                    { 
                        h1_mj_notmetmatch[4]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmetmatch[4]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // mu matching 
                    if(getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0))<1) 
                    { 
                        h1_mj_mumatch[4]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_mumatch[4]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notmumatch[4]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmumatch[4]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // b matching
                    float dRmin=999.;
                    for(int j=0; j<(int)JetPt->size(); j++) 
                    { 
                        if(JetCSV->at(j) < 0.679) continue;
                        float dRtmp = getDR(JetEta->at(j), RA4MusEta->at(0), JetPhi->at(j), RA4MusPhi->at(0)); 
                        if(dRtmp < dRmin) dRmin = dRtmp; 
                    }
                    if(dRmin<1.)
                    { 
                        h1_mj_bmatch[4]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_bmatch[4]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notbmatch[4]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notbmatch[4]->Fill(FatjetPt->at(i), EventWeight);
                    }
                }
                h1_Nfatjet[4]->Fill(Nfatjet_thres, EventWeight);
                h1_Nskinnyjet[4]->Fill(Nskinnyjet, EventWeight);
            }
            else if(Nfatjet_thres>2)
            {
                if(!ChainName.Contains("DATA") && Do3FatjetScale) EventWeight = EventWeight * 1.15; // Scale factor for 3 fatjet  

                h1_muspT[3]->Fill(RA4MusPt->at(0), EventWeight);
                h1_musEta[3]->Fill(RA4MusEta->at(0), EventWeight);
                h1_musPhi[3]->Fill(RA4MusPhi->at(0), EventWeight);
                h1_mT[3]->Fill(mT, EventWeight);
                //h1_genpTttbar[3]->Fill(genpTttbar, EventWeight);
                h1_WpT[3]->Fill(WpT, EventWeight);
                h1_MJ[3]->Fill(MJ, EventWeight);
                h1_HT[3]->Fill(HT, EventWeight);
                h1_MET[3]->Fill(MET, EventWeight);
                h1_METPhi[3]->Fill(METPhi, EventWeight);
                h1_DPhi[3]->Fill(RA4MusPhi->at(0)-METPhi, EventWeight);
                for(int i=0; i<(int)mj->size(); i++) 
                { 
                    if(FatjetPt->at(i)<FatjetpTthres) continue;

                    h1_mj[3]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt[3]->Fill(FatjetPt->at(i), EventWeight);
                    h2_mjvsFatjetPt[3]->Fill(mj->at(i), FatjetPt->at(i), EventWeight);
                    // met matching 
                    if(getDPhi(FatjetPhi->at(i), METPhi)<1)
                    { 
                        h1_mj_metmatch[3]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_metmatch[3]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    else 
                    { 
                        h1_mj_notmetmatch[3]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmetmatch[3]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // b matching
                    float dRmin=999.;
                    for(int j=0; j<(int)JetPt->size(); j++) 
                    { 
                        if(JetCSV->at(j) < 0.679) continue;
                        float dRtmp = getDR(JetEta->at(j), RA4MusEta->at(0), JetPhi->at(j), RA4MusPhi->at(0)); 
                        if(dRtmp < dRmin) dRmin = dRtmp; 
                    }
                    if(dRmin<1.)
                    { 
                        h1_mj_bmatch[3]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_bmatch[3]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notbmatch[3]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notbmatch[3]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // mu matching
                    if(getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0))<1) 
                    { 
                        h1_mj_mumatch[3]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_mumatch[3]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notmumatch[3]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmumatch[3]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    
                }
                h1_Nfatjet[3]->Fill(Nfatjet_thres, EventWeight);
                h1_Nskinnyjet[3]->Fill(Nskinnyjet, EventWeight);
            }
            else if(Nfatjet_thres>1)
            {
                h1_muspT[2]->Fill(RA4MusPt->at(0), EventWeight);
                h1_musEta[2]->Fill(RA4MusEta->at(0), EventWeight);
                h1_musPhi[2]->Fill(RA4MusPhi->at(0), EventWeight);
                h1_mT[2]->Fill(mT, EventWeight);
                //h1_genpTttbar[2]->Fill(genpTttbar, EventWeight);
                h1_WpT[2]->Fill(WpT, EventWeight);
                h1_MJ[2]->Fill(MJ, EventWeight);
                h1_HT[2]->Fill(HT, EventWeight);
                h1_MET[2]->Fill(MET, EventWeight);
                h1_METPhi[2]->Fill(METPhi, EventWeight);
                h1_DPhi[2]->Fill(RA4MusPhi->at(0)-METPhi, EventWeight);
                h2_mjvsmj_2fatjet->Fill(mj->at(0), mj->at(1), EventWeight);
                for(int i=0; i<(int)mj->size(); i++) 
                { 
                    if(FatjetPt->at(i)<FatjetpTthres) continue;

                    h1_mj[2]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt[2]->Fill(FatjetPt->at(i), EventWeight);
                    h2_mjvsFatjetPt[2]->Fill(mj->at(i), FatjetPt->at(i), EventWeight);
                    // met matching 
                    if(getDPhi(FatjetPhi->at(i), METPhi)<1)
                    { 
                        h1_mj_metmatch[2]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_metmatch[2]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    else 
                    { 
                        h1_mj_notmetmatch[2]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmetmatch[2]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // b matching
                    float dRmin=999.;
                    for(int j=0; j<(int)JetPt->size(); j++) 
                    { 
                        if(JetCSV->at(j) < 0.679) continue;
                        float dRtmp = getDR(JetEta->at(j), RA4MusEta->at(0), JetPhi->at(j), RA4MusPhi->at(0)); 
                        if(dRtmp < dRmin) dRmin = dRtmp; 
                    }
                    if(dRmin<1.)
                    { 
                        h1_mj_bmatch[2]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_bmatch[2]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notbmatch[2]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notbmatch[2]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // mu matching
                    if(getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0))<1) 
                    { 
                        h1_mj_mumatch[2]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_mumatch[2]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notmumatch[2]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmumatch[2]->Fill(FatjetPt->at(i), EventWeight);
                    }

                }
                h1_Nfatjet[2]->Fill(Nfatjet_thres, EventWeight);
                h1_Nskinnyjet[2]->Fill(Nskinnyjet, EventWeight);

            }
            else if(Nfatjet_thres>0)
            {
                h1_muspT[1]->Fill(RA4MusPt->at(0), EventWeight);
                h1_musEta[1]->Fill(RA4MusEta->at(0), EventWeight);
                h1_musPhi[1]->Fill(RA4MusPhi->at(0), EventWeight);
                h1_mT[1]->Fill(mT, EventWeight);
                //h1_genpTttbar[1]->Fill(genpTttbar, EventWeight);
                h1_WpT[1]->Fill(WpT, EventWeight);
                h1_MJ[1]->Fill(MJ, EventWeight);
                h1_HT[1]->Fill(HT, EventWeight);
                h1_MET[1]->Fill(MET, EventWeight);
                h1_METPhi[1]->Fill(METPhi, EventWeight);
                h1_DPhi[1]->Fill(RA4MusPhi->at(0)-METPhi, EventWeight);
                for(int i=0; i<(int)mj->size(); i++) 
                { 
                    if(FatjetPt->at(i)<FatjetpTthres) continue;

                    h1_mj[1]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt[1]->Fill(FatjetPt->at(i), EventWeight);
                    h2_mjvsFatjetPt[1]->Fill(mj->at(i), FatjetPt->at(i), EventWeight);
                    // met matching 
                    if(getDPhi(FatjetPhi->at(i), METPhi)<1)
                    { 
                        h1_mj_metmatch[1]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_metmatch[1]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    else 
                    { 
                        h1_mj_notmetmatch[1]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmetmatch[1]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // b matching
                    float dRmin=999.;
                    for(int j=0; j<(int)JetPt->size(); j++) 
                    { 
                        if(JetCSV->at(j) < 0.679) continue;
                        float dRtmp = getDR(JetEta->at(j), RA4MusEta->at(0), JetPhi->at(j), RA4MusPhi->at(0)); 
                        if(dRtmp < dRmin) dRmin = dRtmp; 
                    }
                    if(dRmin<1.)
                    { 
                        h1_mj_bmatch[1]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_bmatch[1]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notbmatch[1]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notbmatch[1]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // mu matching
                    if(getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0))<1) 
                    { 
                        h1_mj_mumatch[1]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_mumatch[1]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notmumatch[1]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmumatch[1]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    
                }
                h1_Nfatjet[1]->Fill(Nfatjet_thres, EventWeight);
                h1_Nskinnyjet[1]->Fill(Nskinnyjet, EventWeight);
            }
            else if(Nfatjet_thres==0)
            {
                h1_muspT[0]->Fill(RA4MusPt->at(0), EventWeight);
                h1_musEta[0]->Fill(RA4MusEta->at(0), EventWeight);
                h1_musPhi[0]->Fill(RA4MusPhi->at(0), EventWeight);
                h1_mT[0]->Fill(mT, EventWeight);
                //h1_genpTttbar[0]->Fill(genpTttbar, EventWeight);
                h1_WpT[0]->Fill(WpT, EventWeight);
                h1_MJ[0]->Fill(MJ, EventWeight);
                h1_HT[0]->Fill(HT, EventWeight);
                h1_MET[0]->Fill(MET, EventWeight);
                h1_METPhi[0]->Fill(METPhi, EventWeight);
                h1_DPhi[0]->Fill(RA4MusPhi->at(0)-METPhi, EventWeight);
                for(int i=0; i<(int)mj->size(); i++) 
                { 
                    if(FatjetPt->at(i)<FatjetpTthres) continue;

                    h1_mj[0]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt[0]->Fill(FatjetPt->at(i), EventWeight);
                    h2_mjvsFatjetPt[0]->Fill(mj->at(i), FatjetPt->at(i), EventWeight);
                    // met matching 
                    if(getDPhi(FatjetPhi->at(i), METPhi)<1)
                    { 
                        h1_mj_metmatch[0]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_metmatch[0]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    else 
                    { 
                        h1_mj_notmetmatch[0]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmetmatch[0]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // b matching
                    float dRmin=999.;
                    for(int j=0; j<(int)JetPt->size(); j++) 
                    { 
                        if(JetCSV->at(j) < 0.679) continue;
                        float dRtmp = getDR(JetEta->at(j), RA4MusEta->at(0), JetPhi->at(j), RA4MusPhi->at(0)); 
                        if(dRtmp < dRmin) dRmin = dRtmp; 
                    }
                    if(dRmin<1.)
                    { 
                        h1_mj_bmatch[0]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_bmatch[0]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notbmatch[0]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notbmatch[0]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    // mu matching
                    if(getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0))<1) 
                    { 
                        h1_mj_mumatch[0]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_mumatch[0]->Fill(FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        h1_mj_notmumatch[0]->Fill(mj->at(i), EventWeight);
                        h1_FatjetPt_notmumatch[0]->Fill(FatjetPt->at(i), EventWeight);
                    }
                    
                }
                h1_Nfatjet[0]->Fill(Nfatjet_thres, EventWeight);
                h1_Nskinnyjet[0]->Fill(Nskinnyjet, EventWeight);
            }

            h1_muspT[5]->Fill(RA4MusPt->at(0), EventWeight);
            h1_musEta[5]->Fill(RA4MusEta->at(0), EventWeight);
            h1_musPhi[5]->Fill(RA4MusPhi->at(0), EventWeight);
            h1_mT[5]->Fill(mT, EventWeight);
            //h1_genpTttbar[5]->Fill(genpTttbar, EventWeight);
            h1_WpT[5]->Fill(WpT, EventWeight);
            h1_MJ[5]->Fill(MJ, EventWeight);
            h1_HT[5]->Fill(HT, EventWeight);
            h1_MET[5]->Fill(MET, EventWeight);
            h1_METPhi[5]->Fill(METPhi, EventWeight);
            h1_DPhi[5]->Fill(RA4MusPhi->at(0)-METPhi, EventWeight);
            for(int i=0; i<(int)mj->size(); i++) 
            { 
                if(FatjetPt->at(i)<FatjetpTthres) continue;

                h1_mj[5]->Fill(mj->at(i), EventWeight);
                h1_FatjetPt[5]->Fill(FatjetPt->at(i), EventWeight);
                h2_mjvsFatjetPt[5]->Fill(mj->at(i), FatjetPt->at(i), EventWeight);
                // met matching 
                if(getDPhi(FatjetPhi->at(i), METPhi)<1)
                { 
                    h1_mj_metmatch[5]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt_metmatch[5]->Fill(FatjetPt->at(i), EventWeight);
                }
                else 
                { 
                    h1_mj_notmetmatch[5]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt_notmetmatch[5]->Fill(FatjetPt->at(i), EventWeight);
                }
                // b matching
                float dRmin=999.;
                for(int j=0; j<(int)JetPt->size(); j++) 
                { 
                    if(JetCSV->at(j) < 0.679) continue;
                    float dRtmp = getDR(JetEta->at(j), RA4MusEta->at(0), JetPhi->at(j), RA4MusPhi->at(0)); 
                    if(dRtmp < dRmin) dRmin = dRtmp; 
                }
                if(dRmin<1.)
                { 
                    h1_mj_bmatch[5]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt_bmatch[5]->Fill(FatjetPt->at(i), EventWeight);
                } 
                else 
                { 
                    h1_mj_notbmatch[5]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt_notbmatch[5]->Fill(FatjetPt->at(i), EventWeight);
                }
                // mu matching
                if(getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0))<1) 
                { 
                    h1_mj_mumatch[5]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt_mumatch[5]->Fill(FatjetPt->at(i), EventWeight);
                } 
                else 
                { 
                    h1_mj_notmumatch[5]->Fill(mj->at(i), EventWeight);
                    h1_FatjetPt_notmumatch[5]->Fill(FatjetPt->at(i), EventWeight);
                }
            }
            h1_Nfatjet[5]->Fill(Nfatjet_thres, EventWeight);
            h1_Nskinnyjet[5]->Fill(Nskinnyjet, EventWeight);
        }
    }
    
    TString HistFileName = ch->GetTitle();
    HistFileName = Form("HistFiles/%s_pT%i.root", HistFileName.Data(), pTR0p5thres);
    cout << "[MJ Analysis] Writing " << HistFileName << endl;
    TFile *HistFile = new TFile(HistFileName, "RECREATE");
    gROOT->cd();
    HistFile->cd();
    // write histograms
    for(int i=0; i<6; i++)  
    {
        h1_MJ[i]->SetDirectory(0);                          h1_MJ[i]->Write();
        h1_HT[i]->SetDirectory(0);                          h1_HT[i]->Write();
        h1_MET[i]->SetDirectory(0);                         h1_MET[i]->Write();
        h1_METPhi[i]->SetDirectory(0);                      h1_METPhi[i]->Write();
        h1_FatjetPt[i]->SetDirectory(0);                    h1_FatjetPt[i]->Write();
        h1_FatjetPt_mumatch[i]->SetDirectory(0);            h1_FatjetPt_mumatch[i]->Write();
        h1_FatjetPt_notmumatch[i]->SetDirectory(0);         h1_FatjetPt_notmumatch[i]->Write();
        h1_FatjetPt_bmatch[i]->SetDirectory(0);             h1_FatjetPt_bmatch[i]->Write();
        h1_FatjetPt_notbmatch[i]->SetDirectory(0);          h1_FatjetPt_notbmatch[i]->Write();
        h1_FatjetPt_metmatch[i]->SetDirectory(0);           h1_FatjetPt_metmatch[i]->Write();
        h1_FatjetPt_notmetmatch[i]->SetDirectory(0);        h1_FatjetPt_notmetmatch[i]->Write();
        h1_DPhi[i]->SetDirectory(0);                        h1_DPhi[i]->Write();
        h1_mj[i]->SetDirectory(0);                          h1_mj[i]->Write();
        h1_mj_mumatch[i]->SetDirectory(0);                  h1_mj_mumatch[i]->Write();
        h1_mj_notmumatch[i]->SetDirectory(0);               h1_mj_notmumatch[i]->Write();
        h1_mj_notbmatch[i]->SetDirectory(0);                h1_mj_notbmatch[i]->Write();
        h1_mj_bmatch[i]->SetDirectory(0);                   h1_mj_bmatch[i]->Write();
        h1_mj_notmetmatch[i]->SetDirectory(0);              h1_mj_notmetmatch[i]->Write();
        h1_mj_metmatch[i]->SetDirectory(0);                 h1_mj_metmatch[i]->Write();
        h1_mT[i]->SetDirectory(0);                          h1_mT[i]->Write();
        //h1_genpTttbar[i]->SetDirectory(0);                  h1_genpTttbar[i]->Write();
        h1_WpT[i]->SetDirectory(0);                         h1_WpT[i]->Write();
        h1_muspT[i]->SetDirectory(0);                       h1_muspT[i]->Write();
        h1_musEta[i]->SetDirectory(0);                      h1_musEta[i]->Write();
        h1_musPhi[i]->SetDirectory(0);                      h1_musPhi[i]->Write();
        h1_Nfatjet[i]->SetDirectory(0);                     h1_Nfatjet[i]->Write();
        h1_Nskinnyjet[i]->SetDirectory(0);                  h1_Nskinnyjet[i]->Write();
        h2_mjvsFatjetPt[i]->SetDirectory(0);                h2_mjvsFatjetPt[i]->Write();
    }
    h2_mjvsmj_2fatjet->SetDirectory(0);                h2_mjvsmj_2fatjet->Write();
    HistFile->Close();

}

//
// Stacks
//
void DrawStack(TString HistName, int pTR0p5=30, int NMergeBins=1) 
{ 
    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_pT%i.root",pTR0p5));
   
    char *var; 
    if(HistName=="MET")                 	var=(char*)"MET [GeV]";
    if(HistName=="METPhi")              	var=(char*)"#phi(MET)";
    if(HistName=="FatjetPt")            	var=(char*)"pT(Fatjet) [GeV]";                   
    if(HistName=="FatjetPt_mumatch")    	var=(char*)"pT(Fatjet)(matched muon) [GeV]";    
    if(HistName=="FatjetPt_notmumatch") 	var=(char*)"pT(Fatjet)(not matched muon) [GeV]";
    if(HistName=="FatjetPt_bmatch")     	var=(char*)"pT(Fatjet)(matched b) [GeV]";        
    if(HistName=="FatjetPt_notbmatch")  	var=(char*)"pT(Fatjet)(not matched b) [GeV]";   
    if(HistName=="FatjetPt_metmatch")   	var=(char*)"pT(Fatjet)(matched met) [GeV]";    
    if(HistName=="FatjetPt_notmetmatch")	var=(char*)"pT(Fatjet)(not matched met) [GeV]"; 
    if(HistName=="DPhi")                	var=(char*)"#Delta(MET, muon)";                  
    if(HistName=="HT")                  	var=(char*)"H_{T} [GeV]";                        
    if(HistName=="MJ")                  	var=(char*)"M_{J} [GeV]";                        
    if(HistName=="mj")                  	var=(char*)"m_{j} [GeV]";                        
    if(HistName=="mj_mumatch")          	var=(char*)"m_{j}(matched muon) [GeV]";         
    if(HistName=="mj_notmumatch")       	var=(char*)"m_{j}(not matched muon) [GeV]";     
    if(HistName=="mj_bmatch")           	var=(char*)"m_{j}(matched b) [GeV]";             
    if(HistName=="mj_notbmatch")        	var=(char*)"m_{j}(not matched b) [GeV]";         
    if(HistName=="mj_metmatch")         	var=(char*)"m_{j}(matched met) [GeV]";           
    if(HistName=="mj_notmetmatch")      	var=(char*)"m_{j}(not matched met) [GeV]";       
    if(HistName=="mT")                  	var=(char*)"m_{T} [GeV]";                        
    if(HistName=="muspT")               	var=(char*)"p_{T}(muon) [GeV]";                  
    if(HistName=="musEta")              	var=(char*)"#eta(muon)";                         
    if(HistName=="musPhi")              	var=(char*)"#phi(muon)";                         
    if(HistName=="Nfatjet")             	var=(char*)"N_{fatjet}";
    if(HistName=="Nskinnyjet")          	var=(char*)"N_{skinny}";
    if(HistName=="WpT")                 	var=(char*)"p_{T}(W) [GeV]";
    //if(HistName=="genpTttbar")          	var=(char*)"gen p_{T}(ttbar) [GeV]";

    TH1F *h1_DATA[6], *h1_T[6], *h1_TT_sl[6], *h1_TT_ll[6], *h1_TT[6], *h1_WJets[6], *h1_DY[6], *h1_MC[6], *h1_Ratio[6]; 
    TH1F *h1_f1400_25[6], *h1_f1400_1000[6];
    TH1F *h1_One[6];
    THStack *st[6];
    TCanvas *c = new TCanvas("c","c",1200,360);  
    c->Divide(4,1);
    for(int i=2; i<6; i++) 
    {
        h1_DATA[i]      = (TH1F*)HistFile->Get(Form("h1_DATA_%s_%ifatjet", HistName.Data(), i)); 
        h1_T[i]         = (TH1F*)HistFile->Get(Form("h1_T_%s_%ifatjet", HistName.Data(), i));
        h1_TT_sl[i]     = (TH1F*)HistFile->Get(Form("h1_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h1_TT_ll[i]     = (TH1F*)HistFile->Get(Form("h1_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h1_WJets[i]     = (TH1F*)HistFile->Get(Form("h1_WJets_%s_%ifatjet", HistName.Data(), i));
//        h1_DY[i]      = (TH1F*)HistFile->Get(Form("h1_DY_%s_%ifatjet", HistName.Data(), i)); 
        h1_One[i]       = InitTH1F(Form("h1_One_%i",i), Form("h1_One_%i",i), 1,
                                    h1_DATA[i]->GetXaxis()->GetBinLowEdge(1),
                                    h1_DATA[i]->GetXaxis()->GetBinUpEdge(h1_DATA[i]->GetXaxis()->GetNbins())  );
        h1_One[i]->SetBinContent(1,1);
        h1_f1400_25[i]      = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1400_25_%s_%ifatjet", HistName.Data(), i)); 
        h1_f1400_1000[i]    = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1400_1000_%s_%ifatjet", HistName.Data(), i)); 

        // merge bins
        h1_DATA[i]->Rebin(NMergeBins);
        h1_T[i]->Rebin(NMergeBins);
        h1_TT_sl[i]->Rebin(NMergeBins);
        h1_TT_ll[i]->Rebin(NMergeBins);
        h1_WJets[i]->Rebin(NMergeBins);
//        h1_DY[i]->Rebin(NMergeBins);
        h1_f1400_25[i]->Rebin(NMergeBins);
        h1_f1400_1000[i]->Rebin(NMergeBins);
        
        h1_MC[i] = (TH1F*)h1_TT_sl[i]->Clone(Form("h1_MC_%s_%ifatjet", HistName.Data(), i));
        h1_MC[i]->Add(h1_TT_ll[i]);
        h1_MC[i]->Add(h1_WJets[i]);
        h1_MC[i]->Add(h1_T[i]);
//        h1_MC[i]->Add(h1_DY[i]);
        h1_Ratio[i] = (TH1F*)h1_DATA[i]->Clone(Form("h1_Ratio_%s_%ifatjet", HistName.Data(), i));
        h1_Ratio[i]->Divide(h1_MC[i]);

        h1_TT[i] = (TH1F*)h1_TT_ll[i]->Clone();
        h1_TT[i]->Add(h1_TT_sl[i]);

        h1cosmetic(h1_DATA[i],          Form("DATA %ifatjet", i),               kBlack, 1, 0,           var);
        h1cosmetic(h1_TT[i],            Form("TT %ifatjet", i),                 kBlack, 1, kYellow,     var);
        h1cosmetic(h1_TT_sl[i],         Form("TT(l) %ifatjet", i),              kBlack, 1, kYellow,     var);
        h1cosmetic(h1_TT_ll[i],         Form("TT(ll) %ifatjet", i),             kBlack, 1, kYellow-1,   var);
        h1cosmetic(h1_T[i],             Form("t+tW %ifatjet", i),               kBlack, 1, kYellow+1,   var);
        h1cosmetic(h1_WJets[i],         Form("WJets %ifatjet", i),              kBlack, 1, kGreen+2,    var);
        h1cosmetic(h1_Ratio[i],         Form(" "),                              kBlack, 1, kBlack,      var);
        h1cosmetic(h1_One[i],           Form(" "),                              kBlue,  1, 0,           var); 
        h1cosmetic(h1_f1400_25[i],      Form("T1tttt(1400,25) %ifatjet", i),    kRed,   1, 0,           var);
        h1cosmetic(h1_f1400_1000[i],    Form("T1tttt(1400,1000) %ifatjet", i),  kBlue,  1, 0,           var);

        c->cd(i-1);

        TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.04);
        pad1->SetRightMargin(0.07);
        pad1->SetLogy(0);
        pad1->Draw();
        pad1->cd();
        pad1->cd()->SetLogy(DoLog);
        TString StackTitle = Form("%i fatjets", i);
        if(i==5) StackTitle = "All fatjets";
        if(i==4) StackTitle = "4+ fatjets";
        st[i] = new THStack( Form("Stack %ifatjet", i), StackTitle);
//        st[i]->Add(h1_DY[i]);
        st[i]->Add(h1_WJets[i]);
        st[i]->Add(h1_T[i]); 
//         st[i]->Add(h1_TT[i]);
         st[i]->Add(h1_TT_ll[i]);
         st[i]->Add(h1_TT_sl[i]);
        st[i]->SetMaximum(h1_DATA[i]->GetMaximum()*(DoLog?50:1.6));
        st[i]->SetMinimum(0.1);
        st[i]->Draw("HIST"); 

        h1_DATA[i]->SetLineColor(kBlack);
        h1_DATA[i]->SetMarkerColor(kBlack);
        h1_DATA[i]->SetMarkerSize(1);
        h1_DATA[i]->SetMarkerStyle(20);
        h1_DATA[i]->SetStats(0);
        h1_DATA[i]->Draw("E SAME");

        TLegend *l1 = new TLegend(0.15, 0.65, 0.85, 0.85);
        l1->SetNColumns(2);
        l1->SetBorderSize(0);
        l1->SetFillColor(0);
        l1->SetFillStyle(0);
        l1->SetTextFont(42);
        l1->SetTextAlign(12);
        l1->SetTextSize(0.035);
        l1->SetFillColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->SetShadowColor(kWhite);
        l1->AddEntry(h1_DATA[i],        " Data",  "lp");
//        l1->AddEntry(h1_TT[i],          " TT",    "f");
        l1->AddEntry(h1_TT_sl[i],        " ",    "");
        l1->AddEntry(h1_TT_sl[i],        " TT(l)",    "f");
        l1->AddEntry(h1_TT_ll[i],        " TT(ll)",    "f");
        l1->AddEntry(h1_T[i],           " t+tW",  "f");
        l1->AddEntry(h1_WJets[i],       " WJets", "f");
//        l1->AddEntry(h1_DY[i],        " DY",    "f");
        l1->AddEntry(h1_f1400_25[i],    "T1tttt[1400,25]x1000",     "l");
        l1->AddEntry(h1_f1400_1000[i],  "T1tttt[1400,1000]x1000",   "l");
        l1->Draw();
        
        h1_f1400_25[i]->Draw("SAME HIST");
        h1_f1400_1000[i]->Draw("SAME HIST");

        c->cd(i-1);
        TPad *pad2 = new TPad("p_pull", "p_pull", 0.0, 0.0, 1.0, 0.3);
        pad2->Draw();
        pad2->cd();
        pad2->SetTopMargin(0.04);
        pad2->SetRightMargin(0.07);
        pad2->SetBottomMargin(0.4);

        h1_One[i]->SetLabelSize(0.1,"XY");
        h1_One[i]->SetTitleSize(0.1,"XY");
        h1_One[i]->SetTitleOffset(1.5);
        h1_One[i]->GetYaxis()->SetNdivisions(3,true);
        h1_One[i]->SetMinimum(0.5);
        h1_One[i]->SetMaximum(1.5);
        h1_One[i]->Draw("HIST");
        //h1_Ratio[i]->SetMinimum(0);
        //h1_Ratio[i]->SetMaximum(2);
        h1_Ratio[i]->SetMarkerStyle(20);
        h1_Ratio[i]->SetMarkerSize(0.5);
        //h1_Ratio[i]->SetLabelSize(0.1,"XY");
        //h1_Ratio[i]->SetTitleSize(0.1,"XY");
        //h1_Ratio[i]->SetTitleOffset(1.5);
        h1_Ratio[i]->Draw("SAME E X0");
    }

    // 
    if(HistName=="mj") HistName="JetMass";
    c->Print( Form("figures_FatjetpT%i/CompareDataMC_%s_pT%i%s%s.pdf", FatjetpTthres,
              HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", DoLog?"_log":"") ); 
    
    // 
    HistFile->Close();
    delete c;
}

//
// 2D histograms 
//
void Draw2D(TString HistName, int pTR0p5=30, int NMergeBins=1) 
{ 
    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_pT%i.root",pTR0p5));
   
    char *Xvar, *Yvar, *Zvar; 
    if(HistName=="mjvsmj")
    {   
        Xvar=(char*)"mj [GeV]";
        Yvar=(char*)"mj [GeV]";
        Zvar=(char*)"Events/bin";
    }
    if(HistName=="mjvsFatjetPt")
    {   
        Xvar=(char*)"mj [GeV]";
        Yvar=(char*)"pT(fatjet) [GeV]";
        Zvar=(char*)"Events/bin";
    }


    TH2F *h2_DATA[6], *h2_T[6], *h2_TT_sl[6], *h2_TT_ll[6], *h2_TT[6], *h2_WJets[6], *h2_DY[6]; 
    TH2F *h2_f1400_25[6], *h2_f1400_1000[6];
    TCanvas *c_TT_sl = new TCanvas("c_TT_sl","c_TT_sl",1600,360);  
    TCanvas *c_TT_ll = new TCanvas("c_TT_ll","c_TT_ll",1600,360);  
    c_TT_sl->Divide(4,1);
    c_TT_ll->Divide(4,1);
    for(int i=2; i<6; i++) 
    {
        if(HistName=="mjvsmj" && i!=2) continue;

        h2_DATA[i]          = (TH2F*)HistFile->Get(Form("h2_DATA_%s_%ifatjet", HistName.Data(), i)); 
        h2_T[i]             = (TH2F*)HistFile->Get(Form("h2_T_%s_%ifatjet", HistName.Data(), i));
        h2_TT_sl[i]         = (TH2F*)HistFile->Get(Form("h2_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h2_TT_ll[i]         = (TH2F*)HistFile->Get(Form("h2_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h2_WJets[i]         = (TH2F*)HistFile->Get(Form("h2_WJets_%s_%ifatjet", HistName.Data(), i));
//        h2_DY[i]          = (TH2F*)HistFile->Get(Form("h2_DY_%s_%ifatjet", HistName.Data(), i)); 
        h2_f1400_25[i]      = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1400_25_%s_%ifatjet", HistName.Data(), i)); 
        h2_f1400_1000[i]    = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1400_1000_%s_%ifatjet", HistName.Data(), i)); 

        h2cosmetic(h2_DATA[i],          Form("DATA %ifatjet", i),               Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_TT_sl[i],         Form("TT(sl) %ifatjet", i),             Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_TT_ll[i],         Form("TT(ll) %ifatjet", i),             Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_T[i],             Form("T %ifatjet", i),                  Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_WJets[i],         Form("WJets %ifatjet", i),              Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_f1400_25[i],      Form("T1tttt(1400,25) %ifatjet", i),    Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_f1400_1000[i],    Form("T1tttt(1400,1000) %ifatjet", i),  Xvar,   Yvar,   Zvar);
    }

    // Drawing
    for(int i=2; i<6; i++) 
    {
        if(HistName=="mjvsmj" && i!=2) continue;
        c_TT_sl->cd(i-1);
        h2_TT_sl[i]->Draw("colz");
        c_TT_ll->cd(i-1);
        h2_TT_ll[i]->Draw("colz");
    }
   
    c_TT_sl->Print( Form("figures_FatjetpT%i/2D_TT_sl_%s_pT%i%s%s.pdf", FatjetpTthres,
                HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", DoLog?"_log":"") ); 
    c_TT_ll->Print( Form("figures_FatjetpT%i/2D_TT_ll_%s_pT%i%s%s.pdf", FatjetpTthres,
                HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", DoLog?"_log":"") ); 
    
    // 
    HistFile->Close();
    delete c_TT_sl;
    delete c_TT_ll;
}


//
// main 
//
void DrawSL() 
{

    gROOT->ProcessLine(".L /Users/jaehyeok/macros/rootlogon.C");
    
    // ----------------------------------------
    //  Define chains  
    // ----------------------------------------
    TChain *ch_data         = new TChain("tree", "DATA");
    TChain *ch_ttbar_sl     = new TChain("tree", "TT_sl");
    TChain *ch_ttbar_ll     = new TChain("tree", "TT_ll");
    TChain *ch_wjets        = new TChain("tree", "WJets");
    TChain *ch_dy           = new TChain("tree", "DY");
    TChain *ch_t            = new TChain("tree", "T");
    TChain *ch_f1400_25     = new TChain("tree", "T1tttt_f1400_25");
    TChain *ch_f1400_1000   = new TChain("tree", "T1tttt_f1400_1000");
   
    TString BabyDir = "../../../babies/HT500MET100/";  
    // Data
    //ch_data->Add(BabyDir+"baby_SingleMu_Run2012A-13Jul2012-v1_f*.root");  // 808 pb if fully processed : 808 * 646/650 = 803 pb 
    //ch_data->Add(BabyDir+"baby_SingleMu_Run2012B-13Jul2012-v1_f*.root");  // 5237-808=4429 pb if fully processed : 4429 * 4247/4294 = 4380 pb  
                                                                            // 803 + 4380 = 5183
    //ch_data->Add(BabyDir+"baby_SingleMu_Run2012*.root");                    // Full SingleMu (19.1 fb-1)
    ch_data->Add(BabyDir+"baby_MuHad_*.root");                              // Full Muhad 
    
    // TT 
    ch_ttbar_sl->Add(BabyDir+"baby_TTJets_SemiLeptMGDecays_8TeV_f*.root");
    ch_ttbar_ll->Add(BabyDir+"baby_TTJets_FullLeptMGDecays_8TeV_f*.root");
    // WJets 
    ch_wjets->Add(BabyDir+"baby_WJetsToLNu_HT-400ToInf_8TeV_f*.root");
    // DY 
    ch_dy->Add(BabyDir+"baby_DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV_f*.root");
    // Singla top 
    ch_t->Add(BabyDir+"baby_*channel*_f*.root");

    // Signal
    ch_f1400_25->Add(BabyDir+"baby_*f1400_25.root");
    ch_f1400_1000->Add(BabyDir+"baby_*f1400_1000.root");
    
    
    // ----------------------------------------
    //  Get number of entries 
    // ----------------------------------------
    cout << "data               : " << ch_data->GetEntries()        << endl;
    cout << "ttbarl             : " << ch_ttbar_sl->GetEntries()    << endl;
    cout << "ttbarll            : " << ch_ttbar_ll->GetEntries()    << endl;
    cout << "wjets              : " << ch_wjets->GetEntries()       << endl;
    cout << "dy                 : " << ch_dy->GetEntries()          << endl;
    cout << "Single top         : " << ch_t->GetEntries()           << endl;
    cout << "T1tttt(1400,25)    : " << ch_f1400_25->GetEntries()    << endl;
    cout << "T1tttt(1400,1000)  : " << ch_f1400_1000->GetEntries()  << endl;
    
    if(!OnlyDraw) 
    {
        // ----------------------------------------
        //  Fill histrograms 
        // ----------------------------------------
        DoOneProcess(ch_data,	    pTR0p5thres); 
        DoOneProcess(ch_ttbar_sl,	pTR0p5thres); 
        DoOneProcess(ch_ttbar_ll,	pTR0p5thres); 
        DoOneProcess(ch_wjets,	    pTR0p5thres); 
        DoOneProcess(ch_dy,	        pTR0p5thres); 
        DoOneProcess(ch_t,	        pTR0p5thres); 
        DoOneProcess(ch_f1400_25,	pTR0p5thres);  
        DoOneProcess(ch_f1400_1000,	pTR0p5thres); 

        // ----------------------------------------
        //  Make the final histogram file
        // ----------------------------------------
        cout << "[MJ Analysis] Merging result files" << endl;
        gSystem->Exec(Form("rm HistFiles/Hist_pT%i.root",pTR0p5thres));
        gSystem->Exec(Form("hadd -f HistFiles/Hist_pT%i.root HistFiles/*_pT%i.root",pTR0p5thres,pTR0p5thres));
    }

    // ----------------------------------------
    //  Draw histograms 
    // ---------------------------------------- 
    DrawStack("muspT"           	);
    DrawStack("musPhi"          	);
    DrawStack("musEta"          	);
    DrawStack("mT"              	);
    DrawStack("mj"              	);
    DrawStack("mj_mumatch"      	);
    DrawStack("mj_notmumatch"   	);
    DrawStack("mj_bmatch"       	);
    DrawStack("mj_notbmatch"    	);
    DrawStack("mj_metmatch"       	);
    DrawStack("mj_notmetmatch"    	);
    DrawStack("MJ"         			);
    DrawStack("HT"         			);
    DrawStack("Nfatjet"    			);
    DrawStack("Nskinnyjet" 			);
    DrawStack("MET"        			);
    DrawStack("METPhi"     			);
    DrawStack("DPhi"       			);
    DrawStack("FatjetPt"   			);
    DrawStack("FatjetPt_mumatch"    );
    DrawStack("FatjetPt_notmumatch" );
    DrawStack("FatjetPt_bmatch"     );
    DrawStack("FatjetPt_notbmatch"  );
    DrawStack("FatjetPt_metmatch"   );
    DrawStack("FatjetPt_notmetmatch");
    DrawStack("WpT"                 );
    //DrawStack("genpTttbar"          );
    Draw2D("mjvsmj"                 );
    Draw2D("mjvsFatjetPt"           );
}

