#include <iostream>
#include <vector>

#include "TCanvas.h"
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
#include "TLorentzVector.h"
#include "TInterpreter.h"
#include "TColor.h"
#include "TStyle.h"

#include "../babytree_manuel.h"

using namespace std;


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
//TH1F initialization
//
TH1F* InitTH1F(const char* Name, const char* Title, int Nbins, double XMin, double XMax)
{
    TH1F *h1 = new TH1F(Name, Title, Nbins, XMin, XMax); 
    h1->Sumw2(); 
    return h1;
}
//
//TH2F initialization
//
TH2F* InitTH2F(const char* Name, const char* Title, int NXbins, double XMin, double XMax, int NYbins, double YMin, double YMax)
{
    TH2F *h2 = new TH2F(Name, Title, NXbins, XMin, XMax, NYbins, YMin, YMax);
    h2->Sumw2();
    return h2;
}

TH2F* InitTH2F(const char* Name, const char* Title, int NXbins, Float_t* xbins, int NYbins, Float_t* ybins){
    TH2F *h2 = new TH2F(Name, Title, NXbins, xbins, NYbins, ybins);
    h2->Sumw2();
    return h2;
}

//
// Fill TH1F
//
void FillTH1F(TH1F* &h1, double var, double weight=1)
{
    if(var >= h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())) 
        var=h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())-0.00001;
    if(var < h1->GetXaxis()->GetBinLowEdge(1)) 
        var=h1->GetXaxis()->GetBinLowEdge(1)+0.00001;
    h1->Fill(var, weight);        
}

//
// Fill TH2F
//
void FillTH2F(TH2F* &h2, float varX, float varY, float weight=1){
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

void h1cosmetic(TH1F* &h1, const char* title, int linecolor=kBlack, int linewidth=1, int fillcolor=0, TString var="", bool LOG=false)
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    if(LOG) h1->SetMinimum(0.1);
    else h1->SetMinimum(0.);
    h1->SetMaximum(h1->GetMaximum()*1.5);
}

//
// h2 cosmetics
//
void h2cosmetic(TH2F* &h2, const char* title, TString Xvar="", TString Yvar="", TString Zvar="")
{
    h2->SetTitle(title);
    h2->SetXTitle(Xvar);
    h2->SetYTitle(Yvar);
    h2->SetZTitle(Zvar);
    h2->SetStats(0);
}

bool goodrun(int run, int lumi)
{ 
    bool isgoodrun=false;
    if(run==251643 && ((lumi>=1 && lumi<=216) || (lumi>=222 && lumi<=606)))  isgoodrun=true;
    if(run==251721 && (lumi>=21 && lumi<=36))  isgoodrun=true;
    return isgoodrun;
}

bool goodevent()
{ 
    return (pass_hbhe_ && pass_cschalo_ && pass_eebadsc_ && pass_goodv_);
}

void RelVal() 
{
    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");  

    // Palette Color scheme -------------------------------------------
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    // ----------------------------------------------------------------

    // 
    TChain *ch_74   = new TChain("tree", "74X");
    TChain *ch_75   = new TChain("tree", "75X");
    //ch_74->Add("babies/small_quick_cfA_746p1_nleps1_trig0ON.root");                            
    //ch_75->Add("babies/small_quick_cfA_751_nleps1_trig0ON.root");                            
    //ch_74->Add("babies/small_quick_cfA_746p1_trig0ON.root");                            
    //ch_75->Add("babies/small_quick_cfA_751_trig0ON.root");                            
    ch_74->Add("babies/small_quick_cfA_746p1_v83_trig0ON_htg400_met_minig150.root");                            
    ch_75->Add("babies/small_quick_cfA_751_v83_trig0ON_htg400_met_minig150.root");                            
    
    // Variables and histograms  
    const int Nvar=26;
    const char*   var[Nvar];     
    int     nbins[Nvar];   
    float   begin[Nvar];   
    float   end[Nvar];     

    var[0]="MET";                               nbins[0]=20;    begin[0]=0;         end[0]=500;
    var[1]="MET #phi"/*"METPhi"*/;              nbins[1]=20;    begin[1]=-3.14;     end[1]=3.14;
    var[2]="HT";                                nbins[2]=20;    begin[2]=0;         end[2]=1500;
    var[3]="MJ";                                nbins[3]=20;    begin[3]=0;         end[3]=500;
    var[4]="m_{T}"/*mt*/;                       nbins[4]=20;    begin[4]=0;         end[4]=300;
    var[5]="n_{jets}"/*"NSkinnyjet"*/;          nbins[5]=10;    begin[5]=-0.5;      end[5]=9.5;
    var[6]="n_{bjets}"/*"NBtagCSVM"*/;          nbins[6]=5;     begin[6]=-0.5;      end[6]=4.5;
    var[7]="n_{fatjets}";                       nbins[7]=10;    begin[7]=-0.5;      end[7]=9.5;
    var[8]="p_{T}(jet1)"/*"JetPt1"*/;           nbins[8]=20;    begin[8]=0;         end[8]=1000;
    var[9]="p_{T}(fatjet1)"/*"FatjetPt1"*/;     nbins[9]=20;    begin[9]=0;         end[9]=1000;
    var[10]="p_{T}(electron)"/*"RA4ElsPt"*/;    nbins[10]=20;   begin[10]=0;        end[10]=300;
    var[11]="#eta(electron)"/*"RA4ElsEta"*/;    nbins[11]=20;   begin[11]=-2.5;     end[11]=2.5;
    var[12]="#phi(electron)"/*"RA4ElsPhi"*/;    nbins[12]=20;   begin[12]=-3.14;    end[12]=3.14;
    var[13]="p_{T}(muon)"/*"RA4MusPt"*/;        nbins[13]=20;   begin[13]=0;        end[13]=300;
    var[14]="#eta(muon)"/*"RA4MusEta"*/;        nbins[14]=20;   begin[14]=-2.5;     end[14]=2.5;
    var[15]="#phi(muon)"/*"RA4MusPhi"*/;        nbins[15]=20;   begin[15]=-3.14;    end[15]=3.14;
    var[16]="#eta(jet1)"/*"JetEta1"*/;          nbins[16]=20;   begin[16]=-2.5;     end[16]=2.5;
    var[17]="#phi(jet1)"/*"JetPhi1"*/;          nbins[17]=20;   begin[17]=-3.14;    end[17]=3.14;
    var[18]="#eta(fatjet1)"/*"FatjetEta1"*/;    nbins[18]=20;   begin[18]=-2.5;     end[18]=2.5;
    var[19]="#phi(fatjet1)"/*"FatjetPhi1"*/;    nbins[19]=20;   begin[19]=-3.14;    end[19]=3.14;
    var[20]="m(J1)"/*"mj1"*/;                   nbins[20]=20;   begin[20]=0;        end[20]=400;
    var[21]="m(jet1)"/*"JetM1"*/;               nbins[21]=20;   begin[21]=0;        end[21]=100;
    var[22]="miniiso(electron)"/*"RA4ElsMinIso"*/;  nbins[22]=20;   begin[22]=0;        end[22]=1;
    var[23]="miniiso(muon)"/*"RA4MusMinIso"*/;      nbins[23]=20;   begin[23]=0;        end[23]=1;
    var[24]="CSV"/*"JetCSV"*/;                      nbins[24]=20;   begin[24]=0;        end[24]=1;
    var[25]="m(J2,...)"/*"mJ"*/;                    nbins[25]=20;   begin[25]=0;        end[25]=200;

    //TH1F* h1_met_mini_74 = InitTH1F("h1_met_mini_74","h1_met_mini_74", 20, 0, 200);
    //TH1F* h1_met_mini_75 = InitTH1F("h1_met_mini_75","h1_met_mini_75", 50, 0, 500);

    TH1F *h1[Nvar][2];
    TH2F *h2[Nvar];
    TH1F *h1diff[Nvar];
    for(int ivar=0; ivar<Nvar; ivar++) 
    { 
            h1[ivar][0] = InitTH1F(Form("h1_%s_74",var[ivar]),Form("h1_%s_74",var[ivar]), nbins[ivar], begin[ivar], end[ivar]);
            h1[ivar][1] = InitTH1F(Form("h1_%s_75",var[ivar]),Form("h1_%s_75",var[ivar]), nbins[ivar], begin[ivar], end[ivar]);
            h2[ivar]    = InitTH2F(Form("h2_%s_74",var[ivar]),Form("h2_%s_74",var[ivar]), nbins[ivar], begin[ivar], end[ivar], nbins[ivar], begin[ivar], end[ivar]);
            if(ivar==5 || ivar==6 || ivar==7) h1diff[ivar] = InitTH1F(Form("h1diff_%s_75",var[ivar]),Form("h1diff_%s_75",var[ivar]), 5, -2.5, 2.5);
            else h1diff[ivar] = InitTH1F(Form("h1diff_%s_75",var[ivar]),Form("h1diff_%s_75",var[ivar]), nbins[ivar], -0.2*(end[ivar]-begin[ivar]), 0.2*(end[ivar]-begin[ivar]) );
    }


    // 
    InitBaby(ch_74); 
    InitBaby(ch_75); 

    //
    Int_t nentries_74 = (Int_t)ch_74->GetEntries();
    Int_t nentries_75 = (Int_t)ch_75->GetEntries(); 
    cout << "Nevt in 74 : " << nentries_74 << endl; 
    cout << "Nevt in 75 : " << nentries_75 << endl;  

    //--------------------------------------
    //  74x loop
    //--------------------------------------
    for(int i74 = 0; i74<nentries_74; i74++)
    {
        if((i74%100)==0) cout << "74x loop : " << i74 << " out of " << nentries_74 << endl;

        ch_74->GetEntry(i74);  
        
        if( !goodrun(run_,lumi_) || HT_<400 || MET_<150 || !goodevent() || NLeps_!=1) continue; 

        FillTH1F(h1[0][0], MET_);
        FillTH1F(h1[1][0], METPhi_);
        FillTH1F(h1[2][0], HT_);
        FillTH1F(h1[3][0], MJ_);
        FillTH1F(h1[4][0], mt_);
        FillTH1F(h1[5][0], Nskinnyjet_);
        FillTH1F(h1[6][0], NBtagCSVM_);
        FillTH1F(h1[7][0], FatjetPt_->size());
        FillTH1F(h1[8][0], JetPt_->at(0));
        FillTH1F(h1[9][0], FatjetPt_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[10][0], RA4ElsPt_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[11][0], RA4ElsEta_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[12][0], RA4ElsPhi_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[13][0], RA4MusPt_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[14][0], RA4MusEta_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[15][0], RA4MusPhi_->at(0));
        FillTH1F(h1[16][0], JetEta_->at(0));
        FillTH1F(h1[17][0], JetPhi_->at(0));
        FillTH1F(h1[18][0], FatjetEta_->at(0));
        FillTH1F(h1[19][0], FatjetPhi_->at(0));
        FillTH1F(h1[20][0], mj_->at(0)); 
        FillTH1F(h1[21][0], JetM_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[22][0], RA4ElsMinIso_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[23][0], RA4MusMinIso_->at(0));
        FillTH1F(h1[24][0], JetCSV_->at(0));
        for(int i=1; i<mj_->size(); i++) FillTH1F(h1[25][0], mj_->at(i)); 
    }
   

    //--------------------------------------
    //  74x loop
    //--------------------------------------
    for(int i75 = 0; i75<nentries_75; i75++)
    {
        if((i75%100)==0) cout << "75x loop : " << i75 << " out of " << nentries_75 << endl;

        ch_75->GetEntry(i75);  

        if( !goodrun(run_,lumi_) || HT_<400 || MET_<150 || !goodevent() || NLeps_!=1 ) continue; 
        //
        FillTH1F(h1[0][1], MET_);
        FillTH1F(h1[1][1], METPhi_);
        FillTH1F(h1[2][1], HT_);
        FillTH1F(h1[3][1], MJ_);
        FillTH1F(h1[4][1], mt_);
        FillTH1F(h1[5][1], Nskinnyjet_);
        FillTH1F(h1[6][1], NBtagCSVM_);
        FillTH1F(h1[7][1], FatjetPt_->size());
        FillTH1F(h1[8][1], JetPt_->at(0));
        FillTH1F(h1[9][1], FatjetPt_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[10][1], RA4ElsPt_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[11][1], RA4ElsEta_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[12][1], RA4ElsPhi_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[13][1], RA4MusPt_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[14][1], RA4MusEta_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[15][1], RA4MusPhi_->at(0));
        FillTH1F(h1[16][1], JetEta_->at(0));
        FillTH1F(h1[17][1], JetPhi_->at(0));
        FillTH1F(h1[18][1], FatjetEta_->at(0));
        FillTH1F(h1[19][1], FatjetPhi_->at(0));
        FillTH1F(h1[20][1], mj_->at(0));
        FillTH1F(h1[21][1], JetM_->at(0));
        if(RA4ElsPt_->size()>0) FillTH1F(h1[22][1], RA4ElsMinIso_->at(0));
        if(RA4MusPt_->size()>0) FillTH1F(h1[23][1], RA4MusMinIso_->at(0));
        FillTH1F(h1[24][1], JetCSV_->at(0));
        for(int i=1; i<mj_->size(); i++) FillTH1F(h1[25][1], mj_->at(i)); 
    }


    //--------------------------------------
    //  Double loops
    //--------------------------------------
    for(int i74 = 0; i74<nentries_74; i74++)
    //for(int i74 = 0; i74<50; i74++)
    {
        if((i74%100)==0) cout << i74 << " out of " << nentries_74 << endl;
        
        //
        // 74X
        //
        ch_74->GetEntry(i74);  
        if( !goodrun(run_,lumi_) || HT_<400 || MET_<150 || !goodevent() || NLeps_!=1) continue; 

        // Temp  
        int event_74x = event_;
        float MET_74x = MET_;
        float METPhi_74x = METPhi_;
        float HT_74x = HT_;
        float MJ_74x = MJ_;
        float mt_74x = mt_;
        int Nskinnyjet_74x = Nskinnyjet_;
        int NBtagCSVM_74x = NBtagCSVM_;
        int NFJ_74x = FatjetPt_->size();
        float JetPt1_74x = JetPt_->at(0);
        float FatjetPt1_74x = FatjetPt_->at(0);
        float RA4ElsPt_74x(-999), RA4ElsEta_74x(-999), RA4ElsPhi_74x(-999), RA4ElsMinIso_74x(-999);  
        if(RA4ElsPt_->size()>0) 
        {
            RA4ElsPt_74x = RA4ElsPt_->at(0);
            RA4ElsEta_74x = RA4ElsEta_->at(0);
            RA4ElsPhi_74x = RA4ElsPhi_->at(0);
            RA4ElsMinIso_74x = RA4ElsMinIso_->at(0);
        } 
        float RA4MusPt_74x(-999), RA4MusEta_74x(-999), RA4MusPhi_74x(-999), RA4MusMinIso_74x(-999);    
        if(RA4MusPt_->size()>0) 
        {
            RA4MusPt_74x = RA4MusPt_->at(0);
            RA4MusEta_74x = RA4MusEta_->at(0);
            RA4MusPhi_74x = RA4MusPhi_->at(0);
            RA4MusMinIso_74x = RA4MusMinIso_->at(0);
        } 
        float JetEta1_74x   = JetEta_->at(0);
        float JetPhi1_74x   = JetPhi_->at(0);
        float FatjetEta1_74x= FatjetEta_->at(0);
        float FatjetPhi1_74x= FatjetPhi_->at(0);
        float FatjetM1_74x  = mj_->at(0);
        float JetM1_74x     = JetM_->at(0);
        float JetCSV1_74x   = JetCSV_->at(0);

        //
        // 75X
        //
        for(int i75 = 0; i75<nentries_75; i75++)
        {

            ch_75->GetEntry(i75);  

            if( !goodrun(run_,lumi_) || HT_<400 || MET_<150 || !goodevent() || NLeps_!=1) continue; 
            if(event_74x != event_)  continue; // select the same event

            //if((JetPt1_74x-JetPt_->at(0))<-50)  cout << event_ << endl;
            //if((JetM1_74x-JetM_->at(0))<-10)    cout << event_ << endl; 
            //if((FatjetPt1_74x-FatjetPt_->at(0))<-100 || (FatjetPt1_74x-FatjetPt_->at(0))>100) cout << event_ << endl;
            //if((FatjetEta1_74x-FatjetEta_->at(0))<-0.3 || (FatjetEta1_74x-FatjetEta_->at(0))>0.3) cout << event_ << endl;
            //if((FatjetPhi1_74x-FatjetPhi_->at(0))>1) cout << event_ << endl;
            //if((FatjetM1_74x-mj_->at(0))<-20 || (FatjetM1_74x-mj_->at(0))>20) cout << event_ << endl; 
            if((MJ_74x-MJ_)<-30 || (MJ_74x-MJ_)>30) cout << event_ << endl; 
            /*
            //DumpEvent
            cout << "74 vs 75 " << endl;
            cout << "pT  : " << JetPt1_74x << " " << JetPt_->at(0) << endl;
            cout << "Eta : " << JetEta1_74x << " " << JetEta_->at(0) << endl;
            cout << "Phi : " << JetPhi1_74x << " " << JetPhi_->at(0) << endl;
            cout << "M   : " << JetM1_74x << " " << JetM_->at(0) << endl;
            cout << "F pT  : " << FatjetPt1_74x << " " << FatjetPt_->at(0) << endl;
            cout << "F Eta : " << FatjetEta1_74x << " " << FatjetEta_->at(0) << endl;
            cout << "F Phi : " << FatjetPhi1_74x << " " << FatjetPhi_->at(0) << endl;
            cout << "F M   : " << FatjetM1_74x << " " << mj_->at(0) << endl;
            */


            //        
            FillTH2F(h2[0], MET_74x,       MET_); 
            FillTH2F(h2[1], METPhi_74x,    METPhi_); 
            FillTH2F(h2[2], HT_74x,        HT_); 
            FillTH2F(h2[3], MJ_74x,        MJ_); 
            FillTH2F(h2[4], mt_74x,        mt_); 
            FillTH2F(h2[5], Nskinnyjet_74x, Nskinnyjet_); 
            FillTH2F(h2[6], NBtagCSVM_74x, NBtagCSVM_); 
            FillTH2F(h2[7], NFJ_74x,       FatjetPt_->size()); 
            if(getDR(JetEta1_74x,JetEta_->at(0),JetPhi1_74x,JetPhi_->at(0))<0.4) FillTH2F(h2[8], JetPt1_74x,   JetPt_->at(0)); 
            if(getDR(FatjetEta1_74x,FatjetEta_->at(0),FatjetPhi1_74x,FatjetPhi_->at(0))<1.0) FillTH2F(h2[9], FatjetPt1_74x,        FatjetPt_->at(0)); 
            if(RA4ElsPt_->size()>0) FillTH2F(h2[10], RA4ElsPt_74x,        RA4ElsPt_->at(0)); 
            if(RA4ElsPt_->size()>0) FillTH2F(h2[11], RA4ElsEta_74x,        RA4ElsEta_->at(0)); 
            if(RA4ElsPt_->size()>0) FillTH2F(h2[12], RA4ElsPhi_74x,        RA4ElsPhi_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH2F(h2[13], RA4MusPt_74x,        RA4MusPt_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH2F(h2[14], RA4MusEta_74x,        RA4MusEta_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH2F(h2[15], RA4MusPhi_74x,        RA4MusPhi_->at(0)); 
            FillTH2F(h2[16], JetEta1_74x,   JetEta_->at(0)); 
            FillTH2F(h2[17], JetPhi1_74x,   JetPhi_->at(0)); 
            FillTH2F(h2[18], FatjetEta1_74x,   FatjetEta_->at(0)); 
            FillTH2F(h2[19], FatjetPhi1_74x,   FatjetPhi_->at(0)); 
            if(getDR(FatjetEta1_74x,FatjetEta_->at(0),FatjetPhi1_74x,FatjetPhi_->at(0))<1.0) FillTH2F(h2[20], FatjetM1_74x,   mj_->at(0)); 
            if(getDR(JetEta1_74x,JetEta_->at(0),JetPhi1_74x,JetPhi_->at(0))<0.4) FillTH2F(h2[21], JetM1_74x,   JetM_->at(0)); 
            if(RA4ElsPt_->size()>0) FillTH2F(h2[22], RA4ElsMinIso_74x,        RA4ElsMinIso_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH2F(h2[23], RA4MusMinIso_74x,        RA4MusMinIso_->at(0)); 
            FillTH2F(h2[24], JetCSV1_74x,   JetCSV_->at(0)); 

            //
            FillTH1F(h1diff[0], MET_74x-MET_); 
            FillTH1F(h1diff[1], METPhi_74x-METPhi_); 
            FillTH1F(h1diff[2], HT_74x-HT_); 
            FillTH1F(h1diff[3], MJ_74x-MJ_); 
            FillTH1F(h1diff[4], mt_74x-mt_); 
            FillTH1F(h1diff[5], Nskinnyjet_74x-Nskinnyjet_); 
            FillTH1F(h1diff[6], NBtagCSVM_74x-NBtagCSVM_); 
            FillTH1F(h1diff[7], NFJ_74x-(int)FatjetPt_->size());  
            if(getDR(JetEta1_74x,JetEta_->at(0),JetPhi1_74x,JetPhi_->at(0))<0.4) FillTH1F(h1diff[8], JetPt1_74x-JetPt_->at(0)); 
            if(getDR(FatjetEta1_74x,FatjetEta_->at(0),FatjetPhi1_74x,FatjetPhi_->at(0))<1.0) FillTH1F(h1diff[9], FatjetPt1_74x-FatjetPt_->at(0));  
            if(RA4ElsPt_->size()>0) FillTH1F(h1diff[10], RA4ElsPt_74x-RA4ElsPt_->at(0)); 
            if(RA4ElsPt_->size()>0) FillTH1F(h1diff[11], RA4ElsEta_74x-RA4ElsEta_->at(0)); 
            if(RA4ElsPt_->size()>0) FillTH1F(h1diff[12], RA4ElsPhi_74x-RA4ElsPhi_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH1F(h1diff[13], RA4MusPt_74x-RA4MusPt_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH1F(h1diff[14], RA4MusEta_74x-RA4MusEta_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH1F(h1diff[15], RA4MusPhi_74x-RA4MusPhi_->at(0)); 
            FillTH1F(h1diff[16], JetEta1_74x-JetEta_->at(0)); 
            FillTH1F(h1diff[17], JetPhi1_74x-JetPhi_->at(0)); 
            FillTH1F(h1diff[18], FatjetEta1_74x-FatjetEta_->at(0)); 
            FillTH1F(h1diff[19], FatjetPhi1_74x-FatjetPhi_->at(0)); 
            if(getDR(FatjetEta1_74x,FatjetEta_->at(0),FatjetPhi1_74x,FatjetPhi_->at(0))<1.0) FillTH1F(h1diff[20], FatjetM1_74x-mj_->at(0)); 
            if(getDR(JetEta1_74x,JetEta_->at(0),JetPhi1_74x,JetPhi_->at(0))<0.4) FillTH1F(h1diff[21], JetM1_74x-JetM_->at(0)); 
            if(RA4ElsPt_->size()>0) FillTH1F(h1diff[22], RA4ElsMinIso_74x-RA4ElsMinIso_->at(0)); 
            if(RA4MusPt_->size()>0) FillTH1F(h1diff[23], RA4MusMinIso_74x-RA4MusMinIso_->at(0)); 
            FillTH1F(h1diff[24], JetCSV1_74x-JetCSV_->at(0)); 
        }
    } 
    
    //
    // Draw 
    //
    for(int ivar=0; ivar<Nvar; ivar++) 
    { 
        TLegend *l1 = new TLegend(0.4, 0.7, 0.87, 0.86);
        l1->SetBorderSize(0);
        l1->SetFillStyle(0);
        l1->SetFillColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->SetShadowColor(kWhite);
        l1->SetTextSize(0.06);
        l1->AddEntry(h1[ivar][0],        Form(" 74X(N=%.0f)", h1[ivar][0]->Integral()),  "l");
        l1->AddEntry(h1[ivar][1],        Form(" 75X(N=%.0f)", h1[ivar][1]->Integral()),  "l");
        l1->AddEntry(h1[ivar][1],        Form(" Common (N=%.0f)", h1diff[ivar]->Integral()),  "");

        h1cosmetic(h1[ivar][0], "Blue:74X, Red:75X", kBlue, 2, 0, var[ivar]);
        h1cosmetic(h1[ivar][1], "Blue:74X, Red:75X", kRed,  2, 0, var[ivar]);
        h2cosmetic(h2[ivar], var[ivar], Form("%s in 74X", var[ivar]), Form("%s in 75X",   var[ivar]), "Entries/bin");
        h1cosmetic(h1diff[ivar], "Difference", kBlue, 2, 0, Form("74X %s - 75X %s", var[ivar], var[ivar]), true);

        TH1F *h1_ratio = (TH1F*)h1[ivar][0]->Clone(); 
        h1_ratio->Divide(h1[ivar][1]);
        h1cosmetic(h1_ratio, "", kBlack, 1, 0, var[ivar]); 
        h1_ratio->SetLineColor(kBlack);
        h1_ratio->SetMarkerStyle(20);
        h1_ratio->SetMarkerSize(1); 
        h1_ratio->SetLabelSize(0.15);
        h1_ratio->SetLabelSize(0.15,"Y");
        h1_ratio->SetTitle("");
        h1_ratio->SetXTitle(var[ivar]);
        h1_ratio->GetXaxis()->SetTitleSize(0.15);
        h1_ratio->GetXaxis()->SetTitleOffset(1);
        h1_ratio->GetYaxis()->SetNdivisions(505);
        h1_ratio->SetMinimum(0);
        h1_ratio->SetMaximum(2);

        TCanvas *c = new TCanvas("c","c",1800,600);  
        c->Divide(3,1); 
        c->cd(1);  
        TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
        TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.33);
        pad1->SetBottomMargin(0.00001);
        pad1->SetBorderMode(0);
        pad2->SetTopMargin(0.00001);
        pad2->SetBottomMargin(0.3);
        pad2->SetBorderMode(0);
        pad1->Draw();
        pad2->Draw();
        pad1->cd();
        h1[ivar][0]->Draw("HIST");
        h1[ivar][1]->Draw("HIST SAME"); 
        l1->Draw("SAME");
        pad2->cd(); 
        h1_ratio->Draw("EP");
        
        c->cd(2);
        if(h2[ivar]->Integral()>0) h2[ivar]->Draw("COLZ");
        
        c->cd(3); 
        c->cd(3)->SetLogy(1); 
        if(h2[ivar]->Integral()>0) h1diff[ivar]->Draw("HIST"); 

        TString strvar = var[ivar];
        strvar.ReplaceAll("_","");  strvar.ReplaceAll(" ","");
        strvar.ReplaceAll("{","");  strvar.ReplaceAll("}","");
        strvar.ReplaceAll("(","");  strvar.ReplaceAll(")","");
        strvar.ReplaceAll("#","");  strvar.ReplaceAll("GeV","");

        c->Print(Form("fig/pdf/RelVal_%s.pdf", strvar.Data())); 
        c->Print(Form("fig/C/RelVal_%s.C", strvar.Data())); 
        delete l1;
        delete c;
        delete h1_ratio;
    }
}

