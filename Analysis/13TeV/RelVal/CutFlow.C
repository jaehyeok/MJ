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

void CutFlow_() 
{ 
    //
    TChain *ch_74   = new TChain("tree", "74X");
    TChain *ch_75   = new TChain("tree", "75X");
    ch_74->Add("babies/small_quick_cfA_746p1.root");                            
    ch_75->Add("babies/small_quick_cfA_751.root");                            
 
    //
    TH1F *h1_74 = new TH1F("h1_74","h1_74",1,-1000,1000); 
    TH1F *h1_75 = new TH1F("h1_75","h1_75",1,-1000,1000); 
    
    //
    // 74X 
    //
    cout << "Cutflow for 74X" << endl;
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock)","goff");
    cout << "|" << setw(30) << "Good run" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1","goff");
    cout << "|" << setw(30) << "Trigger(HT350 MET100)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1 && ht>400","goff");
    cout << "|" << setw(30) << "+ HT>400 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1 && met_mini>150 && ht>400","goff");
    cout << "|" << setw(30) << "+ MET>150 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1 && nleps==1 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nleps=1" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1 && nels==1 && nmus==0 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nels=1" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1 && nmus==1 && nels==0 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nmus=1" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1 && nleps==1 && ht>400 && met_mini>150 && njets>4","goff");
    cout << "|" << setw(30) << "+ njets>=5 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74","goodrun(run,lumiblock) && trig[0]==1 && nleps==1 && ht>400 && met_mini>150 && njets>4 && nbm>0","goff");
    cout << "|" << setw(30) << "+ nb>=1 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 

    //
    // 75X 
    //
    cout << "Cutflow for 75X" << endl;
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock)","goff");
    cout << "|" << setw(30) << "Good run" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1","goff");
    cout << "|" << setw(30) << "Trigger(HT350 MET100)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1 && ht>400","goff");
    cout << "|" << setw(30) << "+ HT>400 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1 && met_mini>150 && ht>400","goff");
    cout << "|" << setw(30) << "+ MET>150 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1 && nleps==1 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nleps=1" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1 && nels==1 && nmus==0 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nels=1" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1 && nmus==1 && nels==0 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nmus=1" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1 && nleps==1 && ht>400 && met_mini>150 && njets>4","goff");
    cout << "|" << setw(30) << "+ njets>=5 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75","goodrun(run,lumiblock) && trig[0]==1 && nleps==1 && ht>400 && met_mini>150 && njets>4 && nbm>0","goff");
    cout << "|" << setw(30) << "+ nb>=1 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 

}

void CutFlow() 
{ 
    //
    TChain *ch_74   = new TChain("tree", "74X");
    TChain *ch_75   = new TChain("tree", "75X");
    //ch_74->Add("babies/small_quick_cfA_746p1.root");                            
    //ch_75->Add("babies/small_quick_cfA_751.root");                            
    ch_74->Add("babies/small_quick_cfA_746p1_v83.root");                            
    ch_75->Add("babies/small_quick_cfA_751_v83.root");                            
 
    //
    TH1F *h1_74 = new TH1F("h1_74","h1_74",1,-1000,1000); 
    TH1F *h1_75 = new TH1F("h1_75","h1_75",1,-1000,1000); 
    
    //
    // 74X 
    // 
    TString goodevent = "(pass_hbhe && pass_cschalo && pass_eebadsc && pass_goodv)";
    TString nels1 = "Sum$((els_pt>20)&&(els_sigid && els_miniso<0.1))==1";
    TString nmus1 = "Sum$((mus_pt>20)&&(mus_sigid && mus_miniso<0.2))==1";
    TString nleps1 = "(Sum$((els_pt>20)&&(els_sigid && els_miniso<0.1))+Sum$((mus_pt>20)&&(mus_sigid && mus_miniso<0.2)))==1";
    TString nels1noiso = "Sum$((els_pt>20)&&(els_sigid && els_miniso<999.1))==1";
    TString nmus1noiso = "Sum$((mus_pt>20)&&(mus_sigid && mus_miniso<999.2))==1";
    TString nleps1noiso = "(Sum$((els_pt>20)&&(els_sigid && els_miniso<999.1))+Sum$((mus_pt>20)&&(mus_sigid && mus_miniso<999.2)))==1";
    cout << "Cutflow for 74X" << endl;
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock)","goff");
    cout << "|" << setw(30) << "Good run" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1","goff");
    cout << "|" << setw(30) << "Trigger(HT350 MET100)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && ht>400","goff");
    cout << "|" << setw(30) << "+ HT>400 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && met_mini>150 && ht>400","goff");
    cout << "|" << setw(30) << "+ MET>150 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nleps1+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nleps=1" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nels1noiso+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nels(without miniiso)=1" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nels1+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nels=1(with miniso)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nmus1noiso+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nmus=1(without miniiso)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nmus1+ "&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nmus=1(with miniso)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nleps1+"&& ht>400 && met_mini>150 && njets>4","goff");
    cout << "|" << setw(30) << "+ njets>=5 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nleps1+"&& ht>400 && met_mini>150 && njets>4 && nbm>0","goff");
    cout << "|" << setw(30) << "+ nb>=1 GeV" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset(); 

    //
    // 75X 
    //
    cout << "Cutflow for 75X" << endl;
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock)","goff");
    cout << "|" << setw(30) << "Good run" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1","goff");
    cout << "|" << setw(30) << "Trigger(HT350 MET100)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && ht>400","goff");
    cout << "|" << setw(30) << "+ HT>400 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && met_mini>150 && ht>400","goff");
    cout << "|" << setw(30) << "+ MET>150 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nleps1+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nleps=1" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nels1noiso+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nels(without miniiso)=1" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nels1+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nels=1(with miniso)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nmus1noiso+"&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nmus=1(without miniiso)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nmus1+ "&& ht>400 && met_mini>150","goff");
    cout << "|" << setw(30) << "+ nmus=1(with miniso)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nleps1+"&& ht>400 && met_mini>150 && njets>4","goff");
    cout << "|" << setw(30) << "+ njets>=5 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 &&"+nleps1+"&& ht>400 && met_mini>150 && njets>4 && nbm>0","goff");
    cout << "|" << setw(30) << "+ nb>=1 GeV" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset(); 


    //
    // Investigation on muons 
    //
    
    // 74X
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1");
    cout << "|" << setw(50) << "Good run + filter + Trigger(HT350 MET100) " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>10)==1","goff");
    cout << "|" << setw(50) << "+ 1 loose muon w/o iso (pT>10)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20)==1","goff");
    cout << "|" << setw(50) << "or 1 loose muon w/o iso (pT>20)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid)==1","goff");
    cout << "|" << setw(50) << "or 1 sig muon w/o iso (pT>20)" << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1","goff");
    cout << "|" << setw(50) << "or 1 signal muon w/ miniso<0.2 (pT>20) " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400","goff");
    cout << "|" << setw(50) << "+ HT>400 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>100","goff");
    cout << "|" << setw(50) << "+ MET>100 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>120","goff");
    cout << "|" << setw(50) << "+ MET>120 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>140","goff");
    cout << "|" << setw(50) << "+ MET>140 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(50) << "+ MET>150 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>160","goff");
    cout << "|" << setw(50) << "+ MET>160 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>180","goff");
    cout << "|" << setw(50) << "+ MET>180 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  
    ch_74->Draw("njets>>h1_74",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>200","goff");
    cout << "|" << setw(50) << "+ MET>200 " << " |"
                << setw(10) << h1_74->Integral() << " |" << endl;
    h1_74->Reset();  

    cout << endl;
    cout << endl;

    // 75X
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1");
    cout << "|" << setw(50) << "Good run + filter + Trigger(HT350 MET100) " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>10)==1","goff");
    cout << "|" << setw(50) << "+ 1 loose muon w/o iso (pT>10)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20)==1","goff");
    cout << "|" << setw(50) << "or 1 loose muon w/o iso (pT>20)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid)==1","goff");
    cout << "|" << setw(50) << "or 1 sig muon w/o iso (pT>20)" << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1","goff");
    cout << "|" << setw(50) << "or 1 signal muon w/ miniso<0.2 (pT>20) " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400","goff");
    cout << "|" << setw(50) << "+ HT>400 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>100","goff");
    cout << "|" << setw(50) << "+ MET>100 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>120","goff");
    cout << "|" << setw(50) << "+ MET>120 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>140","goff");
    cout << "|" << setw(50) << "+ MET>140 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>150","goff");
    cout << "|" << setw(50) << "+ MET>150 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>160","goff");
    cout << "|" << setw(50) << "+ MET>160 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>180","goff");
    cout << "|" << setw(50) << "+ MET>180 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  
    ch_75->Draw("njets>>h1_75",goodevent+" && goodrun(run,lumiblock) && trig[0]==1 && Sum$(mus_pt>20&&mus_sigid&&mus_miniso<0.2)==1 && ht>400 && met_mini>200","goff");
    cout << "|" << setw(50) << "+ MET>200 " << " |"
                << setw(10) << h1_75->Integral() << " |" << endl;
    h1_75->Reset();  

}
