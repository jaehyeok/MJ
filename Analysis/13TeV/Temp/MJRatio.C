
#include <iostream>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TString.h"
#include "TChain.h"
#include "TMarker.h"
#include "TPaletteAxis.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TLorentzVector.h"
#include "TInterpreter.h"

#include "../babytree.h"
#include "../PassSelection.h"

using namespace std;
int pTR0p5          =30;
int FatjetpTthres   =50;
int mjthres         =0;

void myText(Double_t x,Double_t y, const char *text, Color_t red,float tsize) {

    //Double_t tsize=0.05;
    TLatex l; //l.SetTextAlign(12);
    l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextColor(red);
    l.DrawLatex(x,y,text);

}
//
// h1 cosmetics
//
void h1cosmetic(TH1F* &h1, TString title, int linecolor=kBlack, int linewidth=1, int fillcolor=0, TString var="")
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    //h1->SetMinimum(0.1);
}
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

float getDEta(float eta1, float eta2)
{
    return abs(eta1-eta2);
}
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

TH1F* GetProfile(TH2F* h2) 
{ 
    int Nbiny = h2->GetYaxis()->GetNbins(); 
    int Nbinx = h2->GetXaxis()->GetNbins(); 
    TH1F *h1slice = new TH1F("h1slice","h1slice",Nbiny,0,1);
    TH1F *h1      = new TH1F("h1","h1",Nbinx,0,1000);
    
    for(int x=1; x<=h2->GetXaxis()->GetNbins(); x++) 
    { 
        for(int y=1; y<=h2->GetYaxis()->GetNbins(); y++) 
        { 
            h1slice->SetBinContent(y, h2->GetBinContent(x,y));
        }
        cout << h1slice->GetMean() << " +/- " << h1slice->GetMeanError() 
             << " +/- " << h1slice->GetRMS() << endl;
        h1->SetBinContent(x, h1slice->GetMean());
        h1->SetBinError(x, h1slice->GetRMS());
        h1slice->Reset();
    }
    h1->SetMinimum(0.);
    h1->SetMaximum(1.);
    return h1;
}

void MJRatio(bool sig=true, const char* Region="Baseline")
{

    //gInterpreter->ExecuteMacro("~/macros/rootlogon.C");
    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");
    
    // chain      
    TChain *ch        = new TChain("tree");
    TString BabyDir = "/Users/jaehyeok/Research/Tools/fastjet-3.0.6/example/babies/13TeV/Phys14_JetPt20_30Mar2015/";  
    if(!sig) ch->Add(BabyDir+"baby_TTJets*.root");
    if(sig) ch->Add(BabyDir+"baby_*1500*.root");

    TH2F *h2_MJvsMj =  InitTH2F("h2_MJvsMj", "h2_MJvsMj", 10, 0, 1000, 40, 0, 1);

    InitBaby(ch); 
    
    Int_t nentries = (Int_t)ch->GetEntries();
    for(int i = 0; i<nentries; i++)
    {
        ch->GetEntry(i); 

        //if( event_!=111296860) continue; // FIXME   

        if(!PassNLep(1))  continue; // need this upfront because of mT calculation
        
        //
        // HT, Nb and Nskinny with selected jets
        //
        float HT_thres=0.;
        int   Nskinny_thres=0.;
        int   Nbcsvm_thres=0.;
        float Mj=0;
        for(unsigned int ijet=0; ijet<JetPt_->size(); ijet++)
        {
        
            if(JetPt_->at(ijet)<30 || TMath::Abs(JetEta_->at(ijet))>5) continue;
            float JetPz = JetPt_->at(ijet)/TMath::Tan(2*TMath::ATan(TMath::Exp(-JetEta_->at(ijet))));

            Mj = Mj + TMath::Sqrt(JetE_->at(ijet)*JetE_->at(ijet)-JetPt_->at(ijet)*JetPt_->at(ijet)-JetPz*JetPz);
            if(JetPt_->at(ijet)<40 || TMath::Abs(JetEta_->at(ijet))>2.5) continue;
            HT_thres=HT_thres+JetPt_->at(ijet);
            Nskinny_thres++;

            if( JetCSV_->at(ijet)<0.814 ) continue;
            Nbcsvm_thres++;
        }

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
        }
        if(RA4ElsPt_->size()==1) {
            mT  = TMath::Sqrt( 2*MET_*RA4ElsPt_->at(0)*(1-TMath::Cos(METPhi_-RA4ElsPhi_->at(0))) ); 
        }
       
        if( !PassBaselineSelection(HT_thres, MET_, Nbcsvm_thres, Nskinny_thres) ) continue; 
        if( !PassSelection(Region, HT_thres, MET_, Nbcsvm_thres, Nskinny_thres, mT, MJ_thres)) continue;

       
        //cout << MJ_thres << " : " << Mj << endl;
        FillTH2F(h2_MJvsMj, MJ_thres, Mj/MJ_thres, EventWeight_*5000);

        //if(Nfatjet_thres<2) continue;

        //c->Print(savename+".pdf");
        //c->Print(savename+".C");

    }
    TString title = sig?"T1tttt [1500,100]":"t#bar{t}";
    TString outfile = sig?"sig":"ttbar";
    
    float textSize = 0.032;
    TLatex *tex_title = new TLatex(0.5,0.95,title);
    tex_title->SetNDC();
    tex_title->SetTextSize(textSize+0.02);
    tex_title->SetLineWidth(2);
    tex_title->SetTextAlign(22);

    TCanvas *c2 = new TCanvas("c","c",400,400);
    c2->cd(1);
    h2_MJvsMj->RebinY(2);
    h2_MJvsMj->SetXTitle("M_{J}");
    h2_MJvsMj->SetYTitle("<M_{j}/M_{J}>");
    h2_MJvsMj->Draw("colz");
    tex_title->Draw("same");
    //c2->Print(outfile+"_2D.pdf");
    
    TH1F *h1_profile = GetProfile(h2_MJvsMj);
/* 
    TCanvas *c = new TCanvas("c","c",600,400);
    c->cd(1);

    h1cosmetic(h1_profile,          title,    kBlack, 2, 0,           "M_{J}");
    h1_profile->SetMarkerStyle(1); 
    h1_profile->SetYTitle("<M_{j}/M_{J}>");
    
    h1_profile->Draw("e");
    tex_title->Draw("same");
    c->Print(outfile+"_1D.pdf");
*/ 
    cout << h2_MJvsMj->Integral() << endl;
}
