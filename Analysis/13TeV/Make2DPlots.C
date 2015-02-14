#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

bool DoLog          = 1;
int SignalScale     = 1;

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
// 2D histograms 
//
void Make2DPlots(TString HistName, int NMergeBins=1, bool doRatio=false) 
{ 
    TFile* HistFile = TFile::Open("HistFiles/Hist.root");
   
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


    TH2F *h2_DATA[7], *h2_T[7], *h2_TT_sl[7], *h2_TT_ll[7], *h2_TT[7], *h2_WJets[7], *h2_DY[7]; 
    TH2F *h2_Ratio[7], *h2_MC[7]; 
    TH2F *h2_f1500_100[7], *h2_f1200_800[7];
    TCanvas *c_TT_sl = new TCanvas("c_TT_sl","c_TT_sl",1600,360);  
    TCanvas *c_TT_ll = new TCanvas("c_TT_ll","c_TT_ll",1600,360);  
    TCanvas *c_Ratio = new TCanvas("c_Ratio","c_Ratio",1600,360);  
    c_TT_sl->Divide(5,1);
    c_TT_ll->Divide(5,1);
    c_Ratio->Divide(5,1);
    for(int i=2; i<7; i++) 
    {
        if(HistName=="mjvsmj" && i!=2) continue;

        h2_DATA[i]          = (TH2F*)HistFile->Get(Form("h2_DATA_%s_%ifatjet", HistName.Data(), i)); 
        h2_T[i]             = (TH2F*)HistFile->Get(Form("h2_T_%s_%ifatjet", HistName.Data(), i));
        h2_TT_sl[i]         = (TH2F*)HistFile->Get(Form("h2_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h2_TT_ll[i]         = (TH2F*)HistFile->Get(Form("h2_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h2_WJets[i]         = (TH2F*)HistFile->Get(Form("h2_WJets_%s_%ifatjet", HistName.Data(), i));
//        h2_DY[i]          = (TH2F*)HistFile->Get(Form("h2_DY_%s_%ifatjet", HistName.Data(), i)); 
        h2_f1500_100[i]      = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1500_100_%s_%ifatjet", HistName.Data(), i)); 
        h2_f1200_800[i]    = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1200_800_%s_%ifatjet", HistName.Data(), i)); 

        h2cosmetic(h2_DATA[i],          Form("DATA %ifatjet", i),               Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_TT_sl[i],         Form("TT(l) %ifatjet", i),              Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_TT_ll[i],         Form("TT(ll) %ifatjet", i),             Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_T[i],             Form("T %ifatjet", i),                  Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_WJets[i],         Form("WJets %ifatjet", i),              Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_f1500_100[i],      Form("T1tttt(1200,25) %ifatjet", i),    Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_f1200_800[i],    Form("T1tttt(1200,1000) %ifatjet", i),  Xvar,   Yvar,   Zvar);
    }

    // Drawing
    for(int i=2; i<7; i++) 
    {
        if(HistName=="mjvsmj" && i!=2) continue;
        c_TT_sl->cd(i-1);
        h2_TT_sl[i]->Draw("colz");
        c_TT_ll->cd(i-1);
        h2_TT_ll[i]->Draw("colz");
    
        c_TT_sl->Print( Form("fig/2D_TT_sl_%s%s.pdf", 
                    HistName.Data(), DoLog?"_log":"") ); 
        c_TT_ll->Print( Form("fig/2D_TT_ll_%s%s.pdf", 
                    HistName.Data(), DoLog?"_log":"") ); 
    }
  
    if(doRatio) 
    { 
        for(int i=2; i<7; i++) 
        {
            h2_Ratio[i] = (TH2F*) h2_DATA[i]->Clone(); 
            h2_MC[i]    = (TH2F*) h2_TT_sl[i]->Clone();
            h2_MC[i]->Add(h2_TT_ll[i]);
            h2_MC[i]->Add(h2_T[i]);
            h2_MC[i]->Add(h2_WJets[i]);
            h2_Ratio[i]->Divide(h2_MC[i]);
            c_Ratio->cd(i-1);
            h2_Ratio[i]->SetMinimum(0.5);
            h2_Ratio[i]->SetMaximum(1.5);
            h2_Ratio[i]->Draw("colz");
            c_Ratio->Print( Form("fig/2D_Ratio_%s%s.pdf",
                        HistName.Data(), DoLog?"_log":"") ); 
        } 
    }

    // 
    HistFile->Close();
    delete c_TT_sl;
    delete c_TT_ll;
    delete c_Ratio;
}
