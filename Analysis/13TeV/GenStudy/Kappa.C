#include <iostream>

#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "RooStats/RooStatsUtils.h"

using namespace std;

bool doGamma = true;

float Min(float a, float b) { return a <= b ? a : b; }

// Code from http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
// Parameter b could be theta...
double gsl_ran_gamma(const double a, const double b, TRandom3 &rand){
  
  if (a < 1)
  {
    double u = rand.Uniform(1);
    return gsl_ran_gamma(1.0 + a, b, rand) * pow (u, 1.0 / a);
  }

  double x, v, u;
  double d = a - 1.0 / 3.0;
  double c = (1.0 / 3.0) / sqrt (d);
  
  while (1) 
  {
    do 
    {
      x = rand.Gaus(0, 1.0);
      v = 1.0 + c * x;
    }
    while (v <= 0);
      
    v = v * v * v;
    u = rand.Uniform(1);

    if (u < 1 - 0.0331 * x * x * x * x) 
      break;

    if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
      break;
  }
    
  return b * d * v;
}

float GetKappa(TH2F *h2)
{ 
    // kappa = (N4/N3)/(N2/N1) = N4 * N1 / N2 * N3
    float kappa = h2->GetBinContent(2,2) * h2->GetBinContent(1,1) / 
                (h2->GetBinContent(1,2) * h2->GetBinContent(2,1));
    return kappa;
}

//
//
//
TH2F* KappaOneProcess(TChain *ch, float HTcut=400, float METcut=200, int njetscut=4)  
{ 
    TH2F *h2 = new TH2F("h2", "h2", 2, 0, 800, 2, 0, 280);

    float   MET_ = 0;
    float   HT_ = 0;
    float   MJ_ = 0;
    float   mT_ = 0;
    int   njets_ = 0;
    ch->SetBranchAddress("MET",     &MET_);
    ch->SetBranchAddress("HT",      &HT_);
    ch->SetBranchAddress("MJ",      &MJ_);
    ch->SetBranchAddress("mT",      &mT_);
    ch->SetBranchAddress("njets",   &njets_);

    int nentries = ch->GetEntries();
    for(int i=0; i<nentries; i++)  
    { 
        if(i%1000000==0) cout << i << "/" << nentries << endl;
        
        ch->GetEntry(i);        
     
        if(HT_<HTcut || MET_<METcut || njets_<njetscut) continue;

        h2->Fill(Min(MJ_,799.99),Min(mT_,279.99));

    }

    return h2;
}

void Kappa(float HTcut=400, float METcut=200, int njetscut=4)  
{ 

    gStyle->SetPaintTextFormat("6.3f");
    gStyle->SetOptStat(0);

    TChain *ch_tt_sl = new TChain("t"); // 2249 x 50k events
    ch_tt_sl->Add("/Users/jaehyeok/scratch/11Aug2015/small_TTJets_Semi*root");
    TChain *ch_tt_ll = new TChain("t"); // 995 x 50k events
    ch_tt_ll->Add("/Users/jaehyeok/scratch/11Aug2015/small_TTJets_Full*root");

    //
    // with Baseline selection at 8 TeV analysis 
    // N(tt_sl) = 308.5 +/- 4.4
    // N(tt_ll) = 58.5 +/- 1.4 
    // Relative weights calculated with
    //   if(HT_<400 || MET_<200 || njets_<4) continue;
    // so that the yields become above values
    TH2F *h2_tt_sl; 
    TH2F *h2_tt_ll; 
    TH2F *h2_tt; 

    h2_tt_sl = KappaOneProcess(ch_tt_sl, HTcut, METcut, njetscut);
    h2_tt_ll = KappaOneProcess(ch_tt_ll, HTcut, METcut, njetscut);
    
    h2_tt_sl->Sumw2();
    h2_tt_ll->Sumw2();

    h2_tt_sl->Scale(0.000326812);
    h2_tt_ll->Scale(0.000145449); 

    h2_tt = (TH2F*) h2_tt_sl->Clone("h2_tt");
    h2_tt->Add(h2_tt_ll); 
    h2_tt->Sumw2();

    // Get R and kappa
    float kappa = GetKappa(h2_tt);
    cout << "R21 = " << h2_tt->GetBinContent(2,1) / h2_tt->GetBinContent(1,1) << endl; 
    cout << "R43 = " << h2_tt->GetBinContent(2,2) / h2_tt->GetBinContent(1,2) << endl;  
    cout << "Kappa = " << kappa << endl;

    TLatex *tex_cuts = new TLatex(0.2,0.8,Form("HT>%i MET>%i njets#geq%i",(int)HTcut,(int)METcut,(int)njetscut));
    tex_cuts->SetNDC();
    tex_cuts->SetTextSize(0.05);
    tex_cuts->SetLineWidth(2);
    tex_cuts->SetTextAlign(12);
    
    TLatex *tex_R21 = new TLatex(0.2,0.6,Form("R21 = %.3f", h2_tt->GetBinContent(2,1) / h2_tt->GetBinContent(1,1)));
    tex_R21->SetNDC();
    tex_R21->SetTextSize(0.08);
    tex_R21->SetLineWidth(2);
    tex_R21->SetTextAlign(12);
    
    TLatex *tex_R43 = new TLatex(0.2,0.5,Form("R43 = %.3f", h2_tt->GetBinContent(2,2) / h2_tt->GetBinContent(1,2)));
    tex_R43->SetNDC();
    tex_R43->SetTextSize(0.08);
    tex_R43->SetLineWidth(2);
    tex_R43->SetTextAlign(12);
    
    TLatex *tex_kappa = new TLatex(0.2,0.4,Form("kappa = %.3f", kappa));
    tex_kappa->SetNDC();
    tex_kappa->SetTextSize(0.08);
    tex_kappa->SetLineWidth(2);
    tex_kappa->SetTextAlign(12);

    TCanvas *c2d = new TCanvas("c2d", "c2d", 800, 800);  
    c2d->Divide(2,2);
    c2d->cd(1);
    h2_tt_sl->SetTitle("ttbar(1l)");
    h2_tt_sl->SetMarkerSize(3);
    h2_tt_sl->Draw("colz text e");
    c2d->cd(2);
    h2_tt_ll->SetTitle("ttbar(2l)");
    h2_tt_ll->SetMarkerSize(3);
    h2_tt_ll->Draw("colz text e");
    c2d->cd(3);
    h2_tt->SetTitle("ttbar(1l+2l)");
    h2_tt->SetMarkerSize(3);
    h2_tt->Draw("colz text e");
    c2d->cd(4);
    tex_cuts->Draw();
    tex_R21->Draw("same");
    tex_R43->Draw("same");
    tex_kappa->Draw("same");
    c2d->Print(Form("fig/2DMJvsmT_HT%i_MET%i_njets%i.pdf",(int)HTcut,(int)METcut,(int)njetscut));
}
