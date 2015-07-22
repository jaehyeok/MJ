#include <iostream>

#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "RooStats/RooStatsUtils.h"

using namespace std;

bool doGamma = true;
int Npseudo = 100000;

// Code from http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
// Parameter b could be theta...
double gsl_ran_gamma(const double a, const double b, TRandom3 &rand)
{
  
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

//
//
//
void h1cosmetic(TH1F* &h1, char* title, int linecolor=kBlack, int linewidth=1, int fillcolor=0, TString var="")
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.);
}

//
// 
//
void CompareGammaPoisson(float N=10)
{ 

    gStyle->SetOptStat(111111110);
    gStyle->SetStatW(0.3);                
    gStyle->SetStatH(0.25);                

    TRandom1 *frand = new TRandom1();
    TRandom3 rand(1234); // gamma

    int Nbins = (int)(N+TMath::Sqrt(N)*10);

    TH1F *hpoisson    = new TH1F("hpoisson",    "Poisson",   Nbins, -0.5,  Nbins-0.5);
    TH1F *hgamma      = new TH1F("hgamma",      "Gamma",     Nbins, -0.5,  Nbins-0.5);
    hpoisson->Sumw2();
    hgamma->Sumw2();

    for(int i=0; i<Npseudo; i++)  
    { 
        int checkpoint = (int)Npseudo/10;
        //if((i%checkpoint)==0) cout << "Generated " << i << "/" << Npseudo << " experiments " << endl;
        
        float Npoisson      = frand->Poisson(N);
        float Ngamma        = gsl_ran_gamma(N+1,1,rand);;

        hpoisson->Fill(Npoisson);
        hgamma->Fill(Ngamma);

    }
    hpoisson->Scale(1./hpoisson->Integral());
    hgamma->Scale(1./hgamma->Integral());

    h1cosmetic(hpoisson, Form("Poisson(black) and Gamma(red), N=%i",(int)N), kBlack, 3, 0, "");
    h1cosmetic(hgamma,   "Gamma", kRed, 3, 0, "");

    TCanvas *c = new TCanvas("c","c",400,400);
    c->cd(1);
    hpoisson->Draw("HIST");
    hgamma->Draw("HIST SAME");
    c->Print(Form("fig/ComparePoissonGamma_N%i.pdf", (int)N));
}


