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
int Npseudo = 5000;
int NpseudoGamma = 5000;

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

float GetRErr(float Num, float Den, float NumErr, float DenErr)
{
    float R = Num/Den;
    float RErr = R*TMath::Sqrt( DenErr*DenErr/Den/Den + NumErr*NumErr/Num/Num );
    if(Den==0) return -999;
    else return RErr;
}

float GetXErr(float N1, float N2, float N1Err, float N2Err)
{
    return TMath::Sqrt( N1*N1*N2Err*N2Err + N2*N2*N1Err*N1Err );
}


float GetLimitGamma(float Nb, float Nsig) 
{ 
    TRandom3 rand(4321); // gamma
    TRandom3 fNsig_(2121);
    
    TH1F *hbgamma           = new TH1F("hbgamma",        "Nbgamma",       10000, -0.5,  Nb+TMath::Sqrt(Nb)*10);
    TH1F *hsiggamma         = new TH1F("hsiggamma",      "Nsiggamma",   10000, -0.5,  Nb+Nsig+TMath::Sqrt(Nb+Nsig)*10);
    hbgamma->Sumw2();
    hsiggamma->Sumw2();

    for(int i=0; i<NpseudoGamma; i++)   
    { 
        float Nbfluc = gsl_ran_gamma(Nb+1,1,rand);        
        //float Nsigfluc = gsl_ran_gamma(Nsig+1,1,rand);   
        float Nsigfluc      = fNsig_.Poisson(Nsig);
        hbgamma->Fill(Nbfluc);
        hsiggamma->Fill(Nbfluc+Nsigfluc);
    }
    hbgamma->Scale(1./hbgamma->Integral());
    hsiggamma->Scale(1./hsiggamma->Integral());
    
    // Find the bin where there is Nb and calculate p-value 
    int ibinsig = hsiggamma->FindBin(Nb);
    float pvaluesig = hsiggamma->Integral(-1,ibinsig);
    int ibinbkg = hbgamma->FindBin(Nb);
    float pvaluebkg = hbgamma->Integral(-1,ibinbkg); 
    float CLs = pvaluesig / (1-pvaluebkg); 

    cout << "[Gamma] p-value(sig) with Gamma  : " << pvaluesig << endl;
    cout << "[Gamma] p-value(bkg) with Gamma  : " << pvaluebkg << endl;
    cout << "[Gamma] CLs with Gamma           : " << CLs << endl;
   
    if(0) 
    {
        TCanvas *c_gamma = new TCanvas("c_gamma","c_gamma",800,400);
        c_gamma->Divide(2,1);
        c_gamma->cd(1);
        hbgamma->Draw("hist");
        hsiggamma->Draw("same e");
        c_gamma->cd(2);
        //hlimit->Draw("hist");
        c_gamma->Print("fig/test.C");
    }

    delete hbgamma;
    delete hsiggamma;
    
    //return CLs;
    return pvaluesig;
}

//
// 
//
void LimitOneConfig(float Nb, float Nsig)
{ 

    cout << Nb << " " << Nsig << endl;

    gStyle->SetOptStat(111111110);
    gStyle->SetStatW(0.3);                
    gStyle->SetStatH(0.25);                

    TRandom1 *fNb = new TRandom1();
    TRandom2 *fNsig = new TRandom2();
    TRandom3 rand(1234); // gamma

    TH1F *hb        = new TH1F("hb",        "Nb",       10000, -0.5,  Nb+TMath::Sqrt(Nb)*10);
    TH1F *hsig      = new TH1F("hsig",      "Nsig",     10000, -0.5,  Nb+Nsig+TMath::Sqrt(Nb+Nsig)*10);
    TH1F *hlimit    = new TH1F("hlimit",     "Nlimit",  50,     0,  0.5);
    hb->Sumw2();
    hsig->Sumw2();
    hlimit->Sumw2();

    for(int i=0; i<Npseudo; i++)  
    { 
        int checkpoint = (int)Npseudo/10;
        //if((i%checkpoint)==0) cout << "Generated " << i << "/" << Npseudo << " experiments " << endl;
        
        float Nbpoisson      = fNb->Poisson(Nb);
        float Nsigpoisson    = fNsig->Poisson(Nsig);

        float signifgamma = GetLimitGamma(Nbpoisson, Nsig);

        hb->Fill(Nbpoisson);
        hsig->Fill(Nsigpoisson+Nbpoisson);
        hlimit->Fill(signifgamma);

    }
    hb->Scale(1./hb->Integral());
    hsig->Scale(1./hsig->Integral());
  
    // Find the bin where there is Nb and calculate p-value
    int ibinsig = hsig->FindBin(Nb);
    float pvaluesig = hsig->Integral(-1,ibinsig);
    int ibinbkg = hb->FindBin(Nb);
    float pvaluebkg = hb->Integral(-1,ibinbkg);
    float CLs = pvaluesig / (1-pvaluebkg);

    cout << "[Poisson] p-value(sig) with Poisson  : " << pvaluesig << endl;
    cout << "[Poisson] p-value(bkg) with Poisson  : " << pvaluebkg << endl;
    cout << "[Poisson] CLs with Poisson           : " << CLs << endl;

    TCanvas *c = new TCanvas("c","c",800,400);
    c->Divide(2,1);
    c->cd(1);
    hb->Draw("hist");
    hsig->Draw("same e");
    c->cd(2);
    hlimit->Draw("hist");
    c->Print(Form("fig/Limit_Npseudo%i_NpseugoGamma%i_Nb%i_Nsig%i.pdf", Npseudo, NpseudoGamma, (int)Nb, (int)Nsig) );
    c->Print(Form("fig/Limit_Npseudo%i_NpseugoGamma%i_Nb%i_Nsig%i.C", Npseudo, NpseudoGamma, (int)Nb, (int)Nsig) );
   
}

void Limit(float Nb=100, float Nsig=28)
{ 
    LimitOneConfig(20,10);  
}
