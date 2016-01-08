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

void KappaOneConfig( float N1, float N2, float N3, float N4, float scale, bool doPoisson)  
{ 

    cout << N1 << " " << N2 << " " << N3 << " "  << N4 << " "  << scale 
         << endl;
    
    N1=N1*scale;
    N2=N2*scale;
    N3=N3*scale;
    N4=N4*scale;

    gStyle->SetOptStat(111111110);
    gStyle->SetStatW(0.3);                
    gStyle->SetStatH(0.25);                

    TRandom1 *fN1 = new TRandom1();
    TRandom2 *fN2 = new TRandom2();
    TRandom3 *fN3 = new TRandom3();
    TRandom3 *fN4 = new TRandom3(3241);
    TRandom3 rand(1234); // gamma

    TH1F *h1 = new TH1F("h1", "N1", 1000, (N1?N1-TMath::Sqrt(N1)*50:0)>0?N1-TMath::Sqrt(N1)*50:0, N1?TMath::Sqrt(N1)*10+N1:10);
    TH1F *h2 = new TH1F("h2", "N2", 1000, (N2?N2-TMath::Sqrt(N2)*50:0)>0?N2-TMath::Sqrt(N2)*50:0, N2?TMath::Sqrt(N2)*10+N2:10);
    TH1F *h3 = new TH1F("h3", "N3", 1000, (N3?N3-TMath::Sqrt(N3)*50:0)>0?N3-TMath::Sqrt(N3)*50:0, N3?TMath::Sqrt(N3)*10+N3:10);
    TH1F *h4 = new TH1F("h4", "N4", 1000, (N4?N4-TMath::Sqrt(N4)*50:0)>0?N4-TMath::Sqrt(N4)*50:0, N4?TMath::Sqrt(N4)*10+N4:10);
    //TH1F *h  = new TH1F("h",  "N1*N4/N2*N3", 100000, (N3*N2/N1 ? N3*N2/N1*(1-10/TMath::Sqrt(N3)) : 0)>0 ? N3*N2/N1*(1-10/TMath::Sqrt(N3)) : 0, 
    //                                             N3*N2/N1 ? N3*N2/N1*(1+10/TMath::Sqrt(N3)) : 10 );
    TH1F *h  = new TH1F("h",  "N1*N4/N2*N3", 100, 0, 5);

    TH1F *hband = (TH1F*)h->Clone("hband"); 

    for(int i=0; i<Npseudo; i++)  
    { 
        int checkpoint = (int)Npseudo/10;
        if((i%checkpoint)==0) cout << "Generated " << i << "/" << Npseudo << " experiments " << endl;

        float N1poisson   = N1;
        float N2poisson   = N2;
        float N3poisson   = N3;
        float N4poisson   = N4;
     
        if(doPoisson)
        {
            N1poisson   = fN1->Poisson(N1);
            N2poisson   = fN2->Poisson(N2);
            N3poisson   = fN3->Poisson(N3);
            N4poisson   = fN4->Poisson(N4);
        }
        else 
        {
            N1poisson = gsl_ran_gamma(N1+1,1,rand);
            N2poisson = gsl_ran_gamma(N2+1,1,rand);
            N3poisson = gsl_ran_gamma(N3+1,1,rand);
            N4poisson = gsl_ran_gamma(N4+1,1,rand);
        }

        float kappa = 999999.;
        if(N2poisson!=0 && N3poisson!=0) kappa = N1poisson * N4poisson / N2poisson / N3poisson;
        //float kappaData = N3temp;
        
        h1->Fill(N1poisson);
        h2->Fill(N2poisson);
        h3->Fill(N3poisson);
        h4->Fill(N4poisson);

        h->Fill(TMath::Min((Float_t)kappa,(Float_t)4.99));

    }
    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
    h3->Scale(1./h3->Integral());
    h4->Scale(1./h4->Integral());
    h->Scale(1./h->Integral());
  
    cout << "probability for kappa < 1 = "  << h->Integral(-1,20) << endl;
    cout << "probability for kappa > 1 = "  << h->Integral(21,999) << endl;

    //TLatex *tex_N = new TLatex(0.5,0.8,Form("N1=%.1f, N2=%.1f, N3=%.1f, N4=%.1f", N1, N2, N3, N4));
    TLatex *tex_N = new TLatex(0.5,0.8,Form("#frac{N1 #times N4}{N2 #times N3} = #frac{%.1f #times %.1f}{%.1f #times %.1f}", N1, N4, N2, N3));
    tex_N->SetNDC();
    tex_N->SetTextSize(0.08);
    tex_N->SetLineWidth(2);
    tex_N->SetTextAlign(22);
    
    TLatex *tex_problt1 = new TLatex(0.5,0.5,Form("P(x<1) = %.2f",h->Integral(-1,20)));
    tex_problt1->SetNDC();
    tex_problt1->SetTextSize(0.1);
    tex_problt1->SetLineWidth(2);
    tex_problt1->SetTextAlign(22);
    
    TLatex *tex_probgt1 = new TLatex(0.5,0.3,Form("P(x>1) = %.2f",h->Integral(21,999)));
    tex_probgt1->SetNDC();
    tex_probgt1->SetTextSize(0.1);
    tex_probgt1->SetLineWidth(2);
    tex_probgt1->SetTextAlign(22);

    TCanvas *c = new TCanvas("c","c",800,600);
    c->Divide(3,2);
    c->cd(1);
    h1->Draw();
    c->cd(2);
    h2->Draw();
    c->cd(3);
    h3->Draw();
    c->cd(4);
    h4->Draw();
    c->cd(5);
    c->cd(5)->SetLogy(0);
    h->Draw();
    c->cd(6); 
    tex_N->Draw();
    tex_problt1->Draw("same");
    tex_probgt1->Draw("same");

    c->Print(Form("fig/Kappa_Npseudo%i_N1%i_N2%i_N3%i_N4%i_%s.pdf", Npseudo, (int)N1, (int)N2, (int)N3, (int)N4, (doPoisson?"Possion":"Gamma")) );


}

//void Kappa(float N1=100, float N2=10, float N3=40, float N4=4, float scale=1)
void Kappa(float N1=100, float N2=30, float N3=10, float N4=3, float scale=1, bool doPoisson=false)
{ 
/*
    KappaOneConfig(N1,N2,N3,N4,0.25,    doPoisson);  
    KappaOneConfig(N1,N2,N3,N4,0.5,     doPoisson);  
    KappaOneConfig(N1,N2,N3,N4,0.75,    doPoisson);  
    KappaOneConfig(N1,N2,N3,N4,1,       doPoisson);  
    KappaOneConfig(N1,N2,N3,N4,2,       doPoisson);  
    KappaOneConfig(N1,N2,N3,N4,5,       doPoisson);  
    KappaOneConfig(N1,N2,N3,N4,10,      doPoisson);  
*/
    KappaOneConfig(10,10,10,10,10,      doPoisson);  
}
