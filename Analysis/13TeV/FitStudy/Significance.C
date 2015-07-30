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
int Npseudo = 30000;
int NpseudoGamma = 30000;

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


bool IsWithinInterval(float sigma, int ipseudo, float N1true, float N2true, float N3true, int N1poisson, int N2poisson, int N3poisson) 
{ 
    if(N1poisson==0) return false;

    TH1F *h1gamma = new TH1F("h1gamma", "N1poisson", 1000, (N1poisson?N1poisson-TMath::Sqrt(N1poisson)*50:0)>0?N1poisson-TMath::Sqrt(N1poisson)*50:0, 
                        N1poisson?TMath::Sqrt(N1poisson)*10+N1poisson:10);
    TH1F *h2gamma = new TH1F("h2gamma", "N2poisson", 1000, (N2poisson?N2poisson-TMath::Sqrt(N2poisson)*50:0)>0?N2poisson-TMath::Sqrt(N2poisson)*50:0, 
                        N2poisson?TMath::Sqrt(N2poisson)*10+N2poisson:10);
    TH1F *h3gamma = new TH1F("h3gamma", "N3poisson", 1000, (N3poisson?N3poisson-TMath::Sqrt(N3poisson)*50:0)>0?N3poisson-TMath::Sqrt(N3poisson)*50:0, 
                        N3poisson?TMath::Sqrt(N3poisson)*10+N3poisson:10);

    TH1F *hgamma  = new TH1F("hgamma",  "N2poisson*N3poisson/N1poisson", 100000, 
                       (N3poisson*N2poisson/N1poisson ? N3poisson*N2poisson/N1poisson*(1-50/TMath::Sqrt(N3poisson)) : 0)>0 ? N3poisson*N2poisson/N1poisson*(1-50/TMath::Sqrt(N3poisson)) : 0, 
                                                 N3poisson*N2poisson/N1poisson ? N3poisson*N2poisson/N1poisson*(1+50/TMath::Sqrt(N3poisson)) : 10 );
   

//    cout << "IsWithin1Sig : " 
//         << N1true << " " << N2true << " " << N3true << " " 
//         << N1poisson << " " << N2poisson << " " << N3poisson << " " << endl; 

    TRandom3 rand(4321); // gamma
    
    for(int i=0; i<Npseudo; i++)  
    { 
        float N1temp = gsl_ran_gamma(N1poisson+1,1,rand);        
        float N2temp = gsl_ran_gamma(N2poisson+1,1,rand);        
        float N3temp = gsl_ran_gamma(N3poisson+1,1,rand);        

        h1gamma->Fill(N1temp);
        h2gamma->Fill(N2temp);
        h3gamma->Fill(N3temp);
        hgamma->Fill(N1temp==0?999999.:N3temp*N2temp/N1temp);

        //cout << N1temp << " " << N2temp << " " << N3temp << endl;
    }
    hgamma->Scale(1./hgamma->Integral());
    
    //
    // Get 1sigma band
    //
    float   kappa=999999.;
    if(N1poisson!=0) kappa = N3poisson * N2poisson / N1poisson;  
    int     ibinmod  = hgamma->FindBin(kappa); 
    //cout << "i bin mod " << ibinmod << endl; 
   
    float msigma=0;
    float psigma=0;
    int   imsigma=-1;
    int   ipsigma=-1;
    for(int ip=ibinmod; ip<=hgamma->GetXaxis()->GetNbins(); ip++) 
    { 
        psigma = psigma + hgamma->GetBinContent(ip); 
        //cout << psigma << " " << ipsigma << endl;
        if(psigma > 0.5-RooStats::SignificanceToPValue(sigma)) 
        { 
            ipsigma = ip;
            break;
        }
    }
    
    for(int im=ibinmod-1; im>0; im--) 
    { 
        if(im!=0) msigma = msigma + hgamma->GetBinContent(im); 
        //cout << msigma << " " << imsigma << endl;
        if(msigma > 0.5-RooStats::SignificanceToPValue(sigma)) 
        { 
            imsigma = im;
            break;
        }
    }
    
    if(imsigma==-1) 
    { 
        psigma = 0; 
        for(int i=1; i<=hgamma->GetXaxis()->GetNbins(); i++)  
        { 
            psigma = psigma + hgamma->GetBinContent(i); 
            if(psigma > 1-2*RooStats::SignificanceToPValue(sigma)) 
            { 
                ipsigma = i;
                break;
            }
        }
        msigma=0;
    }
   
    if(imsigma==-1) imsigma=1;
 /*   
    TCanvas *c = new TCanvas("c","c",800,600);
    c->Divide(2,2);
    c->cd(1);
    h2->Draw();
    c->cd(2);
    h3->Draw();
    c->cd(3);
    h1->Draw();
    c->cd(4);
    //h->Rebin(100);
    h->Draw();

    c->Print(Form("fig/Gamma_N2N3overN1_pseudo%i_N1poisson%i_N2poisson%i_N3poisson%i.pdf", ipseudo, N1poisson, N2poisson, N3poisson) );

    delete c;
*/
    float msigmagamma = hgamma->GetBinLowEdge(imsigma);
    float psigmagamma = hgamma->GetBinLowEdge(ipsigma+1);

    delete h1gamma; delete h2gamma; delete h3gamma; delete hgamma; 

//    cout << imsigma << " +" << ipsigma << endl;
//    cout << kappa << " -" << msigmagamma << " +" << psigmagamma << endl;

    if(N2true*N3true/N1true<psigmagamma &&  N2true*N3true/N1true>msigmagamma) return true; 
    else return false;
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


float GetSignificangeGamma(float Nb, float Nsig) 
{ 
    TRandom3 rand(4321); // gamma
    
    TH1F *hbgamma        = new TH1F("hbgamma",        "Nbgamma",       100000, -0.5,  Nb+TMath::Sqrt(Nb)*10);
    hbgamma->Sumw2();
    for(int i=0; i<NpseudoGamma; i++)   
    { 
        float Nbgamma = gsl_ran_gamma(Nb+1,1,rand);        
        hbgamma->Fill(Nbgamma);
    }
    hbgamma->Scale(1./hbgamma->Integral());
    
    // Find the bin where there is Nsig and calculate p-value 
    int ibin = hbgamma->FindBin(Nb+Nsig);
    float pvalue = hbgamma->Integral(ibin,1000);
    float significance = RooStats::PValueToSignificance(pvalue);

    //cout << "[Gamma] p-value with Gamma  : " << pvalue << endl;
    //cout << "[Gamma] significance(sigma) : " << significance << endl;

    delete hbgamma;

    return significance;
}

//
// 
//
void SignificanceOneConfig(float Nb, float Nsig)
{ 

    cout << Nb << " " << Nsig << endl;

    gStyle->SetOptStat(111111110);
    gStyle->SetStatW(0.3);                
    gStyle->SetStatH(0.25);                

    TRandom1 *fNb = new TRandom1();
    TRandom2 *fNsig = new TRandom2();
    TRandom3 rand(1234); // gamma

    TH1F *hb        = new TH1F("hb",        "Poisson(x|Nb)",       100000, (Nb-TMath::Sqrt(Nb)*10)<0?-0.5:Nb-TMath::Sqrt(Nb)*10,  Nb+TMath::Sqrt(Nb)*10);
    TH1F *hsignif   = new TH1F("hsignif",   "Significance",  10,     1,  3);
    hb->Sumw2();
    hsignif->Sumw2();

    for(int i=0; i<Npseudo; i++)  
    { 
        int checkpoint = (int)Npseudo/10;
        if((i%checkpoint)==0) cout << "Generated " << i << "/" << Npseudo << " experiments " << endl;

        float Nbpoisson     = Nb;
        float Nsigpoisson   = Nsig;
        
        Nbpoisson      = fNb->Poisson(Nb);
        Nsigpoisson    = fNsig->Poisson(Nsig);

        float signifgamma = GetSignificangeGamma(Nbpoisson, Nsig);

        hb->Fill(Nbpoisson);
        hsignif->Fill(signifgamma);

    }
    hb->Scale(1./hb->Integral());
   
    TLine *lineNsig;
    lineNsig = new TLine(Nsig+Nb, 0, Nsig+Nb, hb->GetMaximum()*1.05);
    lineNsig->SetLineStyle(2);
    lineNsig->SetLineWidth(3);
    lineNsig->SetLineColor(kRed);
    
    TLine *lineSig;
    lineSig = new TLine(2, 0, 2, hsignif->GetMaximum()*1.05);
    lineSig->SetLineStyle(2);
    lineSig->SetLineWidth(3);
    lineSig->SetLineColor(kRed);
    
    TCanvas *c = new TCanvas("c","c",800,400);
    c->Divide(2,1);
    c->cd(1);
    hb->Draw("hist");
    lineNsig->Draw("same");
    c->cd(2);
    hsignif->Draw("hist");
    lineSig->Draw("same");
    
    c->Print(Form("fig/Signif_Npseudo%i_NpseugoGamma%i_Nb%i_Nsig%i.pdf", Npseudo, NpseudoGamma, (int)Nb, (int)Nsig) );
    c->Print(Form("fig/Signif_Npseudo%i_NpseugoGamma%i_Nb%i_Nsig%i.C", Npseudo, NpseudoGamma, (int)Nb, (int)Nsig) );
 

    delete fNb;
    delete fNsig;
    delete hb;
    delete hsignif;
    delete lineNsig;
    delete lineSig;
    delete c;
}

//void Significance(float Nb=10, float Nsig=7)
void Significance(float Nb=100, float Nsig=28)
{ 
//    SignificanceOneConfig(9.8,7.2);  
    SignificanceOneConfig(99.1,20.9);  
    SignificanceOneConfig(799.4,57.6);  
}
