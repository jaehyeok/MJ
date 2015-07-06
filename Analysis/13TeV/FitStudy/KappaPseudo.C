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

using namespace std;

bool doGamma = true;

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


bool IsWithin1Sig(int ipseudo, float N1true, float N2true, float N3true, int N1poisson, int N2poisson, int N3poisson) 
{ 
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
    
    for(int i=0; i<50000; i++)  
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
    float   kappa    = N3poisson * N2poisson / N1poisson;  
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
        if(psigma > /*0.34*/0.4987/*0.475*/) 
        { 
            ipsigma = ip;
            break;
        }
    }
    
    for(int im=ibinmod-1; im>0; im--) 
    { 
        if(im!=0) msigma = msigma + hgamma->GetBinContent(im); 
        //cout << msigma << " " << imsigma << endl;
        if(msigma > /*0.34*/0.4987/*0.475*/) 
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
            if(psigma > /*0.34*2*/0.4987*2/*0.475*2*/) 
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

void KappaPseudoOneConfig( float N1, float N2, float N3,  
                           bool flucN1=false, bool flucN2=true, bool flucN3=true, 
                           float scale1=1., float scale2=1., float scale3=1.) 
{ 

    cout << N1 << " " << N2 << " " << N3 << " " 
         << flucN1 << " " << flucN2 << " " << flucN3 << " " 
         << scale1 << " " << scale2 << " " << scale3 << " " 
         << endl;

    gStyle->SetOptStat(111111110);
    gStyle->SetStatW(0.3);                
    gStyle->SetStatH(0.25);                

    int Npseudo = 50000;
    N1 = N1 * scale1;
    N2 = N2 * scale2;
    N3 = N3 * scale3;
    // Method1
    //float N1=6.1*scale1; 
    //float N2=4.3*scale2; 
    //float N3=0.8*scale3; 
    // Method2
    //float N1=32.7*scale1; 
    //float N2=7.8*scale2; 
    //float N3=3.7*scale3; 
    //float N4=0.4; 
    //float N1=100*scale1; 
    //float N2=100*scale2; 
    //float N3=100*scale3; 
    //float N1err=5.1;
    //float N2err=2.1;
    //float N3err=1.0;
    //float N4err=0.3;
    float N1err=TMath::Sqrt(N1);
    float N2err=TMath::Sqrt(N2);
    float N3err=TMath::Sqrt(N3);

    TRandom1 *fN1 = new TRandom1();
    TRandom2 *fN2 = new TRandom2();
    TRandom3 *fN3 = new TRandom3();
    TRandom3 rand(1234); // gamma

    TH1F *h1 = new TH1F("h1", "N1", 1000, (N1?N1-TMath::Sqrt(N1)*50:0)>0?N1-TMath::Sqrt(N1)*50:0, N1?TMath::Sqrt(N1)*10+N1:10);
    TH1F *h2 = new TH1F("h2", "N2", 1000, (N2?N2-TMath::Sqrt(N2)*50:0)>0?N2-TMath::Sqrt(N2)*50:0, N2?TMath::Sqrt(N2)*10+N2:10);
    TH1F *h3 = new TH1F("h3", "N3", 1000, (N3?N3-TMath::Sqrt(N3)*50:0)>0?N3-TMath::Sqrt(N3)*50:0, N3?TMath::Sqrt(N3)*10+N3:10);
    TH1F *h  = new TH1F("h",  "N2*N3/N1", 100000, (N3*N2/N1 ? N3*N2/N1*(1-10/TMath::Sqrt(N3)) : 0)>0 ? N3*N2/N1*(1-10/TMath::Sqrt(N3)) : 0, 
                                                 N3*N2/N1 ? N3*N2/N1*(1+10/TMath::Sqrt(N3)) : 10 );
    TH1F *hband = (TH1F*)h->Clone("hband"); 

    for(int i=0; i<Npseudo; i++)  
    { 
        if((i%10000)==0) cout << "Generated " << i << " experiments" << endl;

        float N1temp   = N1;
        float N2temp   = N2;
        float N3temp   = N3;
        
        if(flucN1) N1temp   = fN1->Poisson(N1);
        if(flucN2) N2temp   = fN2->Poisson(N2);
        if(flucN3) N3temp   = fN3->Poisson(N3);

 //       cout << N1temp << " " << N2temp << " " << N3temp << endl;

        float kappaData = N3temp * N2temp / N1temp;
        if(N1temp==0) kappaData = 999999.;
        //float kappaData = N3temp;
        
        // Check if N is within 68% interval in Gamma 
        bool Within68 = IsWithin1Sig(i, N1, N2, N3, (int)N1temp, (int)N2temp, (int)N3temp); 
        //cout << Within68 << endl;
/*
        if(doGamma) 
        {
            h1->Fill(N1temp);
            h2->Fill(N2temp);
            h3->Fill(N3temp);
            if(flucN1) N1temp = gsl_ran_gamma(N1temp+1,1,rand); 
            if(flucN2) N2temp = gsl_ran_gamma(N2temp+1,1,rand); 
            if(flucN3) N3temp = gsl_ran_gamma(N3temp+1,1,rand); 
        }
        else 
        { 
            h1->Fill(N1temp);
            h2->Fill(N2temp);
            h3->Fill(N3temp);
        }
*/        
        h1->Fill(N1temp);
        h2->Fill(N2temp);
        h3->Fill(N3temp);

        h->Fill(kappaData);
        if(Within68) hband->Fill(kappaData); 
        //cout << N1temp << " " << N2temp << " " << N3temp << endl;
    }
    hband->Scale(1./h->Integral());
    h->Scale(1./h->Integral());
   
    cout << "h      : " << h->Integral() << endl;
    cout << "hband  : " << hband->Integral() << endl;

    TH1F *hkappa = (TH1F*)h->Clone("hkappa");
    hkappa->SetTitle("(N2 #times N3) / N1");
    TH1F *hkappaband = (TH1F*)hband->Clone("hkappa");

    TCanvas *c = new TCanvas("c","c",800,600);
    c->Divide(2,2);
    c->cd(1);
    h2->Draw();
    c->cd(2);
    h3->Draw();
    c->cd(3);
    h1->Draw();
    c->cd(4);
    hkappa->Draw();
   
    //
    // Get 1sigma band
    //
    float   kappa    = N3 * N2 / N1;  
    //int iN1max = h1->GetMaximumBin();
    //int iN2max = h2->GetMaximumBin();
    //int iN3max = h3->GetMaximumBin();
    //kappa    = h3->GetBinLowEdge(iN3max) * h2->GetBinLowEdge(iN2max) / h1->GetBinLowEdge(iN1max);  
    //kappa    = h->GetBinCenter(h->GetMaximumBin());  
    
    int     ibinmod     = h->FindBin(kappa); 
    //cout << "i bin mod " << ibinmod << endl; 
   
    float msigma=0;
    float psigma=0;
    int   imsigma=-1;
    int   ipsigma=-1;
    for(int ip=ibinmod; ip<=h->GetXaxis()->GetNbins(); ip++) 
    { 
        psigma = psigma + h->GetBinContent(ip); 
        //cout << psigma << " " << ipsigma << endl;
        if(psigma > /*0.34*/0.4987/*0.475*/) 
        { 
            ipsigma = ip;
            break;
        }
    }
    
    for(int im=ibinmod-1; im>0; im--) 
    { 
        if(im!=0) msigma = msigma + h->GetBinContent(im); 
        //cout << msigma << " " << imsigma << endl;
        if(msigma > /*0.34*/0.4987/*0.475*/) 
        { 
            imsigma = im;
            break;
        }
    }
    
    //cout << " ----------------------------- "  << endl;
    //cout << imsigma << " " << ipsigma << endl;
    //cout << msigma << " " << psigma << endl;
    //cout << " ----------------------------- "  << endl;

    if(imsigma==-1) 
    { 
        psigma = 0; 
        for(int i=1; i<=h->GetXaxis()->GetNbins(); i++)  
        { 
            psigma = psigma + h->GetBinContent(i); 
            if(psigma > /*0.34*2*/0.4987*2/*0.475*2*/) 
            { 
                ipsigma = i;
                break;
            }
        }
        msigma=0;
    }
   
    if(imsigma==-1) imsigma=1;
/*
    TH1F *hband = (TH1F*)hkappa->Clone("hband");
    hband->SetTitle("1sig band");
    hband->Reset(); 
    for(int i=imsigma; i<=ipsigma; i++) 
    { 
        hband->SetBinContent(i,hkappa->GetBinContent(i));
    }
*/

    hkappa->Rebin(250);
    hkappaband->Rebin(250);
    
    hkappaband->SetFillColor(kGreen); 

    hkappaband->Draw("HIST SAME"); 

    cout << " ----------------------------- "  << endl;
    cout << imsigma << " " <<  ipsigma << endl;
    cout << msigma << " " << psigma << endl;
    cout << " ----------------------------- "  << endl;

    //cout << kappa << " + " << (h->GetBinLowEdge(ipsigma)-kappa) 
    //              << " - " << (kappa-h->GetBinLowEdge(imsigma)) << endl;
    cout << kappa << " + " << h->GetBinLowEdge(ipsigma+1) 
                     << " - " << h->GetBinLowEdge(imsigma) << endl;
    cout << " ----------------------------- "  << endl;
    
    cout << " Linear approximation -------- "  << endl;
    float N2N3err = GetXErr(N2,N3,N2err,N3err);
    float kappaerr = GetRErr(N2*N3,N1,N2N3err,N1err); 
    cout << N3 * N2 / N1 << " +/- " << kappaerr << endl; 
   
    // print log 
    cout << " LOG : " << Npseudo << "\t" 
                      << N1 << "\t" << N2 << "\t" << N3 << "\t" 
                      << flucN1 << "\t" << flucN2 << "\t" << flucN3 << "\t" 
                      << scale1 << "\t" << scale2 << "\t" << scale3 << "\t" 
                      << Form("%.2f",kappa) << "\t" << Form("%.2f",h->GetBinLowEdge(ipsigma+1)-kappa) << "\t" << Form("%.2f",kappa-h->GetBinLowEdge(imsigma)) 
         << endl; 

    TLine *linekappa;
    linekappa = new TLine(kappa, 0, kappa, hkappa->GetMaximum()*1.05);
    linekappa->SetLineStyle(2);
    linekappa->SetLineWidth(3);
    linekappa->SetLineColor(kRed);
    linekappa->Draw("same");


    //TLatex *tex_kappa = new TLatex(0.7,0.3,Form("#frac{N2 #times N3}{N1}=%.2f^{%.2f}_{%.2f}",kappa,h->GetBinLowEdge(ipsigma)-kappa,kappa-h->GetBinLowEdge(imsigma)));
    TLatex *tex_kappa = new TLatex(0.8,0.3,Form("%.2f^{+%.2f}_{-%.2f}",kappa,h->GetBinLowEdge(ipsigma+1)-kappa,kappa-h->GetBinLowEdge(imsigma)));
    tex_kappa->SetNDC();
    tex_kappa->SetTextSize(0.05);
    tex_kappa->SetLineWidth(2);
    tex_kappa->SetTextAlign(22);
    tex_kappa->SetTextColor(kRed);
    tex_kappa->Draw("same");
    
    //TLatex *tex_kappa_lin = new TLatex(0.8,0.2,Form("%.2f #pm %.2f", N3*N2/N1,kappaerr));
    //tex_kappa_lin->SetNDC();
    //tex_kappa_lin->SetTextSize(0.05);
    //tex_kappa_lin->SetLineWidth(2);
    //tex_kappa_lin->SetTextAlign(22);
    //tex_kappa_lin->Draw("same");
    
    TLatex *tex_coverage = new TLatex(0.8,0.2,Form("Coverage = %.2f", hband->Integral()));
    tex_coverage->SetNDC();
    tex_coverage->SetTextSize(0.05);
    tex_coverage->SetLineWidth(2);
    tex_coverage->SetTextAlign(22);
    tex_coverage->Draw("same");

    TString strN1 = Form("%.1f",N1);
    TString strN2 = Form("%.1f",N2);
    TString strN3 = Form("%.1f",N3); 
    strN1.ReplaceAll(".","p");
    strN2.ReplaceAll(".","p");
    strN3.ReplaceAll(".","p");
    
    TString strscale1 = Form("%.1f",scale1);
    TString strscale2 = Form("%.1f",scale2);
    TString strscale3 = Form("%.1f",scale3); 
    strscale1.ReplaceAll(".","p");
    strscale2.ReplaceAll(".","p");
    strscale3.ReplaceAll(".","p");

    c->Print(Form("fig/N2N3overN1_Npseudo%i_N1%s_N2%s_N3%s_fluct%i%i%i_scale_N1%s_N2%s_N3%s.pdf", Npseudo, 
                                                                 strN1.Data(),strN2.Data(),strN3.Data(), 
                                                                 (int)flucN1, (int)flucN2, (int)flucN3, 
                                                                 strscale1.Data(), strscale2.Data(), strscale3.Data() 
            ));

    delete h1;
    delete h2;
    delete h3;
    delete h;
    delete hkappa;
    delete c;
    delete fN1;
    delete fN2;
    delete fN3;

}

void KappaPseudo(float N1=1, float N2=1, float N3=1)
{ 
    //KappaPseudoOneConfig( N1, N2, N3, true,     true,   true,   1., 1., 1.);

    //KappaPseudoOneConfig(506.19,67.05,43.24,true,true,true,0.3,0.3,0.3);  
    //KappaPseudoOneConfig(227.56,35.59,17.52,true,true,true,0.3,0.3,0.3);
    KappaPseudoOneConfig(151.47,26.76,10.80,true,true,true,0.3,0.3,0.3);
    KappaPseudoOneConfig(103.03,19.62,6.8,true,true,true,0.3,0.3,0.3);
    KappaPseudoOneConfig(71.33,15.13,5.64,true,true,true,0.3,0.3,0.3);
    KappaPseudoOneConfig(47.06,10.58,4.78,true,true,true,0.3,0.3,0.3);
    KappaPseudoOneConfig(32.7,7.8,3.7,true,true,true,0.3,0.3,0.3);

/*
    KappaPseudoOneConfig( N1, N2, N3, true,     true,   false,  1., 1., 1.); 
    KappaPseudoOneConfig( N1, N2, N3, true,     false,  true,   1., 1., 1.); 
    KappaPseudoOneConfig( N1, N2, N3,false,    true,   true,   1., 1., 1.);  

    for(int i1=0; i1<6; i1++) 
    { 
        for(int i2=0; i2<6; i2++) 
        { 
            for(int i3=0; i3<6; i3++) 
            { 
                KappaPseudoOneConfig( true,     true,   true,   1.+i1*0.2, 1.+i2*0.2, 1.+i3*0.2); 
            }
        }
    }
*/
}
