#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip> // for setw()

#include "TH1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TInterpreter.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMath.h"

using namespace std;


//
// h1 cosmetics
//
void h1cosmetic(TH1F* &h1, char* title, int linecolor=kBlack, int linewidth=1, int fillcolor=0, TString var="")
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle("nj");
    h1->SetYTitle(var);
    h1->SetStats(0);
    h1->SetMarkerColor(linecolor);
    h1->SetMarkerSize(1);
    h1->SetMarkerStyle(20);
    h1->SetStats(0);
    //h1->SetMinimum(0.1);
}

float GetRErr(float Num, float Den, float NumErr, float DenErr)
{
    float R = Num/Den;
    float RErr = R*TMath::Sqrt( DenErr*DenErr/Den/Den + NumErr*NumErr/Num/Num );
    if(Num==0) return 0;
    else return RErr;
}

void DrawKappaSummary(int nblow_)
{

    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");
    
    string Bin;
    int METlow  = -999;
    int METhigh = -999;
    int njhigh  = -999;
    int njlow   = -999;
    int nbhigh  = -999;
    int nblow   = -999;
    int MJ      = -999;
    string R;
    string equal;
    float RMC  = -999;
    string MCpm;
    float RMCerr  = -999;
    float Rdata  = -999;
    string datapm;
    float Rdataerr  = -999;
    
    float kappa  = -999;
    string kappapm;
    float kappaerr  = -999;


    // TH1F 
    TH1F *h1_MC     = new TH1F("h1_MC",     "h1_MC",        32,     0,      32);
    TH1F *h1_data   = new TH1F("h1_data",   "h1_data",      32,     0,      32);
    TH1F *h1_ratio  = new TH1F("h1_ratio",  "h1_ratio",     8,      0,      8);
   
    const int numofjbin =  8;
    TLatex *tex_nj[numofjbin]; 
    for(int i=0; i<numofjbin; i++) 
    {
        char *njbin = "n_{j}=5";
        if(nblow_==1) 
        {
           if((i%3)==0) njbin="n_{j}=6"; 
           if((i%3)==1) njbin="n_{j}=7"; 
           if((i%3)==2) njbin="n_{j}#geq8"; 
        }
        if(nblow_==2) 
        {
           if((i%4)==0) njbin="n_{j}=5"; 
           if((i%4)==1) njbin="n_{j}=6"; 
           if((i%4)==2) njbin="n_{j}=7"; 
           if((i%4)==3) njbin="n_{j}#geq8"; 
        }
        float width = 0.7 / (2*(2+nblow_)); 
        tex_nj[i] = new TLatex(0.2+(i+0.5)*width,0.2,njbin);
        tex_nj[i]->SetNDC();
        tex_nj[i]->SetTextSize(0.035);
        tex_nj[i]->SetLineWidth(2);
        tex_nj[i]->SetTextAlign(22);
    }

    TLatex *tex_MET[2]; 
    for(int i=0; i<2; i++) 
    {
        char *METbin = "100<MET<200 GeV";
        if((i%2)==1) METbin = "200<MET<300 GeV";

        float width = 0.7 / 2; 
        tex_MET[i] = new TLatex(0.205+(i+0.5)*width,0.13,METbin);
        tex_MET[i]->SetNDC();
        tex_MET[i]->SetTextSize(0.04);
        tex_MET[i]->SetLineWidth(2);
        tex_MET[i]->SetTextAlign(22);
    }
    
    
    TLatex *tex_selection; 
    tex_selection = new TLatex(0.55,0.95,"n_{lep}=1, H_{T}>500 GeV, MET>100 GeV, n_{j}#geq5, n_{b}#geq1");
    tex_selection->SetNDC();
    tex_selection->SetTextSize(0.045);
    tex_selection->SetLineWidth(2);
    tex_selection->SetTextAlign(22);
    tex_selection->SetTextFont(42);


    string line;
    //ifstream infile ("log_all_v3.txt");
    //ifstream infile ("log_all_01Jun2015_For29May2015slides.txt");
    //ifstream infile ("log_all_01Jun2015_For29May2015slides_NoLepVeto.txt");
    //ifstream infile ("log_all_10Jun2015.txt");
    ifstream infile ("log_all_11Jun2015.txt");
    if (infile.is_open())
    {
        while ( infile.good() )
        {

            // get a line from input file
            getline (infile,line);

            if( line.find("[Bin]")==string::npos )  continue;

            if( line.find("kappa")!=string::npos ) 
            {
                // streaming entry of string to doubles
                stringstream stream(line);
                stream >> Bin;
                stream >> METlow;
                stream >> METhigh;
                stream >> njlow;
                stream >> njhigh;
                stream >> nblow;
                stream >> nbhigh;
                stream >> MJ;
                stream >> R;
                stream >> equal;
                stream >> RMC;
                stream >> MCpm;
                stream >> RMCerr;
                stream >> Rdata;
                stream >> datapm;
                stream >> Rdataerr;
                
                if(MJ!=500) continue; 
                
                if(nblow_==2/*&& nblow==2 && nbhigh>8*/) 
                {
                    if(METlow==100 && METhigh==200) 
                    {
                        if(njlow==5 && njhigh==5)
                        {
                            h1_MC->SetBinContent(1*4-2,RMC);
                            h1_MC->SetBinError(1*4-2,RMCerr);
                            h1_data->SetBinContent(1*4-1,Rdata);
                            h1_data->SetBinError(1*4-1,Rdata==0?0:Rdataerr);
                        }
                        if(njlow==6 && njhigh==6)
                        {
                            h1_MC->SetBinContent(2*4-2,RMC);
                            h1_MC->SetBinError(2*4-2,RMCerr);
                            h1_data->SetBinContent(2*4-1,Rdata);
                            h1_data->SetBinError(2*4-1,Rdata==0?0:Rdataerr);
                        }
                        if(njlow==7 && njhigh==7)
                        {
                            h1_MC->SetBinContent(3*4-2,RMC);
                            h1_MC->SetBinError(3*4-2,RMCerr);
                            h1_data->SetBinContent(3*4-1,Rdata);
                            h1_data->SetBinError(3*4-1,Rdata==0?0:Rdataerr);
                        }
                        if(njlow==8)
                        {
                            h1_MC->SetBinContent(4*4-2,RMC);
                            h1_MC->SetBinError(4*4-2,RMCerr);
                            h1_data->SetBinContent(4*4-1,Rdata);
                            h1_data->SetBinError(4*4-1,Rdata==0?0:Rdataerr);
                        }
                    }
                    if(METlow==200 && METhigh==300) 
                    {
                        if(njlow==5 && njhigh==5)
                        {
                            h1_MC->SetBinContent(5*4-2,RMC);
                            h1_MC->SetBinError(5*4-2,RMCerr);
                            h1_data->SetBinContent(5*4-1,Rdata);
                            h1_data->SetBinError(5*4-1,Rdata==0?0:Rdataerr);
                        }
                        if(njlow==6 && njhigh==6)
                        {
                            h1_MC->SetBinContent(6*4-2,RMC);
                            h1_MC->SetBinError(6*4-2,RMCerr);
                            h1_data->SetBinContent(6*4-1,Rdata);
                            h1_data->SetBinError(6*4-1,Rdata==0?0:Rdataerr);
                        }
                        if(njlow==7 && njhigh==7)
                        {
                            h1_MC->SetBinContent(7*4-2,RMC);
                            h1_MC->SetBinError(7*4-2,RMCerr);
                            h1_data->SetBinContent(7*4-1,Rdata);
                            h1_data->SetBinError(7*4-1,Rdata==0?0:Rdataerr);
                        }
                        if(njlow==8)
                        {
                            h1_MC->SetBinContent(8*4-2,RMC);
                            h1_MC->SetBinError(8*4-2,RMCerr);
                            h1_data->SetBinContent(8*4-1,Rdata);
                            h1_data->SetBinError(8*4-1,Rdata==0?0:Rdataerr);
                        }
                    }
                } //  if(nblow_==2 && nblow==2 && nbhigh==99)

            }   

            if( !infile.good() ) continue;
        }
    }
    infile.close();

    char* selection = "n_{b} = 1";
    if(nblow_==2) selection = "n_{b}#geq1 (#geq2 for n_{j}=5)";
    TString YTitle = Form("#kappa [ %s ]", selection);
    char* Title = Form("nb %i",nblow_); 
    h1cosmetic(h1_MC,       Title,  kRed,     1, 0,   YTitle);
    h1cosmetic(h1_data,     Title,  kBlack,   1, 0,   YTitle);
    h1_MC->SetMarkerStyle(21); 
    h1_MC->SetMarkerSize(0.7); 

    TCanvas *c = new TCanvas("c","c",900,600);
    c->cd();
    TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.25);
    pad1->SetRightMargin(0.1);
    pad1->SetLeftMargin(0.2);
    pad1->SetLogy(0);
    pad1->Draw();
    pad1->cd();

    h1_MC->GetXaxis()->SetLabelSize(0.);
    h1_MC->GetXaxis()->SetNdivisions(8, kFALSE);
    h1_MC->GetXaxis()->SetTitleSize(0.);
    h1_MC->SetMaximum(5.);
    h1_MC->SetMinimum(0.);
    h1_MC->Draw("e");
    h1_data->Draw("same e");
   
    for(int i=0; i<8; i++)  tex_nj[i]->Draw("same");
    for(int i=0; i<2; i++)  tex_MET[i]->Draw("same");
    tex_selection->Draw("same");

    TLine *lineMET[3];  
    for(int i=0; i<3; i++)
    {   
        float width = 2*numofjbin;
        lineMET[i] = new TLine(i*width, -1.2, i*width, 0);  
        lineMET[i]->SetLineStyle(3);
        lineMET[i]->SetLineWidth(2);
        lineMET[i]->Draw("same");
    }

    TLegend *leg = new TLegend(0.5, 0.7, 0.8, 0.85, "");
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);

    leg->AddEntry(h1_data,   " Data (8 TeV, 19.5 fb^{-1})","lp");
    leg->AddEntry(h1_MC,     " MC","lp");
    leg->Draw();


    c->cd();
    TPad *pad2 = new TPad("p_pull", "p_pull", 0.0, 0.0, 1.0, 0.3);
    pad2->SetLeftMargin(0.2);
    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0.1);
    pad2->SetRightMargin(0.1);
    pad2->SetBottomMargin(0.4); 




    // Fill h1_ratio 
    // MC = 2, 6, 10, 14, ...  
    // data = 3, 7, 11, 15, ...   
    for(int i=1; i<=8; i++) 
    { 
        float data      = h1_data->GetBinContent(i*4-1); 
        float dataerr   = h1_data->GetBinError(i*4-1); 
        float mc        = h1_MC->GetBinContent(i*4-2); 
        float mcerr     = h1_MC->GetBinError(i*4-2);
        float ratio     = data/mc; 
        float ratioerr  = GetRErr(data,mc,dataerr,mcerr); 
        h1_ratio->SetBinContent(i, ratio);
        h1_ratio->SetBinError(i, ratioerr);
        cout << i <<  " " 
             << data << " +/- " << dataerr  << " ,,, " 
             << mc << " +/- " << mcerr 
             << " => " << ratio << " +/- " << ratioerr << endl;
    }
    
    TH1F *h1_One       = new TH1F("h1_One", "h1_One", 1, 0, 8);
    h1_One->SetBinContent(1,1);
    h1_One->GetXaxis()->SetNdivisions(16, kFALSE);
    h1_One->GetXaxis()->SetLabelSize(0.);
    h1_One->GetXaxis()->SetTitleSize(0.);
    h1_One->GetYaxis()->SetNdivisions(2, kFALSE);
    h1_One->GetYaxis()->SetLabelSize(0.15);
    h1_One->GetYaxis()->SetTitleSize(0.15);
    h1_One->GetYaxis()->SetTitleOffset(0.6);
    h1_One->SetMinimum(0);
    h1_One->SetMaximum(4);
    
    h1cosmetic(h1_ratio,       Title,  kBlack,     1, 0,  "Data/MC");
    h1cosmetic(h1_One,         " ",    kBlue,      1, 0,  "Data/MC"); 
    
    h1_One->Draw("hist");
    h1_ratio->Draw("same p");


    c->Print("figR/KappaSummary_nb1or2_MJ500.pdf");
    delete c;
}

void DrawRSummary(int nblow_, TString R_)
{

    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");
    
    string Bin;
    int METlow  = -999;
    int METhigh = -999;
    int njhigh  = -999;
    int njlow   = -999;
    int nbhigh  = -999;
    int nblow   = -999;
    int MJ      = -999;
    string R;
    string equal;
    float RMC  = -999;
    string MCpm;
    float RMCerr  = -999;
    float Rdata  = -999;
    string datapm;
    float Rdataerr  = -999;
    
    float kappa  = -999;
    string kappapm;
    float kappaerr  = -999;


    // TH1F 
    TH1F *h1_MC     = new TH1F("h1_MC",     "h1_MC",    4*4*(2+nblow_), 0,  4*4*(2+nblow_));
    TH1F *h1_data   = new TH1F("h1_data",   "h1_data",  4*4*(2+nblow_), 0,  4*4*(2+nblow_));
    TH1F *h1_ratio  = new TH1F("h1_ratio",   "h1_ratio",  4*(2+nblow_),   0,  4*(2+nblow_));
   
    const int numofjbin =  4*(2+nblow_);
    TLatex *tex_nj[numofjbin]; 
    for(int i=0; i<numofjbin; i++) 
    {
        char *njbin = "n_{j}=5";
        if(nblow_==1) 
        {
           if((i%3)==0) njbin="n_{j}=6"; 
           if((i%3)==1) njbin="n_{j}=7"; 
           if((i%3)==2) njbin="n_{j}#geq8"; 
        }
        if(nblow_==2) 
        {
           if((i%4)==0) njbin="n_{j}=5"; 
           if((i%4)==1) njbin="n_{j}=6"; 
           if((i%4)==2) njbin="n_{j}=7"; 
           if((i%4)==3) njbin="n_{j}#geq8"; 
        }
        float width = 0.7 / (4*(2+nblow_)); 
        tex_nj[i] = new TLatex(0.2+(i+0.5)*width,0.2,njbin);
        tex_nj[i]->SetNDC();
        tex_nj[i]->SetTextSize(0.035);
        tex_nj[i]->SetLineWidth(2);
        tex_nj[i]->SetTextAlign(22);
    }

    TLatex *tex_MET[4]; 
    for(int i=0; i<4; i++) 
    {
        char *METbin = "100<MET<200 GeV";
        if((i%2)==1) METbin = "200<MET<300 GeV";

        float width = 0.7 / 4; 
        tex_MET[i] = new TLatex(0.205+(i+0.5)*width,0.13,METbin);
        tex_MET[i]->SetNDC();
        tex_MET[i]->SetTextSize(0.04);
        tex_MET[i]->SetLineWidth(2);
        tex_MET[i]->SetTextAlign(22);
    }
    
    TLatex *tex_MJ[2]; 
    for(int i=0; i<2; i++)  
    { 
        char *MJbin;
        if(R_=="RmT")
        {
            MJbin= "M_{J}<500 GeV";
            if((i%2)==1) MJbin = "M_{J}>500 GeV";
        } 
        if(R_=="RMJ")
        {
            MJbin= "m_{T}<140 GeV";
            if((i%2)==1) MJbin = "m_{T}<140 GeV";
        } 
        
        float width = 0.7 / 2; 
        tex_MJ[i] = new TLatex(0.2+(i+0.5)*width,0.05,MJbin);
        tex_MJ[i]->SetNDC();
        tex_MJ[i]->SetTextSize(0.05);
        tex_MJ[i]->SetLineWidth(2);
        tex_MJ[i]->SetTextAlign(22);

    }
    
    TLatex *tex_selection; 
    tex_selection = new TLatex(0.55,0.95,"n_{lep}=1, H_{T}>500 GeV, MET>100 GeV, n_{j}#geq5, n_{b}#geq1");
    tex_selection->SetNDC();
    tex_selection->SetTextSize(0.045);
    tex_selection->SetLineWidth(2);
    tex_selection->SetTextAlign(22);
    tex_selection->SetTextFont(42);


    string line;
    //ifstream infile ("log_all_v3.txt");
    //ifstream infile ("log_all_01Jun2015_For29May2015slides.txt");
    //ifstream infile ("log_all_01Jun2015_For29May2015slides_NoLepVeto.txt");
    //ifstream infile ("log_all_10Jun2015.txt");
    ifstream infile ("log_all_11Jun2015.txt");
    if (infile.is_open())
    {
        while ( infile.good() )
        {

            // get a line from input file
            getline (infile,line);

            if( line.find("[Bin]")==string::npos )  continue;

            if( line.find("kappa")==string::npos ) 
            {
                // streaming entry of string to doubles
                stringstream stream(line);
                stream >> Bin;
                stream >> METlow;
                stream >> METhigh;
                stream >> njlow;
                stream >> njhigh;
                stream >> nblow;
                stream >> nbhigh;
                stream >> MJ;
                stream >> R;
                stream >> equal;
                stream >> RMC;
                stream >> MCpm;
                stream >> RMCerr;
                stream >> Rdata;
                stream >> datapm;
                stream >> Rdataerr;
                
                if(MJ!=500) continue; 
                
                // low MJ 
                if( (R_=="RMJ" && R=="R21") || (R_=="RmT" && R=="R31") ) 
                { 
                    if(nblow_==1 && nblow==1 && nbhigh==1) 
                    {
                        if(METlow==100 && METhigh==200) 
                        {
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(1*4-2,RMC);
                                h1_MC->SetBinError(1*4-2,RMCerr);
                                h1_data->SetBinContent(1*4-1,Rdata);
                                h1_data->SetBinError(1*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(2*4-2,RMC);
                                h1_MC->SetBinError(2*4-2,RMCerr);
                                h1_data->SetBinContent(2*4-1,Rdata);
                                h1_data->SetBinError(2*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(3*4-2,RMC);
                                h1_MC->SetBinError(3*4-2,RMCerr);
                                h1_data->SetBinContent(3*4-1,Rdata);
                                h1_data->SetBinError(3*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                        if(METlow==200 && METhigh==300) 
                        {
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(4*4-2,RMC);
                                h1_MC->SetBinError(4*4-2,RMCerr);
                                h1_data->SetBinContent(4*4-1,Rdata);
                                h1_data->SetBinError(4*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(5*4-2,RMC);
                                h1_MC->SetBinError(5*4-2,RMCerr);
                                h1_data->SetBinContent(5*4-1,Rdata);
                                h1_data->SetBinError(5*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(6*4-2,RMC);
                                h1_MC->SetBinError(6*4-2,RMCerr);
                                h1_data->SetBinContent(6*4-1,Rdata);
                                h1_data->SetBinError(6*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                    } //  if(nblow_==1 && nblow==1 && nbhigh==1)
                    
                    if(nblow_==2/*&& nblow==2 && nbhigh>8*/) 
                    {
                        if(METlow==100 && METhigh==200) 
                        {
                            if(njlow==5 && njhigh==5)
                            {
                                h1_MC->SetBinContent(1*4-2,RMC);
                                h1_MC->SetBinError(1*4-2,RMCerr);
                                h1_data->SetBinContent(1*4-1,Rdata);
                                h1_data->SetBinError(1*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(2*4-2,RMC);
                                h1_MC->SetBinError(2*4-2,RMCerr);
                                h1_data->SetBinContent(2*4-1,Rdata);
                                h1_data->SetBinError(2*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(3*4-2,RMC);
                                h1_MC->SetBinError(3*4-2,RMCerr);
                                h1_data->SetBinContent(3*4-1,Rdata);
                                h1_data->SetBinError(3*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(4*4-2,RMC);
                                h1_MC->SetBinError(4*4-2,RMCerr);
                                h1_data->SetBinContent(4*4-1,Rdata);
                                h1_data->SetBinError(4*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                        if(METlow==200 && METhigh==300) 
                        {
                            if(njlow==5 && njhigh==5)
                            {
                                h1_MC->SetBinContent(5*4-2,RMC);
                                h1_MC->SetBinError(5*4-2,RMCerr);
                                h1_data->SetBinContent(5*4-1,Rdata);
                                h1_data->SetBinError(5*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(6*4-2,RMC);
                                h1_MC->SetBinError(6*4-2,RMCerr);
                                h1_data->SetBinContent(6*4-1,Rdata);
                                h1_data->SetBinError(6*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(7*4-2,RMC);
                                h1_MC->SetBinError(7*4-2,RMCerr);
                                h1_data->SetBinContent(7*4-1,Rdata);
                                h1_data->SetBinError(7*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(8*4-2,RMC);
                                h1_MC->SetBinError(8*4-2,RMCerr);
                                h1_data->SetBinContent(8*4-1,Rdata);
                                h1_data->SetBinError(8*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                    } //  if(nblow_==2 && nblow==2 && nbhigh==99)
                
                } // R21 
                
                // high MJ 
                if( (R_=="RMJ" && R=="R43") || (R_=="RmT" && R=="R42") ) 
                { 
                    if(nblow_==1 && nblow==1 && nbhigh==1) 
                    {
                        if(METlow==100 && METhigh==200) 
                        {
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(7*4-2,RMC);
                                h1_MC->SetBinError(7*4-2,RMCerr);
                                h1_data->SetBinContent(7*4-1,Rdata);
                                h1_data->SetBinError(7*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(8*4-2,RMC);
                                h1_MC->SetBinError(8*4-2,RMCerr);
                                h1_data->SetBinContent(8*4-1,Rdata);
                                h1_data->SetBinError(8*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(9*4-2,RMC);
                                h1_MC->SetBinError(9*4-2,RMCerr);
                                h1_data->SetBinContent(9*4-1,Rdata);
                                h1_data->SetBinError(9*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                        if(METlow==200 && METhigh==300) 
                        {
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(10*4-2,RMC);
                                h1_MC->SetBinError(10*4-2,RMCerr);
                                h1_data->SetBinContent(10*4-1,Rdata);
                                h1_data->SetBinError(10*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(11*4-2,RMC);
                                h1_MC->SetBinError(11*4-2,RMCerr);
                                h1_data->SetBinContent(11*4-1,Rdata);
                                h1_data->SetBinError(11*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(12*4-2,RMC);
                                h1_MC->SetBinError(12*4-2,RMCerr);
                                h1_data->SetBinContent(12*4-1,Rdata);
                                h1_data->SetBinError(12*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                    } //  if(nblow_==1 && nblow==1 && nbhigh==1)
                    
                    if(nblow_==2/* && nblow==2 && nbhigh>8*/) 
                    {
                        if(METlow==100 && METhigh==200) 
                        {
                            if(njlow==5 && njhigh==5)
                            {
                                h1_MC->SetBinContent(9*4-2,RMC);
                                h1_MC->SetBinError(9*4-2,RMCerr);
                                h1_data->SetBinContent(9*4-1,Rdata);
                                h1_data->SetBinError(9*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(10*4-2,RMC);
                                h1_MC->SetBinError(10*4-2,RMCerr);
                                h1_data->SetBinContent(10*4-1,Rdata);
                                h1_data->SetBinError(10*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(11*4-2,RMC);
                                h1_MC->SetBinError(11*4-2,RMCerr);
                                h1_data->SetBinContent(11*4-1,Rdata);
                                h1_data->SetBinError(11*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(12*4-2,RMC);
                                h1_MC->SetBinError(12*4-2,RMCerr);
                                h1_data->SetBinContent(12*4-1,Rdata);
                                h1_data->SetBinError(12*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                        if(METlow==200 && METhigh==300) 
                        {
                            if(njlow==5 && njhigh==5)
                            {
                                h1_MC->SetBinContent(13*4-2,RMC);
                                h1_MC->SetBinError(13*4-2,RMCerr);
                                h1_data->SetBinContent(13*4-1,Rdata);
                                h1_data->SetBinError(13*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==6 && njhigh==6)
                            {
                                h1_MC->SetBinContent(14*4-2,RMC);
                                h1_MC->SetBinError(14*4-2,RMCerr);
                                h1_data->SetBinContent(14*4-1,Rdata);
                                h1_data->SetBinError(14*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==7 && njhigh==7)
                            {
                                h1_MC->SetBinContent(15*4-2,RMC);
                                h1_MC->SetBinError(15*4-2,RMCerr);
                                h1_data->SetBinContent(15*4-1,Rdata);
                                h1_data->SetBinError(15*4-1,Rdata==0?0:Rdataerr);
                            }
                            if(njlow==8)
                            {
                                h1_MC->SetBinContent(16*4-2,RMC);
                                h1_MC->SetBinError(16*4-2,RMCerr);
                                h1_data->SetBinContent(16*4-1,Rdata);
                                h1_data->SetBinError(16*4-1,Rdata==0?0:Rdataerr);
                            }
                        }
                    } //  if(nblow_==2 && nblow==2 && nbhigh==99)
                } // R43
            }   

            if( !infile.good() ) continue;
        }
    }
    infile.close();

    char* selection = "n_{b} = 1";
    if(nblow_==2) selection = "n_{b}#geq1 (#geq2 for n_{j}=5)";
    TString YTitle = Form("R_{m_{T}} [ %s ]", selection);
    if(R_=="RMJ") YTitle =YTitle =  Form("R_{M_{J}} [ %s ]", selection);
    char* Title = Form("nb %i",nblow_); 
    h1cosmetic(h1_MC,       Title,  kRed,     1, 0,   YTitle);
    h1cosmetic(h1_data,     Title,  kBlack,   1, 0,   YTitle);
    h1_MC->SetMarkerStyle(21); 
    h1_MC->SetMarkerSize(0.7); 

    TCanvas *c = new TCanvas("c","c",900,600);
    c->cd();
    TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.25);
    pad1->SetRightMargin(0.1);
    pad1->SetLeftMargin(0.2);
    pad1->SetLogy(0);
    pad1->Draw();
    pad1->cd();

    float HistMax=0.3;
    if(R_=="RMJ") HistMax=0.5; 

    h1_MC->GetXaxis()->SetLabelSize(0.);
    h1_MC->GetXaxis()->SetNdivisions(16, kFALSE);
    h1_MC->GetXaxis()->SetTitleSize(0.);
    h1_MC->SetMaximum(HistMax);
    h1_MC->SetMinimum(0.);
    h1_MC->Draw("e");
    h1_data->Draw("same e");
    
    for(int i=0; i<4*(2+nblow_); i++)  tex_nj[i]->Draw("same");
    for(int i=0; i<4; i++)  tex_MET[i]->Draw("same");
    for(int i=0; i<2; i++)  tex_MJ[i]->Draw("same");
    tex_selection->Draw("same");

    TLine *lineMET[5];  
    for(int i=0; i<5; i++)
    {   
        float width = numofjbin;
        if(i==0 || i==4) lineMET[i] = new TLine(i*width, -HistMax*0.33, i*width, 0);  
        else if(i==2) lineMET[i] = new TLine(i*width, -HistMax*0.33, i*width, 0.3);  
        else lineMET[i] = new TLine(i*width, -HistMax*0.2, i*width, 0);  
        lineMET[i]->SetLineStyle(3);
        lineMET[i]->SetLineWidth(2);
        lineMET[i]->Draw("same");
    }

    TLegend *leg = new TLegend(0.22, 0.7, 0.5, 0.85, "");
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);

    leg->AddEntry(h1_data,   " Data (8 TeV, 19.5 fb^{-1})","lp");
    leg->AddEntry(h1_MC,     " MC","lp");
    leg->Draw();

    c->cd();
    TPad *pad2 = new TPad("p_pull", "p_pull", 0.0, 0.0, 1.0, 0.3);
    pad2->SetLeftMargin(0.2);
    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0.1);
    pad2->SetRightMargin(0.1);
    pad2->SetBottomMargin(0.4); 




    // Fill h1_ratio 
    // MC = 2, 6, 10, 14, ...  
    // data = 3, 7, 11, 15, ...   
    for(int i=1; i<=16; i++) 
    { 
        float data      = h1_data->GetBinContent(i*4-1); 
        float dataerr   = h1_data->GetBinError(i*4-1); 
        float mc        = h1_MC->GetBinContent(i*4-2); 
        float mcerr     = h1_MC->GetBinError(i*4-2);
        float ratio     = data/mc; 
        float ratioerr  = GetRErr(data,mc,dataerr,mcerr); 
        h1_ratio->SetBinContent(i, ratio);
        h1_ratio->SetBinError(i, ratioerr);
        cout << i <<  " " 
             << data << " +/- " << dataerr  << " ,,, " 
             << mc << " +/- " << mcerr 
             << " => " << ratio << " +/- " << ratioerr << endl;
    }
    
    TH1F *h1_One       = new TH1F("h1_One", "h1_One", 1, 0, 16);
    h1_One->SetBinContent(1,1);
    h1_One->GetXaxis()->SetNdivisions(16, kFALSE);
    h1_One->GetXaxis()->SetLabelSize(0.);
    h1_One->GetXaxis()->SetTitleSize(0.);
    h1_One->GetYaxis()->SetNdivisions(3, kFALSE);
    h1_One->GetYaxis()->SetLabelSize(0.15);
    h1_One->GetYaxis()->SetTitleSize(0.15);
    h1_One->GetYaxis()->SetTitleOffset(0.6);
    h1_One->SetMinimum(0);
    h1_One->SetMaximum(3);
    
    h1cosmetic(h1_ratio,       Title,  kBlack,     1, 0,  "Data/MC");
    h1cosmetic(h1_One,         " ",    kBlue,      1, 0,  "Data/MC"); 
    
    h1_One->Draw("hist");
    h1_ratio->Draw("same p");


    c->Print(Form("figR/%sSummary_nb1or2_MJ500.pdf",R_.Data()));
    delete c;
}

void DrawRone(int METlow_=100, int METhigh_=200, int njlow_=6, int njhigh_=99, int nblow_=2, int nbhigh_=9, int MJ_=500)
{

    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");
    
    string Bin;
    int METlow  = -999;
    int METhigh = -999;
    int njhigh  = -999;
    int njlow   = -999;
    int nbhigh  = -999;
    int nblow   = -999;
    int MJ      = -999;
    string R;
    string equal;
    float RMC  = -999;
    string MCpm;
    float RMCerr  = -999;
    float Rdata  = -999;
    string datapm;
    float Rdataerr  = -999;
    
    float kappa  = -999;
    string kappapm;
    float kappaerr  = -999;



    // TH1F 
    TH1F *h1_R21_MC     = new TH1F("h1_R21_MC",     "h1_R21_MC",    4, 4.5,  8.5);
    TH1F *h1_R21_data   = new TH1F("h1_R21_data",   "h1_R21_data",  4, 4.5,  8.5);
    TH1F *h1_R43_MC     = new TH1F("h1_R43_MC",     "h1_R43_MC",    4, 4.5,  8.5);
    TH1F *h1_R43_data   = new TH1F("h1_R43_data",   "h1_R43_data",  4, 4.5,  8.5);
    TH1F *h1_R31_MC     = new TH1F("h1_R31_MC",     "h1_R31_MC",    4, 4.5,  8.5);
    TH1F *h1_R31_data   = new TH1F("h1_R31_data",   "h1_R31_data",  4, 4.5,  8.5);
    TH1F *h1_R42_MC     = new TH1F("h1_R42_MC",     "h1_R42_MC",    4, 4.5,  8.5);
    TH1F *h1_R42_data   = new TH1F("h1_R42_data",   "h1_R42_data",  4, 4.5,  8.5);
    
    TH1F *h1_kappa_MC   = new TH1F("h1_kappa_MC",   "h1_kappa_MC",    4, 4.5,  8.5);
    TH1F *h1_kappa_Data = new TH1F("h1_kappa_Data", "h1_kappa_Data",  4, 4.5,  8.5);

    string line;
    //ifstream infile ("log_all_v3.txt");
    //ifstream infile ("log_all_01Jun2015_For29May2015slides.txt");
    ifstream infile ("log_all_01Jun2015_For29May2015slides_NoLepVeto.txt");
    if (infile.is_open())
    {
        while ( infile.good() )
        {

            // get a line from input file
            getline (infile,line);

            if( line.find("[Bin]")==string::npos )  continue;

            if( line.find("kappa")==string::npos ) 
            {
                // streaming entry of string to doubles
                stringstream stream(line);
                stream >> Bin;
                stream >> METlow;
                stream >> METhigh;
                stream >> njlow;
                stream >> njhigh;
                stream >> nblow;
                stream >> nbhigh;
                stream >> MJ;
                stream >> R;
                stream >> equal;
                stream >> RMC;
                stream >> MCpm;
                stream >> RMCerr;
                stream >> Rdata;
                stream >> datapm;
                stream >> Rdataerr;

                if(METlow!=METlow_) continue;
                if(nblow!=nblow_) continue;
                if(nbhigh!=nbhigh_) continue;
                if(MJ!=MJ_) continue; 
                if(njlow<njlow_) continue; 
                
                cout << setw(10) << R;
                cout << setw(10) << equal;
                cout << setw(10) << RMC;
                cout << setw(10) << MCpm;
                cout << setw(10) << RMCerr;
                cout << setw(10) << Rdata;
                cout << setw(10) << datapm;
                cout << setw(10) << Rdataerr;
                cout << endl;
                
                if(R=="R21")
                { 
                    h1_R21_MC->SetBinContent(njlow-4,RMC);
                    h1_R21_MC->SetBinError(njlow-4,RMCerr);
                    h1_R21_data->SetBinContent(njlow-4,Rdata);
                    h1_R21_data->SetBinError(njlow-4,Rdata==0?0:Rdataerr);
                }
                if(R=="R43")
                { 
                    h1_R43_MC->SetBinContent(njlow-4,RMC);
                    h1_R43_MC->SetBinError(njlow-4,RMCerr);
                    h1_R43_data->SetBinContent(njlow-4,Rdata);
                    h1_R43_data->SetBinError(njlow-4,Rdata==0?0:Rdataerr);
                }
                if(R=="R31")
                { 
                    h1_R31_MC->SetBinContent(njlow-4,RMC);
                    h1_R31_MC->SetBinError(njlow-4,RMCerr);
                    h1_R31_data->SetBinContent(njlow-4,Rdata);
                    h1_R31_data->SetBinError(njlow-4,Rdata==0?0:Rdataerr);
                }
                if(R=="R42")
                { 
                    h1_R42_MC->SetBinContent(njlow-4,RMC);
                    h1_R42_MC->SetBinError(njlow-4,RMCerr);
                    h1_R42_data->SetBinContent(njlow-4,Rdata);
                    h1_R42_data->SetBinError(njlow-4,Rdata==0?0:Rdataerr);
                }
/*
                cout << setw(10) << Bin;
                cout << setw(10) << METlow;
                cout << setw(10) << METhigh;
                cout << setw(10) << njlow;
                cout << setw(10) << njhigh;
                cout << setw(10) << nblow;
                cout << setw(10) << nbhigh;
                cout << setw(10) << MJ;
                cout << setw(10) << R;
                cout << setw(10) << equal;
                cout << setw(10) << RMC;
                cout << setw(10) << MCpm;
                cout << setw(10) << RMCerr;
                cout << setw(10) << Rdata;
                cout << setw(10) << datapm;
                cout << setw(10) << Rdataerr;
                cout << endl;
*/
            }
            else 
            { 
                stringstream stream(line);
                stream >> Bin;
                stream >> METlow;
                stream >> METhigh;
                stream >> njlow;
                stream >> njhigh;
                stream >> nblow;
                stream >> nbhigh;
                stream >> MJ;
                stream >> R;
                stream >> equal;
                stream >> RMC;
                stream >> MCpm;
                stream >> RMCerr;
                stream >> Rdata;
                stream >> datapm;
                stream >> Rdataerr;
                
                if(METlow!=METlow_) continue;
                if(nblow!=nblow_) continue;
                if(nbhigh!=nbhigh_) continue;
                if(MJ!=MJ_) continue; 
                if(njlow<njlow_) continue; 
 /*               
                cout << setw(10) << Bin;
                cout << setw(10) << METlow;
                cout << setw(10) << METhigh;
                cout << setw(10) << njlow;
                cout << setw(10) << njhigh;
                cout << setw(10) << nblow;
                cout << setw(10) << nbhigh;
                cout << setw(10) << MJ;
                cout << setw(10) << R;
                cout << setw(10) << equal;
                cout << setw(10) << kappa;
                cout << setw(10) << kappapm;
                cout << setw(10) << kappaerr;
                cout << endl;
*/
                h1_kappa_MC->SetBinContent(njlow-4,RMC);
                h1_kappa_MC->SetBinError(njlow-4,RMCerr);
                h1_kappa_Data->SetBinContent(njlow-4,Rdata);
                h1_kappa_Data->SetBinError(njlow-4,Rdata==0?0:Rdataerr);
            }
            
            if( !infile.good() ) continue;
        }
    }
    infile.close();

    char* Title = Form("%i<MET<%i %i<=nj<%i %i<=nb<%i MJ<%i",METlow_,METhigh_,njlow_,njhigh_,nblow_,nbhigh_,MJ_); 
    h1cosmetic(h1_R21_MC,       Title,  kRed,   1, 0,   "R_{M_{J}} (m_{T}<140 GeV)");
    h1cosmetic(h1_R21_data,     Title,  kBlack, 1, 0,   "R_{M_{J}} (m_{T}<140 GeV)");
    h1cosmetic(h1_R43_MC,       Title,  kRed,   1, 0,   "R_{M_{J}} (m_{T}>140 GeV)");
    h1cosmetic(h1_R43_data,     Title,  kBlack, 1, 0,   "R_{M_{J}} (m_{T}>140 GeV)");
    h1cosmetic(h1_R31_MC,       Title,  kRed,   1, 0,   Form("R_{m_{T}} (M_{J}<%i GeV)", MJ_));
    h1cosmetic(h1_R31_data,     Title,  kBlack, 1, 0,   Form("R_{m_{T}} (M_{J}<%i GeV)", MJ_));
    h1cosmetic(h1_R42_MC,       Title,  kRed,   1, 0,   Form("R_{m_{T}} (M_{J}>%i GeV)", MJ_)); 
    h1cosmetic(h1_R42_data,     Title,  kBlack, 1, 0,   Form("R_{m_{T}} (M_{J}>%i GeV)", MJ_)); 
    h1cosmetic(h1_kappa_MC,     Title,  kBlue,  1, 0,   "#kappa"); 
    h1cosmetic(h1_kappa_Data,   Title,  kBlack, 1, 0,   "#kappa"); 
    
    TCanvas *c = new TCanvas("c","c",900,600);
    c->Divide(3,2);
    c->cd(1);
    h1_R21_MC->SetMaximum(1);
    h1_R21_MC->Draw("e");
    h1_R21_data->Draw("same e");
    c->cd(2);
    h1_R43_MC->SetMaximum(1);
    h1_R43_MC->Draw("e");
    h1_R43_data->Draw("same e");
    c->cd(3);
    h1_kappa_MC->SetMinimum(0);
    h1_kappa_MC->SetMaximum(5);
    h1_kappa_MC->Draw("e");
    h1_kappa_Data->Draw("same e");
    c->cd(4);
    h1_R31_MC->SetMaximum(0.3);
    h1_R31_MC->Draw("e");
    h1_R31_data->Draw("same e");
    c->cd(5);
    h1_R42_MC->SetMaximum(0.3);
    h1_R42_MC->Draw("e");
    h1_R42_data->Draw("same e");
    c->Print(Form("figR/R_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i.pdf",METlow_,METhigh_,njlow_,njhigh_,nblow_,nbhigh_,MJ_));
    delete c;
}

void DrawR()
{ 
    //DrawRone(int METlow_=100, int METhigh_=200, int njlow_=6, int njhigh_=99, int nblow_=2, int nbhigh_=9, int MJ_=500)
    /* 
    DrawRone(100,200,5,99,2,9,400);
    DrawRone(200,300,5,99,2,9,400);
    
    DrawRone(100,200,6,99,1,2,400);
    DrawRone(200,300,6,99,1,2,400);
    */
    /*
    DrawRone(100,200,5,99,2,99,500);
    DrawRone(200,300,5,99,2,99,500);
    
    DrawRone(100,200,6,99,1,1,500);
    DrawRone(200,300,6,99,1,1,500);
    */ 
//    DrawRSummary(1);
    DrawRSummary(2, "RMJ");
    DrawRSummary(2, "RmT");
    DrawKappaSummary(2);
}
