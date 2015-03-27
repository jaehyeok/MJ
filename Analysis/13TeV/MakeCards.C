#include <iostream>
#include <fstream>
#include <iomanip> // for setw()

#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TInterpreter.h"
#include "TLatex.h"
#include "TMath.h"

ofstream fout;

//
// Print A line of table
//
void PrintCardOneLine(std::vector<TString> col1to2, std::vector<TString> col3to)
{
    TString line1to2 = Form("%20s%7s",col1to2.at(0).Data(),col1to2.at(1).Data()); 
    TString line3to=line1to2; 
    for(int i=0; i<(int)col3to.size(); i++) line3to = line3to+Form("%10s",col3to.at(i).Data()); 
    cout << line3to << endl;   
}
const char* PrintCardOneLine(const char* col1, const char* col2, const char* col3, const char* col4, const char* col5, const char* col6, const char* col7)
{
    const char* line = Form("%20s%7s%10s%10s%10s%10s%10s", col1,col2,col3,col4,col5,col6,col7); 
    return line;
}
const char* PrintCardOneLine(const char* col1, const char* col2, float col3, float col4, float col5, float col6, float col7)
{
    const char* line = Form("%20s%7s%7.3f%7.3f%7.3f%7.3f%7.3f", col1,col2,col3,col4,col5,col6,col7); 
    return line;
}

void MakeCards(int lepflav=0, const char* Region="")
{ 
    cout << "[MJ Analysis] Make cards for " << Region << endl; 
    TString lepflavname = "all"; 
    if(lepflav==11) lepflavname = "e"; 
    if(lepflav==13) lepflavname = "m"; 

    if(lepflav==0)  cout << "[MJ Analysis] Card for Electron+Muon" << endl;
    if(lepflav==11) cout << "[MJ Analysis] Card for Electron" << endl;
    if(lepflav==13) cout << "[MJ Analysis] Card for Muon" << endl;

    TString HistName="yields";

    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_%s.root", Region));
        
    TH1F *h1_DATA[7], *h1_T[7], *h1_TT_sl[7], *h1_TT_ll[7], *h1_WJets[7], *h1_DY[7], *h1_MC[7];
    TH1F *h1_f1500_100[7], *h1_f1200_800[7];
    for(int i=2; i<7; i++)
    {
        h1_DATA[i]      = (TH1F*)HistFile->Get(Form("h1_DATA_%s_%ifatjet", HistName.Data(), i));
        h1_T[i]         = (TH1F*)HistFile->Get(Form("h1_T_%s_%ifatjet", HistName.Data(), i));
        h1_TT_sl[i]     = (TH1F*)HistFile->Get(Form("h1_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h1_TT_ll[i]     = (TH1F*)HistFile->Get(Form("h1_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h1_WJets[i]     = (TH1F*)HistFile->Get(Form("h1_WJets_%s_%ifatjet", HistName.Data(), i));
        h1_DY[i]        = (TH1F*)HistFile->Get(Form("h1_DY_%s_%ifatjet", HistName.Data(), i));
        h1_f1500_100[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1500_100_%s_%ifatjet", HistName.Data(), i));
        h1_f1200_800[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1200_800_%s_%ifatjet", HistName.Data(), i));

        h1_MC[i] = (TH1F*)h1_TT_sl[i]->Clone("h1_MC");
        h1_MC[i]->Add(h1_TT_ll[i]);
        h1_MC[i]->Add(h1_WJets[i]);
        h1_MC[i]->Add(h1_T[i]);
        h1_MC[i]->Add(h1_DY[i]);
    }   

    //
    // yields in string format 
    //
    const char* Ndata;          
    const char* Nsig;          
    const char* Nttbar;        
    const char* Nsingletop;    
    const char* Nwjets;       
    const char* Ndy;           
    if(lepflav!=11 && lepflav!=13)
    {
        Ndata         = Form("%i",   (int)(h1_MC[6]->Integral(0,1000)+h1_f1500_100[6]->Integral(0,1000))); 
        Nsig          = Form("%.3f", h1_f1500_100[6]->Integral(0,1000)); 
        Nttbar        = Form("%.3f", h1_TT_ll[6]->Integral(0,1000)+h1_TT_sl[6]->Integral(0,1000)); 
        Nsingletop    = Form("%.3f", h1_T[6]->Integral(0,1000)); 
        Nwjets        = Form("%.3f", h1_WJets[6]->Integral(0,1000)); 
        Ndy           = Form("%.3f", h1_DY[6]->Integral(0,1000)); 
    }
    else 
    { 
        int bin = (lepflav-9)/2;
        Ndata         = Form("%i",   (int)(h1_MC[6]->GetBinContent(bin)+h1_f1500_100[6]->GetBinContent(bin))); 
        Nsig          = Form("%.3f", h1_f1500_100[6]->GetBinContent(bin)); 
        Nttbar        = Form("%.3f", h1_TT_ll[6]->GetBinContent(bin)+h1_TT_sl[6]->GetBinContent(bin)); 
        Nsingletop    = Form("%.3f", h1_T[6]->GetBinContent(bin)); 
        Nwjets        = Form("%.3f", h1_WJets[6]->GetBinContent(bin)); 
        Ndy           = Form("%.3f", h1_DY[6]->GetBinContent(bin)); 
    }

    // -------------------------------------
    // Print out on file  
    // -------------------------------------
    fout.open(Form("Cards/%s_%s.dat", Region, lepflavname.Data()));

    fout << "imax 1 number of channels" << endl;
    fout << "jmax * number of background" << endl;
    fout << "kmax * number of nuisance parameters" << endl;
    fout << "Observation " << Ndata << endl;
    fout << PrintCardOneLine("bin",                         "",     Region,     Region,     Region,         Region,     Region  ) << endl;
    fout << PrintCardOneLine("process",                     "",     "T1tttt",   "ttbar",    "singletop",    "DY",       "Wjets" ) << endl;
    fout << PrintCardOneLine("process",                     "",     "0",        "1",        "2",            "3",        "4"     ) << endl;
    fout << PrintCardOneLine("rate",                        "",     Nsig,       Nttbar,     Nsingletop,     Nwjets,     Ndy     ) << endl;
    fout << PrintCardOneLine(Form("ttbar_SF_%s",Region),    "lnN",  "-",       "1.5",       "-",            "-",        "-"     ) << endl;

    fout.close();

}

