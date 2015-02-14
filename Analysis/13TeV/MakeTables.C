#include <iostream>
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

//
// Print A line of table
//
void PrintTableOneLine(TString Process, TH1F* h1[7], int lepflav=0, bool doLatex=true)
{ 
        if(doLatex)
        {
            if(lepflav!=11 && lepflav!=13)
            {
                Double_t error[7];
                for(int i=2; i<7; i++) h1[i]->IntegralAndError(1,10000,error[i]);
                cout << Process << " & " 
                     << Form("$%.2f \\pm %.2f$",h1[2]->Integral(),error[2]) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[3]->Integral(),error[3]) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[4]->Integral(),error[4]) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[5]->Integral(),error[5]) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[6]->Integral(),error[6]) << "\t\\\\" 
                     << endl; 
            }
            else 
            {   
                int bin = (lepflav-9)/2;
                cout << Process << " & " 
                     << Form("$%.2f \\pm %.2f$",h1[2]->GetBinContent(bin),h1[2]->GetBinError(bin)) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[3]->GetBinContent(bin),h1[3]->GetBinError(bin)) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[4]->GetBinContent(bin),h1[4]->GetBinError(bin)) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[5]->GetBinContent(bin),h1[5]->GetBinError(bin)) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[6]->GetBinContent(bin),h1[6]->GetBinError(bin)) << "\t\\\\" 
                     << endl; 
            }
        }
        else 
        {   
            if(lepflav!=11 && lepflav!=13)
            {
                Double_t error[7];
                for(int i=2; i<7; i++) h1[i]->IntegralAndError(1,10000,error[i]);
                cout << "|" << 
                    setw(20) << Process << " |"  <<
                    setw(20) << Form("%.2f +/- %.2f", h1[2]->Integral(), error[2]) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[3]->Integral(), error[3]) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[4]->Integral(), error[4]) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[5]->Integral(), error[5]) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[6]->Integral(), error[6]) << " |" << endl; 
            }
            else 
            {   
                int bin = (lepflav-9)/2;
                cout << "|" << 
                    setw(20) << Process << " |"  <<
                    setw(20) << Form("%.2f +/- %.2f", h1[2]->GetBinContent(bin), h1[2]->GetBinError(bin)) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[3]->GetBinContent(bin), h1[3]->GetBinError(bin)) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[4]->GetBinContent(bin), h1[4]->GetBinError(bin)) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[5]->GetBinContent(bin), h1[5]->GetBinError(bin)) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[6]->GetBinContent(bin), h1[6]->GetBinError(bin)) << " |" << endl; 
            }
        }
        
        // 
        if(Process.Contains("Bkg")) 
        { 
            if(doLatex) 
            {   
                cout << "\\hline \\hline" << endl;  
            }
            else 
            { 
                cout << "|" << 
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "-|" << endl;
            }
        }
}

void MakeTables(int lepflav=0, bool doLatex=false)
{ 
    if(lepflav==0)  cout << "[MJ Table] Yields for Electron+Muon" << endl;
    if(lepflav==11) cout << "[MJ Table] Yields for Electron" << endl;
    if(lepflav==13) cout << "[MJ Table] Yields for Muon" << endl;

    TString HistName="yields";

    TFile* HistFile = TFile::Open("HistFiles/Hist.root");
        
    TH1F *h1_DATA[7], *h1_T[7], *h1_TT_sl[7], *h1_TT_ll[7], *h1_TT[7], *h1_WJets[7], *h1_DY[7], *h1_MC[7];
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

    // print header
    if(doLatex)
    {

        cout << "\\begin{table}[!htb]" << endl;
        cout << "\\centering" << endl;
        cout << "\\begin{tabular}{c | c c c c | c }" << endl;
        cout << "\\hline \\hline" << endl;
        cout << "Process  & " 
             << "$N_{FJ}=2$ \t&" 
             << "$N_{FJ}=3$ \t&" 
             << "$N_{FJ}=4$ \t&" 
             << "$N_{FJ}\\ge 5$ \t&" 
             << "$N_{FJ}\\ge 2$ \t\\\\" 
             << endl; 
        cout << "\\hline" << endl;
    }
    else 
    { 
        cout << "|" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "-|" << endl;
        cout << "|" << 
            setw(20) << "Process" << " |"  <<
            setw(20) << "N_FJ=2" << " |" <<
            setw(20) << "N_FJ=3" << " |" <<
            setw(20) << "N_FJ=4" << " |" <<
            setw(20) << "N_FJ>=5" << " |" <<
            setw(20) << "N_FJ>=2" << " |" << endl;
        cout << "|" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "-|" << endl;
    }

    PrintTableOneLine("TT(l)",              h1_TT_sl,       lepflav,	doLatex);
    PrintTableOneLine("TT(ll)",             h1_TT_ll,       lepflav,	doLatex);
    PrintTableOneLine("t+tW",               h1_T,           lepflav,	doLatex);
    PrintTableOneLine("WJets",              h1_WJets,       lepflav,	doLatex);
    PrintTableOneLine("DY",                 h1_DY,          lepflav,	doLatex);
    PrintTableOneLine("Total Bkg",          h1_MC,          lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1200,800]",   h1_f1200_800,   lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1500,100]",   h1_f1500_100,   lepflav,	doLatex);

    if(doLatex)
    {
        cout << "\\hline \\hline" << endl;
        cout << "\\end{tabular}" << endl;
        //cout << "\\label{tab:exp_sig_0j}" << endl;
        cout << "\\caption{Yields}" << endl;
        cout << "\\end{table}" << endl;
    }
    else 
    { 
        cout << "|" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "-|" << endl;
    }

    cout << endl;
}
