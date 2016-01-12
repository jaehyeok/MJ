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

using namespace std;
ofstream fout;

TString region_title[24] = {
    "R1: 200<MET<400, 6#leq njets#leq8, nb#geq1",	
 	"R2: 200<MET<400, 6#leq njets#leq8, nb=1",	
 	"R2: 200<MET<400, njets#geq9, nb=1",	
 	"R2: 200<MET<400, 6#leq njets#leq8, nb=2",	
 	"R2: 200<MET<400, njets#geq9, nb=2",	
 	"R2: 200<MET<400, 6#leq njets#leq8, nb#geq3",	
 	"R2: 200<MET<400, njets#geq9, nb#geq3",	
 	"R3: 200<MET<400, 6#leq njets#leq8, nb#geq1",	
 	"R4: 200<MET<400, 6#leq njets#leq8, nb=1",	
 	"R4: 200<MET<400, njets#geq9, nb=1",	
 	"R4: 200<MET<400, 6#leq njets#leq8, nb=2",	
 	"R4: 200<MET<400, njets#geq9, nb=2",	
 	"R4: 200<MET<400, 6#leq njets#leq8, nb#geq3",	
 	"R4: 200<MET<400, njets#geq9, nb#geq3",	
 	"R1: MET>400, 6#leq njets#leq8, nb#geq1",	
 	"R2: MET>400, 6#leq njets#leq8, nb=1",	
 	"R2: MET>400, njets#geq9, nb=1",	
 	"R2: MET>400, 6#leq njets#leq8, nb#geq2",	
 	"R2: MET>400, njets#geq9, nb#geq2",	
 	"R3: MET>400, 6#leq njets#leq8, nb#geq1",	
 	"R4: MET>400, 6#leq njets#leq8, nb=1",	
 	"R4: MET>400, njets#geq9, nb=1",	
 	"R4: MET>400, 6#leq njets#leq8, nb#geq2",	
 	"R4: MET>400, njets#geq9, nb#geq2"
};	

//
// Print A line of table
//
void PrintTableOneLine(TString Process, TH1F* h1[7], int lepflav=0, bool doLatex=true, bool PrintFile=true)
{
    // Print out on file 
    if(PrintFile) 
    {
        if(lepflav!=11 && lepflav!=13)
        {
            Double_t error[7];
            for(int i=2; i<7; i++) h1[i]->IntegralAndError(1,10000,error[i]);
            fout << Process << " & " 
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
            fout << Process << " & " 
                 << Form("$%.2f \\pm %.2f$",h1[2]->GetBinContent(bin),h1[2]->GetBinError(bin)) << "\t&" 
                 << Form("$%.2f \\pm %.2f$",h1[3]->GetBinContent(bin),h1[3]->GetBinError(bin)) << "\t&" 
                 << Form("$%.2f \\pm %.2f$",h1[4]->GetBinContent(bin),h1[4]->GetBinError(bin)) << "\t&" 
                 << Form("$%.2f \\pm %.2f$",h1[5]->GetBinContent(bin),h1[5]->GetBinError(bin)) << "\t&" 
                 << Form("$%.2f \\pm %.2f$",h1[6]->GetBinContent(bin),h1[6]->GetBinError(bin)) << "\t\\\\" 
                 << endl; 
        }
        if(Process.Contains("Bkg")) fout << "\\hline \\hline" << endl;  
    }
    else 
    {
        // Print out on screen  
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
        // Print out on screen  
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
}

void MakeTables(int lepflav=0, const char* Selection="", bool doLatex=false, float Lumi=40)
{ 
    if(lepflav==0)  cout << "[MJ Table] Yields for Electron+Muon" << endl;
    if(lepflav==11) cout << "[MJ Table] Yields for Electron" << endl;
    if(lepflav==13) cout << "[MJ Table] Yields for Muon" << endl;

    TString HistName="yields";

    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_%s.root", Selection));
        
    TH1F *h1_DATA[7], *h1_T[7], *h1_TT_sl[7], *h1_TT_ll[7], *h1_WJets[7], *h1_DY[7], *h1_TTV[7], *h1_QCD[7], *h1_Others[7], *h1_MC[7];
    TH1F *h1_f1500_100[7], *h1_f1200_800[7];
    for(int i=2; i<7; i++)
    {
        h1_DATA[i]      = (TH1F*)HistFile->Get(Form("h1_DATA_%s_%ifatjet",   HistName.Data(), i));
        h1_T[i]         = (TH1F*)HistFile->Get(Form("h1_T_%s_%ifatjet",      HistName.Data(), i));
        h1_TT_sl[i]     = (TH1F*)HistFile->Get(Form("h1_TT_sl_%s_%ifatjet",  HistName.Data(), i));
        h1_TT_ll[i]     = (TH1F*)HistFile->Get(Form("h1_TT_ll_%s_%ifatjet",  HistName.Data(), i));
        h1_WJets[i]     = (TH1F*)HistFile->Get(Form("h1_WJets_%s_%ifatjet",  HistName.Data(), i));
        h1_DY[i]        = (TH1F*)HistFile->Get(Form("h1_DY_%s_%ifatjet",     HistName.Data(), i));
        h1_TTV[i]       = (TH1F*)HistFile->Get(Form("h1_TTV_%s_%ifatjet",    HistName.Data(), i));
        h1_QCD[i]       = (TH1F*)HistFile->Get(Form("h1_QCD_%s_%ifatjet",    HistName.Data(), i));
        h1_Others[i]    = (TH1F*)HistFile->Get(Form("h1_Others_%s_%ifatjet", HistName.Data(), i));
        h1_f1500_100[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1500_100_%s_%ifatjet", HistName.Data(), i));
        h1_f1200_800[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1200_800_%s_%ifatjet", HistName.Data(), i));

        h1_MC[i] = (TH1F*)h1_TT_sl[i]->Clone("h1_MC");
        h1_MC[i]->Add(h1_TT_ll[i]);
        h1_MC[i]->Add(h1_WJets[i]);
        h1_MC[i]->Add(h1_T[i]);
        h1_MC[i]->Add(h1_DY[i]);
        h1_MC[i]->Add(h1_TTV[i]);
        h1_MC[i]->Add(h1_QCD[i]);
        h1_MC[i]->Add(h1_Others[i]);
    }   

    // -------------------------------------
    // Print out on screen 
    // -------------------------------------
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

    PrintTableOneLine("TT(l)",              h1_TT_sl,       lepflav,	doLatex, false);
    PrintTableOneLine("TT(ll)",             h1_TT_ll,       lepflav,	doLatex, false);
    PrintTableOneLine("t+tW",               h1_T,           lepflav,	doLatex, false);
    PrintTableOneLine("WJets",              h1_WJets,       lepflav,	doLatex, false);
    PrintTableOneLine("DY",                 h1_DY,          lepflav,	doLatex, false);
    PrintTableOneLine("TTV",                h1_TTV,         lepflav,	doLatex, false);
    PrintTableOneLine("QCD",                h1_QCD,         lepflav,	doLatex, false);
    PrintTableOneLine("Others",             h1_Others,      lepflav,	doLatex, false);
    PrintTableOneLine("Total Bkg",          h1_MC,          lepflav,	doLatex, false);
    PrintTableOneLine("Data",               h1_DATA,        lepflav,	doLatex, false);
    PrintTableOneLine("T1tttt[1200,800]",   h1_f1200_800,   lepflav,	doLatex, false);
    PrintTableOneLine("T1tttt[1500,100]",   h1_f1500_100,   lepflav,	doLatex, false);
        
    if(doLatex)
    {
        cout << "\\hline \\hline" << endl;
        cout << "\\end{tabular}" << endl;
        cout << "\\label{tab:"<< Selection <<"}" << endl;
        cout << "\\caption{" << Selection << "}" << endl;
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
    
    // -------------------------------------
    // Print out on file  
    // -------------------------------------
    fout.open(Form("Tables/TableOfYields_%s.tex", Selection));

    fout << "\\begin{table}[!htb]" << endl;
    fout << "\\centering" << endl;
    fout << "\\begin{tabular}{c | c c c c | c }" << endl;
    fout << "\\hline \\hline" << endl;
    fout << "Process  & " 
         << "$N_{FJ}=2$ \t&" 
         << "$N_{FJ}=3$ \t&" 
         << "$N_{FJ}=4$ \t&" 
         << "$N_{FJ}\\ge 5$ \t&" 
         << "$N_{FJ}\\ge 2$ \t\\\\" 
         << endl; 
    fout << "\\hline" << endl;
    fout << "\\multicolumn{6}{c}{Electron + muon channel} \\\\" << endl; 
    fout << "\\hline" << endl;
    lepflav=0;
    PrintTableOneLine("TT(l)",              h1_TT_sl,       lepflav,	doLatex, true);
    PrintTableOneLine("TT(ll)",             h1_TT_ll,       lepflav,	doLatex, true);
    PrintTableOneLine("t+tW",               h1_T,           lepflav,	doLatex, true);
    PrintTableOneLine("WJets",              h1_WJets,       lepflav,	doLatex, true);
    PrintTableOneLine("DY",                 h1_DY,          lepflav,	doLatex, true);
    PrintTableOneLine("TTV",                h1_TTV,         lepflav,	doLatex, true);
    PrintTableOneLine("QCD",                h1_QCD,         lepflav,	doLatex, true);
    PrintTableOneLine("Others",             h1_Others,      lepflav,	doLatex, true);
    PrintTableOneLine("Total Bkg",          h1_MC,          lepflav,	doLatex, true);
    PrintTableOneLine("Data",               h1_DATA,        lepflav,	doLatex, true);
    PrintTableOneLine("T1tttt[1200,800]",   h1_f1200_800,   lepflav,	doLatex, true);
    PrintTableOneLine("T1tttt[1500,100]",   h1_f1500_100,   lepflav,	doLatex, true);
    fout << "\\hline \\hline" << endl;
    
    lepflav=11;
    fout << "\\multicolumn{6}{c}{Electron channel} \\\\" << endl; 
    fout << "\\hline" << endl;
    PrintTableOneLine("TT(l)",              h1_TT_sl,       lepflav,	doLatex, true);
    PrintTableOneLine("TT(ll)",             h1_TT_ll,       lepflav,	doLatex, true);
    PrintTableOneLine("t+tW",               h1_T,           lepflav,	doLatex, true);
    PrintTableOneLine("WJets",              h1_WJets,       lepflav,	doLatex, true);
    PrintTableOneLine("DY",                 h1_DY,          lepflav,	doLatex, true);
    PrintTableOneLine("TTV",                h1_TTV,         lepflav,	doLatex, true);
    PrintTableOneLine("QCD",                h1_QCD,         lepflav,	doLatex, true);
    PrintTableOneLine("Others",             h1_Others,      lepflav,	doLatex, true);
    PrintTableOneLine("Total Bkg",          h1_MC,          lepflav,	doLatex, true);
    PrintTableOneLine("Data",               h1_DATA,        lepflav,	doLatex, true);
    PrintTableOneLine("T1tttt[1200,800]",   h1_f1200_800,   lepflav,	doLatex, true);
    PrintTableOneLine("T1tttt[1500,100]",   h1_f1500_100,   lepflav,	doLatex, true);
    fout << "\\hline \\hline" << endl;
    
    lepflav=13;
    fout << "\\multicolumn{6}{c}{Muon channel} \\\\" << endl; 
    fout << "\\hline" << endl;
    PrintTableOneLine("TT(l)",              h1_TT_sl,       lepflav,	doLatex, true);
    PrintTableOneLine("TT(ll)",             h1_TT_ll,       lepflav,	doLatex, true);
    PrintTableOneLine("t+tW",               h1_T,           lepflav,	doLatex, true);
    PrintTableOneLine("WJets",              h1_WJets,       lepflav,	doLatex, true);
    PrintTableOneLine("DY",                 h1_DY,          lepflav,	doLatex, true);
    PrintTableOneLine("TTV",                h1_TTV,         lepflav,	doLatex, true);
    PrintTableOneLine("QCD",                h1_QCD,         lepflav,	doLatex, true);
    PrintTableOneLine("Others",             h1_Others,      lepflav,	doLatex, true);
    PrintTableOneLine("Total Bkg",          h1_MC,          lepflav,	doLatex, true);
    PrintTableOneLine("Data",               h1_DATA,        lepflav,	doLatex, true);
    PrintTableOneLine("T1tttt[1200,800]",   h1_f1200_800,   lepflav,	doLatex, true);
    PrintTableOneLine("T1tttt[1500,100]",   h1_f1500_100,   lepflav,	doLatex, true);

    // Print out on file  
    fout << "\\hline \\hline" << endl;
    fout << "\\end{tabular}" << endl;
    fout << "\\label{tab:"<< Selection << "}" << endl;
    fout << "\\caption{Table of yields in " << Selection << ". From top to bottom, electon+muon, electron and muon channels}" << endl;
    fout << "\\end{table}" << endl;
    fout.close();
} 

void MakeTablesAllRegions(int lepflav=0, const char* Selection="", bool doLatex=false, float Lumi=40)
{ 
    if(lepflav==0)  cout << "[MJ Table] Yields for Electron+Muon" << endl;
    if(lepflav==11) cout << "[MJ Table] Yields for Electron" << endl;
    if(lepflav==13) cout << "[MJ Table] Yields for Muon" << endl;

    TString HistName="yields_bins";

    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_%s.root", Selection));
        
    TH2F *h2_DATA[7], *h2_T[7], *h2_TT_sl[7], *h2_TT_ll[7], *h2_WJets[7], *h2_DY[7], *h2_TTV[7], *h2_QCD[7], *h2_Others[7], *h2_MC[7];
    TH2F *h2_f1500_100[7], *h2_f1200_800[7];
    for(int i=0; i<7; i++)
    {
        h2_DATA[i]      = (TH2F*)HistFile->Get(Form("h2_DATA_%s_%ifatjet",   HistName.Data(), i));
        h2_T[i]         = (TH2F*)HistFile->Get(Form("h2_T_%s_%ifatjet",      HistName.Data(), i));
        h2_TT_sl[i]     = (TH2F*)HistFile->Get(Form("h2_TT_sl_%s_%ifatjet",  HistName.Data(), i));
        h2_TT_ll[i]     = (TH2F*)HistFile->Get(Form("h2_TT_ll_%s_%ifatjet",  HistName.Data(), i));
        h2_WJets[i]     = (TH2F*)HistFile->Get(Form("h2_WJets_%s_%ifatjet",  HistName.Data(), i));
        h2_DY[i]        = (TH2F*)HistFile->Get(Form("h2_DY_%s_%ifatjet",     HistName.Data(), i));
        h2_TTV[i]       = (TH2F*)HistFile->Get(Form("h2_TTV_%s_%ifatjet",    HistName.Data(), i));
        h2_QCD[i]       = (TH2F*)HistFile->Get(Form("h2_QCD_%s_%ifatjet",    HistName.Data(), i));
        h2_Others[i]    = (TH2F*)HistFile->Get(Form("h2_Others_%s_%ifatjet", HistName.Data(), i));
        h2_f1500_100[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1500_100_%s_%ifatjet", HistName.Data(), i));
        h2_f1200_800[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1200_800_%s_%ifatjet", HistName.Data(), i));

        h2_MC[i] = (TH2F*)h2_TT_sl[i]->Clone("h2_MC");
        h2_MC[i]->Add(h2_TT_ll[i]);
        h2_MC[i]->Add(h2_WJets[i]);
        h2_MC[i]->Add(h2_T[i]);
        h2_MC[i]->Add(h2_DY[i]);
        h2_MC[i]->Add(h2_TTV[i]);
        h2_MC[i]->Add(h2_QCD[i]);
        h2_MC[i]->Add(h2_Others[i]);
    }   

    // -------------------------------------
    // Print out on screen 
    // -------------------------------------
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
            setw(35) << Form("-----------------------------------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(5)  << Form("-----") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "-|" << endl;
        cout << "|" << 
        /*
            setw(10) << Form("%10s","Region")    << " |"  <<
            setw(10) << Form("%10s","ttbar(sl)") << " |" <<
            setw(10) << Form("%10s","ttbar(ll)") << " |" <<
            setw(10) << Form("%10s","w+jets")    << " |" <<
            setw(10) << Form("%10s","dy")        << " |" <<
            setw(10) << Form("%10s","single top")<< " |" << 
            setw(10) << Form("%10s","ttv")       << " |" << 
            setw(10) << Form("%10s","qcd")       << " |" << 
            setw(10) << Form("%10s","others")    << " |" << 
            setw(10) << Form("%10s","1500,100")  << " |" << 
            setw(10) << Form("%10s","1200,800")  << " |" << endl;
        */
            setw(35) << "Region"    << " |"  <<
            setw(10) << "ttbar(sl)" << " |" <<
            setw(10) << "ttbar(ll)" << " |" <<
            setw(10) << "w+jets"    << " |" <<
            setw(10) << "dy"        << " |" <<
            setw(10) << "single top"<< " |" << 
            setw(10) << "ttv"       << " |" << 
            setw(10) << "qcd"       << " |" << 
            setw(10) << "others"    << " |" << 
            setw(10) << "all bkg"   << " |" << 
            setw(5) << "data"      << " |" << 
            setw(10) << "1500,100"  << " |" << 
            setw(10) << "1200,800"  << " |" << endl;
        cout << "|" <<
            setw(35) << Form("-----------------------------------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(5)  << Form("-----") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "-|" << endl;
    }

    for(int iregion=0; iregion<24; iregion++)  
    {
        cout << "|" <<   
            left << setw(35) << region_title[iregion].ReplaceAll("#geq",">=").ReplaceAll("#leq ","<=").ReplaceAll("#leq","<=").ReplaceAll("MET>400,","MET>400,    ").ReplaceAll("njets>=9,","njets>=9,   ")   << " |"  <<
            left << setw(10) << Form("%.2f",h2_TT_sl[6]->GetBinContent(3,iregion+1))  << " |" << 
            left << setw(10) << Form("%.2f",h2_TT_ll[6]->GetBinContent(3,iregion+1))  << " |" << 
            left << setw(10) << Form("%.2f",h2_WJets[6]->GetBinContent(3,iregion+1))  << " |" << 
            left << setw(10) << Form("%.2f",h2_DY[6]->GetBinContent(3,iregion+1))  << " |" << 
            left << setw(10) << Form("%.2f",h2_T[6]->GetBinContent(3,iregion+1))      << " |" << 
            left << setw(10) << Form("%.2f",h2_TTV[6]->GetBinContent(3,iregion+1))    << " |" << 
            left << setw(10) << Form("%.2f",h2_QCD[6]->GetBinContent(3,iregion+1))    << " |" << 
            left << setw(10) << Form("%.2f",h2_Others[6]->GetBinContent(3,iregion+1)) << " |" << 
            left << setw(10) << Form("%.2f",h2_MC[6]->GetBinContent(3,iregion+1)) << " |" << 
            left << setw(5)  << Form("%.0f",h2_DATA[6]->GetBinContent(3,iregion+1)) << " |" << 
            left << setw(10) << Form("%.2f",h2_f1500_100[6]->GetBinContent(3,iregion+1)) << " |" << 
            left << setw(10) << Form("%.2f",h2_f1200_800[6]->GetBinContent(3,iregion+1)) << " |" << endl;
    } 

    if(doLatex)
    {
        cout << "\\hline \\hline" << endl;
        cout << "\\end{tabular}" << endl;
        cout << "\\label{tab:method2yields" << endl;
        cout << "\\caption{Yields in all regions.} " << endl;
        cout << "\\end{table}" << endl;
    }
    else 
    { 
        cout << "|" <<
            setw(35) << Form("-----------------------------------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(5)  << Form("-----") << "--" <<
            setw(10) << Form("----------") << "--" <<
            setw(10) << Form("----------") << "-|" << endl;
    }
    cout << endl;
    
} 

