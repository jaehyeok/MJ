#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "TROOT.h"
#include "TString.h"
#include "TH1F.h"
#include "THStack.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"

using namespace std;
bool DoRebin        = 0;
bool OneFileTest    = 0;
int pTR0p5     = 30;

//
//TH1F initialization
//
TH1F* InitTH1F(char* Name, char* Title, int Nbins, double XMin, double XMax){
    TH1F *h1 = new TH1F(Name, Title, Nbins, XMin, XMax);
    h1->Sumw2();
    return h1;
}

//
// h1 cosmetics
//
void h1cosmetic(TH1F* &h1, char* title, int linecolor=kBlack, int linewidth=1, int fillcolor=0, TString var="")
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.1);
}

//
// loading PU reweigting histogram
//
double nPUScaleFactor2012(TH1F* h1PU, float npu){
    // Get PU scale factor histogram
    double mynpu = TMath::Min(npu,(float)49.499);
    Int_t npuxbin = h1PU->GetXaxis()->FindBin(mynpu);
    return h1PU->GetBinContent(npuxbin);
}

//
// main
//
void DrawStack(TString HistName, bool doZ, bool RemoveMuon, int pTR0p5) 
{ 
    TString HistDir = "HistFiles";
    if(!doZ) HistDir = HistDir+"_Zveto";
    if(RemoveMuon) HistDir = HistDir+"_MuonRemoved";
    
    TFile* HistFile = TFile::Open(Form("%s/Hist_pT%i.root", HistDir.Data(), pTR0p5));
     
    char *var; 
    bool skip0fatjet = false;
    if(HistName=="MJ")      { var=(char*)"M_{J} [GeV]"; skip0fatjet=true; }
    if(HistName=="mj")      { var=(char*)"m_{j} [GeV]"; skip0fatjet=true; }
    if(HistName=="mll")     var=(char*)"m_{ll} [GeV]";
    if(HistName=="Nfatjet") var=(char*)"N_{fatjet}";

    TH1F *h1_DATA[6], *h1_DY[6], *h1_TT[6], *h1_WJets[6], *h1_QCD[6], *h1_MC[6], *h1_Ratio[6]; 
    TH1F *h1_One[6];
    THStack *st[6];
    TCanvas *c = new TCanvas("c","c",900,600*1.2);  
    c->Divide(3,2);
    for(int i=skip0fatjet?1:0; i<6; i++) 
    {
        h1_DATA[i]  = (TH1F*)HistFile->Get(Form("h1_DATA_%s_%ifatjet", HistName.Data(), i)); 
        h1_DY[i]    = (TH1F*)HistFile->Get(Form("h1_DY_%s_%ifatjet", HistName.Data(), i));
        h1_TT[i]    = (TH1F*)HistFile->Get(Form("h1_TT_%s_%ifatjet", HistName.Data(), i));
        h1_WJets[i] = (TH1F*)HistFile->Get(Form("h1_WJets_%s_%ifatjet", HistName.Data(), i));
        h1_QCD[i]   = (TH1F*)HistFile->Get(Form("h1_QCD_%s_%ifatjet", HistName.Data(), i)); 
        h1_One[i]   = InitTH1F(Form("h1_One_%i",i), Form("h1_One_%i",i), 1, 
                               h1_QCD[i]->GetXaxis()->GetBinLowEdge(1), 
                               h1_QCD[i]->GetXaxis()->GetBinUpEdge(h1_QCD[i]->GetXaxis()->GetNbins())  ); 
        h1_One[i]->SetBinContent(1,1);

        if(DoRebin) {
            int NMergeBins = 2;
            h1_DATA[i]->Rebin(NMergeBins);
            h1_DY[i]->Rebin(NMergeBins);
            h1_TT[i]->Rebin(NMergeBins);
            h1_WJets[i]->Rebin(NMergeBins);
            h1_QCD[i]->Rebin(NMergeBins);
        }
        
        h1_MC[i] = (TH1F*)h1_DY[i]->Clone(Form("h1_MC_%s_%ifatjet", HistName.Data(), i));
        h1_MC[i]->Add(h1_TT[i]);
        h1_MC[i]->Add(h1_WJets[i]);
        h1_MC[i]->Add(h1_QCD[i]);
        h1_Ratio[i] = (TH1F*)h1_DATA[i]->Clone(Form("h1_Ratio_%s_%ifatjet", HistName.Data(), i));
        h1_Ratio[i]->Divide(h1_MC[i]);

        h1cosmetic(h1_DATA[i],  Form("DATA %ifatjet", i),   kBlack, 1, 0,       var);
        h1cosmetic(h1_DY[i],    Form("DY %ifatjet", i),     kBlack, 1, kYellow, var);
        h1cosmetic(h1_TT[i],    Form("TT %ifatjet", i),     kBlack, 1, kOrange, var);
        h1cosmetic(h1_WJets[i], Form("WJets %ifatjet", i),  kBlack, 1, kGreen,  var);
        h1cosmetic(h1_QCD[i],   Form("QCD %ifatjet", i),    kBlack, 1, kBlue,   var);
        h1cosmetic(h1_Ratio[i], Form(""),                   kBlack, 1, kBlack,  var);
        h1cosmetic(h1_One[i],   Form(""),                   kBlue,  1, 0,       var);
  
        c->cd(i+1);
        c->cd(i+1)->SetLogy(1);

        TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.04);
        pad1->SetRightMargin(0.07);
        pad1->SetLogy(1);
        pad1->Draw();
        pad1->cd();
        st[i] = new THStack( Form("Stack %ifatjet", i), 
                             i==5?"All fatjets":Form("%i fatjets; ; ", i));
        st[i]->Add(h1_QCD[i]);
        st[i]->Add(h1_WJets[i]);
        st[i]->Add(h1_TT[i]);
        st[i]->Add(h1_DY[i]); 
        st[i]->SetMaximum(h1_DATA[i]->GetMaximum()*1.5);
        st[i]->SetMinimum(0.1);
        st[i]->Draw("HIST"); 

        h1_DATA[i]->SetLineColor(kBlack);
        h1_DATA[i]->SetMarkerColor(kBlack);
        h1_DATA[i]->SetMarkerSize(1.0);
        h1_DATA[i]->SetMarkerStyle(20);
        h1_DATA[i]->SetStats(0);
        h1_DATA[i]->Draw("E SAME");

        TLegend *l1 = new TLegend(0.7, 0.60, 0.85, 0.85);
        l1->SetFillColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->SetShadowColor(kWhite);
        l1->AddEntry(h1_DATA[i],      " Data",   "lp");
        l1->AddEntry(h1_DY[i],        " DY",     "f");
        l1->AddEntry(h1_TT[i],        " TT",     "f");
        l1->AddEntry(h1_WJets[i],     " WJets",  "f");
        l1->AddEntry(h1_QCD[i],       " QCD",    "f");
        l1->Draw();

        c->cd(i+1);
        TPad *pad2 = new TPad("p_pull", "p_pull", 0.0, 0.0, 1.0, 0.3);
        pad2->Draw();
        pad2->cd();
        pad2->SetTopMargin(0.04);
        pad2->SetRightMargin(0.07);
        pad2->SetBottomMargin(0.4);

        h1_One[i]->SetMinimum(0);
        h1_One[i]->SetMaximum(2);
        h1_One[i]->SetLabelSize(0.1,"XY");
        h1_One[i]->SetTitleSize(0.1,"XY");
        h1_One[i]->SetTitleOffset(1.5);
        h1_One[i]->GetYaxis()->SetNdivisions(3,true);
        h1_One[i]->Draw("HIST");
        h1_Ratio[i]->SetMarkerStyle(20);
        h1_Ratio[i]->SetMarkerSize(0.7);
        h1_Ratio[i]->Draw("SAME E X0");
    }
    
    //
    if(OneFileTest) 
    {
        if(HistName=="mj") HistName="JetMass";
        c->Print( Form("CompareDataMC_%s_pT%i%s%s_OneFileTest.pdf", 
                  HistName.Data(), pTR0p5, doZ?"":"_Zveto", RemoveMuon?"_MuonRemoved":""));
    }
    else 
    { 
        if(HistName=="mj") HistName="JetMass";
        c->Print( Form("CompareDataMC_%s_pT%i%s%s.pdf", 
                  HistName.Data(), pTR0p5, doZ?"":"_Zveto", RemoveMuon?"_MuonRemoved":"")); 
    }
    
    // 
    HistFile->Close();
    delete c;
}


void makeHistFile(bool doZ, bool RemoveMuon) 
{ 
    TString HistDir = "HistFiles";
    if(!doZ) HistDir = HistDir+"_Zveto";
    if(RemoveMuon) HistDir = HistDir+"_MuonRemoved";
   
    cout << "HistDir = " << HistDir << endl;

    vector<TString> DATAFiles;
    vector<TString> DYFiles;
    vector<TString> TTFiles;
    vector<TString> QCDFiles;
    vector<TString> WJetsFiles;
    
    vector<TString> var;
    var.push_back("MJ"); var.push_back("mj");  var.push_back("Nfatjet"); var.push_back("mll"); 
   
    // 
    // 
    // 
    gSystem->Exec("less dataset.txt | grep -v \"#\" > dataset_trimmed.txt");
    string line;
    ifstream dataset("dataset_trimmed.txt");

    string Dataset, Recog, IsData, Lumi, Nfiles; 

    if (dataset.is_open())
    {
        while ( dataset.good() )
        {
            // get a line from input file
            getline (dataset,line);
            stringstream stream(line);
            
            stream >> Dataset;
            stream >> Recog;
            stream >> IsData;
            stream >> Lumi;
            stream >> Nfiles;
            
            if( ! (Recog.find("2012")==string::npos) )  DATAFiles.push_back(Recog); 
            if( ! (Recog.find("DY")==string::npos) )    DYFiles.push_back(Recog); 
            if( ! (Recog.find("TT")==string::npos) )    TTFiles.push_back(Recog); 
            if( ! (Recog.find("W")==string::npos) )     WJetsFiles.push_back(Recog); 
            if( ! (Recog.find("QCD")==string::npos) )   QCDFiles.push_back(Recog); 
        }
    }
    dataset.close();

    if(0) // DEBUG
    {
        cout << "---- DATA " << DATAFiles.size() << endl;
        for(int i=0; i<DATAFiles.size(); i++) cout << DATAFiles.at(i) << endl;
        cout << "---- DY " << DYFiles.size() << endl;
        for(int i=0; i<DYFiles.size(); i++) cout << DYFiles.at(i) << endl;
        cout << "---- TT " << TTFiles.size() << endl;
        for(int i=0; i<TTFiles.size(); i++) cout << TTFiles.at(i) << endl;
        cout << "---- WJets " << WJetsFiles.size() << endl;
        for(int i=0; i<WJetsFiles.size(); i++) cout << WJetsFiles.at(i) << endl;
        cout << "---- QCD " << QCDFiles.size() << endl;
        for(int i=0; i<QCDFiles.size(); i++) cout << QCDFiles.at(i) << endl;
    }


    //
    //
    //
    TH1F *h1_DATA[4][6];
    TH1F *h1_DY[4][6];
    TH1F *h1_TT[4][6];
    TH1F *h1_QCD[4][6];
    TH1F *h1_WJets[4][6];

    for(int ivar=0; ivar<4; ivar++) 
    { 
        for(int ifatjet=0; ifatjet<6; ifatjet++) 
        {
            // DATA
            for(int i=0; i<DATAFiles.size(); i++) 
            {
                
                TString HistName = DATAFiles.at(i); HistName.ReplaceAll("-","_");
                TFile *f = new TFile(Form("%s/%s_pT%i.root", HistDir.Data(), DATAFiles.at(i).Data(), pTR0p5), "READ");
                
                TH1F *h1temp = (TH1F*)f->Get(Form("h1_%s_%s_%ifatjet",HistName.Data(), var.at(ivar).Data(), ifatjet))->Clone();
                
                if(i==0) 
                { 
                    h1_DATA[ivar][ifatjet] = (TH1F*)h1temp->Clone();
                    h1_DATA[ivar][ifatjet]->SetTitle(Form("h1_DATA_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                    h1_DATA[ivar][ifatjet]->SetName(Form("h1_DATA_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                }
                else 
                { 
                    h1_DATA[ivar][ifatjet]->Add(h1temp); 
                } 
            }

            // DY
            for(int i=0; i<DYFiles.size(); i++) 
            {
                
                TString HistName = DYFiles.at(i); HistName.ReplaceAll("-","_");
                TFile *f = new TFile(Form("%s/%s_pT%i.root", HistDir.Data(), DYFiles.at(i).Data(), pTR0p5), "READ");
                
                TH1F *h1temp = (TH1F*)f->Get(Form("h1_%s_%s_%ifatjet",HistName.Data(), var.at(ivar).Data(), ifatjet))->Clone();
                
                if(i==0) 
                { 
                    h1_DY[ivar][ifatjet] = (TH1F*)h1temp->Clone();
                    h1_DY[ivar][ifatjet]->SetTitle(Form("h1_DY_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                    h1_DY[ivar][ifatjet]->SetName(Form("h1_DY_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                }
                else 
                { 
                    h1_DY[ivar][ifatjet]->Add(h1temp); 
                } 
            }

            // TT
            for(int i=0; i<TTFiles.size(); i++) 
            {
                
                TString HistName = TTFiles.at(i); HistName.ReplaceAll("-","_");
                TFile *f = new TFile(Form("%s/%s_pT%i.root", HistDir.Data(), TTFiles.at(i).Data(), pTR0p5), "READ");
                
                TH1F *h1temp = (TH1F*)f->Get(Form("h1_%s_%s_%ifatjet",HistName.Data(), var.at(ivar).Data(), ifatjet))->Clone();
                
                if(i==0) 
                { 
                    h1_TT[ivar][ifatjet] = (TH1F*)h1temp->Clone();
                    h1_TT[ivar][ifatjet]->SetTitle(Form("h1_TT_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                    h1_TT[ivar][ifatjet]->SetName(Form("h1_TT_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                }
                else 
                { 
                    h1_TT[ivar][ifatjet]->Add(h1temp); 
                } 
            }
            
            // WJets
            for(int i=0; i<WJetsFiles.size(); i++) 
            {
                
                TString HistName = WJetsFiles.at(i); HistName.ReplaceAll("-","_");
                TFile *f = new TFile(Form("%s/%s_pT%i.root", HistDir.Data(), WJetsFiles.at(i).Data(), pTR0p5), "READ");
                
                TH1F *h1temp = (TH1F*)f->Get(Form("h1_%s_%s_%ifatjet",HistName.Data(), var.at(ivar).Data(), ifatjet))->Clone();
                
                if(i==0) 
                { 
                    h1_WJets[ivar][ifatjet] = (TH1F*)h1temp->Clone();
                    h1_WJets[ivar][ifatjet]->SetTitle(Form("h1_WJets_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                    h1_WJets[ivar][ifatjet]->SetName(Form("h1_WJets_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                }
                else 
                { 
                    h1_WJets[ivar][ifatjet]->Add(h1temp); 
                } 
            }
            
            // QCD
            for(int i=0; i<QCDFiles.size(); i++) 
            {
                
                TString HistName = QCDFiles.at(i); HistName.ReplaceAll("-","_");
                TFile *f = new TFile(Form("%s/%s_pT%i.root", HistDir.Data(), QCDFiles.at(i).Data(), pTR0p5), "READ");
                
                TH1F *h1temp = (TH1F*)f->Get(Form("h1_%s_%s_%ifatjet",HistName.Data(), var.at(ivar).Data(), ifatjet))->Clone();
                
                if(i==0) 
                { 
                    h1_QCD[ivar][ifatjet] = (TH1F*)h1temp->Clone();
                    h1_QCD[ivar][ifatjet]->SetTitle(Form("h1_QCD_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                    h1_QCD[ivar][ifatjet]->SetName(Form("h1_QCD_%s_%ifatjet", var.at(ivar).Data(), ifatjet));
                }
                else 
                { 
                    h1_QCD[ivar][ifatjet]->Add(h1temp); 
                } 
            }

        }
    }

    //
    //
    //
    TString HistFileName = Form("%s/Hist_pT%i.root", HistDir.Data() ,pTR0p5);
    cout << "[MJ Analysis] Writing " << HistFileName << endl;
    TFile *HistFile = new TFile(HistFileName, "RECREATE");
    gROOT->cd();
    HistFile->cd();
    
    // write histograms
    for(int ivar=0; ivar<4; ivar++) 
    { 
        for(int ifatjet=0; ifatjet<6; ifatjet++) 
        {
            if(ivar<2 && ifatjet==0) continue;
            h1_DATA[ivar][ifatjet]->SetDirectory(0);    h1_DATA[ivar][ifatjet]->Write();
            h1_DY[ivar][ifatjet]->SetDirectory(0);      h1_DY[ivar][ifatjet]->Write();
            h1_TT[ivar][ifatjet]->SetDirectory(0);      h1_TT[ivar][ifatjet]->Write();
            h1_WJets[ivar][ifatjet]->SetDirectory(0);   h1_WJets[ivar][ifatjet]->Write();
            h1_QCD[ivar][ifatjet]->SetDirectory(0);     h1_QCD[ivar][ifatjet]->Write();
        }
    }
    HistFile->Close();

}

    
void DrawDYCR() 
{

    gROOT->ProcessLine(".L ~/macros/JaeStyle.C");
    
    bool doZ=true;
    bool RemoveMuon=true;

    for(int doZ=0; doZ<2; doZ++) 
    { 
        for(int RemoveMuon=0; RemoveMuon<2; RemoveMuon++)  
        {
            // doZ, RemoveMuon
            makeHistFile(doZ, RemoveMuon); 

            //  Draw histograms 
            DrawStack("mll",        doZ, RemoveMuon, pTR0p5);
            DrawStack("MJ",         doZ, RemoveMuon, pTR0p5);
            DrawStack("Nfatjet",    doZ, RemoveMuon, pTR0p5);
            DrawStack("mj",         doZ, RemoveMuon, pTR0p5);
        }
    }

}



