#include <iostream>

#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "THStack.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include "babytree.h"

using namespace std;
bool DoRebin        = 0;
bool OnlyDraw       = 0;
bool OneFileTest    = 0;
bool RemoveMuon     = 1;
int pTR0p5thres     = 30;

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
// per process
//
void DoOneProcess(TChain *ch, int pTR0p5) 
{ 

//    mj
//    MJ
//    FatjetPt 
//    FatjetEta 
//    FatjetPhi

    int             Nfatjet;
    float           MJ;
    vector<float>   *mj;
    vector<float>   *FatjetPt;
    vector<float>   *FatjetEta;
    vector<float>   *FatjetPhi;
    TBranch         *b_MJ;   //!
    TBranch         *b_Nfatjet;   //!
    TBranch         *b_mj;   //!
    TBranch         *b_FatjetPt;   //!
    TBranch         *b_FatjetEta;   //!
    TBranch         *b_FatjetPhi;   //!
    mj        = 0;
    FatjetPt  = 0;
    FatjetEta = 0;
    FatjetPhi = 0;
    ch->SetBranchAddress(Form("MJ_pT%i", pTR0p5), &MJ, &b_MJ);
    ch->SetBranchAddress(Form("Nfatjet_pT%i", pTR0p5), &Nfatjet, &b_Nfatjet);
    ch->SetBranchAddress(Form("mj_pT%i", pTR0p5), &mj, &b_mj);
    ch->SetBranchAddress(Form("FatjetPt_pT%i", pTR0p5), &FatjetPt, &b_FatjetPt);
    ch->SetBranchAddress(Form("FatjetEta_pT%i", pTR0p5), &FatjetEta, &b_FatjetEta);
    ch->SetBranchAddress(Form("FatjetPhi_pT%i", pTR0p5), &FatjetPhi, &b_FatjetPhi);

    TString ChainName = ch->GetTitle();
    cout << "[MJ Analysis] " << ChainName << endl;  

    //
    //
    //
    TH1F *h1_MJ[6], *h1_mll[6], *h1_mj[6], *h1_Nfatjet[6];
    for(int i=0; i<6; i++) 
    {
        h1_MJ[i] = InitTH1F( Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 500);
        h1_mj[i] = InitTH1F( Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 500);
        h1_Nfatjet[i] = InitTH1F( Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  11, -0.5, 10.5);
        h1_mll[i] = InitTH1F( Form("h1_%s_mll_%ifatjet", ch->GetTitle(), i), 
                              Form("h1_%s_mll_%ifatjet", ch->GetTitle(), i), 
                              25, 50, 150);
    }
    // Pile up reweighting hist
    TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_1p317ipb.root");
    TH1F *h1PU = (TH1F*)(fPUFile->Get("puWeights"));
    
    //
    //
    //
    int i_permille_old = 0; 
    TDatime DTStart;
    int StartDate = DTStart.GetDate();
    int StartTime = DTStart.GetTime();
    cout << "[MJ Analysis] Start time : " << (StartTime/10000)%100 << ":"
        << (StartTime/100)%100 << ":" << StartTime%100
        << endl;
   
    //
    //
    //
    InitTree(ch);
    Int_t nentries = (Int_t)ch->GetEntries();
    for(int i = 0; i<nentries; i++)
    {
        ch->GetEntry(i); 

        // Progress indicator begin --------------------------------
        int i_permille = (int)floor(1000 * i / float(nentries));
        TDatime DTCurrent;
        int CurrentDate = DTCurrent.GetDate();
        int CurrentTime = DTCurrent.GetTime();
        int TimeLaps = (CurrentDate-StartDate)*1000000+(CurrentTime-StartTime);
        int TimeToRun = (int)((float)nentries/(float)i)*TimeLaps;
        if (i_permille != i_permille_old)
        {
            // xterm magic from L. Vacavant and A. Cerri
            if (isatty(1))
            {
                printf("\015\033[32m Processed :: \033[1m\033[31m%4.1f %%"
                        "\033[0m\033[32m   Expected processing time :: \033[1m\033[31m%i:%i:%i \033[0m\015",
                        i_permille/10., (TimeToRun/10000)%100<60 ? (TimeToRun/10000)%100 : (TimeToRun/10000)%100-40,
                        (TimeToRun/100)%100<60 ? (TimeToRun/100)%100 : (TimeToRun/100)%100-40,
                        (TimeToRun%100)<60 ? (TimeToRun)%100 : (TimeToRun)%100-40 );
                fflush(stdout);
            }
            i_permille_old = i_permille;
        }
        // Progress indicator end ----------------------------------
      
        bool cut_others = false;
        bool cut_mll    = false;

        // 
        // weights 
        // 
        // Pileup 
        if(!ChainName.Contains("DATA")) 
        {
            EventWeight = EventWeight*nPUScaleFactor2012(h1PU, Npu); // PU   
        }
        // DY reweighting
        if(ChainName.Contains("DY")) 
        { 
            EventWeight = EventWeight * 1177.6*3/2950.;
        }
        // For one file test // FIXME 
        if(OneFileTest && !ChainName.Contains("DATA")) 
        {
            EventWeight = EventWeight * 0.00352; 
        }

        //
        // cuts applied  
        //

        // select only di-muon events
        int Nmuon=0;
        int imus1=-1;
        int imus2=-1;
        for(int ilep=0; ilep<(int)RA4leptonPt->size(); ilep++)
        {
            if(TMath::Abs(RA4leptonId->at(ilep))==13) 
            {
                Nmuon++;
                if(Nmuon==1) imus1=ilep;
                if(Nmuon==2) imus2=ilep;
            }
        } 
        if(RA4leptonId->size()!=2) continue;

        if( Nmuon==2 &&                                        // exactly two muons 
            RA4leptonId->at(imus1)*RA4leptonId->at(imus2)<0    // opposite sign
            //(Npv<15 || Npv>20) 
            ) 
        {   
            cut_mll = true;
            if( TMath::Abs(mll-91.2)<15)                                   // Z mass 
            { 
                cut_others = true;
            }
        }
    
        // Nfatjet counting with threshold 
        int Nfatjet_thres = 0;
        for(int ifatjet=0; ifatjet<(int)FatjetPt->size(); ifatjet++)
        {   
            float FatjetPt_cor = FatjetPt->at(ifatjet); 
            if(RemoveMuon) 
            {
                for(int ilep=0; ilep<(int)RA4leptonPt->size(); ilep++)  
                {
                    float dEta = FatjetEta->at(ifatjet) - RA4leptonEta->at(ilep);
                    float dPhi = FatjetPhi->at(ifatjet) - RA4leptonPhi->at(ilep);
                    if(TMath::Sqrt(dEta*dEta+dPhi*dPhi)<1.0) FatjetPt_cor = FatjetPt_cor - RA4leptonPt->at(ilep);
                }
            } 
            if(FatjetPt_cor>50) Nfatjet_thres++;
        }

        //
        // Fill histograms 
        //
        if(cut_mll) {
            
            if(Nfatjet_thres>3) 
            {
                h1_mll[4]->Fill(mll, EventWeight);
            }
            else if(Nfatjet_thres>2)
            {
                h1_mll[3]->Fill(mll, EventWeight);
            }
            else if(Nfatjet_thres>1)
            {
                h1_mll[2]->Fill(mll, EventWeight);
            }
            else if(Nfatjet_thres>0)
            {
                h1_mll[1]->Fill(mll, EventWeight);
            }
            else if(Nfatjet_thres==0)
            {
                h1_mll[0]->Fill(mll, EventWeight);
            }

            h1_mll[5]->Fill(mll, EventWeight);
        }

        if(cut_others) {
            if(Nfatjet_thres>3) 
            {
                h1_MJ[4]->Fill(MJ, EventWeight);
                for(int imj=0; imj<(int)mj->size(); imj++) h1_mj[4]->Fill(mj->at(imj), EventWeight);
                h1_Nfatjet[4]->Fill(Nfatjet_thres, EventWeight);
            }
            else if(Nfatjet_thres>2)
            {
                h1_MJ[3]->Fill(MJ, EventWeight);
                for(int imj=0; imj<(int)mj->size(); imj++) h1_mj[3]->Fill(mj->at(imj), EventWeight);
                h1_Nfatjet[3]->Fill(Nfatjet_thres, EventWeight);

            }
            else if(Nfatjet_thres>1)
            {
                h1_MJ[2]->Fill(MJ, EventWeight);
                for(int imj=0; imj<(int)mj->size(); imj++) h1_mj[2]->Fill(mj->at(imj), EventWeight);
                h1_Nfatjet[2]->Fill(Nfatjet_thres, EventWeight);

            }
            else if(Nfatjet_thres>0)
            {
                h1_MJ[1]->Fill(MJ, EventWeight);
                for(int imj=0; imj<(int)mj->size(); imj++) h1_mj[1]->Fill(mj->at(imj), EventWeight);
                h1_Nfatjet[1]->Fill(Nfatjet_thres, EventWeight);
            }
            else if(Nfatjet_thres==0)
            {
                h1_MJ[0]->Fill(MJ, EventWeight);
                for(int imj=0; imj<(int)mj->size(); imj++) h1_mj[0]->Fill(mj->at(imj), EventWeight);
                h1_Nfatjet[0]->Fill(Nfatjet_thres, EventWeight);
            }

            h1_MJ[5]->Fill(MJ, EventWeight);
            for(int imj=0; imj<(int)mj->size(); imj++) h1_mj[5]->Fill(mj->at(imj), EventWeight);
            h1_Nfatjet[5]->Fill(Nfatjet_thres, EventWeight);
        }
    }
    
    TString HistFileName = ch->GetTitle();
    HistFileName = Form("HistFiles/%s_pT%i.root", HistFileName.Data(), pTR0p5thres);
    cout << "[MJ Analysis] Writing " << HistFileName << endl;
    TFile *HistFile = new TFile(HistFileName, "RECREATE");
    gROOT->cd();
    HistFile->cd();
    // write histograms
    for(int i=0; i<6; i++)  
    {
        if(i) h1_MJ[i]->SetDirectory(0);       h1_MJ[i]->Write();
        if(i) h1_mj[i]->SetDirectory(0);       h1_mj[i]->Write();
        h1_mll[i]->SetDirectory(0);      h1_mll[i]->Write();
        h1_Nfatjet[i]->SetDirectory(0);  h1_Nfatjet[i]->Write();
    }
    HistFile->Close();

}

//
//
//
void DrawStack(TString HistName, int pTR0p5) 
{ 
    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_pT%i.root",pTR0p5));
   
    char *var; 
    bool skip0fatjet = false;
    if(HistName=="MJ")      { var=(char*)"M_{J} [GeV]"; skip0fatjet=true; }
    if(HistName=="mj")      { var=(char*)"m_{j} [GeV]"; skip0fatjet=true; }
    if(HistName=="mll")     var=(char*)"m_{ll} [GeV]";
    if(HistName=="Nfatjet") var=(char*)"N_{fatjet}";

    TH1F *h1_DATA[6], *h1_DY[6], *h1_TT[6], *h1_WJets[6], *h1_QCD[6], *h1_MC[6], *h1_Ratio[6]; 
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

        h1_Ratio[i]->SetMinimum(0);
        h1_Ratio[i]->SetMaximum(2);
        h1_Ratio[i]->SetMarkerStyle(20);
        h1_Ratio[i]->SetMarkerSize(0.5);
        h1_Ratio[i]->SetLabelSize(0.1,"XY");
        h1_Ratio[i]->SetTitleSize(0.1,"XY");
        h1_Ratio[i]->SetTitleOffset(1.5);
        h1_Ratio[i]->Draw("E X0");
    }
    
    //
    if(OneFileTest) 
    {
        if(HistName=="mj") HistName="JetMass";
        c->Print( Form("CompareDataMC_%s_pT%i%s_OneFileTest.pdf", 
                  HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":""));
    }
    else 
    { 
        if(HistName=="mj") HistName="JetMass";
        c->Print( Form("CompareDataMC_%s_pT%i%s.pdf", 
                  HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"")); 
    }
    
    // 
    HistFile->Close();
    delete c;
}

//
// main 
//
void DrawDYCR() 
{

    gROOT->ProcessLine(".L ~/macros/JaeStyle.C");
    
    // ----------------------------------------
    //  Define chains  
    // ----------------------------------------
    TChain *ch_data  = new TChain("tree", "DATA");
    TChain *ch_dyll  = new TChain("tree", "DY");
    TChain *ch_ttbar = new TChain("tree", "TT");
    TChain *ch_wjets = new TChain("tree", "WJets");
    TChain *ch_qcd   = new TChain("tree", "QCD");
   
    TString BabyDir = "../../baby/";  
    if(OneFileTest) BabyDir = "../../baby_onefile/"; 
    ch_data->Add(BabyDir+"baby_DoubleMu_Run2012D-PromptReco2.root");
    //ch_data->Add(BabyDir+"baby_DoubleMu_Run2012D-PromptReco2DavidSkim.root");
    ch_dyll->Add(BabyDir+"baby_DYJetsToLL_M-50.root");
    ch_ttbar->Add(BabyDir+"baby_TTJets_FullLeptMGDecays.root");
    ch_ttbar->Add(BabyDir+"baby_TTJets_HadronicMGDecays.root");
    ch_ttbar->Add(BabyDir+"baby_TTJets_SemiLeptMGDecays.root");
    ch_wjets->Add(BabyDir+"baby_W2JetsToLNu.root");
    ch_wjets->Add(BabyDir+"baby_W3JetsToLNu.root");
    ch_wjets->Add(BabyDir+"baby_W4JetsToLNu.root");
    ch_qcd->Add(BabyDir+"baby_QCD_HT-1000ToInf.root");
    ch_qcd->Add(BabyDir+"baby_QCD_HT-500To1000.root");
    ch_qcd->Add(BabyDir+"baby_QCD_HT-250To500.root");
    ch_qcd->Add(BabyDir+"baby_QCD_HT-100To250.root");
    
    // ----------------------------------------
    //  Get number of entries 
    // ----------------------------------------
    cout << "data  : " << ch_data->GetEntries() << endl;
    cout << "dyll  : " << ch_dyll->GetEntries() << endl;
    cout << "ttbar : " << ch_ttbar->GetEntries() << endl;
    cout << "wjets : " << ch_wjets->GetEntries() << endl;
    cout << "qcd   : " << ch_qcd->GetEntries() << endl;
    
    if(!OnlyDraw) 
    {
        // ----------------------------------------
        //  Fill histrograms 
        // ----------------------------------------
        DoOneProcess(ch_data,	pTR0p5thres); 
        DoOneProcess(ch_dyll,	pTR0p5thres); 
        DoOneProcess(ch_ttbar,	pTR0p5thres); 
        DoOneProcess(ch_wjets,	pTR0p5thres); 
        DoOneProcess(ch_qcd,	pTR0p5thres); 

        // ----------------------------------------
        //  Make the final histogram file
        // ----------------------------------------
        cout << "[MJ Analysis] Merging result files" << endl;
        gSystem->Exec(Form("rm HistFiles/Hist_pT%i.root",pTR0p5thres));
        gSystem->Exec(Form("hadd -f HistFiles/Hist_pT%i.root HistFiles/*_pT%i.root",pTR0p5thres,pTR0p5thres));
    }

    // ----------------------------------------
    //  Draw histograms 
    // ---------------------------------------- 
    DrawStack("mll",        pTR0p5thres);
    DrawStack("MJ",         pTR0p5thres);
    DrawStack("Nfatjet",    pTR0p5thres);
    DrawStack("mj",         pTR0p5thres);
      

}

