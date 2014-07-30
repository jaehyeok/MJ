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
void DoOneProcess(TString Recognizer, int pTR0p5, bool doZ=true, bool RemoveMuon=false) 
{ 
    TString ChainName = Recognizer;
    ChainName.ReplaceAll("-","_");
    TChain *ch  = new TChain("tree", ChainName.Data());
    TString BabyDir = "/cms24r0/jaehyeok/baby/";  
    //gSystem->Exec(Form("ls %s/baby_%s*.root", BabyDir.Data(), Recognizer.Data())); // DEBIUG
    ch->Add(BabyDir+"baby_"+Recognizer+"*.root");

    TChain *ch_total = new TChain("tree", ChainName.Data());
    ch_total->Add(BabyDir+"baby_"+Recognizer+"*.root");
    
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
    mj          = 0;
    FatjetPt    = 0;
    FatjetEta   = 0;
    FatjetPhi   = 0;
    ch->SetBranchAddress(Form("MJ_pT%i", pTR0p5), &MJ, &b_MJ);
    ch->SetBranchAddress(Form("Nfatjet_pT%i", pTR0p5), &Nfatjet, &b_Nfatjet);
    ch->SetBranchAddress(Form("mj_pT%i", pTR0p5), &mj, &b_mj);
    ch->SetBranchAddress(Form("FatjetPt_pT%i", pTR0p5), &FatjetPt, &b_FatjetPt);
    ch->SetBranchAddress(Form("FatjetEta_pT%i", pTR0p5), &FatjetEta, &b_FatjetEta);
    ch->SetBranchAddress(Form("FatjetPhi_pT%i", pTR0p5), &FatjetPhi, &b_FatjetPhi);

    cout << "[MJ Analysis] " << Recognizer << endl;  

    //
    //
    //
    TH1F *h1_MJ[6], *h1_mll[6], *h1_mj[6], *h1_Nfatjet[6];
    for(int i=0; i<6; i++) 
    {
        h1_MJ[i] = InitTH1F( Form("h1_%s_MJ_%ifatjet", ChainName.Data(), i), 
                             Form("h1_%s_MJ_%ifatjet", ChainName.Data(), i), 
                             20, 0, 500);
        h1_mj[i] = InitTH1F( Form("h1_%s_mj_%ifatjet", ChainName.Data(), i), 
                             Form("h1_%s_mj_%ifatjet", ChainName.Data(), i), 
                             20, 0, 500);
        h1_Nfatjet[i] = InitTH1F( Form("h1_%s_Nfatjet_%ifatjet", ChainName.Data(), i), 
                                  Form("h1_%s_Nfatjet_%ifatjet", ChainName.Data(), i), 
                                  11, -0.5, 10.5);
        if(doZ) 
        {
            h1_mll[i] = InitTH1F( Form("h1_%s_mll_%ifatjet", ChainName.Data(), i), 
                                  Form("h1_%s_mll_%ifatjet", ChainName.Data(), i), 
                                  25, 50, 150); 
        } 
        else if(!doZ) 
        {
            h1_mll[i] = InitTH1F( Form("h1_%s_mll_%ifatjet", ChainName.Data(), i), 
                                  Form("h1_%s_mll_%ifatjet", ChainName.Data(), i), 
                                  25, 56.2, 56.2+25*15); 
        }

    }
    // Pile up reweighting hist
    TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_19p5ifb.root");
    TH1F *h1PU = (TH1F*)(fPUFile->Get("puWeights"));
    
    //
    //
    //
    //int i_permille_old = 0; 
    TDatime DTStart;
    //int StartDate = DTStart.GetDate();
    int StartTime = DTStart.GetTime();
    cout << "[MJ Analysis] Start time : " << (StartTime/10000)%100 << ":"
        << (StartTime/100)%100 << ":" << StartTime%100
        << endl;
   
    //
    //
    //
    InitTree(ch);
    Int_t nentries = (Int_t)ch->GetEntries();
    cout << "Nentries : " << nentries << endl;
    for(int i = 0; i<nentries; i++)
    {
        ch->GetEntry(i); 
/*
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
*/ 
        if((i%1000000)==0) cout << "\015\033[32m[Progress] " << Recognizer << " :: "  << i << " of " << nentries 
                                << "(" << (int)((float)i/(float)nentries*100) << "%) done\033[0m\015" << endl;
        // Progress indicator end ----------------------------------
       

        bool cut_others = false;
        bool cut_mll    = false;

        // 
        // weights 
        // 
        // Pileup 
        if(!Recognizer.Contains("2012")) 
        {
            EventWeight = EventWeight*nPUScaleFactor2012(h1PU, Npu); // PU   
        }
        // DY reweighting
        if(Recognizer.Contains("DY")) 
        { 
            EventWeight = EventWeight * 1177.6*3/2950.;
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
            if( doZ && TMath::Abs(mll-91.2)<15)                         // Z mass 
            { 
                cut_others = true;
            }
            if( !doZ && (TMath::Abs(mll-91.2)>15 && mll>56.2))          // Z veto 
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
                    if(TMath::Sqrt(dEta*dEta+dPhi*dPhi)<0.9) FatjetPt_cor = FatjetPt_cor - RA4leptonPt->at(ilep);
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
    
    TString HistFileName = Form("HistFiles%s%s/%s_pT%i.root", 
                                 doZ ? "":"_Zveto", 
                                 RemoveMuon ? "_MuonRemoved":"", 
                                 Recognizer.Data(), pTR0p5);
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
