
#include <iostream>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TString.h"
#include "TChain.h"
#include "TMarker.h"
using namespace std;

int pTR0p5=30;
int doTT=true;

float getDPhi(float phi1, float phi2)
{
    float absdphi = abs(phi1-phi2);
    if(absdphi < TMath::Pi()) return absdphi;
    else 2*TMath::Pi() - absdphi;
}

void EventDisplayForBaby()
{

    // chain      
    TChain *ch        = new TChain("tree");
    TString BabyDir = "../../../babies/HT500MET100/";  
    if(doTT) ch->Add(BabyDir+"baby_TTJets_SemiLeptMGDecays_8TeV_f*.root");
    if(!doTT) ch->Add(BabyDir+"baby_T1tttt_*25*.root");

    // branch
    int             run;
    int             event;
    int             lumiblock;
    float           EventWeight;
    float           Npu;
    int             Npv;
    int             Nfatjet;
    int             Nskinnyjet;
    int             NBtagCSVM;
    float           MJ;
    float           MET;
    float           HT;
    float           METPhi;
    vector<float>   *mj;
    vector<float>   *FatjetPt;
    vector<float>   *FatjetEta;
    vector<float>   *FatjetPhi;
    vector<float>   *RA4MusPt;
    vector<float>   *RA4MusPhi;
    vector<float>   *RA4MusEta;
    vector<float>   *JetPt;
    vector<float>   *JetEta;
    vector<float>   *JetPhi;
    vector<float>   *JetCSV;
    vector<float>   *RA4MusVetoPt;
    vector<float>   *RA4ElsVetoPt;
    vector<float>   *GenPt;
    vector<float>   *GenPhi;
    vector<float>   *GenEta;
    vector<float>   *GenId;
    vector<float>   *GenMId;
    vector<float>   *GenGMId;

    TBranch         *b_run;   //!
    TBranch         *b_event;   //!
    TBranch         *b_lumiblock;   //!
    TBranch         *b_EventWeight;   //!
    TBranch         *b_Npu;   //!
    TBranch         *b_Npv;   //!
    TBranch         *b_Nfatjet;   //!
    TBranch         *b_Nskinnyjet;   //!
    TBranch         *b_NBtagCSVM;   //!
    TBranch         *b_MJ;   //!
    TBranch         *b_MET;   //!
    TBranch         *b_HT;   //!
    TBranch         *b_METPhi;   //!
    TBranch         *b_mj;   //!
    TBranch         *b_FatjetPt;   //!
    TBranch         *b_FatjetEta;   //!
    TBranch         *b_FatjetPhi;   //!
    TBranch         *b_RA4MusPt;   //!
    TBranch         *b_RA4MusPhi;   //!
    TBranch         *b_RA4MusEta;   //!
    TBranch         *b_JetPt;   //!
    TBranch         *b_JetEta;   //!
    TBranch         *b_JetPhi;   //!
    TBranch         *b_JetCSV;   //!
    TBranch         *b_RA4MusVetoPt;   //!
    TBranch         *b_RA4ElsVetoPt;   //!
    TBranch         *b_GenPt;   //!
    TBranch         *b_GenEta;   //!
    TBranch         *b_GenPhi;   //!
    TBranch         *b_GenId;   //!
    TBranch         *b_GenMId;   //!
    TBranch         *b_GenGMId;   //!
    mj        = 0;
    FatjetPt  = 0;
    FatjetEta = 0;
    FatjetPhi = 0;
    RA4MusPt  = 0;
    RA4MusPhi = 0;
    RA4MusEta = 0;
    JetPt  = 0;
    JetEta = 0;
    JetPhi = 0;
    JetCSV = 0;
    RA4MusVetoPt  = 0;
    RA4ElsVetoPt  = 0;
    GenPt  = 0;
    GenEta = 0;
    GenPhi = 0;
    GenId = 0;
    GenMId = 0;
    GenGMId = 0;

    ch->SetBranchAddress("run", &run, &b_run);
    ch->SetBranchAddress("event", &event, &b_event);
    ch->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
    ch->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
    ch->SetBranchAddress("Npu", &Npu, &b_Npu);
    ch->SetBranchAddress("Npv", &Npv, &b_Npv);
    ch->SetBranchAddress(Form("Nfatjet_pT%i", pTR0p5), &Nfatjet, &b_Nfatjet);
    ch->SetBranchAddress("Nskinnyjet", &Nskinnyjet, &b_Nskinnyjet);
    ch->SetBranchAddress("NBtagCSVM", &NBtagCSVM, &b_NBtagCSVM);
    ch->SetBranchAddress(Form("MJ_pT%i", pTR0p5), &MJ, &b_MJ);
    ch->SetBranchAddress("MET", &MET, &b_MET);
    ch->SetBranchAddress("HT", &HT, &b_HT);
    ch->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
    ch->SetBranchAddress(Form("mj_pT%i", pTR0p5), &mj, &b_mj);
    ch->SetBranchAddress(Form("FatjetPt_pT%i", pTR0p5), &FatjetPt, &b_FatjetPt);
    ch->SetBranchAddress(Form("FatjetEta_pT%i", pTR0p5), &FatjetEta, &b_FatjetEta);
    ch->SetBranchAddress(Form("FatjetPhi_pT%i", pTR0p5), &FatjetPhi, &b_FatjetPhi);
    ch->SetBranchAddress("RA4MusPt", &RA4MusPt, &b_RA4MusPt);
    ch->SetBranchAddress("RA4MusPhi", &RA4MusPhi, &b_RA4MusPhi);
    ch->SetBranchAddress("RA4MusEta", &RA4MusEta, &b_RA4MusEta);
    ch->SetBranchAddress("JetPt",  &JetPt,  &b_JetPt);
    ch->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
    ch->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
    ch->SetBranchAddress("JetCSV", &JetCSV, &b_JetCSV);
    ch->SetBranchAddress("RA4MusVetoPt", &RA4MusVetoPt, &b_RA4MusVetoPt);
    ch->SetBranchAddress("RA4ElsVetoPt", &RA4ElsVetoPt, &b_RA4ElsVetoPt);
    ch->SetBranchAddress("GenPt",  &GenPt,  &b_GenPt);
    ch->SetBranchAddress("GenPhi", &GenPhi, &b_GenPhi);
    ch->SetBranchAddress("GenEta", &GenEta, &b_GenEta);
    ch->SetBranchAddress("GenId", &GenId, &b_GenId);
    ch->SetBranchAddress("GenMId", &GenMId, &b_GenMId);
    ch->SetBranchAddress("GenGMId", &GenGMId, &b_GenGMId);

    
    Int_t nentries = (Int_t)ch->GetEntries();
    //for(int i = 0; i<nentries; i++)
    for(int i = 0; i<1000; i++)
    {
        ch->GetEntry(i); 

        //selection 
        if( RA4MusPt->size()!=1 ) continue;
        if( RA4MusVetoPt->size()>0 ) continue;
        if( RA4ElsVetoPt->size()>0 ) continue;
        
        if( !(HT>500                             && 
             MET>100                             && 
             RA4MusPt->at(0)>20                  &&  
             NBtagCSVM==2) 
            ) continue;
        
        if( Nfatjet!=3) continue; 
        //if( mj->at(0)>40 && mj->at(1)>40) continue; 
        if(getDPhi(RA4MusPhi->at(0),METPhi)>0.5) continue; // FIXME

        //
        TH2F *h2 = new TH2F("h2","h2", 115, -5.0, 5.0, 72, -3.141592, 3.141592);
        //

        // Draw fat jets
        TEllipse *cone[FatjetEta->size()]; 
        for(int ijets=0; ijets<(int)FatjetEta->size(); ijets++){
            h2->Fill(FatjetEta->at(ijets), FatjetPhi->at(ijets), mj->at(ijets));
            cone[ijets] = new TEllipse(FatjetEta->at(ijets), FatjetPhi->at(ijets), 1.2, 1.2);
            cone[ijets]->SetFillStyle(3003);
            cone[ijets]->SetFillColor(kYellow);
    //        cone[ijets]->Draw("P");
        }
       
        // Draw gen particles 
        vector<TMarker> constituents, genpart;
        for(unsigned int imc = 0; imc < GenPt->size(); imc++)
        {
            int id  = GenId->at(imc);
            int mid = GenMId->at(imc);
            int gid = GenGMId->at(imc);
            int marknum=0;
            TString partname;
            if(abs(id) == 1000021){
                marknum=8;
                partname = "gluino";
            }
            if(abs(id) == 1000022){
                marknum=24;
                partname = "LSP";
            }
            if(id == 6 || id==(-6)){
                //top
                marknum=22;
                partname = "top";
            }
            if((id==5 && mid==6) || (id==(-5) && mid==(-6)) ){
                //b from top
                marknum=23;
                partname="b";
            }
            if( (id==24 && mid==6) || (id==(-24) && mid==(-6)) ){
                //W
                marknum=34;
                partname = "W";
            }

            if((id==(-13) && mid==24 && gid==6)||(id==(13) && mid==(-24) && gid==(-6))){
                //mu
                marknum=29;
                partname = "#mu";
            }

            if((id==(14) && mid==24 && gid==6)||(id==(-14) && mid==(-24) && gid==(-6))){
                //nu
                marknum=30;
                partname = "#nu";
            }
            if(abs(id)<5 && abs(mid)==24 && abs(gid)==6) 
            { 
                // Hadronic W
                marknum=3;
                if(abs(id)==1) partname = "d";
                if(abs(id)==2) partname = "u";
                if(abs(id)==3) partname = "s";
                if(abs(id)==4) partname = "c";

            }
            if(marknum!=0){
                TMarker temp = TMarker(GenEta->at(imc),GenPhi->at(imc), marknum);
                if((id>0 && id!=13) || id==(-13)){temp.SetMarkerColor(30); /*leg->AddEntry(&temp,partname,"p");*/} //from t
                else if((id<0 &&id!=(-13)) || id==13){ temp.SetMarkerColor(46); partname = "anti"+partname; /*leg->AddEntry(&temp,partname,"p");*/}//from tbar
                if(abs(id)<5 && abs(mid)==24 && gid==6)  temp.SetMarkerColor(30);
                if(abs(id)<5 && abs(mid)==24 && gid==-6)  temp.SetMarkerColor(46);
                temp.SetMarkerSize(2);
                if(id==6 || id==(-6)) temp.SetMarkerSize(3);
                if(abs(id) > 1000000)  temp.SetMarkerColor(9);
                genpart.push_back(temp);
            }
        }
        
        // Draw MET 
        TMarker met = TMarker(0, METPhi, 27);
        met.SetMarkerSize(3);
        met.SetMarkerColor(kMagenta);
        genpart.push_back(met);


        // jets 
        for(int ij=0;ij<(int)JetPt->size();ij++){
            TMarker jet = TMarker(JetEta->at(ij),JetPhi->at(ij),20);
            if(JetCSV->at(ij)>0.679) jet.SetMarkerColor(38);
            constituents.push_back(jet);
        }
        
        // Draw muons 
        TMarker mu = TMarker(RA4MusEta->at(0), RA4MusPhi->at(0), 27);
        mu.SetMarkerSize(3);
        mu.SetMarkerColor(kRed);

        // Canvas 
        TString cname = Form("Event_%i",event);
        TCanvas *c = new TCanvas(cname,cname,1160,800);
        gPad->SetRightMargin(0.15);
        h2->Draw("colz");
        h2->SetTitle(Form("run=%i lumi=%i event=%i R=%.1f", run, lumiblock, event, 1.2));
        h2->SetMaximum(250); 
        h2->SetStats(0); 
        h2->SetXTitle("#eta"); 
        h2->SetYTitle("#phi"); 
        h2->SetZTitle("m_{j} [GeV]");
        h2->Draw("colz");

        for(int ix=0;ix<(int)genpart.size();ix++)
        {
            genpart[ix].Draw("same");
        }
        for(int iy=0;iy<(int)constituents.size();iy++)
        {
            constituents[iy].Draw("same");
        }
        for(int ifj = 0; ifj< (int)FatjetEta->size(); ifj++)
        {
            cone[ifj]->Draw();
        }
        met.Draw("same");
        mu.Draw("same");
        h2->Draw("colz same");

        //
        c->SaveAs(Form("FiguresEventDisplay/EventDisplay_%s_Run%i_Lumi%i_Event%i_R%.1f_5.pdf", 
                    doTT?"TT_sl":"T1tttt_f1400_25",run, lumiblock, event));
        h2->Reset(); 
        for(int ijets=0; ijets<(int)FatjetEta->size(); ijets++) delete cone[ijets];
   

        // print out 
        cout << "-------------------------------------------------------------------- " << endl;
        cout << Form("Run=%i Lumi=%i Event=%i", run, lumiblock, event) << endl;
        cout << "-- Fat jet " << endl;
        for(int ijets=0; ijets<(int)FatjetEta->size(); ijets++){
            cout << "pT: " << FatjetPt->at(ijets) << " " 
                 << "eta: " << FatjetEta->at(ijets) << " " 
                 << "phi: " << FatjetPhi->at(ijets) << " "
                 << "mj:" << mj->at(ijets) << " "
                 << endl;
        }
        cout << "-- skinny jet " << endl;
        for(int ijets=0; ijets<(int)JetEta->size(); ijets++){
            cout << "pT: " << JetPt->at(ijets) << " " 
                 << "eta: " << JetEta->at(ijets) << " " 
                 << "phi: " << JetPhi->at(ijets) << " "
                 << "csv: " << JetCSV->at(ijets) << " "
                 << endl;
        }
        cout << "-- muon " << endl;
            cout << "pT: " << RA4MusPt->at(0) << " " 
                 << "eta: " << RA4MusEta->at(0) << " " 
                 << "phi: " << RA4MusPhi->at(0) << " "
                 << endl;
        cout << "-- MET " << endl;
            cout << "pT: " << MET << " " 
                 << "phi: " << METPhi << " "
                 << endl;
        cout << "-- Gen " << endl;
        for(unsigned int imc = 0; imc < GenPt->size(); imc++) 
        { 
            int id  = GenId->at(imc);
            if(abs(id) > 25) continue; 
//            if(abs(id) == 21) continue; 
            cout << "Id: "  << id << " " 
                 << "Mother Id: "  << GenMId->at(imc) << " " 
                 << "GMother Id: "  << GenGMId->at(imc) << " " 
                 << "pT: "  << GenPt->at(imc) << " " 
                 << "eta: " << GenEta->at(imc) << " " 
                 << "phi: " << GenPhi->at(imc) << " "
                 << endl;
        }

    }
}
