
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
#include "TPaletteAxis.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TLorentzVector.h"
#include "TInterpreter.h"

#include "../babytree_manuel.h"

using namespace std;

void myText(Double_t x,Double_t y, const char *text, Color_t red,float tsize) {

    //Double_t tsize=0.05;
    TLatex l; //l.SetTextAlign(12);
    l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextColor(red);
    l.DrawLatex(x,y,text);

}

float getDEta(float eta1, float eta2)
{
    return abs(eta1-eta2);
}
float getDPhi(float phi1, float phi2)
{
    float absdphi = abs(phi1-phi2);
    if(absdphi < TMath::Pi()) return absdphi;
    else return (2*TMath::Pi() - absdphi);
}
float getDR(float dphi, float deta)
{
    return TMath::Sqrt(dphi*dphi+deta*deta);
}
float getDR(float eta1, float eta2, float phi1, float phi2)
{
    return getDR(getDPhi(phi1, phi2), eta1-eta2);
}

void EventDisplayForBaby(bool sig=false, bool truth=true, bool fatjet=true, int event=-1)
{

    gInterpreter->ExecuteMacro("~/macros/rootlogon.C");
    
    // chain      
    TChain *ch        = new TChain("tree");
    if(!sig) ch->Add("babies/small_quick_cfA_746p1_nleps1_trig0ON.root");
    if(sig) ch->Add("~/scratch/plots/1d_2015d_13Jan2016/2d/48ipb/*.root");
    
    InitBaby(ch); 
   
    Int_t nentries = (Int_t)ch->GetEntries();
    for(int i = 0; i<nentries; i++)
    {
        ch->GetEntry(i); 

        // apply selections (event==-1 && mj_>xxx && ht_>yyy && ...)
        // or 
        // select an event (event!=-1 && event_=xxxxxx)
        // "event" is one of the arguments  
        
        if(event>0 && event_!=event) continue; 

        // 
        //  Event Display
        // 
        float start=0.62;
        float nextline = start;
        float increment = 0.04;
        float offset=0.02;
        int npert;
        int LSPcol[2] = {kCyan+1,kOrange-4};
        if(!sig) npert = 5;
        else npert=3;
        float xalign = 0.82;
        vector<TMarker> constituents, genpart;
        TEllipse *cone[fjets_pt_->size()]; 
        TH2F *h_fatjets = new TH2F("h_fatjets","h_fatjets", 230, -5.0, 5.0, 144, -3.141592, 3.141592); 
        cout << event_ << endl;
        for(int ifj = 0; ifj< (int)fjets_pt_->size(); ifj++)
        {   if(ifj==0) cout << "... Fatjets info(pT, eta, phi, mj)" << endl;
            h_fatjets->Fill(fjets_eta_->at(ifj),fjets_phi_->at(ifj), fjets_m_->at(ifj));
            cout << fjets_pt_->at(ifj) << " " << fjets_eta_->at(ifj) << " " 
                 <<  fjets_phi_->at(ifj) << " " << fjets_m_->at(ifj) << endl;
        }

        TString cname = Form("Event_%lli_",event_);
        TCanvas *c = new TCanvas(cname,cname,1640,760);
        gPad->SetRightMargin(0.35);
        h_fatjets->SetZTitle("m_{j} [GeV]");
        h_fatjets->GetZaxis()->SetTitleSize(0.05);
        h_fatjets->GetZaxis()->SetTitleOffset(0.55);
        h_fatjets->GetYaxis()->SetTitleOffset(0.6);
        h_fatjets->Draw("colz");
        h_fatjets->GetZaxis()->SetRangeUser(0,1.1*h_fatjets->GetMaximum());
        h_fatjets->Draw("colz");
        gPad->Update();

        TPaletteAxis *palette = (TPaletteAxis*)h_fatjets->GetListOfFunctions()->FindObject("palette");
        for(int ifj = 0; ifj< (int)fjets_pt_->size(); ifj++)
        {
            cone[ifj] = new TEllipse(fjets_eta_->at(ifj),fjets_phi_->at(ifj), 1.2, 1.2);
            //  cone[ifj]->SetFillStyle(3003);
            cone[ifj]->SetFillStyle(0);
            Int_t binx,biny,binz;
            h_fatjets->GetBinXYZ(h_fatjets->FindBin(fjets_eta_->at(ifj),fjets_phi_->at(ifj)),binx,biny,binz);
            //cout<<"bin x bin y "<<binx<<" "<<biny<<endl;
            //cout<<"content "<<h_fatjets->GetBinContent(binx,biny)<<endl;

            //cout<<palette<<endl;
            Int_t ci = palette->GetBinColor(binx,biny);
            if(ifj==0) ci = palette->GetBinColor(binx,biny);
            //cout<<"color"<<ci<<endl;
            cone[ifj]->SetFillColor(kGray);
            cone[ifj]->SetLineColor(ci);
        }


        float lepPt, lepEta, lepPhi; 
        if(nels_==1)  
        {
            for(int i=0; i<(int)els_pt_->size(); i++) 
            {
                if(els_miniso_->at(i)>0.1) continue;
                if(els_sigid_->at(i)!=1) continue;
                if(els_pt_->at(i)<20) continue; 
                lepPt   =els_pt_->at(i);
                lepEta  =els_eta_->at(i);
                lepPhi  =els_phi_->at(i); 
            }
        }
        if(nmus_==1)  
        {
            for(int i=0; i<(int)mus_pt_->size(); i++) 
            {
                if(mus_miniso_->at(i)>0.2) continue;
                if(mus_sigid_->at(i)!=1) continue;
                if(mus_pt_->at(i)<20) continue;
                lepPt   =mus_pt_->at(i);
                lepEta  =mus_eta_->at(i);
                lepPhi  =mus_phi_->at(i); 
            }
        }



        myText(xalign+0.01,0.96,Form("H_{T} = %.0f GeV",ht_),1,0.04);
        myText(xalign+0.01,0.91,Form("M_{J} = %.0f GeV",mj_),1,0.04);
        myText(xalign+0.01,0.86,Form("#slash{E}_{T} = %.0f GeV",met_),1,0.04);
        myText(xalign+0.01,0.81,Form("m_{T} = %.0f GeV",mt_),1,0.04);
        myText(xalign+0.01,0.76,Form("reco %s p_{T} = %.0f GeV", nmus_==1?"#mu":"e", lepPt), kBlack, 0.04); 
        nextline=0.76;
        TMarker recolep = TMarker(lepEta, lepPhi, 27);
        recolep.SetMarkerSize(4);
        recolep.SetMarkerColor(kRed);
        genpart.push_back(recolep);
        TMarker mumark = TMarker(xalign,nextline,27);
        mumark.SetNDC();
        mumark.SetX(xalign);
        mumark.SetY(nextline+0.01);
        mumark.SetMarkerSize(2);
        mumark.SetMarkerColor(kRed);
        genpart.push_back(mumark);
        myText(xalign+0.01,0.71,Form("-- large-R jets --"),1,0.04); 
        myText(xalign+0.01,0.67,Form("(pT, eta, phi, mass)"),1,0.03); 
        nextline=0.635;
        for(int ifj = 0; ifj< (int)fjets_pt_->size(); ifj++)
        {
            myText(xalign+0.01,nextline,Form("%3.0f, %2.1f, %2.1f, %3.0f",
                   fjets_pt_->at(ifj),fjets_eta_->at(ifj),fjets_phi_->at(ifj),fjets_m_->at(ifj)),1,0.03); 
            nextline=nextline-0.035;
        } 
        nextline=nextline-0.015;
        myText(xalign+0.01,nextline,Form("-- AK4 jets --"),1,0.04); 
        nextline=nextline-0.04;
        myText(xalign+0.01,nextline,Form("(pT, eta, phi)"),1,0.03); 
        nextline=nextline-0.035;
        for(int ij = 0; ij< (int)jets_pt_->size(); ij++)
        {            
            if(jets_islep_->at(ij)==1) continue;
            myText(xalign+0.01,nextline,Form("%3.0f, %2.1f, %2.1f",
                   jets_pt_->at(ij),jets_eta_->at(ij),jets_phi_->at(ij)),1,0.03); 
            nextline=nextline-0.035;
        }
        
        if(truth) myText(xalign,0.68,"Gen  p_{T} [GeV]",1,0.04);



        
        if(truth) myText(xalign,0.68,"Gen  p_{T} [GeV]",1,0.04);



        TArrow line1 = TArrow();
        if(truth)line1.DrawLineNDC(xalign-0.02,start+0.035,xalign+0.14,start+0.035);
        TMarker met = TMarker(0, met_phi_, 27);
        met.SetMarkerSize(4);
        met.SetMarkerColor(kMagenta);
        genpart.push_back(met);
        TMarker metmark = TMarker(0.81,0.86,27);
        metmark.SetNDC();
        metmark.SetX(xalign);
        metmark.SetY(0.87);
        metmark.SetMarkerSize(2);
        metmark.SetMarkerColor(kMagenta);
        genpart.push_back(metmark);

        TMarker jetmark = TMarker(0.5,0.5,20);
        jetmark.SetNDC();
        jetmark.SetMarkerSize(1.2);
        jetmark.SetX(0.12);
        jetmark.SetY(0.2);
        genpart.push_back(jetmark);
        myText(0.13,0.19,"AK4 Jets",kBlack,0.04);
        TMarker bjetmark = TMarker(0.5,0.5,20);
        bjetmark.SetNDC();
        bjetmark.SetMarkerSize(1.2);
        bjetmark.SetMarkerColor(8);
        bjetmark.SetX(0.12);
        bjetmark.SetY(0.16);
        genpart.push_back(bjetmark);
        myText(0.13,0.15,"CSVM AK4",kBlack,0.04);
        myText(0.08,0.94,"Ring color indicates FJ mass",kBlack,0.04);

        // if(sig) nextline = start - (4*npert+2)*increment-0.02-offset;
        //else nextline = start - 2*npert*increment-offset;


        for(int ij=0;ij<(int)jets_pt_->size();ij++)
        {
            if(jets_islep_->at(ij)==1) continue;
            if(ij==0) cout << "... skinny jet info (pT eta phi)" << endl;
            cout << jets_pt_->at(ij) << " " << jets_eta_->at(ij) << " " <<  jets_phi_->at(ij) << endl;
                TMarker jet = TMarker(jets_eta_->at(ij),jets_phi_->at(ij),20);
                jet.SetMarkerSize(1.2);
                if(jets_csv_->at(ij)>0.890) jet.SetMarkerColor(8);
                constituents.push_back(jet);
        }

/*
        bool drawn=false;
        for(int ifj = 0; ifj< (int)fjets.size(); ifj++){
            vector<fastjet::PseudoJet> cons = fjets[ifj].constituents();
            for(int ics = 0; ics<(int)cons.size();ics++){
                for(int ifj2=0; ifj2<(int)fjets.size(); ifj2++){
                    if(ifj2==ifj) continue;
                    TArrow first = TArrow(cons[ics].eta(),cons[ics].phi_std(), fjets[ifj].eta(),fjets[ifj].phi_std(),0.04,">");
                    if(deltaR(cons[ics].eta(),cons[ics].phi_std(), fjets[ifj2].eta(),fjets[ifj2].phi_std()) < 1.25){
                        //cout<<"type 1 line: con eta phi FJ eta phi"<<cons[ics].eta()<<" "<<cons[ics].phi_std()<<" "<<fjets[ifj2].eta()<<" "<<fjets[ifj2].phi_std()<<endl;
                        first.DrawArrow(cons[ics].eta(),cons[ics].phi_std(), fjets[ifj].eta(),fjets[ifj].phi_std(),0.0045,"|>");
                        drawn=true;
                        break;
                    }
                    else if(deltaR(cons[ics].eta(),cons[ics].phi_std(), fjets[ifj].eta(),fjets[ifj].phi_std()) > 1.15){
                        //cout<<"type 2 line: con eta phi FJ eta phi"<<cons[ics].eta()<<" "<<cons[ics].phi_std()<<" "<<fjets[ifj].eta()<<" "<<fjets[ifj].phi_std()<<endl;
                        first.DrawArrow(cons[ics].eta(),cons[ics].phi_std(), fjets[ifj].eta(),fjets[ifj].phi_std(),0.0045,"|>");
                        drawn=true;
                        break;}
                }
            }
            cons.clear();
        }
        if(drawn) myText(0.12,0.04,"Arrows indicate ambiguous clustering ownership",kBlack,0.03);

*/

        if(truth)
        {
            int col[4] = {30,38,44,46};
            if(!sig) col[1]=46;
            int LSP,top,b,W,Wda,gl;
            LSP=0;top=0;b=0;W=0;Wda=0;gl=0;
            TLorentzVector vgl[2];
            TLorentzVector vt[4];
            TLorentzVector vLSP[2];

            for(unsigned int imc = 0; imc < mc_id_->size(); imc++)
            {
                int id= (int)mc_id_->at(imc);
                int absid = abs(id);
                int mid = (int)mc_mom_->at(imc);
                int absmid= abs(mid);
                if(id!=mid)
                {
                    int marknum=0;
                    TString partname;
                    int color=0;
                    bool flag=false;
                    if(id==1000021){
                        if(gl<2) vgl[gl].SetPtEtaPhiM(mc_pt_->at(imc), mc_eta_->at(imc), mc_phi_->at(imc), 1500);
                        gl++;
                    }  

                    if(id==1000022)
                    {
                        color = LSPcol[LSP];
                        partname = "#tilde{#chi}^{0}_{1}";
                        if(LSP==0) nextline=start-2*npert*increment;
                        else nextline = start - (4*npert+1)*increment-offset;
                        LSP++;
                        flag=true;
                        marknum=25;//27;

                    }
                    if(absid == 6)
                    {
                        //top
                        marknum=22;
                        partname = "t";
                        flag =true;
                        nextline = start - npert*top*increment;
                        color = col[top];
                        if(!sig) {}
                        if(sig) {if(top>1){ nextline-=(increment+offset);}}
                        if(top<4)vt[top].SetPtEtaPhiM(mc_pt_->at(imc), mc_eta_->at(imc), mc_phi_->at(imc), 172.5); 
                        top++;
                    }

                    if(absid==5  && (mid==6 || mid==(-6)))
                    {
                        //b from top
                        marknum=23;
                        partname="b";
                        flag=true;
                        color = col[b];
                        nextline = start - npert*b*increment - 2*increment;
                        // if(!sig){}
                        if(sig) {if(b>1){ nextline-=(increment+offset);}}
                        b++;

                    }
                    if(absid==24  && (mid==6 || mid==(-6)))
                    {
                        //W
                        marknum=34;
                        partname = "W";
                        flag=true;
                        color = col[W];
                        nextline = start - npert*W*increment - increment;
                        // if(!sig){}
                        if(sig){if(W>1){ nextline-=(increment+offset);}}
                        W++;
                    }


                    if((absid==13 || absid== 11 || absid==15) && absmid==24)
                    {
                        marknum=29;
                        if(absid==13) partname = "#mu";
                        else if(absid==11) partname = "e";
                        else partname = "#tau";
                        if(!sig) flag =true;
                        color = col[TMath::FloorNint(Wda/2)];
                        nextline = start - npert*increment*TMath::FloorNint(Wda/2)-3*increment;
                        if(TMath::FloorNint(Wda/2) != Wda/2) nextline-=increment;
                        Wda++;

                    }

                    if((absid==12 || absid== 14 || absid==16) && absmid==24)
                    {
                        marknum=30;
                        partname = "#nu";     
                        if(!sig) flag =true;
                        color = col[TMath::FloorNint(Wda/2)];
                        nextline = start - npert*increment*TMath::FloorNint(Wda/2)-4*increment;
                        //if(TMath::FloorNint(Wda/2) != Wda/2) nextline-=increment;
                        Wda++;

                    }
                    if((absid>0 && absid<5) && absmid==24)
                    {
                        marknum=3;
                        if(absid==1) partname="d";
                        if(absid==2) partname="u";
                        if(absid==3) partname="s";
                        if(absid==4) partname="c";

                        if(!sig) flag =true;
                        color = col[TMath::FloorNint(Wda/2)];
                        nextline = start - npert*increment*TMath::FloorNint(Wda/2)-3*increment;
                        // if(TMath::FloorNint(Wda/2) != Wda/2) nextline-=increment;
                        if(absid % 2 ==1) nextline-=increment; 
                        Wda++;

                    }

                    if(flag)
                    {
                        TMarker temp = TMarker(mc_eta_->at(imc),mc_phi_->at(imc), marknum);
                        TMarker legmark = TMarker(0.81,nextline,marknum);
                        legmark.SetNDC();
                        legmark.SetX(xalign+0.01);
                        legmark.SetY(nextline+0.011);
                        TString label = partname + "  " + Form("%.0f",mc_pt_->at(imc));
                        myText(xalign+0.03,nextline,partname,color,0.04);
                        myText(xalign+0.07,nextline,Form("%.0f",mc_pt_->at(imc)),color,0.04);
                        temp.SetMarkerColor(color);
                        legmark.SetMarkerColor(color);
                        temp.SetMarkerSize(2);
                        legmark.SetMarkerSize(2);

                        if(absid==6 || absid==13 || absid==1000022){ temp.SetMarkerSize(3);} 
                        genpart.push_back(temp);
                        genpart.push_back(legmark);
                        flag=false;
                    }
                }
            }
            if(sig && LSP==0){
                int color = LSPcol[LSP];
                TString partname = "#tilde{#chi}^{0}_{1}";
                if(LSP==0) nextline=start-2*npert*increment;
                else nextline = start - (4*npert+1)*increment-offset;
                LSP++;
                int  marknum=25;//27;
                // TMarker temp = TMarker(mc_doc_eta->at(imc),mc_doc_phi->at(imc), marknum);
                vLSP[0]= vgl[0]-vt[0]-vt[1];
                TMarker temp = TMarker(vLSP[0].Eta(),vLSP[0].Phi(), marknum);
                TMarker legmark = TMarker(0.81,nextline,marknum);
                legmark.SetNDC();
                legmark.SetX(xalign+0.01);
                legmark.SetY(nextline+0.011);
                TString label = partname + "  " + Form("%.0f",vLSP[0].Pt());
                myText(xalign+0.03,nextline,partname,color,0.04);
                myText(xalign+0.07,nextline,Form("%.0f",vLSP[0].Pt()),color,0.04);
                temp.SetMarkerColor(color);
                legmark.SetMarkerColor(color);
                temp.SetMarkerSize(3);
                legmark.SetMarkerSize(2);


                genpart.push_back(temp);
                genpart.push_back(legmark);

            }
            if(sig && LSP==1){
                int color = LSPcol[LSP];
                TString partname = "#tilde{#chi}^{0}_{1}";
                if(LSP==0) nextline=start-2*npert*increment;
                else nextline = start - (4*npert+1)*increment-offset;
                LSP++;
                int marknum=25;//27;
                // TMarker temp = TMarker(mc_doc_eta->at(imc),mc_doc_phi->at(imc), marknum);
                vLSP[1]= vgl[1]-vt[2]-vt[3];
                TMarker temp = TMarker(vLSP[1].Eta(),vLSP[1].Phi(), marknum);
                TMarker legmark = TMarker(0.81,nextline,marknum);
                legmark.SetNDC();
                legmark.SetX(xalign+0.01);
                legmark.SetY(nextline+0.011);
                TString label = partname + "  " + Form("%.0f",vLSP[1].Pt());
                myText(xalign+0.03,nextline,partname,color,0.04);
                myText(xalign+0.07,nextline,Form("%.0f",vLSP[1].Pt()),color,0.04);
                temp.SetMarkerColor(color);
                legmark.SetMarkerColor(color);
                temp.SetMarkerSize(3);
                legmark.SetMarkerSize(2);


                genpart.push_back(temp);
                genpart.push_back(legmark);

            }
        }

        for(int ix=0;ix<(int)genpart.size();ix++)
        {
            genpart[ix].Draw("same");
        }
        for(int iy=0;iy<(int)constituents.size();iy++)
        {
            constituents[iy].Draw("same");
        }

        for(int ifj = 0; ifj< (int)fjets_pt_->size(); ifj++)
        {
            if(fatjet) cone[ifj]->Draw();
        }
        h_fatjets->Draw("colz same");
        h_fatjets->SetTitle(Form("Event %lli",event_));

        h_fatjets->SetStats(0); 
        h_fatjets->SetXTitle("#eta"); 
        h_fatjets->SetYTitle("#phi"); 

        TString savename;
        if(!sig) savename = "fig/EventDisplay_"+cname+Form("%i_FJ",nfjets_); 
        else  savename = "fig/EventDisplay_"+cname+Form("%i_FJ",nfjets_);
        if(truth) savename+="_truth";
        else savename+="_notruth";	  
        if(fatjet) savename+="_fatjet";
        else savename+="_nofatjet";
        c->Print(savename+".pdf");
        //c->Print(savename+".C");

        c->Close();
        h_fatjets->Delete();
    }
}

void EventDisplay(int event=1) 
{
/*
    // R4 
    event = 112119467; 
    EventDisplayForBaby(true, false, true, event);
    event = 1115852829; 
    EventDisplayForBaby(true, false, true, event);
    event = 226892978; 
    EventDisplayForBaby(true, false, true, event);
    event = 183543170; 
    EventDisplayForBaby(true, false, true, event);
*/
    event = 943292398;
    EventDisplayForBaby(true, false, true, event);
}
