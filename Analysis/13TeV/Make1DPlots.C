#include <iostream>
#include <iomanip> // for setw()

#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"
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
#include "TColor.h"

using namespace std;
//bool DoLog          = 0;
bool doData         = 1;
bool SignalScale    = 1;
bool DrawOnlyAllFJ  = 1;
bool NormToData     = 1;

//
// "OLD(Jack's)" UCSB RA4 Color scheme
//
/*
TColor ucsb_blue(1000, 1/255.,57/255.,166/255.);
TColor ucsb_gold(1001, 255/255.,200/255.,47/255);
TColor penn_red(1002, 149/255.,0/255.,26/255.);
TColor uf_orange(1003, 255/255.,74/255.,0/255.);
TColor uo_green(1004, 0/255.,79/255.,39/255.);
TColor tcu_purple(1005, 52/255.,42/255.,123/255.);
TColor tar_heel_blue(1006, 86/255.,160/255.,211/255.);
TColor sig_teal(1007, 96/255.,159/255.,128/255.);
TColor sig_gold(1008, 215/255.,162/255.,50/255.);
TColor seal_brown(1010, 89/255.,38/255.,11/255.);
*/

TColor light_blue(1011, 153/255.,220/255.,255/255.);
TColor med_blue(1012, 1/255.,148/255.,218/255.);
TColor red(1015, 250/255.,96/255.,1/255.);
TColor skype_green(1018,9/255.,186/255.,1/255.);
// TColor purple(1019, 172/255.,46/255.,135/255.);
TColor purple(1019, 183/255.,66/255.,176/255.);
TColor ucsb_gold(1020, 255/255.,200/255.,47/255);

// New color scheme
enum 
{
    c_tt_1l     = 1012, // ucsb_blue
    c_tt_2l     = 1011, // tar_heel_blue
    c_wjets     = 1018, // ucsb_gold
    c_singlet   = 1015,
    c_zjets     = 1020,
    c_other     = 1019
};

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
    h1->SetMinimum(0.);
}

//
// Stacks
//
void Make1DPlots(TString HistName, char* Selection, int NMergeBins=1, bool DoLog=false, float Lumi=40) 
{ 
    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");

    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_%s.root", Selection));
   
    char *var; 
    if(HistName=="dRlep")                 	var=(char*)"dR(lepton,FJ)";
    if(HistName=="dPhiMET")                	var=(char*)"dPhi(MET,FJ)";
    if(HistName=="dPhiMETlep")             	var=(char*)"dPhi(MET,lep)";
    if(HistName=="dRbmin")                	var=(char*)"min dR(bjet,FJ)";
    if(HistName=="mindPhibb")            	var=(char*)"minimum #Delta#phi(b,b)";
    if(HistName=="MET")                 	var=(char*)"MET [GeV]";
    if(HistName=="METPhi")              	var=(char*)"#phi(MET)";
    if(HistName=="METx")              	    var=(char*)"METx [GeV]";
    if(HistName=="METy")              	    var=(char*)"METy [GeV]";
    if(HistName=="FatjetPt")            	var=(char*)"pT(Fatjet) [GeV]";                   
    if(HistName=="FatjetPt1")            	var=(char*)"pT(Fatjet1) [GeV]";                   
    if(HistName=="FatjetPt2")            	var=(char*)"pT(Fatjet2) [GeV]";                   
    if(HistName=="FatjetPt3")            	var=(char*)"pT(Fatjet3) [GeV]";                   
    if(HistName=="FatjetPt4")            	var=(char*)"pT(Fatjet4) [GeV]";                   
    if(HistName=="mj1")            	        var=(char*)"mj1 [GeV]";                   
    if(HistName=="mj2")            	        var=(char*)"mj2 [GeV]";                   
    if(HistName=="mj3")            	        var=(char*)"mj3 [GeV]";                   
    if(HistName=="mj4")            	        var=(char*)"mj4 [GeV]";                   
    if(HistName=="mj08_1")            	    var=(char*)"m(J_{1}) (R=0.8) [GeV]";                   
    if(HistName=="mj08_2")            	    var=(char*)"m(J_{2}) (R=0.8) [GeV]";                   
    if(HistName=="mj08_3")            	    var=(char*)"m(J_{3}) (R=0.8) [GeV]";                   
    if(HistName=="mj08_4")            	    var=(char*)"m(J_{4}) (R=0.8) [GeV]";                   
    if(HistName=="mj1OverMJ")            	var=(char*)"mj1/MJ";                   
    if(HistName=="mj2OverMJ")            	var=(char*)"mj2/MJ";                   
    if(HistName=="mj3OverMJ")            	var=(char*)"mj3/MJ";                   
    if(HistName=="mj4OverMJ")            	var=(char*)"mj4/MJ";                   
    if(HistName=="N1")            	        var=(char*)"N(constituents)1";                   
    if(HistName=="N2")            	        var=(char*)"N(constituents)2";                   
    if(HistName=="N3")            	        var=(char*)"N(constituents)3";                   
    if(HistName=="N4")            	        var=(char*)"N(constituents)4";                   
    if(HistName=="mjOverPt1")            	var=(char*)"mj1/pT1";                   
    if(HistName=="mjOverPt2")            	var=(char*)"mj2/pT2";                   
    if(HistName=="mjOverPt3")            	var=(char*)"mj3/pT3";                   
    if(HistName=="mjOverPt4")            	var=(char*)"mj4/pT4";                   
    if(HistName=="mj3overmj2")     	        var=(char*)"mj3/mj2 [GeV]";                   
    if(HistName=="mj2overmj1")     	        var=(char*)"mj2/mj1 [GeV]";                   
    if(HistName=="FatjetPhi1")            	var=(char*)"#phi(Fatjet1)";                   
    if(HistName=="FatjetPhi2")            	var=(char*)"#phi(Fatjet2)";                   
    if(HistName=="FatjetPhi3")            	var=(char*)"#phi(Fatjet3)";                   
    if(HistName=="FatjetPhi4")            	var=(char*)"#phi(Fatjet4)";                   
    if(HistName=="FatjetEta1")            	var=(char*)"#eta(Fatjet1)";                   
    if(HistName=="FatjetEta2")            	var=(char*)"#eta(Fatjet2)";                   
    if(HistName=="FatjetEta3")            	var=(char*)"#eta(Fatjet3)";                   
    if(HistName=="FatjetEta4")            	var=(char*)"#eta(Fatjet4)";                   
    if(HistName=="FatjetEta")            	var=(char*)"#eta(Fatjet)";                   
    if(HistName=="DPhi")                	var=(char*)"#Delta(MET, muon)";                  
    if(HistName=="dRFJ")                	var=(char*)"min #Delta R(FJ,FJ)";                  
    if(HistName=="dPhiFJ")                	var=(char*)"min #Delta #phi(FJ,FJ)";                  
    if(HistName=="dEtaFJ")                	var=(char*)"min #Delta #eta(FJ,FJ)";                  
    if(HistName=="HT")                  	var=(char*)"H_{T} [GeV]";                        
    if(HistName=="MJ")                  	var=(char*)"M_{J} [GeV]";                        
    if(HistName=="mj")                  	var=(char*)"m_{j} [GeV]";                        
    if(HistName=="mT")                  	var=(char*)"m_{T} [GeV]";                        
    if(HistName=="muspT")               	var=(char*)"p_{T}(muon) [GeV]";                  
    if(HistName=="muspTminusMET")          	var=(char*)"(MET-p_{T}(muon))/p_{T}(muon)";                  
    if(HistName=="musEta")              	var=(char*)"#eta(muon)";                         
    if(HistName=="musPhi")              	var=(char*)"#phi(muon)";                         
    if(HistName=="elspT")               	var=(char*)"p_{T}(electron) [GeV]";                  
    if(HistName=="elspTminusMET")          	var=(char*)"(MET-p_{T}(electron))/p_{T}(muon)";                  
    if(HistName=="elsEta")              	var=(char*)"#eta(electron)";                         
    if(HistName=="elsPhi")              	var=(char*)"#phi(electron)";                         
    if(HistName=="lepspT")               	var=(char*)"p_{T}(lepton) [GeV]";                  
    if(HistName=="lepsEta")              	var=(char*)"#eta(lepton)";                         
    if(HistName=="lepsPhi")              	var=(char*)"#phi(lepton)";                         
    if(HistName=="Nfatjet")             	var=(char*)"N_{fatjet}";
    if(HistName=="Nskinnyjet")          	var=(char*)"N_{skinny}";
    if(HistName=="Ncsvm")          	        var=(char*)"N_{CSVM}";
    if(HistName=="WpT")                 	var=(char*)"p_{T}(W) [GeV]";
    if(HistName=="mbb")                 	var=(char*)"Maximum m_{bb} [GeV]";

    TH1F *h1_DATA[7], *h1_T[7], *h1_TT_sl[7], *h1_TT_ll[7], *h1_WJets[7], *h1_DY[7], *h1_TTV[7], *h1_QCD[7],*h1_Others[7], *h1_MC[7]; 
    TH1F *h1_One[7], *h1_Ratio[7]; 
    TH1F *h1_f1500_100[7], *h1_f1200_800[7];
    THStack *st[7];
    //TCanvas *c = new TCanvas("c","c",1500,300);  
    //c->Divide(5,1);
    TCanvas *c = new TCanvas("c","c",1200,800);  
    c->Divide(3,2);
    TCanvas *c_AllFJ = new TCanvas("c_AllFJ","c_AllFJ",300,300);
    for(int i=2; i<7; i++) 
    {
        if(i!=6 && DrawOnlyAllFJ)  continue;

        h1_DATA[i]      = (TH1F*)HistFile->Get(Form("h1_DATA_%s_%ifatjet", HistName.Data(), i));  
        h1_T[i]         = (TH1F*)HistFile->Get(Form("h1_T_%s_%ifatjet", HistName.Data(), i));
        h1_TT_sl[i]     = (TH1F*)HistFile->Get(Form("h1_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h1_TT_ll[i]     = (TH1F*)HistFile->Get(Form("h1_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h1_WJets[i]     = (TH1F*)HistFile->Get(Form("h1_WJets_%s_%ifatjet", HistName.Data(), i));
        h1_DY[i]        = (TH1F*)HistFile->Get(Form("h1_DY_%s_%ifatjet", HistName.Data(), i)); 
        h1_TTV[i]       = (TH1F*)HistFile->Get(Form("h1_TTV_%s_%ifatjet", HistName.Data(), i)); 
        h1_QCD[i]       = (TH1F*)HistFile->Get(Form("h1_QCD_%s_%ifatjet", HistName.Data(), i)); 
        h1_Others[i]    = (TH1F*)HistFile->Get(Form("h1_Others_%s_%ifatjet", HistName.Data(), i)); 
        h1_f1500_100[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1500_100_%s_%ifatjet", HistName.Data(), i)); 
        h1_f1200_800[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1200_800_%s_%ifatjet", HistName.Data(), i)); 

        h1_One[i]       = new TH1F( Form("h1_One_%i",i), Form("h1_One_%i",i), 1,
                                    h1_DATA[i]->GetXaxis()->GetBinLowEdge(1),
                                    h1_DATA[i]->GetXaxis()->GetBinUpEdge(h1_DATA[i]->GetXaxis()->GetNbins())  );
        h1_One[i]->SetBinContent(1,1);

        // merge bins
        h1_DATA[i]->Rebin(NMergeBins);
        h1_T[i]->Rebin(NMergeBins);
        h1_TT_sl[i]->Rebin(NMergeBins);
        h1_TT_ll[i]->Rebin(NMergeBins);
        h1_WJets[i]->Rebin(NMergeBins);
        h1_DY[i]->Rebin(NMergeBins);
        h1_TTV[i]->Rebin(NMergeBins);
        h1_QCD[i]->Rebin(NMergeBins);
        h1_Others[i]->Rebin(NMergeBins);
        h1_f1500_100[i]->Rebin(NMergeBins);
        h1_f1200_800[i]->Rebin(NMergeBins);
        
        h1_MC[i] = (TH1F*)h1_TT_sl[i]->Clone(Form("h1_MC_%s_%ifatjet", HistName.Data(), i));
        h1_MC[i]->Add(h1_TT_ll[i]);
        h1_MC[i]->Add(h1_WJets[i]);
        h1_MC[i]->Add(h1_T[i]);
        h1_MC[i]->Add(h1_DY[i]);
        h1_MC[i]->Add(h1_TTV[i]);
        h1_MC[i]->Add(h1_QCD[i]);
        h1_MC[i]->Add(h1_Others[i]);
       
        float DataOverMC = h1_DATA[i]->Integral()/h1_MC[i]->Integral(); 
        if(NormToData) 
        { 
            h1_TT_sl[i]->Scale(DataOverMC);
            h1_TT_ll[i]->Scale(DataOverMC);
            h1_WJets[i]->Scale(DataOverMC);
            h1_T[i]->Scale(DataOverMC);
            h1_DY[i]->Scale(DataOverMC);
            h1_TTV[i]->Scale(DataOverMC);
            h1_QCD[i]->Scale(DataOverMC);
            h1_Others[i]->Scale(DataOverMC);
            h1_MC[i]->Scale(DataOverMC);
        }
        // add DY and QCD to Others 
        h1_Others[i]->Add(h1_DY[i]);
        h1_Others[i]->Add(h1_QCD[i]);

        h1_Ratio[i] = (TH1F*)h1_DATA[i]->Clone(Form("h1_Ratio_%s_%ifatjet", HistName.Data(), i));
        h1_Ratio[i]->Divide(h1_MC[i]);

        h1cosmetic(h1_DATA[i],          Form("DATA %ifatjet", i),               kBlack, 2, 0,           var);
        h1cosmetic(h1_TT_sl[i],         Form("TT(l) %ifatjet", i),              kBlack, 2, c_tt_1l,     var);
        h1cosmetic(h1_TT_ll[i],         Form("TT(ll) %ifatjet", i),             kBlack, 2, c_tt_2l,     var);
        h1cosmetic(h1_T[i],             Form("t+tW %ifatjet", i),               kBlack, 2, c_singlet,   var);
        h1cosmetic(h1_WJets[i],         Form("WJets %ifatjet", i),              kBlack, 2, c_wjets,     var);
        h1cosmetic(h1_Others[i],        Form("Others %ifatjet", i),             kBlack, 2, c_zjets,     var);
        h1cosmetic(h1_TTV[i],           Form("TTV %ifatjet", i),                kBlack, 2, c_other,     var);
        h1cosmetic(h1_f1500_100[i],     Form("T1tttt(1500,100) %ifatjet", i),   kRed,   2, 0,           var);
        h1cosmetic(h1_f1200_800[i],     Form("T1tttt(1200,800) %ifatjet", i),   kBlue,  2, 0,           var);
        h1cosmetic(h1_Ratio[i],         Form(" "),                              kBlack, 2, kBlack,      var);
        h1cosmetic(h1_One[i],           Form(" "),                              kGray,  2, 0,           var);

        bool DoLogOne = (DoLog && h1_MC[i]->Integral()>0);
        if(DrawOnlyAllFJ) c_AllFJ->cd();
        else 
        {
            c->cd(i-1);
        }
        //c->cd(i-1)->SetLeftMargin(0.15);
        //c->cd(i-1)->SetRightMargin(0.07);
        //c->cd(i-1)->SetBottomMargin(0.15);
        //c->cd(i-1)->SetTopMargin(0.1); 

        TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.04);
        pad1->SetRightMargin(0.1);
        pad1->SetLeftMargin(0.2);
        pad1->Draw();
        pad1->cd();
        pad1->cd()->SetLogy(DoLogOne);

        TString StackTitle = Form("%i fatjets", i);
        if(i==6) StackTitle = "All fatjets";
        if(i==5) StackTitle = "5+ fatjets";
        st[i] = new THStack( Form("Stack %ifatjet", i), StackTitle);
        st[i]->Add(h1_TTV[i]);
        st[i]->Add(h1_Others[i]);
        st[i]->Add(h1_WJets[i]);
        st[i]->Add(h1_T[i]);
        st[i]->Add(h1_TT_ll[i]);
        st[i]->Add(h1_TT_sl[i]);
        float HistMax = h1_MC[i]->GetMaximum()>h1_f1500_100[i]->GetMaximum()?h1_MC[i]->GetMaximum():h1_f1500_100[i]->GetMaximum();
        if(h1_DATA[i]->GetMaximum()>HistMax) 
        {  
            int MaxBin = h1_DATA[i]->GetMaximumBin();
            HistMax=h1_DATA[i]->GetBinContent(MaxBin) + h1_DATA[i]->GetBinError(MaxBin); 
        }
        //float HistMax = 4;
        st[i]->SetMaximum(HistMax*(DoLogOne?200:1.7));
        //st[i]->SetMinimum(h1_MC[i]->GetMinimum()*(DoLogOne?1:0));
        st[i]->SetMinimum((DoLogOne?0.005:0));
        st[i]->Draw("HIST"); 

        ///* needs to be in JaeStyle : FIXME  
        //st[i]->GetYaxis()->SetLabelSize(0.07); 
        //st[i]->GetYaxis()->SetTitleOffset(1.1); 
        //st[i]->GetYaxis()->SetTitleSize(0.07); 
        st[i]->GetYaxis()->SetTitle("Events/bin"); 
        st[i]->GetXaxis()->SetLabelSize(0.0); 
        st[i]->GetXaxis()->SetTitle(var); 
        //st[i]->GetXaxis()->SetTitleOffset(1.1); 
        //st[i]->GetXaxis()->SetTitleSize(0.07); 
        //*/
        h1_DATA[i]->SetLineColor(kBlack);
        h1_DATA[i]->SetMarkerColor(kBlack);
        h1_DATA[i]->SetMarkerSize(0.6);
        h1_DATA[i]->SetMarkerStyle(20);
        h1_DATA[i]->SetStats(0);
        if(doData) h1_DATA[i]->Draw("E SAME");

        TLegend *l1 = new TLegend(0.23, 0.60, 0.87, 0.86);
        l1->SetNColumns(2);
        l1->SetBorderSize(0);
        l1->SetFillColor(0);
        l1->SetFillStyle(0);
        l1->SetTextFont(42);
        l1->SetTextAlign(12);
        l1->SetTextSize(0.05);
        l1->SetFillColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->SetShadowColor(kWhite);
        if(doData) l1->AddEntry(h1_DATA[i],        " Data",  "lp");
        if(doData) l1->AddEntry(h1_DATA[i],        " ",  "");
        l1->AddEntry(h1_TT_sl[i],        " t#bar{t}(1#font[12]{l})",    "f");
        l1->AddEntry(h1_TT_ll[i],        " t#bar{t}(2#font[12]{l})",    "f");
        l1->AddEntry(h1_T[i],           " t+tW",  "f");
        l1->AddEntry(h1_WJets[i],       " WJets", "f");
        l1->AddEntry(h1_Others[i],      " DY/QCD/others",    "f");
        l1->AddEntry(h1_TTV[i],         " t#bar{t}V",    "f");
        l1->AddEntry(h1_f1500_100[i],   " T1tttt[1500,100]", "l");
        l1->AddEntry(h1_f1200_800[i],   " T1tttt[1200,800]", "l");
        l1->Draw();
        
        h1_f1500_100[i]->Scale(SignalScale);
        h1_f1200_800[i]->Scale(SignalScale);
        h1_f1500_100[i]->Draw("SAME HIST");
        h1_f1200_800[i]->Draw("SAME HIST");

        // CMS Labels 
        float textSize = 0.05;

        //TLatex *TexEnergyLumi = new TLatex(0.9,0.92,Form("#sqrt[]{s}=13 TeV, L = %.1f pb^{-1}", Lumi));
        TLatex *TexEnergyLumi = new TLatex(0.9,0.92,Form("#font[42]{%.1f fb^{-1} (13 TeV)}", Lumi/1000));
        TexEnergyLumi->SetNDC();
        TexEnergyLumi->SetTextSize(textSize);
        TexEnergyLumi->SetTextAlign (31);
        TexEnergyLumi->SetLineWidth(2);

        TLatex *TexCMS = new TLatex(0.2,0.92,"CMS #font[52]{Preliminary}");
        TexCMS->SetNDC();
        TexCMS->SetTextSize(textSize);
        TexCMS->SetLineWidth(2);
        
        TLatex *SF = new TLatex(0.88,0.84,Form("#font[42]{%s, Data/MC=%.1f%%}", Selection, DataOverMC*100));
        SF->SetNDC();
        SF->SetTextSize(textSize*0.8);
        SF->SetTextAlign (31);
        SF->SetLineWidth(2);
       
        // FIXME : need to add lepton flavor
        TString LabelExt = Form("N_{fatjet} = %i", i);
        if(i==5) LabelExt="N_{fatjet} >= 5";
        TLatex *TexExt = new TLatex(0.85,0.7,LabelExt);
        TexExt->SetTextAlign (31);
        TexExt->SetNDC();
        TexExt->SetTextSize(textSize);
        TexExt->SetLineWidth(2);
        
        TexEnergyLumi->Draw("SAME");
        TexCMS->Draw("SAME");
        if(i!=6) TexExt->Draw("SAME"); 
        SF->Draw("SAME");

        if(DrawOnlyAllFJ) c_AllFJ->cd();
        else 
        {
            c->cd(i-1);
        }
        TPad *pad2 = new TPad("p_pull", "p_pull", 0.0, 0.0, 1.0, 0.3);
        pad2->SetLeftMargin(0.2);
        pad2->Draw();
        pad2->cd();
        pad2->SetTopMargin(0.04);
        pad2->SetRightMargin(0.1);
        pad2->SetBottomMargin(0.4);
        
        h1_One[i]->SetLabelSize(0.16,"XY");
        h1_One[i]->SetTitleSize(0.16,"XY");
        h1_One[i]->SetTitleOffset(1.0);
        h1_One[i]->GetYaxis()->SetNdivisions(3,true);
        h1_One[i]->GetXaxis()->SetNdivisions(5,true);
        h1_One[i]->SetMinimum(0);
        h1_One[i]->SetMaximum(2);
        h1_One[i]->SetYTitle("Data/MC");
        //h1_One[i]->GetYaxis()->SetTitleOffset(1.0);
        h1_One[i]->Draw("HIST");
        //h1_Ratio[i]->SetMinimum(0);
        //h1_Ratio[i]->SetMaximum(2);
        h1_Ratio[i]->SetMarkerStyle(20);
        h1_Ratio[i]->SetMarkerSize(0.6);
        //h1_Ratio[i]->SetLabelSize(0.1,"XY");
        //h1_Ratio[i]->SetTitleSize(0.1,"XY");
        //h1_Ratio[i]->SetTitleOffset(1.5);
        h1_Ratio[i]->Draw("SAME E");
    }

    // 
    if(HistName=="mj") HistName="JetMass";
    if(DrawOnlyAllFJ)
    {
        c_AllFJ->Print( Form("Figures/%s/CompareDataMC_AllFJ_%s_%s%s.pdf", Selection, HistName.Data(), Selection, DoLog?"_log":"") ); 
    }
    else 
    { 
        c->Print( Form("Figures/%s/CompareDataMC_%s_%s%s.pdf", Selection, HistName.Data(), Selection, DoLog?"_log":"") ); 
    } 
    // 
    HistFile->Close();
    delete c; 
    delete c_AllFJ; 

}
