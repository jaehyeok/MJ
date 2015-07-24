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
#include "TPad.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TInterpreter.h"
#include "TLatex.h"
#include "TColor.h"
#include "babytree.h"

using namespace std;
int pTR0p5thres     = 30;
bool doData         = 1;
bool DoLog          = 0;
int SignalScale     = 1000;
bool RemoveMuon     = 0;
bool Do3FatjetScale = 0;
int FatjetpTthres   = 50;

// Btagging prob map function
TFile *f_tageff_;

//
//TH1F initialization
//
TH1F* InitTH1F(char* Name, char* Title, int Nbins, double XMin, double XMax){
    TH1F *h1 = new TH1F(Name, Title, Nbins, XMin, XMax);
    h1->Sumw2();
    return h1;
}
TH1F* InitTH1F(char* Name, char* Title, int Nbins, Float_t* xbins){
    TH1F *h1 = new TH1F(Name, Title, Nbins, xbins);
    h1->Sumw2();
    return h1;
}
//
//TH2F initialization
//
TH2F* InitTH2F(char* Name, char* Title, int NXbins, double XMin, double XMax, int NYbins, double YMin, double YMax){
    TH2F *h2 = new TH2F(Name, Title, NXbins, XMin, XMax, NYbins, YMin, YMax);
    h2->Sumw2();
    return h2;
}
//
// Fill TH1F
//
TH1F* FillTH1F(TH1F* &h1, double var, double weight){
    if(var >= h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())) 
        var=h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins())-0.00001;
    if(var < h1->GetXaxis()->GetBinLowEdge(1)) 
        var=h1->GetXaxis()->GetBinLowEdge(1)+0.00001;
    h1->Fill(var, weight);        
    return h1;
}
//
// Fill TH2F
//
TH2F* FillTH2F(TH2F* &h2, double varX, double varY, double weight){
    if(varX >= h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins())) 
        varX=h2->GetXaxis()->GetBinUpEdge(h2->GetXaxis()->GetNbins())-0.00001;
    if(varY >= h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins())) 
        varY=h2->GetYaxis()->GetBinUpEdge(h2->GetYaxis()->GetNbins())-0.00001;
    if(varX < h2->GetXaxis()->GetBinLowEdge(1)) 
        varX=h2->GetXaxis()->GetBinLowEdge(1)+0.00001;
    if(varY < h2->GetYaxis()->GetBinLowEdge(1)) 
        varY=h2->GetYaxis()->GetBinLowEdge(1)+0.00001;
    h2->Fill(varX, varY, weight);        
    return h2;
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
// h2 cosmetics
//
void h2cosmetic(TH2F* &h2, char* title, TString Xvar="", TString Yvar="", TString Zvar="Events/bin")
{
    h2->SetTitle(title);
    h2->SetXTitle(Xvar);
    h2->SetYTitle(Yvar);
    h2->SetZTitle(Zvar);
    h2->SetStats(0);
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


float GetMuonEff(/*float pt=1,*/ float eta=0)
{ 
   float eff = 1; 
/* 
   // For SingleMu dataset
   if(pt<25) eff = 0;
   
   if(pt>25 && pt<30)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.9;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.81;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.78;
   }
   
   if(pt>30 && pt<35)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.91;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.83;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.80;
   }
   
   if(pt>35 && pt<40)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.92;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.84;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.80;
   }
   
   if(pt>40 && pt<50)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.94;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.85;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.82;
   }
   
   if(pt>50)  
   { 
        if(TMath::Abs(eta)<0.8) eff = 0.94;
        if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.2) eff = 0.86;
        if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<2.1) eff = 0.82;
   }

   if(TMath::Abs(eta)>2.1) eff = 0; 

   eff = 0.88;
*/

    if(TMath::Abs(eta)<0.9) eff = 0.98;
    else eff = 0.84;
   
    return eff;
}

float getDPhi(float phi1, float phi2) 
{ 
    float absdphi = abs(phi1-phi2);
    if(absdphi < TMath::Pi()) return absdphi;
    else return (2*TMath::Pi() - absdphi);
    // Ryan's : http://www.aidansean.com/cheatsheets/?p=105
    //return TMath::Abs(TMath::Abs(TMath::Abs(phi1 - phi2) - TMath::Pi())-TMath::Pi());  

}
float getDEta(float eta1, float eta2)
{
    return abs(eta1-eta2);
}
float getDR(float dphi, float deta)
{
    return TMath::Sqrt(dphi*dphi+deta*deta);
}

float getDR(float eta1, float eta2, float phi1, float phi2)
{
    return getDR(getDPhi(phi1, phi2), eta1-eta2);
}

float getISRSF(float ISRpT)
{ 
    // ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMST2ccMadgraph8TeV
    if( ISRpT>250) return 0.8;
    else if( ISRpT>150 ) return 0.9;
    else if( ISRpT>120 ) return 0.95;
    else return 1;
}

float getMJSF(float MJbin)
{
    if(MJbin>500) return 1.60;
    else if(MJbin>400) return 1.33;
    else if(MJbin>300) return 1.12;
    else if(MJbin>200) return 1.18;
    else if(MJbin>100) return 1.00;
    else return 0.84;
    
    // FNAL workshop
    //if(MJ>500) return 1.8;
    //else if(MJ>400) return 1.5;
    //else if(MJ>300) return 1.2;
    //else if(MJ>200) return 1.3;
    //else if(MJ>100) return 1.1;
    //else return 0.5;
}
/*
//
// Btagging 
//
// here 
TString GetJetTagEffMaps(TString process)
{ 
    if(process=="TT_sl") 
        return "histos_btageff_csvm_TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71.root"; 
    else if(process=="TT_ll") 
        return "histos_btageff_csvm_TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71.root";
    else if(process=="DY") 
        return "histos_btageff_csvm_DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1535_v66.root"; 
    else if(process=="WJets") 
        return "histos_btageff_csvm_DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1535_v66.root"; 
        //return "histos_btageff_csvm_WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1_AODSIM_UCSB1464ra2b_v65.root"; 
    else if(process=="T") 
        return "histos_btageff_csvm_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71.root"; 
    else 
    { 
        cout << "[CAUTION] No btagging mapping file exists!! " << endl;
        return "histos_btageff_csvm_TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71.root";
    }
}

void LoadJetTagEffMaps(TString process) 
{
    int cfAVersion = 71;

    if (cfAVersion<=71) 
    {
        //assert(f_tageff_ ==0);
        TString filename = GetJetTagEffMaps(process); 
        //TString filename = "histos_btageff_csvm_TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71.root";
        filename.Prepend("btagEffMaps/");
        f_tageff_ = new TFile(filename,"READ");
        if (f_tageff_->IsZombie()) 
        {
            //cout<<"Failed to load the b-tag eff map for sample "<<sampleName<<endl;
            //comment this next line to ALLOW JOBS TO RUN WITHOUT BTAG EFF
            //if (!isSampleRealData()) assert(0); //it's just too frustrating to let jobs continue in this state
            delete f_tageff_;
            f_tageff_=0;
        }
        else 
        {
            TH1D * h_btageff = static_cast<TH1D*>(f_tageff_->Get("h_btageff"));
            int nbins=h_btageff->GetNbinsX();
            if (nbins != 17) 
            {
                cout<<" b tag eff map has the wrong number of bins! "<<nbins<<endl;
                assert(0);
            }
            else 
            {
                std::cout << "Successfully loaded b-tag eff maps in " << filename << endl;
            }
        }
    }
}

//get MC btag efficiency
double GetJetTagEff(unsigned int ijet, TH1D* h_btageff, TH1D* h_ctageff, TH1D* h_ltageff) {

    double tageff=0;
    const float pt = JetPt->at(ijet);
    const float eta = fabs(JetEta->at(ijet));
    int flavor = static_cast<int>(JetMCId->at(ijet));

    //x is the pt value that will be used to evaluate the SF
    //the max pt depends on flavor, and, for LF jets, eta
    const double cutoff1=800.; const double cutoff2=700.;
    double x;
    //HF or central LF, max is 800
    if ( abs(flavor)==5 || abs(flavor)==4 || eta<1.5) x = pt > cutoff1 ? cutoff1 : pt;
    //high eta LF, max is 700
    else x = pt > cutoff2 ? cutoff2 : pt;


        //  if (theBTagEffType_ == kBTagEff05 || theBTagEffType_ == kBTagEffup5 || theBTagEffType_ == kBTagEffdown5) { //new  BTV-11-004 prescription 

        if (abs(flavor) ==4 || abs(flavor)==5) { //heavy flavor
            double errFactor = 1;

            if (pt >cutoff1) { //use twice the errors
                errFactor=2;
            }
            if (abs(flavor) == 4)  errFactor=2; //charm has double the errors   "SFc = SFb with twice the quoted uncertainty"
            //not clear to me what to do for charm with pT> cutoff. errFactor of 2 or 4? must be small effect though

            // Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
            double  SFb = 0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)));

            //apply FASTSTIM correction where needed
            //SFb *= GetbJetFastsimSF("value",flavor,pt); // FIXME:for signal 

            // skip syst variations, for now

            // cout<<"jet flavor, pt, SF = "<<abs(flavor)<<" "<<pt<<" "<<SFb<<endl;
            if      (abs(flavor) == 5) tageff = SFb * h_btageff->GetBinContent( h_btageff->FindBin( pt ) );
            else if (abs(flavor) == 4) tageff = SFb * h_ctageff->GetBinContent( h_ctageff->FindBin( pt ) );
            else assert(0);
        } // if heavy flavor
        else { //light flavor [ see https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_Moriond2013.C ]
            //note that these are valid only to a cutoff value, so use 'x' and not 'pt'
            double SF=0;
            double nominalSF=0;
            if ( eta < 0.8 ) {
                nominalSF =  ((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x))); // SF without + or - variation in mistag rate
                SF = nominalSF;
            }
            else if (eta>=0.8 && eta<1.6) {
                nominalSF = ((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)));
                SF = nominalSF;
            }
            else if (eta>=1.6 && eta<=2.4) {
                nominalSF = ((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x)));
                SF = nominalSF;
            }
            //design question -- what do to for jets at eta>2.4? assert? or return tageff=0?
            //i guess tageff 0 makes the most sense, so leave SF = 0

            // double deltaSF = SF - nominalSF; //here is the nominal uncertainty on the fullsim SF

            //5 June 2013 -- new code to apply LF mistags
            //apply FASTSTIM correction where needed
            //SF *= GetbJetFastsimSF("value",flavor,pt);

            // now, skip the extra uncertainties for high pT/eta and Fastsim

            tageff = SF * h_ltageff->GetBinContent( h_ltageff->FindBin( pt )); 
            // cout<<"jet flavor, pt, SF = "<<abs(flavor)<<" "<<pt<<" "<<SF<<endl;

        } // if light flavor

        //    } //if BTV-11-004 prescription
        //  else {      //no longer support old prescriptions
        //  assert(0);
        //    }



    if (tageff<0) tageff=0;
    if (tageff>1) tageff=1;

    return tageff;
}



// This gives you the weight you apply to the MC when you don't cut on b-tags
void CalculateTagProbs(double &Prob0, double &ProbGEQ1, double &Prob1, double &ProbGEQ2,
                       double &Prob2, double &ProbGEQ3, double &Prob3, double &ProbGEQ4, 
                       bool isData) 
{

    TString sampleName="TTJets";

    //must initialize correctly
    Prob2 = 0;
    Prob1 = 0; ProbGEQ1 = 1; Prob0 = 1; ProbGEQ2 = 0;
    Prob3 = 0; ProbGEQ4 = 0;

    if(isData) 
    { // for data, set all to 1--not really needed
        Prob0=1;
        ProbGEQ1=1;
        Prob1=1;
        ProbGEQ2=1;
        Prob2=1;
        ProbGEQ3=1;
        Prob3=1;
        ProbGEQ4=1;
        return;
    }

    if (f_tageff_ == 0) 
    { //if the b-tag eff file is not there make sure it will be obvious to the user -- fill all with zero
        Prob0=0;
        ProbGEQ1=0;
        Prob1=0;
        ProbGEQ2=0;
        Prob2=0;
        ProbGEQ3=0;
        Prob3=0;
        ProbGEQ4=0;
        return;
    }

    char btageffname[200], ctageffname[200], ltageffname[200];
    std::string sbtageff = "h_btageff";  std::string sctageff = "h_ctageff";  std::string sltageff = "h_ltageff";
    sprintf(btageffname,"%s",sbtageff.c_str());   
    sprintf(ctageffname,"%s",sctageff.c_str());   
    sprintf(ltageffname,"%s",sltageff.c_str());   
    TH1D * h_btageff  = static_cast<TH1D*>(f_tageff_->Get(btageffname));
    TH1D * h_ctageff  = static_cast<TH1D*>(f_tageff_->Get(ctageffname));
    TH1D * h_ltageff  = static_cast<TH1D*>(f_tageff_->Get(ltageffname));

    for (unsigned int ijet=0; ijet<JetPt->size(); ++ijet) {
        double subprob1=0;
        double subprob2=0;

        double effi = GetJetTagEff(ijet, h_btageff, h_ctageff, h_ltageff);
        //      cout << "jet: " << ijet << ", effi: " << effi << endl;
        Prob0 = Prob0* ( 1 - effi);

        double product = 1;
        for (unsigned int kjet=0; kjet<JetPt->size(); ++kjet) 
        {
            double effk = GetJetTagEff(kjet, h_btageff, h_ctageff, h_ltageff);
            if(kjet != ijet) product = product*(1-effk);
            if(kjet > ijet)
            {
                double subproduct = 1;
                for (unsigned int jjet=0; jjet<JetPt->size(); ++jjet) 
                {
                    if(jjet != kjet && jjet != ijet)
                    {
                        double effj = GetJetTagEff(jjet, h_btageff, h_ctageff, h_ltageff);
                        subproduct = subproduct*(1-effj);
                        if ( jjet > kjet) 
                        {
                            double subproduct2 = 1;
                            for (unsigned int ljet=0; ljet<JetPt->size(); ++ljet) 
                            {
                                if(ljet != kjet && ljet != ijet && ljet != jjet)
                                {
                                    double effl = GetJetTagEff(ljet, h_btageff, h_ctageff, h_ltageff);
                                    subproduct2 = subproduct2*(1-effl);
                                }
                            } //ljet loop
                            subprob2 += effk*effj*subproduct2;
                        }
                    }
                }//j loop
                subprob1 += effk*subproduct;
            }
        }//k loop

        Prob1 += effi*product;
        Prob2 += effi*subprob1;
        Prob3 += effi*subprob2; 

    }

    //  std::cout << "prob0 = " << Prob0 << ", prob1 = " << Prob1 << ", prob2 = " << Prob2 << std::endl;  

    ProbGEQ1 = 1 - Prob0;
    ProbGEQ2 = 1 - Prob1 - Prob0;
    ProbGEQ3 = 1 - Prob2 - Prob1 - Prob0;
    ProbGEQ4 = 1 - Prob3 - Prob2 - Prob1 - Prob0;

    return;

}
*/
//
// Print A line of table
//
void PrintTableOneLine(TString Process, TH1F* h1[6], int lepflav=0, bool doLatex=true)
{ 
        if(doLatex)
        {
            if(lepflav!=11 && lepflav!=13)
            {
                Double_t error[6];
                for(int i=2; i<6; i++) h1[i]->IntegralAndError(1,10000,error[i]);
                cout << Process << " & " 
                     //<< Form("$%.2f \\pm %.2f$",h1[2]->Integral(),error[2]) << "\t&" 
                     //<< Form("$%.2f \\pm %.2f$",h1[3]->Integral(),error[3]) << "\t&" 
                     //<< Form("$%.2f \\pm %.2f$",h1[4]->Integral(),error[4]) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[5]->Integral(),error[5]) << "\t\\\\" 
                     << endl; 
            }
            else 
            {   
                int bin = (lepflav-9)/2;
                cout << Process << " & " 
                     //<< Form("$%.2f \\pm %.2f$",h1[2]->GetBinContent(bin),h1[2]->GetBinError(bin)) << "\t&" 
                     //<< Form("$%.2f \\pm %.2f$",h1[3]->GetBinContent(bin),h1[3]->GetBinError(bin)) << "\t&" 
                     //<< Form("$%.2f \\pm %.2f$",h1[4]->GetBinContent(bin),h1[4]->GetBinError(bin)) << "\t&" 
                     << Form("$%.2f \\pm %.2f$",h1[5]->GetBinContent(bin),h1[5]->GetBinError(bin)) << "\t\\\\" 
                     << endl; 
            }
        }
        else 
        {   
            if(lepflav!=11 && lepflav!=13)
            {
                Double_t error[6];
                for(int i=2; i<6; i++) h1[i]->IntegralAndError(1,10000,error[i]);
                cout << "|" << 
                    setw(20) << Process << " |"  <<
                    setw(20) << Form("%.2f +/- %.2f", h1[2]->Integral(), error[2]) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[3]->Integral(), error[3]) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[4]->Integral(), error[4]) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[5]->Integral(), error[5]) << " |" << endl; 
            }
            else 
            {   
                int bin = (lepflav-9)/2;
                cout << "|" << 
                    setw(20) << Process << " |"  <<
                    setw(20) << Form("%.2f +/- %.2f", h1[2]->GetBinContent(bin), h1[2]->GetBinError(bin)) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[3]->GetBinContent(bin), h1[3]->GetBinError(bin)) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[4]->GetBinContent(bin), h1[4]->GetBinError(bin)) << " |" <<
                    setw(20) << Form("%.2f +/- %.2f", h1[5]->GetBinContent(bin), h1[5]->GetBinError(bin)) << " |" << endl; 
            }
        }
        
        // 
        if(Process.Contains("Bkgd") || Process.Contains("DATA")) 
        { 
            if(doLatex) 
            {   
                cout << "\\hline" << endl;  
            }
            else 
            { 
                cout << "|" << 
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "--" <<
                    setw(20) << Form("--------------------") << "-|" << endl;
            }
        }
}

//
// per process
//
void DoOneProcess(TChain *ch, int pTR0p5/*, int Nbcut*/, float METlow, float METhigh, int njlow, int njhigh,int nblow, int nbhigh, float MJcut) 
{ 
/*
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
    vector<bool>    *filter;
    vector<float>   *mj;
    vector<float>   *FatjetPt;
    vector<float>   *FatjetEta;
    vector<float>   *FatjetPhi;
    vector<float>   *RA4MusPt;
    vector<float>   *RA4MusPhi;
    vector<float>   *RA4MusEta;
    vector<float>   *RA4ElsPt;
    vector<float>   *RA4ElsPhi;
    vector<float>   *RA4ElsEta;
    vector<float>   *JetPt;
    vector<float>   *JetEta;
    vector<float>   *JetPhi;
    vector<float>   *JetCSV;
    vector<float>   *RA4MusVetoPt;
    vector<float>   *RA4ElsVetoPt;

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
    TBranch         *b_filter;   //!
    TBranch         *b_mj;   //!
    TBranch         *b_FatjetPt;   //!
    TBranch         *b_FatjetEta;   //!
    TBranch         *b_FatjetPhi;   //!
    TBranch         *b_RA4MusPt;   //!
    TBranch         *b_RA4MusPhi;   //!
    TBranch         *b_RA4MusEta;   //!
    TBranch         *b_RA4ElsPt;   //!
    TBranch         *b_RA4ElsPhi;   //!
    TBranch         *b_RA4ElsEta;   //!
    TBranch         *b_JetPt;   //!
    TBranch         *b_JetEta;   //!
    TBranch         *b_JetPhi;   //!
    TBranch         *b_JetCSV;   //!
    TBranch         *b_RA4MusVetoPt;   //!
    TBranch         *b_RA4ElsVetoPt;   //!
    filter        = 0;
    mj        = 0;
    FatjetPt  = 0;
    FatjetEta = 0;
    FatjetPhi = 0;
    RA4MusPt  = 0;
    RA4MusPhi = 0;
    RA4MusEta = 0;
    RA4ElsPt  = 0;
    RA4ElsPhi  = 0;
    RA4ElsEta  = 0;
    JetPt  = 0;
    JetEta = 0;
    JetPhi = 0;
    JetCSV = 0;
    RA4MusVetoPt  = 0;
    RA4ElsVetoPt  = 0;

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
    ch->SetBranchAddress("filter", &filter, &b_filter);
    ch->SetBranchAddress(Form("mj_pT%i", pTR0p5), &mj, &b_mj);
    ch->SetBranchAddress(Form("FatjetPt_pT%i", pTR0p5), &FatjetPt, &b_FatjetPt);
    ch->SetBranchAddress(Form("FatjetEta_pT%i", pTR0p5), &FatjetEta, &b_FatjetEta);
    ch->SetBranchAddress(Form("FatjetPhi_pT%i", pTR0p5), &FatjetPhi, &b_FatjetPhi);
    ch->SetBranchAddress("RA4MusPt", &RA4MusPt, &b_RA4MusPt);
    ch->SetBranchAddress("RA4MusPhi", &RA4MusPhi, &b_RA4MusPhi);
    ch->SetBranchAddress("RA4MusEta", &RA4MusEta, &b_RA4MusEta);
    ch->SetBranchAddress("RA4ElsPt", &RA4ElsPt, &b_RA4ElsPt);
    ch->SetBranchAddress("RA4ElsPhi", &RA4ElsPhi, &b_RA4ElsPhi);
    ch->SetBranchAddress("RA4ElsEta", &RA4ElsEta, &b_RA4ElsEta);
    ch->SetBranchAddress("JetPt",  &JetPt,  &b_JetPt);
    ch->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
    ch->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
    ch->SetBranchAddress("JetCSV", &JetCSV, &b_JetCSV);
    ch->SetBranchAddress("RA4MusVetoPt", &RA4MusVetoPt, &b_RA4MusVetoPt);
    ch->SetBranchAddress("RA4ElsVetoPt", &RA4ElsVetoPt, &b_RA4ElsVetoPt);
*/    
    InitBaby(ch, pTR0p5); 
    TString ChainName = ch->GetTitle();
    cout << "[MJ Analysis] " << ChainName << endl;  

    bool        TrigMuon;
    bool        TrigSingleMuon;
    bool        TrigElectron;
    TBranch     *b_TrigMuon; 
    TBranch     *b_TrigSingleMuon; 
    TBranch     *b_TrigElectron; 
    if(ChainName.Contains("DATA"))  
    { 
        if(ChainName.Contains("MuHad")) 
        {
            ch->SetBranchAddress("TrigMuon", &TrigMuon, &b_TrigMuon);
            ch->SetBranchAddress("TrigSingleMuon", &TrigSingleMuon, &b_TrigSingleMuon);
        }
        if(ChainName.Contains("EHad")) 
        {
            ch->SetBranchAddress("TrigElectron", &TrigElectron, &b_TrigElectron);
        }
    }
    
    float        top1pT;
    float        top1Phi;
    float        top2pT;
    float        top2Phi;
    TBranch     *b_top1pT; 
    TBranch     *b_top1Phi; 
    TBranch     *b_top2pT; 
    TBranch     *b_top2Phi; 
    if(ChainName.Contains("TT") )  
    { 
        ch->SetBranchAddress("top1pT", &top1pT, &b_top1pT);
        ch->SetBranchAddress("top1Phi", &top1Phi, &b_top1Phi);
        ch->SetBranchAddress("top2pT", &top2pT, &b_top2pT);
        ch->SetBranchAddress("top2Phi", &top2Phi, &b_top2Phi);
    }

    //
    //
    //
    TH1F *h1_yields[6];
    TH1F *h1_ttbarpT[6];
    TH1F *h1_MJ[6], *h1_mT[6], *h1_Nskinnyjet[6], *h1_muspT[6], *h1_muspTminusMET[6], *h1_musEta[6], *h1_musPhi[6], 
         *h1_mj[6], *h1_FatjetPt[6], *h1_FatjetEta[6], *h1_mjFJFJ[6], 
         *h1_mjFJFJnotclose[6], *h1_FatjetPtFJFJnotclose[6], *h1_FatjetEtaFJFJnotclose[6],  *h1_FatjetPhiFJFJnotclose[6],
         *h1_FatjetPtFJFJcloseHigherPt[6], *h1_FatjetEtaFJFJcloseHigherPt[6],  *h1_FatjetPhiFJFJcloseHigherPt[6], 
         *h1_FatjetPtFJFJcloseLowerPt[6], *h1_FatjetEtaFJFJcloseLowerPt[6],  *h1_FatjetPhiFJFJcloseLowerPt[6], 
         *h1_FatjetPt1[6],  *h1_FatjetPt2[6],  *h1_FatjetPt3[6],  *h1_FatjetPt4[6],
         *h1_FatjetPhi1[6], *h1_FatjetPhi2[6], *h1_FatjetPhi3[6], *h1_FatjetPhi4[6],
         *h1_FatjetEta1[6], *h1_FatjetEta2[6], *h1_FatjetEta3[6], *h1_FatjetEta4[6],
         *h1_mj_mumatch[6], *h1_FatjetPt_mumatch[6], *h1_mj_notmumatch[6], *h1_FatjetPt_notmumatch[6], 
         *h1_BjetPt_bmatchmu[6], *h1_BjetPt_notbmatchmu[6], *h1_dRbmu_bmatchmu[6], *h1_dRbmu_notbmatchmu[6],  
         *h1_dRFJ[6], *h1_dPhiFJ[6], *h1_dEtaFJ[6],
         *h1_mj_bmatch[6], *h1_FatjetPt_bmatch[6], *h1_mj_notbmatch[6], *h1_FatjetPt_notbmatch[6], 
         *h1_mj_metmatch[6], *h1_FatjetPt_metmatch[6], *h1_mj_notmetmatch[6], *h1_FatjetPt_notmetmatch[6], 
         *h1_HT[6], *h1_MET[6], *h1_MET_mTlt140[6], *h1_MET_mTgt140[6], *h1_dPhi_mTlt140[6], *h1_dPhi_mTgt140[6], *h1_METPhi[6], *h1_METx[6], *h1_METy[6], *h1_DPhi[6], *h1_Nfatjet[6], *h1_WpT[6];
    TH2F *h2_mjvsFatjetPt[6];
    TH2F *h2_mindRFJmjvsmj[6];
    TH2F *h2_MJmT[6];
    TH2F *h2_dPhivsMET[6][4];
    TH2F *h2_dPhiTTvsMET[6][4];

    for(int i=0; i<6; i++) 
    {
        h1_yields[i] = InitTH1F( Form("h1_%s_yields_%ifatjet", ch->GetTitle(), i), 
                                 Form("h1_%s_yields_%ifatjet", ch->GetTitle(), i), 
                                 2, 0, 2); // bin1 for e, bin2 for mu
        Float_t MJbinning[7] = {0,100,200,300,400,500,1000};
        h1_MJ[i] = InitTH1F( Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MJ_%ifatjet", ch->GetTitle(), i), 
                             6, MJbinning);
                             //2, 0, 800);
        h1_mj[i] = InitTH1F( Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mj_%ifatjet", ch->GetTitle(), i), 
                             40, 0, 500);
        h1_ttbarpT[i] = InitTH1F( Form("h1_%s_ttbarpT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_ttbarpT_%ifatjet", ch->GetTitle(), i), 
                             40, 0, 500);
        h1_mjFJFJ[i] = InitTH1F( Form("h1_%s_mjFJFJ_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mjFJFJ_%ifatjet", ch->GetTitle(), i), 
                             80, 0, 800);
        h1_mjFJFJnotclose[i] = InitTH1F( Form("h1_%s_mjFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mjFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                             40, 0, 500);
        h1_FatjetPtFJFJnotclose[i] = InitTH1F( Form("h1_%s_FatjetPtFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_FatjetPtFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                             40, 0, 500);
        h1_FatjetEtaFJFJnotclose[i] = InitTH1F( Form("h1_%s_FatjetEtaFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetEtaFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                            30, -3, 3);
        h1_FatjetPhiFJFJnotclose[i] = InitTH1F( Form("h1_%s_FatjetPhiFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetPhiFJFJnotclose_%ifatjet", ch->GetTitle(), i), 
                            20, -TMath::Pi(), TMath::Pi());
        h1_FatjetEtaFJFJcloseHigherPt[i] = InitTH1F( Form("h1_%s_FatjetEtaFJFJcloseHigherPt_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetEtaFJFJcloseHigherPt_%ifatjet", ch->GetTitle(), i), 
                            30, -3, 3);
        h1_FatjetPhiFJFJcloseHigherPt[i] = InitTH1F( Form("h1_%s_FatjetPhiFJFJcloseHigherPt_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetPhiFJFJcloseHigherPt_%ifatjet", ch->GetTitle(), i), 
                            20, -TMath::Pi(), TMath::Pi());
        h1_FatjetPtFJFJcloseHigherPt[i] = InitTH1F( Form("h1_%s_FatjetPtFJFJcloseHigherPt_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetPtFJFJcloseHigherPt_%ifatjet", ch->GetTitle(), i), 
                            40, 0, 500);
        h1_FatjetEtaFJFJcloseLowerPt[i] = InitTH1F( Form("h1_%s_FatjetEtaFJFJcloseLowerPt_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetEtaFJFJcloseLowerPt_%ifatjet", ch->GetTitle(), i), 
                            30, -3, 3);
        h1_FatjetPhiFJFJcloseLowerPt[i] = InitTH1F( Form("h1_%s_FatjetPhiFJFJcloseLowerPt_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetPhiFJFJcloseLowerPt_%ifatjet", ch->GetTitle(), i), 
                            20, -TMath::Pi(), TMath::Pi());
        h1_FatjetPtFJFJcloseLowerPt[i] = InitTH1F( Form("h1_%s_FatjetPtFJFJcloseLowerPt_%ifatjet", ch->GetTitle(), i), 
                            Form("h1_%s_FatjetPtFJFJcloseLowerPt_%ifatjet", ch->GetTitle(), i), 
                            40, 0, 500);
        h1_mj_mumatch[i] = InitTH1F( Form("h1_%s_mj_mumatch_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_mj_mumatch_%ifatjet", ch->GetTitle(), i), 
                                   40, 0, 500);
        h1_mj_notmumatch[i] = InitTH1F( Form("h1_%s_mj_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                      Form("h1_%s_mj_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                      40, 0, 500);
        h1_mj_bmatch[i] = InitTH1F( Form("h1_%s_mj_bmatch_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_mj_bmatch_%ifatjet", ch->GetTitle(), i), 
                                   40, 0, 500);
        h1_mj_notbmatch[i] = InitTH1F( Form("h1_%s_mj_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                      Form("h1_%s_mj_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                      40, 0, 500);
        h1_mj_metmatch[i] = InitTH1F( Form("h1_%s_mj_metmatch_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_mj_metmatch_%ifatjet", ch->GetTitle(), i), 
                                   40, 0, 500);
        h1_mj_notmetmatch[i] = InitTH1F( Form("h1_%s_mj_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                      Form("h1_%s_mj_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                      40, 0, 500);
        h1_Nfatjet[i] = InitTH1F( Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  Form("h1_%s_Nfatjet_%ifatjet", ch->GetTitle(), i), 
                                  11, -0.5, 10.5);
        h1_Nskinnyjet[i] = InitTH1F( Form("h1_%s_Nskinnyjet_%ifatjet", ch->GetTitle(), i), 
                                     Form("h1_%s_Nskinnyjet_%ifatjet", ch->GetTitle(), i), 
                                     7, 5.5, 12.5);
        Float_t mTbinning[9] = {0,30,60,90,110,140,180,300,500};
        h1_mT[i] = InitTH1F( Form("h1_%s_mT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_mT_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 600);
                             8, mTbinning);
        h1_muspT[i] = InitTH1F( Form("h1_%s_muspT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_muspT_%ifatjet", ch->GetTitle(), i), 
                             20, 0, 200);
        h1_muspTminusMET[i] = InitTH1F( Form("h1_%s_muspTminusMET_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_muspTminusMET_%ifatjet", ch->GetTitle(), i), 
                             //20, 100, 400);
                             20, -2, 8);
        h1_musEta[i] = InitTH1F( Form("h1_%s_musEta_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_musEta_%ifatjet", ch->GetTitle(), i), 
                                20, -2.5, 2.5);
                                //100, -2.1, 2.1);
        h1_musPhi[i] = InitTH1F( Form("h1_%s_musPhi_%ifatjet", ch->GetTitle(), i), 
                                Form("h1_%s_musPhi_%ifatjet", ch->GetTitle(), i), 
                                20, -TMath::Pi(), TMath::Pi());
        h1_HT[i] = InitTH1F( Form("h1_%s_HT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_HT_%ifatjet", ch->GetTitle(), i), 
                             //20, 350, 1350);
                             20, 500, 2000);
        h1_MET[i] = InitTH1F( Form("h1_%s_MET_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MET_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             12, 200, 800);
        h1_MET_mTlt140[i] = InitTH1F( Form("h1_%s_MET_mTlt140_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MET_mTlt140_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             25, 100, 600);
        h1_MET_mTgt140[i] = InitTH1F( Form("h1_%s_MET_mTgt140_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_MET_mTgt140_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             25, 100, 600);
        h1_dPhi_mTlt140[i] = InitTH1F( Form("h1_%s_dPhi_mTlt140_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dPhi_mTlt140_%ifatjet", ch->GetTitle(), i), 
                               32, 0, TMath::Pi());
        h1_dPhi_mTgt140[i] = InitTH1F( Form("h1_%s_dPhi_mTgt140_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dPhi_mTgt140_%ifatjet", ch->GetTitle(), i), 
                               32, 0, TMath::Pi());
        h1_METPhi[i] = InitTH1F( Form("h1_%s_METPhi_%ifatjet", ch->GetTitle(), i), 
                                 Form("h1_%s_METPhi_%ifatjet", ch->GetTitle(), i), 
                                 20, -TMath::Pi(), TMath::Pi());
        h1_METx[i] = InitTH1F( Form("h1_%s_METx_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_METx_%ifatjet", ch->GetTitle(), i), 
                               100, -300, 300);
        h1_METy[i] = InitTH1F( Form("h1_%s_METy_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_METy_%ifatjet", ch->GetTitle(), i), 
                               100, -300, 300);
        h1_DPhi[i] = InitTH1F( Form("h1_%s_DPhi_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_DPhi_%ifatjet", ch->GetTitle(), i), 
                               32, 0, TMath::Pi());
        h1_dRFJ[i] = InitTH1F( Form("h1_%s_dRFJ_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dRFJ_%ifatjet", ch->GetTitle(), i), 
                               20, 1, 5);
        h1_dPhiFJ[i] = InitTH1F( Form("h1_%s_dPhiFJ_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dPhiFJ_%ifatjet", ch->GetTitle(), i), 
                               32, 0, TMath::Pi());
        h1_dEtaFJ[i] = InitTH1F( Form("h1_%s_dEtaFJ_%ifatjet", ch->GetTitle(), i), 
                               Form("h1_%s_dEtaFJ_%ifatjet", ch->GetTitle(), i), 
                               20, 0, 4);
        h1_FatjetPt[i] = InitTH1F( Form("h1_%s_FatjetPt_%ifatjet", ch->GetTitle(), i), 
                                   Form("h1_%s_FatjetPt_%ifatjet", ch->GetTitle(), i), 
                                   20, 0, 800);
        h1_FatjetPt_mumatch[i] = InitTH1F( Form("h1_%s_FatjetPt_mumatch_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_FatjetPt_mumatch_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 800);
        h1_FatjetPt_notmumatch[i] = InitTH1F( Form("h1_%s_FatjetPt_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_FatjetPt_notmumatch_%ifatjet", ch->GetTitle(), i), 
                                            20, 0, 800);
        h1_FatjetPt_bmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_bmatch_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_FatjetPt_bmatch_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 800);
        h1_FatjetPt_notbmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_FatjetPt_notbmatch_%ifatjet", ch->GetTitle(), i), 
                                            20, 0, 800);
        h1_FatjetPt1[i] = InitTH1F( Form("h1_%s_FatjetPt1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt1_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 800);
        h1_FatjetPt2[i] = InitTH1F( Form("h1_%s_FatjetPt2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt2_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 500);
        h1_FatjetPt3[i] = InitTH1F( Form("h1_%s_FatjetPt3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt3_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 300);
        h1_FatjetPt4[i] = InitTH1F( Form("h1_%s_FatjetPt4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPt4_%ifatjet", ch->GetTitle(), i), 
                                    20, 0, 150);
        h1_FatjetPhi1[i] = InitTH1F( Form("h1_%s_FatjetPhi1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi1_%ifatjet", ch->GetTitle(), i), 
                                    32, 0, TMath::Pi());
        h1_FatjetPhi2[i] = InitTH1F( Form("h1_%s_FatjetPhi2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi2_%ifatjet", ch->GetTitle(), i), 
                                    32, 0, TMath::Pi());
        h1_FatjetPhi3[i] = InitTH1F( Form("h1_%s_FatjetPhi3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi3_%ifatjet", ch->GetTitle(), i), 
                                    32, 0, TMath::Pi());
        h1_FatjetPhi4[i] = InitTH1F( Form("h1_%s_FatjetPhi4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetPhi4_%ifatjet", ch->GetTitle(), i), 
                                    32, 0, TMath::Pi());
        h1_FatjetEta1[i] = InitTH1F( Form("h1_%s_FatjetEta1_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta1_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_FatjetEta2[i] = InitTH1F( Form("h1_%s_FatjetEta2_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta2_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_FatjetEta3[i] = InitTH1F( Form("h1_%s_FatjetEta3_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta3_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_FatjetEta4[i] = InitTH1F( Form("h1_%s_FatjetEta4_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta4_%ifatjet", ch->GetTitle(), i), 
                                    30, -3, 3);
        h1_BjetPt_bmatchmu[i] = InitTH1F( Form("h1_%s_BjetPt_bmatchmu_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_BjetPt_bmatchmu_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 800);
        h1_BjetPt_notbmatchmu[i] = InitTH1F( Form("h1_%s_BjetPt_notbmatchmu_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_BjetPt_notbmatchmu_%ifatjet", ch->GetTitle(), i), 
                                            20, 0, 800);
        h1_dRbmu_bmatchmu[i] = InitTH1F( Form("h1_%s_dRbmu_bmatchmu_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_dRbmu_bmatchmu_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 1);
        h1_dRbmu_notbmatchmu[i] = InitTH1F( Form("h1_%s_dRbmu_notbmatchmu_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_dRbmu_notbmatchmu_%ifatjet", ch->GetTitle(), i), 
                                            20, 1, 4);
        h1_FatjetPt_metmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_metmatch_%ifatjet", ch->GetTitle(), i), 
                                         Form("h1_%s_FatjetPt_metmatch_%ifatjet", ch->GetTitle(), i), 
                                         20, 0, 800);
        h1_FatjetPt_notmetmatch[i] = InitTH1F( Form("h1_%s_FatjetPt_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                            Form("h1_%s_FatjetPt_notmetmatch_%ifatjet", ch->GetTitle(), i), 
                                            20, 0, 800);
        h1_FatjetEta[i] = InitTH1F( Form("h1_%s_FatjetEta_%ifatjet", ch->GetTitle(), i), 
                                    Form("h1_%s_FatjetEta_%ifatjet", ch->GetTitle(), i), 
                                    20, -2.5, 2.5);
        h1_WpT[i] = InitTH1F( Form("h1_%s_WpT_%ifatjet", ch->GetTitle(), i), 
                             Form("h1_%s_WpT_%ifatjet", ch->GetTitle(), i), 
                             //20, 0, 500);
                             20, 0, 1000);
        h2_mjvsFatjetPt[i]    =   InitTH2F(Form("h2_%s_mjvsFatjetPt_%ifatjet", ch->GetTitle(), i),
                                           Form("h2_%s_mjvsFatjetPt_%ifatjet", ch->GetTitle(), i), 
                                           20, 0, 500, 20, 0, 800);
        h2_mindRFJmjvsmj[i]    =   InitTH2F(Form("h2_%s_mindRFJmjvsmj_%ifatjet", ch->GetTitle(), i),
                                           Form("h2_%s_mindRFJmjvsmj_%ifatjet", ch->GetTitle(), i), 
                                           30, 0, 300, 20, 0, 100);
        h2_MJmT[i]              =   InitTH2F(Form("h2_%s_MJmT_%ifatjet", ch->GetTitle(), i),
                                            Form("h2_%s_MJmT_%ifatjet", ch->GetTitle(), i), 
                                            //2, 0, 800, 2, 0, 300);
                                            2, 0, 2*MJcut, 2, 0, 280);
        for(int iregion=0; iregion<4; iregion++) 
        {
            h2_dPhivsMET[i][iregion]    =   InitTH2F(Form("h2_%s_dPhivsMET_R%i_%ifatjet", ch->GetTitle(), iregion+1, i),
                                            Form("h2_%s_dPhivsMET_R%i_%ifatjet", ch->GetTitle(), iregion+1, i), 
                                            12, 0, 600, 16, 0, TMath::Pi());
            h2_dPhiTTvsMET[i][iregion]    =   InitTH2F(Form("h2_%s_dPhiTTvsMET_R%i_%ifatjet", ch->GetTitle(), iregion+1, i),
                                            Form("h2_%s_dPhiTTvsMET_R%i_%ifatjet", ch->GetTitle(), iregion+1, i), 
                                            12, 0, 600, 16, 0, TMath::Pi());
        }
    }
    // Pile up reweighting hist
    //TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_1p317ipb.root");
    TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_19p5ifb.root");
    TH1F *h1PU = (TH1F*)(fPUFile->Get("puWeights"));
    
    // Btagging probability
    //LoadJetTagEffMaps(ChainName);
    
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
    //InitTree(ch);
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
            if ( isatty(1) )
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

        // Event Filter for data 
        /*
        0 scrapingVeto_decision
        1 hbhefilter_decision
        2 trackingfailurefilter_decision
        3 cschalofilter_decision
        4 eebadscfilter_decision
        5 ecalTPfilter_decision
        6 ecallaserfilter_decision
        7 trackertoomanyclustersfilter_decision
        8 trackertoomanytripletsfilter_decision
        9 trackertoomanyseedsfilter_decision
        10 ecalBEfilter_decision
        11 greedymuonfilter_decision
        12 inconsistentPFmuonfilter_decision
        13 hcallaserfilter_decision
        14 eenoisefilter_decision
        15 trackercoherentnoisefilter1_decision
        16 trackercoherentnoisefilter2_decision
        
        // Paul's code : /cms1r0/pbgeff/RA4_2012/src/RA4Ana.C
        
        if( hasAGoodVertex
        && scrapingVeto_decision
        && hbhefilter_decision
        && trackingfailurefilter_decision
        && hcallaserfilter_decision
        && cschalofilter_decision
        && eebadscfilter_decision
        && ecalTPfilter_decision
        && passesECALHCALfilter
        && passTrigJSON ) ..
        */
        if( ChainName.Contains("DATA") && 
           !(filter->at(0)*filter->at(1)*filter->at(2)*filter->at(13)*filter->at(3)*filter->at(4)*filter->at(5)) ) continue;

        // 
        // weights 
        // 
        // Temp fixes for wrong event weights 

        // Data : trigger
        if(ChainName.Contains("DATA"))  
        { 
            if(ChainName.Contains("MuHad")) EventWeight = EventWeight * (TrigMuon); 
            if(ChainName.Contains("EHad")) EventWeight = EventWeight * (TrigElectron); 
        }
        // Pileup 
        if(!ChainName.Contains("DATA")) 
        {
            // luminostiry correction   
            EventWeight     = EventWeight*1;//19500*19600;                      
            // PU   
            EventWeight     = EventWeight*nPUScaleFactor2012(h1PU, Npu);        
            // Btagging SF (https://twiki.cern.ch/twiki/pub/CMSPublic/TWikiBTV_Moriond2013/Fig17b.pdf) 
            // More precise numbers : https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53X/SFb-pt_payload_Moriond13.txt 
            if(nblow==1) EventWeight     = EventWeight*0.97;        
            if(nblow==2) EventWeight     = EventWeight*0.97*0.97;        
            if(nblow==3) EventWeight     = EventWeight*0.97*0.97*0.97;        
        }
        // TT reweighting
        // Ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
        if(ChainName.Contains("TT")) 
        { 
            // top pT reweighting 
            if(top1pT>400) top1pT=400;
            if(top2pT>400) top2pT=400;
            float weight_top1pT = TMath::Exp(0.159-0.00141*top1pT);
            float weight_top2pT = TMath::Exp(0.159-0.00141*top2pT);
            EventWeight = EventWeight * TMath::Sqrt(weight_top1pT*weight_top2pT);
            EventWeight = EventWeight * 1.01; // normalization correction calculated from 
                                              // Semileptonic sample => 24953451 / 24687562 = 1.01
            // ISR reweighting  
            // description : https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMST2ccMadgraph8TeV
            //float ISRpx = top1pT*TMath::Cos(top1Phi) + top2pT*TMath::Cos(top2Phi); 
            //float ISRpy = top1pT*TMath::Sin(top1Phi) + top2pT*TMath::Sin(top2Phi); 
            //float ISRpT = TMath::Sqrt(ISRpx*ISRpx+ISRpy*ISRpy);

        }
        
        if(ChainName.Contains("T1tttt")) 
        { 

            vector<float>   *GenPt;
            vector<float>   *GenPhi;
            vector<float>   *GenId;
            TBranch         *b_GenPt;   //!
            TBranch         *b_GenPhi;   //!
            TBranch         *b_GenId;   //!
            GenPt  = 0;
            GenPhi  = 0;
            GenId  = 0;
            ch->SetBranchAddress("GenPt",           &GenPt,         &b_GenPt);
            ch->SetBranchAddress("GenPhi",          &GenPhi,        &b_GenPhi);
            ch->SetBranchAddress("GenId",           &GenId,         &b_GenId);

            float ISRpx=0;
            float ISRpy=0;
            
            for(unsigned int igen=0; igen<GenPt->size(); igen++) 
            { 
                if(GenId->at(igen)!=1000021) continue;
                ISRpx = ISRpx + GenPt->at(igen)*TMath::Cos(GenPhi->at(igen)); 
                ISRpy = ISRpy + GenPt->at(igen)*TMath::Sin(GenPhi->at(igen)); 
            }

            float ISRpT = TMath::Sqrt(ISRpx*ISRpx+ISRpy*ISRpy);
            EventWeight = EventWeight * getISRSF(ISRpT);
        }
        
        // Signal cross section
        // Twiki : SUSY xsec twiki : https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVgluglu 
        // mgluino = 1400 
        // 0.000871201 at Twiki 
        // 0.000256511 in the cfA model param
        // mgluino = 1100 0.0101744 
        // mgluino = 1150 0.00666673 
        // mgluino = 1200 0.00440078 
        // mgluino = 1250 0.00291565 
        // mgluino = 1300 0.00194443
        if(ChainName.Contains("f1150_25")) 
            EventWeight = EventWeight*19600*0.00666673/127220;                  
        if(ChainName.Contains("f1200_25")) 
            //EventWeight = EventWeight*19600*0.000871201/137498;                  
            EventWeight = EventWeight*19600*0.00440078/129178;                  
        if(ChainName.Contains("f1200_1000")) 
            //EventWeight = EventWeight*19600*0.000871201/126811;                      
            EventWeight = EventWeight*19600*0.00440078/96418;                      
        if(ChainName.Contains("f1250_25")) 
            EventWeight = EventWeight*19600*0.00291565/131641;                  
        if(ChainName.Contains("f1300_25")) 
            EventWeight = EventWeight*19600*0.00194443/133438;                  
        
        if(ChainName.Contains("f1100_125"))     EventWeight = EventWeight*19600*0.0101744/125152;                  
        if(ChainName.Contains("f1100_225"))     EventWeight = EventWeight*19600*0.0101744/126179;                  
        if(ChainName.Contains("f1100_325"))     EventWeight = EventWeight*19600*0.0101744/125514;                  
        if(ChainName.Contains("f1100_425"))     EventWeight = EventWeight*19600*0.0101744/125663;                  
        if(ChainName.Contains("f1100_525"))     EventWeight = EventWeight*19600*0.0101744/31521;                  
        if(ChainName.Contains("f1100_625"))     EventWeight = EventWeight*19600*0.0101744/41278;                  

        //
        // cuts applied  
        //

        // select only lepton-muon events
        if( (RA4ElsPt->size()+RA4MusPt->size())!=1)                 continue;
        if( ChainName.Contains("MuHad") && RA4MusPt->size()!=1 )    continue;
        if( ChainName.Contains("EHad")  && RA4ElsPt->size()!=1 )    continue;
        //if( RA4MusVetoPt->size()>0 )                                continue;
        //if( RA4ElsVetoPt->size()>0 )                                continue;
        
        // 
        // Btagging probability 
        // 
        /* 
        double Prob0, ProbGEQ1, Prob1, ProbGEQ2,
               Prob2, ProbGEQ3, Prob3, ProbGEQ4; 
        if(!ChainName.Contains("DATA") && Nbcut!=-1)
        {
            CalculateTagProbs(Prob0, ProbGEQ1, Prob1, ProbGEQ2,
                              Prob2, ProbGEQ3, Prob3, ProbGEQ4, 
                              ChainName.Contains("DATA"));
            if(Nbcut==0) EventWeight = EventWeight * Prob0; 
            if(Nbcut==1) EventWeight = EventWeight * Prob1; 
            if(Nbcut==2) EventWeight = EventWeight * Prob2; 
            if(Nbcut==3) EventWeight = EventWeight * Prob3; 
        //cout << Prob0<< " " << Prob1<< " " << Prob2<< " " << Prob3<< " " << endl; // DEBUG
        }
        */ 

        // 
        // Trigger efficiencies
        // applied after one-lepton selection due to muon trigger efficiency 
        // that depends on pT and Eta of muon 
        if(!ChainName.Contains("DATA")) 
        {
            // Trig Efficiency : HT/MET   
            EventWeight     = EventWeight*0.98;     
            // Trig Efficiency : Muon leg    
           if(RA4MusPt->size()==1)  EventWeight = EventWeight * GetMuonEff(/*RA4MusPt->at(0),*/ RA4MusEta->at(0)); 
            // Trig Efficiency : Electron leg    
            if(RA4ElsPt->size()==1)  EventWeight = EventWeight*0.96;                        
        }

        //
        // MET XY shift correction
        //
        /*
        float metx = MET*TMath::Cos(METPhi);
        float mety = MET*TMath::Sin(METPhi);
        if(ChainName.Contains("DATA")) 
        {
            //2012runAvsNvtx_data
            //px = cms.string("+3.54233e-01 + 2.65299e-01*Nvtx"),
            //py = cms.string("+1.88923e-01 - 1.66425e-01*Nvtx")
            metx -= (+3.54233e-01 + 2.65299e-01*Npv);
            mety -= (+1.88923e-01 - 1.66425e-01*Npv);
        } 
        else 
        {
            //2012runAvsNvtx_mc
            //px = cms.string("-2.99576e-02 - 6.61932e-02*Nvtx"),
            //py = cms.string("+3.70819e-01 - 1.48617e-01*Nvtx")
            metx -= (-2.99576e-02 - 6.61932e-02*Npv);
            mety -= (+3.70819e-01 - 1.48617e-01*Npv);
        }
        MET     = TMath::Sqrt(metx*metx+mety*mety);   
        METPhi  = TMath::ATan2(mety,metx);   
        */

        // Nfatjet counting with threshold 
        int Nfatjet_thres = 0;
        for(int ifatjet=0; ifatjet<(int)FatjetPt->size(); ifatjet++)
        {   
            float FatjetPt_cor = FatjetPt->at(ifatjet); 
//            if(RemoveMuon) 
//            {
//                for(int ilep=0; ilep<(int)RA4MusPt->size(); ilep++)  
//                {
//                    float dEta = FatjetEta->at(ifatjet) - RA4MusEta->at(ilep);
//                    float dPhi = FatjetPhi->at(ifatjet) - RA4MusPhi->at(ilep);
//                    if(TMath::Sqrt(dEta*dEta+dPhi*dPhi)<1.0) FatjetPt_cor = FatjetPt_cor - RA4MusPt->at(ilep);
//                }
//            } 
            if(FatjetPt_cor>FatjetpTthres) Nfatjet_thres++;
        }

        //
        // Fill histograms 
        //
        ///// mT
        float mT=-999;
        float dPhi=-999;
        float WpT=-999;
        if(RA4MusPt->size()==1)
        {
            dPhi = getDPhi(METPhi,RA4MusPhi->at(0));
            mT  = TMath::Sqrt( 2*MET*RA4MusPt->at(0)*(1-TMath::Cos(METPhi-RA4MusPhi->at(0))) ); 
        ///// WpT 
            WpT =  TMath::Sqrt(  
                        (RA4MusPt->at(0)*TMath::Cos(RA4MusPhi->at(0)) + MET*TMath::Cos(METPhi))
                       *(RA4MusPt->at(0)*TMath::Cos(RA4MusPhi->at(0)) + MET*TMath::Cos(METPhi))
                       +(RA4MusPt->at(0)*TMath::Sin(RA4MusPhi->at(0)) + MET*TMath::Sin(METPhi))
                       *(RA4MusPt->at(0)*TMath::Sin(RA4MusPhi->at(0)) + MET*TMath::Sin(METPhi))  ); 
        }
        if(RA4ElsPt->size()==1)
        {
            dPhi = getDPhi(METPhi,RA4ElsPhi->at(0));
            mT  = TMath::Sqrt( 2*MET*RA4ElsPt->at(0)*(1-TMath::Cos(METPhi-RA4ElsPhi->at(0))) ); 
        ///// WpT 
            WpT =  TMath::Sqrt(  
                        (RA4ElsPt->at(0)*TMath::Cos(RA4ElsPhi->at(0)) + MET*TMath::Cos(METPhi))
                       *(RA4ElsPt->at(0)*TMath::Cos(RA4ElsPhi->at(0)) + MET*TMath::Cos(METPhi))
                       +(RA4ElsPt->at(0)*TMath::Sin(RA4ElsPhi->at(0)) + MET*TMath::Sin(METPhi))
                       *(RA4ElsPt->at(0)*TMath::Sin(RA4ElsPhi->at(0)) + MET*TMath::Sin(METPhi))  ); 
        }
        ///// MJ 
        float MJ_tmp=0; 
        for(int imj=0; imj<(int)mj->size(); imj++) 
        {
            if(FatjetPt->at(imj)<FatjetpTthres) continue;
            MJ_tmp = MJ_tmp + mj->at(imj);   
        }
        MJ = MJ_tmp;
      
        Nskinnyjet=0;
        HT=0;
        for(int isj=0; isj<(int)JetPt->size(); isj++) 
        { 
            if(JetPt->at(isj)<30) continue;
            Nskinnyjet++;
            HT = HT + JetPt->at(isj);
        }
        
        //if( ChainName.Contains("TT") ) EventWeight = EventWeight * getMJSF(MJ); // MJ correction 
        //if( Nbcut!=-1 && ChainName.Contains("DATA") && NBtagCSVM!=Nbcut) continue;
        
        //
        // selection
        //
        if( HT>500                              && 
            MET>METlow                             && 
            MET<METhigh                             && 
            //RA4MusPt->at(0)>20                  &&  
            //mT>140                             &&  
            //MJ>500                            &&  
            //NBtagCSVM==Nbcut                        &&      // should change Btagging efficiency numbers as well
            NBtagCSVM>=nblow                        &&      // should change Btagging efficiency numbers as well
            NBtagCSVM<=nbhigh                        &&      // should change Btagging efficiency numbers as well
            //NBtagCSVM>1                        &&      // should change Btagging efficiency numbers as well
            Nskinnyjet>=njlow                  &&
            Nskinnyjet<=njhigh 
           ) cut_others = true;

        if(cut_others) {

            int NFJbin = -1;
            if(Nfatjet_thres>3) NFJbin=4;
            else if(Nfatjet_thres>2) NFJbin=3;
            else if(Nfatjet_thres>1) NFJbin=2;
            else if(Nfatjet_thres>0) NFJbin=1;
            else NFJbin=0;
        
            //if(Nfatjet_thres<2) continue;

            if(NFJbin==-1) 
            {   
                cout << "[MJ Analysis] ERROR : NFJ is cannot be negative" << endl;
                continue;
            }
           
            //
            // Fill histogams 
            //
            // dRmin(FJ) 
            float dRFJmin=999.;
            float dPhiFJmin=999.;
            float dEtaFJmin=999.;
            int idRFJ1=0, idRFJ2=0;
            for(int ifj=0; ifj<(int)mj->size(); ifj++) 
            { 
                if(FatjetPt->at(ifj)<FatjetpTthres) continue;
                for(int jfj=0; jfj<ifj; jfj++) 
                { 
                    float dRFJtmp = getDR(FatjetEta->at(jfj), FatjetEta->at(ifj), FatjetPhi->at(jfj), FatjetPhi->at(ifj)); 
                    if(dRFJtmp<dRFJmin)  
                    {   
                        dRFJmin=dRFJtmp;
                        idRFJ1 = jfj;
                        idRFJ2 = ifj;
                    }
                    float dPhiFJtmp = getDPhi(FatjetPhi->at(jfj), FatjetPhi->at(ifj)); 
                    if(dPhiFJtmp<dPhiFJmin)  dPhiFJmin=dPhiFJtmp;
                    float dEtaFJtmp = getDEta(FatjetEta->at(jfj), FatjetEta->at(ifj)); 
                    if(dEtaFJtmp<dEtaFJmin)  dEtaFJmin=dEtaFJtmp;
                }
            }
            
            /// TEST cuts
            //float dRminFJmu=999.;
            //for(int i=0; i<(int)mj->size(); i++) 
            //{   
            //    if(FatjetPt->at(i)<FatjetpTthres) continue;
            //    //if((i==idRFJ1 || i==idRFJ2)) continue;
            //    if(i!=idRFJ2) continue;
            //    float dRtmp =  getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0));
            //    if(dRtmp<dRminFJmu) dRminFJmu = dRtmp;
            //}
            //if(RA4MusPt->at(0)>50 || RA4MusPt->at(0)<30) continue; // FIXME 
            //if(dRFJmin>2) continue; // FIXME 
            //if(METPhi>0) continue; // FIXME 
            //if(dRminFJmu>1.) continue; // FIXME 
            //if(RA4MusPhi->at(0)>0) continue; // FIXME 
            //if(getDPhi(RA4MusPhi->at(0),METPhi)>0.5) continue; // FIXME 
            //if(Npv<20) continue; // FIXME 
            /// FIXME 
            
           ///* 
            if(ChainName.Contains("TT")) 
            {
                float ISRpx = top1pT*TMath::Cos(top1Phi) + top2pT*TMath::Cos(top2Phi); 
                float ISRpy = top1pT*TMath::Sin(top1Phi) + top2pT*TMath::Sin(top2Phi); 
                float ISRpT = TMath::Sqrt(ISRpx*ISRpx+ISRpy*ISRpy);
                //EventWeight = EventWeight * getISRSF(ISRpT);
                EventWeight = EventWeight * 1.013;
                FillTH1F(h1_ttbarpT[NFJbin], ISRpT , EventWeight);                   FillTH1F(h1_ttbarpT[5], ISRpT, EventWeight);
                
                if(MJ<500 && mT<140) { FillTH2F(h2_dPhivsMET[NFJbin][0], MET, dPhi,  EventWeight);  FillTH2F(h2_dPhivsMET[5][0], MET, dPhi,  EventWeight); }
                if(MJ>500 && mT<140) { FillTH2F(h2_dPhivsMET[NFJbin][1], MET, dPhi,  EventWeight);  FillTH2F(h2_dPhivsMET[5][1], MET, dPhi,  EventWeight); }
                if(MJ<500 && mT>140) { FillTH2F(h2_dPhivsMET[NFJbin][2], MET, dPhi,  EventWeight);  FillTH2F(h2_dPhivsMET[5][2], MET, dPhi,  EventWeight); }
                if(MJ>500 && mT>140) { FillTH2F(h2_dPhivsMET[NFJbin][3], MET, dPhi,  EventWeight);  FillTH2F(h2_dPhivsMET[5][3], MET, dPhi,  EventWeight); }
               
                float dPhiTT = getDPhi(top1Phi,top2Phi);   
                if(MJ<500 && mT<140) { FillTH2F(h2_dPhiTTvsMET[NFJbin][0], MET, dPhiTT,  EventWeight);  FillTH2F(h2_dPhiTTvsMET[5][0], MET, dPhiTT,  EventWeight); }
                if(MJ>500 && mT<140) { FillTH2F(h2_dPhiTTvsMET[NFJbin][1], MET, dPhiTT,  EventWeight);  FillTH2F(h2_dPhiTTvsMET[5][1], MET, dPhiTT,  EventWeight); }
                if(MJ<500 && mT>140) { FillTH2F(h2_dPhiTTvsMET[NFJbin][2], MET, dPhiTT,  EventWeight);  FillTH2F(h2_dPhiTTvsMET[5][2], MET, dPhiTT,  EventWeight); }
                if(MJ>500 && mT>140) { FillTH2F(h2_dPhiTTvsMET[NFJbin][3], MET, dPhiTT,  EventWeight);  FillTH2F(h2_dPhiTTvsMET[5][3], MET, dPhiTT,  EventWeight); }
            }
            //*/
            
            // yields
            if(RA4ElsPt->size()==1) {FillTH1F(h1_yields[NFJbin], 0.5, EventWeight);   FillTH1F(h1_yields[5], 0.5, EventWeight);}
            if(RA4MusPt->size()==1) {FillTH1F(h1_yields[NFJbin], 1.5, EventWeight);   FillTH1F(h1_yields[5], 1.5, EventWeight);}

            if(RA4MusPt->size()==1)
            {
                FillTH1F(h1_muspT[NFJbin], RA4MusPt->at(0), EventWeight);   FillTH1F(h1_muspT[5],RA4MusPt->at(0), EventWeight);
                FillTH1F(h1_musEta[NFJbin], RA4MusEta->at(0), EventWeight); FillTH1F(h1_musEta[5], RA4MusEta->at(0), EventWeight);
                FillTH1F(h1_musPhi[NFJbin], RA4MusPhi->at(0), EventWeight); FillTH1F(h1_musPhi[5], RA4MusPhi->at(0), EventWeight);
                if(getDPhi(RA4MusPhi->at(0),METPhi)<0.4) 
                { 
                    FillTH1F(h1_muspTminusMET[NFJbin], (MET-RA4MusPt->at(0))/RA4MusPt->at(0), EventWeight);   
                    FillTH1F(h1_muspTminusMET[5],(MET-RA4MusPt->at(0))/RA4MusPt->at(0), EventWeight);
                }
                FillTH1F(h1_DPhi[NFJbin], getDPhi(RA4MusPhi->at(0),METPhi), EventWeight);   FillTH1F(h1_DPhi[5], getDPhi(RA4MusPhi->at(0),METPhi), EventWeight);
            } 
            if(RA4ElsPt->size()==1)
            {
                FillTH1F(h1_muspT[NFJbin], RA4ElsPt->at(0), EventWeight);   FillTH1F(h1_muspT[5],RA4ElsPt->at(0), EventWeight);
                FillTH1F(h1_musEta[NFJbin], RA4ElsEta->at(0), EventWeight); FillTH1F(h1_musEta[5], RA4ElsEta->at(0), EventWeight);
                FillTH1F(h1_musPhi[NFJbin], RA4ElsPhi->at(0), EventWeight); FillTH1F(h1_musPhi[5], RA4ElsPhi->at(0), EventWeight);
                if(getDPhi(RA4ElsPhi->at(0),METPhi)<0.4) 
                { 
                    FillTH1F(h1_muspTminusMET[NFJbin], (MET-RA4ElsPt->at(0))/RA4ElsPt->at(0), EventWeight);   
                    FillTH1F(h1_muspTminusMET[5],(MET-RA4ElsPt->at(0))/RA4ElsPt->at(0), EventWeight);
                }
                FillTH1F(h1_DPhi[NFJbin], getDPhi(RA4ElsPhi->at(0),METPhi), EventWeight);   FillTH1F(h1_DPhi[5], getDPhi(RA4ElsPhi->at(0),METPhi), EventWeight);
            } 
            FillTH1F(h1_mT[NFJbin], mT, EventWeight);                   FillTH1F(h1_mT[5], mT, EventWeight);
            FillTH1F(h1_WpT[NFJbin], WpT, EventWeight);                 FillTH1F(h1_WpT[5], WpT, EventWeight);
            FillTH1F(h1_Nfatjet[NFJbin], Nfatjet_thres, EventWeight);   FillTH1F(h1_Nfatjet[5], Nfatjet_thres, EventWeight);
            FillTH1F(h1_Nskinnyjet[NFJbin], Nskinnyjet, EventWeight);   FillTH1F(h1_Nskinnyjet[5], Nskinnyjet, EventWeight);
            FillTH1F(h1_MJ[NFJbin], MJ, EventWeight);                   FillTH1F(h1_MJ[5], MJ, EventWeight);
            FillTH1F(h1_HT[NFJbin], HT, EventWeight);                   FillTH1F(h1_HT[5], HT, EventWeight);
            FillTH1F(h1_MET[NFJbin], MET, EventWeight);                 FillTH1F(h1_MET[5], MET, EventWeight);
            if(mT>140) 
            {   
                FillTH1F(h1_MET_mTgt140[NFJbin], MET, EventWeight);                 FillTH1F(h1_MET_mTgt140[5], MET, EventWeight);
                FillTH1F(h1_dPhi_mTgt140[NFJbin], dPhi, EventWeight);               FillTH1F(h1_dPhi_mTgt140[5], dPhi, EventWeight);
            }
            if(mT<140) 
            {
                FillTH1F(h1_MET_mTlt140[NFJbin], MET, EventWeight);                 FillTH1F(h1_MET_mTlt140[5], MET, EventWeight);
                FillTH1F(h1_dPhi_mTlt140[NFJbin], dPhi, EventWeight);               FillTH1F(h1_dPhi_mTlt140[5], dPhi, EventWeight);
            }
            FillTH1F(h1_METPhi[NFJbin], METPhi, EventWeight);           FillTH1F(h1_METPhi[5], METPhi, EventWeight);
            FillTH1F(h1_METx[NFJbin], MET*TMath::Cos(METPhi), EventWeight);           FillTH1F(h1_METx[5], MET*TMath::Cos(METPhi), EventWeight);
            FillTH1F(h1_METy[NFJbin], MET*TMath::Sin(METPhi), EventWeight);           FillTH1F(h1_METy[5], MET*TMath::Sin(METPhi), EventWeight);
            if(Nfatjet_thres>0) 
            {
                FillTH1F(h1_FatjetPt1[NFJbin], FatjetPt->at(0), EventWeight);   FillTH1F(h1_FatjetPt1[5], FatjetPt->at(0), EventWeight); 
                FillTH1F(h1_FatjetPhi1[NFJbin], FatjetPhi->at(0), EventWeight); FillTH1F(h1_FatjetPhi1[5], FatjetPhi->at(0), EventWeight); 
                FillTH1F(h1_FatjetEta1[NFJbin], FatjetEta->at(0), EventWeight); FillTH1F(h1_FatjetEta1[5], FatjetEta->at(0), EventWeight); 
            }
            if(Nfatjet_thres>1) 
            {
                FillTH1F(h1_FatjetPt2[NFJbin], FatjetPt->at(1), EventWeight);   FillTH1F(h1_FatjetPt2[5], FatjetPt->at(1), EventWeight); 
                FillTH1F(h1_FatjetPhi2[NFJbin], FatjetPhi->at(1), EventWeight); FillTH1F(h1_FatjetPhi2[5], FatjetPhi->at(1), EventWeight); 
                FillTH1F(h1_FatjetEta2[NFJbin], FatjetEta->at(1), EventWeight); FillTH1F(h1_FatjetEta2[5], FatjetEta->at(1), EventWeight); 
            }
            if(Nfatjet_thres>2) 
            {
                FillTH1F(h1_FatjetPt3[NFJbin], FatjetPt->at(2), EventWeight);   FillTH1F(h1_FatjetPt3[5], FatjetPt->at(2), EventWeight); 
                FillTH1F(h1_FatjetPhi3[NFJbin], FatjetPhi->at(2), EventWeight); FillTH1F(h1_FatjetPhi3[5], FatjetPhi->at(2), EventWeight); 
                FillTH1F(h1_FatjetEta3[NFJbin], FatjetEta->at(2), EventWeight); FillTH1F(h1_FatjetEta3[5], FatjetEta->at(2), EventWeight); 
            }
            if(Nfatjet_thres>3) 
            {
                FillTH1F(h1_FatjetPt4[NFJbin], FatjetPt->at(3), EventWeight);   FillTH1F(h1_FatjetPt4[5], FatjetPt->at(3), EventWeight); 
                FillTH1F(h1_FatjetPhi4[NFJbin], FatjetPhi->at(3), EventWeight); FillTH1F(h1_FatjetPhi4[5], FatjetPhi->at(3), EventWeight); 
                FillTH1F(h1_FatjetEta4[NFJbin], FatjetEta->at(3), EventWeight); FillTH1F(h1_FatjetEta4[5], FatjetEta->at(3), EventWeight); 
            }
            if(Nfatjet_thres>4) 
            {
                // add 5th jet to FillTH1F(h1_FatjetPt4
                FillTH1F(h1_FatjetPt4[NFJbin], FatjetPt->at(4), EventWeight);   FillTH1F(h1_FatjetPt4[5], FatjetPt->at(4), EventWeight); 
                FillTH1F(h1_FatjetPhi4[NFJbin], FatjetPhi->at(4), EventWeight); FillTH1F(h1_FatjetPhi4[5], FatjetPhi->at(4), EventWeight); 
                FillTH1F(h1_FatjetEta4[NFJbin], FatjetEta->at(4), EventWeight); FillTH1F(h1_FatjetEta4[5], FatjetEta->at(4), EventWeight); 
            }
            if(Nfatjet_thres>5) 
            {
                // add 6th jet to FillTH1F(h1_FatjetPt4
                FillTH1F(h1_FatjetPt4[NFJbin], FatjetPt->at(5), EventWeight);   FillTH1F(h1_FatjetPt4[5], FatjetPt->at(5), EventWeight); 
                FillTH1F(h1_FatjetPhi4[NFJbin], FatjetPhi->at(5), EventWeight); FillTH1F(h1_FatjetPhi4[5], FatjetPhi->at(5), EventWeight); 
                FillTH1F(h1_FatjetEta4[NFJbin], FatjetEta->at(5), EventWeight); FillTH1F(h1_FatjetEta4[5], FatjetEta->at(5), EventWeight); 
            }
            for(unsigned int imj=0; imj<mj->size(); imj++) 
            { 
                if(FatjetPt->at(imj)<FatjetpTthres) continue;

                FillTH1F(h1_mj[NFJbin], mj->at(imj), EventWeight);                            FillTH1F(h1_mj[5], mj->at(imj), EventWeight);
                FillTH1F(h1_FatjetPt[NFJbin], FatjetPt->at(imj), EventWeight);                FillTH1F(h1_FatjetPt[5], FatjetPt->at(imj), EventWeight);
                FillTH1F(h1_FatjetEta[NFJbin], FatjetEta->at(imj), EventWeight);              FillTH1F(h1_FatjetEta[5], FatjetEta->at(imj), EventWeight);
                FillTH2F(h2_mjvsFatjetPt[NFJbin], mj->at(imj), FatjetPt->at(imj), EventWeight); FillTH2F(h2_mjvsFatjetPt[5], mj->at(imj), FatjetPt->at(imj), EventWeight);
                
                /* Turn off matchings  
                // met matching 
                if(getDPhi(FatjetPhi->at(i), METPhi)<1)
                { 
                    FillTH1F(h1_mj_metmatch[NFJbin], mj->at(i), EventWeight);               FillTH1F(h1_mj_metmatch[5], mj->at(i), EventWeight);
                    FillTH1F(h1_FatjetPt_metmatch[NFJbin], FatjetPt->at(i), EventWeight);   FillTH1F(h1_FatjetPt_metmatch[5], FatjetPt->at(i), EventWeight);
                }
                else 
                { 
                    FillTH1F(h1_mj_notmetmatch[NFJbin], mj->at(i), EventWeight);            FillTH1F(h1_mj_notmetmatch[5], mj->at(i), EventWeight);
                    FillTH1F(h1_FatjetPt_notmetmatch[NFJbin], FatjetPt->at(i), EventWeight);FillTH1F(h1_FatjetPt_notmetmatch[5], FatjetPt->at(i), EventWeight);
                }
                // mu matching 
                if(RA4MusPt->size()==1)
                {
                    if(getDR(FatjetEta->at(i), RA4MusEta->at(0), FatjetPhi->at(i), RA4MusPhi->at(0))<1) 
                    { 
                        FillTH1F(h1_mj_mumatch[NFJbin], mj->at(i), EventWeight);                FillTH1F(h1_mj_mumatch[5], mj->at(i), EventWeight);
                        FillTH1F(h1_FatjetPt_mumatch[NFJbin], FatjetPt->at(i), EventWeight);    FillTH1F(h1_FatjetPt_mumatch[5], FatjetPt->at(i), EventWeight);
                    } 
                    else 
                    { 
                        FillTH1F(h1_mj_notmumatch[NFJbin], mj->at(i), EventWeight);             FillTH1F(h1_mj_notmumatch[5], mj->at(i), EventWeight);
                        FillTH1F(h1_FatjetPt_notmumatch[NFJbin], FatjetPt->at(i), EventWeight); FillTH1F(h1_FatjetPt_notmumatch[5], FatjetPt->at(i), EventWeight);
                    }
                }
                // b matching
                float dRmin=999.;
                for(int j=0; j<(int)JetPt->size(); j++) 
                { 
                    if(JetCSV->at(j) < 0.679) continue;
                    float dRtmp = getDR(JetEta->at(j), FatjetEta->at(i), JetPhi->at(j), FatjetPhi->at(i)); 
                    if(dRtmp < dRmin) dRmin = dRtmp; 
                }
                if(dRmin<1.)
                { 
                    FillTH1F(h1_mj_bmatch[NFJbin], mj->at(i), EventWeight);             FillTH1F(h1_mj_bmatch[5], mj->at(i), EventWeight);
                    FillTH1F(h1_FatjetPt_bmatch[NFJbin], FatjetPt->at(i), EventWeight); FillTH1F(h1_FatjetPt_bmatch[5], FatjetPt->at(i), EventWeight);
                } 
                else 
                { 
                    FillTH1F(h1_mj_notbmatch[NFJbin], mj->at(i), EventWeight);              FillTH1F(h1_mj_notbmatch[5], mj->at(i), EventWeight);
                    FillTH1F(h1_FatjetPt_notbmatch[NFJbin], FatjetPt->at(i), EventWeight);  FillTH1F(h1_FatjetPt_notbmatch[5], FatjetPt->at(i), EventWeight);
                }
                */
            }

            FillTH1F(h1_dRFJ[NFJbin], dRFJmin, EventWeight);                                    FillTH1F(h1_dRFJ[5], dRFJmin, EventWeight);
            FillTH1F(h1_dPhiFJ[NFJbin], dPhiFJmin, EventWeight);                                FillTH1F(h1_dPhiFJ[5], dPhiFJmin, EventWeight);
            FillTH1F(h1_dEtaFJ[NFJbin], dEtaFJmin, EventWeight);                                FillTH1F(h1_dEtaFJ[5], dEtaFJmin, EventWeight);
            FillTH2F(h2_mindRFJmjvsmj[NFJbin], mj->at(idRFJ1), mj->at(idRFJ2), EventWeight);    FillTH2F(h2_mindRFJmjvsmj[5], mj->at(idRFJ1), mj->at(idRFJ2), EventWeight);
            FillTH2F(h2_MJmT[NFJbin], MJ, mT, EventWeight);    FillTH2F(h2_MJmT[5], MJ, mT, EventWeight);
            // mj(FJ1+FJ2) 
            float px1 = FatjetPt->at(idRFJ1)*TMath::Cos(FatjetPhi->at(idRFJ1)); 
            float py1 = FatjetPt->at(idRFJ1)*TMath::Sin(FatjetPhi->at(idRFJ1)); 
            float pz1 = FatjetPt->at(idRFJ1)/TMath::Tan(2*TMath::ATan(TMath::Exp(-FatjetEta->at(idRFJ1)))); 
            float e1  = TMath::Sqrt(mj->at(idRFJ1)*mj->at(idRFJ1)+px1*px1+py1*py1+pz1*pz1);  
            float px2 = FatjetPt->at(idRFJ2)*TMath::Cos(FatjetPhi->at(idRFJ2)); 
            float py2 = FatjetPt->at(idRFJ2)*TMath::Sin(FatjetPhi->at(idRFJ2));
            float pz2 = FatjetPt->at(idRFJ2)/TMath::Tan(2*TMath::ATan(TMath::Exp(-FatjetEta->at(idRFJ2)))); 
            float e2  = TMath::Sqrt(mj->at(idRFJ2)*mj->at(idRFJ2)+px2*px2+py2*py2+pz2*pz2);  
            float mjFJFJ = TMath::Sqrt((e1+e2)*(e1+e2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2)); 
            FillTH1F(h1_mjFJFJ[NFJbin], mjFJFJ, EventWeight);   FillTH1F(h1_mjFJFJ[5], mjFJFJ, EventWeight);
            /*
            // bmatching
            int nbclose=0;
            int nbaway=0;
            for(int i=0; i<(int)mj->size(); i++) 
            {   
                if(FatjetPt->at(i)<FatjetpTthres) continue;
                if(i==idRFJ1 || i==idRFJ2) 
                { 
                    for(int j=0; j<(int)JetPt->size(); j++) 
                    { 
                        if(JetCSV->at(j) < 0.679) continue;
                        float dRtmp = getDR(JetEta->at(j), FatjetEta->at(i), JetPhi->at(j), FatjetPhi->at(i)); 
                        if(dRtmp < 1.0) nbclose++; 
                    }
                }
                else 
                {
                    for(int j=0; j<(int)JetPt->size(); j++) 
                    { 
                        if(JetCSV->at(j) < 0.679) continue;
                        float dRtmp = getDR(JetEta->at(j), FatjetEta->at(i), JetPhi->at(j), FatjetPhi->at(i)); 
                        if(dRtmp < 1.0) nbaway++; 
                    }
                } 
            }
            for(int i=0; i<(int)mj->size(); i++) 
            { 
                if(FatjetPt->at(i)<FatjetpTthres) continue;
                if(i==idRFJ1 || i==idRFJ2) continue; 
                if(nbclose<0) continue; // FIXME 
                if(nbaway<0) continue;  // FIXME 
                FillTH1F(h1_mjFJFJnotclose[NFJbin], mj->at(i), EventWeight);                FillTH1F(h1_mjFJFJnotclose[5], mj->at(i), EventWeight);
                FillTH1F(h1_FatjetPtFJFJnotclose[NFJbin], FatjetPt->at(i), EventWeight);    FillTH1F(h1_FatjetPtFJFJnotclose[5], FatjetPt->at(i), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJnotclose[NFJbin], FatjetEta->at(i), EventWeight);  FillTH1F(h1_FatjetEtaFJFJnotclose[5], FatjetEta->at(i), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJnotclose[NFJbin], FatjetPhi->at(i), EventWeight);  FillTH1F(h1_FatjetPhiFJFJnotclose[5], FatjetPhi->at(i), EventWeight);
            }

            if(FatjetPt->at(idRFJ1)>FatjetPt->at(idRFJ2)) 
            {
                FillTH1F(h1_FatjetPtFJFJcloseHigherPt[NFJbin], FatjetPt->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetPtFJFJcloseHigherPt[5], FatjetPt->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseHigherPt[NFJbin], FatjetEta->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseHigherPt[5], FatjetEta->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseHigherPt[NFJbin], FatjetPhi->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseHigherPt[5], FatjetPhi->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetPtFJFJcloseLowerPt[NFJbin], FatjetPt->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetPtFJFJcloseLowerPt[5], FatjetPt->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseLowerPt[NFJbin], FatjetEta->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseLowerPt[5], FatjetEta->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseLowerPt[NFJbin], FatjetPhi->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseLowerPt[5], FatjetPhi->at(idRFJ2), EventWeight);
            }
            else 
            { 
                FillTH1F(h1_FatjetPtFJFJcloseHigherPt[NFJbin], FatjetPt->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetPtFJFJcloseHigherPt[5], FatjetPt->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseHigherPt[NFJbin], FatjetEta->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseHigherPt[5], FatjetEta->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseHigherPt[NFJbin], FatjetPhi->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseHigherPt[5], FatjetPhi->at(idRFJ2), EventWeight);
                FillTH1F(h1_FatjetPtFJFJcloseLowerPt[NFJbin], FatjetPt->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetPtFJFJcloseLowerPt[5], FatjetPt->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseLowerPt[NFJbin], FatjetEta->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetEtaFJFJcloseLowerPt[5], FatjetEta->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseLowerPt[NFJbin], FatjetPhi->at(idRFJ1), EventWeight);
                FillTH1F(h1_FatjetPhiFJFJcloseLowerPt[5], FatjetPhi->at(idRFJ1), EventWeight);
            }
            // b-mu matching
            if(RA4MusPt->size()==1)
            {
                for(int j=0; j<(int)JetPt->size(); j++) 
                { 
                    if(JetCSV->at(j) < 0.679) continue;
                    float dRtmp = getDR(JetEta->at(j), RA4MusEta->at(0), JetPhi->at(j), RA4MusPhi->at(0)); 
                    if(dRtmp<1.0) 
                    {   
                        FillTH1F(h1_BjetPt_bmatchmu[NFJbin], JetPt->at(j), EventWeight);    FillTH1F(h1_BjetPt_bmatchmu[5], JetPt->at(j), EventWeight);
                        FillTH1F(h1_dRbmu_bmatchmu[NFJbin], dRtmp, EventWeight);            FillTH1F(h1_dRbmu_bmatchmu[5], dRtmp, EventWeight);
                    }
                    else 
                    {
                        FillTH1F(h1_BjetPt_notbmatchmu[NFJbin], JetPt->at(j), EventWeight); FillTH1F(h1_BjetPt_notbmatchmu[5], JetPt->at(j), EventWeight);
                        FillTH1F(h1_dRbmu_notbmatchmu[NFJbin], dRtmp, EventWeight);         FillTH1F(h1_dRbmu_notbmatchmu[5], dRtmp, EventWeight);
                    }
                }
            }
                
            // DEBUG     
            if(NFJbin==-1) 
            {
                cout << "Fatjet" << endl;
                for(int i=0; i<(int)mj->size(); i++)
                {
                    if(FatjetPt->at(i)<FatjetpTthres) continue;
                    if(i==idRFJ1 || i==idRFJ2) cout << "Close :: "; 
                    else cout << "Away :: "; 
                    cout << " FatjetPt->at(i): " << FatjetPt->at(i)
                         << " FatjetEta->at(i): " << FatjetEta->at(i)
                         << " FatjetPhi->at(i): " << FatjetPhi->at(i)  << endl;
                } 

                cout << "Bjet" << endl;
                for(int j=0; j<(int)JetPt->size(); j++)
                {
                    if(JetCSV->at(j) < 0.679) continue;
                    cout << " JetPt->at(j): " << JetPt->at(j)
                         << " JetEta->at(j): " << JetEta->at(j)
                         << " JetPhi->at(j): " << JetPhi->at(j)
                         << " JetCSV->at(j): " << JetCSV->at(j) << endl;
                }
                cout << "-------------------------------------------------------------" << endl;
            } // if(NFJbin==-1) 
            */
        } // if(cut_others) 
    } // for(int i = 0; i<nentries; i++)
    
    TString HistFileName = ch->GetTitle();
    HistFileName = Form("HistFiles/%s_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root", HistFileName.Data(),(int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres);
    cout << "[MJ Analysis] Writing " << HistFileName << endl;
    TFile *HistFile = new TFile(HistFileName, "RECREATE");
    gROOT->cd();
    HistFile->cd();
    // write histograms
    for(int i=0; i<6; i++)  
    {
        h1_yields[i]->SetDirectory(0);                      h1_yields[i]->Write();
        h1_MJ[i]->SetDirectory(0);                          h1_MJ[i]->Write();
        h1_HT[i]->SetDirectory(0);                          h1_HT[i]->Write();
        h1_ttbarpT[i]->SetDirectory(0);                     h1_ttbarpT[i]->Write();
        h1_MET[i]->SetDirectory(0);                         h1_MET[i]->Write();
        h1_MET_mTlt140[i]->SetDirectory(0);                 h1_MET_mTlt140[i]->Write();
        h1_MET_mTgt140[i]->SetDirectory(0);                 h1_MET_mTgt140[i]->Write();
        h1_dPhi_mTlt140[i]->SetDirectory(0);                h1_dPhi_mTlt140[i]->Write();
        h1_dPhi_mTgt140[i]->SetDirectory(0);                h1_dPhi_mTgt140[i]->Write();
        h1_METPhi[i]->SetDirectory(0);                      h1_METPhi[i]->Write();
        h1_METx[i]->SetDirectory(0);                        h1_METx[i]->Write();
        h1_METy[i]->SetDirectory(0);                        h1_METy[i]->Write();
        h1_FatjetPt[i]->SetDirectory(0);                    h1_FatjetPt[i]->Write();
        h1_FatjetPt1[i]->SetDirectory(0);                   h1_FatjetPt1[i]->Write();
        h1_FatjetPt2[i]->SetDirectory(0);                   h1_FatjetPt2[i]->Write();
        h1_FatjetPt3[i]->SetDirectory(0);                   h1_FatjetPt3[i]->Write();
        h1_FatjetPt4[i]->SetDirectory(0);                   h1_FatjetPt4[i]->Write();
        h1_FatjetPhi1[i]->SetDirectory(0);                   h1_FatjetPhi1[i]->Write();
        h1_FatjetPhi2[i]->SetDirectory(0);                   h1_FatjetPhi2[i]->Write();
        h1_FatjetPhi3[i]->SetDirectory(0);                   h1_FatjetPhi3[i]->Write();
        h1_FatjetPhi4[i]->SetDirectory(0);                   h1_FatjetPhi4[i]->Write();
        h1_FatjetEta1[i]->SetDirectory(0);                   h1_FatjetEta1[i]->Write();
        h1_FatjetEta2[i]->SetDirectory(0);                   h1_FatjetEta2[i]->Write();
        h1_FatjetEta3[i]->SetDirectory(0);                   h1_FatjetEta3[i]->Write();
        h1_FatjetEta4[i]->SetDirectory(0);                   h1_FatjetEta4[i]->Write();
        h1_FatjetPt_mumatch[i]->SetDirectory(0);            h1_FatjetPt_mumatch[i]->Write();
        h1_FatjetPt_notmumatch[i]->SetDirectory(0);         h1_FatjetPt_notmumatch[i]->Write();
        h1_FatjetPt_bmatch[i]->SetDirectory(0);             h1_FatjetPt_bmatch[i]->Write();
        h1_FatjetPt_notbmatch[i]->SetDirectory(0);          h1_FatjetPt_notbmatch[i]->Write();
        h1_BjetPt_bmatchmu[i]->SetDirectory(0);             h1_BjetPt_bmatchmu[i]->Write();
        h1_BjetPt_notbmatchmu[i]->SetDirectory(0);          h1_BjetPt_notbmatchmu[i]->Write();
        h1_dRbmu_bmatchmu[i]->SetDirectory(0);              h1_dRbmu_bmatchmu[i]->Write();
        h1_dRbmu_notbmatchmu[i]->SetDirectory(0);           h1_dRbmu_notbmatchmu[i]->Write();
        h1_dRFJ[i]->SetDirectory(0);                        h1_dRFJ[i]->Write();
        h1_dPhiFJ[i]->SetDirectory(0);                      h1_dPhiFJ[i]->Write();
        h1_dEtaFJ[i]->SetDirectory(0);                      h1_dEtaFJ[i]->Write();
        h2_mindRFJmjvsmj[i]->SetDirectory(0);               h2_mindRFJmjvsmj[i]->Write();
        h1_FatjetPt_metmatch[i]->SetDirectory(0);           h1_FatjetPt_metmatch[i]->Write();
        h1_FatjetPt_notmetmatch[i]->SetDirectory(0);        h1_FatjetPt_notmetmatch[i]->Write();
        h1_FatjetEta[i]->SetDirectory(0);                   h1_FatjetEta[i]->Write();
        h1_DPhi[i]->SetDirectory(0);                        h1_DPhi[i]->Write();
        h1_mj[i]->SetDirectory(0);                          h1_mj[i]->Write();
        h1_mjFJFJ[i]->SetDirectory(0);                      h1_mjFJFJ[i]->Write();
        h1_mjFJFJnotclose[i]->SetDirectory(0);              h1_mjFJFJnotclose[i]->Write();
        h1_FatjetPtFJFJnotclose[i]->SetDirectory(0);        h1_FatjetPtFJFJnotclose[i]->Write();
        h1_FatjetEtaFJFJnotclose[i]->SetDirectory(0);       h1_FatjetEtaFJFJnotclose[i]->Write();
        h1_FatjetPhiFJFJnotclose[i]->SetDirectory(0);       h1_FatjetPhiFJFJnotclose[i]->Write();
        h1_FatjetPtFJFJcloseHigherPt[i]->SetDirectory(0);          h1_FatjetPtFJFJcloseHigherPt[i]->Write();
        h1_FatjetEtaFJFJcloseHigherPt[i]->SetDirectory(0);          h1_FatjetEtaFJFJcloseHigherPt[i]->Write();
        h1_FatjetPhiFJFJcloseHigherPt[i]->SetDirectory(0);          h1_FatjetPhiFJFJcloseHigherPt[i]->Write();
        h1_FatjetPtFJFJcloseLowerPt[i]->SetDirectory(0);          h1_FatjetPtFJFJcloseLowerPt[i]->Write();
        h1_FatjetEtaFJFJcloseLowerPt[i]->SetDirectory(0);          h1_FatjetEtaFJFJcloseLowerPt[i]->Write();
        h1_FatjetPhiFJFJcloseLowerPt[i]->SetDirectory(0);          h1_FatjetPhiFJFJcloseLowerPt[i]->Write();
        h1_mj_mumatch[i]->SetDirectory(0);                  h1_mj_mumatch[i]->Write();
        h1_mj_notmumatch[i]->SetDirectory(0);               h1_mj_notmumatch[i]->Write();
        h1_mj_notbmatch[i]->SetDirectory(0);                h1_mj_notbmatch[i]->Write();
        h1_mj_bmatch[i]->SetDirectory(0);                   h1_mj_bmatch[i]->Write();
        h1_mj_notmetmatch[i]->SetDirectory(0);              h1_mj_notmetmatch[i]->Write();
        h1_mj_metmatch[i]->SetDirectory(0);                 h1_mj_metmatch[i]->Write();
        h1_mT[i]->SetDirectory(0);                          h1_mT[i]->Write();
        h1_WpT[i]->SetDirectory(0);                         h1_WpT[i]->Write();
        h1_muspT[i]->SetDirectory(0);                       h1_muspT[i]->Write();
        h1_muspTminusMET[i]->SetDirectory(0);               h1_muspTminusMET[i]->Write();
        h1_musEta[i]->SetDirectory(0);                      h1_musEta[i]->Write();
        h1_musPhi[i]->SetDirectory(0);                      h1_musPhi[i]->Write();
        h1_Nfatjet[i]->SetDirectory(0);                     h1_Nfatjet[i]->Write();
        h1_Nskinnyjet[i]->SetDirectory(0);                  h1_Nskinnyjet[i]->Write();
        h2_mjvsFatjetPt[i]->SetDirectory(0);                h2_mjvsFatjetPt[i]->Write();
        h2_MJmT[i]->SetDirectory(0);                        h2_MJmT[i]->Write();
        for(int iregion=0; iregion<4; iregion++) 
        { 
            h2_dPhivsMET[i][iregion]->SetDirectory(0);                        h2_dPhivsMET[i][iregion]->Write(); 
            h2_dPhiTTvsMET[i][iregion]->SetDirectory(0);                      h2_dPhiTTvsMET[i][iregion]->Write(); 
        }
    }
    HistFile->Close();

}

//
// Stacks
//
void DrawStack(TString HistName, int pTR0p5, int NMergeBins/*, int Nbcut=0*/,float METlow, float METhigh, int njlow, int njhigh,int nblow, int nbhigh, float MJcut) 
{ 

    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");

    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root", (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres));
   
    char *var; 
    if(HistName=="MET")                 	var=(char*)"MET [GeV]";
    if(HistName=="METPhi")              	var=(char*)"#phi(MET)";
    if(HistName=="METx")              	    var=(char*)"METx [GeV]";
    if(HistName=="METy")              	    var=(char*)"METy [GeV]";
    if(HistName=="FatjetPt")            	var=(char*)"pT(Fatjet) [GeV]";                   
    if(HistName=="FatjetPt1")            	var=(char*)"pT(Fatjet1) [GeV]";                   
    if(HistName=="FatjetPt2")            	var=(char*)"pT(Fatjet2) [GeV]";                   
    if(HistName=="FatjetPt3")            	var=(char*)"pT(Fatjet3) [GeV]";                   
    if(HistName=="FatjetPt4")            	var=(char*)"pT(Fatjet4) [GeV]";                   
    if(HistName=="FatjetPhi1")            	var=(char*)"#phi(Fatjet1)";                   
    if(HistName=="FatjetPhi2")            	var=(char*)"#phi(Fatjet2)";                   
    if(HistName=="FatjetPhi3")            	var=(char*)"#phi(Fatjet3)";                   
    if(HistName=="FatjetPhi4")            	var=(char*)"#phi(Fatjet4)";                   
    if(HistName=="FatjetEta1")            	var=(char*)"#eta(Fatjet1)";                   
    if(HistName=="FatjetEta2")            	var=(char*)"#eta(Fatjet2)";                   
    if(HistName=="FatjetEta3")            	var=(char*)"#eta(Fatjet3)";                   
    if(HistName=="FatjetEta4")            	var=(char*)"#eta(Fatjet4)";                   
    if(HistName=="FatjetPt_mumatch")    	var=(char*)"pT(Fatjet)(matched muon) [GeV]";    
    if(HistName=="FatjetPt_notmumatch") 	var=(char*)"pT(Fatjet)(not matched muon) [GeV]";
    if(HistName=="FatjetPt_bmatch")     	var=(char*)"pT(Fatjet)(matched b) [GeV]";        
    if(HistName=="FatjetPt_notbmatch")  	var=(char*)"pT(Fatjet)(not matched b) [GeV]";   
    if(HistName=="BjetPt_bmatchmu")     	var=(char*)"pT(bjet)(matched mu) [GeV]";        
    if(HistName=="BjetPt_notbmatchmu")  	var=(char*)"pT(bjet)(not matched mu) [GeV]";   
    if(HistName=="dRbmu_bmatchmu")     	    var=(char*)"#Delta R(mu,b)(matched mu)";        
    if(HistName=="dRbmu_notbmatchmu")  	    var=(char*)"#Delta R(mu,b)(not matched mu)";   
    if(HistName=="FatjetPt_metmatch")   	var=(char*)"pT(Fatjet)(matched met) [GeV]";    
    if(HistName=="FatjetPt_notmetmatch")	var=(char*)"pT(Fatjet)(not matched met) [GeV]"; 
    if(HistName=="FatjetEta")            	var=(char*)"#eta(Fatjet)";                   
    if(HistName=="DPhi")                	var=(char*)"#Delta(MET, muon)";                  
    if(HistName=="dRFJ")                	var=(char*)"min #Delta R(FJ,FJ)";                  
    if(HistName=="dPhiFJ")                	var=(char*)"min #Delta #phi(FJ,FJ)";                  
    if(HistName=="dEtaFJ")                	var=(char*)"min #Delta #eta(FJ,FJ)";                  
    if(HistName=="HT")                  	var=(char*)"H_{T} [GeV]";                        
    if(HistName=="MJ")                  	var=(char*)"M_{J} [GeV]";                        
    if(HistName=="mj")                  	var=(char*)"m_{j} [GeV]";                        
    if(HistName=="mjFJFJ")               	var=(char*)"m_{j}(FJ+FJ) [GeV]";                        
    if(HistName=="mjFJFJnotclose")        	var=(char*)"m_{j}(FJ+FJ not close) [GeV]";                        
    if(HistName=="FatjetPtFJFJnotclose")   	var=(char*)"pT(Fatjet FJ+FJ not close) [GeV]";                        
    if(HistName=="FatjetEtaFJFJnotclose")   var=(char*)"#eta(Fatjet FJ+FJ not close)";                        
    if(HistName=="FatjetPhiFJFJnotclose")   var=(char*)"#phi(Fatjet FJ+FJ not close)";                        
    if(HistName=="FatjetEtaFJFJclose")      var=(char*)"#eta(Fatjet FJ+FJ close)";                        
    if(HistName=="FatjetEtaFJFJcloseHigherPt")  var=(char*)"#eta(Fatjet FJ+FJ close higher pT)";                        
    if(HistName=="FatjetEtaFJFJcloseLowerPt")   var=(char*)"#eta(Fatjet FJ+FJ close lower pT)";                        
    if(HistName=="FatjetPhiFJFJcloseHigherPt")  var=(char*)"#eta(Fatjet FJ+FJ close higher pT)";                        
    if(HistName=="FatjetPhiFJFJcloseLowerPt")   var=(char*)"#eta(Fatjet FJ+FJ close lower pT)";                        
    if(HistName=="FatjetPtFJFJcloseHigherPt")   var=(char*)"#eta(Fatjet FJ+FJ close higher pT)";                        
    if(HistName=="FatjetPtFJFJcloseLowerPt")    var=(char*)"#eta(Fatjet FJ+FJ close lower pT)";                        
    if(HistName=="mj_mumatch")          	var=(char*)"m_{j}(matched muon) [GeV]";         
    if(HistName=="mj_notmumatch")       	var=(char*)"m_{j}(not matched muon) [GeV]";     
    if(HistName=="mj_bmatch")           	var=(char*)"m_{j}(matched b) [GeV]";             
    if(HistName=="mj_notbmatch")        	var=(char*)"m_{j}(not matched b) [GeV]";         
    if(HistName=="mj_metmatch")         	var=(char*)"m_{j}(matched met) [GeV]";           
    if(HistName=="mj_notmetmatch")      	var=(char*)"m_{j}(not matched met) [GeV]";       
    if(HistName=="mT")                  	var=(char*)"m_{T} [GeV]";                        
    if(HistName=="muspT")               	var=(char*)"p_{T}(muon) [GeV]";                  
    if(HistName=="muspTminusMET")          	var=(char*)"(MET-p_{T}(muon))/p_{T}(muon)";                  
    if(HistName=="musEta")              	var=(char*)"#eta(muon)";                         
    if(HistName=="musPhi")              	var=(char*)"#phi(muon)";                         
    if(HistName=="Nfatjet")             	var=(char*)"N_{fatjet}";
    if(HistName=="Nskinnyjet")          	var=(char*)"n_{j}";
    if(HistName=="WpT")                 	var=(char*)"p_{T}(W) [GeV]";

    TH1F *h1_DATA[6], *h1_DATA_MuHad[6], *h1_DATA_EHad[6], 
         *h1_T[6], *h1_TT_sl[6], *h1_TT_ll[6], *h1_TT[6], *h1_WJets[6], *h1_DY[6], *h1_MC[6], *h1_Ratio[6]; 
    TH1F *h1_f1200_25[6], *h1_f1200_1000[6];
    TH1F *h1_f1100_125[6], *h1_f1100_225[6],*h1_f1100_325[6], *h1_f1100_425[6], *h1_f1100_525[6],*h1_f1100_625[6];
    TH1F *h1_One[6];
    THStack *st[6];
    TCanvas *c = new TCanvas("c","c",1200,360);  
    c->Divide(4,1);
    TCanvas *c_one = new TCanvas("c_one","c_one",400,360);  
    for(int i=2; i<6; i++) 
    {
        
        if(i!=5)  continue; // FIXME 

        h1_DATA_MuHad[i]    = (TH1F*)HistFile->Get(Form("h1_DATA_MuHad_%s_%ifatjet", HistName.Data(), i)); 
        h1_DATA_EHad[i]     = (TH1F*)HistFile->Get(Form("h1_DATA_EHad_%s_%ifatjet", HistName.Data(), i)); 
        // merge MuHad + EHad
        h1_DATA[i]    = (TH1F*)h1_DATA_MuHad[i]->Clone("h1_DATA");  h1_DATA[i]->Add(h1_DATA_EHad[i]);
        h1_T[i]         = (TH1F*)HistFile->Get(Form("h1_T_%s_%ifatjet", HistName.Data(), i));
        h1_TT_sl[i]     = (TH1F*)HistFile->Get(Form("h1_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h1_TT_ll[i]     = (TH1F*)HistFile->Get(Form("h1_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h1_WJets[i]     = (TH1F*)HistFile->Get(Form("h1_WJets_%s_%ifatjet", HistName.Data(), i));
        h1_DY[i]      = (TH1F*)HistFile->Get(Form("h1_DY_%s_%ifatjet", HistName.Data(), i)); 
        h1_One[i]       = InitTH1F(Form("h1_One_%i",i), Form("h1_One_%i",i), 1,
                                    h1_DATA[i]->GetXaxis()->GetBinLowEdge(1),
                                    h1_DATA[i]->GetXaxis()->GetBinUpEdge(h1_DATA[i]->GetXaxis()->GetNbins())  );
        h1_One[i]->SetBinContent(1,1);
        h1_f1200_25[i]      = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1200_25_%s_%ifatjet", HistName.Data(), i)); 
        h1_f1200_1000[i]    = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1200_1000_%s_%ifatjet", HistName.Data(), i)); 

        // merge bins
        h1_DATA[i]->Rebin(NMergeBins);
        h1_T[i]->Rebin(NMergeBins);
        h1_TT_sl[i]->Rebin(NMergeBins);
        h1_TT_ll[i]->Rebin(NMergeBins);
        h1_WJets[i]->Rebin(NMergeBins);
        h1_DY[i]->Rebin(NMergeBins);
        h1_f1200_25[i]->Rebin(NMergeBins);
        h1_f1200_1000[i]->Rebin(NMergeBins);
        
        h1_MC[i] = (TH1F*)h1_TT_sl[i]->Clone(Form("h1_MC_%s_%ifatjet", HistName.Data(), i));
        h1_MC[i]->Add(h1_TT_ll[i]);
        h1_MC[i]->Add(h1_WJets[i]);
        h1_MC[i]->Add(h1_T[i]);
        h1_MC[i]->Add(h1_DY[i]);
        h1_Ratio[i] = (TH1F*)h1_DATA[i]->Clone(Form("h1_Ratio_%s_%ifatjet", HistName.Data(), i));
        h1_Ratio[i]->Divide(h1_MC[i]);

        h1_TT[i] = (TH1F*)h1_TT_ll[i]->Clone();
        h1_TT[i]->Add(h1_TT_sl[i]);

        // UCSB RA4 Color scheme  
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

        h1cosmetic(h1_DATA[i],          Form("DATA %ifatjet", i),               kBlack, 1, 0,           var);
        h1cosmetic(h1_TT[i],            Form("TT %ifatjet", i),                 kBlack, 1, kYellow,     var);
        h1cosmetic(h1_TT_sl[i],         Form("TT(l) %ifatjet", i),              kBlack, 1, 1000/*kYellow*/,     var);
        h1cosmetic(h1_TT_ll[i],         Form("TT(ll) %ifatjet", i),             kBlack, 1, 1006/*kOrange-3*/,   var);
        h1cosmetic(h1_T[i],             Form("t+tW %ifatjet", i),               kBlack, 1, 1004/*kGreen+2*/,   var);
        h1cosmetic(h1_WJets[i],         Form("WJets %ifatjet", i),              kBlack, 1, 1001/*kGray+1*/,    var);
        h1cosmetic(h1_DY[i],            Form("DY %ifatjet", i),                 kBlack, 1, 1005/*kRed+1*/,    var);
        h1cosmetic(h1_Ratio[i],         Form(" "),                              kBlack, 1, kBlack,      var);
        h1cosmetic(h1_One[i],           Form(" "),                              kBlue,  1, 0,           var); 
        h1cosmetic(h1_f1200_25[i],      Form("T1tttt(1200,25) %ifatjet", i),    kRed,   1, 0,           var);
        h1cosmetic(h1_f1200_1000[i],    Form("T1tttt(1200,1000) %ifatjet", i),  kBlue,  1, 0,           var);

        //c->cd(i-1);
        c_one->cd();
        TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.04);
        pad1->SetRightMargin(0.1);
        pad1->SetLeftMargin(0.2);
        pad1->SetLogy(0);
        pad1->Draw();
        pad1->cd();
        pad1->cd()->SetLogy(DoLog);
        TString StackTitle = Form("%i fatjets", i);
        if(i==5) StackTitle = "All fatjets";
        if(i==4) StackTitle = "4+ fatjets";
        st[i] = new THStack( Form("Stack %ifatjet", i), StackTitle);
        st[i]->Add(h1_DY[i]);
        st[i]->Add(h1_WJets[i]);
        st[i]->Add(h1_T[i]); 
//         st[i]->Add(h1_TT[i]);
         st[i]->Add(h1_TT_ll[i]);
         st[i]->Add(h1_TT_sl[i]);
        st[i]->SetMaximum(h1_DATA[i]->GetMaximum()*(DoLog?50:1.6));
        //st[i]->SetMinimum(DoLog?0.001:0.1);
        st[i]->Draw("HIST"); 
        //st[i]->GetYaxis()->SetLabelSize(0.07); 
        st[i]->GetXaxis()->SetLabelSize(0); 
        st[i]->GetYaxis()->SetTitleOffset(1.5); 
        //st[i]->GetYaxis()->SetTitleSize(0.07); 
        st[i]->GetYaxis()->SetTitle("Events/bin"); 

        h1_DATA[i]->SetLineColor(kBlack);
        h1_DATA[i]->SetMarkerColor(kBlack);
        h1_DATA[i]->SetMarkerSize(1);
        h1_DATA[i]->SetMarkerStyle(20);
        h1_DATA[i]->SetStats(0);
        if(doData) h1_DATA[i]->Draw("E SAME");

        TLegend *l1 = new TLegend(0.25, 0.6, 0.90, 0.88);
        l1->SetNColumns(2);
        l1->SetBorderSize(0);
        l1->SetFillColor(0);
        l1->SetFillStyle(0);
        l1->SetTextFont(42);
        l1->SetTextAlign(12);
        l1->SetTextSize(0.06);
        l1->SetFillColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->SetShadowColor(kWhite);
        if(doData) l1->AddEntry(h1_DATA[i],        " Data",  "lp");
//        l1->AddEntry(h1_TT[i],          " TT",    "f");
        l1->AddEntry(h1_TT_sl[i],        " ",    "");
        l1->AddEntry(h1_TT_sl[i],        " t#bar{t}(1#font[12]{l})",    "f");
        l1->AddEntry(h1_TT_ll[i],        " t#bar{t}(2#font[12]{l})",    "f");
        l1->AddEntry(h1_T[i],            " t+tW",  "f");
        l1->AddEntry(h1_WJets[i],        " WJets", "f");
        l1->AddEntry(h1_DY[i],           " Drell-Yan",    "f");
//        l1->AddEntry(h1_f1200_25[i],    Form("T1tttt[1200,25]x%i",SignalScale),     "l");
//        l1->AddEntry(h1_f1200_1000[i],  Form("T1tttt[1200,1000]x%i",SignalScale),   "l");
        l1->Draw();
        
        h1_f1200_25[i]->Scale(SignalScale);
        h1_f1200_1000[i]->Scale(SignalScale);
//        h1_f1200_25[i]->Draw("SAME HIST");
//        h1_f1200_1000[i]->Draw("SAME HIST");

        // CMS Labels
        float textSize = 0.04;

        TLatex *TexEnergyLumi = new TLatex(0.9,0.92,Form("#sqrt{s}=8 TeV, L = %.1f fb^{-1}", 19.5));
        TexEnergyLumi->SetNDC();
        TexEnergyLumi->SetTextSize(textSize);
        TexEnergyLumi->SetTextAlign (31);
        TexEnergyLumi->SetLineWidth(2);

        //TLatex *TexCMS = new TLatex(0.2,0.92,"CMS Preliminary");
        //TLatex *TexCMS = new TLatex(0.2,0.92,"CMS work in progress");
        TLatex *TexCMS = new TLatex(0.2,0.92,"n_{lep}=1, H_{T}>500 GeV, MET>200 GeV, n_{j}#geq6, n_{b}#geq2");
        TexCMS->SetNDC();
        TexCMS->SetTextSize(0.035);
        TexCMS->SetTextFont(42);
        TexCMS->SetLineWidth(2);

        // FIXME : need to add lepton flavor
        TString LabelExt = Form("N_{fatjet} = %i", i);
        if(i==4) LabelExt="N_{fatjet} >= 4";
        if(i==5) LabelExt="N_{fatjet} >= 2";
        TLatex *TexExt = new TLatex(0.85,0.7,LabelExt);
        TexExt->SetTextAlign (31);
        TexExt->SetNDC();
        TexExt->SetTextSize(textSize);
        TexExt->SetLineWidth(2);

        TexEnergyLumi->Draw("SAME");
        TexCMS->Draw("SAME");
        //TexExt->Draw("SAME");

        //c->cd(i-1);
        c_one->cd();
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
        h1_Ratio[i]->SetMarkerSize(0.8);
        //h1_Ratio[i]->SetLabelSize(0.1,"XY");
        //h1_Ratio[i]->SetTitleSize(0.1,"XY");
        //h1_Ratio[i]->SetTitleOffset(1.5);
        h1_Ratio[i]->Draw("SAME E X0");
        
        // For David 
        if(HistName=="mj") HistName="JetMass";
        c_one->Print( Form("fig/CompareDataMC_AllFJ_%s_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i%s.pdf", 
                      HistName.Data(), (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5,DoLog?"_log":"") ); 
       

    }

    // print reweighting factor for MJ
    if(HistName=="MJ")   for(int i=1; i<=h1_Ratio[5]->GetXaxis()->GetNbins(); i++) cout << i << " " << h1_Ratio[5]->GetBinContent(i) << " = " << h1_DATA[5]->GetBinContent(i) << " / " << h1_MC[5]->GetBinContent(i) << endl;                        
    // 
    if(HistName=="mj") HistName="JetMass";
    //c->Print( Form("fig/CompareDataMC_%s_pT%i%sNb%i_%s.pdf", 
    //          HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", Nbcut, DoLog?"_log":"") ); 
//    c->Print( Form("fig/CompareDataMC_%s_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i%s.pdf", 
//              HistName.Data(), (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5,DoLog?"_log":"") ); 
    
    // 
    HistFile->Close();
    delete c; 
    delete c_one; 

}

void PrintTable(float METlow, float METhigh, int njlow, int njhigh,int nblow, int nbhigh, float MJcut, int lepflav=0, bool doLatex=false)
{ 
    if(lepflav==0)  cout << "[MJ Table] Yields for Electron+Muon" << endl;
    if(lepflav==11) cout << "[MJ Table] Yields for Electron" << endl;
    if(lepflav==13) cout << "[MJ Table] Yields for Muon" << endl;

    TString HistName="yields";
    int pTR0p5=30;

    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root", (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres)); 
        
    TH1F *h1_DATA[6], *h1_DATA_MuHad[6], *h1_DATA_EHad[6], 
         *h1_T[6], *h1_TT_sl[6], *h1_TT_ll[6], /*h1_TT[6],*/ *h1_WJets[6], *h1_DY[6], *h1_MC[6];
    TH1F *h1_f1150_25[6], *h1_f1200_25[6], *h1_f1250_25[6], *h1_f1300_25[6], *h1_f1200_1000[6];
    TH1F *h1_f1100_125[6], *h1_f1100_225[6],*h1_f1100_325[6], *h1_f1100_425[6], *h1_f1100_525[6],*h1_f1100_625[6];
    for(int i=2; i<6; i++)
    {

        h1_DATA_MuHad[i]    = (TH1F*)HistFile->Get(Form("h1_DATA_MuHad_%s_%ifatjet", HistName.Data(), i)); 
        h1_DATA_EHad[i]     = (TH1F*)HistFile->Get(Form("h1_DATA_EHad_%s_%ifatjet", HistName.Data(), i)); 
        // merge MuHad + EHad
        h1_DATA[i]      = (TH1F*)h1_DATA_MuHad[i]->Clone("h1_DATA");  h1_DATA[i]->Add(h1_DATA_EHad[i]);
        h1_T[i]         = (TH1F*)HistFile->Get(Form("h1_T_%s_%ifatjet", HistName.Data(), i));
        h1_TT_sl[i]     = (TH1F*)HistFile->Get(Form("h1_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h1_TT_ll[i]     = (TH1F*)HistFile->Get(Form("h1_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h1_WJets[i]     = (TH1F*)HistFile->Get(Form("h1_WJets_%s_%ifatjet", HistName.Data(), i));
        h1_DY[i]        = (TH1F*)HistFile->Get(Form("h1_DY_%s_%ifatjet", HistName.Data(), i));
        h1_f1150_25[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1150_25_%s_%ifatjet", HistName.Data(), i));
        h1_f1200_25[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1200_25_%s_%ifatjet", HistName.Data(), i));
        h1_f1250_25[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1250_25_%s_%ifatjet", HistName.Data(), i));
        h1_f1300_25[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1300_25_%s_%ifatjet", HistName.Data(), i));
        h1_f1200_1000[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1200_1000_%s_%ifatjet", HistName.Data(), i));
        h1_f1100_125[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1100_125_%s_%ifatjet", HistName.Data(), i));
        h1_f1100_225[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1100_225_%s_%ifatjet", HistName.Data(), i));
        h1_f1100_325[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1100_325_%s_%ifatjet", HistName.Data(), i));
        h1_f1100_425[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1100_425_%s_%ifatjet", HistName.Data(), i));
        h1_f1100_525[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1100_525_%s_%ifatjet", HistName.Data(), i));
        h1_f1100_625[i] = (TH1F*)HistFile->Get(Form("h1_T1tttt_f1100_625_%s_%ifatjet", HistName.Data(), i));

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
     //   cout << "\\begin{tabular}{c | c c c | c }" << endl;
        cout << "\\begin{tabular}{c | c }" << endl;
        cout << "\\hline \\hline" << endl;
        cout << "Process  & " 
     //        << "$N_{FJ}=2$ \t&" 
     //        << "$N_{FJ}=3$ \t&" 
     //        << "$N_{FJ}\\ge 4$ \t&" 
     //        << "$N_{FJ}\\ge 2$ \t\\\\" 
             << "Yields \t\\\\" 
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
            setw(20) << Form("--------------------") << "-|" << endl;
        cout << "|" << 
            setw(20) << "Process" << " |"  <<
            setw(20) << "N_FJ=2" << " |" <<
            setw(20) << "N_FJ=3" << " |" <<
            setw(20) << "N_FJ>=4" << " |" <<
            setw(20) << "N_FJ>=2" << " |" << endl;
        cout << "|" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "--" <<
            setw(20) << Form("--------------------") << "-|" << endl;
    }

    PrintTableOneLine("TT(sl)",             h1_TT_sl,       lepflav,	doLatex);
    PrintTableOneLine("TT(ll)",             h1_TT_ll,       lepflav,	doLatex);
    PrintTableOneLine("t+tW",               h1_T,           lepflav,	doLatex);
    PrintTableOneLine("WJets",              h1_WJets,       lepflav,	doLatex);
    PrintTableOneLine("DY",                 h1_DY,          lepflav,	doLatex);
    PrintTableOneLine("Total Bkgd",         h1_MC,          lepflav,	doLatex);
    PrintTableOneLine("DATA",               h1_DATA,        lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1200,1000]",  h1_f1200_1000,  lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1150,25]",    h1_f1150_25,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1200,25]",    h1_f1200_25,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1250,25]",    h1_f1250_25,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1300,25]",    h1_f1300_25,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1100,125]",   h1_f1100_125,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1100,225]",   h1_f1100_225,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1100,325]",   h1_f1100_325,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1100,425]",   h1_f1100_425,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1100,525]",   h1_f1100_525,    lepflav,	doLatex);
    PrintTableOneLine("T1tttt[1100,625]",   h1_f1100_625,    lepflav,	doLatex);

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
            setw(20) << Form("--------------------") << "-|" << endl;
    }

    cout << endl;
}

float GetRErr(float Num, float Den, float NumErr, float DenErr)
{ 
    if(Num==0) Num==1; 
    float R = Num/Den;
    float RErr = R*TMath::Sqrt( DenErr*DenErr/Den/Den + NumErr*NumErr/Num/Num ); 
    if(Den==0) return -999; 
    //else if(Num==0) return 0; 
    else return RErr; 
}

float GetXErr(float N1, float N2, float N1Err, float N2Err)
{   
    if(N1==0)  N1=1;
    return TMath::Sqrt( N1*N1*N2Err*N2Err + N2*N2*N1Err*N1Err ); 
}   

float dividefloat(float num, float den)
{ 
    if (den==0.) return -999;
    else if (den==-999) return -999;
    else if (num==-999) return -999;
    else return num/den;
}

const char* PrintCardOneLine(const char* col1, const char* col2, float col3, float col4)
{
    const char* line = Form("%40s%7s%37.3f%37.3f", col1,col2,col3,col4);
    return line;
}

const char* PrintCardOneLine(const char* col1, const char* col2, const char* col3, const char* col4)
{
    const char* line = Form("%40s%7s%40s%40s", col1,col2,col3,col4);
    return line;
}

void GetPrediction(float METlow, float METhigh, int njlow, int njhigh,int nblow, int nbhigh, float MJcut)
{ 
    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root", (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres));
    TH2F *h2_DATA[6], *h2_DATA_MuHad[6], *h2_DATA_EHad[6], 
         *h2_T[6], *h2_TT_sl[6], *h2_TT_ll[6], /*h2_TT[6],*/ *h2_WJets[6], *h2_DY[6], *h2_MC[6];
    TH2F *h2_f1150_25[6], *h2_f1200_25[6], *h2_f1250_25[6], *h2_f1300_25[6], *h2_f1200_1000[6];
    TH2F *h2_f1100_125[6], *h2_f1100_225[6],*h2_f1100_325[6], *h2_f1100_425[6], *h2_f1100_525[6],*h2_f1100_625[6];
    
    for(int i=2; i<6; i++)
    {
        if(i<5) continue; // only inclusive NFJ bin 
        TString HistName="MJmT";

        h2_DATA_MuHad[i]    = (TH2F*)HistFile->Get(Form("h2_DATA_MuHad_%s_%ifatjet", HistName.Data(), i)); 
        h2_DATA_EHad[i]     = (TH2F*)HistFile->Get(Form("h2_DATA_EHad_%s_%ifatjet", HistName.Data(), i)); 
        // merge MuHad + EHad
        h2_DATA[i]      = (TH2F*)h2_DATA_MuHad[i]->Clone("h2_DATA");  h2_DATA[i]->Add(h2_DATA_EHad[i]);
        h2_T[i]         = (TH2F*)HistFile->Get(Form("h2_T_%s_%ifatjet", HistName.Data(), i));
        h2_TT_sl[i]     = (TH2F*)HistFile->Get(Form("h2_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h2_TT_ll[i]     = (TH2F*)HistFile->Get(Form("h2_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h2_WJets[i]     = (TH2F*)HistFile->Get(Form("h2_WJets_%s_%ifatjet", HistName.Data(), i));
        h2_DY[i]        = (TH2F*)HistFile->Get(Form("h2_DY_%s_%ifatjet", HistName.Data(), i));
        h2_f1150_25[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1150_25_%s_%ifatjet", HistName.Data(), i));
        h2_f1200_25[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1200_25_%s_%ifatjet", HistName.Data(), i));
        h2_f1250_25[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1250_25_%s_%ifatjet", HistName.Data(), i));
        h2_f1300_25[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1300_25_%s_%ifatjet", HistName.Data(), i));
        h2_f1200_1000[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1200_1000_%s_%ifatjet", HistName.Data(), i));
        h2_f1100_125[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1100_125_%s_%ifatjet", HistName.Data(), i));
        h2_f1100_225[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1100_225_%s_%ifatjet", HistName.Data(), i));
        h2_f1100_325[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1100_325_%s_%ifatjet", HistName.Data(), i));
        h2_f1100_425[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1100_425_%s_%ifatjet", HistName.Data(), i));
        h2_f1100_525[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1100_525_%s_%ifatjet", HistName.Data(), i));
        h2_f1100_625[i] = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1100_625_%s_%ifatjet", HistName.Data(), i));

        h2_MC[i] = (TH2F*)h2_TT_sl[i]->Clone("h2_MC");
        h2_MC[i]->Add(h2_TT_ll[i]);
        h2_MC[i]->Add(h2_WJets[i]);
        h2_MC[i]->Add(h2_T[i]);
        h2_MC[i]->Add(h2_DY[i]);
    }  

    // 
    // Yield table 
    //

    float N1MC      =  h2_MC[5]->GetBinContent(1,1);  
    float N2MC      =  h2_MC[5]->GetBinContent(2,1); 
    float N3MC      =  h2_MC[5]->GetBinContent(1,2);  
    float N4MC      =  h2_MC[5]->GetBinContent(2,2); 
    float N1MCErr   =  h2_MC[5]->GetBinError(1,1);  
    float N2MCErr   =  h2_MC[5]->GetBinError(2,1); 
    float N3MCErr   =  h2_MC[5]->GetBinError(1,2);  
    float N4MCErr   =  h2_MC[5]->GetBinError(2,2); 
    float N1Data    =  h2_DATA[5]->GetBinContent(1,1);  
    float N2Data    =  h2_DATA[5]->GetBinContent(2,1); 
    float N3Data    =  h2_DATA[5]->GetBinContent(1,2);  
    float N4Data    =  h2_DATA[5]->GetBinContent(2,2);  
    float N1DataErr =  h2_DATA[5]->GetBinError(1,1);  
    float N2DataErr =  h2_DATA[5]->GetBinError(2,1); 
    float N3DataErr =  h2_DATA[5]->GetBinError(1,2);  
    float N4DataErr =  h2_DATA[5]->GetBinError(2,2);  
    
    float N4T1tttt1150      =  h2_f1150_25[5]->GetBinContent(2,2); 
    float N4T1tttt1200      =  h2_f1200_25[5]->GetBinContent(2,2); 
    float N4T1tttt1250      =  h2_f1250_25[5]->GetBinContent(2,2); 
    float N4T1tttt1300      =  h2_f1300_25[5]->GetBinContent(2,2); 
    float N4T1tttt1150Err   =  h2_f1150_25[5]->GetBinError(2,2); 
    float N4T1tttt1200Err   =  h2_f1200_25[5]->GetBinError(2,2); 
    float N4T1tttt1250Err   =  h2_f1250_25[5]->GetBinError(2,2); 
    float N4T1tttt1300Err   =  h2_f1300_25[5]->GetBinError(2,2); 
    
    float N4T1tttt1100_125      =  h2_f1100_125[5]->GetBinContent(2,2); 
    float N4T1tttt1100_225      =  h2_f1100_225[5]->GetBinContent(2,2); 
    float N4T1tttt1100_325      =  h2_f1100_325[5]->GetBinContent(2,2); 
    float N4T1tttt1100_425      =  h2_f1100_425[5]->GetBinContent(2,2); 
    float N4T1tttt1100_525      =  h2_f1100_525[5]->GetBinContent(2,2); 
    float N4T1tttt1100_625      =  h2_f1100_625[5]->GetBinContent(2,2); 
    float N4T1tttt1100_125Err     =  h2_f1100_125[5]->GetBinError(2,2); 
    float N4T1tttt1100_225Err     =  h2_f1100_225[5]->GetBinError(2,2); 
    float N4T1tttt1100_325Err     =  h2_f1100_325[5]->GetBinError(2,2); 
    float N4T1tttt1100_425Err     =  h2_f1100_425[5]->GetBinError(2,2); 
    float N4T1tttt1100_525Err     =  h2_f1100_525[5]->GetBinError(2,2); 
    float N4T1tttt1100_625Err     =  h2_f1100_625[5]->GetBinError(2,2); 

    //  
    // FIXME : From dilepton events
    // 
    // selection is the following
    //      HT>500                              
    //      MET>200
    //      MJ>(<)500
    //      NBtagCSVM>0
    //      Nskinnyjet>4
    //N3MC = 72.67;   N3MCErr = 2.37;
    //N4MC = 5.24;    N4MCErr = 0.64;
    //N3Data = 89;    N3DataErr = TMath::Sqrt(89);
    //N4Data = 7;     N4DataErr = TMath::Sqrt(7);

    cout << endl; 
    cout << "|" <<
        setw(20) << Form(" MC ----------------") << "--" <<
        setw(20) << Form("--------------------") << "-|" << endl;
    cout << "|" << 
        setw(20) << Form("%.2f +/- %.2f", N3MC, N3MCErr) << " |" <<
        setw(20) << Form("%.2f +/- %.2f", N4MC, N4MCErr) << " |" << endl; 
    cout << "|" << 
        setw(20) << Form("%.2f +/- %.2f", N1MC, N1MCErr) << " |" <<
        setw(20) << Form("%.2f +/- %.2f", N2MC, N2MCErr) << " |" << endl; 
    cout << "|" <<
        setw(20) << Form(" Data --------------") << "--" <<
        setw(20) << Form("--------------------") << "-|" << endl;
    cout << "|" << 
        setw(20) << Form("%.2f +/- %.2f", N3Data, N3DataErr) << " |" <<
        setw(20) << Form("%.2f +/- %.2f", N4Data, N4DataErr) << " |" << endl; 
    cout << "|" << 
        setw(20) << Form("%.2f +/- %.2f", N1Data, N1DataErr) << " |" <<
        setw(20) << Form("%.2f +/- %.2f", N2Data, N2DataErr) << " |" << endl; 
    cout << "|" <<
        setw(20) << Form("--------------------") << "--" <<
        setw(20) << Form("--------------------") << "-|" << endl;

    
    //  
    // Definition of regions 
    //  
    // mT  ----------- 
    //     | R3 | R4 |  
    // 140 ----------- 
    //     | R1 | R2 |  
    //   0 ----------- 
    //         400   MJ 
    
    // 
    // Prediction  
    // 
    // R21MC    = NMC2/NMC1  
    // R43MC    = NMC4/NMC3  
    // R21Data  = NData2/NData1  
    // 
    //  Npred = N3Data * R21Data * (R43MC/R21MC)
    //                             ^^^^^^^^^^^^^ 
    //                                kappa 
    // 


    float R21MC      = dividefloat(N2MC,N1MC);  
    float R21MCErr   = GetRErr(N2MC,N1MC,N2MCErr,N1MCErr);  
    float R43MC      = dividefloat(N4MC,N3MC);  
    float R43MCErr   = GetRErr(N4MC,N3MC,N4MCErr,N3MCErr);  
    float R31MC      = dividefloat(N3MC,N1MC);  
    float R31MCErr   = GetRErr(N3MC,N1MC,N3MCErr,N1MCErr);  
    float R42MC      = dividefloat(N4MC,N2MC);  
    float R42MCErr   = GetRErr(N4MC,N2MC,N4MCErr,N2MCErr);  
    float kappaMC      = dividefloat(R43MC,R21MC);  
    float kappaMCErr   = GetRErr(R43MC,R21MC,R43MCErr,R21MCErr);   // kappa
    
    float R21Data    = dividefloat(N2Data,N1Data);  
    float R21DataErr = GetRErr(N2Data,N1Data,N2DataErr,N1DataErr);   // R21DATA
    float R43Data    = dividefloat(N4Data,N3Data);  
    float R43DataErr = GetRErr(N4Data,N3Data,N4DataErr,N3DataErr);   // R43DATA
    float R31Data    = dividefloat(N3Data,N1Data);  
    float R31DataErr = GetRErr(N3Data,N1Data,N3DataErr,N1DataErr);   // R31DATA
    float R42Data    = dividefloat(N4Data,N2Data);  
    float R42DataErr = GetRErr(N4Data,N2Data,N4DataErr,N2DataErr);   // R42DATA
    float kappaData      = dividefloat(R43Data,R21Data);  
    float kappaDataErr   = GetRErr(R43Data,R21Data,R43DataErr,R21DataErr);   // kappa

    float R21DatakappaMC    = R21Data * kappaMC; 
    float R21DatakappaMCErr = GetXErr(R21Data,kappaMC,R21DataErr,kappaMCErr);  
    float Npred             = N3Data * R21DatakappaMC;  
    float NpredErr          = GetXErr(N3Data,R21DatakappaMC,N3DataErr,R21DatakappaMCErr);  

    cout << "R21MC      = " << Form("%.2f +/- %.2f", R21MC,     R21MCErr)   << endl; 
    cout << "R43MC      = " << Form("%.2f +/- %.2f", R43MC,     R43MCErr)   << endl; 
    cout << "kappaMC    = " << Form("%.2f +/- %.2f", kappaMC,     kappaMCErr)   << endl; 
    cout << "kappaData  = " << Form("%.2f +/- %.2f", kappaData,     kappaDataErr)   << endl; 
    cout << "R21Data    = " << Form("%.2f +/- %.2f", R21Data,   R21DataErr) << endl; 
    cout << "N3Data     = " << Form("%.2f +/- %.2f", N3Data,    N3DataErr)  << endl; 
    cout << "------------ " << Form("--------------------")  << endl; 
    cout << "Npred      = " << Form("%.2f +/- %.2f", Npred,     NpredErr)   << endl; 
    cout << "T1tttt[1150,25] = " << Form("%.2f +/- %.2f", N4T1tttt1150,     N4T1tttt1150Err)   << endl; 
    cout << "T1tttt[1200,25] = " << Form("%.2f +/- %.2f", N4T1tttt1200,     N4T1tttt1200Err)   << endl; 
    cout << "T1tttt[1250,25] = " << Form("%.2f +/- %.2f", N4T1tttt1250,     N4T1tttt1250Err)   << endl; 
    cout << "T1tttt[1300,25] = " << Form("%.2f +/- %.2f", N4T1tttt1300,     N4T1tttt1300Err)   << endl; 
    cout << "------------ " << Form("--------------------")  << endl; 

    cout << "------------ " << Form("--------------------")  << endl; 
    cout << Form("Bin %i %i %i %i %i %i %i", (int)METlow, (int)METhigh, (int)njlow, (int)njhigh, (int)nblow, (int)nbhigh, (int)MJcut)
         << " R        = " << setw(20) << "MC" <<  setw(20) << "Data" << endl; 
    cout << Form("[Bin] %i %i %i %i %i %i %i", (int)METlow, (int)METhigh, (int)njlow, (int)njhigh, (int)nblow, (int)nbhigh, (int)MJcut)
         << " R21      = " << setw(20) << Form("%.3f +/- %.3f", R21MC,     R21MCErr) <<  setw(20) << Form("%.3f +/- %.3f", R21Data,   R21DataErr) << endl; 
    cout << Form("[Bin] %i %i %i %i %i %i %i", (int)METlow, (int)METhigh, (int)njlow, (int)njhigh, (int)nblow, (int)nbhigh, (int)MJcut)
         << " R43      = " << setw(20) << Form("%.3f +/- %.3f", R43MC,     R43MCErr) <<  setw(20) << Form("%.3f +/- %.3f", R43Data,   R43DataErr) << endl; 
    cout << Form("[Bin] %i %i %i %i %i %i %i", (int)METlow, (int)METhigh, (int)njlow, (int)njhigh, (int)nblow, (int)nbhigh, (int)MJcut)
         << " R31      = " << setw(20) << Form("%.3f +/- %.3f", R31MC,     R31MCErr) <<  setw(20) << Form("%.3f +/- %.3f", R31Data,   R31DataErr) << endl; 
    cout << Form("[Bin] %i %i %i %i %i %i %i", (int)METlow, (int)METhigh, (int)njlow, (int)njhigh, (int)nblow, (int)nbhigh, (int)MJcut)
         << " R42      = " << setw(20) << Form("%.3f +/- %.3f", R42MC,     R42MCErr) <<  setw(20) << Form("%.3f +/- %.3f", R42Data,   R42DataErr) << endl; 
    cout << Form("[Bin] %i %i %i %i %i %i %i", (int)METlow, (int)METhigh, (int)njlow, (int)njhigh, (int)nblow, (int)nbhigh, (int)MJcut)
         << " kappa    = " << setw(20) << Form("%.3f +/- %.3f", kappaMC,   kappaMCErr) << setw(20) << Form("%.3f +/- %.3f", kappaData,     kappaDataErr) << endl; 
    cout << "------------ " << Form("--------------------")  << endl; 



    cout << "\\begin{table}[!htb]" << endl;
    cout << "\\centering" << endl;
    cout << "\\begin{tabular}{c|c|c|c|c|c}" << endl;
    cout << "\\hline \\hline" << endl;
    //cout << "Process  & Yields   \\\\" << endl;
    cout << "\\multirow{2}{*}{R3} & MC & " << Form("$%.2f \\pm %.2f$", N3MC, N3MCErr) << " & \\multirow{2}{*}{R4} & MC & " << Form("$%.2f \\pm %.2f$", N4MC, N4MCErr) << "\\\\" << endl;
    cout << " & Data & " << Form("$%i$", (int)N3Data) << " & & Data & " << Form("$%i$", (int)N4Data) << "\\\\" << endl;
    cout << "\\hline \\hline" << endl;
    cout << "\\multirow{2}{*}{R1} & MC & " << Form("$%.2f \\pm %.2f$", N1MC, N1MCErr) << " & \\multirow{2}{*}{R2} & MC & " << Form("$%.2f \\pm %.2f$", N2MC, N2MCErr) << "\\\\" << endl;
    cout << " & Data & " << Form("$%i$", (int)N1Data) << " & & Data & " << Form("$%i$", (int)N2Data) << "\\\\" << endl;
    cout << "\\hline \\hline" << endl;
    cout << "\\multicolumn{2}{c|}{$\\kappa_{MC}$} & " << Form("$%.2f \\pm %.2f$", kappaMC, kappaMCErr) << " & \\multicolumn{2}{c|}{T1tttt [1150, 25]} & " << Form("$%.2f \\pm %.2f$", N4T1tttt1150,     N4T1tttt1150Err) << "\\\\" << endl;
    cout << "\\multicolumn{2}{c|}{$\\kappa_{Data}$} & " << (kappaData==-999?"NaN":Form("$%.2f \\pm %.2f$", kappaData, kappaDataErr)) << " & \\multicolumn{2}{c|}{T1tttt [1200, 25]} & " << Form("$%.2f \\pm %.2f$", N4T1tttt1200,     N4T1tttt1200Err)  << "\\\\" << endl;
    cout << Form("\\multicolumn{2}{c|}{$N_{pred}$} & $%.2f\\pm%.2f$",Npred,NpredErr) << " & \\multicolumn{2}{c|}{T1tttt [1250, 25]} & " << Form("$%.2f \\pm %.2f$", N4T1tttt1250,     N4T1tttt1250Err) << "\\\\" << endl;
    cout << "\\hline \\hline" << endl;
    cout << "\\end{tabular}" << endl;
    cout << Form("\\caption{$%i<MET<%i$ GeV and $%i\\le n_{j} \\le %i$}", (int)METlow, (int)METhigh, (int)njlow, (int)njhigh) << endl;
    cout << "\\end{table}" << endl;


    // -------------------------------------
    // Print out on file
    // -------------------------------------
    TString Region = Form("MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i", (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut);
    
    // 
    // T1tttt 1150_25 
    // 
    ofstream fout_T1tttt_1150_25;
    fout_T1tttt_1150_25.open(Form("Cards/datacard_T1tttt_1150_25_%s.dat", Region.Data()));

    fout_T1tttt_1150_25 << "imax 1 number of channels" << endl;
    fout_T1tttt_1150_25 << "jmax * number of background" << endl;
    fout_T1tttt_1150_25 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1150_25 << "Observation " << N4Data << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine("bin",                         "",     Region,         Region                      ) << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine("process",                     "",     "T1tttt",       "bkg"                       ) << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine("process",                     "",     "0",            "1"                         ) << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1150),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1150_25 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",     "-" ) << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine("T1tttt_1150_25",           "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1150_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    "lnN",  "-",     Form("%.3f",1+NpredErr/Npred) ) << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1150_25 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;

    fout_T1tttt_1150_25.close();

    // 
    // T1tttt 1200_25 
    // 
    ofstream fout_T1tttt_1200_25;
    fout_T1tttt_1200_25.open(Form("Cards/datacard_T1tttt_1200_25_%s.dat", Region.Data()));

    fout_T1tttt_1200_25 << "imax 1 number of channels" << endl;
    fout_T1tttt_1200_25 << "jmax * number of background" << endl;
    fout_T1tttt_1200_25 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1200_25 << "Observation " << N4Data << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine("bin",                         "",     Region,         Region                      ) << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine("process",                     "",     "T1tttt",       "bkg"                       ) << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine("process",                     "",     "0",            "1"                         ) << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1200),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1200_25 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",     "-" ) << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine("T1tttt_1200_25",           "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1200_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    "lnN",  "-",     Form("%.3f",1+NpredErr/Npred) ) << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1200_25 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;

    fout_T1tttt_1200_25.close();

    // 
    // T1tttt 1250_25 
    // 
    ofstream fout_T1tttt_1250_25;
    fout_T1tttt_1250_25.open(Form("Cards/datacard_T1tttt_1250_25_%s.dat", Region.Data()));

    fout_T1tttt_1250_25 << "imax 1 number of channels" << endl;
    fout_T1tttt_1250_25 << "jmax * number of background" << endl;
    fout_T1tttt_1250_25 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1250_25 << "Observation " << N4Data << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine("bin",                         "",     Region,         Region                      ) << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine("process",                     "",     "T1tttt",       "bkg"                       ) << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine("process",                     "",     "0",            "1"                         ) << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1250),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1250_25 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",     "-" ) << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine("T1tttt_1250_25",           "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1250_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    "lnN",  "-",     Form("%.3f",1+NpredErr/Npred) ) << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1250_25 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;

    fout_T1tttt_1250_25.close();

    // 
    // T1tttt 1300_25 
    // 
    ofstream fout_T1tttt_1300_25;
    fout_T1tttt_1300_25.open(Form("Cards/datacard_T1tttt_1300_25_%s.dat", Region.Data()));

    fout_T1tttt_1300_25 << "imax 1 number of channels" << endl;
    fout_T1tttt_1300_25 << "jmax * number of background" << endl;
    fout_T1tttt_1300_25 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1300_25 << "Observation " << N4Data << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine("bin",                                 "",     Region,                      Region                          ) << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine("process",                             "",     "T1tttt",                    "bkg"                           ) << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine("process",                             "",     "0",                         "1"                             ) << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1300),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1300_25 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",                       "-"                             ) << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine("T1tttt_1300_25",                   "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1300_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),     "lnN",  "-",                         Form("%.3f",1+NpredErr/Npred)   ) << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1300_25 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;

    fout_T1tttt_1300_25.close();

    // 
    // T1tttt 1100_125 
    // 
    ofstream fout_T1tttt_1100_125;
    fout_T1tttt_1100_125.open(Form("Cards/datacard_T1tttt_1100_125_%s.dat", Region.Data()));

    fout_T1tttt_1100_125 << "imax 1 number of channels" << endl;
    fout_T1tttt_1100_125 << "jmax * number of background" << endl;
    fout_T1tttt_1100_125 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1100_125 << "Observation " << N4Data << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine("bin",                                 "",     Region,                      Region                          ) << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine("process",                             "",     "T1tttt",                    "bkg"                           ) << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine("process",                             "",     "0",                         "1"                             ) << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1100_125),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1100_125 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",                       "-"                             ) << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine("T1tttt_1100_125",                   "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1100_125 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),     "lnN",  "-",                         Form("%.3f",1+NpredErr/Npred)   ) << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1100_125 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;

    fout_T1tttt_1100_125.close();

    // 
    // T1tttt 1100_225 
    // 
    ofstream fout_T1tttt_1100_225;
    fout_T1tttt_1100_225.open(Form("Cards/datacard_T1tttt_1100_225_%s.dat", Region.Data()));

    fout_T1tttt_1100_225 << "imax 1 number of channels" << endl;
    fout_T1tttt_1100_225 << "jmax * number of background" << endl;
    fout_T1tttt_1100_225 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1100_225 << "Observation " << N4Data << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine("bin",                                 "",     Region,                      Region                          ) << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine("process",                             "",     "T1tttt",                    "bkg"                           ) << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine("process",                             "",     "0",                         "1"                             ) << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1100_225),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1100_225 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",                       "-"                             ) << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine("T1tttt_1100_225",                   "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1100_225 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),     "lnN",  "-",                         Form("%.3f",1+NpredErr/Npred)   ) << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1100_225 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;
    
    fout_T1tttt_1100_225.close();
    
    // 
    // T1tttt 1100_325 
    // 
    ofstream fout_T1tttt_1100_325;
    fout_T1tttt_1100_325.open(Form("Cards/datacard_T1tttt_1100_325_%s.dat", Region.Data()));

    fout_T1tttt_1100_325 << "imax 1 number of channels" << endl;
    fout_T1tttt_1100_325 << "jmax * number of background" << endl;
    fout_T1tttt_1100_325 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1100_325 << "Observation " << N4Data << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine("bin",                                 "",     Region,                      Region                          ) << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine("process",                             "",     "T1tttt",                    "bkg"                           ) << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine("process",                             "",     "0",                         "1"                             ) << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1100_325),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1100_325 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",                       "-"                             ) << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine("T1tttt_1100_325",                   "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1100_325 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),     "lnN",  "-",                         Form("%.3f",1+NpredErr/Npred)   ) << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1100_325 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;
    
    fout_T1tttt_1100_325.close();

    // 
    // T1tttt 1100_425 
    // 
    ofstream fout_T1tttt_1100_425;
    fout_T1tttt_1100_425.open(Form("Cards/datacard_T1tttt_1100_425_%s.dat", Region.Data()));

    fout_T1tttt_1100_425 << "imax 1 number of channels" << endl;
    fout_T1tttt_1100_425 << "jmax * number of background" << endl;
    fout_T1tttt_1100_425 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1100_425 << "Observation " << N4Data << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine("bin",                                 "",     Region,                      Region                          ) << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine("process",                             "",     "T1tttt",                    "bkg"                           ) << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine("process",                             "",     "0",                         "1"                             ) << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1100_425),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1100_425 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",                       "-"                             ) << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine("T1tttt_1100_425",                   "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1100_425 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),     "lnN",  "-",                         Form("%.3f",1+NpredErr/Npred)   ) << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1100_425 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;
    
    fout_T1tttt_1100_425.close();

    // 
    // T1tttt 1100_525 
    // 
    ofstream fout_T1tttt_1100_525;
    fout_T1tttt_1100_525.open(Form("Cards/datacard_T1tttt_1100_525_%s.dat", Region.Data()));

    fout_T1tttt_1100_525 << "imax 1 number of channels" << endl;
    fout_T1tttt_1100_525 << "jmax * number of background" << endl;
    fout_T1tttt_1100_525 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1100_525 << "Observation " << N4Data << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine("bin",                                 "",     Region,                      Region                          ) << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine("process",                             "",     "T1tttt",                    "bkg"                           ) << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine("process",                             "",     "0",                         "1"                             ) << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1100_525),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1100_525 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",                       "-"                             ) << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine("T1tttt_1100_525",                   "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1100_525 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),     "lnN",  "-",                         Form("%.3f",1+NpredErr/Npred)   ) << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1100_525 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;
    
    fout_T1tttt_1100_525.close();

    // 
    // T1tttt 1100_625 
    // 
    ofstream fout_T1tttt_1100_625;
    fout_T1tttt_1100_625.open(Form("Cards/datacard_T1tttt_1100_625_%s.dat", Region.Data()));

    fout_T1tttt_1100_625 << "imax 1 number of channels" << endl;
    fout_T1tttt_1100_625 << "jmax * number of background" << endl;
    fout_T1tttt_1100_625 << "kmax * number of nuisance parameters" << endl;
    fout_T1tttt_1100_625 << "Observation " << N4Data << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine("bin",                                 "",     Region,                      Region                          ) << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine("process",                             "",     "T1tttt",                    "bkg"                           ) << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine("process",                             "",     "0",                         "1"                             ) << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine("rate",                                "",     Form("%.3f",N4T1tttt1100_625),    Form("%.3f",Npred)              ) << endl;
    //fout_T1tttt_1100_625 << PrintCardOneLine(Form("T1tttt_%s",Region.Data()),    "lnN",  "1.3",                       "-"                             ) << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine("T1tttt_1100_625",                   "lnN",  "1.300",     "-" ) << endl;
    //fout_T1tttt_1100_625 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),     "lnN",  "-",                         Form("%.3f",1+NpredErr/Npred)   ) << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine(Form("ttbar_SF_%s",Region.Data()),    Form("gmN %i", (int)N3Data),  "-",     Form("%.3f", R21DatakappaMC) ) << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine(Form("ttbar_uncert_%s",Region.Data()),    "lnN",  "-",     Form("%.3f", (1+R21DatakappaMCErr/R21DatakappaMC)) ) << endl;
    fout_T1tttt_1100_625 << PrintCardOneLine(Form("datacloser_%s",Region.Data()),    "lnN",  "-",     "1.380") << endl;
    
    fout_T1tttt_1100_625.close();


}

//
// 2D histograms 
//
void Draw2D(TString HistName, int pTR0p5/*, int NMergeBins=1*/, bool doRatio/*, int Nbcut=0*/,float METlow, float METhigh, int njlow, int njhigh,int nblow, int nbhigh, float MJcut) 
{ 
    TFile* HistFile = TFile::Open(Form("HistFiles/Hist_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root",(int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres));
   
    char *Xvar, *Yvar, *Zvar; 
    if(HistName=="mjvsmj")
    {   
        Xvar=(char*)"mj [GeV]";
        Yvar=(char*)"mj [GeV]";
        Zvar=(char*)"Events/bin";
    }
    if(HistName=="mjvsFatjetPt")
    {   
        Xvar=(char*)"mj [GeV]";
        Yvar=(char*)"pT(fatjet) [GeV]";
        Zvar=(char*)"Events/bin";
    }


    TH2F *h2_DATA[6], *h2_DATA_MuHad[6], *h2_DATA_EHad[6], *h2_T[6], *h2_TT_sl[6], *h2_TT_ll[6], /*h2_TT[6],*/ *h2_WJets[6], *h2_DY[6]; 
    TH2F *h2_Ratio[6], *h2_MC[6]; 
    TH2F *h2_f1200_25[6], *h2_f1200_1000[6];
    TCanvas *c_TT_sl = new TCanvas("c_TT_sl","c_TT_sl",1600,360);  
    TCanvas *c_TT_ll = new TCanvas("c_TT_ll","c_TT_ll",1600,360);  
    TCanvas *c_Ratio = new TCanvas("c_Ratio","c_Ratio",1600,360);  
    c_TT_sl->Divide(4,1);
    c_TT_ll->Divide(4,1);
    c_Ratio->Divide(4,1);
    for(int i=2; i<6; i++) 
    {
        if(HistName=="mjvsmj" && i!=2) continue;

        h2_DATA_MuHad[i]    = (TH2F*)HistFile->Get(Form("h2_DATA_MuHad_%s_%ifatjet", HistName.Data(), i)); 
        h2_DATA_EHad[i]     = (TH2F*)HistFile->Get(Form("h2_DATA_EHad_%s_%ifatjet", HistName.Data(), i)); 
        // merge MuHad + EHad
        h2_DATA[i]          = (TH2F*)h2_DATA_MuHad[i]->Clone("h2_DATA");  h2_DATA[i]->Add(h2_DATA_EHad[i]);
        h2_T[i]             = (TH2F*)HistFile->Get(Form("h2_T_%s_%ifatjet", HistName.Data(), i));
        h2_TT_sl[i]         = (TH2F*)HistFile->Get(Form("h2_TT_sl_%s_%ifatjet", HistName.Data(), i));
        h2_TT_ll[i]         = (TH2F*)HistFile->Get(Form("h2_TT_ll_%s_%ifatjet", HistName.Data(), i));
        h2_WJets[i]         = (TH2F*)HistFile->Get(Form("h2_WJets_%s_%ifatjet", HistName.Data(), i));
        h2_DY[i]          = (TH2F*)HistFile->Get(Form("h2_DY_%s_%ifatjet", HistName.Data(), i)); 
        h2_f1200_25[i]      = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1200_25_%s_%ifatjet", HistName.Data(), i)); 
        h2_f1200_1000[i]    = (TH2F*)HistFile->Get(Form("h2_T1tttt_f1200_1000_%s_%ifatjet", HistName.Data(), i)); 

        h2cosmetic(h2_DATA[i],          Form("DATA %ifatjet", i),               Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_TT_sl[i],         Form("TT(sl) %ifatjet", i),             Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_TT_ll[i],         Form("TT(ll) %ifatjet", i),             Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_T[i],             Form("T %ifatjet", i),                  Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_WJets[i],         Form("WJets %ifatjet", i),              Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_f1200_25[i],      Form("T1tttt(1200,25) %ifatjet", i),    Xvar,   Yvar,   Zvar);
        h2cosmetic(h2_f1200_1000[i],    Form("T1tttt(1200,1000) %ifatjet", i),  Xvar,   Yvar,   Zvar);

        // MC 
        h2_MC[i] = (TH2F*)h2_TT_sl[i]->Clone(Form("h2_MC_%s_%ifatjet", HistName.Data(), i));
        h2_MC[i]->Add(h2_TT_ll[i]);
        h2_MC[i]->Add(h2_WJets[i]);
        h2_MC[i]->Add(h2_T[i]);
        h2_MC[i]->Add(h2_DY[i]);
    }

    if(HistName=="MJmT") 
    { 
        cout << "-- DATA " << endl; 
        cout << h2_DATA[5]->GetBinContent(1,2) << " " << h2_DATA[5]->GetBinContent(2,2) << endl; 
        cout << h2_DATA[5]->GetBinContent(1,1) << " " << h2_DATA[5]->GetBinContent(2,1) << endl; 
        cout << "-- TT sl " << endl; 
        cout << h2_TT_sl[5]->GetBinContent(1,2) << " " << h2_TT_sl[5]->GetBinContent(2,2) << endl; 
        cout << h2_TT_sl[5]->GetBinContent(1,1) << " " << h2_TT_sl[5]->GetBinContent(2,1) << endl; 
        cout << "-- TT_ll " << endl; 
        cout << h2_TT_ll[5]->GetBinContent(1,2) << " " << h2_TT_ll[5]->GetBinContent(2,2) << endl; 
        cout << h2_TT_ll[5]->GetBinContent(1,1) << " " << h2_TT_ll[5]->GetBinContent(2,1) << endl; 
        cout << "-- T " << endl; 
        cout << h2_T[5]->GetBinContent(1,2) << " " << h2_T[5]->GetBinContent(2,2) << endl; 
        cout << h2_T[5]->GetBinContent(1,1) << " " << h2_T[5]->GetBinContent(2,1) << endl; 
        cout << "-- MC " << endl; 
        cout << h2_MC[5]->GetBinContent(1,2) << " " << h2_MC[5]->GetBinContent(2,2) << endl; 
        cout << h2_MC[5]->GetBinContent(1,1) << " " << h2_MC[5]->GetBinContent(2,1) << endl; 
        cout << "-- f1200_25 " << endl; 
        cout << h2_f1200_25[5]->GetBinContent(1,2) << " " << h2_f1200_25[5]->GetBinContent(2,2) << endl; 
        cout << h2_f1200_25[5]->GetBinContent(1,1) << " " << h2_f1200_25[5]->GetBinContent(2,1) << endl; 
    }

    // Drawing
    for(int i=2; i<6; i++) 
    {
        if(HistName=="mjvsmj" && i!=2) continue;
        c_TT_sl->cd(i-1);
        h2_TT_sl[i]->Draw("colz");
        c_TT_ll->cd(i-1);
        h2_TT_ll[i]->Draw("colz");
    
        //c_TT_sl->Print( Form("fig/2D_TT_sl_%s_pT%i%sNb%i_%s.pdf",
        //            HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", Nbcut, DoLog?"_log":"") ); 
        //c_TT_ll->Print( Form("fig/2D_TT_ll_%s_pT%i%sNb%i_%s.pdf", 
        //            HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", Nbcut, DoLog?"_log":"") ); 
        //c_TT_sl->Print( Form("fig/2D_TT_sl_%s_pT%i%s_%s.pdf",
        //            HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", DoLog?"_log":"") ); 
        //c_TT_ll->Print( Form("fig/2D_TT_ll_%s_pT%i%s_%s.pdf", 
        //            HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", DoLog?"_log":"") ); 
    }
  
    if(doRatio) 
    { 
        for(int i=2; i<6; i++) 
        {
            h2_Ratio[i] = (TH2F*) h2_DATA[i]->Clone(); 
            h2_MC[i]    = (TH2F*) h2_TT_sl[i]->Clone();
            h2_MC[i]->Add(h2_TT_ll[i]);
            h2_MC[i]->Add(h2_T[i]);
            h2_MC[i]->Add(h2_WJets[i]);
            h2_Ratio[i]->Divide(h2_MC[i]);
            c_Ratio->cd(i-1);
            h2_Ratio[i]->SetMinimum(0.5);
            h2_Ratio[i]->SetMaximum(1.5);
            h2_Ratio[i]->Draw("colz");
            //c_Ratio->Print( Form("fig/2D_Ratio_%s_pT%i%sNb%i_%s.pdf", 
            //            HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", Nbcut, DoLog?"_log":"") ); 
            //c_Ratio->Print( Form("fig/2D_Ratio_%s_pT%i%s_%s.pdf", 
            //            HistName.Data(), pTR0p5, RemoveMuon?"_MuonRemoved":"", DoLog?"_log":"") ); 
        } 
    }

    // 
    HistFile->Close();
    delete c_TT_sl;
    delete c_TT_ll;
    delete c_Ratio;
}


//
// main 
//
void DrawSL(float METlow, float METhigh, int njlow, int njhigh,int nblow, int nbhigh, float MJcut, bool OnlyDraw=false/*, int Nbcut=0*/) 
{

    gROOT->ProcessLine(".L /Users/jaehyeok/macros/rootlogon.C");
    
    // ----------------------------------------
    //  Define chains  
    // ----------------------------------------
    TChain *ch_data_muhad   = new TChain("tree", "DATA_MuHad");
    TChain *ch_data_ehad    = new TChain("tree", "DATA_EHad");
    TChain *ch_ttbar_sl     = new TChain("tree", "TT_sl");
    TChain *ch_ttbar_ll     = new TChain("tree", "TT_ll");
    TChain *ch_wjets        = new TChain("tree", "WJets");
    TChain *ch_dy           = new TChain("tree", "DY");
    TChain *ch_t            = new TChain("tree", "T");
    TChain *ch_f1150_25     = new TChain("tree", "T1tttt_f1150_25");
    TChain *ch_f1200_25     = new TChain("tree", "T1tttt_f1200_25");
    TChain *ch_f1250_25     = new TChain("tree", "T1tttt_f1250_25");
    TChain *ch_f1300_25     = new TChain("tree", "T1tttt_f1300_25");
    TChain *ch_f1200_1000   = new TChain("tree", "T1tttt_f1200_1000");
    TChain *ch_f1100_125    = new TChain("tree", "T1tttt_f1100_125");
    TChain *ch_f1100_225    = new TChain("tree", "T1tttt_f1100_225");
    TChain *ch_f1100_325    = new TChain("tree", "T1tttt_f1100_325");
    TChain *ch_f1100_425    = new TChain("tree", "T1tttt_f1100_425");
    TChain *ch_f1100_525    = new TChain("tree", "T1tttt_f1100_525");
    TChain *ch_f1100_625    = new TChain("tree", "T1tttt_f1100_625");
   
    //TString BabyDir = "../../../babies/HT500MET100_17Jan2014/HT500MET100_OR_LL/";  
    TString BabyDir = "../../../babies/HT500MET100/";  
    // Data
    //ch_data->Add(BabyDir+"baby_SingleMu_Run2012A-13Jul2012-v1_f*.root");  // 808 pb if fully processed : 808 * 646/650 = 803 pb 
    //ch_data->Add(BabyDir+"baby_SingleMu_Run2012B-13Jul2012-v1_f*.root");  // 5237-808=4429 pb if fully processed : 4429 * 4247/4294 = 4380 pb  
                                                                            // 803 + 4380 = 5183
    //ch_data->Add(BabyDir+"baby_SingleMu_Run2012*.root");                  // Full SingleMu (19.1 fb-1)
    ch_data_muhad->Add(BabyDir+"baby_MuHad_*.root");                        // Full Muhad 
    ch_data_ehad->Add(BabyDir+"baby_HT_*.root");                            // Full Electronhad 
    ch_data_ehad->Add(BabyDir+"baby_ElectronHad_*.root");                  
    
    // TT 
    //ch_ttbar_sl->Add(BabyDir+"baby_TT_CT*.root");
    ch_ttbar_sl->Add(BabyDir+"baby_TTJets_SemiLeptMGDecays_8TeV_f*.root");
    ch_ttbar_ll->Add(BabyDir+"baby_TTJets_FullLeptMGDecays_8TeV_f*.root");
    // WJets 
    ch_wjets->Add(BabyDir+"baby_WJetsToLNu_HT-400ToInf_8TeV_f*.root");
    //ch_wjets->Add(BabyDir+"baby_W2*.root");
    //ch_wjets->Add(BabyDir+"baby_W3*.root");
    //ch_wjets->Add(BabyDir+"baby_W4*.root");
    // DY 
    ch_dy->Add(BabyDir+"baby_DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV_f*.root");
    // Singla top 
    ch_t->Add(BabyDir+"baby_*channel*_f*.root");

    // Signal
    ch_f1150_25->Add(BabyDir+"baby_*f1150_25.root");
    ch_f1200_25->Add(BabyDir+"baby_*f1200_25.root");
    ch_f1250_25->Add(BabyDir+"baby_*f1250_25.root");
    ch_f1300_25->Add(BabyDir+"baby_*f1300_25.root");
    ch_f1200_1000->Add(BabyDir+"baby_*f1200_1000.root");
    ch_f1100_125->Add(BabyDir+"baby_*f1100_125.root");
    ch_f1100_225->Add(BabyDir+"baby_*f1100_225.root");
    ch_f1100_325->Add(BabyDir+"baby_*f1100_325.root");
    ch_f1100_425->Add(BabyDir+"baby_*f1100_425.root");
    ch_f1100_525->Add(BabyDir+"baby_*f1100_525.root");
    ch_f1100_625->Add(BabyDir+"baby_*f1100_625.root");
    
    
    // ----------------------------------------
    //  Get number of entries 
    // ----------------------------------------
    cout << "data MuHad         : " << ch_data_muhad->GetEntries()        << endl;
    cout << "data ElectronHad   : " << ch_data_ehad->GetEntries()        << endl;
    cout << "ttbarl             : " << ch_ttbar_sl->GetEntries()    << endl;
    cout << "ttbarll            : " << ch_ttbar_ll->GetEntries()    << endl;
    cout << "wjets              : " << ch_wjets->GetEntries()       << endl;
    cout << "dy                 : " << ch_dy->GetEntries()          << endl;
    cout << "Single top         : " << ch_t->GetEntries()           << endl;
    cout << "T1tttt(1150,25)    : " << ch_f1150_25->GetEntries()    << endl;
    cout << "T1tttt(1200,25)    : " << ch_f1200_25->GetEntries()    << endl;
    cout << "T1tttt(1250,25)    : " << ch_f1250_25->GetEntries()    << endl;
    cout << "T1tttt(1300,25)    : " << ch_f1300_25->GetEntries()    << endl;
    cout << "T1tttt(1200,1000)  : " << ch_f1200_1000->GetEntries()  << endl;
    cout << "T1tttt(1100,125)   : " << ch_f1100_125->GetEntries()    << endl;
    cout << "T1tttt(1100,225)   : " << ch_f1100_225->GetEntries()    << endl;
    cout << "T1tttt(1100,325)   : " << ch_f1100_325->GetEntries()    << endl;
    cout << "T1tttt(1100,425)   : " << ch_f1100_425->GetEntries()    << endl;
    cout << "T1tttt(1100,525)   : " << ch_f1100_525->GetEntries()    << endl;
    cout << "T1tttt(1100,625)   : " << ch_f1100_625->GetEntries()    << endl;
    
    if(!OnlyDraw) 
    {
        // ----------------------------------------
        //  Fill histrograms 
        // ----------------------------------------
        DoOneProcess(ch_data_muhad,	    pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_data_ehad,	    pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_ttbar_sl,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_ttbar_ll,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_wjets,	    pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_dy,	        pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_t,	        pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_f1150_25,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1200_25,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1250_25,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1300_25,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1200_1000,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut); 
        DoOneProcess(ch_f1100_125,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1100_225,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1100_325,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1100_425,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1100_525,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  
        DoOneProcess(ch_f1100_625,	pTR0p5thres/*, Nbcut*/, METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);  

        // ----------------------------------------
        //  Make the final histogram file
        // ----------------------------------------
        cout << "[MJ Analysis] Merging result files" << endl;
        gSystem->Exec(Form("rm HistFiles/Hist_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root",(int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres));
        gSystem->Exec(Form("hadd -f HistFiles/Hist_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root HistFiles/*_MET%iTo%i_nj%iTo%i_nb%iTo%i_MJ%i_pT%i.root",
                    (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres,
                    (int)METlow,(int)METhigh,njlow,njhigh,nblow,nbhigh,(int)MJcut,pTR0p5thres));
    }

    // ----------------------------------------
    //  Draw histograms 
    // ---------------------------------------- 
///*
    DrawStack("mT",             30,	1/*,	Nbcut*/,METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);
    DrawStack("mj",             30,	1/*,	Nbcut*/,METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);
    DrawStack("MJ",       	    30,	1/*,	Nbcut*/,METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);
    DrawStack("HT",        	    30,	1/*,	Nbcut*/,METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);
    DrawStack("MET",       	    30,	1/*,	Nbcut*/,METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);
    DrawStack("Nfatjet",        30,	1/*,	Nbcut*/,METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);
    DrawStack("Nskinnyjet",     30,	1/*,	Nbcut*/,METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);
    //Draw2D("MJmT",              30/*,1*/, false/*,	Nbcut*/);
//*/
/*   
    DrawStack("muspT"           	);
    DrawStack("muspTminusMET"  ,30,1);
    DrawStack("musPhi"     	   ,30,1);
    DrawStack("musEta"          	);
    DrawStack("mT"              	);
    DrawStack("mj"             ,30,2);
    DrawStack("mj_mumatch"      	);
    DrawStack("mj_notmumatch"   	);
    DrawStack("mj_bmatch"       	);
    DrawStack("mj_notbmatch"    	);
    DrawStack("mj_metmatch"       	);
    DrawStack("mj_notmetmatch"    	);
    DrawStack("MJ"       	   ,30,2);
    DrawStack("HT"         			);
    DrawStack("Nfatjet"    			);
    DrawStack("Nskinnyjet" 			);
    DrawStack("MET"        			);
    DrawStack("METPhi"     	   ,30,1);
    DrawStack("METx"     	        );
    DrawStack("METy"     	        );
    DrawStack("DPhi"      	   ,30,2);	
    DrawStack("FatjetPt"   			);
    DrawStack("FatjetEta" 			);
    DrawStack("FatjetPt_mumatch"    );
    DrawStack("FatjetPt_notmumatch" );
    DrawStack("FatjetPt_bmatch"     );
    DrawStack("FatjetPt_notbmatch"  );
    DrawStack("FatjetPt_metmatch"   );
    DrawStack("FatjetPt_notmetmatch");
    DrawStack("WpT"                 );
    DrawStack("BjetPt_bmatchmu"   );
    DrawStack("BjetPt_notbmatchmu");
    DrawStack("dRbmu_bmatchmu");
    DrawStack("dRbmu_notbmatchmu");
    DrawStack("dRFJ");
    DrawStack("dPhiFJ");
    DrawStack("dEtaFJ");
    Draw2D("mjvsFatjetPt"           );
   
    
    DrawStack("FatjetPt1"   			);
    DrawStack("FatjetPt2"   			);
    DrawStack("FatjetPt3"   			);
    DrawStack("FatjetPt4"   			);
    DrawStack("FatjetPhi1"   			);
    DrawStack("FatjetPhi2"   			);
    DrawStack("FatjetPhi3"   			);
    DrawStack("FatjetPhi4"   			);
    DrawStack("FatjetEta1"   			);
    DrawStack("FatjetEta2"   			);
    DrawStack("FatjetEta3"   			);
    DrawStack("FatjetEta4"   			);
  
    DrawStack("dRFJ");
    DrawStack("dPhiFJ");
    DrawStack("dEtaFJ");
    DrawStack("mjFJFJ",30,2);
    DrawStack("mjFJFJnotclose",30,1);
    
    DrawStack("FatjetPtFJFJnotclose",30,1);
    DrawStack("FatjetEtaFJFJnotclose",30,1);
    DrawStack("FatjetPhiFJFJnotclose",30,1);

    DrawStack("FatjetPtFJFJcloseHigherPt",30,1);
    DrawStack("FatjetEtaFJFJcloseHigherPt",30,1);
    DrawStack("FatjetPhiFJFJcloseHigherPt",30,1);
    DrawStack("FatjetPtFJFJcloseLowerPt",30,1);
    DrawStack("FatjetEtaFJFJcloseLowerPt",30,1);
    DrawStack("FatjetPhiFJFJcloseLowerPt",30,1);
    
    Draw2D("mindRFJmjvsmj", 30, 1, true );
*/
    // Yield table
    PrintTable(METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut, 0, true);
//    PrintTable(METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut, 11, false);
//    PrintTable(METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut, 13, false);
    
    // Get Prediction 
    GetPrediction(METlow, METhigh, njlow, njhigh, nblow, nbhigh, MJcut);

}

