#include <iostream>
#include <vector>
#include <string>

#include "TROOT.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TLegend.h"
#include "TString.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TLorentzVector.h"

#include "Branch_v71.h"
#include "ObjectSelector_Sync.h"
#include "EventSelector.h"
#include "Utilities.h"
#include "inJSON2012.h"
#include "TOBTECFilter.h"
#include "filters.h"

// include necessary fastjet files
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

//
using namespace std;

//
typedef std::pair<fastjet::PseudoJet,int > FatJetPair;

// Btagging prob map function
TFile *f_tageff_;

//
// Make fat jets 
//
std::vector<FatJetPair> makeFatJet( vector<TLorentzVector> FatJetConstituent, 
                                   double Rparam=1.2, 
                                   int ConstituentpTcut=30, 
                                   float ConstituentEtacut=100) 
{
    std::vector<FatJetPair> FatJets_Pair; 

    // Loop over R=0.5 jets, form into PseudoJets vector
    vector<fastjet::PseudoJet> input_particles;
    double FatjetConstituent_px_tmp, FatjetConstituent_py_tmp, FatjetConstituent_pz_tmp, FatjetConstituent_energy_tmp;

    for(int ijet = 0; ijet<(int)FatJetConstituent.size(); ijet++) 
    { 

        FatjetConstituent_px_tmp        = FatJetConstituent.at(ijet).Px();
        FatjetConstituent_py_tmp        = FatJetConstituent.at(ijet).Py();
        FatjetConstituent_pz_tmp        = FatJetConstituent.at(ijet).Pz();
        FatjetConstituent_energy_tmp    = FatJetConstituent.at(ijet).E();	  
    
//        cout << ijet << " :: " 
//             << FatjetConstituent_px_tmp << " " 
//             << FatjetConstituent_py_tmp << " " 
//             << FatjetConstituent_pz_tmp << " " 
//             << FatjetConstituent_energy_tmp << " " 
//             << endl;

        if(TMath::Sqrt( FatjetConstituent_px_tmp*FatjetConstituent_px_tmp
                       +FatjetConstituent_py_tmp*FatjetConstituent_py_tmp)<ConstituentpTcut) continue;

        if(TMath::Abs(FatJetConstituent.at(ijet).Eta())>ConstituentEtacut) continue;

        input_particles.push_back(fastjet::PseudoJet( FatjetConstituent_px_tmp, FatjetConstituent_py_tmp,
                                                      FatjetConstituent_pz_tmp, FatjetConstituent_energy_tmp));
    }
    
    //
    // Run Fastjet to reconstuct jets 
    //

    // Create an object that represents your choice of jet algorithm and the associated parameters
    fastjet::Strategy strategy = fastjet::Best;
    fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, recomb_scheme, strategy);

    // run the jet clustering with the above jet definition
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);

    // 
    // Get p4 of the reconstructed jets  
    //
    double ptmin = 0.0; // could use 3.0 here, instead of applying later
    vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
    //Sort by pt
    vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(inclusive_jets);
    //fill fastjet output into vectors, continue as original code
    for(int isortjets = 0; isortjets< (int)sorted_jets.size(); isortjets++)
    {
        //store only if pt >3 GeV to match CMS jets
        if(TMath::Sqrt( sorted_jets[isortjets].px()*sorted_jets[isortjets].px()
                       +sorted_jets[isortjets].py()*sorted_jets[isortjets].py())>50) 
        {
            TLorentzVector FatJet_tmp( sorted_jets[isortjets].px(), sorted_jets[isortjets].py(), 
                                       sorted_jets[isortjets].pz(), sorted_jets[isortjets].E());
            FatJets_Pair.push_back(std::make_pair(sorted_jets.at(isortjets),sorted_jets.at(isortjets).constituents().size()));

//            cout << isortjets << " :: "  
//                 << sorted_jets[isortjets].px() << " " 
//                 << sorted_jets[isortjets].py() << " " 
//                 << sorted_jets[isortjets].pz() << " " 
//                 << sorted_jets[isortjets].E() << " " 
//                 << endl;
        }
    }

    return FatJets_Pair;
}


//
// Get mj 
//
float Getmj(double px, double py, double pz, double E)
{

    float mj = TMath::Sqrt(E*E - px*px - py*py - pz*pz); 
    return mj;
}

//
// Get MJ 
//
float GetMJ(vector<float> Vectormj)
{
    
    if(Vectormj.size()>0) 
    {
        float MJ = 0.;
        for(int imj=0; imj<(int)Vectormj.size(); imj++)
        {    
            MJ = MJ + Vectormj.at(imj);
        }
        return MJ;
    } else 
    {
        return -999.;
    }
}

bool hasGoodVertex(int & firstGoodVertex)
{
    bool GoodEvent=false;

    for(unsigned i=0; i<pv_x->size(); i++) {
        if( fabs(pv_isFake->at(i)) < 0.001
                && pv_ndof->at(i) > 4
                // commented out to compare with Robert
                //      && pv_ndof->at(i) >= 4
                && sqrt(pv_x->at(i)*pv_x->at(i)+pv_y->at(i)*pv_y->at(i)) < 2
                && fabs(pv_z->at(i)) < 24)
        {
            if(GoodEvent==false) firstGoodVertex=i;
            GoodEvent=true;
        }
    }
    return GoodEvent;
}

//
// Btag probability from Jack : https://github.com/jbradmil/csa14/blob/master/src/event_handler.cpp
//
string AssembleBTagEffFilename(string sampleName) 
{
    cout << sampleName << endl;
    unsigned found = sampleName.find_last_of("/");
    cout << found << endl;
    std::string truncated_name = sampleName.substr(0,found);
    cout << truncated_name<< endl;
    found = truncated_name.find_last_of("/");
    cout << found << endl;
    truncated_name = truncated_name.substr(found+1);
    cout << truncated_name<< endl;
    cout << "histos_btageff_csvm_"+truncated_name+".root" << endl;
    return "histos_btageff_csvm_"+truncated_name+".root";
}

void LoadJetTagEffMaps(string sampleName) 
{

    int cfAVersion = 71;

    if (cfAVersion<=71) 
    {
        assert(f_tageff_ ==0);
        TString filename = AssembleBTagEffFilename(sampleName);
        cout << filename << endl;
        filename.Prepend("btagEffMaps/");
        f_tageff_ = new TFile(filename,"READ");
        if (f_tageff_->IsZombie()) 
        {
            cout<<"Failed to load the b-tag eff map for sample "<<sampleName<<endl;
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
double GetJetTagEff(unsigned int ijet, TH1D* h_btageff, TH1D* h_ctageff, TH1D* h_ltageff) 
{

    double tageff=0;
    const float pt = jets_AK5PF_pt->at(ijet);
    const float eta = fabs(jets_AK5PF_eta->at(ijet));
    int flavor = static_cast<int>(jets_AK5PF_partonFlavour->at(ijet));

    //x is the pt value that will be used to evaluate the SF
    //the max pt depends on flavor, and, for LF jets, eta
    const double cutoff1=800.; const double cutoff2=700.;
    double x;
    //HF or central LF, max is 800
    if ( abs(flavor)==5 || abs(flavor)==4 || eta<1.5) x = pt > cutoff1 ? cutoff1 : pt;
    //high eta LF, max is 700
    else x = pt > cutoff2 ? cutoff2 : pt;

    //  if (theBTagEffType_ == kBTagEff05 || theBTagEffType_ == kBTagEffup5 || theBTagEffType_ == kBTagEffdown5{ //new  BTV-11-004 prescription 

    if (abs(flavor) ==4 || abs(flavor)==5) 
    { //heavy flavor
        double errFactor = 1;

        if (pt >cutoff1) 
        { //use twice the errors
            errFactor=2;
        }
        if (abs(flavor) == 4)  errFactor=2; //charm has double the errors   "SFc = SFb with twice the quoted uncertainty"
        //not clear to me what to do for charm with pT> cutoff. errFactor of 2 or 4? must be small effect though

        // Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
        double  SFb = 0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)));

        //apply FASTSTIM correction where needed
        //SFb *= GetbJetFastsimSF("value",flavor,pt);

        // skip syst variations, for now

        // cout<<"jet flavor, pt, SF = "<<abs(flavor)<<" "<<pt<<" "<<SFb<<endl;
        if      (abs(flavor) == 5) tageff = SFb * h_btageff->GetBinContent( h_btageff->FindBin( pt ) );
        else if (abs(flavor) == 4) tageff = SFb * h_ctageff->GetBinContent( h_ctageff->FindBin( pt ) );
        else assert(0);
    } // if heavy flavor
    else 
    { //light flavor [ see https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_Moriond2013.C ]
        //note that these are valid only to a cutoff value, so use 'x' and not 'pt'
        double SF=0;
        double nominalSF=0;
        if ( eta < 0.8 ) 
        {
            nominalSF =  ((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x))); // SF without + or - variation in mistag rate
            SF = nominalSF;
        }
        else if (eta>=0.8 && eta<1.6) 
        {
            nominalSF = ((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)));
            SF = nominalSF;
        }
        else if (eta>=1.6 && eta<=2.4) 
        {
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
                       vector<int> GoodJets,  string sampleName) 
{

    //must initialize correctly
    Prob2 = 0;
    Prob1 = 0; ProbGEQ1 = 1; Prob0 = 1; ProbGEQ2 = 0;
    Prob3 = 0; ProbGEQ4 = 0;

    if(sampleName.find("Run2012")!=std::string::npos) 
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

    for(int igoodjet=0; igoodjet<(int)GoodJets.size(); igoodjet++)
    {
        int ijet = GoodJets.at(igoodjet);

        double subprob1=0;
        double subprob2=0;

        double effi = GetJetTagEff(ijet, h_btageff, h_ctageff, h_ltageff);
        //      cout << "jet: " << ijet << ", effi: " << effi << endl;
        Prob0 = Prob0* ( 1 - effi);

        double product = 1;
        for(int kgoodjet=0; kgoodjet<(int)GoodJets.size(); kgoodjet++)
        {
            int kjet = GoodJets.at(kgoodjet);

            double effk = GetJetTagEff(kjet, h_btageff, h_ctageff, h_ltageff);
            if(kjet != ijet) product = product*(1-effk);
            if(kjet > ijet)
            {
                double subproduct = 1;
                for(int jgoodjet=0; jgoodjet<(int)GoodJets.size(); jgoodjet++)
                {
                    int jjet = GoodJets.at(jgoodjet);

                    if(jjet != kjet && jjet != ijet)
                    {
                        double effj = GetJetTagEff(jjet, h_btageff, h_ctageff, h_ltageff);
                        subproduct = subproduct*(1-effj);
                        if ( jjet > kjet) 
                        {
                            double subproduct2 = 1;
                            for(int lgoodjet=0; lgoodjet<(int)GoodJets.size(); lgoodjet++)
                            {
                                int ljet = GoodJets.at(lgoodjet);

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

//
// main function
//
void DoOneProcess(TString InputName, TString ProcessName, int ibegin, int iend, bool isData, float Lumi) 
{
    //
    cout << "[MJ Analysis] ------------------------------------------------------------------------------------"<<endl; 
    cout << "[MJ Analysis] Processing : " << ProcessName << endl;
    cout << "[MJ Analysis] Input dir  : " << InputName << endl;

    TFile *babyFile_ = new TFile(Form("baby_%s_f%iTo%i.root", ProcessName.Data(), ibegin, iend), "RECREATE");
    babyFile_->cd();
    TTree *babyTree_ = new TTree("tree", "A Baby Ntuple");
    
    // 
    // Get tree 
    // 
    TChain * chainA = new TChain("/configurableAnalysis/eventA");   
    TChain * chainB = new TChain("/configurableAnalysis/eventB");  
    
    for(int i=ibegin; i<=iend; i++)  
    {
        gSystem->Exec(Form("ls %s/*f%i_*.root", InputName.Data(), i));
        chainA->Add(Form("%s/*f%i_*.root", InputName.Data(), i));
        chainB->Add(Form("%s/*f%i_*.root", InputName.Data(), i));
        //gSystem->Exec(Form("ls %s/configurableAnalysis_%i_*.root", InputName.Data(), i)); // for non-published cfA samples 
        //chainA->Add(Form("%s/configurableAnalysis_%i_*.root", InputName.Data(), i));
        //chainB->Add(Form("%s/configurableAnalysis_%i_*.root", InputName.Data(), i));
    } 
    TList *l = (TList*)chainA->GetListOfFiles();
    l->Print();

    InitializeA(chainA);
    InitializeB(chainB);
    
    // Get total number of events of a sample 
    int TotalNEntries=1;
    
    //if(!isData && !ProcessName.Contains("Test"))
    if(!isData)
    { 
        TChain * chainATotal = new TChain("/configurableAnalysis/eventA");  
        chainATotal->Add(Form("%s/*.root", InputName.Data()));
        TotalNEntries = (int)chainATotal->GetEntries();
    }

    //
    // Baby variables 
    //
    int run_; 
    int lumiblock_; 
    int event_; 
    vector<bool> filter_;  
    bool TrigMuon_;  
    bool TrigSingleMuon_;  
    bool TrigDiMuon_;  
    bool TrigHTMuon_;  
    bool TrigElectron_;  
    bool TrigSingleElectron_;  
    bool TrigDiElectron_;  
    bool TrigHTElectron_;  
    bool TrigMuEG_;  
    int Nfatjet_pT20_; 
    int Nfatjet_pT30_; 
    int NfatjetCHS_pT20_; 
    int NfatjetCHS_pT30_; 
    int Nskinnyjet_; 
    int NBtagCSVM_; 
    int NskinnyjetCHS_; 
    int NBtagCHSCSVM_; 
    int Npv_; 
    float Npu_; 
    float EventWeight_; 
    float MJ_pT20_; 
    float MJ_pT30_; 
    float MJCHS_pT20_; 
    float MJCHS_pT30_; 
    float MET_;  
    float METPhi_;  
    float HT_;  
    float HTCHS_;  
    float top1pT_;  
    float top1Eta_;  
    float top1Phi_;  
    float top2pT_;  
    float top2Eta_;  
    float top2Phi_;  
    float Prob0_;  
    float Prob1_;  
    float Prob2_;  
    float Prob3_;  
    vector<float> mj_pT20_;
    vector<float> mj_pT30_;
    vector<float> mjCHS_pT20_;
    vector<float> mjCHS_pT30_;
    vector<float> mjLep_pT30_;
    vector<float> mjCHSLep_pT30_;
    vector<float> FatjetPt_pT20_;
    vector<float> FatjetEta_pT20_;
    vector<float> FatjetPhi_pT20_;
    vector<float> FatjetN_pT20_;
    vector<float> FatjetCHSPt_pT20_;
    vector<float> FatjetCHSEta_pT20_;
    vector<float> FatjetCHSPhi_pT20_;
    vector<float> FatjetCHSN_pT20_;
    vector<float> FatjetPt_pT30_;
    vector<float> FatjetEta_pT30_;
    vector<float> FatjetPhi_pT30_;
    vector<float> FatjetN_pT30_;
    vector<float> FatjetCHSPt_pT30_;
    vector<float> FatjetCHSEta_pT30_;
    vector<float> FatjetCHSPhi_pT30_;
    vector<float> FatjetCHSN_pT30_;
    vector<float> FatjetLepPt_pT30_;
    vector<float> FatjetLepEta_pT30_;
    vector<float> FatjetLepPhi_pT30_;
    vector<float> FatjetLepN_pT30_;
    vector<float> FatjetCHSLepPt_pT30_;
    vector<float> FatjetCHSLepEta_pT30_;
    vector<float> FatjetCHSLepPhi_pT30_;
    vector<float> FatjetCHSLepN_pT30_;
    vector<float> RA4ElsPt_;
    vector<float> RA4ElsEta_;
    vector<float> RA4ElsPhi_;
    vector<float> RA4ElsE_;
    vector<float> RA4ElsQ_;
    vector<float> RA4MusPt_;
    vector<float> RA4MusEta_;
    vector<float> RA4MusPhi_;
    vector<float> RA4MusE_;
    vector<float> RA4MusQ_;
    vector<float> RA4ElsVetoPt_;
    vector<float> RA4ElsVetoEta_;
    vector<float> RA4ElsVetoPhi_;
    vector<float> RA4MusVetoPt_;
    vector<float> RA4MusVetoEta_;
    vector<float> RA4MusVetoPhi_;
    vector<float> JetPt_;
    vector<float> JetEta_;
    vector<float> JetPhi_;
    vector<float> JetE_;
    vector<float> JetCSV_;
    vector<float> JetMCId_;
    vector<float> JetCHSPt_;
    vector<float> JetCHSEta_;
    vector<float> JetCHSPhi_;
    vector<float> JetCHSE_;
    vector<float> JetCHSCSV_;
    vector<float> JetCHSMCId_;
    vector<float> JetCHSNCh_;
    vector<float> JetCHSNNeu_;
    vector<float> GenPt_;
    vector<float> GenEta_;
    vector<float> GenPhi_;
    vector<float> GenE_;
    vector<float> GenId_;
    vector<float> GenMId_;
    vector<float> GenGMId_;
    
    babyTree_->Branch("run",            	&run_);    
    babyTree_->Branch("lumiblock",      	&lumiblock_); 
    babyTree_->Branch("event",          	&event_);     
    babyTree_->Branch("filter",          	&filter_);     
    babyTree_->Branch("TrigMuon",          	&TrigMuon_);     
    babyTree_->Branch("TrigSingleMuon",    	&TrigSingleMuon_);     
    babyTree_->Branch("TrigDiMuon",    	    &TrigDiMuon_);     
    babyTree_->Branch("TrigHTMuon",       	&TrigHTMuon_);     
    babyTree_->Branch("TrigElectron",      	&TrigElectron_);     
    babyTree_->Branch("TrigSingleElectron",	&TrigSingleElectron_);     
    babyTree_->Branch("TrigDiElectron", 	&TrigDiElectron_);     
    babyTree_->Branch("TrigHTElectron",    	&TrigHTElectron_);     
    babyTree_->Branch("TrigMuEG",    	    &TrigMuEG_);     
    babyTree_->Branch("Nfatjet_pT20",   	&Nfatjet_pT20_);   
    babyTree_->Branch("Nfatjet_pT30",   	&Nfatjet_pT30_);   
    babyTree_->Branch("NfatjetCHS_pT20", 	&NfatjetCHS_pT20_);   
    babyTree_->Branch("NfatjetCHS_pT30", 	&NfatjetCHS_pT30_);   
    babyTree_->Branch("Nskinnyjet",     	&Nskinnyjet_);
    babyTree_->Branch("NBtagCSVM",     	    &NBtagCSVM_);
    babyTree_->Branch("NskinnyjetCHS",  	&NskinnyjetCHS_);
    babyTree_->Branch("NBtagCHSCSVM",     	&NBtagCHSCSVM_);
    babyTree_->Branch("Npv",            	&Npv_);       
    babyTree_->Branch("Npu",            	&Npu_);       
    babyTree_->Branch("EventWeight",    	&EventWeight_);
    babyTree_->Branch("MJ_pT20",        	&MJ_pT20_);        
    babyTree_->Branch("MJ_pT30",        	&MJ_pT30_);        
    babyTree_->Branch("MJCHS_pT20",        	&MJCHS_pT20_);        
    babyTree_->Branch("MJCHS_pT30",        	&MJCHS_pT30_);        
    babyTree_->Branch("MET",            	&MET_);        
    babyTree_->Branch("METPhi",            	&METPhi_);        
    babyTree_->Branch("HT",             	&HT_);        
    babyTree_->Branch("HTCHS",             	&HTCHS_);        
    babyTree_->Branch("top1pT",          	&top1pT_);        
    babyTree_->Branch("top1Eta",          	&top1Eta_);        
    babyTree_->Branch("top1Phi",          	&top1Phi_);        
    babyTree_->Branch("top2pT",          	&top2pT_);        
    babyTree_->Branch("top2Eta",          	&top2Eta_);        
    babyTree_->Branch("top2Phi",          	&top2Phi_);        
    babyTree_->Branch("Prob0",          	&Prob0_);        
    babyTree_->Branch("Prob1",          	&Prob1_);        
    babyTree_->Branch("Prob2",          	&Prob2_);        
    babyTree_->Branch("Prob3",          	&Prob3_);        
    babyTree_->Branch("mj_pT20",        	&mj_pT20_);     
    babyTree_->Branch("mj_pT30",        	&mj_pT30_);     
    babyTree_->Branch("mjCHS_pT20",        	&mjCHS_pT20_);     
    babyTree_->Branch("mjCHS_pT30",        	&mjCHS_pT30_);     
    babyTree_->Branch("mjLep_pT30",        	&mjLep_pT30_);     
    babyTree_->Branch("mjCHSLep_pT30",        	&mjCHSLep_pT30_);     
    babyTree_->Branch("FatjetPt_pT20",  	&FatjetPt_pT20_); 
    babyTree_->Branch("FatjetEta_pT20", 	&FatjetEta_pT20_);
    babyTree_->Branch("FatjetPhi_pT20",     &FatjetPhi_pT20_);
    babyTree_->Branch("FatjetN_pT20",       &FatjetN_pT20_);
    babyTree_->Branch("FatjetCHSPt_pT20",   &FatjetCHSPt_pT20_); 
    babyTree_->Branch("FatjetCHSEta_pT20",  &FatjetCHSEta_pT20_);
    babyTree_->Branch("FatjetCHSPhi_pT20",  &FatjetCHSPhi_pT20_);
    babyTree_->Branch("FatjetCHSN_pT20",    &FatjetCHSN_pT20_);
    babyTree_->Branch("FatjetPt_pT30",  	&FatjetPt_pT30_); 
    babyTree_->Branch("FatjetEta_pT30", 	&FatjetEta_pT30_);
    babyTree_->Branch("FatjetPhi_pT30",     &FatjetPhi_pT30_);
    babyTree_->Branch("FatjetN_pT30",       &FatjetN_pT30_);
    babyTree_->Branch("FatjetCHSPt_pT30",   &FatjetCHSPt_pT30_); 
    babyTree_->Branch("FatjetCHSEta_pT30",  &FatjetCHSEta_pT30_);
    babyTree_->Branch("FatjetCHSPhi_pT30",  &FatjetCHSPhi_pT30_);
    babyTree_->Branch("FatjetCHSN_pT30",    &FatjetCHSN_pT30_);
    babyTree_->Branch("FatjetLepPt_pT30",   &FatjetLepPt_pT30_); 
    babyTree_->Branch("FatjetLepEta_pT30",  &FatjetLepEta_pT30_);
    babyTree_->Branch("FatjetLepPhi_pT30",  &FatjetLepPhi_pT30_);
    babyTree_->Branch("FatjetLepN_pT30",    &FatjetLepN_pT30_);
    babyTree_->Branch("FatjetCHSLepPt_pT30",   &FatjetCHSLepPt_pT30_); 
    babyTree_->Branch("FatjetCHSLepEta_pT30",  &FatjetCHSLepEta_pT30_);
    babyTree_->Branch("FatjetCHSLepPhi_pT30",  &FatjetCHSLepPhi_pT30_);
    babyTree_->Branch("FatjetCHSLepN_pT30",    &FatjetCHSLepN_pT30_);
    babyTree_->Branch("RA4ElsPt",           &RA4ElsPt_);
    babyTree_->Branch("RA4ElsEta",          &RA4ElsEta_);
    babyTree_->Branch("RA4ElsPhi",          &RA4ElsPhi_);
    babyTree_->Branch("RA4ElsE",            &RA4ElsE_);
    babyTree_->Branch("RA4ElsQ",            &RA4ElsQ_);
    babyTree_->Branch("RA4MusPt",           &RA4MusPt_);
    babyTree_->Branch("RA4MusEta",          &RA4MusEta_);
    babyTree_->Branch("RA4MusPhi",          &RA4MusPhi_);
    babyTree_->Branch("RA4MusE",            &RA4MusE_);
    babyTree_->Branch("RA4MusQ",            &RA4MusQ_);
    babyTree_->Branch("RA4ElsVetoPt",       &RA4ElsVetoPt_);
    babyTree_->Branch("RA4ElsVetoEta",      &RA4ElsVetoEta_);
    babyTree_->Branch("RA4ElsVetoPhi",      &RA4ElsVetoPhi_);
    babyTree_->Branch("RA4MusVetoPt",       &RA4MusVetoPt_);
    babyTree_->Branch("RA4MusVetoEta",      &RA4MusVetoEta_);
    babyTree_->Branch("RA4MusVetoPhi",      &RA4MusVetoPhi_);
    babyTree_->Branch("JetPt",              &JetPt_);
    babyTree_->Branch("JetEta",             &JetEta_);
    babyTree_->Branch("JetPhi",             &JetPhi_);
    babyTree_->Branch("JetE",               &JetE_);
    babyTree_->Branch("JetCSV",             &JetCSV_);
    babyTree_->Branch("JetMCId",            &JetMCId_);
    babyTree_->Branch("JetCHSPt",           &JetCHSPt_);
    babyTree_->Branch("JetCHSEta",          &JetCHSEta_);
    babyTree_->Branch("JetCHSPhi",          &JetCHSPhi_);
    babyTree_->Branch("JetCHSE",            &JetCHSE_);
    babyTree_->Branch("JetCHSCSV",          &JetCHSCSV_);
    babyTree_->Branch("JetCHSMCId",         &JetCHSMCId_);
    babyTree_->Branch("JetCHSNCh",          &JetCHSNCh_);
    babyTree_->Branch("JetCHSNNeu",         &JetCHSNNeu_);
    babyTree_->Branch("GenPt",              &GenPt_);
    babyTree_->Branch("GenEta",             &GenEta_);
    babyTree_->Branch("GenPhi",             &GenPhi_);
    babyTree_->Branch("GenE",               &GenE_);
    babyTree_->Branch("GenId",              &GenId_);
    babyTree_->Branch("GenMId",             &GenMId_);
    babyTree_->Branch("GenGMId",            &GenGMId_);
   
    // 
    // Event weights
    //
    // PU for MC
    TFile *fPUFile = TFile::Open("aux/puWeights_Summer12_53x_True_19p5ifb.root"); 
    TH1F *h1PU = (TH1F*)(fPUFile->Get("puWeights"));
    
    // Get cross section for MC 
    float Xsec=1; 
    if(!isData) Xsec=GetXsec(ProcessName);

    // json for DATA 
    std::vector< std::vector<int> > VRunLumi; VRunLumi.clear();
    if(isData) 
    {
        if(ProcessName.Contains("PromptReco")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("13Jul2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt"<< endl;       
        } else if(ProcessName.Contains("06Aug2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("24Aug2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("11Dec2012")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt" << endl;       
        } else if(ProcessName.Contains("22Jan2013")) 
        { 
            VRunLumi = MakeVRunLumi("json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"); 
            cout << "[MJ Analysis] Running with JSON : " << "Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" << endl;       
        } else 
        {
            cout << "[Error] No proper choice of JSON files!!" << endl;
            return ;
        }
    } else 
    {
        cout << "[MJ Analysis] No JSON files applied because it is MC" << endl;
    }

    //LoadJetTagEffMaps((InputName+"/").Data());
    LoadJetTagEffMaps((InputName).Data());

    //
    // main event loop
    //
    Int_t nentries = (Int_t)chainA->GetEntries();
    cout<<"[MJ Analysis] Number of entries is: "<<nentries<<endl;
    cout<<"[MJ Analysis] Number of entries of total sample is: "<<TotalNEntries<<endl;
    // Progress tracking 
    int i_permille_old = 0; 
    TDatime DTStart;
    int StartDate = DTStart.GetDate(); 
    int StartTime = DTStart.GetTime(); 
    cout << "[MJ Analysis] Start time : " << (StartTime/10000)%100 << ":"
         << (StartTime/100)%100 << ":" << StartTime%100
         << endl; 
    
    for(int i = 0; i<nentries; i++) 
    {
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
//                printf("\015\033[32m Processed :: \033[1m\033[31m%4.1f %%" 
//                       "\033[0m\033[32m   Expected processing time :: \033[1m\033[31m%i:%i:%i \033[0m\015",
//                        i_permille/10., (TimeToRun/10000)%100<60 ? (TimeToRun/10000)%100 : (TimeToRun/10000)%100-40, 
//                                        (TimeToRun/100)%100<60 ? (TimeToRun/100)%100 : (TimeToRun/100)%100-40, 
//                                        (TimeToRun%100)<60 ? (TimeToRun)%100 : (TimeToRun)%100-40 );
//                fflush(stdout);
            }
            i_permille_old = i_permille;
        }
        if(i%10000==0) cout << "Progress : " << i << "/" << nentries << endl;
        // Progress indicator end ----------------------------------
        
        // Access to the event  
        chainA->GetEntry(i);
        chainB->GetEntry(i);
      
        // initialize baby variables
        run_                =   -1;
        lumiblock_          =   -1;
        event_              =   -1;
        TrigMuon_           =   1;
        TrigSingleMuon_     =   1;
        TrigDiMuon_         =   1;
        TrigHTMuon_         =   1;
        TrigElectron_       =   1;
        TrigSingleElectron_ =   1;
        TrigDiElectron_     =   1;
        TrigHTElectron_     =   1;
        TrigMuEG_           =   1;
        Nfatjet_pT20_       =   -1;
        Nfatjet_pT30_       =   -1;
        NfatjetCHS_pT20_    =   -1;
        NfatjetCHS_pT30_    =   -1;
        Nskinnyjet_         =   -1;
        NBtagCSVM_          =   -1;
        NskinnyjetCHS_      =   -1;
        NBtagCHSCSVM_       =   -1;
        Npv_                =   -1;
        Npu_                =   -1;
        EventWeight_        =   1.;
        MJ_pT20_            =-999.;
        MJ_pT30_            =-999.;
        MJCHS_pT20_         =-999.;
        MJCHS_pT30_         =-999.;
        MET_                =-999.;
        METPhi_             =-999.;
        HT_                 =-999.;
        HTCHS_              =-999.;
        top1pT_             =-999.;
        top1Phi_            =-999.;
        top1Eta_            =-999.;
        top2pT_             =-999.;
        top2Eta_            =-999.;
        top2Phi_            =-999.;
        Prob0_              =-999.;
        Prob1_              =-999.;
        Prob2_              =-999.;
        Prob3_              =-999.;
        filter_.clear(); 
        mj_pT20_.clear();
        mj_pT30_.clear();
        mjCHS_pT20_.clear();
        mjCHS_pT30_.clear();
        mjLep_pT30_.clear();
        mjCHSLep_pT30_.clear();
        FatjetPt_pT20_.clear();
        FatjetEta_pT20_.clear();
        FatjetPhi_pT20_.clear();
        FatjetN_pT20_.clear();
        FatjetCHSPt_pT20_.clear();
        FatjetCHSEta_pT20_.clear();
        FatjetCHSPhi_pT20_.clear();
        FatjetCHSN_pT20_.clear();
        FatjetPt_pT30_.clear();
        FatjetEta_pT30_.clear();
        FatjetPhi_pT30_.clear();
        FatjetN_pT30_.clear();
        FatjetCHSPt_pT30_.clear();
        FatjetCHSEta_pT30_.clear();
        FatjetCHSPhi_pT30_.clear();
        FatjetCHSN_pT30_.clear();
        FatjetLepPt_pT30_.clear();
        FatjetLepEta_pT30_.clear();
        FatjetLepPhi_pT30_.clear();
        FatjetLepN_pT30_.clear();
        FatjetCHSLepPt_pT30_.clear();
        FatjetCHSLepEta_pT30_.clear();
        FatjetCHSLepPhi_pT30_.clear();
        FatjetCHSLepN_pT30_.clear();
        RA4ElsPt_.clear();
        RA4ElsEta_.clear();
        RA4ElsPhi_.clear();
        RA4ElsE_.clear();
        RA4ElsQ_.clear();
        RA4MusPt_.clear();
        RA4MusEta_.clear();
        RA4MusPhi_.clear();
        RA4MusE_.clear();
        RA4MusQ_.clear();
        RA4ElsVetoPt_.clear();
        RA4ElsVetoEta_.clear();
        RA4ElsVetoPhi_.clear();
        RA4MusVetoPt_.clear();
        RA4MusVetoEta_.clear();
        RA4MusVetoPhi_.clear();
        JetPt_.clear();
        JetEta_.clear();
        JetPhi_.clear();
        JetE_.clear();
        JetCSV_.clear();
        JetMCId_.clear();
        JetCHSPt_.clear();
        JetCHSEta_.clear();
        JetCHSPhi_.clear();
        JetCHSE_.clear();
        JetCHSCSV_.clear();
        JetCHSMCId_.clear();
        JetCHSNCh_.clear();
        JetCHSNNeu_.clear();
        GenPt_.clear();
        GenEta_.clear();
        GenPhi_.clear();
        GenE_.clear();
        GenId_.clear();
        GenMId_.clear();
        GenGMId_.clear();

        //
        // Core analysis 
        //
        // Get event weight 
        float EventWeight = 1;
        if(!isData) 
        {
            EventWeight = Xsec/TotalNEntries*Lumi; // scale for 1 pb * Lumi
            // need more weights if needed 
        } 
        else 
        {
            if(!inJSON(VRunLumi,run,lumiblock)) continue; // JSON
        }

        // Event quality 
        int firstGoodVertex=0;
        bool GoodVertex = hasGoodVertex(firstGoodVertex);
        if(isData && !GoodVertex) continue;
        
        // filters
        if(isData)
        {
            filter_.push_back(scrapingVeto_decision);
            filter_.push_back(hbhefilter_decision);
            filter_.push_back(trackingfailurefilter_decision);
            filter_.push_back(cschalofilter_decision);
            filter_.push_back(eebadscfilter_decision);
            filter_.push_back(ecalTPfilter_decision);
            filter_.push_back(ecallaserfilter_decision);
            filter_.push_back(trackertoomanyclustersfilter_decision);
            filter_.push_back(trackertoomanytripletsfilter_decision);
            filter_.push_back(trackertoomanyseedsfilter_decision);
            filter_.push_back(ecalBEfilter_decision);
            filter_.push_back(greedymuonfilter_decision);
            filter_.push_back(inconsistentPFmuonfilter_decision);
            filter_.push_back(hcallaserfilter_decision);
            filter_.push_back(eenoisefilter_decision);
            filter_.push_back(trackercoherentnoisefilter1_decision);
            filter_.push_back(trackercoherentnoisefilter2_decision);
            filter_.push_back(TOBTEC_ok());
            
            //filters filterList("Filter/AllBad_HCAL_ECAL_Laser.dat");
            //bool BadHCALECALfilter = filterList.inFilterList(run,lumiblock,event);
            //filter_.push_back(BadHCALECALfilter);
            //cout << TOBTEC_ok() << " " << BadHCALECALfilter << endl;
            //cout << event << " " << TOBTEC_ok() << endl;
        }

        // Get good RA4 muons
        vector<int> RA4MuonVeto; RA4MuonVeto.clear();
        vector<int> RA4Muon = GetRA4Muon(RA4MuonVeto,"", firstGoodVertex);
        // Get good RA4 electrons
        vector<int> RA4ElecVeto; RA4ElecVeto.clear();
        vector<int> RA4Elec = GetRA4Elec(RA4ElecVeto, "", firstGoodVertex, true);
        // Get good skinny jets, HT and B-tagged jets 
        double HT=-999.; 
        vector<int> LooseBJet; 
        vector<int> MediumBJet; 
        vector<int> GoodJets_AK5PFclean = GetJets(RA4Muon,RA4Elec,RA4MuonVeto,RA4ElecVeto,
                                                  HT,LooseBJet,MediumBJet, 
                                                  2.4, 40, 0.3); 
        double HTCHS=-999.; 
        vector<int> LooseBJetCHS; 
        vector<int> MediumBJetCHS; 
        vector<int> GoodJets_AK5PF = GetJetsCHS(RA4Muon,RA4Elec,RA4MuonVeto,RA4ElecVeto,
                                                HTCHS,LooseBJetCHS,MediumBJetCHS, 
                                                2.4, 40, 0.3); 
        //
        // Skim : HT>500 && MET>100 OR dilepton
        //
        // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  
        if( (RA4Muon.size()+RA4Elec.size())==0) continue;
        bool Skim = ((HT>500 || HTCHS>500) && pfTypeImets_et->at(0)>100) || (RA4Muon.size()+RA4Elec.size())>1;
        if(!Skim) continue;
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  

        double Prob0, ProbGEQ1, Prob1, ProbGEQ2,
               Prob2, ProbGEQ3, Prob3, ProbGEQ4;
        CalculateTagProbs(Prob0, ProbGEQ1, Prob1, ProbGEQ2,
                          Prob2, ProbGEQ3, Prob3, ProbGEQ4,
                          GoodJets_AK5PF, InputName.Data());
        //cout << Prob0<< " " << ProbGEQ1<< " " << Prob1<< " " << ProbGEQ2 << " " // DEBUG
        //     << Prob2<< " " << ProbGEQ3<< " " << Prob3<< " " << ProbGEQ4 << " "
        //     << endl; 
        Prob0_ = Prob0;
        Prob1_ = Prob1;
        Prob2_ = Prob2;
        Prob3_ = Prob3;
        
        /*
        // 
        // pT(R=0.5) > 20 GeV
        // 
        
        // Regular jets
        double MJ_pT20=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjet_pT20; 
        vector<int> Vector_GoodFatjet_pT20_Index; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PFclean_R1p2_R0p5pT20_px->size(); ifatjet++) 
        {
            float temp_pT_pT20 = TMath::Sqrt(fastjets_AK5PFclean_R1p2_R0p5pT20_px->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT20_px->at(ifatjet)
                                        +fastjets_AK5PFclean_R1p2_R0p5pT20_py->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT20_py->at(ifatjet));
            if(temp_pT_pT20<50) continue;
            TLorentzVector temp_GoodFatjet_pT20( fastjets_AK5PFclean_R1p2_R0p5pT20_px->at(ifatjet), 
                                            fastjets_AK5PFclean_R1p2_R0p5pT20_py->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT20_pz->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT20_energy->at(ifatjet));
            Vector_GoodFatjet_pT20.push_back(temp_GoodFatjet_pT20);
            Vector_GoodFatjet_pT20_Index.push_back(ifatjet);
        } 
       
        vector<float> Vector_mj_pT20;   // mj
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT20.size(); igoodfatjet++) 
        {
            float temp_mj_pT20 = Getmj(Vector_GoodFatjet_pT20.at(igoodfatjet).Px(), 
                                  Vector_GoodFatjet_pT20.at(igoodfatjet).Py(),
                                  Vector_GoodFatjet_pT20.at(igoodfatjet).Pz(),
                                  Vector_GoodFatjet_pT20.at(igoodfatjet).E());
            Vector_mj_pT20.push_back(temp_mj_pT20);
        }
        MJ_pT20 = GetMJ(Vector_mj_pT20);

        int Nfatjet_pT20 = Vector_GoodFatjet_pT20.size(); 
       
        // CHS jets
        double MJCHS_pT20=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjetCHS_pT20; 
        vector<int> Vector_GoodFatjetCHS_pT20_Index; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PF_R1p2_R0p5pT20_px->size(); ifatjet++) 
        {
            float temp_pT_pT20 = TMath::Sqrt(fastjets_AK5PF_R1p2_R0p5pT20_px->at(ifatjet)*fastjets_AK5PF_R1p2_R0p5pT20_px->at(ifatjet)
                                        +fastjets_AK5PF_R1p2_R0p5pT20_py->at(ifatjet)*fastjets_AK5PF_R1p2_R0p5pT20_py->at(ifatjet));
            if(temp_pT_pT20<50) continue;
            TLorentzVector temp_GoodFatjetCHS_pT20( fastjets_AK5PF_R1p2_R0p5pT20_px->at(ifatjet), 
                                            fastjets_AK5PF_R1p2_R0p5pT20_py->at(ifatjet),
                                            fastjets_AK5PF_R1p2_R0p5pT20_pz->at(ifatjet),
                                            fastjets_AK5PF_R1p2_R0p5pT20_energy->at(ifatjet));
            Vector_GoodFatjetCHS_pT20.push_back(temp_GoodFatjetCHS_pT20);
            Vector_GoodFatjetCHS_pT20_Index.push_back(ifatjet);
        } 
       
        vector<float> Vector_mjCHS_pT20;   // mj
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjetCHS_pT20.size(); igoodfatjet++) 
        {
            float temp_mj_pT20 = Getmj(Vector_GoodFatjetCHS_pT20.at(igoodfatjet).Px(), 
                                       Vector_GoodFatjetCHS_pT20.at(igoodfatjet).Py(),
                                       Vector_GoodFatjetCHS_pT20.at(igoodfatjet).Pz(),
                                       Vector_GoodFatjetCHS_pT20.at(igoodfatjet).E());
            Vector_mjCHS_pT20.push_back(temp_mj_pT20);
        }
        MJCHS_pT20 = GetMJ(Vector_mjCHS_pT20);

        int NfatjetCHS_pT20 = Vector_GoodFatjetCHS_pT20.size(); 
        
        // 
        // pT(R=0.5) > 30 GeV
        // 
        
        // Regular jets
        double MJ_pT30=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjet_pT30; 
        vector<int> Vector_GoodFatjet_pT30_Index; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PFclean_R1p2_R0p5pT30_px->size(); ifatjet++) 
        {
            float temp_pT_pT30 = TMath::Sqrt(fastjets_AK5PFclean_R1p2_R0p5pT30_px->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT30_px->at(ifatjet)
                                        +fastjets_AK5PFclean_R1p2_R0p5pT30_py->at(ifatjet)*fastjets_AK5PFclean_R1p2_R0p5pT30_py->at(ifatjet));
            if(temp_pT_pT30<50) continue;
            TLorentzVector temp_GoodFatjet_pT30( fastjets_AK5PFclean_R1p2_R0p5pT30_px->at(ifatjet), 
                                            fastjets_AK5PFclean_R1p2_R0p5pT30_py->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT30_pz->at(ifatjet),
                                            fastjets_AK5PFclean_R1p2_R0p5pT30_energy->at(ifatjet));
            Vector_GoodFatjet_pT30.push_back(temp_GoodFatjet_pT30);
            Vector_GoodFatjet_pT30_Index.push_back(ifatjet);
        } 
       
        vector<float> Vector_mj_pT30;   // mj
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT30.size(); igoodfatjet++) 
        {
            float temp_mj_pT30 = Getmj(Vector_GoodFatjet_pT30.at(igoodfatjet).Px(), 
                                  Vector_GoodFatjet_pT30.at(igoodfatjet).Py(),
                                  Vector_GoodFatjet_pT30.at(igoodfatjet).Pz(),
                                  Vector_GoodFatjet_pT30.at(igoodfatjet).E());
            Vector_mj_pT30.push_back(temp_mj_pT30);
        }
        MJ_pT30 = GetMJ(Vector_mj_pT30);

        int Nfatjet_pT30 = Vector_GoodFatjet_pT30.size(); 
       
        // CHS jets
        double MJCHS_pT30=-999.; 
        // first, select good fat jets 
        vector<TLorentzVector> Vector_GoodFatjetCHS_pT30; 
        vector<int> Vector_GoodFatjetCHS_pT30_Index; 
        for(int ifatjet=0; ifatjet<(int)fastjets_AK5PF_R1p2_R0p5pT30_px->size(); ifatjet++) 
        {
            float temp_pT_pT30 = TMath::Sqrt(fastjets_AK5PF_R1p2_R0p5pT30_px->at(ifatjet)*fastjets_AK5PF_R1p2_R0p5pT30_px->at(ifatjet)
                                        +fastjets_AK5PF_R1p2_R0p5pT30_py->at(ifatjet)*fastjets_AK5PF_R1p2_R0p5pT30_py->at(ifatjet));
            if(temp_pT_pT30<50) continue;
            TLorentzVector temp_GoodFatjetCHS_pT30( fastjets_AK5PF_R1p2_R0p5pT30_px->at(ifatjet), 
                                            fastjets_AK5PF_R1p2_R0p5pT30_py->at(ifatjet),
                                            fastjets_AK5PF_R1p2_R0p5pT30_pz->at(ifatjet),
                                            fastjets_AK5PF_R1p2_R0p5pT30_energy->at(ifatjet));
            Vector_GoodFatjetCHS_pT30.push_back(temp_GoodFatjetCHS_pT30);
            Vector_GoodFatjetCHS_pT30_Index.push_back(ifatjet);
        } 
       
        vector<float> Vector_mjCHS_pT30;   // mj
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjetCHS_pT30.size(); igoodfatjet++) 
        {
            float temp_mj_pT30 = Getmj(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Px(), 
                                       Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Py(),
                                       Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Pz(),
                                       Vector_GoodFatjetCHS_pT30.at(igoodfatjet).E());
            Vector_mjCHS_pT30.push_back(temp_mj_pT30);
        }
        MJCHS_pT30 = GetMJ(Vector_mjCHS_pT30);

        int NfatjetCHS_pT30 = Vector_GoodFatjetCHS_pT30.size(); 
*/       
       
        //
        // Make fatjet 
        //
        vector<TLorentzVector> FatJetConstituent_pT20; 
        vector<TLorentzVector> FatJetConstituent_pT30; 
        vector<TLorentzVector> FatJetConstituentCHS_pT20; 
        vector<TLorentzVector> FatJetConstituentCHS_pT30; 
        for(unsigned int ijet=0;ijet<jets_AK5PFclean_pt->size();ijet++)
        {
            TLorentzVector tmp(jets_AK5PFclean_px->at(ijet),jets_AK5PFclean_py->at(ijet),
                               jets_AK5PFclean_pz->at(ijet), jets_AK5PFclean_energy->at(ijet));
            if(jets_AK5PFclean_pt->at(ijet)<20) continue;
            FatJetConstituent_pT20.push_back(tmp);
            if(jets_AK5PFclean_pt->at(ijet)<30) continue;
            FatJetConstituent_pT30.push_back(tmp);
        }
        for(unsigned int ijet=0;ijet<jets_AK5PF_pt->size();ijet++)
        {
            TLorentzVector tmp(jets_AK5PF_px->at(ijet), jets_AK5PF_py->at(ijet),
                               jets_AK5PF_pz->at(ijet), jets_AK5PF_energy->at(ijet));
            if(jets_AK5PF_pt->at(ijet)<20) continue;
            FatJetConstituentCHS_pT20.push_back(tmp);
            if(jets_AK5PF_pt->at(ijet)<30) continue;
            FatJetConstituentCHS_pT30.push_back(tmp);
        }

        // make vector of contituents adding RA4 leptons  
        double HT_dummy=-999.; 
        vector<int> LooseBJet_dummy; 
        vector<int> MediumBJet_dummy; 
        vector<int> GoodJets_AK5PFclean_PlusRA4Lep = GetJets(RA4Muon,RA4Elec,RA4MuonVeto,RA4ElecVeto,
                                                             HT_dummy,LooseBJet_dummy,MediumBJet_dummy, 
                                                             5, 30, 0.3); 
        
        vector<TLorentzVector> FatJetConstituent_pT30_PlusRA4Lep; 
        for(int igoodjet=0; igoodjet<(int)GoodJets_AK5PFclean_PlusRA4Lep.size(); igoodjet++) 
        {
            int ijet = GoodJets_AK5PFclean_PlusRA4Lep.at(igoodjet); 
            TLorentzVector tmp(jets_AK5PFclean_px->at(ijet), jets_AK5PFclean_py->at(ijet),
                               jets_AK5PFclean_pz->at(ijet), jets_AK5PFclean_energy->at(ijet));
            FatJetConstituent_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int igoodmus=0; igoodmus<RA4Muon.size(); igoodmus++) 
        {
            int imu = RA4Muon.at(igoodmus); 
            TLorentzVector tmp(mus_px->at(imu), mus_py->at(imu),
                               mus_pz->at(imu), mus_energy->at(imu));
            FatJetConstituent_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int ivetomus=0; ivetomus<RA4MuonVeto.size(); ivetomus++) 
        {
            int imu = RA4MuonVeto.at(ivetomus); 
            TLorentzVector tmp(mus_px->at(imu), mus_py->at(imu),
                               mus_pz->at(imu), mus_energy->at(imu));
            FatJetConstituent_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int igoodels=0; igoodels<RA4Elec.size(); igoodels++) 
        {
            int iel = RA4Elec.at(igoodels); 
            TLorentzVector tmp(els_px->at(iel), els_py->at(iel),
                               els_pz->at(iel), els_energy->at(iel));
            //cout << tmp.Px() << " " << tmp.Py() << " " << tmp.Pz() << " " << tmp.E() << " " << tmp.M() << endl; 
            FatJetConstituent_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int ivetoels=0; ivetoels<RA4ElecVeto.size(); ivetoels++) 
        {
            int iel = RA4ElecVeto.at(ivetoels); 
            TLorentzVector tmp(els_px->at(iel), els_py->at(iel),
                               els_pz->at(iel), els_energy->at(iel));
            FatJetConstituent_pT30_PlusRA4Lep.push_back(tmp);
        }

        //
        double HTCHS_dummy=-999.; 
        vector<int> LooseBJetCHS_dummy; 
        vector<int> MediumBJetCHS_dummy; 
        vector<int> GoodJets_AK5PF_PlusRA4Lep = GetJetsCHS(RA4Muon,RA4Elec,RA4MuonVeto,RA4ElecVeto,
                                                        HTCHS_dummy,LooseBJetCHS_dummy,MediumBJetCHS_dummy, 
                                                        5, 30, 0.3); 

        vector<TLorentzVector> FatJetConstituentCHS_pT30_PlusRA4Lep; 
        for(int igoodjet=0; igoodjet<(int)GoodJets_AK5PF_PlusRA4Lep.size(); igoodjet++) 
        {
            int ijet = GoodJets_AK5PF_PlusRA4Lep.at(igoodjet); 
            TLorentzVector tmp(jets_AK5PF_px->at(ijet), jets_AK5PF_py->at(ijet),
                               jets_AK5PF_pz->at(ijet), jets_AK5PF_energy->at(ijet));
            FatJetConstituentCHS_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int igoodmus=0; igoodmus<RA4Muon.size(); igoodmus++) 
        {
            int imu = RA4Muon.at(igoodmus); 
            TLorentzVector tmp(mus_px->at(imu), mus_py->at(imu),
                               mus_pz->at(imu), mus_energy->at(imu));
            FatJetConstituentCHS_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int ivetomus=0; ivetomus<RA4MuonVeto.size(); ivetomus++) 
        {
            int imu = RA4MuonVeto.at(ivetomus); 
            TLorentzVector tmp(mus_px->at(imu), mus_py->at(imu),
                               mus_pz->at(imu), mus_energy->at(imu));
            FatJetConstituentCHS_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int igoodels=0; igoodels<RA4Elec.size(); igoodels++) 
        {
            int iel = RA4Elec.at(igoodels); 
            TLorentzVector tmp(els_px->at(iel), els_py->at(iel),
                               els_pz->at(iel), els_energy->at(iel));
            //cout << tmp.Px() << " " << tmp.Py() << " " << tmp.Pz() << " " << tmp.E() << " " << tmp.M() << endl; 
            FatJetConstituentCHS_pT30_PlusRA4Lep.push_back(tmp);
        }
        for(int ivetoels=0; ivetoels<RA4ElecVeto.size(); ivetoels++) 
        {
            int iel = RA4ElecVeto.at(ivetoels); 
            TLorentzVector tmp(els_px->at(iel), els_py->at(iel),
                               els_pz->at(iel), els_energy->at(iel));
            FatJetConstituentCHS_pT30_PlusRA4Lep.push_back(tmp);
        }

        std::vector<FatJetPair> FatJet_R1p2_pT20_Eta5                   = makeFatJet(FatJetConstituent_pT20,                1.2, 20, 5);
        std::vector<FatJetPair> FatJet_R1p2_pT30_Eta5                   = makeFatJet(FatJetConstituent_pT30,                1.2, 30, 5);
        std::vector<FatJetPair> FatJetCHS_R1p2_pT20_Eta5                = makeFatJet(FatJetConstituentCHS_pT20,             1.2, 20, 5);
        std::vector<FatJetPair> FatJetCHS_R1p2_pT30_Eta5                = makeFatJet(FatJetConstituentCHS_pT30,             1.2, 30, 5);
        std::vector<FatJetPair> FatJet_R1p2_pT30_Eta5_PlusRA4Lep        = makeFatJet(FatJetConstituent_pT30_PlusRA4Lep,     1.2, 30, 5);
        std::vector<FatJetPair> FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep     = makeFatJet(FatJetConstituentCHS_pT30_PlusRA4Lep,  1.2, 30, 5);

        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT20_Eta5.size();ifj++)
        {
            if(FatJet_R1p2_pT20_Eta5.at(ifj).first.pt()<50) continue;
            mj_pT20_.push_back(FatJet_R1p2_pT20_Eta5.at(ifj).first.m());
            FatjetPt_pT20_.push_back(FatJet_R1p2_pT20_Eta5.at(ifj).first.pt());
            FatjetEta_pT20_.push_back(FatJet_R1p2_pT20_Eta5.at(ifj).first.eta());
            FatjetPhi_pT20_.push_back(FatJet_R1p2_pT20_Eta5.at(ifj).first.phi()<3.141592?FatJet_R1p2_pT20_Eta5.at(ifj).first.phi():FatJet_R1p2_pT20_Eta5.at(ifj).first.phi()-2*3.141592);
            FatjetN_pT20_.push_back(FatJet_R1p2_pT20_Eta5.at(ifj).second);
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT30_Eta5.size();ifj++)
        {
            if(FatJet_R1p2_pT30_Eta5.at(ifj).first.pt()<50) continue;
            mj_pT30_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).first.m());
            FatjetPt_pT30_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).first.pt());
            FatjetEta_pT30_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).first.eta());
            FatjetPhi_pT30_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).first.phi()<3.141592?FatJet_R1p2_pT30_Eta5.at(ifj).first.phi():FatJet_R1p2_pT30_Eta5.at(ifj).first.phi()-2*3.141592);
            FatjetN_pT30_.push_back(FatJet_R1p2_pT30_Eta5.at(ifj).second);
        }
        for(unsigned int ifj=0; ifj<FatJetCHS_R1p2_pT20_Eta5.size();ifj++)
        {
            if(FatJetCHS_R1p2_pT20_Eta5.at(ifj).first.pt()<50) continue;
            mjCHS_pT20_.push_back(FatJetCHS_R1p2_pT20_Eta5.at(ifj).first.m());
            FatjetCHSPt_pT20_.push_back(FatJetCHS_R1p2_pT20_Eta5.at(ifj).first.pt());
            FatjetCHSEta_pT20_.push_back(FatJetCHS_R1p2_pT20_Eta5.at(ifj).first.eta());
            FatjetCHSPhi_pT20_.push_back(FatJetCHS_R1p2_pT20_Eta5.at(ifj).first.phi()<3.141592?FatJetCHS_R1p2_pT20_Eta5.at(ifj).first.phi():FatJetCHS_R1p2_pT20_Eta5.at(ifj).first.phi()-2*3.141592);
            FatjetCHSN_pT20_.push_back(FatJetCHS_R1p2_pT20_Eta5.at(ifj).second);
        }
        for(unsigned int ifj=0; ifj<FatJetCHS_R1p2_pT30_Eta5.size();ifj++)
        {
            if(FatJetCHS_R1p2_pT30_Eta5.at(ifj).first.pt()<50) continue;
            mjCHS_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5.at(ifj).first.m());
            FatjetCHSPt_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5.at(ifj).first.pt());
            FatjetCHSEta_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5.at(ifj).first.eta());
            FatjetCHSPhi_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5.at(ifj).first.phi()<3.141592?FatJetCHS_R1p2_pT30_Eta5.at(ifj).first.phi():FatJetCHS_R1p2_pT30_Eta5.at(ifj).first.phi()-2*3.141592);
            FatjetCHSN_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5.at(ifj).second);
        }
        for(unsigned int ifj=0; ifj<FatJet_R1p2_pT30_Eta5_PlusRA4Lep.size();ifj++)
        {
            if(FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.pt()<50) continue;
            mjLep_pT30_.push_back(FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.m());
            FatjetLepPt_pT30_.push_back(FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.pt());
            FatjetLepEta_pT30_.push_back(FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.eta());
            FatjetLepPhi_pT30_.push_back(FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.phi()<3.141592?FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.phi():FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.phi()-2*3.141592);
            FatjetLepN_pT30_.push_back(FatJet_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).second);
        }
        for(unsigned int ifj=0; ifj<FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.size();ifj++)
        {
            if(FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.pt()<50) continue;
            mjCHSLep_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.m());
            FatjetCHSLepPt_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.pt());
            FatjetCHSLepEta_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.eta());
            FatjetCHSLepPhi_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.phi()<3.141592?FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.phi():FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).first.phi()-2*3.141592);
            FatjetCHSLepN_pT30_.push_back(FatJetCHS_R1p2_pT30_Eta5_PlusRA4Lep.at(ifj).second);
        }

        //
        // Fill the baby variables 
        //
        run_                =   run;
        lumiblock_          =   lumiblock;
        event_              =   event; 
        TrigMuon_           =   PassMuonTrig(); 
        TrigHTMuon_         =   PassHTMuonTrig(); 
        TrigSingleMuon_     =   PassSingleMuonTrig(); 
        TrigDiMuon_         =   PassDiMuonTrig(); 
        TrigElectron_       =   PassElecTrig(); 
        TrigHTElectron_     =   PassHTElecTrig(); 
        TrigSingleElectron_ =   PassSingleElecTrig(); 
        TrigDiElectron_     =   PassDiElecTrig(); 
        TrigMuEG_           =   PassMuEGTrig(); 
        Nfatjet_pT20_       =   mj_pT20_.size();
        Nfatjet_pT30_       =   mj_pT30_.size();
        NfatjetCHS_pT20_    =   mjCHS_pT20_.size();
        NfatjetCHS_pT30_    =   mjCHS_pT30_.size();
        Nskinnyjet_         =   GoodJets_AK5PFclean.size();
        NBtagCSVM_          =   MediumBJet.size();
        NskinnyjetCHS_      =   GoodJets_AK5PF.size();
        NBtagCHSCSVM_       =   MediumBJetCHS.size();
        Npv_                =   Npv;
        if(!isData) Npu_    =   PU_TrueNumInteractions->at(1);
        EventWeight_        =   EventWeight;
        MJ_pT20_            =   GetMJ(mj_pT20_);
        MJ_pT30_            =   GetMJ(mj_pT30_);
        MJCHS_pT20_         =   GetMJ(mjCHS_pT20_);
        MJCHS_pT30_         =   GetMJ(mjCHS_pT30_);
        MET_                =   pfTypeImets_et->at(0);
        METPhi_             =   pfTypeImets_phi->at(0);
        HT_                 =   HT;
        HTCHS_              =   HTCHS;
       /* 
        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT20.size(); igoodfatjet++) 
        { 
            mj_pT20_.push_back(Vector_mj_pT20.at(igoodfatjet));
            FatjetPt_pT20_.push_back(Vector_GoodFatjet_pT20.at(igoodfatjet).Pt());
            FatjetEta_pT20_.push_back(Vector_GoodFatjet_pT20.at(igoodfatjet).Eta());
            FatjetPhi_pT20_.push_back(Vector_GoodFatjet_pT20.at(igoodfatjet).Phi());
            FatjetN_pT20_.push_back(fastjets_AK5PFclean_R1p2_R0p5pT20_nconstituents->at(Vector_GoodFatjet_pT20_Index.at(igoodfatjet)));
        }

        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjetCHS_pT20.size(); igoodfatjet++) 
        { 
            mjCHS_pT20_.push_back(Vector_mjCHS_pT20.at(igoodfatjet));
            FatjetCHSPt_pT20_.push_back(Vector_GoodFatjetCHS_pT20.at(igoodfatjet).Pt());
            FatjetCHSEta_pT20_.push_back(Vector_GoodFatjetCHS_pT20.at(igoodfatjet).Eta());
            FatjetCHSPhi_pT20_.push_back(Vector_GoodFatjetCHS_pT20.at(igoodfatjet).Phi());
            FatjetCHSN_pT20_.push_back(fastjets_AK5PF_R1p2_R0p5pT20_nconstituents->at(Vector_GoodFatjetCHS_pT20_Index.at(igoodfatjet)));
        }

        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjet_pT30.size(); igoodfatjet++) 
        { 
            mj_pT30_.push_back(Vector_mj_pT30.at(igoodfatjet));
            FatjetPt_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Pt());
            FatjetEta_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Eta());
            FatjetPhi_pT30_.push_back(Vector_GoodFatjet_pT30.at(igoodfatjet).Phi());
            FatjetN_pT30_.push_back(fastjets_AK5PFclean_R1p2_R0p5pT30_nconstituents->at(Vector_GoodFatjet_pT30_Index.at(igoodfatjet)));
        }

        for(int igoodfatjet=0; igoodfatjet<(int)Vector_GoodFatjetCHS_pT30.size(); igoodfatjet++) 
        { 
            mjCHS_pT30_.push_back(Vector_mjCHS_pT30.at(igoodfatjet));
            FatjetCHSPt_pT30_.push_back(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Pt());
            FatjetCHSEta_pT30_.push_back(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Eta());
            FatjetCHSPhi_pT30_.push_back(Vector_GoodFatjetCHS_pT30.at(igoodfatjet).Phi());
            FatjetCHSN_pT30_.push_back(fastjets_AK5PF_R1p2_R0p5pT30_nconstituents->at(Vector_GoodFatjetCHS_pT30_Index.at(igoodfatjet)));
        }
       */ 
        for(unsigned int imus=0; imus<RA4Muon.size(); imus++) 
        {
            RA4MusPt_.push_back(mus_pt->at(RA4Muon.at(imus)));
            RA4MusEta_.push_back(mus_eta->at(RA4Muon.at(imus)));
            RA4MusPhi_.push_back(mus_phi->at(RA4Muon.at(imus)));
            RA4MusE_.push_back(mus_energy->at(RA4Muon.at(imus)));
            RA4MusQ_.push_back(mus_charge->at(RA4Muon.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4Elec.size(); iels++) 
        {
            RA4ElsPt_.push_back(els_pt->at(RA4Elec.at(iels)));
            RA4ElsEta_.push_back(els_eta->at(RA4Elec.at(iels)));
            RA4ElsPhi_.push_back(els_phi->at(RA4Elec.at(iels)));
            RA4ElsE_.push_back(els_energy->at(RA4Elec.at(iels)));
            RA4ElsQ_.push_back(els_charge->at(RA4Elec.at(iels)));
        }
        for(unsigned int imus=0; imus<RA4MuonVeto.size(); imus++) 
        {
            RA4MusVetoPt_.push_back(mus_pt->at(RA4MuonVeto.at(imus)));
            RA4MusVetoEta_.push_back(mus_eta->at(RA4MuonVeto.at(imus)));
            RA4MusVetoPhi_.push_back(mus_phi->at(RA4MuonVeto.at(imus)));
        }
        for(unsigned int iels=0; iels<RA4ElecVeto.size(); iels++) 
        {
            RA4ElsVetoPt_.push_back(els_pt->at(RA4ElecVeto.at(iels)));
            RA4ElsVetoEta_.push_back(els_eta->at(RA4ElecVeto.at(iels)));
            RA4ElsVetoPhi_.push_back(els_phi->at(RA4ElecVeto.at(iels)));
        }
        
        for(int igoodjet=0; igoodjet<(int)GoodJets_AK5PFclean.size(); igoodjet++) 
        {
          int ijet = GoodJets_AK5PFclean.at(igoodjet); 
          JetPt_.push_back(jets_AK5PFclean_pt->at(ijet)); 
          JetEta_.push_back(jets_AK5PFclean_eta->at(ijet)); 
          JetPhi_.push_back(jets_AK5PFclean_phi->at(ijet)); 
          JetE_.push_back(jets_AK5PFclean_energy->at(ijet)); 
          JetCSV_.push_back(jets_AK5PFclean_btag_secVertexCombined->at(ijet)); 
          JetMCId_.push_back(jets_AK5PFclean_partonFlavour->at(ijet)); 
        }
        
        for(int igoodjet=0; igoodjet<(int)GoodJets_AK5PF.size(); igoodjet++) 
        {
          int ijet = GoodJets_AK5PF.at(igoodjet); 
          JetCHSPt_.push_back(jets_AK5PF_pt->at(ijet)); 
          JetCHSEta_.push_back(jets_AK5PF_eta->at(ijet)); 
          JetCHSPhi_.push_back(jets_AK5PF_phi->at(ijet)); 
          JetCHSE_.push_back(jets_AK5PF_energy->at(ijet)); 
          JetCHSCSV_.push_back(jets_AK5PF_btag_secVertexCombined->at(ijet)); 
          JetCHSMCId_.push_back(jets_AK5PF_partonFlavour->at(ijet)); 
          JetCHSNCh_.push_back(jets_AK5PF_chg_Mult->at(ijet)); 
          JetCHSNNeu_.push_back(jets_AK5PF_neutral_Mult->at(ijet)); 
        }

        // gen top pT for TTbar samples  
        for(int igen=0; igen<mc_doc_id->size(); igen++) 
        { 
            GenPt_.push_back(   mc_doc_pt->at(igen));  
            GenPhi_.push_back(  mc_doc_phi->at(igen));  
            GenE_.push_back(  mc_doc_energy->at(igen));  
            GenEta_.push_back(  mc_doc_eta->at(igen));  
            GenId_.push_back(   mc_doc_id->at(igen));  
            GenMId_.push_back(  mc_doc_mother_id->at(igen));  
            GenGMId_.push_back( mc_doc_grandmother_id->at(igen));  

            if(ProcessName.Contains("TT")) 
            {
                if(mc_doc_id->at(igen)==6)  
                {
                    top1pT_ = mc_doc_pt->at(igen);  
                    top1Phi_ = mc_doc_phi->at(igen);  
                    top1Eta_ = mc_doc_eta->at(igen);  
                }
                if(mc_doc_id->at(igen)==-6)
                {
                    top2pT_ = mc_doc_pt->at(igen);  
                    top2Phi_ = mc_doc_phi->at(igen);  
                    top2Eta_ = mc_doc_eta->at(igen);  
                }
            }
        }

        babyTree_->Fill(); // Fill all events

        //for(int i=0; i<GoodJets_AK5PFclean.size(); i++) cout << event << " :: " << HT << endl;
        
        // Clean fat jets? 
        // (1) identify if all skinny jets associated with a given fat jet 
        //     are good jets by comparing index of skinny jet in GoodJets_AK5PFclean
        //     and fastjet_.._R1p2pTxx_index
        // (2) if a skinny jet is not in the good jet, its four vector is subtracted 
        //     from the fat jet and mj is calculated again
        // (3) MJ is then calculated


    } // event loop
    cout << endl;
    cout << "[MJ Analysis] Looping over events has been done" << endl;
    
    // clean up  
    fPUFile->Close();

    //
    // Write the baby file 
    //
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();

    // 
    //  
    // 
    TDatime DTEnd;
    int EndTime = DTEnd.GetTime(); 
    cout << "[MJ Analysis] End time : " << (EndTime/10000)%100 << ":"
         << (EndTime/100)%100 << ":" << EndTime%100
         << endl; 
    cout << "[MJ Analysis] Done with " << ProcessName << endl; 
    cout << "[MJ Analysis] ------------------------------------------------------------------------------------"<<endl; 
   
    // cleanup
    delete chainA;
    delete chainB;
}
