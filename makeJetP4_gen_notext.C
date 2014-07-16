// ------------------------------------
// Author : Jae Hyeok Yoo
//          jaehyeok@hep.ucsb.edu
// ------------------------------------
//
// * This code is composed of 3 blocks
//  
// Block 1 : Read information about PF candidates from a cfA ntuple   
//           and write that information on a text file  
// Block 2 : Read the text file from (1) and run Fastjet package. 
//           The information of the recontructed jets will be 
//           written on a text file 
// Block 3 : Read the text file from (2) and make new branches in 
//           the input cfA in (1) to store the information 
// 
// * Pre-requisites
//
// (1) Download(http://fastjet.fr) and set up Fastjet package
//    * version 3.0.6(http://fastjet.fr/repo/fastjet-3.0.6.tar.gz) 
//    * manual : http://fastjet.fr/repo/fastjet-doc-3.0.6.pdf
//      --> Chapter 2 has an instruction to set up the code
// (2) Turn off some of the printouts so that fastjet_example 
//     prints out only information of PF candidates, i.e., 
//     no table header, ... 
// 
// * To run this code, do
//   
//      $root -q -b makeJetP4.C++ 
// 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream> // for stringstream

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"

//RH include necessary fastjet files
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace std;

const bool DEBUG = false;

void makeJetP4_gen_notext(TString InRootFile, double Rparam=1.2, int R0p5JetpTcut=30) { 
    
    cout << " ................................................................." << endl; 
    cout << " ... Processing file = " << InRootFile << endl;
    cout << " ... Rparam = "          << Rparam << endl;
    cout << " ... pT(R=0.5) > "       << R0p5JetpTcut << " GeV" <<  endl;
    cout << " ................................................................." << endl; 
    
    // 
    // Get tree from a cfA ntuple 
    // 
    TFile *f = TFile::Open(InRootFile, "UPDATE");
    TDirectory* dir = f->GetDirectory("configurableAnalysis");
    dir->cd();
    TTree *eventB = (TTree*)dir->Get("eventB");
    
    // 
    // set address of variables to read 
    // 
    UInt_t   event_ = 0;
    eventB->SetBranchAddress("event", &event_);
    UInt_t   run_ = 0;
    eventB->SetBranchAddress("run", &run_);
    UInt_t   lumiblock_ = 0;
    eventB->SetBranchAddress("lumiblock", &lumiblock_);
    vector<float>   *FatjetConstituent_pt_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_gen_pt", &FatjetConstituent_pt_);
    vector<float>   *FatjetConstituent_energy_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_gen_Energy", &FatjetConstituent_energy_); // Note capital E 
    vector<float>   *FatjetConstituent_phi_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_gen_phi", &FatjetConstituent_phi_);
    vector<float>   *FatjetConstituent_eta_ = 0;
    eventB->SetBranchAddress("jets_AK5PFclean_gen_eta", &FatjetConstituent_eta_);
    
    // 
    // Define new variables to write 
    // 
    // need this line because ROOT does not know vector<float>
    // one needs to generate dictionary using rootcint or 
    // just include this line 
    // [ref] http://root.cern.ch/phpBB3/viewtopic.php?t=8467
    gROOT->ProcessLine("#include <vector>"); 
    vector<float>   *Fatjet_px = 0;
    TBranch *pxb =  eventB->Branch(Form("fastjets_AK5PFclean_gen_R1p2_R0p5pT%i_px",R0p5JetpTcut), &Fatjet_px);
    vector<float>   *Fatjet_py = 0;
    TBranch *pyb = eventB->Branch(Form("fastjets_AK5PFclean_gen_R1p2_R0p5pT%i_py",R0p5JetpTcut), &Fatjet_py);
    vector<float>   *Fatjet_pz = 0;
    TBranch *pzb = eventB->Branch(Form("fastjets_AK5PFclean_gen_R1p2_R0p5pT%i_pz",R0p5JetpTcut), &Fatjet_pz);
    vector<float>   *Fatjet_energy = 0;
    TBranch *energyb = eventB->Branch(Form("fastjets_AK5PFclean_gen_R1p2_R0p5pT%i_energy",R0p5JetpTcut), &Fatjet_energy);
    vector<float>   *Fatjet_phi = 0;
    TBranch *phib = eventB->Branch(Form("fastjets_AK5PFclean_gen_R1p2_R0p5pT%i_phi",R0p5JetpTcut), &Fatjet_phi);
    vector<float>   *Fatjet_eta = 0;
    TBranch *etab = eventB->Branch(Form("fastjets_AK5PFclean_gen_R1p2_R0p5pT%i_eta",R0p5JetpTcut), &Fatjet_eta);

    // 
    // Histgrom : to draw eta-phi plot of energy deposit 
    //  (1) Bin size is 0.087x0.087 to mimic the size of hcal tower
    //  (2) Bin Entry is the sum over energies of PF candidates in a given bin  
    // 
    TH2F *h2 = new TH2F("h2","h2", 115, -5.0, 5.0, 72, -3.141592, 3.141592);
    
    // 
    // Loop over entries
    // 
    Int_t nentries = (Int_t)eventB->GetEntries();
    if(DEBUG) nentries = 10;
    cout<<"The number of entries is: "<<nentries<<endl;

    // main event loop
    for(int ib = 0; ib<nentries; ib++) {
      
        // Counting to see progress
        //if(ib%100==0) cout << " ...... " << ib << " events processed "<< endl; 

        // get the entry of event ib
        eventB->GetEntry(ib);
        
        // ---------------------------------------
        //          Block 1 
        // ---------------------------------------

        // Loop over R=0.5 jets, form into PseudoJets vector
        vector<fastjet::PseudoJet> input_particles;
        double FatjetConstituent_px_tmp, FatjetConstituent_py_tmp, FatjetConstituent_pz_tmp, FatjetConstituent_energy_tmp;

        for(int ijet = 0; ijet < (int)FatjetConstituent_pt_->size(); ijet++) { 
            // px: pt*cos phi
            // py: pt*sin phi
            // pz :jets_AK5PFclean_pt/TMath::Tan(2*TMath::ATan(TMath::Exp(-jets_AK5PFclean_eta)))
            FatjetConstituent_px_tmp = FatjetConstituent_pt_->at(ijet)*TMath::Cos(FatjetConstituent_phi_->at(ijet));
            FatjetConstituent_py_tmp = FatjetConstituent_pt_->at(ijet)*TMath::Sin(FatjetConstituent_phi_->at(ijet));
            FatjetConstituent_pz_tmp = FatjetConstituent_pt_->at(ijet)/TMath::Tan(2*TMath::ATan(TMath::Exp(-FatjetConstituent_eta_->at(ijet))));
            FatjetConstituent_energy_tmp = FatjetConstituent_energy_->at(ijet);	  
           
            if(TMath::Sqrt( FatjetConstituent_px_tmp*FatjetConstituent_px_tmp
                           +FatjetConstituent_py_tmp*FatjetConstituent_py_tmp) < R0p5JetpTcut) continue;


            input_particles.push_back(fastjet::PseudoJet(FatjetConstituent_px_tmp,FatjetConstituent_py_tmp,
                                                         FatjetConstituent_pz_tmp,FatjetConstituent_energy_tmp));
            h2->Fill( FatjetConstituent_eta_->at(ijet), FatjetConstituent_phi_->at(ijet), FatjetConstituent_energy_->at(ijet));
        }


        // ---------------------------------------
        //          Block 2 
        // ---------------------------------------

        //
        // Run Fastjet to reconstuct jets 
        //

        // RH create an object that represents your choice of jet algorithm and 
        // the associated parameters
        fastjet::Strategy strategy = fastjet::Best;
        fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, recomb_scheme, strategy);

        // run the jet clustering with the above jet definition
        fastjet::ClusterSequence clust_seq(input_particles, jet_def);

        // ---------------------------------------
        //          Block 3 
        // ---------------------------------------

        // 
        // Get p4 of the reconstructed jets  
        //
        //RH
        double ptmin = 0.0; // could use 3.0 here, instead of applying later 
        vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
        //Sort by pt
        vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(inclusive_jets);
        //fill fastjet output into vectors, continue as original code
        for(int isortjets = 0; isortjets< (int)sorted_jets.size(); isortjets++){
            //store only if pt >3 GeV to match CMS jets
            if(TMath::Sqrt( sorted_jets[isortjets].px()*sorted_jets[isortjets].px()
                           +sorted_jets[isortjets].py()*sorted_jets[isortjets].py())>(DEBUG?30:3)) {
                Fatjet_px->push_back(sorted_jets[isortjets].px());
                Fatjet_py->push_back(sorted_jets[isortjets].py());
                Fatjet_pz->push_back(sorted_jets[isortjets].pz());
                Fatjet_energy->push_back(sorted_jets[isortjets].E());
                Fatjet_eta->push_back(sorted_jets[isortjets].rap());
                Fatjet_phi->push_back( sorted_jets[isortjets].phi()<3.141592?sorted_jets[isortjets].phi():sorted_jets[isortjets].phi()-2*3.141592); 
        
                if(DEBUG) {
                    cout << event_ << " " 
                        << TMath::Sqrt( sorted_jets[isortjets].px()*sorted_jets[isortjets].px()
                                       +sorted_jets[isortjets].py()*sorted_jets[isortjets].py()) << " " 
                        << sorted_jets[isortjets].rap() << " " 
                        //<< sorted_jets[isortjets].phi()<3.141592 ? sorted_jets[isortjets].phi() : (sorted_jets[isortjets].phi()-2*3.141592) << " "  
                        << endl; 
                }
            }
        }
        // Fill the branches 
        pxb->Fill();
        pyb->Fill();
        pzb->Fill();
        energyb->Fill();
        phib->Fill();
        etab->Fill();
        
        // 
        // Draw a lego plot (eta, phi) 
        //
        if(DEBUG) {
            TCanvas *c = new TCanvas();
            c->cd(1);
            h2->Draw("colz");
            h2->SetTitle(Form("run=%i lumi=%i event=%i R=%.1f", run_, lumiblock_, event_, Rparam));
            h2->SetMaximum(50); 
            h2->SetStats(0); 
            h2->SetXTitle("#eta"); 
            h2->SetYTitle("#phi"); 

            //Draw circles around jets
            TEllipse *cone[Fatjet_eta->size()]; 
            for(int ijets=0; ijets<(int)Fatjet_eta->size(); ijets++){
                cone[ijets] = new TEllipse(Fatjet_eta->at(ijets), Fatjet_phi->at(ijets), Rparam, Rparam);
                cone[ijets]->SetFillStyle(3003);
                cone[ijets]->SetFillColor(kYellow);
                cone[ijets]->Draw();
            }

            c->SaveAs(Form("EtaPhiViewPFCand_Run%i_Lumi%i_Event%i_R%.1f_usingR0p5.pdf", 
                        run_, lumiblock_, event_, Rparam));
            h2->Reset(); 
            for(int ijets=0; ijets<(int)Fatjet_eta->size(); ijets++) delete cone[ijets];

        } 
   
        // Clear vectors for the next event 
        Fatjet_px->clear();
        Fatjet_py->clear();
        Fatjet_pz->clear();
        Fatjet_energy->clear();
        Fatjet_eta->clear();
        Fatjet_phi->clear();
        
    } // event loop

    // update the tree and close file
    if(!DEBUG) eventB->Write();
    f->Close();
    
    //
    // cleanup
    //
    delete f;
}
