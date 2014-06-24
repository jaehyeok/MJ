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

using namespace std;

const bool DEBUG = false;
const int R0p5JetpTcut = 30;

void makeJetP4(TString InRootFile, double Rparam=0.5) { 
    
    cout << " ... Processing file = " << InRootFile << endl;
    cout << " ... Rparam = " << Rparam << endl;
    cout << " ................................................................." << endl; 
    
    // 
    // Get tree from a cfA ntuple 
    // 
    TFile *f = TFile::Open("../cfA/"+InRootFile, "UPDATE");
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
    vector<float>   *fastjets_AK5PF_px_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_px", &fastjets_AK5PF_px_);
    vector<float>   *fastjets_AK5PF_py_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_py", &fastjets_AK5PF_py_);
    vector<float>   *fastjets_AK5PF_pz_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_pz", &fastjets_AK5PF_pz_);
    vector<float>   *fastjets_AK5PF_energy_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_energy", &fastjets_AK5PF_energy_);
    vector<float>   *fastjets_AK5PF_phi_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_phi", &fastjets_AK5PF_phi_);
    vector<float>   *fastjets_AK5PF_eta_ = 0;
    eventB->SetBranchAddress("fastjets_AK5PF_eta", &fastjets_AK5PF_eta_);
    
    // 
    // Define new variables to write 
    // 
    // need this line because ROOT does not know vector<float>
    // one needs to generate dictionary using rootcint or 
    // just include this line 
    // [ref] http://root.cern.ch/phpBB3/viewtopic.php?t=8467
    gROOT->ProcessLine("#include <vector>"); 
    vector<float>   *fastjets_AK5PF_R1p2_px = 0;
    TBranch *pxb =  eventB->Branch(Form("fastjets_AK5PF_R1p2_pT%i_px",R0p5JetpTcut), &fastjets_AK5PF_R1p2_px);
    vector<float>   *fastjets_AK5PF_R1p2_py = 0;
    TBranch *pyb = eventB->Branch(Form("fastjets_AK5PF_R1p2_pT%i_py",R0p5JetpTcut), &fastjets_AK5PF_R1p2_py);
    vector<float>   *fastjets_AK5PF_R1p2_pz = 0;
    TBranch *pzb = eventB->Branch(Form("fastjets_AK5PF_R1p2_pT%i_pz",R0p5JetpTcut), &fastjets_AK5PF_R1p2_pz);
    vector<float>   *fastjets_AK5PF_R1p2_energy = 0;
    TBranch *energyb = eventB->Branch(Form("fastjets_AK5PF_R1p2_pT%i_energy",R0p5JetpTcut), &fastjets_AK5PF_R1p2_energy);
    vector<float>   *fastjets_AK5PF_R1p2_phi = 0;
    TBranch *phib = eventB->Branch(Form("fastjets_AK5PF_R1p2_pT%i_phi",R0p5JetpTcut), &fastjets_AK5PF_R1p2_phi);
    vector<float>   *fastjets_AK5PF_R1p2_eta = 0;
    TBranch *etab = eventB->Branch(Form("fastjets_AK5PF_R1p2_pT%i_eta",R0p5JetpTcut), &fastjets_AK5PF_R1p2_eta);

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
        if(ib%100==0) cout << " ...... " << ib << " events processed "<< endl; 

        // get the entry of event ib
        eventB->GetEntry(ib);
        
        // ---------------------------------------
        //          Block 1 
        // ---------------------------------------
        
        // Open a text file to write information about PF candidates  
        ofstream fout;
        fout.open(Form("OneEvent_FastjetsR0p5_tmp_%i.dat", ib));
        
        // loop over Fastjets with R=0.5
        for(int ifastjets = 0; ifastjets < (int)fastjets_AK5PF_px_->size(); ifastjets++) {
            
            if(TMath::Sqrt( fastjets_AK5PF_px_->at(ifastjets)*fastjets_AK5PF_px_->at(ifastjets) 
                           +fastjets_AK5PF_py_->at(ifastjets)*fastjets_AK5PF_py_->at(ifastjets)) < R0p5JetpTcut) continue; 

            fout.width(15); fout << fastjets_AK5PF_px_->at(ifastjets) << "\t";
            fout.width(15); fout << fastjets_AK5PF_py_->at(ifastjets) << "\t";
            fout.width(15); fout << fastjets_AK5PF_pz_->at(ifastjets) << "\t";
            fout.width(15); fout << fastjets_AK5PF_energy_->at(ifastjets) << endl;
           
            float pt = TMath::Sqrt(fastjets_AK5PF_px_->at(ifastjets)*fastjets_AK5PF_px_->at(ifastjets) 
                                  +fastjets_AK5PF_py_->at(ifastjets)*fastjets_AK5PF_py_->at(ifastjets));
            h2->Fill( fastjets_AK5PF_eta_->at(ifastjets), fastjets_AK5PF_phi_->at(ifastjets), pt);
        }
        fout.close();




        // ---------------------------------------
        //          Block 2 
        // ---------------------------------------

        //
        // Run Fastjet to reconstuct jets 
        //
        // Using the example code "fastjet_example.cc" 
        // It takes OneEvent_PFCands_tmp_%i.dat as an input 
        // and writes p4 of reconstructed jets in  OneEvent_PFJets_tmp_%i.dat. 
        // For the details about the code, take a look at fastjet_example.cc.
        gSystem->Exec(Form("../fastjet_example %f < OneEvent_FastjetsR0p5_tmp_%i.dat > OneEvent_FastjetsR1p2_tmp_%i.dat", Rparam, ib, ib));

        // ---------------------------------------
        //          Block 3 
        // ---------------------------------------

        // 
        // Get p4 of the reconstructed jets  
        //
        string line;
        double id, px, py , pz, energy, eta, phi, ncand;
        ifstream fin(Form("OneEvent_FastjetsR1p2_tmp_%i.dat", ib)); 
        if(fin.is_open()) { 
            while(fin.good()){

                // get a line from fin
                getline(fin, line);
                if(line=="") break; // need this to avoid line without entry 
                // Store each element in the line to the defined variables
                stringstream stream(line);
                stream >> id >> px >> py >> pz >> energy >> eta >> phi >> ncand; 
               
                // store only when pT > 3 GeV 
                // (same as CMS jet reconstruction cut)
                if(TMath::Sqrt(px*px+py*py)>(DEBUG?30:3)) {
                    fastjets_AK5PF_R1p2_px->push_back(px);
                    fastjets_AK5PF_R1p2_py->push_back(py);
                    fastjets_AK5PF_R1p2_pz->push_back(pz);
                    fastjets_AK5PF_R1p2_energy->push_back(energy);
                    fastjets_AK5PF_R1p2_eta->push_back(eta);
                    fastjets_AK5PF_R1p2_phi->push_back(phi); 

                    if(DEBUG) {
                        cout << event_ << " " 
                             << id << " " 
                             << TMath::Sqrt(px*px+py*py) << " " 
                             << eta << " " 
                             << phi << " "  
                             << endl; 
                    }
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
            TEllipse *cone[fastjets_AK5PF_R1p2_eta->size()]; 
            for(int ijets=0; ijets<(int)fastjets_AK5PF_R1p2_eta->size(); ijets++){
                cone[ijets] = new TEllipse(fastjets_AK5PF_R1p2_eta->at(ijets), fastjets_AK5PF_R1p2_phi->at(ijets), Rparam, Rparam);
                cone[ijets]->SetFillStyle(3003);
                cone[ijets]->SetFillColor(kYellow);
                cone[ijets]->Draw();
            }

            c->SaveAs(Form("EtaPhiViewPFCand_Run%i_Lumi%i_Event%i_R%.1f_usingR0p5.pdf", 
                        run_, lumiblock_, event_, Rparam));
            h2->Reset(); 
            for(int ijets=0; ijets<(int)fastjets_AK5PF_R1p2_eta->size(); ijets++) delete cone[ijets];

        } 
   
        // Clear vectors for the next event 
        fastjets_AK5PF_R1p2_px->clear();
        fastjets_AK5PF_R1p2_py->clear();
        fastjets_AK5PF_R1p2_pz->clear();
        fastjets_AK5PF_R1p2_energy->clear();
        fastjets_AK5PF_R1p2_eta->clear();
        fastjets_AK5PF_R1p2_phi->clear();
        
        // remove text files
        if(!DEBUG) gSystem->Exec("'rm' OneEvent_Fastjets*.dat");
    
    } // event loop

    // update the tree and close file
    if(!DEBUG) eventB->Write();
    f->Close();
    
    //
    // cleanup
    //
    delete f;
}
