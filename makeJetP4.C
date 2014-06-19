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
//double Rparam = 1.2;        // Distance parameter for jet clustering

void makeJetP4(TString InRootFile, double Rparam=0.5) { 
    
    cout << " ... Processing file = " << InRootFile << endl;
    cout << " ... Rparam = " << Rparam << endl;
    cout << " ................................................................." << endl; 
    
    // 
    // Get tree from a cfA ntuple 
    // 
    //TFile *f = TFile::Open("../cfA/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2029_v71_f1_1_YWH_Test.root", "UPDATE");
    //TFile *f = TFile::Open("../cfA/cfA_JetHT_Run2012D-PromptReco-v1_AOD_UCSB2030_v71_f1_1_fip_Test.root", "UPDATE");
    //TFile *f = TFile::Open("../cfA/NopT10Cut/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB2048_v71_f1_1_DGh_Test.root", "UPDATE");
    //TFile *f = TFile::Open("../cfA/NopT10Cut/cfA_JetHT_Run2012D-PromptReco-v1_AOD_UCSB2049_v71_f1000_1_473_Test.root", "UPDATE");
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
    vector<float>   *pfcand_px_ = 0;
    eventB->SetBranchAddress("pfcand_px", &pfcand_px_);
    vector<float>   *pfcand_py_ = 0;
    eventB->SetBranchAddress("pfcand_py", &pfcand_py_);
    vector<float>   *pfcand_pz_ = 0;
    eventB->SetBranchAddress("pfcand_pz", &pfcand_pz_);
    vector<float>   *pfcand_energy_ = 0;
    eventB->SetBranchAddress("pfcand_energy", &pfcand_energy_);
    vector<float>   *pfcand_phi_ = 0;
    eventB->SetBranchAddress("pfcand_phi", &pfcand_phi_);
    vector<float>   *pfcand_eta_ = 0;
    eventB->SetBranchAddress("pfcand_eta", &pfcand_eta_);
    vector<float>   *pfcand_pt_ = 0;
    eventB->SetBranchAddress("pfcand_pt", &pfcand_pt_);
    
    // 
    // Define new variables to write 
    // 
    // need this line because ROOT does not know vector<float>
    // one needs to generate dictionary using rootcint or 
    // just include this line 
    // [ref] http://root.cern.ch/phpBB3/viewtopic.php?t=8467
    gROOT->ProcessLine("#include <vector>"); 
    vector<float>   *fastjets_AK5PF_px = 0;
    TBranch *pxb =  eventB->Branch("fastjets_AK5PF_px", &fastjets_AK5PF_px);
    vector<float>   *fastjets_AK5PF_py = 0;
    TBranch *pyb = eventB->Branch("fastjets_AK5PF_py", &fastjets_AK5PF_py);
    vector<float>   *fastjets_AK5PF_pz = 0;
    TBranch *pzb = eventB->Branch("fastjets_AK5PF_pz", &fastjets_AK5PF_pz);
    vector<float>   *fastjets_AK5PF_energy = 0;
    TBranch *energyb = eventB->Branch("fastjets_AK5PF_energy", &fastjets_AK5PF_energy);
    vector<float>   *fastjets_AK5PF_phi = 0;
    TBranch *phib = eventB->Branch("fastjets_AK5PF_phi", &fastjets_AK5PF_phi);
    vector<float>   *fastjets_AK5PF_eta = 0;
    TBranch *etab = eventB->Branch("fastjets_AK5PF_eta", &fastjets_AK5PF_eta);

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
        fout.open(Form("OneEvent_PFCands_tmp_%i.dat", ib));
        
        // loop over PF candidates
        for(int ipfcand = 0; ipfcand < (int)pfcand_px_->size(); ipfcand++) {
            fout.width(15); fout << pfcand_px_->at(ipfcand) << "\t";
            fout.width(15); fout << pfcand_py_->at(ipfcand) << "\t";
            fout.width(15); fout << pfcand_pz_->at(ipfcand) << "\t";
            fout.width(15); fout << pfcand_energy_->at(ipfcand) << endl;
            
            h2->Fill( pfcand_eta_->at(ipfcand), 
                      pfcand_phi_->at(ipfcand), 
                      pfcand_pt_->at(ipfcand));
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
        gSystem->Exec(Form("../fastjet_example %f < OneEvent_PFCands_tmp_%i.dat > OneEvent_PFJets_tmp_%i.dat", Rparam, ib, ib));

        // ---------------------------------------
        //          Block 3 
        // ---------------------------------------

        // 
        // Get p4 of the reconstructed jets  
        //
        string line;
        double id, px, py , pz, energy, eta, phi, ncand;
        ifstream fin(Form("OneEvent_PFJets_tmp_%i.dat", ib)); 
        if(fin.is_open()) { 
            while(fin.good()){

                // get a line from fin
                getline(fin, line);
                
                // Store each element in the line to the defined variables
                stringstream stream(line);
                stream >> id >> px >> py >> pz >> energy >> eta >> phi >> ncand; 
               
                // store only when pT > 3 GeV 
                // (same as CMS jet reconstruction cut)
                if(TMath::Sqrt(px*px+py*py)>(DEBUG?30:3)) {
                    fastjets_AK5PF_px->push_back(px);
                    fastjets_AK5PF_py->push_back(py);
                    fastjets_AK5PF_pz->push_back(pz);
                    fastjets_AK5PF_energy->push_back(energy);
                    fastjets_AK5PF_eta->push_back(eta);
                    fastjets_AK5PF_phi->push_back(phi); 

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
            TEllipse *cone[fastjets_AK5PF_eta->size()]; 
            for(int ijets=0; ijets<(int)fastjets_AK5PF_eta->size(); ijets++){
                cone[ijets] = new TEllipse(fastjets_AK5PF_eta->at(ijets), fastjets_AK5PF_phi->at(ijets), Rparam, Rparam);
                cone[ijets]->SetFillStyle(3003);
                cone[ijets]->SetFillColor(kYellow);
                cone[ijets]->Draw();
            }

            c->SaveAs(Form("EtaPhiViewPFCand_Run%i_Lumi%i_Event%i_R%.1f.pdf", 
                        run_, lumiblock_, event_, Rparam));
            h2->Reset(); 
            for(int ijets=0; ijets<(int)fastjets_AK5PF_eta->size(); ijets++) delete cone[ijets];

        } 
   
        // Clear vectors for the next event 
        fastjets_AK5PF_px->clear();
        fastjets_AK5PF_py->clear();
        fastjets_AK5PF_pz->clear();
        fastjets_AK5PF_energy->clear();
        fastjets_AK5PF_eta->clear();
        fastjets_AK5PF_phi->clear();
        
        // remove text files
        if(!DEBUG) gSystem->Exec("'rm' OneEvent_PFCands_tmp*.dat OneEvent_PFJets_tmp*.dat");
    
    } // event loop

    // update the tree and close file
    if(!DEBUG) eventB->Write();
    f->Close();
    
    //
    // cleanup
    //
    delete f;
}
