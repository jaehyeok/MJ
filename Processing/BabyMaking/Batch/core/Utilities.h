//
// Dump of utility functions 
//

//TH1F initialization
TH1F* InitTH1F(char* Name, char* Title, int Nbins, double XMin, double XMax){

    TH1F *h1 = new TH1F(Name, Title, Nbins, XMin, XMax);
    h1->Sumw2();
    return h1;
}

// PU reweigting 
double nPUScaleFactor2012(TH1F* h1PU, float npu){
    // Get PU scale factor histogram 
    double mynpu = TMath::Min(npu,(float)49.499);
    Int_t npuxbin = h1PU->GetXaxis()->FindBin(mynpu);
    return h1PU->GetBinContent(npuxbin);
}

float GetXsec(TString Process) {
    
    // Note that Xsec is in pb
    // Ref : http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2012_350_v6.pdf
    // On my laptop : /Users/jaehyeok/Downloads/AN2012_350_v6.pdf
    // Some are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/CfaNtuple2012 [marked as cfA]
    float Xsec = -999.;
    
    // QCD 
    if(Process=="QCD_HT-100To250_TuneZ2star_8TeV")  Xsec = 10360000.; 
    if(Process=="QCD_HT-250To500_TuneZ2star_8TeV")  Xsec = 276000.;
    if(Process=="QCD_HT-500To1000_TuneZ2star_8TeV") Xsec = 8426.; 
    if(Process=="QCD_HT-1000ToInf_TuneZ2star_8TeV") Xsec = 204.; 
    
    // TT
    if(Process=="TTbar_TuneZ2star_13TeV")           Xsec = 832; 
    if(Process=="TT_CT10_TuneZ2star_8TeV")          Xsec = 234; 
    if(Process=="TTJets_HadronicMGDecays_8TeV")     Xsec = 106.93; 
    if(Process=="TTJets_FullLeptMGDecays_8TeV")     Xsec = 24.56; 
    if(Process=="TTJets_SemiLeptMGDecays_8TeV")     Xsec = 102.5;
   
    // tW 
    if(Process=="T_tW-channel-DR_TuneZ2star_8TeV")      Xsec = 11.1;
    if(Process=="Tbar_tW-channel-DR_TuneZ2star_8TeV")   Xsec = 11.1;
    
    // single top 
    if(Process=="T_t-channel_TuneZ2star_8TeV")          Xsec = 56.4;
    if(Process=="Tbar_t-channel_TuneZ2star_8TeV")       Xsec = 30.7;
    if(Process=="T_s-channel_TuneZ2star_8TeV")          Xsec = 3.79;
    if(Process=="Tbar_s-channel_TuneZ2star_8TeV")       Xsec = 1.76;

    // DY
    if(Process=="DYJetsToLL_M-50_8TeV")                     Xsec = 1177.6*3;  
    if(Process=="DYJetsToLL_HT-200To400_TuneZ2Star_8TeV")   Xsec = 23.43;  
    if(Process=="DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV")   Xsec = 3.36;  

    // W+Jets 
    if(Process=="W2JetsToLNu_8TeV")  Xsec = 1750.;
    if(Process=="W3JetsToLNu_8TeV")  Xsec = 519.; 
    if(Process=="W4JetsToLNu_8TeV")  Xsec = 214.;
    //if(Process=="WJetsToLNu_HT-150To200_8TeV")  Xsec = 1.; // FIXME 
    //if(Process=="WJetsToLNu_HT-200To250_8TeV")  Xsec = 1.; //FIXME 
    if(Process=="WJetsToLNu_HT-250To300_8TeV")  Xsec = 57.26;
    if(Process=="WJetsToLNu_HT-300To400_8TeV")  Xsec = 45.68; 
    if(Process=="WJetsToLNu_HT-400ToInf_8TeV")  Xsec = 30.08;
    
    // ttV 
    if(Process=="TTZJets_8TeV")  Xsec = 0.2057; // cfA
    if(Process=="TTWJets_8TeV")  Xsec = 0.232; // cfA
    
    // Test 
    if(Process=="Test" || Process=="test")  Xsec = 1.; 


    if(Xsec==-999.) { 
        cout << "[Error] Something is wrong with cross section" << endl;
        return Xsec;
    } else { 
        cout << "[MJ Analysis] Cross section =  " << Xsec << " pb" << endl;
        return Xsec;
    }

}

int GetNjet(vector<float> *vec_px, vector<float> *vec_py, float pTthres) {
    int Njet=0;
    for(int ifatjet=0; ifatjet<(int)vec_px->size(); ifatjet++) {
        if(TMath::Sqrt(vec_px->at(ifatjet)*vec_px->at(ifatjet)          
                      +vec_py->at(ifatjet)*vec_py->at(ifatjet))<pTthres) continue;
        Njet++;
    }
    return Njet;
}
