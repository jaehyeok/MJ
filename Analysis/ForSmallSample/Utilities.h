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
    if(Process=="QCD_HT-100To250") Xsec = 10360000.; 
    if(Process=="QCD_HT-250To500") Xsec = 276000.;
    if(Process=="QCD_HT-500To1000") Xsec = 8426.; 
    if(Process=="QCD_HT-1000ToInf") Xsec = 204.; 
    
    // TT
    if(Process=="TTJets_HadronicMGDecays")  Xsec = 106.93; 
    if(Process=="TTJets_FullLeptMGDecays")  Xsec = 24.56; 
    if(Process=="TTJets_SemiLeptMGDecays")  Xsec = 102.5;
    
    // DY
    if(Process=="DYJetsToLL_M-50")  Xsec = 2950.; // cfA 
    // --> Guillelmo's file(~/Migration/scratch/xs.dat) says 1177.6*3 
    // Is 2950.0 with only leading order?

    // W+Jets 
    if(Process=="W2JetsToLNu")  Xsec = 1750.; // cfA
    if(Process=="W3JetsToLNu")  Xsec = 519.; // cfA
    if(Process=="W4JetsToLNu")  Xsec = 214.; // cfA
    
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
