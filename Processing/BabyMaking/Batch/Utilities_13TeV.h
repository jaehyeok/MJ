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

float GetXsec(TString Process) 
{
   
    // Note that xsec is in pb
    // 13 TeV cross section from https://github.com/manuelfs/ra4_code/blob/master/src/utilities.cpp
    float xsec = -999.;
   
    if(Process.Contains("TTJets"))   xsec = 832; 
    
    if(Process.Contains("T_tW-channel-DR"))          xsec = 35.0;
    if(Process.Contains("TToLeptons_t-channel"))     xsec = 45.0;
    if(Process.Contains("TToLeptons_s-channel"))     xsec = 2.0;
    if(Process.Contains("Tbar_tW-channel-DR"))       xsec = 35.0;
    if(Process.Contains("TBarToLeptons_t-channel"))  xsec = 16.9;
    if(Process.Contains("TBarToLeptons_s-channel"))  xsec = 1.0;

    if(Process.Contains("WJetsToLNu_HT-100to200"))   xsec = 1817.0;
    if(Process.Contains("WJetsToLNu_HT-200to400"))   xsec = 471.6;
    if(Process.Contains("WJetsToLNu_HT-400to600"))   xsec = 55.61;
    if(Process.Contains("WJetsToLNu_HT-600toInf"))   xsec = 18.81;
    
    if(Process.Contains("DYJetsToLL_M-50_HT-100to200"))    xsec = 194.3;
    if(Process.Contains("DYJetsToLL_M-50_HT-200to400"))    xsec = 52.24;
    if(Process.Contains("DYJetsToLL_M-50_HT-400to600"))    xsec = 6.546;
    if(Process.Contains("DYJetsToLL_M-50_HT-600toInf"))    xsec = 2.179;

    if(Process.Contains("QCD_Pt-5to10"))             xsec = 80710000000;
    if(Process.Contains("QCD_Pt-10to15"))            xsec = 7528000000;
    if(Process.Contains("QCD_Pt-15to30"))            xsec = 2237000000;
    if(Process.Contains("QCD_Pt-30to50"))            xsec = 161500000;
    if(Process.Contains("QCD_Pt-50to80"))            xsec = 22110000;
    if(Process.Contains("QCD_Pt-80to120"))           xsec = 3000114;
    if(Process.Contains("QCD_Pt-120to170"))          xsec = 493200;
    if(Process.Contains("QCD_Pt-170to300"))          xsec = 120300;
    if(Process.Contains("QCD_Pt-300to470"))          xsec = 7475;
    if(Process.Contains("QCD_Pt-470to600"))          xsec = 587.1;
    if(Process.Contains("QCD_Pt-600to800"))          xsec = 167;
    if(Process.Contains("QCD_Pt-800to1000"))         xsec = 28.25;
    if(Process.Contains("QCD_Pt-1000to1400"))        xsec = 8.195;
    if(Process.Contains("QCD_Pt-1400to1800"))        xsec = 0.7346;
    if(Process.Contains("QCD_Pt-1800to2400"))        xsec = 0.102;
    if(Process.Contains("QCD_Pt-2400to3200"))        xsec = 0.00644;
    if(Process.Contains("QCD_Pt-3200"))              xsec = 0.000163;

    // Test 
    if(Process=="Test" || Process=="test")  xsec = 1.; 


    if(xsec==-999.) 
    { 
        cout << "[Error] Something is wrong with cross section" << endl;
        return xsec;
    } else { 
        cout << "[MJ Analysis] Cross section =  " << xsec << " pb" << endl;
        return xsec;
    }

}

int GetNjet(vector<float> *vec_px, vector<float> *vec_py, float pTthres) {
    int Njet=0;
    for(int ifatjet=0; ifatjet<(int)vec_px->size(); ifatjet++) {
        if(TMath::Sqrt(vec_px->at(ifatjet)*vec_px->at(ifatjet)          
                      +vec_py->at(ifatjet)*vec_py->at(ifatjet)<pTthres) continue;
        Njet++;
    }
    return Njet;
}
