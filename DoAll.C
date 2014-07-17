//
//  main wrapper
//

void DoAll() {

    gROOT->ProcessLine(".L DoOneProcess.C+");
    gSystem->Load("DoOneProcess_C.so");  

    //
    // QCD
    //
    cout << "... Processing QCD ... " << endl; 
    DoOneProcess("/cms26r0/jaehyeok/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1948_v71/cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1948_v71_f1_*.root","QCD_HT-1000ToInf_UCSB1948",1);

    //
    // TTbar 
    //
//    cout << "... Processing TTbar ... " << endl; 
//    DoOneProcess("/cms26r0/jaehyeok/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71/cfA_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_f1_1_Qb0.root", "TTJets_HadronicMGDecays_UCSB1880", 1);

}
