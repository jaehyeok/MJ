//
//  main wrapper
//

void DoAll() {

    // Add a line to remove previously constructed library
    gSystem->Exec("rm DoOneProcess_C.*");
    gROOT->ProcessLine(".L DoOneProcess.C+");
    gSystem->Load("DoOneProcess_C.so");  

    // Set Luminosity 
    //float Lumi = 19500;  // in pb-1
    float Lumi = 1317;  // in pb-1

    // ----------------------------------------
    //              MC
    // ----------------------------------------
    
    //
    // QCD
    //
    DoOneProcess("/cms26r0/jaehyeok/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1948_v71/slim/slim_cfA_QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1948_v71_f11*.root", "QCD_HT-1000ToInf", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1947_v71/slim/slim_cfA_QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1947_v71_f11*.root", "QCD_HT-500To1000", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1946_v71/slim/slim_cfA_QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1946_v71_f11*.root", "QCD_HT-250To500", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1945_v71/slim/slim_cfA_QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1945_v71_f11*.root", "QCD_HT-100To250", false, Lumi);

    //
    // TTbar 
    //
    DoOneProcess("/cms26r0/jaehyeok/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71/cfA_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_f11*.root", "TTJets_HadronicMGDecays", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71/slim/slim_cfA_TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71_f11*.root", "TTJets_FullLeptMGDecays", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71/slim/slim_cfA_TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_f11*.root", "TTJets_SemiLeptMGDecays", false, Lumi);
    
    //
    // W+jets 
    //
    DoOneProcess("/cms26r0/jaehyeok/W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71/cfA_W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71_f11*.root", "W2JetsToLNu", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71/cfA_W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71_f11*.root", "W3JetsToLNu", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71/cfA_W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71_f11*.root", "W4JetsToLNu", false, Lumi);



    //
    // Drell-Yan
    //
    //DoOneProcess("/cms26r0/jaehyeok/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66/cfA_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66_f11*.root", "DYJetsToLL_M-50", false, Lumi);
    //DoOneProcess("/cms26r0/jaehyeok/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66/cfA_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66_f111*.root", "DYJetsToLL_M-50", false, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66/slim/slim_cfA_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66_f11*.root", "DYJetsToLL_M-50", false, Lumi);

    // ----------------------------------------
    //              DATA 
    // ----------------------------------------
    //DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2021_v71/cfA_DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2021_v71_f50_*.root", "DoubleMu_Run2012D-PromptReco", true, Lumi);
    //DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2021_v71/cfA_DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2021_v71_f111*.root", "DoubleMu_Run2012D-PromptReco", true, Lumi);
//    DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012A-13Jul2012-v1_AOD_UCSB2017_v71/cfA_DoubleMu_Run2012A-13Jul2012-v1_AOD_UCSB2017_v71_f*.root", "DoubleMu_Run2012A-13Jul2012", true, Lumi);
//    DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012B-13Jul2012-v4_AOD_UCSB2018_v71/cfA_DoubleMu_Run2012B-13Jul2012-v4_AOD_UCSB2018_v71_f*.root", "DoubleMu_Run2012B-13Jul2012", true, Lumi);
//    DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012C-PromptReco-v2_AOD_UCSB2020_v71/cfA_DoubleMu_Run2012C-PromptReco-v2_AOD_UCSB2020_v71_f*.root", "DoubleMu_Run2012C-PromptReco", true, Lumi);
//    DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012C-24Aug2012-v1_AOD_UCSB2019_v71/cfA_DoubleMu_Run2012C-24Aug2012-v1_AOD_UCSB2019_v71_f*.root", "DoubleMu_Run2012C-24Aug2012", true, Lumi);
    DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2022_v71/slim/slim_cfA_DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2022_v71_f*.root", "DoubleMu_Run2012D-PromptReco2", true, Lumi);
//    DoOneProcess("/cms26r0/jaehyeok/DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2021_v71/cfA_DoubleMu_Run2012D-PromptReco-v1_AOD_UCSB2021_v71_f*.root", "DoubleMu_Run2012D-PromptReco", true, Lumi);
    

    // ----------------------------------------
    //  wMake the final histogram file   
    // ----------------------------------------
    //cout << "[MJ Analysis] Merging result files" << endl;
    //gSystem->Exec("rm ResultRootFiles/Result.root");
    //gSystem->Exec("hadd -f ResultRootFiles/Result.root ResultRootFiles/*_*.root");

}
