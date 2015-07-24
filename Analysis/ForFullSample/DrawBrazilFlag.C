void DrawBrazilFlag() 
{ 


/*
combine -M Asymptotic -v 1 --rMax 2 datacard_T1tttt_1100_125_MET200_nb2To9_MJ400.dat
1100 125 0.1613 0.2323 0.3703 0.6272 1.0520 0.5137
1100 225 0.1681 0.2438 0.3859 0.6536 1.1045 0.5406
1100 325 0.1885 0.2734 0.4328 0.7330 1.2387 0.6059
1100 425 0.2481 0.3578 0.5672 0.9606 1.6113 0.7918
1100 525 0.3716 0.5352 0.8531 1.4312 2.4216 1.1741
1100 625 0.7595 1.1017 1.7437 2.9532 4.9905 2.4525

combine -M Asymptotic -v 1 --rMax 2 datacard_T1tttt_1100_${x}_MET200To9999_nj6To99_nb2To9_MJ500.dat
1100 125 0.2849 0.3964 0.5953 0.9513 1.5083 0.4364
1100 225 0.2991 0.4149 0.6250 0.9912 1.5809 0.4566
1100 325 0.3327 0.4616 0.6953 1.1027 1.7587 0.5086
1100 425 0.4322 0.5996 0.9031 1.4323 2.2843 0.6608
1100 525 0.6370 0.8838 1.3312 2.1113 3.3370 0.9716
1100 625 1.3206 1.8243 2.7375 4.3743 6.9356 2.0143

combine -M Asymptotic -v 1 --rMax 2 datacard_T1tttt_1150_25_MET200_nb2To9_MJ500.dat
1150 25 0.2845 0.4098 0.6531 1.1165 1.8842 0.9530
1200 25 0.3961 0.5705 0.9094 1.5401 2.6026 1.3150
1250 25 0.5770 0.8319 1.3188 2.2544 3.8044 1.9165
1300 25 0.8030 1.1567 1.8438 3.1225 5.3155 2.6713

combine -M Asymptotic -v 1 --rMax 2 datacard_T1tttt_${x}_25_MET200To9999_nj6To99_nb2To9_MJ500.dat > temp.txt
1150 25 0.3494 0.4996 0.7812 1.2982 2.1484 0.6785
1200 25 0.4822 0.6894 1.0781 1.7915 2.9648 0.9355
1250 25 0.7016 1.0031 1.5688 2.6068 4.3140 1.3637
1300 25 0.9819 1.4037 2.1953 3.6479 6.0370 1.9066
*/

    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");

    int Nmass = 6;  

    float m[]           = {125,225,325,425,525,625};
    float ex[]          = {0,0,0,0,0,0};
    
    // method1 
    float exp2sigdown[] = {0.1613,0.1681,0.1885,0.2481,0.3716,0.7595};
    float exp1sigdown[] = {0.2323,0.2438,0.2734,0.3578,0.5352,1.1017};
    float exp[]         = {0.3703,0.3859,0.4328,0.5672,0.8531,1.7437};
    float exp1sigup[]   = {0.6272,0.6536,0.7330,0.9606,1.4312,2.9532};
    float exp2sigup[]   = {1.0520,1.1045,1.2387,1.6113,2.4216,4.9905};
    float obs[]         = {0.5137,0.5406,0.6059,0.7918,1.1741,2.4525};

    // one bin
//    float exp2sigdown[] = {0.2849,0.2991,0.3327,0.4322,0.6370,1.3206};
//    float exp1sigdown[] = {0.3964,0.4149,0.4616,0.5996,0.8838,1.8243};
//    float exp[]         = {0.5953,0.6250,0.6953,0.9031,1.3312,2.7375};
//    float exp1sigup[]   = {0.9513,0.9912,1.1027,1.4323,2.1113,4.3743};
//    float exp2sigup[]   = {1.5083,1.5809,1.7587,2.2843,3.3370,6.9356};
//    float obs[]         = {0.4364,0.4566,0.5086,0.6608,0.9716,2.0143};

/*
    int Nmass = 4;  
    float m[]           = {1150,1200,1250,1300};
    float ex[]          = {0,0,0,0};
    
    // method1 
    float exp2sigdown[] = {0.2845,0.3961,0.5770,0.8030};
    float exp1sigdown[] = {0.4098,0.5705,0.8319,1.1567};
    float exp[]         = {0.6531,0.9094,1.3188,1.8438};
    float exp1sigup[]   = {1.1165,1.5401,2.2544,3.1225};
    float exp2sigup[]   = {1.8842,2.6026,3.8044,5.3155};
    float obs[]         = {0.9530,1.3150,1.9165,2.6713};
    
    // one bin
//    float exp2sigdown[4] = {0.3494,0.4822,0.7016,0.9819};
//    float exp1sigdown[4] = {0.4996,0.6894,1.0031,1.4037};
//    float exp[4]         = {0.7812,1.0781,1.5688,2.1953};
//    float exp1sigup[4]   = {1.2982,1.7915,2.6068,3.6479};
//    float exp2sigup[4]   = {2.1484,2.9648,4.3140,6.0370};
//    float obs[4]         = {0.6785,0.9355,1.3637,1.9066};
*/
    for(int i=0; i<Nmass; i++)
    { 
        exp1sigup[i]    = TMath::Abs(exp1sigup[i]-exp[i]);
        exp1sigdown[i]  = TMath::Abs(exp1sigdown[i]-exp[i]);
        exp2sigup[i]    = TMath::Abs(exp2sigup[i]-exp[i]);
        exp2sigdown[i]  = TMath::Abs(exp2sigdown[i]-exp[i]);
    }

    TGraph* ObsLim = NULL ;
    TGraph* ExpLim = NULL ;
    TGraphAsymmErrors* ExpBand1sig = NULL ;
    TGraphAsymmErrors* ExpBand2sig = NULL ;
   
    TCanvas *cLimit = new TCanvas("cLimit","cLimit",800,500);
    cLimit->cd(1); 

    // 2sig band
    ExpBand2sig = new TGraphAsymmErrors(Nmass,m,exp,ex,ex,exp2sigdown,exp2sigup);
    ExpBand2sig->SetFillColor(90);
    ExpBand2sig->GetYaxis()->SetRangeUser(0,8.);
    ExpBand2sig->GetXaxis()->SetRangeUser(m[0], m[Nmass-1]);
    ExpBand2sig->GetXaxis()->SetTitleOffset(1.3);
    ExpBand2sig->GetYaxis()->SetTitleOffset(1.7);
    ExpBand2sig->GetXaxis()->SetNdivisions(505);
    ExpBand2sig->GetYaxis()->SetNdivisions(505);
    //ExpBand2sig->GetXaxis()->SetTitle("m_{gluino} [GeV]");
    ExpBand2sig->GetXaxis()->SetTitle("m_{LSP} [GeV]");
    ExpBand2sig->GetYaxis()->SetTitle("95% CL limit on #sigma/#sigma_{theory}");
    ExpBand2sig->Draw("A3");
    
    // 1sig band
    ExpBand1sig = new TGraphAsymmErrors(Nmass,m,exp,ex,ex,exp1sigdown,exp1sigup);
    ExpBand1sig->SetFillColor(211);
    ExpBand1sig->Draw("3");

    // exp limit 
    ExpLim = new TGraph(Nmass,m,exp);
    ExpLim->SetLineWidth(2);
    ExpLim->SetLineStyle(2);
    ExpLim->Draw("l");

    // line at 1 
    TLine *l = new TLine(m[0],1,m[Nmass-1],1); 
    l->SetLineWidth(2);
    l->SetLineColor(kRed+1);
    l->Draw("same");

    // obs limit 
    ObsLim = new TGraph(Nmass,m,obs);
    ObsLim->SetMarkerColor(kBlack);
    ObsLim->SetLineWidth(2);
    ObsLim->SetLineColor(kBlack);
    ObsLim->SetMarkerStyle(kFullDotLarge);
    ObsLim->SetMarkerSize(1.0);
    ObsLim->Draw("lp");

    // Legend 
    leg = new TLegend(0.25, 0.6, 0.6, 0.85, "");
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);

    leg->AddEntry(ExpLim,   " median expected","l");
    leg->AddEntry(ExpBand1sig," expected #pm 1#sigma","f");
    leg->AddEntry(ExpBand2sig," expected #pm 2#sigma","f");
    leg->AddEntry(ObsLim,   " observed","lp");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    leg->Draw("same");

    // CMS Labels
    float textSize = 0.04;

    TLatex *TexEnergyLumi = new TLatex(0.9,0.92,Form("#sqrt{s}=8 TeV, L = %.1f fb^{-1}", 19.5));
    TexEnergyLumi->SetNDC();
    TexEnergyLumi->SetTextSize(textSize);
    TexEnergyLumi->SetTextAlign (31);
    TexEnergyLumi->SetLineWidth(2);

    TLatex *TexCMS = new TLatex(0.2,0.92,"CMS work in progress");
    TexCMS->SetNDC();
    TexCMS->SetTextSize(textSize);
    TexCMS->SetLineWidth(2);

    TLatex *TexExt = new TLatex(0.25,0.53,"m_{gluino} = 1100 GeV");
    //TLatex *TexExt = new TLatex(0.25,0.53,"m_{LSP} = 25 GeV");
    TexExt->SetNDC();
    TexExt->SetTextSize(textSize+0.01);
    TexExt->SetLineWidth(2);

    TexEnergyLumi->Draw("SAME");
    TexCMS->Draw("SAME");
    TexExt->Draw("SAME");

    cLimit->Print("Limit_mLSP_mgluino1100_method1.pdf");
    cLimit->Print("Limit_mLSP_mgluino1100_method1.C");
    //cLimit->Print("Limit_mgluino_mLSP25_method1.pdf");
    //cLimit->Print("Limit_mgluino_mLSP25_method1.C");

}
