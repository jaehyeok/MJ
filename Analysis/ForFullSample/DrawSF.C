//
// h1 cosmetics
//
void h1cosmetic(TH1F* &h1, char* title, int linecolor=kBlack, int linewidth=1, int fillcolor=0, TString var="")
{
    h1->SetLineColor(linecolor);
    h1->SetLineWidth(linewidth);
    h1->SetFillColor(fillcolor);
    h1->SetTitle(title);
    h1->SetXTitle(var);
    h1->SetStats(0);
    h1->SetMinimum(0.1);
}

//
// 
//
void DrawSFOne(TString tag, bool ApplySF)
{
    gStyle->SetTextAlign(22);
    gInterpreter->ExecuteMacro("~/macros/JaeStyle.C");

    TH1F *h1master  = new TH1F("h1master",  "h1master", 4,  -0.5,   3.5);
    TH1F *h1data    = new TH1F("h1data",    "h1dat",    4,  -0.5,   3.5);
    TH1F *h1mc      = new TH1F("h1mc",      "h1mc",     4,  -0.5,   3.5);
    TH1F *h1mcerror = new TH1F("h1mcerror", "h1mcerror",4,  -0.5,   3.5);
    TH1F *h1one     = new TH1F("h1one",     "h1one",    1,  -0.5,   3.5);

    int MJ400bin = 9;
    if(tag.Contains("DY")) MJ400bin=14; 


    TLatex *tex_R[4];

    float chi2global = 0;
    // Loop over Nb
    for(int Nbcut=0; Nbcut<4; Nbcut++) 
    {
        TFile* HistFile;
        HistFile = TFile::Open(Form("%s/Hist_Nb%i_pT30.root", tag.Data(), Nbcut));

        TH1F *h1_DATA, *h1_DATA_MuHad, *h1_DATA_EHad, 
             *h1_T, *h1_TT_sl, *h1_TT_ll, *h1_TT, *h1_WJets, *h1_DY, *h1_MC, *h1_Ratio; 

        if(tag.Contains("HistFilesSL")) 
        {
            h1_DATA_MuHad    = (TH1F*)HistFile->Get("h1_DATA_MuHad_MJ_5fatjet"); 
            h1_DATA_EHad     = (TH1F*)HistFile->Get("h1_DATA_EHad_MJ_5fatjet"); 
            // merge MuHad + EHad
            h1_DATA     = (TH1F*)h1_DATA_MuHad->Clone("h1_DATA");  h1_DATA->Add(h1_DATA_EHad);
        }
        if(tag.Contains("HistFilesDY")) 
            h1_DATA    = (TH1F*)HistFile->Get("h1_DATA_MJ_5fatjet"); 

        h1_T        = (TH1F*)HistFile->Get("h1_T_MJ_5fatjet");
        h1_TT_sl    = (TH1F*)HistFile->Get("h1_TT_sl_MJ_5fatjet");
        h1_TT_ll    = (TH1F*)HistFile->Get("h1_TT_ll_MJ_5fatjet");
        h1_WJets    = (TH1F*)HistFile->Get("h1_WJets_MJ_5fatjet");
        h1_DY       = (TH1F*)HistFile->Get("h1_DY_MJ_5fatjet"); 

        h1_MC = (TH1F*)h1_TT_sl->Clone("h1_MC_MJ_5fatjet");
        h1_MC->Add(h1_TT_ll);
        h1_MC->Add(h1_WJets);
        h1_MC->Add(h1_T);
        h1_MC->Add(h1_DY);

        //
        // Get normalization scale factor using MJ<400
        //
        float NdataMJlt400=0;
        float NdataMJlt400Error=0;
        float NmcMJlt400=0;
        float NmcMJlt400Error=0;
        for(int i=1; i<MJ400bin; i++)
        { 
            NdataMJlt400 = NdataMJlt400 + h1_DATA->GetBinContent(i);
            NdataMJlt400Error = TMath::Sqrt( NdataMJlt400Error*NdataMJlt400Error 
                    +h1_DATA->GetBinError(i)*h1_DATA->GetBinError(i));
            NmcMJlt400 = NmcMJlt400 + h1_MC->GetBinContent(i);
            NmcMJlt400Error = TMath::Sqrt( NmcMJlt400Error*NmcMJlt400Error 
                    +h1_MC->GetBinError(i)*h1_MC->GetBinError(i));
        }

        cout << "------------------------------------------------------------------" << endl;
        cout << "Nb  = " << Nbcut << endl;
        cout << "Data : " << NdataMJlt400 <<  " +/- " << NdataMJlt400Error << endl;
        cout << "MC   : " << NmcMJlt400 <<  " +/- " << NmcMJlt400Error << endl;

        float SFMJlt400 = NdataMJlt400/NmcMJlt400; 
        float SFMJlt400Error = NdataMJlt400/NmcMJlt400
                               *TMath::Sqrt(
                    NdataMJlt400Error*NdataMJlt400Error/NdataMJlt400/NdataMJlt400 
                    +NmcMJlt400Error*NmcMJlt400Error/NmcMJlt400/NmcMJlt400 
                    ); 
        cout << "SF : " << SFMJlt400 <<  " +/- " << SFMJlt400Error << endl;

        // 
        // 
        // 
        float NdataMJgt400=0;
        float NdataMJgt400Error=0;
        float NmcMJgt400=0;
        float NmcMJgt400Error=0;
        for(int i=MJ400bin; i<100; i++)
        { 
            NdataMJgt400 = NdataMJgt400 + h1_DATA->GetBinContent(i);
            NdataMJgt400Error = TMath::Sqrt( NdataMJgt400Error*NdataMJgt400Error 
                    +h1_DATA->GetBinError(i)*h1_DATA->GetBinError(i));
            NmcMJgt400 = NmcMJgt400 + h1_MC->GetBinContent(i);
            NmcMJgt400Error = TMath::Sqrt( NmcMJgt400Error*NmcMJgt400Error 
                    +h1_MC->GetBinError(i)*h1_MC->GetBinError(i));
        }
        cout << "... Before renormalization " << endl;  
        cout << "Data : " << NdataMJgt400 <<  " +/- " << NdataMJgt400Error << endl;
        cout << "MC   : " << NmcMJgt400 <<  " +/- " << NmcMJgt400Error << endl;

        // in case you don't want to apply scalee factors
        if(!ApplySF) 
        {
            SFMJlt400=1;
            SFMJlt400Error=0;
        }

        float NmcMJgt400Renorm       = NmcMJgt400*SFMJlt400;
        float NmcMJgt400ErrorRenorm  = TMath::Sqrt( NmcMJgt400*NmcMJgt400*SFMJlt400Error*SFMJlt400Error
                                                   +SFMJlt400*SFMJlt400*NmcMJgt400Error*NmcMJgt400Error
                                                  );

        cout << "... After renormalization " << endl;  
        cout << "MC   : " << NmcMJgt400Renorm <<  " +/- " << NmcMJgt400ErrorRenorm << endl;

        float chi2=0;
        float R = NdataMJgt400/NmcMJgt400Renorm;
        float Rerror = R*TMath::Sqrt( 1/NdataMJgt400
                                     +NmcMJgt400ErrorRenorm*NmcMJgt400ErrorRenorm/NmcMJgt400Renorm/NmcMJgt400Renorm
                                     );
        cout << "R   : " << R <<  " +/- " << Rerror << endl;


        // Fill histograms  
        h1data->SetBinContent(Nbcut+1,NdataMJgt400/NmcMJgt400Renorm);
        h1data->SetBinError(Nbcut+1,NdataMJgt400Error/NmcMJgt400Renorm);
        h1mcerror->SetBinContent(Nbcut+1,NmcMJgt400Renorm/NmcMJgt400Renorm);
        h1mcerror->SetBinError(Nbcut+1,NmcMJgt400ErrorRenorm/NmcMJgt400Renorm);

        chi2 = (NdataMJgt400-NmcMJgt400Renorm)*(NdataMJgt400-NmcMJgt400Renorm)/
                    (NmcMJgt400ErrorRenorm*NmcMJgt400ErrorRenorm+NdataMJgt400);
        cout << "chi2 : " << chi2 << endl;
        chi2global = chi2global + chi2; 
        //h1master->GetXaxis()->SetBinLabel(Nbcut+1, Form("#splitline{Nb=%i}{#chi^{2}=%.2f}",Nbcut, chi2));
        //h1master->GetXaxis()->SetBinLabel(Nbcut+1, Form("#splitline{Nb=%i}{R=%.2f#pm%.2f}",Nbcut, R, Rerror));
        h1master->GetXaxis()->SetBinLabel(Nbcut+1, Form("n_{b-jets}=%i",Nbcut));
        h1master->GetXaxis()->CenterLabels();
        //h1master->GetXaxis()->SetLabelSize(0.05);
        tex_R[Nbcut] = new TLatex(0.2+0.175/2+Nbcut*0.175,0.06,Form("R=%.2f#pm%.2f", R, Rerror));
        tex_R[Nbcut]->SetNDC();
        tex_R[Nbcut]->SetTextSize(0.045);
        tex_R[Nbcut]->SetLineWidth(2);
        tex_R[Nbcut]->SetTextAlign(22);
    }

    chi2global = chi2global / 4;
    cout << "Global chi2 : " <<  chi2global << endl; 

    TString histtitle = "One-lepton events in CR"; 
    if(tag=="HistFilesSL_AK5PFLep") histtitle  = "SLCR"; 
    if(tag=="HistFilesDY_AK5PFLep") histtitle = "DLCR"; 
    h1cosmetic(h1master,          histtitle.Data(),               kBlack, 2, 0,           "");

    //h1master->SetMaximum(h1data->GetMaximum()*2);
    h1master->SetMaximum(3.);
    h1master->SetMinimum(0.);
    h1master->SetYTitle("Data/MC");
    
    h1one->SetBinContent(1,1); 
    
    h1data->SetLineColor(kBlack);
    h1data->SetMarkerColor(kBlack);
    h1data->SetMarkerSize(1);
    h1data->SetMarkerStyle(20);
    h1data->SetStats(0);

    h1mcerror->SetMarkerSize(0);
    h1mcerror->SetFillColor(kBlack);
    h1mcerror->SetLineColor(kWhite);
    h1mcerror->SetFillStyle(3005);
   
    TLegend* l1 = new TLegend (0.25,0.7,0.5,0.85);
    l1->SetBorderSize(0);
    l1->SetFillColor(0);
    l1->SetFillStyle(0);
    l1->SetTextFont(42);
    l1->SetTextAlign(12);
    l1->SetTextSize(0.06);
    l1->AddEntry(h1data,       "Data/MC",        "lp");
    l1->AddEntry(h1mcerror,    "Uncertainty",    "f");
    
    float textSize = 0.032;
    TLatex *tex_renorm1 = new TLatex(0.6,0.8,"MC renormalized");
    tex_renorm1->SetNDC();
    tex_renorm1->SetTextSize(textSize+0.02);
    tex_renorm1->SetLineWidth(2);
    TLatex *tex_renorm2 = new TLatex(0.6,0.75,"using M_{J}<400 GeV");
    tex_renorm2->SetNDC();
    tex_renorm2->SetTextSize(textSize+0.02);
    tex_renorm2->SetLineWidth(2);
    TString title = "1-lepton signal region";
    if(tag=="HistFilesSL_AK5PFLep") title = "1-lepton control region";
    if(tag=="HistFilesDY_AK5PFLep") title = "2-lepton control region";
    TLatex *tex_title = new TLatex(0.5,0.95,title);
    tex_title->SetNDC();
    tex_title->SetTextSize(textSize+0.02);
    tex_title->SetLineWidth(2);
    tex_title->SetTextAlign(22);
    TLatex *tex_chi2 = new TLatex(0.6,0.67,Form("#chi^{2}/ndf = %.2f",chi2global));
    tex_chi2->SetNDC();
    tex_chi2->SetTextSize(textSize+0.02);
    tex_chi2->SetLineWidth(2);

    TCanvas *c = new TCanvas("c","c", 600, 400); 
    h1master->Draw();
    h1one->Draw("SAME HIST");
    h1mcerror->Draw("SAME E2");
    h1data->Draw("SAME E");
    l1->Draw("same");
    if(ApplySF) tex_renorm1->Draw("same");
    if(ApplySF) tex_renorm2->Draw("same");
    tex_title->Draw("same");
    if(ApplySF)tex_chi2->Draw("same");
    for(int i=0; i<4; i++) tex_R[i]->Draw("same");
    TString outfilename = "SLSR"; 
    if(tag=="HistFilesSL_AK5PFLep") outfilename = "SLCR"; 
    if(tag=="HistFilesDY_AK5PFLep") outfilename = "DLCR"; 
    if(ApplySF) outfilename=outfilename+"_AppliedSF";
    else outfilename=outfilename+"_NOTAppliedSF";
    c->Print(outfilename+".pdf");

    delete c;
}

void DrawSF()
{ 
    DrawSFOne("HistFilesSL_AK5PFLep",   true);
    DrawSFOne("HistFilesSL",            true);
    DrawSFOne("HistFilesDY_AK5PFLep",   true);
    
    DrawSFOne("HistFilesSL_AK5PFLep",   false);
    DrawSFOne("HistFilesSL",            false);
    DrawSFOne("HistFilesDY_AK5PFLep",   false);
}
