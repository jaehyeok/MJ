// 
bool PassNLep(unsigned int Nlep)
{
//    if( (RA4MusPt_->size()+RA4ElsPt_->size())==Nlep 
//        //&& RA4MusVetoPt_->size()==0 
//        //&& RA4ElsVetoPt_->size()==0 
//      ) return true;
//    else return false;
    return (nels_+nmus_)==static_cast<int>(Nlep); 
}

// 
bool PassBaselineSelection(float HT, float MET, int Ncsvm, int Nskinny)
{
   return  (HT>500 && MET>200 && Ncsvm>1 && Nskinny>6); 
}

// 
bool PassSelection(TString Region, 
                   float HT, float MET, int Nb, int Njet, float mT, float MJ)
{
    bool passed=false;
   
    if(Region=="Baseline" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > -1 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR0" 
        && HT   > -1 
        && MET  > 400
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > 600 
    )  passed = true;
    
    if(Region=="SR0p1" 
        && HT   > -1 
        && MET  > 400
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > 500 
    )  passed = true;
    
    if(Region=="SR0p2" 
        && HT   > -1 
        && MET  > 400
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > 400 
    )  passed = true;
    
    if(Region=="SR1" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   < 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR2" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR2p1" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > -1 
        && MJ   > 400 
    )  passed = true;
    
    if(Region=="SR2p2" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > 100 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR2p3" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > 200 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR3" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > 7 
        && mT   > -1 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR4" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > 8 
        && mT   > -1 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR5" 
        && HT   > -1 
        && MET  > -1 && MET < 350 
        && Nb   > -1 
        && Njet > -1 
        && mT   > -1 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR6" 
        && HT   > -1 
        && MET  > 350 && MET < 450 
        && Nb   > -1 
        && Njet > -1 
        && mT   > -1 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR7" 
        && HT   > -1 
        && MET  > 450 
        && Nb   > -1 
        && Njet > -1 
        && mT   > -1 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR8" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > 6 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR9" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > 7 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR10" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > 8 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR11" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > 600 
    )  passed = true;
    
    if(Region=="SR12" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > 500 
    )  passed = true;
    
    if(Region=="SR13" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   < 140 
        && MJ   < 400 
    )  passed = true;
    
    if(Region=="SR14" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > 5 
        && mT   > 140 
        && MJ   > 600 
    )  passed = true;
    
    if(Region=="SR15" 
        && HT   > -1 
        && MET  > 400 
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > 600 
    )  passed = true;
    
    if(Region=="SR16" 
        && HT   > -1 
        && MET  > 400 
        && Nb   > -1 
        && Njet > 3 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR17" 
        && HT   > -1 
        && MET  > 400 
        && Nb   > -1 
        && Njet > 4 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR18" 
        && HT   > -1 
        && MET  > 400 
        && Nb   > 1 
        && Njet > 3 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR19" 
        && HT   > -1 
        && MET  > 400 
        && Nb   > 1 
        && Njet > 4 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR20" 
        && HT   > -1 
        && MET  > 400 
        && Nb   > 1 
        && Njet > 5 
        && mT   > 140 
        && MJ   > 400 
    )  passed = true;
    
    if(Region=="SR20p1" 
        && HT   > -1 
        && MET  > 400 
        && Nb   > 1 
        && Njet > 5 
        && mT   > 140 
        && MJ   > 500 
    )  passed = true;
    
    if(Region=="SR21" 
        && HT   > -1 
        && MET  > 350 
        && Nb   > -1 
        && Njet > 5 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR22" 
        && HT   > -1 
        && MET  > 350 
        && Nb   > -1 
        && Njet > 6 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR23" 
        && HT   > -1 
        && MET  > 350 
        && Nb   > -1 
        && Njet > 7 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR24" 
        && HT   > -1 
        && MET  > 350 
        && Nb   > -1 
        && Njet > 8 
        && mT   > 140 
        && MJ   > -1 
    )  passed = true;
    
    if(Region=="SR25" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > 140 
        && MJ   > 300 
    )  passed = true;
    
    if(Region=="TEST" 
        && HT   > 500 
        && MET  > 400 
        && Nb   > 1 
        && Njet > 6 
        && mT   > 140 
        && MJ   > 600 
    )  passed = true;
   
   return passed;
}
