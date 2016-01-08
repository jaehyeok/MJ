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
   //return  (HT>500 && MET>200 && Ncsvm>1 && Nskinny>6); 
   return  true; 
}

// 
bool PassSelection(TString Region, 
                   float HT, float MET, int Nb, int Njet, float mT, float MJ)
{
    bool passed=false;
    
    if(Region=="TEST" 
        && HT   > 400 
        && MET  > 200 
        && Nb   > 0 
        && Njet > 4 
        //&& mT   > 140 
        //&& MJ   > 400 
    )  passed = true;
   
    if(Region=="Baseline" 
        && HT   > -1 
        && MET  > -1 
        && Nb   > -1 
        && Njet > -1 
        && mT   > -1 
        && MJ   > -1 
    )  passed = true;
    
    // 
    // Method 1
    // 
    if(Region=="Method1" 
        && HT   > 500 
        && MET  > 200 
        && Nb   > 1 
        && Njet > 6 
    )  passed = true;
    
    if(Region=="M1R1" 
        && HT   > 500 
        && MET  > 300 
        && Nb   > 1 
        && Njet > 8 
        && mT   < 140 
        && MJ   < 600 
    )  passed = true;
    
    if(Region=="M1R2" 
        && HT   > 500 
        && MET  > 300 
        && Nb   > 1 
        && Njet > 8 
        && mT   < 140 
        && MJ   > 600 
    )  passed = true;
   
    if(Region=="M1R3" 
        && HT   > 500 
        && MET  > 300 
        && Nb   > 1 
        && Njet > 8 
        && mT   > 140 
        && MJ   < 600 
    )  passed = true;
   
    if(Region=="M1R4" 
        && HT   > 500 
        && MET  > 300 
        && Nb   > 1 
        && Njet > 8 
        && mT   > 140 
        && MJ   > 600 
    )  passed = true;
    
    // 
    // Method 2
    // 
    if(Region=="Method2" 
        && HT   > 500 
        && MET  > 200 
        && Nb   > 0 
        && Njet > 6 
    )  passed = true;
    
    if(Region=="M2R1" 
        && HT   > 500 
        && MET  > 250
        && Nb   > 1 
        && Njet > 6 
        && mT   < 140 
        && MJ   < 400 
    )  passed = true;
    
    if(Region=="M2R2" 
        && HT   > 500 
        && MET  > 250 
        && Nb   > 1 
        && Njet > 8 
        && mT   < 140 
        && MJ   > 400 
    )  passed = true;
    
    if(Region=="M2R3" 
        && HT   > 500 
        && MET  > 250
        && Nb   > 1 
        && Njet > 6 
        && mT   > 140 
        && MJ   < 400 
    )  passed = true;
   
    if(Region=="M2R4" 
        && HT   > 500 
        && MET  > 250
        && Nb   > 1 
        && Njet > 8 
        && mT   > 140 
        && MJ   > 400 
    )  passed = true;

   return passed;
}
