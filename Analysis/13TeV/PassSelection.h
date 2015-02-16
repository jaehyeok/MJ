// 
bool PassNLep(unsigned int Nlep)
{
    if( (RA4MusPt_->size()+RA4ElsPt_->size())==Nlep 
        && RA4MusVetoPt_->size()==0 
        && RA4ElsVetoPt_->size()==0 
      ) return true;
    else return false;
       
}

// 
bool PassBaselineSelection()
{
   return  (HT_>750 && MET_>250 && NBtagCSVM_>1 && Nskinnyjet_>5); 
}

// 
bool PassSelection(TString Region, 
                   float HT, float MET, int Nb, int Njet, float mT, float MJ)
{
    bool passed=false;
    
    if(Region=="SR0" 
        && HT > -1 
        && MET > 400
        && Nb > -1 
        && Njet > -1 
        && mT > 150 
        && MJ > 600 
    )  passed = true;
    
    if(Region=="SR1" 
        && HT > 750
        && MET > -1 
        && Nb > -1 
        && Njet > -1 
        && mT < 150 
        && MJ > -1 
    )  passed = true;
    
    return passed;

}
