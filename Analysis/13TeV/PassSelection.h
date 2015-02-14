
bool PassNLep(int Nlep)
{
    if( (RA4MusPt->size()+RA4ElsPt->size())==Nlep 
        && RA4MusVetoPt->size()==0 
        && RA4ElsVetoPt->size()==0 
      ) return true;
    else return false;
       
}

bool PassBaselineSelection()
{
   return  (HT>750 && MET>250 && NBtagCSVM>1 && Nskinnyjet>5); 
}

bool PassSelection(TString Region, 
                   float HT, float MET, int Nb, int Njet, float mT, float MJ)
{
    bool passed=false;
    if(Region=="SR0" 
        && HT > 750
        && MET > 400
        && Nb > 1 
        && Njet > 5 
        && mT > 150 
        && MJ > 600 
    )  passed = true;
    
    return passed;
}
