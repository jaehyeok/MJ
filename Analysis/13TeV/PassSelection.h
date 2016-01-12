//#include "Regions.h"

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

bool PassMethod2Region(TString Region, 
                       float MJ, float mT, float MET, int Nb, int Njet)
{ 

TString region[24] = {
    "r1_lowmet_lownj_allnb",
    "r2_lowmet_lownj_1b",
    "r2_lowmet_highnj_1b",
    "r2_lowmet_lownj_2b",
    "r2_lowmet_highnj_2b",
    "r2_lowmet_lownj_3b",
    "r2_lowmet_highnj_3b",
    "r3_lowmet_lownj_allnb",
    "r4_lowmet_lownj_1b",
    "r4_lowmet_highnj_1b",
    "r4_lowmet_lownj_2b",
    "r4_lowmet_highnj_2b",
    "r4_lowmet_lownj_3b",
    "r4_lowmet_highnj_3b",
    "r1_highmet_lownj_allnb",
    "r2_highmet_lownj_1b",
    "r2_highmet_highnj_1b",
    "r2_highmet_lownj_2b",
    "r2_highmet_highnj_2b",
    "r3_highmet_lownj_allnb",
    "r4_highmet_lownj_1b",
    "r4_highmet_highnj_1b",
    "r4_highmet_lownj_2b",
    "r4_highmet_highnj_2b"
};

    bool passed=false;

    if(Region==region[0]  && MJ<400 && mT<140 && MET>200 && MET<400 && Njet>=6 && Njet<=99 && Nb>=1 && Nb<=99) passed=true;
    if(Region==region[1]  && MJ>400 && mT<140 && MET>200 && MET<400 && Njet>=6 && Njet<=8  && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[2]  && MJ>400 && mT<140 && MET>200 && MET<400 && Njet>=9 && Njet<=99 && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[3]  && MJ>400 && mT<140 && MET>200 && MET<400 && Njet>=6 && Njet<=8  && Nb>=2 && Nb<=2 ) passed=true;
    if(Region==region[4]  && MJ>400 && mT<140 && MET>200 && MET<400 && Njet>=9 && Njet<=99 && Nb>=2 && Nb<=2 ) passed=true;
    if(Region==region[5]  && MJ>400 && mT<140 && MET>200 && MET<400 && Njet>=6 && Njet<=8  && Nb>=3 && Nb<=99) passed=true;
    if(Region==region[6]  && MJ>400 && mT<140 && MET>200 && MET<400 && Njet>=9 && Njet<=99 && Nb>=3 && Nb<=99) passed=true;
    if(Region==region[7]  && MJ<400 && mT>140 && MET>200 && MET<400 && Njet>=6 && Njet<=99 && Nb>=1 && Nb<=99) passed=true;
    if(Region==region[8]  && MJ>400 && mT>140 && MET>200 && MET<400 && Njet>=6 && Njet<=8  && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[9]  && MJ>400 && mT>140 && MET>200 && MET<400 && Njet>=9 && Njet<=99 && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[10] && MJ>400 && mT>140 && MET>200 && MET<400 && Njet>=6 && Njet<=8  && Nb>=2 && Nb<=2 ) passed=true;
    if(Region==region[11] && MJ>400 && mT>140 && MET>200 && MET<400 && Njet>=9 && Njet<=99 && Nb>=2 && Nb<=2 ) passed=true;
    if(Region==region[12] && MJ>400 && mT>140 && MET>200 && MET<400 && Njet>=6 && Njet<=8  && Nb>=3 && Nb<=99) passed=true;
    if(Region==region[13] && MJ>400 && mT>140 && MET>200 && MET<400 && Njet>=9 && Njet<=99 && Nb>=3 && Nb<=99) passed=true;
    if(Region==region[14] && MJ<400 && mT<140 && MET>400            && Njet>=6 && Njet<=99 && Nb>=1 && Nb<=99) passed=true;
    if(Region==region[15] && MJ>400 && mT<140 && MET>400            && Njet>=6 && Njet<=8  && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[16] && MJ>400 && mT<140 && MET>400            && Njet>=9 && Njet<=99 && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[17] && MJ>400 && mT<140 && MET>400            && Njet>=6 && Njet<=8  && Nb>=2 && Nb<=99) passed=true;
    if(Region==region[18] && MJ>400 && mT<140 && MET>400            && Njet>=9 && Njet<=99 && Nb>=2 && Nb<=99) passed=true;
    if(Region==region[19] && MJ<400 && mT>140 && MET>400            && Njet>=6 && Njet<=99 && Nb>=1 && Nb<=99) passed=true;
    if(Region==region[20] && MJ>400 && mT>140 && MET>400            && Njet>=6 && Njet<=8  && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[21] && MJ>400 && mT>140 && MET>400            && Njet>=9 && Njet<=99 && Nb>=1 && Nb<=1 ) passed=true;
    if(Region==region[22] && MJ>400 && mT>140 && MET>400            && Njet>=6 && Njet<=8  && Nb>=2 && Nb<=99) passed=true;
    if(Region==region[23] && MJ>400 && mT>140 && MET>400            && Njet>=9 && Njet<=99 && Nb>=2 && Nb<=99) passed=true;
    
    return passed;
}

// 
bool PassSelection(TString Region, 
                   float HT, float MET, int Nb, int Njet, float mT, float MJ)
{
    bool passed=false;
    
    if(Region=="TEST"  
        //&& mT   < 140 
        //&& MJ   < 400 
        //&& MET  > 200 
        //&& MET  < 400 
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
