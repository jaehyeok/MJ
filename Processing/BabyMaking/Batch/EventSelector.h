/*
Author: Wing To, HEP, UCSB
Org Date: 02-Sep-2010
Helper function for RA4ElecPtlysis Topbox.
*/
#ifndef EVENTSELECTOR_H
#define EVENTSELECTOR_H
#include "TString.h"
#include "TRegexp.h"

bool GetTriggerNames(vector<TString> & names){
	bool result=false;
	
	for(int i=0;i<int(trigger_decision->size());i++){
		for(int j=0;j<int(names.size());j++){
			if(trigger_name->at(i) == names[j]){
				result = true;
			}
		}
		if(!result){
			names.push_back(trigger_name->at(i));
		}
	}
	return result;
}

bool WriteTriggerNames(vector<TString> names){
	for(int i=0;i<int(names.size());i++){
		cout<<names[i]<<endl;
	}
	return true;
}

bool PassMuonTrig(){
	bool result=false;
	for(int i=0;i<int(trigger_decision->size());i++){
		if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){  
			// Trigger on or of all
			vector<TRegexp> MyTrigs;
			//MyTrigs.push_back("HLT_Mu40_PFHT350_v"); //control sample
			//MyTrigs.push_back("HLT_Mu60_PFHT350_v"); //control sample
                        MyTrigs.push_back("HLT_PFHT350_Mu15_PFMET45_v");//
                        MyTrigs.push_back("HLT_PFHT350_Mu15_PFMET50_v");//
                        MyTrigs.push_back("HLT_PFHT400_Mu5_PFMET45_v");//
                        MyTrigs.push_back("HLT_PFHT400_Mu5_PFMET50_v");//
                        MyTrigs.push_back("HLT_PFNoPUHT350_Mu15_PFMET45_v");//
                        MyTrigs.push_back("HLT_PFNoPUHT350_Mu15_PFMET50_v");//
                        MyTrigs.push_back("HLT_PFNoPUHT400_Mu5_PFMET45_v");//
                        MyTrigs.push_back("HLT_PFNoPUHT400_Mu5_PFMET50_v");//
			for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
				if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
					result = true;
					//cout<<"Triggered on "<<trigger_name->at(i)<<endl;
				}
			}			
		}
	}
	return result;
}

bool PassSingleMuonTrig(){
        bool result=false;
        for(int i=0;i<int(trigger_decision->size());i++){
                if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
                        // Trigger on or of all 
                        vector<TRegexp> MyTrigs;
                        MyTrigs.push_back("HLT_IsoMu24_v");//
                        MyTrigs.push_back("HLT_IsoMu24_eta2p1_v");//
                        for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
                                if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
                                        result = true;
                                        //cout<<"Triggered on "<<trigger_name->at(i)<<endl;
                                }       
                        }        
                }       
        }       
        return result;
}  

bool PassHTMuonTrig(){
        bool result=false;
        for(int i=0;i<int(trigger_decision->size());i++){
                if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
                        // Trigger on or of all
                        vector<TRegexp> MyTrigs;
                        MyTrigs.push_back("HLT_Mu40_PFHT350_v"); //control sample
                        MyTrigs.push_back("HLT_Mu60_PFHT350_v"); //control sample
                        MyTrigs.push_back("HLT_Mu40_PFNoPUHT350_v"); //control sample
                        MyTrigs.push_back("HLT_Mu60_PFNoPUHT350_v"); //control sample
                        for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
                                if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
                                        result = true;
                                        //cout<<"Triggered on "<<trigger_name->at(i)<<endl;
                                }
                        }
                }
        }
        return result;
}


bool PassDiMuonTrig(){
        bool result=false;
        for(int i=0;i<int(trigger_decision->size());i++){
                if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
                        // Trigger on or of all
                        vector<TRegexp> MyTrigs;
                        MyTrigs.push_back("HLT_Mu17_Mu8_v");
                        for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
                                if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
                                        result = true;
                                        //cout<<"Triggered on "<<trigger_name->at(i)<<endl;
                                }
                        }
                }
        }
        return result;
}

bool PassElecTrig(){ 
	bool result=false;
	if(run < 130000) {
		return true;
	}
	else{
		for(int i=0;i<int(trigger_decision->size());i++){
			if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
				vector<TRegexp> MyTrigs;
				// 2012 Triggers
                                //MyTrigs.push_back("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v"); //control sample
                                //MyTrigs.push_back("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v"); //control sample
                                MyTrigs.push_back("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v");//
                                MyTrigs.push_back("HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v");//
                                MyTrigs.push_back("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v");//
                                MyTrigs.push_back("HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v");//
                                MyTrigs.push_back("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v");//
                                MyTrigs.push_back("HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v");//
                                MyTrigs.push_back("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v");//
                                MyTrigs.push_back("HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v");//
				for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
					if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
						result = true;
						//cout<<"Triggered on "<<trigger_name->at(i)<<endl;
					}
				}
				
			}
		}
		return result;
	}
}

bool PassSingleElecTrig(){
        bool result=false;
        for(int i=0;i<int(trigger_decision->size());i++){
                if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
                        // Trigger on or of all 
                        vector<TRegexp> MyTrigs;
                        MyTrigs.push_back("HLT_Ele27_WP80_v");//
                        for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
                                if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
                                        result = true;
                                        //cout<<"Triggered on "<<trigger_name->at(i)<<endl;
                                }
                        }
                }
        }
        return result;
}

bool PassHTElecTrig(){
        bool result=false;
        if(run < 130000) {
                return true;
        }       
        else{
                for(int i=0;i<int(trigger_decision->size());i++){
                        if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
                                vector<TRegexp> MyTrigs;
                                // 2012 Triggers
                                MyTrigs.push_back("HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v"); //control sample
                                MyTrigs.push_back("HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v"); //control sample
                                MyTrigs.push_back("HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v"); //control sample
                                MyTrigs.push_back("HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v"); //control sample
                                for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
                                        if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
                                                result = true;
                                                //cout<<"Triggered on "<<trigger_name->at(i)<<endl;
                                        }
                                }

                        }
                }
                return result;
        }
}


bool PassDiElecTrig(){
        bool result=false;
        if(run < 130000) {
                return true;
        }
        else{
                for(int i=0;i<int(trigger_decision->size());i++){
                        if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
                                vector<TRegexp> MyTrigs;
                                // 2012 Triggers
                                MyTrigs.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
                                for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
                                        if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
                                                result = true;
                                                //cout<<"Triggered on "<<trigger_name->at(i)<<endl;
                                        }
                                }

                        }
                }
                return result;
        }
}


bool PassHTTrig(){
	bool result=false;
	if(run < 100000) {
		return true;
	}
	else{for(int i=0;i<int(trigger_decision->size());i++){
			if(trigger_prescalevalue->at(i) == 1 && trigger_decision->at(i)==1){
				vector<TRegexp> MyTrigs;
				// Wing's Triggers
				MyTrigs.push_back("HLT_HT200_v");
				MyTrigs.push_back("HLT_HT250_v");
				MyTrigs.push_back("HLT_HT300_v");
				MyTrigs.push_back("HLT_HT350_v");
				MyTrigs.push_back("HLT_HT400_v");
				MyTrigs.push_back("HLT_HT450_v");
				MyTrigs.push_back("HLT_HT500_v");
				MyTrigs.push_back("HLT_HT550_v");
				MyTrigs.push_back("HLT_HT600_v");
				for(int j=0;(j<int(MyTrigs.size())) && !result;j++){
					if(TString(trigger_name->at(i)).Contains(MyTrigs[j])){
						result = true;
						//cout<<"Triggered on "<<trigger_name->at(i)<<endl;
					}
				}
				
			}
		}
		return result;
	}
}
/*
float GetPVNDof(int & bestPV, int & nGoodPV){
	float BestNDof=0;
	for(unsigned int i=0;i < pv_x->size();i++){ // loop over all Primary Vertices
		if(fabs(pv_isFake->at(i))< 0.001) Hist2D["PVNDoF_Vs_PVNum"]->Fill(i,pv_ndof->at(i));
		if(pv_ndof->at(i) > 4
		&& sqrt(pv_x->at(i)*pv_x->at(i)+pv_y->at(i)*pv_y->at(i)) <= 2
		&& fabs(pv_z->at(i)) < 24){
			nGoodPV++;
			if(pv_ndof->at(i) > BestNDof) { //pick the PV with high ndof.
			BestNDof=pv_ndof->at(i);
			bestPV=i;
		} //pick the PV with highest ndof for later use.
		}
	} // loop over all Primary Vertices
	return BestNDof;
}
*/

#endif //EVENTSELECTOR_H
