#ifndef INJSON2012_H
#define INJSON2012_H
//root function to filter cfA using JSON File.
//Usage: #include "/afs/cern.ch/user/w/wto/public/scripts/inJSON2012.h"
//Create Vector of Run and Lumi by calling vector< vector<int> > VRunLumi = MakeVRunLumi("Golden") before the event loop.
//You can print of a list of Lumiblock by calling CheckVRunLumi(VRunLumi);
//Check if your run is inJSON by calling bool inJSON(VRunLumi,run,lumiblock) in the event loop..
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
vector< vector<int> > MakeVRunLumi(string input){
	ifstream orgJSON;
	if(input == "Golden" || input == "Prompt"){
	  orgJSON.open("json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt");
	}
	else if(input == "MuonPhys" || input == "Muon"){
		orgJSON.open("json/Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt");
	}
	else if(input == "13Jul2012" || input == "13July2012") {
	  orgJSON.open("json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt");
	}
	else if(input == "06Aug2012") {
	  orgJSON.open("json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt");
	}
	else if(input == "24Aug2012") {
	  orgJSON.open("json/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt");
	}
	else if(input == "11Dec2012") {
	  orgJSON.open("json/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt");
	}
	else if(input == "22Jan2013") {
	  orgJSON.open("json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt");
	}
	else if(input == "DCS" || input == "DCSOnly"){
		orgJSON.open("json/DCSOnly/json_DCSONLY.txt");
	}
	else{
		orgJSON.open(input.c_str());
	}
  vector<int> VRunLumi;
	if(orgJSON.is_open()){
  	char inChar;
  	int inInt;
  	string str;
  	while(!orgJSON.eof()){
  	  char next = orgJSON.peek();
  	  if( next == '1' || next == '2' || next == '3' ||
  	      next == '4' || next == '5' || next == '6' ||
  	      next == '7' || next == '8' || next == '9' || 
  	      next == '0'){     
  	    orgJSON >>inInt;
  	    VRunLumi.push_back(inInt);        
  	  }
  	  else if(next == ' '){
  	    getline(orgJSON,str,' ');
  	  }
  	  else{
  	    orgJSON>>inChar;
  	  }
  	}
  }//check if the file opened.
	else{
		cout<<"Invalid JSON File!\n";
	}
	orgJSON.close();
	if(VRunLumi.size() == 0){
		cout<<"No Lumiblock found in JSON file\n";
	}
  vector< vector<int> > VVRunLumi;
  for(int i = 0; i+2 < int(VRunLumi.size());){
    if(VRunLumi[i] > 130000){
      vector<int> RunLumi;
      RunLumi.push_back(VRunLumi[i]);
      while(VRunLumi[i+1] < 130000 && i+1 < int(VRunLumi.size())){
        RunLumi.push_back(VRunLumi[i+1]);
        i++;
      }
      VVRunLumi.push_back(RunLumi);
      i++;
    }
  }
  return VVRunLumi;
}

bool inJSON(vector< vector<int> > VVRunLumi, int Run, int LS){
	bool answer = false;
	if(Run < 120000){
		answer = true;
	}
	else{
  	for(int i = 0; i < int(VVRunLumi.size());i++){
  	  if(Run == VVRunLumi[i][0]){
  	    for(int j = 1; j+1 < int(VVRunLumi[i].size());j=j+2){
  	      if(LS >= VVRunLumi[i][j] && LS <= VVRunLumi[i][j+1]){
  	        answer = true;
  	      }
  	    }
  	  }
  	}
	}
	return answer;
}

void CheckVRunLumi(vector< vector<int> > VVRunLumi){
	for(int i = 0; i < int(VVRunLumi.size());i++){
		cout<<"Run:"<<VVRunLumi[i][0]<<" LS: ";
		for(int j = 1; j+1 < int(VVRunLumi[i].size());j=j+2){
			cout<<VVRunLumi[i][j]<<"-"<<VVRunLumi[i][j+1]<<" ";
		}
		cout<<endl;
	}
}

void CheckVRunLumi2(vector< vector<int> > VVRunLumi){
	for(int i = 0; i < int(VVRunLumi.size());i++){
		for(int j = 1; j+1 < int(VVRunLumi[i].size());j=j+2){
			if(VVRunLumi[i][j] == VVRunLumi[i][j+1]){
				cout<<VVRunLumi[i][0]<<" "<<VVRunLumi[i][j]<<endl;
			}
			else{
				for(int k=VVRunLumi[i][j];k<VVRunLumi[i][j+1]+1;k++){
					cout<<VVRunLumi[i][0]<<" "<<k<<endl;
				}
			}
		}
		cout<<endl;
	}
}


#endif //INJSON2012_H
