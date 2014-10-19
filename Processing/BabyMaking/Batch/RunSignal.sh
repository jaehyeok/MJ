#!/bin/bash

$ROOTSYS/bin/root -b -q DoOneSignal.C++\(\"/net/cms26/cms26r0/jaehyeok/Fatjet/SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-25to500_50GeVX50GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1732reshuf_v68/\",\"T1tttt\",1400,25\) &
$ROOTSYS/bin/root -b -q DoOneSignal.C++\(\"/net/cms26/cms26r0/jaehyeok/Fatjet/SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-25to500_50GeVX50GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1732reshuf_v68/\",\"T1tttt\",1200,25\) &

$ROOTSYS/bin/root -b -q DoOneSignal.C++\(\"/net/cms26/cms26r0/jaehyeok/Fatjet/SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1739reshuf_v68/\",\"T1tttt\",1400,1000\) &
$ROOTSYS/bin/root -b -q DoOneSignal.C++\(\"/net/cms26/cms26r0/jaehyeok/Fatjet/SMS-MadGraph_Pythia6Zstar_8TeV_T1tttt_2J_mGo-1100to1400_mLSP-525to1000_25GeVX25GeV_Binning_Summer12-START52_V9_FSIM-v2_AODSIM_UCSB1739reshuf_v68/\",\"T1tttt\",1200,1000\) &

#echo "... Copying babies to cms26 "
#mv baby_T1ttt* /net/cms26/cms26r0/jaehyeok/baby/Fatjet 
