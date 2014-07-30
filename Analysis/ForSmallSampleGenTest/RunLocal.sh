#!/bin/bash

INPUTDIR=/cms26r0/jaehyeok/genTest/

for FILE in `ls $INPUTDIR/cfA_QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1945_v71_f*.root`; do
    
    echo "[DEBUG] Doing " $FILE 
    echo "[DEBUG] Run macro "
    echo "$ROOTSYS/bin/root -b -q makeJetP4_gen_notext.C++\(\"$FILE\",1.2,pT(R=0.5) cut\)" 
    $ROOTSYS/bin/root -b -q makeJetP4_gen_notext.C++\(\"$FILE\",1.2,10\)  
    $ROOTSYS/bin/root -b -q makeJetP4_gen_notext.C++\(\"$FILE\",1.2,15\)  
    $ROOTSYS/bin/root -b -q makeJetP4_gen_notext.C++\(\"$FILE\",1.2,20\)  
    $ROOTSYS/bin/root -b -q makeJetP4_gen_notext.C++\(\"$FILE\",1.2,25\)  
    $ROOTSYS/bin/root -b -q makeJetP4_gen_notext.C++\(\"$FILE\",1.2,30\)  

done
