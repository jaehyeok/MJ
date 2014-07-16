#!/bin/bash
# note that INPUTDIR and OUTPUTDIR contain DATASET name as well 
INPUTDIR=$1
OUTPUTDIR=$2
FILE=$3
tmpdir=/data2/jaehyeok/$$

echo "[DEBUG] Set up working directory"
mkdir $tmpdir
cp -r /homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/ $tmpdir
cp -r /homes/jaehyeok/Analysis/MJ/fastjet-install/ $tmpdir
cd $tmpdir/fastjet-3.0.6/example/MJ 
ln -sf ../../../fastjet-install/include/fastjet

echo "[DEBUG] Copy input file"
cp $INPUTDIR/$FILE . 

echo "[DEBUG] Run macro "
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,10\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,15\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,20\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,25\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,30\)  

echo "[DEBUG] Copy output file"
cp $FILE $OUTPUTDIR

echo "[DEBUG] Clean up"
cd -
rm -rf $tmpdir
