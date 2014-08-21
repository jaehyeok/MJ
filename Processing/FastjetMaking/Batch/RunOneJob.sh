#!/bin/bash
# note that INPUTDIR and OUTPUTDIR contain DATASET name as well 
INPUTDIR=$1
OUTPUTDIR=$2
FILE=$3
tmpdir=/data2/jaehyeok/$$

echo "[DEBUG] Set up working directory"
mkdir $tmpdir
#if [ ! $(($status)) -eq 0 ] 
#then 
#exit 1
#fi 

echo "[DEBUG] Untar files and set up working directory"
cp /data2/jaehyeok/fastjet-3.0.6_Batch.tar $tmpdir
cp /data2/jaehyeok/fastjet-install.tar $tmpdir 
cd $tmpdir
tar xvf fastjet-3.0.6_Batch.tar
tar xvf fastjet-install.tar
cd $tmpdir/fastjet-3.0.6_Batch/example/MJ 
ln -sf ../../../fastjet-install/include/fastjet

echo "[DEBUG] Copy input file"
cp $INPUTDIR/$FILE . 

echo "[DEBUG] Run macro "

$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,10,\"AK5PFclean\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,15,\"AK5PFclean\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,20,\"AK5PFclean\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,25,\"AK5PFclean\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,30,\"AK5PFclean\"\)  

$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,10,\"AK5PF\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,15,\"AK5PF\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,20,\"AK5PF\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,25,\"AK5PF\"\)  
$ROOTSYS/bin/root -b -q makeJetP4_notext.C++\(\"$FILE\",1.2,30,\"AK5PF\"\)  

echo "[DEBUG] Copy output file"
cp $FILE $OUTPUTDIR

echo "[DEBUG] Clean up"
cd -
rm -rf $tmpdir
