#!/bin/bash
# note that INPUTDIR and OUTPUTDIR contain DATASET name as well 
DIR=$1
RECOGNIZER=$2
BEGINFILE=$3
ENDFILE=$4
ISDATA=$5
LUMI=$6
tmpdir=/data2/jaehyeok/$$

echo "[BabyMaker] Input directory   : $DIR"
echo "[BabyMaker] RECOGNIZER        : $RECOGNIZER"
echo "[BabyMaker] Begin file index  : $BEGINFILE"
echo "[BabyMaker] End file index    : $ENDFILE"
echo "[BabyMaker] ISDATA            : $ISDATA"
echo "[BabyMaker] LUMI              : $LUMI"

echo "[DEBUG] Make temporary directory $tmpdir"
mkdir $tmpdir

echo "[DEBUG] Untar files and set up working directory"
cp /data2/jaehyeok/fastjet-3.0.6_Batch.tar $tmpdir
cp /data2/jaehyeok/fastjet-install.tar $tmpdir
cd $tmpdir
tar xvf fastjet-3.0.6_Batch.tar
tar xvf fastjet-install.tar
cd fastjet-3.0.6_Batch/example/MJ

echo "[BabyMaker] Copy and set up working directory"
cp -r /homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJGit/MJ/Processing/BabyMaking/Batch . 
cd Batch
ln -sf ../../../../fastjet-install/include/fastjet
ln -sf ../../../../fastjet-install/lib/libfastjet.so
pwd

echo "[DEBUG] Run from the file $BEGINFILE for $NFILEPERJOB files"
echo "[DEBUG] Run the macro "

$ROOTSYS/bin/root -b -q DoOneProcess13TeV.C++\(\"$DIR\",\"$RECOGNIZER\",$BEGINFILE,$ENDFILE,$ISDATA,$LUMI\) 
echo "[DEBUG] Copy baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root to /net/cms26/cms26r0/jaehyeok/baby/Fatjet/13TeV"
cp baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root /net/cms26/cms26r0/jaehyeok/baby/Fatjet/13TeV 

#$ROOTSYS/bin/root -b -q DoOneProcess.C++\(\"$DIR\",\"$RECOGNIZER\",$BEGINFILE,$ENDFILE,$ISDATA,$LUMI\) 
#echo "[DEBUG] Copy baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root to /net/cms26/cms26r0/jaehyeok/baby/Fatjet"
#cp baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root /net/cms26/cms26r0/jaehyeok/baby/Fatjet 

echo "[DEBUG] Clean up"
cd $tmpdir 
cd ../
rm -rf $tmpdir
