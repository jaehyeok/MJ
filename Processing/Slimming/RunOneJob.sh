#!/bin/bash
# note that INPUTDIR and OUTPUTDIR contain DATASET name as well 
INPUTDIR=$1
OUTPUTDIR=$2
FILE=$3
OUTFILE=slim_$3
tmpdir=/data2/jaehyeok/$$

echo "[Slimmer] Input directory : $INPUTDIR"
echo "[Slimmer] Input file : $FILE"
echo "[Slimmer] Output directory : $OUTPUTDIR"
echo "[Slimmer] Output file : $OUTFILE"

echo "[Slimmer] Set up working directory"
mkdir $tmpdir
cd $tmpdir 
pwd
cp /homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJ_working/SlimBatch/Slimmer.C .
cp /homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJ_working/SlimBatch/Slimmer.h . 

echo "[DEBUG] Copy input file"
cp $INPUTDIR/$FILE . 

echo "[DEBUG] Run macro "
$ROOTSYS/bin/root -b -q Slimmer.C++\(\"$FILE\",\"$OUTFILE\"\)  

echo "[DEBUG] Copy $OUTFILE to $OUTPUTDIR"
cp $OUTFILE $OUTPUTDIR

echo "[DEBUG] Clean up"
cd -
rm -rf $tmpdir
