#!/bin/bash

INPUTDIR=/net/cms26/cms26r0/jaehyeok
OUTPUTDIR=/net/cms26/cms26r0/jaehyeok
DATASET=$1

# make an output dir 
if [ ! -d "$OUTPUTDIR/$DATASET/slim" ] 
then 
    echo " [DEBUG] Making directory $OUTPUTDIR/$DATASET/slim " 
    mkdir -p $OUTPUTDIR/$DATASET/slim 
fi 

for file in `ls $INPUTDIR/$DATASET/*f*_*.root`; do
    
    FILE="$( cut -d '/' -f 7 <<< "$file" )" 
    echo " [DEBUG] Doing " $FILE 
    JobSubmit.csh ./RunOneJob.sh $INPUTDIR/$DATASET $OUTPUTDIR/$DATASET/slim $FILE  
    #echo "JobSubmit.csh ./RunOneJob.sh $INPUTDIR/$DATASET $OUTPUTDIR/$DATASET $FILE"  

done
