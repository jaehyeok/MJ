#!/bin/bash

INPUTDIR=/net/cms2/cms2r0/cfA              # 13 TeV or official 8 TeV
#INPUTDIR=/net/cms26/cms26r0/jaehyeok/Fatjet # unofficial 8 TeV
DATASET=$1
RECOGNIZER=$2
ISDATA=$3
LUMI=$4
NFILEPERJOB=$5

# make an output dir 
#if [ ! -d "$OUTPUTDIR/$DATASET/slim" ] 
#then 
#    echo " [DEBUG] Making directory $OUTPUTDIR/$DATASET/slim " 
#    mkdir -p $OUTPUTDIR/$DATASET/slim 
#fi 

#### cleaning  
rm filenumber_resubmit.txt filenumbersorted_resubmit.txt

#####
for NONSLIMFILEPATH in `ls $INPUTDIR/$DATASET/*root`; do
    NONSLIMFILE="$( cut -d '/' -f 7 <<< "$NONSLIMFILEPATH" )" # 13 TeV or official 8 TeV 
#    NONSLIMFILE="$( cut -d '/' -f 8 <<< "$NONSLIMFILEPATH" )"  # unofficial 8 TeV 
    NONSLIMFILE=`echo $NONSLIMFILE | sed 's/_f/#/g' `
    NONSLIMFILE=`echo $NONSLIMFILE | sed 's/.root/#root/g' `
    NONSLIMFILE=`echo $NONSLIMFILE | sed 's/configurableAnalysis_/configurableAnalysis#/g' ` # non-published samples
    NONSLIMINDEX="$( cut -d '#' -f 2 <<< "$NONSLIMFILE" )" 
    NONSLIMINDEX=`echo $NONSLIMINDEX | cut -d "_" -f 1 `
    echo $NONSLIMINDEX >> filenumber_resubmit.txt
done 

#####
sort -g filenumber_resubmit.txt > filenumbersorted_resubmit.txt
NUMLINE=`wc -l filenumbersorted_resubmit.txt | awk '{print $1}'`
LASTFILENUM=`tail -n +$NUMLINE filenumbersorted_resubmit.txt`  
echo "The Last file number : $LASTFILENUM"

####
COUNTER=1
NUMOUTPUTFILES=0
while [ $COUNTER -lt $(($LASTFILENUM)) ]; do
    #echo The counter is $COUNTER
    if [ $(($LASTFILENUM-$COUNTER)) -lt $(($NFILEPERJOB)) ] 
    then 
        BEGINFILE=$COUNTER
        ENDFILE=$(($LASTFILENUM))
        if [ ! -e /net/cms26/cms26r0/jaehyeok/baby/Fatjet/13TeV/baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root ] 
#        if [ ! -e /net/cms26/cms26r0/jaehyeok/baby/Fatjet/baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root ] 
        then
            echo -e "\033[5;34m JobSubmit.csh ./RunOneJob.sh  $INPUTDIR/$DATASET/ $RECOGNIZER $COUNTER $(($LASTFILENUM)) $ISDATA $LUMI \033[0m"
#            JobSubmit.csh ./RunOneJob.sh  $INPUTDIR/$DATASET/ $RECOGNIZER $COUNTER $(($LASTFILENUM)) $ISDATA $LUMI
        fi 
    else  
        BEGINFILE=$COUNTER
        ENDFILE=$(($COUNTER+$NFILEPERJOB-1))
        if [ ! -e /net/cms26/cms26r0/jaehyeok/baby/Fatjet/13TeV/baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root ] 
#        if [ ! -e /net/cms26/cms26r0/jaehyeok/baby/Fatjet/baby_${RECOGNIZER}_f${BEGINFILE}To${ENDFILE}.root ] 
        then
            echo -e "\033[5;34m JobSubmit.csh ./RunOneJob.sh $INPUTDIR/$DATASET/ $RECOGNIZER $COUNTER $(($COUNTER+$NFILEPERJOB-1)) $ISDATA $LUMI \033[0m" 
#            JobSubmit.csh ./RunOneJob.sh $INPUTDIR/$DATASET/ $RECOGNIZER $COUNTER $(($COUNTER+$NFILEPERJOB-1)) $ISDATA $LUMI 
        fi
    fi
    let NUMOUTPUTFILES=NUMOUTPUTFILES+1
    let COUNTER=COUNTER+$NFILEPERJOB
done

NUMEXISTOUTPUTFILES=`ls /net/cms26/cms26r0/jaehyeok/baby/Fatjet/13TeV/baby_${RECOGNIZER}_f*To*.root | wc -l`
#NUMEXISTOUTPUTFILES=`ls /net/cms26/cms26r0/jaehyeok/baby/Fatjet/baby_${RECOGNIZER}_f*To*.root | wc -l`
echo "Number of target output files : $NUMOUTPUTFILES"
echo "Number of existing output files : $NUMEXISTOUTPUTFILES"

# 
# merge babies
# 
if [ $(($NUMEXISTOUTPUTFILES)) -eq $(($NUMOUTPUTFILES)) ]
then 
    echo "Babyies for $DATASET are complete" 
    echo -e "\033[5;34m No need for resubmission.\033[0m"
#    echo "hadd /net/cms26/cms26r0/jaehyeok/baby/Fatjet/baby_${RECOGNIZER}.root /net/cms26/cms26r0/jaehyeok/baby/Fatjet/baby_${RECOGNIZER}_f*To*.root"
fi

