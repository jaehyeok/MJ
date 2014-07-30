#!/bin/bash

#########################################
### !! This script is NOT COMPLETE !! ###
#########################################
## developing version is at : 
##    /homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJ_working/SlimBatch 
## There is an issue with Matched being alwasy false 
########################################
Dataset=$1
InputDir=/net/cms26/cms26r0/jaehyeok/$Dataset
OutputDir=/net/cms26/cms26r0/jaehyeok/$Dataset/slim

# make resubmit dir 
if [ ! -d "Resubmit" ] 
then 
    echo " [DEBUG] Making directory Resubmit " 
    mkdir -p Resubmit 
fi 


for NonSlimFilePATH in `ls $InputDir/*f*_*.root`; do
#for nonslimfile in `ls /net/cms26/cms26r0/jaehyeok/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66/cfA_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1572_v66_f1_1_QVN.root`; do
    NonSlimFile="$( cut -d '/' -f 7 <<< "$NonSlimFilePATH" )" 
    NonSlimFile=`echo $NonSlimFile | sed 's/_f/#/g' `
    NonSlimIndex="$( cut -d '#' -f 2 <<< "$NonSlimFile" )" 
    NonSlimIndex=`echo $NonSlimIndex | cut -d "_" -f 1 `
    echo "Non-slim index : $NonSlimIndex"
    
    for SlimFilePATH in `ls $OutputDir/slim*f*_*.root`; do
        SlimFile="$( cut -d '/' -f 8 <<< "$SlimFilePATH" )" 
        SlimFile=`echo $SlimFile | sed 's/_f/#/g' `
        SlimIndex="$( cut -d '#' -f 2 <<< "$SlimFile" )" 
        SlimIndex=`echo $SlimIndex | cut -d "_" -f 1 `
        echo "Slim index : $SlimIndex"
       
        # Here, do comparison 
        Matched=0 
        if [ $(($NonSlimIndex)) -eq $(($SlimIndex)) ]
        then 
            Matched=1
        fi
        # if there is no matching, submit that job 
        if [ $(($Matched)) -eq 0 ] 
        then 
            echo JobSubmit.csh ./RunOneJob.sh $INPUTDIR/$DATASET $OUTPUTDIR/$DATASET/slim $FILE  
            #JobSubmit.csh ./RunOneJob.sh $INPUTDIR/$DATASET $OUTPUTDIR/$DATASET/slim $FILE  
        fi

    done
done

