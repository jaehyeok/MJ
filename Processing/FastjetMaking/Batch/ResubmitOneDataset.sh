#!/bin/bash

InputDir=/net/cms26/cms26r0/jaehyeok
OutputDir=/net/cms26/cms26r0/jaehyeok/Fatjet
#OutputDir=/net/cms26/cms26r0/jaehyeok/CHS
Dataset=$1

# make resubmit dir 
if [ ! -d "Resubmit" ] 
then 
    echo " [DEBUG] Making directory Resubmit " 
    mkdir -p Resubmit 
fi 


if [ ! -e Resubmit/cfA_${Dataset}.txt ]
then 
    echo "[Resubmit] Making a list of input cfA files"
    for cfAPath in `ls $InputDir/$Dataset/configurableAnalysis_*.root`; do
        cfAFile="$( cut -d '/' -f 7 <<< "$cfAPath" )" 
        cfAFile=`echo $cfAFile | sed 's/_/#/g' `
        cfAFile=`echo $cfAFile | sed 's/\./#/g' `
        cfAIndex="$( cut -d '#' -f 2 <<< "$cfAFile" )" 
        #cfAIndex=`echo $cfAIndex | cut -d "_" -f 1 `
        echo $cfAIndex >> Resubmit/cfA_${Dataset}.txt
    done
fi 

echo "[Resubmit] Making a list of existing output files"
[ -e Resubmit/fatjet_${Dataset}.txt ] && rm Resubmit/fatjet_${Dataset}.txt
for FatjetFilePATH in `ls $OutputDir/$Dataset/configurableAnalysis_*.root`; do
    FatjetFile="$( cut -d '/' -f 8 <<< "$FatjetFilePATH" )"  # should fix this number 7 or 8 depending on the file path 
    FatjetFile=`echo $FatjetFile | sed 's/_/#/g' `
    FatjetFile=`echo $FatjetFile | sed 's/\./#/g' `
    FatjetIndex="$( cut -d '#' -f 2 <<< "$FatjetFile" )" 
    #FatjetIndex=`echo $FatjetIndex | cut -d "_" -f 1 `
    echo $FatjetIndex >> Resubmit/fatjet_${Dataset}.txt
done 

cat Resubmit/cfA_${Dataset}.txt | while read cfAIndex 
do

    echo "0" > Resubmit/matched_tmp.txt
    
    cat Resubmit/fatjet_${Dataset}.txt | while read FatjetIndex 
    do
        # Here, do comparison 
        if [ $(($cfAIndex)) -eq $(($FatjetIndex)) ]
        then 
            echo "1" > Resubmit/matched_tmp.txt
            break
        fi
    done
    
    Matched=`cat Resubmit/matched_tmp.txt`
     
    # if there is no matching, submit that job 
    if [ $Matched  -eq "0" ] 
    then 
        FilePath=`ls $InputDir/$Dataset/configurableAnalysis_${cfAIndex}.root`
        File="$( cut -d '/' -f 7 <<< "$FilePath" )"
        #echo "[debug] Submit ${cfAIndex}"
        echo "JobSubmit.csh ./RunOneJob.sh $InputDir/$Dataset $OutputDir/$Dataset $File" 
        JobSubmit.csh ./RunOneJob.sh $InputDir/$Dataset $OutputDir/$Dataset $File
    fi

done

rm Resubmit/matched_tmp.txt
