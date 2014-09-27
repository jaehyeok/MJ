#!/bin/bash

InputDir=/net/cms26/cms26r0/jaehyeok
OutputDir=/net/cms26/cms26r0/jaehyeok/Fatjet
#OutputDir=/net/cms26/cms26r0/jaehyeok/CHS
Dataset=$1
mass=$2

# make resubmit dir 
if [ ! -d "Resubmit" ] 
then 
    echo " [DEBUG] Making directory Resubmit " 
    mkdir -p Resubmit 
fi 


echo "[Resubmit] Making a list of input cfA files"
[ -e Resubmit/cfA_${Dataset}_f${mass}.txt ] && rm Resubmit/cfA_${Dataset}_f${mass}.txt
for cfAPath in `ls $InputDir/$Dataset/*f${mass}_*.root`; do
    cfAFile="$( cut -d '/' -f 7 <<< "$cfAPath" )" 
    cfAFile=`echo $cfAFile | sed 's/_f/#/g' `
    cfAFile=`echo $cfAFile | sed 's/\./#/g' `
    cfAIndex="$( cut -d '#' -f 2 <<< "$cfAFile" )" 
    cfAIndex1=`echo $cfAIndex | cut -d "_" -f 3 `
    cfAIndex2=`echo $cfAIndex | cut -d "_" -f 4 `
    echo $cfAIndex1 $cfAIndex2 >> Resubmit/cfA_${Dataset}_f${mass}.txt
done

echo "[Resubmit] Making a list of existing output files"
[ -e Resubmit/fatjet_${Dataset}_f${mass}.txt ] && rm Resubmit/fatjet_${Dataset}_f${mass}.txt
for FatjetFilePATH in `ls $OutputDir/$Dataset/*f${mass}_*.root`; do
    FatjetFile="$( cut -d '/' -f 8 <<< "$FatjetFilePATH" )"  # should fix this number 7 or 8 depending on the file path 
    FatjetFile=`echo $FatjetFile | sed 's/_f/#/g' `
    FatjetFile=`echo $FatjetFile | sed 's/\./#/g' `
    FatjetIndex="$( cut -d '#' -f 2 <<< "$FatjetFile" )" 
    FatjetIndex1=`echo $FatjetIndex | cut -d "_" -f 3 `
    FatjetIndex2=`echo $FatjetIndex | cut -d "_" -f 4 `
    echo $FatjetIndex1 $FatjetIndex2 >> Resubmit/fatjet_${Dataset}_f${mass}.txt
done 

cat Resubmit/cfA_${Dataset}_f${mass}.txt | while read cfAIndex1 cfAIndex2
do

    echo "0" > Resubmit/matched_tmp.txt
    
    cat Resubmit/fatjet_${Dataset}_f${mass}.txt | while read FatjetIndex1 FatjetIndex2 
    do
        # Here, do comparison 
        if [ $(($cfAIndex1)) -eq $(($FatjetIndex1)) ]
        then 
            if [ $(($cfAIndex2)) -eq $(($FatjetIndex2)) ]
            then 
                echo "1" > Resubmit/matched_tmp.txt
                break
            fi
        fi
    done
    
    Matched=`cat Resubmit/matched_tmp.txt`
     
    # if there is no matching, submit that job 
    if [ $Matched  -eq "0" ] 
    then 
        FilePath=`ls $InputDir/$Dataset/*_f${mass}_${cfAIndex1}_${cfAIndex2}.root`
        File="$( cut -d '/' -f 7 <<< "$FilePath" )"
        echo "[debug] Submit ${cfAIndex1}_${cfAIndex2}"
        echo "JobSubmit.csh ./RunOneJob.sh $InputDir/$Dataset $OutputDir/$Dataset $File" 
        JobSubmit.csh ./RunOneJob.sh $InputDir/$Dataset $OutputDir/$Dataset $File
    fi

done

rm Resubmit/matched_tmp.txt
