#!/bin/bash

cat dataset/datasetPhys14.txt | grep -v "\#" | while read line  
do 

    echo "-----------------------------------------------------------------------------"
    echo "./ResubmitOneDataset.sh $line"
    ./ResubmitOneDataset.sh $line

    DATASET=`echo $line | awk '{print $1}'`
    echo "Done with resubmitting $DATASET"
    #echo "Let's get some sleep for 1 hour"
    #sleep 3600
    echo "-----------------------------------------------------------------------------"

done
