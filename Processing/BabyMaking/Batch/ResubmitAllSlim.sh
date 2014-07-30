#!/bin/bash

#for ARG in `less dataset.txt | grep -v "#"`; do
    
cat dataset_slim.txt | grep -v "\#" | while read line  
do 

    echo "-----------------------------------------------------------------------------"
    echo "./ResubmitOneDatasetSlim.sh $line"
    ./ResubmitOneDatasetSlim.sh $line

    DATASET=`echo $line | awk '{print $1}'`
    echo "Done with resubmitting $DATASET"
    echo "-----------------------------------------------------------------------------"

done
