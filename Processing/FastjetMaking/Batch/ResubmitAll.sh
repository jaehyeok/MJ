#!/bin/bash

for DATASET in `cat datasetCSA.txt | grep -v "#"`; do
    
    echo " [DEBUG] Submitting " $DATASET 
    ./ResubmitOneDataset.sh $DATASET  

done
