#!/bin/bash

for DATASET in `cat dataset.txt | grep -v "#"`; do
    
    echo " [DEBUG] Submitting " $DATASET 
    ./ResubmitOneDataset.sh $DATASET  

done
