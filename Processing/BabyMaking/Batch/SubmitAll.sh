#!/bin/bash


cat dataset/datasetPhys14.txt | grep -v "\#" | while read line  
do 
    
    echo "./RunOneDataset.sh $line"
    ./RunOneDataset.sh $line

    DATASET=`echo $line | awk '{print $1}'`
    echo "Done with submitting $DATASET"
    echo "Now wait 90 minutes before submitting jobs for next dataset "
    sleep 5400 # wait for 90 min before submitting jobs for next dataset 
    echo -e "\033[5;34m-----------------------------------------------------------------\033[0m"
    echo -e "\033[5;34m-----------------------------------------------------------------\033[0m"

done
