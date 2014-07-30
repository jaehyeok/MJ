#!/bin/bash

cat dataset_slim.txt | grep -v "\#" | while read line  
do 
    echo "./RunOneDatasetSlim.sh $line"
    ./RunOneDatasetSlim.sh $line

    DATASET=`echo $line | awk '{print $1}'`
    echo "Done with submitting $DATASET"
    echo "Now wait 10 minutes before submitting jobs for next dataset "
    echo -e "\033[5;34m-----------------------------------------------------------------\033[0m"
    echo -e "\033[5;34m-----------------------------------------------------------------\033[0m"
    #sleep 600 # wait for 10 min before submitting jobs for next dataset 

done
