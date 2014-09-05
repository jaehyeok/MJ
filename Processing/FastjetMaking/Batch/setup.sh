#!/bin/sh 
#
# This script makes tar files of fastjet package and copy 
# them to each node to save time for copying the package 
# each time a job runs
#

# This script should run in one directory above fastjet package.
# You should change the path according to your directory structure.
cd ../../../../../../ 

# Making tar files
tar cvf fastjet-3.0.6_Batch.tar fastjet-3.0.6_Batch
tar cvf fastjet-install.tar fastjet-install

# Copying tar files to all nodes.
# You should change jaehyeok to your account name
for x in {0..33}; do 
    echo "copying files to cms$x"
    scp fast*tar jaehyeok@cms${x}:/data2/jaehyeok
done

# back to the original directory 
cd -
