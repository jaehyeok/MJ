# This script is to tar fastjec package and copy it to 
# the ucsb nodes in order to save time to copy files 
# from home direcory to nodes for every job
#
# This script needs to be copied to the top of fastjet directory
#
tar cvf fastjet-3.0.6_Batch.tar fastjet-3.0.6_Batch
tar cvf fastjet-install.tar fastjet-install

for x in {0..33}; do 
    echo "copying files to cms$x"
    scp fast*tar jaehyeok@cms${x}:/data2/jaehyeok
done
