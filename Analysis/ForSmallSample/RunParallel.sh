#
# This script is meant to run on UCSB clusters
# 
# Directories for each process must be made first 
# then DoAll.C needs to be copied to each direcory 
# with only correponding datasets are turned on
#

echo " -------------------------"
echo " begin time"
echo " -------------------------"
date

cd TT
pwd
cp ../DoOneProcess.C .
cp ../Utilities.h .
cp ../inJSON2012.h .
cp ../ObjectSelector_Sync.h .
cp ../slimmedBranch.h .
ln -sf ../aux/
ln -sf ../json/
$ROOTSYS/bin/root -b -q DoAll.C &
cd - 

cd WJets 
pwd
cp ../DoOneProcess.C .
cp ../Utilities.h .
cp ../inJSON2012.h .
cp ../ObjectSelector_Sync.h .
cp ../slimmedBranch.h .
ln -sf ../aux/
ln -sf ../json/
$ROOTSYS/bin/root -b -q DoAll.C &
cd - 

cd QCD 
pwd
cp ../DoOneProcess.C .
cp ../Utilities.h .
cp ../inJSON2012.h .
cp ../ObjectSelector_Sync.h .
cp ../slimmedBranch.h .
ln -sf ../aux/
ln -sf ../json/
$ROOTSYS/bin/root -b -q DoAll.C &
cd - 

cd DY 
pwd
cp ../DoOneProcess.C .
cp ../Utilities.h .
cp ../inJSON2012.h .
cp ../ObjectSelector_Sync.h .
cp ../slimmedBranch.h .
ln -sf ../aux/
ln -sf ../json/
$ROOTSYS/bin/root -b -q DoAll.C &
cd - 

cd DATA 
pwd
cp ../DoOneProcess.C .
cp ../Utilities.h .
cp ../inJSON2012.h .
cp ../ObjectSelector_Sync.h .
cp ../slimmedBranch.h .
ln -sf ../aux/
ln -sf ../json/
$ROOTSYS/bin/root -b -q DoAll.C &
cd - 

## copy 
#cp */baby/*root baby
