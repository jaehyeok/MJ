BEGINFILE=1
ENDFILE=10
DIR=/A/B/C
echo "[BabyMaker] Copy input cfA files to this directory"
for index in `seq $BEGINFILE $ENDFILE`; do
    ls /net/cms26/cms26r0/jaehyeok/Fatjet/SingleMu_Run2012D-PromptReco-v1_AOD_1/configurableAnalysis_${index}_* 
    echo "cp /net/cms26/cms26r0/jaehyeok/Fatjet/SingleMu_Run2012D-PromptReco-v1_AOD_1/configurableAnalysis_${index}_* ."
done
