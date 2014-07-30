#!/bin/bash 

# set up configuration  
NumOfMaxJob=7

for RemoveMuon in {0,1}; do
	for DoZ in {0,1}; do
		for jetpTthres in {30,20,10}; do
		    
		    echo "------------[ Jet pT > $jetpTthres GeV ]--------------------------------------------------"
		    
		    cat dataset.txt | grep -v "\#" | while read line 
		    do  
		        # 
		        # set up directory and copy files 
		        # 
		        Recognizer=`echo $line | awk '{print $2}'`  
		        if [ ! -d "$Recognizer" ] 
		        then 
		            echo "making directory : $Recognizer" 
		            mkdir $Recognizer  
		        fi
		        cd $Recognizer  
		        cp ../babytree.h .
		        cp ../DoOneProcess.C .
		        ln -sf ../../../MJ_working/aux/
		        
		        # 
		        # make histogram dir
		        # 
		        # Decide diretory name
		        HistDirName="HistFiles"
		        if [ $DoZ -eq 1 ] 
		        then
		            if [ $RemoveMuon -eq 1 ] 
		            then
		                HistDirName="HistFiles_MuonRemoved"
		            else 
		                HistDirName="HistFiles"
		            fi
		        else 
		            if [ $RemoveMuon -eq 1 ] 
		            then
		                HistDirName="HistFiles_Zveto_MuonRemoved"
		            else 
		                HistDirName="HistFiles_Zveto"
		            fi
		        fi
		
		        # make diretory
		        if [ ! -d $HistDirName ] 
		        then 
		            echo "making direcotyr : $HistDirName"
		            mkdir $HistDirName 
		        fi
		        
		        #
		        # run the macro controlling the maximum number of jobs(<= $NumOfMaxJob)
		        #
		        while true 
		        do
		            numjobs=`ps | grep root.exe | wc -l`
		            echo -e "\033[5;34m *** Number of jobs : $(($numjobs)) *** \033[0m"
		            if [ $(($numjobs)) -ge $(($NumOfMaxJob)) ]
		            then
		               echo "[pT>$jetpTthres GeV / DoZ=$DoZ / RemoveMuon=$RemoveMuon] running maximum number of jobs($NumOfMaxJob) now" 
		                sleep 10
		            else 
		                $ROOTSYS/bin/root -q -b DoOneProcess.C++\(\"$Recognizer\",$jetpTthres,$DoZ,$RemoveMuon\) &
		                echo "[pT>$jetpTthres GeV / DoZ=$DoZ / RemoveMuon=$RemoveMuon] $ROOTSYS/bin/root -q -b DoOneProcess.C++\(\"$Recognizer\",$jetpTthres,$DoZ,$RemoveMuon\) &"
		                cd -
		                sleep 10
		                break
		            fi
		            sleep 10
		            echo "[pT>$jetpTthres GeV / DoZ=$DoZ / RemoveMuon=$RemoveMuon] Wait 10s for current jobs to finish"
		        done  # while true
		    done # while line
		done # for jetpT
	done # for DoZ 
done # for RemoveMuon 
