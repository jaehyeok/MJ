#!/bin/bash 
[ ! -d "HistFiles" ] && echo "making directory : HistFiles" && mkdir HistFiles
[ ! -d "HistFiles_MuonRemoved" ] && echo "making directory : HistFiles_MuonRemoved" && mkdir HistFiles_MuonRemoved
[ ! -d "HistFiles_Zveto" ] && echo "making directory : HistFiles_Zveto" && mkdir HistFiles_Zveto
[ ! -d "HistFiles_Zveto_MuonRemoved_Zveto" ] && echo "making directory : HistFiles_Zveto_MuonRemoved" && mkdir HistFiles_Zveto_MuonRemoved


scp jaehyeok@cms2:/homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJTest/Analysis/*/HistFiles/* HistFiles
scp jaehyeok@cms2:/homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJTest/Analysis/*/HistFiles_MuonRemoved/* HistFiles_MuonRemoved
scp jaehyeok@cms2:/homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJTest/Analysis/*/HistFiles_Zveto/* HistFiles_Zveto
scp jaehyeok@cms2:/homes/jaehyeok/Analysis/MJ/fastjet-3.0.6/example/MJTest/Analysis/*/HistFiles_Zveto_MuonRemoved/* HistFiles_Zveto_MuonRemoved
