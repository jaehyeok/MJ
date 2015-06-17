MJ Analysis code
=== 

Analysis code for MJ analysis at 13 TeV. Now it runs over Jae's baby ntuples. The contents can be found here : https://github.com/jaehyeok/MJ/blob/master/Processing/BabyMaking/Batch/DoOneProcess13TeV.C. The code will be updated so that it can run over Manuel Adams ntuples.

#### Structure of code  

The codes are composed of several ROOT macros and a couple header files. The ROOT macros are capsulized according to tasks. The reason why ROOT macros are chosen over header files, is to use them independently in case you want to run a part of the code. For example, if you want to re-make only MJ plot, you can just run `Make1DPlot.C` without running the whole machinery again and this can save much time depending on the size of baby ntuples. Following is the list of scripts and what they do :

    babytree.h : sets branch addresses 
    PassSelection.h : defines selection 
    DoAnalysis.C : is the main wrapper that runs following macros     
    MakeHists.C : loops over events and make histograms to be used for making plots and tables 
    Make1DPlots.C : makes 1D plots from the histrogram file genereated by MakeHists.C 
    Make2DPlots.C : makes 2D plots from the histrogram file genereated by MakeHists.C 
    MakeTables.C : makes table of yields from the histrogram file genereated by MakeHists.C. It can print out in two formats, plain text or Latex. 
    MakeCards.C : makes cards for combine tool. It generates cards for electron, muon and inclusive channel. [* Right now, it is in a very preliminary form, but once we have real systematics, then more realistic cards will be made.]    


#### How to run the code  

Before running the code, you need to set the directory of babies appropriately. Modify this line in `DoAnalysis.C` 

```
TString BabyDir = "/Users/jaehyeok/Research/Tools/fastjet-3.0.6/example/babies/13TeV/HT750MET250/";
```

In addition, make sure you have following directories 

```
HistFiles
Figures
Tables
Cards
```

If not, please make them. 

Now, you are ready to run the code. Do 

```
root -q -b DoAnalysis.C
```

Then, you will have histogram files in `HistFiles`, plots in pdf format in `Figures` and tables of yields in Latex format in `Tables`.  

#### Possible improvements 

One obvious thing is to write a separate code to do background estimation 

```
Calculation of k1 and k2, and their uncertainties (which syst uncert to include needs to be determined)
Application of these 
```

