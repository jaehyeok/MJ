MJ Analysis

-----------------------------
 Some comments
-----------------------------
* All of these work provided that fastjet package is set up 

-----------------------------
 Structure of directory 
-----------------------------

/Analysis
    |------ ForFullSample
        - Full 8 TeV sample analysis
        - Make histograms from babies per data sample and merge them 
        - Making histograms from babies done in parallel on UCSB HEP cluster
          (preferably on cms24/26 where babies are)
    |------ ForSmallSample
        - Analysis for a chunk of 8 TeV data (1.317 /fb)
        - Everything done locally, parallelized if necessary
    |------ ForSmallSampleGenTest
        - MC test using gen jets 
        - Purpose is to study the correlation between MJ(vec) and MJ(scalar)

/Dev 

/Processing
    |------ FastjetMaking
        |------ Local
        |------ Batch
            - For batch submission to UCSB HEP cluster
            - Be careful about the directory structure because the whole Fastjet
              directory is sent to a node
            - Needs a script that does resubmission of failed jobs
    |------ Slimming
        - This is Batch submission 
        - Resubmission script is under development 
    |------ BabyMaking
        |------ Batch
            - Now there are two sets of files (for slimmed and non-slimmed cfAs) 
              because some of the necessary cfA samples were not slimmed
            - Has resubmission script 

