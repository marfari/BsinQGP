****************************
          README
****************************

area setup:  
source /cvmfs/cms.cern.ch/cmsset_default.sh  
export SCRAM_ARCH=slc6_amd64_gcc700  
cmsrel CMSSW_10_3_1_patch3  
cd CMSSW_10_3_1_patch3/src  
mkdir UserCode  
git clone [repo]  
mkdir prefiltered_trees  

to compile: scram b  

First run selection_skimmer.cc!  

It runs over data and mc samples in lstore  
right now for B+  
to run over other samples, change inputs  
It loops over the data or mc samples, flattening the tree, applying the selection, and writing a new tree with less variables  

running example: selection_skimmer --colsyst PbPb --mc 0 --channel 1  

Then you can run sideban_subtraction.cc over the outputs of selection_skimmer that are in prefiltered_trees  
It applies a sideband subtraction method to the data and then compares the result with the monte carlo distributions  

First do:  
mkdir mc_validation  
mkdir mc_validation/sideband_sub  

To run:  
sideband_sub

Some important git commands:


1. git add [files]  

2. git commit -m ["message"]  

3. git push origin master (vai online)  

Create remote repository and import files  

1. git remote add [repository name] [repository address]  

2. git pull [repository name] 



