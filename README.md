# MiniAODanalyzer
For Xbb double btag trigger study  
Simple, standalone miniAOD analyzer
1. Set-up CMSSW area, compile the code
```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
mkdir bbtoDijet
cd bbtoDijet
git clone https://github.com/kakwok/MiniAODanalyzer.git
cd MiniAODanalyzer
scram b -j8
```
2. Run analyzer on miniAOD input
For one file  
```
cmsenv
cd python 
cmsRun ConfFile_cfg.py inputFile=MY_MINIAOD.root maxEvent=100    
```
For multipile input file, edit data path in `runSignal.py`, then
```
eosmount MYEOSMOUNTPOINT
python runSignal.py
```
