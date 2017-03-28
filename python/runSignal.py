import os
import glob

def getCommaString(List):
    return ','.join(map(str,List))

eos="/afs/cern.ch/user/k/kakwok/eos/cms/store/mc/RunIISummer16MiniAODv2"
#dataset="/Spin0_ggPhi12j_g1_75_PseudoScalar_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/"
#dataset="/Spin0_ggPhi12j_g1_250_PseudoScalar_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/"
dataset="/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/"
inputFiles = glob.glob(eos+dataset+"*/*.root")

inputFile  = getCommaString(inputFiles)
outputFile = "%s_nTuple.root"% dataset.split("/")[1]

cmd="cmsRun ConfFile_cfg.py inputFile=%s outputFile=%s maxEvent=-1"%(inputFile,outputFile)
print cmd
os.system(cmd)
