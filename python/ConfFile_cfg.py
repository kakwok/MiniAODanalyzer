import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HLTMiniAOD")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

#------------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------------

options = VarParsing.VarParsing()
options.register('inputFile',
                "input_miniAOD.root",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "Input miniAOD")
options.register('outputFile',
                "ntuple_output.root",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "filename of output root file")

options.register('maxEvent',
                100,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                "max event to process")


options.parseArguments()

inputFileList = options.inputFile.split(",")
readFiles = cms.untracked.vstring()
for f in inputFileList:
    f = "file:"+f
    readFiles.extend([f])
inputFile  = options.inputFile
outputFile = options.outputFile

print " Going to use input    = %s" % inputFile
print " Going to write output = %s" % outputFile

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvent) )


process.source = cms.Source("PoolSource",
     fileNames = readFiles 
)

process.TFileService=cms.Service("TFileService",
        fileName=cms.string(options.outputFile),
        closeFileFast = cms.untracked.bool(True)
)
process.mini = cms.EDAnalyzer('MiniAODanalyzer',
  triggerTag = cms.InputTag("TriggerResults","","HLT"),
  jetTag =  cms.InputTag("slimmedJets"),
  DEBUG = cms.untracked.bool(False)
)


process.p = cms.Path(process.mini)
