import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HLTMiniAOD")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.MessageLogger.destinations = ['cout', 'cerr']

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
options.register('reportEvery',
                1,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                "Report every N events for cout/cerr")
options.register('debug',
                False,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.bool,
                "Print debug info")

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
print " Debug mode            = %s" % options.debug

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvent) )


process.source = cms.Source("PoolSource",
     fileNames = readFiles 
)

process.TFileService=cms.Service("TFileService",
        fileName=cms.string(options.outputFile),
        closeFileFast = cms.untracked.bool(True)
)
process.mini = cms.EDAnalyzer('MiniAODanalyzer',
  triggerTag             = cms.InputTag("TriggerResults","","HLT"),
  jetTag                 = cms.InputTag("slimmedJets"),
  prunedGenParticles     = cms.InputTag("prunedGenParticles"),
  packedGenParticles     = cms.InputTag("packedGenParticles"),
  DEBUG = cms.untracked.bool(options.debug)
)


process.p = cms.Path(process.mini)
