import FWCore.ParameterSet.Config as cms

## import the function to create the process
from JMEAnalysis.JMEValidator.FrameworkConfiguration import createProcess

## parse some parameters from external line                                                                                                                                     
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('globalTag',   "MCRUN2_74_V7",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'input global tag to be used');
options.register ('isMC'     ,   True,            VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate data or MC');

options.parseArguments()

## create the process instance
process = createProcess(options.isMC, options.globalTag)

if len(options.inputFiles) == 0 and options.isMC == True:
    options.inputFiles.append('/store/relval/CMSSW_7_4_1/RelValFS_TTbar_13_PUAVE35/MINIAODSIM/PU25ns_MCRUN2_74_V9_FastSim-v1/00000/1868AA47-19ED-E411-9D57-0025905A6080.root');


## set input files
process.source.fileNames = cms.untracked.vstring(options.inputFiles);
## set max events 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
## output name
if options.outputFile == "" or options.outputFile == "output.root" :
    process.TFileService.fileName = cms.string('output_mc.root') if options.isMC else cms.string('output_data.root')
else:    
    process.TFileService.fileName = options.outputFile

## logger
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
