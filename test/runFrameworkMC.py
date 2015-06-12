import FWCore.ParameterSet.Config as cms

## import the function to create the process
from JMEAnalysis.JMEValidator.FrameworkConfiguration import createProcess

## parse some parameters from external line                                                                                                                                     
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('globalTag',   "MCRUN2_74_V9",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'input global tag to be used');
options.register ('isMC'     ,   True,            VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate data or MC');
options.register ('runPuppiMuonIso',   True,      VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate to run or not puppi iso for mons');
options.register ('muonIsoCone',     0.4,         VarParsing.multiplicity.singleton, VarParsing.varType.float,  'value to be used for muon isolation cone');
options.register ('muonCollection',  "slimmedMuons",  VarParsing.multiplicity.singleton, VarParsing.varType.string,  'default benchmark of muons to be considered');
options.register ('electronCollection',  "slimmedElectrons",  VarParsing.multiplicity.singleton, VarParsing.varType.string,  'default benchmark of electrons to be considered');
options.register ('tauCollection',  "slimmedTaus",  VarParsing.multiplicity.singleton, VarParsing.varType.string,  'default benchmark of taus to be considered');
options.register ('dropAnamyzerDumpEDM', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'do not run the analyzer and store an edm file');
options.parseArguments()

## create the process instance
process = createProcess(options.isMC, options.globalTag, options.muonCollection, options.runPuppiMuonIso, options.muonIsoCone, 
                        options.electronCollection, options.tauCollection,options.dropAnamyzerDumpEDM)

if len(options.inputFiles) == 0 and options.isMC == True:
    options.inputFiles.append('/store/relval/CMSSW_7_4_1/RelValTTbarLepton_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/7A4C3E6D-E7EC-E411-8CBE-0025905A60FE.root');


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

#! Output and Log                                                                                                                                                            
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(True)
 

if options.dropAnamyzerDumpEDM :

    process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string('output.root'),
                                      outputCommands = cms.untracked.vstring('keep *_slimmedMuons*_*_*',
                                                                             'keep *_slimmedElectrons*_*_*',
                                                                             'keep *_patJets*Cleaned*_*_*',
                                                                             'keep *_slimmed*Tau*_*_*',
                                                                             'keep *_slimmed*MET*_*_*', 
                                                                             'keep *_mvaPUPPET_Z_*',
                                                                             'keep *_mvaPUPPET_recoilsForMvaPUPPET_*',
                                                                             'keep *_mvaPUPPET__*') )
    
    process.out = cms.EndPath(process.output)
    
