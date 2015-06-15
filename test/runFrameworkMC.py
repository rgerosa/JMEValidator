import FWCore.ParameterSet.Config as cms
## parse some parameters from external line                                                                                                                                     
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('globalTag',        "MCRUN2_74_V9",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'input global tag to be used');
options.register ('isMC'     ,        True,            VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate data or MC');
options.register ('runPuppiMuonIso',  True,            VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate to run or not puppi iso for mons');
options.register ('muonIsoCone',      0.4,             VarParsing.multiplicity.singleton, VarParsing.varType.float,  'value to be used for muon isolation cone');
options.register ('dropAnalyzerDumpEDM', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'do not run the analyzer and store an edm file');
options.register ('runMVAPUPPETAnalysis', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'run a specific analysis for MVA MET : Z->LL events');
options.register ('muonTypeID',       "Tight",         VarParsing.multiplicity.singleton, VarParsing.varType.string, 'muon ID to be considered for MVA PUPPET analysis ');
options.register ('electronTypeID',   "Medium",        VarParsing.multiplicity.singleton, VarParsing.varType.string, 'electron ID to be considered for MVA PUPPET analysis ');
options.register ('tauTypeID',        "Loose",         VarParsing.multiplicity.singleton, VarParsing.varType.string, 'tau ID to be considered for MVA PUPPET analysis ');
options.parseArguments()

## import the function to create the process
from JMEAnalysis.JMEValidator.FrameworkConfiguration import createProcess

process = createProcess(options.isMC, options.globalTag, options.muonTypeID, options.runPuppiMuonIso, options.muonIsoCone, 
                        options.electronTypeID, options.tauTypeID,options.dropAnalyzerDumpEDM, options.runMVAPUPPETAnalysis)

if len(options.inputFiles) == 0 and options.isMC == True:
    options.inputFiles.append('/store/relval/CMSSW_7_4_4/RelValZMM_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_38Tbis-v1/00000/64DC187C-5D09-E511-A14C-0025905964B2.root');
#    options.inputFiles.append('/store/relval/CMSSW_7_4_4/RelValZEE_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_38Tbis-v1/00000/D4F6F957-4809-E511-AD3B-00261894386F.root');


## set input files
process.source.fileNames = cms.untracked.vstring(options.inputFiles);
## set max events 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
## output name
if not options.dropAnalyzerDumpEDM :
   if options.outputFile == "" or options.outputFile == "output.root" :
       process.TFileService.fileName = cms.string('output_mc.root') if options.isMC else cms.string('output_data.root')
   else:    
       process.TFileService.fileName = options.outputFile

## logger
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 50

#! Output and Log                                                                                                                                                            
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(True)
 
if options.dropAnalyzerDumpEDM :
    process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string('output.root'),
                                      outputCommands = cms.untracked.vstring('keep *_slimmedMuons*_*_*',
                                                                             'keep *_slimmedElectrons*_*_*',
                                                                             'keep *_patJets*Cleaned*_*_*',
                                                                             'keep *_slimmed*Tau*_*_*',
                                                                             'keep *_slimmed*MET*_*_*',
                                                                             'keep *_*mvaPUPPET*_*_*',
                                                                             'keep *_*recoil*_*_*',
                                                                             'keep *_*ZtagBoson*_*_*')
                                      )
    
    process.out = cms.EndPath(process.output)
    
