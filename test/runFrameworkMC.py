import sys
if not hasattr(sys, 'argv'):
  sys.argv = ["cmsRun", "runFrameworkMC.py"]

import FWCore.ParameterSet.Config as cms
## parse some parameters from external line                                                                                                                                     
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC
options.register ('isMC',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to indicate data or MC');
## processName
options.register ('processName','JRA',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');
## conditions
options.register ('globalTag',"MCRUN2_74_V9",VarParsing.multiplicity.singleton,VarParsing.varType.string,'input global tag to be used');
options.register ('muonIsoCone',0.4,VarParsing.multiplicity.singleton, VarParsing.varType.float,  'value to be used for muon isolation cone');

## Lepton ID
options.register ('muonTypeID',    "Tight",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'muon ID to be considered for MVA PUPPET analysis ');
options.register ('electronTypeID',"Medium", VarParsing.multiplicity.singleton, VarParsing.varType.string, 'electron ID to be considered for MVA PUPPET analysis ');
options.register ('tauTypeID',     "Loose",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'tau ID to be considered for MVA PUPPET analysis ');
## selections
options.register ('applyZSelections',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'apply selection for Zll events when runMVAPUPPETAnalysis is true');
options.register ('jetPtCut',1,VarParsing.multiplicity.singleton, VarParsing.varType.float, 'apply a jet pt cut for mva met input');
## JEC
options.register ('useJECFromDB',         False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'read JEC from the database for special JEC not in GT');
options.parseArguments()

## import the function to create the process
from JMEAnalysis.JMEValidator.FrameworkConfiguration import createProcess


## create the process with all the information
process = createProcess(options.isMC, ## MC or data
                        options.processName,
                        options.globalTag, ## GT
                        options.muonTypeID, options.muonIsoCone,## muons
                        options.electronTypeID, ## electrons
                        options.tauTypeID,## taus
                        options.applyZSelections,
                        options.jetPtCut,
                        options.useJECFromDB ## JEC
                        );

####### files
inputFiles = []
if options.isMC == True:
     inputFiles.append('/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/78750D2F-726D-E511-A7F7-0025905C2CEA.root')
elif options.isMC == False:
     inputFiles.append('/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/244/00000/E42FEF61-6E27-E511-B93A-02163E0143C0.root')


## set input files
process.source = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring(inputFiles);

## count the number of events
process.AllEvents = cms.EDFilter("PassFilter",
    srcGenEventInfo   = cms.InputTag("generator"),
    isMC      = cms.bool(options.isMC)                              
  )
process.counterPath = cms.Path(process.AllEvents)
## set max events 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
## output name
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName = cms.string('output.root')
process.TFileService.closeFileFast = cms.untracked.bool(True)


## logger
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 50

#! Output and Log                                                                                                                                                            
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(True)
 
"""
    process.output = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string('output_particles.root'),
                                      outputCommands = cms.untracked.vstring('keep *_slimmedMuons'+options.muonTypeID+'*_*_*',
                                                                             'keep *_slimmedElectrons'+options.electronTypeID+'*_*_*',
                                                                             'keep *_*selectedPatJets*_*_*',
                                                                             'keep *_*selectedPatJets*Cleaned*_*_*',
                                                                             'keep *_slimmed*Taus*'+options.tauTypeID+'*_*_*',
                                                                             'keep *_slimmed*MET*_*_*',
                                                                             'keep *_*mvaPUPPET*_*_*',
                                                                             'keep *_*recoil*_*_*',
                                                                             'keep *_*ZdiLepton*_*_*',
                                                                             'keep *_*LeptonMerge*_*_*',
                                                                             'keep *_*ZtagBoson*_*_*',
                                                                             'keep *_*ak4GenJetsNoNu*_*_*',
                                                                             'keep *_*offlineSlimmedPrimaryVertices*_*_*',
                                                                             'keep *_*packedGenLeptons*_*_*',
                                                                             'keep *_*puppi*_*_*',
                                                                             'keep *_*pfAllChargedParticlesPuppi*_*_*',
                                                                             'keep *_*pfPileUpIso*_*_*',
                                                                             'keep *_*pfNoPileUpIso*_*_*',
                                                                             'keep *_*pfAllNeutralParticlesPuppi*_*_*',
                                                                             'keep *_*pupuppi*_*_*',
                                                                             'keep *_*pfPuppi*_*_*',
                                                                             'keep *_*pfPUPuppi*_*_*',
                                                                             'keep *_*pfPUPuppiCharge*_*_*',
                                                                             'keep *_*pfAllNeutralParticlesPuppiPU*_*_*',
                                                                             'keep *_*pfChargedPV*_*_*',
                                                                             'keep *_*pfNeutrals*_*_*',
                                                                             'keep *_*pfChargedPU*_*_*',
                                                                             'keep *_*pfCandidatesForMET*_*_*',
                                                                             'keep *_*pfCandidatesForMETCH*_*_*',
                                                                             'keep *_*packed*Candidates*_*_*',
                                                                             'keep *_*Met*_*_*',
                                                                             'keep *_*TriggerResults*_*_*',
                                                                             'keep *_*neutral*_*_*',
                                                                             ),        
                                      SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring('p'))
                                      )
    
    process.out = cms.EndPath(process.output)
    

"""
processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
