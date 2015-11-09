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
## iPUPPI options
options.register ('runPuppiMuonIso',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate to run or not puppi iso for mons');
options.register ('muonIsoCone',0.4,VarParsing.multiplicity.singleton, VarParsing.varType.float,  'value to be used for muon isolation cone');
## PUPPET analysis
options.register ('runMVAPUPPETAnalysis',True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'run a specific analysis for MVA MET : Z->LL events');
## store and edm to debug
options.register ('dropAnalyzerDumpEDM',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'do not run the analyzer and store an edm file');
## Lepton ID
options.register ('muonTypeID',    "Tight",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'muon ID to be considered for MVA PUPPET analysis ');
options.register ('electronTypeID',"Medium", VarParsing.multiplicity.singleton, VarParsing.varType.string, 'electron ID to be considered for MVA PUPPET analysis ');
options.register ('tauTypeID',     "Loose",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'tau ID to be considered for MVA PUPPET analysis ');
## selections
options.register ('applyZSelections',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'apply selection for Zll events when runMVAPUPPETAnalysis is true');
options.register ('applyWSelections',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'apply selection for Wlnu events when runMVAPUPPETAnalysis is true');
options.register ('jetPtCut',0.,VarParsing.multiplicity.singleton, VarParsing.varType.float, 'apply a jet pt cut for mva met input');
## JEC
options.register ('isRunningOn25ns',      False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'true when running on 25ns and JEC from DB should be red');
options.register ('useJECFromDB',         False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'read JEC from the database for special JEC not in GT');
options.register ('applyJECtoPuppiJets',  False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'apply or not JEC on puppi jets');
## Puppi particles
options.register ('runPuppiDiagnostics',  False,VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'run Puppi diagnostic and store in the output');
## Puppi special setup
options.register ('etaCutForMetDiagnostic', 10.0,VarParsing.multiplicity.singleton, VarParsing.varType.float, 'introduce a cut for the diagnostic of the MET');
options.register ('ptNeutralCut',           [],     VarParsing.multiplicity.list, VarParsing.varType.float, 'ptNeutral cut in each eta bin');
options.register ('ptNeutralCutSlope',      [],  VarParsing.multiplicity.list, VarParsing.varType.float, 'ptNeutral cut in each eta bin');
options.register ('etaBinPuppi',            [],     VarParsing.multiplicity.list, VarParsing.varType.float, 'eta bin for puppi algo');
options.register ('puppiUseCharge',         [], VarParsing.multiplicity.list, VarParsing.varType.bool, 'use charge constraint in puppi algo');
options.register ('ptThresholdForTypeIPuppi', 20.,    VarParsing.multiplicity.singleton, VarParsing.varType.float, 'pt threshold for typeI puppi Met');
options.register ('runPUPPINoLeptons', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'remove charged identified leptons before running puppi');
options.parseArguments()


## import the function to create the process
from JMEAnalysis.JMEValidator.FrameworkConfiguration import createProcess

if options.applyWSelections and options.applyZSelections :
      sys.exit("both options applyZSelections and applyWSelections are set to true --> please check")

## create the process with all the information
process = createProcess(options.isMC, ## MC or data
                        options.processName,
                        options.globalTag, ## GT
                        options.muonTypeID, options.runPuppiMuonIso, options.muonIsoCone, ## muons
                        options.electronTypeID, ## electrons
                        options.tauTypeID,## taus
                        options.dropAnalyzerDumpEDM, ## debug in EDM file
                        options.applyZSelections, options.applyWSelections, ## special flags for PUPPI analysis
                        options.jetPtCut,
                        options.applyJECtoPuppiJets, ## JEC for puppi
                        options.useJECFromDB ## JEC
                        );

####### files
if len(options.inputFiles) == 0 and options.isMC == True:
      options.inputFiles.append('/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/78750D2F-726D-E511-A7F7-0025905C2CEA.root')
elif len(options.inputFiles) == 0 and options.isMC == False:
      options.inputFiles.append('/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/244/00000/E42FEF61-6E27-E511-B93A-02163E0143C0.root')


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
    

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
