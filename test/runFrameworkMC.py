import FWCore.ParameterSet.Config as cms
## parse some parameters from external line                                                                                                                                     
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('globalTag',        "MCRUN2_74_V9",  VarParsing.multiplicity.singleton, VarParsing.varType.string, 'input global tag to be used');
options.register ('isMC'     ,        True,            VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate data or MC');
options.register ('runPuppiMuonIso',  False,           VarParsing.multiplicity.singleton, VarParsing.varType.bool,   'flag to indicate to run or not puppi iso for mons');
options.register ('muonIsoCone',      0.4,             VarParsing.multiplicity.singleton, VarParsing.varType.float,  'value to be used for muon isolation cone');
options.register ('dropAnalyzerDumpEDM',  False,       VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'do not run the analyzer and store an edm file');
options.register ('runMVAPUPPETAnalysis', True,        VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'run a specific analysis for MVA MET : Z->LL events');
options.register ('muonTypeID',       "Tight",         VarParsing.multiplicity.singleton, VarParsing.varType.string, 'muon ID to be considered for MVA PUPPET analysis ');
options.register ('electronTypeID',   "Medium",        VarParsing.multiplicity.singleton, VarParsing.varType.string, 'electron ID to be considered for MVA PUPPET analysis ');
options.register ('tauTypeID',        "Loose",         VarParsing.multiplicity.singleton, VarParsing.varType.string, 'tau ID to be considered for MVA PUPPET analysis ');
options.register ('applyZSelections'  , True,          VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'apply selection for Zll events when runMVAPUPPETAnalysis is true');
options.parseArguments()

## import the function to create the process
from JMEAnalysis.JMEValidator.FrameworkConfiguration import createProcess

process = createProcess(options.isMC, options.globalTag, options.muonTypeID, options.runPuppiMuonIso, options.muonIsoCone, 
                        options.electronTypeID, options.tauTypeID,options.dropAnalyzerDumpEDM, options.runMVAPUPPETAnalysis, options.applyZSelections)

if len(options.inputFiles) == 0 and options.isMC == True:

     #options.inputFiles.append('root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/AsymptNoPURawReco_MCRUN2_74_V9A-v3/10000/263601E1-AB15-E511-B132-3417EBE4E882.root');
     options.inputFiles.append('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/38D1C54C-0F02-E511-A54E-AC853D9F5256.root')
     #options.inputFiles.append('/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/02DE3B74-6C08-E511-ABE3-0025905A60D0.root')

#     options.inputFiles.append('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/AsymptNoPUbx25Reco_MCRUN2_74_V9-v3/00000/02FAF8EE-3608-E511-AFC8-0025905A612C.root')
#      options.inputFiles.append('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/AsymptFlat0to50bx50Reco_MCRUN2_74_V9A-v3/00000/023F427F-0E08-E511-A813-0025905A60EE.root')
     #options.inputFiles.append('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/AsymptFlat0to50bx25Reco_MCRUN2_74_V9-v3/10000/0031CCC7-B007-E511-A963-0025905964CC.root')



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
                                      fileName = cms.untracked.string('output_50ns_Zee.root'),
                                      outputCommands = cms.untracked.vstring('keep *_slimmedMuons'+options.muonTypeID+'*_*_*',
                                                                             'keep *_slimmedElectrons'+options.electronTypeID+'*_*_*',
                                                                             'keep *_*selectedPatJets*_*_*',
                                                                             'keep *_*selectedPatJets*Cleaned*_*_*',
                                                                             'keep *_slimmed*Taus*'+options.tauTypeID+'*_*_*',
                                                                             'keep *_slimmed*MET*_*_*',
                                                                             'keep *_*mvaPUPPET*_*_*',
                                                                             'keep *_*recoil*_*_*',
                                                                             'keep *_slimmedJets_*_*',
                                                                             'keep *_slimmedJetsPuppi_*_*',
                                                                             'keep *_*ZdiLepton*_*_*',
                                                                             'keep *_*LeptonMerge*_*_*',
                                                                             'keep *_*ZtagBoson*_*_*',
                                                                             'keep *_*ak4GenJetsNoNu*_*_*',
                                                                             'keep *_*offlineSlimmedPrimaryVertices*_*_*',
                                                                             'keep *_*packedGenLeptons*_*_*'),
                                      SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring('p'))
                                      )
    
    process.out = cms.EndPath(process.output)
    

#processDumpFile = open('processDump.py', 'w')
#print >> processDumpFile, process.dumpPython()
