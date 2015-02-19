#! Text To ASCII from http://patorjk.com/software/taag/#p=display&f=Big&t=MUON%20ISOLATION

import FWCore.ParameterSet.Config as cms

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Jet and reference kinematics
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PileupNtupleMakerParameters = cms.PSet(
    # record flavor information, consider both RefPt and JetPt
    doComposition   = cms.bool(True),
    doFlavor        = cms.bool(True),
    doRefPt         = cms.bool(True),
    doJetPt         = cms.bool(True),
    # MATCHING MODE: deltaR(ref,jet)
    deltaRMax       = cms.double(99.9),
    # deltaR(ref,parton) IF doFlavor is True
    deltaRPartonMax = cms.double(0.25),
    # consider all matched references
    nJetMax         = cms.uint32(0),
)
 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Process
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process = cms.Process("JRA")


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Conditions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "PHYS14_25_V2::All"


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
                                    
dyFiles = cms.untracked.vstring(
#######
# QCD #
#######
	'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/1020E374-B26B-E411-8F91-E0CB4E29C513.root',
#	'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/1EF51024-986B-E411-A6F6-20CF300E9EAF.root',
#	'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/5AE5B0FC-986B-E411-ACED-20CF3027A57B.root',
###########
# DY Jets #
###########
	# '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
	# '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root',
    )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.source = cms.Source("PoolSource", fileNames = dyFiles )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('test.root')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Run PUPPI, make some new jet collections
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# from CommonTools.PileupAlgos.Puppi_cff import puppi
process.load('CommonTools.PileupAlgos.Puppi_cff');
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

process.load('JMEAnalysis.JMEValidator.makeJets_cff')
puppi_onMiniAOD = cms.Sequence(process.puppi * process.pfCHS * process.AK4GenJets * process.AK8GenJets * process.AK4PFchsJets * process.AK8PFchsJets * process.AK4PFJetsPuppi * process.AK8PFJetsPuppi)
setattr(process,'puppi_onMiniAOD',puppi_onMiniAOD)

#! convert the PUPPI jets into pat::jets
from JMEAnalysis.JMEValidator.convertPFToPATJet_cff import convertPFToPATJet
convertPFToPATJet(process,'AK4PFchsJets','AK4PFchsJets','ak4',0.4,'AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
convertPFToPATJet(process,'AK8PFchsJets','AK8PFchsJets','ak8',0.8,'AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
convertPFToPATJet(process,'AK4PFJetsPuppi','AK4PFJetsPuppi','ak4',0.4,'AK4PFchs', [])
convertPFToPATJet(process,'AK8PFJetsPuppi','AK8PFJetsPuppi','ak8',0.8,'AK8PFchs', [])
conversion_sequence = cms.Sequence(process.patJetsAK4PFchsJets*process.patJetsAK8PFchsJets*process.patJetsAK4PFJetsPuppi*process.patJetsAK8PFJetsPuppi)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! JME stuff (analyzer)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

jetCollections = [];
correctionLevels = [];
jetSrcName = [];

jetCollections.append('AK4PFchs');
correctionLevels.append(['L1FastJet']);
jetSrcName.append('patJetsAK4PFchsJets');

jetCollections.append('AK8PFchs');
correctionLevels.append(['L1FastJet']);
jetSrcName.append('patJetsAK8PFchsJets');

jetCollections.append('AK4PUPPI');
correctionLevels.append([]);
jetSrcName.append('patJetsAK4PFJetsPuppi');

jetCollections.append('AK8PUPPI');
correctionLevels.append([]);
jetSrcName.append('patJetsAK8PFJetsPuppi');

validator_sequence = cms.Sequence()
setattr(process,"validator_sequence",validator_sequence)
for i in range(len(jetCollections)):
	pnm = cms.EDAnalyzer('validatorTreeMaker',
	                    PileupNtupleMakerParameters,
						JetCorLabel       = cms.string(jetCollections[i]),
						JetCorLevels      = cms.vstring(correctionLevels[i]),
						srcJet            = cms.InputTag(jetSrcName[i]),
						srcRho            = cms.InputTag('fixedGridRhoAll'),
						srcVtx            = cms.InputTag('offlineSlimmedPrimaryVertices'),
						srcMuons          = cms.InputTag('selectedPatMuons')
			 )

	#process.myseq = cms.Sequence(process.pnm)
	setattr(process,'nt_'+jetCollections[i],pnm)
	validator_sequence = cms.Sequence(validator_sequence*pnm)

process.p = cms.Path( puppi_onMiniAOD * conversion_sequence * validator_sequence );


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
                                  #outputCommands = cms.untracked.vstring('drop *','keep *_puppi_*_*'),
                                  outputCommands = cms.untracked.vstring('keep *'),
                                  fileName       = cms.untracked.string ("Output.root")                                                                                                                   
)
# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.output) 

#!
#! THAT'S ALL! CAN YOU BELIEVE IT? :-D
#!
