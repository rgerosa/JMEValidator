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
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.GlobalTag.globaltag = "MCRUN2_74_V7::All"


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

infiles = cms.untracked.vstring(
#######
# QCD #
#######
    '/store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/MCRUN2_74_V7-v1/00000/5CCC0483-C9D1-E411-BD1E-0025905A60DE.root',
###########
# DY Jets #
###########
	# '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
	# '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root',
    )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource", fileNames = infiles )

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Run PUPPI, make some new jet collections
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from JMEAnalysis.JetToolbox.jetToolbox_cff import *
jetToolbox( process, 'ak4', 'ak4JetSubs', 'out', PUMethod='Puppi', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
jetToolbox( process, 'ak4', 'ak4JetSubs', 'out', PUMethod='SK', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
#jetToolbox( process, 'ak4', 'ak4JetSubs', 'out', PUMethod='CS', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
jetToolbox( process, 'ak4', 'ak4JetSubs', 'out', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']) # CHS jets?
jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', PUMethod='Puppi', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', PUMethod='SK', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
#jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', PUMethod='CS', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']) # CHS jets?

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('test.root')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! JME stuff (analyzer)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

jetCollections = [];
correctionLevels = [];
jetSrcName = [];

jetCollections.append('AK4PFchs');
correctionLevels.append(['L1FastJet']);
jetSrcName.append('selectedPatJetsAK4PFCHS');

jetCollections.append('AK8PFchs');
correctionLevels.append(['L1FastJet']);
jetSrcName.append('selectedPatJetsAK8PFCHS');

jetCollections.append('AK4PUPPI');
correctionLevels.append([]);
jetSrcName.append('selectedPatJetsAK4PFPuppi');

jetCollections.append('AK8PUPPI');
correctionLevels.append([]);
jetSrcName.append('selectedPatJetsAK8PFPuppi');


jetCollections.append('AK4SK');
correctionLevels.append([]);
jetSrcName.append('selectedPatJetsAK4PFSK');

jetCollections.append('AK8SK');
correctionLevels.append([]);
jetSrcName.append('selectedPatJetsAK8PFSK');

#jetCollections.append('AK4CS');
#correctionLevels.append([]);
#jetSrcName.append('selectedPatJetsAK4PFCS');

#jetCollections.append('AK8CS');
#correctionLevels.append([]);
#jetSrcName.append('selectedPatJetsAK8PFCS');


validator_sequence = cms.Sequence()
setattr(process,"validator_sequence",validator_sequence)
for i in range(len(jetCollections)):
	pnm = cms.EDAnalyzer('validatorTreeMaker',
			     PileupNtupleMakerParameters,
			     JetCorLabel       = cms.string(jetCollections[i]),
			     JetCorLevels      = cms.vstring(correctionLevels[i]),
			     srcJet            = cms.InputTag(jetSrcName[i]),
			     srcRho            = cms.InputTag('fixedGridRhoAllFastjet'),
			     srcVtx            = cms.InputTag('offlineSlimmedPrimaryVertices'),
			     srcMuons          = cms.InputTag('selectedPatMuons')
			     )

	#process.myseq = cms.Sequence(process.pnm)
	setattr(process,'nt_'+jetCollections[i],pnm)
	validator_sequence = cms.Sequence(validator_sequence*pnm)

# process.p = cms.Path( puppi_onMiniAOD * corrservices_sequence * conversion_sequence * validator_sequence );

process.puppiReader = cms.EDAnalyzer("puppiReader")

process.p = cms.Path( process.puppiReader + validator_sequence )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
                                  #outputCommands = cms.untracked.vstring('drop *','keep *_puppi_*_*'),
                                  outputCommands = cms.untracked.vstring('keep *'),
                                  fileName       = cms.untracked.string ("Output.root")                                                                                                                   
)
# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.out) 

#!
#! THAT'S ALL! CAN YOU BELIEVE IT? :-D
#!
