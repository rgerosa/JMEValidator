import FWCore.ParameterSet.Config as cms

# Common parameters used in all modules
CommonParameters = cms.PSet(
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
 
process = cms.Process("JRA")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.GlobalTag.globaltag = "MCRUN2_74_V7::All"


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
                                    
inputFiles = cms.untracked.vstring(
        '/store/relval/CMSSW_7_4_1/RelValFS_TTbar_13_PUAVE35/MINIAODSIM/PU25ns_MCRUN2_74_V9_FastSim-v1/00000/1868AA47-19ED-E411-9D57-0025905A6080.root'
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource", fileNames = inputFiles )

# Services
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName = cms.string('output.root')

process.out = cms.OutputModule("PoolOutputModule",
        outputCommands  = cms.untracked.vstring(),
        fileName       = cms.untracked.string("output_edm.root")
        )

# Create all needed jets collections

# jetsCollections is a dictionnary containing all the informations needed for creating a new jet collection. The format used is :
#  "name": {
#      "algo": string ; the jet algorithm to use
#      "pu_methods" : array of strings ; which PU method to use
#      "jec_payloads" : array of strings ; which JEC payload to use for making the JEC. The size must match the size of pu_methods
#      "jec_levels" : array of strings ; which JEC levels to apply
#  }

jetsCollections = {
        'AK4': {
            'algo': 'ak4',
            'pu_methods': ['Puppi', 'CHS', ''],
            'jec_payloads': ['AK4PFPUPPI', 'AK4PFchs', 'AK4PF'],
            'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
            },

        'AK8': {
            'algo': 'ak8',
            'pu_methods': ['Puppi', 'CHS', ''],
            'jec_payloads': ['AK8PFPUPPI', 'AK8PFchs', 'AK8PF'],
            'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
            },
        }

from JMEAnalysis.JetToolbox.jetToolbox_cff import *

for name, params in jetsCollections.items():
    for index, pu_method in enumerate(params['pu_methods']):
        # Add the jet collection
        jetToolbox(process, params['algo'], 'dummy', 'out', PUMethod = pu_method, JETCorrPayload = params['jec_payloads'][index], JETCorrLevels = params['jec_levels'])


# Configure the analyzers

process.jmfw_analyzers = cms.Sequence()

for name, params in jetsCollections.items():
    for index, pu_method in enumerate(params['pu_methods']):

        algo = params['algo'].upper()
        jetCollection = 'selectedPatJets%sPF%s' % (algo, pu_method)

        print('Adding analyzer for jets collection \'%s\'' % jetCollection)

        analyzer = cms.EDAnalyzer('JetMETAnalyzer',
                CommonParameters,
                JetCorLabel   = cms.string(params['jec_payloads'][index]),
                JetCorLevels  = cms.vstring(params['jec_levels']),
                srcJet        = cms.InputTag(jetCollection),
                srcRho        = cms.InputTag('fixedGridRhoAllFastjet'),
                srcVtx        = cms.InputTag('offlineSlimmedPrimaryVertices'),
                srcMuons      = cms.InputTag('selectedPatMuons')
                )

        setattr(process, 'jmfw_%s' % params['jec_payloads'][index], analyzer)
        process.jmfw_analyzers += analyzer

process.puppiReader = cms.EDAnalyzer("puppiAnalyzer",
                                        treeName = cms.string("puppiTree"),
										maxEvents = cms.int32(1000),
                                        nAlgos = cms.InputTag("puppi", "PuppiNAlgos", "JRA"),
                                        rawAlphas = cms.InputTag("puppi", "PuppiRawAlphas", "JRA"),
                                        alphas = cms.InputTag("puppi", "PuppiAlphas", "JRA"),
                                        alphasMed = cms.InputTag("puppi", "PuppiAlphasMed", "JRA"),
                                        alphasRms = cms.InputTag("puppi", "PuppiAlphasRms", "JRA"),
                                        packedPFCandidates = cms.InputTag("packedPFCandidates", "", "PAT")
									)

process.p = cms.Path( process.puppiReader + process.jmfw_analyzers )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.out) 

#!
#! THAT'S ALL! CAN YOU BELIEVE IT? :-D
#!
