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
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.GlobalTag.globaltag = "PHYS14_25_V2::All"


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
                                    
inputFiles = cms.untracked.vstring(
        # 'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/1020E374-B26B-E411-8F91-E0CB4E29C513.root'
        'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
        # 'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/Neutrino_Pt-2to20_gun/MINIAODSIM/AVE20BX25_tsg_PHYS14_25_V3-v1/00000/007B75A0-538F-E411-B87D-00259059649C.root',
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource", fileNames = inputFiles )

# Services
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
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

from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMetPuppi = pfMet.clone();
process.pfMetPuppi.src = cms.InputTag('puppi')

process.selectedMuonsForZ = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
   (pfIsolationR04().sumChargedHadronPt+
    max(0.,pfIsolationR04().sumNeutralHadronEt+
    pfIsolationR04().sumPhotonEt-
    0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))

process.leptonsAndMET = cms.EDAnalyzer("LeptonsAndMETAnalyzer",
                                        srcIsoMuons = cms.InputTag("selectedMuonsForZ"),
                                        srcMET = cms.InputTag("slimmedMETs"),
                                        srcPUPPET = cms.InputTag("pfMetPuppi")                                                                                
                                    )

process.p = cms.Path( process.pfMetPuppi + process.selectedMuonsForZ + process.puppiReader + process.leptonsAndMET + process.jmfw_analyzers )


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
