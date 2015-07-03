import FWCore.ParameterSet.Config as cms
<<<<<<< HEAD
=======


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
>>>>>>> origin
 
def runMuonIsolation(process):
	
	process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
	process.load("Configuration.EventContent.EventContent_cff")
<<<<<<< HEAD
	process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
	process.load('Configuration.StandardSequences.MagneticField_38T_cff')
	process.GlobalTag.globaltag = "MCRUN2_74_V7::All"
=======
	process.load('Configuration.StandardSequences.Geometry_cff')
	process.load('Configuration.StandardSequences.MagneticField_38T_cff')
	process.GlobalTag.globaltag = "PHYS14_25_V2::All"
>>>>>>> origin
	
	### load default PAT sequence
	process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
	process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
	process.load("JMEAnalysis.JMEValidator.convertPackedCandToRecoCand_cff")
	
	process.packedPFCandidatesWoMuon  = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV>=2 && abs(pdgId)!=13 " ) )
	convertedPackedPFCandidatesWoMuon = cms.EDProducer('convertCandToRecoCand',
							   src = cms.InputTag('packedPFCandidatesWoMuon')
							   )
	
	setattr( process, 'convertedPackedPFCandidatesWoMuon', convertedPackedPFCandidatesWoMuon )
<<<<<<< HEAD

	# presequences needed for PUPPI and PF-Weighting
	process.patseq = cms.Sequence(process.convertedPackedPFCandidates *
				      process.convertedPackedPFCandidatesWoMuon *
				      process.pfParticleSelectionForIsoSequence *
				      process.selectedPatMuons
				      )

	process.p = cms.Path(process.patseq)
	
	# # change the input collections
	process.particleFlowPtrs.src = 'convertedPackedPFCandidates'
	process.pfPileUpIsoPFBRECO.Vertices = 'offlineSlimmedPrimaryVertices'
	process.pfPileUpPFBRECO.Vertices    = 'offlineSlimmedPrimaryVertices'
=======
	process.patseq = cms.Sequence(process.convertedPackedPFCandidates *
				      convertedPackedPFCandidatesWoMuon *
				      process.patCandidates * process.selectedPatCandidates)
	process.p = cms.Path(process.patseq)
	
	# change the input collections
	process.particleFlowPtrs.src = 'convertedPackedPFCandidates'
	process.pfPileUpIso.Vertices = 'offlineSlimmedPrimaryVertices'
	process.pfPileUp.Vertices    = 'offlineSlimmedPrimaryVertices'
	
	# remove unnecessary PAT modules
	process.p.remove(process.makePatElectrons)
	process.p.remove(process.makePatPhotons)
	process.p.remove(process.makePatJets)
	process.p.remove(process.makePatTaus)
	process.p.remove(process.makePatMETs)
	process.p.remove(process.patCandidateSummary)
	process.p.remove(process.selectedPatElectrons)
	process.p.remove(process.selectedPatPhotons)
	process.p.remove(process.selectedPatJets)
	process.p.remove(process.selectedPatTaus)
	process.p.remove(process.selectedPatCandidateSummary)
>>>>>>> origin
	
	### muon selection
	process.selectedPatMuons.src = 'slimmedMuons'
	process.selectedPatMuons.cut = 'pt>10 && abs(eta)<2.4'
	
	# load user-defined particle collections (e.g. PUPPI)
	
	# -- PF-Weighted
	process.load('CommonTools.ParticleFlow.deltaBetaWeights_cff')
<<<<<<< HEAD
	process.pfWeightedPhotons.src = 'pfAllPhotonsPFBRECO'
	process.pfWeightedPhotons.chargedFromPV = 'pfAllChargedParticlesPFBRECO'
	process.pfWeightedPhotons.chargedFromPU = 'pfPileUpAllChargedParticlesPFBRECO'
	process.pfWeightedNeutralHadrons.src = 'pfAllNeutralHadronsPFBRECO'
	process.pfWeightedNeutralHadrons.chargedFromPV = 'pfAllChargedParticlesPFBRECO'
	process.pfWeightedNeutralHadrons.chargedFromPU = 'pfPileUpAllChargedParticlesPFBRECO'

=======
	
>>>>>>> origin
	# -- PUPPI
	from JMEAnalysis.JMEValidator.pfPUPPISequence_cff import *
	load_pfPUPPI_sequence(process, 'pfPUPPISequence', algo = 'PUPPI',
	  src_puppi = 'pfAllHadronsAndPhotonsForPUPPI',
	  cone_puppi_central = 0.5
	)
	
	# change the input collections
	process.pfAllHadronsAndPhotonsForPUPPI.src = 'convertedPackedPFCandidates'
	process.particleFlowPUPPI.candName = 'packedPFCandidates'
	process.particleFlowPUPPI.vertexName = 'offlineSlimmedPrimaryVertices'
	
	# -- PUPPI isolation calculation without muon
	load_pfPUPPI_sequence(process, 'pfNoMuonPUPPISequence', algo = 'NoMuonPUPPI',
	  src_puppi = 'pfAllHadronsAndPhotonsForNoMuonPUPPI',
	  cone_puppi_central = 0.5
	)
	process.pfAllHadronsAndPhotonsForNoMuonPUPPI.src = 'convertedPackedPFCandidatesWoMuon'
	process.particleFlowNoMuonPUPPI.candName         = 'packedPFCandidatesWoMuon'
	process.particleFlowNoMuonPUPPI.vertexName       = 'offlineSlimmedPrimaryVertices'
	
<<<<<<< HEAD
	process.ParticleIsoSequences = cms.Sequence(process.pfDeltaBetaWeightingSequence * 
												process.pfPUPPISequence * 
												process.pfNoMuonPUPPISequence
												)
=======
	from JMEAnalysis.JMEValidator.makePUPPIJets_cff import *
	load_PUPPIJet_sequence(process,"PUPPIJetSequence",[0.4,0.8])
	
	process.p.replace(
	  process.pfParticleSelectionSequence,
	  process.pfParticleSelectionSequence  *
	  process.pfDeltaBetaWeightingSequence *
	  process.pfPUPPISequence *
	  process.pfNoMuonPUPPISequence *
	  process.PUPPIJetSequence *
	  process.patJetPartonMatchAK4PUPPIJets *
	  process.patJetGenJetMatchAK4PUPPIJets *
	  process.patJetsAK4PUPPIJets *
	  process.patJetPartonMatchAK8PUPPIJets *
	  process.patJetGenJetMatchAK8PUPPIJets *
	  process.patJetsAK8PUPPIJets
	)
	
>>>>>>> origin
	
	from JMEAnalysis.JMEValidator.MuonPFIsolationSequence_cff import *
	muon_src, cone_size = 'selectedPatMuons', 0.4
	
	process.pfCHLVForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV>=2 && abs(charge) > 0"))
	process.pfCHPUForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV<=1 && abs(charge) > 0"))
	process.pfPhotonsForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("pdgId==22"))
	process.pfNHForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("pdgId!=22 && abs(charge) == 0" ))
	
	load_muonPFiso_sequence(process, 'MuonPFIsoSequenceSTAND', algo = 'R04STAND',
	  src = muon_src,
	  src_charged_hadron = 'pfCHLVForIso',
	  src_neutral_hadron = 'pfNHForIso',
	  src_photon         = 'pfPhotonsForIso',
	  src_charged_pileup = 'pfCHPUForIso',
	  coneR = cone_size
	)
	
	load_muonPFiso_sequence(process, 'MuonPFIsoSequencePFWGT', algo = 'R04PFWGT',
	  src = muon_src,
	  src_neutral_hadron = 'pfWeightedNeutralHadrons',
	  src_photon         = 'pfWeightedPhotons',
	  coneR = cone_size
	)
	
	load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPI', algo = 'R04PUPPI',
	  src = muon_src,
	  src_charged_hadron = 'pfPUPPIChargedHadrons',
	  src_neutral_hadron = 'pfPUPPINeutralHadrons',
	  src_photon         = 'pfPUPPIPhotons',
	  coneR = cone_size
	)
	
	load_muonPFiso_sequence(process, 'MuonPFIsoSequenceNoMuonPUPPI', algo = 'R04NOMUONPUPPI',
	  src = muon_src,
	  src_charged_hadron = 'pfNoMuonPUPPIChargedHadrons',
	  src_neutral_hadron = 'pfNoMuonPUPPINeutralHadrons',
	  src_photon         = 'pfNoMuonPUPPIPhotons',
	  coneR = cone_size
	)
	
<<<<<<< HEAD
	# process.muPFIsoDepositCharged.src = 'slimmedMuons'
=======
	process.muPFIsoDepositCharged.src = 'slimmedMuons'
>>>>>>> origin
	process.muonMatch.src = 'slimmedMuons'
	process.muonMatch.matched = 'prunedGenParticles'
	process.patMuons.pvSrc = 'offlineSlimmedPrimaryVertices'
	process.p.remove(process.patMuons)
	
	process.MuonPFIsoSequences = cms.Sequence(
	  process.MuonPFIsoSequenceSTAND *
	  process.MuonPFIsoSequencePFWGT *
	  process.MuonPFIsoSequencePUPPI *
	  process.MuonPFIsoSequenceNoMuonPUPPI
	)
	
	process.p.replace(
	  process.selectedPatMuons,
	  process.selectedPatMuons *
<<<<<<< HEAD
	  process.ParticleIsoSequences *
	  process.MuonPFIsoSequences
	)
	
=======
	  process.MuonPFIsoSequences
	)
	
	process.leptonsAndMET = cms.EDAnalyzer("LeptonsAndMETAnalyzer",
	                                       srcIsoMuons = cms.InputTag("selectedMuonsForZ"),
	                                       srcMET = cms.InputTag("slimmedMETs"),
	                                       srcPUPPET = cms.InputTag("pfMetPuppi"),
	                                       srcVtx            = cms.InputTag('offlineSlimmedPrimaryVertices'),
	                                       srcMuons          = cms.InputTag( muon_src ),
	                                       srcVMCHSTAND      = cms.InputTag('muPFIsoValueCHR04STAND'),
	                                       srcVMNHSTAND      = cms.InputTag('muPFIsoValueNHR04STAND'),
	                                       srcVMPhSTAND      = cms.InputTag('muPFIsoValuePhR04STAND'),
	                                       srcVMPUSTAND      = cms.InputTag('muPFIsoValuePUR04STAND'),
	                                       srcVMCHPFWGT      = cms.InputTag('muPFIsoValueCHR04PFWGT'),
	                                       srcVMNHPFWGT      = cms.InputTag('muPFIsoValueNHR04PFWGT'),
	                                       srcVMPhPFWGT      = cms.InputTag('muPFIsoValuePhR04PFWGT'),
	                                       srcVMCHPUPPI      = cms.InputTag('muPFIsoValueCHR04PUPPI'),
	                                       srcVMNHPUPPI      = cms.InputTag('muPFIsoValueNHR04PUPPI'),
	                                       srcVMPhPUPPI      = cms.InputTag('muPFIsoValuePhR04PUPPI'),
	                                       srcVMCHNOMUONPUPPI      = cms.InputTag('muPFIsoValueCHR04NOMUONPUPPI'),
	                                       srcVMNHNOMUONPUPPI      = cms.InputTag('muPFIsoValueNHR04NOMUONPUPPI'),
	                                       srcVMPhNOMUONPUPPI      = cms.InputTag('muPFIsoValuePhR04NOMUONPUPPI')
	                                       )
	
>>>>>>> origin
