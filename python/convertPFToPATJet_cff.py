import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

#from PhysicsTools.PatAlgos.recoLayer0.jetCorrections_cff import *
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
#from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *


def convertPFToPATJet(proc, inputColl, outputColl, alg, r, corrAlgo, corrLevels):

	addJetCollection(
		proc,
		labelName = outputColl,
		jetSource = cms.InputTag(inputColl),
		algo = alg,
		rParam = r,
		# jetCorrections = None,
		trackSource = cms.InputTag('unpackedTracksAndVertices'),
		pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
		jetCorrections = (corrAlgo, cms.vstring(corrLevels), 'None')
		#btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
	)

	for module in [getattr(proc,'patJets'+outputColl)]:
		module.addJetCharge = False
		module.addBTagInfo = False
		module.getJetMCFlavour = False
		module.addAssociatedTracks = False
		module.addGenPartonMatch = False
		module.addGenJetMatch = True
		
		if 'Puppi' in inputColl: module.addJetCorrFactors = False
		else: module.addJetCorrFactors = True

	getattr(proc, 'patJetGenJetMatch'+inputColl).matched = cms.InputTag(alg.upper()+'GenJets');

