import FWCore.ParameterSet.Config as cms

convertedRecoPFCandidates = cms.EDProducer('convertRecoCandToPackedCand',
                                           src = cms.InputTag('PFCandidates'),
                                           srcFromPVLoose = cms.InputTag("pfNoPileUpJME"),
                                           srcVtx = cms.InputTag("offlineSlimmedPrimaryVertices")
                                           )
