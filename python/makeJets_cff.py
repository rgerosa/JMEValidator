import FWCore.ParameterSet.Config as cms
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

#PFCHS
pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

# make GEN jets
AK4GenJets = ak4GenJets.clone( rParam = 0.4 )
AK4GenJets.src = cms.InputTag('packedGenParticles');

AK8GenJets = ak4GenJets.clone( rParam = 0.8 )
AK8GenJets.src = cms.InputTag('packedGenParticles');

AK4PFchsJets = ak4PFJets.clone();
AK4PFchsJets.src = cms.InputTag('pfCHS','','JRA');

AK8PFchsJets = ak4PFJets.clone( rParam = 0.8 );
AK8PFchsJets.src = cms.InputTag('pfCHS','','JRA');

AK4PFJetsPuppi = ak4PFJetsPuppi.clone( )
AK4PFJetsPuppi.src =  cms.InputTag('puppi','','JRA') #PFJetParameters

AK8PFJetsPuppi = ak4PFJetsPuppi.clone( rParam = 0.8 )
AK8PFJetsPuppi.src =  cms.InputTag('puppi','','JRA') #PFJetParameters
