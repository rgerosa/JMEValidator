import FWCore.ParameterSet.Config as cms
import sys

from PhysicsTools.PatAlgos.tools.tauTools import *

def createProcess(isMC, ## isMC flag
                  processName,
                  globalTag, ## global tag
                  muonTypeID, muonIsoCone, ## muons
                  electronTypeID, ## electrons
                  tauTypeID, ## taus
                  applyZSelections, 
                  jetPtCut,
                  useJECFromLocalDB
                  ):

    process = cms.Process(processName)

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')

    process.GlobalTag.globaltag = globalTag

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Input
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ###########################
    ## Electrons and photons ##
    ###########################
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, switchOnVIDPhotonIdProducer, DataFormat, setupAllVIDIdsInModule, setupVIDElectronSelection, setupVIDPhotonSelection

    if not hasattr(process,"egmGsfElectronIDs"):
        electronIdModules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                             'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

        switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

        for idMod in electronIdModules:
            setupAllVIDIdsInModule(process, idMod, setupVIDElectronSelection)

    if not hasattr(process,"VersionedPhotonIdProducer"):
        switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

        photonIdModules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

        for idMod in photonIdModules:
            setupAllVIDIdsInModule(process, idMod, setupVIDPhotonSelection)
     ## create the Path
    process.jmfw_analyzers = cms.Sequence()
    process.p = cms.Path(process.jmfw_analyzers)
 
    from JMEAnalysis.JMEValidator.runMVAMET_cff import runMVAMET
 
    runMVAMET( process, 
                  processName,
                  isMC,
                  srcMuons = "slimmedMuons", 
                  muonTypeID = "Tight", 
                  iso_map_muons = [], 
                  typeIsoMuons = "dBeta",
                  srcElectrons = "slimmedElectrons", 
                  electronTypeID = "Tight", 
                  electronID_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight',
                  iso_map_electrons = [], 
                  typeIsoElectrons = "rhoCorr",
                  srcTaus = "slimmedTaus", 
                  tauTypeID = tauTypeID, 
                  doTauCleaning = True,
                  jetCollectionPF    = "slimmedJets",
                  dRCleaning = 0.3, 
                  jetPtCut = jetPtCut, 
                  jetEtaCut = 5.,
                  etaCutForMetDiagnostic = 1000,
                  cleanGenJets = True,
                  applyZSelections = applyZSelections
                  )
 
 
    # Run
    if isMC:
        process.run = cms.EDAnalyzer('RunAnalyzer')
        process.jmfw_analyzers += process.run
    


    setattr(process, "MVAMETAnalyzer", 
            cms.EDAnalyzer('METAnalyzer',
                           isMC      = cms.bool(isMC),
                           srcJet    = cms.InputTag("slimmedJetsCleaned"),
                           srcJetPF  = cms.InputTag("slimmedJetsCleaned"),
                           srcVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                           srcZboson = cms.InputTag("MVAMET","ZtagBoson"),
                           srcLeptons = cms.InputTag("LeptonMerge"),
                           srcGenMet = cms.InputTag("slimmedMETs","","PAT"),
                           srcGenJets          = cms.InputTag("slimmedGenJets","","PAT"),
                           srcGenJetsCleaned   = cms.InputTag("slimmedGenJets", "", "PAT"),
                           srcGenParticles     = cms.InputTag("prunedGenParticles","","PAT"),
                           srcGenEventInfo     = cms.InputTag("generator"),

                           srcMETs             = cms.VInputTag(
                                                 cms.InputTag("slimmedMETs"),
                                                 cms.InputTag("patpfMET"),
                                                 cms.InputTag("patpfTrackMET"),
                                                 cms.InputTag("patpfNoPUMET"),
                                                 cms.InputTag("patpfPUCorrectedMET"),
                                                 cms.InputTag("patpfPUMET"),
                                                 cms.InputTag("slimmedMETsPuppi"),
                                                 cms.InputTag("MVAMET", "MVAMET") 
                                                 ),
                           srcRecoils          = cms.VInputTag(
                                                               cms.InputTag("MVAMET", "recoilslimmedMETs"),
                                                               cms.InputTag("MVAMET", "recoilpatpfMET"),
                                                               cms.InputTag("MVAMET", "recoilpatpfTrackMET"),
                                                               cms.InputTag("MVAMET", "recoilpatpfNoPUMET"),
                                                               cms.InputTag("MVAMET", "recoilpatpfPUCorrectedMET"),
                                                               cms.InputTag("MVAMET", "recoilpatpfPUMET"),
                                                               cms.InputTag("MVAMET", "recoilslimmedMETsPuppi"),
                                                               cms.InputTag("MVAMET", "recoilMVAMET")
                           ),
                           dRgenMatching = cms.double(0.3)
                           )
            )

    process.jmfw_analyzers += getattr(process,"MVAMETAnalyzer")

    return process
