import FWCore.ParameterSet.Config as cms
import sys

from PhysicsTools.PatAlgos.tools.tauTools import *

def createProcess(isMC, ## isMC flag
                  processName,
                  globalTag, ## global tag
                  muonTypeID, muonIsoCone, ## muons
                  electronTypeID, ## electrons
                  tauTypeID, ## taus
                  applyZSelections, ## special settings for PUPPET
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
    
            
    # Create METs from CHS and PUPPI

    ## Raw PF METs
    etaCutForMetDiagnostic = 100
    ## make a selection on eta according to the value defined
    process.pfCandidatesForMET = cms.EDFilter("CandPtrSelector",
                                              src = cms.InputTag("packedPFCandidates"),
                                              cut = cms.string("abs(eta) < %f"%etaCutForMetDiagnostic))


    ## CHS pat MET; raw PF is the slimmedMet in miniAOD + typeI correction
    from RecoMET.METProducers.PFMET_cfi import pfMet
    process.pfCandidatesForMETCHS = cms.EDFilter("CandPtrSelector",
                                                 src = cms.InputTag("chs"),
                                                 cut = cms.string("abs(eta) < %f"%etaCutForMetDiagnostic))


    process.pfMetCHS        = pfMet.clone()
    process.pfMetCHS.src    = cms.InputTag("pfCandidatesForMETCHS") ## packed candidates without fromPV < 1
    process.pfMetCHS.alias  = cms.string('pfMetCHS')

    ## Type 1 corrections
    from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
    from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

    ######
    ### Standard TypeI correction
    ######

     ## create the Path
    process.jmfw_analyzers = cms.Sequence()
    process.p = cms.Path(process.jmfw_analyzers)
 
    ## re run HBHE filter
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
 
    process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
                                                        inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
                                                        reverseDecision = cms.bool(False)
                                                        )
 
    process.jmfw_analyzers += process.HBHENoiseFilterResultProducer
    process.jmfw_analyzers += process.ApplyBaselineHBHENoiseFilter 
 
    from JMEAnalysis.JMEValidator.runMVAPUPPET_cff import runMVAPUPPET
 
    runMVAPUPPET( process, 
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
                  etaCutForMetDiagnostic = etaCutForMetDiagnostic,
                  cleanGenJets = True,
                  applyZSelections = applyZSelections
                  )
 
 
    # Run
    if isMC:
        process.run = cms.EDAnalyzer('RunAnalyzer')
        process.jmfw_analyzers += process.run
    


    setattr(process, "PUPPET", 
            cms.EDAnalyzer('PUPPETAnalyzer',
                           isMC      = cms.bool(isMC),
                           srcJet    = cms.InputTag("slimmedJetsCleaned"),
                           srcJetPF  = cms.InputTag("slimmedJetsCleaned"),
                           srcVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                           srcZboson = cms.InputTag("mvaMET","ZtagBoson"),
                           srcLeptons = cms.InputTag("LeptonMerge"),
                           srcGenMet = cms.InputTag("slimmedMETs","","PAT"),
                           srcGenJets          = cms.InputTag("slimmedGenJets","","PAT"),
                           srcGenJetsCleaned   = cms.InputTag("slimmedGenJets", "", "PAT"),
                           srcGenParticles     = cms.InputTag("prunedGenParticles","","PAT"),
                           srcGenEventInfo     = cms.InputTag("generator"),

                           srcMETs             = cms.VInputTag(cms.InputTag("patpfMET"),
                                                 cms.InputTag("patpfTrackMET"),
                                                 cms.InputTag("patpfNoPUMET"),
                                                 cms.InputTag("patpfPUCorrectedMET"),
                                                 cms.InputTag("patpfPUMET"),
                                                 cms.InputTag("slimmedMETsPuppi"),
                                                 cms.InputTag("slimmedMETs"),
                                                 cms.InputTag("mvaMET", "mvaMET") 
),
                           srcRecoils          = cms.VInputTag(cms.InputTag("mvaMET", "recoilpatpfMET"),
                                                               cms.InputTag("mvaMET", "recoilpatpfTrackMET"),
                                                               cms.InputTag("mvaMET", "recoilpatpfNoPUMET"),
                                                               cms.InputTag("mvaMET", "recoilpatpfPUCorrectedMET"),
                                                               cms.InputTag("mvaMET", "recoilpatpfPUMET"),
                                                               cms.InputTag("mvaMET", "recoilslimmedMETsPuppi"),
                                                               cms.InputTag("mvaMET", "recoilslimmedMETs"),
                                                               cms.InputTag("mvaMET", "recoilmvaMET")
                           ),
                           dRgenMatching = cms.double(0.3),
                           srcMetFiltersBits = cms.InputTag("TriggerResults","","PAT"),
                           srcTriggerBits = cms.InputTag("TriggerResults","","HLT"),
                           srcTriggerPrescales = cms.InputTag('patTrigger')
                           )
            )
    
    if not isMC:
        getattr(process,"PUPPET").srcMetFiltersBits = cms.InputTag("TriggerResults","","RECO")

    process.jmfw_analyzers += getattr(process,"PUPPET")

    return process
