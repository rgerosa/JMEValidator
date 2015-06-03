def createProcess(isMC, globalTag):

    import FWCore.ParameterSet.Config as cms

    # Common parameters used in all modules
    JetAnalyserCommonParameters = cms.PSet(
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

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
    process.load("Configuration.EventContent.EventContent_cff")
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')

    process.GlobalTag.globaltag = globalTag


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Input
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
    process.source = cms.Source("PoolSource")

    # Services
    process.load('FWCore.MessageLogger.MessageLogger_cfi')
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
    process.load('CommonTools.UtilAlgos.TFileService_cfi')
    process.TFileService.fileName = cms.string('output_mc.root') if isMC else cms.string('output_data.root')

    # Create all needed jets collections

    # jetsCollections is a dictionnary containing all the informations needed for creating a new jet collection. The format used is :
    #  "name": {
    #      "algo": string ; the jet algorithm to use
    #      "pu_methods" : array of strings ; which PU method to use
    #      "jec_payloads" : array of strings ; which JEC payload to use for making the JEC. The size must match the size of pu_methods
    #      "jec_levels" : array of strings ; which JEC levels to apply
    #      "pu_jet_id": run the pu jet id or not. Very time consuming
    #  }

    # Jet corrections
    process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

    jetsCollections = {
            'AK1': {
                'algo': 'ak1',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK1PFPUPPI', 'AK1PFchs', 'AK1PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK2': {
                'algo': 'ak2',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK2PFPUPPI', 'AK2PFchs', 'AK2PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK3': {
                'algo': 'ak3',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK3PFPUPPI', 'AK3PFchs', 'AK3PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK4': {
                'algo': 'ak4',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK4PFPUPPI', 'AK4PFchs', 'AK4PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': True,
                'qg_tagger': True,
                },

            'AK5': {
                'algo': 'ak5',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK5PFPUPPI', 'AK5PFchs', 'AK5PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK6': {
                'algo': 'ak6',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK6PFPUPPI', 'AK6PFchs', 'AK6PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK7': {
                'algo': 'ak7',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK7PFPUPPI', 'AK7PFchs', 'AK7PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK8': {
                'algo': 'ak8',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK8PFPUPPI', 'AK8PFchs', 'AK8PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK9': {
                'algo': 'ak9',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK9PFPUPPI', 'AK9PFchs', 'AK9PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },

            'AK10': {
                'algo': 'ak10',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK10PFPUPPI', 'AK10PFchs', 'AK10PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': False,
                },
            }

    from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
    from PhysicsTools.PatAlgos.tools.helpers import loadWithPostfix, applyPostfix

    process.load('RecoJets.JetProducers.QGTagger_cfi')

    for name, params in jetsCollections.items():
        for index, pu_method in enumerate(params['pu_methods']):
            # Add the jet collection
            jetToolbox(process, params['algo'], 'dummy', 'out', PUMethod = pu_method, JETCorrPayload = params['jec_payloads'][index], JETCorrLevels = params['jec_levels'], addPUJetID = False)

            algo = params['algo'].upper()
            jetCollection = '%sPFJets%s' % (params['algo'], pu_method)
            postfix = '%sPF%s' % (algo, pu_method)

            # FIXME: PU Jet id is not working with puppi jets
            if params['pu_jet_id'] and pu_method != 'Puppi':

                # PU jet Id
                loadWithPostfix(process, 'RecoJets.JetProducers.pileupjetidproducer_cfi', postfix)
                applyPostfix(process, "pileupJetIdEvaluator", postfix).jets = cms.InputTag(jetCollection)
                applyPostfix(process, "pileupJetIdCalculator", postfix).jets = cms.InputTag(jetCollection)
                applyPostfix(process, "pileupJetIdEvaluator", postfix).rho = cms.InputTag("fixedGridRhoFastjetAll")
                applyPostfix(process, "pileupJetIdEvaluator", postfix).vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
                applyPostfix(process, "pileupJetIdCalculator", postfix).rho = cms.InputTag("fixedGridRhoFastjetAll")
                applyPostfix(process, "pileupJetIdCalculator", postfix).vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")

                # Add informations as userdata: easily accessible
                applyPostfix(process, 'patJets', postfix).userData.userFloats.src += ['pileupJetIdEvaluator%s:fullDiscriminant' % postfix]
                applyPostfix(process, 'patJets', postfix).userData.userInts.src += ['pileupJetIdEvaluator%s:cutbasedId' % postfix, 'pileupJetIdEvaluator%s:fullId' % postfix]

            # Quark / gluon discriminator
            # FIXME: Puppi needs some love
            if 'qg_tagger' in params and params['qg_tagger'] and pu_method != 'Puppi':

                taggerPayload = 'QGL_%sPF%s' % (algo, pu_method.lower())

                setattr(process, 'QGTagger%s' % postfix, process.QGTagger.clone(
                        srcJets = cms.InputTag(jetCollection),
                        jetsLabel = cms.string(taggerPayload)
                    ))

                applyPostfix(process, "patJets", postfix).userData.userFloats.src += ['QGTagger%s:qgLikelihood' % postfix]

    # Compute PF-weighted and PUPPI-weighted isolation
    # See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonIsolationForRun2 for details

    muon_src, cone_size = 'slimmedMuons', 0.4

    ## Create PF candidate collections from packed PF candidates
    ### Using CHS
    process.pfPileUpIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV <= 1"))
    process.pfNoPileUpIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV > 1"))

    process.pfAllPhotons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("pfNoPileUpIso"), cut = cms.string("pdgId == 22"))
    process.pfAllNeutralHadrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("pfNoPileUpIso"), cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
    process.pfAllChargedParticles = cms.EDFilter("CandPtrSelector", src = cms.InputTag("pfNoPileUpIso"), cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212 || pdgId == 11 || pdgId == -11 || pdgId == 13 || pdgId == -13"))
    process.pfAllChargedHadrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("pfNoPileUpIso"), cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))

    process.pfPileUpAllChargedParticles = process.pfAllChargedParticles.clone(
            src = 'pfPileUpIso'
            )

    ### Using puppi
    process.puppiR05 = process.puppi.clone()
    process.puppiR05.algos[0].puppiAlgos[0].cone = 0.5

    process.pfAllPhotonsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 22"))
    process.pfAllNeutralHadronsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
    process.pfAllChargedHadronsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))

    ### Using puppi, but without muons
    ### FIXME: Reference code [1] excludes particles no coming from PV. It leads to an inconsistency between the two puppi collections (one is done on all pf candidates, the other only on
    ### candidates coming from PV)
    ### [1] https://github.com/cms-jet/JMEValidator/blob/a61ebd818c82dc9eab9d47b616ea85136488e77c/python/runMuonIsolation_cff.py#L16
    process.packedPFCandidatesNoMuon = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV > 1 && abs(pdgId) != 13"))
    process.puppiR05NoMu = process.puppiR05.clone(
            candName = 'packedPFCandidatesNoMuon'
            )

    process.pfAllPhotonsPuppiNoMuon = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05NoMu"), cut = cms.string("pdgId == 22"))
    process.pfAllNeutralHadronsPuppiNoMuon = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05NoMu"), cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
    process.pfAllChargedHadronsPuppiNoMuon = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05NoMu"), cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))


    ## Create pf weighted collections
    process.load('CommonTools.ParticleFlow.deltaBetaWeights_cff')

    ## Create isoDeposits with the newly created pf particles collections.
    from JMEAnalysis.JMEValidator.MuonIsolationTools import load_muonPFiso_sequence

    ### PF weighted isolation
    load_muonPFiso_sequence(process, 'MuonPFIsoSequencePFWGT', algo = 'R04PFWGT',
            src = muon_src,
            src_neutral_hadron = 'pfWeightedNeutralHadrons',
            src_photon         = 'pfWeightedPhotons',
            coneR = cone_size
            )

    ### PUPPI weighted isolation
    load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPI', algo = 'R04PUPPI',
            src = muon_src,
            src_charged_hadron = 'pfAllChargedHadronsPuppi',
            src_neutral_hadron = 'pfAllNeutralHadronsPuppi',
            src_photon         = 'pfAllPhotonsPuppi',
            coneR = cone_size
            )

    ### PUPPI weighted isolation without muons
    load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPINoMu', algo = 'R04PUPPINoMu',
            src = muon_src,
            src_charged_hadron = 'pfAllChargedHadronsPuppiNoMuon',
            src_neutral_hadron = 'pfAllNeutralHadronsPuppiNoMuon',
            src_photon         = 'pfAllPhotonsPuppiNoMuon',
            coneR = cone_size
            )

    # Compute electrons and photons IDs
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, switchOnVIDPhotonIdProducer, DataFormat, setupAllVIDIdsInModule, setupVIDElectronSelection, setupVIDPhotonSelection
    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

    electronIdModules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

    photonIdModules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

    for idMod in electronIdModules:
        setupAllVIDIdsInModule(process, idMod, setupVIDElectronSelection)

    for idMod in photonIdModules:
        setupAllVIDIdsInModule(process, idMod, setupVIDPhotonSelection)

    # Create METs from CHS and PUPPI
    from PhysicsTools.PatAlgos.tools.metTools import addMETCollection

    ## Gen MET

    ### Copied from https://github.com/cms-sw/cmssw/blob/2b75137e278b50fc967f95929388d430ef64710b/RecoMET/Configuration/python/GenMETParticles_cff.py#L37
    process.genParticlesForMETAllVisible = cms.EDProducer(
            "InputGenJetsParticleSelector",
            src = cms.InputTag("prunedGenParticles"),
            partonicFinalState = cms.bool(False),
            excludeResonances = cms.bool(False),
            excludeFromResonancePids = cms.vuint32(),
            tausAsJets = cms.bool(False),

            ignoreParticleIDs = cms.vuint32(
                1000022,
                1000012, 1000014, 1000016,
                2000012, 2000014, 2000016,
                1000039, 5100039,
                4000012, 4000014, 4000016,
                9900012, 9900014, 9900016,
                39, 12, 14, 16
                )
            )
    process.load('RecoMET.METProducers.genMetTrue_cfi')

    ## Raw PF METs
    process.load('RecoMET.METProducers.PFMET_cfi')
    process.pfMetCHS = process.pfMet.clone()
    process.pfMetCHS.src = cms.InputTag("chs")
    process.pfMetCHS.alias = cms.string('pfMetCHS')
    addMETCollection(process, labelName='patPFMetCHS', metSource='pfMetCHS') # RAW CHS MET
    process.patPFMetCHS.addGenMET = False

    process.pfMetPuppi = process.pfMet.clone()
    process.pfMetPuppi.src = cms.InputTag("puppi")
    process.pfMetPuppi.alias = cms.string('pfMetPuppi')
    addMETCollection(process, labelName='patPFMetPuppi', metSource='pfMetPuppi') # RAW puppi MET
    process.patPFMetPuppi.addGenMET = False

    ## Type 1 corrections
    from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
    from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

    ### PUPPI
    if not hasattr(process, 'ak4PFJetsPuppi'):
        print("WARNING: No AK4 puppi jets produced. Type 1 corrections for puppi MET are not available.")
    else:
        process.corrPfMetType1Puppi = corrPfMetType1.clone(
            src = 'ak4PFJetsPuppi',
            jetCorrLabel = 'ak4PFCHSL2L3Corrector', #FIXME: Use PUPPI corrections when available?
        )
        del process.corrPfMetType1Puppi.offsetCorrLabel # no L1 for PUPPI jets
        process.pfMetT1Puppi = pfMetT1.clone(
            src = 'pfMetPuppi',
            srcCorrections = [ cms.InputTag("corrPfMetType1Puppi","type1") ]
        )

        addMETCollection(process, labelName='patMETPuppi', metSource='pfMetT1Puppi') # T1 puppi MET
        process.patMETPuppi.addGenMET = False

    ### CHS
    if not hasattr(process, 'ak4PFJetsCHS'):
        print("WARNING: No AK4 CHS jets produced. Type 1 corrections for CHS MET are not available.")
    else:
        process.corrPfMetType1CHS = corrPfMetType1.clone(
            src = 'ak4PFJetsCHS',
            jetCorrLabel = 'ak4PFL1FastL2L3Corrector',
            offsetCorrLabel = 'ak4PFCHSL1FastjetCorrector'
        )
        process.pfMetT1CHS = pfMetT1.clone(
            src = 'pfMetCHS',
            srcCorrections = [ cms.InputTag("corrPfMetType1CHS","type1") ]
        )

        addMETCollection(process, labelName='patMETCHS', metSource='pfMetT1CHS') # T1 puppi MET
        process.patMETCHS.addGenMET = False


    ## Slimmed METs
    from PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi import slimmedMETs
    #### CaloMET is not available in MiniAOD
    del slimmedMETs.caloMET

    ### PUPPI
    process.slimmedMETsPuppi = slimmedMETs.clone()

    if hasattr(process, "patMETPuppi"):
        # Create MET from Type 1 PF collection
        process.patMETPuppi.addGenMET = True
        process.slimmedMETsPuppi.src = cms.InputTag("patMETPuppi")
        process.slimmedMETsPuppi.rawUncertainties = cms.InputTag("patPFMetPuppi") # only central value
    else:
        # Create MET from RAW PF collection
        process.patPFMetPuppi.addGenMET = True
        process.slimmedMETsPuppi.src = cms.InputTag("patPFMetPuppi")
        del process.slimmedMETsPuppi.rawUncertainties # not available

    del process.slimmedMETsPuppi.type1Uncertainties # not available
    del process.slimmedMETsPuppi.type1p2Uncertainties # not available

    ### CHS
    process.slimmedMETsCHS = slimmedMETs.clone()
    if hasattr(process, "patMETCHS"):
        # Create MET from Type 1 PF collection
        process.patMETCHS.addGenMET = True
        process.slimmedMETsCHS.src = cms.InputTag("patMETCHS")
        process.slimmedMETsCHS.rawUncertainties = cms.InputTag("patPFMetCHS") # only central value
    else:
        # Create MET from RAW PF collection
        process.patPFMetCHS.addGenMET = True
        process.slimmedMETsCHS.src = cms.InputTag("patPFMetCHS")
        del process.slimmedMETsCHS.rawUncertainties # not available

    del process.slimmedMETsCHS.type1Uncertainties # not available
    del process.slimmedMETsCHS.type1p2Uncertainties # not available

    # Configure the analyzers

    process.jmfw_analyzers = cms.Sequence()

    # Run
    process.run = cms.EDAnalyzer('RunAnalyzer')
    process.jmfw_analyzers += process.run

    # Event
    from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets
    process.kt6PFJetsRhos = kt6PFJets.clone(
            src = cms.InputTag('packedPFCandidates'),
            doFastJetNonUniform = cms.bool(True),
            puCenters = cms.vdouble(5,-4,-3,-2,-1,0,1,2,3,4,5),
            puWidth = cms.double(.8),
            nExclude = cms.uint32(2))

    process.event = cms.EDAnalyzer('EventAnalyzer',
            rho        = cms.InputTag('fixedGridRhoFastjetAll'),
            rhos       = cms.InputTag('kt6PFJetsRhos', 'rhos'),
            vertices   = cms.InputTag('offlineSlimmedPrimaryVertices')
            )

    process.jmfw_analyzers += process.event

    # HLTs
    process.hlt = cms.EDAnalyzer('HLTAnalyzer',
            src = cms.InputTag('TriggerResults', '', 'HLT'),
            prescales = cms.InputTag('patTrigger'),
            objects = cms.InputTag("selectedPatTrigger"),
            )

    process.jmfw_analyzers += process.hlt

    # Vertices
    process.vertex = cms.EDAnalyzer('VertexAnalyzer',
            src = cms.InputTag('offlineSlimmedPrimaryVertices')
            )

    process.jmfw_analyzers += process.vertex

    # Muons
    process.muons = cms.EDAnalyzer('MuonAnalyzer',
            src = cms.InputTag('slimmedMuons'),
            vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
            rho = cms.InputTag('fixedGridRhoFastjetAll'),
            isoValue_NH_pfWeighted_R04 = cms.InputTag('muPFIsoValueNHR04PFWGT'),
            isoValue_Ph_pfWeighted_R04 = cms.InputTag('muPFIsoValuePhR04PFWGT'),

            isoValue_CH_puppiWeighted_R04 = cms.InputTag('muPFIsoValueCHR04PUPPI'),
            isoValue_NH_puppiWeighted_R04 = cms.InputTag('muPFIsoValueNHR04PUPPI'),
            isoValue_Ph_puppiWeighted_R04 = cms.InputTag('muPFIsoValuePhR04PUPPI'),

            isoValue_CH_puppiNoMuonWeighted_R04 = cms.InputTag('muPFIsoValueCHR04PUPPINoMu'),
            isoValue_NH_puppiNoMuonWeighted_R04 = cms.InputTag('muPFIsoValueNHR04PUPPINoMu'),
            isoValue_Ph_puppiNoMuonWeighted_R04 = cms.InputTag('muPFIsoValuePhR04PUPPINoMu'),
            )

    process.jmfw_analyzers += process.muons

    # Electrons
    process.electrons = cms.EDAnalyzer('ElectronAnalyzer',
            src = cms.InputTag('slimmedElectrons'),
            vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
            conversions = cms.InputTag('reducedEgamma:reducedConversions'),
            beamspot = cms.InputTag('offlineBeamSpot'),
            rho = cms.InputTag('fixedGridRhoFastjetAll'),
            ids = cms.VInputTag('egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto', 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose', 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium', 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight', 'egmGsfElectronIDs:heepElectronID-HEEPV51')
            )

    process.jmfw_analyzers += process.electrons

    # Photons
    process.photons = cms.EDAnalyzer('PhotonAnalyzer',
            src = cms.InputTag('slimmedPhotons'),
            electrons = cms.InputTag('slimmedElectrons'),
            vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
            conversions = cms.InputTag('reducedEgamma:reducedConversions'),
            beamspot = cms.InputTag('offlineBeamSpot'),
            rho = cms.InputTag('fixedGridRhoFastjetAll'),
            phoChargedHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
            phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
            phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
            effAreaChHadFile = cms.FileInPath("EgammaAnalysis/PhotonTools/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
            effAreaNeuHadFile = cms.FileInPath("EgammaAnalysis/PhotonTools/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
            effAreaPhoFile = cms.FileInPath("EgammaAnalysis/PhotonTools/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt"),
            ids = cms.VInputTag('egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose', 'egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium', 'egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight')
            )

    process.jmfw_analyzers += process.photons

    # Jets
    for name, params in jetsCollections.items():
        for index, pu_method in enumerate(params['pu_methods']):

            algo = params['algo'].upper()
            jetCollection = 'selectedPatJets%sPF%s' % (algo, pu_method)

            print('Adding analyzer for jets collection \'%s\'' % jetCollection)

            analyzer = cms.EDAnalyzer('JMEJetAnalyzer',
                    JetAnalyserCommonParameters,
                    JetCorLabel   = cms.string(params['jec_payloads'][index]),
                    JetCorLevels  = cms.vstring(params['jec_levels']),
                    srcJet        = cms.InputTag(jetCollection),
                    srcRho        = cms.InputTag('fixedGridRhoFastjetAll'),
                    srcVtx        = cms.InputTag('offlineSlimmedPrimaryVertices'),
                    srcMuons      = cms.InputTag('selectedPatMuons')
                    )

            setattr(process, params['jec_payloads'][index], analyzer)

            process.jmfw_analyzers += analyzer

    # MET
    process.met_chs = cms.EDAnalyzer('JMEMETAnalyzer',
            src = cms.InputTag('slimmedMETsCHS', '', 'JRA'),
            caloMET = cms.InputTag('slimmedMETs', '', 'PAT')
            )
    process.jmfw_analyzers += process.met_chs

    process.met_puppi = cms.EDAnalyzer('JMEMETAnalyzer',
            src = cms.InputTag('slimmedMETsPuppi', '', 'JRA'),
            caloMET = cms.InputTag('slimmedMETsPuppi', '', 'PAT')
            )
    process.jmfw_analyzers += process.met_puppi

    # Puppi ; only for the first 1000 events of the job
    ## Turn on diagnostic
    process.puppi.puppiDiagnostics = cms.bool(True)
    process.puppiReader = cms.EDAnalyzer("puppiAnalyzer",
                                            treeName = cms.string("puppiTree"),
                                            maxEvents = cms.int32(1000),
                                            nAlgos = cms.InputTag("puppi", "PuppiNAlgos", "JRA"),
                                            rawAlphas = cms.InputTag("puppi", "PuppiRawAlphas", "JRA"),
                                            alphas = cms.InputTag("puppi", "PuppiAlphas", "JRA"),
                                            alphasMed = cms.InputTag("puppi", "PuppiAlphasMed", "JRA"),
                                            alphasRms = cms.InputTag("puppi", "PuppiAlphasRms", "JRA"),
                                            packedPFCandidates = cms.InputTag("packedPFCandidates")
                                        )

    process.jmfw_analyzers += process.puppiReader

    process.p = cms.Path(process.jmfw_analyzers)


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Output and Log
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
    process.options.allowUnscheduled = cms.untracked.bool(True)

    return process

    #!
    #! THAT'S ALL! CAN YOU BELIEVE IT? :-D
    #!
