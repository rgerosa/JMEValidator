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
    #      "pu_jet_id": run the pu jet id or not. Very time consuming
    #  }

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
            jetToolbox(process, params['algo'], 'dummy', 'out', PUMethod = pu_method, JETCorrPayload = params['jec_payloads'][index], JETCorrLevels = params['jec_levels'])

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

    # Muons
    process.muons = cms.EDAnalyzer('MuonAnalyzer',
            src = cms.InputTag('slimmedMuons'),
            vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
            rho = cms.InputTag('fixedGridRhoFastjetAll'),
            )

    process.jmfw_analyzers += process.muons

    # Electrons
    process.electrons = cms.EDAnalyzer('ElectronAnalyzer',
            src = cms.InputTag('slimmedElectrons'),
            vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
            conversions = cms.InputTag('reducedEgamma:reducedConversions'),
            beamspot = cms.InputTag('offlineBeamSpot'),
            rho = cms.InputTag('fixedGridRhoFastjetAll'),
            )

    process.jmfw_analyzers += process.electrons

    # Photons
    process.photons = cms.EDAnalyzer('PhotonAnalyzer',
            src = cms.InputTag('slimmedPhotons'),
            vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
            conversions = cms.InputTag('reducedEgamma:reducedConversions'),
            beamspot = cms.InputTag('offlineBeamSpot'),
            rho = cms.InputTag('fixedGridRhoFastjetAll'),
            phoChargedHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
            phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
            phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
            effAreaChHadFile = cms.FileInPath("EgammaAnalysis/PhotonTools/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
            effAreaNeuHadFile = cms.FileInPath("EgammaAnalysis/PhotonTools/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
            effAreaPhoFile = cms.FileInPath("EgammaAnalysis/PhotonTools/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt")
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

            setattr(process, 'jmfw_%s' % params['jec_payloads'][index], analyzer)

            process.jmfw_analyzers += analyzer

    # MET
    process.met_chs = cms.EDAnalyzer('JMEMETAnalyzer',
            src = cms.InputTag('slimmedMETs')
            )
    process.jmfw_analyzers += process.met_chs

    process.met_puppi = cms.EDAnalyzer('JMEMETAnalyzer',
            src = cms.InputTag('slimmedMETsPuppi')
            )
    process.jmfw_analyzers += process.met_puppi

    # Puppi ; only for the first 1000 events of the job
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

    process.jmfw_analyzers += process.puppiReader

    process.p = cms.Path(process.jmfw_analyzers)


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Output and Log
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
    process.options.allowUnscheduled = cms.untracked.bool(True)

    process.out.outputCommands  = cms.untracked.vstring('drop *')

    # schedule definition
    process.outpath  = cms.EndPath(process.out)

    return process

    #!
    #! THAT'S ALL! CAN YOU BELIEVE IT? :-D
    #!
