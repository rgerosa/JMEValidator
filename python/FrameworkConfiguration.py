def createProcess(isMC, globalTag):

    import FWCore.ParameterSet.Config as cms

    # Common parameters used in all modules
    JetMETAnalyserCommonParameters = cms.PSet(
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

    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
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
    #  }

    jetsCollections = {
            'AK1': {
                'algo': 'ak1',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK1PFPUPPI', 'AK1PFchs', 'AK1PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK2': {
                'algo': 'ak2',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK2PFPUPPI', 'AK2PFchs', 'AK2PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK3': {
                'algo': 'ak3',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK3PFPUPPI', 'AK3PFchs', 'AK3PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK4': {
                'algo': 'ak4',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK4PFPUPPI', 'AK4PFchs', 'AK4PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK5': {
                'algo': 'ak5',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK5PFPUPPI', 'AK5PFchs', 'AK5PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK6': {
                'algo': 'ak6',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK6PFPUPPI', 'AK6PFchs', 'AK6PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK7': {
                'algo': 'ak7',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK7PFPUPPI', 'AK7PFchs', 'AK7PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK8': {
                'algo': 'ak8',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK8PFPUPPI', 'AK8PFchs', 'AK8PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK9': {
                'algo': 'ak9',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK9PFPUPPI', 'AK9PFchs', 'AK9PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },

            'AK10': {
                'algo': 'ak10',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK10PFPUPPI', 'AK10PFchs', 'AK10PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
                },
            }

    from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

    for name, params in jetsCollections.items():
        for index, pu_method in enumerate(params['pu_methods']):
            # Add the jet collection
            jetToolbox(process, params['algo'], 'dummy', 'out', PUMethod = pu_method, JETCorrPayload = params['jec_payloads'][index], JETCorrLevels = params['jec_levels'])


    # Configure the analyzers

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

    process.jmfw_analyzers = cms.Sequence()

    for name, params in jetsCollections.items():
        for index, pu_method in enumerate(params['pu_methods']):

            algo = params['algo'].upper()
            jetCollection = 'selectedPatJets%sPF%s' % (algo, pu_method)

            print('Adding analyzer for jets collection \'%s\'' % jetCollection)

            analyzer = cms.EDAnalyzer('JetMETAnalyzer',
                    JetMETAnalyserCommonParameters,
                    JetCorLabel   = cms.string(params['jec_payloads'][index]),
                    JetCorLevels  = cms.vstring(params['jec_levels']),
                    srcJet        = cms.InputTag(jetCollection),
                    srcRho        = cms.InputTag('fixedGridRhoFastjetAll'),
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

    process.p = cms.Path( process.puppiReader + process.event + process.jmfw_analyzers )


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Output and Log
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
    process.options.allowUnscheduled = cms.untracked.bool(True)

    process.out.outputCommands  = cms.untracked.vstring('drop *')

    # schedule definition
    process.outpath  = cms.EndPath(process.out)

    return process

    #!
    #! THAT'S ALL! CAN YOU BELIEVE IT? :-D
    #!
