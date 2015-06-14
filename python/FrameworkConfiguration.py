import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.tauTools import *

def createProcess(isMC, globalTag, muonCollection, runPuppiMuonIso, muonIsoCone, electronCollection, tauCollection, dropAnalyzerDumpEDM, runMVAPUPPETAnalysis):

    process = cms.Process("JRA")

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')

    process.GlobalTag.globaltag = globalTag

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
        nJetMax         = cms.uint32( 0),
    )


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Input
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    process.source = cms.Source("PoolSource")

    process.load('CommonTools.UtilAlgos.TFileService_cfi')
    process.TFileService.fileName = cms.string('output_mc.root') if isMC else cms.string('output_data.root')

    # Jet corrections
    process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
    # QG tagger
    process.load('RecoJets.JetProducers.QGTagger_cfi')
    # tool box
    from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
    # manipulate post fix
    from PhysicsTools.PatAlgos.tools.helpers   import loadWithPostfix, applyPostfix


    #######################
    ### JET COLLECTIONS ###
    #######################

    # jetsCollections is a dictionnary containing all the informations needed for creating a new jet collection. The format used is :
    #  "name": {
    #      "algo": string ; the jet algorithm to use
    #      "pu_methods" : array of strings ; which PU method to use
    #      "jec_payloads" : array of strings ; which JEC payload to use for making the JEC. The size must match the size of pu_methods
    #      "jec_levels" : array of strings ; which JEC levels to apply
    #      "pu_jet_id": run the pu jet id or not. Very time consuming
    #  }

    jetsCollections = {
             'AK4': {
                'algo': 'ak4',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK4PFPUPPI', 'AK4PFchs', 'AK4PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': True,
                'qg_tagger': True,
                },
            }
    

    ## loop on the jet collections : generic container just defined by clustering algorithm and cone dimension
    for name, params in jetsCollections.items():
        ## loop on the pileup methos
        for index, pu_method in enumerate(params['pu_methods']):

            algo          = params['algo'].upper()
            jetCollection = '%sPFJets%s' % (params['algo'], pu_method)
            postfix       = '%sPF%s' % (algo, pu_method)

            # Add the jet collection via the jetToolBox
            jetToolbox(process, params['algo'], postfix+"Sequence", 'out', PUMethod = pu_method, JETCorrPayload = params['jec_payloads'][index], JETCorrLevels = params['jec_levels'], addPUJetID = False)

            # FIXME: PU Jet id is not working with puppi jets or SK jets
            if params['pu_jet_id'] and pu_method != 'Puppi' and pu_method != 'SK':

                # PU jet Id  .. the pileup jet id is run at posteriori since it does not work in the jet tool box for puppi and SK jets
                loadWithPostfix(process, 'RecoJets.JetProducers.pileupjetidproducer_cfi', postfix)
                applyPostfix(process, "pileupJetIdEvaluator", postfix).jets      = cms.InputTag(jetCollection)
                applyPostfix(process, "pileupJetIdCalculator", postfix).jets     = cms.InputTag(jetCollection)
                applyPostfix(process, "pileupJetIdEvaluator", postfix).rho       = cms.InputTag("fixedGridRhoFastjetAll")
                applyPostfix(process, "pileupJetIdEvaluator", postfix).vertexes  = cms.InputTag("offlineSlimmedPrimaryVertices")
                applyPostfix(process, "pileupJetIdCalculator", postfix).rho      = cms.InputTag("fixedGridRhoFastjetAll")
                applyPostfix(process, "pileupJetIdCalculator", postfix).vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")

                # Add informations as userdata: easily accessible
                applyPostfix(process, 'patJets', postfix).userData.userFloats.src += ['pileupJetIdEvaluator%s:fullDiscriminant' % postfix]
                applyPostfix(process, 'patJets', postfix).userData.userInts.src   += ['pileupJetIdEvaluator%s:cutbasedId' % postfix, 'pileupJetIdEvaluator%s:fullId' % postfix]

            # Quark / gluon discriminator
            # FIXME: Puppi needs some love
            # FIXME: So does SK
            if 'qg_tagger' in params and params['qg_tagger'] and pu_method != 'Puppi' and pu_method != 'SK':

                taggerPayload = 'QGL_%sPF%s' % (algo, pu_method.lower())

                setattr(process, 'QGTagger%s' % postfix, process.QGTagger.clone(
                        srcJets = cms.InputTag(jetCollection),
                        jetsLabel = cms.string(taggerPayload)))

                applyPostfix(process, "patJets", postfix).userData.userFloats.src += ['QGTagger%s:qgLikelihood' % postfix]

    ######################
    ### MUONS ISOLATION ##
    ######################

    # Compute PF-weighted and PUPPI-weighted isolation
    # See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonIsolationForRun2 for details

    ## Create PF candidate collections from packed PF candidates Using CHS
    process.pfPileUpIso   = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV <= 1")) ## cut away PV particles
    process.pfNoPileUpIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV > 1"))  ## take particles from PV only

    process.pfAllPhotons        = cms.EDFilter("CandPtrSelector", src = cms.InputTag("pfNoPileUpIso"), cut = cms.string("pdgId == 22"))
    process.pfAllNeutralHadrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("pfNoPileUpIso"), 
                                               cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
    
    process.pfAllChargedParticles = cms.EDFilter("CandPtrSelector",src = cms.InputTag("pfNoPileUpIso"),
                                                 cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212 || pdgId == 11 || pdgId == -11 || pdgId == 13 || pdgId == -13"))
    process.pfAllChargedHadrons   = cms.EDFilter("CandPtrSelector",src = cms.InputTag("pfNoPileUpIso"), 
                                                 cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))
    process.pfPileUpAllChargedParticles = process.pfAllChargedParticles.clone( src = 'pfPileUpIso')
    
    ### Using puppi with R05 for muons isolation, is not the once of jets but the puppi cone in the barrel region
    if runPuppiMuonIso :

        process.puppiR05 = process.puppi.clone()
        process.puppiR05.algos[0].puppiAlgos[0].cone = 0.5

        process.pfAllPhotonsPuppi        = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 22"))
        process.pfAllNeutralHadronsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
        process.pfAllChargedHadronsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))

        ### Using puppi, but without muons
        ### FIXME: Reference code [1] excludes particles no coming from PV. It leads to an inconsistency between the two puppi collections (one is done on all pf candidates, the other only on candidates coming from PV)
        ### [1] https://github.com/cms-jet/JMEValidator/blob/a61ebd818c82dc9eab9d47b616ea85136488e77c/python/runMuonIsolation_cff.py#L16
        process.packedPFCandidatesNoMuon = cms.EDFilter("CandPtrSelector", 
                                                        src = cms.InputTag("packedPFCandidates"), 
                                                        cut = cms.string("fromPV > 1 && abs(pdgId) != 13"))
        process.puppiR05NoMu = process.puppiR05.clone(
            candName = 'packedPFCandidatesNoMuon'
            )

        process.pfAllPhotonsPuppiNoMuon        = cms.EDFilter("CandPtrSelector", 
                                                              src = cms.InputTag("puppiR05NoMu"), 
                                                              cut = cms.string("pdgId == 22"))
        process.pfAllNeutralHadronsPuppiNoMuon = cms.EDFilter("CandPtrSelector", 
                                                              src = cms.InputTag("puppiR05NoMu"), 
                                                              cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
        process.pfAllChargedHadronsPuppiNoMuon = cms.EDFilter("CandPtrSelector", 
                                                              src = cms.InputTag("puppiR05NoMu"), 
                                                              cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))


    ## Create pf weighted collections
    process.load('CommonTools.ParticleFlow.deltaBetaWeights_cff')

    ## Create isoDeposits with the newly created pf particles collections.
    from JMEAnalysis.JMEValidator.MuonIsolationTools import load_muonPFiso_sequence

    ### PF weighted isolation
    load_muonPFiso_sequence(process, 
                            'MuonPFIsoSequencePFWGT', 
                            algo = 'R04PFWGT',
                            src = muonCollection,
                            src_neutral_hadron = 'pfWeightedNeutralHadrons',
                            src_photon         = 'pfWeightedPhotons',
                            coneR = muonIsoCone
                            )
    if runPuppiMuonIso :

        ### PUPPI weighted isolation
        load_muonPFiso_sequence(process, 
                                'MuonPFIsoSequencePUPPI', 
                                algo = 'R04PUPPI',
                                src =  muonCollection,
                                src_charged_hadron = 'pfAllChargedHadronsPuppi',
                                src_neutral_hadron = 'pfAllNeutralHadronsPuppi',
                                src_photon         = 'pfAllPhotonsPuppi',
                                coneR = muonIsoCone
                                )
    
        ### PUPPI weighted isolation without muons
        load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPINoMu', algo = 'R04PUPPINoMu',
                                src =  muonCollection,
                                src_charged_hadron = 'pfAllChargedHadronsPuppiNoMuon',
                                src_neutral_hadron = 'pfAllNeutralHadronsPuppiNoMuon',
                                src_photon         = 'pfAllPhotonsPuppiNoMuon',
                                coneR = muonIsoCone
                                )

    
    

    if runMVAPUPPETAnalysis :

        ##############
        ## apply muon ID : label is used to form the name, tight indicates the id to be applied from POGs, iso map are non standard iso info, typeIsoval is related to what apply: 0 no PU correction, 1 dBeta correction, 2 is rho*Area, 3 is PFWeighted correction, 4 is puppi
        ##############

        from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyMuonID
        
        applyMuonID(process, label = "Tight", src  = muonCollection, type = 'tightID',
                    iso_map_charged_hadron  = '', 
                    iso_map_neutral_hadron  = '', 
                    iso_map_photon          = '',                 
                    rho = 'fixedGridRhoFastjetAll',
                    typeIsoVal = 1
                    )

        applyMuonID(process, label = "TightDBeta", src  = muonCollection, type = 'tightID',
                    iso_map_charged_hadron  = '', 
                    iso_map_neutral_hadron  = 'muPFIsoValueNHR04PFWGT', 
                    iso_map_photon          = 'muPFIsoValuePhR04PFWGT',                 
                    rho = 'fixedGridRhoFastjetAll',
                    typeIsoVal = 3
                    )

        if runPuppiMuonIso :
            
            applyMuonID(process, label = "TightPuppiNoMu", src  = muonCollection, type = 'tightID',
                        iso_map_charged_hadron  = 'muPFIsoValueCHR04PUPPINoMu', 
                        iso_map_neutral_hadron  = 'muPFIsoValueNHR04PUPPINoMu', 
                        iso_map_photon          = 'muPFIsoValuePhR04PUPPINoMu',                 
                        rho = 'fixedGridRhoFastjetAll',
                        typeIsoVal = 4
                        )

    ###########################
    ## Electrons and photons ##
    ###########################

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

    if runMVAPUPPETAnalysis:

        from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyElectronID
        
        #applyElectronID(process, label = "Tight", src  = electronCollection, 
        #                iso_map_charged_hadron  = '', 
        #                iso_map_neutral_hadron  = '', 
        #                iso_map_photon          = '',                 
        #                rho = 'fixedGridRhoFastjetAll',
        #                typeIsoVal = 1,
        #                electron_id_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight'
        #                )
        
        applyElectronID(process, label = "Medium", src  = electronCollection, 
                        iso_map_charged_hadron  = '', 
                        iso_map_neutral_hadron  = '', 
                        iso_map_photon          = '',                 
                        rho = 'fixedGridRhoFastjetAll',
                        typeIsoVal = 1,
                        electron_id_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium'
                        )

    ##################
    #### TAU ID ######
    ##################
    muonIDLabelForCleaning     = ["Tight"]                    
    electronIDLabelForCleaning = ["Medium"]
    tauIDLabelForCleaning      = ["Loose"]

    if runMVAPUPPETAnalysis :
        from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyTauID
        
        for muon in muonIDLabelForCleaning :
            for electron in electronIDLabelForCleaning:

                applyTauID( process, label = "Loose", src = tauCollection,
                            muonCollection     = muonCollection+muon,
                            electronCollection = electronCollection+electron)

        #applyTauID( process, label = "Tight", src = tauCollection,
        #            muonCollection     = "slimmedMuonsTight",
        #            electronCollection = "slimmedElectronsTight")
        

    ##################
    ### clean jets ###
    ##################
        
    if runMVAPUPPETAnalysis :

        from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import cleanJetsFromLeptons

        for name, params in jetsCollections.items():
        ## loop on the pileup methos                                                                                                                                            
            for index, pu_method in enumerate(params['pu_methods']):
                postfix       = '%sPF%s' % (algo, pu_method)

                for muon in muonIDLabelForCleaning :
                    for electron in electronIDLabelForCleaning :
                        for tau in tauIDLabelForCleaning :
                            if tau == "" :
                                cleanJetsFromLeptons(process,"Cleaned"+"Mu"+muon+"Ele"+electron,                                 
                                                     jetCollection      = "selectedPatJets"+postfix,
                                                     muonCollection     = muonCollection+muon,
                                                     electronCollection = electronCollection+electron,
                                                     tauCollection      = tauCollection+tau,
                                                     jetPtCut  = 15.,
                                                     jetEtaCut = 5.,
                                                     dRCleaning = 0.3) 
                            else:
                                cleanJetsFromLeptons(process,"Cleaned"+"Mu"+muon+"Ele"+electron+"Tau"+tau,                                 
                                                     jetCollection      = "selectedPatJets"+postfix,
                                                     muonCollection     = muonCollection+muon,
                                                     electronCollection = electronCollection+electron,
                                                     tauCollection      = "",
                                                     jetPtCut  = 15.,
                                                     jetEtaCut = 5.,
                                                     dRCleaning = 0.3) 

    

    ##################################                    
    ##### run pu puppi for PUPPET ####
    ##################################
    if runMVAPUPPETAnalysis:                            
        process.pupuppi = process.puppi.clone()
        process.pupuppi.invertPuppi = True


    # Create METs from CHS and PUPPI
    from PhysicsTools.PatAlgos.tools.metTools import addMETCollection

    ## Gen MET ###
    ### Copied from https://github.com/cms-sw/cmssw/blob/2b75137e278b50fc967f95929388d430ef64710b/RecoMET/Configuration/python/GenMETParticles_cff.py#L37

    process.load('RecoMET.METProducers.genMetTrue_cfi')

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

    ## CHS pat MET; raw PF is the slimmedMet in miniAOD + typeI correction
    from RecoMET.METProducers.PFMET_cfi import pfMet
    process.pfMetCHS        = pfMet.clone()
    process.pfMetCHS.src    = cms.InputTag("chs") ## packed candidates without fromPV < 1
    process.pfMetCHS.alias  = cms.string('pfMetCHS')
    addMETCollection(process, labelName='patPFMetCHS', metSource='pfMetCHS') # Convert the CHS PFMet in PAT MET
    process.patPFMetCHS.addGenMET = False

    ### puppi raw met
    process.pfMetPuppi       = pfMet.clone()
    process.pfMetPuppi.src   = cms.InputTag("puppi")
    process.pfMetPuppi.alias = cms.string('pfMetPuppi')
    addMETCollection(process, labelName='patPFMetPuppi', metSource='pfMetPuppi') # RAW puppi MET
    process.patPFMetPuppi.addGenMET = False

    ## Type 1 corrections
    from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
    from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

    ### PUPPI TypeI corrected
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
        ## new PAT Met correction
        addMETCollection(process, labelName='patMETPuppi', metSource='pfMetT1Puppi') # T1 puppi MET
        process.patMETPuppi.addGenMET = False

    ### CHS
    if not hasattr(process, 'ak4PFJetsCHS'):
        print("WARNING: No AK4 CHS jets produced. Type 1 corrections for CHS MET are not available.")
    else:
        process.corrPfMetType1CHS = corrPfMetType1.clone(
            src             = 'ak4PFJetsCHS',
            jetCorrLabel    = 'ak4PFL1FastL2L3Corrector',
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

    ### PUPPI : make slimmed METs in order to embed both corrected and not corrected one after TypeI
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


    process.jmfw_analyzers = cms.Sequence()
    process.p = cms.Path(process.jmfw_analyzers)


    if runMVAPUPPETAnalysis :

        ######## add other specific set of particles and MET collections, in this case not TypeI corrected, but we will still use the same workflow    
        ## all charge particles from PUPPI : hadrons + leptons (e,mu,tau) --> trak met
        process.pfAllChargedParticlesPuppi = cms.EDFilter("CandPtrSelector",
                                                          src = cms.InputTag("puppi"),
                                                          cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 22\
12 || pdgId == -2212 || pdgId == 11 || pdgId == -11 || pdgId == 13 || pdgId ==-13 || pdgId == 15 || pdgId == -15"))
        
        process.pfMetPuppiChargedPV       = pfMet.clone()
        process.pfMetPuppiChargedPV.src   = cms.InputTag("pfAllChargedParticlesPuppi") ## packed candidates without fromPV < 1
        process.pfMetPuppiChargedPV.alias = cms.string('pfMetPuppiChargedPV')
        addMETCollection(process, labelName='patPFMetPuppiChargedPV', metSource='pfMetPuppiChargedPV') # Convert the CHS PFMet in PAT MET
        process.patPFMetPuppiChargedPV.addGenMET = True

        process.slimmedMETsPuppiChargedPV = slimmedMETs.clone()
        process.slimmedMETsPuppiChargedPV.src = cms.InputTag("patPFMetPuppiChargedPV")
        process.slimmedMETsPuppiChargedPV.rawUncertainties = cms.InputTag("patPFMetPuppiChargedPV") # only central value
        del process.slimmedMETsPuppiChargedPV.type1Uncertainties # not available                                                                                            
        del process.slimmedMETsPuppiChargedPV.type1p2Uncertainties # not available                                                                                 

        ### charge candidate from PU and build the Met ( we can take the one defined for the isolation)
        process.pfMetChargedPU       = pfMet.clone() 
        process.pfMetChargedPU.src   = cms.InputTag("pfPileUpIso")
        process.pfMetChargedPU.alias = cms.string('pfMetPuppiChargedPU')
        addMETCollection(process, labelName='patPFMetChargedPU', metSource='pfMetChargedPU') 
        process.patPFMetChargedPU.addGenMET = True

        process.slimmedMETsPuppiChargedPU = slimmedMETs.clone()
        process.slimmedMETsPuppiChargedPU.src = cms.InputTag("patPFMetChargedPU")
        process.slimmedMETsPuppiChargedPU.rawUncertainties = cms.InputTag("patPFMetChargedPU") # only central value
        
        del process.slimmedMETsPuppiChargedPU.type1Uncertainties # not available                                                                                           
        del process.slimmedMETsPuppiChargedPU.type1p2Uncertainties # not available                                                                                      
        
        ## neutral particles from PV starting from puppi particles
        process.pfAllNeutralParticlesPuppi  = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppi"), 
                                                           cut = cms.string("pdgId == 22 || pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
        
        
        process.pfMetPuppiNeutralPV       = pfMet.clone() 
        process.pfMetPuppiNeutralPV.src   = cms.InputTag("pfAllNeutralParticlesPuppi")
        process.pfMetPuppiNeutralPV.alias = cms.string('pfMetPuppiNeutralPV')
        addMETCollection(process, labelName='patPFMetPuppiNeutralPV', metSource='pfMetPuppiNeutralPV') 
        process.patPFMetPuppiNeutralPV.addGenMET = True
        
        process.slimmedMETsPuppiNeutralPV = slimmedMETs.clone()
        process.slimmedMETsPuppiNeutralPV.src = cms.InputTag("patPFMetPuppiNeutralPV")
        process.slimmedMETsPuppiNeutralPV.rawUncertainties = cms.InputTag("patPFMetPuppiNeutralPV") # only central value
        
        del process.slimmedMETsPuppiNeutralPV.type1Uncertainties # not available                                                                                                
        del process.slimmedMETsPuppiNeutralPV.type1p2Uncertainties # not available                                                                                          
    

        ## neutral particles from PU starting from inverted puppi particles
        process.pfAllNeutralParticlesPuppiPU  = cms.EDFilter("CandPtrSelector", src = cms.InputTag("pupuppi"), 
                                                             cut = cms.string("pdgId == 22 || pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
    
        
        process.pfMetPuppiNeutralPU       = pfMet.clone() 
        process.pfMetPuppiNeutralPU.src   = cms.InputTag("pfAllNeutralParticlesPuppiPU")
        process.pfMetPuppiNeutralPU.alias = cms.string('pfMetPuppiNeutralPU')
        addMETCollection(process, labelName='patPFMetPuppiNeutralPU', metSource='pfMetPuppiNeutralPU') 
        process.patPFMetPuppiNeutralPU.addGenMET = True

        process.slimmedMETsPuppiNeutralPU = slimmedMETs.clone()
        process.slimmedMETsPuppiNeutralPU.src = cms.InputTag("patPFMetPuppiNeutralPU")
        process.slimmedMETsPuppiNeutralPU.rawUncertainties = cms.InputTag("patPFMetPuppiNeutralPU") # only central value

        del process.slimmedMETsPuppiNeutralPU.type1Uncertainties # not available                                                                                                
        del process.slimmedMETsPuppiNeutralPU.type1p2Uncertainties # not available                                                                                           

        ## merge all the lepton collection into a single one
        for muon in muonIDLabelForCleaning :
            for electron in electronIDLabelForCleaning :
                for tau in tauIDLabelForCleaning :                    

                  ## look for a good ZLL candidate i.e.: two tight muons within Mz, two tight electron or taus in Mz
                  setattr(process,"diMuon"+muon,cms.EDProducer("CandViewCombiner",
                                                             decay       = cms.string(muonCollection+muon+"@+ "+muonCollection+muon+"@-"),
                                                             checkCharge = cms.bool(True),
                                                             cut         = cms.string("mass > 70 && mass < 110 & charge=0"),                                              
                                                             ))

                  process.jmfw_analyzers += getattr(process,"diMuon"+muon);

                  setattr(process,"diElectron"+electron,cms.EDProducer("CandViewCombiner",
                                                             decay       = cms.string(electronCollection+electron+"@+ "+electronCollection+electron+"@-"),
                                                             checkCharge = cms.bool(True),
                                                             cut         = cms.string("mass > 70 && mass < 110 & charge=0"),                                       
                                                             ))

                  process.jmfw_analyzers += getattr(process,"diElectron"+electron);

                  setattr(process,"diTau"+tau,cms.EDProducer("CandViewCombiner",
                                                           decay       = cms.string(tauCollection+tau+"@+ "+tauCollection+tau+"@-"),
                                                           checkCharge = cms.bool(True),
                                                           cut         = cms.string("mass > 70 && mass < 110 & charge=0"),                                               
                                                           ))

                  process.jmfw_analyzers += getattr(process,"diTau"+tau);

                  ## merge all the Z canddates and ask for only one candidate per event
                  setattr(process,"diLeptonMergeMu"+muon+"Ele"+electron+"Tau"+tau, cms.EDProducer("CandViewMerger",
                                                                                                  src = cms.VInputTag("diMuon"+muon,"diElectron"+electron,"diTau"+tau)
                                                                                                  ))

                  process.jmfw_analyzers += getattr(process,"diLeptonMergeMu"+muon+"Ele"+electron+"Tau"+tau);

                  setattr(process,"diLeptonMergeMu"+muon+"Ele"+electron+"Tau"+tau+"Filter",cms.EDFilter("PATCandViewCountFilter",
                                                                                                     minNumber = cms.uint32(1),
                                                                                                     maxNumber = cms.uint32(1),
                                                                                                     src = cms.InputTag("diLeptonMergeMu"+muon+"Ele"+electron+"Tau"+tau)
                                                                                                     ))
                                                                                                                         
                  process.jmfw_analyzers += getattr(process,"diLeptonMergeMu"+muon+"Ele"+electron+"Tau"+tau+"Filter");

                  ## now merge all the previous leptons in one single collection and ask for no more than 2 tight leptons
                  setattr(process,"leptonMergeMu"+muon+"Ele"+electron+"Tau"+tau, cms.EDProducer("CandViewMerger",
                                                                                                  src = cms.VInputTag(muonCollection+muon,electronCollection+electron,tauCollection+tau)
                                                                                                  ))

                  process.jmfw_analyzers += getattr(process,"leptonMergeMu"+muon+"Ele"+electron+"Tau"+tau);

                  setattr(process,"leptonMergeMu"+muon+"Ele"+electron+"Tau"+tau+"Filter",cms.EDFilter("PATCandViewCountFilter",
                                                                                                     minNumber = cms.uint32(2),
                                                                                                     maxNumber = cms.uint32(2),
                                                                                                     src = cms.InputTag("leptonMergeMu"+muon+"Ele"+electron+"Tau"+tau)
                                                                                                     ))
        
                  process.jmfw_analyzers += getattr(process,"leptonMergeMu"+muon+"Ele"+electron+"Tau"+tau+"Filter");

                  ### in this way we are sure to have a Z candidate with no more than 2 leptons in the event (after selections)
                  setattr(process,"mvaPUPPETMu"+muon+"Ele"+electron+"Tau"+tau, cms.EDProducer("mvaPUPPET",
                                   srcMETs      = cms.VInputTag("slimmedMETs","slimmedMETsCHS", "slimmedMETsPuppi","slimmedMETsPuppiChargedPV","slimmedMETsPuppiChargedPU","slimmedMETsPuppiNeutralPV","slimmedMETsPuppiNeutralPU"),
                                   referenceMET   = cms.InputTag("slimmedMETsPuppi"),
                                   srcVertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                   srcJets        = cms.InputTag("selectedPatJetsAK4PFPuppiCleaned"+"Mu"+muon+"Ele"+electron+"Tau"+tau),
                                   inputFileNames = cms.PSet(
                                           #PhiCorrectionWeightFile = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrmet_7X_BX50_Jan2015.root'),
                                           #RecoilCorrectionWeightFile  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7X_BX50_Jan2015.root')
                                   ),
                                   srcLeptons = cms.VInputTag("leptonMergeMu"+muon+"Ele"+electron+"Tau"+tau) ))

                  process.jmfw_analyzers += getattr(process,"mvaPUPPETMu"+muon+"Ele"+electron+"Tau"+tau);


    if dropAnalyzerDumpEDM:        
        return process


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
    
    # Muons : tight muons, DBWeight and puppiNoMu corrected
    if not runMVAPUPPETAnalysis :

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
                                       isoValue_Ph_puppiNoMuonWeighted_R04 = cms.InputTag('muPFIsoValuePhR04PUPPINoMu'))
                                       
        process.jmfw_analyzers += process.muons

        # Electrons
        process.electrons = cms.EDAnalyzer('ElectronAnalyzer',
                                                src = cms.InputTag('slimmedElectrons'),
                                                vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                                conversions = cms.InputTag('reducedEgamma:reducedConversions'),
                                                beamspot = cms.InputTag('offlineBeamSpot'),
                                                rho = cms.InputTag('fixedGridRhoFastjetAll'),
                                                ids = cms.VInputTag()
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

        # Jets : store the cleaned jet collection for our case
        for name, params in jetsCollections.items():
            for index, pu_method in enumerate(params['pu_methods']):

                algo = params['algo'].upper()
                jetCollection = 'selectedPatJets%sPF%s' % (algo, pu_method)

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
        
    else:
        
        for muon in muonIDLabelForCleaning :
            for electron in electronIDLabelForCleaning :
                for tau in tauIDLabelForCleaning :           
                    setattr(process, "puppet"+"Mu"+muon+"Ele"+electron+"Tau"+tau, 
                            cms.EDAnalyzer('PUPPETAnalyzer',
                                           srcJet    = cms.InputTag("selectedPatJetsAK4PFPuppiCleaned"+"Mu"+muon+"Ele"+electron+"Tau"+tau),
                                           srcVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                           srcZboson = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau+":ZtagBoson"),
                                           srcRecoilPFMet = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau,"recoilslimmedMETs"),
                                           srcRecoilPFCHSMet = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau,"recoilslimmedMETsCHS"),
                                           srcRecoilPFPuppiMet = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau,"recoilslimmedMETsPuppi"),
                                           srcRecoilPFPuppiMet_ChargedPV = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau,"recoilslimmedMETsPuppiChargedPV"),
                                           srcRecoilPFPuppiMet_ChargedPU = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau,"recoilslimmedMETsPuppiChargedPU"),
                                           srcRecoilPFPuppiMet_NeutralPV = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau,"recoilslimmedMETsPuppiNeutralPV"),
                                           srcRecoilPFPuppiMet_NeutralPU = cms.InputTag("mvaPUPPET"+"Mu"+muon+"Ele"+electron+"Tau"+tau,"recoilslimmedMETsPuppiNeutralPU")))
                                           
                    process.jmfw_analyzers += getattr(process,"puppetMu"+muon+"Ele"+electron+"Tau"+tau)
                                           

    # Puppi ; only for the first 1000 events of the job
    ## Turn on diagnostic
#    process.puppi.puppiDiagnostics = cms.bool(True)
#    process.puppiReader = cms.EDAnalyzer("puppiAnalyzer",
#                                            treeName = cms.string("puppiTree"),
#                                            maxEvents = cms.int32(1000),
#                                            nAlgos = cms.InputTag("puppi", "PuppiNAlgos", "JRA"),
#                                            rawAlphas = cms.InputTag("puppi", "PuppiRawAlphas", "JRA"),
#                                            alphas = cms.InputTag("puppi", "PuppiAlphas", "JRA"),
#                                            alphasMed = cms.InputTag("puppi", "PuppiAlphasMed", "JRA"),
#                                            alphasRms = cms.InputTag("puppi", "PuppiAlphasRms", "JRA"),
#                                            packedPFCandidates = cms.InputTag("packedPFCandidates")
#                                        )

#    process.jmfw_analyzers += process.puppiReader
    
    return process

    #!
    #! THAT'S ALL! CAN YOU BELIEVE IT? :-D
    #!

