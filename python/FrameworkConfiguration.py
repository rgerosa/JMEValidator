import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.tauTools import *

def get_cone_size(algo):
    import re
    cone_size = re.search('(\d+)$', algo)
    if cone_size is None:
        raise ValueError('Cannot extract cone size from algorithm name')
    
    return int(cone_size.group(1))

def get_jec_payload(algo, pu_method):
    
    # FIXME: Until PUPPI and SK payloads are in the GT, use CHS corrections
    jec_payloads = {
                'Puppi': 'AK%dPFchs',
                'CHS': 'AK%dPFchs',
                'SK': 'AK%dPFchs',
                '': 'AK%dPF',
                }
    
    
    cone_size = get_cone_size(algo)
    
    if not pu_method in jec_payloads:
        print('WARNING: JEC payload not found for method %r. Using default one.' % pu_method)
        return 'None'
    
    return jec_payloads[pu_method] % cone_size

def get_jec_levels(pu_method):
    
    jec_levels = {
        'Puppi': ['L2Relative', 'L3Absolute'],
        'CHS': ['L1FastJet', 'L2Relative', 'L3Absolute'],
        'SK': ['L2Relative', 'L3Absolute'],
        '': ['L1FastJet', 'L2Relative', 'L3Absolute'],
        }
    
    
    if not pu_method in jec_levels:
        print('WARNING: JEC levels not found for method %r. Using default ones.' % pu_method)
        return ['None']
    
    return jec_levels[pu_method]


def useJECFromDB(process, db):
    process.load("CondCore.DBCommon.CondDBCommon_cfi")

    process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(0)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(),

            connect = cms.string('sqlite:%s' % db)
         
            )

    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

def checkForTag(db_file, tag):
    import sqlite3

    db_file = db_file.replace('sqlite:', '')

    connection = sqlite3.connect(db_file)
    
    res = connection.execute('select TAG_NAME from IOV where TAG_NAME=?', tag).fetchall()

    return len(res) != 0

def appendJECToDB(process, payload, prefix):

    for set in process.jec.toGet:
        if set.label == payload:
            return

    tag = 'JetCorrectorParametersCollection_%s_%s' % (prefix, payload)
    if not checkForTag(process.jec.connect.value(), (tag,)):
        print("WARNING: The JEC payload %r is not present in the database you want to use. Corrections for this payload will be loaded from the Global Tag" % payload)
        return

    process.jec.toGet += [cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string(tag),
            label  = cms.untracked.string(payload)
            )]




def createProcess(isMC, globalTag, muonTypeID, runPuppiMuonIso, muonIsoCone, electronTypeID, tauTypeID, dropAnalyzerDumpEDM, runMVAPUPPETAnalysis,applyZSelections,applyWSelections,applyJECtoPuppiJets,runPuppiDiagnostics):

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

    process = cms.Process("JRA")

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
    process.load("Configuration.EventContent.EventContent_cff")
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')

    process.GlobalTag.globaltag = globalTag

    run_on_25ns = True
    
    # Custom JEC
    readJECFromDB = False

    jec_database = 'PY8_RunIISpring15DR74_bx25_MC.db'
    if not run_on_25ns:
        jec_database = 'PY8_RunIISpring15DR74_bx50_MC.db'

    jec_db_prefix = 'PY8_RunIISpring15DR74_bx25_MC'
    if not run_on_25ns:
        jec_db_prefix = 'PY8_RunIISpring15DR74_bx50_MC'

    if readJECFromDB:
        useJECFromDB(process, jec_database)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Input
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    process.source = cms.Source("PoolSource")

    process.load('CommonTools.UtilAlgos.TFileService_cfi')
    process.TFileService.fileName = cms.string('output_mc.root') if isMC else cms.string('output_data.root')
    process.TFileService.closeFileFast = cms.untracked.bool(True)

    # Jet corrections
    process.load('JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff')

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

    # Jet corrections
    process.load('JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff')

    ## loop on the jet collections : generic container just defined by clustering algorithm and cone dimension
    for name, params in jetsCollections.items():
        ## loop on the pileup methos
        for index, pu_method in enumerate(params['pu_methods']):
            # Add the jet collection

            jec_payload = get_jec_payload(params['algo'], pu_method)
            jec_levels  = get_jec_levels(pu_method)

            if readJECFromDB:
                appendJECToDB(process, jec_payload, jec_db_prefix)

            jetToolbox(process, params['algo'], 'dummy', 'out', PUMethod = pu_method, JETCorrPayload = jec_payload, JETCorrLevels = jec_levels, addPUJetID = False)

            algo          = params['algo'].upper()
            jetCollection = '%sPFJets%s' % (params['algo'], pu_method)
            postfix       = '%sPF%s' % (algo, pu_method)

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


            if applyJECtoPuppiJets == False:
                applyPostfix(process, 'patJets', postfix).addJetCorrFactors = cms.bool(False)
                applyPostfix(process, 'patJets', postfix).jetCorrFactorsSource = cms.VInputTag(cms.InputTag(""))

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
                            src = "slimmedMuons",
                            src_neutral_hadron = 'pfWeightedNeutralHadrons',
                            src_photon         = 'pfWeightedPhotons',
                            coneR = muonIsoCone
                            )
    if runPuppiMuonIso :

        ### PUPPI weighted isolation
        load_muonPFiso_sequence(process, 
                                'MuonPFIsoSequencePUPPI', 
                                algo = 'R04PUPPI',
                                src =  "slimmedMuons",
                                src_charged_hadron = 'pfAllChargedHadronsPuppi',
                                src_neutral_hadron = 'pfAllNeutralHadronsPuppi',
                                src_photon         = 'pfAllPhotonsPuppi',
                                coneR = muonIsoCone
                                )
    
        ### PUPPI weighted isolation without muons
        load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPINoMu', algo = 'R04PUPPINoMu',
                                src =  "slimmedMuons",
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

        if muonTypeID == "Tight" :

            applyMuonID(process, label = "Tight", src  = "slimmedMuons", type = 'tightID',
                        iso_map_charged_hadron  = '', 
                        iso_map_neutral_hadron  = '', 
                        iso_map_photon          = '',                 
                        rho = 'fixedGridRhoFastjetAll',
                        typeIsoVal = 1
                        )
        elif muonTypeID == "TightDBeta" :    

            applyMuonID(process, label = "TightDBeta", src  = "slimmedMuons", type = 'tightID',
                        iso_map_charged_hadron  = '', 
                        iso_map_neutral_hadron  = 'muPFIsoValueNHR04PFWGT', 
                        iso_map_photon          = 'muPFIsoValuePhR04PFWGT',                 
                        rho = 'fixedGridRhoFastjetAll',
                        typeIsoVal = 3
                        )

        elif muonTypeID == "TightPuppiNoMu" :    
            if not runPuppiMuonIso :
                print "no available puppi isolation, please check runPuppiMuonIso option";
                exit();
            
            applyMuonID(process, label = "TightPuppiNoMu", src  = "slimmedMuons", type = 'tightID',
                        iso_map_charged_hadron  = 'muPFIsoValueCHR04PUPPINoMu', 
                        iso_map_neutral_hadron  = 'muPFIsoValueNHR04PUPPINoMu', 
                        iso_map_photon          = 'muPFIsoValuePhR04PUPPINoMu',                 
                        rho = 'fixedGridRhoFastjetAll',
                        typeIsoVal = 4
                        )

        else:
                exit("not recognized muon type ID --> exit");

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

        if electronTypeID == "Tight" :

            applyElectronID(process, label = "Tight", src  = "slimmedElectrons", 
                            iso_map_charged_hadron  = '', 
                            iso_map_neutral_hadron  = '', 
                            iso_map_photon          = '',                 
                            rho = 'fixedGridRhoFastjetAll',
                            typeIsoVal = 1,
                            electron_id_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight'
                            )

        elif electronTypeID == "Medium" :

            applyElectronID(process, label = "Medium", src  = "slimmedElectrons", 
                            iso_map_charged_hadron  = '', 
                            iso_map_neutral_hadron  = '', 
                            iso_map_photon          = '',                 
                            rho = 'fixedGridRhoFastjetAll',
                            typeIsoVal = 1,
                            electron_id_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium'
                        )

        else:
            exit("not recognized electron type ID --> exit");
            
    
    ##################
    #### TAU ID ######
    ##################

    if runMVAPUPPETAnalysis :
        from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyTauID

        if tauTypeID != "":
            
            applyTauID( process, label = tauTypeID, src = "slimmedTaus",
                        muonCollection     = "slimmedMuons"+muonTypeID,
                        electronCollection = "slimmedElectrons"+electronTypeID)
    
    ##################
    ### clean jets ###
    ##################

    ## define gen leptons (muons and electrons)\        
    process.packedGenLeptons = cms.EDFilter("CandPtrSelector",
                                            cut = cms.string('(abs(pdgId) = 11 || abs(pdgId) = 13) && pt > 10'),
                                            src = cms.InputTag("packedGenParticles")
                                            )
    

        
    if runMVAPUPPETAnalysis :

        from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import cleanJetsFromLeptons, cleanGenJetsFromGenLeptons

        for name, params in jetsCollections.items():
        ## loop on the pileup methos                                                                                                                                            
            for index, pu_method in enumerate(params['pu_methods']):
                postfix       = '%sPF%s' % (algo, pu_method)

                setattr(getattr(process,"selectedPatJets"+postfix),"cut",cms.string('pt > 30'))
                
                cleanJetsFromLeptons(process,"Cleaned"+"Mu"+muonTypeID+"Ele"+electronTypeID,                                 
                                     jetCollection      = "selectedPatJets"+postfix,
                                     muonCollection     = "slimmedMuons"+muonTypeID,
                                     electronCollection = "slimmedElectrons"+electronTypeID,
                                     tauCollection      = "",
                                     jetPtCut   = 30.,
                                     jetEtaCut  = 5.,
                                     dRCleaning = 0.3) 

            ## clean also the related Gen Jet Collection
            genAlgo = algo.replace("AK","ak")
            setattr(process,"pat"+genAlgo+"GenJetsNoNu", cms.EDProducer("PATJetProducer",
                                                                        jetSource = cms.InputTag(genAlgo+"GenJetsNoNu"),
                                                                        addJetCorrFactors = cms.bool(False),
                                                                        addJetCharge = cms.bool(False),
                                                                        addGenJetMatch = cms.bool(False),
                                                                        embedGenJetMatch = cms.bool(False),
                                                                        addAssociatedTracks = cms.bool(False),
                                                                        addBTagInfo = cms.bool(False),
                                                                        partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
                                                                        addGenPartonMatch = cms.bool(False),
                                                                        addTagInfos = cms.bool(False),
                                                                        addPartonJetMatch = cms.bool(False),
                                                                        embedGenPartonMatch = cms.bool(False),
                                                                        useLegacyJetMCFlavour = cms.bool(False),
                                                                        addEfficiencies = cms.bool(False),
                                                                        embedPFCandidates = cms.bool(False),
                                                                        addJetFlavourInfo = cms.bool(False),
                                                                        addResolutions = cms.bool(False),
                                                                        getJetMCFlavour = cms.bool(False),
                                                                        addDiscriminators = cms.bool(False),
                                                                        addJetID = cms.bool(False)
             ))

            setattr(process,"selectedPat"+genAlgo+"GenJetsNoNu",cms.EDFilter("PATJetSelector",
                                                                             cut = cms.string('pt > 30'),
                                                                             src = cms.InputTag("pat"+genAlgo+"GenJetsNoNu")
                                                                             ))

            cleanGenJetsFromGenLeptons (process,
                                        jetCollection       = "selectedPat"+genAlgo+"GenJetsNoNu",
                                        genLeptonCollection = "packedGenLeptons",
                                        jetPtCut         = 30,
                                        jetEtaCut        = 5,
                                        dRCleaning       = 0.3)


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

    ## Raw PF METs
    process.load('RecoMET.METProducers.PFMET_cfi')
    process.pfMet.src = cms.InputTag('packedPFCandidates')
    addMETCollection(process, labelName='patPFMet', metSource='pfMet') # RAW MET
    process.patPFMet.addGenMET = False

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
    ### Standard
    if not hasattr(process, 'ak4PFJets'):
        print("WARNING: No AK4 jets produced. Type 1 corrections for MET are not available.")
    else:
        process.corrPfMetType1 = corrPfMetType1.clone(
            src = 'ak4PFJets',
            jetCorrLabel = 'ak4PFL1FastL2L3Corrector',
            offsetCorrLabel = 'ak4PFL1FastjetCorrector'
        )
        process.pfMetT1 = pfMetT1.clone(
            src = 'pfMet',
            srcCorrections = [ cms.InputTag("corrPfMetType1","type1") ]
        )

        addMETCollection(process, labelName='patMET', metSource='pfMetT1') # T1 MET
        process.patMET.addGenMET = False

    ### PUPPI
    if not hasattr(process, 'ak4PFJetsPuppi'):
        print("WARNING: No AK4 puppi jets produced. Type 1 corrections for puppi MET are not available.")
    else:
        if applyJECtoPuppiJets :
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

    ### CHS TypeI corrected
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

        addMETCollection(process, labelName='patMETCHS', metSource='pfMetT1CHS') # T1 CHS MET
        process.patMETCHS.addGenMET = False


    ## Slimmed METs
    from PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi import slimmedMETs
    #### CaloMET is not available in MiniAOD
    del slimmedMETs.caloMET

    ### PUPPI : make slimmed METs in order to embed both corrected and not corrected one after TypeI
    ### Standard
    process.slimmedMETs = slimmedMETs.clone()
    if hasattr(process, 'patMET'):
        # Create MET from Type 1 PF collection
        process.patMET.addGenMET = True
        process.slimmedMETs.src = cms.InputTag("patMET")
        process.slimmedMETs.rawUncertainties = cms.InputTag("patPFMet") # only central value
    else:
        # Create MET from RAW PF collection
        process.patPFMet.addGenMET = True
        process.slimmedMETs.src = cms.InputTag("patPFMet")
        del process.slimmedMETs.rawUncertainties # not available

    del process.slimmedMETs.type1Uncertainties # not available
    del process.slimmedMETs.type1p2Uncertainties # not available

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
        process.slimmedMETsPuppi.rawUncertainties = cms.InputTag("patPFMetPuppi") # only central value

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


    ## create the Path
    process.jmfw_analyzers = cms.Sequence()
    process.p = cms.Path(process.jmfw_analyzers)


    if runMVAPUPPETAnalysis :

        ######## add other specific set of particles and MET collections, in this case not TypeI corrected, but we will still use the same workflow    
        ## all charge particles from PUPPI : hadrons + leptons (e,mu,tau) --> trak met

        ## particles for chs 
        process.pfChargedPV = cms.EDFilter("CandPtrSelector",
                                           src = cms.InputTag("chs"),
                                           cut = cms.string("pt > 0  && charge!=0"))

        process.pfNeutrals = cms.EDFilter("CandPtrSelector",
                                          src = cms.InputTag("chs"),
                                          cut = cms.string("pt > 0 && charge == 0"))


        process.pfChargedPU = cms.EDFilter("CandPtrSelector",
                                           cut = cms.string('!fromPV'),
                                           src = cms.InputTag("packedPFCandidates")
                                           )

        #### puppi particles
        process.pfPuppiAll = cms.EDFilter("CandPtrSelector",
                                          src = cms.InputTag("puppi"),
                                          cut = cms.string("pt > 0"))

        process.pfAllChargedParticlesPuppi = cms.EDFilter("CandPtrSelector",
                                                          src = cms.InputTag("puppi"),
                                                          cut = cms.string("charge !=0 && pt > 0"))

        process.pfAllNeutralParticlesPuppi  = cms.EDFilter("CandPtrSelector", 
                                                           src = cms.InputTag("puppi"), 
                                                           cut = cms.string("charge == 0 && pt > 0"))

        ## inverted puppi
        process.pfPUPuppi = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag("pupuppi"),
                                         cut = cms.string("pt > 0"))

        process.pfPUPuppiCharge = cms.EDFilter("CandPtrSelector",
                                               src = cms.InputTag("pupuppi"),
                                               cut = cms.string("pt > 0 && charge != 0"))

        process.pfAllNeutralParticlesPuppiPU  = cms.EDFilter("CandPtrSelector", 
                                                             src = cms.InputTag("pupuppi"), 
                                                             cut = cms.string("charge=0 && pt > 0"))
    

        #### Missing energies        
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
        process.pfMetChargedPU.src   = cms.InputTag("pfChargedPU")
        process.pfMetChargedPU.alias = cms.string('pfMetPuppiChargedPU')
        addMETCollection(process, labelName='patPFMetChargedPU', metSource='pfMetChargedPU') 
        process.patPFMetChargedPU.addGenMET = True

        process.slimmedMETsPuppiChargedPU = slimmedMETs.clone()
        process.slimmedMETsPuppiChargedPU.src = cms.InputTag("patPFMetChargedPU")
        process.slimmedMETsPuppiChargedPU.rawUncertainties = cms.InputTag("patPFMetChargedPU") # only central value
        
        del process.slimmedMETsPuppiChargedPU.type1Uncertainties # not available                                                                                           
        del process.slimmedMETsPuppiChargedPU.type1p2Uncertainties # not available                                                                                      
        
        ## neutral particles from PV starting from puppi particles                
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

        
        ## look for a good ZLL candidate i.e.: two tight muons within Mz, two tight electron or taus in Mz
        if not applyWSelections :
            setattr(process,"ZdiMuon"+muonTypeID,cms.EDProducer("CandViewCombiner",
                                                                decay       = cms.string("slimmedMuons"+muonTypeID+"@+ "+"slimmedMuons"+muonTypeID+"@-"),
                                                                checkCharge = cms.bool(True),
                                                                cut         = cms.string("mass > 70 && mass < 110 & charge=0"),                                              
                                                                ))
        
            process.jmfw_analyzers += getattr(process,"ZdiMuon"+muonTypeID);
        
            setattr(process,"ZdiElectron"+electronTypeID,cms.EDProducer("CandViewCombiner",
                                                                       decay       = cms.string("slimmedElectrons"+electronTypeID+"@+ "+"slimmedElectrons"+electronTypeID+"@-"),
                                                                       checkCharge = cms.bool(True),
                                                                       cut         = cms.string("mass > 70 && mass < 110 & charge=0"),                                       
                                                                       ))
        
            process.jmfw_analyzers += getattr(process,"ZdiElectron"+electronTypeID);

            if tauTypeID != "":

                setattr(process,"ZdiTau"+tauTypeID,cms.EDProducer("CandViewCombiner",
                                                                  decay       = cms.string("slimmedTaus"+tauTypeID+"Cleaned@+ "+"slimmedTaus"+tauTypeID+"Cleaned@-"),
                                                                  checkCharge = cms.bool(True),
                                                                  cut         = cms.string("mass > 70 && mass < 110 & charge=0"),                                               
                                                                  ))
            
                process.jmfw_analyzers += getattr(process,"ZdiTau"+tauTypeID);


                ## merge all the Z canddates and ask for only one candidate per event
                setattr(process,"ZdiLepton", cms.EDProducer("CandViewMerger",
                                                            src = cms.VInputTag("ZdiMuon"+muonTypeID,"ZdiElectron"+electronTypeID,"ZdiTau"+tauTypeID)
                                                            ))
            else:


                ## merge all the Z canddates and ask for only one candidate per event
                setattr(process,"ZdiLepton", cms.EDProducer("CandViewMerger",
                                                            src = cms.VInputTag("ZdiMuon"+muonTypeID,"ZdiElectron"+electronTypeID)
                                                            ))

                process.jmfw_analyzers += getattr(process,"ZdiLepton");
    
            ## filter number of Z candidates
            if applyZSelections :

                setattr(process,"ZdiLeptonFilter",cms.EDFilter("PATCandViewCountFilter",
                                                               minNumber = cms.uint32(1),
                                                               maxNumber = cms.uint32(1),
                                                               src = cms.InputTag("ZdiLepton")
                                                               ))

                process.jmfw_analyzers += getattr(process,"ZdiLeptonFilter");

            ## now merge all the previous leptons in one single collection and ask for no more than 2 tight leptons
            if tauTypeID != "" :
                setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                              src = cms.VInputTag("slimmedMuons"+muonTypeID,"slimmedElectrons"+electronTypeID,"slimmedTaus"+tauTypeID+"Cleaned")
                                                              ))

            else:
                setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                              src = cms.VInputTag("slimmedMuons"+muonTypeID,"slimmedElectrons"+electronTypeID)
                                                              ))

            process.jmfw_analyzers += getattr(process,"LeptonMerge");
            
            ## apply selections on Nleptons
            if applyZSelections :
                setattr(process,"LeptonMergeFilter",cms.EDFilter("PATCandViewCountFilter",
                                                                 minNumber = cms.uint32(2),
                                                                 maxNumber = cms.uint32(2),
                                                                 src = cms.InputTag("LeptonMerge")
                                                                 ))
                process.jmfw_analyzers += getattr(process,"LeptonMergeFilter");
                

        ### apply simpler selection for W+jets events
        else:

            if tauTypeID != "" :
                setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                              src = cms.VInputTag("slimmedMuons"+muonTypeID,"slimmedElectrons"+electronTypeID,"slimmedTaus"+tauTypeID+"Cleaned")
                                                              ))

            else:
                setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                              src = cms.VInputTag("slimmedMuons"+muonTypeID,"slimmedElectrons"+electronTypeID)
                                                              ))

            process.jmfw_analyzers += getattr(process,"LeptonMerge");
            
            setattr(process,"LeptonMergeFilter",cms.EDFilter("PATCandViewCountFilter",
                                                             minNumber = cms.uint32(1),
                                                             maxNumber = cms.uint32(1),
                                                             src = cms.InputTag("LeptonMerge")
                                                             ))
        
            process.jmfw_analyzers += getattr(process,"LeptonMergeFilter");
            

        ### in this way we are sure to have a Z candidate with no more than 2 leptons in the event (after selections)
        jetColl = "selectedPatJetsAK4PFPuppiCleaned"+"Mu"+muonTypeID+"Ele"+electronTypeID;

        setattr(process,"mvaPUPPET", cms.EDProducer("mvaPUPPET",
                                                    srcMETs      = cms.VInputTag("slimmedMETs","slimmedMETsCHS", "slimmedMETsPuppi","slimmedMETsPuppiChargedPV","slimmedMETsPuppiChargedPU","slimmedMETsPuppiNeutralPV","slimmedMETsPuppiNeutralPU"),
                                                    referenceMET   = cms.InputTag("slimmedMETsPuppi"),
                                                    srcVertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                    srcJets        = cms.InputTag(jetColl),
                                                    inputFileNames = cms.PSet(
                    #PhiCorrectionWeightFile = cms.FileInPath('JMEAnalysis/JMEValidator/data/PhiCor_13TeV.root'),
                    #RecoilCorrectionWeightFile  = cms.FileInPath('JMEAnalysis/JMEValidator/data/RecoilCor_13TeV.root')
                   ),
                                                    srcLeptons = cms.VInputTag("LeptonMerge") ))

        process.jmfw_analyzers += getattr(process,"mvaPUPPET");
        

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

        # Jets
        for name, params in jetsCollections.items():
            for index, pu_method in enumerate(params['pu_methods']):

                algo = params['algo'].upper()
                jetCollection = 'selectedPatJets%sPF%s' % (algo, pu_method)
                
                print('Adding analyzer for jets collection \'%s\'' % jetCollection)
                
                 # FIXME: Remove once PUPPI and SK payloads are in the GT
                jec_payload = get_jec_payload(algo, pu_method)
                jec_levels = get_jec_levels(pu_method)

                if jec_payload == 'None':
                    jec_payload = '';
                    
                if jec_levels == ['None']:
                    jec_levels = []
                    
                analyzer = cms.EDAnalyzer('JMEJetAnalyzer',
                                          JetAnalyserCommonParameters,
                                          JetCorLabel   = cms.string(jec_payload),
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

        setattr(process, "PUPPET", 
                cms.EDAnalyzer('PUPPETAnalyzer',
                               srcJet    = cms.InputTag("selectedPatJetsAK4PFPuppiCleaned"+"Mu"+muonTypeID+"Ele"+electronTypeID),
                               srcVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               srcZboson = cms.InputTag("mvaPUPPET","ZtagBoson"),
                               srcGenMet = cms.InputTag("slimmedMETs","","PAT"),
                               srcGenJets          = cms.InputTag("slimmedGenJets","","PAT"),
                               srcGenJetsCleaned   = cms.InputTag("selectedPatak4GenJetsNoNuCleaned"),
                               srcGenParticles     = cms.InputTag("prunedGenParticles","","PAT"),
                               srcRecoilPFMet      = cms.InputTag("mvaPUPPET","recoilslimmedMETs"),
                               srcRecoilPFCHSMet   = cms.InputTag("mvaPUPPET","recoilslimmedMETsCHS"),
                               srcRecoilPFPuppiMet = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppi"),
                               srcRecoilPFPuppiMet_ChargedPV = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiChargedPV"),
                               srcRecoilPFPuppiMet_ChargedPU = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiChargedPU"),
                               srcRecoilPFPuppiMet_NeutralPV = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiNeutralPV"),
                               srcRecoilPFPuppiMet_NeutralPU = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiNeutralPU"),
                               srcMVAMet     = cms.InputTag("mvaPUPPET","mvaMET"),
                               dRgenMatching = cms.double(0.3)))
        

        process.jmfw_analyzers += getattr(process,"PUPPET")
                                           

    # Puppi ; only for the first 1000 events of the job
    ## Turn on diagnostic
    if runPuppiDiagnostics :
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

    return process

    #!
    #! THAT'S ALL! CAN YOU BELIEVE IT? :-D
    #!

