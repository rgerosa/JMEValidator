import FWCore.ParameterSet.Config as cms
import sys

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
                'Puppi': 'AK%dPFPuppi',
                'CHS': 'AK%dPFchs',
                'SK': 'AK%dPFchs',
                '': 'AK%dPF',
                }
    
    
    cone_size = get_cone_size(algo)
    
    if not pu_method in jec_payloads:
        print('WARNING: JEC payload not found for method %r. Using default one.' % pu_method)
        return 'None'
    
    return jec_payloads[pu_method] % cone_size

def get_jec_levels(pu_method, isMC = True, useJECFromDB = False):

    if isMC :

            jec_levels = {
                'Puppi': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'CHS': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'SK': ['L2Relative', 'L3Absolute'],
                '': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                }

    else:

        jec_levels = {
            'Puppi': ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
            'CHS': ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
            'SK': ['L2Relative', 'L3Absolute','L2L3Residual'],
            '': ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
            }
    
        
    if not pu_method in jec_levels:
        print('WARNING: JEC levels not found for method %r. Using default ones.' % pu_method)
        return ['None']
    
    return jec_levels[pu_method]


def useJECFromDB(process, db, postfix = ""):

    process.load("CondCore.DBCommon.CondDBCommon_cfi")

    setattr(process,"jec"+postfix, cms.ESSource("PoolDBESSource",
                                                DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
                                                timetype = cms.string('runnumber'),
                                                toGet = cms.VPSet(),                                        
                                                connect = cms.string('sqlite:%s' % db)))

    setattr(process,"es_prefer_jec"+postfix, cms.ESPrefer('PoolDBESSource','jec'+postfix))

def checkForTag(db_file, tag):

    import sqlite3

    db_file = db_file.replace('sqlite:', '')

    connection = sqlite3.connect(db_file)
    
    res = connection.execute('select TAG_NAME from IOV where TAG_NAME=?', tag).fetchall()

    return len(res) != 0

def appendJECToDB(process, payload, prefix, postfix=""):

    instance = getattr(process,"jec"+postfix);

    for set in instance.toGet:
        if set.label == payload:
            return

    tag = 'JetCorrectorParametersCollection_%s_%s' % (prefix, payload)
    tag = tag.replace('.db','')

    if not checkForTag(instance.connect.value(), (tag,)):
        print("WARNING: The JEC payload %r is not present in the database you want to use. Corrections for this payload will be loaded from the Global Tag" % payload)
        return


    instance.toGet += [cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string(tag),
            label  = cms.untracked.string(payload)
            )]


def createProcess(isMC, ## isMC flag
                  globalTag, ## global tag
                  muonTypeID, runPuppiMuonIso, muonIsoCone, ## muons
                  electronTypeID, ## electrons
                  tauTypeID, ## taus
                  dropAnalyzerDumpEDM, ## debug EDM 
                  runMVAPUPPETAnalysis,applyZSelections,applyWSelections, ## special settings for PUPPET
                  jetPtCut,
                  applyJECtoPuppiJets,
                  runPuppiDiagnostics,
                  isRunningOn25ns,
                  useJECFromLocalDB,
                  etaCutForMetDiagnostic,
                  ptNeutralCut,
                  ptNeutralCutSlope,
                  etaBinPuppi,
                  puppiCone,
                  puppiUseCharge,
                  ptThresholdForTypeIPuppi,
                  runPUPPINoLeptons):

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

    process.GlobalTag.globaltag = globalTag

    jec_database_PF = 'Summer15_50nsV2_DATA.db'
    if not isRunningOn25ns:
        jec_database_PF = 'Summer15_50nsV2_DATA.db'

    jec_database_Puppi = 'Summer15_50nsV2_DATA.db'

    if not isRunningOn25ns:
        jec_database_Puppi = 'Summer15_50nsV2_DATA.db'


    if useJECFromLocalDB:
        useJECFromDB(process, jec_database_PF)
        useJECFromDB(process, jec_database_Puppi,"_puppi")

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Input
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    process.source = cms.Source("PoolSource")

    process.load('CommonTools.UtilAlgos.TFileService_cfi')
    process.TFileService.fileName = cms.string('output_mc.root') if isMC else cms.string('output_data.root')
    process.TFileService.closeFileFast = cms.untracked.bool(True)

    
    ## count the number of events
    process.AllEvents = cms.EDFilter("PassFilter",
        srcGenEventInfo     = cms.InputTag("generator"),
        isMC      = cms.bool(isMC),
      )
    process.counterPath = cms.Path(process.AllEvents)

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

    if isMC:
        
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
        
    else:

        jetsCollections = {
            'AK4': {
                'algo': 'ak4',
                'pu_methods': ['Puppi', 'CHS', ''],
                'jec_payloads': ['AK4PFPUPPI', 'AK4PFchs', 'AK4PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
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
            jec_levels  = get_jec_levels(pu_method,isMC,useJECFromLocalDB)

            if useJECFromLocalDB:
                if pu_method == "Puppi":
                    appendJECToDB(process, jec_payload, jec_database_Puppi.replace("_PUPPI","").replace("_Puppi",""),"_puppi")
                else:
                    appendJECToDB(process, jec_payload, jec_database_PF.replace("_PFCHS","").replace("_PF",""))

            jetToolbox(process, params['algo'], 'dummy', 'out', runOnMC=isMC, PUMethod = pu_method, JETCorrPayload = jec_payload, JETCorrLevels = jec_levels, addPUJetID = False)

            if useJECFromLocalDB and pu_method == "Puppi":
                getattr(process,"patJetCorrFactors"+name+"PFPuppi").payload = cms.string(jec_payload)

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


            if applyJECtoPuppiJets == False and pu_method == 'Puppi':
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
    

            

    # Create METs from CHS and PUPPI
    from PhysicsTools.PatAlgos.tools.metTools import addMETCollection

    ## Gen MET ###
    ### Copied from https://github.com/cms-sw/cmssw/blob/2b75137e278b50fc967f95929388d430ef64710b/RecoMET/Configuration/python/GenMETParticles_cff.py#L37

    if isMC :
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

    ## make a selection on eta according to the value defined
    process.pfCandidatesForMET = cms.EDFilter("CandPtrSelector",
                                              src = cms.InputTag("packedPFCandidates"),
                                              cut = cms.string("abs(eta) < %f"%etaCutForMetDiagnostic))


    process.load('RecoMET.METProducers.PFMET_cfi')
    process.pfMet.src = cms.InputTag('pfCandidatesForMET')
    addMETCollection(process, labelName='patPFMet', metSource='pfMet') # RAW MET
    process.patPFMet.addGenMET = False

    ## CHS pat MET; raw PF is the slimmedMet in miniAOD + typeI correction
    from RecoMET.METProducers.PFMET_cfi import pfMet
    process.pfCandidatesForMETCHS = cms.EDFilter("CandPtrSelector",
                                                 src = cms.InputTag("chs"),
                                                 cut = cms.string("abs(eta) < %f"%etaCutForMetDiagnostic))


    process.pfMetCHS        = pfMet.clone()
    process.pfMetCHS.src    = cms.InputTag("pfCandidatesForMETCHS") ## packed candidates without fromPV < 1
    process.pfMetCHS.alias  = cms.string('pfMetCHS')
    addMETCollection(process, labelName='patPFMetCHS', metSource='pfMetCHS') # Convert the CHS PFMet in PAT MET
    process.patPFMetCHS.addGenMET = False

    ## Type 1 corrections
    from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
    from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

    ######
    ### Standard TypeI correction
    ######

    from CommonTools.RecoAlgos.pfJetSelector_cfi import pfJetSelector

    if not hasattr(process, 'ak4PFJets'):
        print("WARNING: No AK4 jets produced. Type 1 corrections for MET are not available.")
    else:        

        
        process.ak4PFJetsForTypeI = pfJetSelector.clone(
            src = cms.InputTag( "ak4PFJets" ),
            cut = cms.string( "abs(eta)< %f"%(etaCutForMetDiagnostic) )
                )
        
        if isMC :

            process.corrPfMetType1 = corrPfMetType1.clone(
                src = 'ak4PFJetsForTypeI',
                jetCorrLabel = 'ak4PFL1FastL2L3Corrector',
                offsetCorrLabel = 'ak4PFL1FastjetCorrector',
                )
        else:

            process.corrPfMetType1 = corrPfMetType1.clone(
                src = 'ak4PFJetsForTypeI',
                jetCorrLabel = 'ak4PFL1FastL2L3ResidualCorrector',
                offsetCorrLabel = 'ak4PFL1FastjetCorrector',
                )

        
        
        process.pfMetT1 = pfMetT1.clone(
            src = 'pfMet',
            srcCorrections = [ cms.InputTag("corrPfMetType1","type1") ]
            )
        
        addMETCollection(process, labelName='patMET', metSource='pfMetT1') # T1 MET
        process.patMET.addGenMET = False

            
    ### CHS TypeI corrected
    if not hasattr(process, 'ak4PFJetsCHS'):
        print("WARNING: No AK4 CHS jets produced. Type 1 corrections for CHS MET are not available.")
    else:

        process.ak4PFJetsCHSForTypeI = pfJetSelector.clone(
            src = cms.InputTag( "ak4PFJetsCHS" ),
            cut = cms.string( "abs(eta)<%f"%(etaCutForMetDiagnostic) )
            )

        if isMC:
            process.corrPfMetType1CHS = corrPfMetType1.clone(
                src             = 'ak4PFJetsCHSForTypeI',
                jetCorrLabel    = 'ak4PFCHSL1FastL2L3Corrector',
                offsetCorrLabel = 'ak4PFCHSL1FastjetCorrector'
                )
        else:
            process.corrPfMetType1CHS = corrPfMetType1.clone(
                src             = 'ak4PFJetsCHSForTypeI',
                jetCorrLabel    = 'ak4PFCHSL1FastL2L3ResidualCorrector',
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
        process.patMET.addGenMET = isMC
        process.slimmedMETs.src = cms.InputTag("patMET")
        process.slimmedMETs.rawUncertainties = cms.InputTag("patPFMet") # only central value
    else:
        # Create MET from RAW PF collection
        process.patPFMet.addGenMET = isMC
        process.slimmedMETs.src = cms.InputTag("patPFMet")
        del process.slimmedMETs.rawUncertainties # not available

    del process.slimmedMETs.type1Uncertainties # not available
    del process.slimmedMETs.type1p2Uncertainties # not available


    ### CHS
    process.slimmedMETsCHS = slimmedMETs.clone()
    if hasattr(process, "patMETCHS"):
        # Create MET from Type 1 PF collection
        process.patMETCHS.addGenMET = isMC
        process.slimmedMETsCHS.src = cms.InputTag("patMETCHS")
        process.slimmedMETsCHS.rawUncertainties = cms.InputTag("patPFMetCHS") # only central value
    else:
        # Create MET from RAW PF collection
        process.patPFMetCHS.addGenMET = isMC
        process.slimmedMETsCHS.src = cms.InputTag("patPFMetCHS")
        del process.slimmedMETsCHS.rawUncertainties # not available

    del process.slimmedMETsCHS.type1Uncertainties # not available
    del process.slimmedMETsCHS.type1p2Uncertainties # not available


    ## create the Path
    process.jmfw_analyzers = cms.Sequence()
    process.p = cms.Path(process.jmfw_analyzers)

    if runMVAPUPPETAnalysis :

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
                      jetCollectionPuppi = "selectedPatJetsAK4PFPuppi", 
                      jetCollectionPF    = "selectedPatJetsAK4PF", 
                      dRCleaning = 0.3, 
                      jetPtCut = jetPtCut, 
                      jetEtaCut = 5.,
                      etaCutForMetDiagnostic = etaCutForMetDiagnostic,
                      genJetCollection = "ak4GenJetsNoNu",
                      cleanGenJets = True,
                      applyTypeICorrection = applyJECtoPuppiJets, 
                      useJECFromLocalDB = useJECFromLocalDB,                      
                      applyZSelections = applyZSelections, 
                      applyWSelections = applyWSelections,
                      runPUPPINoLeptons = runPUPPINoLeptons,
                      )
    

        ## change cone and use charge for the tracking region
        if len(ptNeutralCut) !=3 or len(ptNeutralCutSlope)!=3 or len(etaBinPuppi)!=3 or len(puppiCone)!=3 or len(puppiUseCharge)!=3 :
            sys.exit("puppi parameters not corrected --> please check")

        process.puppi.producePackedCollection = cms.bool(True)

        for iBin in range(len(ptNeutralCut)):
            if iBin == 0 :
                process.puppi.algos[iBin].etaMin = cms.double(0.);
            else:
                process.puppi.algos[iBin].etaMin = cms.double(etaBinPuppi[iBin-1]);
 
            process.puppi.algos[iBin].etaMax            = cms.double(etaBinPuppi[iBin]);
            process.puppi.algos[iBin].MinNeutralPt      = cms.double(ptNeutralCut[iBin]);
            process.puppi.algos[iBin].MinNeutralPtSlope        = cms.double(ptNeutralCutSlope[iBin]);
            process.puppi.algos[iBin].puppiAlgos[0].cone       = cms.double(puppiCone[iBin])
            process.puppi.algos[iBin].puppiAlgos[0].useCharged = cms.bool(puppiUseCharge[iBin])


        ######## add other specific set of particles and MET collections, in this case not TypeI corrected, but we will still use the same workflow    
        ## all charge particles from PUPPI : hadrons + leptons (e,mu,tau) --> trak met
        ## particles for chs 

        process.pfChargedPV = cms.EDFilter("CandPtrSelector",
                                           src = cms.InputTag("chs"),
                                           cut = cms.string("pt > 0  && charge!=0 && abs(eta) < %f"%etaCutForMetDiagnostic))

        process.pfNeutrals = cms.EDFilter("CandPtrSelector",
                                          src = cms.InputTag("chs"),
                                          cut = cms.string("pt > 0 && charge == 0 && abs(eta) < %f"%etaCutForMetDiagnostic))


        process.pfChargedPU = cms.EDFilter("CandPtrSelector",
                                           cut = cms.string('!fromPV && abs(eta) < %f'%etaCutForMetDiagnostic),
                                           src = cms.InputTag("packedPFCandidates")
                                           )


        ## inverted puppi
        process.pfPuppiAll = cms.EDFilter("CandPtrSelector",
                                          src = cms.InputTag("puppi"),
                                          cut = cms.string("pt > 0 && abs(eta) < %f"%(etaCutForMetDiagnostic)))



        ## inverted puppi                                                                                                                                                      
        process.pfPUPuppi = cms.EDFilter("CandPtrSelector",
                                     src = cms.InputTag("pupuppi"),
                                     cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

        process.pfPUPuppiCharge = cms.EDFilter("CandPtrSelector",
                                               src = cms.InputTag("pupuppi"),
                                               cut = cms.string("pt > 0 && charge != 0 && abs(eta) < %f"%etaCutForMetDiagnostic))
        

        process.jmfw_analyzers += getattr(process,"mvaPUPPET");
        
        
    if dropAnalyzerDumpEDM:        
        return process
   
    # Run
    if isMC and not runMVAPUPPETAnalysis:
        process.run = cms.EDAnalyzer('RunAnalyzer')
        process.jmfw_analyzers += process.run
    

    if not runMVAPUPPETAnalysis:

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
                               isMC      = cms.bool(isMC),
                               srcJet    = cms.InputTag("selectedPatJetsAK4PFPuppiCleaned"),
                               srcJetPF  = cms.InputTag("selectedPatJetsAK4PFCleaned"),
                               srcVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               srcZboson = cms.InputTag("mvaPUPPET","ZtagBoson"),
                               srcLeptons = cms.InputTag("LeptonMerge"),
                               srcGenMet = cms.InputTag("slimmedMETs","","PAT"),
                               srcGenJets          = cms.InputTag("slimmedGenJets","","PAT"),
                               srcGenJetsCleaned   = cms.InputTag("selectedPatak4GenJetsNoNuCleaned"),
                               srcGenParticles     = cms.InputTag("prunedGenParticles","","PAT"),
                               srcGenEventInfo     = cms.InputTag("generator"),
                               srcPFMet            = cms.InputTag("slimmedMETs","","JRA"),
                               srcPFCHSMet         = cms.InputTag("slimmedMETsCHS","","JRA"),
                               srcPFPuppiMet       = cms.InputTag("slimmedMETsPuppi","","JRA"),
                               srcRecoilPFMet      = cms.InputTag("mvaPUPPET","recoilslimmedMETs"),
                               srcRecoilPFCHSMet   = cms.InputTag("mvaPUPPET","recoilslimmedMETsCHS"),
                               srcRecoilPFPuppiMet = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppi"),
                               srcRecoilPFPuppiMet_ChargedPV = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiChargedPV"),
                               srcRecoilPFPuppiMet_ChargedPU = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiChargedPU"),
                               srcRecoilPFPuppiMet_NeutralPV = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiNeutralPV"),
                               srcRecoilPFPuppiMet_NeutralPU = cms.InputTag("mvaPUPPET","recoilslimmedMETsPuppiNeutralPU"),
                               srcRecoilPFChargedPVNeutralPVPUJetID = cms.InputTag("mvaPUPPET","recoilslimmedMETsChargedPVNeutralPVPUJetID"),
                               srcRecoilPFChargedPUNeutralPUPUJetID = cms.InputTag("mvaPUPPET","recoilslimmedMETsChargedPUNeutralPUPUJetID"),
                               srcRecoilPFChargedPVNeutralPV        = cms.InputTag("mvaPUPPET","recoilslimmedMETsChargedPVNeutralPV"),
                               srcMVAMet     = cms.InputTag("mvaPUPPET","mvaMET"),
                               dRgenMatching = cms.double(0.3),
                               srcMetFiltersBits = cms.InputTag("TriggerResults","","PAT"),
                               srcTriggerBits = cms.InputTag("TriggerResults","","HLT"),
                               srcTriggerPrescales = cms.InputTag('patTrigger')
                               )
                )
        
        if not isMC:
            getattr(process,"PUPPET").srcMetFiltersBits = cms.InputTag("TriggerResults","","RECO")

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
    
