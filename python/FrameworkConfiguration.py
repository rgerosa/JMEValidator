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
                #'Puppi': 'AK%dPFPuppi',
                'CHS': 'AK%dPFchs',
                #'SK': 'AK%dPFchs',
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
                #'Puppi': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'CHS': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                #'SK': ['L2Relative', 'L3Absolute'],
                '': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                }

    else:

        jec_levels = {
            #'Puppi': ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
            'CHS': ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
            #'SK': ['L2Relative', 'L3Absolute','L2L3Residual'],
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
                  processName,
                  globalTag, ## global tag
                  muonTypeID, muonIsoCone, ## muons
                  electronTypeID, ## electrons
                  tauTypeID, ## taus
                  applyZSelections,applyWSelections, ## special settings for PUPPET
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

    if isMC:
        
        jetsCollections = {
            'AK4': {
                'algo': 'ak4',
                'pu_methods': ['CHS', ''],
                'jec_payloads': ['AK4PFchs', 'AK4PF'],
                'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute'],
                'pu_jet_id': True,
                'qg_tagger': True,
                },
            }
        
    else:

        jetsCollections = {
            'AK4': {
                'algo': 'ak4',
                'pu_methods': ['CHS', ''],
                'jec_payloads': ['AK4PFchs', 'AK4PF'],
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
                appendJECToDB(process, jec_payload, jec_database_PF.replace("_PFCHS","").replace("_PF",""))

            jetToolbox(process, params['algo'], 'dummy', 'out', runOnMC=isMC, PUMethod = pu_method, JETCorrPayload = jec_payload, JETCorrLevels = jec_levels, addPUJetID = True)

            algo          = params['algo'].upper()
            jetCollection = '%sPFJets%s' % (params['algo'], pu_method)
            postfix       = '%sPF%s' % (algo, pu_method)
            if params['pu_jet_id']:

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
            if 'qg_tagger' in params and params['qg_tagger']:

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
    #from PhysicsTools.PatAlgos.tools.metTools import addMETCollection

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
                  jetCollectionPF    = "selectedPatJetsAK4PF", 
                  dRCleaning = 0.3, 
                  jetPtCut = jetPtCut, 
                  jetEtaCut = 5.,
                  etaCutForMetDiagnostic = etaCutForMetDiagnostic,
                  genJetCollection = "ak4GenJetsNoNu",
                  cleanGenJets = True,
                  applyZSelections = applyZSelections, 
                  applyWSelections = applyWSelections
                  )
 
 
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

    # Run
    if isMC:
        process.run = cms.EDAnalyzer('RunAnalyzer')
        process.jmfw_analyzers += process.run
    


    setattr(process, "PUPPET", 
            cms.EDAnalyzer('PUPPETAnalyzer',
                           isMC      = cms.bool(isMC),
                           srcJet    = cms.InputTag("selectedPatJetsAK4PFCleaned"),
                           srcJetPF  = cms.InputTag("selectedPatJetsAK4PFCleaned"),
                           srcVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                           srcZboson = cms.InputTag("mvaMET","ZtagBoson"),
                           srcLeptons = cms.InputTag("LeptonMerge"),
                           srcGenMet = cms.InputTag("slimmedMETs","","PAT"),
                           srcGenJets          = cms.InputTag("slimmedGenJets","","PAT"),
                           srcGenJetsCleaned   = cms.InputTag("selectedPatak4GenJetsNoNuCleaned"),
                           srcGenParticles     = cms.InputTag("prunedGenParticles","","PAT"),
                           srcGenEventInfo     = cms.InputTag("generator"),
                           srcPFMet            = cms.InputTag("slimmedMETs"),
                           srcPFCHSMet         = cms.InputTag("patPFMetCHS"),
                           srcPFPuppiMet       = cms.InputTag("slimmedMETsPuppi"),

                           srcRecoilPFPuppiMet       = cms.InputTag("mvaMET", "recoilslimmedMETsPuppi"),
                           srcRecoilPFMet      = cms.InputTag("mvaMET","recoilslimmedMETs"),
                           srcRecoilPFCHSMet   = cms.InputTag("mvaMET","recoilpatPFMetCHS"),
                           srcRecoilPFChargedPVNeutralPVPUJetID = cms.InputTag("mvaMET","recoilpatPFMetChargedPVNeutralPVPUJetID"),
                           srcRecoilPFChargedPUNeutralPUPUJetID = cms.InputTag("mvaMET","recoilpatPFMetChargedPUNeutralPUPUJetID"),
                           srcRecoilPFChargedPVNeutralPV        = cms.InputTag("mvaMET","recoilpatPFMetChargedPVNeutralPV"),

                           srcMVAMet     = cms.InputTag("mvaMET","mvaMET"),
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
