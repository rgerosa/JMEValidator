import FWCore.ParameterSet.Config as cms
import sys

from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyMuonID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyElectronID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyTauID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import cleanJetsFromLeptons
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import cleanGenJetsFromGenLeptons
#from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
from RecoMET.METProducers.PFMET_cfi import pfMet
from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

def runMVAPUPPET(process,
                 processName,
                 isMC,
                 srcMuons =  "slimmedMuons", ## inputMuonCollection
                 muonTypeID    = "Tight", ## type of muon ID to be applied                                                                                                   
                 iso_map_muons = [], ## isolation maps in case they have been re-computed (charged, neutral, photon)                                                         
                 typeIsoMuons  = "dBeta", ## isolation type to be used for muons                                                                                               
                 relativeIsoCutMuons = 0.12,
                 srcElectrons = "slimmedElectrons", 
                 electronTypeID= "Tight", 
                 electronID_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium',
                 electronID_map_loose = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose',
                 iso_map_electrons = [], 
                 typeIsoElectrons = "rhoCorr",
                 relativeIsoCutEletrons = 0.12,
                 srcTaus = "slimmedTaus",
                 tauTypeID = "Loose",
                 doTauCleaning = True,
                 jetCollectionPF    = "selectedPatJetsAK4PF" , 
                 dRCleaning= 0.3, 
                 jetPtCut = 0., 
                 jetEtaCut =5.,
                 genJetCollection = "",
                 cleanGenJets = False,
                 etaCutForMetDiagnostic = 10.,
                 applyTypeICorrection = True, 
                 applyZSelections = True,
                 applyWSelections = False,
                 putRecoilInsideEvent = True
                 ):

    relativeIsoCutMuonsLoose = relativeIsoCutMuons+0.05;
    relativeIsoCutEletronsLoose = relativeIsoCutEletrons+0.05;    

    ### run Muon ID
    if len(iso_map_muons) < 3 :
        
        applyMuonID(process, 
                    src   = srcMuons,
                    label = muonTypeID, 
                    iso_map_charged_hadron  = '',
                    iso_map_neutral_hadron  = '',
                    iso_map_photon          = '',
                    typeIso                 = typeIsoMuons,
                    relativeIsolationCutVal = relativeIsoCutMuons
                    )
    else:

        applyMuonID(process, 
                    src   = srcMuons, 
                    label = muonTypeID, 
                    iso_map_charged_hadron  = iso_map_muons[0],
                    iso_map_neutral_hadron  = iso_map_muons[1],
                    iso_map_photon          = iso_map_muons[2],
                    rho = 'fixedGridRhoFastjetAll',
                    typeIso                 = typeIsoMuons,
                    relativeIsolationCutVal = relativeIsoCutMuons
                )


    ### run Electron ID
    if len(iso_map_electrons) < 3 :
        applyElectronID(process, 
                        label = electronTypeID, 
                        src   = srcElectrons,
                        iso_map_charged_hadron  = '',
                        iso_map_neutral_hadron  = '',
                        iso_map_photon          = '',
                        typeIso = typeIsoElectrons,
                        electron_id_map = electronID_map,
                        relativeIsolationCutVal = relativeIsoCutEletrons
                        )
    else:
        applyElectronID(process, 
                        label = electronTypeID, 
                        src   = srcElectrons,
                        iso_map_charged_hadron  = iso_map_electrons[0],
                        iso_map_neutral_hadron  = iso_map_electrons[1],
                        iso_map_photon          = iso_map_electrons[2],
                        typeIso = typeIsoElectrons,
                        electron_id_map = electronID_map,
                        relativeIsolationCutVal = relativeIsoCutEletrons
                        )



    ### run tau ID                                        
    if doTauCleaning :
        applyTauID( process, 
                    label = tauTypeID, 
                    src = srcTaus,
                    muonCollection     = srcMuons+muonTypeID,
                    electronCollection = srcElectrons+electronTypeID)

    else:
        applyTauID( process, label = tauTypeID, 
                    src = srcTaus,
                    muonCollection     = "",
                    electronCollection = "")


    ############
    ############

    if applyZSelections and applyWSelections:
        sys.exit("Z and W selections cannot be applied at the same time --> exit");

    if not applyZSelections and not applyWSelections:
        sys.exit("Z and W selections cannot be false at the same time --> exit");

    ##### apply Z selection ######
    if applyZSelections :
        
        setattr(process,"ZdiMuon"+muonTypeID,cms.EDProducer("CandViewCombiner",
                                                            decay       = cms.string(srcMuons+muonTypeID+"@+ "+srcMuons+muonTypeID+"@-"),
                                                            checkCharge = cms.bool(True),
                                                            cut         = cms.string("mass > 70 && mass < 110 & charge=0")))
        
        setattr(process,"ZdiElectron"+electronTypeID,cms.EDProducer("CandViewCombiner",
                                                                    decay       = cms.string(srcElectrons+electronTypeID+"@+ "+srcElectrons+electronTypeID+"@-"),
                                                                    checkCharge = cms.bool(True),
                                                                    cut         = cms.string("mass > 70 && mass < 110 & charge=0"),
                                                                    ))

        process.jmfw_analyzers += getattr(process,"ZdiMuon"+muonTypeID);
        process.jmfw_analyzers += getattr(process,"ZdiElectron"+electronTypeID);

        if tauTypeID != "":
                
            setattr(process,"ZdiTau"+tauTypeID,cms.EDProducer("CandViewCombiner",
                                                              decay       = cms.string(srcTaus+tauTypeID+"Cleaned@+ "+srcTaus+tauTypeID+"Cleaned@-"),
                                                              checkCharge = cms.bool(True),
                                                              cut         = cms.string("mass > 70 && mass < 110 & charge=0")))

            process.jmfw_analyzers += getattr(process,"ZdiTau"+tauTypeID);

            ## merge all the Z canddates and ask for only one candidate per event                                                                                            
            setattr(process,"ZdiLepton", cms.EDProducer("CandViewMerger",
                                                        src = cms.VInputTag("ZdiMuon"+muonTypeID,"ZdiElectron"+electronTypeID,"ZdiTau"+tauTypeID)))
 
        else:

                ## merge all the Z canddates and ask for only one candidate per event                                                                                           
                setattr(process,"ZdiLepton", cms.EDProducer("CandViewMerger",
                                                            src = cms.VInputTag("ZdiMuon"+muonTypeID,"ZdiElectron"+electronTypeID)))

        process.jmfw_analyzers += getattr(process,"ZdiLepton");

        
        setattr(process,"ZdiLeptonFilter",cms.EDFilter("PATCandViewCountFilter",
                                                       minNumber = cms.uint32(1),
                                                       maxNumber = cms.uint32(1),
                                                       src = cms.InputTag("ZdiLepton")))
        
        process.jmfw_analyzers += getattr(process,"ZdiLeptonFilter");

        ### count the number of leptons        
        if tauTypeID != "" :
            setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                          src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID,srcTaus+tauTypeID+"Cleaned")))

        else:
            setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                          src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID)))
            
        process.jmfw_analyzers += getattr(process,"LeptonMerge");
            
        setattr(process,"LeptonMergeFilter",cms.EDFilter("PATCandViewCountFilter",
                                                         minNumber = cms.uint32(2),
                                                         maxNumber = cms.uint32(2),
                                                         src = cms.InputTag("LeptonMerge")
                                                         ))
        process.jmfw_analyzers += getattr(process,"LeptonMergeFilter");

    if applyWSelections:

        if tauTypeID != "" :
            setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                          src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID,srcTaus+tauTypeID+"Cleaned")))

        else:
            setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                          src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID)))

        process.jmfw_analyzers += getattr(process,"LeptonMerge");
            
        setattr(process,"LeptonMergeFilter",cms.EDFilter("PATCandViewCountFilter",
                                                         minNumber = cms.uint32(1),
                                                         maxNumber = cms.uint32(1),
                                                         src = cms.InputTag("LeptonMerge")
                                                         ))

        process.jmfw_analyzers += getattr(process,"LeptonMergeFilter");
            

        #### run loose muon selection
        if len(iso_map_muons) < 3 :
        
            applyMuonID(process, 
                        src   = srcMuons,
                        label = "Loose", 
                        iso_map_charged_hadron  = '',
                        iso_map_neutral_hadron  = '',
                        iso_map_photon          = '',
                        typeIso                 = typeIsoMuons,
                        relativeIsolationCutVal = relativeIsoCutMuonsLoose
                    )
        else:

            applyMuonID(process, 
                        src   = srcMuons, 
                        label = "Loose", 
                        iso_map_charged_hadron  = iso_map_muons[0],
                        iso_map_neutral_hadron  = iso_map_muons[1],
                        iso_map_photon          = iso_map_muons[2],
                        rho = 'fixedGridRhoFastjetAll',
                        typeIso                 = typeIsoMuons,
                        relativeIsolationCutVal = relativeIsoCutMuonsLoose
                        )


        ### run Electron ID
        if len(iso_map_electrons) < 3 :
            applyElectronID(process, 
                            label = "Loose", 
                            src   = srcElectrons,
                            iso_map_charged_hadron  = '',
                            iso_map_neutral_hadron  = '',
                            iso_map_photon          = '',
                            typeIso = typeIsoElectrons,
                            electron_id_map = electronID_map_loose,
                            relativeIsolationCutVal = relativeIsoCutEletronsLoose
                            )
        else:
            applyElectronID(process, 
                            label = "Loose", 
                            src   = srcElectrons,
                            iso_map_charged_hadron  = iso_map_electrons[0],
                            iso_map_neutral_hadron  = iso_map_electrons[1],
                            iso_map_photon          = iso_map_electrons[2],
                            typeIso = typeIsoElectrons,
                            electron_id_map = electronID_map_loose,
                            relativeIsolationCutVal = relativeIsoCutEletronsLoose
                            )



        #### loose lepton veto
        if tauTypeID != "" :
            setattr(process,"LeptonMergeLoose", cms.EDProducer("CandViewMerger",
                                                               src = cms.VInputTag(srcMuons+"Loose",srcElectrons+"Loose",srcTaus+tauTypeID+"Cleaned")))

        else:
            setattr(process,"LeptonMergeLoose", cms.EDProducer("CandViewMerger",
                                                               src = cms.VInputTag(srcMuons+"Loose",srcElectrons+"Loose")))

        process.jmfw_analyzers += getattr(process,"LeptonMergeLoose");
            
        setattr(process,"LeptonMergeFilterLoose",cms.EDFilter("PATCandViewCountFilter",
                                                              minNumber = cms.uint32(1),
                                                              maxNumber = cms.uint32(1),
                                                              src = cms.InputTag("LeptonMerge")
                                                              ))

        process.jmfw_analyzers += getattr(process,"LeptonMergeFilterLoose");
        


    ## jet lepton cleaning
    setattr(getattr(process,jetCollectionPF),"cut",cms.string('pt > %f'%jetPtCut))

    cleanJetsFromLeptons(process,
                         label = "Cleaned",
                         jetCollection      = jetCollectionPF,
                         muonCollection     = srcMuons+muonTypeID,
                         electronCollection = srcElectrons+electronTypeID,
                         tauCollection      = srcTaus+tauTypeID+"Cleaned",
                         jetPtCut   = jetPtCut,
                         jetEtaCut  = jetEtaCut,
                         dRCleaning = dRCleaning)



    ## clean also the related Gen Jet Collection                                                                                                                        
    if isMC and cleanGenJets:
        process.packedGenLeptons = cms.EDFilter("CandPtrSelector",
                                                cut = cms.string('(abs(pdgId) = 11 || abs(pdgId) = 13) && pt > 10'),
                                                src = cms.InputTag("packedGenParticles")
                                                )
 

        setattr(process,"pat"+genJetCollection, cms.EDProducer("PATJetProducer",
                                                            jetSource = cms.InputTag(genJetCollection),
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

        setattr(process,"selectedPat"+genJetCollection,cms.EDFilter("PATJetSelector",
                                                                    cut = cms.string('pt > %f'%jetPtCut),
                                                                    src = cms.InputTag("pat"+genJetCollection)
                                                                    ))
        
        cleanGenJetsFromGenLeptons (process,
                                    jetCollection       = "selectedPat"+genJetCollection,
                                    genLeptonCollection = "packedGenLeptons",
                                    jetPtCut         = jetPtCut,
                                    jetEtaCut        = jetEtaCut,
                                    dRCleaning       = dRCleaning)

    ### puppi setup of PUPPI MET
    process.packedPFCandidatesNoLepton = cms.EDProducer("packedCandidateFilterParticles",
                                                        src      = cms.InputTag("packedPFCandidates"),
                                                        srcMuons = cms.InputTag(srcMuons+muonTypeID),
                                                        srcElectrons = cms.InputTag(srcElectrons+electronTypeID),
                                                        srcTaus = cms.InputTag(""))

    

    ### produce the collection neutrals in and out jets passing PU jet id
    process.neutralInJets = cms.EDProducer("neutralCandidatePUIDJets",
                                           srcJets = cms.InputTag(jetCollectionPF+"Cleaned"),
                                           srcCandidates = cms.InputTag("packedPFCandidates"),
                                           neutralParticlesPVJetsLabel = cms.string("neutralPassingPUIDJets"),
                                           neutralParticlesPUJetsLabel = cms.string("neutralFailingPUIDJets"),
                                           neutraParticlesPVLabel = cms.string("neutralParticlesPV"),
                                           jetPUDIWP = cms.string("user"),
                                           jetPUIDMapLabel = cms.string("fullDiscriminant"))

    #### puppi charged particles == charged PV                                                                                                                             
    process.pfAllChargedParticles = cms.EDFilter("CandPtrSelector",
                                                      src = cms.InputTag("packedPFCandidates"),
                                                      cut = cms.string("charge !=0 && pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

    #### charged particles PU                                                                                                                             
    process.pfChargedPU = cms.EDFilter("CandPtrSelector",
                                       cut = cms.string('!fromPV && abs(eta) < %f'%etaCutForMetDiagnostic),
                                       src = cms.InputTag("packedPFCandidates")
                                       )


    ### neutrals puppi candidate
    process.pfAllNeutralParticlesPuppi  = cms.EDFilter("CandPtrSelector",
                                                       src = cms.InputTag("puppi"),
                                                       cut = cms.string("charge == 0 && pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

    ### particles charged PV + neutrals in jet passing PUJET 
    process.pfChargedPVNeutralsPVPUJetIDMerge = cms.EDProducer("CandViewMerger",
                                                               src = cms.VInputTag("pfAllChargedParticles",cms.InputTag("neutralInJets","neutralPassingPUIDJets"))
                                                               )

    process.pfChargedPVNeutralsPVPUJetID = cms.EDFilter("CandPtrSelector",
                                                       src = cms.InputTag("pfChargedPVNeutralsPVPUJetIDMerge"),
                                                       cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

    process.pfChargedPUNeutralsPUPUJetIDMerge = cms.EDProducer("CandViewMerger",
                                                               src = cms.VInputTag("pfChargedPU",cms.InputTag("neutralInJets","neutralFailingPUIDJets"))
                                                               )

    process.pfChargedPUNeutralsPUPUJetID = cms.EDFilter("CandPtrSelector",
                                                        src = cms.InputTag("pfChargedPUNeutralsPUPUJetIDMerge"),
                                                        cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

    ## chargePV + all neutrals - PU neutrals
    process.pfChargedPVNeutralsPVMerge = cms.EDProducer("CandViewMerger",
                                                        src = cms.VInputTag(cms.InputTag("neutralInJets","neutralParticlesPV"),"pfAllChargedParticles"))


    process.pfChargedPVNeutralsPV = cms.EDFilter("CandPtrSelector",
                                                 src = cms.InputTag("pfChargedPVNeutralsPVMerge"),
                                                 cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))


    ### last inputs for standard MVA met
    process.pfMetChargedPVNeutralPVPUJetID       = pfMet.clone()
    process.pfMetChargedPVNeutralPVPUJetID.src   = cms.InputTag("pfChargedPVNeutralsPVPUJetID")
    process.pfMetChargedPVNeutralPVPUJetID.alias = cms.string('pfMetChargedPVNeutralPVPUJetID')
    from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
    patMETsForMVA = patMETs.clone()
    patMETsForMVA.computeMETSignificance = cms.bool(True)
    patMETsForMVA.addGenMET = cms.bool(False)
    patMETsForMVA.srcJets = cms.InputTag("selectedPatJetsAK4PF")
    patMETsForMVA.srcLeptons = cms.InputTag("selectedPatJetsAK4PF")
    setattr(patMETsForMVA,"srcLeptons", cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID,srcTaus+tauTypeID+"Cleaned"))
    
    process.patPFMetChargedPVNeutralPVPUJetID = patMETsForMVA.clone()
    process.patPFMetChargedPVNeutralPVPUJetID.metSource = cms.InputTag("pfMetChargedPVNeutralPVPUJetID")

    ### last inputs for standard MVA met
    process.pfMetChargedPUNeutralPUPUJetID       = pfMet.clone()
    process.pfMetChargedPUNeutralPUPUJetID.src   = cms.InputTag("pfChargedPUNeutralsPUPUJetID")
    process.pfMetChargedPUNeutralPUPUJetID.alias = cms.string('pfMetChargedPUNeutralPUPUJetID')
    process.patPFMetChargedPUNeutralPUPUJetID = patMETsForMVA.clone()
    process.patPFMetChargedPUNeutralPUPUJetID.metSource = cms.InputTag("pfMetChargedPUNeutralPUPUJetID")
    
    ### last inputs for standard MVA met
    process.pfMetChargedPVNeutralPV       = pfMet.clone()
    process.pfMetChargedPVNeutralPV.src   = cms.InputTag("pfChargedPVNeutralsPV")
    process.pfMetChargedPVNeutralPV.alias = cms.string('pfMetChargedPVNeutralPV')
    process.patPFMetChargedPVNeutralPV = patMETsForMVA.clone()
    process.patPFMetChargedPVNeutralPV.metSource = cms.InputTag("pfMetChargedPVNeutralPV")

    process.patPFMetCHS = patMETsForMVA.clone()
    process.patPFMetCHS.metSource = cms.InputTag("pfMetCHS")
    
    ### MVA PUPPET
    setattr(process,"mvaMET", cms.EDProducer("mvaPUPPET",                                                
                                                referenceMET = cms.InputTag("slimmedMETs"),
                                                debug = cms.bool(True),
                                                srcMETs      = cms.VInputTag(cms.InputTag("slimmedMETs"),
                                                                             cms.InputTag("patPFMetCHS"),
                                                                             cms.InputTag("patPFMetChargedPVNeutralPVPUJetID"),
                                                                             cms.InputTag("patPFMetChargedPUNeutralPUPUJetID"),
                                                                             cms.InputTag("patPFMetChargedPVNeutralPV"),
                                                                             cms.InputTag("slimmedMETsPuppi")
 ),
                                                inputMETFlags = cms.vint32(1,1,1,0,1,0,1),
                                                srcJets        = cms.InputTag(jetCollectionPF+"Cleaned"),
                                                srcVertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                srcTaus        = cms.InputTag(srcTaus+tauTypeID+"Cleaned"),
                                                srcMuons       = cms.InputTag(srcMuons+muonTypeID),
                                                srcPuppiWeights     = cms.InputTag("puppi"),
                                                inputFileNames = cms.PSet(
                #PhiCorrectionWeightFile = cms.FileInPath('JMEAnalysis/JMEValidator/data/PhiCorrection_PUPPI.root'),                                    
                #RecoilCorrectionWeightFile  = cms.FileInPath('JMEAnalysis/JMEValidator/data/RecoilCorrection_PUPPI.root')                           
                                    ),
                                                srcLeptons  = cms.VInputTag("LeptonMerge"),
                                                mvaMETLabel = cms.string("mvaMET"),
                                                ZbosonLabel = cms.string("ZtagBoson"),
                                                produceRecoils = cms.bool(putRecoilInsideEvent)
                                                ))
