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

    #### Input definitions like in classic MVA MET
    #### tracks from PV
    process.pfChargedPV = cms.EDFilter("CandPtrSelector",
                                       cut = cms.string('fromPV && charge !=0'),
                                       src = cms.InputTag("packedPFCandidates")
                                       )
    #### tracks not from PV
    process.pfChargedPU = cms.EDFilter("CandPtrSelector",
                                       cut = cms.string('!fromPV && charge !=0'),
                                       src = cms.InputTag("packedPFCandidates")
                                       )
    #### Neutrals in Jets passing PU Jet ID
    #### and Neutrals in Jets not passing PU Jet ID
    ### TODO: unclustered Neutrals
    process.neutralInJets = cms.EDProducer("neutralCandidatePUIDJets",
                                           srcJets = cms.InputTag(jetCollectionPF+"Cleaned"),
                                           srcCandidates = cms.InputTag("packedPFCandidates"),
                                           neutralParticlesPVJetsLabel = cms.string("neutralPassingPUIDJets"),
                                           neutralParticlesPUJetsLabel = cms.string("neutralFailingPUIDJets"),
                                           neutralParticlesUnclustered = cms.string("neutralParticlesUnclustered"),
                                           jetPUDIWP = cms.string("user"),
                                           jetPUIDMapLabel = cms.string("fullDiscriminant"))
    

    #### Merge collections to produce corresponding METs
    #### PF MET
    #process.pfMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag("pfAllChargedParticles",cms.InputTag("packedPFCandidates"))
    process.pfMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(cms.InputTag("packedPFCandidates")))
    #### Track MET
    process.pfTrackMETCands = process.pfChargedPV.clone() #cms.EDProducer("CandViewMerger", src = cms.VInputTag("pfAllChargedParticles",cms.InputTag("pfChargedPV"))
    ## No-PU MET
    process.pfNoPUMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag("pfChargedPV",cms.InputTag("neutralInJets","neutralPassingPUIDJets")))
    ## PU corrected MET
    process.pfPUCorrectedMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag("pfChargedPV",cms.InputTag("neutralInJets", "neutralPassingPUIDJets"), cms.InputTag("neutralInJets","neutralParticlesUnclustered")))
    ## PU MET
    process.pfPUMETCands = cms.EDProducer("CandViewMerger", src = cms.VInputTag(cms.VInputTag("pfChargedPU",cms.InputTag("neutralInJets","neutralFailingPUIDJets"))))
                                                              
    from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
    patMETsForMVA = patMETs.clone()
    patMETsForMVA.computeMETSignificance = cms.bool(True)
    patMETsForMVA.addGenMET = cms.bool(False)
    patMETsForMVA.srcJets = cms.InputTag("selectedPatJetsAK4PF")
    patMETsForMVA.srcLeptons = cms.InputTag("selectedPatJetsAK4PF")
    setattr(patMETsForMVA,"srcLeptons", cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID,srcTaus+tauTypeID+"Cleaned"))

    for met in ["pfMET", "pfTrackMET", "pfNoPUMET", "pfPUCorrectedMET", "pfPUMET"]:
        setattr(process, met, pfMet.clone())
        setattr(getattr(process, met), "src", cms.InputTag(met+"Cands"))
        setattr(getattr(process, met), "alias", cms.string(met))
        setattr(process, "pat"+met, patMETsForMVA.clone())
        setattr(getattr(process, "pat"+met), "metSource", cms.InputTag(met))

    ### MVA PUPPET
    setattr(process,"mvaMET", cms.EDProducer("mvaPUPPET",                                                
                                                referenceMET = cms.InputTag("slimmedMETs"),
                                                debug = cms.bool(True),
                                                srcMETs      = cms.VInputTag(cms.InputTag("patpfMET"),
                                                                             cms.InputTag("patpfTrackMET"),
                                                                             cms.InputTag("patpfNoPUMET"),
                                                                             cms.InputTag("patpfPUCorrectedMET"),
                                                                             cms.InputTag("patpfPUMET"),
                                                                             cms.InputTag("slimmedMETsPuppi"),
                                                                             cms.InputTag("slimmedMETs")
 ),
                                                inputMETFlags = cms.vint32(1,1,1,1,1,0,1,1),
                                                srcJets        = cms.InputTag(jetCollectionPF+"Cleaned"),
                                                srcVertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                srcTaus        = cms.InputTag(srcTaus+tauTypeID+"Cleaned"),
                                                srcMuons       = cms.InputTag(srcMuons+muonTypeID),
                                                inputFileNames = cms.PSet(
                #PhiCorrectionWeightFile = cms.FileInPath('JMEAnalysis/JMEValidator/data/PhiCorrection_PUPPI.root'),                                    
                #RecoilCorrectionWeightFile  = cms.FileInPath('JMEAnalysis/JMEValidator/data/RecoilCorrection_PUPPI.root')                           
                                    ),
                                                srcLeptons  = cms.VInputTag("LeptonMerge"),
                                                mvaMETLabel = cms.string("mvaMET"),
                                                ZbosonLabel = cms.string("ZtagBoson"),
                                                produceRecoils = cms.bool(putRecoilInsideEvent)
                                                ))
