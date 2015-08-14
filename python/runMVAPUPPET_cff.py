import FWCore.ParameterSet.Config as cms
import sys

from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyMuonID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyElectronID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import applyTauID
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import cleanJetsFromLeptons
from JMEAnalysis.JMEValidator.LeptonSelectionTools_cff import cleanGenJetsFromGenLeptons
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
from RecoMET.METProducers.PFMET_cfi import pfMet
from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
from PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi import slimmedMETs
from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1

def runMVAPUPPET(process,
                 isMC,
                 srcMuons =  "slimmedMuons", ## inputMuonCollection
                 muonTypeID    = "Tight", ## type of muon ID to be applied                                                                                                   
                 iso_map_muons = [], ## isolation maps in case they have been re-computed (charged, neutral, photon)                                                         
                 typeIsoMuons  = "dBeta", ## isolation type to be used for muons                                                                                               
                 relativeIsoCutMuons = 0.12,
                 srcElectrons = "slimmedElectrons", 
                 electronTypeID= "Tight", 
                 electronID_map = 'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight',
                 iso_map_electrons = [], 
                 typeIsoElectrons = "rhoCorr",
                 relativeIsoCutEletrons = 0.12,
                 srcTaus = "slimmedTaus",
                 tauTypeID = "Loose",
                 doTauCleaning = True,
                 jetCollectionPuppi = "selectedPatJetsAK4PFPuppi" , 
                 jetCollectionPF    = "selectedPatJetsAK4PF" , 
                 dRCleaning= 0.3, 
                 jetPtCut = 0., 
                 jetEtaCut =5.,
                 genJetCollection = "",
                 cleanGenJets = False,
                 etaCutForMetDiagnostic = 10.,
                 applyTypeICorrection = True, 
                 ptThresholdForTypeIPuppi = 20, 
                 useJECFromLocalDB = True,                  
                 applyZSelections = True,
                 applyWSelections = False,
                 putRecoilInsideEvent = True
                 ):


    ## check if puppi is setup or not
    if not hasattr(process, 'puppi'):
        process.load('CommonTools.PileupAlgos.Puppi_cff')
        puppi.candName = cms.InputTag('packedPFCandidates')
        puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
        

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



    ##### apply W or Z selections
    if not applyWSelections :
        setattr(process,"ZdiMuon"+muonTypeID,cms.EDProducer("CandViewCombiner",
                                                            decay       = cms.string(srcMuons+muonTypeID+"@+ "+srcMuons+muonTypeID+"@-"),
                                                            checkCharge = cms.bool(True),
                                                            cut         = cms.string("mass > 70 && mass < 110 & charge=0"),
                                                            ))
        
        process.jmfw_analyzers += getattr(process,"ZdiMuon"+muonTypeID);

        setattr(process,"ZdiElectron"+electronTypeID,cms.EDProducer("CandViewCombiner",
                                                                    decay       = cms.string(srcElectrons+electronTypeID+"@+ "+srcElectrons+electronTypeID+"@-"),
                                                                    checkCharge = cms.bool(True),
                                                                    cut         = cms.string("mass > 70 && mass < 110 & charge=0"),
                                                                    ))

        process.jmfw_analyzers += getattr(process,"ZdiElectron"+electronTypeID);

        if tauTypeID != "":
                
            setattr(process,"ZdiTau"+tauTypeID,cms.EDProducer("CandViewCombiner",
                                                              decay       = cms.string(srcTaus+tauTypeID+"Cleaned@+ "+srcTaus+tauTypeID+"Cleaned@-"),
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
                                                          src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID,srcTaus+tauTypeID+"Cleaned")
                                                          ))

        else:
            setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                          src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID)
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
                                                          src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID,srcTaus+tauTypeID+"Cleaned")
                                                              ))

        else:
            setattr(process,"LeptonMerge", cms.EDProducer("CandViewMerger",
                                                              src = cms.VInputTag(srcMuons+muonTypeID,srcElectrons+electronTypeID)
                                                              ))

        process.jmfw_analyzers += getattr(process,"LeptonMerge");
            
        setattr(process,"LeptonMergeFilter",cms.EDFilter("PATCandViewCountFilter",
                                                         minNumber = cms.uint32(1),
                                                         maxNumber = cms.uint32(1),
                                                         src = cms.InputTag("LeptonMerge")
                                                         ))

        process.jmfw_analyzers += getattr(process,"LeptonMergeFilter");
            

    ## jet lepton cleaning
    setattr(getattr(process,jetCollectionPuppi),"cut",cms.string('pt > %f'%jetPtCut))
    setattr(getattr(process,jetCollectionPF),"cut",cms.string('pt > %f'%jetPtCut))

    cleanJetsFromLeptons(process,
                         label = "Cleaned",
                         jetCollection      = jetCollectionPuppi,
                         muonCollection     = srcMuons+muonTypeID,
                         electronCollection = srcElectrons+electronTypeID,
                         tauCollection      = srcTaus+tauTypeID+"Cleaned",
                         jetPtCut   = jetPtCut,
                         jetEtaCut  = jetEtaCut,
                         dRCleaning = dRCleaning)

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

    process.puppi.candName = cms.InputTag('packedPFCandidatesNoLepton')
    
    process.pupuppi = process.puppi.clone()
    process.pupuppi.invertPuppi = True
    

    ### produce the collection neutrals in and out jets passing PU jet id
    process.neutralInJets = cms.EDProducer("neutralCandidatePUIDJets",
                                           srcJets = cms.InputTag(jetCollectionPF+"Cleaned"),
                                           srcCandidates = cms.InputTag("packedPFCandidates"),
                                           neutralParticlesPVJetsLabel = cms.string("neutralPassingPUIDJets"),
                                           neutralParticlesPUJetsLabel = cms.string("neutralFailingPUIDJets"),
                                           neutraParticlesPVLabel = cms.string("neutralParticlesPV"),
                                           jetPUDIWP = cms.string("user"),
                                           jetPUIDMapLabel = cms.string("fullDiscriminant"))

    ### puppi raw met                                                                                                                                             
    process.pfCandidatesForPuppiMET = cms.EDFilter("CandPtrSelector",
                                                   src = cms.InputTag("puppi"),
                                                   cut = cms.string("abs(eta) < %f"%etaCutForMetDiagnostic))
    process.pfMetPuppi       = pfMet.clone()
    process.pfMetPuppi.src   = cms.InputTag("pfCandidatesForPuppiMET")
    process.pfMetPuppi.alias = cms.string('pfMetPuppi')
    addMETCollection(process, labelName='patPFMetPuppi', metSource='pfMetPuppi') # RAW puppi MET                                                                   
    process.patPFMetPuppi.addGenMET = False
    
    ## Type 1 corrections                                                                                                                                                       
    if not hasattr(process, 'ak4PFJetsPuppi'):
        print("WARNING: No AK4 puppi jets produced. Type 1 corrections for puppi MET are not available.")
    else:
        if applyTypeICorrection :

            from CommonTools.RecoAlgos.pfJetSelector_cfi import pfJetSelector

            process.ak4PFJetsPuppiForTypeI = pfJetSelector.clone(
                src = cms.InputTag( "ak4PFJetsPuppi" ),
                cut = cms.string( "abs(eta)<%f"%(etaCutForMetDiagnostic) )
                )

            if useJECFromLocalDB :
                
                process.ak4PuppiL1FastjetCorrector = process.ak4PFCHSL1FastjetCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
                process.ak4PuppiL2RelativeCorrector = process.ak4PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
                process.ak4PuppiL3AbsoluteCorrector = process.ak4PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
                process.ak4PuppiL1FastL2L3Corrector = process.ak4PFL1FastL2L3Corrector.clone(
                    correctors = cms.VInputTag("ak4PuppiL1FastjetCorrector", "ak4PuppiL2RelativeCorrector", "ak4PuppiL3AbsoluteCorrector")
                    )
                process.ak4PuppiResidualCorrector = process.ak4PFResidualCorrector.clone( algorithm = 'AK4PFPuppi' )
                process.ak4PuppiL1FastL2L3ResidualCorrector = process.ak4PFL1FastL2L3ResidualCorrector.clone( 
                    correctors = cms.VInputTag("ak4PuppiL1FastjetCorrector", "ak4PuppiL2RelativeCorrector", "ak4PuppiL3AbsoluteCorrector", "ak4PuppiResidualCorrector")
                    )

                if isMC:

                    process.corrPfMetType1Puppi = corrPfMetType1.clone(
                        src = 'ak4PFJetsPuppiForTypeI',
                        jetCorrLabel = 'ak4PuppiL1FastL2L3Corrector',
                        offsetCorrLabel = 'ak4PuppiL1FastjetCorrector',
                        type1JetPtThreshold = cms.double(ptThresholdForTypeIPuppi)
                        )

                else:
                    process.corrPfMetType1Puppi = corrPfMetType1.clone(
                        src = 'ak4PFJetsPuppiForTypeI',
                        jetCorrLabel = 'ak4PuppiL1FastL2L3ResidualCorrector',
                        offsetCorrLabel = 'ak4PuppiL1FastjetCorrector',
                        type1JetPtThreshold = cms.double(ptThresholdForTypeIPuppi)
                        )

            else :
                
                if isMC :
                    process.corrPfMetType1Puppi = corrPfMetType1.clone(
                        src = 'ak4PFJetsPuppiForTypeI',
                        jetCorrLabel = 'ak4PFCHSL1FastL2L3Corrector', #FIXME: Use PUPPI corrections when available?                                                             
                        offsetCorrLabel = 'ak4PFCHSL1FastjetCorrector'
                        )
                else:           
                    process.corrPfMetType1Puppi = corrPfMetType1.clone(
                        src = 'ak4PFJetsPuppiForTypeI',
                        jetCorrLabel = 'ak4PFCHSL1FastL2L3ResidualCorrector', #FIXME: Use PUPPI corrections when available?                                                 
                        offsetCorrLabel = 'ak4PFCHSL1FastjetCorrector'
                        )
                

            process.pfMetT1Puppi = pfMetT1.clone(
                src = 'pfMetPuppi',
                srcCorrections = [ cms.InputTag("corrPfMetType1Puppi","type1") ]
                )

            ## new PAT Met correction                                                                                                                                           
            addMETCollection(process, labelName='patMETPuppi', metSource='pfMetT1Puppi') # T1 puppi MET                                                                     
            process.patMETPuppi.addGenMET = False


    ### PUPPI                                                                                                                                                             
    process.slimmedMETsPuppi = slimmedMETs.clone()
    if hasattr(process, "patMETPuppi"):
        # Create MET from Type 1 PF collection                                                                                                                                  
        process.patMETPuppi.addGenMET = isMC
        process.slimmedMETsPuppi.src = cms.InputTag("patMETPuppi")
        process.slimmedMETsPuppi.rawUncertainties = cms.InputTag("patPFMetPuppi") # only central value                                                                          
    else:
        # Create MET from RAW PF collection                                                                                                                                     
        process.patPFMetPuppi.addGenMET = isMC
        process.slimmedMETsPuppi.src = cms.InputTag("patPFMetPuppi")
        process.slimmedMETsPuppi.rawUncertainties = cms.InputTag("patPFMetPuppi") # only central value                                                                          

    del process.slimmedMETsPuppi.type1Uncertainties # not available                                                                                                             
    del process.slimmedMETsPuppi.type1p2Uncertainties # not available                                                                                                          

    
    #### puppi charged particles == charged PV                                                                                                                             
    process.pfAllChargedParticlesPuppi = cms.EDFilter("CandPtrSelector",
                                                      src = cms.InputTag("puppi"),
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
                                                               src = cms.VInputTag("pfAllChargedParticlesPuppi",cms.InputTag("neutralInJets","neutralPassingPUIDJets"))
                                                               )

    process.pfChargedPVNeutralsPVPUJetID = cms.EDFilter("CandPtrSelector",
                                                       src = cms.InputTag("pfChargedPVNeutralsPVPUJetIDMerge"),
                                                       cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

    ## inverted puppi : neutrals PU and charged puppi                                                                                                                        
    process.pfPUPuppi = cms.EDFilter("CandPtrSelector",
                                     src = cms.InputTag("pupuppi"),
                                     cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))
                                     
    process.pfPUPuppiCharge = cms.EDFilter("CandPtrSelector",
                                           src = cms.InputTag("pupuppi"),
                                           cut = cms.string("pt > 0 && charge != 0 && abs(eta) < %f"%etaCutForMetDiagnostic))
    
    process.pfAllNeutralParticlesPuppiPU  = cms.EDFilter("CandPtrSelector",
                                                     src = cms.InputTag("pupuppi"),
                                                     cut = cms.string("charge=0 && pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))


    ### Charged PU + neutrals in jet failing PUID
    process.pfChargedPUNeutralsPUPUJetIDMerge = cms.EDProducer("CandViewMerger",
                                                               src = cms.VInputTag("pfChargedPU",cms.InputTag("neutralInJets","neutralFailingPUIDJets"))
                                                               )

    process.pfChargedPUNeutralsPUPUJetID = cms.EDFilter("CandPtrSelector",
                                                        src = cms.InputTag("pfChargedPUNeutralsPUPUJetIDMerge"),
                                                        cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

    
    ## chargePV + all neutrals - PU neutrals
    process.pfChargedPVNeutralsPVMerge = cms.EDProducer("CandViewMerger",
                                                        src = cms.VInputTag(cms.InputTag("neutralInJets","neutralParticlesPV"),"pfAllChargedParticlesPuppi"))


    process.pfChargedPVNeutralsPV = cms.EDFilter("CandPtrSelector",
                                                 src = cms.InputTag("pfChargedPVNeutralsPVMerge"),
                                                 cut = cms.string("pt > 0 && abs(eta) < %f"%etaCutForMetDiagnostic))

    #### Missing energies                                                                                                                                                    
    process.pfMetPuppiChargedPV       = pfMet.clone()
    process.pfMetPuppiChargedPV.src   = cms.InputTag("pfAllChargedParticlesPuppi") ## packed candidates without fromPV < 1                                                   
    process.pfMetPuppiChargedPV.alias = cms.string('pfMetPuppiChargedPV')
    addMETCollection(process, labelName='patPFMetPuppiChargedPV', metSource='pfMetPuppiChargedPV') # Convert the CHS PFMet in PAT MET                                        
    process.patPFMetPuppiChargedPV.addGenMET = isMC

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
    process.patPFMetChargedPU.addGenMET = isMC
    
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
    process.patPFMetPuppiNeutralPV.addGenMET = isMC
    
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
    process.patPFMetPuppiNeutralPU.addGenMET = isMC
    
    process.slimmedMETsPuppiNeutralPU = slimmedMETs.clone()
    process.slimmedMETsPuppiNeutralPU.src = cms.InputTag("patPFMetPuppiNeutralPU")
    process.slimmedMETsPuppiNeutralPU.rawUncertainties = cms.InputTag("patPFMetPuppiNeutralPU") # only central value                                                         
    
    del process.slimmedMETsPuppiNeutralPU.type1Uncertainties # not available                                                                                                
    del process.slimmedMETsPuppiNeutralPU.type1p2Uncertainties # not available                                                                                               
    

    ### last inputs for standard MVA met
    process.pfMetChargedPVNeutralPVPUJetID       = pfMet.clone()
    process.pfMetChargedPVNeutralPVPUJetID.src   = cms.InputTag("pfChargedPVNeutralsPVPUJetID")
    process.pfMetChargedPVNeutralPVPUJetID.alias = cms.string('pfMetChargedPVNeutralPVPUJetID')
    addMETCollection(process, labelName='patPFMetChargedPVNeutralPVPUJetID', metSource='pfMetChargedPVNeutralPVPUJetID')
    process.patPFMetChargedPVNeutralPVPUJetID.addGenMET = isMC
    
    process.slimmedMETsChargedPVNeutralPVPUJetID = slimmedMETs.clone()
    process.slimmedMETsChargedPVNeutralPVPUJetID.src = cms.InputTag("patPFMetChargedPVNeutralPVPUJetID")
    process.slimmedMETsChargedPVNeutralPVPUJetID.rawUncertainties = cms.InputTag("patPFMetChargedPVNeutralPVPUJetID") # only central value                                   
    
    del process.slimmedMETsChargedPVNeutralPVPUJetID.type1Uncertainties # not available                                                                                        
    del process.slimmedMETsChargedPVNeutralPVPUJetID.type1p2Uncertainties # not available                                                                                       

    ### last inputs for standard MVA met
    process.pfMetChargedPUNeutralPUPUJetID       = pfMet.clone()
    process.pfMetChargedPUNeutralPUPUJetID.src   = cms.InputTag("pfChargedPUNeutralsPUPUJetID")
    process.pfMetChargedPUNeutralPUPUJetID.alias = cms.string('pfMetChargedPUNeutralPUPUJetID')
    addMETCollection(process, labelName='patPFMetChargedPUNeutralPUPUJetID', metSource='pfMetChargedPUNeutralPUPUJetID')
    process.patPFMetChargedPUNeutralPUPUJetID.addGenMET = isMC
    
    process.slimmedMETsChargedPUNeutralPUPUJetID = slimmedMETs.clone()
    process.slimmedMETsChargedPUNeutralPUPUJetID.src = cms.InputTag("patPFMetChargedPUNeutralPUPUJetID")
    process.slimmedMETsChargedPUNeutralPUPUJetID.rawUncertainties = cms.InputTag("patPFMetChargedPUNeutralPUPUJetID") # only central value                                   
    
    del process.slimmedMETsChargedPUNeutralPUPUJetID.type1Uncertainties # not available                                                                                        
    del process.slimmedMETsChargedPUNeutralPUPUJetID.type1p2Uncertainties # not available                                                                                       

    ### last inputs for standard MVA met
    process.pfMetChargedPVNeutralPV       = pfMet.clone()
    process.pfMetChargedPVNeutralPV.src   = cms.InputTag("pfChargedPVNeutralsPV")
    process.pfMetChargedPVNeutralPV.alias = cms.string('pfMetChargedPVNeutralPV')
    addMETCollection(process, labelName='patPFMetChargedPVNeutralPV', metSource='pfMetChargedPVNeutralPV')
    process.patPFMetChargedPVNeutralPV.addGenMET = isMC
    
    process.slimmedMETsChargedPVNeutralPV = slimmedMETs.clone()
    process.slimmedMETsChargedPVNeutralPV.src = cms.InputTag("patPFMetChargedPVNeutralPV")
    process.slimmedMETsChargedPVNeutralPV.rawUncertainties = cms.InputTag("patPFMetChargedPVNeutralPV") # only central value                                   
    
    del process.slimmedMETsChargedPVNeutralPV.type1Uncertainties # not available                                                                                        
    del process.slimmedMETsChargedPVNeutralPV.type1p2Uncertainties # not available                                                                                              

    ### MVA PUPPET
    setattr(process,"mvaPUPPET", cms.EDProducer("mvaPUPPET",
                                                referenceMET = cms.InputTag("slimmedMETsPuppi"),
                                                srcMETs      = cms.VInputTag(cms.InputTag("slimmedMETs","","JRA"),
                                                                             cms.InputTag("slimmedMETsCHS"),
                                                                             cms.InputTag("slimmedMETsPuppi","","JRA"),
                                                                             cms.InputTag("slimmedMETsPuppiChargedPU"),
                                                                             cms.InputTag("slimmedMETsPuppiChargedPV"),
                                                                             cms.InputTag("slimmedMETsPuppiNeutralPV"),
                                                                             cms.InputTag("slimmedMETsPuppiNeutralPU"),
                                                                             cms.InputTag("slimmedMETsChargedPVNeutralPVPUJetID"),
                                                                             cms.InputTag("slimmedMETsChargedPUNeutralPUPUJetID"),
                                                                             cms.InputTag("slimmedMETsChargedPVNeutralPV")),
                                                inputMETFlags = cms.vint32(1,1,0,0,0,0,0,1,0,1),
                                                srcJets        = cms.InputTag(jetCollectionPuppi+"Cleaned"),
                                                srcVertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                srcTaus        = cms.InputTag(srcTaus+tauTypeID+"Cleaned"),
                                                srcMuons       = cms.InputTag(srcMuons+muonTypeID),
                                                srcPuppiWeights     = cms.InputTag("puppi"),
                                                inputFileNames = cms.PSet(
                                    PhiCorrectionWeightFile = cms.FileInPath('JMEAnalysis/JMEValidator/data/PhiCorrection_PUPPI.root'),                                            
                                    RecoilCorrectionWeightFile  = cms.FileInPath('JMEAnalysis/JMEValidator/data/RecoilCorrection_PUPPI.root')                                      
                    ),
                                                srcLeptons  = cms.VInputTag("LeptonMerge"),
                                                mvaMETLabel = cms.string("mvaMET"),
                                                ZbosonLabel = cms.string("ZtagBoson"),
                                                produceRecoils = cms.bool(putRecoilInsideEvent)
                                                ))
