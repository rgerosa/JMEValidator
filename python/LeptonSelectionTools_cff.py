import FWCore.ParameterSet.Config as cms

##### to apply muon ID
def applyMuonID(process, 
                label = "Tight",  ##add to the process module
                src   = 'slimmedMuons',  ## input collection
                type  = 'tightID',  ## string to identify the type of the muon ID                                                                                             
                iso_map_charged_hadron = '', ## isolation value map from charged hadrons                                                                          
                iso_map_neutral_hadron = '', ## isolation value map from neutral hadrons                                                                                
                iso_map_photon         = '', ## isolation value map from photons
                rho     = 'fixedGridRhoFastjetAll',                                                                                                                      
                vertex  = 'offlineSlimmedPrimaryVertices',
                ptVal   = 10., ## pt Cut 
                etaVal  = 2.4, ## eta Cut
                typeIsoVal = 0 ## inter number to identify the type of isolation to be used in the lepton ID
                ):

  if "Tight" in label or "tight" in label :
    relativeIsolationCutVal = 0.13
  else:
    relativeIsolationCutVal = 0.20

  setattr(process,src+label, cms.EDProducer("patMuonIDIsoSelector",
                                            src = cms.InputTag(src),
                                            rho = cms.InputTag(rho),
                                            vertex = cms.InputTag(vertex),
                                            charged_hadron_iso = cms.InputTag(iso_map_charged_hadron),
                                            neutral_hadron_iso = cms.InputTag(iso_map_neutral_hadron),
                                            photon_iso         = cms.InputTag(iso_map_photon),                                                       
                                            relativeIsolationCut = cms.double(relativeIsolationCutVal),
                                            typeID = cms.string(type),
                                            ptCut  = cms.double(ptVal),
                                            etaCut = cms.double(etaVal),
                                            typeIso = cms.int32(typeIsoVal)
                                            ))

##### to apply electron ID
def applyElectronID(process, 
                    label = "Tight", ##add to the process module 
                    src  = 'slimmedElectrons', ## input collection
                    type = 'tightID', ## string to identify the type of the electron ID                                                                            
                    iso_map_charged_hadron = '', ## isolation value map from charged hadrons                                                                              
                    iso_map_neutral_hadron = '', ## isolation value map from neutral hadrons                                                                         
                    iso_map_photon         = '', ## isolation value map from photons
                    electron_id_map        = '', ## value map for electron id 
                    rho     = 'fixedGridRhoFastjetAll',                                                                                                                      
                    vertex  = 'offlineSlimmedPrimaryVertices',
                    ptVal   = 10., ## pt Cut
                    etaVal  = 2.4, ## eta Cut
                    typeIsoVal = 0 ## inter number to identify the type of isolation to be used in the lepton ID 
                    ):

  if "Tight" in label or "tight" in label :
    relativeIsolationCutVal = 0.13
  else:
    relativeIsolationCutVal = 0.20

  setattr(process,src+label, cms.EDProducer("patElectronIDIsoSelector",
                                            src = cms.InputTag(src),
                                            rho = cms.InputTag(rho),
                                            vertex = cms.InputTag(vertex),
                                            charged_hadron_iso = cms.InputTag(iso_map_charged_hadron),
                                            neutral_hadron_iso = cms.InputTag(iso_map_neutral_hadron),
                                            photon_iso         = cms.InputTag(iso_map_photon),                                                       
                                            relativeIsolationCut = cms.double(relativeIsolationCutVal),
                                            electron_id = cms.InputTag(electron_id_map),
                                            ptCut    = cms.double(ptVal),
                                            etaCut   = cms.double(etaVal),
                                            typeIso  = cms.int32(typeIsoVal)
                                            ))


##### to apply tau ID + lepton cleaning
def    applyTauID( process, 
                   label = "Loose", 
                   src    = "slimmedTaus", ## input collection
                   ptCut  = 20., ## ptCut
                   etaCut = 2.4, ## eta Cut
                   muonCollection     = "", ## muon collection for cleaning
                   electronCollection = "", ## electron collection for cleaning
                   dRCleaning = 0.3 ## dR for cleaning
                   ):


  if label == "Loose" :

    ## apply tau loose ID
    setattr(process,src+label,cms.EDFilter("PATTauRefSelector",
                                           src = cms.InputTag(src),
                                           filter = cms.bool(False),
                                           cut  = cms.string('pt > %f & abs(eta) < %f & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") < 1.5 & tauID("againstMuonLoose3") > 0.5 & tauID("againstElectronLooseMVA5") > 0.5 '%(ptCut,etaCut))
                                           ));


  ## medium tau id
  elif label == "Medium" :  

    setattr(process,src+label, cms.EDFilter("PATTauRefSelector",
                                            src = cms.InputTag(src),
                                            filter = cms.bool(False),
                                            cut  = cms.string('pt > %f & abs(eta) < %f & tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") < 1.5 & tauID("againstMuonMedium3") > 0.5 & tauID("againstElectronMediumMVA5") > 0.5 '%(ptCut,etaCut)),
                                           ));
  ## tau tight ID
  elif label == "Tight" :  

    setattr(process,src+label, cms.EDFilter("PATTauRefSelector",
                                           src = cms.InputTag(src),
                                           filter = cms.bool(False),
                                           cut  = cms.string('pt > %f & abs(eta) < %f & tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") < 1.5 & tauID("againstMuonTight3") > 0.5 & tauID("againstElectronTightMVA5") > 0.5 '%(ptCut,etaCut))
                                            ));

  else :
    print "not known string for tau ID .. please check!! "
    exit(1);


  ## cleaning taus from other leptons
  tausNotOverlappingWithLeptons =  cms.EDProducer("PATTauCleaner",
                                                  src = cms.InputTag(src+label),
                                                  checkOverlaps = cms.PSet(),
                                                  finalCut = cms.string(''),
                                                  preselection = cms.string('')
                                                  )


  if muonCollection != '':
    setattr(tausNotOverlappingWithLeptons.checkOverlaps,"muons",cms.PSet( src = cms.InputTag(muonCollection),
                                                                          algorithm = cms.string("byDeltaR"),
                                                                          preselection        = cms.string(""),
                                                                          deltaR              = cms.double(dRCleaning),
                                                                          checkRecoComponents = cms.bool(False),
                                                                          pairCut             = cms.string(""),
                                                                          requireNoOverlaps   = cms.bool(True))) 
    if electronCollection != '':
      setattr(tausNotOverlappingWithLeptons.checkOverlaps,"electrons",cms.PSet( src = cms.InputTag(electronCollection),
                                                                                algorithm = cms.string("byDeltaR"),
                                                                                preselection        = cms.string(""),
                                                                                deltaR              = cms.double(dRCleaning),
                                                                                checkRecoComponents = cms.bool(False),
                                                                                pairCut             = cms.string(""),
                                                                                requireNoOverlaps   = cms.bool(True))) 
      
    ## copy in the process
    setattr(process,src+label+"Cleaned",tausNotOverlappingWithLeptons.clone());


## to clean hets from selected leptons
def cleanJetsFromLeptons (process, 
                          label = "Cleaned", ## label to be added to output collection
                          jetCollection       = '', ## input jet collection
                          muonCollection      = '', ## input muon collection
                          electronCollection  = '', ## input electron collection
                          tauCollection       = '', ## input tau collection
                          jetPtCut       = 15., ## pt threshold on jets
                          jetEtaCut      = 5.0, ## eta cut
                          dRCleaning     = 0.3  ## cleaning cone
                          ):

  from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import cleanPatJets 

  jetsNotOverlappingWithLeptons =  cms.EDProducer("PATJetCleaner",
                                                  src = cms.InputTag(jetCollection),
                                                  preselection = cms.string(('pt > %f && abs(eta) < %f')%(jetPtCut,jetEtaCut)),
                                                  checkOverlaps = cms.PSet(),
                                                  finalCut = cms.string('')
                                                  )

  if muonCollection != '':
    setattr(jetsNotOverlappingWithLeptons.checkOverlaps,"muons",cms.PSet( src = cms.InputTag(muonCollection),
                                                                          algorithm = cms.string("byDeltaR"),
                                                                          preselection        = cms.string(""),
                                                                          deltaR              = cms.double(dRCleaning),
                                                                          checkRecoComponents = cms.bool(False),
                                                                          pairCut             = cms.string(""),
                                                                          requireNoOverlaps   = cms.bool(True))) 
  if electronCollection != '':
    setattr(jetsNotOverlappingWithLeptons.checkOverlaps,"electrons",cms.PSet( src = cms.InputTag(electronCollection),
                                                                              algorithm = cms.string("byDeltaR"),
                                                                              preselection        = cms.string(""),
                                                                              deltaR              = cms.double(dRCleaning),
                                                                              checkRecoComponents = cms.bool(False),
                                                                              pairCut             = cms.string(""),
                                                                              requireNoOverlaps   = cms.bool(True))) 

  if tauCollection != '':  
    setattr(jetsNotOverlappingWithLeptons.checkOverlaps,"taus",cms.PSet( src = cms.InputTag(tauCollection),
                                                                         algorithm = cms.string("byDeltaR"),
                                                                         preselection        = cms.string(""),
                                                                         deltaR              = cms.double(dRCleaning),
                                                                         checkRecoComponents = cms.bool(False),
                                                                         pairCut             = cms.string(""),
                                                                         requireNoOverlaps   = cms.bool(True))) 

  setattr(process,jetCollection+label,jetsNotOverlappingWithLeptons);
