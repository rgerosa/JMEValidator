import FWCore.ParameterSet.Config as cms

def applyMuonID(process, label = "Tight", 
                src  = 'slimmedMuons', 
                type = 'tightID',                                                                                              
                iso_map_charged_hadron  = '',                                                                                                                                   
                iso_map_neutral_hadron = '',                                                                                                                                    
                iso_map_photon         = '',                                                                                                                                    
                rho     = 'fixedGridRhoFastjetAll',                                                                                                                      
                vertex  = 'offlineSlimmedPrimaryVertices',
                ptVal   = 5.,
                etaVal  = 2.4,
                typeIsoVal = 0
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


def applyElectronID(process, label = "Tight", 
                    src  = 'slimmedElectrons', 
                    type = 'tightID',                                                                                              
                    iso_map_charged_hadron = '',                                                                               
                    iso_map_neutral_hadron = '',                                                                                                                       
                    iso_map_photon         = '',                                                                                                                                
                    electron_id_map        = '',
                    rho     = 'fixedGridRhoFastjetAll',                                                                                                                      
                    vertex  = 'offlineSlimmedPrimaryVertices',
                    ptVal   = 5.,
                    etaVal  = 2.4,
                    typeIsoVal = 0
                    ):

  if "Tight" in label or "tight" in label :
    relativeIsolationCutVal = 0.09
  else:
    relativeIsolationCutVal = 0.15

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




def cleanJetsFromLeptons (process, 
                          label = "Cleaned",
                          jetCollection  = '',
                          muonCollection  = '',
                          electronCollection  = '',
                          tauCollection  = '',
                          jetPtCut       = 0.,
                          jetEtaCut      = 99,
                          dRCleaning     = 0.3
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
