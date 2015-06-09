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


#applyElectronID(process, label = "Veto",  src  = electronCollection,
#                identification_map     = '',
#                iso_map_chared_hadron  = '',                                                                                                                                   
#                iso_map_neutral_hadron = '',                                                                                                                                   
#                iso_map_photon         = '',                                                                                                                                   # 
#                rho = 'fixedGridRhoFastjetAll'):



#  setattr(process,electronCollection+label, cms.EDProducer("patElectronIDIsoSelector",
#                                                           src = cms.InputTag(electronCollection),
#                                                           rho = cms.InputTag(rho),
#                                                           charged_hadron_iso = cms.ItputTag(iso_map_chared_hadron),
#                                                           neutral_hadron_iso = cms.InputTag(iso_map_neutral_hadron),
#                                                           photon_iso         = cms.InputTag(iso_map_photon),
#                                                           id_map             = cms.InputTag(identification_map),
#                                                           relativeIsolationCut = cms.double(0.1),
#                                                           typeID = cms.string(identification_map)
#                                                           )

