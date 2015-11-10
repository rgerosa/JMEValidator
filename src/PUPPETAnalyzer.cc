////////////////////////////////////////////////////////////////////////////////
//
// PUPPETAnalyzer
// ------------------
//
//                        05/2015    Sébastien Brochet <sebastien.brochet@cern.ch>
////////////////////////////////////////////////////////////////////////////////

#include "JMEAnalysis/JMEValidator/interface/PUPPETAnalyzer.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include <vector>
#include <iostream>
#include <regex>
#include <string>

////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
PUPPETAnalyzer::PUPPETAnalyzer(const edm::ParameterSet& iConfig):
  JME::Analyzer(iConfig){
  
  if (iConfig.existsAs<edm::InputTag>("srcJet"))
    srcJet_ = iConfig.getParameter<edm::InputTag>("srcJet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input jet collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcJetPF"))
    srcJetPF_ = iConfig.getParameter<edm::InputTag>("srcJetPF");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PF jet collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcVertex"))
    srcVertex_ = iConfig.getParameter<edm::InputTag>("srcVertex");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input vertex collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcZboson"))
    srcZboson_ = iConfig.getParameter<edm::InputTag>("srcZboson");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input Zboson collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcLeptons"))
    srcLeptons_ = iConfig.getParameter<edm::InputTag>("srcLeptons");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input Leptons collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenJets"))
    srcGenJets_ = iConfig.getParameter<edm::InputTag>("srcGenJets");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input gen jet collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenJetsCleaned"))
    srcGenJetsCleaned_ = iConfig.getParameter<edm::InputTag>("srcGenJetsCleaned");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input gen jet cleaned collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenMet"))
    srcGenMet_ = iConfig.getParameter<edm::InputTag>("srcGenMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input gen MET location not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenParticles"))
    srcGenParticles_ = iConfig.getParameter<edm::InputTag>("srcGenParticles");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input gen particle collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenEventInfo"))
    srcGenEventInfo_ = iConfig.getParameter<edm::InputTag>("srcGenEventInfo");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input no gen event info \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFMet"))
    srcRecoilPFMet_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFMet recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcPFMet"))
    srcPFMet_ = iConfig.getParameter<edm::InputTag>("srcPFMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFMet not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcPFCHSMet"))
    srcPFCHSMet_ = iConfig.getParameter<edm::InputTag>("srcPFCHSMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFCHSMet not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcPFPuppiMet"))
    srcPFPuppiMet_ = iConfig.getParameter<edm::InputTag>("srcPFPuppiMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFPuppiMet not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFCHSMet"))
    srcRecoilPFCHSMet_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFCHSMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFCHSMet recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFPuppiMet"))
    srcRecoilPFPuppiMet_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFPuppiMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFPuppiMet recoil not given \n";

/*
  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFPuppiMet_ChargedPV"))
    srcRecoilPFPuppiMet_ChargedPV_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFPuppiMet_ChargedPV");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFPuppiMet_ChargedPV recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFPuppiMet_ChargedPU"))
    srcRecoilPFPuppiMet_ChargedPU_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFPuppiMet_ChargedPU");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFPuppiMet_ChargedPU recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFPuppiMet_NeutralPV"))
    srcRecoilPFPuppiMet_NeutralPV_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFPuppiMet_NeutralPV");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFPuppiMet_NeutralPV recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFPuppiMet_NeutralPU"))
    srcRecoilPFPuppiMet_NeutralPU_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFPuppiMet_NeutralPU");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFPuppiMet_NeutralPU recoil not given \n";
*/
  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFChargedPVNeutralPVPUJetID"))
   srcRecoilPFChargedPVNeutralPVPUJetID_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFChargedPVNeutralPVPUJetID");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFChargedPVNeutralPVPUJetID recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFChargedPUNeutralPUPUJetID"))
   srcRecoilPFChargedPUNeutralPUPUJetID_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFChargedPUNeutralPUPUJetID");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFChargedPUNeutralPUPUJetID recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFChargedPVNeutralPV"))
   srcRecoilPFChargedPVNeutralPV_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFChargedPVNeutralPV");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFChargedPVNeutralPV recoil not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcMVAMet"))
    srcMVAMet_ = iConfig.getParameter<edm::InputTag>("srcMVAMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input MVAMet recoil not given \n";

  if (iConfig.existsAs<double>("dRgenMatching"))
    dRgenMatching_ = iConfig.getParameter<double>("dRgenMatching");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input dR for gen jet matching not given \n";

  if (iConfig.existsAs<bool>("isMC"))
    isMC_ = iConfig.getParameter<bool>("isMC");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input isMC \n";

  if (iConfig.existsAs<edm::InputTag >("srcMetFiltersBits"))
    srcMetFiltersBits_ = iConfig.getParameter<edm::InputTag >("srcMetFiltersBits");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] missing met filter bits  \n";

  if (iConfig.existsAs<edm::InputTag >("srcTriggerBits"))
    srcTriggerBits_ = iConfig.getParameter<edm::InputTag >("srcTriggerBits");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] missing trigger bits  \n";

  if (iConfig.existsAs<edm::InputTag >("srcTriggerPrescales"))
    srcTriggerPrescales_ = iConfig.getParameter<edm::InputTag >("srcTriggerBits");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] missing trigger prescales  \n";
  
  if(!(srcJet_ == edm::InputTag("")))
    srcJetToken_ = consumes<pat::JetCollection>(srcJet_);

  if(!(srcJetPF_ == edm::InputTag("")))
    srcJetPFToken_ = consumes<pat::JetCollection>(srcJetPF_);

  if(!(srcZboson_ == edm::InputTag("")))
    srcZbosonToken_ = consumes<std::vector<reco::Particle>>(srcZboson_);

  if(!(srcLeptons_ == edm::InputTag("")))
    srcLeptonsToken_ = consumes<reco::CandidateView>(srcLeptons_);

  if(!(srcVertex_ == edm::InputTag("")))
    srcVertexToken_ = consumes<reco::VertexCollection>(srcVertex_);

  if(!(srcGenJets_ == edm::InputTag("")) and isMC_)
    srcGenJetsToken_ = consumes<reco::GenJetCollection>(srcGenJets_);

  if(!(srcGenJetsCleaned_ == edm::InputTag("")) and isMC_)
    srcGenJetsCleanedToken_ = consumes<pat::JetCollection>(srcGenJetsCleaned_);

  if(!(srcGenMet_ == edm::InputTag("")) and isMC_)
    srcGenMetToken_ = consumes<pat::METCollection>(srcGenMet_);

  if(!(srcGenParticles_ == edm::InputTag("")) and isMC_)
    srcGenParticlesToken_ = consumes<reco::GenParticleCollection>(srcGenParticles_);

  if(!(srcGenEventInfo_ == edm::InputTag("")) and isMC_)
    srcGenEventInfoToken_ = consumes<GenEventInfoProduct>(srcGenEventInfo_);

  if(!(srcRecoilPFMet_ == edm::InputTag("")))
    srcRecoilPFMetToken_ = consumes<pat::METCollection>(srcRecoilPFMet_);

  if(!(srcPFMet_ == edm::InputTag("")))
    srcPFMetToken_ = consumes<pat::METCollection>(srcPFMet_);

  if(!(srcRecoilPFCHSMet_ == edm::InputTag("")))
    srcRecoilPFCHSMetToken_ = consumes<pat::METCollection>(srcRecoilPFCHSMet_);

  if(!(srcPFCHSMet_ == edm::InputTag("")))
    srcPFCHSMetToken_ = consumes<pat::METCollection>(srcPFCHSMet_);

  if(!(srcRecoilPFPuppiMet_ == edm::InputTag("")))
    srcRecoilPFPuppiMetToken_ = consumes<pat::METCollection>(srcRecoilPFPuppiMet_);

  if(!(srcPFPuppiMet_ == edm::InputTag("")))
    srcPFPuppiMetToken_ = consumes<pat::METCollection>(srcPFPuppiMet_);

  if(!(srcRecoilPFChargedPVNeutralPVPUJetID_ == edm::InputTag("")))
    srcRecoilPFChargedPVNeutralPVPUJetIDToken_ = consumes<pat::METCollection>(srcRecoilPFChargedPVNeutralPVPUJetID_);

  if(!(srcRecoilPFChargedPUNeutralPUPUJetID_ == edm::InputTag("")))
    srcRecoilPFChargedPUNeutralPUPUJetIDToken_ = consumes<pat::METCollection>(srcRecoilPFChargedPUNeutralPUPUJetID_);

  if(!(srcRecoilPFChargedPVNeutralPV_ == edm::InputTag("")))
    srcRecoilPFChargedPVNeutralPVToken_ = consumes<pat::METCollection>(srcRecoilPFChargedPVNeutralPV_);

  if(!(srcMVAMet_ == edm::InputTag("")))
    srcMVAMetToken_ = consumes<pat::METCollection>(srcMVAMet_);

  if(!(srcMetFiltersBits_ == edm::InputTag("")))
      srcMetFiltersBitsToken_ = consumes<edm::TriggerResults>(srcMetFiltersBits_);

  if(!(srcTriggerBits_ == edm::InputTag("")))
      srcTriggerBitsToken_ = consumes<edm::TriggerResults>(srcTriggerBits_);

  if(!(srcTriggerPrescales_ == edm::InputTag("")))
    srcTriggerPrescalesToken_ = consumes<pat::PackedTriggerPrescales>(srcTriggerPrescales_);

}


//______________________________________________________________________________
PUPPETAnalyzer::~PUPPETAnalyzer(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________

//______________________________________________________________________________
void PUPPETAnalyzer::analyze(const edm::Event& iEvent,
			     const edm::EventSetup& iSetup){


  // met filters
  edm::Handle<edm::TriggerResults> triggerBit;
  iEvent.getByToken(srcMetFiltersBitsToken_,triggerBit);
  const edm::TriggerNames &metNames = iEvent.triggerNames(*triggerBit);
  for(size_t ibit = 0; ibit < triggerBit->size(); ibit++){
    if( metNames.triggerName(ibit) == "Flag_trackingFailureFilter")
      flag_trackingFailureFilter_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) ==  "Flag_goodVertices")
      flag_goodVertices_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_CSCTightHaloFilter")
      flag_CSCTightHaloFilter_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_trkPOGFilters") 
      flag_trkPOGFilters_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_trkPOG_logErrorTooManyClusters") 
      flag_trkPOG_logErrorTooManyClusters_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_EcalDeadCellTriggerPrimitiveFilter") 
      flag_EcalDeadCellTriggerPrimitiveFilter_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_ecalLaserCorrFilter") 
      flag_ecalLaserCorrFilter_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_trkPOG_manystripclus53X") 
      flag_trkPOG_manystripclus53X_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_eeBadScFilter") 
      flag_eeBadScFilter_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_METFilters") 
      flag_METFilters_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_trkPOG_toomanystripclus53X") 
      flag_trkPOG_toomanystripclus53X_ = triggerBit->accept(ibit);
    else if(metNames.triggerName(ibit) == "Flag_hcalLaserEventFilter") 
      flag_hcalLaserEventFilter_ = triggerBit->accept(ibit);
    else if( metNames.triggerName(ibit) == "Flag_HBHENoiseFilter") 
      flag_HBHENoiseFilter_ = triggerBit->accept(ibit);
  }

  // trigger bits 
  iEvent.getByToken(srcTriggerBitsToken_,triggerBit);
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerBit);
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(srcTriggerPrescalesToken_, triggerPrescales);

  for (size_t iBit = 0 ; iBit < triggerBit->size(); iBit++) {
    if (triggerPrescales.isValid()){
      if(triggerPrescales->getPrescaleForIndex(iBit) !=1) 
	continue;
    }

    if((TString(triggerNames.triggerName(iBit)).Contains("HLT_DoubleMu") or TString(triggerNames.triggerName(iBit)).Contains("HLT_Mu")) and triggerBit->accept(iBit)){
      DoubleMuPaths_.push_back(triggerNames.triggerName(iBit));
    }
    else if((TString(triggerNames.triggerName(iBit)).Contains("HLT_DoubleEle") or TString(triggerNames.triggerName(iBit)).Contains("HLT_Ele")) and triggerBit->accept(iBit)){
      DoubleElePaths_.push_back(triggerNames.triggerName(iBit));
    }
    else if((TString(triggerNames.triggerName(iBit)).Contains("Double") and TString(triggerNames.triggerName(iBit)).Contains("Tau")) and triggerBit->accept(iBit)){
      DoubleTauPaths_.push_back(triggerNames.triggerName(iBit));
    }
  }


  // take gen level info weight  
  int ijet = 0;

  if(isMC_){
    edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
    iEvent.getByToken(srcGenEventInfoToken_,GenEventInfoHandle);
    eventMCWeight = GenEventInfoHandle->weight()/fabs(GenEventInfoHandle->weight());

    // find generator level Z kinematics
    edm::Handle<reco::GenParticleCollection> GenParticlesHandle;
    iEvent.getByToken(srcGenParticlesToken_, GenParticlesHandle);
    reco::GenParticle GenBoson;

    for(auto aGenParticle : *GenParticlesHandle){
      if(aGenParticle.pdgId() == 23 or abs(aGenParticle.pdgId()) == 24 or aGenParticle.pdgId() == 22){
	if(aGenParticle.numberOfDaughters() !=2) continue;
	bool goodBoson = false;
	int numberOfLeptonDaughter   = 0;
	int leptonDaughterPdgId      = 0;
	int numberOfNeutrinoDaughter = 0;
	for(unsigned int i0 = 0; i0 < aGenParticle.numberOfDaughters(); i0++) {
	  const reco::GenParticle *daughter = aGenParticle.daughterRef(i0).get();
	  if(abs(daughter->pdgId()) == 11 or abs(daughter->pdgId()) == 13 or abs(daughter->pdgId()) == 15){
	    numberOfLeptonDaughter ++;
	    leptonDaughterPdgId = abs(daughter->pdgId());
	}
	  else if(abs(daughter->pdgId()) == 12 or abs(daughter->pdgId()) == 14 or abs(daughter->pdgId()) == 16)
	    numberOfNeutrinoDaughter++;
	}
	
	if(numberOfLeptonDaughter == 1 and numberOfNeutrinoDaughter == 1)
	  goodBoson = true;

	if(numberOfLeptonDaughter == 2 and numberOfNeutrinoDaughter == 0)
	  goodBoson = true;
	
	if(goodBoson == false) continue;

	GenBoson_Pt_  = aGenParticle.pt();
	GenBoson_Eta_ = aGenParticle.eta();
	GenBoson_Phi_ = aGenParticle.phi();
	GenBoson_M_   = aGenParticle.mass();
	GenBoson.setP4(aGenParticle.p4());
	GenBoson_daughter_ = leptonDaughterPdgId;
      }
    }
  

    // store gen jets and gen jets cleaned info
    edm::Handle<reco::GenJetCollection> GenJetsHandle;
    iEvent.getByToken(srcGenJetsToken_, GenJetsHandle);
    
    edm::Handle<pat::JetCollection> GenJetsCleanedHandle;
    iEvent.getByToken(srcGenJetsCleanedToken_, GenJetsCleanedHandle);
    
    NGenJets_        = GenJetsHandle->size();
    NGenJetsCleaned_ = GenJetsHandle->size();

    for( auto GenJet : *GenJetsHandle){     
      if(ijet == 0){
	GenLeadingJet_Pt_  = GenJet.pt();
	GenLeadingJet_Eta_ = GenJet.eta();
	GenLeadingJet_Phi_ = GenJet.phi();
	GenLeadingJet_M_   = GenJet.mass();
	ijet++;
	continue;
      }
      else if(ijet == 1){
	GenTrailingJet_Pt_  = GenJet.pt();
	GenTrailingJet_Eta_ = GenJet.eta();
	GenTrailingJet_Phi_ = GenJet.phi();
	GenTrailingJet_M_   = GenJet.mass();
	break;
      }
    }
    
    ijet = 0;
    
    for( auto GenJet : *GenJetsCleanedHandle){     
      if(ijet == 0){
	GenLeadingJetCleaned_Pt_  = GenJet.pt();
	GenLeadingJetCleaned_Eta_ = GenJet.eta();
	GenLeadingJetCleaned_Phi_ = GenJet.phi();
	GenLeadingJetCleaned_M_   = GenJet.mass();
	ijet++;
	continue;
      }
      else if(ijet == 1){
	GenTrailingJetCleaned_Pt_  = GenJet.pt();
	GenTrailingJetCleaned_Eta_ = GenJet.eta();
	GenTrailingJetCleaned_Phi_ = GenJet.phi();
	GenTrailingJetCleaned_M_   = GenJet.mass();
	break;
      }
    }


    edm::Handle<std::vector<pat::MET>> GenMetHandle;
    iEvent.getByToken(srcGenMetToken_, GenMetHandle);
    
    const reco::GenMET* genMET_ = GenMetHandle->at(0).genMET();
    GenRecoil_sumEt_ = genMET_->sumEt();
    
    pat::MET genRecoil;
    genRecoil.setP4(-GenBoson.p4()-genMET_->p4());
    GenRecoil_Pt_    = genRecoil.pt();
    GenRecoil_Phi_   = genRecoil.phi();
    RecoilVec.SetMagPhi(GenRecoil_Pt_,reco::deltaPhi(GenRecoil_Phi_,GenBoson.phi()));
    GenRecoil_PerpZ_ = RecoilVec.Py();
    GenRecoil_LongZ_ = RecoilVec.Px();
    
    BosonVec.SetMagPhi(GenBoson.pt(),reco::deltaPhi(GenBoson.phi(),GenRecoil_Phi_));
    GenBoson_PerpU_ = BosonVec.Py();
    GenBoson_LongU_ = BosonVec.Px() - GenBoson.pt();
  }


  // store jet info
  edm::Handle<std::vector<pat::Jet>> jetHandle;
  iEvent.getByToken(srcJetToken_, jetHandle);

  // store jet PF info
  edm::Handle<std::vector<pat::Jet>> jetPFHandle;
  iEvent.getByToken(srcJetPFToken_, jetPFHandle);

  NCleanedJets_  = jetHandle->size();
  NCleanedJetsPF_  = jetPFHandle->size();

  ijet = 0;
  for( auto jet : *jetHandle){     
    if(ijet == 0){
      LeadingJet_Pt_  = jet.pt();
      LeadingJet_Eta_ = jet.eta();
      LeadingJet_Phi_ = jet.phi();
      LeadingJet_M_   = jet.mass();
      ijet++;
      continue;
    }
    else if(ijet == 1){
      TrailingJet_Pt_  = jet.pt();
      TrailingJet_Eta_ = jet.eta();
      TrailingJet_Phi_ = jet.phi();
      TrailingJet_M_   = jet.mass();
      break;
    }
  }


  ijet = 0;
  for( auto jet : *jetPFHandle){     
    if(ijet == 0){
      LeadingJetPF_Pt_  = jet.pt();
      LeadingJetPF_Eta_ = jet.eta();
      LeadingJetPF_Phi_ = jet.phi();
      LeadingJetPF_M_   = jet.mass();
      ijet++;
      continue;
    }
    else if(ijet == 1){
      TrailingJetPF_Pt_  = jet.pt();
      TrailingJetPF_Eta_ = jet.eta();
      TrailingJetPF_Phi_ = jet.phi();
      TrailingJetPF_M_   = jet.mass();
      break;
    }
  }

  // vertexes
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(srcVertexToken_, vertexHandle);
  NVertex_  = vertexHandle->size();

  // Zboson
  edm::Handle<std::vector<reco::Particle>> BosonHandle;
  iEvent.getByToken(srcZbosonToken_, BosonHandle);

  const reco::Particle& Boson = BosonHandle->at(0);
  Boson_Pt_  = Boson.pt();
  Boson_Eta_ = Boson.eta();
  Boson_Phi_ = Boson.phi();
  Boson_M_   = Boson.p4().M();
  Boson_daughter_ = Boson.pdgId();

  // Leptons
  edm::Handle<reco::CandidateView> LeptonHandle;
  iEvent.getByToken(srcLeptonsToken_, LeptonHandle);

  int iLep = 0;
  for(reco::CandidateView::const_iterator lepton = LeptonHandle->begin(); lepton != LeptonHandle->end(); ++lepton){

    if(iLep == 0){
      LeadingLepton_Pt_  = (*lepton).pt();
      LeadingLepton_Eta_ = (*lepton).eta();
      LeadingLepton_Phi_ = (*lepton).phi();
      LeadingLepton_M_   = (*lepton).p4().M();
    }
    else if(iLep == 1){

      if((*lepton).pt() > LeadingLepton_Pt_){

	TrailingLepton_Pt_  = LeadingLepton_Pt_;
	TrailingLepton_Eta_ = LeadingLepton_Eta_;
	TrailingLepton_Phi_ = LeadingLepton_Phi_;
	TrailingLepton_M_   = LeadingLepton_M_;
	
	LeadingLepton_Pt_  = (*lepton).pt();
	LeadingLepton_Eta_ = (*lepton).eta();
	LeadingLepton_Phi_ = (*lepton).phi();
	LeadingLepton_M_   = (*lepton).p4().M();
      }
      else {
	TrailingLepton_Pt_  = (*lepton).pt();
	TrailingLepton_Eta_ = (*lepton).eta();
	TrailingLepton_Phi_ = (*lepton).phi();
	TrailingLepton_M_   = (*lepton).p4().M();
      }
    }
    iLep++;
  }
    
  // reco recoils
  edm::Handle<std::vector<pat::MET>> RecoilPFMetHandle;
  iEvent.getByToken(srcRecoilPFMetToken_, RecoilPFMetHandle);

  edm::Handle<std::vector<pat::MET>> PFMetHandle;
  iEvent.getByToken(srcPFMetToken_, PFMetHandle);

  const pat::MET& MetPF = PFMetHandle->at(0);
  PFMet_sumEt_ = MetPF.sumEt();
  PFMet_Pt_    = MetPF.pt();
  PFMet_Phi_   = MetPF.phi();
//  PFMet_uncorrected_sumEt_ = MetPF.uncorrectedSumEt();
//  PFMet_uncorrected_Pt_    = MetPF.uncorrectedPt();
//  PFMet_uncorrected_Phi_   = MetPF.uncorrectedPhi();

  const pat::MET& recoilMetPF = RecoilPFMetHandle->at(0);
  recoilPFMet_sumEt_ = recoilMetPF.sumEt();
  recoilPFMet_Pt_    = recoilMetPF.pt();
  recoilPFMet_Phi_   = recoilMetPF.phi();

  RecoilVec.SetMagPhi(recoilPFMet_Pt_,reco::deltaPhi(recoilPFMet_Phi_ ,Boson_Phi_));
  recoilPFMet_PerpZ_ = RecoilVec.Py();
  recoilPFMet_LongZ_ = RecoilVec.Px();
  BosonVec.SetMagPhi(Boson_Pt_, reco::deltaPhi(Boson_Phi_, recoilPFMet_Phi_ + TMath::Pi()));
  recoilPFMet_Boson_PerpU_ = - BosonVec.Py();
  recoilPFMet_Boson_LongU_ = BosonVec.Px() - recoilPFMet_Pt_; 
//  pat::MET recoilMetPF_uncorrected = getUncorrectedRecoil(recoilMetPF);
//  recoilPFMet_uncorrected_sumEt_ = recoilMetPF_uncorrected.sumEt();
//  recoilPFMet_uncorrected_Pt_    = recoilMetPF_uncorrected.pt();
//  recoilPFMet_uncorrected_Phi_   = recoilMetPF_uncorrected.phi();
//  RecoilVec.SetMagPhi(recoilPFMet_uncorrected_Pt_,reco::deltaPhi(recoilPFMet_uncorrected_Phi_,Boson_Phi_));
//  recoilPFMet_uncorrected_PerpZ_ = RecoilVec.Py();
//  recoilPFMet_uncorrected_LongZ_ = RecoilVec.Px();
  // fix this
 // BosonVec.SetMagPhi(Boson_Pt_, reco::deltaPhi(Boson_Phi_, recoilPFMet_uncorrected_Phi_ + TMath::Pi()));
//  recoilPFMet_uncorrected_Boson_PerpU_ = - BosonVec.Py();
//  recoilPFMet_uncorrected_Boson_LongU_ = BosonVec.Px() - recoilPFMet_uncorrected_Pt_; 

  // reco recoil CHS met
  edm::Handle<std::vector<pat::MET>> RecoilPFMetCHSHandle;
  iEvent.getByToken(srcRecoilPFCHSMetToken_, RecoilPFMetCHSHandle);

  edm::Handle<std::vector<pat::MET>> PFCHSMetHandle;
  iEvent.getByToken(srcPFCHSMetToken_, PFCHSMetHandle);

  const pat::MET& MetCHSPF = PFCHSMetHandle->at(0);
  PFCHSMet_sumEt_ = MetCHSPF.sumEt();
  PFCHSMet_Pt_    = MetCHSPF.pt();
  PFCHSMet_Phi_   = MetCHSPF.phi();
/*
  PFCHSMet_uncorrected_sumEt_ = MetCHSPF.uncorrectedSumEt();
  PFCHSMet_uncorrected_Pt_    = MetCHSPF.uncorrectedPt();
  PFCHSMet_uncorrected_Phi_   = MetCHSPF.uncorrectedPhi();
*/
  const pat::MET& recoilMetPFCHS = RecoilPFMetCHSHandle->at(0);
  recoilPFCHSMet_sumEt_ = recoilMetPFCHS.sumEt();
  recoilPFCHSMet_Pt_    = recoilMetPFCHS.pt();
  recoilPFCHSMet_Phi_   = recoilMetPFCHS.phi();
  RecoilVec.SetMagPhi(recoilPFCHSMet_Pt_,reco::deltaPhi(recoilPFCHSMet_Phi_,Boson_Phi_));
  recoilPFCHSMet_PerpZ_ = RecoilVec.Py();
  recoilPFCHSMet_LongZ_ = RecoilVec.Px();

//  pat::MET recoilMetPFCHS_uncorrected = getUncorrectedRecoil( recoilMetPFCHS);
  //recoilPFCHSMet_uncorrected_sumEt_ = recoilMetPFCHS_uncorrected.sumEt();
  //recoilPFCHSMet_uncorrected_Pt_    = recoilMetPFCHS_uncorrected.pt();
  //recoilPFCHSMet_uncorrected_Phi_   = recoilMetPFCHS_uncorrected.phi();
 // RecoilVec.SetMagPhi(recoilPFCHSMet_uncorrected_Pt_,reco::deltaPhi(recoilPFCHSMet_uncorrected_Phi_,Boson_Phi_));
  //recoilPFCHSMet_uncorrected_PerpZ_ = RecoilVec.Py();
  //recoilPFCHSMet_uncorrected_LongZ_ = RecoilVec.Px();

  // reco recoil puppi met
  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppiHandle;
  iEvent.getByToken(srcRecoilPFPuppiMetToken_, RecoilPFMetPuppiHandle);

  edm::Handle<std::vector<pat::MET>> PFPuppiMetHandle;
  iEvent.getByToken(srcPFPuppiMetToken_, PFPuppiMetHandle);

  const pat::MET& MetPuppiPF = PFPuppiMetHandle->at(0);
  PFPuppiMet_sumEt_ = MetPuppiPF.sumEt();
  PFPuppiMet_Pt_    = MetPuppiPF.pt();
  PFPuppiMet_Phi_   = MetPuppiPF.phi();
/*
  PFPuppiMet_uncorrected_sumEt_ = MetPuppiPF.uncorrectedSumEt();
  PFPuppiMet_uncorrected_Pt_    = MetPuppiPF.uncorrectedPt();
  PFPuppiMet_uncorrected_Phi_   = MetPuppiPF.uncorrectedPhi();
*/
  const pat::MET& recoilMetPFPuppi = RecoilPFMetPuppiHandle->at(0);
  recoilPFPuppiMet_sumEt_ = recoilMetPFPuppi.sumEt();
  recoilPFPuppiMet_Pt_    = recoilMetPFPuppi.pt();
  recoilPFPuppiMet_Phi_   = recoilMetPFPuppi.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_Pt_,reco::deltaPhi(recoilPFPuppiMet_Phi_,Boson_Phi_));
  recoilPFPuppiMet_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_LongZ_ = RecoilVec.Px();
  BosonVec.SetMagPhi(Boson_Pt_, reco::deltaPhi(Boson_Phi_, recoilPFPuppiMet_Phi_ + TMath::Pi()));
  recoilPFPuppiMet_Boson_PerpU_ = - BosonVec.Py();
  recoilPFPuppiMet_Boson_LongU_ = BosonVec.Px() - recoilPFPuppiMet_Pt_; 

  // charged PV + neutrals in jet passing PUJet ID
  edm::Handle<std::vector<pat::MET>> RecoilPFMet_ChargedPVNeutralPVPUJetIDHandle;
  iEvent.getByToken(srcRecoilPFChargedPVNeutralPVPUJetIDToken_, RecoilPFMet_ChargedPVNeutralPVPUJetIDHandle);

  const pat::MET& recoilMetPF_ChargedPVNeutralPVPUJetID = RecoilPFMet_ChargedPVNeutralPVPUJetIDHandle->at(0);
  recoilPFMet_ChargedPVNeutralPVPUJetID_sumEt_ = recoilMetPF_ChargedPVNeutralPVPUJetID.sumEt();
  recoilPFMet_ChargedPVNeutralPVPUJetID_Pt_    = recoilMetPF_ChargedPVNeutralPVPUJetID.pt();
  recoilPFMet_ChargedPVNeutralPVPUJetID_Phi_   = recoilMetPF_ChargedPVNeutralPVPUJetID.phi();
  RecoilVec.SetMagPhi(recoilPFMet_ChargedPVNeutralPVPUJetID_Pt_,reco::deltaPhi(recoilPFMet_ChargedPVNeutralPVPUJetID_Phi_,Boson_Phi_));
  recoilPFMet_ChargedPVNeutralPVPUJetID_PerpZ_ = RecoilVec.Py();
  recoilPFMet_ChargedPVNeutralPVPUJetID_LongZ_ = RecoilVec.Px();

  // charged PU + neutrals in jet failing PUJet ID
  edm::Handle<std::vector<pat::MET>> RecoilPFMet_ChargedPUNeutralPUPUJetIDHandle;
  iEvent.getByToken(srcRecoilPFChargedPUNeutralPUPUJetIDToken_, RecoilPFMet_ChargedPUNeutralPUPUJetIDHandle);

  const pat::MET& recoilMetPF_ChargedPUNeutralPUPUJetID = RecoilPFMet_ChargedPUNeutralPUPUJetIDHandle->at(0);
  recoilPFMet_ChargedPUNeutralPUPUJetID_sumEt_ = recoilMetPF_ChargedPUNeutralPUPUJetID.sumEt();
  recoilPFMet_ChargedPUNeutralPUPUJetID_Pt_    = recoilMetPF_ChargedPUNeutralPUPUJetID.pt();
  recoilPFMet_ChargedPUNeutralPUPUJetID_Phi_   = recoilMetPF_ChargedPUNeutralPUPUJetID.phi();
  RecoilVec.SetMagPhi(recoilPFMet_ChargedPUNeutralPUPUJetID_Pt_,reco::deltaPhi(recoilPFMet_ChargedPUNeutralPUPUJetID_Phi_,Boson_Phi_));
  recoilPFMet_ChargedPUNeutralPUPUJetID_PerpZ_ = RecoilVec.Py();
  recoilPFMet_ChargedPUNeutralPUPUJetID_LongZ_ = RecoilVec.Px();

  // charged PV + neutrals PV as all neutrals -  neutrals in jet failing PVJet ID
  edm::Handle<std::vector<pat::MET>> RecoilPFMet_ChargedPVNeutralPVHandle;
  iEvent.getByToken(srcRecoilPFChargedPVNeutralPVToken_, RecoilPFMet_ChargedPVNeutralPVHandle);

  const pat::MET& recoilMetPF_ChargedPVNeutralPV = RecoilPFMet_ChargedPVNeutralPVHandle->at(0);
  recoilPFMet_ChargedPVNeutralPV_sumEt_ = recoilMetPF_ChargedPVNeutralPV.sumEt();
  recoilPFMet_ChargedPVNeutralPV_Pt_    = recoilMetPF_ChargedPVNeutralPV.pt();
  recoilPFMet_ChargedPVNeutralPV_Phi_   = recoilMetPF_ChargedPVNeutralPV.phi();
  RecoilVec.SetMagPhi(recoilPFMet_ChargedPVNeutralPV_Pt_,reco::deltaPhi(recoilPFMet_ChargedPVNeutralPV_Phi_,Boson_Phi_));
  recoilPFMet_ChargedPVNeutralPV_PerpZ_ = RecoilVec.Py();
  recoilPFMet_ChargedPVNeutralPV_LongZ_ = RecoilVec.Px();

  // MVA met
  edm::Handle<std::vector<pat::MET>> MVAMetHandle;
  iEvent.getByToken(srcMVAMetToken_, MVAMetHandle);

  const pat::MET& MVAMet = MVAMetHandle->at(0);
  MVAMet_sumEt_ = MVAMet.sumEt();
  MVAMet_Pt_    = MVAMet.pt();
  MVAMet_Phi_   = MVAMet.phi();
  RecoilVec.SetMagPhi(MVAMet_Pt_,reco::deltaPhi(MVAMet_Phi_,Boson_Phi_));
  MVAMet_PerpZ_ = RecoilVec.Py();
  MVAMet_LongZ_ = RecoilVec.Px();


  // dump all jet info
  AllJets_Pt_.clear();
  AllJets_Eta_.clear();
  AllJets_Phi_.clear();
  AllJets_M_.clear();
  GenMatchedJets_Pt_.clear();
  GenMatchedJets_Eta_.clear();
  GenMatchedJets_Phi_.clear();
  GenMatchedJets_M_.clear();

  for(auto jet : *jetHandle){
    AllJets_Pt_.push_back(jet.pt());
    AllJets_Eta_.push_back(jet.eta());
    AllJets_Phi_.push_back(jet.phi());
    AllJets_M_.push_back(jet.mass());
    if(isMC_){
      edm::Handle<reco::GenJetCollection> GenJetsHandle;
      iEvent.getByToken(srcGenJetsToken_, GenJetsHandle);
      for(auto GenJet : *GenJetsHandle){
	if(reco::deltaR(jet.eta(),jet.phi(),GenJet.eta(),GenJet.phi()) < dRgenMatching_){
	  GenMatchedJets_Pt_.push_back(GenJet.pt());
	  GenMatchedJets_Eta_.push_back(GenJet.eta());
	  GenMatchedJets_Phi_.push_back(GenJet.phi());
	  GenMatchedJets_M_.push_back(GenJet.mass());
	  break;
	}
      }
    }
  }

  NGenMatchedJets_ = GenMatchedJets_Pt_.size();

  // dump all PF jet info
  AllJetsPF_Pt_.clear();
  AllJetsPF_Eta_.clear();
  AllJetsPF_Phi_.clear();
  AllJetsPF_M_.clear();
  GenMatchedJetsPF_Pt_.clear();
  GenMatchedJetsPF_Eta_.clear();
  GenMatchedJetsPF_Phi_.clear();
  GenMatchedJetsPF_M_.clear();

  for(auto jet : *jetPFHandle){
    AllJetsPF_Pt_.push_back(jet.pt());
    AllJetsPF_Eta_.push_back(jet.eta());
    AllJetsPF_Phi_.push_back(jet.phi());
    AllJetsPF_M_.push_back(jet.mass());
    if(isMC_){
      edm::Handle<reco::GenJetCollection> GenJetsHandle;
      iEvent.getByToken(srcGenJetsToken_, GenJetsHandle);
      for(auto GenJet : *GenJetsHandle){
	if(reco::deltaR(jet.eta(),jet.phi(),GenJet.eta(),GenJet.phi()) < dRgenMatching_){
	  GenMatchedJetsPF_Pt_.push_back(GenJet.pt());
	  GenMatchedJetsPF_Eta_.push_back(GenJet.eta());
	  GenMatchedJetsPF_Phi_.push_back(GenJet.phi());
	  GenMatchedJetsPF_M_.push_back(GenJet.mass());
	  break;
	}
      }
    }
  }

  NGenMatchedJetsPF_ = GenMatchedJetsPF_Pt_.size();

  tree.fill();
  
}

pat::MET PUPPETAnalyzer::getUncorrectedRecoil(const pat::MET& input)
{
  // turn recoil in direction of MET
  TVector2 helperVector2(input.px(), input.py());
  helperVector2 = helperVector2.Rotate(M_PI);
  reco::Candidate::LorentzVector helperVector4( helperVector2.Px(), helperVector2.Py(), 0, input.sumEt());
  pat::MET negRecoil(input);
  negRecoil.setP4(helperVector4);

  //retrieve uncorrected values and rotate them back by PI
  //helperVector2.SetMagPhi( negRecoil.uncorrectedPt(), negRecoil.uncorrectedPhi());
  helperVector2 = helperVector2.Rotate(- M_PI);
  helperVector4.SetXYZT(helperVector2.Px(), helperVector2.Py(), 0, negRecoil.sumEt());
  pat::MET uncorrectedRecoil(input);
  uncorrectedRecoil.setP4(helperVector4);
  //uncorrectedRecoil.setSumEt(negRecoil.uncorrectedSumEt());
  return uncorrectedRecoil;
}
////////////////////////////////////////////////////////////////////////////////
// define PUPPETAnalyzer as a plugin
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PUPPETAnalyzer);
