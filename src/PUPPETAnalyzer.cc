////////////////////////////////////////////////////////////////////////////////
//
// PUPPETAnalyzer
// ------------------
//
//                        05/2015    Sébastien Brochet <sebastien.brochet@cern.ch>
////////////////////////////////////////////////////////////////////////////////

#include "JMEAnalysis/JMEValidator/interface/PUPPETAnalyzer.h"

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

  if (iConfig.existsAs<edm::InputTag>("srcVertex"))
    srcVertex_ = iConfig.getParameter<edm::InputTag>("srcVertex");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input vertex collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcZboson"))
    srcZboson_ = iConfig.getParameter<edm::InputTag>("srcZboson");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input Zboson collection not given \n";

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

  if (iConfig.existsAs<edm::InputTag>("srcMVAMet"))
    srcMVAMet_ = iConfig.getParameter<edm::InputTag>("srcMVAMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input MVAMet recoil not given \n";

  if (iConfig.existsAs<double>("dRgenMatching"))
    dRgenMatching_ = iConfig.getParameter<double>("dRgenMatching");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input dR for gen jet matching not given \n";

  if(!(srcJet_ == edm::InputTag("")))
    srcJetToken_ = consumes<pat::JetCollection>(srcJet_);

  if(!(srcZboson_ == edm::InputTag("")))
    srcZbosonToken_ = consumes<std::vector<reco::Particle>>(srcZboson_);

  if(!(srcVertex_ == edm::InputTag("")))
    srcVertexToken_ = consumes<reco::VertexCollection>(srcVertex_);

  if(!(srcGenJets_ == edm::InputTag("")))
    srcGenJetsToken_ = consumes<reco::GenJetCollection>(srcGenJets_);

  if(!(srcGenJetsCleaned_ == edm::InputTag("")))
    srcGenJetsCleanedToken_ = consumes<pat::JetCollection>(srcGenJetsCleaned_);

  if(!(srcGenMet_ == edm::InputTag("")))
    srcGenMetToken_ = consumes<pat::METCollection>(srcGenMet_);

  if(!(srcGenParticles_ == edm::InputTag("")))
    srcGenParticlesToken_ = consumes<reco::GenParticleCollection>(srcGenParticles_);

  if(!(srcGenEventInfo_ == edm::InputTag("")))
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

  if(!(srcRecoilPFPuppiMet_ChargedPV_ == edm::InputTag("")))
    srcRecoilPFPuppiMet_ChargedPVToken_ = consumes<pat::METCollection>(srcRecoilPFPuppiMet_ChargedPV_);

  if(!(srcRecoilPFPuppiMet_ChargedPU_ == edm::InputTag("")))
    srcRecoilPFPuppiMet_ChargedPUToken_ = consumes<pat::METCollection>(srcRecoilPFPuppiMet_ChargedPU_);

  if(!(srcRecoilPFPuppiMet_NeutralPV_ == edm::InputTag("")))
    srcRecoilPFPuppiMet_NeutralPVToken_ = consumes<pat::METCollection>(srcRecoilPFPuppiMet_NeutralPV_);

  if(!(srcRecoilPFPuppiMet_NeutralPU_ == edm::InputTag("")))
    srcRecoilPFPuppiMet_NeutralPUToken_ = consumes<pat::METCollection>(srcRecoilPFPuppiMet_NeutralPU_);

  if(!(srcMVAMet_ == edm::InputTag("")))
    srcMVAMetToken_ = consumes<pat::METCollection>(srcMVAMet_);

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

  // take gen level info weight
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

  int ijet = 0;
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

  // store jet info
  edm::Handle<std::vector<pat::Jet>> jetHandle;
  iEvent.getByToken(srcJetToken_, jetHandle);

  NCleanedJets_  = jetHandle->size();

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

  // Gen Recoils
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

  // reco recoils
  edm::Handle<std::vector<pat::MET>> RecoilPFMetHandle;
  iEvent.getByToken(srcRecoilPFMetToken_, RecoilPFMetHandle);

  edm::Handle<std::vector<pat::MET>> PFMetHandle;
  iEvent.getByToken(srcPFMetToken_, PFMetHandle);

  const pat::MET& MetPF = PFMetHandle->at(0);
  PFMet_sumEt_ = MetPF.sumEt();
  PFMet_Pt_    = MetPF.pt();
  PFMet_Phi_   = MetPF.phi();
  PFMet_uncorrected_sumEt_ = MetPF.uncorrectedSumEt();
  PFMet_uncorrected_Pt_    = MetPF.uncorrectedPt();
  PFMet_uncorrected_Phi_   = MetPF.uncorrectedPhi();

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
  recoilPFMet_uncorrected_sumEt_ = recoilMetPF.uncorrectedSumEt();
  recoilPFMet_uncorrected_Pt_    = recoilMetPF.uncorrectedPt();
  recoilPFMet_uncorrected_Phi_   = recoilMetPF.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFMet_uncorrected_Pt_,reco::deltaPhi(recoilPFMet_uncorrected_Phi_,Boson_Phi_));
  recoilPFMet_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFMet_uncorrected_LongZ_ = RecoilVec.Px();
  BosonVec.SetMagPhi(Boson_Pt_, reco::deltaPhi(Boson_Phi_, recoilPFMet_uncorrected_Phi_ + TMath::Pi()));
  recoilPFMet_Boson_PerpU_ = - BosonVec.Py();
  recoilPFMet_Boson_LongU_ = BosonVec.Px() - recoilPFMet_uncorrected_Pt_; 

  // reco recoil CHS met
  edm::Handle<std::vector<pat::MET>> RecoilPFMetCHSHandle;
  iEvent.getByToken(srcRecoilPFCHSMetToken_, RecoilPFMetCHSHandle);

  edm::Handle<std::vector<pat::MET>> PFCHSMetHandle;
  iEvent.getByToken(srcPFCHSMetToken_, PFCHSMetHandle);

  const pat::MET& MetCHSPF = PFCHSMetHandle->at(0);
  PFCHSMet_sumEt_ = MetCHSPF.sumEt();
  PFCHSMet_Pt_    = MetCHSPF.pt();
  PFCHSMet_Phi_   = MetCHSPF.phi();
  PFCHSMet_uncorrected_sumEt_ = MetCHSPF.uncorrectedSumEt();
  PFCHSMet_uncorrected_Pt_    = MetCHSPF.uncorrectedPt();
  PFCHSMet_uncorrected_Phi_   = MetCHSPF.uncorrectedPhi();

  const pat::MET& recoilMetPFCHS = RecoilPFMetCHSHandle->at(0);
  recoilPFCHSMet_sumEt_ = recoilMetPFCHS.sumEt();
  recoilPFCHSMet_Pt_    = recoilMetPFCHS.pt();
  recoilPFCHSMet_Phi_   = recoilMetPFCHS.phi();
  RecoilVec.SetMagPhi(recoilPFCHSMet_Pt_,reco::deltaPhi(recoilPFCHSMet_Phi_,Boson_Phi_));
  recoilPFCHSMet_PerpZ_ = RecoilVec.Py();
  recoilPFCHSMet_LongZ_ = RecoilVec.Px();

  recoilPFCHSMet_uncorrected_sumEt_ = recoilMetPFCHS.uncorrectedSumEt();
  recoilPFCHSMet_uncorrected_Pt_    = recoilMetPFCHS.uncorrectedPt();
  recoilPFCHSMet_uncorrected_Phi_   = recoilMetPFCHS.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFCHSMet_uncorrected_Pt_,reco::deltaPhi(recoilPFCHSMet_uncorrected_Phi_,Boson_Phi_));
  recoilPFCHSMet_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFCHSMet_uncorrected_LongZ_ = RecoilVec.Px();

  // reco recoil puppi met
  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppiHandle;
  iEvent.getByToken(srcRecoilPFPuppiMetToken_, RecoilPFMetPuppiHandle);

  edm::Handle<std::vector<pat::MET>> PFPuppiMetHandle;
  iEvent.getByToken(srcPFPuppiMetToken_, PFPuppiMetHandle);

  const pat::MET& MetPuppiPF = PFPuppiMetHandle->at(0);
  PFPuppiMet_sumEt_ = MetPuppiPF.sumEt();
  PFPuppiMet_Pt_    = MetPuppiPF.pt();
  PFPuppiMet_Phi_   = MetPuppiPF.phi();
  PFPuppiMet_uncorrected_sumEt_ = MetPuppiPF.uncorrectedSumEt();
  PFPuppiMet_uncorrected_Pt_    = MetPuppiPF.uncorrectedPt();
  PFPuppiMet_uncorrected_Phi_   = MetPuppiPF.uncorrectedPhi();


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

  recoilPFPuppiMet_uncorrected_sumEt_ = recoilMetPFPuppi.uncorrectedSumEt();
  recoilPFPuppiMet_uncorrected_Pt_    = recoilMetPFPuppi.uncorrectedPt();
  recoilPFPuppiMet_uncorrected_Phi_   = recoilMetPFPuppi.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_uncorrected_Pt_,reco::deltaPhi(recoilPFPuppiMet_uncorrected_Phi_,Boson_Phi_));
  recoilPFPuppiMet_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_uncorrected_LongZ_ = RecoilVec.Px();
  BosonVec.SetMagPhi(Boson_Pt_, reco::deltaPhi(Boson_Phi_, recoilPFPuppiMet_uncorrected_Phi_ + TMath::Pi()));
  recoilPFPuppiMet_uncorrected_Boson_PerpU_ = - BosonVec.Py();
  recoilPFPuppiMet_uncorrected_Boson_LongU_ = BosonVec.Px() - recoilPFPuppiMet_uncorrected_Pt_; 

  // track met
  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_ChargedPVHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_ChargedPVToken_, RecoilPFMetPuppi_ChargedPVHandle);

  const pat::MET& recoilMetPFPuppi_ChargedPV = RecoilPFMetPuppi_ChargedPVHandle->at(0);
  recoilPFPuppiMet_ChargedPV_sumEt_ = recoilMetPFPuppi_ChargedPV.sumEt();
  recoilPFPuppiMet_ChargedPV_Pt_    = recoilMetPFPuppi_ChargedPV.pt();
  recoilPFPuppiMet_ChargedPV_Phi_   = recoilMetPFPuppi_ChargedPV.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_ChargedPV_Pt_,reco::deltaPhi(recoilPFPuppiMet_ChargedPV_Phi_ ,Boson_Phi_));
  recoilPFPuppiMet_ChargedPV_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_ChargedPV_LongZ_ = RecoilVec.Px();

  // charged from PU
  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_ChargedPUHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_ChargedPUToken_, RecoilPFMetPuppi_ChargedPUHandle);

  const pat::MET& recoilMetPFPuppi_ChargedPU = RecoilPFMetPuppi_ChargedPUHandle->at(0);
  recoilPFPuppiMet_ChargedPU_sumEt_ = recoilMetPFPuppi_ChargedPU.sumEt();
  recoilPFPuppiMet_ChargedPU_Pt_    = recoilMetPFPuppi_ChargedPU.pt();
  recoilPFPuppiMet_ChargedPU_Phi_   = recoilMetPFPuppi_ChargedPU.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_ChargedPU_Pt_,reco::deltaPhi(recoilPFPuppiMet_ChargedPU_Phi_,Boson_Phi_));
  recoilPFPuppiMet_ChargedPU_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_ChargedPU_LongZ_ = RecoilVec.Px();

  // neutral puppi from PV
  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_NeutralPVHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_NeutralPVToken_, RecoilPFMetPuppi_NeutralPVHandle);

  const pat::MET& recoilMetPFPuppi_NeutralPV = RecoilPFMetPuppi_NeutralPVHandle->at(0);
  recoilPFPuppiMet_NeutralPV_sumEt_ = recoilMetPFPuppi_NeutralPV.sumEt();
  recoilPFPuppiMet_NeutralPV_Pt_    = recoilMetPFPuppi_NeutralPV.pt();
  recoilPFPuppiMet_NeutralPV_Phi_   = recoilMetPFPuppi_NeutralPV.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_NeutralPV_Pt_,reco::deltaPhi(recoilPFPuppiMet_NeutralPV_Phi_,Boson_Phi_));
  recoilPFPuppiMet_NeutralPV_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_NeutralPV_LongZ_ = RecoilVec.Px();

  // neutral puppi from PU
  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_NeutralPUHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_NeutralPUToken_, RecoilPFMetPuppi_NeutralPUHandle);

  const pat::MET& recoilMetPFPuppi_NeutralPU = RecoilPFMetPuppi_NeutralPUHandle->at(0);
  recoilPFPuppiMet_NeutralPU_sumEt_ = recoilMetPFPuppi_NeutralPU.sumEt();
  recoilPFPuppiMet_NeutralPU_Pt_    = recoilMetPFPuppi_NeutralPU.pt();
  recoilPFPuppiMet_NeutralPU_Phi_   = recoilMetPFPuppi_NeutralPU.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_NeutralPU_Pt_,reco::deltaPhi(recoilPFPuppiMet_NeutralPU_Phi_,Boson_Phi_));
  recoilPFPuppiMet_NeutralPU_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_NeutralPU_LongZ_ = RecoilVec.Px();

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

  NGenMatchedJets_ = GenMatchedJets_Pt_.size();
  tree.fill();
  
}

////////////////////////////////////////////////////////////////////////////////
// define PUPPETAnalyzer as a plugin
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PUPPETAnalyzer);
