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

  if (iConfig.existsAs<edm::InputTag>("srcGenMet"))
    srcGenMet_ = iConfig.getParameter<edm::InputTag>("srcGenMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input gen MET location not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcGenParticles"))
    srcGenParticles_ = iConfig.getParameter<edm::InputTag>("srcGenParticles");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input gen particle collection not given \n";

  if (iConfig.existsAs<edm::InputTag>("srcRecoilPFMet"))
    srcRecoilPFMet_ = iConfig.getParameter<edm::InputTag>("srcRecoilPFMet");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input PFMet recoil not given \n";

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

  if(!(srcGenMet_ == edm::InputTag("")))
    srcGenMetToken_ = consumes<pat::METCollection>(srcGenMet_);

  if(!(srcGenParticles_ == edm::InputTag("")))
    srcGenParticlesToken_ = consumes<reco::GenParticleCollection>(srcGenParticles_);

  if(!(srcRecoilPFMet_ == edm::InputTag("")))
    srcRecoilPFMetToken_ = consumes<pat::METCollection>(srcRecoilPFMet_);

  if(!(srcRecoilPFCHSMet_ == edm::InputTag("")))
    srcRecoilPFCHSMetToken_ = consumes<pat::METCollection>(srcRecoilPFCHSMet_);

  if(!(srcRecoilPFPuppiMet_ == edm::InputTag("")))
    srcRecoilPFPuppiMetToken_ = consumes<pat::METCollection>(srcRecoilPFPuppiMet_);

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

    edm::Handle<reco::GenParticleCollection> GenParticlesHandle;
  iEvent.getByToken(srcGenParticlesToken_, GenParticlesHandle);

  for(auto GenParticle : *GenParticlesHandle){
    if(GenParticle.pdgId() == 23){
        GenZ_Pt_  = GenParticle.pt();
        GenZ_Eta_ = GenParticle.eta();
        GenZ_Phi_ = GenParticle.phi();
        GenZ_M_   = GenParticle.mass();
    }
  }  

  edm::Handle<reco::GenJetCollection> GenJetsHandle;
  iEvent.getByToken(srcGenJetsToken_, GenJetsHandle);

  NGenJets_ = GenJetsHandle->size();

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
  edm::Handle<std::vector<reco::Particle>> ZbosonHandle;
  iEvent.getByToken(srcZbosonToken_, ZbosonHandle);

  const reco::Particle& Zboson = ZbosonHandle->at(0);
  Zboson_Pt_  = Zboson.pt();
  Zboson_Eta_ = Zboson.eta();
  Zboson_Phi_ = Zboson.phi();
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> p4(Zboson.pt(), Zboson.eta(), Zboson.phi(), Zboson.energy());
  Zboson_M_   = p4.M();
  Zboson_daughter_ = Zboson.pdgId();

  // Recoils

  edm::Handle<std::vector<pat::MET>> GenMetHandle;
  iEvent.getByToken(srcGenMetToken_, GenMetHandle);

  const reco::GenMET* genMET_ = GenMetHandle->at(0).genMET();
  GenRecoil_sumEt_ = genMET_->sumEt() - Zboson.et();
  pat::MET genMETplusZ;
  genMETplusZ.setP4(genMET_->p4() + Zboson.p4());
  GenRecoil_Pt_    = genMETplusZ.pt();
  RecoilVec.SetMagPhi(GenRecoil_Pt_,genMETplusZ.phi() - TMath::Pi());
  GenRecoil_Phi_   = genMETplusZ.phi();
  RecoilVec.SetMagPhi(GenRecoil_Pt_,GenRecoil_Phi_ - Zboson_Phi_);
  GenRecoil_PerpZ_ = RecoilVec.Py();
  GenRecoil_LongZ_ = RecoilVec.Px();


  edm::Handle<std::vector<pat::MET>> RecoilPFMetHandle;
  iEvent.getByToken(srcRecoilPFMetToken_, RecoilPFMetHandle);

  const pat::MET& metPF = RecoilPFMetHandle->at(0);
  recoilPFMet_sumEt_ = metPF.sumEt();
  recoilPFMet_Pt_    = metPF.pt();
  recoilPFMet_Phi_   = metPF.phi();
  RecoilVec.SetMagPhi(recoilPFMet_Pt_,recoilPFMet_Phi_ - Zboson_Phi_);
  recoilPFMet_PerpZ_ = RecoilVec.Py();
  recoilPFMet_LongZ_ = RecoilVec.Px();

  recoilPFMet_uncorrected_sumEt_ = metPF.uncorrectedSumEt();
  recoilPFMet_uncorrected_Pt_    = metPF.uncorrectedPt();
  recoilPFMet_uncorrected_Phi_   = metPF.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFMet_uncorrected_Pt_,recoilPFMet_uncorrected_Phi_ - Zboson_Phi_);
  recoilPFMet_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFMet_uncorrected_LongZ_ = RecoilVec.Px();

  edm::Handle<std::vector<pat::MET>> RecoilPFMetCHSHandle;
  iEvent.getByToken(srcRecoilPFCHSMetToken_, RecoilPFMetCHSHandle);

  const pat::MET& metPFCHS = RecoilPFMetCHSHandle->at(0);
  recoilPFCHSMet_sumEt_ = metPFCHS.sumEt();
  recoilPFCHSMet_Pt_    = metPFCHS.pt();
  recoilPFCHSMet_Phi_   = metPFCHS.phi();
  RecoilVec.SetMagPhi(recoilPFCHSMet_Pt_,recoilPFCHSMet_Phi_ - Zboson_Phi_);
  recoilPFCHSMet_PerpZ_ = RecoilVec.Py();
  recoilPFCHSMet_LongZ_ = RecoilVec.Px();

  recoilPFCHSMet_uncorrected_sumEt_ = metPFCHS.uncorrectedSumEt();
  recoilPFCHSMet_uncorrected_Pt_    = metPFCHS.uncorrectedPt();
  recoilPFCHSMet_uncorrected_Phi_   = metPFCHS.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFCHSMet_uncorrected_Pt_,recoilPFCHSMet_uncorrected_Phi_ - Zboson_Phi_);
  recoilPFCHSMet_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFCHSMet_uncorrected_LongZ_ = RecoilVec.Px();

  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppiHandle;
  iEvent.getByToken(srcRecoilPFPuppiMetToken_, RecoilPFMetPuppiHandle);

  const pat::MET& metPFPuppi = RecoilPFMetPuppiHandle->at(0);
  recoilPFPuppiMet_sumEt_ = metPFPuppi.sumEt();
  recoilPFPuppiMet_Pt_    = metPFPuppi.pt();
  recoilPFPuppiMet_Phi_   = metPFPuppi.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_Pt_,recoilPFPuppiMet_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_LongZ_ = RecoilVec.Px();

  recoilPFPuppiMet_uncorrected_sumEt_ = metPFPuppi.uncorrectedSumEt();
  recoilPFPuppiMet_uncorrected_Pt_    = metPFPuppi.uncorrectedPt();
  recoilPFPuppiMet_uncorrected_Phi_   = metPFPuppi.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_uncorrected_Pt_,recoilPFPuppiMet_uncorrected_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_uncorrected_LongZ_ = RecoilVec.Px();

  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_ChargedPVHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_ChargedPVToken_, RecoilPFMetPuppi_ChargedPVHandle);

  const pat::MET& metPFPuppi_ChargedPV = RecoilPFMetPuppi_ChargedPVHandle->at(0);
  recoilPFPuppiMet_ChargedPV_sumEt_ = metPFPuppi_ChargedPV.sumEt();
  recoilPFPuppiMet_ChargedPV_Pt_    = metPFPuppi_ChargedPV.pt();
  recoilPFPuppiMet_ChargedPV_Phi_   = metPFPuppi_ChargedPV.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_ChargedPV_Pt_,recoilPFPuppiMet_ChargedPV_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_ChargedPV_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_ChargedPV_LongZ_ = RecoilVec.Px();

  recoilPFPuppiMet_ChargedPV_uncorrected_sumEt_ = metPFPuppi_ChargedPV.uncorrectedSumEt();
  recoilPFPuppiMet_ChargedPV_uncorrected_Pt_    = metPFPuppi_ChargedPV.uncorrectedPt();
  recoilPFPuppiMet_ChargedPV_uncorrected_Phi_   = metPFPuppi_ChargedPV.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_ChargedPV_uncorrected_Pt_,recoilPFPuppiMet_ChargedPV_uncorrected_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_ChargedPV_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_ChargedPV_uncorrected_LongZ_ = RecoilVec.Px();

  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_ChargedPUHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_ChargedPUToken_, RecoilPFMetPuppi_ChargedPUHandle);

  const pat::MET& metPFPuppi_ChargedPU = RecoilPFMetPuppi_ChargedPUHandle->at(0);
  recoilPFPuppiMet_ChargedPU_sumEt_ = metPFPuppi_ChargedPU.sumEt();
  recoilPFPuppiMet_ChargedPU_Pt_    = metPFPuppi_ChargedPU.pt();
  recoilPFPuppiMet_ChargedPU_Phi_   = metPFPuppi_ChargedPU.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_ChargedPU_Pt_,recoilPFPuppiMet_ChargedPU_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_ChargedPU_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_ChargedPU_LongZ_ = RecoilVec.Px();

  recoilPFPuppiMet_ChargedPU_uncorrected_sumEt_ = metPFPuppi_ChargedPU.uncorrectedSumEt();
  recoilPFPuppiMet_ChargedPU_uncorrected_Pt_    = metPFPuppi_ChargedPU.uncorrectedPt();
  recoilPFPuppiMet_ChargedPU_uncorrected_Phi_   = metPFPuppi_ChargedPU.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_ChargedPU_uncorrected_Pt_,recoilPFPuppiMet_ChargedPU_uncorrected_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_ChargedPU_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_ChargedPU_uncorrected_LongZ_ = RecoilVec.Px();

  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_NeutralPVHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_NeutralPVToken_, RecoilPFMetPuppi_NeutralPVHandle);

  const pat::MET& metPFPuppi_NeutralPV = RecoilPFMetPuppi_NeutralPVHandle->at(0);
  recoilPFPuppiMet_NeutralPV_sumEt_ = metPFPuppi_NeutralPV.sumEt();
  recoilPFPuppiMet_NeutralPV_Pt_    = metPFPuppi_NeutralPV.pt();
  recoilPFPuppiMet_NeutralPV_Phi_   = metPFPuppi_NeutralPV.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_NeutralPV_Pt_,recoilPFPuppiMet_NeutralPV_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_NeutralPV_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_NeutralPV_LongZ_ = RecoilVec.Px();

  recoilPFPuppiMet_NeutralPV_uncorrected_sumEt_ = metPFPuppi_NeutralPV.uncorrectedSumEt();
  recoilPFPuppiMet_NeutralPV_uncorrected_Pt_    = metPFPuppi_NeutralPV.uncorrectedPt();
  recoilPFPuppiMet_NeutralPV_uncorrected_Phi_   = metPFPuppi_NeutralPV.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_NeutralPV_uncorrected_Pt_,recoilPFPuppiMet_NeutralPV_uncorrected_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_NeutralPV_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_NeutralPV_uncorrected_LongZ_ = RecoilVec.Px();

  edm::Handle<std::vector<pat::MET>> RecoilPFMetPuppi_NeutralPUHandle;
  iEvent.getByToken(srcRecoilPFPuppiMet_NeutralPUToken_, RecoilPFMetPuppi_NeutralPUHandle);

  const pat::MET& metPFPuppi_NeutralPU = RecoilPFMetPuppi_NeutralPUHandle->at(0);
  recoilPFPuppiMet_NeutralPU_sumEt_ = metPFPuppi_NeutralPU.sumEt();
  recoilPFPuppiMet_NeutralPU_Pt_    = metPFPuppi_NeutralPU.pt();
  recoilPFPuppiMet_NeutralPU_Phi_   = metPFPuppi_NeutralPU.phi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_NeutralPU_Pt_,recoilPFPuppiMet_NeutralPU_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_NeutralPU_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_NeutralPU_LongZ_ = RecoilVec.Px();

  recoilPFPuppiMet_NeutralPU_uncorrected_sumEt_ = metPFPuppi_NeutralPU.uncorrectedSumEt();
  recoilPFPuppiMet_NeutralPU_uncorrected_Pt_    = metPFPuppi_NeutralPU.uncorrectedPt();
  recoilPFPuppiMet_NeutralPU_uncorrected_Phi_   = metPFPuppi_NeutralPU.uncorrectedPhi();
  RecoilVec.SetMagPhi(recoilPFPuppiMet_NeutralPU_uncorrected_Pt_,recoilPFPuppiMet_NeutralPU_uncorrected_Phi_ - Zboson_Phi_);
  recoilPFPuppiMet_NeutralPU_uncorrected_PerpZ_ = RecoilVec.Py();
  recoilPFPuppiMet_NeutralPU_uncorrected_LongZ_ = RecoilVec.Px();

  edm::Handle<std::vector<pat::MET>> MVAMetHandle;
  iEvent.getByToken(srcMVAMetToken_, MVAMetHandle);

  const pat::MET& MVAMet = MVAMetHandle->at(0);
  MVAMet_sumEt_ = MVAMet.sumEt();
  MVAMet_Pt_    = MVAMet.pt();
  MVAMet_Phi_   = MVAMet.phi();
  RecoilVec.SetMagPhi(MVAMet_Pt_,MVAMet_Phi_ - Zboson_Phi_);
  MVAMet_PerpZ_ = RecoilVec.Py();
  MVAMet_LongZ_ = RecoilVec.Px();

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
