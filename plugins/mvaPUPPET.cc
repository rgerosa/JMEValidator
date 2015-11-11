#include <iostream>

#include "JMEAnalysis/JMEValidator/plugins/mvaPUPPET.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/Particle.h"

typedef std::vector<reco::Particle> ParticleCollection;

mvaPUPPET::mvaPUPPET(const edm::ParameterSet& cfg){

  // get MET that the mva is applied on
  referenceMET_      = consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("referenceMET"));
  referenceMET_name_ = cfg.getParameter<edm::InputTag>("referenceMET").label();

  if (cfg.existsAs<bool>("produceRecoils"))
    produceRecoils_ = cfg.getParameter<bool>("produceRecoils");
  else
    produceRecoils_ = false;


  // get tokens for input METs and prepare for saving the corresponding recoils to the event
  srcMETTags_   = cfg.getParameter<vInputTag>("srcMETs");
  for(vInputTag::const_iterator it=srcMETTags_.begin();it!=srcMETTags_.end();it++) {
    srcMETs_.push_back( consumes<pat::METCollection >( *it ) );
    if(produceRecoils_)
      produces<pat::METCollection>( "recoil"+(*it).label() );
  }

  // take flags for the met
  if (cfg.existsAs<std::vector<int> >("inputMETFlags"))
    srcMETFlags_ = cfg.getParameter<std::vector<int>>("inputMETFlags");
  
  if(srcMETFlags_.size() != srcMETTags_.size()+1)
    throw cms::Exception("mvaPUPPET::mvaPUPPET") << " Failed to load MET flags   !!\n";

  //get leptons to calculate Z vector and save it as a reco::candidate back to the event
  vInputTag srcLeptonsTags = cfg.getParameter<vInputTag>("srcLeptons");
  for(vInputTag::const_iterator it=srcLeptonsTags.begin();it!=srcLeptonsTags.end();it++) {
    srcLeptons_.push_back( consumes<reco::CandidateView >( *it ) );
  }
    
  // get debug flags
  if (cfg.existsAs<bool>("debug"))
    debug_ = cfg.getParameter<bool>("debug");
  else
    debug_ = false;

  // get vertices
  srcVertices_  = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("srcVertices"));
  // get jets 
  srcJets_      = consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("srcJets"));
  // get taus
  srcTaus_      = consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("srcTaus"));
  // get muons
  srcMuons_     = consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("srcMuons"));

  // load weight files
  edm::ParameterSet cfgInputFileNames = cfg.getParameter<edm::ParameterSet>("inputFileNames");
  
  if(cfgInputFileNames.existsAs<edm::FileInPath>("PhiCorrectionWeightFile")){
    inputFileNamePhiCorrection_ = cfgInputFileNames.getParameter<edm::FileInPath>("PhiCorrectionWeightFile");
    mvaReaderPhiCorrection_     = loadMVAfromFile(inputFileNamePhiCorrection_, variablesForPhiTraining_, "PhiCorrectedRecoil");
  }
  
  if(cfgInputFileNames.existsAs<edm::FileInPath>("RecoilCorrectionWeightFile")){
    inputFileNameRecoilCorrection_ = cfgInputFileNames.getParameter<edm::FileInPath>("RecoilCorrectionWeightFile");
    mvaReaderRecoilCorrection_  = loadMVAfromFile(inputFileNameRecoilCorrection_, variablesForRecoilTraining_, "LongZCorrectedRecoil");
  }
  
  // prepare for saving the final mvaMET to the event
  if(cfg.existsAs<std::string>("mvaMETLabel"))
    mvaMETLabel_ = cfg.getParameter<std::string>("mvaMETLabel");
  else
    mvaMETLabel_ = "mvaMET";

  produces<pat::METCollection>(mvaMETLabel_);

  // produce Zboson tag
  if(cfg.existsAs<std::string>("ZbosonLabel"))
    ZbosonLabel_ = cfg.getParameter<std::string>("ZbosonLabel");
  else
    ZbosonLabel_ = "ZtagBoson";
  
  if(produceRecoils_)
    produces<ParticleCollection>(ZbosonLabel_);
}

mvaPUPPET::~mvaPUPPET(){}

void mvaPUPPET::produce(edm::Event& evt, const edm::EventSetup& es){

  var_.clear();
  
  // calculate Z for recoil
  reco::Particle Z;
  Z.setP4(reco::Candidate::LorentzVector(0, 0, 0, 0));
  
  // take the tau lepton collection
  edm::Handle<pat::TauCollection> tauCollectionHandle;
  evt.getByToken(srcTaus_, tauCollectionHandle);
  const pat::TauCollection tauCollection = *(tauCollectionHandle.product());

  // take mu lepton collection
  edm::Handle<pat::MuonCollection> muCollectionHandle;
  evt.getByToken(srcMuons_, muCollectionHandle);
  const pat::MuonCollection muCollection = *(muCollectionHandle.product());

  
  // vector of packed candidates to store charged and neutral particles belonging to the tau jets
  std::vector<reco::CandidatePtr> chargedTauJetCandidates;
  std::vector<reco::CandidatePtr> neutralTauJetCandidates;

  // loop on the indentified leptons
  float sumEt_Leptons = 0;
	for ( std::vector<edm::EDGetTokenT<reco::CandidateView > >::const_iterator srcLeptons_i = srcLeptons_.begin(); srcLeptons_i != srcLeptons_.end(); ++srcLeptons_i )
	{
		edm::Handle<reco::CandidateView> leptons;		
		evt.getByToken(*srcLeptons_i, leptons);
		for ( reco::CandidateView::const_iterator lepton = leptons->begin();
			lepton != leptons->end(); ++lepton )
		{
				math::PtEtaPhiELorentzVectorD p4Photon;      
				for (auto muon : muCollection)
				{
					if( muon.p4() == lepton->p4())
					{
						p4Photon.SetPt(muon.pfEcalEnergy()/TMath::CosH(muon.p4().eta()));
						p4Photon.SetEta(muon.p4().eta());
						p4Photon.SetPhi(muon.p4().phi());
						p4Photon.SetE(muon.pfEcalEnergy());
					}
				}

				if(p4Photon.E() > 0 )
				{
					Z.setP4(Z.p4()+ lepton->p4()+p4Photon);
					sumEt_Leptons += lepton->p4().Et()+p4Photon.Et();
				}
			else
			{
				Z.setP4(Z.p4()+ lepton->p4());
				sumEt_Leptons += lepton->p4().Et();
			}

			Z.setPdgId(abs(lepton->pdgId()));
			for(std::vector<pat::Tau>::const_iterator tau = tauCollection.begin(); tau!= tauCollection.end(); ++tau)
			{
				if(lepton->p4() != tau->p4()) 				
					continue;
					
					// take the PF constituents of tau lepton
				for(auto candidate : tau->signalCands())
				{
					if(abs(candidate->pdgId()) > 11 and abs(candidate->pdgId()) < 16) continue;
					if(candidate->charge() !=0)
					{
						chargedTauJetCandidates.push_back(candidate);
					}
					else
					{
						neutralTauJetCandidates.push_back(candidate);
					}
				}
			}      
		}
	}

	reco::METCovMatrix rotateToZFrame;
	rotateToZFrame(0,0) = rotateToZFrame(1,1) = std::cos(- Z.p4().Phi());
	rotateToZFrame(0,1) =   std::sin(- Z.p4().Phi());
	rotateToZFrame(1,0) = - std::sin(- Z.p4().Phi());

  // add the Z-boson kinematic information
  // var_["z_pT"]  = Z.pt();
  // var_["z_Phi"] = Z.phi();
  // var_["z_m"]   = Z.mass();

  // and save the Z back to the event 
  if(produceRecoils_){
    std::auto_ptr<ParticleCollection> recoZParticleCollection(new ParticleCollection());
    recoZParticleCollection->push_back(Z);
    evt.put(recoZParticleCollection,ZbosonLabel_);
  }

  // get reference MET and calculate its recoil
  edm::Handle<pat::METCollection> referenceMETs;
  evt.getByToken(referenceMET_, referenceMETs);
  assert((*referenceMETs).size() == 1);
  auto referenceMET = (*referenceMETs)[0];
  reco::Candidate::LorentzVector referenceRecoil;
  if(srcMETFlags_.at(0))
    referenceRecoil = - referenceMET.p4() - Z.p4();
  else
    referenceRecoil = - referenceMET.p4();

  std::string reference = "recoilPFPuppiMet";
  
  // calculate the recoils and save them to MET objects
  int i = 0;
  float sumEt_TauJetCharge  = 0;
  float sumEt_TauJetNeutral = 0;
  std::vector<int>::const_iterator itMETFlags = srcMETFlags_.begin();
  itMETFlags++;
  for ( std::vector<edm::EDGetTokenT<pat::METCollection> >::const_iterator srcMET = srcMETs_.begin(); srcMET != srcMETs_.end() && itMETFlags!=srcMETFlags_.end(); ++srcMET, ++itMETFlags )
  {
    //get inputs
    std::string collection_name = srcMETTags_[i++].label();
    if(debug_) 
      std::cout << "colname: " << collection_name << std::endl;

    edm::Handle<pat::METCollection> MET;
    evt.getByToken(*srcMET, MET);
    assert((*MET).size() == 1 );

    if(debug_)
      std::cout << "tomap: " << collection_name << std::endl;

    reco::Particle tauJetSpouriousComponents;
    tauJetSpouriousComponents.setP4(reco::Candidate::LorentzVector(0, 0, 0, 0));

    if(collection_name.find("ChargedPV") != std::string::npos)
    {
      for( auto particle : neutralTauJetCandidates)
      {
        tauJetSpouriousComponents.setP4(tauJetSpouriousComponents.p4()+particle->p4());            
        sumEt_TauJetNeutral += particle->p4().Et();	
      }
    }    
    else if(collection_name.find("NeutralPV") != std::string::npos)
    {
      for( auto particle : chargedTauJetCandidates)
      {
        tauJetSpouriousComponents.setP4(tauJetSpouriousComponents.p4()+particle->p4());            
        sumEt_TauJetCharge += particle->p4().Et();
      }
    }

    // calculate recoil    
    pat::MET Recoil((*MET)[0]); 

    if((*itMETFlags)){
      if(collection_name.find("ChargedPV") != std::string::npos){
	Recoil.setP4(-Z.p4() +tauJetSpouriousComponents.p4() - (*MET)[0].p4());
        Recoil.setSumEt((*MET)[0].sumEt()-sumEt_TauJetCharge-sumEt_Leptons);    
      }
      else{
	Recoil.setP4(-Z.p4() - (*MET)[0].p4());
        Recoil.setSumEt((*MET)[0].sumEt()-sumEt_Leptons);    
      }
    }
    else{
      if(collection_name.find("NeutralPV") != std::string::npos){
	Recoil.setP4(tauJetSpouriousComponents.p4() - (*MET)[0].p4());      
        Recoil.setSumEt((*MET)[0].sumEt()-sumEt_TauJetNeutral);    
      }
      else{
	Recoil.setP4(-(*MET)[0].p4());		
	Recoil.setSumEt((*MET)[0].sumEt());
      }
    }

    if(produceRecoils_){
      std::auto_ptr<pat::METCollection> patMETRecoilCollection(new pat::METCollection());
      patMETRecoilCollection->push_back(Recoil);
      evt.put(patMETRecoilCollection, "recoil"+collection_name);
    }

    // This only does the PU and PV stuff here
    if (TString(collection_name).Contains(referenceMET_name_) and collection_name != referenceMET_name_){
      TString tempName = Form("%s",collection_name.c_str());
      tempName.ReplaceAll(referenceMET_name_,"");
      collection_name = tempName;
      reco::METCovMatrix rotatedCovMatrix = rotateToZFrame * Recoil.getSignificanceMatrix();
      addToMap(Recoil.p4(), Recoil.sumEt(), collection_name, reference, 1, rotatedCovMatrix);
    }
    else {
      TString tempName = Form("%s",collection_name.c_str());
      tempName.ReplaceAll("slimmedMETs","recoilPF");
      tempName = tempName + "Met";
      collection_name = tempName;
      reco::METCovMatrix rotatedCovMatrix = rotateToZFrame * Recoil.getSignificanceMatrix();
      addToMap(Recoil.p4(), Recoil.sumEt(), "", collection_name, 1, rotatedCovMatrix);
    }
  }

  edm::Handle<pat::JetCollection> jets;
  evt.getByToken(srcJets_, jets);
  size_t jetsSize = jets->size();
  for( size_t iJet = 0; iJet <= 1; iJet++)
  {
    var_["Jet" + std::to_string(iJet)+ "_Pt"]  = (jetsSize > iJet) ? jets->at(iJet).p4().pt() : 0;
    var_["Jet" + std::to_string(iJet)+ "_Eta"] = (jetsSize > iJet) ? jets->at(iJet).p4().Eta() : 0;
    var_["Jet" + std::to_string(iJet)+ "_Phi"] = (jetsSize > iJet) ? jets->at(iJet).p4().Phi() : 0;
    var_["Jet" + std::to_string(iJet)+ "_M"]   = (jetsSize > iJet) ? jets->at(iJet).p4().M() : 0;
  }

  var_["NCleanedJets"] = countJets(*jets, 5);

  // treat other collections and save to map
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(srcVertices_, vertices);
  var_["NVertex"] = countVertices(*vertices);
  

  // print whole map
  for(auto entry : var_){
    if(debug_)
      std::cout << "map " << entry.first << "/" << entry.second << std::endl;
  }
  
  // evaluate phi training and apply angular correction
  Float_t PhiAngle = 0.;
  if(inputFileNamePhiCorrection_.fullPath() != "")
    PhiAngle = GetResponse(mvaReaderPhiCorrection_, variablesForPhiTraining_);

  auto refRecoil = TVector2(referenceRecoil.px(), referenceRecoil.py());
  refRecoil = refRecoil.Rotate(PhiAngle);
  reco::Candidate::LorentzVector phiCorrectedRecoil(refRecoil.Px(), refRecoil.Py(), 0, referenceMET.sumEt());
  // addToMap(phiCorrectedRecoil, referenceMET.sumEt(), "", reference); //, referenceMET.sumEt());
  
  var_["PhiCorrectedRecoil_Phi"] = TVector2::Phi_mpi_pi(refRecoil.Phi());

  // evaluate second training and apply recoil correction
  Float_t RecoilCorrection = 1.0;
  if(inputFileNameRecoilCorrection_.fullPath() != "")
    RecoilCorrection = GetResponse(mvaReaderRecoilCorrection_, variablesForRecoilTraining_);
  refRecoil *= RecoilCorrection;
  
  // calculate new mvaMET
  pat::MET mvaMET(referenceMET);
  reco::Candidate::LorentzVector recoilP4(refRecoil.Px(), refRecoil.Py(), 0, referenceMET.sumEt());
  reco::Candidate::LorentzVector metP4 = - Z.p4() + recoilP4;
  mvaMET.setP4(metP4);
  
  //// save results to event
  std::auto_ptr<pat::METCollection> patMETCollection(new pat::METCollection());
  patMETCollection->push_back(mvaMET);
  evt.put(patMETCollection,"mvaMET");

}

//void mvaPUPPET::addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &name, const std::string &type){
//	addToMap(p4, sumEt, name, type, 1);
//}

void mvaPUPPET::addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &name, const std::string &type, double divisor, reco::METCovMatrix &covMatrix){
  addToMap(p4, sumEt, name, type, divisor);
  if(name == "")
  {
    var_[type + "_cov00" ] = covMatrix(0,0);
    var_[type + "_cov11" ] = covMatrix(1,1);
  }
  else
  {
    var_[type + "_" + name + "_cov00" ] = covMatrix(0,0);
    var_[type + "_" + name + "_cov11" ] = covMatrix(1,1);
  }
}

void mvaPUPPET::addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &name, const std::string &type, double divisor){
  if(name == "")
  {
    var_[type + "_Pt" ] = p4.pt();
    var_[type + "_Phi" ] = p4.phi();
    var_[type + "_sumEt" ] = sumEt/divisor;
  }
  else
  {
    var_[type + "_" + name + "_Pt" ] = p4.pt();
    var_[type + "_" + name + "_Phi" ] = p4.phi();
    var_[type + "_" + name + "_sumEt" ] = sumEt/divisor;
  }
}

unsigned int mvaPUPPET::countVertices(const reco::VertexCollection& vertices){
  return vertices.size();
}

unsigned int mvaPUPPET::countJets(const pat::JetCollection& jets, const float maxPt){
  int nJets = 0;
  for(auto jet : jets){
    if(jet.pt() > maxPt)
      nJets++;
  }
  return nJets;
}

const Float_t mvaPUPPET::GetResponse(const GBRForest * Reader, std::vector<std::string> &variableNames ){

  Float_t * mvaInputVector = createFloatVector(variableNames);
  double result = Reader->GetResponse(mvaInputVector);
  delete mvaInputVector;
  return result;
}


const GBRForest* mvaPUPPET::loadMVAfromFile(const edm::FileInPath& inputFileName, std::vector<std::string>& trainingVariableNames, std::string mvaName){
  
  if ( inputFileName.location()==edm::FileInPath::Unknown ) 
    throw cms::Exception("PFMETAlgorithmMVA::loadMVA") << " Failed to find File = " << inputFileName << " !!\n";
  
  TFile* inputFile = new TFile(inputFileName.fullPath().data());
  
  std::vector<std::string> *lVec = (std::vector<std::string>*)inputFile->Get("varlist");

  for(unsigned int i=0; i< lVec->size();++i)
  {
    trainingVariableNames.push_back(lVec->at(i));
    if(debug_) 
      std::cout << "training variable " << i << ":" << lVec->at(i) << std::endl;
  }
  const GBRForest* mva = (GBRForest*)inputFile->Get(mvaName.data());
  if ( !mva )
    throw cms::Exception("PFMETAlgorithmMVA::loadMVA") << " Failed to load MVA from file = " << inputFileName.fullPath().data() << " !!\n";
  
  delete inputFile;
  
  return mva;
}

Float_t* mvaPUPPET::createFloatVector(std::vector<std::string> variableNames){
  Float_t* floatVector = new Float_t[variableNames.size()];
  std::cout << "creating Float Vector: " << std::endl;
  for(size_t i = 0; i < variableNames.size(); ++i){
      floatVector[i] = var_[variableNames[i]];
      std::cout << variableNames[i] << " = " << floatVector[i] << std::endl;
  }
  std::cout << "--------" << std::endl;
  return floatVector;
}


DEFINE_FWK_MODULE(mvaPUPPET);
