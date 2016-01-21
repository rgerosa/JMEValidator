#include <iostream>

#include "JMEAnalysis/JMEValidator/plugins/MVAMET.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/Particle.h"

typedef std::vector<reco::Particle> ParticleCollection;

MVAMET::MVAMET(const edm::ParameterSet& cfg){

  // get MET that the mva is applied on
  referenceMET_      = consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("referenceMET"));

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
    throw cms::Exception("MVAMET::MVAMET") << " Failed to load MET flags   !!\n";

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
  edm::FileInPath weightFile; 
  if(cfg.existsAs<edm::FileInPath>("weightFile"))
    weightFile = cfg.getParameter<edm::FileInPath>("weightFile");
  mvaReaderPhiCorrection_     = loadMVAfromFile(weightFile, variablesForPhiTraining_, "PhiCorrectedRecoil"); 
  mvaReaderRecoilCorrection_  = loadMVAfromFile(weightFile, variablesForRecoilTraining_, "LongZCorrectedRecoil");
  mvaReaderCovU1_             = loadMVAfromFile(weightFile, variablesForCovU1_, "CovU1");
  mvaReaderCovU2_             = loadMVAfromFile(weightFile, variablesForCovU2_, "CovU2");

  // prepare for saving the final mvaMET to the event
  if(cfg.existsAs<std::string>("MVAMETLabel"))
    mvaMETLabel_ = cfg.getParameter<std::string>("MVAMETLabel");
  else
    mvaMETLabel_ = "MVAMET";

  produces<pat::METCollection>(mvaMETLabel_);
  if(produceRecoils_)
      produces<pat::METCollection>( "recoil"+mvaMETLabel_ );

  // produce Zboson tag
  if(cfg.existsAs<std::string>("ZbosonLabel"))
    ZbosonLabel_ = cfg.getParameter<std::string>("ZbosonLabel");
  else
    ZbosonLabel_ = "ZtagBoson";
  
  if(produceRecoils_)
    produces<ParticleCollection>(ZbosonLabel_);
}

MVAMET::~MVAMET(){}

void MVAMET::calculateRecoil(edm::Handle<pat::METCollection> MET, reco::Particle Z, reco::Particle tauJetSpouriousComponents, float sumEt_TauJetCharge, float sumEt_TauJetNeutral, float sumEt_Leptons , int METFlag, edm::Event& evt, std::string collection_name, float divisor)
{

    reco::METCovMatrix rotateToZFrame;
    rotateToZFrame(0,0) = rotateToZFrame(1,1) = std::cos(- Z.p4().Phi());
    rotateToZFrame(0,1) =   std::sin(- Z.p4().Phi());
    rotateToZFrame(1,0) = - std::sin(- Z.p4().Phi());

    pat::MET Recoil((*MET)[0]); 

    Recoil.setP4(tauJetSpouriousComponents.p4() - (*MET)[0].p4());
    Recoil.setSumEt((*MET)[0].sumEt());    
    // subtract Z-Boson if contained in MET to get recoil
    if (!(METFlag==2))
    {
      Recoil.setP4(Recoil.p4()-Z.p4());
      Recoil.setSumEt(Recoil.sumEt()-sumEt_TauJetCharge-sumEt_Leptons);    
    }
 
    if (!(METFlag==1))
    {
      Recoil.setSumEt(Recoil.sumEt()-sumEt_TauJetNeutral);    
    }

    reco::METCovMatrix rotatedCovMatrix = rotateToZFrame * Recoil.getSignificanceMatrix();
    Recoil.setSignificanceMatrix( rotatedCovMatrix );

    if(produceRecoils_){
      std::auto_ptr<pat::METCollection> patMETRecoilCollection(new pat::METCollection());
      patMETRecoilCollection->push_back(Recoil);
      evt.put(patMETRecoilCollection, "recoil"+collection_name);
    }

    addToMap(Recoil.p4(), Recoil.sumEt(), "recoil"+ collection_name, divisor, rotatedCovMatrix);

}

void MVAMET::produce(edm::Event& evt, const edm::EventSetup& es){

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
  pat::MET referenceRecoil(referenceMET);
  if ( srcMETFlags_.at(0) == 2 )
  {
    referenceRecoil.setP4(- referenceMET.p4());
  }
  else
  {
    referenceRecoil.setP4( - referenceMET.p4() - Z.p4() );
    referenceRecoil.setSumEt( referenceMET.sumEt() - sumEt_Leptons);
  }
  
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
    // charged components
    if ( ((*itMETFlags)==0) or ((*itMETFlags)==1) )
    {
      for( auto particle : chargedTauJetCandidates)
      {
        tauJetSpouriousComponents.setP4(tauJetSpouriousComponents.p4()+particle->p4());            
        sumEt_TauJetCharge += particle->p4().Et();
      }
    }
    // neutral components
    if ( ((*itMETFlags)==0) or ((*itMETFlags)==2) )
    {
      for( auto particle : neutralTauJetCandidates)
      {
        tauJetSpouriousComponents.setP4(tauJetSpouriousComponents.p4()+particle->p4());            
        sumEt_TauJetNeutral += particle->p4().Et();	
      }
    }    

    // calculate recoil   
    calculateRecoil(MET, Z, tauJetSpouriousComponents, sumEt_TauJetCharge, sumEt_TauJetNeutral, sumEt_Leptons, (*itMETFlags), evt, collection_name, referenceRecoil.sumEt());
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

  auto refRecoil = TVector2(referenceRecoil.p4().px(), referenceRecoil.p4().py());
  refRecoil = refRecoil.Rotate(PhiAngle);
  reco::Candidate::LorentzVector phiCorrectedRecoil(refRecoil.Px(), refRecoil.Py(), 0, referenceMET.sumEt());
  addToMap(phiCorrectedRecoil, referenceMET.sumEt(), "PhiCorrectedRecoil", 1); //, referenceMET.sumEt());
  
  var_["PhiCorrectedRecoil_Phi"] = TVector2::Phi_mpi_pi(refRecoil.Phi());

  // evaluate second training and apply recoil correction
  Float_t RecoilCorrection = 1.0;
  RecoilCorrection = GetResponse(mvaReaderRecoilCorrection_, variablesForRecoilTraining_);
  refRecoil *= RecoilCorrection;


  // create pat::MET from recoil
  pat::MET recoilmvaMET(referenceMET);
  reco::Candidate::LorentzVector recoilP4(refRecoil.Px(), refRecoil.Py(), 0, referenceMET.sumEt());
  recoilmvaMET.setP4(recoilP4);
  addToMap(recoilP4, referenceMET.sumEt(), "LongZCorrectedRecoil", 1); //, referenceMET.sumEt());

  // evaluate covariance matrix regression
  Float_t CovU1 = 0;
  Float_t CovU2 = 0;
  CovU1 = GetResponse(mvaReaderCovU1_, variablesForCovU1_) * refRecoil.Mod();
  CovU2 = GetResponse(mvaReaderCovU2_, variablesForCovU2_) * refRecoil.Mod();
  reco::METCovMatrix recoilmvaMETCov;
  recoilmvaMETCov(0, 0) =  std::pow(CovU1, 2);
  recoilmvaMETCov(0, 1) =  recoilmvaMETCov(1, 0) = 0;
  recoilmvaMETCov(1, 1) = std::pow(CovU2, 2);
  recoilmvaMET.setSignificanceMatrix(recoilmvaMETCov);
  
  //// save results to event
  std::auto_ptr<pat::METCollection> recoilpatMETCollection(new pat::METCollection());
  recoilpatMETCollection->push_back(recoilmvaMET);
  evt.put(recoilpatMETCollection,"recoilMVAMET");

  // calculate new mvaMET
  pat::MET mvaMET(referenceMET);
  reco::Candidate::LorentzVector metP4 = - Z.p4() + recoilP4;
  mvaMET.setP4(metP4);

  reco::METCovMatrix mvaMETCov;
  double cosPhi =  std::cos(refRecoil.Phi());
  double sinPhi =  std::sin(refRecoil.Phi());
  mvaMETCov(0, 0) =  std::pow(CovU1, 2)*cosPhi*cosPhi + std::pow(CovU2, 2)*sinPhi*sinPhi;
  mvaMETCov(0, 1) = -std::pow(CovU1, 2) * sinPhi*cosPhi + std::pow(CovU2, 2) * sinPhi*cosPhi;
  mvaMETCov(1, 0) =  mvaMETCov(0, 1);
  mvaMETCov(1, 1) =  std::pow(CovU1, 2) * sinPhi*sinPhi + std::pow(CovU2, 2) * cosPhi*cosPhi;
  mvaMET.setSignificanceMatrix(mvaMETCov);

  //// save results to event
  std::auto_ptr<pat::METCollection> patMETCollection(new pat::METCollection());
  patMETCollection->push_back(mvaMET);
  evt.put(patMETCollection,"MVAMET");

}

void MVAMET::addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &type, double divisor, reco::METCovMatrix &covMatrix)
{
  addToMap(p4, sumEt, type, divisor);
  var_[type +  "_Cov00" ] = covMatrix(0,0);
  var_[type +  "_Cov11" ] = covMatrix(1,1);
}

void MVAMET::addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &type, double divisor)
{
  var_[type + "_Pt" ] = p4.pt();
  var_[type + "_Phi" ] = p4.phi();
  var_[type + "_sumEt" ] = sumEt;
  var_[type + "_sumEtFraction" ] = sumEt/divisor;
}

unsigned int MVAMET::countVertices(const reco::VertexCollection& vertices){
  return vertices.size();
}

unsigned int MVAMET::countJets(const pat::JetCollection& jets, const float maxPt){
  int nJets = 0;
  for(auto jet : jets){
    if(jet.pt() > maxPt)
      nJets++;
  }
  return nJets;
}

const Float_t MVAMET::GetResponse(const GBRForest * Reader, std::vector<std::string> &variableNames ){

  Float_t * mvaInputVector = createFloatVector(variableNames);
  double result = Reader->GetResponse(mvaInputVector);
  delete mvaInputVector;
  return result;
}


const GBRForest* MVAMET::loadMVAfromFile(const edm::FileInPath& inputFileName, std::vector<std::string>& trainingVariableNames, std::string mvaName){
  
  if ( inputFileName.location()==edm::FileInPath::Unknown ) 
    throw cms::Exception("PFMETAlgorithmMVA::loadMVA") << " Failed to find File = " << inputFileName << " !!\n";
  
  TFile* inputFile = new TFile(inputFileName.fullPath().data());
  std::string variableListName = mvaName + "varlist";
  std::vector<std::string> *lVec = (std::vector<std::string>*)inputFile->Get(variableListName.c_str());

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

Float_t* MVAMET::createFloatVector(std::vector<std::string> variableNames){
  Float_t* floatVector = new Float_t[variableNames.size()];
  for(size_t i = 0; i < variableNames.size(); ++i){
      floatVector[i] = var_[variableNames[i]];
  }
  return floatVector;
}


DEFINE_FWK_MODULE(MVAMET);
