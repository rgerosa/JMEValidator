#include "JMEAnalysis/JMEValidator/plugins/mvaPUPPET.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/Particle.h"

typedef std::vector<reco::Particle> ParticleCollection;

mvaPUPPET::mvaPUPPET(const edm::ParameterSet& cfg){

  // get tokens for input METs and prepare for saving the corresponding recoils to the event
  srcMETTags_   = cfg.getParameter<vInputTag>("srcMETs");
  for(vInputTag::const_iterator it=srcMETTags_.begin();it!=srcMETTags_.end();it++) {
    srcMETs_.push_back( consumes<pat::METCollection >( *it ) );
    produces<pat::METCollection>( "recoil"+(*it).label() );
  }
  
  
  // get debug flags
  if (cfg.existsAs<bool>("debug"))
    debug_ = cfg.getParameter<bool>("debug");
  else
    debug_ = false;
  
  // get MET that the mva is applied on
  referenceMET_ = consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("referenceMET"));
  referenceMET_name_ = cfg.getParameter<edm::InputTag>("referenceMET").label();
  for(auto itag : srcMETTags_){
    if(debug_)
      std::cout << "itag: " << itag.label() << std::endl;
  }
  
  // get vertices
  srcVertices_ = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("srcVertices"));

  // get jets
  srcJets_     = consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("srcJets"));
  
  //get leptons to calculate Z vector and save it as a reco::candidate back to the event
  vInputTag srcLeptonsTags = cfg.getParameter<vInputTag>("srcLeptons");
  for(vInputTag::const_iterator it=srcLeptonsTags.begin();it!=srcLeptonsTags.end();it++) {
    srcLeptons_.push_back( consumes<reco::CandidateView >( *it ) );
  }
  
  // load weight files
  edm::ParameterSet cfgInputFileNames = cfg.getParameter<edm::ParameterSet>("inputFileNames");
  
  if(cfgInputFileNames.existsAs<edm::FileInPath>("PhiCorrectionWeightFile")){
    inputFileNamePhiCorrection_ = cfgInputFileNames.getParameter<edm::FileInPath>("PhiCorrectionWeightFile");
    mvaReaderPhiCorrection_     = loadMVAfromFile(inputFileNamePhiCorrection_, variablesForPhiTraining_, "PhiCor");
  }
  
  if(cfgInputFileNames.existsAs<edm::FileInPath>("RecoilCorrectionWeightFile")){
    inputFileNameRecoilCorrection_ = cfgInputFileNames.getParameter<edm::FileInPath>("RecoilCorrectionWeightFile");
    mvaReaderRecoilCorrection_  = loadMVAfromFile(inputFileNameRecoilCorrection_, variablesForRecoilTraining_, "RecoilCor");
  }
  
  // prepare for saving the final mvaPUPPET to the event
  produces<pat::METCollection>("mvaMET");
  // produce Zboson tag
  produces<ParticleCollection>("ZtagBoson");
}

mvaPUPPET::~mvaPUPPET(){}

void mvaPUPPET::produce(edm::Event& evt, const edm::EventSetup& es){

	var_.clear();

	// calculate Z for recoil
	//reco::Particle Z;
	reco::Particle Z;
	Z.setP4(reco::Candidate::LorentzVector(0, 0, 0, 0));
	for ( std::vector<edm::EDGetTokenT<reco::CandidateView > >::const_iterator srcLeptons_i = srcLeptons_.begin();
		srcLeptons_i != srcLeptons_.end(); ++srcLeptons_i ){

	  edm::Handle<reco::CandidateView> leptons;		
	  evt.getByToken(*srcLeptons_i, leptons);
	  for ( reco::CandidateView::const_iterator lepton = leptons->begin();
		lepton != leptons->end(); ++lepton ){
	    Z.setP4(Z.p4()+ lepton->p4());
	    Z.setPdgId(abs(lepton->pdgId()));
	  }
	}
	var_["z_pT"] = Z.pt();
	var_["z_Phi"] = Z.phi();
	var_["z_m"] = Z.mass();

	// and save the Z back to the event 
	std::auto_ptr<ParticleCollection> recoZParticleCollection(new ParticleCollection());
	recoZParticleCollection->push_back(Z);
	evt.put(recoZParticleCollection, "ZtagBoson");

	// get reference MET and calculate its recoil
	edm::Handle<pat::METCollection> referenceMETs;
	evt.getByToken(referenceMET_, referenceMETs);
	assert((*referenceMETs).size() == 1);
	auto referenceMET = (*referenceMETs)[0];
	reco::Candidate::LorentzVector referenceNegRecoil = - Z.p4() - referenceMET.p4();
	std::string reference = "referenceNegRecoil";
	addToMap(referenceNegRecoil, referenceMET.sumEt(), reference, referenceMET_name_);

	// calculate the recoils and save them to MET objects
	int i = 0;
	for ( std::vector<edm::EDGetTokenT<pat::METCollection> >::const_iterator srcMET = srcMETs_.begin();
		srcMET != srcMETs_.end(); ++srcMET ){

		//get inputs
		std::string collection_name = srcMETTags_[i++].label();
		if(debug_) 
		  std::cout << "colname: " << collection_name << std::endl;
		edm::Handle<pat::METCollection> MET;
		evt.getByToken(*srcMET, MET);
		assert((*MET).size() == 1 );
		std::string string_input = "negRecoil"; 
		if(debug_)
		  std::cout << "tomap: " << collection_name << std::endl;

		// calculate recoil
		pat::MET negRecoil((*MET)[0]); 
		negRecoil.setP4(-Z.p4() - (*MET)[0].p4());

		addToMap(negRecoil.p4(), (*MET)[0].sumEt(), string_input, collection_name, referenceMET.sumEt());
		std::auto_ptr<pat::METCollection> patMETRecoilCollection(new pat::METCollection());
		patMETRecoilCollection->push_back(negRecoil);
		evt.put(patMETRecoilCollection, "recoil"+collection_name);
	}

	// print whole map
	for(auto entry : var_){
	  if(debug_)
	    std::cout << "map" << entry.first << "/" << entry.second << std::endl;
	}
	
	// treat other collections and save to map
	edm::Handle<reco::VertexCollection> vertices;
	evt.getByToken(srcVertices_, vertices);
	var_["nVertices"] = countVertices(*vertices);

	edm::Handle<pat::JetCollection> jets;
	evt.getByToken(srcJets_, jets);
	var_["nJets"] = countJets(*jets, 30);

	// evaluate phi training and apply angular correction
	Float_t PhiAngle = 0.;
	if(inputFileNamePhiCorrection_.fullPath() != "")
	  PhiAngle = GetResponse(mvaReaderPhiCorrection_, variablesForPhiTraining_);

	
	auto refNegRecoil = TVector2(referenceNegRecoil.px(), referenceNegRecoil.py());
	refNegRecoil = refNegRecoil.Rotate(PhiAngle);
	reco::Candidate::LorentzVector phiCorrectedNegRecoil(refNegRecoil.Px(), refNegRecoil.Py(), 0, referenceMET.sumEt());
	addToMap(phiCorrectedNegRecoil, referenceMET.sumEt(), reference, referenceMET_name_);

	// evaluate second training and apply recoil correction
	Float_t RecoilCorrection = 1.0;
	if(inputFileNameRecoilCorrection_.fullPath() != "")
	 RecoilCorrection = GetResponse(mvaReaderRecoilCorrection_, variablesForRecoilTraining_);
	refNegRecoil *= RecoilCorrection;
	
	// calculate new mvaPUPPET
	pat::MET mvaMET(referenceMET);
	reco::Candidate::LorentzVector recoilP4(refNegRecoil.Px(), refNegRecoil.Py(), 0, referenceMET.sumEt());
	reco::Candidate::LorentzVector metP4 = - Z.p4() - recoilP4;
	mvaMET.setP4(metP4);

	//// save results to event
	// mvaPUPPET
	std::auto_ptr<pat::METCollection> patMETCollection(new pat::METCollection());
	patMETCollection->push_back(mvaMET);
	evt.put(patMETCollection,"mvaMET");
	
}

void mvaPUPPET::addToMap(reco::Candidate::LorentzVector p4, double sumEt, std::string &name, std::string &type)
{
	addToMap(p4, sumEt, name, type, 1);
}

void mvaPUPPET::addToMap(reco::Candidate::LorentzVector p4, double sumEt, std::string &name, std::string &type, double divisor)
{
  if(debug_)
    std::cout << "hier" << name << ", " << type << std::endl;

  var_[type + "_" + name + "_pt" ] = p4.pt();
  var_[type + "_" + name + "_phi" ] = p4.phi();
  var_[type + "_" + name + "_sumEt" ] = sumEt/divisor;
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


const GBRForest* mvaPUPPET::loadMVAfromFile(const edm::FileInPath& inputFileName, std::vector<std::string>& trainingVariableNames, std::string mvaName)
{
  if ( inputFileName.location()==edm::FileInPath::Unknown ) 
    throw cms::Exception("PFMETAlgorithmMVA::loadMVA") 
      << " Failed to find File = " << inputFileName << " !!\n";

  TFile* inputFile = new TFile(inputFileName.fullPath().data());
  
  std::vector<std::string> *lVec = (std::vector<std::string>*)inputFile->Get("varlist");
  for(unsigned int i=0; i< lVec->size();++i)
    trainingVariableNames.push_back(lVec->at(i));
  const GBRForest* mva = (GBRForest*)inputFile->Get(mvaName.data());
  if ( !mva )
    throw cms::Exception("PFMETAlgorithmMVA::loadMVA")
      << " Failed to load MVA from file = " << inputFileName.fullPath().data() << " !!\n";
  
  delete inputFile;
  
  return mva;
}

Float_t* mvaPUPPET::createFloatVector(std::vector<std::string> variableNames)
{
    Float_t* floatVector = new Float_t[variableNames.size()];
    for(size_t i = 0; i < variableNames.size(); ++i)
    {
        floatVector[i] = var_[variableNames[i]];
    }
    return floatVector;
}


DEFINE_FWK_MODULE(mvaPUPPET);
