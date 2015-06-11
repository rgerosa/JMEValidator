#include "JMEAnalysis/JMEValidator/plugins/mvaPUPPET.h"
#include "FWCore/Framework/interface/MakerMacros.h"



mvaPUPPET::mvaPUPPET(const edm::ParameterSet& cfg)
{
	std::cout << "hallo welt" << std::endl;
	srcMETTags_   = cfg.getParameter<vInputTag>("srcMETs");
	referenceMET_name_ = cfg.getParameter<edm::InputTag>("referenceMET").label();
	for(vInputTag::const_iterator it=srcMETTags_.begin();it!=srcMETTags_.end();it++) {
		srcMETs_.push_back( consumes<pat::METCollection >( *it ) );
	}

	referenceMET_ = consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("referenceMET"));
	for(auto itag : srcMETTags_)
		std::cout << "itag: " << itag.label() << std::endl;
	srcVertices_ = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("srcVertices"));
	srcJets_     = consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("srcJets"));

	edm::ParameterSet cfgInputFileNames = cfg.getParameter<edm::ParameterSet>("inputFileNames");

//	edm::FileInPath inputFileNamePhiCorrection = cfgInputFileNames.getParameter<edm::FileInPath>("PhiCorrectionWeightFile");
//	mvaReaderPhiCorrection_     = loadMVAfromFile(inputFileNamePhiCorrection, variablesForPhiTraining_);

//	edm::FileInPath inputFileNameRecoilCorrection = cfgInputFileNames.getParameter<edm::FileInPath>("RecoilCorrectionWeightFile");
//	mvaReaderRecoilCorrection_  = loadMVAfromFile(inputFileNameRecoilCorrection, variablesForRecoilTraining_);
	std::cout << "init done" << std::endl;
}

mvaPUPPET::~mvaPUPPET()
{
	std::cout << "destructor" << std::endl;
}


void mvaPUPPET::produce(edm::Event& evt, const edm::EventSetup& es)
{
	std::cout << "produce" << std::endl;
	var_.clear();


	edm::Handle<pat::METCollection> referenceMETs;
	evt.getByToken(referenceMET_, referenceMETs);
	assert((*referenceMETs).size() == 1);
	auto referenceMET = (*referenceMETs)[0];
	std::string reference = "reference";
	addToMap(referenceMET, reference, referenceMET_name_);

	std::cout << "map" << std::endl;

	int i = 0;
	for ( std::vector<edm::EDGetTokenT<pat::METCollection> >::const_iterator srcMET = srcMETs_.begin();
		srcMET != srcMETs_.end(); ++srcMET )
	{
		std::string collection_name = srcMETTags_[i].label();
		std::cout << "colname: " << collection_name << std::endl;
		edm::Handle<pat::METCollection> MET;
		evt.getByToken(*srcMET, MET);
		assert((*MET).size() == 1 );
		std::string string_input = "input"; 
		std::cout << "tomap: " << collection_name << std::endl;
		for(auto met: (*MET))
			addToMap(met, string_input, collection_name, referenceMET.sumEt());
		++i;
	}
	
	for(auto entry : var_)
		std::cout << "map" << entry.first << "/" << entry.second << std::endl;


	edm::Handle<reco::VertexCollection> vertices;
	evt.getByToken(srcVertices_, vertices);
	var_["nVertices"] = countVertices(*vertices);

	edm::Handle<pat::JetCollection> jets;
	evt.getByToken(srcJets_, jets);
	var_["nJets"] = countJets(*jets, 30);


	// evaluate training
	Float_t PhiAngle = GetResponse(mvaReaderPhiCorrection_, variablesForPhiTraining_);

	//calculate phi corrected MET and store it in new pat::MET
	auto metDir = TVector2(referenceMET.px(), referenceMET.py());
	metDir.Rotate(PhiAngle);

	pat::MET phiCorrectedMET(referenceMET);
	reco::Candidate::LorentzVector lv(metDir.Px(), metDir.Py(), 0, referenceMET.sumEt());
	phiCorrectedMET.setP4(lv);

	std::cout << "dphi: " << referenceMET.phi() - phiCorrectedMET.phi() << std::endl;

	Float_t result2 = GetResponse(mvaReaderRecoilCorrection_, variablesForRecoilTraining_);
	std::cout << result2 << std::endl;
	// evaluieren des GBR Trees
	// funktion : string zu float Wert aus PAT::MET bzw. andere grÃsse
	// Name: input_collection_value, input_ak4pfjets_
	// oder  phiCor_collection_value
	// oder  mva_collection_value
	// final: patPFMVAMetPuppi
	// korrektur anwenden, nochmal evaluieren
	// anwenden
	// neues MET Objekt erstellen und als Collection event speichern

}

void mvaPUPPET::addToMap(pat::MET &met, std::string &name, std::string &type)
{
	addToMap(met, name, type, 1);
}

void mvaPUPPET::addToMap(pat::MET &met, std::string &name, std::string &type, double divisor)
{
	std::cout << "hier" << name << ", " << type << std::endl;
	var_[type + "_" + name + "_pt" ] = met.pt();
	var_[type + "_" + name + "_phi" ] = met.phi();
	var_[type + "_" + name + "_sumEt" ] = met.sumEt()/divisor;
	std::cout << "end" << std::endl;
}

unsigned int mvaPUPPET::countVertices(const reco::VertexCollection& vertices)
{
	return vertices.size();
}

unsigned int mvaPUPPET::countJets(const pat::JetCollection& jets, const float maxPt)
{
	int nJets = 0;
	for(auto jet : jets)
	{
		if(jet.pt() > maxPt)
			nJets++;
	}
	return nJets;
}

const Float_t mvaPUPPET::GetResponse(const GBRForest * Reader, std::vector<std::string> &variableNames )
{
	return 1.0; // return dummy until weightfile is there
    Float_t * mvaInputVector = createFloatVector(variableNames);
    double result = Reader->GetResponse(mvaInputVector);
    delete mvaInputVector;
    return result;
}


const GBRForest* mvaPUPPET::loadMVAfromFile(const edm::FileInPath& inputFileName, std::vector<std::string>& trainingVariableNames)
{
	if ( inputFileName.location()==edm::FileInPath::Unknown ) 
	throw cms::Exception("PFMETAlgorithmMVA::loadMVA") 
	<< " Failed to find File = " << inputFileName << " !!\n";
	TFile* inputFile = new TFile(inputFileName.fullPath().data());

	std::vector<std::string> *lVec = (std::vector<std::string>*)inputFile->Get("varlist");
	for(unsigned int i=0; i< lVec->size();++i)
		trainingVariableNames.push_back(lVec->at(i));

	const GBRForest* mva = (GBRForest*)inputFile->Get("Forest");
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
