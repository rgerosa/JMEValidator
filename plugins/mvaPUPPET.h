#ifndef JMEValidator_mvaPUPPET_h
#define JMEValidator_mvaPUPPET_h

/** \class mvaPUPPET
 *  *
 *  * Apply recoil corrections and improve PUPPI missing Et (a.k.a. PUPPET) 
 *  *
 * */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include <iostream>
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <TFile.h>
#include <TVector2.h>
#include <TNtuple.h>

class mvaPUPPET : public edm::stream::EDProducer<>
{
public:
	mvaPUPPET(const edm::ParameterSet&);
	~mvaPUPPET();



private:
	void produce(edm::Event&, const edm::EventSetup&);
	std::map<std::string, Float_t> var_;
	edm::EDGetTokenT<reco::VertexCollection> srcVertices_;
	typedef std::vector<edm::InputTag> vInputTag;
	vInputTag srcMETTags_;
	std::vector<edm::EDGetTokenT<pat::METCollection > > srcMETs_;
	edm::EDGetTokenT<pat::METCollection> referenceMET_;
	std::string referenceMET_name_;
	edm::EDGetTokenT<pat::JetCollection> srcJets_;
	std::vector<edm::EDGetTokenT<reco::CandidateView > > srcLeptons_;


	unsigned int countVertices(const reco::VertexCollection& vertices);
	unsigned int countJets(const pat::JetCollection& jets, const float maxPt);
	Float_t* createFloatVector(std::vector<std::string> variableNames);

	std::vector<std::string> variablesForPhiTraining_ = {"nVertices", "nJets" };
	std::vector<std::string> variablesForRecoilTraining_ = {"nVertices", "nJets" };

	const GBRForest* loadMVAfromFile(const edm::FileInPath& inputFileName, std::vector<std::string>& trainingVariableNames);
	const Float_t GetResponse(const GBRForest * Reader, std::vector<std::string> &variableNames );
	void addToMap(reco::Candidate::LorentzVector p4, double sumEt, std::string &name, std::string &type);
	void addToMap(reco::Candidate::LorentzVector p4, double sumEt, std::string &name, std::string &type, double divisor);

	const GBRForest* mvaReaderPhiCorrection_;
	const GBRForest* mvaReaderRecoilCorrection_;

	bool writeNtuple_ = true;
	std::vector<std::string> varForSkim_;
	TFile* skimOutputFile_;
	TNtuple* skimNTuple_;
	void writeSkim();
}; 
#endif
