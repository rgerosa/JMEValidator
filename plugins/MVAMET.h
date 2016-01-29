#ifndef JMEValidator_MVAMET_h
#define JMEValidator_MVAMET_h

/** \class MVAMET
 *  * Apply recoil corrections and improve PUPPI missing Et (a.k.a. PUPPET) 
 * */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include <iostream>
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <TFile.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TMath.h>

class metPlus : public pat::MET {
public: 
  float sumEt_TauJetCharge  = 0;
  float sumEt_TauJetNeutral = 0;
  int METFlag;
  std::string collection_name;
  reco::Particle tauJetSpouriousComponents;
  // datentyp: sumEt_TauJetCharge, TauJetNeutral, metflat, collectionname
};

class recoilingBoson : public reco::Particle {
  public: 
  recoilingBoson()
  {
    sumEt_Leptons = 0;
  }
  std::vector<reco::CandidatePtr> chargedTauJetCandidates;
  std::vector<reco::CandidatePtr> neutralTauJetCandidates;
  float sumEt_Leptons;
  std::vector<int> pdgIds;
  std::vector<const reco::Candidate*> leptons;
};

class MVAMET : public edm::stream::EDProducer<> {

 public:

  // basic constructor from parameter set
  MVAMET(const edm::ParameterSet&);
  ~MVAMET();
  
  void produce(edm::Event&, const edm::EventSetup&);
  typedef std::vector<edm::InputTag> vInputTag;

  // create a vector given input variables
  Float_t* createFloatVector(std::vector<std::string> variableNames);

  unsigned int countVertices(const reco::VertexCollection& vertices);
  unsigned int countJets(const pat::JetCollection& jets, const float maxPt);

  // load MVA file produced in the training
  const GBRForest* loadMVAfromFile(const edm::FileInPath& inputFileName, std::vector<std::string>& trainingVariableNames, std::string mvaName);
  // read the response
  const Float_t GetResponse(const GBRForest * Reader,std::vector<std::string> &variableNames );

  // to correctly create the map of regression input vriables
  //void addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &name, const std::string &type);
  void addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &type, double divisor);
  void addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &type, double divisor, reco::METCovMatrix &covMatrix);


  void calculateRecoil(metPlus* MET, recoilingBoson &Z, edm::Event& evt, float divisor);
private:
  void doCombinations(int offset, int k);
  void saveMap(edm::Event& evt);
  void calculateRecoilingObjects(edm::Event& evt, const pat::MuonCollection&, const pat::TauCollection& );
  void cleanLeptonsFromSS();
  void handleMuons(const reco::Candidate* lepton, recoilingBoson& Z, const pat::MuonCollection& );
  void handleTaus(const reco::Candidate* lepton, recoilingBoson& Z, const pat::TauCollection& );
  void fillEventInformation(edm::Event&);
  std::string mvaMETLabel_;
  std::string ZbosonLabel_;


  vInputTag srcMETTags_;
  
  std::vector<edm::EDGetTokenT<pat::METCollection > >  srcMETs_;
  edm::EDGetTokenT<pat::METCollection>                 referenceMET_;
  edm::EDGetTokenT<reco::VertexCollection>             srcVertices_;
  edm::EDGetTokenT<pat::JetCollection>                 srcJets_;
  std::vector<edm::EDGetTokenT<reco::CandidateView > > srcLeptons_;
  edm::EDGetTokenT<pat::TauCollection>                 srcTaus_;
  edm::EDGetTokenT<pat::MuonCollection>                srcMuons_;
//  edm::EDGetTokenT<pat::ElectronCollection>            srcElectrons_;
  
  std::string referenceMET_name_;
  
  std::vector<int> srcMETFlags_;
  
  std::map<std::string, Float_t> var_;
  
  edm::FileInPath inputFileNamePhiCorrection_;
  edm::FileInPath inputFileNameRecoilCorrection_;
 
  std::vector<const reco::Candidate*> allLeptons_;
  std::vector<std::vector<const reco::Candidate*>> combinations_;
  std::vector<const reco::Candidate*> combination_;
  
  std::vector<std::string> variablesForPhiTraining_  = {};
  std::vector<std::string> variablesForRecoilTraining_  = {};
  std::vector<std::string> variablesForCovU1_  = {};
  std::vector<std::string> variablesForCovU2_  = {};
  edm::Ptr<reco::Candidate> leptonPtr_;

  const GBRForest* mvaReaderPhiCorrection_;
  const GBRForest* mvaReaderRecoilCorrection_;
  const GBRForest* mvaReaderCovU1_;
  const GBRForest* mvaReaderCovU2_;

  bool debug_;
  bool saveMap_;
  bool produceRecoils_;
  std::vector<recoilingBoson> Bosons_;
}; 
#endif
