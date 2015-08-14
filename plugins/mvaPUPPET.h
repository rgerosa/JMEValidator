#ifndef JMEValidator_mvaPUPPET_h
#define JMEValidator_mvaPUPPET_h

/** \class mvaPUPPET
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

class mvaPUPPET : public edm::stream::EDProducer<> {

 public:

  // basic constructor from parameter set
  mvaPUPPET(const edm::ParameterSet&);
  ~mvaPUPPET();
  
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
  void addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &name, const std::string &type);
  void addToMap(reco::Candidate::LorentzVector p4, double sumEt, const std::string &name, const std::string &type, double divisor);


private:

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
  
  edm::InputTag srcPuppiWeights_;
  edm::EDGetTokenT<edm::ValueMap<float> > puppiWeights_;
  
  std::string referenceMET_name_;
  
  std::vector<int> srcMETFlags_;
  
  std::map<std::string, Float_t> var_;
  
  edm::FileInPath inputFileNamePhiCorrection_;
  edm::FileInPath inputFileNameRecoilCorrection_;
  
  
  std::vector<std::string> variablesForPhiTraining_  = {};/*    = {"recoilPFPuppiMet_Pt", */
	/* 						  "recoilPFPuppiMet_Phi", */
	/* 						  "recoilPFPuppiMet_sumEt", */
	/* 						  "recoilPFPuppiMet_ChargedPU_Pt", */
	/* 						  "recoilPFPuppiMet_ChargedPU_Phi", */
	/* 						  "recoilPFPuppiMet_ChargedPU_sumEt", */
	/* 						  "recoilPFPuppiMet_ChargedPV_Pt", */
	/* 						  "recoilPFPuppiMet_ChargedPV_Phi", */
	/* 						  "recoilPFPuppiMet_ChargedPV_sumEt", */
	/* 						  "recoilPFPuppiMet_NeutralPV_Pt", */
	/* 						  "recoilPFPuppiMet_NeutralPV_Phi", */
	/* 						  "recoilPFPuppiMet_NeutralPV_sumEt", */
	/* 						  "recoilPFPuppiMet_NeutralPU_Pt", */
	/* 						  "recoilPFPuppiMet_NeutralPU_Phi", */
	/* 						  "recoilPFPuppiMet_NeutralPU_sumEt", */
	/* 						  "LeadingJet_Phi", */
	/* 						  "LeadingJet_Eta", */
	/* 						  "LeadingJet_M", */
	/* 						  "LeadingJet_Pt", */
	/* 						  "TrailingJet_Phi", */
	/* 						  "TrailingJet_Eta", */
	/* 						  "TrailingJet_M", */
	/* 						  "TrailingJet_Pt", */
	/* 						  "NVertex", */
	/* 						  "NCleanedJets" */
        /* }; */


  std::vector<std::string> variablesForRecoilTraining_  = {};/*   = {"recoilPFPuppiMet_Pt", */
	/* 						     "recoilPFPuppiMet_Phi", */
	/* 						     "recoilPFPuppiMet_sumEt", */
	/* 						     "recoilPFPuppiMet_ChargedPU_Pt", */
	/* 						     "recoilPFPuppiMet_ChargedPU_Phi", */
	/* 						     "recoilPFPuppiMet_ChargedPU_sumEt", */
	/* 						     "recoilPFPuppiMet_ChargedPV_Pt", */
	/* 						     "recoilPFPuppiMet_ChargedPV_Phi", */
	/* 						     "recoilPFPuppiMet_ChargedPV_sumEt", */
	/* 						     "recoilPFPuppiMet_NeutralPV_Pt", */
	/* 						     "recoilPFPuppiMet_NeutralPV_Phi", */
	/* 						     "recoilPFPuppiMet_NeutralPV_sumEt", */
	/* 						     "recoilPFPuppiMet_NeutralPU_Pt", */
	/* 						     "recoilPFPuppiMet_NeutralPU_Phi", */
	/* 						     "recoilPFPuppiMet_NeutralPU_sumEt", */
	/* 						     "LeadingJet_Phi", */
	/* 						     "LeadingJet_Eta", */
	/* 						     "LeadingJet_M", */
	/* 						     "LeadingJet_Pt", */
	/* 						     "TrailingJet_Phi", */
	/* 						     "TrailingJet_Eta", */
	/* 						     "TrailingJet_M", */
	/* 						     "TrailingJet_Pt", */
	/* 						     "NVertex", */
	/* 						     "NCleanedJets" */
        /* }; */
	
  const GBRForest* mvaReaderPhiCorrection_;
  const GBRForest* mvaReaderRecoilCorrection_;

  bool debug_;
  bool produceRecoils_;
}; 
#endif
