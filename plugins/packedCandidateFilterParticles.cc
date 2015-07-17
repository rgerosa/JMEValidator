#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
 
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include <iostream>
#include <memory>

using namespace std;
using namespace edm;
using namespace pat;


////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
class packedCandidateFilterParticles : public edm::stream::EDProducer<> {

public:
  // construction/destruction
  packedCandidateFilterParticles(const edm::ParameterSet& iConfig);
  ~packedCandidateFilterParticles() {;}
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();

private:
  // input particle to be considered : packed candidate
  edm::InputTag src_;
  // leptons to be cleaned 
  edm::InputTag srcMuons_;
  edm::InputTag srcElectrons_;
  edm::InputTag srcTaus_;
  // vertex collection
  edm::InputTag vertex_;


  // tokens
  edm::EDGetTokenT<pat::PackedCandidateCollection>    srcToken_ ;
  edm::EDGetTokenT<pat::MuonCollection>               srcMuonsToken_ ;
  edm::EDGetTokenT<pat::ElectronCollection>           srcElectronsToken_ ;
  edm::EDGetTokenT<pat::TauCollection>                srcTausToken_ ;

};




////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
packedCandidateFilterParticles::packedCandidateFilterParticles(const edm::ParameterSet& iConfig){

  if(iConfig.existsAs<edm::InputTag >("src")){
    src_ = iConfig.getParameter<edm::InputTag>("src");
    if(src_ == edm::InputTag(""))
      throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] no packed candidates --> at least one \n"; 
  }
  else throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] no packed candidates is given \n";  

  if(iConfig.existsAs<edm::InputTag >("srcMuons")){
    srcMuons_ = iConfig.getParameter<edm::InputTag>("srcMuons");
  }
  else throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] no muons \n";  

  if(iConfig.existsAs<edm::InputTag >("srcElectrons")){
    srcElectrons_ = iConfig.getParameter<edm::InputTag>("srcElectrons");
  }
  else throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] no electrons \n";  

  if(iConfig.existsAs<edm::InputTag >("srcTaus")){
    srcTaus_ = iConfig.getParameter<edm::InputTag>("srcTaus");
  }
  else throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] no Taus \n";  

  // tokens
  if(!(src_ == edm::InputTag(""))) 
    srcToken_ = consumes<pat::PackedCandidateCollection>(src_);

  if(!(srcMuons_ == edm::InputTag(""))) 
    srcMuonsToken_ = consumes<pat::MuonCollection>(srcMuons_);

  if(!(srcElectrons_ == edm::InputTag(""))) 
    srcElectronsToken_ = consumes<pat::ElectronCollection>(srcElectrons_);

  if(!(srcTaus_ == edm::InputTag(""))) 
    srcTausToken_ = consumes<pat::TauCollection>(srcTaus_);

   produces<pat::PackedCandidateCollection>();

}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void packedCandidateFilterParticles::produce(edm::Event& iEvent,const edm::EventSetup& iSetup){

  pat::PackedCandidateCollection packedCandidateColl ;
  edm::Handle<pat::PackedCandidateCollection> packedCandidateHandle;
  if(!(src_ == edm::InputTag(""))){
    iEvent.getByToken(srcToken_,packedCandidateHandle);
    if(packedCandidateHandle.failedToGet()){
      throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] failed to get packed candidate collection \n";
    }
    packedCandidateColl = *(packedCandidateHandle.product());
  }
 
  edm::Handle<pat::MuonCollection> MuonCollectionHandle;
  pat::MuonCollection muonColl;
  if(!(srcMuons_ == edm::InputTag(""))){  
    iEvent.getByToken(srcMuonsToken_,MuonCollectionHandle);
    if(MuonCollectionHandle.failedToGet()){
      throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] failed to get muon collection \n"; 
    }
    muonColl = *(MuonCollectionHandle.product());
  }    

  edm::Handle<pat::ElectronCollection> ElectronCollectionHandle;
  pat::ElectronCollection electronColl;
  if(!(srcElectrons_ == edm::InputTag(""))){  
    iEvent.getByToken(srcElectronsToken_,ElectronCollectionHandle);
    if(ElectronCollectionHandle.failedToGet()){
      throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] failed to get electron collection \n"; 
    }
    electronColl = *(ElectronCollectionHandle.product());
  }    

  edm::Handle<pat::TauCollection> TauCollectionHandle;
  pat::TauCollection tauColl;
  if(!(srcTaus_ == edm::InputTag(""))){  
    iEvent.getByToken(srcTausToken_,TauCollectionHandle);
    if(TauCollectionHandle.failedToGet()){
      throw cms::Exception("Configuration")<<"[packedCandidateFilterParticles] failed to get tau collection \n"; 
    }
    tauColl = *(TauCollectionHandle.product());
  }    

  auto_ptr<pat::PackedCandidateCollection> filteredPackedCandidates(new vector<pat::PackedCandidate>);

  for( auto particle : packedCandidateColl){

    if( fabs(particle.pdgId()) !=13 and fabs(particle.pdgId()) !=15 and fabs(particle.pdgId()) !=11){
      filteredPackedCandidates->push_back(particle);
      continue; 
    }

    // check muons 
    if(fabs(particle.pdgId()) == 13) {
      bool goodMuon = true ;
      for( auto muon : muonColl){
	if( reco::deltaR(particle.p4(),muon.p4()) < 0.01){
	  goodMuon = false;
	  break;
	}
      }

      if(goodMuon and fabs(particle.pdgId()) == 13){
	filteredPackedCandidates->push_back(particle);	
	continue;
      }
      else continue;
    }

    // check electrons 
    if(fabs(particle.pdgId()) == 11) {
      bool goodElectron = true ;
      for( auto ele : electronColl){
	if( reco::deltaR(particle.p4(),ele.p4()) < 0.01){
	  goodElectron = false;
	  break;
	}
      }

      if(goodElectron and fabs(particle.pdgId()) == 11){
	filteredPackedCandidates->push_back(particle);	
	continue;
      }
      else continue;
    }

    // check taus
    if(fabs(particle.pdgId()) == 15) {
      bool goodTau = true ;
      for( auto tau : tauColl){
	if( reco::deltaR(particle.p4(),tau.p4()) < 0.01){
	  goodTau = false;
	  break;
	}
      }

      if(goodTau and fabs(particle.pdgId()) == 15){
	filteredPackedCandidates->push_back(particle);	
	continue;
      }
      else continue;
    }

  }

  iEvent.put(filteredPackedCandidates);

}


//______________________________________________________________________________
void packedCandidateFilterParticles::endJob(){}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(packedCandidateFilterParticles);
