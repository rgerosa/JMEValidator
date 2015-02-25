////////////////////////////////////////////////////////////////////////////////
//
// convertRecoCandToPackedCand
// ---------------------------
//
//                          02/24/2015 Alexx Perloff <aperloff@physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <iostream>
#include <memory>

using namespace std;
using namespace edm;
using namespace reco;

typedef edm::View<reco::PFCandidate> PFCandidateView;
enum { dimension = 3 };
typedef math::Error<dimension>::type Error;
typedef math::Error<dimension>::type CovarianceMatrix;

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
class convertRecoCandToPackedCand : public edm::EDProducer
{
public:
  // construction/destruction
  convertRecoCandToPackedCand(const edm::ParameterSet& iConfig);
  ~convertRecoCandToPackedCand() {;}
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();

private:
  // member data
  edm::InputTag src_;
  edm::InputTag srcFromPVLoose_;
  edm::InputTag srcVtx_;

  edm::Handle<reco::CandidateView> PFCands_;
  edm::Handle<reco::PFCandidateFwdPtrVector> PFCandsFromPVLoose_;
  edm::Handle<reco::VertexCollection> vtxs_;

  std::string  moduleName_;

};




////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
convertRecoCandToPackedCand::convertRecoCandToPackedCand(const edm::ParameterSet& iConfig)
  : src_(iConfig.getParameter<InputTag>("src"))
  , srcFromPVLoose_(iConfig.getParameter<InputTag>("srcFromPVLoose"))
  , srcVtx_(iConfig.getParameter<InputTag>("srcVtx"))
  , moduleName_(iConfig.getParameter<string>("@module_label"))
{
  produces<pat::PackedCandidateCollection>();
}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void convertRecoCandToPackedCand::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  auto_ptr<pat::PackedCandidateCollection> packedCands(new vector<pat::PackedCandidate>);
  
  iEvent.getByLabel(src_,PFCands_);
  iEvent.getByLabel(srcFromPVLoose_,PFCandsFromPVLoose_);
  iEvent.getByLabel(srcVtx_,vtxs_);
  
  size_t nCands = (size_t)PFCands_->size();

  std::vector<pat::PackedCandidate::PVAssoc> fromPV(PFCands_->size(), pat::PackedCandidate::NoPV);
  for (const reco::PFCandidateFwdPtr &ptr : *PFCandsFromPVLoose_) {
     if (ptr.ptr().id() == PFCands_.id()) {
        fromPV[ptr.ptr().key()] = pat::PackedCandidate::PVLoose;
     } else if (ptr.backPtr().id() == PFCands_.id()) {
        fromPV[ptr.backPtr().key()] = pat::PackedCandidate::PVLoose;
     } else {
        throw cms::Exception("Configuration", "The elements from 'inputCollectionFromPVLoose' don't point to 'inputCollection'\n");
     }
  }
  
  reco::VertexRef PV(vtxs_.id());
  if (!vtxs_->empty()) {
     PV = reco::VertexRef(vtxs_, 0);
  }  
    
  for(unsigned int iCand = 0; iCand<nCands; iCand++) {
     reco::CandidateBaseRef rCand = PFCands_->refAt(iCand);
     //const Error cm = PFCands_->at(iCand).vertexCovariance();
     //const reco::VertexRef vertex(PFCands_->at(iCand).vertex(),cm,
     //                             PFCands_->at(iCand).vertexChi2(),PFCands_->at(iCand).vertexNdof(),
     //                             size_t(100000));
     //pat::PackedCandidate intPackedCand(*rCand,vertex);
     pat::PackedCandidate intPackedCand(*rCand,PV);
     packedCands->push_back(intPackedCand);
     packedCands->back().setFromPV(fromPV[iCand]);
  }

  iEvent.put(packedCands);
}


//______________________________________________________________________________
void convertRecoCandToPackedCand::endJob()
{
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(convertRecoCandToPackedCand);
