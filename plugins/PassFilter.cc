#ifndef PassFilter_h
#define PassFilter_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"



class PassFilter : public edm::EDFilter {
 
public:
 
  //! ctor
  explicit PassFilter (const edm::ParameterSet&);
 
  //! dtor 
  ~PassFilter();
 
private:
 
  //! the actual filter method 
  bool filter(edm::Event&, const edm::EventSetup&);
 
private:
 
  TH1F* _totalEvents; 
};

#endif

//! ctor
PassFilter::PassFilter(const edm::ParameterSet& iConfig) {
  edm::Service<TFileService> fs;
  _totalEvents = fs -> make<TH1F>("totalEvents", "totalEvents", 1,  0., 1.);
}

// ----------------------------------------------------------------

//! dtor
PassFilter::~PassFilter()
{}

// ----------------------------------------------------------------


//! loop over the reco particles and count leptons
bool PassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  _totalEvents -> Fill(0.5);
  return true;
}

DEFINE_FWK_MODULE(PassFilter);
