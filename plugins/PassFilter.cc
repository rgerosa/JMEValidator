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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
  TH1F* _weightEvents; 
  bool isMC_;
  edm::EDGetTokenT<GenEventInfoProduct> srcGenEventInfoToken_;
  edm::InputTag srcGenEventInfo_;
};

#endif

//! ctor
PassFilter::PassFilter(const edm::ParameterSet& iConfig) {
  edm::Service<TFileService> fs;
  _totalEvents = fs -> make<TH1F>("totalEvents",   "totalEvents", 1,  0., 1.);
  _weightEvents = fs -> make<TH1F>("weightEvents", "weightEvents", 2,  -2., 2.);

  if (iConfig.existsAs<edm::InputTag>("srcGenEventInfo"))
    srcGenEventInfo_ = iConfig.getParameter<edm::InputTag>("srcGenEventInfo");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input no gen event info \n";
  if (iConfig.existsAs<bool>("isMC"))
    isMC_ = iConfig.getParameter<bool>("isMC");
  else throw cms::Exception("Configuration")<<"[PUPPETAnalyzer] input isMC \n";
  if(!(srcGenEventInfo_ == edm::InputTag("")) and isMC_)
    srcGenEventInfoToken_ = consumes<GenEventInfoProduct>(srcGenEventInfo_);
}

// ----------------------------------------------------------------

//! dtor
PassFilter::~PassFilter()
{};

// ----------------------------------------------------------------


//! loop over the reco particles and count leptons
bool PassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  _totalEvents -> Fill(0.5);
  if(isMC_){
    edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
    iEvent.getByToken(srcGenEventInfoToken_,GenEventInfoHandle);
    _weightEvents -> Fill(GenEventInfoHandle->weight()/fabs(GenEventInfoHandle->weight()));
  }
  return true;
}

DEFINE_FWK_MODULE(PassFilter);
