// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  jing li
//         Created:  Mon, 04 May 2015 17:00:27 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TString.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <Math/VectorUtil.h>
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace std;
using namespace reco;
//
// class declaration
//

class genjetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit genjetAnalyzer(const edm::ParameterSet&);
      ~genjetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
	  edm::InputTag srcGenJets_;
	  edm::Service<TFileService> fs_;
	  std::vector<float> *genjetpt_,*genjeteta_,*genjetphi_,*genjetenergy_;
	  TTree * outTree_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
genjetAnalyzer::genjetAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   srcGenJets_ = iConfig.getParameter<edm::InputTag>        ("genjets");

}


genjetAnalyzer::~genjetAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
genjetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   genjetpt_->clear();
   genjeteta_->clear();
   genjetphi_->clear();
   genjetenergy_->clear();

   edm::Handle<reco::GenJetCollection> genjets;
   iEvent.getByLabel(srcGenJets_,genjets);

   for(reco::GenJetCollection::const_iterator igen = genjets->begin();igen != genjets->end(); ++igen) {
	   genjetpt_->push_back(igen->pt());
	   genjeteta_->push_back(igen->eta());
	   genjetphi_->push_back(igen->phi());
	   genjetenergy_->push_back(igen->energy());
   }
   outTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
genjetAnalyzer::beginJob()
{
	outTree_ = fs_->make<TTree>("t","t");
	genjetpt_        = new std::vector<float>;
	genjeteta_       = new std::vector<float>;
	genjetphi_       = new std::vector<float>;
	genjetenergy_    = new std::vector<float>;
	outTree_->Branch("genjetpt"    ,"vector<float>" ,&genjetpt_);
	outTree_->Branch("genjeteta"   ,"vector<float>" ,&genjeteta_);
	outTree_->Branch("genjetphi"   ,"vector<float>" ,&genjetphi_);
	outTree_->Branch("genjetEnergy","vector<float>" ,&genjetenergy_);
//	cout<<"Begin job finished"<<endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
genjetAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
genjetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
genjetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
genjetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
genjetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
genjetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(genjetAnalyzer);
