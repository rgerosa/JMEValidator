////////////////////////////////////////////////////////////////////////////////
//
// puppiReader
// -----------
//
//                        01/07/2014 Alexx Perloff   <aperloff@physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////

#include "JMEAnalysis/JMEValidator/interface/puppi_Ntuple.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/View.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////

class puppiReader : public edm::EDAnalyzer
{
public:
  // construction/destruction
  explicit puppiReader(const edm::ParameterSet& iConfig);
  virtual ~puppiReader();

private:
  // member functions
  void beginJob();
  void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob(){;}

private:
  // member data
  std::string   moduleLabel_;
  std::string   JetCorLabel_;
  std::vector<std::string> JetCorLevels_;
  
  // tree
  TTree*        tree_;
  puppiNtuple*  Ntuple_;
};


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
puppiReader::puppiReader(const edm::ParameterSet& iConfig)
{

}


//______________________________________________________________________________
puppiReader::~puppiReader()
{

}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void puppiReader::beginJob()
{
  edm::Service<TFileService> fs;
  if (!fs) throw edm::Exception(edm::errors::Configuration,
				"TFileService missing from configuration!");

  tree_=fs->make<TTree>("puppiTree","puppiTree");
  Ntuple_ = new puppiNtuple(tree_,true);
}


//______________________________________________________________________________
void puppiReader::analyze(const edm::Event& iEvent,
                                  const edm::EventSetup& iSetup)
{

  edm::Handle<double> nalgos;
  edm::InputTag src_nalgos = edm::InputTag("puppi","PuppiNAlgos","JRA");
  iEvent.getByLabel(src_nalgos,nalgos);

  edm::Handle<std::vector<double>> alphas;
  edm::InputTag src_alphas = edm::InputTag("puppi","PuppiRawAlphas","JRA");
  iEvent.getByLabel(src_alphas,alphas);

  // edm::Handle< std::vector< pat::PackedCandidate > > packedCands;
  edm::InputTag src_packedCands = edm::InputTag("packedPFCandidates","","PAT");
  // iEvent.getByLabel(src_packedCands,packedCands);
  edm::Handle<reco::CandidateView> hPFProduct;
  iEvent.getByLabel(src_packedCands,hPFProduct);
  const reco::CandidateView *pfCol = hPFProduct.product();

  //EVENT INFORMATION
  Ntuple_->run = iEvent.id().run();
  Ntuple_->lumi = iEvent.id().luminosityBlock();
  Ntuple_->evt = iEvent.id().event();

  // std::cout << "nalgos = " << *nalgos << std::endl;
  Ntuple_->nalgos = float(*nalgos);

  Ntuple_->px->clear();
  Ntuple_->py->clear();
  Ntuple_->pz->clear();
  Ntuple_->e->clear();
  Ntuple_->alphas->clear();
  Ntuple_->id->clear();
  Ntuple_->charge->clear();
  Ntuple_->fromPV->clear();

  // std::cout << "reading candidates..." << std::endl;
  for(reco::CandidateView::const_iterator itPF = pfCol->begin(); itPF!=pfCol->end(); itPF++) {
    //const reco::PFCandidate *pPF = dynamic_cast<const reco::PFCandidate*>(&(*itPF));
    Ntuple_->px->push_back( itPF->px() );
    Ntuple_->py->push_back( itPF->py() );
    Ntuple_->pz->push_back( itPF->pz() );
    Ntuple_->e->push_back( itPF->energy() );
    Ntuple_->id->push_back( float(itPF->pdgId()) );
    Ntuple_->charge->push_back( float(itPF->charge()) );
    
    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&(*itPF));
    Ntuple_->fromPV->push_back( float(lPack->fromPV()) );
  }


  // std::cout << "reading alphas..." << std::endl;
  for(unsigned int i = 0; i < alphas->size(); i++){
    Ntuple_->alphas->push_back( (*alphas)[i] );
  }

  tree_->Fill();
  
  return;
}


////////////////////////////////////////////////////////////////////////////////
// define puppiReader as a plugin
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(puppiReader);
