////////////////////////////////////////////////////////////////////////////////
//
// LeptonsAndMETAnalyzer
// -----------
//
//                        01/07/2014 Alexx Perloff   <aperloff@physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////

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

#include "JMEAnalysis/JMEValidator/interface/LeptonsAndMETAnalyzer.h"

#include <vector>

////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
LeptonsAndMETAnalyzer::LeptonsAndMETAnalyzer(const edm::ParameterSet& iConfig)
    : JME::Analyzer(iConfig)
    , srcIsoMuons_   (consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("srcIsoMuons")))
    , srcMET_        (consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET")))
    , srcPUPPET_     (consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("srcPUPPET")))
{

}


//______________________________________________________________________________
LeptonsAndMETAnalyzer::~LeptonsAndMETAnalyzer()
{

}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void LeptonsAndMETAnalyzer::beginJob()
{
  JME::Analyzer::beginJob();
  eventCounter_ = 0;
}


//______________________________________________________________________________
void LeptonsAndMETAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<reco::Candidate> >           muonsForZ;
  edm::Handle<std::vector<pat::MET> >            mets;
  edm::Handle<std::vector<reco::PFMET> >         puppet;

  //EVENT INFORMATION
  run = iEvent.id().run();
  lumiBlock = iEvent.id().luminosityBlock();
  event = iEvent.id().event();

  iEvent.getByToken(srcMET_,mets);
  assert(mets.isValid());
  const pat::MET  &inPFMET = mets.product()->front();
  pfMET      = inPFMET.pt();
  pfMETphi   = inPFMET.phi();

  iEvent.getByToken(srcPUPPET_,puppet);
  assert(mets.isValid());
  const reco::PFMET &inPuppET = puppet.product()->front();
  puppET      = inPuppET.pt();
  puppETphi   = inPuppET.phi();

  iEvent.getByToken(srcIsoMuons_, muonsForZ);
  // std::cout << muonsForZ->size() << std::endl;
  int nZ = 0;
  TLorentzVector theZCand;
  for(size_t i = 0, n = muonsForZ->size(); i < n; ++i) {
    edm::Ptr<reco::Candidate> muPtr1 = muonsForZ->ptrAt(i);
    for(size_t j = 1, n = muonsForZ->size(); j < n; ++j) {
      edm::Ptr<reco::Candidate> muPtr2 = muonsForZ->ptrAt(j);
      if (muPtr1->charge()*muPtr2->charge() < 1){
        // std::cout << " oppo " << std::endl;
        TLorentzVector Zcand = TLorentzVector(muPtr1->px(),muPtr1->py(),muPtr1->pz(),muPtr1->energy()) + TLorentzVector(muPtr2->px(),muPtr2->py(),muPtr2->pz(),muPtr2->energy());
        float Zmass = Zcand.M();
        // std::cout << "Zmass = " << Zmass << std::endl;
        float Zdiff = fabs(Zmass - 91.2);
        if (Zdiff < 15){
          nZ += 1;
          theZCand = Zcand;
        }
      }
    }
  }
  

  if (nZ > 0){
    nZcands = nZ;
    ZpT = theZCand.Pt();
    Zphi = theZCand.Phi();
    Zmass = theZCand.M();
    pfMET_uPara = 0.;
    pfMET_uPerp = 0.;
    puppET_uPara = 0.;
    puppET_uPerp = 0.;
    recoilComputation( pfMET, pfMETphi, ZpT, Zphi, pfMET_uPara, pfMET_uPerp);
    recoilComputation( puppET, puppETphi, ZpT, Zphi, puppET_uPara, puppET_uPerp);
  }
  else{
    nZcands = 0;
    ZpT = -1;
    Zphi = -1;
    Zmass = -1;
    pfMET_uPara = -1;
    pfMET_uPerp = -1;
    puppET_uPara = -1;
    puppET_uPerp = -1;
  }

  tree.fill();  
  eventCounter_++;
}

void LeptonsAndMETAnalyzer::recoilComputation( float &met, float &metPhi, float &Zpt, float &Zphi, float &upara, float &uperp){
  
  double p1UX  = met*cos(metPhi) + Zpt*cos(Zphi);
  double p1UY  = met*sin(metPhi) + Zpt*sin(Zphi);
  double p1U   = sqrt(p1UX*p1UX+p1UY*p1UY);
  double p1Cos = - (p1UX*cos(Zphi) + p1UY*sin(Zphi))/p1U;
  double p1Sin =   (p1UX*sin(Zphi) - p1UY*cos(Zphi))/p1U;
  upara   = p1U*p1Cos;//*(pU1*(iGenPt > 10) + (iGenPt > 10)*((1.-iGenPt/10.)*(pU1-1.)+1.));
  uperp   = p1U*p1Sin;

}


////////////////////////////////////////////////////////////////////////////////
// define LeptonsAndMETAnalyzer as a plugin
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LeptonsAndMETAnalyzer);
