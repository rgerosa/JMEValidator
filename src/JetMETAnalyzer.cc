////////////////////////////////////////////////////////////////////////////////
//
// JetMETAnalyzer
// ------------------
//
//                        01/07/2014 Alexx Perloff   <aperloff@physics.tamu.edu>
//                        04/2015    Sébastien Brochet <sebastien.brochet@cern.ch>
////////////////////////////////////////////////////////////////////////////////

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

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

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

#include "JMEAnalysis/JMEValidator/interface/JetMETAnalyzer.h"

#include <vector>
#include <iostream>
#include <regex>
#include <string>

////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
JetMETAnalyzer::JetMETAnalyzer(const edm::ParameterSet& iConfig)
  : JME::PhysicsObjectAnalyzer(iConfig)
  , JetCorLabel_   (iConfig.getParameter<std::string>("JetCorLabel"))
  , JetCorLevels_  (iConfig.getParameter<std::vector<std::string>>("JetCorLevels"))
  , srcJet_        (consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJet")))
  , srcVtx_        (consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcVtx")))
  , srcMuons_      (consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuons")))
  , doComposition_ (iConfig.getParameter<bool>("doComposition"))
  , doFlavor_      (iConfig.getParameter<bool>("doFlavor"))
  , nJetMax_       (iConfig.getParameter<unsigned int>("nJetMax"))
  , deltaRMax_(0.0)
  , deltaPhiMin_(3.141)
  , deltaRPartonMax_(0.0)
  , jetCorrector_(0)
{
  if (iConfig.exists("deltaRMax")) {
    deltaRMax_=iConfig.getParameter<double>("deltaRMax");
  }
  else
    throw cms::Exception("MissingParameter")<<"Set *either* deltaRMax (matching)"
					    <<" *or* deltaPhiMin (balancing)";

  // Jet CORRECTOR
  // if(!JetCorLevels_.empty()) {
  //   vector<JetCorrectorParameters> vPar;
  //   string jetCorPar = "../data/PHYS14_V2_MC_L1FastJet_"+JetCorLabel_.substr(0,JetCorLabel_.size()-2)+".txt";
  //   cout << "Getting JEC from file " << jetCorPar  << " ... ";
  //   vPar.push_back(JetCorrectorParameters(jetCorPar));
  //   jetCorrector_ = new FactorizedJetCorrector(vPar);
  //   cout << "DONE" << endl;
  // }
  
  std::cout << "|---- JetMETAnalyzer: Initialyzing..." << std::endl;
  std::cout << "|---- JetMETAnalyzer: Applying these jet corrections: ( " << JetCorLabel_;
  for (const std::string& level: JetCorLevels_)
     std::cout << ", " << level;
  std::cout << " )" << std::endl;

  std::cout << "|---- JetMETAnalyzer: RUNNING ON " << moduleLabel_ << " FOR "
       << JetCorLabel_.substr(0,3) << " JETS";
  if      (JetCorLabel_.find("chs") != std::string::npos)   std::cout << " USING CHS" << std::endl;
  else if (JetCorLabel_.find("PUPPI") != std::string::npos) std::cout << " USING PUPPI" << std::endl;
  else                                                      std::cout << std::endl;
}


//______________________________________________________________________________
JetMETAnalyzer::~JetMETAnalyzer()
{

}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________

//______________________________________________________________________________
void JetMETAnalyzer::analyze(const edm::Event& iEvent,
                                  const edm::EventSetup& iSetup)
{

  // // EVENT DATA HANDLES
  edm::Handle<reco::CandidateView>               refs;
  edm::Handle<std::vector<pat::Jet> >            jets;
  edm::Handle<std::vector<reco::Vertex>>         vtx;
  edm::Handle<edm::View<pat::Muon> >             muons;

  iEvent.getByToken(srcVtx_, vtx);

  // REFERENCES & RECOJETS
  iEvent.getByToken(srcJet_, jets);
  
  //loop over the jets and fill the ntuple
  size_t nJet = (nJetMax_ == 0) ? jets->size() : std::min(nJetMax_, (unsigned int) jets->size());
  for (size_t iJet = 0; iJet < nJet; iJet++) {

     pat::Jet const & jet = jets->at(iJet);
     if (jet.pt() < 5)
         continue;

     const reco::GenJet* ref = jet.genJet();

     if (ref) {
         isMatched.push_back(true);
         refdrjt.push_back(reco::deltaR(jet, *ref));
         refpdgid.push_back(ref->pdgId());
         refe.push_back(ref->energy());
         refpt.push_back(ref->pt());
         refeta.push_back(ref->eta());
         refphi.push_back(ref->phi());
         refm.push_back(ref->mass());
         refy.push_back(ref->rapidity());
         refarea.push_back(ref->jetArea());
     }
     else {
         isMatched.push_back(false);
         refdrjt.push_back(0);
         refpdgid.push_back(0.);
         refe.push_back(0.);
         refpt.push_back(0.);
         refeta.push_back(0.);
         refphi.push_back(0.);
         refm.push_back(0.);
         refy.push_back(0.);
         refarea.push_back(0.);
     }

     refrank.push_back(nref);

     // New jet flavor informations
     // See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
     partonFlavor.push_back(jet.partonFlavour());
     hadronFlavor.push_back(jet.hadronFlavour());

     // b-tagging discriminators
     btagDiscri = jet.getPairDiscri();

     // PU Jet Id
     for (const std::string& userFloatName: jet.userFloatNames()) {
        // Look for a string starting with 'pileupJetIdEvaluator'
        if (userFloatName.find("pileupJetIdEvaluator") == 0)
            pujetid_fulldiscriminant.push_back(jet.userFloat(userFloatName));
     }

     for (const std::string& userIntName: jet.userIntNames()) {
         static std::regex cutbasedIdRegex("pileupJetIdEvaluator(.*):cutbasedId");
         static std::regex fullIdRegex("pileupJetIdEvaluator(.*):fullId");

         if (std::regex_match(userIntName, cutbasedIdRegex))
             pujetid_cutbasedid.push_back(jet.userInt(userIntName));

         if (std::regex_match(userIntName, fullIdRegex))
             pujetid_fullid.push_back(jet.userInt(userIntName));
     }

     extractBasicProperties(jet);

     jtarea.push_back(jet.jetArea());
     jtjec.push_back(jet.jecFactor(0));

     computeBetaStar(jet, *vtx);

     nref++;
  }

  tree.fill();
}

void JetMETAnalyzer::computeBetaStar(const pat::Jet& jet, const std::vector<reco::Vertex>& vertices) {

    int nCh_tmp(0), nNeutrals_tmp(0);
    float sumTkPt(0.0);
    float beta_tmp(0.0), betaStar_tmp(0.0), betaStarClassic_tmp(0.0), betaClassic_tmp(0.0);
    float pTMax(0.0), dZ2(-999);
    float sumW(0.0), sumW2(0.0), sumWdR2(0.0);

    const size_t nRings = 9;
    float sum_rings[nRings] = {0};

    for (size_t j = 0; j < jet.numberOfDaughters(); j++) {
        const auto& part = jet.daughterPtr(j);
        if (!(part.isAvailable() && part.isNonnull()) ){
            continue;
        }

        if (fabs(part->charge()) > 0)
            nCh_tmp++;
        else 
            nNeutrals_tmp++;

        float dR = reco::deltaR(*part, jet);
        float weight = part->pt();
        float weight2 = weight * weight;
        sumWdR2      += weight2 * dR * dR;
        sumW         += weight;
        sumW2        += weight2;

        size_t index = (int) (dR * 10);
        if (index > nRings - 1)
            index = nRings - 1;
        sum_rings[index] += weight;

        reco::CandidatePtr pfJetConstituent = jet.sourceCandidatePtr(j);
        const reco::Candidate* icand = pfJetConstituent.get();
        const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>( icand );
        if (lPack) {

            if (fabs(lPack->charge()) > 0) {

                if (lPack->pt() > pTMax) {
                    pTMax = lPack->pt();
                    dZ2 = lPack->dz();
                }

                float tkpt = lPack->pt();
                sumTkPt += tkpt;
                bool inVtx0 = (lPack->fromPV() == 3);
                bool inAnyOther = (lPack->fromPV() == 0);
                double dZ0 = lPack->dz();
                double dZ_tmp = dZ0;

                for (const auto& iv: vertices) {
                    if (iv.isFake())
                        continue;
                    if (fabs(lPack->dz(iv.position())) < fabs(dZ_tmp)) {
                        dZ_tmp = lPack->dz(iv.position());
                    }
                }

                if (inVtx0) {
                    betaClassic_tmp += tkpt;
                }
                else if (inAnyOther) {
                    betaStarClassic_tmp += tkpt;
                }
                if (fabs(dZ0) < 0.2) {
                    beta_tmp += tkpt;
                }
                else if (fabs(dZ_tmp) < 0.2) {
                    betaStar_tmp += tkpt;
                }
            }
        }
    }

    if (sumW > 0) {
        DRweighted.push_back(sumWdR2 / sumW2);
        fRing0.push_back(sum_rings[0] / sumW);
        fRing1.push_back(sum_rings[1] / sumW);
        fRing2.push_back(sum_rings[2] / sumW);
        fRing3.push_back(sum_rings[3] / sumW);
        fRing4.push_back(sum_rings[4] / sumW);
        fRing5.push_back(sum_rings[5] / sumW);
        fRing6.push_back(sum_rings[6] / sumW);
        fRing7.push_back(sum_rings[7] / sumW);
        fRing8.push_back(sum_rings[8] / sumW);
        ptD.push_back(sqrt(sumW2) / sumW);
    }
    else{
        DRweighted.push_back(-999);
        fRing0.push_back(-999);
        fRing1.push_back(-999);
        fRing2.push_back(-999);
        fRing3.push_back(-999);
        fRing4.push_back(-999);
        fRing5.push_back(-999);
        fRing6.push_back(-999);
        fRing7.push_back(-999);
        fRing8.push_back(-999);
        ptD.push_back(-999);
    }
    if (sumTkPt > 0) {
        beta.push_back(beta_tmp/sumTkPt);
        betaStar.push_back(betaStar_tmp/sumTkPt);
        betaClassic.push_back(betaClassic_tmp/sumTkPt);
        betaStarClassic.push_back(betaStarClassic_tmp/sumTkPt);
    }
    else{
        beta.push_back(-999);
        betaStar.push_back(-999);
        betaClassic.push_back(-999);
        betaStarClassic.push_back(-999);
    }

    dZ.push_back(dZ2);
    nCh.push_back(nCh_tmp);
    nNeutrals.push_back(nNeutrals_tmp);
}


////////////////////////////////////////////////////////////////////////////////
// define JetMETAnalyzer as a plugin
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetMETAnalyzer);
