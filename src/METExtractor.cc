////////////////////////////////////////////////////////////////////////////////
//
// METAnalyzer
// ------------------
//
//                        05/2015    Sébastien Brochet <sebastien.brochet@cern.ch>
////////////////////////////////////////////////////////////////////////////////

#include "JMEAnalysis/JMEValidator/interface/METAnalyzer.h"

#include <vector>
#include <iostream>
#include <regex>
#include <string>

////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
METAnalyzer::METAnalyzer(const edm::ParameterSet& iConfig): JME::PhysicsObjectAnalyzer(iConfig),
    src_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("src")))
{
}


//______________________________________________________________________________
METAnalyzer::~METAnalyzer()
{

}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________

//______________________________________________________________________________
void METAnalyzer::analyze(const edm::Event& iEvent,
        const edm::EventSetup& iSetup)
{

    edm::Handle<std::vector<pat::MET>> metHandle;
    iEvent.getByToken(src_, metHandle);

    const pat::MET& met = metHandle->at(0);

    extractBasicProperties(met);
    sumEt.push_back(met.sumEt());
    significance.push_back(met.significance());

    uncorrectedPt.push_back(met.uncorrectedPt());
    uncorrectedPhi.push_back(met.uncorrectedPhi());
    uncorrectedSumEt.push_back(met.uncorrectedSumEt());

    extractGenProperties(met.genMET());
    if (met.genMET()) {
        gen_sumEt.push_back(met.genMET()->sumEt());
        gen_significance.push_back(met.genMET()->significance());
    } else {
        gen_sumEt.push_back(0);
        gen_significance.push_back(0);
    }

    // FIXME: Calo met will be removed soon
    caloMET_pt.push_back(met.caloMETPt());
    caloMET_phi.push_back(met.caloMETPhi());
    caloMET_sumEt.push_back(met.caloMETSumEt());

    tree.fill();
}

////////////////////////////////////////////////////////////////////////////////
// define METAnalyzer as a plugin
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METAnalyzer);
