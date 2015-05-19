#include "JMEAnalysis/JMEValidator/interface/MuonAnalyzer.h"

#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig): JME::LeptonAnalyzer(iConfig),
    muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho")))
{
    // Empty
}

MuonAnalyzer::~MuonAnalyzer() {
    // Empty
}

void MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<pat::MuonCollection> muonsHandle;
    iEvent.getByToken(muons_, muonsHandle);

    edm::Handle<reco::VertexCollection> verticesHandle;
    iEvent.getByToken(vertices_, verticesHandle);
    const auto& primaryVertex = verticesHandle->at(0);

    edm::Handle<double> rhoHandle;
    iEvent.getByToken(rhoToken_, rhoHandle);
    double rho = *rhoHandle;

    // Loop over muons
    for (const pat::Muon& muon: *muonsHandle) {
        extractBasicProperties(muon);
        extractGenProperties(muon.genLepton());

        reco::MuonPFIsolation pfIso = muon.pfIsolationR03();
        computeRelativeIsolationR03(muon, pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt, pfIso.sumPhotonEt, pfIso.sumPUPt, muon.eta(), rho);

        pfIso = muon.pfIsolationR04();
        computeRelativeIsolationR04(muon, pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt, pfIso.sumPhotonEt, pfIso.sumPUPt, muon.eta(), rho);

        isLoose_.push_back(muon.isLooseMuon());
        isSoft_.push_back(muon.isSoftMuon(primaryVertex));
        isTight_.push_back(muon.isTightMuon(primaryVertex));
        isHighPt_.push_back(muon.isHighPtMuon(primaryVertex));

    }

    tree.fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonAnalyzer);
