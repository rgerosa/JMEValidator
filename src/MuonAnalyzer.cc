#include "JMEAnalysis/JMEValidator/interface/MuonAnalyzer.h"

MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig): JME::PhysicsObjectAnalyzer(iConfig),
    muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
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

    // Loop over muons
    for (const pat::Muon& muon: *muonsHandle) {
        extractBasicProperties(muon);
        extractGenProperties(muon.genLepton());

        isLoose_.push_back(muon.isLooseMuon());
        isSoft_.push_back(muon.isSoftMuon(primaryVertex));
        isTight_.push_back(muon.isTightMuon(primaryVertex));
        isHighPt_.push_back(muon.isHighPtMuon(primaryVertex));

        // Isolation
        chargedHadronIso_.push_back(muon.chargedHadronIso());
        neutralHadronIso_.push_back(muon.neutralHadronIso());
        photonIso_.push_back(muon.photonIso());
        puChargedHadronIso_.push_back(muon.puChargedHadronIso());
        relativeIso.push_back(muon.chargedHadronIso() + muon.neutralHadronIso() + muon.photonIso());
        relativeIso_deltaBeta.push_back((muon.chargedHadronIso() + std::max((muon.neutralHadronIso() + muon.photonIso()) - 0.5 * muon.puChargedHadronIso(), 0.0)) / muon.pt());
    }

    tree.fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonAnalyzer);
