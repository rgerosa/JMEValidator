#include "JMEAnalysis/JMEValidator/interface/MuonAnalyzer.h"

#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig): JME::LeptonAnalyzer(iConfig),
    muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho")))
{
    if (iConfig.existsAs<std::vector<edm::InputTag>>("ids")) {
        const std::vector<edm::InputTag>& ids = iConfig.getParameter<std::vector<edm::InputTag>>("ids");

        for (const edm::InputTag& id: ids) {
            idTokens_.push_back(std::make_pair(id.instance(), consumes<edm::ValueMap<bool>>(id)));
        }
    }
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

    std::vector<std::pair<std::string, edm::Handle<edm::ValueMap<bool>>>> idHandles;
    for (auto& token: idTokens_) {
        edm::Handle<edm::ValueMap<bool>> idHandle;
        iEvent.getByToken(token.second, idHandle);
        idHandles.push_back(std::make_pair(token.first, idHandle));
    }

    // Loop over muons
    size_t index = 0;
    for (const pat::Muon& muon: *muonsHandle) {
        extractBasicProperties(muon);
        extractGenProperties(muon.genLepton());

        reco::MuonPFIsolation pfIso = muon.pfIsolationR03();
        computeRelativeIsolationR03(muon, pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt, pfIso.sumPhotonEt, pfIso.sumPUPt, muon.eta(), rho);

        pfIso = muon.pfIsolationR04();
        computeRelativeIsolationR04(muon, pfIso.sumChargedHadronPt, pfIso.sumNeutralHadronEt, pfIso.sumPhotonEt, pfIso.sumPUPt, muon.eta(), rho);

        isLoose_.push_back(muon::isLooseMuon(muon));
        isMedium_.push_back(muon::isMediumMuon(muon));
        isSoft_.push_back(muon::isSoftMuon(muon, primaryVertex));
        isTight_.push_back(muon::isTightMuon(muon, primaryVertex));
        isHighPt_.push_back(muon::isHighPtMuon(muon, primaryVertex));

        // IDs
        const pat::MuonRef muonRef(muonsHandle, index++);
        std::map<std::string, bool> ids;
        for (auto& idHandle: idHandles) {
            ids[idHandle.first] = (*(idHandle.second))[muonRef];
        }
        ids_.push_back(ids);

    }

    tree.fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonAnalyzer);
