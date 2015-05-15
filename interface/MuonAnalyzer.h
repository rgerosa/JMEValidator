#pragma once

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"


class MuonAnalyzer: public JME::PhysicsObjectAnalyzer {
    public:
        explicit MuonAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~MuonAnalyzer();

        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    private:
        edm::EDGetTokenT<pat::MuonCollection> muons_;
        edm::EDGetTokenT<reco::VertexCollection> vertices_;

    private:
        std::vector<bool>& isLoose_ = tree["isLoose"].write<std::vector<bool>>();
        std::vector<bool>& isSoft_ = tree["isSoft"].write<std::vector<bool>>();
        std::vector<bool>& isTight_ = tree["isTight"].write<std::vector<bool>>();
        std::vector<bool>& isHighPt_ = tree["isHighPt"].write<std::vector<bool>>();

        // Isolation
        std::vector<float>& chargedHadronIso_ = tree["chargedHadronIso"].write<std::vector<float>>();
        std::vector<float>& neutralHadronIso_ = tree["neutralHadronIso"].write<std::vector<float>>();
        std::vector<float>& photonIso_ = tree["photonIso"].write<std::vector<float>>();
        std::vector<float>& puChargedHadronIso_ = tree["puChargedHadronIso"].write<std::vector<float>>();
        std::vector<float>& relativeIso = tree["relativeIso"].write<std::vector<float>>();
        std::vector<float>& relativeIso_deltaBeta = tree["relativeIso_deltaBeta"].write<std::vector<float>>();
};

