#pragma once

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "JMEAnalysis/JMEValidator/interface/IsolatedPhysicsObjectAnalyzer.h"


class MuonAnalyzer: public JME::IsolatedPhysicsObjectAnalyzer {
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
};

