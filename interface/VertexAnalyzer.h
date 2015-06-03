#pragma once

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

#include <vector>
#include <string>

class VertexAnalyzer: public JME::Analyzer {
    public:
        explicit VertexAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~VertexAnalyzer();

        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    private:
        edm::EDGetTokenT<reco::VertexCollection> srcToken_;

    private:
        // Tree members
        std::vector<float>& normalizedChi2_ = tree["normalizedChi2"].write<std::vector<float>>();
        std::vector<float>& ndof_ = tree["ndof"].write<std::vector<float>>();
        std::vector<bool>& isFake_ = tree["isFake"].write<std::vector<bool>>();
        std::vector<bool>& isValid_ = tree["isValid"].write<std::vector<bool>>();
        std::vector<reco::Vertex::Point>& position_ = tree["position"].write<std::vector<reco::Vertex::Point>>();
        std::vector<reco::Vertex::CovarianceMatrix>& covariance_ = tree["covariance"].write<std::vector<reco::Vertex::CovarianceMatrix>>();
};
