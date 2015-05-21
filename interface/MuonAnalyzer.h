#pragma once

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "JMEAnalysis/JMEValidator/interface/LeptonAnalyzer.h"


class MuonAnalyzer: public JME::LeptonAnalyzer {
    public:
        explicit MuonAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~MuonAnalyzer();

        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    private:
        virtual float getEffectiveArea(float eta, JME::ConeSize coneSize) override {
            eta = fabs(eta);
            switch (coneSize) {
                case JME::ConeSize::R03:
                    if (eta < 0.8)
                        return 0.0913;
                    else if (eta < 1.3)
                        return 0.0765;
                    else if (eta < 2.0)
                        return 0.0546;
                    else if (eta < 2.2)
                        return 0.0728;
                    else
                        return 0.1177;

                    break;

                case JME::ConeSize::R04:
                    if (eta < 0.8)
                        return 0.1564;
                    else if (eta < 1.3)
                        return 0.1325;
                    else if (eta < 2.0)
                        return 0.0913;
                    else if (eta < 2.2)
                        return 0.1212;
                    else
                        return 0.2085;

                    break;
            }

            return 1;
        }

        edm::EDGetTokenT<pat::MuonCollection> muons_;
        edm::EDGetTokenT<reco::VertexCollection> vertices_;
        edm::EDGetTokenT<double> rhoToken_;
        std::vector<std::pair<std::string, edm::EDGetTokenT<edm::ValueMap<bool>>>> idTokens_;

    private:
        std::vector<bool>& isLoose_ = tree["isLoose"].write<std::vector<bool>>();
        std::vector<bool>& isMedium_ = tree["isMedium"].write<std::vector<bool>>();
        std::vector<bool>& isSoft_ = tree["isSoft"].write<std::vector<bool>>();
        std::vector<bool>& isTight_ = tree["isTight"].write<std::vector<bool>>();
        std::vector<bool>& isHighPt_ = tree["isHighPt"].write<std::vector<bool>>();

        std::vector<std::map<std::string, bool>>& ids_ = tree["ids"].write<std::vector<std::map<std::string, bool>>>();
};

