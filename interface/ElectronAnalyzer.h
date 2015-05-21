#pragma once

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "JMEAnalysis/JMEValidator/interface/LeptonAnalyzer.h"


class ElectronAnalyzer: public JME::LeptonAnalyzer {
    public:
        explicit ElectronAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~ElectronAnalyzer();

        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    private:
        virtual float getEffectiveArea(float eta, JME::ConeSize coneSize) override {
            eta = fabs(eta);
            switch (coneSize) {
                case JME::ConeSize::R03:
                    if (eta < 0.8)
                        return 0.1013;
                    else if (eta < 1.3)
                        return 0.0988;
                    else if (eta < 2.0)
                        return 0.0572;
                    else if (eta < 2.2)
                        return 0.0842;
                    else
                        return 0.1530;

                    break;

                case JME::ConeSize::R04:
                    if (eta < 0.8)
                        return 0.1830;
                    else if (eta < 1.3)
                        return 0.1734;
                    else if (eta < 2.0)
                        return 0.1077;
                    else if (eta < 2.2)
                        return 0.1565;
                    else
                        return 0.2680;

                    break;
            }

            return 1;
        }

        edm::EDGetTokenT<pat::ElectronCollection> electrons_;
        edm::EDGetTokenT<reco::VertexCollection> vertices_;
        edm::EDGetTokenT<reco::ConversionCollection> conversions_;
        edm::EDGetTokenT<reco::BeamSpot> beamspot_;
        edm::EDGetTokenT<double> rhoToken_;

        std::vector<std::pair<std::string, edm::EDGetTokenT<edm::ValueMap<bool>>>> idTokens_;

    private:
        std::vector<float>& supercluster_eta = tree["supercluster_eta"].write<std::vector<float>>();
        std::vector<float>& supercluster_phi = tree["supercluster_phi"].write<std::vector<float>>();
        std::vector<float>& dEtaIn = tree["dEtaIn"].write<std::vector<float>>();
        std::vector<float>& dPhiIn = tree["dPhiIn"].write<std::vector<float>>();
        std::vector<float>& hOverE = tree["hOverE"].write<std::vector<float>>();
        std::vector<float>& full5x5_sigmaIetaIeta = tree["full5x5_sigmaIetaIeta"].write<std::vector<float>>();
        std::vector<float>& ooEmooP = tree["ooEmooP"].write<std::vector<float>>();
        std::vector<float>& d0 = tree["d0"].write<std::vector<float>>();
        std::vector<float>& dz = tree["dz"].write<std::vector<float>>();
        std::vector<int>& expectedMissingInnerHits = tree["expected_missing_inner_hits"].write<std::vector<int>>();
        std::vector<bool>& passConversionVeto = tree["pass_conversion_veto"].write<std::vector<bool>>();

        std::vector<std::map<std::string, bool>>& ids_ = tree["ids"].write<std::vector<std::map<std::string, bool>>>();


};

