#pragma once

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"

namespace JME {

    class IsolatedPhysicsObjectAnalyzer : public PhysicsObjectAnalyzer {
        public:

            // construction/destruction
            explicit IsolatedPhysicsObjectAnalyzer(const edm::ParameterSet& iConfig)
                : PhysicsObjectAnalyzer(iConfig) {
                    // Empty;
                }

            virtual ~IsolatedPhysicsObjectAnalyzer() {
                // Empty
            }

        protected:

            template<typename T>
                void computeIsolation(const T& object) {
                    chargedHadronIso_.push_back(object.chargedHadronIso());
                    neutralHadronIso_.push_back(object.neutralHadronIso());
                    photonIso_.push_back(object.photonIso());
                    puChargedHadronIso_.push_back(object.puChargedHadronIso());
                    relativeIso.push_back(object.chargedHadronIso() + object.neutralHadronIso() + object.photonIso());
                    relativeIso_deltaBeta.push_back((object.chargedHadronIso() + std::max((object.neutralHadronIso() + object.photonIso()) - 0.5 * object.puChargedHadronIso(), 0.0)) / object.pt());
                }

        protected:
            std::vector<float>& chargedHadronIso_ = tree["chargedHadronIso"].write<std::vector<float>>();
            std::vector<float>& neutralHadronIso_ = tree["neutralHadronIso"].write<std::vector<float>>();
            std::vector<float>& photonIso_ = tree["photonIso"].write<std::vector<float>>();
            std::vector<float>& puChargedHadronIso_ = tree["puChargedHadronIso"].write<std::vector<float>>();
            std::vector<float>& relativeIso = tree["relativeIso"].write<std::vector<float>>();
            std::vector<float>& relativeIso_deltaBeta = tree["relativeIso_deltaBeta"].write<std::vector<float>>();
    };
}
