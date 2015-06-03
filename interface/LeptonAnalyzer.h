#pragma once

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"

namespace JME {

    enum class ConeSize {
        R03,
        R04
    };

    class LeptonAnalyzer : public PhysicsObjectAnalyzer {
        public:

            // construction/destruction
            explicit LeptonAnalyzer(const edm::ParameterSet& iConfig)
                : PhysicsObjectAnalyzer(iConfig) {
                    // Empty;
                }

            virtual ~LeptonAnalyzer() {
                // Empty
            }

        protected:

            virtual float getEffectiveArea(float eta, ConeSize) {
                return 1;
            }

            template<typename T>
                void computeRelativeIsolationR03(const T& object, float chargedHadronIso, float neutralHadronIso, float photonIso, float puChargedHadronIso, float eta, float rho) {
                    chargedHadronIsoR03_.push_back(chargedHadronIso);
                    neutralHadronIsoR03_.push_back(neutralHadronIso);
                    photonIsoR03_.push_back(photonIso);
                    puChargedHadronIsoR03_.push_back(puChargedHadronIso);
                    relativeIsoR03_.push_back((chargedHadronIso + neutralHadronIso + photonIso) / object.pt());
                    relativeIsoR03_deltaBeta_.push_back((chargedHadronIso + std::max((neutralHadronIso + photonIso) - 0.5f * puChargedHadronIso, 0.0f)) / object.pt());

                    float EA = getEffectiveArea(eta, ConeSize::R03);
                    relativeIsoR03_withEA_.push_back((chargedHadronIso + std::max((neutralHadronIso + photonIso) - rho * EA, 0.0f)) / object.pt());
                }

            template<typename T>
                void computeRelativeIsolationR04(const T& object, float chargedHadronIso, float neutralHadronIso, float photonIso, float puChargedHadronIso, float eta, float rho) {
                    chargedHadronIsoR04_.push_back(chargedHadronIso);
                    neutralHadronIsoR04_.push_back(neutralHadronIso);
                    photonIsoR04_.push_back(photonIso);
                    puChargedHadronIsoR04_.push_back(puChargedHadronIso);
                    relativeIsoR04_.push_back((chargedHadronIso + neutralHadronIso + photonIso) / object.pt());
                    relativeIsoR04_deltaBeta_.push_back((chargedHadronIso + std::max((neutralHadronIso + photonIso) - 0.5f * puChargedHadronIso, 0.0f)) / object.pt());

                    float EA = getEffectiveArea(eta, ConeSize::R04);
                    relativeIsoR04_withEA_.push_back((chargedHadronIso + std::max((neutralHadronIso + photonIso) - rho * EA, 0.0f)) / object.pt());
                }

        protected:
            std::vector<float>& chargedHadronIsoR03_ = tree["chargedHadronIsoR03"].write<std::vector<float>>();
            std::vector<float>& neutralHadronIsoR03_ = tree["neutralHadronIsoR03"].write<std::vector<float>>();
            std::vector<float>& photonIsoR03_ = tree["photonIsoR03"].write<std::vector<float>>();
            std::vector<float>& puChargedHadronIsoR03_ = tree["puChargedHadronIsoR03"].write<std::vector<float>>();
            std::vector<float>& relativeIsoR03_ = tree["relativeIsoR03"].write<std::vector<float>>();
            std::vector<float>& relativeIsoR03_deltaBeta_ = tree["relativeIsoR03_deltaBeta"].write<std::vector<float>>();
            std::vector<float>& relativeIsoR03_withEA_ = tree["relativeIsoR03_withEA"].write<std::vector<float>>();

            std::vector<float>& chargedHadronIsoR04_ = tree["chargedHadronIsoR04"].write<std::vector<float>>();
            std::vector<float>& neutralHadronIsoR04_ = tree["neutralHadronIsoR04"].write<std::vector<float>>();
            std::vector<float>& photonIsoR04_ = tree["photonIsoR04"].write<std::vector<float>>();
            std::vector<float>& puChargedHadronIsoR04_ = tree["puChargedHadronIsoR04"].write<std::vector<float>>();
            std::vector<float>& relativeIsoR04_ = tree["relativeIsoR04"].write<std::vector<float>>();
            std::vector<float>& relativeIsoR04_deltaBeta_ = tree["relativeIsoR04_deltaBeta"].write<std::vector<float>>();
            std::vector<float>& relativeIsoR04_withEA_ = tree["relativeIsoR04_withEA"].write<std::vector<float>>();
    };
}
