#pragma once

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

#include <Math/Vector4D.h>


namespace JME {

    class PhysicsObjectAnalyzer : public Analyzer {
        public:

            typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> PtEtaPhiEVector;

            // construction/destruction
            explicit PhysicsObjectAnalyzer(const edm::ParameterSet& iConfig)
                : Analyzer(iConfig) {
                    // Empty;
                }

            virtual ~PhysicsObjectAnalyzer() {
                // Empty
            }

        protected:

            template<typename T>
                void extractBasicProperties(const T& object) {
                    PtEtaPhiEVector p(object.pt(), object.eta(), object.phi(), object.energy());
                    p4.push_back(p);
                    pt.push_back(p.Pt());
                    eta.push_back(p.Eta());
                    phi.push_back(p.Phi());
                    e.push_back(p.E());

                    m.push_back(object.mass());
                    y.push_back(object.rapidity());
                    charge.push_back(object.charge());
                }


        protected:
            std::vector<PtEtaPhiEVector>& p4 = tree["p4"].write<std::vector<PtEtaPhiEVector>>();
            std::vector<float>& pt = tree["pt"].write<std::vector<float>>();
            std::vector<float>& eta = tree["eta"].write<std::vector<float>>();
            std::vector<float>& phi = tree["phi"].write<std::vector<float>>();
            std::vector<float>& e = tree["e"].write<std::vector<float>>();
            std::vector<float>& m = tree["m"].write<std::vector<float>>();
            std::vector<float>& y = tree["y"].write<std::vector<float>>();
            std::vector<int>& charge = tree["charge"].write<std::vector<int>>();
    };
}
