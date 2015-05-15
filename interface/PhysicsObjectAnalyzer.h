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

            template<typename T>
                void extractGenProperties(const T& genObject) {
                    PtEtaPhiEVector p(genObject.pt(), genObject.eta(), genObject.phi(), genObject.energy());
                    gen_p4.push_back(p);
                    gen_pt.push_back(p.Pt());
                    gen_eta.push_back(p.Eta());
                    gen_phi.push_back(p.Phi());
                    gen_e.push_back(p.E());

                    gen_m.push_back(genObject.mass());
                    gen_y.push_back(genObject.rapidity());
                    gen_charge.push_back(genObject.charge());
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

            std::vector<PtEtaPhiEVector>& gen_p4 = tree["gen_p4"].write<std::vector<PtEtaPhiEVector>>();
            std::vector<float>& gen_pt = tree["gen_pt"].write<std::vector<float>>();
            std::vector<float>& gen_eta = tree["gen_eta"].write<std::vector<float>>();
            std::vector<float>& gen_phi = tree["gen_phi"].write<std::vector<float>>();
            std::vector<float>& gen_e = tree["gen_e"].write<std::vector<float>>();
            std::vector<float>& gen_m = tree["gen_m"].write<std::vector<float>>();
            std::vector<float>& gen_y = tree["gen_y"].write<std::vector<float>>();
            std::vector<int>& gen_charge = tree["gen_charge"].write<std::vector<int>>();
    };
}
