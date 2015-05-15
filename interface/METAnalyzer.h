#pragma once

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"

#include "DataFormats/PatCandidates/interface/MET.h"

class METAnalyzer : public JME::PhysicsObjectAnalyzer
{
    public:
        // construction/destruction
        explicit METAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~METAnalyzer();

    private:

        // member functions
        void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);

    private:
        edm::EDGetTokenT<std::vector<pat::MET>> src_;

        // Tree branches
        std::vector<float>& sumEt = tree["sumEt"].write<std::vector<float>>();
        std::vector<float>& significance = tree["significance"].write<std::vector<float>>();

        std::vector<float>& uncorrectedPt = tree["uncorrectedPt"].write<std::vector<float>>();
        std::vector<float>& uncorrectedPhi = tree["uncorrectedPhi"].write<std::vector<float>>();
        std::vector<float>& uncorrectedSumEt = tree["uncorrectedSumEt"].write<std::vector<float>>();

        std::vector<float>& gen_sumEt = tree["gen_sumEt"].write<std::vector<float>>();
        std::vector<float>& gen_significance = tree["gen_significance"].write<std::vector<float>>();

        // Calo met
        std::vector<float>& caloMET_pt = tree["caloMET_pt"].write<std::vector<float>>();
        std::vector<float>& caloMET_phi = tree["caloMET_phi"].write<std::vector<float>>();
        std::vector<float>& caloMET_sumEt = tree["caloMET_sumEt"].write<std::vector<float>>();
};
