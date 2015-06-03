#pragma once

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

#include <vector>
#include <string>

class HLTAnalyzer: public JME::Analyzer {
    public:
        explicit HLTAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~HLTAnalyzer();

        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    private:
        edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
        edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;

    private:
        // Tree members
        std::vector<std::string>& paths = tree["paths"].write<std::vector<std::string>>();
        std::vector<int>& prescales = tree["prescales"].write<std::vector<int>>();

        // HLT objects
        std::vector<float>& objects_pt = tree["objects_pt"].write<std::vector<float>>();
        std::vector<float>& objects_eta = tree["objects_eta"].write<std::vector<float>>();
        std::vector<float>& objects_phi = tree["objects_phi"].write<std::vector<float>>();
        std::vector<float>& objects_e = tree["objects_e"].write<std::vector<float>>();

        std::vector<std::vector<int>>& objects_paths = tree["objects_paths"].write<std::vector<std::vector<int>>>();
        std::vector<std::vector<bool>>& objects_paths_isl3 = tree["objects_paths_isl3"].write<std::vector<std::vector<bool>>>();
        std::vector<std::vector<bool>>& objects_paths_islast = tree["objects_paths_islast"].write<std::vector<std::vector<bool>>>();
};
