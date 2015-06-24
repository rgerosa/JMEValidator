#pragma once

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

class EventAnalyzer: public JME::Analyzer {
    public:
        explicit EventAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~EventAnalyzer();

        virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    private:
        edm::EDGetTokenT<double> rhoToken_;
        edm::EDGetTokenT<std::vector<double>> rhosToken_;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfoToken_;
        edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
        edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;

    private:
        // Tree members
        float& rho_ = tree["rho"].write<float>();
        std::vector<float>& rhos_ = tree["rhos"].write<std::vector<float>>();
        ULong64_t& npv = tree["npv"].write<ULong64_t>();
        ULong64_t& run = tree["run"].write<ULong64_t>();
        ULong64_t& lumiBlock = tree["lumi"].write<ULong64_t>();
        ULong64_t& event = tree["evt"].write<ULong64_t>();
        std::vector<int>& npus = tree["npus"].write<std::vector<int>>();
        std::vector<float>& tnpus = tree["tnpus"].write<std::vector<float>>();
        std::vector<int>& bxns = tree["bxns"].write<std::vector<int>>();
        std::vector<int>& bunchSpacings_ = tree["bunch_spacings"].write<std::vector<int>>();
        float& pthat_ = tree["pthat"].write<float>();
        float& weight_ = tree["weight"].write<float>();
        std::vector<float>& sumpt_lowpt_ = tree["pu_sumpt_lowpt"].write<std::vector<float>>();
        std::vector<float>& sumpt_highpt_ = tree["pu_sumpt_highpt"].write<std::vector<float>>();
        std::vector<int>& ntrks_lowpt_ = tree["pu_ntrks_lowpt"].write<std::vector<int>>();
        std::vector<int>& ntrks_highpt_ = tree["pu_ntrks_highpt"].write<std::vector<int>>();
        std::vector<std::vector<float>>& pu_zpositions_ = tree["pu_zpositions"].write<std::vector<std::vector<float>>>();
        int& npuIT = tree["npuIT"].write<int>();
        int& npuOOT = tree["npuOOT"].write<int>();
        int& nTrueInt = tree["nTrueInt"].write<int>();

        float& alphaQCD_ = tree["alphaQCD"].write<float>();
        float& alphaQED_ = tree["alphaQED"].write<float>();
        float& qScale_ = tree["qScale"].write<float>();
        std::pair<int, int>& pdfID_ = tree["pdf_id"].write<std::pair<int, int>>();
        std::pair<float, float>& pdfX_ = tree["pdf_x"].write<std::pair<float, float>>();
        int& nMEPartons_ = tree["nMEPartons"].write<int>();
        int& nMEPartonsFiltered_ = tree["nMEPartonsFiltered"].write<int>();
};
