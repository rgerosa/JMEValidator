#pragma once

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"

class JMEJetAnalyzer : public JME::PhysicsObjectAnalyzer
{
    public:
        // construction/destruction
        explicit JMEJetAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~JMEJetAnalyzer();

    private:

        // member functions
        void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
        void endJob(){;}

        void computeBetaStar(const pat::Jet& jet, const std::vector<reco::Vertex>& vertices);

    private:
        // member data
        std::string   JetCorLabel_;
        std::vector<std::string> JetCorLevels_;

        edm::EDGetTokenT<std::vector<pat::Jet>> srcJet_;
        edm::EDGetTokenT<std::vector<reco::Vertex>> srcVtx_;
        edm::EDGetTokenT<std::vector<pat::Muon>> srcMuons_;

        bool          doComposition_;
        bool          doFlavor_;
        unsigned int  nJetMax_;
        double        deltaRMax_;
        double        deltaPhiMin_;
        double        deltaRPartonMax_;
        FactorizedJetCorrector* jetCorrector_;

        // Tree branches
        std::vector<int>& refpdgid = tree["refpdgid"].write<std::vector<int>>();
        std::vector<float>& refdrjt = tree["refdrjt"].write<std::vector<float>>();
        std::vector<float>& refarea = tree["refarea"].write<std::vector<float>>();
        std::vector<int8_t>&   partonFlavor = tree["partonFlavor"].write<std::vector<int8_t>>();
        std::vector<int8_t>&   hadronFlavor = tree["hadronFlavor"].write<std::vector<int8_t>>();
        std::vector<float>& jtarea = tree["jtarea"].write<std::vector<float>>();
        std::vector<float>& jtjec = tree["jtjec"].write<std::vector<float>>();
        std::vector<float>& beta = tree["beta"].write<std::vector<float>>();
        std::vector<float>& betaStar = tree["betaStar"].write<std::vector<float>>();
        std::vector<float>& betaClassic = tree["betaClassic"].write<std::vector<float>>();
        std::vector<float>& betaStarClassic = tree["betaStarClassic"].write<std::vector<float>>();
        std::vector<float>& dZ = tree["dZ"].write<std::vector<float>>();
        std::vector<float>& DRweighted = tree["DRweighted"].write<std::vector<float>>();
        std::vector<float>& fRing0 = tree["fRing0"].write<std::vector<float>>();
        std::vector<float>& fRing1 = tree["fRing1"].write<std::vector<float>>();
        std::vector<float>& fRing2 = tree["fRing2"].write<std::vector<float>>();
        std::vector<float>& fRing3 = tree["fRing3"].write<std::vector<float>>();
        std::vector<float>& fRing4 = tree["fRing4"].write<std::vector<float>>();
        std::vector<float>& fRing5 = tree["fRing5"].write<std::vector<float>>();
        std::vector<float>& fRing6 = tree["fRing6"].write<std::vector<float>>();
        std::vector<float>& fRing7 = tree["fRing7"].write<std::vector<float>>();
        std::vector<float>& fRing8 = tree["fRing8"].write<std::vector<float>>();
        std::vector<uint16_t>& nCh = tree["nCh"].write<std::vector<uint16_t>>();
        std::vector<uint16_t>& nNeutrals = tree["nNeutrals"].write<std::vector<uint16_t>>();
        std::vector<float>& ptD = tree["ptD"].write<std::vector<float>>();
        std::vector<float>& pujetid_fulldiscriminant = tree["PUJetId_fullDiscriminant"].write<std::vector<float>>();
        std::vector<int>& pujetid_cutbasedid = tree["PUJetId_cutBasedId"].write<std::vector<int>>();
        std::vector<int>& pujetid_fullid = tree["PUJetId_fullId"].write<std::vector<int>>();
        std::vector<float>& qg_tagger = tree["QGTagger_qgLikelihood"].write<std::vector<float>>();

        std::vector<float>& chargedEmEnergyFraction = tree["chargedEmEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& chargedHadronEnergyFraction = tree["chargedHadronEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& chargedMuEnergyFraction = tree["chargedMuEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& electronEnergyFraction = tree["electronEnergyFraction"].write<std::vector<float>>();


        std::vector<float>& HFEMEnergyFraction = tree["HFEMEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& HFHadronEnergyFraction = tree["HFHadronEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& hoEnergyFraction = tree["hoEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& muonEnergyFraction = tree["muonEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& neutralEmEnergyFraction = tree["neutralEmEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& neutralHadronEnergyFraction = tree["neutralHadronEnergyFraction"].write<std::vector<float>>();
        std::vector<float>& photonEnergyFraction = tree["photonEnergyFraction"].write<std::vector<float>>();

        // Calo Jet
        //std::vector<float>& emEnergyFraction = tree["emEnergyFraction"].write<std::vector<float>>();
        //std::vector<float>& emEnergyInEB = tree["emEnergyInEB"].write<std::vector<float>>();
        //std::vector<float>& emEnergyInEE = tree["emEnergyInEE"].write<std::vector<float>>();
        //std::vector<float>& emEnergyInHF = tree["emEnergyInHF"].write<std::vector<float>>();
        //std::vector<float>& energyFractionHadronic = tree["energyFractionHadronic"].write<std::vector<float>>();
};
