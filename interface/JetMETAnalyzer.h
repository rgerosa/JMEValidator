#pragma once

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"

class JetMETAnalyzer : public JME::PhysicsObjectAnalyzer
{
    public:
        // construction/destruction
        explicit JetMETAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~JetMETAnalyzer();

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
        int& nref = tree["nref"].write<int>();
        std::vector<int>& refrank = tree["refrank"].write<std::vector<int>>();
        std::vector<int>& refpdgid = tree["refpdgid"].write<std::vector<int>>();
        std::vector<float>& refdrjt = tree["refdrjt"].write<std::vector<float>>();
        std::vector<float>& refe = tree["refe"].write<std::vector<float>>();
        std::vector<float>& refpt = tree["refpt"].write<std::vector<float>>();
        std::vector<float>& refeta = tree["refeta"].write<std::vector<float>>();
        std::vector<float>& refphi = tree["refphi"].write<std::vector<float>>();
        std::vector<float>& refm = tree["refm"].write<std::vector<float>>();
        std::vector<float>& refy = tree["refy"].write<std::vector<float>>();
        std::vector<float>& refarea = tree["refarea"].write<std::vector<float>>();
        std::vector<int>&   partonFlavor = tree["partonFlavor"].write<std::vector<int>>();
        std::vector<int>&   hadronFlavor = tree["hadronFlavor"].write<std::vector<int>>();
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
        std::vector<float>& nCh = tree["nCh"].write<std::vector<float>>();
        std::vector<float>& nNeutrals = tree["nNeutrals"].write<std::vector<float>>();
        std::vector<float>& ptD = tree["ptD"].write<std::vector<float>>();
        std::vector<bool>& isMatched = tree["isMatched"].write<std::vector<bool>>();
        std::vector<std::pair<std::string, float>>& btagDiscri = tree["btagDiscriminators"].write<std::vector<std::pair<std::string, float>>>();
        std::vector<float>& pujetid_fulldiscriminant = tree["PUJetId_fullDiscriminant"].write<std::vector<float>>();
        std::vector<int>& pujetid_cutbasedid = tree["PUJetId_cutBasedId"].write<std::vector<int>>();
        std::vector<int>& pujetid_fullid = tree["PUJetId_fullId"].write<std::vector<int>>();
        std::vector<float>& qg_tagger = tree["QGTagger_qgLikelihood"].write<std::vector<float>>();
};
