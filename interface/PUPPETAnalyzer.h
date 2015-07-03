#ifndef PUPPETAnalyzer_H
#define PUPPETAnalyzer_H

#include "JMEAnalysis/JMEValidator/interface/PhysicsObjectAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

class PUPPETAnalyzer : public JME::Analyzer
{
    public:
        // construction/destruction
        explicit PUPPETAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~PUPPETAnalyzer();

    private:

        // member functions
        void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);

    private:

        double dRgenMatching_;
        TVector2 RecoilVec;

	edm::InputTag srcRecoilPFMet_;
	edm::InputTag srcRecoilPFCHSMet_;
	edm::InputTag srcRecoilPFPuppiMet_;
	edm::InputTag srcRecoilPFPuppiMet_ChargedPV_;
	edm::InputTag srcRecoilPFPuppiMet_ChargedPU_;
	edm::InputTag srcRecoilPFPuppiMet_NeutralPV_;
	edm::InputTag srcRecoilPFPuppiMet_NeutralPU_;
	edm::InputTag srcJet_;
	edm::InputTag srcZboson_;
	edm::InputTag srcVertex_;
	edm::InputTag srcMVAMet_;
        edm::InputTag srcGenMet_;
        edm::InputTag srcGenJets_;
        edm::InputTag srcGenJetsCleaned_;
        edm::InputTag srcGenParticles_;

        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFCHSMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_ChargedPVToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_ChargedPUToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_NeutralPVToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_NeutralPUToken_;
        edm::EDGetTokenT<std::vector<pat::Jet>> srcJetToken_;
        edm::EDGetTokenT<std::vector<reco::Particle>>   srcZbosonToken_;
        edm::EDGetTokenT<reco::VertexCollection>        srcVertexToken_;
	edm::EDGetTokenT<std::vector<pat::MET>>         srcMVAMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>>         srcGenMetToken_;
        edm::EDGetTokenT<std::vector<reco::GenJet>>     srcGenJetsToken_;
        edm::EDGetTokenT<pat::JetCollection>            srcGenJetsCleanedToken_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> srcGenParticlesToken_;

        // Tree branches
        float& GenZ_Pt_  = tree["GenZ_Pt"].write<float>();
        float& GenZ_Eta_ = tree["GenZ_Eta"].write<float>();
        float& GenZ_Phi_ = tree["GenZ_Phi"].write<float>();
        float& GenZ_M_   = tree["GenZ_M"].write<float>();
        int& GenZ_daughter_ = tree["GenZ_daughter"].write<int>();

        int& NGenJets_        = tree["NGenJets"].write<int>();
        int& NGenJetsCleaned_ = tree["NGenJetsCleaned"].write<int>();
        int& NGenMatchedJets_ = tree["NGenMatchedJets"].write< int>();

        float& GenLeadingJet_Pt_  = tree["GenLeadingJet_Pt"].write<float>();
        float& GenLeadingJet_Eta_ = tree["GenLeadingJet_Eta"].write<float>();
        float& GenLeadingJet_Phi_ = tree["GenLeadingJet_Phi"].write<float>();
        float& GenLeadingJet_M_   = tree["GenLeadingJet_M"].write<float>();

        float& GenLeadingJetCleaned_Pt_  = tree["GenLeadingJetCleaned_Pt"].write<float>();
        float& GenLeadingJetCleaned_Eta_ = tree["GenLeadingJetCleaned_Eta"].write<float>();
        float& GenLeadingJetCleaned_Phi_ = tree["GenLeadingJetCleaned_Phi"].write<float>();
        float& GenLeadingJetCleaned_M_   = tree["GenLeadingJetCleaned_M"].write<float>();

        float& GenTrailingJet_Pt_  = tree["GenTrailingJet_Pt"].write<float>();
        float& GenTrailingJet_Eta_ = tree["GenTrailingJet_Eta"].write<float>();
        float& GenTrailingJet_Phi_ = tree["GenTrailingJet_Phi"].write<float>();
        float& GenTrailingJet_M_   = tree["GenTrailingJet_M"].write<float>();

        float& GenTrailingJetCleaned_Pt_  = tree["GenTrailingJetCleaned_Pt"].write<float>();
        float& GenTrailingJetCleaned_Eta_ = tree["GenTrailingJetCleaned_Eta"].write<float>();
        float& GenTrailingJetCleaned_Phi_ = tree["GenTrailingJetCleaned_Phi"].write<float>();
        float& GenTrailingJetCleaned_M_   = tree["GenTrailingJetCleaned_M"].write<float>();

        float& GenRecoil_sumEt_ = tree["GenRecoil_sumEt"].write<float>();
        float& GenRecoil_Pt_    = tree["GenRecoil_Pt"].write<float>();
        float& GenRecoil_Phi_   = tree["GenRecoil_Phi"].write<float>();
        float& GenRecoil_PerpZ_ = tree["GenRecoil_PerpZ"].write<float>();
        float& GenRecoil_LongZ_ = tree["GenRecoil_LongZ"].write<float>();

        float& recoilPFMet_sumEt_ = tree["recoilPFMet_sumEt"].write<float>();
        float& recoilPFMet_Pt_    = tree["recoilPFMet_Pt"].write<float>();
        float& recoilPFMet_Phi_   = tree["recoilPFMet_Phi"].write<float>();
        float& recoilPFMet_PerpZ_ = tree["recoilPFMet_PerpZ"].write<float>();
        float& recoilPFMet_LongZ_ = tree["recoilPFMet_LongZ"].write<float>();

        float& recoilPFMet_uncorrected_sumEt_ = tree["recoilPFMet_uncorrected_sumEt"].write<float>();
        float& recoilPFMet_uncorrected_Pt_    = tree["recoilPFMet_uncorrected_Pt"].write<float>();
        float& recoilPFMet_uncorrected_Phi_   = tree["recoilPFMet_uncorrected_Phi"].write<float>();
        float& recoilPFMet_uncorrected_PerpZ_ = tree["recoilPFMet_uncorrected_PerpZ"].write<float>();
        float& recoilPFMet_uncorrected_LongZ_ = tree["recoilPFMet_uncorrected_LongZ"].write<float>();

        float& recoilPFCHSMet_sumEt_ = tree["recoilPFCHSMet_sumEt"].write<float>();
        float& recoilPFCHSMet_Pt_    = tree["recoilPFCHSMet_Pt"].write<float>();
        float& recoilPFCHSMet_Phi_   = tree["recoilPFCHSMet_Phi"].write<float>();
        float& recoilPFCHSMet_PerpZ_ = tree["recoilPFCHSMet_PerpZ"].write<float>();
        float& recoilPFCHSMet_LongZ_ = tree["recoilPFCHSMet_LongZ"].write<float>();

        float& recoilPFCHSMet_uncorrected_sumEt_ = tree["recoilPFCHSMet_uncorrected_sumEt"].write<float>();
        float& recoilPFCHSMet_uncorrected_Pt_    = tree["recoilPFCHSMet_uncorrected_Pt"].write<float>();
        float& recoilPFCHSMet_uncorrected_Phi_   = tree["recoilPFCHSMet_uncorrected_Phi"].write<float>();
        float& recoilPFCHSMet_uncorrected_PerpZ_ = tree["recoilPFCHSMet_uncorrected_PerpZ"].write<float>();
        float& recoilPFCHSMet_uncorrected_LongZ_ = tree["recoilPFCHSMet_uncorrected_LongZ"].write<float>();

        float& recoilPFPuppiMet_sumEt_ = tree["recoilPFPuppiMet_sumEt"].write<float>();
        float& recoilPFPuppiMet_Pt_    = tree["recoilPFPuppiMet_Pt"].write<float>();
        float& recoilPFPuppiMet_Phi_   = tree["recoilPFPuppiMet_Phi"].write<float>();
        float& recoilPFPuppiMet_PerpZ_ = tree["recoilPFPuppiMet_PerpZ"].write<float>();
        float& recoilPFPuppiMet_LongZ_ = tree["recoilPFPuppiMet_LongZ"].write<float>();

        float& recoilPFPuppiMet_uncorrected_sumEt_ = tree["recoilPFPuppiMet_uncorrected_sumEt"].write<float>();
        float& recoilPFPuppiMet_uncorrected_Pt_    = tree["recoilPFPuppiMet_uncorrected_Pt"].write<float>();
        float& recoilPFPuppiMet_uncorrected_Phi_   = tree["recoilPFPuppiMet_uncorrected_Phi"].write<float>();
        float& recoilPFPuppiMet_uncorrected_PerpZ_ = tree["recoilPFPuppiMet_uncorrected_PerpZ"].write<float>();
        float& recoilPFPuppiMet_uncorrected_LongZ_ = tree["recoilPFPuppiMet_uncorrected_LongZ"].write<float>();

        float& recoilPFPuppiMet_ChargedPV_sumEt_ = tree["recoilPFPuppiMet_ChargedPV_sumEt"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_Pt_    = tree["recoilPFPuppiMet_ChargedPV_Pt"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_Phi_   = tree["recoilPFPuppiMet_ChargedPV_Phi"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_PerpZ_ = tree["recoilPFPuppiMet_ChargedPV_PerpZ"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_LongZ_ = tree["recoilPFPuppiMet_ChargedPV_LongZ"].write<float>();

        float& recoilPFPuppiMet_ChargedPV_uncorrected_sumEt_ = tree["recoilPFPuppiMet_ChargedPV_uncorrected_sumEt"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_uncorrected_Pt_    = tree["recoilPFPuppiMet_ChargedPV_uncorrected_Pt"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_uncorrected_Phi_   = tree["recoilPFPuppiMet_ChargedPV_uncorrected_Phi"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_uncorrected_PerpZ_ = tree["recoilPFPuppiMet_ChargedPV_uncorrected_PerpZ"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_uncorrected_LongZ_ = tree["recoilPFPuppiMet_ChargedPV_uncorrected_LongZ"].write<float>();

        float& recoilPFPuppiMet_ChargedPU_sumEt_ = tree["recoilPFPuppiMet_ChargedPU_sumEt"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_Pt_    = tree["recoilPFPuppiMet_ChargedPU_Pt"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_Phi_   = tree["recoilPFPuppiMet_ChargedPU_Phi"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_PerpZ_ = tree["recoilPFPuppiMet_ChargedPU_PerpZ"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_LongZ_ = tree["recoilPFPuppiMet_ChargedPU_LongZ"].write<float>();

        float& recoilPFPuppiMet_ChargedPU_uncorrected_sumEt_ = tree["recoilPFPuppiMet_ChargedPU_uncorrected_sumEt"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_uncorrected_Pt_    = tree["recoilPFPuppiMet_ChargedPU_uncorrected_Pt"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_uncorrected_Phi_   = tree["recoilPFPuppiMet_ChargedPU_uncorrected_Phi"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_uncorrected_PerpZ_ = tree["recoilPFPuppiMet_ChargedPU_uncorrected_PerpZ"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_uncorrected_LongZ_ = tree["recoilPFPuppiMet_ChargedPU_uncorrected_LongZ"].write<float>();

        float& recoilPFPuppiMet_NeutralPV_sumEt_ = tree["recoilPFPuppiMet_NeutralPV_sumEt"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_Pt_    = tree["recoilPFPuppiMet_NeutralPV_Pt"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_Phi_   = tree["recoilPFPuppiMet_NeutralPV_Phi"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_PerpZ_ = tree["recoilPFPuppiMet_NeutralPV_PerpZ"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_LongZ_ = tree["recoilPFPuppiMet_NeutralPV_LongZ"].write<float>();

        float& recoilPFPuppiMet_NeutralPV_uncorrected_sumEt_ = tree["recoilPFPuppiMet_NeutralPV_uncorrected_sumEt"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_uncorrected_Pt_    = tree["recoilPFPuppiMet_NeutralPV_uncorrected_Pt"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_uncorrected_Phi_   = tree["recoilPFPuppiMet_NeutralPV_uncorrected_Phi"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_uncorrected_PerpZ_ = tree["recoilPFPuppiMet_NeutralPV_uncorrected_PerpZ"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_uncorrected_LongZ_ = tree["recoilPFPuppiMet_NeutralPV_uncorrected_LongZ"].write<float>();

        float& recoilPFPuppiMet_NeutralPU_sumEt_ = tree["recoilPFPuppiMet_NeutralPU_sumEt"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_Pt_    = tree["recoilPFPuppiMet_NeutralPU_Pt"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_Phi_   = tree["recoilPFPuppiMet_NeutralPU_Phi"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_PerpZ_ = tree["recoilPFPuppiMet_NeutralPU_PerpZ"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_LongZ_ = tree["recoilPFPuppiMet_NeutralPU_LongZ"].write<float>();

        float& recoilPFPuppiMet_NeutralPU_uncorrected_sumEt_ = tree["recoilPFPuppiMet_NeutralPU_uncorrected_sumEt"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_uncorrected_Pt_    = tree["recoilPFPuppiMet_NeutralPU_uncorrected_Pt"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_uncorrected_Phi_   = tree["recoilPFPuppiMet_NeutralPU_uncorrected_Phi"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_uncorrected_PerpZ_ = tree["recoilPFPuppiMet_NeutralPU_uncorrected_PerpZ"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_uncorrected_LongZ_ = tree["recoilPFPuppiMet_NeutralPU_uncorrected_LongZ"].write<float>();

	float& MVAMet_sumEt_ = tree["MVAMet_sumEt"].write<float>();
	float& MVAMet_Pt_    = tree["MVAMet_Pt"].write<float>();
	float& MVAMet_Phi_   = tree["MVAMet_Phi"].write<float>();
	float& MVAMet_PerpZ_ = tree["MVAMet_PerpZ"].write<float>();
	float& MVAMet_LongZ_ = tree["MVAMet_LongZ"].write<float>();

	int& NCleanedJets_ =  tree["NCleanedJets"].write<int>();
	int& NVertex_      =  tree["NVertex"].write<int>();

	float& LeadingJet_Pt_  = tree["LeadingJet_Pt"].write<float>();
	float& LeadingJet_Eta_ = tree["LeadingJet_Eta"].write<float>();
	float& LeadingJet_Phi_ = tree["LeadingJet_Phi"].write<float>();
	float& LeadingJet_M_   = tree["LeadingJet_M"].write<float>();

	float& TrailingJet_Pt_  = tree["TrailingJet_Pt"].write<float>();
	float& TrailingJet_Eta_ = tree["TrailingJet_Eta"].write<float>();
	float& TrailingJet_Phi_ = tree["TrailingJet_Phi"].write<float>();
	float& TrailingJet_M_   = tree["TrailingJet_M"].write<float>();

	float& Zboson_Pt_     =  tree["Zboson_Pt"].write<float>();
	float& Zboson_Phi_    =  tree["Zboson_Phi"].write<float>();
	float& Zboson_Eta_    =  tree["Zboson_Eta"].write<float>();
	float& Zboson_M_      =  tree["Zboson_M"].write<float>();
	int& Zboson_daughter_ =  tree["Zboson_daughter"].write<int>();

        std::vector<float>& AllJets_Pt_  = tree["AllJets_Pt"].write<std::vector<float>>();
        std::vector<float>& AllJets_Eta_ = tree["AllJets_Eta"].write<std::vector<float>>();
        std::vector<float>& AllJets_Phi_ = tree["AllJets_Phi"].write<std::vector<float>>();
        std::vector<float>& AllJets_M_   = tree["AllJets_M"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_Pt_  = tree["GenMatchedJets_Pt"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_Eta_ = tree["GenMatchedJets_Eta"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_Phi_ = tree["GenMatchedJets_Phi"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_M_   = tree["GenMatchedJets_M"].write<std::vector<float>>();

};

#endif
