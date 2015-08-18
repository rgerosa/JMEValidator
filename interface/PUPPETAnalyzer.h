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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

class PUPPETAnalyzer : public JME::Analyzer {
    public:
        // construction/destruction
        explicit PUPPETAnalyzer(const edm::ParameterSet& iConfig);
        virtual ~PUPPETAnalyzer();

    private:

        // member functions
        void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
        pat::MET getUncorrectedRecoil(const pat::MET& input);

    private:

        double dRgenMatching_;
	bool   isMC_;
        TVector2 RecoilVec;
        TVector2 BosonVec;

	edm::InputTag srcPFMet_;
	edm::InputTag srcPFPuppiMet_;
	edm::InputTag srcPFCHSMet_;

	edm::InputTag srcRecoilPFMet_;
	edm::InputTag srcRecoilPFCHSMet_;
	edm::InputTag srcRecoilPFPuppiMet_;
	edm::InputTag srcRecoilPFPuppiMet_ChargedPV_;
	edm::InputTag srcRecoilPFPuppiMet_ChargedPU_;
	edm::InputTag srcRecoilPFPuppiMet_NeutralPV_;
	edm::InputTag srcRecoilPFPuppiMet_NeutralPU_;
	edm::InputTag srcRecoilPFChargedPVNeutralPVPUJetID_;
	edm::InputTag srcRecoilPFChargedPUNeutralPUPUJetID_;
	edm::InputTag srcRecoilPFChargedPVNeutralPV_;

	edm::InputTag srcJet_;
	edm::InputTag srcJetPF_;

	edm::InputTag srcZboson_;
	edm::InputTag srcLeptons_;
	edm::InputTag srcVertex_;
	edm::InputTag srcMVAMet_;

        edm::InputTag srcGenMet_;

        edm::InputTag srcGenJets_;
        edm::InputTag srcGenJetsCleaned_;

        edm::InputTag srcGenParticles_;
	edm::InputTag srcGenEventInfo_;


	edm::InputTag srcMetFiltersBits_;
	edm::InputTag srcTriggerBits_;
	edm::InputTag srcTriggerPrescales_;

	edm::EDGetTokenT<std::vector<pat::MET>> srcPFMetToken_;
	edm::EDGetTokenT<std::vector<pat::MET>> srcPFCHSMetToken_;
	edm::EDGetTokenT<std::vector<pat::MET>> srcPFPuppiMetToken_;

        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFCHSMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_ChargedPVToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_ChargedPUToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_NeutralPVToken_;
        edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFPuppiMet_NeutralPUToken_;
	edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFChargedPVNeutralPVPUJetIDToken_;
	edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFChargedPUNeutralPUPUJetIDToken_;
	edm::EDGetTokenT<std::vector<pat::MET>> srcRecoilPFChargedPVNeutralPVToken_;

        edm::EDGetTokenT<std::vector<pat::Jet>> srcJetToken_;
        edm::EDGetTokenT<std::vector<pat::Jet>> srcJetPFToken_;
        edm::EDGetTokenT<std::vector<reco::Particle>>   srcZbosonToken_;
        edm::EDGetTokenT<reco::CandidateView>           srcLeptonsToken_;
        edm::EDGetTokenT<reco::VertexCollection>        srcVertexToken_;

	edm::EDGetTokenT<std::vector<pat::MET>>         srcMVAMetToken_;
        edm::EDGetTokenT<std::vector<pat::MET>>         srcGenMetToken_;

        edm::EDGetTokenT<std::vector<reco::GenJet>>     srcGenJetsToken_;
        edm::EDGetTokenT<pat::JetCollection>            srcGenJetsCleanedToken_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> srcGenParticlesToken_;
	edm::EDGetTokenT<GenEventInfoProduct> srcGenEventInfoToken_;

	edm::EDGetTokenT<edm::TriggerResults>  srcMetFiltersBitsToken_;
	edm::EDGetTokenT<edm::TriggerResults>  srcTriggerBitsToken_;
	edm::EDGetTokenT<pat::PackedTriggerPrescales>  srcTriggerPrescalesToken_;

        // MC weight factor
	float& eventMCWeight = tree["eventMCWeight"].write<float>(); 
	
	// Generator level boson 
        float& GenBoson_Pt_  = tree["GenBoson_Pt"].write<float>();
        float& GenBoson_Eta_ = tree["GenBoson_Eta"].write<float>();
        float& GenBoson_Phi_ = tree["GenBoson_Phi"].write<float>();
        float& GenBoson_M_   = tree["GenBoson_M"].write<float>();
        int&   GenBoson_daughter_ = tree["GenBoson_daughter"].write<int>();

	// Generator level jets
        int& NGenJets_        = tree["NGenJets"].write<int>();
        int& NGenJetsCleaned_ = tree["NGenJetsCleaned"].write<int>();
        int& NGenMatchedJets_ = tree["NGenMatchedJets"].write< int>();
        int& NGenMatchedJetsPF_ = tree["NGenMatchedJetsPF"].write< int>();

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

	// Generator level recoil

        float& GenRecoil_sumEt_ = tree["GenRecoil_sumEt"].write<float>();
        float& GenRecoil_Pt_    = tree["GenRecoil_Pt"].write<float>();
        float& GenRecoil_Phi_   = tree["GenRecoil_Phi"].write<float>();
        float& GenRecoil_PerpZ_ = tree["GenRecoil_PerpZ"].write<float>();
        float& GenRecoil_LongZ_ = tree["GenRecoil_LongZ"].write<float>();

        float& GenBoson_PerpU_  = tree["GenBoson_PerpZ"].write<float>();
        float& GenBoson_LongU_  = tree["GenBoson_LongZ"].write<float>();

	// recoils
        float& recoilPFMet_sumEt_       = tree["recoilPFMet_sumEt"].write<float>();
        float& recoilPFMet_Pt_          = tree["recoilPFMet_Pt"].write<float>();
        float& recoilPFMet_Phi_         = tree["recoilPFMet_Phi"].write<float>();
        float& recoilPFMet_PerpZ_       = tree["recoilPFMet_PerpZ"].write<float>();
        float& recoilPFMet_LongZ_       = tree["recoilPFMet_LongZ"].write<float>();
        float& recoilPFMet_Boson_PerpU_ = tree["recoilPFMet_Boson_PerpU"].write<float>();
        float& recoilPFMet_Boson_LongU_ = tree["recoilPFMet_Boson_LongU"].write<float>();

	float& PFMet_Pt_    = tree["PFMet_Pt"].write<float>();
	float& PFMet_Phi_   = tree["PFMet_Phi"].write<float>();
	float& PFMet_sumEt_ = tree["PFMet_SumEt"].write<float>();

	float& PFMet_uncorrected_Pt_    = tree["PFMet_uncorrected_Pt"].write<float>();
	float& PFMet_uncorrected_Phi_   = tree["PFMet_uncorrected_Phi"].write<float>();
	float& PFMet_uncorrected_sumEt_ = tree["PFMet_uncorrected_SumEt"].write<float>();

        float& recoilPFMet_uncorrected_sumEt_ = tree["recoilPFMet_uncorrected_sumEt"].write<float>();
        float& recoilPFMet_uncorrected_Pt_    = tree["recoilPFMet_uncorrected_Pt"].write<float>();
        float& recoilPFMet_uncorrected_Phi_   = tree["recoilPFMet_uncorrected_Phi"].write<float>();
        float& recoilPFMet_uncorrected_PerpZ_ = tree["recoilPFMet_uncorrected_PerpZ"].write<float>();
        float& recoilPFMet_uncorrected_LongZ_ = tree["recoilPFMet_uncorrected_LongZ"].write<float>();
        float& recoilPFMet_uncorrected_Boson_PerpU_ = tree["recoilPFMet_uncorrected_Boson_PerpU"].write<float>();
        float& recoilPFMet_uncorrected_Boson_LongU_ = tree["recoilPFMet_uncorrected_Boson_LongU"].write<float>();

        float& recoilPFCHSMet_sumEt_ = tree["recoilPFCHSMet_sumEt"].write<float>();
        float& recoilPFCHSMet_Pt_    = tree["recoilPFCHSMet_Pt"].write<float>();
        float& recoilPFCHSMet_Phi_   = tree["recoilPFCHSMet_Phi"].write<float>();
        float& recoilPFCHSMet_PerpZ_ = tree["recoilPFCHSMet_PerpZ"].write<float>();
        float& recoilPFCHSMet_LongZ_ = tree["recoilPFCHSMet_LongZ"].write<float>();

	float& PFCHSMet_Pt_    = tree["PFCHSMet_Pt"].write<float>();
	float& PFCHSMet_Phi_   = tree["PFCHSMet_Phi"].write<float>();
	float& PFCHSMet_sumEt_ = tree["PFCHSMet_SumEt"].write<float>();

	float& PFCHSMet_uncorrected_Pt_    = tree["PFCHSMet_uncorrected_Pt"].write<float>();
	float& PFCHSMet_uncorrected_Phi_   = tree["PFCHSMet_uncorrected_Phi"].write<float>();
	float& PFCHSMet_uncorrected_sumEt_ = tree["PFCHSMet_uncorrected_SumEt"].write<float>();

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
        float& recoilPFPuppiMet_Boson_PerpU_ = tree["recoilPFPuppiMet_Boson_PerpU"].write<float>();
        float& recoilPFPuppiMet_Boson_LongU_ = tree["recoilPFPuppiMet_Boson_LongU"].write<float>();

	float& PFPuppiMet_Pt_    = tree["PFPuppiMet_Pt"].write<float>();
	float& PFPuppiMet_Phi_   = tree["PFPuppiMet_Phi"].write<float>();
	float& PFPuppiMet_sumEt_ = tree["PFPuppiMet_SumEt"].write<float>();

	float& PFPuppiMet_uncorrected_Pt_    = tree["PFPuppiMet_uncorrected_Pt"].write<float>();
	float& PFPuppiMet_uncorrected_Phi_   = tree["PFPuppiMet_uncorrected_Phi"].write<float>();
	float& PFPuppiMet_uncorrected_sumEt_ = tree["PFPuppiMet_uncorrected_SumEt"].write<float>();

        float& recoilPFPuppiMet_uncorrected_sumEt_ = tree["recoilPFPuppiMet_uncorrected_sumEt"].write<float>();
        float& recoilPFPuppiMet_uncorrected_Pt_    = tree["recoilPFPuppiMet_uncorrected_Pt"].write<float>();
        float& recoilPFPuppiMet_uncorrected_Phi_   = tree["recoilPFPuppiMet_uncorrected_Phi"].write<float>();
        float& recoilPFPuppiMet_uncorrected_PerpZ_ = tree["recoilPFPuppiMet_uncorrected_PerpZ"].write<float>();
        float& recoilPFPuppiMet_uncorrected_LongZ_ = tree["recoilPFPuppiMet_uncorrected_LongZ"].write<float>();
        float& recoilPFPuppiMet_uncorrected_Boson_PerpU_ = tree["recoilPFPuppiMet_uncorrected_Boson_PerpU"].write<float>();
        float& recoilPFPuppiMet_uncorrected_Boson_LongU_ = tree["recoilPFPuppiMet_uncorrected_Boson_LongU"].write<float>();

        float& recoilPFPuppiMet_ChargedPV_sumEt_ = tree["recoilPFPuppiMet_ChargedPV_sumEt"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_Pt_    = tree["recoilPFPuppiMet_ChargedPV_Pt"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_Phi_   = tree["recoilPFPuppiMet_ChargedPV_Phi"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_PerpZ_ = tree["recoilPFPuppiMet_ChargedPV_PerpZ"].write<float>();
        float& recoilPFPuppiMet_ChargedPV_LongZ_ = tree["recoilPFPuppiMet_ChargedPV_LongZ"].write<float>();

        float& recoilPFPuppiMet_ChargedPU_sumEt_ = tree["recoilPFPuppiMet_ChargedPU_sumEt"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_Pt_    = tree["recoilPFPuppiMet_ChargedPU_Pt"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_Phi_   = tree["recoilPFPuppiMet_ChargedPU_Phi"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_PerpZ_ = tree["recoilPFPuppiMet_ChargedPU_PerpZ"].write<float>();
        float& recoilPFPuppiMet_ChargedPU_LongZ_ = tree["recoilPFPuppiMet_ChargedPU_LongZ"].write<float>();

        float& recoilPFPuppiMet_NeutralPV_sumEt_ = tree["recoilPFPuppiMet_NeutralPV_sumEt"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_Pt_    = tree["recoilPFPuppiMet_NeutralPV_Pt"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_Phi_   = tree["recoilPFPuppiMet_NeutralPV_Phi"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_PerpZ_ = tree["recoilPFPuppiMet_NeutralPV_PerpZ"].write<float>();
        float& recoilPFPuppiMet_NeutralPV_LongZ_ = tree["recoilPFPuppiMet_NeutralPV_LongZ"].write<float>();

        float& recoilPFPuppiMet_NeutralPU_sumEt_ = tree["recoilPFPuppiMet_NeutralPU_sumEt"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_Pt_    = tree["recoilPFPuppiMet_NeutralPU_Pt"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_Phi_   = tree["recoilPFPuppiMet_NeutralPU_Phi"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_PerpZ_ = tree["recoilPFPuppiMet_NeutralPU_PerpZ"].write<float>();
        float& recoilPFPuppiMet_NeutralPU_LongZ_ = tree["recoilPFPuppiMet_NeutralPU_LongZ"].write<float>();

        float& recoilPFMet_ChargedPVNeutralPVPUJetID_sumEt_ = tree["recoilPFMet_ChargedPVNeutralPVPUJetID_sumEt"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPVPUJetID_Pt_    = tree["recoilPFMet_ChargedPVNeutralPVPUJetID_Pt"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPVPUJetID_Phi_   = tree["recoilPFMet_ChargedPVNeutralPVPUJetID_Phi"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPVPUJetID_PerpZ_ = tree["recoilPFMet_ChargedPVNeutralPVPUJetID_PerpZ"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPVPUJetID_LongZ_ = tree["recoilPFMet_ChargedPVNeutralPVPUJetID_LongZ"].write<float>();

        float& recoilPFMet_ChargedPUNeutralPUPUJetID_sumEt_ = tree["recoilPFMet_ChargedPUNeutralPUPUJetID_sumEt"].write<float>();
        float& recoilPFMet_ChargedPUNeutralPUPUJetID_Pt_    = tree["recoilPFMet_ChargedPUNeutralPUPUJetID_Pt"].write<float>();
        float& recoilPFMet_ChargedPUNeutralPUPUJetID_Phi_   = tree["recoilPFMet_ChargedPUNeutralPUPUJetID_Phi"].write<float>();
        float& recoilPFMet_ChargedPUNeutralPUPUJetID_PerpZ_ = tree["recoilPFMet_ChargedPUNeutralPUPUJetID_PerpZ"].write<float>();
        float& recoilPFMet_ChargedPUNeutralPUPUJetID_LongZ_ = tree["recoilPFMet_ChargedPUNeutralPUPUJetID_LongZ"].write<float>();

        float& recoilPFMet_ChargedPVNeutralPV_sumEt_ = tree["recoilPFMet_ChargedPVNeutralPV_sumEt"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPV_Pt_    = tree["recoilPFMet_ChargedPVNeutralPV_Pt"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPV_Phi_   = tree["recoilPFMet_ChargedPVNeutralPV_Phi"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPV_PerpZ_ = tree["recoilPFMet_ChargedPVNeutralPV_PerpZ"].write<float>();
        float& recoilPFMet_ChargedPVNeutralPV_LongZ_ = tree["recoilPFMet_ChargedPVNeutralPV_LongZ"].write<float>();


	float& MVAMet_sumEt_ = tree["MVAMet_sumEt"].write<float>();
	float& MVAMet_Pt_    = tree["MVAMet_Pt"].write<float>();
	float& MVAMet_Phi_   = tree["MVAMet_Phi"].write<float>();
	float& MVAMet_PerpZ_ = tree["MVAMet_PerpZ"].write<float>();
	float& MVAMet_LongZ_ = tree["MVAMet_LongZ"].write<float>();

	// jets and boson
	int& NCleanedJets_ =  tree["NCleanedJets"].write<int>();
	int& NCleanedJetsPF_ =  tree["NCleanedJetsPF"].write<int>();
	int& NVertex_      =  tree["NVertex"].write<int>();

	float& LeadingJet_Pt_  = tree["LeadingJet_Pt"].write<float>();
	float& LeadingJet_Eta_ = tree["LeadingJet_Eta"].write<float>();
	float& LeadingJet_Phi_ = tree["LeadingJet_Phi"].write<float>();
	float& LeadingJet_M_   = tree["LeadingJet_M"].write<float>();

	float& TrailingJet_Pt_  = tree["TrailingJet_Pt"].write<float>();
	float& TrailingJet_Eta_ = tree["TrailingJet_Eta"].write<float>();
	float& TrailingJet_Phi_ = tree["TrailingJet_Phi"].write<float>();
	float& TrailingJet_M_   = tree["TrailingJet_M"].write<float>();

	float& LeadingJetPF_Pt_  = tree["LeadingJet_Pt"].write<float>();
	float& LeadingJetPF_Eta_ = tree["LeadingJet_Eta"].write<float>();
	float& LeadingJetPF_Phi_ = tree["LeadingJet_Phi"].write<float>();
	float& LeadingJetPF_M_   = tree["LeadingJet_M"].write<float>();

	float& TrailingJetPF_Pt_  = tree["TrailingJet_Pt"].write<float>();
	float& TrailingJetPF_Eta_ = tree["TrailingJet_Eta"].write<float>();
	float& TrailingJetPF_Phi_ = tree["TrailingJet_Phi"].write<float>();
	float& TrailingJetPF_M_   = tree["TrailingJet_M"].write<float>();

	float& Boson_Pt_     =  tree["Boson_Pt"].write<float>();
	float& Boson_Phi_    =  tree["Boson_Phi"].write<float>();
	float& Boson_Eta_    =  tree["Boson_Eta"].write<float>();
	float& Boson_M_      =  tree["Boson_M"].write<float>();
	int& Boson_daughter_ =  tree["Boson_daughter"].write<int>();

	float& LeadingLepton_Pt_     =  tree["LeadingLepton_Pt"].write<float>();
	float& LeadingLepton_Phi_    =  tree["LeadingLepton_Phi"].write<float>();
	float& LeadingLepton_Eta_    =  tree["LeadingLepton_Eta"].write<float>();
	float& LeadingLepton_M_      =  tree["LeadingLepton_M"].write<float>();

	float& TrailingLepton_Pt_     =  tree["TrailingLepton_Pt"].write<float>();
	float& TrailingLepton_Phi_    =  tree["TrailingLepton_Phi"].write<float>();
	float& TrailingLepton_Eta_    =  tree["TrailingLepton_Eta"].write<float>();
	float& TrailingLepton_M_      =  tree["TrailingLepton_M"].write<float>();

        std::vector<float>& AllJets_Pt_  = tree["AllJets_Pt"].write<std::vector<float>>();
        std::vector<float>& AllJets_Eta_ = tree["AllJets_Eta"].write<std::vector<float>>();
        std::vector<float>& AllJets_Phi_ = tree["AllJets_Phi"].write<std::vector<float>>();
        std::vector<float>& AllJets_M_   = tree["AllJets_M"].write<std::vector<float>>();

        std::vector<float>& AllJetsPF_Pt_  = tree["AllJetsPF_Pt"].write<std::vector<float>>();
        std::vector<float>& AllJetsPF_Eta_ = tree["AllJetsPF_Eta"].write<std::vector<float>>();
        std::vector<float>& AllJetsPF_Phi_ = tree["AllJetsPF_Phi"].write<std::vector<float>>();
        std::vector<float>& AllJetsPF_M_   = tree["AllJetsPF_M"].write<std::vector<float>>();

        std::vector<float>& GenMatchedJets_Pt_  = tree["GenMatchedJets_Pt"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_Eta_ = tree["GenMatchedJets_Eta"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_Phi_ = tree["GenMatchedJets_Phi"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJets_M_   = tree["GenMatchedJets_M"].write<std::vector<float>>();

        std::vector<float>& GenMatchedJetsPF_Pt_  = tree["GenMatchedJets_Pt"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJetsPF_Eta_ = tree["GenMatchedJets_Eta"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJetsPF_Phi_ = tree["GenMatchedJets_Phi"].write<std::vector<float>>();
        std::vector<float>& GenMatchedJetsPF_M_   = tree["GenMatchedJets_M"].write<std::vector<float>>();

	/// met filter
	int& flag_HBHENoiseFilter_      =  tree["flag_HBHENoiseFilter"].write<int>();
	int& flag_CSCTightHaloFilter_   =  tree["flag_CSCTightHaloFilter"].write<int>();
        int& flag_hcalLaserEventFilter_ =  tree["flag_hcalLaserEventFilter"].write<int>();
	int& flag_EcalDeadCellTriggerPrimitiveFilter_ = tree["flag_EcalDeadCellTriggerPrimitiveFilter"].write<int>();
	int& flag_goodVertices_          = tree["flag_goodVertices"].write<int>();
	int& flag_trackingFailureFilter_ = tree["flag_trackingFailureFilter"].write<int>();
	int& flag_eeBadScFilter_         = tree["flag_eeBadScFilter"].write<int>();
	int& flag_ecalLaserCorrFilter_   = tree["flag_ecalLaserCorrFilter"].write<int>();
	int& flag_trkPOGFilters_         =  tree["flag_trkPOGFilters"].write<int>();  
	int& flag_trkPOG_manystripclus53X_    = tree["flag_trkPOG_manystripclus53X"].write<int>();
	int& flag_trkPOG_toomanystripclus53X_ = tree["flag_trkPOG_toomanystripclus53X"].write<int>();
	int& flag_trkPOG_logErrorTooManyClusters_ = tree["flag_trkPOG_logErrorTooManyClusters"].write<int>();
	int& flag_METFilters_ = tree["flag_METFilters"].write<int>();


	std::vector<std::string>& DoubleMuPaths_    = tree["DoubleMuPaths"].write<std::vector<std::string> >();
	std::vector<std::string>& DoubleElePaths_   = tree["DoubleElePaths"].write<std::vector<std::string> >();
	std::vector<std::string>& DoubleTauPaths_   = tree["DoubleTauPaths"].write<std::vector<std::string> >();

};

#endif
