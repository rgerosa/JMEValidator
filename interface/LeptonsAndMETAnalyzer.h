#pragma once

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

#include <vector>

class LeptonsAndMETAnalyzer : public JME::Analyzer
{
public:
  // construction/destruction
  explicit LeptonsAndMETAnalyzer(const edm::ParameterSet& iConfig);
  virtual ~LeptonsAndMETAnalyzer();

private:
  // member functions
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup) override;
  virtual void endJob(){;}
  void recoilComputation( float &met, float &metPhi, float &Zpt, float &Zphi, float &upara, float &uperp);

private:
  // member data
  int eventCounter_; 
  int eventLimit_; 

  // Tokens
  edm::EDGetTokenT<edm::View<reco::Candidate>> srcIsoMuons_;
  edm::EDGetTokenT<std::vector<pat::MET>> srcMET_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> srcPUPPET_;

  edm::InputTag srcVtx_;
  edm::InputTag srcMuons_;
  edm::InputTag srcVMCHSTAND_;
  edm::InputTag srcVMNHSTAND_;
  edm::InputTag srcVMPhSTAND_;
  edm::InputTag srcVMPUSTAND_;
  edm::InputTag srcVMNHPFWGT_;
  edm::InputTag srcVMPhPFWGT_;
  edm::InputTag srcVMCHPUPPI_;
  edm::InputTag srcVMNHPUPPI_;
  edm::InputTag srcVMPhPUPPI_;
  edm::InputTag srcVMCHNOMUONPUPPI_;
  edm::InputTag srcVMNHNOMUONPUPPI_;
  edm::InputTag srcVMPhNOMUONPUPPI_;

  // Tree branches
  ULong64_t& run = tree["run"].write<ULong64_t>();
  ULong64_t& lumiBlock = tree["lumi"].write<ULong64_t>();
  ULong64_t& event = tree["evt"].write<ULong64_t>();
  float& pfMET = tree["pfMET"].write<float>();
  float& pfMETphi = tree["pfMETphi"].write<float>();
  float& puppET = tree["puppET"].write<float>();
  float& puppETphi = tree["puppETphi"].write<float>();

  // drellyan analysis for MET
  float& nZcands = tree["nZcands"].write<float>();
  float& ZpT = tree["ZpT"].write<float>();
  float& Zphi = tree["Zphi"].write<float>();
  float& Zmass = tree["Zmass"].write<float>();
  float& pfMET_uPara = tree["pfMET_uPara"].write<float>();
  float& pfMET_uPerp = tree["pfMET_uPerp"].write<float>();
  float& puppET_uPara = tree["puppET_uPara"].write<float>();
  float& puppET_uPerp = tree["puppET_uPerp"].write<float>();

  // Muon isolation analysis
  std::vector<float> & mupt  = tree["mupt" ].write<std::vector<float> >();
  std::vector<float> & mueta = tree["mueta"].write<std::vector<float> >();
  std::vector<float> & muphi = tree["muphi"].write<std::vector<float> >();
  std::vector<float> & mue   = tree["mue"  ].write<std::vector<float> >();
  std::vector<float> & muIsoRAW   = tree["muIsoRAW"]  .write<std::vector<float> >();
  std::vector<float> & muIsoSTAND = tree["muIsoSTAND"].write<std::vector<float> >();
  std::vector<float> & muIsoPFWGT = tree["muIsoPFWGT"].write<std::vector<float> >();
  std::vector<float> & muIsoPUPPI = tree["muIsoPUPPI"].write<std::vector<float> >();
  std::vector<float> & muIsoNOMUONPUPPI = tree["muIsoNOMUONPUPPI"].write<std::vector<float> >();
  std::vector<float> & muIso_CH     	 = tree["muIso_CH"           ].write<std::vector<float> >();
  std::vector<float> & muIso_NU     	 = tree["muIso_NU"	   ].write<std::vector<float> >();     
  std::vector<float> & muIso_PH     	 = tree["muIso_PH"	   ].write<std::vector<float> >();     
  std::vector<float> & muIso_PU     	 = tree["muIso_PU"	   ].write<std::vector<float> >();     
  std::vector<float> & muIso_NUPFW  	 = tree["muIso_NUPFW"        ].write<std::vector<float> >();
  std::vector<float> & muIso_PHPFW  	 = tree["muIso_PHPFW"        ].write<std::vector<float> >();
  std::vector<float> & muIso_CHPUPPI	 = tree["muIso_CHPUPPI"      ].write<std::vector<float> >();
  std::vector<float> & muIso_NUPUPPI	 = tree["muIso_NUPUPPI"      ].write<std::vector<float> >();
  std::vector<float> & muIso_PHPUPPI	 = tree["muIso_PHPUPPI"      ].write<std::vector<float> >();
  std::vector<float> & muIso_CHNOMUONPUPPI = tree["muIso_CHNOMUONPUPPI"].write<std::vector<float> >();
  std::vector<float> & muIso_NUNOMUONPUPPI = tree["muIso_NUNOMUONPUPPI"].write<std::vector<float> >();
  std::vector<float> & muIso_PHNOMUONPUPPI = tree["muIso_PHNOMUONPUPPI"].write<std::vector<float> >();


};
