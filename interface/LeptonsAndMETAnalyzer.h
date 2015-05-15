#pragma once

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

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

};
