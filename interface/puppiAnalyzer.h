#pragma once

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

class puppiAnalyzer : public JME::Analyzer
{
public:
  // construction/destruction
  explicit puppiAnalyzer(const edm::ParameterSet& iConfig);
  virtual ~puppiAnalyzer();

private:
  // member functions
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup) override;
  virtual void endJob(){;}

private:
  // member data
  int eventCounter_; 
  int eventLimit_; 

  // Tokens
  edm::EDGetTokenT<double> nAlgosToken_;
  edm::EDGetTokenT<std::vector<double>> rawAlphasToken_;
  edm::EDGetTokenT<std::vector<double>> alphasToken_;
  edm::EDGetTokenT<std::vector<double>> alphasMedToken_;
  edm::EDGetTokenT<std::vector<double>> alphasRmsToken_;
  edm::EDGetTokenT<reco::CandidateView> packedPFCandidatesToken_;

  // Tree branches
  float& nalgos = tree["nalgos"].write<float>();
  std::vector<float>& px = tree["px"].write<std::vector<float>>();
  std::vector<float>& py = tree["py"].write<std::vector<float>>();
  std::vector<float>& pz = tree["pz"].write<std::vector<float>>();
  std::vector<float>& e = tree["e"].write<std::vector<float>>();
  std::vector<float>& alphas = tree["alphas"].write<std::vector<float>>();
  std::vector<float>& thealphas = tree["thealphas"].write<std::vector<float>>();
  std::vector<float>& thealphasmed = tree["thealphasmed"].write<std::vector<float>>();
  std::vector<float>& thealphasrms = tree["thealphasrms"].write<std::vector<float>>();
  std::vector<float>& id = tree["id"].write<std::vector<float>>();
  std::vector<float>& charge = tree["charge"].write<std::vector<float>>();
  std::vector<float>& fromPV = tree["fromPV"].write<std::vector<float>>();
};
