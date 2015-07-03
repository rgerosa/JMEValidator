#pragma once

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
<<<<<<< HEAD
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "JMEAnalysis/TreeWrapper/interface/TreeWrapper.h"

namespace JME {
  class Analyzer : public edm::EDAnalyzer {
=======
#include "JMEAnalysis/TreeWrapper/interface/TreeWrapper.h"


namespace JME {
    class Analyzer : public edm::EDAnalyzer {
>>>>>>> origin
        public:
            // construction/destruction
            explicit Analyzer(const edm::ParameterSet& iConfig);
            virtual ~Analyzer();

        protected:

            // member functions
<<<<<<< HEAD
            virtual void beginJob() ;
=======
            virtual void beginJob() override;
>>>>>>> origin
            virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) = 0;

        protected:

            // member data
            std::string moduleLabel_;

            // tree
            std::string treeName_;
            ROOT::TreeWrapper tree;
    };
}
