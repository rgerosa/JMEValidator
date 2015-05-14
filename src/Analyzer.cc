#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "JMEAnalysis/JMEValidator/interface/Analyzer.h"

#include <TTree.h>

JME::Analyzer::Analyzer(const edm::ParameterSet& iConfig)
    : moduleLabel_(iConfig.getParameter<std::string>("@module_label")) {

        if (iConfig.existsAs<std::string>("treeName"))
            treeName_ = iConfig.getParameter<std::string>("treeName");
        else
            treeName_ = "t";
    }

JME::Analyzer::~Analyzer() {
    // Empty
}

void JME::Analyzer::beginJob()
{
    edm::Service<TFileService> fs;
    if (!fs)
        throw edm::Exception(edm::errors::Configuration, "TFileService missing from configuration!");

    TTree* tree_ = fs->make<TTree>(treeName_.c_str(), treeName_.c_str());
    tree.init(tree_);
}
