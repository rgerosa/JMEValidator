#include "TTree.h"
#include "../interface/GBRTrainer.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
//#include "Cintex/Cintex.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <cassert>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "../interface/applyphi.h"

using namespace std;




int main(int argc, char* argv[] ) {
	// config file einlesen
  boost::property_tree::ptree pt;
  try 
  {
    std::ifstream in(argv[1]);
    std::stringstream ss;
    ss << in.rdbuf();      
    boost::property_tree::json_parser::read_json(ss, pt);
  } catch (std::exception const& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  std::vector<boost::property_tree::ptree> trainingProperties;
  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("activeTrainings"))
  {
    assert(v.first.empty());
    trainingProperties.push_back( pt.get_child(v.second.data()));
  }

  std::string inputFilename = pt.get<std::string>("inputFile");

  TFile *inputFile = TFile::Open(inputFilename.c_str());
  TTree *inputTree = (TTree*)(inputFile->Get("PUPPET/t"));


  for(size_t iTrain = 0; iTrain < trainingProperties.size(); ++iTrain)
  {
    applyTraining *user = new applyTraining(trainingProperties[iTrain], inputTree);
    user->getResults();
    delete user;
  }
}

