#include "TFile.h"
#include "TTree.h"
#include "../interface/GBRTrainer.h"
#include "../../../CondFormats/EgammaObjects/interface/GBRForest.h"
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
#include "../interface/flattenDistribution.h"

using namespace std;

void doTraining(boost::property_tree::ptree &pt, TTree* lRegress)
{

  std::vector<std::string> *friends = new std::vector<std::string>;
  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("friends"))
  {
    assert(v.first.empty());
    friends->push_back( v.second.data() );
  }
  for (int i=0; i<int(friends->size()); ++i) {
    lRegress->AddFriend(friends->at(i).c_str(), (friends->at(i)+".root").c_str());
  }
  std::string weight = pt.get<std::string>("weight");
  GBRTrainer *train = new GBRTrainer;
  train->AddTree(lRegress);
  train->SetTrainingCut(weight);
  train->SetMinEvents( pt.get<int>("minEvents") );
  train->SetShrinkage( pt.get<float>("shrinkage") );
  train->SetMinCutSignificance( pt.get<float>("minCutSignificance") );  
  cout << " ===> " << weight << endl;  
  train->SetTargetVar( pt.get<std::string>("targetVariable") );

  std::vector<std::string> *lVec = new std::vector<std::string>;
  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("trainingVariables"))
  {
    assert(v.first.empty());
    lVec->push_back( v.second.data() );
  }

  for (int i=0; i<int(lVec->size()); ++i) {
     train->AddInputVar(lVec->at(i));
  }
  
  std::cout << "Training Forest with " << pt.get<int>("nTrees")<< " Trees" << std::endl; 
  std::cout << pt.get<std::string>("desc") << std::endl;
  const GBRForest *forest = train->TrainForest( pt.get<int>("nTrees"));
  //ROOT::Cintex::Cintex::Enable();
  TFile *fout = new TFile(pt.get<std::string>("weightfilename").c_str(),"RECREATE");    
  std::string mvaResponseName = pt.get<std::string>("name");
  fout->WriteObject(lVec, "varlist");
  fout->WriteObject(forest, mvaResponseName.c_str());
  fout->Close();
}


int main(int argc, char* argv[] ) {
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

  std::vector<boost::property_tree::ptree> reweightingProperties;
  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("activeReweightings"))
  {
    assert(v.first.empty());
    reweightingProperties.push_back( pt.get_child(v.second.data()));
  }

  std::string inputFilename = pt.get<std::string>("inputFile");

  TFile *inputFile = TFile::Open(inputFilename.c_str());
  TTree *inputTree = (TTree*)(inputFile->Get("Flat"));


  for(size_t iTrain = 0; iTrain < reweightingProperties.size(); ++iTrain)
  {
    std::string friendFilename;
    std::string friendTreename;
    distributionFlatter *flatter = new distributionFlatter(reweightingProperties[iTrain], inputTree);
    flatter->calculateWeights();
    flatter->writeWeightToTree(friendFilename, friendTreename); 
    inputTree->AddFriend(friendTreename.c_str(), friendFilename.c_str());
  }


  for(size_t iTrain = 0; iTrain < trainingProperties.size(); ++iTrain)
  {
    int mode = trainingProperties[iTrain].get<int>("mode");
    if( mode > 0)
    {
      doTraining(trainingProperties[iTrain], inputTree);
    }

    std::string friendFilename;
    std::string friendTreename;
    {
    applyTraining user = applyTraining(trainingProperties[iTrain], inputTree, friendFilename, friendTreename);
    user.getResults();
    }
    inputTree->AddFriend(friendTreename.c_str(), friendFilename.c_str());
  }
  //merge inputTree and its friends to one file
  //std::string outputFilename = "output.root";
  //TFile *outputFile = TFile::Open(outputFilename.c_str(), "RECREATE");
  //TTree *outputTree = inputTree->CloneTree();
  //outputTree->Write();
  //outputFile->Close();
}
