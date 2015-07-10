#include "TFile.h"
#include "TTree.h"
#include "../interface/GBRTrainer.h"
#include "../interface/GBRTrainer2D.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "CondFormats/EgammaObjects/interface/GBRForest2D.h"

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

  // Reading in from config file 
  std::vector<std::string> *friends = new std::vector<std::string>;
  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("friends"))
  {
    assert(v.first.empty());
    friends->push_back( v.second.data() );
  }
  for (int i=0; i<int(friends->size()); ++i) {
    lRegress->AddFriend(friends->at(i).c_str(), (friends->at(i)+".root").c_str());
  }
  std::vector<std::string> *targetVariables = new std::vector<std::string>;
  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("targetVariables"))
  {
    assert(v.first.empty());
    targetVariables->push_back( v.second.data() );
  }

  std::vector<std::string> *lVec = new std::vector<std::string>;
  BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("trainingVariables"))
  {
    assert(v.first.empty());
    lVec->push_back( v.second.data() );
  }

  std::string weight = pt.get<std::string>("weight");

  TFile *fout = new TFile(pt.get<std::string>("weightfilename").c_str(),"RECREATE");

  std::string mvaResponseName = pt.get<std::string>("name");

  if( targetVariables->size() == 1)
  {
    // 1-dim training
    GBRTrainer *train = new GBRTrainer;
    train->AddTree(lRegress);
    train->SetTrainingCut(weight);
    train->SetMinEvents( pt.get<int>("minEvents") );
    train->SetShrinkage( pt.get<float>("shrinkage") );
    train->SetMinCutSignificance( pt.get<float>("minCutSignificance") );  
    cout << " ===> " << weight << endl;  
    train->SetTargetVar( targetVariables->at(0) );
    for (int i=0; i<int(lVec->size()); ++i) {
     train->AddInputVar(lVec->at(i));
    }
    std::cout << "Training Forest with " << pt.get<int>("nTrees")<< " Trees" << std::endl; 
    std::cout << pt.get<std::string>("desc") << std::endl;
    const GBRForest *forest = train->TrainForest( pt.get<int>("nTrees"));
    fout->WriteObject(forest, mvaResponseName.c_str());
  }
  else if (targetVariables->size() == 2)
  {
   // 2-dim training
    GBRTrainer2D *train = new GBRTrainer2D;
    train->SetTree(lRegress);
    train->SetTargetXVar(targetVariables->at(0));
    train->SetTargetYVar(targetVariables->at(1));
    train->SetTrainingCut(weight);
    train->SetMinEvents( pt.get<int>("minEvents") );
    train->SetShrinkage( pt.get<float>("shrinkage") );
    for (int i=0; i<int(lVec->size()); ++i) {
     train->AddInputVar(lVec->at(i));
    }
    const GBRForest2D *forest = train->TrainForest( pt.get<int>("nTrees"));
    fout->WriteObject(forest, mvaResponseName.c_str());
  } 

  fout->WriteObject(lVec, "varlist");
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
  std::string Treename = pt.get<std::string>("Folder");
  TTree *inputTree = (TTree*)(inputFile->Get(Treename.c_str()));

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

    int nTargetVars = 0;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, trainingProperties[iTrain].get_child("targetVariables"))
    {
      nTargetVars++;
    }
    if(nTargetVars == 1)
    {
      applyTraining1D user = applyTraining1D(trainingProperties[iTrain], inputTree, friendFilename, friendTreename);
      user.getResults();
    }
    if(nTargetVars == 2)
    {
      applyTraining2D user = applyTraining2D(trainingProperties[iTrain], inputTree, friendFilename, friendTreename);
      user.getResults();
    }
    inputTree->AddFriend(friendTreename.c_str(), friendFilename.c_str());
  }
}
