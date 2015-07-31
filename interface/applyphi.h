#ifndef APPLYTRAINING
#define APPLYTRAINING

#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "GBRTrainer.h"
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



class applyTraining {
  public:
    applyTraining(boost::property_tree::ptree &pt, TTree *inputTree);
    applyTraining(boost::property_tree::ptree &pt, TTree *inputTree, std::string &friendFilename, std::string &friendTreename);
    void registerUpdatedFourVector();
    void registerUpdatedMET();
    void registerUpdatedCovMatrix();
    void calculateUpdatedFourVector();
    void calculateUpdatedMET();
    void calculateUpdatedMETCovMatrix();
    void getResults();
    virtual ~applyTraining();

  protected:
  int _mode = 0;
  void wireInputs();
  virtual void eventLoop();
  TLorentzVector _z;
  float _z_pT                = 0;
  float _z_Phi               = 0;
  float _mvaCov1, _mvaCov2;

  TLorentzVector _oldU;
  float _old_U;
  float _old_UPhi;

  float _mvaResponse         = 0;

  TLorentzVector _newU;
  float _new_U;
  float _new_LongZ;
  float _new_PerpZ;
  float _new_UPhi;

  TLorentzVector _MET;
  float _new_met;
  float _new_metphi;

  float _Cov11, _Cov12, _Cov21, _Cov22;

  std::string _iTrain;
  std::string _iName;
  std::string _applyMVAto;
  std::string _mvaResponseName;
  TFile *_lFForest;
  const std::vector<std::string> *_lVars;
  int _lN = 0;
  TTreeFormula **_lFVars;
  Float_t *_lVals;

  TFile *_lInput;
  TTree *_lTree;

  int _lNEvents = 0;

  std::string _outputFilename;
  TFile *_lOFile;
  TTree *_lOTree;

  double _xResult;
  double _yResult;
};

class applyTraining1D : public applyTraining {
  public:
    applyTraining1D(boost::property_tree::ptree &pt, TTree *inputTree, std::string &friendFilename, std::string &friendTreename) : 
      applyTraining::applyTraining(pt, inputTree, friendFilename, friendTreename),
        _lForest((_mode>0)  ? (GBRForest*)_lFForest->Get(_mvaResponseName.c_str()) : NULL)
{
};


  private:
    const GBRForest * _lForest;
    virtual void eventLoop();
};

class applyTraining2D : public applyTraining {
  public:
    applyTraining2D(boost::property_tree::ptree &pt, TTree *inputTree, std::string &friendFilename, std::string &friendTreename):
      applyTraining::applyTraining(pt, inputTree, friendFilename, friendTreename),
        _lForest( (GBRForest2D*)_lFForest->Get(_mvaResponseName.c_str())) {};

  private:
    GBRForest2D * _lForest;
    virtual void eventLoop();
};

#endif
