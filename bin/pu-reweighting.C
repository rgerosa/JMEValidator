#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TTreeFormula.h"
#include <iostream>
#include <regex>
//#include "Cintex/Cintex.h"
//
using namespace std;




int main(int argc, char* argv[] ) {

  
  TFile *histoFile = TFile::Open((argv[1]));
  TH1D *weightHisto = (TH1D*) histoFile->Get("Data");
  TH1D *MCHisto = (TH1D*) histoFile->Get("MC");

//  std::cout << dataHisto->GetNbinsX() << "/" << MCHisto->GetNbinsX() << std::endl;
  weightHisto->Divide(weightHisto, MCHisto, 1, 1);

  for(int i = 0; i < weightHisto->GetNbinsX(); i++)
  {
    std::cout << "bin " << i << " : " << weightHisto->GetBinContent(i) << std::endl;
  }

  TFile *lInput = TFile::Open((argv[2]));
  TTree *lTree = (TTree*) lInput->Get("PUPPET/t");
// npv.root open
// open Histogram data
// print contents
// loop over tree, add variable puWeight
  std::string inputstring1(argv[2]);
  std::string outputfilename(inputstring1.substr(0, inputstring1.find_last_of(".")) + "_pu_weight.root");

  TFile *lOFile = new TFile(outputfilename.c_str(), "RECREATE");
  lOFile->mkdir("PUPPET");
  lOFile->cd("PUPPET");
  TTree *lOTree = lTree->CloneTree(0);
  int lNEvents = lTree->GetEntries();

  float puWeight;
  lOTree->Branch("puWeight", &puWeight, "puWeight/F");
  TTreeFormula *npv = new TTreeFormula("npv", "NVertex", lTree);
  for (Long64_t i0=0; i0<lNEvents;i0++)
  {
    lTree->GetEntry(i0);
    puWeight = weightHisto->GetBinContent(int(npv->EvalInstance())+1);
 //  std::cout << "Event: " << i0 << ", npv: " << npv->EvalInstance() << ", int: " << int(npv->EvalInstance())<< ", weight: "  <<puWeight<< ", check: " << weightHisto->GetBinLowEdge(int(npv->EvalInstance())) << std::endl; 
    if(i0%1000000==0)
        std::cout << "Event " << i0 << "/" << lNEvents << std::endl;
        lOTree->Fill();
  }

  lOTree->Write();
  lOFile->Close();
  std::cout << "Writing output file " << outputfilename << ".... Done!" << std::endl;
  return 0;
}
