#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <iostream>
#include <regex>
//#include "Cintex/Cintex.h"
//
using namespace std;




int main(int argc, char* argv[] ) {

//read in
//evaluate cut
//write out
  std::cout << "Selecting from file " << argv[1]<< " only events with trigger '" << argv[2] << "'" << std::endl;
  TFile *lInput = TFile::Open((argv[1]));
  TTree *lTree = (TTree*) lInput->Get("PUPPET/t");
  std::string inputstring1(argv[1]);
  std::string inputstring2(argv[2]);
  std::string outputfilename(inputstring1.substr(0, inputstring1.find_last_of(".")) + "_" + inputstring2 + ".root");
  std::regex expression(argv[2]); //"(nLeptons==2)*(lepton1_charge!=0)*(lepton1_charge==-lepton2_charge)*(z_m>60)*(z_m<120)";
  //TTreeFormula *Cut = new TTreeFormula("name", expression, lTree);

  TFile *lOFile = new TFile(outputfilename.c_str(), "RECREATE");
  TTree *lOTree = lTree->CloneTree(0);
  int lNEvents = lTree->GetEntries();
  std::vector<std::string> *triggers= 0;
  lTree->SetBranchAddress("DoubleMuPaths", &triggers);
  for (Long64_t i0=0; i0<lNEvents;i0++)
  {
    if(i0%1000000==0)
        std::cout << "Event " << i0 << "/" << lNEvents << std::endl;
    lTree->GetEntry(i0);
    std::cout << "Trigger: " << triggers->size() << std::endl;
    for(size_t i = 0; i < triggers->size(); i++)
    {
    std::cout << "i: " << i << ": " << triggers->at(i) << std::endl;
      //if(triggers->at(i).find(expression)!=std::string::npos)
      if(std::regex_match(triggers->at(i), expression))
      {
        std::cout << "match" << std::endl;
        lOTree->Fill();
        //break;
      }
    }
    // std::cout << triggers->at(i) << ", ";
    //std::cout << std::endl;
    //float result = (float)Cut->EvalInstance();
  }

  lOTree->Write();
  lOFile->Close();
  std::cout << "Writing output file " << outputfilename << ".... Done!" << std::endl;
  return 0;
}
