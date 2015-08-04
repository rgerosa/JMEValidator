#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TTreeFormula.h"
#include <iostream>
#include <regex>
#include "TVector2.h"
//#include "Cintex/Cintex.h"
//
using namespace std;




int main(int argc, char* argv[] ) {

  // input file öffnen (argv[1])
  // recoilPFPuppiMet Vector füllen
  // Boson_Pt vector füllen
  // PFPuppiMet2_ füllen
  
  TFile *lInput = TFile::Open((argv[1]));
  TTree *lTree = (TTree*) lInput->Get("PUPPET/t");
// npv.root open
// open Histogram data
// print contents
// loop over tree, add variable puWeight
  std::string inputstring1(argv[1]);
  std::string outputfilename(inputstring1.substr(0, inputstring1.find_last_of(".")) + "_puppiMET.root");

  TFile *lOFile = new TFile(outputfilename.c_str(), "RECREATE");
  lOFile->mkdir("PUPPET");
  lOFile->cd("PUPPET");
  TTree *lOTree = lTree->CloneTree(0);
  int lNEvents = lTree->GetEntries();

  float recoilPFPuppiMet_Pt, recoilPFPuppiMet_Phi;
  float recoilPFPuppiMet_uncorrected_Pt, recoilPFPuppiMet_uncorrected_Phi;
  float Boson_Pt, Boson_Phi;

  lTree->SetBranchAddress("recoilPFPuppiMet_Pt", &recoilPFPuppiMet_Pt);
  lTree->SetBranchAddress("recoilPFPuppiMet_Phi", &recoilPFPuppiMet_Phi);
  lTree->SetBranchAddress("recoilPFPuppiMet_uncorrected_Pt", &recoilPFPuppiMet_uncorrected_Pt);
  lTree->SetBranchAddress("recoilPFPuppiMet_uncorrected_Phi", &recoilPFPuppiMet_uncorrected_Phi);
  lTree->SetBranchAddress("Boson_Pt", &Boson_Pt);
  lTree->SetBranchAddress("Boson_Phi", &Boson_Phi);

  float PFPuppiMet2_Pt, PFPuppiMet2_Phi;
  float PFPuppiMet2_uncorrected_Pt, PFPuppiMet2_uncorrected_Phi;

  lOTree->Branch("PFPuppiMet2_Pt",  &PFPuppiMet2_Pt,  "PFPuppiMet2_Pt/F");
  lOTree->Branch("PFPuppiMet2_Phi", &PFPuppiMet2_Phi, "PFPuppiMet2_Phi/F");

  lOTree->Branch("PFPuppiMet2_uncorrected_Pt",  &PFPuppiMet2_uncorrected_Pt,  "PFPuppiMet2_uncorrected_Pt/F");
  lOTree->Branch("PFPuppiMet2_uncorrected_Phi", &PFPuppiMet2_uncorrected_Phi, "PFPuppiMet2_uncorrected_Phi/F");

  TVector2 puppiRecoil, puppiRecoil_uncorrected, Boson, puppiMet, puppiMet_uncorrected;


  for (Long64_t i0=0; i0<lNEvents;i0++)
  {
    lTree->GetEntry(i0);

    puppiRecoil.SetMagPhi(recoilPFPuppiMet_Pt, recoilPFPuppiMet_Phi);
    puppiRecoil_uncorrected.SetMagPhi(recoilPFPuppiMet_uncorrected_Pt, recoilPFPuppiMet_uncorrected_Phi);
    Boson.SetMagPhi(Boson_Pt, Boson_Phi);

    puppiMet = Boson + puppiRecoil;
    puppiMet = puppiMet.Rotate(TMath::Pi());

    puppiMet_uncorrected = Boson + puppiRecoil_uncorrected;
    puppiMet_uncorrected = puppiMet_uncorrected.Rotate(TMath::Pi());

    PFPuppiMet2_Pt  = puppiMet.Mod();
    PFPuppiMet2_Phi = TVector2::Phi_mpi_pi(puppiMet.Phi());

    PFPuppiMet2_uncorrected_Pt  = puppiMet_uncorrected.Mod();
    PFPuppiMet2_uncorrected_Phi = TVector2::Phi_mpi_pi(puppiMet_uncorrected.Phi());
    //std::cout << "puppi recoil" << recoilPFPuppiMet_Pt << "/" << recoilPFPuppiMet_Phi << ", boson: " << Boson_Pt << "/" << Boson_Phi << ", met: " << PFPuppiMet2_Pt << ", metphi: " << PFPuppiMet2_Phi << std::endl;
    if(i0%1000000==0)
        std::cout << "Event " << i0 << "/" << lNEvents << std::endl;
    lOTree->Fill();
  }

  lOTree->Write();
  lOFile->Close();
  std::cout << "Writing output file " << outputfilename << ".... Done!" << std::endl;
  return 0;
}
