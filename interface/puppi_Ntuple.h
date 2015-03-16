//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 13 19:38:57 2015 by ROOT version 5.34/18
// from TTree puppiTree/puppiTree
// found on file: ../test/test.root
//////////////////////////////////////////////////////////

#ifndef puppiNtuple_h
#define puppiNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

using std::vector;

class puppiNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Long64_t        run;
   Long64_t        lumi;
   Long64_t        evt;
   Float_t         nalgos;
   vector<float>   *px;
   vector<float>   *py;
   vector<float>   *pz;
   vector<float>   *e;
   vector<float>   *alphas;
   vector<float>   *id;
   vector<float>   *charge;
   vector<float>   *fromPV;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_nalgos;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_e;   //!
   TBranch        *b_alphas;   //!
   TBranch        *b_id;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_fromPV;   //!

   puppiNtuple(TTree *tree=0, bool newTree = false);
   virtual ~puppiNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     MakeTree(TTree *tree);
};

#endif

#ifdef puppiNtuple_cxx
puppiNtuple::puppiNtuple(TTree *tree, bool newTree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (newTree) {
      MakeTree(tree);
   }
   else {
      Init(tree);
   }
}

puppiNtuple::~puppiNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t puppiNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t puppiNtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void puppiNtuple::MakeTree(TTree *tree)
{
   // Set object pointer
   px     = new vector<float>;
   py     = new vector<float>;
   pz     = new vector<float>;
   e      = new vector<float>;
   alphas = new vector<float>;
   id     = new vector<float>;
   charge = new vector<float>;
   fromPV = new vector<float>;   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->Branch("run", &run, "run/L");
   fChain->Branch("lumi", &lumi, "lumi/L");
   fChain->Branch("evt", &evt, "evt/L");
   fChain->Branch("nalgos", &nalgos, "nalgos/F");
   fChain->Branch("px",     "vector<float>", &px);
   fChain->Branch("py",     "vector<float>", &py);
   fChain->Branch("pz",     "vector<float>", &pz);
   fChain->Branch("e",      "vector<float>", &e);
   fChain->Branch("alphas", "vector<float>", &alphas);
   fChain->Branch("id",     "vector<float>", &id);
   fChain->Branch("charge", "vector<float>", &charge);
   fChain->Branch("fromPV", "vector<float>", &fromPV);
   Notify();
}

void puppiNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   px = 0;
   py = 0;
   pz = 0;
   e = 0;
   alphas = 0;
   id = 0;
   charge = 0;
   fromPV = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("nalgos", &nalgos, &b_nalgos);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("e", &e, &b_e);
   fChain->SetBranchAddress("alphas", &alphas, &b_alphas);
   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("fromPV", &fromPV, &b_fromPV);
   Notify();
}

Bool_t puppiNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void puppiNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t puppiNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef puppiNtuple_cxx
