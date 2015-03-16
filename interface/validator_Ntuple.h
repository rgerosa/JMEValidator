//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 16 15:35:42 2014 by ROOT version 5.32/00
// from TTree t/t
//////////////////////////////////////////////////////////

#ifndef validatorNtuple_h
#define validatorNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

using std::vector;

class validatorNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<int>*    npus;
   vector<float>*  tnpus;
   vector<int>*    bxns;
   Float_t         rho;
   Float_t         beta;
   Float_t         betaStar;
   Long64_t        npv;
   Long64_t        run;
   Long64_t        lumi;
   Long64_t        evt;
   UInt_t          nref;
   vector<int>*    refrank;   //[nref]
   vector<int>*    refpdgid;   //[nref]
   vector<int>*    refpdgid_algorithmicDef;   //[nref]
   vector<int>*    refpdgid_physicsDef;   //[nref]
   vector<float>*  refe;   //[nref]
   vector<float>*  refpt;   //[nref]
   vector<float>*  refeta;   //[nref]
   vector<float>*  refphi;   //[nref]
   vector<float>*  refy;   //[nref]
   vector<float>*  refdrjt;   //[nref]
   vector<float>*  refarea;   //[nref]
   vector<float>*  jte;   //[nref]
   vector<float>*  jtpt;   //[nref]
   vector<float>*  jteta;   //[nref]
   vector<float>*  jtphi;   //[nref]
   vector<float>*  jty;   //[nref]
   vector<float>*  jtjec;   //[nref]
   vector<float>*  jtarea;   //[nref]
   vector<float>*  jtchf;   //[nref]
   vector<float>*  jtnhf;   //[nref]
   vector<float>*  jtnef;   //[nref]
   vector<float>*  jtcef;   //[nref]
   vector<float>*  jtmuf;   //[nref]
   vector<float>*  jthfhf;   //[nref]
   vector<float>*  jthfef;   //[nref]
   UChar_t         nmu;   //[nmu]
   Float_t         mupt[92];   //[nmu]
   Float_t         mueta[92];   //[nmu]
   Float_t         muphi[92];   //[nmu]
   Float_t         mue[92];   //[nmu]
   Float_t         muIsoRAW[92];   //[nmu]
   Float_t         muIsoSTAND[92];   //[nmu]
   Float_t         muIsoPFWGT[92];   //[nmu]
   Float_t         muIsoPUPPI[92];   //[nmu]

   // List of branches
   TBranch        *b_npus;   //!
   TBranch        *b_tnpus;   //!
   TBranch        *b_bxns;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_betaStar;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_refrank;   //!
   TBranch        *b_refpdgid;   //!
   TBranch        *b_refpdgid_algorithmicDef;   //!
   TBranch        *b_refpdgid_physicsDef;   //!
   TBranch        *b_refe;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refarea;   //!
   TBranch        *b_jte;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtjec;   //!
   TBranch        *b_jtarea;   //!
   TBranch        *b_jtchf;   //!
   TBranch        *b_jtnhf;   //!
   TBranch        *b_jtnef;   //!
   TBranch        *b_jtcef;   //!
   TBranch        *b_jtmuf;   //!
   TBranch        *b_jthfhf;   //!
   TBranch        *b_jthfef;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_mupt;   //!
   TBranch        *b_mueta;   //!
   TBranch        *b_muphi;   //!
   TBranch        *b_mue;   //!
   TBranch        *b_muIsoRAW;   //!
   TBranch        *b_muIsoSTAND;   //!
   TBranch        *b_muIsoPFWGT;   //!
   TBranch        *b_muIsoPUPPI;   //!

   validatorNtuple(TTree *tree=0, bool newTree = false);
   virtual ~validatorNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     MakeTree(TTree *tree);

   int itIndex();
   double sumEOOT();
   double sumLOOT();
};

#endif

#ifdef validatorNtuple_cxx
validatorNtuple::validatorNtuple(TTree *tree, bool newTree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if(newTree) {
      MakeTree(tree);
   }
   else {
      Init(tree);
   }
}

validatorNtuple::~validatorNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t validatorNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t validatorNtuple::LoadTree(Long64_t entry)
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

void validatorNtuple::MakeTree(TTree *tree)
{
   // Set object pointer
   npus         = new vector<int>;
   tnpus        = new vector<float>;
   bxns         = new vector<int>;
   refrank      = new vector<int>;
   refpdgid     = new vector<int>;
   refpdgid_algorithmicDef = new vector<int>;
   refpdgid_physicsDef = new vector<int>;
   refe         = new vector<float>;
   refpt        = new vector<float>;
   refeta       = new vector<float>;
   refphi       = new vector<float>;
   refy         = new vector<float>;
   refdrjt      = new vector<float>;
   refarea      = new vector<float>;
   jte          = new vector<float>;
   jtpt         = new vector<float>;
   jteta        = new vector<float>;
   jtphi        = new vector<float>;
   jty          = new vector<float>;
   jtjec        = new vector<float>;
   jtarea       = new vector<float>;
   jtchf        = new vector<float>;
   jtnhf        = new vector<float>;
   jtnef        = new vector<float>;
   jtcef        = new vector<float>;
   jtmuf        = new vector<float>;
   jthfhf       = new vector<float>;
   jthfef       = new vector<float>;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->Branch("npus", "vector<Int_t>", &npus);
   fChain->Branch("tnpus", "vector<Float_t>", &tnpus);
   fChain->Branch("bxns", "vector<Int_t>", &bxns);
   fChain->Branch("rho", &rho, "rho/F");
   fChain->Branch("beta", &beta, "beta/F");
   fChain->Branch("betaStar", &betaStar, "betaStar/F");
   fChain->Branch("npv", &npv, "npv/L");
   fChain->Branch("run", &run, "run/L");
   fChain->Branch("lumi", &lumi, "lumi/L");
   fChain->Branch("evt", &evt, "evt/L");
   
   fChain->Branch("nref", &nref, "nref/I");
   fChain->Branch("refrank", "vector<Int_t>", &refrank);
   fChain->Branch("refpdgid", "vector<Int_t>", &refpdgid);
   fChain->Branch("refpdgid_algorithmicDef", "vector<Int_t>", &refpdgid_algorithmicDef);
   fChain->Branch("refpdgid_physicsDef", "vector<Int_t>", &refpdgid_physicsDef);
   fChain->Branch("refe", "vector<Float_t>", &refe);      
   fChain->Branch("refpt", "vector<Float_t>", &refpt);    
   fChain->Branch("refeta", "vector<Float_t>", &refeta);  
   fChain->Branch("refphi", "vector<Float_t>", &refphi);  
   fChain->Branch("refy", "vector<Float_t>", &refy);      
   fChain->Branch("refdrjt", "vector<Float_t>", &refdrjt);
   fChain->Branch("refarea", "vector<Float_t>", &refarea);
   fChain->Branch("jte", "vector<Float_t>", &jte);        
   fChain->Branch("jtpt", "vector<Float_t>", &jtpt);      
   fChain->Branch("jteta", "vector<Float_t>", &jteta);    
   fChain->Branch("jtphi", "vector<Float_t>", &jtphi);    
   fChain->Branch("jty", "vector<Float_t>", &jty);        
   fChain->Branch("jtjec", "vector<Float_t>", &jtjec);
   fChain->Branch("jtarea", "vector<Float_t>", &jtarea);
   fChain->Branch("jtchf", "vector<Float_t>", &jtchf);
   fChain->Branch("jtnhf", "vector<Float_t>", &jtnhf);
   fChain->Branch("jtnef", "vector<Float_t>", &jtnef);
   fChain->Branch("jtcef", "vector<Float_t>", &jtcef);
   fChain->Branch("jtmuf", "vector<Float_t>", &jtmuf);
   fChain->Branch("jthfhf", "vector<Float_t>", &jthfhf);
   fChain->Branch("jthfef", "vector<Float_t>", &jthfef);

   fChain->Branch("nmu", &nmu, "nmu/b");
   // fChain->Branch("mupt", mupt, "mupt[nmu]/F");
   // fChain->Branch("mueta", mueta, "mueta[nmu]/F");
   // fChain->Branch("muphi", muphi, "muphi[nmu]/F");
   // fChain->Branch("mue", mue, "mue[nmu]/F");
   // fChain->Branch("muIsoRAW", muIsoRAW, "muIsoRAW[nmu]/F");
   // fChain->Branch("muIsoSTAND", muIsoSTAND, "muIsoSTAND[nmu]/F");
   // fChain->Branch("muIsoPFWGT", muIsoPFWGT, "muIsoPFWGT[nmu]/F");
   // fChain->Branch("muIsoPUPPI", muIsoPUPPI, "muIsoPUPPI[nmu]/F");
   Notify();
}

void validatorNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   npus         = new vector<int>;
   tnpus        = new vector<float>;
   bxns         = new vector<int>;
   refrank      = new vector<int>;
   refpdgid     = new vector<int>;
   refpdgid_algorithmicDef = new vector<int>;
   refpdgid_physicsDef = new vector<int>;
   refe         = new vector<float>;
   refpt        = new vector<float>;
   refeta       = new vector<float>;
   refphi       = new vector<float>;
   refy         = new vector<float>;
   refdrjt      = new vector<float>;
   refarea      = new vector<float>;
   jte          = new vector<float>;
   jtpt         = new vector<float>;
   jteta        = new vector<float>;
   jtphi        = new vector<float>;
   jty          = new vector<float>;
   jtjec        = new vector<float>;
   jtarea       = new vector<float>;
   jtchf        = new vector<float>;
   jtnhf        = new vector<float>;
   jtnef        = new vector<float>;
   jtcef        = new vector<float>;
   jtmuf        = new vector<float>;
   jthfhf       = new vector<float>;
   jthfef       = new vector<float>;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npus", &npus, &b_npus);
   fChain->SetBranchAddress("tnpus", &tnpus, &b_tnpus);
   fChain->SetBranchAddress("bxns", &bxns, &b_bxns);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("beta", &beta, &b_beta);
   fChain->SetBranchAddress("betaStar", &betaStar, &b_betaStar);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);

   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("refrank", &refrank, &b_refrank);
   fChain->SetBranchAddress("refpdgid", &refpdgid, &b_refpdgid);
   fChain->SetBranchAddress("refpdgid_algorithmicDef", &refpdgid_algorithmicDef, &b_refpdgid_algorithmicDef);
   fChain->SetBranchAddress("refpdgid_physicsDef", &refpdgid_physicsDef, &b_refpdgid_physicsDef);
   fChain->SetBranchAddress("refe", &refe, &b_refe);
   fChain->SetBranchAddress("refpt", &refpt, &b_refpt);
   fChain->SetBranchAddress("refeta", &refeta, &b_refeta);
   fChain->SetBranchAddress("refphi", &refphi, &b_refphi);
   fChain->SetBranchAddress("refy", &refy, &b_refy);
   fChain->SetBranchAddress("refdrjt", &refdrjt, &b_refdrjt);
   fChain->SetBranchAddress("refarea", &refarea, &b_refarea);
   fChain->SetBranchAddress("jte", &jte, &b_jte);
   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("jty", &jty, &b_jty);
   fChain->SetBranchAddress("jtjec", &jtjec, &b_jtjec);
   fChain->SetBranchAddress("jtarea", &jtarea, &b_jtarea);
   fChain->SetBranchAddress("jtchf", &jtchf, &b_jtchf);
   fChain->SetBranchAddress("jtnhf", &jtnhf, &b_jtnhf);
   fChain->SetBranchAddress("jtnef", &jtnef, &b_jtnef);
   fChain->SetBranchAddress("jtcef", &jtcef, &b_jtcef);
   fChain->SetBranchAddress("jtmuf", &jtmuf, &b_jtmuf);
   fChain->SetBranchAddress("jthfhf", &jthfhf, &b_jthfhf);
   fChain->SetBranchAddress("jthfef", &jthfef, &b_jthfef);
   fChain->SetBranchAddress("nmu", &nmu, &b_nmu);
   // fChain->SetBranchAddress("mupt", mupt, &b_mupt);
   // fChain->SetBranchAddress("mueta", mueta, &b_mueta);
   // fChain->SetBranchAddress("muphi", muphi, &b_muphi);
   // fChain->SetBranchAddress("mue", mue, &b_mue);
   // fChain->SetBranchAddress("muIsoRAW", muIsoRAW, &b_muIsoRAW);
   // fChain->SetBranchAddress("muIsoSTAND", muIsoSTAND, &b_muIsoSTAND);
   // fChain->SetBranchAddress("muIsoPFWGT", muIsoPFWGT, &b_muIsoPFWGT);
   // fChain->SetBranchAddress("muIsoPUPPI", muIsoPUPPI, &b_muIsoPUPPI);
   Notify();
}

Bool_t validatorNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void validatorNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t validatorNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef validatorNtuple_cxx
