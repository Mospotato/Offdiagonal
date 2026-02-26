#ifndef ReadUrQMD_h
#define ReadUrQMD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include "TH1.h"
#include "TH2.h"
#include "TProfile3D.h"
#include "TLorentzVector.h"
// Header file for the classes stored in the TTree if any.
#define Max_Num 9999
class ReadUrQMD : public TSelector
{
public:
   TTree *fChain;
   Int_t mul;
   Float_t b;
   Int_t Npart;
   Int_t pid[Max_Num];
   Float_t px[Max_Num];
   Float_t py[Max_Num];
   Float_t pz[Max_Num];

   // List of branches
   TBranch *b_mul;
   TBranch *b_b;
   TBranch *b_Npart;
   TBranch *b_pid;
   TBranch *b_px; 
   TBranch *b_py; 
   TBranch *b_pz; 

   ReadUrQMD(TTree * /*tree*/ = 0) : fChain(0) {FXTFlag = false;}
   ReadUrQMD(Int_t Energy, Bool_t flag, std::string);
   virtual ~ReadUrQMD() {}
   virtual Int_t Version() const { return 2; }
   virtual void Begin(TTree *tree);
   virtual void SlaveBegin(TTree *tree);
   virtual void Init(TTree *tree);
   virtual Bool_t Notify();
   virtual Bool_t Process(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void SetOption(const char *option) { fOption = option; }
   virtual void SetObject(TObject *obj) { fObject = obj; }
   virtual void SetInputList(TList *input) { fInput = input; }
   virtual TList *GetOutputList() const { return fOutput; }
   virtual void SlaveTerminate();
   virtual void Terminate();
   void GetMoment(Int_t EPDMult, const std::vector<std::map<Int_t, Int_t>> &ParticleMap);
   void getMoment(const std::map<Int_t, Int_t> &ParticleMap, std::vector<Double_t>&);
   Bool_t IsEtaValid(TLorentzVector &vector);
   Double_t GetRapidity(TLorentzVector &vector);
private:
   Bool_t FXTFlag;
   TVector3 BoostFactor;
   Double_t BeamRapidity;
   std::string JobId;
   TFile *OutProfile;
   TH1D *hNpart;
   TH1D *hRefMult;
   TH2D *hMultNpart;
   TH2D *hMultImpact;
   std::vector<std::vector<std::vector<Double_t>>> ConstituentMap;
   ClassDef(ReadUrQMD, 0);
};

#endif

#ifdef ReadUrQMD_cxx
void ReadUrQMD::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mul", &mul, &b_mul);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("Npart", &Npart, &b_Npart);
   fChain->SetBranchAddress("pid", pid, &b_pid);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
}

Bool_t ReadUrQMD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ReadTree_cxx
