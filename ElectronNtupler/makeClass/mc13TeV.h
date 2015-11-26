//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep  9 07:51:59 2015 by ROOT version 6.02/05
// from TTree ElectronTree/Electron data
// found on file: ../DYJetsToLL_M50.root
//////////////////////////////////////////////////////////

#ifndef mc13TeV_h
#define mc13TeV_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class mc13TeV {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        RunNo;
   Double_t        EvtNo;
   Double_t        Lumi;
   Double_t        Bunch;
   Int_t           weightsP;
   Int_t           weightsN;
   Double_t        theWeight;
   Int_t           pvNTracks;
   Int_t           nPV;
   Int_t           nPU;
   Int_t           nPUTrue;
   Float_t         rho;
   Bool_t          singleElectron_mc;
   Bool_t          singleElectron_data;
   Bool_t          doubleElectron;
   vector<double>  *pt_leg1;
   vector<double>  *eta_leg1;
   vector<double>  *phi_leg1;
   vector<double>  *pt_leg2;
   vector<double>  *eta_leg2;
   vector<double>  *phi_leg2;
   vector<double>  *pt_SE_MC;
   vector<double>  *eta_SE_MC;
   vector<double>  *phi_SE_MC;
   vector<double>  *pt_SE_Data;
   vector<double>  *eta_SE_Data;
   vector<double>  *phi_SE_Data;
   Int_t           nEle;
   Int_t           nGenEle;
   vector<float>   *ZMass;
   vector<float>   *ZPt;
   vector<float>   *ZEta;
   vector<float>   *ZRap;
   vector<float>   *ZPhi;
   vector<float>   *gPre_energy;
   vector<float>   *gPre_px;
   vector<float>   *gPre_py;
   vector<float>   *gPre_pz;
   vector<float>   *gPre_pt;
   vector<float>   *gPre_eta;
   vector<float>   *gPre_rap;
   vector<float>   *gPre_phi;
   vector<float>   *gPost_energy;
   vector<float>   *gPost_px;
   vector<float>   *gPost_py;
   vector<float>   *gPost_pz;
   vector<float>   *gPost_pt;
   vector<float>   *gPost_eta;
   vector<float>   *gPost_rap;
   vector<float>   *gPost_phi;
   vector<float>   *pt;
   vector<float>   *eta;
   vector<float>   *rap;
   vector<float>   *phi;
   vector<float>   *energy;
   vector<float>   *mass;
   vector<float>   *charge;
   vector<float>   *enSC;
   vector<float>   *preEnSC;
   vector<float>   *rawEnSC;
   vector<float>   *etSC;
   vector<float>   *etaSC;
   vector<float>   *phiSC;
   vector<float>   *full5x5_sigmaIetaIeta;
   vector<float>   *E1x5;
   vector<float>   *E2x5;
   vector<float>   *E5x5;
   vector<float>   *hOverE;
   vector<float>   *etaScWidth;
   vector<float>   *phiScWidth;
   vector<float>   *r9;
   vector<float>   *dEtaIn;
   vector<float>   *dPhiIn;
   vector<float>   *isoChargedHadrons;
   vector<float>   *isoNeutralHadrons;
   vector<float>   *isoPhotons;
   vector<float>   *isoChargedFromPU;
   vector<float>   *isoDeltaBeta;
   vector<float>   *isoRho;
   vector<float>   *ooEmooP;
   vector<float>   *d0;
   vector<float>   *dz;
   vector<int>     *expectedMissingInnerHits;
   vector<int>     *passConversionVeto;
   vector<float>   *brem;
   vector<int>     *isTrue;
   vector<int>     *passVetoId;
   vector<int>     *passLooseId;
   vector<int>     *passMediumId;
   vector<int>     *passTightId;
   vector<float>   *mvaVal;
   vector<int>     *mvaCat;
   vector<int>     *eleEcalDrivenSeed;

   // List of branches
   TBranch        *b_RunNo;   //!
   TBranch        *b_EvtNo;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_Bunch;   //!
   TBranch        *b_weighstP;   //!
   TBranch        *b_weightsN;   //!
   TBranch        *b_theWeight;   //!
   TBranch        *b_pvNTracks;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nPUTrue;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_singleElectron_mc;   //!
   TBranch        *b_singleElectron_data;   //!
   TBranch        *b_doubleElectron;   //!
   TBranch        *b_pt_leg1;   //!
   TBranch        *b_eta_leg1;   //!
   TBranch        *b_phi_leg1;   //!
   TBranch        *b_pt_leg2;   //!
   TBranch        *b_eta_leg2;   //!
   TBranch        *b_phi_leg2;   //!
   TBranch        *b_pt_SE_MC;   //!
   TBranch        *b_eta_SE_MC;   //!
   TBranch        *b_phi_SE_MC;   //!
   TBranch        *b_pt_SE_Data;   //!
   TBranch        *b_eta_SE_Data;   //!
   TBranch        *b_phi_SE_Data;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_nGenEle;   //!
   TBranch        *b_ZMass;   //!
   TBranch        *b_ZPt;   //!
   TBranch        *b_ZEta;   //!
   TBranch        *b_ZRap;   //!
   TBranch        *b_ZPhi;   //!
   TBranch        *b_gPre_energy;   //!
   TBranch        *b_gPre_px;   //!
   TBranch        *b_gPre_py;   //!
   TBranch        *b_gPre_pz;   //!
   TBranch        *b_gPre_pt;   //!
   TBranch        *b_gPre_eta;   //!
   TBranch        *b_gPre_rap;   //!
   TBranch        *b_gPre_phi;   //!
   TBranch        *b_gPost_energy;   //!
   TBranch        *b_gPost_px;   //!
   TBranch        *b_gPost_py;   //!
   TBranch        *b_gPost_pz;   //!
   TBranch        *b_gPost_pt;   //!
   TBranch        *b_gPost_eta;   //!
   TBranch        *b_gPost_rap;   //!
   TBranch        *b_gPost_phi;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_rap;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_enSC;   //!
   TBranch        *b_preEnSC;   //!
   TBranch        *b_rawEnSC;   //!
   TBranch        *b_etSC;   //!
   TBranch        *b_etaSC;   //!
   TBranch        *b_phiSC;   //!
   TBranch        *b_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_E1x5;   //!
   TBranch        *b_E2x5;   //!
   TBranch        *b_E5x5;   //!
   TBranch        *b_hOverE;   //!
   TBranch        *b_etaScWidth;   //!
   TBranch        *b_phiScWidth;   //!
   TBranch        *b_r9;   //!
   TBranch        *b_dEtaIn;   //!
   TBranch        *b_dPhiIn;   //!
   TBranch        *b_isoChargedHadrons;   //!
   TBranch        *b_isoNeutralHadrons;   //!
   TBranch        *b_isoPhotons;   //!
   TBranch        *b_isoChargedFromPU;   //!
   TBranch        *b_isoDeltaBeta;   //!
   TBranch        *b_isoRho;   //!
   TBranch        *b_ooEmooP;   //!
   TBranch        *b_d0;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_expectedMissingInnerHits;   //!
   TBranch        *b_passConversionVeto;   //!
   TBranch        *b_brem;   //!
   TBranch        *b_isTrue;   //!
   TBranch        *b_passVetoId;   //!
   TBranch        *b_passLooseId;   //!
   TBranch        *b_passMediumId;   //!
   TBranch        *b_passTightId;   //!
   TBranch        *b_mvaVal;   //!
   TBranch        *b_mvaCat;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!

   mc13TeV(TTree *tree=0);
   virtual ~mc13TeV();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Double_t deltaPhi(Double_t phientry1 , Double_t phientry2 );
   virtual Double_t deltaEta(Double_t etaentry1 , Double_t etaentry2);
   virtual Double_t deltaR(Double_t etaentry1 , Double_t phientry1 , Double_t etaentry2 , Double_t phientry2);
};

#endif

#ifdef mc13TeV_cxx
mc13TeV::mc13TeV(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DYJetsToLL_M50.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DYJetsToLL_M50.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ntupler;1");
      dir->GetObject("ElectronTree",tree);

   }
   Init(tree);
}

mc13TeV::~mc13TeV()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mc13TeV::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mc13TeV::LoadTree(Long64_t entry)
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

void mc13TeV::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pt_leg1 = 0;
   eta_leg1 = 0;
   phi_leg1 = 0;
   pt_leg2 = 0;
   eta_leg2 = 0;
   phi_leg2 = 0;
   pt_SE_MC = 0;
   eta_SE_MC = 0;
   phi_SE_MC = 0;
   pt_SE_Data = 0;
   eta_SE_Data = 0;
   phi_SE_Data = 0;
   ZMass = 0;
   ZPt = 0;
   ZEta = 0;
   ZRap = 0;
   ZPhi = 0;
   gPre_energy = 0;
   gPre_px = 0;
   gPre_py = 0;
   gPre_pz = 0;
   gPre_pt = 0;
   gPre_eta = 0;
   gPre_rap = 0;
   gPre_phi = 0;
   gPost_energy = 0;
   gPost_px = 0;
   gPost_py = 0;
   gPost_pz = 0;
   gPost_pt = 0;
   gPost_eta = 0;
   gPost_rap = 0;
   gPost_phi = 0;
   pt = 0;
   eta = 0;
   rap = 0;
   phi = 0;
   energy = 0;
   mass = 0;
   charge = 0;
   enSC = 0;
   preEnSC = 0;
   rawEnSC = 0;
   etSC = 0;
   etaSC = 0;
   phiSC = 0;
   full5x5_sigmaIetaIeta = 0;
   E1x5 = 0;
   E2x5 = 0;
   E5x5 = 0;
   hOverE = 0;
   etaScWidth = 0;
   phiScWidth = 0;
   r9 = 0;
   dEtaIn = 0;
   dPhiIn = 0;
   isoChargedHadrons = 0;
   isoNeutralHadrons = 0;
   isoPhotons = 0;
   isoChargedFromPU = 0;
   isoDeltaBeta = 0;
   isoRho = 0;
   ooEmooP = 0;
   d0 = 0;
   dz = 0;
   expectedMissingInnerHits = 0;
   passConversionVeto = 0;
   brem = 0;
   isTrue = 0;
   passVetoId = 0;
   passLooseId = 0;
   passMediumId = 0;
   passTightId = 0;
   mvaVal = 0;
   mvaCat = 0;
   eleEcalDrivenSeed = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
   fChain->SetBranchAddress("EvtNo", &EvtNo, &b_EvtNo);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("Bunch", &Bunch, &b_Bunch);
   fChain->SetBranchAddress("weightsP", &weightsP, &b_weighstP);
   fChain->SetBranchAddress("weightsN", &weightsN, &b_weightsN);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("pvNTracks", &pvNTracks, &b_pvNTracks);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nPUTrue", &nPUTrue, &b_nPUTrue);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("singleElectron_mc", &singleElectron_mc, &b_singleElectron_mc);
   fChain->SetBranchAddress("singleElectron_data", &singleElectron_data, &b_singleElectron_data);
   fChain->SetBranchAddress("doubleElectron", &doubleElectron, &b_doubleElectron);
   fChain->SetBranchAddress("pt_leg1", &pt_leg1, &b_pt_leg1);
   fChain->SetBranchAddress("eta_leg1", &eta_leg1, &b_eta_leg1);
   fChain->SetBranchAddress("phi_leg1", &phi_leg1, &b_phi_leg1);
   fChain->SetBranchAddress("pt_leg2", &pt_leg2, &b_pt_leg2);
   fChain->SetBranchAddress("eta_leg2", &eta_leg2, &b_eta_leg2);
   fChain->SetBranchAddress("phi_leg2", &phi_leg2, &b_phi_leg2);
   fChain->SetBranchAddress("pt_SE_MC", &pt_SE_MC, &b_pt_SE_MC);
   fChain->SetBranchAddress("eta_SE_MC", &eta_SE_MC, &b_eta_SE_MC);
   fChain->SetBranchAddress("phi_SE_MC", &phi_SE_MC, &b_phi_SE_MC);
   fChain->SetBranchAddress("pt_SE_Data", &pt_SE_Data, &b_pt_SE_Data);
   fChain->SetBranchAddress("eta_SE_Data", &eta_SE_Data, &b_eta_SE_Data);
   fChain->SetBranchAddress("phi_SE_Data", &phi_SE_Data, &b_phi_SE_Data);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("nGenEle", &nGenEle, &b_nGenEle);
   fChain->SetBranchAddress("ZMass", &ZMass, &b_ZMass);
   fChain->SetBranchAddress("ZPt", &ZPt, &b_ZPt);
   fChain->SetBranchAddress("ZEta", &ZEta, &b_ZEta);
   fChain->SetBranchAddress("ZRap", &ZRap, &b_ZRap);
   fChain->SetBranchAddress("ZPhi", &ZPhi, &b_ZPhi);
   fChain->SetBranchAddress("gPre_energy", &gPre_energy, &b_gPre_energy);
   fChain->SetBranchAddress("gPre_px", &gPre_px, &b_gPre_px);
   fChain->SetBranchAddress("gPre_py", &gPre_py, &b_gPre_py);
   fChain->SetBranchAddress("gPre_pz", &gPre_pz, &b_gPre_pz);
   fChain->SetBranchAddress("gPre_pt", &gPre_pt, &b_gPre_pt);
   fChain->SetBranchAddress("gPre_eta", &gPre_eta, &b_gPre_eta);
   fChain->SetBranchAddress("gPre_rap", &gPre_rap, &b_gPre_rap);
   fChain->SetBranchAddress("gPre_phi", &gPre_phi, &b_gPre_phi);
   fChain->SetBranchAddress("gPost_energy", &gPost_energy, &b_gPost_energy);
   fChain->SetBranchAddress("gPost_px", &gPost_px, &b_gPost_px);
   fChain->SetBranchAddress("gPost_py", &gPost_py, &b_gPost_py);
   fChain->SetBranchAddress("gPost_pz", &gPost_pz, &b_gPost_pz);
   fChain->SetBranchAddress("gPost_pt", &gPost_pt, &b_gPost_pt);
   fChain->SetBranchAddress("gPost_eta", &gPost_eta, &b_gPost_eta);
   fChain->SetBranchAddress("gPost_rap", &gPost_rap, &b_gPost_rap);
   fChain->SetBranchAddress("gPost_phi", &gPost_phi, &b_gPost_phi);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("rap", &rap, &b_rap);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("enSC", &enSC, &b_enSC);
   fChain->SetBranchAddress("preEnSC", &preEnSC, &b_preEnSC);
   fChain->SetBranchAddress("rawEnSC", &rawEnSC, &b_rawEnSC);
   fChain->SetBranchAddress("etSC", &etSC, &b_etSC);
   fChain->SetBranchAddress("etaSC", &etaSC, &b_etaSC);
   fChain->SetBranchAddress("phiSC", &phiSC, &b_phiSC);
   fChain->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta, &b_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("E1x5", &E1x5, &b_E1x5);
   fChain->SetBranchAddress("E2x5", &E2x5, &b_E2x5);
   fChain->SetBranchAddress("E5x5", &E5x5, &b_E5x5);
   fChain->SetBranchAddress("hOverE", &hOverE, &b_hOverE);
   fChain->SetBranchAddress("etaScWidth", &etaScWidth, &b_etaScWidth);
   fChain->SetBranchAddress("phiScWidth", &phiScWidth, &b_phiScWidth);
   fChain->SetBranchAddress("r9", &r9, &b_r9);
   fChain->SetBranchAddress("dEtaIn", &dEtaIn, &b_dEtaIn);
   fChain->SetBranchAddress("dPhiIn", &dPhiIn, &b_dPhiIn);
   fChain->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
   fChain->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
   fChain->SetBranchAddress("isoPhotons", &isoPhotons, &b_isoPhotons);
   fChain->SetBranchAddress("isoChargedFromPU", &isoChargedFromPU, &b_isoChargedFromPU);
   fChain->SetBranchAddress("isoDeltaBeta", &isoDeltaBeta, &b_isoDeltaBeta);
   fChain->SetBranchAddress("isoRho", &isoRho, &b_isoRho);
   fChain->SetBranchAddress("ooEmooP", &ooEmooP, &b_ooEmooP);
   fChain->SetBranchAddress("d0", &d0, &b_d0);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits, &b_expectedMissingInnerHits);
   fChain->SetBranchAddress("passConversionVeto", &passConversionVeto, &b_passConversionVeto);
   fChain->SetBranchAddress("brem", &brem, &b_brem);
   fChain->SetBranchAddress("isTrue", &isTrue, &b_isTrue);
   fChain->SetBranchAddress("passVetoId", &passVetoId, &b_passVetoId);
   fChain->SetBranchAddress("passLooseId", &passLooseId, &b_passLooseId);
   fChain->SetBranchAddress("passMediumId", &passMediumId, &b_passMediumId);
   fChain->SetBranchAddress("passTightId", &passTightId, &b_passTightId);
   fChain->SetBranchAddress("mvaVal", &mvaVal, &b_mvaVal);
   fChain->SetBranchAddress("mvaCat", &mvaCat, &b_mvaCat);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   Notify();
}

Bool_t mc13TeV::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mc13TeV::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mc13TeV::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mc13TeV_cxx
