//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 26 07:31:57 2016 by ROOT version 6.02/05
// from TTree ElectronTree/Electron data
// found on file: /tmp/rchawla/SingleElectron_Run2015D.root
//////////////////////////////////////////////////////////

#ifndef dataEff_vs_Mass_h
#define dataEff_vs_Mass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class dataEff_vs_Mass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        RunNo;
   Double_t        EvtNo;
   Double_t        Lumi;
   Double_t        Bunch;
   Double_t        theWeight;
   Int_t           pvNTracks;
   Int_t           nPV;
   Int_t           nPU;
   Int_t           nPUTrue;
   Float_t         rho;
   Bool_t          singleEle;
   Bool_t          singleMuon;
   Bool_t          doubleElectron;
   Bool_t          doubleEMu_17_8;
   Bool_t          doubleEMu_12_17;
   Bool_t          doubleEMu_23_8;
   vector<double>  *pt_leg1;
   vector<double>  *eta_leg1;
   vector<double>  *phi_leg1;
   vector<double>  *pt_leg2;
   vector<double>  *eta_leg2;
   vector<double>  *phi_leg2;
   vector<double>  *pt_Ele;
   vector<double>  *eta_Ele;
   vector<double>  *phi_Ele;
   Int_t           tauFlag;
   vector<float>   *gen_ptTau;
   vector<float>   *gen_etaTau;
   vector<float>   *gen_phiTau;
   Int_t           nEle;
   Int_t           nGenEle;
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
   vector<int>     *eleEcalDrivenSeed;
   vector<int>     *isGLBmuon;
   vector<int>     *isPFmuon;
   Int_t           nMuons;
   vector<bool>    *isLoose;
   vector<bool>    *isTight;
   vector<double>  *ptMuon;
   vector<double>  *etaMuon;
   vector<double>  *phiMuon;
   vector<double>  *energyMuon;
   vector<double>  *chargeMuon;
   vector<double>  *isoChargedHadronPfR04Muon;
   vector<double>  *isoNeutralHadronPfR04Muon;
   vector<double>  *isoGammaPfR04Muon;
   vector<double>  *isoChargedFromPUMuon;
   vector<double>  *isoPFMuon;

   // List of branches
   TBranch        *b_RunNo;   //!
   TBranch        *b_EvtNo;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_Bunch;   //!
   TBranch        *b_theWeight;   //!
   TBranch        *b_pvNTracks;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nPUTrue;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_singleEle;   //!
   TBranch        *b_singleMuon;   //!
   TBranch        *b_doubleElectron;   //!
   TBranch        *b_doubleEMu_17_8;   //!
   TBranch        *b_doubleEMu_12_17;   //!
   TBranch        *b_doubleEMu_23_8;   //!
   TBranch        *b_pt_leg1;   //!
   TBranch        *b_eta_leg1;   //!
   TBranch        *b_phi_leg1;   //!
   TBranch        *b_pt_leg2;   //!
   TBranch        *b_eta_leg2;   //!
   TBranch        *b_phi_leg2;   //!
   TBranch        *b_pt_Ele;   //!
   TBranch        *b_eta_Ele;   //!
   TBranch        *b_phi_Ele;   //!
   TBranch        *b_tauFlag;   //!
   TBranch        *b_gen_ptTau;   //!
   TBranch        *b_gen_etaTau;   //!
   TBranch        *b_gen_phiTau;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_nGenEle;   //!
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
   TBranch        *b_eleEcalDrivenSeed;   //!
   TBranch        *b_isGLBmuon;   //!
   TBranch        *b_isPFmuon;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_isLoose;   //!
   TBranch        *b_isTight;   //!
   TBranch        *b_ptMuon;   //!
   TBranch        *b_etaMuon;   //!
   TBranch        *b_phiMuon;   //!
   TBranch        *b_energyMuon;   //!
   TBranch        *b_chargeMuon;   //!
   TBranch        *b_isoChargedHadronPfR04Muon;   //!
   TBranch        *b_isoNeutralHadronPfR04Muon;   //!
   TBranch        *b_isoGammaPfR04Muon;   //!
   TBranch        *b_isoChargedFromPUMuon;   //!
   TBranch        *b_isoPFMuon;   //!

   dataEff_vs_Mass(TTree *tree=0);
   virtual ~dataEff_vs_Mass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef dataEff_vs_Mass_cxx
dataEff_vs_Mass::dataEff_vs_Mass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/rchawla/SingleElectron_Run2015D.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/tmp/rchawla/SingleElectron_Run2015D.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/tmp/rchawla/SingleElectron_Run2015D.root:/ntupler");
      dir->GetObject("ElectronTree",tree);

   }
   Init(tree);
}

dataEff_vs_Mass::~dataEff_vs_Mass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dataEff_vs_Mass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dataEff_vs_Mass::LoadTree(Long64_t entry)
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

void dataEff_vs_Mass::Init(TTree *tree)
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
   pt_Ele = 0;
   eta_Ele = 0;
   phi_Ele = 0;
   gen_ptTau = 0;
   gen_etaTau = 0;
   gen_phiTau = 0;
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
   eleEcalDrivenSeed = 0;
   isGLBmuon = 0;
   isPFmuon = 0;
   isLoose = 0;
   isTight = 0;
   ptMuon = 0;
   etaMuon = 0;
   phiMuon = 0;
   energyMuon = 0;
   chargeMuon = 0;
   isoChargedHadronPfR04Muon = 0;
   isoNeutralHadronPfR04Muon = 0;
   isoGammaPfR04Muon = 0;
   isoChargedFromPUMuon = 0;
   isoPFMuon = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
   fChain->SetBranchAddress("EvtNo", &EvtNo, &b_EvtNo);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("Bunch", &Bunch, &b_Bunch);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("pvNTracks", &pvNTracks, &b_pvNTracks);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nPUTrue", &nPUTrue, &b_nPUTrue);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("singleEle", &singleEle, &b_singleEle);
   fChain->SetBranchAddress("singleMuon", &singleMuon, &b_singleMuon);
   fChain->SetBranchAddress("doubleElectron", &doubleElectron, &b_doubleElectron);
   fChain->SetBranchAddress("doubleEMu_17_8", &doubleEMu_17_8, &b_doubleEMu_17_8);
   fChain->SetBranchAddress("doubleEMu_12_17", &doubleEMu_12_17, &b_doubleEMu_12_17);
   fChain->SetBranchAddress("doubleEMu_23_8", &doubleEMu_23_8, &b_doubleEMu_23_8);
   fChain->SetBranchAddress("pt_leg1", &pt_leg1, &b_pt_leg1);
   fChain->SetBranchAddress("eta_leg1", &eta_leg1, &b_eta_leg1);
   fChain->SetBranchAddress("phi_leg1", &phi_leg1, &b_phi_leg1);
   fChain->SetBranchAddress("pt_leg2", &pt_leg2, &b_pt_leg2);
   fChain->SetBranchAddress("eta_leg2", &eta_leg2, &b_eta_leg2);
   fChain->SetBranchAddress("phi_leg2", &phi_leg2, &b_phi_leg2);
   fChain->SetBranchAddress("pt_Ele", &pt_Ele, &b_pt_Ele);
   fChain->SetBranchAddress("eta_Ele", &eta_Ele, &b_eta_Ele);
   fChain->SetBranchAddress("phi_Ele", &phi_Ele, &b_phi_Ele);
   fChain->SetBranchAddress("tauFlag", &tauFlag, &b_tauFlag);
   fChain->SetBranchAddress("gen_ptTau", &gen_ptTau, &b_gen_ptTau);
   fChain->SetBranchAddress("gen_etaTau", &gen_etaTau, &b_gen_etaTau);
   fChain->SetBranchAddress("gen_phiTau", &gen_phiTau, &b_gen_phiTau);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("nGenEle", &nGenEle, &b_nGenEle);
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
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("isGLBmuon", &isGLBmuon, &b_isGLBmuon);
   fChain->SetBranchAddress("isPFmuon", &isPFmuon, &b_isPFmuon);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("isLoose", &isLoose, &b_isLoose);
   fChain->SetBranchAddress("isTight", &isTight, &b_isTight);
   fChain->SetBranchAddress("ptMuon", &ptMuon, &b_ptMuon);
   fChain->SetBranchAddress("etaMuon", &etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", &phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("energyMuon", &energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("chargeMuon", &chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("isoChargedHadronPfR04Muon", &isoChargedHadronPfR04Muon, &b_isoChargedHadronPfR04Muon);
   fChain->SetBranchAddress("isoNeutralHadronPfR04Muon", &isoNeutralHadronPfR04Muon, &b_isoNeutralHadronPfR04Muon);
   fChain->SetBranchAddress("isoGammaPfR04Muon", &isoGammaPfR04Muon, &b_isoGammaPfR04Muon);
   fChain->SetBranchAddress("isoChargedFromPUMuon", &isoChargedFromPUMuon, &b_isoChargedFromPUMuon);
   fChain->SetBranchAddress("isoPFMuon", &isoPFMuon, &b_isoPFMuon);
   Notify();
}

Bool_t dataEff_vs_Mass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dataEff_vs_Mass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t dataEff_vs_Mass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef dataEff_vs_Mass_cxx
