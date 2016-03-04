//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 30 11:36:46 2016 by ROOT version 6.02/05
// from TTree ElectronTree/Electron data
// found on file: /tmp/rchawla/SingleElectron_Run2015D.root
//////////////////////////////////////////////////////////

#ifndef data13TeV_h
#define data13TeV_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "TClonesArray.h"

class data13TeV {
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
   Bool_t          single_Ele23;
   Bool_t          single_Ele27;
   Bool_t          singleMuon;
   Bool_t          doubleElectron;
   Bool_t          doubleEMu_17_8;
   Bool_t          doubleEMu_12_17;
   Bool_t          doubleEMu_23_8;
   vector<double>  *et_Photon;
   Int_t           singlePhoton;
   Int_t           prescalePhoton;
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
   vector<float>   *gen_postFSR_ene;
   vector<float>   *gen_postFSR_px;
   vector<float>   *gen_postFSR_py;
   vector<float>   *gen_postFSR_pz;
   vector<float>   *gen_postFSR_pt;
   vector<float>   *gen_postFSR_eta;
   vector<float>   *gen_postFSR_rap;
   vector<float>   *gen_postFSR_phi;
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
   vector<float>   *eleInBarrel;
   vector<float>   *eleInEndcap;
   vector<int>     *isGLBmuon;
   vector<int>     *isPFmuon;
   Int_t           nMuons;
   vector<bool>    *isLoose;
   vector<bool>    *isTight;
   vector<bool>    *isHighPt;
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
   vector<double>  *isoTrkMuon;
   vector<double>  *metPt;
   vector<double>  *metPhi;
   vector<double>  *metSumEt;
   Int_t           nPhotons;
   vector<float>   *ptPhoton;
   vector<float>   *etaPhoton;
   vector<float>   *phiPhoton;
   TClonesArray    *gen_preFSR;

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
   TBranch        *b_single_Ele23;   //!
   TBranch        *b_single_Ele27;   //!
   TBranch        *b_singleMuon;   //!
   TBranch        *b_doubleElectron;   //!
   TBranch        *b_doubleEMu_17_8;   //!
   TBranch        *b_doubleEMu_12_17;   //!
   TBranch        *b_doubleEMu_23_8;   //!
   TBranch        *b_et_Photon;   //!
   TBranch        *b_singlePhoton;   //!
   TBranch        *b_prescalePhoton;   //!
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
   TBranch        *b_gen_postFSR_ene;   //!
   TBranch        *b_gen_postFSR_px;   //!
   TBranch        *b_gen_postFSR_py;   //!
   TBranch        *b_gen_postFSR_pz;   //!
   TBranch        *b_gen_postFSR_pt;   //!
   TBranch        *b_gen_postFSR_eta;   //!
   TBranch        *b_gen_postFSR_rap;   //!
   TBranch        *b_gen_postFSR_phi;   //!
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
   TBranch        *b_eleInBarrel;   //!
   TBranch        *b_eleInEndcap;   //!
   TBranch        *b_isGLBmuon;   //!
   TBranch        *b_isPFmuon;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_isLoose;   //!
   TBranch        *b_isTight;   //!
   TBranch        *b_isHighPt;   //!
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
   TBranch        *b_isoTrkMuon;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_metSumEt;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_ptPhoton;   //!
   TBranch        *b_etaPhoton;   //!
   TBranch        *b_phiPhoton;   //!
   TBranch        *b_gen_preFSR;   //!

   data13TeV(TTree *tree=0);
   virtual ~data13TeV();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef data13TeV_cxx
data13TeV::data13TeV(TTree *tree) : fChain(0) 
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

data13TeV::~data13TeV()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t data13TeV::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t data13TeV::LoadTree(Long64_t entry)
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

void data13TeV::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   et_Photon = 0;
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
   gen_postFSR_ene = 0;
   gen_postFSR_px = 0;
   gen_postFSR_py = 0;
   gen_postFSR_pz = 0;
   gen_postFSR_pt = 0;
   gen_postFSR_eta = 0;
   gen_postFSR_rap = 0;
   gen_postFSR_phi = 0;
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
   eleInBarrel = 0;
   eleInEndcap = 0;
   isGLBmuon = 0;
   isPFmuon = 0;
   isLoose = 0;
   isTight = 0;
   isHighPt = 0;
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
   isoTrkMuon = 0;
   metPt = 0;
   metPhi = 0;
   metSumEt = 0;
   ptPhoton = 0;
   etaPhoton = 0;
   phiPhoton = 0;
   gen_preFSR = 0;
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
   fChain->SetBranchAddress("single_Ele23", &single_Ele23, &b_single_Ele23);
   fChain->SetBranchAddress("single_Ele27", &single_Ele27, &b_single_Ele27);
   fChain->SetBranchAddress("singleMuon", &singleMuon, &b_singleMuon);
   fChain->SetBranchAddress("doubleElectron", &doubleElectron, &b_doubleElectron);
   fChain->SetBranchAddress("doubleEMu_17_8", &doubleEMu_17_8, &b_doubleEMu_17_8);
   fChain->SetBranchAddress("doubleEMu_12_17", &doubleEMu_12_17, &b_doubleEMu_12_17);
   fChain->SetBranchAddress("doubleEMu_23_8", &doubleEMu_23_8, &b_doubleEMu_23_8);
   fChain->SetBranchAddress("et_Photon", &et_Photon, &b_et_Photon);
   fChain->SetBranchAddress("singlePhoton", &singlePhoton, &b_singlePhoton);
   fChain->SetBranchAddress("prescalePhoton", &prescalePhoton, &b_prescalePhoton);
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
   fChain->SetBranchAddress("gen_postFSR_ene", &gen_postFSR_ene, &b_gen_postFSR_ene);
   fChain->SetBranchAddress("gen_postFSR_px", &gen_postFSR_px, &b_gen_postFSR_px);
   fChain->SetBranchAddress("gen_postFSR_py", &gen_postFSR_py, &b_gen_postFSR_py);
   fChain->SetBranchAddress("gen_postFSR_pz", &gen_postFSR_pz, &b_gen_postFSR_pz);
   fChain->SetBranchAddress("gen_postFSR_pt", &gen_postFSR_pt, &b_gen_postFSR_pt);
   fChain->SetBranchAddress("gen_postFSR_eta", &gen_postFSR_eta, &b_gen_postFSR_eta);
   fChain->SetBranchAddress("gen_postFSR_rap", &gen_postFSR_rap, &b_gen_postFSR_rap);
   fChain->SetBranchAddress("gen_postFSR_phi", &gen_postFSR_phi, &b_gen_postFSR_phi);
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
   fChain->SetBranchAddress("eleInBarrel", &eleInBarrel, &b_eleInBarrel);
   fChain->SetBranchAddress("eleInEndcap", &eleInEndcap, &b_eleInEndcap);
   fChain->SetBranchAddress("isGLBmuon", &isGLBmuon, &b_isGLBmuon);
   fChain->SetBranchAddress("isPFmuon", &isPFmuon, &b_isPFmuon);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("isLoose", &isLoose, &b_isLoose);
   fChain->SetBranchAddress("isTight", &isTight, &b_isTight);
   fChain->SetBranchAddress("isHighPt", &isHighPt, &b_isHighPt);
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
   fChain->SetBranchAddress("isoTrkMuon", &isoTrkMuon, &b_isoTrkMuon);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("metSumEt", &metSumEt, &b_metSumEt);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("ptPhoton", &ptPhoton, &b_ptPhoton);
   fChain->SetBranchAddress("etaPhoton", &etaPhoton, &b_etaPhoton);
   fChain->SetBranchAddress("phiPhoton", &phiPhoton, &b_phiPhoton);
   fChain->SetBranchAddress("gen_preFSR", &gen_preFSR, &b_gen_preFSR);
   Notify();
}

Bool_t data13TeV::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void data13TeV::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t data13TeV::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef data13TeV_cxx
