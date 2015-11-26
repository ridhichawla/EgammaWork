#define mc13TeV_cxx
#include "mc13TeV.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
//#include "DataFormats/Math/interface/deltaR.h"

Double_t mc13TeV::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

Double_t mc13TeV::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t mc13TeV::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

void mc13TeV::Loop()
{
  TFile *f1 = TFile::Open("/afs/cern.ch/work/r/rchawla/private/Analysis_ee_13TeV/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/PileUp/dyM50_pileup.root");
  //TFile *f2 = TFile::Open("/afs/cern.ch/work/r/rchawla/private/CMSSW_5_3_22/src/Analysis/El_analyzer/Macros/New/MCEE/WTSfiles/MCDist_M50_13TeV.root");

  // data histogram 
  TH1F *DATA_puDist = (TH1F*)f1->Get("DATA_nPV");
  //DATA_puDist->Scale(1/DATA_puDist->Integral());

  // mc histogram 
  TH1F *MC_puDist = (TH1F*)f1->Get("MC_nPV");
  //MC_puDist->Scale(1/MC_puDist->Integral());

  TH1F *weights_ = (TH1F*)DATA_puDist->Clone("weights_");
  weights_->Divide(MC_puDist);

  TFile *file = new TFile("dyJtoLL_M50.root", "recreate");
  //std::ofstream wgts;
  //wgts.open ("Weights.txt");

  //double Mz = 91.1876;

  int kin, id;
  bool isCharge, isMass;
  int good_elec, gen_elec_pre, gen_elec_post;
  int n_ele, n_trig, n_kin, n_id, n_kin_id, n_charge, n_all, n_all_cuts, nZ;
  int n_dR;
  double dR;
  double sum_weights;
  double Z_Mass, Z_Y, Z_Rap, Z_Eta, Z_Pt, Z_Phi, ZMass_G_pre, ZMass_G_post, ZMass_R_pre, ZMass_R_post, ZMass_di_gPre, ZMass_di_gPost;
  TLorentzVector ele1,ele2,dielectron;
  TLorentzVector reco1_pre,reco2_pre,direco_pre,reco1_post,reco2_post,direco_post;
  TLorentzVector gen1_pre,gen2_pre,digen_pre,gen1_post,gen2_post,digen_post;
  TLorentzVector gPre1,gPre2,di_gPre;
  TLorentzVector gPost1,gPost2,di_gPost;

  vector <int> indx;
  
  vector <double> newscEt;
  vector <double> newscEta;
  vector <double> newscPhi;
  vector <double> newscEnr;
  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> neweleEnr;
  vector <double> newelePhi;
  vector <double> neweleCharge;

  vector <double> newgPre_Pt;
  vector <double> newgPre_Eta;
  vector <double> newgPre_Enr;
  vector <double> newgPre_Phi;

  vector <double> newgPost_Pt;
  vector <double> newgPost_Eta;
  vector <double> newgPost_Enr;
  vector <double> newgPost_Phi;

  vector <double> recomatchedGen;
  vector <double> Recomatchedgen;
  vector <double> genmatchedReco;
  vector <double> Genmatchedreco;

  vector <double> recoPre_Pt;
  vector <double> recoPre_Eta;
  vector <double> recoPre_Enr;
  vector <double> recoPre_Phi;

  vector <double> genPre_Pt;
  vector <double> genPre_Eta;
  vector <double> genPre_Enr;
  vector <double> genPre_Phi;

  vector <double> recoPost_Pt;
  vector <double> recoPost_Eta;
  vector <double> recoPost_Enr;
  vector <double> recoPost_Phi;

  vector <double> genPost_Pt;
  vector <double> genPost_Eta;
  vector <double> genPost_Enr;
  vector <double> genPost_Phi;

  TH1F *MC_nPV = new TH1F("MC_nPV", "MC_nPV", 50, 0, 50);
  MC_nPV->Sumw2();

  //TH1F *MC_puDistW = new TH1F("MC_puDistW","MC_puDistW",50,0,50);
  //MC_puDistW->Sumw2();

  TH1F *pt_ = new TH1F("pt_", "pt_", 100, 0, 700);
  TH1F *eta_ = new TH1F("eta_", "eta_", 50, -2.5, 2.5);
  TH1F *et_sc_ = new TH1F("et_sc_", "et_sc_", 100, 0 ,700);
  TH1F *eta_sc_ = new TH1F("eta_sc_", "eta_sc_", 50, -2.5, 2.5);

  TH1F *pt_lead = new TH1F("pt_lead", "pt_lead", 100, 0, 700);
  TH1F *eta_lead = new TH1F("eta_lead", "eta_lead", 50, -2.5, 2.5);
  TH1F *pt_slead = new TH1F("pt_slead", "pt_slead", 100, 0, 700);
  TH1F *eta_slead = new TH1F("eta_slead", "eta_slead", 48, -2.4, 2.4);

  TH1F *et_sc_lead = new TH1F("et_sc_lead", "et_sc_lead", 105, 0 ,700);
  TH1F *eta_sc_lead = new TH1F("eta_sc_lead", "eta_sc_lead", 48, -2.4, 2.4);
  TH1F *et_sc_slead = new TH1F("et_sc_slead", "et_sc_slead", 105, 0 ,700);
  TH1F *eta_sc_slead = new TH1F("eta_sc_slead", "eta_sc_slead", 48, -2.4, 2.4);

  pt_->Sumw2(); eta_->Sumw2(); et_sc_->Sumw2(); eta_sc_->Sumw2();
  pt_lead->Sumw2(); eta_lead->Sumw2(); pt_slead->Sumw2(); eta_slead->Sumw2();
  et_sc_lead->Sumw2(); eta_sc_lead->Sumw2(); et_sc_slead->Sumw2(); eta_sc_slead->Sumw2();

  TH1F *B_DeltaEta = new TH1F("B_DeltaEta", "B_deltaEta", 100, -0.01, 0.01);
  TH1F *B_DeltaPhi = new TH1F("B_DeltaPhi", "B_deltaPhi", 80, -0.08, 0.08);
  TH1F *B_H_E = new TH1F("B_H_E", "B_H_E", 80, 0, 0.1);
  TH1F *B_SigmaEta = new TH1F("B_SigmaEta", "B_sigmaEta", 100, 0, 0.1);
  TH1F *B_E_P = new TH1F("B_E_P", "B_E_P", 100, 0, 0.5);
  TH1F *B_D0 = new TH1F("B_D0", "B_D0", 80, 0, 0.2);
  TH1F *B_Dz = new TH1F("B_Dz", "B_Dz", 80, 0, 0.2);
  TH1F *B_MissHits = new TH1F("B_MissHits", "B_MissHits", 4, 0, 4);
  TH1F *B_ConVeto = new TH1F("B_ConVeto", "B_ConVeto", 2, 0, 2);
  TH1F *B_fBrem = new TH1F("B_fBrem", "B_fBrem", 40, 0, 1);
  TH1F *B_r9 = new TH1F("B_r9", "B_r9", 40, 0, 1);
  TH1F *B_RelRhoIso = new TH1F("B_RelRhoIso", "B_RelRhoIso", 80, 0, 0.6);
  //TH1F *B_RelEMIso = new TH1F("B_RelEMIso", "B_RelEMIso", 60, 0, 0.6);
  //TH1F *B_RelNeutralIso = new TH1F("B_RelNeutralIso", "B_RelNeutralIso", 60, 0, 0.6);

  B_DeltaEta->Sumw2(); B_DeltaPhi->Sumw2(); B_H_E->Sumw2(); B_SigmaEta->Sumw2(); B_E_P->Sumw2(); B_D0->Sumw2(); B_Dz->Sumw2();
  B_MissHits->Sumw2(); B_ConVeto->Sumw2(); B_fBrem->Sumw2(); B_r9->Sumw2(); B_RelRhoIso->Sumw2();

  TH1F *E_DeltaEta = new TH1F("E_DeltaEta", "E_DeltaEta", 100, -0.01, 0.01);
  TH1F *E_DeltaPhi = new TH1F("E_DeltaPhi", "E_DeltaPhi", 80, -0.08, 0.08);
  TH1F *E_H_E = new TH1F("E_H_E", "E_H_E", 80, 0, 0.1);
  TH1F *E_SigmaEta = new TH1F("E_SigmaEta", "E_SigmaEta", 100, 0.0, 0.1);
  TH1F *E_E_P = new TH1F("E_E_P", "E_E_P", 100, 0, 0.5);
  TH1F *E_D0 = new TH1F("E_D0", "E_D0", 80, 0, 0.2);
  TH1F *E_Dz = new TH1F("E_Dz", "E_Dz", 80, 0, 0.2);
  TH1F *E_MissHits = new TH1F("E_MissHits", "E_MissHits", 4, 0, 4);
  TH1F *E_ConVeto = new TH1F("E_ConVeto", "E_ConVeto", 2, 0, 2);
  TH1F *E_fBrem = new TH1F("E_fBrem", "E_fBrem", 40, 0, 1);
  TH1F *E_r9 = new TH1F("E_r9", "E_r9", 40, 0.0, 1);
  TH1F *E_RelRhoIso = new TH1F("E_RelRhoIso", "E_RelRhoIso", 80, 0, 0.6);
  //TH1F *E_RelEMIso = new TH1F("E_RelEMIso", "E_RelEMIso", 60, 0, 0.6);
  //TH1F *E_RelNeutralIso = new TH1F("E_RelNeutralIso", "E_RelNeutralIso", 60, 0, 0.6);

  E_DeltaEta->Sumw2(); E_DeltaPhi->Sumw2(); E_H_E->Sumw2(); E_SigmaEta->Sumw2(); E_E_P->Sumw2(); E_D0->Sumw2(); E_Dz->Sumw2();
  E_MissHits->Sumw2(); E_ConVeto->Sumw2(); E_fBrem->Sumw2(); E_r9->Sumw2(); E_RelRhoIso->Sumw2();

  //int n = 4;
  //const double xbins[32] = {45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510};
  //const double xbins[19] = {45., 50., 55., 60., 64., 68., 72., 76., 81., 86., 91., 96., 101., 106., 110., 115., 120., 126., 133.};
  
  Double_t xbins[42] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 1000, 1500, 2000};
  TH1D *genPre_Mass_ = new TH1D("genPre_Mass_", "genPre_Mass_", 41, xbins);
  //TH1D *recoPre_Mass = new TH1D("recoPre_Mass", "recoPre_Mass", 41, xbins);
  TH2D *hRM_pre_RG = new TH2D("Response_Matrix_pre_RG", "Response_Matrix_pre_RG", 41, xbins, 41, xbins);
  TH2D *hRM_pre_GR = new TH2D("Response_Matrix_pre_GR", "Response_Matrix_pre_GR", 41, xbins, 41, xbins);
  TH2D *hRM_noWeight_pre_RG = new TH2D("Response_Matrix_noWeight_pre_RG", "Response_Matrix_noWeight_pre_RG", 41, xbins, 41, xbins);
  TH2D *hRM_noWeight_pre_GR = new TH2D("Response_Matrix_noWeight_pre_GR", "Response_Matrix_noWeight_pre_GR", 41, xbins, 41, xbins);

  TH1D *genPost_Mass_ = new TH1D("genPost_Mass_", "genPost_Mass_", 41, xbins);
  //TH1D *recoPost_Mass = new TH1D("recoPost_Mass", "recoPost_Mass", 41, xbins);
  TH2D *hRM_post_RG = new TH2D("Response_Matrix_post_RG", "Response_Matrix_post_RG", 41, xbins, 41, xbins);
  TH2D *hRM_post_GR = new TH2D("Response_Matrix_post_GR", "Response_Matrix_post_GR", 41, xbins, 41, xbins);
  TH2D *hRM_noWeight_post_RG = new TH2D("Response_Matrix_noWeight_post_RG", "Response_Matrix_noWeight_post_RG", 41, xbins, 41, xbins);
  TH2D *hRM_noWeight_post_GR = new TH2D("Response_Matrix_noWeight_post_GR", "Response_Matrix_noWeight_post_GR", 41, xbins, 41, xbins);

  TH1D *ZMass_ = new TH1D("ZMass_", "ZMass_", 41, xbins);
  TH1D *ZMass_0to400 = new TH1D("ZMass_0to400", "ZMass_0to400", 100, 0, 400);
  TH1D *ZMass_60to120 = new TH1D("ZMass_60to120", "ZMass_60to120", 60, 60, 120);
  TH1D *ZPt_ = new TH1D("ZPt_", "ZPt_", 100, 0, 700);
  TH1D *ZY_ = new TH1D("ZY_", "ZY_", 150, -15, 15);
  TH1D *ZRap_ = new TH1D("ZRap_", "ZRap_", 60, -3, 3);
  TH1D *ZEta_ = new TH1D("ZEta_", "ZEta_", 60, -8, 8);
  TH1D *ZPhi_ = new TH1D("ZPhi_", "ZPhi_", 50, -3.5, 3.5);

  /*TH1D *ZMass_gPre_ = new TH1D("ZMass_gPre_", "ZMass_gPre_", 31, xbins);
    TH1D *ZPt_gPre_ = new TH1D("ZPt_gPre_", "ZPt_gPre_", 100, 0, 700);
  //TH1D *ZY_gPre_ = new TH1D("ZY_gPre_", "ZY_gPre_", 150, -15, 15);
  TH1D *ZRap_gPre_ = new TH1D("ZRap_gPre_", "ZRap_gPre_", 60, -3, 3);
  TH1D *ZEta_gPre_ = new TH1D("ZEta_gPre_", "ZEta_gPre_", 60, -8, 8);
  TH1D *ZPhi_gPre_ = new TH1D("ZPhi_gPre_", "ZPhi_gPre_", 50, -3.5, 3.5);*/

  //TH1D *ZMass_di_gPre_ = new TH1D("ZMass_di_gPre_", "ZMass_di_gPre_", 31, xbins);
  TH1D *pt_gPre_ = new TH1D("pt_gPre_", "pt_gPre_", 100, 0, 400);
  TH1D *eta_gPre_ = new TH1D("eta_gPre_", "eta_gPre_", 50, -2.5, 2.5);

  //TH1D *ZMass_di_gPost_ = new TH1D("ZMass_di_gPost_", "ZMass_di_gPost_", 31, xbins);
  TH1D *pt_gPost_ = new TH1D("pt_gPost_", "pt_gPost_", 100, 0, 400);
  TH1D *eta_gPost_ = new TH1D("eta_gPost_", "eta_gPost_", 50, -2.5, 2.5);

  TH1D *genPreSize_  = new TH1D("genPreSize_", "genPreSize_", 10, 0, 10);
  TH1D *genPostSize_ = new TH1D("genPostSize_", "genPostSize_", 10, 0, 10);

  hRM_pre_RG->Sumw2(); hRM_noWeight_pre_RG->Sumw2(); hRM_pre_GR->Sumw2(); hRM_noWeight_pre_GR->Sumw2(); genPre_Mass_->Sumw2();
  hRM_post_RG->Sumw2(); hRM_noWeight_post_RG->Sumw2(); hRM_post_GR->Sumw2(); hRM_noWeight_post_GR->Sumw2(); genPost_Mass_->Sumw2();
  pt_gPre_->Sumw2(); eta_gPre_->Sumw2(); 
  pt_gPost_->Sumw2(); eta_gPost_->Sumw2(); 
  genPreSize_->Sumw2(); genPostSize_->Sumw2(); 

  ZMass_->Sumw2(); ZMass_0to400->Sumw2(); ZMass_60to120->Sumw2(); ZPt_->Sumw2(); ZY_->Sumw2(); ZRap_->Sumw2(); ZEta_->Sumw2(); ZPhi_->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 50;
  cout<<"entries: "<<nentries<<endl;

  n_ele = 0;
  n_trig = 0;
  n_kin = 0;
  n_id = 0;
  n_charge = 0;
  n_all = 0;
  n_dR = 0;
  n_all_cuts = 0;
  nZ = 0;
  sum_weights = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    kin = 0;
    id = 0;
    isCharge = false;
    isMass = false;

    indx.clear();

    int index1[gPre_pt->size()];
    float pt1[gPre_pt->size()];

    for(unsigned int a=0; a<gPre_pt->size(); a++)
    { 
      pt1[a]=gPre_pt->at(a);
    }

    int sizem = sizeof(pt1)/sizeof(pt1[0]);
    TMath::Sort(sizem,pt1,index1,true);

    int index2[gPost_pt->size()];
    float pt2[gPost_pt->size()];

    for(unsigned int b=0; b<gPost_pt->size(); b++)
    {
      pt2[b]=gPost_pt->at(b);
    }

    int sizen = sizeof(pt2)/sizeof(pt2[0]);
    TMath::Sort(sizen,pt2,index2,true);

    //cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //cout<<""<<endl;
    Z_Mass=0.; Z_Y=0.; Z_Rap=0.; Z_Eta=0.; Z_Pt=0.; Z_Phi=0.; ZMass_G_pre=0.; ZMass_R_pre=0.; ZMass_G_post=0.; ZMass_R_post=0.; ZMass_di_gPre=0.; ZMass_di_gPost=0.;
    good_elec = 0;
    gen_elec_pre = 0;
    gen_elec_post = 0;
    dR = 0;
    int bin = 0;
    double puWeights = 0.0;
    unsigned int matchGen1 = 0;
    unsigned int matchGen2 = 0;
    unsigned int matchReco1 = 0;
    unsigned int matchReco2 = 0;

    bin = weights_->GetXaxis()->FindBin(nPV);
    puWeights = weights_->GetBinContent(bin);
    //cout<<"puWeights: "<<puWeights<<endl;
    //MC_puDistW->Fill(nPV,puWeights);

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear(); newelePt.clear();
    neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    newgPre_Pt.clear(); newgPre_Eta.clear(); newgPre_Enr.clear(); newgPre_Phi.clear(); 
    newgPost_Pt.clear(); newgPost_Eta.clear(); newgPost_Enr.clear(); newgPost_Phi.clear();

    recomatchedGen.clear(); Recomatchedgen.clear(); genmatchedReco.clear(); Genmatchedreco.clear();

    recoPre_Pt.clear(); recoPre_Eta.clear(); recoPre_Enr.clear(); recoPre_Phi.clear();
    genPre_Pt.clear();  genPre_Eta.clear(); genPre_Enr.clear(); genPre_Phi.clear();

    recoPost_Pt.clear(); recoPost_Eta.clear(); recoPost_Enr.clear(); recoPost_Phi.clear(); 
    genPost_Pt.clear(); genPost_Eta.clear(); genPost_Enr.clear(); genPost_Phi.clear();

    MC_nPV->Fill(nPV,puWeights);
    //cout<<"NLO Weights: "<<theWeight<<endl;
    sum_weights = sum_weights+theWeight;

    if(nEle>=2) {
      ++n_ele;
      if(singleElectron_mc) {
	++n_trig;
	if(singleElectron_mc == 0) cout<<"Wrong"<<endl;

	for(int k=0;k<nEle;k++){

	  //bool isBarrel = (fabs(etaSC->at(k)) <= 1.479);
	  //bool isEndcap = (fabs(etaSC->at(k)) > 1.479 && fabs(etaSC->at(k)) < 2.5);

	  if(fabs(eta->at(k)) < 2.5 && pt->at(k) > 25. && !(fabs(etaSC->at(k)) > 1.4442 && fabs(etaSC->at(k)) < 1.566)){
	    kin++;

	    if(passMediumId->at(k) == 1 && eleEcalDrivenSeed->at(k) == 1){
	      id++;
	      indx.push_back(k);

	      if(passMediumId->at(k) == 0) cout<<"Wrong ID: "<<endl;
	      if(eleEcalDrivenSeed->at(k) == 0) cout<<"Wrong ECAL ID: "<<endl;

	      good_elec = good_elec + 1;

	      newscEt.push_back(etSC->at(k));
	      newscEta.push_back(etaSC->at(k));
	      newscPhi.push_back(phiSC->at(k));
	      newscEnr.push_back(enSC->at(k));
	      newelePt.push_back(pt->at(k));
	      neweleEta.push_back(eta->at(k));
	      neweleEnr.push_back(energy->at(k));
	      newelePhi.push_back(phi->at(k));
	      neweleCharge.push_back(charge->at(k));
	    }
	  }
	}
      }
    }

    if(good_elec>=2){
      if(neweleCharge[0] * neweleCharge[1] == -1){
	isCharge = true;
	++n_all_cuts;
	
	for(unsigned int i=0; i<newelePt.size(); i++)
	{
	  pt_->Fill(newelePt[i],theWeight*puWeights);
	  eta_->Fill(neweleEta[i],theWeight*puWeights);
	  et_sc_->Fill(newscEt[i],theWeight*puWeights);
	  eta_sc_->Fill(newscEta[i],theWeight*puWeights);
	}
	
	pt_lead->Fill(newelePt[0],theWeight*puWeights);
	eta_lead->Fill(neweleEta[0],theWeight*puWeights);
	pt_slead->Fill(newelePt[1],theWeight*puWeights);
	eta_slead->Fill(neweleEta[1],theWeight*puWeights);

	et_sc_lead->Fill(newscEt[0],theWeight*puWeights);
	et_sc_slead->Fill(newscEt[1],theWeight*puWeights);
	eta_sc_lead->Fill(newscEta[0],theWeight*puWeights);
	eta_sc_slead->Fill(newscEta[1],theWeight*puWeights);

	ele1.SetPtEtaPhiE(newscEt[0],newscEta[0],newscPhi[0],newscEnr[0]);
	ele2.SetPtEtaPhiE(newscEt[1],newscEta[1],newscPhi[1],newscEnr[1]);

	dielectron=ele1+ele2;
	Z_Mass = dielectron.M();

	if(Z_Mass >= 60 && Z_Mass <= 120){
	  isMass = true;
	  ++nZ;
	}
	
	Z_Pt = dielectron.Pt();
	Z_Y = dielectron.Y();
	Z_Rap = dielectron.Rapidity();
	Z_Eta = dielectron.Eta();
	Z_Phi = dielectron.Phi();

	ZMass_0to400->Fill(Z_Mass,theWeight*puWeights);
	ZMass_60to120->Fill(Z_Mass,theWeight*puWeights);
	ZMass_->Fill(Z_Mass,theWeight*puWeights);
	ZPt_->Fill(Z_Pt,theWeight*puWeights);
	ZY_->Fill(Z_Y,theWeight*puWeights);
	ZRap_->Fill(Z_Rap,theWeight*puWeights);
	ZEta_->Fill(Z_Eta,theWeight*puWeights);
	ZPhi_->Fill(Z_Phi,theWeight*puWeights);
      } 
    }

    if(kin>=2){
      n_kin++;
      if(id>=2){
	n_id++;
	if(isCharge){
	  n_charge++;
	  if(isMass){
	    n_all++;
	  } //isMass
	} // isCharge
      } // id
    } // kin


    //cout<<"event: "<<jentry<<endl;
    //cout<<"#####################################################"<<endl;

    genPreSize_->Fill(gPre_eta->size());
    genPostSize_->Fill(gPost_eta->size());

    if(gPre_pt->size()>0.){
      for(unsigned int l=0;l<gPre_eta->size();l++){

	pt_gPre_->Fill(gPre_pt->at(index1[l]));
	eta_gPre_->Fill(gPre_eta->at(index1[l]));

	if(fabs(gPre_eta->at(index1[l])) < 2.5 && gPre_pt->at(index1[l]) > 25. && !(fabs(gPre_eta->at(index1[l])) > 1.4442 && fabs(gPre_eta->at(index1[l])) < 1.566)){

	  gen_elec_pre = gen_elec_pre + 1;

	  newgPre_Pt.push_back(gPre_pt->at(index1[l]));
	  newgPre_Eta.push_back(gPre_eta->at(index1[l]));
	  newgPre_Enr.push_back(gPre_energy->at(index1[l]));
	  newgPre_Phi.push_back(gPre_phi->at(index1[l]));
	}
      }
    }
    //cout<<"2"<<endl;
    /*if(gPre_pt->size()>=2){
      cout<<"entry: "<<jentry<<endl;
      cout<<""<<endl;
      cout<<"gPre_pt->at(0): "<<gPre_pt->at(index1[0])<<"   "<<"gPre_eta->at(0): "<<gPre_eta->at(index1[0])<<endl;
      cout<<"gPre_pt->at(1): "<<gPre_pt->at(index1[1])<<"   "<<"gPre_eta->at(1): "<<gPre_eta->at(index1[1])<<endl;
      cout<<"newgPre_Pt size: "<<newgPre_Pt.size()<<endl;}*/

    if(gen_elec_pre>=2.){
      gPre1.SetPtEtaPhiE(newgPre_Pt.at(0),newgPre_Eta.at(0),newgPre_Phi.at(0),newgPre_Enr.at(0));
      gPre2.SetPtEtaPhiE(newgPre_Pt.at(1),newgPre_Eta.at(1),newgPre_Phi.at(1),newgPre_Enr.at(1));

      di_gPre=gPre1+gPre2;
      ZMass_di_gPre=di_gPre.M();
      genPre_Mass_->Fill(ZMass_di_gPre,theWeight*puWeights);
    }

    //cout<<"****************************************************"<<endl;
    //cout<<""<<endl;
    //cout<<"3"<<endl;
    if(gPost_pt->size()>0.){   
      for(unsigned int m=0;m<gPost_eta->size();m++){

	pt_gPost_->Fill(gPost_pt->at(index2[m]));
	eta_gPost_->Fill(gPost_eta->at(index2[m]));

	if(fabs(gPost_eta->at(index2[m])) < 2.5 && gPost_pt->at(index2[m]) > 25. && !(fabs(gPost_eta->at(index2[m])) > 1.4442 && fabs(gPost_eta->at(index2[m])) < 1.566)){

	  gen_elec_post = gen_elec_post + 1;

	  newgPost_Pt.push_back(gPost_pt->at(index2[m]));
	  newgPost_Eta.push_back(gPost_eta->at(index2[m]));
	  newgPost_Enr.push_back(gPost_energy->at(index2[m]));
	  newgPost_Phi.push_back(gPost_phi->at(index2[m]));
	}
      }
    }       
    //cout<<"4"<<endl;
    /*if(gPost_pt->size()>=2){
    //cout<<"entry: "<<jentry<<endl;
    cout<<""<<endl;
    cout<<"gen_elec_post: "<<gen_elec_post<<endl;
    cout<<"gPost_pt size: "<<gPost_pt->size()<<endl;
    cout<<"gPost_pt->at(0): "<<gPost_pt->at(index2[0])<<"   "<<"gPost_eta->at(0): "<<gPost_eta->at(index2[0])<<endl;
    cout<<"gPost_pt->at(1): "<<gPost_pt->at(index2[1])<<"   "<<"gPost_eta->at(1): "<<gPost_eta->at(index2[1])<<endl;
    //cout<<"gPost_pt->at(2): "<<gPost_pt->at(index2[2])<<"   "<<"gPost_eta->at(2): "<<gPost_eta->at(index2[2])<<endl;
    cout<<"newgPost_Pt size: "<<newgPost_Pt.size()<<endl;}*/

    //cout<<"****************************************************"<<endl;
    //cout<<""<<endl;

    if(gen_elec_post>=2.){
      gPost1.SetPtEtaPhiE(newgPost_Pt.at(0),newgPost_Eta.at(0),newgPost_Phi.at(0),newgPost_Enr.at(0));
      gPost2.SetPtEtaPhiE(newgPost_Pt.at(1),newgPost_Eta.at(1),newgPost_Phi.at(1),newgPost_Enr.at(1));

      di_gPost=gPost1+gPost2;
      ZMass_di_gPost=di_gPost.M();
      genPost_Mass_->Fill(ZMass_di_gPost,theWeight*puWeights);
    }
    //cout<<"5"<<endl;

    //if(good_elec<2) continue;
    //if(gen_elec_pre<2) continue;
    //if(gen_elec_post<2) continue;

    if(gen_elec_pre>=2){
      for(unsigned int ireco = 0; ireco < neweleEta.size(); ireco++){
	double dR_comp_pre = 1000.;
	for(unsigned int igen = 0; igen < newgPre_Eta.size(); igen++){
	  dR = deltaR(neweleEta[ireco], newelePhi[ireco], newgPre_Eta[igen], newgPre_Phi[igen]);
	  //cout<<"event: "<<jentry<<"   "<<"Reco : "<<ireco<<"   "<<"Gen : "<<igen<<" Delta R : "<<dR<<endl;

	  if(dR < 0.1  ){
	    if (dR < dR_comp_pre)
	    {
	      dR_comp_pre = dR;
	      matchGen1 = igen ; matchReco1 = ireco ; //iter[ireco]=igen;
	      recomatchedGen.push_back(igen);
	      Recomatchedgen.push_back(ireco);
	    }

	    //cout<<"newPre_Pt: "<<newgPre_Pt[0]<<"   "<<"newPre_Eta: "<<newgPre_Eta[0]<<endl;
	    //cout<<"elePre_Pt: "<<newelePt[0]<<"   "<<"elePre_Eta: "<<neweleEta[0]<<endl;

	    genPre_Pt.push_back(newgPre_Pt[matchGen1]); genPre_Eta.push_back(newgPre_Eta[matchGen1]); genPre_Enr.push_back(newgPre_Enr[matchGen1]); genPre_Phi.push_back(newgPre_Phi[matchGen1]);
	    recoPre_Pt.push_back(newelePt[matchReco1]); recoPre_Eta.push_back(neweleEta[matchReco1]); recoPre_Enr.push_back(neweleEnr[matchReco1]); recoPre_Phi.push_back(newelePhi[matchReco1]);
	  }
	}
      }
    }

    if(recoPre_Pt.size()>=2)
    {
      //cout<<"recoPre_Pt: "<<recoPre_Pt[0]<<"   "<<"recoPre_Eta: "<<recoPre_Eta[0]<<endl;

      reco1_pre.SetPtEtaPhiE(recoPre_Pt[0],recoPre_Eta[0],recoPre_Phi[0],recoPre_Enr[0]);
      reco2_pre.SetPtEtaPhiE(recoPre_Pt[1],recoPre_Eta[1],recoPre_Phi[1],recoPre_Enr[1]);
      direco_pre=reco1_pre+reco2_pre;
      ZMass_R_pre = direco_pre.M();
    }

    if(genPre_Pt.size()>=2)
    {
      //cout<<"genPre_Pt: "<<genPre_Pt[0]<<"   "<<"genPre_Eta: "<<genPre_Eta[0]<<endl;

      gen1_pre.SetPtEtaPhiE(genPre_Pt[0],genPre_Eta[0],genPre_Phi[0],genPre_Enr[0]);
      gen2_pre.SetPtEtaPhiE(genPre_Pt[1],genPre_Eta[1],genPre_Phi[1],genPre_Enr[1]);
      digen_pre=gen1_pre+gen2_pre;
      ZMass_G_pre = digen_pre.M();
    }

    //genPre_Mass->Fill(ZMass_G_pre,theWeight*puWeights);
    //recoPre_Mass->Fill(ZMass_R_pre,theWeight*puWeights);
    hRM_pre_RG->Fill(ZMass_R_pre,ZMass_G_pre,theWeight*puWeights);
    hRM_pre_GR->Fill(ZMass_G_pre,ZMass_R_pre,theWeight*puWeights);
    
    hRM_noWeight_pre_RG->Fill(ZMass_R_pre,ZMass_G_pre);
    hRM_noWeight_pre_GR->Fill(ZMass_G_pre,ZMass_R_pre);

    if(gen_elec_post>=2){
      for(unsigned int ireco = 0; ireco < neweleEta.size(); ireco++){
	double dR_comp_post = 1000.;
	for(unsigned int igen = 0; igen < newgPost_Eta.size(); igen++){
	  dR = deltaR(neweleEta[ireco], newelePhi[ireco], newgPost_Eta[igen], newgPost_Phi[igen]);
	  //cout<<"event: "<<jentry<<"   "<<"Reco : "<<ireco<<"   "<<"Gen : "<<igen<<" Delta R : "<<dR<<endl;

	  if(dR < 0.1  ){
	    if (dR < dR_comp_post)
	    {
	      dR_comp_post = dR;
	      matchGen2 = igen ; matchReco2 = ireco ; //iter[ireco]=igen;
	      recomatchedGen.push_back(igen);
	      Recomatchedgen.push_back(ireco);
	    }

	    genPost_Pt.push_back(newgPost_Pt[matchGen2]); genPost_Eta.push_back(newgPost_Eta[matchGen2]); genPost_Enr.push_back(newgPost_Enr[matchGen2]); genPost_Phi.push_back(newgPost_Phi[matchGen2]);
	    recoPost_Pt.push_back(newelePt[matchReco2]); recoPost_Eta.push_back(neweleEta[matchReco2]); recoPost_Enr.push_back(neweleEnr[matchReco2]); recoPost_Phi.push_back(newelePhi[matchReco2]);
	  }
	}
      }
    }

    if(recoPost_Pt.size()>=2)
    {
      reco1_post.SetPtEtaPhiE(recoPost_Pt[0],recoPost_Eta[0],recoPost_Phi[0],recoPost_Enr[0]);
      reco2_post.SetPtEtaPhiE(recoPost_Pt[1],recoPost_Eta[1],recoPost_Phi[1],recoPost_Enr[1]);
      direco_post=reco1_post+reco2_post;
      ZMass_R_post = direco_post.M();
    }

    if(genPost_Pt.size()>=2)
    {
      gen1_post.SetPtEtaPhiE(genPost_Pt[0],genPost_Eta[0],genPost_Phi[0],genPost_Enr[0]);
      gen2_post.SetPtEtaPhiE(genPost_Pt[1],genPost_Eta[1],genPost_Phi[1],genPost_Enr[1]);
      digen_post=gen1_post+gen2_post;
      ZMass_G_post = digen_post.M();
    }

    //genPost_Mass->Fill(ZMass_G_post,theWeight*puWeights);
    //recoPost_Mass->Fill(ZMass_R_post,theWeight*puWeights);
    hRM_post_RG->Fill(ZMass_R_post,ZMass_G_post,theWeight*puWeights);
    hRM_post_GR->Fill(ZMass_G_post,ZMass_R_post,theWeight*puWeights);
    
    hRM_noWeight_post_RG->Fill(ZMass_R_post,ZMass_G_post);
    hRM_noWeight_post_GR->Fill(ZMass_G_post,ZMass_R_post);

  } // event

  //hRM->Draw("colz");
  //hRM->SaveAs("Response.C");
  //hRM->SaveAs("Response.png");

  cout<<"n_ele: "<<n_ele<<"   "<<"n_trig: "<<n_trig<<"   "<<"n_kin: "<<n_kin<<"   "<<"n_id: "<<n_id<<"   "<<"n_charge: "<<n_charge<<"   "<<"n_all: "<<n_all<<endl;
  cout<<"n_all_cuts: "<<n_all_cuts<<"   "<<"nZ: "<<nZ<<endl;
  printf ("sum_weights: %f \n", sum_weights);

  //MC_nPV->Scale(1/MC_nPV->Integral());
  file->Write();
  file->Close();
}
