#define dataWZGroup_cxx
#include "dataWZGroup.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
//#include "DataFormats/Math/interface/deltaR.h"

void dataWZGroup::Loop()
{
  TFile *file = new TFile("WZGroup.root", "recreate");
  //TFile *file = new TFile("doubleEG.root", "recreate");

  //double Mz = 91.1876;

  int good_elec;
  int n_ele, n_trig, n_kin, n_id, n_charge, n_all;
  bool isKin, isID;
  int n_dR;
  double ZMass, ZY, ZRap, ZEta, ZPt, ZPhi;
  //double deltaR1, deltaR2;
  TLorentzVector ele1,ele2,dielectron,ele12,ele22,dielectron2;

  vector <double> newscEt;
  vector <double> newscEta;
  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> neweleEnr;
  vector <double> newelePhi;
  vector <double> neweleCharge;
  vector <double> neweleIso;

  TH1F *eta_sc = new TH1F("eta_sc", "eta_sc", 50, -2.5, 2.5);
  TH1F *et_sc = new TH1F("et_sc", "et_sc", 100, 0 ,700);
  TH1F *et_sc_rebin = new TH1F("et_sc_rebin", "et_sc_rebin", 700, 0 ,700);

  TH1F *pt_lead = new TH1F("pt_lead", "pt_lead", 100, 0, 700);
  TH1F *eta_lead = new TH1F("eta_lead", "eta_lead", 50, -2.5, 2.5);
  TH1F *pt_slead = new TH1F("pt_slead", "pt_slead", 100, 0, 700);
  TH1F *eta_slead = new TH1F("eta_slead", "eta_slead", 48, -2.4, 2.4);

  TH1F *et_sc_lead = new TH1F("et_sc_lead", "et_sc_lead", 105, 0 ,700);
  TH1F *et_sc_slead = new TH1F("et_sc_slead", "et_sc_slead", 105, 0 ,700);
  TH1F *eta_sc_lead = new TH1F("eta_sc_lead", "eta_sc_lead", 48, -2.4, 2.4);
  TH1F *eta_sc_slead = new TH1F("eta_sc_slead", "eta_sc_slead", 48, -2.4, 2.4);

  /*TH1F *B_delta_eta = new TH1F("B_delta_eta", "B_delta_eta", 260, -0.01, 0.01);
    TH1F *B_delta_phi = new TH1F("B_delta_phi", "B_delta_phi", 160, -0.08, 0.08);
    TH1F *B_H_E = new TH1F("B_H_E", "B_H_E", 100, 0, 0.1);
    TH1F *B_sigma_eta = new TH1F("B_sigma_eta", "B_sigma_eta", 100, 0, 0.1);
    TH1F *B_E_P = new TH1F("B_E_P", "B_E_P", 400, 0, 0.5);

    TH1F *E_delta_eta = new TH1F("E_delta_eta", "E_delta_eta", 260, -0.01, 0.01);
    TH1F *E_delta_phi = new TH1F("E_delta_phi", "E_delta_phi", 160, -0.08, 0.08);
    TH1F *E_H_E = new TH1F("E_H_E", "E_H_E", 100, 0, 0.1);
    TH1F *E_sigma_eta = new TH1F("E_sigma_eta", "E_sigma_eta", 100, 0, 0.1);
    TH1F *E_E_P = new TH1F("E_E_P", "E_E_P", 400, 0, 0.5);*/

  TH1F *Data_nPV = new TH1F("Data_nPV", "Data_nPV", 50, 0, 50);
  TH1F *goodeleNum = new TH1F("goodeleNum", "goodeleNum", 7, 0, 7);
  
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

  //Double_t xbins[7] = {15.,30.,60.,120.,240.,600.,2000.};

  //TH1D *z_mass_gen = new TH1D("z_mass_gen", "z_mass_gen", 6, xbins);
  //TH1D *z_mass_gen_cut = new TH1D("z_mass_gen_cut", "z_mass_gen_cut", 6, xbins);
  //TH1D *z_mass_gen = new TH1D("z_mass_gen", "z_mass_gen", 100, 0., 200.);
  //TH1D *z_mass_gen_cut = new TH1D("z_mass_gen_cut", "z_mass_gen_cut", 100, 0., 200.);
  TH1D *z_mass_1 = new TH1D("z_mass_1", "z_mass_1", 400, 0, 400);
  TH1D *z_mass_2 = new TH1D("z_mass_2", "z_mass_2", 80, 0, 400);
  TH1D *z_mass = new TH1D("z_mass", "z_mass", 60, 60, 120);
  TH1D *z_pt = new TH1D("z_pt", "z_pt", 100, 0, 700);
  TH1D *z_y = new TH1D("z_y", "z_y", 150, -15, 15);
  TH1D *z_rap = new TH1D("z_rap", "z_rap", 60, -3, 3);
  TH1D *z_eta = new TH1D("z_eta", "z_eta", 60, -8, 8);
  TH1D *z_phi = new TH1D("z_phi", "z_phi", 50, -3.5, 3.5);

  //TH2D *hRM = new TH2D("Response Matrix", "Response Matrix", 6, xbins, 6, xbins);
  //TH2D *hRM_match = new TH2D("Response Matrix_match", "Response Matrix_match", 6, xbins, 6, xbins);

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //int nentries = 200;
  cout<<"entries: "<<nentries<<endl;

  n_ele = 0;
  n_trig = 0;
  n_kin = 0;
  n_id = 0;
  n_charge = 0;
  n_all = 0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //cout<<"entry: "<<jentry<<endl;
    //cout<<"No. of electrons: "<<nEle<<endl;

    good_elec = 0;

    newscEt.clear();
    newscEta.clear();
    newelePt.clear();
    neweleEta.clear();
    neweleEnr.clear();
    newelePhi.clear();
    neweleCharge.clear();
    neweleIso.clear();

    Data_nPV->Fill(nPV);

    //if(nEle<2) continue;
    //++n_ele;
    
    if(!singleElectron_data) continue;
    ++n_trig;
    
    for(int k=0;k<nEle;k++){

      if(fabs(etaSC->at(k)) >= 2.5) continue;
      if(etSC->at(k) <= 25.) continue;
      if(passMediumId->at(k) == 0) continue;

      good_elec = good_elec + 1;
      newscEt.push_back(etSC->at(k));
      newscEta.push_back(etaSC->at(k));
      newelePt.push_back(pt->at(k));
      neweleEta.push_back(eta->at(k));
      neweleEnr.push_back(energy->at(k));
      newelePhi.push_back(phi->at(k));
      neweleCharge.push_back(charge->at(k));
      neweleIso.push_back(isoRho->at(k));
      //goodeleNum->Fill(nEle);
      //}

    }
    //if(isKin && isID){
    //goodeleNum->Fill(nEle);
    //++n_id;}

    if(good_elec<2) continue;
    if(good_elec>2)
    cout<<"good elec "<<good_elec<<endl;
    //cout<<"entry: "<<jentry<<endl;}

    //if(jentry==11171 || jentry==39817 || jentry==54637 || jentry==77355 || jentry==182372 || jentry==186789)
    //cout<<"ele pt[3] "<<newscEt[3]<<endl;

    //goodeleNum->Fill(good_elec);

    //if(fabs(newscEta[0]) >= 2.5 && fabs(neweleEta[1]) >= 2.5) continue;// && newscEt[0] <= 25. && newscEt[1] <= 25. && newscEt[2] <= 25.) continue;
    ++n_id;
    if(neweleCharge[0]==neweleCharge[1]) continue;
    ++n_charge;

    for(unsigned int i=0; i<newscEt.size(); i++)
    {
      eta_sc->Fill(newscEta[i]);
      et_sc->Fill(newscEt[i]);
      et_sc_rebin->Fill(newscEt[i]);
    }

    pt_lead->Fill(newelePt[0]);
    eta_lead->Fill(neweleEta[0]);
    pt_slead->Fill(newelePt[1]);
    eta_slead->Fill(neweleEta[1]);

    et_sc_lead->Fill(newscEt[0]);
    et_sc_slead->Fill(newscEt[1]);
    eta_sc_lead->Fill(newscEta[0]);
    eta_sc_slead->Fill(newscEta[1]);

    ele1.SetPtEtaPhiE(newelePt[0],neweleEta[0],newelePhi[0],neweleEnr[0]);
    ele2.SetPtEtaPhiE(newelePt[1],neweleEta[1],newelePhi[1],neweleEnr[1]);

    dielectron=ele1+ele2;
    ZMass = dielectron.M();
    ZPt = dielectron.Pt();
    ZY = dielectron.Y();
    ZRap = dielectron.Rapidity();
    ZEta = dielectron.Eta();
    ZPhi = dielectron.Phi();

    //cout<<"Event: "<<jentry<<"   Z mass gen: "<<Z_mass[0]<<"  Z mass reco: "<<ZMass<<endl;

    z_mass_1->Fill(ZMass);
    z_mass_2->Fill(ZMass);
    z_mass->Fill(ZMass);
    z_pt->Fill(ZPt);
    z_y->Fill(ZY);
    z_rap->Fill(ZRap);
    z_eta->Fill(ZEta);
    z_phi->Fill(ZPhi);
    ++n_all;

  } // event

  //cout<<"n_ele: "<<n_ele<<"   "<<"n_trig: "<<n_trig<<endl;//"   "<<"n_kin: "<<n_kin<<"   "
  cout<<"n_id: "<<n_id<<"   "<<"n_charge: "<<n_charge<<endl;//"   "<<"n_all: "<<n_all<<endl;
  cout<<"n_all: "<<n_all<<endl;
  cout<<"60-120: "<<z_mass->Integral()<<endl;
  file->Write();
  file->Close();
}
