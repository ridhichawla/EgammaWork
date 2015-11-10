#define mcTTBar_cxx
#include "mcTTBar.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
//#include "DataFormats/Math/interface/deltaR.h"

void mcTTBar::Loop()
{
  TFile *file = new TFile("ttbar.root", "recreate");

  int good_elec, good_muon;
  int n_all_cuts, nZ, n_dR;
  double ZMass, ZY, ZRap, ZEta, ZPt, ZPhi;
  TLorentzVector ele1,ele2,dielectron;

  vector <double> newscEt;
  vector <double> newscEta;
  vector <double> newscPhi;
  vector <double> newscEnr;
  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> neweleEnr;
  vector <double> newelePhi;
  vector <double> neweleCharge;

  vector <double> newmuonPt;
  vector <double> newmuonEta;
  vector <double> newmuonEnr;
  vector <double> newmuonPhi;
  vector <double> newmuonCharge;

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

  //const double xbins[38] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510};
  //const double xbins[19] = {45., 50., 55., 60., 64., 68., 72., 76., 81., 86., 91., 96., 101., 106., 110., 115., 120., 126., 133.};
  //const double xbins[5] = {10.,30.,60.,120.,500.};
  //Double_t xbins[7] = {15.,30.,60.,120.,240.,600.,2000.};

  Double_t xbins[42] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 1000, 1500, 2000};
  TH1D *ZMass_ = new TH1D("ZMass_", "ZMass_", 41, xbins);
  TH1D *ZMass_0to400 = new TH1D("ZMass_0to400", "ZMass_0to400", 100, 0, 400);
  TH1D *ZMass_60to120 = new TH1D("ZMass_60to120", "ZMass_60to120", 60, 60, 120);
  TH1D *ZPt_ = new TH1D("ZPt_", "ZPt_", 100, 0, 700);
  TH1D *ZY_ = new TH1D("ZY_", "ZY_", 150, -15, 15);
  TH1D *ZRap_ = new TH1D("ZRap_", "ZRap_", 60, -3, 3);
  TH1D *ZEta_ = new TH1D("ZEta_", "ZEta_", 60, -8, 8);
  TH1D *ZPhi_ = new TH1D("ZPhi_", "ZPhi_", 50, -3.5, 3.5);

  ZMass_->Sumw2();
  ZMass_0to400->Sumw2(); ZMass_60to120->Sumw2(); ZPt_->Sumw2(); ZY_->Sumw2(); ZRap_->Sumw2(); ZEta_->Sumw2(); ZPhi_->Sumw2();

  //TH2D *hRM = new TH2D("Response Matrix", "Response Matrix", 6, xbins, 6, xbins);
  //TH2D *hRM_match = new TH2D("Response Matrix_match", "Response Matrix_match", 6, xbins, 6, xbins);

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 30;
  cout<<"entries: "<<nentries<<endl;

  n_dR = 0;
  n_all_cuts = 0;
  nZ = 0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    good_elec = 0;
    good_muon = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();
    newmuonPt.clear(); newmuonEta.clear(); newmuonEnr.clear(); newmuonPhi.clear(); newmuonCharge.clear();

    Data_nPV->Fill(nPV);

    if(nEle>=1){
      for(int k=0;k<nEle;k++){

	//bool isBarrel = (fabs(etaSC->at(k)) <= 1.479);
	//bool isEndcap = (fabs(etaSC->at(k)) > 1.479 && fabs(etaSC->at(k)) < 2.5);

	if(fabs(eta->at(k)) < 2.5 && pt->at(k) > 25. && !(fabs(etaSC->at(k)) > 1.4442 && fabs(etaSC->at(k)) < 1.566)){

	  if(passMediumId->at(k) == 1 && eleEcalDrivenSeed->at(k) == 1){
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

    //cout<<"event: "<<jentry<<endl;
    if(nMuons>=1){

      for(int l=0;l<nMuons;l++){

	if(fabs(etaMuon->at(l)) < 2.4 && ptMuon->at(l) > 20.){
	  if(isTight->at(l)){
	    if(isoPFMuon->at(l) < 0.15){
	      good_muon = good_muon+1;

	      newmuonPt.push_back(ptMuon->at(l));
	      newmuonEta.push_back(etaMuon->at(l));
	      newmuonPhi.push_back(phiMuon->at(l));
	      newmuonEnr.push_back(energyMuon->at(l));
	      newmuonCharge.push_back(chargeMuon->at(l));

	    }
	  }
	}
      }
    }
    //cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //cout<<""<<endl;
    /*
       if(good_elec>=2){

       if(neweleCharge[0] * neweleCharge[1] == -1){
       ++n_all_cuts;

       for(unsigned int i=0; i<newelePt.size(); i++)
       {
       pt_->Fill(newelePt[i]);
       eta_->Fill(neweleEta[i]);
       et_sc_->Fill(newscEt[i]);
       eta_sc_->Fill(newscEta[i]);
       }

       pt_lead->Fill(newelePt[0]);
       eta_lead->Fill(neweleEta[0]);
       pt_slead->Fill(newelePt[1]);
       eta_slead->Fill(neweleEta[1]);

       et_sc_lead->Fill(newscEt[0]);
       et_sc_slead->Fill(newscEt[1]);
       eta_sc_lead->Fill(newscEta[0]);
       eta_sc_slead->Fill(newscEta[1]);

       ele1.SetPtEtaPhiE(newscEt[0],newscEta[0],newscPhi[0],newscEnr[0]);
       ele2.SetPtEtaPhiE(newscEt[1],newscEta[1],newscPhi[1],newscEnr[1]);

       dielectron=ele1+ele2;
       ZMass = dielectron.M();

       if(ZMass >= 60 && ZMass <= 120){
       ++nZ;
       }

       ZPt = dielectron.Pt();
       ZY = dielectron.Y();
       ZRap = dielectron.Rapidity();
       ZEta = dielectron.Eta();
       ZPhi = dielectron.Phi();

       ZMass_0to400->Fill(ZMass);
       ZMass_60to120->Fill(ZMass);
       ZMass_->Fill(ZMass);
       ZPt_->Fill(ZPt);
       ZY_->Fill(ZY);
       ZRap_->Fill(ZRap);
       ZEta_->Fill(ZEta);
       ZPhi_->Fill(ZPhi);
       }
       }
       */
    //cout<<"new Pt size: "<<newelePt.size()<<endl;
    /*bool leg1_e1 = false;
      bool leg1_e2 = false;
      bool leg2_e1 = false;
      bool leg2_e2 = false;

      for(unsigned int j = 0; j < pt_f1->size(); j++){
      double dR1 = deltaR(neweleEta[0],newelePhi[0],eta_f1->at(j),phi_f1->at(j));
      double dR2 = deltaR(neweleEta[1],newelePhi[1],eta_f1->at(j),phi_f1->at(j));
      if (dR1 < 0.1) leg1_e1 = true;
      if (dR2 < 0.1) leg1_e2 = true;
      }

      for(unsigned int k = 0; k < pt_f2->size(); k++){
      double dR3 = deltaR(neweleEta[0],newelePhi[0],eta_f2->at(k),phi_f2->at(k));
      double dR4 = deltaR(neweleEta[1],newelePhi[1],eta_f2->at(k),phi_f2->at(k));
      if (dR3 < 0.1) leg2_e1 = true;
      if (dR4 < 0.1) leg2_e2 = true;
      }*/

    //if((leg1_e1 && leg2_e2) || (leg1_e2 && leg2_e1)){

    //++n_dR;
    //} // dR matching

    //hRM->Fill(Mee_post,ZMass);
    //if(deltaR1 < 0.3 && deltaR2 < 0.3) hRM_match->Fill(Mee_post,ZMass);
    //hRM->Draw("colz");
    //cout<<"Event: "<<jentry<<"   Z mass gen: "<<Mee<<"  Z mass reco: "<<ZMass<<endl;

  } // event

  cout<<"n_all_cuts: "<<n_all_cuts<<"   "<<"nZ: "<<nZ<<endl;

  //Data_nPV->Scale(1/Data_nPV->Integral());
  file->Write();
  file->Close();
}
