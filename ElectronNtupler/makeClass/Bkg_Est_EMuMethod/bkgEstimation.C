#define bkgEstimation_cxx
#include "bkgEstimation.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void bkgEstimation::Loop()
{
  TFile *file = new TFile("muonEG_Run2015D.root", "recreate");

  int good_elec, good_muon;
  double sum_weights;
  double ZMass, EMuMass, Wjet_Mass;
  TLorentzVector ele1,ele2,diZMass;
  TLorentzVector emu1,emu2,diEMuMass;
  TLorentzVector wjets1,wjets2,diWjets_Mass;

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

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  TH1D *Z_Mass= new TH1D("Z_Mass", "Z_Mass", 45, xbins);
  TH1D *EMu_Mass= new TH1D("EMu_Mass", "EMu_Mass", 45, xbins);
  TH1D *Wjets_Mass= new TH1D("Wjets_Mass", "Wjets_Mass", 45, xbins);

  TH1F *elePt  = new TH1F("elePt", "elePt", 100, 0, 700);
  TH1F *eleEta = new TH1F("eleEta", "eleEta", 50, -2.5, 2.5);
  TH1F *elePhi = new TH1F("elePhi", "elePhi", 50, -3.5, 3.5);

  TH1F *muonPt  = new TH1F("muonPt", "muonPt", 100, 0, 700);
  TH1F *muonEta = new TH1F("muonEta", "muonEta", 50, -2.5, 2.5);
  TH1F *muonPhi = new TH1F("muonPhi", "muonPhi", 50, -3.5, 3.5);

  Z_Mass->Sumw2(); EMu_Mass->Sumw2(); Wjets_Mass->Sumw2();

  elePt->Sumw2(); eleEta->Sumw2(); elePhi->Sumw2();
  muonPt->Sumw2(); muonEta->Sumw2(); muonPhi->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 5000;
  cout<<"entries: "<<nentries<<endl;

  sum_weights = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    int index1[pt->size()];
    float ptnew1[pt->size()];

    for(unsigned int el=0; el<pt->size(); el++) {
      ptnew1[el]=pt->at(el); }

    int size1 = sizeof(ptnew1)/sizeof(ptnew1[0]);
    TMath::Sort(size1,ptnew1,index1,true);

    int index2[ptMuon->size()];
    float ptnew2[ptMuon->size()];

    for(unsigned int mu=0; mu<ptMuon->size(); mu++) {
      ptnew2[mu]=ptMuon->at(mu); }

    int size2 = sizeof(ptnew2)/sizeof(ptnew2[0]);
    TMath::Sort(size2,ptnew2,index2,true);

    good_elec = 0;
    good_muon = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();
    newmuonPt.clear(); newmuonEta.clear(); newmuonEnr.clear(); newmuonPhi.clear(); newmuonCharge.clear();

    //sum_weights = sum_weights+theWeight;

    //if(tauFlag){
    for(int k=0;k<nEle;k++){

      if(fabs(eta->at(index1[k])) < 2.5 && !(fabs(etaSC->at(index1[k])) > 1.4442 && fabs(etaSC->at(index1[k])) < 1.566)){

	if(passMediumId->at(index1[k]) == 1 && eleEcalDrivenSeed->at(index1[k]) == 1){
	  if(passMediumId->at(index1[k]) == 0) cout<<"Wrong ID: "<<endl;
	  if(eleEcalDrivenSeed->at(index1[k]) == 0) cout<<"Wrong ECAL ID: "<<endl;

	  good_elec = good_elec + 1;

	  newscEt.push_back(etSC->at(index1[k]));
	  newscEta.push_back(etaSC->at(index1[k]));
	  newscPhi.push_back(phiSC->at(index1[k]));
	  newscEnr.push_back(enSC->at(index1[k]));
	  newelePt.push_back(pt->at(index1[k]));
	  neweleEta.push_back(eta->at(index1[k]));
	  neweleEnr.push_back(energy->at(index1[k]));
	  newelePhi.push_back(phi->at(index1[k]));
	  neweleCharge.push_back(charge->at(index1[k]));
	} // ID
      } // eta
    } // nEle

    for(int l=0;l<nMuons;l++){

      if(fabs(etaMuon->at(index2[l])) < 2.4){
	if(isHEEP->at(index2[l])){
	  if(isoTrkMuon->at(index2[l]) < 0.1){
	    good_muon = good_muon+1;

	    newmuonPt.push_back(ptMuon->at(index2[l]));
	    newmuonEta.push_back(etaMuon->at(index2[l]));
	    newmuonPhi.push_back(phiMuon->at(index2[l]));
	    newmuonEnr.push_back(energyMuon->at(index2[l]));
	    newmuonCharge.push_back(chargeMuon->at(index2[l]));

	  } // Isolation
	} // ID
      } // eta
    } //nMuons


    if(good_elec==2){
      if(newelePt.at(0) > 25. && newelePt.at(1) > 15.){
	if(neweleCharge.at(0)*neweleCharge.at(1) == -1){

	  ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	  diZMass=ele1+ele2;
	  ZMass = diZMass.M();

	  Z_Mass->Fill(ZMass);
	  //Z_Mass->Fill(ZMass,theWeight);

	} // opposite charge
      } // pt cut
    } // only two electrons

    if(good_elec==1 && good_muon==1){
      if(newelePt.at(0) > 25. && newmuonPt.at(0) > 15.){
	if(neweleCharge.at(0)*newmuonCharge.at(0) == -1){

	  if(neweleCharge.at(0)*newmuonCharge.at(0) ==1) cout<<"Opposite charge condition is not satisfied"<<endl;

	  elePt->Fill(newelePt.at(0));
	  eleEta->Fill(neweleEta.at(0));
	  elePhi->Fill(newelePhi.at(0));

	  muonPt->Fill(newmuonPt.at(0));
	  muonEta->Fill(newmuonEta.at(0));
	  muonPhi->Fill(newmuonPhi.at(0));

	  /*elePt->Fill(newelePt.at(0),theWeight);
	    eleEta->Fill(neweleEta.at(0),theWeight);
	    elePhi->Fill(newelePhi.at(0),theWeight);

	    muonPt->Fill(newmuonPt.at(0),theWeight);
	    muonEta->Fill(newmuonEta.at(0),theWeight);
	    muonPhi->Fill(newmuonPhi.at(0),theWeight);*/

	  emu1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  emu2.SetPtEtaPhiE(newmuonPt.at(0),newmuonEta.at(0),newmuonPhi.at(0),newmuonEnr.at(0));

	  diEMuMass=emu1+emu2;
	  EMuMass = diEMuMass.M();

	  EMu_Mass->Fill(EMuMass);
	  //EMu_Mass->Fill(EMuMass,theWeight);

	} // opposite charge
      } // pt cut
    } // one ele and one muon

    if(good_elec==1 && good_muon==1){
      if(newelePt.at(0) > 25. && newmuonPt.at(0) > 15.){
	if(neweleCharge.at(0)*neweleCharge.at(0) == 1){
	  wjets1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  wjets2.SetPtEtaPhiE(newmuonPt.at(0),newmuonEta.at(0),newmuonPhi.at(0),newmuonEnr.at(0));

	  diWjets_Mass=wjets1+wjets2;
	  Wjet_Mass = diWjets_Mass.M();

	  Wjets_Mass->Fill(Wjet_Mass);
	  //Wjets_Mass->Fill(Wjet_Mass,theWeight);

	} // same charge
      } // pt cut
    } // one ele ans one muon

    //} // tauFlag

  } // event

  //cout<<"n_all_cuts: "<<n_all_cuts<<"   "<<"nZ: "<<nZ<<endl;

  //printf ("sum_weights: %f \n", sum_weights);

  file->Write();
  file->Close();
}
