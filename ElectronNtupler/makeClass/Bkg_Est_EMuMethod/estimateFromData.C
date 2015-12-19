#define estimateFromData_cxx
#include "estimateFromData.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
//#include "DataFormats/Math/interface/deltaR.h"

void estimateFromData::Loop()
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

  TH1F *elePt_z  = new TH1F("elePt_z", "elePt_z", 30, 0, 150);
  TH1F *elePt  = new TH1F("elePt", "elePt", 100, 0, 700);
  TH1F *eleEta = new TH1F("eleEta", "eleEta", 50, -2.5, 2.5);
  TH1F *elePhi = new TH1F("elePhi", "elePhi", 50, -3.5, 3.5);

  TH1F *muonPt_z  = new TH1F("muonPt_z", "muonPt_z", 30, 0, 150);
  TH1F *muonPt  = new TH1F("muonPt", "muonPt", 100, 0, 700);
  TH1F *muonEta = new TH1F("muonEta", "muonEta", 50, -2.5, 2.5);
  TH1F *muonPhi = new TH1F("muonPhi", "muonPhi", 50, -3.5, 3.5);

  Z_Mass->Sumw2(); EMu_Mass->Sumw2(); Wjets_Mass->Sumw2();

  elePt_z->Sumw2(); elePt->Sumw2(); eleEta->Sumw2(); elePhi->Sumw2();
  muonPt_z->Sumw2(); muonPt->Sumw2(); muonEta->Sumw2(); muonPhi->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 500;
  cout<<"entries: "<<nentries<<endl;

  sum_weights = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    int index[pt->size()];
    float ptnew[pt->size()];

    for(unsigned int el=0; el<pt->size(); el++)
    {
      ptnew[el]=pt->at(el);
    }

    int size = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(size,ptnew,index,true);

    good_elec = 0;
    good_muon = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();
    newmuonPt.clear(); newmuonEta.clear(); newmuonEnr.clear(); newmuonPhi.clear(); newmuonCharge.clear();

    sum_weights = sum_weights+theWeight;

    if(ptMuon->size() >=2){
      if(ptMuon->at(0) < ptMuon->at(1)) {cout<<"Sorting required"<<endl;}}

    if(!doubleEMu_17_8) continue;
    if(!doubleEMu_17_8) cout<<"Trigger not applied"<<endl;
    
    //if(doubleEMu_17_8){
      for(int k=0;k<nEle;k++){

	if(fabs(eta->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){

	  if(passMediumId->at(index[k]) == 1 && eleEcalDrivenSeed->at(index[k]) == 1){
	    if(passMediumId->at(index[k]) == 0) cout<<"Wrong ID: "<<endl;
	    if(eleEcalDrivenSeed->at(index[k]) == 0) cout<<"Wrong ECAL ID: "<<endl;

	    good_elec = good_elec + 1;

	    newscEt.push_back(etSC->at(index[k]));
	    newscEta.push_back(etaSC->at(index[k]));
	    newscPhi.push_back(phiSC->at(index[k]));
	    newscEnr.push_back(enSC->at(index[k]));
	    newelePt.push_back(pt->at(index[k]));
	    neweleEta.push_back(eta->at(index[k]));
	    neweleEnr.push_back(energy->at(index[k]));
	    newelePhi.push_back(phi->at(index[k]));
	    neweleCharge.push_back(charge->at(index[k]));

	  } // ID
	} // eta
      } // nEle
    //}

    //if(doubleEMu_17_8){
      for(int l=0;l<nMuons;l++){

	if(fabs(etaMuon->at(l)) < 2.4){
	  if(isTight->at(l)){
	    if(isoPFMuon->at(l) < 0.15){
	      good_muon = good_muon+1;

	      newmuonPt.push_back(ptMuon->at(l));
	      newmuonEta.push_back(etaMuon->at(l));
	      newmuonPhi.push_back(phiMuon->at(l));
	      newmuonEnr.push_back(energyMuon->at(l));
	      newmuonCharge.push_back(chargeMuon->at(l));

	    } // Isolation
	  } // ID
	} // eta
      } // nMuons
    //} // trigger

    //cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    //cout<<""<<endl;

    if(good_elec==2){
      //cout<<"event: "<<jentry<<"   "<<"ele 1: "<<newelePt.at(0)<<"   "<<"ele2: "<<newelePt.at(1)<<endl;
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

	  elePt_z->Fill(newelePt.at(0));
	  elePt->Fill(newelePt.at(0));
	  eleEta->Fill(neweleEta.at(0));
	  elePhi->Fill(newelePhi.at(0));

	  muonPt_z->Fill(newmuonPt.at(0));
	  muonPt->Fill(newmuonPt.at(0));
	  muonEta->Fill(newmuonEta.at(0));
	  muonPhi->Fill(newmuonPhi.at(0));

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


  } // event

  //cout<<"n_all_cuts: "<<n_all_cuts<<"   "<<"nZ: "<<nZ<<endl;

  //printf ("sum_weights: %f \n", sum_weights);

  file->Write();
  file->Close();
}
