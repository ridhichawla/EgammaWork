#define data13TeV_cxx
#include "data13TeV.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
//#include "DataFormats/Math/interface/deltaR.h"

void data13TeV::Loop()
{
  TFile *file = new TFile("single_electron.root", "recreate");

  int good_elec, n_all_cuts;
  double Z_Mass, Z_Rap, Z_Eta, Z_Pt, Z_Phi;
  TLorentzVector ele1,ele2,dielectron;

  vector <double> newscEt; vector <double> newscEta; vector <double> newscPhi; vector <double> newscEnr;
  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

  TH1F *Data_nPV = new TH1F("Data_nPV", "Data_nPV", 50, 0, 50);
  Data_nPV->Sumw2();

  TH1F *elePt  = new TH1F("elePt", "elePt", 100, 0, 700);
  TH1F *eleEta = new TH1F("eleEta", "eleEta", 50, -2.5, 2.5);
  TH1F *scEt   = new TH1F("scEt", "scEt", 100, 0 ,700);
  TH1F *scEta  = new TH1F("scEta", "scEta", 50, -2.5, 2.5);

  TH1F *elePt_lead   = new TH1F("elePt_lead", "elePt_lead", 100, 0, 700);
  TH1F *eleEta_lead  = new TH1F("eleEta_lead", "eleEta_lead", 50, -2.5, 2.5);
  TH1F *elePt_slead  = new TH1F("elePt_slead", "elePt_slead", 100, 0, 700);
  TH1F *eleEta_slead = new TH1F("eleEta_slead", "eleEta_slead", 50, -2.5, 2.5);

  TH1F *scEt_lead   = new TH1F("scEt_lead", "scEt_lead", 100, 0 ,700);
  TH1F *scEta_lead  = new TH1F("scEta_lead", "scEta_lead", 50, -2.5, 2.5);
  TH1F *scEt_slead  = new TH1F("scEt_slead", "scEt_slead", 100, 0 ,700);
  TH1F *scEta_slead = new TH1F("scEta_slead", "scEta_slead", 50, -2.5, 2.5);

  elePt->Sumw2(); eleEta->Sumw2(); scEt->Sumw2(); scEta->Sumw2();
  elePt_lead->Sumw2(); eleEta_lead->Sumw2(); elePt_slead->Sumw2(); eleEta_slead->Sumw2();
  scEt_lead->Sumw2(); scEta_lead->Sumw2(); scEt_slead->Sumw2(); scEta_slead->Sumw2();

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  TH1D *ZMass = new TH1D("ZMass", "ZMass", 45, xbins);
  TH1D *ZMass_0to500 = new TH1D("ZMass_0to500", "ZMass_0to500", 100, 0, 500);
  TH1D *ZMass_60to120 = new TH1D("ZMass_60to120", "ZMass_60to120", 60, 60, 120);
  TH1D *ZPt = new TH1D("ZPt", "ZPt", 100, 0, 700);
  TH1D *ZRapidity = new TH1D("ZRapidity", "ZRapidity", 60, -3, 3);
  TH1D *ZEta = new TH1D("ZEta", "ZEta", 60, -8, 8);
  TH1D *ZPhi = new TH1D("ZPhi", "ZPhi", 50, -3.5, 3.5);

  ZMass->Sumw2(); ZMass_0to500->Sumw2(); ZMass_60to120->Sumw2(); ZPt->Sumw2(); ZRapidity->Sumw2(); ZEta->Sumw2(); ZPhi->Sumw2();
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 500;
  cout<<"entries: "<<nentries<<endl;

  n_all_cuts = 0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    int index[pt->size()];
    float ptnew[pt->size()];

    for(unsigned int el=0; el<pt->size(); el++)
    {
      ptnew[el]=pt->at(el);
    }

    int sizer = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(sizer,ptnew,index,true);

    Z_Mass=0.; Z_Rap=0.; Z_Eta=0.; Z_Pt=0.; Z_Phi=0.;
    good_elec = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    Data_nPV->Fill(nPV);

    if(!single_Ele23) continue;
    if(nEle>=2) {

      for(int k=0;k<nEle;k++){

	if(passMediumId->at(index[k]) == 1){
	  if(fabs(etaSC->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){

	    if(passMediumId->at(index[k]) == 0) cout<<"Wrong ID: "<<endl;

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

	  } // eta
	} // ID
      } // nEle

      if(good_elec==2){

	if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){
	  ++n_all_cuts;

	  for(unsigned int i=0; i<newelePt.size(); i++)
	  {
	    elePt->Fill(newelePt[i]);
	    eleEta->Fill(neweleEta[i]);
	    scEt->Fill(newscEt[i]);
	    scEta->Fill(newscEta[i]);
	  }

	  elePt_lead->Fill(newelePt.at(0));
	  eleEta_lead->Fill(neweleEta.at(0));
	  elePt_slead->Fill(newelePt.at(1));
	  eleEta_slead->Fill(neweleEta.at(1));

	  scEt_lead->Fill(newscEt.at(0));
	  scEt_slead->Fill(newscEt.at(1));
	  scEta_lead->Fill(newscEta.at(0));
	  scEta_slead->Fill(newscEta.at(1));

	  ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	  dielectron=ele1+ele2;
	  Z_Mass = dielectron.M();

	  Z_Pt = dielectron.Pt();
	  Z_Rap = dielectron.Rapidity();
	  Z_Eta = dielectron.Eta();
	  Z_Phi = dielectron.Phi();

	  ZMass_0to500->Fill(Z_Mass);
	  ZMass_60to120->Fill(Z_Mass);
	  ZMass->Fill(Z_Mass);
	  ZPt->Fill(Z_Pt);
	  ZRapidity->Fill(Z_Rap);
	  ZEta->Fill(Z_Eta);
	  ZPhi->Fill(Z_Phi);
	} // pt
      } // good electrons
    } // nEle>=2
  } // event

  cout<<"n_all_cuts: "<<n_all_cuts<<endl;

  file->Write();
  file->Close();
}
