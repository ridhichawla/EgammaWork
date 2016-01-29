#define mc13TeV_cxx
#include "mc13TeV.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void mc13TeV::Loop()
{

  TFile *f1 = TFile::Open("../../PileUp/dataPUDist.root");
  TFile *f2 = TFile::Open("../../PileUp/PileUp_MC.root");

  // data histogram 
  TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
  DATA_puDist->Scale(1/DATA_puDist->Integral());

  // mc histogram 
  TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");

  TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
  weights->Divide(MC_puDist);

  TFile *file = new TFile("dyee_M-10to50.root", "recreate");

  int good_elec, n_all_cuts;
  double sum_weights;
  double Z_Mass, Z_Y, Z_Rap, Z_Eta, Z_Pt, Z_Phi;
  TLorentzVector ele1,ele2,dielectron;

  vector <double> newscEt; vector <double> newscEta; vector <double> newscPhi; vector <double> newscEnr;
  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

  TH1F *MC_nPV = new TH1F("MC_nPV", "MC_nPV", 50, 0, 50);
  MC_nPV->Sumw2();

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

  TH1D *ZMass_ = new TH1D("ZMass_", "ZMass_", 45, xbins);
  TH1D *ZMass_0to500 = new TH1D("ZMass_0to500", "ZMass_0to500", 100, 0, 500);
  TH1D *ZMass_60to120 = new TH1D("ZMass_60to120", "ZMass_60to120", 60, 60, 120);
  TH1D *ZPt_ = new TH1D("ZPt_", "ZPt_", 100, 0, 700);
  TH1D *ZY_ = new TH1D("ZY_", "ZY_", 150, -15, 15);
  TH1D *ZRap_ = new TH1D("ZRap_", "ZRap_", 60, -3, 3);
  TH1D *ZEta_ = new TH1D("ZEta_", "ZEta_", 60, -8, 8);
  TH1D *ZPhi_ = new TH1D("ZPhi_", "ZPhi_", 50, -3.5, 3.5);

  ZMass_->Sumw2(); ZMass_0to500->Sumw2(); ZMass_60to120->Sumw2(); ZPt_->Sumw2(); ZY_->Sumw2(); ZRap_->Sumw2(); ZEta_->Sumw2(); ZPhi_->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 100;
  cout<<"entries: "<<nentries<<endl;

  sum_weights = 0.0;
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

    Z_Mass=0.; Z_Y=0.; Z_Rap=0.; Z_Eta=0.; Z_Pt=0.; Z_Phi=0.;
    good_elec = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    int bin = 0;
    double puWeights = 0.0;

    //sum_weights = sum_weights+theWeight;

    bin = weights->GetXaxis()->FindBin(nPUTrue);
    puWeights = weights->GetBinContent(bin);
    //cout<<"puWeights: "<<puWeights<<endl;

    MC_nPV->Fill(nPV,puWeights*theWeight);

    //cout<<"1"<<endl;
/*
    if(!tauFlag){
      if(single_Ele27){

	for(int k=0;k<nEle;k++){

	  if(passMediumId->at(index[k]) == 1){
	    if(fabs(eta->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){

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
      } // trigger

      //cout<<"2"<<endl;

      if(good_elec==2){

	if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){
	  ++n_all_cuts;

	  for(unsigned int i=0; i<newelePt.size(); i++)
	  {
	    elePt->Fill(newelePt[i],theWeight);
	    eleEta->Fill(neweleEta[i],theWeight);
	    scEt->Fill(newscEt[i],theWeight);
	    scEta->Fill(newscEta[i],theWeight);
	  }

	  elePt_lead->Fill(newelePt.at(0),theWeight);
	  eleEta_lead->Fill(neweleEta.at(0),theWeight);
	  elePt_slead->Fill(newelePt.at(1),theWeight);
	  eleEta_slead->Fill(neweleEta.at(1),theWeight);

	  scEt_lead->Fill(newscEt.at(0),theWeight);
	  scEt_slead->Fill(newscEt.at(1),theWeight);
	  scEta_lead->Fill(newscEta.at(0),theWeight);
	  scEta_slead->Fill(newscEta.at(1),theWeight);

	  ele1.SetPtEtaPhiE(newscEt.at(0),newscEta.at(0),newscPhi.at(0),newscEnr.at(0));
	  ele2.SetPtEtaPhiE(newscEt.at(1),newscEta.at(1),newscPhi.at(1),newscEnr.at(1));

	  dielectron=ele1+ele2;
	  Z_Mass = dielectron.M();

	  Z_Pt = dielectron.Pt();
	  Z_Y = dielectron.Y();
	  Z_Rap = dielectron.Rapidity();
	  Z_Eta = dielectron.Eta();
	  Z_Phi = dielectron.Phi();

	  ZMass_0to500->Fill(Z_Mass,theWeight);
	  ZMass_60to120->Fill(Z_Mass,theWeight);
	  ZMass_->Fill(Z_Mass,theWeight);
	  ZPt_->Fill(Z_Pt,theWeight);
	  ZY_->Fill(Z_Y,theWeight);
	  ZRap_->Fill(Z_Rap,theWeight);
	  ZEta_->Fill(Z_Eta,theWeight);
	  ZPhi_->Fill(Z_Phi,theWeight);

	} // pt
      } // good electrons

    }*/ // tauFlag
  } // event

  //cout<<"all cuts: "<<n_all_cuts<<endl;
  //cout<<" Events with positive weights: "<<evt_wgts<<endl;

  //printf ("sum_weights: %f \n", sum_weights);

  //MC_nPV->Scale(1/MC_nPV->Integral());
  file->Write();
  file->Close();
}
