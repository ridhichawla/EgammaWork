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

  TFile *file = new TFile("dyee_M50.root", "recreate");

  int good_elec, n_all_cuts;
  double sum_weights;
  double Z_Mass, Z_Rap, Z_Eta, Z_Pt, Z_Phi;
  TLorentzVector ele1,ele2,dielectron;

  vector <double> newscEt; vector <double> newscEta; vector <double> newscPhi; vector <double> newscEnr;
  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

  TH1F *MC_nPV     = new TH1F("MC_nPV", "MC_nPV", 50., 0., 50.);
  TH1F *MC_nPV_wts = new TH1F("MC_nPV_wts", "MC_nPV_wts", 50., 0., 50.);
  MC_nPV->Sumw2(); MC_nPV_wts->Sumw2();

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
  TH1D *ZPt = new TH1D("ZPt", "ZPt", 100, 0, 700);
  TH1D *ZRapidity = new TH1D("ZRapidity", "ZRapidity", 60, -3, 3);
  TH1D *ZEta = new TH1D("ZEta", "ZEta", 60, -8, 8);
  TH1D *ZPhi = new TH1D("ZPhi", "ZPhi", 50, -3.5, 3.5);

  ZMass->Sumw2(); ZPt->Sumw2(); ZRapidity->Sumw2(); ZEta->Sumw2(); ZPhi->Sumw2();

  if (fChain == 0) return;

  //Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries = 5000;
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

    int index[ptElec->size()];
    float ptElecnew[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++)
    {
      ptElecnew[el]=ptElec->at(el);
    }

    int sizer = sizeof(ptElecnew)/sizeof(ptElecnew[0]);
    TMath::Sort(sizer,ptElecnew,index,true);

    Z_Mass=0.; Z_Rap=0.; Z_Eta=0.; Z_Pt=0.; Z_Phi=0.;
    good_elec = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    int bin = 0;
    double puWeights = 1.0;

    sum_weights = sum_weights+theWeight;

    //double weight = DATA_puDist->GetBinContent(DATA_puDist->FindBin(nPUTrue))/MC_puDist->GetBinContent(MC_puDist->FindBin(nPUTrue));
    //double weight = DATA_puDist->GetBinContent(nPUTrue+1)/MC_puDist->GetBinContent(nPUTrue+1);

    bin = weights->GetXaxis()->FindBin(nPUTrue);
    puWeights = weights->GetBinContent(bin);

    //if(tauFlag){

    MC_nPV->Fill(nPV,theWeight);
    MC_nPV_wts->Fill(nPV,puWeights*theWeight);

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

	    newelePt.push_back(ptElec->at(index[k]));
	    neweleEta.push_back(etaElec->at(index[k]));
	    neweleEnr.push_back(energyElec->at(index[k]));
	    newelePhi.push_back(phiElec->at(index[k]));
	    neweleCharge.push_back(chargeElec->at(index[k]));

	  } // eta
	} // ID
      } // nEle

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

	  ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	  dielectron=ele1+ele2;
	  Z_Mass = dielectron.M();

	  Z_Pt = dielectron.Pt();
	  Z_Rap = dielectron.Rapidity();
	  Z_Eta = dielectron.Eta();
	  Z_Phi = dielectron.Phi();

	  ZMass->Fill(Z_Mass,theWeight);
	  ZPt->Fill(Z_Pt,theWeight);
	  ZRapidity->Fill(Z_Rap,theWeight);
	  ZEta->Fill(Z_Eta,theWeight);
	  ZPhi->Fill(Z_Phi,theWeight);

	} // pt
      } // good electrons
    } // nEle>=2
    //}  // tauFlag
  } // event

  cout<<"all cuts: "<<n_all_cuts<<endl;

  printf ("sum_weights: %f \n", sum_weights);

  file->Write();
  file->Close();
}
