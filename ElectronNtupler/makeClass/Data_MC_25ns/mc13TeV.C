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
  /*TFile *f1 = TFile::Open("/afs/cern.ch/work/r/rchawla/private/Analysis_ee_13TeV/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/PileUp/dyM50_pileup.root");
  //TFile *f2 = TFile::Open("/afs/cern.ch/work/r/rchawla/private/CMSSW_5_3_22/src/Analysis/El_analyzer/Macros/New/MCEE/WTSfiles/MCDist_M50_13TeV.root");

  // data histogram 
  TH1F *DATA_puDist = (TH1F*)f1->Get("DATA_nPV");
  //DATA_puDist->Scale(1/DATA_puDist->Integral());

  // mc histogram 
  TH1F *MC_puDist = (TH1F*)f1->Get("MC_nPV");
  //MC_puDist->Scale(1/MC_puDist->Integral());

  TH1F *weights_ = (TH1F*)DATA_puDist->Clone("weights_");
  weights_->Divide(MC_puDist);*/

  TFile *file = new TFile("wjets.root", "recreate");
  //std::ofstream wgts;
  //wgts.open ("Weights.txt");

  int good_elec, gen_elec_pre, gen_elec_post;
  int n_notaus, n_ele, n_trig, n_kin, n_id, n_all_cuts;
  int kin, id;
  double sum_weights;
  double Z_Mass, Z_Y, Z_Rap, Z_Eta, Z_Pt, Z_Phi;
  TLorentzVector ele1,ele2,dielectron;

  vector <double> newscEt; vector <double> newscEta; vector <double> newscPhi; vector <double> newscEnr;
  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

  //TH1F *MC_nPV = new TH1F("MC_nPV", "MC_nPV", 50, 0, 50);
  //MC_nPV->Sumw2();

  //TH1F *MC_puDistW = new TH1F("MC_puDistW","MC_puDistW",50,0,50);
  //MC_puDistW->Sumw2();

  TH1D *weights = new TH1D("weights", "weights", 60000, -300000, 300000);
  
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
  TH1D *ZMass_0to400 = new TH1D("ZMass_0to400", "ZMass_0to400", 100, 0, 400);
  TH1D *ZMass_60to120 = new TH1D("ZMass_60to120", "ZMass_60to120", 60, 60, 120);
  TH1D *ZPt_ = new TH1D("ZPt_", "ZPt_", 100, 0, 700);
  TH1D *ZY_ = new TH1D("ZY_", "ZY_", 150, -15, 15);
  TH1D *ZRap_ = new TH1D("ZRap_", "ZRap_", 60, -3, 3);
  TH1D *ZEta_ = new TH1D("ZEta_", "ZEta_", 60, -8, 8);
  TH1D *ZPhi_ = new TH1D("ZPhi_", "ZPhi_", 50, -3.5, 3.5);

  ZMass_->Sumw2(); ZMass_0to400->Sumw2(); ZMass_60to120->Sumw2(); ZPt_->Sumw2(); ZY_->Sumw2(); ZRap_->Sumw2(); ZEta_->Sumw2(); ZPhi_->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 50;
  cout<<"entries: "<<nentries<<endl;

  n_notaus = 0; n_ele = 0; n_trig = 0; n_kin = 0; n_id = 0; n_all_cuts = 0;
  sum_weights = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    kin = 0; id = 0; //isCharge = false; isMass = false;

    int index[pt->size()];
    float ptnew[pt->size()];

    for(unsigned int r=0; r<pt->size(); r++)
    {
      ptnew[r]=pt->at(r);
    }

    int sizer = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(sizer,ptnew,index,true);

    Z_Mass=0.; Z_Y=0.; Z_Rap=0.; Z_Eta=0.; Z_Pt=0.; Z_Phi=0.;
    good_elec = 0;

    //int bin = 0;
    //double puWeights = 0.0;

    /*bin = weights_->GetXaxis()->FindBin(nPV);
      puWeights = weights_->GetBinContent(bin);*/
    //cout<<"puWeights: "<<puWeights<<endl;
    //MC_puDistW->Fill(nPV,puWeights);

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    //MC_nPV->Fill(nPV,puWeights);
    //cout<<"NLO Weights: "<<theWeight<<endl;
    sum_weights = sum_weights+theWeight;
    weights->Fill(theWeight);

    //cout<<"1"<<endl;
/*
    //if(!tauFlag){
      ++n_notaus;
      if(nEle>=2) {
	++n_ele;
	if(doubleElectron) {
	  ++n_trig;
	  if(doubleElectron == 0) cout<<"Wrong"<<endl;

	  for(int k=0;k<nEle;k++){

	    if(passMediumId->at(index[k]) == 1 && eleEcalDrivenSeed->at(index[k]) == 1){
	      id++;
	      if(fabs(eta->at(index[k])) < 2.5 && pt->at(index[k]) > 10. && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){
		kin++;

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
	      }
	    }
	  }
	}
      }
    //}

    //cout<<"2"<<endl;

    if(good_elec>=2){
      if(newelePt.at(0) < newelePt.at(1)) cout<<"event: "<<jentry<<"   "<<"Sorting not proper: "<<"   "<<"reco pt lead: "<<newelePt.at(0)<<"   "<<"reco pt sublead: "<<newelePt.at(1)<<endl;

      if(newelePt.at(0) > 20. && newelePt.at(1) > 10.){
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

	ele1.SetPtEtaPhiE(newscEt.at(0),newscEta.at(0),newscPhi.at(0),newscEnr.at(0));
	ele2.SetPtEtaPhiE(newscEt.at(1),newscEta.at(1),newscPhi.at(1),newscEnr.at(1));

	dielectron=ele1+ele2;
	Z_Mass = dielectron.M();

	//if(Z_Mass >= 60 && Z_Mass <= 120){ isMass = true; ++nZ;}

	Z_Pt = dielectron.Pt();
	Z_Y = dielectron.Y();
	Z_Rap = dielectron.Rapidity();
	Z_Eta = dielectron.Eta();
	Z_Phi = dielectron.Phi();

	ZMass_0to400->Fill(Z_Mass);
	ZMass_60to120->Fill(Z_Mass);
	ZMass_->Fill(Z_Mass);
	ZPt_->Fill(Z_Pt);
	ZY_->Fill(Z_Y);
	ZRap_->Fill(Z_Rap);
	ZEta_->Fill(Z_Eta);
	ZPhi_->Fill(Z_Phi);
      }
    }

    if(id>=2){
      n_id++;
      if(kin>=2){
	n_kin++;
      } // id
    } // kin
*/
  } // event

  //cout<<"no taus: "<<n_notaus<<"   "<<"ele: "<<n_ele<<"   "<<"trigger: "<<n_trig<<"   "<<"ID: "<<n_id<<"   "<<"kin: "<<n_kin<<"   "<<"all: "<<n_all_cuts<<endl;
  printf ("sum_weights: %f \n", sum_weights);

  //MC_nPV->Scale(1/MC_nPV->Integral());
  file->Write();
  file->Close();
}
