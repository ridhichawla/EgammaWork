#define jetRate_Est_cxx
#include "jetRate_Est.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void jetRate_Est::Loop()
{
  TFile *file = new TFile("single_photon.root", "recreate");

  int isNum, isDen;

  vector <double> newscEt;
  vector <double> newscEta;
  vector <double> newscPhi;
  vector <double> newscEnr;

  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> neweleEnr;
  vector <double> newelePhi;
  vector <double> neweleCharge;

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  TH1F *numPt = new TH1F("numPt", "numPt", 100, 0, 500);
  TH1F *denPt = new TH1F("denPt", "denPt", 100, 0, 500);

  TH1F *eleEta = new TH1F("eleEta", "eleEta", 50, -2.5, 2.5);
  TH1F *elePhi = new TH1F("elePhi", "elePhi", 50, -3.5, 3.5);

  numPt->Sumw2(); denPt->Sumw2();
  eleEta->Sumw2(); elePhi->Sumw2();

  isNum = 0;
  isDen = 0;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 500;
  cout<<"entries: "<<nentries<<endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    int index[pt->size()];
    float ptnew[pt->size()];

    for(unsigned int el=0; el<pt->size(); el++) {
      ptnew[el]=pt->at(el); }

    int size1 = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(size1,ptnew,index,true);

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    for(int k=0;k<nEle;k++){

      if(nEle <=1 ){
	if(metPt->at(0) < 10.){
	  if(expectedMissingInnerHits->at(index[k]) == 0){

	    isDen++;

	    newelePt.push_back(pt->at(index[k]));
	    neweleEta.push_back(eta->at(index[k]));
	    neweleEnr.push_back(energy->at(index[k]));
	    newelePhi.push_back(phi->at(index[k]));
	    neweleCharge.push_back(charge->at(index[k]));

	    denPt->Fill(pt->at(index[k]));
	  }
	}
      }
    }

    //cout<<"pt size: "<<newelePt.size()<<endl;
    //cout<<"Den: "<<isDen<<endl;
    //isDen = isDen + newelePt.size();
    
    for(unsigned int m=0;m<newelePt.size();m++){

      if(passMediumId->at(index[m]) == 1 && eleEcalDrivenSeed->at(index[m]) == 1){
	if(newelePt.at(m) > 20.){

	  //cout<<"Num: "<<isNum<<endl;
	  isNum++;
	  numPt->Fill(newelePt.at(m));
	}
      }
    }

  } // event

  cout<<"Numerator: "<<isNum<<"   "<<"Denominator: "<<isDen<<endl;

  file->Write();
  file->Close();
}
