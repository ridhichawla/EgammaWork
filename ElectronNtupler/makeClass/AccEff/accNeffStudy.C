#define accNeffStudy_cxx
#include "accNeffStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void accNeffStudy::Loop()
{

  TFile *file = new TFile("singleElectron.root", "recreate");

  double sum_weights;
  double num_Mass, den_Mass;
  TLorentzVector numEle1, numEle2, numDielectron, denEle1, denEle2, denDielectron;

  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleMediumId;

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  TH1D *numMass = new TH1D("numMass", "numMass", 45, xbins);
  TH1D *denMass = new TH1D("denMass", "denMass", 45, xbins);

  numMass->Sumw2(); denMass->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 10;
  cout<<"entries: "<<nentries<<endl;

  sum_weights = 0.0;

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

    num_Mass=0.; den_Mass=0.;
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleMediumId.clear();

    sum_weights = sum_weights+theWeight;

    //if(!tauFlag){

    for(int k=0;k<nEle;k++){
      if(fabs(etaSC->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){

	newelePt.push_back(pt->at(index[k]));
	neweleEta.push_back(eta->at(index[k]));
	neweleEnr.push_back(energy->at(index[k]));
	newelePhi.push_back(phi->at(index[k]));
	neweleMediumId.push_back(passMediumId->at(index[k]));

      } // eta
    } //nEle

    //cout<<"2"<<endl;

    if(newelePt.size()>=2){

      if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){

	denEle1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	denEle2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	denDielectron=denEle1+denEle2;
	den_Mass = denDielectron.M();

	denMass->Fill(den_Mass);

	if(neweleMediumId.at(0)==1 && neweleMediumId.at(1)==1){
	  numEle1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  numEle2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	  numDielectron=numEle1+numEle2;
	  num_Mass = numDielectron.M();

	  numMass->Fill(num_Mass);
	}

      } // pt
    } // good electrons
    //}  // tauFlag
  } // event

  printf ("sum_weights: %f \n", sum_weights);

  file->Write();
  file->Close();
}
