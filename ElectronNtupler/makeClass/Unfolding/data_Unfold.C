#define data_Unfold_cxx
#include "data_Unfold.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void data_Unfold::Loop()
{

  TFile *file = new TFile("single_electron.root", "recreate"); //dyJtoLL_M50.root", "recreate");
  
  int good_elec;
  double Z_Mass;
  
  TLorentzVector ele1,ele2,dielectron;

  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

  Double_t xbins[48] = {0,10,15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  TH1D *reco_ZMass = new TH1D("ZMass", "ZMass", 47, xbins);
  reco_ZMass->Sumw2();

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

    for(unsigned int r=0; r<pt->size(); r++)
    {
      ptnew[r]=pt->at(r);
    }

    int sizer = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(sizer,ptnew,index,true);
    
    //cout<<""<<endl;
    Z_Mass=0.;
    good_elec = 0;

    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    //cout<<"1"<<endl;

    if(nEle>=2) {
      if(single_Ele23) {

	for(int k=0;k<nEle;k++){

	  if(passMediumId->at(index[k]) == 1){
	    if(fabs(eta->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){

	      if(passMediumId->at(index[k]) == 0) cout<<"Wrong ID: "<<endl;

	      good_elec = good_elec + 1;

	      newelePt.push_back(pt->at(index[k]));
	      neweleEta.push_back(eta->at(index[k]));
	      neweleEnr.push_back(energy->at(index[k]));
	      newelePhi.push_back(phi->at(index[k]));
	      neweleCharge.push_back(charge->at(index[k]));

	    } // kin cuts
	  } // ID
	} // nEle
      } // trigger
    } // nEle >=2

    //cout<<"2"<<endl;

    if(good_elec==2){
      if(newelePt.at(0) < newelePt.at(1)) cout<<"event: "<<jentry<<"   "<<"Sorting not proper: "<<"   "<<"reco pt lead: "<<newelePt.at(0)<<"   "<<"reco pt sublead: "<<newelePt.at(1)<<endl;

      if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){

	ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	dielectron=ele1+ele2;
	Z_Mass = dielectron.M();
	reco_ZMass->Fill(Z_Mass);
      }
    }

    //cout<<"3"<<endl;

  } // event

  file->Write();
  file->Close();
}
