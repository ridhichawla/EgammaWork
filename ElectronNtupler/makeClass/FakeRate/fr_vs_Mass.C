#define fr_vs_Mass_cxx
#include "fr_vs_Mass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void fr_vs_Mass::Loop()
{

  TFile *file = new TFile("SPhoton_forFakeRate.root", "recreate");
  TTree *tree = new TTree("tree"," after preselections tree");

  // Branch variable declaration
  double Ele1PT;
  double Ele2PT;
  double Ele1Eta;
  double Ele2Eta;
  double diEle_Mass;

  // Branch declaration
  tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D"); 
  tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");         
  tree->Branch("Ele1Eta", &Ele1Eta, "Ele1Eta/D");         
  tree->Branch("Ele2Eta", &Ele2Eta, "Ele2Eta/D");
  tree->Branch("diEle_Mass", &diEle_Mass, "diEle_Mass/D");

  int count, good_elec;
  bool passID;
  bool passKin1;
  bool passKin2;
  TLorentzVector ele1,ele2,diMass;

  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi;

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 100;
  cout<<"entries: "<<nentries<<endl;

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

    good_elec = 0;
    count = 0;
    passID = false;
    passKin1 = false;
    passKin2 = false;

    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear();

    if(metPt->at(0) > 10.) continue;

    if(singlePhoton){
      for(int j=0;j<nEle;j++){

	passID = passMediumId->at(index[j]);
	passKin1 = pt->at(index[j]) > 20.;
	passKin2 = fabs(eta->at(index[j])) < 2.5;

	if(passID && passKin1 && passKin2) count++;

      } // nEle

      if(count <= 1){
	for(int k=0;k<nEle;k++){

	  if(expectedMissingInnerHits->at(index[k]) == 0){

	    good_elec = good_elec + 1;

	    newelePt.push_back(pt->at(index[k]));
	    neweleEta.push_back(eta->at(index[k]));
	    newelePhi.push_back(phi->at(index[k]));
	    neweleEnr.push_back(energy->at(index[k]));

	  } // missing hits
	} // nEle
      } // count

      for(unsigned int l=0;l<newelePt.size();l++){

	if(good_elec==2){
	  ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	  diMass=ele1+ele2;
	  diEle_Mass = diMass.M();

	  Ele1PT = newelePt.at(0);
	  Ele1Eta = neweleEta.at(0);
	  Ele2PT = newelePt.at(1);
	  Ele2Eta = neweleEta.at(1);
	}

	tree->Fill();

      } // newelePt
    } // trigger
  } // event

  file->Write();
  file->Close();
}
