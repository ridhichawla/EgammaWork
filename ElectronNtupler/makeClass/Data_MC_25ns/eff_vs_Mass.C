#define eff_vs_Mass_cxx
#include "eff_vs_Mass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void eff_vs_Mass::Loop()
{

  TFile *file = new TFile("DY_forEff.root", "recreate");
  TTree *tree = new TTree("tree"," after preselections tree");

  // Branch variable declaration
  double Ele1PT;
  double Ele2PT;
  double Ele1EtaSC;
  double Ele2EtaSC;
  double ZMass;
  double genWeights;

  // Branch declaration
  tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D"); 
  tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");         
  tree->Branch("Ele1EtaSC", &Ele1EtaSC, "Ele1EtaSC/D");         
  tree->Branch("Ele2EtaSC", &Ele2EtaSC, "Ele2EtaSC/D");
  tree->Branch("ZMass", &ZMass, "ZMass/D");
  tree->Branch("genWeights", &genWeights, "genWeights/D");

  int good_elec;
  TLorentzVector ele1,ele2,dielectron;

  vector <double> newscEt; vector <double> newscEta; vector <double> newscPhi; vector <double> newscEnr;
  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

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

    //Z_Mass=0.;
    good_elec = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    //cout<<"1"<<endl;

    genWeights = theWeight;

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

	  Ele1PT = newelePt.at(0);
	  Ele1EtaSC = newscEta.at(0);
	  Ele2PT = newelePt.at(1);
	  Ele2EtaSC = newscEta.at(1);

	  ele1.SetPtEtaPhiE(newscEt.at(0),newscEta.at(0),newscPhi.at(0),newscEnr.at(0));
	  ele2.SetPtEtaPhiE(newscEt.at(1),newscEta.at(1),newscPhi.at(1),newscEnr.at(1));

	  dielectron=ele1+ele2;
	  ZMass = dielectron.M();

	  tree->Fill();

	} // pt
      } // good electrons

    } // tauFlag
  } // event

  file->Write();
  file->Close();
}
