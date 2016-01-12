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

  int count;
  bool passID, passECAL;

  int isNum_BRL, isDen_BRL;
  int isNum_ECAP1, isDen_ECAP1;
  int isNum_ECAP2, isDen_ECAP2;

  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> newelePassMedium;
  vector <double> newelePassECAL;

  TH1F *numPt_BRL = new TH1F("numPt_BRL", "numPt_BRL", 100, 0, 500);
  TH1F *numPt_ECAP1 = new TH1F("numPt_ECAP1", "numPt_ECAP1", 100, 0, 500);
  TH1F *numPt_ECAP2 = new TH1F("numPt_ECAP2 ", "numPt_ECAP2 ", 100, 0, 500);

  TH1F *denPt_BRL = new TH1F("denPt_BRL", "denPt_BRL", 100, 0, 500);
  TH1F *denPt_ECAP1 = new TH1F("denPt_ECAP1", "denPt_ECAP1", 100, 0, 500);
  TH1F *denPt_ECAP2  = new TH1F("denPt_ECAP2 ", "denPt_ECAP2 ", 100, 0, 500);

  numPt_BRL->Sumw2(); denPt_BRL->Sumw2();
  numPt_ECAP1->Sumw2(); denPt_ECAP1->Sumw2();
  numPt_ECAP2->Sumw2(); denPt_ECAP2->Sumw2();

  isNum_BRL = 0; isDen_BRL = 0;
  isNum_ECAP1 = 0; isDen_ECAP1 = 0;
  isNum_ECAP2 = 0; isDen_ECAP2 = 0;

  if (fChain == 0) return;

  //Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries = 5000;
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

    count = 0;
    newelePt.clear(); neweleEta.clear(); newelePassMedium.clear(); newelePassECAL.clear();

    //if(!singlePhoton_175) continue;
    if(metPt->at(0) > 10.) continue;

    for(int j=0;j<nEle;j++){

      passID = passMediumId->at(index[j]);
      passECAL = eleEcalDrivenSeed->at(index[j]);
      if(passID && passECAL) count++;

    } //nEle

    if(count <= 1){
      for(int k=0;k<nEle;k++){

	if(expectedMissingInnerHits->at(index[k]) == 0){
	  if(fabs(eta->at(index[k])) < 2.5){

	    newelePt.push_back(pt->at(index[k]));
	    neweleEta.push_back(eta->at(index[k]));
	    newelePassMedium.push_back(passMediumId->at(index[k]));
	    newelePassECAL.push_back(eleEcalDrivenSeed->at(index[k]));

	  } // eta
	} // missing hits
      } // nEle
    } // count

    // // event

    for(unsigned int l=0;l<newelePt.size();l++){
      if(fabs(neweleEta.at(l)) < 1.4442){
	isDen_BRL++;
	denPt_BRL->Fill(newelePt.at(l));
      }

      if(fabs(neweleEta.at(l)) > 1.566 && fabs(neweleEta.at(l)) < 2.0){
	isDen_ECAP1++;
	denPt_ECAP1->Fill(newelePt.at(l));
      }

      if(fabs(neweleEta.at(l)) > 2.0 && fabs(neweleEta.at(l)) < 2.5){
	isDen_ECAP2++;
	denPt_ECAP2->Fill(newelePt.at(l));
      }

      if(newelePassMedium.at(l) == 1 && newelePassECAL.at(l) == 1){
	if(newelePt.at(l) > 20.){

	  if(fabs(neweleEta.at(l)) < 1.4442){
	    isNum_BRL++;
	    numPt_BRL->Fill(newelePt.at(l));
	  }

	  if(fabs(neweleEta.at(l)) > 1.566 && fabs(neweleEta.at(l)) < 2.0){
	    isNum_ECAP1++;
	    numPt_ECAP1->Fill(newelePt.at(l));
	  }

	  if(fabs(neweleEta.at(l)) > 2.0 && fabs(neweleEta.at(l)) < 2.5){
	    isNum_ECAP2++;
	    numPt_ECAP2->Fill(newelePt.at(l));
	  }

	} // pt
      } // pass Medium & ECAL driven
    } // l

  } // event

  cout<<"Numerator Barrel: "<<isNum_BRL<<"   "<<"Denominator Barrel: "<<isDen_BRL<<endl;
  cout<<"Numerator Endcap 1: "<<isNum_ECAP1<<"   "<<"Denominator Endcap 1: "<<isDen_ECAP1<<endl;
  cout<<"Numerator Endcap 2: "<<isNum_ECAP2<<"   "<<"Denominator Endcap 2: "<<isDen_ECAP2<<endl;

  file->Write();
  file->Close();
}
