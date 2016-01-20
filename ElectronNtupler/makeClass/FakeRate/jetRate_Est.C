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
  bool passID;
  //bool passECAL;
  //bool passKin1;
  bool passKin2;

  int isNum, isDen;
  int isNum_BRL, isDen_BRL;
  int isNum_ECAP1, isDen_ECAP1;
  int isNum_ECAP2, isDen_ECAP2;

  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> newelePassMedium;
  vector <double> newelePassECAL;

  Double_t xbins[16] = {0,5,10,15,20,25,30,35,40,50,60,70,90,150,200,300};

  /*TH1F *numPt_BRL = new TH1F("numPt_BRL", "numPt_BRL", 100, 0, 500);
    TH1F *numPt_ECAP1 = new TH1F("numPt_ECAP1", "numPt_ECAP1", 100, 0, 500);
    TH1F *numPt_ECAP2 = new TH1F("numPt_ECAP2 ", "numPt_ECAP2 ", 100, 0, 500);

    TH1F *denPt_BRL = new TH1F("denPt_BRL", "denPt_BRL", 100, 0, 500);
    TH1F *denPt_ECAP1 = new TH1F("denPt_ECAP1", "denPt_ECAP1", 100, 0, 500);
    TH1F *denPt_ECAP2  = new TH1F("denPt_ECAP2 ", "denPt_ECAP2 ", 100, 0, 500);*/

  TH1F *etPhoton = new TH1F("etPhoton", "etPhoton", 100, 0, 700);
  TH1F *etPhoton_preScale = new TH1F("etPhoton_preScale", "etPhoton_preScale", 100, 0, 700);

  TH1F *numPt     = new TH1F("numPt", "numPt", 15, xbins);
  TH1F *numPt_BRL = new TH1F("numPt_BRL", "numPt_BRL", 15, xbins);
  TH1F *numPt_ECAP1 = new TH1F("numPt_ECAP1", "numPt_ECAP1", 15, xbins);
  TH1F *numPt_ECAP2 = new TH1F("numPt_ECAP2 ", "numPt_ECAP2 ", 15, xbins);

  TH1F *denPt = new TH1F("denPt", "denPt", 15, xbins);
  TH1F *denPt_BRL = new TH1F("denPt_BRL", "denPt_BRL", 15, xbins);
  TH1F *denPt_ECAP1 = new TH1F("denPt_ECAP1", "denPt_ECAP1", 15, xbins);
  TH1F *denPt_ECAP2  = new TH1F("denPt_ECAP2 ", "denPt_ECAP2 ", 15, xbins);

  etPhoton->Sumw2(); etPhoton_preScale->Sumw2();
  numPt->Sumw2(); denPt->Sumw2();
  numPt_BRL->Sumw2(); denPt_BRL->Sumw2();
  numPt_ECAP1->Sumw2(); denPt_ECAP1->Sumw2();
  numPt_ECAP2->Sumw2(); denPt_ECAP2->Sumw2();

  isNum = 0; isDen = 0;
  isNum_BRL = 0; isDen_BRL = 0;
  isNum_ECAP1 = 0; isDen_ECAP1 = 0;
  isNum_ECAP2 = 0; isDen_ECAP2 = 0;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 50000;
  cout<<"entries: "<<nentries<<endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%500000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    int index[pt->size()];
    float ptnew[pt->size()];

    for(unsigned int el=0; el<pt->size(); el++) {
      ptnew[el]=pt->at(el); }

    int size1 = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(size1,ptnew,index,true);

    count = 0;
    passID = false;
    //passECAL = false;
    //passKin1 = false;
    passKin2 = false;
    newelePt.clear(); neweleEta.clear(); newelePassMedium.clear(); newelePassECAL.clear();

    if(metPt->at(0) > 10.) continue;

    //if(singlePhoton_30 || singlePhoton_36 || singlePhoton_50 || singlePhoton_75 || singlePhoton_90 || singlePhoton_120 || singlePhoton_175)
    if(singlePhoton){

      etPhoton->Fill(et_Photon->at(0));
      etPhoton_preScale->Fill(et_Photon->at(0),prescalePhoton);

      for(int j=0;j<nEle;j++){

	passID = passMediumId->at(index[j]);
	//passECAL = eleEcalDrivenSeed->at(index[j]);
	//passKin1 = pt->at(index[j]) > 20.;
	passKin2 = fabs(eta->at(index[j])) < 2.5;

	if(passID && passKin2) count++;

      } //nEle

      if(count <= 1){
	for(int k=0;k<nEle;k++){

	  if(expectedMissingInnerHits->at(index[k]) == 0){
	    //if(pt->at(index[k]) > 20. && fabs(eta->at(index[k])) < 2.5){

	    newelePt.push_back(pt->at(index[k]));
	    neweleEta.push_back(eta->at(index[k]));
	    newelePassMedium.push_back(passMediumId->at(index[k]));
	    newelePassECAL.push_back(eleEcalDrivenSeed->at(index[k]));

	    //} // eta
	  } // missing hits
	} // nEle
      } // count

      // // event

      for(unsigned int l=0;l<newelePt.size();l++){
	
	isDen++;
	denPt->Fill(newelePt.at(l),prescalePhoton);

	if(fabs(neweleEta.at(l)) < 1.4442){
	  isDen_BRL++;
	  denPt_BRL->Fill(newelePt.at(l),prescalePhoton);
	}

	if(fabs(neweleEta.at(l)) > 1.566 && fabs(neweleEta.at(l)) < 2.0){
	  isDen_ECAP1++;
	  denPt_ECAP1->Fill(newelePt.at(l),prescalePhoton);
	}

	if(fabs(neweleEta.at(l)) > 2.0 && fabs(neweleEta.at(l)) < 2.5){
	  isDen_ECAP2++;
	  denPt_ECAP2->Fill(newelePt.at(l),prescalePhoton);
	}

	if(newelePassMedium.at(l) == 1){
	  //if(newelePt.at(l) > 20.){

	  isNum++;
	  numPt->Fill(newelePt.at(l),prescalePhoton);
	  
	  if(fabs(neweleEta.at(l)) < 1.4442){
	    isNum_BRL++;
	    numPt_BRL->Fill(newelePt.at(l),prescalePhoton);
	  }

	  if(fabs(neweleEta.at(l)) > 1.566 && fabs(neweleEta.at(l)) < 2.0){
	    isNum_ECAP1++;
	    numPt_ECAP1->Fill(newelePt.at(l),prescalePhoton);
	  }

	  if(fabs(neweleEta.at(l)) > 2.0 && fabs(neweleEta.at(l)) < 2.5){
	    isNum_ECAP2++;
	    numPt_ECAP2->Fill(newelePt.at(l),prescalePhoton);
	  }

	  //} // pt
	} // pass Medium
      } // l
    } // trigger
  } // event

  cout<<"Numerator: "<<isNum<<"   "<<"Denominator: "<<isDen<<endl;
  cout<<"Numerator Barrel: "<<isNum_BRL<<"   "<<"Denominator Barrel: "<<isDen_BRL<<endl;
  cout<<"Numerator Endcap 1: "<<isNum_ECAP1<<"   "<<"Denominator Endcap 1: "<<isDen_ECAP1<<endl;
  cout<<"Numerator Endcap 2: "<<isNum_ECAP2<<"   "<<"Denominator Endcap 2: "<<isDen_ECAP2<<endl;

  file->Write();
  file->Close();
}
