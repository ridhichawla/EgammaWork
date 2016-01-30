#define jet_FakeRate_cxx
#include "jet_FakeRate.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void jet_FakeRate::Loop()
{
  TFile *file = new TFile("single_photon.root", "recreate");

  int count, good_elec;
  bool passID;
  bool passKin1;
  bool passKin2;

  int isNum, isDen;
  int isNum_BRL, isDen_BRL;
  int isNum_ECAP1, isDen_ECAP1;
  int isNum_ECAP2, isDen_ECAP2;

  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> newelePhi;
  vector <double> neweleEnr;
  vector <double> newelePassMedium;

  Double_t x1bin[18] = {0,5,10,15,20,25,30,35,40,50,60,70,90,150,200,300,400,500};
  Double_t x2bin[6] = {0.,0.5,1.0,1.5,2.0,2.5};

  TH1F *et_HLT = new TH1F("et_HLT", "et_HLT", 100, 0, 700);
  TH1F *et_HLT_preScale = new TH1F("et_HLT_preScale", "et_HLT_preScale", 100, 0, 700);

  TH1F *numPt     = new TH1F("numPt", "numPt", 17, x1bin);
  TH2F *numPt_Eta = new TH2F("numPt_Eta", "numPt_Eta", 17, x1bin, 5, x2bin);
  TH1F *numPt_BRL = new TH1F("numPt_BRL", "numPt_BRL", 17, x1bin);
  TH1F *numPt_ECAP1 = new TH1F("numPt_ECAP1", "numPt_ECAP1", 17, x1bin);
  TH1F *numPt_ECAP2 = new TH1F("numPt_ECAP2 ", "numPt_ECAP2 ", 17, x1bin);

  TH1F *denPt = new TH1F("denPt", "denPt", 17, x1bin);
  TH2F *denPt_Eta = new TH2F("denPt_Eta", "denPt_Eta", 17, x1bin, 5, x2bin);
  TH1F *denPt_BRL = new TH1F("denPt_BRL", "denPt_BRL", 17, x1bin);
  TH1F *denPt_ECAP1 = new TH1F("denPt_ECAP1", "denPt_ECAP1", 17, x1bin);
  TH1F *denPt_ECAP2  = new TH1F("denPt_ECAP2 ", "denPt_ECAP2 ", 17, x1bin);

  TH1F *numEta     = new TH1F("numEta", "numEta", 5, x2bin);
  TH1F *denEta     = new TH1F("denEta", "denEta", 5, x2bin);

  //TH1F *metWJets   = new TH1F("metWJets", "metWJets", 20, 0, 200);

  et_HLT->Sumw2(); et_HLT_preScale->Sumw2();

  numPt->Sumw2(); denPt->Sumw2();
  numPt_Eta->Sumw2(); denPt_Eta->Sumw2();
  numPt_BRL->Sumw2(); denPt_BRL->Sumw2();
  numPt_ECAP1->Sumw2(); denPt_ECAP1->Sumw2();
  numPt_ECAP2->Sumw2(); denPt_ECAP2->Sumw2();
  numEta->Sumw2(); denEta->Sumw2();

  //metWJets->Sumw2();

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

    good_elec = 0;
    count = 0;
    passID = false;
    passKin1 = false;
    passKin2 = false;
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); newelePassMedium.clear();

    //metWJets->Fill(metPt->at(0),theWeight);

    if(metPt->at(0) > 10.) continue;

    if(singlePhoton){

      et_HLT->Fill(et_Photon->at(0));
      et_HLT_preScale->Fill(et_Photon->at(0),prescalePhoton);

      for(int j=0;j<nEle;j++){

	passID = passMediumId->at(index[j]);
	passKin1 = pt->at(index[j]) > 20.;
	passKin2 = fabs(eta->at(index[j])) < 2.5;

	if(passID && passKin1 && passKin2) count++;

      } //nEle

      if(count <= 1){
	for(int k=0;k<nEle;k++){

	  if(expectedMissingInnerHits->at(index[k]) == 0){
	    good_elec = good_elec + 1;

	    newelePt.push_back(pt->at(index[k]));
	    neweleEta.push_back(eta->at(index[k]));
	    newelePhi.push_back(phi->at(index[k]));
	    neweleEnr.push_back(energy->at(index[k]));
	    newelePassMedium.push_back(passMediumId->at(index[k]));

	  } // missing hits
	} // nEle
      } // count

      for(unsigned int l=0;l<newelePt.size();l++){

	isDen++;
	denPt->Fill(newelePt.at(l),prescalePhoton);
	denEta->Fill(neweleEta.at(l),prescalePhoton);
	denPt_Eta->Fill(newelePt.at(l),neweleEta.at(l),prescalePhoton);

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

	  isNum++;
	  numPt->Fill(newelePt.at(l),prescalePhoton);
	  numEta->Fill(neweleEta.at(l),prescalePhoton);
	  numPt_Eta->Fill(newelePt.at(l),neweleEta.at(l),prescalePhoton);

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

	} // pass Medium
      } // newelePt size
    } // trigger
  } // event

  cout<<"Numerator: "<<isNum<<"   "<<"Denominator: "<<isDen<<endl;
  cout<<"Numerator Barrel: "<<isNum_BRL<<"   "<<"Denominator Barrel: "<<isDen_BRL<<endl;
  cout<<"Numerator Endcap 1: "<<isNum_ECAP1<<"   "<<"Denominator Endcap 1: "<<isDen_ECAP1<<endl;
  cout<<"Numerator Endcap 2: "<<isNum_ECAP2<<"   "<<"Denominator Endcap 2: "<<isDen_ECAP2<<endl;

  file->Write();
  file->Close();
}
