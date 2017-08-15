#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void jet_FakeRate() {

  TFile f1("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_07082017/Data/SE_2015.root");
  TTree *T1 = (TTree*)f1.Get("ntupler/ElectronTree");

  Bool_t          Ele8_PFJet;
  Bool_t          Ele12_PFJet;
  Bool_t          Ele23_PFJet;
  Bool_t          Ele33_PFJet;
  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
  vector<float>   *etaSC;
  vector<float>   *full5x5_sigmaIetaIeta;
  vector<float>   *dEtaIn;
  vector<float>   *dPhiIn;
  vector<float>   *hOverE;
  vector<int>     *passMediumId;

  ptElec = 0;
  etaElec = 0;
  phiElec = 0;
  energyElec = 0;
  etaSC = 0;
  full5x5_sigmaIetaIeta = 0;
  dEtaIn = 0;
  dPhiIn = 0;
  hOverE = 0;
  passMediumId = 0; 

  T1->SetBranchStatus("*", 0);
  T1->SetBranchStatus("Ele8_PFJet", 1);
  T1->SetBranchStatus("Ele12_PFJet", 1);
  T1->SetBranchStatus("Ele23_PFJet", 1);
  T1->SetBranchStatus("Ele33_PFJet", 1);
  T1->SetBranchStatus("ptElec", 1);
  T1->SetBranchStatus("etaElec", 1);
  T1->SetBranchStatus("phiElec", 1);
  T1->SetBranchStatus("energyElec", 1);
  T1->SetBranchStatus("etaSC", 1);
  T1->SetBranchStatus("full5x5_sigmaIetaIeta", 1);
  T1->SetBranchStatus("dEtaIn", 1);
  T1->SetBranchStatus("dPhiIn", 1);
  T1->SetBranchStatus("hOverE", 1);
  T1->SetBranchStatus("passMediumId", 1);

  T1->SetBranchAddress("Ele8_PFJet", &Ele8_PFJet);
  T1->SetBranchAddress("Ele12_PFJet", &Ele12_PFJet);
  T1->SetBranchAddress("Ele23_PFJet", &Ele23_PFJet);
  T1->SetBranchAddress("Ele33_PFJet", &Ele33_PFJet);
  T1->SetBranchAddress("ptElec", &ptElec);
  T1->SetBranchAddress("etaElec", &etaElec);
  T1->SetBranchAddress("phiElec", &phiElec);
  T1->SetBranchAddress("energyElec", &energyElec);
  T1->SetBranchAddress("etaSC", &etaSC);
  T1->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta);
  T1->SetBranchAddress("dEtaIn", &dEtaIn);
  T1->SetBranchAddress("dPhiIn", &dPhiIn);
  T1->SetBranchAddress("hOverE", &hOverE);
  T1->SetBranchAddress("passMediumId", &passMediumId);

  TFile *file = new TFile("FakeRate_SE15.root", "recreate");

  int mediumId_barrel, mediumId_endcap;
  double sigma_ieta, dEta, dPhi, hovere;
  vector<double> isden_barrel; vector<double> isden_endcap;

  Double_t x1bin[6] = {10,20,30,40,50,10000};
  //Double_t x1bin[10] = {10,15,20,25,30,40,50,70,100,10000};
  int nbins = 5;

  TH1D *numPt_barrel  = new TH1D("numPt_barrel", "numPt_barrel", nbins, x1bin);
  TH1D *numPt_endcap = new TH1D("numPt_endcap", "numPt_endcap", nbins, x1bin);

  TH1D *denPt_barrel  = new TH1D("denPt_barrel", "denPt_barrel", nbins, x1bin);
  TH1D *denPt_endcap = new TH1D("denPt_endcap", "denPt_endcap", nbins, x1bin);

  TH1D *isden_Barrel = new TH1D("isden_Barrel", "isden_Barrel", 5, 0, 5);
  TH1D *isden_Endcap = new TH1D("isden_Endcap", "isden_Endcap", 5, 0, 5);

  numPt_barrel->Sumw2(); denPt_barrel->Sumw2();
  numPt_endcap->Sumw2(); denPt_endcap->Sumw2();
  isden_Barrel->Sumw2(); isden_Endcap->Sumw2();

  //int nentries = T1->GetEntries();
  int nentries = 100000;
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int jentry=0; jentry < nentries; jentry++) {
    T1->GetEntry(jentry);

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    // Sorting
    int index[ptElec->size()];
    float pt[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++) {
      pt[el]=ptElec->at(el); }

    int size = sizeof(pt)/sizeof(pt[0]);
    TMath::Sort(size,pt,index,true);

    mediumId_barrel = 0; mediumId_endcap = 0;
    sigma_ieta = 0.; dEta = 0.; dPhi = 0.; hovere = 0.;
    isden_barrel.clear(); isden_endcap.clear();

    //cout<<"Entry = "<<jentry<<endl;
    //cout<<"Nele = "<<ptElec->size()<<endl;
    //cout<<"Ele8_PFJet = "<<Ele8_PFJet<<"   Ele12_PFJet = "<<Ele12_PFJet<<"   Ele23_PFJet = "<<Ele23_PFJet<<"   Ele33_PFJet = "<<Ele33_PFJet<<endl;

    if(Ele8_PFJet || Ele12_PFJet || Ele23_PFJet || Ele33_PFJet){

      cout<<"Entry = "<<jentry<<endl;
      cout<<"Nele = "<<ptElec->size()<<endl;
      
      for(int i=0;i<ptElec->size();i++){

	cout<<"Pt = "<<ptElec->at(index[i])<<"   SC Eta = "<<etaSC->at(index[i])<<endl;

	if(ptElec->at(index[i]) > 10. && fabs(etaSC->at(index[i])) < 2.5 && !(fabs(etaSC->at(index[i])) > 1.4442 && fabs(etaSC->at(index[i])) < 1.566)){

	  double sigma_ieta = full5x5_sigmaIetaIeta->at(index[i]);
	  double dEta = dEtaIn->at(index[i]);
	  double dPhi = dPhiIn->at(index[i]);
	  double hovere = hOverE->at(index[i]);

	  cout<<"sigma_ieta = "<<sigma_ieta<<"   dEta = "<<dEta<<"   dPhi = "<<dPhi<<"   hovere = "<<hovere<<endl;

	  if(fabs(etaSC->at(index[i])) < 1.4442 && sigma_ieta < 0.011 && fabs(dEta) < 0.01 && fabs(dPhi) < 0.15 && hovere < 0.1) isden_barrel.push_back(index[i]);
	  if(fabs(etaSC->at(index[i])) > 1.566 && sigma_ieta < 0.031 && hovere < 0.12) isden_endcap.push_back(index[i]);
	}
      }

      cout<<"isden_barrel size = "<<isden_barrel.size()<<"   isden_endcap size = "<<isden_endcap.size()<<endl;

      isden_Barrel->Fill(isden_barrel.size());
      isden_Endcap->Fill(isden_endcap.size());
      
      for(int j=0; j<isden_barrel.size(); j++){

	cout<<"isden_barrel["<<j<<"] = "<<j<<endl;

	cout<<"****************Denominator Barrel****************"<<endl;
	cout<<"Pt["<<j<<"] = "<<ptElec->at(j)<<"   Eta["<<j<<"] = "<<etaSC->at(j)<<endl;
	cout<<"Medium Id["<<j<<"] = "<<passMediumId->at(j)<<endl;

	denPt_barrel->Fill(ptElec->at(j));

	int mediumId_barrel = passMediumId->at(j);

	if(mediumId_barrel){

	  cout<<"*****************Numerator Barrel*****************"<<endl;
	  cout<<"Pt["<<j<<"] = "<<ptElec->at(j)<<"   Eta["<<j<<"] = "<<etaSC->at(j)<<endl;
	  numPt_barrel->Fill(ptElec->at(j));
	}
      }

      for(int k=0; k<isden_endcap.size(); k++){

	cout<<"isden_endcap["<<k<<"] = "<<k<<endl;

	cout<<"****************Denominator Endcap****************"<<endl;
	cout<<"Pt["<<k<<"] = "<<ptElec->at(k)<<"   Eta["<<k<<"] = "<<etaSC->at(k)<<endl;
	cout<<"Medium Id["<<k<<"] = "<<passMediumId->at(k)<<endl;

	denPt_endcap->Fill(ptElec->at(k));

	int mediumId_endcap = passMediumId->at(k);

	if(mediumId_endcap){

	  cout<<"*****************Numerator Endcap*****************"<<endl;
	  cout<<"Pt["<<k<<"] = "<<ptElec->at(k)<<"   Eta["<<k<<"] = "<<etaSC->at(k)<<endl;
	  numPt_endcap->Fill(ptElec->at(k));
	}
      }

      cout<<""<<endl;

    } // trigger

    //cout<<""<<endl;
  } // event

  file->Write();
  file->Close();
}
