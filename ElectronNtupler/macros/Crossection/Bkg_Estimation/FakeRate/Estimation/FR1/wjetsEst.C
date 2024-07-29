#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void wjetsEst_fromData(){

  TFile f1("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Data/SE_2015.root");
  TTree *T1 = (TTree*)f1.Get("ntupler/ElectronTree");

  Bool_t          Ele23_WPLoose;
  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
  vector<float>   *etaSC;
  vector<int>     *passMediumId;

  ptElec = 0;
  etaElec = 0;
  phiElec = 0;
  energyElec = 0;
  etaSC = 0;
  passMediumId = 0;

  T1->SetBranchStatus("*",0);
  T1->SetBranchStatus("Ele23_WPLoose", 1);
  T1->SetBranchStatus("ptElec", 1);
  T1->SetBranchStatus("etaElec", 1);
  T1->SetBranchStatus("phiElec", 1);
  T1->SetBranchStatus("energyElec", 1);
  T1->SetBranchStatus("etaSC", 1);
  T1->SetBranchStatus("passMediumId", 1);

  T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
  T1->SetBranchAddress("ptElec", &ptElec);
  T1->SetBranchAddress("etaElec", &etaElec);
  T1->SetBranchAddress("phiElec", &phiElec);
  T1->SetBranchAddress("energyElec", &energyElec);
  T1->SetBranchAddress("etaSC", &etaSC);
  T1->SetBranchAddress("passMediumId", &passMediumId);

  TFile *file = new TFile("WJets_fromData_ControlRegion.root", "recreate");
  TTree *tree = new TTree("tree"," after preselections tree");

  int count, count1;
  int mediumId;
  bool passKin;
  bool ptcut;
  TLorentzVector ele1,ele2,diReco;

  // Branch variable declaration
  double ElePT_Fail, ElePT_Pass;
  double EleEta_Fail, EleEta_Pass;
  double wjetEstMass;
  bool BB, BE, EE;

  // Branch declaration
  tree->Branch("ElePT_Fail", &ElePT_Fail, "ElePT_Fail/D");
  tree->Branch("EleEta_Fail", &EleEta_Fail, "EleEta_Fail/D");
  tree->Branch("ElePT_Pass", &ElePT_Pass, "ElePT_Pass/D");
  tree->Branch("EleEta_Pass", &EleEta_Pass, "EleEta_Pass/D");
  tree->Branch("wjetEstMass", &wjetEstMass, "wjetEstMass/D");
  tree->Branch("BB", &BB, "BB/B");
  tree->Branch("BE", &BE, "BE/B");
  tree->Branch("EE", &EE, "EE/B");

  vector<int> idx; vector<int> idx1;

  int nentries = T1->GetEntries();
  //int nentries = 1000;
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

    count = 0; count1 = 0;
    mediumId = 0;
    passKin = false;
    ptcut = false;
    bool etacut = false;
    idx.clear();
    idx1.clear();
    BB=false; BE=false; EE=false;

    if(!Ele23_WPLoose) continue;        // trigger not satisfied

    for(int i=0;i<ptElec->size();i++){

      mediumId = passMediumId->at(index[i]);
      //passKin = (fabs(etaSC->at(index[i])) < 2.5 && !(fabs(etaSC->at(index[i])) > 1.4442 && fabs(etaSC->at(index[i])) < 1.566));

      if(mediumId){
	count++;
	idx.push_back(index[i]);
      }

      if(!mediumId){
	count1++;
	idx1.push_back(index[i]);
      }

    } // pt size

    if(count == 1 && count1 == 1){

      //cout<<"idx size: "<<idx.size()<<"   "<<"idx1 size: "<<idx1.size()<<endl;

      if((idx[0] < idx1[0]) && ptElec->at(idx[0]) > 30 && ptElec->at(idx1[0]) > 10.) ptcut = true;
      if((idx[0] > idx1[0]) && ptElec->at(idx1[0]) > 30 && ptElec->at(idx[0]) > 10.) ptcut = true; 

      if((fabs(etaSC->at(idx[0])) < 2.5 && !(fabs(etaSC->at(idx[0])) > 1.4442 && fabs(etaSC->at(idx[0])) < 1.566)) && (fabs(etaSC->at(idx1[0])) < 2.5 && !(fabs(etaSC->at(idx1[0])) > 1.4442 && fabs(etaSC->at(idx1[0])) < 1.566))) etacut = true; 

      if(ptcut && etacut) {

	ElePT_Fail  = ptElec->at(idx1[0]);
	EleEta_Fail = etaElec->at(idx1[0]);
	ElePT_Pass  = ptElec->at(idx[0]);
	EleEta_Pass = etaElec->at(idx[0]);

	ele1.SetPtEtaPhiE(ptElec->at(idx[0]),etaElec->at(idx[0]),phiElec->at(idx[0]),energyElec->at(idx[0]));
	ele2.SetPtEtaPhiE(ptElec->at(idx1[0]),etaElec->at(idx1[0]),phiElec->at(idx1[0]),energyElec->at(idx1[0]));

	diReco=ele1+ele2;
	wjetEstMass = diReco.M();

	if(fabs(etaSC->at(idx[0])) < 1.4442 && fabs(etaSC->at(idx1[0])) < 1.4442) BB = true;
	if((fabs(etaSC->at(idx[0])) < 1.4442 && fabs(etaSC->at(idx1[0])) > 1.566) || (fabs(etaSC->at(idx[0])) > 1.566 && fabs(etaSC->at(idx1[0])) < 1.4442)) BE =true;
	if(fabs(etaSC->at(idx[0])) > 1.566 && fabs(etaSC->at(idx1[0])) > 1.566) EE =true;

	tree->Fill();
      }
    }// count==1

  } // event

  file->Write();
  file->Close();
}
