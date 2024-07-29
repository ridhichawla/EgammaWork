#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void dijetEst_fromData(){

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

  TFile *file = new TFile("DiJet_fromData_ControlRegion.root", "recreate");
  TTree *tree = new TTree("tree"," after preselections tree");

  int count, cnt1, cnt2, cnt3;
  int mediumId;
  bool passKin;
  TLorentzVector ele1,ele2,diReco;

  // Branch variable declaration
  double Ele1PT; double Ele2PT;
  double Ele1Eta; double Ele2Eta;
  double qcdEstMass;

  // Branch declaration
  tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D");
  tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");
  tree->Branch("Ele1Eta", &Ele1Eta, "Ele1Eta/D");
  tree->Branch("Ele2Eta", &Ele2Eta, "Ele2Eta/D");
  tree->Branch("qcdEstMass", &qcdEstMass, "qcdEstMass/D");

  vector<int> idx;
  cnt1 = 0.;
  cnt2 = 0.;
  cnt3 = 0.;

  int nentries = T1->GetEntries();
  //int nentries = 20;
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

    count = 0;
    mediumId = 0;
    passKin = false;
    bool etacut = false;
    idx.clear();

    //cout<<"Trigger = "<<Ele23_WPLoose<<endl;
    //cout<<"No. of electrons = "<<ptElec->size()<<endl;

    if(!Ele23_WPLoose) continue;        // trigger not satisfied
    //if(Ele23_WPLoose) continue;        // trigger satisfied

    cnt1++;
    
    for(int i=0;i<ptElec->size();i++){

      //if(ptElec->at(index[i]) > 1000.) cout<<"Ele PT: "<<ptElec->at(index[i])<<endl;
      mediumId = passMediumId->at(index[i]);
      //passKin = (fabs(etaSC->at(index[i])) < 2.5 && !(fabs(etaSC->at(index[i])) > 1.4442 && fabs(etaSC->at(index[i])) < 1.566));

      if(!mediumId){
        //cout<<"Event = "<<jentry<<"   No. of electrons = "<<ptElec->size()<<endl;
	count++;
	idx.push_back(index[i]);
      }
    }

    //cout<<""<<endl;

    //cout<<"No. of electrons = "<<ptElec->size()<<endl;
    if(count == 2) cnt2++;
    if(count > 2) cnt3++;
    
    if(count == 2){

      //if(idx.size() != 2) cout<<"idx size: "<<idx.size()<<endl;
      //cout<<"idx[0]: "<<idx[0]<<"   "<<"idx[1]: "<<idx[1]<<endl;
      //cout<<""<<endl;
      
      if((fabs(etaSC->at(idx[0])) < 2.5 && !(fabs(etaSC->at(idx[0])) > 1.4442 && fabs(etaSC->at(idx[0])) < 1.566)) && (fabs(etaSC->at(idx[1])) < 2.5 && !(fabs(etaSC->at(idx[1])) > 1.4442 && fabs(etaSC->at(idx[1])) < 1.566))) etacut = true;

      if(ptElec->at(idx[0]) > 30 && ptElec->at(idx[1]) > 10 && etacut){
	Ele1PT  = ptElec->at(idx[0]);
	Ele1Eta = etaElec->at(idx[0]);

	//if(ptElec->at(idx[0]) > 1000. || ptElec->at(idx[1]) > 1000.) cout<<"Ele1 PT: "<<ptElec->at(idx[0])<<"   "<<"Ele2 PT: "<<ptElec->at(idx[1])<<endl;

	Ele2PT  = ptElec->at(idx[1]);
	Ele2Eta = etaElec->at(idx[1]);

	ele1.SetPtEtaPhiE(ptElec->at(idx[0]),etaElec->at(idx[0]),phiElec->at(idx[0]),energyElec->at(idx[0]));
	ele2.SetPtEtaPhiE(ptElec->at(idx[1]),etaElec->at(idx[1]),phiElec->at(idx[1]),energyElec->at(idx[1]));

	diReco=ele1+ele2;
	qcdEstMass = diReco.M();

	tree->Fill();

      } // pt
    } // count
  } // event

  cout<<"Total Events = "<<cnt1<<"   Events with ele == 2: "<<cnt2<<"   Events with ele > 2: "<<cnt3<<endl;

  file->Write();
  file->Close();
}
