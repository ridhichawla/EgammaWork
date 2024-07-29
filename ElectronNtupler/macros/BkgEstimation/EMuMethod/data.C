#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void data() {

  TFile f1("/tmp/rchawla/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Data/MuEG_2015.root");
  TTree *T1 = (TTree*)f1.Get("ntupler/ElectronTree");

  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *rapElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
  vector<float>   *chargeElec;
  vector<float>   *etaSC;
  vector<int>     *passVetoId;
  vector<float>   *ptMuon;
  vector<float>   *etaMuon;
  vector<float>   *phiMuon;
  vector<float>   *energyMuon;
  vector<float>   *chargeMuon;
  vector<float>   *isoPFMuon;
  vector<bool>    *isTightMuon;
  vector<int>     *passMediumId;
  bool            Mu8_Ele17;
  Bool_t          Ele23_WPLoose;

  ptElec = 0;
  etaElec = 0;
  rapElec = 0;
  phiElec = 0;
  energyElec = 0;
  chargeElec = 0;
  etaSC = 0;
  passVetoId = 0;
  passMediumId = 0;
  ptMuon = 0;
  etaMuon = 0;
  phiMuon = 0;
  energyMuon = 0;
  chargeMuon = 0;
  isoPFMuon = 0;
  isTightMuon = 0;

  T1->SetBranchStatus("*",0);
  T1->SetBranchStatus("ptElec", 1);
  T1->SetBranchStatus("etaElec", 1);
  T1->SetBranchStatus("phiElec", 1);
  T1->SetBranchStatus("energyElec", 1);
  T1->SetBranchStatus("etaSC", 1);
  T1->SetBranchStatus("passMediumId", 1);
  T1->SetBranchStatus("chargeElec", 1);
  T1->SetBranchStatus("ptMuon", 1);
  T1->SetBranchStatus("etaMuon", 1);
  T1->SetBranchStatus("phiMuon", 1);
  T1->SetBranchStatus("energyMuon", 1);
  T1->SetBranchStatus("chargeMuon", 1);
  T1->SetBranchStatus("isoPFMuon", 1);
  T1->SetBranchStatus("isTightMuon", 1);
  T1->SetBranchStatus("Mu8_Ele17", 1);
  T1->SetBranchStatus("Ele23_WPLoose", 1);

  T1->SetBranchAddress("ptElec", &ptElec);
  T1->SetBranchAddress("etaElec", &etaElec);
  T1->SetBranchAddress("phiElec", &phiElec);
  T1->SetBranchAddress("energyElec", &energyElec);
  T1->SetBranchAddress("chargeElec", &chargeElec);
  T1->SetBranchAddress("etaSC", &etaSC);
  T1->SetBranchAddress("passMediumId", &passMediumId);
  T1->SetBranchAddress("ptMuon", &ptMuon);
  T1->SetBranchAddress("etaMuon", &etaMuon);
  T1->SetBranchAddress("phiMuon", &phiMuon);
  T1->SetBranchAddress("energyMuon", &energyMuon);
  T1->SetBranchAddress("chargeMuon", &chargeMuon);
  T1->SetBranchAddress("isoPFMuon", &isoPFMuon);
  T1->SetBranchAddress("isTightMuon", &isTightMuon);
  T1->SetBranchAddress("Mu8_Ele17", &Mu8_Ele17);
  T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);

  TFile *file = new TFile("muonEG.root", "recreate");
  TTree *tree = new TTree("tree"," after preselections tree");

  int mediumId;
  TLorentzVector ele1,ele2,diMass;
  TLorentzVector emu1,emu2,diEMuMass;
  TLorentzVector wjet1,wjet2,diWjetMass;

  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;
  vector <double> newmuonPt; vector <double> newmuonEta; vector <double> newmuonEnr; vector <double> newmuonPhi; vector <double> newmuonCharge;

  // Branch variable declaration
  double ElePT, EleEta, ElePhi, MuonPT, MuonEta, MuonPhi;
  double EEMass, EMuMass, WjetMass;

  // Branch declaration
  tree->Branch("ElePT", &ElePT, "ElePT/D");
  tree->Branch("EleEta", &EleEta, "EleEta/D");
  tree->Branch("ElePhi", &ElePhi, "ElePhi/D");
  tree->Branch("MuonPT", &MuonPT, "MuonPT/D");
  tree->Branch("MuonEta", &MuonEta, "MuonEta/D");
  tree->Branch("MuonPhi", &MuonPhi, "MuonPhi/D");
  tree->Branch("EEMass", &EEMass, "EEMass/D");
  tree->Branch("EMuMass", &EMuMass, "EMuMass/D");
  tree->Branch("WjetMass", &WjetMass, "WjetMass/D");

  int nentries = T1->GetEntries();
  //int nentries = 100000;
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int jentry=0; jentry < nentries; jentry++) {
    T1->GetEntry(jentry);

    ElePT = -999.; EleEta = -999.; ElePhi = -999.; MuonPT = -999.; MuonEta = -999.; MuonPhi = -999.;
    EEMass = -999.; EMuMass = -999.; WjetMass = -999.;

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    mediumId = 0;
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();
    newmuonPt.clear(); newmuonEta.clear(); newmuonEnr.clear(); newmuonPhi.clear(); newmuonCharge.clear();

    // Sorting for Electron
    int index1[ptElec->size()];
    float pt1[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++) {
      pt1[el]=ptElec->at(el); }

    int size1 = sizeof(pt1)/sizeof(pt1[0]);
    TMath::Sort(size1,pt1,index1,true);

    // Sorting for muons
    int index2[ptMuon->size()];
    float pt2[ptMuon->size()];

    for(unsigned int mu=0; mu<ptMuon->size(); mu++) {
      pt2[mu]=ptMuon->at(mu); }

    int size2 = sizeof(pt2)/sizeof(pt2[0]);
    TMath::Sort(size2,pt2,index2,true);

    for(int i=0;i<ptElec->size();i++){

      mediumId = passMediumId->at(index1[i]);

      if(mediumId) {
	if(fabs(etaSC->at(index1[i])) < 2.5 && !(fabs(etaSC->at(index1[i])) > 1.4442 && fabs(etaSC->at(index1[i])) < 1.566)){

	  newelePt.push_back(ptElec->at(index1[i]));
	  neweleEta.push_back(etaElec->at(index1[i]));
	  neweleEnr.push_back(energyElec->at(index1[i]));
	  newelePhi.push_back(phiElec->at(index1[i]));
	  neweleCharge.push_back(chargeElec->at(index1[i]));

	} // eta
      } // ID
    } // size

    for(int j=0;j<ptMuon->size();j++){

      if(fabs(etaMuon->at(index2[j])) < 2.4){
	if(isTightMuon->at(index2[j])){
	  if(isoPFMuon->at(index2[j]) < 0.15){

	    newmuonPt.push_back(ptMuon->at(index2[j]));
	    newmuonEta.push_back(etaMuon->at(index2[j]));
	    newmuonPhi.push_back(phiMuon->at(index2[j]));
	    newmuonEnr.push_back(energyMuon->at(index2[j]));
	    newmuonCharge.push_back(chargeMuon->at(index2[j]));

	  } // Isolation
	} // ID
      } // eta
    } // size

    // EE Mass
    if(Ele23_WPLoose){
      if(newelePt.size()==2){
	if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){
	  if(neweleCharge.at(0)*neweleCharge.at(1) == -1){

	    ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    diMass=ele1+ele2;
	    EEMass = diMass.M();

	  } // opposite charge
	} // pt cut
      } // only two electrons
    }

    if(Mu8_Ele17){
      if(newelePt.size()==1 && newmuonPt.size()==1){
	if(newelePt.at(0) > 25. && newmuonPt.at(0) > 15.){
	  if(neweleCharge.at(0)*newmuonCharge.at(0) == -1){

	    ElePT = newelePt.at(0);
	    EleEta = neweleEta.at(0);
	    ElePhi = newelePhi.at(0);

	    MuonPT = newmuonPt.at(0);
	    MuonEta = newmuonEta.at(0);
	    MuonPhi = newmuonPhi.at(0);

	    emu1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    emu2.SetPtEtaPhiE(newmuonPt.at(0),newmuonEta.at(0),newmuonPhi.at(0),newmuonEnr.at(0));

	    diEMuMass=emu1+emu2;
	    EMuMass = diEMuMass.M();

	  } // opposite charge
	} // pt cut
      } // one ele and one muon
    }

    if(Mu8_Ele17){
      if(newelePt.size()==1 && newmuonPt.size()==1){
	if(newelePt.at(0) > 25. && newmuonPt.at(0) > 15.){
	  if(neweleCharge.at(0)*newmuonCharge.at(0) == 1){

	    wjet1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    wjet2.SetPtEtaPhiE(newmuonPt.at(0),newmuonEta.at(0),newmuonPhi.at(0),newmuonEnr.at(0));

	    diWjetMass=wjet1+wjet2;
	    WjetMass = diWjetMass.M();

	  } // opposite charge
	} // pt cut
      } // one ele and one muon
    }

    tree->Fill();

  } // event

  file->Write();
  file->Close();
} 
