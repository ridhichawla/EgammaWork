#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

Double_t deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t dEta = deltaEta(eta1, eta2);
  Double_t dPhi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(dEta*dEta + dPhi*dPhi);
  return dr;
}

void EffStudy() {
  TFile f1("/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/nTuples/analysis_4March/DYtoEE_M50.root");

  TTree *T1 = (TTree*)f1.Get("ntupler/ElectronTree");

  vector<float>   *genPostFSR_Pt;
  vector<float>   *genPostFSR_Eta;
  vector<float>   *genPostFSR_Rap;
  vector<float>   *genPostFSR_Phi;
  vector<float>   *genPostFSR_En;
  vector<int>     *passVetoId;
  vector<int>     *passLooseId;
  vector<int>     *passMediumId;
  vector<int>     *passTightId;
  vector<int>     *passHEEPId;
  vector<int>     *isPassMedium_NoPt;
  vector<int>     *isPassMedium_NoScEta;
  vector<int>     *isPassMedium_NoDEta;
  vector<int>     *isPassMedium_NoDPhi;
  vector<int>     *isPassMedium_NoSigmaEtaEta;
  vector<int>     *isPassMedium_NoHOverE;
  vector<int>     *isPassMedium_NoDxy;
  vector<int>     *isPassMedium_NoDz;
  vector<int>     *isPassMedium_NoEInvP;
  vector<int>     *isPassMedium_NoPFIso;
  vector<int>     *isPassMedium_NoConVeto;
  vector<int>     *isPassMedium_NoMissHits;
  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *rapElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
  vector<float>   *etaSC;
  Int_t tauFlag; 

  genPostFSR_Pt = 0;
  genPostFSR_Eta = 0;
  genPostFSR_Rap = 0;
  genPostFSR_Phi = 0;
  genPostFSR_En = 0;
  ptElec = 0;
  etaElec = 0;
  rapElec = 0;
  phiElec = 0;
  energyElec = 0;
  etaSC = 0;
  passVetoId = 0;
  passLooseId = 0;
  passMediumId = 0;
  passTightId = 0;
  passHEEPId = 0;
  isPassMedium_NoPt = 0;
  isPassMedium_NoScEta = 0;
  isPassMedium_NoDEta = 0;
  isPassMedium_NoDPhi = 0;
  isPassMedium_NoSigmaEtaEta = 0;
  isPassMedium_NoHOverE = 0;
  isPassMedium_NoDxy = 0;
  isPassMedium_NoDz = 0;
  isPassMedium_NoEInvP = 0;
  isPassMedium_NoPFIso = 0;
  isPassMedium_NoConVeto = 0;
  isPassMedium_NoMissHits = 0;

  T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
  T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
  T1->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap);
  T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
  T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
  T1->SetBranchAddress("ptElec", &ptElec);
  T1->SetBranchAddress("etaElec", &etaElec);
  T1->SetBranchAddress("rapElec", &rapElec);
  T1->SetBranchAddress("phiElec", &phiElec);
  T1->SetBranchAddress("energyElec", &energyElec);
  T1->SetBranchAddress("passVetoId", &passVetoId);
  T1->SetBranchAddress("passLooseId", &passLooseId);
  T1->SetBranchAddress("passMediumId", &passMediumId);
  T1->SetBranchAddress("passTightId", &passTightId);
  T1->SetBranchAddress("passHEEPId", &passHEEPId);
  T1->SetBranchAddress("isPassMedium_NoPt", &isPassMedium_NoPt);
  T1->SetBranchAddress("isPassMedium_NoScEta", &isPassMedium_NoScEta);
  T1->SetBranchAddress("isPassMedium_NoDEta", &isPassMedium_NoDEta);
  T1->SetBranchAddress("isPassMedium_NoDPhi", &isPassMedium_NoDPhi);
  T1->SetBranchAddress("isPassMedium_NoSigmaEtaEta", &isPassMedium_NoSigmaEtaEta);
  T1->SetBranchAddress("isPassMedium_NoHOverE", &isPassMedium_NoHOverE);
  T1->SetBranchAddress("isPassMedium_NoDxy", &isPassMedium_NoDxy);
  T1->SetBranchAddress("isPassMedium_NoDz", &isPassMedium_NoDz);
  T1->SetBranchAddress("isPassMedium_NoEInvP", &isPassMedium_NoEInvP);
  T1->SetBranchAddress("isPassMedium_NoPFIso", &isPassMedium_NoPFIso);
  T1->SetBranchAddress("isPassMedium_NoConVeto", &isPassMedium_NoConVeto);
  T1->SetBranchAddress("isPassMedium_NoMissHits", &isPassMedium_NoMissHits);
  T1->SetBranchAddress("tauFlag", &tauFlag);
  T1->SetBranchAddress("etaSC", &etaSC);

  TFile *file = new TFile("DY_50toinf.root", "recreate");

  double dR, massGen;
  TLorentzVector gen1,gen2,diGen;
  TLorentzVector numEle1, numEle2, numDiele, numEle1HEEP, numEle2HEEP, numDieleHEEP, denEle1, denEle2, denDiele;
  TLorentzVector numDieleNoDEta, numEle1NoDEta, numEle2NoDEta;
  TLorentzVector numDieleNoDPhi, numEle1NoDPhi, numEle2NoDPhi;
  TLorentzVector numDieleNoSigmaEtaEta, numEle1NoSigmaEtaEta, numEle2NoSigmaEtaEta;
  TLorentzVector numDieleNoHOverE, numEle1NoHOverE, numEle2NoHOverE;
  TLorentzVector numDieleNoDxy, numEle1NoDxy, numEle2NoDxy;
  TLorentzVector numDieleNoDz, numEle1NoDz, numEle2NoDz;
  TLorentzVector numDieleNoEInvP, numEle1NoEInvP, numEle2NoEInvP;
  TLorentzVector numDieleNoPFIso, numEle1NoPFIso, numEle2NoPFIso;
  TLorentzVector numDieleNoConVeto, numEle1NoConVeto, numEle2NoConVeto;
  TLorentzVector numDieleNoMissHits, numEle1NoMissHits, numEle2NoMissHits;

  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi;
  vector <double> neweleMediumId; vector <double> neweleHEEPId; 
  vector<double> neweleMediumNoDEta; vector<double> neweleMediumNoDPhi; vector<double> neweleMediumNoSigmaEtaEta;
  vector<double> neweleMediumNoHOverE; vector<double> neweleMediumNoDxy; vector<double> neweleMediumNoDz; vector<double> neweleMediumNoEInvP;
  vector<double> neweleMediumNoPFIso; vector<double> neweleMediumNoConVeto; vector<double> neweleMediumNoMissHits;
  
  vector <double> newgenPt; vector <double> newgenEta; vector <double> newgenEnr; vector <double> newgenPhi;
  vector <double> elePt; vector <double> eleEta; vector <double> eleEnr; vector <double> elePhi;
  vector <double> eleMediumId; vector <double> eleHEEPId;
  vector <double> eleMediumNoDEta; vector <double> eleMediumNoDPhi; vector <double> eleMediumNoSigmaEtaEta; vector <double> eleMediumNoHOverE;
  vector <double> eleMediumNoDxy; vector <double> eleMediumNoDz; vector <double> eleMediumNoEInvP; vector <double> eleMediumNoPFIso;
  vector <double> eleMediumNoConVeto; vector <double> eleMediumNoMissHits;

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};
  
  TH1D *numMass = new TH1D("numMass", "numMass", 45, xbins);
  TH1D *numMassHEEP = new TH1D("numMassHEEP", "numMassHEEP", 45, xbins);
  TH1D *denMass = new TH1D("denMass", "denMass", 45, xbins);
  
  TH1D *numMassNoDEta        = new TH1D("numMassNoDEta", "numMassNoDEta", 45, xbins);
  TH1D *numMassNoDPhi        = new TH1D("numMassNoDPhi", "numMassNoDPhi", 45, xbins);
  TH1D *numMassNoSigmaEtaEta = new TH1D("numMassNoSigmaEtaEta", "numMassNoSigmaEtaEta", 45, xbins);
  TH1D *numMassNoHOverE      = new TH1D("numMassNoHOverE", "numMassNoHOverE", 45, xbins);
  TH1D *numMassNoDxy         = new TH1D("numMassNoDxy", "numMassNoDxy", 45, xbins);
  TH1D *numMassNoDz          = new TH1D("numMassNoDz", "numMassNoDz", 45, xbins);
  TH1D *numMassNoEInvP       = new TH1D("numMassNoEInvP", "numMassNoEInvP", 45, xbins);
  TH1D *numMassNoPFIso       = new TH1D("numMassNoPFIso", "numMassNoPFIso", 45, xbins);
  TH1D *numMassNoConVeto     = new TH1D("numMassNoConVeto", "numMassNoConVeto", 45, xbins);
  TH1D *numMassNoMissHits    = new TH1D("numMassNoMissHits", "numMassNoMissHits", 45, xbins);

  numMass->Sumw2(); numMassHEEP->Sumw2(); denMass->Sumw2();
  numMassNoDEta->Sumw2(); numMassNoDPhi->Sumw2(); numMassNoSigmaEtaEta->Sumw2(); numMassNoHOverE->Sumw2(); numMassNoDxy->Sumw2();
  numMassNoDz->Sumw2(); numMassNoEInvP->Sumw2(); numMassNoPFIso->Sumw2(); numMassNoConVeto->Sumw2(); numMassNoMissHits->Sumw2();

  //double xsec = 6024.;
  //double procEvents = 451498239235.;
  //double lumiWeight = 6024./451498239235.;

  int nentries = T1->GetEntries();
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int j=0; j < nentries; j++) {
    T1->GetEntry(j);

    if(j%1000000 == 0){
      cout << "Events Processed :  " << j << endl;
    }
    
    int index1[ptElec->size()];
    float pt1[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++)
    {
      pt1[el]=ptElec->at(el);
    }

    int sizer = sizeof(pt1)/sizeof(pt1[0]);
    TMath::Sort(sizer,pt1,index1,true);

    int index2[genPostFSR_Pt->size()];
    float pt2[genPostFSR_Pt->size()];

    for(unsigned int b=0; b<genPostFSR_Pt->size(); b++)
    {
      pt2[b]=genPostFSR_Pt->at(b);
    }

    int sizen = sizeof(pt2)/sizeof(pt2[0]);
    TMath::Sort(sizen,pt2,index2,true);

    dR = 0.; massGen = 0.0;
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleMediumId.clear(); neweleHEEPId.clear();
    neweleMediumNoDEta.clear(); neweleMediumNoDPhi.clear(); neweleMediumNoSigmaEtaEta.clear(); neweleMediumNoHOverE.clear();
    neweleMediumNoDxy.clear(); neweleMediumNoDz.clear(); neweleMediumNoEInvP.clear(); neweleMediumNoPFIso.clear();
    neweleMediumNoConVeto.clear(); neweleMediumNoMissHits.clear();

    newgenPt.clear(); newgenEta.clear(); newgenEnr.clear(); newgenPhi.clear();
    elePt.clear(); eleEta.clear(); eleEnr.clear(); elePhi.clear(); eleMediumId.clear(); eleHEEPId.clear();
    eleMediumNoDEta.clear(); eleMediumNoDPhi.clear(); eleMediumNoSigmaEtaEta.clear(); eleMediumNoHOverE.clear();
    eleMediumNoDxy.clear(); eleMediumNoDz.clear(); eleMediumNoEInvP.clear(); eleMediumNoPFIso.clear();
    eleMediumNoConVeto.clear(); eleMediumNoMissHits.clear();

    unsigned int matchGen = 0;
    unsigned int matchReco = 0;

    if(genPostFSR_Pt->size() > 2) cout<<"size > 2"<<endl;

    if(genPostFSR_Pt->size() == 2){
      gen1.SetPtEtaPhiE(genPostFSR_Pt->at(0),genPostFSR_Eta->at(0),genPostFSR_Phi->at(0),genPostFSR_En->at(0));
      gen2.SetPtEtaPhiE(genPostFSR_Pt->at(1),genPostFSR_Eta->at(1),genPostFSR_Phi->at(1),genPostFSR_En->at(1));

      diGen=gen1+gen2;
      massGen=diGen.M();

    }

    if(massGen > 200.) continue;
    if(!tauFlag){

      if(ptElec->size()>=2.){
	for(unsigned int i=0;i<ptElec->size();i++){
	  if(fabs(etaSC->at(index1[i])) < 2.5 && !(fabs(etaSC->at(index1[i])) > 1.4442 && fabs(etaSC->at(index1[i])) < 1.566)){

	    newelePt.push_back(ptElec->at(index1[i]));
	    neweleEta.push_back(etaElec->at(index1[i]));
	    neweleEnr.push_back(energyElec->at(index1[i]));
	    newelePhi.push_back(phiElec->at(index1[i]));
	    neweleMediumId.push_back(passMediumId->at(index1[i]));
	    neweleHEEPId.push_back(passHEEPId->at(index1[i]));

	    neweleMediumNoDEta.push_back(isPassMedium_NoDEta->at(index1[i]));
	    neweleMediumNoDPhi.push_back(isPassMedium_NoDPhi->at(index1[i]));
	    neweleMediumNoSigmaEtaEta.push_back(isPassMedium_NoSigmaEtaEta->at(index1[i]));
	    neweleMediumNoHOverE.push_back(isPassMedium_NoHOverE->at(index1[i]));
	    neweleMediumNoDxy.push_back(isPassMedium_NoDxy->at(index1[i]));
	    neweleMediumNoDz.push_back(isPassMedium_NoDz->at(index1[i]));
	    neweleMediumNoEInvP.push_back(isPassMedium_NoEInvP->at(index1[i]));
	    neweleMediumNoPFIso.push_back(isPassMedium_NoPFIso->at(index1[i]));
	    neweleMediumNoConVeto.push_back(isPassMedium_NoConVeto->at(index1[i]));
	    neweleMediumNoMissHits.push_back(isPassMedium_NoMissHits->at(index1[i]));

	  } // eta
	} //nEle
      }
      if(genPostFSR_Pt->size() >= 2.){
	for(unsigned int j=0;j<genPostFSR_Eta->size();j++){
	  if(fabs(genPostFSR_Eta->at(index2[j])) < 2.5 && !(fabs(genPostFSR_Eta->at(index2[j])) > 1.4442 && fabs(genPostFSR_Eta->at(index2[j])) < 1.566)){

	    newgenPt.push_back(genPostFSR_Pt->at(index2[j]));
	    newgenEta.push_back(genPostFSR_Eta->at(index2[j]));
	    newgenEnr.push_back(genPostFSR_En->at(index2[j]));
	    newgenPhi.push_back(genPostFSR_Phi->at(index2[j]));
	  }
	}
      }

      for(unsigned int ireco = 0; ireco < neweleEta.size(); ireco++){
	double dR_comp_post = 1000.;
	for(unsigned int igen = 0; igen < newgenEta.size(); igen++){
	  dR = deltaR(neweleEta[ireco], newelePhi[ireco], newgenEta[igen], newgenPhi[igen]);
	  //cout<<"dR: "<<dR<<endl;
	  if(dR < 0.1){
	    if (dR < dR_comp_post)
	    {
	      dR_comp_post = dR;
	      matchGen = igen ; matchReco = ireco ; //iter[ireco]=igen;
	    }
	    elePt.push_back(newelePt[matchReco]);
	    eleEta.push_back(neweleEta[matchReco]);
	    eleEnr.push_back(neweleEnr[matchReco]);
	    elePhi.push_back(newelePhi[matchReco]);
	    eleMediumId.push_back(neweleMediumId[matchReco]);
	    eleHEEPId.push_back(neweleHEEPId[matchReco]);

	    eleMediumNoDEta.push_back(neweleMediumNoDEta[matchReco]);
	    eleMediumNoDPhi.push_back(neweleMediumNoDPhi[matchReco]);
	    eleMediumNoSigmaEtaEta.push_back(neweleMediumNoSigmaEtaEta[matchReco]);
	    eleMediumNoHOverE.push_back(neweleMediumNoHOverE[matchReco]);
	    eleMediumNoDxy.push_back(neweleMediumNoDxy[matchReco]);
	    eleMediumNoDz.push_back(neweleMediumNoDz[matchReco]);
	    eleMediumNoEInvP.push_back(neweleMediumNoEInvP[matchReco]);
	    eleMediumNoPFIso.push_back(neweleMediumNoPFIso[matchReco]);
	    eleMediumNoConVeto.push_back(neweleMediumNoConVeto[matchReco]);
	    eleMediumNoMissHits.push_back(neweleMediumNoMissHits[matchReco]);
	  }
	}
      }          

      if(elePt.size()==2){
	if(elePt.at(0) > 30. && elePt.at(1) > 10.){

	  denEle1.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	  denEle2.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	  denDiele=denEle1+denEle2;
	  denMass->Fill(denDiele.M());

	  if(eleMediumId.at(0)==1 && eleMediumId.at(1)==1){
	    numEle1.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDiele=numEle1+numEle2;
	    numMass->Fill(numDiele.M());
	  }

	  if(eleHEEPId.at(0)==1 && eleHEEPId.at(1)==1){
	    numEle1HEEP.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2HEEP.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleHEEP=numEle1HEEP+numEle2HEEP;
	    numMassHEEP->Fill(numDiele.M());
	  }

	  if(eleMediumNoDEta.at(0)==1 && eleMediumNoDEta.at(1)==1){
	    numEle1NoDEta.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoDEta.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoDEta=numEle1NoDEta+numEle2NoDEta;
	    numMassNoDEta->Fill(numDieleNoDEta.M());
	  }

	  if(eleMediumNoDPhi.at(0)==1 && eleMediumNoDPhi.at(1)==1){
	    numEle1NoDPhi.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoDPhi.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoDPhi=numEle1NoDPhi+numEle2NoDPhi;
	    numMassNoDPhi->Fill(numDieleNoDPhi.M());
	  }

	  if(eleMediumNoSigmaEtaEta.at(0)==1 && eleMediumNoSigmaEtaEta.at(1)==1){
	    numEle1NoSigmaEtaEta.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoSigmaEtaEta.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoSigmaEtaEta=numEle1NoSigmaEtaEta+numEle2NoSigmaEtaEta;
	    numMassNoSigmaEtaEta->Fill(numDieleNoSigmaEtaEta.M());
	  }

	  if(eleMediumNoHOverE.at(0)==1 && eleMediumNoHOverE.at(1)==1){
	    numEle1NoHOverE.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoHOverE.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoHOverE=numEle1NoHOverE+numEle2NoHOverE;
	    numMassNoHOverE->Fill(numDieleNoHOverE.M());
	  }

	  if(eleMediumNoDxy.at(0)==1 && eleMediumNoDxy.at(1)==1){
	    numEle1NoDxy.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoDxy.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoDxy=numEle1NoDxy+numEle2NoDxy;
	    numMassNoDxy->Fill(numDieleNoDxy.M());
	  }

	  if(eleMediumNoDz.at(0)==1 && eleMediumNoDz.at(1)==1){
	    numEle1NoDz.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoDz.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoDz=numEle1+numEle2;
	    numMassNoDz->Fill(numDieleNoDz.M());
	  }

	  if(eleMediumNoEInvP.at(0)==1 && eleMediumNoEInvP.at(1)==1){
	    numEle1NoEInvP.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoEInvP.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoEInvP=numEle1NoEInvP+numEle2NoEInvP;
	    numMassNoEInvP->Fill(numDieleNoDz.M());
	  }

	  if(eleMediumNoPFIso.at(0)==1 && eleMediumNoPFIso.at(1)==1){
	    numEle1NoPFIso.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoPFIso.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoPFIso=numEle1NoPFIso+numEle2NoPFIso;
	    numMassNoPFIso->Fill(numDieleNoPFIso.M());
	  }

	  if(eleMediumNoConVeto.at(0)==1 && eleMediumNoConVeto.at(1)==1){
	    numEle1NoConVeto.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoConVeto.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoConVeto=numEle1NoConVeto+numEle2NoConVeto;
	    numMassNoConVeto->Fill(numDieleNoConVeto.M());
	  }

	  if(eleMediumNoMissHits.at(0)==1 && eleMediumNoMissHits.at(1)==1){
	    numEle1NoMissHits.SetPtEtaPhiE(elePt.at(0),eleEta.at(0),elePhi.at(0),eleEnr.at(0));
	    numEle2NoMissHits.SetPtEtaPhiE(elePt.at(1),eleEta.at(1),elePhi.at(1),eleEnr.at(1));

	    numDieleNoMissHits=numEle1NoMissHits+numEle2NoMissHits;
	    numMassNoMissHits->Fill(numDieleNoMissHits.M());
	  }

	} // pt
      } // good electrons
    }  // tauFlag

  } //event loop
  file->Write();
  file->Close();
}
