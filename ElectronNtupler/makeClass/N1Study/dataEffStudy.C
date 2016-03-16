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

void dataEffStudy() {
  TFile f1("/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/nTuples/analysis_4March/SE_Run2015.root");
  //TFile f1("../../data_check.root");


  TTree *T1 = (TTree*)f1.Get("ntupler/ElectronTree");

  Bool_t          Ele23_WPLoose;
  vector<double>  *pt_Ele23;
  vector<double>  *eta_Ele23;
  vector<double>  *phi_Ele23;
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

  pt_Ele23 = 0;
  eta_Ele23 = 0;
  phi_Ele23 = 0;
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

  T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
  T1->SetBranchAddress("pt_Ele23", &pt_Ele23);
  T1->SetBranchAddress("eta_Ele23", &eta_Ele23);
  T1->SetBranchAddress("phi_Ele23", &phi_Ele23);
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
  T1->SetBranchAddress("etaSC", &etaSC);

  TFile *file = new TFile("single_electron.root", "recreate");

  bool trigMatch_lead, trigMatch_slead;
  double dR, dR1, dR2;
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

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  //TH1I *singleEle = new TH1I("singleEle", "singleEle", );

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

  int nentries = T1->GetEntries();
  //int nentries = 5000;
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int i=1; i <= nentries; i++) {
    T1->GetEntry(i);

    if(i%1000000 == 0){
      cout << "Events Processed :  " << i << endl;
    }

    int index1[ptElec->size()];
    float pt1[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++)
    {
      pt1[el]=ptElec->at(el);
    }

    int size1 = sizeof(pt1)/sizeof(pt1[0]);
    TMath::Sort(size1,pt1,index1,true);

    int index2[pt_Ele23->size()];
    float pt2[pt_Ele23->size()];

    for(unsigned int trig=0; trig<pt_Ele23->size(); trig++)
    {
      pt2[trig]=pt_Ele23->at(trig);
    }

    int size2 = sizeof(pt2)/sizeof(pt2[0]);
    TMath::Sort(size2,pt2,index2,true);

    trigMatch_lead = false; trigMatch_slead = false;
    dR = 0.; dR1 = 0.; dR2 = 0.;
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleMediumId.clear(); neweleHEEPId.clear();
    neweleMediumNoDEta.clear(); neweleMediumNoDPhi.clear(); neweleMediumNoSigmaEtaEta.clear(); neweleMediumNoHOverE.clear();
    neweleMediumNoDxy.clear(); neweleMediumNoDz.clear(); neweleMediumNoEInvP.clear(); neweleMediumNoPFIso.clear();
    neweleMediumNoConVeto.clear(); neweleMediumNoMissHits.clear();

    if(!Ele23_WPLoose) continue;
    if(!Ele23_WPLoose) cout<<"Wrong trigger result"<<endl;

    if(ptElec->size()>=2.){
      for(unsigned int j=0;j<ptElec->size();j++){
	if(fabs(etaSC->at(index1[j])) < 2.5 && !(fabs(etaSC->at(index1[j])) > 1.4442 && fabs(etaSC->at(index1[j])) < 1.566)){

	  newelePt.push_back(ptElec->at(index1[j]));
	  neweleEta.push_back(etaElec->at(index1[j]));
	  neweleEnr.push_back(energyElec->at(index1[j]));
	  newelePhi.push_back(phiElec->at(index1[j]));
	  neweleMediumId.push_back(passMediumId->at(index1[j]));
	  neweleHEEPId.push_back(passHEEPId->at(index1[j]));

	  neweleMediumNoDEta.push_back(isPassMedium_NoDEta->at(index1[j]));
	  neweleMediumNoDPhi.push_back(isPassMedium_NoDPhi->at(index1[j]));
	  neweleMediumNoSigmaEtaEta.push_back(isPassMedium_NoSigmaEtaEta->at(index1[j]));
	  neweleMediumNoHOverE.push_back(isPassMedium_NoHOverE->at(index1[j]));
	  neweleMediumNoDxy.push_back(isPassMedium_NoDxy->at(index1[j]));
	  neweleMediumNoDz.push_back(isPassMedium_NoDz->at(index1[j]));
	  neweleMediumNoEInvP.push_back(isPassMedium_NoEInvP->at(index1[j]));
	  neweleMediumNoPFIso.push_back(isPassMedium_NoPFIso->at(index1[j]));
	  neweleMediumNoConVeto.push_back(isPassMedium_NoConVeto->at(index1[j]));
	  neweleMediumNoMissHits.push_back(isPassMedium_NoMissHits->at(index1[j]));

	} // eta
      } // ptElec size
    } // electrons >=2

    if(newelePt.size()==2){

      //cout<<"Medium ID: "<<neweleMediumId.at(0)<<"   "<<neweleMediumId.at(1)<<endl;
      //cout<<"entry: "<<i<<"   "<<"pt: "<<newelePt.at(0)<<"   "<<newelePt.at(1)<<"   "<<"eta: "<<neweleEta.at(0)<<"   "<<neweleEta.at(1)<<"   "<<"phi: "<<newelePhi.at(0)<<"   "<<newelePhi.at(1)<<endl;

      for(unsigned int k = 0; k < pt_Ele23->size(); k++){
	//cout<<"entry: "<<i<<"   "<<"Trig pt: "<<pt_Ele23->at(k)<<"   "<<"eta: "<<eta_Ele23->at(k)<<"   "<<"Trig phi: "<<phi_Ele23->at(k)<<endl;
	double dR1_comp = 1000.;
	double dR2_comp = 1000.;
	
	dR1 = deltaR(neweleEta[0], newelePhi[0], eta_Ele23->at(k), phi_Ele23->at(k));
	dR2 = deltaR(neweleEta[1], newelePhi[1], eta_Ele23->at(k), phi_Ele23->at(k));

	if(dR1 < 0.1){
	  if (dR1 < dR1_comp)
	  {
	    dR1_comp = dR1;
	    trigMatch_lead = true;
	    //cout<<"dR1: "<<dR1<<endl;
	  }
	} // dR1

	if(dR2 < 0.1){
	  if (dR2 < dR2_comp)
	  {
	    dR2_comp = dR2;
	    trigMatch_slead = true; 
	    //cout<<"dR2: "<<dR2<<endl;
	  }
	} // dR2
      } // pt_Ele23 

      if(trigMatch_lead || trigMatch_slead){
	
	if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){

	  denEle1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  denEle2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	  denDiele=denEle1+denEle2;
	  denMass->Fill(denDiele.M());

	  if(neweleMediumId.at(0)==1 && neweleMediumId.at(1)==1){
	    numEle1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDiele=numEle1+numEle2;
	    numMass->Fill(numDiele.M());
	  }

	  if(neweleHEEPId.at(0)==1 && neweleHEEPId.at(1)==1){
	    numEle1HEEP.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2HEEP.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleHEEP=numEle1HEEP+numEle2HEEP;
	    numMassHEEP->Fill(numDiele.M());
	  }

	  if(neweleMediumNoDEta.at(0)==1 && neweleMediumNoDEta.at(1)==1){
	    numEle1NoDEta.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoDEta.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoDEta=numEle1NoDEta+numEle2NoDEta;
	    numMassNoDEta->Fill(numDieleNoDEta.M());
	  }

	  if(neweleMediumNoDPhi.at(0)==1 && neweleMediumNoDPhi.at(1)==1){
	    numEle1NoDPhi.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoDPhi.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoDPhi=numEle1NoDPhi+numEle2NoDPhi;
	    numMassNoDPhi->Fill(numDieleNoDPhi.M());
	  }

	  if(neweleMediumNoSigmaEtaEta.at(0)==1 && neweleMediumNoSigmaEtaEta.at(1)==1){
	    numEle1NoSigmaEtaEta.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoSigmaEtaEta.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoSigmaEtaEta=numEle1NoSigmaEtaEta+numEle2NoSigmaEtaEta;
	    numMassNoSigmaEtaEta->Fill(numDieleNoSigmaEtaEta.M());
	  }

	  if(neweleMediumNoHOverE.at(0)==1 && neweleMediumNoHOverE.at(1)==1){
	    numEle1NoHOverE.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoHOverE.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoHOverE=numEle1NoHOverE+numEle2NoHOverE;
	    numMassNoHOverE->Fill(numDieleNoHOverE.M());
	  }

	  if(neweleMediumNoDxy.at(0)==1 && neweleMediumNoDxy.at(1)==1){
	    numEle1NoDxy.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoDxy.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoDxy=numEle1NoDxy+numEle2NoDxy;
	    numMassNoDxy->Fill(numDieleNoDxy.M());
	  }

	  if(neweleMediumNoDz.at(0)==1 && neweleMediumNoDz.at(1)==1){
	    numEle1NoDz.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoDz.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoDz=numEle1+numEle2;
	    numMassNoDz->Fill(numDieleNoDz.M());
	  }

	  if(neweleMediumNoEInvP.at(0)==1 && neweleMediumNoEInvP.at(1)==1){
	    numEle1NoEInvP.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoEInvP.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoEInvP=numEle1NoEInvP+numEle2NoEInvP;
	    numMassNoEInvP->Fill(numDieleNoDz.M());
	  }

	  if(neweleMediumNoPFIso.at(0)==1 && neweleMediumNoPFIso.at(1)==1){
	    numEle1NoPFIso.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoPFIso.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoPFIso=numEle1NoPFIso+numEle2NoPFIso;
	    numMassNoPFIso->Fill(numDieleNoPFIso.M());
	  }

	  if(neweleMediumNoConVeto.at(0)==1 && neweleMediumNoConVeto.at(1)==1){
	    numEle1NoConVeto.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoConVeto.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoConVeto=numEle1NoConVeto+numEle2NoConVeto;
	    numMassNoConVeto->Fill(numDieleNoConVeto.M());
	  }

	  if(neweleMediumNoMissHits.at(0)==1 && neweleMediumNoMissHits.at(1)==1){
	    numEle1NoMissHits.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	    numEle2NoMissHits.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	    numDieleNoMissHits=numEle1NoMissHits+numEle2NoMissHits;
	    numMassNoMissHits->Fill(numDieleNoMissHits.M());
	  }

	} // pt
      } // trigger match
    } // exactly two electrons

  } //event loop

  file->Write();
  file->Close();
}
