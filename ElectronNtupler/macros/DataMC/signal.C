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

void signal() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;

  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

  double xsec[11] = {18610./3,5870./3,226./3,7.67/3,0.423/3,0.24/3,0.035/3,0.03/3,0.016/3,0.002/3,0.00054/3};
  double sumofWts[11] = {771413889185.162476,144505031098.323120,219889705.060318,7008766.904321,122987.746342,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

  workdir = "/tmp/rchawla/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";

  InputFiles_signal_DY.clear();

  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_10to50.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_50toInf.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_100to200.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_200to400.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_400to500.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_500to700.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_700to800.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_800to1000.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1000to1500.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1500to2000.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_2000to3000.root"));

  int nsample = InputFiles_signal_DY.size();
  TFile* file[11]; 

  for(unsigned int jentry = 0; jentry < nsample; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("../dataPUDist.root");
    TFile *f2 = TFile::Open("../PileUp_MC.root");

    // data histogram 
    TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1D *MC_puDist = (TH1D*)f2->Get("pileup_MC");
    TH1D *weights = (TH1D*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    vector<float>   *genPostFSR_Pt;
    vector<float>   *genPostFSR_Eta;
    vector<float>   *genPostFSR_Phi;
    vector<float>   *genPostFSR_En;
    vector<float>   *genPreFSR_Pt; 
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<int>     *passMediumId;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etSC;
    vector<float>   *etaSC;
    vector<float>   *full5x5_sigmaIetaIeta;
    vector<float>   *dEtaIn;
    vector<float>   *dPhiIn;
    vector<float>   *hOverE;
    vector<float>   *isoRho;
    vector<float>   *ooEmooP;
    vector<float>   *d0;
    vector<float>   *dz;
    vector<int>     *expectedMissingInnerHits;
    vector<int>     *passConversionVeto;
    Int_t           tauFlag;
    Double_t        theWeight;
    Bool_t          Ele23_WPLoose;
    vector<double>  *pt_Ele23;
    vector<double>  *eta_Ele23;
    vector<double>  *phi_Ele23;
    Int_t           nPV;
    Int_t           nPUTrue;

    genPostFSR_Pt = 0;
    genPostFSR_Eta = 0;
    genPostFSR_Phi = 0;
    genPostFSR_En = 0;
    genPreFSR_Pt = 0; 
    genPreFSR_Eta = 0; 
    genPreFSR_Phi = 0; 
    genPreFSR_En = 0;
    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etSC = 0;
    etaSC = 0;
    passMediumId = 0;
    full5x5_sigmaIetaIeta = 0;
    dEtaIn = 0;
    dPhiIn = 0;
    hOverE = 0;
    isoRho = 0;
    ooEmooP = 0;
    d0 = 0;
    dz = 0;
    expectedMissingInnerHits = 0;
    passConversionVeto = 0;
    pt_Ele23 = 0;
    eta_Ele23 = 0;
    phi_Ele23 = 0;

    T1->SetBranchStatus("*",0);
    T1->SetBranchStatus("genPostFSR_Pt", 1);
    T1->SetBranchStatus("genPostFSR_Eta", 1);
    T1->SetBranchStatus("genPostFSR_Phi", 1);
    T1->SetBranchStatus("genPostFSR_En", 1);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("nPV", 1);
    T1->SetBranchStatus("nPUTrue", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("rapElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etSC", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("full5x5_sigmaIetaIeta", 1);
    T1->SetBranchStatus("dEtaIn", 1);
    T1->SetBranchStatus("dPhiIn", 1);
    T1->SetBranchStatus("hOverE", 1);
    T1->SetBranchStatus("isoRho", 1);
    T1->SetBranchStatus("ooEmooP", 1);
    T1->SetBranchStatus("d0", 1);
    T1->SetBranchStatus("dz", 1);
    T1->SetBranchStatus("expectedMissingInnerHits", 1);
    T1->SetBranchStatus("passConversionVeto", 1);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("pt_Ele23", 1);
    T1->SetBranchStatus("eta_Ele23", 1);
    T1->SetBranchStatus("phi_Ele23", 1);

    T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
    T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
    T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
    T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt); 
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta); 
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi); 
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etSC", &etSC);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta);
    T1->SetBranchAddress("dEtaIn", &dEtaIn);
    T1->SetBranchAddress("dPhiIn", &dPhiIn);
    T1->SetBranchAddress("hOverE", &hOverE);
    T1->SetBranchAddress("isoRho", &isoRho);
    T1->SetBranchAddress("ooEmooP", &ooEmooP);
    T1->SetBranchAddress("d0", &d0);
    T1->SetBranchAddress("dz", &dz);
    T1->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits);
    T1->SetBranchAddress("passConversionVeto", &passConversionVeto);
    T1->SetBranchAddress("passConversionVeto", &passConversionVeto);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("pt_Ele23", &pt_Ele23);
    T1->SetBranchAddress("eta_Ele23", &eta_Ele23);
    T1->SetBranchAddress("phi_Ele23", &phi_Ele23);
    T1->SetBranchAddress("nPV", &nPV);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("ROOTFiles/latest_3Apr/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    int mediumId;
    double dR1, dR2;
    double massGPre;
    TLorentzVector genPre1,genPre2,diGen;
    TLorentzVector gen1_NoAcc,gen2_NoAcc,digen_NoAcc;
    TLorentzVector gen1,gen2,digen;
    TLorentzVector ele1,ele2,dielectron;

    // Branch variable declaration
    double sigmaEta1, dEta1, dPhi1, he1, pfIso1, ep1, D01, Dz1;
    double sigmaEta2, dEta2, dPhi2, he2, pfIso2, ep2, D02, Dz2;
    bool missHits1, conVeto1, missHits2, conVeto2;
    double primVtx;
    double Ele1PT, Ele2PT, Ele1Eta, Ele2Eta, SC1Et, SC2Et, SC1Eta, SC2Eta, Ele1Phi, Ele2Phi, Ele1Enr, Ele2Enr;
    double ZMass, ZPT, ZEta, ZRapidity, ZPhi;
    double lumiWeight, genWeight, PUWeight;
    bool BB, BE, EE;
    bool Ele1Barrel, Ele2Barrel, Ele1Endcap, Ele2Endcap;

    // Branch declaration
    tree->Branch("sigmaEta1", &sigmaEta1, "sigmaEta1/D");
    tree->Branch("sigmaEta2", &sigmaEta2, "sigmaEta2/D");
    tree->Branch("dEta1", &dEta1, "dEta1/D");
    tree->Branch("dEta2", &dEta2, "dEta2/D");
    tree->Branch("dPhi1", &dPhi1, "dPhi1/D");
    tree->Branch("dPhi2", &dPhi2, "dPhi2/D");
    tree->Branch("he1", &he1, "he1/D");
    tree->Branch("he2", &he2, "he2/D");
    tree->Branch("pfIso1", &pfIso1, "pfIso1/D");
    tree->Branch("pfIso2", &pfIso2, "pfIso2/D");
    tree->Branch("ep1", &ep1, "ep1/D");
    tree->Branch("ep2", &ep2, "ep2/D");
    tree->Branch("D01", &D01, "D01/D");
    tree->Branch("D02", &D02, "D02/D");
    tree->Branch("Dz1", &Dz1, "Dz1/D");
    tree->Branch("Dz2", &Dz2, "Dz2/D");
    tree->Branch("missHits1", &missHits1, "missHits1/B");
    tree->Branch("missHits2", &missHits2, "missHits2/B");
    tree->Branch("conVeto1", &conVeto1, "conVeto1/B");
    tree->Branch("conVeto2", &conVeto2, "conVeto2/B");
    tree->Branch("primVtx", &primVtx, "primVtx/D");
    tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D");
    tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");
    tree->Branch("Ele1Eta", &Ele1Eta, "Ele1Eta/D");
    tree->Branch("Ele2Eta", &Ele2Eta, "Ele2Eta/D");
    tree->Branch("SC1Et", &SC1Et, "SC1Et/D");
    tree->Branch("SC2Et", &SC2Et, "SC2Et/D");
    tree->Branch("SC1Eta", &SC1Eta, "SC1Eta/D");
    tree->Branch("SC2Eta", &SC2Eta, "SC2Eta/D");
    tree->Branch("Ele1Phi", &Ele1Phi, "Ele1Phi/D");
    tree->Branch("Ele2Phi", &Ele2Phi, "Ele2Phi/D");
    tree->Branch("Ele1Enr", &Ele1Enr, "Ele1Enr/D");
    tree->Branch("Ele2Enr", &Ele2Enr, "Ele2Enr/D");
    tree->Branch("ZMass", &ZMass, "ZMass/D");
    tree->Branch("ZPT", &ZPT, "ZPT/D");
    tree->Branch("ZEta", &ZEta, "ZEta/D");
    tree->Branch("ZRapidity", &ZRapidity, "ZRapidity/D");
    tree->Branch("ZPhi", &ZPhi, "ZPhi/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/D");
    tree->Branch("BB", &BB, "BB/B");
    tree->Branch("BE", &BE, "BE/B");
    tree->Branch("EE", &EE, "EE/B");
    tree->Branch("Ele1Barrel", &Ele1Barrel, "Ele1Barrel/B");
    tree->Branch("Ele2Barrel", &Ele2Barrel, "Ele2Barrel/B");
    tree->Branch("Ele1Endcap", &Ele1Endcap, "Ele1Endcap/B");
    tree->Branch("Ele2Endcap", &Ele2Endcap, "Ele2Endcap/B");

    vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi;
    vector <double> newscEta; vector <double> newscEt;
    vector <double> neweleSigmaEta; vector <double> neweleDEta; vector <double> neweleDPhi; vector <double> neweleHE; vector <double> newelePFIso;
    vector <double> neweleEP; vector <double> neweleD0; vector <double> neweleDz; vector <double> neweleMissHits; vector <double> neweleConVeto;
    vector <double> newgenPt; vector <double> newgenEta; vector <double> newgenEnr; vector <double> newgenPhi;

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    int nentries = T1->GetEntries();
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      // Sorting Reco level
      int index[ptElec->size()];
      float pt1[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++)
      {
	pt1[el]=ptElec->at(el);
      }

      int size = sizeof(pt1)/sizeof(pt1[0]);
      TMath::Sort(size,pt1,index,true);

      massGPre=0.;
      dR1 = 0.; dR2 = 0.;
      BB=false; BE=false; EE=false;
      Ele1Barrel=false; Ele2Barrel=false; Ele1Endcap=false; Ele2Endcap=false;
      mediumId = 0;
      newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); newscEta.clear(); newscEt.clear();
      neweleSigmaEta.clear(); neweleDEta.clear(); neweleDPhi.clear(); neweleHE.clear(); newelePFIso.clear();
      neweleEP.clear(); neweleD0.clear(); neweleDz.clear(); neweleMissHits.clear(); neweleConVeto.clear();
      newgenPt.clear(); newgenEta.clear(); newgenEnr.clear(); newgenPhi.clear();

      // PU Weight
      int bin = 0;
      double puWeights = 1.0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      PUWeight = weights->GetBinContent(bin);

      if(genPreFSR_Pt->size() == 2){
	genPre1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	genPre2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=genPre1+genPre2;
	massGPre=diGen.M();
      }

      if(jentry==1 && massGPre >= 100.) continue;        // Gen Mass cut ----- for 50 to inf sample
      if(tauFlag) continue;

      primVtx = nPV;

      if(!Ele23_WPLoose) continue;

      for(int j=0;j<ptElec->size();j++){

	mediumId = passMediumId->at(index[j]);

	if(mediumId) {

	  if(fabs(etaSC->at(index[j])) < 2.5 && !(fabs(etaSC->at(index[j])) > 1.4442 && fabs(etaSC->at(index[j])) < 1.566)){

	    newelePt.push_back(ptElec->at(index[j]));
	    neweleEta.push_back(etaElec->at(index[j]));
	    neweleEnr.push_back(energyElec->at(index[j]));
	    newelePhi.push_back(phiElec->at(index[j]));
	    newscEta.push_back(etaSC->at(index[j]));
	    newscEt.push_back(etSC->at(index[j]));

	    neweleSigmaEta.push_back(full5x5_sigmaIetaIeta->at(index[j]));
	    neweleDEta.push_back(dEtaIn->at(index[j]));
	    neweleDPhi.push_back(dPhiIn->at(index[j]));
	    neweleHE.push_back(hOverE->at(index[j]));
	    newelePFIso.push_back(isoRho->at(index[j]));
	    neweleEP.push_back(ooEmooP->at(index[j]));
	    neweleD0.push_back(d0->at(index[j]));
	    neweleDz.push_back(dz->at(index[j]));
	    neweleMissHits.push_back(expectedMissingInnerHits->at(index[j]));
	    neweleConVeto.push_back(passConversionVeto->at(index[j]));

	  } // eta
	} // ID
      } // pt size

      if(newelePt.size()==2){

	int countTrig = 0;

	for(unsigned int j = 0; j < pt_Ele23->size(); j++){

	  dR1 = deltaR(neweleEta.at(0), newelePhi.at(0), eta_Ele23->at(j), phi_Ele23->at(j));
	  dR2 = deltaR(neweleEta.at(1), newelePhi.at(1), eta_Ele23->at(j), phi_Ele23->at(j));

	  if(dR1 < 0.1 || dR2 < 0.1) countTrig++;

	} // pt_Ele23

	if(newelePt.at(0) > 30. && newelePt.at(1) > 10. && countTrig != 0){

	  sigmaEta1 = neweleSigmaEta.at(0);
	  sigmaEta2 = neweleSigmaEta.at(1);
	  dEta1 = neweleDEta.at(0);
	  dEta2 = neweleDEta.at(1);
	  dPhi1 = neweleDPhi.at(0);
	  dPhi2 = neweleDPhi.at(1);
	  he1 = neweleHE.at(0);
	  he2 = neweleHE.at(1);
	  pfIso1 = newelePFIso.at(0);
	  pfIso2 = newelePFIso.at(1);
	  ep1 = neweleEP.at(0);
	  ep2 = neweleEP.at(1);
	  D01 = neweleD0.at(0);
	  D02 = neweleD0.at(1);
	  Dz1 = neweleDz.at(0);
	  Dz2 = neweleDz.at(1);
	  missHits1 = neweleMissHits.at(0);
	  missHits2 = neweleMissHits.at(1);
	  conVeto1 = neweleConVeto.at(0);
	  conVeto2 = neweleConVeto.at(1);

	  Ele1PT  = newelePt.at(0);
	  Ele2PT  = newelePt.at(1);
	  Ele1Eta = neweleEta.at(0);
	  Ele2Eta = neweleEta.at(1);
	  SC1Et   = newscEt.at(0);
	  SC2Et   = newscEt.at(1);
	  SC1Eta  = newscEta.at(0);
	  SC2Eta  = newscEta.at(1);
	  Ele1Phi = newelePhi.at(0);
	  Ele2Phi = newelePhi.at(1);
	  Ele1Enr = neweleEnr.at(0);
	  Ele2Enr = neweleEnr.at(1);

	  ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	  ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));
	  dielectron=ele1+ele2;

	  ZMass = dielectron.M();
	  ZPT = dielectron.Pt();
	  ZEta = dielectron.Eta();
	  ZRapidity = dielectron.Rapidity();
	  ZPhi = dielectron.Phi();

	  if(fabs(newscEta.at(0)) < 1.4442) Ele1Barrel = true;
	  if(fabs(newscEta.at(1)) < 1.4442) Ele2Barrel = true;
	  if(fabs(newscEta.at(0)) > 1.566)  Ele1Endcap = true;
	  if(fabs(newscEta.at(1)) > 1.566)  Ele2Endcap = true;

	  if(fabs(newscEta.at(0)) < 1.4442 && fabs(newscEta.at(1)) < 1.4442) BB = true;
	  if((fabs(newscEta.at(0)) < 1.4442 && fabs(newscEta.at(1)) > 1.566) || (fabs(newscEta.at(0)) > 1.566 && fabs(newscEta.at(1)) < 1.4442)) BE =true;
	  if(fabs(newscEta.at(0)) > 1.566 && fabs(newscEta.at(1)) > 1.566) EE =true;

	  lumiWeight = lumi_Weight;
	  genWeight  = theWeight;

	  tree->Fill();

	} // pt
      } // good electrons
    } // event Loop

    file[jentry]->Write();
    file[jentry]->Close();

    cout<<""<<endl;

  } // file Loop
}
