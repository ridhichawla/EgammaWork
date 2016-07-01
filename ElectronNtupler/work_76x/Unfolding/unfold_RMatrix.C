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

//bool unfold_RMatrix::mySortFnx(Double_t i, Double_t j) { return i>j; }

void unfold_RMatrix() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};
  //int mass[9] = {10,50,120,200,400,800,1400,2300,3500};

  double xsec[11] = {18609.9/3,5789./3,226./3,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
  double sumofWts[11] = {758771667841.549683,150453917158.292908,219884886.214768,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};
  
  //double xsec[8] = {18609.9/3,1975.,19.32,2.731,0.241,0.01678,0.00139,0.00008948};
  //double noEvts[8] = {631905.,2984400.,100000.,100000.,100000.,100000.,100000.,100000.};

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

  /*InputFiles_signal_DY.push_back(TFile::Open(workdir+"MG_M5_50.root"));
  //InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_10to50.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"Powheg_50_120.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"Powheg_120_200.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"Powheg_200_400.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"Powheg_400_800.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"Powheg_800_1400.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"Powheg_1400_2300.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"Powheg_2300_3500.root"));*/

  int nsample = InputFiles_signal_DY.size();
  TFile* file[11];
  //TFile* file[8];

  for(unsigned int jentry = 0; jentry < nsample; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("../dataPUDist.root");
    TFile *f2 = TFile::Open("../PileUp_MC.root");

    // data histogram 
    TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");
    TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
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
    vector<float>   *etaSC;
    Int_t           tauFlag;
    Double_t        theWeight;
    Bool_t          Ele23_WPLoose;
    Int_t           nPUTrue;
    vector<double>  *eta_Ele23;
    vector<double>  *phi_Ele23;

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
    etaSC = 0;
    passMediumId = 0;
    eta_Ele23 = 0;
    phi_Ele23 = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("genPostFSR_Pt", 1);
    T1->SetBranchStatus("genPostFSR_Eta", 1);
    T1->SetBranchStatus("genPostFSR_Phi", 1);
    T1->SetBranchStatus("genPostFSR_En", 1);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("nPUTrue", 1);
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
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);
    T1->SetBranchAddress("eta_Ele23", &eta_Ele23);
    T1->SetBranchAddress("phi_Ele23", &phi_Ele23);

    file[jentry] = new TFile(Form("detRes_Unfold/Systematics/Mass_%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    //double lumi_Weight = xsec[jentry]/noEvts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    int mediumId, dzMaskId, passId;
    bool trigMatch_lead, trigMatch_slead;
    double massGPre;
    double dR, dR1, dR2;
    TLorentzVector ele1,ele2,dielectron;
    TLorentzVector gen1,gen2,digen;
    TLorentzVector genPre1,genPre2,diGen;

    // Branch variable declaration
    double massReco, massGen;
    bool isReco, isGen;
    double lumiWeight, genWeight, PUWeight;
    bool BB, BE, EE;

    // Branch declaration
    tree->Branch("massReco", &massReco, "massReco/D");
    tree->Branch("massGen", &massGen, "massGen/D");
    tree->Branch("isReco", &isReco, "isReco/B");
    tree->Branch("isGen", &isGen, "isGen/B");
    tree->Branch("BB", &BB, "BB/B");
    tree->Branch("BE", &BE, "BE/B");
    tree->Branch("EE", &EE, "EE/B");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/D");

    vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> newscEta;
    vector <double> newgenPt; vector <double> newgenEta; vector <double> newgenEnr; vector <double> newgenPhi;

    int nentries = T1->GetEntries();
    //int nentries = 50000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      // Sorting Reco level
      int index1[ptElec->size()];
      float pt1[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++)
      {
	pt1[el]=ptElec->at(el);
      }

      int size = sizeof(pt1)/sizeof(pt1[0]);
      TMath::Sort(size,pt1,index1,true);

      // Sorting Post Gen level
      int index2[genPostFSR_Pt->size()];
      float pt2[genPostFSR_Pt->size()];

      for(unsigned int gn=0; gn<genPostFSR_Pt->size(); gn++)
      {
	pt2[gn]=genPostFSR_Pt->at(gn);
      }

      int sizen = sizeof(pt2)/sizeof(pt2[0]);
      TMath::Sort(sizen,pt2,index2,true);

      trigMatch_lead = false; trigMatch_slead = false;
      massGPre=0.;massReco=-999.; massGen=-999.;
      dR = 0.; dR1 = 0.; dR2 = 0.;
      isReco = false; isGen = false;
      BB = false; BE = false; EE = false;
      mediumId = 0; dzMaskId = 0; passId = 0;
      newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); newscEta.clear();
      newgenPt.clear(); newgenEta.clear(); newgenEnr.clear(); newgenPhi.clear();

      // PU Weight
      int bin = 0;
      //double puWeights = 1.0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      double puweight = weights->GetBinContent(bin);

      if(genPreFSR_Pt->size() == 2){
	genPre1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	genPre2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=genPre1+genPre2;
	massGPre=diGen.M();

      }

      if(jentry==1 && massGPre > 100.) continue;        // Gen Mass cut ----- for 50 to inf sample

      if(!tauFlag && genPreFSR_Pt->size() == 2) {
      //if(jentry==0 && tauFlag && genPreFSR_Pt->size() != 2) continue;

	if(Ele23_WPLoose) {
	  if(ptElec->size()>=2.) {

	    for(int j=0;j<ptElec->size();j++){

	      mediumId = passMediumId->at(index1[j]);

	      if(mediumId) {

		if(fabs(etaSC->at(index1[j])) < 2.5 && !(fabs(etaSC->at(index1[j])) > 1.4442 && fabs(etaSC->at(index1[j])) < 1.566)){

		  newelePt.push_back(ptElec->at(index1[j]));
		  neweleEta.push_back(etaElec->at(index1[j]));
		  newscEta.push_back(etaSC->at(index1[j]));
		  neweleEnr.push_back(energyElec->at(index1[j]));
		  newelePhi.push_back(phiElec->at(index1[j]));

		} // eta
	      } // ID
	    } // pt size
	  } // pt size >= 2
	} // trigger



	if(newelePt.size() == 2){

	  for(unsigned int k = 0; k < eta_Ele23->size(); k++){

	    double dR1_comp = 1000.;
	    double dR2_comp = 1000.;

	    dR1 = deltaR(neweleEta.at(0), newelePhi.at(0), eta_Ele23->at(k), phi_Ele23->at(k));
	    dR2 = deltaR(neweleEta.at(1), newelePhi.at(1), eta_Ele23->at(k), phi_Ele23->at(k));

	    if(dR1 < 0.1){
	      if (dR1 < dR1_comp)
	      { 
		dR1_comp = dR1;
		trigMatch_lead = true;
	      }
	    } // dR1

	    if(dR2 < 0.1){
	      if (dR2 < dR2_comp)
	      { 
		dR2_comp = dR2; 
		trigMatch_slead = true;
	      }
	    } // dR2
	  } // eta_Ele23

	  if(trigMatch_lead || trigMatch_slead){
	    if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){

	      isReco = true;
	      ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	      ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	      dielectron=ele1+ele2;
	      massReco = dielectron.M();

	      if(fabs(newscEta.at(0)) < 1.4442 && fabs(newscEta.at(1)) < 1.4442) BB = true;
	      if((fabs(newscEta.at(0)) < 1.4442 && fabs(newscEta.at(1)) > 1.566) || (fabs(newscEta.at(0)) > 1.566 && fabs(newscEta.at(1)) < 1.4442)) BE = true;
	      if(fabs(newscEta.at(0)) > 1.566 && fabs(newscEta.at(1)) > 1.566) EE = true;
	    } // pt on reco
	  }


	  if(genPostFSR_Pt->size() >= 2.){  
	    for(unsigned int k=0;k<genPostFSR_Pt->size();k++){

	      if(fabs(genPostFSR_Eta->at(index2[k])) < 2.5 && !(fabs(genPostFSR_Eta->at(index2[k])) > 1.4442 && fabs(genPostFSR_Eta->at(index2[k])) < 1.566)){

		newgenPt.push_back(genPostFSR_Pt->at(index2[k]));
		newgenEta.push_back(genPostFSR_Eta->at(index2[k]));
		newgenEnr.push_back(genPostFSR_En->at(index2[k]));
		newgenPhi.push_back(genPostFSR_Phi->at(index2[k]));
	      } // eta on gen
	    } // post FSR size
	  } // post FSR size == 2

	  if(newgenPt.size() == 2){
	    if(newgenPt.at(0) > 30. && newgenPt.at(1) > 10.){

	      isGen = true;
	      gen1.SetPtEtaPhiE(newgenPt.at(0),newgenEta.at(0),newgenPhi.at(0),newgenEnr.at(0));
	      gen2.SetPtEtaPhiE(newgenPt.at(1),newgenEta.at(1),newgenPhi.at(1),newgenEnr.at(1));

	      digen=gen1+gen2;
	      massGen=digen.M();
	    } // acceptance on gen
	  } // newgen size == 2
	} // newele size == 2

	lumiWeight = lumi_Weight;
	genWeight = theWeight;
	PUWeight = puweight;

	tree->Fill();
      } // pre FSR == 2
    } // event

    file[jentry]->Write();
    file[jentry]->Close();

  } //file Loop
}
