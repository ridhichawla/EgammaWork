#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void subtractEWK_2L_tautau() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  
  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

  double xsec[11] = {6016.88,1873.52,76.2401,2.67606,0.139728,0.0792496,0.0123176,0.01042,0.00552772,0.000741613,0.000178737};
  //double xsec[11] = {18610./3,5870./3,226./3,7.67/3,0.423/3,0.24/3,0.035/3,0.03/3,0.016/3,0.002/3,0.00054/3};
  double sumofWts[11] = {771413889185.162476,144505031098.323120,219889705.060318,7008766.904321,122987.746342,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701}; 

  workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";
  
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

    TFile *f1 = TFile::Open("../../../../../dataPUDist.root");
    TFile *f2 = TFile::Open("../../../../../PileUp_MC.root");

    // data histogram 
    TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1D *MC_puDist = (TH1D*)f2->Get("pileup_MC");
    TH1D *weights = (TH1D*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    vector<int>     *passMediumId;
    Int_t           tauFlag;
    Double_t        theWeight;
    Bool_t          Ele23_WPLoose;
    Int_t           nPUTrue;

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

    T1->SetBranchStatus("*", 0);
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

    file[jentry] = new TFile(Form("2Loose/DYTT_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    //file[jentry] = new TFile(Form("Test/DYTT_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    int count;
    int mediumId;
    bool passKin;
    double massGen;
    TLorentzVector ele1,ele2,diReco;
    TLorentzVector gen1,gen2,diGen;

    // Branch variable declaration
    double Ele1PT; double Ele2PT;
    double Ele1Eta; double Ele2Eta;
    double qcdEstMass;
    double lumiWeight, genWeight, PUWeight;

    // Branch declaration
    tree->Branch("Ele1PT", &Ele1PT, "Ele1PT/D");
    tree->Branch("Ele2PT", &Ele2PT, "Ele2PT/D");
    tree->Branch("Ele1Eta", &Ele1Eta, "Ele1Eta/D");
    tree->Branch("Ele2Eta", &Ele2Eta, "Ele2Eta/D");
    tree->Branch("qcdEstMass", &qcdEstMass, "qcdEstMass/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/D");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    vector<int> idx;

    int nentries = T1->GetEntries();
    //int nentries = 50000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      int index[ptElec->size()];
      float pt[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++) {
	pt[el]=ptElec->at(el); }

      int size = sizeof(pt)/sizeof(pt[0]);
      TMath::Sort(size,pt,index,true);

      count = 0;
      mediumId = 0;
      bool etacut = false;
      passKin = false;
      massGen = 0.0;
      idx.clear();

      // PU Weight
      int bin = 0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      double puweight = weights->GetBinContent(bin);

      if(genPreFSR_Pt->size() > 2) cout<<"size > 2"<<endl;

      if(genPreFSR_Pt->size() == 2){
	gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=gen1+gen2;
	massGen=diGen.M();

      }

      if(!Ele23_WPLoose) continue;
      //if(Ele23_WPLoose) continue;
      
      if(jentry==1 && massGen >= 100.) continue;         // Gen Mass cut ----- for 50 to inf sample
      if(!tauFlag) continue;                             // taus

      for(int j=0;j<ptElec->size();j++){

	mediumId = passMediumId->at(index[j]);
	//passKin = (fabs(etaSC->at(index[j])) < 2.5 && !(fabs(etaSC->at(index[j])) > 1.4442 && fabs(etaSC->at(index[j])) < 1.566));

	if(!mediumId){
	  count++;
	  idx.push_back(index[j]);
	}

      } //pt size > 2.

      if(count == 2){

	if(idx.size() != 2) cout<<"idx size: "<<idx.size()<<endl;

	if((fabs(etaSC->at(idx[0])) < 2.5 && !(fabs(etaSC->at(idx[0])) > 1.4442 && fabs(etaSC->at(idx[0])) < 1.566)) && (fabs(etaSC->at(idx[1])) < 2.5 && !(fabs(etaSC->at(idx[1])) > 1.4442 && fabs(etaSC->at(idx[1])) < 1.566))) etacut = true;

	if(ptElec->at(idx[0]) > 30 && ptElec->at(idx[1]) > 10 && etacut){
	  Ele1PT  = ptElec->at(idx[0]);
	  Ele1Eta = etaElec->at(idx[0]);

	  Ele2PT  = ptElec->at(idx[1]);
	  Ele2Eta = etaElec->at(idx[1]);

	  ele1.SetPtEtaPhiE(ptElec->at(idx[0]),etaElec->at(idx[0]),phiElec->at(idx[0]),energyElec->at(idx[0]));
	  ele2.SetPtEtaPhiE(ptElec->at(idx[1]),etaElec->at(idx[1]),phiElec->at(idx[1]),energyElec->at(idx[1]));

	  diReco=ele1+ele2;
	  qcdEstMass = diReco.M();

	  lumiWeight = lumi_Weight;
	  genWeight  = theWeight;
	  PUWeight = puweight;

	  tree->Fill();

	} //pt
      }// count==1

    } // event

    file[jentry]->Write();
    file[jentry]->Close();
    
    cout<<""<<endl;

  } // file Loop
}
