#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void treeRead() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

  double xsec[11] = {18609.9/3,5789./3,226./3,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
  double sumofWts[11] = {758771667841.549683,150453917158.292908,219884886.214768,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

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

    vector<float>   *genPostFSR_Pt;
    vector<float>   *genPostFSR_Eta;
    vector<float>   *genPostFSR_Rap;
    vector<float>   *genPostFSR_Phi;
    vector<float>   *genPostFSR_En;
    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Rap;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<float>   *genPhoton_Pt;
    vector<float>   *genPhoton_Eta;
    vector<float>   *genPhoton_Rap;
    vector<float>   *genPhoton_Phi;
    vector<float>   *genPhoton_En;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *rapElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *chargeElec;
    vector<float>   *etaSC;
    vector<int>     *passMediumId;
    vector<double>  *pt_Ele23;
    vector<double>  *eta_Ele23;
    vector<double>  *phi_Ele23;
    Bool_t          Ele23_WPLoose;
    Int_t           tauFlag;
    Int_t           nPV;
    Int_t           nPUTrue;
    Double_t        theWeight;
    vector<float>   *ptMuon;
    vector<float>   *etaMuon;
    vector<float>   *phiMuon;
    vector<float>   *energyMuon;
    vector<float>   *chargeMuon;
    vector<float>   *isoPFMuon;
    vector<bool>    *isTightMuon;
    bool            Mu8_Ele17;

    genPostFSR_Pt = 0;
    genPostFSR_Eta = 0;
    genPostFSR_Rap = 0;
    genPostFSR_Phi = 0;
    genPostFSR_En = 0;
    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Rap = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;
    genPhoton_Pt = 0;
    genPhoton_Eta = 0;
    genPhoton_Rap = 0;
    genPhoton_Phi = 0;
    genPhoton_En = 0;
    ptElec = 0;
    etaElec = 0;
    rapElec = 0;
    phiElec = 0;
    energyElec = 0;
    chargeElec = 0;
    etaSC = 0;
    passMediumId = 0;
    pt_Ele23 = 0;
    eta_Ele23 = 0;
    phi_Ele23 = 0;
    ptMuon = 0;
    etaMuon = 0;
    phiMuon = 0;
    energyMuon = 0;
    chargeMuon = 0;
    isoPFMuon = 0;
    isTightMuon = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("genPostFSR_Pt", 1);
    T1->SetBranchStatus("genPostFSR_Eta", 1);
    T1->SetBranchStatus("genPostFSR_Rap", 1);
    T1->SetBranchStatus("genPostFSR_Phi", 1);
    T1->SetBranchStatus("genPostFSR_En", 1);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Rap", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("genPhoton_Pt", 1);
    T1->SetBranchStatus("genPhoton_Eta", 1);
    T1->SetBranchStatus("genPhoton_Rap", 1);
    T1->SetBranchStatus("genPhoton_Phi", 1);
    T1->SetBranchStatus("genPhoton_En", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("rapElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("chargeElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("pt_Ele23", 1);
    T1->SetBranchStatus("eta_Ele23", 1);
    T1->SetBranchStatus("phi_Ele23", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("nPV", 1);
    T1->SetBranchStatus("nPUTrue", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("ptMuon", 1);
    T1->SetBranchStatus("etaMuon", 1);
    T1->SetBranchStatus("phiMuon", 1);
    T1->SetBranchStatus("energyMuon", 1);
    T1->SetBranchStatus("chargeMuon", 1);
    T1->SetBranchStatus("isoPFMuon", 1);
    T1->SetBranchStatus("isTightMuon", 1);
    T1->SetBranchStatus("Mu8_Ele17", 1);


    T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
    T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
    T1->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap);
    T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
    T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("genPhoton_Pt", &genPhoton_Pt);
    T1->SetBranchAddress("genPhoton_Eta", &genPhoton_Eta);
    T1->SetBranchAddress("genPhoton_Rap", &genPhoton_Rap);
    T1->SetBranchAddress("genPhoton_Phi", &genPhoton_Phi);
    T1->SetBranchAddress("genPhoton_En", &genPhoton_En);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("rapElec", &rapElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("chargeElec", &chargeElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("pt_Ele23", &pt_Ele23);
    T1->SetBranchAddress("eta_Ele23", &eta_Ele23);
    T1->SetBranchAddress("phi_Ele23", &phi_Ele23);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("nPV", &nPV);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("ptMuon", &ptMuon);
    T1->SetBranchAddress("etaMuon", &etaMuon);
    T1->SetBranchAddress("phiMuon", &phiMuon);
    T1->SetBranchAddress("energyMuon", &energyMuon);
    T1->SetBranchAddress("chargeMuon", &chargeMuon);
    T1->SetBranchAddress("isoPFMuon", &isoPFMuon);
    T1->SetBranchAddress("isTightMuon", &isTightMuon);
    T1->SetBranchAddress("Mu8_Ele17", &Mu8_Ele17);

    file[jentry] = new TFile(Form("Calibrated/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    // Branch declaration
    tree->Branch("genPostFSR_Pt", &genPostFSR_Pt);
    tree->Branch("genPostFSR_Eta", &genPostFSR_Eta);
    tree->Branch("genPostFSR_Rap", &genPostFSR_Rap);
    tree->Branch("genPostFSR_Phi", &genPostFSR_Phi);
    tree->Branch("genPostFSR_En", &genPostFSR_En);
    tree->Branch("genPreFSR_Pt", &genPreFSR_Pt);
    tree->Branch("genPreFSR_Eta", &genPreFSR_Eta);
    tree->Branch("genPreFSR_Rap", &genPreFSR_Rap);
    tree->Branch("genPreFSR_Phi", &genPreFSR_Phi);
    tree->Branch("genPreFSR_En", &genPreFSR_En);
    tree->Branch("genPhoton_Pt", &genPhoton_Pt);
    tree->Branch("genPhoton_Eta", &genPhoton_Eta);
    tree->Branch("genPhoton_Rap", &genPhoton_Rap);
    tree->Branch("genPhoton_Phi", &genPhoton_Phi);
    tree->Branch("genPhoton_En", &genPhoton_En);
    tree->Branch("ptElec", &ptElec);
    tree->Branch("etaElec", &etaElec);
    tree->Branch("rapElec", &rapElec);
    tree->Branch("phiElec", &phiElec);
    tree->Branch("energyElec", &energyElec);
    tree->Branch("chargeElec", &chargeElec);
    tree->Branch("etaSC", &etaSC);
    tree->Branch("passMediumId", &passMediumId);
    tree->Branch("pt_Ele23", &pt_Ele23);
    tree->Branch("eta_Ele23", &eta_Ele23);
    tree->Branch("phi_Ele23", &phi_Ele23);
    tree->Branch("tauFlag", &tauFlag, "tauFlag/I");
    tree->Branch("nPV", &nPV, "nPV/I");
    tree->Branch("nPUTrue", &nPUTrue, "nPUTrue/I");
    tree->Branch("theWeight", &theWeight, "theWeight/D");
    tree->Branch("Ele23_WPLoose", &Ele23_WPLoose, "Ele23_WPLoose/B");
    tree->Branch("ptMuon", &ptMuon);
    tree->Branch("etaMuon", &etaMuon);
    tree->Branch("phiMuon", &phiMuon);
    tree->Branch("energyMuon", &energyMuon);
    tree->Branch("chargeMuon", &chargeMuon);
    tree->Branch("isoPFMuon", &isoPFMuon);
    tree->Branch("isTightMuon", &isTightMuon);
    tree->Branch("Mu8_Ele17", &Mu8_Ele17, "Mu8_Ele17/B");

    int nentries = T1->GetEntries();
    //int nentries = 50000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      /*for(int k=0; k<genPostFSR_Pt->size(); k++) {
	cout<<"entry: "<<i<<"   "<<"gen Post pt: "<<genPostFSR_Pt->at(k)<<endl;
	}*/

      tree->Fill();
    } // event

    file[jentry]->Write();
    file[jentry]->Close();

  } // file Loop
}
