#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
#include <TFile.h>
#include <TTree.h>

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

    vector<float>   *genPostFSR_Pt;
    vector<float>   *genPostFSR_Eta;
    vector<float>   *genPostFSR_Phi;
    vector<float>   *genPostFSR_En;
    vector<float>   *genPreFSR_Pt; 
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    Int_t           tauFlag;
    Double_t        theWeight;

    genPostFSR_Pt = 0;
    genPostFSR_Eta = 0;
    genPostFSR_Phi = 0;
    genPostFSR_En = 0;
    genPreFSR_Pt = 0; 
    genPreFSR_Eta = 0; 
    genPreFSR_Phi = 0; 
    genPreFSR_En = 0;

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
    T1->SetBranchStatus("theWeight", 1);

    T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
    T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
    T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
    T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt); 
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta); 
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi); 
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);

    file[jentry] = new TFile(Form("DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    double massPreGen,massPostGen;
    TLorentzVector pregen1,pregen2,diPreGen;
    TLorentzVector postgen1,postgen2,diPostGen;

    // Branch variable declaration
    double lumiWeight, genWeight;
    double preFSR_ZMass;
    double postFSR_ZMass;

    // Branch declaration
    tree->Branch("preFSR_ZMass", &preFSR_ZMass, "preFSR_ZMass/D");
    tree->Branch("postFSR_ZMass", &postFSR_ZMass, "postFSR_ZMass/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    int nentries = T1->GetEntries();
    //int nentries = 50000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      if(tauFlag) continue;

      if(genPreFSR_Pt->size() == 2){
	pregen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	pregen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diPreGen=pregen1+pregen2;
	
	if(jentry==1 && diPreGen.M() >= 100.) continue;
	preFSR_ZMass=diPreGen.M();

	postgen1.SetPtEtaPhiE(genPostFSR_Pt->at(0),genPostFSR_Eta->at(0),genPostFSR_Phi->at(0),genPostFSR_En->at(0));
	postgen2.SetPtEtaPhiE(genPostFSR_Pt->at(1),genPostFSR_Eta->at(1),genPostFSR_Phi->at(1),genPostFSR_En->at(1));

	diPostGen=postgen1+postgen2;
	postFSR_ZMass=diPostGen.M();

	lumiWeight = lumi_Weight;
	genWeight  = theWeight;

	tree->Fill();
      }

    } // event Loop

    file[jentry]->Write();
    file[jentry]->Close();

    cout<<""<<endl;
  } // file Loop
}
