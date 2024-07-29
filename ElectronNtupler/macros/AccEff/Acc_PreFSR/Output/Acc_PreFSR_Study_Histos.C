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

Bool_t isPassAccCondition_GenLepton_ECALGAP(double leadpt, double leadeta, double subleadpt, double subleadeta)
{
  Bool_t isPassAcc = kFALSE;

  if( leadpt > 30. && fabs(leadeta) < 2.5 && !( fabs(leadeta) > 1.4442 && fabs(leadeta) < 1.566 ) &&
      subleadpt  > 10.  && fabs(subleadeta) < 2.5 && !( fabs(subleadeta) > 1.4442 && fabs(subleadeta) < 1.566 ) )
    isPassAcc = 1;

  return isPassAcc;
}

Bool_t isPassAccCondition_Electron(double leadpt, double leadeta, double subleadpt, double subleadeta)
{
  Bool_t isPassAcc = kFALSE;
  if( leadpt > 30 && fabs(leadeta) < 2.5 && !( fabs(leadeta) > 1.4442 && fabs(leadeta) < 1.566 ) &&
      subleadpt  > 10  && fabs(subleadeta)  < 2.5 && !( fabs(subleadeta) > 1.4442 && fabs(subleadeta) < 1.566 ) )
    isPassAcc = kTRUE;

  return isPassAcc;
}

void Acc_PreFSR_Study_Histos() {

  ofstream outfile;
  outfile.open("log_EventInfo.txt");
  
  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;

  //int mass[14] = {10,10,50,50,100,200,400,500,700,800,1000,1500,2000,3000};
  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

  //double xsec[13] = {18610./3,18610./3,5870./3,5870./3,226./3,7.67/3,0.423/3,0.24/3,0.035/3,0.03/3,0.016/3,0.002/3,0.00054/3};
  //double sumofWts[13] = {771413889185.162476,771413889185.162476,144505031098.323120,144505031098.323120,219889705.060318,7008766.904321,122987.746342,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};
  
  double xsec[11] = {18610./3,5870./3,226./3,7.67/3,0.423/3,0.24/3,0.035/3,0.03/3,0.016/3,0.002/3,0.00054/3};
  double sumofWts[11] = {771413889185.162476,144505031098.323120,219889705.060318,7008766.904321,122987.746342,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

  /*workdir = "/eos/cms/store/group/phys_smp/rchawla/DY_76X_SmearSyst/";
  InputFiles_signal_DY.clear();

  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_10to50_part1.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_10to50_part2.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_50toInf_part1.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_50toInf_part2.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_100to200.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_200to400.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_400to500.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_500to700.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_700to800.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_800to1000.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1000to1500.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1500to2000.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_2000to3000.root"));*/
  
  //workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";
  workdir = "/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_24072017/DY_Signal/";
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

  for(unsigned int j = 2; j < 3; ++j ) {

    TFile *f1 = TFile::Open("../../dataPUDist.root");
    TFile *f2 = TFile::Open("../../PileUp_MC.root");

    //data PU histo
    TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    //mc PU histo
    TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");
    TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    TTree * tmpTree = (TTree*)InputFiles_signal_DY.at(j)->Get("ntupler/ElectronTree");

    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Px;
    vector<float>   *genPreFSR_Py;
    vector<float>   *genPreFSR_Pz;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<float>   *genPostFSR_Pt;
    vector<float>   *genPostFSR_Px;
    vector<float>   *genPostFSR_Py;
    vector<float>   *genPostFSR_Pz;
    vector<float>   *genPostFSR_Eta;
    vector<float>   *genPostFSR_Phi;
    vector<float>   *genPostFSR_En;
    vector<float>   *genPhoton_Pt;
    vector<float>   *genPhoton_Px;
    vector<float>   *genPhoton_Py;
    vector<float>   *genPhoton_Pz;
    vector<float>   *genPhoton_Eta;
    vector<float>   *genPhoton_Phi;
    vector<float>   *genPhoton_En;
    Int_t           tauFlag; 
    Double_t        theWeight;
    Int_t           nPUTrue;
    Double_t        RunNo;
    Double_t        Lumi;
    Double_t        EvtNo;

    genPreFSR_Pt = 0;
    genPreFSR_Px = 0;
    genPreFSR_Py = 0;
    genPreFSR_Pz = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;
    genPostFSR_Pt = 0;
    genPostFSR_Px = 0;
    genPostFSR_Py = 0;
    genPostFSR_Pz = 0;
    genPostFSR_Eta = 0;
    genPostFSR_Phi = 0;
    genPostFSR_En = 0;
    genPhoton_Pt = 0;
    genPhoton_Px = 0;
    genPhoton_Py = 0;
    genPhoton_Pz = 0;
    genPhoton_Eta = 0;
    genPhoton_Phi = 0;
    genPhoton_En = 0;

    tmpTree->SetBranchStatus("*", 0);
    tmpTree->SetBranchStatus("genPreFSR_Pt", 1);
    tmpTree->SetBranchStatus("genPreFSR_Px", 1);
    tmpTree->SetBranchStatus("genPreFSR_Py", 1);
    tmpTree->SetBranchStatus("genPreFSR_Pz", 1);
    tmpTree->SetBranchStatus("genPreFSR_Eta", 1);
    tmpTree->SetBranchStatus("genPreFSR_Phi", 1);
    tmpTree->SetBranchStatus("genPreFSR_En", 1);
    tmpTree->SetBranchStatus("genPostFSR_Pt", 1);
    tmpTree->SetBranchStatus("genPostFSR_Px", 1);
    tmpTree->SetBranchStatus("genPostFSR_Py", 1);
    tmpTree->SetBranchStatus("genPostFSR_Pz", 1);
    tmpTree->SetBranchStatus("genPostFSR_Eta", 1);
    tmpTree->SetBranchStatus("genPostFSR_Phi", 1);
    tmpTree->SetBranchStatus("genPostFSR_En", 1);
    tmpTree->SetBranchStatus("genPhoton_Pt", 1);
    tmpTree->SetBranchStatus("genPhoton_Px", 1);
    tmpTree->SetBranchStatus("genPhoton_Py", 1);
    tmpTree->SetBranchStatus("genPhoton_Pz", 1);
    tmpTree->SetBranchStatus("genPhoton_Eta", 1);
    tmpTree->SetBranchStatus("genPhoton_Phi", 1);
    tmpTree->SetBranchStatus("genPhoton_En", 1);
    tmpTree->SetBranchStatus("tauFlag", 1);
    tmpTree->SetBranchStatus("theWeight", 1);
    tmpTree->SetBranchStatus("nPUTrue", 1);
    tmpTree->SetBranchStatus("RunNo", 1);
    tmpTree->SetBranchStatus("Lumi", 1);
    tmpTree->SetBranchStatus("EvtNo", 1);

    tmpTree->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    tmpTree->SetBranchAddress("genPreFSR_Px", &genPreFSR_Px);
    tmpTree->SetBranchAddress("genPreFSR_Py", &genPreFSR_Py);
    tmpTree->SetBranchAddress("genPreFSR_Pz", &genPreFSR_Pz);
    tmpTree->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    tmpTree->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    tmpTree->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    tmpTree->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
    tmpTree->SetBranchAddress("genPostFSR_Px", &genPostFSR_Px);
    tmpTree->SetBranchAddress("genPostFSR_Py", &genPostFSR_Py);
    tmpTree->SetBranchAddress("genPostFSR_Pz", &genPostFSR_Pz);
    tmpTree->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
    tmpTree->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
    tmpTree->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
    tmpTree->SetBranchAddress("genPhoton_Pt", &genPhoton_Pt);
    tmpTree->SetBranchAddress("genPhoton_Px", &genPhoton_Px);
    tmpTree->SetBranchAddress("genPhoton_Py", &genPhoton_Py);
    tmpTree->SetBranchAddress("genPhoton_Pz", &genPhoton_Pz);
    tmpTree->SetBranchAddress("genPhoton_Eta", &genPhoton_Eta);
    tmpTree->SetBranchAddress("genPhoton_Phi", &genPhoton_Phi);
    tmpTree->SetBranchAddress("genPhoton_En", &genPhoton_En);
    tmpTree->SetBranchAddress("tauFlag", &tauFlag);
    tmpTree->SetBranchAddress("theWeight", &theWeight);
    tmpTree->SetBranchAddress("nPUTrue", &nPUTrue);
    tmpTree->SetBranchAddress("RunNo", &RunNo);
    tmpTree->SetBranchAddress("Lumi", &Lumi);
    tmpTree->SetBranchAddress("EvtNo", &EvtNo);

    file[j] = new TFile(Form("Check/DYEE_M%dto%d.root",mass[j],mass[j+1]),"RECREATE");

    const Int_t nMassBin = 43;

    Double_t MassBinEdges[44] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1500,3000};

    TH1D *h_mass_AccTotal = new TH1D("h_mass_AccTotal", "", nMassBin, MassBinEdges);
    TH1D *h_mass_AccPass = new TH1D("h_mass_AccPass", "", nMassBin, MassBinEdges);
    TH1D *h_mass_EffTotal = new TH1D("h_mass_EffTotal", "", nMassBin, MassBinEdges);
    TH1D *h_mass_EffPass = new TH1D("h_mass_EffPass", "", nMassBin, MassBinEdges);

    h_mass_AccTotal->Sumw2(); h_mass_AccPass->Sumw2(); h_mass_EffTotal->Sumw2(); h_mass_EffPass->Sumw2();

    int count;
    double dR, dR1, massGen;
    double postFSR_Mass;
    double preFSR_Mass;
    TLorentzVector gen_preFSR,gen_preFSR1;
    TLorentzVector fourmom,fourmom1;
    TLorentzVector SumPhotonMom,SumPhotonMom1;
    TLorentzVector pre1,pre2,diPre;
    TLorentzVector post1,post2,diPost;
    TLorentzVector gen1,gen2,diGen;

    //for sorted electrons
    vector <double> newgenPt; vector <double> newgenEta; vector <double> newgenPhi; vector <double> newgenEn;
    vector <double> newphoPt; vector <double> newphoEta; vector <double> newphoPhi; vector <double> newphoEn;
    vector <double> newgenPx; vector <double> newgenPy; vector <double> newgenPz;
    vector <double> newphoPx; vector <double> newphoPy; vector <double> newphoPz;

    vector <double> gPreFSR_Pt; vector <double> gPreFSR_Eta; vector <double> gPreFSR_Phi; vector <double> gPreFSR_En;
    vector <double> gPreFSR_Px; vector <double> gPreFSR_Py; vector <double> gPreFSR_Pz; vector <double> gPreFSR_Pt1;

    count = 0.;

    TH1D *h_Pt_Lead_NoAcc_PreFSR = new TH1D("h_Pt_Lead_NoAcc_PreFSR", "h_Pt_Lead_NoAcc_PreFSR", 10000, 0,  10000);
    TH1D *h_Pt_Lead_Acc_PreFSR   = new TH1D("h_Pt_Lead_Acc_PreFSR", "h_Pt_Lead_Acc_PreFSR", 10000, 0,  10000);
    TH1D *h_Pt_SubLead_NoAcc_PreFSR = new TH1D("h_Pt_SubLead_NoAcc_PreFSR", "h_Pt_SubLead_NoAcc_PreFSR", 10000, 0,  10000);
    TH1D *h_Pt_SubLead_Acc_PreFSR   = new TH1D("h_Pt_SubLead_Acc_PreFSR", "h_Pt_SubLead_Acc_PreFSR", 10000, 0,  10000);

    TH1D *h_Eta_Lead_NoAcc_PreFSR = new TH1D("h_Eta_Lead_NoAcc_PreFSR", "h_Eta_Lead_NoAcc_PreFSR", 2000, -10, 10);
    TH1D *h_Eta_Lead_Acc_PreFSR   = new TH1D("h_Eta_Lead_Acc_PreFSR", "h_Eta_Lead_Acc_PreFSR", 600, -3, 3);
    TH1D *h_Eta_SubLead_NoAcc_PreFSR = new TH1D("h_Eta_SubLead_NoAcc_PreFSR", "h_Eta_SubLead_NoAcc_PreFSR", 2000, -10, 10);
    TH1D *h_Eta_SubLead_Acc_PreFSR   = new TH1D("h_Eta_SubLead_Acc_PreFSR", "h_Eta_SubLead_Acc_PreFSR", 600, -3, 3);

    TH1D *h_Pt_Lead_NoAcc_PostFSR = new TH1D("h_Pt_Lead_NoAcc_PostFSR", "h_Pt_Lead_NoAcc_PostFSR", 10000, 0,  10000);
    TH1D *h_Pt_SubLead_NoAcc_PostFSR = new TH1D("h_Pt_SubLead_NoAcc_PostFSR", "h_Pt_SubLead_NoAcc_PostFSR", 10000, 0,  10000);
    TH1D *h_Eta_Lead_NoAcc_PostFSR = new TH1D("h_Eta_Lead_NoAcc_PostFSR", "h_Eta_Lead_NoAcc_PostFSR", 2000, -10, 10);
    TH1D *h_Eta_SubLead_NoAcc_PostFSR = new TH1D("h_Eta_SubLead_NoAcc_PostFSR", "h_Eta_SubLead_NoAcc_PostFSR", 2000, -10, 10);

    h_Pt_Lead_NoAcc_PreFSR->Sumw2(); h_Pt_Lead_Acc_PreFSR->Sumw2(); h_Pt_SubLead_NoAcc_PreFSR->Sumw2(); h_Pt_SubLead_Acc_PreFSR->Sumw2();
    h_Eta_Lead_NoAcc_PreFSR->Sumw2(); h_Eta_Lead_Acc_PreFSR->Sumw2(); h_Eta_SubLead_NoAcc_PreFSR->Sumw2(); h_Eta_SubLead_Acc_PreFSR->Sumw2();
    
    h_Pt_Lead_NoAcc_PostFSR->Sumw2(); h_Pt_SubLead_NoAcc_PostFSR->Sumw2();
    h_Eta_Lead_NoAcc_PostFSR->Sumw2(); h_Eta_SubLead_NoAcc_PostFSR->Sumw2();

    double lumiWeight = (xsec[j]/sumofWts[j])*2258.066;
    cout<<"DY Sample: "<<mass[j]<<"to"<<mass[j+1]<<endl;
    cout<<"Xsec: "<<xsec[j]<<"   Sum of weights: "<<sumofWts[j]<<endl;
    //cout << "Lumiweight = " << lumiWeight << endl;

    int nentries = tmpTree->GetEntries();
    //int nentries = 1000;
    cout<<"entries: "<<nentries<<endl;

    for (unsigned int k=0; k < nentries; k++){

      tmpTree->GetEntry(k);

      if(k%1000000 == 0){
	cout << "Events Processed :  " << k << endl;
      }

      //PUWeight
      int bin = 0;
      double puWeights = 1.0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      puWeights = weights->GetBinContent(bin);

      // Sorting Gen Photons
      int index1[genPhoton_Pt->size()];
      float pt1[genPhoton_Pt->size()];

      for(unsigned int ph=0; ph<genPhoton_Pt->size(); ph++)
      {
	pt1[ph]=genPhoton_Pt->at(ph);
      }

      int size1 = sizeof(pt1)/sizeof(pt1[0]);
      TMath::Sort(size1,pt1,index1,true);

      // sorting of gen electrons : postFSR
      int index2[genPostFSR_Pt->size()];
      float pt2[genPostFSR_Pt->size()];

      for(unsigned int b=0; b<genPostFSR_Pt->size(); b++)
      {
	pt2[b]=genPostFSR_Pt->at(b);
      }
      int size2 = sizeof(pt2)/sizeof(pt2[0]);
      TMath::Sort(size2,pt2,index2,true);

      //clearing of vectors
      dR = 0.; dR1 = 0.; massGen = 0.0;
      postFSR_Mass = -999.; preFSR_Mass = -999.;
      newgenPt.clear(); newgenPt.clear(); newgenEta.clear(); newgenPhi.clear(); newgenEn.clear();
      newphoPt.clear(); newgenPt.clear(); newphoEta.clear(); newphoPhi.clear(); newphoEn.clear();
      newgenPx.clear(); newgenPy.clear(); newgenPz.clear();
      newphoPx.clear(); newphoPy.clear(); newphoPz.clear();
      gPreFSR_Pt.clear(); gPreFSR_Eta.clear(); gPreFSR_Phi.clear(); gPreFSR_En.clear();
      gPreFSR_Px.clear(); gPreFSR_Py.clear(); gPreFSR_Pz.clear(); gPreFSR_Pt1.clear();

      if(genPreFSR_Pt->size() == 2){

	double E1 = sqrt(genPreFSR_Px->at(0)*genPreFSR_Px->at(0) + genPreFSR_Py->at(0)*genPreFSR_Py->at(0) + genPreFSR_Pz->at(0)*genPreFSR_Pz->at(0) + 0);
	double E2 = sqrt(genPreFSR_Px->at(1)*genPreFSR_Px->at(1) + genPreFSR_Py->at(1)*genPreFSR_Py->at(1) + genPreFSR_Pz->at(1)*genPreFSR_Pz->at(1) + 0);
	gen1.SetPxPyPzE(genPreFSR_Px->at(0),genPreFSR_Py->at(0),genPreFSR_Pz->at(0),E1);
	gen2.SetPxPyPzE(genPreFSR_Px->at(1),genPreFSR_Py->at(1),genPreFSR_Pz->at(1),E2);

	//gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	//gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=gen1+gen2;
	massGen=diGen.M();

	//cout<<"Pre-FSR mass with SetPxPyPzE = "<<massGen<<endl;
	//cout<<"pre FSR Mass with SetPtEtaPhiE = "<<massGen<<endl;
      }

      //if((j==2 || j==3) && massGen >= 100) continue;
      if(j==1 && massGen >= 100) continue;

      if(!tauFlag && genPreFSR_Pt->size() == 2){

	for(unsigned int j=0;j<genPhoton_Pt->size();j++){

	  newphoPt.push_back(genPhoton_Pt->at(index1[j]));
	  newphoPx.push_back(genPhoton_Px->at(index1[j]));
	  newphoPy.push_back(genPhoton_Py->at(index1[j]));
	  newphoPz.push_back(genPhoton_Pz->at(index1[j]));
	  newphoEta.push_back(genPhoton_Eta->at(index1[j]));
	  newphoPhi.push_back(genPhoton_Phi->at(index1[j]));
	  newphoEn.push_back(genPhoton_En->at(index1[j]));
	}

	for(unsigned int i=0;i<genPostFSR_Pt->size();i++){

	  newgenPt.push_back(genPostFSR_Pt->at(index2[i]));
	  newgenPx.push_back(genPostFSR_Px->at(index2[i]));
	  newgenPy.push_back(genPostFSR_Py->at(index2[i]));
	  newgenPz.push_back(genPostFSR_Pz->at(index2[i]));
	  newgenEta.push_back(genPostFSR_Eta->at(index2[i]));
	  newgenPhi.push_back(genPostFSR_Phi->at(index2[i]));
	  newgenEn.push_back(genPostFSR_En->at(index2[i]));
	}

	double E3 = sqrt(newgenPx.at(0)*newgenPx.at(0) + newgenPy.at(0)*newgenPy.at(0) + newgenPz.at(0)*newgenPz.at(0) + 0);
	double E4 = sqrt(newgenPx.at(1)*newgenPx.at(1) + newgenPy.at(1)*newgenPy.at(1) + newgenPz.at(1)*newgenPz.at(1) + 0);

	post1.SetPxPyPzE(newgenPx.at(0),newgenPy.at(0),newgenPz.at(0),E3);
	post2.SetPxPyPzE(newgenPx.at(1),newgenPy.at(1),newgenPz.at(1),E4);

	//post1.SetPtEtaPhiE(newgenPt.at(0),newgenEta.at(0),newgenPhi.at(0),newgenEn.at(0));
	//post2.SetPtEtaPhiE(newgenPt.at(1),newgenEta.at(1),newgenPhi.at(1),newgenEn.at(1));

	diPost=post1+post2;
	postFSR_Mass=diPost.M();

	//cout<<"Post-FSR mass = "<<postFSR_Mass<<endl;

	for(unsigned int igen1 = 0; igen1 < newgenPt.size(); igen1++){
	  SumPhotonMom1.SetPxPyPzE(0.,0.,0.,0.);

	  if(newphoPt.size() >= 1.){

	    for(unsigned int ipho1 = 0; ipho1 < newphoPt.size(); ipho1++){

	      double E51 = sqrt(newgenPx.at(igen1)*newgenPx.at(igen1) + newgenPy.at(igen1)*newgenPy.at(igen1) + newgenPz.at(igen1)*newgenPz.at(igen1) + 0);
	      gen_preFSR1.SetPxPyPzE(newgenPx.at(igen1), newgenPy.at(igen1), newgenPz.at(igen1), E51);
	      fourmom1.SetPxPyPzE(0.,0.,0.,0.);

	      dR1 = deltaR(newphoEta.at(ipho1), newphoPhi.at(ipho1), newgenEta.at(igen1), newgenPhi.at(igen1));

	      if(dR1 < 0.1){

		double E61 = sqrt(newphoPx.at(ipho1)*newphoPx.at(ipho1) + newphoPy.at(ipho1)*newphoPy.at(ipho1) + newphoPz.at(ipho1)*newphoPz.at(ipho1) + 0);
		fourmom1.SetPxPyPzE(newphoPx.at(ipho1), newphoPy.at(ipho1), newphoPz.at(ipho1), E61);
		SumPhotonMom1 = SumPhotonMom1 + fourmom1;

	      }
	    }

	    gen_preFSR1 = gen_preFSR1 + SumPhotonMom1;
	    gPreFSR_Pt1.push_back(gen_preFSR1.Pt());

	  }

	  else {
	    gPreFSR_Pt1.push_back(newgenPt.at(igen1));
	  }
	}

	// sorting of gen electrons : dressed-level
	int index3[gPreFSR_Pt1.size()];
	float pt3[gPreFSR_Pt1.size()];

	for(unsigned int c=0; c<gPreFSR_Pt1.size(); c++)
	{
	  pt3[c]=gPreFSR_Pt1.at(c);
	}
	int size3 = sizeof(pt3)/sizeof(pt3[0]);
	TMath::Sort(size3,pt3,index3,true);

	//if(EvtNo == 43936){

	  //cout<<"Size = "<<newgenPt.size()<<"   "<<gPreFSR_Pt1.size()<<endl;
	  //cout<<"Post-FSR     Pt Lead = "<<newgenPt.at(0)<<"   Sub-lead = "<<newgenPt.at(1)<<endl;
	  //cout<<"Dressed      Pt Lead = "<<gPreFSR_Pt1.at(0)<<"   Sub-lead = "<<gPreFSR_Pt1.at(1)<<endl;
	  //cout<<""<<endl;

	  if(count >= 1000.) continue;

	  if(gPreFSR_Pt1.size() == 2){

	    if(gPreFSR_Pt1.at(index3[1]) > 40. && gPreFSR_Pt1.at(index3[1]) < 60.) {

	      count++;

	      outfile<<"========================================================================================="<<endl;
	      outfile<<setprecision(0)<<fixed<<"[run, lumi, event] = ("<<RunNo<<", "<<Lumi<<", "<<EvtNo<<")"<<endl;
	      outfile<<""<<endl;

	      outfile<<setprecision(3)<<fixed<<"[genlep_postFSR1] (Pt, eta, phi) = (   "<<newgenPt.at(0)<<",    "<<newgenEta.at(0)<<",    "<<newgenPhi.at(0)<<")"<<endl;
	      outfile<<setprecision(3)<<fixed<<"[genlep_postFSR2] (Pt, eta, phi) = (   "<<newgenPt.at(1)<<",    "<<newgenEta.at(1)<<",    "<<newgenPhi.at(1)<<")"<<endl;
	      outfile<<""<<endl;

	      for(unsigned int igen = 0; igen < newgenPt.size(); igen++){
		SumPhotonMom.SetPxPyPzE(0.,0.,0.,0.);
		//SumPhotonMom.SetPtEtaPhiE(0.,0.,0.,0.);

		outfile<<"##### Dressing post-FSR"<<igen+1<<" #####"<<endl;
		outfile<<"[post-FSR] (Pt, eta, phi) = (   "<<newgenPt.at(igen)<<",    "<<newgenEta.at(igen)<<",    "<<newgenPhi.at(igen)<<")"<<endl;

		if(newphoPt.size() >= 1.){

		  for(unsigned int ipho = 0; ipho < newphoPt.size(); ipho++){

		    double E5 = sqrt(newgenPx.at(igen)*newgenPx.at(igen) + newgenPy.at(igen)*newgenPy.at(igen) + newgenPz.at(igen)*newgenPz.at(igen) + 0);
		    gen_preFSR.SetPxPyPzE(newgenPx.at(igen), newgenPy.at(igen), newgenPz.at(igen), E5);
		    //gen_preFSR.SetPtEtaPhiE(newgenPt.at(igen), newgenEta.at(igen), newgenPhi.at(igen), newgenEn.at(igen));

		    fourmom.SetPxPyPzE(0.,0.,0.,0.);
		    //fourmom.SetPtEtaPhiE(0.,0.,0.,0.);

		    dR = deltaR(newphoEta.at(ipho), newphoPhi.at(ipho), newgenEta.at(igen), newgenPhi.at(igen));

		    if(dR < 0.1){

		      //cout<<"dR = "<<dR<<endl;

		      double E6 = sqrt(newphoPx.at(ipho)*newphoPx.at(ipho) + newphoPy.at(ipho)*newphoPy.at(ipho) + newphoPz.at(ipho)*newphoPz.at(ipho) + 0);
		      fourmom.SetPxPyPzE(newphoPx.at(ipho), newphoPy.at(ipho), newphoPz.at(ipho), E6);
		      //fourmom.SetPtEtaPhiE(newphoPt.at(ipho), newphoEta.at(ipho), newphoPhi.at(ipho), newphoEn.at(ipho));

		      SumPhotonMom = SumPhotonMom + fourmom;

		      outfile<<"        [Photon (dR="<<dR<<")] (Pt, eta, phi) = (    "<<fourmom.Pt()<<",    "<<fourmom.Eta()<<",    "<<fourmom.Phi()<<") -> Sum of photon momentum: (Pt, eta, phi) = (    "<<fourmom.Pt()<<",    "<<fourmom.Eta()<<",    "<<fourmom.Phi()<<endl;
		      //outfile<<"        [Photon (dR="<<dR<<")] (Pt, eta, phi) = (    "<<SumPhotonMom.Pt()<<",    "<<SumPhotonMom.Eta()<<",    "<<SumPhotonMom.Phi()<<") -> Sum of photon momentum: (Pt, eta, phi) = (    "<<SumPhotonMom.Pt()<<",    "<<SumPhotonMom.Eta()<<",    "<<SumPhotonMom.Phi()<<endl;

		    }

		    //outfile<<"        [Photon (dR="<<dR<<")] (Pt, eta, phi) = (    "<<SumPhotonMom.Pt()<<",    "<<SumPhotonMom.Eta()<<",    "<<SumPhotonMom.Phi()<<") -> Sum of photon momentum: (Pt, eta, phi) = (    "<<SumPhotonMom.Pt()<<",    "<<SumPhotonMom.Eta()<<",    "<<SumPhotonMom.Phi()<<endl;
		  }

		  gen_preFSR = gen_preFSR + SumPhotonMom;

		  outfile<<"[dressed lepton] (Pt, eta, phi) = (   "<<gen_preFSR.Pt()<<",    "<<gen_preFSR.Eta()<<",    "<<gen_preFSR.Phi()<<")"<<endl;
		  outfile<<""<<endl;

		  gPreFSR_Pt.push_back(gen_preFSR.Pt());
		  gPreFSR_Px.push_back(gen_preFSR.Px());
		  gPreFSR_Py.push_back(gen_preFSR.Py());
		  gPreFSR_Pz.push_back(gen_preFSR.Pz());
		  gPreFSR_Eta.push_back(gen_preFSR.Eta());
		  gPreFSR_Phi.push_back(gen_preFSR.Phi());
		  gPreFSR_En.push_back(gen_preFSR.Energy());

		  //cout<<"Px = "<<gen_preFSR.Px()<<endl;
		}

		else {

		  outfile<<"[dressed lepton] (Pt, eta, phi) = (   "<<gen_preFSR.Pt()<<",    "<<gen_preFSR.Eta()<<",    "<<gen_preFSR.Phi()<<")"<<endl;
		  outfile<<""<<endl;

		  gPreFSR_Pt.push_back(newgenPt.at(igen));
		  gPreFSR_Px.push_back(newgenPx.at(igen));
		  gPreFSR_Py.push_back(newgenPy.at(igen));
		  gPreFSR_Pz.push_back(newgenPz.at(igen));
		  gPreFSR_Eta.push_back(newgenEta.at(igen));
		  gPreFSR_Phi.push_back(newgenPhi.at(igen));
		  gPreFSR_En.push_back(newgenEn.at(igen));
		}
	      }

	      // sorting of gen electrons : dressed-level
	      int index4[gPreFSR_Pt.size()];
	      float pt4[gPreFSR_Pt.size()];

	      for(unsigned int d=0; d<gPreFSR_Pt.size(); d++)
	      {                       
		pt4[d]=gPreFSR_Pt.at(d);
	      }                                     
	      int size4 = sizeof(pt4)/sizeof(pt4[0]);     
	      TMath::Sort(size4,pt4,index4,true);

	      outfile<<""<<endl;
	      outfile<<"[dressed lepton1] (Pt, eta, phi) = (   "<<gPreFSR_Pt.at(index4[0])<<",    "<<gPreFSR_Eta.at(index4[0])<<",    "<<gPreFSR_Phi.at(index4[0])<<")"<<endl;
	      outfile<<"[dressed lepton2] (Pt, eta, phi) = (   "<<gPreFSR_Pt.at(index4[1])<<",    "<<gPreFSR_Eta.at(index4[1])<<",    "<<gPreFSR_Phi.at(index4[1])<<")"<<endl;
	      outfile<<"========================================================================================="<<endl;
	      outfile<<""<<endl;

	      double E7 = sqrt(gPreFSR_Px.at(0)*gPreFSR_Px.at(0) + gPreFSR_Py.at(0)*gPreFSR_Py.at(0) + gPreFSR_Pz.at(0)*gPreFSR_Pz.at(0) + 0);
	      double E8 = sqrt(gPreFSR_Px.at(1)*gPreFSR_Px.at(1) + gPreFSR_Py.at(1)*gPreFSR_Py.at(1) + gPreFSR_Pz.at(1)*gPreFSR_Pz.at(1) + 0);

	      //cout<<"E7 = "<<E7<<"   E8 = "<<E8<<endl;

	      pre1.SetPxPyPzE(gPreFSR_Px.at(0),gPreFSR_Py.at(0),gPreFSR_Pz.at(0),E7);
	      pre2.SetPxPyPzE(gPreFSR_Px.at(1),gPreFSR_Py.at(1),gPreFSR_Pz.at(1),E8);

	      //pre1.SetPtEtaPhiE(gPreFSR_Pt.at(0),gPreFSR_Eta.at(0),gPreFSR_Phi.at(0),gPreFSR_En.at(0));
	      //pre2.SetPtEtaPhiE(gPreFSR_Pt.at(1),gPreFSR_Eta.at(1),gPreFSR_Phi.at(1),gPreFSR_En.at(1));

	      diPre=pre1+pre2;
	      preFSR_Mass=diPre.M();

	      //cout<<"Mass = "<<preFSR_Mass<<endl;

	      Bool_t Flag_PassAcc = kFALSE;

	      Flag_PassAcc = isPassAccCondition_GenLepton_ECALGAP(gPreFSR_Pt.at(0),gPreFSR_Eta.at(0),gPreFSR_Pt.at(1),gPreFSR_Eta.at(1));

	      h_Pt_Lead_NoAcc_PostFSR->Fill(newgenPt.at(0),theWeight*lumiWeight);
	      h_Pt_SubLead_NoAcc_PostFSR->Fill(newgenPt.at(1),theWeight*lumiWeight);
	      h_Eta_Lead_NoAcc_PostFSR->Fill(newgenEta.at(0),theWeight*lumiWeight);
	      h_Eta_SubLead_NoAcc_PostFSR->Fill(newgenEta.at(1),theWeight*lumiWeight);

	      h_Pt_Lead_NoAcc_PreFSR->Fill(gPreFSR_Pt.at(0),theWeight*lumiWeight);
	      h_Pt_SubLead_NoAcc_PreFSR->Fill(gPreFSR_Pt.at(1),theWeight*lumiWeight);
	      h_Eta_Lead_NoAcc_PreFSR->Fill(gPreFSR_Eta.at(0),theWeight*lumiWeight);
	      h_Eta_SubLead_NoAcc_PreFSR->Fill(gPreFSR_Eta.at(1),theWeight*lumiWeight);

	      if(Flag_PassAcc) {
		h_Pt_Lead_Acc_PreFSR->Fill(gPreFSR_Pt.at(0),theWeight*lumiWeight);
		h_Pt_SubLead_Acc_PreFSR->Fill(gPreFSR_Pt.at(1),theWeight*lumiWeight);
		h_Eta_Lead_Acc_PreFSR->Fill(gPreFSR_Eta.at(0),theWeight*lumiWeight);
		h_Eta_SubLead_Acc_PreFSR->Fill(gPreFSR_Eta.at(1),theWeight*lumiWeight);

		h_mass_AccTotal->Fill(preFSR_Mass,theWeight*lumiWeight);
		h_mass_AccPass->Fill(preFSR_Mass,theWeight*lumiWeight);
	      }

	      else
	      {
		h_mass_AccTotal->Fill(preFSR_Mass,theWeight*lumiWeight);
	      }

	    }
	  }
	//}
      }  // tauFlag
    } //event loop

    file[j]->Write();
    cout<<"Events = "<<count<<endl;
    cout<<""<<endl;

  }//file loop
}
