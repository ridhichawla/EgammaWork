#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TMath.h>
#include <TText.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TColor.h>
#include <TLatex.h>
#include <TEfficiency.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
using namespace std;

bool acceptance(double pt, double eta) {
  if( pt > 10. && fabs(eta) < 2.5 && (fabs(eta) < 1.4442 || fabs(eta) > 1.566) ) return true;
  else return false;
}

double FR_template(double pt, double eta){

  double fakerate = -999;

  if(fabs(eta) < 1.4442){

    double FR[] = {0.127899,0.141108,0.130491,0.131567,0.15511,0.149678,0.13468};
    if(pt < 25) fakerate = FR[0];
    else if(pt >= 25 && pt < 35) fakerate = FR[1];
    else if(pt >= 35 && pt < 45) fakerate = FR[2];
    else if(pt >= 45 && pt < 55) fakerate = FR[3];
    else if(pt >= 55 && pt < 100) fakerate = FR[4];
    else if(pt >= 100 && pt < 200) fakerate = FR[5];
    else fakerate = FR[6];
  }

  else{

    double FR[] = {0.10106,0.188843,0.189266,0.208474,0.218081,0.256777,0.36591};
    if(pt < 25) fakerate = FR[0];
    else if(pt >= 25 && pt < 35) fakerate = FR[1];
    else if(pt >= 35 && pt < 45) fakerate = FR[2];
    else if(pt >= 45 && pt < 55) fakerate = FR[3];
    else if(pt >= 55 && pt < 100) fakerate = FR[4];
    else if(pt >= 100 && pt < 200) fakerate = FR[5];
    else fakerate = FR[6];
  }

  return fakerate;
}

double FR_ratio(double pt, double eta){

  double fakerate = -999;

  if(fabs(eta) < 1.4442){

    double FR[] = {0.0615396,0.101094,0.12104,0.122836,0.142145,0.136914,0.106336};
    if(pt < 25) fakerate = FR[0];
    else if(pt >= 25 && pt < 35) fakerate = FR[1];
    else if(pt >= 35 && pt < 45) fakerate = FR[2];
    else if(pt >= 45 && pt < 55) fakerate = FR[3];
    else if(pt >= 55 && pt < 100) fakerate = FR[4];
    else if(pt >= 100 && pt < 200) fakerate = FR[5];
    else fakerate = FR[5];
  }

  else{

    double FR[] = {0.160235,0.194553,0.18197,0.182569,0.18552,0.186755,0.233082};
    if(pt < 25) fakerate = FR[0];
    else if(pt >= 25 && pt < 35) fakerate = FR[1];
    else if(pt >= 35 && pt < 45) fakerate = FR[2];
    else if(pt >= 45 && pt < 55) fakerate = FR[3];
    else if(pt >= 55 && pt < 100) fakerate = FR[4];
    else if(pt >= 100 && pt < 200) fakerate = FR[5];
    else fakerate = FR[6];
  }

  return fakerate;
}

void applyFR(int index) {

  cout<<"Chain"<<endl;
  bool mc = false;
  //PhysicsEvent* event = new PhysicsEvent();
  TChain* chain = new TChain("ntupler/ElectronTree");

  if(index==-1) {
    chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Data/SE_2015.root");
  }

  else {
    mc = true;
    if(index==1) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/DY_10to50.root");
    else if(index==2) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/DY_50toInf.root");
    else if(index==3) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/TTbar.root");
    else if(index==4) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/WJetsToLNu.root");
    else if(index==7) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/diBoson_WW.root");
    else if(index==8) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/diBoson_WZ.root");
    else if(index==9) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/diBoson_ZZ.root");
    else if(index==10) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/GammaJets_15_6000.root");
    else if(index==19) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/Single_antiTop.root");
    else if(index==20) chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/SingleTop.root");

    else {
      cout<<"Wrong input"<<endl;
      return;
    }
  }

  TFile *f1 = TFile::Open("/afs/cern.ch/user/r/rchawla/dataPUDist.root");
  TFile *f2 = TFile::Open("/afs/cern.ch/user/r/rchawla/PileUp_MC.root");

  TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
  DATA_puDist->Scale(1/DATA_puDist->Integral());

  TH1D *MC_puDist = (TH1D*)f2->Get("pileup_MC");
  TH1D *weights = (TH1D*)DATA_puDist->Clone("weights");
  weights->Divide(MC_puDist);

  Int_t           nPV;
  Int_t           nPUTrue;
  Double_t        theWeight;
  Double_t        EvtNo;
  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
  vector<float>   *chargeElec;
  vector<float>   *etaSC;
  vector<float>   *full5x5_sigmaIetaIeta;
  vector<float>   *isoRho;
  vector<int>     *passMediumId;
  vector<int>     *isPassMedium_NoSigmaEtaEta;
  vector<int>     *isPassMedium_NoPFIso;

  ptElec = 0;
  etaElec = 0;
  phiElec = 0;
  energyElec = 0;
  chargeElec = 0;
  etaSC = 0;
  full5x5_sigmaIetaIeta = 0;
  isoRho = 0;
  passMediumId = 0;
  isPassMedium_NoSigmaEtaEta = 0;
  isPassMedium_NoPFIso = 0;

  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("nPV", 1);
  chain->SetBranchStatus("nPUTrue", 1);
  chain->SetBranchStatus("theWeight", 1);
  chain->SetBranchStatus("EvtNo", 1);
  chain->SetBranchStatus("ptElec", 1);
  chain->SetBranchStatus("etaElec", 1);
  chain->SetBranchStatus("phiElec", 1);
  chain->SetBranchStatus("energyElec", 1);
  chain->SetBranchStatus("chargeElec", 1);
  chain->SetBranchStatus("etaSC", 1);
  chain->SetBranchStatus("full5x5_sigmaIetaIeta", 1);
  chain->SetBranchStatus("isoRho", 1);
  chain->SetBranchStatus("passMediumId", 1);
  chain->SetBranchStatus("isPassMedium_NoSigmaEtaEta", 1);
  chain->SetBranchStatus("isPassMedium_NoPFIso", 1);

  chain->SetBranchAddress("nPV", &nPV);
  chain->SetBranchAddress("nPUTrue", &nPUTrue);
  chain->SetBranchAddress("theWeight", &theWeight);
  chain->SetBranchAddress("EvtNo", &EvtNo);
  chain->SetBranchAddress("ptElec", &ptElec);
  chain->SetBranchAddress("etaElec", &etaElec);
  chain->SetBranchAddress("phiElec", &phiElec);
  chain->SetBranchAddress("chargeElec", &chargeElec);
  chain->SetBranchAddress("energyElec", &energyElec);
  chain->SetBranchAddress("etaSC", &etaSC);
  chain->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta);
  chain->SetBranchAddress("isoRho", &isoRho);
  chain->SetBranchAddress("passMediumId", &passMediumId);
  chain->SetBranchAddress("isPassMedium_NoSigmaEtaEta", &isPassMedium_NoSigmaEtaEta);
  chain->SetBranchAddress("isPassMedium_NoPFIso", &isPassMedium_NoPFIso);

  if(index==-1) index=6;
  TFile* file = new TFile("histograms/nocorr/fake"+TString::Itoa(index,10)+".root","RECREATE");

  int binnum = 43;
  double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

  TH1D *histDijetPFIso1_barrel = new TH1D("histDijetPFIso1_barrel", "", 500, 0, 1.0);
  TH1D *histDijetPFIso2_barrel = new TH1D("histDijetPFIso2_barrel", "", 500, 0, 1.0);
  TH1D *histDijetSigie1_barrel = new TH1D("histDijetSigie1_barrel", "", 200, 0, 0.1);
  TH1D *histDijetSigie2_barrel = new TH1D("histDijetSigie2_barrel", "", 200, 0, 0.1);

  TH1D *histDijetPFIso1_endcap = new TH1D("histDijetPFIso1_endcap", "", 500, 0, 1.0);
  TH1D *histDijetPFIso2_endcap = new TH1D("histDijetPFIso2_endcap", "", 500, 0, 1.0);
  TH1D *histDijetSigie1_endcap = new TH1D("histDijetSigie1_endcap", "", 200, 0, 0.1);
  TH1D *histDijetSigie2_endcap = new TH1D("histDijetSigie2_endcap", "", 200, 0, 0.1);

  TH1D* histDijet1 = new TH1D("histDijet1","",binnum,bins);
  TH1D* histDijet2 = new TH1D("histDijet2","",binnum,bins);
  TH1D* histDijet1_barrel = new TH1D("histDijet1_barrel","",binnum,bins);
  TH1D* histDijet2_barrel = new TH1D("histDijet2_barrel","",binnum,bins);
  TH1D* histDijet1_barend = new TH1D("histDijet1_barend","",binnum,bins);
  TH1D* histDijet2_barend = new TH1D("histDijet2_barend","",binnum,bins);
  TH1D* histDijet1_endcap = new TH1D("histDijet1_endcap","",binnum,bins);
  TH1D* histDijet2_endcap = new TH1D("histDijet2_endcap","",binnum,bins);
  TH1D* histSameDijet1 = new TH1D("histSameDijet1","",binnum,bins);
  TH1D* histSameDijet2 = new TH1D("histSameDijet2","",binnum,bins);

  TH1D* fitDijet1 = new TH1D("fitDijet1","",37,15,200);
  TH1D* fitDijet2 = new TH1D("fitDijet2","",37,15,200);
  TH1D* fitSameDijet1 = new TH1D("fitSameDijet1","",37,15,200);
  TH1D* fitSameDijet2 = new TH1D("fitSameDijet2","",37,15,200);

  TH1D* rapDijet1 = new TH1D("rapDijet1","",50,-2.5,2.5);
  TH1D* rapDijet2 = new TH1D("rapDijet2","",50,-2.5,2.5);
  TH1D* rapSameDijet1 = new TH1D("rapSameDijet1","",50,-2.5,2.5);
  TH1D* rapSameDijet2 = new TH1D("rapSameDijet2","",50,-2.5,2.5);

  histDijetPFIso1_barrel->Sumw2();
  histDijetPFIso2_barrel->Sumw2();
  histDijetSigie1_barrel->Sumw2();
  histDijetSigie2_barrel->Sumw2();

  histDijetPFIso1_endcap->Sumw2();
  histDijetPFIso2_endcap->Sumw2();
  histDijetSigie1_endcap->Sumw2();
  histDijetSigie2_endcap->Sumw2();

  histDijet1->Sumw2();
  histDijet2->Sumw2();
  histDijet1_barrel->Sumw2();
  histDijet2_barrel->Sumw2();
  histDijet1_barend->Sumw2();
  histDijet2_barend->Sumw2();
  histDijet1_endcap->Sumw2();
  histDijet2_endcap->Sumw2();
  histSameDijet1->Sumw2();
  histSameDijet2->Sumw2();

  fitDijet1->Sumw2();
  fitDijet2->Sumw2();
  fitSameDijet1->Sumw2();
  fitSameDijet2->Sumw2();

  rapDijet1->Sumw2();
  rapDijet2->Sumw2();
  rapSameDijet1->Sumw2();
  rapSameDijet2->Sumw2();

  TH1D *histWjetsPFIso1_barrel = new TH1D("histWjetsPFIso1_barrel", "", 50, 0, 0.1);
  TH1D *histWjetsPFIso2_barrel = new TH1D("histWjetsPFIso2_barrel", "", 50, 0, 0.1);
  TH1D *histWjetsSigie1_barrel = new TH1D("histWjetsSigie1_barrel", "", 200, 0, 0.1);
  TH1D *histWjetsSigie2_barrel = new TH1D("histWjetsSigie2_barrel", "", 200, 0, 0.1);

  TH1D *histWjetsPFIso1_endcap = new TH1D("histWjetsPFIso1_endcap", "", 50, 0, 0.1);
  TH1D *histWjetsPFIso2_endcap = new TH1D("histWjetsPFIso2_endcap", "", 50, 0, 0.1);
  TH1D *histWjetsSigie1_endcap = new TH1D("histWjetsSigie1_endcap", "", 200, 0, 0.1);
  TH1D *histWjetsSigie2_endcap = new TH1D("histWjetsSigie2_endcap", "", 200, 0, 0.1);

  TH1D* histWJets1 = new TH1D("histWJets1","",binnum,bins);
  TH1D* histWJets2 = new TH1D("histWJets2","",binnum,bins);
  TH1D* histWJets1_barrel = new TH1D("histWJets1_barrel","",binnum,bins);
  TH1D* histWJets2_barrel = new TH1D("histWJets2_barrel","",binnum,bins);
  TH1D* histWJets1_barend = new TH1D("histWJets1_barend","",binnum,bins);
  TH1D* histWJets2_barend = new TH1D("histWJets2_barend","",binnum,bins);
  TH1D* histWJets1_endcap = new TH1D("histWJets1_endcap","",binnum,bins);
  TH1D* histWJets2_endcap = new TH1D("histWJets2_endcap","",binnum,bins);
  TH1D* histSameWJets1 = new TH1D("histSameWJets1","",binnum,bins);
  TH1D* histSameWJets2 = new TH1D("histSameWJets2","",binnum,bins);

  TH1D* fitWJets1 = new TH1D("fitWJets1","",37,15,200);
  TH1D* fitWJets2 = new TH1D("fitWJets2","",37,15,200);
  TH1D* fitSameWJets1 = new TH1D("fitSameWJets1","",37,15,200);
  TH1D* fitSameWJets2 = new TH1D("fitSameWJets2","",37,15,200);

  TH1D* rapWJets1 = new TH1D("rapWJets1","",48,-2.4,2.4);
  TH1D* rapWJets2 = new TH1D("rapWJets2","",48,-2.4,2.4);
  TH1D* rapSameWJets1 = new TH1D("rapSameWJets1","",48,-2.4,2.4);
  TH1D* rapSameWJets2 = new TH1D("rapSameWJets2","",48,-2.4,2.4);

  TH1D* histPass = new TH1D("histPass","",10,0,10);
  TH1D* histFail = new TH1D("histFail","",10,0,10);

  histWjetsPFIso1_barrel->Sumw2();
  histWjetsPFIso2_barrel->Sumw2();
  histWjetsSigie1_barrel->Sumw2();
  histWjetsSigie2_barrel->Sumw2();

  histWjetsPFIso1_endcap->Sumw2();
  histWjetsPFIso2_endcap->Sumw2();
  histWjetsSigie1_endcap->Sumw2();
  histWjetsSigie2_endcap->Sumw2();

  histWJets1->Sumw2();
  histWJets2->Sumw2();
  histWJets1_barrel->Sumw2();
  histWJets2_barrel->Sumw2();
  histWJets1_barend->Sumw2();
  histWJets2_barend->Sumw2();
  histWJets1_endcap->Sumw2();
  histWJets2_endcap->Sumw2();
  histSameWJets1->Sumw2();
  histSameWJets2->Sumw2();

  fitWJets1->Sumw2();
  fitWJets2->Sumw2();
  fitSameWJets1->Sumw2();
  fitSameWJets2->Sumw2();

  rapWJets1->Sumw2();
  rapWJets2->Sumw2();
  rapSameWJets1->Sumw2();
  rapSameWJets2->Sumw2();

  TLorentzVector elec1, elec2, dielec;
  TLorentzVector recoElec;
  vector<double> passingElec; vector<double> failingElec;

  bool ptcut;

  double pt = 0;
  double eta = 0;
  double wt = 1.0;
  double genwt = 1.0;
  double wtsum = 0;
  bool leading = false;

  double FR1_template;
  double FR2_template;
  double FR1_ratio;
  double FR2_ratio;
  double weight_template;
  double weight_ratio;
  double mass;
  double sign;
  double rap;

  int nPass = 0;
  int nFail = 0;
  int count = 0;
  int total = 0;
  int oppsign = 0;
  int samesign = 0;

  int nentries = chain->GetEntries();
  //int nentries = 20000;
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int jentry=0; jentry < nentries; jentry++) {
    chain->GetEntry(jentry);

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    total++;

    // Sorting
    int index[ptElec->size()];
    float pt[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++) {
      pt[el]=ptElec->at(el); }

    int size = sizeof(pt)/sizeof(pt[0]);
    TMath::Sort(size,pt,index,true);

    double PUWeight = 1.0;
    int bin = 0;
    bin = weights->GetXaxis()->FindBin(nPUTrue);
    PUWeight = weights->GetBinContent(bin);

    if(mc) {
      genwt = theWeight;
      wt = theWeight*PUWeight;
      wtsum += genwt;
    }

    else wt = 1.0;

    ptcut = false;
    leading = false; 
    passingElec.clear();
    failingElec.clear();

    //printf("Event  = %f\n", EvtNo);
    //cout<<"Nele = "<<ptElec->size()<<endl;

    for(int i=0;i<ptElec->size();i++){

      Bool_t isAcc = kFALSE;
      //cout<<"Pt["<<index[i]<<"] = "<<ptElec->at(index[i])<<"   Eta["<<index[i]<<"] = "<<etaSC->at(index[i])<<endl;
      isAcc = acceptance(ptElec->at(index[i]), etaSC->at(index[i]));
      //cout<<"isAcc = "<<isAcc<<endl;

      //cout<<"isPassMedium_NoSigmaEtaEta["<<index[i]<<"] = "<<isPassMedium_NoSigmaEtaEta->at(index[i])<<"   isPassMedium_NoPFIso["<<index[i]<<"] = "<<isPassMedium_NoPFIso->at(index[i])<<endl;
      //cout<<"passMediumId["<<index[i]<<"] = "<<passMediumId->at(index[i])<<endl;

      //if(isPassMedium_NoSigmaEtaEta->at(index[i]) == 1 && isAcc)
      if((isPassMedium_NoSigmaEtaEta->at(index[i]) == 1 || isPassMedium_NoPFIso->at(index[i]) == 1) && isAcc){
	if(ptElec->at(index[i]) > 30.) leading = true;

	if(passMediumId->at(index[i]) == 1) passingElec.push_back(index[i]);
	else failingElec.push_back(index[i]);
      }
    }

    //cout<<"leading = "<<leading<<endl;

    //if(!leading) continue;

    histPass->Fill(passingElec.size(),wt);
    histFail->Fill(failingElec.size(),wt);

    nPass += passingElec.size();
    nFail += failingElec.size();

    //cout<<"passingElec = "<<passingElec.size()<<endl;
    //cout<<"failingElec = "<<failingElec.size()<<endl;

    if(failingElec.size() > 1){

      //cout<<"Nele = "<<ptElec->size()<<endl;
      //cout<<"failingElec = "<<failingElec.size()<<endl;
      //cout<<"passingElec = "<<passingElec.size()<<endl;

      if(failingElec.size() >= 2){

	//cout<<"ptElec["<<failingElec[0]<<"] = "<<ptElec->at(failingElec[0])<<"   ptElec["<<failingElec[1]<<"] = "<<ptElec->at(failingElec[1])<<endl;
	if((failingElec[0] < failingElec[1]) && ptElec->at(failingElec[0]) > 30 && ptElec->at(failingElec[1]) > 10) ptcut = true;
	else if((failingElec[0] > failingElec[1]) && ptElec->at(failingElec[1]) > 30 && ptElec->at(failingElec[0]) > 10) ptcut = true;

	if(ptcut){

	  //printf("Event  = %f\n", EvtNo);
	  count++;

	  FR1_template = FR_template(ptElec->at(failingElec[0]), etaElec->at(failingElec[0]));
	  FR2_template = FR_template(ptElec->at(failingElec[1]), etaElec->at(failingElec[1]));
	  FR1_ratio = FR_ratio(ptElec->at(failingElec[0]), etaElec->at(failingElec[0]));
	  FR2_ratio = FR_ratio(ptElec->at(failingElec[1]), etaElec->at(failingElec[1]));

	  //cout<<"weight = "<<wt<<endl;
	  //cout<<"FR1_template = "<<FR1_template<<"   "<<"FR2_template = "<<FR2_template<<endl;
	  //cout<<"FR1_ratio = "<<FR1_ratio<<"   "<<"FR2_ratio = "<<FR2_ratio<<endl;

	  weight_template = wt*FR1_template*FR2_template/((1-FR1_template)*(1-FR2_template));
	  weight_ratio = wt*FR1_ratio*FR2_ratio/((1-FR1_ratio)*(1-FR2_ratio));

	  //cout<<"weight_template = "<<weight_template<<endl;
	  //cout<<"weight_ratio = "<<weight_ratio<<endl;

	  elec1.SetPtEtaPhiE(ptElec->at(failingElec[0]), etaElec->at(failingElec[0]), phiElec->at(failingElec[0]), energyElec->at(failingElec[0]));
	  elec2.SetPtEtaPhiE(ptElec->at(failingElec[1]), etaElec->at(failingElec[1]), phiElec->at(failingElec[1]), energyElec->at(failingElec[1]));

	  dielec = elec1+elec2;
	  mass = dielec.M();
	  rap = dielec.Rapidity();
	  sign = chargeElec->at(failingElec[0])*chargeElec->at(failingElec[1]);

	  //cout<<"mass = "<<mass<<"   "<<"rapidity = "<<rap<<"   "<<"sign = "<<sign<<endl;

	  if( sign < 0 ) {
	    //oppsign++;
	    if( mass > 15 && mass < 3000) {
	      oppsign++;
	      histDijet1->Fill(mass, weight_template);
	      histDijet2->Fill(mass, weight_ratio);

	      if(etaSC->at(failingElec[0]) < 1.4442 && etaSC->at(failingElec[1]) < 1.4442){
		histDijet1_barrel->Fill(mass, weight_template);
		histDijet2_barrel->Fill(mass, weight_ratio);

		histDijetPFIso1_barrel->Fill(isoRho->at(failingElec[0]), wt);
		histDijetPFIso2_barrel->Fill(isoRho->at(failingElec[1]), wt);
		histDijetSigie1_barrel->Fill(full5x5_sigmaIetaIeta->at(failingElec[0]), wt);
		histDijetSigie2_barrel->Fill(full5x5_sigmaIetaIeta->at(failingElec[1]), wt);
	      }
	      else if((etaSC->at(failingElec[0]) < 1.4442 && etaSC->at(failingElec[1]) > 1.566) || (etaSC->at(failingElec[0]) > 1.566 && etaSC->at(failingElec[1]) < 1.4442)){
		histDijet1_barend->Fill(mass, weight_template);
		histDijet2_barend->Fill(mass, weight_ratio);
	      }
	      else if(etaSC->at(failingElec[0]) > 1.566 && etaSC->at(failingElec[1]) > 1.566){
		histDijet1_endcap->Fill(mass, weight_template);
		histDijet2_endcap->Fill(mass, weight_ratio);

		histDijetPFIso1_endcap->Fill(isoRho->at(failingElec[0]), wt);
		histDijetPFIso2_endcap->Fill(isoRho->at(failingElec[1]), wt);
		histDijetSigie1_endcap->Fill(full5x5_sigmaIetaIeta->at(failingElec[0]), wt);
		histDijetSigie2_endcap->Fill(full5x5_sigmaIetaIeta->at(failingElec[1]), wt);
	      }

	      fitDijet1->Fill(mass, weight_template);
	      fitDijet2->Fill(mass, weight_ratio);
	      rapDijet1->Fill(rap, weight_template);
	      rapDijet2->Fill(rap, weight_ratio);
	    }
	  }
	  else {
	    //samesign++;
	    if( mass > 15 && mass < 3000) {
	      samesign++;
	      histSameDijet1->Fill(mass, weight_template);
	      histSameDijet2->Fill(mass, weight_ratio);
	      fitSameDijet1->Fill(mass, weight_template);
	      fitSameDijet2->Fill(mass, weight_ratio);
	      rapSameDijet1->Fill(rap, weight_template);
	      rapSameDijet2->Fill(rap, weight_ratio);
	    }
	  }
	}
      }
    }

    else if(failingElec.size() == 1 && passingElec.size() == 1){

      //cout<<"Nele = "<<ptElec->size()<<endl;
      //cout<<"failingElec = "<<failingElec.size()<<endl;
      //cout<<"passingElec = "<<passingElec.size()<<endl;

      if((failingElec[0] < passingElec[0]) && ptElec->at(failingElec[0]) > 30 && ptElec->at(passingElec[0]) > 10) ptcut = true;
      else if((failingElec[0] > passingElec[0]) && ptElec->at(passingElec[0]) > 30 && ptElec->at(failingElec[0]) > 10) ptcut = true;

      //cout<<"ptElec["<<failingElec[0]<<"] = "<<ptElec->at(failingElec[0])<<"   ptElec["<<passingElec[0]<<"] = "<<ptElec->at(passingElec[0])<<endl;
      if(ptcut){

	FR1_template = FR_template(ptElec->at(failingElec[0]), etaElec->at(failingElec[0]));
	FR1_ratio = FR_ratio(ptElec->at(failingElec[0]), etaElec->at(failingElec[0]));

	//cout<<"weight = "<<wt<<endl;
	//cout<<"FR1_template = "<<FR1_template<<endl;
	//cout<<"FR1_ratio = "<<FR1_ratio<<endl;

	weight_template = wt*FR1_template/(1-FR1_template);
	weight_ratio = wt*FR1_ratio/(1-FR1_ratio);

	//cout<<"weight_template = "<<weight_template<<endl;
	//cout<<"weight_ratio = "<<weight_ratio<<endl;

	elec1.SetPtEtaPhiE(ptElec->at(failingElec[0]), etaElec->at(failingElec[0]), phiElec->at(failingElec[0]), energyElec->at(failingElec[0]));
	elec2.SetPtEtaPhiE(ptElec->at(passingElec[0]), etaElec->at(passingElec[0]), phiElec->at(passingElec[0]), energyElec->at(passingElec[0]));

	dielec = elec1+elec2;
	mass = dielec.M();
	rap = dielec.Rapidity();
	sign = chargeElec->at(failingElec[0])*chargeElec->at(passingElec[0]);

	//cout<<"mass = "<<mass<<"   "<<"rapidity = "<<rap<<"   "<<"sign = "<<sign<<endl;

	if( sign < 0 ) {
	  if( mass > 15 && mass < 3000) {
	    histWJets1->Fill(mass, weight_template);
	    histWJets2->Fill(mass, weight_ratio);

	    if(etaSC->at(failingElec[0]) < 1.4442 && etaSC->at(passingElec[0]) < 1.4442){
	      histWJets1_barrel->Fill(mass, weight_template);
	      histWJets2_barrel->Fill(mass, weight_ratio);

	      histWjetsPFIso1_barrel->Fill(isoRho->at(failingElec[0]), wt);
	      histWjetsPFIso2_barrel->Fill(isoRho->at(passingElec[0]), wt);
	      histWjetsSigie1_barrel->Fill(full5x5_sigmaIetaIeta->at(failingElec[0]), wt);
	      histWjetsSigie2_barrel->Fill(full5x5_sigmaIetaIeta->at(passingElec[0]), wt);
	    }
	    else if((etaSC->at(failingElec[0]) < 1.4442 && etaSC->at(passingElec[0]) > 1.566) || (etaSC->at(failingElec[0]) > 1.566 && etaSC->at(passingElec[0]) < 1.4442)){
	      histWJets1_barend->Fill(mass, weight_template);
	      histWJets2_barend->Fill(mass, weight_ratio);
	    }
	    else if(etaSC->at(failingElec[0]) > 1.566 && etaSC->at(passingElec[0]) > 1.566){
	      histWJets1_endcap->Fill(mass, weight_template);
	      histWJets2_endcap->Fill(mass, weight_ratio);

	      histWjetsPFIso1_endcap->Fill(isoRho->at(failingElec[0]), wt);
	      histWjetsPFIso2_endcap->Fill(isoRho->at(passingElec[0]), wt);
	      histWjetsSigie1_endcap->Fill(full5x5_sigmaIetaIeta->at(failingElec[0]), wt);
	      histWjetsSigie2_endcap->Fill(full5x5_sigmaIetaIeta->at(passingElec[0]), wt);
	    }

	    fitWJets1->Fill(mass, weight_template);
	    fitWJets2->Fill(mass, weight_ratio);
	    rapWJets1->Fill(rap, weight_template);
	    rapWJets2->Fill(rap, weight_ratio);
	  }
	}

	else {
	  if( mass > 15 && mass < 3000) {
	    histSameWJets1->Fill(mass, weight_template);
	    histSameWJets2->Fill(mass, weight_ratio);
	    fitSameWJets1->Fill(mass, weight_template);
	    fitSameWJets2->Fill(mass, weight_ratio);
	    rapSameWJets1->Fill(rap, weight_template);
	    rapSameWJets2->Fill(rap, weight_ratio);
	  }
	}
      }
    }

    //cout<<""<<endl;

  } // event

  cout<<"# of passing electrons = "<<nPass<<endl;
  cout<<"# of failing electrons = "<<nFail<<endl;
  cout<<endl;
  //cout<<"# of passing electrons per event = "<<nPass/wtsum<<endl;
  //cout<<"# of failing electrons per event = "<<nFail/wtsum<<endl;
  //cout<<endl;

  cout<<"Success"<<endl;
  cout<<wtsum<<endl;

  cout<<"# of events with opposite sign electrons passing dijet selection = "<<oppsign<<endl;
  cout<<"# of events with same sign electrons passing dijet selection = "<<samesign<<endl;
  cout<<"# of events passing dijet selection = "<<count<<endl;
  cout<<"# of events = "<<total<<endl;

  file->Write();
  file->Close();

} // applyFR
