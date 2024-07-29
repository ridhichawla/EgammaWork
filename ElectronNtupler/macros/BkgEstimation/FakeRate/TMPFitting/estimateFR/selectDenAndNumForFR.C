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

void selectDenAndNumForFR(int index) {

  cout<<"Chain"<<endl;
  bool mc = false;

  TChain* chain = new TChain("ntupler/ElectronTree");

  if(index==-1) {
    chain->Add("/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Data/Photon_2015.root");
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
    else if(index==11) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt15to20.root");
    else if(index==12) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt20to30.root");
    else if(index==13) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt30to50.root");
    else if(index==14) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt50to80.root");
    else if(index==15) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt80to120.root");
    else if(index==16) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt120to170.root");
    else if(index==17) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt170to300.root");
    else if(index==18) chain->Add("/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_16082017/Backgrounds/QCD_Pt300toInf.root");

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

  Double_t        EvtNo;
  Int_t           singlePhoton;
  Int_t           prescalePhoton;
  Int_t           nPV;
  Int_t           nPUTrue;
  Double_t        theWeight;
  vector<float>   *ptElec;
  vector<float>   *etaElec;
  vector<float>   *phiElec;
  vector<float>   *energyElec;
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
  etaSC = 0;
  full5x5_sigmaIetaIeta = 0;
  isoRho = 0;
  passMediumId = 0;
  isPassMedium_NoSigmaEtaEta = 0;
  isPassMedium_NoPFIso = 0;

  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("EvtNo", 1);
  chain->SetBranchStatus("singlePhoton", 1);
  chain->SetBranchStatus("prescalePhoton", 1);
  chain->SetBranchStatus("nPV", 1);
  chain->SetBranchStatus("nPUTrue", 1);
  chain->SetBranchStatus("theWeight", 1);
  chain->SetBranchStatus("ptElec", 1);
  chain->SetBranchStatus("etaElec", 1);
  chain->SetBranchStatus("phiElec", 1);
  chain->SetBranchStatus("energyElec", 1);
  chain->SetBranchStatus("etaSC", 1);
  chain->SetBranchStatus("full5x5_sigmaIetaIeta", 1);
  chain->SetBranchStatus("isoRho", 1);
  chain->SetBranchStatus("passMediumId", 1);
  chain->SetBranchStatus("isPassMedium_NoSigmaEtaEta", 1);
  chain->SetBranchStatus("isPassMedium_NoPFIso", 1);

  chain->SetBranchAddress("EvtNo", &EvtNo);
  chain->SetBranchAddress("singlePhoton", &singlePhoton);
  chain->SetBranchAddress("prescalePhoton", &prescalePhoton);
  chain->SetBranchAddress("nPV", &nPV);
  chain->SetBranchAddress("nPUTrue", &nPUTrue);
  chain->SetBranchAddress("theWeight", &theWeight);
  chain->SetBranchAddress("ptElec", &ptElec);
  chain->SetBranchAddress("etaElec", &etaElec);
  chain->SetBranchAddress("phiElec", &phiElec);
  chain->SetBranchAddress("energyElec", &energyElec);
  chain->SetBranchAddress("etaSC", &etaSC);
  chain->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta);
  chain->SetBranchAddress("isoRho", &isoRho);
  chain->SetBranchAddress("passMediumId", &passMediumId);
  chain->SetBranchAddress("isPassMedium_NoSigmaEtaEta", &isPassMedium_NoSigmaEtaEta);
  chain->SetBranchAddress("isPassMedium_NoPFIso", &isPassMedium_NoPFIso);

  if(index==-1) index=6;
  TFile* file = new TFile("histograms/hist"+TString::Itoa(index,10)+".root","RECREATE");

  double Pt, Eta, SCEta, Sigma_eta, Isolation, genwt, wt, wtsum;
  vector <double> passingElectron;

  //Double_t x1bin[6] = {10,20,30,40,50,10000};
  Double_t x1bin[8] = {10,25,35,45,55,100,200,1000};
  int nbins = 7;

  TH1D *numerator_pt          = new TH1D("numerator_pt", "numerator_pt", nbins, x1bin);
  TH1D *denominator_pt        = new TH1D("denominator_pt", "denominator_pt", nbins, x1bin);
  TH1D *numerator_pt_barrel   = new TH1D("numerator_pt_barrel", "numerator_pt_barrel", nbins, x1bin);
  TH1D *denominator_pt_barrel = new TH1D("denominator_pt_barrel", "denominator_pt_barrel", nbins, x1bin);
  TH1D *numerator_pt_endcap   = new TH1D("numerator_pt_endcap", "numerator_pt_endcap", nbins, x1bin);
  TH1D *denominator_pt_endcap = new TH1D("denominator_pt_endcap", "denominator_pt_endcap", nbins, x1bin);

  TH1D *numerator_eta      = new TH1D("numerator_eta", "numerator_eta", 60, -3, 3);
  TH1D *denominator_eta    = new TH1D("denominator_eta", "denominator_eta", 60, -3, 3);

  TH1D *numerator          = new TH1D("numerator", "numerator", 200, 0, 0.1);
  TH1D *denominator        = new TH1D("denominator", "denominator", 200, 0, 0.1);
  TH1D *numerator_barrel   = new TH1D("numerator_barrel", "numerator_barrel", 100, 0, 0.05);
  TH1D *denominator_barrel = new TH1D("denominator_barrel", "denominator_barrel", 100, 0, 0.05);
  TH1D *numerator_endcap   = new TH1D("numerator_endcap", "numerator_endcap", 100, 0.01, 0.06);
  TH1D *denominator_endcap = new TH1D("denominator_endcap", "denominator_endcap", 100, 0.01, 0.06);
  //TH1D *numerator_barrel   = new TH1D("numerator_barrel", "numerator_barrel", 100, 0, 0.1);
  //TH1D *denominator_barrel = new TH1D("denominator_barrel", "denominator_barrel", 100, 0, 0.1);
  //TH1D *numerator_endcap   = new TH1D("numerator_endcap", "numerator_endcap", 100, 0, 0.1);
  //TH1D *denominator_endcap = new TH1D("denominator_endcap", "denominator_endcap", 100, 0, 0.1);

  numerator_pt->Sumw2(); denominator_pt->Sumw2(); numerator_eta->Sumw2(); denominator_eta->Sumw2();
  numerator_pt_barrel->Sumw2(); denominator_pt_barrel->Sumw2();
  numerator_pt_endcap->Sumw2(); denominator_pt_endcap->Sumw2();

  numerator->Sumw2(); numerator_barrel->Sumw2(); numerator_endcap->Sumw2();
  denominator->Sumw2(); denominator_barrel->Sumw2(); denominator_endcap->Sumw2();

  wtsum = 0.;
  
  int nentries = chain->GetEntries();
  //int nentries = 100;
  cout<<"entries: "<<nentries<<endl;
  for (unsigned int jentry=0; jentry < nentries; jentry++) {
    chain->GetEntry(jentry);

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    // Sorting
    int index[ptElec->size()];
    float pt[ptElec->size()];

    for(unsigned int el=0; el<ptElec->size(); el++) {
      pt[el]=ptElec->at(el); }

    int size = sizeof(pt)/sizeof(pt[0]);
    TMath::Sort(size,pt,index,true);

    Pt = 0.;
    Eta = 0.;
    SCEta = 0.;
    Sigma_eta = 0.;
    Isolation = 0.;
    genwt = 0.;
    wt = 0.;
    passingElectron.clear();

    double PUWeight = 1.0;
    int bin = 0;
    bin = weights->GetXaxis()->FindBin(nPUTrue);
    PUWeight = weights->GetBinContent(bin);

    if(mc) {
      genwt = theWeight;
      wt = theWeight*PUWeight;
      wtsum += genwt;

      //cout<<"gen weight = "<<genwt<<"   "<<"sum = "<<wtsum<<endl;
    }
    else wt = prescalePhoton;

    if(!singlePhoton) continue;

    //printf("Event  = %f\n",EvtNo);
    //cout<<"Nele = "<<ptElec->size()<<endl;

    for(int i=0;i<ptElec->size();i++){

      Bool_t isAcc = kFALSE;
      //cout<<"Pt["<<i<<"] = "<<ptElec->at(index[i])<<"   Eta["<<i<<"] = "<<etaSC->at(index[i])<<endl;
      isAcc = acceptance(ptElec->at(index[i]), etaSC->at(index[i]));
      //cout<<"isAcc = "<<isAcc<<endl;

      if(isAcc){
	passingElectron.push_back(index[i]);
      }
    }

    //cout<<"Npassing = "<<passingElectron.size()<<endl;

    for(int j=0;j<passingElectron.size();j++){

      //cout<<"passingElectron ["<<j<<"] = "<<passingElectron.at(j)<<endl;

      Pt = ptElec->at(passingElectron[j]);
      Eta = etaElec->at(passingElectron[j]);
      SCEta = etaSC->at(passingElectron[j]);
      Sigma_eta = full5x5_sigmaIetaIeta->at(passingElectron[j]);
      Isolation = isoRho->at(passingElectron[j]);

      //cout<<"Pt["<<j<<"] = "<<Pt<<"   SC Eta["<<j<<"] = "<<SCEta<<"   Sigma_eta["<<j<<"] = "<<Sigma_eta<<"   Isolation["<<j<<"] = "<<Isolation<<endl;

      //cout<<"Medium ID with sigma eta masked["<<j<<"] = "<<isPassMedium_NoSigmaEtaEta->at(passingElectron[j])<<"   Medium ID with isolation masked["<<j<<"] = "<<isPassMedium_NoPFIso->at(passingElectron[j])<<endl;

      if(isPassMedium_NoSigmaEtaEta->at(passingElectron[j]) == 1 || isPassMedium_NoPFIso->at(passingElectron[j]) == 1){

	//cout<<"Denominator   Pt["<<j<<"] = "<<Pt<<"   Eta["<<j<<"] = "<<Eta<<"   Sigma_eta["<<j<<"] = "<<Sigma_eta<<"   Isolation["<<j<<"] = "<<Isolation<<endl; 
	denominator_pt->Fill(Pt,wt);
	denominator_eta->Fill(Eta,wt);
	denominator->Fill(Sigma_eta,wt);

	if(fabs(SCEta) < 1.4442){
	  denominator_pt_barrel->Fill(Pt,wt);
	  denominator_barrel->Fill(Sigma_eta,wt);
	}

	else{
	  denominator_pt_endcap->Fill(Pt,wt);
	  denominator_endcap->Fill(Sigma_eta,wt);
	}

	//cout<<"Medium ID["<<j<<"] = "<<passMediumId->at(passingElectron[j])<<endl;

	if(passMediumId->at(passingElectron[j]) == 0) continue;

	//cout<<"Numerator   Pt["<<j<<"] = "<<Pt<<"   SC Eta["<<j<<"] = "<<SCEta<<"   Sigma_eta["<<j<<"] = "<<Sigma_eta<<endl;

	numerator_pt->Fill(Pt,wt);
	numerator_eta->Fill(Eta,wt);
	numerator->Fill(Sigma_eta,wt);

	if(fabs(SCEta) < 1.4442){
	  numerator_pt_barrel->Fill(Pt,wt);
	  numerator_barrel->Fill(Sigma_eta,wt);
	}

	else{
	  numerator_pt_endcap->Fill(Pt,wt);
	  numerator_endcap->Fill(Sigma_eta,wt);
	}

      }

    }

    //cout<<""<<endl;

  }

  file->Write();
  file->Close();

  cout<<"Sum of weights = "<<wtsum<<endl;
}
