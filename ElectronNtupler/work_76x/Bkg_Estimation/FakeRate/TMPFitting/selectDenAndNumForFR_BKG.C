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

bool acceptance(double pt_, double eta_) {
  if( pt_ > 10. && fabs(eta_) < 2.5 && (fabs(eta_) < 1.4442 || fabs(eta_) > 1.566) ) return true;
  else return false;
}

void selectDenAndNumForFR_BKG() {

  TString workdir;
  std::vector<TFile*> InputFiles_bkg;

  const char *bkg[8] = {"WJetsToLNu", "TTbar", "diBoson_WW", "diBoson_WZ", "diBoson_ZZ", "Single_antiTop", "SingleTop", "QCD"};

  double xsec[8] = {61526.7,831.76,118.7,47.13,16.523,35.6,35.6,162060000.};
  double noEvts[8] = {3731926637458.121094,97994304.,988416.,999996.,985598.,999399.,999999.,37835400.};

  workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/";
  //workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated_19042017/";

  InputFiles_bkg.clear();

  InputFiles_bkg.push_back(TFile::Open(workdir+"WJetsToLNu.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"TTbar.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WW.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_ZZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"Single_antiTop.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"SingleTop.root"));
  //InputFiles_bkg.push_back(TFile::Open(workdir+"QCD.root"));

  int nsample = InputFiles_bkg.size();
  TFile* file[8];

  for(unsigned int jentry = 0; jentry < 1; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_bkg.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("/afs/cern.ch/user/r/rchawla/dataPUDist.root");
    TFile *f2 = TFile::Open("/afs/cern.ch/user/r/rchawla/PileUp_MC.root");

    TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    TH1D *MC_puDist = (TH1D*)f2->Get("pileup_MC");
    TH1D *weights = (TH1D*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    Int_t           singlePhoton;
    Int_t           prescalePhoton;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    vector<float>   *full5x5_sigmaIetaIeta;
    vector<int>     *passMediumId;
    vector<int>     *isPassMedium_NoSigmaEtaEta;
    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Rap;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    Int_t           tauFlag;
    Double_t        theWeight;
    Int_t           nPV;
    Int_t           nPUTrue;

    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    full5x5_sigmaIetaIeta = 0;
    passMediumId = 0;
    isPassMedium_NoSigmaEtaEta = 0;
    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Rap = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("singlePhoton", 1);
    T1->SetBranchStatus("prescalePhoton", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("full5x5_sigmaIetaIeta", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("isPassMedium_NoSigmaEtaEta", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Rap", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("nPV", 1);
    T1->SetBranchStatus("nPUTrue", 1);
    T1->SetBranchStatus("theWeight", 1);

    T1->SetBranchAddress("singlePhoton", &singlePhoton);
    T1->SetBranchAddress("prescalePhoton", &prescalePhoton);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("isPassMedium_NoSigmaEtaEta", &isPassMedium_NoSigmaEtaEta);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("nPV", &nPV);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("%s.root",bkg[jentry]),"RECREATE");

    double lumi_Weight = xsec[jentry]/noEvts[jentry];
    cout<<"Background Sample: "<<bkg[jentry]<<endl;

    double Pt, Eta, Sigma_eta;
    vector <double> passingElectron;

    Double_t x1bin[6] = {10,20,30,40,50,10000};
    int nbins = 5;

    TH1D *numerator_pt      = new TH1D("numerator_pt", "numerator_pt", nbins, x1bin);
    TH1D *numerator_pt_barrel  = new TH1D("numerator_pt_barrel", "numerator_pt_barrel", nbins, x1bin);
    TH1D *numerator_pt_endcap = new TH1D("numerator_pt_endcap", "numerator_pt_endcap", nbins, x1bin);

    TH1D *denominator_pt      = new TH1D("denominator_pt", "denominator_pt", nbins, x1bin);
    TH1D *denominator_pt_barrel  = new TH1D("denominator_pt_barrel", "denominator_pt_barrel", nbins, x1bin);
    TH1D *denominator_pt_endcap = new TH1D("denominator_pt_endcap", "denominator_pt_endcap", nbins, x1bin);

    TH1D *numerator_eta     = new TH1D("numerator_eta", "numerator_eta", 60, -3, 3);
    TH1D *denominator_eta     = new TH1D("denominator_eta", "denominator_eta", 60, -3, 3);

    TH1D *numerator        = new TH1D("numerator", "numerator", 1000, 0, 1);
    TH1D *numerator_barrel    = new TH1D("numerator_barrel", "numerator_barrel", 1000, 0, 1);
    TH1D *numerator_endcap   = new TH1D("numerator_endcap", "numerator_endcap", 1000, 0, 1);

    TH1D *denominator        = new TH1D("denominator", "denominator", 1000, 0, 1);
    TH1D *denominator_barrel    = new TH1D("denominator_barrel", "denominator_barrel", 1000, 0, 1);
    TH1D *denominator_endcap   = new TH1D("denominator_endcap", "denominator_endcap", 1000, 0, 1);

    numerator_pt->Sumw2(); denominator_pt->Sumw2(); numerator_eta->Sumw2(); denominator_eta->Sumw2();
    numerator_pt_barrel->Sumw2(); denominator_pt_barrel->Sumw2();
    numerator_pt_endcap->Sumw2(); denominator_pt_endcap->Sumw2();

    numerator->Sumw2(); numerator_barrel->Sumw2(); numerator_endcap->Sumw2();
    denominator->Sumw2(); denominator_barrel->Sumw2(); denominator_endcap->Sumw2();

    int nentries = T1->GetEntries();
    //int nentries = 2000000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int jentry=908406; jentry < 908420; jentry++) {
      T1->GetEntry(jentry);

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
      Sigma_eta = 0.;
      passingElectron.clear();

      double PUWeight = 1.0;

      int bin = 0;
      double puWeights = 1.0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      PUWeight = weights->GetBinContent(bin);

      if(!singlePhoton) continue;

      cout<<"Nele = "<<ptElec->size()<<endl;

      for(int i=0;i<ptElec->size();i++){

	Bool_t isAcc = kFALSE;
	cout<<"Pt["<<i<<"] = "<<ptElec->at(index[i])<<"   Eta["<<i<<"] = "<<etaElec->at(index[i])<<endl;
	isAcc = acceptance(ptElec->at(index[i]), etaElec->at(index[i]));
	cout<<"isAcc = "<<isAcc<<endl;

	if(isAcc){
	  passingElectron.push_back(index[i]);
	}
      }

      cout<<"Npassing = "<<passingElectron.size()<<endl;

      for(int j=0;j<passingElectron.size();j++){

	cout<<"passingElectron ["<<j<<"] = "<<passingElectron.at(j)<<endl;

	Pt = ptElec->at(passingElectron[j]);
	Eta = etaElec->at(passingElectron[j]);
	Sigma_eta = full5x5_sigmaIetaIeta->at(passingElectron[j]);

	cout<<"Pt["<<j<<"] = "<<Pt<<"   Eta["<<j<<"] = "<<Eta<<"   sigma_ieta_ieta["<<j<<"] = "<<Sigma_eta<<endl;

	//cout<<"Medium ID["<<j<<"] = "<<passMediumId->at(passingElectron[j])<<endl;
	cout<<"Medium ID with sigma eta masked["<<j<<"] = "<<isPassMedium_NoSigmaEtaEta->at(passingElectron[j])<<endl;

	if(isPassMedium_NoSigmaEtaEta->at(passingElectron[j]) == 0) continue;

	cout<<"Denominator-> Pt["<<j<<"] = "<<Pt<<"   Eta["<<j<<"] = "<<Eta<<"   sigma_ieta_ieta["<<j<<"] = "<<Sigma_eta<<endl; 
	denominator_pt->Fill(Pt,lumi_Weight*PUWeight*2258.066*theWeight);
	denominator_eta->Fill(Eta,lumi_Weight*PUWeight*2258.066*theWeight);
	denominator->Fill(Sigma_eta,lumi_Weight*PUWeight*2258.066*theWeight);

	if(fabs(Eta) < 1.4442){
	  denominator_pt_barrel->Fill(Pt,lumi_Weight*PUWeight*2258.066*theWeight);
	  denominator_barrel->Fill(Sigma_eta,lumi_Weight*PUWeight*2258.066*theWeight);
	}

	else{
	  denominator_pt_endcap->Fill(Pt,lumi_Weight*PUWeight*2258.066*theWeight);
	  denominator_endcap->Fill(Sigma_eta,lumi_Weight*PUWeight*2258.066*theWeight);
	}

	cout<<"Medium ID with sigma eta unmasked["<<j<<"] = "<<passMediumId->at(passingElectron[j])<<endl;

	if(passMediumId->at(passingElectron[j]) == 0) continue;

	cout<<"Numerator->   Pt["<<j<<"] = "<<Pt<<"   Eta["<<j<<"] = "<<Eta<<"   sigma_ieta_ieta["<<j<<"] = "<<Sigma_eta<<endl;

	numerator_pt->Fill(Pt,lumi_Weight*PUWeight*2258.066*theWeight);
	numerator_eta->Fill(Eta,lumi_Weight*PUWeight*2258.066*theWeight);
	numerator->Fill(Sigma_eta,lumi_Weight*PUWeight*2258.066*theWeight);

	if(fabs(Eta) < 1.4442){
	  numerator_pt_barrel->Fill(Pt,lumi_Weight*PUWeight*2258.066*theWeight);
	  numerator_barrel->Fill(Sigma_eta,lumi_Weight*PUWeight*2258.066*theWeight);

	  if(Sigma_eta >= 0.0101) cout<<"Entry = "<<jentry<<"   sigma_ieta_ieta = "<<Sigma_eta<<"   eta = "<<Eta<<endl;
	  //cout<<""<<endl;
	}

	else{
	  numerator_pt_endcap->Fill(Pt,lumi_Weight*PUWeight*2258.066*theWeight);
	  numerator_endcap->Fill(Sigma_eta,lumi_Weight*PUWeight*2258.066*theWeight);
	}

      }

      cout<<""<<endl;

    }

    file[jentry]->Write();
    file[jentry]->Close();

    cout<<""<<endl;
  } // file Loop  
}
