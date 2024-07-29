#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void BKG() {

  TString workdir;
  std::vector<TFile*> InputFiles_bkg;
  const char *bkg[3] = {"GammaJets", "WJetsToLNu", "WGamma"};

  //workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated_19042017/";
  workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated_18072017/";
  InputFiles_bkg.clear();

  InputFiles_bkg.push_back(TFile::Open(workdir+"GammaJets_15_6000.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"WJetsToLNu.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"WGamma.root"));

  int nsample = InputFiles_bkg.size();

  for(unsigned int jentry = 2; jentry < nsample; ++jentry){
    TTree * T1 = (TTree*)InputFiles_bkg.at(jentry)->Get("ntupler/ElectronTree");

    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    Double_t        theWeight;
    Int_t           tauFlag;

    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("tauFlag", 1);

    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("tauFlag", &tauFlag);

    double sum1_weights;
    double sum2_weights;
    double sum3_weights;
    double sum4_weights;
    double sum5_weights;
    double sum6_weights;
    double sum7_weights;
    double sum8_weights;

    sum1_weights = 0.;
    sum2_weights = 0.;
    sum3_weights = 0.;
    sum4_weights = 0.;
    sum5_weights = 0.;
    sum6_weights = 0.;
    sum7_weights = 0.;
    sum8_weights = 0.;

    cout<<"Background Sample: "<<bkg[jentry]<<endl;
    double massGPre, diLHE_Mass;
    TLorentzVector genPre1,genPre2,diGen;
    TLorentzVector lhe1,lhe2,diLHE;
    vector <int> idx1;

    int nentries = T1->GetEntries();
    //int nentries = 7000000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      idx1.clear();
      diLHE_Mass=0.;
      massGPre=0.;

      // Sum of weights for nominal samples
      if(genPreFSR_Pt->size() == 2){
	genPre1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	genPre2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=genPre1+genPre2;
	massGPre=diGen.M();
      }

      sum1_weights = sum1_weights + theWeight;
      if(jentry==1 && massGPre < 100.) sum2_weights = sum2_weights + theWeight;
      if(jentry==1 && massGPre < 100. && genPreFSR_Pt->size() == 2) sum3_weights = sum3_weights + theWeight;
      if(jentry==1 && massGPre < 100. && !tauFlag) sum4_weights = sum4_weights + theWeight;
      if(jentry==1 && massGPre < 100. && !tauFlag && genPreFSR_Pt->size() == 2) sum5_weights = sum5_weights + theWeight;
      if(genPreFSR_Pt->size() == 2) sum6_weights = sum6_weights + theWeight;
      if(!tauFlag) sum7_weights = sum7_weights + theWeight;
      if(!tauFlag && genPreFSR_Pt->size() == 2) sum8_weights = sum8_weights + theWeight;


    } // event Loop

    printf ("sum of weights No condition: %f \n",sum1_weights);
    printf ("sum of weights with mass condition: %f \n",sum2_weights);
    printf ("sum of weights with mass+isHardProcess condition: %f \n",sum3_weights);
    printf ("sum of weights with mass+noTau condition: %f \n",sum4_weights);
    printf ("sum of weights with mass+isHardProcess+noTau condition: %f \n",sum5_weights);
    printf ("sum of weights with isHardProcess condition: %f \n",sum6_weights);
    printf ("sum of weights with noTau condition: %f \n",sum7_weights);
    printf ("sum of weights with isHardProcess+noTau condition: %f \n",sum8_weights);
    cout<<""<<endl;

  } // file Loop
}
