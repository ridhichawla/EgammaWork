#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void SumofWeights() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  int mass[13] = {10,50,50,100,200,400,500,700,800,1000,1500,2000,3000};
  //int mass[8] = {10,50,50,50,50,50,50,3000};
  //int mass[7] = {1,2,3,4,5,6,7};

  //workdir = "/tmp/rchawla/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";
  //workdir = "/tmp/rchawla/eos/cms/store/user/arun/DY_SmearSyst/";
  //workdir = "/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/PHOTOS_FSRCorr/";
  workdir = "/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/PHOTOS_FSRCorr/v201611208_1st_ChangeMassScale_HighMassSample/";
  InputFiles_signal_DY.clear();

  //InputFiles_signal_DY.push_back(TFile::Open(workdir+"MG_M5_50.root"));
  /*InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_10to50.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_50toInf.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_100to200.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_200to400.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_400to500.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_500to700.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_700to800.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_800to1000.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1000to1500.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_1500to2000.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_2000to3000.root"));*/

  /*InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_1.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_2.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_3.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_4.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_5.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_6.root"));
    InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M50toInf_Photos.root"));*/

  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M10to50/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M50toInf_part1/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M50toInf_part2/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M100to200/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M200to400/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M400to500/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M500to700/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M700to800/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M800to1000/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M1000to1500/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M1500to2000/ntuple_skim.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"M2000to3000/ntuple_skim.root"));

  int nsample = InputFiles_signal_DY.size();

  for(unsigned int jentry = 0; jentry < nsample; ++jentry){
    //TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("ntupler/ElectronTree");
    //TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("DYTree");
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("recoTree/DYTree");

    /*vector<float>   *genPreFSR_Pt;
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
      T1->SetBranchStatus("tauFlag", 1);
      T1->SetBranchStatus("theWeight", 1);

      T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
      T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
      T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
      T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
      T1->SetBranchAddress("theWeight", &theWeight);
      T1->SetBranchAddress("tauFlag", &tauFlag);*/

    //cout<<"1"<<endl;

    static const int MPSIZE = 2000;

    int GENnPair;
    int nLHEParticle;
    double LHELepton_Px[MPSIZE];
    double LHELepton_Py[MPSIZE];
    double LHELepton_Pz[MPSIZE];
    double LHELepton_E[MPSIZE];
    int LHELepton_ID[MPSIZE];
    int LHELepton_status[MPSIZE];
    double GENLepton_pT[MPSIZE];
    int GENLepton_ID[MPSIZE];
    int GENLepton_fromHardProcessFinalState[MPSIZE];
    double GENEvt_weight;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("GENnPair", 1);
    T1->SetBranchStatus("nLHEParticle", 1);
    T1->SetBranchStatus("LHELepton_Px", 1);
    T1->SetBranchStatus("LHELepton_Py", 1);
    T1->SetBranchStatus("LHELepton_Pz", 1);
    T1->SetBranchStatus("LHELepton_E", 1);
    T1->SetBranchStatus("LHELepton_ID", 1);
    T1->SetBranchStatus("LHELepton_status", 1);
    T1->SetBranchStatus("GENLepton_pT", 1);
    T1->SetBranchStatus("GENLepton_ID", 1);
    T1->SetBranchStatus("GENLepton_fromHardProcessFinalState", 1);
    T1->SetBranchStatus("GENEvt_weight", 1);

    T1->SetBranchAddress("GENnPair", &GENnPair);
    T1->SetBranchAddress("nLHEParticle", &nLHEParticle);
    T1->SetBranchAddress("LHELepton_Px", &LHELepton_Px);
    T1->SetBranchAddress("LHELepton_Py", &LHELepton_Py);
    T1->SetBranchAddress("LHELepton_Pz", &LHELepton_Pz);
    T1->SetBranchAddress("LHELepton_E", &LHELepton_E);
    T1->SetBranchAddress("LHELepton_ID", &LHELepton_ID);
    T1->SetBranchAddress("LHELepton_status", &LHELepton_status);
    T1->SetBranchAddress("GENLepton_pT", &GENLepton_pT);
    T1->SetBranchAddress("GENLepton_ID", &GENLepton_ID);
    T1->SetBranchAddress("GENLepton_fromHardProcessFinalState", &GENLepton_fromHardProcessFinalState);
    T1->SetBranchAddress("GENEvt_weight", &GENEvt_weight);

    double sum1_weights;
    double sum2_weights;
    double sum3_weights;

    sum1_weights = 0.;
    sum2_weights = 0.;
    sum3_weights = 0.;

    cout<<"Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;
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
      
      // Sum of weights for PHOTOS high mass samples
      for(int i=0; i<nLHEParticle; i++){
	if(abs(LHELepton_ID[i]) == 11) idx1.push_back(i);
      }

      if(idx1.size() == 2){

	if(abs(LHELepton_ID[idx1.at(0)]) == 11 && abs(LHELepton_ID[idx1.at(1)]) == 11 && LHELepton_status[idx1.at(0)] == 1 && LHELepton_status[idx1.at(1)] == 1){
	  lhe1.SetPxPyPzE(LHELepton_Px[idx1.at(0)],LHELepton_Py[idx1.at(0)],LHELepton_Pz[idx1.at(0)],LHELepton_E[idx1.at(0)]);
	  lhe2.SetPxPyPzE(LHELepton_Px[idx1.at(1)],LHELepton_Py[idx1.at(1)],LHELepton_Pz[idx1.at(1)],LHELepton_E[idx1.at(1)]);

	  diLHE=lhe1+lhe2;
	  diLHE_Mass=diLHE.M();
	}
      }

      sum1_weights = sum1_weights + GENEvt_weight;
      if(idx1.size() == 2) sum2_weights = sum2_weights + GENEvt_weight;
      if((jentry==1 || jentry==2) && diLHE_Mass < 100. && idx1.size() == 2) sum3_weights = sum3_weights + GENEvt_weight;


      // Sum of weights for PHOTOS samples (M10-50 and M50-Inf)
      /*sum1_weights = sum1_weights + GENEvt_weight;

	for(int i=0; i<GENnPair; i++){
	if(abs(GENLepton_ID[i]) == 11 && GENLepton_fromHardProcessFinalState[i] == 1) idx1.push_back(i);
	}

	if(idx1.size() == 2){
	sum2_weights = sum2_weights + GENEvt_weight;
	}*/


      // Sum of weights for nominal samples
      /*if(genPreFSR_Pt->size() == 2){
	genPre1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	genPre2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=genPre1+genPre2;
	massGPre=diGen.M();
	}

	sum1_weights = sum1_weights + theWeight;

	if(!tauFlag && genPreFSR_Pt->size() == 2) {
	sum2_weights = sum2_weights + theWeight;

	if(jentry==1 && massGPre < 100.) sum3_weights = sum3_weights + theWeight;
	}*/

    } // event Loop

    printf ("sum of weights No condition: %f \n",sum1_weights);
    printf ("sum of weights with condition: %f \n",sum2_weights);
    printf ("sum of weights with mass condition: %f \n",sum3_weights);
    cout<<""<<endl;

  } // file Loop
}
