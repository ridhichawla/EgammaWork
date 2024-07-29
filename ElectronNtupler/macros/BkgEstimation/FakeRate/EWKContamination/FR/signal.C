#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

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

    TFile *f1 = TFile::Open("../../../../../dataPUDist.root");
    TFile *f2 = TFile::Open("../../../../../PileUp_MC.root");

    TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    TH1D *MC_puDist = (TH1D*)f2->Get("pileup_MC");
    TH1D *weights = (TH1D*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    Bool_t          Ele23_WPLoose;
    vector<int>     *expectedMissingInnerHits;
    vector<double>  *metPt;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    vector<int>     *passMediumId;
    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Rap;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    Int_t           tauFlag;
    Int_t           singlePhoton;
    Double_t        theWeight;
    Int_t           nPV;
    Int_t           nPUTrue;

    expectedMissingInnerHits = 0;
    metPt = 0;
    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    passMediumId = 0; 
    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Rap = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("expectedMissingInnerHits", 1);
    T1->SetBranchStatus("metPt", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Rap", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("singlePhoton", 1);
    T1->SetBranchStatus("nPV", 1);
    T1->SetBranchStatus("nPUTrue", 1);
    T1->SetBranchStatus("theWeight", 1);

    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits);
    T1->SetBranchAddress("metPt", &metPt);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("singlePhoton", &singlePhoton);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("nPV", &nPV);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("ROOTFiles/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");

    int count;
    bool passKin;
    int mediumId, passId;
    double massGen;
    TLorentzVector gen1,gen2,diGen;

    vector <double> newelePt; vector <double> neweleEta;  vector <double> newscEta; vector <double> neweleMediumId;

    Double_t x1bin[10] = {10,15,20,25,30,40,50,70,100,10000};
    int nbins = 9;
    
    TH1D *numPt      = new TH1D("numPt", "numPt", nbins, x1bin);
    TH1D *numPt_BRL  = new TH1D("numPt_BRL", "numPt_BRL", nbins, x1bin);
    TH1D *numPt_ECAP = new TH1D("numPt_ECAP", "numPt_ECAP", nbins, x1bin);

    TH1D *denPt      = new TH1D("denPt", "denPt", nbins, x1bin);
    TH1D *denPt_BRL  = new TH1D("denPt_BRL", "denPt_BRL", nbins, x1bin);
    TH1D *denPt_ECAP = new TH1D("denPt_ECAP", "denPt_ECAP", nbins, x1bin);

    numPt->Sumw2(); denPt->Sumw2();
    numPt_BRL->Sumw2(); denPt_BRL->Sumw2();
    numPt_ECAP->Sumw2(); denPt_ECAP->Sumw2();

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    int nentries = T1->GetEntries();
    //int nentries = 50;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      // Sorting
      int index[ptElec->size()];
      float pt[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++) {
	pt[el]=ptElec->at(el); }

      int size = sizeof(pt)/sizeof(pt[0]);
      TMath::Sort(size,pt,index,true);

      count = 0;
      passKin = false;
      mediumId = 0; passId = 0;
      massGen = 0.0;
      double PUWeight = 1.0;

      newelePt.clear(); neweleEta.clear(); newscEta.clear(); neweleMediumId.clear();

      if(genPreFSR_Pt->size() == 2){
	gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=gen1+gen2;
	massGen=diGen.M();
      }

      if(!Ele23_WPLoose) continue;
      //if(!singlePhoton) continue;
      
      //if(metPt->at(0) >= 20.) continue;
      
      if(jentry==1 && massGen >= 100.) continue;
      if(tauFlag) continue;

      int bin = 0;
      double puWeights = 1.0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      PUWeight = weights->GetBinContent(bin);

      for(int j=0;j<ptElec->size();j++){

	mediumId = passMediumId->at(index[j]);
	passKin = (fabs(etaSC->at(index[j])) < 2.5 && !(fabs(etaSC->at(index[j])) > 1.4442 && fabs(etaSC->at(index[j])) < 1.566));

	if(mediumId && passKin) count++;
      }

      if(count <= 1){
	for(unsigned int k=0;k<ptElec->size();k++){
	  if(fabs(etaSC->at(index[k])) < 2.5){

	    if(etaSC->at(index[k]) > 2.5) cout<<"eta > 2.5   "<<"entry: "<<jentry<<"   "<<"no. of electrons: "<<ptElec->size()<<"   "<<etaSC->at(index[k])<<endl;

	    newelePt.push_back(ptElec->at(index[k]));
	    neweleEta.push_back(etaElec->at(index[k]));
	    newscEta.push_back(etaSC->at(index[k]));
	    neweleMediumId.push_back(passMediumId->at(index[k]));
	  }
	}
      }

      for(unsigned int l=0;l<newelePt.size();l++){

	denPt->Fill(newelePt.at(l),lumi_Weight*PUWeight*2316.969*theWeight);

	if(fabs(newscEta.at(l)) < 1.4442) denPt_BRL->Fill(newelePt.at(l),lumi_Weight*PUWeight*2316.969*theWeight);
	else if(fabs(newscEta.at(l)) > 1.566 && fabs(newscEta.at(l)) < 2.5) denPt_ECAP->Fill(newelePt.at(l),lumi_Weight*PUWeight*2316.969*theWeight);

	passId = neweleMediumId.at(l);
	if(passId){

	  numPt->Fill(newelePt.at(l),lumi_Weight*PUWeight*2316.969*theWeight);

	  if(fabs(newscEta.at(l)) < 1.4442) numPt_BRL->Fill(newelePt.at(l),lumi_Weight*PUWeight*2316.969*theWeight);
	  else if(fabs(newscEta.at(l)) > 1.566 && fabs(newscEta.at(l)) < 2.5) numPt_ECAP->Fill(newelePt.at(l),lumi_Weight*PUWeight*2316.969*theWeight);

	} // ID
      }// elePt size 
    } // event

    file[jentry]->Write();
    file[jentry]->Close();
    
    cout<<""<<endl;

  } // file Loop
}
