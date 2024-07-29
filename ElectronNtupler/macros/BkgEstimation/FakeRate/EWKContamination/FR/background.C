#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void background() {

  TString workdir;
  std::vector<TFile*> InputFiles_bkg;
  
  const char *bkg[8] = {"GammaJets", "WJetsToLNu", "TTbar", "diBoson_WW", "diBoson_WZ", "diBoson_ZZ", "Single_antiTop", "SingleTop"};
  
  //double xsec[8] = {365896.,61526.7,831.76,118.7,66.1,15.4,35.6,35.6};
  double xsec[8] = {365896.,61526.7,831.76,118.7,47.13,16.523,35.6,35.6};
  double noEvts[8] = {54136.012704,3731926637458.121094,85849690.,988418.,1000000.,985600.,999400.,1000000.};
  
  workdir = "/tmp/rchawla/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/Backgrounds/";
  
  InputFiles_bkg.clear();

  InputFiles_bkg.push_back(TFile::Open(workdir+"GammaJets_15_6000.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"WJetsToLNu.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"TTbar.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WW.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_ZZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"Single_antiTop.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"SingleTop.root"));

  int nsample = InputFiles_bkg.size();
  TFile* file[8];

  for(unsigned int jentry = 0; jentry < nsample; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_bkg.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("../../../../../dataPUDist.root");
    TFile *f2 = TFile::Open("../../../../../PileUp_MC.root");

    // data histogram 
    TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
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
    Int_t           singlePhoton;
    Double_t        theWeight;
    Int_t           nPUTrue;

    expectedMissingInnerHits = 0;
    metPt = 0;
    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    passMediumId = 0; 

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
    T1->SetBranchStatus("singlePhoton", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("nPUTrue", 1);

    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits);
    T1->SetBranchAddress("metPt", &metPt);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("singlePhoton", &singlePhoton);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("ROOTFiles/%s.root",bkg[jentry]),"RECREATE");
    
    double lumiWeight = xsec[jentry]/noEvts[jentry];
    cout<<"Background Sample: "<<bkg[jentry]<<endl;

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

    int nentries = T1->GetEntries();
    //int nentries = 500;
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

      // PU Weight
      int bin = 0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      PUWeight = weights->GetBinContent(bin);

      if(!Ele23_WPLoose) continue;
      //if(singlePhoton) continue;

      //if(metPt->at(0) >= 20.) continue;
      
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

	denPt->Fill(newelePt.at(l),lumiWeight*theWeight*PUWeight*2316.969);

	if(fabs(newscEta.at(l)) < 1.4442) denPt_BRL->Fill(newelePt.at(l),lumiWeight*theWeight*PUWeight*2316.969);
	else if(fabs(newscEta.at(l)) > 1.566 && fabs(newscEta.at(l)) < 2.5) denPt_ECAP->Fill(newelePt.at(l),lumiWeight*theWeight*PUWeight*2316.969);

	passId = neweleMediumId.at(l);
	if(passId){

	  numPt->Fill(newelePt.at(l),lumiWeight*theWeight*PUWeight*2316.969);

	  if(fabs(newscEta.at(l)) < 1.4442) numPt_BRL->Fill(newelePt.at(l),lumiWeight*theWeight*PUWeight*2316.969);
	  else if(fabs(newscEta.at(l)) > 1.566 && fabs(newscEta.at(l)) < 2.5) numPt_ECAP->Fill(newelePt.at(l),lumiWeight*theWeight*PUWeight*2316.969);

	} // ID
      }// elePt size
    } // event

    file[jentry]->Write();
    file[jentry]->Close();
    
    cout<<""<<endl;

  } // file Loop
}
