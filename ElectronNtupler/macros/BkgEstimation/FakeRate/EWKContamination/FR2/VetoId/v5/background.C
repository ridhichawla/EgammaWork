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
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

void background() {

  TString workdir;
  std::vector<TFile*> InputFiles_bkg;

  const char *bkg[9] = {"GammaJets", "WJetsToLNu", "TTbar", "diBoson_WW", "diBoson_WZ", "diBoson_ZZ", "Single_antiTop", "SingleTop", "WGamma"};

  //double xsec[8] = {365896.,61526.7,831.76,118.7,66.1,15.4,35.6,35.6};
  double xsec[9] = {365896.,61526.7,831.76,118.7,47.13,16.523,35.6,35.6,489.};
  double noEvts[9] = {52839.384531,3731926637458.121094,97994304.,988416.,999996.,985598.,999399.,999999.,4155331783.583204};

  workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated_19042017/";

  InputFiles_bkg.clear();

  InputFiles_bkg.push_back(TFile::Open(workdir+"GammaJets_15_6000.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"WJetsToLNu.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"TTbar.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WW.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_WZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"diBoson_ZZ.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"Single_antiTop.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"SingleTop.root"));
  InputFiles_bkg.push_back(TFile::Open(workdir+"WGamma.root"));

  int nsample = InputFiles_bkg.size();
  TFile* file[9];

  for(unsigned int jentry = 0; jentry < nsample; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_bkg.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("/afs/cern.ch/user/r/rchawla/dataPUDist.root");
    TFile *f2 = TFile::Open("/afs/cern.ch/user/r/rchawla/PileUp_MC.root");

    // data histogram 
    TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1D *MC_puDist = (TH1D*)f2->Get("pileup_MC");
    TH1D *weights = (TH1D*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    Bool_t          Ele23_WPLoose;
    vector<double>  *etSPhoHLT;
    vector<double>  *etaSPhoHLT;
    vector<double>  *phiSPhoHLT;
    Int_t           singlePhoton;
    Int_t           prescalePhoton;
    vector<int>     *expectedMissingInnerHits;
    vector<double>  *metPt;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    vector<int>     *passVetoId;
    vector<int>     *passMediumId;
    Double_t        theWeight;
    Int_t           nPUTrue;

    etSPhoHLT = 0;
    etaSPhoHLT = 0;
    phiSPhoHLT = 0;
    expectedMissingInnerHits = 0;
    metPt = 0;
    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    passVetoId = 0; 
    passMediumId = 0; 

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("etSPhoHLT", 1);
    T1->SetBranchStatus("etaSPhoHLT", 1);
    T1->SetBranchStatus("phiSPhoHLT", 1);
    T1->SetBranchStatus("singlePhoton", 1);
    T1->SetBranchStatus("prescalePhoton", 1);
    T1->SetBranchStatus("expectedMissingInnerHits", 1);
    T1->SetBranchStatus("metPt", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passVetoId", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("nPUTrue", 1);

    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("etSPhoHLT", &etSPhoHLT);
    T1->SetBranchAddress("etaSPhoHLT", &etaSPhoHLT);
    T1->SetBranchAddress("phiSPhoHLT", &phiSPhoHLT);
    T1->SetBranchAddress("singlePhoton", &singlePhoton);
    T1->SetBranchAddress("prescalePhoton", &prescalePhoton);
    T1->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits);
    T1->SetBranchAddress("metPt", &metPt);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passVetoId", &passVetoId);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("%s.root",bkg[jentry]),"RECREATE");

    double lumiWeight = xsec[jentry]/noEvts[jentry];
    cout<<"Background Sample: "<<bkg[jentry]<<endl;

    int count;
    bool passKin;
    int mediumId, vetoId, passId;
    double dR; vector <bool> isJet; vector <double> idx;
    double massGen;
    TLorentzVector gen1,gen2,diGen;

    vector <double> newelePt; vector <double> neweleEta; vector <double> newelePhi; vector <double> newscEta; vector <double> neweleMediumId;
    vector <double> elePt; vector <double> eleEta; vector <double> scEta; vector <double> eleMediumId;

    Double_t x1bin[6] = {10,20,30,40,50,10000};
    //Double_t x1bin[10] = {10,15,20,25,30,40,50,70,100,10000};
    int nbins = 5;

    TH1D *numPt      = new TH1D("numPt", "numPt", nbins, x1bin);
    TH1D *numPt_BRL  = new TH1D("numPt_BRL", "numPt_BRL", nbins, x1bin);
    TH1D *numPt_ECAP = new TH1D("numPt_ECAP", "numPt_ECAP", nbins, x1bin);

    TH1D *denPt      = new TH1D("denPt", "denPt", nbins, x1bin);
    TH1D *denPt_BRL  = new TH1D("denPt_BRL", "denPt_BRL", nbins, x1bin);
    TH1D *denPt_ECAP = new TH1D("denPt_ECAP", "denPt_ECAP", nbins, x1bin);

    TH1D *numEta     = new TH1D("numEta", "numEta", 60, -3, 3);
    TH1D *denEta     = new TH1D("denEta", "denEta", 60, -3, 3);

    TH1D *delta_R1    = new TH1D("delta_R1", "delta_R1", 1000, 0, 10);
    TH1D *delta_R2    = new TH1D("delta_R2", "delta_R2", 1000, 0, 10);

    numPt->Sumw2(); denPt->Sumw2(); numEta->Sumw2(); denEta->Sumw2();
    numPt_BRL->Sumw2(); denPt_BRL->Sumw2();
    numPt_ECAP->Sumw2(); denPt_ECAP->Sumw2();
    delta_R1->Sumw2(); delta_R2->Sumw2();

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
      mediumId = 0; vetoId = 0; passId = 0;
      dR = 0.;
      massGen = 0.0;
      double PUWeight = 1.0;
      newelePt.clear(); neweleEta.clear(); newelePhi.clear(); newscEta.clear(); neweleMediumId.clear();
      elePt.clear(); eleEta.clear(); scEta.clear(); eleMediumId.clear();

      // PU Weight
      int bin = 0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      PUWeight = weights->GetBinContent(bin);

      if(!singlePhoton) continue;
      //if(!Ele23_WPLoose) continue;

      if(metPt->at(0) >= 20.) continue;

      for(int j=0;j<ptElec->size();j++){

	mediumId = passMediumId->at(index[j]);

	if(mediumId) count++;
      }

      for(unsigned int k=0;k<ptElec->size();k++){

	vetoId = passVetoId->at(index[k]);
	passKin = (ptElec->at(index[k]) > 10. && fabs(etaSC->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566));

	if(vetoId && passKin){

	  if(etaSC->at(index[k]) > 2.5) cout<<"eta > 2.5   "<<"entry: "<<jentry<<"   "<<"no. of electrons: "<<ptElec->size()<<"   "<<etaSC->at(index[k])<<endl;

	  newelePt.push_back(ptElec->at(index[k]));
	  neweleEta.push_back(etaElec->at(index[k]));
	  newelePhi.push_back(phiElec->at(index[k]));
	  newscEta.push_back(etaSC->at(index[k]));
	  neweleMediumId.push_back(passMediumId->at(index[k]));
	}
      }

      isJet.clear(); idx.clear();

      for(unsigned int m=0; m<etaSPhoHLT->size(); m++){
	for(unsigned int n=0;n<newelePt.size();n++){

	  dR = deltaR(etaSPhoHLT->at(m), phiSPhoHLT->at(m), neweleEta.at(n), newelePhi.at(n));
	  delta_R1->Fill(dR);

	  if(dR > 0.05){
	    delta_R2->Fill(dR);
	    isJet.push_back(true);
	    idx.push_back(n);
	  }
	}
      }

      if(isJet.size() > 0){
	for(unsigned int s=0; s<idx.size(); s++){
	  elePt.push_back(newelePt.at(idx[s]));
	  eleEta.push_back(neweleEta.at(idx[s]));
	  scEta.push_back(newscEta.at(idx[s]));
	  eleMediumId.push_back(neweleMediumId.at(idx[s]));
	}
      }

      for(unsigned int l=0;l<elePt.size();l++){

	denPt->Fill(elePt.at(l),lumiWeight*theWeight*PUWeight*2258.066);
	denEta->Fill(eleEta.at(l),lumiWeight*theWeight*PUWeight*2258.066);

	if(fabs(scEta.at(l)) < 1.4442) denPt_BRL->Fill(elePt.at(l),lumiWeight*theWeight*PUWeight*2258.066);
	else if(fabs(scEta.at(l)) > 1.566 && fabs(scEta.at(l)) < 2.5) denPt_ECAP->Fill(elePt.at(l),lumiWeight*theWeight*PUWeight*2258.066);

	passId = eleMediumId.at(l);
	if(passId){

	  numPt->Fill(elePt.at(l),lumiWeight*theWeight*PUWeight*2258.066);
	  numEta->Fill(eleEta.at(l),lumiWeight*theWeight*PUWeight*2258.066);

	  if(fabs(scEta.at(l)) < 1.4442) numPt_BRL->Fill(elePt.at(l),lumiWeight*theWeight*PUWeight*2258.066);
	  else if(fabs(scEta.at(l)) > 1.566 && fabs(scEta.at(l)) < 2.5) numPt_ECAP->Fill(elePt.at(l),lumiWeight*theWeight*PUWeight*2258.066);

	} // ID
      }// elePt size
    } // event

    file[jentry]->Write();
    file[jentry]->Close();

    cout<<""<<endl;

  } // file Loop
}
