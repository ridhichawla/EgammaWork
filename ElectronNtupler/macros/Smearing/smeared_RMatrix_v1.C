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

//bool unfold_RMatrix::mySortFnx(Double_t i, Double_t j) { return i>j; }

void smeared_RMatrix_v1() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  int mass[14] = {10,10,50,50,100,200,400,500,700,800,1000,1500,2000,3000};

  double xsec[13] = {18609.9/3,18609.9/3,5789./3,5789./3,226./3,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
  //double sumofWts[13] = {771365896995.030029,771365896995.030029,144503601061.654709,144503601061.654709,219884886.214768,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};
  double sumofWts[13] = {771365896802.749878,771365896802.749878,144503601049.458435,144503601049.458435,219884886.214768,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

  workdir = "/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/DY_76X_SmearSyst/";
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
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"DY_2000to3000.root"));

  int nsample = InputFiles_signal_DY.size();
  TFile* file[13];

  for(unsigned int jentry = 0; jentry < 2; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("ntupler/ElectronTree");

    TFile *f1 = TFile::Open("../dataPUDist.root");
    TFile *f2 = TFile::Open("../PileUp_MC.root");

    // data histogram 
    TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1F *MC_puDist = (TH1F*)f2->Get("pileup_MC");
    TH1F *weights = (TH1F*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    vector<float>   *genPostFSR_Pt;
    vector<float>   *genPostFSR_Eta;
    vector<float>   *genPostFSR_Phi;
    vector<float>   *genPostFSR_En;
    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<int>     *passMediumId;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *etaSC;
    Int_t           tauFlag;
    Double_t        theWeight;
    Bool_t          Ele23_WPLoose;
    Int_t           nPUTrue;
    vector<double>  *eta_Ele23;
    vector<double>  *phi_Ele23;
    vector<double>  *error_scale;
    vector<float>  *sigma;
    vector<float>  *sigma1;
    vector<float>  *sigma2;
    vector<float>  *sigma3;
    vector<float>  *sigma4;
    double EvtNo;

    genPostFSR_Pt = 0;
    genPostFSR_Eta = 0;
    genPostFSR_Phi = 0;
    genPostFSR_En = 0;
    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;
    ptElec = 0;
    etaElec = 0;
    phiElec = 0;
    energyElec = 0;
    etaSC = 0;
    passMediumId = 0;
    eta_Ele23 = 0;
    phi_Ele23 = 0;
    error_scale = 0;
    sigma = 0;
    sigma1 = 0;
    sigma2 = 0;
    sigma3 = 0;
    sigma4 = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("genPostFSR_Pt", 1);
    T1->SetBranchStatus("genPostFSR_Eta", 1);
    T1->SetBranchStatus("genPostFSR_Phi", 1);
    T1->SetBranchStatus("genPostFSR_En", 1);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("nPUTrue", 1);
    T1->SetBranchStatus("eta_Ele23", 1);
    T1->SetBranchStatus("phi_Ele23", 1);
    T1->SetBranchStatus("EvtNo", 1);
    T1->SetBranchStatus("error_scale", 1);
    T1->SetBranchStatus("sigma", 1);
    T1->SetBranchStatus("sigma1", 1);
    T1->SetBranchStatus("sigma2", 1);
    T1->SetBranchStatus("sigma3", 1);
    T1->SetBranchStatus("sigma4", 1);

    T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
    T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
    T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
    T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);
    T1->SetBranchAddress("eta_Ele23", &eta_Ele23);
    T1->SetBranchAddress("phi_Ele23", &phi_Ele23);
    T1->SetBranchAddress("EvtNo", &EvtNo);
    T1->SetBranchAddress("error_scale", &error_scale);
    T1->SetBranchAddress("sigma", &sigma);
    T1->SetBranchAddress("sigma1", &sigma1);
    T1->SetBranchAddress("sigma2", &sigma2);
    T1->SetBranchAddress("sigma3", &sigma3);
    T1->SetBranchAddress("sigma4", &sigma4);

    file[jentry] = new TFile(Form("files_29Jan/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    int n = 50;
    int nSmearMap = 100;

    int mediumId;
    double massGPre;
    double dR1, dR2, dR3, dR4, dR5, dR6, dR7, dR8, dR9, dR10, dR11, dR12;
    double events;
    TLorentzVector ele1,ele2,dielectron,ele11,ele12,dielectron1,ele21,ele22,dielectron2,ele31,ele32,dielectron3,ele41,ele42,dielectron4,ele51,ele52,dielectron5;
    TLorentzVector gen1,gen2,digen;
    TLorentzVector gen1_NoAcc,gen2_NoAcc,digen_NoAcc;
    TLorentzVector genPre1,genPre2,diGen;

    // Branch variable declaration
    double massReco, massReco1, massReco2, massReco3, massReco4, massReco5, massGen;
    bool isReco, isReco1, isReco2, isReco3, isReco4, isReco5, isGen;
    double lumiWeight, genWeight, PUWeight;
    bool BB, BE, EE, BB1, BE1, EE1, BB2, BE2, EE2, BB3, BE3, EE3, BB4, BE4, EE4, BB5, BE5, EE5;

    vector <double> newelePt; vector <double> newelePt1; vector <double> newelePt2; vector <double> newelePt3; vector <double> newelePt4;
    vector <double> newelePt5; vector <double> neweleEta; vector <double> neweleEnr; vector <double> neweleEnr1; vector <double> neweleEnr2;
    vector <double> neweleEnr3; vector <double> neweleEnr4; vector <double> neweleEnr5;
    vector <double> newelePhi; vector <double> newscEta;
    vector <double> newgenPt; vector <double> newgenEta; vector <double> newgenEnr; vector <double> newgenPhi;
    vector <double> scale; vector <double> sig1; vector <double> sig2; vector <double> sig3; vector <double> sig4; vector <double> sig5;

    int nentries = T1->GetEntries();
    //int nentries = 20;
    cout<<"entries: "<<nentries<<endl;

    TTree* tree[nSmearMap];
    char name[nSmearMap];

    for(Int_t i_map=0; i_map<nSmearMap; i_map++) {

      sprintf(name, "%s%d", "tree", i_map);
      tree[i_map] = new TTree(name, "after preselections tree");

      tree[i_map]->Branch("massReco", &massReco, "massReco/D");
      tree[i_map]->Branch("massReco1", &massReco1, "massReco1/D");
      tree[i_map]->Branch("massReco2", &massReco2, "massReco2/D");
      tree[i_map]->Branch("massReco3", &massReco3, "massReco3/D");
      tree[i_map]->Branch("massReco5", &massReco5, "massReco5/D");
      tree[i_map]->Branch("massGen", &massGen, "massGen/D");
      tree[i_map]->Branch("isReco", &isReco, "isReco/B");
      tree[i_map]->Branch("isReco1", &isReco1, "isReco1/B");
      tree[i_map]->Branch("isReco2", &isReco2, "isReco2/B");
      tree[i_map]->Branch("isReco3", &isReco3, "isReco3/B");
      tree[i_map]->Branch("isReco5", &isReco5, "isReco5/B");
      tree[i_map]->Branch("isGen", &isGen, "isGen/B");
      tree[i_map]->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
      tree[i_map]->Branch("genWeight", &genWeight, "genWeight/D");
      tree[i_map]->Branch("PUWeight", &PUWeight, "PUWeight/D");
      tree[i_map]->Branch("BB", &BB, "BB/B"); tree[i_map]->Branch("BE", &BE, "BE/B"); tree[i_map]->Branch("EE", &EE, "EE/B");
      tree[i_map]->Branch("BB1", &BB1, "BB1/B"); tree[i_map]->Branch("BE1", &BE1, "BE1/B"); tree[i_map]->Branch("EE1", &EE1, "EE1/B");
      tree[i_map]->Branch("BB2", &BB2, "BB2/B"); tree[i_map]->Branch("BE2", &BE2, "BE2/B"); tree[i_map]->Branch("EE2", &EE2, "EE2/B");
      tree[i_map]->Branch("BB3", &BB3, "BB3/B"); tree[i_map]->Branch("BE3", &BE3, "BE3/B"); tree[i_map]->Branch("EE3", &EE3, "EE3/B");
      tree[i_map]->Branch("BB5", &BB5, "BB5/B"); tree[i_map]->Branch("BE5", &BE5, "BE5/B"); tree[i_map]->Branch("EE5", &EE5, "EE5/B");
    }

    for(int ismr=0;ismr<nSmearMap;ismr++){

      cout<<"Smearing Map: "<<ismr<<endl;

      for (unsigned int i=0; i < nentries; i++) {
	T1->GetEntry(i);

	if(i%1000000 == 0){
	  //if(i%10 == 0)
	  cout << "Events Processed :  " << i << endl;
	}

	// Sorting Post Gen level
	int idx[genPostFSR_Pt->size()];
	float ptn[genPostFSR_Pt->size()];

	for(unsigned int gn=0; gn<genPostFSR_Pt->size(); gn++)
	{
	  ptn[gn]=genPostFSR_Pt->at(gn);
	}

	int sizen = sizeof(ptn)/sizeof(ptn[0]);
	TMath::Sort(sizen,ptn,idx,true);

	massGPre=0.; massGen=-999.; massReco=-999.; massReco1=-999.; massReco2=-999.; massReco3=-999.; massReco4=-999.; massReco5=-999.;
	dR1 = 0.; dR2 = 0.; dR3 = 0.; dR4 = 0.; dR5 = 0.; dR6 = 0.; dR7 = 0.; dR8 = 0.; dR9 = 0.; dR10 = 0.; dR11 = 0.; dR12 = 0.;
	isReco = false; isReco1 = false; isReco2 = false; isReco3 = false; isReco4 = false; isReco5 = false; isGen = false;
	BB=false; BE=false; EE=false; BB1=false; BE1=false; EE1=false; BB2=false; BE2=false; EE2=false; BB3=false; BE3=false; EE3=false;
	BB4=false; BE4=false; EE4=false; BB5=false; BE5=false; EE5=false;
	mediumId = 0; 
	newelePt.clear(); newelePt1.clear(); newelePt2.clear(); newelePt3.clear(); newelePt4.clear(); newelePt5.clear(); neweleEta.clear();
	neweleEnr.clear(); neweleEnr1.clear(); neweleEnr2.clear(); neweleEnr3.clear(); neweleEnr4.clear(); neweleEnr5.clear();
	newelePhi.clear(); newscEta.clear();
	newgenPt.clear(); newgenEta.clear(); newgenEnr.clear(); newgenPhi.clear();
	scale.clear(); sig1.clear(); sig2.clear(); sig3.clear(); sig4.clear(); sig5.clear();

	// PU Weight
	int bin = 0;
	bin = weights->GetXaxis()->FindBin(nPUTrue);
	float puweight = weights->GetBinContent(bin);

	if(genPreFSR_Pt->size() == 2){
	  genPre1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	  genPre2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	  diGen=genPre1+genPre2;
	  massGPre=diGen.M();
	}

	//printf ("Event  Number = %f \n",EvtNo);

	if((jentry==2 || jentry==3) && massGPre >= 100.) continue;        // Gen Mass cut ----- for 50 to inf sample

	if(tauFlag) continue;
	if(genPreFSR_Pt->size() != 2) continue;

	// -------------------------- Gen Information with acceptance ------------------------- //
	for(unsigned int j=0;j<genPostFSR_Pt->size();j++){

	  if(fabs(genPostFSR_Eta->at(idx[j])) < 2.5 && !(fabs(genPostFSR_Eta->at(idx[j])) > 1.4442 && fabs(genPostFSR_Eta->at(idx[j])) < 1.566)){

	    newgenPt.push_back(genPostFSR_Pt->at(idx[j]));
	    newgenEta.push_back(genPostFSR_Eta->at(idx[j]));
	    newgenEnr.push_back(genPostFSR_En->at(idx[j]));
	    newgenPhi.push_back(genPostFSR_Phi->at(idx[j]));
	  } // eta on gen
	} // post FSR size

	if(newgenPt.size() == 2){
	  if(newgenPt.at(0) > 30. && newgenPt.at(1) > 10.){

	    isGen = true;
	    gen1.SetPtEtaPhiE(newgenPt.at(0),newgenEta.at(0),newgenPhi.at(0),newgenEnr.at(0));
	    gen2.SetPtEtaPhiE(newgenPt.at(1),newgenEta.at(1),newgenPhi.at(1),newgenEnr.at(1));

	    digen=gen1+gen2;
	    massGen=digen.M();
	    //cout<<"mass gen: "<<massGen<<endl;
	  } // acceptance on gen
	} // newgen size == 2

	TRandom3 eran;
	eran.SetSeed(0);

	//TRandom3* eran;
	//eran = new TRandom3(0);

	if(Ele23_WPLoose){

	  for(int k=0;k<ptElec->size();k++){

	    mediumId = passMediumId->at(k);
	    if(mediumId) {

	      if(fabs(etaSC->at(k)) < 2.5 && !(fabs(etaSC->at(k)) > 1.4442 && fabs(etaSC->at(k)) < 1.566)){

		//printf ("Event = %f \n",EvtNo);
		//cout<<"electron = "<<k<<" eta = "<<etaElec->at(k)<<" eta SC = "<<etaSC->at(k)<<endl;
		//cout<<"error_scale = "<<error_scale->at(k)<<endl;
		//cout<<"Nominal sigma = "<<sigma->at(k)<<"   "<<"sigma2 = "<<sigma1->at(k)<<"   "<<"sigma3 = "<<sigma2->at(k)<<"   "<<"sigma4 = "<<sigma3->at(k)<<"   "<<"sigma5 = "<<sigma4->at(k)<<endl;

		double shift1  = eran.Gaus(1 + error_scale->at(k), sigma->at(k));
		double shift2  = eran.Gaus(1 + error_scale->at(k), sigma1->at(k));
		double shift3  = eran.Gaus(1 + error_scale->at(k), sigma2->at(k));
		double shift5  = eran.Gaus(1 + error_scale->at(k), sigma4->at(k));

		//cout<<shift1<<"   "<<shift2<<"   "<<shift3<<"   "<<shift4<<"   "<<shift5<<endl;
		//cout<<"Unsmeared energy = "<<energyElec->at(k)<<"   Unsmeared PT = "<<ptElec->at(k)<<endl;
		//cout<<"Smeared energies = "<<energyElec->at(k)*shift1<<"   "<<energyElec->at(k)*shift2<<"   "<<energyElec->at(k)*shift3<<"   "<<energyElec->at(k)*shift4<<"   "<<energyElec->at(k)*shift5<<endl;
		//cout<<"Smeared PT       = "<<ptElec->at(k)*shift1<<"   "<<ptElec->at(k)*shift2<<"   "<<ptElec->at(k)*shift3<<"   "<<ptElec->at(k)*shift4<<"   "<<ptElec->at(k)*shift5<<endl;
		//cout<<""<<endl;

		newelePt.push_back(ptElec->at(k));
		newelePt1.push_back(ptElec->at(k)*shift1);
		newelePt2.push_back(ptElec->at(k)*shift2);
		newelePt3.push_back(ptElec->at(k)*shift3);
		newelePt5.push_back(ptElec->at(k)*shift5);

		neweleEta.push_back(etaElec->at(k));
		newscEta.push_back(etaSC->at(k));
		newelePhi.push_back(phiElec->at(k));

		neweleEnr.push_back(energyElec->at(k));
		neweleEnr1.push_back((energyElec->at(k))*shift1);
		neweleEnr2.push_back((energyElec->at(k))*shift2);
		neweleEnr3.push_back((energyElec->at(k))*shift3);
		neweleEnr5.push_back((energyElec->at(k))*shift5);

	      } // eta
	    } // ID
	  } // pt size

	  // Sorting Reco level
	  int index[newelePt.size()];
	  float pt[newelePt.size()];

	  for(unsigned int el=0; el<newelePt.size(); el++)
	  {
	    pt[el]=newelePt.at(el);
	  }

	  int size = sizeof(pt)/sizeof(pt[0]);
	  TMath::Sort(size,pt,index,true);

	  // Sorting Reco1 level
	  int index1[newelePt1.size()];
	  float pt1[newelePt1.size()];

	  for(unsigned int el1=0; el1<newelePt1.size(); el1++)
	  { 
	    pt1[el1]=newelePt1.at(el1);
	  }

	  int size1 = sizeof(pt1)/sizeof(pt1[0]);
	  TMath::Sort(size1,pt1,index1,true);

	  // Sorting Reco2 level
	  int index2[newelePt2.size()];
	  float pt2[newelePt2.size()];

	  for(unsigned int el2=0; el2<newelePt2.size(); el2++)
	  { 
	    pt2[el2]=newelePt2.at(el2);
	  }

	  int size2 = sizeof(pt2)/sizeof(pt2[0]);
	  TMath::Sort(size2,pt2,index2,true);

	  // Sorting Reco3 level
	  int index3[newelePt3.size()];
	  float pt3[newelePt3.size()];

	  for(unsigned int el3=0; el3<newelePt3.size(); el3++)
	  { 
	    pt3[el3]=newelePt3.at(el3);
	  }

	  int size3 = sizeof(pt3)/sizeof(pt3[0]);
	  TMath::Sort(size3,pt3,index3,true);

	  // Sorting Reco5 level
	  int index5[newelePt5.size()];
	  float pt5[newelePt5.size()];

	  for(unsigned int el5=0; el5<newelePt5.size(); el5++)
	  { 
	    pt5[el5]=newelePt5.at(el5);
	  }

	  int size5 = sizeof(pt5)/sizeof(pt5[0]);
	  TMath::Sort(size5,pt5,index5,true);

	  // Select two highest Pt electrons // --------------------------------------------
	  if(newelePt.size() == 2){
	    int countTrig = 0;
	    for(unsigned int l = 0; l < eta_Ele23->size(); l++){
	      dR1 = deltaR(neweleEta.at(index[0]), newelePhi.at(index[0]), eta_Ele23->at(l), phi_Ele23->at(l));
	      dR2 = deltaR(neweleEta.at(index[1]), newelePhi.at(index[1]), eta_Ele23->at(l), phi_Ele23->at(l));
	      if(dR1 < 0.1 || dR2 < 0.1) countTrig++;
	    } // eta_Ele23

	    if(newelePt.at(index[0]) > 30. && newelePt.at(index[1]) > 10. && countTrig != 0){
	      isReco = true;
	      ele1.SetPtEtaPhiE(newelePt.at(index[0]),neweleEta.at(index[0]),newelePhi.at(index[0]),neweleEnr.at(index[0]));
	      ele2.SetPtEtaPhiE(newelePt.at(index[1]),neweleEta.at(index[1]),newelePhi.at(index[1]),neweleEnr.at(index[1]));

	      dielectron=ele1+ele2;
	      massReco=dielectron.M();

	      if(fabs(newscEta.at(index[0])) < 1.4442 && fabs(newscEta.at(index[1])) < 1.4442) BB = true;
	      if((fabs(newscEta.at(index[0])) < 1.4442 && fabs(newscEta.at(index[1])) > 1.566) || (fabs(newscEta.at(index[0])) > 1.566 && fabs(newscEta.at(index[1])) < 1.4442)) BE =true;
	      if(fabs(newscEta.at(index[0])) > 1.566 && fabs(newscEta.at(index[1])) > 1.566) EE =true;
	    }
	  }

	  if(newelePt1.size() == 2){
	    int countTrig1 = 0;
	    for(unsigned int m = 0; m < eta_Ele23->size(); m++){
	      dR3 = deltaR(neweleEta.at(index1[0]), newelePhi.at(index1[0]), eta_Ele23->at(m), phi_Ele23->at(m));
	      dR4 = deltaR(neweleEta.at(index1[1]), newelePhi.at(index1[1]), eta_Ele23->at(m), phi_Ele23->at(m));
	      if(dR3 < 0.1 || dR4 < 0.1) countTrig1++;
	    }

	    if(newelePt1.at(index1[0]) > 30. && newelePt1.at(index1[1]) > 10. && countTrig1 != 0){
	      isReco1 = true;
	      ele11.SetPtEtaPhiE(newelePt1.at(index1[0]),neweleEta.at(index1[0]),newelePhi.at(index1[0]),neweleEnr1.at(index1[0]));
	      ele12.SetPtEtaPhiE(newelePt1.at(index1[1]),neweleEta.at(index1[1]),newelePhi.at(index1[1]),neweleEnr1.at(index1[1]));

	      dielectron1=ele11+ele12;
	      massReco1=dielectron1.M();

	      if(fabs(newscEta.at(index1[0])) < 1.4442 && fabs(newscEta.at(index1[1])) < 1.4442) BB1 = true;
	      if((fabs(newscEta.at(index1[0])) < 1.4442 && fabs(newscEta.at(index1[1])) > 1.566) || (fabs(newscEta.at(index1[0])) > 1.566 && fabs(newscEta.at(index1[1])) < 1.4442)) BE1 = true;
	      if(fabs(newscEta.at(index1[0])) > 1.566 && fabs(newscEta.at(index1[1])) > 1.566) EE1 = true;

	    }
	  }

	  if(newelePt2.size() == 2){
	    int countTrig2 = 0;
	    for(unsigned int n = 0; n < eta_Ele23->size(); n++){
	      dR5 = deltaR(neweleEta.at(index2[0]), newelePhi.at(index2[0]), eta_Ele23->at(n), phi_Ele23->at(n));
	      dR6 = deltaR(neweleEta.at(index2[1]), newelePhi.at(index2[1]), eta_Ele23->at(n), phi_Ele23->at(n));
	      if(dR5 < 0.1 || dR6 < 0.1) countTrig2++;
	    } // eta_Ele23 

	    if(newelePt2.at(index2[0]) > 30. && newelePt2.at(index2[1]) > 10. && countTrig2 != 0){
	      isReco2 = true;
	      ele21.SetPtEtaPhiE(newelePt2.at(index2[0]),neweleEta.at(index2[0]),newelePhi.at(index2[0]),neweleEnr2.at(index2[0]));
	      ele22.SetPtEtaPhiE(newelePt2.at(index2[1]),neweleEta.at(index2[1]),newelePhi.at(index2[1]),neweleEnr2.at(index2[1]));

	      dielectron2=ele21+ele22;
	      massReco2=dielectron2.M();

	      if(fabs(newscEta.at(index2[0])) < 1.4442 && fabs(newscEta.at(index2[1])) < 1.4442) BB2 = true;
	      if((fabs(newscEta.at(index2[0])) < 1.4442 && fabs(newscEta.at(index2[1])) > 1.566) || (fabs(newscEta.at(index2[0])) > 1.566 && fabs(newscEta.at(index2[1])) < 1.4442)) BE2 = true;
	      if(fabs(newscEta.at(index2[0])) > 1.566 && fabs(newscEta.at(index2[1])) > 1.566) EE2 = true;

	    } 
	  }

	  if(newelePt3.size() == 2){
	    int countTrig3 = 0;
	    for(unsigned int p = 0; p < eta_Ele23->size(); p++){
	      dR7 = deltaR(neweleEta.at(index3[0]), newelePhi.at(index3[0]), eta_Ele23->at(p), phi_Ele23->at(p));
	      dR8 = deltaR(neweleEta.at(index3[1]), newelePhi.at(index3[1]), eta_Ele23->at(p), phi_Ele23->at(p));
	      if(dR7 < 0.1 || dR8 < 0.1) countTrig3++;
	    } // eta_Ele23

	    if(newelePt3.at(index3[0]) > 30. && newelePt3.at(index3[1]) > 10. && countTrig3 != 0){
	      isReco3 = true;
	      ele31.SetPtEtaPhiE(newelePt3.at(index3[0]),neweleEta.at(index3[0]),newelePhi.at(index3[0]),neweleEnr3.at(index3[0]));
	      ele32.SetPtEtaPhiE(newelePt3.at(index3[1]),neweleEta.at(index3[1]),newelePhi.at(index3[1]),neweleEnr3.at(index3[1]));

	      dielectron3=ele31+ele32;
	      massReco3=dielectron3.M();

	      if(fabs(newscEta.at(index3[0])) < 1.4442 && fabs(newscEta.at(index3[1])) < 1.4442) BB3 = true;
	      if((fabs(newscEta.at(index3[0])) < 1.4442 && fabs(newscEta.at(index3[1])) > 1.566) || (fabs(newscEta.at(index3[0])) > 1.566 && fabs(newscEta.at(index3[1])) < 1.4442)) BE3 = true;
	      if(fabs(newscEta.at(index3[0])) > 1.566 && fabs(newscEta.at(index3[1])) > 1.566) EE3 = true;

	    }
	  }

	  if(newelePt5.size() == 2){
	    int countTrig5 = 0;
	    for(unsigned int r = 0; r < eta_Ele23->size(); r++){
	      dR11 = deltaR(neweleEta.at(index5[0]), newelePhi.at(index5[0]), eta_Ele23->at(r), phi_Ele23->at(r));
	      dR12 = deltaR(neweleEta.at(index5[1]), newelePhi.at(index5[1]), eta_Ele23->at(r), phi_Ele23->at(r));
	      if(dR11 < 0.1 || dR12 < 0.1) countTrig5++;
	    } // eta_Ele23

	    if(newelePt5.at(index5[0]) > 30. && newelePt5.at(index5[1]) > 10. && countTrig5 != 0){
	      isReco5 = true;
	      ele51.SetPtEtaPhiE(newelePt5.at(index5[0]),neweleEta.at(index5[0]),newelePhi.at(index5[0]),neweleEnr5.at(index5[0]));
	      ele52.SetPtEtaPhiE(newelePt5.at(index5[1]),neweleEta.at(index5[1]),newelePhi.at(index5[1]),neweleEnr5.at(index5[1]));

	      dielectron5=ele51+ele52;
	      massReco5=dielectron5.M();

	      if(fabs(newscEta.at(index5[0])) < 1.4442 && fabs(newscEta.at(index5[1])) < 1.4442) BB5 = true;
	      if((fabs(newscEta.at(index5[0])) < 1.4442 && fabs(newscEta.at(index5[1])) > 1.566) || (fabs(newscEta.at(index5[0])) > 1.566 && fabs(newscEta.at(index5[1])) < 1.4442)) BE5 = true;
	      if(fabs(newscEta.at(index5[0])) > 1.566 && fabs(newscEta.at(index5[1])) > 1.566) EE5 = true;

	    }
	  }

	} // trigger fire

	lumiWeight = lumi_Weight;
	genWeight = theWeight;
	PUWeight = puweight;

	//cout<<"mass reco: "<<massReco<<endl;
	tree[ismr]->Fill();
      } // event

      cout<<""<<endl;
    } // smeared maps*/

    file[jentry]->Write();
    file[jentry]->Close();

    //cout<<""<<endl;
  } //file Loop
}
