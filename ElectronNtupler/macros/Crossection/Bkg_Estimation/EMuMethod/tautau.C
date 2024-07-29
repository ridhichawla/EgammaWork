#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

void tautau() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;

  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};
  
  double xsec[11] = {6016.88,1873.52,76.2401,2.67606,0.139728,0.0792496,0.0123176,0.01042,0.00552772,0.000741613,0.000178737};
  //double xsec[11] = {18610./3,5870./3,226./3,7.67/3,0.423/3,0.24/3,0.035/3,0.03/3,0.016/3,0.002/3,0.00054/3};
  double sumofWts[11] = {771413889185.162476,144505031098.323120,219889705.060318,7008766.904321,122987.746342,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

  workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";
  
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

    TFile *f1 = TFile::Open("../../../dataPUDist.root");
    TFile *f2 = TFile::Open("../../../PileUp_MC.root");

    // data histogram 
    TH1D *DATA_puDist = (TH1D*)f1->Get("pileup");
    DATA_puDist->Scale(1/DATA_puDist->Integral());

    // mc histogram 
    TH1D *MC_puDist = (TH1D*)f2->Get("pileup_MC");
    TH1D *weights = (TH1D*)DATA_puDist->Clone("weights");
    weights->Divide(MC_puDist);

    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Rap;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<float>   *ptElec;
    vector<float>   *etaElec;
    vector<float>   *rapElec;
    vector<float>   *phiElec;
    vector<float>   *energyElec;
    vector<float>   *chargeElec;
    vector<float>   *etaSC;
    vector<float>   *ptMuon;
    vector<float>   *etaMuon;
    vector<float>   *phiMuon;
    vector<float>   *energyMuon;
    vector<float>   *chargeMuon;
    vector<float>   *isoPFMuon;
    vector<bool>    *isTightMuon;
    vector<int>     *passMediumId;
    Int_t           tauFlag;
    Double_t        theWeight;
    bool            Mu8_Ele17;
    bool            Ele23_WPLoose;
    Int_t           nPUTrue;

    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Rap = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;
    ptElec = 0;
    etaElec = 0;
    rapElec = 0;
    phiElec = 0;
    energyElec = 0;
    chargeElec = 0;
    etaSC = 0;
    passMediumId = 0;
    ptMuon = 0;
    etaMuon = 0;
    phiMuon = 0;
    energyMuon = 0;
    chargeMuon = 0;
    isoPFMuon = 0;
    isTightMuon = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Rap", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("ptElec", 1);
    T1->SetBranchStatus("etaElec", 1);
    T1->SetBranchStatus("phiElec", 1);
    T1->SetBranchStatus("energyElec", 1); 
    T1->SetBranchStatus("chargeElec", 1);
    T1->SetBranchStatus("etaSC", 1);
    T1->SetBranchStatus("passMediumId", 1);
    T1->SetBranchStatus("ptMuon", 1);
    T1->SetBranchStatus("etaMuon", 1);
    T1->SetBranchStatus("phiMuon", 1); 
    T1->SetBranchStatus("energyMuon", 1);
    T1->SetBranchStatus("chargeMuon", 1);
    T1->SetBranchStatus("isoPFMuon", 1);
    T1->SetBranchStatus("isTightMuon", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("theWeight", 1);
    T1->SetBranchStatus("Mu8_Ele17", 1);
    T1->SetBranchStatus("Ele23_WPLoose", 1);
    T1->SetBranchStatus("nPUTrue", 1);

    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("ptElec", &ptElec);
    T1->SetBranchAddress("etaElec", &etaElec);
    T1->SetBranchAddress("phiElec", &phiElec);
    T1->SetBranchAddress("energyElec", &energyElec);
    T1->SetBranchAddress("chargeElec", &chargeElec);
    T1->SetBranchAddress("etaSC", &etaSC);
    T1->SetBranchAddress("passMediumId", &passMediumId);
    T1->SetBranchAddress("ptMuon", &ptMuon);
    T1->SetBranchAddress("etaMuon", &etaMuon);
    T1->SetBranchAddress("phiMuon", &phiMuon);
    T1->SetBranchAddress("energyMuon", &energyMuon);
    T1->SetBranchAddress("chargeMuon", &chargeMuon);
    T1->SetBranchAddress("isoPFMuon", &isoPFMuon);
    T1->SetBranchAddress("isTightMuon", &isTightMuon);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("theWeight", &theWeight);
    T1->SetBranchAddress("Mu8_Ele17", &Mu8_Ele17);
    T1->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);

    file[jentry] = new TFile(Form("ROOTFiles/DYTT_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    double massGen;
    int mediumId;
    TLorentzVector gen1,gen2,diGen;
    TLorentzVector ele1,ele2,diMass;
    TLorentzVector emu1,emu2,diEMuMass;
    TLorentzVector wjet1,wjet2,diWjetMass;

    vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge; vector <double> newscEta;
    vector <double> newmuonPt; vector <double> newmuonEta; vector <double> newmuonEnr; vector <double> newmuonPhi; vector <double> newmuonCharge;

    // Branch variable declaration
    double ElePT, EleEta, ElePhi, MuonPT, MuonEta, MuonPhi;
    double EEMass, EMuMass, WjetMass;
    double lumiWeight, genWeight, PUWeight;

    // Branch declaration
    tree->Branch("ElePT", &ElePT, "ElePT/D");
    tree->Branch("EleEta", &EleEta, "EleEta/D");
    tree->Branch("ElePhi", &ElePhi, "ElePhi/D");
    tree->Branch("MuonPT", &MuonPT, "MuonPT/D");
    tree->Branch("MuonEta", &MuonEta, "MuonEta/D");
    tree->Branch("MuonPhi", &MuonPhi, "MuonPhi/D");
    tree->Branch("EEMass", &EEMass, "EEMass/D");
    tree->Branch("EMuMass", &EMuMass, "EMuMass/D");
    tree->Branch("WjetMass", &WjetMass, "WjetMass/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/D");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    int nentries = T1->GetEntries();
    //int nentries = 5000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {

      ElePT = -999.; EleEta = -999.; ElePhi = -999.; MuonPT = -999.; MuonEta = -999.; MuonPhi = -999.;
      EEMass = -999.; EMuMass = -999.; WjetMass = -999.;

      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      // Sorting for Electron
      int index1[ptElec->size()];
      float pt1[ptElec->size()];

      for(unsigned int el=0; el<ptElec->size(); el++) {
	pt1[el]=ptElec->at(el); }

      int size1 = sizeof(pt1)/sizeof(pt1[0]);
      TMath::Sort(size1,pt1,index1,true);

      mediumId = 0;
      newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear(); newscEta.clear();
      newmuonPt.clear(); newmuonEta.clear(); newmuonEnr.clear(); newmuonPhi.clear(); newmuonCharge.clear();

      // Sorting for muons
      int index2[ptMuon->size()];
      float pt2[ptMuon->size()];

      for(unsigned int mu=0; mu<ptMuon->size(); mu++) {
	pt2[mu]=ptMuon->at(mu); }

      int size2 = sizeof(pt2)/sizeof(pt2[0]);
      TMath::Sort(size2,pt2,index2,true);

      // PU Weight
      int bin = 0;
      bin = weights->GetXaxis()->FindBin(nPUTrue);
      double puweight = weights->GetBinContent(bin);

      if(genPreFSR_Pt->size() > 2) cout<<"size > 2"<<endl;

      if(genPreFSR_Pt->size() == 2){
	gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));

	diGen=gen1+gen2;
	massGen=diGen.M();

      }

      if(jentry==1 && massGen >= 100.) continue;        // Gen Mass cut ----- for 50 to inf sample
      
      if(!tauFlag) continue;

      for(int j=0;j<ptElec->size();j++){

	mediumId = passMediumId->at(index1[j]);

	if(mediumId) {

	  if(fabs(etaSC->at(index1[j])) < 2.5 && !(fabs(etaSC->at(index1[j])) > 1.4442 && fabs(etaSC->at(index1[j])) < 1.566)){

	    newelePt.push_back(ptElec->at(index1[j]));
	    neweleEta.push_back(etaElec->at(index1[j]));
	    newscEta.push_back(etaSC->at(index1[j]));
	    neweleEnr.push_back(energyElec->at(index1[j]));
	    newelePhi.push_back(phiElec->at(index1[j]));
	    neweleCharge.push_back(chargeElec->at(index1[j]));

	  } // eta
	} // ID
      } // size

      for(int k=0;k<ptMuon->size();k++){

	if(fabs(etaMuon->at(index2[k])) < 2.4){
	  if(isTightMuon->at(index2[k])){
	    if(isoPFMuon->at(index2[k]) < 0.15){

	      newmuonPt.push_back(ptMuon->at(index2[k]));
	      newmuonEta.push_back(etaMuon->at(index2[k]));
	      newmuonPhi.push_back(phiMuon->at(index2[k]));
	      newmuonEnr.push_back(energyMuon->at(index2[k]));
	      newmuonCharge.push_back(chargeMuon->at(index2[k]));

	    } // Isolation
	  } // ID
	} // eta
      } // size

      if(Ele23_WPLoose){
	if(newelePt.size()==2){
	  if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){
	    if(neweleCharge.at(0)*neweleCharge.at(1) == -1){

	      ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	      ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	      diMass=ele1+ele2;
	      EEMass = diMass.M();

	    } // opposite charge
	  } // pt cut
	} // only two electrons
      }

      if(Mu8_Ele17){
	if(newelePt.size()==1 && newmuonPt.size()==1){
	  if(newelePt.at(0) > 25. && newmuonPt.at(0) > 15.){
	    if(neweleCharge.at(0)*newmuonCharge.at(0) == -1){

	      ElePT = newelePt.at(0);
	      EleEta = neweleEta.at(0);
	      ElePhi = newelePhi.at(0);

	      MuonPT = newmuonPt.at(0);
	      MuonEta = newmuonEta.at(0);
	      MuonPhi = newmuonPhi.at(0);

	      emu1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	      emu2.SetPtEtaPhiE(newmuonPt.at(0),newmuonEta.at(0),newmuonPhi.at(0),newmuonEnr.at(0));

	      diEMuMass=emu1+emu2;
	      EMuMass = diEMuMass.M();

	    } // opposite charge
	  } // pt cut
	} // one ele and one muon
      }

      if(Mu8_Ele17){
	if(newelePt.size()==1 && newmuonPt.size()==1){
	  if(newelePt.at(0) > 25. && newmuonPt.at(0) > 15.){
	    if(neweleCharge.at(0)*newmuonCharge.at(0) == 1){

	      wjet1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	      wjet2.SetPtEtaPhiE(newmuonPt.at(0),newmuonEta.at(0),newmuonPhi.at(0),newmuonEnr.at(0));

	      diWjetMass=wjet1+wjet2;
	      WjetMass = diWjetMass.M();

	    } // opposite charge
	  } // pt cut
	} // one ele and one muon
      }

      lumiWeight = lumi_Weight;
      genWeight  = theWeight;
      PUWeight = puweight;

      tree->Fill();
      // } // pre FSR == 2
    } // event

    file[jentry]->Write();
    file[jentry]->Close();
    
    cout<<""<<endl;

  } // file Loop
}
