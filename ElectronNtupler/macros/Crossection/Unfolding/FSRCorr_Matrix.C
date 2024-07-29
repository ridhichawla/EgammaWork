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

void FSRCorr_Matrix() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  
  int mass[12] = {10,50,100,200,400,500,700,800,1000,1500,2000,3000};

  double xsec[11] = {6016.88,1873.52,76.2401,2.67606,0.139728,0.0792496,0.0123176,0.01042,0.00552772,0.000741613,0.000178737}; 
  //double xsec[11] = {18610./3,5870./3,226./3,7.67/3,0.423/3,0.24/3,0.035/3,0.03/3,0.016/3,0.002/3,0.00054/3};
  double sumofWts[11] = {771413889185.162476,144505031098.323120,219889705.060318,7008766.904321,122987.746342,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

  workdir = "/eos/cms/store/group/phys_higgs/cmshww/arun/DYAnalysis_76X_Calibrated/DY_Signal/";
  //workdir = "/eos/cms/store/user/arun/DYAnalysis_76X_Calibrated_24072017/DY_Signal/";
  
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

    vector<float>   *genPostFSR_Pt;
    vector<float>   *genPostFSR_Eta;
    vector<float>   *genPostFSR_Rap;
    vector<float>   *genPostFSR_Phi;
    vector<float>   *genPostFSR_En;
    vector<float>   *genPreFSR_Pt;
    vector<float>   *genPreFSR_Eta;
    vector<float>   *genPreFSR_Rap;
    vector<float>   *genPreFSR_Phi;
    vector<float>   *genPreFSR_En;
    vector<float>   *genPhoton_Pt;
    vector<float>   *genPhoton_Eta;
    vector<float>   *genPhoton_Phi;
    vector<float>   *genPhoton_En;
    Int_t           tauFlag;
    Int_t           nPUTrue;
    Double_t        theWeight;

    genPostFSR_Pt = 0;
    genPostFSR_Eta = 0;
    genPostFSR_Rap = 0;
    genPostFSR_Phi = 0;
    genPostFSR_En = 0;
    genPreFSR_Pt = 0;
    genPreFSR_Eta = 0;
    genPreFSR_Rap = 0;
    genPreFSR_Phi = 0;
    genPreFSR_En = 0;
    genPhoton_Pt = 0;
    genPhoton_Eta = 0;
    genPhoton_Phi = 0;
    genPhoton_En = 0;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("genPostFSR_Pt", 1);
    T1->SetBranchStatus("genPostFSR_Eta", 1);
    T1->SetBranchStatus("genPostFSR_Rap", 1);
    T1->SetBranchStatus("genPostFSR_Phi", 1);
    T1->SetBranchStatus("genPostFSR_En", 1);
    T1->SetBranchStatus("genPreFSR_Pt", 1);
    T1->SetBranchStatus("genPreFSR_Eta", 1);
    T1->SetBranchStatus("genPreFSR_Rap", 1);
    T1->SetBranchStatus("genPreFSR_Phi", 1);
    T1->SetBranchStatus("genPreFSR_En", 1);
    T1->SetBranchStatus("genPhoton_Pt", 1);
    T1->SetBranchStatus("genPhoton_Eta", 1);
    T1->SetBranchStatus("genPhoton_Phi", 1);
    T1->SetBranchStatus("genPhoton_En", 1);
    T1->SetBranchStatus("tauFlag", 1);
    T1->SetBranchStatus("nPUTrue", 1);
    T1->SetBranchStatus("theWeight", 1);

    T1->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt);
    T1->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta);
    T1->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap);
    T1->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi);
    T1->SetBranchAddress("genPostFSR_En", &genPostFSR_En);
    T1->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt);
    T1->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta);
    T1->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap);
    T1->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi);
    T1->SetBranchAddress("genPreFSR_En", &genPreFSR_En);
    T1->SetBranchAddress("genPhoton_Pt", &genPhoton_Pt);
    T1->SetBranchAddress("genPhoton_Eta", &genPhoton_Eta);
    T1->SetBranchAddress("genPhoton_Phi", &genPhoton_Phi);
    T1->SetBranchAddress("genPhoton_En", &genPhoton_En);
    T1->SetBranchAddress("tauFlag", &tauFlag);
    T1->SetBranchAddress("nPUTrue", &nPUTrue);
    T1->SetBranchAddress("theWeight", &theWeight);

    file[jentry] = new TFile(Form("FSRCorr/default/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    //file[jentry] = new TFile(Form("FSRCorr/default/check_24072017_ntuples/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");

    TTree *tree = new TTree("tree"," after preselections tree");
    //TTree *tree = new TTree("tree"," after preselections tree");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    double dR, massGPre;
    double postFSR_Mass;
    double preFSR_Mass;
    double lumiWeight, genWeight;
    TLorentzVector gen_preFSR;
    TLorentzVector fourmom;
    TLorentzVector SumPhotonMom;
    TLorentzVector pre1,pre2,diPre;
    TLorentzVector post1,post2,diPost;
    TLorentzVector gen1,gen2,diGen;

    // Branch declaration
    tree->Branch("postFSR_Mass", &postFSR_Mass, "postFSR_Mass/D");
    tree->Branch("preFSR_Mass", &preFSR_Mass, "preFSR_Mass/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");

    vector <double> newgenPt; vector <double> newgenEta; vector <double> newgenPhi; vector <double> newgenEn;
    vector <double> newphoPt; vector <double> newphoEta; vector <double> newphoPhi; vector <double> newphoEn;
    vector <double> gPreFSR_Pt; vector <double> gPreFSR_Eta; vector <double> gPreFSR_Phi; vector <double> gPreFSR_En;

    int nentries = T1->GetEntries();
    //int nentries = 100000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i=0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      // Sorting Gen Photons
      int index1[genPhoton_Pt->size()];
      float pt1[genPhoton_Pt->size()];

      for(unsigned int ph=0; ph<genPhoton_Pt->size(); ph++)
      { 
	pt1[ph]=genPhoton_Pt->at(ph);
      }

      int size1 = sizeof(pt1)/sizeof(pt1[0]);
      TMath::Sort(size1,pt1,index1,true);

      // Sorting Post Gen level
      int index2[genPostFSR_Pt->size()];
      float pt2[genPostFSR_Pt->size()];

      for(unsigned int gn=0; gn<genPostFSR_Pt->size(); gn++)
      {
	pt2[gn]=genPostFSR_Pt->at(gn);
      }

      int sizen = sizeof(pt2)/sizeof(pt2[0]);
      TMath::Sort(sizen,pt2,index2,true);

      dR = 0.;
      postFSR_Mass = -999.; preFSR_Mass = -999.;
      massGPre = 0.;

      newgenPt.clear(); newgenPt.clear(); newgenEta.clear(); newgenPhi.clear(); newgenEn.clear();
      newphoPt.clear(); newgenPt.clear(); newphoEta.clear(); newphoPhi.clear(); newphoEn.clear();
      gPreFSR_Pt.clear(); gPreFSR_Eta.clear(); gPreFSR_Phi.clear(); gPreFSR_En.clear();

      if(genPreFSR_Pt->size() == 2){
	gen1.SetPtEtaPhiE(genPreFSR_Pt->at(0),genPreFSR_Eta->at(0),genPreFSR_Phi->at(0),genPreFSR_En->at(0));
	gen2.SetPtEtaPhiE(genPreFSR_Pt->at(1),genPreFSR_Eta->at(1),genPreFSR_Phi->at(1),genPreFSR_En->at(1));
	diGen=gen1+gen2;
	massGPre=diGen.M();
      }

      if(jentry==1 && massGPre >= 100.) continue;

      if(genPreFSR_Pt->size() == 2 && !tauFlag) {

	for(unsigned int j=0;j<genPhoton_Pt->size();j++){

	  newphoPt.push_back(genPhoton_Pt->at(index1[j]));
	  newphoEta.push_back(genPhoton_Eta->at(index1[j]));
	  newphoPhi.push_back(genPhoton_Phi->at(index1[j]));
	  newphoEn.push_back(genPhoton_En->at(index1[j]));
	}

	for(unsigned int i=0;i<genPostFSR_Pt->size();i++){

	  newgenPt.push_back(genPostFSR_Pt->at(index2[i]));
	  newgenEta.push_back(genPostFSR_Eta->at(index2[i]));
	  newgenPhi.push_back(genPostFSR_Phi->at(index2[i]));
	  newgenEn.push_back(genPostFSR_En->at(index2[i]));
	}

	post1.SetPtEtaPhiE(newgenPt.at(0),newgenEta.at(0),newgenPhi.at(0),newgenEn.at(0));
	post2.SetPtEtaPhiE(newgenPt.at(1),newgenEta.at(1),newgenPhi.at(1),newgenEn.at(1));
	diPost=post1+post2;
	postFSR_Mass=diPost.M();

	for(unsigned int igen = 0; igen < newgenPt.size(); igen++){
	  SumPhotonMom.SetPtEtaPhiE(0.,0.,0.,0.);

	  if(newphoPt.size() >= 1.){
	    for(unsigned int ipho = 0; ipho < newphoPt.size(); ipho++){

	      gen_preFSR.SetPtEtaPhiE(newgenPt.at(igen), newgenEta.at(igen), newgenPhi.at(igen), newgenEn.at(igen));
	      fourmom.SetPtEtaPhiE(0.,0.,0.,0.);

	      dR = deltaR(newphoEta.at(ipho), newphoPhi.at(ipho), newgenEta.at(igen), newgenPhi.at(igen));

	      if(dR < 0.1){

		fourmom.SetPtEtaPhiE(newphoPt.at(ipho), newphoEta.at(ipho), newphoPhi.at(ipho), newphoEn.at(ipho));
		SumPhotonMom = SumPhotonMom + fourmom;

	      }
	    }

	    gen_preFSR = gen_preFSR + SumPhotonMom;

	    gPreFSR_Pt.push_back(gen_preFSR.Pt());
	    gPreFSR_Eta.push_back(gen_preFSR.Eta());
	    gPreFSR_Phi.push_back(gen_preFSR.Phi());
	    gPreFSR_En.push_back(gen_preFSR.Energy());
	  }

	  else {

	    gPreFSR_Pt.push_back(newgenPt.at(igen));
	    gPreFSR_Eta.push_back(newgenEta.at(igen));
	    gPreFSR_Phi.push_back(newgenPhi.at(igen));
	    gPreFSR_En.push_back(newgenEn.at(igen));
	  }
	}

	pre1.SetPtEtaPhiE(gPreFSR_Pt.at(0),gPreFSR_Eta.at(0),gPreFSR_Phi.at(0),gPreFSR_En.at(0));
	pre2.SetPtEtaPhiE(gPreFSR_Pt.at(1),gPreFSR_Eta.at(1),gPreFSR_Phi.at(1),gPreFSR_En.at(1));
	diPre=pre1+pre2;
	preFSR_Mass=diPre.M();

	lumiWeight = lumi_Weight;
	genWeight = theWeight;

	tree->Fill();
      } // pre FSR == 2
    } // event

    file[jentry]->Write();
    file[jentry]->Close();
    
    cout<<""<<endl;

  } // file Loop
}
