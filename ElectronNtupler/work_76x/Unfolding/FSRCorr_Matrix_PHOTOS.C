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

void SysFSRCorr_Matrix() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  int mass[8] = {1,2,3,4,5,6,7,8};

  double xsec[7] = {18609.9/3, 18609.9/3, 18609.9/3, 18609.9/3, 18609.9/3, 18609.9/3, 6024./3};
  double sumofWts[7] = {136500365161.982367, 136500365161.982367, 136500365161.982367, 136500365161.982367, 136500365161.982367, 136500365161.982367, 48855157880.975021};

  workdir = "/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/PHOTOS_FSRCorr/";
  //workdir = "/afs/cern.ch/user/k/kplee/work/public/v20160922_1st_GENInfo_FSRfromPHOTOS/";
  InputFiles_signal_DY.clear();

  InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_1.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_2.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_3.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_4.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_5.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M10to50_Photos_6.root"));
  InputFiles_signal_DY.push_back(TFile::Open(workdir+"ROOTFile_FlatNutple_DYLL_M50toInf_Photos.root"));

  int nsample = InputFiles_signal_DY.size();
  TFile* file[8];

  for(unsigned int jentry = 6; jentry < nsample; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("DYTree");

    static const int MPSIZE = 2000;

    int evtNum;
    int GENnPair;
    int nGenOthers;
    double GENLepton_phi[MPSIZE];
    double GENLepton_eta[MPSIZE];
    double GENLepton_pT[MPSIZE];
    double GENLepton_Px[MPSIZE];
    double GENLepton_Py[MPSIZE];
    double GENLepton_Pz[MPSIZE];
    int GENLepton_ID[MPSIZE];
    int GENLepton_isHardProcess[MPSIZE];
    int GENLepton_fromHardProcessFinalState[MPSIZE];
    int GENLepton_fromHardProcessDecayed[MPSIZE];
    double GenOthers_phi[MPSIZE];
    double GenOthers_eta[MPSIZE];
    double GenOthers_pT[MPSIZE];
    double GenOthers_Px[MPSIZE];
    double GenOthers_Py[MPSIZE];
    double GenOthers_Pz[MPSIZE];
    double GenOthers_mother[MPSIZE];
    double GENEvt_weight;

    T1->SetBranchStatus("*", 0);
    T1->SetBranchStatus("evtNum", 1);
    T1->SetBranchStatus("GENnPair", 1);
    T1->SetBranchStatus("nGenOthers", 1);
    T1->SetBranchStatus("GENLepton_pT", 1);
    T1->SetBranchStatus("GENLepton_eta", 1);
    T1->SetBranchStatus("GENLepton_phi", 1);
    T1->SetBranchStatus("GENLepton_ID", 1);
    T1->SetBranchStatus("GENLepton_isHardProcess", 1);
    T1->SetBranchStatus("GENLepton_fromHardProcessFinalState", 1);
    T1->SetBranchStatus("GenOthers_pT", 1);
    T1->SetBranchStatus("GenOthers_eta", 1);
    T1->SetBranchStatus("GenOthers_phi", 1);
    T1->SetBranchStatus("GenOthers_mother", 1);
    T1->SetBranchStatus("GENEvt_weight", 1);

    T1->SetBranchAddress("evtNum", &evtNum);
    T1->SetBranchAddress("GENnPair", &GENnPair);
    T1->SetBranchAddress("nGenOthers", &nGenOthers);
    T1->SetBranchAddress("GENLepton_pT", &GENLepton_pT);
    T1->SetBranchAddress("GENLepton_eta", &GENLepton_eta);
    T1->SetBranchAddress("GENLepton_phi", &GENLepton_phi);
    T1->SetBranchAddress("GENLepton_ID", &GENLepton_ID);
    T1->SetBranchAddress("GENLepton_isHardProcess", &GENLepton_isHardProcess);
    T1->SetBranchAddress("GENLepton_fromHardProcessFinalState", &GENLepton_fromHardProcessFinalState);
    T1->SetBranchAddress("GenOthers_pT", &GenOthers_pT);
    T1->SetBranchAddress("GenOthers_eta", &GenOthers_eta);
    T1->SetBranchAddress("GenOthers_phi", &GenOthers_phi);
    T1->SetBranchAddress("GenOthers_mother", &GenOthers_mother);
    T1->SetBranchAddress("GENEvt_weight", &GENEvt_weight);

    //file[jentry] = new TFile(Form("FSRCorr/test_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    file[jentry] = new TFile(Form("FSRCorr/altSignal/Checks/DYEE_%d_5.root",mass[jentry]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    double dR;
    int nGenPhotons, nGenPho_dR;
    double gen1Post_Pt, gen2Post_Pt, gen1Post_Eta, gen2Post_Eta, gen1Post_Phi, gen2Post_Phi, postFSR_Pt, postFSR_Rap, postFSR_Mass;
    double gen1Pre_Pt, gen2Pre_Pt, gen1Pre_Eta, gen2Pre_Eta, gen1Pre_Phi, gen2Pre_Phi, preFSR_Pt, preFSR_Rap, preFSR_Mass;
    double photon_Pt, photonPt_dR, photonEn_dR;
    double lumiWeight, genWeight;
    vector<double> delta_R;
    TLorentzVector gen_preFSR;
    TLorentzVector fourmom;
    TLorentzVector SumPhotonMom;
    TLorentzVector pre1,pre2,diPre;
    TLorentzVector post1,post2,diPost;
    TLorentzVector gen1,gen2,diGen;

    // Branch declaration
    tree->Branch("gen1Post_Pt", &gen1Post_Pt, "gen1Post_Pt/D");
    tree->Branch("gen2Post_Pt", &gen2Post_Pt, "gen2Post_Pt/D");
    tree->Branch("gen1Post_Eta", &gen1Post_Eta, "gen1Post_Eta/D");
    tree->Branch("gen2Post_Eta", &gen2Post_Eta, "gen2Post_Eta/D");
    tree->Branch("gen1Post_Phi", &gen1Post_Phi, "gen1Post_Phi/D");
    tree->Branch("gen2Post_Phi", &gen2Post_Phi, "gen2Post_Phi/D");
    tree->Branch("postFSR_Pt", &postFSR_Pt, "postFSR_Pt/D");
    tree->Branch("postFSR_Rap", &postFSR_Rap, "postFSR_Rap/D");
    tree->Branch("postFSR_Mass", &postFSR_Mass, "postFSR_Mass/D");
    tree->Branch("gen1Pre_Pt", &gen1Pre_Pt, "gen1Pre_Pt/D");
    tree->Branch("gen2Pre_Pt", &gen2Pre_Pt, "gen2Pre_Pt/D");
    tree->Branch("gen1Pre_Eta", &gen1Pre_Eta, "gen1Pre_Eta/D");
    tree->Branch("gen2Pre_Eta", &gen2Pre_Eta, "gen2Pre_Eta/D");
    tree->Branch("gen1Pre_Phi", &gen1Pre_Phi, "gen1Pre_Phi/D");
    tree->Branch("gen2Pre_Phi", &gen2Pre_Phi, "gen2Pre_Phi/D");
    tree->Branch("preFSR_Pt", &preFSR_Pt, "preFSR_Pt/D");
    tree->Branch("preFSR_Rap", &preFSR_Rap, "preFSR_Rap/D");
    tree->Branch("preFSR_Mass", &preFSR_Mass, "preFSR_Mass/D");

    tree->Branch("delta_R", &delta_R);
    tree->Branch("nGenPhotons", &nGenPhotons, "nGenPhotons/I");
    tree->Branch("nGenPho_dR", &nGenPho_dR, "nGenPho_dR/I");
    tree->Branch("photon_Pt", &photon_Pt, "photon_Pt/D");
    //tree->Branch("photon_En", &photon_En, "photon_En/D");
    tree->Branch("photonPt_dR", &photonPt_dR, "photonPt_dR/D");
    tree->Branch("photonEn_dR", &photonEn_dR, "photonEn_dR/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");

    vector <int> idx1; vector <int> idx2; vector <int> photon_dR;
    vector <double> genPreFSR_Pt; vector <double> genPreFSR_Eta; vector <double> genPreFSR_Phi;

    int nentries = T1->GetEntries();
    //int nentries = 10;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int jentry=8000000; jentry < nentries; jentry++) {
      T1->GetEntry(jentry);

      if(jentry%1000000 == 0){
	cout << "Events Processed :  " << jentry << endl;
      }

      dR = 0.;
      nGenPhotons = -999; nGenPho_dR = -999;
      gen1Post_Pt = -999.; gen2Post_Pt = -999.; gen1Post_Eta = -999.; gen2Post_Eta = -999.; gen1Post_Phi = -999.; gen2Post_Phi = -999.; postFSR_Pt = -999.; postFSR_Rap = -999.; postFSR_Mass = -999.;
      gen1Pre_Pt = -999.; gen2Pre_Pt = -999.; gen1Pre_Eta = -999.; gen2Pre_Eta = -999.; gen1Pre_Phi = -999.; gen2Pre_Phi = -999.; preFSR_Pt = -999.; preFSR_Rap = -999.; preFSR_Mass = -999.;
      photon_Pt = 0.; photonPt_dR = 0.; photonEn_dR = 0.;
      idx1.clear(); idx2.clear(); photon_dR.clear(); delta_R.clear();
      genPreFSR_Pt.clear(); genPreFSR_Eta.clear(); genPreFSR_Phi.clear();

      for(int i=0; i<GENnPair; i++){
	if(abs(GENLepton_ID[i]) == 11 && GENLepton_fromHardProcessFinalState[i] == 1) idx1.push_back(i);
      }

      for(int j=0; j<nGenOthers; j++){
	if(abs(GenOthers_mother[j]) == 11) idx2.push_back(j);
      }	

      if(idx1.size() == 2){

	//nGenPhotons = nGenOthers;
	nGenPhotons = idx2.size();

	gen1Post_Pt = GENLepton_pT[idx1.at(0)];
	gen2Post_Pt = GENLepton_pT[idx1.at(1)];
	gen1Post_Eta = GENLepton_eta[idx1.at(0)];
	gen2Post_Eta = GENLepton_eta[idx1.at(1)];
	gen1Post_Phi = GENLepton_phi[idx1.at(0)];
	gen2Post_Phi = GENLepton_phi[idx1.at(1)];

	post1.SetPtEtaPhiM(GENLepton_pT[idx1.at(0)],GENLepton_eta[idx1.at(0)],GENLepton_phi[idx1.at(0)],0);
	post2.SetPtEtaPhiM(GENLepton_pT[idx1.at(1)],GENLepton_eta[idx1.at(1)],GENLepton_phi[idx1.at(1)],0);

	diPost=post1+post2;
	postFSR_Mass=diPost.M();
	postFSR_Pt=diPost.Pt();
	postFSR_Rap=diPost.Rapidity();

	//cout<<"size: "<<nGenOthers<<endl;
	for(int k=0; k<idx2.size(); k++){   // nGenOthers
	  //cout<<"photon pt: "<<photon_Pt<<"+"<<GenOthers_pT[j]<<endl;

	  photon_Pt = photon_Pt + GenOthers_pT[k];
	}

	//cout<<"post FSR Mass: "<<postFSR_Mass<<endl;
	//cout<<"Post FSR PT: "<<GENLepton_pT[idx1.at(0)]<<"   "<<GENLepton_pT[idx1.at(1)]<<endl;
	//cout<<"No. of photons: "<<nGenOthers<<endl;

	//cout<<"No. of Photons: "<<nGenOthers<<endl;
	for(unsigned int igen = 0; igen < idx1.size(); igen++){
	  SumPhotonMom.SetPtEtaPhiM(0.,0.,0.,0.);

	  //if(nGenOthers >= 1.){
	  if(idx2.size() >= 1.){

	    //for(unsigned int ipho = 0; ipho < nGenOthers; ipho++){
	    for(unsigned int ipho = 0; ipho < idx2.size(); ipho++){

	      gen_preFSR.SetPtEtaPhiM(GENLepton_pT[idx1.at(igen)],GENLepton_eta[idx1.at(igen)],GENLepton_phi[idx1.at(igen)],0);
	      //cout<<"ipho: "<<ipho<<"   Pre FSR PT: "<<gen_preFSR.Pt()<<endl;
	      fourmom.SetPtEtaPhiM(0.,0.,0.,0.);

	      dR = deltaR(GenOthers_eta[idx2.at(ipho)],GenOthers_phi[idx2.at(ipho)],GENLepton_eta[idx1.at(igen)],GENLepton_phi[idx1.at(igen)]);
	      //cout<<"dR: "<<dR<<endl;

	      delta_R.push_back(dR);

	      if(dR < 0.1){

		//cout<<"Photons in dR cone: "<<ipho<<endl;
		photon_dR.push_back(ipho);

		fourmom.SetPtEtaPhiM(GenOthers_pT[idx2.at(ipho)],GenOthers_eta[idx2.at(ipho)],GenOthers_phi[idx2.at(ipho)],0);
		//cout<<"entry: "<<jentry<<endl;
		//cout<<"mother ID: "<<GenOthers_mother[ipho]<<"   dR: "<<dR<<endl;
		//cout<<"Photon PT: "<<fourmom.Pt()<<endl;

		SumPhotonMom = SumPhotonMom + fourmom;

		photonPt_dR = photonPt_dR + fourmom.Pt();
		photonEn_dR = photonEn_dR + fourmom.Energy();
		//cout<<"SumPhoton PT: "<<SumPhotonMom.Pt()<<endl;

	      }
	    }

	    gen_preFSR = gen_preFSR + SumPhotonMom;
	    //cout<<"Pre FSR PT added with Photon PT: "<<gen_preFSR.Pt()<<endl;
	    //cout<<""<<endl;

	    genPreFSR_Pt.push_back(gen_preFSR.Pt());
	    genPreFSR_Eta.push_back(gen_preFSR.Eta());
	    genPreFSR_Phi.push_back(gen_preFSR.Phi());
	  }

	  else {

	    genPreFSR_Pt.push_back(GENLepton_pT[idx1.at(igen)]);
	    genPreFSR_Eta.push_back(GENLepton_eta[idx1.at(igen)]);
	    genPreFSR_Phi.push_back(GENLepton_phi[idx1.at(igen)]);
	    //cout<<"igen: "<<igen<<"   Pre FSR PT not added with Photon PT: "<<GENLepton_pT[idx1.at(igen)]<<endl;
	    //cout<<""<<endl;
	  }
	}

	nGenPho_dR = photon_dR.size();
	//cout<<"entry: "<<jentry<< "   No. of photons within dR cone: "<<nGenPho_dR<<endl;

	pre1.SetPtEtaPhiM(genPreFSR_Pt.at(0),genPreFSR_Eta.at(0),genPreFSR_Phi.at(0),0);
	pre2.SetPtEtaPhiM(genPreFSR_Pt.at(1),genPreFSR_Eta.at(1),genPreFSR_Phi.at(1),0);

	gen1Pre_Pt  = genPreFSR_Pt.at(0);
	gen2Pre_Pt  = genPreFSR_Pt.at(1);
	gen1Pre_Eta = genPreFSR_Eta.at(0);
	gen2Pre_Eta = genPreFSR_Eta.at(1);
	gen1Pre_Phi = genPreFSR_Phi.at(0);
	gen2Pre_Phi = genPreFSR_Phi.at(1);

	diPre=pre1+pre2;
	preFSR_Mass=diPre.M();
	preFSR_Pt=diPre.Pt();
	preFSR_Rap=diPre.Rapidity();
	//cout<<"pre FSR Mass: "<<preFSR_Mass<<endl;
	//cout<<""<<endl;
      }

      lumiWeight = lumi_Weight;
      genWeight = GENEvt_weight;

      tree->Fill();
    } // event

    file[jentry]->Write();
    file[jentry]->Close();

  } // file Loop
}
