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

  double xsec[11] = {18609.9/3,5789./3,226./3,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
  double sumofWts[11] = {771365896802.749878,144503601049.458435,219884886.214768,7008726.024242,122990.185607,69998.861766,11080.956088,9415.627747,4893.463037,664.033593,159.615701};

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
    TTree *tree = new TTree("tree"," after preselections tree");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    double dR, massGPre;
    int nGenPhotons, nGenPho_dR;
    double gen1Post_Pt, gen2Post_Pt, gen1Post_Eta, gen2Post_Eta, gen1Post_Phi, gen2Post_Phi, postFSR_Pt, postFSR_Rap, postFSR_Mass;
    double gen1Pre_Pt, gen2Pre_Pt, gen1Pre_Eta, gen2Pre_Eta, gen1Pre_Phi, gen2Pre_Phi, preFSR_Pt, preFSR_Rap, preFSR_Mass;
    double photon_Pt, photon_En, photonPt_dR, photonEn_dR;
    double lumiWeight, genWeight;
    vector<double> delta_R;
    TLorentzVector gen_preFSR;
    TLorentzVector fourmom;
    TLorentzVector SumPhotonMom;
    TLorentzVector pre1,pre2,diPre;
    TLorentzVector post1,post2,diPost;
    TLorentzVector gen1,gen2,diGen;
    int phoindexele1, phoindexele2;

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
    tree->Branch("photon_En", &photon_En, "photon_En/D");
    tree->Branch("photonPt_dR", &photonPt_dR, "photonPt_dR/D");
    tree->Branch("photonEn_dR", &photonEn_dR, "photonEn_dR/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");

    vector <int> idx_gen; vector <int> idx_pho; vector <int> photon_dR;
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
      nGenPhotons = -999; nGenPho_dR = -999;
      gen1Post_Pt = -999.; gen2Post_Pt = -999.; gen1Post_Eta = -999.; gen2Post_Eta = -999.; gen1Post_Phi = -999.; gen2Post_Phi = -999.; postFSR_Pt = -999.; postFSR_Rap = -999.; postFSR_Mass = -999.;
      gen1Pre_Pt = -999.; gen2Pre_Pt = -999.; gen1Pre_Eta = -999.; gen2Pre_Eta = -999.; gen1Pre_Phi = -999.; gen2Pre_Phi = -999.; preFSR_Pt = -999.; preFSR_Rap = -999.; preFSR_Mass = -999.;
      photon_Pt = 0.; photon_En = 0.; photonPt_dR = 0.; photonEn_dR = 0.;
      massGPre = 0.;
      //phoindexele1 = 0; phoindexele2 = 0;
      
      idx_gen.clear(); idx_pho.clear(); photon_dR.clear(); delta_R.clear();
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

      /*if(genPostFSR_Pt->size() == 2){
	if(genPreFSR_Pt->size() != 2) cout<<"check the event: "<<i<<endl;
      }*/

      //if(!tauFlag == 0) cout<<"check the event with taus: "<<i<<endl;

      if(genPreFSR_Pt->size() == 2 && !tauFlag) {

	nGenPhotons = genPhoton_Pt->size();
	//cout<<"entry: "<<i<<"   No of photons: "<<genPhoton_Pt->size()<<endl;

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

	gen1Post_Pt  = newgenPt.at(0);
	gen2Post_Pt  = newgenPt.at(1);
	gen1Post_Eta = newgenEta.at(0);
	gen2Post_Eta = newgenEta.at(1);
	gen1Post_Phi = newgenPhi.at(0);
	gen2Post_Phi = newgenPhi.at(1);

	post1.SetPtEtaPhiE(newgenPt.at(0),newgenEta.at(0),newgenPhi.at(0),newgenEn.at(0));
	post2.SetPtEtaPhiE(newgenPt.at(1),newgenEta.at(1),newgenPhi.at(1),newgenEn.at(1));
	diPost=post1+post2;
	postFSR_Mass=diPost.M();
	postFSR_Pt=diPost.Pt();
	postFSR_Rap=diPost.Rapidity();

	//cout<<"size: "<<newphoPt.size()<<endl;
	for(unsigned int k=0;k<newphoPt.size();k++){

	  //cout<<"photon pt: "<<photon_Pt<<"+"<<newphoPt.at(k)<<endl;
	  
	  photon_Pt = photon_Pt + newphoPt.at(k);
	  photon_En = photon_En + newphoEn.at(k);
	
	}

	//cout<<"Total photon pt: "<<photon_Pt<<endl;
	//cout<<""<<endl;
	
	//if(newgenPt.size() != 2) cout<<"Wrong "<<"   size: "<<newgenPt.size()<<endl;

	for(unsigned int igen = 0; igen < newgenPt.size(); igen++){
	  SumPhotonMom.SetPtEtaPhiE(0.,0.,0.,0.);

	  phoindexele1 = -999; phoindexele2 = -999;
	  //if(igen==0) cout<<"I Loop : "<<igen<<endl;
	  //if(igen==1) cout<<"II Loop : "<<igen<<endl;

	  if(newphoPt.size() >= 1.){
	    for(unsigned int ipho = 0; ipho < newphoPt.size(); ipho++){

	      gen_preFSR.SetPtEtaPhiE(newgenPt.at(igen), newgenEta.at(igen), newgenPhi.at(igen), newgenEn.at(igen));
	      //cout<<"ipho: "<<ipho<<"   Pre FSR PT: "<<gen_preFSR.Pt()<<endl;
	      fourmom.SetPtEtaPhiE(0.,0.,0.,0.);

	      dR = deltaR(newphoEta.at(ipho), newphoPhi.at(ipho), newgenEta.at(igen), newgenPhi.at(igen));
	      delta_R.push_back(dR);

	      //phoindex.push_back(ipho);
	      //genindex.push_back(igen);

	      //cout<<"delta R: "<<dR<<endl;
	      //cout<<""<<endl;

	      if(dR < 0.1){

		//cout<<"Photons in dR cone: "<<ipho<<endl;
		
		photon_dR.push_back(ipho);

		idx_gen.push_back(igen);
		idx_pho.push_back(ipho);
		fourmom.SetPtEtaPhiE(newphoPt.at(ipho), newphoEta.at(ipho), newphoPhi.at(ipho), newphoEn.at(ipho));
		//cout<<"Photon PT: "<<fourmom.Pt()<<endl;
		SumPhotonMom = SumPhotonMom + fourmom;
		//cout<<"SumPhoton PT: "<<SumPhotonMom.Pt()<<endl;

		photonPt_dR = photonPt_dR + fourmom.Pt();
		photonEn_dR = photonEn_dR + fourmom.Energy();

		//cout<<"photonPt_dR: "<<photonPt_dR<<endl;
		//cout<<""<<endl;

		//cout<<"photon index: "<<ipho<<endl;
		if(igen==0) phoindexele1 = ipho;
		if(igen==1) phoindexele2 = ipho;
		//if(phoindexele1==phoindexele2) cout<<"index photon: "<<ipho<<endl;
	      }
	    }

	    gen_preFSR = gen_preFSR + SumPhotonMom;
	    //cout<<"Pre FSR PT added with Photon PT: "<<gen_preFSR.Pt()<<endl;
	    //cout<<""<<endl;

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

	//cout<<"********************************"<<endl;

	nGenPho_dR = photon_dR.size();

	pre1.SetPtEtaPhiE(gPreFSR_Pt.at(0),gPreFSR_Eta.at(0),gPreFSR_Phi.at(0),gPreFSR_En.at(0));
	pre2.SetPtEtaPhiE(gPreFSR_Pt.at(1),gPreFSR_Eta.at(1),gPreFSR_Phi.at(1),gPreFSR_En.at(1));

	gen1Pre_Pt  = gPreFSR_Pt.at(0);
	gen2Pre_Pt  = gPreFSR_Pt.at(1);
	gen1Pre_Eta = gPreFSR_Eta.at(0);
	gen2Pre_Eta = gPreFSR_Eta.at(1);
	gen1Pre_Phi = gPreFSR_Phi.at(0);
	gen2Pre_Phi = gPreFSR_Phi.at(1);

	diPre=pre1+pre2;
	preFSR_Mass=diPre.M();
	preFSR_Pt=diPre.Pt();
	preFSR_Rap=diPre.Rapidity();

	//cout<<"pre FSR Mass: "<<preFSR_Mass<<endl;
	//cout<<""<<endl;


	lumiWeight = lumi_Weight;
	genWeight = theWeight;

	tree->Fill();
      } // pre FSR == 2
    } // event

    file[jentry]->Write();
    file[jentry]->Close();

  } // file Loop
}
