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

void FSRCorr_Matrix_PHOTOS() {

  TString workdir;
  std::vector<TFile*> InputFiles_signal_DY;
  int mass[13] = {10,50,50,100,200,400,500,700,800,1000,1500,2000,3000};

  //double xsec[13] = {18609.9/3,5789./3,5789./3,226./3,7.68/3,0.423/3,0.24/3,0.036/3,0.03/3,0.0159/3,0.00201/3,0.00054/3};
  double xsec[13] = {18610./3,5870./3,5870./3,226./3,7.67/3,0.423/3,0.24/3,0.035/3,0.03/3,0.016/3,0.002/3,0.00054/3};
  double sumofWts[12] = {109933846165.897476,50495259821.712626,50495259821.712626,180250845.335368,1431007.165396,88508.734153,41819.599379,9360.877373,4418.256088,3922.345749,464.860914,137.147596};
  //double sumofWts[12] = {329083310509.273560,155672868766.721511,155672868766.721511,540691818.559381,4323599.795268,265931.120467,126228.177749,27794.464027,13222.642455,11737.389431,1394.885489,406.643229};

  workdir = "/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/PHOTOS_FSRCorr/v201611208_1st_ChangeMassScale_HighMassSample/";
  InputFiles_signal_DY.clear();

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
  TFile* file[12];

  for(unsigned int jentry = 0; jentry < nsample; ++jentry) {
    TTree * T1 = (TTree*)InputFiles_signal_DY.at(jentry)->Get("recoTree/DYTree");

    static const int MPSIZE = 2000;

    int evtNum;
    int GENnPair;
    int nGenOthers;
    int nLHEParticle;
    double LHELepton_Px[MPSIZE];
    double LHELepton_Py[MPSIZE];
    double LHELepton_Pz[MPSIZE];
    double LHELepton_E[MPSIZE];
    int LHELepton_ID[MPSIZE];
    int LHELepton_status[MPSIZE];
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
    T1->SetBranchStatus("nLHEParticle", 1);
    T1->SetBranchStatus("LHELepton_Px", 1);
    T1->SetBranchStatus("LHELepton_Py", 1);
    T1->SetBranchStatus("LHELepton_Pz", 1);
    T1->SetBranchStatus("LHELepton_E", 1);
    T1->SetBranchStatus("LHELepton_ID", 1);
    T1->SetBranchStatus("LHELepton_status", 1);
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
    T1->SetBranchAddress("nLHEParticle", &nLHEParticle);
    T1->SetBranchAddress("LHELepton_Px", &LHELepton_Px);
    T1->SetBranchAddress("LHELepton_Py", &LHELepton_Py);
    T1->SetBranchAddress("LHELepton_Pz", &LHELepton_Pz);
    T1->SetBranchAddress("LHELepton_E", &LHELepton_E);
    T1->SetBranchAddress("LHELepton_ID", &LHELepton_ID);
    T1->SetBranchAddress("LHELepton_status", &LHELepton_status);
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

    file[jentry] = new TFile(Form("FSRCorr/altSignal/v2_HighMassSamples/DYEE_M%dto%d.root",mass[jentry],mass[jentry+1]),"RECREATE");
    //file[jentry] = new TFile(Form("/tmp/rchawla/eos/cms/store/group/phys_smp/rchawla/PHOTOS_FSRCorr/v201611208_1st_ChangeMassScale_HighMassSample/Output/DYEE_M%dto%d_part2.root",mass[jentry],mass[jentry+1]),"RECREATE");
    TTree *tree = new TTree("tree"," after preselections tree");

    double lumi_Weight = xsec[jentry]/sumofWts[jentry];
    cout<<"DY Sample: "<<mass[jentry]<<"to"<<mass[jentry+1]<<endl;

    double dR;
    double diLHE_Mass;
    double postFSR_Mass;
    double preFSR_Mass;
    double lumiWeight, genWeight;
    TLorentzVector gen_preFSR;
    TLorentzVector fourmom;
    TLorentzVector SumPhotonMom;
    TLorentzVector lhe1,lhe2,diLHE;
    TLorentzVector pre1,pre2,diPre;
    TLorentzVector post1,post2,diPost;
    TLorentzVector gen1,gen2,diGen;

    // Branch declaration
    tree->Branch("postFSR_Mass", &postFSR_Mass, "postFSR_Mass/D");
    tree->Branch("preFSR_Mass", &preFSR_Mass, "preFSR_Mass/D");
    tree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
    tree->Branch("genWeight", &genWeight, "genWeight/D");

    vector <int> idx1; vector <int> idx2; vector <int> photon_dR;
    vector <double> genPreFSR_Pt; vector <double> genPreFSR_Eta; vector <double> genPreFSR_Phi;

    int nentries = T1->GetEntries();
    //int nentries = 7000000;
    cout<<"entries: "<<nentries<<endl;
    for (unsigned int i = 0; i < nentries; i++) {
      T1->GetEntry(i);

      if(i%1000000 == 0){
	cout << "Events Processed :  " << i << endl;
      }

      dR = 0.;
      postFSR_Mass = -999.; preFSR_Mass = -999.;
      idx1.clear(); idx2.clear();
      genPreFSR_Pt.clear(); genPreFSR_Eta.clear(); genPreFSR_Phi.clear();

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

	for(int j=0; j<nGenOthers; j++){
	  if(abs(GenOthers_mother[j]) == 11) idx2.push_back(j);
	}	

	if((jentry==1 || jentry==2) && diLHE_Mass >= 100.) continue;

	post1.SetPtEtaPhiM(GENLepton_pT[idx1.at(0)],GENLepton_eta[idx1.at(0)],GENLepton_phi[idx1.at(0)],0);
	post2.SetPtEtaPhiM(GENLepton_pT[idx1.at(1)],GENLepton_eta[idx1.at(1)],GENLepton_phi[idx1.at(1)],0);
	diPost=post1+post2;
	postFSR_Mass=diPost.M();

	for(unsigned int igen = 0; igen < 2; igen++){
	  SumPhotonMom.SetPtEtaPhiM(0.,0.,0.,0.);

	  if(nGenOthers >= 1.){
	    //if(idx2.size() >= 1.)

	    for(unsigned int ipho = 0; ipho < nGenOthers; ipho++){
	      //for(unsigned int ipho = 0; ipho < idx2.size(); ipho++)

	      gen_preFSR.SetPtEtaPhiM(GENLepton_pT[idx1.at(igen)],GENLepton_eta[idx1.at(igen)],GENLepton_phi[idx1.at(igen)],0);
	      fourmom.SetPtEtaPhiM(0.,0.,0.,0.);

	      dR = deltaR(GenOthers_eta[ipho],GenOthers_phi[ipho],GENLepton_eta[idx1.at(igen)],GENLepton_phi[idx1.at(igen)]);

	      if(dR < 0.1){

		fourmom.SetPtEtaPhiM(GenOthers_pT[ipho],GenOthers_eta[ipho],GenOthers_phi[ipho],0);
		SumPhotonMom = SumPhotonMom + fourmom;

	      }
	    }

	    gen_preFSR = gen_preFSR + SumPhotonMom;

	    genPreFSR_Pt.push_back(gen_preFSR.Pt());
	    genPreFSR_Eta.push_back(gen_preFSR.Eta());
	    genPreFSR_Phi.push_back(gen_preFSR.Phi());
	  }

	  else {

	    genPreFSR_Pt.push_back(GENLepton_pT[idx1.at(igen)]);
	    genPreFSR_Eta.push_back(GENLepton_eta[idx1.at(igen)]);
	    genPreFSR_Phi.push_back(GENLepton_phi[idx1.at(igen)]);
	  }
	}

	pre1.SetPtEtaPhiM(genPreFSR_Pt.at(0),genPreFSR_Eta.at(0),genPreFSR_Phi.at(0),0);
	pre2.SetPtEtaPhiM(genPreFSR_Pt.at(1),genPreFSR_Eta.at(1),genPreFSR_Phi.at(1),0);
	diPre=pre1+pre2;
	preFSR_Mass=diPre.M();
      } // idx1.size() == 2

      lumiWeight = lumi_Weight;
      genWeight = GENEvt_weight;

      tree->Fill();
    } // event

    file[jentry]->Write();
    file[jentry]->Close();

    cout<<""<<endl;

  } // file Loop
}
