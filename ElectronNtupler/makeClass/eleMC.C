#define eleMC_cxx
#include "eleMC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
//#include "DataFormats/Math/interface/deltaR.h"

void eleMC::Loop()
{
  /*
     TFile *f1 = TFile::Open("/afs/cern.ch/work/r/rchawla/private/CMSSW_5_3_22/src/RecoLuminosity/LumiDB/PU_DATA/PileupHist_Run2012.root");
     TFile *f2 = TFile::Open("/afs/cern.ch/work/r/rchawla/private/CMSSW_5_3_22/src/Analysis/El_analyzer/Macros/New/MCEE/WTSfiles/MCDist_M50_13TeV.root");

  // data histogram
  TH1F *DATA_puDist = (TH1F*)f1->Get("pileup");
  DATA_puDist->Scale(1/DATA_puDist->Integral());

  // mc histogram
  TH1F *MC_puDist = (TH1F*)f2->Get("MC_puDist");
  MC_puDist->Scale(1/MC_puDist->Integral());

  TH1F *weights_ = (TH1F*)DATA_puDist->Clone("weights_");
  weights_->Divide(MC_puDist);
  */
  TFile *file = new TFile("dyJtoLL_M50.root", "recreate");

  //double Mz = 91.1876;

  int good_elec;
  int n_dR;
  //int good_jets;
  double ZMass, ZY, ZRap, ZEta, ZPt, ZPhi;
  //double deltaR1, deltaR2;
  TLorentzVector ele1,ele2,dielectron,ele12,ele22,dielectron2;

  vector <double> newscEt;
  vector <double> newscEta;
  vector <double> newelePt;
  vector <double> neweleEta;
  vector <double> neweleEnr;
  vector <double> newelePhi;
  vector <double> neweleIso;

  //vector <int>    jetmarker;

  //TH1F *MC_puDistW = new TH1F("MC_puDistW","MC_puDistW",50,0,50);

  TH1F *eta_sc = new TH1F("eta_sc", "eta_sc", 50, -2.5, 2.5);
  TH1F *et_sc = new TH1F("et_sc", "et_sc", 100, 0 ,700);
  TH1F *et_sc_rebin = new TH1F("et_sc_rebin", "et_sc_rebin", 700, 0 ,700);

  TH1F *pt_lead = new TH1F("pt_lead", "pt_lead", 100, 0, 700);
  TH1F *eta_lead = new TH1F("eta_lead", "eta_lead", 50, -2.5, 2.5);
  TH1F *pt_slead = new TH1F("pt_slead", "pt_slead", 100, 0, 700);
  TH1F *eta_slead = new TH1F("eta_slead", "eta_slead", 48, -2.4, 2.4);

  TH1F *et_sc_lead = new TH1F("et_sc_lead", "et_sc_lead", 105, 0 ,700);
  TH1F *et_sc_slead = new TH1F("et_sc_slead", "et_sc_slead", 105, 0 ,700);
  TH1F *eta_sc_lead = new TH1F("eta_sc_lead", "eta_sc_lead", 48, -2.4, 2.4);
  TH1F *eta_sc_slead = new TH1F("eta_sc_slead", "eta_sc_slead", 48, -2.4, 2.4);

  /*TH1F *B_delta_eta = new TH1F("B_delta_eta", "B_delta_eta", 260, -0.01, 0.01);
    TH1F *B_delta_phi = new TH1F("B_delta_phi", "B_delta_phi", 160, -0.08, 0.08);
    TH1F *B_H_E = new TH1F("B_H_E", "B_H_E", 100, 0, 0.1);
    TH1F *B_sigma_eta = new TH1F("B_sigma_eta", "B_sigma_eta", 100, 0, 0.1);
    TH1F *B_E_P = new TH1F("B_E_P", "B_E_P", 400, 0, 0.5);

    TH1F *E_delta_eta = new TH1F("E_delta_eta", "E_delta_eta", 260, -0.01, 0.01);
    TH1F *E_delta_phi = new TH1F("E_delta_phi", "E_delta_phi", 160, -0.08, 0.08);
    TH1F *E_H_E = new TH1F("E_H_E", "E_H_E", 100, 0, 0.1);
    TH1F *E_sigma_eta = new TH1F("E_sigma_eta", "E_sigma_eta", 100, 0, 0.1);
    TH1F *E_E_P = new TH1F("E_E_P", "E_E_P", 400, 0, 0.5);*/

  TH1F *h_nPV = new TH1F("h_nPV", "h_nPV", 50, 0, 50);
  
  TH1F *B_DeltaEta = new TH1F("B_DeltaEta", "B_deltaEta", 20, -0.01, 0.01);
  TH1F *B_DeltaPhi = new TH1F("B_DeltaPhi", "B_deltaPhi", 20, -0.08, 0.08);
  TH1F *B_H_E = new TH1F("B_H_E", "B_H_E", 20, 0, 0.1);
  TH1F *B_SigmaEta = new TH1F("B_SigmaEta", "B_sigmaEta", 40, 0, 0.1);
  TH1F *B_E_P = new TH1F("B_E_P", "B_E_P", 20, 0, 0.5);
  TH1F *B_D0 = new TH1F("B_D0", "B_D0", 40, 0, 0.2);
  TH1F *B_Dz = new TH1F("B_Dz", "B_Dz", 40, 0, 0.2);
  TH1F *B_MissHits = new TH1F("B_MissHits", "B_MissHits", 4, 0, 4);
  TH1F *B_ConVeto = new TH1F("B_ConVeto", "B_ConVeto", 2, 0, 2);
  TH1F *B_fBrem = new TH1F("B_fBrem", "B_fBrem", 40, 0, 1);
  TH1F *B_r9 = new TH1F("B_r9", "B_r9", 40, 0, 1);
  TH1F *B_RelRhoIso = new TH1F("B_RelRhoIso", "B_RelRhoIso", 30, 0, 0.6);
  //TH1F *B_RelEMIso = new TH1F("B_RelEMIso", "B_RelEMIso", 60, 0, 0.6);
  //TH1F *B_RelNeutralIso = new TH1F("B_RelNeutralIso", "B_RelNeutralIso", 60, 0, 0.6);

  TH1F *E_DeltaEta = new TH1F("E_DeltaEta", "E_DeltaEta", 20, -0.01, 0.01);
  TH1F *E_DeltaPhi = new TH1F("E_DeltaPhi", "E_DeltaPhi", 20, -0.08, 0.08);
  TH1F *E_H_E = new TH1F("E_H_E", "E_H_E", 20, 0, 0.1);
  TH1F *E_SigmaEta = new TH1F("E_SigmaEta", "E_SigmaEta", 40, 0.0, 0.1);
  TH1F *E_E_P = new TH1F("E_E_P", "E_E_P", 20, 0, 0.5);
  TH1F *E_D0 = new TH1F("E_D0", "E_D0", 40, 0, 0.2);
  TH1F *E_Dz = new TH1F("E_Dz", "E_Dz", 40, 0, 0.2);
  TH1F *E_MissHits = new TH1F("E_MissHits", "E_MissHits", 4, 0, 4);
  TH1F *E_ConVeto = new TH1F("E_ConVeto", "E_ConVeto", 2, 0, 2);
  TH1F *E_fBrem = new TH1F("E_fBrem", "E_fBrem", 40, 0, 1);
  TH1F *E_r9 = new TH1F("E_r9", "E_r9", 40, 0.0, 1);
  TH1F *E_RelRhoIso = new TH1F("E_RelRhoIso", "E_RelRhoIso", 30, 0, 0.6);
  //TH1F *E_RelEMIso = new TH1F("E_RelEMIso", "E_RelEMIso", 60, 0, 0.6);
  //TH1F *E_RelNeutralIso = new TH1F("E_RelNeutralIso", "E_RelNeutralIso", 60, 0, 0.6);

  //Double_t xbins[7] = {15.,30.,60.,120.,240.,600.,2000.};

  //TH1D *z_mass_gen = new TH1D("z_mass_gen", "z_mass_gen", 6, xbins);
  //TH1D *z_mass_gen_cut = new TH1D("z_mass_gen_cut", "z_mass_gen_cut", 6, xbins);
  //TH1D *z_mass_gen = new TH1D("z_mass_gen", "z_mass_gen", 100, 0., 200.);
  //TH1D *z_mass_gen_cut = new TH1D("z_mass_gen_cut", "z_mass_gen_cut", 100, 0., 200.);
  TH1D *z_mass_1 = new TH1D("z_mass_1", "z_mass_1", 1000, 0, 1000);
  TH1D *z_mass_2 = new TH1D("z_mass_2", "z_mass_2", 500, 0, 1000);
  TH1D *z_mass = new TH1D("z_mass", "z_mass", 60, 60, 120);
  TH1D *z_pt = new TH1D("z_pt", "z_pt", 100, 0, 700);
  TH1D *z_y = new TH1D("z_y", "z_y", 150, -15, 15);
  TH1D *z_rap = new TH1D("z_rap", "z_rap", 60, -3, 3);
  TH1D *z_eta = new TH1D("z_eta", "z_eta", 60, -8, 8);
  TH1D *z_phi = new TH1D("z_phi", "z_phi", 50, -3.5, 3.5);

  //TH2D *hRM = new TH2D("Response Matrix", "Response Matrix", 6, xbins, 6, xbins);
  //TH2D *hRM_match = new TH2D("Response Matrix_match", "Response Matrix_match", 6, xbins, 6, xbins);

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"entries: "<<nentries<<endl;

  n_dR = 0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    //cout<<"entry: "<<jentry<<endl;
    //cout<<"No. of electrons: "<<nEle<<endl;

    good_elec = 0;
    //int bin = 0;
    //double Weights = 0.0;

    //bin = weights_->GetXaxis()->FindBin(PUI);
    //Weights = weights_->GetBinContent(bin);
    //MC_puDistW->Fill(PUI,Weights);

    newscEt.clear();
    newscEta.clear();
    newelePt.clear();
    neweleEta.clear();
    neweleEnr.clear();
    newelePhi.clear();
    neweleIso.clear();

    h_nPV->Fill(pVtx);

    for(int k=0;k<nEle;k++){

      bool isBarrel = (fabs(etaSC->at(k)) <= 1.479);
      bool isEndcap = (fabs(etaSC->at(k)) > 1.479 && fabs(etaSC->at(k)) < 2.5);

      /*bool isEtaB = false;
	bool isPhiB = false;
	bool isSigmaB = false;
	bool isH_EB = false;
	bool isE_PB = false;

	bool isEtaE = false;
	bool isPhiE = false;
	bool isSigmaE = false;
	bool isH_EE = false;
	bool isE_PE = false;

	if(isBarrel && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dPhiIn->at(k)) < 0.035973 && full5x5_sigmaIetaIeta->at(k) < 0.009996 && hOverE->at(k) < 0.050537 && abs(d0->at(k)) < 0.012235 && abs(dz->at(k)) < 0.042020 && ooEmooP->at(k) < 0.091942 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.107587) isEtaB = true;
	if(isEtaB) B_delta_eta->Fill(dEtaIn->at(k)); 

	if(isEndcap && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dPhiIn->at(k)) < 0.067879 && full5x5_sigmaIetaIeta->at(k) < 0.030135 && hOverE->at(k) < 0.086782 && abs(d0->at(k)) < 0.036719 && abs(dz->at(k)) < 0.138142 && ooEmooP->at(k) < 0.100683 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.113254) isEtaE = true;
	if(isEtaE) E_delta_eta->Fill(dEtaIn->at(k));

	if(isBarrel && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.008925 && full5x5_sigmaIetaIeta->at(k) < 0.009996 && hOverE->at(k) < 0.050537 && abs(d0->at(k)) < 0.012235 && abs(dz->at(k)) < 0.042020 && ooEmooP->at(k) < 0.091942 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.107587) isPhiB = true;
	if(isPhiB) B_delta_phi->Fill(dPhiIn->at(k));

	if(isEndcap && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.007429 && full5x5_sigmaIetaIeta->at(k) < 0.030135 && hOverE->at(k) < 0.086782 && abs(d0->at(k)) < 0.036719 && abs(dz->at(k)) < 0.138142 && ooEmooP->at(k) < 0.100683  && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.113254) isPhiE = true;
	if(isPhiE) E_delta_phi->Fill(dPhiIn->at(k));

	if(isBarrel && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.008925 && abs(dPhiIn->at(k)) < 0.035973 && hOverE->at(k) < 0.050537 && abs(d0->at(k)) < 0.012235 && abs(dz->at(k)) < 0.042020 && ooEmooP->at(k) < 0.091942 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.107587) isSigmaB = true; 
	if(isSigmaB) B_sigma_eta->Fill(full5x5_sigmaIetaIeta->at(k));

	if(isEndcap && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.007429 && abs(dPhiIn->at(k)) < 0.067879 && hOverE->at(k) < 0.086782  && abs(d0->at(k)) < 0.036719 && abs(dz->at(k)) < 0.138142 && ooEmooP->at(k) < 0.100683 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.113254) isSigmaE = true;     
	if(isSigmaE) E_sigma_eta->Fill(full5x5_sigmaIetaIeta->at(k));

	if(isBarrel && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.008925 && abs(dPhiIn->at(k)) < 0.035973 && full5x5_sigmaIetaIeta->at(k) < 0.009996 && abs(d0->at(k)) < 0.012235 && abs(dz->at(k)) < 0.042020 && ooEmooP->at(k) < 0.091942 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.107587) isH_EB = true;
	if(isH_EB) B_H_E->Fill(hOverE->at(k));

	if(isEndcap && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.007429 && abs(dPhiIn->at(k)) < 0.067879 && full5x5_sigmaIetaIeta->at(k) < 0.030135 && abs(d0->at(k)) < 0.036719 && abs(dz->at(k)) < 0.138142 && ooEmooP->at(k) < 0.100683 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.113254) isH_EE = true;
	if(isH_EE) E_H_E->Fill(hOverE->at(k));

	if(isBarrel && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.008925 && abs(dPhiIn->at(k)) < 0.035973 && full5x5_sigmaIetaIeta->at(k) < 0.009996 && hOverE->at(k) < 0.050537 && abs(d0->at(k)) < 0.012235 && abs(dz->at(k)) < 0.042020 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.107587) isE_PB = true;
	if(isE_PB) B_E_P->Fill(ooEmooP->at(k));

	if(isEndcap && pt->at(k) > 15. && etaSC->at(k) <2.4 && abs(dEtaIn->at(k)) < 0.007429 && abs(dPhiIn->at(k)) < 0.067879 && full5x5_sigmaIetaIeta->at(k) < 0.030135 && hOverE->at(k) < 0.086782 && abs(d0->at(k)) < 0.036719 && abs(dz->at(k)) < 0.138142 && expectedMissingInnerHits->at(k) <= 1 && passConversionVeto->at(k) == 1 && isoRho->at(k) < 0.113254) isE_PE = true;
	if(isE_PE) E_E_P->Fill(ooEmooP->at(k));*/

      //cout<<"pt: "<<pt->at(k)<<endl;
      //cout<<" ID: "<<passMediumId->at(k)<<endl;

      //bool isBarrel = (fabs(sc_eta->at(index[k])) <= 1.479);
      //bool isEndcap = (fabs(sc_eta->at(index[k])) > 1.479 && fabs(sc_eta->at(index[k])) < 2.5);

      if(pt->at(k) < 15.) continue;
      //if(passMediumId->at(k) == 0) continue;

      //if(doubleElectron){
	if(isBarrel){
	  B_DeltaEta->Fill(dEtaIn->at(k));
	  B_DeltaPhi->Fill(dPhiIn->at(k));
	  B_SigmaEta->Fill(full5x5_sigmaIetaIeta->at(k));
	  B_H_E->Fill(hOverE->at(k));
	  B_E_P->Fill(ooEmooP->at(k));
	  B_D0->Fill(d0->at(k));
	  B_Dz->Fill(dz->at(k));
	  B_MissHits->Fill(expectedMissingInnerHits->at(k));
	  B_ConVeto->Fill(passConversionVeto->at(k));
	  B_fBrem->Fill(brem->at(k));
	  B_r9->Fill(r9->at(k));
	  B_RelRhoIso->Fill(isoRho->at(k));
	  //B_RelEMIso->Fill(isoNeutralHadrons->at(k));
	  //B_RelNeutralIso->Fill(isoPhotons->at(k));

	}

	if(isEndcap){
	  E_DeltaEta->Fill(dEtaIn->at(k));
	  E_DeltaPhi->Fill(dPhiIn->at(k));
	  E_SigmaEta->Fill(full5x5_sigmaIetaIeta->at(k));
	  E_H_E->Fill(hOverE->at(k));
	  E_E_P->Fill(ooEmooP->at(k));
	  E_D0->Fill(d0->at(k));
	  E_Dz->Fill(dz->at(k));
	  E_MissHits->Fill(expectedMissingInnerHits->at(k));
	  E_ConVeto->Fill(passConversionVeto->at(k));
	  E_fBrem->Fill(brem->at(k));
	  E_r9->Fill(r9->at(k));
	  E_RelRhoIso->Fill(isoRho->at(k));
	  //E_RelEMIso->Fill(isoNeutralHadrons->at(k));
	  //E_RelNeutralIso->Fill(isoPhotons->at(k));
	}
      //}
      
      //bool real_electron = true;
      //      if(real_electron)
      //      {
      
      if(passMediumId->at(k) == 0) continue;

      good_elec = good_elec + 1;

      newscEt.push_back(etSC->at(k));
      newscEta.push_back(etaSC->at(k));
      newelePt.push_back(pt->at(k));
      neweleEta.push_back(eta->at(k));
      neweleEnr.push_back(energy->at(k));
      newelePhi.push_back(phi->at(k));
      neweleIso.push_back(isoRho->at(k));
      //      }
    }

    if(good_elec<2) continue;

    if(fabs(neweleEta[0]) >= 2.4 && fabs(neweleEta[1]) >= 2.4 && newelePt[0] <= 20.) continue;
    if(!doubleElectron) continue;
    /*
       if(fabs(newscEta[0]) <= 1.479 && (neweleIso[0]/newelePt[0]) >= 0.15) continue;
       if(fabs(newscEta[1]) <= 1.479 && (neweleIso[1]/newelePt[1]) >= 0.15) continue;

       if((fabs(newscEta[0]) > 1.479 && fabs(newscEta[0]) < 2.5) && newelePt[0] > 20. && (neweleIso[0]/newelePt[0]) >= 0.15) continue;
       if((fabs(newscEta[1]) > 1.479 && fabs(newscEta[1]) < 2.5) && newelePt[1] > 20. && (neweleIso[1]/newelePt[1]) >= 0.15) continue;

       if((fabs(newscEta[0]) > 1.479 && fabs(newscEta[0]) < 2.5) && newelePt[0] < 20. && (neweleIso[0]/newelePt[0]) >= 0.10) continue;
       if((fabs(newscEta[1]) > 1.479 && fabs(newscEta[1]) < 2.5) && newelePt[1] < 20. && (neweleIso[1]/newelePt[1]) >= 0.10) continue;*/

    //cout<<"new Pt size: "<<newelePt.size()<<endl;
    /*bool leg1_e1 = false;
      bool leg1_e2 = false;
      bool leg2_e1 = false;
      bool leg2_e2 = false;

      for(unsigned int j = 0; j < pt_f1->size(); j++){
      double dR1 = deltaR(neweleEta[0],newelePhi[0],eta_f1->at(j),phi_f1->at(j));
      double dR2 = deltaR(neweleEta[1],newelePhi[1],eta_f1->at(j),phi_f1->at(j));
      if (dR1 < 0.1) leg1_e1 = true;
      if (dR2 < 0.1) leg1_e2 = true;
      }

      for(unsigned int k = 0; k < pt_f2->size(); k++){
      double dR3 = deltaR(neweleEta[0],newelePhi[0],eta_f2->at(k),phi_f2->at(k));
      double dR4 = deltaR(neweleEta[1],newelePhi[1],eta_f2->at(k),phi_f2->at(k));
      if (dR3 < 0.1) leg2_e1 = true;
      if (dR4 < 0.1) leg2_e2 = true;
      }*/

    //if((leg1_e1 && leg2_e2) || (leg1_e2 && leg2_e1)){

    //++n_dR;
    for(unsigned int i=0; i<newelePt.size(); i++)
    {
      eta_sc->Fill(newscEta[i]);
      et_sc->Fill(newscEt[i]);
      et_sc_rebin->Fill(newscEt[i]);
    }

    pt_lead->Fill(newelePt[0]);
    eta_lead->Fill(neweleEta[0]);
    pt_slead->Fill(newelePt[1]);
    eta_slead->Fill(neweleEta[1]);

    et_sc_lead->Fill(newscEt[0]);
    et_sc_slead->Fill(newscEt[1]);
    eta_sc_lead->Fill(newscEta[0]);
    eta_sc_slead->Fill(newscEta[1]);

    ele1.SetPtEtaPhiE(newelePt[0],neweleEta[0],newelePhi[0],neweleEnr[0]);
    ele2.SetPtEtaPhiE(newelePt[1],neweleEta[1],newelePhi[1],neweleEnr[1]);

    dielectron=ele1+ele2;
    ZMass = dielectron.M();
    ZPt = dielectron.Pt();
    ZY = dielectron.Y();
    ZRap = dielectron.Rapidity();
    ZEta = dielectron.Eta();
    ZPhi = dielectron.Phi();

    //cout<<"Event: "<<jentry<<"   Z mass gen: "<<Z_mass[0]<<"  Z mass reco: "<<ZMass<<endl;

    z_mass_1->Fill(ZMass);
    z_mass_2->Fill(ZMass);
    z_mass->Fill(ZMass);
    z_pt->Fill(ZPt);
    z_y->Fill(ZY);
    z_rap->Fill(ZRap);
    z_eta->Fill(ZEta);
    z_phi->Fill(ZPhi);
    //} // dR matching

    //hRM->Fill(Mee_post,ZMass);
    //if(deltaR1 < 0.3 && deltaR2 < 0.3) hRM_match->Fill(Mee_post,ZMass);
    //hRM->Draw("colz");
    //cout<<"Event: "<<jentry<<"   Z mass gen: "<<Mee<<"  Z mass reco: "<<ZMass<<endl;

  } // event

  //cout<<"n_dR: "<<n_dR<<endl;
  file->Write();
  file->Close();
}
