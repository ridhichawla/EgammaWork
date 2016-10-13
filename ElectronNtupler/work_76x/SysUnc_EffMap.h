#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
#include "/home/ridhi/RooUnfold/src/RooUnfoldBayes.h"
#include "/home/ridhi/RooUnfold/src/RooUnfold.h"

#define nEtaBin 10
#define nPtBin 5
#define nEffMap 500

class SysUncTools_Base
{
public:
	Int_t nMassBin;
	Double_t MassBinEdges[44];

	SysUncTools_Base()
	{
		// -- Setting the number of mass bin & mass bin edges -- //
		nMassBin = 43;
  		Double_t temp[43+1] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,
							91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,
							220,243,273,320,380,440,510,600,700,830,1000,1500,3000};
		for(Int_t i=0; i<44; i++)
			MassBinEdges[i] = temp[i];

	}

	void Calculation_DiffXsec(TH1D* h_yield, RooUnfoldResponse* UnfoldRes1,
					TH1D* h_mass_SFCorr, TH1D* h_AccEff,
					RooUnfoldResponse* UnfoldRes3, TH1D *h_xSec_dM_FSRCorr)
   	{


		// -- Unfolding correction for detector resolution: apply on the data before applying Acc*Eff correction -- //
        	RooUnfoldBayes *UnfoldBayes1 = new RooUnfoldBayes(UnfoldRes1, h_yield, 15); //Unfolding Data Histo
		h_yield_Unfolded = (TH1D*)UnfoldBayes1->Hreco();

		
		// -- Apply Efficiency correction factors to each Yield -- //
        	ApplyFinal_EffCorr(h_yield_Unfolded, h_mass_SFCorr, h_yield_Unfolded_SFCorr); //Apply SF Corrections
	
		// -- Acc*Eff Correction -- //
		ApplyFinal_EffCorr(h_yield_Unfolded_SFCorr, h_AccEff, h_yield_Unfolded_AccEff); //Apply Acc*Eff Corrections
	
		// -- FSR Correction -- //
		RooUnfoldBayes *UnfoldBayes3 = new RooUnfoldBayes(UnfoldRes3, h_yield_Unfolded_AccEff, 15);
		TH1D *h_FSRCorrected = (TH1D*)UnfoldBayes3->Hreco();
	
		// -- Obtain differential cross section -- //
		h_xSec_dM_FSRCorr->Sumw2();
		Obtain_dSigma_dM(h_FSRCorrected, h_xSec_dM_FSRCorr); // -- differential cross section -- //
		h_xSec_dM_FSRCorr->Scale(1/2316.969 ); // -- cross section -- //
	}


	//======================== Function to compute Bkg Subtracted Histogram
	void ObtainYieldHistogram(TH1D *h_Data, TH1D *h_EMu, TH1D *h_WZ, TH1D *h_ZZ, TH1D *h_Dijet, TH1D *h_Wjet, TH1D *Final) {

		for(int i=1;i<44; i++){
	
			double data_org = h_Data->GetBinContent(i);
			double data_org_err = h_Data->GetBinError(i);
		
			double emu = h_EMu->GetBinContent(i);
			double emu_err = h_EMu->GetBinError(i);

			double wz  = h_WZ->GetBinContent(i);
			double wz_err  = h_WZ->GetBinError(i);

			double zz  = h_ZZ->GetBinContent(i);
			double zz_err  = h_ZZ->GetBinError(i);
		
			double qcd = h_Dijet->GetBinContent(i);
			double qcd_err = h_Dijet->GetBinError(i);
	
			double wjet = h_Wjet->GetBinContent(i);
			double wjet_err = h_Wjet->GetBinError(i);
		
			double data_bkgsub = data_org - (emu+wz+zz+wjet+qcd);
			double data_bkgsub_err = sqrt((data_org_err*data_org_err) + (emu_err*emu_err) + (wz_err*wz_err) + (zz_err*zz_err) + (wjet_err*wjet_err) + (qcd_err*qcd_err));
	
			Final->SetBinContent(i, data_bkgsub);
			Final->SetBinError(i, data_bkgsub_err);

		}
	}

	//==================== Function to apply Final corrections
	void ApplyFinal_EffCorr(TH1D *unfold, TH1D *effCorr, TH1D *effCorrFinal) {

		for(int i=1; i<44; i++){
			double unfoldData = unfold->GetBinContent(i);
			double effCorrec  = effCorr->GetBinContent(i);
			//cout<<"eff Correction: "<<effCorrec<<endl;

			double unfoldDataErr = unfold->GetBinError(i);
			double effCorrecErr  = effCorr->GetBinError(i);

			//if(unfoldData > 0. && effCorrec > 0.) {
			double FinalCorr = unfoldData/effCorrec;
			double FinalCorrErr = unfoldData/effCorrec * sqrt((unfoldDataErr/unfoldData)*(unfoldDataErr/unfoldData)+(effCorrecErr/effCorrec)*(effCorrecErr/effCorrec));

			effCorrFinal->SetBinContent(i,FinalCorr);
			effCorrFinal->SetBinError(i,FinalCorrErr);

			//}
		}
	}

	void Calc_AccEff(TH1D *AccTotal, TH1D *AccPass, TH1D *EffTotal, TH1D *EffPass, TH1D *AccEff) {

		for(int i=1; i<44; i++){
			double low  = AccPass->GetBinLowEdge(i);
			double high = AccPass->GetBinLowEdge(i+1);

			double numAcc = AccPass->GetBinContent(i);
			double denAcc = AccTotal->GetBinContent(i);
			double denerr = (AccTotal->GetBinContent(i)*AccTotal->GetBinContent(i))/(AccTotal->GetBinError(i)*AccTotal->GetBinError(i));
			double Acc = numAcc/denAcc;
			double errAcc = sqrt(Acc*(1-Acc)/denerr);

			double numEff = EffPass->GetBinContent(i);
			double denEff = EffTotal->GetBinContent(i);
			double denerr1 = (EffTotal->GetBinContent(i)*EffTotal->GetBinContent(i))/(EffTotal->GetBinError(i)*EffTotal->GetBinError(i));
			double Eff = numEff/denEff;
			double errEff = sqrt(Eff*(1-Eff)/denerr1);

			double totalEff = Acc * Eff;
			double totalErr = totalEff * sqrt((errAcc/Acc)*(errAcc/Acc) + (errEff/Eff)*(errEff/Eff));

			AccEff->SetBinContent(i, totalEff);
			AccEff->SetBinError(i, totalErr);

		}
	}

	//==================== Function to calculate diff XSec
	void Obtain_dSigma_dM(TH1D *h, TH1D *h1) {

		for(int i=1; i<44; i++){
			double BinWidth = h->GetBinWidth(i);

			double xSec = h->GetBinContent(i);
			double xSec_dM = xSec/BinWidth;

			double error_before = h->GetBinError(i);
			double error_after = error_before/BinWidth;
	
			h1->SetBinContent(i, xSec_dM);
			h1->SetBinError(i, error_after);
		}
	}


};


class SysUncTools_EffCorr : public SysUncTools_Base
{
public:

	// -- Reco + ID Efficiency & error -- //
	Double_t Eff_Reco_Data[nEtaBin][nPtBin];
	Double_t EffErr_Tot_Reco_Data[nEtaBin][nPtBin];

	// -- ID Efficiency & error -- //
	Double_t Eff_ID_Data[nEtaBin][nPtBin];
	Double_t EffErr_Tot_ID_Data[nEtaBin][nPtBin];
	
	// -- Trigger Efficiency & error -- //
	Double_t Eff_Trigger_Data[nEtaBin][nPtBin];
	Double_t EffErr_Tot_Trigger_Data[nEtaBin][nPtBin];

	// -- Reco Efficiency & error -- //
	Double_t Eff_Reco_MC[nEtaBin][nPtBin];
	Double_t EffErr_Tot_Reco_MC[nEtaBin][nPtBin];

	// -- ID Efficiency & error -- //
	Double_t Eff_ID_MC[nEtaBin][nPtBin];
	Double_t EffErr_Tot_ID_MC[nEtaBin][nPtBin];
	
	// -- Trigger Efficiency & error -- //
	Double_t Eff_Trigger_MC[nEtaBin][nPtBin];
	Double_t EffErr_Tot_Trigger_MC[nEtaBin][nPtBin];

	// -- Reco Efficiency smeared map -- //
	Double_t Eff_Reco_Data_Smeared[nEffMap][nEtaBin][nPtBin];
	Double_t Eff_Reco_MC_Smeared[nEffMap][nEtaBin][nPtBin];
	////////////////////////////////////////////

	// -- Isolation Efficiency smeared map -- //
	Double_t Eff_ID_Data_Smeared[nEffMap][nEtaBin][nPtBin];
	Double_t Eff_ID_MC_Smeared[nEffMap][nEtaBin][nPtBin];
	////////////////////////////////////////////

	// -- Trigger Efficiency smeared map -- //
	Double_t Eff_Trigger_Data_Smeared[nEffMap][nEtaBin][nPtBin];
	Double_t Eff_Trigger_MC_Smeared[nEffMap][nEtaBin][nPtBin];
	
	// -- Differential X-section: Calculated in "CalcXsec_AllMap" Method -- //
	TH1D *h_xSec_dM_FSRCorr_CV;
	h_xSec_dM_FSRCorr_CV->Sumw2();
	TH1D *h_xSec_dM_FSRCorr_Smeared[nEffMap];
	h_xSec_dM_FSRCorr_Smeared[nEffMap]->Sumw2();
	TH1D *h_RelDiff_massBin[44];

	TH1D *h_mass_SFCorr_CV;
	h_mass_SFCorr_CV->Sumw2();
	TH1D *h_mass_SFCorr_Smeared[nEffMap];
	h_mass_SFCorr_Smeared[nEffMap]->Sumw2();

	SysUncTools_EffCorr()
	{
		// -- Setting the differential x-section histograms -- //
		h_xSec_dM_FSRCorr_CV = new TH1D("h_xSec_dM_FSRCorr_CV", "", nMassBin, MassBinEdges);
		for(Int_t i=0; i<nEffMap; i++)
			h_xSec_dM_FSRCorr_Smeared[i] = new TH1D("h_xSec_dM_FSRCorr_Smeared_"+TString::Format("%d", i), "", nMassBin, MassBinEdges);

		// -- Setting the h_RelDiff_massBin -- //
		for(Int_t i=1; i<44; i++)
			h_RelDiff_massBin[i] = new TH1D("h_RelDiff_massBin_"+TString::Format("%d", i), "", 4000, -2, 2);

		this->SetupCentralValueStatError();
		this->MakeSmearedEffMap();
	}



	//======================Calculating Efficiency Maps
	void SetupCentralValueStatError() {

		TFile *file1 = new TFile("RECO_2D.root");
		TFile *file2 = new TFile("Medium_2D.root");
		TFile *file3 = new TFile("Trig_2D.root");

		TH2D *h_Reco_Data    = (TH2D*)file1->Get("hdeff");
		TH2D *h_Reco_MC      = (TH2D*)file1->Get("hmeff");

		TH2D *h_ID_Data   = (TH2D*)file2->Get("hdeff");
		TH2D *h_ID_MC     = (TH2D*)file2->Get("hmeff");

		TH2D *h_Trigger_Data = (TH2D*)file3->Get("hdeff");
		TH2D *h_Trigger_MC   = (TH2D*)file3->Get("hmeff");

		Int_t nEtaBins = h_Reco_Data->GetNbinsX();
		Int_t nPtBins = h_Reco_Data->GetNbinsY();

		for(Int_t iter_x = 0; iter_x < nEtaBins; iter_x++)
		{
			for(Int_t iter_y = 0; iter_y < nPtBins; iter_y++)
			{
				Int_t i_etabin = iter_x + 1;
				Int_t i_ptbin = iter_y + 1;

				// -- Reco cetral value & error -- // 
				Double_t Reco_Data = h_Reco_Data->GetBinContent(i_etabin, i_ptbin);
				Double_t Reco_Data_Error = h_Reco_Data->GetBinError(i_etabin, i_ptbin);

				Eff_Reco_Data[iter_x][iter_y] = Reco_Data;
				EffErr_Tot_Reco_Data[iter_x][iter_y] = Reco_Data_Error;

				Double_t Reco_MC = h_Reco_MC->GetBinContent(i_etabin, i_ptbin);
				Double_t Reco_MC_Error = h_Reco_MC->GetBinError(i_etabin, i_ptbin);

				Eff_Reco_MC[iter_x][iter_y] = Reco_MC;
				EffErr_Tot_Reco_MC[iter_x][iter_y] = Reco_MC_Error;
				//////////////////////////////////////////////////////////////////////////////

				// -- ID central value & error -- //
				Double_t ID_Data = h_ID_Data->GetBinContent(i_etabin, i_ptbin);
				Double_t ID_Data_Error = h_ID_Data->GetBinError(i_etabin, i_ptbin);

				Eff_ID_Data[iter_x][iter_y] = ID_Data;
				EffErr_Tot_ID_Data[iter_x][iter_y] = ID_Data_Error;

				Double_t ID_MC = h_ID_MC->GetBinContent(i_etabin, i_ptbin);
				Double_t ID_MC_Error = h_ID_MC->GetBinError(i_etabin, i_ptbin);

				Eff_ID_MC[iter_x][iter_y] = ID_MC;
				EffErr_Tot_ID_MC[iter_x][iter_y] = ID_MC_Error;
				//////////////////////////////////////////////////////////////////////////////

				// -- Ele23 Trigger central value & error -- //
				Double_t Trigger_Data = h_Trigger_Data->GetBinContent(i_etabin, i_ptbin);
				Double_t Trigger_Data_Error = h_Trigger_Data->GetBinError(i_etabin, i_ptbin);

				Eff_Trigger_Data[iter_x][iter_y] = Trigger_Data;
				EffErr_Tot_Trigger_Data[iter_x][iter_y] = Trigger_Data_Error;
  
				Double_t Trigger_MC = h_Trigger_MC->GetBinContent(i_etabin, i_ptbin);
				Double_t Trigger_MC_Error = h_Trigger_MC->GetBinError(i_etabin, i_ptbin);

				Eff_Trigger_MC[iter_x][iter_y] = Trigger_MC;
				EffErr_Tot_Trigger_MC[iter_x][iter_y] = Trigger_MC_Error;

			}
		}
     
		cout << "========================================================================" << endl;
		cout << "[Setting for efficiency cetral value and total error is completed]" << endl;
		cout << "========================================================================" << endl;
		cout << endl;
	}


	void MakeSmearedEffMap() {
		
		TRandom3 eran;
		eran.SetSeed(0);

		for(Int_t i_map=0; i_map<nEffMap; i_map++)
		{
			for(Int_t iter_x=0; iter_x<nEtaBin; iter_x++)
			{
				for(Int_t iter_y=0; iter_y<nPtBin; iter_y++)
				{
	      
					Eff_Reco_Data_Smeared[i_map][iter_x][iter_y] = Eff_Reco_Data[iter_x][iter_y] + eran.Gaus(0.0, 1.0) * EffErr_Tot_Reco_Data[iter_x][iter_y];
					Eff_Reco_MC_Smeared[i_map][iter_x][iter_y] = Eff_Reco_MC[iter_x][iter_y] + eran.Gaus(0.0, 1.0) * EffErr_Tot_Reco_MC[iter_x][iter_y];

					Eff_ID_Data_Smeared[i_map][iter_x][iter_y] = Eff_ID_Data[iter_x][iter_y] + eran.Gaus(0.0, 1.0) * EffErr_Tot_ID_Data[iter_x][iter_y];
					Eff_ID_MC_Smeared[i_map][iter_x][iter_y] = Eff_ID_MC[iter_x][iter_y] + eran.Gaus(0.0, 1.0) * EffErr_Tot_ID_MC[iter_x][iter_y];

					Eff_Trigger_Data_Smeared[i_map][iter_x][iter_y] = Eff_Trigger_Data[iter_x][iter_y] + eran.Gaus(0.0, 1.0) * EffErr_Tot_Trigger_Data[iter_x][iter_y];
					Eff_Trigger_MC_Smeared[i_map][iter_x][iter_y] = Eff_Trigger_MC[iter_x][iter_y] + eran.Gaus(0.0, 1.0) * EffErr_Tot_Trigger_MC[iter_x][iter_y];


				} // -- end of for(Int_t iter_y=0; iter_y<nPtBin; iter_y++) -- //

			} // -- end of for(Int_t iter_x=0; iter_x<nEtaBin; iter_x++) -- //

		} // -- end of for(Int_t i_map=0; i_map<nEffMap; i_map++) -- //

		cout << "======================================" << endl;
		cout << "[" << nEffMap << " Smeared Maps are produced]" << endl;
		cout << "======================================" << endl;
		cout << endl;

	}

	Double_t EfficiencySF_EventWeight(float pt1, float eta1, float pt2, float eta2) {

		Double_t weight1 = -999;

		// *************************** Electron 1**************************
		Int_t ptbin1 = FindPtBin(pt1);
		Int_t etabin1 = FindEtaBin(eta1);

		Double_t Eff_ele1_Data = Eff_Reco_Data[etabin1][ptbin1] * Eff_ID_Data[etabin1][ptbin1];
		Double_t Eff_ele1_MC = Eff_Reco_MC[etabin1][ptbin1] * Eff_ID_MC[etabin1][ptbin1];

		// *************************** Electron 2**************************
		Int_t ptbin2 = FindPtBin(pt2);
		Int_t etabin2 = FindEtaBin(eta2);

		Double_t Eff_ele2_Data = Eff_Reco_Data[etabin2][ptbin2] * Eff_ID_Data[etabin2][ptbin2];
		Double_t Eff_ele2_MC = Eff_Reco_MC[etabin2][ptbin2] * Eff_ID_MC[etabin2][ptbin2];

		Double_t Eff_EventTrig_Data = 0;
		Double_t Eff_EventTrig_MC = 0;

		Double_t Eff_Trig_ele1_Data = Eff_Trigger_Data[etabin1][ptbin1];
		Double_t Eff_Trig_ele2_Data = Eff_Trigger_Data[etabin2][ptbin2];
		Eff_EventTrig_Data = Eff_Trig_ele1_Data + Eff_Trig_ele2_Data - Eff_Trig_ele1_Data * Eff_Trig_ele2_Data;

		Double_t Eff_Trig_ele1_MC = Eff_Trigger_MC[etabin1][ptbin1];
		Double_t Eff_Trig_ele2_MC = Eff_Trigger_MC[etabin2][ptbin2];
		Eff_EventTrig_MC = Eff_Trig_ele1_MC + Eff_Trig_ele2_MC - Eff_Trig_ele1_MC * Eff_Trig_ele2_MC;

		Double_t Eff_Data_all = Eff_ele1_Data * Eff_ele2_Data * Eff_EventTrig_Data;
		Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC * Eff_EventTrig_MC;

		weight1 = Eff_Data_all / Eff_MC_all;

		return weight1;
	}


	Double_t EfficiencySF_Smeared_EventWeight(Int_t i_map, float pt1, float eta1, float pt2, float eta2) {

		Double_t weight2 = -999;

		// *************************** Electron 1**************************
		Int_t ptbin1 = FindPtBin(pt1);
		Int_t etabin1 = FindEtaBin(eta1);
	
		Double_t Eff_ele1_Data = Eff_Reco_Data_Smeared[i_map][etabin1][ptbin1] * Eff_ID_Data_Smeared[i_map][etabin1][ptbin1];
		Double_t Eff_ele1_MC = Eff_Reco_MC_Smeared[i_map][etabin1][ptbin1] * Eff_ID_MC_Smeared[i_map][etabin1][ptbin1];
	
		// *************************** Electron 2**************************
		Int_t ptbin2 = FindPtBin(pt2);
		Int_t etabin2 = FindEtaBin(eta2);
	
		Double_t Eff_ele2_Data = Eff_Reco_Data_Smeared[i_map][etabin2][ptbin2] * Eff_ID_Data_Smeared[i_map][etabin2][ptbin2];
		Double_t Eff_ele2_MC = Eff_Reco_MC_Smeared[i_map][etabin2][ptbin2] * Eff_ID_MC_Smeared[i_map][etabin2][ptbin2];
	
		Double_t Eff_EventTrig_Data = 0;
		Double_t Eff_EventTrig_MC = 0;
	
		Double_t Eff_Trig_ele1_Data = Eff_Trigger_Data_Smeared[i_map][etabin1][ptbin1];
		Double_t Eff_Trig_ele2_Data = Eff_Trigger_Data_Smeared[i_map][etabin2][ptbin2];
		Eff_EventTrig_Data = 1 + Eff_Trig_ele1_Data + Eff_Trig_ele2_Data - Eff_Trig_ele1_Data * Eff_Trig_ele2_Data;
	
		Double_t Eff_Trig_ele1_MC = Eff_Trigger_MC_Smeared[i_map][etabin1][ptbin1];
		Double_t Eff_Trig_ele2_MC = Eff_Trigger_MC_Smeared[i_map][etabin2][ptbin2];
		Eff_EventTrig_MC = 1 + Eff_Trig_ele1_MC + Eff_Trig_ele2_MC - Eff_Trig_ele1_MC * Eff_Trig_ele2_MC;
	
		Double_t Eff_Data_all = Eff_ele1_Data * Eff_ele2_Data * Eff_EventTrig_Data;
		Double_t Eff_MC_all = Eff_ele1_MC * Eff_ele2_MC * Eff_EventTrig_MC;
	
		// cout << "Eff_Data_all: " << Eff_Data_all << ", Eff_MC_all: " << Eff_MC_all << endl;
		weight2 = Eff_Data_all / Eff_MC_all;

		return weight2;
	}


	void CorrectedEff_AllMap() {

		cout << "===============================================================================" << endl;
		cout << "[Start the calculation of *Corrected* MC-truth efficiency for each Smeared map]" << endl;
		cout << "===============================================================================" << endl;
		cout << endl;
	
		TFile *f = TFile::Open("DY_forEff_M10to3000_13June.root");
	
		double Ele1PT, Ele2PT, Ele1Eta, Ele2Eta, ZMass, lumiWeights, genWeights;
	
		TTree *t;
		t = (TTree*)f->Get("tree");
		Long64_t entries = t->GetEntries();
		cout<<"Total entries: "<<entries<<endl;
		cout<<""<<endl;
	
		t->SetBranchAddress("Ele1PT",&Ele1PT);
		t->SetBranchAddress("Ele1Eta",&Ele1Eta);
		t->SetBranchAddress("Ele2PT",&Ele2PT);
		t->SetBranchAddress("Ele2Eta",&Ele2Eta);
		t->SetBranchAddress("ZMass",&ZMass);
		t->SetBranchAddress("lumiWeights",&lumiWeights);
		t->SetBranchAddress("genWeights",&genWeights);
	
		TH1D *h_mass_EffTotal = new TH1D("h_mass_SFTotal", "", nMassBin, MassBinEdges);
		h_mass_EffTotal->Sumw2();
		TH1D *h_mass_EffPass_Corr_CV = new TH1D("h_mass_EffPass_Corr_CentralValue", "", nMassBin, MassBinEdges);
		h_mass_EffPass_Corr_CV->Sumw2();
	
		TH1D *h_mass_EffPass_Corr[nEffMap];	
		for(Int_t i_map=0; i_map<nEffMap; i_map++) {
			h_mass_EffPass_Corr[i_map] = new TH1D("h_mass_EffPass_Corr_"+TString::Format("%d", i_map), "", nMassBin, 	MassBinEdges);
			h_mass_EffPass_Corr[i_map]->Sumw2();
		}

		Double_t Eff_SF_CV = -999; // -- Efficiency correction factor for the central value -- //
		Double_t Eff_SF[nEffMap] = {-999}; // -- Efficiency correction factor for each smeared map -- //

	    	for(Long64_t j=0; j<10000; j++){
	      		t->GetEntry(j);

			if(j%1000 == 0) cout<<"entry: "<<j<<endl;

			h_mass_EffTotal->Fill(ZMass, lumiWeights * genWeights);
	
			Eff_SF_CV = EfficiencySF_EventWeight(Ele1PT, Ele1Eta, Ele2PT, Ele2Eta);

			h_mass_EffPass_Corr_CV->Fill(ZMass, lumiWeights * genWeights * Eff_SF_CV);

			for(Int_t i_map=0; i_map<nEffMap; i_map++)
			{
				Eff_SF[i_map] = EfficiencySF_Smeared_EventWeight(i_map, Ele1PT, Ele1Eta, Ele2PT, Ele2Eta);
				h_mass_EffPass_Corr[i_map]->Fill(ZMass, lumiWeights * genWeights * Eff_SF[i_map]);
	
			}

		}

		h_mass_SFCorr_CV = (TH1D*) h_mass_EffPass_Corr_CV->Clone();
		h_mass_SFCorr_CV->Divide(h_mass_EffTotal);

		for(Int_t i_map=0; i_map<nEffMap; i_map++)
		{
			h_mass_SFCorr_Smeared[i_map] = (TH1D*)h_mass_EffPass_Corr[i_map]->Clone();
			h_mass_SFCorr_Smeared[i_map]->Divide(h_mass_EffTotal);
		}

	}
	
	Int_t FindPtBin(Double_t Pt) {
	
		const Int_t nPtBins = 5;
		Double_t PtBinEdges[nPtBins+1] = {10, 20, 30, 40, 50, 2000};

		Int_t ptbin = 9999;

		// -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
		if( Pt >= PtBinEdges[nPtBins] )
			ptbin = nPtBins-1;
		else
		{
			for(Int_t i=0; i<nPtBins; i++)
			{
				if( Pt >= PtBinEdges[i] && Pt < PtBinEdges[i+1] )
				{
					ptbin = i;
					break;
				}
			}
		}
	
		return ptbin;
	}

	Int_t FindEtaBin(Double_t eta) {

		const Int_t nEtaBins = 10;
		Double_t EtaBinEdges[nEtaBins+1] = {-2.5, -2, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2, 2.5};

		Int_t etabin = 9999;

		for(Int_t i=0; i<nEtaBins; i++)
		{
			if( eta >= EtaBinEdges[i] && eta < EtaBinEdges[i+1] )
			{
				etabin = i;
				break;
			}
		}
	
		return etabin;
	}


	void CalcXsec_AllMap(TString version) {

		gSystem->Load("libRooUnfold");
	
		//============ Acceptance & Efficiency
		TFile *f1 = TFile::Open("DY_10to3000_29Aug.root");

		TH1D *h_AccPass = (TH1D*)f1->Get("h_mass_AccPass");	TH1D *h_AccTotal = (TH1D*)f1->Get("h_mass_AccTotal");
		TH1D *h_EffTotal = (TH1D*)f1->Get("h_mass_EffTotal");	TH1D *h_EffPass = (TH1D*)f1->Get("h_mass_EffPass");

		//============ Data & Data-Driven background
		TFile *f2 = TFile::Open("SingleElectron_UnScaleCorr_new_case2.root");
		TFile *f3 = TFile::Open("estimated_EEMass_v2.root");
		TFile *f4 = TFile::Open("Est_fmMC_diBoson.root");
		TFile *f5 = TFile::Open("dijet_wjet_BkgFR.root");
  
		TH1D *hData  = (TH1D*)f2->Get("dieleMass");
		TH1D *hEMu   = (TH1D*)f3->Get("estimated_EEMass");
		TH1D *hWZ    = (TH1D*)f4->Get("h_wz_corr_Mass");
		TH1D *hZZ    = (TH1D*)f4->Get("h_zz_corr_Mass");
		TH1D *hDijet = (TH1D*)f5->Get("DiJet_fromData_ControlRegion");
		TH1D *hWjet  = (TH1D*)f5->Get("WJets_fromData_ControlRegion");

		//============ Unfolding
		TFile *f6 = TFile::Open("RespObj_detUnfolding.root");
		RooUnfoldResponse *UnfoldRes1 = (RooUnfoldResponse *)f6->Get("Unfold_DetectorRes1");
		RooUnfoldResponse *UnfoldRes3 = (RooUnfoldResponse *)f6->Get("Unfold_FSRCorr");

		// Histograms declaration 					
		//**********************************************************************************

		TH1D *h_AccEff  = new TH1D("h_AccEff", "", nMassBin, MassBinEdges);

		TH1D *h_yield = new TH1D("h_yield", "", nMassBin, MassBinEdges);
		TH1D *h_yield_Unfolded = new TH1D("h_yield_Unfolded","", nMassBin, MassBinEdges);

		TH1D *h_Efficiency_SFCorr = new TH1D("Efficiency_SFCorr", "", nMassBin, MassBinEdges);
		TH1D *h_AccEffCorr = new TH1D("h_AccEffCorr", "", nMassBin, MassBinEdges);

		TH1D *h_yield_Unfolded_SFCorr = new TH1D("h_yield_Unfolded_SFCorr", "", nMassBin, MassBinEdges);
		TH1D *h_yield_Unfolded_AccEff = new TH1D("h_yield_Unfolded_AccEff", "", nMassBin, MassBinEdges);
		TH1D *h_FSRCorrected = new TH1D("h_FSRCorrected", "", nMassBin, MassBinEdges);


		// Call Functions
		//****************************************************************************************
		Calc_AccEff(h_AccTotal, h_AccPass, h_EffTotal, h_EffPass, h_AccEff); //Computing Acceptance * Eff

		//ObtainYieldHistogram(hData, hEMu, hDijet, hWjet, h_yield); //Background Subtracted Data Histo
		ObtainYieldHistogram(hData, hEMu, hWZ, hZZ, hDijet, hWjet, h_yield);

		Calculation_DiffXsec(h_yield, UnfoldRes1, h_mass_SFCorr_CV, h_AccEff, 
					UnfoldRes3, h_xSec_dM_FSRCorr_CV);

		for(Int_t i_map=0; i_map<nEffMap; i_map++) {
			Calculation_DiffXsec(h_yield, UnfoldRes1, h_mass_SFCorr_Smeared[i_map], h_AccEff, 
					UnfoldRes3, h_xSec_dM_FSRCorr_Smeared[i_map]);


			// -- Calculate the difference with the central value -- //
			TH1D *h_diff_CV_Smeared = (TH1D*)h_xSec_dM_FSRCorr_Smeared[i_map]->Clone();
			
			h_diff_CV_Smeared->Scale( -1 );
			h_diff_CV_Smeared->Add(h_xSec_dM_FSRCorr_CV, 1);


			// -- Insert the difference between central value and smeared value in the histogram -- //
			for(Int_t i=1; i<44; i++)
			{

				Double_t diff = h_diff_CV_Smeared->GetBinContent(i);
				Double_t CentralValue = h_xSec_dM_FSRCorr_CV->GetBinContent(i);
				h_RelDiff_massBin[i]->Fill( diff / CentralValue );

				printf("\t[%d mass bin] RelDiff = %9.6e\n", i, diff / CentralValue);
			}
		}


	} // SysUncEffMaps_v1() ends
};
