#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TEfficiency.h>

#include <vector>

#include <SysUnc_EffMap.h>

void SysUnc_EffCorr()
{

	TFile *f_output = new TFile("ROOTFile_SysUnc_EffCorr.root", "RECREATE");

	SysUncTools_EffCorr *tools = new SysUncTools_EffCorr();

	tools->CorrectedEff_AllMap();
	f_output->cd();

	tools->h_mass_SFCorr_CV->Write();
	for(Int_t i=0; i<nEffMap; i++)
		tools->h_mass_SFCorr_Smeared[i]->Write();

	tools->CalcXsec_AllMap("v20161007_UpdateTnPResults");
	f_output->cd();

	tools->h_xSec_dM_FSRCorr_CV->Write();
	for(Int_t i=0; i<nEffMap; i++)
		tools->h_xSec_dM_FSRCorr_Smeared[i]->Write();

	for(Int_t i=1; i<44; i++)
		tools->h_RelDiff_massBin[i]->Write();

}
