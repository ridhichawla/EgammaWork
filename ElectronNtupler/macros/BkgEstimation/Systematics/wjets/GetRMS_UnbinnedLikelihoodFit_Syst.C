#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"

using namespace RooFit;

#define nMassBin 43
#define nMap 1000

void Set_FitInitValues(Double_t BinCenter, Double_t &RangeMax, Double_t &Sigma_Init);
TTree* makeTTree(TFile *f_input, Int_t i_bin);

void GetRMS_UnbinnedLikelihoodFit_Syst()
{

  TFile *f1 = TFile::Open("ROOTFile_SystUnc_BkgEst.root");
  TFile *f2 = new TFile("RelUnc_BkgEst_SystUnc_UnbinnedFit.root","RECREATE");

  Double_t MassBinEdges[nMassBin+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
    64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
    110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
    200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
    830, 1000, 1500, 3000}; // -- Merging high-mass bins -- //

  TH1D *h_RMS = new TH1D("h_RMS", "", nMassBin, MassBinEdges);
  TH1D *h_mean = new TH1D("h_mean", "", nMassBin, MassBinEdges);

  for(Int_t i=0; i<nMassBin; i++)
  {
    Int_t i_bin = i+1;

    TTree *tree = makeTTree(f1, i_bin);

    Double_t BinCenter = ( MassBinEdges[i] + MassBinEdges[i+1] ) / 2.0;
    Double_t RangeMax = 0.1;
    Double_t Sigma_Init = 0.01;
    Set_FitInitValues(BinCenter, RangeMax, Sigma_Init);

    RooRealVar RelDiff("RelDiff","(#sigma_{Smeared} - #sigma_{CV}) / #sigma_{CV}", (-1)*RangeMax, RangeMax);
    RooDataSet data("data","data", tree, RelDiff);

    // --- Make plot of binned dataset showing Poisson error bars (RooFit default)
    RooPlot* frame = RelDiff.frame( Title(TString::Format("%.0lf < M < %.0lf (%02d bin)", MassBinEdges[i], MassBinEdges[i+1], i_bin)) );

    RooRealVar mean("mean", "mean", 0, -2, 2) ;
    RooRealVar sigma("sigma", "sigma", Sigma_Init, 0.0001, 2);
    RooGaussian gauss("gauss", "gauss", RelDiff, mean, sigma);
    gauss.fitTo(data);

    data.plotOn(frame, Binning(50));
    gauss.plotOn(frame);
    gauss.paramOn(frame,Layout(0.6, 0.9, 0.9));
    frame->getAttText()->SetTextSize(0.02);

    //TString CanvasName = TString::Format("c_RelDiff_Bin%02d", i_bin);
    TString CanvasName = "c_RelDiff_BkgEst_Syst_" + TString::Format("_%d", i_bin) + "bin";
    //if( FileName.Contains("FpoF") )
    //CanvasName.ReplaceAll("RelDiff_", "FpoF_RelDiff_");

    //TCanvas *c = new TCanvas(CanvasName, "", 500, 500); c->cd();
    //frame->Draw();

    //c->SaveAs(CanvasName+".pdf");

    f2->cd();
    //c->Write();

    Double_t RMS = sigma.getVal();
    h_RMS->SetBinContent(i_bin, RMS);
    h_RMS->SetBinError(i_bin, 0); // -- how to get the error of the sigma? -- //

    Double_t Mean = mean.getVal();
    h_mean->SetBinContent(i_bin, Mean);
    h_mean->SetBinError(i_bin, 0);
  }

  f2->cd();
  h_RMS->Write();
  h_mean->Write();

}

void Set_FitInitValues( Double_t BinCenter, Double_t &RangeMax, Double_t &Sigma_Init )
{
  /*if( BinCenter > 1500 )
  {
    RangeMax = 4;
    Sigma_Init = 2;
  }*/
  if( BinCenter > 830 )
  {
    RangeMax = 4;
    Sigma_Init = 0.5;
  }
  else if( BinCenter > 380 )
  {
    RangeMax = 0.5;
    Sigma_Init = 0.02;
  }
  else if( BinCenter > 133 )
  {
    RangeMax = 0.3;
    Sigma_Init = 0.01;
  }
  else
  {
    RangeMax = 0.1;
    Sigma_Init = 0.01;
  }
}

TTree* makeTTree(TFile *f_input, Int_t i_bin) 
{
  TTree* tree = new TTree("tree","tree");

  Double_t* RelDiff = new Double_t;

  tree->Branch("RelDiff", RelDiff, "RelDiff/D");

  TString HistName = "h_xSec_dM_FSRCorr"; // h_DiffXsec
  TString FileName = f_input->GetName();
  if( FileName.Contains("FpoF") )
    HistName = "h_FpoF_DiffXsec";

  f_input->cd();
  TH1D *h_DiffXsec_CV = (TH1D*)f_input->Get(HistName+"_CV")->Clone();
  Double_t value_CV = h_DiffXsec_CV->GetBinContent(i_bin);

  for(Int_t i=0; i<nMap; i++)
  {
    TH1D *h_DiffXsec_Smeared = (TH1D*)f_input->Get( HistName+"_Smeared_"+TString::Format("%d", i) );
    Double_t value_Smeared = h_DiffXsec_Smeared->GetBinContent(i_bin);

    *RelDiff = ( value_CV - value_Smeared ) / value_CV;

    // printf("[\t%d th RelDiff = %lf]\n", i, *RelDiff);

    tree->Fill();
  }

  return tree;
}
