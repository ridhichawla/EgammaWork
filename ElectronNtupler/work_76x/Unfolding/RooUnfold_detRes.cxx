//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 348 2014-08-08 22:18:23Z T.J.Adye@rl.ac.uk $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#endif


//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfold_detRes()
{

  TFile *f1 = TFile::Open("/home/ridhi/Work/Analysis/76X/Unfolding/detRes_Unfold/DYEE_final_new_test.root");
  //TFile *f1 = TFile::Open("/home/ridhi/Work/Analysis/76X/Unfolding/detRes_Unfold/DY_M10to3500_new.root");

  TFile *f2 = TFile::Open("/home/ridhi/Work/Analysis/76X/DYSpectrum/MediumId/histo_bkgsub_datadriven.root");

  //TFile *file1 = new TFile ("/home/ridhi/Work/Analysis/76X/Unfolding/detRes_Unfold/RespObj_detRes.root","RECREATE");
  //TFile *file1 = new TFile ("/home/ridhi/Work/Analysis/76X/Unfolding/detRes_Unfold/RespObj_detRes_AltSig2.root","RECREATE");

  //TH2D *h2Resp = (TH2D*)f1->Get("RespMatrix");
  TH2D *h2Resp = (TH2D*)f1->Get("resp_DetRes");
  TH1D *hTrue  = (TH1D*)f1->Get("gen_EEMass");
  TH1D *hMeas  = (TH1D*)f1->Get("reco_EEMass");

  TH1D *hData  = (TH1D*)f2->Get("dieleMass");

  //*********************************************************************************************
  TCanvas *c1 = new TCanvas("c1","",141,211,546,489);
  c1->Draw();
  c1->cd();

  c1_2 = new TPad("c1_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c1_2->SetBottomMargin(0.1);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetLogx();
  c1_2->SetLogy();
  c1_2->SetGridx();
  c1_2->SetGridy();

  Double_t xbins[44] ={15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,
320,380,440,510,600,700,830,1000,1500,3000};

  /*string file;
  std::ofstream outfile;
  outfile.open ("/home/ridhi/Work/Analysis/76X/Unfolding/Inputs1/nEvents_Iter14.txt");
  outfile<<std::fixed;*/

  RooUnfoldResponse response (hMeas, hTrue, h2Resp);
  response.UseOverflow(false);

  cout << "==================================== UNFOLD ===================================" << endl;

  RooUnfoldBayes   unfold (&response, hMeas, 14);    // OR
  //RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
  //RooUnfoldTUnfold unfold (&response, hMeas);

  //hData->Print("all");

  TH1D* hUnfold= (TH1D*) unfold.Hreco();
  cout<<"Fakes : "<<response.FakeEntries()<<endl;
  //cout<<"nEvents unfolded 60-120: "<<hUnfold->Integral(10,22)<<endl;

  int nbins = hUnfold->GetNbinsX();
  cout<<"nbins: "<<nbins<<endl;

  for(int i=1; i<nbins+1; i++) {

   double low = hUnfold->GetBinLowEdge(i);
   double high = hUnfold->GetBinLowEdge(i+1);

   double binContent = hUnfold->GetBinContent(i);
   double binError   = hUnfold->GetBinError(i);
   //cout<<low<<"   "<<high<<endl;
   //cout<<"Unfolded data: "<<binContent<<"   "<<binError<<endl;
   //outfile<<setprecision(2)<<binContent<<endl;
   //cout<<""<<endl;
  }

  hTrue->SetStats(0);
  hUnfold->SetStats(0);

  hTrue->SetLineColor(kOrange-2);
  hTrue->SetFillColor(kOrange-2);

  /*hMeas->SetLineColor(kBlue);
  hTrue->SetLineColor(kRed);
  hTrue->SetMarkerColor(kRed);
  hTrue->SetMarkerStyle(24);
  hTrue->SetMarkerSize(0.71);*/

  hUnfold->SetLineColor(kBlack);
  hUnfold->SetMarkerColor(kBlack);
  hUnfold->SetMarkerStyle(20);
  hUnfold->SetMarkerSize(0.7);

  hData->SetLineColor(kRed);
  hData->SetMarkerColor(kRed);
  hData->SetMarkerStyle(22);
  hData->SetMarkerSize(0.9);

  hTrue->SetTitle("");
  hTrue->GetXaxis()->SetTitle(""); 
  hTrue->GetXaxis()->SetTickLength(0.03);
  hTrue->GetXaxis()->SetTitleOffset(1.05);
  hTrue->GetXaxis()->SetLabelSize(0.03);
  hTrue->GetXaxis()->SetLabelOffset(999);
  hTrue->GetXaxis()->SetRangeUser(15.,3000.);
 
  hTrue->GetYaxis()->SetTitle("Number of Events");
  hTrue->GetYaxis()->SetLabelFont(42);
  hTrue->GetYaxis()->SetLabelSize(0.03);
  hTrue->GetYaxis()->SetTitleSize(0.05);
  hTrue->GetYaxis()->SetTitleOffset(0.90);
  hTrue->GetYaxis()->SetTitleFont(42);
  hTrue->GetYaxis()->SetRangeUser(0.1,1000000.);

  hTrue->Draw("hist][");
  hUnfold->Draw("esame");
  hData->Draw("esame");
  //hMeas->Draw("hist same");
  //hUnfold->Print("all");

  TLegend *legend1=new TLegend(0.485166,0.7123843,0.8847336,0.8630624);
  legend1->SetTextSize(0.04);
  legend1->SetFillColor(0);
  legend1->SetLineColor(1);
  //legend1->SetFillStyle(0);
  legend1->AddEntry(hData,"Raw (Data, PreUnfolded)");
  legend1->AddEntry(hUnfold,"Unfolded (Data)");
  //legend1->AddEntry(hMeas,"Reco Level (Measured)","l");
  legend1->AddEntry(hTrue,"Gen Level (Post FSR)","f");
  legend1->Draw();

  c1->cd();

  TH1F *hRatio = (TH1F *)hUnfold->Clone();
  //TH1F *hRatio = (TH1F *)hData->Clone("hRatio");
  hRatio->Divide(hTrue);

  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerSize(0.6);

  hRatio->SetTitle("  ");  
  hRatio->GetXaxis()->SetTitle("Dielectron Invariant Mass [GeV]");
  hRatio->GetXaxis()->SetLabelFont(42);
  hRatio->GetXaxis()->SetLabelSize(0.13);
  hRatio->GetXaxis()->SetLabelOffset(0.01);
  hRatio->GetXaxis()->SetTitleSize(0.16);
  hRatio->GetXaxis()->SetTitleOffset(0.9);
  hRatio->GetXaxis()->SetRangeUser(15.,3000.);
  hRatio->GetXaxis()->SetMoreLogLabels();
  hRatio->GetXaxis()->SetNoExponent();

  hRatio->GetYaxis()->SetTitle("Unfolded/Gen");
  hRatio->GetYaxis()->SetTitleSize(0.15);
  hRatio->GetYaxis()->SetTitleOffset(0.3);
  hRatio->GetYaxis()->SetLabelFont(42);
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetRangeUser(0.9,1.1);
  hRatio->GetYaxis()->SetNdivisions(5);

  c1_1 = new TPad("c1_1", "newpad",0.008064516,0.0116071,0.9899194,0.2299107);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetLogx();
  c1_1->SetGridx();
  c1_1->SetGridy();
  
  c1_1->Range(-85.9335,-19.83656,785.9335,21.48034);
  c1_1->SetFillColor(0);
  c1_1->SetBorderMode(0);
  c1_1->SetBorderSize(1);
  c1_1->SetTopMargin(0.03067478);
  c1_1->SetBottomMargin(0.3047036);
  c1_1->SetFrameBorderMode(0);
  c1_1->SetFrameBorderMode(0);

  hRatio->Draw("");
  TLine l1(15.,1.0,3000.,1.0);

  //c1->SaveAs("/home/ridhi/Work/Analysis/76X/Plots/13TeV/Unfolding/unfold_Data_AltSig_Madgraph.png");
  //c1->SaveAs("/home/ridhi/Work/Analysis/76X/Plots/13TeV/Unfolding/unfoldMC_closure.png");
  c1->Draw();
  c1->Modified();
  c1->Update();

  //response.Write("UnfoldRes_DetectorRes");
  //file1->Close();

}

#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
