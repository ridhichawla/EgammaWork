#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TMath.h>
#include <TText.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TColor.h>
#include <TLatex.h>
#include <TEfficiency.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;

void setDataHist(TH1D* hist);
void setMCHist(TH1D* hist, const int& color);
void setFRHist(TH1D* hist);
TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator);
TH1D* FRBytRatio(TH1D** numerator, TH1D** denominator);

void estimateFR() {

  const int ptbinnum_endcap = 7;
  double ptbin_endcap[ptbinnum_endcap+1] = {10,25,35,45,55,100,200,1000};
  const int ptbinnum = 7;
  double ptbin[ptbinnum+1] = {10,25,35,45,55,100,200,1000};
  //const int etabinnum = 2;
  //double etabin[etabinnum+1] = {0,1.2,2.4};

  TFile* f[11];
  f[0]  = new TFile("histograms/hist0.root","READ");
  f[3]  = new TFile("histograms/hist3.root","READ");
  f[4]  = new TFile("histograms/hist4.root","READ");
  f[5]  = new TFile("histograms/hist5.root","READ");
  f[6]  = new TFile("histograms/hist6.root","READ");
  f[7]  = new TFile("histograms/hist7.root","READ");
  f[8]  = new TFile("histograms/hist8.root","READ");
  f[9]  = new TFile("histograms/hist9.root","READ");
  f[10] = new TFile("histograms/hist10.root","READ");

  TH1D* denominator_pt_fit_barrel[11];
  TH1D* denominator_pt_fit_endcap[11];
  TH1D* denominator_pt_xsec_barrel[11];
  TH1D* denominator_pt_xsec_endcap[11];

  TH1D* numerator_pt_barrel[11];
  TH1D* numerator_pt_endcap[11];

  TH1D* denominator_barrel[11];
  TH1D* denominator_endcap[11];

  //cout<<"1"<<endl;
  for(int i=0;i<11;i++) {
    if(i==1||i==2) continue;

    denominator_pt_fit_barrel[i] = (TH1D*)f[i]->Get("denominator_pt_barrel")->Clone("denominator_pt_fit_barrel"+TString::Itoa(i,10));
    denominator_pt_xsec_barrel[i] = (TH1D*)f[i]->Get("denominator_pt_barrel")->Clone("denominator_pt_xsec_barrel"+TString::Itoa(i,10));
    denominator_barrel[i] = (TH1D*)f[i]->Get("denominator_barrel")->Clone("denominator_barrel"+TString::Itoa(i,10));
    numerator_pt_barrel[i] = (TH1D*)f[i]->Get("numerator_pt_barrel")->Clone("numerator_pt_barrel"+TString::Itoa(i,10));
    denominator_pt_fit_endcap[i] = (TH1D*)f[i]->Get("denominator_pt_endcap")->Clone("denominator_pt_fit_endcap"+TString::Itoa(i,10));
    denominator_pt_xsec_endcap[i] = (TH1D*)f[i]->Get("denominator_pt_endcap")->Clone("denominator_pt_xsec_endcap"+TString::Itoa(i,10));
    denominator_endcap[i] = (TH1D*)f[i]->Get("denominator_endcap")->Clone("denominator_endcap"+TString::Itoa(i,10));
    numerator_pt_endcap[i] = (TH1D*)f[i]->Get("numerator_pt_endcap")->Clone("numerator_pt_endcap"+TString::Itoa(i,10));

    if(i==6) {
      setDataHist( denominator_pt_fit_barrel[i] );
      setDataHist( denominator_pt_xsec_barrel[i] );
      setDataHist( denominator_barrel[i] );
      setDataHist( numerator_pt_barrel[i] );
      setDataHist( denominator_pt_fit_endcap[i] );
      setDataHist( denominator_pt_xsec_endcap[i] );
      setDataHist( denominator_endcap[i] );
      setDataHist( numerator_pt_endcap[i] );
    }
    else {
      setMCHist( denominator_pt_fit_barrel[i], i );
      setMCHist( denominator_pt_xsec_barrel[i], i );
      setMCHist( denominator_barrel[i], i );
      setMCHist( numerator_pt_barrel[i], i );
      setMCHist( denominator_pt_fit_endcap[i], i );
      setMCHist( denominator_pt_xsec_endcap[i], i );
      setMCHist( denominator_endcap[i], i );
      setMCHist( numerator_pt_endcap[i], i );
    }

  }

  //DYMuMu, ttbar, WJets, WW, tautau, QCD
  double nEvts[11] = {8.1236e+07, 0, 0, 85849572, 3.73193e+12, 2.23076e+06, 0, 988416, 999996, 985598, 54136};
  double xsec[11] = {8.1236e+07, 0, 0, 831.76, 61526.7, 2.23076e+06, 0, 118.7, 47.13, 16.523, 365896};
  double norm_xsec[11];
  double norm_fit_barrel[11];
  double norm_fit_endcap[11];
  double lumi = 2258.066;

  for(int i=0;i<11;i++) {
    norm_xsec[i] = lumi*xsec[i]/nEvts[i];
  }

  // GJets is set to float within 5% //
  /*norm_fit_barrel[0]  = 5.3644e+06/denominator_barrel[0]->Integral();
  norm_fit_barrel[3]  = 4.6221e+05/denominator_barrel[3]->Integral();
  norm_fit_barrel[4]  = 3.0477e+07/denominator_barrel[4]->Integral();
  norm_fit_barrel[5]  = 2.5396e+07/denominator_barrel[5]->Integral();
  norm_fit_barrel[7]  = 4.1481e+04/denominator_barrel[7]->Integral();
  norm_fit_barrel[8]  = 1.4544e+04/denominator_barrel[8]->Integral();
  norm_fit_barrel[9]  = 4.5428e+03/denominator_barrel[9]->Integral();
  norm_fit_barrel[10] = 1.7324e+06/denominator_barrel[10]->Integral();

  norm_fit_endcap[0]  = 1.4493e+06/denominator_endcap[0]->Integral();
  norm_fit_endcap[3]  = 5.7606e+04/denominator_endcap[3]->Integral();
  norm_fit_endcap[4]  = 3.0460e+06/denominator_endcap[4]->Integral();
  norm_fit_endcap[5]  = 5.8361e+07/denominator_endcap[5]->Integral();
  norm_fit_endcap[7]  = 9.6170e+03/denominator_endcap[7]->Integral();
  norm_fit_endcap[8]  = 3.1885e+03/denominator_endcap[8]->Integral();
  norm_fit_endcap[9]  = 1.0079e+03/denominator_endcap[9]->Integral();
  norm_fit_endcap[10] = 2.2278e+06/denominator_endcap[10]->Integral();*/

  // GJets is set to float within 50% //
  norm_fit_barrel[0]  = 5.1504e+06/denominator_barrel[0]->Integral();
  norm_fit_barrel[3]  = 4.4012e+05/denominator_barrel[3]->Integral();
  norm_fit_barrel[4]  = 3.1387e+07/denominator_barrel[4]->Integral();
  norm_fit_barrel[5]  = 2.5415e+07/denominator_barrel[5]->Integral();
  norm_fit_barrel[7]  = 3.9496e+04/denominator_barrel[7]->Integral();
  norm_fit_barrel[8]  = 1.3851e+04/denominator_barrel[8]->Integral();
  norm_fit_barrel[9]  = 4.2715e+03/denominator_barrel[9]->Integral();
  norm_fit_barrel[10] = 1.0030e+06/denominator_barrel[10]->Integral();

  norm_fit_endcap[0]  = 1.4493e+06/denominator_endcap[0]->Integral();
  norm_fit_endcap[3]  = 6.3669e+04/denominator_endcap[3]->Integral();
  norm_fit_endcap[4]  = 3.0460e+06/denominator_endcap[4]->Integral();
  norm_fit_endcap[5]  = 5.9362e+07/denominator_endcap[5]->Integral();
  norm_fit_endcap[7]  = 9.6170e+03/denominator_endcap[7]->Integral();
  norm_fit_endcap[8]  = 3.1885e+03/denominator_endcap[8]->Integral();
  norm_fit_endcap[9]  = 1.0079e+03/denominator_endcap[9]->Integral();
  norm_fit_endcap[10] = 1.1725e+06/denominator_endcap[10]->Integral();

  for(int i=0; i<11; i++) {
    if(i==1||i==2||i==6) continue;

    cout<<"i = "<<i<<endl;
    cout<<"norm_fit_barrel["<<i<<"]"<<"   "<<norm_fit_barrel[i]<<endl;
    cout<<"norm_fit_endcap["<<i<<"]"<<"   "<<norm_fit_endcap[i]<<endl;
    cout<<""<<endl;

    denominator_barrel[i]->Scale(norm_fit_barrel[i]);
    denominator_pt_fit_barrel[i]->Scale(norm_fit_barrel[i]);
    denominator_pt_xsec_barrel[i]->Scale(norm_xsec[i]);
    numerator_pt_barrel[i]->Scale(norm_xsec[i]);

    denominator_endcap[i]->Scale(norm_fit_endcap[i]);
    denominator_pt_fit_endcap[i]->Scale(norm_fit_endcap[i]);
    denominator_pt_xsec_endcap[i]->Scale(norm_xsec[i]);
    numerator_pt_endcap[i]->Scale(norm_xsec[i]);

  }

  TH1D* FR_template_barrel = (TH1D*)FRByTemplate(numerator_pt_barrel, denominator_pt_fit_barrel);
  TH1D* FR_template_endcap = (TH1D*)FRByTemplate(numerator_pt_endcap, denominator_pt_fit_endcap);

  TH1D* FR_xsec_barrel = (TH1D*)FRBytRatio(numerator_pt_barrel, denominator_pt_xsec_barrel);
  TH1D* FR_xsec_endcap = (TH1D*)FRBytRatio(numerator_pt_endcap, denominator_pt_xsec_endcap);

  int W = 1200;
  int H = 1200;

  int H_ref = 1200;
  int W_ref = 1200;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  /*lumi_13TeV = "2759 pb^{-1}";
    writeExtraText = true;
    extraText = "Preliminary";
    drawLogo = false;
    tdrGrid(true);
    lumiTextSize = 0.5;
    cmsTextSize = 0.75;*/

  TH1D* ptFrame = new TH1D("ptFrame","",4,10,1000);
  ptFrame->SetStats(kFALSE);
  ptFrame->GetXaxis()->SetTitle("p_{T}[GeV]");
  ptFrame->GetYaxis()->SetTitle("Fake Rate");

  ptFrame->SetMinimum(0.);
  ptFrame->SetMaximum(1.0);
  ptFrame->GetXaxis()->SetTitleOffset(1);
  ptFrame->GetYaxis()->SetTitleOffset(1.1);
  ptFrame->GetXaxis()->SetTitleSize(0.05);
  ptFrame->GetYaxis()->SetTitleSize(0.05);  
  ptFrame->GetXaxis()->SetLabelSize(0.035);
  ptFrame->GetYaxis()->SetLabelSize(0.035); 
  ptFrame->GetXaxis()->SetMoreLogLabels(); 

  TCanvas* canv = new TCanvas("canv","",600,600);
  canv->SetFillColor(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );

  FR_template_barrel->SetMarkerSize(1.5);
  FR_template_endcap->SetMarkerSize(1.5);
  FR_xsec_barrel->SetMarkerSize(1.5);
  FR_xsec_endcap->SetMarkerSize(1.5);

  FR_template_barrel->SetLineColor(1);
  FR_template_barrel->SetLineWidth(1);
  FR_template_barrel->SetMarkerStyle(20);
  FR_template_barrel->SetMarkerColor(1);

  FR_xsec_barrel->SetLineColor(2);
  FR_xsec_barrel->SetLineWidth(1);
  FR_xsec_barrel->SetMarkerStyle(21);
  FR_xsec_barrel->SetMarkerColor(2);

  FR_template_endcap->SetLineColor(1);
  FR_template_endcap->SetLineWidth(1);
  FR_template_endcap->SetMarkerStyle(20);
  FR_template_endcap->SetMarkerColor(1);

  FR_xsec_endcap->SetLineColor(2);
  FR_xsec_endcap->SetLineWidth(1);
  FR_xsec_endcap->SetMarkerStyle(21);
  FR_xsec_endcap->SetMarkerColor(2);

  canv->cd();
  canv->SetLogx();

  TLegend* legend2 = new TLegend(.45,.65,.75,.89);
  legend2->AddEntry(FR_template_barrel,"Template fitting");
  legend2->AddEntry(FR_xsec_barrel,"Ratio method");
  legend2->SetBorderSize(0);

  ptFrame->Draw();
  //CMS_lumi(canv,4,11);
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  FR_template_barrel->Draw("EPSAME");
  FR_xsec_barrel->Draw("EPSAME");
  legend2->Draw("SAME");
  canv->Print("print/FR_Barrel.png");
  canv->Print("print/FR_Barrel.pdf");

  canv->Clear();
  ptFrame->Draw();
  //CMS_lumi(canv,4,11);
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  FR_template_endcap->Draw("EPSAME");
  FR_xsec_endcap->Draw("EPSAME");
  legend2->Draw("SAME");
  canv->Print("print/FR_Endcap.png");
  canv->Print("print/FR_Endcap.pdf");

  TFile* g = new TFile("result/fakerate.root","RECREATE");
  FR_template_barrel->Write();
  FR_template_endcap->Write();
  FR_xsec_barrel->Write();
  FR_xsec_endcap->Write();
  g->Close();

  cout<<"Template Barrel"<<endl;
  for(int i=1; i<ptbinnum+1; i++) {
    if(i!=1) cout<<",";
    cout<<FR_template_barrel->GetBinContent(i);
  }
  cout<<endl;
  cout<<"Template Endcap"<<endl;
  for(int i=1; i<ptbinnum_endcap+1; i++) {
    if(i!=1) cout<<",";
    cout<<FR_template_endcap->GetBinContent(i);
  }
  cout<<endl;
  cout<<"Ratio Barrel"<<endl;
  for(int i=1; i<ptbinnum+1; i++) {
    if(i!=1) cout<<",";
    cout<<FR_xsec_barrel->GetBinContent(i);
  }
  cout<<endl;
  cout<<"Ratio Endcap"<<endl;
  for(int i=1; i<ptbinnum_endcap+1; i++) {
    if(i!=1) cout<<",";
    cout<<FR_xsec_endcap->GetBinContent(i);
  }
  cout<<endl;

}

void setDataHist(TH1D* hist) {
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(33);
  hist->SetMarkerSize(3);
  hist->SetStats(kFALSE);
  //hist->Sumw2();
}

void setMCHist(TH1D* hist, const int& color) {
  hist->SetFillColor(color+2);
  hist->SetStats(kFALSE);
  //hist->Sumw2();
}

TH1D* FRByTemplate(TH1D** numerator, TH1D** denominator) {

  TString name = ( ((TString)(denominator[5]->GetName())).Contains("barrel") ) ? "FR_template_barrel" : "FR_template_endcap";

  TH1D* num = (TH1D*)numerator[5]->Clone(name);
  TH1D* den = (TH1D*)numerator[0]->Clone(name+"_");

  num->Multiply(numerator[6]);

  den->Add(numerator[3]);
  den->Add(numerator[4]);
  den->Add(numerator[5]);
  den->Add(numerator[7]);
  den->Add(numerator[8]);
  den->Add(numerator[9]);
  den->Add(numerator[10]);
  den->Multiply(denominator[5]);

  num->Divide(den);

  delete den;
  return num;
}

TH1D* FRBytRatio(TH1D** numerator, TH1D** denominator) {

  TString name = ( ((TString)(denominator[5]->GetName())).Contains("barrel") ) ? "FR_xsec_barrel" : "FR_xsec_endcap";

  TH1D* num = (TH1D*)denominator[0]->Clone(name);
  TH1D* den = (TH1D*)numerator[0]->Clone(name+"_");

  num->Add(denominator[3]);
  num->Add(denominator[4]);
  num->Add(denominator[5]);
  num->Add(denominator[7]);
  num->Add(denominator[8]);
  num->Add(denominator[9]);
  num->Add(denominator[10]);
  num->Multiply(numerator[5]);
  num->Multiply(numerator[6]);

  den->Add(numerator[3]);
  den->Add(numerator[4]);
  den->Add(numerator[5]);
  den->Add(numerator[7]);
  den->Add(numerator[8]);
  den->Add(numerator[9]);
  den->Add(numerator[10]);
  den->Multiply(denominator[5]);
  den->Multiply(denominator[6]);

  num->Divide(den);

  delete den;
  return num;
}
