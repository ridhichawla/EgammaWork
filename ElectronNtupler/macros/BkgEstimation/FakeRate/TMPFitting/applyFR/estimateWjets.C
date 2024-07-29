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

void estimateWjets() {

  int W = 1200;
  int H = 1200;

  int H_ref = 1200;
  int W_ref = 1200;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  // UPDATED IN 2017
  //lumi_13TeV = "2759 pb^{-1}";
  //lumiTextSize = 0.5;
  //writeExtraText = true;
  //extraText = "Preliminary";
  //drawLogo = false;

  const int binnum = 43;
  double bins[binnum+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1500, 3000};

  TH1D* massFrame = new TH1D("massFrame","",38,15,3000);
  massFrame->SetMinimum(0.001);
  massFrame->SetMaximum(1000000);
  massFrame->SetStats(kFALSE);
  massFrame->GetXaxis()->SetTitle("Mass[GeV]");
  //massFrame->GetXaxis()->CenterTitle(kTRUE);
  //massFrame->GetYaxis()->CenterTitle(kTRUE);
  massFrame->GetYaxis()->SetTitleOffset(1);
  massFrame->GetYaxis()->SetTitle("Number of events");
  massFrame->GetXaxis()->SetTitleSize(0);
  massFrame->GetYaxis()->SetTitleSize(0.05);
  massFrame->GetXaxis()->SetLabelSize(0);
  //massFrame->GetYaxis()->SetLabelSize(0.025); 
  massFrame->GetXaxis()->SetMoreLogLabels();

  TFile* f[14];
  f[0] = new TFile("histograms/nocorr/fake2.root","READ");
  f[3] = new TFile("histograms/nocorr/fake1.root","READ");
  f[1] = new TFile("histograms/nocorr/fake3.root","READ");
  f[6] = new TFile("histograms/nocorr/fake6.root","READ");
  
  //f[7] = new TFile("result/nocorr/dijetcorr.root","READ");

  TH1D* dijet_template[14];
  TH1D* dijetSS_template[14];
  TH1D* dijet_ratio[14];
  TH1D* dijetSS_ratio[14];
  TH1D* wjets_template[14];
  TH1D* wjetsSS_template[14];
  TH1D* wjets_ratio[14];
  TH1D* wjetsSS_ratio[14];

  TFile* gg = new TFile("result/nocorr/wjets.root","RECREATE");
  
  dijet_template[5] = (TH1D*)f[6]->Get("histDijet1")->Clone();
  //dijet_template[5]->SetFillColor(7);
  dijet_template[5]->SetStats(kFALSE);
  //dijet_template[5]->Sumw2();

  dijetSS_template[5] = (TH1D*)f[6]->Get("histSameDijet1")->Clone();
  //dijetSS_template[5]->SetLineColor(2);
  //dijetSS_template[5]->SetMarkerColor(2);
  //dijetSS_template[5]->SetMarkerStyle(22);
  //dijetSS_template[5]->SetMarkerSize(3);
  dijetSS_template[5]->SetStats(kFALSE);
  //dijetSS_template[5]->Sumw2();

  dijet_template[0] = (TH1D*)f[0]->Get("histDijet1")->Clone();
  dijet_template[3] = (TH1D*)f[3]->Get("histDijet1")->Clone();
  //dijet_template[0]->SetFillColor(2);
  dijet_template[0]->SetStats(kFALSE);
  //dijet_template[0]->Sumw2();

  dijet_template[1] = (TH1D*)f[1]->Get("histDijet1")->Clone();
  //dijet_template[1]->SetFillColor(3);
  dijet_template[1]->SetStats(kFALSE);
  //dijet_template[1]->Sumw2();

  dijet_ratio[5] = (TH1D*)f[6]->Get("histDijet2")->Clone();
  //dijet_ratio[5]->SetFillColor(7);
  dijet_ratio[5]->SetStats(kFALSE);
  //dijet_ratio[5]->Sumw2();

  dijetSS_ratio[5] = (TH1D*)f[6]->Get("histSameDijet2")->Clone();
  //dijetSS_ratio[5]->SetLineColor(2);
  //dijetSS_ratio[5]->SetMarkerColor(2);
  //dijetSS_ratio[5]->SetMarkerStyle(22);
  //dijetSS_ratio[5]->SetMarkerSize(3);
  dijetSS_ratio[5]->SetStats(kFALSE);
  //dijetSS_ratio[5]->Sumw2();

  dijet_ratio[0] = (TH1D*)f[0]->Get("histDijet2")->Clone();
  dijet_ratio[3] = (TH1D*)f[3]->Get("histDijet2")->Clone();
  //dijet_ratio[0]->SetFillColor(2);
  dijet_ratio[0]->SetStats(kFALSE);
  //dijet_ratio[0]->Sumw2();

  dijet_ratio[1] = (TH1D*)f[1]->Get("histDijet2")->Clone();
  //dijet_ratio[1]->SetFillColor(3);
  dijet_ratio[1]->SetStats(kFALSE);
  //dijet_ratio[1]->Sumw2();

  wjets_template[5] = (TH1D*)f[6]->Get("histWJets1")->Clone();
  //wjets_template[5]->SetFillColor(9);
  wjets_template[5]->SetStats(kFALSE);
  //wjets_template[5]->Sumw2();

  wjetsSS_template[5] = (TH1D*)f[6]->Get("histSameWJets1")->Clone();
  wjetsSS_template[1] = (TH1D*)f[1]->Get("histSameWJets1")->Clone();
  wjetsSS_template[0] = (TH1D*)f[0]->Get("histSameWJets1")->Clone();
  wjetsSS_template[3] = (TH1D*)f[3]->Get("histSameWJets1")->Clone();
  //wjetsSS_template[5]->SetFillColor(9);
  //wjetsSS_template[5]->SetLineColor(9);
  //wjetsSS_template[5]->SetMarkerColor(2);
  //wjetsSS_template[5]->SetMarkerStyle(22);
  //wjetsSS_template[5]->SetMarkerSize(4);
  wjetsSS_template[5]->SetStats(kFALSE);
  //wjetsSS_template[5]->Sumw2();

  wjets_template[0] = (TH1D*)f[0]->Get("histWJets1")->Clone();
  wjets_template[3] = (TH1D*)f[3]->Get("histWJets1")->Clone();
  //wjets_template[0]->SetFillColor(2);
  wjets_template[0]->SetStats(kFALSE);
  //wjets_template[0]->Sumw2();

  wjets_template[1] = (TH1D*)f[1]->Get("histWJets1")->Clone();
  //wjets_template[1]->SetFillColor(3);
  wjets_template[1]->SetStats(kFALSE);
  //wjets_template[1]->Sumw2();

  wjets_ratio[5] = (TH1D*)f[6]->Get("histWJets2")->Clone();
  //wjets_ratio[5]->SetFillColor(9);
  wjets_ratio[5]->SetStats(kFALSE);
  //wjets_ratio[5]->Sumw2();

  wjetsSS_ratio[5] = (TH1D*)f[6]->Get("histSameWJets2")->Clone();
  wjetsSS_ratio[1] = (TH1D*)f[1]->Get("histSameWJets2")->Clone();
  wjetsSS_ratio[0] = (TH1D*)f[0]->Get("histSameWJets2")->Clone();
  wjetsSS_ratio[3] = (TH1D*)f[3]->Get("histSameWJets2")->Clone();
  //wjetsSS_ratio[5]->SetFillColor(9);
  //wjetsSS_ratio[5]->SetLineColor(9);
  //wjetsSS_ratio[5]->SetMarkerColor(2);
  //wjetsSS_ratio[5]->SetMarkerStyle(22);
  //wjetsSS_ratio[5]->SetMarkerSize(4);
  wjetsSS_ratio[5]->SetStats(kFALSE);
  //wjetsSS_ratio[5]->Sumw2();

  wjets_ratio[0] = (TH1D*)f[0]->Get("histWJets2")->Clone();
  wjets_ratio[3] = (TH1D*)f[3]->Get("histWJets2")->Clone();
  //wjets_ratio[0]->SetFillColor(2);
  wjets_ratio[0]->SetStats(kFALSE);
  //wjets_ratio[0]->Sumw2();

  wjets_ratio[1] = (TH1D*)f[1]->Get("histWJets2")->Clone();
  //wjets_ratio[1]->SetFillColor(3);
  wjets_ratio[1]->SetStats(kFALSE);
  //wjets_ratio[1]->Sumw2();

  //cout<<"2"<<endl;

  //DYMuMu, ttbar, WJets, WW, tautau, QCD
  double nEvts[14] = {4.5151e+11, 85849572, 3.73193e+12, 2.30899e+12, 3268361, 2.23076e+06, 0 , 0, 0, 0, 0, 988416, 999996,985598};
  double xsec[14] = {2008.4*3, 831.76, 61526.7, 18610, 1915, 2.23076e+06, 0, 0, 0, 0, 0, 118.7, 47.13, 16.523};
  double norm[14];
  double norm1[14];
  double norm2[14];
  // UPDATED IN 2017
  double lumi = 2258.066;
  // double lumi = 2832.673;

  for(int i=0;i<14;i++) {
    if(i>5&&i<11) continue;
    norm[i] = (xsec[i]*lumi)/nEvts[i];
    //cout<<norm[i]<<endl;
  }

  double n_DYJets = 3.6456e+04;
  double n_QCD    = 9.4631e+02;
  double n_WJets  = 6.6929e+02;
  double n_ttbar  = 4.3185e+02;

  double nn_DYJets = 3.0771e+04;
  double nn_QCD    = 4.8031e+02;
  double nn_WJets  = 6.7655e+02;
  double nn_ttbar  = 5.5031e+02;

  dijet_template[0]->Scale(norm[0]);
  dijet_template[3]->Scale(norm[3]);
  dijet_template[0]->Add(dijet_template[3]);
  dijet_ratio[0]->Scale(norm[0]);
  dijet_ratio[3]->Scale(norm[3]);
  dijet_ratio[0]->Add(dijet_ratio[3]);
  dijet_template[1]->Scale(norm[1]);
  dijet_ratio[1]->Scale(norm[1]);

  wjets_template[1]->Scale(norm[1]);
  wjets_ratio[1]->Scale(norm[1]);
  wjets_template[0]->Scale(norm[0]);
  wjets_template[3]->Scale(norm[3]);
  wjets_template[0]->Add(wjets_template[3]);
  wjets_ratio[0]->Scale(norm[0]);
  wjets_ratio[3]->Scale(norm[3]);
  wjets_ratio[0]->Add(wjets_ratio[3]);

  wjetsSS_template[1]->Scale(norm[1]);
  wjetsSS_ratio[1]->Scale(norm[1]);
  wjetsSS_template[0]->Scale(norm[0]);
  wjetsSS_template[3]->Scale(norm[3]);
  wjetsSS_template[0]->Add(wjetsSS_template[3]);
  wjetsSS_ratio[0]->Scale(norm[0]);
  wjetsSS_ratio[3]->Scale(norm[3]);
  wjetsSS_ratio[0]->Add(wjetsSS_ratio[3]);

  dijet_template[5]->Add(dijet_template[0],-1.0);
  dijet_template[5]->Add(dijet_template[1],-1.0);
  dijet_ratio[5]->Add(dijet_ratio[0],-1.0);
  dijet_ratio[5]->Add(dijet_ratio[1],-1.0);

  wjetsSS_template[5]->Write("dataSS_template");
  wjetsSS_ratio[5]->Write("dataSS_ratio");

  cout<<"data(template): "<<wjetsSS_template[5]->Integral()<<endl;
  cout<<"data(ratio): "<<wjetsSS_ratio[5]->Integral()<<endl;

  wjetsSS_template[5]->Add(dijetSS_template[5],-2.0);
  wjetsSS_template[5]->Add(wjetsSS_template[1],-1.0);
  wjetsSS_template[5]->Add(wjetsSS_template[0],-1.0);
  wjetsSS_ratio[5]->Add(dijetSS_ratio[5],-2.0);
  wjetsSS_ratio[5]->Add(wjetsSS_ratio[1],-1.0);
  wjetsSS_ratio[5]->Add(wjetsSS_ratio[0],-1.0);

  cout<<"dijet(template): "<<dijetSS_template[5]->Integral()<<endl;
  cout<<"dijet(ratio): "<<dijetSS_ratio[5]->Integral()<<endl;
  cout<<"ttbar(template): "<<wjetsSS_template[1]->Integral()<<endl;
  cout<<"ttbar(ratio): "<<wjetsSS_ratio[1]->Integral()<<endl;
  cout<<"DY(template): "<<wjetsSS_template[0]->Integral()<<endl;
  cout<<"DY(ratio): "<<wjetsSS_ratio[0]->Integral()<<endl;

  cout<<"Wjets(template) before: "<<wjetsSS_template[5]->Integral()<<endl;
  cout<<"Wjets(ratio) before: "<<wjetsSS_ratio[5]->Integral()<<endl;

  wjetsSS_template[0]->Write("dySS_template");
  wjetsSS_template[1]->Write("ttbarSS_template");
  dijetSS_template[5]->Write("dijetSS_template");
  wjetsSS_template[5]->Write("wjetsSS_template");

  wjetsSS_ratio[0]->Write("dySS_ratio");
  wjetsSS_ratio[1]->Write("ttbarSS_ratio");
  dijetSS_ratio[5]->Write("dijetSS_ratio");
  wjetsSS_ratio[5]->Write("wjetsSS_ratio");

  //wjets_template[1]->Scale(1/norm[1]);
  //wjets_ratio[1]->Scale(1/norm[1]);

  norm1[0] = n_DYJets/wjets_template[0]->Integral();
  norm1[1] = n_ttbar/wjets_template[1]->Integral();
  norm1[2] = n_WJets/wjetsSS_template[5]->Integral();
  norm1[5] = n_QCD/dijet_template[5]->Integral();
  norm2[0] = nn_DYJets/wjets_ratio[0]->Integral();
  norm2[1] = nn_ttbar/wjets_ratio[1]->Integral();
  norm2[2] = nn_WJets/wjetsSS_ratio[5]->Integral();
  norm2[5] = nn_QCD/dijet_ratio[5]->Integral();

  wjets_template[6] = (TH1D*)wjets_template[5]->Clone();
  wjets_ratio[6] = (TH1D*)wjets_ratio[5]->Clone();

  wjets_template[6]->Write("data_template");
  wjets_ratio[6]->Write("data_ratio");

  cout<<"data(template): "<<wjets_template[6]->Integral()<<endl;
  cout<<"data(ratio): "<<wjets_ratio[6]->Integral()<<endl;

  wjets_template[5] = (TH1D*)wjetsSS_template[5]->Clone(); 
  wjets_ratio[5] = (TH1D*)wjetsSS_ratio[5]->Clone(); 

  wjets_template[0]->Scale(norm1[0]);
  wjets_template[1]->Scale(norm1[1]);
  dijet_template[5]->Scale(norm1[5]);
  wjets_template[5]->Scale(norm1[2]);
  wjets_ratio[0]->Scale(norm2[0]);
  wjets_ratio[1]->Scale(norm2[1]);
  dijet_ratio[5]->Scale(norm2[5]);
  wjets_ratio[5]->Scale(norm2[2]);

  wjets_template[0]->Write("dy_template");
  wjets_template[1]->Write("ttbar_template");
  dijet_template[5]->Write("dijet_template");
  wjets_template[5]->Write("wjets_template");

  wjets_ratio[0]->Write("dy_ratio");
  wjets_ratio[1]->Write("ttbar_ratio");
  dijet_ratio[5]->Write("dijet_ratio");
  wjets_ratio[5]->Write("wjets_ratio");

  cout<<"dijet(template): "<<dijet_template[5]->Integral()<<endl;
  cout<<"dijet(ratio): "<<dijet_ratio[5]->Integral()<<endl;
  cout<<"ttbar(template): "<<wjets_template[1]->Integral()<<endl;
  cout<<"ttbar(ratio): "<<wjets_ratio[1]->Integral()<<endl;
  cout<<"DY(template): "<<wjets_template[0]->Integral()<<endl;
  cout<<"DY(ratio): "<<wjets_ratio[0]->Integral()<<endl;
  cout<<"Wjets(template) after: "<<wjets_template[5]->Integral()<<endl;
  cout<<"Wjets(ratio) after: "<<wjets_ratio[5]->Integral()<<endl;

  //setTDRStyle();
  //tdrGrid(true);

  //lumiTextSize = 0.6;
  //cmsTextSize = 1.0;

  for(int i=1; i<44; i++) {
    if(wjets_template[5]->GetBinContent(i) < 0) {
      wjets_template[5]->SetBinContent(i,0.0);
      wjets_template[5]->SetBinError(i,0.0);
    }

    if(wjets_ratio[5]->GetBinContent(i) < 0) {
      wjets_ratio[5]->SetBinContent(i,0.0);
      wjets_ratio[5]->SetBinError(i,0.0);
    }

    if( wjetsSS_template[5]->GetBinContent(i) < 0 ) {
      wjetsSS_template[5]->SetBinContent(i,0.0);
      wjetsSS_template[5]->SetBinError(i,0.0);
    }
  }

  //Smooth
  wjets_template[5]->Smooth();
  wjets_ratio[5]->Smooth();
  wjetsSS_template[5]->Smooth();
  wjetsSS_ratio[5]->Smooth();

  wjets_template[5]->Write("wjets_template_smooth");
  wjetsSS_template[5]->Write("wjetsSS_template_smooth");
  wjets_ratio[5]->Write("wjets_ratio_smooth");
  wjetsSS_ratio[5]->Write("wjetsSS_ratio_smooth");

  /*TCanvas* canv = new TCanvas("canv","",1200,1200);
    canv->SetFillColor(0);
    canv->SetLeftMargin( L/W );
    canv->SetRightMargin( R/W );
    canv->SetTopMargin( T/H );
    canv->SetBottomMargin( B/H );

    wjets_template[5]->Draw("HIST");
  //CMS_lumi(canv,4,11);
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  canv->SetLogx();
  canv->Print("print/wjets_template.pdf");
  wjetsSS_template[5]->Draw("HISTPSAME");
  legg->Draw("SAME");
  canv->Print("print/wjetsBoth_template.pdf");
  canv->Clear();
  wjetsSS_template[5]->SetFillColor(9);
  wjetsSS_template[5]->Draw("HIST");
  //CMS_lumi(canv,4,11);
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  canv->SetLogx();
  canv->Print("print/wjetsSS_template.pdf");
  canv->Clear();

  wjets_ratio[5]->Draw("HIST");
  //CMS_lumi(canv,4,11);
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  canv->SetLogx();
  canv->Print("print/wjets_ratio.pdf");
  legg->Draw("SAME");
  wjetsSS_ratio[5]->Draw("HISTPSAME");
  canv->Print("print/wjetsBoth_ratio.pdf");
  canv->Clear();
  wjetsSS_ratio[5]->SetFillColor(9);
  wjetsSS_ratio[5]->Draw("HIST");
  //CMS_lumi(canv,4,11);
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  canv->SetLogx();
  canv->Print("print/wjetsSS_ratio.pdf");
  canv->Clear();*/

  double error = 0;
  wjets_template[5]->IntegralAndError(1,43,error);
  cout<<"WJets(template) OS: "<<wjets_template[5]->Integral(1,43)<<"+-"<<error<<endl;
  error = 0;
  wjetsSS_template[5]->IntegralAndError(1,43,error);
  cout<<"WJets(template) SS: "<<wjetsSS_template[5]->Integral(1,43)<<"+-"<<error<<endl;
  error = 0;
  wjets_ratio[5]->IntegralAndError(1,43,error);
  cout<<"WJets(ratio) OS: "<<wjets_ratio[5]->Integral(1,43)<<"+-"<<error<<endl;
  error = 0;
  wjetsSS_ratio[5]->IntegralAndError(1,43,error);
  cout<<"WJets(ratio) SS: "<<wjetsSS_ratio[5]->Integral(1,43)<<"+-"<<error<<endl;
  error = 0;

  TH1D* wjets = (TH1D*)wjets_template[5]->Clone();
  TH1D* wjets_control = (TH1D*)wjets_ratio[5]->Clone();
  wjets->SetName("wjets");

  TH1D* wjets_total      = new TH1D("wjets_total","",binnum,bins);
  TH1D* wjets_systematic = new TH1D("wjets_systematic","",binnum,bins);
  TH1D* wjets_stat       = new TH1D("wjets_stat","",binnum,bins);

  for(int i=1; i<binnum+1; i++) {
    double systematic = fabs(wjets->GetBinContent(i) - wjets_control->GetBinContent(i));
    double stat = wjets->GetBinError(i);
    double total = sqrt( systematic*systematic + stat*stat );
    if(wjets->GetBinContent(i)==0) {
      systematic = 0;
      stat = 0;
      total = 0;
    }

    wjets_systematic->SetBinContent(i,systematic);
    wjets_stat->SetBinContent(i,stat);
    wjets_total->SetBinContent(i,total);

    wjets->SetBinError(i,total);
  }

  //TFile* gg = new TFile("result/wjets.root","RECREATE");
  wjets->Write();
  wjets_systematic->Write();
  wjets_stat->Write();
  gg->Close();

  //wjets_systematic->Divide(wjets);
  /*wjets_systematic->Divide(wjets);
    wjets_stat->Divide(wjets);
    wjets_total->Divide(wjets);

    wjets_total->SetMarkerStyle(20);
    wjets_total->SetMarkerSize(3);
    wjets_total->SetMarkerColor(1);

    wjets_systematic->SetMarkerStyle(22);
    wjets_systematic->SetMarkerSize(3);
    wjets_systematic->SetMarkerColor(2);

    wjets_stat->SetMarkerStyle(21);
    wjets_stat->SetMarkerSize(3);
    wjets_stat->SetMarkerColor(4);

    TLegend* leg = new TLegend(.8,.15,.95,.3);
    leg->AddEntry(wjets_total,"Total","P");
    leg->AddEntry(wjets_systematic,"Sys.","P");
    leg->AddEntry(wjets_stat,"Stat.","P");

    massFrame->GetYaxis()->SetTitle("Unceratinty");
    massFrame->SetMinimum(0);
    massFrame->SetMaximum(1);
    massFrame->GetXaxis()->SetLabelSize(0.025);
    massFrame->GetXaxis()->SetTitleSize(0.05);

    massFrame->Draw();
    wjets_total->Draw("HISTPSAME");
    wjets_systematic->Draw("HISTPSAME");
    wjets_stat->Draw("HISTPSAME");
    massFrame->Draw("AXISSAME");
    leg->Draw("SAME");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    canv->Print("print/wjets_uncertainty.pdf");*/
}
