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

void estimateDijet() {

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
  //lumi_13TeV = "2833 pb^{-1}";
  //lumiTextSize = 0.5;
  //writeExtraText = true;
  //extraText = "Preliminary";
  //drawLogo = false;

  int binnum = 43;
  double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

  int binnum1 = 42;
  double bins1[43] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

  TH1D* massFrame = new TH1D("massFrame","",38,15,3000);
  massFrame->SetMinimum(0.001);
  massFrame->SetMaximum(1000000);
  massFrame->SetStats(kFALSE);
  massFrame->GetXaxis()->SetTitle("Mass[GeV]");
  massFrame->GetYaxis()->SetTitleOffset(1);
  massFrame->GetYaxis()->SetTitle("Number of events");
  massFrame->GetXaxis()->SetTitleSize(0);
  massFrame->GetYaxis()->SetTitleSize(0.05);
  massFrame->GetXaxis()->SetLabelSize(0);
  massFrame->GetXaxis()->SetMoreLogLabels();

  TFile* f[15];
  f[0]  = new TFile("histograms/singlebin/fake2.root","READ");
  f[3]  = new TFile("histograms/singlebin/fake1.root","READ");
  f[1]  = new TFile("histograms/singlebin/fake3.root","READ");
  f[6]  = new TFile("histograms/singlebin/fake6.root","READ");
  f[14] = new TFile("histograms/singlebin/fake10.root","READ");

  TFile* f1 = new TFile("ROOTFile_EffSF_NEW.root","READ");
  TFile* gg = new TFile("result/singlebin/dijet.root","RECREATE");

  TH1D* wjets_template[15]; // just for draw MC histograms
  TH1D* dijet_template[15];
  TH1D* dijetSS_template[15];
  TH1D* dijet_ratio[15];
  TH1D* dijetSS_ratio[15];

  dijet_template[5] = (TH1D*)f[6]->Get("histDijet1_81_101_1");
  //dijet_template[5]->SetFillColor(7);
  dijet_template[5]->SetStats(kFALSE);

  dijetSS_template[5] = (TH1D*)f[6]->Get("histSameDijet1");
  dijetSS_template[0] = (TH1D*)f[0]->Get("histSameDijet1");
  dijetSS_template[3] = (TH1D*)f[3]->Get("histSameDijet1");
  dijetSS_template[1] = (TH1D*)f[1]->Get("histSameDijet1");
  dijetSS_template[14] = (TH1D*)f[14]->Get("histSameDijet1");
  dijetSS_template[5]->SetLineColor(1);
  dijetSS_template[5]->SetMarkerColor(1);
  dijetSS_template[5]->SetMarkerStyle(22);
  dijetSS_template[5]->SetMarkerSize(1.5);
  dijetSS_template[5]->SetStats(kFALSE);

  cout<<"data(template) OS: "<<dijet_template[5]->Integral()<<endl;
  cout<<"data(template) SS: "<<dijetSS_template[5]->Integral()<<endl;
  dijet_template[5]->Write("data_templates");
  dijetSS_template[5]->Write("dataSS_templates");

  /////////
  //wjets_template[0] = (TH1D*)f[0]->Get("histWJets1");
  //wjets_template[3] = (TH1D*)f[3]->Get("histWJets1");
  //wjets_template[0]->SetFillColor(2);
  //wjets_template[0]->SetStats(kFALSE);

  //wjets_template[1] = (TH1D*)f[1]->Get("histWJets1");
  //wjets_template[1]->SetFillColor(3);
  //wjets_template[1]->SetStats(kFALSE);

  /////////
  dijet_template[0] = (TH1D*)f[0]->Get("histDijet1_81_101_1");
  dijet_template[3] = (TH1D*)f[3]->Get("histDijet1_81_101_1");
  //dijet_template[0]->SetFillColor(2);
  dijet_template[0]->SetStats(kFALSE);

  dijet_template[1] = (TH1D*)f[1]->Get("histDijet1_81_101_1");
  //dijet_template[1]->SetFillColor(3);
  dijet_template[1]->SetStats(kFALSE);

  dijet_template[14] = (TH1D*)f[14]->Get("histDijet1_81_101_1");

  dijet_ratio[5] = (TH1D*)f[6]->Get("histDijet2_81_101_1");
  //dijet_ratio[5]->SetFillColor(7);
  dijet_ratio[5]->SetStats(kFALSE);

  dijetSS_ratio[5] = (TH1D*)f[6]->Get("histSameDijet2");
  dijetSS_ratio[0] = (TH1D*)f[0]->Get("histSameDijet2");
  dijetSS_ratio[3] = (TH1D*)f[3]->Get("histSameDijet2");
  dijetSS_ratio[1] = (TH1D*)f[1]->Get("histSameDijet2");
  dijetSS_ratio[14] = (TH1D*)f[14]->Get("histSameDijet2");
  dijetSS_ratio[5]->SetLineColor(1);
  dijetSS_ratio[5]->SetMarkerColor(1);
  dijetSS_ratio[5]->SetMarkerStyle(22);
  dijetSS_ratio[5]->SetMarkerSize(1.5);
  dijetSS_ratio[5]->SetStats(kFALSE);

  cout<<"data(ratio) OS: "<<dijet_ratio[5]->Integral()<<endl;
  cout<<"data(ratio) SS: "<<dijetSS_ratio[5]->Integral()<<endl;
  dijet_ratio[5]->Write("data_ratio");
  dijetSS_ratio[5]->Write("dataSS_ratio");

  dijet_ratio[0] = (TH1D*)f[0]->Get("histDijet2_81_101_1");
  dijet_ratio[3] = (TH1D*)f[3]->Get("histDijet2_81_101_1");
  //dijet_ratio[0]->SetFillColor(2);
  dijet_ratio[0]->SetStats(kFALSE);

  dijet_ratio[1] = (TH1D*)f[1]->Get("histDijet2_81_101_1");
  //dijet_ratio[1]->SetFillColor(3);
  dijet_ratio[1]->SetStats(kFALSE);

  dijet_ratio[14] = (TH1D*)f[14]->Get("histDijet2_81_101_1");

  //DYMuMu, ttbar, WJets, WW, tautau, QCD
  double nEvts[15] = {4.5151e+11, 85849572, 3.73193e+12, 2.30899e+12, 3268361, 2.23076e+06, 0, 0, 0, 0, 0, 988416, 999996, 985598, 54136};
  double xsec[15] = {2008.4*3, 831.76, 61526.7, 18610, 1915, 2.23076e+06, 0, 0, 0, 0, 0, 118.7, 47.13, 16.523, 365896};
  double norm[15];
  // UPDATED IN 2017
  double lumi = 2258.066;
  // double lumi = 2832.673;

  for(int i=0;i<15;i++) {
    if(i>5&&i<11) continue;
    norm[i] = (xsec[i]*lumi)/nEvts[i];
    //cout<<norm[i]<<endl;
  }

  ////
  //wjets_template[0]->Scale(norm[0]);
  //wjets_template[3]->Scale(norm[3]);
  //wjets_template[0]->Add(wjets_template[3]);
  //wjets_template[1]->Scale(norm[1]);

  ////
  dijet_template[0]->Scale(norm[0]);
  dijet_template[3]->Scale(norm[3]);
  dijet_template[0]->Add(dijet_template[3]);
  dijet_template[1]->Scale(norm[1]);
  dijet_template[14]->Scale(norm[14]);

  dijetSS_template[0]->Scale(norm[0]);
  dijetSS_template[3]->Scale(norm[3]);
  dijetSS_template[0]->Add(dijetSS_template[3]);
  dijetSS_template[1]->Scale(norm[1]);
  dijetSS_template[14]->Scale(norm[14]);

  dijet_ratio[0]->Scale(norm[0]);
  dijet_ratio[3]->Scale(norm[3]);
  dijet_ratio[0]->Add(dijet_ratio[3]);
  dijet_ratio[1]->Scale(norm[1]);
  dijet_ratio[14]->Scale(norm[14]);

  dijetSS_ratio[0]->Scale(norm[0]);
  dijetSS_ratio[3]->Scale(norm[3]);
  dijetSS_ratio[0]->Add(dijetSS_ratio[3]);
  dijetSS_ratio[1]->Scale(norm[1]);
  dijetSS_ratio[14]->Scale(norm[14]);

  TH1D *h_sf = (TH1D*)f1->Get("h_EffSF2");

  for(int i=1; i<43; i++){

    //cout<<dijet_template[5]->GetBinLowEdge(i)<<"   "<<dijet_template[5]->GetBinLowEdge(i+1)<<endl;

    double scalefactor = h_sf->GetBinContent(i);
    double SFErr = h_sf->GetBinError(i);

    double dy1 = dijet_template[0]->GetBinContent(i);
    double dy1err = dijet_template[0]->GetBinError(i);
    double ttbar1 = dijet_template[1]->GetBinContent(i);
    double tt1err = dijet_template[1]->GetBinError(i);
    double gjets1 = dijet_template[14]->GetBinContent(i);
    double gj1err = dijet_template[14]->GetBinError(i);

    double dy2 = dijet_ratio[0]->GetBinContent(i);
    double dy2err = dijet_ratio[0]->GetBinError(i);
    double ttbar2 = dijet_ratio[1]->GetBinContent(i);
    double tt2err = dijet_ratio[1]->GetBinError(i);
    double gjets2 = dijet_ratio[14]->GetBinContent(i);
    double gj2err = dijet_ratio[14]->GetBinError(i);

    double dy1Err = dy1*scalefactor * sqrt( (dy1err*dy1err)/(dy1*dy1) + (SFErr*SFErr)/(scalefactor*scalefactor) );
    double tt1Err = ttbar1*scalefactor * sqrt( (tt1err*tt1err)/(ttbar1*ttbar1) + (SFErr*SFErr)/(scalefactor*scalefactor) );
    double gj1Err = gjets1*scalefactor * sqrt( (gj1err*gj1err)/(gjets1*gjets1) + (SFErr*SFErr)/(scalefactor*scalefactor) );

    double dy2Err = dy2*scalefactor * sqrt( (dy2err*dy2err)/(dy2*dy2) + (SFErr*SFErr)/(scalefactor*scalefactor) );
    double tt2Err = ttbar2*scalefactor * sqrt( (tt2err*tt2err)/(ttbar2*ttbar2) + (SFErr*SFErr)/(scalefactor*scalefactor) );
    double gj2Err = gjets2*scalefactor * sqrt( (gj2err*gj2err)/(gjets2*gjets2) + (SFErr*SFErr)/(scalefactor*scalefactor) );

    //cout<<"Before Content = "<<dy1<<"   "<<ttbar1<<"   "<<gjets1<<endl;
    //cout<<"Before Error   = "<<dy1err<<"   "<<tt1err<<"   "<<gj1err<<endl;

    dijet_template[0]->SetBinContent(i, scalefactor*dy1);
    if(fabs(dy1) > 0) dijet_template[0]->SetBinError(i, dy1Err);
    else dijet_template[0]->SetBinError(i, 0);

    dijet_template[1]->SetBinContent(i, scalefactor*ttbar1);
    if(fabs(ttbar1) > 0) dijet_template[1]->SetBinError(i, tt1Err);
    else dijet_template[1]->SetBinError(i, 0);

    dijet_template[14]->SetBinContent(i, scalefactor*gjets1);
    if(fabs(gjets1) > 0) dijet_template[14]->SetBinError(i, gj1Err);
    else dijet_template[14]->SetBinError(i, 0);

    dijet_ratio[0]->SetBinContent(i, scalefactor*dy2);
    if(fabs(dy2) > 0) dijet_ratio[0]->SetBinError(i, dy2Err);
    else dijet_ratio[0]->SetBinError(i, 0);

    dijet_ratio[1]->SetBinContent(i, scalefactor*ttbar2);
    if(fabs(ttbar2) > 0) dijet_ratio[1]->SetBinError(i, tt2Err);
    else dijet_ratio[1]->SetBinError(i, 0);

    dijet_ratio[14]->SetBinContent(i, scalefactor*gjets2);
    if(fabs(gjets2) > 0) dijet_ratio[14]->SetBinError(i, gj2Err);
    else dijet_ratio[14]->SetBinError(i, 0);

    //cout<<"After Content  = "<<dijet_template[0]->GetBinContent(i)<<"   "<<dijet_template[1]->GetBinContent(i)<<"   "<<dijet_template[14]->GetBinContent(i)<<endl;
    //cout<<"After Error    = "<<dijet_template[0]->GetBinError(i)<<"   "<<dijet_template[1]->GetBinError(i)<<"   "<<dijet_template[14]->GetBinError(i)<<endl;
    //cout<<""<<endl;
  }

  /*TH1D *h_num = (TH1D*)f1->Get("hist_mass_num");
    TH1D *h_den = (TH1D*)f1->Get("hist_mass_den");

    for(int i=1; i<44; i++){

  //cout<<dijet_template[5]->GetBinLowEdge(i)<<"   "<<dijet_template[5]->GetBinLowEdge(i+1)<<endl;

  double numEff = h_num->GetBinContent(i);
  double numErr = h_num->GetBinError(i);
  double denEff = h_den->GetBinContent(i);
  double denErr = h_den->GetBinError(i);

  double scalefactor = numEff/denEff;
  double denerr1 = (h_den->GetBinContent(i)*h_den->GetBinContent(i))/(h_den->GetBinError(i)*h_den->GetBinError(i));
  double SFErr = sqrt(scalefactor*(1-scalefactor)/denerr1);

  //cout<<"SF = "<<scalefactor<<endl;

  double dy1 = dijet_template[0]->GetBinContent(i);
  double dy1err = dijet_template[0]->GetBinError(i);
  double ttbar1 = dijet_template[1]->GetBinContent(i);
  double tt1err = dijet_template[1]->GetBinError(i);
  double gjets1 = dijet_template[14]->GetBinContent(i);
  double gj1err = dijet_template[14]->GetBinError(i);

  double dy2 = dijet_ratio[0]->GetBinContent(i);
  double dy2err = dijet_ratio[0]->GetBinError(i);
  double ttbar2 = dijet_ratio[1]->GetBinContent(i);
  double tt2err = dijet_ratio[1]->GetBinError(i);
  double gjets2 = dijet_ratio[14]->GetBinContent(i);
  double gj2err = dijet_ratio[14]->GetBinError(i);

  double dy1Err = dy1*scalefactor * sqrt( (dy1err*dy1err)/(dy1*dy1) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double tt1Err = ttbar1*scalefactor * sqrt( (tt1err*tt1err)/(ttbar1*ttbar1) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double gj1Err = gjets1*scalefactor * sqrt( (gj1err*gj1err)/(gjets1*gjets1) + (SFErr*SFErr)/(scalefactor*scalefactor) );

  double dy2Err = dy2*scalefactor * sqrt( (dy2err*dy2err)/(dy2*dy2) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double tt2Err = ttbar2*scalefactor * sqrt( (tt2err*tt2err)/(ttbar2*ttbar2) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double gj2Err = gjets2*scalefactor * sqrt( (gj2err*gj2err)/(gjets2*gjets2) + (SFErr*SFErr)/(scalefactor*scalefactor) );

  //cout<<"Before Content = "<<dy<<"   "<<ttbar<<endl;
  //cout<<"Before Error   = "<<dyerr<<"   "<<tterr<<endl;

  dijet_template[0]->SetBinContent(i, scalefactor*dy1);
  if(fabs(dy1) > 0) dijet_template[0]->SetBinError(i, dy1Err);
  else dijet_template[0]->SetBinError(i, 0);

  dijet_template[1]->SetBinContent(i, scalefactor*ttbar1);
  if(fabs(ttbar1) > 0) dijet_template[1]->SetBinError(i, tt1Err);
  else dijet_template[1]->SetBinError(i, 0);

  dijet_template[14]->SetBinContent(i, scalefactor*gjets1);
  if(fabs(gjets1) > 0) dijet_template[14]->SetBinError(i, gj1Err);
  else dijet_template[14]->SetBinError(i, 0);

  dijet_ratio[0]->SetBinContent(i, scalefactor*dy2);
  if(fabs(dy2) > 0) dijet_ratio[0]->SetBinError(i, dy2Err);
  else dijet_ratio[0]->SetBinError(i, 0);

  dijet_ratio[1]->SetBinContent(i, scalefactor*ttbar2);
  if(fabs(ttbar2) > 0) dijet_ratio[1]->SetBinError(i, tt2Err);
  else dijet_ratio[1]->SetBinError(i, 0);

  dijet_ratio[14]->SetBinContent(i, scalefactor*gjets2);
  if(fabs(gjets2) > 0) dijet_ratio[14]->SetBinError(i, gj2Err);
  else dijet_ratio[14]->SetBinError(i, 0);

  //cout<<"After Content  = "<<dijet_template[0]->GetBinContent(i)<<"   "<<dijet_template[1]->GetBinContent(i)<<endl;
  //cout<<"After Error    = "<<dijet_template[0]->GetBinError(i)<<"   "<<dijet_template[1]->GetBinError(i)<<endl;
  //cout<<""<<endl;

}*/

cout<<"DY(template) OS: "<<dijet_template[0]->Integral()<<endl;
cout<<"DY(template) SS: "<<dijetSS_template[0]->Integral()<<endl;
cout<<"DY(ratio) OS: "<<dijet_ratio[0]->Integral()<<endl;
cout<<"DY(ratio) SS: "<<dijetSS_ratio[0]->Integral()<<endl;

cout<<"ttbar(template) OS: "<<dijet_template[1]->Integral()<<endl;
cout<<"ttbar(template) SS: "<<dijetSS_template[1]->Integral()<<endl;
cout<<"ttbar(ratio) OS: "<<dijet_ratio[1]->Integral()<<endl;
cout<<"ttbar(ratio) SS: "<<dijetSS_ratio[1]->Integral()<<endl;

cout<<"GJets(template) OS: "<<dijet_template[14]->Integral()<<endl;
cout<<"GJets(template) SS: "<<dijetSS_template[14]->Integral()<<endl;
cout<<"GJets(ratio) OS: "<<dijet_ratio[14]->Integral()<<endl;
cout<<"GJets(ratio) SS: "<<dijetSS_ratio[14]->Integral()<<endl;

dijet_template[0]->Write("dy_template");
dijet_template[1]->Write("ttbar_template");
dijet_template[14]->Write("gjets_template");
dijetSS_template[0]->Write("dySS_template");
dijetSS_template[1]->Write("ttbarSS_template");
dijetSS_template[14]->Write("gjetsSS_template");

dijet_ratio[0]->Write("dy_ratio");
dijet_ratio[1]->Write("ttbar_ratio");
dijet_ratio[14]->Write("gjets_ratio");
dijetSS_ratio[0]->Write("dySS_ratio");
dijetSS_ratio[1]->Write("ttbarSS_ratio");
dijetSS_ratio[14]->Write("gjetsSS_ratio");

for(int i=1; i<44; i++) {
  if(dijet_template[0]->GetBinContent(i) < 0) {
    dijet_template[0]->SetBinContent(i,0.0);
    dijet_template[0]->SetBinError(i,0.0);
  }
}
dijet_template[5]->Add(dijet_template[0],-1.0);
dijet_template[5]->Add(dijet_template[1],-1.0);
dijet_template[5]->Add(dijet_template[14],-1.0);

dijet_ratio[5]->Add(dijet_ratio[0],-1.0);
dijet_ratio[5]->Add(dijet_ratio[1],-1.0);
dijet_ratio[5]->Add(dijet_ratio[14],-1.0);

for(int i=1; i<44; i++) {
  if(dijet_template[5]->GetBinContent(i) < 0) {
    dijet_template[5]->SetBinContent(i,0.0);
    dijet_template[5]->SetBinError(i,0.0);
  }
  if(dijet_ratio[5]->GetBinContent(i) < 0) {
    dijet_ratio[5]->SetBinContent(i,0.0);
    dijet_ratio[5]->SetBinError(i,0.0);
  }
}

double error = 0;
dijet_template[5]->IntegralAndError(1,43,error);
cout<<"QCD(template) OS = "<<dijet_template[5]->Integral(1,43)<<"+-"<<error<<endl;
error = 0;
dijetSS_template[5]->IntegralAndError(1,43,error);
cout<<"QCD(template) SS = "<<dijetSS_template[5]->Integral(1,43)<<"+-"<<error<<endl;
error = 0;
dijet_ratio[5]->IntegralAndError(1,43,error);
cout<<"QCD(ratio) OS = "<<dijet_ratio[5]->Integral(1,43)<<"+-"<<error<<endl;
error = 0;
dijetSS_ratio[5]->IntegralAndError(1,43,error);
cout<<"QCD(ratio) SS = "<<dijetSS_ratio[5]->Integral(1,43)<<"+-"<<error<<endl;

dijet_template[5]->Write("dijet_template");
dijetSS_template[5]->Write("dijetSS_template");
dijet_ratio[5]->Write("dijet_ratio");
dijetSS_ratio[5]->Write("dijetSS_ratio");
//dijet->Write();
//dijet_systematic->Write();
//dijet_stat->Write();

dijet_template[5]->Smooth();
dijetSS_template[5]->Smooth();
dijet_ratio[5]->Smooth();
dijetSS_ratio[5]->Smooth();

dijet_template[5]->Write("dijet_template_smooth");
dijetSS_template[5]->Write("dijetSS_template_smooth");
dijet_ratio[5]->Write("dijet_ratio_smooth");
dijetSS_ratio[5]->Write("dijetSS_ratio_smooth");

gg->Close();

}
