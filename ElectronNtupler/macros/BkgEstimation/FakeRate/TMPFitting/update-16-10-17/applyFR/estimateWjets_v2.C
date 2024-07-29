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

void estimateWjets_v2() {

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
  f[0]  = new TFile("histograms/nocorr/fake2.root","READ");
  f[3]  = new TFile("histograms/nocorr/fake1.root","READ");
  f[1]  = new TFile("histograms/nocorr/fake3.root","READ");
  f[6]  = new TFile("histograms/nocorr/fake6.root","READ");
  f[14] = new TFile("histograms/nocorr/fake10.root","READ");

  TFile* f1 = new TFile("ROOTFile_EffSF.root","READ");
  TFile* gg = new TFile("result/nocorr/wjets_v2.root","RECREATE");

  TH1D* wjets_template[15];
  TH1D* dijet_template[15];

  dijet_template[5] = (TH1D*)f[6]->Get("histDijet1");
  dijet_template[0] = (TH1D*)f[0]->Get("histDijet1");
  dijet_template[3] = (TH1D*)f[3]->Get("histDijet1");
  dijet_template[1] = (TH1D*)f[1]->Get("histDijet1");
  dijet_template[14] = (TH1D*)f[14]->Get("histDijet1");

  /////////
  wjets_template[5] = (TH1D*)f[6]->Get("histWJets1");
  wjets_template[0] = (TH1D*)f[0]->Get("histWJets1");
  wjets_template[3] = (TH1D*)f[3]->Get("histWJets1");
  wjets_template[1] = (TH1D*)f[1]->Get("histWJets1");
  wjets_template[14] = (TH1D*)f[14]->Get("histWJets1")->Clone();

  cout<<"data(template) OS: "<<wjets_template[5]->Integral()<<endl;
  wjets_template[5]->Write("data_templates");

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
  wjets_template[0]->Scale(norm[0]);
  wjets_template[3]->Scale(norm[3]);
  wjets_template[0]->Add(wjets_template[3]);
  wjets_template[1]->Scale(norm[1]);
  wjets_template[14]->Scale(norm[4]);

  ////
  dijet_template[0]->Scale(norm[0]);
  dijet_template[3]->Scale(norm[3]);
  dijet_template[0]->Add(dijet_template[3]);
  dijet_template[1]->Scale(norm[1]);
  dijet_template[14]->Scale(norm[14]);

  TH1D *h_num = (TH1D*)f1->Get("hist_mass_num");
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

  double dy2 = wjets_template[0]->GetBinContent(i);
  double dy2err = wjets_template[0]->GetBinError(i);
  double ttbar2 = wjets_template[1]->GetBinContent(i);
  double tt2err = wjets_template[1]->GetBinError(i);
  double gjets2 = wjets_template[14]->GetBinContent(i);
  double gj2err = wjets_template[14]->GetBinError(i);

  double dy1Err = dy1*scalefactor * sqrt( (dy1err*dy1err)/(dy1*dy1) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double tt1Err = ttbar1*scalefactor * sqrt( (tt1err*tt1err)/(ttbar1*ttbar1) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double gj1Err = gjets1*scalefactor * sqrt( (gj1err*gj1err)/(gjets1*gjets1) + (SFErr*SFErr)/(scalefactor*scalefactor) );

  double dy2Err = dy2*scalefactor * sqrt( (dy2err*dy2err)/(dy2*dy2) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double tt2Err = ttbar2*scalefactor * sqrt( (tt2err*tt2err)/(ttbar2*ttbar2) + (SFErr*SFErr)/(scalefactor*scalefactor) );
  double gj2Err = gjets2*scalefactor * sqrt( (gj2err*gj2err)/(gjets2*gjets2) + (SFErr*SFErr)/(scalefactor*scalefactor) );

  //cout<<"Before Content = "<<dy<<"   "<<ttbar<<endl;
  //cout<<"Before Error   = "<<dyerr<<"   "<<tterr<<endl;

  // Dijet
  dijet_template[0]->SetBinContent(i, scalefactor*dy1);
  if(fabs(dy1) > 0) dijet_template[0]->SetBinError(i, dy1Err);
  else dijet_template[0]->SetBinError(i, 0);

  dijet_template[1]->SetBinContent(i, scalefactor*ttbar1);
  if(fabs(ttbar1) > 0) dijet_template[1]->SetBinError(i, tt1Err);
  else dijet_template[1]->SetBinError(i, 0);

  dijet_template[14]->SetBinContent(i, scalefactor*gjets1);
  if(fabs(gjets1) > 0) dijet_template[14]->SetBinError(i, gj1Err);
  else dijet_template[14]->SetBinError(i, 0);

  // W+jets
  wjets_template[0]->SetBinContent(i, scalefactor*dy2);
  if(fabs(dy2) > 0) wjets_template[0]->SetBinError(i, dy2Err);
  else wjets_template[0]->SetBinError(i, 0);

  wjets_template[1]->SetBinContent(i, scalefactor*ttbar2);
  if(fabs(ttbar2) > 0) wjets_template[1]->SetBinError(i, tt2Err);
  else wjets_template[1]->SetBinError(i, 0);

  wjets_template[14]->SetBinContent(i, scalefactor*gjets2);
  if(fabs(gjets2) > 0) wjets_template[14]->SetBinError(i, gj2Err);
  else wjets_template[14]->SetBinError(i, 0);

  //cout<<"After Content  = "<<dijet_template[0]->GetBinContent(i)<<"   "<<dijet_template[1]->GetBinContent(i)<<endl;
  //cout<<"After Error    = "<<dijet_template[0]->GetBinError(i)<<"   "<<dijet_template[1]->GetBinError(i)<<endl;
  //cout<<""<<endl;

  }

  cout<<"DY(template) OS: "<<wjets_template[0]->Integral()<<endl;
  cout<<"ttbar(template) OS: "<<wjets_template[1]->Integral()<<endl;
  cout<<"GJets(template) OS: "<<wjets_template[14]->Integral()<<endl;

  wjets_template[0]->Write("dy_template");
  wjets_template[1]->Write("ttbar_template");
  wjets_template[14]->Write("gjets_template");

  // Dijet
  /*for(int i=1; i<44; i++) {
    if(dijet_template[0]->GetBinContent(i) < 0) {
      dijet_template[0]->SetBinContent(i,0.0);
      dijet_template[0]->SetBinError(i,0.0);
    }
  }*/

  //dijet_template[5]->Add(dijet_template[0],-1.0);
  //dijet_template[5]->Add(dijet_template[1],-1.0);
  //dijet_template[5]->Add(dijet_template[14],-1.0);

  dijet_template[5]->SetBinContent(1,4.971472);
  dijet_template[5]->SetBinContent(2,4.519224);
  dijet_template[5]->SetBinContent(3,3.971289);
  dijet_template[5]->SetBinContent(4,4.993672);
  dijet_template[5]->SetBinContent(5,8.349676);
  dijet_template[5]->SetBinContent(6,13.45228);
  dijet_template[5]->SetBinContent(7,15.51968);
  dijet_template[5]->SetBinContent(8,14.46188);
  dijet_template[5]->SetBinContent(9,14.18202);
  dijet_template[5]->SetBinContent(10,12.6324);
  dijet_template[5]->SetBinContent(11,11.5465);
  dijet_template[5]->SetBinContent(12,9.929977);
  dijet_template[5]->SetBinContent(13,5.892304);
  dijet_template[5]->SetBinContent(14,7.517796);
  dijet_template[5]->SetBinContent(15,8.987323);
  dijet_template[5]->SetBinContent(16,8.987323);
  dijet_template[5]->SetBinContent(17,8.987323);
  dijet_template[5]->SetBinContent(18,8.987323);
  dijet_template[5]->SetBinContent(19,7.761368);
  dijet_template[5]->SetBinContent(20,4.971742);
  dijet_template[5]->SetBinContent(21,6.162208);
  dijet_template[5]->SetBinContent(22,5.395173);
  dijet_template[5]->SetBinContent(23,6.818448);
  dijet_template[5]->SetBinContent(24,7.419263);
  dijet_template[5]->SetBinContent(25,7.32303);
  dijet_template[5]->SetBinContent(26,7.489979);
  dijet_template[5]->SetBinContent(27,6.341876);
  dijet_template[5]->SetBinContent(28,6.527589);
  dijet_template[5]->SetBinContent(29,6.238975);
  dijet_template[5]->SetBinContent(30,6.564028);
  dijet_template[5]->SetBinContent(31,6.506175);
  dijet_template[5]->SetBinContent(32,6.228585);
  dijet_template[5]->SetBinContent(33,6.724893);
  dijet_template[5]->SetBinContent(34,7.021114);
  dijet_template[5]->SetBinContent(35,5.805506);
  dijet_template[5]->SetBinContent(36,3.256544);
  dijet_template[5]->SetBinContent(37,1.552025);
  dijet_template[5]->SetBinContent(38,1.7438);
  dijet_template[5]->SetBinContent(39,1.287416);
  dijet_template[5]->SetBinContent(40,0.8415104);
  dijet_template[5]->SetBinContent(41,0.6650521);
  dijet_template[5]->SetBinContent(42,0.45739);
  dijet_template[5]->SetBinContent(43,0.330468);
  dijet_template[5]->SetBinError(1,0.7263963);
  dijet_template[5]->SetBinError(2,0.6007138);
  dijet_template[5]->SetBinError(3,0.5490095);
  dijet_template[5]->SetBinError(4,0.651571);
  dijet_template[5]->SetBinError(5,0.6764778);
  dijet_template[5]->SetBinError(6,0.7757187);
  dijet_template[5]->SetBinError(7,1.031641);
  dijet_template[5]->SetBinError(8,0.8987661);
  dijet_template[5]->SetBinError(9,1.007484);
  dijet_template[5]->SetBinError(10,0.9025045);
  dijet_template[5]->SetBinError(11,0.9308435);
  dijet_template[5]->SetBinError(12,0.9497515);
  dijet_template[5]->SetBinError(13,1.287991);
  dijet_template[5]->SetBinError(14,1.34219);
  dijet_template[5]->SetBinError(15,4.28108);
  dijet_template[5]->SetBinError(16,4.28108);
  dijet_template[5]->SetBinError(17,4.28108);
  dijet_template[5]->SetBinError(18,4.28108);
  dijet_template[5]->SetBinError(19,0.8039756);
  dijet_template[5]->SetBinError(20,0.5832915);
  dijet_template[5]->SetBinError(21,0.6232218);
  dijet_template[5]->SetBinError(22,1.603471);
  dijet_template[5]->SetBinError(23,0.5873509);
  dijet_template[5]->SetBinError(24,0.6100289);
  dijet_template[5]->SetBinError(25,0.6069087);
  dijet_template[5]->SetBinError(26,0.600678);
  dijet_template[5]->SetBinError(27,0.5358132);
  dijet_template[5]->SetBinError(28,0.5423073);
  dijet_template[5]->SetBinError(29,0.5580924);
  dijet_template[5]->SetBinError(30,0.5754691);
  dijet_template[5]->SetBinError(31,0.6225296);
  dijet_template[5]->SetBinError(32,0.5758989);
  dijet_template[5]->SetBinError(33,0.6068115);
  dijet_template[5]->SetBinError(34,0.6575865);
  dijet_template[5]->SetBinError(35,0.6496811);
  dijet_template[5]->SetBinError(36,0.4756398);
  dijet_template[5]->SetBinError(37,0.3328225);
  dijet_template[5]->SetBinError(38,0.3713172);
  dijet_template[5]->SetBinError(39,0.4008877);
  dijet_template[5]->SetBinError(40,0.303381);
  dijet_template[5]->SetBinError(41,0.2523161);
  dijet_template[5]->SetBinError(42,0.2338805);
  dijet_template[5]->SetBinError(43,0.3330084);

  /*for(int i=1; i<44; i++) {
    if(dijet_template[5]->GetBinContent(i) < 0) {
      dijet_template[5]->SetBinContent(i,0.0);
      dijet_template[5]->SetBinError(i,0.0);
    }
  }*/

  // W+jets
  for(int i=1; i<44; i++) {
    if(wjets_template[0]->GetBinContent(i) < 0) {
      wjets_template[0]->SetBinContent(i,0.0);
      wjets_template[0]->SetBinError(i,0.0);
    }
  }

  wjets_template[5]->Add(dijet_template[5],-2.0);
  wjets_template[5]->Add(wjets_template[0],-1.0);
  wjets_template[5]->Add(wjets_template[1],-1.0);
  wjets_template[5]->Add(wjets_template[14],-1.0);

  /*for(int i=1; i<44; i++) {
    if(wjets_template[5]->GetBinContent(i) < 0) {
      wjets_template[5]->SetBinContent(i,0.0);
      wjets_template[5]->SetBinError(i,0.0);
    } 
  }*/

  double error = 0;
  dijet_template[5]->IntegralAndError(1,43,error);
  cout<<"QCD(template) OS = "<<dijet_template[5]->Integral(1,43)<<"+-"<<error<<endl;
  error = 0;
  wjets_template[5]->IntegralAndError(1,43,error);
  cout<<"W+Jets(template) OS = "<<wjets_template[5]->Integral(1,43)<<"+-"<<error<<endl;

  dijet_template[5]->Write("dijet_template");
  wjets_template[5]->Write("wjets_template");

  dijet_template[5]->Smooth();
  dijet_template[5]->Write("dijet_template_smooth");
  wjets_template[5]->Smooth();
  wjets_template[5]->Write("wjets_template_smooth");

  gg->Close();

}
