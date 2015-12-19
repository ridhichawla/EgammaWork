{
  TFile *f0 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/muonEG_Run2015D.root");
  TFile *f1 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/dyll_M10-50.root");  
  TFile *f2 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/dyll_M-50.root");
  TFile *f3 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/tt2L2Nu.root");
  TFile *f4 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/ww2L2Nu.root");
  TFile *f5 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/wz3LNu.root");
  TFile *f6 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/zz4L.root");
  TFile *f7 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/st_sch_4f.root");
  TFile *f8 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/st_tch_4f.root");
  TFile *f9 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/st_tch_atop_4f.root");
  TFile *f10 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/st_tch_top_4f.root");
  TFile *f11 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/st_tW_atop_5f.root");
  TFile *f12 = TFile::Open("/home/ridhi/Work/Analysis/Data-MC/Tuples/13TeV/BkgEstimate/st_tW_top_5f.root");

  TFile *file = TFile::Open("eMumass.root","recreate");

  //double w1 = (2076.132*18610)/905268829598;
  //double w2 = (2076.132*6025.2)/451498239235;
  double w1 = (2076.132*6203.3)/905268829598;
  double w2 = (2076.132*2008)/451498239235;
  double w3 = (2076.132*87.31)/4997000;
  double w4 = (2076.132*10.481)/1965200;
  double w5 = (2076.132*4.42965)/1980800;
  double w6 = (2076.132*1.256)/6665004;
  double w7 = (2076.132*10.3)/3318690.136;
  double w8 = (2076.132*70.3)/1340168885;
  double w9 = (2076.132*26.2)/1680200;
  double w10 = (2076.132*44.1)/3299800;
  double w11 = (2076.132*35.6)/988500;
  double w12 = (2076.132*35.6)/995600;

//***************************************e mu Mass****************************************

  TCanvas *c1 =  new TCanvas("c1","",66,164,377,413);
  c1->Draw();
  c1->cd();
  
  c1_2 = new TPad("c1_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c1_2->SetBottomMargin(0.1);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetLogx();
  c1_2->SetLogy();
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);
  c1_2->SetGridx();
  c1_2->SetGridy();

  TH1F *h_data_EMu_Mass       = (TH1F*)f0->Get("EMu_Mass");
  TH1F *h_dyll10to50_EMu_Mass = (TH1F*)f1->Get("EMu_Mass");
  TH1F *h_dyll50_EMu_Mass     = (TH1F*)f2->Get("EMu_Mass");
  TH1F *h_tt_EMu_Mass         = (TH1F*)f3->Get("EMu_Mass");
  TH1F *h_ww_EMu_Mass         = (TH1F*)f4->Get("EMu_Mass");
  TH1F *h_wz_EMu_Mass         = (TH1F*)f5->Get("EMu_Mass");
  TH1F *h_zz_EMu_Mass         = (TH1F*)f6->Get("EMu_Mass");
  TH1F *h_s_4f_EMu_Mass       = (TH1F*)f7->Get("EMu_Mass");
  //TH1F *h_t_4f_EMu_Mass       = (TH1F*)f8->Get("EMu_Mass");
  TH1F *h_st_at_4f_EMu_Mass   = (TH1F*)f9->Get("EMu_Mass");
  TH1F *h_st_t_4f_EMu_Mass    = (TH1F*)f10->Get("EMu_Mass");
  TH1F *h_tW_at_5f_EMu_Mass   = (TH1F*)f11->Get("EMu_Mass");
  TH1F *h_tW_t_5f_EMu_Mass    = (TH1F*)f12->Get("EMu_Mass");
 
  cout<<"DY to tau tau: "<<h_dyll50_EMu_Mass->GetEntries()<<endl;
  cout<<"ttbar: "<<h_tt_EMu_Mass->GetEntries()<<endl;

  h_dyll10to50_EMu_Mass->Scale(w1);
  h_dyll50_EMu_Mass->Scale(w2);
  h_tt_EMu_Mass->Scale(w3);
  h_ww_EMu_Mass->Scale(w4);
  h_wz_EMu_Mass->Scale(w5);
  h_zz_EMu_Mass->Scale(w6);
  h_s_4f_EMu_Mass->Scale(w7);
  //h_t_4f_EMu_Mass->Scale(w8);
  h_st_at_4f_EMu_Mass->Scale(w9);
  h_st_t_4f_EMu_Mass->Scale(w10);
  h_tW_at_5f_EMu_Mass->Scale(w11);
  h_tW_t_5f_EMu_Mass->Scale(w12);

  TH1F *h_dyll_EMu_Mass = (TH1F*)h_dyll10to50_EMu_Mass->Clone("h_dyll_EMu_Mass");
  h_dyll_EMu_Mass->Add(h_dyll50_EMu_Mass);

  TH1F *h_st_at_EMu_Mass = (TH1F*)h_st_at_4f_EMu_Mass->Clone("h_st_at_EMu_Mass");
  h_st_at_EMu_Mass->Add(h_tW_at_5f_EMu_Mass);

  TH1F *h_st_t_EMu_Mass = (TH1F*)h_st_t_4f_EMu_Mass->Clone("h_st_t_EMu_Mass");
  h_st_t_EMu_Mass->Add(h_tW_t_5f_EMu_Mass);
  h_st_t_EMu_Mass->Add(h_s_4f_EMu_Mass);
  //h_st_t_EMu_Mass->Add(h_t_4f_EMu_Mass);
    
  h_data_EMu_Mass->SetMarkerColor(1);
  h_data_EMu_Mass->SetMarkerStyle(20);
  h_data_EMu_Mass->SetMarkerSize(0.5); 
  h_data_EMu_Mass->SetLineColor(1);

  h_dyll_EMu_Mass->SetLineColor(kBlue+1);
  h_dyll_EMu_Mass->SetFillColor(kBlue+1);
  h_tt_EMu_Mass->SetLineColor(kGreen+2);
  h_tt_EMu_Mass->SetFillColor(kGreen+2);
  h_ww_EMu_Mass->SetLineColor(kCyan-7);  //Purple
  h_ww_EMu_Mass->SetFillColor(kCyan-7);
  h_wz_EMu_Mass->SetLineColor(kPink-9);
  h_wz_EMu_Mass->SetFillColor(kPink-9);
  h_zz_EMu_Mass->SetLineColor(kOrange-8);
  h_zz_EMu_Mass->SetFillColor(kOrange-8);
  h_st_at_EMu_Mass->SetLineColor(kOrange-2);
  h_st_at_EMu_Mass->SetFillColor(kOrange-2);
  h_st_t_EMu_Mass->SetLineColor(2);
  h_st_t_EMu_Mass->SetFillColor(2);
   
  h_data_EMu_Mass->SetStats(0);
  h_dyll_EMu_Mass->SetStats(0);
  THStack *hs_EMu_Mass = new THStack("hs","");

  hs_EMu_Mass->Add(h_tt_EMu_Mass);
  hs_EMu_Mass->Add(h_st_at_EMu_Mass);
  hs_EMu_Mass->Add(h_st_t_EMu_Mass);
  hs_EMu_Mass->Add(h_ww_EMu_Mass);
  hs_EMu_Mass->Add(h_wz_EMu_Mass);
  hs_EMu_Mass->Add(h_zz_EMu_Mass);
  hs_EMu_Mass->Add(h_dyll_EMu_Mass);

  /*cout<<"Data: "<<h_data_EMu_Mass->Integral()<<endl;
  cout<<"DY to tau tau: "<<h_dyll_EMu_Mass->Integral()<<endl;
  cout<<"TT: "<<h_tt_EMu_Mass->Integral()<<endl;
  cout<<"ZZ: "<<h_zz_EMu_Mass->Integral()<<endl;
  cout<<"WZ: "<<h_wz_EMu_Mass->Integral()<<endl;
  cout<<"WW: "<<h_ww_EMu_Mass->Integral()<<endl;*/

  hs_EMu_Mass->Draw("hist");
  h_data_EMu_Mass->Draw("same p");

  hs_EMu_Mass->SetTitle("");
  hs_EMu_Mass->GetXaxis()->SetTitle("");
  hs_EMu_Mass->GetXaxis()->SetTickLength(0.03);
  hs_EMu_Mass->GetXaxis()->SetTitleOffset(1.05);
  hs_EMu_Mass->GetXaxis()->SetLabelSize(0.03);
  hs_EMu_Mass->GetXaxis()->SetLabelOffset(999);
  hs_EMu_Mass->GetXaxis()->SetRangeUser(15.,3000.);
  hs_EMu_Mass->GetXaxis()->SetMoreLogLabels();
  hs_EMu_Mass->GetXaxis()->SetNoExponent();
 
  hs_EMu_Mass->GetYaxis()->SetTitle("Number of Events");
  hs_EMu_Mass->GetYaxis()->SetLabelFont(42);
  hs_EMu_Mass->GetYaxis()->SetLabelSize(0.03);
  hs_EMu_Mass->GetYaxis()->SetTitleSize(0.04);
  hs_EMu_Mass->GetYaxis()->SetTitleOffset(1.15);
  hs_EMu_Mass->GetYaxis()->SetTitleFont(42);  
  hs_EMu_Mass->SetMinimum(0.3);
  hs_EMu_Mass->SetMaximum(10000);
  
  //hs->Print("all");  
  
  TLegend *leg = new TLegend(0.635783,0.6164841,0.8271252,0.8837705,NULL,"brNDC");
  leg->AddEntry(h_data_EMu_Mass,"Data","lp");
  leg->AddEntry(h_tt_EMu_Mass,"TT","f");
  leg->AddEntry(h_st_t_EMu_Mass,"tW","f");
  leg->AddEntry(h_st_at_EMu_Mass,"#bar{t}W","f");
  leg->AddEntry(h_dyll_EMu_Mass,"DY->#tau#tau","f");
  leg->AddEntry(h_ww_EMu_Mass,"WW","f");
  leg->AddEntry(h_wz_EMu_Mass,"WZ","f");
  leg->AddEntry(h_zz_EMu_Mass,"ZZ","f");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.02888889);
  leg->Draw();

  c1->cd();

  TH1F *h_mc_EMu_Mass = h_dyll_EMu_Mass->Clone("h_mc_EMu_Mass");
  h_mc_EMu_Mass->Add(h_tt_EMu_Mass);
  h_mc_EMu_Mass->Add(h_ww_EMu_Mass);
  h_mc_EMu_Mass->Add(h_wz_EMu_Mass);
  h_mc_EMu_Mass->Add(h_zz_EMu_Mass);
  h_mc_EMu_Mass->Add(h_st_at_EMu_Mass);
  h_mc_EMu_Mass->Add(h_st_t_EMu_Mass);
  cout<<"Total: "<<h_mc_EMu_Mass->Integral()<<endl;

  Int_t nbins = h_mc_EMu_Mass->GetSize();
  //cout<<"nbins: "<<nbins<<endl; 

  Double_t error;
  Double_t error1 = 0;
  for (Int_t i=1;i<nbins-1;i++) {
  error1 += (h_mc_EMu_Mass->GetBinError(i)) * (h_mc_EMu_Mass->GetBinError(i));
  }

  error = TMath::Sqrt(error1);
  cout<<"total error: "<<error<<endl;

  TH1F *hratio_EMu_Mass = h_data_EMu_Mass->Clone("hratio_EMu_Mass");
  hratio_EMu_Mass->Divide(h_mc_EMu_Mass);
     
  hratio_EMu_Mass->SetMarkerColor(kBlack);
  hratio_EMu_Mass->SetMarkerStyle(20);
  hratio_EMu_Mass->SetLineColor(kBlack);
  hratio_EMu_Mass->SetMarkerSize(0.4);

  hratio_EMu_Mass->SetTitle("  ");  
  hratio_EMu_Mass->GetXaxis()->SetTitle("e#mu mass GeV");
  hratio_EMu_Mass->GetXaxis()->SetLabelFont(42);
  hratio_EMu_Mass->GetXaxis()->SetLabelSize(0.15);
  hratio_EMu_Mass->GetXaxis()->SetTitleSize(0.14);
  hratio_EMu_Mass->GetXaxis()->SetTitleOffset(1);
  hratio_EMu_Mass->GetXaxis()->SetRangeUser(15.,3000.);
  hratio_EMu_Mass->GetXaxis()->SetMoreLogLabels();
  hratio_EMu_Mass->GetXaxis()->SetNoExponent();

  hratio_EMu_Mass->GetYaxis()->SetTitle("Data/MC");
  hratio_EMu_Mass->GetYaxis()->SetTitleSize(0.15);
  hratio_EMu_Mass->GetYaxis()->SetTitleOffset(0.3);
  hratio_EMu_Mass->GetYaxis()->SetLabelFont(42);
  hratio_EMu_Mass->GetYaxis()->SetLabelSize(0.1);
  hratio_EMu_Mass->GetYaxis()->SetRangeUser(0.,2.);
  hratio_EMu_Mass->GetYaxis()->SetNdivisions(5);

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

  hratio_EMu_Mass->Draw();
  TLine l1(15.,1.0,3000.,1.0);
  l1->SetLineWidth(1);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(5);
  l1->Draw("same");
  c1->Draw();
  c1->Update();
  
  c1->SaveAs("/home/ridhi/Work/Analysis/Data-MC/plots/13TeV/BKGEstimate/emu_Mass.png");
  //c1->SaveAs("/home/ridhi/Work/Analysis/Data-MC/plots/13TeV/BKGEstimate/emu_Mass.pdf");
  
  //c1->Write();
  h_dyll_EMu_Mass->Delete();
  hratio_EMu_Mass->Delete();

  file->Write();
  file->Close();
}
