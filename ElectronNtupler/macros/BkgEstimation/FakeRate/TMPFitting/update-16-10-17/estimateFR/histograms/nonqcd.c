void nonqcd() {

  TFile *f1 = TFile::Open("hist1.root");
  TFile *f2 = TFile::Open("hist2.root");
  TFile *f3 = TFile::Open("hist3.root");
  TFile *f4 = TFile::Open("hist4.root");
  TFile *f7 = TFile::Open("hist7.root");
  TFile *f8 = TFile::Open("hist8.root");
  TFile *f9 = TFile::Open("hist9.root");

  TH1D *h1_num_pt_barrel = (TH1D*)f1->Get("numerator_pt_barrel");
  TH1D *h2_num_pt_barrel = (TH1D*)f2->Get("numerator_pt_barrel");
  TH1D *h3_num_pt_barrel = (TH1D*)f3->Get("numerator_pt_barrel");
  TH1D *h4_num_pt_barrel = (TH1D*)f4->Get("numerator_pt_barrel");
  TH1D *h7_num_pt_barrel = (TH1D*)f7->Get("numerator_pt_barrel");
  TH1D *h8_num_pt_barrel = (TH1D*)f8->Get("numerator_pt_barrel");
  TH1D *h9_num_pt_barrel = (TH1D*)f9->Get("numerator_pt_barrel");

  TH1D *h1_num_pt_endcap = (TH1D*)f1->Get("numerator_pt_endcap");
  TH1D *h2_num_pt_endcap = (TH1D*)f2->Get("numerator_pt_endcap");
  TH1D *h3_num_pt_endcap = (TH1D*)f3->Get("numerator_pt_endcap");
  TH1D *h4_num_pt_endcap = (TH1D*)f4->Get("numerator_pt_endcap");
  TH1D *h7_num_pt_endcap = (TH1D*)f7->Get("numerator_pt_endcap");
  TH1D *h8_num_pt_endcap = (TH1D*)f8->Get("numerator_pt_endcap");
  TH1D *h9_num_pt_endcap = (TH1D*)f9->Get("numerator_pt_endcap");

  TH1D *h1_den_pt_barrel = (TH1D*)f1->Get("denominator_pt_barrel");
  TH1D *h2_den_pt_barrel = (TH1D*)f2->Get("denominator_pt_barrel");
  TH1D *h3_den_pt_barrel = (TH1D*)f3->Get("denominator_pt_barrel");
  TH1D *h4_den_pt_barrel = (TH1D*)f4->Get("denominator_pt_barrel");
  TH1D *h7_den_pt_barrel = (TH1D*)f7->Get("denominator_pt_barrel");
  TH1D *h8_den_pt_barrel = (TH1D*)f8->Get("denominator_pt_barrel");
  TH1D *h9_den_pt_barrel = (TH1D*)f9->Get("denominator_pt_barrel");

  TH1D *h1_den_pt_endcap = (TH1D*)f1->Get("denominator_pt_endcap");
  TH1D *h2_den_pt_endcap = (TH1D*)f2->Get("denominator_pt_endcap");
  TH1D *h3_den_pt_endcap = (TH1D*)f3->Get("denominator_pt_endcap");
  TH1D *h4_den_pt_endcap = (TH1D*)f4->Get("denominator_pt_endcap");
  TH1D *h7_den_pt_endcap = (TH1D*)f7->Get("denominator_pt_endcap");
  TH1D *h8_den_pt_endcap = (TH1D*)f8->Get("denominator_pt_endcap");
  TH1D *h9_den_pt_endcap = (TH1D*)f9->Get("denominator_pt_endcap");

  TH1D *h1_num_barrel = (TH1D*)f1->Get("numerator_barrel");
  TH1D *h2_num_barrel = (TH1D*)f2->Get("numerator_barrel");
  TH1D *h3_num_barrel = (TH1D*)f3->Get("numerator_barrel");
  TH1D *h4_num_barrel = (TH1D*)f4->Get("numerator_barrel");
  TH1D *h7_num_barrel = (TH1D*)f7->Get("numerator_barrel");
  TH1D *h8_num_barrel = (TH1D*)f8->Get("numerator_barrel");
  TH1D *h9_num_barrel = (TH1D*)f9->Get("numerator_barrel");

  TH1D *h1_num_endcap = (TH1D*)f1->Get("numerator_endcap");
  TH1D *h2_num_endcap = (TH1D*)f2->Get("numerator_endcap");
  TH1D *h3_num_endcap = (TH1D*)f3->Get("numerator_endcap");
  TH1D *h4_num_endcap = (TH1D*)f4->Get("numerator_endcap");
  TH1D *h7_num_endcap = (TH1D*)f7->Get("numerator_endcap");
  TH1D *h8_num_endcap = (TH1D*)f8->Get("numerator_endcap");
  TH1D *h9_num_endcap = (TH1D*)f9->Get("numerator_endcap");

  TH1D *h1_den_barrel = (TH1D*)f1->Get("denominator_barrel");
  TH1D *h2_den_barrel = (TH1D*)f2->Get("denominator_barrel");
  TH1D *h3_den_barrel = (TH1D*)f3->Get("denominator_barrel");
  TH1D *h4_den_barrel = (TH1D*)f4->Get("denominator_barrel");
  TH1D *h7_den_barrel = (TH1D*)f7->Get("denominator_barrel");
  TH1D *h8_den_barrel = (TH1D*)f8->Get("denominator_barrel");
  TH1D *h9_den_barrel = (TH1D*)f9->Get("denominator_barrel");

  TH1D *h1_den_endcap = (TH1D*)f1->Get("denominator_endcap");
  TH1D *h2_den_endcap = (TH1D*)f2->Get("denominator_endcap");
  TH1D *h3_den_endcap = (TH1D*)f3->Get("denominator_endcap");
  TH1D *h4_den_endcap = (TH1D*)f4->Get("denominator_endcap");
  TH1D *h7_den_endcap = (TH1D*)f7->Get("denominator_endcap");
  TH1D *h8_den_endcap = (TH1D*)f8->Get("denominator_endcap");
  TH1D *h9_den_endcap = (TH1D*)f9->Get("denominator_endcap");

  //double lumi = 2258.066;
  double lumi = 1;

  double w1 = (18610*lumi)/2.30899e+12;
  double w2 = (2008.4*3*lumi)/4.5151e+11;
  double w3 = (831.76*lumi)/85849572;
  double w4 = (61526.7*lumi)/3.73193e+12;
  double w7 = (118.7*lumi)/988416;
  double w8 = (47.13*lumi)/999996;
  double w9 = (16.523*lumi)/985598;

  h1_num_pt_barrel->Scale(w1);
  h2_num_pt_barrel->Scale(w2);
  h3_num_pt_barrel->Scale(w3);
  h4_num_pt_barrel->Scale(w4);
  h7_num_pt_barrel->Scale(w7);
  h8_num_pt_barrel->Scale(w8);
  h9_num_pt_barrel->Scale(w9);
  
  h1_num_pt_endcap->Scale(w1);
  h2_num_pt_endcap->Scale(w2);
  h3_num_pt_endcap->Scale(w3);
  h4_num_pt_endcap->Scale(w4);
  h7_num_pt_endcap->Scale(w7);
  h8_num_pt_endcap->Scale(w8);
  h9_num_pt_endcap->Scale(w9);
  
  h1_num_barrel->Scale(w1);
  h2_num_barrel->Scale(w2);
  h3_num_barrel->Scale(w3);
  h4_num_barrel->Scale(w4);
  h7_num_barrel->Scale(w7);
  h8_num_barrel->Scale(w8);
  h9_num_barrel->Scale(w9);
  
  h1_num_endcap->Scale(w1);
  h2_num_endcap->Scale(w2);
  h3_num_endcap->Scale(w3);
  h4_num_endcap->Scale(w4);
  h7_num_endcap->Scale(w7);
  h8_num_endcap->Scale(w8);
  h9_num_endcap->Scale(w9);

  h1_den_pt_barrel->Scale(w1);
  h2_den_pt_barrel->Scale(w2);
  h3_den_pt_barrel->Scale(w3);
  h4_den_pt_barrel->Scale(w4);
  h7_den_pt_barrel->Scale(w7);
  h8_den_pt_barrel->Scale(w8);
  h9_den_pt_barrel->Scale(w9);
  
  h1_den_pt_endcap->Scale(w1);
  h2_den_pt_endcap->Scale(w2);
  h3_den_pt_endcap->Scale(w3);
  h4_den_pt_endcap->Scale(w4);
  h7_den_pt_endcap->Scale(w7);
  h8_den_pt_endcap->Scale(w8);
  h9_den_pt_endcap->Scale(w9);
  
  h1_den_barrel->Scale(w1);
  h2_den_barrel->Scale(w2);
  h3_den_barrel->Scale(w3);
  h4_den_barrel->Scale(w4);
  h7_den_barrel->Scale(w7);
  h8_den_barrel->Scale(w8);
  h9_den_barrel->Scale(w9);
  
  h1_den_endcap->Scale(w1);
  h2_den_endcap->Scale(w2);
  h3_den_endcap->Scale(w3);
  h4_den_endcap->Scale(w4);
  h7_den_endcap->Scale(w7);
  h8_den_endcap->Scale(w8);
  h9_den_endcap->Scale(w9);

  TH1D *h_num_pt_barrel = (TH1D*)h1_num_pt_barrel->Clone("h_num_pt_barrel");
  h_num_pt_barrel->Add(h2_num_pt_barrel);
  h_num_pt_barrel->Add(h3_num_pt_barrel);
  h_num_pt_barrel->Add(h4_num_pt_barrel);
  h_num_pt_barrel->Add(h7_num_pt_barrel);
  h_num_pt_barrel->Add(h8_num_pt_barrel);
  h_num_pt_barrel->Add(h9_num_pt_barrel);

  TH1D *h_num_pt_endcap = (TH1D*)h1_num_pt_endcap->Clone("h_num_pt_endcap");
  h_num_pt_endcap->Add(h2_num_pt_endcap);
  h_num_pt_endcap->Add(h3_num_pt_endcap);
  h_num_pt_endcap->Add(h4_num_pt_endcap);
  h_num_pt_endcap->Add(h7_num_pt_endcap);
  h_num_pt_endcap->Add(h8_num_pt_endcap);
  h_num_pt_endcap->Add(h9_num_pt_endcap);

  TH1D *h_num_barrel = (TH1D*)h1_num_barrel->Clone("h_num_barrel");
  h_num_barrel->Add(h2_num_barrel);
  h_num_barrel->Add(h3_num_barrel);
  h_num_barrel->Add(h4_num_barrel);
  h_num_barrel->Add(h7_num_barrel);
  h_num_barrel->Add(h8_num_barrel);
  h_num_barrel->Add(h9_num_barrel);

  TH1D *h_num_endcap = (TH1D*)h1_num_endcap->Clone("h_num_endcap");
  h_num_endcap->Add(h2_num_endcap);
  h_num_endcap->Add(h3_num_endcap);
  h_num_endcap->Add(h4_num_endcap);
  h_num_endcap->Add(h7_num_endcap);
  h_num_endcap->Add(h8_num_endcap);
  h_num_endcap->Add(h9_num_endcap);

  TH1D *h_den_pt_barrel = (TH1D*)h1_den_pt_barrel->Clone("h_den_pt_barrel");
  h_den_pt_barrel->Add(h2_den_pt_barrel);
  h_den_pt_barrel->Add(h3_den_pt_barrel);
  h_den_pt_barrel->Add(h4_den_pt_barrel);
  h_den_pt_barrel->Add(h7_den_pt_barrel);
  h_den_pt_barrel->Add(h8_den_pt_barrel);
  h_den_pt_barrel->Add(h9_den_pt_barrel);

  TH1D *h_den_pt_endcap = (TH1D*)h1_den_pt_endcap->Clone("h_den_pt_endcap");
  h_den_pt_endcap->Add(h2_den_pt_endcap);
  h_den_pt_endcap->Add(h3_den_pt_endcap);
  h_den_pt_endcap->Add(h4_den_pt_endcap);
  h_den_pt_endcap->Add(h7_den_pt_endcap);
  h_den_pt_endcap->Add(h8_den_pt_endcap);
  h_den_pt_endcap->Add(h9_den_pt_endcap);

  TH1D *h_den_barrel = (TH1D*)h1_den_barrel->Clone("h_den_barrel");
  h_den_barrel->Add(h2_den_barrel);
  h_den_barrel->Add(h3_den_barrel);
  h_den_barrel->Add(h4_den_barrel);
  h_den_barrel->Add(h7_den_barrel);
  h_den_barrel->Add(h8_den_barrel);
  h_den_barrel->Add(h9_den_barrel);

  TH1D *h_den_endcap = (TH1D*)h1_den_endcap->Clone("h_den_endcap");
  h_den_endcap->Add(h2_den_endcap);
  h_den_endcap->Add(h3_den_endcap);
  h_den_endcap->Add(h4_den_endcap);
  h_den_endcap->Add(h7_den_endcap);
  h_den_endcap->Add(h8_den_endcap);
  h_den_endcap->Add(h9_den_endcap);

  TFile *file = TFile::Open("nonqcd.root", "RECREATE");

  h_num_pt_barrel->Write("numerator_pt_barrel");
  h_num_pt_endcap->Write("numerator_pt_endcap");
  h_num_barrel->Write("numerator_barrel");
  h_num_endcap->Write("numerator_endcap");
  h_den_pt_barrel->Write("denominator_pt_barrel");
  h_den_pt_endcap->Write("denominator_pt_endcap");
  h_den_barrel->Write("denominator_barrel");
  h_den_endcap->Write("denominator_endcap");

  //cout<<h_num_barrel->Integral()<<endl;

  //file->Write();
  file->Close();

}
