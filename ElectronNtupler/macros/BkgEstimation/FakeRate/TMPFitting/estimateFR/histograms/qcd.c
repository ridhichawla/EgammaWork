void qcd() {

  TFile *f11 = TFile::Open("hist11.root");
  TFile *f12 = TFile::Open("hist12.root");
  TFile *f13 = TFile::Open("hist13.root");
  TFile *f14 = TFile::Open("hist14.root");
  TFile *f15 = TFile::Open("hist15.root");
  TFile *f16 = TFile::Open("hist16.root");
  TFile *f17 = TFile::Open("hist17.root");
  TFile *f18 = TFile::Open("hist18.root");

  TH1D *h1_num_pt_barrel = (TH1D*)f11->Get("numerator_pt_barrel");
  TH1D *h2_num_pt_barrel = (TH1D*)f12->Get("numerator_pt_barrel");
  TH1D *h3_num_pt_barrel = (TH1D*)f13->Get("numerator_pt_barrel");
  TH1D *h4_num_pt_barrel = (TH1D*)f14->Get("numerator_pt_barrel");
  TH1D *h5_num_pt_barrel = (TH1D*)f15->Get("numerator_pt_barrel");
  TH1D *h6_num_pt_barrel = (TH1D*)f16->Get("numerator_pt_barrel");
  TH1D *h7_num_pt_barrel = (TH1D*)f17->Get("numerator_pt_barrel");
  TH1D *h8_num_pt_barrel = (TH1D*)f18->Get("numerator_pt_barrel");

  TH1D *h1_num_pt_endcap = (TH1D*)f11->Get("numerator_pt_endcap");
  TH1D *h2_num_pt_endcap = (TH1D*)f12->Get("numerator_pt_endcap");
  TH1D *h3_num_pt_endcap = (TH1D*)f13->Get("numerator_pt_endcap");
  TH1D *h4_num_pt_endcap = (TH1D*)f14->Get("numerator_pt_endcap");
  TH1D *h5_num_pt_endcap = (TH1D*)f15->Get("numerator_pt_endcap");
  TH1D *h6_num_pt_endcap = (TH1D*)f16->Get("numerator_pt_endcap");
  TH1D *h7_num_pt_endcap = (TH1D*)f17->Get("numerator_pt_endcap");
  TH1D *h8_num_pt_endcap = (TH1D*)f18->Get("numerator_pt_endcap");

  TH1D *h1_den_pt_barrel = (TH1D*)f11->Get("denominator_pt_barrel");
  TH1D *h2_den_pt_barrel = (TH1D*)f12->Get("denominator_pt_barrel");
  TH1D *h3_den_pt_barrel = (TH1D*)f13->Get("denominator_pt_barrel");
  TH1D *h4_den_pt_barrel = (TH1D*)f14->Get("denominator_pt_barrel");
  TH1D *h5_den_pt_barrel = (TH1D*)f15->Get("denominator_pt_barrel");
  TH1D *h6_den_pt_barrel = (TH1D*)f16->Get("denominator_pt_barrel");
  TH1D *h7_den_pt_barrel = (TH1D*)f17->Get("denominator_pt_barrel");
  TH1D *h8_den_pt_barrel = (TH1D*)f18->Get("denominator_pt_barrel");

  TH1D *h1_den_pt_endcap = (TH1D*)f11->Get("denominator_pt_endcap");
  TH1D *h2_den_pt_endcap = (TH1D*)f12->Get("denominator_pt_endcap");
  TH1D *h3_den_pt_endcap = (TH1D*)f13->Get("denominator_pt_endcap");
  TH1D *h4_den_pt_endcap = (TH1D*)f14->Get("denominator_pt_endcap");
  TH1D *h5_den_pt_endcap = (TH1D*)f15->Get("denominator_pt_endcap");
  TH1D *h6_den_pt_endcap = (TH1D*)f16->Get("denominator_pt_endcap");
  TH1D *h7_den_pt_endcap = (TH1D*)f17->Get("denominator_pt_endcap");
  TH1D *h8_den_pt_endcap = (TH1D*)f18->Get("denominator_pt_endcap");

  TH1D *h1_num_barrel = (TH1D*)f11->Get("numerator_barrel");
  TH1D *h2_num_barrel = (TH1D*)f12->Get("numerator_barrel");
  TH1D *h3_num_barrel = (TH1D*)f13->Get("numerator_barrel");
  TH1D *h4_num_barrel = (TH1D*)f14->Get("numerator_barrel");
  TH1D *h5_num_barrel = (TH1D*)f15->Get("numerator_barrel");
  TH1D *h6_num_barrel = (TH1D*)f16->Get("numerator_barrel");
  TH1D *h7_num_barrel = (TH1D*)f17->Get("numerator_barrel");
  TH1D *h8_num_barrel = (TH1D*)f18->Get("numerator_barrel");

  TH1D *h1_num_endcap = (TH1D*)f11->Get("numerator_endcap");
  TH1D *h2_num_endcap = (TH1D*)f12->Get("numerator_endcap");
  TH1D *h3_num_endcap = (TH1D*)f13->Get("numerator_endcap");
  TH1D *h4_num_endcap = (TH1D*)f14->Get("numerator_endcap");
  TH1D *h5_num_endcap = (TH1D*)f15->Get("numerator_endcap");
  TH1D *h6_num_endcap = (TH1D*)f16->Get("numerator_endcap");
  TH1D *h7_num_endcap = (TH1D*)f17->Get("numerator_endcap");
  TH1D *h8_num_endcap = (TH1D*)f18->Get("numerator_endcap");

  TH1D *h1_den_barrel = (TH1D*)f11->Get("denominator_barrel");
  TH1D *h2_den_barrel = (TH1D*)f12->Get("denominator_barrel");
  TH1D *h3_den_barrel = (TH1D*)f13->Get("denominator_barrel");
  TH1D *h4_den_barrel = (TH1D*)f14->Get("denominator_barrel");
  TH1D *h5_den_barrel = (TH1D*)f15->Get("denominator_barrel");
  TH1D *h6_den_barrel = (TH1D*)f16->Get("denominator_barrel");
  TH1D *h7_den_barrel = (TH1D*)f17->Get("denominator_barrel");
  TH1D *h8_den_barrel = (TH1D*)f18->Get("denominator_barrel");

  TH1D *h1_den_endcap = (TH1D*)f11->Get("denominator_endcap");
  TH1D *h2_den_endcap = (TH1D*)f12->Get("denominator_endcap");
  TH1D *h3_den_endcap = (TH1D*)f13->Get("denominator_endcap");
  TH1D *h4_den_endcap = (TH1D*)f14->Get("denominator_endcap");
  TH1D *h5_den_endcap = (TH1D*)f15->Get("denominator_endcap");
  TH1D *h6_den_endcap = (TH1D*)f16->Get("denominator_endcap");
  TH1D *h7_den_endcap = (TH1D*)f17->Get("denominator_endcap");
  TH1D *h8_den_endcap = (TH1D*)f18->Get("denominator_endcap");

  //double lumi = 2258.066;
  double lumi = 1;

  double w1 = (2526000*lumi)/5471093;
  double w2 = (4833200*lumi)/9260312;
  double w3 = (6850000*lumi)/4693590;
  double w4 = (1900000*lumi)/22463790;
  double w5 = (478520*lumi)/36029601;
  double w6 = (68592*lumi)/36202369;
  double w7 = (20859*lumi)/11518861;
  double w8 = (1350*lumi)/7340288;

  h1_num_pt_barrel->Scale(w1);  h2_num_pt_barrel->Scale(w2);  h3_num_pt_barrel->Scale(w3);  h4_num_pt_barrel->Scale(w4);
  h5_num_pt_barrel->Scale(w5);  h6_num_pt_barrel->Scale(w6);  h7_num_pt_barrel->Scale(w7);  h8_num_pt_barrel->Scale(w8);

  h1_num_pt_endcap->Scale(w1);  h2_num_pt_endcap->Scale(w2);  h3_num_pt_endcap->Scale(w3);  h4_num_pt_endcap->Scale(w4);
  h5_num_pt_endcap->Scale(w5);  h6_num_pt_endcap->Scale(w6);  h7_num_pt_endcap->Scale(w7);  h8_num_pt_endcap->Scale(w8);

  h1_num_barrel->Scale(w1);  h2_num_barrel->Scale(w2);  h3_num_barrel->Scale(w3);  h4_num_barrel->Scale(w4);
  h5_num_barrel->Scale(w5);  h6_num_barrel->Scale(w6);  h7_num_barrel->Scale(w7);  h8_num_barrel->Scale(w8);

  h1_num_endcap->Scale(w1);  h2_num_endcap->Scale(w2);  h3_num_endcap->Scale(w3);  h4_num_endcap->Scale(w4);
  h5_num_endcap->Scale(w5);  h6_num_endcap->Scale(w6);  h7_num_endcap->Scale(w7);  h8_num_endcap->Scale(w8);

  h1_den_pt_barrel->Scale(w1);  h2_den_pt_barrel->Scale(w2);  h3_den_pt_barrel->Scale(w3);  h4_den_pt_barrel->Scale(w4);
  h5_den_pt_barrel->Scale(w5);  h6_den_pt_barrel->Scale(w6);  h7_den_pt_barrel->Scale(w7);  h8_den_pt_barrel->Scale(w8);

  h1_den_pt_endcap->Scale(w1);  h2_den_pt_endcap->Scale(w2);  h3_den_pt_endcap->Scale(w3);  h4_den_pt_endcap->Scale(w4);
  h5_den_pt_endcap->Scale(w5);  h6_den_pt_endcap->Scale(w6);  h7_den_pt_endcap->Scale(w7);  h8_den_pt_endcap->Scale(w8);

  h1_den_barrel->Scale(w1);  h2_den_barrel->Scale(w2);  h3_den_barrel->Scale(w3);  h4_den_barrel->Scale(w4);
  h5_den_barrel->Scale(w5);  h6_den_barrel->Scale(w6);  h7_den_barrel->Scale(w7);  h8_den_barrel->Scale(w8);

  h1_den_endcap->Scale(w1);  h2_den_endcap->Scale(w2);  h3_den_endcap->Scale(w3);  h4_den_endcap->Scale(w4);
  h5_den_endcap->Scale(w5);  h6_den_endcap->Scale(w6);  h7_den_endcap->Scale(w7);  h8_den_endcap->Scale(w8);

  TH1D *h_num_pt_barrel = (TH1D*)h1_num_pt_barrel->Clone("h_num_pt_barrel");
  h_num_pt_barrel->Add(h2_num_pt_barrel);
  h_num_pt_barrel->Add(h3_num_pt_barrel);
  h_num_pt_barrel->Add(h4_num_pt_barrel);
  h_num_pt_barrel->Add(h5_num_pt_barrel);
  h_num_pt_barrel->Add(h6_num_pt_barrel);
  h_num_pt_barrel->Add(h7_num_pt_barrel);
  h_num_pt_barrel->Add(h8_num_pt_barrel);

  TH1D *h_num_pt_endcap = (TH1D*)h1_num_pt_endcap->Clone("h_num_pt_endcap");
  h_num_pt_endcap->Add(h2_num_pt_endcap);
  h_num_pt_endcap->Add(h3_num_pt_endcap);
  h_num_pt_endcap->Add(h4_num_pt_endcap);
  h_num_pt_endcap->Add(h5_num_pt_endcap);
  h_num_pt_endcap->Add(h6_num_pt_endcap);
  h_num_pt_endcap->Add(h7_num_pt_endcap);
  h_num_pt_endcap->Add(h8_num_pt_endcap);

  TH1D *h_num_barrel = (TH1D*)h1_num_barrel->Clone("h_num_barrel");
  h_num_barrel->Add(h2_num_barrel);
  h_num_barrel->Add(h3_num_barrel);
  h_num_barrel->Add(h4_num_barrel);
  h_num_barrel->Add(h5_num_barrel);
  h_num_barrel->Add(h6_num_barrel);
  h_num_barrel->Add(h7_num_barrel);
  h_num_barrel->Add(h8_num_barrel);

  TH1D *h_num_endcap = (TH1D*)h1_num_endcap->Clone("h_num_endcap");
  h_num_endcap->Add(h2_num_endcap);
  h_num_endcap->Add(h3_num_endcap);
  h_num_endcap->Add(h4_num_endcap);
  h_num_endcap->Add(h5_num_endcap);
  h_num_endcap->Add(h6_num_endcap);
  h_num_endcap->Add(h7_num_endcap);
  h_num_endcap->Add(h8_num_endcap);

  TH1D *h_den_pt_barrel = (TH1D*)h1_den_pt_barrel->Clone("h_den_pt_barrel");
  h_den_pt_barrel->Add(h2_den_pt_barrel);
  h_den_pt_barrel->Add(h3_den_pt_barrel);
  h_den_pt_barrel->Add(h4_den_pt_barrel);
  h_den_pt_barrel->Add(h5_den_pt_barrel);
  h_den_pt_barrel->Add(h6_den_pt_barrel);
  h_den_pt_barrel->Add(h7_den_pt_barrel);
  h_den_pt_barrel->Add(h8_den_pt_barrel);

  TH1D *h_den_pt_endcap = (TH1D*)h1_den_pt_endcap->Clone("h_den_pt_endcap");
  h_den_pt_endcap->Add(h2_den_pt_endcap);
  h_den_pt_endcap->Add(h3_den_pt_endcap);
  h_den_pt_endcap->Add(h4_den_pt_endcap);
  h_den_pt_endcap->Add(h5_den_pt_endcap);
  h_den_pt_endcap->Add(h6_den_pt_endcap);
  h_den_pt_endcap->Add(h7_den_pt_endcap);
  h_den_pt_endcap->Add(h8_den_pt_endcap);

  TH1D *h_den_barrel = (TH1D*)h1_den_barrel->Clone("h_den_barrel");
  h_den_barrel->Add(h2_den_barrel);
  h_den_barrel->Add(h3_den_barrel);
  h_den_barrel->Add(h4_den_barrel);
  h_den_barrel->Add(h5_den_barrel);
  h_den_barrel->Add(h6_den_barrel);
  h_den_barrel->Add(h7_den_barrel);
  h_den_barrel->Add(h8_den_barrel);

  TH1D *h_den_endcap = (TH1D*)h1_den_endcap->Clone("h_den_endcap");
  h_den_endcap->Add(h2_den_endcap);
  h_den_endcap->Add(h3_den_endcap);
  h_den_endcap->Add(h4_den_endcap);
  h_den_endcap->Add(h5_den_endcap);
  h_den_endcap->Add(h6_den_endcap);
  h_den_endcap->Add(h7_den_endcap);
  h_den_endcap->Add(h8_den_endcap);

  TFile *file = TFile::Open("hist5.root", "RECREATE");

  h_num_pt_barrel->Write("numerator_pt_barrel");
  h_num_pt_endcap->Write("numerator_pt_endcap");
  h_num_barrel->Write("numerator_barrel");
  h_num_endcap->Write("numerator_endcap");
  h_den_pt_barrel->Write("denominator_pt_barrel");
  h_den_pt_endcap->Write("denominator_pt_endcap");
  h_den_barrel->Write("denominator_barrel");
  h_den_endcap->Write("denominator_endcap");

  cout<<h_num_barrel->Integral()<<endl;

  //file->Write();
  file->Close();

}
