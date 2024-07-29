void dy() {

  TFile *f1 = TFile::Open("hist1.root");
  TFile *f2 = TFile::Open("hist2.root");

  TH1D *h1_num_pt_barrel = (TH1D*)f1->Get("numerator_pt_barrel");
  TH1D *h2_num_pt_barrel = (TH1D*)f2->Get("numerator_pt_barrel");

  TH1D *h1_num_pt_endcap = (TH1D*)f1->Get("numerator_pt_endcap");
  TH1D *h2_num_pt_endcap = (TH1D*)f2->Get("numerator_pt_endcap");

  TH1D *h1_den_pt_barrel = (TH1D*)f1->Get("denominator_pt_barrel");
  TH1D *h2_den_pt_barrel = (TH1D*)f2->Get("denominator_pt_barrel");

  TH1D *h1_den_pt_endcap = (TH1D*)f1->Get("denominator_pt_endcap");
  TH1D *h2_den_pt_endcap = (TH1D*)f2->Get("denominator_pt_endcap");

  TH1D *h1_num_barrel = (TH1D*)f1->Get("numerator_barrel");
  TH1D *h2_num_barrel = (TH1D*)f2->Get("numerator_barrel");

  TH1D *h1_num_endcap = (TH1D*)f1->Get("numerator_endcap");
  TH1D *h2_num_endcap = (TH1D*)f2->Get("numerator_endcap");

  TH1D *h1_den_barrel = (TH1D*)f1->Get("denominator_barrel");
  TH1D *h2_den_barrel = (TH1D*)f2->Get("denominator_barrel");

  TH1D *h1_den_endcap = (TH1D*)f1->Get("denominator_endcap");
  TH1D *h2_den_endcap = (TH1D*)f2->Get("denominator_endcap");

  //double lumi = 2258.066;
  double lumi = 1;

  double w1 = (18610*lumi)/2.30899e+12;
  double w2 = (2008.4*3*lumi)/4.5151e+11;

  h1_num_pt_barrel->Scale(w1);  h2_num_pt_barrel->Scale(w2);
  h1_num_pt_endcap->Scale(w1);  h2_num_pt_endcap->Scale(w2);
  h1_num_barrel->Scale(w1);  h2_num_barrel->Scale(w2);
  h1_num_endcap->Scale(w1);  h2_num_endcap->Scale(w2);

  h1_den_pt_barrel->Scale(w1);  h2_den_pt_barrel->Scale(w2);
  h1_den_pt_endcap->Scale(w1);  h2_den_pt_endcap->Scale(w2);
  h1_den_barrel->Scale(w1);  h2_den_barrel->Scale(w2);
  h1_den_endcap->Scale(w1);  h2_den_endcap->Scale(w2);

  TH1D *h_num_pt_barrel = (TH1D*)h1_num_pt_barrel->Clone("h_num_pt_barrel");
  h_num_pt_barrel->Add(h2_num_pt_barrel);

  TH1D *h_num_pt_endcap = (TH1D*)h1_num_pt_endcap->Clone("h_num_pt_endcap");
  h_num_pt_endcap->Add(h2_num_pt_endcap);

  TH1D *h_num_barrel = (TH1D*)h1_num_barrel->Clone("h_num_barrel");
  h_num_barrel->Add(h2_num_barrel);

  TH1D *h_num_endcap = (TH1D*)h1_num_endcap->Clone("h_num_endcap");
  h_num_endcap->Add(h2_num_endcap);

  TH1D *h_den_pt_barrel = (TH1D*)h1_den_pt_barrel->Clone("h_den_pt_barrel");
  h_den_pt_barrel->Add(h2_den_pt_barrel);

  TH1D *h_den_pt_endcap = (TH1D*)h1_den_pt_endcap->Clone("h_den_pt_endcap");
  h_den_pt_endcap->Add(h2_den_pt_endcap);

  TH1D *h_den_barrel = (TH1D*)h1_den_barrel->Clone("h_den_barrel");
  h_den_barrel->Add(h2_den_barrel);

  TH1D *h_den_endcap = (TH1D*)h1_den_endcap->Clone("h_den_endcap");
  h_den_endcap->Add(h2_den_endcap);

  TFile *file = TFile::Open("hist0.root", "RECREATE");

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
