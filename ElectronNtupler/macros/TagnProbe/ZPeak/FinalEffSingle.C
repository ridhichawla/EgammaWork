#include "LeptonScaleLookup1.cc"

void FinalEffSingle()
{
  gStyle->SetPaintTextFormat("1.3f");
  TFile *f = new TFile("../input_eleVar/DY_forEff_M10to3000.root");
  TTree *t;
  t = (TTree*)f->Get("tree");

  TFile *file = TFile::Open("ROOTFile_EffSF.root","recreate");
	
  double Ele1PT, Ele1Eta, Ele2PT, Ele2Eta,ZMass, genWeights,lumiWeights;

  t->SetBranchAddress("Ele1PT",&Ele1PT);
  t->SetBranchAddress("Ele1Eta",&Ele1Eta);
  t->SetBranchAddress("Ele2PT",&Ele2PT);
  t->SetBranchAddress("Ele2Eta",&Ele2Eta);
  t->SetBranchAddress("ZMass",&ZMass);
  t->SetBranchAddress("genWeights",&genWeights);
  t->SetBranchAddress("lumiWeights",&lumiWeights);

  Double_t xbins[32] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

  double nMassBin = 31;

  TH1D *hist_mass_num = new TH1D("hist_mass_num","Post FSR Mass",nMassBin,xbins);
  TH1D *hist_mass_den = new TH1D("hist_mass_den","Post FSR Mass",nMassBin,xbins);
  
  hist_mass_num->Sumw2();
  hist_mass_den->Sumw2();

  float eff_mc_RECO, eff_mc_ID, eff_mc_Trig;
  float eff_data_RECO, eff_data_ID, eff_data_Trig;
  float eff_mc, eff_data, eff_SF;
  int binX;
  int binY;
  int nentries = (int)t->GetEntries();

  cout << "Number of Entries = " << nentries << endl;

  string mc   = "mc";
  string data = "data";


  for(int i=0; i<nentries; i++)
  {
    if(i%100000==0)cout<<"OUT:"<<i<<endl;
    LeptonScaleLookup a("RECO_Stat.root");
    LeptonScaleLookup b("Medium_Stat.root");
    LeptonScaleLookup c("Trig_Stat.root");
    t->GetEntry(i);

    eff_mc_RECO   = a.GetExpectedEfficiency(Ele1Eta,Ele1PT,Ele2Eta,Ele2PT,mc);
    eff_data_RECO = a.GetExpectedEfficiency(Ele1Eta,Ele1PT,Ele2Eta,Ele2PT,data);

    eff_mc_ID     = b.GetExpectedEfficiency(Ele1Eta,Ele1PT,Ele2Eta,Ele2PT,mc);
    eff_data_ID   = b.GetExpectedEfficiency(Ele1Eta,Ele1PT,Ele2Eta,Ele2PT,data);

    eff_mc_Trig   = c.GetExpectedTriggerEfficiency(Ele1Eta,Ele1PT,Ele2Eta,Ele2PT,mc);
    eff_data_Trig = c.GetExpectedTriggerEfficiency(Ele1Eta,Ele1PT,Ele2Eta,Ele2PT,data);

    eff_mc   = eff_mc_RECO*eff_mc_ID*eff_mc_Trig;
    eff_data = eff_data_RECO*eff_data_ID*eff_data_Trig;
    eff_SF   = eff_data/eff_mc;
    if(eff_data == 0 || eff_mc == 0) cout<<"entry: "<<i<<"   "<<eff_data<<"   "<<eff_mc<<endl;

    /*cout<<"entry: "<<i<<endl;
    cout<<"ZMass: "<<ZMass<<endl;
    cout<<"pt1: "<<Ele1PT<<"   "<<"eta1: "<<Ele1Eta<<endl;
    cout<<"pt2: "<<Ele2PT<<"   "<<"eta2: "<<Ele2Eta<<endl;
    cout<<"RECO       eff MC: " <<eff_mc_RECO<<endl;//"   "<<"eff Data: "<<eff_data_RECO<<endl;
    cout<<"ID         eff MC: " <<eff_mc_ID<<"   "<<"eff Data: "<<eff_data_ID<<endl;
    cout<<"Trig       eff MC: " <<eff_mc_Trig<<"   "<<"eff Data: "<<eff_data_Trig<<endl;
    cout<<"Total eff MC: "<<eff_mc<<"   "<<"Total eff Data: "<<eff_data<<endl;
    cout<<"Eff SF CV: "<<eff_SF<<endl;
    cout<<""<<endl;*/

    hist_mass_num->Fill(ZMass,genWeights*lumiWeights*eff_SF);
    hist_mass_den->Fill(ZMass,genWeights*lumiWeights);

  }

  file->Write();
  file->Close();

}
