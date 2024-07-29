#include "LeptonScaleLookup1.h"
#include <string.h>
#pragma once
LeptonScaleLookup::LeptonScaleLookup(std::string filename)
{
  file_ = new TFile(filename.c_str(), "READ");
  
  h2_eff_mc   = (TH2F*)file_->Get("hmeff");
  assert(h2_eff_mc);

  h2_eff_data   = (TH2F*)file_->Get("hdeff");
  assert(h2_eff_data);

}

LeptonScaleLookup::~LeptonScaleLookup()
{
  file_->Close();
  delete file_;
}

float LeptonScaleLookup::GetEfficiency(float eta, float pt, TH2F *hist) 
{
  // look up the efficiency
  const Int_t nPtBins = 6;
  Int_t binX, binY;
  
  binX = hist->GetXaxis()->FindFixBin(eta);

  // -- if Pt is larger than the largest Pt bin edge, SF is same with the value for the last bin -- // 
  if( pt >= 2000 ) binY = nPtBins-1;
  else binY = hist->GetYaxis()->FindFixBin(pt);
  //cout<<"binX: "<<binX<<"   "<<"binY: "<<binY<<endl;
  //cout<<""<<endl;

  return hist->GetBinContent(binX, binY);
}

float LeptonScaleLookup::GetError(float eta, float pt, TH2F *hist)
{
  
  // look up the error
  Int_t binX = hist->GetXaxis()->FindFixBin(eta);
  Int_t binY = hist->GetYaxis()->FindFixBin(pt);

  return hist->GetBinError(binX, binY);
}

float LeptonScaleLookup::GetExpectedTriggerEfficiency(float eta1, float pt1, float eta2, float pt2, string TrigType)
{
  
  float eff_mc_1, eff_mc_2, eff_data_1, eff_data_2, evt_eff;

  // get individual leg efficiencies
  string mc = "mc";
  string data = "data";  

  eff_mc_1 = GetEfficiency(eta1, pt1, h2_eff_mc);
  eff_mc_2 = GetEfficiency(eta2, pt2, h2_eff_mc);
  
  eff_data_1 = GetEfficiency(eta1, pt1, h2_eff_data);
  eff_data_2 = GetEfficiency(eta2, pt2, h2_eff_data);

  //cout<<"ptbin1: "<<eff_mc_1<<"   "<<

  if (mc.compare(TrigType)== 0)
  {
    evt_eff = 1 - (1-eff_mc_1)*(1-eff_mc_2);
  }

  if (data.compare(TrigType)== 0)
  {
    evt_eff = 1 - (1-eff_data_1)*(1-eff_data_2);
  }
  
  return evt_eff;
  
}


float LeptonScaleLookup::GetExpectedEfficiency(float eta1, float pt1, float eta2, float pt2, string TrigType)
{
  
  float eff_mc_1, eff_mc_2, eff_data_1, eff_data_2, evt_eff;
    
  // get individual efficiencies
  string mc   = "mc";
  string data = "data";  


  eff_mc_1 = GetEfficiency(eta1, pt1, h2_eff_mc);
  eff_mc_2 = GetEfficiency(eta2, pt2, h2_eff_mc);
  //cout<<"eff MC 1: "<<eff_mc_1<<"   "<<"eff MC 2: "<<eff_mc_2<<endl;
  
  eff_data_1 = GetEfficiency(eta1, pt1, h2_eff_data);
  eff_data_2 = GetEfficiency(eta2, pt2, h2_eff_data);
  //cout<<"eff Data 1: "<<eff_data_1<<"   "<<"eff Data 2: "<<eff_data_2<<endl;

  if (mc.compare(TrigType)== 0)
  {
    evt_eff = (eff_mc_1)*(eff_mc_2);
  }

  if(data.compare(TrigType)== 0)
  {
    evt_eff = (eff_data_1)*(eff_data_2);
  }  
  
  return evt_eff;
  
}
