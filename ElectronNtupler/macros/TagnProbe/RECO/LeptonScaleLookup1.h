#ifndef LEPTONSCALELOOKUP_H
#define LEPTONSCALELOOKUP_H

#include "TFile.h"
#include <iostream>
#include <string>
#include "TH2F.h"
#include <cmath>

using namespace std;

class LeptonScaleLookup {

 public:
  
  LeptonScaleLookup(std::string filename);
  ~LeptonScaleLookup();
  
  float GetEfficiency(float eta, float pt, TH2F *hist);
  float GetError(float eta, float pt, TH2F *hist);
  
  float GetExpectedTriggerEfficiency(float eta1, float pt1, float eta2, float pt2,string TrigType);
  float GetEfficiencyReference(float pt1, float pt2);
  float GetExpectedError(float eta1, float pt1, float eta2, float pt2);
  void printTable(std::string name);
  float GetExpectedEfficiency(float eta1, float pt1, float eta2, float pt2,string TrigType);
  
 private:
  TFile *file_;
  TH2F *h2_eff_data;
  TH2F *h2_eff_mc ;
  
};

void LoopupAll(std::string file);

#endif
