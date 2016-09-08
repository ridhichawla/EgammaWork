//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 348 2014-08-08 22:18:23Z T.J.Adye@rl.ac.uk $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#endif


//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfold()
{

  TFile *f1 = TFile::Open("/home/ridhi/Work/Analysis/76X/Unfolding/detResponse/DYEE_new.root");
  TFile *f2 = TFile::Open("/home/ridhi/Work/Analysis/76X/Unfolding/detResponse/DYEE_altSig_new.root");

  TFile *f3 = TFile::Open("/home/ridhi/Work/Analysis/76X/Unfolding/FSRCorr/DYEE_new.root");

  TFile *file = new TFile ("/home/ridhi/Work/Analysis/76X/FinalCorr_diffXsec/RespObj_detUnfolding.root","RECREATE");

  TH2D *h2Resp1 = (TH2D*)f1->Get("resp1_DetRes");
  TH1D *hTrue1  = (TH1D*)f1->Get("gen_EEMass");
  TH1D *hMeas1  = (TH1D*)f1->Get("reco_EEMass");

  TH2D *h2Resp2 = (TH2D*)f2->Get("resp1_DetRes");
  TH1D *hTrue2  = (TH1D*)f2->Get("gen_EEMass");
  TH1D *hMeas2  = (TH1D*)f2->Get("reco_EEMass");

  /*TH2D *h2Resp3 = (TH2D*)f1->Get("resp3_DetRes");
  TH1D *hTrue3  = (TH1D*)f1->Get("gen1_EEMass");
  TH1D *hMeas3  = (TH1D*)f1->Get("reco1_EEMass");*/

  TH2D *h2Resp4 = (TH2D*)f3->Get("resp1_FSRCorr");
  TH1D *hTrue4  = (TH1D*)f3->Get("preFSR_Mass");
  TH1D *hMeas4  = (TH1D*)f3->Get("postFSR_Mass");

  RooUnfoldResponse response1 (hMeas1, hTrue1, h2Resp1); // default sample
  response1.UseOverflow(true);

  RooUnfoldResponse response2 (hMeas2, hTrue2, h2Resp2); // alternate sample
  response2.UseOverflow(true);

  //RooUnfoldResponse response3 (hMeas3, hTrue3, h2Resp3); // alternate: reweighted sample
  //response3.UseOverflow(true);

  RooUnfoldResponse response4 (hMeas4, hTrue4, h2Resp4); // FSR
  response4.UseOverflow(true);

  // ******************************** Fakes ************************************
  cout<<"Fakes1 : "<<response1.FakeEntries()<<endl;
  cout<<"Fakes2 : "<<response2.FakeEntries()<<endl;
  //cout<<"Fakes3 : "<<response3.FakeEntries()<<endl;
  cout<<"Fakes4 : "<<response4.FakeEntries()<<endl;

  response1.Write("Unfold_DetectorRes1");
  response2.Write("Unfold_DetectorRes2");
  //response3.Write("Unfold_DetectorRes3");
  response4.Write("Unfold_FSRCorr");

  file->Close();

}

#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
