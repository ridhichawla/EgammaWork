#define cutFlow_cxx
#include "cutFlow.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <math.h>
//#include "DataFormats/Math/interface/deltaR.h"

void cutFlow::Loop()
{
  TFile *file = new TFile("output.root", "recreate");
  TH1F * wts = new TH1F("wts", "wts", 30000, -30000, 30000);

  int kin, id;
  int n_ele, n_trig, n_kin, n_id, n_all;// n_charge;
  bool isCharge, isMass;
  int n_dR;
  //double ZMass, ZY, ZRap, ZEta, ZPt, ZPhi;
  //double deltaR1, deltaR2;
  TLorentzVector ele1,ele2,dielectron,ele12,ele22,dielectron2;

  std::ofstream out;
  out.open ("Output.txt");

  vector <int> indx;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //int nentries = 50;
  cout<<"entries: "<<nentries<<endl;

  n_ele = 0;
  n_trig = 0;
  n_kin = 0;
  n_id = 0;
  n_charge = 0;
  n_all = 0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //cout<<"Weights: "<<theWeight<<endl;
    if(theWeight>0) wts->Fill(theWeight);
    
    kin = 0;
    id = 0;
    isCharge = false;
    isMass = false;

    indx.clear();

    if(nEle>=2){
      //++n_ele;

      if(sEleMC){
	//if(!doubleElectron) continue;
	//++n_trig;

	for(int k=0;k<nEle;k++){
	  //cout<<"Before: "<<"   "<<"pt: "<<etSC->at(k)<<"   "<<"eta: "<<etaSC->at(k)<<endl;
	  /*
	     if(jentry==451114)
	     cout<<"Event: "<<jentry<<"   "<<"Electrons: "<<nEle<<"   "<<"pt: "<<pt->at(k)<<"   "<<"eta: "<<eta->at(k)<<"   "<<"etSC: "<<etSC->at(k)<<"   "<<"etaSC: "<<etaSC->at(k)<<endl;
	     if(jentry==630371)
	     cout<<"Event: "<<jentry<<"   "<<"Electrons: "<<nEle<<"   "<<"pt: "<<pt->at(k)<<"   "<<"eta: "<<eta->at(k)<<"   "<<"etSC: "<<etSC->at(k)<<"   "<<"etaSC: "<<etaSC->at(k)<<endl;
	     if(jentry==788438)
	     cout<<"Event: "<<jentry<<"   "<<"Electrons: "<<nEle<<"   "<<"pt: "<<pt->at(k)<<"   "<<"eta: "<<eta->at(k)<<"   "<<"etSC: "<<etSC->at(k)<<"   "<<"etaSC: "<<etaSC->at(k)<<endl;
	     if(jentry==1256064)
	     cout<<"Event: "<<jentry<<"   "<<"Electrons: "<<nEle<<"   "<<"pt: "<<pt->at(k)<<"   "<<"eta: "<<eta->at(k)<<"   "<<"etSC: "<<etSC->at(k)<<"   "<<"etaSC: "<<etaSC->at(k)<<endl;
	     if(jentry==1860409)
	     cout<<"Event: "<<jentry<<"   "<<"Electrons: "<<nEle<<"   "<<"pt: "<<pt->at(k)<<"   "<<"eta: "<<eta->at(k)<<"   "<<"etSC: "<<etSC->at(k)<<"   "<<"etaSC: "<<etaSC->at(k)<<endl;
	     */

	  //if(etSC->at(k) < 10.) continue;
	  //if(etaSC->at(k) > 2.7) continue;
	  //if(etSC->at(k) < 10.) cout<<"Wrong pt: "<<endl;
	  //if(etaSC->at(k) > 2.7) cout<<"Wrong eta: "<<endl;
	  //cout<<"After: "<<"   "<<"pt: "<<pt->at(k)<<"   "<<"eta: "<<etaSC->at(k)<<endl;

	  //isKin_ID = false;
	  //if(fabs(etaSC->at(k)) < 2.5 && etSC->at(k) > 25. && !(fabs(etaSC->at(k)) > 1.4442 && fabs(etaSC->at(k)) < 1.566))
	  if(fabs(eta->at(k)) < 2.5 && pt->at(k) > 25. && !(fabs(etaSC->at(k)) > 1.4442 && fabs(etaSC->at(k)) < 1.566)){
	    //if(fabs(etaSC->at(k)) > 2.5) cout<<"Wrong eta: "<<endl;
	    //if(etSC->at(k) < 25) cout<<"Wrong pT: "<<endl;
	    //if((fabs(etaSC->at(k)) > 1.4442 && fabs(etaSC->at(k)) < 1.566)) entry<<"Gap region included: "<<endl;

	    //if(pt->at(k) < 10.) cout<<"Event: "<<jentry<<"   "<<"pT less than 10"<<"   "<<pt->at(k)<<"   "<<"etSC: "<<etSC->at(k)<<endl;
	    //if(fabs(eta->at(k)) > 2.7) cout<<"Event: "<<jentry<<"   "<<"eta greater than 2.7"<<"   "<<eta->at(k)<<"   "<<"etaSC: "<<etaSC->at(k)<<endl;
	    kin++;
	    if(passMediumId->at(k) == 1 && eleEcalDrivenSeed->at(k) == 1){
	      //if(passMediumId->at(k) == 0) entry<<"Wrong ID: "<<endl;
	      //if(eleEcalDrivenSeed->at(k) == 0) entry<<"Wrong ECAL ID: "<<endl;
	      id++;
	      indx.push_back(k);
	      /*if(id==2){
		ele1.SetPtEtaPhiE((*etSC)[indx[0]],(*etaSC)[indx[0]],(*phiSC)[indx[0]],(*enSC)[indx[0]]);
		ele2.SetPtEtaPhiE((*etSC)[indx[1]],(*etaSC)[indx[1]],(*phiSC)[indx[1]],(*enSC)[indx[1]]);

		dielectron=ele1+ele2;
		ZMass = dielectron.M();
		z_mass_ID->Fill(ZMass);
		}*/

	      //if(id==2 && (*charge)[indx[0]] * (*charge)[indx[1]]== -1){
		//isCharge = true;

		if(indx.size()>2) cout<<"Event: "<<jentry<<"   "<<"No. of electrons exceeds 2"<<endl;
		//cout<<"Entry: "<<jentry<<"   "<<"Electrons: "<<nEle<<endl;
		//newelePt.push_back((*pt)[indx]);
		//cout<<"Event: "<<jentry<<"   "<<"indx size: "<<indx.size()<<endl;
		//for(unsigned int i=0;i<indx.size();i++)

		ele1.SetPtEtaPhiE((*etSC)[indx[0]],(*etaSC)[indx[0]],(*phiSC)[indx[0]],(*enSC)[indx[0]]);
		ele2.SetPtEtaPhiE((*etSC)[indx[1]],(*etaSC)[indx[1]],(*phiSC)[indx[1]],(*enSC)[indx[1]]);

		dielectron=ele1+ele2;
		ZMass = dielectron.M();

		if(ZMass >= 60 && ZMass <= 120){
		  isMass = true;
		} // mass

	      //} //charge
	    } // id
	  } // kinematic
	} // electron
      } // trigger
    } // nEle

    if(nEle>=2){
      n_ele++;
      if(sEleMC){
	n_trig++;
	if(kin>=2){
	  n_kin++;
	  if(id>=2){
	    n_id++;
	    //if(isCharge){
	      //n_charge++;
	      if(isMass){
		n_all++;
	      } //isMass
	    //} // isCharge
	  } // id
	} // kin
      } //trigger
    } // nEle

  } // event

  cout<<"n_ele: "<<n_ele<<"   "<<"n_trig: "<<n_trig<<"   "<<"n_kin: "<<n_kin<<"   "<<"n_id: "<<n_id<<"   "<<"n_charge: "<<n_charge<<"   "<<"n_all: "<<n_all<<endl;
  file->Write();
  file->Close();
}
