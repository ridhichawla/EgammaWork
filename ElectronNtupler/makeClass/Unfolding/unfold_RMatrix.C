#define unfold_RMatrix_cxx
#include "unfold_RMatrix.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

Double_t unfold_RMatrix::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

Double_t unfold_RMatrix::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t unfold_RMatrix::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

//bool unfold_RMatrix::mySortFnx(Double_t i, Double_t j) { return i>j; }

void unfold_RMatrix::Loop()
{

  TFile *file = new TFile("dyll_M-50.root", "recreate");
  
  int good_elec, gen_elec_pre, gen_elec_post;
  double dR;
  double sum_weights;
  double Z_Mass, Z_Y, Z_Rap, Z_Eta, Z_Pt, Z_Phi, ZMass_G_pre, ZMass_G_post, ZMass_R_pre, ZMass_R_post, ZMass_di_gPre, ZMass_di_gPost;
  double ZMass_R_post_barrel, ZMass_R_post_endcap, ZMass_R_post_b_e;
  
  TLorentzVector ele1,ele2,dielectron;
  TLorentzVector reco1_pre,reco2_pre,direco_pre,reco1_post,reco2_post,direco_post;
  
  TLorentzVector reco1_post_barrel,reco2_post_barrel,direco_post_barrel;
  TLorentzVector reco1_post_endcap,reco2_post_endcap,direco_post_endcap;
  TLorentzVector reco1_post_b_e,reco2_post_b_e,direco_post_b_e;

  TLorentzVector gen1_pre,gen2_pre,digen_pre,gen1_post,gen2_post,digen_post;
  TLorentzVector gPre1,gPre2,di_gPre;
  TLorentzVector gPost1,gPost2,di_gPost;

  vector <double> newscEt; vector <double> newscEta; vector <double> newscPhi; vector <double> newscEnr;
  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

  vector <double> newgPre_Pt; vector <double> newgPre_Eta; vector <double> newgPre_Enr; vector <double> newgPre_Phi;
  vector <double> newgPost_Pt; vector <double> newgPost_Eta; vector <double> newgPost_Enr; vector <double> newgPost_Phi;

  vector <double> recomatchedGen; vector <double> Recomatchedgen; vector <double> genmatchedReco; vector <double> Genmatchedreco;

  vector <double> recoPre_Pt; vector <double> recoPre_Eta; vector <double> recoPre_Enr; vector <double> recoPre_Phi;
  vector <double> genPre_Pt; vector <double> genPre_Eta; vector <double> genPre_Enr; vector <double> genPre_Phi;

  vector <double> recoPost_Pt; vector <double> recoPost_Eta; vector <double> recoPost_Enr; vector <double> recoPost_Phi;
  vector <double> genPost_Pt; vector <double> genPost_Eta; vector <double> genPost_Enr; vector <double> genPost_Phi;

  Double_t xbins[46] = {15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  TH1D *genPreFSR_Mass = new TH1D("genPreFSR_Mass", "genPreFSR_Mass", 45, xbins);
  TH2D *responsePre = new TH2D("RM_preFSR", "RM_preFSR", 45, xbins, 45, xbins);

  TH1D *genPostFSR_Mass = new TH1D("genPostFSR_Mass", "genPostFSR_Mass", 45, xbins);
  TH2D *responsePost = new TH2D("RM_postFSR", "RM_postFSR", 45, xbins, 45, xbins);
  TH2D *responsePost_barrel = new TH2D("RM_postFSR_BB", "RM_postFSR_BB", 45, xbins, 45, xbins);
  TH2D *responsePost_b_e = new TH2D("RM_postFSR_EB", "RM_postFSR_EB", 45, xbins, 45, xbins);
  TH2D *responsePost_endcap = new TH2D("RM_postFSR_EE", "RM_postFSR_EE", 45, xbins, 45, xbins);

  TH1D *ZMass_ = new TH1D("ZMass_", "ZMass_", 45, xbins);

  responsePre->Sumw2(); genPreFSR_Mass->Sumw2();
  responsePost->Sumw2(); genPostFSR_Mass->Sumw2();
  ZMass_->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 500;
  cout<<"entries: "<<nentries<<endl;

  sum_weights = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%500000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    int index[pt->size()];
    float ptnew[pt->size()];

    for(unsigned int r=0; r<pt->size(); r++)
    {
      ptnew[r]=pt->at(r);
    }

    int sizer = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(sizer,ptnew,index,true);

    /*int index1[gPre_pt->size()];
    float pt1[gPre_pt->size()];

    for(unsigned int a=0; a<gPre_pt->size(); a++)
    { 
      pt1[a]=gPre_pt->at(a);
    }

    int sizem = sizeof(pt1)/sizeof(pt1[0]);
    TMath::Sort(sizem,pt1,index1,true);*/

    int index2[gen_postFSR_pt->size()];
    float pt2[gen_postFSR_pt->size()];

    for(unsigned int b=0; b<gen_postFSR_pt->size(); b++)
    {
      pt2[b]=gen_postFSR_pt->at(b);
    }

    int sizen = sizeof(pt2)/sizeof(pt2[0]);
    TMath::Sort(sizen,pt2,index2,true);

    Z_Mass=0.; Z_Y=0.; Z_Rap=0.; Z_Eta=0.; Z_Pt=0.; Z_Phi=0.; ZMass_G_pre=0.; ZMass_R_pre=0.; ZMass_G_post=0.; ZMass_R_post=0.; ZMass_di_gPre=0.; ZMass_di_gPost=0.;
    ZMass_R_post_barrel=0.; ZMass_R_post_endcap=0.; ZMass_R_post_b_e=0.;
    good_elec = 0;
    gen_elec_pre = 0;
    gen_elec_post = 0;
    dR = 0;

    unsigned int matchGen1 = 0;
    unsigned int matchGen2 = 0;
    unsigned int matchReco1 = 0;
    unsigned int matchReco2 = 0;

    newscEt.clear(); newscEta.clear(); newscPhi.clear(); newscEnr.clear();
    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    newgPre_Pt.clear(); newgPre_Eta.clear(); newgPre_Enr.clear(); newgPre_Phi.clear(); 
    newgPost_Pt.clear(); newgPost_Eta.clear(); newgPost_Enr.clear(); newgPost_Phi.clear();

    recomatchedGen.clear(); Recomatchedgen.clear(); genmatchedReco.clear(); Genmatchedreco.clear();

    recoPre_Pt.clear(); recoPre_Eta.clear(); recoPre_Enr.clear(); recoPre_Phi.clear();
    genPre_Pt.clear();  genPre_Eta.clear(); genPre_Enr.clear(); genPre_Phi.clear();

    recoPost_Pt.clear(); recoPost_Eta.clear(); recoPost_Enr.clear(); recoPost_Phi.clear(); 
    genPost_Pt.clear(); genPost_Eta.clear(); genPost_Enr.clear(); genPost_Phi.clear();

    //sum_weights = sum_weights+theWeight;

    //cout<<"1"<<endl;

    if(!tauFlag){
      if(nEle>=2) {
	if(doubleElectron) {
	  if(doubleElectron == 0) cout<<"Trigger not applied correctly"<<endl;

	  for(int k=0;k<nEle;k++){

	    if(passMediumId->at(index[k]) == 1 && eleEcalDrivenSeed->at(index[k]) == 1){
	      if(fabs(eta->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){

		if(passMediumId->at(index[k]) == 0) cout<<"Wrong ID: "<<endl;
		if(eleEcalDrivenSeed->at(index[k]) == 0) cout<<"Not ECAL Driven: "<<endl;

		good_elec = good_elec + 1;

		newscEt.push_back(etSC->at(index[k]));
		newscEta.push_back(etaSC->at(index[k]));
		newscPhi.push_back(phiSC->at(index[k]));
		newscEnr.push_back(enSC->at(index[k]));
		newelePt.push_back(pt->at(index[k]));
		neweleEta.push_back(eta->at(index[k]));
		neweleEnr.push_back(energy->at(index[k]));
		newelePhi.push_back(phi->at(index[k]));
		neweleCharge.push_back(charge->at(index[k]));

	      } //kin
	    } // ID
	  } // k loop
	} //trigger
      } // nEle
    } // tau flag

    //cout<<"2"<<endl;

    if(good_elec==2){
      if(newelePt.at(0) < newelePt.at(1)) cout<<"event: "<<jentry<<"   "<<"Sorting not proper: "<<"   "<<"reco pt lead: "<<newelePt.at(0)<<"   "<<"reco pt sublead: "<<newelePt.at(1)<<endl;

      if(newelePt.at(0) > 25. && newelePt.at(1) > 20.){

	ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	dielectron=ele1+ele2;
	Z_Mass = dielectron.M();

	ZMass_->Fill(Z_Mass,theWeight);
      }
    }

    //cout<<"3"<<endl;
    /*
       if(gPre_pt->size()>0.){
       for(unsigned int l=0;l<gPre_eta->size();l++){

       if(fabs(gPre_eta->at(index1[l])) < 2.5 && !(fabs(gPre_eta->at(index1[l])) > 1.4442 && fabs(gPre_eta->at(index1[l])) < 1.566)){

       gen_elec_pre = gen_elec_pre + 1;

       newgPre_Pt.push_back(gPre_pt->at(index1[l]));
       newgPre_Eta.push_back(gPre_eta->at(index1[l]));
       newgPre_Enr.push_back(gPre_energy->at(index1[l]));
       newgPre_Phi.push_back(gPre_phi->at(index1[l]));
       }
       }
       }

    //cout<<"4"<<endl;

    if(gen_elec_pre>=2.){

    if(newgPre_Pt.at(0) < newgPre_Pt.at(1)) cout<<"event: "<<jentry<<"   "<<"Sorting not properly done"<<"   "<<"new gen pre pt lead: "<<newgPre_Pt.at(0)<<"   "<<"new gen pre pt sublead: "<<newgPre_Pt.at(1)<<endl;

    if(newgPre_Pt.at(0) > 20. && newgPre_Pt.at(1) > 10.){

    gPre1.SetPtEtaPhiE(newgPre_Pt.at(0),newgPre_Eta.at(0),newgPre_Phi.at(0),newgPre_Enr.at(0));
    gPre2.SetPtEtaPhiE(newgPre_Pt.at(1),newgPre_Eta.at(1),newgPre_Phi.at(1),newgPre_Enr.at(1));

    di_gPre=gPre1+gPre2;
    ZMass_di_gPre=di_gPre.M();
    genPreFSR_Mass->Fill(ZMass_di_gPre,theWeight);
    }
    }
    */
    //cout<<"5"<<endl;

    if(gen_postFSR_pt->size()>=2.){   
      for(unsigned int m=0;m<gen_postFSR_eta->size();m++){

	if(fabs(gen_postFSR_eta->at(index2[m])) < 2.5 && !(fabs(gen_postFSR_eta->at(index2[m])) > 1.4442 && fabs(gen_postFSR_eta->at(index2[m])) < 1.566)){

	  gen_elec_post = gen_elec_post + 1;

	  newgPost_Pt.push_back(gen_postFSR_pt->at(index2[m]));
	  newgPost_Eta.push_back(gen_postFSR_eta->at(index2[m]));
	  newgPost_Enr.push_back(gen_postFSR_ene->at(index2[m]));
	  newgPost_Phi.push_back(gen_postFSR_phi->at(index2[m]));
	}
      }
    }

    if(gen_elec_post==2){

      if(newgPost_Pt.at(0) < newgPost_Pt.at(1)) cout<<"entry: "<<jentry<<"   "<<"Sorting not done properly"<<"   "<<"new gen post pt lead: "<<newgPost_Pt.at(0)<<"   "<<"new gen post pt sublead: "<<newgPost_Pt.at(1)<<endl;

      if(newgPost_Pt.at(0) > 25. && newgPost_Pt.at(1) > 20.){

	gPost1.SetPtEtaPhiE(newgPost_Pt.at(0),newgPost_Eta.at(0),newgPost_Phi.at(0),newgPost_Enr.at(0));
	gPost2.SetPtEtaPhiE(newgPost_Pt.at(1),newgPost_Eta.at(1),newgPost_Phi.at(1),newgPost_Enr.at(1));

	di_gPost=gPost1+gPost2;
	ZMass_di_gPost=di_gPost.M();
	genPostFSR_Mass->Fill(ZMass_di_gPost,theWeight);
      }
    }

    //cout<<"7"<<endl;
    /*
       if(gen_elec_pre>=2){
       for(unsigned int ireco = 0; ireco < neweleEta.size(); ireco++){
       double dR_comp_pre = 1000.;
       for(unsigned int igen = 0; igen < newgPre_Eta.size(); igen++){
       dR = deltaR(neweleEta[ireco], newelePhi[ireco], newgPre_Eta[igen], newgPre_Phi[igen]);

       if(dR < 0.1  ){
       if (dR < dR_comp_pre)
       {
       dR_comp_pre = dR;
       matchGen1 = igen ; matchReco1 = ireco ; //iter[ireco]=igen;
       recomatchedGen.push_back(igen);
       Recomatchedgen.push_back(ireco);
       }

       genPre_Pt.push_back(newgPre_Pt[matchGen1]); genPre_Eta.push_back(newgPre_Eta[matchGen1]); genPre_Enr.push_back(newgPre_Enr[matchGen1]); genPre_Phi.push_back(newgPre_Phi[matchGen1]);
       recoPre_Pt.push_back(newelePt[matchReco1]); recoPre_Eta.push_back(neweleEta[matchReco1]); recoPre_Enr.push_back(neweleEnr[matchReco1]); recoPre_Phi.push_back(newelePhi[matchReco1]);
       }
       }
       }
       }

    //cout<<"8"<<endl;

    if(recoPre_Pt.size()>=2)
    {
    reco1_pre.SetPtEtaPhiE(recoPre_Pt.at(0),recoPre_Eta.at(0),recoPre_Phi.at(0),recoPre_Enr.at(0));
    reco2_pre.SetPtEtaPhiE(recoPre_Pt.at(1),recoPre_Eta.at(1),recoPre_Phi.at(1),recoPre_Enr.at(1));
    direco_pre=reco1_pre+reco2_pre;
    ZMass_R_pre = direco_pre.M();
    }

    //cout<<"9"<<endl;

    if(genPre_Pt.size()>=2)
    {
    gen1_pre.SetPtEtaPhiE(genPre_Pt.at(0),genPre_Eta.at(0),genPre_Phi.at(0),genPre_Enr.at(0));
    gen2_pre.SetPtEtaPhiE(genPre_Pt.at(1),genPre_Eta.at(1),genPre_Phi.at(1),genPre_Enr.at(1));
    digen_pre=gen1_pre+gen2_pre;
    ZMass_G_pre = digen_pre.M();
    }

    responsePre->Fill(ZMass_G_pre,ZMass_R_pre,theWeight);
    */

    if(gen_elec_post==2){
      for(unsigned int ireco = 0; ireco < neweleEta.size(); ireco++){
	double dR_comp_post = 1000.;
	for(unsigned int igen = 0; igen < newgPost_Eta.size(); igen++){
	  dR = deltaR(neweleEta[ireco], newelePhi[ireco], newgPost_Eta[igen], newgPost_Phi[igen]);

	  if(dR < 0.1  ){
	    if (dR < dR_comp_post)
	    {
	      dR_comp_post = dR;
	      matchGen2 = igen ; matchReco2 = ireco ; //iter[ireco]=igen;
	      recomatchedGen.push_back(igen);
	      Recomatchedgen.push_back(ireco);
	    }

	    genPost_Pt.push_back(newgPost_Pt[matchGen2]); genPost_Eta.push_back(newgPost_Eta[matchGen2]); genPost_Enr.push_back(newgPost_Enr[matchGen2]); genPost_Phi.push_back(newgPost_Phi[matchGen2]);
	    recoPost_Pt.push_back(newelePt[matchReco2]); recoPost_Eta.push_back(neweleEta[matchReco2]); recoPost_Enr.push_back(neweleEnr[matchReco2]); recoPost_Phi.push_back(newelePhi[matchReco2]);
	  }
	}
      }
    }

    if(recoPost_Pt.size()==2)
    {
      reco1_post.SetPtEtaPhiE(recoPost_Pt.at(0),recoPost_Eta.at(0),recoPost_Phi.at(0),recoPost_Enr.at(0));
      reco2_post.SetPtEtaPhiE(recoPost_Pt.at(1),recoPost_Eta.at(1),recoPost_Phi.at(1),recoPost_Enr.at(1));

      if(recoPost_Eta.at(0) < 1.4442 && recoPost_Eta.at(1) < 1.4442){
	reco1_post_barrel.SetPtEtaPhiE(recoPost_Pt.at(0),recoPost_Eta.at(0),recoPost_Phi.at(0),recoPost_Enr.at(0));
	reco2_post_barrel.SetPtEtaPhiE(recoPost_Pt.at(1),recoPost_Eta.at(1),recoPost_Phi.at(1),recoPost_Enr.at(1));
      }

      if((recoPost_Eta.at(0) > 1.566 && recoPost_Eta.at(0) < 2.5) && (recoPost_Eta.at(1) > 1.566 && recoPost_Eta.at(1) < 2.5)){
	reco1_post_endcap.SetPtEtaPhiE(recoPost_Pt.at(0),recoPost_Eta.at(0),recoPost_Phi.at(0),recoPost_Enr.at(0));
	reco2_post_endcap.SetPtEtaPhiE(recoPost_Pt.at(1),recoPost_Eta.at(1),recoPost_Phi.at(1),recoPost_Enr.at(1));
      }

      if((recoPost_Eta.at(0) < 1.4442 && (recoPost_Eta.at(1) > 1.566 && recoPost_Eta.at(1) < 2.5)) || ((recoPost_Eta.at(0) > 1.566 && recoPost_Eta.at(0) < 2.5) && recoPost_Eta.at(1) < 1.4442)) {
	reco1_post_b_e.SetPtEtaPhiE(recoPost_Pt.at(0),recoPost_Eta.at(0),recoPost_Phi.at(0),recoPost_Enr.at(0));
	reco2_post_b_e.SetPtEtaPhiE(recoPost_Pt.at(1),recoPost_Eta.at(1),recoPost_Phi.at(1),recoPost_Enr.at(1));
      }

      direco_post=reco1_post+reco2_post;
      ZMass_R_post = direco_post.M();

      direco_post_barrel=reco1_post_barrel+reco2_post_barrel;
      ZMass_R_post_barrel = direco_post_barrel.M();

      direco_post_endcap=reco1_post_endcap+reco2_post_endcap;
      ZMass_R_post_endcap = direco_post_endcap.M();

      direco_post_b_e=reco1_post_b_e+reco2_post_b_e;
      ZMass_R_post_b_e = direco_post_b_e.M();
    }

    if(genPost_Pt.size()==2)
    {
      gen1_post.SetPtEtaPhiE(genPost_Pt.at(0),genPost_Eta.at(0),genPost_Phi.at(0),genPost_Enr.at(0));
      gen2_post.SetPtEtaPhiE(genPost_Pt.at(1),genPost_Eta.at(1),genPost_Phi.at(1),genPost_Enr.at(1));
      digen_post=gen1_post+gen2_post;
      ZMass_G_post = digen_post.M();
    }

    responsePost->Fill(ZMass_G_post,ZMass_R_post,theWeight);
    responsePost_barrel->Fill(ZMass_G_post,ZMass_R_post_barrel,theWeight);
    responsePost_endcap->Fill(ZMass_G_post,ZMass_R_post_endcap,theWeight);
    responsePost_b_e->Fill(ZMass_G_post,ZMass_R_post_b_e,theWeight);

  } // event

  file->Write();
  file->Close();
}
