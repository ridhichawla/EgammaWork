//********************************* Unfolding Detector Resolution*********************

const int nBins=43;  
Double_t xbin[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

int insideMassRange(double m)
{
  return ((m<xbin[0]) || (m>xbin[nBins])) ? 0 : 1;
}

void treeRead1() 
{
  const char *file1[2] = {"DYEE_final","DY_AltSig_M10to3500"};

  char name[200];
  double massReco, massGen, lumiWeight, genWeight, PUWeight;
  char isReco, isGen, BB, BE, EE;

  //const int nBins=43;
  //Double_t xbin[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
  //Double_t xbin[32] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

  TH1D *reco_EEMass[2]; TH1D *gen_EEMass[2]; TH2D *RespMatrix[2]; TH2D *isMiss_RespMatrix[2]; TH2D *isFake_RespMatrix[2];
  TH2D * resp_DetRes[2];

  TFile *file2[2];

  for(int i=0;i<2;i++){

    sprintf(name, "/home/ridhi/Work/Analysis/76X/Tuples/13TeV/Unfolding/detRes_Unfold/%s.root", file1[i]);
    TFile *f = TFile::Open(name);

    sprintf(name, "detRes_Unfold/%s_new_test.root", file1[i]);
    //sprintf(name, "detRes_Unfold/ZPeak_xsection/%s_new.root", file1[i]);
    file2[i] = new TFile(name,"RECREATE");
 
    TTree *t;
    t = (TTree*)f->Get("tree");
    Long64_t entries = t->GetEntries();
    cout<<"entry: "<<entries<<endl;

    t->SetBranchAddress("massReco",&massReco);
    t->SetBranchAddress("massGen",&massGen);
    t->SetBranchAddress("lumiWeight",&lumiWeight);
    t->SetBranchAddress("genWeight",&genWeight);
    t->SetBranchAddress("PUWeight",&PUWeight);
    t->SetBranchAddress("isReco",&isReco);
    t->SetBranchAddress("isGen",&isGen);
    t->SetBranchAddress("BB",&BB);
    t->SetBranchAddress("BE",&BE);
    t->SetBranchAddress("EE",&EE);

    sprintf(name, "%s",file1[i]);
    reco_EEMass[i] = new TH1D("reco_EEMass", "reco_EEMass", 43.0, xbin);
    gen_EEMass[i]  = new TH1D("gen_EEMass", "gen_EEMass", 43.0, xbin);
    RespMatrix[i]  = new TH2D("RespMatrix", "RespMatrix", 43.0, xbin, 43.0, xbin);
    isMiss_RespMatrix[i] = new TH2D("isMiss_RespMatrix", "isMass_RespMatrix", 43.0, xbin, 43.0, xbin);
    isFake_RespMatrix[i] = new TH2D("isFake_RespMatrix", "isFake_RespMatrix", 43.0, xbin, 43.0, xbin);

    resp_DetRes[i] = new TH2D("resp_DetRes", "", 43.0, xbin, 43.0, xbin);
    
    reco_EEMass[i]->Sumw2(); gen_EEMass[i]->Sumw2(); RespMatrix[i]->Sumw2(); isMiss_RespMatrix[i]->Sumw2(); isFake_RespMatrix[i]->Sumw2();

    for(Long64_t j=0;j<entries;j++){
      t->GetEntry(j);

      double scale = lumiWeight*genWeight*PUWeight*2316.969;

      if(massReco > -999.) reco_EEMass[i]->Fill(massReco,scale);
      if(massGen > -999.) gen_EEMass[i]->Fill(massGen,scale);

      if(isReco && isGen) RespMatrix[i]->Fill(massReco,massGen,scale);
      if(!isReco && isGen) isMiss_RespMatrix[i]->Fill(massReco,massGen,scale);
      if(isReco && !isGen) isFake_RespMatrix[i]->Fill(massReco,massGen,scale);

      if(massReco > -999. && massGen > -999.) {
      	if (insideMassRange(massReco) || insideMassRange(massGen)) {
      		//if (!insideMassRange(massReco) &&  insideMassRange(massGen)) resp_DetRes->Fill(massReco, massGen, scale);
      		//else if ( insideMassRange(MassReco) && !insideMassRange(massGen)) resp_DetRes->Fill(massReco, massGen, scale);

		if(insideMassRange(massReco) && insideMassRange(massGen)) resp_DetRes[i]->Fill(massReco, massGen, scale);
        }
      }
    }

    file2[i]->Write();
    file2[i]->Close();
  }
}
