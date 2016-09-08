//********************************* Unfolding Detector Resolution*********************

//TFile *file = TFile::Open("/home/ridhi/Work/Analysis/76X/Unfolding/Ratio_dataMC.root");

const int nBins=43;
Double_t xbin[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

//TH1D *h_scale = (TH1D*)file->Get("h_ratio");

int insideMassRange(double m)
{
  return ((m<xbin[0]) || (m>xbin[nBins])) ? 0 : 1;
}

void treeRead1() 
{
  const char *file1[2] = {"DYEE","DYEE_altSig"};
  //const char *file1[2] = {"DYEE_1_M10to3500","DYEE_2_M10to3500"};

  char name[200];
  double massReco, massGen, lumiWeight, genWeight, PUWeight;
  char isReco, isGen, BB, BE, EE;

  TH1D *reco_EEMass[2]; TH1D *gen_EEMass[2]; //TH1D *reco1_EEMass[2]; TH1D *gen1_EEMass[2];
  TH2D * resp1_DetRes[2]; TH2D * resp2_DetRes[2]; //TH2D * resp3_DetRes[2];
  TH2D * resp_DetRes_isMiss[2]; TH2D * resp_DetRes_isFake[2];

  TFile *file2[2];

  for(int i=0;i<2;i++){

    sprintf(name, "/home/ridhi/Work/Analysis/76X/Unfolding/detResponse/%s.root", file1[i]);
    //sprintf(name, "/home/ridhi/Work/Analysis/76X/Unfolding/unfolding_Checks/%s.root", file1[i]);
    TFile *f = TFile::Open(name);

    sprintf(name, "detResponse/%s_new.root", file1[i]);
    //sprintf(name, "unfolding_Checks/%s_new.root", file1[i]);
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

    //reco1_EEMass[i] = new TH1D("reco1_EEMass", "reco1_EEMass", 43.0, xbin);
    //gen1_EEMass[i]  = new TH1D("gen1_EEMass", "gen1_EEMass", 43.0, xbin);

    resp1_DetRes[i] = new TH2D("resp1_DetRes", "", 43.0, xbin, 43.0, xbin);
    resp2_DetRes[i] = new TH2D("resp2_DetRes", "", 43.0, xbin, 43.0, xbin);
    //resp3_DetRes[i] = new TH2D("resp3_DetRes", "", 43.0, xbin, 43.0, xbin);
    resp_DetRes_isMiss[i] = new TH2D("resp_DetRes_isMiss", "", 43.0, xbin, 43.0, xbin);
    resp_DetRes_isFake[i] = new TH2D("resp_DetRes_isFake", "", 43.0, xbin, 43.0, xbin);
    
    reco_EEMass[i]->Sumw2(); gen_EEMass[i]->Sumw2();
    resp1_DetRes[i]->Sumw2(); resp2_DetRes[i]->Sumw2(); resp_DetRes_isMiss[i]->Sumw2(); resp_DetRes_isFake[i]->Sumw2(); 

    for(Long64_t j=0;j<entries;j++){
      t->GetEntry(j);

      //double RatioValue = h_scale->GetBinContent(h_scale->GetXaxis()->FindFixBin(massReco));

      double scale = lumiWeight*genWeight*PUWeight*2316.969;

      if(massReco > -999.) reco_EEMass[i]->Fill(massReco,scale);
      if(massGen > -999.) gen_EEMass[i]->Fill(massGen,scale);

      //if(massReco > -999.) reco1_EEMass[i]->Fill(massReco,scale*RatioValue);
      //if(massGen > -999.) gen1_EEMass[i]->Fill(massGen,scale*RatioValue);

      //if(isReco && isGen) RespMatrix[i]->Fill(massReco,massGen,scale);
      //if(!isReco && isGen) isMiss_RespMatrix[i]->Fill(massReco,massGen,scale);
      //if(isReco && !isGen) isFake_RespMatrix[i]->Fill(massReco,massGen,scale);

      if(massReco > -999. && massGen > -999.) {
      	if (insideMassRange(massReco) || insideMassRange(massGen)) {
      		if (!insideMassRange(massReco) &&  insideMassRange(massGen)) resp_DetRes_isMiss[i]->Fill(massReco, massGen, scale);
      		else if ( insideMassRange(massReco) && !insideMassRange(massGen)) resp_DetRes_isFake[i]->Fill(massReco, massGen, scale);

		else (insideMassRange(massReco) && insideMassRange(massGen)) {
			resp1_DetRes[i]->Fill(massReco, massGen, scale);
			resp2_DetRes[i]->Fill(massGen, massReco, scale);
			//resp3_DetRes[i]->Fill(massReco, massGen, scale*RatioValue);
		}
        }
      }
    }

    file2[i]->Write();
    file2[i]->Close();
  }
}
