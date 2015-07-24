// -*- C++ -*-
//
// Package:    EgammaWork/ElectronNtupler
// Class:      SimpleElectronNtupler
// 
/**\class SimpleElectronNtupler SimpleElectronNtupler.cc EgammaWork/ElectronNtupler/plugins/SimpleElectronNtupler.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefHolder.h"
#include "DataFormats/Common/interface/RefVectorHolder.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "/afs/cern.ch/work/r/rchawla/private/CMSSW_7_4_0/src/EgammaWork/ElectronNtupler/plugins/utils.h"
//#include "FWCore/ParameterSet/interface/FileInPath.h"
//#include "/afs/cern.ch/work/r/rchawla/private/CMSSW_7_4_0/src/EgammaWork/ElectronNtupler/plugins/FileInPath.h"
//
// class declaration
//

class SimpleElectronNtupler : public edm::EDAnalyzer {
  public:
    explicit SimpleElectronNtupler(const edm::ParameterSet&);
    ~SimpleElectronNtupler();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0, 
      TRUE_PROMPT_ELECTRON, 
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    int matchToTruth(const edm::Ptr<reco::GsfElectron> el,  
	const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);
    void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);

    // ----------member data ---------------------------
    // Data members that are the same for AOD and miniAOD
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

    // AOD case data members
    edm::EDGetToken electronsToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

    // MiniAOD case data members
    edm::EDGetToken electronsMiniAODToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;

    // ID decisions objects
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

    // MVA values and categories (optional)
    edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

    edm::Service<TFileService> fs;
    TTree *electronTree_;

    // If MC
    bool misMC;

    // Histograms
    TH1F* nevents_;

    // General
    Double_t RunNo_, EvtNo_, Lumi_, Bunch_;

    // Vars for PVs
    Int_t pvNTracks_;

    // Vars for pile-up
    Int_t nPUTrue_;    // true pile-up
    Int_t nPU_;        // generated pile-up
    Int_t nPV_;        // number of reconsrtucted primary vertices
    Float_t rho_;      // the rho variable

    // Trigger
    std::string string_singleEle;
    std::string string_doubleEle;
    std::vector<std::string> hlNames_;
    std::vector<std::string> single_electron_triggers_in_run;
    std::vector<std::string> double_electron_triggers_in_run;
    bool singleElectron;
    bool doubleElectron;

    std::vector<int> idx_singleElectron;
    std::vector<int> idx_doubleElectron;

    // all electron variables
    Int_t nElectrons_;
    Int_t nGenElectrons_;

    std::vector<Float_t> ZMass_;
    std::vector<Float_t> ZPt_;
    std::vector<Float_t> ZEta_;
    std::vector<Float_t> ZRap_;
    std::vector<Float_t> ZPhi_;

    std::vector<Float_t> gPre_energy_;
    std::vector<Float_t> gPre_px_;
    std::vector<Float_t> gPre_py_;
    std::vector<Float_t> gPre_pz_;
    std::vector<Float_t> gPre_pt_;
    std::vector<Float_t> gPre_eta_;
    std::vector<Float_t> gPre_rap_;
    std::vector<Float_t> gPre_phi_;

    std::vector<Float_t> gPost_energy_;
    std::vector<Float_t> gPost_px_;
    std::vector<Float_t> gPost_py_;
    std::vector<Float_t> gPost_pz_;
    std::vector<Float_t> gPost_pt_;
    std::vector<Float_t> gPost_eta_;
    std::vector<Float_t> gPost_rap_;
    std::vector<Float_t> gPost_phi_;

    std::vector<Float_t> pt_;
    std::vector<Float_t> eta_;
    std::vector<Float_t> rap_;
    std::vector<Float_t> phi_;
    std::vector<Float_t> energy_;
    std::vector<Float_t> mass_;
    std::vector<Float_t> charge_;

    std::vector<Float_t> enSC_;
    std::vector<Float_t> preEnSC_;
    std::vector<Float_t> rawEnSC_;
    std::vector<Float_t> etSC_;
    std::vector<Float_t> etaSC_;
    std::vector<Float_t> phiSC_;

    std::vector<Float_t> full5x5_sigmaIetaIeta_;
    std::vector<Float_t> E1x5_;
    std::vector<Float_t> E2x5_;
    std::vector<Float_t> E5x5_;
    std::vector<Float_t> hOverE_;
    std::vector<Float_t> etaScWidth_;
    std::vector<Float_t> phiScWidth_;
    std::vector<Float_t> r9_;

    std::vector<Float_t> dEtaIn_;
    std::vector<Float_t> dPhiIn_;
    std::vector<Float_t> isoChargedHadrons_;
    std::vector<Float_t> isoNeutralHadrons_;
    std::vector<Float_t> isoPhotons_;
    std::vector<Float_t> isoChargedFromPU_;
    std::vector<Float_t> isoDeltaBeta_;
    std::vector<Float_t> isoRho_;
    std::vector<Float_t> ooEmooP_;
    std::vector<Float_t> d0_;
    std::vector<Float_t> dz_;
    std::vector<Int_t>   expectedMissingInnerHits_;
    std::vector<Int_t>   passConversionVeto_;
    std::vector<Float_t> brem_;
    std::vector<Int_t>   isTrue_;

    std::vector<Float_t> eleInBarrel_;
    std::vector<Float_t> eleInEndcap_;

    std::vector<Int_t> passVetoId_;
    std::vector<Int_t> passLooseId_;
    std::vector<Int_t> passMediumId_;
    std::vector<Int_t> passTightId_;

    std::vector<Float_t> mvaValue_;
    std::vector<Int_t>   mvaCategory_;

    std::vector<double> pt_leg1;
    std::vector<double> eta_leg1;
    std::vector<double> phi_leg1;

    std::vector<double> pt_leg2;
    std::vector<double> eta_leg2;
    std::vector<double> phi_leg2;

    std::vector<double> pt_SE;
    std::vector<double> eta_SE;
    std::vector<double> phi_SE;

    //std::vector<pat::TriggerObjectStandAlone> leg1triggerObj;
    //std::vector<pat::TriggerObjectStandAlone> leg2triggerObj;

    //std::vector<std::string> tagFilterName_;
    //std::vector<std::string> probeFilterName_;

    //std::vector<bool> passleg1Trigger;
    //std::vector<bool> passleg2Trigger;

    double DeltaR(const pat::Electron& e, std::vector<pat::TriggerObjectStandAlone> object);
};

SimpleElectronNtupler::SimpleElectronNtupler(const edm::ParameterSet& iConfig):
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{
  nevents_ = fs->make<TH1F>("nevents_","nevents_",2,0,2);
  string_singleEle = "HLT_Ele27_eta2p1_WPLoose_Gsf_v";
  string_doubleEle = "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  misMC            = iConfig.getUntrackedParameter<bool>("isMC");

  // Prepare tokens for all input collections and objects

  // Trigger
  triggerToken_  = mayConsume<edm::TriggerResults>
    (iConfig.getParameter<edm::InputTag>
     ("trigger"));

  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>
    (iConfig.getParameter<edm::InputTag>
     ("objects"));

  /*tagFilterName_ = (iConfig.getParameter<std::vector<std::string>>
    ("tagFilterName"));

    probeFilterName_ = (iConfig.getParameter<std::vector<std::string>>
    ("probeFilterName"));

    if (tagFilterName_.size() != probeFilterName_.size()) {
    std::cout << "Need to specify the same numbers of tag and probe filters." << std::endl;
    abort();
    }*/

  // Universal tokens for AOD and miniAOD  
  pileupToken_ = consumes<edm::View<PileupSummaryInfo> >
    (iConfig.getParameter <edm::InputTag>
     ("pileup"));

  rhoToken_    = consumes<double> 
    (iConfig.getParameter <edm::InputTag>
     ("rho"));

  beamSpotToken_    = consumes<reco::BeamSpot> 
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));

  // AOD tokens
  electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));

  vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));

  // MiniAOD tokens
  // For electrons, use the fact that pat::Electron can be cast into GsfElectron
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

  vtxMiniAODToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("verticesMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

  conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversionsMiniAOD"));

  //
  // Set up the ntuple structure
  //
  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

  electronTree_->Branch("RunNo", &RunNo_, "RunNo/D");
  electronTree_->Branch("EvtNo", &EvtNo_, "EvtNo/D");
  electronTree_->Branch("Lumi", &Lumi_, "Lumi/D");
  electronTree_->Branch("Bunch", &Bunch_, "Bunch/D");

  electronTree_->Branch("pvNTracks"    ,  &pvNTracks_ , "pvNTracks/I");

  electronTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  electronTree_->Branch("nPU"        ,  &nPU_     , "nPU/I");
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue_ , "nPUTrue/I");
  electronTree_->Branch("rho"        ,  &rho_ , "rho/F");

  electronTree_->Branch("singleElectron"    ,  &singleElectron    );
  electronTree_->Branch("doubleElectron"    ,  &doubleElectron    );

  electronTree_->Branch("pt_leg1"    ,  &pt_leg1    );
  electronTree_->Branch("eta_leg1"   ,  &eta_leg1   );
  electronTree_->Branch("phi_leg1"   ,  &phi_leg1   );
  electronTree_->Branch("pt_leg2"    ,  &pt_leg2    );
  electronTree_->Branch("eta_leg2"   ,  &eta_leg2   );
  electronTree_->Branch("phi_leg2"   ,  &phi_leg2   );
  electronTree_->Branch("pt_SE"    ,  &pt_SE    );
  electronTree_->Branch("eta_SE"   ,  &eta_SE   );
  electronTree_->Branch("phi_SE"   ,  &phi_SE   );

  electronTree_->Branch("nEle"    ,  &nElectrons_ , "nEle/I");
  electronTree_->Branch("nGenEle"    ,  &nGenElectrons_ , "nGenEle/I");
  electronTree_->Branch("ZMass"    ,  &ZMass_    );
  electronTree_->Branch("ZPt"    ,  &ZPt_    );
  electronTree_->Branch("ZEta"    ,  &ZEta_    );
  electronTree_->Branch("ZRap"    ,  &ZRap_    );
  electronTree_->Branch("ZPhi"    ,  &ZPhi_    );
  electronTree_->Branch("gPre_energy"    ,  &gPre_energy_    );
  electronTree_->Branch("gPre_px"    ,  &gPre_px_    );
  electronTree_->Branch("gPre_py"    ,  &gPre_py_    );
  electronTree_->Branch("gPre_pz"    ,  &gPre_pz_    );
  electronTree_->Branch("gPre_pt"    ,  &gPre_pt_    );
  electronTree_->Branch("gPre_eta"    ,  &gPre_eta_    );
  electronTree_->Branch("gPre_rap"    ,  &gPre_rap_    );
  electronTree_->Branch("gPre_phi"    ,  &gPre_phi_    );
  electronTree_->Branch("gPost_energy"    ,  &gPost_energy_    );
  electronTree_->Branch("gPost_px"    ,  &gPost_px_    );
  electronTree_->Branch("gPost_py"    ,  &gPost_py_    );
  electronTree_->Branch("gPost_pz"    ,  &gPost_pz_    );
  electronTree_->Branch("gPost_pt"    ,  &gPost_pt_    );
  electronTree_->Branch("gPost_eta"    ,  &gPost_eta_    );
  electronTree_->Branch("gPost_rap"    ,  &gPost_rap_    );
  electronTree_->Branch("gPost_phi"    ,  &gPost_phi_    );
  electronTree_->Branch("pt"    ,  &pt_    );
  electronTree_->Branch("eta"    ,  &eta_    );
  electronTree_->Branch("rap"    ,  &rap_    );
  electronTree_->Branch("phi"    ,  &phi_    );
  electronTree_->Branch("energy"    ,  &energy_    );
  electronTree_->Branch("mass"    ,  &mass_    );
  electronTree_->Branch("charge"    ,  &charge_    );
  electronTree_->Branch("enSC" ,  &enSC_ );
  electronTree_->Branch("preEnSC" ,  &preEnSC_ );
  electronTree_->Branch("rawEnSC" ,  &rawEnSC_ );
  electronTree_->Branch("etSC" ,  &etSC_ );
  electronTree_->Branch("etaSC" ,  &etaSC_ );
  electronTree_->Branch("phiSC" ,  &phiSC_ );
  electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_);
  electronTree_->Branch("E1x5" ,  &E1x5_ );
  electronTree_->Branch("E2x5" ,  &E2x5_ );
  electronTree_->Branch("E2x5" ,  &E2x5_ );
  electronTree_->Branch("hOverE",  &hOverE_);
  electronTree_->Branch("etaScWidth" ,  &etaScWidth_ );
  electronTree_->Branch("phiScWidth" ,  &phiScWidth_ );
  electronTree_->Branch("r9",  &r9_);

  electronTree_->Branch("dEtaIn",  &dEtaIn_);
  electronTree_->Branch("dPhiIn",  &dPhiIn_);
  electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  electronTree_->Branch("isoPhotons"             , &isoPhotons_);
  electronTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  electronTree_->Branch("isoDeltaBeta"     , &isoDeltaBeta_);
  electronTree_->Branch("isoRho"     , &isoRho_);
  electronTree_->Branch("ooEmooP", &ooEmooP_);
  electronTree_->Branch("d0"     , &d0_);
  electronTree_->Branch("dz"     , &dz_);
  electronTree_->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_);
  electronTree_->Branch("passConversionVeto", &passConversionVeto_);
  electronTree_->Branch("brem",  &brem_);
  electronTree_->Branch("isTrue"    , &isTrue_);

  electronTree_->Branch("passVetoId"  ,  &passVetoId_ );
  electronTree_->Branch("passLooseId"  ,  &passLooseId_ );
  electronTree_->Branch("passMediumId" ,  &passMediumId_ );
  electronTree_->Branch("passTightId"  ,  &passTightId_ );
  electronTree_->Branch("mvaVal" ,  &mvaValue_ );
  electronTree_->Branch("mvaCat" ,  &mvaCategory_ );

}


SimpleElectronNtupler::~SimpleElectronNtupler()
{
}

// member functions

// ------------ method called for each event  ------------
  void
SimpleElectronNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  //bool trigResult = false;

  nevents_->Fill(1);
  RunNo_ = iEvent.id().run();
  EvtNo_ = iEvent.id().event();
  Lumi_  = iEvent.luminosityBlock();
  Bunch_ = iEvent.bunchCrossing();

  //cout<<"1"<<endl;

  Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByLabel("generator", genEvtInfo);

  std::vector<double> evtWeights = genEvtInfo->weights();
  double theWeight = genEvtInfo->weight();

  Handle<LHEEventProduct> EvtHandle;
  iEvent.getByLabel("externalLHEProducer", EvtHandle);

  int whichWeight = 1;
  theWeight *= EvtHandle->weights()[whichWeight].wgt/EvtHandle->originalXWGTUP(); 
  
  cout<<"Final Weights: "<<theWeight<<endl;
  
  
  // Get Triggers
  Handle<edm::TriggerResults> triggerHandle;
  iEvent.getByToken(triggerToken_, triggerHandle);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerHandle);

  //std::vector<std::string> filterLabels;
  std::string DETagFilter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
  std::string DEProbeFilter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
  std::string SEFilter("hltEle27WPLooseGsfTrackIsoFilter");

  /*bool trigResult = false;

  for (unsigned int i=0; i<triggerHandle->size(); i++)
  {
    std::string trigName = triggerNames.triggerName(i);
    trigResult = triggerHandle->accept(i);
    cout<<"Name of Trigger = "<<trigName<<"   Trigger Result = "<<trigResult<<"   Trigger Number = "<<i<<endl;
  }*/

  //bool tagPass = false;
  //bool probePass = false;

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackPathNames(triggerNames);
    for (unsigned j = 0; j < obj.filterLabels().size(); ++j){
      //cout<<obj.filterLabels()[j]<<endl;
      if((DETagFilter.compare(obj.filterLabels()[j]))==0){
	//tagPass = true;
	pt_leg1.push_back(obj.pt());
	eta_leg1.push_back(obj.eta());
	phi_leg1.push_back(obj.phi());
      }

      if((DEProbeFilter.compare(obj.filterLabels()[j]))==0){
	//probePass = true;
	pt_leg2.push_back(obj.pt());
	eta_leg2.push_back(obj.eta());
	phi_leg2.push_back(obj.phi());
      }

      if((SEFilter.compare(obj.filterLabels()[j]))==0){
	//probePass = true;
	pt_SE.push_back(obj.pt());
	eta_SE.push_back(obj.eta());
	phi_SE.push_back(obj.phi());
      }

    }
  } 

  hlNames_ = triggerNames.triggerNames();
  int ntriggers = hlNames_.size();
  //cout<<"ntriggers: "<<ntriggers<<endl;
  Int_t hsize = Int_t(triggerHandle->size());
  //cout<<"hsize: "<<hsize<<endl;

  for (int itrigger=0; itrigger<ntriggers; itrigger++)
  {
    std::string hltname(triggerNames.triggerName(itrigger));

    size_t found_singleEle = hltname.find(string_singleEle);
    size_t found_doubleEle = hltname.find(string_doubleEle);

    if(found_singleEle !=string::npos){
      single_electron_triggers_in_run.push_back(hltname);
    }

    if(found_doubleEle !=string::npos){
      double_electron_triggers_in_run.push_back(hltname);
    }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_electron_triggers_in_run.size(); itrigger++){
    idx_singleElectron.push_back(triggerNames.triggerIndex(single_electron_triggers_in_run[itrigger]));
    if(idx_singleElectron.size()>0)
      if(idx_singleElectron[itrigger] < hsize){
	singleElectron = (triggerHandle->accept(idx_singleElectron[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)double_electron_triggers_in_run.size(); itrigger++){
    idx_doubleElectron.push_back(triggerNames.triggerIndex(double_electron_triggers_in_run[itrigger]));
    if(idx_doubleElectron.size()>0)
      if(idx_doubleElectron[itrigger] < hsize){
	doubleElectron = (triggerHandle->accept(idx_doubleElectron[itrigger]));
      }
  }

  /*for (int itrigger=0; itrigger<ntriggers; itrigger++)
  {
    std::string hltname(triggerNames.triggerName(itrigger));
    //cout<<"hltname: "<<hltname<<endl;
    size_t found_doubleEle = hltname.find(string_doubleEle);
    //cout<<"found_doubleEle: "<<found_doubleEle<<endl;
    if(found_doubleEle !=string::npos){
      double_electron_triggers_in_run.push_back(hltname);
    }
  }

  //cout<<"hltname: "<<hltname<<endl;
  //cout<<"found_doubleEle: "<<found_doubleEle<<endl;
  //cout<<"size: "<<(int)double_electron_triggers_in_run.size()<<endl;
  for ( int itrigger = 0 ; itrigger < (int)double_electron_triggers_in_run.size(); itrigger++){

    idx_doubleElectron.push_back(triggerNames.triggerIndex(double_electron_triggers_in_run[itrigger]));
    //cout<<"idx_doubleElectron: "<<idx_doubleElectron[itrigger]<<endl;
    //cout<<"idx_doubleElectron.size(): "<<idx_doubleElectron.size()<<endl;
    if(idx_doubleElectron.size()>0)
      //cout<<"idx_doubleElectron.size(): "<<idx_doubleElectron.size()<<endl;
      if(idx_doubleElectron[itrigger] < hsize){
	//cout<<"idx_doubleElectron[itrigger]: "<<idx_doubleElectron[itrigger]<<endl;
	doubleElectron = (triggerHandle->accept(idx_doubleElectron[itrigger]));
      }
  } 
*/
  // Get Pileup info
  Handle<edm::View<PileupSummaryInfo> > pileupHandle;
  iEvent.getByToken(pileupToken_, pileupHandle);
  if(misMC){
    for( auto & puInfoElement : *pileupHandle){
      if( puInfoElement.getBunchCrossing() == 0 ){
	nPU_    = puInfoElement.getPU_NumInteractions();
	nPUTrue_= puInfoElement.getTrueNumInteractions();
      }
    }
  }

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  

  // We use exactly the same handle for AOD and miniAOD formats since pat::Electron objects can be recast as reco::GsfElectron objects.
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  bool isAOD = true;
  iEvent.getByToken(electronsToken_, electrons);
  if( !electrons.isValid() ){
    isAOD = false;
    iEvent.getByToken(electronsMiniAODToken_,electrons);
  }

  // Get the MC collection

  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  if(misMC){
    nGenElectrons_ = 0;

    for(size_t i = 0; i < genParticles->size(); ++i){
      const GenParticle &e = (*genParticles)[i];
      int id = e.pdgId();
      int st = e.status();
      const Candidate * mother = e.mother();

      if (fabs(id)==11 && st==1) {
	gPost_energy_.push_back(e.energy());
	gPost_px_.push_back(e.px());
	gPost_py_.push_back(e.py());
	gPost_pz_.push_back(e.pz());
	gPost_pt_.push_back(e.pt());
	gPost_eta_.push_back(e.eta());
	gPost_rap_.push_back(e.rapidity());
	gPost_phi_.push_back(e.phi());
      }

      //cout<<"particle no. = "<<i<<"   ID = "<<id<<"   STATUS = "<<st<<endl;
      //  if (id==23 && st==22){cout<<"id: "<<id<<endl;}

      if (fabs(id)==11 && st==23 && mother->pdgId()==23){  
	nGenElectrons_++;

	ZMass_.push_back(mother->mass());
	ZPt_.push_back(mother->pt());
	ZEta_.push_back(mother->eta());
	ZRap_.push_back(mother->rapidity());
	ZPhi_.push_back(mother->phi());

	gPre_energy_.push_back(e.energy());
	gPre_px_.push_back(e.px());
	gPre_py_.push_back(e.py());
	gPre_pz_.push_back(e.pz());
	gPre_pt_.push_back(e.pt());
	gPre_eta_.push_back(e.eta());
	gPre_rap_.push_back(e.rapidity());
	gPre_phi_.push_back(e.phi());

      }
    }
  }

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  if( isAOD )
    iEvent.getByToken(vtxToken_, vertices);
  else
    iEvent.getByToken(vtxMiniAODToken_, vertices);

  if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &pv = vertices->front();
  nPV_    = vertices->size();

  // Find the first vertex in the collection that passes  good quality criteria
  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
      vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = vtx->isFake();
    if( !isAOD )
      isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if ( !isFake
	&&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	&& fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs

  // Seems always zero. Not stored in miniAOD...?
  pvNTracks_ = firstGoodVertex->nTracks();

  // Get the conversions collection
  edm::Handle<reco::ConversionCollection> conversions;
  if(isAOD)
    iEvent.getByToken(conversionsToken_, conversions);
  else
    iEvent.getByToken(conversionsMiniAODToken_, conversions);

  // Get the electron ID data from the event stream.
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
  iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

  nElectrons_ = 0;
  //cout<<"Electrons: "<<electrons->size()<<"      trigger: "<<doubleElectron<<endl;

  // Loop over electrons
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);

    // Kinematics
    if( el->pt() < 10 ) // keep only electrons above 10 GeV
      continue;

    nElectrons_++;
    //cout<<"pt: "<<el->pt()<<endl;
    pt_.push_back( el->pt() );
    eta_.push_back( el->eta() );
    rap_.push_back( el->rapidity() );
    phi_.push_back( el->phi() );
    energy_.push_back( el->energy() );
    mass_.push_back( el->mass() );
    charge_.push_back( el->charge() );

    double R = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y() +el->superCluster()->z()*el->superCluster()->z());
    double Rt = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y());

    enSC_.push_back(el->superCluster()->energy());
    preEnSC_.push_back(el->superCluster()->preshowerEnergy());
    rawEnSC_.push_back(el->superCluster()->rawEnergy());
    etSC_.push_back( (el->superCluster()->energy())*(Rt/R) );
    etaSC_.push_back( el->superCluster()->eta() );
    phiSC_.push_back( el->superCluster()->phi() );

    // ECAL
    full5x5_sigmaIetaIeta_.push_back( el->full5x5_sigmaIetaIeta() );
    E1x5_.push_back(el->e1x5());
    E2x5_.push_back(el->e2x5Max());
    E5x5_.push_back(el->e5x5());
    hOverE_.push_back( el->hcalOverEcal() );
    etaScWidth_.push_back(el->superCluster()->etaWidth());
    phiScWidth_.push_back(el->superCluster()->phiWidth());
    r9_.push_back(el->r9());

    // ECAL + Track
    dEtaIn_.push_back( el->deltaEtaSuperClusterTrackAtVtx() );
    dPhiIn_.push_back( el->deltaPhiSuperClusterTrackAtVtx() );
    // |1/E-1/p| = |1/E - EoverPinner/E| is computed below. The if protects against ecalEnergy == inf or zero
    // (always the case for miniAOD for electrons <5 GeV)
    if( el->ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP_.push_back( 1e30 );
    }else if( !std::isfinite(el->ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP_.push_back( 1e30 );
    }else{
      ooEmooP_.push_back( fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) );
    }

    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();

    // Compute isolation with delta beta correction for PU
    isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
    isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
    isoPhotons_.push_back( pfIso.sumPhotonEt );
    isoChargedFromPU_.push_back( pfIso.sumPUPt );

    float abseta = fabs(el->superCluster()->eta());

    // The effective areas constants file in the local release or default CMSSW, whichever is found
    edm::FileInPath eaConstantsFile("RecoEgamma/ElectronIdentification/data/PHYS14/effAreaElectrons_cone03_pfNeuHadronsAndPhotons.txt");
    EffectiveAreas effectiveAreas(eaConstantsFile.fullPath());
    float eA = effectiveAreas.getEffectiveArea(abseta);

    isoDeltaBeta_.push_back((pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt))/(el->pt()));
    isoRho_.push_back((pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho_ * eA))/(el->pt()));

    // Track - Impact Parameter, Conversion rejection, Converted
    reco::GsfTrackRef theTrack = el->gsfTrack();
    d0_.push_back( (-1) * theTrack->dxy(firstGoodVertex->position() ) );
    dz_.push_back( theTrack->dz( firstGoodVertex->position() ) );
    expectedMissingInnerHits_.push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );
    bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, theBeamSpot->position());
    passConversionVeto_.push_back( (int) passConvVeto );
    brem_.push_back(el->fbrem());

    // Match to generator level truth
    if(misMC) isTrue_.push_back( matchToTruth( el, genParticles) );

    // ID
    bool isPassVeto  = (*veto_id_decisions)[el];
    bool isPassLoose  = (*loose_id_decisions)[el];
    bool isPassMedium = (*medium_id_decisions)[el];
    bool isPassTight  = (*tight_id_decisions)[el];
    passVetoId_.push_back  ( (int)isPassVeto  );
    passLooseId_.push_back ( (int)isPassLoose );
    passMediumId_.push_back( (int)isPassMedium);
    passTightId_.push_back ( (int)isPassTight );

    mvaValue_.push_back( (*mvaValues)[el] );
    mvaCategory_.push_back( (*mvaCategories)[el] );

    eleInBarrel_.push_back(el->isEB());
    eleInEndcap_.push_back(el->isEE());

  }

  // Save this electron's info
  electronTree_->Fill();

  hlNames_.clear();
  single_electron_triggers_in_run.clear();
  double_electron_triggers_in_run.clear();
  idx_singleElectron.clear();
  idx_doubleElectron.clear();

  // Clear vectors
  pt_leg1.clear();
  eta_leg1.clear();
  phi_leg1.clear();
  pt_leg2.clear();
  eta_leg2.clear();
  phi_leg2.clear();

  pt_SE.clear();
  eta_SE.clear();
  phi_SE.clear();

  if(misMC){
    ZMass_.clear();
    ZPt_.clear();
    ZEta_.clear();
    ZRap_.clear();
    ZPhi_.clear();
    gPre_energy_.clear();
    gPre_pt_.clear();
    gPre_px_.clear();
    gPre_py_.clear();
    gPre_pz_.clear();
    gPre_eta_.clear();
    gPre_rap_.clear();
    gPre_phi_.clear();
    gPost_energy_.clear();
    gPost_pt_.clear();
    gPost_px_.clear();
    gPost_py_.clear();
    gPost_pz_.clear();
    gPost_eta_.clear();
    gPost_rap_.clear();
    gPost_phi_.clear();

    isTrue_.clear();
  }

  pt_.clear();
  eta_.clear();
  rap_.clear();
  phi_.clear();
  energy_.clear();
  mass_.clear();
  charge_.clear();

  enSC_.clear();
  preEnSC_.clear();
  rawEnSC_.clear();
  etSC_.clear();
  etaSC_.clear();
  phiSC_.clear();
  dEtaIn_.clear();
  dPhiIn_.clear();
  full5x5_sigmaIetaIeta_.clear();
  E1x5_.clear();
  E2x5_.clear();
  E5x5_.clear();
  hOverE_.clear();
  etaScWidth_.clear();
  phiScWidth_.clear();
  r9_.clear();

  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoChargedFromPU_.clear();
  isoDeltaBeta_.clear();
  isoRho_.clear();
  ooEmooP_.clear();
  d0_.clear();
  dz_.clear();
  expectedMissingInnerHits_.clear();
  passConversionVeto_.clear();
  brem_.clear();

  passVetoId_.clear();
  passLooseId_.clear();
  passMediumId_.clear();
  passTightId_.clear();
  mvaCategory_.clear();
  mvaValue_.clear();

}
// ------------ method called once each job just before starting event loop  ------------
  void 
SimpleElectronNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
SimpleElectronNtupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   SimpleElectronNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   SimpleElectronNtupler::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   SimpleElectronNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   SimpleElectronNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleElectronNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int SimpleElectronNtupler::matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
    const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("SimpleElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }

  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void SimpleElectronNtupler::findFirstNonElectronMother(const reco::Candidate *particle,
    int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("SimpleElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleElectronNtupler);
