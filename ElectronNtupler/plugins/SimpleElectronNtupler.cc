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
#include <regex>

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
#include "TParticle.h"
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
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//#include "/afs/cern.ch/work/r/rchawla/private/Analysis_ee_13TeV/CMSSW_7_4_15/src/EgammaWork/ElectronNtupler/plugins/utils.h"

//#include "HLTrigger/HLTanalyzers/interface/JetUtil.h"
//#include "/afs/cern.ch/work/r/rchawla/private/Analysis_ee_13TeV/CMSSW_7_4_0/src/EgammaWork/ElectronNtupler/plugins/utils.h"
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
    edm::EDGetTokenT<pat::PackedTriggerPrescales>triggerPrescale_; 
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;

    // AOD case data members
    edm::EDGetToken electronsToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

    //edm::EDGetToken VertexToken_;
    
    // MiniAOD case data members
    edm::EDGetToken muonsMiniAODToken_;
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
    //edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
    //edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

    virtual Double_t deltaPhi(Double_t phientry1 , Double_t phientry2 );
    virtual Double_t deltaEta(Double_t etaentry1 , Double_t etaentry2);
    virtual Double_t deltaR(Double_t etaentry1 , Double_t phientry1 , Double_t etaentry2 , Double_t phientry2);
    
    edm::Service<TFileService> fs;
    TTree *electronTree_;

    //edm::LumiReWeighting LumiWeights_;

    // If MC
    bool misMC;
    bool misNLO;
    bool misSIG;

    // Histograms
    TH1F* nevents_;

    // Weights for MC@NLO
    double theWeight;

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
    std::string string_singleMuon;
    std::string string_doubleEle;

    std::string string_emu_17_8;
    std::string string_emu_12_17;
    std::string string_emu_23_8;

    std::string string_photon30;
    std::string string_photon36;
    std::string string_photon50;
    std::string string_photon75;
    std::string string_photon90;
    std::string string_photon120;
    std::string string_photon175;

    // Triggers for Fake-Rate method
    std::regex Photon_30;
    std::regex Photon_36;
    std::regex Photon_50;
    std::regex Photon_75;
    std::regex Photon_90;
    std::regex Photon_120;
    std::regex Photon_175;

    int prescalePhoton_30;
    int prescalePhoton_36;
    int prescalePhoton_50;
    int prescalePhoton_75;
    int prescalePhoton_90;
    int prescalePhoton_120;
    int prescalePhoton_175;

    std::vector<std::string> hlNames_;
    std::vector<std::string> single_electron_triggers_in_run;
    std::vector<std::string> single_muon_triggers_in_run;
    std::vector<std::string> double_electron_triggers_in_run;
    std::vector<std::string> double_emu_17_8_triggers_in_run;
    std::vector<std::string> double_emu_12_17_triggers_in_run;
    std::vector<std::string> double_emu_23_8_triggers_in_run;

    std::vector<std::string> single_photon_30_triggers_in_run;
    std::vector<std::string> single_photon_36_triggers_in_run;
    std::vector<std::string> single_photon_50_triggers_in_run;
    std::vector<std::string> single_photon_75_triggers_in_run;
    std::vector<std::string> single_photon_90_triggers_in_run;
    std::vector<std::string> single_photon_120_triggers_in_run;
    std::vector<std::string> single_photon_175_triggers_in_run;

    bool singleEle;
    bool singleMuon;
    bool doubleElectron;
    bool doubleEMu_17_8;
    bool doubleEMu_12_17;
    bool doubleEMu_23_8;
    
    bool singlePhoton_30;
    bool singlePhoton_36;
    bool singlePhoton_50;
    bool singlePhoton_75;
    bool singlePhoton_90;
    bool singlePhoton_120;
    bool singlePhoton_175;

    std::vector<int> idx_singleEle;
    std::vector<int> idx_singleMuon;
    std::vector<int> idx_doubleElectron;
    std::vector<int> idx_doubleEMu_17_8;
    std::vector<int> idx_doubleEMu_12_17;
    std::vector<int> idx_doubleEMu_23_8;

    std::vector<int> idx_singlePhoton_30;
    std::vector<int> idx_singlePhoton_36;
    std::vector<int> idx_singlePhoton_50;
    std::vector<int> idx_singlePhoton_75;
    std::vector<int> idx_singlePhoton_90;
    std::vector<int> idx_singlePhoton_120;
    std::vector<int> idx_singlePhoton_175;

    // tau variables
    Int_t tauFlag;
    std::vector<Float_t> gen_ptTau_;
    std::vector<Float_t> gen_etaTau_;
    std::vector<Float_t> gen_phiTau_;

    // all electron variables
    Int_t nElectrons_;
    Int_t nGenElectrons_;
    Double_t SumPhotonMom;

    std::vector<Float_t> gen_preFSR_ene_;
    std::vector<Float_t> gen_preFSR_px_;
    std::vector<Float_t> gen_preFSR_py_;
    std::vector<Float_t> gen_preFSR_pz_;
    std::vector<Float_t> gen_preFSR_pt_;
    std::vector<Float_t> gen_preFSR_eta_;
    std::vector<Float_t> gen_preFSR_rap_;
    std::vector<Float_t> gen_preFSR_phi_;

    std::vector<Float_t> gen_postFSR_ene_;
    std::vector<Float_t> gen_postFSR_px_;
    std::vector<Float_t> gen_postFSR_py_;
    std::vector<Float_t> gen_postFSR_pz_;
    std::vector<Float_t> gen_postFSR_pt_;
    std::vector<Float_t> gen_postFSR_eta_;
    std::vector<Float_t> gen_postFSR_rap_;
    std::vector<Float_t> gen_postFSR_phi_;

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

    //std::vector<Float_t> mvaValue_;
    //std::vector<Int_t>   mvaCategory_;

    std::vector<Int_t> eleEcalDrivenSeed_;

    std::vector<double> pt_leg1;
    std::vector<double> eta_leg1;
    std::vector<double> phi_leg1;

    std::vector<double> pt_leg2;
    std::vector<double> eta_leg2;
    std::vector<double> phi_leg2;

    std::vector<double> pt_Ele;
    std::vector<double> eta_Ele;
    std::vector<double> phi_Ele;

    // all muon variables
    Int_t nMuons_;

    std::vector<bool>   isLoose_;
    std::vector<bool>   isTight_;
    std::vector<bool>   isHEEP_;
    std::vector<int>    isGLBmuon_;
    std::vector<int>    isPFmuon_;

    std::vector<double> ptMuon_;
    std::vector<double> etaMuon_;
    std::vector<double> phiMuon_;
    std::vector<double> energyMuon_;
    std::vector<double> chargeMuon_;

    std::vector<double> isoChargedHadronPfR04Muon_;
    std::vector<double> isoNeutralHadronPfR04Muon_;
    std::vector<double> isoGammaPfR04Muon_;
    std::vector<double> isoChargedFromPUMuon_;
    std::vector<double> isoPFMuon_;
    std::vector<double> isoTrkMuon_;

    std::vector<double> metPt_;
    std::vector<double> metPhi_;
    std::vector<double> metSumEt_;

    double DeltaR(const pat::Electron& e, std::vector<pat::TriggerObjectStandAlone> object);
};

SimpleElectronNtupler::SimpleElectronNtupler(const edm::ParameterSet& iConfig):
  //edm::EDAnalyzer(),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap")))
  //mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  //mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{
  nevents_ = fs->make<TH1F>("nevents_","nevents_",2,0,2);
  string_singleEle      = "HLT_Ele23_WPLoose_Gsf_v";
  string_singleMuon     = "HLT_IsoMu20_v";
  string_doubleEle      = "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  string_emu_17_8       = "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v";
  string_emu_12_17      = "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
  string_emu_23_8       = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";

  // triggers for Fake-Rate method
  Photon_30  = "(HLT_Photon30_v)(.*)";
  Photon_36  = "(HLT_Photon36_v)(.*)";
  Photon_50  = "(HLT_Photon50_v)(.*)";
  Photon_75  = "(HLT_Photon75_v)(.*)";
  Photon_90  = "(HLT_Photon90_v)(.*)";
  Photon_120 = "(HLT_Photon120_v)(.*)";
  Photon_175 = "(HLT_Photon175_v)(.*)";

  misMC            = iConfig.getUntrackedParameter<bool>("isMC");
  misSIG           = iConfig.getUntrackedParameter<bool>("isSIG");
  misNLO           = iConfig.getUntrackedParameter<bool>("isNLO");

  // initialize 1-D reweighting
  /*LumiWeights_ = edm::LumiReWeighting("dyJtoLL_M50.root", 
    "MyDataPileupHistogram.root", 
    "MC_nPU", 
    "pileup");*/

  // Prepare tokens for all input collections and objects

  // Trigger
  triggerToken_  = mayConsume<edm::TriggerResults>
    (iConfig.getParameter<edm::InputTag>
     ("trigger"));

  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>
    (iConfig.getParameter<edm::InputTag>
     ("objects"));

  triggerPrescale_ = consumes<pat::PackedTriggerPrescales>
    (iConfig.getParameter<edm::InputTag>
     ("prescale"));

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
  muonsMiniAODToken_        = mayConsume<edm::View<pat::Muon> >
    (iConfig.getParameter<edm::InputTag>
     ("muonsMiniAOD"));

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

  metToken_ = mayConsume<pat::METCollection>
    (iConfig.getParameter<edm::InputTag>
     ("metsMiniAOD"));

  //
  // Set up the ntuple structure
  //
  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

  electronTree_->Branch("RunNo", &RunNo_, "RunNo/D");
  electronTree_->Branch("EvtNo", &EvtNo_, "EvtNo/D");
  electronTree_->Branch("Lumi", &Lumi_, "Lumi/D");
  electronTree_->Branch("Bunch", &Bunch_, "Bunch/D");

  electronTree_->Branch("theWeight", &theWeight, "theWeight/D");

  electronTree_->Branch("pvNTracks"    ,  &pvNTracks_ , "pvNTracks/I");

  electronTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  electronTree_->Branch("nPU"        ,  &nPU_     , "nPU/I");
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue_ , "nPUTrue/I");
  electronTree_->Branch("rho"        ,  &rho_ , "rho/F");

  electronTree_->Branch("singleEle"     ,  &singleEle    );
  electronTree_->Branch("singleMuon"    ,  &singleMuon   );
  electronTree_->Branch("doubleElectron"    ,  &doubleElectron    );
  electronTree_->Branch("doubleEMu_17_8"    ,  &doubleEMu_17_8    );
  electronTree_->Branch("doubleEMu_12_17"    ,  &doubleEMu_12_17  );
  electronTree_->Branch("doubleEMu_23_8"    ,  &doubleEMu_23_8    );

  electronTree_->Branch("prescalePhoton_30", &prescalePhoton_30);
  electronTree_->Branch("prescalePhoton_36", &prescalePhoton_36);
  electronTree_->Branch("prescalePhoton_50", &prescalePhoton_50);
  electronTree_->Branch("prescalePhoton_75", &prescalePhoton_75);
  electronTree_->Branch("prescalePhoton_90", &prescalePhoton_90);
  electronTree_->Branch("prescalePhoton_120", &prescalePhoton_120);
  electronTree_->Branch("prescalePhoton_175", &prescalePhoton_175);

  electronTree_->Branch("singlePhoton_30", &singlePhoton_30);
  electronTree_->Branch("singlePhoton_36", &singlePhoton_36);
  electronTree_->Branch("singlePhoton_50", &singlePhoton_50);
  electronTree_->Branch("singlePhoton_75", &singlePhoton_75);
  electronTree_->Branch("singlePhoton_90", &singlePhoton_90);
  electronTree_->Branch("singlePhoton_120", &singlePhoton_120);
  electronTree_->Branch("singlePhoton_175", &singlePhoton_175);

  electronTree_->Branch("pt_leg1"    ,  &pt_leg1    );
  electronTree_->Branch("eta_leg1"   ,  &eta_leg1   );
  electronTree_->Branch("phi_leg1"   ,  &phi_leg1   );
  electronTree_->Branch("pt_leg2"    ,  &pt_leg2    );
  electronTree_->Branch("eta_leg2"   ,  &eta_leg2   );
  electronTree_->Branch("phi_leg2"   ,  &phi_leg2   );
  electronTree_->Branch("pt_Ele"    ,  &pt_Ele    );
  electronTree_->Branch("eta_Ele"   ,  &eta_Ele   );
  electronTree_->Branch("phi_Ele"   ,  &phi_Ele   );

  electronTree_->Branch("tauFlag", &tauFlag, "tauFlag/I");
  electronTree_->Branch("gen_ptTau", &gen_ptTau_);
  electronTree_->Branch("gen_etaTau", &gen_etaTau_);
  electronTree_->Branch("gen_phiTau", &gen_phiTau_);

  electronTree_->Branch("nEle"    ,  &nElectrons_ , "nEle/I");
  electronTree_->Branch("nGenEle"    ,  &nGenElectrons_ , "nGenEle/I");
  electronTree_->Branch("gen_preFSR_ene"    ,  &gen_preFSR_ene_    );
  electronTree_->Branch("gen_preFSR_px"    ,  &gen_preFSR_px_    );
  electronTree_->Branch("gen_preFSR_py"    ,  &gen_preFSR_py_    );
  electronTree_->Branch("gen_preFSR_pz"    ,  &gen_preFSR_pz_    );
  electronTree_->Branch("gen_preFSR_pt"    ,  &gen_preFSR_pt_    );
  electronTree_->Branch("gen_preFSR_eta"    ,  &gen_preFSR_eta_    );
  electronTree_->Branch("gen_preFSR_rap"    ,  &gen_preFSR_rap_    );
  electronTree_->Branch("gen_preFSR_phi"    ,  &gen_preFSR_phi_    );

  electronTree_->Branch("gen_postFSR_ene"    ,  &gen_postFSR_ene_    );
  electronTree_->Branch("gen_postFSR_px"    ,  &gen_postFSR_px_    );
  electronTree_->Branch("gen_postFSR_py"    ,  &gen_postFSR_py_    );
  electronTree_->Branch("gen_postFSR_pz"    ,  &gen_postFSR_pz_    );
  electronTree_->Branch("gen_postFSR_pt"    ,  &gen_postFSR_pt_    );
  electronTree_->Branch("gen_postFSR_eta"    ,  &gen_postFSR_eta_    );
  electronTree_->Branch("gen_postFSR_rap"    ,  &gen_postFSR_rap_    );
  electronTree_->Branch("gen_postFSR_phi"    ,  &gen_postFSR_phi_    );

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
  electronTree_->Branch("E5x5" ,  &E5x5_ );
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
  //electronTree_->Branch("mvaVal" ,  &mvaValue_ );
  //electronTree_->Branch("mvaCat" ,  &mvaCategory_ );
  electronTree_->Branch("eleEcalDrivenSeed", &eleEcalDrivenSeed_);
  electronTree_->Branch("eleInBarrel", &eleInBarrel_);
  electronTree_->Branch("eleInEndcap", &eleInEndcap_);

  electronTree_->Branch("isGLBmuon", &isGLBmuon_);
  electronTree_->Branch("isPFmuon", &isPFmuon_);

  electronTree_->Branch("nMuons", &nMuons_);
  electronTree_->Branch("isLoose", &isLoose_);
  //electronTree_->Branch("isMedium", &isMedium);
  electronTree_->Branch("isTight", &isTight_);
  electronTree_->Branch("isHEEP", &isHEEP_);

  electronTree_->Branch("ptMuon", &ptMuon_);
  electronTree_->Branch("etaMuon", &etaMuon_);
  electronTree_->Branch("phiMuon", &phiMuon_);
  electronTree_->Branch("energyMuon", &energyMuon_);
  electronTree_->Branch("chargeMuon", &chargeMuon_);

  electronTree_->Branch("isoChargedHadronPfR04Muon", &isoChargedHadronPfR04Muon_);
  electronTree_->Branch("isoNeutralHadronPfR04Muon", &isoNeutralHadronPfR04Muon_);
  electronTree_->Branch("isoGammaPfR04Muon", &isoGammaPfR04Muon_);
  electronTree_->Branch("isoChargedFromPUMuon", &isoChargedFromPUMuon_);
  electronTree_->Branch("isoPFMuon", &isoPFMuon_);
  electronTree_->Branch("isoTrkMuon", &isoTrkMuon_);

  electronTree_->Branch("metPt", &metPt_);
  electronTree_->Branch("metPhi", &metPhi_);
  electronTree_->Branch("metSumEt", &metSumEt_);

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

  nevents_->Fill(1);
  RunNo_ = iEvent.id().run();
  EvtNo_ = iEvent.id().event();
  Lumi_  = iEvent.luminosityBlock();
  Bunch_ = iEvent.bunchCrossing();

  //cout<<"1"<<endl;

  if(misMC){// && misNLO){
    Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator", genEvtInfo);

    //std::vector<double> evtWeights = genEvtInfo->weights();

    if (genEvtInfo.isValid()) {
      theWeight = genEvtInfo->weight();
    }
    //cout<<"Weight: "<<theWeight<<endl;
  }

  // MET
  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByToken(metToken_, metHandle);
  
  const pat::MET &met = metHandle->front();
  //std::cout<<"met pt: "<<met.pt()<< "   "<<"met phi: "<<"   "<<met.phi()<<"   "<<"met sum et: "<<met.sumEt()<<endl;
  metPt_.push_back(met.pt());
  metPhi_.push_back(met.phi());
  metSumEt_.push_back(met.sumEt());

  // Get Triggers
  Handle<edm::TriggerResults> triggerHandle;
  iEvent.getByToken(triggerToken_, triggerHandle);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescale_, triggerPrescales);

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerHandle);

  //std::vector<std::string> filterLabels;
  std::string DETagFilter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
  std::string DEProbeFilter("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
  std::string SEFilter("hltEle23WPLooseGsfTrackIsoFilter");
  //std::string SEDataFilter("hltSingleEle22WPLooseGsfTrackIsoFilter");//hltEle23WP75GsfTrackIsoFilter");


  //bool trigResult = false;

  for (unsigned int i=0; i<triggerHandle->size(); i++)
  {
    std::string trigName = triggerNames.triggerName(i);
    //trigResult = triggerHandle->accept(i);
    //cout<<"Name of Trigger = "<<trigName<<"   Trigger Result = "<<trigResult<<"   Trigger Number = "<<i<<endl;

    if(std::regex_match(triggerNames.triggerName(i),Photon_30)) prescalePhoton_30 = triggerPrescales->getPrescaleForIndex(i);
    if(std::regex_match(triggerNames.triggerName(i),Photon_36)) prescalePhoton_36 = triggerPrescales->getPrescaleForIndex(i);
    if(std::regex_match(triggerNames.triggerName(i),Photon_50)) prescalePhoton_50 = triggerPrescales->getPrescaleForIndex(i);
    if(std::regex_match(triggerNames.triggerName(i),Photon_75)) prescalePhoton_75 = triggerPrescales->getPrescaleForIndex(i);
    if(std::regex_match(triggerNames.triggerName(i),Photon_90)) prescalePhoton_90 = triggerPrescales->getPrescaleForIndex(i);
    if(std::regex_match(triggerNames.triggerName(i),Photon_120)) prescalePhoton_120 = triggerPrescales->getPrescaleForIndex(i);
    if(std::regex_match(triggerNames.triggerName(i),Photon_175)) prescalePhoton_175 = triggerPrescales->getPrescaleForIndex(i);

    //cout<<"Trigger:"<<trigName<<"   "<<"prescale:"<<triggerPrescales->getPrescaleForIndex(i)<<"   "<<(triggerHandle->accept(i) ? "PASS" : "fail (or not run)")<<endl;
  }

  //if(prescalePhoton_30 != 1) cout<<"prescales photon 30: "<<prescalePhoton_30<<endl;
  //cout<<"prescales photon: "<<prescalePhoton_30<<"   "<<prescalePhoton_36<<"   "<<prescalePhoton_50<<"   "<<prescalePhoton_75<<"   "<<prescalePhoton_90<<"   "<<prescalePhoton_120<<"   "<<prescalePhoton_175<<endl;
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackPathNames(triggerNames);
    //cout<<"Trigger: "<<obj.unpackPathNames(triggerNames)<<endl;
    for (unsigned j = 0; j < obj.filterLabels().size(); ++j){
      //cout<<obj.filterLabels()[j]<<endl;
      if((DETagFilter.compare(obj.filterLabels()[j]))==0){
	pt_leg1.push_back(obj.pt());
	eta_leg1.push_back(obj.eta());
	phi_leg1.push_back(obj.phi());
      }

      if((DEProbeFilter.compare(obj.filterLabels()[j]))==0){
	pt_leg2.push_back(obj.pt());
	eta_leg2.push_back(obj.eta());
	phi_leg2.push_back(obj.phi());
      }

      if((SEFilter.compare(obj.filterLabels()[j]))==0){
	pt_Ele.push_back(obj.pt());
	eta_Ele.push_back(obj.eta());
	phi_Ele.push_back(obj.phi());
      }

    } // obj.filterLabels().size()
  } // triggerObjects

  hlNames_ = triggerNames.triggerNames();
  int ntriggers = hlNames_.size();
  Int_t hsize = Int_t(triggerHandle->size());

  for (int itrigger=0; itrigger<ntriggers; itrigger++)
  {
    std::string hltname(triggerNames.triggerName(itrigger));
    size_t found_singleEle   = hltname.find(string_singleEle);
    size_t found_singleMuon = hltname.find(string_singleMuon);
    size_t found_doubleEle = hltname.find(string_doubleEle);
    size_t found_emu_17_8 = hltname.find(string_emu_17_8);
    size_t found_emu_12_17 = hltname.find(string_emu_12_17);
    size_t found_emu_23_8 = hltname.find(string_emu_23_8);

    size_t found_photon_30 = hltname.find(string_photon30);
    size_t found_photon_36 = hltname.find(string_photon36);
    size_t found_photon_50 = hltname.find(string_photon50);
    size_t found_photon_75 = hltname.find(string_photon75);
    size_t found_photon_90 = hltname.find(string_photon90);
    size_t found_photon_120 = hltname.find(string_photon120);
    size_t found_photon_175 = hltname.find(string_photon175);

    if(found_singleEle !=string::npos) single_electron_triggers_in_run.push_back(hltname);
    if(found_singleMuon !=string::npos) single_muon_triggers_in_run.push_back(hltname);
    if(found_doubleEle !=string::npos) double_electron_triggers_in_run.push_back(hltname);
    if(found_emu_17_8 !=string::npos) double_emu_17_8_triggers_in_run.push_back(hltname);
    if(found_emu_12_17 !=string::npos) double_emu_12_17_triggers_in_run.push_back(hltname);
    if(found_emu_23_8 !=string::npos) double_emu_23_8_triggers_in_run.push_back(hltname);

    if(found_photon_30 !=string::npos) single_photon_30_triggers_in_run.push_back(hltname);
    if(found_photon_36 !=string::npos) single_photon_36_triggers_in_run.push_back(hltname);
    if(found_photon_50 !=string::npos) single_photon_50_triggers_in_run.push_back(hltname);
    if(found_photon_75 !=string::npos) single_photon_75_triggers_in_run.push_back(hltname);
    if(found_photon_90 !=string::npos) single_photon_90_triggers_in_run.push_back(hltname);
    if(found_photon_120 !=string::npos) single_photon_120_triggers_in_run.push_back(hltname);
    if(found_photon_175 !=string::npos) single_photon_175_triggers_in_run.push_back(hltname);
  }

  for ( int itrigger = 0 ; itrigger < (int)single_electron_triggers_in_run.size(); itrigger++){
    idx_singleEle.push_back(triggerNames.triggerIndex(single_electron_triggers_in_run[itrigger]));
    if(idx_singleEle.size()>0)
      if(idx_singleEle[itrigger] < hsize){
	singleEle = (triggerHandle->accept(idx_singleEle[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)double_electron_triggers_in_run.size(); itrigger++){
    idx_doubleElectron.push_back(triggerNames.triggerIndex(double_electron_triggers_in_run[itrigger]));
    if(idx_doubleElectron.size()>0)
      if(idx_doubleElectron[itrigger] < hsize){
	doubleElectron = (triggerHandle->accept(idx_doubleElectron[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_muon_triggers_in_run.size(); itrigger++){
    idx_singleMuon.push_back(triggerNames.triggerIndex(single_muon_triggers_in_run[itrigger]));
    if(idx_singleMuon.size()>0)
      if(idx_singleMuon[itrigger] < hsize){
	singleMuon = (triggerHandle->accept(idx_singleMuon[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)double_emu_17_8_triggers_in_run.size(); itrigger++){
    idx_doubleEMu_17_8.push_back(triggerNames.triggerIndex(double_emu_17_8_triggers_in_run[itrigger]));
    if(idx_doubleEMu_17_8.size()>0)
      if(idx_doubleEMu_17_8[itrigger] < hsize){
	doubleEMu_17_8 = (triggerHandle->accept(idx_doubleEMu_17_8[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)double_emu_12_17_triggers_in_run.size(); itrigger++){
    idx_doubleEMu_12_17.push_back(triggerNames.triggerIndex(double_emu_12_17_triggers_in_run[itrigger]));
    if(idx_doubleEMu_12_17.size()>0)
      if(idx_doubleEMu_12_17[itrigger] < hsize){
	doubleEMu_12_17 = (triggerHandle->accept(idx_doubleEMu_12_17[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)double_emu_23_8_triggers_in_run.size(); itrigger++){
    idx_doubleEMu_23_8.push_back(triggerNames.triggerIndex(double_emu_23_8_triggers_in_run[itrigger]));
    if(idx_doubleEMu_23_8.size()>0)
      if(idx_doubleEMu_23_8[itrigger] < hsize){
	doubleEMu_23_8 = (triggerHandle->accept(idx_doubleEMu_23_8[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_photon_30_triggers_in_run.size(); itrigger++){
    idx_singlePhoton_30.push_back(triggerNames.triggerIndex(single_photon_30_triggers_in_run[itrigger]));
    if(idx_singlePhoton_30.size()>0)
      if(idx_singlePhoton_30[itrigger] < hsize){
	singlePhoton_30 = (triggerHandle->accept(idx_singlePhoton_30[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_photon_36_triggers_in_run.size(); itrigger++){
    idx_singlePhoton_36.push_back(triggerNames.triggerIndex(single_photon_36_triggers_in_run[itrigger]));
    if(idx_singlePhoton_36.size()>0)
      if(idx_singlePhoton_36[itrigger] < hsize){
	singlePhoton_36 = (triggerHandle->accept(idx_singlePhoton_36[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_photon_50_triggers_in_run.size(); itrigger++){
    idx_singlePhoton_50.push_back(triggerNames.triggerIndex(single_photon_50_triggers_in_run[itrigger]));
    if(idx_singlePhoton_50.size()>0)
      if(idx_singlePhoton_50[itrigger] < hsize){
	singlePhoton_50 = (triggerHandle->accept(idx_singlePhoton_50[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_photon_75_triggers_in_run.size(); itrigger++){
    idx_singlePhoton_75.push_back(triggerNames.triggerIndex(single_photon_75_triggers_in_run[itrigger]));
    if(idx_singlePhoton_75.size()>0)
      if(idx_singlePhoton_75[itrigger] < hsize){
	singlePhoton_75 = (triggerHandle->accept(idx_singlePhoton_75[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_photon_90_triggers_in_run.size(); itrigger++){
    idx_singlePhoton_90.push_back(triggerNames.triggerIndex(single_photon_90_triggers_in_run[itrigger]));
    if(idx_singlePhoton_90.size()>0)
      if(idx_singlePhoton_90[itrigger] < hsize){
	singlePhoton_90 = (triggerHandle->accept(idx_singlePhoton_90[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_photon_120_triggers_in_run.size(); itrigger++){
    idx_singlePhoton_120.push_back(triggerNames.triggerIndex(single_photon_120_triggers_in_run[itrigger]));
    if(idx_singlePhoton_120.size()>0)
      if(idx_singlePhoton_120[itrigger] < hsize){
	singlePhoton_120 = (triggerHandle->accept(idx_singlePhoton_120[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_photon_175_triggers_in_run.size(); itrigger++){
    idx_singlePhoton_175.push_back(triggerNames.triggerIndex(single_photon_175_triggers_in_run[itrigger]));
    if(idx_singlePhoton_175.size()>0)
      if(idx_singlePhoton_175[itrigger] < hsize){
	singlePhoton_175 = (triggerHandle->accept(idx_singlePhoton_175[itrigger]));
      }
  }

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

  //double MyWeight = LumiWeights_.weight(nPUTrue_);
  //cout<<"MyWeight: "<<MyWeight<<endl;

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

  if(misMC && misSIG){
    nGenElectrons_ = 0;
    tauFlag = 0;
    SumPhotonMom = 0.;

    for(size_t i = 0; i < genParticles->size(); ++i){
      const GenParticle &genlep = (*genParticles)[i];

      int id = genlep.pdgId();

      TLorentzVector fourmom;
      fourmom.SetPxPyPzE(genlep.px(), genlep.py(), genlep.pz(), genlep.energy());
      //TVector3 p = fourmom.Vect();

      /*int n = genlep.numberOfDaughters();
	if(fabs(id)==11){
	for(int j = 0; j < n; ++j) {
	const Candidate * d = genlep.daughter(j);
	int dauId = d->pdgId();
      //cout<<"ID: "<<id<<"   "<<"Daughter ID: "<<dauId<<endl;

      if(fabs(dauId) == 22){
      Double_t dR = deltaR(d->eta(), d->phi(), genlep.eta(), genlep.phi());
      cout<<"dR: "<<dR<<"   "<<"ID: "<<id<<"   "<<"Daughter ID: "<<dauId<<endl;

      // Sum of all photon's momentum near the post-FSR electron
      if(dR<0.1)
      {
      SumPhotonMom = SumPhotonMom + fourmom.Momentum;
      }
      }
      }
      }*/

      //const Candidate * Mother = genlep.mother(i);
      //int momId = Mother->pdgId();
      //cout<<"Mother  ID:"<<momId<<endl;

      if(fabs(id)==11 && genlep.fromHardProcessFinalState()==1){ // post FSR
	gen_postFSR_ene_.push_back(genlep.energy());
	gen_postFSR_px_.push_back(genlep.px());
	gen_postFSR_py_.push_back(genlep.py());
	gen_postFSR_pz_.push_back(genlep.pz());
	gen_postFSR_pt_.push_back(genlep.pt());
	gen_postFSR_eta_.push_back(genlep.eta());
	gen_postFSR_rap_.push_back(genlep.rapidity());
	gen_postFSR_phi_.push_back(genlep.phi());

      }

      // Only for the photons whose mother is electron or anti-electron
      /*if( fabs(id) == 22 && fabs(Mother->pdgId()) == 11)
	{

	Double_t dR = deltaR(gen_postFSR_eta_[i], gen_postFSR_phi_[i], genlep.eta(), genlep.phi());
      //cout<<"dR: "<<dR<<endl;

      // Sum of all photon's momentum near the post-FSR electron
      if(dR<0.1)
      {
      //SumPhotonMom = SumPhotonMom + p;
      }

      }*/

      if(fabs(id)==11 && genlep.fromHardProcessBeforeFSR()==1){
	nGenElectrons_++;

	gen_preFSR_ene_.push_back(genlep.energy());
	gen_preFSR_px_.push_back(genlep.px());
	gen_preFSR_py_.push_back(genlep.py());
	gen_preFSR_pz_.push_back(genlep.pz());
	gen_preFSR_pt_.push_back(genlep.pt());
	gen_preFSR_eta_.push_back(genlep.eta());
	gen_preFSR_rap_.push_back(genlep.rapidity());
	gen_preFSR_phi_.push_back(genlep.phi());

      }

      // Separating the taus coming from Z decay
      if(abs(id)==15 && genlep.fromHardProcessDecayed()==1){
	gen_ptTau_.push_back(genlep.pt());
	gen_etaTau_.push_back(genlep.eta());
	gen_phiTau_.push_back(genlep.phi());
      }

    }

    //genlep_preFSR->Momentum = genlep_Mom_postFSR + SumPhotonMom;

    // Select the events containing 2 taus from hard-process
    if(gen_ptTau_.size()==2){
      tauFlag = 1;
    }

  }

  //cout<<"gen post FSR: "<<gen_postFSR_pt_.size()<<endl;

  //cout<<"2"<<endl;

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  if( isAOD )
    iEvent.getByToken(vtxToken_, vertices);
  else
    iEvent.getByToken(vtxMiniAODToken_, vertices);

  if (vertices->empty()) return; // skip the event if no PV found
  nPV_ = vertices->size();
  int firstGoodVertexIdx = 0;

  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    firstGoodVertex = vtx;
    break;
  }

  if (firstGoodVertex==vertices->end())
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
  //edm::Handle<edm::ValueMap<float> > mvaValues;
  //edm::Handle<edm::ValueMap<int> > mvaCategories;
  //iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  //iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

  nElectrons_ = 0;
  //cout<<"Event: "<<iEvent.id().event()<<"   "<<"Electrons: "<<electrons->size()<<endl;//"      trigger: "<<doubleElectron<<endl;

  // Loop over electrons
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);

    // Kinematics

    if(el->pt() > 10. && el->eta() < 2.7){

      nElectrons_++;
      //cout<<"ele pt: "<<el->pt()<<"   "<<"ele eta: "<<el->eta()<<endl;
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
	//printf("Electron energy is zero!\n");
	ooEmooP_.push_back( 1e30 );
      }else if( !std::isfinite(el->ecalEnergy())){
	//printf("Electron energy is not finite!\n");
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

      //mvaValue_.push_back( (*mvaValues)[el] );
      //mvaCategory_.push_back( (*mvaCategories)[el] );

      eleInBarrel_.push_back(el->isEB());
      eleInEndcap_.push_back(el->isEE());

      // ECAL driven
      eleEcalDrivenSeed_.push_back(el->ecalDrivenSeed());

    }
  }

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonsMiniAODToken_,muons);

  //cout<<"Muons: "<<muons->size()<<endl;
  nMuons_ = 0;

  for(unsigned i=0; i < muons->size();++i ) {
    const auto mu = muons->ptrAt(i);

    if(mu->pt() > 10. && mu->eta() < 2.7){
      nMuons_++;

      isGLBmuon_.push_back(mu->isGlobalMuon());
      isPFmuon_.push_back(mu->isPFMuon());

      // Kinematics
      ptMuon_.push_back(mu->pt());
      etaMuon_.push_back(mu->eta());
      phiMuon_.push_back(mu->phi());
      energyMuon_.push_back(mu->energy());
      chargeMuon_.push_back(mu->charge());

      // isLoose
      bool isLooseId = muons->at(i).isLooseMuon();
      isLoose_.push_back(isLooseId);

      // isTight
      bool isTightId = muons->at(i).isTightMuon(vertices->at(0));
      isTight_.push_back(isTightId);

      // is HEEP
      bool isHEEPId = muons->at(i).isHighPtMuon(vertices->at(0));
      isHEEP_.push_back(isHEEPId);

      // pf isolation 
      isoChargedHadronPfR04Muon_.push_back(mu->pfIsolationR04().sumChargedHadronPt);
      isoNeutralHadronPfR04Muon_.push_back(mu->pfIsolationR04().sumNeutralHadronEt);
      isoGammaPfR04Muon_.push_back(mu->pfIsolationR04().sumPhotonEt);
      isoChargedFromPUMuon_.push_back(mu->pfIsolationR04().sumPUPt);

      isoPFMuon_.push_back((mu->pfIsolationR04().sumChargedHadronPt + max<float>(0.0, mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5 * (mu->pfIsolationR04().sumPUPt)))/(mu->pt()));

      // tracker isolation
      isoTrkMuon_.push_back((mu->isolationR03().sumPt)/mu->pt());

    }
  }

  //cout<<"4"<<endl;

  // Save this electron's info
  electronTree_->Fill();

  hlNames_.clear();
  single_electron_triggers_in_run.clear();
  single_muon_triggers_in_run.clear();
  double_electron_triggers_in_run.clear();
  double_emu_17_8_triggers_in_run.clear();
  double_emu_12_17_triggers_in_run.clear();
  double_emu_23_8_triggers_in_run.clear();

  single_photon_30_triggers_in_run.clear();
  single_photon_36_triggers_in_run.clear();
  single_photon_50_triggers_in_run.clear();
  single_photon_75_triggers_in_run.clear();
  single_photon_90_triggers_in_run.clear();
  single_photon_120_triggers_in_run.clear();
  single_photon_175_triggers_in_run.clear();

  idx_singleEle.clear();
  idx_singleMuon.clear();
  idx_doubleElectron.clear();
  idx_doubleEMu_17_8.clear();
  idx_doubleEMu_12_17.clear();
  idx_doubleEMu_23_8.clear();

  idx_singlePhoton_30.clear();
  idx_singlePhoton_36.clear();
  idx_singlePhoton_50.clear();
  idx_singlePhoton_75.clear();
  idx_singlePhoton_90.clear();
  idx_singlePhoton_120.clear();
  idx_singlePhoton_175.clear();

  // Clear vectors
  pt_leg1.clear();
  eta_leg1.clear();
  phi_leg1.clear();
  pt_leg2.clear();
  eta_leg2.clear();
  phi_leg2.clear();

  pt_Ele.clear();
  eta_Ele.clear();
  phi_Ele.clear();

  if(misMC){
    gen_ptTau_.clear();
    gen_etaTau_.clear();
    gen_phiTau_.clear();

    gen_preFSR_ene_.clear();
    gen_preFSR_pt_.clear();
    gen_preFSR_px_.clear();
    gen_preFSR_py_.clear();
    gen_preFSR_pz_.clear();
    gen_preFSR_eta_.clear();
    gen_preFSR_rap_.clear();
    gen_preFSR_phi_.clear();

    gen_postFSR_ene_.clear();
    gen_postFSR_pt_.clear();
    gen_postFSR_px_.clear();
    gen_postFSR_py_.clear();
    gen_postFSR_pz_.clear();
    gen_postFSR_eta_.clear();
    gen_postFSR_rap_.clear();
    gen_postFSR_phi_.clear();

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
  //mvaCategory_.clear();
  //mvaValue_.clear();
  eleEcalDrivenSeed_.clear();

  eleInBarrel_.clear();
  eleInEndcap_.clear();

  isLoose_.clear();
  isTight_.clear();
  isHEEP_.clear();
  isGLBmuon_.clear();
  isPFmuon_.clear();

  ptMuon_.clear();
  etaMuon_.clear();
  phiMuon_.clear();
  energyMuon_.clear();
  chargeMuon_.clear();

  isoChargedHadronPfR04Muon_.clear();
  isoNeutralHadronPfR04Muon_.clear();
  isoGammaPfR04Muon_.clear();
  isoChargedFromPUMuon_.clear();
  isoPFMuon_.clear();
  isoTrkMuon_.clear();

  metPt_.clear();
  metPhi_.clear();
  metSumEt_.clear();

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

  Double_t SimpleElectronNtupler::deltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t pi = 3.1415927;
    Double_t dphi = fabs(phi1 - phi2);
    if(dphi >= pi) dphi = 2. * pi - dphi;
    return dphi;
  }

  Double_t SimpleElectronNtupler::deltaEta(Double_t eta1, Double_t eta2)
  {
    Double_t deta = fabs(eta1 - eta2);
    return deta;
  }

  Double_t SimpleElectronNtupler::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
  {
    Double_t deta = deltaEta(eta1, eta2);
    Double_t dphi = deltaPhi(phi1, phi2);
    Double_t dr = sqrt(deta*deta + dphi*dphi);
    return dr;
  }

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
