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

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

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
#include "TClonesArray.h"

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

#include "/afs/cern.ch/work/r/rchawla/public/utils.h"

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

    // MiniAOD case data members
    edm::EDGetToken muonsMiniAODToken_;
    edm::EDGetToken electronsMiniAODToken_;
    edm::EDGetToken photonsMiniAODToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;

    // ID decisions objects
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

    edm::Service<TFileService> fs;
    TTree *electronTree_;

    //edm::LumiReWeighting LumiWeights_;

    // If MC
    bool misMC;
    bool misSIG;

    // Weights for MC@NLO
    double theWeight;

    // General
    Double_t RunNo, EvtNo, Lumi, Bunch;

    // Vars for PVs
    Int_t pvNTracks;

    // Vars for pile-up
    Int_t nPUTrue;    // true pile-up
    Int_t nPU;        // generated pile-up
    Int_t nPV;        // number of reconsrtucted primary vertices
    Float_t rho;      // the rho variable
    //Double_t PUWeight_;

    // Trigger
    std::string string_singleEle23;
    std::string string_singleEle27;
    std::string string_singleMuon;
    std::string string_doubleEle;

    std::string string_emu_17_8;
    std::string string_emu_12_17;
    std::string string_emu_23_8;

    std::vector<std::string> hlNames_;
    std::vector<std::string> single_electron_23_triggers_in_run;
    std::vector<std::string> single_electron_27_triggers_in_run;
    std::vector<std::string> single_muon_triggers_in_run;
    std::vector<std::string> double_electron_triggers_in_run;
    std::vector<std::string> double_emu_17_8_triggers_in_run;
    std::vector<std::string> double_emu_12_17_triggers_in_run;
    std::vector<std::string> double_emu_23_8_triggers_in_run;

    bool single_Ele23;
    bool single_Ele27;
    bool singleMuon;
    bool doubleElectron;
    bool doubleEMu_17_8;
    bool doubleEMu_12_17;
    bool doubleEMu_23_8;

    std::vector<int> idx_singleEle23;
    std::vector<int> idx_singleEle27;
    std::vector<int> idx_singleMuon;
    std::vector<int> idx_doubleElectron;
    std::vector<int> idx_doubleEMu_17_8;
    std::vector<int> idx_doubleEMu_12_17;
    std::vector<int> idx_doubleEMu_23_8;

    // tau variables
    Int_t tauFlag;
    Int_t nGenTaus;

    // all electron variables
    Int_t nGenElectrons;

    std::vector<Int_t>   genlep_Id_;
    std::vector<Int_t>   genlepMother_Id_;
    std::vector<Bool_t>  fromHProcessFinalState_;
    std::vector<Bool_t>  fromHProcessDecayed_;

    std::vector<Float_t> genlep_Px_;
    std::vector<Float_t> genlep_Py_;
    std::vector<Float_t> genlep_Pz_;
    std::vector<Float_t> genlep_Pt_;
    std::vector<Float_t> genlep_Eta_;
    std::vector<Float_t> genlep_Rap_;
    std::vector<Float_t> genlep_Phi_;
    std::vector<Float_t> genlep_En_;

    std::vector<Float_t> genPostFSR_Pt_;
    std::vector<Float_t> genPostFSR_Eta_;
    std::vector<Float_t> genPostFSR_Phi_;
    std::vector<Float_t> genPostFSR_En_;
    
    std::vector<Float_t> genPhoton_Pt_;
    std::vector<Float_t> genPhoton_Eta_;
    std::vector<Float_t> genPhoton_Phi_;
    std::vector<Float_t> genPhoton_En_;

    Int_t nElectrons;
    
    std::vector<Float_t> ptElec_;
    std::vector<Float_t> etaElec_;
    std::vector<Float_t> rapElec_;
    std::vector<Float_t> phiElec_;
    std::vector<Float_t> energyElec_;
    std::vector<Float_t> massElec_;
    std::vector<Float_t> chargeElec_;

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
    std::vector<Int_t> eleEcalDrivenSeed_;

    // Trigger objects
    std::vector<double> pt_Ele23;
    std::vector<double> eta_Ele23;
    std::vector<double> phi_Ele23;

    // Triggers for Fake-Rate method
    std::regex Photon_22;
    std::regex Photon_30;
    std::regex Photon_36;
    std::regex Photon_50;
    std::regex Photon_75;
    std::regex Photon_90;
    std::regex Photon_120;
    std::regex Photon_175;

    int singlePhoton;
    int prescalePhoton;

    bool photon22;
    bool photon30;
    bool photon36;
    bool photon50;
    bool photon75;
    bool photon90;
    bool photon120;
    bool photon175;

    std::vector<double> etPhoton_;

    // all muon variables
    Int_t nMuons;

    std::vector<bool>   isLooseMuon_;
    std::vector<bool>   isTightMuon_;
    std::vector<bool>   isHighPtMuon_;

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

    // MET variables
    std::vector<double> metPt_;
    std::vector<double> metPhi_;
    std::vector<double> metSumEt_;

    // photon variables
    Int_t nPhotons;

    std::vector<Float_t> ptPhoton_;
    std::vector<Float_t> etaPhoton_;
    std::vector<Float_t> phiPhoton_;

    double DeltaR(const pat::Electron& e, std::vector<pat::TriggerObjectStandAlone> object);

};

SimpleElectronNtupler::SimpleElectronNtupler(const edm::ParameterSet& iConfig):
  
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap")))

{

  string_singleEle23    = "HLT_Ele23_WPLoose_Gsf_v";
  string_singleEle27    = "HLT_Ele27_WP85_Gsf_v";
  string_singleMuon     = "HLT_IsoMu20_v";
  string_doubleEle      = "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  string_emu_17_8       = "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v";
  string_emu_12_17      = "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
  string_emu_23_8       = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";

  // triggers for Fake-Rate method
  Photon_22  = "(HLT_Photon22_v)(.*)";
  Photon_30  = "(HLT_Photon30_v)(.*)";
  Photon_36  = "(HLT_Photon36_v)(.*)";
  Photon_50  = "(HLT_Photon50_v)(.*)";
  Photon_75  = "(HLT_Photon75_v)(.*)";
  Photon_90  = "(HLT_Photon90_v)(.*)";
  Photon_120 = "(HLT_Photon120_v)(.*)";
  Photon_175 = "(HLT_Photon175_v)(.*)";

  misMC            = iConfig.getUntrackedParameter<bool>("isMC");
  misSIG           = iConfig.getUntrackedParameter<bool>("isSIG");

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

  photonsMiniAODToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));

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

  electronTree_->Branch("RunNo", &RunNo, "RunNo/D");
  electronTree_->Branch("EvtNo", &EvtNo, "EvtNo/D");
  electronTree_->Branch("Lumi", &Lumi, "Lumi/D");
  electronTree_->Branch("Bunch", &Bunch, "Bunch/D");

  electronTree_->Branch("theWeight", &theWeight, "theWeight/D");

  electronTree_->Branch("pvNTracks"    ,  &pvNTracks , "pvNTracks/I");

  electronTree_->Branch("nPV"        ,  &nPV     , "nPV/I");
  electronTree_->Branch("nPU"        ,  &nPU     , "nPU/I");
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue , "nPUTrue/I");
  electronTree_->Branch("rho"        ,  &rho , "rho/F");
  //electronTree_->Branch("PUWeight", &PUWeight_, "PUWeight/D");

  electronTree_->Branch("single_Ele23"     ,  &single_Ele23    );
  electronTree_->Branch("single_Ele27"     ,  &single_Ele27    );
  electronTree_->Branch("singleMuon"    ,  &singleMuon   );
  electronTree_->Branch("doubleElectron"    ,  &doubleElectron    );
  electronTree_->Branch("doubleEMu_17_8"    ,  &doubleEMu_17_8    );
  electronTree_->Branch("doubleEMu_12_17"    ,  &doubleEMu_12_17  );
  electronTree_->Branch("doubleEMu_23_8"    ,  &doubleEMu_23_8    );

  electronTree_->Branch("etPhoton", &etPhoton_);
  electronTree_->Branch("singlePhoton", &singlePhoton);
  electronTree_->Branch("prescalePhoton", &prescalePhoton);

  electronTree_->Branch("pt_Ele23"    ,  &pt_Ele23    );
  electronTree_->Branch("eta_Ele23"   ,  &eta_Ele23   );
  electronTree_->Branch("phi_Ele23"   ,  &phi_Ele23   );

  electronTree_->Branch("nEle"    ,  &nElectrons , "nEle/I");
  electronTree_->Branch("nGenEle"    ,  &nGenElectrons , "nGenEle/I");
  electronTree_->Branch("nGenTau"    ,  &nGenTaus , "nGenTau/I");
  electronTree_->Branch("tauFlag", &tauFlag, "tauFlag/I");

  electronTree_->Branch("genlep_Id", &genlep_Id_);
  electronTree_->Branch("genlepMother_Id", &genlepMother_Id_);
  electronTree_->Branch("fromHProcessFinalState", &fromHProcessFinalState_);
  electronTree_->Branch("fromHProcessDecayed", &fromHProcessDecayed_);
  
  electronTree_->Branch("genlep_Px"    ,  &genlep_Px_    );
  electronTree_->Branch("genlep_Py"    ,  &genlep_Py_    );
  electronTree_->Branch("genlep_Pz"    ,  &genlep_Pz_    );
  electronTree_->Branch("genlep_Pt"    ,  &genlep_Pt_    );
  electronTree_->Branch("genlep_Eta"   ,  &genlep_Eta_    );
  electronTree_->Branch("genlep_Rap"   ,  &genlep_Rap_    );
  electronTree_->Branch("genlep_Phi"   ,  &genlep_Phi_    );
  electronTree_->Branch("genlep_En"    ,  &genlep_En_    );

  electronTree_->Branch("genPostFSR_Pt"    ,  &genPostFSR_Pt_    );
  electronTree_->Branch("genPostFSR_Eta"   ,  &genPostFSR_Eta_    );
  electronTree_->Branch("genPostFSR_Phi"   ,  &genPostFSR_Phi_    );
  electronTree_->Branch("genPostFSR_En"    ,  &genPostFSR_En_    );
  
  electronTree_->Branch("genPhoton_Pt"    ,  &genPhoton_Pt_    );
  electronTree_->Branch("genPhoton_Eta"   ,  &genPhoton_Eta_    );
  electronTree_->Branch("genPhoton_Phi"   ,  &genPhoton_Phi_    );
  electronTree_->Branch("genPhoton_En"    ,  &genPhoton_En_    );

  electronTree_->Branch("ptElec"    ,  &ptElec_    );
  electronTree_->Branch("etaElec"   ,  &etaElec_    );
  electronTree_->Branch("rapElec"   ,  &rapElec_    );
  electronTree_->Branch("phiElec"   ,  &phiElec_    );
  electronTree_->Branch("energyElec",  &energyElec_    );
  electronTree_->Branch("massElec"  ,  &massElec_    );
  electronTree_->Branch("chargeElec",  &chargeElec_    );
  electronTree_->Branch("enSC"    ,  &enSC_ );
  electronTree_->Branch("preEnSC" ,  &preEnSC_ );
  electronTree_->Branch("rawEnSC" ,  &rawEnSC_ );
  electronTree_->Branch("etSC"  ,  &etSC_ );
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
  electronTree_->Branch("eleEcalDrivenSeed", &eleEcalDrivenSeed_);
  electronTree_->Branch("eleInBarrel", &eleInBarrel_);
  electronTree_->Branch("eleInEndcap", &eleInEndcap_);

  electronTree_->Branch("nMuons", &nMuons);
  electronTree_->Branch("isLooseMuon", &isLooseMuon_);
  electronTree_->Branch("isTightMuon", &isTightMuon_);
  electronTree_->Branch("isHighPtMuon", &isHighPtMuon_);

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

  electronTree_->Branch("nPhotons",  &nPhotons , "nPhotons/I");
  electronTree_->Branch("ptPhoton"  ,  &ptPhoton_    );
  electronTree_->Branch("etaPhoton" ,  &etaPhoton_ );
  electronTree_->Branch("phiPhoton" ,  &phiPhoton_ );

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

  RunNo = iEvent.id().run();
  EvtNo = iEvent.id().event();
  Lumi  = iEvent.luminosityBlock();
  Bunch = iEvent.bunchCrossing();

  //cout<<"1"<<endl;

  if(misMC){

    Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator", genEvtInfo);

    if (genEvtInfo.isValid()) {
      theWeight = genEvtInfo->weight();
    }
  }

  // Get Pileup info

  if(misMC){
    Handle<edm::View<PileupSummaryInfo> > pileupHandle;
    iEvent.getByToken(pileupToken_, pileupHandle);
    for( auto & puInfoElement : *pileupHandle){
      if( puInfoElement.getBunchCrossing() == 0 ){
	nPU     = puInfoElement.getPU_NumInteractions();
	nPUTrue = puInfoElement.getTrueNumInteractions();
      }
    }
  }

  //cout<<"2"<<endl;
  
  // Get Triggers
  Handle<edm::TriggerResults> triggerHandle;
  iEvent.getByToken(triggerToken_, triggerHandle);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescale_, triggerPrescales);

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerHandle);

  // Filter Labels
  std::string SEFilter("hltEle23WPLooseGsfTrackIsoFilter");

  std::string photon22Filter("hltEG22HEFilter");
  std::string photon30Filter("hltEG30HEFilter");
  std::string photon36Filter("hltEG36HEFilter");
  std::string photon50Filter("hltEG50HEFilter");
  std::string photon75Filter("hltEG75HEFilter");
  std::string photon90Filter("hltEG90HEFilter");
  std::string photon120Filter("hltEG120HEFilter");
  std::string photon175Filter("hltEG175HEFilter");

  /*bool trigResult = false;
    for (unsigned int i=0; i<triggerHandle->size(); i++)
    {
    std::string trigName = triggerNames.triggerName(i);
  //trigResult = triggerHandle->accept(i);
  //cout<<"Name of Trigger = "<<trigName<<"   Trigger Result = "<<trigResult<<"   Trigger Number = "<<i<<endl;
  }*/

  singlePhoton = 0;
  prescalePhoton = 0;

  photon22 = false; photon30 = false; photon36 = false; photon50 = false; photon75 = false; photon90 = false; photon120 = false; photon175 = false;

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackPathNames(triggerNames);
    //cout<<"Trigger: "<<obj.unpackPathNames(triggerNames)<<endl;
    for (unsigned j = 0; j < obj.filterLabels().size(); ++j){
      
      if((SEFilter.compare(obj.filterLabels()[j]))==0){
	pt_Ele23.push_back(obj.pt());
	eta_Ele23.push_back(obj.eta());
	phi_Ele23.push_back(obj.phi());
      }

      if((photon22Filter.compare(obj.filterLabels()[j]))==0 || (photon30Filter.compare(obj.filterLabels()[j]))==0 || (photon36Filter.compare(obj.filterLabels()[j]))==0 || (photon50Filter.compare(obj.filterLabels()[j]))==0 || (photon75Filter.compare(obj.filterLabels()[j]))==0 || (photon90Filter.compare(obj.filterLabels()[j]))==0 || (photon120Filter.compare(obj.filterLabels()[j]))==0 || (photon175Filter.compare(obj.filterLabels()[j]))==0){

	if(obj.et() >= 22. && obj.et() < 30.){
	  if((photon22Filter.compare(obj.filterLabels()[j]))==0){
	    photon22 = true;
	    etPhoton_.push_back(obj.et());
	  }
	}

	if(obj.et() >= 30. && obj.et() < 36.){
	  if((photon30Filter.compare(obj.filterLabels()[j]))==0){
	    photon30 = true;
	    etPhoton_.push_back(obj.et()); 
	  }
	}

	if(obj.et() >=36. && obj.et() < 50.){
	  if((photon36Filter.compare(obj.filterLabels()[j]))==0){
	    photon36 = true;
	    etPhoton_.push_back(obj.et());
	  }
	}

	if(obj.et() >= 50. && obj.et() < 75.){
	  if((photon50Filter.compare(obj.filterLabels()[j]))==0){
	    photon50 = true;
	    etPhoton_.push_back(obj.et());
	  }
	}

	if(obj.et() >= 75. && obj.et() < 90.){
	  if((photon75Filter.compare(obj.filterLabels()[j]))==0){
	    photon75 = true;
	    etPhoton_.push_back(obj.et());
	  }
	}

	if(obj.et() >= 90. && obj.et() < 120.){
	  if((photon90Filter.compare(obj.filterLabels()[j]))==0){
	    photon90 = true;
	    etPhoton_.push_back(obj.et());
	  }
	}

	if(obj.et() >= 120. && obj.et() < 175.){
	  if((photon120Filter.compare(obj.filterLabels()[j]))==0){
	    photon120 = true;
	    etPhoton_.push_back(obj.et());
	  }
	}

	if(obj.et() >= 175.){
	  if((photon175Filter.compare(obj.filterLabels()[j]))==0){
	    photon175 = true;
	    etPhoton_.push_back(obj.et());
	  }
	}

      } // OR of triggers
    } // obj.filterLabels().size()
  } // triggerObjects

  for (unsigned int i=0; i<triggerHandle->size(); i++)
  {
    if(photon22){
      if(std::regex_match(triggerNames.triggerName(i),Photon_22)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon30){
      if(std::regex_match(triggerNames.triggerName(i),Photon_30)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon36){
      if(std::regex_match(triggerNames.triggerName(i),Photon_36)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon50){
      if(std::regex_match(triggerNames.triggerName(i),Photon_50)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon75){
      if(std::regex_match(triggerNames.triggerName(i),Photon_75)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon90){
      if(std::regex_match(triggerNames.triggerName(i),Photon_90)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon120){
      if(std::regex_match(triggerNames.triggerName(i),Photon_120)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon175){
      if(std::regex_match(triggerNames.triggerName(i),Photon_175)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

  }

  //cout<<"singlePhoton: "<<singlePhoton<<"   "<<"prescalePhoton: "<<prescalePhoton<<endl;

  //cout<<"3"<<endl;
  
  hlNames_ = triggerNames.triggerNames();
  int ntriggers = hlNames_.size();
  Int_t hsize = Int_t(triggerHandle->size());

  for (int itrigger=0; itrigger<ntriggers; itrigger++)
  {
    std::string hltname(triggerNames.triggerName(itrigger));
    size_t found_singleEle23   = hltname.find(string_singleEle23);
    size_t found_singleEle27   = hltname.find(string_singleEle27);
    size_t found_singleMuon = hltname.find(string_singleMuon);
    size_t found_doubleEle = hltname.find(string_doubleEle);
    size_t found_emu_17_8 = hltname.find(string_emu_17_8);
    size_t found_emu_12_17 = hltname.find(string_emu_12_17);
    size_t found_emu_23_8 = hltname.find(string_emu_23_8);

    if(found_singleEle23 !=string::npos) single_electron_23_triggers_in_run.push_back(hltname);
    if(found_singleEle27 !=string::npos) single_electron_27_triggers_in_run.push_back(hltname);
    if(found_singleMuon !=string::npos) single_muon_triggers_in_run.push_back(hltname);
    if(found_doubleEle !=string::npos) double_electron_triggers_in_run.push_back(hltname);
    if(found_emu_17_8 !=string::npos) double_emu_17_8_triggers_in_run.push_back(hltname);
    if(found_emu_12_17 !=string::npos) double_emu_12_17_triggers_in_run.push_back(hltname);
    if(found_emu_23_8 !=string::npos) double_emu_23_8_triggers_in_run.push_back(hltname);

  }

  for ( int itrigger = 0 ; itrigger < (int)single_electron_23_triggers_in_run.size(); itrigger++){
    idx_singleEle23.push_back(triggerNames.triggerIndex(single_electron_23_triggers_in_run[itrigger]));
    if(idx_singleEle23.size()>0)
      if(idx_singleEle23[itrigger] < hsize){
	single_Ele23 = (triggerHandle->accept(idx_singleEle23[itrigger]));
      }
  }

  for ( int itrigger = 0 ; itrigger < (int)single_electron_27_triggers_in_run.size(); itrigger++){
    idx_singleEle27.push_back(triggerNames.triggerIndex(single_electron_27_triggers_in_run[itrigger]));
    if(idx_singleEle27.size()>0)
      if(idx_singleEle27[itrigger] < hsize){
	single_Ele27 = (triggerHandle->accept(idx_singleEle27[itrigger]));
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

  //cout<<"4"<<endl;
  
  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho = *rhoH;

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

  //cout<<"##################################"<<endl;

  // Get the MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  if(misMC && misSIG){

    nGenElectrons = 0;
    nGenTaus = 0;
    tauFlag  = 0;

    for(size_t i = 0; i < genParticles->size(); ++i){
      const GenParticle &genlep = (*genParticles)[i];

      genlep_Id_.push_back(genlep.pdgId());
      if(abs(genlep.pdgId())==22) genlepMother_Id_.push_back(genlep.mother(0)->pdgId());
      
      fromHProcessFinalState_.push_back(genlep.fromHardProcessFinalState());
      fromHProcessDecayed_.push_back(genlep.fromHardProcessDecayed());

      nGenElectrons++;
      
      genlep_Px_.push_back(genlep.px());
      genlep_Py_.push_back(genlep.py());
      genlep_Pz_.push_back(genlep.pz());
      genlep_Pt_.push_back(genlep.pt());
      genlep_Eta_.push_back(genlep.eta());
      genlep_Rap_.push_back(genlep.rapidity());
      genlep_Phi_.push_back(genlep.phi());
      genlep_En_.push_back(genlep.energy());

      // gen Post FSR
      if(abs(genlep.pdgId())==11 && genlep.fromHardProcessFinalState()==1){
	
	genPostFSR_Pt_.push_back(genlep.pt());
	genPostFSR_Eta_.push_back(genlep.eta());
	genPostFSR_Phi_.push_back(genlep.phi());
	genPostFSR_En_.push_back(genlep.energy());

      }

      if(abs(genlep.pdgId())==22 && abs(genlep.mother(0)->pdgId())==11){

	genPhoton_Pt_.push_back(genlep.pt());
	genPhoton_Eta_.push_back(genlep.eta());
	genPhoton_Phi_.push_back(genlep.phi());
	genPhoton_En_.push_back(genlep.energy());
      }

      // Separating the taus coming from Z decay
      if(abs(genlep.pdgId())==15 && genlep.fromHardProcessDecayed()==1) nGenTaus++;
    }

    // Select the events containing 2 taus from hard-process
    if(nGenTaus == 2) tauFlag = 1;

  } // isMC && isSIG

  //cout<<"5"<<endl;

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  if( isAOD )
    iEvent.getByToken(vtxToken_, vertices);
  else
    iEvent.getByToken(vtxMiniAODToken_, vertices);

  if (vertices->empty()) return; // skip the event if no PV found
  nPV = vertices->size();
  int firstGoodVertexIdx = 0;

  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    firstGoodVertex = vtx;
    break;
  }

  if (firstGoodVertex==vertices->end())
    return; // skip event if there are no good PVs

  // Seems always zero. Not stored in miniAOD...?
  pvNTracks = firstGoodVertex->nTracks();

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

  nElectrons = 0;

  // Loop over electrons
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);

    // Kinematics

    if(el->pt() > 10. && el->eta() < 2.7){

      nElectrons++;
      ptElec_.push_back( el->pt() );
      etaElec_.push_back( el->eta() );
      rapElec_.push_back( el->rapidity() );
      phiElec_.push_back( el->phi() );
      energyElec_.push_back( el->energy() );
      massElec_.push_back( el->mass() );
      chargeElec_.push_back( el->charge() );

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
      isoRho_.push_back((pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * eA))/(el->pt()));

      // Track - Impact Parameter, Conversion rejection, Converted
      reco::GsfTrackRef theTrack = el->gsfTrack();
      d0_.push_back( (-1) * theTrack->dxy(firstGoodVertex->position() ) );
      dz_.push_back( theTrack->dz( firstGoodVertex->position() ) );
      expectedMissingInnerHits_.push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );
      bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, theBeamSpot->position());
      passConversionVeto_.push_back( (int) passConvVeto );
      brem_.push_back(el->fbrem());

      // Electron ID
      bool isPassVeto  = (*veto_id_decisions)[el];
      bool isPassLoose  = (*loose_id_decisions)[el];
      bool isPassMedium = (*medium_id_decisions)[el];
      bool isPassTight  = (*tight_id_decisions)[el];
      passVetoId_.push_back  ( (int)isPassVeto  );
      passLooseId_.push_back ( (int)isPassLoose );
      passMediumId_.push_back( (int)isPassMedium);
      passTightId_.push_back ( (int)isPassTight );

      eleInBarrel_.push_back(el->isEB());
      eleInEndcap_.push_back(el->isEE());

      // ECAL driven
      eleEcalDrivenSeed_.push_back(el->ecalDrivenSeed());

    }
  }

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonsMiniAODToken_,muons);

  nMuons = 0;

  for(unsigned i=0; i < muons->size();++i ) {
    const auto mu = muons->ptrAt(i);

    if(mu->pt() > 10. && mu->eta() < 2.7){
      nMuons++;

      // Kinematics
      ptMuon_.push_back(mu->pt());
      etaMuon_.push_back(mu->eta());
      phiMuon_.push_back(mu->phi());
      energyMuon_.push_back(mu->energy());
      chargeMuon_.push_back(mu->charge());

      // isLoose
      bool isLooseId = muons->at(i).isLooseMuon();
      isLooseMuon_.push_back(isLooseId);

      // isTight
      bool isTightId = muons->at(i).isTightMuon(vertices->at(0));
      isTightMuon_.push_back(isTightId);

      // is High Pt
      bool isHighPt_Id = muons->at(i).isHighPtMuon(vertices->at(0));
      isHighPtMuon_.push_back(isHighPt_Id);

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

  // MET
  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  const pat::MET &met = metHandle->front();
  metPt_.push_back(met.pt());
  metPhi_.push_back(met.phi());
  metSumEt_.push_back(met.sumEt());

  edm::Handle<edm::View<reco::Photon> > photons;
  iEvent.getByToken(photonsMiniAODToken_,photons);

  nPhotons = 0;

  // Loop over photons
  for (size_t i = 0; i < photons->size(); ++i){
    const auto pho = photons->ptrAt(i);

    // Kinematics
    if(pho->pt() < 15.){
      nPhotons++;

      // Kinematics
      ptPhoton_  .push_back(pho->pt());
      etaPhoton_ .push_back(pho->superCluster()->eta());
      phiPhoton_ .push_back(pho->superCluster()->phi());
    }
  }

  //cout<<"4"<<endl;

  // Save this electron's info
  electronTree_->Fill();

  hlNames_.clear();
  single_electron_23_triggers_in_run.clear();
  single_electron_27_triggers_in_run.clear();
  single_muon_triggers_in_run.clear();
  double_electron_triggers_in_run.clear();
  double_emu_17_8_triggers_in_run.clear();
  double_emu_12_17_triggers_in_run.clear();
  double_emu_23_8_triggers_in_run.clear();

  idx_singleEle23.clear();
  idx_singleEle27.clear();
  idx_singleMuon.clear();
  idx_doubleElectron.clear();
  idx_doubleEMu_17_8.clear();
  idx_doubleEMu_12_17.clear();
  idx_doubleEMu_23_8.clear();

  // Clear vectors
  pt_Ele23.clear();
  eta_Ele23.clear();
  phi_Ele23.clear();

  etPhoton_.clear();

  if(misMC && misSIG){

    genlep_Id_.clear();
    genlepMother_Id_.clear();
    fromHProcessFinalState_.clear();
    fromHProcessDecayed_.clear();

    genlep_Pt_.clear();
    genlep_Px_.clear();
    genlep_Py_.clear();
    genlep_Pz_.clear();
    genlep_Eta_.clear();
    genlep_Rap_.clear();
    genlep_Phi_.clear();
    genlep_En_.clear();

    genPostFSR_Pt_.clear();
    genPostFSR_Eta_.clear();
    genPostFSR_Phi_.clear();
    genPostFSR_En_.clear();
    
    genPhoton_Pt_.clear();
    genPhoton_Eta_.clear();
    genPhoton_Phi_.clear();
    genPhoton_En_.clear();
  }

  ptElec_.clear();
  etaElec_.clear();
  rapElec_.clear();
  phiElec_.clear();
  energyElec_.clear();
  massElec_.clear();
  chargeElec_.clear();

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
  eleEcalDrivenSeed_.clear();

  eleInBarrel_.clear();
  eleInEndcap_.clear();

  isLooseMuon_.clear();
  isTightMuon_.clear();
  isHighPtMuon_.clear();

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

  ptPhoton_.clear();
  etaPhoton_.clear();
  phiPhoton_.clear();

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

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleElectronNtupler);
