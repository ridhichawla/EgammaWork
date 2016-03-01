import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'                                             #MC
#process.GlobalTag.globaltag = 'GR_E_V49::All'                                                 #Data

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = cms.string("74X_dataRun2_v4")
#process.GlobalTag.globaltag = cms.string("74X_dataRun2_reMiniAOD_v0")                                     #Data
#process.GlobalTag.globaltag = cms.string("74X_dataRun2_Prompt_v4")
process.GlobalTag.globaltag = cms.string("74X_mcRun2_asymptotic_v2")                                      #MC

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM
    '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
    
    #'/store/data/Run2015D/SingleElectron/MINIAOD/PromptReco-v4/000/258/159/00000/0EC56452-186C-E511-8158-02163E0146D5.root'
    #'/store/data/Run2015D/SingleElectron/MINIAOD/05Oct2015-v1/10000/00991D45-4E6F-E511-932C-0025905A48F2.root'
    
    #'/store/data/Run2015D/MuonEG/MINIAOD/PromptReco-v4/000/258/159/00000/64914E6C-F26B-E511-B0C8-02163E0142D1.root'
    #'/store/data/Run2015D/MuonEG/MINIAOD/05Oct2015-v2/60000/00D43A12-C573-E511-8F4C-0025905A60E0.root'
    
    #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/2699AF1F-9F6D-E511-BDD8-0025901D4AF0.root'
    '/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/041705A6-6F6F-E511-AC9C-001E6757F1D4.root'
    #'/store/mc/RunIISpring15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/008865AA-596D-E511-92F1-0025905A6110.root'
    #'/store/mc/RunIISpring15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/00087FEB-236E-E511-9ACB-003048FF86CA.root'
    
    )

#
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False
if useAOD == True :
    inputFiles = inputFilesAOD
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
	         #'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#
# Configure the ntupler module
#

process.ntupler = cms.EDAnalyzer('SimpleElectronNtupler',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 isMC     = cms.untracked.bool(True),
				 isSIG    = cms.untracked.bool(True),
				 trigger  = cms.InputTag("TriggerResults::HLT"),
				 prescale = cms.InputTag("patTrigger"),
				 #pileup   = cms.InputTag("addPileupInfo"),
                                 pileup   = cms.InputTag("slimmedAddPileupInfo"),
				 rho      = cms.InputTag("fixedGridRhoFastjetAll"),
                                 beamSpot = cms.InputTag('offlineBeamSpot'),
                                 #
                                 # Objects specific to AOD format
                                 #
                                 electrons    = cms.InputTag("gedGsfElectrons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 vertices     = cms.InputTag("offlinePrimaryVertices"),
                                 conversions  = cms.InputTag('allConversions'),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 muonsMiniAOD    = cms.InputTag("slimmedMuons"),
				 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
				 photonsMiniAOD      = cms.InputTag("slimmedPhotons"),
				 metsMiniAOD         = cms.InputTag("slimmedMETs"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
				 #
                                 # ID decisions (common to all formats)
                                 #
				 eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                 eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                 eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                 eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
				 eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
				 phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
                                 phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
                                 phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
				 #
                                 # ValueMaps with MVA results
                                 #
                                 #mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues"),
                                 #mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories"),
				 # This is a fairly verbose mode if switched on, with full cut flow 
				 # diagnostics for each candidate. Use it in a low event count test job.
				 eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium",),
				 eleIdVerbose = cms.bool(False),
				 
				 objects = cms.InputTag('selectedPatTrigger')
				 #tagFilterName   = cms.vstring("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter"),
				 #probeFilterName = cms.vstring("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter")
                                 )

process.primaryVertexFilter  = cms.EDFilter("VertexSelector",
      src = cms.InputTag('offlineSlimmedPrimaryVertices'),
      cut = cms.string('!isFake && ndof > 4.0 && position.Rho < 2.0 && abs(z) < 24'),
      filter = cms.bool(True)  ## otherwise it won't filter the events, just produce an empty vertex collection.
      )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('DY_check.root')
				   #fileName = cms.string('SingleElectron_Run2015C_05Oct2015-v1.root')
                                   )


process.p = cms.Path(process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.ntupler)
#process.p = cms.Path(process.egmPhotonIDSequence * process.ntupler)
