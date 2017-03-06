import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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
#process.GlobalTag.globaltag = cms.string("76X_mcRun2_asymptotic_v12")
process.GlobalTag.globaltag = cms.string("76X_dataRun2_v15")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
        calibratedPatElectrons = cms.PSet(
	initialSeed = cms.untracked.uint32(1),
	engineName = cms.untracked.string('TRandom3')
	),
	calibratedElectrons = cms.PSet(
	initialSeed = cms.untracked.uint32(1),
	engineName = cms.untracked.string('TRandom3')
	),
	                                           )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM
    '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM

       #'/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/40000/1EE52939-48D8-E511-8D0A-20CF3019DF17.root',
       
       '/store/data/Run2015D/SingleElectron/MINIAOD/16Dec2015-v1/20001/86E53624-32A7-E511-A28D-0025905A60D0.root',
       #'/store/data/Run2015C_25ns/SinglePhoton/MINIAOD/16Dec2015-v1/50000/02E07C2B-90B2-E511-96BB-0CC47A78A3F4.root'
    
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

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

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

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons')
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')

#
# Configure the ntupler module
#
process.selectedElectrons = cms.EDFilter("PATElectronSelector",
src = cms.InputTag("slimmedElectrons"),
cut = cms.string("pt > 10 && abs(eta)<2.5")
)

# should have no effect on slimmedElectrons
#process.calibratedPatElectrons.isMC = cms.bool(True)

process.ntupler = cms.EDAnalyzer('SimpleElectronNtupler',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 isMC     = cms.untracked.bool(False),
				 isSIG    = cms.untracked.bool(False),
				 trigger  = cms.InputTag("TriggerResults::HLT"),
				 prescale = cms.InputTag("patTrigger"),
				 #pileup   = cms.InputTag("addPileupInfo"),
                                 pileup   = cms.InputTag("slimmedAddPileupInfo"),
				 rho      = cms.InputTag("fixedGridRhoFastjetAll"),
                                 beamSpot = cms.InputTag('offlineBeamSpot'),
                                 eventWeight   = cms.InputTag("generator"),
                                 #
                                 # Objects specific to AOD format
                                 #
                                 #electrons    = cms.InputTag("gedGsfElectrons"),
                                 electrons    = cms.InputTag("calibratedElectrons"),
				 genParticles = cms.InputTag("genParticles"),
                                 vertices     = cms.InputTag("offlinePrimaryVertices"),
                                 conversions  = cms.InputTag('allConversions'),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 muonsMiniAOD    = cms.InputTag("slimmedMuons"),
				 #electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
				 #electronsMiniAOD    = cms.InputTag("selectedElectrons"),
				 electronsMiniAOD    = cms.InputTag("calibratedPatElectrons"),
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
				 # This is a fairly verbose mode if switched on, with full cut flow 
				 # diagnostics for each candidate. Use it in a low event count test job.
				 eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium",),
				 eleIdVerbose = cms.bool(False),
				 objects = cms.InputTag('selectedPatTrigger')
                                 )

process.primaryVertexFilter  = cms.EDFilter("VertexSelector",
      src = cms.InputTag('offlineSlimmedPrimaryVertices'),
      cut = cms.string('!isFake && ndof > 4.0 && position.Rho < 2.0 && abs(z) < 24'),
      filter = cms.bool(True)  ## otherwise it won't filter the events, just produce an empty vertex collection.
      )

process.hcalDDDRecConstants = cms.ESProducer( "HcalDDDRecConstantsESModule",
  appendToDataLabel = cms.string( "" )
)
process.hcalDDDSimConstants = cms.ESProducer( "HcalDDDSimConstantsESModule",
  appendToDataLabel = cms.string( "" )
)

#from PhysicsTools.PatAlgos.preselection_eltau_cff import *
#process.leadingElectron = process.slimmedElectrons.clone(
#    src = 'slimmedElectrons', cut = 'pt>20'
#    )
#process.load("PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi")
#process.leadingElectronRequirement = process.countPatElectrons.clone(minNumber = 1, src = 'leadingElectron')
#process.leadingleptonsequence = cms.Sequence(process.leadingElectron+process.leadingElectronRequirement)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Data_76X.root')
				   #fileName = cms.string('Data.root')
                                   )


process.p = cms.Path(process.selectedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.primaryVertexFilter + process.ntupler)
#process.p = cms.Path(process.selectedElectrons * process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.ntupler)
#process.p = cms.Path(process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.ntupler)
