import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'                                             #MC
#process.GlobalTag.globaltag = 'GR_E_V49::All'                                                 #Data

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = cms.string("74X_dataRun2_Prompt_v0")                            #Data
process.GlobalTag.globaltag = cms.string("MCRUN2_74_V9A")                                      #MC

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM
    '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root',
    '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/040D9AF7-FB6B-E411-8106-0025907DBA06.root',
    '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04442001-036C-E411-9C90-0025901D42C0.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
    
    #'/DYToEE_M-50_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/MINIAODSIM'
    #'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM' 
    #'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/02DE3B74-6C08-E511-ABE3-0025905A60D0.root'
    #'/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'
    #'/DoubleEG/Run2015B-PromptReco-v1/MINIAOD'
    #'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/244/00000/6A0A8868-4B27-E511-B3F8-02163E011BD1.root'
    #'/SingleElectron/Run2015B-PromptReco-v1/MINIAOD'
    #'/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/244/00000/084C9A66-9227-E511-91E0-02163E0133F0.root'
    #'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'
    #'/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/08B2593F-FAFD-E411-AC92-485B39800C3B.root'
    #'/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM'
    #'/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
    #'/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
    '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/MINIAODSIM'
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
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
	                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff',
			 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

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
				 trigger  = cms.InputTag("TriggerResults::HLT"),
				 pileup   = cms.InputTag("addPileupInfo"),
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
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
				 #
                                 # ID decisions (common to all formats)
                                 #
				 eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
                                 eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
                                 eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
                                 eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
				 #
                                 # ValueMaps with MVA results
                                 #
                                 mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues"),
                                 mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories"),
				 
				 objects = cms.InputTag('selectedPatTrigger')
				 #tagFilterName   = cms.vstring("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter"),
				 #probeFilterName = cms.vstring("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter")
                                 )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('TT.root')
#                                   fileName = cms.string('DoubleEG.root')
#                                   fileName = cms.string('WJetsToLNu.root')
#                                   fileName = cms.string('check.root')
                                   )


process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntupler)
#process.p = cms.Path(process.ntupler)
