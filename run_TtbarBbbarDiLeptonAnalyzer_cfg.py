import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarBbbarDiLeptonAnalyzer")

#process.Tracer = cms.Service("Tracer") 
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames = ['file:/xrootd/store/group/CAT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/v8-0-4_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170113_045423/0000/catTuple_1.root',]
#process.source.fileNames = ['file:/xrootd/store/group/CAT/DoubleMuon/v8-0-4_Run2016D-23Sep2016-v1/170113_134243/0000/catTuple_1.root',]
#process.source.fileNames = ['file:/xrootd/store/group/CAT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_103557/0000/catTuple_1.root']
#process.source.fileNames = ['file:/xrootd/store/group/CAT/chanwook/private/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/private_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180807_064829/0000/catTuple_1.root']
process.source.fileNames = ['file:/xrootd/store/group/CAT/DoubleMuon/v8-0-6_Run2016B-03Feb2017_ver2-v2/171207_054751/0000/catTuple_2.root']
#process.source.fileNames = ['file:/xrootd/store/group/CAT/v8-0-8/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/171126_161209/0000/catTuple_1.root',]
#process.source.fileNames = ['file:/xrootd/store/group/CAT/DoubleMuon/v8-0-6_Run2016D-03Feb2017-v1/170303_095442/0000/catTuple_1.root']
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-5/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-5/DoubleEG_Run2015D-16Dec2015-v2.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-5/DoubleMuon_Run2015D-16Dec2015-v1.root',]
#process.source.fileNames = ['/store/user/jhgoh/CATTools/sync/v7-6-5/MuonEG_Run2015D-16Dec2015-v1.root',]

#import os
#useGold = True
#isRunData = False
lumiMask = 'lumiMask'
#if useGold:
#    catmet = 'catMETs'
#    if isRunData:
#        #lumiFile = 'Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#        lumiFile = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#        from FWCore.PythonUtilities.LumiList import LumiList
#        lumiList = LumiList(os.environ["CMSSW_BASE"]+'/src/CATTools/CatProducer/prod/LumiMask/'+lumiFile)
#        process.source.lumisToProcess = lumiList.getVLuminosityBlockRange()
    
#    process.load("CATTools.CatProducer.pileupWeight_cff")
#    from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
#    process.pileupWeight.weightingMethod = "RedoWeight"
#    process.pileupWeight.pileupRD = pileupWeightMap["Run2015_25nsV1"]
#    process.pileupWeight.pileupUp = pileupWeightMap["Run2015Up_25nsV1"]
#    process.pileupWeight.pileupDn = pileupWeightMap["Run2015Dn_25nsV1"]

#process.load("CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionAlgos_cff")


#import os
#lumiFile = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
#from FWCore.PythonUtilities.LumiList import LumiList
#lumiList = LumiList(os.environ["CMSSW_BASE"]+'/src/CATTools/CatProducer/data/LumiMask/'+lumiFile)
#process.source.lumisToProcess = lumiList.getVLuminosityBlockRange()


####for running genTop on the fly. however it is running slowly.
#process.load("CATTools.CatProducer.genTopProducer_cfi")
#from CATTools.CatProducer.Tools.tools import genHFTool
#genHFTool(process,True)

process.load("CATTools.CatAnalyzer.filters_cff")

##for only ttbar signal mc sample
process.load("CATTools.CatProducer.mcTruthTop.partonTop_cfi")
process.load("CATTools.CatAnalyzer.topPtWeightProducer_cfi")
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")

from CATTools.CatAnalyzer.leptonSF_cff import *

## Redo the pileup weight - necessary for v765 production
'''
process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeight = process.pileupWeight.clone()
from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
process.redoPileupWeight.weightingMethod = "RedoWeight"
process.redoPileupWeight.pileupMC = pileupWeightMap["2016_25ns_SpringMC"]
process.redoPileupWeight.pileupRD = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON"]
process.redoPileupWeight.pileupUp = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Up"]
process.redoPileupWeight.pileupDn = pileupWeightMap["Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Dn"]
'''
process.cattree = cms.EDAnalyzer("TtbarBbbarDiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    recoFiltersMC = cms.InputTag("filterRECOMC"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweights = cms.InputTag("flatGenWeights", "pdf"),
    scaleupweights = cms.InputTag("flatGenWeights", "scaleup"),
    scaledownweights = cms.InputTag("flatGenWeights", "scaledown"),
    topPtWeight = cms.InputTag("topPtWeight"),

    lumiSelection = cms.InputTag(lumiMask),
    puweight = cms.InputTag("pileupWeight"),
    puweightUp = cms.InputTag("pileupWeight","up"),
    #puweightDown = cms.InputTag("pileupWeight","dn"),
    puweightDown = cms.InputTag("pileupWeight","dn"), ## Redo the pileup weight for downward variation in v765 prod
    #isRunH = cms.bool(True),
    isRunH = cms.bool(False),

    #DiLepton and SingleLepton
    trigMUELRunBtoG = cms.InputTag("filterTrigMUELRunBtoG"),
    trigMUELRunH = cms.InputTag("filterTrigMUELRunH"),
    trigMUEL = cms.InputTag("filterTrigMUELRunBtoG"),

    trigMUMURunBtoG = cms.InputTag("filterTrigMUMURunBtoG"),
    trigMUMURunH = cms.InputTag("filterTrigMUMURunH"),
    trigMUMU = cms.InputTag("filterTrigMUMURunBtoG"),

    trigELEL = cms.InputTag("filterTrigELEL"),

    #DiLepton 
    trigMUELRunBtoGDiLepton = cms.InputTag("filterTrigMUELRunBtoGDiLepton"),
    trigMUELRunHDiLepton = cms.InputTag("filterTrigMUELRunHDiLepton"),
    trigMUELDiLepton = cms.InputTag("filterTrigMUELRunBtoGDiLepton"),

    trigMUMURunBtoGDiLepton = cms.InputTag("filterTrigMUMURunBtoGDiLepton"),
    trigMUMURunHDiLepton = cms.InputTag("filterTrigMUMURunHDiLepton"),
    trigMUMUDiLepton = cms.InputTag("filterTrigMUMURunBtoGDiLepton"),

    trigELELDiLepton = cms.InputTag("filterTrigELELDiLepton"),

    #SigleLepton 
    trigMUELRunBtoGSingleLepton = cms.InputTag("filterTrigMUELRunBtoGSingleLepton"),
    trigMUELRunHSingleLepton = cms.InputTag("filterTrigMUELRunHSingleLepton"),
    trigMUELSingleLepton = cms.InputTag("filterTrigMUELRunBtoGSingleLepton"),

    trigMUELRunBtoGSingleElectron = cms.InputTag("filterTrigMUELRunBtoGSingleElectron"),
    trigMUELRunHSingleElectron = cms.InputTag("filterTrigMUELRunHSingleElectron"),
    trigMUELSingleElectron = cms.InputTag("filterTrigMUELRunBtoGSingleElectron"),

    trigMUELRunBtoGSingleMuon = cms.InputTag("filterTrigMUELRunBtoGSingleMuon"),
    trigMUELRunHSingleMuon = cms.InputTag("filterTrigMUELRunHSingleMuon"),
    trigMUELSingleMuon = cms.InputTag("filterTrigMUELRunBtoGSingleMuon"),

    trigMUMURunBtoGSingleLepton = cms.InputTag("filterTrigMUMURunBtoGSingleLepton"),
    trigMUMURunHSingleLepton = cms.InputTag("filterTrigMUMURunHSingleLepton"),
    trigMUMUSingleLepton = cms.InputTag("filterTrigMUMURunBtoGSingleLepton"),

    trigELELSingleLepton = cms.InputTag("filterTrigELELSingleLepton"),

    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),

    #elecSF = electronSFWP90,
    elecSF = combineSF(electronSFRecoOnly,electronSFCutBasedIDTightWPIdOnly),

    muonSF = combineSF(combineSF(muonSFTightIdOnly, muonSFTightIsoOnly), muonSFTrackingOnly),

    muonSFGH = combineSF(combineSF(muonSFTightGHIdOnly, muonSFTightGHIsoOnly), muonSFTrackingGHOnly),


    genTtbarId = cms.InputTag("GenTtbarCategories", "genTtbarId"),
    #genTtbarId30 = cms.InputTag("GenTtbarCategories30", "genTtbarId"),
    #genTtbarId40 = cms.InputTag("GenTtbarCategories40", "genTtbarId"),

    #GenJets = cms.InputTag("slimmedGenJets"),
    #GenParticles = cms.InputTag("prunedGenParticles"),
    GenTop = cms.InputTag("catGenTops"),
    #GenTop = cms.InputTag("catGenTops","","TtbarDiLeptonAnalyzer"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree.root"
))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000


