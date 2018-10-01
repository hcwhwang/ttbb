import FWCore.ParameterSet.Config as cms
process = cms.Process("plotVsData")

#process.Tracer = cms.Service("Tracer") 
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames = ['file:/xrootd/store/group/CAT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/v8-0-4_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170113_045423/0000/catTuple_1.root',]
#process.source.fileNames = ['file:/xrootd/store/group/CAT/DoubleMuon/v8-0-4_Run2016D-23Sep2016-v1/170113_134243/0000/catTuple_1.root',]
process.source.fileNames = ['file:/xrootd/store/group/CAT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/v8-0-6_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170303_103557/0000/catTuple_1.root']
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
process.load("CATTools.CatAnalyzer.topPtWeightProducer_cfi")
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")

from CATTools.CatAnalyzer.leptonSF_cff import *

## Redo the pileup weight - necessary for v765 production
process.load("CATTools.CatProducer.pileupWeight_cff")
from CATTools.CatProducer.pileupWeight_cff import pileupWeightMap
process.redoPileupWeightB = process.pileupWeight.clone()
process.redoPileupWeightB.weightingMethod = "RedoWeight"
process.redoPileupWeightB.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.redoPileupWeightB.pileupRD = pileupWeightMap["2016B"]
process.redoPileupWeightB.pileupUp = pileupWeightMap["2016B_Up"]
process.redoPileupWeightB.pileupDn = pileupWeightMap["2016B_Dn"]

process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeightC = process.pileupWeight.clone()
process.redoPileupWeightC.weightingMethod = "RedoWeight"
process.redoPileupWeightC.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.redoPileupWeightC.pileupRD = pileupWeightMap["2016C"]
process.redoPileupWeightC.pileupUp = pileupWeightMap["2016C_Up"]
process.redoPileupWeightC.pileupDn = pileupWeightMap["2016C_Dn"]

process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeightD = process.pileupWeight.clone()
process.redoPileupWeightD.weightingMethod = "RedoWeight"
process.redoPileupWeightD.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.redoPileupWeightD.pileupRD = pileupWeightMap["2016D"]
process.redoPileupWeightD.pileupUp = pileupWeightMap["2016D_Up"]
process.redoPileupWeightD.pileupDn = pileupWeightMap["2016D_Dn"]

process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeightE = process.pileupWeight.clone()
process.redoPileupWeightE.weightingMethod = "RedoWeight"
process.redoPileupWeightE.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.redoPileupWeightE.pileupRD = pileupWeightMap["2016E"]
process.redoPileupWeightE.pileupUp = pileupWeightMap["2016E_Up"]
process.redoPileupWeightE.pileupDn = pileupWeightMap["2016E_Dn"]

process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeightF = process.pileupWeight.clone()
process.redoPileupWeightF.weightingMethod = "RedoWeight"
process.redoPileupWeightF.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.redoPileupWeightF.pileupRD = pileupWeightMap["2016F"]
process.redoPileupWeightF.pileupUp = pileupWeightMap["2016F_Up"]
process.redoPileupWeightF.pileupDn = pileupWeightMap["2016F_Dn"]

process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeightG = process.pileupWeight.clone()
process.redoPileupWeightG.weightingMethod = "RedoWeight"
process.redoPileupWeightG.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.redoPileupWeightG.pileupRD = pileupWeightMap["2016G"]
process.redoPileupWeightG.pileupUp = pileupWeightMap["2016G_Up"]
process.redoPileupWeightG.pileupDn = pileupWeightMap["2016G_Dn"]

process.load("CATTools.CatProducer.pileupWeight_cff")
process.redoPileupWeightH = process.pileupWeight.clone()
process.redoPileupWeightH.weightingMethod = "RedoWeight"
process.redoPileupWeightH.pileupMC = pileupWeightMap["2016_25ns_Moriond17MC"]
process.redoPileupWeightH.pileupRD = pileupWeightMap["2016H"]
process.redoPileupWeightH.pileupUp = pileupWeightMap["2016H_Up"]
process.redoPileupWeightH.pileupDn = pileupWeightMap["2016H_Dn"]

process.cattree = cms.EDAnalyzer("plotVsData",
    recoFilters = cms.InputTag("filterRECO"),
    recoFiltersMC = cms.InputTag("filterRECOMC"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweights = cms.InputTag("flatGenWeights", "pdf"),
    scaleupweights = cms.InputTag("flatGenWeights", "scaleup"),
    scaledownweights = cms.InputTag("flatGenWeights", "scaledown"),
    topPtWeight = cms.InputTag("topPtWeight"),

    lumiSelection = cms.InputTag(lumiMask),
    puweight = cms.InputTag("redoPileupWeightB"),
    puweightUp = cms.InputTag("redoPileupWeightB","up"),
    #puweightDown = cms.InputTag("pileupWeight","dn"),
    puweightDown = cms.InputTag("redoPileupWeightB","dn"), ## Redo the pileup weight for downward variation in v765 prod
    puB = cms.InputTag("redoPileupWeightB"),
    puBUp = cms.InputTag("redoPileupWeightB","up"),
    puBDown = cms.InputTag("redoPileupWeightB","dn"),

    puC = cms.InputTag("redoPileupWeightC"),
    puCUp = cms.InputTag("redoPileupWeightC","up"),
    puCDown = cms.InputTag("redoPileupWeightC","dn"),

    puD = cms.InputTag("redoPileupWeightD"),
    puDUp = cms.InputTag("redoPileupWeightD","up"),
    puDDown = cms.InputTag("redoPileupWeightD","dn"),

    puE = cms.InputTag("redoPileupWeightE"),
    puEUp = cms.InputTag("redoPileupWeightE","up"),
    puEDown = cms.InputTag("redoPileupWeightE","dn"),

    puF = cms.InputTag("redoPileupWeightF"),
    puFUp = cms.InputTag("redoPileupWeightF","up"),
    puFDown = cms.InputTag("redoPileupWeightF","dn"),

    puG = cms.InputTag("redoPileupWeightG"),
    puGUp = cms.InputTag("redoPileupWeightG","up"),
    puGDown = cms.InputTag("redoPileupWeightG","dn"),

    puH = cms.InputTag("redoPileupWeightH"),
    puHUp = cms.InputTag("redoPileupWeightH","up"),
    puHDown = cms.InputTag("redoPileupWeightH","dn"),

    trigMUEL = cms.InputTag("filterTrigMUEL"),
    trigMUMU = cms.InputTag("filterTrigMUMU"),
    trigELEL = cms.InputTag("filterTrigELEL"),
    trigMUELMC = cms.InputTag("filterTrigMUELMC"),

    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),

    #elecSF = electronSFWP90,
    elecSF = combineSF(electronSFRecoOnly,electronSFCutBasedIDMediumWPIdOnly),

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


