#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

#analysis = 'h2muAnalyzer'
#analysis = 'ttbbAnalyzer'
#analysis = 'triggerEfficiency'
#analysis = 'TTBBGenAnalyzer'
analysis = 'TtbarLeptonJetsAnalyzer'
#analysis = 'plotVsData'

pythonCfg = 'run_'+analysis+'_cfg.py'
#analysis=analysis+'Silver'

import os,json
datadir = os.environ["CMSSW_BASE"]+'/src/CATTools/CatAnalyzer/data/dataset/'
dataset_json = datadir + 'dataset.json'
#dataList = ['TT_powheg', 'WJets','WW','WWW','WWZ','WZ','WZZ','ZZ','ZZZ','DYJets','DYJets_10to50','SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','ttW','ttZ']
#dataList = ['TT_powheg', 'WJets','WW','WWW','WWZ','WZ','WZZ','ZZ','ZZZ','DYJets','DYJets_10to50', 'SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW']
#dataList = ['DYJets_10to50']
dataList = ['SingleElectron_Run2016B','SingleElectron_Run2016C','SingleElectron_Run2016D','SingleElectron_Run2016E','SingleElectron_Run2016F','SingleElectron_Run2016G','SingleElectron_Run2016H_v2','SingleElectron_Run2016H_v3','SingleMuon_Run2016B','SingleMuon_Run2016C','SingleMuon_Run2016D','SingleMuon_Run2016E','SingleMuon_Run2016F','SingleMuon_Run2016G','SingleMuon_Run2016H_v2','SingleMuon_Run2016H_v3','TT_powheg', 'TT_powheg_UEup','TT_powheg_UEdown','TT_powheg_FSRup','TT_powheg_FSRdown','TT_powheg_ISRup','TT_powheg_ISRdown','TTJets_aMC', 'TT_powheg_herwig', 'WJets','WW','WZ','ZZ','DYJets','DYJets_10to50','SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','QCD_Pt-15to20_MuEnriched', 'QCD_Pt-20to30_MuEnriched', 'QCD_Pt-30to50_MuEnriched', 'QCD_Pt-50to80_MuEnriched', 'QCD_Pt-80to120_MuEnriched', 'QCD_Pt-120to170_MuEnriched', 'QCD_Pt-170to300_MuEnriched', 'QCD_Pt-300to470_MuEnriched', 'QCD_Pt-470to600_MuEnriched', 'QCD_Pt-600to800_MuEnriched', 'QCD_Pt-800to1000_MuEnriched', 'QCD_Pt-1000toInf_MuEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched', 'QCD_Pt-80to120_EMEnriched', 'QCD_Pt-120to170_EMEnriched', 'QCD_Pt-170to300_EMEnriched', 'QCD_Pt-300toInf_EMEnriched']
#dataList = ['SingleElectron_Run2016B','SingleElectron_Run2016C','SingleElectron_Run2016D','SingleElectron_Run2016E','SingleElectron_Run2016F','SingleElectron_Run2016G','SingleElectron_Run2016H_v2','SingleElectron_Run2016H_v3','SingleMuon_Run2016B','SingleMuon_Run2016C','SingleMuon_Run2016D','SingleMuon_Run2016E','SingleMuon_Run2016F','SingleMuon_Run2016G','SingleMuon_Run2016H_v2','SingleMuon_Run2016H_v3','MuonEG_Run2016B','MuonEG_Run2016C','MuonEG_Run2016D','MuonEG_Run2016E','MuonEG_Run2016F','MuonEG_Run2016G','MuonEG_Run2016H_v2','MuonEG_Run2016H_v3']
#dataList += ['WWW','WWZ','WZZ','ZZZ']
#dataList = ['TT_powheg','WJets','WW','WZ','ZZ','DYJets','DYJets_10to50','SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','ttW','ttZ']
#dataList = ['DoubleEG_Run2016G', 'DoubleEG_Run2016H_v2', 'SingleTbar_tW', 'SingleTop_tW']
#dataList = ['MET_Run2016B','MET_Run2016C','MET_Run2016D','MET_Run2016E','MET_Run2016F','MET_Run2016G','MET_Run2016H_v2','MET_Run2016H_v3']
baseDir = '/store/user/chanwook/ttbb/LepJets'
#dataList = ['DoubleEG_Run2016D']
#baseDir = '.'
#os.system('rm -rf /xrootd/%s/*' % (baseDir))
with open(dataset_json) as data_file:    
    data = json.load(data_file)
    for i in data:
        #print data[0]
        datasetName = i['name']
        if datasetName not in dataList: continue
        #if 'DoubleMuon' not in datasetName and 'DoubleEG' not in datasetName and 'MuonEG' not in datasetName: continue
        #if "QCD" in datasetName:
        #    continue
        if "ttH" in datasetName:
            continue       
        transDest = '%s/%s' % (baseDir, datasetName)
        #os.system('mkdir /xrootd/%s' % (transDest))
        fileList = datadir + 'dataset_' + datasetName + '.txt'
        jobName = analysis+'_'+datasetName
        createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 10 --transferDest %s"%(pythonCfg, jobName, fileList, transDest)
        #createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 10 "%(pythonCfg, jobName, fileList)
        print createbatch
        os.system(createbatch)
        
