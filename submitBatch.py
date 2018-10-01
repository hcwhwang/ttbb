#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

#analysis = 'h2muAnalyzer'
#analysis = 'ttbbAnalyzer'
#analysis = 'triggerEfficiency'
#analysis = 'TTBBGenAnalyzer'
analysis = 'TtbarBbbarDiLeptonJESAnalyzer'
#analysis = 'TtbarBbbarDiLeptonAnalyzer'
#analysis = 'TtbarBbbarDiLeptonDoubleRDAnalyzer'
#analysis = 'TtbarBbbarDiLeptonSingleMuonAnalyzer'
#analysis = 'TtbarBbbarDiLeptonSingleElectronAnalyzer'
#analysis = 'TtbarBbbarDiLeptonPlotAnalyzer'
#analysis = 'plotVsData'

pythonCfg = 'run_'+analysis+'_cfg.py'
#analysis=analysis+'Silver'

import os,json
datadir = os.environ["CMSSW_BASE"]+'/src/CATTools/CatAnalyzer/data/dataset/'
dataset_json = datadir + 'dataset.json'
#dataList = ['TT_powheg', 'WJets','WW','WWW','WWZ','WZ','WZZ','ZZ','ZZZ','DYJets','DYJets_10to50','SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','ttW','ttZ']
#dataList = ['TT_powheg', 'WJets','WW','WWW','WWZ','WZ','WZZ','ZZ','ZZZ','DYJets','DYJets_10to50', 'SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW']
#dataList = ['DYJets_10to50']
#dataList = ['ttH_bb', 'TT_powheg', 'TT_powheg_up_part1', 'TT_powheg_up_part2','TT_powheg_down_part1','TT_powheg_down_part2','TT_powheg_fsrup_part1','TT_powheg_fsrup_part2','TT_powheg_fsrdown_part1','TT_powheg_fsrdown_part2','TT_powheg_isrup','TT_powheg_isrdown_part1','TT_powheg_isrdown_part2','TT_powheg_hdampup_part1','TT_powheg_hdampup_part2','TT_powheg_hdampdown_part1','TT_powheg_hdampdown_part2','TT_powheg_erdon_part1','TT_powheg_erdon_part2','TT_powheg_qcderdon_part1','TT_powheg_qcderdon_part2','TT_powheg_gluonmove','WJets_part1','WJets_part2','WW_part1','WW_part2','WZ_part1','WZ_part2','ZZ_part1','ZZ_part2','DYJets','DYJets_10to50_part1', 'DYJets_10to50_part2', 'DYJets_10to50_part3', 'SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','TTW_part1','TTW_part2','TTZ_part1','TTZ_part2','TTZ_part3','WWW','WWZ','WZZ','ZZZ','TTLL_powheg']
#dataList = ['SingleElectron_Run2016B','SingleElectron_Run2016C','SingleElectron_Run2016D','SingleElectron_Run2016E','SingleElectron_Run2016F','SingleElectron_Run2016G','SingleMuon_Run2016B','SingleMuon_Run2016C','SingleMuon_Run2016D','SingleMuon_Run2016E','SingleMuon_Run2016F','SingleMuon_Run2016G','DoubleEG_Run2016B','DoubleEG_Run2016C','DoubleEG_Run2016D','DoubleEG_Run2016E','DoubleEG_Run2016F','DoubleEG_Run2016G','DoubleMuon_Run2016B','DoubleMuon_Run2016C','DoubleMuon_Run2016D','DoubleMuon_Run2016E','DoubleMuon_Run2016F','DoubleMuon_Run2016G','MuonEG_Run2016B','MuonEG_Run2016C','MuonEG_Run2016D','MuonEG_Run2016E','MuonEG_Run2016F','MuonEG_Run2016G']
#dataList = ['DoubleMuon_Run2016H_v2', 'DoubleMuon_Run2016H_v3', 'MuonEG_Run2016H_v2', 'MuonEG_Run2016H_v3','DoubleEG_Run2016H_v2', 'DoubleEG_Run2016H_v3','SingleElectron_Run2016H_v2','SingleElectron_Run2016H_v3','SingleMuon_Run2016H_v2','SingleMuon_Run2016H_v3',]
dataList = ['TT_powheg','TTLL_powheg', 'WJets_part1','WJets_part2','WW_part1','WW_part2','WZ_part1','WZ_part2','ZZ_part1','ZZ_part2','DYJets','DYJets_10to50_part1', 'DYJets_10to50_part2', 'DYJets_10to50_part3', 'SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','TTW_part1','TTW_part2','TTZ_part1','TTZ_part2','TTZ_part3','WWW','WWZ','WZZ','ZZZ']
#dataList = ['TT_powheg_ISRdown','TTJets_aMC']
#dataList = ['TT_powheg','WJets','WW','WZ','ZZ','DYJets','DYJets_10to50','SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','ttW','ttZ']
#dataList = ['DoubleEG_Run2016G', 'DoubleEG_Run2016H_v2', 'SingleTbar_tW', 'SingleTop_tW']
#dataList = ['MET_Run2016B','MET_Run2016C','MET_Run2016D','MET_Run2016E','MET_Run2016F','MET_Run2016G','MET_Run2016H_v2','MET_Run2016H_v3']
#dataList = ['TTLL_powheg','TT_powheg']
#dataList = ['DoubleEG_Run2016B','DoubleEG_Run2016C','DoubleEG_Run2016D','DoubleEG_Run2016E','DoubleEG_Run2016F','DoubleEG_Run2016G','DoubleMuon_Run2016B','DoubleMuon_Run2016C','DoubleMuon_Run2016D','DoubleMuon_Run2016E','DoubleMuon_Run2016F','DoubleMuon_Run2016G','MuonEG_Run2016B','MuonEG_Run2016C','MuonEG_Run2016D','MuonEG_Run2016E','MuonEG_Run2016F','MuonEG_Run2016G']
#dataList = ['DoubleMuon_Run2016H_v2', 'DoubleMuon_Run2016H_v3', 'MuonEG_Run2016H_v2', 'MuonEG_Run2016H_v3','DoubleEG_Run2016H_v2', 'DoubleEG_Run2016H_v3']
#dataList = ['SingleElectron_Run2016B','SingleElectron_Run2016C','SingleElectron_Run2016D','SingleElectron_Run2016E','SingleElectron_Run2016F','SingleElectron_Run2016G']
#dataList = ['SingleElectron_Run2016H_v2','SingleElectron_Run2016H_v3']
#dataList = ['SingleMuon_Run2016B','SingleMuon_Run2016C','SingleMuon_Run2016D','SingleMuon_Run2016E','SingleMuon_Run2016F','SingleMuon_Run2016G']
#dataList = ['SingleMuon_Run2016H_v2','SingleMuon_Run2016H_v3']
#dataList = ['TT_powheg_up_part1', 'TT_powheg_up_part2','TT_powheg_down_part1','TT_powheg_down_part2','TT_powheg_fsrup_part1','TT_powheg_fsrup_part2','TT_powheg_fsrdown_part1','TT_powheg_fsrdown_part2','TT_powheg_isrup','TT_powheg_isrdown_part1','TT_powheg_isrdown_part2','TT_powheg_hdampup_part1','TT_powheg_hdampup_part2','TT_powheg_hdampdown_part1','TT_powheg_hdampdown_part2','TT_powheg_erdon_part1','TT_powheg_erdon_part2','TT_powheg_qcderdon_part1','TT_powheg_qcderdon_part2','TT_powheg_gluonmove']
#baseDir = '/store/user/chanwook/ttbb/final/preapproval/test/'
baseDir = '/xrootd_user/chanwook/xrootd/ttbb/final/preapproval/test/'
#baseDir = '/store/user/chanwook/ttbb/final/preapproval/JES/'
#baseDir = '/store/user/chanwook/ttbb/v808'
#dataList = ['DoubleEG_Run2016D']
#baseDir = '.'
#os.system('rm -rf /xrootd/%s/*' % (baseDir))
with open(dataset_json) as data_file:    
    data = json.load(data_file)
    for i in data:
        #print data[0]
        datasetName = i['name']
        if datasetName not in dataList: continue
        #if 'powheg' not in datasetName: continue
        #if 'DoubleMuon' not in datasetName and 'DoubleEG' not in datasetName and 'MuonEG' not in datasetName: continue
        if "QCD" in datasetName:
            continue
        transDest = '%s/%s' % (baseDir, datasetName)
        #os.system('mkdir /xrootd/%s' % (transDest))
        fileList = datadir + 'dataset_' + datasetName + '.txt'
        jobName = analysis+'_'+datasetName
        #if 'TT' in datasetName : createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 1 --transferDest %s"%(pythonCfg, jobName, fileList, transDest)
        createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 3 --transferDest %s"%(pythonCfg, jobName, fileList, transDest)
        #createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 10 "%(pythonCfg, jobName, fileList)
        print createbatch
        os.system(createbatch)
        
