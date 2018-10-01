import os

analysis = 'TtbarBbbarDiLeptonAnalyzer'
dataList = ['MuonEG_Run2016', 'DoubleEG_Run2016', 'DoubleMuon_Run2016', 'DoubleEG_Run2016B','DoubleEG_Run2016C','DoubleEG_Run2016D','DoubleEG_Run2016E','DoubleEG_Run2016F','DoubleEG_Run2016G','DoubleEG_Run2016H_v2','DoubleEG_Run2016H_v3','DoubleMuon_Run2016B','DoubleMuon_Run2016C','DoubleMuon_Run2016D','DoubleMuon_Run2016E','DoubleMuon_Run2016F','DoubleMuon_Run2016G','DoubleMuon_Run2016H_v2','DoubleMuon_Run2016H_v3','MuonEG_Run2016B','MuonEG_Run2016C','MuonEG_Run2016D','MuonEG_Run2016E','MuonEG_Run2016F','MuonEG_Run2016G','MuonEG_Run2016H_v2','MuonEG_Run2016H_v3', 'TT_powheg', 'WJets','WW','WZ','ZZ','DYJets','DYJets_10to50','SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','ttW','ttZ']

datadir = '%s/src/CATTools/CatAnalyzer/test/' % os.environ["CMSSW_BASE"]

for d in dataList:
  fileName = analysis + '_' + d + '.root'
  fullFileName = datadir + fileName
  if os.path.exists(fileName): continue
  os.system("cp %s ." % fullFileName)
