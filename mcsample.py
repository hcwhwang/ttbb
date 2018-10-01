import json, os, sys
from histoHelper import *

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ["CMSSW_BASE"]))

MCList =  ["TT_powheg", "TTLL_powheg", "WJets",'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW', 'WWZ', 'WZZ', 'ZZZ', 'ttW','ttZ', 'DYJets', 'DYJets_10to50']
RDList = [['MuonEG_Run2016B','MuonEG_Run2016C','MuonEG_Run2016D','MuonEG_Run2016E','MuonEG_Run2016F','MuonEG_Run2016G','MuonEG_Run2016H_v2','MuonEG_Run2016H_v3'],['DoubleEG_Run2016B','DoubleEG_Run2016C','DoubleEG_Run2016D','DoubleEG_Run2016E','DoubleEG_Run2016F','DoubleEG_Run2016G','DoubleEG_Run2016H_v2','DoubleEG_Run2016H_v3'],['DoubleMuon_Run2016B','DoubleMuon_Run2016C','DoubleMuon_Run2016D','DoubleMuon_Run2016E','DoubleMuon_Run2016F','DoubleMuon_Run2016G','DoubleMuon_Run2016H_v2','DoubleMuon_Run2016H_v3']]
runRange = ["272007-275376","275657-276283", "276315-276811", "276831-277420", "277772-278808", "278820-280385", "280919-284044",""] 
for i,mc in enumerate(MCList):
  data = findDataSet(mc, datasets)
  fullName = data["DataSetName"].split("/")
  name = fullName[1].replace("_","\_")
  name = "/"+name
  print name, "&", data["xsec"], "\\\\"
print ""
for i,rd in enumerate(RDList):
  for j,run in enumerate(rd):
    data = findDataSet(run, datasets)
    print runRange[j], "&", data["DataSetName"].replace("_","\_"), "\\\\"
  print "\\hline"
print ""
print "dataset & run range & trigger selection"
print "MuonEG & 272007-284044 & Run B-G: HLT\_Mu23\_TrkIsoVVL\_Ele12\_CaloIdL\_TrackIdL\_IsoVL\_v* \\\\"
print "& & Run B-G: HLT\_Mu8\_TrkIsoVVL\_Ele23\_CaloIdL\_TrackIdL\_IsoVL\_v* \\\\" 
print "& & Run H: HLT\_Mu23\_TrkIsoVVL\_Ele12\_CaloIdL\_TrackIdL\_IsoVL\_DZ\_v* \\\\"
print "& & Run H: HLT\_Mu8\_TrkIsoVVL\_Ele23\_CaloIdL\_TrackIdL\_IsoVL\_DZ\_v* \\\\\hline"
print "DoubleEG & 272007-284044 & HLT\_Ele23\_Ele12\_CaloIdL\_TrackIdL\_IsoVL\_DZ\_v* \\\\\hline"
print "DoubleMuon & 272007-284044 & Run B-G: HLT\_Mu17\_TrkIsoVVL\_Mu8\_TrkIsoVVL\_v* \\\\"
print "& & Run B-G: HLT\_Mu17\_TrkIsoVVL\_TkMu8\_TrkIsoVVL\_v* \\\\"
print "& & Run H: HLT\_Mu17\_TrkIsoVVL\_Mu8\_TrkIsoVVL\_DZ\_v* \\\\"
print "& & Run H: HLT\_Mu17\_TrkIsoVVL\_TkMu8\_TrkIsoVVL\_DZ\_v* \\\\"
