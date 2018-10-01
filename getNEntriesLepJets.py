#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from histoHelper import *
from ROOT import *
from sysWeight_cfi import *

def mAND(aaa,bbb):
  return "(" +aaa+ " && "+bbb+")"
def op_(aaa):
  return "!(" + aaa + ")"
def weightCut(a, b):
  cut = "(%s)*(%s)" % (a, b)
  return cut

analysis = "TtbarLeptonJetsAnalyzer_"


datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ["CMSSW_BASE"]))
#rootdir = "%s/src/CATTools/CatAnalyzer/test" % os.environ["CMSSW_BASE"]
rootdir = "/xrootd/store/user/chanwook/ttbb/"
datalumi = 36.5*1000
ChList = ['El','Mu']
RDList = ['SingleElectron', 'SingleMuon']
#MCList = [ "WJets",'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'ttW','ttZ']
MCList = [ 'WJets', 'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW','WWZ','WZZ','ZZZ','ttW','ttZ']
#TTList = ['TT_powheg']
treeN = "cattree/nom"
tree = "cattree/nom2"
weight = 'tri*weight*puweight*mueffweight*eleffweight*csvweights2[0]'
weight1 = 'tri*weight*puweight*mueffweight*eleffweight'
weightN = 'weight*puweight'
output = []
dyscale = json.load(open("drellyanresult.json"))
for step in range(1,4):
  for i,d in enumerate(RDList):
    print "step = ", step
    print "channel = ", ChList[i]
    nameList = []
    nList = []
    errList = []
    chCut = "filtered == 1 && tri>0 && channel==%s" % (i+1)
    stepCut = "step>=%s" % step
    cut = mAND(chCut, stepCut)
    if step<2: totW = weight1
    else: totW = weight
    weightedCut = "(%s)*(%s)" % (cut, totW)
   
    print RDList[i]
    rdfile = rootdir + RDList[i]+'_Run2016.root'
    n = getWeightedEntries(rdfile, tree, "filtered", cut)
    err = getWeightedEntriesError(rdfile, tree, "filtered", cut)
    print n, " ", err
    print ""
    nameList.append(RDList[i])
    nList.append(n)
    errList.append(err)

    for j,mc in enumerate(MCList):
      #if i!=0: continue
      data = findDataSet(mc, datasets)
      print mc
      mcfile = rootdir + MCList[j]+'.root'
      mcname = MCList[j]+'.root'
      if mcname not in sorted(os.listdir(rootdir)): continue
      totN = getWeightedEntries(mcfile, treeN, "filtered", weightN)
      scale = datalumi*data["xsec"]/totN
      n = getWeightedEntries(mcfile, tree, "filtered", weightedCut, scale)
      err = getWeightedEntriesError(mcfile, tree, "filtered", weightedCut, scale)
      print n, " ", err
      print ""
      nameList.append(mc)
      nList.append(n)
      errList.append(err)
    '''
    for m,dy in enumerate(DYList):
      data = findDataSet(dy, datasets)
      print dy
      dyfile = rootdir + DYList[m]+'.root'
      dyname = DYList[m]+'.root'
      if dyname not in sorted(os.listdir(rootdir)): continue
      totN = getWeightedEntries(dyfile, treeN, "filtered", weightN)
      scale = datalumi*data["xsec"]/totN
      for dych in dyscale:
        if dych["channel"] == ChList[i]: dyS = dych
      scale *= dyS["step%s" % step] 
      n = getWeightedEntries(dyfile, tree, "filtered", weightedCut, scale)
      err = getWeightedEntriesError(dyfile, tree, "filtered", weightedCut, scale)
      print n, " ", err
      print ""
      nameList.append(dy)
      nList.append(n)
      errList.append(err)
'''
    outputSub = {"ChStep":"%s%s"%(ChList[i],step), "name":nameList, "entries":nList, "error":errList}    
    output.append(outputSub)
    print ""
    print ""
print output
file = open("entries.json", "w")
jsonoutput = json.dumps(output, indent=4)
file.write(jsonoutput)
file.close()
