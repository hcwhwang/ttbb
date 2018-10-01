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
def getEntries(step, rootfile, tree, cut):
  if step <= 3: return getWeightedEntries(rootfile, tree, "filtered", cut) 
  else: return getWeightedEntries2D(rootfile, tree, "jets_bDiscriminatorCSV[csvd_jetid[2]]:jets_bDiscriminatorCSV[csvd_jetid[3]]", cut) 
analysis = "TtbarBbbarDiLeptonAnalyzer_"

#lepton1_pt>20 && lepton2_pt>20lepton1_pt>20 && lepton2_pt>20visible="(NJets20>=6 && NbJets20>=2 && (((lepton1_pt>20 && abs(lepton1_eta)<2.4)) || (lepton2_pt>20 && abs(lepton2_eta)<2.4)))"
visible="(NJets30>=4 && NbJets30>=2 && ((lepton1_pt>25 && lepton2_pt>20)||(lepton1_pt>20 && lepton2_pt>25)) && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"

ttbb = mAND("(NbJets30>=4)",visible)
ttbj = mAND("(NbJets30==3)",visible)
ttcc = mAND("((NcJets30>=2) && !(NbJets30>=3))",visible)
ttlf = mAND("(!(NbJets30>=4) && !(NbJets30==3) && !(NcJets30>=2))",visible)
ttothers = op_(visible)
full ="(NaddJets30 >= 2)"
ttbb_full = "(NaddbJets30 >= 2)"
ttbj_full = "(NaddJets30 >= 2 && NaddbJets30 == 1)"
ttcc_full = "(NaddJets30 >= 2 && NaddcJets30 >= 2 && NaddbJets30==0)"
ttlf_full = "( !"+ttbb_full+" && !"+ttbj_full+" && !"+ttcc_full+"  && NaddJets30 >= 2)"
ttothers_full = op_(full)
ttjj_cut = [ttbb, ttbj, ttcc, ttlf, ttothers, visible]
#ttjj_cut = [ttbb_full, ttbj_full, ttcc_full, ttlf_full, ttothers_full, full]
ttjj_title = ["ttbb", "ttbj", "ttcc", "ttlf", "ttothers", "ttjj"]


datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/totdataset.json" % os.environ["CMSSW_BASE"]))
#rootdir = "%s/src/CATTools/CatAnalyzer/test" % os.environ["CMSSW_BASE"]
#rootdir = "/xrootd/store/user/chanwook/ttbb/final/preapproval/"
rootdir = "./"
datalumi = 35.82*1000
ChList = ['MuEl','ElEl','MuMu']
RDList = ['MuonEG', 'DoubleEG', 'DoubleMuon']
#MCList = [ "WJets",'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'ttW','ttZ']
DYList = ['DYJets', 'DYJets_10to50']
MCList = [ 'WJets', 'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW','WWZ','WZZ','ZZZ','TTW','TTZ','ttH_bb']
TTList = ['TT_powheg_tot']
treeN = "cattree/nom"
tree = "cattree/nom2"
weight = 'tri*weight*puweight*mueffweight*eleffweight*csvweights2[0]'
weight1 = 'tri*weight*puweight*mueffweight*eleffweight'
weightN = 'weight'
var = "jets_bDiscriminatorCSV[csvd_jetid[2]]:jets_bDiscriminatorCSV[csvd_jetid[3]]"
output = []
dyscale = json.load(open("drellyanresult.json"))
for step in range(1,6):
  for i,d in enumerate(RDList):
    print "step = ", step
    print "channel = ", ChList[i]
    nameList = []
    nList = []
    errList = []
    if step == 5: chCut = "filtered == 1 && tri>0 && channel==%s && nbjetM30>=2" % (i+1)
    else :chCut = "filtered == 1 && tri>0 && channel==%s" % (i+1)
    stepCut = "step>=%s" % step
    cut = mAND(chCut, stepCut)
    if step<4: totW = weight1
    else: totW = weight
    weightedCut = "(%s)*(%s)" % (cut, totW)
   
    print RDList[i]
    rdfile = rootdir + RDList[i]+'_Run2016_tot.root'
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

    for k,tt in enumerate(TTList):
      data = findDataSet('TT_powheg', datasets)
      ttfile = rootdir + TTList[k]+'.root'
      wentries1 = getWeightedEntries(rootdir+"TT_powheg.root", treeN, "filtered", weightN)
      wentries2 = getWeightedEntries(rootdir+"TTLL_powheg.root", treeN, "filtered", weightN)

      scale = datalumi*data["xsec"]/(wentries1 + wentries2*(1./0.11))
      #totN = getWeightedEntries(ttfile, treeN, "filtered", weightN)
      #scale = datalumi*data["xsec"]/totN
      for l,cate in enumerate(ttjj_title):
        print cate
        ttCut = mAND(ttjj_cut[l], cut)
        ttWeightedCut = "(%s)*(%s)" % (ttCut, weight)
        n = getWeightedEntries(ttfile, tree, "filtered", ttWeightedCut, scale)
        err = getWeightedEntriesError(ttfile, tree, "filtered", ttWeightedCut, scale)
        print n, " ", err
        print ""
        nameList.append(cate)
        nList.append(n)
        errList.append(err)
    outputSub = {"ChStep":"%s%s"%(ChList[i],step), "name":nameList, "entries":nList, "error":errList}    
    output.append(outputSub)
    print ""
    print ""
print output
file = open("entriesPA.json", "w")
jsonoutput = json.dumps(output, indent=4)
file.write(jsonoutput)
file.close()
