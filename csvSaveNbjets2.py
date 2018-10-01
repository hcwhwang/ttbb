#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from histoHelper import *
from ROOT import *
#from sysWeight_cfi import *
from sysWeightPDF_cfi import *

def mAND(aaa,bbb):
  return "(" +aaa+ " && "+bbb+")"
def op_(aaa):
  return "!(" + aaa + ")"
def weightCut(a, b):
  cut = "(%s)*(%s)" % (a, b)
  return cut

#analysis = "TtbarBbbarDiLeptonAnalyzer_"
analysis = ""
overallCut = "tri>0 && filtered==1 && step>=5 && nbjetM30>=2"
#overallCut = "filtered==1 && step>=5"
channelMCCut = "(channel==1 || channel==2 || channel==3)"

#visible="(NJets20>=6 && NbJets20>=2 && (((lepton1_pt>20 && abs(lepton1_eta)<2.4)) || (lepton2_pt>20 && abs(lepton2_eta)<2.4)))"
visible="(NJets20>=4 && NbJets20>=2 && lepton1_pt>20 && lepton2_pt>20 && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"

ttbb = mAND("(NbJets20>=4)",visible)
ttbj = mAND("(NbJets20==3)",visible)
ttcc = mAND("((NcJets20>=2) && !(NbJets20>=3))",visible)
ttlf = mAND("(!(NbJets20>=4) && !(NbJets20==3) && !(NcJets20>=2))",visible)
ttothers = op_(visible)

#full ="(diLeptonicM1==1 && NaddJets20 >= 2)"
#ttbb_full = "(NaddbJets20 >= 2 && diLeptonicM1==1)"
#ttb_full = "(NaddJets20 >= 2 && NaddbJets20 == 1 && diLeptonicM1==1 && !(genTtbarId%100==52))"
#tt2b_full = "(NaddJets20 >= 2 && NaddbJets20 == 1 && diLeptonicM1==1 && (genTtbarId%100==52))"
#ttcc_full = "(NaddJets20 >= 2 && NaddcJets20 >= 2 && NaddbJets20==0 && diLeptonicM1==1)"
#ttlf_full = "( !"+ttbb+" && !"+ttb+" && !"+ttcc+"  && NaddJets20 >= 2 && diLeptonicM1==1)"
full ="(NaddJets20 >= 2)"
ttbb_full = "(NaddbJets20 >= 2)"
ttbj_full = "(NaddJets20 >= 2 && NaddbJets20 == 1)"
ttcc_full = "(NaddJets20 >= 2 && NaddcJets20 >= 2 && NaddbJets20==0 )"
ttlf_full = "( !"+ttbb+" && !"+ttbj+" && !"+ttcc+"  && NaddJets20 >= 2 )"
ttothers_full = op_(full)

ttjj_cut = [ttbb, ttbj, ttcc, ttlf, ttothers]
ttjj_full_cut = [ttbb_full, ttbj_full, ttcc_full, ttlf_full, ttothers_full]
ttjj_title = ["ttbb", "ttbj", "ttcc", "ttlf", "ttothers"]

name = "Nbjets2_pdf"
rootFile = 'csv%s.root' % name
outputJson = 'output%s.json' % name

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ["CMSSW_BASE"]))
#rootdir = "%s/src/CATTools/CatAnalyzer/test" % os.environ["CMSSW_BASE"]
rootdir = "/xrootd/store/user/chanwook/ttbb/"

datalumi = 36.8*1000
RDList = ['MuonEG_Run2016', 'DoubleEG_Run2016', 'DoubleMuon_Run2016']
MCList = [ "WJets",'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW', 'WZZ', 'ZZZ', 'ttW','ttZ']
DYList = ['DYJets', 'DYJets_10to50']#MCList = [ 'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW','WWZ','WZZ','ZZZ','ttW','ttZ']
TTList = ['TT_powheg', 'TTJets_aMC', 'TT_powheg_herwig', 'TT_powheg_FSRdown', 'TT_powheg_FSRup', 'TT_powheg_ISRdown', 'TT_powheg_ISRup', 'TT_powheg_UEup', 'TT_powheg_UEdown']
ChList = ['MuEl','ElEl','MuMu']
ChCutList = ['channel==1','channel==2','channel==3']
nomtreepath = "cattree/nom"
output = []
dyscale = json.load(open("drellyanresult.json"))
#binning = [10, 0, 1., 10, 0.5426, 1.]
#binning = [10, 0, 1., 10, 0.8484, 1.]
binning = [10, 0, 1., 10, 0, 1.]
if os.path.exists(rootFile): os.system("rm -f %s" % rootFile)
for i,d in enumerate(mceventweight):
  print d["name"]
  if d["name"]=="PW_Up": PUweight = "weight*puweightUp"
  elif d["name"]=="PW_Down": PUweight = "weight*puweightDown"
  else: PUweight = "weight*puweight"
  tree = d["tree"]
  treepath = "cattree/"+tree
  weight = d["var"] 
  bkgCut = mAND(overallCut, channelMCCut)
  weightedBkgCut = "(%s)*(%s)" % (bkgCut, weight)
  rdCut = [mAND(overallCut, "channel==1"), mAND(overallCut, "channel==2"), mAND(overallCut, "channel==3")]

  var = "jets_bDiscriminatorCSV[csvd_jetid[2]]:jets_bDiscriminatorCSV[csvd_jetid[3]]"

  hBkgAllName = "Bkg_" + tree + "_" + d["name"]
  hBkgAll = TH2F(hBkgAllName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])

  hRDAllName = "Data_" + tree + "_" + d["name"]
  hRDAll = TH2F(hRDAllName,"",binning[0],binning[1], binning[2], binning[3], binning[4], binning[5])
  if i==0:
    for i,file in enumerate(RDList):
      histName = RDList[i]
      print histName
      fileName = rootdir + RDList[i]+'.root'
      f = TFile.Open(fileName)
      t = f.Get(treepath)
      h = TH2F(histName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
      t.Project(histName, var, rdCut[i])
      print h.Integral()
      #if (i==0): hRDAll = h.Clone()
      #else : hRDAll.Add(h)
      hRDAll.Add(h)
  print "data tot = ", hRDAll.Integral()
  hList = []
  for i,file in enumerate(MCList):
    #if i!=0: continue
    histName = MCList[i]+ "_" + tree + "_" + d["name"]
    data = findDataSet(MCList[i], datasets)
    fileName = rootdir + MCList[i]+'.root'
    #if fileName not in sorted(os.listdir(rootdir)): continue
    f = TFile.Open(fileName,"READ")
    t = f.Get(treepath)
    h = TH2F(histName,histName,binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
    weightedBkgCut = weightCut(bkgCut, weight)
    t.Project(histName, var, weightedBkgCut)
    wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight)
    scale = datalumi*data["xsec"]/wentries
    h.Scale(scale)
    if d["name"]=="Bkg_Up" and "SingleT" in MCList[i]: h.Scale(1.5)
    if d["name"]=="Bkg_Down" and "SingleT" in MCList[i]: h.Scale(0.5)
    print histName, data["xsec"], wentries, scale
    print h.Integral()
    hh = h.Clone()
    hList.append(hh)
    hBkgAll.Add(hh)
  for i,file in enumerate(DYList):
    histName = DYList[i]+ "_" + tree + "_" + d["name"]
    data = findDataSet(DYList[i], datasets)
    fileName = rootdir + DYList[i]+'.root'
    #if fileName not in sorted(os.listdir(rootdir)): continue
    f = TFile.Open(fileName,"READ")
    t = f.Get(treepath)
    wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight)
    scale = datalumi*data["xsec"]/wentries
    hDY = TH2F(histName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])    
    for j,ch in enumerate(ChList):
      dyCut = mAND(overallCut,ChCutList[j])
      weightedDYCut = "(%s)*(%s)" % (dyCut, weight)
      histDYChName = DYList[i]+ "_" + ch + "_" + tree + "_" + d["name"]
      h = TH2F(histDYChName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
      t.Project(histDYChName, var, weightedDYCut)    
      for dych in dyscale:
        if dych["channel"] == ch:dy=dych
      print scale
      scale *= dy["step4"]
      print scale
      h.Scale(scale)
      print h.Integral()
      hh = h.Clone()
      hDY.Add(hh)
    print histName, data["xsec"], wentries, scale
    print hDY.Integral()
    hList.append(hDY)
    hBkgAll.Add(hDY)
  print "bkg tot = ", hBkgAll.Integral()
  if d["name"].split("_")[0] == 'MuF' or d["name"].split("_")[0] == 'MuR':
    scaleW = d["var"].strip(")").split('*')[-1]
    PUweight += "*%s" % scaleW
    print PUweight
  if "pdf" in d["name"]:
    scaleW = d["var"].strip(")").split('*')[-1]
    PUweight += "*%s" % scaleW
    print PUweight

  hTTList = []
  for i,ttName in enumerate(TTList):
    if (ttName!='TT_powheg' and ttName!=d['name']) or (ttName=='TT_powheg' and 'TT' in d['name']):continue 
    data = findDataSet(TTList[i], datasets)
    fileName = rootdir + TTList[i]+'.root'
    wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight) 
    print "test"
    scale = datalumi*data["xsec"]/wentries

    f = TFile.Open(fileName)
    t = f.Get(treepath)
    #nttjjV = getWeightedEntries(fileName, treepath, "filtered", weightCut(visible,weight))
    #nttbbV = getWeightedEntries(fileName, treepath, "filtered", weightCut(ttbb,weight))
    #nttjjF = getWeightedEntries(fileName, treepath, "filtered", weightCut(full,weight))
    #nttbbF = getWeightedEntries(fileName, treepath, "filtered", weightCut(ttbb_full,weight))

    nttjjV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(visible,PUweight),scale)
    nttbbV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttbb,PUweight),scale)
    nttjjF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(full,PUweight),scale)
    nttbbF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttbb_full,PUweight),scale)
    nttjjR = getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(visible,bkgCut),weight),scale)
    nttbbR = getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(ttbb,bkgCut),weight),scale)
    if d["name"]=="Ttcc_Fraction_Up":
      nttjjV += getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc,PUweight),0.5*scale)
      nttjjF += getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc_full,PUweight),0.5*scale) 
      nttjjR += getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(ttcc,bkgCut),weight),0.5*scale)
    elif d["name"]=="Ttcc_Fraction_Down":
      nttjjV -= getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc,PUweight),0.5*scale)
      nttjjF -= getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc_full,PUweight),0.5*scale) 
      nttjjR -= getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(ttcc,bkgCut),weight),0.5*scale)
    out = {'name':d['name'], 'tot':datalumi*data["xsec"], 'nttjjF':nttjjF, 'nttbbF':nttbbF, 'nttjjV':nttjjV, 'nttbbV':nttbbV, 'nttjjR':nttjjR, 'nttbbR':nttbbR}
    output.append(out)
    print out 
    for j,cate in enumerate(ttjj_title):
      histName = TTList[i] + "_" + ttjj_title[j] + "_" + tree + "_" + d["name"]
      ttCut = mAND(bkgCut,ttjj_cut[j])
      ttCutV = ttjj_cut[j]
      ttCutF = ttjj_full_cut[j]
      ttWeightedCut = weightCut(ttCut, weight)
      h = TH2F(histName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
      t.Project(histName, var, ttWeightedCut)
      h.Scale(scale)
      print histName, data["xsec"], wentries, scale
      print histName, " = ", h.Integral()
      hh = h.Clone()
      hTTList.append(hh)
  fCSV = TFile.Open(rootFile, "UPDATE")

  hBkgAll.Write()
  hRDAll.Write()

  for h in hList:
    #print hh.Title()
    h.Write()
  for h in hTTList:
    h.Write()

  fCSV.Close()
  print ""
  print ""
file = open(outputJson, "w")
jsonoutput = json.dumps(output, indent=4)
file.write(jsonoutput)
file.close()
