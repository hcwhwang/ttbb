#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from histoHelper import *
from ROOT import *
from sysWeight_cfi import *
#from sysWeightPDF_cfi import *
def mEntryToZero(h2D, nEvents):
  n2DXBin = h2D.GetXaxis().GetNbins() 
  n2DYBin = h2D.GetYaxis().GetNbins() 
  for i in range(n2DXBin):
    for j in range(n2DYBin):
      content = h2D.GetBinContent(i+1,j+1)
      if content < 0:
        h2D.SetBinContent(i+1,j+1,0)
        h2D.SetBinError(i+1,j+1,0)
        nEvents -= content
  return h2D, nEvents
          
def mAND(aaa,bbb):
  return "(" +aaa+ " && "+bbb+")"
def op_(aaa):
  return "!(" + aaa + ")"
def weightCut(a, b):
  cut = "(%s)*(%s)" % (a, b)
  return cut

#analysis = "TtbarBbbarDiLeptonAnalyzer_"
analysis = ""
#overallCut = "tri>0 && filtered==1 && step>=5 && nbjetL30>=3"
overallCut = "tri>0 && filtered==1 && step>=4"
#overallCut = "filtered==1 && step>=5"



name = "ttLepJets"
rootFile = '%s.root' % name
outputJson = 'output%s.json' % name

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ["CMSSW_BASE"]))
#rootdir = "%s/src/CATTools/CatAnalyzer/test" % os.environ["CMSSW_BASE"]
rootdir = "/xrootd/store/user/heewon/Ttbar/LepJets/"

datalumi = 36.8*1000
RDList = ['SingleElectron_Run2016', 'SingleMuon_Run2016']
MCList = [ 'WJets','WW','WZ','ZZ','DYJets','DYJets_10to50','SingleTop_s','SingleTop_t','SingleTbar_t','SingleTop_tW','SingleTbar_tW','QCD_Pt-15to20_MuEnriched', 'QCD_Pt-20to30_MuEnriched', 'QCD_Pt-30to50_MuEnriched', 'QCD_Pt-50to80_MuEnriched', 'QCD_Pt-80to120_MuEnriched', 'QCD_Pt-120to170_MuEnriched', 'QCD_Pt-170to300_MuEnriched', 'QCD_Pt-300to470_MuEnriched', 'QCD_Pt-470to600_MuEnriched', 'QCD_Pt-600to800_MuEnriched', 'QCD_Pt-800to1000_MuEnriched', 'QCD_Pt-1000toInf_MuEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched', 'QCD_Pt-80to120_EMEnriched', 'QCD_Pt-120to170_EMEnriched', 'QCD_Pt-170to300_EMEnriched', 'QCD_Pt-300toInf_EMEnriched'] 
TTList = ['TT_powheg', 'TT_powheg_FSR_Up', 'TT_powheg_FSR_Down', 'TT_powheg_ISR_Down', 'TT_powheg_ISR_Up', 'TT_powheg_UE_Up', 'TT_powheg_UE_Down', 'TT_powheg_herwig', 'TTJets_aMC']
ChList = ['ElJets','MuJets']
ChCutList = ['channel==1','channel==2']
nomtreepath = "cattree/nom"
output = []
nRDList = []
#binning = [10, 0, 1., 10, 0.5426, 1.]
#binning = [10, 0, 1., 10, 0.8484, 1.]


if os.path.exists(rootFile): os.system("rm -f %s" % rootFile)
nSys = 0
for sys,d in enumerate(mceventweight):
  print d["name"]
  #if d["name"]!='csvweight': continue
  #if "FSR" not in d["name"]: continue
  nSys += 1
  if d["name"]=="PW_Up": PUweight = "weight*puweightUp"
  elif d["name"]=="PW_Down": PUweight = "weight*puweightDown"
  else: PUweight = "weight*puweight"
  tree = d["tree"]
  treepath = "cattree/"+tree
  weight = d["var"] 

  if d["name"] == "csvweight": datacard = "rate"
  for ch in range(2):
    totCut = mAND(overallCut, ChCutList[ch])
    weightedCut = "(%s)*(%s)" % (totCut, weight)

    if d["name"] == "csvweight":
      hBkgName = ChList[ch] + "_bkg"

    elif "Up" or "Down" in d["name"]:
      sysName = d["name"].replace("_","")
      hBkgName = ChList[ch] + "_bkg" + "_" + sysName

    hBkg = TH1D(hBkgName,"",1,0,1)
    nBkg = 0
    for i,file in enumerate(MCList):
      #if i!=0: continue
      histName = MCList[i]+ "_" + tree + "_" + d["name"]
      data = findDataSet(MCList[i], datasets)
      fileName = rootdir + MCList[i]+'.root'
      if MCList[i]+'.root' not in sorted(os.listdir(rootdir)): continue
      h = TH1D(histName,histName,1,0,1)
      wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight)
      scale = datalumi*data["xsec"]/wentries
      nBkgTemp = getWeightedEntries(fileName, treepath, "filtered", weightedCut,scale)
      nBkgTempErr = getWeightedEntriesError(fileName, treepath, "filtered", weightedCut,scale)
      h.SetBinContent(1,nBkgTemp)
      h.SetBinError(1,nBkgTempErr) 
      if d["name"]=="Bkg_Up" and "SingleT" in MCList[i]: h.Scale(1.5)
      if d["name"]=="Bkg_Down" and "SingleT" in MCList[i]: h.Scale(0.5)

      hAndN = mEntryToZero(h,nBkgTemp)      
      print histName, data["xsec"], wentries, scale
      print hAndN[0].Integral(), hAndN[1] 
      hh = hAndN[0].Clone()
      hBkg.Add(hh)
      nBkg += hAndN[1]
    print "bkg tot = ", hBkg.Integral()
    if d["name"].split("_")[0] == 'MuF' or d["name"].split("_")[0] == 'MuR':
      scaleW = d["var"].strip(")").split('*')[-1]
      PUweight += "*%s" % scaleW
      print PUweight
    if "pdf" in d["name"]:
      scaleW = d["var"].strip(")").split('*')[-1]
      PUweight += "*%s" % scaleW
      print PUweight

    hTTList = []
    isTT = False
    for i,ttName in enumerate(TTList):
      if 'TT' not in ttName: continue
      #if ttName != 'TT_powheg': continue
      if (ttName!='TT_powheg' and ttName!=d['name']) or (ttName=='TT_powheg' and 'TT' in d['name']):continue 
      if "FSR" in d['name'] or "ISR" in d['name'] or "UE" in d['name']: 
        ttName = ttName.replace('_Up','up')
        ttName = ttName.replace('_Down','down')
 
      data = findDataSet(ttName, datasets)
      fileName = rootdir + ttName+'.root'
      if ttName+'.root' not in sorted(os.listdir(rootdir)): continue
      isTT = True
      wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight) 
      scale = datalumi*data["xsec"]/wentries

      #nttjjV = getWeightedEntries(fileName, treepath, "filtered", weightCut(visible,weight))
      #nttbbV = getWeightedEntries(fileName, treepath, "filtered", weightCut(ttbb,weight))
      #nttjjF = getWeightedEntries(fileName, treepath, "filtered", weightCut(full,weight))
      #nttbbF = getWeightedEntries(fileName, treepath, "filtered", weightCut(ttbb_full,weight))

      nttF = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight,scale)
      nttR = getWeightedEntries(fileName, treepath, "filtered", weightedCut,scale)
      nttRErr = getWeightedEntriesError(fileName, treepath, "filtered", weightedCut,scale)
      if d["name"] == "csvweight":
        hTTName = ChList[ch] + "_" + "tt"

      elif "Up" or "Down" in d["name"]:
        sysName = d["name"].replace("_","")
        hTTName = ChList[ch] + "_" + "tt" + "_" + sysName

      h = TH1D(hTTName,"",1,0,1)
      h.SetBinContent(1,nttR)
      h.SetBinError(1,nttRErr)
      print hTTName, data["xsec"], nttR
      hAndN = mEntryToZero(h,h.Integral())      
      print hAndN[0].Integral(), hAndN[1] 

      hh = hAndN[0].Clone()
      hh = h.Clone()
      hTTList.append(hh)
      #nTTJJCh = getWeightedEntries2D(fileName, treepath, var, ttjjWeightedCut,scale)
      #nTTCh = getWeightedEntries2D(fileName, treepath, var, ttWeightedCut,ttScale)
      if d["name"] == "csvweight": datacard += " %f" % hh.Integral()
    hRDName = ChList[ch] + "_data_obs"
    if nSys == 1:
    #if sys==0:
      print hRDName
      fileName = rootdir + RDList[ch]+'.root'
      nRD = getWeightedEntries(fileName, nomtreepath, "filtered", totCut)
      nRDErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", totCut)
      hRD = TH1D(hRDName,"",1,0,1)
      hRD.SetBinContent(1,nRD)
      hRD.SetBinError(1,nRDErr)
      nRDList.append(nRD)    
    print "data = ", nRDList[ch]
    print "nBkg = ", nBkg
    if d["name"] == "csvweight": datacard += " %f" % nBkg
    fCSV = TFile.Open(rootFile, "UPDATE")

    if hBkg: hBkg.Write()
    if isTT: hTTList[0].Write()
    if sys==0: hRD.Write()

    fCSV.Close()
    print ""
    print ""
    if isTT: 
      outputSub = {"sys":d["name"],"Channel":ChList[ch],"tt":nttR,"ttErr":nttRErr,"ttF":nttF,"bkg":nBkg,"data":nRDList[ch]}
      print outputSub
      output.append(outputSub)
if datacard: print datacard
file = open(outputJson, "a")
jsonoutput = json.dumps(output, indent=4)
file.write(jsonoutput)
file.close()
