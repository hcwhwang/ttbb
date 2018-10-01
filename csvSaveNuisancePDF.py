#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from histoHelper import *
from ROOT import *
from sysWeightPDF_cfi import *
#from sysWeightPDF_cfi import *
def unroll(h2D, name):
  n2DXBin = h2D.GetXaxis().GetNbins() 
  n2DYBin = h2D.GetYaxis().GetNbins() 
  nBin = h2D.GetXaxis().GetNbins()*h2D.GetYaxis().GetNbins()
  h = TH1F(name,name,nBin,0,nBin)
  h.Sumw2()
  for i in range(n2DXBin):
    for j in range(n2DYBin):
      content = h2D.GetBinContent(i+1,j+1)
      error = h2D.GetBinError(i+1,j+1)
      binNumber = i*n2DYBin + j + 1
      h.SetBinContent(binNumber,content)
      h.SetBinError(binNumber,error)
  return h

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
overallCut = "tri>0 && filtered==1 && step>=5"
#overallCut = "filtered==1 && step>=5"

channelMCCut = "(channel==1 || channel==2 || channel==3)"
#visible="(NJets20>=6 && NbJets20>=2 && (((lepton1_pt>20 && abs(lepton1_eta)<2.4)) || (lepton2_pt>20 && abs(lepton2_eta)<2.4)))"
full ="(NaddJets20 >= 2)"
ttbb_full = "(NaddbJets20 >= 2)"



visible="(NJets20>=4 && NbJets20>=2 && lepton1_pt>20 && lepton2_pt>20 && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"

ttbb = mAND("(NbJets20>=4)",visible)
ttbj = mAND("(NbJets20==3)",visible)
ttcclf = mAND("!(NbJets20>=3)",visible)
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
ttcclf_full = "( !"+ttbb+" && !"+ttbj+" && NaddJets20 >= 2 )"
ttothers_full = op_(full)

ttjj_cut = [ttbb, ttbj, ttcclf, ttcclf, ttothers]
ttjj_full_cut = [ttbb_full, ttbj_full, ttcclf_full, ttcclf_full, ttothers_full]
ttjj_title = ["ttbb", "ttbj", "ttcclf", "ttcclf2","ttothers"]

name = "Nbjets2VisNestPDF"
#name = "test1"
rootFile = 'csv%s_nuisance.root' % name
outputJson = 'output%s_nuisance.json' % name

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ["CMSSW_BASE"]))
#rootdir = "%s/src/CATTools/CatAnalyzer/test" % os.environ["CMSSW_BASE"]
rootdir = "/xrootd/store/user/chanwook/ttbb/v808/"

datalumi = 36.8*1000
RDList = ['MuonEG_Run2016', 'DoubleEG_Run2016', 'DoubleMuon_Run2016']
MCList = [ 'WJets_part2','SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ_part2', 'WW_part2', 'WZ_part2', 'WWW', 'WZZ', 'ZZZ', 'TTW_part2', 'TTZ_part2',]
DYList = ['DYJets', 'DYJets_10to50_part2']#MCList = [ 'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW','WWZ','WZZ','ZZZ','ttW','ttZ']
TTList = ['TT_powheg', 'FSR_Up', 'FSR_Down', 'ISR_Down', 'ISR_Up', 'UE_Up', 'UE_Down']#, 'TT_aMC', 'TT_powheg_evtgen']
ChList = ['MuEl','ElEl','MuMu']
ChCutList = ['channel==1','channel==2','channel==3']
nomtreepath = "cattree/nom"
output = []
dyscale = json.load(open("drellyanresult.json"))
#binning = [10, 0, 1., 10, 0.5426, 1.]
#binning = [10, 0, 1., 10, 0.8484, 1.]


binning = [10, 0, 1., 10, 0, 1.]
#if os.path.exists(rootFile): os.system("rm -f %s" % rootFile)
for sys,d in enumerate(mceventweight):
  print d["name"],sys
  #if sys < 33: continue
  if d["name"]!='pdf101': continue
  tree = d["tree"]
  treepath = "cattree/"+tree
  weight = d["var"] 
  var = "jets_bDiscriminatorCSV[csvd_jetid[2]]:jets_bDiscriminatorCSV[csvd_jetid[3]]"
  var2 = "jets_bDiscriminatorCSV[csvd_jetid[3]]"

  if d["name"] == "csvweight": datacard = "rate"
  for ch in range(3):
    if d["name"]=="PW_Up": PUweight = "weight*puweightUp"
    elif d["name"]=="PW_Down": PUweight = "weight*puweightDown"
    else: PUweight = "weight*puweight"
    print PUweight
    totCut = mAND(overallCut, ChCutList[ch])
    weightedCut = "(%s)*(%s)" % (totCut, weight)

    if d["name"] == "csvweight":
      hBkgName = ChList[ch] + "_bkg"


    elif "Up" or "Down" in d["name"]:
      sysName = d["name"].replace("_","")
      hBkgName = ChList[ch] + "_bkg" + "_" + sysName

    hBkg = TH2F(hBkgName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
    nBkg = 0
    for i,file in enumerate(MCList):
      #if i!=0: continue
      histName = MCList[i]+ "_" + tree + "_" + d["name"]
      data = findDataSet(MCList[i], datasets)
      fileName = rootdir + MCList[i]+'.root'
      #if fileName not in sorted(os.listdir(rootdir)): continue
      f = TFile.Open(fileName,"READ")
      t = f.Get(treepath)
      h = TH2F(histName,histName,binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
      t.Project(histName, var, weightedCut)
      wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight)
      scale = datalumi*data["xsec"]/wentries
      h.Scale(scale)
      nBkgTemp = getWeightedEntries2D(fileName, treepath, var, weightedCut,scale)
      if d["name"]=="Bkg_Up" and "SingleT" in MCList[i]: h.Scale(1.5)
      if d["name"]=="Bkg_Down" and "SingleT" in MCList[i]: h.Scale(0.5)

      hAndN = mEntryToZero(h,nBkgTemp)      
      print histName, data["xsec"], wentries, scale
      print hAndN[0].Integral(), hAndN[1] 
      hh = hAndN[0].Clone()
      hBkg.Add(hh)
      nBkg += hAndN[1]
    for i,file in enumerate(DYList):
      histName = DYList[i]+ "_" + tree + "_" + d["name"]
      data = findDataSet(DYList[i], datasets)
      fileName = rootdir + DYList[i]+'.root'
      #if fileName not in sorted(os.listdir(rootdir)): continue
      f = TFile.Open(fileName,"READ")
      t = f.Get(treepath)
      wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight)
      scale = datalumi*data["xsec"]/wentries
      dyCut = mAND(overallCut,ChCutList[ch])
      weightedDYCut = "(%s)*(%s)" % (dyCut, weight)
      histDYChName = DYList[i]+ "_" + ChList[ch] + "_" + tree + "_" + d["name"]
      h = TH2F(histDYChName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
      t.Project(histDYChName, var, weightedDYCut)    
      for dych in dyscale:
        if dych["channel"] == ChList[ch]:dy=dych
      print scale
      scale *= dy["step4"]
      print scale
      nBkgTemp = getWeightedEntries2D(fileName, treepath, var, weightedCut,scale)
      h.Scale(scale)

      hAndN = mEntryToZero(h,nBkgTemp)      
      print histName, data["xsec"], wentries, scale
      print hAndN[0].Integral(), hAndN[1] 

      hh = hAndN[0].Clone()
      hBkg.Add(hh)
      nBkg += hAndN[1]
    print "bkg tot = ", hBkg.Integral()

    if d["name"].split("_")[0] == 'MuF' or d["name"].split("_")[0] == 'MuR'or d["name"].split("_")[0] == 'ME':
      scaleW = d["var"].strip(")").split('*')[-1]
      PUweight += "*%s" % scaleW
      print PUweight


    hTTList = []
    for i,ttName in enumerate(TTList):
      #if ttName != 'TT_powheg': continue
      if (ttName!='TT_powheg' and ttName!=d['name']) or (ttName=='TT_powheg' and 'TT' in d['name']):continue 
      if ttName == 'FSR_Down': 
        ttName = 'TT_powheg_fsrdown'
        sysName = 'FSRDown' 
      elif ttName == 'FSR_Up': 
        ttName = 'TT_powheg_fsrup'
        sysName = 'FSRUp' 
      elif ttName == 'ISR_Down': 
        ttName = 'TT_powheg_isrdown'
        sysName = 'ISRDown' 
      elif ttName == 'ISR_Up': 
        ttName = 'TT_powheg_isrup'
        sysName = 'ISRUp' 
      elif ttName == 'UE_Down': 
        ttName = 'TT_powheg_down'
        sysName = 'UEDown' 
      elif ttName == 'UE_Up': 
        ttName = 'TT_powheg_up'
        sysName = 'UEUp' 

 
      data = findDataSet(ttName, datasets)
      fileName = rootdir + ttName+'.root'
      wentries = getWeightedEntries(fileName, nomtreepath, "filtered", PUweight) 
      print "test"
      scale = datalumi*data["xsec"]/wentries

      f = TFile.Open(fileName)
      t = f.Get(treepath)
      #nttjjV = getWeightedEntries(fileName, treepath, "filtered", weightCut(visible,weight))
      #nttbbV = getWeightedEntries(fileName, treepath, "filtered", weightCut(ttbb,weight))
      #nttjjF = getWeightedEntries(fileName, treepath, "filtered", weightCut(full,weight))
      #nttbbF = getWeightedEntries(fileName, treepath, "filtered", weightCut(ttbb_full,weight))

      nttjjF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(full,PUweight),scale)
      nttbbF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttbb_full,PUweight),scale)
      nttjjV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(visible,PUweight),scale)
      nttbbV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttbb,PUweight),scale)
      nttjjFErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(full,PUweight),scale)
      nttbbFErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(ttbb_full,PUweight),scale)
      nttjjVErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(visible,PUweight),scale)
      nttbbVErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(ttbb,PUweight),scale)
      #nttbjV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttbj,PUweight),scale)
      nttcclfV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcclf,PUweight),scale)
      #nttothersV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttothers,PUweight),scale)
      #nttjjF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(mAND(full,ChCutList[ch]),PUweight),scale)
      #nttbbF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(mAND(ttbb_full,ChCutList[ch]),PUweight),scale)
      nttjjR = getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(visible,overallCut),ChCutList[ch]),weight),scale)
      nttbbR = getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttbb,overallCut),ChCutList[ch]),weight),scale)
      nttbjR = getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttbj,overallCut),ChCutList[ch]),weight),scale)
      nttcclfR = getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttcclf,overallCut),ChCutList[ch]),weight),scale)
      #nttothersR = getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttothers,overallCut),ChCutList[ch]),weight),scale)
      '''
      if d["name"]=="Ttcc_Fraction_Up":
        nttjjV += getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc,PUweight),0.5*scale)
        nttjjF += getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc_full,PUweight),0.5*scale) 
        nttjjR += getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttcc,overallCut),ChCutList[ch]),weight),0.5*scale)
        nttcclfR += getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttcc,overallCut),ChCutList[ch]),weight),0.5*scale)
      elif d["name"]=="Ttcc_Fraction_Down":
        nttjjV -= getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc,PUweight),0.5*scale)
        nttjjF -= getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttcc_full,PUweight),0.5*scale) 
        nttjjR -= getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttcc,overallCut),ChCutList[ch]),weight),0.5*scale)
        nttjjR -= getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttcc,overallCut),ChCutList[ch]),weight),0.5*scale)
        nttcclfR -= getWeightedEntries2D(fileName, treepath, var, weightCut(mAND(mAND(ttcc,overallCut),ChCutList[ch]),weight),0.5*scale)
      '''
      out = {'name':d['name'], 'nttjjF':nttjjF, 'nttbbF':nttbbF, 'nttjjV':nttjjV, 'nttbbV':nttbbV, 'nttjjFErr':nttjjFErr, 'nttbbFErr':nttbbFErr, 'nttjjVErr':nttjjVErr, 'nttbbVErr':nttbbVErr}
      if ch==0:output.append(out)
      #print out 
      print "init r_vis_ttbb = ", nttbbV/nttjjV
      print "init r_vis_ttcclf = ", nttcclfV/nttjjV
      print "init xsec_vis_ttjj = ", nttjjV/datalumi
      print "init xsec_vis_ttbb = ", nttbbV/datalumi
      for j,cate in enumerate(ttjj_title):
        if d["name"] == "csvweight":
          hTTName = ChList[ch] + "_" + ttjj_title[j]

        elif "Up" or "Down" in d["name"]:
          hTTName = ChList[ch] + "_" + ttjj_title[j] + "_" + sysName

        ttCut = mAND(overallCut,ttjj_cut[j])
        ttWeightedCut = weightCut(mAND(ttCut,ChCutList[ch]), weight)

        ttjjCut = mAND(overallCut,visible)
        ttjjWeightedCut = weightCut(mAND(ttCut,ChCutList[ch]), weight)

        h = TH2F(hTTName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
        t.Project(hTTName, var, ttWeightedCut)
        '''
        if cate == "ttbb": ttScale = scale*datalumi/nttbbV
        elif cate == "ttbj": ttScale = scale*datalumi/nttbbV
        elif cate == "ttcclf": ttScale = scale*datalumi/nttcclfV
        else: ttScale = scale*datalumi/nttjjV
        '''
        if cate == "ttbb": ttScale = scale*datalumi/nttbbV
        elif cate == "ttbj": ttScale = scale*datalumi/nttbbV
        elif cate == "ttcclf": ttScale = scale*datalumi/nttcclfR*(nttbbR+nttbjR)/nttbbV
        elif cate == "ttcclf2": ttScale = scale*datalumi/nttjjV*nttjjR/nttcclfR
        else: ttScale = scale*datalumi/nttjjV
        print "test integral ", h.Integral()
        #print "test getWeigted ", getWeightedEntries2D(fileName, treepath, var, ttjjWeightedCut)
        h.Scale(ttScale)
        print hTTName, data["xsec"], wentries, ttScale
        hAndN = mEntryToZero(h,h.Integral())      
        print hAndN[0].Integral(), hAndN[1] 

        hh = hAndN[0].Clone()
        hh = h.Clone()
        hTTList.append(hh)
        #nTTJJCh = getWeightedEntries2D(fileName, treepath, var, ttjjWeightedCut,scale)
        #nTTCh = getWeightedEntries2D(fileName, treepath, var, ttWeightedCut,ttScale)
        if cate in ["ttbb","ttbj","ttcclf"]: print "n%s = " % cate, hh.Integral()
        else: print "n%s = " % cate, hh.Integral()
        if d["name"] == "csvweight": datacard += " %f" % hh.Integral()
    hRDName = ChList[ch] + "_data_obs"
    if sys==0:
      print hRDName
      fileName = rootdir + RDList[ch]+'.root'
      f = TFile.Open(fileName)
      t = f.Get(treepath)
      hRD = TH2F(hRDName,"",binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
      t.Project(hRDName, var, totCut)
      print hRD.Integral()
      #if (i==0): hRDAll = h.Clone()
      #else : hRDAll.Add(h)
    print "nBkg = ", nBkg
    if d["name"] == "csvweight": datacard += " %f" % nBkg
    fCSV = TFile.Open(rootFile, "UPDATE")

    #unroll(hBkg,ChList[ch] + "_bkg" + "_" + sysName).Write()
    #unroll(hBkg,ChList[ch] + "_bkg").Write()
    unroll(hBkg,hBkg.GetName()).Write()
    if sys==0: unroll(hRD,hRD.GetName()).Write()

    for h in hTTList:
      unroll(h,h.GetName()).Write()
    fCSV.Close()
    print ""
    print ""
if datacard: print datacard
print output
file = open(outputJson, "w")
jsonoutput = json.dumps(output, indent=4)
file.write(jsonoutput)
file.close()
