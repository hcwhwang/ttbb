#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from histoHelper import *
from ROOT import *
from sysWeight_cfi import *
from math import *

gROOT.SetBatch(True)
def mAND(aaa,bbb):
  return "(" +aaa+ " && "+bbb+")"
def op_(aaa):
  return "!(" + aaa + ")"
def weightCut(a, b):
  cut = "(%s)*(%s)" % (a, b)
  return cut

datalumi = 36.5
CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV "%(datalumi)

#analysis = "TtbarBbbarDiLeptonAnalyzer_"
analysis = ""
store = "/xrootd/store/user/chanwook/ttbb/"
overallCut = "tri>0 && filtered==1 && step>=5 && (channel==1 || channel==2 || channel==3) && nbjetL30>=3"
#overallCut = "filtered==1 && step>=5"
channelMCCut = "(channel==1 || channel==2 || channel==3)"

#visible="(NJets20>=6 && NbJets20>=2 && (((lepton1_pt>20 && abs(lepton1_eta)<2.4)) || (lepton2_pt>20 && abs(lepton2_eta)<2.4)))"
visible="(NJets20>=4 && NbJets20>=2 && lepton1_pt>20 && lepton2_pt>20 && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"

ttbb = mAND("(NbJets20>=4)",visible)
#ttb = mAND("(NbJets20==3 && !(genTtbarId%100==52))",visible)
#tt2b = mAND("(NbJets20==3 && (genTtbarId%100==52))",visible)
ttbj = mAND("(NbJets20==3)",visible)
ttcc = mAND("((NcJets20>=2) && !(NbJets20>=3))",visible)
ttlf = mAND("(!(NbJets20>=4) && !(NbJets20==3) && !(NcJets20>=2))",visible)
ttothers = op_(visible)

ttbarid_cut = [ttbb,ttbj,ttcc,ttlf]
ttbarid_title = ["t#bar{t}b#bar{b}", "t#bar{t}bj","t#bar{t}c#bar{c}", "t#bar{t}lf"]
ttbarid_name = ["ttbb","ttbj","ttcc","ttlf"]
ttbarid_linecolor = [2,28,3,4]


#sysList = ["PW","JER","JES","HF","HF_Stats1","HF_Stats2","LF","LF_Stats1","LF_Stats2","CQ_Err1","CQ_Err2","Mu_Eff","El_Eff","Mu_Pt","El_Pt","Trig"]
sysList = ["HF","HF_Stats1","HF_Stats2","LF","LF_Stats1","LF_Stats2","CQ_Err1","CQ_Err2"]
orderName = ["1st jet","2nd jet","3rd jet","4th jet"]

f = TFile.Open(store + "TT_powheg.root")
la = TLatex()

sysIndex = ["","_Up","_Down"]
binning = [10, 0.5426, 1., 10, 0., 1.]
#binning = [10, 0., 1., 10, 0., 1.]
#binning = [10, 0, 1.]
var = "jets_bDiscriminatorCSV[csvd_jetid[3]]:jets_bDiscriminatorCSV[csvd_jetid[2]]"
for i,sys in enumerate(sysList):  
  #if i!=0:continue
  print sys
  #if j!=2:continue
  for l,cut in enumerate(ttbarid_cut): 
    print ttbarid_name[l]
    hList = []
    for k in range(3):
      if k==0: dataVal = findDataSet("csvweight",mceventweight)
      else : dataVal = findDataSet(sys+sysIndex[k],mceventweight)
      tree =  "cattree/"+dataVal["tree"]
      hName = "%s_%s%s" % (ttbarid_name[l],sys,sysIndex[k]) 
      h = ROOT.TH2F(hName,"",binning[0],binning[1],binning[2],binning[3],binning[4],binning[5])
      ttcut = mAND(overallCut, ttbarid_cut[l]) 
      cut = "(%s)*(%s)" % (ttcut, dataVal["var"])
      t = f.Get(tree)
      t.Project(hName, var, cut)
      h.Scale(1./h.Integral())
      hList.append(h)
    Lup, Ldown = 0., 0.
    ndof = 0.
    for x in range(binning[0]):
      for y in range(binning[3]):
        nCenBin = hList[0].GetBinContent(x+1,y+1) 
        nUpBin = hList[1].GetBinContent(x+1,y+1) 
        nDownBin = hList[2].GetBinContent(x+1,y+1)
        if nCenBin == 0: continue
        ndof += 1
        Lup += (nCenBin-nUpBin)**2/nCenBin  
        Ldown += (nCenBin-nDownBin)**2/nCenBin 
    print "Chi^2/ndof = ", Lup/ndof, Ldown/ndof
    print " "
  print " " 
