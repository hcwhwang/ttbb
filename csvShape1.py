#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from histoHelper import *
from ROOT import *
from sysWeight_cfi import *

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
overallCut = "tri>0 && filtered==1 && step>=5 && (channel==1 || channel==2 || channel==3)"
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


sysList = ["PW","JER","JES","HF","HF_Stats1","HF_Stats2","LF","LF_Stats1","LF_Stats2","CQ_Err1","CQ_Err2","Mu_Eff","El_Eff","Mu_Pt","El_Pt","Trig"]
sysList = ["HF","LF"]
sysIndex = ["","_Up","_Down"]
orderName = ["1st jet","2nd jet","3rd jet","4th jet"]

f = TFile.Open(store + "TT_powheg.root")
la = TLatex()

for i,sys in enumerate(sysList):  
  #if i!=0:continue
  print sys
  leg = TLegend(0.65, 0.7, 0.92, 0.92)
  leg.SetFillStyle(0)
  leg.SetBorderSize(0)
  var1 = "jets_bDiscriminatorCSV[csvd_jetid[2]]"
  var2 = "jets_bDiscriminatorCSV[csvd_jetid[3]]"
  var = "jets_bDiscriminatorCSV[csvd_jetid[3]]/jets_bDiscriminatorCSV[csvd_jetid[2]]"
  #canName = "csv_%s_test_inverse" % (sys)
  canName = "csv_%s_test_ratio" % (sys)
  c = TCanvas(canName,"",500,500)
  hFrame = TH1F(canName,";b jet discriminator;Normalized Entries",10,0,1) 
  hList = []
  for k in range(3):
    for l,cut in enumerate(ttbarid_cut): 
      if k==0:data = findDataSet("csvweight",mceventweight)
      else: data = findDataSet(sys+sysIndex[k],mceventweight)
      hName1 = "3rd_%s_%s%s" % (ttbarid_name[l],sys,sysIndex[k]) 
      hName2 = "4th_%s_%s%s" % (ttbarid_name[l],sys,sysIndex[k]) 
      h1 = TH1F(hName1,";b jet discriminator;Normalized Entries",10,0,1)
      h2 = TH1F(hName2,";b jet discriminator;Normalized Entries",10,0,1)
      ttcut = mAND(overallCut, ttbarid_cut[l]) 
      cut = "(%s)*(%s)" % (ttcut, data["var"])
      tree =  "cattree/"+data["tree"]
      t = f.Get(tree)
      #t.Project(hName1, var, cut)
      t.Project(hName1, var1, cut)
      t.Project(hName2, var2, cut)
      h1.SetLineColor(ttbarid_linecolor[l])
      h1.SetLineStyle(k+1)
      h1.SetLineWidth(2)
      for b in range(h1.GetXaxis().GetNbins()):
        #h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h2.GetBinContent(h1.GetXaxis().GetNbins()-b))
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h2.GetBinContent(b+1))
      hList.append(h1)
      if k==0:leg.AddEntry(h1,ttbarid_title[l],"l")
  hFrame.SetMinimum(0)
  hMax = max(h.GetMaximum()/h.Integral() for h in hList)
  hFrame.SetMaximum(hMax*1.4)
  hFrame.Draw()
  for h in hList:
    h.DrawNormalized("sameHIST")
  leg.Draw()
  CMS_lumi.CMS_lumi(c,0,11)
  la.DrawLatex(.1,hFrame.GetMaximum()*0.8,sys)
  c.SaveAs("%s.png"%c.GetTitle())
