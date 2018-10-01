#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from histoHelper import *
from ROOT import *
from sysWeight_cfi import *

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
def mAND(aaa,bbb):
  return "(" +aaa+ " && "+bbb+")"
def op_(aaa):
  return "!(" + aaa + ")"
def weightCut(a, b):
  cut = "(%s)*(%s)" % (a, b)
  return cut

datalumi = 36.8
#CMS_lumi.lumi_sqrtS = "%.2f fb^{-1}, #sqrt{s} = 13 TeV "%(datalumi)

#analysis = "TtbarBbbarDiLeptonAnalyzer_"
analysis = ""
store = "/xrootd/store/user/chanwook/ttbb/final/preapproval/"
#overallCut = "tri>0 && filtered==1 && step>=5 && (channel==1 || channel==2 || channel==3) && nbjetL30>=3"
overallCut = "tri>0 && filtered==1 && step>=5 && (channel==1 || channel==2 || channel==3)"
#overallCut = "filtered==1 && step>=5"
channelMCCut = "(channel==1 || channel==2 || channel==3)"

#visible="(NJets20>=6 && NbJets20>=2 && (((lepton1_pt>20 && abs(lepton1_eta)<2.4)) || (lepton2_pt>20 && abs(lepton2_eta)<2.4)))"
visible="(NJets30>=4 && NbJets30>=2 && ((lepton1_pt>20 && lepton2_pt>25) || (lepton1_pt>25 && lepton2_pt>20)) && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"

ttbb = mAND("(NbJets30>=4)",visible)
#ttb = mAND("(NbJets30==3 && !(genTtbarId%100==52))",visible)
#tt2b = mAND("(NbJets30==3 && (genTtbarId%100==52))",visible)
ttbj = mAND("(NbJets30==3)",visible)
ttcc = mAND("((NcJets30>=2) && !(NbJets30>=3))",visible)
ttlf = mAND("(!(NbJets30>=4) && !(NbJets30==3) && !(NcJets30>=2))",visible)
ttothers = op_(visible)

ttbarid_cut = [ttbb,ttbj,ttcc,ttlf]
ttbarid_title = ["t#bar{t}b#bar{b}", "t#bar{t}bj","t#bar{t}c#bar{c}", "t#bar{t}LF"]
ttbarid_name = ["ttbb","ttbj","ttcc","ttLF"]
ttbarid_linecolor = [2,28,3,4]


#sysList = ["PW","JER","JES","HF","HF_Stats1","HF_Stats2","LF","LF_Stats1","LF_Stats2","CQ_Err1","CQ_Err2","Mu_Eff","El_Eff","Mu_Pt","El_Pt","Trig"]
orderName = ["1st jet","2nd jet","3rd jet","4th jet"]

f = TFile.Open(store + "TT_powheg_tot.root")

binning = [10, 0.5426, 1.]
binning2D = [10, 0.5426, 1, 10, 0.5426, 1]

#binning = [10, 0, 1.]
#binning2D = [10, 0, 1, 10, 0, 1]

binning = [5, 0, 1.]
binning2D = [5, 0, 1, 5, 0, 1]

name = "Nbjets2Bin"
var2D = "jets_bDiscriminatorCSV[csvd_jetid[3]]:jets_bDiscriminatorCSV[csvd_jetid[2]]"
for i in range(4):
  la = TLatex()
  leg = TLegend(0.65, 0.7, 0.92, 0.92)
  leg.SetFillStyle(0)
  leg.SetBorderSize(0)
  var = "jets_bDiscriminatorCSV[csvd_jetid[%i]]"%i
  canName = "csv_%i_%s" % ((i+1),name)
  c = TCanvas(canName,"",500,500)
  hFrame = TH1F(canName,";b jet discriminator;Normalized Entries",binning[0],binning[1],binning[2]) 
  hList = []
  for l,cut in enumerate(ttbarid_cut): 
    dataVal = findDataSet("csvweight",mceventweight)
    tree =  "cattree/"+dataVal["tree"]
    hName = "%s" % (ttbarid_name[l]) 
    h = ROOT.TH1F(hName,";b jet discriminator;Normalized Entries",binning[0],binning[1],binning[2])
    ttcut = mAND(overallCut, ttbarid_cut[l]) 
    cut = "(%s)*(%s)" % (ttcut, dataVal["var"])
    t = f.Get(tree)
    t.Project(hName, var, cut)
    h.SetLineColor(ttbarid_linecolor[l])
    h.SetLineWidth(2)
    hList.append(h)
    leg.AddEntry(h,ttbarid_title[l],"l")
    hFrame.SetMinimum(0)
    hMax = max(h.GetMaximum()/h.Integral() for h in hList)
    hFrame.SetMaximum(hMax*1.4)
    hFrame.Draw()
    
    if i==0:
      la2D= TLatex()
      c2D = TCanvas("csv_2D_%s_%s" % (ttbarid_name[l],name),"",500,500)
      h2D = ROOT.TH2F(hName+"2D",";b discriminator of 3rd jet;b discriminator of 4th jet",binning2D[0],binning2D[1],binning2D[2],binning2D[3],binning2D[4],binning2D[5])
      t.Project(hName+"2D", var2D, cut)
      h2D.DrawNormalized("COLZ")
      la2D.DrawLatex((binning2D[2]-binning2D[1])*0.1+binning2D[1],(binning2D[5]-binning2D[4])*0.8+binning2D[4],ttbarid_title[l])
      CMS_lumi.CMS_lumi(c2D,0,11)
      c2D.SaveAs("%s.png" % c2D.GetTitle())
      c2D.SaveAs("%s.pdf" % c2D.GetTitle())

  c.cd()
  for h in hList:
    h.DrawNormalized("sameHIST")
  leg.Draw()
  CMS_lumi.CMS_lumi(c,0,11)
  la.DrawLatex((binning[2]-binning[1])*0.1+binning[1],hFrame.GetMaximum()*0.8,orderName[i])
  c.SaveAs("%s.png" % c.GetTitle())
  c.SaveAs("%s.pdf" % c.GetTitle())

