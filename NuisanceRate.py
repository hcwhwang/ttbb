from ROOT import *
import os

rootFile = "%s/src/CATTools/CatAnalyzer/test/csvNbjets2VisNest_nuisance.root" % os.environ["CMSSW_BASE"]

sysList = ["FSR" ,"ISR", "UE", "MuR", "MuF", "PDF"]
chList = ["MuEl","ElEl","MuMu"]
nameList = [ "ttbb", "ttbj","ttcclf", "ttcclf2", "ttothers", "dy", "vv", "vvv", "ttv", "st"]
f = TFile.Open(rootFile)
print "ttbb ttbj ttcclf ttcclf2 ttothers dy vv vvv ttv st"
for i,sys in enumerate(sysList):
  datacard = "%s rate " % sys
  for j,name in enumerate(nameList):
    for k,ch in enumerate(chList):
      hCen = f.Get("%s_%s" % (ch, name))
      hUp = f.Get("%s_%s_%sUp" % (ch, name, sys))
      hDown = f.Get("%s_%s_%sDown" % (ch, name, sys))
      datacard += (" %f/%f" % (hCen.Integral()/hDown.Integral(), hCen.Integral()/hUp.Integral()))
  print datacard
f.Close()
