from ROOT import *
import json, os

info = json.load(open("%s/src/CATTools/CatAnalyzer/test/outputttLepJets.json" % os.environ["CMSSW_BASE"]))
lumi = 36.8*1000

prevSys, prevCh = "",""
totSys = 0.
tempSys = 0.
prevSysName = ""
for i,sys in enumerate(info):
  if sys["sys"]==prevSys and sys["Channel"]!=prevCh:
    nRD += sys["data"]
    nBkg += sys["bkg"]
    nTT += sys["tt"]
    xsec = (nRD-nBkg)/lumi/(nTT/nTTTot)
    if sys["sys"]=="csvweight": 
      xsecCen = xsec
      xsecCenStat = (pow(nRD,0.5))/lumi/(nTT/nTTTot)
    else:
      sysName = sys["sys"].replace("_Up","").replace("_Down","").replace("_NOM","")
      print sysName
      if sysName == prevSysName:
        if abs(xsec-xsecCen)>tempSys: tempSys = abs(xsec-xsecCen)
      else: 
        totSys = pow(totSys*totSys + tempSys*tempSys, 0.5)
        tempSys = abs(xsec-xsecCen)  
        print totSys, tempSys, xsec
        prevSysName = sysName
  elif sys["sys"]!=prevSys: 
    nRD = sys["data"]
    nBkg = sys["bkg"]
    nTT = sys["tt"]
    nTTTot = sys["ttF"]
  prevSys = sys["sys"]
  prevCh = sys["Channel"]

print "xsec = ", xsecCen, " +- ", xsecCenStat, "(Stat.) +- ", totSys, "(Sys.)" 
sysErr = round(totSys/xsecCen*100,3) 
statErr = round(xsecCenStat/xsecCen*100,3)
totErr = round(pow(sysErr*sysErr+statErr*statErr,0.5),3)
print "sys err = ", sysErr, "%" 
print "stat err = ", statErr, "%" 
print "tot err = ", totErr, "%" 
