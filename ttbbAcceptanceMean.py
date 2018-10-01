from ROOT import *
import json, os

info = json.load(open("%s/src/CATTools/CatAnalyzer/test/outputNbjets2VisNest_nuisance.json" % os.environ["CMSSW_BASE"]))
lumi = 36.8*1000

RVis = 0.018012
xsecTtjjVis = 4.1968
RVisStatErr = 0.0007842
RVisSysErr = 0.001119
xsecTtjjVisStatErr = 0.03397
xsecTtjjVisSysErr = 0.30165
def bOverAErr(b,bErr,a,aErr):
  return  pow(aErr*aErr*b*b/a/a + bErr*bErr, 0.5)/a
def aTimesBErr(a,aErr,b,bErr):
  return pow(b*b*aErr*aErr+a*a*bErr*bErr,0.5)

xsecTtbbVis = RVis*xsecTtjjVis
xsecTtbbVisStatErr = aTimesBErr(RVis,RVisStatErr,xsecTtjjVis,xsecTtjjVisStatErr)
xsecTtbbVisSysErr = aTimesBErr(RVis,RVisSysErr,xsecTtjjVis,xsecTtjjVisSysErr)

print "visible ttbb/ttjj = ", RVis
print "stat err", RVisStatErr, round(100*RVisStatErr/RVis,2), "%"
print "sys err", RVisSysErr, round(100*RVisSysErr/RVis,2), "%"
print " "

print "visible ttjj xsec = ", xsecTtjjVis
print "stat err", xsecTtjjVisStatErr, round(100*xsecTtjjVisStatErr/xsecTtjjVis,2), "%"
print "sys err", xsecTtjjVisSysErr, round(100*xsecTtjjVisSysErr/xsecTtjjVis,2), "%"
print " "

print "visible ttbb xsec = ", xsecTtbbVis
print "stat err", xsecTtbbVisStatErr, round(100*xsecTtbbVisStatErr/xsecTtbbVis,2), "%"
print "sys err", xsecTtbbVisSysErr, round(100*xsecTtbbVisSysErr/xsecTtbbVis,2), "%"
print " "

sysList = [["csvweight"],["PU_Up","PU_Down"],["JER_Up","JER_Down"],["LF_Up","LF_Down"],["HF_Up","HF_Down"],["HF_Stats1_Up","HF_Stats1_Down"],["HF_Stats2_Up","HF_Stats2_Down"],["LF_Stats1_Up","LF_Stats1_Down"],["LF_Stats2_Up","LF_Stats2_Down"],["CQ_Err1_Up","CQ_Err1_Down"],["CQ_Err2_Up","CQ_Err2_Down"],["Mu_Eff_Up","Mu_Eff_Down"],["El_Eff_Up","El_Eff_Down"],["Mu_Scale_Up","Mu_Scale_Down"],["El_Scale_Up","El_Scale_Down"],["Trigger_Up","Trigger_Down"],["Bkg_Up","Bkg_Down"],["ISR_Up","ISR_Down"],["FSR_Up","FSR_Down"],["UE_Up","UE_Down"],["MuF_Up","MuR_Up","MuF_MuR_Up","MuF_Down","MuR_Down","MuF_MuR_Down"]]

AttjjSysErr = 0.
AttbbSysErr = 0.
AccRSysErr = 0.
Attbb, Attjj = 0, 0
totN = 0
for i,sysGroup in enumerate(sysList):
  for j,sysSub in enumerate(sysGroup):
    for k,sys in enumerate(info):
      sysName = sys["name"]

      nttjjF = sys["nttjjF"]
      nttbbF = sys["nttbbF"]
      
      nttjjV = sys["nttjjV"]
      nttbbV = sys["nttbbV"]

              
      if sysSub == sysName:
        tempAttjj = nttjjV/nttjjF
        tempAttbb = nttbbV/nttbbF
        Attjj += tempAttjj
        Attbb += tempAttbb
        totN += 1
Attjj /= totN
Attbb /= totN
        
AccR = Attbb/Attjj

RFull = RVis/AccR
RFullStatErr = RVisStatErr/AccR

xsecTtjjFull = xsecTtjjVis/Attjj
xsecTtjjFullStatErr = xsecTtjjVisStatErr/Attjj

xsecTtbbFull = xsecTtbbVis/Attbb
print Attbb
xsecTtbbFullStatErr = xsecTtbbVisStatErr/Attbb

for i,sysGroup in enumerate(sysList):
  maxAttjjDiff = 0.
  maxAttbbDiff = 0.
  maxAccRDiff = 0.
  for j,sysSub in enumerate(sysGroup):
    for k,sys in enumerate(info):
      sysName = sys["name"]

      nttjjF = sys["nttjjF"]
      nttbbF = sys["nttbbF"]
      
      nttjjV = sys["nttjjV"]
      nttbbV = sys["nttbbV"]

              
      if sysSub == sysName != "csvweight":
        sysAttjj = nttjjV/nttjjF
        sysAttbb = nttbbV/nttbbF
        sysAccR = sysAttbb/sysAttjj
        
        #print sysAccR
        tempAttjjDiff = abs(Attjj-sysAttjj)
        tempAttbbDiff = abs(Attbb-sysAttbb)
        tempAccRDiff = abs(AccR-sysAccR)

        if tempAttjjDiff > maxAttjjDiff: maxAttjjDiff = tempAttjjDiff
        if tempAttbbDiff > maxAttbbDiff: maxAttbbDiff = tempAttbbDiff
        if tempAccRDiff > maxAccRDiff: maxAccRDiff = tempAccRDiff
  print sysSub, round(abs(maxAttjjDiff/Attjj)*100,2)
  #print sysSub, round(abs(maxAttbbDiff/Attbb)*100,2)
  AttjjSysErr = pow(AttjjSysErr*AttjjSysErr+maxAttjjDiff*maxAttjjDiff,0.5)
  AttbbSysErr = pow(AttbbSysErr*AttbbSysErr+maxAttbbDiff*maxAttbbDiff,0.5)
  AccRSysErr = pow(AccRSysErr*AccRSysErr+maxAccRDiff*maxAccRDiff,0.5)

RFullSysErr = bOverAErr(RVis,RVisSysErr,AccR,AccRSysErr)
xsecTtjjFullSysErr = bOverAErr(xsecTtjjVis,xsecTtjjVisSysErr,Attjj,AttjjSysErr)
xsecTtbbFullSysErr = bOverAErr(xsecTtbbVis,xsecTtbbVisSysErr,Attbb,AttbbSysErr)

print "full ttbb/ttjj = ", RFull
print "stat err", RFullStatErr, round(100*RFullStatErr/RFull,2), "%"
print "sys err", RFullSysErr, round(100*RFullSysErr/RFull,2), "%"
print " "

print "full ttjj xsec = ", xsecTtjjFull
print "stat err", xsecTtjjFullStatErr, round(100*xsecTtjjFullStatErr/xsecTtjjFull,2), "%"
print "sys err", xsecTtjjFullSysErr, round(100*xsecTtjjFullSysErr/xsecTtjjFull,2), "%"
print " "

print "full ttbb xsec = ", xsecTtbbFull
print "stat err", xsecTtbbFullStatErr, round(100*xsecTtbbFullStatErr/xsecTtbbFull,2), "%"
print "sys err", xsecTtbbFullSysErr, round(100*xsecTtbbFullSysErr/xsecTtbbFull,2), "%"
print " "


