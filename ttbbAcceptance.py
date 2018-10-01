from ROOT import *
import json, os

info = json.load(open("%s/src/CATTools/CatAnalyzer/test/outputPA_nuisance.json" % os.environ["CMSSW_BASE"]))
lumi = 35.8*1000

RVis = 0.018
xsecTtjjVis = 2.171
RVisStatErr = 0.0008555
RVisSysErr = 0.001375
xsecTtjjVisStatErr = 0.0186
xsecTtjjVisSysErr = 0.15
rkCor = 0.4
def bOverAErr(b,bErr,a,aErr):
  return  pow(aErr*aErr*b*b/a/a + bErr*bErr, 0.5)/a
def aTimesBErr(a,aErr,b,bErr,corAB=0.):
  return pow(b*b*aErr*aErr+a*a*bErr*bErr +2*a*b*corAB*aErr*bErr,0.5)

xsecTtbbVis = RVis*xsecTtjjVis
xsecTtbbVisStatErr = aTimesBErr(RVis,RVisStatErr,xsecTtjjVis,xsecTtjjVisStatErr,rkCor)
xsecTtbbVisSysErr = aTimesBErr(RVis,RVisSysErr,xsecTtjjVis,xsecTtjjVisSysErr,rkCor)

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

sysList = [["csvweight"],["ISR_Up","ISR_Down"],["FSR_Up","FSR_Down"],["UE_Up","UE_Down"],["MuF_Up","MuF_Down"],["MuR_Up","MuR_Down"],["erdon_Up","erdon_Down"],["qcderdon_Up","qcderdon_Down"],["gluonmove_Up","gluonmove_Down"],["hdamp_Up","hdamp_Down"],["topPt_Up","topPt_Down"],["PDF_Up","PDF_Down"],["PDFAlphaS_Up","PDFAlphaS_Down"]]

AttjjSysErr = 0.
AttbbSysErr = 0.
AccRSysErr = 0.
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

      nttjjR = sys["nttjjR"]
      nttbbR = sys["nttbbR"]

              
      if sysSub == sysName == "csvweight":
        Attjj = nttjjV/nttjjF
        Attbb = nttbbV/nttbbF
        Ettjj = nttjjR/nttjjV
        Ettbb = nttbbR/nttbbV
        
        
        AccR = Attbb/Attjj

        RFull = RVis/AccR

        xsecTtjjFull = xsecTtjjVis/Attjj

        xsecTtbbFull = xsecTtbbVis/Attbb
        #print Attbb
        print "nominal & ", round(Attbb,4), " & ", round(Attjj,4), "\\\\"
        #print "nominal & ", round(Ettbb,3), " & ", round(Ettjj,3), "\\\\"
      if i==0 and abs(Attbb-nttbbV/nttbbF)>0.0001: print sysName, " & ", round((Attbb-nttbbV/nttbbF)/Attbb*100,1), "\% & ", round((Attjj-nttjjV/nttjjF)/Attjj*100,1), "\% \\\\"
      #if i==0 and abs(Ettbb-nttbbR/nttbbV)>0.01: print sysName, " & ", round((Ettbb-nttbbR/nttbbV)/Ettbb*100,1), "\% & ", round((Ettjj-nttjjR/nttjjV)/Ettjj*100,1), "\% \\\\"
      elif sysSub == sysName:
        sysAttjj = nttjjV/nttjjF
        sysAttbb = nttbbV/nttbbF
        sysAccR = sysAttbb/sysAttjj
        
        tempAttjjDiff = abs(Attjj-sysAttjj)
        tempAttbbDiff = abs(Attbb-sysAttbb)
        tempAccRDiff = abs(AccR-sysAccR)

        if tempAttjjDiff > maxAttjjDiff: maxAttjjDiff = tempAttjjDiff
        if tempAttbbDiff > maxAttbbDiff: maxAttbbDiff = tempAttbbDiff
        if tempAccRDiff > maxAccRDiff: maxAccRDiff = tempAccRDiff

  print " "      
  print sysSub, round(abs(maxAttjjDiff/Attjj)*100,2)
  #print sysSub, round(abs(maxAttbbDiff/Attbb)*100,2)
  AttjjSysErr = pow(AttjjSysErr*AttjjSysErr+maxAttjjDiff*maxAttjjDiff,0.5)
  AttbbSysErr = pow(AttbbSysErr*AttbbSysErr+maxAttbbDiff*maxAttbbDiff,0.5)
  AccRSysErr = pow(AccRSysErr*AccRSysErr+maxAccRDiff*maxAccRDiff,0.5)

RFullSysErr = bOverAErr(RVis,RVisSysErr,AccR,AccRSysErr)
xsecTtjjFullSysErr = bOverAErr(xsecTtjjVis,xsecTtjjVisSysErr,Attjj,AttjjSysErr)
xsecTtbbFullSysErr = bOverAErr(xsecTtbbVis,xsecTtbbVisSysErr,Attbb,AttbbSysErr)

print "full ttbb/ttjj = ", RFull
print "stat err", RVisStatErr/(Attbb/Attjj), round(100*(RVisStatErr/(Attbb/Attjj))/RFull,2), "%"
print "sys err", RFullSysErr, round(100*RFullSysErr/RFull,2), "%"
print " "

print "full ttjj xsec = ", xsecTtjjFull
print "stat err", xsecTtjjVisStatErr/Attjj, round(100*xsecTtjjVisStatErr/Attjj/xsecTtjjFull,2), "%"
print "sys err", xsecTtjjFullSysErr, round(100*xsecTtjjFullSysErr/xsecTtjjFull,2), "%"
print " "

print "full ttbb xsec = ", xsecTtbbFull
print "stat err", xsecTtbbVisStatErr/Attbb, round(100*xsecTtbbVisStatErr/Attbb/xsecTtbbFull,2), "%"
print "sys err", xsecTtbbFullSysErr, round(100*xsecTtbbFullSysErr/xsecTtbbFull,2), "%"
print " "


