import json, os

def findDataSet(name, datasets):
    for data in datasets:
        if data["ChStep"] == name:
            return data
    return None

inputJson = "entries.json"
title = ["$t\overline{t}+Jets$","$Single Top$","$VV$","$W\\rightarrow l\\nu$","$Z/\gamma^{\star} \\rightarrow ll$","QCD","$Total$","$Data$"]
ChList = ["El","Mu"]


qcdList = ['QCD_Pt-15to20_MuEnriched', 'QCD_Pt-20to30_MuEnriched', 'QCD_Pt-30to50_MuEnriched', 'QCD_Pt-50to80_MuEnriched', 'QCD_Pt-80to120_MuEnriched', 'QCD_Pt-120to170_MuEnriched', 'QCD_Pt-170to300_MuEnriched', 'QCD_Pt-300to470_MuEnriched', 'QCD_Pt-470to600_MuEnriched', 'QCD_Pt-600to800_MuEnriched', 'QCD_Pt-800to1000_MuEnriched', 'QCD_Pt-1000toInf_MuEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched', 'QCD_Pt-80to120_EMEnriched', 'QCD_Pt-120to170_EMEnriched', 'QCD_Pt-170to300_EMEnriched', 'QCD_Pt-300toInf_EMEnriched']
info = json.load(open("%s/src/CATTools/CatAnalyzer/test/%s" % (os.environ["CMSSW_BASE"], inputJson)))
nChList = []
errChList = []
for c in ChList:
  nStepList = []
  errStepList = []
  for s in range(1,6):
    nData, nWJets, nSingleTop, nVV, nDYJets, nQCD = 0,0,0,0,0,0
    nTT = 0
    nDataErr, nWJetsErr, nSingleTopErr, nVVErr, nDYJetsErr, nQCDErr = 0,0,0,0,0,0
    nTT = 0
    nTot, nTotErr = 0,0

    name = "%s%i" % (c,s)
    data = findDataSet(name, info)
    for i,n in enumerate(data["name"]):
      if n in ["SingleElectron","SingleMuon"]:
        nData += data["entries"][i]    
        print nData
        nDataErr += data["error"][i]**2
      if n == "WJets":
        nWJets += data["entries"][i]    
        nWJetsErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2
      if "SingleT" in n:
        nSingleTop += data["entries"][i]    
        nSingleTopErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2
      if n in ["ZZ","WW","ZZ"]:
        nVV += data["entries"][i]    
        nVVErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2
      if "DY" in n:
        nDYJets += data["entries"][i]    
        nDYJetsErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n == "TT_Powheg":
        nTT += data["entries"][i]
        nTTErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n in qcdList:
        nQCD += data["entries"][i]
        nQCDErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
    
    nList = [nTT,nSingleTop, nVV, nWJets, nDYJets, nQCD,nTot, nData]
    errList = [nTTErr,nSingleTopErr, nVVErr, nWJetsErr, nDYJetsErr, nQCDErr,  nTotErr, nDataErr]
    nStrList = []
    nStepList.append(nList)        
    errStepList.append(errList)
  nChList.append(nStepList)
  errChList.append(errStepList)

nAllList = []
errAllList = []
for s in range(5):
  nAllStepList = [0,0,0,0,0,0,0,0]
  errAllStepList = [0,0,0,0,0,0,0,0]
  for c in range(3):
    for i,v in enumerate(nChList[c][s]):  
      nAllStepList[i] += v
      errAllStepList[i] += errChList[c][s][i]
  nAllList.append(nAllStepList) 
  errAllList.append(errAllStepList) 
nChList.append(nAllList)
errChList.append(errAllList)

for i in range(4): 
  for j,t in enumerate(title):
    sList = [0,0,0,0]
    sErrList = [0,0,0,0]
    for k,v in enumerate(nChList[i]):
      sList[k] = round(v[j],1)
      sErrList[k] = round(errChList[i][k][j]**0.5,1)
    if t=="$Data$":
      print t,"&",sList[0],"&",sList[1],"&",sList[2],"&",sList[3],"&",sList[4],"\\\\\hline"
    else:
      print t,"&","$%s \pm %s$" % (sList[0],sErrList[0]),"&","$%s \pm %s$" % (sList[1],sErrList[1]),"&","$%s \pm %s$" % (sList[2],sErrList[2]),"&","$%s \pm %s$" % (sList[3],sErrList[3]),"&","$%s \pm %s$" % (sList[4],sErrList[4]),"\\\\\hline"
  print ""
  print ""
