import json, os

def findDataSet(name, datasets):
    for data in datasets:
        if data["ChStep"] == name:
            return data
    return None

title = ["$t\\bar{t}+b\\bar{b}$","$t\\bar{t}+bj$","$t\\bar{t}+c\\bar{c}$","$t\\bar{t}$LF","$t\\bar{t}$ others","$t\\bar{t}$+V","Single Top","VV","VVV","$W\\rightarrow l\\nu$","$Z+jets$", "$t\\bar{t}+H \\rightarrow t\\bar{t}+b\\bar{b}$","Total","Data"]
ChList = ["MuEl","ElEl","MuMu"]

info = json.load(open("%s/src/CATTools/CatAnalyzer/test/entriesPA.json" % os.environ["CMSSW_BASE"]))
nChList = []
errChList = []
for c in ChList:
  nStepList = []
  errStepList = []
  for s in range(1,6):
    nData, nWJets, nSingleTop, nVV, nVVV, nttV, nDYJets, nttH = 0,0,0,0,0,0,0,0
    nttbb, nttbj, nttcc, nttlf, nttot = 0,0,0,0,0
    nDataErr, nWJetsErr, nSingleTopErr, nVVErr, nVVVErr, nttVErr, nDYJetsErr, nttHErr = 0,0,0,0,0,0,0,0
    nttbbErr, nttbjErr, nttccErr, nttlfErr, nttotErr = 0,0,0,0,0
    nTot, nTotErr = 0,0

    name = "%s%i" % (c,s)
    data = findDataSet(name, info)
    for i,n in enumerate(data["name"]):
      if n in ["MuonEG","DoubleEG","DoubleMuon"]:
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
      if n in ["ZZ","WW","WZ"]:
        nVV += data["entries"][i]    
        nVVErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2
      if n in ["ZZZ","WZZ","WWZ","WWW"]:
        nVVV += data["entries"][i]    
        nVVVErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2
      
      if n in ["TTW","TTZ"]:
        nttV += data["entries"][i]    
        nttVErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2
      if "DY" in n:
        nDYJets += data["entries"][i]    
        nDYJetsErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n == "ttH_bb":
        nttH += data["entries"][i]    
        nttHErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n == "ttbb":
        nttbb += data["entries"][i]
        nttbbErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n == "ttbj":
        nttbj += data["entries"][i]
        nttbjErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n == "ttcc":
        nttcc += data["entries"][i]
        nttccErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n == "ttlf":
        nttlf += data["entries"][i]
        nttlfErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
      if n == "ttothers":
        nttot += data["entries"][i]
        nttotErr += data["error"][i]**2
        nTot += data["entries"][i]
        nTotErr += data["error"][i]**2  
    nList = [nttbb, nttbj, nttcc, nttlf, nttot, nttV, nSingleTop, nVV, nVVV, nWJets, nDYJets, nttH, nTot, nData]
    errList = [nttbbErr, nttbjErr, nttccErr, nttlfErr, nttotErr, nttVErr, nSingleTopErr, nVVErr, nVVVErr, nWJetsErr, nDYJetsErr, nttHErr, nTotErr, nDataErr]
    nStrList = []
    nStepList.append(nList)        
    errStepList.append(errList)
  nChList.append(nStepList)
  errChList.append(errStepList)

nAllList = []
errAllList = []
for s in range(5):
  nAllStepList = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  errAllStepList = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  for c in range(3):
    for i,v in enumerate(nChList[c][s]):  
      nAllStepList[i] += v
      errAllStepList[i] += errChList[c][s][i]
  nAllList.append(nAllStepList) 
  errAllList.append(errAllStepList) 
nChList.append(nAllList)
errChList.append(errAllList)

def Round(a,b):
  if a < 0 or a==0.0:
    a,b = 0, 0
  elif a >= 100:
    a = int(round(a))
    b = int(round(b))
  elif a >=10:
    a = round(a,1)
    b = round(b,1)
  else:
    a = round(a,2)
    b = round(b,2)
  text = "%s $\pm$ %s" % (a,b)
  return text    

for i in range(4): 
  for j,t in enumerate(title):
    sList = [0,0,0,0,0]
    sErrList = [0,0,0,0,0]
    for k,v in enumerate(nChList[i]):
      sList[k] = round(v[j],1)
      sErrList[k] = round(errChList[i][k][j]**0.5,1)
    if t=="Data":
      print t,"&",int(sList[0]),"&",int(sList[1]),"&",int(sList[2]),"&",int(sList[3]),"&",int(sList[4]),"\\\\"
    else:
      print t,"&",Round(sList[0],sErrList[0]),"&",Round(sList[1],sErrList[1]),"&",Round(sList[2],sErrList[2]),"&",Round(sList[3],sErrList[3]),"&",Round(sList[4],sErrList[4]),"\\\\"
    if j in [0,9,10,11]:
      print "\hline" 
  print ""
  print ""

for i,t in enumerate(title):
  if t=="Data":
    print t,"&",int(nChList[0][4][i]),"&",int(nChList[1][4][i]),"&",int(nChList[2][4][i]),"&",int(nChList[3][4][i]),"\\\\"
  else:
    print t,"&",Round(nChList[0][4][i],errChList[0][4][i]**0.5),"&",Round(nChList[1][4][i],errChList[1][4][i]**0.5),"&",Round(nChList[2][4][i],errChList[2][4][i]**0.5),"&",Round(nChList[3][4][i],errChList[3][4][i]**0.5),"\\\\"
  if i in [0,9,10,11]:
    print "\hline" 
