from ROOT import *
import json, os

def bOverAErr(b,bErr,a,aErr):
  return  pow(aErr*aErr*b*b/a/a + bErr*bErr, 0.5)/a
def findDataSet(name, datasets):
    for data in datasets:
        if data["name"] == name:
            return data
    return None


info = json.load(open("%s/src/CATTools/CatAnalyzer/test/outputNbjets2VisNestTot_nuisance.json" % os.environ["CMSSW_BASE"]))

lumi = 35.8*1000
scale = ["MuF", "MuR"]
for i,sub in enumerate(info):
  if sub["name"]=='csvweight':
    visR = sub["nttbbV"]/sub["nttjjV"]
    visTTJJ = sub["nttjjV"]/lumi
    visTTBB = sub["nttbbV"]/lumi
    fullR = sub["nttbbF"]/sub["nttjjF"]
    fullTTJJ = sub["nttjjF"]/lumi
    fullTTBB = sub["nttbbF"]/lumi
visRErr, visTTJJErr, visTTBBErr = 0, 0, 0
fullRErr, fullTTJJErr, fullTTBBErr = 0, 0, 0
for sys in scale:
  dataUp = findDataSet(sys+"_Up", info)
  dataDown = findDataSet(sys+"_Down", info)

  tempVisRErr = max(abs(visR-dataUp["nttbbV"]/dataUp["nttjjV"]), abs(visR-dataDown["nttbbV"]/dataDown["nttjjV"]))
  visRErr += tempVisRErr*tempVisRErr
  
  tempVisTTJJErr = max(abs(visTTJJ-dataUp["nttjjV"]/lumi), abs(visTTJJ-dataDown["nttjjV"]/lumi))
  visTTJJErr += tempVisTTJJErr*tempVisTTJJErr

  tempVisTTBBErr = max(abs(visTTBB-dataUp["nttbbV"]/lumi), abs(visTTBB-dataDown["nttbbV"]/lumi))
  visTTBBErr += tempVisTTBBErr*tempVisTTBBErr

  tempFisRErr = max(abs(fullR-dataUp["nttbbF"]/dataUp["nttjjF"]), abs(fullR-dataDown["nttbbF"]/dataDown["nttjjF"]))
  fullRErr += tempFisRErr*tempFisRErr
  
  tempFisTTJJErr = max(abs(fullTTJJ-dataUp["nttjjF"]/lumi), abs(fullTTJJ-dataDown["nttjjF"]/lumi))
  fullTTJJErr += tempFisTTJJErr*tempFisTTJJErr

  tempFisTTBBErr = max(abs(fullTTBB-dataUp["nttbbF"]/lumi), abs(fullTTBB-dataDown["nttbbF"]/lumi))
  fullTTBBErr += tempFisTTBBErr*tempFisTTBBErr

print visR, pow(visRErr,0.5)
print visTTJJ, pow(visTTJJErr,0.5)
print visTTBB, pow(visTTBBErr,0.5)
print fullR, pow(fullRErr,0.5)
print fullTTJJ, pow(fullTTJJErr,0.5)
print fullTTBB, pow(fullTTBBErr,0.5)
