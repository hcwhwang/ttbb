from ROOT import *
import json, os

  

#info = json.load(open("%s/src/HiggsAnalysis/CombinedLimit/test/impactsJESBTag.json" % os.environ["CMSSW_BASE"]))
#info = json.load(open("%s/src/HiggsAnalysis/CombinedLimit/test/impactsRate.json" % os.environ["CMSSW_BASE"]))
info = json.load(open("%s/src/HiggsAnalysis/CombinedLimit/test/impacts.json" % os.environ["CMSSW_BASE"]))
accInfo = json.load(open("/cms/ldap_home/chanwook/work/CatProduction/cattools/src/CATTools/CatAnalyzer/test/outputNbjets2VisNestTot_nuisance.json"))
rootFileName1 = 'higgsCombine_paramFit_Test_'
rootFileName2 = '.MultiDimFit.mH120.root'
nameList = [['BtagLF','b-tag light flavour'],['BtagHF','b-tag heavy flavour'],['BtagLFStats1','b-tag light flavour stat. 1'],['BtagHFStats1','b-tag heavy flavour stat. 1'],['BtagLFStats2','b-tag light flavour stat. 2'],['BtagHFStats2','b-tag heavy flavour stat. 2'],['BtagCQErr1','b-tag c flavour stat. 1'],['BtagCQErr2','b-tag c flavour stat. 2'],['PU','PileUp'],['JER','Jet Energy Resolution'],['Trigger','Trigger Efficiency'],['MuEff','Muon Efficiency'],['ElEff','Electron Efficiency'],['MuScale','Muon Energy Scale'],['ElScale','Electron Energy Scale'],['FSR','FSR'],['ISR','ISR'],['PDF','PDF'],['PDFAlphaS','PDF - ${\\alpha}_{s}$'],['erdon', 'Color Reconnection - erdon'],['qcderdon', 'Color Reconnection - qcderdon'],['gluonmove', 'Color Reconnection - gluonmove'],['hdamp','ME-PS Matching'],['topPt', 'Top ${p}_{T}$ Reweighting'],['UE','Underlying Event'],['MuR','${Q}^{2}$ Scale - Renormalization'],['MuF','${Q}^{2}$ Scale - Factorization'],['Lumi','Luminosity'],['st','Background Modeling - Single Top'],['vv','Background Modeling - VV'],['vvv','Background Modeling - VVV'],['ttv','Background Modeling - tt+V'],['dy','Background Modeling - Drell-Yan']]
def findInfo(sys):
  for i,info in enumerate(accInfo):
    if info["name"]==sys: return info
totImpactR = 0.
totImpactK = 0.
totImpactRUp = 0.
totImpactRDown = 0.
bestR = 0.
bestK = 0.

sysList = ""
mcstatRErr = 0.
btagRErr = 0.
jesRErr = 0.
mcstatKErr = 0.
btagKErr = 0.
jesKErr = 0.

for i,subInfo in enumerate(info["POIs"]):
  if subInfo["name"]=="r": bestR = subInfo["fit"][1]
  elif subInfo["name"]=="k": bestK = subInfo["fit"][1]


rVis = bestR
ttjjVis = bestK
ttbbVis = bestR*bestK
statVisR = bestR*4.6/100
statVisTtjj = bestK*0.82/100
statVisTtbb = pow(rVis*rVis*statVisTtjj*statVisTtjj+ttjjVis*ttjjVis*statVisR*statVisR,0.5)

cenInfo = findInfo("csvweight")
cenAttbb = cenInfo["nttbbV"]/cenInfo["nttbbF"]
cenAttjj = cenInfo["nttjjV"]/cenInfo["nttjjF"]

rFull = rVis/cenAttbb*cenAttjj
ttjjFull = ttjjVis/cenAttjj
ttbbFull = ttbbVis/cenAttbb

statFullR = statVisR/cenAttbb*cenAttjj
statFullTtjj = statVisTtjj/cenAttjj
statFullTtbb = statVisTtbb/cenAttbb

impactArray = []
tempSysFullR = 0.
tempSysFullTtjj = 0.
tempSysFullTtbb = 0.
for i,subInfo in enumerate(info["params"]):
  if findInfo(subInfo["name"]+"_Up"):
    upInfo = findInfo(subInfo["name"]+"_Up")
    upAccTtbb = upInfo["nttbbV"]/upInfo["nttbbF"]
    upAccTtjj = upInfo["nttjjV"]/upInfo["nttjjF"]

    upRVis = subInfo["r"][2]
    upTtjjVis = subInfo["k"][2]
    upTtbbVis = upRVis*upTtjjVis

    upRFull = upRVis/upAccTtbb*upAccTtjj
    upTtjjFull = upTtjjVis/upAccTtjj
    upTtbbFull = upTtbbVis/upAccTtbb

    downInfo = findInfo(subInfo["name"]+"_Down")
    downAccTtbb = downInfo["nttbbV"]/downInfo["nttbbF"]
    downAccTtjj = downInfo["nttjjV"]/downInfo["nttjjF"]

    downRVis = subInfo["r"][0]
    downTtjjVis = subInfo["k"][0]
    downTtbbVis = downRVis*downTtjjVis

    downRFull = downRVis/downAccTtbb*downAccTtjj
    downTtjjFull = downTtjjVis/downAccTtjj
    downTtbbFull = downTtbbVis/downAccTtbb

  else: 
    upRVis = subInfo["r"][2]
    upTtjjVis = subInfo["k"][2]
    upTtbbVis = upRVis*upTtjjVis

    upRFull = upRVis/cenAttbb*cenAttjj
    upTtjjFull = upTtjjVis/cenAttjj
    upTtbbFull = upTtbbVis/cenAttbb

    downRVis = subInfo["r"][0]
    downTtjjVis = subInfo["k"][0]
    downTtbbVis = downRVis*downTtjjVis

    downRFull = downRVis/cenAttbb*cenAttjj
    downTtjjFull = downTtjjVis/cenAttjj
    downTtbbFull = downTtbbVis/cenAttbb

  tempSubSysFullR = max(abs(rFull-upRFull),abs(rFull-downRFull))
  tempSubSysFullTtjj = max(abs(ttjjFull-upTtjjFull),abs(ttjjFull-downTtjjFull))
  tempSubSysFullTtbb = max(abs(ttbbFull-upTtbbFull),abs(ttbbFull-downTtbbFull))

  print tempSubSysFullR, tempSubSysFullTtjj, tempSubSysFullTtbb

  tempSysFullR = pow(tempSysFullR*tempSysFullR+tempSubSysFullR*tempSubSysFullR,0.5)
  tempSysFullTtjj = pow(tempSysFullTtjj*tempSysFullTtjj+tempSubSysFullTtjj*tempSubSysFullTtjj,0.5)
  tempSysFullTtbb = pow(tempSysFullTtbb*tempSysFullTtbb+tempSubSysFullTtbb*tempSubSysFullTtbb,0.5)

  impactArray.append((subInfo["name"],subInfo["impact_r"]))
  pull = round(subInfo["fit"][1],2)
  if pull == -0.0: pull = 0.0
  pullSigP = abs(round(subInfo["fit"][2]-pull,2))
  pullSigM = abs(round(subInfo["fit"][0]-pull,2))
  valueRPercent = round(subInfo["impact_r"]/bestR*100.,2)
  valueKPercent = round(subInfo["impact_k"]/bestK*100.,2)
  isInList = False
  for name in nameList:
    if subInfo["name"]==name[0]: 
      print name[1], '&', pull, '& -%.2f/+%.2f' % (pullSigM, pullSigP), '&', valueRPercent, '&', valueKPercent, '&', round(pow(valueRPercent*valueRPercent+valueKPercent*valueKPercent,0.5),2),'\\\\'
      isInList = True
  if not isInList and 'prop' not in subInfo["name"]: print subInfo["name"], '&', pull, '& -%.2f/+%.2f' % (pullSigM, pullSigP), '&', valueRPercent, '&', valueKPercent, '\\\\'
  if 'prop' in subInfo["name"]: 
    mcstatRErr = pow(mcstatRErr*mcstatRErr + valueRPercent*valueRPercent, 0.5)
    mcstatKErr = pow(mcstatKErr*mcstatKErr + valueKPercent*valueKPercent, 0.5)
  if 'Btag' in subInfo["name"]: 
    btagRErr = pow(btagRErr*btagRErr + valueRPercent*valueRPercent, 0.5)
    btagKErr = pow(btagKErr*btagKErr + valueKPercent*valueKPercent, 0.5)
  if 'JES' in subInfo["name"]: 
    jesRErr = pow(jesRErr*jesRErr + valueRPercent*valueRPercent, 0.5)
    jesKErr = pow(jesKErr*jesKErr + valueKPercent*valueKPercent, 0.5)
print 'MC Statistics & - & - & ', round(mcstatRErr,2), '&', round(mcstatKErr,2), '&', round(pow(mcstatRErr*mcstatRErr + mcstatKErr*mcstatKErr,0.5),2),'\\\\'  
print 'b-tag total & - & - & ', round(btagRErr,2), '&', round(btagKErr,2), '\\\\'  
print 'Jet Energy Scale & - & - & ', round(jesRErr,2), '&', round(jesKErr,2), '&', round(pow(jesRErr*jesRErr + jesKErr*jesKErr,0.5),2),'\\\\'

for i,subInfo in enumerate(info["params"]):
  pull = round(subInfo["fit"][1],2)
  if pull == -0.0: pull = 0.0
  pullSigP = abs(round(subInfo["fit"][2]-pull,2))
  pullSigM = abs(round(subInfo["fit"][0]-pull,2))
  valueRPercent = round(subInfo["impact_r"]/bestR*100.,2)
  valueKPercent = round(subInfo["impact_k"]/bestK*100.,2)
  isInList = False
  for name in nameList:
    if subInfo["name"]==name[0]: 
      print name[1], '&', round(pow(valueRPercent*valueRPercent+valueKPercent*valueKPercent,0.5),2) ,'&', valueKPercent, '&', valueRPercent, '\\\\'
      isInList = True
  if not isInList and 'prop' not in subInfo["name"]: print subInfo["name"], '&', round(pow(valueRPercent*valueRPercent+valueKPercent*valueKPercent,0.5),2), '&', valueKPercent, '&', valueRPercent, '\\\\'
  if 'prop' in subInfo["name"]: 
    mcstatRErr = pow(mcstatRErr*mcstatRErr + valueRPercent*valueRPercent, 0.5)
    mcstatKErr = pow(mcstatKErr*mcstatKErr + valueKPercent*valueKPercent, 0.5)
  if 'Btag' in subInfo["name"]: 
    btagRErr = pow(btagRErr*btagRErr + valueRPercent*valueRPercent, 0.5)
    btagKErr = pow(btagKErr*btagKErr + valueKPercent*valueKPercent, 0.5)
  if 'JES' in subInfo["name"]: 
    jesRErr = pow(jesRErr*jesRErr + valueRPercent*valueRPercent, 0.5)
    jesKErr = pow(jesKErr*jesKErr + valueKPercent*valueKPercent, 0.5)
print 'MC Statistics & ', round(pow(mcstatRErr*mcstatRErr + mcstatKErr*mcstatKErr,0.5),2), '&',round(mcstatKErr,2), '&', round(mcstatRErr,2),  '\\\\'  
print 'b-tag total & ', round(btagRErr,2), '&', round(btagKErr,2), '\\\\'  
print 'Jet Energy Scale & ', round(pow(jesRErr*jesRErr + jesKErr*jesKErr,0.5),2), '&', round(jesKErr,2), '&', round(jesRErr,2) ,'\\\\'  

impactArray.sort(key=lambda x: x[1], reverse = True)
rankList = []
for i in range(31):
  rankList.append(impactArray[i][0])
print rankList

print "rFull = ", rFull, statFullR, tempSysFullR
print "ttjjFull = ", ttjjFull, statFullTtjj, tempSysFullTtjj 
print "ttbbFull = ", ttbbFull, statFullTtbb, tempSysFullTtbb 
