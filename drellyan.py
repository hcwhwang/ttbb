#!/usr/bin/env python
import ROOT, CMS_lumi, json, os, getopt, sys
from histoHelper import *
ROOT.gROOT.SetBatch(True)

datalumi = 35.82
#datalumi = 5.92
datalumi *= 1000

#rootfileDir = "/xrootd/store/user/chanwook/ttbb/final/preapproval/" 
rootfileDir = "./"
rootfileDirJES = "/xrootd/store/user/chanwook/ttbb/final/preapproval/JES/" 
#rootfileDir = "/xrootd/store/user/chanwook/ntuple/ReReco/" 
#rootfileDir = ""
analysis = "TtbarBbbarDiLeptonAnalyzer_"
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/totdataset.json" % os.environ['CMSSW_BASE']))

dataset = ['MuonEG_Run2016_tot','DoubleEG_Run2016_tot','DoubleMuon_Run2016_tot']
chlist = ["MuEl", "ElEl", "MuMu"]
stepList = [2, 3, 4,5]

#cut = '1'
cut = 'filtered==1 && tri>0'
weight = 'tri*weight*mueffweight*eleffweight*puweight'
weight1 = 'tri*weight*mueffweight*eleffweight*puweight*csvweights2[0]'

enttname = "cattree/nom"
tname = "cattree/nom2"
dyratio = [[1. for x in range(6)] for x in range(4)]
text = ""
muelesti = {}
elelesti = {}
mumuesti = {}
dyesti = [muelesti, elelesti, mumuesti]
for i,esti in enumerate(dyesti):
    esti["channel"] = chlist[i]
for i,step in enumerate(stepList):
    scale = 1.
    dycut = ""
    if step == 2: dycut = "(step1==1)*"
    if step == 3: dycut = "(step1==1)*(step3==1)*"
    if step == 4: dycut = "(step1==1)*(step3==1)*(step4==1)*"
    if step == 5: dycut = "(step1==1)*(step3==1)*(step4==1)*(step5==1)*"

    #rfname = rootfileDir + analysis + 'DYJets' +".root"
    rfname = rootfileDir + 'DYJets' +".root"
    data = findDataSet('DYJets', datasets)
    scale = datalumi*data["xsec"]
    if step<4: dyweight = weight
    else : dyweight = weight1
    wentries = getWeightedEntries(rfname, enttname, "filtered", "weight")
    scale = scale/wentries

   
    mc_ee_in = makeTH1(rfname,tname,"mc_ee_in", [2,0,2], "filtered", '('+dycut+'(%s && channel==2 && step2==0))*(%s)'%(cut,dyweight), scale)
    mc_mm_in = makeTH1(rfname,tname,"mc_mm_in", [2,0,2], "filtered", '('+dycut+'(%s && channel==3 && step2==0))*(%s)'%(cut,dyweight), scale)
    mc_ee_out = makeTH1(rfname,tname,"mc_ee_out", [2,0,2], "filtered", '('+dycut+'(%s && channel==2 && step2==1))*(%s)'%(cut,dyweight), scale)
    mc_mm_out = makeTH1(rfname,tname,"mc_mm_out", [2,0,2], "filtered", '('+dycut+'(%s && channel==3 && step2==1))*(%s)'%(cut,dyweight), scale)
    #rfname = rootfileDir + analysis + 'DYJets_10to50' +".root"
    rfname = rootfileDir + 'DYJets_10to50' +".root"
    data = findDataSet('DYJets_10to50', datasets)
    scale = datalumi*data["xsec"]

    wentries = getWeightedEntries(rfname, enttname, "filtered", "weight")    
    scale = scale/wentries

    mc_ee_in.Add(makeTH1(rfname,tname,"mc_ee_in", [2,0,2], "filtered", '('+dycut+'(%s && channel==2 && step2==0))*(%s)'%(cut,dyweight), scale))
    mc_mm_in.Add(makeTH1(rfname,tname,"mc_mm_in", [2,0,2], "filtered", '('+dycut+'(%s && channel==3 && step2==0))*(%s)'%(cut,dyweight), scale))
    mc_ee_out.Add(makeTH1(rfname,tname,"mc_ee_out", [2,0,2], "filtered", '('+dycut+'(%s && channel==2 && step2==1))*(%s)'%(cut,dyweight), scale))
    mc_mm_out.Add(makeTH1(rfname,tname,"mc_mm_out", [2,0,2], "filtered", '('+dycut+'(%s && channel==3 && step2==1))*(%s)'%(cut,dyweight), scale))

    rfname = rootfileDir + dataset[1-1]+".root"
    rd_em_in = makeTH1(rfname, tname,'rd_em_in', [2,0,2], "filtered", dycut+'(%s && channel==1 && ((ll_m > 76) && (ll_m < 106)))'%(cut))
    rfname = rootfileDir + dataset[2-1] +".root"
    rd_ee_in = makeTH1(rfname, tname,'rd_ee_in', [2,0,2], "filtered", dycut+'(%s && channel==2 && step2 ==0)'%(cut))
    rfname = rootfileDir + dataset[3-1] +".root"
    rd_mm_in = makeTH1(rfname, tname,'rd_mm_in', [2,0,2], "filtered", dycut+'(%s && channel==3 && step2 ==0)'%(cut))

    rfname = rootfileDir + dataset[2-1] +".root"
    rd_ee_in_loose = makeTH1(rfname, tname,'rd_ee_in_loose', [2,0,2], "filtered", '(step1==1)*(%s && channel==2 && step2 ==0)'%(cut))
    rfname = rootfileDir + dataset[3-1] +".root"
    rd_mm_in_loose = makeTH1(rfname, tname,'rd_mm_in_loose', [2,0,2], "filtered", '(step1==1)*(%s && channel==3 && step2 ==0)'%(cut))

    hList = [mc_ee_in, mc_ee_out, mc_mm_in, mc_mm_out, rd_ee_in, rd_mm_in, rd_em_in, rd_ee_in_loose, rd_mm_in_loose]
    hNList = []
    for h in hList:
        hNList.append(h.Integral())
        print h.GetName(), h.Integral()
    kEE = math.sqrt(hNList[7]/hNList[8])/2.
    kMM = math.sqrt(hNList[8]/hNList[7])/2.

    rMC_ee = hNList[1]/hNList[0]
    rMC_mm = hNList[3]/hNList[2]

    nOutEst_ee = rMC_ee*(hNList[4] - hNList[6]*kEE)
    nOutEst_mm = rMC_mm*(hNList[5] - hNList[6]*kMM)

    scale_ee = nOutEst_ee/hNList[1]
    scale_mm = nOutEst_mm/hNList[3]
   
    print "DY estimation for", step, "ee =",scale_ee, "mm =",scale_mm
    print "kMM = ", kMM, "kEE = ", kEE
    dyratio[1][step] = pow(scale_ee*scale_mm, 0.5)
    dyratio[2][step] = scale_ee
    dyratio[3][step] = scale_mm

for ch,esti in enumerate(dyesti):
    for st in range(1,6):
        esti["step%s"%st] = dyratio[ch+1][st]
result = [muelesti, elelesti, mumuesti]
import json
with open("drellyanresult.json","w") as outfile:
    json.dump(result,outfile)    
print result 
