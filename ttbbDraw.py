#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from ROOT import *
from histoHelper import *
from sysWeightPlot_cfi import *
ROOT.gROOT.SetBatch(True)
'''
topDraw.py -a 1 -s 1 -c 'tri==1&&filtered==1' -b [40,0,40] -p nvertex -x 'no. vertex' &
topDraw.py -a 1 -s 1 -b [100,-3,3] -p lep1_eta,lep2_eta -x '#eta' &
'''
#datalumi = 1.56
datalumi = 35.8
#datalumi = 2.65
#datalumi = 4.35
#datalumi = 3.05
#datalumi = 8.53
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1} ( 13 TeV ) "%(datalumi)
datalumi = datalumi*1000 # due to fb

analysis = 'TtbarBbbarDiLeptonAnalyzer_'
#mcfilelist = ['TT_powheg', "WJets",'SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW','WWZ','WZZ','ZZZ','ttW','ttZ','DYJets', 'DYJets_10to50']
#mcfilelist = ['TT_powheg','SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s','ZZ', 'WW', 'WZ', 'WWW','WWZ','WZZ','ZZZ','ttW','ttZ','DYJets','DYJets_10to50']
mcfilelist = ['TT_powheg', 'WJets','SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s', 'TTW', 'TTZ', 'ZZ', 'WW', 'WZ', 'WWW', 'WWZ','WZZ', 'ZZZ', 'DYJets','DYJets_10to50']
#mcfilelist = ['DYJets','DYJets_10to50']
#rdfilelist = ['MuonEG_Run2016','DoubleEG_Run2016','DoubleMuon_Run2016']

rdfilelist = ['MuonEG_Run2016','DoubleEG_Run2016','DoubleMuon_Run2016']
rootfileDir = "/xrootd/store/user/chanwook/ttbb/final/final/"
#rootfileDir = ""
channel_name = ['MuEl', 'ElEl', 'MuMu', "All"]

def mAND(aaa,bbb):
  return "(" +aaa+ " && "+bbb+")"
def mAND2(aaa):
  bbb=""
  for i,ii in enumerate(aaa):
    if i==0 : 
      bbb+= ii
    else : 
      bbb=mAND(ii,bbb)
  return bbb
def op_(aaa):
  return "!(" + aaa + ")"




visible="(NJets20>=4 && NbJets20>=2 && lepton1_pt>20 && lepton2_pt>20 && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"
ttbb = mAND("(NbJets20>=4)",visible)
ttbj = mAND("(NbJets20==3)",visible)
ttcc = mAND("((NcJets20>=2) && !(NbJets20>=3))",visible)
ttlf = mAND("(!(NbJets20>=4) && !(NbJets20==3) && !(NcJets20>=2))",visible)

#full phase
#fullphase ="(diLeptonicM1==1 && NaddJets20 >= 2)"
#TTJJ = "(NaddJets20 >= 2 && diLeptonicM1==1)"
#ttbb = "(NaddbJets20 >= 2 && diLeptonicM1==1)"
#ttb = "(NaddJets20 >= 2 && NaddbJets20 == 1 && diLeptonicM1==1 && !(genTtbarId%100==52))"
#tt2b = "(NaddJets20 >= 2 && NaddbJets20 == 1 && diLeptonicM1==1 && (genTtbarId%100==52))"
#ttcc = "(NaddJets20 >= 2 && NaddcJets20 >= 2 && NaddbJets20==0 && diLeptonicM1==1)"
#ttlf = "( !"+ttbb+" && !"+ttb+" && !"+ttcc+"  && NaddJets20 >= 2 && diLeptonicM1==1)"
ttothers = op_(visible)

ttbarid_cut = [ttbb,ttbj,ttcc,ttlf,ttothers]
#ttbarid_cut = ['&&(genTtbarId%100>52&&genTtbarId%100<56)','&&(genTtbarId%100>50&&genTtbarId%100<53)', '&&(genTtbarId%100>42&&genTtbarId%100<46)', '&&(genTtbarId%100>40&&genTtbarId%100<43)', '&&genTtbarId%100==0', '&&((genTtbarId%100>0&&genTtbarId%100<41)||(genTtbarId%100>45&&genTtbarId%100<51)||(genTtbarId%100>55))']
ttbarid_title = ["t#bar{t}b#bar{b}", "t#bar{t}bj", "t#bar{t}c#bar{c}", "t#bar{t}lf", "t#bar{t}others"]
ttbarid_fillcolor = ["#660000", "#ffcc00", "#cc6600", "#ff0000", "#ff6565"]
ttbarid_linecolor = ["#000000", "#000000", "#000000", "#000000", "#000000"]
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/totdataset.json" % os.environ["CMSSW_BASE"]))
plotdir = "%s/src/CATTools/CatAnalyzer/test/ControlPlots/" % os.environ["CMSSW_BASE"]
#defalts
step = 1
#channel = 1

#cut = 'tri!=0&&filtered==1'
cut = '1'
#weight = 'weight*puweight'
weight = '1'
binning = [60, 20, 320]
plotvar = 'll_m'
x_name = 'mass [GeV]'
y_name = 'Number of Events'
dolog = False
plottitle = 'invariant_mass'
try:
  opts, args = getopt.getopt(sys.argv[1:],"hdc:w:b:p:x:y:s:t:",["cut","weight","binning","plotvar","x_name","y_name","dolog","step","title"])
except getopt.GetoptError:          
  print 'Usage : ./topDraw.py.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -t <title>'
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print 'Usage : ./topDraw.py.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -t <title>'
    sys.exit()
  elif opt in ("-c", "--cut"):
    cut = arg
  #elif opt in ("-a", "--channel"):
  #    channel = int(arg)
  elif opt in ("-s", "--step"):
    step = int(arg)
  elif opt in ("-w", "--weight"):
    weight = arg
  elif opt in ("-b", "--binning"):
    print arg
    binning = eval(arg)
  elif opt in ("-p", "--plotvar"):
    plotvar = arg
  elif opt in ("-x", "--x_name"):
    x_name = arg
  elif opt in ("-y", "--y_name"):
    y_name = arg
  elif opt in ("-d", "--dolog"):
    dolog = True
  elif opt in ("-t", "--title"):
    plottitle = arg
if plotvar == 'ttbar_deltaphi':
  plotvar = 'acos(cos(top1_phi-top2_phi))'
fitFile = "fitresultNbjetsL3.txt"
  
if 'puweight' not in weight: sysweights = mceventweightNOPU
elif 'csv' in weight: sysweights = mceventweightCSV
else : sysweights = mceventweight

print sysweights
tnameEnt = "cattree/nom"
tname = "cattree/nom2"
mchistList = []

#if channel == 1: ttother_tcut = "!(parton_channel==2 && ((parton_mode1==1 && parton_mode2==2) || (parton_mode1==2 && parton_mode2==1)))"
#elif channel == 2: ttother_tcut = "!(parton_channel==2 && (parton_mode1==2 && parton_mode2==2))"
#elif channel == 3: ttother_tcut = "!(parton_channel==2 && (parton_mode1==1 && parton_mode2==1))"

#ttother_tcut = '(%s&&%s&&%s)*%s'%(stepch_tcut,cut,ttother_tcut,weight)
def xyname(x_name, y_name, step, channel, binning):
  print binning
  #x_name = "Step "+str(step)+" "+channel_name[channel]+" "+x_name
  if len(binning) <= 3:
    num = (binning[2]-binning[1])/float(binning[0])
    if num != 1:
      if x_name.endswith(']'):
          unit = "["+x_name.split('[')[1]
      else: unit = ""
      y_name = y_name + " / %g%s"%(num,unit)
  return x_name, y_name

print plotvar
print plottitle
var = plotvar.split(',')[0]
for i in range(4):
  if "[%s]"%i in var:
    var = var[:len(var)-3]
    if i ==0: var = "1st_"+var
    if i ==1: var = "2nd_"+var
    if i ==2: var = "3rd_"+var
    if i ==3: var = "4th_"+var


dyscale = json.load(open("drellyanresult.json"))

mchistallList = []
rdhistallList = []
mcsysallList = []

mchisttot = []
tthisttot = []
mcsystot = []
if step==0: scut = "1"
else:
  for s in range(1,step+1):
    if (s==1): scut = "step1==1"
    else : scut = mAND(scut, "step%s==1" % (s))

summchistList = []
if 'afterfit' in plottitle:
  print "afterfit"
  fitF = open(fitFile,"r")
  line = fitF.readline()
  v = line.split(" ")
  rttbb = float(v[1])
  rttbbErr = float(v[2])
  K = float(v[3])
  KErr = float(v[4])
  rttbj = float(v[5])
  preRttbb = float(v[6])
  preRttbj = float(v[7]) 
  fitF.close()
for ch in range(3):
  name = xyname(x_name, y_name, step, ch, binning)
  for dych in dyscale:
    if dych["channel"] == channel_name[ch]:
      dy = dych   

  mchistList = []
  tthistList = []
  mcscaleList = []
  stepch_tcut =  mAND(scut,'channel==%i'%(ch+1))
  hbkgall = ROOT.TH1D("hbkgall%i_%i" % (ch+1,s),"",binning[0],binning[1],binning[2])
  for i, mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    title = data["title"]
    if 'DYJets' in mcname and step >= 2: 
      if step <= 4: scale = scale*dy["step%s"%step]
      else: scale = scale*dy["step4"]
    #rfname = rootfileDir + analysis+mcname +".root"
    rfname = rootfileDir +mcname +".root"
    if mcname+".root" not in sorted(os.listdir(rootfileDir)): continue
    wentries = getWeightedEntries(rfname, tnameEnt, "filtered", 'weight')
    print mcname, " ", wentries
  
    scale = scale/wentries
    mcscaleList.append(scale)
    if mcname == 'TT_powheg':
      for i,ttbarId in enumerate(ttbarid_cut):
        ttbarcut = mAND(cut,ttbarid_cut[i])
        tcut = '(%s&&%s)*%s'%(stepch_tcut,ttbarcut,weight)
        if 'afterfit' in plottitle and i==0:
          scale *= K*rttbb/preRttbb
        elif 'afterfit' in plottitle and i==1:
          scale *= K*rttbj/preRttbj
        elif 'afterfit' in plottitle and (i==2 or i==3):
          scale *= K*(1-rttbb-rttbj)/(1-preRttbb-preRttbj)
        elif 'afterfit' in plottitle and i==4:
          scale *= K
        mchist = makeTH1(rfname, tname, ttbarid_title[i], binning, plotvar, tcut, scale)
        lC = TColor.GetColor(ttbarid_linecolor[i])
        fC = TColor.GetColor(ttbarid_fillcolor[i])

        mchist.SetLineColor(lC)
        mchist.SetFillColor(fC)
        mchistList.append(mchist)
    elif 'tt' not in mcname:
      tcut = '(%s&&%s)*%s'%(stepch_tcut,cut,weight)
      mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
      mchist.SetLineColor(colour)
      mchist.SetFillColor(colour)
      mchistList.append(mchist)
      hbkgall.Add(mchist)
    else :
      tcut = '(%s&&%s)*%s'%(stepch_tcut,cut,weight)
      if 'afterfit' in plottitle:
        fitF = open(fitFile,"r")
        while True:
          line = fitF.readline()
          if not line: break
          v = line.split(" ")
          sysName = str(v[0])
          if sysName == mcname: 
            KSys = float(v[3])
            scale *= KSys
        fitF.close()
      if mcname == "TT_aMC":
        tthist = makeTH1(rfname, tname, "aMC@NLO + PHYTHIA8", binning, plotvar, tcut, scale)
        tthist.SetLineColor(kBlue)
        tthist.SetLineStyle(2)
        tthistList.append(tthist)
      elif mcname == "TT_powheg_herwig":
        tthist = makeTH1(rfname, tname, "POWHEG + HERWIG", binning, plotvar, tcut, scale)
        tthist.SetLineColor(kGreen)
        tthist.SetLineStyle(3)
        tthistList.append(tthist)
      
  mcsysList = []

  for s,systErr in enumerate(sysweights):
    if 'afterfit' in plottitle:
      fitF = open(fitFile,"r")
      while True:
        line = fitF.readline()
        if not line: break
        v = line.split(" ")
        sysName = str(v[0])
        if sysName == systErr["name"]: 
          KSys = float(v[3])
          mcscaleList[0] *= KSys  
      fitF.close()
    htitle = "hmcall%i_%i" % (ch+1,s) 
    hmcall = ROOT.TH1D(htitle,"",binning[0],binning[1],binning[2])
    print "sys : " + systErr["name"]
    for m, mc in enumerate(mcfilelist):
      if mc in ['TT_aMC', 'TT_powheg_herwig']: continue
      hmctitle = htitle+mc
      #fname = rootfileDir + analysis + mc +".root"
      fname = rootfileDir + mc +".root"
      tree = "cattree/" + systErr["tree"]
      wcut = '(%s&&%s)*%s'%(stepch_tcut,cut,systErr["var"])                 
      mch = makeTH1(fname, tree, hmctitle, binning, plotvar, wcut, mcscaleList[m])
      hmcall.Add(mch)
    if s!=0: hmcall.Add(mcsysList[0],-1)
    mcsysList.append(hmcall)
  #rfname = rootfileDir + analysis + rdfilelist[ch] +".root"
  rfname = rootfileDir +  rdfilelist[ch] +".root"
  rdcut = '(%s&&%s)'%(stepch_tcut,cut)
  rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, rdcut)
  rdhist.SetLineColor(1)
  #var = ''.join(i for i in var if not i.isdigit())
  outfile = "%s_s%d_%s"%(channel_name[ch],step,var)
  plotname = plotdir + "%s_s%d_%s"%(channel_name[ch],step,plottitle)
  drawTH1(outfile, plotname, CMS_lumi, mchistList, tthistList, hbkgall, rdhist, name[0], name[1], mcsysList,dolog)
  if ch==0:
    rdhisttot = copy.deepcopy(rdhist)
    for i in range(len(mchistList)):
      mchist = copy.deepcopy(mchistList[i])
      mchisttot.append(mchist)
    for j in range(len(mcsysList)):
      mcsys = copy.deepcopy(mcsysList[j])
      mcsystot.append(mcsys)
    for k in range(len(tthistList)):
      tthist = copy.deepcopy(tthistList[k])
      tthisttot.append(tthist)
    hbkgalltot = copy.deepcopy(hbkgall) 
  else:
    rdhisttot.Add(rdhist)
    for i in range(len(mchistList)):
      mchisttot[i].Add(mchistList[i])
    for j in range(len(mcsysList)):
      mcsystot[j].Add(mcsysList[j])
    for k in range(len(tthistList)):
      tthisttot[k].Add(tthistList[k])
    hbkgalltot.Add(hbkgall)
outfile = "%s_s%d_%s"%("All",step,var)
plotname = plotdir + "All_s%d_%s"%(step,plottitle)
name = xyname(x_name, y_name, step, 3, binning)
drawTH1(outfile, plotname, CMS_lumi, mchisttot, tthisttot, hbkgalltot, rdhisttot, name[0], name[1], mcsystot, dolog)

