#!/usr/bin/env/python
import CMS_lumi, json, os, getopt, sys, copy
from ROOT import *
from histoHelper import *
from sysWeightPlot_cfi import *
#from sysWeightPlot_cfi import *
ROOT.gROOT.SetBatch(True)
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
redraw = False
plottitle = 'invariant_mass'
opts, arv = getopt.getopt(sys.argv[1:],"hdrc:w:b:p:x:y:s:t:",["cut","weight","binning","plotvar","x_name","y_name","dolog","step","title","redraw"])
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
  elif opt in ("-r", "--redraw"):
    redraw = True
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

#datalumi = 1.56
datalumi = 35.8
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1} ( 13 TeV ) "%(datalumi)
datalumi = datalumi*1000 # due to fb

analysis = 'TtbarBbbarDiLeptonAnalyzer_'
mcfilelist = ['TT_powheg', 'WJets','SingleTbar_tW', 'SingleTbar_t', 'SingleTop_t','SingleTop_tW', 'SingleTop_s', 'TTW', 'TTZ', 'ZZ', 'WW', 'WZ', 'WWW', 'WWZ','WZZ', 'ZZZ', 'DYJets','DYJets_10to50']
rdfilelist = ['MuonEG_Run2016_tot','DoubleEG_Run2016_tot','DoubleMuon_Run2016_tot']
#rootfileDir = "/xrootd/store/user/chanwook/ttbb/final/preapproval/"
#rootfileDir = "/cms/ldap_home/chanwook/work/CatProduction/cattools/src/CATTools/CatAnalyzer/test/rootFile/"
rootfileDir = "./"
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



visible="(NJets30>=4 && NbJets30>=2 && ((lepton1_pt>25 && lepton2_pt>20)||(lepton1_pt>20 && lepton2_pt>25)) && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"
ttbb = mAND("(NbJets30>=4)",visible)
ttbj = mAND("(NbJets30==3)",visible)
ttcc = mAND("((NcJets30>=2) && !(NbJets30>=3))",visible)
ttlf = mAND("(!(NbJets30>=4) && !(NbJets30==3) && !(NcJets30>=2))",visible)

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
ttbarid_title = ["t#bar{t}b#bar{b}", "t#bar{t}bj", "t#bar{t}c#bar{c}", "t#bar{t}LF", "t#bar{t} others"]
ttbarid_plottitle = ["ttbb","ttbj","ttcc","ttLF","ttothers"]
ttbarid_fillcolor = ["#660000", "#ffcc00", "#cc6600", "#ff0000", "#ff6565"]
ttbarid_linecolor = ["#000000", "#000000", "#000000", "#000000", "#000000"]
datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/totdataset.json" % os.environ["CMSSW_BASE"]))
#plotdir = "%s/src/CATTools/CatAnalyzer/test/preapproval/" % os.environ["CMSSW_BASE"]
plotdir = "%s/src/CATTools/CatAnalyzer/test/preapproval/test/" % os.environ["CMSSW_BASE"]
 
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
plotTitle = str(step)+plottitle
plotName = str(step)+x_name.replace(" ","")
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
  K = 2.171/2.411
  initVisR = 0.01318
  fitVisR = 0.018
  initRttbb = [0.0291, 0.0254, 0.0280]
  initRttbj = [0.0886, 0.0836, 0.0908]
  effttbb = [0.195, 0.041, 0.084]
  effttjj = [0.088, 0.022, 0.040]
  fitRttbb = [fitVisR*effttbb[0]/effttjj[0], fitVisR*effttbb[1]/effttjj[1], fitVisR*effttbb[2]/effttjj[2]]
  fitRttbj = [fitRttbb[0]*initRttbj[0]/initRttbb[0], fitRttbb[1]*initRttbj[1]/initRttbb[1],fitRttbb[2]*initRttbj[2]/initRttbb[2]]
  scalettbb = [K*fitRttbb[0]/initRttbb[0], K*fitRttbb[1]/initRttbb[1],K*fitRttbb[2]/initRttbb[2]]
  scalettbj = [K*fitRttbj[0]/initRttbj[0], K*fitRttbj[1]/initRttbj[1], K*fitRttbj[2]/initRttbj[2]]
  scalettcclf = [K*(1-fitRttbb[0]-fitRttbj[0])/(1-initRttbb[0]-initRttbj[0]), K*(1-fitRttbb[1]-fitRttbj[1])/(1-initRttbb[1]-initRttbj[1]), K*(1-fitRttbb[2]-fitRttbj[2])/(1-initRttbb[2]-initRttbj[2])]  
  scalettothers = [K, K, K]

  fitScaleList = [scalettbb, scalettbj, scalettcclf, scalettcclf, scalettothers]  
  print fitScaleList
  '''
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
  '''
if len(binning)==3: is1D=True
elif len(binning)==6: is1D=False

if redraw == False:
  rootFile = TFile.Open("%s.root"%plotTitle, "recreate")
  rootFile.Close()
for ch in range(3):
  if redraw == True:continue 
  #if ch!=2: continue
  for dych in dyscale:
    if dych["channel"] == channel_name[ch]:
      dy = dych   

  mchistList = []
  tthistList = []
  mcscaleList = []
  stepch_tcut =  mAND(scut,'channel==%i'%(ch+1))
  
  rfname = rootfileDir +  rdfilelist[ch] +".root"
  rdcut = '(%s&&%s)'%(stepch_tcut,cut)
  if is1D:rdhist = makeTH1(rfname, tname, 'data', binning, plotvar, rdcut)
  else: rdhist = makeTH2toTH1(rfname, tname, 'data', binning, plotvar, rdcut)
  rdhist.SetLineColor(1)
  rdhist.SetName(channel_name[ch])
  print channel_name[ch]
  histRoot = TFile.Open("%s.root" % plotTitle, "update")
  rdhist.Write()
  histRoot.Close()
  if is1D: hbkgall = ROOT.TH1D("hbkgall%i_%i" % (ch+1,s),"",binning[0],binning[1],binning[2])
  else: hbkgall = ROOT.TH1D("hbkgall%i_%i" % (ch+1,s),"",binning[0]*binning[3],0,binning[0]*binning[3])
  for i, mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    scale = datalumi*data["xsec"]
    colour = data["colour"]
    title = data["title"]
    if 'DYJets' in mcname and step >= 2: 
      scale = scale*dy["step%s"%step]
    #rfname = rootfileDir + analysis+mcname +".root"
    rfname = rootfileDir +mcname +".root"
    if mcname+".root" not in sorted(os.listdir(rootfileDir)): continue
    wentries = getWeightedEntries(rfname, tnameEnt, "filtered", 'weight')
  
    scale = scale/wentries
    fitScale = 1.
    if mcname == 'TT_powheg':
      wentries1 = getWeightedEntries(rootfileDir+"TT_powheg.root", tnameEnt, "filtered", 'weight')
      wentries2 = getWeightedEntries(rootfileDir+"TTLL_powheg.root", tnameEnt, "filtered", 'weight')
      scale = datalumi*data['xsec']/(wentries1 + wentries2*(1./0.1049))
      mcscaleList.append(scale)
      rfname = rootfileDir + 'TT_powheg_tot.root'
      for i,ttbarId in enumerate(ttbarid_cut):
        ttbarcut = mAND(cut,ttbarid_cut[i])
        tcut = '(%s&&%s)*%s'%(stepch_tcut,ttbarcut,weight)
        if 'afterfit' in plottitle: fitScale = fitScaleList[i][ch]
        if is1D: mchist = makeTH1(rfname, tname, ttbarid_title[i], binning, plotvar, tcut, scale*fitScale)
        else: mchist = makeTH2toTH1(rfname, tname, ttbarid_title[i], binning, plotvar, tcut, scale*fitScale)
        ttPlotName = "%s_%s" % (channel_name[ch], ttbarid_plottitle[i])
        print ttPlotName
        mchist.SetName(ttPlotName)
        #mchist.SetTitle(ttPlotName)
        histRoot = TFile.Open("%s.root" % plotTitle, "update")
        mchist.Write()
        histRoot.Close()
    #elif 'TT' not in mcname:
    else:
      mcscaleList.append(scale)
      tcut = '(%s&&%s)*%s'%(stepch_tcut,cut,weight)
      if is1D: mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
      else: mchist = makeTH2toTH1(rfname, tname, title, binning, plotvar, tcut, scale)
      mcPlotName = "%s_%s" % (channel_name[ch], mcname)
      print mcPlotName
      mchist.SetName(mcPlotName)
      #mchist.SetTitle(mcPlotName)
      histRoot = TFile.Open("%s.root" % plotTitle, "update")
      mchist.Write()
      histRoot.Close()
      hbkgall.Add(mchist)

  hbkgall.SetName("%s_bkg" % channel_name[ch])
  histRoot = TFile.Open("%s.root" % plotTitle, "update")
  hbkgall.Write()
  histRoot.Close()    
  mcsysList = []

  #if 'afterfitt' not in plottitle: 
  if 'csv' not in plottitle: 
    for s,systErr in enumerate(sysweights):
      htitle = "hmcall%i_%i" % (ch+1,s) 
      if is1D: hmcall = ROOT.TH1D(htitle,"",binning[0],binning[1],binning[2])
      else: hmcall = ROOT.TH1D(htitle,"",binning[0]*binning[3],0,binning[0]*binning[3])
      for m, mc in enumerate(mcfilelist):
        if mc in ['TT_aMC', 'TT_powheg_herwig']: continue
        hmctitle = htitle+mc
        #fname = rootfileDir + analysis + mc +".root"
        if mc=='TT_powheg':
          if systErr["name"] == 'FSR_Down':
            mc = 'TT_powheg_fsrdown'
          elif systErr["name"] == 'FSR_Up':
            mc = 'TT_powheg_fsrup'
          elif systErr["name"] == 'ISR_Down':
            mc = 'TT_powheg_isrdown'
          elif systErr["name"] == 'ISR_Up':
            mc = 'TT_powheg_isrup'
          elif systErr["name"] == 'UE_Down':
            mc = 'TT_powheg_down'
          elif systErr["name"] == 'UE_Up':
            mc = 'TT_powheg_up'
          elif systErr["name"] == 'hdamp_Down':
            mc = 'TT_powheg_hdampdown'
          elif systErr["name"] == 'hdamp_Up':
            mc = 'TT_powheg_hdampup'
          elif systErr["name"] == 'erdon_Up':
            mc = 'TT_powheg_erdon'
          #elif systErr["name"] == 'erdon_Down':
          #  mc = 'TT_powheg_erdon'
          elif systErr["name"] == 'qcderdon_Up':
            mc = 'TT_powheg_qcderdon'
          #elif systErr["name"] == 'qcderdon_Down':
          #  mc = 'TT_powheg_qcderdon'
          elif systErr["name"] == 'gluonmove_Up':
            mc = 'TT_powheg_gluonmove'
          #elif systErr["name"] == 'gluonmove_Down':
          #  mc = 'TT_powheg_gluonmove'
          else: mc = 'TT_powheg_tot'

        fname = rootfileDir + mc +".root"
        tree = "cattree/" + systErr["tree"]
        wcut = '(%s&&%s)*%s'%(stepch_tcut,cut,systErr["var"])                 
        if 'MuR' in systErr["name"] or 'MuF' in systErr["name"] or 'PDF' in systErr['name']: 
          if 'TT_powheg' in mc: mc = 'TT_powheg'
          data = findDataSet(mc, datasets)
          scale = datalumi*data["xsec"]
          scaleW = systErr["var"].strip(")").split('*')[-1]
          scaleW = 'weight*%s' % scaleW

          if mc == 'TT_powheg':         
            wentries1 = getWeightedEntries(rootfileDir+"TT_powheg.root", tnameEnt, "filtered", scaleW)
            wentries2 = getWeightedEntries(rootfileDir+"TTLL_powheg.root", tnameEnt, "filtered", scaleW)
            scale = scale/(wentries1 + wentries2*(1./0.1049))
          else:
            wentries = getWeightedEntries(fname, tnameEnt, "filtered", "weight")
            scale = scale/wentries
            wcut = '(%s&&%s)*%s'%(stepch_tcut,cut,weight)
          if is1D: mch = makeTH1(fname, tree, hmctitle, binning, plotvar, wcut, scale)
          else: mch = makeTH2toTH1(fname, tree, hmctitle, binning, plotvar, wcut, scale)
        elif 'TT_powheg_' in mc and 'tot' not in mc:
          data = findDataSet('TT_powheg', datasets)
          scale = datalumi*data["xsec"]
          wentries = getWeightedEntries(fname, tnameEnt, "filtered", 'weight')
          ttScale = scale/wentries
          if is1D: mch = makeTH1(fname, tree, hmctitle, binning, plotvar, wcut, ttScale)
          else: mch = makeTH2toTH1(fname, tree, hmctitle, binning, plotvar, wcut, ttScale)
        else:
          if is1D: mch = makeTH1(fname, tree, hmctitle, binning, plotvar, wcut, mcscaleList[m])
          else: mch = makeTH2toTH1(fname, tree, hmctitle, binning, plotvar, wcut, mcscaleList[m])
        hmcall.Add(mch)
      if s==0: hCenter = copy.deepcopy(hmcall)
      if s!=0: hmcall.Add(hCenter,-1)
      mcSysName = "%s_sys_%s" % (channel_name[ch], systErr["name"])
      print mcSysName
      hmcall.SetName(mcSysName)
      histRoot = TFile.Open("%s.root" % plotTitle, "update")
      hmcall.Write()
      histRoot.Close()

histRoot = TFile.Open("%s.root" % plotTitle, "read")
mchisttot = []
tthisttot = []
for ch in range(3):
  mchistList = []
  tthistList = []

  for i, mcname in enumerate(mcfilelist):
    data = findDataSet(mcname, datasets)
    colour = data["colour"]
    title = data["title"]

    if mcname == 'TT_powheg':
      for i, ttbarId in enumerate(ttbarid_cut):
        h = histRoot.Get("%s_%s" % (channel_name[ch],ttbarid_plottitle[i])) 
        lC = TColor.GetColor(ttbarid_linecolor[i])
        fC = TColor.GetColor(ttbarid_fillcolor[i])
        h.SetLineColor(lC) 
        h.SetFillColor(fC) 
        mchistList.append(h)
    else:
      h = histRoot.Get("%s_%s" % (channel_name[ch], mcname))
      h.SetLineColor(colour)
      h.SetFillColor(colour)
      mchistList.append(h) 
  mcsysList = []
  i = 0
  for h in histRoot.GetListOfKeys():
    if h.GetClassName()!='TH1D': continue
    if 'sys' not in h.GetName(): continue
    if channel_name[ch] not in h.GetName(): continue
    #if 'csvweight' not in h.GetName(): continue
    #if 'Btag' not in h.GetName() and 'csvweight' not in h.GetName(): continue
    sysHist = histRoot.Get(h.GetName())
    mcsysList.append(sysHist)
    #if i<=30:mcsysList.append(sysHist)
    i+=1
  print mcsysList
  rdhist = histRoot.Get(channel_name[ch])
  rdhist.SetLineColor(1) 
  plotDirName = plotdir + plotTitle + channel_name[ch]
  hbkgall = histRoot.Get("%s_bkg" % channel_name[ch])
  name = xyname(x_name, y_name, step, ch, binning)
  drawTH1(plotName+channel_name[ch], plotDirName, CMS_lumi, mchistList, tthistList, hbkgall, rdhist, name[0], name[1], mcsysList,dolog)
  #drawTH1(outfile, plotname, CMS_lumi, mchistList, tthistList, hbkgall, rdhist, name[0], name[1], dolog)
  if ch==0:
    histList = []
    rdhisttot = copy.deepcopy(rdhist)
    histList.append(rdhisttot)
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
    histList.append(hbkgalltot)
  else:
  #elif ch==1 or ch==2:
    rdhisttot.Add(copy.deepcopy(rdhist))
    for i in range(len(mchistList)):
      mchisttot[i].Add(copy.deepcopy(mchistList[i]))
    for j in range(len(mcsysList)):
      mcsystot[j].Add(copy.deepcopy(mcsysList[j]))
    for k in range(len(tthistList)):
      tthisttot[k].Add(copy.deepcopy(tthistList[k]))
    hbkgalltot.Add(copy.deepcopy(hbkgall))
plotDirName = plotdir + plotTitle + "All" 
name = xyname(x_name, y_name, step, 3, binning)
drawTH1(plotName, plotDirName, CMS_lumi, mchisttot, tthisttot, hbkgalltot, rdhisttot, name[0], name[1], mcsystot, dolog)
histRoot.Close()
