import json, os, CMS_lumi
from histoHelper import *
from ROOT import * 
from array import array
from sysWeightPDF_cfi import *

gROOT.SetBatch(True)
def mAND(aaa,bbb):
  return "(" +aaa+ " && "+bbb+")"
def weightCut(a, b):
  cut = "(%s)*(%s)" % (a, b)
  return cut
def bOverAErr(b,bErr,a,aErr):
  return  pow(aErr*aErr*b*b/a/a + bErr*bErr, 0.5)/a
def abOverCErr(a,aErr,b,bErr,c,cErr):
  return pow(b*b/c/c*aErr*aErr + a*a/c/c*bErr*bErr + a*a*b*b/c/c/c/c*cErr*cErr,0.5)
lumi = 36.8
lumiText = "%.2f fb^{-1}, #sqrt{s} = 13 TeV "%(lumi)
lumi *= 1000
ttXsec = 831.76
ttXsecErr = 45.62
nbjets="Nbjets2"
#nbjets="NbjetsL3"
#jsonFile = "outputNbjets%i.json"%nbjets
#resultFile = "fitresultNbjets%i.txt"%nbjets
jsonFile = "output%s.json"%nbjets
pdfjsonFile = "output%s_pdf.json"%nbjets
#pdfjsonFile = "outputNbjetsL3_pdf.json"
resultFile = "fitresult%s.txt"%nbjets
resultPDFFile = "fitresult%s_pdf.txt"%nbjets
#resultPDFFile = "fitresultNbjetsL3_pdf.txt"
mcFile = "mcresultsNbjets3.json"
recalculateMC = False 
entries = json.load(open("%s/src/CATTools/CatAnalyzer/test/%s" % (os.environ["CMSSW_BASE"], jsonFile)))
pdfentries = json.load(open("%s/src/CATTools/CatAnalyzer/test/%s" % (os.environ["CMSSW_BASE"], pdfjsonFile)))
totSys = 0.
sysNameList = ["Pileup","JES \\& JER", "b tag (b quark flavour)", "b tag (c quark flavour)", "b tag (light flavour)", "Ratio of $t\\overline{t}b\\overline{b}$ and $t\\overline{t}bj$", "Background modelling", "$t\\overline{t}c\\overline{c}$ fraction in the fit", "Lepton", "MC generator", "Scale $(\\mu_{F}$ and $\\mu_{R})$", "PS scale", "PDF", "Efficiency $(t\\overline{t}c\\overline{c} fraction)$", "Luminosity","Trigger"]
#sysGroupList = [["PW_Up","PW_Down"],["JER_NOM","JER_Up","JER_Down","JES_Up","JES_Down"],["HF_Up","HF_Down","HF_Stats1_Up","HF_Stats1_Down","HF_Stats2_Up","HF_Stats2_Down"],["CQ_Err1_Up","CQ_Err1_Down","CQ_Err2_Up","CQ_Err2_Down"],["LF_Up","LF_Down","LF_Stats1_Up","LF_Stats1_Down","LF_Stats2_Up","LF_Stat2_Down"],["Free_ttbj"],["Bkg_Up","Bkg_Down"],[],["Mu_Eff_Up","Mu_Eff_Down","Mu_Pt_Up","Mu_Pt_Down","El_Eff_Up","El_Eff_Down","El_Pt_Up","El_Pt_Down"],["TT_powheg_herwig","TTJets_aMC"],[],[],[],[],[],["Trig_Up","Trig_Down"]]
sysGroupList = [[["PW_Up","PW_Down"]],[["JER_NOM","JER_Up","JER_Down"],["JES_Up","JES_Down"]],[["HF_Up","HF_Down"],["HF_Stats1_Up","HF_Stats1_Down"],["HF_Stats2_Up","HF_Stats2_Down"]],[["CQ_Err1_Up","CQ_Err1_Down"],["CQ_Err2_Up","CQ_Err2_Down"]],[["LF_Up","LF_Down"],["LF_Stats1_Up","LF_Stats1_Down"],["LF_Stats2_Up","LF_Stat2_Down"]],[["Free_ttbj"]],[["Bkg_Up","Bkg_Down"]],["Ttcc_Fraction_Fit_Up","Ttcc_Fraction_Fit_Down"],[["Mu_Eff_Up","Mu_Eff_Down"],["Mu_Pt_Up","Mu_Pt_Down"],["El_Eff_Up","El_Eff_Down"],["El_Pt_Up","El_Pt_Down"]],[],[["MuF_Up","MuR_Up","MuF_MuR_Up","MuF_Down","MuR_Down","MuF_MuR_Down"]],[["TT_powheg_FSRdown","TT_powheg_FSRup"],["TT_powheg_ISRdown","TT_powheg_ISRup"]],[],[["Ttcc_Fraction_Up","Ttcc_Fraction_Down"]],[],[["Trig_Up","Trig_Down"]]]
pdfsys = []
for pdf in mceventweight:
   pdfsys.append(pdf["name"])
sysGroupList[-4].append(pdfsys)
print sysGroupList

visible="(NJets20>=4 && NbJets20>=2 && lepton1_pt>20 && lepton2_pt>20 && abs(lepton1_eta)<2.4 && abs(lepton2_eta)<2.4)"
ttbb = mAND("(NbJets20>=4)",visible)

full ="NaddJets20 >= 2"
ttbb_full = "NaddbJets20 >= 2"

rootdir = "/xrootd/store/user/chanwook/ttbb/"
TTList = ['TT_powheg']
weight = "weight*puweight"
nomtreepath = "cattree/nom"
mcXsecList = []

if recalculateMC:
  for i,ttName in enumerate(TTList):
      fileName = rootdir + TTList[i]+'.root'
      ntot = getWeightedEntries(fileName, nomtreepath, "filtered", weight)
      nttjjV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(visible,weight))
      nttbbV = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttbb,weight))
      nttjjF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(full,weight))
      nttbbF = getWeightedEntries(fileName, nomtreepath, "filtered", weightCut(ttbb_full,weight))

      ntotErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weight)
      nttjjVErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(visible,weight))
      nttbbVErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(ttbb,weight))
      nttjjFErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(full,weight))
      nttbbFErr = getWeightedEntriesError(fileName, nomtreepath, "filtered", weightCut(ttbb_full,weight))


      rVMC = nttbbV/nttjjV
      rVMCErr = bOverAErr(nttbbV,nttbbVErr,nttjjV,nttjjVErr)
      xttjjVMC = ttXsec*nttjjV/ntot
      #xttjjVMCErr = ttXsec*bOverAErr(nttjjV,nttjjVErr,ntot,ntotErr)
      xttjjVMCErr = abOverCErr(ttXsec,ttXsecErr,nttjjV,nttjjVErr,ntot,ntotErr) 
      xttbbVMC = ttXsec*nttbbV/ntot
      #xttbbVMCErr = ttXsec*bOverAErr(nttbbV,nttbbVErr,ntot,ntotErr)
      xttbbVMCErr = abOverCErr(ttXsec,ttXsecErr,nttbbV,nttbbVErr,ntot,ntotErr) 

      rFMC = nttbbF/nttjjF
      rFMCErr = bOverAErr(nttbbF,nttbbFErr,nttjjF,nttjjFErr)
      xttjjFMC = ttXsec*nttjjF/ntot
      #xttjjFMCErr = ttXsec*bOverAErr(nttjjF,nttjjFErr,ntot,ntotErr)
      xttjjFMCErr = abOverCErr(ttXsec,ttXsecErr,nttjjF,nttjjFErr,ntot,ntotErr) 
      xttbbFMC = ttXsec*nttbbF/ntot
      #xttbbFMCErr = ttXsec*bOverAErr(nttbbF,nttbbFErr,ntot,ntotErr)
      xttbbFMCErr = abOverCErr(ttXsec,ttXsecErr,nttbbF,nttbbFErr,ntot,ntotErr) 

      mcXsec = [rVMC, rVMCErr, xttjjVMC, xttjjVMCErr, xttbbVMC, xttbbVMCErr, rFMC, rFMCErr, xttjjFMC, xttjjFMCErr, xttbbFMC, xttbbFMCErr]
      out = {'name':ttName, 'values':mcXsec}
      mcXsecList.append(out)
  file = open(mcFile, "w")
  file.write(json.dumps(mcXsecList, indent=4))
  file.close()
else:
  file = open(mcFile).read()
  mcXsecList = json.loads(file)
mc = findDataSet("TT_powheg",mcXsecList)["values"]
f = open(resultFile,"r")
while True:
  line = f.readline()
  if not line: break
  v = line.split(" ")
  name = str(v[0])
  if name=="csvweight": 
    r = float(v[1])
    rErr = float(v[2])
    k = float(v[3])
    kErr = float(v[4])
    sys = findDataSet(name,entries)
    rEff = sys["nttjjR"]/sys["nttbbR"]*sys["nttbbV"]/sys["nttjjV"]
    rAccEff = sys["nttjjR"]/sys["nttbbR"]*sys["nttbbF"]/sys["nttjjF"]
    rV = rEff*r
    rF = rAccEff*r 
    rVErr = rEff*rErr
    rFErr = rAccEff*rErr

    rk = r*k
    xttjjV = sys["nttjjV"]*k/lumi
    xttjjVErr = sys["nttjjV"]*kErr/lumi
    xttbbV = sys["nttjjR"]*k*r/lumi/(sys["nttbbR"]/sys["nttbbV"])
    xttbbVErr = sys["nttjjR"]*pow(r*r*kErr*kErr+k*k*rErr*rErr,0.5)/lumi/(sys["nttbbR"]/sys["nttbbV"])

    xttjjF = sys["nttjjF"]*k/lumi
    xttjjFErr = sys["nttjjF"]*kErr/lumi
    xttbbF = sys["nttjjR"]*k*r/lumi/(sys["nttbbR"]/sys["nttbbF"])
    xttbbFErr = sys["nttjjR"]*pow(r*r*kErr*kErr+k*k*rErr*rErr,0.5)/lumi/(sys["nttbbR"]/sys["nttbbF"])

    cttjjR = sys["nttjjR"] 
    cttbbR = sys["nttbbR"]
    cttjjV = sys["nttjjV"]
    cttbbV = sys["nttbbV"]
    cttjjF = sys["nttjjF"]
    cttbbF = sys["nttbbF"]
    print "Process & Acceptance ($\mathcal{A}$) & Efficiency ($\epsilon$) \\\\\hline"
    print "$t\overline{t}b\overline{b}$ &", round(sys["nttbbV"]/sys["nttbbF"]*100,1), "\\% &", round(sys["nttbbR"]/sys["nttbbV"]*100), "\\%","\\\\\hline"
    print "$t\overline{t}j\overline{j}$ &", round(sys["nttjjV"]/sys["nttjjF"]*100,1), "\\% &", round(sys["nttjjR"]/sys["nttjjV"]*100), "\\%","\\\\\hline"
print r, " ", rErr, " ", rV, " ", rF
print ""

for p in ["V","F"]:
  totSysCate = [0.,0.,0.]
  for i,n in enumerate(sysGroupList):
    sysCate = [0.,0.,0.]
    for j in range(3): 
      for k,m in enumerate(n):
        tempList = []
        for l,ud in enumerate(m):
          if sysNameList[i]=="PDF": f = open(resultPDFFile,"r") 
          else : f = open(resultFile,"r")
          while True:
            line = f.readline()
            if not line: break
            v = line.split(" ")
            name = str(v[0])
            sysValue = float(v[1])
            kValue = float(v[3])
            if sysNameList[i]=="PDF": sys = findDataSet(v[0],pdfentries)
            else : sys = findDataSet(v[0],entries)
            if not sys: continue
            if ud==name:
              
              if p=="V" and j==0:
                val = sys["nttjjR"]*kValue*sysValue/lumi/(sys["nttbbR"]/sys["nttbbV"])
                V = xttbbV
              elif p=="V" and j==1:
                val = sys["nttjjV"]*kValue/lumi
                V = xttjjV 
              elif p=="V" and j==2:
                val = sys["nttjjR"]/sys["nttbbR"]*sys["nttbbV"]/sys["nttjjV"]*sysValue
                V = rV
              elif p=="F" and j==0:
                val = sys["nttjjR"]*kValue*sysValue/lumi/(sys["nttbbR"]/sys["nttbbF"])
                #val = cttjjR/(cttbbR/cttbbF)/lumi*kValue*sysValue
                V = xttbbF 
              elif p=="F" and j==1:
                val = sys["nttjjF"]*kValue/lumi
                #val = cttjjF*kValue/lumi
                V = xttjjF
              elif p=="F" and j==2:
                val = sys["nttjjR"]/sys["nttbbR"]*sys["nttbbF"]/sys["nttjjF"]*sysValue
                #val = cttjjR/cttbbR*cttbbF/cttjjF*sysValue
                V = rF
              if sysNameList[i]=="PDF": tempList.append(V-val)
              else : tempList.append((V-val)**2)
        if not tempList: continue
        if sysNameList[i]=="PDF": subSys = (sum(tempList)/len(tempList))**2
        else : subSys = max(tempList)  
        #if j==2:print ud, pow(subSys, 0.5)*1000
        sysCate[j] += subSys     
    for c in range(3):
      totSysCate[c] += sysCate[c]
    if p=="V":
      if n: print sysNameList[i],"& ",round(abs(pow(sysCate[0],0.5)/xttbbV*100.),1), "\\% & ", round(abs(pow(sysCate[1],0.5)/xttjjV*100.),1), "\\% & ",round(abs(pow(sysCate[2],0.5)/rV*100.),1), "\\%","\\\\\hline"
      else: print sysNameList[i],"& Will be updated & Will be updated & Will be updated \\\\\hline" 
    elif p=="F":
      if n: print sysNameList[i],"& ",round(abs(pow(sysCate[0],0.5)/xttbbF*100.),1), "\\% & ", round(abs(pow(sysCate[1],0.5)/xttjjF*100.),1), "\\% & ",round(abs(pow(sysCate[2],0.5)/rF*100.),1), "\\%","\\\\\hline" 
      else: print sysNameList[i],"& Will be updated & Will be updated & Will be updated \\\\\hline" 
  if p=="V": print "Total Uncertainty & ",round(abs(pow(totSysCate[0],0.5)/xttbbV*100.),1), "\\% & ", round(abs(pow(totSysCate[1],0.5)/xttjjV*100.),1), "\\% & ",round(abs(pow(totSysCate[2],0.5)/rV*100.),1), "\\%","\\\\\hline" 
  if p=="F": print "Total Uncertainty & ",round(abs(pow(totSysCate[0],0.5)/xttbbF*100.),1), "\\% & ", round(abs(pow(totSysCate[1],0.5)/xttjjF*100.),1), "\\% & ",round(abs(pow(totSysCate[2],0.5)/rF*100.),1), "\\%","\\\\\hline" 
  print ""
  print "" 
  if p=="V":
    obs = [xttbbV, xttbbVErr, pow(xttbbVErr**2 + totSysCate[0], 0.5), xttjjV, xttjjVErr, pow(xttbbVErr**2 + totSysCate[1], 0.5), 100*rV, 100*rVErr, 100*pow(rVErr**2 + totSysCate[2], 0.5)]
    exp = [mc[4], mc[5], mc[2], mc[3], 100*mc[0], 100*mc[1]]
    print "Measurement & $", round(xttbbV,3), "\\pm", round(pow(xttbbVErr**2+totSysCate[0],0.5),3), " $ & $", round(xttjjV,2), "\\pm", round(pow(xttjjVErr**2+totSysCate[1],0.5),2), "$ & $", round(rV,4), "\\pm", round(pow(rVErr**2+totSysCate[2],0.5),4), "$ \\\\"
    print " "
    print "Measurement & $", round(xttbbV,3), "\\pm", round(xttbbVErr,3), "\\pm", round(pow(totSysCate[0],0.5),3), "$ & $", round(xttjjV,2), "\\pm", round(xttjjVErr,2), "\\pm", round(pow(totSysCate[1],0.5),2), "$ & $", round(rV,3), "\\pm", round(rVErr,3), "\\pm", round(pow(totSysCate[2],0.5),3), "$ \\\\"
    print "SM (POWHEG) & $", round(mc[4],3), "\\pm", round(mc[5],3), "$ & $", round(mc[2],2), "\\pm", round(mc[3],2), "$ & $", round(mc[0],3), "\\pm", round(mc[1],4), "$ \\\\\hline" 
    print " "  

    print "Measurement ", round(xttbbV,3), "+-", round(pow(xttbbVErr**2+totSysCate[0],0.5),3), ",  ", round(xttjjV,2), "+-", round(pow(xttjjVErr**2+totSysCate[1],0.5),2), ",  ", round(rV,3), "+-", round(pow(rVErr**2+totSysCate[2],0.5),3)
    print " "
    print "Measurement ", round(xttbbV,3), "+-", round(xttbbVErr,3), "+-", round(pow(totSysCate[0],0.5),3), ",  ", round(xttjjV,2), "+-", round(xttjjVErr,2), "+-", round(pow(totSysCate[1],0.5),2), ",  ", round(rV,4), "+-", round(rVErr,4), "+-", round(pow(totSysCate[2],0.5),4)
    print "SM (POWHEG) ", round(mc[4],3), "+-", round(mc[5],3), ",  ", round(mc[2],2), "+-", round(mc[3],2), ",  ", round(mc[0],3), "+-", round(mc[1],4) 
    print " "  

  elif p=="F":
    obs = [xttbbF, xttbbFErr, pow(xttbbFErr**2 + totSysCate[0], 0.5), xttjjF, xttjjFErr, pow(xttbbFErr**2 + totSysCate[1], 0.5), 100*rF, 100*rFErr, 100*pow(rFErr**2 + totSysCate[2], 0.5)]
    exp = [mc[10], mc[11], mc[8], mc[9], 100*mc[6], 100*mc[7]]
    print "Measurement & $", round(xttbbF,2), "\\pm", round(pow(xttbbFErr**2+totSysCate[0],0.5),2),"$ & $", int(round(xttjjF)), "\\pm", int(round(pow(xttjjFErr**2+totSysCate[1],0.5))),"$ & $", round(rF,3), "\\pm", round(pow(rFErr**2+totSysCate[2],0.5),3), "$ \\\\"
    print " "
    print "Measurement & $", round(xttbbF,2), "\\pm", round(xttbbFErr,2), "\\pm", round(pow(totSysCate[0],0.5),2), "$ & $", int(round(xttjjF)), "\\pm", int(round(xttjjFErr)), "\\pm", int(round(pow(totSysCate[1],0.5))), "$ & $", round(rF,3), "\\pm", round(rFErr,3), "\\pm", round(pow(totSysCate[2],0.5),3), "$ \\\\"
    print "SM (POWHEG) & $", round(mc[10],2), "\\pm", round(mc[11],2), "$ & $", int(round(mc[8])), "\\pm", int(round(mc[9])), "$ & $", round(mc[6],3), "\\pm", round(mc[7],4), "$ \\\\\hline"  
    print " "

    print "Measurement ", round(xttbbF,2), "+-", round(pow(xttbbFErr**2+totSysCate[0],0.5),2),",  ", int(round(xttjjF)), "+-", int(round(pow(xttjjFErr**2+totSysCate[1],0.5))),",  ", round(rF,3), "+-", round(pow(rFErr**2+totSysCate[2],0.5),3)
    print " "
    print "Measurement ", round(xttbbF,2), "+-", round(xttbbFErr,2), "+-", round(pow(totSysCate[0],0.5),2), ",  ", int(round(xttjjF)), "+-", int(round(xttjjFErr)), "+-", int(round(pow(totSysCate[1],0.5))), ",  ", round(rF,4), "+-", round(rFErr,4), "+-", round(pow(totSysCate[2],0.5),4)
    print "SM (POWHEG) ", round(mc[10],2), "+-", round(mc[11],2), ",  ", int(round(mc[8])), "+-", int(round(mc[9])), ",  ", round(mc[6],3), "+-", round(mc[7],4)
  print ""
  print ""

  c = TCanvas(p,p,1000,500)
  name = ["t#bar{t}b#bar{b} #sigma (pb)","t#bar{t}jj #sigma (pb)","t#bar{t}b#bar{b}/t#bar{t}jj (%)"]
  objList = []
  for i in range(3):
    pad = TPad(name[i],name[i],0.25*i,0,0.25*(i+1),1)
    pad.Draw()
    pad.cd()

    gStatSys = TGraphErrors()#"gStatSys"+p,"gStatSys"+p)
    gStatSys.SetName(name[i]+"StatSys"+p)
    gStatSys.SetPoint(0,obs[3*i],1)
    gStatSys.SetPointError(0,obs[3*i+2],1)
    gStatSys.SetFillColor(2)
    gStatSys.SetFillStyle(3001)
    gStatSys.GetXaxis().SetLimits(obs[3*i]-obs[3*i+2]*2, obs[3*i]+obs[3*i+2]*3)
    gStatSys.GetXaxis().SetLabelSize(0.035)
    gStatSys.GetXaxis().SetTitle(name[i])
    gStatSys.GetYaxis().SetTickSize(0)
    gStatSys.GetYaxis().SetTickLength(0)
    gStatSys.GetYaxis().SetLabelSize(0)
    gStatSys.Draw("a2")

    gStat = TGraphErrors()#"gStat"+p,"gStat"+p)
    gStat.SetName(name[i]+"Stat"+p)
    gStat.SetPoint(0,obs[3*i],1)
    gStat.SetPointError(0,obs[3*i+1],1)
    gStat.SetFillColor(3)
    gStat.SetFillStyle(3001)
    gStat.Draw("2")

    gExp = TGraphErrors()
    gExp.SetName(name[i]+"Exp"+p)
    gExp.SetPoint(0,exp[2*i],1)
    #gExp.SetPointError(0,exp[2*i+1],0)
    gExp.SetMarkerStyle(21)
    gExp.SetMarkerSize(1)
    gExp.SetMarkerColor(kBlue)
    gExp.Draw("p")  

    l = TLine(obs[3*i],0,obs[3*i],2)
    #l.SetName(name[i]+p)
    l.SetLineStyle(2)
    l.SetLineWidth(2)
    l.Draw()

    objList.extend([gStatSys, gStat, gExp,l])
    pad.Update()
    c.cd()
  
  padInfo = TPad("3","3",0.75,0,1,1)
  padInfo.Draw()
  padInfo.cd()

  la1 = TLatex()
  la1.SetTextSize(0.1)
  la1.SetTextFont(61)
  la1.DrawLatex(.1,.9,"CMS")

  la2 = TLatex()
  la2.SetTextSize(0.075)
  la2.SetTextFont(52)
  la2.DrawLatex(.35,.9,"Preliminary")
  la2.DrawLatex(.1,.8,lumiText)

  leg = TLegend(0.1, 0.5, 0.9, 0.75)
  leg.SetTextSize(0.075)
  leg.SetFillStyle(0)
  leg.SetBorderSize(0)
  leg.AddEntry(l,"Measurement","l")
  leg.AddEntry(gStat,"Stat.","f")
  leg.AddEntry(gStatSys,"Stat. #oplus Sys.","f")
  leg.AddEntry(gExp, "POWHEG", "p")
  leg.Draw()
  c.cd() 
  c.Update() 
  c.Draw()
  c.SaveAs("obsVsExp%s_%s.png"%(p,nbjets))
#print "\\item ${\\sigma}_{t\\overline{t}b\\overline{b}}/{\\sigma}_{t\\overline{t}j\\overline{j}}(Full)=",round(rF,4), "\\pm", round(rErr*rAccEff,4), "(stat.)", "\\pm", round(totSysF,4),"(syst.)$"
#print texName,"&" " ", round(abs((r-float(v[1])*accAndEff)/r*100.),1),"\\%","\\\\\hline"
