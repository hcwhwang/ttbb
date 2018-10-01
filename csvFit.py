from ROOT import *
import sys, json, os
from sysWeight_cfi import *
from array import array
gROOT.SetBatch(True)
def addLegendCMS():
  #tex2 = TLatex(0.3715952,0.9146667,"Preliminary")
  tex2 = TLatex(-20.,50.,"Preliminary")
  tex2.SetNDC()
  tex2.SetTextAlign(12)
  tex2.SetX(0.25)
  tex2.SetY(0.97)
  tex2.SetTextColor(2)
  tex2.SetTextFont(42)
  tex2.SetTextSize(0.05)
  tex2.SetTextSizePixels(24)
  #tex2.Draw()

  return tex2

def make_legend(xmin,ymin,xmax,ymax):
  #leg = TLegend(0.65,0.7, 0.89,0.89)
  leg = TLegend(xmin,ymin,xmax,ymax)
  leg.SetFillColor(0)
  leg.SetLineColor(1)
  leg.SetTextFont(62)
  leg.SetTextSize(0.03)

  leg.SetBorderSize(1)
  leg.SetLineStyle(1)
  leg.SetLineWidth(1)
  leg.SetLineColor(0)

  return leg
name = "Nbjets2"
result = "fitresult%s.txt" % name
rootFile = "csv%s.root" % name
if os.path.exists(result): os.system("rm -f %s" % result)
f = TFile.Open(rootFile)

for i,d in enumerate(mceventweight):
  #if i!=0: continue
  mcname = mceventweight[i]["name"]
  print mcname 

  if mcname not in ['TT_powheg_FSRdown', 'TT_powheg_FSRup', 'TT_powheg_ISRdown', 'TT_powheg_ISRup', 'TT_powheg_UEdown', 'TT_powheg_UEup', 'TT_powheg_herwig', 'TTJets_aMC']:
    httbb = f.Get("TT_powheg_ttbb_%s_%s" % (mceventweight[i]["tree"],mceventweight[i]["name"]))
    httbj =  f.Get("TT_powheg_ttbj_%s_%s" % (mceventweight[i]["tree"],mceventweight[i]["name"]))
    httcc = f.Get("TT_powheg_ttcc_%s_%s" % (mceventweight[i]["tree"],mceventweight[i]["name"]))
    httlf = f.Get("TT_powheg_ttlf_%s_%s" % (mceventweight[i]["tree"],mceventweight[i]["name"]))
    httot = f.Get("TT_powheg_ttothers_%s_%s" % (mceventweight[i]["tree"],mceventweight[i]["name"]))
  else:
    httbb = f.Get("%s_ttbb_%s_%s" % (mcname, mceventweight[i]["tree"],mcname))
    httbj =  f.Get("%s_ttbj_%s_%s" % (mcname, mceventweight[i]["tree"],mcname))
    httcc = f.Get("%s_ttcc_%s_%s" % (mcname, mceventweight[i]["tree"],mcname))
    httlf = f.Get("%s_ttlf_%s_%s" % (mcname, mceventweight[i]["tree"],mcname))
    httot = f.Get("%s_ttothers_%s_%s" % (mcname, mceventweight[i]["tree"],mcname))
  hbkg = f.Get("Bkg_%s_%s" % (mceventweight[i]["tree"],mceventweight[i]["name"]))

  httcclf = httcc.Clone()
  httcclf.Add(httlf)
  hTTList = [httbb, httbj, httcc, httlf, httot, hbkg]
  hrd = f.Get("Data_nom2_csvweight")

  nList = []
  for j,h in enumerate(hTTList):
    nList.append(h.Integral())

  n_ttbb = nList[0]
  n_ttbj = nList[1]
  n_ttcc = nList[2]
  n_ttlf = nList[3]
  n_ttot = nList[4]
  n_bkg  = nList[5]

  n_ttjj = n_ttbb+n_ttbj+n_ttcc+n_ttlf
  n_ttbar = n_ttjj+n_ttot

  rttbb = n_ttbb/n_ttjj
  print "init r = ", rttbb
  rttbj  = n_ttbj/n_ttjj
  rttcc = (n_ttcc)/n_ttjj

  x = RooRealVar("x","x",0,1)
  y = RooRealVar("y","y",0,1)

  RttbbReco=RooRealVar("RttbbReco","RttbbReco",rttbb,rttbb,rttbb);
  RttbjReco =RooRealVar("RttbjReco", "RttbjReco", rttbj, rttbj, rttbj);
  RttccReco=RooRealVar("RttccReco","RttccReco",rttcc,rttcc,rttcc);

  #fttbb   =RooRealVar(    "fttbb",                "fttbb",           rttbb, 0.0, 0.1)
  fttbb   =RooRealVar(    "fttbb",                "fttbb",           rttbb, 0.0, 0.9)
  fttbjcon  =RooFormulaVar("fttbjcon",          "fttbj","@0/@1*@2",RooArgList(fttbb,RttbbReco,RttbjReco) )  # constraint fttb with fttbb
  fttbj  =RooRealVar(   "fttbj",                "fttbj",          rttbj, 0.0, 0.9)  # free fttbj
  fttcc =RooRealVar(  "fttcc",              "fttcc",           rttcc, 0.0, 0.9)  # free fttbbcc
  k      =RooRealVar(       "k","normalization factor",           1, 0.7, 1.3)
  nttjj =RooRealVar(    "nttjj","number of ttjj events",                            n_ttjj , n_ttjj, n_ttjj)
  knttjj=RooFormulaVar("knttjj","number of ttjj events after fitting","k*nttjj",    RooArgList(k,nttjj) )
  nttot =RooRealVar(    "nttot","number of ttot events",                            n_ttot , n_ttot, n_ttot)
  knttot=RooFormulaVar("knttot","number of ttot events after fitting","k*nttot",    RooArgList(k,nttot) )
  nbkg  =RooRealVar(     "nbkg","number of background events",                      n_bkg , n_bkg, n_bkg)

  nttcc =RooRealVar(   "nttcc","number of ttcc events",                         n_ttcc , n_ttcc, n_ttcc)
  knttcc=RooFormulaVar("knttcc","number of ttcc events after fitting","k*nttcc",RooArgList(k,nttcc) )

  xArg = RooArgList(x, y)
  data    = RooDataHist("data",    "data set with (x)",   xArg, hrd)
  ttbb    = RooDataHist("ttbb",    "ttbb set with (x)",   xArg, httbb)
  ttbj    = RooDataHist("ttbj",    "ttbj set with (x)",   xArg, httbj)
  ttcc    = RooDataHist("ttcc",    "ttcc set with (x)",   xArg, httcc)
  ttlf    = RooDataHist("ttlf",    "ttlf set with (x)",   xArg, httlf)
  ttcclf  = RooDataHist("ttcclf",  "ttcclf set with (x)", xArg, httcclf)
  ttot    = RooDataHist("ttot",    "ttot set with (x)",   xArg, httot)
  bkg     = RooDataHist("bkg",     "bkg set with (x)",    xArg, hbkg)

  ttbbpdf      = RooHistPdf("ttbbpdf",     "ttbbpdf",      RooArgSet(RooArgList(x,y)), ttbb)
  ttbjpdf      = RooHistPdf("ttbjpdf",     "ttbjpdf",      RooArgSet(RooArgList(x,y)), ttbj)
  ttccpdf      = RooHistPdf("ttccpdf",     "ttccpdf",      RooArgSet(RooArgList(x,y)), ttcc)
  ttlfpdf      = RooHistPdf("ttlfpdf",     "ttlfpdf",      RooArgSet(RooArgList(x,y)), ttlf)
  ttcclfpdf    = RooHistPdf("ttcclfpdf",   "ttcclfpdf",    RooArgSet(RooArgList(x,y)), ttcclf)
  ttotpdf      = RooHistPdf("ttotpdf",     "ttotpdf",      RooArgSet(RooArgList(x,y)), ttot)
  bkgpdf      = RooHistPdf("bkgpdf",     "bkgpdf",      RooArgSet(RooArgList(x,y)), bkg)

  #freeTTB = False
  
  freeTTB = False
  freeTTCC = False
  if mceventweight[i]["name"]=="Free_ttbj":
    freeTTB = True
    freeTTCC = False 
  if mceventweight[i]["name"]=="Free_ttcc":
    freeTTB = False
    freeTTCC = True 
  if freeTTB and not freeTTCC  : model  = RooAddPdf("model",   "model",RooArgList( ttbbpdf, ttbjpdf, ttcclfpdf), RooArgList(fttbb,fttbj))
  #if freeTTB and not freeTTCC  : model  = RooAddPdf("model",   "model",RooArgList( ttbbpdf, ttbjpdf, ttcclfpdf), RooArgList(fttbb,fttbj))
  elif not freeTTB and freeTTCC: model  = RooAddPdf("model",   "model",RooArgList( ttbbpdf, ttbjpdf,  ttccpdf, ttlfpdf), RooArgList(fttbb,fttbjcon,fttcc))
  elif freeTTB and freeTTCC    : model  = RooAddPdf("model",   "model",RooArgList( ttbbpdf, ttbjpdf, ttccpdf, ttlfpdf), RooArgList(fttbb,fttbj, fttcc))
  else                         : model  = RooAddPdf("model",   "model",RooArgList( ttbbpdf, ttbjpdf,  ttcclfpdf), RooArgList(fttbb,fttbjcon))
  #else                         : model  = RooAddPdf("model",   "model",RooArgList( ttbbpdf, ttbjpdf,ttcclfpdf), RooArgList(fttbb,fttbjcon))

  model2 = RooAddPdf("model2", "model2",RooArgList( model, ttotpdf, bkgpdf),              RooArgList(knttjj,knttot,nbkg)) # k*bkg
  #model2 = RooAddPdf("model2", "model2",RooArgList( model),              RooArgList(knttjj)) # k*bkg
  #model2 = RooAddPdf("model2", "model2",RooArgList( model, ttotpdf, bkgpdf),              RooArgList(knttjj,knttot,nbkg)) # fixing bkg
  model2.fitTo(data)

  fttbb.Print()
  k.Print()

  print mceventweight[i]["name"]
  file = open(result,"a")

  if not freeTTB :line = "%s %f %f %f %f %f %f %f\n" % (d["name"], fttbb.getVal(), fttbb.getError(), k.getVal(), k.getError(), fttbjcon.getVal(), rttbb, rttbj)
  else :line = "%s %f %f %f %f %f %f %f %f\n" % (d["name"], fttbb.getVal(), fttbb.getError(), k.getVal(), k.getError(), fttbj.getVal(), rttbb, rttbj, fttbj.getError())
  file.write(line)
  file.close()
  #continue
  if mceventweight[i]["name"] not in ["csvweight","Free_ttbj","HF_Down"]: continue
  cR10 = TCanvas("R10", "R", 1)#500, 500)

  nll = model2.createNLL(data)
  #nll = model2.createNLL(ttlf)
  #RooMinuit(nk1).migrad() 
  RFrame = fttbb.frame()
  nll.plotOn(RFrame,RooFit.ShiftToZero())
  RFrame.SetMaximum(4.);RFrame.SetMinimum(0)
  RFrame.GetXaxis().SetTitle("R_{Reco}=ttbb/ttjj")
  RFrame.SetTitle("")
  RFrame.Draw()
  pt = addLegendCMS()

  lineKKK = TLine(0,0,0,0)
  lineKKK.SetLineColor(kBlue)
  lineKKK.SetLineWidth(3)

  line = TLine(RFrame.GetXaxis().GetXmin() ,0.5,RFrame.GetXaxis().GetXmax(),0.5)
  line.SetLineColor(kRed)
  line.Draw()

  lineTbb = TLine(rttbb,RFrame.GetMaximum(),rttbb,0)
  lineTbb.SetLineStyle(2)
  lineTbb.Draw()

  l1 = make_legend(0.49,0.76,0.93,0.88)
  l1.AddEntry(lineTbb,"prefit: R="+str(round(rttbb*10000)/10000),"l")
  recoR = fttbb.getVal()
  recoRerror = fttbb.getError()
  recoRP = fttbj.getVal()
  l1.AddEntry(lineKKK,"fit: R="+str(round(recoR*10000)/10000)+" #pm "+str(round(recoRerror*10000)/10000)+"","l")
  l1.SetTextSize(0.04)
  l1.SetFillColor(0)
  l1.SetLineColor(0)
  l1.Draw()

  pt.Draw()
  cR10.SaveAs("ControlPlots/fttbbNLL_%s_%s.png"%(mceventweight[i]["name"],name))

  cR00 = TCanvas("R00", "R", 1)#500, 500)

  nll = model2.createNLL(data)
  #nll = model2.createNLL(ttlf)
  #RooMinuit(nk1).migrad() 
  RFrameK = k.frame()
  nll.plotOn(RFrameK,RooFit.ShiftToZero())
  RFrameK.SetMaximum(4.);RFrameK.SetMinimum(0)
  RFrameK.GetXaxis().SetTitle("K")
  RFrameK.SetTitle("")
  RFrameK.Draw()


  lineK = TLine(RFrameK.GetXaxis().GetXmin() ,0.5,RFrameK.GetXaxis().GetXmax(),0.5)
  lineK.SetLineColor(kRed)
  lineK.Draw()

  lineTbbK = TLine(1,RFrameK.GetMaximum(),1,0)
  lineTbbK.SetLineStyle(2)
  lineTbbK.Draw()

  l1K = make_legend(0.49,0.76,0.93,0.88)
  l1K.AddEntry(lineTbbK,"prefit: k=1.0","l")
  kVal = k.getVal()
  kValError = k.getError()
  l1K.AddEntry(lineKKK,"fit: k="+str(round(kVal*10000)/10000)+" #pm "+str(round(kValError*10000)/10000)+"","l")
  l1K.SetTextSize(0.04)
  l1K.SetFillColor(0)
  l1K.SetLineColor(0)
  l1K.Draw()

  pt.Draw()
  cR00.SaveAs("ControlPlots/kNLL_%s_%s.png"%(mceventweight[i]["name"],name))

  cR11 = TCanvas("R11", "R", 1)# 500, 500)
  xframe = x.frame()
  data.plotOn(xframe, RooFit.DataError(RooAbsData.SumW2) )
  #model2.paramOn(xframe, RooFit.Layout(0.65,0.9,0.9) )
  model2.plotOn(xframe)
  chi2 = xframe.chiSquare(2)
  ndof = xframe.GetNbinsX()
  print "chi2 = "+ str(chi2)
  print "ndof = "+ str(ndof)
  xframe.SetMaximum(xframe.GetMaximum()*1.5)
  xframe.Draw()

  pt.Draw()
  cR11.SaveAs("ControlPlots/csv4th_%s_%s.png"%(mceventweight[i]["name"],name))

  cR12 = TCanvas("R12", "R", 1)# 500, 500)
  yframe = y.frame()
  data.plotOn(yframe, RooFit.DataError(RooAbsData.SumW2) )
  #model2.paramOn(yframe, RooFit.Layout(0.65,0.9,0.9) )
  model2.plotOn(yframe)
  chi2 = yframe.chiSquare(2)
  ndof = yframe.GetNbinsX()
  print "chi2 = "+ str(chi2)
  print "ndof = "+ str(ndof)
  yframe.SetMaximum(yframe.GetMaximum()*1.5)
  yframe.Draw()

  pt.Draw()
  cR12.SaveAs("ControlPlots/csv3rd_%s_%s.png"%(mceventweight[i]["name"],name))
  if freeTTB:
    cR13 = TCanvas("R13", "R", 1)# 500, 500)
    nll22 = model2.createNLL(data)
    m=RooMinuit(nll22)
    frameNLLContour = m.contour(fttbb, fttbj,1,2,3)
    frameNLLContour.GetXaxis().SetTitle("R (ttbb/ttjj)")
    frameNLLContour.GetYaxis().SetTitle("R' (ttbj/ttjj)")
    frameNLLContour.SetMarkerStyle(21)
    frameNLLContour.Draw()
  
    preM = TMarker(rttbb,rttbj,20)
    preM.SetMarkerColor(kRed)
    preM.Draw()
    preM2 = TMarker(rttbb,rttbj,20)
    preM2.SetMarkerColor(kBlack)

    
    preM23 = TMarker(0.0421,0.1064,20)
    preM23.SetMarkerColor(kBlue)
    #preM23.Draw()

    pt.Draw()
  
    l2 = make_legend(0.38,0.63,0.94,0.90)
    l2.AddEntry(preM,"prefit: R="+str(round(rttbb,4)),"p")
    l2.AddEntry(preM,"prefit: R'="+str(round(rttbj,4)),"p")
    #l2.AddEntry(preM23,"fit(cen.): R="+str((recoR)) ,"p")
    #l2.AddEntry(preM23,"fit(cen.): R2="+str((recoRP)),"p")
  
    l2.AddEntry(preM2,"fit: R="+str(round(fttbb.getVal(),5))+" #pm "+str(round(fttbb.getError(),5))+"","p")
    l2.AddEntry(preM2,"fit: R'="+str(round(fttbj.getVal(),5))+" #pm "+str(round(fttbj.getError(),5))+"","p")
    l2.SetTextSize(0.04)
    l2.SetFillColor(0)
    l2.SetLineColor(0)
    l2.Draw()    
    cR13.SaveAs("ControlPlots/ttbbVsTtbj_%s_%s.png"%(mceventweight[i]["name"],name))
f.Close()
