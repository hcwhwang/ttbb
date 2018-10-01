from ROOT import *
import copy

f = TFile.Open("csvNbjets2VisNestPDF_nuisance.root","UPDATE")

ChList = ["MuEl", "ElEl", "MuMu"]
PdfList = ["ttbb","ttbj","ttcclf","ttcclf2","ttothers","bkg"]
#VaryList = ["MuRUp","MuFUp","MuRMuFUp","MuRDown","MuFDown","MuRMuFDown"]
sysName = "PDF"
for i,ch in enumerate(ChList):
  #if i!=0: continue
  for j,pdf in enumerate(PdfList):
    #if j!=0: continue
    hUpName = ch + "_" + pdf + "_PDFUp"
    hDownName = ch + "_" + pdf + "_PDFDown"
    hPDFUp = TH1F(hUpName, hUpName, 100, 0, 100)
    hPDFDown = TH1F(hDownName, hDownName, 100, 0, 100)
    for b in range(100):
      #if b!=4: continue
      binList = []
      errList = []
      for k in range(100):
        hVary = f.Get(ch + "_" + pdf + "_pdf" + str(k))
        print hVary.GetBinContent(b+1), hVary.GetBinError(b+1)
        binList.append(hVary.GetBinContent(b+1))
        errList.append(hVary.GetBinError(b+1))
      mean = sum(binList)/float(len(binList)) 
      if mean==0.: 
        hPDFUp.SetBinContent(b+1,0)
        hPDFUp.SetBinError(b+1,0) 
        hPDFDown.SetBinContent(b+1,0)
        hPDFDown.SetBinError(b+1,0) 
      else:
        rms = 0.
        for x, xx in enumerate(binList):
          rms += (mean-xx)*(mean-xx)
        rms = rms/float(len(binList))
        rms = pow(rms,0.5)
        binErrUp, binErrDown = 0., 0.
        binErrPart1, binErrPart2 = 0., 0.
        hPDFUp.SetBinContent(b+1,mean+rms)  
        for v, value in enumerate(binList):
          binErrUp += pow(errList[v]*(rms+value-mean)/float(len(binList))/rms,2)
          binErrDown += pow(errList[v]*(rms-value+mean)/float(len(binList))/rms,2)
          if binErrDown<0: binErrDown=0.
        binErrUp = pow(binErrUp, 0.5)
        binErrDown = pow(binErrDown, 0.5)
        hPDFUp.SetBinError(b+1, binErrUp)
        print mean, rms, binErrUp, binErrDown
        if mean>rms: 
          hPDFDown.SetBinContent(b+1, mean-rms)
          hPDFDown.SetBinError(b+1, binErrDown)
        else:
          hPDFDown.SetBinContent(b+1,0)
          hPDFDown.SetBinError(b+1,0)
    f1 = TFile.Open("csvNbjets2VisNestPDFIndi1_nuisance.root","UPDATE")
    hPDFUp.Write()
    hPDFDown.Write()
    f1.Close()
f.Close()


