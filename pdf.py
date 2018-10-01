from ROOT import *
import copy

f = TFile.Open("csvNbjets2VisNestPDF_nuisance.root","UPDATE")

ChList = ["MuEl", "ElEl", "MuMu"]
PdfList = ["ttbb","ttbj","ttcclf","ttcclf2","ttothers","bkg"]
#VaryList = ["MuRUp","MuFUp","MuRMuFUp","MuRDown","MuFDown","MuRMuFDown"]
sysName = "PDF"
for i,ch in enumerate(ChList):
  for j,pdf in enumerate(PdfList):
    for k in range(102):
      hVary = f.Get(ch + "_" + pdf + "_pdf" + str(k))

      hPDF = copy.deepcopy(hVary)
      if k%2==0: hName = ch + "_" + pdf + "_PDF" + str(((k+1)//2)+1) + "Up"
      elif k%2==1: hName = ch + "_" + pdf + "_PDF" + str(((k+1)//2)+1) + "Down" 
      print k, hName
      hPDF.SetName(hName)
      hPDF.SetTitle(hName)

      f1 = TFile.Open("csvNbjets2VisNestPDFIndi1_nuisance.root","UPDATE")
      hPDF.Write()
      f1.Close()
f.Close()


