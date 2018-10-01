from ROOT import *
import copy

f = TFile.Open("csvNbjets2VisNestPDF_nuisance.root","UPDATE")

ChList = ["MuEl", "ElEl", "MuMu"]
PdfList = ["ttbb","ttbj","ttcclf","ttcclf2","ttothers","bkg"]
#VaryList = ["MuRUp","MuFUp","MuRMuFUp","MuRDown","MuFDown","MuRMuFDown"]
sysName = "PDF"
for i,ch in enumerate(ChList):
  for j,pdf in enumerate(PdfList):
    hCen = f.Get(ch + "_" + pdf)
    #hName = ch + "_" + pdf + "_" + sysName
    #hUp = TH1F(hName+"Up", hName+"Up",100,0,100)
    #hDown = TH1F(hName+"Down", hName+"Down",100,0,100)
    for b in range(1,101):
      print ch, pdf, b
      upValue, downValue = 0, 9999999999
      upValueErr, downValueErr = 0, 0
      bb = 10
      '''
      if hCen.GetBinContent(b)==0:
        hUp.SetBinContent(b,0.)
        hUp.SetBinError(b,0.)
        hDown.SetBinContent(b,0.)
        hDown.SetBinError(b,0.)
        continue
      '''
      #for k,vary in enumerate(VaryList):
      for k in range(102):
        if b!=1: continue
        #hVary = f.Get(ch + "_" + pdf + "_" + vary)
        hVary = f.Get(ch + "_" + pdf + "_pdf" + str(k))
        hName = ch + "_" + pdf + "_pdf" + str(k)#case2

        hUp = copy.deepcopy(hVary)#case2
        hDown = copy.deepcopy(hVary)#case2
        hUp.SetName(hName+"Up")
        hUp.SetTitle(hName+"Up")
        hDown.SetName(hName+"Down")
        hDown.SetTitle(hName+"Down")

        f1 = TFile.Open("csvNbjets2VisNestPDFIndi_nuisance.root","UPDATE")
        hUp.Write()#case2
        hDown.Write()#case2
        f1.Close()
        '''
        tempValue = hVary.GetBinContent(b)
        tempValueErr = hVary.GetBinError(b)
        if tempValue > upValue:
          upValue = tempValue
          upValueErr = tempValueErr
          bb = k
        if tempValue < downValue:
          downValue = tempValue
          downValueErr = tempValueErr
      hUp.SetBinContent(b, upValue)
      hUp.SetBinError(b, upValueErr)
      hDown.SetBinContent(b, downValue)
      hDown.SetBinError(b, downValueErr)
    hUp.Write()
    hDown.Write()
    '''
f.Close()

'''
c1 = TCanvas("c1","c1",500,500)
hUp.Draw()
c2 = TCanvas("c2","c2",500,500)
hDown.Draw()
c3 = TCanvas("c3","c3",500,500)
hCen.Draw()
'''
'''
for i in range(1,101):
  print hUp.GetBinContent(i)
  print hCen.GetBinContent(i)
  print hDown.GetBinContent(i)
  print ""
'''
