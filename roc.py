from ROOT import *

f = TFile.Open("TMVA.root")

algoList =["Cuts", "Likelihood", "LD", "BDT","BDT","BDT"]
algoSubList =["CutsPCA", "Likelihood", "LD", "BDT","BDTG","BDTB"]
nameList = ["Cut Optimisation", "1-D Likelihood", "Linear Discriminant Analysis",  "Neural Network","K-Nearest Neighbor","Boosted Decision Tree"]
colList = [2,3,4,5,6,7,8]

c = TCanvas("c","c",500,500)
#gFrame = TGraph()
#gFrame.SetTitle(";Sensitivity;Specificit")
#gFrame.SetName("ROC")

leg = TLegend(0.5, 0.6, 0.9, 0.9)
leg.SetFillStyle(0)
leg.SetBorderSize(0)

gStyle.SetOptStat(0)
for i,algo in enumerate(algoList):
  #if i!=0: continue
  h = f.Get("Method_%s/%s/MVA_%s_rejBvsS" % (algo, algoSubList[i], algoSubList[i]))
  print nameList[i], h.Integral()
  #g = TGraph(algoSubList[i],";Sensitivity;Specificity")
  globals()[algoSubList[i]]= TGraph()
  globals()[algoSubList[i]].SetName(algoSubList[i])
  if i==0:globals()[algoSubList[i]].SetTitle("ROC;Sensitivity;Specificity")
  else:globals()[algoSubList[i]].SetTitle(";Sensitivity;Specificity")
  globals()[algoSubList[i]].SetLineColor(colList[i])
  globals()[algoSubList[i]].SetLineWidth(2)
  for j in range(h.GetNbinsX()):
    globals()[algoSubList[i]].SetPoint(j,0.005+0.01*j,h.GetBinContent(j+1))
  leg.AddEntry(globals()[algoSubList[i]], nameList[i], "l")
  if i==0:globals()[algoSubList[i]].Draw("AL")
  else:globals()[algoSubList[i]].Draw("same")
leg.Draw()
f.Close()

c.SaveAs("ROC.pdf")
