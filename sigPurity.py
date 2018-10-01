from ROOT import *
from array import array

gROOT.SetBatch(True)
c = TCanvas("c","c",900,300)
c.Divide(3)
col = ["#660000","#ffcc00","#cc6600","#ff0000","#ff6565",6]
nbjets2 = [534.1, 1328.8, 542.1, 12809.4, 1809.3, 183.7+511.5+2.6+0.6+0.0+5.0]
nbjetsL3 = [446.9, 887.6, 334.5, 3736.4, 685.1, 84.7+188.8+1.3+0.4+0.0+0]
nbjets3 = [279.7, 425.3, 103.6, 556.3, 114.3, 26.2+42.5+0.0+0.1+0.0+0.0]
entry = [nbjets2,nbjetsL3,nbjets3]
label = ["t#bar{t}b#bar{b}","t#bar{t}bj","t#bar{t}c#bar{c}","t#bar{t}lf","t#bar{t}others", "bkg"]
pie1 = TPie("pie1","nMbjets#geq2",6)
pie2 = TPie("pie2","nMbjets#geq2 & nLbjets#geq3",6)
pie3 = TPie("pie3","nMbjets#geq3",6)

pieList = [pie1,pie2,pie3]
for i in range(3):
  c.cd(i+1) 
  pieList[i].SetRadius(0.35)
  pieList[i].SetLabelFormat("#splitline{%txt}{%perc}")
  for j in range(6):
    pieList[i].SetEntryVal(j,entry[i][j])
    pieList[i].SetEntryFillColor(j,TColor.GetColor(col[j]))
    pieList[i].SetEntryLabel(j,label[j])

  pieList[i].SetEntryRadiusOffset(0,0.05)
  pieList[i].Draw()
  c.Update()
c.SaveAs("sigPurity.png")
