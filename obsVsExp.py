from ROOT import *
from array import array

#gStyle.SetPadTickY(kFALSE)
line = TLine(0,0,0.5,0.5)

c = TCanvas("c","c",500,500)
ax = [3]
axx = [3]
ay = [4]
aexl = [1]
aexll = [0.5]
aeyl = [4]
gae = TGraphErrors(1, array("d",ax), array("d",ay), array("d",aexl), array("d",aeyl))
gae.SetFillColor(2)
gae.SetFillStyle(3001)
gae.GetYaxis().SetTickSize(0)
gae.GetYaxis().SetLabelSize(0)
gae.Draw("a2")

ga = TGraphErrors(1, array("d",axx), array("d",ay), array("d",aexll), array("d",aeyl))
ga.SetFillColor(3)
ga.SetFillStyle(3001)
ga.Draw("2")

