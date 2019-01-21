from ROOT import *
import CMS_lumi, tdrstyle
#tdrstyle.setTDRStyle()

#sysList = ["FSR", "ISR", "UE", "qcderdon", "gluonmove", "hdamp", "topPt", "MuR", "PDFAlphaS", "BtagLF", "BtagHF", "BtagHFStats1", "BtagHFStats2", "BtagLFStats2", "BtagCQErr1", "BtagCQErr2", "JER", "PU", "ElScale", "ElEff", "dy", "ttv", "Lumi", "JESFlavorQCD", "JESRelativePtEC1", "JESRelativeBal", "JESPileUpDataMC", "JESPileUpPtRef", "JESPileUpPtEC1","r","k"]
sysList = [u'BtagLF', u'BtagHF', u'gluonmove', u'FSR', u'BtagCQErr1', u'BtagHFStats2', u'BtagLFStats2', u'JESFlavorQCD', u'hdamp', u'erdon', u'ISR', u'Lumi', u'PDFAlphaS', u'Trigger', u'BtagCQErr2', u'JESRelativeBal', u'prop_binMuEl_bin99', u'st', u'topPt', u'MuR', u'BtagHFStats1', u'JESRelativeFSR', u'prop_binMuEl_bin19', u'JESPileUpPtRef', u'UE', u'JESPileUpDataMC', u'prop_binMuEl_bin29', u'prop_binMuEl_bin39', u'prop_binElEl_bin79', u'prop_binMuEl_bin79', u'prop_binMuEl_bin49']
sysTotList = []
gStyle.SetOptStat(0)
#phase='full'
phase='vis'
#f = TFile.Open("fitDiagnostics%s.root"%phase)
f = TFile.Open("fitDiagnostics.root")
h = f.Get("covariance_fit_s")
h.GetXaxis().SetLabelSize(0.03)
h.GetYaxis().SetLabelSize(0.03)
h.GetZaxis().SetLabelSize(0.03)
p = h.GetListOfFunctions().FindObject("palette")
p.SetX1NDC(0.92)
p.SetX2NDC(0.95)
p.SetY1NDC(0.1)
p.SetY2NDC(0.9)
h1 = TH1D("h1","nuisance-r;Correlation Factor;Arbitrary Unit",50,-1,1)
h1.SetLineColor(kRed)
h1.SetLineWidth(2)

h2 = TH1D("h2","nuisance-k;Correlation Factor;Arbitrary Unit",50,-1,1)
h2.SetLineColor(kBlue)
h2.SetLineWidth(2)

h3 = TH1D("h3","nuisance-nuisance;Correlation Factor;Arbitrary Unit",50,-1,1)
h3.SetLineColor(kBlack)
h3.SetLineWidth(2)

hCor = TH2D("hCor", "Correlations between NPs and POIs", len(sysList), 0, len(sysList), len(sysList), 0, len(sysList))
leg = TLegend(0.65, 0.7, 0.92, 0.92)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(h1,h1.GetTitle(),"l")
leg.AddEntry(h2,h2.GetTitle(),"l")
leg.AddEntry(h3,h3.GetTitle(),"l")

for i in range(h.GetXaxis().GetNbins()):
  sysTotList.append(h.GetXaxis().GetBinLabel(i+1))
for i in range(h.GetXaxis().GetNbins()):
  for j in range(h.GetYaxis().GetNbins()):
    xTitle = h.GetXaxis().GetBinLabel(i+1)
    yTitle = h.GetXaxis().GetBinLabel(j+1)
    if h.GetXaxis().GetBinLabel(i+1)=='r' and h.GetYaxis().GetBinLabel(j+1)=='k':
      print "correlation r and k", h.GetBinContent(i+1,j+1)
    if xTitle in sysList and yTitle in sysList:
      hCor.GetXaxis().SetBinLabel(sysList.index(xTitle)+1, xTitle)
      hCor.GetYaxis().SetBinLabel(sysList.index(yTitle)+1, yTitle)
      hCor.SetBinContent(sysList.index(xTitle)+1, sysList.index(yTitle)+1, h.GetBinContent(i+1, h.GetYaxis().GetNbins()-j))
    if j>=-i+h.GetYaxis().GetNbins()-1: continue
    if h.GetXaxis().GetBinLabel(i+1)=='r' and h.GetYaxis().GetBinLabel(j+1)=='k': continue
    if h.GetXaxis().GetBinLabel(i+1)=='k' and h.GetYaxis().GetBinLabel(j+1)=='r': continue
    if h.GetYaxis().GetBinLabel(j+1)=='r': h1.Fill(h.GetBinContent(i+1,j+1))
    elif h.GetYaxis().GetBinLabel(j+1)=='k': h2.Fill(h.GetBinContent(i+1,j+1))
    else: h3.Fill(h.GetBinContent(i+1,j+1))
    #print h.GetXaxis().GetBinLabel(i+1), h.GetYaxis().GetBinLabel(j+1), h.GetBinContent(i+1,j+1)
    #if h.GetBinContent(i+1,j+1)<0.05: h.SetBinContent(i+1,j+1,0)
'''
c = TCanvas("c","c",1000,700)
h.Draw("colz")
'''
for h in [h1,h2,h3]:
  h.Scale(1./h.Integral())

hFrame = TH1D("hFrame",";Correlation Factor;Arbitrary Unit",50,-1,1)
hFrame.SetMaximum(max(h.GetMaximum() for h in [h1,h2,h3])*1.2) 
c1 = TCanvas("c1","c1",700,700)
hFrame.Draw()
h1.Draw("same")
h2.Draw("same")
h3.Draw("same")
leg.Draw()
c1.SaveAs("correlation1D%s.png"%phase)
c1.SaveAs("correlation1D%s.pdf"%phase)
c2 = TCanvas("c2","c2",1000,700)
hCor.GetXaxis().SetLabelSize(0.025)
hCor.GetXaxis().SetLabelOffset(0.001)
hCor.GetYaxis().SetLabelSize(0.025)
hCor.GetYaxis().SetLabelOffset(0.001)
hCor.GetZaxis().SetLabelSize(0.025)
hCor.GetZaxis().SetRangeUser(-1, 1)
hCor.GetXaxis().LabelsOption("v")
hCor.Draw("colz")
c2.SaveAs("correlation%s.png"%phase)
c2.SaveAs("correlation%s.pdf"%phase)

