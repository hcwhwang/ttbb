from ROOT import *
gStyle.SetOptStat(0)
c = TCanvas("c","c",500,500)


f = TFile.Open("DoubleEG_Run2016B.root")
t = f.Get("cattree/nom")
h = TH1F("h","DoubleEG_Run2016B;nvertex;events",40,0,40)
h.SetLineColor(kRed)
t.Project("h","nvertex","nvertex<40")
h.Draw()

ff = TFile.Open("../ICHEP/DoubleEG_Run2016B.root")
#f = TFile.Open("/cms/home/chanwook/work/test/CMSSW_8_0_12/src/CATTools/CatAnalyzer/test/cattree.root")
tt = ff.Get("cattree/nom")
h1 = TH1F("h1","h1",50,0,50)
h1.SetLineColor(kBlue)
tt.Project("h1","nvertex","nvertex<40")
h1.Draw("same")

l = TLegend(0.6,0.6,0.9,0.9)
l.SetFillStyle(0)
l.SetBorderSize(0)
l.AddEntry(h1,"ReReco","l")
l.AddEntry(h,"PromptReco","l")
l.Draw()



#nbin = h.GetXaxis().GetNbins()
#array = ""
#for i in range(nbin):
#  print format(h.GetBinContent(i+1),'10e'),',',
#h.Draw()
