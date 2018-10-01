from ROOT import *
import copy
'''
fJES = TFile.Open("csvNbjets2VisNestPDFIndi1_nuisance.root")
for i in fJES.GetListOfKeys():
  if i.GetClassName()!= 'TH1F': continue
  if 'data' in i.GetTitle() : continue
  #if 'FSR' not in i.GetTitle() and 'ISR' not in i.GetTitle(): continue
  #if 'PS' in i.GetTitle(): continue
  if 'PDF' not in i.GetTitle(): continue
  print i.GetTitle()
  h = fJES.Get(i.GetTitle())
  hh = copy.deepcopy(h)
  f = TFile.Open("csvNbjets2VisNest_nuisance.root","UPDATE")
  h.Write()
  f.Close()
'''
f = TFile.Open("csvNbjets2VisNestFullTot_nuisance.root","UPDATE")
for i in f.GetListOfKeys():
  if i.GetClassName()!= 'TH1F': continue
  #if 'UE' not in i.GetTitle() and 'FSR' not in i.GetTitle() and 'ISR' not in i.GetTitle(): continue
  #if 'JES' not in i.GetTitle(): continue
  if 'JES' in i.GetTitle():gDirectory.Delete("%s;1" % i.GetTitle())
'''
  if 'PDF' in i.GetTitle():
    hTitle = i.GetTitle().split("_")[0] + "_" + i.GetTitle().split("_")[1]
    h = f.Get(hTitle)
    h1 = f.Get(i.GetTitle())
    print i.GetTitle(), h.Integral()-h1.Integral(), h1.Integral()
    if i.GetTitle().endswith("Down"): print " "
'''
f.Close()
#gDirectory.Delete("MuMu_data_obs;2")
