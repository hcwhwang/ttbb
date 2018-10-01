from ROOT import *

r = TRandom()
r1 = TRandom()

h = TH1F("h","h", 100, -50, 50)
h1 = TH1F("h1","h1", 100, -50, 50)
hNew = TH1F()
for i in range(1000):
  rN = r.Gaus(10,10)
  h.SetLineColor(kRed)
   
  h.Fill(rN)

for i in range(10000):
  rN1 = r1.Uniform(-50,50)
  h1.SetLineColor(kBlue)
   
  h1.Fill(rN1)

hNew.Add(h)
hNew.Add(h1)
hNew.Draw()
