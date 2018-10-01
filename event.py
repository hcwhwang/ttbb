from ROOT import *

f = TFile("TtbarBbbarDiLeptonAnalyzer_SingleMuon_Run2016H_v3.root","r")
t = f.Get("cattree/nom2")

f1 = open("event.txt","w")
arr = []
for i in range(t.GetEntries()):
#for i in range(1000000):
  t.GetEntry(i)
  if i % 10000 == 0: print i
  if t.channel!=1: continue
  #if t.tri_single<=0.: continue
  #if t.tri_di<=0.: continue
  if t.tri_di>0: continue
  if t.tri_single<=0: continue
  idN = t.event*pow(10,len(str(t.run))) + t.run
  #print idN
  if idN in arr: print "duplicated ", idN
  else: 
    arr.append(idN)
    f1.write(str(idN)+'\n')
    #print str(idN)
f1.close()
