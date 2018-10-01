from ROOT import *

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


f = TFile("TtbarBbbarDiLeptonAnalyzer_SingleMuon_Run2016H_v3.root","r")
t = f.Get("cattree/nom2")
fDi = open("eventDoubleDi.txt", "w")
fSi = open("eventDoubleSi.txt", "w")
fBoth = open("eventDoubleBoth.txt", "w")

arr = []
f1 = open("eventDu.txt", "r")
while True:
  value = f1.readline()
  if isfloat(value):
    line = long(float(value))
    print line
  else:
    print value
    break
  arr.append(line)
  if not value: break
f1.close()
for i in range(t.GetEntries()):
  t.GetEntry(i)
  if i % 10000 == 0: print i
  if t.channel!=1: continue
  if t.tri_di>0: continue
  idN = t.event*pow(10,len(str(t.run))) + t.run
  if idN in arr:continue
  if t.tri_singleM>0 and t.tri_singleE<=0: fDi.write(str(idN)+'\n')
  elif t.tri_singleM<=0 and t.tri_singleE>0: fSi.write(str(idN)+'\n')
  elif t.tri_singleM>0 and t.tri_singleE>0: fBoth.write(str(idN)+'\n')
fDi.close()
fSi.close()
fBoth.close()

