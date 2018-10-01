from ROOT import *

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

#f = TFile("/xrootd/store/user/chanwook/ttbb/final/preapproval/test/SingleMuon_Run2016H_v3.root","r")
f = TFile("TtbarBbbarDiLeptonAnalyzer_SingleElectron_Run2016H_v3.root","r")
t = f.Get("cattree/nom2")


fDi = open("eventDi.txt", "w")
fSi = open("eventSi.txt", "w")
fBoth = open("eventBoth.txt", "w")

fDu = open("eventDu.txt", "w")
fNotTriDu = open("eventNotTriDu.txt", "w")
fDuDi = open("eventDuDi.txt", "w")
fDuSi = open("eventDuSi.txt", "w")
fDuBoth = open("eventDuBoth.txt", "w")

arr = []
f1 = open("event.txt", "r")
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
  if t.tri_single<=0.: continue
  if t.tri_di>0.: continue
  idN = t.event*pow(10,len(str(t.run))) + t.run
  if idN in arr:
    if t.tri_single>0: fDu.write(str(idN)+'\n') 
    else: fNotTriDu.write(str(idN)+'\n')
    if t.tri_singleM>0 and t.tri_singleE<=0: fDuDi.write(str(idN)+'\n')
    elif t.tri_singleM<=0 and t.tri_singleE>0: fDuSi.write(str(idN)+'\n')
    elif t.tri_singleM>0 and t.tri_singleE>0: fDuBoth.write(str(idN)+'\n')
  else:
    if t.tri_singleM>0 and t.tri_singleE<=0: fDi.write(str(idN)+'\n')
    elif t.tri_singleM<=0 and t.tri_singleE>0: fSi.write(str(idN)+'\n')
    elif t.tri_singleM>0 and t.tri_singleE>0: fBoth.write(str(idN)+'\n')
  
fDi.close()
fSi.close()
fBoth.close()
fDu.close()
fNotTriDu.close()
fDuDi.close()
fDuSi.close()
fDuBoth.close()
