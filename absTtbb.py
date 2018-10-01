from ROOT import *

def aTimesb(a, aErr, b, bErr, corAB):
  err = b*b*aErr*aErr + a*a*bErr*bErr + 2*a*b*corAB*aErr*bErr
  return pow(err,0.5)/(a*b)*100, pow(err,0.5)

rVis = 0.018
kVis = 2.176
rVisStat = rVis*4.92/100
rVisSys = rVis *7.58/100
kVisStat = kVis*0.87/100
kVisSys = kVis*6.68/100
corVisRK = 0.325

rFull = 0.020
kFull = 145.36
rFullStat = rFull*5.01/100
rFullSys = rFull *7.34/100
kFullStat = kFull*0.87/100
kFullSys = kFull*7.34/100
corFullRK = 0.45

print "vis ttbb = ", rVis*kVis
print "vis ttbb stat = " , aTimesb(rVis, rVisStat, kVis, kVisStat, corVisRK)  
print "vis ttbb sys = " , aTimesb(rVis, rVisSys, kVis, kVisSys, corVisRK) 
print " " 
print "full ttbb = ", rFull*kFull
print "full ttbb stat = " , aTimesb(rFull, rFullStat, kFull, kFullStat, corFullRK) 
print "full ttbb sys = " , aTimesb(rFull, rFullSys, kFull, kFullSys, corFullRK) 
