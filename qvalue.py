from ROOT import *
from math import *

s = 10.
b = 10.
tau = 1.
mu = 0.
muTest = 0.
def muHat(n,m,tau,s):
  return (n-m/tau)/s

def bHat(m,tau):
  return m/tau

def bHatHat(n,m,tau,s,mu):
  a = (n+m-(1+tau)*mu*s)/(2*(1+tau))
  b = ((n+m-(1+tau)*mu*s)**2+4*(1+tau)*m*mu*s)/(4*(1+tau)**2)
  return a+pow(b,0.5)
def L(a,b):
  #print a**b, exp(-a), factorial(b), (a**b)*exp(-a)/factorial(b)
  return (a**b)*exp(-a)/factorial(b) 

h = TH1F("h","h",25,0,40)
h.SetMaximum(10)
h.SetMinimum(pow(10,-8))
for i in range(10000000):
  p1 = TRandom(i)
  p2 = TRandom(i+1)
  pN = p1.Poisson(mu*s+b)
  pM = p2.Poisson(tau*b)
  muH = muHat(pN,pM,tau,s)
  bH = bHat(pM,tau)
  bHH = bHatHat(pN,pM,tau,s,muTest)
  L1 = L(muTest*s+bHH,pN)*L(tau*bHH,pM)
  #print L(muTest*s+bHH,pN), L(tau*bHH,pM)
  L2 = L(muH*s+bH,pN)*L(tau*bH,pM)
  if muH >= 0: value = -2*log(L1/L2)
  else: value = 0
  h.Fill(value)
c = TCanvas("c","c",500,500)
c.SetLogy()
h.Scale(1./h.Integral())
h.Draw()

