from ROOT import *

jesList = ["AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeFSR", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF"]

#for jes in jesList:
for pdf in range(102):
  #print '{"name":"JES_%s_Up", "tree":"JES_%s_up", "var":"("+baseWeight+"*csvweights2[0])"},' % (jes, jes)
  #print '{"name":"JES_%s_Down", "tree":"JES_%s_dw", "var":"("+baseWeight+"*csvweights2[0])"},' % (jes, jes)
  print 'pdf%i shape2N  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1' % pdf
