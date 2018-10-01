import sys, os
from ROOT import *

#analysis = '0000'
analysis = 'TtbarBbbarDiLeptonAnalyzer'
#analysis = 'TtbarBbbarLeptonJetsAnalyzer'
#analysis = 'CSVTemplateGetter'
#analysis = 'TTBBGenAnalyzer'
#rdlist = ['MuonEG_Run2015', 'DoubleEG_Run2015', 'DoubleMuon_Run2015']
for diritem in sorted(os.listdir(".")):
    #if not diritem.startswith(analysis): continue
    #if diritem in skiplist: continue
    if diritem == 'JES': continue
    if '.root' not in diritem: continue
    if 'part1' not in diritem: continue
    name = diritem[0:-11]
    commandhadd = "hadd %s.root " % (name)
    for diritem1 in sorted(os.listdir(".")):
        if '.root' not in diritem1: continue
        if name != diritem1[0:-11]: continue
        commandhadd += diritem1 
        commandhadd += " "
    if not os.path.exists("%s.root" % name): os.system(commandhadd)
    #for rootitem in sorted(os.listdir(diritem)):
    #    if rootitem.endswith(".root"):
    #        commandhadd += " %s/%s" % (diritem, rootitem)
    #os.system("rm -rf %s" % diritem)
    print "@@ %s created" % diritem
