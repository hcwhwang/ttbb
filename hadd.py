import sys, os
from ROOT import *

#analysis = '0000'
analysis = ['TtbarBbbarDiLeptonDoubleRDAnalyzer', 'TtbarBbbarDiLeptonSingleMuonAnalyzer', 'TtbarBbbarDiLeptonSingleElectronAnalyzer', 'TtbarBbbarDiLeptonAnalyzer', 'TtbarBbbarDiLeptonJESAnalyzer']
#analysis = 'TtbarBbbarLeptonJetsAnalyzer'
#analysis = 'CSVTemplateGetter'
#analysis = 'TTBBGenAnalyzer'
#rdlist = ['MuonEG_Run2015', 'DoubleEG_Run2015', 'DoubleMuon_Run2015']
rdlist = ['DoubleEG_Run2016', 'DoubleMuon_Run2016', 'MuonEG_Run2016', 'SingleElectron_Run2016', 'SingleMuon_Run2016']
skiplist = ['TT_powheg_UEup','TT_powheg_UEdown','TT_powheg_FSRup','TT_powheg_FSRdown','TT_powheg_ISRup','TT_powheg_ISRdown']
partlist = [['DYJets.root','DYJets_10to50_part1.root','DYJets_10to50_part2.root','DYJets_10to50_part3.root'],['TTW_part1.root','TTW_part2.root'],['TTZ_part1.root','TTZ_part2.root','TTZ_part3.root'],['WJets_part1.root','WJets_part2.root'],['WW_part1.root','WW_part2.root'],['WZ_part1.root','WZ_part2.root'],['ZZ_part1.root','ZZ_part2.root']]
for diritem in sorted(os.listdir(".")):
    #if analysis not in diritem: continue
    if not any(ana in diritem for ana in analysis): continue
    #if not diritem.startswith(analysis): continue
    #if diritem in skiplist: continue
    newname = diritem[len(diritem.split("_")[0])+1:]
    if diritem == 'JES': continue
    if not os.path.isdir(diritem): continue
    if diritem.endswith(".py") or diritem.endswith(".root"): continue
    if os.path.exists("%s.root" % newname): continue
    os.system("cd %s" % diritem)
    #dataset = diritem[len(analysis)+1:]
    commandhadd = "hadd %s/%s.root" % (diritem, newname) 
    for rootitem in sorted(os.listdir(diritem)):
        if rootitem.endswith(".root"):
            commandhadd += " %s/%s" % (diritem, rootitem)
    os.system(commandhadd)
    #os.system("rm -f %s/%s.root" % (diritem, newname))
    os.system("mv -f %s/%s.root ." % (diritem, newname))
    #os.system("rm -rf %s" % diritem)
    print "@@ %s.root created" % diritem
'''
for i,rd in enumerate(rdlist):
    #rd = analysis + "_" + rd
    tot = rd + '.root'
    b = rd + 'B.root'
    c = rd + 'C.root'
    d = rd + 'D.root'
    e = rd + 'E.root'
    f = rd + 'F.root'
    g = rd + 'G.root'
    h2 = rd + 'H_v2.root'
    h3 = rd + 'H_v3.root'
    if not os.path.exists(tot): os.system("hadd %s %s %s %s %s %s %s %s %s" % (tot,b,c,d,e,f,g,h2,h3))
    #os.system("rm -f %s %s " % (c,d))
'''
print "@@ completed!!"
