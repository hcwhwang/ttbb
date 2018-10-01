import sys, os
import ROOT

step = [1, 2, 3, 4, 5, 6]
var = ['nvertex', 'll_m', 'met', 'njet30', 'nbjetM30', 'nbjetT30']
binning = ['[40,0,40]', '[60,20,320]', '[30,0,150]', '[10,0,10]', '[10,0,10]', '[10,0,10]']
x_name = ['no. vertex', 'll mass (GeV/c^{2})', 'MET (GeV/c)', 'no. jet', 'no. M b-jet', 'no. T b-jet']


weight = "tri*weight*puweight*mueffweight*eleffweight*csvweights2[0]"
#weight1 = "'weight*puweight*mueffweight*eleffweight*tri'"
weight1 = "tri*weight*puweight*mueffweight*eleffweight"
weight2 = "tri*weight*mueffweight*eleffweight"

os.system("python ttLepJetsDraw.py  -s 1  -b [40,0,40] -c 'filtered==1&&tri>0' -p nvertex -x 'no. vertex' -w %s  -d -t 'puweight'" % weight1)
os.system("python ttLepJetsDraw.py  -s 1  -b [40,0,40] -c 'filtered==1&&tri>0' -p nvertex -x 'no. vertex' -w %s -t 'no_puweight' " % weight2)
os.system("python ttLepJetsDraw.py  -s 1 -b [60,0,300] -c 'filtered==1&&tri>0' -p lep_pt -x 'lepton p_{T} (GeV/c)' -w %s -d -t 'lep1_pt'" % weight1)
os.system("python ttLepJetsDraw.py  -s 1 -b [60,0,300] -c 'filtered==1&&tri>0' -p met -x 'Missing E_{T} (GeV/c)' -w %s -d -t 'met'" % weight1)
#os.system("python ttLepJetsDraw.py  -s 1ered==1' -b [40,0,40] -p nvertex -x 'no. vertex' -w 1 " )
os.system("python ttLepJetsDraw.py  -s 2 -b [60,0,300] -c 'filtered==1&&tri>0' -p met -x 'Missing E_{T} (GeV/c)' -w %s -d -t 'met'" % weight1)

os.system("python ttLepJetsDraw.py  -s 1 -b [10,0,10] -c 'filtered==1&&tri>0' -p njet25 -x 'Jet25 Multiplicity' -w %s -d -t 'njet'" % weight1 )
os.system("python ttLepJetsDraw.py  -s 2 -b [10,0,10] -c 'filtered==1&&tri>0' -p njet25 -x 'Jet25 Multiplicity' -w %s -d -t 'njet'" % weight )
os.system("python ttLepJetsDraw.py  -s 3 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetM25 -x 'b Jet25 Midium Multiplicity' -w %s -d -t 'nbjetM'" % weight)
os.system("python ttLepJetsDraw.py  -s 4 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetM25 -x 'b Jet25 Midium Multiplicity' -w %s -d -t 'nbjetM'" % weight)
#os.system("python ttLepJetsDraw.py  -s 3 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x 'CSVv2 3rd' -w %s -d -t 'csv_3rd_weight'" % weight)
#os.system("python ttLepJetsDraw.py  -s 3 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x 'CSVv2 4th' -w %s -d -t 'csv_4th_weight'" % weight)
#os.system("python ttLepJetsDraw.py  -s 3 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x 'CSVv2 3rd' -w %s -d -t 'csv_3rd_noweight'" % weight1)
#os.system("python ttLepJetsDraw.py  -s 3 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x 'CSVv2 4th' -w %s -d -t 'csv_4th_noweight'" % weight1)
#os.system("python ttLepJetsDraw.py  -s 5 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetT30 -x 'b Jet30 Tight Multiplicity' -w weight*puweight -d " )
#os.system("python ttLepJetsDraw.py  -s 6 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetM30 -x 'b Jet30 Midium Multiplicity' -w weight*puweight -d " )
#os.system("python ttLepJetsDraw.py  -s 6 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetT30 -x 'b Jet30 Tight Multiplicity' -w weight*puweight -d " )


