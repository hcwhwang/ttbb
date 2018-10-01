import sys, os
import ROOT

weight = "tri*weight*puweight*mueffweight*eleffweight*csvweights2[0]"
#weight1 = "'weight*puweight*mueffweight*eleffweight*tri'"
weight1 = "tri*weight*puweight*mueffweight*eleffweight"
weight2 = "tri*weight*mueffweight*eleffweight"

#os.system("python csv1DUnroll.py  -s 3 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[0]] -x 'CSVv2 1st' -w %s -d -t 'csv_1st_weight'" % weight)
#os.system("python csv1DUnroll.py  -s 3 -b '[7,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[0]] -x 'CSVv2 1st' -w %s -d -t 'csv_1st_weight_bin7'" % weight)
#os.system("python csv1DUnroll.py  -s 1 -b '[50,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV -x 'CSVv2' -w %s  -d -t 'csv_all_noweight_bin50'" % weight1)
#os.system("python csv1DUnroll.py -s 5 -b '[10,0,1,10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]]:jets_bDiscriminatorCSV[csvd_jetid[3]] -x 'CSVv2 3rd & 4th' -w %s  -d -t 'csv_3rdVs4th_afterfit'" % weight)
#os.system("python csv1DUnroll.py -s 5 -b '[10,0,1,10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]]:jets_bDiscriminatorCSV[csvd_jetid[3]] -x '3rd & 4th b-discriminator' -w %s  -d -t 'csv_3rdVs4th'" % weight)
#os.system("python csv1DUnroll.py  -s 5 -b [20,0,200] -c 'filtered==1&&tri>0' -p jets_pt -x 'Jet p_{T} (GeV)' -w %s -d -t 'jet_pt'" % weight)
#os.system("python csv1DUnroll.py  -s 5 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x '3rd b-discriminator' -w %s -d -t 'csv_3rd_weight'" % weight)
#os.system("python csv1DUnroll.py  -s 5 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x '4th b-discriminator' -w %s -d -t 'csv_4th_weight'" % weight)
#os.system("python csv1DUnroll.py  -s 5 -b '[6,0,6]' -c 'filtered==1&&tri>0' -p nbjetM30 -x 'b-tagged Jet Multiplicity' -w %s -d -t 'nbjetM'" % weight)
os.system("python csv1DUnroll.py  -s 4 -b '[10,0,10]' -c 'filtered==1&&tri>0' -p njet30 -x 'Jet Multiplicity' -w %s -d -t 'njet_notrigger'" % weight)
#os.system("python csv1DUnroll.py  -s 5 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x '3rd b-discriminator' -w %s -d -t 'csv_3rd_afterfit'" % weight)
#os.system("python csv1DUnroll.py  -s 5 -b [10,0,1] -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x '4th b-discriminator' -w %s -d -t 'csv_4th_afterfit'" % weight)
