import sys, os
import ROOT

weight = "tri*weight*puweight*mueffweight*eleffweight*csvweights2[0]"
#weight1 = "'weight*puweight*mueffweight*eleffweight*tri'"
weight1 = "tri*weight*puweight*mueffweight*eleffweight"
weight2 = "tri*weight*mueffweight*eleffweight"

#os.system("python histSaver.py  -s 5 -b '[10,0,1.]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x '3rd b-discriminator' -w %s -d -t 'csv_3rd_afterfit'" % weight)
#os.system("python histSaver.py  -s 5 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x '4rd b-discriminator' -w %s -d -t 'csv_4th_afterfit'" % weight)
os.system("python histSaver.py  -s 1  -b '[40,0,40]' -c 'filtered==1&&tri>0 ' -p nvertex -x 'primary vertex' -w %s  -d -t 'puweight'" % weight1)
os.system("python histSaver.py  -s 1  -b '[40,0,40]' -c 'filtered==1&&tri>0' -p nvertex -x 'primary vertex' -w %s -d -t 'no_puweight' " % weight2)
os.system("python histSaver.py  -s 1 -b '[30,0,300]' -c 'filtered==1&&tri>0 ' -p ll_m -x 'dilepton invariant mass (GeV)' -w %s -d  -t 'll_m'" % weight1 )
os.system("python histSaver.py  -s 2 -b '[30,0,300]' -c 'filtered==1&&tri>0 ' -p ll_m -x 'dilepton invariant mass (GeV)' -w %s -d -t 'll_m'" % weight1 )
os.system("python histSaver.py  -s 3 -b '[10,0,10]' -c 'filtered==1&&tri>0' -p njet30 -x 'jet multiplicity' -w %s -d -t 'njet' " % weight1 )
os.system("python histSaver.py  -s 2 -b '[30,0,300]' -c 'filtered==1&&tri>0' -p met -x 'missing E_{T} (GeV)' -w %s -d -t 'met'" % weight1)
#os.system("python histSaver.py  -s 1ered==1' -b [40,0,40] -p nvertex -x 'no. vertex' -w 1 " )
os.system("python histSaver.py  -s 3 -b '[30,0,300]' -c 'filtered==1&&tri>0' -p met -x 'missing E_{T} (GeV)' -w %s -d -t 'met'" % weight1)

os.system("python histSaver.py  -s 3 -b '[10,0,10]' -c 'filtered==1&&tri>0' -p njet30 -x 'jet multiplicity' -w %s -d -t 'njet'" % weight1 )
os.system("python histSaver.py  -s 4 -b '[10,0,10]' -c 'filtered==1&&tri>0' -p njet30 -x 'jet multiplicity' -w %s -d -t 'njet'" % weight )
os.system("python histSaver.py  -s 4 -b '[6,0,6]' -c 'filtered==1&&tri>0' -p nbjetM30 -x 'b-tagged jet multiplicity' -w %s -d -t 'nbjet'" % weight)
os.system("python histSaver.py  -s 5 -b '[6,0,6]' -c 'filtered==1&&tri>0' -p nbjetM30 -x 'b-tagged jet multiplicity' -w %s -d -t 'nbjet'" % weight)
os.system("python histSaver.py  -s 4 -b '[10,0,1.]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x '3rd b-discriminator' -w %s -d -t 'csv_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x '4rd b-discriminator' -w %s -d -t 'csv_4th_noweight'" % weight1)
#os.system("python histSaver.py  -s 5 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetT30 -x 'b Jet30 Tight Multiplicity' -w weight*puweight -d " )
#os.system("python histSaver.py  -s 6 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetM30 -x 'b Jet30 Midium Multiplicity' -w weight*puweight -d " )
#os.system("python histSaver.py  -s 6 -b [10,0,10] -c 'filtered==1&&tri>0' -p nbjetT30 -x 'b Jet30 Tight Multiplicity' -w weight*puweight -d " )


os.system("python histSaver.py  -s 3 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[0]] -x 'CSVv2 1st' -w %s -d -t 'csv_1st_noweight'" % weight1)
os.system("python histSaver.py  -s 3 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[1]] -x 'CSVv2 2nd' -w %s -d -t 'csv_2nd_noweight'" % weight1)
os.system("python histSaver.py  -s 3 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[0]] -x 'CSVv2 1st' -w %s -d -t 'csv_1st_weight'" % weight)
os.system("python histSaver.py  -s 3 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[1]] -x 'CSVv2 2nd' -w %s -d -t 'csv_2nd_weight'" % weight)
os.system("python histSaver.py  -s 4 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[0]] -x 'CSVv2 1st' -w %s -d -t 'csv_1st_weight'" % weight)
os.system("python histSaver.py  -s 4 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[1]] -x 'CSVv2 2nd' -w %s -d -t 'csv_2nd_weight'" % weight)
os.system("python histSaver.py  -s 4 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x 'CSVv2 3rd' -w %s -d -t 'csv_3rd_weight'" % weight)
os.system("python histSaver.py  -s 4 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x 'CSVv2 4th' -w %s -d -t 'csv_4th_weight'" % weight)
os.system("python histSaver.py  -s 5 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x 'CSVv2 3rd' -w %s -d -t 'csv_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x 'CSVv2 4th' -w %s -d -t 'csv_4th_noweight'" % weight1)

os.system("python histSaver.py  -s 5 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[2]] -x 'CSVv2 3rd' -w %s -d -t 'csv_3rd_weight'" % weight)
os.system("python histSaver.py  -s 5 -b '[10,0,1]' -c 'filtered==1&&tri>0' -p jets_bDiscriminatorCSV[csvd_jetid[3]] -x 'CSVv2 4th' -w %s -d -t 'csv_4th_weight'" % weight)
os.system("python histSaver.py  -s 1 -b '[15,0,150]' -c 'filtered==1&&tri>0' -p lep1_pt -x '1st leading lepton p_{T} (GeV)' -w %s -d -t 'lep1_pt'" % weight1)
os.system("python histSaver.py  -s 1 -b '[15,0,150]' -c 'filtered==1&&tri>0' -p lep2_pt -x '2nd leading lepton p_{T} (GeV)' -w %s -d -t 'lep2_pt'" % weight1)

os.system("python histSaver.py  -s 1 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p lep1_eta -x '1st leading lepton #eta' -w %s -d -t 'lep1_eta'" % weight1)
os.system("python histSaver.py  -s 1 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p lep2_eta -x '2nd leading lepton #eta' -w %s -d -t 'lep2_eta'" % weight1)

os.system("python histSaver.py  -s 1 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p lep1_phi -x '1st leading lepton #phi' -w %s -d -t 'lep1_phi'" % weight1)
os.system("python histSaver.py  -s 1 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p lep2_phi -x '2nd leading lepton #phi' -w %s -d -t 'lep2_phi'" % weight1)

os.system("python histSaver.py  -s 4 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[0]] -x '1st CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_1st_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[1]] -x '2nd CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_2nd_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[2]] -x '3rd CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[3]] -x '4th CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_4th_noweight'" % weight1)

os.system("python histSaver.py  -s 4 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[0]] -x '1st CSV leading jet #eta' -w %s -d -t 'jet_eta_1st_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[1]] -x '2nd CSV leading jet #eta' -w %s -d -t 'jet_eta_2nd_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[2]] -x '3rd CSV leading jet #eta' -w %s -d -t 'jet_eta_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[3]] -x '4th CSV leading jet #eta' -w %s -d -t 'jet_eta_4th_noweight'" % weight1)

os.system("python histSaver.py  -s 4 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[0]] -x '1st CSV leading jet #phi' -w %s -d -t 'jet_phi_1st_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[1]] -x '2nd CSV leading jet #phi' -w %s -d -t 'jet_phi_2nd_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[2]] -x '3rd CSV leading jet #phi' -w %s -d -t 'jet_phi_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 4 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[3]] -x '4th CSV leading jet #phi' -w %s -d -t 'jet_phi_4th_noweight'" % weight1)

os.system("python histSaver.py  -s 5 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[0]] -x '1st CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_1st_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[1]] -x '2nd CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_2nd_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[2]] -x '3rd CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[20,0,200]' -c 'filtered==1&&tri>0' -p jets_pt[csvd_jetid[3]] -x '4th CSV leading jet p_{T} (GeV)' -w %s -d -t 'jet_pt_4th_noweight'" % weight1)

os.system("python histSaver.py  -s 5 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[0]] -x '1st CSV leading jet #eta' -w %s -d -t 'jet_eta_1st_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[1]] -x '2nd CSV leading jet #eta' -w %s -d -t 'jet_eta_2nd_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[2]] -x '3rd CSV leading jet #eta' -w %s -d -t 'jet_eta_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[10,-3,3]' -c 'filtered==1&&tri>0' -p jets_eta[csvd_jetid[3]] -x '4th CSV leading jet #eta' -w %s -d -t 'jet_eta_4th_noweight'" % weight1)

os.system("python histSaver.py  -s 5 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[0]] -x '1st CSV leading jet #phi' -w %s -d -t 'jet_phi_1st_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[1]] -x '2nd CSV leading jet #phi' -w %s -d -t 'jet_phi_2nd_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[2]] -x '3rd CSV leading jet #phi' -w %s -d -t 'jet_phi_3rd_noweight'" % weight1)
os.system("python histSaver.py  -s 5 -b '[10,0,3.14]' -c 'filtered==1&&tri>0' -p jets_phi[csvd_jetid[3]] -x '4th CSV leading jet #phi' -w %s -d -t 'jet_phi_4th_noweight'" % weight1)
