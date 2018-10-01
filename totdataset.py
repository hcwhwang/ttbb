from ROOT import *
import json, os, sys

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ["CMSSW_BASE"]))
output = []
for i,name in enumerate(datasets):
  print name["name"]
  if 'part1' in name["name"]:
    suboutput = {"colour":name["colour"], "name":name["name"][0:-6], "title":name["title"], "xsec":name["xsec"]}
  elif 'part' not in name["name"]:
    suboutput = {"colour":name["colour"], "name":name["name"], "title":name["title"], "xsec":name["xsec"]}
  else: continue
  output.append(suboutput)
file = open("%s/src/CATTools/CatAnalyzer/data/dataset/totdataset.json" % os.environ["CMSSW_BASE"], "w")
file.write(json.dumps(output, indent=4))
file.close()

