import json, os, sys

mceventweight = []
baseWeight = "tri*weight*puweight*mueffweight*eleffweight*csvweights2[0]"
for i in range(101):
  out = {"name":"pdf%i"%i, "tree":"nom2", "var":"("+baseWeight+"*pdfWeights[%i])"%i}
  mceventweight.append(out)
file = open("sysWeightPDF_cfi.py", "w")
jsonoutput = json.dumps(mceventweight)
file.write(jsonoutput)
file.close()
