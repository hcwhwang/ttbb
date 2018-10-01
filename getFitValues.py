import os 
from sysWeight_cfi import *
from ROOT import *

for i,d in enumerate(mceventweight):
  os.system("python -i csvFit.py %i" % i)
