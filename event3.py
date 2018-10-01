from ROOT import *

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
arr = []
f1 = open("eventSiM.txt", "r")
while True:
  value = f1.readline()
  if isfloat(value):
    line = long(float(value))
    print line
  else:
    print value
    break
  arr.append(line)
  if not value: break
f1.close()

f2 = open("eventSiE.txt", "r")
nDu = 0
while True:
  value = f2.readline()
  if isfloat(value):
    line = long(float(value))
    print line
  else:
    print value
    break
  if line in arr:
    nDu += 1
    print line
  if not value: break
print nDu
f2.close()



