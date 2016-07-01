# sensitivityRead.py
#
# Parse sensitivtyTestOutput.py files
# find highest success rate vs lowest number of instances
#
# LKS, September 2015
#
import numpy as np
fname='sensitivityTestOutput.txt'

with open(fname) as f:
    content = f.readlines()

successRate=[]
instances=[]
TperpC=[]
TparC=[]
PressureC=[]
DensityC=[]

for i in range(len(content)):
    if content[i][0:7]=='Success':
        successRate.append(float(content[i][38:-2]))
    elif content[i][0:5]=='total':
        instances.append(float(content[i][20:-1]))
    elif content[i][0:3]=='for':
        temp=content[i].split()
        TperpC.append(float(temp[1][6:-1]))
        TparC.append(float(temp[2][5:-1]))
        PressureC.append(float(temp[3][9:-1]))
        DensityC.append(float(temp[4][8:]))
successRate=np.array(successRate)
high=np.where(successRate > .8)[0]
instances=np.array(instances)
highIn=np.min(instances[high])
rightIndex=np.where((instances==highIn)&(successRate > 0.8))[0][0]

print "Tperp = "+ str(TperpC[rightIndex])
print "Tpar = " + str(TparC[rightIndex])
print "Pressure = " + str(PressureC[rightIndex])
print "DensityC = " + str(DensityC[rightIndex])


        
    
