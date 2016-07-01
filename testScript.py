# testScript.py
#
# figure out dates of anamolous dayside negative charging event
#
# LKS, Halloween 2015! Boo!
#
#
import numpy as np
from spacepy import pycdf
import os
import glob
import datetime
import matplotlib.dates as dt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from spacepy import datamodel as dm
from numpy import ma

#
#
os.chdir('EFW_L3_A')
files=glob.glob('*.cdf')
for i in range(len(files)):
    try:
        pyf=pycdf.CDF(files[i])
        potential=-1*np.array(pyf['Vavg'][...])
        epoch=np.array(pyf['epoch'][...])
        MLT=np.array(np.swapaxes(pyf['mlt_lshell_mlat'],1 ,0)[0])
        lowEvents=np.where(potential<-10)[0]
        lEpoch=epoch[lowEvents]
        lMLT=MLT[lowEvents]
        MLTday=np.where((lMLT > 8) & (lMLT < 16))[0]
        print lEpoch[MLTday]
    except(OSError):
        print "Bad data"

    
