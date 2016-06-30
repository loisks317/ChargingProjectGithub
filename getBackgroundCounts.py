# getBackgroundCounts.py
#
# average the counts in each energy channel across 6 months to remove the background
#
# September 2015 LKS
#
#

import os
import glob 
import numpy as np
from matplotlib import pyplot as plt
import datetime
import matplotlib.dates as dates
import spacepy.pybats.kyoto as spk
from spacepy import pycdf
import itertools as itert
import math
from numpy import ma
import pandas as pd

#
# Steps
# 1. read in the file
# 2. sort by L shell (don't worry about MLT), maybe .1 resolution
# 3. take median of each bin in each energy channel
# 4. save that as the 'background' value
#
dateStart='20130201'
dateEnd='20150401'
#
#
EnBins=72
Lbins=56
LbinsArr=np.linspace(1, 6.5,Lbins)
SortedLists=[[ [] for i in range(Lbins)] for j in range(EnBins)]
#
#
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
while DT != endDt:
     date=datetime.datetime.strftime(DT, '%Y%m%d')
     os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_A')
     f=glob.glob('*'+date+'*')
     try:
         pyf=pycdf.CDF(f[0])
         #
         # put into pandas arrays
         energy=np.swapaxes(pyf['HOPE_ENERGY_Ion'][...],1,0)
         L=pyf['L_Ion'][...]
         counts=pyf['Counts_P_Omni'][...]

         for iL in range(Lbins-1):
             for iEn in range(EnBins):
                 temp=np.where((np.array(L) >= LbinsArr[iL]) & (np.array(L) < LbinsArr[iL+1]))[0]
                 try:
                     # make sure this is correct and that temp and iEn are in
                     # the correct order
                     SortedLists[iEn][iL]+=list(counts[temp][iEn])
                 except(IndexError):
                     continue
     except:
         print 'Bad date: ' + date
     DT=DT+datetime.timedelta(days=1)
#
# get the median
print "getting median"
medData=[[ [] for i in range(Lbins)] for j in range(EnBins)]
for iL in range(Lbins-1):
    for iEn in range(EnBins):
        medData[iEn][iL]=np.nanmedian(np.array(SortedLists[iEn][iL]))
import pickle
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
dir_name='BackgroundCounts'
if not os.path.exists(dir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(dir_name, 0777) 
os.chdir(dir_name)#
with open("BgCounts.p", "wb") as f:
      pickle.dump(medData, f)
