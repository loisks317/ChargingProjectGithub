# MagLatNonEclipse.py
#
# exploring what is happening at high magnetic latitude, positive charging
# times that are NOT during eclipse
#
# LKS, August 2015
#
# imports
import h5py
import numpy as np
from spacepy import pycdf
import os
import glob
import datetime
import matplotlib.dates as dt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
# data
stuff=['density', 'epoch', 'L', 'MLAT', 'MLT', 'potential', 'eclipseFlags']
var=['electronFluxPA', 'electronFlux', 'electronEnergy', 'epoch', 'L', 'MLT', 'eclipseFlags']
hourGap=5
#
dataEFW={}
# open the files
# get the EFW data
os.chdir('Combined_EFW')
for iStuff in stuff:
    dataEFW[iStuff]=h5py.File('combined_EFW_L3_A'+'_'+iStuff+'.h5', 'r')[iStuff]
os.chdir('..')
#
# combined variables
ECL0=np.where(np.array(dataEFW['eclipseFlags'])[:,0]==1)[0]
ECL1=np.where(np.array(dataEFW['eclipseFlags'])[:,1]==1)[0]
ECL2=np.where(np.array(dataEFW['eclipseFlags'])[:,2]==1)[0]
#
#
MLAT=np.array(dataEFW['MLAT'])
potential=-1.0*np.array(dataEFW['potential'])
#
# nan the bad indexes
MLAT[ECL0]=np.nan; MLAT[ECL1]=np.nan; MLAT[ECL2]=np.nan
potential[ECL0]=np.nan; potential[ECL1]=np.nan; potential[ECL2]=np.nan
#
# now looks for where MLAT > 5 and potential > 10
pos=np.where((MLAT > 5) & (potential>10))[0]
#
# get the unique dates from this
dates=dt.num2date(np.array(dataEFW['epoch'])[pos])
listOfDates=[datetime.datetime.strftime(i,'%Y%m%d') for i in dates]
uniqueDates=list(set(listOfDates))
uniqueDates.sort()
