#compareChargingTimes.py
#
# directly compare times of charging in HOPE vs EFW
# put error bars on the HOPE times for width of the energy channel
#
# LKS, October 2015
#
# imports 
import numpy as np
from spacepy import pycdf
import glob
import os
import datetime
import pickle 
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#
# load in the times
EFWSC=pickle.load(open("EFWSC.p", "rb"))
EFWTimes=pickle.load(open("EFWTimes.p","rb"))
IonSC=pickle.load(open("IonAlgoLine.p","rb"))
IonTimes=pickle.load(open("IonAlgoTime.p", "rb"))
#
# get energy channels energy
# then get energy channel widths
os.chdir('HOPE_L3_A')
f=glob.glob('rbspa_rel02_ect-hope-PA-L3_20130329_v5.0.1.cdf')
pyf=pycdf.CDF(f[0])
energyChannels=pyf['HOPE_ENERGY_Ion'][...][0] # typical energy channel
delta=pyf['ENERGY_Ion_DELTA'][...][0]
os.chdir('..')
#
# make sure they are arrays
IonSC=np.array(IonSC)
IonTimes=np.array(IonTimes)
EFWSC=np.array(EFWSC)
EFWTimes=np.array(EFWTimes)
#
# create error bar array
ErrorBars=np.zeros(len(IonSC))
for iEn in range(1,len(energyChannels)-1):
    indices=np.where((IonSC > energyChannels[iEn-1]) & (IonSC <= energyChannels[iEn+1]))[0]
    ErrorBars[indices]=delta[iEn]
#
#


fig=plt.figure()
ax1=fig.add_subplot(111)
plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
ax1.set_ylabel('| Potential [V] |', fontsize=25, fontweight='bold')
ax1.set_xlabel('Date', fontsize=25, fontweight='bold')
ax1.set_yscale('log')
ax1.set_ylim(1e0,1e3)
# 10260 is the last charging event of 2013 before 2014
ax1.set_xlim(EFWTimes[0], EFWTimes[10260])
font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 22}
plt.rc('font', **font)
plt.scatter(EFWTimes, EFWSC, s=30, c='green', edgecolors='none')
plt.scatter(IonTimes[IonSC>10], IonSC[IonSC>10], s=30, c='gold', edgecolors='none')
plt.errorbar(IonTimes[IonSC>10],IonSC[IonSC>10], yerr=ErrorBars[IonSC>10], linestyle='none')
hfmt = dt.DateFormatter('%Y \n -%m')
ax1.set_xlim(EFWTimes[0], EFWTimes[10260])
ax1.xaxis.set_major_locator(dt.MonthLocator(interval=1))
ax1.xaxis.set_major_formatter(hfmt)
subdir_name='ChargingScatter'
if not os.path.exists(subdir_name):
     os.umask(0) # unmask if necessary
     os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig('compareIonEFW.png')
plt.close(fig)
os.chdir('..')

