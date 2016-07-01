# AEpotential.py
#
# compare AE with s/c potential on a long term basis
# choose -25 V or less events 
#
# LKS, August 2015
#
# imports
import h5py
import numpy as np
from spacepy import pycdf
from spacepy import datamodel as dm
import os
import glob
import datetime
import matplotlib.dates as dates
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.interpolate import interp1d
import scipy
import scipy.stats
#
# start the date loop
dateStart='20130101'
dateEnd='20141231'
# data
hourGap=5
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
sat=['A', 'B']
aeList=[]
alList=[]
aeTimeList=[]
potentialList=[]
dateArr=[]
threshold=-100
#
# get all the AE data
# first get 2013
os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/AEindex/AE_121314/2013')
f=glob.glob('*')
for iFile in range(len(f)):
    data=dm.readJSONheadedASCII(f[iFile])
    aeList+=list(data['AE'])
    aeTimeList+=list(data['DateTime'])
    alList+=list(data['AL'])
os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/AEindex/AE_121314/2014')
f=glob.glob('*')
for iFile in range(len(f)):
    data=dm.readJSONheadedASCII(f[iFile])
    aeList+=list(data['AE'])
    aeTimeList+=list(data['DateTime'])
    alList+=list(data['AL'])
os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/')
#
AETimes=np.zeros(len(aeTimeList))
# now have complete lists of times and indexes
for iTime in range(len(aeTimeList)):
    AETimes[iTime]=dates.date2num(datetime.datetime.strptime(aeTimeList[iTime], '%Y-%m-%dT%H:%M:%S'))
#
# now have all the AE stuff
# get the s/c potential stuf
potentialLow=[]
potentialTimes=[]
sats=['A', 'B']
while DT != endDt:
  date=datetime.datetime.strftime(DT, '%Y%m%d')
  DT=datetime.datetime.strptime(date, '%Y%m%d')
  for iSat in range(len(sats)):
    os.chdir('EFW_L3_'+sats[iSat])
    f=glob.glob('*'+date+'*')
    try:
        pyf=pycdf.CDF(f[0])
    except:
        continue
    potential=-1.0*np.array(pyf['Vavg'][...])
    epochE=dates.date2num(pyf['epoch'][...])
    lowP=np.where(potential < threshold)[0]
    try:
        potentialLow+=list(potential[lowP])
        potentialTimes+=list(epochE[lowP])
    except:
        potentialLow+=[]
        potentialTimes+=[]
    os.chdir('..')
  DT=DT+datetime.timedelta(days=1)

# function to find nearest points
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx # return the index

AElow=np.zeros(len(potentialTimes))
ALlow=np.zeros(len(potentialTimes))
print "at the Find Nearest"
print "total points: " + str(len(potentialTimes))
# find nearest points in AE
for iTime in range(len(potentialTimes)):
    nearest=find_nearest(AETimes, potentialTimes[iTime])
    AElow[iTime]=aeList[nearest]
    ALlow[iTime]=alList[nearest]

#
# now correlate these two things
seriesB = np.transpose(AElow)
seriesA = np.transpose(potentialLow)
ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
print("Spearman's R is {0} (p={1})".format(spear,spear_p))
fig2=plt.figure()
ax2=fig2.add_subplot(111)
ax2.plot(seriesA, seriesB, '.', color='k')
xx=np.arange(-200,-25)
ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
plt.draw()
font = {'family' : 'normal',
           'weight' : 'bold',
           'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
ax2.set_ylabel('AE', fontweight='bold')
ax2.set_xlabel('$\phi$',fontweight='bold')
ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(-50,500 ), xycoords='data', color='k')
subdir_name='SubstormCorr'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)
plt.savefig('AEcorr_'+str(threshold)+'.pdf')
plt.close()
os.chdir('..')

#
# now correlate these two things
seriesB = np.transpose(ALlow)
seriesA = np.transpose(potentialLow)
ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
print("Spearman's R is {0} (p={1})".format(spear,spear_p))
fig2=plt.figure()
ax2=fig2.add_subplot(111)
ax2.plot(seriesA, seriesB, '.', color='k')
xx=np.arange(-200,-25)
ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
plt.draw()
font = {'family' : 'normal',
           'weight' : 'bold',
           'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
ax2.set_ylabel('AL', fontweight='bold')
ax2.set_xlabel('$\phi$',fontweight='bold')
ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(-50,500 ), xycoords='data', color='k')
subdir_name='SubstormCorr'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)
plt.savefig('ALcorr_'+str(threshold)+'.pdf')
plt.close()
os.chdir('..')
