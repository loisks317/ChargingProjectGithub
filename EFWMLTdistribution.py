# EFWMLTdistribution.py
#
# look at MLT distribution of EFW data
#
# LKS, August 2015
#
#
import h5py
import numpy as np
from spacepy import pycdf
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
#
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
    dataEFW[iStuff]=h5py.File('combined_EFW_L3_A_'+iStuff+'.h5', 'r')[iStuff]
os.chdir('..')

nMLT=48
nL=25
nMLAT=41
MLTbins=np.linspace(0.25, 23.75, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
MLATarr=[ [] for i in range(nMLAT)]
Lbins=np.linspace(0.25, 6.25, nL)
MLATbins=np.linspace(-20,10, nMLAT)

# set the right index
base=datetime.datetime.strptime('20150401', '%Y%m%d')
b2=dates.date2num(base)
diff=np.array(dataEFW['epoch'][:])-b2
minIndex=np.where(diff==np.min(np.abs(diff)))[0]

dataEFW['MLT']=dataEFW['MLT'][:minIndex]
dataEFW['MLAT']=dataEFW['MLAT'][:minIndex]
dataEFW['epoch']=dataEFW['epoch'][:minIndex]
eclipseTimes0=np.where(np.array(dataEFW['eclipseFlags'][0][:minIndex]) ==1)
eclipseTimes1=np.where(np.array(dataEFW['eclipseFlags'][1][:minIndex]) ==1)
eclipseTimes2=np.where(np.array(dataEFW['eclipseFlags'][2][:minIndex]) ==1)
#
# set the eclipse times
dataEFW['MLT'][eclipseTimes0]=np.nan
dataEFW['MLT'][eclipseTimes1]=np.nan
dataEFW['MLT'][eclipseTimes2]=np.nan
dataEFW['MLAT'][eclipseTimes0]=np.nan
dataEFW['MLAT'][eclipseTimes1]=np.nan
dataEFW['MLAT'][eclipseTimes2]=np.nan


# # screen for eclipse times
dataEFW['MLT']

for iMLT in range(nMLT):
    try:
        MLTindex=np.where((dataEFW['MLT'] >= iMLT*0.5) & (dataEFW['MLT'] < (iMLT+1)*0.5))[0]
        MLTarr[iMLT]=len(MLTindex)
    except:
        MLTarr[iMLT]=0
for iMLAT in range(nMLAT):
    try:
               MLATindex=np.where((dataEFW['MLAT'] >= iMLAT*1-10) & (dataEFW['MLAT'] < (iMLAT+1)-10))[0]
               MLATarr[iMLAT]=len(MLATindex)
    except:
               MLATarr[iMLAT]=0
# plot bar for each bin size
fig2=plt.figure()
ax2=fig2.add_subplot(111)
ax2.set_yscale('log')
ax2.set_xlim(0,24)
ax2.bar(MLTbins, MLTarr, 0.5, color='b')
ax2.set_title('MLT distribution')
font = {'family' : 'normal',
          'weight' : 'bold',
          'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
ax2.set_ylabel('Number of Occurrences', fontweight='bold')
ax2.set_xlabel('MLT ',fontweight='bold')
os.chdir('MLT_barplots')              
plt.savefig('MLT_all_distribution.pdf')
os.chdir('..')

# plot bar for each bin size
fig2=plt.figure()
ax2=fig2.add_subplot(111)
#ax2.set_yscale('linear')
ax2.set_xlim(-20,10)
ax2.bar(MLATbins, MLATarr, 1, color='lightseagreen')
ax2.set_title('MLAT distribution')
font = {'family' : 'normal',
          'weight' : 'bold',
          'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.25, right=0.8, top=0.85, bottom=0.15)
ax2.set_ylabel('Fraction of Total', fontweight='bold')
ax2.set_xlabel('MLAT ',fontweight='bold')
os.chdir('MLAT_barplots')              
plt.savefig('MLAT_all_distribution.pdf')
os.chdir('..')
