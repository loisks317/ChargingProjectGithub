# HOPEMLTdistribution.py
#
# HOPE's MLT distribution
#
# LKS, August 2015
#
# get the HOPE ephemeris data
#
#
import numpy as np
from spacepy import pycdf
import glob
from matplotlib import pyplot as plt
import os
import datetime
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
import matplotlib.collections as collections
import matplotlib.dates as dt
import spacepy.datamodel as dm
#
# define parameters
startDate='20130101'
endDate='20150501'
dt=datetime.datetime.strptime(startDate,'%Y%m%d')
dtEnd=datetime.datetime.strptime(endDate,'%Y%m%d')
sats=['a','b']
cSats=['A', 'B']
#
# array presets
nMLT=48
nL=25
nMLAT=41
MLTbins=np.linspace(0.25, 23.75, nMLT)
MLTarr=[ 0 for i in range(nMLT)]
Larr=[ 0 for i in range(nL)]
MLATarr=[ 0 for i in range(nMLAT)]
Lbins=np.linspace(0.25, 6.25, nL)
MLATbins=np.linspace(-20,20, nMLAT)
#
#
while dt != dtEnd:
 for iSat in range(len(sats)):
    date=datetime.datetime.strftime(dt, '%Y%m%d')
    fem='rbsp'+sats[iSat]+'_def_MagEphem_OP77Q_'+date+'*'
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMPHEMERIS_'+cSats[iSat])
    files=glob.glob(fem)[0]
    try:
        pyfem=dm.readJSONheadedASCII(files)



        
    MLT=pyfem['CDMAG_MLT']
    MLT[MLT<0]=np.nan
    L=pyfem['L']; L[L<0]=np.nan
    m_L = np.ma.masked_array(L,np.isnan(L))
    L=np.mean(m_L, axis=1)
    MLAT=pyfem['CDMAG_MLAT'] # centered dipole MLAT                             
    MLAT[MLAT<-1000]=np.nan
    # add to arrays
    for iMLT in range(nMLT):
     try:
        MLTindex=np.where((MLT >= iMLT*0.5) & (MLT < (iMLT+1)*0.5))[0]
        MLTarr[iMLT]+=len(MLTindex)
     except:
        MLTarr[iMLT]+=0
    for iMLAT in range(nMLAT):
     try:
               MLATindex=np.where((MLAT >= iMLAT*1-10) & (MLAT < (iMLAT+1)-10))[0]
               MLATarr[iMLAT]+=len(MLATindex)
     except:
               MLATarr[iMLAT]+=0
    except:
        print "Didn't Work :( on " + str(datetime.datetime.strftime(dt,'%Y%m%d')) 
    dt=dt+datetime.timedelta(days=1)
# plot bar for each bin size
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
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
os.chdir('DistributionPlots')              
plt.savefig('MLT_all_distribution.pdf')
os.chdir('..')

# plot bar for each bin size
fig2=plt.figure()
ax2=fig2.add_subplot(111)
ax2.set_yscale('log')
ax2.set_xlim(-25,15)
ax2.bar(MLATbins, MLATarr, 1, color='r')
ax2.set_title('MLAT distribution')
font = {'family' : 'normal',
          'weight' : 'bold',
          'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
ax2.set_ylabel('Number of Occurrences', fontweight='bold')
ax2.set_xlabel('MLAT ',fontweight='bold')
os.chdir('DistributionPlots')              
plt.savefig('MLAT_all_distribution.pdf')
os.chdir('..')
