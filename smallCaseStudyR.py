# smallCaseStudyR.py
#
# try a select interval and see if there is a high correlation coefficient
# within that interval
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
from scipy.interpolate import interp1d
import scipy
import scipy.stats
#
#
date='20130301'
# data
hourGap=5
#
# load in EFW data
# open the files
os.chdir('EFW_L3_A')
f=glob.glob('*'+date+'*')
pyf=pycdf.CDF(f[0])
density=(list(pyf['density'][...]))
potential=-1*pyf['Vavg'][...]
epochE=dt.date2num(pyf['epoch'][...])
MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
#
# now get HOPE data
#
os.chdir('..')
os.chdir('HOPE_L3_A')
f=glob.glob('*'+date+'*')
pyf=pycdf.CDF(f[0])

electronFlux=pyf['FEDO'][...]
electronEnergy=pyf['HOPE_ENERGY_Ele'][...]
electronFlux=electronFlux*electronEnergy
epoch=pyf['Epoch_Ele'][...]
L=pyf['L_Ele'][...]
MLT=pyf['MLT_Ele'][...]
os.chdir('..')
mTime=dt.date2num(epoch)

# define a charging event
lowPotential=np.where(potential < -5)[0]
lowPotential[0]=-10
# now sequence
potentialEvents={}
count=-1
for i in range(1,len(lowPotential)):
    if lowPotential[i] != lowPotential[i-1]+1:
        # out of sequence
        count+=1
        potentialEvents[count]=[lowPotential[i]]
    else:
        potentialEvents[count].extend([lowPotential[i]])

# Now calculate a correlation coefficient in the ones left

# now test to make sure each event is long enough
for iCount in range(count):
    interpFlux=[[] for i in range(72)]
    if len(potentialEvents[count]) < 100:
        print "Too small of an event"
        potentialEvents.pop(count,0) # this should remove the item
        # doesn't change the count though, will need to account
        # for that later

#
# Now calculate a correlation coefficient in the ones left
    else:
        f=interp1d(epochE[potentialEvents[iCount]],potential[potentialEvents[iCount]])
        # now find where mTime is greater than epochE[startEvent]
        AdeleIsAwesome=np.where(mTime > epochE[potentialEvents[iCount]][0])[0][0]
        AdeleIsAwesomer=np.where(mTime > epochE[potentialEvents[iCount]][-1])[0][0]-1
        try:
                   tempC=f(mTime[AdeleIsAwesome:AdeleIsAwesomer])
        except(ValueError):
                   # gah fuck it
                   print "you fucked it up girl" 
                   mTime[mTime<epochE[potentialEvents[iCount]][0]]=epochE[potentialEvents[iCount]][0]
                   mTime[mTime>epochE[potentialEvents[iCount]][-1]]=epochE[potentialEvents[iCount]][-1]
                   tempC=f(mTime[AdeleIsAwesome:AdeleIsAwesomer])
        allPotential=list(tempC)
        for iEn in range(72):
            temp=np.swapaxes(electronFlux,1,0)[iEn][AdeleIsAwesome:AdeleIsAwesomer]
            interpFlux[iEn]+=list(temp)
            try:
                interpFlux[iEn]=np.array(interpFlux[iEn][~np.isnan(interpFlux[iEn])])
            except:
                interpFlux[iEn]=np.array(interpFlux[iEn])
            seriesB = np.transpose(np.log10(interpFlux[iEn]))
            seriesA = np.transpose(allPotential)
            ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
            tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
            spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
            print("For event count " + str(iCount) + " Kendall's tau for (Energy={0}) is {1} (p={2})".format(iEn,tau,p_value))
            print("Spearman's R for (Energy={0}) is {1} (p={2})".format(iEn,spear,spear_p))
            fig2=plt.figure()
            ax2=fig2.add_subplot(111)
            ax2.plot(seriesA, seriesB, '.', color='k')
            xx=np.arange(-200,0)
            ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
            plt.draw()
            font = {'family' : 'normal',
                      'weight' : 'bold',
                      'size'   : 22}
            plt.rc('font', **font)
            plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
            ax2.set_ylabel('Log(Energy Flux)', fontweight='bold')
            ax2.set_xlabel('$\phi$',fontweight='bold')
            ax2.set_ylim(9,12)
            ax2.set_xlim(-75,0)
            ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(-50, 12), xycoords='data', color='k')
            os.chdir('CorrelationPlots_SpecificEvent')
            plt.savefig("event="+str(iCount)+"_energy="+str(iEn)+'.pdf')
            plt.close()
            os.chdir('..')
