# integratedFluence.py
#
# see if there is a link between certain energy electron fluxes
# and spacecraft charging
# fluence over an hour/time scale with multiple correlation variables

#
#
#
# LKS, July 2015
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
# threshold
threshold=-0
maxVal=-200
timeWindow = 0.04 # approximately 1 hour
windowStr="1.0hours_extremeCharging"
scatter='on'
hitTest='off'

# do a -50 V + and - test 
#
# unique dates for less than 50V
uniqueDates=['20130108',
 '20130113',
 '20130114',
 '20130117',
 '20130120',
 '20130126',
 '20130202',
 '20130208',
 '20130213',
 '20130214',
 '20130216',
 '20130217',
 '20130220',
 '20130222',
 '20130223',
 '20130226',
 '20130301',
 '20130316',
 '20130317',
 '20130320',
 '20130321',
 '20130323',
 '20130327',
 '20130328',
 '20130329',
 '20130404',
 '20130407',
 '20130424',
 '20130426',
 '20130501']

#
# open up the files
# interpolate HOPE onto EFW s/c
# compare by energy

TotalFluence=[[] for i in range(72)]
allPotential=[]
IntegratedFluence=[]
#
# now load in HOPE data for this

for iDate in range(len(uniqueDates)):
       date=uniqueDates[iDate]
       # first get EFW
       os.chdir('EFW_L3_A')
       f=glob.glob('*'+uniqueDates[iDate]+'*')
       pyf=pycdf.CDF(f[0])
       density=(list(pyf['density'][...]))
       potential=-1*pyf['Vavg'][...]
       epochE=dt.date2num(pyf['epoch'][...])
       MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
       LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
       MLATE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
       # identify points of charging
       chargingIndex=np.where(potential < threshold)[0]
       chargingTimes=epochE[chargingIndex]
       charging=potential[chargingIndex]


       
       #
       # now get HOPE data
       os.chdir('..')
       os.chdir('HOPE_L3_A')
       f=glob.glob('*'+uniqueDates[iDate]+'*')
       pyf=pycdf.CDF(f[0])
       
       electronFlux=pyf['FEDO'][...]
       electronEnergy=pyf['HOPE_ENERGY_Ele'][...]
       electronFlux=electronFlux*electronEnergy
       epoch=pyf['Epoch_Ele'][...]
       L=pyf['L_Ele'][...]
       MLT=pyf['MLT_Ele'][...]
       os.chdir('..')
       mTime=dt.date2num(epoch)

       
  
      # so take the time window of previous fluence
       for iTime in range(len(chargingTimes)):
          if chargingTimes[iTime]-np.floor(chargingTimes[iTime])>timeWindow:
          #start=np.where(chargingTimes>=(chargingTimes[iTime]-timeWindow))[0][0]
           allPotential.append(charging[iTime])
           HopeStart=np.where(mTime>=(chargingTimes[iTime]-timeWindow))[0][0]
           HopeEnd=np.where(mTime>=chargingTimes[iTime])[0][0]
           tempSums=[]
           for iEn in range(72):
              fluTemp=np.swapaxes(electronFlux, 1, 0)[iEn][HopeStart:HopeEnd]
              # interp to the same number of points
              xvals=np.linspace(mTime[HopeStart], mTime[HopeEnd], 200)
              summedFluence=np.sum(np.interp(xvals, mTime[HopeStart:HopeEnd], fluTemp))
              TotalFluence[iEn].append(summedFluence)
              tempSums.append(summedFluence)
           IntegratedFluence.append(np.sum(np.array(tempSums)))
        
            
        
      

        # find correlation between the two
if scatter == 'on':
 # for iEn in range(72):
    seriesB = np.transpose(np.log10(IntegratedFluence))
    seriesA = np.transpose(allPotential)
    ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
    tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
    spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
    print("Kendall's tau for is {0} (p={1})".format(tau,p_value))
    print("Spearman's R for is {0} (p={1})".format(spear,spear_p))
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(seriesA, seriesB, '.', color='k')
    xx=np.arange(-200,0)
    ax2.set_ylim(14.0, 15.5)
    ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
    plt.draw()
    font = {'family' : 'normal',
                      'weight' : 'bold',
                      'size'   : 22}
    plt.rc('font', **font)
    plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
    ax2.set_ylabel('Log(Total Fluence)', fontweight='bold')
    ax2.set_xlabel('$\phi$',fontweight='bold')
    ax2.set_xlim(maxVal,threshold)
    ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(-100,13.5), xycoords='data', color='k')
    subdir_name='IntegratedFluence'
    if not os.path.exists(subdir_name):
               os.umask(0) # unmask if necessary
               os.makedirs(subdir_name, 0777) 
    os.chdir(subdir_name)
    plt.savefig('totalFlux_threshold='+str(threshold)+'_window='+windowStr+'.pdf')
    plt.close()
    os.chdir('..')

