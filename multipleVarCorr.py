# multipleVarCorr.py

# correlate multiple variables
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
#
# threshold
threshold=-50
maxVal=-200
scatter='on'



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
#
# now load in HOPE data for this

for n1 in range(71):
  for n2 in range(71):
   interpFlux=[[] for i in range(72)]
   allPotential=[]
   allMLT=[]
   allL=[]
   allMLAT=[]
   TotalFluence=[[] for i in range(72)]
   allPotential=[]
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
       #MLTEchar=MLTE[chargingIndex]
       #LEchar=LE[chargingIndex]
       #MLATEchar=MLTE[chargingIndex]
       
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
       f=interp1d(epochE,potential)
       #fMLT=interp1d(epochE,MLTE)
       # continue to work on this here tomorrow
       try:
           temp=f(mTime)
       except(ValueError):
           mTime[mTime<epochE[0]]=epochE[0]
           mTime[mTime>epochE[-1]]=epochE[-1]
           temp=f(mTime)
       null=np.where(temp<threshold)[0] # if it's below the threshold, don't include
       temp[null]=np.nan
       allPotential+=list(temp)
       #allMLT+=list(MLTEchar)
       for iEn in range(72):
           temp=np.swapaxes(electronFlux,1,0)[iEn]
           temp[null]=np.nan
           interpFlux[iEn]+=list(temp)
        

        # find correlation between the two
   if scatter == 'on':
        allPotential=np.array(allPotential)[~np.isnan(allPotential)]
        fluxLow=np.log10(np.array(interpFlux[n1])[~np.isnan(interpFlux[n1])])
        fluxHigh=np.log10(np.array(interpFlux[n2])[~np.isnan(interpFlux[n2])])
        rx1=scipy.stats.spearmanr(fluxLow, allPotential)[0]
        rx2=scipy.stats.spearmanr(fluxHigh, allPotential)[0]
        rx3=scipy.stats.spearmanr(fluxLow,fluxHigh)[0]
        R=np.sqrt((rx1**2 + rx2**2 - 2*rx1*rx2*rx3)/(1-rx3**2))
        print 'n1: ' + str(n1)+ ' n2: ' + str(n2)+ ' r = '+ str(R)
                    
#                     
#   fig2=plt.figure()
#   ax2=fig2.add_subplot(111)
#   ax2.plot(seriesA, seriesB, '.', color='k')
##   xx=np.arange(-10,10)
##   ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
#   plt.draw()
#   font = {'family' : 'normal',
#                     'weight' : 'bold',
#                     'size'   : 22}
#   plt.rc('font', **font)
#   plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
#   ax2.set_ylabel('Log(Energy Flux)', fontweight='bold')
#   ax2.set_xlabel('$\phi$',fontweight='bold')
#   ax2.set_ylim(9,12)
#   ax2.set_xlim(-10,10)
#   ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(10, 12), xycoords='data', color='k')
#   os.chdir('CorrelationPlots')
#   plt.savefig("lowEnergy="+str(iEn)+'_threshold='+str(threshold)+'.pdf')
#   plt.close()
#   os.chdir('..')

  
