c# trending.py
#
# see if there is a link between certain energy electron fluxes
# and spacecraft charging
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
import pandas as pd

#
# threshold
threshold=-10
scatter='on'
hitTest='off'
#

#
# open up the files
# interpolate HOPE onto EFW s/c
# compare by energy

interpFlux=[[] for i in range(72)]
allPotential=[]
ExcludeDates=['20130206', '20130227','20130314','20130321', '20130327','20130328','20130407', '20130417', '20130508']
dateStart='20130201'
dateEnd='20130601'
allTime=[]
dt1=datetime.datetime.strptime(dateStart, '%Y%m%d')
date=dateStart
DT=datetime.datetime.strptime(date, '%Y%m%d')
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
#
# now load in HOPE data for this

while DT != endDt:
       # approximately from january to may
 date=datetime.datetime.strftime(DT, '%Y%m%d')
 DT=datetime.datetime.strptime(date, '%Y%m%d')
 check=0
 for iDate in ExcludeDates:
     if iDate == date:
      check=1
 if check == 0:

       # first get EFW
       os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_A')
       f=glob.glob('*'+date+'*')
       try:
        pyf=pycdf.CDF(f[0])
       except():
              continue
       density=(list(pyf['density'][...]))
       potential=-1*pyf['Vavg'][...]
       epochE=dt.date2num(pyf['epoch'][...])
       MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
       LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
       MLATE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
       #
       # now get HOPE data
       os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_A')
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


       # resample
       HOPEDF={}
       for iChan in range(72):
          temp=pd.DataFrame(np.swapaxes(electronFlux,1,0)[iChan], index=epoch, columns=['flux'])
          HOPEDF[iChan]=temp['flux'].resample('1Min')
       temp=pd.DataFrame(potential, index=dt.num2date(epochE), columns=['potential'])
       EFWDF=temp['potential'].resample('1Min')
       

       # find where EFW potential is negative
       neg=np.where(np.array(EFWDF)<threshold)[0]

       allPotential+=list(np.array(EFWDF[neg]))
       allTime+=list(np.array(EFWDF.index[neg]))
       
       for iEn in range(72):
           interpFlux[iEn]+=list(np.array(HOPEDF[iChan][neg]))
       DT=DT+datetime.timedelta(days=1)
 else:
        DT=DT+datetime.timedelta(days=1)

        # find correlation between the two
if scatter == 'on':
  allPotential=np.array(allPotential)[~np.isnan(allPotential)]
  for iEn in range(72):
    interpFlux[iEn]=np.array(interpFlux[iEn])[~np.isnan(interpFlux[iEn])]
    seriesB = np.transpose(np.log10(interpFlux[iEn]))
    identify=np.where(np.array(seriesB)<10)[0]
    seriesA = np.transpose(allPotential)
    print np.array(allTime)[identify]
    ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
    tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
    spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
    print("Kendall's tau for (Energy={0}) is {1} (p={2})".format(iEn,tau,p_value))
    print("Spearman's R for (Energy={0}) is {1} (p={2})".format(iEn,spear,spear_p))
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(seriesA, seriesB, '.', color='k')
    #xx=np.arange(-10,10)
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
    ax2.set_xlim(0,-200)
    ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(10, 12), xycoords='data', color='k')
    os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/ScatterCaseStudies')
    plt.savefig("lowEnergy="+str(iEn)+'_threshold='+str(threshold)+'.pdf')
    plt.close()
    os.chdir('..')

#
# Hit Test
# set threshold for potential and for electron fluxes for each energy channel
#
fluxThresh=1e11
#potentialThresh=-80
if hitTest == 'on':
  allPotential=np.array(allPotential)[~np.isnan(allPotential)]
  hitRatio=[ [] for i in range(72)]
  for iEn in range(72):
    interpFlux[iEn]=np.array(interpFlux[iEn])[~np.isnan(interpFlux[iEn])]
    for iHit in range(len(interpFlux[iEn])):
        # 1 = double yes
        if (interpFlux[iEn][iHit] >= fluxThresh) & (allPotential[iHit] <= potentialThresh):
            hitRatio[iEn].append(1)
        # 2 = flux yes, potential no
        if (interpFlux[iEn][iHit] >= fluxThresh) & (allPotential[iHit] >= potentialThresh):
            hitRatio[iEn].append(2)
        # 3 = flux no, potential yes
        if (interpFlux[iEn][iHit] <= fluxThresh) & (allPotential[iHit] <= potentialThresh):
            hitRatio[iEn].append(3)        
        # 4 = flux no, potential no
        if (interpFlux[iEn][iHit] <= fluxThresh) & (allPotential[iHit] >= potentialThresh):
            hitRatio[iEn].append(4)
    # percentages
    hitRatio[iEn]=np.array(hitRatio[iEn])
    yesyes=1.0*len(np.where(hitRatio[iEn]==1)[0])/len(hitRatio[iEn])
    yesno=1.0*len(np.where(hitRatio[iEn]==2)[0])/len(hitRatio[iEn])
    noyes=1.0*len(np.where(hitRatio[iEn]==3)[0])/len(hitRatio[iEn])
    nono=1.0*len(np.where(hitRatio[iEn]==4)[0])/len(hitRatio[iEn])
    print "for energy="+str(iEn)
    print "High Flux, High Charging: " + str(yesyes)
    print "High Flux, Low Charging: " + str(yesno)
    print "Low Flux, High Charging: " + str(noyes)
    print "Low Flux, Low Charging: " + str(nono)
    
