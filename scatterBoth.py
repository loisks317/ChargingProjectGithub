# scatterBoth.py
#
# scatter electron fluxes with ion line and s/c potential from EFW
# to see if there is an 'arm' in either of the datasets
#
# LKS, September 2015
#
#
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
from scipy.interpolate import interp1d
import scipy
import pickle
import scipy.stats
import pandas as pd
# start the date loop
dateStart='20130201'
dateEnd='20150401'
#dateEnd='20130214'
EnBins=72
Lbins=56
LbinsArr=np.linspace(1, 6.5,Lbins)
endDT=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(dateStart, '%Y%m%d')
sat=['A','B']
EI=[]
IonTimes=[]
EFWTimes=[]
EFWSC=[]
EFluxesIon=[ [] for i in range(72)]
EFluxesE=[ [] for i in range(72)]
for iSat in range(len(sat)):
     print "through a satellite"
     while DT != endDT:
        date=datetime.datetime.strftime(DT, '%Y%m%d')
        DT=datetime.datetime.strptime(date, '%Y%m%d')
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+sat[iSat])
        f=glob.glob('*'+date+'*')
        try:
            pyf=pycdf.CDF(f[0])
        except:
            print "bad data"
            DT=DT+datetime.timedelta(days=1)
            continue
        #
        # resample EFW potential
        efwTime=pd.DatetimeIndex(pyf['epoch'][...])
        dfepotential=pd.DataFrame(pyf['Vavg'][...], index=efwTime, columns=['potential'])
        rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
        rng=rt[::1]
        EFWpotential=dfepotential['potential'].resample('1min', how='median').reindex(index=rng, fill_value=np.nan)
        #
        # get HOPE electron fluxes
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+sat[iSat])
        f=glob.glob('*'+date+'*')
        try:
            pyf=pycdf.CDF(f[0])
        except:
              DT=DT+datetime.timedelta(days=1)
              print 'bad data'
              continue
        electronFlux=np.swapaxes(pyf['FEDO'][...],1,0)
        electronEnergy=np.swapaxes(pyf['HOPE_ENERGY_Ele'][...], 1, 0)
        energyFlux=np.array(electronFlux)*np.array(electronEnergy)
        epoch=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
        EF={}
        EFlux={}
        for iEn in range(len(energyFlux)):
            EF[iEn]=pd.DataFrame(energyFlux[iEn], index=epoch, columns=['EF'])
            EFlux[iEn]=EF[iEn]['EF'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)

        #
        # now get HOPE ion line
        counts=pyf['Counts_P_Omni'][...]
        epochI=pyf['Epoch_Ion'][...]
        HIE=pyf['HOPE_ENERGY_Ion'][...]
        Lion=pyf['L_Ion'][...]
        #
        # now get the background data
#        os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/BackgroundCounts')
#        bgCounts=pickle.load(open('bgCounts.p','rb'))
#        for iTime in range(len(epochI)):
#            for iEN in range(72):
#               idx = np.argmin(np.abs(LbinsArr - np.array(Lion)[iTime]))
#               try:
#                   counts[iTime][iEN]=counts[iTime][iEN]-bgCounts[iEN][idx]
#               except(ValueError):
#                   counts[iTime][iEN]=counts[iTime][iEN]
        #
        # extraction
#        os.chdir('..')
#        chargingLine=[]
#        chargingTimes=[]
#        for iTime in range(len(epochI)):
#            # find the gradient
#          check=0
#          if Lion[iTime] > 4.5:
#            temp=np.where(counts[iTime]>0)[0][0:5] # should be clear
#            for iTemp in range(len(temp)):
#                # is it above like 6 eV
#                if (temp[iTemp]>10) and (temp[iTemp]<45):
#                    # is it above 50 counts at hte ion line
#                    try:
#                     if (counts[iTime][temp[iTemp]]-counts[iTime][temp[iTemp+2]])>10:
#                      # is it only a few thick under these same requirements
#                     #  if temp[iTemp+3] > temp[iTemp]+3:
#                        # this is a charging line
#                        chargingLine.append(HIE[iTime][temp[iTemp]])
#                        chargingTimes.append(epochI[iTime])
#                        check=1
#                    except(IndexError):
#                        break
#          if check!=1:
#              chargingLine.append(0)
#              chargingTimes.append(epochI[iTime])            
                         
         #
         # resample the Ion Line data
        #epochIon=pd.DatetimeIndex(chargingTimes)
        #IonLinet=pd.DataFrame(chargingLine, index=epochIon, columns=['Ion'])
        #IonLine=IonLinet['Ion'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
        #IonLine=np.array(IonLine)
        #
        # now I have everything, add to big array
        for iEn in range(72):
            EFluxesE[iEn]+=list(EFlux[iEn][EFWpotential>1])
            #EFluxesIon[iEn]+=list(EFlux[iEn][IonLine>1])
        #IonTimes+=list(np.array(rt)[IonLine>1])
        #EI+=list(IonLine[IonLine>1])
        EFWTimes+=list(rt[EFWpotential>1])
        EFWSC+=list(EFWpotential[EFWpotential>1])
        print date
        DT=DT+datetime.timedelta(days=1)
         
     date=dateStart
     DT=datetime.datetime.strptime(date, '%Y%m%d')
# loop it
print "now plotting"
os.chdir('ScatterBoth')
#
# pickle the times here
#with open('IonTimes.p','wb') as f:
#     pickle.dump(IonTimes, f)
#with open('IonLine.p','wb') as f:
#     pickle.dump(EI,f)
with open('EfluxesSC.', 'wb') as f:
     pickle.dump(EFluxesE,f)
with open('EFWTimes.p','wb') as f:
     pickle.dump(EFWTimes, f)
with open('EFWSC.p','wb') as f:
         pickle.dump(EFWSC,f)
for iEn in range(72):
#    fig2=plt.figure()
#    ax2=fig2.add_subplot(111)
#    ax2.plot(EFluxesIon[iEn],EI, '.', color='k')
#    ax2.set_yscale('log')
#    ax2.set_xscale('log')
#    ax2.set_xlim(1e8,1e12)
#    ax2.set_ylim(1e0,1e4)
#    plt.draw()
#    font = {'family' : 'normal',
#                      'weight' : 'bold',
#                      'size'   : 22}
#    plt.rc('font', **font)
#    plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
#    ax2.set_ylabel('Spacecraft Potential [V]', fontweight='bold')
#    ax2.set_xlabel('Energy Flux',fontweight='bold')
#    subdir_name='ScatterBoth'
#    if not os.path.exists(subdir_name):
#               os.umask(0) # unmask if necessary
#               os.makedirs(subdir_name, 0777) 
#    os.chdir(subdir_name)
#    plt.savefig('IonLine_Energy='+str(electronEnergy[iEn][0])+'_eV.png')
#    plt.close()

    #
    # plot EFW stats
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(EFluxesE[iEn],EFWSC, '.', color='k')
    ax2.set_xlim(1e9,1e12)
    ax2.set_ylim(1e1,1e3)
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    plt.draw()
    font = {'family' : 'normal',
                      'weight' : 'bold',
                      'size'   : 22}
    plt.rc('font', **font)
    plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
    ax2.set_ylabel('Spacecraft Potential [V]', fontweight='bold')
    ax2.set_xlabel('Energy Flux',fontweight='bold')
    subdir_name='ScatterBoth'
    if not os.path.exists(subdir_name):
               os.umask(0) # unmask if necessary
               os.makedirs(subdir_name, 0777) 
    os.chdir(subdir_name)
    plt.savefig('EFWSC_Energy='+str(electronEnergy[iEn][0])+'_eV.png')
    plt.close()
    os.chdir('..')
