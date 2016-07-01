# getEclipseTimes.py
#
# get the charging line from HOPE during eclipse
#
# LKS, September 2015
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
#
#
# start the date loop
dateStart='20130201'
dateEnd='20150401'
#
nLbins=11
nmlt_bins=48
endDT=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(dateStart, '%Y%m%d')
sat=['A','B']
Lbins=56
LbinsArr=np.linspace(1, 6.5,Lbins)
eFluxEclipse=[[] for i in range(72)]
scIonEclipse=[]
EclipseTimes=[]
for iSat in range(len(sat)):
     print "through a satellite"
     while DT != endDT:
  
        date=datetime.datetime.strftime(DT, '%Y%m%d')
        DT=datetime.datetime.strptime(date, '%Y%m%d')
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+sat[iSat])
        f=glob.glob('*'+date+'*')

        #
        try:
            pyf=pycdf.CDF(f[0])
        except:
            print "bad data"
            DT=DT+datetime.timedelta(days=1)
            continue
        #
        # resample EFW potential
        efwTime=pyf['epoch'][...]
        chargingFlagsTemp=pyf['flags_charging_bias_eclipse'][...]
        # identify any times of eclipse as 1 or 0
        chargingFlags=np.zeros(len(efwTime))
        for iTime in range(len(efwTime)):
            if any(np.array(chargingFlagsTemp[iTime]) >0):
                chargingFlags[iTime]=1
            else:
                chargingFlags[iTime]=0
        # great, now get the HOPE data
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
        epochEle=pyf['Epoch_Ele'][...]
        
        #
        # now get HOPE ion line
        counts=pyf['Counts_P_Omni'][...]
        epochI=pyf['Epoch_Ion'][...]
        HIE=pyf['HOPE_ENERGY_Ion'][...]
        Lion=pyf['L_Ion'][...]
        os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/BackgroundCounts')
        bgCounts=pickle.load(open('bgCounts.p','rb'))
        for iTime in range(len(epochI)):
            for iEN in range(72):
               idx = np.argmin(np.abs(LbinsArr - np.array(Lion)[iTime]))
               try:
                   counts[iTime][iEN]=counts[iTime][iEN]-bgCounts[iEN][idx]
               except(ValueError):
                   counts[iTime][iEN]=counts[iTime][iEN]
               #
        # extraction
        os.chdir('..')
        chargingLine=[]
        chargingTimes=[]
        for iTime in range(len(epochI)):
            # find the gradient
          check=0
          if Lion[iTime] > 4.5:
            temp=np.where(counts[iTime]>0)[0][0:5] # should be clear
            for iTemp in range(len(temp)):
                # is it above like 6 eV
                if (temp[iTemp]>10) and (temp[iTemp]<45):
                    # is it above 50 counts at hte ion line
                    try:
                     if (counts[iTime][temp[iTemp]]-counts[iTime][temp[iTemp+2]])>10:
                      # is it only a few thick under these same requirements
                     #  if temp[iTemp+3] > temp[iTemp]+3:
                        # this is a charging line
                        chargingLine.append(HIE[iTime][temp[iTemp]])
                        chargingTimes.append(epochI[iTime])
                        check=1
                    except(IndexError):
                        break
          if check!=1:
              chargingLine.append(0)
              chargingTimes.append(epochI[iTime])
        #
        # resample everything and then take times greater than 0 for
        # the charging flags
        
        epochEflux=pd.DatetimeIndex(epochEle)
        epochIflux=pd.DatetimeIndex(chargingTimes)
        epochPot=pd.DatetimeIndex(efwTime)

        rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
        Efluxes={}
        for iEN in range(72):
            temp=pd.DataFrame(energyFlux[iEN], index=epochEflux, columns=['Ele'])
            Efluxes[iEN]=np.array(temp['Ele'].resample('1min',how='median').reindex(index=rt, fill_value=np.nan))
        IonSCtemp=pd.DataFrame(chargingLine, index=epochIflux, columns=['Ion'])
        IonSC=np.array(IonSCtemp['Ion'].resample('1min',how='median').reindex(index=rt, fill_value=np.nan))
        ChargeFlagsTemp=pd.DataFrame(chargingFlags, index=epochPot, columns=['Flags'])
        ChargeFlags=np.array(ChargeFlagsTemp['Flags'].resample('1min',how='median').reindex(index=rt, fill_value=np.nan))

        eclipseTimes=np.where(ChargeFlags>0)[0]
        for iEN in range(72):
            eFluxEclipse[iEN]+=list(Efluxes[iEN][eclipseTimes])
        scIonEclipse+=list(IonSC[eclipseTimes])
        EclipseTimes+=list(rt[eclipseTimes])
        DT=DT+datetime.timedelta(days=1)
         
     date=dateStart
     DT=datetime.datetime.strptime(date, '%Y%m%d')        
                                    
        #
        # 
os.chdir('ScatterBoth')
with open('IonEclipseTimes.p','wb') as f:
     pickle.dump(EclipseTimes, f)
with open('IonLineEclipse.p','wb') as f:
     pickle.dump(scIonEclipse,f)
with open('EFluxEclipse.p','wb') as f:
     pickle.dump(eFluxEclipse,f)
     
