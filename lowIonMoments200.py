# lowIonMoments.py
#
#  for spacecraft charging moments


import os
import glob 
import numpy as np
from matplotlib import pyplot as plt
import datetime
import matplotlib.dates as dates
import spacepy.pybats.kyoto as spk
from spacepy import pycdf
import itertools as itert
import math
from numpy import ma
import pandas as pd
import pickle
# 
# SETTINGS
n1=0
n2=35 # 210 eV
PAs=11 # 11 pitch angle bins on HOPE
dateStart='20130201'           # starting date
#dateEnd='20130301'
dateEnd='20130701'
#dateEnd='20150401'             # ending date
nLbins=11                       # for the number of L shells, spaced by .25
LbinsArr=np.linspace(1.5, 4,nLbins)
nmlt_bins=48                   # 30 minute time resolution in MLT
MLTbinsArr=np.linspace(0.25, 23.75, nmlt_bins)
satellite=['A', 'B']
name_species=[ 'FPDO', 'FHEDO', 'FODO']
spe=['P']#,'He','O']
mass=[2]#, 4, 16]
#
# Data arrays

# prep the start times
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
monthCur=DT.month
for iSpe in range(len(name_species)):
 dfs={} # for a and b data storage
 dfsT={}
 for iSat in range(len(satellite)):
  fplasmaDensity=[]
  plasmaTime=[]
  DT=datetime.datetime.strptime(dateStart, '%Y%m%d') 
  while DT != endDt:
     date=datetime.datetime.strftime(DT, '%Y%m%d')
     DT=datetime.datetime.strptime(date, '%Y%m%d')
     # load in HOPE files
     os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+satellite[iSat])
     f=glob.glob('*'+date+'*')
     try:
         pyf=pycdf.CDF(f[0])
         #
         # put into pandas arrays
         energy=np.swapaxes(pyf['HOPE_ENERGY_Ion'][...],1,0)
         epoch=pyf['Epoch_Ion'][...]
         MLT=pyf['MLT_Ion'][...]
         L=pyf['L_Ion'][...]
         EnergyDelta=np.swapaxes(pyf['ENERGY_Ion_DELTA'][...],1,0)[n1:n2]
         FLUX=pyf[name_species[iSpe]][...]
         FLUXl=np.swapaxes(FLUX,1,0)[n1:n2]
         dataFrames={}
         kyoto=spk.fetch('kp', DT, DT)
         day=datetime.datetime.strftime(DT,'%d') # need to get day of month
         kp_arr=kyoto['kp'][(int(day)-1)*8: int(day)*8]
         kpC=0
         kpM=0
         nKP=np.zeros(1440)
         for iKP in range(1440):
              nKP[iKP]=kp_arr[kpM]
              kpC+=1
              if kpC == 180:
                 kpC=0
                 kpM+=1
              else:
                 continue

         rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
         rng=rt[::1]
         kyotodf2=pd.DataFrame({'kp':nKP}, index=rng)
         #kyotodf2=kyotodf['kp'].resample('5min',how='median').reindex(index=rng)
         COUNTS=np.swapaxes(pyf['Counts_'+spe[iSpe]+'_Omni'][...],1,0)[n1:n2]
         #
         # each dataFrame is set up to be by the energy channel number (so 0 = 0.9 eV)
         for iN in range(n1,n2):
            index=np.where(energy[iN] > 20)[0]
            energy[iN][index]=np.nan
            FLUXl[iN][index]=np.nan # bad energies
            COUNTS[iN][index]=np.nan # bad energies
            offFlux=np.where((FLUXl[iN] <= 0) & (FLUXl > 1e12))[0]
            FLUXl[iN][offFlux]=np.nan
            COUNTS[iN][offFlux]=np.nan
            dataFrames[str(iN)]=pd.DataFrame({'MLT':MLT, 'L':L, 'flux':FLUXl[iN], 'energy':energy[iN], 'delta':EnergyDelta[iN]}, index=epoch)
            # now resample the data frame
            #
            # NEED TO WRITE MY OWN CODE TO RESAMPLE FOR 5 MINUTES BUT FACTOR IN WEIGHTED AVERAGE
            def weightedAvg(window,df,weight,DT):
               # window = time window
               # df = data frame to avg over
               # weight = array of weights
               # DT = starting date time
               cutoff=1440/window # number of windows in a day by minutes
               newDF={'MLT':np.zeros(cutoff),'L':np.zeros(cutoff), 'flux':np.zeros(cutoff), 'energy':np.zeros(cutoff), 'delta':np.zeros(cutoff) }
               weight=np.array(weight)
               for iTime in range( cutoff-1):
                  sI=df.index.searchsorted(DT+iTime*datetime.timedelta(minutes=1))
                  eI=df.index.searchsorted(DT+(iTime+1)*datetime.timedelta(minutes=1))
                  # we have five minute window now get the weights
                  totalC=1.0*np.nansum(weight[sI:eI])
                  newDF['MLT'][iTime]=np.nansum(np.array(df['MLT'])[sI:eI]*weight[sI:eI])/totalC
                  newDF['L'][iTime]=np.nansum(np.array(df['L'])[sI:eI]*weight[sI:eI])/totalC
                  newDF['flux'][iTime]=np.nansum(np.array(df['flux'])[sI:eI]*weight[sI:eI])/totalC
                  newDF['energy'][iTime]=np.nansum(np.array(df['energy'])[sI:eI]*weight[sI:eI])/totalC
                  newDF['delta'][iTime]=np.nansum(np.array(df['delta'])[sI:eI]*weight[sI:eI])/totalC
               newDF['MLT'][-1]=df['MLT'][-1]
               newDF['L'][-1]=df['L'][-1]
               newDF['flux'][-1]=df['flux'][-1]
               newDF['energy'][-1]=df['energy'][-1]
               newDF['delta'][-1]=df['delta'][-1]
               timeArr=pd.date_range(DT,DT+datetime.timedelta(days=1), freq='1Min')[:-1]
               DF=pd.DataFrame(newDF, index=timeArr)
               # 5 minute weighted average
               return DF, timeArr

            temp=weightedAvg(1, dataFrames[str(iN)], COUNTS[iN], DT)
            dataFrames[str(iN)]=temp[0]
            tempTime=temp[1]
            #dataFrames[str(iN)]['kp']=kyotodf2 # add resampled kp onto the resampled flux
            #
            # now sum for plasma density
            Energy_Ev=dataFrames[str(iN)]['energy']*(1.602*1e-19)
            mass_adj=mass[iSpe]*(1.67*1e-27) 
            dataFrames[str(iN)]['density']=4.0*np.pi*(1.0/np.sqrt(2.0*Energy_Ev/mass_adj))*np.array(dataFrames[str(iN)]['delta'])*(1.0e-5)*np.array(dataFrames[str(iN)]['flux'])


         # sum the plasma density
         count=0
         for iN in range(n1,n2):
             if count ==0:
                         temp2=np.array(dataFrames[str(iN)]['density'])
                         temp2[np.isnan(temp2)]=0
                         densT=temp2
                         count=1
             else:
                         # add on to get total density
                         temp2=np.array(dataFrames[str(iN)]['density'])
                         temp2[np.isnan(temp2)]=0
                         densT=densT+temp2
         if len(tempTime) == len(densT):
           plasmaTime+=list(tempTime)
           fplasmaDensity+=list(densT)
                     
     except:
        print 'bad date: ' + str(date)
     yrCur=DT.year
     DT=DT+datetime.timedelta(days=1)

     # save by month
     if DT.month != monthCur:
     # save this file
        os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/PlasmaDensity200eV')
        with open(spe[iSpe]+'_sat=rbsp'+str(satellite[iSat])+'_'+str(yrCur)+'_'+str(monthCur), 'wb') as f:
              pickle.dump(fplasmaDensity, f)
        with open('time_'+spe[iSpe]+'_sat=rbsp'+str(satellite[iSat])+'_'+str(yrCur)+'_'+str(monthCur), 'wb') as f:
              pickle.dump(plasmaTime, f)          
       
        monthCur=DT.month
        # clean plasmadensity
        fplasmaDensity=[]
        plasmaTime=[]
  




         
         

 












