# PotentialCorrectedDensity.py
#
# calculate plasma density with spacecraft potential from
# EFW factored in. Bin accordingly. 1 minute resamples.
#
# LKS - November, 2015 for GRL 
# 
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
# 
# SETTINGS
n1=0
n2=36
PAs=11 # 11 pitch angle bins on HOPE
dateStart='20130201'           # starting date
dateEnd='20130601'
nLbins=11                       # for the number of L shells, spaced by .25
LbinsArr=np.linspace(1.5, 4,nLbins)
nmlt_bins=48                   # 30 minute time resolution in MLT
MLTbinsArr=np.linspace(0.25, 23.75, nmlt_bins)
HopeBins=[9.84899998e-01,   1.19595003e+00,   1.33664989e+00,
         1.54769993e+00,   1.82909989e+00,   2.18085003e+00,
         2.53260016e+00,   2.95469975e+00,   3.37679982e+00,
         3.93959999e+00,   4.64310026e+00,   5.34659958e+00,
         6.26114988e+00,   7.31639957e+00,   8.51235008e+00,
         9.91934967e+00,   1.15373993e+01,   1.34368496e+01,
         1.56880503e+01,   1.82909985e+01,   2.12456989e+01,
         2.47631989e+01,   2.89138508e+01,   3.36273003e+01,
         3.91849518e+01,   4.57275009e+01,   5.32549515e+01,
         6.20486946e+01,   7.23197937e+01,   8.42089539e+01,
         9.81382446e+01,   1.14318748e+02,   1.33242905e+02,
         1.55262451e+02,   1.80869843e+02,   2.10768585e+02]
satellite=['A', 'B']
name_species=[ 'FPDO', 'FHEDO', 'FODO']
spe=['P','He','O']
mass=[2, 4, 16]
#
# Data arrays

# prep the start times
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
for iSpe in range(len(name_species)):
 dfs={} # for a and b data storage
 for iSat in range(len(satellite)):
  fplasmaDensity=[]
  fPotential=[]
  fTime=[]
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
         EnergyDelta=np.swapaxes(pyf['ENERGY_Ion_DELTA'][...],1,0)[n1:n2]
         FLUX=pyf[name_species[iSpe]][...]
         FLUXl=np.swapaxes(FLUX,1,0)[n1:n2]
         dataFrames={}
         Fluxdf={}
         Energydf={}
         Deltadf={}
         Densitydf=np.zeros(1440) # number of minutes in day
         #
         # get the spacecraft potential
         os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+satellite[iSat])
         f=glob.glob('*'+date+'*')
         pyfe=pycdf.CDF(f[0])
         potential=-1*pyfe['Vavg'][...]
         
         rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
         rng=rt[::1]
         #
         pEpoch=pd.DatetimeIndex(pyfe['epoch'][...])
         df=pd.DataFrame(potential, index=pEpoch, columns=['potential'])
         Phi=np.array(df['potential'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
         fPotential+=list(Phi)
         fTime+=list(rng)
         Phi=np.absolute(Phi)
         #
         COUNTS=np.swapaxes(pyf['Counts_'+spe[iSpe]+'_Omni'][...],1,0)[n1:n2]
         #
         # each dataFrame is set up to be by the energy channel number (so 0 = 0.9 eV)
         for iN in range(n1,n2):
            badDataScreen=np.where(energy[0]>20)[0] # perfect!
            #index=np.where(energy[iN] > 20)[0]
            #energy[iN][index]=np.nan
            #FLUXl[iN][index]=np.nan # bad energies
            #COUNTS[iN][index]=np.nan # bad energies
            offFlux=np.where((FLUXl[iN] < 0) & (FLUXl > 1e12))[0]
            phiHigh=np.where(Phi>10)[0]
            FLUXl[iN][offFlux]=np.nan
            COUNTS[iN][offFlux]=np.nan
            FLUXl[iN][phiHigh]=np.nan
            COUNTS[iN][phiHigh]=np.nan
            dataFrames[str(iN)]=pd.DataFrame({ 'flux':FLUXl[iN], 'energy':energy[iN], 'delta':EnergyDelta[iN]}, index=epoch)
            # now resample the data frame
            #

            Fluxdf[str(iN)]=np.array(dataFrames[str(iN)]['flux'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
            Deltadf[str(iN)]=np.array(dataFrames[str(iN)]['delta'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
            Energydf[str(iN)]=np.array(dataFrames[str(iN)]['energy'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
            # now everything is resampled
            # add resampled kp onto the resampled flux
            #
            # change the energy of the fluxes
            # add the s/c potential
            Energydf[str(iN)]=np.array(Energydf[str(iN)])+np.array(Phi)
            #
            # now sum for plasma density
            Energy_Ev=np.array(Energydf[str(iN)]*(1.602*1e-19))
            mass_adj=mass[iSpe]*(1.67*1e-27) 
            Densitydf=np.array(Densitydf)+4.0*np.pi*(1.0/np.sqrt(2.0*Energy_Ev/mass_adj))*np.array(Deltadf[str(iN)])*(1.0e-5)*np.array(Fluxdf[str(iN)])
         #print(Densitydf)
         fplasmaDensity+=list(np.array(Densitydf))
         
     except:
        print('bad date: ' + str(date))
     DT=DT+datetime.timedelta(days=1)
     #
     #
  try:
     fplasmaDensity=np.array(fplasmaDensity)
     fTime=np.array(fTime)
     fPotential=np.array(fPotential)
  except(TypeError):
     fplasmaDensity=fplasmaDensity
     print("error, probably no density from this day")
     print(str(date))

     # save this file
  os.chdir('/Users/loisks/Desktop/ResearchProjects/HOPEPlasmaDensities/')
  #
  #
  import pickle
  with open('1_210eV_corrected_plasma_density_'+str(n1)+'_'+str(n2)+'_sat='+satellite[iSat]+'_species='+spe[iSpe], 'wb') as f1:
      pickle.dump(fplasmaDensity, f1)
  with open('1_210eV_corrected_time_'+str(n1)+'_'+str(n2)+'_sat='+satellite[iSat]+'_species='+spe[iSpe], 'wb') as f2:
      pickle.dump(fTime,f2)
  with open('1_210eV_potential_'+str(n1)+'_'+str(n2)+'_sat='+satellite[iSat]+'_species='+spe[iSpe], 'wb') as f3:
      pickle.dump(fPotential, f3)
      



   
         
         

        
