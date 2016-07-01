# electronfluxesbinned.py
#
# create a probability distribution function for e Fluxes
# and see how it compares with the charging distributions
#
# LKS, October 2015
#
# imports
import h5py
import numpy as np
from spacepy import pycdf
import os
import glob
import matplotlib.dates as dt
logBins=np.logspace(5,9,40)
sats=['A', 'B']
electronFlux=[[] for i in range(72)]
EFE=[[0 for j in range(40)] for i in range(72)]
for isat in range(len(sats)):
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/')
    #
    # empty lists to append
    os.chdir('HOPE_L3_'+sats[isat])
    #
    # get all the files
    listAll=glob.glob('*')
    #
    # loop
    for iCDF in range(len(listAll)):
     try:
        pyf=pycdf.CDF(listAll[iCDF])
        for i in range(72):
            # bin here
            Fluxes=np.array(np.swapaxes(pyf['FEDO'][...],1,0)[i])/1000.0* np.array(np.swapaxes(pyf['HOPE_ENERGY_Ele'][...], 1, 0)[i])
            #L=np.array(pyf['L_Ele'])
            #Lscreen=np.where(L>4)[0]
            for iF in range(39):
                temp=np.where((Fluxes>= logBins[iF]) & (Fluxes<logBins[iF+1]))[0]
                EFE[i][iF]=EFE[i][iF]+len(list(temp))            
     except:
         print listAll[iCDF]
    os.chdir('..')
#
# now bin
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
import pickle
pickle.dump(EFE, open('binnedelectronfluxes.p', 'wb'))
