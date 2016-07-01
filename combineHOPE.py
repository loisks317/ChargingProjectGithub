# combineHOPE.py
#
# this file combines the HOPE pitch angle resolved data from Jan 2013 to
# to July 2015 to make it more accessible
#
# LKS, July 2015
#
# imports
import h5py
import numpy as np
from spacepy import pycdf
import os
import glob
import matplotlib.dates as dt
#
# satellites
sats=['HOPE_L3_A', 'HOPE_L3_B']
for isat in sats:
    os.chdir(isat)
    #
    # empty lists to append
    electronFlux=[[] for i in range(72)]
    electronFluxPA=[ [ [] for i in range(11)] for j in range(72) ]
    electronEnergy=[ [] for i in range(72)]
    epoch=[]
    MLT=[]
    L=[]
    #
    # get all the files
    listAll=glob.glob('*')
    #
    # loop
    for iCDF in range(len(listAll)):
        pyf=pycdf.CDF(listAll[iCDF])
        for i in range(72):
            for j in range(11):
                electronFluxPA[i][j]+=list(np.swapaxes(pyf['FEDU'][...],2,0)[i][j])
            electronFlux[i]+=list(np.swapaxes(pyf['FEDO'][...],1,0)[i])
            electronEnergy[i]+=list(np.swapaxes(pyf['HOPE_ENERGY_Ele'][...],1,0)[i])
        epoch+=list(pyf['Epoch_Ele'][...])
        L+=list(pyf['L_Ele'][...])
        MLT+=list(pyf['MLT_Ele'][...])
    os.chdir('..')
    if not os.path.exists('Combined_HOPE'):
        os.umask(0) # unmask if necessary
        os.makedirs('Combined_HOPE', 0777) 
    os.chdir('Combined_HOPE')
    epoch=dt.date2num(epoch)
    var=['electronFluxPA', 'electronFlux', 'electronEnergy', 'epoch', 'L', 'MLT']
    dataVar=[electronFluxPA, electronFlux, electronEnergy, epoch, L, MLT]
    for iFile in range(len(dataVar)):
        FILE='combined_'+isat+'_'+var[iFile]+'.h5'
        f=h5py.File(FILE, "w")
        f.create_dataset(var[iFile], data=dataVar[iFile])
        f.close()
    os.chdir('..')
