# combineEFW.py
#
# this file combines the spacecraft potential from EFW from Jan 2013 to
# July 2015 to make it more accessible
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
#
# satellites
sats=['EFW_L3_A', 'EFW_L3_B']
for isat in sats:
    os.chdir(isat)
    #
    # empty lists to append
    potential=[]
    density=[]
    epoch=[]
    MLT=[]
    L=[]
    MLAT=[]
    eclipseFlags=[[] for i in range(3)]
    #
    # get all the files
    listAll=glob.glob('*')
    #
    # loop
    for iCDF in range(len(listAll)):
        pyf=pycdf.CDF(listAll[iCDF])
        potential+=(list(pyf['Vavg'][...]))
        density+=(list(pyf['density'][...]))
        epoch+=(list(pyf['epoch'][...]))
        MLT+=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
        L+=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
        MLAT+=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
        eclipseFlags[0]+=list(pyf['flags_charging_bias_eclipse'][...][:,0])
        eclipseFlags[1]+=list(pyf['flags_charging_bias_eclipse'][...][:,1])
        eclipseFlags[2]+=list(pyf['flags_charging_bias_eclipse'][...][:,2])
        

    os.chdir('..')
    if not os.path.exists('Combined_EFW'):
        os.umask(0) # unmask if necessary
        os.makedirs('Combined_EFW', 0777) 
    os.chdir('Combined_EFW')
    var=['potential', 'density', 'epoch', 'MLT', 'L', 'MLAT', 'eclipseFlags']
    epoch=dt.date2num(epoch)
    dataVar=[potential, density, epoch, MLT, L, MLAT, eclipseFlags]
    for iFile in range(len(dataVar)):
        FILE='combined_'+isat+'_'+var[iFile]+'.h5'
        f=h5py.File(FILE, "w")  
        f.create_dataset(var[iFile], data=dataVar[iFile])
        f.close()
    os.chdir('..')
    
