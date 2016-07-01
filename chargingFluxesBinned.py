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
EFE=[[0 for j in range(30)] for i in range(72)]


# eclipseEvents.py
#
# determine negative charging during eclipse vs not during
#
# LKS, October 2015
#
#
import h5py
import numpy as np
from spacepy import pycdf
import os
import glob
import datetime
import matplotlib.dates as dt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from spacepy import datamodel as dm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
import pandas as pd
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

#
# create a half circle for density maps
def dual_half_circle(center, radius, angle=90, ax=None, colors=('w','k'),
                     **kwargs):
    from matplotlib.patches import Wedge
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0],transform=ax.transData._b, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], transform=ax.transData._b,**kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]

iDate='20130201'
dt=datetime.datetime.strptime(iDate,'%Y%m%d')
endDate='20150401'
satellites=['A', 'B']
#
# MLT and L-Shell
nMLT=48
nL=24
MLTbins=np.linspace(0.25, 23.75, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
Lbins=np.linspace(1.5, 7.25, nL)
electronFluxes=[ [] for iEn in range(72)]
totalFluxes=[ [] for iEn in range(72)]
lenMLTL=[[ 0 for i in range(nL)] for j in range(nMLT)]
#
chargingArr=[]
RE=6378000

for iSat in range(len(satellites)):
 iDate='20130201'
 dt=datetime.datetime.strptime(iDate,'%Y%m%d')
 for isat in range(len(satellites)):
   while iDate != endDate:
    try:
       iDate=datetime.datetime.strftime(dt,'%Y%m%d')
       dt=datetime.datetime.strptime(iDate, '%Y%m%d')
       if (iDate=='20140428') or (iDate=='20140429'):
           # bad data!
           skip
       # get the EFW data
       os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+satellites[isat])
       f=glob.glob('*'+iDate+'*')
       pyf=pycdf.CDF(f[0])
       potential=-1*pyf['Vavg'][...]
       rng = pd.period_range(iDate,periods=1440, freq='T').to_timestamp()
       pEpoch=pd.DatetimeIndex(pyf['epoch'][...])
       df=pd.DataFrame(potential, index=pEpoch, columns=['potential'])
       Phi=np.array(df['potential'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
       negCharge=np.where(Phi<0)[0]

       #
       # now get the ephemeris data
       #eclipseFlag=np.zeros(len(rng))
       os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+satellites[isat])
       files=glob.glob('*'+iDate+'*')
       pyf=pycdf.CDF(files[0])
       eEpoch=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
       eFluxes=np.swapaxes(np.array(pyf['FEDO'])*np.array(pyf['HOPE_ENERGY_Ele']), 1, 0)

       for iEner in range(72):
           df=pd.DataFrame(eFluxes[iEner], index=eEpoch, columns=['Eflux'])
           eEF=np.array(df['Eflux'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
           electronFluxes[iEner]+=list(eEF[negCharge])
           totalFluxes[iEner]+=list(eEF)
       chargingArr+=list(Phi[negCharge])
       
    except:
           print iDate
           print len(chargingArr)
           print len(electronFluxes[10])
    dt=dt+datetime.timedelta(days=1)
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/')


os.chdir('/Users/loisks/Documents/Functions/')
import pickling as pick
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/')
pick.hdf5_data_save1(totalFluxes, 'AllEFluxes', 'ChargingElectronFluxes', 72)
pick.hdf5_data_save1(electronFluxes,'EFluxesE', 'ChargingElectronFluxes', 72)
#pick.hdf5_data_save1(chargingArr, 'EFWSC',  'ChargingElectronFluxes', 72)
#import pickle
#pickle.dump(electronFluxes, open('EFluxesE.p', 'wb'))
pickle.dump(chargingArr, open('EFWSC.p', 'wb'))
#pickle.dump(totalFluxes, open('AllEFluxes.p', 'wb'))


