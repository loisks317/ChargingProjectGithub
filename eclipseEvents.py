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
endDate='20150501'
satellites=['A', 'B']
#
# MLT and L-Shell
nMLT=48
nL=24
MLTbins=np.linspace(0.25, 23.75, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
Lbins=np.linspace(1.5, 7.25, nL)
MLTLBins=[[ [] for i in range(nL)] for j in range(nMLT)]
lenMLTL=[[ 0 for i in range(nL)] for j in range(nMLT)]
#
chargingArr=[]
chargingTimes=[]
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
       Phi=df['potential'].resample('1min', how='median').reindex(rng, fill_value=np.nan)
       #
       # now get the ephemeris data
       eclipseFlag=np.zeros(len(rng))
       os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMPHEMERIS_'+satellites[isat])
       files=glob.glob('*'+iDate+'*')
       pyfem=dm.readJSONheadedASCII(files[0])
       gseCoords=np.swapaxes(pyfem["Rgse"], 1, 0)
       for iGSE in range(len(rng)):
          # subtract out an RE to get difference above the surface
          if ((gseCoords[1][iGSE])**2 + (gseCoords[2][iGSE])**2) < 1:
              if gseCoords[0][iGSE] < 0:
                  # satellite is in Earth's shadow
                  eclipseFlag[iGSE]=1
                  print 'eclipse flag'
       LShell=np.swapaxes(pyfem['L'],1,0)[0] # weird multiple rows with L
       MLT=pyfem['CDMAG_MLT']

      # now let's get dates and times where there is charging and it's
      # not in eclipse
       keyTimes=np.where((np.array(Phi) < 0) & (eclipseFlag != 1))[0]
       chargingYa=np.array(Phi[keyTimes])
       LYa=np.array(LShell[keyTimes])
       MLTYa=np.array(MLT[keyTimes])
       for iMLT in range(nMLT):
               for iL in range(nL):
                   temp=np.where((MLTYa >= MLTbins[iMLT]-0.25) & (MLTYa < MLTbins[iMLT]+0.25))[0]
                   temp2=np.where((LYa >= Lbins[iL]-.125) & (LYa < Lbins[iL]+0.125))[0]
                   # get the set
                   matches=list(set(temp) & set(temp2))
                   MLTLBins[iMLT][iL]+=list(chargingYa[matches])

                
       chargingArr+=list(Phi[keyTimes])
       chargingTimes+=list(rng[keyTimes])
    except:
           print iDate
    dt=dt+datetime.timedelta(days=1)
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/')

import pickle
pickle.dump(chargingArr, open('nonEclipseArr.p', 'wb'))
pickle.dump(chargingTimes, open('nonEclipseTimes.p', 'wb'))
for iMlt in range(nMLT):
    for iL in range(nL):
            lenMLTL[iMLT][iL]+=len(MLTLBins[iMLT][iL])
pickle.dump(MLTLBins, open('nonEclipseSorted.p', 'wb'))

pickle.dump(lenMLTL, open('nonEclipseSortedLen.p', 'wb'))


                   
        # calculate for each location if it is in eclipse or not
#

# have to get ephemeris and EFW data, resample EFW to one minute intervals
# use Michelle's formula to determine if in shadow or not
#

# need EFW data
