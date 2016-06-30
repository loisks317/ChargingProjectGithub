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
MLTbins=np.linspace(0, 24, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
Lbins=np.linspace(1.5, 7.25, nL)
MLTLBins=[[ [] for i in range(nL)] for j in range(nMLT)]
lenMLTL=[[ 0 for i in range(nL)] for j in range(nMLT)]
#

import pickle
#chargingArr=pickle.load(open('nonEclipseArr.p', 'rb'))
#chargingTimes=pickle.load( open('nonEclipseTimes.p', 'rb'))
MLTLbins=pickle.load( open('nonEclipseSorted.p', 'rb'))
#lenMLTL=pickle.load(open('nonEclipseSortedLen.p', 'rb'))
for iMLT in range(nMLT):
    for iL in range(nL):
        lenMLTL[iMLT][iL] = len(MLTLbins[iMLT][iL])

# # plot?
# intitialize the figure
fig=plt.figure()
ax=fig.add_subplot(111, polar=True)
datah_m=ma.masked_invalid(np.array(lenMLTL).transpose())
mlt=np.array(MLTbins)*(15*np.pi/180.)
X,Y=np.meshgrid(mlt, Lbins)
ax.set_ylim(0, 7.25)
dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
vmin=1
vmax=3
col=ax.pcolormesh( X, Y, np.log10(datah_m), cmap='jet', vmin=vmin, vmax=vmax )
#col=ax.pcolormesh( X, Y, datah_m, cmap='jet', vmin=vmin, vmax=vmax )
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
#cb=plt.colorbar(col,cax=cbaxes)
plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
xT2=['', '2', '', '4', '', '6', '', '8']
xL2=['00','', '06','', '12','', '18']
cb.set_label("Occurrences", fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=25)    
cb.ax.tick_params(labelsize=35) 
ax.set_yticklabels(xT2, fontsize=30)
ax.set_xticklabels(xL2, fontsize=30)
ax.grid(True)
plt.draw()
subdir_name='MLTL_Dialplots'
if not os.path.exists(subdir_name):
   os.umask(0) # unmask if necessary
   os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig('MLTL_Dialplots_NonEclipseNegCharging'+'.pdf')
plt.close(fig)
os.chdir('..')  
