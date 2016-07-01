# EFWscatterP.py
#
# check the positive potential points real quick
#
# LKS August 2015
#
#
# imports
import h5py
import numpy as np
from spacepy import pycdf
import os
import glob
import datetime
import matplotlib.dates as dt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#
# read in correct file
date='20130316'
hourGap=5
os.chdir('EFW_L3_A')
f=glob.glob('*'+date+'*')
pyf=pycdf.CDF(f[0])
density=(list(pyf['density'][...]))
potential=pyf['Vavg'][...]
epochE=dt.date2num(pyf['epoch'][...])
MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
os.chdir('..')
#
#
fig=plt.figure()
ax1=fig.add_subplot(111)
plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
ax1.set_ylabel('Potential', fontsize=25, fontweight='bold')
ax1.set_xlabel('Time', fontsize=25, fontweight='bold')
font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 22}
plt.rc('font', **font)
plt.scatter(epochE, -1*potential, s=40, c='b')
hfmt = dt.DateFormatter('%H:%M')
ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
ax1.xaxis.set_major_formatter(hfmt)
subdir_name='ChargingScatter'
if not os.path.exists(subdir_name):
     os.umask(0) # unmask if necessary
     os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig(date+'_chargingScatter.pdf')
plt.close(fig)
os.chdir('..')
