# thermalCurrent.py
#
# calculate the thermal current and compare to threshold of 70 nA/m^2
# John Wygant says this formula is ne sqrt[Te]
#
# October 2015, LKS
#
#
# imports
import h5py
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
#
#
date='20130214'
kelvin=11604.505
hourGap=4
#
#

os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_MOM_A')
f=glob.glob('*'+date+'*')
pyf=pycdf.CDF(f[0])
hd=np.array(pyf['Dens_e_200'][...])
tpt=np.array(pyf['Tpar_e_200'][...])
tpp=np.array(pyf['Tperp_e_200'][...])
eT=dates.date2num(pyf['Epoch_Ele'][...])
#
# do conversion here to get micro amps over m^2
# this is n* sqrt(T) with n in cm^-3 and T in eV. There needs to be
# a 0.0642 conversion factor in here too
Te=(tpt+2*tpp)/3.0

tCur=hd*np.sqrt(Te)*0.094
#
# compare with potential from EFW
os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_A')
f=glob.glob('*'+date+'*')
pyf=pycdf.CDF(f[0])
Vavg=-1*np.array(pyf['Vavg'][...])
epoch=dates.date2num(pyf['epoch'][...])
#
# plot both of these
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/')
fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.set_ylabel('Current ($\mu$A/m$^{2}$)', fontsize=22, fontweight='bold')
ax1.set_xlabel('Time ', fontsize=22, fontweight='bold')
plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
plt.rc('font', **font)
ax1.set_ylim(0,30)
p1,=ax1.plot(eT, tCur, c='lightseagreen', lw=2, label="Thermal Current")
ax1.yaxis.label.set_color('lightseagreen')
#p2,=ax1.plot(eT, np.zeros(len(eT))+70, c='purple', lw=2, label='Photoemission')
ax2=ax1.twinx()
ax2.set_ylabel("$\phi$ [-V]", fontsize = 22, fontweight='bold')
ax2.yaxis.label.set_color('darkviolet')
ax2.plot(epoch, -1*np.array(Vavg), c='darkviolet', lw=2, label="S/C")
#plt.legend([p1,p2], bbox_to_anchor=[1.2, 0.7])
hfmt = dates.DateFormatter('%H:%M')
ax1.xaxis.set_major_locator(dates.HourLocator(interval=hourGap))
ax1.xaxis.set_major_formatter(hfmt)

# save the figure
subdir_name='ThermalCurrent'
if not os.path.exists(subdir_name):
     os.umask(0) # unmask if necessary
     os.makedirs(subdir_name) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig(date+'_'+'ThermalCurrent.pdf')
plt.close(fig)
os.chdir('..')
