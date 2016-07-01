# magLatTest.py
#
# Test and see if magnetic latitude is a factor in s/c charging
# use the rainbow plots from rad belt work
#
# LKS, August 2015
#
#

# get the data
import h5py
import numpy as np
from spacepy import pycdf
import os
import glob
import datetime
import matplotlib.dates as dt
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#
#
hourGap=5
sat=['A','B']
isat=0
uniqueDates=['20130126',
 '20130217',
 '20130309',
 '20130310',
 '20130327',
 '20131127',
 '20141223',
 '20141226',
 '20150106',
 '20150121',
 '20150126',
 '20150130',
 '20150131',
 '20150201',
 '20150202',
 '20150204',
 '20150205',
 '20150207',
 '20150208',
 '20150210',
 '20150217',
 '20150222',
 '20150223',
 '20150224',
 '20150225',
 '20150227',
 '20150228',
 '20150301',
 '20150302',
 '20150303',
 '20150305',
 '20150307',
 '20150308',
 '20150309',
 '20150310',
 '20150311',
 '20150312',
 '20150313',
 '20150314',
 '20150315',
 '20150317',
 '20150318',
 '20150319',
 '20150320',
 '20150321',
 '20150322',
 '20150323',
 '20150324',
 '20150325',
 '20150326',
 '20150327',
 '20150328',
 '20150329',
 '20150330',
 '20150331',
 '20150401',
 '20150402',
 '20150403',
 '20150404',
 '20150405',
 '20150411',
 '20150412',
 '20150415',
 '20150417',
 '20150505',
 '20150510']
for iDate in range(len(uniqueDates)):
     try:
       date=uniqueDates[iDate]
       print date
       # first get EFW
       os.chdir('EFW_L3_A')
       f=glob.glob('*'+uniqueDates[iDate]+'*')
       pyf=pycdf.CDF(f[0])
       density=(list(pyf['density'][...]))
       potential=pyf['Vavg'][...]
       epochE=dt.date2num(pyf['epoch'][...])
       MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
       LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
       MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
       #
       # now get HOPE data
       os.chdir('..')
       os.chdir('HOPE_L3_A')
       f=glob.glob('*'+uniqueDates[iDate]+'*')
       pyf=pycdf.CDF(f[0])
       
       electronFlux=pyf['FEDO'][...]
       electronEnergy=pyf['HOPE_ENERGY_Ele'][...]
       electronFlux=electronFlux*electronEnergy
       epoch=pyf['Epoch_Ele'][...]
       L=pyf['L_Ele'][...]
       MLT=pyf['MLT_Ele'][...]
       #MLAT=pyf['MLAT_Ele'][...]
       os.chdir('..')
       mTime=dt.date2num(epoch)
       #
       # now we have both
       fig=plt.figure()
       cmap=plt.cm.jet
       cbmin=8
       cbmax=11
       ax1=fig.add_subplot(111)
       plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.32)
       ax1.set_ylabel('Energy (eV)', fontsize=25, fontweight='bold')
       ax1.set_xlabel('Time', fontsize=25, fontweight='bold')
       ax1.set_yscale('log')
       mEnergies=np.mean(electronEnergy, axis=0)
       X,Y=np.meshgrid(mTime, mEnergies)
       mdata=ma.masked_invalid(np.log10(electronFlux).transpose())
       col=ax1.pcolormesh(X,Y,mdata,cmap='jet', vmax=cbmax, vmin=cbmin)
       font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}
       plt.rc('font', **font)
       cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
       cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(cbmin,cbmax))
       cb.set_label('Energy Flux', fontsize=25, fontweight='bold')
       ax1.tick_params(axis='both', which='major', labelsize=22)
       cb.ax.tick_params(labelsize=30)
       hfmt = dt.DateFormatter('%H:%M')
       ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
       ax1.xaxis.set_major_formatter(hfmt)
       ax1.set_ylim(1e2, 50000)

       
       # also overplot S/C Potential
       ax2=ax1.twinx()
       #ax2.plot(epochE,-1*potential , lw=1, c='Black')
       #ax2.plot(epochE,-1*potential, lw=2, c=MLAT, cmap=cm.hot)
       points=np.array([np.array(epochE),np.array(-1*potential)]).T.reshape(-1,1,2)
       segments=np.concatenate([points[:-1], points[1:]], axis=1)
       lc=LineCollection(segments, cmap=plt.get_cmap('gray'),norm=plt.Normalize(-20,10))
       lc2=LineCollection(segments, cmap=plt.get_cmap('gray'))
       lc2.set_array(np.array(np.zeros(len(MLAT))))
       lc2.set_linewidth(6)
       ax2.add_collection(lc2)
       ax2.set_ylim(np.min(-1*potential), np.max(-1*potential))
       lc.set_array(np.array(MLAT))
       lc.set_linewidth(3)
       ax2.add_collection(lc)
       ax2.set_ylabel('S/C Potential', fontsize=25, fontweight='bold')
   
       # Add L LABEL
       ax3 = fig.add_axes(ax1.get_position())
       ax3.patch.set_visible(False)
       ax3.yaxis.set_visible(False)
       ax3.set_yscale('log')
       for spinename, spine in ax3.spines.iteritems():
           if spinename != 'bottom':
               spine.set_visible(False)
       ax3.spines['bottom'].set_position(('outward', 65))
       L_conc=[round(x, 2) for x in L]
       p3=ax3.plot(mTime,electronFlux  , 'g', alpha=0)
       # I hate doing this but
       spacing = len(L_conc)/7
       spacedArray=L_conc[::spacing]
       ax3.set_xticks(epoch[::spacing])
       nans=np.where(np.array(L_conc[::spacing]) < 0)
       if len(nans[0])>0:
           try:
               newSpot = L_conc[spacing*7+5]
               print newSpot
               if newSpot < 0:
                   newSpot=L_conc[spacing*7-5]
                   print newSpot
               spacedArray[nans[0][0]]=newSpot
           except(IndexError):
               newSpot=L_conc[spacing*7-5]
               print newSpot
               spacedArray[nans[0][0]]=newSpot
       ax3.set_xticklabels(spacedArray)
       ax3.set_xlabel('L-Shell', fontsize=22, fontweight='bold')
       #
       # Add MLT Label
       ax5 = fig.add_axes(ax1.get_position())
       ax5.patch.set_visible(False)
       ax5.yaxis.set_visible(False)
       ax5.set_yscale('log')
       for spinename, spine in ax3.spines.iteritems():
           if spinename != 'bottom':
               spine.set_visible(False)
       ax5.spines['bottom'].set_position(('outward', 125))
       MLT_conc=[round(x, 1) for x in MLT]
       p4=ax5.plot(mTime,  electronFlux, 'purple', alpha=0)
       spacing = len(MLT_conc)/8
       ax5.set_xticks(mTime[::spacing])
       ax5.set_xticklabels(MLT_conc[::spacing])
       ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
       

       subdir_name='MagLatChargingPlots'
       if not os.path.exists(subdir_name):
            os.umask(0) # unmask if necessary
            os.makedirs(subdir_name, 0777) 
       os.chdir(subdir_name)#
       fig.set_size_inches(13,9)
       plt.savefig(date+'_'+sat[isat]+'_spectrogram.pdf')
       plt.close(fig)
       os.chdir('..')
     except(IOError):
        print 'bad date: ' + date  


