# caseStudy.py
#
# see if there is a link between certain energy electron fluxes
# and spacecraft charging
# fluence over an hour/time scale with multiple correlation variables
# try scikit elastic net approach
# look at protons too
# look at temperatures too

#
#
#
# LKS, August 2015
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
from scipy.interpolate import interp1d
import scipy
import scipy.stats

#
# threshold
threshold=0
maxVal=-200
scatter='on'
hitTest='off'

# do a -50 V + and - test 
#
# unique dates for greater than 50V
uniqueDates=['20130208',
 '20130213',
 '20130214',
 '20130220',
 '20130301',
 '20130309',
 '20130310',
 '20130311',
 '20130313',
 '20130315',
 '20130316',
 '20130317',
 '20130320',
 '20130321',
 '20130322',
 '20130323',
 '20130324',
 '20130328',
 '20130329',
 '20130331',
 '20130401',
 '20130402',
 '20130404',
 '20130406',
 '20130407',
 '20130424',
 '20130501',
 '20130518',
 '20130525',
 '20130629',
 '20130714',
 '20140903',
 '20140912',
 '20141010',
 '20141011',
 '20141014',
 '20141015',
 '20141018',
 '20141020',
 '20141021',
 '20141022',
 '20141027',
 '20141028',
 '20141104',
 '20141105',
 '20141107',
 '20141110',
 '20141112',
 '20141115',
 '20141116',
 '20141117',
 '20141118',
 '20141119',
 '20150107',
 '20150302',
 '20150307',
 '20150308',
 '20150309',
 '20150311',
 '20150312',
 '20150313',
 '20150314',
 '20150315',
 '20150316',
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
 '20150402',
 '20150403',
 '20150404',
 '20150410',
 '20150411',
 '20150414',
 '20150415',
 '20150416',
 '20150417']

#
# open up the files
# interpolate HOPE onto EFW s/c
# compare by energy

interpFlux=[[] for i in range(72)]
allPotential=[]
#
# now load in HOPE data for this

for iDate in range(len(uniqueDates)):
       date=uniqueDates[iDate]
       # first get EFW
       os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_A')
       f=glob.glob('*'+uniqueDates[iDate]+'*')
       pyf=pycdf.CDF(f[0])
       density=(list(pyf['density'][...]))
       potential=-1*pyf['Vavg'][...]
       epochE=dt.date2num(pyf['epoch'][...])
       MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
       LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
       MLATE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
       #
       # now get HOPE data
       os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_A')
       f=glob.glob('*'+uniqueDates[iDate]+'*')
       pyf=pycdf.CDF(f[0])
       
       electronFlux=pyf['FEDO'][...]
       electronEnergy=pyf['HOPE_ENERGY_Ele'][...]
       electronFlux=electronFlux*electronEnergy
       epoch=pyf['Epoch_Ele'][...]
       L=pyf['L_Ele'][...]
       MLT=pyf['MLT_Ele'][...]
       os.chdir('..')
       mTime=dt.date2num(epoch)

       # now interpolate onto EFW times

       f=interp1d(epochE,potential)
       try:
           temp=f(mTime)
       except(ValueError):
           mTime[mTime<epochE[0]]=epochE[0]
           mTime[mTime>epochE[-1]]=epochE[-1]
           temp=f(mTime)
       null=np.where(temp>threshold)[0]
       temp[null]=np.nan
       nullMax=np.where(temp<maxVal)[0]
       temp[nullMax]=np.nan
       allPotential=list(temp)
       allPotential=np.array(allPotential)[~np.isnan(allPotential)]
       os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
       for iEn in range(72):
           temp=np.swapaxes(electronFlux,1,0)[iEn]
           temp[null]=np.nan
           temp[nullMax]=np.nan
           interpFlux[iEn]=list(temp)
       # determine correlation coefficient for each day    
           interpFlux[iEn]=np.array(interpFlux[iEn])[~np.isnan(interpFlux[iEn])]
           seriesB = np.transpose(np.log10(interpFlux[iEn]))
           seriesA = np.transpose(allPotential)
           ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
           tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
           spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
           print("Kendall's tau for (Energy={0}) is {1} (p={2})".format(iEn,tau,p_value))
           print("Spearman's R for (Energy={0}) is {1} (p={2})".format(iEn,spear,spear_p))
           fig2=plt.figure()
           ax2=fig2.add_subplot(111)
           ax2.plot(seriesA, seriesB, '.', color='k')
           xx=np.arange(-200,0)
           ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
           plt.draw()
           font = {'family' : 'normal',
                      'weight' : 'bold',
                      'size'   : 22}
           plt.rc('font', **font)
           plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
           ax2.set_ylabel('Log(Energy Flux)', fontweight='bold')
           ax2.set_xlabel('$\phi$',fontweight='bold')
           ax2.set_ylim(9,12)
           ax2.set_xlim(maxVal,threshold)
           ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(-50, 12), xycoords='data', color='k')
           subdir_name='ScatterCaseStudies'
           if not os.path.exists(subdir_name):
               os.umask(0) # unmask if necessary
               os.makedirs(subdir_name, 0777) 
           os.chdir(subdir_name)
           plt.savefig('date='+str(uniqueDates[iDate])+"_energy="+str(iEn)+'_threshold='+str(threshold)+'.pdf')
           plt.close()
           os.chdir('..')



