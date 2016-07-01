# LTbins.py

# bin by MLT and L-Shell for charging
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
import matplotlib.dates as dates
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
#
#
# establish bin sizes
binWidth=10
maxCharging=-200
nX=25
chargingBins=np.linspace(maxCharging,50,nX)
nMLT=48
nL=25
nMLAT=41
MLTbins=np.linspace(0.25, 23.75, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
MLATarr=[ [] for i in range(nMLAT)]
Lbins=np.linspace(0.25, 6.25, nL)
MLATbins=np.linspace(-20,20, nMLAT)

#
# unique dates for less than 50V
# important to normalize this


# now load in HOPE data for this
interpFlux=[[] for i in range(72)]
allPotential=[]
allMLT=[]
allL=[]
allMLAT=[]
TotalFluence=[[] for i in range(72)]
allPotential=[]
iDate='20130101'
dt=datetime.datetime.strptime(iDate,'%Y%m%d')
endDate='20150501'
while iDate != endDate:
#for iDate in range(len(uniqueDates)):
   try:
       iDate=datetime.datetime.strftime(dt,'%Y%m%d')
       dt=datetime.datetime.strptime(iDate, '%Y%m%d')

       # first get EFW
       os.chdir('EFW_L3_A')
       f=glob.glob('*'+iDate+'*')
       pyf=pycdf.CDF(f[0])
       density=(list(pyf['density'][...]))
       potential=-1*pyf['Vavg'][...]
       epochE=dates.date2num(pyf['epoch'][...])
       MLTE=np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]
       LE=np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]
       MLATE=np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]
       # identify points of charging
       #
       # now bin everything
       for i in range(nMLT):
           try:
               MLTindex=np.where((MLTE >= i*0.5) & (MLTE < (i+1)*0.5))[0]
               MLTarr[i].extend(list(potential[MLTindex]))
           except:
               MLTarr[i].extend([])
       for i in range(nL):
           try:
               Lindex=np.where((LE >= i*0.25) & (LE < (i+1)*0.25))[0]
               Larr[i].extend(list(potential[Lindex]))
           except:
               Larr[i].extend([])
       for i in range(nMLAT):
           try:
               MLATindex=np.where((MLATE >= i*1-10) & (MLATE < (i+1)-10))[0]
               MLATarr[i].extend(list(potential[MLATindex]))
           except:
               MLATarr[i].extend([])
       os.chdir('..')
       dt=dt+datetime.timedelta(days=1)
   except(IndexError):
       os.chdir('..')
       dt=dt+datetime.timedelta(days=1)           
       #
       # Now sort by threshold
# get the totals
tMLTarr=[ [] for i in range(nMLT)]
tLarr=[ [] for i in range(nL)]
tMLATarr=[ [] for i in range(nMLAT)]
for iMLT in range(nMLT):
    try:
        tMLTarr[iMLT]=len(np.array(MLTarr[iMLT]))*1.0
    except:
        tMLTarr[iMLT]=0
for iL in range(nL):
    try:
        tLarr[iL]=len(np.array(Larr[iL]))*1.0
    except:
        tLarr[iL]=0
for iMLAT in range(nMLAT):
       try:
                tMLATarr[iMLAT]=len(np.array(MLATarr[iMLAT]))*1.0
       except:
                tMLATarr[iMLAT]=0        

prevMLTarr=[ 0 for i in range(nMLT)]
prevLarr=[ 0 for i in range(nL)]
prevMLATarr=[ 0 for i in range(nMLAT)]
for iBar in range(nX):
        sMLTarr=[ [] for i in range(nMLT)]
        sLarr=[ [] for i in range(nL)]
        sMLATarr=[ [] for i in range(nMLAT)]
        if maxCharging+iBar*binWidth >= 0:
            prevMLTarr=[ 0 for i in range(nMLT)]
            prevLarr=[ 0 for i in range(nL)]
            prevMLATarr=[ 0 for i in range(nMLAT)]           
        for iMLT in range(nMLT):
            indexes=np.where((np.array(MLTarr[iMLT])>=maxCharging+iBar*binWidth) & (np.array(MLTarr[iMLT]) < maxCharging+(1+iBar)*binWidth))[0]
            try:
                sMLTarr[iMLT]=len(np.array(MLTarr[iMLT])[indexes])/tMLTarr[iMLT]
                prevMLTarr[iMLT]+=sMLTarr[iMLT]
            except:
                sMLTarr[iMLT]=0
        for iL in range(nL):
            indexes=np.where((np.array(Larr[iL])>=maxCharging+iBar*binWidth) & (np.array(Larr[iL]) < maxCharging+(1+iBar)*binWidth))[0]
            try:
                sLarr[iL]=len(np.array(Larr[iL])[indexes])/tLarr[iL]
                prevLarr[iL]+=sLarr[iL]
            except:
                sLarr[iL]=0
        for iMLAT in range(nMLAT):
            indexes=np.where((np.array(MLATarr[iMLAT])>=maxCharging+iBar*binWidth) & (np.array(MLATarr[iMLAT]) < maxCharging+(1+iBar)*binWidth))[0]
            try:
                sMLATarr[iMLAT]=len(np.array(MLATarr[iMLAT])[indexes])/tMLATarr[iMLAT]
                prevMLATarr[iMLAT]+=sMLATarr[iMLAT]
            except:
                sMLATarr[iMLAT]=0

                
        # plot bar for each bin size
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.set_xlim(0,24)
        ax2.bar(MLTbins, prevMLTarr, 0.5, color='b')
        ax2.set_title("Normalized Charging less than " + str(maxCharging+(1+iBar)*binWidth))
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Fraction of Total', fontweight='bold')
        ax2.set_xlabel('MLT ',fontweight='bold')
        os.chdir('MLT_barplots')              
        plt.savefig("Norm_MLT_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth))
        os.chdir('..')

        # plot bar for L
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.bar(Lbins, prevLarr, 0.5, color='b')
        ax2.set_title("Normalized Charging less than " + str(maxCharging+(1+iBar)*binWidth))
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Fraction of Total', fontweight='bold')
        ax2.set_xlabel('L-Shell ',fontweight='bold')
        os.chdir('L_barplots')              
        plt.savefig("Norm_L_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth))
        os.chdir('..')

        # plot bar for MLAT
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.bar(MLATbins, prevMLATarr, 1, color='b')
        ax2.set_title("Normalized Charging less than " + str(maxCharging+(1+iBar)*binWidth))
        #ax2.set_title("Normalized Charging From " + str(maxCharging+iBar*binWidth) + ' to ' + str(maxCharging+(1+iBar)*binWidth))
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Fraction of Total', fontweight='bold')
        ax2.set_xlabel('MLAT ',fontweight='bold')
        os.chdir('MLAT_barplots')              
        plt.savefig("Norm_MLAT_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth))
        os.chdir('..')

  
