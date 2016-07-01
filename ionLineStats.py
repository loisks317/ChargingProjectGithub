# ionLineStats.py
#
# compare ion line with electron fluxes long term
#
# LKS, September 2015
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
from collections import deque,Counter
from bisect import insort, bisect_left
from itertools import islice
import pickle
import pandas as pd


def runningMedian(seq, M):
    """
     Purpose: Find the median for the points in a sliding window (odd number in size) 
              as it is moved from left to right by one point at a time.
      Inputs:
            seq -- list containing items for which a running median (in a sliding window) 
                   is to be calculated
              M -- number of items in window (window size) -- must be an integer > 1
      Otputs:
         medians -- list of medians with size N - M + 1
       Note:
         1. The median of a finite list of numbers is the "center" value when this list
            is sorted in ascending order. 
         2. If M is an even number the two elements in the window that
            are close to the center are averaged to give the median (this
            is not by definition)
    """   
    seq = iter(seq)
    s = []   
    m = M // 2

    # Set up list s (to be sorted) and load deque with first window of seq
    s = [item for item in islice(seq,M)]    
    d = deque(s)

    # Simple lambda function to handle even/odd window sizes    
    median = lambda : s[m] if bool(M&1) else (s[m-1]+s[m])*0.5

    # Sort it in increasing order and extract the median ("center" of the sorted window)
    s.sort()    
    medians = [median()]   

    # Now slide the window by one point to the right for each new position (each pass through 
    # the loop). Stop when the item in the right end of the deque contains the last item in seq
    for item in seq:
        old = d.popleft()          # pop oldest from left
        d.append(item)             # push newest in from right
        del s[bisect_left(s, old)] # locate insertion point and then remove old 
        insort(s, item)            # insert newest such that new sort is not required        
        medians.append(median())  
    return medians

#
# data
sat=['A', 'B']


## setup
## extract the ion line for each day
## get the electron fluxes by energy channel
## use pandas to interpolate the data
## pickle all in folder 

EnBins=72
Lbins=56
LbinsArr=np.linspace(1, 6.5,Lbins)
dateStart='20130201'
dateEnd='20150401'
#dateEnd='20130601'
dfs={}
endDT=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(dateStart, '%Y%m%d')

#
IonLineC=[]
EFluxC=[ [] for i in range(72)]
for isat in range(len(sat)):
    dataEFW={}
    dataHOPE={}


    IonLineData=[]
    ElectronFluxData=[[] for i in range(72)]
    # open the files
#   
#   # now load in HOPE data for this
    while DT != endDT:
     try:

        date=datetime.datetime.strftime(DT, '%Y%m%d')
        DT=datetime.datetime.strptime(date, '%Y%m%d')
        # first get EFW
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+sat[isat])
        f=glob.glob('*'+date+'*')
        pyf=pycdf.CDF(f[0])
        EFWdensity=(list(pyf['density'][...]))
        EFWpotential=np.array(pyf['Vavg'][...])*-1.0
        epochEFW=pd.DatetimeIndex(pyf['epoch'][...])
        # now get HOPE data
        os.chdir('..')
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+sat[isat])
        f=glob.glob('*'+date+'*')
        pyf=pycdf.CDF(f[0])
        
        
        counts=pyf['Counts_P_Omni'][...]
        HIE=pyf['HOPE_ENERGY_Ion'][...]
        EE=pyf['HOPE_ENERGY_Ele'][...]
        IonFlux=pyf['FPDO'][...]
        ElectronFlux=np.swapaxes(pyf['FEDO'][...],1,0)*np.swapaxes(EE,1,0)
        
        epochI=pyf['Epoch_Ion'][...]
        epochE=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
        epoch=epochI
        Lion=pyf['L_Ion'][...]
        MLT=pyf['MLT_Ele'][...]
        Lele=pyf['L_Ele'][...]



        os.chdir('..')
        mTime=dt.date2num(epochI)
        #
        # now get the background data
        os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/BackgroundCounts')
        bgCounts=pickle.load(open('bgCounts.p','rb'))
        for iTime in range(len(epochI)):
            for iEN in range(72):
               idx = np.argmin(np.abs(LbinsArr - np.array(Lion)[iTime]))
               try:
                   counts[iTime][iEN]=counts[iTime][iEN]-bgCounts[iEN][idx]
               except(ValueError):
                   counts[iTime][iEN]=counts[iTime][iEN]
        #
        # extract an ion line
        os.chdir('..')
        chargingLine=[]
        chargingTimes=[]

        gradient=np.zeros((len(epoch),55))
        for iTime in range(len(epoch)):
            cCount=-1
            if Lion[iTime]>4:
                for iEn in range(55):
                 gradient[iTime][iEn]=((counts[iTime][iEn+1]-counts[iTime][iEn])/(counts[iTime][iEn]*1.0))
                 a=np.where(gradient[iTime] > 10)[0]
                 try:
                     test=len(a)
                     cCount=a[0]
                 except:
                     cCount=0
                 if cCount >= 1:
                     sc=HIE[iTime][cCount]
                 else:
                     sc=0
            chargingLine.append(sc)
            chargingTimes.append(epoch[iTime])

        # not sure if smoothing is working                  
        #cLine=runningMedian(chargingLine, 10)
        #chargingTimes=runningMedian(dt.date2num(chargingTimes),10)
        #chargingTimes=dt.num2date(chargingTimes)
        cLine=chargingLine
        chargingTimes=chargingTimes

        
        # set up the data frames
        dfpotential=pd.DataFrame(EFWpotential, index=epochEFW, columns=['potential'])
        dfcharging=pd.DataFrame(cLine,index=pd.DatetimeIndex(chargingTimes),columns=['charging'])
        dfElectronFlux={}
        for iEn in range(72):
            temp=ElectronFlux[iEn][Lele>4]
            dfElectronFlux[iEn]=pd.DataFrame(temp,index=epochE[Lele>4],columns=['flux'])
        rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
        rng=rt[::5]

        # interpolate the data frames
        EFWpot=dfpotential['potential'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
        IonLine=dfcharging['charging'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
        EFLUX={}
        for iEn in range(72):
            EFLUX[iEn]=dfElectronFlux[iEn]['flux'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)

        IonLineData+=list(IonLine)
        for iEn in range(72):
            ElectronFluxData[iEn]+=list(EFLUX[iEn])

        DT=DT+datetime.timedelta(days=1)

     except:
         print date
         DT=DT+datetime.timedelta(days=1)
         continue
#
    # combine both satellites 
    for iEn in range(72):
        EFluxC[iEn]+=ElectronFluxData[iEn]
    IonLineC+=IonLineData

# plot the results
for i in range(72):
           left, width = 0.1, 0.65
           bottom, height = 0.1, 0.65
           bottom_h = left_h = left+width+0.02

           rect_scatter = [left, bottom, width, height]
           rect_histx = [left, bottom_h, width, 0.2]
           rect_histy = [left_h, bottom, 0.2, height]
           fig=plt.figure(1, figsize=(8,8))
           axScatter = plt.axes(rect_scatter)

           # the scatter plot:
           axScatter.plot(EFluxC[i],IonLineC,'.', c='blue')
           
           axScatter.set_ylabel('S/C Potential [V]')
           axScatter.set_xlabel('Energy Flux')
           axScatter.set_yscale('log')
           axScatter.set_xscale('log')
           axScatter.set_xlim(1e7,1e11)
           font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}
           plt.rc('font', **font)
           #cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
           #cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(0,4))
           #cb.set_label('Number of Points', fontsize=25, fontweight='bold')
           axScatter.tick_params(axis='both', which='major', labelsize=22)
   
           subdir_name='IonStats'
           if not os.path.exists(subdir_name):
              os.umask(0) # unmask if necessary
              os.makedirs(subdir_name, 0777) 
           os.chdir(subdir_name)#
           fig.set_size_inches(13,9)
           print 'saving figure'
           plt.savefig('ionLine_Eflux='+str(i)+'.png')
           plt.close(fig)
           os.chdir('..')
            
               
