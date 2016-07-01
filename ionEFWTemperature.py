# ionEFWTemperature.py
#
# compare ion line with electron temperature long term
# indices for chargingStats are 23,38, 53, 68
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
import bisect
from itertools import islice
import pickle
import scipy
import scipy.stats
import pandas as pd
#os.chdir('/Users/loisks/Desktop/Functions/')
#import pickling
#os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/')


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
dateEnd='20130601'
#dateEnd='20150401'
#dateEnd='20130601'
dfs={}
endDT=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(dateStart, '%Y%m%d')
bins1=np.logspace(5,9,40)
#
# 

#print('here')
#IonLineC=[]
#IonChargingTimesC=[]
#ElectronTemperatureC=[]
#os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/PlasmaDensity200eV')
#lowIonsData=[]
#lowIonsTime=[]
#datesArr=['2013_2','2013_3','2013_4','2013_5']
#for iDate2 in range(len(datesArr)):
#   lowIonsData+=list(pickle.load(open("P_sat=rbspA_"+datesArr[iDate2], "rb")))
#   lowIonsTime+=list(pickle.load(open("time_P_sat=rbspA_"+datesArr[iDate2], "rb")))
#pdTime=pd.to_datetime(np.array(lowIonsTime))
#lowITime=[ datetime.datetime.strptime(str(i), '%Y-%m-%d %H:%M:%S') for i in pdTime]
#for isat in range(len(sat)):
#    dataEFW={}
#    dataHOPE={}
#
#    IonChargingTimes=[]
#    IonLineData=[]
#    ElectronTemperature=[]
#    # open the files
##   
##   # now load in HOPE data for this
#    while DT != endDT:
#     try:
#
#        date=datetime.datetime.strftime(DT, '%Y%m%d')
#        DT=datetime.datetime.strptime(date, '%Y%m%d')
#        # first get EFW
#        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+sat[isat])
#        f=glob.glob('*'+date+'*')
#        pyf=pycdf.CDF(f[0])
#        EFWdensity=(list(pyf['density'][...]))
#        EFWpotential=np.array(pyf['Vavg'][...])*-1.0
#        epochEFW=pd.DatetimeIndex(pyf['epoch'][...])
#        # now get HOPE data
#        os.chdir('..')
#        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+sat[isat])
#        f=glob.glob('*'+date+'*')
#        pyf=pycdf.CDF(f[0])
#        counts=pyf['Counts_P_Omni'][...]
#        HIE=np.average(pyf['HOPE_ENERGY_Ion'][...], axis=0)
#        IonFlux=pyf['FPDO'][...]
#        epochI=pyf['Epoch_Ion'][...]
#        epoch=epochI
#        Lion=pyf['L_Ion'][...]
#
#
#        # NEED TO METHODICALLY GO THROUGH THIS CODE AND DO GOOD JOB
#
#        os.chdir('..')
#        mTime=dt.date2num(epochI)
#        #
#        
#        # now get the background data
#        os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/BackgroundCounts')
#        bgCounts=pickle.load(open('bgCounts.p','rb'))
#        for iTime in range(len(epochI)):
#            for iEN in range(72):
#               idx = np.argmin(np.abs(LbinsArr - np.array(Lion)[iTime]))
#               try:
#                   counts[iTime][iEN]=counts[iTime][iEN]-bgCounts[iEN][idx]
#               except(ValueError):
#                   counts[iTime][iEN]=counts[iTime][iEN]
#        #
#        # extract an ion line
#        os.chdir('..')
#
#        # Now get Electron temperature
#
#        # # everything is loaded in 
#        # now put into pandas and do 1 minute averages
#        
#        dfpotential=pd.DataFrame(EFWpotential, index=epochEFW, columns=['potential'])
#      #  dfcFlags=pd.DataFrame(chargingFlags,index=epochEFW,columns=['charging'])
#        dfLion=pd.DataFrame(Lion, index=epochI, columns=['Lion'])
#
#        dfIonCounts={}
#        for iEn in range(72):
#            #temp=ElectronFlux[iEn]
#            tempI=np.swapaxes(counts,1,0)[iEn]
#            dfIonCounts[iEn]=pd.DataFrame(tempI,index=epochI,columns=['counts'])
#        rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
#        rng=rt[::1]
#        #
#        # Electron Temperatures
#        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_MOM_'+sat[isat])
#        f=glob.glob('*'+date+'*')
#        try: 
#            pyf=pycdf.CDF(f[0])
#        except:
#            continue
#
#
#        tpt=pyf['Tpar_e_200'][...]
#        tpp=pyf['Tperp_e_200'][...]
#        tpt[tpt>1e24]=np.nan
#        tpp[tpp>1e24]=np.nan
#        Tboth=np.array(tpt)*(2/3.)+np.array(tpp)*(1/3.)
#        htime=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
#        dfT=pd.DataFrame(Tboth, index=htime, columns=['T'])
#        T=dfT['T'].resample('1min', how='median').reindex(index=rng)
#        #
#        # now get T avg
#        try:
#            ionBiRight=bisect.bisect_right(lowITime, DT)
#            ionBiLeft=bisect.bisect_left(lowITime,DT+datetime.timedelta(days=1))
#            ionDensitydf=pd.DataFrame(lowIonsData[ionBiRight:ionBiLeft], index=lowITime[ionBiRight:ionBiLeft], columns=['density'])
#            ionDensity=ionDensitydf['density'].resample('1min', how='median').reindex(index=rng)
#        except(ValueError):
#            print("No Hope Data for " + str(DT))
#            continue 
#        # interpolate the data frames
#        EFWpot=dfpotential['potential'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
#        chargingTimes=np.where((np.array(EFWpot) < 0))[0]
#        LionC=np.array(dfLion['Lion'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan))[chargingTimes]
#        #cFlags=dfcFlags['charging'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
#        hd=pyf['Dens_e_200'][...]
#        hd[hd<0]=np.nan
#        hopeData={}
#        hopeData['epoch']=htime
#        hopeData['density']=hd*1.0e6
#       
#        dfdensity=pd.DataFrame(hopeData['density'], index=hopeData['epoch'], columns=['density'])
#        density=dfdensity['density'].resample('1min', how='median').reindex(index=rng)
#        Tavg=((np.array(ionDensity)*200) + (np.array(density)*np.array(T)))/( np.array(ionDensity)+np.array(density))
#    #EFLUX=[ [] for iNN in range(72)]
#        IonCounts=[ [] for iNN in range(EnBins)]
#
#        
#        for iEn in range(72):
#        #    EFLUX[iEn]=list(np.array(dfElectronFlux[iEn]['flux'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan))[chargingTimes])
#            IonCounts[iEn]+=list(np.array(dfIonCounts[iEn]['counts'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan))[chargingTimes])
#            #ElectronFluxData[iEn]+=list(EFLUX[iEn])
#        
#            
#        # now get the ion charging line
#        
#        chargingLine=[]
#        chargingTimesIon=[]
#        temperature=[]
#        #ElectronFluxT=[ [] for iNN in range(72)]
#        counts=np.swapaxes(IonCounts,1,0)
#
#        
#        for iTime in range(len(chargingTimes)):
#          if LionC[iTime] > 5: # must be near apogee
#            check=0
#            temp=np.where((counts[iTime]>0))[0][0:5] # should be clear
#            for iTemp in range(len(temp)):
#                # is it above like 6 eV
#
#                if ((temp[iTemp]>15) & (temp[iTemp]<50)):
#                    # is it above 50 counts at hte ion line
#                    try:
#                       if (counts[iTime][temp[iTemp]]-counts[iTime][temp[iTemp+2]])>10:
#                      # is it only a few thick under these same requirements
#                     #  if temp[iTemp+3] > temp[iTemp]+3:
#                        # this is a charging line
#                         chargingLine.append(HIE[temp[iTemp]])
#                         chargingTimesIon.append(rng[chargingTimes][iTime])
#                         temperature.append(Tavg[chargingTimes][iTime])
#                         check=1
#                    except(IndexError):
#                        # I have it deliberately error after 3
#                        break
#            if check!=1:
#              # nothing
#              continue
#
#
## make sure we get the right values here for ion extracted data 
#          
#        IonLineData+=list(chargingLine)
#        IonChargingTimes+=list(chargingTimesIon)
#        
#        ElectronTemperature+=list(temperature)
#        DT=DT+datetime.timedelta(days=1)
#
#     except:
#         print(date)
#         print(len(ElectronTemperature))
#         print(len(IonLineData))
#         DT=DT+datetime.timedelta(days=1)
#         continue
##
#    # combine both satellites 
#    ElectronTemperatureC+=ElectronTemperature
#    IonLineC+=IonLineData
#    IonChargingTimesC+=IonChargingTimes
#
#
#import pickle
#
#        
#os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
#IonLineC=np.array(IonLineC)
#tA=np.where((IonLineC < 10) | (IonLineC > 1000))
#IonLineC[tA]=np.nan
#pickle.dump(IonLineC, open('IonAlgoLineTavg.p','wb'))
#pickle.dump(IonChargingTimesC,open('IonAlgoTimeTavg.p','wb'))
#pickle.dump(ElectronTemperatureC, open('ElectronTavgC.p','wb'))


#
#
import pickle
os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/')
IonLineC=pickle.load(open('IonAlgoLineTAvg.p','rb'))
IonLineC=np.array(IonLineC)
tA=np.where((IonLineC < 10) | (IonLineC > 1000))
IonLineC[tA]=np.nan
#EFluxC=np.array(pickle.load(open('IonAlgoEflux.p','rb')))/1000.0
ETemp=np.array(pickle.load(open('ElectronTavgC.p', 'rb')))

#os.chdir('ChargingElectronFluxes')
EFW_T=pickle.load(open('EFW_T.p', 'rb'))
EFWSC=pickle.load(open('EFWSCT.p', 'rb'))
#os.chdir('..')
#Numbers=pickle.load(open('binnedelectronfluxes.p', 'rb'))
#dataNumbers=[]
#for iNum in range(len(Numbers)):
#      dataNumbers+=[bins1[iNum]]*Numbers[iNum]
#
#
# stats on ion
#seriesB = np.transpose(np.log10(ETemp))
#seriesA = np.transpose(np.log10(IonLineC))
#ts_res = scipy.stats.theilslopes(seriesA, seriesB, 0.95)
#tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
#spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
#print("Kendall's tau for  is {0} (p={1})".format(tau,p_value))
#print("Spearman's R for  is {0} (p={1})".format(spear,spear_p))
  

#
# rebin by logspace

left, width = 0.15, 0.8
bottom, height = 0.15, 0.8
left_h = left+width+0.03
bottom_h = bottom+height+0.03

rect_scatter = [left, bottom, width, height]
#rect_histx = [left, bottom_h, width, 0.18]
#rect_histy = [left_h+0.015, bottom, 0.18, height]
fig=plt.figure(1, figsize=(8,8))
axScatter = plt.axes(rect_scatter)
#axHistx = plt.axes(rect_histx)
#axHisty = plt.axes(rect_histy)

axScatter.plot(np.array(EFW_T),-1*np.array(EFWSC), '.', color='darkorange')
axScatter.plot(np.array(ETemp),IonLineC, '.', color='red')
axScatter.set_xlim(1e3,1e4)
axScatter.set_ylim(1e0,1e3)
axScatter.set_yscale('log')
axScatter.set_xscale('log')
plt.draw()
font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
axScatter.set_ylabel('Spacecraft Potential [-V]', fontweight='bold', fontsize=22)
axScatter.set_xlabel('<T$_{e}$> [eV]',fontweight='bold', fontsize=22)
#xx=np.arange(np.nanmin(seriesB),np.nanmax(seriesB))
#axScatter.plot(10**np.array(xx), 10**np.array(ts_res[1] + ts_res[0] * xx),'-', linewidth=2, color='k')
#axHistx.hist(np.array(dataNumbers[iEn])/1000.0, bins=bins1, stacked=True, color='blue')
#axHistx.hist(np.array(EFW_T), bins=bins1,stacked=True, color='blue')
#axHistx.set_xlim( axScatter.get_xlim() )
#axHistx.set_yscale('log')
#axHistx.set_xscale('log')
axScatter.tick_params(axis='both', which='major', labelsize=22)
#axHistx.tick_params(axis='x', bottom='off', top='off', labelbottom='off', labeltop='off')    
    #
# plot EFW stats
#fig2=plt.figure()
#ax2=fig2.add_subplot(111)
#ax2.plot(np.array(EFluxesE[iEn])/1000.0,EFWSC, '.', color='lightseagreen')
#ax2.set_xlim(1e6,1e9)
#ax2.set_ylim(1e0,1e3)
#ax2.set_yscale('log')
#ax2.set_xscale('log')
#plt.draw()
#font = {'family' : 'normal',
#                  'weight' : 'bold',
#                  'size'   : 22}
#plt.rc('font', **font)
#plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
#ax2.set_ylabel('Spacecraft Potential [V]', fontweight='bold')
#ax2.set_xlabel('Energy Flux',fontweight='bold')
subdir_name='NewIonStats'
if not os.path.exists(subdir_name):
           os.umask(0) # unmask if necessary
           os.makedirs(subdir_name) 
os.chdir(subdir_name)
plt.savefig('lowBoth_EFWT_Energy.png')
plt.close()
os.chdir('..')
