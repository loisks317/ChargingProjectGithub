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
from itertools import islice
import pickle
import pandas as pd
import scipy
import scipy.stats
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
#dateEnd='20130401'
dateEnd='20150401'
#dateEnd='20130601'
dfs={}
endDT=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(dateStart, '%Y%m%d')
bins1=np.logspace(5,9,40)
#
# 

print 'here'
#IonLineC=[]
#IonChargingTimesC=[]
#ElectronPressureC=[]
#for isat in range(len(sat)):
#    dataEFW={}
#    dataHOPE={}
#
#    IonChargingTimes=[]
#    IonLineData=[]
#    ElectronPressure=[]
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
#        chargingFlagsTemp=pyf['flags_charging_bias_eclipse'][...]
#        chargingFlags=np.zeros(len(epochEFW))
#        for iTime in range(len(epochEFW)):
#            if any(np.array(chargingFlagsTemp[iTime]) >0):
#                chargingFlags[iTime]=1
#            else:
#                chargingFlags[iTime]=0
#        # now get HOPE data
#        os.chdir('..')
#        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+sat[isat])
#        f=glob.glob('*'+date+'*')
#        pyf=pycdf.CDF(f[0])
#        
#        
#        counts=pyf['Counts_P_Omni'][...]
#        HIE=np.average(pyf['HOPE_ENERGY_Ion'][...], axis=0)
#        #EE=pyf['HOPE_ENERGY_Ele'][...]
#        IonFlux=pyf['FPDO'][...]
#        #ElectronFlux=np.swapaxes(pyf['FEDO'][...],1,0)*np.swapaxes(EE,1,0)
#        
#        epochI=pyf['Epoch_Ion'][...]
#        #epochE=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
#        epoch=epochI
#        Lion=pyf['L_Ion'][...]
#        #MLT=pyf['MLT_Ele'][...]
#        #Lele=pyf['L_Ele'][...]
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
#        dfcFlags=pd.DataFrame(chargingFlags,index=epochEFW,columns=['charging'])
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
#        kelvin=11604.505
#        boltzmann=1.38*1e-23
#        tpt=np.array(pyf['Tpar_e_200'][...])*kelvin
#        tpp=np.array(pyf['Tperp_e_200'][...])*kelvin
#        tpt[tpt>1e24]=np.nan
#        tpp[tpp>1e24]=np.nan
#        hd=np.array(pyf['Dens_e_200'][...])*1e6
#        hd[hd<0]=np.nan
#        Tboth=np.array(tpt)*(2/3.)+np.array(tpp)*(1/3.)
#        Pressure=hd*boltzmann*Tboth/1.0e-9
#        htime=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
#        dfP=pd.DataFrame(Pressure, index=htime, columns=['P'])
#        P=dfP['P'].resample('1min', how='median').reindex(index=rng)
#
#        # interpolate the data frames
#        EFWpot=dfpotential['potential'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
#        LionCL=dfLion['Lion'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
#        cFlags=dfcFlags['charging'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
#        #EFLUX=[ [] for iNN in range(72)]
#        IonCounts=[ [] for iNN in range(EnBins)]
#
#        #
#        # identify times of charging or eclipse
#        # where cFlags > 0 or EFW charging < 0
#        chargingTimes=np.where((np.array(cFlags)>0) | (np.array(EFWpot) < 0))[0]
#        LionC=LionCL[chargingTimes]
#        
#        for iEn in range(72):
#        #    EFLUX[iEn]=list(np.array(dfElectronFlux[iEn]['flux'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan))[chargingTimes])
#            IonCounts[iEn]+=list(np.array(dfIonCounts[iEn]['counts'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan))[chargingTimes])
#            #ElectronFluxData[iEn]+=list(EFLUX[iEn])
#        cTimes=np.zeros(len(rng))
#        cTimes[chargingTimes]=1
#            
#        # now get the ion charging line
#        
#        chargingLine=[]
#        chargingTimesIon=[]
#        pressure=[]
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
#                         pressure.append(P[chargingTimes][iTime])
#                         check=1
#                    except(IndexError):
#                        break
#            if check!=1:
#              # nothing
#              continue
#
#        IonLineData+=list(chargingLine)
#        IonChargingTimes+=list(chargingTimesIon)
#        
#        ElectronPressure+=list(pressure)
#
#        DT=DT+datetime.timedelta(days=1)
#
#     except:
#         print date
#         print len(ElectronPressure)
#         print len(IonLineData)
#         DT=DT+datetime.timedelta(days=1)
#         continue
##
#    # combine both satellites 
#    ElectronPressureC+=ElectronPressure
#    IonLineC+=IonLineData
#    IonChargingTimesC+=IonChargingTimes
#
#
##stop
#import pickle
#
#        
#os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
## plot the results
#IonLineC=np.array(IonLineC)
#tA=np.where((IonLineC < 10) | (IonLineC > 1000))
#IonLineC[tA]=np.nan
#pickle.dump(IonLineC, open('IonAlgoLine.p','wb'))
#pickle.dump(IonChargingTimesC,open('IonAlgoTime.p','wb'))
#pickle.dump(ElectronPressureC, open('ElectronPressureC.p','wb'))
#
#
#stop
IonLineC=pickle.load(open('IonAlgoLine.p','rb'))
IonLineC=np.array(IonLineC)
tA=np.where((IonLineC < 10) | (IonLineC > 1000))
IonLineC[tA]=np.nan
#EFluxC=np.array(pickle.load(open('IonAlgoEflux.p','rb')))/1000.0
ETemp=np.array(pickle.load(open('ElectronPressureC.p', 'rb')))

#os.chdir('ChargingElectronFluxes')
EFW_P=pickle.load(open('EFW_P.p', 'rb'))
kelvin=11604.505
EFW_P=np.array(EFW_P)*kelvin
EFWSC=pickle.load(open('EFWSC_P.p', 'rb'))
#os.chdir('..')
#Numbers=pickle.load(open('binnedelectronfluxes.p', 'rb'))
#dataNumbers=[]
#for iNum in range(len(Numbers)):
#      dataNumbers+=[bins1[iNum]]*Numbers[iNum]
#
# run a correlation coefficient
#
# remove nans and combine HOPE + EFW
cP=list(np.array(EFW_P[~np.isnan(EFW_P)])) + list(ETemp[~np.isnan(IonLineC)])
cSC=list(-1*np.array(EFWSC)[~np.isnan(EFW_P)]) + list(IonLineC[~np.isnan(IonLineC)])
#
# remove unphysical values

#cP=list(ETemp[~np.isnan(IonLineC)])
#cSC=list(IonLineC[~np.isnan(IonLineC)])


#cSC=cSC[cP>0]; cP=cP[cP > 0]
# determine best fit line
#cSC=np.array(cSC); cP=np.array(cP)
#cSC[cP>100]=np.nan
#cP[cP>100]=np.nan
#cSC[cP<.3]=np.nan
#cP[cP<.3]=np.nan
#cP[cSC<1]=np.nan
#cSC[cSC<1]=np.nan
#cSC=np.array(cSC); cP=np.array(cP)*100
#cP=np.log10(cP); cSC=np.log10(cSC)
#cSC=cSC[~np.isnan(cP)]; cP=cP[~np.isnan(cP)]
#cP=cP[~np.isnan(cSC)]; cSC=cSC[~np.isnan(cSC)]
#par = np.polyfit(cP, cSC, 1, full=True)
#
#slope=par[0][0]
#intercept=par[0][1]
#print slope
#print intercept
#xl = [min(cP),max(cP)]
#yl = [slope*xx + intercept  for xx in xl]
#variance = np.var(cSC)
#residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(cP,cSC)])
#Rsqr = np.round(1-residuals/variance, decimals=2)

#seriesB = np.transpose(np.log10(cP[cSC>1]))
#seriesA = np.transpose(np.log10(cSC[cSC>1]))
#
#ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
#tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
#spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
#print("Kendall's tau for ( {0} (p={1})".format(tau,p_value))
#print("Spearman's R for ( is {0} (p={1})".format(spear,spear_p))
  




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

axScatter.plot(np.array(EFW_P)*100,-1*np.array(EFWSC), '.', color='darkorange')
axScatter.plot(np.array(ETemp)*100,IonLineC, '.', color='red')
axScatter.set_xlim(1e0,1e3)
axScatter.set_ylim(1e0,1e3)
axScatter.set_yscale('log')
axScatter.set_xscale('log')
plt.draw()
font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
axScatter.set_ylabel('Spacecraft Potential [V]', fontweight='bold', fontsize=22)
axScatter.set_xlabel('Pressure [nPa]',fontweight='bold', fontsize=22)
#axHistx.hist(np.array(dataNumbers[iEn])/1000.0, bins=bins1, stacked=True, color='blue')
#axHistx.hist(np.array(EFW_T), bins=bins1,stacked=True, color='blue')
#axHistx.set_xlim( axScatter.get_xlim() )
#axHistx.set_yscale('log')
#axHistx.set_xscale('log')
axScatter.tick_params(axis='both', which='major', labelsize=22)
#xx=np.logspace(-2, 2, 100)
#axScatter.plot(10**np.array(xl), 10**np.array(yl), '-', lw=2, c='k')
#axScatter.annotate('r='+str(np.round(spear*100)/100.0), xy=( 1, 10000), xycoords='data', color='k')
#print Rsqr
#plt.text(.9*(10**max(cP))+.1*(10**min(cP)),.9*(10**max(cSC))+.1*10**(min(cSC)),'$R^2 = %0.2f$'% Rsqr, fontsize=30)
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
           os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)
plt.savefig('Both_EFWPressure_Energy.png')
plt.close()
os.chdir('..')
