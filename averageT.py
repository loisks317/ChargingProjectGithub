# averageT.py
#
# average temperature plots over 6 month period 
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
import scipy.stats
import pandas as pd
import bisect
#
# start the date loop
dateStart='20130201'
#dateEnd='20130118'
dateEnd='20130601'
#dateEnd='20150401'
ExcludeDates=['20130206', '20130227', '20130327', '20130417', '20130508']
# data
hourGap=5
nx=101
threshold1=-1
threshold2=-25
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
sat=['A']

Tavg=[]
pot=[]

# get the low ion Plasma data
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/PlasmaDensity200eV')
lowIonsData=[]
lowIonsTime=[]
datesArr=['2013_2','2013_3','2013_4','2013_5']
for iDate2 in range(len(datesArr)):
   lowIonsData+=list(pickle.load(open("P_sat=rbspA_"+datesArr[iDate2], "rb")))
   lowIonsTime+=list(pickle.load(open("time_P_sat=rbspA_"+datesArr[iDate2], "rb")))
pdTime=pd.to_datetime(np.array(lowIonsTime))
lowITime=[ datetime.datetime.strptime(str(i), '%Y-%m-%d %H:%M:%S') for i in pdTime]
print "converted the ion Plasma Times"
while DT != endDt:
 date=datetime.datetime.strftime(DT, '%Y%m%d')
 DT=datetime.datetime.strptime(date, '%Y%m%d')
 check=0
 for iDate in ExcludeDates:
     if iDate == date:
      check=1
 if check == 0:   
  for iSat in range(len(sat)):    
    #
    # load in EFW data
    # open the files
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+sat[iSat])
    f=glob.glob('*'+date+'*')
    try:
        pyf=pycdf.CDF(f[0])
    except:
        continue
    efwData={}
    dd=(list(pyf['density'][...]))
    dd[dd<0]=np.nan
    efwData['density']=dd
    efwData['potential']=pyf['Vavg'][...]
    efwData['epochE']=pd.DatetimeIndex(pyf['epoch'][...])
    dfepotential=pd.DataFrame( efwData['potential'], index=efwData['epochE'], columns=['potential'])
    rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
    rng=rt[::1]

        #
    # now get HOPE data
    #
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_MOM_'+sat[iSat])
    f=glob.glob('*'+date+'*')
    try: 
        pyf=pycdf.CDF(f[0])
    except:
        continue

    hopeData={}
    tpt=pyf['Tpar_e_200'][...]
    tpt[tpt>1e24]=np.nan
    hopeData['Tpar']=tpt
    tpp=pyf['Tperp_e_200'][...]
    tpp[tpp>1e24]=np.nan
    hopeData['Tperp']=tpp
    Tpre=(1/3.*hopeData['Tpar'] + 2/3.*hopeData['Tperp'])

    
    hd=pyf['Dens_e_200'][...]
    hd[hd<0]=np.nan
    hopeData['density']=hd*1.0e6
    hopeData['epoch']=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
    dfT=pd.DataFrame(Tpre, index=hopeData['epoch'], columns=['T'])
    dfdensity=pd.DataFrame(hopeData['density'], index=hopeData['epoch'], columns=['density'])
    T=dfT['T'].resample('1min', how='median').reindex(index=rng)
    density=dfdensity['density'].resample('1min', how='median').reindex(index=rng)
    potential=np.array(dfepotential['potential'].resample('1min', how='median').reindex(index=rng))



    ## one more, get low ion density
    # get the right times for Ion Density
    try:
      ionBiRight=bisect.bisect_right(lowITime, DT)
      ionBiLeft=bisect.bisect_left(lowITime,DT+datetime.timedelta(days=1))
      ionDensitydf=pd.DataFrame(lowIonsData[ionBiRight:ionBiLeft], index=lowITime[ionBiRight:ionBiLeft], columns=['density'])
      #for iDen in range(len(ionDensitydf['density'])):
      # try:
      #  a=len(ionDensitydf['density'][iDen])
      #  ionDensitydf['density'][iDen]=np.nan
      # except(TypeError):
      #  ionDensitydf['density'][iDen]=ionDensitydf['density'][iDen].astype(float)    
      ionDensity=ionDensitydf['density'].resample('1min', how='median').reindex(index=rng)
    except(ValueError):
      print "No Hope Data for " + str(DT)
      continue 

    Ttemp=((np.array(ionDensity)*200) + (np.array(density)*np.array(T)))/( np.array(ionDensity)+np.array(density))
    neg=np.where(np.array(potential)>0)
    pot+=list(potential[neg])
    Tavg+=list(Ttemp[neg])

                  
  
  DT=DT+datetime.timedelta(days=1)
 else:
  DT=DT+datetime.timedelta(days=1)

# pickle import here 


# now make scatter plots
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/')
import pickle
IonSC=pickle.load(open('IonAlgoLineTavg.p', 'rb'))
ITavg=pickle.load(open('ElectronTavgC.p', 'rb'))

# stats on ion
IonSC=np.array(IonSC)
ITavg=np.array(ITavg)
seriesB = np.transpose(np.log10(ITavg[IonSC>10]))
seriesA = np.transpose(np.log10(IonSC[IonSC>10]))
ts_res = scipy.stats.theilslopes(seriesA, seriesB, 0.95)
tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
print("Kendall's tau for  is {0} (p={1})".format(tau,p_value))
print("Spearman's R for  is {0} (p={1})".format(spear,spear_p))



left, width = 0.15, 0.8
bottom, height = 0.15, 0.8
left_h = left+width+0.03
bottom_h = bottom+height+0.03

rect_scatter = [left, bottom, width, height]
fig=plt.figure(1, figsize=(8,8))
axScatter = plt.axes(rect_scatter)
axScatter.plot(np.array(Tavg),pot, '.', color='darkorange')
axScatter.plot(np.array(ITavg), IonSC, '.', color='red')
xx=np.array([np.nanmin(seriesB),4.25])
#axScatter.plot(10**np.array(xx), 10**np.array(ts_res[1] + ts_res[0] * xx),'-', linewidth=2, color='k')
axScatter.set_xlim(1e3,1e4)
axScatter.set_ylim(1e-1,1e3)
axScatter.set_yscale('log')
axScatter.set_xscale('log')
plt.draw()
font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
axScatter.set_ylabel('Spacecraft Potential [V]', fontweight='bold', fontsize=22)
axScatter.set_xlabel('<T$_{e}$>',fontweight='bold', fontsize=22)
axScatter.tick_params(axis='both', which='major', labelsize=22)
subdir_name='NewIonStats'
if not os.path.exists(subdir_name):
           os.umask(0) # unmask if necessary
           os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)
plt.savefig('Both_AveT.png')
plt.close()
os.chdir('..')




            
               

           
        
         
### DEBUG THIS IN MORNING
