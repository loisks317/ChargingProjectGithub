# scatterPT3d.py
#
# three variable scatter plots with > -25 and < -25 V of charging by color
#
#
# LKS, September  2015
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
import pylab
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
from mpl_toolkits.mplot3d import Axes3D

dateStart='20130201'
dateEnd='20130210'
#dateEnd='20130118'
#dateEnd='20130601'
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
lowPotential=[]
lowTpar=[]
lowTperp=[]
lowPressure=[]
lowDensity=[]
lowEden=[]
lowMLT=[]
lowL=[]
lowMLAT=[]
lowIDen=[]
mPotential=[]
mTpar=[]
mTperp=[]
mPressure=[]
mDensity=[]
mEden=[]
mIDen=[]
mMLT=[]
mL=[]
mMLAT=[]
highPotential=[]
highTpar=[]
highTperp=[]
highPressure=[]
highDensity=[]
highEden=[]
highMLT=[]
highL=[]
highMLAT=[]
highIDen=[]
dateArr=[]

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
    efwData['potential']=-1.0*pyf['Vavg'][...]
    efwData['epochE']=pd.DatetimeIndex(pyf['epoch'][...])
    dfepotential=pd.DataFrame( efwData['potential'], index=efwData['epochE'], columns=['potential'])
    dfedensity=pd.DataFrame( efwData['density'], index=efwData['epochE'], columns=['density'])
    dfeMLT=pd.DataFrame(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0], index=efwData['epochE'], columns=['MLT'])
    dfeMLAT=pd.DataFrame(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2], index=efwData['epochE'], columns=['MLAT'])
    dfeL=pd.DataFrame(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1], index=efwData['epochE'], columns=['L'])
    rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
    rng=rt[::1]
    eDensity=dfedensity['density'].resample('1min', how='median').reindex(index=rng)
    potential=dfepotential['potential'].resample('1min', how='median').reindex(index=rng, fill_value=np.nan)
    MLT=dfeMLT['MLT'].resample('1min', how='median').reindex(index=rng, fill_value=np.nan)
    MLAT=dfeMLAT['MLAT'].resample('1min', how='median').reindex(index=rng, fill_value=np.nan)
    L=dfeL['L'].resample('1min', how='median').reindex(index=rng, fill_value=np.nan)
        #
    # now get HOPE data
    #
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_MOM_'+sat[iSat])
    f=glob.glob('*'+date+'*')
    try: 
        pyf=pycdf.CDF(f[0])
    except:
        continue

    kelvin=11604.505
    boltzmann=1.38*1e-23
    hopeData={}
    tpt=pyf['Tpar_e_200'][...]*kelvin
    tpt[tpt>1e24]=np.nan
    hopeData['Tpar']=tpt
    tpp=pyf['Tperp_e_200'][...]*kelvin
    tpp[tpp>1e24]=np.nan
    hopeData['Tperp']=tpp
    hd=pyf['Dens_e_200'][...]
    hd[hd<0]=np.nan
    hopeData['density']=hd*1.0e6
    hopeData['epoch']=pd.DatetimeIndex(pyf['Epoch_Ele'][...])
    hopeData['pressure']=(hopeData['density']*boltzmann*(1/3.*hopeData['Tpar'] + 2/3.*hopeData['Tperp']))/1.0e-9 # check this
    dfTpar=pd.DataFrame(hopeData['Tpar'], index=hopeData['epoch'], columns=['Tpar'])
    dfTperp=pd.DataFrame(hopeData['Tperp'], index=hopeData['epoch'], columns=['Tperp'])
    dfdensity=pd.DataFrame(hopeData['density'], index=hopeData['epoch'], columns=['density'])
    dfpressure=pd.DataFrame(hopeData['pressure'], index=hopeData['epoch'], columns=['pressure'])
    Tpar=dfTpar['Tpar'].resample('1min', how='median').reindex(index=rng)
    Tperp=dfTperp['Tperp'].resample('1min', how='median').reindex(index=rng)
    density=dfdensity['density'].resample('1min', how='median').reindex(index=rng)
    pressure=dfpressure['pressure'].resample('1min', how='median').reindex(index=rng)



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

    lowIndex=np.where(np.array(potential)>threshold1)[0]
    mIndex=np.where(np.array(potential)<=threshold1)[0]
    highIndex=np.where(np.array(potential)<=threshold2)[0]
    
    lowPotential+=list(potential[lowIndex])
    lowTpar+=list(Tpar[lowIndex])
    lowTperp+=list(Tperp[lowIndex])
    lowPressure+=list(pressure[lowIndex])
    lowDensity+=list(density[lowIndex])
    lowEden+=list(eDensity[lowIndex])
    lowMLT+=list(MLT[lowIndex])
    lowMLAT+=list(MLAT[lowIndex])
    lowL+=list(L[lowIndex])
    lowIDen+=list(ionDensity[lowIndex])

    mPotential+=list(potential[mIndex])
    mTpar+=list(Tpar[mIndex])
    mTperp+=list(Tperp[mIndex])
    mPressure+=list(pressure[mIndex])
    mDensity+=list(density[mIndex])
    mEden+=list(eDensity[mIndex])
    mMLT+=list(MLT[mIndex])
    mMLAT+=list(MLAT[mIndex])
    mL+=list(L[mIndex])
    mIDen+=list(ionDensity[mIndex])

    highPotential+=list(potential[highIndex])
    highTpar+=list(Tpar[highIndex])
    highTperp+=list(Tperp[highIndex])
    highPressure+=list(pressure[highIndex])
    highDensity+=list(density[highIndex])
    highEden+=list(eDensity[highIndex])
    highMLT+=list(MLT[highIndex])
    highMLAT+=list(MLAT[highIndex])
    highL+=list(L[highIndex])
    highIDen+=list(ionDensity[highIndex])
                  
  
  DT=DT+datetime.timedelta(days=1)
 else:
  DT=DT+datetime.timedelta(days=1)


fig=pylab.figure()
ax=Axes3D(fig)
ax.scatter(lowEden[lowEden>0], lowTpar[lowTpar>0], lowTperp[lowTperp>0], c='b')
ax.scatter(mEden[mEden>0], mTpar[mTpar>0], mTperp[mTperp>0], c='magenta')
ax.scatter(highEden[highEden>0], highTpar[highTpar>0], highTperp[highTperp>0], c='gold')
#ax.set_xlim(1e4,1e7 )
#ax.set_ylim(1e6,1e8)
#ax.set_zlim(1e6,1e8)
ax.xaxis.set_scale('log')
ax.yaxis.set_scale('log')
ax.zaxis.set_scale('log')
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/DoubleScatter3d')
pylab.savefig('3var_Compare.png')
os.chdir('..')

