# scatterPT.py
#
# two variable scatter plots with > -25 and < -25 V of charging by color
# maybe do roxanne style plot normalized by color
#
#
# LKS, August 2015
#
# imports
from __future__ import print_function
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
import pandas as pd

#
# start the date loop
dateStart='20130101'
#dateEnd='20130118'
dateEnd='20130601'
#dateEnd='20150401'
# data
hourGap=5
nx=101
threshold1=-1
threshold2=-25

sat=['A', 'B']

#
# thresholds

densityT=[1.0*1e5, 1.5*1e5, 2.0*1e5,2.5*1e5, 3.0*1e5, 3.5*1e5, 4.0*1e5]
LT=5
MLTlowT=0
MLThighT=6
TperpT=[.9*1e7, 1.0*1e7, 1.5*1e7, 2.0*1e7, 2.5*1e7]
pressureT=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
TparT=[.9*1e7, 1.0*1e7, 1.5*1e7, 2.0*1e7, 2.5*1e7]



dateArr=[]
for iDen in range(len(densityT)):
    for iTperpT in range(len(TperpT)):
        for iTparT in range(len(TparT)):
         for iPress in range(len(pressureT)):
                            
             potentialm=[]
             checkP=[]
             potentialhigh=[]
             potentialall=[]
             date=dateStart
             endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
             DT=datetime.datetime.strptime(date, '%Y%m%d')
             while DT != endDt:
               date=datetime.datetime.strftime(DT, '%Y%m%d')
               DT=datetime.datetime.strptime(date, '%Y%m%d')
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
                 pos=np.where(np.array(efwData['potential']) > 10)[0] # eclipse
                 efwData['potential'][pos]=np.nan
                 dfepotential=pd.DataFrame( efwData['potential'], index=efwData['epochE'], columns=['potential'])
                 dfedensity=pd.DataFrame( efwData['density'], index=efwData['epochE'], columns=['density'])
                 dfeMLT=pd.DataFrame(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0], index=efwData['epochE'], columns=['MLT'])
                 dfeMLAT=pd.DataFrame(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2], index=efwData['epochE'], columns=['MLAT'])
                 dfeL=pd.DataFrame(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1], index=efwData['epochE'], columns=['L'])
                 rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
                 rng=rt[::5]
                 eDensity=dfedensity['density'].resample('5min', how='median').reindex(index=rng)
                 potential=dfepotential['potential'].resample('5min', how='median').reindex(index=rng, fill_value=np.nan)
                 MLT=dfeMLT['MLT'].resample('5min', how='median').reindex(index=rng, fill_value=np.nan)
                 MLAT=dfeMLAT['MLAT'].resample('5min', how='median').reindex(index=rng, fill_value=np.nan)
                 L=dfeL['L'].resample('5min', how='median').reindex(index=rng, fill_value=np.nan)
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
                 Tpar=dfTpar['Tpar'].resample('5min', how='median').reindex(index=rng)
                 Tperp=dfTperp['Tperp'].resample('5min', how='median').reindex(index=rng)
                 density=dfdensity['density'].resample('5min', how='median').reindex(index=rng)
                 pressure=dfpressure['pressure'].resample('5min', how='median').reindex(index=rng)
                 #
             
                 
             
                 # check for thresholds
                 MLT=np.array(MLT)
                 L=np.array(L)
                 density=np.array(density)
                 Tpar=np.array(Tpar)
                 Tperp=np.array(Tperp)
                 pressure=np.array(pressure)
                 MLTi=np.where((MLT > MLTlowT) & (MLT < MLThighT))[0]
                 Li=np.where(L > LT)[0]
                 Di=np.where(density > densityT[iDen])[0]
                 Pi=np.where(pressure > pressureT[iPress])[0]
                 Tpari=np.where(Tpar > TparT[iTparT])[0]
                 Tperpi=np.where(Tperp > TperpT[iTperpT])[0]
                 #
                 #matches=[i for i, j,k,l,m in zip(MLTi,Di,Tpari,Tperpi,Pi) if i==j==k==l==m]
                 matches=list(set(MLTi) & set(Di) & set(Tpari) & set(Tperpi) & set(Pi) & set(Li))
                 potentialm+=list(potential[matches])
             
                 # now check how we are doing in catch the extreme charging events with these parameters
                 phi=np.where(np.array(potential) < -25)[0]
                 kai=np.where(np.array(potential) < -1][0]
                 potentialall+=list(potential[kai])
                 potentialhigh+=list(potential[phi])
             
                 MLT=np.array(MLT)
                 L=np.array(L)
                 density=np.array(density)
                 Tpar=np.array(Tpar)
                 Tperp=np.array(Tperp)
                 pressure=np.array(pressure)
                 MLTe=np.where((MLT[phi] > MLTlowT) & (MLT[phi] < MLThighT))[0]
                 Le=np.where(L[phi] > LT)[0]
                 De=np.where(density[phi] > densityT[iDen])[0]
                 Pe=np.where(pressure[phi] > pressureT[iPress])[0]
                 Tpare=np.where(Tpar[phi] > TparT[iTparT])[0]
                 Tperpe=np.where(Tperp[phi] > TperpT[iTperpT])[0]
                 #
                 #matches=[i for i, j,k,l,m in zip(MLTi,Di,Tpari,Tperpi,Pi) if i==j==k==l==m]
                 matches=list(set(MLTe) & set(De) & set(Tpare) & set(Tperpe) & set(Pe) & set(Le))
                 checkP+=list(potential[phi][matches])

                               
               
               DT=DT+datetime.timedelta(days=1)
             os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
             print("got to here, check file" )
             log=open('sensitivityTestOutput.txt', 'a')
             print( "for Tperp="+str(TperpT[iTperpT])+", Tpar="+str(TparT[iTparT])+", Pressure="+str(pressureT[iPress])+", Density="+str(densityT[iDen]), file=log)
             print( "total identifiers = " + str(len(potentialm)), file=log)
             try:
                print ("Success on high potential catches is: " + str(len(checkP)/(1.0*len(potentialhigh))), file=log)
             except:
                print ('bad',file=log)
                                            
             
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')

