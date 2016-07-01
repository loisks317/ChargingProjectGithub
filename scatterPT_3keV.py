# scatterPT.py
#
# two variable scatter plots with > -25 and < -25 V of charging by color
# maybe do roxanne style plot normalized by color
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
#
eC='3keV'
eChan=46 # 46 = 3 keV
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
lowPotential=[]
lowTpar=[]
lowTperp=[]
lowPressure=[]
lowDensity=[]
lowEden=[]
lowMLT=[]
lowTC=[]
lowL=[]
lowMLAT=[]
lowIDen=[]
mPotential=[]
mTpar=[]
mTperp=[]
mTC=[]
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
highTC=[]
highL=[]
highMLAT=[]
highIDen=[]
dateArr=[]

# get the low ion Plasma data
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/PlasmaDensity200eV')
lowIonsData=[]
lowIonsTime=[]
datesArr=['2013_2','2013_3','2013_4','2013_5']
os.chdir('/Users/loisks/Desktop/ResearchProjects/HOPEPlasmaDensities')
lowIonsData+=list(pickle.load(open("1_210eV_corrected_plasma_density_0_36_sat=A_species=P", "rb")))
lowIonsTime+=list(pickle.load(open("1_210eV_corrected_time_0_36_sat=A_species=P", "rb")))   
pdTime=pd.to_datetime(np.array(lowIonsTime))
lowITime=[ datetime.datetime.strptime(str(i), '%Y-%m-%d %H:%M:%S') for i in pdTime]
print("converted the ion Plasma Times")
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
    dd=np.array(pyf['density'][...])
    dd[dd<0]=np.nan
    efwData['density']=list(dd)
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
      print("No Hope Data for " + str(DT))
      continue

# now determine where the 3 keV electron fluxes are above a certain threshold
# 3* 1e10
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+sat[iSat])
    f=glob.glob('*'+date+'*')
    try: 
        pyfE=pycdf.CDF(f[0])
    except:
        print("No Hope Data for " + str(DT) )      
        continue
    epochE=pd.DatetimeIndex(pyfE['Epoch_Ele'][...])
    EE=pyfE['HOPE_ENERGY_Ele'][...]
    # 53 = 3.3 keV
    ElectronFlux=np.array(np.swapaxes(pyfE['FEDO'][...],1,0)*np.swapaxes(EE,1,0))[eChan]
    # resample
    

    dfEflux=pd.DataFrame(ElectronFlux, index=epochE, columns=['Eflux'])
    Eflux=np.array(dfEflux['Eflux'].resample('1min', how='median').reindex(index=rng))
    highFlux=np.where(np.array(Eflux) > 3*1e10)[0]   
    
    lowIndex=set(np.where(np.array(potential)>threshold1)[0]) & set(highFlux)
    mIndex=set(np.where(np.array(potential)<=threshold1)[0]) & set(highFlux)
    highIndex=set(np.where(np.array(potential)<=threshold2)[0]) & set(highFlux)


    lowPotential+=list(potential[lowIndex])
    lowTpar+=list(Tpar[lowIndex])
    lowTperp+=list(Tperp[lowIndex])
    lowTC+=list(((1/3.0)*Tpar[lowIndex]) + ((2/3.0)*Tperp[lowIndex]))
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
    mTC+=list(((1/3.0)*Tpar[mIndex]) + ((2/3.0)*Tperp[mIndex]))
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
    highTC+=list(((1/3.0)*Tpar[highIndex]) + ((2/3.0)*Tperp[highIndex]))
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

# binall the data
#binNo=500
#TparBinsLow=np.logspace(np.log10(np.nanmin(lowTpar)), np.log10(np.nanmax(lowTpar)), binNo)
#TperpBinsLow=np.logspace(np.log10(np.nanmin(lowTperp)), np.log10(np.nanmax(lowTperp)), binNo)
#PressureBinsLow=np.logspace(np.log10(np.nanmin(lowPressure)), np.log10(np.nanmax(lowPressure)), binNo)
#DensityBinsLow=np.logspace(np.log10(np.nanmin(lowDensity)), np.log10(np.nanmax(lowDensity)), binNo)
#TparBinsHigh=np.logspace(np.log10(np.nanmin(highTpar)), np.log10(np.nanmax(highTpar)), binNo)
#TperpBinsHigh=np.logspace(np.log10(np.nanmin(highTperp)), np.log10(np.nanmax(highTperp)), binNo)
#PressureBinsHigh=np.logspace(np.log10(np.nanmin(highPressure)), np.log10(np.nanmax(highPressure)), binNo)
#DensityBinsHigh=np.logspace(np.log10(np.nanmin(highDensity)), np.log10(np.nanmax(highDensity)), binNo)

#TPAlow=np.zeros(binNo)
#TPPlow=np.zeros(binNo)
#Plow=np.zeros(binNo)
#Dlow=np.zeros(binNo)
#TPAhigh=np.zeros(binNo)




# now make scatter plots
os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/')
lowDensity=np.array(lowDensity)/1.0e6
mDensity=np.array(mDensity)/1.0e6
highDensity=np.array(highDensity)/1.0e6
aList=[ lowIDen,lowTpar, lowTperp,lowTC, lowPressure, lowDensity, lowMLT, lowMLAT, lowL]
bList=[highIDen,highTpar, highTperp,highTC, highPressure, highDensity, highMLT, highMLAT, highL]
cList=[mIDen, mTpar, mTperp,mTC,  mPressure, mDensity,mMLT,mMLAT,mL]
#aList=[TparBinsLow, TperpBinsLow, PressureBinsLow, DensityBinsLow]
#bList=[TparBinsHigh, TperpBinsHigh, PressureBinsHigh, DensityBinsHigh]
#aC=[TPAlow, TPPlow, Plow, Dlow]
#xbC=[TPAhigh, TPPhigh, Phigh, Dhigh]
xlimsAlow=[1e-3,1e6, 1e6, 1e6, 1e-3, 1e-2,0, -20,1]
xlimsBhigh=[1e1,1e9, 1e9,1e9, 1e2, 1e2,24,10,7]
names=['Low Ion Density [cm$^{-3}$]','Tpar [K]', 'Tperp [K]', 'T [K]', 'Pressure [nPa]', 'Density [cm$^{-3}$]', 'MLT', 'MLAT [Degrees]', 'L']
for i in range(len(aList)):
    for j in range(len(aList)):
        # define bins
           left, width = 0.1, 0.65
           bottom, height = 0.1, 0.65
           bottom_h = left_h = left+width+0.03

           rect_scatter = [left, bottom, width, height]
           rect_histx = [left, bottom_h, width, 0.18]
           rect_histy = [left_h+0.015, bottom, 0.18, height]
           fig=plt.figure(1, figsize=(8,8))
           axScatter = plt.axes(rect_scatter)
           axHistx = plt.axes(rect_histx)
           axHisty = plt.axes(rect_histy)


           # the scatter plot:
           axScatter.plot(aList[i], aList[j], '.', c='blue', label='low')
           axScatter.plot(cList[i], cList[j], '.', c='magenta')
           axScatter.plot(bList[i], bList[j], '.', c='gold', label='low')


           # now determine nice limits by hand:
           xymax = np.nanmax( [np.nanmax((aList[i])), np.nanmax((aList[j]))] )
           xymin = np.nanmin( [np.nanmin((aList[i])), np.nanmin((aList[j]))] )
           


           axScatter.set_xlabel(names[i], fontsize=22, fontweight='bold')
           axScatter.set_ylabel(names[j], fontsize=22, fontweight='bold')
           axScatter.set_xlim(xlimsAlow[i], xlimsBhigh[i])
           axScatter.set_ylim(xlimsAlow[j], xlimsBhigh[j])
           if i < 6:
             bins1 = np.logspace(np.log10(xlimsAlow[i]), np.log10(xlimsBhigh[i]),50)
             axScatter.set_xscale('log')
           else:
             bins1 = np.linspace(xlimsAlow[i], xlimsBhigh[i],50)
             axScatter.set_xscale('linear')
           if j <6:
             bins2 = np.logspace(np.log10(xlimsAlow[j]), np.log10(xlimsBhigh[j]),50)
             axScatter.set_yscale('log')
           else:
             bins2=np.linspace(xlimsAlow[j], xlimsBhigh[j], 50)
             axScatter.set_yscale('linear')
           axHistx.hist(aList[i], bins=bins1, stacked=True, color='blue')
           axHistx.hist(cList[i], bins=bins1, stacked=True, color='magenta')
           axHistx.hist(bList[i], bins=bins1, stacked=True, color='gold')
           axHisty.hist(aList[j], bins=bins2, stacked=True,orientation='horizontal', color='blue')
           axHisty.hist(cList[j], bins=bins2, stacked=True, orientation='horizontal', color='magenta')
           axHisty.hist(bList[j], bins=bins2, stacked=True, orientation='horizontal', color='gold')
           axHistx.set_xlim( axScatter.get_xlim() )
           axHisty.set_ylim( axScatter.get_ylim() )
           if i < 6:
             axHistx.set_xscale('log') 
           if j < 6:
             axHisty.set_yscale('log')
           #if (i < 4) & (j<4):
           # should be for the labels
           axHisty.set_xscale('log')
           axHistx.set_yscale('log')
          
           # work on this later!

           plt.legend()
           font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}
           plt.rc('font', **font)
           #cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
           #cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(0,4))
           #cb.set_label('Number of Points', fontsize=25, fontweight='bold')
           axScatter.tick_params(axis='both', which='major', labelsize=22)
           axHistx.tick_params(axis='x', bottom='off', top='off', labelbottom='off', labeltop='off')
           axHisty.tick_params(axis='y', bottom='off', top='off', labelbottom='off', labeltop='off', labelleft='off')
       
           #cb.ax.tick_params(labelsize=30)
           subdir_name='3keVDoubleScatterNew'
           if not os.path.exists(subdir_name):
              os.umask(0) # unmask if necessary
              os.makedirs(subdir_name) 
           os.chdir(subdir_name)#
           fig.set_size_inches(13,9)
           plt.savefig(names[i]+ '_'+names[j]+'_threshold='+str(threshold1)+'_'+str(threshold2)+'_hist+'+eC+'.png')
           plt.close(fig)
           os.chdir('..')
            
               

           
        
         
### DEBUG THIS IN MORNING
