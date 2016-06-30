# twoVarScatsLow.py
#
# two variable scatter plots with > -25 of charging by color
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
import scipy.stats
#
# start the date loop
dateStart='20130101'
#dateEnd='20130118'
dateEnd='20130601'
# data
hourGap=5
nx=101
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
sat=['A', 'B']
overallPotential=[]
overallTpar=[]
overallTperp=[]
overallPressure=[]
overallDensity=[]
dateArr=[]
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
    density=(list(pyf['density'][...]))
    potential=-1.0*pyf['Vavg'][...]
    epochE=dates.date2num(pyf['epoch'][...])
    MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
    LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
    MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
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
    Tpar=pyf['Tpar_e_200'][...]*kelvin
    Tperp=pyf['Tperp_e_200'][...]*kelvin
    density=pyf['Dens_e_200'][...]
    density[density<0]=np.nan
    Tperp[Tperp>1e24]=np.nan
    Tpar[Tpar>1e24]=np.nan
    epoch=pyf['Epoch_Ele'][...]
    pressure=density*boltzmann*(1/3.*Tpar + 2/3.*Tperp) # check this
    

    # interpolate onto potential times
    # calculate correlation coefficient for s/c charging and temperature
    os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
    mTime=dates.date2num(epoch)

                                   # define a charging event
    lowPotential=np.where(potential > -25)[0]
    try:
     lowPotential[0]=-100000 # choose an absurdly low value so it doesn't match
     # now sequence
     potentialEvents={}
     count=-1
     for i in range(1,len(lowPotential)):
        if lowPotential[i] != lowPotential[i-1]+1:
            # out of sequence
            count+=1
            potentialEvents[count]=[lowPotential[i]]
        else:
            potentialEvents[count].extend([lowPotential[i]])
    
     # Now calculate a correlation coefficient in the ones left
    
     # now test to make sure each event is long enough
     iCount=0
     while iCount < count:
        interpPressure=[]
        interpTpar=[]
        interpTperp=[]
        interpDensity=[]
        if len(potentialEvents[iCount]) < 10:
            print "Too small of an event"
            potentialEvents.pop(iCount,0) # this should remove the item
            # doesn't change the count though, will need to account
            # for that later
            
    
    #
    # Now calculate a correlation coefficient in the ones left
        else:
            f=interp1d(epochE[potentialEvents[iCount]],potential[potentialEvents[iCount]])
            # now find where mTime is greater than epochE[startEvent]
            AdeleIsAwesome=np.where(mTime > epochE[potentialEvents[iCount]][0])[0][0]
            AdeleIsAwesomer=np.where(mTime > epochE[potentialEvents[iCount]][-1])[0][0]-1
            try:
                       tempC=f(mTime[AdeleIsAwesome:AdeleIsAwesomer])
            except(ValueError):
                       # gah fuck it
                       print "you fucked it up girl" 
                       mTime[mTime<epochE[potentialEvents[iCount]][0]]=epochE[potentialEvents[iCount]][0]
                       mTime[mTime>epochE[potentialEvents[iCount]][-1]]=epochE[potentialEvents[iCount]][-1]
                       tempC=f(mTime[AdeleIsAwesome:AdeleIsAwesomer])
            allPotential=list(tempC)
            dateArr.append(date)
            interpTpar=list(Tpar[AdeleIsAwesome:AdeleIsAwesomer])
            interpTperp=list(Tperp[AdeleIsAwesome:AdeleIsAwesomer])
            interpP=list(pressure[AdeleIsAwesome:AdeleIsAwesomer])
            interpDensity=list(density[AdeleIsAwesome:AdeleIsAwesomer])
            try:
                    interpTpar=list(interpTpar[~np.isnan(interpTpar)])
                    interpTperp=list(interpTperp[~np.isnan(interpTperp)])
                    interpP=list(interpP[~np.isnan(interpP)])
                    interpDensity=list(interpDensity[~np.isnan(interpDensity)])
            except:
                    interpTpar=list(interpTpar)
                    interpTperp=list(interpTperp)
                    interpP=list(interpP)
                    interpDensity=list(interpDensity)
            # add them all
            overallTpar+=interpTpar
            overallTperp+=interpTperp
            overallPressure+=interpP
            overallDensity+=interpDensity

        iCount+=1        
    except(IndexError):
        # nothing happens
        print 'no events on ' + date
  DT=DT+datetime.timedelta(days=1)

# now make scatter plots
aList=[overallTpar, overallTperp, overallPressure, overallDensity]
names=['Tpar', 'Tperp', 'Pressure', 'Density']
for i in range(len(aList)):
    for j in range(len(aList)):
        # define bins
           boxA=np.logspace(np.nanmin(np.log10(aList[i])), np.nanmax(np.log10(aList[i])), nx)
           boxB=np.logspace(np.nanmin(np.log10(aList[j])), np.nanmax(np.log10(aList[j])), nx)
           diffA=np.nanmax(aList[i])- np.nanmin(aList[i])
           diffB=np.nanmax(aList[j])-np.nanmin(aList[j])
           data=[[0 for s in range(nx-1)] for r in range(nx-1)]
           #dataB=[0 for s in range(1000)] 
           for ii in range(nx-1):
            for jj in range(nx-1):
               try:
                   index=np.where((aList[i] >= boxA[ii]) & (aList[i] < boxA[ii+1]))[0]
               except(ValueError):
                   index = []
               try:
                   indexj=np.where((aList[j] >= boxB[jj]) & (aList[j] < boxB[jj+1]))[0]
               except(ValueError):
                   indexj=[]
               #dataA[ii]=len(index)
               data[ii][jj]=len(index)+len(indexj)
                   
           #dataA=np.array(dataA)
           #dataB=np.array(dataB)
           #data=np.vstack((dataA,dataB)) # THIS IS NOT RIGHT
           #
           fig=plt.figure()
           ax1=fig.add_subplot(111)
           plt.subplots_adjust(left=0.11, right=0.75, top=0.92, bottom=0.4)
           X,Y=np.meshgrid(boxA, boxB)
           ax1.set_xlim(boxA[0], boxA[-1])
           ax1.set_ylim(boxB[0], boxB[-1])
           ax1.set_xscale('log')
           ax1.set_yscale('log')
           ax1.set_xlabel(names[i])
           ax1.set_ylabel(names[j])
           col=ax1.pcolormesh(X[:nx-1],Y[:nx-1],np.log10(np.array(data)),cmap='jet', vmin=0, vmax=6)
           font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}
           plt.rc('font', **font)
           cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
           cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(0,6))
           cb.set_label('Number of Points', fontsize=25, fontweight='bold')
           ax1.tick_params(axis='both', which='major', labelsize=22)
           cb.ax.tick_params(labelsize=30)
           subdir_name='DoubleScatter'
           if not os.path.exists(subdir_name):
              os.umask(0) # unmask if necessary
              os.makedirs(subdir_name, 0777) 
           os.chdir(subdir_name)#
           fig.set_size_inches(13,9)
           plt.savefig(names[i]+ '_'+names[j]+'_low.pdf')
           plt.close(fig)
           os.chdir('..')
            
               

           
        
         
### DEBUG THIS IN MORNING
