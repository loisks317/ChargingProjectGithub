# TPstats.py
#
# look at temperature and pressure statistics relative to s/c charging events
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
dateEnd='20150401'
# data
hourGap=5
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
sat=['A', 'B']
corrCoeffTpar=[]
corrCoeffTperp=[]
corrCoeffP=[]
corrCoeffMultiple=[]
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
    epoch=pyf['Epoch_Ele'][...]
    pressure=density*boltzmann*(1/3.*Tpar + 2/3.*Tperp) # check this
    

    # interpolate onto potential times
    # calculate correlation coefficient for s/c charging and temperature
    os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
    mTime=dates.date2num(epoch)

                                   # define a charging event
    lowPotential=np.where(potential < -25)[0]
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
            try:
                    interpTpar=np.array(interpTpar[~np.isnan(interpTpar)])
                    interpTperp=np.array(interpTperp[~np.isnan(interpTperp)])
                    interpP=np.array(interpP[~np.isnan(interpP)])
            except:
                    interpTpar=np.array(interpTpar)
                    interpTperp=np.array(interpTperp)
                    interpP=np.array(interpP)
            # TPar
            try:
                 seriesA = np.transpose(interpTpar)
                 seriesT = np.transpose(allPotential)
                 ts_res = scipy.stats.theilslopes(seriesA, seriesT, 0.95)
                 tau, p_value = scipy.stats.kendalltau(seriesT, seriesA)
                 spear, spear_p = scipy.stats.kendalltau(seriesT, seriesA)
                 corrCoeffTpar.append(spear)

            except(ValueError):
                    print 'Weird Error, continue'
            # Tperp
            try:
                 seriesB = np.transpose(interpTperp)
                 seriesT = np.transpose(allPotential)
                 ts_res = scipy.stats.theilslopes(seriesB, seriesT, 0.95)
                 tau, p_value = scipy.stats.kendalltau(seriesT, seriesB)
                 spear, spear_p = scipy.stats.kendalltau(seriesT, seriesB)
                 corrCoeffTperp.append(spear)

            except(ValueError):
                    print 'Weird Error, continue'

            # Pressure
            try:
                 seriesC = np.transpose(interpP)
                 seriesT = np.transpose(allPotential)
                 ts_res = scipy.stats.theilslopes(seriesC, seriesT, 0.95)
                 tau, p_value = scipy.stats.kendalltau(seriesT, seriesC)
                 spear, spear_p = scipy.stats.kendalltau(seriesT, seriesC)
                 corrCoeffP.append(spear)

            except(ValueError):
                    print 'Weird Error, continue'
            #
            # multiple correlation
            rx1=scipy.stats.spearmanr(interpTpar, np.array(allPotential))[0]
            rx2=scipy.stats.spearmanr(interpP, np.array(allPotential))[0]
            rx3=scipy.stats.spearmanr(interpTpar,interpP)[0]
            R=np.sqrt((rx1**2 + rx2**2 - 2*rx1*rx2*rx3)/(1-rx3**2))
            corrCoeffMultiple.append(R)


            #                     
            fig2=plt.figure()
            ax2=fig2.add_subplot(111)
            ax2.plot(seriesA, seriesT, '.', color='b')
            ax3=ax2.twinx()
            ax3.plot(np.array(seriesC)*1.0e-9, seriesT, '.', color='gold')
            ax3.set_ylabel('Pressure [nPa]')
        #    xx=np.arange(-10,10)
        #    ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
            plt.draw()
            font = {'family' : 'normal',
                              'weight' : 'bold',
                              'size'   : 22}
            plt.rc('font', **font)
            plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
            ax2.set_ylabel('T par [K]', fontweight='bold')
            ax2.set_xlabel('$\phi$',fontweight='bold')
           # ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(10, 12), xycoords='data', color='k')
           
            os.chdir('CorrelationPlots')
            plt.savefig('potential_tpar_pressure_scatter.pdf')
            plt.close()
            os.chdir('..')

            # MAKE SURE THIS WORKS WHEN I GET BACK!!!
            # CHECK SYM H 


        iCount+=1        
    except(IndexError):
        # nothing happens
        print 'no events on ' + date
  DT=DT+datetime.timedelta(days=1)


# pickle the arrays
import pickle

pickle.dump(corrCoeffTpar, open( "TparCorr.p", "wb" ) )
pickle.dump(corrCoeffTperp, open( "TperpCorr.p", "wb" ) )
pickle.dump(corrCoeffP, open( "PCorr.p", "wb" ) )
pickle.dump(corrCoeffMultiple,open("PTparCorr.p", "wb"))
    
print "median Tpar R: " + str(np.nanmedian(corrCoeffTpar))
print "median Tperp R: " + str(np.nanmedian(corrCoeffTperp))
print "median Pressure R: " + str(np.nanmedian(corrCoeffP))
print "median Multiple Corr R: " + str(np.nanmedian(corrCoeffMultiple))
                               
