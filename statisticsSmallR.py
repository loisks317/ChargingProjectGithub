# statisticsSmallR.py
#
# calculate the correlation coefficient of -5 V or less charging events
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
import matplotlib.dates as dt
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
dateStart='20130101'
dateEnd='20150401'
# data
hourGap=5
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
sat=['A','B']
ExcludeDates=['20130206', '20130227', '20130327', '20130417', '20130508']
typicalE=[  1.49845505e+01,   1.68136501e+01,   1.88537998e+01,
         2.11753502e+01,   2.37079487e+01,   2.65922985e+01,
         2.98283978e+01,   3.34866028e+01,   3.75669022e+01,
         4.21396484e+01,   4.72751961e+01,   5.29735489e+01,
         5.94457474e+01,   6.66917953e+01,   7.47820511e+01,
         8.38572006e+01,   9.40579453e+01,   1.05524994e+02,
         1.18328697e+02,   1.32680099e+02,   1.48860596e+02,
         1.66940536e+02,   1.87201355e+02,   2.09994751e+02,
         2.35531784e+02,   2.64164246e+02,   2.96243835e+02,
         3.32263031e+02,   3.72714294e+02,   4.18019684e+02,
         4.68812378e+02,   5.25795898e+02,   5.89744080e+02,
         6.61430664e+02,   7.41840698e+02,   8.32099792e+02,
         9.33263123e+02,   1.04666724e+03,   1.17393042e+03,
         1.31667065e+03,   1.47678711e+03,   1.65632043e+03,
         1.85766199e+03,   2.08355591e+03,   2.33688623e+03,
         2.62095947e+03,   2.93964502e+03,   3.29702295e+03,
         3.69787744e+03,   4.14748389e+03,   4.65175293e+03,
         5.21729639e+03,   5.85164258e+03,   6.56309180e+03,
         7.36107178e+03,   8.25599512e+03,   9.25981836e+03,
         1.03856299e+04,   1.16482715e+04,   1.30644863e+04,
         1.46528506e+04,   1.64343926e+04,   1.84324746e+04,
         2.06735430e+04,   2.31870078e+04,   2.60061426e+04,
         2.91679551e+04,   3.27142266e+04,   3.66916719e+04,
         4.11526406e+04,   4.61560000e+04,   5.17676797e+04]
corrCoeff=[[] for i in range(72)]
dateArr=[]
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
    density=(list(pyf['density'][...]))
    potential=-1.0*pyf['Vavg'][...]
    epochE=dt.date2num(pyf['epoch'][...])
    MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
    LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
    MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
    #
    # now get HOPE data
    #
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+sat[iSat])
    f=glob.glob('*'+date+'*')
    try: 
        pyf=pycdf.CDF(f[0])
    except:
        continue
    
    electronFlux=pyf['FEDO'][...]
    electronEnergy=pyf['HOPE_ENERGY_Ele'][...]
    electronFlux=electronFlux*electronEnergy
    epoch=pyf['Epoch_Ele'][...]
    L=pyf['L_Ele'][...]
    MLT=pyf['MLT_Ele'][...]
    os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
    mTime=dt.date2num(epoch)
    
    # define a charging event
    lowPotential=np.where(potential < -1)[0]
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
        interpFlux=[[] for i in range(72)]
        if len(potentialEvents[iCount]) < 50:
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
            for iEn in range(72):
                temp=np.swapaxes(electronFlux,1,0)[iEn][AdeleIsAwesome:AdeleIsAwesomer]
                interpFlux[iEn]+=list(temp)
                try:
                    interpFlux[iEn]=np.array(interpFlux[iEn][~np.isnan(interpFlux[iEn])])
                except:
                    interpFlux[iEn]=np.array(interpFlux[iEn])
                try:
                 seriesB = np.transpose(np.log10(interpFlux[iEn]))
                 seriesA = np.transpose(allPotential)
                 ts_res = scipy.stats.theilslopes(seriesB, seriesA, 0.95)
                 tau, p_value = scipy.stats.kendalltau(seriesA, seriesB)
                 spear, spear_p = scipy.stats.kendalltau(seriesA, seriesB)
                 corrCoeff[iEn].append(spear)

                except(ValueError):
                    print 'Weird Error, continue'
                print("For event count " + str(iCount) + " Kendall's tau for (Energy={0}) is {1} (p={2})".format(iEn,tau,p_value))
                print("Spearman's R for (Energy={0}) is {1} (p={2})".format(iEn,spear,spear_p))
                fig2=plt.figure()
                ax2=fig2.add_subplot(111)
                ax2.plot(seriesA, seriesB, '.', color='k')
                xx=np.arange(-200,0)
                ax2.plot(xx, ts_res[1] + ts_res[0] * xx,'-', linewidth=2, color='r')
                plt.draw()
                font = {'family' : 'normal',
                          'weight' : 'bold',
                          'size'   : 22}
                plt.rc('font', **font)
                plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
                ax2.set_ylabel('Log(Energy Flux)', fontweight='bold')
                ax2.set_xlabel('$\phi$',fontweight='bold')
                ax2.set_ylim(9,12)
                ax2.set_xlim(np.min(seriesA),0)
                ax2.annotate('r='+str(np.round(spear*100)/100.0), xy=(-10, 12), xycoords='data', color='k')
                os.chdir('CorrelationPlots_SpecificEvent')
                plt.savefig("date="+str(date)+"_sat="+str(sat[iSat])+"_event="+str(iCount)+"_energy="+str(iEn)+'.pdf')
                plt.close()
                os.chdir('..')
        iCount+=1        
    except(IndexError):
        # nothing happens
        print 'no events on ' + date
  DT=DT+datetime.timedelta(days=1)
 else:
  DT=DT+datetime.timedelta(days=1)
      

# now make a plot of the coff coeff var
stdCorr=[[] for i in range(72)]
medCorr=[[] for i in range(72)]
for iStd in range(72):
    stdCorr[iStd]=np.std(np.array(corrCoeff[iStd])[~np.isnan(np.array(corrCoeff[iStd]))])
    medCorr[iStd]=np.nanmedian(corrCoeff[iStd])

# pickle the arrays
import pickle

pickle.dump(dateArr, open( "dateArr25V.p", "wb" ) )
pickle.dump(corrCoeff, open( "corrCoeff25V.p", "wb" ) )
    
fig=plt.figure()
ax=fig.add_subplot(111)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
ax.errorbar(typicalE, medCorr, yerr=stdCorr, fmt='o')
ax.set_xscale('log')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('R')
font = {'family' : 'normal',
          'weight' : 'bold',
          'size'   : 22}
plt.rc('font', **font)
plt.savefig('dateStart='+dateStart+'_dateEnd='+dateEnd+'_corrCoeff25V.pdf')
plt.close()
