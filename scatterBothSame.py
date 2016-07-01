# scatterBothDiff.py
#
# get the difference between EFW and HOPE Charging Line
# also include as a function of percentage, October 2015
#
# LKS, September 2015
#
# import list 
import os
import glob 
import numpy as np
from matplotlib import pyplot as plt
import datetime
import matplotlib.dates as dt
import pickle
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator, MonthLocator
import pandas as pd
from spacepy import pycdf
from scipy import stats
#
#
nLbins=11
nmlt_bins=48
#
#
def depickle(name):
    with open(name, 'rb') as f:
        u = pickle._Unpickler(f)
        u.encoding = 'latin1'
        outf = u.load()
        return outf
os.chdir('HOPE_L3_A')
f=glob.glob('rbspa_rel02_ect-hope-PA-L3_20130524_v5.0.1.cdf')
pyf=pycdf.CDF(f[0])
energyChannels=pyf['HOPE_ENERGY_Ion'][...][0] # typical energy channel
delta=pyf['ENERGY_Ion_DELTA'][...][0]
os.chdir('..')

os.chdir('ScatterBoth')
    
IonTime=depickle('IonTime.p')
IonLine=np.array(depickle('IonLine.p'))
EFWSC=np.array(depickle('EFWSC.p'))
EFWTimes=depickle('EFWTimes.p')
EFluxesE=depickle('EFluxesE.p')
EFluxesIon=depickle('EfluxesIon.p')

IT=np.zeros(len(IonTime))
ET=np.zeros(len(EFWTimes))
for iTime in range(len(IT)):
    IT[iTime]=dt.date2num(datetime.datetime.strptime(str(pd.to_datetime(IonTime)[iTime]), '%Y-%m-%d %H:%M:%S'))
for iTime in range(len(ET)):
    ET[iTime]=dt.date2num(datetime.datetime.strptime(str(pd.to_datetime(EFWTimes)[iTime]), '%Y-%m-%d %H:%M:%S'))
#
#
print("begining other sort")
ionMatches=[]
EFWMatches=[]
#
# find datetime matches
for iT in range(len(IT)):
    for eT in range(len(ET)):
        if IT[iT]==ET[eT]:
            ionMatches.append(iT)
            EFWMatches.append(eT)
    
os.chdir('..')
subdir_name='ScatterCompare'
if not os.path.exists(subdir_name):
   os.umask(0) # unmask if necessary
   os.makedirs(subdir_name) 
os.chdir(subdir_name)#

IL=IonLine[ionMatches]
ErrorBars=np.zeros(len(IL))
for iEn in range(1,len(energyChannels)-1):
    indices=np.where((IL > energyChannels[iEn-1]) & (IL <= energyChannels[iEn+1]))[0]
    ErrorBars[indices]=delta[iEn]

fig=plt.figure()
ax=fig.add_subplot(111)
#ax.set_xlim(-100,100)
plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)
ax.scatter(EFWSC[EFWMatches], IL,s=60, c='LightSeaGreen', alpha=0.35, edgecolors='none')
plt.errorbar(EFWSC[EFWMatches],IL, yerr=ErrorBars, linestyle='none', color='k', alpha=0.1)

ax.plot(np.linspace(0,250,250), ls ='--', c='DarkViolet', lw=3)
#ax.plot(range(-100,100),np.zeros(200)+10, '-.', lw=2, c='k')
#ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylim(0,450)
ax.set_xlim(0,250)
ax.set_xlabel('EFW $\phi$ [-V]') 
ax.set_ylabel('Ion Line $\phi$ [-V]')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
plt.savefig('DiffScatter.png')
plt.close()



fig=plt.figure()
ax=fig.add_subplot(111)
#ax.set_xlim(-100,100)
plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)
ax.scatter(EFWSC[EFWMatches], IL,s=60, c='LightSeaGreen', alpha=0.35, edgecolors='none')
plt.errorbar(EFWSC[EFWMatches],IL, yerr=ErrorBars, linestyle='none', color='k', alpha=0.1)
# screen so only 1 std from line\
#

limD=np.zeros(len(EFWSC[EFWMatches]))
limY=np.zeros(len(EFWSC[EFWMatches]))
for iEN in range(len(EFWSC[EFWMatches])):
    if ((IL[iEN] > EFWSC[EFWMatches][iEN]+EFWSC[EFWMatches][iEN]*.6) or (IL[iEN]< EFWSC[EFWMatches][iEN]*.4)):
        limD[iEN]=np.nan
        limY[iEN]=np.nan
    else:
        limD[iEN]=IL[iEN]
        limY[iEN]=EFWSC[EFWMatches][iEN]
limD=limD[~np.isnan(limD)]
limY=limY[~np.isnan(limY)]

from scipy.optimize import curve_fit
newX = np.logspace(0, 3, base=10) 
def myExpFunc(x, a, b):
    return a * np.power(x, b)
popt, pcov = curve_fit(myExpFunc,limD, limY)
plt.plot(newX, myExpFunc(newX, *popt), ls='--', color='DarkViolet', lw=3, 
         label="({0:.3f}*x**{1:.3f})".format(*popt))
print("Exponential Fit: y = (a*(x**b))")
print("\ta = popt[0] = {0}\n\tb = popt[1] = {1}".format(*popt))




#
#slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(limY), np.log10(limD))
#plt.plot(np.logspace(0,3,100), (slope)*(np.logspace(0,3,100))+(10**intercept), ls='--',lw=3, c='DarkViolet')



#print(r_value)
#print(slope)

#ax.plot(np.linspace(0,1000,1000), ls ='--', c='DarkViolet', lw=3)
#ax.plot(range(-100,100),np.zeros(200)+10, '-.', lw=2, c='k')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e0, 1e3)
ax.set_xlim(1e0,1e3)
ax.set_xlabel('EFW $\phi$ [-V]') 
ax.set_ylabel('Ion Line $\phi$ [-V]')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
plt.savefig('DiffScatterlog.png')
plt.close()
#


fig=plt.figure()
ax=fig.add_subplot(111)
#ax.set_xlim(-100,100)
plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)
ax.scatter(EFWSC[EFWMatches], IL,s=60, c='LightSeaGreen', alpha=0.35, edgecolors='none')
#plt.errorbar(EFWSC[EFWMatches],IL, yerr=ErrorBars, linestyle='none', color='k', alpha=0.1)

ax.plot(np.linspace(0,1000,1000), ls ='--', c='DarkViolet', lw=3)
#ax.plot(range(-100,100),np.zeros(200)+10, '-.', lw=2, c='k')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e0, 1e3)
ax.set_xlim(1e0,1e3)
ax.set_xlabel('EFW $\phi$ [-V]') 
ax.set_ylabel('Ion Line $\phi$ [-V]')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
plt.savefig('DiffScatterlognoErr.png')
plt.close()

#
#
os.chdir('..')



# energy
#HIE=[ 1.49845505e+01,   1.68136501e+01,   1.88537998e+01,
#         2.11753502e+01,   2.37079487e+01,   2.65922985e+01,
#         2.98283978e+01,   3.34866028e+01,   3.75669022e+01,
#         4.21396484e+01,   4.72751961e+01,   5.29735489e+01,
#         5.94457474e+01,   6.66917953e+01,   7.47820511e+01,
#         8.38572006e+01,   9.40579453e+01,   1.05524994e+02,
#         1.18328697e+02,   1.32680099e+02,   1.48860596e+02,
#         1.66940536e+02,   1.87201355e+02,   2.09994751e+02,
#         2.35531784e+02,   2.64164246e+02,   2.96243835e+02,
#         3.32263031e+02,   3.72714294e+02,   4.18019684e+02,
#         4.68812378e+02,   5.25795898e+02,   5.89744080e+02,
#         6.61430664e+02,   7.41840698e+02,   8.32099792e+02,
#         9.33263123e+02,   1.04666724e+03,   1.17393042e+03,
#         1.31667065e+03,   1.47678711e+03,   1.65632043e+03,
#         1.85766199e+03,   2.08355591e+03,   2.33688623e+03,
#         2.62095947e+03,   2.93964502e+03,   3.29702295e+03,
#         3.69787744e+03,   4.14748389e+03,   4.65175293e+03,
#         5.21729639e+03,   5.85164258e+03,   6.56309180e+03,
#         7.36107178e+03,   8.25599512e+03,   9.25981836e+03,
#         1.03856299e+04,   1.16482715e+04,   1.30644863e+04,
#         1.46528506e+04,   1.64343926e+04,   1.84324746e+04,
#         2.06735430e+04,   2.31870078e+04,   2.60061426e+04,
#         2.91679551e+04,   3.27142266e+04,   3.66916719e+04,
#         4.11526406e+04,   4.61560000e+04,   5.17676797e+04]
#
#
#
##
##
#IT=np.zeros(len(IonTime))
#ET=np.zeros(len(EFWTimes))
#for iTime in range(len(IT)):
#    IT[iTime]=dt.date2num(datetime.datetime.strptime(str(pd.to_datetime(IonTime)[iTime]), '%Y-%m-%d %H:%M:%S'))
#for iTime in range(len(ET)):
#    ET[iTime]=dt.date2num(datetime.datetime.strptime(str(pd.to_datetime(EFWTimes)[iTime]), '%Y-%m-%d %H:%M:%S'))
##
## make a general line plot
#fig=plt.figure()
#months = MonthLocator(interval=5) 
#months2=MonthLocator(interval=3)
#daysFmt = DateFormatter('%Y- \n %m')
#
#fig.gca().xaxis.set_major_locator(months)
#fig.gca().xaxis.set_major_formatter(daysFmt)
#fig.gca().xaxis.set_minor_locator(months2)
#ax=fig.add_subplot(111)
#plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.25)
##ax.scatter(np.array(IT), np.array(IonLine), c='b')
#ET=np.array(ET)
#sET=[x for (y,x) in sorted(zip(ET,EFWSC))]
#sortedET=sorted(ET)
#ax.plot(sortedET,-1*np.array(sET),c='blueviolet', lw=2)
#
##ax.set_yscale('log')
#ax.set_ylim(-205, 5)
#ax.set_ylabel('Spacecraft Potential [V]')
#ax.set_xlabel('Date')
##plt.legend(['Ion Line', 'EFW'])
#
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#plt.rc('font', **font)
#plt.savefig('EFWmission.png')
#plt.close()
#
##
## now separate by when they both align
#matchesIon=[i for i, item in enumerate(IT) if item in set(ET)]
#matchesE=[i for i, item in enumerate(ET) if item in set(IT)]
#
#
#
##
### make a general line plot of matches
#fig=plt.figure()
#months = MonthLocator(interval=6) 
#months2=MonthLocator(interval=3)
#daysFmt = DateFormatter('%Y- \n %m')
#
#fig.gca().xaxis.set_major_locator(months)
#fig.gca().xaxis.set_major_formatter(daysFmt)
#fig.gca().xaxis.set_minor_locator(months2)
#ax=fig.add_subplot(111)
#plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.2)
##ax.scatter(np.array(IT)[matchesIon], np.array(IonLine)[matchesIon], c='b')
#ax.plot(np.array(ET),np.array(EFWSC),c='purple', lw=2)
#ax.set_yscale('log')
#ax.set_ylim(1e0,1e4)
#ax.set_ylabel('Spacecraft Potential [V]')
#ax.set_xlabel('Time')
#plt.legend(['Ion Line', 'EFW'], loc='upper right',prop={'size':16})
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#plt.rc('font', **font)
#plt.savefig('comparisonMatches.png')
#plt.close()
##
## get the average differences and plot
#tempI=list(IT[matchesIon])
#tempE=list(ET[matchesE])
#IonTimeM=[i for i, x in enumerate(tempI) if tempI.count(x) > 1]
#ETimeM=[i for i, x in enumerate(tempE) if tempE.count(x) > 1]
## remove every other index
#listRMI=IonTimeM[1::2]
#listRME=ETimeM[1::2]
#
#IonLineR=[i for j, i in enumerate(np.array(IonLine)[matchesIon]) if j not in listRMI]
#EFWSCR=[i for j, i in enumerate(np.array(EFWSC)[matchesE]) if j not in listRME]
#IT=tempI
#ET=tempE
#
#
#diff=np.array(IonLineR) - np.array(EFWSCR)
#IonLineR=np.array(IonLineR)
#EFWSCR=np.array(EFWSCR)
#percent=100*(IonLineR - EFWSCR)/((np.abs(IonLineR + EFWSCR))/2.0)
#fig=plt.figure()
#ax=fig.add_subplot(111)
#binwidth=50
#plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)
#plt.hist(diff[diff<0],bins=range(-200, 1000 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DarkViolet'])
#plt.hist(diff[diff>0],bins=range(-200, 1000 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DeepPink'])
#ax.plot(range(-200,1000),np.zeros(1200)+10, '-.', lw=2, c='k')
#ax.set_yscale('log')
#ax.set_xlabel('Ion $\phi$ - EFW $\phi$ [V]')
#ax.set_ylabel('Occurrences')
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#plt.rc('font', **font)
#plt.savefig('DiffcomparisonMatches.png')
#plt.close()
#
#
#### PERCENTAGE
#print "median percentage difference is: " + str(np.nanmedian(percent))
#fig=plt.figure()
#ax=fig.add_subplot(111)
#binwidth=25
#plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)
#plt.hist(percent[percent<0],bins=range(-250, 250 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DarkViolet'])
#plt.hist(percent[percent>0],bins=range(-250, 250 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DeepPink'])
#ax.plot(range(-250,250),np.zeros(500), '-.', lw=2, c='k')
#ax.set_yscale('log')
#ax.set_xlabel('Percentage Difference ($\%$)')
#ax.set_ylabel('Occurrences')
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 20}
#plt.rc('font', **font)
#plt.savefig('PercentageMatches.png')
#plt.close()
#
#
#
##
## ZOOMED IN
## 
#diff=np.array(IonLineR) - np.array(EFWSCR)
#fig=plt.figure()
#ax=fig.add_subplot(111)
#binwidth=5
#ax.set_xlim(-100,100)
#plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)
#plt.hist(diff[diff<0],bins=range(-100, 100 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DarkViolet'])
#plt.hist(diff[diff>0],bins=range(-100, 100 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DeepPink'])
#ax.plot(range(-100,100),np.zeros(200)+10, '-.', lw=2, c='k')
#ax.set_yscale('log')
#ax.set_xlabel('Ion $\phi$ - EFW $\phi$ [V]')
#ax.set_ylabel('Occurrences')
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#plt.rc('font', **font)
#plt.savefig('DiffcomparisonMatchesZoomedIn.png')
#plt.close()
#
#
##
## ZOOMED IN PERCENTAGES
## 
#fig=plt.figure()
#ax=fig.add_subplot(111)
#binwidth=5
#ax.set_xlim(-100,100)
#plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.2)
#plt.hist(percent[percent<0],bins=range(-100, 100 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DarkViolet'])
#plt.hist(percent[percent>0],bins=range(-100, 100 + binwidth, binwidth), histtype='bar', rwidth=1, color=['DeepPink'])
#ax.plot(range(-100,100),np.zeros(200)+10, '-.', lw=2, c='k')
#ax.set_yscale('log')
#ax.set_xlabel('Perecent Difference [$\%$]')
#ax.set_ylabel('Occurrences')
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#plt.rc('font', **font)
#plt.savefig('PercentagecomparisonMatchesZoomedIn.png')
#plt.close()
#
#
#
##
## scatter plot of matching times
#for iEn in range(72):
#    fig2=plt.figure()
#    ax2=fig2.add_subplot(111)
#    ax2.plot(np.array(EFluxesIon[iEn])[matchesIon],np.array(IonLine)[matchesIon], '.', color='k')
#    ax2.set_yscale('log')
#    ax2.set_xscale('log')
#    ax2.set_xlim(1e8,1e12)
#    ax2.set_ylim(1e0,1e4)
#    plt.draw()
#    font = {'family' : 'normal',
#                      'weight' : 'bold',
#                      'size'   : 22}
#    plt.rc('font', **font)
#    plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
#    ax2.set_ylabel('Spacecraft Potential [V]', fontweight='bold')
#    ax2.set_xlabel('Energy Flux',fontweight='bold')
#    #subdir_name='ScatterBoth'
#    #if not os.path.exists(subdir_name):
#    #           os.umask(0) # unmask if necessary
#    #           os.makedirs(subdir_name, 0777) 
#    #os.chdir(subdir_name)
#    plt.savefig('Matches_IonLine_Energy='+str(HIE[iEn])+'_eV.png')
#    plt.close()
#
#    #
#    # plot EFW stats
#    fig2=plt.figure()
#    ax2=fig2.add_subplot(111)
#    ax2.plot(np.array(EFluxesE[iEn])[matchesE],np.array(EFWSC)[matchesE], '.', color='k')
#    ax2.set_xlim(1e8,1e12)
#    ax2.set_ylim(1e0,1e4)
#    ax2.set_yscale('log')
#    ax2.set_xscale('log')
#    plt.draw()
#    font = {'family' : 'normal',
#                      'weight' : 'bold',
#                      'size'   : 22}
#    plt.rc('font', **font)
#    plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
#    ax2.set_ylabel('Spacecraft Potential [V]', fontweight='bold')
#    ax2.set_xlabel('Energy Flux',fontweight='bold')
#    #subdir_name='ScatterBoth'
#    #if not os.path.exists(subdir_name):
#    #           os.umask(0) # unmask if necessary
#    #           os.makedirs(subdir_name, 0777) 
#    #os.chdir(subdir_name)
#    plt.savefig('Matches_EFWSC_Energy='+str(HIE[iEn])+'_eV.png')
#    plt.close()


# get the 
