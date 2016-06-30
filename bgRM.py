# bgRM.py
#
# spectograms with removed background fluxes
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
sat=['A']
hourGap=5
uniqueDates=['20130208','20130317']
EnBins=72
Lbins=56
LbinsArr=np.linspace(1, 6.5,Lbins)
#
for isat in range(len(sat)):
    dataEFW={}
    dataHOPE={}
    # open the files
#   
#   # now load in HOPE data for this
    for iDate in range(len(uniqueDates)):
    # try:
        date=uniqueDates[iDate]
        print date
        # first get EFW
        os.chdir('EFW_L3_'+sat[isat])
        f=glob.glob('*'+uniqueDates[iDate]+'*')
        pyf=pycdf.CDF(f[0])
        EFWdensity=(list(pyf['density'][...]))
        EFWpotential=pyf['Vavg'][...]
        epochE=dt.date2num(pyf['epoch'][...])
        # now get HOPE data
        os.chdir('..')
        os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/HOPE_L3_'+sat[isat])
        f=glob.glob('*'+uniqueDates[iDate]+'*')
        pyf=pycdf.CDF(f[0])
        
        
        counts=pyf['Counts_P_Omni'][...]
        HIE=pyf['HOPE_ENERGY_Ion'][...]
        IonFlux=pyf['FPDO'][...]
        epochI=pyf['Epoch_Ion'][...]
        epoch=epochI
        Lion=pyf['L_Ion'][...]
        MLT=pyf['MLT_Ion'][...]

        os.chdir('..')
        mTime=dt.date2num(epochI)
        #
        # now get the background data
        os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject')
        bgCounts=pickle.load(open('bgCounts.p','rb'))
        for iTime in range(len(epochI)):
            for iEN in range(72):
                try:
                    idx = np.argmin(np.abs(LbinsArr - np.array(Lion)[iTime]))
                    counts[iTime][iEN]=counts[iTime][iEN]-bgCounts[iEN][idx]
                except(ValueError):
                    counts[iTime][iEN]=np.nan
        #
        # extract an ion line
        chargingLine=[]
        chargingTimes=[]
        for iTime in range(len(epochI)):
            # find the gradient
          check=0
          if Lion[iTime] > 5:
            temp=np.where(counts[iTime]>0)[0][0:5] # should be clear
            for iTemp in range(len(temp)):
                # is it above like 6 eV and less than 3100 eV
                if (temp[iTemp]>10) and (temp[iTemp]<45):
                    # is it above 50 counts at hte ion line
                    try:
                     if (counts[iTime][temp[iTemp]]-counts[iTime][temp[iTemp+2]])>10:
                      # is it only a few thick under these same requirements
                     #  if temp[iTemp+3] > temp[iTemp]+3:
                        # this is a charging line
                        chargingLine.append(HIE[iTime][temp[iTemp]])
                        chargingTimes.append(epochI[iTime])
                        check=1
                    except(IndexError):
                        break
          if check!=1:
              chargingLine.append(0)
              chargingTimes.append(epochI[iTime])
                            
 
         
        chargingLine=runningMedian(chargingLine, 10)
        chargingTimes=runningMedian(dt.date2num(chargingTimes),10)
        # now we have both
        fig=plt.figure()
        cmap=plt.cm.jet
        cbmin=-50
        cbmax=500
        ax1=fig.add_subplot(111)
        plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
        ax1.set_ylabel('Energy (eV)', fontsize=25, fontweight='bold')
        ax1.set_xlabel('Time', fontsize=25, fontweight='bold')
        ax1.set_yscale('log')
        mEnergies=np.nanmedian(HIE, axis=0)
        X,Y=np.meshgrid(mTime, mEnergies)
        mdata=ma.masked_invalid(counts.transpose())
        col=ax1.pcolormesh(X,Y,mdata,cmap='jet', vmax=cbmax, vmin=cbmin)
        font = {'family' : 'normal',
                 'weight' : 'bold',
                 'size'   : 22}
        plt.rc('font', **font)
        cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
        cb=plt.colorbar(col,cax=cbaxes)
        #cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(cbmin,cbmax))
        cb.set_label('# of Counts', fontsize=25, fontweight='bold')
        ax1.tick_params(axis='both', which='major', labelsize=22)
        cb.ax.tick_params(labelsize=30)
        hfmt = dt.DateFormatter('%H:%M')
        ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
        ax1.xaxis.set_major_formatter(hfmt)
        ax1.set_ylim(1, 45000)
        
        
        # also overplot S/C Potential
        ax2=ax1.twinx()
        ax2.plot(chargingTimes,-1*np.array(chargingLine) , lw=6, c='Silver')
        ax2.plot(chargingTimes,-1*np.array(chargingLine), lw=2, c='Black')
        #ax2.scatter(chargingTimes,-1*np.array(chargingLine), s=40, c='Silver')
        #ax2.scatter(chargingTimes,-1*np.array(chargingLine), s=35, c='Black')
        ax2.plot(epochE,-1* EFWpotential, lw=3, c='red')
        ax2.set_ylabel('S/C Potential', fontsize=25, fontweight='bold')
        ax2.set_ylim(-1000,0)
        
        # Add L LABEL
        ax3 = fig.add_axes(ax1.get_position())
        ax3.patch.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.set_yscale('log')
        for spinename, spine in ax3.spines.iteritems():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax3.spines['bottom'].set_position(('outward', 65))
        L_conc=[round(x, 2) for x in Lion]
        p3=ax3.plot(mTime,np.swapaxes(counts,1,0)[0]  , 'g', alpha=0)
        # I hate doing this but
        spacing = len(L_conc)/7
        spacedArray=L_conc[::spacing]
        ax3.set_xticks(epoch[::spacing])
        nans=np.where(np.array(L_conc[::spacing]) < 0)
        if len(nans[0])>0:
            try:
                newSpot = L_conc[spacing*7+5]
                print newSpot
                if newSpot < 0:
                    newSpot=L_conc[spacing*7-5]
                    print newSpot
                spacedArray[nans[0][0]]=newSpot
            except(IndexError):
                newSpot=L_conc[spacing*7-5]
                print newSpot
                spacedArray[nans[0][0]]=newSpot
        ax3.set_xticklabels(spacedArray)
        ax3.set_xlabel('L-Shell', fontsize=22, fontweight='bold')
        #
        # Add MLT Label
        ax5 = fig.add_axes(ax1.get_position())
        ax5.patch.set_visible(False)
        ax5.yaxis.set_visible(False)
        ax5.set_yscale('log')
        for spinename, spine in ax3.spines.iteritems():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax5.spines['bottom'].set_position(('outward', 125))
        MLT_conc=[round(x, 1) for x in MLT]
        p4=ax5.plot(mTime,  np.swapaxes(counts,1,0)[0], 'purple', alpha=0)
        spacing = len(MLT_conc)/8
        ax5.set_xticks(mTime[::spacing])
        ax5.set_xticklabels(MLT_conc[::spacing])
        ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
        
        
        subdir_name='/Users/loisks/Documents/ResearchProjects/ChargingProject/PosChargingPlots'
        if not os.path.exists(subdir_name):
             os.umask(0) # unmask if necessary
             os.makedirs(subdir_name, 0777) 
        os.chdir(subdir_name)#
        fig.set_size_inches(13,9)
        plt.savefig('counts_'+date+'_ion_'+str(sat[isat])+'_spectrogram.png')
        plt.close(fig)
        os.chdir('..')
   #  except(OSError):
   #     print 'bad date: ' + date  
#
        # PLOT DIFF # FLUX

        # now we have both
        fig=plt.figure()
        cmap=plt.cm.jet
        cbmin=4
        cbmax=8
        ax1=fig.add_subplot(111)
        plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
        ax1.set_ylabel('Energy (eV)', fontsize=25, fontweight='bold')
        ax1.set_xlabel('Time', fontsize=25, fontweight='bold')
        ax1.set_yscale('log')
        mEnergies=np.nanmedian(HIE, axis=0)
        X,Y=np.meshgrid(mTime, mEnergies)
        mdata=ma.masked_invalid(np.log10(IonFlux).transpose())
        col=ax1.pcolormesh(X,Y,mdata,cmap='jet', vmax=cbmax, vmin=cbmin)
        font = {'family' : 'normal',
                 'weight' : 'bold',
                 'size'   : 22}
        plt.rc('font', **font)
        cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
        cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(cbmin,cbmax+1))
        cb.set_label('Differential # Flux', fontsize=25, fontweight='bold')
        ax1.tick_params(axis='both', which='major', labelsize=22)
        cb.ax.tick_params(labelsize=30)
        hfmt = dt.DateFormatter('%H:%M')
        ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
        ax1.xaxis.set_major_formatter(hfmt)
        ax1.set_ylim(1, 45000)
        
        
        # also overplot S/C Potential
        ax2=ax1.twinx()
        ax2.plot(chargingTimes,-1*np.array(chargingLine) , lw=6, c='Silver')
        ax2.plot(chargingTimes,-1*np.array(chargingLine), lw=2, c='Black')
        #ax2.scatter(chargingTimes,-1*np.array(chargingLine), s=40, c='Silver')
        #ax2.scatter(chargingTimes,-1*np.array(chargingLine), s=35, c='Black')
        ax2.plot(epochE,-1* EFWpotential, lw=3, c='red')
        ax2.set_ylabel('S/C Potential', fontsize=25, fontweight='bold')
        ax2.set_ylim(-400,0)
        
        # Add L LABEL
        ax3 = fig.add_axes(ax1.get_position())
        ax3.patch.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.set_yscale('log')
        for spinename, spine in ax3.spines.iteritems():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax3.spines['bottom'].set_position(('outward', 65))
        L_conc=[round(x, 2) for x in Lion]
        p3=ax3.plot(mTime,np.swapaxes(counts,1,0)[0]  , 'g', alpha=0)
        # I hate doing this but
        spacing = len(L_conc)/7
        spacedArray=L_conc[::spacing]
        ax3.set_xticks(epoch[::spacing])
        nans=np.where(np.array(L_conc[::spacing]) < 0)
        if len(nans[0])>0:
            try:
                newSpot = L_conc[spacing*7+5]
                print newSpot
                if newSpot < 0:
                    newSpot=L_conc[spacing*7-5]
                    print newSpot
                spacedArray[nans[0][0]]=newSpot
            except(IndexError):
                newSpot=L_conc[spacing*7-5]
                print newSpot
                spacedArray[nans[0][0]]=newSpot
        ax3.set_xticklabels(spacedArray)
        ax3.set_xlabel('L-Shell', fontsize=22, fontweight='bold')
        #
        # Add MLT Label
        ax5 = fig.add_axes(ax1.get_position())
        ax5.patch.set_visible(False)
        ax5.yaxis.set_visible(False)
        ax5.set_yscale('log')
        for spinename, spine in ax3.spines.iteritems():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax5.spines['bottom'].set_position(('outward', 125))
        MLT_conc=[round(x, 1) for x in MLT]
        p4=ax5.plot(mTime,  np.swapaxes(counts,1,0)[0], 'purple', alpha=0)
        spacing = len(MLT_conc)/8
        ax5.set_xticks(mTime[::spacing])
        ax5.set_xticklabels(MLT_conc[::spacing])
        ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
        
        
        subdir_name='/Users/loisks/Documents/ResearchProjects/ChargingProject/PosChargingPlots'
        if not os.path.exists(subdir_name):
             os.umask(0) # unmask if necessary
             os.makedirs(subdir_name, 0777) 
        os.chdir(subdir_name)#
        fig.set_size_inches(13,9)
        plt.savefig('flux_'+date+'_ion_'+str(sat[isat])+'_spectrogram.png')
        plt.close(fig)
        os.chdir('..')
