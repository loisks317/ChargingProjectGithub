# readFiles.py
#
# read the EFW and HOPE combined files
# interpolate the times
# identify peak times of spacecraft charging
# associate with electron fluxes
#
# LKS, July 2015
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

def RunningMedian(seq, M):
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
hourGap=1
hour1=12
hour2=16
Lbins=56
LbinsArr=np.linspace(1, 6.5,Lbins)
uniqueDates=['20130208']
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
        print(date)
        datetime1=datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]), int(hour1), 0, 0)
        datetime2=datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]), int(hour2), 0, 0)
        # first get EFW
        os.chdir('EFW_L3_'+sat[isat])
        f=glob.glob('*'+uniqueDates[iDate]+'*')
        pyf=pycdf.CDF(f[0])
        EFWdensity=(list(pyf['density'][...]))
        EFWpotential=pyf['Vavg'][...]
        epochE=dt.date2num(pyf['epoch'][...])
        #MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
        #LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
        #MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
        #
        # now get HOPE data
        os.chdir('..')
        os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/HOPE_L3_'+sat[isat])
        f=glob.glob('*'+uniqueDates[iDate]+'*')
        pyf=pycdf.CDF(f[0])
        
        electronFlux=pyf['FEDO'][...]
        counts=pyf['Counts_P_Omni'][...]
        HIE=np.average(pyf['HOPE_ENERGY_Ion'][...],axis=0)
        electronEnergy=pyf['HOPE_ENERGY_Ele'][...]
        electronFlux=electronFlux*electronEnergy
        epoch=pyf['Epoch_Ele'][...]
        epochI=pyf['Epoch_Ion'][...]
        L=pyf['L_Ele'][...]
        Lion=pyf['L_Ion'][...]
        MLT=pyf['MLT_Ele'][...]
        #MLAT=pyf['MLAT_Ele'][...]
        os.chdir('..')
        mTime=dt.date2num(epoch)
        startEFW=np.where(epoch > datetime1)[0][0]
        endEFW=np.where(epoch > datetime2)[0][0]
        #
        # now get the background data
        os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/BackgroundCounts')
        # python 3.4 work around for compatibility with 2.7
        with open('bgCounts.p', 'rb') as f:
            u=pickle._Unpickler(f)
            u.encoding='latin1'
            bgCounts=u.load()
        #bgCounts=pickle.load(open('bgCounts.p','rb'))
        for iTime in range(len(epochI)):
            for iEN in range(72):
               idx = np.argmin(np.abs(LbinsArr - np.array(Lion)[iTime]))
               try:
                   counts[iTime][iEN]=counts[iTime][iEN]-bgCounts[iEN][idx]
               except(ValueError):
                   counts[iTime][iEN]=counts[iTime][iEN]
        #
        # extract an ion line
        os.chdir('..')
        #
        # GET THE ION LINE
        chargingLine=[]
        chargingTimes=[]
        gradient=np.zeros((len(epochI),55))
        for iTime in range(len(epochI)):
          if Lion[iTime] > 5: # must be near apogee
            check=0
            temp=np.where((np.array(counts)[iTime]>0))[0][0:5] # should be clear
            for iTemp in range(len(temp)):
                # is it above like 6 eV

                if ((temp[iTemp]>10) & temp[iTemp]<45):
                    # is it above 50 counts at hte ion line
                    try:
                       if (counts[iTime][temp[iTemp]]-counts[iTime][temp[iTemp+2]])>10:
                      # is it only a few thick under these same requirements
                     #  if temp[iTemp+3] > temp[iTemp]+3:
                        # this is a charging line
                         chargingLine.append(HIE[temp[iTemp]])
                         chargingTimes.append(epochI[iTime])
                         check=1
                   # except(IndexError):
                    except(IndexError):
                        break
                        print('here')
            if check!=1:
                    # nothing
                    continue
        


        #chargingLine=RunningMedian(chargingLine,40)
        #chargingTimes=RunningMedian(dt.date2num(chargingTimes),40)
        #stop
    
       
 #   # now we have both
 #   fig=plt.figure()
 #   cmap=plt.cm.jet
 #   cbmin=8
 #   cbmax=11
 #   ax1=fig.add_subplot(111)
 #   plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
 #   ax1.set_ylabel('Energy (eV)', fontsize=25, fontweight='bold')
 #   ax1.set_xlabel('Time', fontsize=25, fontweight='bold')
 #   ax1.set_yscale('log')
 #   mEnergies=np.mean(electronEnergy, axis=0)
 #   X,Y=np.meshgrid(mTime, mEnergies)
 #   mdata=ma.masked_invalid(np.log10(electronFlux).transpose())
 #   col=ax1.pcolormesh(X,Y,mdata,cmap='jet', vmax=cbmax, vmin=cbmin)
 #   font = {'family' : 'normal',
 #            'weight' : 'bold',
 #            'size'   : 22}
 #   plt.rc('font', **font)
 #   cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
 #   cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(cbmin,cbmax))
 #   cb.set_label('Energy Flux', fontsize=25, fontweight='bold')
 #   ax1.tick_params(axis='both', which='major', labelsize=22)
 #   cb.ax.tick_params(labelsize=30)
 #   hfmt = dt.DateFormatter('%H:%M')
 #   ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
 #   ax1.xaxis.set_major_formatter(hfmt)
 #   ax1.set_ylim(1e2, 50000)
 #
 #   
 #   # also overplot S/C Potential
 #   ax2=ax1.twinx()
 #   ax2.plot(chargingTimes,-1*np.array(chargingLine) , lw=6, c='Silver')
 #   ax2.plot(chargingTimes,-1*np.array(chargingLine), lw=2, c='Black')
 #   ax2.set_ylabel('S/C Potential', fontsize=25, fontweight='bold')
 #   ax2.set_ylim(-300,0)
 #
 #   # Add L LABEL
 #   ax3 = fig.add_axes(ax1.get_position())
 #   ax3.patch.set_visible(False)
 #   ax3.yaxis.set_visible(False)
 #   ax3.set_yscale('log')
 #   for spinename, spine in ax3.spines.iteritems():
 #       if spinename != 'bottom':
 #           spine.set_visible(False)
 #   ax3.spines['bottom'].set_position(('outward', 65))
 #   L_conc=[round(x, 2) for x in L]
 #   p3=ax3.plot(mTime,electronFlux  , 'g', alpha=0)
 #   # I hate doing this but
 #   spacing = len(L_conc)/7
 #   spacedArray=L_conc[::spacing]
 #   ax3.set_xticks(epoch[::spacing])
 #   nans=np.where(np.array(L_conc[::spacing]) < 0)
 #   if len(nans[0])>0:
 #       try:
 #           newSpot = L_conc[spacing*7+5]
 #           print newSpot
 #           if newSpot < 0:
 #               newSpot=L_conc[spacing*7-5]
 #               print newSpot
 #           spacedArray[nans[0][0]]=newSpot
 #       except(IndexError):
 #           newSpot=L_conc[spacing*7-5]
 #           print newSpot
 #           spacedArray[nans[0][0]]=newSpot
 #   ax3.set_xticklabels(spacedArray)
 #   ax3.set_xlabel('L-Shell', fontsize=22, fontweight='bold')
 #   #
 #   # Add MLT Label
 #   ax5 = fig.add_axes(ax1.get_position())
 #   ax5.patch.set_visible(False)
 #   ax5.yaxis.set_visible(False)
 #   ax5.set_yscale('log')
 #   for spinename, spine in ax3.spines.iteritems():
 #       if spinename != 'bottom':
 #           spine.set_visible(False)
 #   ax5.spines['bottom'].set_position(('outward', 125))
 #   MLT_conc=[round(x, 1) for x in MLT]
 #   p4=ax5.plot(mTime,  electronFlux, 'purple', alpha=0)
 #   spacing = len(MLT_conc)/8
 #   ax5.set_xticks(mTime[::spacing])
 #   ax5.set_xticklabels(MLT_conc[::spacing])
 #   ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
 #   
 #   # Add MLAT Label
 #  # ax6 = fig.add_axes(ax1.get_position())
 #  # ax6.patch.set_visible(False)
 #  # ax6.yaxis.set_visible(False)
 #  # ax6.set_yscale('log')
 #  # for spinename, spine in ax3.spines.iteritems():
 #  #     if spinename != 'bottom':
 #  #         spine.set_visible(False)
 #  # ax6.spines['bottom'].set_position(('outward', 185))
 #  # MLAT_conc=[round(x, 1) for x in MLAT]
 #  # p5=ax6.plot(mTime,  electronFlux, 'purple', alpha=0)
 #  # spacing = len(MLAT_conc)/8
 #  # ax6.set_xticks(mTime[::spacing])
 #  # ax6.set_xticklabels(MLT_conc[::spacing])
 #  # ax6.set_xlabel('MLAT', fontsize=22, fontweight='bold')
 #
 #   subdir_name='/Users/loisks/Documents/ResearchProjects/ChargingProject/PosChargingPlots'
 #   if not os.path.exists(subdir_name):
 #        os.umask(0) # unmask if necessary
 #        os.makedirs(subdir_name, 0777) 
 #   os.chdir(subdir_name)#
 #   fig.set_size_inches(13,9)
 #   plt.savefig('IonLine_'+date+'_'+str(sat[isat])+'_spectrogram.pdf')
 #   plt.close(fig)
 #   os.chdir('..')
 # except(OSError):
 #       print 'bad date: ' + date  
#
#
# generate ion spectrogram
# IONSSSSSS
       #try:
       # now remove non ordered chargingtimes
        for iT in range(len(chargingTimes)-2):
           if (chargingTimes[iT+1]-chargingTimes[iT]).total_seconds() > 3600:
               if (chargingTimes[iT+2]-chargingTimes[iT+1]).total_seconds() > 3600:
                   # this is a bad value
                   chargingLine[iT]=np.nan
        if (chargingTimes[-1]-chargingTimes[-2]).total_seconds() > 3600:
           chargingLine[-1]=np.nan

        date=uniqueDates[iDate]
        print(date)
       # first get EFW
       #os.chdir('EFW_L3_'+sat[isat])
       #f=glob.glob('*'+uniqueDates[iDate]+'*')
       #pyf=pycdf.CDF(f[0])
       #density=(list(pyf['density'][...]))
       #potential=pyf['Vavg'][...]
       #epochE=dt.date2num(pyf['epoch'][...])
       #MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
       #LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
       ##
       #MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
       #
       # now get HOPE data
        os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/HOPE_L3_'+sat[isat])
        f=glob.glob('*'+uniqueDates[iDate]+'*')
        pyf=pycdf.CDF(f[0])
        
        ionFlux=pyf['FPDO'][...]
        electronEnergy=pyf['HOPE_ENERGY_Ion'][...]
        #ionFlux=electronFlux*electronEnergy
        epoch=pyf['Epoch_Ion'][...]
        L=pyf['L_Ion'][...]
        MLT=pyf['MLT_Ion'][...]
        #MLAT=pyf['MLAT_Ele'][...]
        os.chdir('..')
        mTime=dt.date2num(epoch)
        startHope=np.where(epoch > datetime1)[0][0]
        endHope=np.where(epoch > datetime2)[0][0]
        #
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
        mEnergies=HIE
        X,Y=np.meshgrid(mTime, mEnergies)
        mdata=ma.masked_invalid(np.log10(ionFlux).transpose())
        col=ax1.pcolormesh(X,Y,mdata,cmap='jet', vmax=cbmax, vmin=cbmin)
        font = {'family' : 'normal',
                 'weight' : 'bold',
                 'size'   : 22}
        plt.rc('font', **font)
        cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
        cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(cbmin,cbmax+1))
        cb.set_label('# Flux', fontsize=25, fontweight='bold')
        ax1.tick_params(axis='both', which='major', labelsize=22)
        cb.ax.tick_params(labelsize=30)
        hfmt = dt.DateFormatter('%H:%M')
        ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
        ax1.xaxis.set_major_formatter(hfmt)
        ax1.set_ylim(3, 45000)
        ax1.set_xlim(epoch[startHope], epoch[endHope])
        
        # also overplot S/C Potential
        ax2=ax1.twinx()
        #ax2.plot(chargingTimes,-1*np.array(chargingLine) , lw=6, c='Silver')
        #ax2.plot(chargingTimes,-1*np.array(chargingLine), lw=2, c='Black')
        #ax2.scatter(chargingTimes,-1*np.array(chargingLine), s=40, c='Black')
        ax2.plot(epochE,-1* EFWpotential, lw=3, c='red')
        ax2.set_ylabel('S/C Potential [V]', fontsize=25, fontweight='bold')
        ax2.set_ylim(-1000,0)
        ax2.set_xlim(epoch[startHope], epoch[endHope])
        
        # Add L LABEL
        ax3 = fig.add_axes(ax1.get_position())
        ax3.patch.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.set_yscale('log')
        for spinename, spine in ax3.spines.items():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax3.spines['bottom'].set_position(('outward', 65))
        L_conc=[round(x, 2) for x in L]
        p3=ax3.plot(mTime,np.swapaxes(ionFlux,1,0)[0]  , 'g', alpha=0)
        # I hate doing this but
        spacing =int( len(L_conc)/7)
        spacedArray=L_conc[::spacing]
        ax3.set_xticks(epoch[::spacing])
        nans=np.where(np.array(L_conc[::spacing]) < 0)
        if len(nans[0])>0:
            try:
                newSpot = L_conc[spacing*7+5]
                print(newSpot)
                if newSpot < 0:
                    newSpot=L_conc[spacing*7-5]
                    print(newSpot)
                spacedArray[nans[0][0]]=newSpot
            except(IndexError):
                newSpot=L_conc[spacing*7-5]
                print(newSpot)
                spacedArray[nans[0][0]]=newSpot
        ax3.set_xticklabels(spacedArray)
        ax3.set_xlabel('L', fontsize=22, fontweight='bold')
        #
        # Add MLT Label
        ax5 = fig.add_axes(ax1.get_position())
        ax5.patch.set_visible(False)
        ax5.yaxis.set_visible(False)
        ax5.set_yscale('log')
        for spinename, spine in ax3.spines.items():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax5.spines['bottom'].set_position(('outward', 125))
        MLT_conc=[round(x, 1) for x in MLT]
        p4=ax5.plot(mTime,  np.swapaxes(ionFlux,1,0)[0], 'purple', alpha=0)
        spacing = int(len(MLT_conc)/8)
        ax5.set_xticks(mTime[::spacing])
        ax5.set_xticklabels(MLT_conc[::spacing])
        ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
        
        # Add MLAT Label
     #   ax6 = fig.add_axes(ax1.get_position())
     #   ax6.patch.set_visible(False)
     #   ax6.yaxis.set_visible(False)
     #   ax6.set_yscale('log')
     #   for spinename, spine in ax3.spines.iteritems():
     #       if spinename != 'bottom':
     #           spine.set_visible(False)
     #   ax6.spines['bottom'].set_position(('outward', 185))
     #   MLAT_conc=[round(x, 1) for x in MLAT]
     #   p5=ax6.plot(mTime,  electronFlux, 'purple', alpha=0)
     #   spacing = len(MLAT_conc)/8
     #   ax6.set_xticks(mTime[::spacing])
     #   ax6.set_xticklabels(MLT_conc[::spacing])
     #   ax6.set_xlabel('MLAT', fontsize=22, fontweight='bold')
        
        subdir_name='/Users/loisks/Documents/ResearchProjects/ChargingProject/PosChargingPlots'
        if not os.path.exists(subdir_name):
             os.umask(0) # unmask if necessary
             os.makedirs(subdir_name) 
        os.chdir(subdir_name)#
        fig.set_size_inches(13,9)
        plt.savefig('IonLine_'+date+'_ion_'+str(sat[isat])+'_spectrogram.png')
        plt.close(fig)
        os.chdir('..')
   #  except(OSError):
   #     print 'bad date: ' + date  
#
#
#
#
# COUNTS spectrogram
# generate ion spectrogram
# IONSSSSSS
       #try:
        date=uniqueDates[iDate]
        print(date)
       # first get EFW
       #os.chdir('EFW_L3_'+sat[isat])
       #f=glob.glob('*'+uniqueDates[iDate]+'*')
       #pyf=pycdf.CDF(f[0])
       #density=(list(pyf['density'][...]))
       #potential=pyf['Vavg'][...]
       #epochE=dt.date2num(pyf['epoch'][...])
       #MLTE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]))
       #LE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]))
       ##
       #MLAT=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
       #
       # now get HOPE data
        os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/HOPE_L3_'+sat[isat])
        f=glob.glob('*'+uniqueDates[iDate]+'*')
        pyf=pycdf.CDF(f[0])
        
        counts=pyf['Counts_P_Omni'][...]
        electronEnergy=pyf['HOPE_ENERGY_Ion'][...]
        width=np.nanmean(pyf['ENERGY_Ion_DELTA'][...],axis=0)
        err=np.zeros(len(chargingLine))
        for iChar in range(len(chargingLine)):
            err[iChar]=min(width, key=lambda x:abs(x-chargingLine[iChar])) 
        #ionFlux=electronFlux*electronEnergy
        epoch=pyf['Epoch_Ion'][...]
        L=pyf['L_Ion'][...]
        MLT=pyf['MLT_Ion'][...]
        #MLAT=pyf['MLAT_Ele'][...]
        os.chdir('..')
        mTime=dt.date2num(epoch)


        
        #
        # now we have both
        fig=plt.figure()
        cmap=plt.cm.jet
        cbmin=0
        cbmax=500
        ax1=fig.add_subplot(111)
        plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
        ax1.set_ylabel('Energy (eV)', fontsize=25, fontweight='bold')
        ax1.set_xlabel('Time', fontsize=25, fontweight='bold')
        ax1.set_yscale('log')
        mEnergies=HIE
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
        cb.set_label('Counts', fontsize=25, fontweight='bold')
        ax1.tick_params(axis='both', which='major', labelsize=22)
        cb.ax.tick_params(labelsize=30)
        hfmt = dt.DateFormatter('%H:%M')
        ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
        ax1.xaxis.set_major_formatter(hfmt)
        ax1.set_ylim(3, 45000)
        ax1.set_xlim(epoch[startHope], epoch[endHope])
        
        
        # also overplot S/C Potential
        ax2=ax1.twinx()
        ax2.plot(chargingTimes,-1*np.array(chargingLine) , lw=6, c='Silver')
        ax2.plot(chargingTimes,-1*np.array(chargingLine), lw=2, c='Black')
       # ax2.errorbar(chargingTimes, np.array(err), fmt='o')
        ax2.plot(epochE,-1* EFWpotential, lw=3, c='red')
        ax2.set_ylabel('S/C Potential [V]', fontsize=25, fontweight='bold')
        #ax2.set_yscale('log')
        ax2.set_ylim(-1000,0)
        ax2.set_xlim(epoch[startHope], epoch[endHope])
        
        # Add L LABEL
        ax3 = fig.add_axes(ax1.get_position())
        ax3.patch.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.set_yscale('log')
        for spinename, spine in ax3.spines.items():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax3.spines['bottom'].set_position(('outward', 65))
        L_conc=[round(x, 2) for x in L]
        p3=ax3.plot(mTime[startHope:endHope],np.swapaxes(ionFlux,1,0)[0][startHope:endHope]  , 'g', alpha=0)
        # I hate doing this but
        spacing = int(len(L_conc[startHope:endHope])/7)
        spacedArray=L_conc[startHope:endHope][::spacing]
        ax3.set_xticks(epoch[startHope:endHope][::spacing])
        nans=np.where(np.array(L_conc[startHope:endHope][::spacing]) < 0)
        if len(nans[0])>0:
            try:
                newSpot = L_conc[startHope:endHope][spacing*7+5]
                print(newSpot)
                if newSpot < 0:
                    newSpot=L_conc[startHope:endHope][spacing*7-5]
                    print(newSpot)
                spacedArray[nans[0][0]]=newSpot
            except(IndexError):
                newSpot=L_conc[startHope:endHope][spacing*7-5]
                print(newSpot)
                spacedArray[nans[0][0]]=newSpot
        ax3.set_xticklabels(spacedArray)
        ax3.set_xlabel('L', fontsize=22, fontweight='bold')
        #
        # Add MLT Label
        ax5 = fig.add_axes(ax1.get_position())
        ax5.patch.set_visible(False)
        ax5.yaxis.set_visible(False)
        ax5.set_yscale('log')
        for spinename, spine in ax3.spines.items():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax5.spines['bottom'].set_position(('outward', 125))
        MLT_conc=[round(x, 1) for x in MLT]
        p4=ax5.plot(mTime[startHope:endHope],  np.swapaxes(ionFlux,1,0)[0][startHope:endHope], 'purple', alpha=0)
        spacing = int(len(MLT_conc[startHope:endHope])/8)
        ax5.set_xticks(mTime[startHope:endHope][::spacing])
        ax5.set_xticklabels(MLT_conc[startHope:endHope][::spacing])
        ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
        
        subdir_name='/Users/loisks/Documents/ResearchProjects/ChargingProject/PosChargingPlots'
        if not os.path.exists(subdir_name):
             os.umask(0) # unmask if necessary
             os.makedirs(subdir_name) 
        os.chdir(subdir_name)#
        fig.set_size_inches(13,9)
        plt.savefig('IonLine_counts_'+date+'_ion_'+str(sat[isat])+'_spectrogram.png')
        plt.close(fig)
        os.chdir('..')
#
#
#
# INVERT THE S/C CHARGING
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
        mEnergies=HIE
        X,Y=np.meshgrid(mTime, mEnergies)
        mdata=ma.masked_invalid(np.log10(ionFlux).transpose())
        col=ax1.pcolormesh(X,Y,mdata,cmap='jet', vmax=cbmax, vmin=cbmin)
        font = {'family' : 'normal',
                 'weight' : 'bold',
                 'size'   : 22}
        plt.rc('font', **font)
        cbaxes = fig.add_axes([0.82, 0.27, 0.03, 0.65])
        #cb=plt.colorbar(col,cax=cbaxes)
        cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(cbmin,cbmax+1))
        cb.set_label('cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$', fontsize=25, fontweight='bold')
        ax1.tick_params(axis='both', which='major', labelsize=22)
        cb.ax.tick_params(labelsize=30)
        hfmt = dt.DateFormatter('%H:%M')
        ax1.xaxis.set_major_locator(dt.HourLocator(interval=hourGap))
        ax1.xaxis.set_major_formatter(hfmt)
        ax1.set_ylim(3, 45000)
        ax1.set_xlim(epoch[startHope], epoch[endHope])
        
        
        # also overplot S/C Potential
        ax2=ax1.twinx()
        ax2.plot(chargingTimes,np.array(chargingLine) , lw=6, color='Silver')
        ax2.plot(chargingTimes,np.array(chargingLine), lw=2, color='Black')
        #ax2.errorbar(chargingTimes,np.array(err), fmt='o')
        ax2.plot(epochE, EFWpotential, lw=3, color='red')
        ax2.set_ylabel('Negative S/C Potential [V]', fontsize=25, fontweight='bold')
        ax2.set_yscale('log')
        ax2.set_ylim(3,45000)
        ax2.set_xlim(epoch[startHope], epoch[endHope])
        #ax2.set_ylim(-1000,0)
        
        # Add L LABEL
        ax3 = fig.add_axes(ax1.get_position())
        ax3.patch.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.set_yscale('log')
        for spinename, spine in ax3.spines.items():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax3.spines['bottom'].set_position(('outward', 65))
        L_conc=[round(x, 2) for x in L]
        p3=ax3.plot(mTime[startHope:endHope],np.swapaxes(ionFlux,1,0)[0][startHope:endHope]  , 'g', alpha=0)
        # I hate doing this but
        spacing = int(len(L_conc[startHope:endHope])/7)
        spacedArray=L_conc[startHope:endHope][::spacing]
        ax3.set_xticks(epoch[startHope:endHope][::spacing])
        nans=np.where(np.array(L_conc[startHope:endHope][::spacing]) < 0)
        if len(nans[0])>0:
            try:
                newSpot = L_conc[startHope:endHope][spacing*7+5]
                print(newSpot)
                if newSpot < 0:
                    newSpot=L_conc[startHope:endHope][spacing*7-5]
                    print(newSpot)
                spacedArray[nans[0][0]]=newSpot
            except(IndexError):
                newSpot=L_conc[startHope:endHope][spacing*7-5]
                print(newSpot)
                spacedArray[nans[0][0]]=newSpot
        ax3.set_xticklabels(spacedArray)
        ax3.set_xlabel('L', fontsize=22, fontweight='bold')
        #
        # Add MLT Label
        ax5 = fig.add_axes(ax1.get_position())
        ax5.patch.set_visible(False)
        ax5.yaxis.set_visible(False)
        ax5.set_yscale('log')
        for spinename, spine in ax3.spines.items():
            if spinename != 'bottom':
                spine.set_visible(False)
        ax5.spines['bottom'].set_position(('outward', 125))
        MLT_conc=[round(x, 1) for x in MLT]
        p4=ax5.plot(mTime[startHope:endHope],  np.swapaxes(np.log10(ionFlux),1,0)[0][startHope:endHope], 'purple', alpha=0)
        spacing = int(len(MLT_conc[startHope:endHope])/8)
        ax5.set_xticks(mTime[startHope:endHope][::spacing])
        ax5.set_xticklabels(MLT_conc[startHope:endHope][::spacing])
        ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
        
        subdir_name='/Users/loisks/Documents/ResearchProjects/ChargingProject/PosChargingPlots'
        if not os.path.exists(subdir_name):
             os.umask(0) # unmask if necessary
             os.makedirs(subdir_name) 
        os.chdir(subdir_name)#
        fig.set_size_inches(13,9)
        plt.savefig('IonLine_inverse_'+date+'_ion_'+str(sat[isat])+'_spectrogram_hr1='+str(hour1)+'_'+str(hour2)+'.png')
        plt.close(fig)
        os.chdir('..')
        print( np.min(chargingLine))
        print( np.nanmax(chargingLine))
