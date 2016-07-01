# LTbins.py

# bin by MLT and L-Shell for charging
# October 2015 - bin into 1 V bins in the 0 to V bin
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
#
# create a half circle for density maps
def dual_half_circle(center, radius, angle=90, ax=None, colors=('w','k'),
                     **kwargs):
    from matplotlib.patches import Wedge
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0],transform=ax.transData._b, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], transform=ax.transData._b,**kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]
#
# establish bin sizes
binWidth=1
maxCharging=0
nX=10
chargingBins=np.linspace(maxCharging,9,nX)
nMLT=49
nL=29
nMLAT=41
MLTbins=np.linspace(0, 24, nMLT)
MLTarr=[ [] for i in range(nMLT)]
Larr=[ [] for i in range(nL)]
MLATarr=[ [] for i in range(nMLAT)]
Lbins=np.linspace(0.25, 7.25, nL)
MLATbins=np.linspace(-20,20, nMLAT)
MLTLBins=[[ [] for i in range(nL)] for j in range(nMLT)]
lenMLTL=[[ 0 for i in range(nL)] for j in range(nMLT)]

#
# unique dates for less than 50V
# important to normalize this


# now load in HOPE data for this
interpFlux=[[] for i in range(72)]
allPotential=[]
allMLT=[]
allL=[]
allMLAT=[]
TotalFluence=[[] for i in range(72)]
allPotential=[]
binsF=np.linspace(-200,40,25)
binnedP=[ 0 for i in range(25)]
iDate='20130101'
dt=datetime.datetime.strptime(iDate,'%Y%m%d')
endDate='20150501'
while iDate != endDate:
#for iDate in range(len(uniqueDates)):
   try:
       iDate=datetime.datetime.strftime(dt,'%Y%m%d')
       dt=datetime.datetime.strptime(iDate, '%Y%m%d')

       # first get EFW
       os.chdir('EFW_L3_A')
       f=glob.glob('*'+iDate+'*')
       pyf=pycdf.CDF(f[0])
       density=(list(pyf['density'][...]))
       potential=-1*pyf['Vavg'][...]
       epochE=dates.date2num(pyf['epoch'][...])
       MLTE=np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[0]
       LE=np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[1]
       MLATE=np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]
       eclipseFlags0=pyf['flags_charging_bias_eclipse'][...][:,0]
       eclipseFlags1=pyf['flags_charging_bias_eclipse'][...][:,1]
       eclipseFlags2=pyf['flags_charging_bias_eclipse'][...][:,2]
       eclT0=np.where(eclipseFlags0==1)[0]
       eclT1=np.where(eclipseFlags1==1)[0]
       eclT2=np.where(eclipseFlags2==1)[0]

       #
       # nan the eclipse times
       MLTE[eclT0]=np.nan; MLTE[eclT1]=np.nan; MLTE[eclT2]=np.nan
       LE[eclT0]=np.nan; LE[eclT1]=np.nan; LE[eclT2]=np.nan
       MLATE[eclT0]=np.nan; MLATE[eclT1]=np.nan; MLATE[eclT2]=np.nan
      #
      # make a simple bar plot of the number events in
      # each charging window
      #
       p=np.array(potential)
       p[eclT0]=np.nan
       p[eclT1]=np.nan
       p[eclT2]=np.nan
       for iP in range(9):
         temp=np.where( (p >= iP) & (p < iP+1))[0]
         binnedP[iP]+=len(temp)

       
       # identify points of charging
       #
       # now bin everything
       for i in range(nMLT):
           try:
               MLTindex=np.where((MLTE >= i*0.5) & (MLTE < (i+1)*0.5))[0]
               MLTarr[i].extend(list(potential[MLTindex]))
           except:
               MLTarr[i].extend([])
       for i in range(nL):
           try:
               Lindex=np.where((LE >= i*0.25) & (LE < (i+1)*0.25))[0]
               Larr[i].extend(list(potential[Lindex]))
           except:
               Larr[i].extend([])
       for i in range(nMLAT):
           try:
               MLATindex=np.where((MLATE >= i*1-10) & (MLATE < (i+1)-10))[0]
               MLATarr[i].extend(list(potential[MLATindex]))
           except:
               MLATarr[i].extend([])
       LE=np.array(LE); MLTE=np.array(MLTE)
       for k in range(nL):
          for j in range(nMLT):
             try:
               multindex=np.where((MLTE >= j*0.5) & (MLTE < (j+1)*0.5) & (LE >= k*0.25) & (LE < (k+1)*0.25) )[0]
               MLTLBins[j][k].extend(list(potential[multindex]))
               #lenMLTL[j][k]=len(list(potential[multindex]))
               
             except:
               MLTLBins[j][k].extend([])
               #lenMLTL[j][k]=0
             
                  
       os.chdir('..')
       dt=dt+datetime.timedelta(days=1)
   except(IndexError):
       os.chdir('..')
       dt=dt+datetime.timedelta(days=1)           
       #
       # Now sort by threshold
# get the totals
tMLTarr=[ [] for i in range(nMLT)]
tLarr=[ [] for i in range(nL)]
tMLATarr=[ [] for i in range(nMLAT)]
for iMLT in range(nMLT):
    try:
        tMLTarr[iMLT]=len(np.array(MLTarr[iMLT]))*1.0
    except:
        tMLTarr[iMLT]=0
for iL in range(nL):
    try:
        tLarr[iL]=len(np.array(Larr[iL]))*1.0
    except:
        tLarr[iL]=0
for iMLAT in range(nMLAT):
       try:
                tMLATarr[iMLAT]=len(np.array(MLATarr[iMLAT]))*1.0
       except:
                tMLATarr[iMLAT]=0

for iMLT in range(nMLT):
   for iL in range(nL):
      try:
         lenMLTL[iMLT][iL]=len(np.array(MLTLBins[iMLT][iL]))*1.0
      except(OSError):
         lenMLTL[iMLT][iL]=0
                      
       
for iBar in range(nX):
        sMLTarr=[ [] for i in range(nMLT)]
        sLarr=[ [] for i in range(nL)]
        sMLATarr=[ [] for i in range(nMLAT)]
        oMLTarr=[ [] for i in range(nMLT)]
        oLarr=[ [] for i in range(nL)]
        oMLATarr=[ [] for i in range(nMLAT)]
        sMLTL=[[ [] for i in range(nL)] for j in range(nMLT)]
        for iMLT in range(nMLT):
            indexes=np.where((np.array(MLTarr[iMLT])>=maxCharging+iBar*binWidth) & (np.array(MLTarr[iMLT]) < maxCharging+(1+iBar)*binWidth))[0]
            try:
                sMLTarr[iMLT]=len(np.array(MLTarr[iMLT])[indexes])/tMLTarr[iMLT]
                oMLTarr[iMLT]=len(np.array(MLTarr[iMLT])[indexes])
            except:
                sMLTarr[iMLT]=0
                oMLTarr[iMLT]=0
        for iL in range(nL):
            indexes=np.where((np.array(Larr[iL])>=maxCharging+iBar*binWidth) & (np.array(Larr[iL]) < maxCharging+(1+iBar)*binWidth))[0]
            try:
                sLarr[iL]=len(np.array(Larr[iL])[indexes])/tLarr[iL]
                oLarr[iL]=len(np.array(Larr[iL])[indexes])
            except:
                sLarr[iL]=0
                oLarr[iL]=0
        for iMLAT in range(nMLAT):
            indexes=np.where((np.array(MLATarr[iMLAT])>=maxCharging+iBar*binWidth) & (np.array(MLATarr[iMLAT]) < maxCharging+(1+iBar)*binWidth))[0]
            try:
                sMLATarr[iMLAT]=len(np.array(MLATarr[iMLAT])[indexes])/tMLATarr[iMLAT]
                oMLATarr[iMLAT]=len(np.array(MLATarr[iMLAT])[indexes])
            except:
                sMLATarr[iMLAT]=0
                oMLATarr[iMLAT]=0
        for iMLT in range(nMLT):
           for iL in range(nL):
            indexes=np.where((np.array(MLTLBins[iMLT][iL])>=maxCharging+iBar*binWidth) & (np.array(MLTLBins[iMLT][iL]) < maxCharging+(1+iBar)*binWidth) )[0]
            try:
                sMLTL[iMLT][iL]=len(np.array(MLTLBins[iMLT][iL])[indexes])/lenMLTL[iMLT][iL]
                
            except:
                sMLTL[iMLT][iL]=0

                
        # plot bar for each bin size
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.set_xlim(0,24)
        ax2.bar(MLTbins, sMLTarr, 0.5, color='b')
        ax2.set_title(str(maxCharging+iBar*binWidth) + ' V to ' + str(maxCharging+(1+iBar)*binWidth) +' V')
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Fraction of Total', fontweight='bold')
        ax2.set_xlabel('MLT ',fontweight='bold')
        subdir_name='MLT_barplots_0_10V'
        if not os.path.exists(subdir_name):
               os.umask(0) # unmask if necessary
               os.makedirs(subdir_name, 0777) 
        os.chdir('MLT_barplots_0_10V')              
        plt.savefig("Norm_MLT_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        os.chdir('..')

        # plot bar for L
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.bar(Lbins, sLarr, 0.25, color='b')
        ax2.set_title(str(maxCharging+iBar*binWidth) + ' V to ' + str(maxCharging+(1+iBar)*binWidth)+' V')
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Fraction of Total', fontweight='bold')
        ax2.set_xlabel('L-Shell ',fontweight='bold')
        subdir_name='L_barplots_0_10V'
        if not os.path.exists(subdir_name):
               os.umask(0) # unmask if necessary
               os.makedirs(subdir_name, 0777) 
        os.chdir('L_barplots_0_10V')              
        plt.savefig("Norm_L_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        os.chdir('..')

        # plot bar for MLAT
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        #ax2.set_yscale('log')
        ax2.set_xlim(-20,10)
        #ax2.set_ylim(0.9,1)
        ax2.bar(MLATbins, sMLATarr, 1, color='lightseagreen')
        ax2.set_title(str(maxCharging+iBar*binWidth) + ' V to ' + str(maxCharging+(1+iBar)*binWidth)+' V')
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Fraction of Total', fontweight='bold')
        ax2.set_xlabel('MLAT [Degrees]',fontweight='bold')
        subdir_name='MLAT_barplots_0_10V'
        if not os.path.exists(subdir_name):
               os.umask(0) # unmask if necessary
               os.makedirs(subdir_name, 0777) 
        os.chdir('MLAT_barplots_0_10V')              
        plt.savefig("Norm_MLAT_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        os.chdir('..')

  


        ####### Non normalized
                # plot bar for each bin size
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.set_xlim(0,24)
        ax2.bar(MLTbins, oMLTarr, 0.5, color='b')
        ax2.set_title(" Charging From " + str(maxCharging+iBar*binWidth) + ' to ' + str(maxCharging+(1+iBar)*binWidth))
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Occurrences', fontweight='bold')
        ax2.set_xlabel('MLT ',fontweight='bold')
        os.chdir('MLT_barplots_0_10V')              
        plt.savefig("MLT_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        os.chdir('..')

        # plot bar for L
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.bar(Lbins, oLarr, 0.5, color='b')
        ax2.set_title("Charging From " + str(maxCharging+iBar*binWidth) + ' to ' + str(maxCharging+(1+iBar)*binWidth))
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Occurrences', fontweight='bold')
        ax2.set_xlabel('L-Shell ',fontweight='bold')
        os.chdir('L_barplots_0_10V')              
        plt.savefig("L_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        os.chdir('..')

        # plot bar for MLAT
        fig2=plt.figure()
        ax2=fig2.add_subplot(111)
        ax2.set_yscale('log')
        ax2.bar(MLATbins, oMLATarr, 1, color='lightseagreen')
        ax2.set_title("Charging From " + str(maxCharging+iBar*binWidth) + ' to ' + str(maxCharging+(1+iBar)*binWidth))
        font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
        plt.rc('font', **font)
        plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
        ax2.set_ylabel('Occurrences', fontweight='bold')
        ax2.set_xlabel('MLAT ',fontweight='bold')
        os.chdir('MLAT_barplots_0_10V')              
        plt.savefig("MLAT_ChargingWindow_" + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        os.chdir('..')



       # intitialize the figure
        fig=plt.figure()
        ax=fig.add_subplot(111, polar=True)
       # datah_m=ma.masked_invalid(np.array(sMLTL).transpose())
        datah_m=ma.masked_invalid(np.array(sMLTL).transpose())
        mlt=np.array(MLTbins)*(15*np.pi/180.)
        X,Y=np.meshgrid(mlt, Lbins)
        ax.set_ylim(0, 7.25)
        dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
        cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
        vmin=-2
        vmax=0
        col=ax.pcolormesh( X, Y, np.log10(datah_m), cmap='jet', vmin=vmin, vmax=vmax )
        #col=ax.pcolormesh( X, Y, datah_m, cmap='jet', vmin=vmin, vmax=vmax )
        cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
        #cb=plt.colorbar(col,cax=cbaxes)
        plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
        xT2=['', '2', '', '4', '', '6', '', '8']
        xL2=['00','', '06','', '12','', '18']
        cb.set_label("Fraction of Total", fontsize=30)
        ax.tick_params(axis='both', which='major', labelsize=25)    
        cb.ax.tick_params(labelsize=35) 
        ax.set_yticklabels(xT2, fontsize=30)
        ax.set_xticklabels(xL2, fontsize=30)
        ax.grid(True)
        plt.draw()
        subdir_name='MLTL_Dialplots_0_10V'
        if not os.path.exists(subdir_name):
           os.umask(0) # unmask if necessary
           os.makedirs(subdir_name, 0777) 
        os.chdir(subdir_name)#
        fig.set_size_inches(13,9)
        plt.savefig('MLTL_Dialplots_ChargingWindow_' + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        plt.close(fig)
        os.chdir('..')  

        #
        # make a line plot
        #
        fig3=plt.figure()
        ax2=fig3.add_subplot(111)
        #ax2.set_ylim(1e-,1e0)
        ax2.set_yscale('log')
        ax2.set_ylabel('Fraction of Total')
        ax2.set_xlabel('MLT')
        ax2.plot(np.linspace(0,24,nMLT), np.swapaxes(sMLTL,1,0)[9], lw=2)
        subdir_name='MLTL_LinePlots'
        if not os.path.exists(subdir_name):
           os.umask(0) # unmask if necessary
           os.makedirs(subdir_name, 0777) 
        os.chdir(subdir_name)#
        fig3.set_size_inches(13,9)
        plt.savefig('MLTL_Dialplots_ChargingWindow_' + str(maxCharging+iBar*binWidth) + '_' + str(maxCharging+(1+iBar)*binWidth)+'.png')
        plt.close(fig3)
        os.chdir('..')          
   
fig2=plt.figure()
ax2=fig2.add_subplot(111)
ax2.set_yscale('log')
ax2.set_xlim(-200,50)
ax2.set_ylim(1e0,1e5)
ax2.bar(binsF,binnedP, 10, color='darkviolet')
ax2.set_title("Charging Magnitudes")
font = {'family' : 'normal',
          'weight' : 'bold',
          'size'   : 22}
plt.rc('font', **font)
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
ax2.set_ylabel('Occurrences', fontweight='bold')
ax2.set_xlabel('Charging Window Lower Bound [V]',fontweight='bold')
os.chdir('MLTL_Dialplots_0_10V')              
plt.savefig("ChargingMagnitudes.pdf")
os.chdir('..')


