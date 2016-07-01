# LTbins.py

# bin by MLT and L-Shell for charging
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
binWidth=10
maxCharging=-200
nX=25
chargingBins=np.linspace(maxCharging,50,nX)
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
iDate='20130201'
dt=datetime.datetime.strptime(iDate,'%Y%m%d')
endDate='20150501'
satellites=['A', 'B']

MLTLP=[[ [[] for i in range(nL)] for j in range(nMLT)] for k in range(4)]
tMLTLP=[[ [] for j in range(nL)] for k in range(nMLT)]
for iSat in range(len(satellites)):
 iDate='20130201'
 dt=datetime.datetime.strptime(iDate,'%Y%m%d')
 while iDate != endDate:
#for iDate in range(len(uniqueDates)):
   try:
       iDate=datetime.datetime.strftime(dt,'%Y%m%d')
       dt=datetime.datetime.strptime(iDate, '%Y%m%d')
       if (iDate=='20140428') or (iDate=='20140429'):
           # this is bad data
           skip

       # first get EFW
       os.chdir('EFW_L3_'+satellites[iSat])
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
      #MLTE[eclT0]=np.nan; MLTE[eclT1]=np.nan; MLTE[eclT2]=np.nan
      #LE[eclT0]=np.nan; LE[eclT1]=np.nan; LE[eclT2]=np.nan
      #MLATE[eclT0]=np.nan; MLATE[eclT1]=np.nan; MLATE[eclT2]=np.nan
      #
      # make a simple bar plot of the number events in
      # each charging window
      #
       p=np.array(potential)
      #p[eclT0]=np.nan
      #p[eclT1]=np.nan
      #p[eclT2]=np.nan

       #
       # create four bins, sort by MLT L
       for iMLT in range(nMLT):
           for iL in range(nL):
               MLTindex=np.where((MLTE >= iMLT*0.5) & (MLTE < (iMLT+1)*0.5))[0]
               Lindex=np.where((LE >= iL*0.25) & (LE < (iL+1)*0.25))[0]
               A1=np.where(p < -10)[0]
               set1=list(set(MLTindex) & set(Lindex) & set(A1))
               if len(set1) > 0:
                   MLTLP[0][iMLT][iL]+=(list([set1]))
               A2=np.where((p>= -10) & (p < 0))[0]
               set2=list(set(MLTindex) & set(Lindex) & set(A2))
               MLTLP[1][iMLT][iL]+=(list(p[set2]))              
               A3=np.where((p>=0) & (p< 10))[0]
               set3=list(set(MLTindex) & set(Lindex) & set(A3))
               MLTLP[2][iMLT][iL]+=(list(p[set3]))                 
               A4=np.where(p>=10)[0]
               set4=list(set(MLTindex) & set(Lindex) & set(A4))
               MLTLP[3][iMLT][iL]+=(list(p[set4]))
               tMLTLP[iMLT][iL]+=(list(set(MLTindex) & set(Lindex)))
               
       os.chdir('..')
       dt=dt+datetime.timedelta(days=1)
   except:
       os.chdir('..')
       dt=dt+datetime.timedelta(days=1)           
       #
       # Now sort by threshold
# get the lengths
os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject')
nMLTLP=[[[ 0 for i in range(nL)] for j in range(nMLT)] for k in range(4)]
for iMLT in range(nMLT):
    for iL in range(nL):
        for iN in range(4):
            if len(tMLTLP[iMLT][iL]) > 0:
                nMLTLP[iN][iMLT][iL]=len(MLTLP[iN][iMLT][iL])/(1.0*len(tMLTLP[iMLT][iL]))

Labels=["< -10", "-10 to 0", "0 to 10", "> 10"]
for iN in range(4):
       # intitialize the figure
        fig=plt.figure()
        ax=fig.add_subplot(111, polar=True)
        datah_m=ma.masked_invalid(np.array(nMLTLP[iN]).transpose())
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
        subdir_name='MLTL_Dialplots'
        if not os.path.exists(subdir_name):
           os.umask(0) # unmask if necessary
           os.makedirs(subdir_name, 0777) 
        os.chdir(subdir_name)#
        fig.set_size_inches(13,9)
        plt.savefig('MLTL_Dialplots_ChargingWindow_' + Labels[iN]+'.pdf')
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
        ax2.plot(np.linspace(0,24,nMLT), np.swapaxes(nMLTLP[1],1,0)[9], lw=2)
        subdir_name='MLTL_LinePlots'
        if not os.path.exists(subdir_name):
           os.umask(0) # unmask if necessary
           os.makedirs(subdir_name, 0777) 
        os.chdir(subdir_name)#
        fig3.set_size_inches(13,9)
        plt.savefig('ChargingWindow_' + Labels[1]+'.png')
        plt.close(fig3)
        os.chdir('..')  



   
#fig2=plt.figure()
#ax2=fig2.add_subplot(111)
#ax2.set_yscale('log')
#ax2.set_xlim(-200,50)
#ax2.set_ylim(1e0,1e5)
#ax2.bar(binsF,binnedP, 10, color='darkviolet')
#ax2.set_title("Charging Magnitudes \n Across Mission")
#font = {'family' : 'normal',
#          'weight' : 'bold',
#          'size'   : 22}
#plt.rc('font', **font)
#plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.15)
#ax2.set_ylabel('Occurrences', fontweight='bold')
#ax2.set_xlabel('Charging Window Lower Bound',fontweight='bold')
#os.chdir('MLTL_Dialplots')              
#plt.savefig("ChargingMagnitudes.pdf")
#os.chdir('..')
