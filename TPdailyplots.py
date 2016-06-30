# TPdailyplots.py
#
# look at temperature and pressure values on a daily scale compared to s/c
# charging
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
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
#
# start the date loop
uniqueDates=['20130214']
# '20130113',
# '20130114',
# '20130117',
# '20130120',
# '20130126',
# '20130202',
# '20130208',
# '20130213',
# '20130214',
# '20130216',
# '20130217',
# '20130220',
# '20130222',
# '20130223',
# '20130226',
# '20130301',
# '20130316',
# '20130317',
# '20130320',
# '20130321',
# '20130323',
# '20130327',
# '20130328',
# '20130329',
# '20130404',
# '20130407',
# '20130424',
# '20130426',
# '20130501']
# data
hourGap=5
sat=['A', 'B']
corrCoeffTpar=[]
corrCoeffTperp=[]
corrCoeffP=[]
dateArr=[]
for date in uniqueDates:
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
    MLATE=(list(np.swapaxes(pyf['mlt_lshell_mlat'][...],1,0)[2]))
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
    epoch=dates.date2num(pyf['Epoch_Ele'][...])
    density[density<0]=np.nan
    Tperp[Tperp>1e24]=np.nan
    pressure=density*1e6*boltzmann*(1/3.*Tpar + 2/3.*Tperp) # check this
    temperature = (1/3.*Tpar + 2/3.*Tperp)

    # interpolate onto potential times
    # calculate correlation coefficient for s/c charging and temperature
    os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject')
    #
    # now plot electron density, Tpar, s/c charging
    fig=plt.figure()
    ax1= host_subplot(111, axes_class=AA.Axes)
    ax1.set_ylabel('Temperature [K]')
    ax1.set_yscale('log')
    plt.subplots_adjust(left=0.11, right=0.70, top=0.92, bottom=0.4)
    font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}
    plt.rc('font', **font)
    #p1,=ax1.plot(epoch, Tpar, c='r', lw=2, label="Tpar")
    p1,=ax1.plot(epoch, temperature, c='r', lw=2, label="Temperature")
    hfmt = dates.DateFormatter('%H:%M')
    ax1.xaxis.set_major_locator(dates.HourLocator(interval=hourGap))
    ax1.xaxis.set_major_formatter(hfmt)
    ax2=ax1.twinx()
    ax2.set_yscale('log')
    ax2.set_ylabel('Pressure [nPa]')
    p2,=ax2.plot(epoch,np.array( pressure)/1.0e-9, c='g', lw=2, label="P")
    axE=ax1.twinx()
    axE.set_ylabel('Potential [-V]')
    p3,=axE.plot(epochE, -1*potential, c='Orange', lw=2, label='$\phi$')
    offset = 100
    new_fixed_axis = axE.get_grid_helper().new_fixed_axis
    axE.axis["right"] = new_fixed_axis(loc="right",
                                        axes=axE,
                                        offset=(offset, 0))
    axn=ax1.twinx()
    axn.set_ylabel('Density [cm$^-3$]')
    axn.set_yscale('log')
    p6,=axn.plot(epoch, density, c='Blue', lw=2, label='density')
    offset = 200
    new_fixed_axis = axn.get_grid_helper().new_fixed_axis
    axn.axis["right"] = new_fixed_axis(loc="right",
                                        axes=axn,
                                        offset=(offset, 0))

    axE.axis["right"].toggle(all=True)
    ax1.axis["left"].label.set_color(p1.get_color())
    ax2.axis["right"].label.set_color(p2.get_color())
    axE.axis["right"].label.set_color(p3.get_color())
    axn.axis["right"].toggle(all=True)
    axn.axis["right"].label.set_color(p6.get_color())

    # line plot should be ready to go
    

    
    
    #
    # Add L LABEL
    ax3 = fig.add_axes(ax1.get_position())
    ax3.patch.set_visible(False)
    ax3.yaxis.set_visible(False)
    ax3.set_yscale('log')
    for spinename, spine in ax3.spines.items():
        if spinename != 'bottom':
            spine.set_visible(False)
    ax3.spines['bottom'].set_position(('outward', 65))
    L_conc=[round(x, 2) for x in LE]
    p3=ax3.plot(epochE,-1*potential  , 'g', alpha=0)
    # I hate doing this but
    spacing = int(len(L_conc)/7)
    spacedArray=L_conc[::spacing]
    ax3.set_xticks(epochE[::spacing])
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
    MLT_conc=[round(x, 1) for x in MLTE]
    p4=ax5.plot(epochE,  potential, 'purple', alpha=0)
    spacing = int(len(MLT_conc)/8)
    ax5.set_xticks(epochE[::spacing])
    ax5.set_xticklabels(MLT_conc[::spacing])
    ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
    
    # Add MLAT Label
    ax6 = fig.add_axes(ax1.get_position())
    ax6.patch.set_visible(False)
    ax6.yaxis.set_visible(False)
    ax6.set_yscale('log')
    for spinename, spine in ax3.spines.items():
        if spinename != 'bottom':
            spine.set_visible(False)
    ax6.spines['bottom'].set_position(('outward', 185))
    MLAT_conc=[round(x, 1) for x in MLATE]
    p5=ax6.plot(epochE,  potential, 'purple', alpha=0)
    spacing = int(len(MLAT_conc)/8)
    ax6.set_xticks(epochE[::spacing])
    ax6.set_xticklabels(MLAT_conc[::spacing])
    ax6.set_xlabel('MLAT [Degrees]', fontsize=22, fontweight='bold')
    #
    # save the figure
    subdir_name='Temp_P_Compare'
    if not os.path.exists(subdir_name):
         os.umask(0) # unmask if necessary
         os.makedirs(subdir_name) 
    os.chdir(subdir_name)#
    fig.set_size_inches(13,9)
    plt.savefig(date+'_'+str(sat[iSat])+'_tempPcompare.pdf')
    plt.close(fig)
    os.chdir('..')
