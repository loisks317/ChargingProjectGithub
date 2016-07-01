# ionLineHOPE.py
#
# determine the ion line and HOPE to calculate spacecraft charging
#
# LKS, September 2015
#
import os
import glob 
import numpy as np
from matplotlib import pyplot as plt
import datetime
import matplotlib.dates as dates
from spacepy import pycdf
import itertools as itert
import math
from numpy import ma
import pandas as pd
import pickle
#
#
dateStart='20130208'
dateEnd='20130209'
dt1=datetime.datetime.strptime(dateStart, '%Y%m%d')
date=dateStart
DT=datetime.datetime.strptime(date, '%Y%m%d')
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
#
# now load in HOPE data for this
chargingLine=[]
chargingTimes=[]
while DT != endDt:
       # approximately from january to may
  date=datetime.datetime.strftime(DT, '%Y%m%d')
  DT=datetime.datetime.strptime(date, '%Y%m%d')
  #
  # get the HOPE data
  #os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_A')
  os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/HOPE_L3_A')
  f=glob.glob('*'+date+'*')
  try:
     pyf=pycdf.CDF(f[0])
  
     counts=pyf['Counts_P_Omni'][...]
     epoch=pyf['Epoch_Ion'][...]
     HIE=pyf['HOPE_ENERGY_Ion'][...]
     Lion=pyf['L_Ion'][...]

     
     # figure out if there is an ion line
     gradient=np.zeros((len(epoch),55))
     for iTime in range(len(epoch)):
       cCount=-1
       if Lion[iTime]>4:
     #   for iEn in range(55):
#
#            if counts[iTime][iEn]==0:
#                cCount+=1
#            elif(counts[iTime][iEn]>=100): # more than 10 counts in the line bin
#                print counts[iTime][iEn]
#                print iEn
#                print epoch[iTime]
#                break

             # try to impose a gradient test

        for iEn in range(55):
                 gradient[iTime][iEn]=((counts[iTime][iEn+1]-counts[iTime][iEn])/(counts[iTime][iEn]*1.0))
        a=np.where(gradient[iTime] > 10)[0]
        try:
                     test=len(a)
                     cCount=a[0]
        except:
                     cCount=0
             


        if cCount >= 1:
            sc=HIE[iTime][cCount]
        else:
            sc=0
        chargingLine.append(sc)
        chargingTimes.append(epoch[iTime])
  except:
      print 'bad date' + str(date)
  DT=DT+datetime.timedelta(days=1)

#os.chdir('/Users/loisks/Desktop/ResearchProjects/ChargingProject/')
os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/')
with open('chargingLine.p','wb') as f:
       pickle.dump(  chargingLine, f)
with open('chargingTimes.p','wb') as f:
       pickle.dump(chargingTimes,f)
