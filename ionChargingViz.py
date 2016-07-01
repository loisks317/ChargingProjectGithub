# ionChargingViz.py
#
# visualize the HOPE ion line
#
# LKS, September 2015
#
import os
import glob 
import numpy as np
from matplotlib import pyplot as plt
import datetime
import matplotlib.dates as dt
import spacepy.pybats.kyoto as spk
from spacepy import pycdf
import itertools as itert
import math
from numpy import ma
import pandas as pd
import pickle
#
with open('chargingLine.p', 'rb') as f:
    data=pickle.load(f)
with open('chargingTimes.p', 'rb') as f:
    time=pickle.load(f)

data=np.array(data)
time=np.array(time)
#
# screen for bad data
time[data>50000]=np.nan
data[data>50000]=np.nan
time=time[~np.isnan(data)]
data=data[~np.isnan(data)]
#
# plot
fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.set_yscale('log')
ax1.set_ylim(1e1,1e5)
ax1.scatter(dt.date2num(time), data, s=30)
hfmt = dt.DateFormatter('%m-%d')
ax1.xaxis.set_major_locator(dt.DayLocator(interval=10))
ax1.xaxis.set_major_formatter(hfmt)
plt.show()
