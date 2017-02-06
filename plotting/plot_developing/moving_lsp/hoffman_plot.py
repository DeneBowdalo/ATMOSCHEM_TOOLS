#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
#from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import stats
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import modules
from collections import OrderedDict
import operator
from netCDF4 import Dataset
from scipy.odr import Model, RealData, ODR, Data
from pylab import *
import pandas as pd
import matplotlib.dates as dates
from collections import Counter
from itertools import groupby

#-----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()

species = 'O3'
start_year = 1990
end_year = 2014

years = np.arange(start_year,end_year+1)

fig, (ax)  = plt.subplots(1,figsize=(19,12))
fig.patch.set_facecolor('white')

#find constant refs in all files
obs_data = []

obs_fname = '/work/home/db876/observations/surface/%s/process/EMEP_SURFACE_%s_%s_%s_gapsallowed.nc'%(species,species,start_year,end_year+1)
obs_ts_grp = Dataset(obs_fname)
obs_site_group = obs_ts_grp.groups['no0042g'] 
data = obs_site_group.variables['%s'%(species.lower())][:]

#get average seasonal waveform

s = start_year
e = start_year+5

#all_seasonal_waveforms = np.empty((2010-1990,8766))

#count = 0
#while e < end_year+2:
#    print s,e
#    s_group = Dataset('obs_sig_periods_%s_%s.nc'%(s,e))
#    s_data = s_group.groups['no0042g']
#    seasonal_waveform = s_data.variables['seasonal_waveform'][:]
#    s+=1
#    e+=1

#    all_seasonal_waveforms[count,:] = seasonal_waveform
#    count+=1

#ave_seasonal_waveform = np.average(all_seasonal_waveforms,axis=0)

s = 0
#test = data < 0
#data[test] = np.nan

for c in range(len(years)):
    year = years[c]

    d0 = datetime.datetime(year,1,1,0,0)
    d1 = datetime.datetime(year+1,1,1,0,0)
    delta = d1 - d0
    n_days = delta.days
    n_hours = n_days*24

    e = s+n_hours
    year_data = np.array([data[s:e]])


    x_axis = np.arange(0,n_hours+1,1)                                                                                                                                                                       
    y_axis = np.array([year,year+1])

    X, Y = np.meshgrid(x_axis, y_axis)

    print year_data.shape

    #remove average / or ave seasonal waveform
#    if year_data.shape[1] == 8760:
#        asw = ave_seasonal_waveform[:8760]
#        year_data[0,:] = year_data[0,:]-asw
#    if year_data.shape[1] == 8784:
#        asw = np.append(ave_seasonal_waveform,[ave_seasonal_waveform[-1]]*18)    
#        year_data[0,:] = year_data[0,:]-asw

    pl = plt.pcolor(X,Y,year_data,cmap=plt.get_cmap('coolwarm'),vmin=0, vmax=60)
    pl.cmap.set_under('white')    

    s = e

plt.xlim(0,8784)    
plt.ylim(years[0],years[-1]+1)
plt.colorbar()

plt.tight_layout(pad = 3.08)

plt.show()
