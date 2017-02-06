import modules
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import signal
from netCDF4 import Dataset
import datetime
import pandas as pd

data_type = raw_input('amp, ph or ave?\n')
if data_type != 'ave':
    period = raw_input('\nd or s?\n')
else:
    period = ''

typ = raw_input('\nabs or pc?\n')


data = Dataset('ACCMIP_hourly_%s%s_std.nc'%(period,data_type))
abs_std = data.variables['absolute_std'][:]
frac_std = data.variables['fractional_std'][:]
lat_c = data.variables['lat_centre'][:]
lon_c = data.variables['lon_centre'][:]
lat_e = data.variables['lat_edges'][:]
lon_e = data.variables['lon_edges'][:]
cesmcam = np.ravel(data.variables['cesmcam'][:])
cmam = np.ravel(data.variables['cmam'][:])
geosccm = np.ravel(data.variables['geosccm'][:])
geoschem = np.ravel(data.variables['geoschem'][:])
gfdl = np.ravel(data.variables['gfdl'][:])
giss = np.ravel(data.variables['giss'][:])
mirocchem = np.ravel(data.variables['mirocchem'][:])

if typ == 'abs':
    z = abs_std
elif typ == 'pc':
    z = frac_std

all_lat_c = [[i]*len(lon_c) for i in lat_c]
all_lat_c = [item for sublist in all_lat_c for item in sublist]
all_lon_c = [lon_c] * len(lat_c) 
all_lon_c = [item for sublist in all_lon_c for item in sublist]

anthfile = '/work/home/db876/plotting_tools/core_tools/anthro2_a2000.nc'
anthload = Dataset(anthfile)
class_result,class_name = modules.anthrome_classify(anthload,all_lat_c,all_lon_c)

areas = ['ANT','OC','S_O','AF','SE_US','S_US','W_US','N_US','NE_US','W_CAN','E_CAN','S_EU','C_EU','NW_EU','N_EU','E_EU','AS','N_O','ARC']

plot_type = raw_input('\nd, s or full?\n')

if plot_type == 'd':
    obs_datetimes = obs_datetimes[:24]
    model_datetimes = model_datetimes[:24]
if plot_type == 's':
    obs_datetimes = obs_datetimes[:8766]
    model_datetimes = model_datetimes[:8766]

obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')



#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111)

models = ['CESMCAM','CMAM','GEOSCCM','GEOSCHEM','GFDL','GISS','MIROCCHEM']
facecolors = ['#7570b3','#771b9e','#1b429e','#9e1b42','#429e1b','#1b9e77','#9e771b']
whiskercolors = ['#7570b3','#7570b3','#771b9e','#771b9e','#1b429e','#1b429e','#9e1b42','#9e1b42','#429e1b','#429e1b','#1b9e77','#1b9e77','#9e771b','#9e771b']

start_n = 1
end_n = 8

for c in classes:
    test = class_name == c
    cesmcam_cut = cesmcam[test]
    cmam_cut = cmam[test]
    geosccm_cut = geosccm[test]
    geoschem_cut = geoschem[test]
    gfdl_cut = gfdl[test]
    giss_cut = giss[test]
    mirocchem_cut = mirocchem[test]
    #n, bins, patches = plt.hist(cesmcam_cut, 50, normed = 1, histtype='step') 
    
    all_data = [cesmcam_cut,cmam_cut,geosccm_cut,geoschem_cut,gfdl_cut,giss_cut,mirocchem_cut]
    
    pos = np.arange(start_n,end_n,1)
    bp = ax.boxplot(all_data,patch_artist=True,widths=1,positions = pos)

    count = 0    
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = facecolors[count])
        count+=1

    count=0
    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color=whiskercolors[count], linewidth=2)
        count+=1

    count=0
    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color=whiskercolors[count], linewidth=2)
        count+=1

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5,markersize=1)
        
    start_n+=8
    end_n+=8
    
    
#-------------------------------------------

mids = np.arange(4,end_n-8,8)
ax.set_xlim([0,end_n-8])
ax.set_xticks(mids)
ax.set_xticklabels(classes)

h1, = ax.plot([1,1],color=facecolors[0])
h2, = ax.plot([1,1],color=facecolors[1])
h3, = ax.plot([1,1],color=facecolors[2])
h4, = ax.plot([1,1],color=facecolors[3])
h5, = ax.plot([1,1],color=facecolors[4])
h6, = ax.plot([1,1],color=facecolors[5])
h7, = ax.plot([1,1],color=facecolors[6])

plt.legend((h1,h2,h3,h4,h5,h6,h7),models)
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)
h4.set_visible(False)
h5.set_visible(False)
h6.set_visible(False)
h7.set_visible(False)

plt.show()

