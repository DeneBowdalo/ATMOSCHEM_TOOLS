#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid
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
import cmath

period_type = raw_input('d or s?\n')

ts_files =  glob.glob('*HN.nc')
lsp_files = glob.glob('obs_sig*')

ts_files = [ts_files[35]]
lsp_files = [lsp_files[35]]

present_dir = os.getcwd()
paths = present_dir.split('/')
species = paths[-2]

#set up plot
latlower_setup = [20,30,20,-90]
latupper_setup = [80,72,55,90]
lonwest_setup = [-170,-15,115,-180]
loneast_setup = [-50,35,155,180]
label_out = ['NA','EU','AS','ZZ']
label = ['NA','EU','AS','ROW']

fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')

ax = plt.subplot2grid((2,2), (0,0), colspan=2, rowspan=2)
ax1 = plt.subplot2grid((2,2), (0,0), colspan=1)
ax2 = plt.subplot2grid((2,2), (0,1), colspan=1)
ax3 = plt.subplot2grid((2,2), (1,0), colspan=1)
ax4 = plt.subplot2grid((2,2), (1,1), colspan=1)
ax_list = [ax1,ax2,ax3,ax4]

fig.patch.set_facecolor('white')

#set colorbar for average
cmap=plt.cm.jet
#set max val for colorbar
max_cb = 50

a = Basemap(projection='cyl',llcrnrlat=latlower_setup[0],urcrnrlat=latupper_setup[0],llcrnrlon=lonwest_setup[0],urcrnrlon=loneast_setup[0],resolution='c',ax = ax1)    
b = Basemap(projection='cyl',llcrnrlat=latlower_setup[1],urcrnrlat=latupper_setup[1],llcrnrlon=lonwest_setup[1],urcrnrlon=loneast_setup[1],resolution='c',ax = ax2)
c = Basemap(projection='cyl',llcrnrlat=latlower_setup[2],urcrnrlat=latupper_setup[2],llcrnrlon=lonwest_setup[2],urcrnrlon=loneast_setup[2],resolution='c',ax = ax3)
d = Basemap(projection='cyl',llcrnrlat=latlower_setup[3],urcrnrlat=latupper_setup[3],llcrnrlon=lonwest_setup[3],urcrnrlon=loneast_setup[3],resolution='c',ax = ax4)

a.drawcoastlines()
a.drawmapboundary()
b.drawcoastlines()
b.drawmapboundary()
c.drawcoastlines()
c.drawmapboundary()
d.drawcoastlines()
d.drawmapboundary()

count = 0
for ax in ax_list:
    
    if (count == 0) or (count == 3):
        ax.text(0.03, 0.18, label[count], transform=ax.transAxes,fontsize=30, va='top')
    if (count == 1) or (count == 2):
        ax.text(0.03, 0.97, label[count], transform=ax.transAxes,fontsize=30, va='top')

    count+=1


for f in range(len(ts_files)):
    obs_ts_fname = ts_files[f]
    last_split = obs_ts_fname.split('/')[-1]
    print last_split
    start_year = last_split.split('_')[3] 
    end_year = last_split.split('_')[4]
    obs_lsp_fname = 'obs_sig_periods_%s_%s.nc'%(start_year,end_year) 

    print ts_files[f]
    print obs_lsp_fname
    print start_year,end_year

    obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_gap_inds = modules.read_obs_all(obs_ts_fname,species,start_year,end_year)

    if len(obs_refs) == 0:
        continue

    tags = modules.get_tags(obs_refs)

    obs_daily_mag,obs_daily_mag_spring,obs_daily_mag_summer,obs_daily_mag_autumn,obs_daily_mag_winter,obs_daily_phase\
    ,obs_daily_phase_spring,obs_daily_phase_spring,obs_daily_phase_summer,obs_daily_phase_autumn,obs_daily_phase_winter\
    ,obs_daily_phase_min,obs_daily_phase_min_spring,obs_daily_phase_min_summer,obs_daily_phase_min_autumn,obs_daily_phase_min_winter\
    ,obs_ave,obs_ave_spring,obs_ave_summer,obs_ave_autumn,obs_ave_winter,obs_ave_day,obs_ave_night,obs_seasonal_mag,obs_seasonal_mag_day\
    ,obs_seasonal_mag_night,obs_seasonal_phase,obs_seasonal_phase_day,obs_seasonal_phase_night,obs_seasonal_phase_min\
    ,obs_seasonal_phase_min_day,obs_seasonal_phase_min_night,obs_daily_ff,obs_daily_ff_spring,obs_daily_ff_summer,obs_daily_ff_autumn\
    ,obs_daily_ff_winter,obs_seasonal_ff,obs_seasonal_ff_day,obs_seasonal_ff_night,obs_daily_waveform,obs_daily_waveform_spring\
    ,obs_daily_waveform_summer,obs_daily_waveform_autumn,obs_daily_waveform_winter,obs_seasonal_waveform,obs_seasonal_waveform_day\
    ,obs_seasonal_waveform_night,obs_full_waveform,obs_periodic_variance_daily,obs_periodic_variance_daily_spring,obs_periodic_variance_daily_summer\
    ,obs_periodic_variance_daily_autumn,obs_periodic_variance_daily_winter,obs_periodic_variance_seasonal,obs_periodic_variance_seasonal_day\
    ,obs_periodic_variance_seasonal_night,obs_periodic_variance_all = modules.get_obs_specific(obs_lsp_fname,obs_refs)

    print obs_seasonal_phase

    #obs_seasonal_phase = obs_seasonal_phase[:4]
    #obs_seasonal_mag = obs_seasonal_mag[:4]
    #obs_refs = obs_refs[:4]
    #obs_lons = obs_lons[:4]
    #obs_lats = obs_lats[:4]
    #tags = tags[:4]
    #print obs_seasonal_phase
    #print tags
    
    #obs_seasonal_phase = np.array([6.])
    #obs_seasonal_mag = np.array([100.])
    #obs_refs = obs_refs[:1]
    #obs_lons = obs_lons[:1]
    #obs_lats = obs_lats[:1]
    #tags = tags[:1]
    
    #scale all seasonal mags by max 

    #convert phase from time to radians (0 to 2pi)
    #obs_seasonal_phase = ((2.*np.pi)/12.)*obs_seasonal_phase
    
    #phase vectors visually start from 3 on a clockface and go anti-clockwise
    #correct to make go from what is 12 on clockface, and go clockwise
    #test1 = (obs_seasonal_phase >= 0) & (obs_seasonal_phase < np.pi)
    #test2 = (obs_seasonal_phase >= np.pi) & (obs_seasonal_phase < (np.pi*2.))
    #obs_seasonal_phase[test1] = ((np.pi*2.) - obs_seasonal_phase[test1])
    #obs_seasonal_phase[test2] = ((np.pi*2.) - obs_seasonal_phase[test2]) 
    #obs_seasonal_phase = obs_seasonal_phase + np.pi/2.
    
    
    
    
    #plot
    #---------------------------------------
    count = 0

    for ax in ax_list:

        if count == 3:
            test = []
            other_tags = ['AF','ANT','ARC','OC','O','SA']
            for i in tags:
                if i in other_tags:
                    test.append(True)
                else:
                    test.append(False)
            test = np.array(test)
        else:
            current_tag = label_out[count] 
            test = tags == current_tag
    
        current_lons = obs_lons[test]
        current_lats = obs_lats[test]
        
        if count == 0:
            X,Y = a(current_lons,current_lats)
            for i in range(len(current_lons)):
                #convert amplitude and phase to complex real and imag parts
                obs_complex = cmath.rect(obs_seasonal_mag[i], obs_seasonal_phase[i])
                all1 = a.quiver(X[i],Y[i],obs_complex.real,obs_complex.imag,width=0.005,color=cmap((1./max_cb)*obs_ave[i]))
                
        if count == 1:
            X,Y = b(current_lons,current_lats)
            for i in range(len(current_lons)):
                #convert amplitude and phase to complex real and imag parts
                obs_complex = cmath.rect(obs_seasonal_mag[i], obs_seasonal_phase[i])
                all1 = b.quiver(X[i],Y[i],obs_complex.real,obs_complex.imag,width=0.005,color=cmap((1./max_cb)*obs_ave[i]))
                
        if count == 2:
            X,Y = c(current_lons,current_lats)
            for i in range(len(current_lons)):
                #convert amplitude and phase to complex real and imag parts
                obs_complex = cmath.rect(obs_seasonal_mag[i], obs_seasonal_phase[i])
                all1 = c.quiver(X[i],Y[i],obs_complex.real,obs_complex.imag,width=0.005,color=cmap((1./max_cb)*obs_ave[i]))
                
        if count == 3:
            X,Y = d(current_lons,current_lats)
            for i in range(len(current_lons)):
                #convert amplitude and phase to complex real and imag parts
                obs_complex = cmath.rect(obs_seasonal_mag[i], obs_seasonal_phase[i])
                all1 = d.quiver(X[i],Y[i],obs_complex.real,obs_complex.imag,width=0.005,color=cmap((1./max_cb)*obs_ave[i]))

        count+=1
    
plt.tight_layout(pad = 1.08)
fig.subplots_adjust(wspace=0.01)
fig.subplots_adjust(hspace=0.01)

plt.show()