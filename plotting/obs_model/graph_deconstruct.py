import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from cmath import *
from math import *
import csv
import datetime
from scipy import signal
import multiprocessing
import datetime
import time
import modules
import random
import numpy.fft
from cmath import *
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import scipy
import numpy.fft as FFT
from netCDF4 import Dataset
from matplotlib import dates

pi2 = np.pi*2

first_year = 2006
last_year = 2011

#read in obs data
root_grp = Dataset('GLOBAL_SURFACE_O3_2006_2012.nc')
valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

site_ref = raw_input('Choose site from list. Sites with full set of yearly files between %i & %i are:\n%s\n'%(first_year,last_year+1,'   '.join(i for i in valid_refs)))

#read in specific site data
site_group = root_grp.groups[site_ref]

#read in variables for site
obs_var = site_group.variables['o3'][:]
full_obs_var = obs_var[:]
full_obs_var_mask = np.ma.masked_where(full_obs_var<=0,full_obs_var)
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude

obs_times = modules.date_process(obs_date,obs_time)
obs_times = np.array(obs_times)
obs_times_full = obs_times[:]

##cut out invalid obs data
obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
valids = obs_var > 0
obs_var = obs_var[valids]
obs_times = obs_times[valids]

obs_ave = np.average(obs_var)

year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change model times to datetimes
obs_date = obs_date.astype('str')
obs_time = obs_time.astype('str')

for date in obs_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in obs_time:
    if time == '0':
        hour_val.append(0)
        minute_val.append(0)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

full_obs_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]
full_obs_datetimes = np.array(full_obs_datetimes)
obs_datetimes = full_obs_datetimes[valids]

#read model vals

#lat lon edges for 2x2.5 grid

lat_c = np.arange(-88.,89.,2)
lat_c = np.insert(lat_c,0,-89.5)
lat_c = np.append(lat_c,89.5)

lon_c = np.arange(-180,178,2.5)

lat_e = np.arange(-89.,90,2)
lat_e = np.insert(lat_e,0,-90.)
lat_e = np.append(lat_e,90.)

lon_e = np.arange(-181.25,179,2.5)

#return gridbox number required from lat lon of obs

#convert obs_lons to same grid as model 
if obs_lon > 178.75:
    diff = obs_lon - 178.75
    obs_lon = -181.25 + diff

lat_indices = []
lon_indices = []
gridbox_num = []

lat_i = np.searchsorted(lat_e,obs_lat,side='left')
lon_i = np.searchsorted(lon_e,obs_lon,side='left')

lat_index = lat_i-1
lon_index = lon_i-1

#find gridbox number of model for obs data
if lat_index > 0:
    lat_index = 144*lat_index

gridbox_num = np.append(gridbox_num, lat_index + lon_index)

f = '/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_O3.npy'
read = np.load(f)
model_var =  read[gridbox_num::13104]
model_var = model_var*1e9
print model_var
model_ave = np.average(model_var)



daily_obs_mags = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/obs/obs_magnitudes/obs_daily_magnitudes.npy')                                                                                                   
ha_obs_mags = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/obs/obs_magnitudes/obs_half_annual_magnitudes.npy')
annual_obs_mags = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/obs/obs_magnitudes/obs_annual_magnitudes.npy')
daily_obs_phases = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/obs/obs_phases/obs_daily_phases.npy')                                                                                                   
ha_obs_phases = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/obs/obs_phases/obs_half_annual_phases.npy')
annual_obs_phases = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/obs/obs_phases/obs_annual_phases.npy')

daily_model_mags = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/2x2.5_same_gaps/model_magnitudes/model_daily_magnitudes.npy')
ha_model_mags = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/2x2.5_same_gaps/model_magnitudes/model_half_annual_magnitudes.npy')
annual_model_mags = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/2x2.5_same_gaps/model_magnitudes/model_annual_magnitudes.npy')
daily_model_phases = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/2x2.5_same_gaps/model_phases/model_daily_phases.npy')
ha_model_phases = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/2x2.5_same_gaps/model_phases/model_half_annual_phases.npy')
annual_model_phases = np.load('/work/home/db876/xy/O3/MERGE_2006_2012/2x2.5_same_gaps/model_phases/model_annual_phases.npy')

valid_refs = valid_refs.tolist()
site_index = valid_refs.index(site_ref)

daily_obs_mag = daily_obs_mags[site_index]
ha_obs_mag = ha_obs_mags[site_index]
annual_obs_mag = annual_obs_mags[site_index]
daily_obs_phase = daily_obs_phases[site_index]
ha_obs_phase = ha_obs_phases[site_index]
annual_obs_phase = annual_obs_phases[site_index]

daily_model_mag = daily_model_mags[site_index]
ha_model_mag = ha_model_mags[site_index] 
annual_model_mag = annual_model_mags[site_index]
daily_model_phase = daily_model_phases[site_index]
ha_model_phase = ha_model_phases[site_index]
annual_model_phase = annual_model_phases[site_index]

if annual_obs_phase >= 6:
    peak_1 = ha_obs_phase+6
    peak_2 = ha_obs_phase
    diff_1 = np.abs(annual_obs_phase - peak_1)
    diff_2 = np.abs(annual_obs_phase - peak_2)
    if diff_2 > diff_1:
        p_1 = peak_1
        p_2 = peak_2
    else:
        p_1 = peak_2                                                                                                                                                                                                              
        p_2 = peak_1
else:
    peak_1 = ha_obs_phase
    peak_2 = ha_obs_phase+6
    diff_1 = np.abs(annual_obs_phase - peak_1)
    diff_2 = np.abs(annual_obs_phase - peak_2)
    if diff_2 > diff_1:
        p_1 = peak_1
        p_2 = peak_2
    else:
        p_1 = peak_2
        p_2 = peak_1 

print 'HA Peak = ', ha_obs_phase
print 'Annual Peak = ', annual_obs_phase
print 'Peak 1 = ', p_1
print 'Peak 2 = ', p_2
ra = 100./annual_obs_mag
print ha_obs_mag
print annual_obs_mag
print 'HA % of Annual Peak = ', ra*ha_obs_mag

#convert phases to radians
calc = pi2/24.
daily_obs_phase = daily_obs_phase * calc
daily_model_mag = daily_model_mag * calc
calc = pi2/6.
ha_obs_phase = ha_obs_phase * calc
ha_model_phase = ha_model_phase * calc
calc = pi2/12.
annual_obs_phase = annual_obs_phase * calc
annual_model_phase = annual_model_phase * calc


#obs waveform
daily_obs_wave = daily_obs_mag*(np.cos((pi2*obs_times_full/1.)-(daily_obs_phase)))
ha_obs_wave = ha_obs_mag*(np.cos((pi2*obs_times_full/(365.25/2.))-(ha_obs_phase)))
print annual_obs_mag,annual_obs_phase
annual_obs_wave = annual_obs_mag*(np.cos((pi2*obs_times_full/365.25)-(annual_obs_phase)))

big_obs_wave = daily_obs_wave+ha_obs_wave+annual_obs_wave
big_obs_wave = big_obs_wave+obs_ave

daily_obs_wave = daily_obs_wave+obs_ave
ha_obs_wave = ha_obs_wave+obs_ave
annual_obs_wave = annual_obs_wave+obs_ave

big_obs_wave = np.ma.masked_where(np.isnan(big_obs_wave) ,big_obs_wave)

print type(big_obs_wave)
print big_obs_wave
print np.ma.average(big_obs_wave)
print np.min(big_obs_wave)
print np.max(big_obs_wave)

#calc variance captured
obs_var_raw = np.ma.var(full_obs_var_mask)
obs_var_daily = np.ma.var(daily_obs_wave)
obs_var_ha = np.ma.var(ha_obs_wave)
obs_var_annual = np.ma.var(annual_obs_wave)

obs_var_lsp = np.ma.var(big_obs_wave)

daily_obs_var_cap = np.around((100./obs_var_raw)*obs_var_daily,2)
ha_obs_var_cap = np.around((100./obs_var_raw)*obs_var_ha,2)
annual_obs_var_cap = np.around((100./obs_var_raw)*obs_var_annual,2)
total_obs_var_cap = np.around((100./obs_var_raw)*obs_var_lsp,2)

obs_g = [daily_obs_var_cap,ha_obs_var_cap,annual_obs_var_cap,total_obs_var_cap]

#model waveform
daily_model_wave = daily_model_mag*(np.cos((pi2*obs_times_full/1.)-(daily_model_phase)))
ha_model_wave = ha_model_mag*(np.cos((pi2*obs_times_full/(365.25/2.))-(ha_model_phase)))
annual_model_wave = annual_model_mag*(np.cos((pi2*obs_times_full/365.25)-(annual_model_phase)))

big_model_wave = daily_model_wave+ha_model_wave+annual_model_wave
big_model_wave = big_model_wave+model_ave

daily_model_wave = daily_model_wave+model_ave
ha_model_wave = ha_model_wave+model_ave
annual_model_wave = annual_model_wave+model_ave

big_model_wave = np.ma.masked_where(np.isnan(big_model_wave) ,big_model_wave)

#calc variance captured
model_var_raw = np.ma.var(model_var)
model_var_daily = np.ma.var(daily_model_wave)
model_var_ha = np.ma.var(ha_model_wave)
model_var_annual = np.ma.var(annual_model_wave)

model_var_lsp = np.ma.var(big_model_wave)

daily_model_var_cap = np.around((100./model_var_raw)*model_var_daily,2)
ha_model_var_cap = np.around((100./model_var_raw)*model_var_ha,2)
annual_model_var_cap = np.around((100./model_var_raw)*model_var_annual,2)
total_model_var_cap = np.around((100./model_var_raw)*model_var_lsp,2)

model_g = [daily_model_var_cap,ha_model_var_cap,annual_model_var_cap,total_model_var_cap]

#plot all
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')

ind = np.arange(4)
width = 0.325

labels = ['Diurnal','Half-Annual','Annual','Total Variance Captured']

ax = fig.add_subplot(2,1,1)

ax.bar(ind+width, obs_g, width, color=['b','g','r','black'])
#ax.set_xticks(ind+(width*1.5))
ax.set_ylim(0,100)
ax.set_xticks([0.5,1.5,2.5,3.5])
ax.set_xticklabels(labels,fontsize=15)
ax.grid(True)

ax1 = fig.add_subplot(2,1,2)
#ax1.set_xticks(ind+(width*1.5))
ax1.bar(ind+width, model_g, width, color=['b','g','r','black'])
ax1.set_ylim(0,100)
ax1.set_xticks([0.5,1.5,2.5,3.5])
ax1.set_xticklabels(labels,fontsize=15)
ax1.grid(True)

ax.annotate('Obs.', xy=(.01,.93), xycoords='axes fraction', fontsize =20)
ax1.annotate('GEOS 2x2.5', xy=(.01,.93), xycoords='axes fraction', fontsize =20)

ax.set_title('%s'%(site_ref),fontsize=25)

plt.show()



