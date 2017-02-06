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
import redfit
import random
import numpy.fft
from cmath import *
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import scipy
import numpy.fft as FFT
from netCDF4 import Dataset

pi2 = np.pi*2


#read in obs data
root_grp = Dataset('../GAW_O3/GAW_SURFACE_O3_2006_2012.nc')

#read in specific site data
site_group = root_grp.groups['RCV']
#site_group = root_grp.groups['mnm']

#read in variables for site
obs_var = site_group.variables['o3'][:]
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude

obs_times = modules.date_process(obs_date,obs_time)
obs_times = np.array(obs_times)
obs_times_full = obs_times[:]

##cut out invalid obs data
#valids = obs_var > 0
#obs_var = obs_var[valids]
#obs_times = obs_times[valids]
#obs_date = obs_date[valids]
#obs_time = obs_time[valids]

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

full_model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]
full_model_datetimes = np.array(full_model_datetimes)
#obs_datetimes = full_obs_datetimes[valids]


#process model
#----------------------------------------
#read in model data
model_dict = {'4x5':'../binary_logs/4x5_GRID_O3.npy','2x2.5':'../binary_logs/2x2.5_GRID_O3.npy'}
model_version = raw_input('\nChoose Model Version.\n%s\n'%('   '.join(i for i in model_dict)))
model_f = model_dict[model_version]

model_data = read = np.load(model_f)
if model_version == '2x2.5':
    model_time = np.arange(0,2191,1./24)
else:
    model_time = np.arange(0,2190,1./24)
    full_model_datetimes = full_model_datetimes[1:]

#get model grid dims. for sim. type
lat_c,lat_e,lon_c,lon_e = modules.model_grids(model_version)
gridbox_count = len(lat_c)*len(lon_c)

#get model gridbox for obs site
gridbox_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)

model_var = model_data[gridbox_n::gridbox_count]
model_var = model_var*1e9

model_ave = np.mean(model_var)

#make time start from 0    
model_times_from0 = modules.phase_start_correct(model_time)

periods, mag, ph,fr,fi,amp_corr= modules.take_lomb(model_times_from0,model_var,4,1./24)

zoomfact = 1000

daily_amp,daily_phase = modules.periodic_interp(fr,fi,zoomfact,periods,1.,len(obs_var),amp_corr)
ha_amp,ha_phase = modules.periodic_interp(fr,fi,zoomfact,periods,365.25/2.,len(obs_var),amp_corr)
annual_amp,annual_phase = modules.periodic_interp(fr,fi,zoomfact,periods,365.25,len(obs_var),amp_corr)

#correct for phase shift from sites where raw times do not start from 0
daily_phase = modules.phase_start_point_correct(1.,daily_phase,obs_times)
ha_phase = modules.phase_start_point_correct(365.25/2.,ha_phase,obs_times)
annual_phase = modules.phase_start_point_correct(365.25,annual_phase,obs_times)

lon_step_time  = 24./360.

#convert site_lon to 0 to 360 degs
if obs_lon < 0:
    obs_lon = 360-np.abs(obs_lon)

#transform from UTC time to solar time 
sun_time = lon_step_time*obs_lon
time_diff = sun_time - 0
if time_diff > 12:
    time_diff = time_diff-24

#convert daily phase from UTC to solar time
daily_phase = daily_phase + ((pi2/24.)*time_diff)
if daily_phase >= np.pi:
    daily_phase = -np.pi+(daily_phase-np.pi)
if daily_phase < -np.pi:
    daily_phase = np.pi-(np.abs(daily_phase)-np.pi)

#convert phase to time
daily_phase_time = modules.convert_phase_units_actual_single(daily_phase,24)
ha_phase_time = modules.convert_phase_units_actual_single(ha_phase,6)
annual_phase_time = modules.convert_phase_units_actual_single(annual_phase,12)

daily_wave = daily_amp*(np.cos((pi2*model_time/1.)-(daily_phase)))
ha_wave = ha_amp*(np.cos((pi2*model_time/(365.25/2.))-(ha_phase)))
annual_wave = annual_amp*(np.cos((pi2*model_time/365.25)-(annual_phase)))

big_wave = daily_wave+ha_wave+annual_wave
big_wave = big_wave+model_ave

daily_wave = daily_wave+model_ave
ha_wave = ha_wave+model_ave
annual_wave = annual_wave+model_ave


#calc variance captured
var_raw = np.var(model_var)
var_lsp = np.var(big_wave)

print 'var raw = ',var_raw
print 'var lsp = ',var_lsp

var_cap = (100./var_raw)*var_lsp
print 'var cap = ',var_cap

print daily_phase_time

fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

ax.grid(True)
ax.plot(full_model_datetimes,model_var,alpha=0.3,color='black')
ax.plot(full_model_datetimes,big_wave,alpha=0.5,label='Diurnal + Half-Annual + Annual',color='yellow')
ax.plot(full_model_datetimes,daily_wave,alpha = 0.6,linestyle = ':',color='green',linewidth=8,label='Diurnal')
ax.plot(full_model_datetimes,ha_wave,linestyle = '-.',color='blue',linewidth=8,label='Half-Annual')
ax.plot(full_model_datetimes,annual_wave,linestyle = '--',color='red',linewidth=8,label='Annual')
ax.set_xlabel('Time',fontsize = 15)
ax.set_ylabel('Concentration (ppb)',fontsize = 15)
ax.tick_params(axis='both', which='major', labelsize=15)
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(11.)

leg=ax.legend(loc=0, prop={'size':20})
leg.get_frame().set_alpha(0.4)

fig.tight_layout()

plt.show()



