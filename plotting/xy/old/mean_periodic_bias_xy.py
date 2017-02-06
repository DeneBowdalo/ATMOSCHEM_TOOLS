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

def interactive(event):
    modules.clicker_interactive_map_obsmodel(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_daily_waveform,obs_seasonal_waveform,obs_full_waveform,model_daily_waveform,model_seasonal_waveform,model_full_waveform,fig,all_m)

#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

model_ts_grp = Dataset(model_fname)
model_var = model_ts_grp.variables[species.lower()][:]
model_date = model_ts_grp.variables['date'][:]
model_time = model_ts_grp.variables['time'][:]
lat_e = model_ts_grp.variables['lat_edges'][:]
lon_e = model_ts_grp.variables['lon_edges'][:]
lat_c = model_ts_grp.variables['lat_centre'][:]
lon_c = model_ts_grp.variables['lon_centre'][:]
grid_size = model_ts_grp.variables['grid_size'][:]
grid_size = grid_size[0]

#process model dates and model times to datetimes, then process pandas objects

year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change model times to datetimes
model_date = model_date.astype('str')
model_time = model_time.astype('str')

for date in model_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in model_time:
    if np.float64(time) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(time) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

model_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#------------------------------------------------------------
#read in obs time series data
obs_ts_grp = Dataset(obs_fname)
obs_refs_dict = obs_ts_grp.groups

obs_refs = []
obs_lats = []
obs_lons = []
obs_alt = []

for i in obs_refs_dict.keys():
    i = i.encode('ascii')
    obs_refs = np.append(obs_refs,i)

for ref in obs_refs:
    obs_ts_group = obs_ts_grp.groups[ref] 
    obs_lats = np.append(obs_lats,obs_ts_group.latitude)
    obs_lons = np.append(obs_lons,obs_ts_group.longitude)
    obs_alt = np.append(obs_alt,obs_ts_group.altitude)
    obs_date = obs_ts_group.variables['date'][:]
    obs_time = obs_ts_group.variables['time'][:]
    
for i in range(len(obs_refs)):
    obs_refs[i] = obs_refs[i].lower()
    
#process obs dates and obs times to datetimes, then process pandas objects
year_val = []
month_val = []
day_val = []
hour_val = []
minute_val = []

#change obs times to datetimes
obs_date = obs_date.astype('str')
obs_time = obs_time.astype('str')

for date in obs_date:
    year_val.append(int(date[0:4]))
    month_val.append(int(date[4:6]))
    day_val.append(int(date[6:8]))

for time in obs_time:
    if np.float64(time) == 0:
        hour_val.append(0)
        minute_val.append(0)
    elif np.float64(time) == 30:
        hour_val.append(0)
        minute_val.append(30)
    else:
        hour_val.append(int(time[0:-2]))
        minute_val.append(int(time[-2:]))

obs_datetimes = [datetime.datetime(year = year_val[i], month = month_val[i], day = day_val[i], hour = hour_val[i], minute = minute_val[i]) for i in range(len(year_val))]

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#--------------------------------------------------------
#load in periodic lsp data
obs_period_grp = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
model_period_grp = Dataset('model_sig_periods.nc')

model_daily_mag = []
model_daily_phase = []
obs_daily_mag = []
obs_daily_phase = []
model_ave = []
obs_ave = []
model_seasonal_mag = []
model_seasonal_phase = []
model_seasonal_phase_min = []
obs_seasonal_mag = []
obs_seasonal_phase = []
obs_seasonal_phase_min = [] 
obs_daily_waveform = []
model_daily_waveform = []
obs_seasonal_waveform = []                                                                                                                                                                                                                       
model_seasonal_waveform = []
obs_full_waveform = []                                                                                                                                                                                                                       
model_full_waveform = []
model_daily_ff = []
model_seasonal_ff = []
obs_daily_ff = []
obs_seasonal_ff = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_mod = model_period_grp.groups[ref]
    
    model_daily_mag = np.append(model_daily_mag,site_group_mod.daily_amplitude)
    model_daily_phase = np.append(model_daily_phase,site_group_mod.daily_phase)
    model_ave = np.append(model_ave,site_group_mod.average)
    model_seasonal_mag = np.append(model_seasonal_mag,site_group_mod.seasonal_amplitude)
    model_seasonal_phase = np.append(model_seasonal_phase,site_group_mod.seasonal_phase)
    model_seasonal_phase_min = np.append(model_seasonal_phase_min,site_group_mod.seasonal_phase_min)   
    model_daily_waveform.append(site_group_mod.variables['daily_waveform'][:])
    model_seasonal_waveform.append(site_group_mod.variables['seasonal_waveform'][:]) 
    model_full_waveform.append(site_group_mod.variables['all_waveform'][:]) 
    model_daily_ff = np.append(model_daily_ff,site_group_mod.daily_ff) 
    model_seasonal_ff = np.append(model_seasonal_ff,site_group_mod.seasonal_ff)

    obs_daily_mag = np.append(obs_daily_mag,site_group_obs.daily_amplitude)
    obs_daily_phase = np.append(obs_daily_phase,site_group_obs.daily_phase)
    obs_ave = np.append(obs_ave,site_group_obs.average)
    obs_seasonal_mag = np.append(obs_seasonal_mag,site_group_obs.seasonal_amplitude)
    obs_seasonal_phase = np.append(obs_seasonal_phase,site_group_obs.seasonal_phase)
    obs_seasonal_phase_min = np.append(obs_seasonal_phase_min,site_group_obs.seasonal_phase_min)
    obs_daily_waveform.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveform.append(site_group_obs.variables['seasonal_waveform'][:]) 
    obs_full_waveform.append(site_group_obs.variables['all_waveform'][:])
    obs_daily_ff = np.append(obs_daily_ff,site_group_obs.daily_ff) 
    obs_seasonal_ff = np.append(obs_seasonal_ff,site_group_obs.seasonal_ff)


tags = modules.get_tags(obs_refs)
loc_dict = {'ANT':'Antarctica','ARC':'Arctic','AF':'Africa','AS':'Asia','EU':'Europe','OC':'Oceania','O':'Oceanic Sites','NA':'North America','SA':'South America'}
loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','AS':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}

max_ave_diff = np.abs(model_ave - obs_ave)/np.max(np.abs(model_ave - obs_ave))
max_mag_daily_diff = np.abs(model_daily_mag - obs_daily_mag)/np.max(np.abs(model_daily_mag - obs_daily_mag))
max_mag_seasonal_diff = np.abs(model_seasonal_mag - obs_seasonal_mag)/np.max(np.abs(model_seasonal_mag - obs_seasonal_mag))
max_phmax_daily_diff = np.abs(model_daily_phase - obs_daily_phase)/np.max(np.abs(model_daily_phase - obs_daily_phase))
max_phmax_seasonal_diff = np.abs(model_seasonal_phase - obs_seasonal_phase)/np.max(np.abs(model_seasonal_phase - obs_seasonal_phase))
max_ff_daily_diff = np.abs(model_daily_ff - obs_daily_ff)/np.max(np.abs(model_daily_ff - obs_daily_ff))
max_ff_seasonal_diff = np.abs(model_seasonal_ff - obs_seasonal_ff)/np.max(np.abs(model_seasonal_ff - obs_seasonal_ff))

daily_cs = np.average([max_mag_daily_diff,max_phmax_daily_diff,max_ff_daily_diff],axis=0)
seasonal_cs = np.average([max_mag_seasonal_diff,max_phmax_seasonal_diff,max_ff_seasonal_diff],axis=0) 

for i in range(len(tags)):
    plt.plot(max_ave_diff[i],seasonal_cs[i],color=loc_colors[tags[i]],marker='o',linestyle='None')
plt.show()
