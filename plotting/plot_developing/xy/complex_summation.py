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

#-----------------------------------------------------
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
    obs_site_group = obs_ts_grp.groups[ref] 
    obs_lats = np.append(obs_lats,obs_site_group.latitude)
    obs_lons = np.append(obs_lons,obs_site_group.longitude)
    obs_alt = np.append(obs_alt,obs_site_group.altitude)
    obs_date = obs_site_group.variables['date'][:]
    obs_time = obs_site_group.variables['time'][:]
    
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

obs_dh3_amp = []
obs_dh2_amp = []
obs_dh1_amp = []
obs_df_amp = []
obs_dh3_ph = []
obs_dh2_ph = []
obs_dh1_ph = []
obs_df_ph = []
obs_sh3_amp = []
obs_sh2_amp = []
obs_sh1_amp = []
obs_sf_amp = []
obs_sh3_ph = []
obs_sh2_ph = []
obs_sh1_ph = []
obs_sf_ph = []
obs_ave = []
obs_d_mag = []
obs_d_ph = []
obs_s_mag = []
obs_s_ph = []
model_dh3_amp = []
model_dh2_amp = []
model_dh1_amp = []
model_df_amp = []
model_dh3_ph = []
model_dh2_ph = []
model_dh1_ph = []
model_df_ph = []
model_sh3_amp = []
model_sh2_amp = []
model_sh1_amp = []
model_sf_amp = []
model_sh3_ph = []
model_sh2_ph = []
model_sh1_ph = []
model_sf_ph = []
model_ave = []
model_d_mag = []
model_d_ph = []
model_s_mag = []
model_s_ph = []


for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_mod = model_period_grp.groups[ref]

    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_dh3_amp.append(obs_site_group.daily_harmonic3_amplitude)
    obs_dh2_amp.append(obs_site_group.daily_harmonic2_amplitude)
    obs_dh1_amp.append(obs_site_group.daily_harmonic1_amplitude)
    obs_df_amp.append(obs_site_group.original_daily_amplitude)
    obs_dh3_ph.append(obs_site_group.daily_harmonic3_phase)
    obs_dh2_ph.append(obs_site_group.daily_harmonic2_phase)
    obs_dh1_ph.append(obs_site_group.daily_harmonic1_phase)
    obs_df_ph.append(obs_site_group.original_daily_phase)
    obs_sh3_amp.append(obs_site_group.seasonal_harmonic3_amplitude)
    obs_sh2_amp.append(obs_site_group.seasonal_harmonic2_amplitude)
    obs_sh1_amp.append(obs_site_group.seasonal_harmonic1_amplitude)
    obs_sf_amp.append(obs_site_group.annual_amplitude)
    obs_sh3_ph.append(obs_site_group.seasonal_harmonic3_phase)
    obs_sh2_ph.append(obs_site_group.seasonal_harmonic2_phase)
    obs_sh1_ph.append(obs_site_group.seasonal_harmonic1_phase)
    obs_sf_ph.append(obs_site_group.annual_phase)
    
    obs_ave.append(obs_site_group.average)
    obs_d_mag.append(obs_site_group.daily_amplitude)
    obs_d_ph.append(obs_site_group.daily_phase)
    obs_s_mag.append(obs_site_group.seasonal_amplitude)
    obs_s_ph.append(obs_site_group.seasonal_phase)
    
    model_site_group.append(model_period_grp.groups[ref]
    model_dh3_amp.append(model_site_group.daily_harmonic3_amplitude)
    model_dh2_amp.append(model_site_group.daily_harmonic2_amplitude)
    model_dh1_amp.append(model_site_group.daily_harmonic1_amplitude)
    model_df_amp.append(model_site_group.original_daily_phase)
    model_dh3_ph.append(model_site_group.daily_harmonic3_phase)
    model_dh2_ph.append(model_site_group.daily_harmonic2_phase)
    model_dh1_ph.append(model_site_group.daily_harmonic1_phase)
    model_df_ph.append(model_site_group.original_daily_phase)
    model_sh3_amp.append(model_site_group.seasonal_harmonic3_amplitude)
    model_sh2_amp.append(model_site_group.seasonal_harmonic2_amplitude)
    model_sh1_amp.append(model_site_group.seasonal_harmonic1_amplitude)
    model_sf_amp.append(model_site_group.annual_amplitude)
    model_sh3_ph.append(model_site_group.seasonal_harmonic3_phase)
    model_sh2_ph.append(model_site_group.seasonal_harmonic2_phase)
    model_sh1_ph.append(model_site_group.seasonal_harmonic1_phase)
    model_sf_ph.append(model_site_group.annual_phase)

    model_d_mag.append(model_site_group.daily_amplitude)
    model_d_ph.append(model_site_group.daily_phase)
    model_s_mag.append(model_site_group.seasonal_amplitude)
    model_s_ph.append(model_site_group.seasonal_phase)
    model_ave.append(model_site_group.average)
    
obs_sh3_ph,model_sh3_ph = convert_phase_rect(obs_sh3_ph,model_sh3_ph,12.)
obs_sh2_ph,model_sh2_ph = convert_phase_rect(obs_sh2_ph,model_sh2_ph,12.)
obs_sh1_ph,model_sh1_ph = convert_phase_rect(obs_sh1_ph,model_sh1_ph,12.)
obs_sf_ph,model_sf_ph = convert_phase_rect(obs_sf_ph,model_sf_ph,12.)

for j in range(len(obs_sh3_ph))
    obs_sh3_rect = cmath.rect(obs_sh3_amp[j], obs_sh3_ph[j])
    obs_sh2_rect = cmath.rect(obs_sh2_amp[j], obs_sh2_ph[j])
    obs_sh1_rect = cmath.rect(obs_sh1_amp[j], obs_sh1_ph[j])
    obs_sf_rect = cmath.rect(obs_sf_amp[j], obs_sf_ph[j])
    model_sh3_rect = cmath.rect(model_sh3_amp[j], model_sh3_ph[j])
    model_sh2_rect = cmath.rect(model_sh2_amp[j], model_sh2_ph[j])
    model_sh1_rect = cmath.rect(model_sh1_amp[j], model_sh1_ph[j])
    model_sf_rect = cmath.rect(model_sf_amp[j], model_sf_ph[j])
    
    cs1 = model_sh3_rect+model_sh2_rect
    cs2 = cs1+model_sh1_rect
    cs3 = cs2+model_sf_rect
    
    print cs1,cs2,cs3
    
    
    
