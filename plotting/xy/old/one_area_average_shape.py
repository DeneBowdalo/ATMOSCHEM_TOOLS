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
    
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)


#--------------------------------------------------------
#load in periodic lsp data
obs_period_grp = Dataset('../obs_%s_%s/obs_sig_periods.nc'%(vres,timeres))
model_period_grp = Dataset('model_sig_periods.nc')

obs_daily_waveforms = []
obs_seasonal_waveforms = []
obs_full_waveforms = []
model_daily_waveforms = []
model_seasonal_waveforms = []
model_full_waveforms = []

for ref in obs_refs:
    site_group_obs = obs_period_grp.groups[ref]
    site_group_mod = model_period_grp.groups[ref]
    
    model_daily_waveforms.append(site_group_mod.variables['daily_waveform'][:])
    model_seasonal_waveforms.append(site_group_mod.variables['seasonal_waveform'][:])
    model_full_waveforms.append(site_group_mod.variables['all_waveform'][:])
    
    obs_daily_waveforms.append(site_group_obs.variables['daily_waveform'][:])
    obs_seasonal_waveforms.append(site_group_obs.variables['seasonal_waveform'][:])
    obs_full_waveforms.append(site_group_obs.variables['all_waveform'][:])


obs_daily_waveforms = np.array(obs_daily_waveforms)
obs_seasonal_waveforms = np.array(obs_seasonal_waveforms)
obs_full_waveforms = np.array(obs_full_waveforms)	
model_daily_waveforms = np.array(model_daily_waveforms)
model_seasonal_waveforms = np.array(model_seasonal_waveforms)
model_full_waveforms = np.array(model_full_waveforms)

#-----------------------------------
#get area
area = raw_input('ANT, S_O, S_NA, S_EU, S_AS, SW_NA, C_EU, SE_AS, SE_NA, E_EU, C_AS, CE_NA, NW_EU, NE_AS, C_NA, N_EU, OC, NW_NA, AF, NE_NA, SA, AL, N_O or ARC?\n')
plot_type = raw_input('\nd or s?\n')

ave_or_no = raw_input('\nave? y or n?\n')

if plot_type == 'd':
    obs_datetimes = obs_datetimes[:24]
    model_datetimes = model_datetimes[:24]
if plot_type == 's':
    obs_datetimes = obs_datetimes[:8766]
    model_datetimes = model_datetimes[:8766]

obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')

area_boundaries,area_tags,area_labels = modules.area_dicts()

#get time series and waveforms in cut region

area_grid = area_boundaries[area]
area_tag = area_tags[area]
area_label = area_labels[area]

cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)

test = obs_daily_waveforms < 0
obs_daily_waveforms[test] = np.nan

obs_daily_waveforms = np.nanmean(obs_daily_waveforms[cut_test,:],axis=0)
obs_seasonal_waveforms = np.nanmean(obs_seasonal_waveforms[cut_test,:],axis=0)

model_daily_waveforms = np.nanmean(model_daily_waveforms[cut_test,:],axis=0)
model_seasonal_waveforms = np.nanmean(model_seasonal_waveforms[cut_test,:],axis=0)

if plot_type == 'd':
    ave_obs_waveform = obs_daily_waveforms 
    ave_model_waveform = model_daily_waveforms    
if plot_type == 's':
    ave_obs_waveform = obs_seasonal_waveforms
    ave_model_waveform = model_seasonal_waveforms
    
if ave_or_no == 'n':
    ave_obs_waveform = ave_obs_waveform - np.average(ave_obs_waveform)
    ave_model_waveform = ave_model_waveform - np.average(ave_model_waveform)

#--------------------------------------------------------------------
#set up plot
fig, (ax) = plt.subplots(1,figsize=(19,13))
fig.patch.set_facecolor('white')

ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = 'black',linestyle='_',linewidth=3,markeredgecolor='None')
ax.plot_date(model_time_pd.to_pydatetime(),ave_model_waveform,color = 'red',linestyle='_',linewidth=3,markeredgecolor='None')
ax.set_xlabel('Time',fontsize=25)
ax.set_ylabel('Concentration (ppb)',fontsize=25)


ax.text(0.9, 0.9, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=35)

plt.show()
