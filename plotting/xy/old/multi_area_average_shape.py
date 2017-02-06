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

#-----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year = modules.get_model_info(present_dir)

all_split = model_fname.split('/')
all_split2 = all_split[-1].split('_')
model_name = all_split2[0]

obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

print obs_refs

#get observational location tags 
#EU = europe, AF = africa, NA  = north america, SA = south america, ANT = antarctica, ARC = arctic, O = oceanic, OC = oceania, AS = asia
tags = modules.get_tags(obs_refs)

#if 'rig' in obs_refs:
#    obs_refs[obs_refs.index('rig')] = 'rig_photo'

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
areas = ['ANT','S_O','OC','AF','SA','NE_NA','CE_NA','SE_NA','S_NA','SW_NA','NW_NA','C_NA','S_EU','C_EU','NW_EU','N_EU','E_EU','NE_AS','SE_AS','C_AS','S_AS','N_O','AL','ARC']
plot_type = raw_input('\nd or s?\n')

ave_or_no = raw_input('\nave? y or n?\n')

if plot_type == 'd':
    obs_datetimes = obs_datetime_time[:24]
    model_datetimes = model_datetime_time[:24]
if plot_type == 's':
    obs_datetimes = obs_datetime_time[:8766]
    model_datetimes = model_datetime_time[:8766]

obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')

area_boundaries,area_tags,area_labels = modules.area_dicts()

fig, axes = plt.subplots(nrows=5, ncols=5,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0
for ax in axes.flat:
    try:
        area = areas[count]
    except:
        ax.axis('off')
        continue

    print area

    area_grid = area_boundaries[area]
    area_tag = area_tags[area]
    area_label = area_labels[area]

    cut_test = modules.area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag)

    obs_d_w = np.nanmean(obs_daily_waveforms[cut_test,:],axis=0)
    obs_s_w = np.nanmean(obs_seasonal_waveforms[cut_test,:],axis=0)

    model_d_w = np.nanmean(model_daily_waveforms[cut_test,:],axis=0)
    model_s_w = np.nanmean(model_seasonal_waveforms[cut_test,:],axis=0)

    if plot_type == 'd':
        ave_obs_waveform = obs_d_w
        ave_model_waveform = model_d_w    
    if plot_type == 's':
        ave_obs_waveform = obs_s_w
        ave_model_waveform = model_s_w

    if ave_or_no == 'n':
        ave_obs_waveform = ave_obs_waveform - np.average(ave_obs_waveform)
        ave_model_waveform = ave_model_waveform - np.average(ave_model_waveform)

    ax.plot_date(obs_time_pd.to_pydatetime(),ave_obs_waveform,color = 'black',linestyle='_',linewidth=1,markersize=1,markeredgecolor='None')
    ax.plot_date(obs_time_pd.to_pydatetime(),ave_model_waveform,color = 'red',linestyle='_',linewidth=1,markersize=1,markeredgecolor='None')
    #ax.set_xlabel('Time',fontsize=25)
    #ax.set_ylabel('Concentration (ppb)',fontsize=25)
    ax.text(0.87, 0.91, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    ax.fill_between(obs_time_pd.to_pydatetime(),ave_obs_waveform, ave_model_waveform,where=ave_model_waveform>ave_obs_waveform, facecolor='yellow', interpolate=True)
    ax.fill_between(obs_time_pd.to_pydatetime(),ave_obs_waveform, ave_model_waveform,where=ave_obs_waveform>ave_model_waveform, facecolor='blue', interpolate=True)
    if plot_type == 'd':
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H'))
    if plot_type == 's':
        ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    
    tst = False
    for t in range(len(ave_obs_waveform)):
        b = np.isnan(ave_obs_waveform[t])
        if b == False:
            tst= True
    if tst == False:
        ax.plot_date(obs_time_pd.to_pydatetime(),len(obs_time_pd)*[1],markersize=0.000001)
    
    count+=1
    
plt.tight_layout(pad = 3.08)

h1, = ax.plot([1,1],color='black',marker='o',linestyle='None',markersize=10)
h2, = ax.plot([1,1],color='red',marker='o',linestyle='None',markersize=10)

plt.legend((h1,h2),['Observations',model_name],loc='lower left',prop={'size':20},fancybox=True,ncol=1,markerscale=1,bbox_to_anchor=(-0.25,0))
h1.set_visible(False)
h2.set_visible(False)

plt.show()
